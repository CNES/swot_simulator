# Copyright (c) 2020 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Settings handling
-----------------
"""
from typing import Any, Dict, Iterator, Tuple
import contextlib
import importlib
import os
import logging
import traceback
import types
from . import math
from .plugins import ssh

#: Default working directory
DEFAULT_WORKING_DIRECTORY = os.path.join(os.path.expanduser('~'),
                                         "swot_simulator")

#: Module logger
LOGGER = logging.getLogger(__name__)


def execfile_(filepath: str, _globals: Any) -> None:
    """Executes a Python code defined in a file"""
    with open(filepath, 'rb') as stream:
        source = stream.read()

    code = compile(source, filepath, 'exec')
    exec(code, _globals)


@contextlib.contextmanager
def cd(target_dir: str) -> Iterator[None]:
    """Moves to a directory and returns to the working directory"""
    cwd = os.getcwd()
    try:
        os.chdir(target_dir)
        yield
    finally:
        os.chdir(cwd)


def eval_config_file(filename: str) -> Dict:
    """Evaluate a config file."""
    path = os.path.abspath(filename)
    dirname = os.path.dirname(path)

    namespace = dict(__file__=path)

    with cd(dirname):
        # during executing config file, current dir is changed to ``confdir``.
        try:
            execfile_(filename, namespace)
        except SyntaxError as err:
            raise RuntimeError(
                f"There is a syntax error in your configuration file: {err}\n")
        except SystemExit:
            raise RuntimeError(
                "The configuration file (or one of the modules it imports) "
                "called sys.exit()")
        except Exception:
            raise RuntimeError(
                "There is a programmable error in your configuration "
                f"file:\n\n{traceback.format_exc()}")

    return namespace


def error_classes() -> Iterator[str]:
    """Get the list of classes implementing random error generation."""
    module = importlib.import_module(".error", package="swot_simulator")
    for item in dir(module):
        if isinstance(getattr(module, item), type):
            yield item


class NumberOfBeams(int):
    def __init__(self, value):
        if value not in [1, 2]:
            raise ValueError("nbeam must be in [1, 2]")


class Parameters:
    """
    Simulator parameter management.

    The simulator parameters are defined in a Python file. The expected
    parameters are described in the code below.

    .. literalinclude:: ../settings.py

    Args:
        path (str, optional): Path to the configuration file used to override
            the default settings.
    """
    #: Known parameters.
    CONFIG_VALUES: Dict[str, Tuple[Any, Any]] = dict(
        area=(None, [float, 4]),
        beam_position=((-20, 20), [float, 2]),
        complete_product=(False, bool),
        cycle_duration=(20.86455, float),
        delta_ac=(2.0, float),
        delta_al=(2.0, float),
        ephemeris_cols=(None, [int, 3]),
        ephemeris=(None, str),
        error_spectrum=(None, str),
        corrected_roll_phase_dataset=(None, str),
        half_gap=(10.0, float),
        half_swath=(60.0, float),
        height=(891000, float),
        karin_noise=(None, str),
        len_repeat=(20000, float),
        nadir=(False, bool),
        nbeam=(2, NumberOfBeams),
        noise=(None, [str, -1]),
        nrand_karin=(1000, int),
        nseed=(0, int),
        shift_lon=(None, float),
        shift_time=(None, float),
        sigma=(6, float),
        ssh_plugin=(None, ssh.Interface),
        swath=(True, bool),
        swh=(2, int),
        working_directory=(DEFAULT_WORKING_DIRECTORY, str),
    )

    #: Arguments that must be defined by the user.
    REQUIRED = ["ephemeris", "error_spectrum", "karin_noise"]

    def __init__(self, overrides: Dict[str, Any]):
        for required in self.REQUIRED:
            if required not in overrides:
                raise TypeError(f"missing required argument: {required!r}")
        self._init_user_parameters(overrides)

        noise = getattr(self, "noise")
        if noise is not None:
            if "corrected_roll_phase" in noise:
                if "roll_phase" in noise:
                    raise TypeError(
                        "option 'corrected_roll_phase' not allowed with option "
                        "'roll_phase'")
                if "corrected_roll_phase_dataset" not in overrides:
                    raise TypeError("missing required argument: "
                                    "'corrected_roll_phase_dataset'")
            noise = [
                "".join(word.capitalize() for word in item.split("_"))
                for item in noise
            ]
            unknowns = set(noise) - set(error_classes())
            if unknowns:
                raise ValueError(
                    f"Unknown error generators: {', '.join(unknowns)}")
            setattr(self, "noise", noise)

    def _convert_overrides(self, name: str, value: Any) -> Any:
        expected_type = self.CONFIG_VALUES[name][1]
        try:
            if isinstance(expected_type, list):
                if not isinstance(value, list):
                    raise ValueError
                length = expected_type[1]
                if length == -1:
                    length = len(value)
                if len(value) != length:
                    raise ValueError
                expected_type = expected_type[0]
                for idx, item in enumerate(value):
                    value[idx] = expected_type(item)
            else:
                mro = expected_type.__mro__[-2]
                if mro == ssh.Interface:
                    value = ssh.Plugin.register(value)
                else:
                    value = expected_type(value)
        except ValueError:
            raise ValueError("invalid value %r for config value %r" %
                             (value, name))
        return value

    def _init_user_parameters(self, overrides: Dict[str, Any]):
        settings = dict(
            (key, value[0]) for key, value in self.CONFIG_VALUES.items())
        for name, value in overrides.items():
            if name in ["__file__", "__builtins__"] or isinstance(
                    value, types.ModuleType) or isinstance(
                        value, types.FunctionType):
                continue
            try:
                if name not in settings:
                    LOGGER.warning(
                        'unknown config value %r in override, ignoring', name)
                    continue
                settings[name] = self._convert_overrides(name, value)
            except ValueError as exc:
                LOGGER.warning("%s", exc)
        self.__dict__.update((key, value) for key, value in settings.items())

    @property
    def box(self) -> math.Box:
        """Get the geographical area of interest."""
        area = self.__dict__["area"]
        if area is None:
            return math.Box()
        return math.Box(math.Point(*area[:2]), math.Point(*area[-2:]))
