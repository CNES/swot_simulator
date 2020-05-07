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


class Parameters:
    """
    Simulator parameter management.

    Args:
        path (str, optional): Path to the configuration file used to override
        the default settings
    """
    CONFIG_VALUES: Dict[str, Tuple[Any, Any]] = dict(
        area=(None, [float, float, float, float]),
        cycle_duration=(20.86455, float),
        delta_al=(2.0, float),
        delta_ac=(2.0, float),
        ephemeris=(None, str),
        ephemeris_cols=(None, [int, int, int]),
        complete_product=(True, bool),
        half_gap=(10.0, float),
        half_swath=(60.0, float),
        height=(891000, float),
        nadir=(False, bool),
        swath=(True, bool),
        ssh_plugin=(None, ssh.Interface),
        shift_lon=(None, float),
        shift_time=(None, float),
        noise=(True, bool),
        len_repeat=(20000, float),
        save_signal=(True, bool),
        roll_phase_file=(None, str),
        error_spectrum_file=(None, str),
        karin_file=(None, str),
        swh=(2, int),
        nrand_karin=(1000, int),
        nbeam=(2, int),
        sigma=(6, float),
        beam_position=((-20, 20), list),
        nseed=(0, int),
        working_directory=(DEFAULT_WORKING_DIRECTORY, str),
    )

    def __init__(self, overrides: Dict[str, Any]):
        if "ephemeris" not in overrides:
            raise TypeError("missing required argument: 'ephemeris'")
        self._init_user_parameters(overrides)

    def _convert_overrides(self, name: str, value: Any) -> Any:
        expected_type = self.CONFIG_VALUES[name][1]
        try:
            if isinstance(expected_type, list):
                if not isinstance(value,
                                  list) or len(value) != len(expected_type):
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

        if self.__dict__['nbeam'] not in [1, 2, 12]:
            raise ValueError(
                f"nbeam should be either 1 or 2 or 12, not {self.__dict__['nbeam']}"
            )

    @property
    def box(self) -> math.Box:
        area = self.__dict__["area"]
        if area is None:
            return math.Box()
        return math.Box(math.Point(*area[:2]), math.Point(*area[-2:]))
