# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Settings handling
-----------------
"""
from typing import Any, Dict, Iterator, Tuple
import contextlib
import copy
import importlib
import logging
import os
import re
import pathlib
import traceback
import types
import numpy as np
from . import math
from . import plugins
from . import PRODUCT_TYPE

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


def error_keywords() -> Iterator[str]:
    """Get the list of the error keywords"""
    for item in error_classes():
        item = item[:1].lower() + item[1:]
        yield re.sub(r'([A-Z])', r'_\1', item).lower()


class NumberOfBeams(int):
    """Handle the number of beams"""
    def __new__(cls, value, *args, **kwargs):
        if value not in [1, 2]:
            raise ValueError("nbeam must be in [1, 2]")
        return super().__new__(cls, value, *args, **kwargs)  # type: ignore


class TimeDelta:
    """Handle a time delta in seconds in the configuration file."""
    def __call__(self, value):
        return np.timedelta64(value, "s")


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
    CONFIG_VALUES: Dict[str, Tuple[Any, Any, str]] = dict(
        area=(None, [float, 4],
              ("Geographical area to simulate defined by the minimum and "
               "maximum corner point :lon_min, lat_min, lon_max, lat_max. "
               "Default: -180, -90, 180, 90")),
        beam_position=((-20, 20), [float, 2],
                       ("Number of beam used to correct wet troposphere "
                        "signal (1, 2 or 'both')")),
        central_pixel=(False, bool,
                       ("If true, the swath, in the final dataset, will "
                        "contain a center pixel divided in half by the "
                        "reference ground track")),
        complete_product=(False, bool,
                          ("If true, the generated netCDF file will be the "
                           "complete product compliant with SWOT's Product "
                           "Description Document (PDD), otherwise only the "
                           "calculated variables will be written to the "
                           "netCDF file")),
        cycle_duration=(20.86455, float,
                        "Duration of a cycle in number of fractional days"),
        delta_ac=(2.0, float,
                  ("Distance, in km, between two points across track "
                   "direction")),
        delta_al=(2.0, float,
                  ("Distance, in km, between two points along track "
                   "direction")),
        ephemeris_cols=(None, [int, 3],
                        ("Index of columns to read in the ephemeris file "
                         "containing, respectively, longitude in degrees, "
                         "latitude in degrees and the number of seconds "
                         "elapsed since the start of the orbit. Default: "
                         "[1, 2, 0]")),
        ephemeris=(None, str,
                   ("Ephemeris file to read containing the satellite's "
                    "orbit.")),
        error_spectrum=(None, str,
                        "File containing spectrum of instrument error"),
        corrected_roll_phase_dataset=(None, str,
                                      ("Estimated roll phase dataset. "
                                       "Default: None")),
        half_gap=(10.0, float,
                  ("Distance, in km, between the nadir and the center of the "
                   "first pixel of the swath")),
        half_swath=(60.0, float,
                    ("Distance, in km, between the nadir and the center of "
                     "the last pixel of the swath")),
        height=(891000.0, float, "Satellite altitude (m)"),
        karin_noise=(None, str,
                     "KaRIN file containing spectrum for several SWH"),
        len_repeat=(20000.0, float, "Repeat length"),
        nadir=(False, bool, "True to generate Nadir products"),
        nbeam=(2, NumberOfBeams,
               ("Number of beam used to correct wet troposphere signal "
                "(1, 2 or 'both')")),
        noise=(None, [str, -1],
               ("The calculation of roll errors can be simulated, option "
                "\"roll_phase\", or interpolated, option "
                "\"corrected_roll_phase\", from the dataset specified by the "
                "option \"roll_phase_dataset\". Therefore, these two options "
                "are mutually exclusive. In other words, if the "
                "\"roll_phase\" option is present, the "
                "\"corrected_roll_phase\" option must be omitted, and vice "
                "versa")),
        nseed=(0, int, ("Seed for RandomState. Must be convertible to 32 bit "
                        "unsigned integers")),
        product_type=("expert", str,
                      ("Type of SWOT product to be generated. Possible "
                       "products are \"expert\", \"basic\" and \"wind_wave\". "
                       "Default to expert")),
        requirement_bounds=(None, [float, 2],
                            ("Limits of SWOT swath requirements. Measurements "
                             "outside the span will be set with fill values")),
        shift_lon=(None, float, "Orbit shift in longitude (degrees)"),
        shift_time=(None, TimeDelta(), "Orbit shift in time (seconds)"),
        sigma=(6.0, float, "Gaussian footprint of sigma (km)"),
        ssh_plugin=(None, plugins.Interface,
                    ("The plug-in handling the SSH interpolation under the "
                     "satellite swath")),
        swh_plugin=(None, plugins.Interface,
                    "SWH plugin to interpolate model SWH on the SWOT grid"),
        swath=(True, bool, "True to generate swath products"),
        swh=(2.0, float, "SWH for the region"),
        working_directory=(DEFAULT_WORKING_DIRECTORY, str,
                           ("The working directory. By default, files are "
                            "generated in the user's root directory")))

    #: Arguments that must be defined by the user.
    REQUIRED = ["ephemeris", "error_spectrum", "karin_noise"]

    def __init__(self, overrides: Dict[str, Any]):
        for required in self.REQUIRED:
            if required not in overrides:
                raise TypeError(f"missing required argument: {required!r}")
        self._init_user_parameters(overrides)

        product_type = getattr(self, "product_type")
        if product_type not in PRODUCT_TYPE:
            raise ValueError(f"Unknown product type: {product_type}")

        if product_type == "wind_wave":
            if getattr(self, "ssh_plugin") is not None:
                raise ValueError("The wind/wave product cannot store SSH.")

            if getattr(self, "noise") is not None:
                raise ValueError("The wind/wave product cannot store errors.")

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
        else:
            noise = []
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
                if mro == plugins.Interface:
                    value = plugins.Plugin.register(value)
                else:
                    value = expected_type(value)
        except ValueError:
            raise ValueError("invalid value %r for config value %r" %
                             (value, name))
        return value

    def _init_user_parameters(self, overrides: Dict[str, Any]):
        # To avoid side effects, default values are copied.
        settings = dict((key, copy.copy(value[0]))
                        for key, value in self.CONFIG_VALUES.items())
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


def required_settings():
    """Get the path to the simulator data."""
    data = pathlib.Path(__file__).parent.joinpath("..", "data")
    result = dict()
    for item in Parameters.REQUIRED:
        result[item] = list()
    for item in os.listdir(data):
        for data_type, values in result.items():
            if item.startswith(data_type):
                values.append(str(data.joinpath(item).resolve()))
    return result


def text_template() -> str:
    """Get the string representing the default configuration template of the
    simulator."""
    def wrap(s: str) -> Iterator[str]:
        """Cuts the line into several lines of 80 characters maximum"""
        def line(items):
            """Returns the line built"""
            return "# " + " ".join(items)

        words = []
        for item in s.split():
            if sum((len(item) for item in words)) + (len(words) - 1) * 2 > 78:
                yield line(words)
                words.clear()
            words.append(item)
        if words:
            yield line(words)

    required = required_settings()

    result = [
        "import swot_simulator.plugins.ssh",
        "import swot_simulator.plugins.swh", ""
    ]
    for key, (default, _type, help) in Parameters.CONFIG_VALUES.items():
        result += list(wrap(help))
        if key in required:
            # If there are several values for a parameter, only the first one is
            # enabled, the others are disabled.
            enable = True
            for item in required[key]:
                result += [("" if enable else "# ") + key + " = " + repr(item)]
                enable = False
            result.append("")
        elif key == "noise":
            result += [
                key + " = " + repr(
                    list(
                        set(error_keywords()) -
                        set(("corrected_roll_phase", )))), " "
            ]
        else:
            result += [key + " = " + repr(default), ""]
    # Remove the last carriage return
    result.pop()
    return "\n".join(result)


def code_template() -> Dict[str, Any]:
    """Get the code representing the default configuration template of the
    simulator."""
    required = required_settings()
    result = dict()
    for key, (default, _type, _help) in Parameters.CONFIG_VALUES.items():
        if key in required:
            result[key] = required[key][0]
        elif key == "noise":
            result[key] = list(
                set(error_keywords()) - set(("corrected_roll_phase", )))
        else:
            result[key] = default
    return result
