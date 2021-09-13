# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Handle the product specification
================================
"""
from typing import Any, Dict, Iterator, List, Optional, Tuple, Union
import copy
import enum
import os
import xml.etree.ElementTree as xt

#
import numpy as np
import numpy.typing as npt
import xarray as xr

#
from . import EXPERT, PRODUCT_TYPE, UNSMOOTHED
from . import math


class Side(enum.Enum):
    """Represents both sides of the swath."""
    LEFT = "left"
    RIGHT = "right"


def _find(element: xt.Element, tag: str) -> xt.Element:
    """Find a tag in the xml format specifcation file"""
    result = element.find(tag)
    if result is None:
        raise RuntimeError("The XML tag '" + tag + "' doesn't exist")
    return result


def _parse_type(dtype, width, signed) -> Union[npt.DTypeLike, np.generic]:
    """Parse type from xml format specification file. """
    if dtype == "real":
        return getattr(np, "float" + width)
    if dtype == "integer":
        return getattr(np, ("u" if not signed else "") + "int" + width)
    if dtype == "string":
        return np.str_
    if dtype == "char":
        return np.dtype(f"S{width}")  # type: ignore
    raise ValueError("Data type '" + dtype + "' is not recognized.")


def cast_to_dtype(attr_value: Union[int, float], properties: Dict[str, str]):
    """Cast the attribute value to the numpy data type required"""
    return getattr(np, properties["dtype"])(attr_value)


def _strtobool(value: str) -> bool:
    """Return the boolean value encoded in the string"""
    value = value.lower()
    if value == "true":
        return True
    if value == "false":
        return False
    raise ValueError(f"invalid truth value {value!r}")


def _strip_shape(name: str) -> str:
    """Strip the dimension name."""
    if name[0] == "/":
        return name[1:]
    return name


def _parser(tree: xt.ElementTree):
    """Parse variables, attributes and shapes from xml format specification
    file"""
    variables = dict()
    attributes = dict()
    shapes = dict()

    for item in tree.getroot().findall("shape"):
        dims = tuple(
            _strip_shape(dim.attrib["name"])
            for dim in item.findall("dimension"))
        if dims:
            shapes[item.attrib["name"]] = dims

    for item in _find(_find(tree.getroot(), 'science'), 'nodes'):
        dtype = _parse_type(
            item.tag, item.attrib["width"],
            _strtobool(item.attrib["signed"])
            if "signed" in item.attrib else None)
        if not isinstance(dtype, np.dtype):
            dtype = dtype.__name__  # type: ignore
        annotation = item.find("annotation")
        if annotation is None:
            continue
        varname = item.attrib["name"]
        if varname.startswith("/@"):
            attributes[varname[2:]] = dict(attrs=annotation.attrib,
                                           dtype=dtype)
        else:
            del annotation.attrib["app"]
            variables[varname[1:]] = dict(attrs=annotation.attrib,
                                          dtype=dtype,
                                          shape=shapes[item.attrib["shape"]])

    return variables, attributes


def parse_specification_file(path: str) -> Tuple:
    """Parse the XML specification file"""
    return _parser(xt.parse(path))


def encode_fill_value(properties: Dict[str, Any]) -> np.ndarray:
    """Returns the fill value encoded with the data type specified"""
    dtype = properties["dtype"]
    attrs = properties["attrs"]
    fill_value = attrs["_FillValue"] if "_FillValue" in attrs else 0
    if isinstance(dtype, str):
        return getattr(np, dtype)(fill_value)
    return np.array(fill_value, dtype)


def build_array(name: str, variables: Dict[str, Dict[str, Any]],
                data: np.ndarray):
    """Builds the array from the data and the variables"""
    properties = variables[name]
    attrs = copy.deepcopy(properties["attrs"])

    # Reading the storage properties of the variable
    encoding: Dict[str, Any] = dict(dtype=properties["dtype"])

    # If the variable defines a fill value.
    if "_FillValue" in attrs:
        encoding["_FillValue"] = encode_fill_value(properties)
        del attrs["_FillValue"]

    # Some values read from the XML files must be decoded
    # TODO(fbriol): The type of these attributes should be determined
    # from their type, but at the moment this is not possible.
    for item in ["add_offset", "scale_factor"]:
        if item in attrs:
            attrs[item] = float(attrs[item])
    for item in ["valid_range", "valid_min", "valid_max"]:
        if item in attrs:
            attrs[item] = cast_to_dtype(attrs[item], properties)
    if "flag_values" in attrs:
        items = attrs["flag_values"].split()
        attrs["flag_values"] = np.array(
            [cast_to_dtype(item, properties) for item in items],
            properties["dtype"]) if len(items) != 1 else cast_to_dtype(
                float(attrs["flag_values"]), properties)
    # if "scale_factor" in attrs and "add_offset" not in attrs:
    #     attrs["add_offset"] = 0.0
    # if "add_offset" in attrs and "scale_factor" not in attrs:
    #     attrs["scale_factor"] = 1.0
    return {
        name: encoding
    }, xr.DataArray(data=data,
                    dims=properties["shape"],
                    name=name,
                    attrs=attrs)


class ProductSpecification:
    """Parse and load into memory the product specification"""
    SPECIFICATION = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "l2b-ssh.xml")

    def __init__(self, product_type: Optional[str]):
        product_type = product_type or EXPERT
        self.unsmoothed = product_type == UNSMOOTHED
        self.specification = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            PRODUCT_TYPE[product_type])
        self.variables, self.attributes = parse_specification_file(
            self.specification)

    def __contains__(self, item):
        if self.unsmoothed:
            item = f"{Side.LEFT.value}/{item}"
        return item in self.variables

    def ssh_karin_name(self) -> str:
        """Get the name of the variable containing the SSH."""
        return "ssh_karin" + ("_2" if self.unsmoothed else "")

    def _variable(self, name: str):
        """Get the name of the template variable."""
        return self.variables[f"{Side.LEFT.value}/{name}" if self.
                              unsmoothed else name]

    def _shape(self, swath: bool = True):
        """Get the shape of the variable"""
        if swath:
            return self._variable("longitude")["shape"]
        return self._variable("time")["shape"]

    def num_sides(self) -> Dict[str, int]:
        """Get the name of the dimensions representing the number of sides."""
        if self.unsmoothed:
            return {
                f"{Side.LEFT.value}/num_sides": 1,
                f"{Side.RIGHT.value}/num_sides": 1
            }
        return {"num_sides": 2}

    def _names(self, name: str, split: bool = True) -> Tuple[str, ...]:
        """Get the names of the variables in the dataset."""
        if self.unsmoothed and split:
            return (f"{Side.LEFT.value}/{name}", f"{Side.RIGHT.value}/{name}")
        return (name, )

    def _data_array(
            self,
            name: str,
            data: np.ndarray,
            split: bool = True) -> Optional[Tuple[Dict, List[xr.DataArray]]]:
        """Returns a tuple containing the netCDF encoding of the variable and
        the data array."""
        result = {}, []

        if split and self.unsmoothed:
            if len(data.shape) == 2:
                middle = data.shape[1] // 2
                swaths = (data[:, :middle], data[:, middle:])
            else:
                swaths = (data, data)
            names = self._names(name)
        else:
            swaths = (data, )
            names = self._names(name, split=False)

        for ix, varname in enumerate(names):
            if varname not in self.variables:
                return None
            encoding_array, array = build_array(varname, self.variables,
                                                swaths[ix])
            result[0].update(encoding_array)
            result[1].append(array)
        return result if result[1] else None

    def time(self, time: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Gets the time axis"""
        return self._data_array("time", time)  # type: ignore

    def x_ac(self,
             x_ac: np.ndarray) -> Optional[Tuple[Dict, List[xr.DataArray]]]:
        """Returns the properties of the variable describing the cross track
        distance"""
        return self._data_array("cross_track_distance", x_ac * 1000)  # km -> m

    def lon(self, lon: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the longitudes of
        the swath"""
        # Longitude must be in [0, 360.0[
        result = self._data_array("longitude",
                                  math.normalize_longitude(lon, 0))
        assert result is not None
        return result

    def lon_nadir(
            self, lon_nadir: np.ndarray
    ) -> Optional[Tuple[Dict, List[xr.DataArray]]]:
        """Returns the properties of the variable describing the longitudes of
        the reference ground track"""
        # Longitude must be in [0, 360.0[
        return self._data_array("longitude_nadir",
                                math.normalize_longitude(lon_nadir, 0))

    def lat(self, lat: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the latitudes of
        the swath"""
        result = self._data_array("latitude", lat)
        assert result is not None
        return result

    def lat_nadir(
            self, lat_nadir: np.ndarray
    ) -> Optional[Tuple[Dict, List[xr.DataArray]]]:
        """Returns the properties of the variable describing the latitudes of
        the reference ground track"""
        return self._data_array("latitude_nadir", lat_nadir)

    def ssh_karin(self, ssh: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the SSH measured
        by KaRIn"""
        result = self._data_array(self.ssh_karin_name(), ssh)
        assert result is not None
        return result

    def ssh_nadir(self, array: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the SSH to
        nadir."""
        name = "ssh_nadir"
        for item in self._names(name, split=False):
            self.variables[item] = dict(attrs=dict(
                _FillValue=2147483647,
                coordinates='longitude latitude',
                long_name='sea surface height',
                scale_factor=0.0001,
                standard_name='seaurface height above reference ellipsoid',
                units='m',
                valid_min=np.int32(-15000000),
                valid_max=np.int32(150000000)),
                                        dtype='int32',
                                        shape=self._shape(swath=False))
        return self._data_array(name, array)  # type: ignore

    def simulated_true_ssh_nadir(
            self, array: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the SSH to nadir
        free of measurement errors."""
        name = "simulated_true_ssh_nadir"
        for item in self._names(name, split=False):
            self.variables[item] = dict(attrs=dict(
                _FillValue=2147483647,
                coordinates='time',
                long_name='sea surface height',
                scale_factor=0.0001,
                standard_name='sea surface height above reference ellipsoid',
                units='m',
                valid_min=np.int32(-15000000),
                valid_max=np.int32(150000000),
                comment='Height of the sea surface free of measurement errors.'
            ),
                                        dtype='int32',
                                        shape=self._shape(swath=False))
        return self._data_array(name, array)  # type: ignore

    def swh_karin(self, swh: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the SWH measured
        by KaRIn"""
        result = self._data_array("swh_karin", swh)
        assert result is not None
        return result

    def swh_nadir(self, array: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the SWH to nadir
        free of measurement errors."""
        name = "swh_nadir"
        for item in self._names(name, split=False):
            self.variables[item] = dict(attrs=dict(
                _FillValue=2147483647,
                coordinates='time',
                long_name='Significant Wave Height',
                scale_factor=0.0001,
                standard_name='Sigificant Wave Height',
                units='m',
                valid_min=np.int32(-15000000),
                valid_max=np.int32(150000000),
            ),
                                        dtype='int32',
                                        shape=self._shape(swath=False))
        return self._data_array(name, array)  # type: ignore

    def simulated_true_ssh_karin(
            self, array: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the SSH KaRIn free
        of measurement errors."""
        name = "simulated_true_ssh_karin"
        for item in self._names(name):
            self.variables[item] = dict(attrs=dict(
                _FillValue=2147483647,
                coordinates='longitude latitude',
                long_name='sea surface height',
                scale_factor=0.0001,
                standard_name='sea surface height above reference ellipsoid',
                units='m',
                valid_min=np.int32(-15000000),
                valid_max=np.int32(150000000),
                comment='Height of the sea surface free of measurement errors.'
            ),
                                        dtype='int32',
                                        shape=self._shape())
        return self._data_array(name, array)  # type: ignore

    def simulated_error_baseline_dilation(
            self, array: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the error due to
        baseline mast dilation"""
        name = "simulated_error_baseline_dilation"
        for item in self._names(name):
            self.variables[item] = dict(attrs=dict(
                _FillValue=2147483647,
                long_name='Error due to baseline mast dilation',
                scale_factor=0.0001,
                units='m',
                coordinates='longitude latitude'),
                                        dtype='int32',
                                        shape=self._shape())
        return self._data_array(name, array)  # type: ignore

    def simulated_error_roll(
            self, array: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the error due to
        roll"""
        name = "simulated_error_roll"
        for item in self._names(name):
            self.variables[item] = dict(attrs=dict(
                _FillValue=2147483647,
                long_name='Error due to roll',
                units='m',
                scale_factor=0.0001,
                coordinates='longitude latitude'),
                                        dtype='int32',
                                        shape=self._shape())
        return self._data_array(name, array)  # type: ignore

    def simulated_error_phase(
            self, array: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the error due to
        phase"""
        name = "simulated_error_phase"
        for item in self._names(name):
            self.variables[item] = dict(attrs=dict(
                _FillValue=2147483647,
                long_name='Error due to phase',
                units='m',
                scale_factor=0.0001,
                coordinates='longitude latitude'),
                                        dtype='int32',
                                        shape=self._shape())
        return self._data_array(name, array)  # type: ignore

    def simulated_error_orbital(
            self, array: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the error due to
        orbital perturbations."""
        name = "simulated_error_orbital"
        for item in self._names(name):
            self.variables[item] = dict(attrs=dict(
                _FillValue=2147483647,
                long_name='Error due to orbital perturbations',
                units='m',
                scale_factor=0.0001,
                coordinates='longitude latitude'),
                                        dtype='int32',
                                        shape=self._shape())
        return self._data_array(name, array)  # type: ignore

    def simulated_roll_phase_estimate(
            self, array: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the roll phase
        correction estimated"""
        name = "simulated_roll_phase_estimate"
        for item in self._names(name):
            self.variables[item] = dict(attrs=dict(
                _FillValue=2147483647,
                long_name='Error after estimation of roll phase',
                units='m',
                scale_factor=0.0001,
                coordinates='longitude latitude'),
                                        dtype='int32',
                                        shape=self._shape())
        return self._data_array(name, array)  # type: ignore

    def simulated_error_karin(
            self, array: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the KaRIn error"""
        name = "simulated_error_karin"
        for item in self._names(name):
            self.variables[item] = dict(attrs=dict(
                _FillValue=2147483647,
                long_name='KaRIn error',
                units='m',
                scale_factor=0.0001,
                coordinates='longitude latitude'),
                                        dtype='int32',
                                        shape=self._shape())
        return self._data_array(name, array)  # type: ignore

    def simulated_error_timing(
            self, array: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the timing
        error"""
        name = "simulated_error_timing"
        for item in self._names(name):
            self.variables[item] = dict(attrs=dict(
                _FillValue=2147483647,
                long_name='Timing error',
                units='m',
                scale_factor=0.0001,
                coordinates='longitude latitude'),
                                        dtype='int32',
                                        shape=self._shape())
        return self._data_array(name, array)  # type: ignore

    def simulated_error_troposphere(
            self, array: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the error due to
        wet troposphere path delay"""
        name = "simulated_error_troposphere"
        for item in self._names(name):
            self.variables[item] = dict(attrs=dict(
                _FillValue=2147483647,
                long_name='Error due to wet troposphere path delay',
                units='m',
                scale_factor=0.0001,
                coordinates='longitude latitude'),
                                        dtype='int32',
                                        shape=self._shape())
        return self._data_array(name, array)  # type: ignore

    def simulated_error_troposphere_nadir(
            self, array: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the error due to
        wet troposphere path delay to nadir."""
        name = "simulated_error_troposphere_nadir"
        for item in self._names(name, split=False):
            self.variables[item] = dict(attrs=dict(
                _FillValue=2147483647,
                long_name='Error due to wet troposphere path delay',
                units='m',
                scale_factor=0.0001,
                coordinates='longitude latitude'),
                                        dtype='int32',
                                        shape=self._shape(swath=False))
        return self._data_array(name, array)  # type: ignore

    def simulated_error_altimeter(
            self, array: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        """Returns the properties of the variable describing the altimeter
        error"""
        name = "simulated_error_altimeter"
        for item in self._names(name, split=False):
            self.variables[item] = dict(attrs=dict(
                _FillValue=2147483647,
                long_name='Altimeter error',
                standard_name='',
                units='m',
                scale_factor=0.0001,
                coordinates='longitude latitude'),
                                        dtype='int32',
                                        shape=self._shape(swath=False))
        return self._data_array(name, array)  # type: ignore

    def fill_variables(self, variables,
                       shape) -> Iterator[Tuple[Dict, List[xr.DataArray]]]:
        """Returns the properties of variables present in the official format
        of the SWOT product, but not calculated by this software."""
        for item, properties in self.variables.items():
            if item in variables:
                continue
            fill_value = encode_fill_value(properties)
            data = np.full(tuple(shape[dim] for dim in properties["shape"]),
                           fill_value, properties["dtype"])
            variable = self._data_array(item, data, split=False)
            if variable is not None:
                yield variable
        return StopIteration
