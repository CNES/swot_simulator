# Copyright (c) 2020 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Parse/Load the product specification
====================================
"""
from typing import Any, Dict, Iterator, List, Optional, Tuple, Union
import datetime
import copy
import collections
import logging
import os
import pathlib
import xml.etree.ElementTree as xt
#
import netCDF4
import numpy as np
import xarray as xr
from . import orbit_propagator
from . import math

LOGGER = logging.getLogger(__name__)

REFERENCE = "Gaultier, L., C. Ubelmann, and L.-L. Fu, 2016: The " \
    "Challenge of Using Future SWOT Data for Oceanic Field Reconstruction." \
    " J. Atmos. Oceanic Technol., 33, 119-126, doi:10.1175/jtech-d-15-0160" \
    ".1. http://dx.doi.org/10.1175/JTECH-D-15-0160.1."

GROUP = collections.OrderedDict(
    basic=dict(
        description="Basic SSH measurement data and related information for "
        "the full swath."),
    error=dict(description="Simulated measurement errors."),
    windwave=dict(
        description="Wind and wave measurement data and related information "
        "for the full swath."),
    expert=dict(
        description="Detailed contextual information, for the full swath, on "
        "the SWOT measurements; this information is intended to allow expert "
        "users to perform advanced analyses."))


def _find(element: xt.Element, tag: str) -> xt.Element:
    """ Find a tag in the xml format specifcation file"""
    result = element.find(tag)
    if result is None:
        raise RuntimeError("The XML tag '" + tag + "' doesn't exist")
    return result


def _parse_type(dtype, width, signed):
    """ Parse type from xml format specification file. """
    if dtype == "real":
        return getattr(np, "float" + width)
    elif dtype == "integer":
        return getattr(np, ("u" if not signed else "") + "int" + width)
    elif dtype == "string":
        return np.str
    elif dtype == "char":
        return np.dtype(f"S{width}")
    raise ValueError("Data type '" + dtype + "' is not recognized.")


def _cast_to_dtype(attr_value: Union[int, float], properties: Dict[str, str]):
    return getattr(np, properties["dtype"])(attr_value)


def global_attributes(attributes: Dict[str, Dict[str, Any]], cycle_number: int,
                      pass_number: int, date: np.ndarray, lng: np.ndarray,
                      lat: np.ndarray) -> Dict[str, Any]:
    def _iso_date(date: np.datetime64) -> str:
        return datetime.datetime.utcfromtimestamp(
            date.astype("datetime64[us]").astype("int64") *
            1e-6).isoformat() + "Z"

    def _iso_duration(timedelta: np.timedelta64) -> str:
        seconds = timedelta.astype("timedelta64[s]").astype("int64")
        hours = seconds // 3600
        seconds -= hours * 3600
        minutes = seconds // 60
        seconds -= minutes * 60

        result = ""
        if hours:
            result += f"{hours}H"
        if minutes or result:
            result += f"{minutes}M"
        result += f"{seconds}S"
        return "P" + result

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S : Creation")

    ellipsoid_semi_major_axis = _cast_to_dtype(
        1, attributes["ellipsoid_semi_major_axis"])
    ellipsoid_flattening = _cast_to_dtype(
        0, attributes["ellipsoid_semi_major_axis"])

    result = collections.OrderedDict({
        "Conventions":
        "CF-1.7",
        "title":
        attributes["title"]["attrs"]["description"],
        "institution":
        "CNES/JPL",
        "source":
        "Simulate product",
        "history":
        now,
        "platform":
        "SWOT",
        "references":
        REFERENCE,
        "reference_document":
        "D-56407_SWOT_Product_Description_L2_LR_SSH",
        "contact":
        "CNES aviso@altimetry.fr, JPL podaac@podaac.jpl.nasa.gov",
        "cycle_number":
        _cast_to_dtype(cycle_number, attributes["cycle_number"]),
        "pass_number":
        _cast_to_dtype(pass_number, attributes["pass_number"]),
        "time_coverage_start":
        _iso_date(date[0]),
        "time_coverage_end":
        _iso_date(date[-1]),
        "time_coverage_duration":
        _iso_duration(date[-1] - date[0]),
        "time_coverage_resolution":
        "P1S",
        "geospatial_lon_min":
        lng.min(),
        "geospatial_lon_max":
        lng.max(),
        "geospatial_lat_min":
        lat.min(),
        "geospatial_lat_max":
        lat.max()
    })
    if len(lng.shape) == 2:
        result.update({
            "left_first_longitude": lng[0, 0],
            "left_first_latitude": lat[0, 0],
            "left_last_longitude": lng[-1, 0],
            "left_last_latitude": lat[-1, 0],
            "right_first_longitude": lng[0, -1],
            "right_first_latitude": lat[0, -1],
            "right_last_longitude": lng[-1, -1],
            "right_last_latitude": lat[-1, -1],
        })
    result.update({
        "wavelength":
        _cast_to_dtype(0.008385803020979, attributes["wavelength"]),
        "orbit_solution":
        "POE",
    })
    for item in attributes:
        if item.startswith("xref_input"):
            result[item] = "N/A"
    result.update({
        "ellipsoid_semi_major_axis": ellipsoid_semi_major_axis,
        "ellipsoid_flattening": ellipsoid_flattening
    })
    return result


def _strtobool(value: str) -> bool:
    value = value.lower()
    if value == "true":
        return True
    elif value == "false":
        return False
    raise ValueError(f"invalid truth value {value!r}")


def _parser(tree: xt.ElementTree):
    """ Parse variables, attributes and shapes from xml format specification
    file"""
    variables = dict()
    attributes = dict()
    shapes = dict()

    for item in tree.getroot().findall("shape"):
        dims = [dim.attrib["name"] for dim in item.findall("dimension")]
        if dims:
            shapes[item.attrib["name"]] = tuple(dims)

    for item in _find(_find(tree.getroot(), 'science'), 'nodes'):
        dtype = _parse_type(
            item.tag, item.attrib["width"],
            _strtobool(item.attrib["signed"])
            if "signed" in item.attrib else None)
        if not isinstance(dtype, np.dtype):
            dtype = dtype.__name__
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


def _create_variable_args(encoding: Dict[str, Dict], name: str,
                          variable: xr.Variable) -> Tuple[str, Dict[str, Any]]:
    """Initiation of netCDF4.Dataset.createVariable method parameters from
    user-defined encoding information.
    """
    kwargs = dict()
    keywords = encoding[name] if name in encoding else dict()
    if "_FillValue" in keywords:
        keywords["fill_value"] = keywords.pop("_FillValue")
    dtype = keywords.pop("dtype", variable.dtype)
    for key, value in dict(zlib=True,
                           complevel=4,
                           shuffle=True,
                           fletcher32=False,
                           contiguous=False,
                           chunksizes=None,
                           endian='native',
                           least_significant_digit=None,
                           fill_value=None).items():
        kwargs[key] = keywords.pop(key, value)
    return dtype, kwargs


def _create_variable(xr_dataset: xr.Dataset, nc_dataset: netCDF4.Dataset,
                     encoding: Dict[str, Dict[str, Dict[str, Any]]], name: str,
                     unlimited_dims: Optional[List[str]],
                     variable: xr.Variable) -> None:
    """Creation and writing of the NetCDF variable"""
    unlimited_dims = unlimited_dims or list()

    variable.attrs.pop("_FillValue", None)
    # Encode datetime64 to float64
    if np.issubdtype(variable.dtype, np.datetime64):
        # 946684800000000 number of microseconds between 2000-01-01 and
        # 1970-01-01
        variable.values = (
            variable.values.astype("datetime64[us]").astype("int64") -
            946684800000000) * 1e-6
        assert (
            variable.attrs["units"] == "seconds since 2000-01-01 00:00:00.0")
    dtype, kwargs = _create_variable_args(encoding, name, variable)

    parts = name.split("/")
    name = parts.pop()
    group = parts.pop() if parts else None

    if group is not None:
        if group not in nc_dataset.groups:
            # Creating the group, dimensions and attributes
            nc_dataset = nc_dataset.createGroup(group)
            for dim_name, size in xr_dataset.dims.items():
                nc_dataset.createDimension(
                    dim_name, None if dim_name in unlimited_dims else size)

            if group in GROUP:
                nc_dataset.setncatts(GROUP[group])
        else:
            nc_dataset = nc_dataset.groups[group]
    ncvar = nc_dataset.createVariable(name, dtype, variable.dims, **kwargs)
    ncvar.setncatts(variable.attrs)
    values = variable.values
    if kwargs['fill_value'] is not None:
        if values.dtype.kind == "f" and np.any(np.isnan(values)):
            values[np.isnan(values)] = kwargs['fill_value']
        values = np.ma.array(values, mask=values == kwargs['fill_value'])
    nc_dataset[name][:] = values


def to_netcdf(dataset: xr.Dataset,
              path: Union[str, pathlib.Path],
              encoding: Optional[Dict[str, Dict]] = None,
              unlimited_dims: Optional[List[str]] = None,
              **kwargs):
    """Write dataset contents to a netCDF file"""
    encoding = encoding or dict()

    if isinstance(path, str):
        path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with netCDF4.Dataset(path, **kwargs) as stream:
        stream.setncatts(dataset.attrs)

        for name, variable in dataset.coords.items():
            _create_variable(dataset, stream, encoding, name, unlimited_dims,
                             variable)

        for name, variable in dataset.data_vars.items():
            _create_variable(dataset, stream, encoding, name, unlimited_dims,
                             variable)


class ProductSpecification:
    """Parse and load into memory the product specification"""
    SPECIFICATION = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "l2b-expert.xml")

    def __init__(self):
        self.variables, self.attributes = _parser(xt.parse(self.SPECIFICATION))

    @staticmethod
    def fill_value(properties):
        dtype = properties["dtype"]
        if isinstance(dtype, str):
            return getattr(np, dtype)(properties["attrs"]["_FillValue"])
        return np.array(properties["attrs"]["_FillValue"], dtype)

    def time(self, time: np.ndarray) -> Tuple[Dict, List[xr.DataArray]]:
        encoding, data_array = self._data_array("basic/time", time)
        return {"basic/time": encoding}, [data_array]

    def _data_array(self, name: str,
                    data: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        properties = self.variables[name]
        attrs = copy.deepcopy(properties["attrs"])

        # The fill value is casted to the target value of the variable
        fill_value = self.fill_value(properties)
        del attrs["_FillValue"]

        # Reading the storage properties of the variable ()
        encoding: Dict[str, Any] = dict(_FillValue=fill_value,
                                        dtype=properties["dtype"])

        # Some values read from the XML files must be decoded
        # TODO(fbriol): The type of these attributes should be determined from
        # their type, but at the moment this is not possible.
        for item in ["add_offset", "scale_factor"]:
            if item in attrs:
                attrs[item] = float(attrs[item])
        for item in ["valid_range", "valid_min", "valid_max"]:
            if item in attrs:
                attrs[item] = _cast_to_dtype(attrs[item], properties)
        if "flag_values" in attrs:
            items = attrs["flag_values"].split()
            attrs["flag_values"] = np.array(
                [_cast_to_dtype(item, properties) for item in items],
                properties["dtype"]) if len(items) != 1 else _cast_to_dtype(
                    float(attrs["flag_values"]), properties)
        # if "scale_factor" in attrs and "add_offset" not in attrs:
        #     attrs["add_offset"] = 0.0
        # if "add_offset" in attrs and "scale_factor" not in attrs:
        #     attrs["scale_factor"] = 1.0
        return encoding, xr.DataArray(data=data,
                                      dims=properties["shape"],
                                      name=name,
                                      attrs=attrs)

    def x_ac(self, x_ac: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return self._data_array("expert/cross_track_distance", x_ac)

    def lon(self, lon: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        # Longitude must be in [0, 360.0[
        return self._data_array("basic/longitude",
                                math.normalize_longitude(lon, 0))

    def lon_nadir(self, lon_nadir: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        # Longitude must be in [0, 360.0[
        return self._data_array("expert/longitude_nadir",
                                math.normalize_longitude(lon_nadir, 0))

    def lat(self, lat: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return self._data_array("basic/latitude", lat)

    def lat_nadir(self, lat_nadir: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return self._data_array("expert/latitude_nadir", lat_nadir)

    def ssh_karin(self, ssh: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return self._data_array("basic/ssh_karin", ssh)

    def ssh_karin_uncert(self, array: np.ndarray) -> None:
        self._data_array("ssh_karin_uncert", array)

    def ssh_nadir(self, array: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return {
            '_FillValue': 2147483647,
            'dtype': 'int32'
        }, xr.DataArray(data=array,
                        dims=self.variables["basic/time"]["shape"],
                        name="basic/ssh_nadir",
                        attrs={
                            'coordinates': 'longitude latitude',
                            'long_name': 'sea surface height',
                            'scale_factor': 0.0001,
                            'standard_name':
                            'sea surface height above reference ellipsoid',
                            'units': 'm',
                            'valid_min': np.int32(-15000000),
                            'valid_max': np.int32(150000000)
                        })

    def ssh_nadir_error(self, array: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return {
            '_FillValue': 2147483647,
            'dtype': 'int32'
        }, xr.DataArray(
            data=array,
            dims=self.variables["basic/time"]["shape"],
            name="error/ssh_nadir",
            attrs={
                'coordinates': 'longitude latitude',
                'long_name': 'sea surface height',
                'scale_factor': 0.0001,
                'standard_name':
                'sea surface height above reference ellipsoid',
                'units': 'm',
                'valid_min': np.int32(-15000000),
                'valid_max': np.int32(150000000),
                'comment':
                'Height of the sea surface free of measurement errors.'
            })

    def ssh_karin_error(self, array: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return {
            '_FillValue': 2147483647,
            'dtype': 'int32'
        }, xr.DataArray(
            data=array,
            dims=self.variables["basic/ssh_karin"]["shape"],
            name="error/ssh_karin",
            attrs={
                'coordinates': 'longitude latitude',
                'long_name': 'sea surface height',
                'scale_factor': 0.0001,
                'standard_name':
                'sea surface height above reference ellipsoid',
                'units': 'm',
                'valid_min': np.int32(-15000000),
                'valid_max': np.int32(150000000),
                'comment':
                'Height of the sea surface free of measurement errors.'
            })

    def baseline_dilation(self,
                          array: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return {
            '_FillValue': 2147483647,
            'dtype': 'int32'
        }, xr.DataArray(data=array,
                        dims=self.variables["basic/ssh_karin"]["shape"],
                        name="error/baseline_dilation",
                        attrs={
                            'long_name': 'Error due to baseline mast dilation',
                            'scale_factor': 0.0001,
                            '_FillValue': 2147483647,
                            'units': 'm',
                            'coordinates': 'longitude latitude'
                        })

    def roll(self, array: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return {
            '_FillValue': 2147483647,
            'dtype': 'int32'
        }, xr.DataArray(data=array,
                        dims=self.variables["basic/ssh_karin"]["shape"],
                        name="error/roll",
                        attrs={
                            'long_name': 'Error due to roll',
                            'units': 'm',
                            'scale_factor': 0.0001,
                            'coordinates': 'longitude latitude'
                        })

    def phase(self, array: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return {
            '_FillValue': 2147483647,
            'dtype': 'int32'
        }, xr.DataArray(data=array,
                        dims=self.variables["basic/ssh_karin"]["shape"],
                        name="error/phase",
                        attrs={
                            'long_name': 'Error due to phase',
                            'units': 'm',
                            'scale_factor': 0.0001,
                            'coordinates': 'longitude latitude'
                        })

    # def roll_phase_est(self, array: np.ndarray) -> Tuple[Dict, xr.DataArray]:
    #     return {
    #         '_FillValue': 2147483647,
    #         'dtype': 'int32'
    #     }, xr.DataArray(data=array,
    #                     dims=self.variables["basic/ssh_karin"]["shape"],
    #                     name="basic/roll_phase_est",
    #                     attrs={
    #                         'long_name': 'Error after estimaton of roll phase',
    #                         'units': 'm',
    #                         'scale_factor': 0.0001,
    #                         'valid_min': -15000000.0,
    #                         'valid_max': 150000000.0,
    #                         'coordinates': 'longitude latitude'
    #                     })

    def karin(self, array: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return {
            '_FillValue': 2147483647,
            'dtype': 'int32'
        }, xr.DataArray(data=array,
                        dims=self.variables["basic/ssh_karin"]["shape"],
                        name="error/karin",
                        attrs={
                            'long_name': 'KaRIn error',
                            'units': 'm',
                            'scale_factor': 0.0001,
                            'coordinates': 'longitude latitude'
                        })

    def timing(self, array: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return {
            '_FillValue': 2147483647,
            'dtype': 'int32'
        }, xr.DataArray(data=array,
                        dims=self.variables["basic/ssh_karin"]["shape"],
                        name="error/timing",
                        attrs={
                            'long_name': 'Timing error',
                            'units': 'm',
                            'scale_factor': 0.0001,
                            'coordinates': 'longitude latitude'
                        })

    def wet_troposphere(self, array: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return {
            '_FillValue': 2147483647,
            'dtype': 'int32',
            'scale_factor': 0.0001
        }, xr.DataArray(data=array,
                        dims=self.variables["basic/ssh_karin"]["shape"],
                        name="error/wet_troposphere",
                        attrs={
                            'long_name':
                            'Error due to wet troposphere path delay',
                            'units': 'm',
                            'scale_factor': 0.0001,
                            'coordinates': 'longitude latitude'
                        })

    def wet_troposphere_nadir(self,
                              array: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return {
            '_FillValue': 2147483647,
            'dtype': 'int32'
        }, xr.DataArray(data=array,
                        dims=self.variables["basic/time"]["shape"],
                        name="error/wet_troposphere_nadir",
                        attrs={
                            'long_name':
                            'Error due to wet troposphere path delay',
                            'units': 'm',
                            'scale_factor': 0.0001,
                            'coordinates': 'longitude latitude'
                        })

    def altimeter(self, array: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return {
            '_FillValue': 2147483647,
            'dtype': 'int32'
        }, xr.DataArray(data=array,
                        dims=self.variables["basic/time"]["shape"],
                        name="error/altimeter",
                        attrs={
                            'long_name': 'Altimeter error',
                            'standard_name': '',
                            'units': 'm',
                            'scale_factor': 0.0001,
                            'coordinates': 'longitude latitude'
                        })

    def fill_variables(self, variables,
                       shape) -> Iterator[Tuple[Dict, xr.DataArray]]:
        for item in self.variables:
            if item in variables:
                continue
            properties = self.variables[item]
            fill_value = self.fill_value(properties)
            data = np.full(tuple(shape[dim] for dim in properties["shape"]),
                           fill_value, properties["dtype"])
            yield self._data_array(item, data)


class Nadir:
    def __init__(self,
                 track: orbit_propagator.Pass,
                 standalone: Optional[bool] = True):
        self.standalone = standalone
        self.product_spec = ProductSpecification()
        self.num_lines = track.time.size
        self.encoding, self.data_vars = self.product_spec.time(track.time)
        self._data_array("lon_nadir",
                         track.lon_nadir)._data_array("lat_nadir",
                                                      track.lat_nadir)

    def _data_array(self, attr, data: np.ndarray):
        encoding, array = getattr(self.product_spec, attr)(data)
        if self.standalone:
            array.name = array.name.replace("_nadir", "")
        self.encoding[array.name] = encoding
        self.data_vars.append(array)
        return self

    def ssh(self, array: np.ndarray) -> None:
        self._data_array("ssh_nadir", array)

    def ssh_error(self, array: np.ndarray) -> None:
        self._data_array("ssh_nadir_error", array)

    def update_noise_errors(self, noise_errors: Dict[str, np.ndarray]) -> None:
        for k, v in noise_errors.items():
            if v.ndim == 1:
                self._data_array(k, v)

    def to_netcdf(self, cycle_number: int, pass_number: int, path: str,
                  complete_product: bool) -> None:
        LOGGER.info("write %s", path)
        data_vars = dict((item.name, item) for item in self.data_vars)
        # Variables that are not calculated are filled in in order to have a
        # product compatible with the PDD SWOT. Longitude is used as a
        # template.
        if complete_product and "basic/longitude" in data_vars:
            item = data_vars["basic/longitude"]
            shape = dict(zip(item.dims, item.shape))
            shape["num_sides"] = 2
            for encoding, array in self.product_spec.fill_variables(
                    data_vars.keys(), shape):
                self.encoding[array.name] = encoding
                self.data_vars.append(array)
        if "basic/longitude" in data_vars:
            lng = data_vars["basic/longitude"]
            lat = data_vars["basic/latitude"]
        elif self.standalone:
            lng = data_vars["expert/longitude"]
            lat = data_vars["expert/latitude"]
        else:
            lng = data_vars["expert/longitude_nadir"]
            lat = data_vars["expert/latitude_nadir"]

        # Variables must be written in the declaration order of the XML file
        unordered_vars = dict((item.name, item) for item in self.data_vars)
        data_vars = collections.OrderedDict()
        for item in self.product_spec.variables:
            if item in unordered_vars:
                data_vars[item] = unordered_vars[item]

        # Add variables that are not defined in the XML file
        for item in set(unordered_vars.keys()) - set(data_vars.keys()):
            data_vars[item] = unordered_vars[item]

        dataset = xr.Dataset(data_vars=data_vars,
                             attrs=global_attributes(
                                 self.product_spec.attributes, cycle_number,
                                 pass_number, self.data_vars[0].values, lng,
                                 lat))
        to_netcdf(dataset, path, self.encoding, mode="w")


class Swath(Nadir):
    def __init__(self, track: orbit_propagator.Pass) -> None:
        super().__init__(track, False)
        self.num_pixels = track.x_ac.size
        _x_ac = np.full((track.time.size, track.x_ac.size),
                        track.x_ac,
                        dtype=track.x_ac.dtype)
        self._data_array("x_ac", _x_ac)._data_array("lon",
                                                    track.lon)._data_array(
                                                        "lat", track.lat)

    def ssh(self, array: np.ndarray) -> None:
        self._data_array("ssh_karin", array)

    def ssh_error(self, array: np.ndarray) -> None:
        self._data_array("ssh_karin_error", array)

    def update_noise_errors(self, noise_errors: Dict[str, np.ndarray]) -> None:
        for k, v in noise_errors.items():
            if v.ndim == 2:
                self._data_array(k, v)
