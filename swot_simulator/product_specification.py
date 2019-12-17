# Copyright (c) 2019 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Parse/Load the product specification
====================================
"""
from typing import Any, Dict, Optional, Tuple
import datetime
import copy
import logging
import os
import xml.etree.ElementTree as xt
import numpy as np
import xarray as xr
from . import orbit_propagator
from . import math

LOGGER = logging.getLogger(__name__)

REFERENCE = "Gaultier, L., C. Ubelmann, and L.-L. Fu, 2016: The " \
    "Challenge of Using Future SWOT Data for Oceanic Field Reconstruction." \
    " J. Atmos. Oceanic Technol., 33, 119-126, doi:10.1175/jtech-d-15-0160" \
    ".1. http://dx.doi.org/10.1175/JTECH-D-15-0160.1."


def _find(element: xt.Element, tag: str) -> xt.Element:
    result = element.find(tag)
    if result is None:
        raise RuntimeError("The XML tag '" + tag + "' doesn't exist")
    return result


def _parse_type(dtype, width, signed):
    if dtype == "real":
        return getattr(np, "float" + width)
    elif dtype == "integer":
        return getattr(np, "u" if not signed else "" + "int" + width)
    elif dtype == "string":
        return np.str
    elif dtype == "char":
        assert int(width) == 1
        return np.int8
    raise ValueError("Data type '" + dtype + "' is not recognized.")


def global_attributes(attributes: Dict[str, Dict[str, str]], cycle_number: int,
                      pass_number: int, date: np.ndarray) -> Dict[str, Any]:
    def _encode(attr_value: str, properties: Dict[str, str]):
        return getattr(np, properties["dtype"])(attr_value)

    def _iso_date(date: np.datetime64) -> str:
        return datetime.datetime.utcfromtimestamp(
            date.astype("datetime64[us]").astype("int64") * 1e-6).isoformat()

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

    return {
        # "contact": "TODO",
        "Conventions": "CF-1.6",
        "cycle_number": _encode(cycle_number, attributes["cycle_number"]),
        "pass_number": _encode(pass_number, attributes["pass_number"]),
        "institution": "CNES/JPL",
        "mission_name": "SWOT",
        "references": REFERENCE,
        "standard_name_vocabulary": "CF Standard Name Table vNN",
        "source": "Simulate product",
        "time_coverage_start": _iso_date(date[0]),
        "time_coverage_end": _iso_date(date[-1]),
        "time_coverage_duration": _iso_duration(date[-1] - date[0]),
        "time_coverage_resolution": "P1S",
        "history": now
    }


def _parser(tree: xt.ElementTree):
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
            bool(item.attrib["signed"]) if "signed" in item.attrib else None)
        annotation = item.find("annotation")
        if annotation is None:
            continue
        varname = item.attrib["name"]
        if varname.startswith("/@"):
            attributes[varname[2:]] = dict(attrs=annotation.attrib,
                                           dtype=dtype.__name__)
        else:
            del annotation.attrib["app"]
            variables[varname[1:]] = dict(attrs=annotation.attrib,
                                          dtype=dtype.__name__,
                                          shape=shapes[item.attrib["shape"]])

    return variables, attributes


class ProductSpecification:
    """Parse and load into memory the product specification"""
    SPECIFICATION = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "l2b-expert.xml")

    def __init__(self):
        self.variables, self.attributes = _parser(xt.parse(self.SPECIFICATION))

    def time(self, time: np.ndarray) -> xr.DataArray:
        return xr.DataArray(data=time,
                            dims=self.variables["basic/time"]["shape"],
                            name="time",
                            attrs=dict(long_name="time in UTC",
                                       standard_name="time",
                                       axis="T"))

    def _data_array(self,
                    name: str,
                    data: np.ndarray,
                    axis: Optional[bool] = False) -> Tuple[Dict, xr.DataArray]:
        properties = self.variables[name]
        attrs = copy.deepcopy(properties["attrs"])

        # The fill value is casted to the target value of the variable
        if axis is not None:
            fill_value = getattr(np, properties["dtype"])(attrs["_FillValue"])
        else:
            fill_value = None
        del attrs["_FillValue"]

        # Reading the storage properties of the variable ()
        encoding = dict(_FillValue=fill_value, dtype=properties["dtype"])
        for item in ["add_offset", "scale_factor"]:
            if item in attrs:
                encoding[item] = float(attrs[item])
                del attrs[item]

        # Some values read from the XML files must be decoded
        for item in ["valid_min", "valid_max"]:
            if item in attrs:
                attrs[item] = float(attrs[item])
        return encoding, xr.DataArray(data=data,
                                      dims=properties["shape"],
                                      name=name.split("/")[-1],
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

    def ssh_nadir(self, ssh: np.ndarray) -> Tuple[Dict, xr.DataArray]:
        return {
            '_FillValue': 2147483647,
            'dtype': 'int32',
            'scale_factor': 0.0001
        }, xr.DataArray(data=ssh,
                        dims=self.variables["basic/time"]["shape"],
                        name="ssh_nadir",
                        attrs={
                            'long_name': 'sea surface height',
                            'standard_name':
                            'sea surface height above reference ellipsoid',
                            'units': 'm',
                            'valid_min': -15000000.0,
                            'valid_max': 150000000.0,
                            'coordinates': 'longitude latitude'
                        })


class Nadir:
    def __init__(self,
                 track: orbit_propagator.Pass,
                 standalone: Optional[bool] = True):
        self.standalone = standalone
        self.product_spec = ProductSpecification()
        self.num_lines = track.time.size
        self.data_vars = [
            self.product_spec.time(track.time),
        ]
        self.encoding = dict(
            (item.name, {
                "_FillValue": None,
                "units": "microseconds since 2000-01-01 00:00:00+00:00"
            }) for item in self.data_vars)
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

    def to_netcdf(self, cycle_number: int, pass_number: int,
                  path: str) -> None:
        dataset = xr.Dataset(data_vars=dict(
            (item.name, item) for item in self.data_vars),
                             attrs=global_attributes(
                                 self.product_spec.attributes, cycle_number,
                                 pass_number, self.data_vars[0].values))
        LOGGER.info("write %s", path)
        dataset.to_netcdf(path, encoding=self.encoding)


class Swath(Nadir):
    def __init__(self, track: orbit_propagator.Pass) -> None:
        super().__init__(track, False)
        self.num_pixels = track.x_ac.size
        self._data_array(
            "x_ac",
            np.full((track.time.size, track.x_ac.size),
                    track.x_ac,
                    dtype=track.x_ac.dtype))._data_array(
                        "lon", track.lon)._data_array("lat", track.lat)

    def ssh(self, array: np.ndarray) -> None:
        self._data_array("ssh_karin", array)
