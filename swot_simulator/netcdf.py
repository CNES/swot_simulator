# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Write SWOT data to netCDF4 files.
=================================
"""
from typing import Any, Dict, Hashable, List, Optional, Tuple, Union
import collections
import copy
import datetime
import logging
import os
import pathlib

#
import netCDF4
import numpy as np
import xarray as xr

#
from . import __version__
from . import math
from . import orbit_propagator
from . import product_specification

#: Logger for this module
LOGGER = logging.getLogger(__name__)

#: The reference document of this simulator.
REFERENCE = "Gaultier, L., C. Ubelmann, and L.-L. Fu, 2016: The " \
    "Challenge of Using Future SWOT Data for Oceanic Field Reconstruction." \
    " J. Atmos. Oceanic Technol., 33, 119-126, doi:10.1175/jtech-d-15-0160" \
    ".1. http://dx.doi.org/10.1175/JTECH-D-15-0160.1."

#: Product specification file for the nadir product
NADIR_SPECIFICATION = pathlib.Path(__file__).parent / "nadir-gpr.xml"


def global_attributes(attributes: Dict[str, Dict[str, Any]], cycle_number: int,
                      pass_number: int, date: np.ndarray, lng: np.ndarray,
                      lat: np.ndarray, lon_eq: float,
                      time_eq: np.datetime64) -> Dict[str, Any]:
    """Calculates the global attributes of the pass"""
    def _iso_date(date: np.datetime64) -> str:
        """Return the time formatted according to ISO."""
        epoch = date.astype("datetime64[us]").astype(
            "int64") * 1e-6  # type: ignore
        return datetime.datetime.utcfromtimestamp(epoch).isoformat() + "Z"

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%SZ : Creation")

    ellipsoid_semi_major_axis = product_specification.cast_to_dtype(
        6378137, attributes["ellipsoid_semi_major_axis"])
    ellipsoid_flattening = product_specification.cast_to_dtype(
        1 / 298.25722356, attributes["ellipsoid_flattening"])
    result = collections.OrderedDict(
        Conventions="CF-1.7",
        title=attributes["title"]["attrs"]["description"],
        institution="CNES/JPL",
        source="Simulate product",
        history=now,
        platform="SWOT",
        product_version=__version__,
        references=REFERENCE,
        reference_document="D-56407_SWOT_Product_Description_L2_LR_SSH",
        contact="CNES aviso@altimetry.fr, JPL podaac@podaac.jpl.nasa.gov",
        cycle_number=product_specification.cast_to_dtype(
            cycle_number, attributes["cycle_number"]),
        pass_number=product_specification.cast_to_dtype(
            pass_number, attributes["pass_number"]),
        equator_time=_iso_date(time_eq),
        equator_longitude=math.normalize_longitude(lon_eq, 0),  # type: ignore
        time_coverage_start=_iso_date(date[0]),
        time_coverage_end=_iso_date(date[-1]),
        geospatial_lon_min=lng.min(),
        geospatial_lon_max=lng.max(),
        geospatial_lat_min=lat.min(),
        geospatial_lat_max=lat.max())
    if len(lng.shape) == 2:
        result.update([
            ("left_first_longitude", lng[0, 0]),
            ("left_first_latitude", lat[0, 0]),
            ("left_last_longitude", lng[-1, 0]),
            ("left_last_latitude", lat[-1, 0]),
            ("right_first_longitude", lng[0, -1]),
            ("right_first_latitude", lat[0, -1]),
            ("right_last_longitude", lng[-1, -1]),
            ("right_last_latitude", lat[-1, -1]),
        ])
    result.update([
        ("wavelength",
         product_specification.cast_to_dtype(0.008385803020979,
                                             attributes["wavelength"])),
        ("orbit_solution", "POE")
    ])
    for item in attributes:
        if item.startswith("xref_input"):
            result[item] = "N/A"
    result.update([("ellipsoid_semi_major_axis", ellipsoid_semi_major_axis),
                   ("ellipsoid_flattening", ellipsoid_flattening)])
    return result


def _split_group_name(name):
    """Get the name of the group and the variable deduced from the name."""
    if "/" in name:
        return tuple(name.rsplit("/", 1))
    return None, name


def _create_variable_args(
        encoding: Dict[Hashable, Dict], name: Hashable,
        variable: xr.DataArray) -> Tuple[Hashable, Dict[str, Any]]:
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


def _group_attributes(side: product_specification.Side) -> Dict[str, str]:
    """Gets the attributes of a group (unsmoothed products)."""
    if side == product_specification.Side.LEFT.value:
        return dict(description="Unsmoothed SSH measurement data and related "
                    "information for the left half swath.")
    return dict(description="Unsmoothed SSH measurement data and related "
                "information for the right half swath.")


def _create_variable(xr_dataset: xr.Dataset, nc_dataset: netCDF4.Dataset,
                     encoding: Dict[str, Dict[str, Dict[str, Any]]],
                     name: Hashable, unlimited_dims: Optional[List[str]],
                     variable: xr.DataArray) -> None:
    """Creation and writing of the NetCDF variable"""
    unlimited_dims = unlimited_dims or list()

    variable.attrs.pop("_FillValue", None)
    # Encode datetime64 to float64
    if np.issubdtype(variable.dtype, np.datetime64):
        # 946684800000000 number of microseconds between 2000-01-01 and
        # 1970-01-01
        values = (variable.values.astype("datetime64[us]").astype("int64") -
                  946684800000000) * 1e-6
        if variable.name in xr_dataset.coords:
            xr_dataset = xr_dataset.assign_coords(
                coords={variable.name: values})
            attrs = variable.attrs
            variable = xr_dataset[variable.name]
            variable.attrs.update(attrs)
        else:
            variable.values = values
        assert (
            variable.attrs["units"] == "seconds since 2000-01-01 00:00:00.0")
    dtype, kwargs = _create_variable_args(encoding, name, variable)

    group, name = _split_group_name(name)

    if group is not None:
        if group not in nc_dataset.groups:
            nc_dataset = nc_dataset.createGroup(group)
            if group in ["left", "right"]:
                nc_dataset.setncatts(
                    _group_attributes(
                        getattr(product_specification.Side,
                                group.upper()).value))
        else:
            nc_dataset = nc_dataset.groups[group]

    # If the dimensions doesn't exist then we have to create them.
    if not nc_dataset.dimensions:
        for dim_name, size in xr_dataset.dims.items():
            dim_group, dim_name = _split_group_name(dim_name)
            if dim_group == group:
                nc_dataset.createDimension(
                    dim_name, None if dim_name in unlimited_dims else size)

    ncvar = nc_dataset.createVariable(
        name, dtype,
        tuple(_split_group_name(item)[-1] for item in variable.dims), **kwargs)
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


def _write_nadir_product(ds: xr.Dataset, path: str,
                         complete_product: bool) -> None:
    """Write the nadir product to a netCDF file."""
    variables, _attributes = product_specification.parse_specification_file(
        str(NADIR_SPECIFICATION))
    variables["data_01/ku/simulated_error_altimeter"] = dict(
        attrs=dict(_FillValue=2147483647,
                   long_name='Altimeter error',
                   standard_name='',
                   units='m',
                   scale_factor=0.0001,
                   coordinates='/data_01/longitude /data_01/latitude'),
        dtype='int32',
        shape=("data_01/time", ))

    ds = ds.rename_dims({"num_lines": "data_01/time"})
    ds = ds.rename_vars(
        dict(time="data_01/time",
             latitude_nadir="data_01/latitude",
             longitude_nadir="data_01/longitude",
             simulated_error_altimeter="data_01/ku/simulated_error_altimeter"))
    # Rename the simulated variables to match the specification pattern (These
    # variables do not exist in the official product. Only the SSHA is present.)
    for name in ["simulated_true_ssh_nadir", "ssh_nadir", "swh_nadir"]:
        if name in ds.variables:
            group_path = f"data_01/ku/{name.replace('_nadir', '')}"
            ds = ds.rename_vars({name: group_path})
            attrs = ds.variables[group_path].attrs.copy()
            attrs["coordinates"] = "/data_01/longitude /data_01/latitude"
            attrs["_FillValue"] = 2147483647
            variables[group_path] = dict(attrs=attrs,
                                         dtype='int32',
                                         shape=("data_01/time", ))
    data_vars = {}
    encoding = {}
    for name in variables:
        if name not in ds.data_vars:
            continue
        encoding_array, array = product_specification.build_array(
            name, variables, ds.variables[name].values)
        data_vars[array.name] = array
        encoding.update(encoding_array)

    shape = ds.variables["data_01/time"].shape

    if complete_product:
        for item, properties in variables.items():
            if item in data_vars:
                continue
            data = np.full(shape,
                           product_specification.encode_fill_value(properties),
                           properties["dtype"])
            encoding_array, array = product_specification.build_array(
                item, variables, data)
            data_vars[array.name] = array
            encoding.update(encoding_array)

    # import pdb; pdb.set_trace()
    ds.attrs["title"] = "GDR - Reduced dataset"
    ds = xr.Dataset(data_vars=data_vars, attrs=ds.attrs)
    to_netcdf(ds, path, encoding=encoding, mode="w")


class Nadir:
    """Handle the nadir measurements dataset.

    Args:
        track (orbit_propagator.Pass): Properties of the half-orbit to write
    """
    def __init__(self,
                 track: orbit_propagator.Pass,
                 product_type: Optional[str] = None):
        self.product_spec = product_specification.ProductSpecification(
            product_type)
        self.num_lines = track.time.size
        self.num_pixels = 1
        self.encoding, self.data_vars = self.product_spec.time(track.time)
        self._data_array("lon_nadir", track.lon_nadir)
        self._data_array("lat_nadir", track.lat_nadir)
        self._time_eq = track.time_eq
        self._lon_eq = track.lon_eq
        self.ascending = np.all(np.diff(track.lat_nadir) > 0)

    def _data_array(self,
                    attr: str,
                    data: np.ndarray,
                    fill_value: Optional[np.ndarray] = None) -> None:
        """Get the properties of a data array to be inserted in the dataset.

        Args:
            attr (dict): Name of the method of the `ProductSpec` class defining
                the variable to be inserted.
            data (np.array): data to be recorded
            fill_value (np.array, optional): Fill value of the reference
                track vector to be inserted in the centre of the swath.
        """
        # Is it necessary to insert a central pixel?
        if len(data.shape) == 2 and self.num_pixels % 2 == 1:
            # If the data is a grid, then the defined fill values are inserted
            # or a Nan vector if the fill values are undefined.
            if fill_value is None:
                fill_value = np.full((data.shape[0], ),
                                     np.nan,
                                     dtype=np.float64)
            middle = data.shape[1] // 2
            data = np.c_[data[:, :middle], fill_value[:, np.newaxis],
                         data[:, middle:]]
        variable = getattr(self.product_spec, attr)(data)
        if variable is not None:
            encoding, array = variable
            self.encoding.update(encoding)
            self.data_vars += array

    def ssh(self, array: np.ndarray) -> None:
        """Sets the variable describing the SSH to nadir.

        Args:
            array (np.ndarray): Data to be recorded
        """
        self._data_array("ssh_nadir", array)

    def swh(self, array: np.ndarray) -> None:
        """Sets the variable describing the SWH to nadir.

        Args:
            array (np.ndarray): Data to be recorded
        """
        self._data_array("swh_nadir", array)

    def simulated_true_ssh(self, array: np.ndarray) -> None:
        """Sets the variable describing the SSH to nadir free of measurement
        errors.

        Args:
            array (np.ndarray): Data to be recorded
        """
        self._data_array("simulated_true_ssh_nadir", array)

    def update_noise_errors(self, noise_errors: Dict[str, np.ndarray]) -> None:
        """Sets the values of the simulated errors.

        Args:
            noise_errors (dict): Simulated errors to be recorded
        """
        for k, v in noise_errors.items():
            if v.ndim == 1:
                self._data_array(k, v)

    def to_xarray(self, cycle_number: int, pass_number: int,
                  complete_product: bool) -> xr.Dataset:
        """Converts this instance into a xarray dataset.

        Args:
            cycle_number (int): Cycle number.
            pass_number (int): Pass number.
            complete_product (bool): True if you want to obtain a complete
                SWOT dataset, i.e. containing all the variables of the
                official dataset, even those not calculated by the simulator.

        Returns:
            xr.Dataset: xarray dataset
        """
        data_vars = dict((item.name, item) for item in self.data_vars)

        if "longitude" in data_vars:
            lng = data_vars["longitude"]
            lat = data_vars["latitude"]
        elif f"{product_specification.Side.LEFT.value}/longitude" in data_vars:
            lng = np.c_[
                data_vars[f"{product_specification.Side.LEFT.value}/longitude"],
                data_vars[
                    f"{product_specification.Side.RIGHT.value}/longitude"]]
            lat = np.c_[
                data_vars[f"{product_specification.Side.LEFT.value}/latitude"],
                data_vars[f"{product_specification.Side.RIGHT.value}/latitude"]]
        else:
            lng = data_vars["longitude_nadir"]
            lat = data_vars["latitude_nadir"]

        # Variables that are not calculated are filled in in order to have a
        # product compatible with the PDD SWOT. Longitude is used as a
        # template.
        pdd_vars = copy.copy(self.data_vars)
        if complete_product and len(lng.shape) == 2:

            # Creation of an associative dictionary between the name of the
            # dimensions of the dataset and their values.
            dimensions = []
            {
                dimensions.extend(list(zip(item.dims, item.shape)))
                for item in data_vars.values()
            }
            shape = dict(dimensions)

            # The simulator does not handle variables using the "num_sides"
            # dimension.
            shape.update(self.product_spec.num_sides())
            for encoding, array in self.product_spec.fill_variables(
                    data_vars.keys(), shape):
                self.encoding.update(encoding)
                pdd_vars += array

        # Variables must be written in the declaration order of the XML file
        unordered_vars = dict((item.name, item) for item in pdd_vars)
        data_vars = collections.OrderedDict()
        for item in self.product_spec.variables:
            if item in unordered_vars:
                data_vars[item] = unordered_vars[item]

        # Add variables that are not defined in the XML file
        for item in set(unordered_vars.keys()) - set(data_vars.keys()):
            data_vars[item] = unordered_vars[item]

        return xr.Dataset(data_vars=data_vars,
                          attrs=global_attributes(self.product_spec.attributes,
                                                  cycle_number, pass_number,
                                                  self.data_vars[0].values,
                                                  np.asarray(lng),
                                                  np.asarray(lat),
                                                  self._lon_eq, self._time_eq))

    def _write_netcdf(self, path: str, cycle_number: int, pass_number: int,
                      complete_product: bool) -> None:
        """Writes the Nadir dataset to a netCDF file.

        This method must be specialized to write Nadir and Karin products.
        Indeed, these two products have different structures and cannot be
        processed in the same way.
        """
        dataset = self.to_xarray(cycle_number, pass_number, False)
        _write_nadir_product(dataset, path, complete_product)

    def to_netcdf(self, cycle_number: int, pass_number: int, path: str,
                  complete_product: bool) -> None:
        """Writes the dataset in a netCDF file.

        Args:
            cycle_number (int): Cycle number.
            pass_number (int): Pass number.
            complete_product (bool): True if you want to obtain a complete
                SWOT dataset, i.e. containing all the variables of the
                official dataset, even those not calculated by the simulator.
        """
        LOGGER.info("write %s", path)
        try:
            self._write_netcdf(path, cycle_number, pass_number,
                               complete_product)
        except:
            if os.path.exists(path):
                os.unlink(path)
            raise


class Swath(Nadir):
    """Handle the KaRIn measurements dataset.

    Args:
        track (orbit_propagator.Pass): Properties of the half-orbit to write
        central_pixel (bool): Inserts or not in the swath, a central pixel
            divided in two by the reference ground track.
    """
    def __init__(self,
                 track: orbit_propagator.Pass,
                 central_pixel: bool = False,
                 product_type: Optional[str] = None) -> None:
        super().__init__(track, product_type)
        self.num_pixels = track.x_ac.size + int(central_pixel)
        _x_ac = np.full((track.time.size, track.x_ac.size),
                        track.x_ac,
                        dtype=track.x_ac.dtype)
        self._data_array("x_ac", _x_ac,
                         np.full((track.time.size, ), 0, dtype=np.float64))
        self._data_array("lon", track.lon, track.lon_nadir)
        self._data_array("lat", track.lat, track.lat_nadir)

    def ssh(self, array: np.ndarray) -> None:
        """Sets the variable describing the KaRIn SSH.

        Args:
            array (np.ndarray): Data to be recorded
        """
        self._data_array("ssh_karin", array)

    def swh(self, array: np.ndarray) -> None:
        """Sets the variable describing the KaRIn SWH.

        Args:
            array (np.ndarray): Data to be recorded
        """
        if not self.product_spec.unsmoothed:
            self._data_array("swh_karin", array)

    def simulated_true_ssh(self, array: np.ndarray) -> None:
        """Sets the variable describing the KaRIn SSH free of measurement
        errors.

        Args:
            array (np.ndarray): Data to be recorded
        """
        self._data_array("simulated_true_ssh_karin", array)

    def update_noise_errors(self, noise_errors: Dict[str, np.ndarray]) -> None:
        """Sets the values of the simulated errors.

        Args:
            noise_errors (dict): Simulated errors to be recorded
        """
        if self.product_spec.ssh_karin_name() in self.product_spec:
            for k, v in noise_errors.items():
                if v.ndim == 2:
                    self._data_array(k, v)

    def _write_netcdf(self, path: str, cycle_number: int, pass_number: int,
                      complete_product: bool) -> None:
        """Writes the swath dataset to a netCDF file."""
        dataset = self.to_xarray(cycle_number, pass_number, complete_product)
        to_netcdf(dataset, path, self.encoding, mode="w")
