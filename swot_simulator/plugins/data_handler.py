# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Data handler
============

Set of classes that helps developping plugins for ssh and swh. This includes a default
netcdf reader for loading a dataset, and default interpolators for both regular and
irregular grids.
"""
import abc
import datetime
import logging
import os
import re
import time

import dask.array as da
import numba as nb
import numpy as np
import pyinterp
import pyinterp.backends.xarray
import xarray as xr

from . import Interface

LOGGER = logging.getLogger(__name__)


class DatasetLoader:
    """
    Interface that specializes in loading data for the plugin. This is helpful to
    separate the data loading and interpolator definition. The dataset loader has only
    one job : given a time range for interpolation, loads the model dataset that is
    needed to do the said interpolation and transform it to have canonical variable
    names.
    """
    @abc.abstractmethod
    def load_dataset(self, first_date: np.datetime64,
                     last_date: np.datetime64) -> xr.Dataset:
        """
        Loads the data under the form of a xr.Dataset. The loaded dataset should
        contain values that allow interpolating first_date and last_date. This means
        its time interval is a little large than [first_date, last_date].

        Moreover, the dataset should refer to the longitude, latitude, time and sea
        surface height using canonical names: lon, lat, time, ssh

        Args:
            first_date (np.datetime64): The first date that needs to be interpolated
            last_date (np.datetime64): The last date that needs to be interpolated

        Returns:
            xr.Dataset: dataset containing lon, lat, time and ssh variables, with
            canonical names.

        See also:
            _shift_date
        """
        ...

    @staticmethod
    def _shift_date(date: np.datetime64, shift: int,
                    time_delta: np.timedelta64) -> np.datetime64:
        """
        Shift the input date using the time_delta of original data. This is
        useful to generate a time interval for which we need an original value.

        Args:
            date (np.datetime64): interpolation date
            shift (int): 1 for a later date, -1 for an earlier one
            time_delta (np.timedelta64): delta specifying the time resolution of the
            model data.

        Returns:
            The input date if it is the input date is a multiple of time_delta (meaning
            it is on the model time axis). Else, the output is shifted.

        Example:
            If we have data on [t0, t1, dt], and we want an interpolation over [T0, T1],
            then we must make sure that t0 <= T0 - dt and t1 >= T1 + dt. If this
            condition is satisfied, interpolation at T0 and T1 will be possible. If this
            condition is not satisfied, interpolation becomes extrapolation.
        """
        # Before comparing the date and timedelta, ensure they have the same unit
        date_same_unit = date + time_delta - time_delta
        time_delta_same_unit = date + time_delta - date

        if date_same_unit.astype("int64") % time_delta_same_unit.astype(
                "int64") != 0:  # type: ignore
            return date + time_delta * shift
        return date

    @staticmethod
    def _calculate_time_delta(dates: xr.DataArray) -> np.timedelta64:
        """Calculation of the delta T between two consecutive grids."""
        frequency = np.diff(dates)
        try:
            if not np.all(frequency == frequency[0]):
                raise RuntimeError(
                    "Time series does not have a constant step between two "
                    f"grids: {set(frequency)} seconds")
            return np.timedelta64(frequency[0], "ns")
        except IndexError as exc:
            raise RuntimeError(
                "Check that your list of data is not empty") from exc


@nb.njit("(float32[::1])(int64[::1], float32[:, ::1], int64[::1])",
         cache=True,
         nogil=True)  # type:ignore
def _time_interp(xp: np.ndarray, yp: np.ndarray, xi: np.ndarray) -> np.ndarray:
    """Time interpolation for the different spatial grids interpolated on the
    SWOT data"""
    xp_diff = np.diff(xp)

    assert xp.shape[0] == yp.shape[0] and yp.shape[1] == xi.shape[0]
    assert np.all(xp_diff == xp_diff[0])

    result = np.empty(yp.shape[1:], dtype=yp.dtype)

    step = 1.0 / xp_diff[0]
    size = xp.size

    for ix in range(yp.shape[0]):
        index = int(np.around((xi[ix] - xp[0]) * step))
        assert index >= 0 or index <= size
        if index == size - 1:
            i0 = index - 1
            i1 = index
        else:
            i0 = index
            i1 = index + 1
        t0 = xp[i0]
        t1 = xp[i1]
        dt = t1 - t0
        w0 = (t1 - xi[ix]) / dt
        w1 = (xi[ix] - t0) / dt

        for jx in range(yp.shape[1]):
            result[jx] = (w0 * yp[i0, jx] + w1 * yp[i1, jx]) / (w0 + w1)

    return result


class NetcdfLoader(DatasetLoader):
    """
    Plugin that implements a netcdf reader. The netcdf reader works on files whose
    names have the date in it. A pattern (ex. P(?<date>.*).nc), associated with a
    date formatter (ex. %Y%m%d) is used to get build the time series.

    Netcdf files can be expensive to concatenate if there are a lot of files. This
    loader avoid loading too much files by building a dictionary matching file paths
    to their time. During the interpolation, where only a given time period is
    needed, only the files that cover the time period are loaded in the dataset.
    """
    def __init__(self,
                 path: str,
                 date_fmt: str,
                 lon_name: str = "lon",
                 lat_name: str = "lat",
                 ssh_name: str = "ssh",
                 time_name: str = "time",
                 pattern: str = ".nc"):
        """
        Initialization of the netcdf loader.

        Args:
            path (str): Folder containing the netcdf files
            date_fmt (str): date formatter
            lon_name (str): longitude name in the netcdf files. Defaults to
                'lon'
            lat_name (str): latitude name in the netcdf files. Defaults to 'lat'
            ssh_name (str): sea surface height name in the netcdf files.
                Defaults to 'ssh'
            time_name (str): time name in the netcdf files. Defaults to 'time'
            pattern (str): Pattern for the NetCDF file names. It should contain
                the P(?<date>) group to retrieve the time

        Example:
            If we have netcdf files whose names are model_20120305_12h.nc, we must
            define the following to retrieve the time::

                loader = NetcdfLoader(
                    '.',
                    pattern='model_P(?<date>\\w+).nc',
                    date_fmt='%Y%m%d_%Hh'
                )
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"{path!r}")

        self.lon_name = lon_name
        self.lat_name = lat_name
        self.ssh_name = ssh_name
        self.time_name = time_name

        self.regex = re.compile(pattern).search
        self.date_fmt = date_fmt
        self.time_series = self._walk_netcdf(path)
        self.time_delta = self._calculate_time_delta(self.time_series["date"])

    def _walk_netcdf(self, path: str) -> np.ndarray:
        # Walks a netcdf folder and finds data files in it
        items = []
        length = -1

        for dir_path, _, filenames in os.walk(path):
            for filename in filenames:
                match = self.regex(filename)
                if match:
                    time_counter = np.datetime64(
                        datetime.datetime.strptime(match.group("date"),
                                                   self.date_fmt))
                    filepath = os.path.join(dir_path, filename)
                    items.append((time_counter, filepath))
                    length = max(length, len(filepath))

        # The time series is encoded in a structured Numpy array containing
        # the date and path to the file.
        time_series = np.array(
            items,
            dtype={
                "names": ("date", "path"),
                "formats": ("datetime64[s]", f"U{length}"),
            },
        )
        time_series = time_series[np.argsort(time_series["date"])]
        return time_series

    def select_netcdf_files(self, first_date: np.datetime64,
                            last_date: np.datetime64) -> np.ndarray:
        first_date = self._shift_date(first_date.astype("datetime64[ns]"), -1,
                                      self.time_delta)
        last_date = self._shift_date(last_date.astype("datetime64[ns]"), 1,
                                     self.time_delta)
        if first_date < self.time_series["date"][
                0] or last_date > self.time_series["date"][-1]:
            raise IndexError(
                f"period [{first_date}, {last_date}] is out of range: "
                f"[{self.time_series['date'][0]}, {self.time_series['date'][-1]}]"
            )

        selected = np.logical_and(self.time_series["date"] >= first_date,
                                  self.time_series["date"] < last_date)
        return selected

    def load_dataset(self, first_date: np.datetime64,
                     last_date: np.datetime64):
        LOGGER.debug("fetch data for %s, %s", first_date, last_date)
        selected = self.select_netcdf_files(first_date, last_date)
        dataset = xr.open_mfdataset(self.time_series["path"][selected],
                                    concat_dim=self.time_name,
                                    combine="nested")

        if self.time_name not in dataset.coords:
            LOGGER.debug(
                "Time coordinate %s was not found, assigning "
                "axis with time from file names", self.time_name)
            dataset = dataset.assign_coords(
                {self.time_name: self.time_series["dates"][selected]})

        LOGGER.debug(
            "Renaming dataset with canonical names lon, lat, time, ssh")
        return dataset.rename({
            self.lon_name: "lon",
            self.lat_name: "lat",
            self.ssh_name: "ssh",
            self.time_name: "time"
        })


class IrregularGridHandler(Interface):
    """
    Default interpolator for an irregular grid. First, uses an RTree to do the spatial
    interpolation of all model grid, then do the time interpolation with a simple
    weighting of two grid
    """
    def __init__(self, dataset_loader: DatasetLoader):
        self.dataset_loader = dataset_loader

    def interpolate(self, lon: np.ndarray, lat: np.ndarray, dates: np.ndarray):
        """Interpolate the SSH for the given coordinates"""
        # Spatial interpolation of the SSH on the different selected grids.
        dataset = self.dataset_loader.load_dataset(
            dates.min(),  # type: ignore
            dates.max())  # type: ignore
        dates_p = dataset.time.load().data
        lon_p = dataset.lon.data
        lat_p = dataset.lat.data
        ssh_p = dataset.ssh.data

        start_time = time.time()
        layers = []
        for index in range(len(ssh_p)):
            layers.append(
                self._spatial_interp(ssh_p[index, :], lon_p, lat_p, lon, lat))
        layers = np.stack(layers)
        LOGGER.debug("interpolation completed in %.2fs for period %s, %s",
                     time.time() - start_time, dates.min(), dates.max())

        # Time interpolation of the SSH.
        return _time_interp(
            dates_p.astype("int64"),
            layers,
            dates.astype("datetime64[us]").astype("int64"),
        )

    @staticmethod
    def _spatial_interp(z_model: da.Array, x_model: da.Array,
                        y_model: da.Array, x_sat: np.ndarray,
                        y_sat: np.ndarray) -> np.ndarray:
        mesh = pyinterp.RTree()
        mesh.packing(
            np.vstack((x_model.compute(), y_model.compute())).T,
            z_model.compute())

        z, _ = mesh.radial_basis_function(
            np.vstack((x_sat, y_sat)).T.astype("float32"),
            within=True,
            k=11,
            rbf="thin_plate",
            num_threads=1,
        )
        return z.astype("float32")


class CartesianGridHandler(Interface):
    """
    Default interpolator for regular grid.
    Uses pyinterp.backends.xarray.Grid3D.trivariate interpolator
    """
    def __init__(self, dataset_loader: DatasetLoader):
        self.dataset_loader = dataset_loader

    def interpolate(self, lon: np.ndarray, lat: np.ndarray,
                    dates: np.ndarray) -> np.ndarray:
        """Interpolate the SSH to the required coordinates"""
        dataset = self.dataset_loader.load_dataset(
            dates.min(),  # type: ignore
            dates.max())  # type: ignore

        interpolator = pyinterp.backends.xarray.Grid3D(dataset.ssh)
        ssh = interpolator.trivariate(dict(longitude=lon,
                                           latitude=lat,
                                           time=dates),
                                      bounds_error=True,
                                      interpolator="bilinear")
        return ssh
