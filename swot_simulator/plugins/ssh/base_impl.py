# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
import logging
import os
import re
from datetime import datetime
import time
from typing import Optional

import dask.array as da
import numba as nb
import numpy as np
import pyinterp
import pyinterp.backends.xarray
import xarray as xr

from .. import Interface

LOGGER = logging.getLogger(__name__)


class DatasetLoader:
    def load_dataset(
        self, first_date: np.datetime64, last_date: np.datetime64
    ) -> xr.Dataset:
        """
        Loads the data under the form of a xarray.Dataset. The loaded dataset should
        contain values that allow interpolating first_date and last_date. This means
        its time interval is a little large than [first_date, last_date].

        Moreover, the dataset should refer to the longitude, latitude, time and sea
        surface height using canonical names: lon, lat, time, ssh

        See also
        --------
        _shift_date

        :param first_date: The first date that needs to be interpolated
        :param last_date: The last date that needs to be interpolated
        :return: An xr.Dataset containing lon, lat, time and ssh
        """
        raise RuntimeError("Not implemented")

    @staticmethod
    def _shift_date(
        date: np.datetime64, shift: int, time_delta: np.datetime64
    ) -> np.datetime64:
        """
        Shift the input date using the time_delta of original data. This is
        useful to generate a time interval for which we need an original value.

        Example: if we have data on [t0, t1, dt], and we want an interpolation over
        [T0, T1], then we must make sure that t0 <= T0 - dt and t1 >= T1 + dt. If this
        condition is satisfied, interpolation at T0 and T1 will be possible. If this
        condition is not satisfied, interpolation becomes extrapolation.
        """
        if date.astype("int64") % time_delta.astype("int64") != 0:
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
                    f"grids: {set(frequency)} seconds"
                )
            return np.timedelta64(frequency[0], "s")
        except IndexError:
            raise RuntimeError("Check that your list of data is not empty")


@nb.njit(
    "(float32[::1])(int64[::1], float32[:, ::1], int64[::1])", cache=True, nogil=True
)
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

    def __init__(
        self,
        path: str,
        lon_name: str = "lon",
        lat_name: str = "lat",
        ssh_name: str = "ssh",
        time_name: str = "time",
        pattern: str = ".nc",
        date_fmt: Optional[str] = None,
    ):
        """
        Initialization of the netcdf loader.

        Example
        -------
        If we have netcdf files whose names are
        model_20120305_12h.nc, we must define the following to retrieve the time:
        pattern='model_P(?<date>\w+).nc'
        date_fmt='%Y%m%d_%Hh'

        :param path: Folder containing the netcdf files
        :param lon_name: longitude name in the netcdf files. Defaults to 'lon'
        :param lat_name: latitude name in the netcdf files. Defaults to 'lat'
        :param ssh_name: sea surface height name in the netcdf files. Defaults to 'ssh'
        :param time_name: time name in the netcdf files. Defaults to 'time'
        :param pattern: Pattern for the NetCDF file names. It should contain the
        P(?<date>) group to retrieve the time
        :param date_fmt: date formatter
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"{path!r}")

        self.lon_name = lon_name
        self.lat_name = lat_name
        self.ssh_name = ssh_name
        self.time_name = time_name

        self.regex = re.compile(pattern).search
        self.date_fmt = date_fmt
        self.ts = self._walk_netcdf(path)
        self.time_delta = self._calculate_time_delta(self.ts["date"])

    def _walk_netcdf(self, path: str) -> np.array:
        # Walks a netcdf folder and finds data files in it
        items = []
        length = -1

        for dir_path, _, filenames in os.walk(path):
            for filename in filenames:
                match = self.regex(filename)
                if match:
                    time_counter = np.datetime64(
                        datetime.strptime(match.group("date"), self.date_fmt)
                    )
                    filepath = os.path.join(dir_path, filename)
                    items.append((time_counter, filepath))
                    length = max(length, len(filepath))

        # The time series is encoded in a structured Numpy array containing
        # the date and path to the file.
        ts = np.array(
            items,
            dtype={
                "names": ("date", "path"),
                "formats": ("datetime64[s]", f"U{length}"),
            },
        )
        ts = ts[np.argsort(ts["date"])]
        return ts

    def select_netcdf_files(
        self, first_date: np.datetime64, last_date: np.datetime64
    ) -> np.array:
        first_date = self._shift_date(first_date, -1, self.time_delta)
        last_date = self._shift_date(last_date, 1, self.time_delta)
        if first_date < self.ts["date"][0] or last_date > self.ts["date"][-1]:
            raise IndexError(
                f"period [{first_date}, {last_date}] is out of range: "
                f"[{self.ts['date'][0]}, {self.ts['date'][-1]}]"
            )

        selected = np.logical_and(self.ts["date"] >= first_date, self.ts["date"] < last_date)
        return selected

    def load_dataset(self, first_date: np.datetime64, last_date: np.datetime64):
        LOGGER.debug(f"fetch data for {first_date}, {last_date}")
        selected = self.select_netcdf_files(first_date, last_date)
        ds = xr.open_mfdataset(
            self.ts["path"][selected], concat_dim=self.time_name, combine="nested"
        )

        if self.time_name not in ds.coords:
            LOGGER.debug(
                f"Time coordinate {self.time_name} was not found, assigning "
                f"axis with time from file names"
            )
            ds = ds.assign_coords({self.time_name: self.ts["dates"][selected]})

        LOGGER.debug(f"Renaming dataset with canonical names lon, lat, time, ssh")
        return ds.rename(
            {
                self.lon_name: "lon",
                self.lat_name: "lat",
                self.ssh_name: "ssh",
                self.time_name: "time",
            }
        )


class IrregularGridHandler(Interface):
    def __init__(self, dataset_loader: DatasetLoader):
        self.dataset_loader = dataset_loader

    def interpolate(self, lon, lat, times):
        """Interpolate the SSH for the given coordinates"""
        # Spatial interpolation of the SSH on the different selected grids.
        ds = self.dataset_loader.load_dataset(times.min(), times.max())
        dates_p = ds.time.load().data
        lon_p = ds.lon.load().data
        lat_p = ds.lat.load().data
        ssh_p = ds.ssh

        start_time = time.time()
        layers = []
        for index in range(len(ssh_p)):
            layers.append(self._spatial_interp(ssh_p[index, :], lon_p, lat_p, lon, lat))
        layers = np.stack(layers)
        LOGGER.debug(
            "interpolation completed in %.2fs for period %s, %s",
            time.time() - start_time,
            times.min(),
            times.max(),
        )

        # Time interpolation of the SSH.
        return _time_interp(
            dates_p.astype("int64"),
            layers,
            times.astype("datetime64[us]").astype("int64"),
        )

    @staticmethod
    def _spatial_interp(
        z_model: da.array,
        x_model: da.array,
        y_model: da.array,
        x_sat: np.ndarray,
        y_sat: np.ndarray,
    ) -> np.ndarray:
        mesh = pyinterp.RTree()
        mesh.packing(np.vstack((x_model, y_model)).T, z_model)

        z, _ = mesh.radial_basis_function(
            np.vstack((x_sat, y_sat)).T.astype("float32"),
            within=True,
            k=11,
            rbf="linear",
            num_threads=1,
        )
        return z.astype("float32")


class CartesianGridHandler(Interface):
    def __init__(self, dataset_loader: DatasetLoader):
        self.dataset_loader = dataset_loader

    def interpolate(
        self, lon: np.ndarray, lat: np.ndarray, dates: np.ndarray
    ) -> np.ndarray:
        """Interpolate the SSH to the required coordinates"""
        ds = self.dataset_loader.load_dataset(dates.min(), dates.max())

        interpolator = pyinterp.backends.xarray.Grid3D(ds.ssh)
        ssh = interpolator.trivariate(
            dict(longitude=lon, latitude=lat, time=dates.astype("datetime64[ns]")),
            bounds_error=True,
            interpolator="bilinear",
        )
        return ssh
