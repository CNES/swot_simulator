# Copyright (c) 2020 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Interpolation of the SSH HYCOM
==============================
"""
import os
import re
import numpy as np
import pyinterp
import pyinterp.backends.xarray
import xarray as xr

from . import detail


class HYCOM(detail.CartesianGridHandler):
    """
    Interpolation of the SSH HYCOM
    """

    #: Decode the product date encoded from the file name.
    #hycom_GLBu0.08_191_2012031900_t012.nc
    PATTERN = re.compile(
        r"hycom_GLBu0.08_191_(\d{4})(\d{2})(\d{2})(\d{2})_t\d{3}.nc").search

    def load_ts(self):
        """Loading in memory the time axis of the time series"""
        items = []
        length = -1

        for root, _, files in os.walk(self.path):
            for item in files:
                match = self.PATTERN(item)
                if match is not None or True:
                    filename = os.path.join(root, item)
                    items.append(
                        (np.datetime64(f"{match.group(1)}-{match.group(2)}-"
                                       f"{match.group(3)}"), filename))
                    length = max(length, len(filename))

        # The time series is encoded in a structured Numpy array containing
        # the date and path to the file.
        ts = np.array(items,
                      dtype={
                          'names': ('date', 'path'),
                          'formats': ('datetime64[s]', f'U{length}')
                      })
        self.ts = ts[np.argsort(ts["date"])]

        # The frequency between the grids must be constant.
        frequency = set(np.diff(self.ts["date"].astype(np.int64)))
        if len(frequency) > 1:
            raise RuntimeError(
                "Time series does not have a constant step between two "
                f"grids: {frequency} seconds")
        elif len(frequency) != 1:
            raise RuntimeError("Check that your list of data is not empty")
        # The frequency is stored in order to load the grids required to
        # interpolate the SSH.
        self.dt = np.timedelta64(frequency.pop(), 's')

    def load_dataset(
            self, first_date: np.datetime64,
            last_date: np.datetime64) -> pyinterp.backends.xarray.Grid3D:
        """Loads the 3D cube describing the SSH in time and space."""
        if first_date < self.ts["date"][0] or last_date > self.ts["date"][-1]:
            raise IndexError(
                f"period [{first_date}, {last_date}] is out of range: "
                f"[{self.ts['date'][0]}, {self.ts['date'][-1]}]")
        first_date -= self.dt
        last_date += self.dt

        selected = self.ts["path"][:]

        ds = xr.open_mfdataset(selected,
                               concat_dim="time",
                               combine="nested",
                               decode_times=False)
        ds = ds.sel(depth=0)

        x_axis = pyinterp.Axis(ds.variables["lon"][:], is_circle=True)
        y_axis = pyinterp.Axis(ds.variables["lat"][:])
        hours = (ds.variables['time'][:].data *
                 3600000000).astype('timedelta64[us]')
        time = np.datetime64('2000') + hours
        z_axis = pyinterp.TemporalAxis(time)
        var = ds.surf_el[:].T
        return pyinterp.Grid3D(x_axis, y_axis, z_axis, var)

    def interpolate(self, lon: np.ndarray, lat: np.ndarray,
                    time: np.ndarray) -> np.ndarray:
        """Interpolate the SSH to the required coordinates"""
        interpolator = self.load_dataset(time.min(), time.max())
        ssh = pyinterp.trivariate(interpolator,
                                  lon.flatten(),
                                  lat.flatten(),
                                  time,
                                  bounds_error=True,
                                  interpolator='bilinear').reshape(lon.shape)
        return ssh
