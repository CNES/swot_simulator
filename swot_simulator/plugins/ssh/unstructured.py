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

from .. import UnstructuredHandler


class unstructured(UnstructuredHandler):
    """
    Interpolation of the SSH HYCOM
    """

    #: Decode the product date encoded from the file name.
    #hycom_GLBu0.08_191_2012031900_t012.nc
    PATTERN = re.compile(
        r"schout_baltic_(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2}).nc").search

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

        selected = self.ts["path"][(self.ts["date"] >= first_date)
                                   & (self.ts["date"] < last_date)]
        ds = xr.open_mfdataset(selected, concat_dim="time", combine="nested")
                               #decode_times=False)
        lon = ds.SCHISM_hgrid_node_x.values
        lat = ds.SCHISM_hgrid_node_y.values
        ssh = ds.elev

        time = [np.datetime64(x) for x in ds.time.values]
        list_mesh = []
        for t in range(np.shape(ssh)[0]):
            mesh = pyinterp.RTree()
            mesh.packing(np.vstack((lon[t, :], lat[t, :])).T, ssh[t, :])
            list_mesh.append(mesh)
        return list_mesh, time

    def interpolate(self, lon: np.ndarray, lat: np.ndarray,
                    time: np.ndarray) -> np.ndarray:
        """Interpolate the SSH to the required coordinates"""

        interpolator, mtime = self.load_dataset(time.min(), time.max())
        middle_time = time[int(np.shape(time)[0]/2)]
        _ind = np.argmin(abs(np.array(mtime) - middle_time))
        print('##index', _ind)
        mesh = interpolator[_ind]
        print(mesh)
        ssh, _ = mesh.inverse_distance_weighting(np.vstack((lon.flatten(),
                                                 lat.flatten())).T,
                                                 within=False,  # Extrapolation is forbidden
                                                 k=5,  # We are looking for at most 11 neighbours
                                                 num_threads=0)
        ssh = ssh.reshape(lon.shape)
        return ssh
