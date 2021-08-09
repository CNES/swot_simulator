# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Interpolation of the SSH NEMO
==============================
"""
import numpy as np
import pyinterp.backends.xarray
import pyinterp
import time
import logging

from .. import data_handler

LOGGER = logging.getLogger(__name__)


class WW3(data_handler.IrregularGridHandler):
    """
    Interpolation of the SSH NEMO (CMEMS L4 products).
    """
    def __init__(self, path: str):
        loader = data_handler.NetcdfLoader(
            path,
            date_fmt="%Y%m%d",
            ssh_name="sossheig",
            lon_name="nav_lon",
            lat_name="nav_lat",

            time_name="time_counter",
            pattern=r"PSY4V3R1_1hAV_(?P<date>\w+)_\d{8}_gridT.nc")
        super().__init__(loader)


    def interpolate(self, lon: np.ndarray, lat: np.ndarray, dates: np.ndarray):
        """Interpolate the SSH for the given coordinates"""
        # Spatial interpolation of the SSH on the different selected grids.
        dataset = self.dataset_loader.load_dataset(
            dates.min(),  # type: ignore
            dates.max())  # type: ignore
        dataset = dataset.sel(deptht=0)
        dates_p = dataset.time.load().data
        lon_p = dataset.lon.data
        lat_p = dataset.lat.data
        ssh_p = dataset.ssh.data

        start_time = time.time()
        layers = []
        for index in range(len(ssh_p)):
            layers.append(
                self._spatial_interp(ssh_p[index, :, :], lon_p, lat_p, lon, lat))
        layers = np.stack(layers)
        LOGGER.debug("interpolation completed in %.2fs for period %s, %s",
                     time.time() - start_time, dates.min(), dates.max())
        print('bouhtime', dates_p)
        dates_p = dates_p.astype("int64") * 10**6
        ssh_interp = data_handler._time_interp(dates_p.astype("int64"),layers,
                                   dates.astype("datetime64[us]").astype("int64"),)
        print('bouhtime', ssh_interp)
        LOGGER.debug('bouhssh')
        # Time interpolation of the SSH.
        return data_handler._time_interp(
            dates_p.astype("int64"),
            layers,
            dates.astype("datetime64[us]").astype("int64"),)

