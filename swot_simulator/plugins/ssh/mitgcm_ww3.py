# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Interpolation of the SSH from MITGCM interpolated model for WW3
===============================================================
"""
import logging

import numpy as np
import xarray as xr

from swot_simulator.plugins.ssh.base_impl import NetcdfLoader, CartesianGridHandler

LOGGER = logging.getLogger(__name__)


class MITGCM_WW3(CartesianGridHandler):
    """
    Interpolation of the SSH from MITGCM interpolated for WWW3.
    """

    def __init__(self, path: str):
        loader = MITGCM_WW3.OverriddenNetcdfLoader(
            path,
            ssh_name="wlv",
            pattern=r"ww3.(?P<date>\d{8})_wlv.nc",
            date_fmt="%Y%m%d",
        )
        super().__init__(loader)

    class OverriddenNetcdfLoader(NetcdfLoader):
        def load_dataset(self, first_date: np.datetime64, last_date: np.datetime64):
            LOGGER.debug(f"fetch data for {first_date}, {last_date}")
            selected = self.select_netcdf_files(first_date, last_date)
            ds = xr.open_mfdataset(
                self.ts["path"][selected],
                concat_dim=self.time_name,
                combine="nested",
                decode_times=True,
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
