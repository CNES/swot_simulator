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

from .. import data_handler

LOGGER = logging.getLogger(__name__)


class MITGCM_WW3(data_handler.CartesianGridHandler):
    """
    Interpolation of the SSH from MITGCM interpolated for WWW3.
    """
    def __init__(self, path: str):
        loader = MITGCM_WW3.OverriddenNetcdfLoader(
            path,
            date_fmt="%Y%m%d",
            ssh_name="wlv",
            pattern=r"ww3.(?P<date>\d{8})_wlv.nc")
        super().__init__(loader)

    class OverriddenNetcdfLoader(data_handler.NetcdfLoader):
        def load_dataset(self, first_date: np.datetime64,
                         last_date: np.datetime64):
            LOGGER.debug("fetch data for %s, %s", first_date, last_date)
            selected = self.select_netcdf_files(first_date, last_date)
            dataset = xr.open_mfdataset(self.time_series["path"][selected],
                                        concat_dim=self.time_name,
                                        combine="nested",
                                        decode_times=True)

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
