# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Interpolation of the SSH from MITGCM interpolated model for WW3
===============================================================
"""
import numpy as np
import xarray as xr

from swot_simulator.plugins.ssh.base_impl import NetcdfLoader, CartesianGridHandler


class MITGCM_WW3(CartesianGridHandler):
    """
    Interpolation of the SSH from MITGCM interpolated for WWW3.
    """

    def __init__(self, path: str):
        loader = MITGCM_WW3.OverriddenNetcdfLoader(
            path,
            ssh_name="wlv",
            pattern=r"ww3.(?P<date>\w+)_wlv.nc",
            date_fmt="%Y%m%d",
        )
        super().__init__(loader)

    class OverriddenNetcdfLoader(NetcdfLoader):
        def load_dataset(self, first_date: np.datetime64, last_date: np.datetime64):
            selected = self.select_netcdf_files(first_date, last_date)
            ds = xr.open_mfdataset(
                self.ts["path"][selected],
                concat_dim=self.time_name,
                combine="nested",
                decode_times=True,
            )
            if self.time_name not in ds.coords:
                return ds.assign_coords({self.time_name: self.ts["dates"][selected]})
            return ds
