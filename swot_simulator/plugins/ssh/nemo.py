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

from .. import data_handler


class NEMO(data_handler.CartesianGridHandler):
    """
    Interpolation of the SSH NEMO (CMEMS L4 products).
    """
    def __init__(self, path: str):
        loader = data_handler.NetcdfLoader(
            path,
            date_fmt="%Y%m%d_PGS_%H",
            ssh_name="sossheigh_corrected",
            time_name="time_counter",
            pattern=r"ORCA12-T.*_1h_grid_T_(?P<date>\w+).nc")
        super().__init__(loader)

    def interpolate(self, lon: np.ndarray, lat: np.ndarray,
                    dates: np.ndarray) -> np.ndarray:
        """Interpolate the SSH to the required coordinates"""
        dataset = self.dataset_loader.load_dataset(
            dates.min(),  # type: ignore
            dates.max())  # type: ignore

        interpolator = pyinterp.backends.xarray.Grid3D(dataset.ssh)
        ssh = interpolator.trivariate(
            dict(longitude=lon, latitude=lat, time_counter=dates),
            interpolator="bilinear",
        )
        return ssh
