# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Interpolation of the SSH AVISO
==============================
"""
import numpy as np
import pyinterp.backends.xarray

from .. import data_handler


class AVISO(data_handler.CartesianGridHandler):
    """
    Interpolation of the SSH AVISO (CMEMS L4 products).
    """
    def __init__(self, path: str):
        loader = data_handler.NetcdfLoader(
            path,
            date_fmt="%Y%m%d",
            lon_name='longitude',
            lat_name='latitude',
            ssh_name="adt",
            pattern=r"phy_l4_(?P<date>\d{8})_\d{8}.nc")
        super().__init__(loader)

    def interpolate(self, lon: np.ndarray, lat: np.ndarray,
                    dates: np.ndarray) -> np.ndarray:
        """Interpolate the SSH to the required coordinates"""
        dataset = self.dataset_loader.load_dataset(
            dates.min(),  # type: ignore
            dates.max())  # type: ignore

        interpolator = pyinterp.backends.xarray.Grid3D(dataset.ssh)
        ssh = interpolator.trivariate(
            dict(lon=lon, lat=lat, time=dates),
            interpolator="bilinear",
        )
        return ssh
