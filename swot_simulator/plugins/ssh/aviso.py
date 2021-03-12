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

from swot_simulator.plugins.ssh.base_impl import CartesianGridHandler, NetcdfLoader


class AVISO(CartesianGridHandler):
    """
    Interpolation of the SSH AVISO (CMEMS L4 products).
    """

    def __init__(self, path: str):
        loader = NetcdfLoader(
            path,
            ssh_name="adt",
            pattern=r"phy_l4_(?P<date>\w+).nc",
            date_fmt="%Y%m%d_%H%M%S%f",
        )
        super().__init__(loader)

    def interpolate(
        self, lon: np.ndarray, lat: np.ndarray, times: np.ndarray
    ) -> np.ndarray:
        """Interpolate the SSH to the required coordinates"""
        ds = self.dataset_loader.load_dataset(times.min(), times.max())

        interpolator = pyinterp.backends.xarray.Grid3D(ds.ssh)
        ssh = interpolator.trivariate(
            dict(longitude=lon, latitude=lat, time=times.astype("datetime64[ns]")),
            interpolator="bilinear",
        )
        return ssh
