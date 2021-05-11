# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Interpolation of the SWH from WW3
=================================
"""
import numpy as np
import pyinterp.backends.xarray

from .. import data_handler


class WW3(data_handler.CartesianGridHandler):
    """
    Interpolation of the SWH from WW3.
    """
    def __init__(self, path: str):
        loader = data_handler.NetcdfLoader(path,
                                           date_fmt="%Y%m%d",
                                           ssh_name="hs",
                                           pattern=r"ww3.(?P<date>\w+)_hs\.nc")
        super().__init__(loader)

    def interpolate(self, lon: np.ndarray, lat: np.ndarray,
                    dates: np.ndarray) -> np.ndarray:
        """Interpolate the SSH to the required coordinates"""
        dataset = self.dataset_loader.load_dataset(
            dates.min(),  # type: ignore
            dates.max())  # type: ignore

        interpolator = pyinterp.backends.xarray.Grid3D(dataset.ssh)
        swh = interpolator.trivariate(
            dict(longitude=lon, latitude=lat, time=dates),
            interpolator="bilinear",
        )
        return swh
