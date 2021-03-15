# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Interpolation of the SSH HYCOM
==============================
"""
import numpy as np

from .base_impl import NetcdfLoader, CartesianGridHandler


class HYCOM(CartesianGridHandler):
    """
    Interpolation of the SSH HYCOM
    """

    def __init__(self, path: str):
        #: Decode the product date encoded from the file name.
        # hycom_GLBu0.08_191_2012031900_t012.nc
        loader = HYCOM.OverriddenNetcdfLoader(
            path,
            ssh_name="surf_el",
            pattern=r"hycom_GLBu0.08_191_(?P<date>\d{8})\d{2}_t\d{3}.nc",
            date_fmt="%Y%m%d",
        )
        super().__init__(loader)

    class OverriddenNetcdfLoader(NetcdfLoader):
        def load_dataset(self, first_date: np.datetime64, last_date: np.datetime64):
            ds = super().load_dataset(first_date, last_date)
            return ds.sel(depth=0)
