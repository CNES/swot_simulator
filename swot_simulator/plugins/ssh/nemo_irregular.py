# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Interpolation of the SSH SYMPHONIE
==================================
"""
from .. import data_handler


class NEMO_IRREGULAR(data_handler.IrregularGridHandler):
    """
    Interpolation of the SSH SYMPHONIE

    Args:
        path (str): path to the folder containing the SYMPHONIE files.
    """

    def __init__(self, path: str):
        print(path)
        loader = data_handler.NetcdfLoader(
            path,
            date_fmt="%Y%m%d",
            lon_name="nav_lon",
            lat_name="nav_lat",
            time_name="time_counter",
            ssh_name="sossheig",
            pattern=r"ORCA36-T426a_1hAV_(?P<date>\d{8})-(?P<date2>\d{8})_gridT.nc")
        super().__init__(loader)
