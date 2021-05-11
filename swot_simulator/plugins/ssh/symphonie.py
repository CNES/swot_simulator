# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Interpolation of the SSH SYMPHONIE
==================================
"""
from .. import data_handler


class SYMPHONIE(data_handler.IrregularGridHandler):
    """
    Interpolation of the SSH SYMPHONIE
    """
    def __init__(self, path: str):
        loader = data_handler.NetcdfLoader(
            path,
            date_fmt="%Y%m%d_%Hh",
            lon_name="lon_t",
            lat_name="lat_t",
            ssh_name="ssh_ib_detided",
            pattern=r"(?P<date>\d{8}_\d{2}h).detided.symphonie_bobshelf.nc")
        super().__init__(loader)
