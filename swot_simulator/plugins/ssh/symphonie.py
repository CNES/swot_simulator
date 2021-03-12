# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Interpolation of the SSH Symphonie
==================================
"""
from swot_simulator.plugins.ssh.base_impl import NetcdfLoader, IrregularGridHandler


class Symphonie(IrregularGridHandler):
    """
    Interpolation of the SSH NEMO (CMEMS L4 products).
    """

    def __init__(self, path: str):
        loader = NetcdfLoader(
            path,
            lon_name="lon_t",
            lat_name="lat_t",
            ssh_name="ssh_ib_detided",
            pattern=r"(?P<date>\w+).detided.symphonie_bobshelf.nc",
            date_fmt="%Y%m%d_%Hh",
        )
        super().__init__(loader)
