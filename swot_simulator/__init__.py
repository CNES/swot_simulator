# Copyright (c) 2020 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Constants that are defined for SWOT simulator
---------------------------------------------
"""
# EARTH CONSTANTS

#: Volumetric mean radius (km)
VOLUMETRIC_MEAN_RADIUS: float = 6371.0

#: Convert degree to km
DEG_2_KM: float = 111.11

#: CELERITY
CELERITY: float = 299800000.0

#: Baseline (m)
BASELINE = 10

#: KA frequency (in Hz)
F_KA = 35750000000.0

#: Product types
PRODUCT_TYPE = {
    "expert": "l2b-expert.xml",
    "basic": "l2b-ssh.xml",
    "wind_wave": "l2b-windwave.xml"
}
