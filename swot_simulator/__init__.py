# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Constants that are defined for SWOT simulator
---------------------------------------------
"""
from . import version

#: Module Version
__version__ = version.release()

#: Module release date
__date__ = version.date()

del version

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

#: Basic product
BASIC = "basic"

#: Expert product
EXPERT = "expert"

#: High resolution
UNSMOOTHED = "unsmoothed"

#: Wind-wave product
WIND_WAVE = "wind_wave"

#: Product types
PRODUCT_TYPE = {
    BASIC: "l2b-ssh.xml",
    EXPERT: "l2b-expert.xml",
    UNSMOOTHED: "l2a-highres.xml",
    WIND_WAVE: "l2b-windwave.xml"
}
