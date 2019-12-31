# Copyright (c) 2020 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
from typing import Any, Dict, Optional, Tuple
from . import math
from .plugins import ssh


def eval_config_file(filename: str) -> Dict:
    ...


class Parameters:
    area: Optional[Tuple[float, float, float, float]]
    cycle_duration: float
    delta_al: float
    delta_ac: float
    ephemeris: Optional[str]
    ephemeris_cols: Optional[Tuple[int, int, int]]
    complete_product: bool
    half_gap: float
    half_swath: float
    height: float
    nadir: bool
    swath: bool
    ssh_plugin: Optional[ssh.Interface]
    shift_lon: Optional[float]
    shift_time: Optional[float]
    working_directory: str

    def __init__(self, override: Dict[str, Any]) -> None:
        ...

    @property
    def box(self) -> math.Box:
        ...
