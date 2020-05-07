# Copyright (c) 2020 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
from typing import Any, Dict, List, Optional, Tuple
from . import math
from .plugins import ssh


def eval_config_file(filename: str) -> Dict:
    ...


class Parameters:
    area: Optional[Tuple[float, float, float, float]]
    beam_position: List[float]
    complete_product: bool
    cycle_duration: float
    delta_ac: float
    delta_al: float
    ephemeris_cols: Optional[Tuple[int, int, int]]
    ephemeris: Optional[str]
    error_spectrum_file: Optional[str]
    half_gap: float
    half_swath: float
    height: float
    karin_file: Optional[str]
    len_repeat: float
    nadir: bool
    nbeam: int
    noise: bool
    nrand_karin: int
    nseed: int
    roll_phase_file: str
    save_signal: bool
    shift_lon: Optional[float]
    shift_time: Optional[float]
    sigma: float
    ssh_plugin: Optional[ssh.Interface]
    swath: bool
    swh: int
    working_directory: str

    def __init__(self, override: Dict[str, Any]) -> None:
        ...

    @property
    def box(self) -> math.Box:
        ...
