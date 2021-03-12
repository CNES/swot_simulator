# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
from typing import Any, Dict, List, Optional, Tuple, Union
from . import math
from .plugins import ssh
from .plugins import swh as _swh
import numpy as np


def eval_config_file(filename: str) -> Dict:
    ...


class Parameters:
    area: Optional[Tuple[float, float, float, float]]
    beam_position: List[float]
    central_pixel: bool
    complete_product: bool
    corrected_roll_phase_dataset: Optional[str]
    cycle_duration: float
    delta_ac: float
    delta_al: float
    ephemeris_cols: Optional[Tuple[int, int, int]]
    ephemeris: Optional[str]
    error_spectrum: Optional[str]
    half_gap: float
    half_swath: float
    height: float
    karin_noise: Optional[str]
    len_repeat: float
    nadir: bool
    nbeam: int
    noise: List[str]
    nseed: int
    product_type: str
    requirement_bounds: Optional[Tuple[float, float]]
    shift_lon: Optional[float]
    shift_time: Optional[float]
    sigma: float
    ssh_plugin: Optional[ssh.Interface]
    swh_plugin: Optional[_swh.Interface]
    swath: bool
    swh: float
    hierarchical_groups: bool
    working_directory: str

    def __init__(self, override: Dict[str, Any]) -> None:
        ...

    @staticmethod
    def load_default() -> 'Parameters':
        ...

    @property
    def box(self) -> math.Box:
        ...

    def rng(self) -> np.random.Generator:
        ...


def template(python: bool = False) -> Union[str, Dict[str, Any]]:
    ...