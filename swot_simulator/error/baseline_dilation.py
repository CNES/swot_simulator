# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Baseling dilation errors
------------------------
"""
from typing import Dict
import logging

import numpy as np

from .. import random_signal
from .. import settings
from .. import VOLUMETRIC_MEAN_RADIUS, BASELINE

#: Logger of this module
LOGGER = logging.getLogger(__name__)


class BaselineDilation:
    """
    Baseline dilation errors

    Args:
        parameters (settings.Parameters): Simulation settings
        dilation_psd (numpy.ndarray): Power spectral density of the baseline
            dilation.
        spatial_frequency (numpy.ndarray): Spatial frequency
    """

    def __init__(self, parameters: settings.Parameters,
                 dilation_psd: np.ndarray,
                 spatial_frequency: np.ndarray) -> None:
        LOGGER.info("Initialize baseline dilation error")
        delta_al = 2 * parameters.delta_al
        assert parameters.height is not None
        height = parameters.height * 1e-3
        self.conversion_factor = -((1 + height / VOLUMETRIC_MEAN_RADIUS) /
                                   (height * BASELINE)) * 1e-3

        self.signal = random_signal.Signal1D(spatial_frequency,
                                             dilation_psd,
                                             rng=parameters.rng(),
                                             fmin=1 / parameters.len_repeat,
                                             fmax=1 / delta_al,
                                             alpha=10)

    def _generate_1d(self, x_al: np.ndarray) -> np.ndarray:
        # Generate 1d baseline dilation using the power spectrum:
        dil = self.signal(x_al)

        # Compute the associated baseline dilation error on the swath in m
        return self.conversion_factor * dil

    def generate(self, x_al: np.ndarray,
                 x_ac: np.ndarray) -> Dict[str, np.ndarray]:
        """Generate the baseline dilation error

        Args:
            x_al (numpy.ndarray): Along track distance
            x_ac (numpy.ndarray): Across track distance

        Returns:
            dict: variable name and errors simulated.
        """
        baseline_dilation_1d = self._generate_1d(x_al)
        return {
            "simulated_error_baseline_dilation":
            x_ac**2 * baseline_dilation_1d[:, np.newaxis]
        }
