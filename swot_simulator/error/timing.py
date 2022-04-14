# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Timing errors
-------------
"""
from typing import Dict
import logging

import numpy as np

from .. import CELERITY, random_signal, settings

#: Logger of this module
LOGGER = logging.getLogger(__name__)


class Timing:
    """Timing errors.

    Args:
        parameters (settings.Parameters): Simulation settings
        timing_psd (numpy.ndarray): Power spectral density the timing error
        spatial_frequency (numpy.ndarray): Spatial frequency
    """
    CONVERSION_FACTOR = CELERITY * 5e-13

    def __init__(self, parameters: settings.Parameters, timing_psd: np.ndarray,
                 spatial_frequency: np.ndarray) -> None:
        LOGGER.info("Initialize timing error")
        # Store the generation parameters of the random signal.
        delta_al = 2 * parameters.delta_al
        self.timing_l = random_signal.Signal1D(spatial_frequency,
                                               timing_psd,
                                               rng=parameters.rng(),
                                               fmin=1 / parameters.len_repeat,
                                               fmax=1 / delta_al,
                                               alpha=10)
        self.timing_r = random_signal.Signal1D(spatial_frequency,
                                               timing_psd,
                                               rng=parameters.rng(),
                                               fmin=1 / parameters.len_repeat,
                                               fmax=1 / delta_al,
                                               alpha=10)

    def _generate_1d(self, x_al: np.ndarray) -> np.ndarray:
        # Generate 1d timing using the power spectrum:
        timing_l = self.timing_l(x_al)
        timing_r = self.timing_r(x_al)
        # Compute the corresponding timing error on the swath in m
        return np.array([
            self.CONVERSION_FACTOR * timing_l,
            self.CONVERSION_FACTOR * timing_r
        ]).T

    def generate(self, x_al: np.ndarray,
                 x_ac: np.ndarray) -> Dict[str, np.ndarray]:
        """Generate timing errors.

        Args:
            x_al (numpy.ndarray): Along track distance
            x_ac (numpy.ndarray): Across track distance

        Returns:
            dict: variable name and errors simulated.
        """
        timing_1d = self._generate_1d(x_al)
        num_pixels = x_ac.shape[0]
        swath_center = num_pixels // 2
        num_lines = timing_1d.shape[0]
        timing = np.empty((num_lines, num_pixels))
        ones_ac = np.ones((swath_center, ))
        timing[:, :swath_center] = ones_ac * timing_1d[:, 0, np.newaxis]
        timing[:, swath_center:] = ones_ac * timing_1d[:, 1, np.newaxis]

        return {"simulated_error_timing": timing}
