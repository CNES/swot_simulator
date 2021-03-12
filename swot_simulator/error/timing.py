# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Timing errors
-------------
"""
from typing import Dict
import numpy as np
from .. import random_signal

from .. import CELERITY
from .. import settings


class Timing:
    """
    Timing errors

    Args:
        parameters (settings.Parameters): Simulation settings
        timing_psd (numpy.ndarray): Power spectral density the timing error
        spatial_frequency (numpy.ndarray): Spatial frequency
    """
    CONVERSION_FACTOR = CELERITY * 5e-13

    def __init__(self, parameters: settings.Parameters, timing_psd: np.ndarray,
                 spatial_frequency: np.ndarray) -> None:
        # Store the generation parameters of the random signal.
        self.rng_l = parameters.rng()
        self.rng_r = parameters.rng()
        self.len_repeat = parameters.len_repeat
        self.delta_al = parameters.delta_al

        # Get baseline dilation power spectrum
        self.timing_psd = timing_psd
        self.spatial_frequency = spatial_frequency

    def _generate_1d(self, x_al: np.ndarray) -> np.ndarray:
        # Generate 1d timing using the power spectrum:
        timing_l = random_signal.gen_signal_1d(self.spatial_frequency,
                                               self.timing_psd,
                                               x_al,
                                               rng=self.rng_l,
                                               fmin=1 / self.len_repeat,
                                               fmax=1 / (2 * self.delta_al),
                                               alpha=10)
        timing_r = random_signal.gen_signal_1d(self.spatial_frequency,
                                               self.timing_psd,
                                               x_al,
                                               rng=self.rng_r,
                                               fmin=1 / self.len_repeat,
                                               fmax=1 / (2 * self.delta_al),
                                               alpha=10)
        # Compute the corresponding timing error on the swath in m
        return np.array([
            self.CONVERSION_FACTOR * timing_l,
            self.CONVERSION_FACTOR * timing_r
        ]).T

    def generate(self, x_al: np.ndarray,
                 x_ac: np.ndarray) -> Dict[str, np.ndarray]:
        """Generate timing errors

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
