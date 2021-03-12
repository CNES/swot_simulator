# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Baseling dilation errors
------------------------
"""
from typing import Dict
import numpy as np

from .. import random_signal
from .. import settings
from .. import VOLUMETRIC_MEAN_RADIUS, BASELINE


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
        # Store the generation parameters of the random signal.
        self.rng = parameters.rng()
        self.len_repeat = parameters.len_repeat
        self.delta_al = parameters.delta_al

        # Get baseline dilation power spectrum
        self.psbd = dilation_psd
        self.freq = spatial_frequency

        # TODO
        height = parameters.height * 1e-3
        self.conversion_factor = -((1 + height / VOLUMETRIC_MEAN_RADIUS) /
                                   (height * BASELINE)) * 1e-3

    def _generate_1d(self, x_al: np.ndarray) -> np.ndarray:
        # Generate 1d baseline dilation using the power spectrum:
        dil = random_signal.gen_signal_1d(self.freq,
                                          self.psbd,
                                          x_al,
                                          rng=self.rng,
                                          fmin=1 / self.len_repeat,
                                          fmax=1 / (2 * self.delta_al),
                                          alpha=10)

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
