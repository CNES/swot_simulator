# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Roll errors
-----------
"""
from typing import Dict, Optional, Tuple
import numpy as np

from .. import random_signal
from .. import settings
from .. import VOLUMETRIC_MEAN_RADIUS, CELERITY, F_KA, BASELINE


class RollPhase:
    """
    Roll errors

    Args:
        parameters (settings.Parameters): Simulation settings
        roll_psd (numpy.ndarray): Power spectral density the roll control angle
        gyro_psd (numpy.ndarray): Power spectral density the roll error
            knowledge (also called gyro error)
        phase_psd (numpy.ndarray): Power spectral density the error angle
        spatial_frequency (numpy.ndarray): Spatial frequency
    """
    def __init__(self, parameters: settings.Parameters, roll_psd: np.ndarray,
                 gyro_psd: np.ndarray, phase_psd: np.ndarray,
                 spatial_frequency: np.ndarray) -> None:
        # Store the generation parameters of the random signal.
        self.rng_theta = parameters.rng()
        self.rng_theta_l = parameters.rng()
        self.rng_theta_r = parameters.rng()
        self.len_repeat = parameters.len_repeat
        self.delta_al = parameters.delta_al

        # Get baseline dilation power spectrum
        self.roll_psd = roll_psd + gyro_psd
        indexes = np.where(spatial_frequency <= 1 / 40000)[0]
        self.roll_psd[indexes] = self.roll_psd[indexes][-1]
        self.phase_psd = phase_psd
        self.spatial_frequency = spatial_frequency

        # TODO
        assert parameters.height is not None
        height = parameters.height * 1e-3
        self.phase_conversion_factor = (
            1 / (F_KA * 2 * np.pi / CELERITY * BASELINE) *
            (1 + height / VOLUMETRIC_MEAN_RADIUS) * np.pi / 180 * 1e3)
        self.roll_conversion_factor = (
            (1 + height / VOLUMETRIC_MEAN_RADIUS) * np.pi / 180 / 3600) * 1e3

    def _generate_1d(self, x_al: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:

        # Compute roll angle using the power spectrum
        # Compute left and right phase angles the power spectrum
        theta = random_signal.gen_signal_1d(self.spatial_frequency,
                                            self.roll_psd,
                                            x_al,
                                            rng=self.rng_theta,
                                            fmin=1 / self.len_repeat,
                                            fmax=1 / (2 * self.delta_al),
                                            alpha=10)
        theta_l = random_signal.gen_signal_1d(self.spatial_frequency,
                                              self.phase_psd,
                                              x_al,
                                              rng=self.rng_theta_l,
                                              fmin=1 / self.len_repeat,
                                              fmax=1 / (2 * self.delta_al),
                                              alpha=10)
        theta_r = random_signal.gen_signal_1d(self.spatial_frequency,
                                              self.phase_psd,
                                              x_al,
                                              rng=self.rng_theta_r,
                                              fmin=1 / self.len_repeat,
                                              fmax=1 / (2 * self.delta_al),
                                              alpha=10)
        # Compute the associated roll  error on the swath in m
        roll = self.roll_conversion_factor * theta
        # Compute the associated phase error on the swath in m
        phase = np.array([
            self.phase_conversion_factor * theta_l.T,
            self.phase_conversion_factor * theta_r.T
        ])
        return roll, phase.T

    def generate(
        self,
        x_al: np.ndarray,
        x_ac: np.ndarray,
    ) -> Dict[str, np.ndarray]:
        """Generate roll and phase errors

        Args:
            x_al (numpy.ndarray): Along track distance
            x_ac (numpy.ndarray): Across track distance

        Returns:
            dict: variable name and errors simulated.
        """
        roll_1d, phase_1d = self._generate_1d(x_al)
        num_pixels = x_ac.shape[0]
        swath_center = num_pixels // 2
        ac_l = x_ac[:swath_center]
        ac_r = x_ac[swath_center:]

        phase = np.full((phase_1d.shape[0], num_pixels), np.nan)
        phase[:, :swath_center] = ac_l * phase_1d[:, 0, np.newaxis]
        phase[:, swath_center:] = ac_r * phase_1d[:, 1, np.newaxis]

        # rollphase_est_1d = np.zeros((np.shape(phase_1d.T)))
        # rollphase_est = np.full((np.shape(rollphase_est_1d)[0], nac), np.nan)
        # rollphase_est[:, :mid_swath] = np.mat(rollphase_est_1d[:, 0]).T * ac_l
        # rollphase_est[:, mid_swath:] = np.mat(rollphase_est_1d[:, 1]).T * ac_r
        return {
            "simulated_error_roll": x_ac * roll_1d[:, np.newaxis],
            "simulated_error_phase": phase
        }
