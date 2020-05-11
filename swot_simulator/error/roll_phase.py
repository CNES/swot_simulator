from typing import Iterator, Tuple
import numpy as np

from . import utils
from .. import settings
from .. import VOLUMETRIC_MEAN_RADIUS, CELERITY, F_KA, BASELINE


class RollPhase:
    def __init__(self, parameters: settings.Parameters, roll_psd: np.ndarray,
                 phase_psd: np.ndarray, spatial_frequency: np.ndarray) -> None:
        # Store the generation parameters of the random signal.
        self.nseed = parameters.nseed + 2
        self.len_repeat = parameters.len_repeat
        self.delta_al = parameters.delta_al

        # Get baseline dilation power spectrum
        self.roll_psd = roll_psd
        self.phase_psd = phase_psd
        self.spatial_frequency = spatial_frequency

        # TODO
        height = parameters.height * 1e-3
        self.phase_conversion_factor = (
            1 / (F_KA * 2 * np.pi / CELERITY * BASELINE) *
            (1 + height / VOLUMETRIC_MEAN_RADIUS) * np.pi / 180 * 1e3)
        self.roll_conversion_factor = (
            (1 + height / VOLUMETRIC_MEAN_RADIUS) * np.pi / 180 / 3600) * 1e3

    def _generate_1d(self, x_al: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:

        # Compute roll angle using the power spectrum
        # Compute left and right phase angles the power spectrum
        theta = utils.gen_signal_1d(self.spatial_frequency,
                                    self.roll_psd,
                                    x_al,
                                    nseed=self.nseed,
                                    fmin=1 / self.len_repeat,
                                    fmax=1 / (2 * self.delta_al),
                                    alpha=10)
        theta_l = utils.gen_signal_1d(self.spatial_frequency,
                                      self.phase_psd,
                                      x_al,
                                      nseed=self.nseed + 100,
                                      fmin=1 / self.len_repeat,
                                      fmax=1 / (2 * self.delta_al),
                                      alpha=10)
        theta_r = utils.gen_signal_1d(self.spatial_frequency,
                                      self.phase_psd,
                                      x_al,
                                      nseed=self.nseed + 200,
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
    ) -> Iterator[Tuple[str, np.ndarray]]:
        """TODO"""
        roll_1d, phase_1d = self._generate_1d(x_al)
        num_pixels = x_ac.shape[0]
        swath_center = num_pixels // 2
        ac_l = x_ac[:swath_center]
        ac_r = x_ac[swath_center:]

        phase = np.full((phase_1d.shape[0], num_pixels), np.nan)
        phase[:, :swath_center] = ac_l * phase_1d[:, 0, np.newaxis]
        phase[:, swath_center:] = ac_r * phase_1d[:, 1, np.newaxis]
        yield ("phase", phase)

        # rollphase_est_1d = np.zeros((np.shape(phase_1d.T)))
        # rollphase_est = np.full((np.shape(rollphase_est_1d)[0], nac), np.nan)
        # rollphase_est[:, :mid_swath] = np.mat(rollphase_est_1d[:, 0]).T * ac_l
        # rollphase_est[:, mid_swath:] = np.mat(rollphase_est_1d[:, 1]).T * ac_r
        yield ("roll", x_ac * roll_1d[:, np.newaxis])
