import numpy as np
import xarray as xr
from . import utils

from .. import CELERITY
from .. import settings


class Timing:
    CONVERSION_FACTOR = CELERITY * 5e-13

    def __init__(self, parameters: settings.Parameters, timing_psd: np.ndarray,
                 spatial_frequency: np.ndarray) -> None:
        # Store the generation parameters of the random signal.
        self.nseed = parameters.nseed + 3
        self.len_repeat = parameters.len_repeat
        self.delta_al = parameters.delta_al

        # Get baseline dilation power spectrum
        self.timing_psd = timing_psd
        self.spatial_frequency = spatial_frequency

    def _generate_1d(self, x_al: np.ndarray) -> np.ndarray:
        # Generate 1d timing using the power spectrum:
        timing_l = utils.gen_signal_1d(self.spatial_frequency,
                                       self.timing_psd,
                                       x_al,
                                       nseed=self.nseed,
                                       fmin=1 / self.len_repeat,
                                       fmax=1 / (2 * self.delta_al),
                                       alpha=10)
        timing_r = utils.gen_signal_1d(self.spatial_frequency,
                                       self.timing_psd,
                                       x_al,
                                       nseed=self.nseed + 100,
                                       fmin=1 / self.len_repeat,
                                       fmax=1 / (2 * self.delta_al),
                                       alpha=10)
        # Compute the corresponding timing error on the swath in m
        return np.array([
            self.CONVERSION_FACTOR * timing_l,
            self.CONVERSION_FACTOR * timing_r
        ]).T

    def generate(self, x_al: np.ndarray, x_ac: np.ndarray) -> xr.DataArray:
        '''Reconstruct 2D errors from 1D instrumental error simulation'''
        timing_1d = self._generate_1d(x_al)
        num_pixels = x_ac.shape[0]
        swath_center = num_pixels // 2
        num_lines = timing_1d.shape[0]
        timing = np.empty((num_lines, num_pixels))
        ones_ac = np.ones((swath_center, ))
        timing[:, :swath_center] = ones_ac * timing_1d[:, 0, np.newaxis]
        timing[:, swath_center:] = ones_ac * timing_1d[:, 1, np.newaxis]

        return xr.DataArray(timing,
                            dims=("num_lines", "num_pixels"),
                            name="timing")
