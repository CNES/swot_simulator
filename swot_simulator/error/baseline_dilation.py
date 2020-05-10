import numpy as np
import xarray as xr

from . import utils
from .. import settings
from .. import VOLUMETRIC_MEAN_RADIUS, BASELINE


class BaselineDilation:
    def __init__(self, parameters: settings.Parameters,
                 dilation_psd: np.ndarray,
                 spatial_frequency: np.ndarray) -> None:
        # Store the generation parameters of the random signal.
        self.nseed = parameters.nseed + 1
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
        dil = utils.gen_signal_1d(self.freq,
                                  self.psbd,
                                  x_al,
                                  nseed=self.nseed,
                                  fmin=1 / self.len_repeat,
                                  fmax=1 / (2 * self.delta_al),
                                  alpha=10)

        # Compute the associated baseline dilation error on the swath in m
        return self.conversion_factor * dil

    def generate(self, x_al: np.ndarray, x_ac: np.ndarray) -> xr.DataArray:
        """TODO"""
        baseline_dilation_1d = self._generate_1d(x_al)
        return xr.DataArray(
            x_ac**2 * baseline_dilation_1d[:, np.newaxis],
            dims=("num_lines", "num_pixels"),
            name="baseline_dilation")
