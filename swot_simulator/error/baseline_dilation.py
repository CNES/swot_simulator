import numpy as np
import xarray as xr

from . import ErrorStat as Base
from . import utils
from .. import settings
from .. import VOLUMETRIC_MEAN_RADIUS, BASELINE


class ErrorStat(Base):
    def __init__(self, ds: xr.Dataset,
                 parameters: settings.Parameters) -> None:
        super().__init__(parameters)
        self.nseed += 1
        # Get baseline dilation power spectrum
        self.psbd = ds['dilationPSD'].data
        self.freq = ds['spatial_frequency'].data

    def make_error(self, x_al: np.ndarray) -> np.ndarray:
        # Generate 1d baseline dilation using the power spectrum:
        dil = utils.gen_signal_1d(self.freq,
                                 self.psbd,
                                 x_al,
                                 nseed=self.nseed,
                                 fmin=1 / self.len_repeat,
                                 fmax=1 / (2 * self.delta_al),
                                 alpha=10)

        # - Compute the associated baseline dilation  error on the swath in m
        _to_km = -((1 + self.height / VOLUMETRIC_MEAN_RADIUS) /
                   (self.height * BASELINE)) * 10**-3
        return _to_km * dil[:]

    def reconstruct_2d_error(self, x_ac: np.ndarray,
                             baseline_dilation_1d: np.ndarray) -> np.ndarray:
        """Reconstruct 2D errors from 1D instrumental error simulation"""
        ac2 = np.mat((x_ac)**2)
        baseline_dilation1d = np.mat(baseline_dilation_1d)
        baseline_dilation = baseline_dilation1d.T * ac2
        return baseline_dilation
