import numpy as np
import xarray as xr
from . import utils

from . import ErrorStat as Base
from .. import CELERITY
from .. import settings


class ErrorStat(Base):
    TO_M = CELERITY / 2 * 1e-12

    def __init__(self, ds: xr.Dataset,
                 parameters: settings.Parameters) -> None:
        super().__init__(parameters)
        self.nseed += 3
        # Get baseline dilation power spectrum
        self.pstim = ds.timingPSD.data
        self.freq = ds.spatial_frequency.data

    def make_error(self, x_al: np.ndarray) -> np.ndarray:
        # Generate 1d timing using the power spectrum:
        tim_l = utils.gen_signal_1d(self.freq,
                                   self.pstim,
                                   x_al,
                                   nseed=self.nseed,
                                   fmin=1 / self.len_repeat,
                                   fmax=1 / (2 * self.delta_al),
                                   alpha=10)
        tim_r = utils.gen_signal_1d(self.freq,
                                   self.pstim,
                                   x_al,
                                   nseed=self.nseed + 100,
                                   fmin=1 / self.len_repeat,
                                   fmax=1 / (2 * self.delta_al),
                                   alpha=10)
        # - Compute the corresponding timing error on the swath in m
        timing_1d = np.concatenate(
            ([self.TO_M * tim_l[:].T], [self.TO_M * tim_r[:].T]), axis=0)
        return timing_1d.T

    def reconstruct_2d_error(self, x_ac: np.ndarray,
                             timing_1d: np.ndarray) -> np.array:
        ''' Reconstruct 2D errors from 1D instrumental error simulation '''
        nac = np.shape(x_ac)[0]
        mid_swath = nac // 2
        nal = np.shape(timing_1d)[0]
        timing = np.zeros((nal, nac))
        ones_ac = np.mat(np.ones((mid_swath)))
        timing[:, :mid_swath] = np.mat(timing_1d[:, 0]).T * ones_ac
        timing[:, mid_swath:] = np.mat(timing_1d[:, 1]).T * ones_ac
        return timing
