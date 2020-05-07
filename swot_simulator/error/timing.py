from typing import IO, Optional
import numpy as np
import xarray as xr
from . import utils

from .. import VOLUMETRIC_MEAN_RADIUS


class ErrorStat:
    def __init__(self, ds: xr.Dataset) -> None:
        # Get baseline dilation power spectrum
        self.pstim = ds.timingPSD.data
        self.freq = ds.spatial_frequency.data
        self.timing1d: Optional[np.ndarray] = None

    def make_error(self,
                   x_al: np.ndarray,
                   dal: float,
                   sat_const: dict,
                   len_repeat: float,
                   nseed: int = 0):
        # Generate 1d timing using the power spectrum:
        tim_l = utils.gen_signal1d(self.freq,
                                   self.pstim,
                                   x_al,
                                   nseed=nseed,
                                   fmin=1 / len_repeat,
                                   fmax=1 / (2 * dal),
                                   alpha=10)
        tim_r = utils.gen_signal1d(self.freq,
                                   self.pstim,
                                   x_al,
                                   nseed=nseed + 100,
                                   fmin=1 / len_repeat,
                                   fmax=1 / (2 * dal),
                                   alpha=10)
        # - Compute the corresponding timing error on the swath in m
        _to_m = sat_const['C'] / 2 * 10**(-12)
        timing1d = np.concatenate(([_to_m * tim_l[:].T], [_to_m * tim_r[:].T]),
                                  axis=0)
        self.timing1d = timing1d.T

    def reconstruct_2d_error(self, x_ac: np.ndarray) -> np.array:
        ''' Reconstruct 2D errors from 1D instrumental error simulation '''
        nac = np.shape(x_ac)[0]
        mid_swath = nac // 2
        timing1d = self.timing1d.copy()
        nal = np.shape(timing1d)[0]
        timing = np.zeros((nal, nac))
        ones_ac = np.mat(np.ones((mid_swath)))
        timing[:, :mid_swath] = np.mat(timing1d[:, 0]).T * ones_ac
        timing[:, mid_swath:] = np.mat(timing1d[:, 1]).T * ones_ac
        return timing
