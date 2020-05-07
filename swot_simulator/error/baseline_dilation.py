from typing import Optional
import numpy as np

from . import utils
from .. import VOLUMETRIC_MEAN_RADIUS


class ErrorStat:
    def __init__(self, ds) -> None:
        # Get baseline dilation power spectrum
        self.psbd = ds['dilationPSD'].data
        self.freq = ds['spatial_frequency'].data
        self.baseline_dilation1d: Optional[np.ndarray] = None

    def make_error(self,
                   x_al: np.ndarray,
                   dal: float,
                   sat_const: dict,
                   len_repeat: float,
                   nseed: int = 0):
        sat_elev = sat_const['height']
        # Fka = sat_const['Fka']
        # Generate 1d baseline dilation using the power spectrum:
        dil = utils.gen_signal1d(self.freq,
                                 self.psbd,
                                 x_al,
                                 nseed=nseed,
                                 fmin=1. / len_repeat,
                                 fmax=1. / (2 * dal),
                                 alpha=10)

        # - Compute the associated baseline dilation  error on the swath in m
        _to_km = -((1 + sat_elev / VOLUMETRIC_MEAN_RADIUS) /
                   (sat_elev * sat_const['B'])) * 10**-3
        self.baseline_dilation1d = _to_km * dil[:]

    def reconstruct_2d_error(self, x_ac: np.ndarray):
        """Reconstruct 2D errors from 1D instrumental error simulation"""
        ac2 = np.mat((x_ac)**2)
        baseline_dilation1d = np.mat(self.baseline_dilation1d)
        baseline_dilation = baseline_dilation1d.T * ac2
        return baseline_dilation
