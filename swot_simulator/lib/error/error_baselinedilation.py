import numpy as np
from . import utils
import sys
import logging
from typing import IO, Optional


# Define logging level for debug purposes
logger = logging.getLogger(__name__)

VOLUMETRIC_MEAN_RADIUS = 6371

class error_stat():
    def __init__(self, ds) -> None:
        # Get baseline dilation power spectrum
        self.PSbd = ds.dilationPSD.data
        self.freq = ds.spatial_frequency.data


    def make_error(self, x_al: np.ndarray, dal: float, sat_const: dict,
                   len_repeat:float, nseed: Optional[int]=0):
        Rearth = VOLUMETRIC_MEAN_RADIUS
        sat_elev = sat_const['height']
        Fka = sat_const['Fka']
        # Generate 1d baseline dilation using the power spectrum:
        dil = utils.gen_signal1d(self.freq, self.PSbd, x_al, nseed=nseed,
                           fmin=1./len_repeat, fmax=1./(2*dal), alpha=10)
        # - Compute the associated baseline dilation  error on the swath in m
        _to_km = -((1 + sat_elev / Rearth) #* 10**-6 #* 10**6
                  /(sat_elev * sat_const['B'])) * 10**-3
        self.baseline_dilation1d = _to_km * dil[:]
        del dil

    def reconstruct_2D_error(self, x_ac: np.ndarray):
        """ Reconstruct 2D errors from 1D instrumental error simulation """
        nac = np.shape(x_ac)[0]
        mid_swath = int(nac / 2)
        ac_l = np.mat(x_ac[:mid_swath])
        ac_r = np.mat(x_ac[mid_swath:])
        ac2 = np.mat((x_ac)**2)
        baseline_dilation1d =  np.mat(self.baseline_dilation1d)
        baseline_dilation = baseline_dilation1d.T * ac2
        return baseline_dilation
