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
        self.PStim = ds.timingPSD.data
        self.freq = ds.spatial_frequency.data


    def init_error_savesignal(self, delta_al:float, lambda_max: float,
                               npseudoper: int, len_repeat:int)-> None:
        """Compute random coefficients using the power spectrum """
        gencoef = utils.gen_rcoeff_signal1d(self.freq, self.PStim,
                                            2 * delta_al, lambda_max,
                                            npseudoper, len_repeat)
        self.A_tim_l, self.phi_tim_l = gencoef
        gencoef = utils.gen_rcoeff_signal1d(self.freq, self.PStim,
                                            2 * delta_al, lambda_max,
                                            npseudoper, len_repeat)
        self.A_tim_r, self.phi_tim_r = gencoef

    def init_error_gensignal(ncomp1d: int)-> None:
        gencoef = utils.gen_coeff_signal1d(self.freq, self.PStim, ncomp1d)
        self.A_tim_l, self.phi_tim_l, self.fr_tim_l = gencoef
        gencoef = utils.gen_coeff_signal1d(self.freq, self.PStim, ncomp1d)
        self.A_tim_r, self.phi_tim_r, self.fr_tim_r = gencoef


    def make_error(self, x_al: np.array, al_cycle: float,
                   cycle: int, dal: float,
                   npseudoper: int, sat_const: dict, len_repeat:float,
                   lmax: Optional[float]=20000,
                   savesignal: Optional[bool]=True)-> None:
        Rearth = VOLUMETRIC_MEAN_RADIUS
        sat_elev = sat_const['height']
        # - Compute timing delay using random coefficients or signals
        # previously initialized with the power spectrum
        if savesignal is True:
            xx = (np.float64(x_al[:]) + float(cycle * al_cycle)) % (len_repeat)
            tim_l = utils.gen_signal1d(xx, self.A_tim_l, self.phi_tim_l,
                                       2 * dal, lmax, npseudoper)
            tim_r = utils.gen_signal1d(xx, self.A_tim_r, self.phi_tim_r,
                                           2 * dal, lmax, npseudoper)
        else:
            tim_l = np.zeros((nal))
            tim_r = np.zeros((nal))
            # - Compute the associated phase error on the swath in m
            for comp in range(0, ncomp1d):
                phase_x_al = (2. * np.pi * float(self.fr_tim_r[comp])
                              * (np.float64(x_al[:])
                              + float(cycle * al_cycle))) % (2.*np.pi)
                tim_r[:] = (tim_r[:] + 2 * self.A_tim_r[comp]
                            * np.cos(phase_x_al[:]
                            + self.phi_tim_r[comp]))
                phase_x_al = (2. * np.pi * float(self.fr_tim_l[comp])
                              * (np.float64(x_al[:])
                              + float(cycle * al_cycle))) % (2.*np.pi)
                tim_l[:] = (tim_l[:] + 2 * self.A_tim_l[comp]
                            * numpy.cos(phase_x_al[:]
                            + self.phi_tim_l[comp]))
        # - Compute the corresponding timing error on the swath in m
        _to_m = sat_const['C']/2 * 10**(-12)
        timing1d = np.concatenate(([_to_m * tim_l[:].T],
                                     [_to_m * tim_r[:].T]), axis=0)
        self.timing1d = timing1d.T


    def reconstruct_2D_error(self, x_ac: np.ndarray)-> np.array:
        ''' Reconstruct 2D errors from 1D instrumental error simulation '''
        nac = np.shape(x_ac)[0]
        mid_swath = int(nac / 2)
        ac_l = np.mat(x_ac[:mid_swath])
        ac_r = np.mat(x_ac[mid_swath:])
        timing1d = self.timing1d
        nal = np.shape(timing1d)[0]
        timing = np.zeros((nal, nac))
        ones_ac = np.mat(np.ones((mid_swath)))
        timing[:, :mid_swath] = np.mat(timing1d[:, 0]).T * ones_ac
        timing[:, mid_swath:] = np.mat(timing1d[:, 1]).T * ones_ac
        return timing
