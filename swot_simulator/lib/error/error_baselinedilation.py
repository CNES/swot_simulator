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


    def init_error_savesignal(self, delta_al:float, lambda_max:float,
                              npseudoper: int, len_repeat: int):
        """Compute random coefficients using the power spectrum """
        gencoef = utils.gen_rcoeff_signal1d(self.freq, self.PSbd,
                                            2 * delta_al, lambda_max,
                                            npseudoper, len_repeat)
        self.A_bd, self.phi_bd = gencoef

    def init_error_gensignal(self, ncomp1d: int):
        """Compute random signal using the power spectrum """
        gencoef = utils.gen_coeff_signal1d(self.freq, self.PSbd, ncomp1d)
        self.A_bd, self.phi_bd, self.fr_bd = gencoef


    '''
    def load_error(self, handler, parameters):
        self.A_bd = np.array(fid.variables['A_bd'][:]).squeeze()
        if np.shape(self.A_bd)[0] != parameters.ncomp1d:
            logger.error('{1} dimensions are different from ncomp1d={2}\n'
                         'remove {1} or adjust ncomp1d number in parameter'
                         'file'.format(parameters.file_coeff,
                                       parameters.ncomp1d))
            sys.exit(1)
        self.phi_bd = np.array(fid.variables['phi_bd'][:]).squeeze()
        self.fr_bd = np.array(fid.variables['fr_bd'][:]).squeeze()
    '''

    def make_error(self, x_al: np.array, al_cycle: float,
                   cycle: int, dal: float,
                   npseudoper: int, sat_const: dict, len_repeat:float,
                   lmax: Optional[float]=20000,
                   savesignal: Optional[bool]=True):
        Rearth = VOLUMETRIC_MEAN_RADIUS
        sat_elev = sat_const['height']
        Fka = sat_const['Fka']
        # - Compute baseline dilation using random coefficients or signals
        # previously initialized with the power spectrum
        if savesignal is True:
            xx = (np.float64(x_al[:]) + float(cycle * al_cycle)) % (len_repeat)
            dil = utils.gen_signal1d(xx, self.A_bd, self.phi_bd,
                                         2 * dal, lmax, npseudoper)

        else:
            dil = np.zeros((nal))
            for comp in range(0, self.ncomp1d):
                phase_x_al = (2. * np.pi * float(self.fr_bd[comp])
                              * (np.float64(x_al[:])
                              + float(cycle * al_cycle))) % (2.*np.pi)
                dil[:] = (dil[:] + 2 * self.A_bd[comp]
                      * np.cos(phase_x_al[:] + self.phi_bd[comp]))
        # - Compute the associated baseline dilation  error
        #   on the swath in m
        #TODO Check formulae
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
