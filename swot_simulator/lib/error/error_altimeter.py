import numpy as np
from . import utils
import sys
import logging
from typing import IO, Optional


# Define logging level for debug purposes
logger = logging.getLogger(__name__)


class error_stat():
    '''Class errornadir defines the error on the nadir.
    Random realisation of errors can be initialized using init_error.
    If the random realisations have already been computed and stored in file
    file_coeff, the random realisations are read directly using load_coeff.
    The correspondg errors on a swath can be computed using make_error. '''
    def __init__(self, ncomp1d: int, dal: float) -> None:
        self.ncomp1d = ncomp1d
        # - Define the sepctrum of the nadir instrument error
        f = np.arange(1./3000., 1./float(2.*dal), 1./3000.)
        PSD = 8 + 1.05 * 10**(-4) * f**(-2.2)
        indf = np.where(f < 0.00023627939582672978)
        PSD[indf] = 10**4
        # Convert spectrum in m2/cy
        self.PSD = PSD * 10**(-4)
        self.freq = f


    def init_error_savesignal(self, dal: int, lmax:float, npseudoper:int,
                              len_repeat: int,
                              savesignal: Optional[bool]=True) -> None:
        '''Initialization of errors: Random realisation of errors are computed
        using a known power spectrum.
        The outputs are the amplitude, the phase and the frequency of each
        random realisation.
        By default, there are ncomp2d=2000 random realisations for the
        wet tropo and ncomp1d=2000 random realisations for the nadir 1d
        spectrum error.'''
        if savesignal is True:
            gencoef = utils.gen_rcoeff_signal1d(self.freq, self.PSD, 2 * dal,
                                                lmax, npseudoper, len_repeat)
            self.A, self.phi = gencoef
        return None

    def init_error_gensignal(self, ncomp1d: int) -> None:
        gencoef = utils.gen_coeff_signal1d(self.freq, self.PSD, ncomp1d)
        self.A, self.phi, self.f = gencoef
        return None


    def make_error(self, x_al: np.array, al_cycle: float,
                   cycle: int, dal: float,
                   npseudoper:int, len_repeat:float,
                   lmax: Optional[float]=20000,
                   savesignal: Optional[bool]=True):
        ''' Build errors corresponding to each selected noise
        among the effect of the wet_tropo, and the instrumental error
        '''
        nal = np.shape(x_al)[0]
        # - Compute random noise of 10**2 cm**2/(km/cycle)
        # - Compute the correspond error on the nadir in m
        if savesignal is True:
            xx = (np.float64(x_al[:]) + float(cycle
                  * al_cycle)) % len_repeat
            errnadir = utils.gen_signal1d(xx, self.A, self.phi,
                                          2 * dal, lmax, npseudoper)
        else:
            errnadir = np.zeros((nal))
            for comp in range(0, self.ncomp1d):
                phase_x_al = (2. * np.pi * float(self.freq[comp])
                              * (np.float64(x_al[:])
                              + float(cycle * al_cycle))) % (2.*np.pi)
                errnadir[:] = (errnadir[:] + 2 * self.A[comp]
                               * np.cos(phase_x_al[:] + self.phi[comp]))
        # - Compute the correspond timing error on the swath in m
        self.nadir = errnadir[:]
        return None

