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
    def __init__(self, dal: float) -> None:
        # - Define the sepctrum of the nadir instrument error
        f = np.arange(1./3000., 1./float(2.*dal), 1./3000.)
        PSD = 8 + 1.05 * 10**(-4) * f**(-2.2)
        indf = np.where(f < 0.00023627939582672978)
        PSD[indf] = 10**4
        # Convert spectrum in m2/cy
        self.PSD = PSD * 10**(-4)
        self.freq = f



    def make_error(self, x_al: np.array, dal: float,
                  len_repeat:float, nseed: Optional[int]=0):
        ''' Build errors corresponding to each selected noise
        among the effect of the wet_tropo, and the instrumental error
        '''
        # - Compute random noise of 10**2 cm**2/(km/cycle)
        # - Compute the correspond error on the nadir in m
        errnadir = utils.gen_signal1d(self.freq, self.PSD, x_al, nseed=nseed,
                           fmin=1./len_repeat, fmax=1./(2*dal), alpha=10)

        # - Compute the correspond timing error on the swath in m
        self.nadir = errnadir[:]
        return None

