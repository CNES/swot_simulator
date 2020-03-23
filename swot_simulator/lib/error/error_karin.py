import numpy as np
from . import utils
import sys
import logging
from typing import IO, Optional


# Define logging level for debug purposes
logger = logging.getLogger(__name__)

VOLUMETRIC_MEAN_RADIUS = 6371

class error_stat():
    '''Karin instrumental error computed from random realization '''
    def __init__(self, karin_file: str, nrand=2000) -> None:
        # Get baseline dilation power spectrum
        self.nrand = nrand
        self.karin_file = karin_file


    def init_error(self, nac: int, nseed:int):
        np.random.seed(nseed)
        # generate random noise for left and right part of the mast
        self.A_karin_l = np.random.normal(0.0, np.float64(1),
                                          (self.nrand, nac))
        np.random.seed(nseed+ 1)
        self.A_karin_r = np.random.normal(0.0, np.float64(1),
                                          (self.nrand, nac))

    '''
    def load_error(self, handler, parameters):
        _A_karin = handler.variables['A_karin_l'][:, :]
        self.A_karin_l = np.array(_A_karin).squeeze()
        _A_karin = handler.variables['A_karin_r'][:, :]
        self.A_karin_r = np.array(_A_karin).squeeze()
        if np.shape(self.A_karin_l)[0] != parameters.nrandkarin:
            logger.error('{1} dimensions are different from nrandkarin={2}'
                         '\n remove {1} or adjust nrandkarin number in '
                         'parameter file'.format(parameters.file_coeff,
                                                 parameters.nrandkarin))
            sys.exit(1)
    '''

    def make_error(self, x_al: np.array, x_ac: np.array, al_cycle:float,
                   cycle: int, dal: int, dac: int, swh) -> None:
        # - Formula of karin noise as a function of x_ac (smile shape)
        # - Load Karin noise from file:
        karin_x_ac, hsdt = utils.read_file_karin(self.karin_file, swh)
        nal = np.shape(x_al)[0]
        nac = np.shape(x_ac)[0]
        self.karin = np.zeros((nal, nac))
        if isinstance(swh, (list, np.ndarray)):
            for ial in range(nal):
                sigma_karin = np.interp(np.abs(x_ac), karin_x_ac, hsdt[ial, :])
                size_grid = np.sqrt(np.float64(dal * dac))
                sigma_karin = sigma_karin / size_grid
                # - Compute random karin error
                Ai = (((np.float64(x_al[ial]) + float(cycle * al_cycle))
                      / dal) % self.nrand).astype('int')
                for j in range(0, nac):
                    self.karin[ial, j] = (sigma_karin[j]) * self.A_karin_r[Ai, j]

        else:
            sigma_karin = np.interp(np.abs(x_ac), karin_x_ac, hsdt)
            size_grid = np.sqrt(np.float64(dal * dac))
            sigma_karin = sigma_karin / size_grid
            # - Compute random karin error
            Ai = (((np.float64(x_al) + float(cycle * al_cycle))
                  / dal) % self.nrand).astype('int')
            for j in range(nac):
                self.karin[:, j] = (sigma_karin[j]) * self.A_karin_r[Ai, j]


    def reconstruct_2D_error(self):
        ''' Reconstruct 2D errors from 1D instrumental error simulation '''
        return self.karin

