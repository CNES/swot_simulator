from typing import IO, Optional
import numpy as np
from . import utils


class ErrorStat():
    '''Karin instrumental error computed from random realization '''
    def __init__(self, karin_file: str, nrand=2000) -> None:
        # Get baseline dilation power spectrum
        self.nrand = nrand
        self.karin_file = karin_file
        self.a_karin_l: Optional[np.ndarray] = None
        self.a_karin_r: Optional[np.ndarray] = None
        self.karin: Optional[np.ndarray] = None

    def init_error(self, nac: int, nseed: int):
        np.random.seed(nseed)
        # generate random noise for left and right part of the mast
        self.a_karin_l = np.random.normal(0, 1, (self.nrand, nac))
        np.random.seed(nseed + 1)
        self.a_karin_r = np.random.normal(0, 1, (self.nrand, nac))

    def make_error(self, x_al: np.array, x_ac: np.array, al_cycle: float,
                   cycle: int, delta_al: float, delta_ac: float, swh) -> None:
        # - Formula of karin noise as a function of x_ac (smile shape)
        # - Load Karin noise from file:
        karin_x_ac, hsdt = utils.read_file_karin(self.karin_file, swh)
        nal = np.shape(x_al)[0]
        nac = np.shape(x_ac)[0]
        self.karin = np.zeros((nal, nac))
        if isinstance(swh, (list, np.ndarray)):
            for ial in range(nal):
                sigma_karin = np.interp(np.abs(x_ac), karin_x_ac, hsdt[ial, :])
                size_grid = np.sqrt(np.float64(delta_al * delta_ac))
                sigma_karin = sigma_karin / size_grid
                # - Compute random karin error
                ai = ((
                    (np.float64(x_al[ial]) + float(cycle * al_cycle)) / delta_al) %
                      self.nrand).astype('int')
                for jx in range(0, nac):
                    self.karin[ial, jx] = (
                        sigma_karin[jx]) * self.a_karin_r[ai, jx]

        else:
            sigma_karin = np.interp(np.abs(x_ac), karin_x_ac, hsdt)
            size_grid = np.sqrt(np.float64(delta_al * delta_ac))
            sigma_karin = sigma_karin / size_grid
            # - Compute random karin error
            ai = (((np.float64(x_al) + float(cycle * al_cycle)) / delta_al) %
                  self.nrand).astype('int')
            for jx in range(nac):
                self.karin[:, jx] = (sigma_karin[jx]) * self.a_karin_r[ai, jx]

    def reconstruct_2d_error(self):
        ''' Reconstruct 2D errors from 1D instrumental error simulation '''
        return self.karin
