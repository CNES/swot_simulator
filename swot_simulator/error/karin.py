import numpy as np

from . import ErrorStat as Base
from . import utils
from .. import settings


class ErrorStat(Base):
    """Karin instrumental error computed from random realization"""
    def __init__(self, parameters: settings.Parameters, nac: int) -> None:
        super().__init__(parameters)
        # Get baseline dilation power spectrum
        self.nrand = parameters.nrand_karin
        assert parameters.karin_file is not None
        self.karin_file = parameters.karin_file
        self.swh = parameters.swh

        # generate random noise for left and right part of the mast
        np.random.seed(parameters.nseed)
        self.a_karin_l = np.random.normal(0, 1, (self.nrand, nac))

        np.random.seed(parameters.nseed + 1)
        self.a_karin_r = np.random.normal(0, 1, (self.nrand, nac))

    def make_error(self, x_al: np.array, x_ac: np.array,
                   curvilinear_distance: float, cycle: int) -> np.ndarray:
        # - Formula of karin noise as a function of x_ac (smile shape)
        # - Load Karin noise from file:
        karin_x_ac, hsdt = utils.read_file_karin(self.karin_file, self.swh)
        nal = np.shape(x_al)[0]
        nac = np.shape(x_ac)[0]
        karin = np.zeros((nal, nac))

        sigma_karin = np.interp(np.abs(x_ac), karin_x_ac, hsdt)
        size_grid = np.sqrt(self.delta_al * self.delta_ac)
        sigma_karin = sigma_karin / size_grid
        # - Compute random karin error
        ai = (((np.float64(x_al) + float(cycle * curvilinear_distance)) /
               self.delta_al) % self.nrand).astype('int')
        for jx in range(nac):
            karin[:, jx] = (sigma_karin[jx]) * self.a_karin_r[ai, jx]

        return karin
