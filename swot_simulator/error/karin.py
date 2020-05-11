from typing import Iterator, Tuple
import numpy as np

from . import utils
from .. import settings


class Karin:
    """Karin instrumental error computed from random realization"""
    def __init__(self, parameters: settings.Parameters) -> None:
        assert parameters.karin_noise is not None

        # Store the generation parameters of the random signal.
        self.x_ac, self.hsdt = utils.read_file_karin(parameters.karin_noise,
                                                     parameters.swh)

        self.delta_ac = parameters.delta_ac
        self.delta_al = parameters.delta_al
        self.nrand_karin = parameters.nrand_karin
        self.nseed = parameters.nseed

    def generate(self, x_al: np.array, x_ac: np.array,
                 curvilinear_distance: float,
                 cycle: int) -> Iterator[Tuple[str, np.ndarray]]:
        num_pixels = x_ac.shape[0]

        # Generate random noise for left and right part of the mast
        np.random.seed(self.nseed + 1)
        a_karin = np.random.normal(0, 1, (self.nrand_karin, num_pixels))

        # Formula of karin noise as a function of x_ac (smile shape)
        sigma_karin = np.interp(np.abs(x_ac), self.x_ac, self.hsdt)
        size_grid = np.sqrt(self.delta_al * self.delta_ac)
        sigma_karin = sigma_karin / size_grid

        # Compute random karin error
        ai = (((x_al + cycle * curvilinear_distance) / self.delta_al) %
              self.nrand_karin).astype(np.uint64)
        yield ("karin", sigma_karin * a_karin[ai, :])
