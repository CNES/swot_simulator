import numpy as np
import xarray as xr

from . import utils
from .. import settings


class Karin:
    """Karin instrumental error computed from random realization"""
    def __init__(self, parameters: settings.Parameters,
                 num_pixels: int) -> None:
        assert parameters.karin_file is not None

        # Store the generation parameters of the random signal.
        self.x_ac, self.hsdt = utils.read_file_karin(parameters.karin_file,
                                                     parameters.swh)

        self.delta_al = parameters.delta_al
        self.delta_ac = parameters.delta_ac
        self.nrand_karin = parameters.nrand_karin

        # Generate random noise for left and right part of the mast
        np.random.seed(parameters.nseed)
        self.a_karin_l = np.random.normal(0, 1, (self.nrand_karin, num_pixels))

        np.random.seed(parameters.nseed + 1)
        self.a_karin_r = np.random.normal(0, 1, (self.nrand_karin, num_pixels))

    def generate(self, x_al: np.array, x_ac: np.array,
                 curvilinear_distance: float, cycle: int) -> xr.DataArray:
        # Formula of karin noise as a function of x_ac (smile shape)
        sigma_karin = np.interp(np.abs(x_ac), self.x_ac, self.hsdt)
        size_grid = np.sqrt(self.delta_al * self.delta_ac)
        sigma_karin = sigma_karin / size_grid

        # Compute random karin error
        ai = (((x_al + cycle * curvilinear_distance) / self.delta_al) %
              self.nrand_karin).astype(np.uint64)
        return xr.DataArray(sigma_karin * self.a_karin_r[ai, :],
                            dims=("num_lines", "num_pixels"),
                            name="karin")
