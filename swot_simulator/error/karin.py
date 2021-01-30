# Copyright (c) 2020 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Karin noise
-----------
"""
from typing import Dict
import numpy as np

from . import utils
from .. import settings


class Karin:
    """Karin instrumental error computed from random realization

    Args:
        parameters (settings.Parameters): Simulation settings
    """
    def __init__(self, parameters: settings.Parameters) -> None:
        assert parameters.karin_noise is not None

        # Store the generation parameters of the random signal.
        self.hsdt, self.x_ac, self.swh = utils.read_file_karin(
                                                     parameters.karin_noise)

        self.delta_ac = parameters.delta_ac
        self.delta_al = parameters.delta_al
        self.nrand_karin = parameters.nrand_karin
        self.nseed = parameters.nseed

    def generate(self, x_al: np.array, x_ac: np.array,
                 swh: np.array) -> Dict[str, np.ndarray]:
        """Generate the karin noise

        Args:
            x_al (numpy.ndarray): Along track distance
            x_ac (numpy.ndarray): Across track distance
            swh (numpy.ndarray): Significant wave height. Used to modulate
                instrumental noise as a function of sea state.

        Returns:
            dict: variable name and errors simulated.
        """
        num_pixels = x_ac.shape[0]

        # Generate random noise for left and right part of the mast
        np.random.seed(self.nseed + 1)
        a_karin = np.random.normal(0, 1, (self.nrand_karin, num_pixels))

        # Formula of karin noise as a function of x_ac (smile shape)
        sigma_karin =  utils.interpolate_file_karin(swh, x_ac, self.hsdt,
                                                    self.x_ac, self.swh)
        size_grid = np.sqrt(self.delta_al * self.delta_ac)
        sigma_karin = sigma_karin / size_grid

        # Compute random karin error
        ai = ((x_al / self.delta_al) % self.nrand_karin).astype(np.uint64)
        return {"simulated_error_karin": sigma_karin * a_karin[ai, :]}
