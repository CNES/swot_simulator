# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Karin noise
-----------
"""
from typing import Dict
import numpy as np

from .. import random_signal
from .. import settings


class Karin:
    """Karin instrumental error computed from random realization

    Args:
        parameters (settings.Parameters): Simulation settings
    """
    def __init__(self, parameters: settings.Parameters) -> None:
        assert parameters.karin_noise is not None

        # Store the generation parameters of the random signal.
        self.hsdt, self.x_ac, self.swh = random_signal.read_file_karin(
            parameters.karin_noise)

        self.size_grid = (parameters.delta_ac * parameters.delta_al)**0.5

        # Hack for unsmoothed products at high resolution
        if self.size_grid < 1:
            self.size_grid *= 8 / (40**.5)

    def generate(self, seed: int, x_al: np.ndarray, x_ac: np.ndarray,
                 swh: np.ndarray) -> Dict[str, np.ndarray]:
        """Generate the karin noise

        Args:
            seed (int): Random seed used to initialize the pseudo-random
                number generator.
            x_al (numpy.ndarray): Along track distance
            x_ac (numpy.ndarray): Across track distance
            swh (numpy.ndarray): Significant wave height. Used to modulate
                instrumental noise as a function of sea state.

        Returns:
            dict: variable name and errors simulated.
        """
        num_pixels = x_ac.shape[0]
        num_lines = x_al.shape[0]

        # Generate random noise for left and right part of the mast
        rng = np.random.default_rng(seed=seed)
        a_karin = rng.normal(0, 1, (num_lines, num_pixels))

        # Formula of karin noise as a function of x_ac (smile shape)
        sigma_karin = random_signal.interpolate_file_karin(
            swh, x_ac, self.hsdt, self.x_ac, self.swh) / self.size_grid

        # Compute random karin error
        return {"simulated_error_karin": sigma_karin * a_karin}
