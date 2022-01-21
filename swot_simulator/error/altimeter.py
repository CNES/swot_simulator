# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Altimeter instrument error
--------------------------
"""
from typing import Dict
import logging

import numpy as np

from .. import random_signal
from .. import settings

#: Logger of this module
LOGGER = logging.getLogger(__name__)


class Altimeter:
    """Altimeter instrument error.

    Args:
        parameters (settings.Parameters): Simulation settings
    """

    def __init__(self, parameters: settings.Parameters) -> None:
        LOGGER.info("Initialize altimeter instrument error")
        # Store the generation parameters of the random signal.
        delta_al = 2 * parameters.delta_al
        # Define the sepctrum of the nadir instrument error
        freq = np.arange(1 / 3000, 1 / delta_al, 1 / 3000)
        psd = 8 + 1.05 * 1e-4 * freq**(-2.2)
        psd[freq < 0.00023627939582672978] = 1e4

        self.signal = random_signal.Signal1D(
            freq,
            # Convert spectrum in m2/cy
            psd * 1e-4,
            parameters.rng(),
            fmin=1 / parameters.len_repeat,
            fmax=1 / delta_al,
            alpha=10)

    def generate(self, x_al: np.ndarray) -> Dict[str, np.ndarray]:
        """Generate altimeter instrument error.

        Args:
            x_al (numpy.ndarray): Along track distance

        Returns:
            dict: variable name and errors simulated.
        """
        # Compute random noise of 10**2 cm**2/(km/cycle)
        # Compute the correspond error on the nadir in m
        error = self.signal(x_al)
        return {"simulated_error_altimeter": error}
