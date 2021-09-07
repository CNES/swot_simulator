# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Orbital error
-------------
"""
from typing import Dict, Tuple

#
import dask.array as da
import numpy as np

#
from .. import random_signal
from .. import settings
from .. import VOLUMETRIC_MEAN_RADIUS

#: Signal amplitude of the orbital error in micro-radians
AMPLITUDE = 100

#: Delta T of the spatial sampling in seconds
DT = 60


def _orbital_error_spectrum(
        orbit_duration: np.timedelta64,
        rng: np.random.Generator) -> Tuple[np.ndarray, float]:
    """Calculate orbital error spectrum

    Args:
        orbit_duration (float): Orbit duration in fractional days
        rng (np.random.Generator): Random number generator

    Returns:
        tuple: (yg, fmaxr)
    """
    df = 1 / (1000 * 86400)
    spatial_frequency = np.arange(df, 1 / DT, df)
    orbital_frequency = 1 / float(
        orbit_duration.astype("timedelta64[us]").astype("float64") * 1e-6)
    sigma_peak = orbital_frequency / 1000
    ps_orbital = np.exp(-0.5 * (spatial_frequency - orbital_frequency)**2 /
                        sigma_peak**2)
    ps_orbital[ps_orbital < 1 / 1000] = 0.
    ps_orbital /= np.sum(ps_orbital * df)
    ps_orbital *= AMPLITUDE**2
    return random_signal.gen_psd_1d(spatial_frequency,
                                    ps_orbital,
                                    rng,
                                    alpha=10)


class Orbital:
    """
    Simulate the error orbital

    Args:
        parameters (Parameters): Simulation parameters.
        orbit_duration (np.timedelta64): Orbit duration.
    """
    def __init__(self, parameters: settings.Parameters,
                 orbit_duration: np.timedelta64) -> None:
        yg, self.fmaxr = _orbital_error_spectrum(orbit_duration,
                                                 parameters.rng())
        self.yg = da.from_array(yg, name="orbital_error").persist()

        assert parameters.height is not None
        height = parameters.height * 1e-3
        self.conversion_factor = (1 + height / VOLUMETRIC_MEAN_RADIUS) * 1e-3

    def generate(
        self,
        time: np.ndarray,
        x_ac: np.ndarray,
    ) -> Dict[str, np.ndarray]:
        """Generate orbital error

        Args:
            time (np.ndarray): time vector

        Returns:
            np.ndarray: orbital error
        """
        time = time.astype("datetime64[us]").astype("float64") * 1e-6
        xg = np.linspace(0, 0.5 / self.fmaxr * self.yg.shape[0],
                         self.yg.shape[0])
        error_orbital = np.interp(np.mod(time, xg.max()), xg,
                                  self.yg.compute())
        return {
            "simulated_error_orbital":
            x_ac * error_orbital[:, np.newaxis] * self.conversion_factor,
        }
