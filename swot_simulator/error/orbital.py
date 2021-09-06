# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Orbital error
-------------
"""
from typing import NamedTuple, Tuple
import logging

#
import numpy as np

#
from .. import random_signal

#: Module logger
LOGGER = logging.getLogger(__name__)

#: Signal amplitude of the orbital error in micro-radians
AMPLITUDE = 100

#: Delta T of the spatial sampling in seconds
DT = 60


class OrbitalErrorSpectrum(NamedTuple):
    """Orbital Error Spectrum"""
    yg: np.ndarray
    fmaxr: float


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
    ps_orbital = np.exp(-0.5 *
                        (spatial_frequency - orbital_frequency)**2 /
                        sigma_peak**2)
    ps_orbital[ps_orbital < 1 / 1000] = 0.
    ps_orbital /= np.sum(ps_orbital * df)
    ps_orbital *= AMPLITUDE**2
    return random_signal.gen_psd_1d(spatial_frequency,
                                    ps_orbital,
                                    rng,
                                    alpha=10)


class Model:
    """
    Simulate the orbital error

    Args:
        orbit_duration (np.timedelta64): Orbit duration
        rng (np.random.Generator): Random number generator
    """
    def __init__(self, orbit_duration: np.timedelta64,
                 rng: np.random.Generator) -> None:
        LOGGER.info("Generating orbital error spectrum")
        yg, self.fmaxr = _orbital_error_spectrum(orbit_duration, rng)
        self.yg = yg.astype(np.float32)

    def generate(
        self,
        time: np.ndarray,
    ) -> np.ndarray:
        """Generate orbital error

        Args:
            time (np.ndarray): time vector

        Returns:
            np.ndarray: orbital error
        """
        time = time.astype("datetime64[us]").astype("float64") * 1e-6
        xg = np.linspace(0, 0.5 / self.fmaxr * self.yg.shape[0],
                         self.yg.shape[0])
        return np.interp(np.mod(time, xg.max()), xg, self.yg)
