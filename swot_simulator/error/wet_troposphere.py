# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Wet troposphere errors
----------------------
"""
from typing import Dict, List, Tuple
import numba as nb
import numba.typed
import numpy as np
import scipy.ndimage.filters

from .. import random_signal
from .. import settings
from .. import F_KA, VOLUMETRIC_MEAN_RADIUS, CELERITY, BASELINE


@nb.njit(cache=True, nogil=True)
def _meshgrid(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    nx = x.size
    ny = y.size
    mx = np.empty((ny, nx))
    my = np.empty((ny, nx))
    for idx in range(ny):
        mx[idx, :] = x
    for idx in range(nx):
        my[:, idx] = y
    return mx, my


@nb.njit(cache=True, nogil=True)
def _calculate_path_delay_lr(beam_positions: List[float], sigma: float,
                             radio_r: np.ndarray, radio_l: np.ndarray,
                             x_al: np.ndarray, x_ac_large: np.ndarray,
                             wt_large: np.ndarray):
    beam_r = np.empty((x_al.shape[0], ))
    beam_l = np.empty((x_al.shape[0], ))
    # Find righ and leftacross track indices in the gaussian
    # footprint of 2.*p.sigma
    ind_r = x_ac_large + beam_positions[1]
    indac_r = np.where((ind_r < 2 * sigma) & (ind_r > -2 * sigma))[0]
    ind_l = x_ac_large + beam_positions[0]
    indac_l = np.where((ind_l < 2 * sigma) & (ind_l > -2 * sigma))[0]

    factor = 1 / (2 * np.pi * sigma**2)

    for idx, xal in enumerate(x_al):
        # Find along track indices in the gaussian footprint
        # of 2.*p.sigma
        delta_x_al = x_al - xal
        indal = np.where((delta_x_al <= (2 * sigma))
                         & (delta_x_al > (-2 * sigma)))[0]
        slice_al = slice(indal[0], indal[-1] + 1)
        slice_acr = slice(indac_r[0], indac_r[-1] + 1)
        slice_acl = slice(indac_l[0], indac_l[-1] + 1)
        x, y = _meshgrid(x_ac_large[slice_acr], x_al[slice_al] - xal)
        # Compute path delay on left and right gaussian footprint
        g = factor * np.exp(-(x**2 + y**2) / (2 * sigma**2))
        beam_r[idx] = np.sum(
            g * wt_large[slice_al, slice_acr]) / np.sum(g) + radio_r[idx]
        x, y = _meshgrid(x_ac_large[slice_acl], x_al[slice_al] - xal)
        g = factor * np.exp(-(x**2 + y**2) / (2 * sigma**2))
        beam_l[idx] = np.sum(
            g * wt_large[slice_al, slice_acl]) / np.sum(g) + radio_l[idx]
    return beam_r, beam_l


@nb.njit(cache=True, nogil=True)
def _calculate_path_delay(sigma: float, radio: np.ndarray, x_al: np.ndarray,
                          x_ac_large: np.ndarray, wt_large: np.ndarray):
    beam = np.empty((x_al.shape[0], ))
    # Find across track indices in the gaussian footprint of
    # 2. * sigma
    indac = np.where((x_ac_large < 2 * sigma) & (x_ac_large > -2 * sigma))[0]
    factor = 1 / (2 * np.pi * sigma**2)
    for idx, xal in enumerate(x_al):
        delta_x_al = x_al[:] - xal
        indal = np.where((delta_x_al <= (2 * sigma))
                         & (delta_x_al > (-2 * sigma)))[0]
        slice_al = slice(indal[0], indal[-1] + 1)
        slice_ac = slice(indac[0], indac[-1] + 1)
        x, y = _meshgrid(x_ac_large[slice_ac], x_al[slice_al] - xal)
        # Compute path delay on gaussian footprint
        g = factor * np.exp(-(x**2 + y**2) / (2 * sigma**2))
        beam[idx] = np.sum(
            g * wt_large[slice_al, slice_ac]) / np.sum(g) + radio[idx]
    return beam


class WetTroposphere:
    """Wet troposphere errors

    Args:
        parameters (settings.Parameters): Simulation settings
    """
    ALPHA = 10
    LC_MAX = 500
    F_MAX = 0.05

    def __init__(self, parameters: settings.Parameters) -> None:
        # Store the generation parameters of the random signal.
        self.beam_positions = parameters.beam_position
        self.delta_ac = parameters.delta_ac
        self.delta_al = parameters.delta_al
        self.len_repeat = parameters.len_repeat
        self.nbeam = parameters.nbeam
        self.rng = parameters.rng()
        self.rng_radio_l = parameters.rng()
        self.rng_radio_r = parameters.rng()
        self.sigma = parameters.sigma
        # TODO
        self.conversion_factor = (
            1 / (F_KA * 2 * np.pi / CELERITY * BASELINE) *
            (1 + (parameters.height * 1e-3) / VOLUMETRIC_MEAN_RADIUS) * np.pi /
            180 * 1e3)
        # Define power spectrum of error in path delay due to wet tropo
        freq = np.arange(1 / 3000, 1 / (2 * self.delta_al), 1 / 3000)
        # Global mean wet tropo power spectrum in cm**2/(cycle/km)
        # for L >= 100 km
        pswt = 3.156 * 1e-05 * freq**(-8 / 3)
        # Wet tropo power spectrum in cm**2/(cycle/km) for L < 100 km
        mask = freq > 1e-2
        pswt[mask] = 1.4875 * 1e-4 * freq[mask]**(-2.33)
        self.pswt = pswt
        self.freq = freq
        self.fminx = 1 / self.len_repeat
        self.ps2d, self.f = random_signal.gen_ps2d(freq,
                                                   pswt,
                                                   fminx=self.fminx,
                                                   fminy=1 / self.LC_MAX,
                                                   fmax=self.F_MAX,
                                                   alpha=self.ALPHA,
                                                   lf_extpl=True,
                                                   hf_extpl=True)

    def _radiometer_error(self,
                          x_al: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        # Define radiometer error power spectrum for a beam
        # High frequencies are cut to filter the associated error:
        # during the reconstruction of the wet trop signal
        psradio = 9.5 * 1e-5 * self.freq**-1.79
        psradio[self.freq < 1e-3] = 9.5 * 1e-5 * (1e-3)**-1.79
        mask = (self.freq > 0.0023) & (self.freq <= 0.0683)
        psradio[mask] = 0.036 * self.freq[mask]**-0.814
        psradio[self.freq > 0.0683] = 0.32
        # Compute random coefficients (1D) for the radiometer error
        # power spectrum for right and left beams
        hrad = random_signal.gen_signal_1d(self.freq,
                                           psradio,
                                           x_al,
                                           fmin=1 / self.len_repeat,
                                           fmax=1 / (2 * self.delta_al),
                                           alpha=10,
                                           rng=self.rng_radio_r,
                                           hf_extpl=True,
                                           lf_extpl=True)
        radio_r = hrad * 1e-2
        hrad = random_signal.gen_signal_1d(self.freq,
                                           psradio,
                                           x_al,
                                           fmin=1 / self.len_repeat,
                                           fmax=1 / (2 * self.delta_al),
                                           alpha=10,
                                           rng=self.rng_radio_l,
                                           hf_extpl=True,
                                           lf_extpl=True)
        radio_l = hrad * 1e-2

        return radio_r, radio_l

    def generate(self, x_al: np.array,
                 x_ac: np.array) -> Dict[str, np.ndarray]:
        """
        Generate wet troposphere errors

        Args:
            x_al (numpy.ndarray): Along track distance
            x_ac (numpy.ndarray): Across track distance

        Returns:
            dict: variable name and errors simulated.
        """
        num_lines = x_al.shape[0]
        num_pixels = x_ac.shape[0]

        # Initialization of radiometer error in right and left beam
        radio_r, radio_l = self._radiometer_error(x_al)
        # Initialization of swath matrices and large swath matrices (which
        # include wet tropo data around the nadir and outside the swath)
        start_x = -2 * self.sigma / self.delta_ac + x_ac[0]
        stop_x = (2 * self.sigma / self.delta_ac + x_ac[-1] + self.delta_ac)
        # x_ac_large and wt_large are necessary to compute the gaussian
        # footprint of a beam on the nadir or near the edge of the swath
        x_ac_large = np.arange(start_x, stop_x, self.delta_ac)
        naclarge = np.shape(x_ac_large)[0]
        # Compute path delay error due to wet tropo and radiometer error
        # using random coefficient initialized with power spectrums
        wt = random_signal.gen_signal_2d_rectangle(self.ps2d,
                                                   self.f,
                                                   x_al,
                                                   x_ac,
                                                   fminx=self.fminx,
                                                   fminy=1 / self.LC_MAX,
                                                   fmax=self.F_MAX,
                                                   alpha=self.ALPHA,
                                                   rng=self.rng)
        wt = wt.T * 1e-2
        wt_large = random_signal.gen_signal_2d_rectangle(self.ps2d,
                                                         self.f,
                                                         x_al,
                                                         x_ac_large,
                                                         fminx=self.fminx,
                                                         fminy=1 / self.LC_MAX,
                                                         fmax=self.F_MAX,
                                                         alpha=self.ALPHA,
                                                         rng=self.rng)
        wt_large = wt_large.T * 1e-2

        # Compute Residual path delay error after a 1-beam radiometer
        # correction
        if self.nbeam == 1:
            beam = _calculate_path_delay(self.sigma, radio_l, x_al, x_ac_large,
                                         wt_large)
            beam = scipy.ndimage.filters.gaussian_filter(
                beam, 30. / self.delta_al)
            beam2d = np.vstack(num_pixels * (beam, )).T
            # Compute residual path delay
            wet_tropo = wt - beam2d
            wet_tropo_nadir = wt_large[:, naclarge // 2] - beam

        # Compute Residual path delay error after a 2-beams radiometer
        # correction
        elif self.nbeam == 2:
            beam_r, beam_l = _calculate_path_delay_lr(
                numba.typed.List(self.beam_positions), self.sigma, radio_r,
                radio_l, x_al, x_ac_large, wt_large)
            # Filtering beam signal to cut frequencies higher than 125 km
            beam_r = scipy.ndimage.filters.gaussian_filter(
                beam_r, 30 / self.delta_al)
            beam_l = scipy.ndimage.filters.gaussian_filter(
                beam_l, 30 / self.delta_al)
            # Compute residual path delay (linear combination of left
            # and right path delay)
            polyfit = np.polynomial.polynomial.polyfit
            pol = polyfit([self.beam_positions[0], self.beam_positions[1]],
                          [beam_l, beam_r], 1)
            beam = (np.array(num_pixels * [pol[0]]).T +
                    np.array(num_lines * [x_ac]) *
                    np.array(num_pixels * [pol[1]]).T)
            wet_tropo = wt - beam
            wet_tropo_nadir = wt_large[:, naclarge //
                                       2] - beam[:, num_pixels // 2]
        else:
            raise ValueError("nbeam must be in [1, 2]")

        # wt_nadir = wt_large[:, naclarge // 2]

        return {
            "simulated_error_troposphere": wet_tropo,
            "simulated_error_troposphere_nadir": wet_tropo_nadir
        }
