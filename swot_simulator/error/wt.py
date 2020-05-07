from typing import Optional
import numpy as np
import scipy.ndimage.filters
from . import utils

from .. import VOLUMETRIC_MEAN_RADIUS


class ErrorStat:
    def __init__(self, delta_al: float) -> None:
        # - Define power spectrum of error in path delay
        #   due to wet tropo
        f = np.arange(1. / 3000., 1. / float(2. * delta_al), 1. / 3000.)
        # - Global mean wet tropo power spectrum in cm**2/(cycle/km)
        #   for L >= 100 km
        pswt = 3.156 * 10**-5 * f**(-8. / 3.)  # *10**-4
        # - Wet tropo power spectrum in cm**2/(cycle/km) for L < 100 km
        indf = np.where(f > 10**-2)
        pswt[indf] = 1.4875 * 10**-4 * f[indf]**(-2.33)  # *10**-4
        self.pswt = pswt
        self.freq = f
        self.radio_r: Optional[np.ndarray] = None
        self.radio_l: Optional[np.ndarray] = None
        self.wet_tropo1: Optional[np.ndarray] = None
        self.wet_tropo2: Optional[np.ndarray] = None
        self.wt: Optional[np.ndarray] = None
        self.wet_tropo1nadir: Optional[np.ndarray] = None
        self.wet_tropo2nadir: Optional[np.ndarray] = None
        self.wtnadir: Optional[np.ndarray] = None

    def init_radio(self,
                   x_al: np.ndarray,
                   dal: float,
                   len_repeat: float,
                   nseed: int = 0) -> None:
        # - Define radiometer error power spectrum for a beam
        #   High frequencies are cut to filter the associated error:
        #   during the reconstruction of the wet trop signal
        psradio = 9.5 * 10**-5 * self.freq**-1.79
        psradio[np.where(
            (self.freq < 1. / 1000.))] = 9.5 * 10**-5 * (10**-3)**-1.79
        indf = np.where((self.freq > 0.0023) & (self.freq <= 0.0683))
        psradio[indf] = 0.036 * self.freq[indf]**-0.814
        psradio[np.where(self.freq > 0.0683)] = 0.32
        # - Compute random coefficients (1D) for the radiometer error
        #   power spectrum for right and left beams
        _hrad = utils.gen_signal1d(self.freq,
                                   psradio,
                                   x_al,
                                   fmin=1 / len_repeat,
                                   fmax=1 / (2 * dal),
                                   alpha=10,
                                   nseed=nseed + 100,
                                   hf_extpl=True,
                                   lf_extpl=True)
        self.radio_r = _hrad * 10**(-2)
        _hrad = utils.gen_signal1d(self.freq,
                                   psradio,
                                   x_al,
                                   fmin=1 / len_repeat,
                                   fmax=1 / (2 * dal),
                                   alpha=10,
                                   nseed=nseed + 200,
                                   hf_extpl=True,
                                   lf_extpl=True)
        self.radio_l = _hrad * 10**(-2)

    def make_error(self,
                   x_al: np.array,
                   x_ac: np.array,
                   dal: float,
                   dac: float,
                   nbeam: int,
                   sigma: float,
                   beam_pos: list,
                   sat_const: dict,
                   len_repeat,
                   nseed: int = 0,
                   lac_max: float = 500) -> None:
        if nbeam not in [1, 2, 12]:
            raise ValueError(
                f"nbeam should be either 1 or 2 or 12, not {nbeam}")

        fka = sat_const['Fka']
        sat_elev = sat_const['height']
        nal = np.shape(x_al)[0]
        nac = np.shape(x_ac)[0]
        _to_km = (1 / (fka * 2 * np.pi / sat_const['C'] * sat_const['B']) *
                  (1 + sat_elev / VOLUMETRIC_MEAN_RADIUS) * np.pi / 180 *
                  10**3)
        # - Initialization of radiometer error in right and left beam
        self.init_radio(x_al, dal, len_repeat, nseed=nseed)
        # - Initialization of swath matrices and large swath matrices
        #   (which include wet tropo data around the nadir and outside
        #   the swath)
        #   x_ac_large and wt_large are necessary to compute the gaussian
        #   footprint of a beam on the nadir or near the edge of the swath
        start_x = -2 * sigma / float(dac) + x_ac[0]
        stop_x = (2 * sigma / float(dac) + x_ac[-1] + dac)
        x_ac_large = np.arange(start_x, stop_x, dac)
        naclarge = np.shape(x_ac_large)[0]
        # - Compute path delay error due to wet tropo and radiometer error
        #   using random coefficient initialized with power spectrums
        self.wet_tropo1 = np.zeros((nal, nac))
        self.wet_tropo2 = np.zeros((nal, nac))
        fminx = 1. / len_repeat
        fminy = 1. / lac_max
        wt = utils.gen_signal2d_rectangle(self.freq,
                                          self.pswt,
                                          x_al,
                                          x_ac,
                                          fminx=fminx,
                                          fminy=fminy,
                                          fmax=0.05,
                                          alpha=10,
                                          nseed=nseed,
                                          lf_extpl=True)
        self.wt = wt.transpose() * 10**(-2)
        wt_large = utils.gen_signal2d_rectangle(self.freq,
                                                self.pswt,
                                                x_al,
                                                x_ac_large,
                                                fminx=fminx,
                                                fminy=fminy,
                                                fmax=0.05,
                                                alpha=10,
                                                nseed=nseed,
                                                lf_extpl=True)
        wt_large = wt_large.transpose() * 10**(-2)
        # - Compute Residual path delay error after a 1-beam radiometer
        #   correction
        if nbeam in [1, 12]:
            beam = np.zeros((nal))
            diff_h1 = np.zeros((nal, nac))
            # - Find across track indices in the gaussian footprint of
            #   2. * sigma
            indac = np.where((x_ac_large < 2 * sigma)
                             & (x_ac_large > -2 * sigma))[0]
            for idx in range(0, nal):
                # - Find along track indices in the gaussian footprint of
                #   2.*p.sigma
                delta_x_al = x_al[:] - x_al[idx]
                indal = np.where((delta_x_al <= (2 * sigma))
                                 & (delta_x_al > (-2 * sigma)))[0]
                slice_al = slice(min(indal), max(indal) + 1)
                slice_ac = slice(min(indac), max(indac) + 1)
                x, y = np.meshgrid(x_ac_large[slice_al],
                                   x_al[slice_ac] - x_al[idx])
                # - Compute path delay on gaussian footprint
                g = 1 / (2 * np.pi * sigma**2) * np.exp(-(x**2 + y**2) /
                                                        (2 * sigma**2))
                beam[idx] = (np.sum(np.sum(g * wt_large[slice_al, slice_ac])) /
                             np.sum(np.sum(g)) + self.radio_l[idx])
            # - Filtering beam signal to cut frequencies higher than 125 km
            beam = scipy.ndimage.filters.gaussian_filter(beam, 30. / dal)
            beam2d = np.array(nac * [beam]).T
            # - Compute residual path delay
            diff_h1 = self.wt - beam2d
            self.wet_tropo1 = +diff_h1  # en 2d
            self.wet_tropo1nadir = wt_large[:, naclarge // 2] - beam[:]
        # - Compute Residual path delay error after a 2-beams radiometer
        #   correction
        if nbeam in [2, 12]:
            beam_r = np.zeros((nal))
            beam_l = np.zeros((nal))
            diff_h2 = np.zeros((nal, nac))
            diff_h2nadir = np.zeros((nal))
            # - Find righ and leftacross track indices in the gaussian
            #   footprint of 2.*p.sigma
            ind_r = x_ac_large + beam_pos[1]
            indac_r = np.where((ind_r < 2 * sigma) & (ind_r > -2 * sigma))[0]
            ind_l = x_ac_large + beam_pos[0]
            indac_l = np.where((ind_l < 2 * sigma) & (ind_l > -2 * sigma))[0]
            for idx in range(nal):
                # - Find along track indices in the gaussian footprint
                #   of 2.*p.sigma
                delta_x_al = x_al[:] - x_al[idx]
                indal = np.where((delta_x_al <= (2 * sigma))
                                 & (delta_x_al > (-2 * sigma)))[0]
                slice_al = slice(min(indal), max(indal) + 1)
                slice_acr = slice(min(indac_r), max(indac_r) + 1)
                slice_acl = slice(min(indac_l), max(indac_l) + 1)
                x, y = np.meshgrid(x_ac_large[slice_acr],
                                   x_al[slice_al] - x_al[idx])
                # - Compute path delay on left and right gaussian footprint
                g = 1. / (2. * np.pi * sigma**2) * np.exp(-(x**2. + y**2.) /
                                                          (2. * sigma**2))
                beam_r[idx] = (
                    np.sum(np.sum(g * wt_large[slice_al, slice_acr])) /
                    np.sum(np.sum(g)) + self.radio_r[idx])
                x, y = np.meshgrid(x_ac_large[slice_acl],
                                   x_al[slice_al] - x_al[idx])
                g = 1. / (2. * np.pi * sigma**2) * np.exp(-(x**2. + y**2.) /
                                                          (2. * sigma**2))
                beam_l[idx] = (
                    np.sum(np.sum(g * wt_large[slice_al, slice_acl])) /
                    np.sum(np.sum(g)) + self.radio_l[idx])
            # - Filtering beam signal to cut frequencies higher than 125 km
            beam_r = scipy.ndimage.filters.gaussian_filter(beam_r, 30 / dal)
            beam_l = scipy.ndimage.filters.gaussian_filter(beam_l, 30 / dal)
            # - Compute residual path delay (linear combination of left
            #   and right path delay)
            #pol = np.polyfit([beam_pos[0], beam_pos[1]],
            #                    [beam_l[i], beam_r[i]], 1)
            polyfit = np.polynomial.polynomial.polyfit
            pol = polyfit([beam_pos[0], beam_pos[1]], [beam_l, beam_r], 1)
            beam = (np.array(nac * [pol[0]]).T +
                    np.array(nal * [x_ac]) * np.array(nac * [pol[1]]).T)
            #pbeam = np.polyval(pol, -x_ac)
            diff_h2 = self.wt - beam
            diff_h2nadir = wt_large[:, naclarge // 2] - beam[:, nac // 2]
            self.wet_tropo2 = diff_h2  # en 2d
            self.wet_tropo2nadir = diff_h2nadir  # en 1d

        self.wtnadir = wt_large[:, naclarge // 2]

    def reconstruct_2d_error(self) -> np.array:
        return self.wt, self.wet_tropo2, self.wet_tropo1,

    def reconstruct_nadir(self) -> np.array:
        return self.wtnadir, self.wet_tropo2nadir, self.wet_tropo1nadir
