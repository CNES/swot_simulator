import numpy as np
from . import utils
import sys
import logging
from typing import IO, Optional
from scipy.ndimage.filters import gaussian_filter
import warnings
warnings.simplefilter('ignore', np.RankWarning)


# Define logging level for debug purposes
logger = logging.getLogger(__name__)

VOLUMETRIC_MEAN_RADIUS = 6371

class error_stat():
    def __init__(self, dal: int) -> None:
        # - Define power spectrum of error in path delay
        #   due to wet tropo
        f = np.arange(1./3000., 1./float(2. * dal), 1./3000.)
        # - Global mean wet tropo power spectrum in cm**2/(cycle/km)
        #   for L >= 100 km
        PSwt = 3.156 * 10**-5 * f**(-8./3.)  # *10**-4
        # - Wet tropo power spectrum in cm**2/(cycle/km) for L < 100 km
        indf = np.where(f > 10**-2)
        PSwt[indf] = 1.4875 * 10**-4 * f[indf]**(-2.33)  # *10**-4
        self.PSwt = PSwt
        self.freq = f


    def init_radio(self, x_al:np.ndarray, dal: float, len_repeat: float,
                             nseed: Optional[int]=0)-> None:
        # - Compute random coefficients in 2D using the previously
        #   defined power spectrum
        gencoef = utils.gen_coeff_signal2d(self.freq, self.PSwt, ncomp2d,
                                           nseed + 100)
        self.A_wt, self.phi_wt, self.frx_wt, self.fry_wt = gencoef
        # - Define radiometer error power spectrum for a beam
        #   High frequencies are cut to filter the associated error:
        #   during the reconstruction of the wet trop signal
        PSradio = 9.5 * 10**-5 * self.freq**-1.79
        PSradio[np.where((self.freq < 1./1000.))] = 9.5 * 10**-5 * (10**-3)**-1.79
        indf = np.where((self.freq > 0.0023) & (self.freq <= 0.0683))
        PSradio[indf] = 0.036 * self.freq[indf]**-0.814
        PSradio[np.where(self.freq > 0.0683)] = 0.32
        # - Compute random coefficients (1D) for the radiometer error
        #   power spectrum for right and left beams
        _hrad = gen_signal1d(self.freq, PSradio, x_al, fmin=1./len_repeat,
                             fmax=1./(2*dal), alpha=10, seed=nseed + 100,
                             hf_extpl=True, lf_extpl=True)
        self.radio_r = _hrad
        _hrad = gen_signal1d(self.freq, PSradio, x_al, fmin=1./len_repeat,
                             fmax=1./(2*dal), alpha=10, seed=nseed + 200,
                             hf_extpl=True, lf_extpl=True)
        self.radio_l = _hrad


    def make_error(self, x_al: np.array, x_ac: np.array, dal: float,
                   dac:float, nbeam:int, sigma:float,
                   beam_pos: list, sat_const: dict, len_repeat,
                   nseed: Optional[int]=0,
                   lac_max: Optional[float]=500) -> None:
        Rearth = VOLUMETRIC_MEAN_RADIUS
        Fka = sat_const['Fka']
        sat_elev = sat_const['height']
        nal = np.shape(x_al)[0]
        nac = np.shape(x_ac)[0]
        _to_km = (1 / (Fka * 2 * np.pi / sat_const['C'] * sat_const['B'])
                  * (1 + sat_elev/Rearth)*np.pi/180. * 10**3)
        # - Initialization of radiometer error in right and left beam
        self.init_radio(self, x_al, dal, len_repeat, nseed=nseed)
        # - Initialization of swath matrices and large swath matrices
        #   (which include wet tropo data around the nadir and outside
        #   the swath)
        #   x_ac_large and wt_large are necessary to compute the gaussian
        #   footprint of a beam on the nadir or near the edge of the swath
        start_x = -2. * sigma / float(dac) + x_ac[0]
        stop_x = (2. * sigma / float(dac) +x_ac[-1] + dac)
        x_ac_large = np.arange(start_x, stop_x, dac)
        naclarge = np.shape(x_ac_large)[0]
        # - Compute path delay error due to wet tropo and radiometer error
        #   using random coefficient initialized with power spectrums
        self.wet_tropo1 = np.zeros((nal, nac))
        self.wet_tropo2 = np.zeros((nal, nac))
        fminx = 1./len_repeat
        fminy = 1./lac_max
        wt = utils.gen_signal2d_rectangle(self.freq, self.PSwt, x_al, x_ac,
                                            fminx=fminx, fminy=fminy,
                                            fmax=1./20, alpha=10, nseed=nseed,
                                            lf_extpl=True)
        wt_large = utils.gen_signal2d_rectangle(self.freq, self.PSwt, x_al,
                                                x_ac_large, fminx=fminx,
                                                fminy=fminy, fmax=1./20,
                                                alpha=10, nseed=nseed,
                                                lf_extpl=True)
        # - Compute Residual path delay error after a 1-beam radiometer
        #   correction
        if nbeam == 1 or nbeam == 12:
            beam = np.zeros((nal))
            diff_h1 = np.zeros((nal, nac))
            # - Find across track indices in the gaussian footprint of
            #   2. * sigma
            indac = np.where((x_ac_large < 2.*sigma)
                             & (x_ac_large > -2.*sigma))[0]
            for i in range(0, nal):
                # - Find along track indices in the gaussian footprint of
                #   2.*p.sigma
                delta_x_al = x_al[:] - x_al[i]
                indal = np.where((delta_x_al <= (2 * sigma))
                                    & (delta_x_al > (-2 * sigma)))[0]
                slice_al = slice(min(indal), max(indal) + 1)
                slice_ac = slice(min(indac), max(indac) + 1)
                x, y = np.meshgrid(x_ac_large[slice_al],
                                   x_al[slice_ac] - x_al[i])
                # - Compute path delay on gaussian footprint
                G = 1./(2.*np.pi*sigma**2) * np.exp(-(x**2. + y**2.)
                                                      / (2.*sigma**2))
                beam[i] = (sum(sum(G * wt_large[slice_al, slice_ac]))
                           / sum(sum(G)) + self.radio_l[i])
            # - Filtering beam signal to cut frequencies higher than 125 km
            beam = gaussian_filter(beam, 30. / dal)
            beam2d = np.array(nac*[beam]).T
            # - Compute residual path delay
            diff_h1 = self.wt - beam2d
            self.wet_tropo1 = + diff_h1  # en 2d
            self.wet_tropo1nadir = wt_large[:, int(naclarge/2.)] - beam[:]
        # - Compute Residual path delay error after a 2-beams radiometer
        #   correction
        if nbeam == 2 or nbeam == 12:
            beam_r = np.zeros((nal))
            beam_l = np.zeros((nal))
            diff_h2 = np.zeros((nal, nac))
            diff_h2nadir = np.zeros((nal))
            # - Find righ and leftacross track indices in the gaussian
            #   footprint of 2.*p.sigma
            ind_r = x_ac_large + beam_pos[1]
            indac_r = np.where((ind_r < 2.*sigma) & (ind_r > -2.*sigma))[0]
            ind_l = x_ac_large + beam_pos[0]
            indac_l = np.where((ind_l < 2.*sigma) & (ind_l > -2.*sigma))[0]
            for i in range(nal):
                # - Find along track indices in the gaussian footprint
                #   of 2.*p.sigma
                delta_x_al = x_al[:] - x_al[i]
                indal = np.where((delta_x_al <= (2*sigma))
                                 & (delta_x_al > (-2*sigma)))[0]
                slice_al = slice(min(indal), max(indal) + 1)
                slice_acr = slice(min(indac_r), max(indac_r) + 1)
                slice_acl = slice(min(indac_l), max(indac_l) + 1)
                x, y = np.meshgrid(x_ac_large[slice_acr],
                                      x_al[slice_al] - x_al[i])
                # - Compute path delay on left and right gaussian footprint
                G = 1. / (2.*np.pi*sigma**2) * np.exp(-(x**2. + y**2.)
                                                        / (2.*sigma**2))
                beam_r[i] = (sum(sum(G*wt_large[slice_al, slice_acr]))
                             / sum(sum(G))+self.radio_r[i])
                beam_l[i] = (sum(sum(G*wt_large[slice_al, slice_acl]))
                             / sum(sum(G)) + self.radio_l[i])
            # - Filtering beam signal to cut frequencies higher than 125 km
            beam_r = gaussian_filter(beam_r, 30. / dal)
            beam_l = gaussian_filter(beam_l, 30. / dal)
                # - Compute residual path delay (linear combination of left
                #   and right path delay)
                #pol = np.polyfit([beam_pos[0], beam_pos[1]],
                #                    [beam_l[i], beam_r[i]], 1)
            polyfit = np.polynomial.polynomial.polyfit
            pol = polyfit([beam_pos[0], beam_pos[1]], [beam_l, beam_r], 1)
            beam = (np.array(nac*[pol[0]]).T
                    + np.array(nal*[x_ac]) * np.array(nac*[pol[1]]).T)
            #pbeam = np.polyval(pol, -x_ac)
            diff_h2 = self.wt - beam
            diff_h2nadir = wt_large[:, int(naclarge/2.)] - beam[:,int(nac/2.)]
            self.wet_tropo2 = + diff_h2  # en 2d
            self.wet_tropo2nadir = + diff_h2nadir  # en 1d
        self.wtnadir = + wt_large[:, int(naclarge/2.)]
        self.wt = wt
        if (not nbeam == 1) and (not nbeam == 2) \
           and (not nbeam == 12):
            logger.error("\n nbeam = {} \n".format(nbeam))
            logger.error("wrong number of beam, nbeam should be either 1"
                         "or 2 or 12")
            sys.exit(1)

    def reconstruct_2D_error(self)-> np.array:
        return self.wt, self.wet_tropo2, self.wet_tropo1,

    def reconstruct_nadir(self)-> np.array:
        return self.wtnadir, self.wet_tropo2nadir, self.wet_tropo1nadir

