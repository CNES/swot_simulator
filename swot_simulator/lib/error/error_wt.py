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
    def __init__(self, dal: int, ncomp2d: int) -> None:
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


    def init_error_gensignal(self, ncomp2d: int)-> None:
        # - Compute random coefficients in 2D using the previously
        #   defined power spectrum
        gencoef = utils.gen_coeff_signal2d(self.freq, self.PSwt, ncomp2d)
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
        gencoef = utils.gen_coeff_signal1d(self.freq, PSradio, ncomp2d)
        self.A_radio_r, self.phi_radio_r, self.fr_radio_r = gencoef
        gencoef = utils.gen_coeff_signal1d(self.freq, PSradio, ncomp2d)
        self.A_radio_l, self.phi_radio_l, self.fr_radio_l = gencoef


    def make_error(self, x_al: np.array, x_ac: np.array, al_cycle: float,
                   cycle:int, dal: float, dac:float, nbeam:int, sigma:float,
                   beam_pos: list, sat_const: dict) -> None:
        Rearth = VOLUMETRIC_MEAN_RADIUS
        Fka = sat_const['Fka']
        sat_elev = sat_const['height']
        nal = np.shape(x_al)[0]
        nac = np.shape(x_ac)[0]
        _to_km = (1 / (Fka * 2 * np.pi / sat_const['C'] * sat_const['B'])
                  * (1 + sat_elev/Rearth)*np.pi/180. * 10**3)
        # - Initialization of radiometer error in right and left beam
        err_radio_r = np.zeros((nal))
        err_radio_l = np.zeros((nal))
        # - Initialization of swath matrices and large swath matrices
        #   (which include wet tropo data around the nadir and outside
        #   the swath)
        #   x_ac_large and wt_large are necessary to compute the gaussian
        #   footprint of a beam on the nadir or near the edge of the swath
        x, y = np.meshgrid(x_ac, x_al + float(cycle * al_cycle))
        start_x = -2. * sigma / float(dac) + x_ac[0]
        stop_x = (2. * sigma / float(dac) +x_ac[-1] + dac)
        x_ac_large = np.arange(start_x, stop_x, dac)
        naclarge = np.shape(x_ac_large)[0]
        wt_large = np.zeros((nal, naclarge))
        y_al = np.float64(x_al) + float(cycle * al_cycle)
        x_large, y_large = np.meshgrid(x_ac_large, y_al)
        ncomp2d = np.shape(self.frx_wt)[0]
        # - Compute path delay error due to wet tropo and radiometer error
        #   using random coefficient initialized with power spectrums
        self.wt = np.zeros((nal, nac))
        self.wet_tropo1 = np.zeros((nal, nac))
        self.wet_tropo2 = np.zeros((nal, nac))
        for comp in range(ncomp2d):
            phase_x_al = (2. * np.pi * (float(self.frx_wt[comp]) * x
                          + float(self.fry_wt[comp]) * y)) % (2.*np.pi)
            self.wt = (self.wt + self.A_wt[comp]
                       * np.cos(phase_x_al + self.phi_wt[comp])*10**-2)
            #import pdb; pdb.set_trace()
            phase_x_al_large = (2. * np.pi * (float(self.frx_wt[comp])
                                * x_large + float(self.fry_wt[comp])
                                * y_large)) % (2.*np.pi)
            wt_large = (wt_large + self.A_wt[comp] * np.cos(phase_x_al_large
                        + self.phi_wt[comp])*10**-2)
            phase_x_al = (2. * np.pi * float(self.fr_radio_r[comp])
                          * (x_al + float(cycle * al_cycle))) % (2.*np.pi)
            err_radio_r = (err_radio_r + 2*self.A_radio_r[comp]
                           * np.cos(phase_x_al+ self.phi_radio_r[comp])*10**-2)
            phase_x_al = (2. * np.pi * float(self.fr_radio_l[comp])
                          * (x_al + float(cycle * al_cycle))) % (2.*np.pi)
            err_radio_l = (err_radio_l + 2*self.A_radio_l[comp]
                           * np.cos(phase_x_al+ self.phi_radio_l[comp])*10**-2)
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
                           / sum(sum(G)) + err_radio_l[i])
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
                             / sum(sum(G))+err_radio_r[i])
                beam_l[i] = (sum(sum(G*wt_large[slice_al, slice_acl]))
                             / sum(sum(G)) + err_radio_l[i])
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

"""

class errorwtnadir():
    '''Class errornadir defines the error on the nadir.
    Random realisation of errors can be initialized using init_error.
    If the random realisations have already been computed and stored in file
    file_coeff, the random realisations are read directly using load_coeff.
    The correspondg errors on a swath can be computed using make_error. '''
    def __init__(self):
        return

    def init_error(self, parameters):
        '''Initialization of errors: Random realisation of errors are computed
        using a known power spectrum.
        The outputs are the amplitude, the phase and the frequency of each
        random realisation.
        By default, there are ncomp2d=2000 random realisations for the
        wet tropo and ncomp1d=2000 random realisations for the nadir 1d
        spectrum error.'''
        # - Define power spectrum of error in path delay due to wet tropo
        f = numpy.arange(1./3000., 1./float(2.*parameters.delta_al), 1./3000.)
        # - Global mean wet tropo power spectrum in cm**2/(cycle/km)
        #   for L >= 100 km
        PSwt = 3.156 * 10**-5 * f**(-8./3.)  # *10**-4
        # - Wet tropo power spectrum in cm**2/(cycle/km) for L < 100 km
        indf = numpy.where(f > 10**-2)
        PSwt[indf] = 1.4875 * 10**-4 * f[indf]**(-2.33)  # *10**-4
        # - Compute random coefficients in 2D using the previously defined
        #   power spectrum
        gencoef = mod_tools.gen_coeff_signal2d(f, PSwt, parameters.ncomp2d)
        self.A_wt, self.phi_wt, self.frx_wt, self.fry_wt = gencoef
        # - Define radiometer error power spectrum for a beam
        #   High frequencies are cut to filter the associated error during
        #   the reconstruction of the wet trop signal
        # f=numpy.arange(1./3000.,1./float(20.),1./3000.)
        PSradio = 9.5 * 10**-5 * f**-1.79
        PSradio[numpy.where((f < 1./1000.))] = 9.5 * 10**-5*(10**-3)**-1.79
        indf = numpy.where((f > 0.0023) & (f <= 0.0683))
        PSradio[indf] = 0.036 * f[indf]**-0.814
        PSradio[numpy.where(f > 0.0683)] = 0.32
        # - Compute random coefficients (1D) for the radiometer error power
        #   spectrum for right and left beams
        gencoef = mod_tools.gen_coeff_signal1d(f, PSradio, parameters.ncomp2d)
        self.A_radio, self.phi_radio, self.fr_radio = gencoef
        return None

    def load_coeff(self, p):
        '''Load existing random realisations that has been stored in
        nadir+file_coeff. The outputs are the amplitude,
        the phase and the frequency of each random realisation.
        There are ncomp random realisations.'''
        try:
            ifile = '{}_nadir.nc'.format(parameters.file_coeff[:-3])
            fid = netCDF4.Dataset(ifile, 'r')
        except IOError:
            logger.error('There was an error opening the file nadir '
                         '{}'.format(ifile))
            sys.exit(1)
        self.A_wt = numpy.array(fid.variables['A_wt'][:]).squeeze()
        if numpy.shape(self.A_wt)[0] != self.ncomp2d:
            logger.error('{1} dimensions are different from ncomp2d={2}\n'
                         'remove {1} or adjust ncomp2d number in parameter'
                         'file'.format(p.file_coeff, self.ncomp2d))
            sys.exit(1)
        self.phi_wt = numpy.array(fid.variables['phi_wt'][:]).squeeze()
        self.frx_wt = numpy.array(fid.variables['frx_wt'][:]).squeeze()
        self.fry_wt = numpy.array(fid.variables['fry_wt'][:]).squeeze()
        self.A_radio = numpy.array(fid.variables['A_radio'][:]).squeeze()
        _tmpradio = fid.variables['phi_radio'][:]
        self.phi_radio = numpy.array(_tmpradio).squeeze()
        _tmpradio = fid.variables['fr_radio'][:]
        self.fr_radio = numpy.array(_tmpradio).squeeze()
        fid.close()
        return None

    def make_error(self, orb, cycle, SSH_true, p):
        ''' Build errors corresponding to each selected noise
        among the effect of the wet_tropo, and the instrumental error
        '''
        nal = numpy.shape(ssh_true)[0]
        # - compute random noise of 10**2 cm**2/(km/cycle)
        # - compute the correspond error on the nadir in m
        # - Initialization of radiometer error in right and left beam
        err_radio = numpy.zeros((nal))
        # - Initialization of swath matrices and large swath matrices
        #   (which include wet tropo data around the nadir and outside the
        #   swath)
        #   x_ac_large and wt_large are necessary to compute the gaussian
        #   footprint of a beam on the nadir or near the edge of the swath
        x_ac_large = numpy.arange(-2. * p.sigma/float(p.delta_ac),
                                  2.*p.sigma/float(p.delta_ac)+p.delta_ac,
                                  p.delta_ac)
        wt_large = numpy.zeros((numpy.shape(orb.x_al[:])[0],
                               numpy.shape(x_ac_large)[0]))
        y_al = numpy.float64(orb.x_al[:]) + float(cycle * orb.al_cycle)
        x_large, y_large = numpy.meshgrid(x_ac_large, y_al)
        # - Compute path delay error due to wet tropo and radiometer error
        #   using random coefficient initialized with power spectrums
        for comp in range(0, self.ncomp2d):
            phase_x_al_large = (2. * pi * (float(self.frx_wt[comp])
                                * (numpy.float64(x_large))
                                + float(self.fry_wt[comp])
                                * numpy.float64(y_large))) % (2.*pi)
            wt_large = (wt_large + self.A_wt[comp]
                        * numpy.cos(phase_x_al_large
                        + self.phi_wt[comp])*10**-2)
            phase_x_al = (2. * pi * float(self.fr_radio[comp])
                          * (numpy.float64(orb.x_al[:])
                          + float(cycle*orb.al_cycle))) % (2.*pi)
            err_radio = (err_radio + 2*self.A_radio[comp]
                         * numpy.cos(phase_x_al + self.phi_radio[comp])
                         * 10**-2)
        # - Compute Residual path delay error after a 1-beam radiometer
        #   correction
        beam = numpy.zeros((nal))
        diff_h1 = numpy.zeros((nal))
        # - Find across track indices in the gaussian footprint of
        #   2.*p.sigma
        indac = numpy.where((x_ac_large < 2.*p.sigma)
                            & (x_ac_large > -2.*p.sigma))[0]
        for i in range(0, nal):
            # - Find along track indices in the gaussian footprint of
            #   2.*p.sigma
            delta_al = orb.x_al[:] - orb.x_al[i]
            indal = numpy.where((delta_al <= (2*p.sigma))
                                & (delta_al > -2*p.sigma))[0]
            slice_ac = slice(min(indac), max(indac) + 1)
            slice_al = slice(min(indal), max(indal) + 1)
            x, y = numpy.meshgrid(x_ac_large[slice_ac],
                                  orb.x_al[slice_al] - orb.x_al[i])
            # - Compute path delay on gaussian footprint
            G = 1. / (2.*pi*p.sigma**2) * numpy.exp(-(x**2.+y**2.)
                                                    / (2.*p.sigma**2))
            beam[i] = (sum(sum(G*wt_large[slice_al, slice_ac]))
                       / sum(sum(G)) + err_radio[i])
        # - Filtering beam signal to cut frequencies higher than 125 km
        beam = gaussian_filter(beam, 30./p.delta_al)
        # - Compute residual path delay
        diff_h1 = wt_large[:, int(numpy.shape(wt_large)[1]/2.)] - beam
        self.wet_tropo1 = diff_h1
        self.wt = wt_large[:, int(numpy.shape(wt_large)[1]/2.)]
        return None
"""
