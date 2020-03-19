import numpy as np
import swotsimulator.rw_data as rw_data
import swotsimulator.const as const
import swotsimulator.mod_tools as mod_tools
from math import pi, sqrt
import netCDF4
import sys
from scipy.ndimage.filters import gaussian_filter
import logging

# Define logging level for debug purposes
logger = logging.getLogger(__name__)


def generate_randreal(parameters: settings.Parameters):
    '''Run reprodictible: generate or load nrand random numbers:
         - If one instrumental error needs to be computed, load frequencies
       and spectrum from the instrumental netcdf file:
    '''
    file_instr = read_file_instr(file=parameters.file_inst_error)
    spatial_frequency = []
    file_instr.read_var(spatial_frequency=spatial_frequency)
    # - Cut frequencies larger than Nyquist frequency and cut long
    #   wavelengths (larger than p.lambda_max)
    cut_min = 1. / float(2 * parameters.delta_al)
    cut_max = 1. / parameters.lambda_max
    ind = numpy.where((file_instr.spatial_frequency < cut_min)
                              & (file_instr.spatial_frequency > 0)
                              & (file_instr.spatial_frequency > cut_max))[0]
    freq = file_instr.spatial_frequency[ind]
    f_cut = 1. / parameters.lambda_cut
    ind_cut = numpy.where(freq < f_cut)
    min_ind_cut = numpy.min(numpy.where(freq >= f_cut))
    return file_instr, ind_cut, min_ind_cut, freq




class error_altimeter():
    '''Class errornadir defines the error on the nadir.
    Random realisation of errors can be initialized using init_error.
    If the random realisations have already been computed and stored in file
    file_coeff, the random realisations are read directly using load_coeff.
    The correspondg errors on a swath can be computed using make_error. '''
    def __init__(self, p, nadir=None, wet_tropo1=None, wt=None):
        self.nadir = nadir
        self.wet_tropo1 = wet_tropo1
        self.wt = wt
        self.nrand = getattr(p, 'nrandkarin', 1000)
        p.nrandkarin = self.nrand
        self.ncomp2d = getattr(p, 'ncomp2d', 2000)
        p.ncomp2d = self.ncomp2d
        self.ncomp1d = getattr(p, 'ncomp1d', 2000)
        p.ncomp1d = self.ncomp1d

    def init_error(self, parameters):
        '''Initialization of errors: Random realisation of errors are computed
        using a known power spectrum.
        The outputs are the amplitude, the phase and the frequency of each
        random realisation.
        By default, there are ncomp2d=2000 random realisations for the
        wet tropo and ncomp1d=2000 random realisations for the nadir 1d
        spectrum error.'''
        # Run reprodictible: generate or load nrand random numbers:
        # - Compute random coefficients in 1D for the nadir error
        wnoise = getattr(p, 'wnoise', 100)
        p.wnoise = wnoise
        dal = parameters.delta_al
        lmax = parameters.lambda_max
        npseudoper = parameters.npseudoper
        # - Define the sepctrum of the nadir instrument error
        # self.A=numpy.random.normal(0.0,sqrt(p.wnoise)
        # /numpy.float64(sqrt(2*p.delta_al)), (self.nrand))*0.01
        f = numpy.arange(1./3000., 1./float(2.*p.delta_al), 1./3000.)
        PSD = 8 + 1.05 * 10**(-4) * f**(-2.2)
        indf = numpy.where(f < 0.00023627939582672978)
        PSD[indf] = 10**4
        # Convert spectrum in m2/cy
        PSD = PSD * 10**(-4)
        if parameters.savesignal is True:
            gencoef = mod_tools.gen_rcoeff_signal1d(f, PSD, 2 * dal, lmax,
                                                    npseudoper,
                                                    parameters.len_repeat)
            self.A, self.phi = gencoef
        else:
            gencoef = mod_tools.gen_coeff_signal1d(f, PSD, parameters.ncomp1d)
            self.A, self.phi, self.f = gencoef
        return None

    def load_coeff(self, p):
        '''Load existing random realisations that has been stored in
        nadir+file_coeff. The outputs are the amplitude,
        the phase and the frequency of each random realisation.
        There are ncomp random realisations.'''
        try:
            fid = netCDF4.Dataset('{}_nadir.nc'.format(p.file_coeff[:-3]), 'r')
        except IOError:
            logger.error('There was an error opening the file nadir '
                         '{}_nadir.nc'.format(p.file_coeff[:-3]))
            sys.exit(1)
        self.A = numpy.array(fid.variables['A'][:]).squeeze()
        self.f = numpy.array(fid.variables['f'][:]).squeeze()
        self.phi = numpy.array(fid.variables['phi'][:]).squeeze()
        return None

    def make_error(self, orb, cycle, SSH_true, p):
        ''' Build errors corresponding to each selected noise
        among the effect of the wet_tropo, and the instrumental error
        '''
        nal = numpy.shape(SSH_true)[0]
        # - Compute random noise of 10**2 cm**2/(km/cycle)
        # - Compute the correspond error on the nadir in m
        if p.savesignal is True:
            xx = (numpy.float64(orb.x_al[:]) + float(cycle
                  * orb.al_cycle)) % p.len_repeat
            errnadir = mod_tools.gen_signal1d(xx, self.A, self.phi,
                                              2 * p.delta_al, p.lambda_max,
                                              p.npseudoper)
        else:
            errnadir = numpy.zeros((nal))
            for comp in range(0, self.ncomp1d):
                phase_x_al = (2. * pi * float(self.f[comp])
                              * (numpy.float64(orb.x_al[:])
                              + float(cycle*orb.al_cycle))) % (2.*pi)
                errnadir[:] = (errnadir[:] + 2*self.A[comp]
                               * numpy.cos(phase_x_al[:]+self.phi[comp]))
        # - Compute the correspond timing error on the swath in m
        self.nadir = errnadir[:]
        return None

    def save_coeff(self, p):
        '''Save random realisations to enable runs to be reproducible.
        The ncomp1d random phase phi, amplitude A and frequency fr for
        1D spectrum and ncomp2d random phase phi, amplitude A and frequencies
        frx and fry for 2D spectrum are saved in nadirfile_coeff for each error
        and can be loaded using load_coeff.
        '''
        # - Open Netcdf file in write mode
        fid = netCDF4.Dataset('{}_nadir.nc'.format(p.file_coeff[:-3]), 'w')
        fid.description = "Random coefficients from orbit simulator"

        # - Create dimensions
        fid.createDimension('nrand1d', self.ncomp1d)
        fid.createDimension('nrand2d', self.ncomp2d)
        # - Create and write Variables
        var = fid.createVariable('A', 'f4', ('nrand1d', ))
        var[:] = self.A
        var = fid.createVariable('f', 'f4', ('nrand1d', ))
        var[:] = self.f
        var = fid.createVariable('phi', 'f4', ('nrand1d', ))
        var[:] = self.phi

        # var = fid.createVariable('phi', 'f4', ('ninstr',))
        # var[:] = self.phi
        # var = fid.createVariable('fr', 'f4', ('ninstr',))
        # var[:] = self.fr
        fid.close()
        return None
