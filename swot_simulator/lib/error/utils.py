import numpy as np
import sys
import logging
import xarray as xr
# Define logging level for debug purposes
logger = logging.getLogger(__name__)


def read_file_instr(file_instr: str, delta_al: float, lambda_max:float): # -> xr.core.dataset.Dataset:
    """ Retrieve power spectrum from intrumental noise file provided by
    """
    #TODO proof errors if issue reading the file
    try:
        ds = xr.open_dataset(file_instr)
    except IOError:
        logger.error('There was an error opening the file'
                     '{}'.format(file_instr))
        sys.exit(1)
    # Set spatial frequency to spatial coordinate  
    ds.coords['nfreq'] = ds.spatial_frequency
    # Extract spatial frequency relevant to our sampling
    # (Cut frequencies larger than Nyquist frequency and long wavelength)
    cut_min = 1. / (2 * delta_al)
    cut_max = max(0, 1. / lambda_max)
    # TODO remove ugly trick from slice
    da = ds.sel(nfreq=slice(cut_max + 0.000001, cut_min - 0.000001))
    return da


def read_swh(swh_file):
    i = 0
    #TODO proof errors if issue reading the file
    try:
        ds = xr.open_dataset(swh_file)
    except IOError:
        logger.error('There was an error opening the file'
                     '{}'.format(swh_file))
        sys.exit(1)
    lon = ds.longitude[:]
    lat = ds.latitude[:]
    swh = ds.hs[i, :, :]
    swh = np.ma.masked_invalid(swh)
    swh[swh.mask] = numpy.nan

    fid.close()
    return lon, lat, swh



def read_file_karin(file_karin: str, swh): # -> xr.core.dataset.Dataset:
    """ Retrieve power spectrum from intrumental noise file provided by
    """
    #TODO proof errors if issue reading the file
    try:
        ds = xr.open_dataset(file_karin)
    except IOError:
        logger.error('There was an error opening the file'
                     '{}'.format(file_karin))
        sys.exit(1)
    file_hsdt = ds['height_sdt'].data
    x_ac = ds['cross_track'].data
    swh_tmp = ds['SWH'].data

    if isinstance(swh, (list, np.ndarray)):
        hsdt = np.full((len(swh), len(x_ac)), numpy.nan)
        for ial in range(len(swh)):
            # TODO take into account across track variability
            if not np.isfinite(swh[ial]):
                swh[ial] = 1
            i = np.argmin(abs(swh_tmp - swh[ial]))
            if swh_tmp[i] > swh[ial]:
                i += 1
            if np.max(swh_tmp) <= swh[ial]:
                _hsdt = hsdt[-1, :]
            else:
                rswh = swh[ial] - swh_tmp[i]
                _hsdt = file_hsdt[i, :] * (1 - rswh) + rswh * file_hsdt[i+1, :]
            hsdt[ial, :] = _hsdt
    else:
        i = np.argmin(abs(swh_tmp - swh))
        if swh_tmp[i] > swh:
            i += 1
        if np.max(swh_tmp) <= swh:
            hsdt = hsdt[-1, :]
            logger.warn('WARNING: swh={} is greater than the maximum value'
                        ' in {}, therefore swh is set to the file maximum'
                        'value '.format(swh, file_karin, np.max(swh_tmp)))
        else:
            rswh = swh - swh_tmp[i]
            hsdt = file_hsdt[i, :] * (1 - rswh) + rswh * file_hsdt[i + 1, :]
    return x_ac, hsdt


def gen_rcoeff_signal1d(f:np.ndarray, PS:np.ndarray,
                        lambda_min:float, lambda_max:float, npseudoper:int,
                        repeat:int, nseed:int): # ->Tuple[np.ndarray, np.ndarray]:
    '''Generate nc random coefficient from a spectrum PS
    with frequencies f. \n
    Return Amplitude, phase and frequency of nc realisations'''

    logffl = np.arange(np.log10(1. / lambda_max),
                       np.log10(1. / lambda_min + 1. / lambda_max),
                       np.log10(1 + 1. / npseudoper))
    ffl = 10**(logffl)
    logf = np.log10(f)
    logPS = np.log10(PS)
    logPSl = np.interp(logffl, logf, logPS)
    PSl = 10**(logPSl)
    A = np.sqrt(2 * PSl * (ffl / npseudoper))
    phi = []
    for k in range(len(ffl)):
        np.random.seed(nseed + k*1000)
        phi.append(2 * np.pi * np.random.random(int(2 * repeat * ffl[k]
                                                  / npseudoper + 3)))
    return A, phi


def gen_coeff_signal1d(f: np.ndarray, PS: np.ndarray, nc: int, nseed:int): # -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''Generate nc random coefficient from a spectrum PS
    with frequencies f. \n
    Return Amplitude, phase and frequency of nc realisations'''

    f1 = np.min(f)
    f2 = np.max(f)
    logf = np.log10(f)
    logf1 = np.log10(f1)
    logf2 = np.log10(f2)
    np.random.seed(nseed)
    # ''' Compute nc random vectors in [logf1, logf2] '''
    logfr = (logf2 - logf1)*np.random.random(nc) + logf1
    fr = 10**(logfr)
    # ''' Compute amplitude for random vectors '''
    logPS = np.log10(PS)
    logPSr = np.interp(logfr, logf, logPS)
    PSr = 10**(logPSr)
    A = np.sqrt(0.5 * PSr * fr * ((f2/f1)**(1./nc)-1))
    # ''' Compute nc phases in (0,2*pi)'''
    np.random.seed(nseed + 1000)
    phi = 2 * np.pi * np.random.random(nc)
    return A, phi, fr


def gen_signal1d(xx: np.ndarray, A:np.ndarray, phi:np.ndarray,
                 lmin:float, lmax:float, npseudoper:float)-> np.ndarray:
    '''Generate 1d random noise signal from coefficent computed using
     gen_rcoeff_signal1d. \n
    Return The random noise signal'''
    S = np.full(np.shape(xx), 0.)
    logffl = np.arange(np.log10(1. / lmax),
                       np.log10(1. / lmin + 1. / lmax),
                       np.log10(1 + 1. / npseudoper))
    ffl = 10**(logffl)
    for k in range(len(ffl)):
        ka = 2 * (xx * ffl[k] / npseudoper).astype(int) + 1
        Cka = np.abs(np.sin(2 * np.pi * xx * ffl[k] / npseudoper / 2.))
        kb = (2 * ((xx + npseudoper / 2. / ffl[k]) * ffl[k]
              / npseudoper).astype(int))
        Ckb = np.abs(np.sin(2 * np.pi * (xx + npseudoper / 2 / ffl[k])
                     * ffl[k] / npseudoper / 2))
        S = (S + A[k] * np.cos(2 * np.pi * ffl[k] * xx + phi[k][ka]) * Cka
             + A[k] * np.cos(2 * np.pi * ffl[k] * xx + phi[k][kb]) * Ckb)
    return S


def gen_coeff_signal2d(f: np.ndarray, PS: np.ndarray, nc:int, nseed:int): # -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''Generate nc random coefficient from a spectrum PS
    with frequencies f. \n
    Inputs are: frequency [f], spectrum [PS], number of realisation [nc]
    Return Amplitude, phase and frequency in 2D (frx, fry) of nc
    realisations'''

    # ''' Compute nc random vectors in an annular
    # (radius is between (min(f), max(f)) '''
    f1 = np.min(f)
    f2 = np.max(f)
    logf = np.log10(f)
    np.random.seed(nseed)
    fr = (f2 - f1) * np.random.random(nc) + f1
    # logfr = numpy.log10(fr)
    direction = 2. * np.pi * np.random.random(nc)
    frx = fr * np.cos(direction)
    fry = fr * np.sin(direction)

    # ''' Compute coefficients corresponding to random vectors '''
    logPS = np.log10(PS)
    logPSr = np.interp(np.log10(fr), logf, logPS)
    PSr = 10.**logPSr
    A = np.sqrt(0.5 * PSr * 2 * np.pi * (f2 - f1) / (nc))

    # ''' Compute nc phases in (0,2*pi)'''
    np.random.seed(nseed + 1000)
    phi = 2 * np.pi * np.random.random(nc)
    return A, phi, frx, fry

