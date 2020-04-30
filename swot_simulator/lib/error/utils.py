import numpy as np
from scipy import interpolate, fftpack
import sys
import logging
import xarray as xr
from typing import IO, Optional
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


def gen_signal1d(fi: np.ndarray, PSi: np.ndarray, x:np.ndarray,
                 nseed: Optional[int]=0, fmin: Optional[float]=None,
                 fmax: Optional[float]=None, alpha: Optional[int]=10,
                 lf_extpl: Optional[bool]=False,
                 hf_extpl: Optional[bool]=False) -> np.ndarray:
    ''' Generate 1d random signal using Fouriner coefficient '''
    # Make sure fi, PSi does not contain the zero frequency:
    PSi = PSi[fi>0]
    fi = fi[fi>0]

    interpolator = interpolate.interp1d
    # Interpolation function for the non-zero part of the spectrum
    finterp = interpolator(np.log(fi[PSi>0]), np.log(PSi[PSi>0]),
                           bounds_error=False, fill_value="extrapolate")

    # Adjust fmin and fmax to fi bounds if not specified:
    if fmin == None:
         fmin = fi[0]
    if fmax == None:
         fmax = fi[-1]
    # Go alpha times further in frequency to avoid interpolation aliasing.
    fmaxr = alpha * fmax

    f = np.arange(fmin, fmaxr + fmin, fmin)
    PS = np.exp(finterp(np.log(f)))

    # lf_extpl=True prolongates the PSi as a plateau below min(fi).
    # Otherwise, we consider zeros values. same for hf
    if lf_extpl is True:
         PS[f<fi[0]] = PSi[0]
    else:
         PS[f<fi[0]] = 0.
    if hf_extpl is True:
         PS[f>fi[-1]] = PSi[-1]
    else:
         PS[f>fi[-1]] = 0.
    PS[f>fmax]=0.

    # Detect the sections (if any) where PSi==0 and apply it to PS
    finterp_mask = interpolator(fi, PSi, bounds_error=False,
                                              fill_value="extrapolate")
    PSmask = finterp_mask(f)
    PS[PSmask==0.] = 0.

    phase = np.empty((2 * len(f) + 1))
    np.random.seed(nseed)
    phase[1:(len(f) + 1)] = np.random.random(len(f)) * 2 * np.pi
    phase[0] = 0.
    phase[-len(f):] = -phase[1:(len(f) + 1)][::-1]

    FFT1A = np.concatenate(([0], 0.5 * PS, 0.5 * PS[::-1]), axis=0)
    FFT1A = np.sqrt(FFT1A) * np.exp(1j * phase) / fmin**0.5

    yg = 2 * fmaxr * np.real(fftpack.ifft(FFT1A))
    xg = np.linspace(0, 0.5 / fmaxr * yg.shape[0], yg.shape[0])

    finterp = interpolator(xg, yg)
    y = finterp(np.mod(x, xg.max()))

    return y


def gen_signal2d_rectangle(fi: np.ndarray, PSi: np.ndarray, x: np.ndarray,
                           y:np.ndarray, fminx: Optional[float]=None,
                           fminy: Optional[float]=None,
                           fmax: Optional[float]=None,
                           alpha: Optional[int]=10, nseed: Optional[int]=0,
                           lf_extpl: Optional[bool]=False,
                           hf_extpl: Optional[bool]=False) -> np.ndarray:

    revert=False
    if fminy < fminx:
        revert = True
        fmin = +fminy
        fminy = +fminx
        x,y = y,x
    else:
        fmin = + fminx
    # Go alpha times further in frequency to avoid interpolation aliasing.
    fmaxr = alpha * fmax

    # Make sure fi, PSi does not contain the zero frequency:
    PSi = PSi[fi>0]
    fi = fi[fi>0]

    # Interpolation function for the non-zero part of the spectrum
    interp1 = interpolate.interp1d
    finterp = interp1(np.log(fi[PSi>0]), np.log(PSi[PSi>0]),
                      bounds_error=False, fill_value="extrapolate")


    f = np.arange(fmin, fmaxr + fmin, fmin)


    # lf_extpl=True prolongates the PSi as a plateau below min(fi).
    # Otherwise, we consider zeros values. same for hf
    if lf_extpl is True:
         PS[f<fi[0]] = PSi[0]
    else:
         PS[f<fi[0]] = 0.
    if hf_extpl is True:
         PS[f>fi[-1]] = PSi[-1]
    else:
         PS[f>fi[-1]] = 0.
    PS[f>fmax]=0.

    # Detect the sections (if any) where PSi==0 and apply it to PS
    finterp_mask = interp1(fi, PSi, bounds_error=False,
                           fill_value="extrapolate")
    PSmask = finterp_mask(f)
    PS[PSmask==0.] = 0.
    PS1D = PS

    # Build the 2D PSD following the given 1D PSD
    fx = np.concatenate(([0], f))
    fy = np.concatenate(([0], np.arange(fminy, fmaxr + fminy, fminy)))
    fx2, fy2 = np.meshgrid(fx, fy)
    f2 = np.sqrt((fx2**2 + fy2**2))
    dfx = fmin
    dfy = fminy

    PS2D = np.zeros(np.shape(f2))

    for iff in range(len(f)):
        fc = f2[:, (-iff-1)]
        ind1 = np.where((f2 >= (f[-iff-1]-dfx/2)) & (f2<(f[-iff-1]+dfx/2)))
        S = np.sum(PS2D[:,-iff-1]) * dfx * dfy
        MISS = PS1D[-iff-1] * dfx - S
        if MISS<=0:
            PS2D[ind1] = 0.
        else:
            PS2D[ind1] = MISS / len(ind1) / (dfx * dfy)

    PS2D[f2>fmax] = 0
    np.random.seed(nseed)
    phase = np.random.random((2*len(fy) - 1, len(fx))) * 2 * np.pi
    phase[0, 0] = 0.
    phase[-len(fy)+1:, 0] = -phase[1: len(fy), 0][::-1]

    FFT2A = np.concatenate((0.25 * PS2D, 0.25 * PS2D[1:, :][::-1, :]), axis=0)
    FFT2A = np.sqrt(FFT2A) * np.exp(1j * phase) / np.sqrt((dfx*dfy))
    FFT2 = np.zeros((2*len(fy)-1, 2*len(fx)-1), dtype=complex)
    FFT2[:, :len(fx)] = FFT2A
    FFT2[1:, -len(fx)+1:] = FFT2A[1:,1:].conj()[::-1, ::-1]
    FFT2[0, -len(fx)+1:] = FFT2A[0,1:].conj()[::-1]

    sg = (4 * fy[-1] * fx[-1]) * np.real(fftpack.ifft2(FFT2))
    xg = np.linspace(0, 1./fmin, sg.shape[1])
    yg = np.linspace(0, 1./fminy, sg.shape[0])
    xgmax = xg.max()
    ygmax = yg.max()
    finterp = interpolate.interp2d(xg, yg, sg)

    yl = y - y[0]
    yl = yl[yl<yg.max()]
    xl = x - x[0]
    xl = xl[xl<xg.max()]
    rectangle = finterp(xl,yl)
    x_n, x_r = np.divmod(x.max() - x[0], xgmax)
    y_n, y_r = np.divmod(y.max() - y[0], ygmax)

    signal = np.zeros((len(y), len(x)))

    for i_x_n in range(int(x_n + 1)):
        ix0 = np.where(((x - x[0]) >= (i_x_n * xgmax))
                       & ((x-x[0]) < ((i_x_n+1)*xgmax)))[0]
        for i_y_n in range(int(y_n + 1)):
            iy0 = np.where(((y - y[0]) >= (i_y_n * ygmax))
                           & ((y - y[0]) < ((i_y_n+1)*ygmax)))[0]
            _rect = rectangle[: len(iy0), : len(ix0)]
            signal[iy0[0]: iy0[-1] + 1, ix0[0]: ix0[-1] + 1] = _rect

    if revert is True:
        return signal.transpose()
    else:
        return signal
