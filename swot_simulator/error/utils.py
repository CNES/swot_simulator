from typing import Optional
import logging
import numpy as np
import scipy.interpolate
import scipy.fftpack
import xarray as xr

LOGGER = logging.getLogger(__name__)


def read_file_instr(file_instr: str, delta_al: float,
                    lambda_max: float) -> xr.Dataset:
    """ Retrieve power spectrum from intrumental noise file provided by
    """
    ds = xr.load_dataset(file_instr)

    # Set spatial frequency to spatial coordinate
    ds.coords['nfreq'] = ds.spatial_frequency

    # Extract spatial frequency relevant to our sampling
    # (Cut frequencies larger than Nyquist frequency and long wavelength)
    cut_min = 1 / (2 * delta_al)
    cut_max = max(0, 1 / lambda_max)

    # TODO remove ugly trick from slice
    da = ds.sel(nfreq=slice(cut_max + 0.000001, cut_min - 0.000001))
    return da


# def read_swh(swh_file):
#     ds = xr.open_dataset(swh_file)

#     lon = ds.longitude[:]
#     lat = ds.latitude[:]
#     swh = ds.hs[0, :, :]
#     swh = np.ma.masked_invalid(swh)
#     swh[swh.mask] = numpy.nan

#     fid.close()
#     return lon, lat, swh


def read_file_karin(file_karin: str, swh):  # -> xr.core.dataset.Dataset:
    """ Retrieve power spectrum from intrumental noise file provided by
    """
    ds = xr.open_dataset(file_karin)

    file_hsdt = ds['height_sdt'].data
    x_ac = ds['cross_track'].data
    swh_tmp = ds['SWH'].data

    if isinstance(swh, (list, np.ndarray)):
        hsdt = np.full((len(swh), len(x_ac)), np.nan)
        for idx, item in enumerate(swh):
            # TODO take into account across track variability
            if not np.isfinite(item):
                item = 1
            indices = np.argmin(abs(swh_tmp - item))
            if swh_tmp[indices] > item:
                indices += 1
            if np.max(swh_tmp) <= item:
                _hsdt = hsdt[-1, :]
            else:
                rswh = item - swh_tmp[indices]
                _hsdt = file_hsdt[indices, :] * (
                    1 - rswh) + rswh * file_hsdt[indices + 1, :]
            hsdt[idx, :] = _hsdt
    else:
        indices = np.argmin(abs(swh_tmp - swh))
        if swh_tmp[indices] > swh:
            indices += 1
        if np.max(swh_tmp) <= swh:
            hsdt = np.full((len(swh), len(x_ac)), np.nan)
            hsdt = hsdt[-1, :]
            LOGGER.warning(
                'WARNING: swh=%r is greater than the maximum value'
                ' in %s, therefore swh is set to the file maximum'
                'value ', swh, file_karin)
        else:
            rswh = swh - swh_tmp[indices]
            hsdt = file_hsdt[indices, :] * (
                1 - rswh) + rswh * file_hsdt[indices + 1, :]
    return x_ac, hsdt


def gen_signal1d(fi: np.ndarray,
                 psi: np.ndarray,
                 x: np.ndarray,
                 nseed: int = 0,
                 fmin: Optional[float] = None,
                 fmax: Optional[float] = None,
                 alpha: int = 10,
                 lf_extpl: Optional[bool] = False,
                 hf_extpl: Optional[bool] = False) -> np.ndarray:
    """Generate 1d random signal using Fouriner coefficient"""
    # Make sure fi, PSi does not contain the zero frequency:
    psi = psi[fi > 0]
    fi = fi[fi > 0]

    interpolator = scipy.interpolate.interp1d
    # Interpolation function for the non-zero part of the spectrum
    finterp = interpolator(np.log(fi[psi > 0]),
                           np.log(psi[psi > 0]),
                           bounds_error=False,
                           fill_value="extrapolate")

    # Adjust fmin and fmax to fi bounds if not specified:
    fmin = fmin or fi[0]
    fmax = fmax or fi[-1]

    # Go alpha times further in frequency to avoid interpolation aliasing.
    fmaxr = alpha * fmax

    f = np.arange(fmin, fmaxr + fmin, fmin)
    ps = np.exp(finterp(np.log(f)))

    # lf_extpl=True prolongates the PSi as a plateau below min(fi).
    # Otherwise, we consider zeros values. same for hf
    ps[f < fi[0]] = psi[0] if lf_extpl else 0
    ps[f > fi[-1]] = psi[-1] if hf_extpl else 0
    ps[f > fmax] = 0

    # Detect the sections (if any) where PSi==0 and apply it to PS
    finterp_mask = interpolator(fi,
                                psi,
                                bounds_error=False,
                                fill_value="extrapolate")
    psmask = finterp_mask(f)
    ps[psmask == 0.] = 0.

    phase = np.empty((2 * len(f) + 1))
    np.random.seed(nseed)
    phase[1:(len(f) + 1)] = np.random.random(len(f)) * 2 * np.pi
    phase[0] = 0.
    phase[-len(f):] = -phase[1:(len(f) + 1)][::-1]

    fft1a = np.concatenate(([0], 0.5 * ps, 0.5 * ps[::-1]), axis=0)
    fft1a = np.sqrt(fft1a) * np.exp(1j * phase) / fmin**0.5

    yg = 2 * fmaxr * np.real(scipy.fftpack.ifft(fft1a))
    xg = np.linspace(0, 0.5 / fmaxr * yg.shape[0], yg.shape[0])

    finterp = interpolator(xg, yg)
    y = finterp(np.mod(x, xg.max()))

    return y


def gen_signal2d_rectangle(fi: np.ndarray,
                           psi: np.ndarray,
                           x: np.ndarray,
                           y: np.ndarray,
                           fminx: Optional[float] = None,
                           fminy: Optional[float] = None,
                           fmax: Optional[float] = None,
                           alpha: int = 10,
                           nseed: int = 0,
                           lf_extpl: bool = False,
                           hf_extpl: bool = False) -> np.ndarray:

    revert = False
    if fminy < fminx:
        revert = True
        fmin = fminy
        fminy = fminx
        x, y = y, x
    else:
        fmin = fminx

    # Go alpha times further in frequency to avoid interpolation aliasing.
    fmaxr = alpha * fmax

    # Make sure fi, PSi does not contain the zero frequency:
    psi = psi[fi > 0]
    fi = fi[fi > 0]

    # Interpolation function for the non-zero part of the spectrum
    interp1 = scipy.interpolate.interp1d
    finterp = interp1(np.log(fi[psi > 0]),
                      np.log(psi[psi > 0]),
                      bounds_error=False,
                      fill_value="extrapolate")
    f = np.arange(fmin, fmaxr + fmin, fmin)
    ps = np.exp(finterp(np.log(f)))

    # lf_extpl=True prolongates the PSi as a plateau below min(fi).
    # Otherwise, we consider zeros values. same for hf
    ps[f < fi[0]] = psi[0] if lf_extpl else 0
    ps[f > fi[-1]] = psi[-1] if hf_extpl else 0
    ps[f > fmax] = 0

    # Detect the sections (if any) where PSi==0 and apply it to PS
    finterp_mask = interp1(fi,
                           psi,
                           bounds_error=False,
                           fill_value="extrapolate")
    psmask = finterp_mask(f)
    ps[psmask == 0] = 0
    ps1d = ps

    # Build the 2D PSD following the given 1D PSD
    fx = np.concatenate(([0], f))
    fy = np.concatenate(([0], np.arange(fminy, fmaxr + fminy, fminy)))
    fx2, fy2 = np.meshgrid(fx, fy)
    f2 = np.sqrt((fx2**2 + fy2**2))
    dfx = fmin
    dfy = fminy

    ps2d = np.zeros(np.shape(f2))

    for iff in range(len(f)):
        ind1 = np.where((f2 >= (f[-iff - 1] - dfx / 2))
                        & (f2 < (f[-iff - 1] + dfx / 2)))
        s = np.sum(ps2d[:, -iff - 1]) * dfx * dfy
        miss = ps1d[-iff - 1] * dfx - s
        if miss <= 0:
            ps2d[ind1] = 0.
        else:
            ps2d[ind1] = miss / len(ind1) / (dfx * dfy)

    ps2d[f2 > fmax] = 0
    np.random.seed(nseed)
    phase = np.random.random((2 * len(fy) - 1, len(fx))) * 2 * np.pi
    phase[0, 0] = 0.
    phase[-len(fy) + 1:, 0] = -phase[1:len(fy), 0][::-1]

    fft2a = np.concatenate((0.25 * ps2d, 0.25 * ps2d[1:, :][::-1, :]), axis=0)
    fft2a = np.sqrt(fft2a) * np.exp(1j * phase) / np.sqrt((dfx * dfy))
    fft2 = np.zeros((2 * len(fy) - 1, 2 * len(fx) - 1), dtype=complex)
    fft2[:, :len(fx)] = fft2a
    fft2[1:, -len(fx) + 1:] = fft2a[1:, 1:].conj()[::-1, ::-1]
    fft2[0, -len(fx) + 1:] = fft2a[0, 1:].conj()[::-1]

    sg = (4 * fy[-1] * fx[-1]) * np.real(scipy.fftpack.ifft2(fft2))
    xg = np.linspace(0, 1. / fmin, sg.shape[1])
    yg = np.linspace(0, 1. / fminy, sg.shape[0])
    xgmax = xg.max()
    ygmax = yg.max()
    finterp = scipy.interpolate.interp2d(xg, yg, sg)

    yl = y - y[0]
    yl = yl[yl < yg.max()]
    xl = x - x[0]
    xl = xl[xl < xg.max()]
    rectangle = finterp(xl, yl)
    x_n, _ = np.divmod(x.max() - x[0], xgmax)
    y_n, _ = np.divmod(y.max() - y[0], ygmax)

    signal = np.zeros((len(y), len(x)))

    for i_x_n in range(int(x_n + 1)):
        ix0 = np.where(((x - x[0]) >= (i_x_n * xgmax))
                       & ((x - x[0]) < ((i_x_n + 1) * xgmax)))[0]
        for i_y_n in range(int(y_n + 1)):
            iy0 = np.where(((y - y[0]) >= (i_y_n * ygmax))
                           & ((y - y[0]) < ((i_y_n + 1) * ygmax)))[0]
            _rect = rectangle[:len(iy0), :len(ix0)]
            signal[iy0[0]:iy0[-1] + 1, ix0[0]:ix0[-1] + 1] = _rect

    return signal.transpose() if revert else signal
