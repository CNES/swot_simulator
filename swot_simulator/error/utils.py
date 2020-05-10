from typing import Optional, Tuple
import warnings
import numba as nb
import numpy as np
import scipy.interpolate
import xarray as xr
try:
    import mkl_fft
    IFFT = mkl_fft.ifft
    IFFT2 = mkl_fft.ifft2
except ImportError:
    IFFT = np.ifft
    IFFT2 = np.ifft2


def read_file_instr(file_instr: str, delta_al: float,
                    lambda_max: float) -> xr.Dataset:
    """ Retrieve power spectrum from intrumental noise file provided by
    """
    dataset = xr.load_dataset(file_instr)

    # Set spatial frequency to spatial coordinate
    dataset = dataset.swap_dims(dict(nfreq="spatial_frequency"))

    # Extract spatial frequency relevant to our sampling
    # (Cut frequencies larger than Nyquist frequency and long wavelength)
    cut_min = 1 / (2 * delta_al)
    cut_max = max(0, 1 / lambda_max)

    # TODO remove ugly trick from slice
    return dataset.sel(spatial_frequency=slice(cut_max + 0.000001, cut_min -
                                               0.000001))


# def read_swh(swh_file):
#     ds = xr.open_dataset(swh_file)

#     lon = ds.longitude[:]
#     lat = ds.latitude[:]
#     swh = ds.hs[0, :, :]
#     swh = np.ma.masked_invalid(swh)
#     swh[swh.mask] = numpy.nan

#     fid.close()
#     return lon, lat, swh


def read_file_karin(path: str,
                    threshold: float) -> Tuple[np.ndarray, np.ndarray]:
    """Retrieve power spectrum from intrumental noise file provided by
    """
    with xr.open_dataset(path) as dataset:
        height_sdt = dataset['height_sdt'].data
        cross_track = dataset['cross_track'].data
        swh = dataset['SWH'].data

    indices = np.argmin(np.abs(swh - threshold))
    if swh[indices] > threshold:
        indices += 1
    if swh.max() <= threshold:
        hsdt = np.full((len(cross_track), ), np.nan)
        warnings.warn(
            f'swh={threshold} is greater than the maximum value '
            f'in {path!r}, therefore swh is set to the file maximum '
            'value', RuntimeWarning)
    else:
        rswh = threshold - swh[indices]
        hsdt = height_sdt[indices, :] * (
            1 - rswh) + rswh * height_sdt[indices + 1, :]
    return cross_track, hsdt


def gen_signal_1d(fi: np.ndarray,
                  psi: np.ndarray,
                  x: np.ndarray,
                  nseed: int = 0,
                  fmin: Optional[float] = None,
                  fmax: Optional[float] = None,
                  alpha: int = 10,
                  lf_extpl: bool = False,
                  hf_extpl: bool = False) -> np.ndarray:
    """Generate 1d random signal using Fouriner coefficient"""
    # Make sure fi, PSi does not contain the zero frequency:
    psi = psi[fi > 0]
    fi = fi[fi > 0]

    # Adjust fmin and fmax to fi bounds if not specified
    fmin = fmin or fi[0]
    fmax = fmax or fi[-1]

    # Go alpha times further in frequency to avoid interpolation aliasing.
    fmaxr = alpha * fmax

    # Interpolation of the non-zero part of the spectrum
    f = np.arange(fmin, fmaxr + fmin, fmin)
    mask = psi > 0
    ps = np.exp(np.interp(np.log(f), np.log(fi[mask]), np.log(psi[mask])))

    # lf_extpl=True prolongates the PSi as a plateau below min(fi).
    # Otherwise, we consider zeros values. same for hf
    ps[f < fi[0]] = psi[0] if lf_extpl else 0
    ps[f > fi[-1]] = psi[-1] if hf_extpl else 0
    ps[f > fmax] = 0

    # Detect the sections (if any) where PSi==0 and apply it to PS
    mask = np.interp(f, fi, psi)
    ps[mask == 0] = 0

    f_size = f.size
    phase = np.empty((2 * f_size + 1))
    np.random.seed(nseed)
    phase[1:(f_size + 1)] = np.random.random(f_size) * 2 * np.pi
    phase[0] = 0
    phase[-f_size:] = -phase[1:(f_size + 1)][::-1]

    fft1a = np.concatenate((np.array([0]), 0.5 * ps, 0.5 * ps[::-1]), axis=0)
    fft1a = np.sqrt(fft1a) * np.exp(1j * phase) / fmin**0.5

    yg = 2 * fmaxr * np.real(IFFT(fft1a))
    xg = np.linspace(0, 0.5 / fmaxr * yg.shape[0], yg.shape[0])

    return np.interp(np.mod(x, xg.max()), xg, yg)


@nb.njit("(float64[:, ::1])"
         "(float64[::1], float64[:, ::1], float64[::1], float64, float64)",
         cache=True)
def _calculate_ps2d(f: np.ndarray, f2: np.ndarray, ps1d: np.ndarray,
                    dfx: np.ndarray, dfy: np.ndarray) -> np.ndarray:
    result = np.zeros(f2.shape)
    view = result.ravel()
    dfx_2 = dfx * 0.5
    dfx_y = dfx * dfy
    for idx in range(-1, -f.size - 1, -1):
        item = f[idx]
        mask = (f2 >= (item - dfx_2)) & (f2 < (item + dfx_2))
        amount = np.sum(result[:, idx]) * dfx_y
        miss = ps1d[idx] * dfx - amount
        view[mask.ravel()] = 0 if miss <= 0 else miss * 0.5 / dfx_y
    return result


@nb.njit("(float64[:, ::1])"
         "(float64[:, ::1], float64[::1], float64[::1], float64, float64)",
         cache=True)
def _calculate_signal(rectangle, x, y, xgmax, ygmax):
    x_n = (x.max() - x[0]) // xgmax
    y_n = (y.max() - y[0]) // ygmax
    result = np.zeros((len(y), len(x)))

    for i_x_n in range(int(x_n + 1)):
        ix0 = np.where(((x - x[0]) >= (i_x_n * xgmax))
                       & ((x - x[0]) < ((i_x_n + 1) * xgmax)))[0]
        for i_y_n in range(int(y_n + 1)):
            iy0 = np.where(((y - y[0]) >= (i_y_n * ygmax))
                           & ((y - y[0]) < ((i_y_n + 1) * ygmax)))[0]
            _rect = rectangle[:len(iy0), :len(ix0)]
            result[iy0[0]:iy0[-1] + 1, ix0[0]:ix0[-1] + 1] = _rect
    return result


def gen_signal_2d_rectangle(fi: np.ndarray,
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
    f = np.arange(fmin, fmaxr + fmin, fmin)
    ps = np.exp(np.interp(np.log(f), np.log(fi[psi > 0]),
                          np.log(psi[psi > 0])))

    # lf_extpl=True prolongates the PSi as a plateau below min(fi).
    # Otherwise, we consider zeros values. same for hf
    ps[f < fi[0]] = psi[0] if lf_extpl else 0
    ps[f > fi[-1]] = psi[-1] if hf_extpl else 0
    ps[f > fmax] = 0

    # Detect the sections (if any) where PSi==0 and apply it to PS
    psmask = np.interp(f, fi, psi)
    ps[psmask == 0] = 0
    ps1d = ps

    # Build the 2D PSD following the given 1D PSD
    fx = np.concatenate(([0], f))
    fy = np.concatenate(([0], np.arange(fminy, fmaxr + fminy, fminy)))
    fx2, fy2 = np.meshgrid(fx, fy)
    f2 = np.sqrt((fx2**2 + fy2**2))
    dfx = fmin
    dfy = fminy
    ps2d = _calculate_ps2d(f, f2, ps1d, dfx, dfy)
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

    sg = (4 * fy[-1] * fx[-1]) * np.real(IFFT2(fft2))
    xg = np.linspace(0, 1 / fmin, sg.shape[1])
    yg = np.linspace(0, 1 / fminy, sg.shape[0])
    xgmax = xg.max()
    ygmax = yg.max()

    yl = y - y[0]
    yl = yl[yl < yg.max()]
    xl = x - x[0]
    xl = xl[xl < xg.max()]
    rectangle = np.ascontiguousarray(
        scipy.interpolate.interp2d(xg, yg, sg)(xl, yl))
    signal = _calculate_signal(rectangle, x, y, xgmax, ygmax)

    return signal.transpose() if revert else signal
