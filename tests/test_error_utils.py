import os
import pickle
import pytest
import xarray as xr

import swot_simulator.error.utils as utils

ROOT = os.path.dirname(os.path.abspath(__file__))


def test_gen_signal_1d():
    with open(os.path.join(ROOT, "data", "gen_signal_1d.bin"), "rb") as stream:
        (fi, psi, x, nseed, fmin, fmax, alpha, lf_extpl, hf_extpl,
         expected) = pickle.load(stream)

    result = utils.gen_signal_1d(fi, psi, x, nseed, fmin, fmax, alpha,
                                 lf_extpl, hf_extpl)
    assert (result - expected).mean() < 1e-12


def test_gen_signal_2d_rectangle():
    with open(os.path.join(ROOT, "data", "gen_signal_2d_rectangle.bin"),
              "rb") as stream:
        (fi, psi, x, y, fminx, fminy, fmax, alpha, nseed, lf_extpl, hf_extpl,
         expected) = pickle.load(stream)

    ps2d, f = utils.gen_ps2d(fi, psi, fminx, fminy, fmax, alpha, lf_extpl,
                             hf_extpl)
    result = utils.gen_signal_2d_rectangle(ps2d, f, x, y, fminx, fminy, fmax,
                                           alpha, nseed)
    assert (result - expected).mean() < 1e-12


def test_read_file_karin():
    cross_track, swh = utils.read_file_karin(
        os.path.join(ROOT, "..", "data", "karin_noise_v2.nc"), 2)
    assert cross_track.shape == swh.shape

    with pytest.warns(RuntimeWarning):
        cross_track, swh = utils.read_file_karin(
            os.path.join(ROOT, "..", "data", "karin_noise_v2.nc"), 200)


def test_read_file_instr():
    dataset = utils.read_file_instr(
        os.path.join(ROOT, "..", "data", "error_spectrum.nc"), 2.0, 20000)
    assert isinstance(dataset, xr.Dataset)
