import os
import pickle
import numpy as np
import xarray as xr

import swot_simulator
import swot_simulator.random_signal as random_signal

ROOT = os.path.dirname(os.path.abspath(__file__))


def test_gen_signal_1d():
    with open(os.path.join(ROOT, "data", "gen_signal_1d.bin"), "rb") as stream:
        (fi, psi, x, nseed, fmin, fmax, alpha, lf_extpl, hf_extpl,
         _expected) = pickle.load(stream)

    rng = np.random.default_rng(seed=nseed)
    result = random_signal.gen_signal_1d(fi, psi, x, rng, fmin, fmax, alpha,
                                         lf_extpl, hf_extpl)
    assert result.mean() < 1


def test_gen_signal_2d_rectangle():
    with open(os.path.join(ROOT, "data", "gen_signal_2d_rectangle.bin"),
              "rb") as stream:
        (fi, psi, x, y, fminx, fminy, fmax, alpha, nseed, lf_extpl, hf_extpl,
         _expected) = pickle.load(stream)

    ps2d, f = random_signal.gen_ps2d(fi, psi, fminx, fminy, fmax, alpha,
                                     lf_extpl, hf_extpl)
    rng = np.random.default_rng(seed=nseed)
    result = random_signal.gen_signal_2d_rectangle(ps2d, f, x, y, rng, fminx,
                                                   fminy, fmax, alpha)
    assert result.mean() < 1


def test_read_file_karin():
    height_sdt, cross_track, swh = random_signal.read_file_karin(
        str(swot_simulator.DATA.joinpath("karin_noise_v2.nc")))


def test_read_file_instr():
    dataset = random_signal.read_file_instr(
        str(swot_simulator.DATA.joinpath("error_spectrum.nc")), 2.0)
    assert isinstance(dataset, xr.Dataset)
