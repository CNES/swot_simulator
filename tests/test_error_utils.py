import os
import pickle

import swot_simulator.error.utils as utils

ROOT = os.path.dirname(os.path.abspath(__file__))


def test_gen_signal_1d():
    with open(os.path.join(ROOT, "gen_signal_1d.bin"), "rb") as stream:
        (fi, psi, x, nseed, fmin, fmax, alpha, lf_extpl, hf_extpl,
         expected) = pickle.load(stream)

    result = utils.gen_signal_1d(fi, psi, x, nseed, fmin, fmax, alpha,
                                 lf_extpl, hf_extpl)
    assert (result - expected).mean() < 1e-16


def test_gen_signal_2d_rectangle():
    with open(os.path.join(ROOT, "gen_signal_2d_rectangle.bin"),
              "rb") as stream:
        (fi, psi, x, y, fminx, fminy, fmax, alpha, nseed, lf_extpl, hf_extpl,
         expected) = pickle.load(stream)

    result = utils.gen_signal_2d_rectangle(fi, psi, x, y, fminx, fminy, fmax,
                                           alpha, nseed, lf_extpl, hf_extpl)
    assert (result - expected).mean() < 1e-16
