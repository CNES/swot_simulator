import os
import pickle
import numpy as np
import swot_simulator.error.altimeter
import swot_simulator.error.baseline_dilation
import swot_simulator.error.karin
import swot_simulator.error.roll_phase
import swot_simulator.error.timing
import swot_simulator.error.wet_troposphere
import swot_simulator.random_signal
import swot_simulator.settings

ROOT = os.path.dirname(os.path.abspath(__file__))


def load_parameters():
    return swot_simulator.settings.Parameters(
        dict(ephemeris=os.path.join(ROOT, "..", "data",
                                    "ephemeris_calval_june2015_ell.txt"),
             nadir=True,
             swath=True,
             error_spectrum=os.path.join(ROOT, "..", "data",
                                         "global_sim_instrument_error.nc"),
             karin_noise=os.path.join(ROOT, "..", "data", "karin_noise_v2.nc"),
             cycle_duration=0.99349,
             height=857244))


def load_data(error=None):
    with open(os.path.join(ROOT, "data", "errors.bin"), "rb") as stream:
        (x_al, x_ac, curvilinear_distance, cycle_number,
         errors) = pickle.load(stream)
    return (x_al, x_ac, curvilinear_distance, cycle_number,
            errors[error] if error is not None else errors)


def test_error_karin():
    parameters = load_parameters()
    (x_al, x_ac, _curvilinear_distance, _cycle_number,
     _signal) = load_data('karin')
    error = swot_simulator.error.karin.Karin(parameters)
    swh = np.array([
        0,
    ] * len(x_ac))

    generated1 = error.generate(0, x_al, x_ac, swh)
    generated2 = error.generate(0, x_al, x_ac, swh)
    assert np.all(generated1['simulated_error_karin'] ==
                  generated2['simulated_error_karin'])
    generated2 = error.generate(1, x_al, x_ac, swh)
    assert np.all(generated1['simulated_error_karin'] !=
                  generated2['simulated_error_karin'])


def test_baseline_dilation():
    parameters = load_parameters()
    error_spectrum = swot_simulator.random_signal.read_file_instr(
        os.path.join(ROOT, "..", "data", "error_spectrum.nc"), 2.0)

    (x_al, x_ac, _, _, _) = load_data('baseline_dilation')
    error = swot_simulator.error.baseline_dilation.BaselineDilation(
        parameters, error_spectrum['dilationPSD'].data,
        error_spectrum['spatial_frequency'].data)
    generated = error.generate(x_al, x_ac)
    assert generated['simulated_error_baseline_dilation'].mean() < 1


def test_roll_phase():
    parameters = load_parameters()
    (x_al, x_ac, _, _, _) = load_data()
    error_spectrum = swot_simulator.random_signal.read_file_instr(
        os.path.join(ROOT, "..", "data", "error_spectrum.nc"), 2.0)

    error = swot_simulator.error.roll_phase.RollPhase(
        parameters, error_spectrum['rollPSD'].data,
        error_spectrum['gyroPSD'].data, error_spectrum['phasePSD'].data,
        error_spectrum['spatial_frequency'].data)
    generated = error.generate(x_al, x_ac)  # type: ignore
    assert generated['simulated_error_phase'].mean() < 1
    assert generated['simulated_error_roll'].mean() < 1


def test_timing():
    parameters = load_parameters()
    error_spectrum = swot_simulator.random_signal.read_file_instr(
        os.path.join(ROOT, "..", "data", "error_spectrum.nc"), 2.0)

    (x_al, x_ac, _, _, _) = load_data('timing')
    error = swot_simulator.error.timing.Timing(
        parameters, error_spectrum['timingPSD'].data,
        error_spectrum['spatial_frequency'].data)
    generated = error.generate(x_al, x_ac)
    assert generated['simulated_error_timing'].mean() < 1


def test_wet_troposphere():
    parameters = load_parameters()
    (x_al, x_ac, _, _, _) = load_data()

    parameters.nbeam = 1
    error = swot_simulator.error.wet_troposphere.WetTroposphere(parameters)
    generated = error.generate(x_al, x_ac)

    # assert abs(wt.data - expected['wt']['one']['wt']).mean() < 1e-12
    assert generated['simulated_error_troposphere'].mean() < 1

    # assert abs(wt_nadir.data -
    #            expected['wt_nadir']['one']['wt']).mean() < 1e-12
    assert generated['simulated_error_troposphere_nadir'].mean() < 1

    parameters.nbeam = 2
    error = swot_simulator.error.wet_troposphere.WetTroposphere(parameters)
    generated = error.generate(x_al, x_ac)

    # assert abs(wt.data - expected['wt']['two']['wt']).mean() < 1e-12
    assert generated['simulated_error_troposphere'].mean() < 1

    # assert abs(wt_nadir.data -
    #            expected['wt_nadir']['two']['wt']).mean() < 1e-12
    assert generated['simulated_error_troposphere_nadir'].mean() < 1


def test_altimeter():
    parameters = load_parameters()
    (x_al, _, _, _, _) = load_data('altimeter')

    error = swot_simulator.error.altimeter.Altimeter(parameters)
    generated = error.generate(x_al)
    assert generated['simulated_error_altimeter'].mean() < 1
