import os
import pickle
import swot_simulator.settings
import swot_simulator.error.altimeter
import swot_simulator.error.baseline_dilation
import swot_simulator.error.karin
import swot_simulator.error.roll_phase
import swot_simulator.error.timing
import swot_simulator.error.utils
import swot_simulator.error.wet_troposphere

ROOT = os.path.dirname(os.path.abspath(__file__))


def load_parameters():
    return swot_simulator.settings.Parameters(
        dict(ephemeris=os.path.join(ROOT, "..", "data",
                                    "ephem_calval_june2015_ell.txt"),
             nadir=True,
             swath=True,
             error_spectrum=os.path.join(ROOT, "..", "data",
                                         "global_sim_instrument_error.nc"),
             karin_noise=os.path.join(ROOT, "..", "data",
                                      "karin_noise_v2.nc")))


def load_data(error=None):
    with open(os.path.join(ROOT, "data", "errors.bin"), "rb") as stream:
        (x_al, x_ac, curvilinear_distance, cycle_number,
         errors) = pickle.load(stream)
    return (x_al, x_ac, curvilinear_distance, cycle_number,
            errors[error] if error is not None else errors)


def test_error_karin():
    parameters = load_parameters()
    (x_al, x_ac, curvilinear_distance, cycle_number,
     expected) = load_data('karin')
    error = swot_simulator.error.karin.Karin(parameters)

    generated = error.generate(x_al, x_ac, curvilinear_distance, cycle_number)
    assert abs(expected - generated['karin']).mean() < 1e-12


def test_baseline_dilation():
    parameters = load_parameters()
    error_spectrum = swot_simulator.error.utils.read_file_instr(
        os.path.join(ROOT, "..", "data", "error_spectrum.nc"), 2.0, 20000)

    (x_al, x_ac, _, _, expected) = load_data('baseline_dilation')
    error = swot_simulator.error.baseline_dilation.BaselineDilation(
        parameters, error_spectrum['dilationPSD'].data,
        error_spectrum['spatial_frequency'].data)
    generated = error.generate(x_al, x_ac)
    assert abs(expected - generated['baseline_dilation']).mean() < 1e-12


def test_roll_phase():
    parameters = load_parameters()
    (x_al, x_ac, _, _, expected) = load_data()
    error_spectrum = swot_simulator.error.utils.read_file_instr(
        os.path.join(ROOT, "..", "data", "error_spectrum.nc"), 2.0, 20000)

    error = swot_simulator.error.roll_phase.RollPhase(
        parameters, error_spectrum['rollPSD'].data,
        error_spectrum['phasePSD'].data,
        error_spectrum['spatial_frequency'].data)
    generated = error.generate(x_al, x_ac)
    assert abs(expected['phase'] - generated['phase']).mean() < 1e-12
    assert abs(expected['roll'] - generated['roll']).mean() < 1e-12


def test_timing():
    parameters = load_parameters()
    error_spectrum = swot_simulator.error.utils.read_file_instr(
        os.path.join(ROOT, "..", "data", "error_spectrum.nc"), 2.0, 20000)

    (x_al, x_ac, _, _, expected) = load_data('timing')
    error = swot_simulator.error.timing.Timing(
        parameters, error_spectrum['timingPSD'].data,
        error_spectrum['spatial_frequency'].data)
    generated = error.generate(x_al, x_ac)
    assert abs(expected - generated['timing']).mean() < 1e-12


def test_wet_troposphere():
    parameters = load_parameters()
    (x_al, x_ac, _, _, expected) = load_data()

    parameters.nbeam = 1
    error = swot_simulator.error.wet_troposphere.WetTroposphere(parameters)
    generated = error.generate(x_al, x_ac)

    # assert abs(wt.data - expected['wt']['one']['wt']).mean() < 1e-12
    assert abs(generated['wet_troposphere'] -
               expected['wt']['one']['wet_tropo']).mean() < 1e-12

    # assert abs(wt_nadir.data -
    #            expected['wt_nadir']['one']['wt']).mean() < 1e-12
    assert abs(generated['wet_troposphere_nadir'] -
               expected['wt_nadir']['one']['wet_tropo']).mean() < 1e-12

    parameters.nbeam = 2
    error = swot_simulator.error.wet_troposphere.WetTroposphere(parameters)
    generated = error.generate(x_al, x_ac)

    # assert abs(wt.data - expected['wt']['two']['wt']).mean() < 1e-12
    assert abs(generated['wet_troposphere'] -
               expected['wt']['two']['wet_tropo']).mean() < 1e-12

    # assert abs(wt_nadir.data -
    #            expected['wt_nadir']['two']['wt']).mean() < 1e-12
    assert abs(generated['wet_troposphere_nadir'].data -
               expected['wt_nadir']['two']['wet_tropo']).mean() < 1e-12


def test_altimeter():
    parameters = load_parameters()
    (x_al, _, _, _, expected) = load_data('altimeter')

    error = swot_simulator.error.altimeter.Altimeter(parameters)
    generated = error.generate(x_al)
    assert abs(generated['altimeter'] - expected).mean() < 1e-12
