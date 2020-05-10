import os
import pickle
import swot_simulator.settings
import swot_simulator.error.altimeter
import swot_simulator.error.baseline_dilation
import swot_simulator.error.karin
import swot_simulator.error.roll_phase
import swot_simulator.error.timing
import swot_simulator.error.utils
import swot_simulator.error.wt

ROOT = os.path.dirname(os.path.abspath(__file__))


def load_parameters():
    return swot_simulator.settings.Parameters(
        dict(ephemeris=os.path.join(ROOT, "..", "data",
                                    "ephem_calval_june2015_ell.txt"),
             nadir=True,
             swath=True,
             error_spectrum_file=os.path.join(
                 ROOT, "..", "data", "global_sim_instrument_error.nc"),
             karin_file=os.path.join(ROOT, "..", "data", "karin_noise_v2.nc")))


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
    error = swot_simulator.error.karin.Karin(parameters, x_ac.shape[0])

    generated = error.generate(x_al, x_ac, curvilinear_distance, cycle_number)
    assert abs(expected - generated.data).mean() < 1e-20


def test_baseline_dilation():
    parameters = load_parameters()
    error_spectrum = swot_simulator.error.utils.read_file_instr(
        os.path.join(ROOT, "..", "data", "error_spectrum_file.nc"), 2.0, 20000)

    (x_al, x_ac, _, _, expected) = load_data('baseline_dilation')
    error = swot_simulator.error.baseline_dilation.BaselineDilation(
        parameters, error_spectrum['dilationPSD'].data,
        error_spectrum['spatial_frequency'].data)
    generated = error.generate(x_al, x_ac)
    assert abs(expected - generated.data).mean() < 1e-17


def test_roll_phase():
    parameters = load_parameters()
    (x_al, x_ac, _, _, expected) = load_data()
    error_spectrum = swot_simulator.error.utils.read_file_instr(
        os.path.join(ROOT, "..", "data", "error_spectrum_file.nc"), 2.0, 20000)

    error = swot_simulator.error.roll_phase.RollPhase(
        parameters, error_spectrum['rollPSD'].data,
        error_spectrum['phasePSD'].data,
        error_spectrum['spatial_frequency'].data)
    generated = error.generate(x_al, x_ac)
    assert abs(expected['phase'] - generated[0].data).mean() < 1e-16
    assert abs(expected['roll'] - generated[1].data).mean() < 1e-17


def test_timing():
    parameters = load_parameters()
    error_spectrum = swot_simulator.error.utils.read_file_instr(
        os.path.join(ROOT, "..", "data", "error_spectrum_file.nc"), 2.0, 20000)

    (x_al, x_ac, _, _, expected) = load_data('timing')
    error = swot_simulator.error.timing.Timing(
        parameters, error_spectrum['timingPSD'].data,
        error_spectrum['spatial_frequency'].data)
    generated = error.generate(x_al, x_ac)
    assert abs(expected - generated).mean() < 1e-17


def test_wet_troposphere():
    parameters = load_parameters()
    (x_al, x_ac, _, _, expected) = load_data()

    parameters.nbeam = 1
    error = swot_simulator.error.wt.WetTroposphere(parameters)
    wt, wet_tropo, wt_nadir, wet_tropo_nadir = error.generate(x_al, x_ac)

    assert abs(wt.data - expected['wt']['one']['wt']).mean() < 1e-16
    assert abs(wet_tropo.data -
               expected['wt']['one']['wet_tropo']).mean() < 1e-16

    assert abs(wt_nadir.data -
               expected['wt_nadir']['one']['wt']).mean() < 1e-16
    assert abs(wet_tropo_nadir.data -
               expected['wt_nadir']['one']['wet_tropo']).mean() < 1e-16

    parameters.nbeam = 2
    error = swot_simulator.error.wt.WetTroposphere(parameters)
    wt, wet_tropo, wt_nadir, wet_tropo_nadir = error.generate(x_al, x_ac)

    assert abs(wt.data - expected['wt']['two']['wt']).mean() < 1e-16
    assert abs(wet_tropo.data -
               expected['wt']['two']['wet_tropo']).mean() < 1e-16

    assert abs(wt_nadir.data -
               expected['wt_nadir']['two']['wt']).mean() < 1e-16
    assert abs(wet_tropo_nadir.data -
               expected['wt_nadir']['two']['wet_tropo']).mean() < 1e-16


def test_altimeter():
    parameters = load_parameters()
    (x_al, _, _, _, expected) = load_data('altimeter')

    error = swot_simulator.error.altimeter.Altimeter(parameters)
    generated = error.generate(x_al)
    assert abs(generated.data - expected).mean() < 1e-17
