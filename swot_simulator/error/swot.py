import xarray as xr
import numpy as np
from . import baseline_dilation
from . import roll_phase
from . import karin
from . import timing
from . import wt
from . import utils

#from swot_simulator.lib.error import utils as utils

sat_const = {}
# to retrieve from parameter file
sat_const['height'] = 891  #*10**3
# - Baseline (m)
sat_const['B'] = 10
# (in Hz)
sat_const['Fka'] = 35.75 * 10**9
# Number of days of one SWOT cycle
sat_const['tcycle'] = 20.86455
# - Light speed (m/s)
sat_const['C'] = 2.998 * 10**8
# swh = 2
#TODO read swh


def make_error(x_al: np.array, x_ac: np.array, al_cycle: float, time: np.array,
               cycle_number: int, delta_al: float, delta_ac: float,
               first_date: np.datetime64, par_err: dict) -> dict:
    dic_err = {}
    ns = par_err['nseed']
    # KaRIN random error
    err = karin.ErrorStat(par_err['karin_file'], par_err['nrand_karin'])
    err.init_error(np.shape(x_ac)[0], ns)
    err.make_error(x_al, x_ac, al_cycle, cycle_number, delta_al, delta_ac,
                   par_err['swh'])
    dic_err['karin'] = err.karin
    ds = utils.read_file_instr(par_err['error_spectrum_file'], delta_al,
                               par_err['len_repeat'])
    err = baseline_dilation.ErrorStat(ds)
    err.make_error(x_al, delta_al, sat_const, par_err['len_repeat'], nseed=(ns + 1))
    dic_err['baseline_dilation'] = err.reconstruct_2d_error(x_ac)
    err = roll_phase.ErrorStat(ds)
    err.make_error(time,
                   x_al,
                   delta_al,
                   sat_const,
                   par_err['len_repeat'],
                   nseed=(ns + 2),
                   roll_phase_file=par_err['roll_phase_file'],
                   first_date=first_date)
    result = err.reconstruct_2d_error(x_ac)
    dic_err['phase'], dic_err['roll'], dic_err['roll_phase_est'] = result
    err = timing.ErrorStat(ds)
    err.make_error(x_al, delta_al, sat_const, par_err['len_repeat'], nseed=(ns + 3))
    dic_err['timing'] = err.reconstruct_2d_error(x_ac)
    err = wt.ErrorStat(delta_al)
    err.make_error(x_al,
                   x_ac,
                   delta_al,
                   delta_ac,
                   par_err['nbeam'],
                   par_err['sigma'],
                   par_err['beam_position'],
                   sat_const,
                   par_err['len_repeat'],
                   nseed=(ns + 4),
                   lac_max=500)

    dic_err['wt'], dic_err['wet_tropo2'], _ = err.reconstruct_2d_error()
    return dic_err


def add_error(dic_err: dict, ssh_true: np.array) -> np.array:
    ssh = +ssh_true
    list_key_error = [
        'timing', 'karin', 'roll_phase_est', 'wet_tropo2', 'baseline_dilation'
    ]
    for key in dic_err.keys():
        if key in list_key_error:
            ssh = ssh + dic_err[key]
    return ssh
