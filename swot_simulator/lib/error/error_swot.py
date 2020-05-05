from . import error_baselinedilation
from . import error_roll_phase
from . import error_karin
from . import error_timing
from . import error_wt
from . import utils
import xarray as xr
import numpy as np

#from swot_simulator.lib.error import utils as utils

sat_const = {}
# to retrieve from parameter file
sat_const['height'] = 891 #*10**3
# - Baseline (m)
sat_const['B'] = 10
# (in Hz)
sat_const['Fka'] = 35.75 * 10**9
# Number of days of one SWOT cycle
sat_const['tcycle'] = 20.86455
# - Light speed (m/s)
sat_const['C'] = 2.998*10**8
# swh = 2
#TODO read swh


def make_error(x_al:np.array, x_ac:np.array, al_cycle:float, time:np.array,
               cycle_number:int, dal:int, dac:int, first_date: np.datetime64,
               par_err:dict)->dict:
    dic_err = {}
    ns = par_err['nseed']
    # KaRIN random error
    karin = error_karin.error_stat(par_err['karin_file'],
                                   par_err['nrand_karin'])
    karin.init_error(np.shape(x_ac)[0], ns)
    karin.make_error(x_al, x_ac, al_cycle, cycle_number, dal, dac,
                     par_err['swh'])
    dic_err['karin'] = karin.karin
    ds = utils.read_file_instr(par_err['error_spectrum_file'], dal,
                               par_err['len_repeat'])
    bd = error_baselinedilation.error_stat(ds)
    bd.make_error(x_al, dal, sat_const, par_err['len_repeat'], nseed=(ns + 1))
    dic_err['baseline_dilation'] = bd.reconstruct_2D_error(x_ac)
    rp = error_roll_phase.error_stat(ds)
    rp.make_error(time, x_al, dal, sat_const, par_err['len_repeat'],
                  nseed=(ns + 2), roll_phase_file=par_err['roll_phase_file'],
                  first_date=first_date)
    result = rp.reconstruct_2D_error(x_ac)
    dic_err['phase'], dic_err['roll'], dic_err['roll_phase_est'] = result
    ti = error_timing.error_stat(ds)
    ti.make_error(x_al, dal, sat_const, par_err['len_repeat'], nseed=(ns + 3))
    dic_err['timing'] = ti.reconstruct_2D_error(x_ac)
    wt = error_wt.error_stat(dal)
    wt.make_error(x_al, x_ac, dal, dac, par_err['nbeam'], par_err['sigma'],
                  par_err['beam_position'], sat_const, par_err['len_repeat'],
                  nseed=(ns + 4), lac_max=500)

    result = wt.reconstruct_2D_error()
    dic_err['wt'], dic_err['wet_tropo2'], var = result
    return dic_err

def add_error(dic_err:dict, ssh_true:np.array) -> np.array:
    ssh = + ssh_true
    list_key_error = ['timing', 'karin', 'roll_phase_est', 'wet_tropo2',
                      'baseline_dilation']
    for key in dic_err.keys():
        if key in list_key_error:
            ssh = ssh + dic_err[key]
    return ssh
