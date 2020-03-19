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
    ds = utils.read_file_instr(par_err['error_spectrum_file'], dal,
                               par_err['lambda_max'])
    bd = error_baselinedilation.error_stat(ds)
    rp = error_roll_phase.error_stat(ds)
    karin = error_karin.error_stat(par_err['karin_file'],
                                   par_err['nrand_karin'])
    ti = error_timing.error_stat(ds)
    wt = error_wt.error_stat(dal, par_err['ncomp2d'])
    if par_err['save_signal'] is True:
        bd.init_error_savesignal(dal, par_err['lambda_max'], par_err['npseudoper'], par_err['len_repeat'])
        rp.init_error_savesignal(dal, par_err['lambda_max'], par_err['npseudoper'], par_err['len_repeat'])
        ti.init_error_savesignal(dal, par_err['lambda_max'], par_err['npseudoper'], par_err['len_repeat'])
    else:
        bd.init_error_gensignal(par_err['ncomp1d'])
        rp.init_error_gensignal(par_err['ncomp1d'])
        ti.init_error_gensignal(par_err['ncomp1d'])
    wt.init_error_gensignal(par_err['ncomp2d'])
    karin.init_error(np.shape(x_ac)[0])
    bd.make_error(x_al, al_cycle, cycle_number, dal, par_err['npseudoper'], sat_const,
                  par_err['len_repeat'], lmax=par_err['lambda_max'], savesignal=par_err['save_signal'])
    rp.make_error(time, x_al, al_cycle, cycle_number, dal, par_err['npseudoper'],
                  sat_const, par_err['len_repeat'], lmax=par_err['lambda_max'],
                  savesignal=par_err['save_signal'],
                  roll_phase_file=par_err['roll_phase_file'],
                  first_date=first_date)
    ti.make_error(x_al, al_cycle, cycle_number, dal, par_err['npseudoper'], sat_const,
                  par_err['len_repeat'], lmax=par_err['lambda_max'], savesignal=par_err['save_signal'])
    wt.make_error(x_al, x_ac, al_cycle, cycle_number, dal, dac,
                  par_err['nbeam'], par_err['sigma'],
                  par_err['beam_position'], sat_const)

    karin.make_error(x_al, x_ac, al_cycle, cycle_number, dal, dac,
                     par_err['swh'])
    dic_err['baseline_dilation'] = bd.reconstruct_2D_error(x_ac)
    result = rp.reconstruct_2D_error(x_ac)
    dic_err['phase'], dic_err['roll'], dic_err['roll_phase_est'] = result
    dic_err['karin'] = karin.karin
    dic_err['timing'] = ti.reconstruct_2D_error(x_ac)
    result = wt.reconstruct_2D_error()
    dic_err['wt'], dic_err['wet_tropo2'], var = result

    return dic_err

