from . import error_altimeter
from . import utils
import xarray as xr
import numpy as np


def make_error(x_al:np.array, al_cycle:float, cycle_number:int, dal:int,
               par_err:dict)->dict:
    dic_err = {}
    err = error_altimeter.error_stat(par_err['ncomp1d'], dal)
    if par_err['save_signal'] is True:
        err.init_error_savesignal(dal, par_err['lambda_max'],
                                  par_err['npseudoper'], par_err['len_repeat'],
                                  nseed=par_err['nseed'])
    else:
        err_nad.init_error_gensignal(par_err['ncomp1d'],
                                     nseed=par_err['nseed'])
    err.make_error(x_al, al_cycle, cycle_number, dal, par_err['npseudoper'],
                   par_err['len_repeat'], lmax=par_err['lambda_max'],
                   savesignal=par_err['save_signal'])
    dic_err['err_altimeter'] = err.nadir
    return dic_err

def add_error(dic_err:dict, ssh_true:np.array) -> np.array:
    ssh = + ssh_true
    for key in dic_err.keys():
        ssh = ssh + dic_err[key]
    return ssh
