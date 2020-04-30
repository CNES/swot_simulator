from . import error_altimeter
from . import utils
import numpy as np


def make_error(x_al:np.ndarray, al_cycle:float, cycle_number:int, dal:int,
               par_err:dict)->dict:
    dic_err = {}
    err = error_altimeter.error_stat(par_err['ncomp1d'], dal)
    err.make_error(x_al, dal, par_err['len_repeat'], nseed=par_err['nseed'])
    dic_err['err_altimeter'] = err.nadir
    return dic_err

def add_error(dic_err:dict, ssh_true:np.ndarray) -> np.ndarray:
    ssh = + ssh_true
    for key in dic_err.keys():
        ssh = ssh + dic_err[key]
    return ssh
