from typing import Dict
import numpy as np
from . import altimeter


def make_error(x_al: np.ndarray, dal: float,
               par_err: dict) -> Dict[str, np.ndarray]:
    err = altimeter.ErrorStat(dal)
    err.make_error(x_al, dal, par_err['len_repeat'], nseed=par_err['nseed'])
    return {'err_altimeter': err.nadir}


def add_error(dic_err: Dict[str, np.ndarray],
              ssh_true: np.ndarray) -> np.ndarray:
    ssh = ssh_true.copy()
    for key in dic_err:
        ssh += dic_err[key]
    return ssh
