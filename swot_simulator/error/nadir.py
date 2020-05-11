# from typing import Dict
# import numpy as np
# from . import altimeter

# from .. import settings


# def make_error(x_al: np.ndarray, parameters: settings.Parameters) -> Dict[str, np.ndarray]:
#     err = altimeter.ErrorStat(parameters)
#     return {'err_altimeter': err.make_error(x_al)}


# def add_error(dic_err: Dict[str, np.ndarray],
#               ssh_true: np.ndarray) -> np.ndarray:
#     ssh = ssh_true.copy()
#     for key in dic_err:
#         ssh += dic_err[key]
#     return ssh
