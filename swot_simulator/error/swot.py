# import numpy as np
# from . import baseline_dilation
# from . import karin
# from . import roll_phase
# from . import timing
# from . import utils
# from . import wt
# from .. import settings


# def make_error(x_al: np.array, x_ac: np.array, curvilinear_distance: float,
#                time: np.array, cycle_number: int, first_date: np.datetime64,
#                parameters: settings.Parameters) -> dict:
#     dic_err = {}
#     err = karin.ErrorStat(parameters, np.shape(x_ac)[0])
#     dic_err['karin'] = err.make_error(x_al, x_ac, curvilinear_distance,
#                                       cycle_number)

#     assert parameters.error_spectrum_file is not None
#     ds = utils.read_file_instr(parameters.error_spectrum_file,
#                                parameters.delta_al, parameters.len_repeat)
#     err = baseline_dilation.ErrorStat(ds, parameters)
#     dic_err['baseline_dilation'] = err.make_error(x_al)

#     err = roll_phase.ErrorStat(ds, parameters)
#     for k, v in err.make_error(time,
#                                x_al,
#                                roll_phase_file=parameters.roll_phase_file,
#                                first_date=first_date):
#         dic_err[k] = v

#     err = timing.ErrorStat(ds, parameters)
#     dic_err['timing'] = err.reconstruct_2d_error(x_ac, err.make_error(x_al))

#     err = wt.ErrorStat(x_al, parameters)
#     errors = err.make_error(x_al, x_ac, lac_max=500)
#     dic_err['wt'], dic_err['wet_tropo2'] = errors.wt, errors.wet_tropo2

#     return dic_err


# def add_error(dic_err: dict, ssh_true: np.array) -> np.array:
#     ssh = +ssh_true
#     list_key_error = [
#         'timing', 'karin', 'roll_phase_est', 'wet_tropo2', 'baseline_dilation'
#     ]
#     for key in dic_err.keys():
#         if key in list_key_error:
#             ssh = ssh + dic_err[key]
#     return ssh
