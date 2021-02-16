# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Estimated roll errors
---------------------
"""
from typing import Dict, Tuple
import numpy as np
import xarray as xr

from .. import settings


class CorrectedRollPhase:
    """
    Corrected roll errors

    Args:
        parameters (settings.Parameters): Simulation settings
        first_date (numpy.datetime64): Date of the first simulated
            measurement.
    """
    def __init__(self, parameters: settings.Parameters,
                 first_date: np.datetime64) -> None:
        with xr.open_dataset(parameters.corrected_roll_phase_dataset) as ds:
            time_date = first_date + (ds.time[:].astype(np.float32) * 86400 *
                                      1000000).astype("timedelta64[us]")
            self.time_date = time_date.astype('datetime64[us]').astype(
                'float64') * 0.001
            self.proll_err = ds.proll_err
            self.p1phase_err = ds.p1phase_err
            self.p2phase_err = ds.p2phase_err
            self.slope1_est = ds.slope1_est
            self.slope2_est = ds.slope2_est
            self.slope1_err = ds.slope1_err
            self.slope2_err = ds.slope2_err

    def _generate_1d(
            self,
            time: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        time = time.astype('datetime64[us]').astype('float64')
        phase_l = np.interp(time, self.time_date, self.p1phase_err)
        phase_r = np.interp(time, self.time_date, self.p2phase_err)
        phase1d = np.concatenate(([phase_l.T], [phase_r.T]), axis=0)
        roll = np.interp(time, self.time_date, self.proll_err)

        est_l = np.interp(time, self.time_date, self.slope1_est)
        err_l = np.interp(time, self.time_date, self.slope1_err)
        rem_l = est_l - err_l
        est_r = np.interp(time, self.time_date, self.slope2_est)
        err_r = np.interp(time, self.time_date, self.slope2_err)
        rem_r = est_r - err_r
        theta2 = np.concatenate(([rem_l.T], [rem_r.T]), axis=0)

        return roll * 1e-3, phase1d.T * 1e-3, theta2.T * 1e-3

    def generate(
        self,
        time: np.ndarray,
        x_ac: np.ndarray,
    ) -> Dict[str, np.ndarray]:
        """Interpolate roll and phase and errors

        Args:
            time (numpy.ndarray): Date of measurements.
            x_ac (numpy.ndarray): Across track distance.

        Returns:
            dict: variable name and errors simulated.
        """
        roll_1d, phase_1d, rollphase_est_1d = self._generate_1d(time)
        num_pixels = x_ac.shape[0]
        swath_center = num_pixels // 2
        ac_l = x_ac[:swath_center]
        ac_r = x_ac[swath_center:]

        phase = np.full((phase_1d.shape[0], num_pixels), np.nan)
        phase[:, :swath_center] = ac_l * phase_1d[:, 0, np.newaxis]
        phase[:, swath_center:] = ac_r * phase_1d[:, 1, np.newaxis]

        rollphase_est = np.full((phase_1d.shape[0], num_pixels), np.nan)
        rollphase_est[:, :swath_center] = np.mat(rollphase_est_1d[:,
                                                                  0]).T * ac_l
        rollphase_est[:,
                      swath_center:] = np.mat(rollphase_est_1d[:, 1]).T * ac_r
        return {
            "simulated_error_roll": x_ac * roll_1d[:, np.newaxis],
            "simulated_error_phase": phase,
            "simulated_roll_phase_estimate": rollphase_est
        }
