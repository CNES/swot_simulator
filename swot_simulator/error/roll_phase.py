from typing import Optional, Tuple
import collections
import numpy as np
import xarray as xr

from . import ErrorStat as Base
from . import utils
from .. import settings
from .. import VOLUMETRIC_MEAN_RADIUS, CELERITY, F_KA, BASELINE


def read_roll_phase(roll_phase_file: str, first_date: np.datetime64) -> Tuple:
    """Retrieve estimated roll phase from roll phase file provided by CLS"""
    ds = xr.load_dataset(roll_phase_file)
    time_date = first_date + (ds.time[:].astype(np.float32) * 86400 *
                              1000000).astype("timedelta64[us]")
    time_date = time_date.astype('datetime64[us]').astype('float64') / 1000
    remaining1 = ds.slope1_err[:] - ds.slope1_est[:]
    remaining2 = ds.slope2_err[:] - ds.slope2_est[:]
    return ds, time_date, remaining1, remaining2


#: Simulated Error
Error = collections.namedtuple('Error',
                               ['phase1d', 'rollphase_est1d', 'roll1d'])


class ErrorStat(Base):
    def __init__(self, ds: xr.Dataset,
                 parameters: settings.Parameters) -> None:
        super().__init__(parameters)
        self.nseed += 2
        # Get baseline dilation power spectrum
        self.psroll = ds.rollPSD.data
        self.psphase = ds.phasePSD.data
        self.freq = ds.spatial_frequency.data

    def make_error(self,
                   time: np.ndarray,
                   x_al: np.ndarray,
                   roll_phase_file: Optional[str] = None,
                   first_date: Optional[np.datetime64] = None) -> Error:
        # Rearth = VOLUMETRIC_MEAN_RADIUS
        _to_km = (1 / (F_KA * 2 * np.pi / CELERITY * BASELINE) *
                  (1 + self.height / VOLUMETRIC_MEAN_RADIUS) * np.pi / 180 *
                  10**3)
        if roll_phase_file is not None:
            # - Compute the associated phase error on the swath in m
            results = read_roll_phase(roll_phase_file, first_date)
            roll_phase, time_date, _, _ = results
            time = time.astype('datetime64[us]').astype('float64')

            theta = np.interp(time, time_date, roll_phase.proll_err)
            phase_l = np.interp(time, time_date, roll_phase.p1phase_err)
            phase_r = np.interp(time, time_date, roll_phase.p2phase_err)
            phase1d = np.concatenate(([phase_l.T], [phase_r.T]), axis=0)
            roll = np.interp(time, time_date, roll_phase.proll_err)
            phase1d = phase1d.T * 10**(-3)

            est_l = np.interp(time, time_date, roll_phase.slope1_est)
            err_l = np.interp(time, time_date, roll_phase.slope1_err)
            rem_l = est_l - err_l
            est_r = np.interp(time, time_date, roll_phase.slope2_est)
            err_r = np.interp(time, time_date, roll_phase.slope2_err)
            rem_r = est_r - err_r
            theta2 = np.concatenate(([rem_l.T], [rem_r.T]), axis=0)
            rollphase_est1d = theta2.T * 10**(-3)
            roll1d = roll[:] * 10**(-3)
        else:
            # - Compute roll angle using the power spectrum
            # - Compute left and right phase angles the power spectrum
            theta = utils.gen_signal_1d(self.freq,
                                       self.psroll,
                                       x_al,
                                       nseed=self.nseed,
                                       fmin=1 / self.len_repeat,
                                       fmax=1 / (2 * self.delta_al),
                                       alpha=10)
            theta_l = utils.gen_signal_1d(self.freq,
                                         self.psphase,
                                         x_al,
                                         nseed=self.nseed + 100,
                                         fmin=1 / self.len_repeat,
                                         fmax=1 / (2 * self.delta_al),
                                         alpha=10)
            theta_r = utils.gen_signal_1d(self.freq,
                                         self.psphase,
                                         x_al,
                                         nseed=self.nseed + 200,
                                         fmin=1 / self.len_repeat,
                                         fmax=1 / (2 * self.delta_al),
                                         alpha=10)
            # - Compute the associated roll  error on the swath in m
            roll1d = ((1 + self.height / VOLUMETRIC_MEAN_RADIUS) * theta[:] *
                      np.pi / 180. / 3600.) * 10**3
            # - Compute the associated phase error on the swath in m
            phase1d = np.concatenate(
                ([_to_km * theta_l[:].T], [_to_km * theta_r[:].T]), axis=0)
            phase1d = phase1d.T
            rollphase_est1d = np.zeros((np.shape(phase1d.T)))

        return Error(roll1d, phase1d, rollphase_est1d)

    def reconstruct_2d_error(self, x_ac: np.ndarray, error: Error
                             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Reconstruct 2D errors from 1D instrumental error simulation"""
        nac = np.shape(x_ac)[0]
        mid_swath = nac // 2
        ac_l = np.mat(x_ac[:mid_swath])
        ac_r = np.mat(x_ac[mid_swath:])

        phase1d = error.phase1d
        phase = np.full((np.shape(phase1d)[0], nac), np.nan)
        phase[:, :mid_swath] = np.mat(phase1d[:, 0]).T * ac_l
        phase[:, mid_swath:] = np.mat(phase1d[:, 1]).T * ac_r

        phaseroll1d = error.rollphase_est1d
        rollphase_est = np.full((np.shape(phaseroll1d)[0], nac), np.nan)
        rollphase_est[:, :mid_swath] = np.mat(phaseroll1d[:, 0]).T * ac_l
        rollphase_est[:, mid_swath:] = np.mat(phaseroll1d[:, 1]).T * ac_r

        roll = np.mat(error.roll1d).T * np.mat(x_ac)
        return phase, roll, rollphase_est
