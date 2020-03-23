import numpy as np
from . import utils
import sys
import xarray as xr
import logging
from typing import IO, Optional


# Define logging level for debug purposes
logger = logging.getLogger(__name__)

VOLUMETRIC_MEAN_RADIUS = 6371


def read_roll_phase(roll_phase_file:str, first_date: np.datetime64):
    ''' Retrieve estimated roll phase from roll phase file provided by CLS'''
    #TODO proof errors if issue reading the file
    try:
        ds = xr.open_dataset(roll_phase_file)
        time_date = first_date + (ds.time[:].astype(np.float32) * 86400
                                     * 1000000).astype("timedelta64[us]")
        time_date = time_date.astype('datetime64[us]').astype('float64') / 1000
        remaining1 = ds.slope1_err[:] - ds.slope1_est[:]
        remaining2 = ds.slope2_err[:] - ds.slope2_est[:]
    except IOError:
        logger.error('There was an error opening the file'
                     ' {}'.format(roll_phase_file))
        return None
    return ds, time_date, remaining1, remaining2


class error_stat():
    def __init__(self, ds,  # delta_al:float, lambda_max:float,
                 ) -> None:
        # Get baseline dilation power spectrum
        self.PSroll = ds.rollPSD.data
        self.PSphase = ds.phasePSD.data
        self.freq = ds.spatial_frequency.data

    def init_error_savesignal(self, delta_al:float, lambda_max:float,
                              npseudoper: int, len_repeat: int,
                              nseed: Optional[int]=0)-> None:
        """Compute random coefficients using the power spectrum """
        gencoef = utils.gen_rcoeff_signal1d(self.freq, self.PSroll,
                                            2 * delta_al, lambda_max,
                                            npseudoper, len_repeat, nseed)
        self.A_roll, self.phi_roll = gencoef
        gencoef = utils.gen_rcoeff_signal1d(self.freq, self.PSphase,
                                            2 * delta_al, lambda_max,
                                            npseudoper, len_repeat, nseed)
        self.A_phase_l, self.phi_phase_l = gencoef
        gencoef = utils.gen_rcoeff_signal1d(self.freq, self.PSphase,
                                            2 * delta_al, lambda_max,
                                            npseudoper, len_repeat, nseed)

        self.A_phase_r, self.phi_phase_r = gencoef

    def init_error_gensignal(self, ncomp1d: int,
                             nseed: Optional[int]=0)-> None:
        """Compute random signal using the power spectrum """
        gencoef = utils.gen_coeff_signal1d(self.freq, self.PSroll, ncomp1d,
                                           nseed)
        self.A_roll, self.phi_roll, self.fr_roll = gencoef
        gencoef = utils.gen_coeff_signal1d(self.freq, self.PSphase, ncomp1d,
                                           nseed)
        self.A_phase_l, self.phi_phase_l, self.fr_phase_l = gencoef
        gencoef = utils.gen_coeff_signal1d(self.freq, self.PSphase, ncomp1d,
                                           nseed)
        self.A_phase_r, self.phi_phase_r, self.fr_phase_r = gencoef


    def make_error(self, time: np.array, x_al: np.array, al_cycle: float,
                   cycle: int,  dal: float,
                   npseudoper: int, sat_const: dict, len_repeat:float,
                   lmax: Optional[float]=20000,
                   savesignal: Optional[bool]=True,
                   roll_phase_file: Optional[str]=None,
                   first_date: Optional[np.datetime64]=None):
        Rearth = VOLUMETRIC_MEAN_RADIUS
        sat_elev = sat_const['height']
        Fka = sat_const['Fka']
        _to_km = (1 / (Fka * 2*np.pi / sat_const['C'] * sat_const['B'])
                  * (1 + sat_elev/Rearth)*np.pi/180. * 10**3)
        if not (roll_phase_file is None):
         # - Compute the associated phase error on the swath in m
             results = read_roll_phase(roll_phase_file, first_date)
             roll_phase, time_date, remaining1, remaining2 = results
             time = time.astype('datetime64[us]').astype('float64')
             #_time = np.mod(time, np.max(time_date))
             theta = np.interp(time, time_date, roll_phase.proll_err)
             phase_l = np.interp(time, time_date, roll_phase.p1phase_err)
             phase_r = np.interp(time, time_date, roll_phase.p2phase_err)
             phase1d = np.concatenate(([phase_l.T], [phase_r.T]), axis=0)
             roll = np.interp(time, time_date, roll_phase.proll_err)
             self.phase1d = phase1d.T * 10**(-3)
             #rem_l = np.interp(time, time_date, remaining1)
             #rem_r = np.interp(time, time_date, remaining2)

             est_l = np.interp(time, time_date, roll_phase.slope1_est)
             err_l = np.interp(time, time_date, roll_phase.slope1_err)
             rem_l = est_l - err_l
             est_r = np.interp(time, time_date, roll_phase.slope2_est)
             err_r = np.interp(time, time_date, roll_phase.slope2_err)
             rem_r = est_r - err_r
             theta2 = np.concatenate(([rem_l.T], [rem_r.T]), axis=0)
             self.rollphase_est1d = theta2.T *10**(-3)
             self.roll1d = roll[:] * 10**(-3)
        else:
            # - Compute roll angle using random coefficients or signals
            # previously initialized with the power spectrum
            # - Compute left and right phase angles using random
            #   coefficients or signals previously initialized
            if savesignal is True:
                xx = (np.float64(x_al[:]) + float(cycle * al_cycle)
                     ) % (len_repeat)
                # roll
                theta = utils.gen_signal1d(xx, self.A_roll, self.phi_roll,
                                           2 * dal, lmax, npseudoper)
                theta_l = utils.gen_signal1d(xx, self.A_phase_l,
                                             self.phi_phase_l,
                                             2 * dal, lmax, npseudoper)
                theta_r = utils.gen_signal1d(xx, self.A_phase_r,
                                             self.phi_phase_r, 2 * dal, lmax,
                                             npseudoper)
            else:
                theta_l = np.zeros((nal))
                theta_r = np.zeros((nal))
                theta = np.zeros((nal))
                for comp in range(0, self.ncomp1d):
                    phase_x_al = (2. * np.pi * float(self.fr_phase_l[comp])
                                  * (np.float64(x_al[:]) +
                                  + float(cycle * al_cycle)))%(2.*np.pi)
                    theta_l[:] = (theta_l[:] + 2 * self.A_phase_l[comp]
                                  * np.cos(phase_x_al[:]
                                  + self.phi_phase_l[comp]))
                    phase_x_al = (2. * np.pi * float(self.fr_phase_r[comp])
                                  * (np.float64(x_al[:])
                                  + float(cycle * al_cycle))) % (2.*np.pi)
                    theta_r[:] = (theta_r[:] + 2 * self.A_phase_r[comp]
                                  * np.cos(phase_x_al[:]
                                  + self.phi_phase_r[comp]))
                    phase_x_al = (2. * np.pi * float(self.fr_roll[comp])
                                  * (np.float64(sgrid.x_al[:])
                                  + float(cycle * al_cycle))) % (2.*np.pi)
                    theta[:] = (theta[:] + 2 * self.A_roll[comp]
                                * np.cos(phase_x_al[:]
                                + self.phi_roll[comp]))
            # - Compute the associated roll  error on the swath in m
            self.roll1d = ((1 + sat_elev / Rearth) * theta[:]
                           * np.pi/180./3600.) *10**3
            # - Compute the associated phase error on the swath in m
            phase1d = np.concatenate(([_to_km * theta_l[:].T],
                                     [_to_km * theta_r[:].T]), axis=0)
            self.phase1d = phase1d.T
            self.rollphase_est1d = np.zeros((np.shape(phase1d.T)))
            del theta_l, theta_r
            del theta

    def reconstruct_2D_error(self, x_ac: np.ndarray):
        ''' Reconstruct 2D errors from 1D instrumental error simulation '''
        nac = np.shape(x_ac)[0]
        mid_swath = int(nac / 2)
        ac_l = np.mat(x_ac[:mid_swath])
        ac_r = np.mat(x_ac[mid_swath:])
        phase1d = self.phase1d
        phase = np.full((np.shape(phase1d)[0], nac), np.nan)
        phase[:, :mid_swath] = np.mat(phase1d[:, 0]).T * ac_l
        phase[:, mid_swath:] = np.mat(phase1d[:, 1]).T * ac_r
        phaseroll1d = self.rollphase_est1d
        rollphase_est = np.full((np.shape(phaseroll1d)[0], nac), np.nan)
        rollphase_est[:, :mid_swath] = np.mat(phaseroll1d[:, 0]).T * ac_l
        rollphase_est[:, mid_swath:] = np.mat(phaseroll1d[:, 1]).T * ac_r
        ac = np.mat(x_ac)
        roll1d = self.roll1d
        roll = np.mat(roll1d).T * ac
        return phase, roll, rollphase_est


