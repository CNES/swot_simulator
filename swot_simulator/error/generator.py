# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Generate instrumental errors
----------------------------
"""
from typing import Dict, Union
import dask.distributed
import numpy as np
from . import (Altimeter, BaselineDilation, CorrectedRollPhase, Karin, Orbital,
               RollPhase, Timing, WetTroposphere)
from .. import random_signal
from .. import settings


class Generator:
    """Instrumental error generator.

    Args:
        parameters (settings.Parameters): Simulation settings
        first_date (numpy.datetime64): Date of the first simulated
            measurement.
        orbit_duration (numpy.timedelta64, optional): Orbit duration.
    """
    def __init__(self, parameters: settings.Parameters,
                 first_date: np.datetime64, orbit_duration: np.timedelta64):
        #: The list of user-defined error generators
        self.generators = []

        assert parameters.error_spectrum is not None
        error_spectrum = random_signal.read_file_instr(
            parameters.error_spectrum, parameters.delta_al)

        for item in parameters.noise:
            if item == Altimeter.__name__:
                self.generators.append(Altimeter(parameters))
            elif item == BaselineDilation.__name__:
                self.generators.append(
                    BaselineDilation(parameters,
                                     error_spectrum['dilationPSD'].data,
                                     error_spectrum['spatial_frequency'].data))
            elif item == CorrectedRollPhase.__name__:
                self.generators.append(
                    CorrectedRollPhase(parameters, first_date))
            elif item == Karin.__name__:
                self.generators.append(Karin(parameters))
            elif item == Orbital.__name__:
                self.generators.append(Orbital(parameters, orbit_duration))
            elif item == RollPhase.__name__:
                self.generators.append(
                    RollPhase(parameters, error_spectrum['rollPSD'].data,
                              error_spectrum['gyroPSD'].data,
                              error_spectrum['phasePSD'].data,
                              error_spectrum['spatial_frequency'].data))
            elif item == Timing.__name__:
                self.generators.append(
                    Timing(parameters, error_spectrum['timingPSD'].data,
                           error_spectrum['spatial_frequency'].data))
            elif item == WetTroposphere.__name__:
                self.generators.append(WetTroposphere(parameters))
            else:
                # A new error generation class has been implemented but it is
                # not handled by this object.
                raise ValueError(f"unknown error generation class: {item}")

    def generate(self, cycle_number: int, pass_number: int,
                 curvilinear_distance: float, time: np.ndarray,
                 x_al: np.ndarray, x_ac: np.ndarray,
                 swh: Union[np.ndarray, float]) -> Dict[str, np.ndarray]:
        """Generate errors

        Args:
            cycle_number (int): Cycle number.
            pass_number (int): Pass number.
            curvilinear_distance (float): Curvilinear distance covered by the
                satellite during a complete cycle.
            time (numpy.ndarray): Date of measurements.
            x_al (numpy.ndarray): Along track distance.
            x_ac (numpy.ndarray): Across track distance.
            swh (numpy.ndarray, float): Significant wave height. Used to
                modulate instrumental noise as a function of sea state.

        Returns:
            dict: Associative array between error variables and simulated
            values.
        """
        # The variable "x_al" contains the distance a along track from the
        # beginning of the cycle. In order not to reproduce the same error from
        # one cycle to another, we transform it into the distance covered since
        # the first cycle.
        x_al += curvilinear_distance * (cycle_number - 1)

        result = {}
        if not self.generators or x_al.shape[0] == 0:
            return result

        futures = []
        with dask.distributed.worker_client() as client:
            for item in self.generators:
                if isinstance(item, Altimeter):
                    futures.append(client.submit(item.generate, x_al))
                elif isinstance(item, BaselineDilation):
                    futures.append(client.submit(item.generate, x_al, x_ac))
                elif isinstance(item, CorrectedRollPhase):
                    futures.append(client.submit(item.generate, time, x_ac))
                elif isinstance(item, Karin):
                    if isinstance(swh, float):
                        swh = np.full((x_ac.size, ), swh)
                    futures.append(
                        client.submit(item.generate,
                                      cycle_number * 10000 + pass_number, x_al,
                                      x_ac, swh))
                elif isinstance(item, Orbital):
                    futures.append(client.submit(item.generate, time, x_ac))
                elif isinstance(item, RollPhase):
                    futures.append(client.submit(item.generate, x_al, x_ac))
                elif isinstance(item, Timing):
                    futures.append(client.submit(item.generate, x_al, x_ac))
                elif isinstance(item, WetTroposphere):
                    futures.append(client.submit(item.generate, x_al, x_ac))

            for future in dask.distributed.as_completed(futures):
                result.update(future.result())
        return result
