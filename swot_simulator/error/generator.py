from typing import Dict, Optional
import dask.distributed
import numpy as np
from .. import settings
from . import (Altimeter, BaselineDilation, Karin, RollPhase, Timing,
               WetTroposphere)
from . import utils


class Generator:
    def __init__(self, parameters: settings.Parameters):
        #: The list of user-defined error generators
        self.generators = []

        assert parameters.error_spectrum is not None
        error_spectrum = utils.read_file_instr(parameters.error_spectrum,
                                               parameters.delta_al,
                                               parameters.len_repeat)

        for item in parameters.noise:
            if item == Altimeter.__name__:
                self.generators.append(Altimeter(parameters))
            elif item == BaselineDilation.__name__:
                self.generators.append(
                    BaselineDilation(parameters,
                                     error_spectrum['dilationPSD'].data,
                                     error_spectrum['spatial_frequency'].data))
            elif item == Karin.__name__:
                self.generators.append(Karin(parameters))
            elif item == RollPhase.__name__:
                self.generators.append(
                    RollPhase(parameters, error_spectrum['rollPSD'].data,
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

    def generate(self, cycle_number: int, curvilinear_distance: float,
                 x_ac: np.ndarray, x_al: np.ndarray) -> Dict[str, np.ndarray]:
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
                elif isinstance(item, Karin):
                    futures.append(
                        client.submit(item.generate, x_al, x_ac,
                                      curvilinear_distance, cycle_number))
                elif isinstance(item, RollPhase):
                    futures.append(client.submit(item.generate, x_al, x_ac))
                elif isinstance(item, Timing):
                    futures.append(client.submit(item.generate, x_al, x_ac))
                elif isinstance(item, WetTroposphere):
                    futures.append(client.submit(item.generate, x_al, x_ac))

            for future in dask.distributed.as_completed(futures):
                result.update(dict(future.result()))
        return result
