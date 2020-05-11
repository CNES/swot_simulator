from typing import Iterator, Tuple
import numpy as np

from . import utils
from .. import settings


class Altimeter:
    """Class errornadir defines the error on the nadir."""
    def __init__(self, parameters: settings.Parameters) -> None:
        # Store the generation parameters of the random signal.
        self.delta_al = 2 * parameters.delta_al
        self.nseed = parameters.nseed + 1
        self.len_repeat = parameters.len_repeat
        # Define the sepctrum of the nadir instrument error
        freq = np.arange(1 / 3000, 1 / self.delta_al, 1 / 3000)
        psd = 8 + 1.05 * 1e-4 * freq**(-2.2)
        psd[freq < 0.00023627939582672978] = 1e4

        # Convert spectrum in m2/cy
        self.psd = psd * 1e-4
        self.freq = freq

    def generate(self, x_al: np.array) -> Iterator[Tuple[str, np.ndarray]]:
        """Build errors corresponding to each selected noise
        among the effect of the wet_tropo, and the instrumental error
        """
        # Compute random noise of 10**2 cm**2/(km/cycle)
        # Compute the correspond error on the nadir in m
        error = utils.gen_signal_1d(self.freq,
                                    self.psd,
                                    x_al,
                                    nseed=self.nseed,
                                    fmin=1 / self.len_repeat,
                                    fmax=1 / self.delta_al,
                                    alpha=10)
        yield ("altimeter", error)
