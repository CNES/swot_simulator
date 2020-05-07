import numpy as np

from . import ErrorStat as Base
from . import utils
from .. import settings


class ErrorStat(Base):
    """Class errornadir defines the error on the nadir."""
    def __init__(self, parameters: settings.Parameters) -> None:
        super().__init__(parameters)
        # Define the sepctrum of the nadir instrument error
        freq = np.arange(1 / 3000, 1 / (2 * parameters.delta_al), 1 / 3000)
        psd = 8 + 1.05 * 10**(-4) * freq**(-2.2)
        indf = np.where(freq < 0.00023627939582672978)
        psd[indf] = 10**4

        # Convert spectrum in m2/cy
        self.psd = psd * 10**(-4)
        self.freq = freq

    def make_error(self, x_al: np.array) -> np.ndarray:
        """Build errors corresponding to each selected noise
        among the effect of the wet_tropo, and the instrumental error
        """
        # Compute random noise of 10**2 cm**2/(km/cycle)
        # Compute the correspond error on the nadir in m
        return utils.gen_signal1d(self.freq,
                                  self.psd,
                                  x_al,
                                  nseed=self.nseed,
                                  fmin=1 / self.len_repeat,
                                  fmax=1 / (2 * self.delta_al),
                                  alpha=10)
