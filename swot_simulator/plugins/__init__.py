# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Plug-in handlers
----------------
"""
import numpy as np


class Interface:
    """Interface of a plugin"""
    @classmethod
    def interpolate(cls, lon: np.ndarray, lat: np.ndarray,
                    dates: np.ndarray) -> np.ndarray:
        """Interpolate the geophysical field for the given coordinates"""
        raise RuntimeError("You must register a plugin")


class Puppet(Interface):
    """Interpolation routine used for testing."""
    def interpolate(self, lon, _lat, _time):
        """Returns an array filled with zeros"""
        return np.full_like(lon, 0)


class Plugin:
    """Plug-in to interpolate a geophysical field"""
    def __init__(self):
        self.plugin = Interface()

    def interpolate(self, lon: np.ndarray, lat: np.ndarray,
                    time: np.ndarray) -> np.ndarray:
        """Interpolate the geophysical field for the given coordinates."""
        return self.plugin.interpolate(lon, lat, time)

    @classmethod
    def register(cls, plugin):
        """Register the user plugin"""
        if not isinstance(plugin, Interface):
            raise TypeError("plugin must be a sub-class of "
                            f"{Interface.__class__.__name__}")  # type: ignore
        result = cls()
        result.plugin = plugin
        return result
