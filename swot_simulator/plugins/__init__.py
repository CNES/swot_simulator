# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Plug-in handlers
----------------
"""
import abc
import os
import numpy as np
import pyinterp.backends.xarray


class Interface:
    """Interface of a plugin"""
    @classmethod
    def interpolate(cls, lon: np.ndarray, lat: np.ndarray,
                    time: np.ndarray) -> np.ndarray:
        """Interpolate the geophysical field for the given coordinates"""
        raise RuntimeError("You must register a plugin")


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
                            f"{Interface.__class__.__name__}")
        result = cls()
        result.plugin = plugin
        return result


class CartesianGridHandler(Interface):
    """Abstract class of the interpolation of a series of grids."""
    def __init__(self, path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"{path!r}")
        self.path = path
        self.ts = None
        self.dt = None
        self.load_ts()

    @abc.abstractmethod
    def load_ts(self):
        """Loading in memory the time axis of the time series"""
        ...

    @abc.abstractmethod
    def load_dataset(
            self, first_date: np.datetime64,
            last_date: np.datetime64) -> pyinterp.backends.xarray.Grid3D:
        """Loads the 3D cube describing the SSH in time and space."""
        ...
