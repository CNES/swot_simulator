# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Mathematical routines
---------------------
"""
from typing import Optional, Tuple
import collections
import numpy as np
import numba as nb


@nb.njit(cache=True, nogil=True)
def normalize_longitude(lon: np.ndarray,
                        lon_min: Optional[float] = -180.0) -> np.ndarray:
    """Normalize longitudes between [lon_min, lon_min + 360]"""
    return ((lon - lon_min) % 360) + lon_min


@nb.njit('UniTuple(float64[:], 3)(float64[:], float64[:])',
         cache=True,
         nogil=True)
def spher2cart(lon: np.ndarray,
               lat: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert spherical coordinates to cartesian coordinates"""
    rlon = np.radians(lon)
    rlat = np.radians(lat)
    x = np.cos(rlon) * np.cos(rlat)
    y = np.sin(rlon) * np.cos(rlat)
    z = np.sin(rlat)
    return x, y, z


@nb.njit('UniTuple(float64[:], 2)(float64[:], float64[:], float64[:])',
         cache=True,
         nogil=True)
def cart2spher(x: np.ndarray, y: np.ndarray,
               z: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Convert cartesian coordinates to spherical coordinates"""
    indexes = np.where((x == 0) & (y == 0))[0]
    if indexes.size:
        x[indexes] = np.nan
        y[indexes] = np.nan
    lat = np.arctan2(z, np.sqrt(x * x + y * y))
    lon = np.arctan2(y, x)
    if indexes.size:
        lon[indexes] = 0
        lat[indexes] = np.pi * 0.5 * np.sign(z[indexes])
    return np.degrees(lon), np.degrees(lat)


@nb.njit('float64[:, ::1](float64, float64[::1])', cache=True, nogil=True)
def rotation_3d_matrix(theta: float, axis: np.ndarray) -> np.ndarray:
    """Creates a rotation matrix: Slow method.

    Inputs are rotation angle theta and rotation axis axis.
    The rotation matrix correspond to a rotation of angle theta
    with respect to axis axis."""
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta * 0.5)
    b, c, d = -axis * np.sin(theta * 0.5)
    result = np.empty((3, 3), dtype=axis.dtype)
    result[0, 0] = a * a + b * b - c * c - d * d
    result[0, 1] = 2 * (b * c - a * d)
    result[0, 2] = 2 * (b * d + a * c)
    result[1, 0] = 2 * (b * c + a * d)
    result[1, 1] = a * a + c * c - b * b - d * d
    result[1, 2] = 2 * (c * d - a * b)
    result[2, 0] = 2 * (b * d - a * c)
    result[2, 1] = 2 * (c * d + a * b)
    result[2, 2] = a * a + d * d - b * b - c * c
    return result


@nb.njit('float64[::1](float64[::1], float64[::1], float64)',
         cache=True,
         nogil=True)
def curvilinear_distance(lon: np.ndarray, lat: np.ndarray,
                         radius: float) -> np.ndarray:
    """Calculating the curvilinear distance."""
    lat1 = np.radians(lat[0:-1])
    lat2 = np.radians(lat[1:])
    dlon = np.radians(np.diff(lon))

    distance = np.zeros_like(lon)
    distance[1:] = np.arccos(
        np.sin(lat1) * np.sin(lat2) +
        np.cos(lat1) * np.cos(lat2) * np.cos(dlon)) * radius
    return np.cumsum(distance)


@nb.njit('float64[:, :](float64[:, :])', cache=True, nogil=True)
def satellite_direction(location: np.ndarray) -> np.ndarray:
    """Calculate satellite direction"""
    direction = np.empty_like(location)

    denominator = np.sqrt(location[1:-1, 0]**2 + location[1:-1, 1]**2 +
                          location[1:-1, 2]**2)
    direction[1:-1, 0] = (location[2:, 0] - location[:-2, 0]) / denominator
    direction[1:-1, 1] = (location[2:, 1] - location[:-2, 1]) / denominator
    direction[1:-1, 2] = (location[2:, 2] - location[:-2, 2]) / denominator

    direction[0, :] = direction[1, :]
    direction[0, :] = direction[1, :]
    direction[0, :] = direction[1, :]
    direction[-1, :] = direction[-2, :]
    direction[-1, :] = direction[-2, :]
    direction[-1, :] = direction[-2, :]

    return direction


@nb.njit('UniTuple(float64, 2)(float64, float64, float64)',
         cache=True,
         nogil=True)
def __cartesian2spherical(x: float, y: float, z: float) -> Tuple[float, float]:
    """Convert cartesian coordinates to spherical coordinates for scalar
    values. This function is for internal use only"""
    if x == 0 and y == 0:
        return 0, np.degrees(np.pi * 0.5 * np.sign(z))
    lat = np.arctan2(z, np.sqrt(x * x + y * y))
    lon = np.arctan2(y, x)
    return np.degrees(lon), np.degrees(lat)


@nb.njit('UniTuple(float64[:, ::1], 2)(float64, float64, int64, float64, '
         'float64[:, ::1], float64[:, ::1])',
         cache=True,
         nogil=True)
def calculate_swath(delta_ac: float, half_gap: float, half_swath: int,
                    radius: float, location: np.ndarray,
                    direction: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Calculation of swath geometry"""
    lon = np.empty((location.shape[0], 2 * half_swath), dtype=np.float64)
    lat = np.empty((location.shape[0], 2 * half_swath), dtype=np.float64)

    for ix in range(len(location)):
        for jx in range(0, int(half_swath)):
            rotation = rotation_3d_matrix(-(jx * delta_ac + half_gap) / radius,
                                          direction[ix, :])

            loc = np.dot(rotation, location[ix])
            kx = half_swath + jx
            lon[ix,
                kx], lat[ix,
                         kx] = __cartesian2spherical(loc[0], loc[1], loc[2])

            loc = np.dot(rotation.T, location[ix])
            kx = half_swath - jx - 1
            lon[ix,
                kx], lat[ix,
                         kx] = __cartesian2spherical(loc[0], loc[1], loc[2])
    return lon, lat


Point = collections.namedtuple('Point', ['lon', 'lat'])
Point.__doc__ = """Defines a 2D Point"""


class Box:
    """Defines a box made of two describing points.
    """
    def __init__(self,
                 min_corner: Optional[Point] = None,
                 max_corner: Optional[Point] = None):
        self.min_corner = min_corner or Point(-180, -90)
        self.max_corner = max_corner or Point(180, 90)

    def __repr__(self):
        return "%s(%s, %s)" % (self.__class__.__name__, str(
            self.min_corner), str(self.max_corner))

    def within(self, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
        """Returns true if the coordinates are located in the box."""
        lon = normalize_longitude(lon, lon_min=self.min_corner.lon)
        lon_is_in_range = (lon >= self.min_corner.lon) | (
            lon <= self.max_corner.lon
        ) if self.max_corner.lon < self.min_corner.lon else (
            lon >= self.min_corner.lon) & (lon <= self.max_corner.lon)

        return (lat >= self.min_corner.lat) & (
            lat <= self.max_corner.lat) & lon_is_in_range
