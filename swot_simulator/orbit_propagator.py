# Copyright (c) 2020 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Orbit Propagator
----------------
"""
from typing import Dict, Iterable, Optional, TextIO, Tuple
import logging
import numpy as np
from . import math_func
from . import settings
from . import VOLUMETRIC_MEAN_RADIUS

LOGGER = logging.getLogger(__name__)


def load_ephemeris(
        stream: TextIO, cols: Optional[Iterable[int]] = None
) -> Tuple[Dict[str, float], Tuple]:
    """Loads a tabular file describing a satellite orbit."""
    if cols is None:
        cols = (1, 2, 0)
    LOGGER.info("loading ephemeris:  %s", stream.name)

    comments = []
    data = []
    for item in stream:
        item = item.strip()
        comments.append(item[1:]) if item.startswith("#") \
            else data.append(list(float(value) for value in item.split()))
    data = np.asarray(data, dtype=np.float64)

    def to_dict(comments) -> Dict[str, float]:
        """Returns a dictionary describing the parameters of the orbit."""
        result = dict()
        for item in comments:
            key, value = item.split("=")
            result[key.strip()] = float(value)
        return result

    return to_dict(comments), tuple(data[:, item] for item in cols)


def interpolate(lon: np.ndarray, lat: np.ndarray,
                time: np.ndarray) -> np.ndarray:
    """Interpolate the given orbit at high resolution (0.5 seconds)"""
    x, y, z = math_func.spher2cart(lon, lat)
    time_hr = np.arange(time[0], time[-1], 0.5, dtype=time.dtype)

    x = np.interp(time_hr, time, x)
    y = np.interp(time_hr, time, y)
    z = np.interp(time_hr, time, z)

    lon, lat = math_func.cart2spher(x, y, z)

    return lon, lat, time_hr


def rearrange_orbit(cycle_duration: float, lon: np.ndarray, lat: np.ndarray,
                    time: np.ndarray
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Rearrange orbit starting from pass 1

    Detect the beginning of pass 1 in the ephemeris. By definition, it is
    the first passage at southernmost latitude."""
    dy = np.roll(lat, 1) - lat
    indexes = np.where((dy < 0) & (np.roll(dy, 1) >= 0))[0]

    # Shift coordinates, so that the first point of the orbit is the beginning
    # of pass 1
    shift = indexes[-1]

    lon = np.hstack([lon[shift:], lon[:shift]])
    lat = np.hstack([lat[shift:], lat[:shift]])
    time = np.hstack([time[shift:], time[:shift]])
    time = (time - time[0]) % cycle_duration
    if time[np.where(time < 0)]:
        LOGGER.warn('there are negative times in your orbit')

    # Compute the initial time of each pass
    dy = np.roll(lat, 1) - lat
    indexes = np.where(((dy < 0) & (np.roll(dy, 1) >= 0))
                       | ((dy > 0)
                          & (np.roll(dy, 1) <= 0)))
    return lon, lat, time, time[indexes[0]]


def select_box(box: math_func.Box, lon: np.ndarray, lat: np.ndarray,
               time: np.ndarray
               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Selects the orbit in the defined box."""
    mask = box.within(lon, lat)
    x_al = math_func.curvilinear_distance(lon, lat, VOLUMETRIC_MEAN_RADIUS)
    return lon[mask], lat[mask], time[mask], x_al[mask]


class Orbit:
    """Properties of one orbit"""
    def __init__(self, height: float, lat: np.ndarray, lon: np.ndarray,
                 pass_time: np.ndarray, time: np.ndarray, x_al: np.ndarray,
                 al_cycle, shift_time: Optional[float]):
        self.height = height
        self.lat = lat
        self.lon = lon
        self.pass_time = pass_time
        self.shift_time = shift_time
        self.time = time
        self.x_al = x_al
        self.al_cycle = al_cycle

    def passes_per_cycle(self):
        return len(self.pass_time)

    def pass_duration(self, number: int) -> np.timedelta64:
        if number == self.passes_per_cycle():
            return np.timedelta64(
                int((self.time[-1] - self.pass_time[-1]) * 1e6), 'us')
        return np.timedelta64(
            int((self.pass_time[number] - self.pass_time[number - 1]) * 1e6),
            'us')

    def decode_absolute_pass_number(self, number: int) -> Tuple[int, int]:
        """Calculate the cycle, pass  number from a given absolute pass number.

        Return:
            tuple: cycle, track number
        """
        number -= 1
        return (int(number / self.passes_per_cycle()) + 1,
                (number % self.passes_per_cycle()) + 1)


class Pass:
    """Handle one pass"""
    def __init__(self, lat_nadir: np.ndarray, lat: np.ndarray,
                 lon_nadir: np.ndarray, lon: np.ndarray, time: np.ndarray,
                 x_ac: np.ndarray, x_al: np.ndarray, al_cycle:float):
        self.lat_nadir = lat_nadir
        self.lat = lat
        self.lon_nadir = lon_nadir
        self.lon = lon
        self.timedelta = ((time * 1e6).astype(
            np.int64)).astype("timedelta64[us]")
        self._time = None
        self.x_ac = x_ac
        self.x_al = x_al
        self.al_cycle = al_cycle

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, date: np.datetime64):
        self._time = self.timedelta + date


def calculate_orbit(parameters: settings.Parameters,
                    ephemeris: TextIO) -> Orbit:
    """Computes the orbit nadir on a subdomain.

    The path of the satellite is given by the orbit file and the subdomain
    corresponds to the one in the model. Note that a subdomain can be manually
    added in the parameters file."""
    properties, (lon, lat,
                 time) = load_ephemeris(ephemeris,
                                        cols=parameters.ephemeris_cols)

    if np.any(lon > 360) or np.any(np.abs(lat) > 90):
        raise RuntimeError("The definition of ephemeris is incorrect. Is the "
                           "'ephemeris_cols' parameter well defined?")

    # Put longitude in [-180, 180[
    lon = math_func.normalize_longitude(lon)

    for key, value in properties.items():
        if key == 'cycle_duration':
            parameters.cycle_duration = value
        elif key == 'height':
            parameters.height = value
        else:
            raise RuntimeError(f"The {key!r} parameter defined in the "
                               "ephemeris is not handled.")

    # If orbit is at low resolution, interpolate the orbit provided
    if np.mean(np.diff(time)) > 0.5:
        lon, lat, time = interpolate(lon, lat, time)

    # Cut orbit if more than an orbit cycle is provided
    indexes = np.where(time < parameters.cycle_duration * 86400)[0]
    lon = lon[indexes]
    lat = lat[indexes]
    time = time[indexes]

    # Get cycle period.
    cycle_duration = time[-1] + time[1] - time[0]

    # shift time if the user needs to shift the time of the orbit
    if parameters.shift_time is not None:
        indexes = np.where(time >= parameters.shift_time)[0]
        lon = np.hstack([lon[indexes[0]:], lon[:indexes[0]]])
        lat = np.hstack([lat[indexes[0]:], lat[:indexes[0]]])
    del indexes

    if parameters.shift_lon is not None:
        lon = math_func.normalize_longitude(lon + parameters.shift_lon)

    # Rearrange orbit starting from pass 1
    lon, lat, time, pass_time = rearrange_orbit(cycle_duration, lon, lat, time)

    # Get the selected domain
    lon, lat, time, distance = select_box(parameters.box, lon, lat, time)

    # Interpolate the final orbit according the given along track resolution
    x_al = np.arange(distance[0],
                     distance[-2],
                     parameters.delta_al,
                     dtype=distance.dtype)
    time = np.interp(x_al, distance[:-1], time[:-1])
    lon = np.interp(x_al, distance[:-1], lon[:-1])
    lat = np.interp(x_al, distance[:-1], lat[:-1])
    al_cycle = distance[-1]
    return Orbit(parameters.height, lat, lon, np.sort(pass_time), time, x_al,
                 al_cycle, parameters.shift_time)


def calculate_pass(pass_number: int, orbit: Orbit,
                   parameters: settings.Parameters) -> Optional[Pass]:
    index = pass_number - 1
    # Selected indexes corresponding to the current pass
    if index == len(orbit.pass_time) - 1:
        indexes = np.where(orbit.time >= orbit.pass_time[-1])[0]
    else:
        indexes = np.where((orbit.time >= orbit.pass_time[index])
                           & (orbit.time < orbit.pass_time[index + 1]))[0]

    if len(indexes) < 5:
        return

    lon_nadir = orbit.lon[indexes]
    lat_nadir = orbit.lat[indexes]
    time = orbit.time[indexes]
    x_al = orbit.x_al[indexes]
    al_cycle = orbit.al_cycle

    # Compute accross track distances from nadir
    # Number of points in half of the swath
    half_swath = int((parameters.half_swath - parameters.half_gap) /
                     parameters.delta_ac) + 1
    x_ac = np.arange(half_swath) * parameters.delta_al + parameters.half_gap
    x_ac = np.hstack((-np.flip(x_ac), x_ac))

    location = np.ascontiguousarray(
        np.vstack(math_func.spher2cart(lon_nadir, lat_nadir)).T)
    satellite_direction = math_func.satellite_direction(location)
    lon, lat = math_func.calculate_swath(parameters.delta_ac, parameters.half_gap,
                                    half_swath, VOLUMETRIC_MEAN_RADIUS,
                                    location, satellite_direction)

    return Pass(lat_nadir, lat, lon_nadir, lon, time, x_ac, x_al, al_cycle)
