# Copyright (c) 2020 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Orbit Propagator
----------------
"""
from typing import Dict, Iterator, Optional, TextIO, Tuple
import datetime
import logging
import numpy as np
from . import math
from . import settings
from . import VOLUMETRIC_MEAN_RADIUS

LOGGER = logging.getLogger(__name__)


def load_ephemeris(stream: TextIO, cols: Optional[Tuple[int, int, int]] = None
                   ) -> Tuple[Dict[str, float], Tuple]:
    """Loads a tabular file describing a satellite orbit.

    Args:
        stream (typing.TextIO): Stream to the text file to be processed.
        cols (tuple): A tuple of 3 elements describing the columns
            representing the longitude, latitude and number of seconds elapsed
            since the beginning of the orbit.

    Returns:
        tuple: A tuple of 2 elements containing the properties of the orbit
        and the positions of the satellite.
    """
    if cols is None:
        cols = (1, 2, 0)
    LOGGER.info("loading ephemeris: %s", stream.name)

    comments = []
    data = []
    for item in stream:
        item = item.strip()
        if item.startswith("#"):
            comments.append(item[1:])
        else:
            data.append(list(float(value) for value in item.split()))
    data = np.asarray(data, dtype=np.float64)

    def to_dict(comments) -> Dict[str, float]:
        """Returns a dictionary describing the parameters of the orbit."""
        result = dict()
        for item in comments:
            key, value = item.split("=")
            result[key.strip()] = float(value)
        return result

    return to_dict(comments), tuple(data[:, item] for item in cols)


def interpolate(lon: np.ndarray, lat: np.ndarray, time: np.ndarray) -> Tuple:
    """Interpolate the given orbit at high resolution (0.5 seconds)

    Args:
        lon (numpy.ndarray): Longitudes (in degrees)
        lat (numpy.ndarray): Latitudes (in degrees)
        time (np.ndarray): Date of the positions (in seconds).

    Returns:
        tuple: lon, lat, and time interpolated
    """
    x, y, z = math.spher2cart(lon, lat)
    time_hr = np.arange(time[0], time[-1], 0.5, dtype=time.dtype)

    x = np.interp(time_hr, time, x)
    y = np.interp(time_hr, time, y)
    z = np.interp(time_hr, time, z)

    lon, lat = math.cart2spher(x, y, z)

    return lon, lat, time_hr


def rearrange_orbit(cycle_duration: float, lon: np.ndarray, lat: np.ndarray,
                    time: np.ndarray
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rearrange orbit starting from pass 1

    Detect the beginning of pass 1 in the ephemeris. By definition, it is
    the first passage at southernmost latitude.

    Args:
        cycle_duration (float): Cycle time in seconds.
        lon (numpy.ndarray): Longitudes (in degrees)
        lat (numpy.ndarray): Latitudes (in degrees)
        time (np.ndarray): Date of the positions (in seconds).

    Returns:
        tuple: lon, lat, and time rearranged.
    """
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
        LOGGER.warning('there are negative times in your orbit')
    return lon, lat, time


def calculate_pass_time(lat: np.ndarray, time: np.ndarray) -> np.ndarray:
    """Compute the initial time of each pass

    Args:
        lat (numpy.ndarray): Latitudes (in degrees)
        time (np.ndarray): Date of the latitudes (in seconds).

    Returns:
        numpy.ndarray: Start date of half-orbits.
    """
    dy = np.roll(lat, 1) - lat
    indexes = np.where(((dy < 0) & (np.roll(dy, 1) >= 0))
                       | ((dy > 0)
                          & (np.roll(dy, 1) <= 0)))
    return time[indexes[0]]


def select_box(box: math.Box, lon: np.ndarray, lat: np.ndarray,
               time: np.ndarray
               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Selects the orbit in the defined box.

    Args:
        box (math.Box): Geographical selection
        lon (numpy.ndarray): Longitudes (in degrees)
        lat (numpy.ndarray): Latitudes (in degrees)
        time (np.ndarray): Date of the positions (in seconds).

    Returns:
        tuple: the selected positions and associated dates.
    """
    mask = box.within(lon, lat)
    x_al = math.curvilinear_distance(lon, lat, VOLUMETRIC_MEAN_RADIUS)
    return lon[mask], lat[mask], time[mask], x_al[mask]


class Orbit:
    """Properties of one orbit

    Args:
        heigth (float): Satellite height (in m)
        lat (numpy.ndarray): Latitudes (in degrees)
        lon (numpy.ndarray): Longitudes (in degrees)
        pass_time (np.ndarray): Start date of half-orbits.
        time (np.ndarray): Date of the positions (in seconds).
        x_al (np.ndarray): Along track distance
        curvilinear_distance (float): Distance covered by the satellite during
            a complete cycle.
        shift_time (float): Time shift to be applied.
    """
    def __init__(self, height: float, lat: np.ndarray, lon: np.ndarray,
                 pass_time: np.ndarray, time: np.ndarray, x_al: np.ndarray,
                 curvilinear_distance: float, shift_time: Optional[float]):
        self.height = height
        self.lat = lat
        self.lon = lon
        self.pass_time = pass_time
        self.shift_time = shift_time
        self.time = time
        self.x_al = x_al
        self.curvilinear_distance = curvilinear_distance

    def cycle_duration(self) -> np.timedelta64:
        """Get the cycle duration"""
        return ((self.time[-1]).astype(np.int64) *
                1e6).astype("timedelta64[us]")

    def passes_per_cycle(self):
        """Get the number of passes per cycle"""
        return len(self.pass_time)

    def pass_duration(self, number: int) -> np.timedelta64:
        """Get the duration of a given pass.

        Args:
            number (int): track number (must be in [1, passes_per_cycle()])

        Returns:
            numpy.datetime64: track duration
        """
        if number == self.passes_per_cycle():
            return np.timedelta64(
                int((self.time[-1] - self.pass_time[-1]) * 1e6), 'us')
        return np.timedelta64(
            int((self.pass_time[number] - self.pass_time[number - 1]) * 1e6),
            'us')

    def decode_absolute_pass_number(self, number: int) -> Tuple[int, int]:
        """Calculate the cycle and pass number from a given absolute pass
        number.

        Args:
            number (int): absolute pass number

        Returns:
            tuple: cycle and pass number
        """
        number -= 1
        return (int(number / self.passes_per_cycle()) + 1,
                (number % self.passes_per_cycle()) + 1)

    def encode_absolute_pass_number(self, cycle_number: int,
                                    pass_number: int) -> int:
        """Calculate the absolute pass number for a given half-orbit.

        Args:
            cycle_number (int): Cycle number
            pass_number (int): Pass number

        Returns:
            int: Absolute pass number
        """
        passes_per_cycle = self.passes_per_cycle()
        if not 1 <= pass_number <= passes_per_cycle:
            raise ValueError(f"pass_number must be in [1, {passes_per_cycle}")
        return (cycle_number - 1) * self.passes_per_cycle() + pass_number

    def iterate(self,
                first_date: Optional[np.datetime64] = None,
                last_date: Optional[np.datetime64] = None,
                absolute_pass_number: int = 1
                ) -> Iterator[Tuple[int, int, np.datetime64]]:
        """Obtain all half-orbits within the defined time interval.

        Args:
            first_date (numpy.datetime64): First date of the period to be
                considered.
            last_date (numpy.datetime64): Last date of the period to be
                considered.
            absolute_pass_number (int, optional): Absolute number of the first
                pass to be returned.

        Returns:
            iterator: An iterator for all passes in the interval pointing to
            the cycle number, pass number and start date of the half-orbit.
        """
        date = first_date or np.datetime64(datetime.datetime.now())
        last_date = last_date or first_date + self.cycle_duration()
        while date <= last_date:
            cycle_number, pass_number = self.decode_absolute_pass_number(
                absolute_pass_number)

            yield cycle_number, pass_number, date

            # Shift the date of the duration of the generated pass
            date += self.pass_duration(pass_number)

            # Update of the number of the next pass to be generated
            absolute_pass_number += 1
        return StopIteration


class Pass:
    """Handle one pass

    Args:
        lat_nadir (np.ndarray): Latitudes at nadir
        lat (np.ndarray): Latitudes of the swath
        lon_nadir (np.ndarray): Longitudes at nadir
        lon (np.ndarray): Longitudes of the swath
        x_ac (np.ndarray): Across track distance.
        x_al (np.ndarray): Along track distance.
    """
    def __init__(self, lat_nadir: np.ndarray, lat: np.ndarray,
                 lon_nadir: np.ndarray, lon: np.ndarray, time: np.ndarray,
                 x_ac: np.ndarray, x_al: np.ndarray):
        self.lat_nadir = lat_nadir
        self.lat = lat
        self.lon_nadir = lon_nadir
        self.lon = lon
        self.timedelta = ((time * 1e6).astype(
            np.int64)).astype("timedelta64[us]")
        self._time = None
        self.x_ac = x_ac
        self.x_al = x_al

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, date: np.datetime64) -> None:
        self._time = date + (self.timedelta - self.timedelta[0])


def calculate_orbit(parameters: settings.Parameters,
                    ephemeris: TextIO) -> Orbit:
    """Computes the orbit nadir on a subdomain.

    The path of the satellite is given by the orbit file and the subdomain
    corresponds to the one in the model. Note that a subdomain can be manually
    added in the parameters file.

    Args:
        parameters (settings.Parameters): Simulation parameters.
        ephemeris (typing.TextIO): Stream to the ephemeris CSV file to be
            processed.

    Returns:
        Orbit: The orbit's loaded into memory.
    """
    properties, (lon, lat,
                 time) = load_ephemeris(ephemeris,
                                        cols=parameters.ephemeris_cols)

    if np.any(lon > 360) or np.any(np.abs(lat) > 90):
        raise RuntimeError("The definition of ephemeris is incorrect. Is the "
                           "'ephemeris_cols' parameter well defined?")

    # Put longitude in [-180, 180[
    lon = math.normalize_longitude(lon)

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
        lon = math.normalize_longitude(lon + parameters.shift_lon)

    # Rearrange orbit starting from pass 1
    lon, lat, time = rearrange_orbit(cycle_duration, lon, lat, time)

    # Calculates the along track distance
    distance = math.curvilinear_distance(lon, lat, VOLUMETRIC_MEAN_RADIUS)

    # Interpolate the final orbit according the given along track resolution
    x_al = np.arange(distance[0],
                     distance[-2],
                     parameters.delta_al,
                     dtype=distance.dtype)
    time = np.interp(x_al, distance[:-1], time[:-1])

    x, y, z = math.spher2cart(lon, lat)
    x = np.interp(x_al, distance[:-1], x[:-1])
    y = np.interp(x_al, distance[:-1], y[:-1])
    z = np.interp(x_al, distance[:-1], z[:-1])
    lon, lat = math.cart2spher(x, y, z)

    return Orbit(parameters.height, lat, lon,
                 np.sort(calculate_pass_time(lat, time)), time, x_al,
                 distance[-1], parameters.shift_time)


def calculate_pass(pass_number: int, orbit: Orbit,
                   parameters: settings.Parameters) -> Optional[Pass]:
    """Get the properties of an half-orbit

    Args:
        pass_number (int): Pass number
        orbit (Orbit): Orbit describing the pass to be calculated.
        parameters (settings.Parameters): Simulation parameters.

    Returns:
        Pass, optional: The properties of the pass or None, if the pas is
        missing in the selected area.
    """
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

    # Selects the orbit in the defined box
    mask = parameters.box.within(lon_nadir, lat_nadir)
    if np.all(~mask):
        return
    lon_nadir = lon_nadir[mask]
    lat_nadir = lat_nadir[mask]
    time = time[mask]
    x_al = x_al[mask]

    # Compute accross track distances from nadir
    # Number of points in half of the swath
    half_swath = int((parameters.half_swath - parameters.half_gap) /
                     parameters.delta_ac) + 1
    x_ac = np.arange(half_swath) * parameters.delta_al + parameters.half_gap
    x_ac = np.hstack((-np.flip(x_ac), x_ac))

    location = np.ascontiguousarray(
        np.vstack(math.spher2cart(lon_nadir, lat_nadir)).T)
    satellite_direction = math.satellite_direction(location)
    lon, lat = math.calculate_swath(parameters.delta_ac, parameters.half_gap,
                                    half_swath, VOLUMETRIC_MEAN_RADIUS,
                                    location, satellite_direction)

    return Pass(lat_nadir, lat, lon_nadir, lon, time, x_ac, x_al)
