# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Interpolate SSH from MIT/GCM model
==================================
"""
import logging
import time
import dask.array as da
import numba as nb
import numpy as np
import pyinterp
import xarray as xr
from .. import Interface

LOGGER = logging.getLogger(__name__)


@nb.njit('(float32[::1])(int64[::1], float32[:, ::1], int64[::1])',
         cache=True,
         nogil=True)
def _time_interp(xp: np.ndarray, yp: np.ndarray, xi: np.ndarray) -> np.ndarray:
    """Time interpolation for the different spatial grids interpolated on the
    SWOT data"""
    xp_diff = np.diff(xp)

    assert xp.shape[0] == yp.shape[0] and yp.shape[1] == xi.shape[0]
    assert np.all(xp_diff == xp_diff[0])

    result = np.empty(yp.shape[1:], dtype=yp.dtype)

    step = 1.0 / xp_diff[0]
    size = xp.size

    for ix in range(yp.shape[0]):
        index = int(np.around((xi[ix] - xp[0]) * step))
        assert index >= 0 or index <= size
        if index == size - 1:
            i0 = index - 1
            i1 = index
        else:
            i0 = index
            i1 = index + 1
        t0 = xp[i0]
        t1 = xp[i1]
        dt = t1 - t0
        w0 = (t1 - xi[ix]) / dt
        w1 = (xi[ix] - t0) / dt

        for jx in range(yp.shape[1]):
            result[jx] = (w0 * yp[i0, jx] + w1 * yp[i1, jx]) / (w0 + w1)

    return result


def _spatial_interp(z_model: da.array, x_model: da.array, y_model: da.array,
                    x_sat: np.ndarray, y_sat: np.ndarray):
    mesh = pyinterp.RTree(dtype="float32")
    x, y, z = (), (), ()

    start_time = time.time()

    for face in range(13):
        x_face = x_model[face, :].compute()
        y_face = y_model[face, :].compute()

        # We test if the face covers the satellite positions.
        ix0, ix1 = x_face.min(), x_face.max()
        iy0, iy1 = y_face.min(), y_face.max()

        box = pyinterp.geodetic.Box2D(pyinterp.geodetic.Point2D(ix0, iy0),
                                      pyinterp.geodetic.Point2D(ix1, iy1))
        mask = box.covered_by(x_sat, y_sat)
        if not np.any(mask == 1):
            continue
        del box, mask

        # The undefined values are filtered
        z_face = z_model[face, :].compute()
        defined = ~np.isnan(z_face)
        x += (x_face[defined].flatten(), )
        y += (y_face[defined].flatten(), )
        z += (z_face[defined].flatten(), )

    # The tree is built and the interpolation is calculated
    x = np.concatenate(x)
    y = np.concatenate(y)
    coordinates = np.vstack((x, y)).T
    del x, y

    z = np.concatenate(z)
    LOGGER.debug("loaded %d MB in %.2fs",
                 (coordinates.nbytes + z.nbytes) // 1024**2,
                 time.time() - start_time)
    start_time = time.time()
    mesh.packing(coordinates, z)
    LOGGER.debug("mesh build in %.2fs", time.time() - start_time)

    del coordinates, z

    start_time = time.time()
    z, _ = mesh.radial_basis_function(np.vstack(
        (x_sat, y_sat)).T.astype("float32"),
                                      within=True,
                                      k=11,
                                      radius=55000,
                                      rbf="thin_plate",
                                      num_threads=1)
    LOGGER.debug("interpolation done in %.2fs", time.time() - start_time)
    del mesh
    return z.astype("float32")


class MITGCM(Interface):
    def __init__(self, xc: xr.DataArray, yc: xr.DataArray, eta: xr.DataArray):
        self.lon = xc.data
        self.lat = yc.data
        self.ssh = eta.data
        self.ts = eta.time.data.astype("datetime64[us]")
        self.dt = self._calculate_dt(self.ts)

    @staticmethod
    def _calculate_dt(dates: xr.DataArray):
        """Calculation of the delta T between two consecutive grids"""
        frequency = np.diff(dates)
        if not np.all(frequency == frequency[0]):
            raise RuntimeError(
                "Time series does not have a constant step between two "
                f"grids: {set(frequency)} seconds")
        return frequency[0]

    def _grid_date(self, date: np.datetime64, shift: int):
        """Calculates the grid date immediately before or after the date
        provided"""
        if date.astype("int64") % self.dt.astype("int64") != 0:
            return date + self.dt * shift
        return date

    def interpolate(self, lon: np.ndarray, lat: np.ndarray,
                    dates: np.ndarray) -> np.ndarray:
        """Interpolate the SSH for the given coordinates"""
        first_date = self._grid_date(dates[0], -1)
        last_date = self._grid_date(dates[-1], 1)

        if first_date < self.ts[0] or last_date > self.ts[-1]:
            raise IndexError(
                f"period [{first_date}, {last_date}] is out of range: "
                f"[{self.ts[0]}, {self.ts[-1]}]")

        # Mask for selecting data covering the time period provided.
        mask = (self.ts >= first_date) & (self.ts <= last_date)

        LOGGER.debug("fetch data for %s, %s", first_date, last_date)

        # 4D cube representing the data necessary for interpolation.
        frame = self.ssh[mask]

        # Spatial interpolation of the SSH on the different selected grids.
        start_time = time.time()
        layers = []
        for index in range(len(frame)):
            layers.append(
                _spatial_interp(frame[index, :], self.lon, self.lat, lon, lat))

        # Time interpolation of the SSH.
        layers = np.stack(layers)
        LOGGER.debug("interpolation completed in %.2fs for period %s, %s",
                     time.time() - start_time, first_date, last_date)
        return _time_interp(self.ts[mask].astype("int64"), layers,
                            dates.astype("datetime64[us]").astype("int64"))
