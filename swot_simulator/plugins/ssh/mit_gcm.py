import logging
import dask.distributed
import dask.array as da
import numba as nb
import numpy as np
import pyinterp
import xarray as xr
from . import detail

LOGGER = logging.getLogger(__name__)


@nb.njit('(float32[::1])(int64[::1], float32[:, ::1], int64[::1])', cache=True)
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

        # The undefined values are filtered
        z_face = z_model[face, :].compute()
        defined = ~np.isnan(z_face)
        x += (x_face[defined].flatten(), )
        y += (y_face[defined].flatten(), )
        z += (z_face[defined].flatten(), )

    # The tree is built and the interpolation is calculated
    mesh.packing(
        np.vstack((np.concatenate(x), np.concatenate(y))).T, np.concatenate(z))
    del x, y, z
    z, _ = mesh.inverse_distance_weighting(np.vstack((x_sat, y_sat)).T,
                                           within=True,
                                           k=11,
                                           radius=55000,
                                           num_threads=1)
    return z


class MITGCM(detail.Interface):
    def __init__(self, xc: xr.DataArray, yc: xr.DataArray, eta: xr.DataArray):
        self.lon = xc.data
        self.lat = yc.data
        self.ssh = eta.data
        self.ts = eta.time.data.astype("datetime64[us]")
        self.dt = self._calculate_dt(self.ts)

    @staticmethod
    def _calculate_dt(time: xr.DataArray):
        """Calculation of the delta T between two consecutive grids"""
        frequency = np.diff(time)
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
                    time: np.ndarray) -> np.ndarray:
        """Interpolate the SSH for the given coordinates"""
        first_date = self._grid_date(time[0], -1)
        last_date = self._grid_date(time[-1], 1)

        if first_date < self.ts[0] or last_date > self.ts[-1]:
            raise IndexError(
                f"period [{first_date}, {last_date}] is out of range: "
                f"[{self.ts[0]}, {self.ts[-1]}]")

        # Mask for selecting data covering the time period provided.
        mask = (self.ts >= first_date) & (self.ts <= last_date)

        LOGGER.debug(f"fetch data for {first_date}, {last_date}")

        # 4D cube representing the data necessary for interpolation.
        frame = self.ssh[mask]

        # Spatial interpolation of the SSH on the different selected grids.
        with dask.distributed.worker_client() as client:
            x_sat = client.scatter(lon)
            y_sat = client.scatter(lat)
            x_model = client.scatter(self.lon)
            y_model = client.scatter(self.lat)

            futures = []
            for index in range(len(frame)):
                z_model = client.scatter(frame[index, :])
                futures.append(
                    client.submit(_spatial_interp, z_model, x_model, y_model,
                                  x_sat, y_sat))

            spatial_interp = client.gather(futures)
            client.cancel([x_sat, y_sat, x_model, y_model])

        # Time interpolation of the SSH.
        return _time_interp(self.ts[mask].astype("int64"),
                            np.stack(spatial_interp),
                            time.astype("datetime64[us]").astype("int64"))
