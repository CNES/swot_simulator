# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Interpolation of the SSH on the NZ grid.
========================================
"""
from typing import Optional
import pathlib

#
import numpy as np
import pyinterp
import pyinterp.backends.xarray
import xarray as xr

#
from .. import Interface


class Base(Interface):
    """Base class for NZ interface.

    Args:
        path: Path to the time series.
        ssh: Name of the variable containing the SSH.
    """
    def __init__(self, path: str, ssh: Optional[str] = None):
        self.ssh = ssh or "ssh"
        # Delta time between two time steps.
        self._dt = np.timedelta64(3, "h")
        # Load time series.
        self._ts = self._load_time_series(path)
        # Have we found the right maps to process?
        if len(self._ts) == 0:
            raise ValueError(f"No files found in {path}")

    def _load_time_series(self, path: str) -> np.ndarray:
        """Load time series from a directory"""
        items = []
        previous = None
        for item in sorted(pathlib.Path(path).glob("*.nc")):
            with xr.open_dataset(item) as ds:
                current = ds.ocean_time.values[0].astype("datetime64[M]")
                if (previous is not None
                        and (current - previous != np.timedelta64(1, "M"))):
                    raise ValueError(f"Time series not continuous")
                items.append((current, str(item)))
                previous = current
        length = max(len(item[1]) for item in items)
        return np.array(
            items,
            dtype={
                "names": ("date", "path"),
                "formats": ("datetime64[M]", f"U{length}"),
            },
        )

    @staticmethod
    def _floor_to_dt(value: np.datetime64) -> np.datetime64:
        """Floor a datetime64 to the nearest dt"""
        integral = int(value.astype("<M8[h]").astype("int64") /
                       3)  # type: ignore
        return np.datetime64(integral * 3, "h")

    def _select_ds(self, first_date: np.datetime64,
                   last_date: np.datetime64) -> xr.Dataset:
        """Select the time series to process"""
        first_ts = self._floor_to_dt(first_date)
        last_ts = self._floor_to_dt(last_date) + self._dt
        first_month = first_ts.astype("M8[M]")
        last_month = last_ts.astype("M8[M]")
        ts = self._ts["date"]
        if first_month < ts[0] or last_month > ts[-1]:
            upper_limit = (ts[-1] + np.timedelta64(1, 'M')).astype("M8[s]")
            raise IndexError(
                f"period [{first_date}, {last_date}] is out of range: "
                f"[{ts[0].astype('M8[s]')}, {upper_limit}[")
        mask = (ts >= first_month) & (ts <= last_month)

        paths = self._ts["path"][mask]
        ds = xr.open_dataset(paths[0]).isel(ocean_time=slice(0, -1, None))
        z0 = ds.sel(ocean_time=slice(first_ts, last_ts))

        if len(paths) > 1:
            ds = xr.open_dataset(paths[1]).isel(ocean_time=slice(0, -1, None))
            z1 = ds.sel(ocean_time=slice(first_ts, last_ts))
            return xr.concat([z0, z1], dim="ocean_time")
        return z0


class NZCartesian(Base):
    """Handle the interpolation on NZ cartesian grids

    Args:
        path: Path to the time series.
    """
    def __init__(self, path: str):
        super().__init__(path, ssh="ssh")

    def interpolate(self, lon: np.ndarray, lat: np.ndarray,
                    dates: np.ndarray) -> np.ndarray:
        """Interpolate the SSH to the required coordinates"""
        ds = self._select_ds(
            dates.min(),  # type: ignore
            dates.max())  # type: ignore
        assert np.all(np.diff(ds.ocean_time.values) == self._dt)
        interpolator = pyinterp.backends.xarray.RegularGridInterpolator(
            ds[self.ssh])
        return interpolator(dict(lat_rho=lat.ravel(),
                                 lon_rho=lon.ravel(),
                                 ocean_time=dates.ravel()),
                            method="bilinear",
                            bounds_error=False).reshape(lon.shape)


class NZMesh(Base):
    """Handle the interpolation on NZ mesh grids
    """
    def __init__(self, path: str):
        super().__init__(path, ssh="zeta")

    def interpolate(self, lon: np.ndarray, lat: np.ndarray,
                    dates: np.ndarray) -> np.ndarray:
        """Interpolate the SSH to the required coordinates"""
        ds = self._select_ds(
            dates.min(),  # type: ignore
            dates.max())  # type: ignore

        assert np.all(np.diff(ds.ocean_time.values) == self._dt)
        assert np.all(np.diff(ds.lon_rho.values, axis=0) < 1e-10)
        assert np.all(np.diff(ds.lat_rho.values, axis=1) < 1e-10)

        t_axis = pyinterp.TemporalAxis(ds.ocean_time.values)

        grid3d = pyinterp.Grid3D(
            pyinterp.Axis(ds.lon_rho.values[0, :], is_circle=True),
            pyinterp.Axis(ds.lat_rho.values[:, 0]), t_axis,
            ds[self.ssh].values.T)

        ssh = pyinterp.trivariate(grid3d,
                                  lon.ravel(),
                                  lat.ravel(),
                                  t_axis.safe_cast(dates.ravel()),
                                  num_threads=1).reshape(lon.shape)
        return ssh
