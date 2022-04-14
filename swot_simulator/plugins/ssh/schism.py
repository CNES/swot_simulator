# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Interpolation of the SSH NEMO
==============================
"""
import pathlib
import re

import numpy as np
import pyinterp
import xarray as xr

from .. import Interface, data_handler


class SCHISM(Interface):
    """Interpolation of the SSH NEMO."""

    def __init__(self, path: str):
        # Delta time between two time steps.
        self._dt = np.timedelta64(3, "h")
        # Load time series.
        self._ts = self._load_time_series(path)
        # Have we found the right maps to process?
        if len(self._ts) == 0:
            raise ValueError(f"No files found in {path}")

    def _load_time_series(self, path: str) -> np.ndarray:
        """Load the time series of the SSH."""
        items = []
        pattern = re.compile(r"schism_(?P<date>\d{8}T\d{6})_\d{8}T\d{6}.nc")
        for item in pathlib.Path(path).iterdir():
            match = pattern.match(item.name)
            if match:
                date_str = match.group("date")
                date = np.datetime64(
                    date_str[:4] + "-" + date_str[4:6] + "-" + date_str[6:8] +
                    "T" + date_str[9:11] + ":" + date_str[11:13] + ":" +
                    date_str[13:15], "s")
                for ix in range(8):
                    items.append((ix, date, str(item)))
                    date += self._dt
        length = max(len(item[2]) for item in items)
        return np.array(
            items,
            dtype={
                "names": ("index", "date", "path"),
                "formats": ("uint64", "datetime64[s]", f"U{length}"),
            },
        )

    def _select_ds(self, first_date: np.datetime64,
                   last_date: np.datetime64) -> np.ndarray:
        """Select the time series to process."""
        first_date -= self._dt
        last_date += self._dt
        ts = self._ts["date"]
        if first_date < ts[0] or last_date > ts[-1]:
            raise IndexError(
                f"period [{first_date}, {last_date}] is out of range: "
                f"[{ts[0]}, {ts[-1]}]")
        mask = (ts >= first_date) & (ts <= last_date)
        return self._ts[mask]

    @staticmethod
    def _rtree(ds: xr.Dataset, index: np.uint64) -> pyinterp.RTree:
        """Constructs the interpolator."""
        ssh = ds.variables["elev"].isel(dict(time=index)).values

        # Filter out NaN values.
        mask = ~np.isnan(ssh)
        lon = ds.variables["SCHISM_hgrid_node_x"].values[mask]
        lat = ds.variables["SCHISM_hgrid_node_y"].values[mask]
        ssh = ssh[mask]

        coordinates = np.vstack((lon, lat)).T

        mesh = pyinterp.RTree(dtype=np.dtype("float32"))
        mesh.packing(coordinates, ssh)
        return mesh

    def interpolate(self, lon: np.ndarray, lat: np.ndarray,
                    dates: np.ndarray) -> np.ndarray:
        """Interpolate the SSH to the required coordinates."""
        selected = self._select_ds(
            dates.min(),  # type: ignore
            dates.max())  # type: ignore

        layers = []
        xp = []

        for ix in range(len(selected)):
            index = selected[ix]["index"]
            ds = xr.open_dataset(selected[ix]["path"])
            # Checks that's the time delta between two time steps is the
            # expected one.
            assert np.all(np.diff(ds.time.values) == self._dt)
            xp.append(np.asarray(ds.time[index]).astype("datetime64[us]"))
            mesh = self._rtree(ds, index)
            ssh, _ = mesh.radial_basis_function(
                np.vstack((lon.ravel(), lat.ravel())).T.astype("float32"),
                within=True,
                k=11,
                radius=10000,  # Spatial resolution in meters.
                rbf="thin_plate",
                num_threads=1,
            )
            layers.append(ssh)

        layers = np.stack(layers)
        xp = np.array(xp)

        # Time interpolation of the SSH.
        return data_handler.time_interp(
            xp.astype("int64"),
            layers,
            dates.astype("datetime64[us]").astype("int64"),
        ).reshape(lon.shape)
