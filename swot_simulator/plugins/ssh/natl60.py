# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Interpolate SSH from NATL60 model
=================================
"""
import logging
import time

import dask.array as da
import numpy as np
import pyinterp
import pyinterp.geodetic
import xarray as xr

from .. import data_handler

LOGGER = logging.getLogger(__name__)


class NATL60(data_handler.IrregularGridHandler):
    def __init__(self, path: str):
        loader = NATL60.ZarrLoader(path)
        super().__init__(loader)

    class ZarrLoader(data_handler.DatasetLoader):
        def __init__(self, path: str):
            with xr.open_zarr(path) as ds:
                lon = ds.nav_lon
                lat = ds.nav_lat
                ssh = ds.sossheig
                time = ds.time_counter.data.astype("datetime64[us]")

            self.dataset = xr.Dataset(
                dict(ssh=ssh),
                coords=dict(lon=lon, lat=lat, time=time)
            )
            self.time_delta = self._calculate_time_delta(self.dataset.time)

        def load_dataset(self, first_date: np.datetime64,
                         last_date: np.datetime64):
            first_date = self._shift_date(first_date.astype("datetime64[ns]"),
                                          -1, self.time_delta)
            last_date = self._shift_date(last_date.astype("datetime64[ns]"), 1,
                                         self.time_delta)

            if first_date < self.dataset.time[
                    0] or last_date > self.dataset.time[-1]:
                raise IndexError(
                    f"period [{first_date}, {last_date}] is out of range: "
                    f"[{self.dataset.time[0]}, {self.dataset.time[-1]}]")

            # Mask for selecting data covering the time period provided.
            mask = (self.dataset.time.data >=
                    first_date) & (self.dataset.time.data <= last_date)
            return self.dataset.isel(time=np.argwhere(mask).squeeze())

    @staticmethod
    def _spatial_interp(
        z_model: da.Array,
        x_model: da.Array,
        y_model: da.Array,
        x_sat: np.ndarray,
        y_sat: np.ndarray,
    ):
        mesh = pyinterp.RTree(dtype=np.dtype("float32"))

        start_time = time.time()

        ssh = z_model.compute()
        defined = ~np.isnan(ssh)
        ssh = ssh[defined]

        lon = x_model[defined].compute()
        lat = y_model[defined].compute()

        # The tree is built and the interpolation is calculated
        coordinates = np.vstack((lon, lat)).T
        del lon, lat

        LOGGER.debug(
            "loaded %d MB in %.2fs",
            (coordinates.nbytes + ssh.nbytes) // 1024**2,
            time.time() - start_time,
        )

        start_time = time.time()
        mesh.packing(coordinates, ssh)
        del coordinates, ssh
        LOGGER.debug("mesh build in %.2fs", time.time() - start_time)

        start_time = time.time()
        z_sat, _ = mesh.radial_basis_function(
            np.vstack((x_sat, y_sat)).T.astype("float32"),
            within=True,
            k=11,
            radius=8000,
            rbf="thin_plate",
            num_threads=1,
        )
        LOGGER.debug("interpolation done in %.2fs", time.time() - start_time)
        del mesh
        return z_sat.astype("float32")
