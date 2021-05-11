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
import numpy as np
import pyinterp
import pyinterp.geodetic
import xarray as xr

from .. import data_handler

LOGGER = logging.getLogger(__name__)


class MITGCM(data_handler.IrregularGridHandler):
    def __init__(self, grid_path: str, eta_path: str):
        loader = MITGCM.ZarrLoader(grid_path, eta_path)
        super().__init__(loader)

    class ZarrLoader(data_handler.DatasetLoader):
        def __init__(self, grid_path: str, eta_path: str):
            dataset = xr.merge(
                [xr.open_zarr(grid_path),
                 xr.open_zarr(eta_path)]).rename({
                     "XC": "lon",
                     "YC": "lat",
                     "Eta": "ssh"
                 })
            self.dataset = dataset[["ssh"]]
            self.dataset.dtime.load()
            self.time_delta = self._calculate_time_delta(self.dataset.dtime)

        def load_dataset(self, first_date: np.datetime64,
                         last_date: np.datetime64):
            first_date = self._shift_date(first_date.astype("datetime64[ns]"),
                                          -1, self.time_delta)
            last_date = self._shift_date(last_date.astype("datetime64[ns]"), 1,
                                         self.time_delta)

            if first_date < self.dataset.dtime[
                    0] or last_date > self.dataset.dtime[-1]:
                raise IndexError(
                    f"period [{first_date}, {last_date}] is out of range: "
                    f"[{self.dataset.dtime[0]}, {self.dataset.dtime[-1]}]")

            # Mask for selecting data covering the time period provided.
            mask = (self.dataset.dtime.data >=
                    first_date) & (self.dataset.dtime.data <= last_date)
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
        LOGGER.debug(
            "loaded %d MB in %.2fs",
            (coordinates.nbytes + z.nbytes) // 1024**2,
            time.time() - start_time,
        )
        start_time = time.time()
        mesh.packing(coordinates, z)
        LOGGER.debug("mesh build in %.2fs", time.time() - start_time)

        del coordinates, z

        start_time = time.time()
        z, _ = mesh.radial_basis_function(
            np.vstack((x_sat, y_sat)).T.astype("float32"),
            within=True,
            k=11,
            radius=55000,
            rbf="thin_plate",
            num_threads=1,
        )
        LOGGER.debug("interpolation done in %.2fs", time.time() - start_time)
        del mesh
        return z.astype("float32")
