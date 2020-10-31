# Copyright (c) 2020 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Main program
------------

This module defines the main :func:`function <launch>` handling the simulation
of SWOT products as well as the entry point of the main program.
"""
from typing import Dict, Optional, Tuple
import argparse
import datetime
import logging
import os
import sys
import traceback
import dask.distributed
import dask.bag
import dateutil.parser
import numpy as np
from . import exception
from . import logbook
from . import product_specification
from . import orbit_propagator
from . import settings
from .error import generator

#: Logger of this module
LOGGER = logging.getLogger(__name__)

#: Tasks launched
TASKS = {}


def datetime_type(value: str) -> np.datetime64:
    """The option should define a datetime

    Args:
        value (str): Value to parse

    Returns:
        numpy.datetime64: The date parsed

    Raises:
        ValueError: if the option provided does not set a valid date.
    """
    # Check if value is a datetime
    try:
        result = dateutil.parser.parse(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            f"invalid date time {value!r}: {error!s}")
    return np.datetime64(result)


def writable_directory(value: str) -> str:
    """The option should define a writable directory

    Args:
        value (str): Value to parse

    Returns:
        str: The path to the directory

    Raises:
        argparse.ArgumentTypeError: If the option provided does not define a
        writable directory.
    """
    # TODO: create directory if it does not exist?
    # Check if value exists
    if not os.path.exists(value):
        _err = f"no such file or directory: {value!r}"
        raise argparse.ArgumentTypeError(_err)
    # Check if value is a directory
    if not os.path.isdir(value):
        raise argparse.ArgumentTypeError(f"{value!r} is not a directory")
    # Check if directory value is writable
    if not os.access(value, os.W_OK):
        raise argparse.ArgumentTypeError(
            f"cannot open directory {value!r} for writing")
    return value


def usage() -> argparse.Namespace:
    """Parse the options provided on the command line.

    Returns:
        argparse.Namespace: The parameters provided on the command line.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("settings",
                        type=argparse.FileType('r'),
                        help="Path to the parameters file")
    group = parser.add_argument_group("General", "Simulation general settings")
    group.add_argument("--first-date",
                       help="The first date to be processed. "
                       "Default to the current date",
                       type=datetime_type,
                       default=np.datetime64(datetime.date.today()))
    group.add_argument("--last-date",
                       help="The last date to be processed. "
                       "Default to the last date allowing to cover an entire "
                       "cycle.",
                       type=datetime_type)
    group = parser.add_argument_group("Execution",
                                      "Runtime parameters options ")
    group.add_argument("--debug",
                       action="store_true",
                       help="Put swot simulator in debug mode")
    group.add_argument("--log",
                       metavar='PATH',
                       help="Path to the logbook to use",
                       type=argparse.FileType("w"))
    group.add_argument("--scheduler-file",
                       help="Path to a file with scheduler information to "
                       "launch swot simulator on a cluster. By "
                       "default, use a local cluster.",
                       metavar='PATH',
                       type=argparse.FileType("r"))
    group = parser.add_argument_group("LocalCluster",
                                      "Dask local cluster option")
    group.add_argument("--n-workers",
                       help="Number of workers to start (Default to 1)",
                       type=int,
                       metavar='N',
                       default=1)
    group.add_argument("--processes",
                       help="Whether to use processes (True) or threads "
                       "(False).  Defaults to False",
                       action="store_true")
    group.add_argument("--threads-per-worker",
                       help="Number of threads per each worker. "
                       "(Default to 1)",
                       type=int,
                       metavar='N',
                       default=1)
    namespace = argparse.Namespace()
    namespace, _ = parser._parse_known_args(sys.argv[1:], namespace)
    if "scheduler_file" in namespace:
        for item in ["n_workers", "processes", "threads_per_worker"]:
            if item in namespace:
                item = item.replace("_", "-")
                raise RuntimeError(
                    f"--{item}: not allowed with argument --scheduler-file")
    return parser.parse_args()


def file_path(first_date: np.datetime64,
              last_date: np.datetime64,
              cycle_number: int,
              pass_number: int,
              parameters: settings.Parameters,
              nadir: bool = False) -> str:
    """Get the absolute path of the file to be created.

    Args:
        first_date (numpy.datetime64): Date of the first simulated measurement.
        last_date (numpy.datetime64): Date of the last simulated measurement.
        cycle_number (int): Cycle number of the file to be created.
        pass_number (int): Pass number of the file to be created.
        parameters (settings.Parameters): Simulation parameters.
        nadir (bool, optional): True if the product to be created contains the
            nadir measurements, false if it is a product containing the swath
            measurements.

    Returns:
        str: The path to the file to be created.
    """
    first_date = first_date.astype(datetime.datetime)
    last_date = last_date.astype(datetime.datetime)
    product_type = "nadir" if nadir else "karin"
    dirname = os.path.join(parameters.working_directory, product_type,
                           first_date.strftime("%Y"))
    os.makedirs(dirname, exist_ok=True)
    if nadir:
        filename = (f"SWOT_GPN_2P1P{cycle_number:03d}_{pass_number:03d}_"
                    f"{first_date:%Y%m%d}_{first_date:%H%M%S}_"
                    f"{last_date:%Y%m%d}_{last_date:%H%M%S}.nc")
    else:
        product_type = "".join(
            [item.capitalize() for item in parameters.product_type.split("_")])
        filename = (f"SWOT_L2_LR_SSH_{product_type}_"
                    f"{cycle_number:03d}_{pass_number:03d}_"
                    f"{first_date:%Y%m%dT%H%M%S}_{last_date:%Y%m%dT%H%M%S}_"
                    "DG10_01.nc")
    return os.path.join(dirname, filename)


def sum_error(errors: Dict[str, np.ndarray], swath: bool = True) -> np.ndarray:
    """Calculate the sum of errors

    Args:
        errors (dict): Simulated errors.
        swath (bool, optional): True if the measurements processed concern the
            swath.

    Returns:
        numpy.ndarray: The sum of errors
    """
    dims = 2 if swath else 1
    return np.add.reduce(
        [item for item in errors.values() if len(item.shape) == dims])


def simulate(args: Tuple[int, int, np.datetime64],
             error_generator: generator.Generator,
             orbit: orbit_propagator.Orbit,
             parameters: settings.Parameters,
             logging_server: Optional[Tuple[str, int, int]] = None) -> None:
    """Simulate a pass.

    Args:
        simulation arguments (tuple): cycle & pass number, date of the first
            half-orbit measurement
        error_generator (generator.Generator): Measurement error generator.
        orbit (orbit_propagator.Orbit): Orbit propagator.
        parameters (settings.Parameters): Simulation parameters.
        logging_server (tuple, optional): Log server connection settings.
    """
    cycle_number, pass_number, date = args

    # Initialize this worker's logger.
    if logging_server is not None:
        logbook.setup_worker_logging(logging_server)

    # Calculation of the end date of the track.
    last_date = date + orbit.pass_shift(pass_number)

    # Paths of products to be generated.
    swath_path = None
    nadir_path = None

    # Generation of the names of the files to be created.
    if parameters.swath:
        swath_path = file_path(date, last_date, cycle_number, pass_number,
                               parameters)

        # If the product has already been produced, the generation of this
        # half orbit is disabled.
        if os.path.exists(swath_path):
            swath_path = None

    if parameters.nadir:
        nadir_path = file_path(date,
                               last_date,
                               cycle_number,
                               pass_number,
                               parameters,
                               nadir=True)
        if os.path.exists(nadir_path):
            nadir_path = None

    # To continue, there must be at least one task left for this pass.
    if swath_path is None and nadir_path is None:
        return

    # Compute the spatial/temporal position of the satellite
    track = orbit_propagator.calculate_pass(pass_number, orbit, parameters)
    if track is None:
        return

    # Set the simulated date
    track.time = date

    # Mask to set the measurements outside the requirements of the mission to
    # NaN.
    mask = track.mask()
    # Interpolation of the SWH if the user wishes.
    if parameters.swh_plugin is not None:
        # The nadir and swath data are concatenated to process the
        # interpolation of the SSH in one time (swath and nadir).
        lon = np.c_[track.lon, track.lon_nadir[:, np.newaxis]]
        lat = np.c_[track.lat, track.lat_nadir[:, np.newaxis]]
        swath_time = np.repeat(track.time, lon.shape[1]).reshape(lon.shape)
        swh = parameters.swh_plugin.interpolate(lon.flatten(), lat.flatten(),
                                                swath_time.flatten())
        swh_all = swh.reshape(lon.shape)
        swh = +swh_all[:, :-1]
    else:
        swh_all = None
        swh = np.array([
            parameters.swh,
        ] * len(track.x_ac))

    # Calculation of instrumental errors
    noise_errors = error_generator.generate(cycle_number,
                                            orbit.curvilinear_distance,
                                            track.time, track.x_al, track.x_ac,
                                            swh)
    for error in noise_errors.values():
        # Only the swaths must be masked
        if len(error.shape) == 2:
            error *= mask

    # Interpolation of the SSH if the user wishes.
    if parameters.ssh_plugin is not None:
        # The nadir and swath data are concatenated to process the
        # interpolation of the SSH in one time (swath and nadir).
        lon = np.c_[track.lon, track.lon_nadir[:, np.newaxis]]
        lat = np.c_[track.lat, track.lat_nadir[:, np.newaxis]]
        swath_time = np.repeat(track.time, lon.shape[1]).reshape(lon.shape)
        ssh = parameters.ssh_plugin.interpolate(lon.flatten(), lat.flatten(),
                                                swath_time.flatten())
        ssh = ssh.reshape(lon.shape)
    else:
        ssh = None

    if swath_path:
        LOGGER.info("generate swath %d/%d [%s, %s]", cycle_number, pass_number,
                    track.time[0], track.time[-1])

        # Create the swath dataset
        product = product_specification.Swath(track, parameters.central_pixel,
                                              parameters.product_type)

        if ssh is not None:
            product.ssh((ssh[:, :-1] * mask) + sum_error(noise_errors))
            if noise_errors:
                product.simulated_true_ssh(ssh[:, :-1])
        if swh_all is not None:
            product.swh((swh * mask))

        product.update_noise_errors(noise_errors)
        product.to_netcdf(cycle_number, pass_number, swath_path,
                          parameters.complete_product)

    # Create the nadir dataset
    if nadir_path:
        LOGGER.info("generate nadir %d/%d [%s, %s]", cycle_number, pass_number,
                    track.time[0], track.time[-1])

        product = product_specification.Nadir(track,
                                              standalone=not parameters.swath)

        # Interpolation of the SSH if the user wishes.
        if ssh is not None:
            product.ssh(ssh[:, -1] + sum_error(noise_errors, swath=False))
            if noise_errors:
                product.simulated_true_ssh(ssh[:, -1])
        if swh_all is not None:
            product.swh(swh_all[:, -1])
        product.update_noise_errors(noise_errors)
        product.to_netcdf(cycle_number, pass_number, nadir_path,
                          parameters.complete_product)


def launch(client: dask.distributed.Client,
           parameters: settings.Parameters,
           logging_server: Optional[Tuple[str, int, int]],
           first_date: Optional[np.datetime64] = None,
           last_date: Optional[np.datetime64] = None):
    """Executes the simulation set to the selected period.

    Args:
        client (dask.distributed.Client): Client connected to the Dask
            cluster.
        parameters (settings.Parameters): Simulation parameters.
        logging_server (tuple, optional): Log server connection settings.
        first_date (numpy.datetime64): First date of the simulation.
        last_date (numpy.datetime64): Last date of the simulation.
    """
    # Displaying Dask client information.
    LOGGER.info(client)

    assert parameters.ephemeris is not None

    # Calculation of the properties of the orbit to be processed.
    with open(parameters.ephemeris, "r") as stream:  # type: TextIO
        orbit = orbit_propagator.calculate_orbit(parameters,
                                                 stream)  # type: ignore

    # Initialization of measurement error generators
    error_generator = generator.Generator(parameters, first_date)

    # Scatter data into distributed memory
    _error_generator = client.scatter(error_generator)
    _parameters = client.scatter(parameters)
    _orbit = client.scatter(orbit)

    bag = dask.bag.from_sequence(orbit.iterate(first_date, last_date),
                                 npartitions=len(
                                     client.scheduler_info()['workers']))
    bag.map(simulate,
            error_generator=_error_generator,
            orbit=_orbit,
            parameters=_parameters,
            logging_server=logging_server).compute()


def main():
    """main function"""
    args = usage()

    # Setup log
    logger, logging_server = logbook.setup(args.log, args.debug)

    client = dask.distributed.Client(
        dask.distributed.LocalCluster(
            protocol="tcp",
            n_workers=args.n_workers,
            processes=args.processes,
            threads_per_worker=args.threads_per_worker)
    ) if args.scheduler_file is None else dask.distributed.Client(
        scheduler_file=args.scheduler_file.name)

    client.wait_for_workers(1)

    try:
        parameters = settings.eval_config_file(args.settings.name)
        launch(client, settings.Parameters(parameters), logging_server,
               args.first_date, args.last_date)

        client.close()
        logger.info("End of processing.")
        return 0
    #: pylint: disable=broad-except
    # All exceptions are caught in the main function to display it in the log
    # and return the appropriate error code.
    except Exception as exc:
        # Clients are stopped before writing the error message to the log to
        # ensure that the log will end with the exception captured.
        client.close()
        logger.error(
            exception.structured_traceback(
                exc, traceback.extract_tb(sys.exc_info()[2])))
        logger.error("End of processing.")
    #: pylint: enable=broad-except
    return 1


if __name__ == '__main__':
    sys.exit(main())
