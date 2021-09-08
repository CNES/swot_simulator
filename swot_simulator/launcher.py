# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Main program
------------

This module defines the main :func:`function <launch>` handling the simulation
of SWOT products as well as the entry point of the main program.
"""
from typing import Dict, List, Optional, Tuple
import argparse
import datetime
import importlib
import logging
import os
import pathlib
import platform
import sys
import traceback

import dask.distributed
import dateutil.parser
import numpy as np

from . import (dispatch, exception, logbook, netcdf, orbit_propagator,
               settings, __date__, __version__)
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


def usage() -> argparse.Namespace:
    """Parse the options provided on the command line.

    Returns:
        argparse.Namespace: The parameters provided on the command line.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-h',
                        '--help',
                        action='store_true',
                        help='show this help message and exit')
    group = parser.add_argument_group("General", "Simulation general settings")
    group.add_argument("--first-date",
                       help="The first date to be processed. "
                       "Default to the current date",
                       type=datetime_type,
                       default=np.datetime64("now"))
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
    group = parser.add_argument_group("Configuration")
    group.add_argument("--template",
                       help="Writes the default configuration of the "
                       "simulator into the file and ends the program.",
                       metavar="PATH",
                       type=argparse.FileType("w"))
    namespace = argparse.Namespace()
    namespace, _ = parser._parse_known_args(sys.argv[1:], namespace)

    def add_settings(parser):
        """Added the argument defining the settings of the simulator."""
        parser.add_argument("settings",
                            type=argparse.FileType('r'),
                            help="Path to the parameters file")

    # Displays help and ends the program.
    if "help" in namespace:
        add_settings(parser)
        parser.print_help()
        parser.exit(0)

    # Checking exclusive options.
    if "scheduler_file" in namespace:
        for item in ["n_workers", "processes", "threads_per_worker"]:
            if item in namespace:
                item = item.replace("_", "-")
                raise RuntimeError(
                    f"--{item}: not allowed with argument --scheduler-file")

    # Write the template configurayion file and ends the programm
    if "template" in namespace:
        namespace.template.write(settings.template())
        sys.stdout.write(f"""
The template has been written in the file: {namespace.template.name!r}.
""")
        parser.exit(0)

    # The partial analysis of the command line arguments is finished, the last
    # argument is added and parsed one last time.
    add_settings(parser)

    return parser.parse_args()


def file_path(first_date: np.datetime64,
              last_date: np.datetime64,
              cycle_number: int,
              pass_number: int,
              parameters: settings.Parameters,
              nadir: bool = False) -> Optional[str]:
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
        str, optional: The path to the file to be created or None if the file
        already exists
    """
    first_datetime: datetime.datetime = first_date.astype(
        datetime.datetime)  # type: ignore
    last_datetime: datetime.datetime = last_date.astype(
        datetime.datetime)  # type: ignore
    product_type = "nadir" if nadir else "karin"
    dirname = os.path.join(parameters.working_directory, product_type,
                           first_datetime.strftime("%Y"))
    os.makedirs(dirname, exist_ok=True)
    if nadir:
        filename = (f"SWOT_GPR_2P1P{cycle_number:03d}_{pass_number:03d}_"
                    f"{first_datetime:%Y%m%d}_{first_datetime:%H%M%S}_"
                    f"{last_datetime:%Y%m%d}_{last_datetime:%H%M%S}.nc")
    else:
        product_type = "".join(
            [item.capitalize() for item in parameters.product_type.split("_")])
        filename = (f"SWOT_L2_LR_SSH_{product_type}_"
                    f"{cycle_number:03d}_{pass_number:03d}_"
                    f"{first_datetime:%Y%m%dT%H%M%S}_"
                    f"{last_datetime:%Y%m%dT%H%M%S}_"
                    "DG10_01.nc")
    result = os.path.join(dirname, filename)
    # If the product has already been produced, the generation of this
    # half orbit is disabled.
    if os.path.exists(result):
        LOGGER.debug("%s already generated: %d/%d",
                     "nadir" if nadir else "swath", cycle_number, pass_number)
        return None
    return result


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

    if parameters.nadir:
        nadir_path = file_path(date,
                               last_date,
                               cycle_number,
                               pass_number,
                               parameters,
                               nadir=True)

    # To continue, there must be at least one task left for this pass.
    if swath_path is None and nadir_path is None:
        return

    # Compute the spatial/temporal position of the satellite
    track = orbit_propagator.calculate_pass(pass_number, orbit, parameters)
    if track is None:
        return

    # Set the simulated date
    track.set_simulated_date(date)

    # Mask to set the measurements outside the requirements of the mission to
    # NaN.
    mask = track.mask()

    # The nadir and swath data are concatenated to process the interpolation
    # in one time (swath and nadir).
    lon = np.c_[track.lon, track.lon_nadir[:, np.newaxis]]
    lat = np.c_[track.lat, track.lat_nadir[:, np.newaxis]]
    swath_time = np.repeat(track.time, lon.shape[1]).reshape(lon.shape)

    # Interpolation of the SWH if the user wishes.
    if parameters.swh_plugin is not None:
        swh = parameters.swh_plugin.interpolate(lon.flatten(), lat.flatten(),
                                                swath_time.flatten())
        swh = swh.reshape(lon.shape)
    else:
        swh = None

    # Calculation of instrumental errors
    noise_errors = error_generator.generate(
        cycle_number, pass_number, orbit.curvilinear_distance, track.time,
        track.x_al, track.x_ac,
        swh[:, :-1] if swh is not None else parameters.swh)
    for error in noise_errors.values():
        # Only the swaths must be masked
        if len(error.shape) == 2:
            error *= mask

    # Interpolation of the SSH if the user wishes.
    if parameters.ssh_plugin is not None:
        ssh = parameters.ssh_plugin.interpolate(lon.flatten(), lat.flatten(),
                                                swath_time.flatten())
        ssh = ssh.reshape(lon.shape)
    else:
        ssh = None

    if swath_path:
        LOGGER.info("generate swath %d/%d [%s, %s]", cycle_number, pass_number,
                    track.time[0], track.time[-1])

        # Create the swath dataset
        product = netcdf.Swath(track, parameters.central_pixel,
                               parameters.product_type)

        if ssh is not None:
            product.ssh((ssh[:, :-1] * mask) + sum_error(noise_errors))
            if noise_errors:
                product.simulated_true_ssh(ssh[:, :-1])
        if swh is not None:
            product.swh((swh[:, :-1] * mask))

        product.update_noise_errors(noise_errors)
        product.to_netcdf(cycle_number, pass_number, swath_path,
                          parameters.complete_product)

    # Create the nadir dataset
    if nadir_path:
        LOGGER.info("generate nadir %d/%d [%s, %s]", cycle_number, pass_number,
                    track.time[0], track.time[-1])

        product = netcdf.Nadir(track)

        if ssh is not None:
            product.ssh(ssh[:, -1] + sum_error(noise_errors, swath=False))
            if noise_errors:
                product.simulated_true_ssh(ssh[:, -1])
        if swh is not None:
            product.swh(swh[:, -1])
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
        first_date (numpy.datetime64, optional): First date of the simulation.
        last_date (numpy.datetime64, optional): Last date of the simulation.
    """
    first_date = first_date or np.datetime64("now")

    # Displaying Dask client information.
    LOGGER.info(client)

    assert parameters.ephemeris is not None

    # Calculation of the properties of the orbit to be processed.
    with open(parameters.ephemeris, "r") as stream:  # type: TextIO
        orbit = orbit_propagator.calculate_orbit(parameters,
                                                 stream)  # type: ignore

    # Initialization of measurement error generators
    error_generator = generator.Generator(parameters, first_date,
                                          orbit.orbit_duration())

    # Scatter data into distributed memory
    _error_generator = client.scatter(error_generator)
    _parameters = client.scatter(parameters)
    _orbit = client.scatter(orbit)

    dispatch.compute(client,
                     simulate,
                     orbit.iterate(first_date, last_date),
                     error_generator=_error_generator,
                     orbit=_orbit,
                     parameters=_parameters,
                     logging_server=logging_server)


def software_dependencies() -> List[Tuple[str, str]]:
    """Returns the software dependencies of this release"""
    deps = [
        ("conda", lambda moduleule: moduleule.__version__),
        ("dask", lambda module: module.__version__),
        ("distributed", lambda module: module.__version__),
        ("netCDF4", lambda module: module.__version__),
        ("numba", lambda module: module.__version__),
        ("numpy", lambda module: module.__version__),
        ("pip", lambda module: module.__version__),
        ("pyinterp", lambda module: module.__version__),
        ("scipy", lambda module: module.__version__),
        ("setuptools", lambda module: module.__version__),
        ("xarray", lambda module: module.__version__),
        ("swot_simulator", lambda module: module.__version__),
    ]
    result: List[Tuple[str, str]] = [("python", platform.python_version())]
    for module_name, version_getter in deps:
        try:
            if module_name in sys.modules:
                mod = sys.modules[module_name]
            else:
                mod = importlib.import_module(module_name)
            version = version_getter(mod)
            result.append((module_name, version))
        except ImportError:
            pass
    return result


def copy_parameters(parameters: settings.Parameters,
                    path: pathlib.Path) -> None:
    """Copies the overriden parameters into the working directory.

    Args:
        parameters (settings.Parameters): Simulation parameters.
        filename (pathlib.Path): Path to the overriden parameters file used.
    """
    target = pathlib.Path(parameters.working_directory).joinpath(path.name)
    target.parent.mkdir(parents=True, exist_ok=True)
    version = 1
    now = datetime.datetime.now().isoformat()
    while target.exists():
        target = pathlib.Path(str(target.with_suffix("")) + f" ({version}).py")
        version += 1

    with target.open("w") as target_stream:
        with path.open("r") as source_stream:
            target_stream.write(f"""
# This file has been automatically generated by the SWOT simulator on {now}.
# It contains the overriden parameters used to run the simulation.
# Please do not modify it manually.
#
# Dependency versions
# -------------------
""")
            for dependency, version in sorted(software_dependencies()):
                target_stream.write(f"# * {dependency}: {version}\n")
            target_stream.write(source_stream.read())


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
        parameters = settings.Parameters(
            settings.eval_config_file(args.settings.name))

        # Keep track of the overriden parameters
        copy_parameters(parameters, pathlib.Path(args.settings.name))
        launch(client, parameters, logging_server, args.first_date,
               args.last_date)

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
