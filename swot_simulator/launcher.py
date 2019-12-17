# Copyright (c) 2019 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Main program
------------
"""
from typing import IO, Optional
import argparse
import datetime
import logging
import os
import sys
import traceback
import dask.distributed
import dateutil
import numpy as np
from . import exception
from . import logbook
from . import product_specification
from . import orbit_propagator
from . import settings

#: Logger of this module
LOGGER = logging.getLogger(__name__)


def datetime_type(value):
    """The option should define a datetime"""
    try:
        value = dateutil.parser.parse(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            f"invalid date time {value!r}: {error!s}")
    return np.datetime64(value)


def writable_directory(value):
    """The option should define a writable directory"""
    if not os.path.exists(value):
        raise argparse.ArgumentTypeError(
            f"no such file or directory: {value!r}")
    if not os.path.isdir(value):
        raise argparse.ArgumentTypeError(f"{value!r} is not a directory")
    if not os.access(value, os.W_OK):
        raise argparse.ArgumentTypeError(
            f"cannot open directory {value!r} for writing")
    return value


def usage() -> argparse.Namespace:
    """Command syntax"""
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
    group.add_argument("--threads-per-worker",
                       help="Number of threads per each worker. "
                       f"Defaults to {os.cpu_count()}",
                       type=int,
                       metavar='N',
                       default=os.cpu_count())
    group.add_argument("--scheduler-file",
                       help="Path to a file with scheduler information to "
                       "launch swot simulator on a cluster. By "
                       "default, use a local cluster.",
                       metavar='PATH',
                       type=argparse.FileType("r"))

    return parser.parse_args()


def file_path(date: np.datetime64,
              cycle_number: int,
              pass_number: int,
              working_directory: str,
              nadir: Optional[bool] = False):
    """Get the absolute path of the file to be created."""
    date = datetime.datetime.utcfromtimestamp(
        date.astype("datetime64[s]").astype("int64"))
    product_type = "nadir" if nadir else "karin"
    dirname = os.path.join(working_directory, product_type,
                           date.strftime("%Y"))
    os.makedirs(dirname, exist_ok=True)

    return os.path.join(
        dirname,
        f"swot_{product_type}_c{cycle_number:03d}_p{pass_number:03d}.nc")


def simulate(cycle_number: int, pass_number: int, date: np.datetime64,
             orbit: orbit_propagator.Orbit,
             parameters: settings.Parameters) -> None:
    LOGGER.info(f"generate pass {cycle_number}/{pass_number}")

    # Compute the spatial/temporal position of the satellite
    track = orbit_propagator.calculate_pass(pass_number, orbit, parameters)
    track.time = date

    if parameters.swath:
        # Create the swath dataset
        product = product_specification.Swath(track)

        # If the file has already been generated, the other operations are
        # ignored
        path = file_path(date, cycle_number, pass_number,
                         parameters.working_directory)
        if not os.path.exists(path):

            # Interpolation of the SSH if the user wishes.
            if parameters.ssh_plugin is not None:
                swath_time = np.repeat(track.time, track.lon.shape[1]).reshape(
                    track.lon.shape)
                product.ssh(
                    parameters.ssh_plugin.interpolate(
                        track.lon.flatten(), track.lat.flatten(),
                        swath_time.flatten()).reshape(track.lon.shape))

            product.to_netcdf(cycle_number, pass_number, path)

    # Create the nadir dataset
    if parameters.nadir:
        product = product_specification.Nadir(track)

        # If the file has already been generated, the other operations are
        # ignored
        path = file_path(date,
                         cycle_number,
                         pass_number,
                         parameters.working_directory,
                         nadir=True)
        if not os.path.exists(path):

            # Interpolation of the SSH if the user wishes.
            if parameters.ssh_plugin is not None:
                product.ssh(
                    parameters.ssh_plugin.interpolate(track.lon_nadir,
                                                      track.lat_nadir,
                                                      track.time))

            product.to_netcdf(cycle_number, pass_number, path)


def launch(client: dask.distributed.Client,
           parameters: settings.Parameters,
           first_date: Optional[datetime.datetime] = None,
           last_date: Optional[datetime.datetime] = None):
    # Displaying Dask client information.
    LOGGER.info(client)

    #: pylint: disable=no-member
    # The attribute "ephemerid" is set dynamically

    # Calculation of the properties of the orbit to be processed.
    with open(parameters.ephemeris) as stream:
        orbit = orbit_propagator.calculate_orbit(parameters, stream)
    #: pylint: enable=no-member

    if last_date is None:
        # By default, only one cycle is processed.
        last_date = first_date + ((orbit.time[-1]).astype(np.int64) *
                                  1000000).astype("timedelta64[us]")
    else:
        last_date = last_date

    # For the moment, we start the processing on the first pass
    absolute_track = 1

    # For the entire period to be generated, the generation of orbits is
    # assigned to dask.
    futures = []

    _parameters = client.scatter(parameters)
    _orbit = client.scatter(orbit)

    date = first_date or datetime.datetime.now()
    while date <= last_date:
        cycle, track = orbit.decode_absolute_pass_number(absolute_track)

        # Generation of the simulated product.
        futures.append(
            client.submit(simulate, cycle, track, date, _orbit, _parameters))

        # Shift the date of the duration of the generated pass
        date += orbit.pass_duration(track)

        # Update of the number of the next pass to be generated
        absolute_track += 1

    client.gather(futures)


def main():
    """main function"""
    args = usage()

    # Setup log
    logger = logbook.setup(args.log, args.debug)

    client = dask.distributed.Client(
        dask.distributed.LocalCluster(
            n_workers=1,
            processes=False,
            threads_per_worker=args.threads_per_worker)
    ) if args.scheduler_file is None else dask.distributed.Client(
        scheduler_file=args.scheduler_file.name)

    try:
        parameters = settings.eval_config_file(args.settings.name)
        launch(client, settings.Parameters(parameters), args.first_date,
               args.last_date)

        client.close()
        logger.info("End of processing.")
        return 0
    except Exception as exc:
        # Clients are stopped before writing the error message to the log to
        # ensure that the log will end with the exception captured.
        client.close()
        logger.error(
            exception.structured_traceback(
                exc, traceback.extract_tb(sys.exc_info()[2])))
        logger.error("End of processing.")
    return 1


if __name__ == '__main__':
    sys.exit(main())
