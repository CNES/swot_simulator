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
from typing import Dict, List, Tuple
import argparse
import datetime
import importlib
import pathlib
import platform
import sys
import traceback

import dask.distributed
import dateutil.parser
import numpy as np

from .. import exception, launcher, logbook, settings


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


def software_dependencies() -> List[Tuple[str, str]]:
    """Returns the software dependencies of this release"""
    deps = [
        ("conda", lambda module: module.__version__),
        ("dask", lambda module: module.__version__),
        ("distributed", lambda module: module.__version__),
        ("netCDF4", lambda module: module.__version__),
        ("numba", lambda module: module.__version__),
        ("numpy", lambda module: module.__version__),
        ("pyinterp", lambda module: module.__version__),
        ("scipy", lambda module: module.__version__),
        ("xarray", lambda module: module.__version__),
        ("swot_simulator", lambda module: module.__version__),
        # Setup
        ("setuptools", lambda module: module.__version__),
        ("pip", lambda module: module.__version__),
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
    """Copies the overridden parameters into the working directory.

    Args:
        parameters (settings.Parameters): Simulation parameters.
        filename (pathlib.Path): Path to the overridden parameters file used.
    """
    target = pathlib.Path(parameters.working_directory).joinpath(path.name)
    target.parent.mkdir(parents=True, exist_ok=True)
    version = 1
    now = datetime.datetime.now().isoformat()
    name = path.with_suffix("").name
    while target.exists():
        target = target.parent.joinpath(f"{name} ({version}).py")
        version += 1

    with target.open("w") as target_stream:
        with path.open("r") as source_stream:
            target_stream.write(f"""
# This file has been automatically generated by the SWOT simulator on {now}.
# It contains the overridden parameters used to run the simulation.
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
        launcher.launch(client, parameters, logging_server, args.first_date,
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
