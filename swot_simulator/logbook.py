# Copyright (c) 2019 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Logging handlers
----------------
"""
from typing import IO
import logging
import pathlib
import sys


def _config_logger(stream: IO, level: int, name: str) -> logging.Logger:
    """Configures logbook handler"""
    logger = logging.getLogger(name)
    logger.propagate = True
    formatter = logging.Formatter(
        '[%(levelname)1.1s - %(asctime)s] %(message)s',
        datefmt='%b %d %H:%M:%S')
    handler = logging.StreamHandler(stream) if not isinstance(
        stream, logging.StreamHandler) else stream
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    if level:
        logger.setLevel(level)
    return logger


def setup(stream: IO, debug: bool) -> logging.Logger:
    """Setup the logging system"""
    if stream is None:
        stream = sys.stdout
    level = logging.DEBUG if debug else logging.INFO
    # Capture dask.distributed
    _config_logger(stream, level=level, name="distributed")
    return _config_logger(stream,
                          level=level,
                          name=pathlib.Path(__file__).absolute().parent.name)
