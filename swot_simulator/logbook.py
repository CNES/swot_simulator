# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Logging handlers
----------------
"""
from typing import Awaitable, IO, Optional, Tuple, Union
import logging
import logging.handlers
import pathlib
import pickle
import socket
import struct
import sys
import threading
import tornado.ioloop
import tornado.iostream
import tornado.log
import tornado.tcpserver

#: Synchronize logs for workers
LOCK = threading.RLock()


class LogRecordSocketReceiver(tornado.tcpserver.TCPServer):
    """
    Simple TCP socket-based logging receiver.

    Args:
        name (str, optional): Logger name.
        **kwargs: Base class keyword arguments.
    """
    def __init__(self, name: Optional[str] = None, **kwargs):
        super().__init__(**kwargs)
        self.logname = name

    async def handle_stream(self, stream: tornado.iostream.IOStream,
                            address: Tuple) -> Optional[Awaitable[None]]:
        """Override to handle a new IOStream from an incoming connection.

        Args:
            stream (tornado.iostream.IOStream): Incoming stream.
            address (tuple): IP address and TCP port.
        """
        while True:
            try:
                chunk = await stream.read_bytes(4)
                if len(chunk) < 4:
                    break
                slen = struct.unpack('>L', chunk)[0]
                chunk = await stream.read_bytes(slen)
                while len(chunk) < slen:
                    chunk += await stream.read_bytes(slen - len(chunk))
            except tornado.iostream.StreamClosedError:
                break
            obj = pickle.loads(chunk)
            record = logging.makeLogRecord(obj)
            setattr(record, "ip", address[0])
            self.handle_log_record(record)

    def handle_log_record(self, record: logging.LogRecord) -> None:
        """Handle an incoming logging reccord

        Args:
            record (logging.LogRecord): The event being logged.
        """
        # if a name is specified, we use the named logger rather than the one
        # implied by the record.
        if self.logname is not None:
            name = self.logname
        else:
            name = record.name
        logger = logging.getLogger(name)

        # N.B. EVERY record gets logged. This is because Logger.handle
        # is normally called AFTER logger-level filtering. If you want
        # to do filtering, do it at the client end to save wasting
        # cycles and network bandwidth!
        logger.handle(record)


class LogServer:
    """Handle the log server

    Args:
        hostname (str, optional): Specifies the IP address we want to listen
            to.
        port (int, optional): Specifies the port we want to listen to.
    """
    def __init__(
            self,
            hostname: Optional[str] = None,
            port: int = logging.handlers.DEFAULT_TCP_LOGGING_PORT) -> None:
        hostname = hostname or socket.gethostname()
        self.ip = socket.gethostbyname(hostname)
        self.port = port

        server = LogRecordSocketReceiver()
        server.listen(self.port, self.ip)

        ioloop = tornado.ioloop.IOLoop.current()
        self.thread = threading.Thread(target=ioloop.start)

    def start(self, daemon: bool = False) -> None:
        """Start the server

        Args:
            daemon (bool, optional): Indicates whether this thread is a daemon
                thread (True) or not (False)
        """
        if daemon:
            self.thread.daemon = True
        self.thread.start()

    def __iter__(self):
        yield self.ip
        yield self.port


class LogFormatter(tornado.log.LogFormatter):
    """Inserts the IP address of the worker or scheduler.

    Args:
        *args: Base class arguments.
        **kwarsg: Base class keyword arguments.
    """
    def __init__(self, *args, **kwargs):
        hostname = socket.gethostname()
        self._ip = socket.gethostbyname(hostname)
        super().__init__(*args, **kwargs)

    def format(self, record: logging.LogRecord) -> str:
        """Do formatting for a record.

        Args:
            record (logging.Record): information to the event being logged.

        Returns:
            str: The record formatted.
        """
        if "ip" not in record.__dict__:
            record.__dict__['ip'] = self._ip
        return super().format(record)


def _config_logger(stream: Union[IO[str], logging.Handler], level: int,
                   name: str) -> logging.Logger:
    """Configures logbook handler"""
    logger = logging.getLogger(name)
    logger.propagate = True
    formatter = LogFormatter(
        '%(color)s[%(levelname)1.1s - %(ip)s - %(asctime)s - %(module)s]'
        '%(end_color)s %(message)s',
        datefmt='%b %d %H:%M:%S')
    handler = logging.StreamHandler(stream) if not isinstance(
        stream, logging.Handler) else stream
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    if level:
        logger.setLevel(level)
    return logger


def setup(stream: IO[str],
          debug: bool) -> Tuple[logging.Logger, Tuple[str, int, int]]:
    """Setup the logging system

    Args:
        stream (io): Flux used to write in the log.
        debug (bool): True if the log should record debugging traces.

    Returns:
        tuple: The channel logging and log server settings.
    """
    logging_server = LogServer()
    logging_server.start(True)

    if stream is None:
        stream = sys.stdout
    level = logging.DEBUG if debug else logging.INFO
    # Capture dask.distributed
    _config_logger(stream, level=level, name="distributed")
    return _config_logger(
        stream,
        level=level,
        name=pathlib.Path(__file__).absolute().parent.name), (
            logging_server.ip, logging_server.port, level)


def setup_worker_logging(logging_server: Tuple[str, int, int]):
    """Setup the logging server to log worker calculations

    Args:
        logging_server (tuple): Log server connection settings.
    """
    with LOCK:
        for name in [
                "distributed",
                pathlib.Path(__file__).absolute().parent.name
        ]:
            logger = logging.getLogger(name)
            # If this logger is already initialized we do nothing
            if logger.handlers:
                continue
            stream = logging.handlers.SocketHandler(logging_server[0],
                                                    logging_server[1])
            _config_logger(stream, level=logging_server[2], name=name)
