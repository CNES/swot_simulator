# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Dispatch task on free workers
=============================
"""
from typing import Any, Callable, Iterator, List, Set
import time
import dask.distributed


def _available_workers(client: dask.distributed.Client) -> Set[str]:
    """Get the list of available workers.

    Args:
        client (dask.distributed.Client): Client connected to the Dask
        cluster.

    Returns:
        list: The list of available workers.
    """
    while True:
        info = client.scheduler_info()
        result = set(info['workers']) - set(
            k for k, v in client.processing().items() if v)
        if result:
            return result
        time.sleep(0.1)


def compute(client: dask.distributed.Client, func: Callable, seq: Iterator,
            *args, **kwargs) -> List[Any]:
    """Distribute the execution of functions to free workers, i.e. those who
    do not perform any tasks.

    Args:
        client (dask.distributed.Client): Client connected to the Dask
            cluster.
        func (callable): Function to execute
        seq (iterable): The sequence of arguments handled by ``func``.
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Returns:
        list: The result of the execution of the functions.
    """
    completed = dask.distributed.as_completed()
    workers = set()
    result = []
    iterate = True
    # As long as there is data to traverse in the iterator.
    while iterate:
        # As long as there are free workers
        while completed.count() < len(client.scheduler_info()['workers']):
            try:
                if not workers:
                    workers = _available_workers(client)
                item = next(seq)
                completed.add(
                    client.submit(func,
                                  item,
                                  *args,
                                  workers=workers.pop(),
                                  allow_other_workers=False,
                                  **kwargs))
            except StopIteration:
                iterate = False
                break
        # The computation queue is full, we consume the finished jobs to be
        # able to continue.
        if iterate:
            try:
                result += client.gather(completed.next_batch())
            except StopIteration:
                pass
    result += [item.result() for item in completed]
    return result
