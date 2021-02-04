# Copyright (c) 2020 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Dispatch task on free workers
=============================
"""
from typing import Any, Callable, Iterator, List, Tuple
import time
import dask.distributed


def _available_workers(client: dask.distributed.Client) -> Tuple[int, List]:
    """Get the list of available workers.

    Args:
        client (dask.distributed.Client): Client connected to the Dask
        cluster.

    Returns:
        tuple: The number of workers and the list of available workers.
    """
    while True:
        info = client.scheduler_info()
        executing = [
            k for k, v in info['workers'].items()
            if v['metrics']['executing'] == 0
        ]
        if executing:
            workers = len(info['workers'])
            return workers, executing
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
        *args, **kwargs : any
            Extra arguments and keyword arguments to pass to ``func``.

    Returns:
        list: The result of the execution of the functions.
    """
    completed = dask.distributed.as_completed()
    result = []
    iterate = True
    workers, free_workers = _available_workers(client)
    # As long as there is data to traverse in the iterator.
    while iterate:
        # As long as there are free workers
        while completed.count() < workers:
            try:
                if not free_workers:
                    workers, free_workers = _available_workers(client)
                item = next(seq)
                completed.add(
                    client.submit(func,
                                  item,
                                  *args,
                                  workers=free_workers.pop(),
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
