Getting started 
---------------

To run the simulator:

.. code-block:: bash

   swot_simulator --first-date '<date>' --last-date \
     '<date>' --threads-per-worker <nproc> \
     <parameter_file> 

.. note::

    The format of the dates expected by the ``first-date``  and ``last-date``
    parameters are quite free. You may enter all strings recognized by
    `dateutil.parser <https://dateutil.readthedocs.io/en/stable/parser.html>`_
    parsing. For example:

        * 2021-01-01
        * 2021-12-31T00:00:00

Example:

.. code-block:: bash

   swot_simulator --first-date '2019-08-02T12:00:00' --last-date \
     '2019-08-06T00:00:00' --threads-per-worker 3  params_aviso.py

You can activate the debug mode by using `--debug` option in the command line.

.. note::

    It's possible to use Dask on a cluster to distribute the calculation.
    For example:

    .. code-block:: bash

        swot_simulator --first-date '2019-08-02T12:00:00' --last-date \
            '2019-08-06T00:00:00' params_aviso.py \
            --scheduler-file scheduler_info.json
