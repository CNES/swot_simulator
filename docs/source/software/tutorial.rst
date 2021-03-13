Getting started
===============

.. _main_program:

Main program
------------

The software is structured like a library. It is driven by a
:py:func:`main program <swot_simulator.launcher.main>` that you can launch to
generate your simulations. When this Python module is installed, the main
program is called ``swot_simulator``. The program accepts options to alter its
default behavior and a :ref:`configuration file <params-file>`. The ``--help``
option allows to display the online help of the command:

.. code-block:: text

    usage: swot_simulator [-h] [--first-date FIRST_DATE] [--last-date LAST_DATE]
                          [--debug] [--log PATH] [--scheduler-file PATH]
                          [--n-workers N] [--processes] [--threads-per-worker N]
                        settings

    positional arguments:
        settings              Path to the parameters file

    optional arguments:
        -h, --help            show this help message and exit

    General:
        Simulation general settings

        --first-date FIRST_DATE
                            The first date to be processed. Default to the
                            current date
        --last-date LAST_DATE
                            The last date to be processed. Default to the last
                            date allowing to cover an entire cycle.

    Execution:
        Runtime parameters options

        --debug             Put swot simulator in debug mode
        --log PATH          Path to the logbook to use
        --scheduler-file PATH
                            Path to a file with scheduler information to launch
                            swot simulator on a cluster. By default, use a local
                            cluster.

    LocalCluster:
        Dask local cluster option

        --n-workers N       Number of workers to start (Default to 1)
        --processes         Whether to use processes (True) or threads (False).
                            Defaults to False
        --threads-per-worker N
                            Number of threads per each worker. (Default to 1)

    Configuration:
        --template PATH     Writes the default configuration of the simulator
                            into the file and ends the program.

The first step is to create a configuration file for the simulator. To do this,
execute the following command:

.. code-block::

    swot_simulator --template conf.py

The program will create the configuration file ``conf.py`` and ends the program.
The  :py:const:`ephemeris <settings.ephemeris>`,
:py:const:`error_spectrum <settings.error_spectrum>`,
:py:const:`karin_noise <settings.karin_noise>` parameters are filled using the
data delivered with this simulator in order to simulate the swaths on the
scientific orbit without interpolating SSH. In other words, this setting is
sufficient to launch the simulator and generate the various errors, but the
generated products will not contain any interpolated SSH under the satellite
swath.

To start the simulator execute the following command.

.. code-block:: bash

   swot_simulator conf.py

The program will launch and generate a product cycle for the selected orbit. The
first half orbit will begin on the date the simulator is launched. You can
change this behavior by setting a start and end date for the simulation using
the ``first-date`` and ``last-date`` parameters.

.. note::

    The format of the dates expected by the ``first-date``  and ``last-date``
    parameters are quite free. You may enter all strings recognized by
    `dateutil.parser <https://dateutil.readthedocs.io/en/stable/parser.html>`_
    parsing. For example:

    * 2021-01-01
    * 2021-12-31T00:00:00

The data will be written in the directory specified by the
:py:const:`working_directory <settings.working_directory>` parameter.

By default, the program generates the products on a single thread and process.
You can change this behavior by modifying the program's ``n-workers``,
``processes`` and ``threads-per-worker`` options.

If you can run the program on an HPC, you can generate the products using a dask
cluster. To do this, you need to start the dask cluster, write the cluster
properties to a JSON file, for example the ``scheduler.json`` file, and then run
the program by specifying this file to simulator using the ``scheduler-file``
option.

Library
-------

It is also possible to generate swaths using the Python module directly. The
following `link
<https://binder.pangeo.io/v2/gh/CNES/swot_simulator/master?filepath=notebooks>`_
allows you to launch a notebok on the `Binder <https://jupyter.org/binder>`_
`Pangeo <https://pangeo.io/>`_ to show you how to simulate SWOT data in a
Jupyter notebook.
