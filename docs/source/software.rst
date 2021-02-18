The software
=============

The software is written in Python3, and is mainly tested on python3.8 and 3.9.
All the parameters that can be modified by the user are read in a params file
(e.g. params.py) specified by the user.

The software is divided in 6 main modules: 

+-------------------------------------------+----------------------------------+
| Python module                             | Description                      |
+===========================================+==================================+
|:mod:`swot_simulator.launcher`             | Entry point of the main program. |
+-------------------------------------------+----------------------------------+
|:mod:`swot_simulator.orbit_propagator`     | Generates the SWOT swath.        |
+-------------------------------------------+----------------------------------+
|:mod:`swot_simulator.error.generator`      | Generates all the errors on the  |
|                                           | swath.                           |
+-------------------------------------------+----------------------------------+
|:mod:`swot_simulator.product_specification`| Defines the format to save output|
|                                           | in L2 SWOT-like data.            |
+-------------------------------------------+----------------------------------+
|:mod:`swot_simulator.settings`             | Handles the simulation           |
|                                           | parameters.                      |
+-------------------------------------------+----------------------------------+
|:mod:`swot_simulator.plugins`              | Modules containing plugins to    |
|                                           | readand interpolate SSH and SWH  |
|                                           | models.                          |
+-------------------------------------------+----------------------------------+

Inputs
-------

To read and interpolate SSH and SWH inputs on the SWOT grid, you have to provide
your own reader as a plugin. Several readers are already available in the
`plugins/ssh
<https://github.com/CNES/swot_simulator/tree/master/swot_simulator/plugins/ssh>`_
directory for the SSH:

+--------------------------------------------+----------------------------------+
| Python module                              | Description                      |
+============================================+==================================+
|:mod:`swot_simulator.plugins.ssh.aviso`     | AVISO/CMEMS data, which are on   |
|                                            | regular grid with one timestep   |
|                                            | per file.                        |
+--------------------------------------------+----------------------------------+
|:mod:`swot_simulator.plugins.ssh.hycom`     | HYCOM regional data.             |
+--------------------------------------------+----------------------------------+
|:mod:`swot_simulator.plugins.ssh.mitgcm`    | LLC4320 MITGCM zarr data         |
|                                            | (format compatible with the one  |
|                                            | available on the CNES HAL        |
|                                            | supercomputer)                   |
+--------------------------------------------+----------------------------------+
|:mod:`swot_simulator.plugins.ssh.mitgcm_ww3`| LLC4320 MITGCM data converted in |
|                                            | netCDF on Datarmor IFREMER       |
|                                            | supercomputer and compatible with|
|                                            | WW3 format.                      |
+--------------------------------------------+----------------------------------+

One reader is available in the `plugins/swh
<https://github.com/CNES/swot_simulator/tree/master/swot_simulator/plugins/swh>`_
directory for the SWH:

+-------------------------------------+----------------------------------------+
| Python module                       | Description                            |
+=====================================+========================================+
|:mod:`swot_simulator.plugins.swh.ww3`| WW3 outputs (available on Datarmor     |
|                                     | IFREMER supercomputer).                |
+-------------------------------------+----------------------------------------+

The plugin for the ssh is specified in the :py:const:`ssh_plugin
<settings.ssh_plugin>` key, and the one for the swh in the :py:const:`swh_plugin
<settings.swh_plugin>` key.

To disable interpolation of SWH or SSH models, the related parameters must be
set to ``None``.

 
Generation of the SWOT grid
----------------------------


:mod:`swot_simulator.orbit_propagator` module generates the SWOT grid. The
:py:const:`ephemeris <settings.ephemeris>` keyword sets the path to the used
orbit file. This file defines three different columns (longitude, latitude,
and the corresponding time) the position of the satellite at each point of the
orbit. You can specify column order with the parameters :py:const:`ephemeris_col
<settings.ephemeris_col>` (default is ``[1, 2, 0]``). The orbit is interpolated
at the along-track resolution specified by the user (parameter
:py:const:`delta_al <settings.delta_al>`). At each point of the orbit, the swath
is calculated with an across-track resolution defined by the :py:const:`delta_ac
<settings.delta_ac>` parameter.


The width of the swath (:py:const:`half_swath <settings.half_swath>`) and the
gap between the nadir and the swath (:py:const:`half_gap <settings.half_gap>`)
can also be defined according to :ref:`Fig. 2 <Fig2>`

By default, all measurements on the swath are valid. In reality, measurements
outside the requirements of the SWOT swath are subject to errors. To simulate
this, you can define the limits of the swath validity with the
:py:const:`half_gap <settings.requirement_bounds>` parameter. Measurements
outside the span will be set with NaN.


The generation of SWOT products is global, i.e., the products cover the globe's
entire surface. With the help of the parameter :py:const:`area <settings.area>`,
you can define a more restricted area of interest. The expected syntax of the
parameter to define a bounding box is: ``area = [lon_min, lat_min, lon_max,
lat_max]``.

.. seealso::
    :ref:`ProposedSWOTorbits`

Sampled SSH and error fields
----------------------------

At each pass, for each cycle, an output netCDF file containing the SSH
interpolated from the model (if :py:const:`ssh_plugin <settings.ssh_plugin>` is
not set to `None`), and the different errors. The naming and the file format are
compliant with SWOT's PDD (Product Description Document), set
:py:const:`complete_product <settings.complete_product>` key to ``True`` to
store all the variables and ``False`` to store only the one compute by the
simulator.

You can also select the output type :py:const:`product_type
<settings.product_type>` among ``expert``, ``ssh``, and ``windwave``. Default is
``expert`` and details are in the corresponding xml file.

The output directory is defined in :py:const:`working_directory
<settings.working_directory>` key, default is the user`s root. The output file
names are stored in the output directory in ``karin/<year>`` directory for SWOT
and ``nadir/<year>`` directory for the nadir. The naming follows the pattern
``SWOT_L2_LR_Expert_[cycle]_[pass]_[start_time]_[stop_time]_DG10_01.nc`` for the
Karin products and ``SWOT_GPN_2P1P_[cycle]_[pass]_[start_time]_[stop_time].nc``
for the nadir products. The interpolation method is specified in your SSH
plugin. `pangeo-pyinterp <https://github.com/CNES/pangeo-pyinterp>`_ module is
used in the examples.

The :py:const:`nadir <settings.nadir>` and :py:const:`karin <settings.karin>`
parameters enable (parameter set to ``True``) or disable (parameter set to
``False``) the generation of these products.

Computation of :py:const:`errors <settings.noise>` are specified as a list:
``noise = ['altimeter', 'baseline_dilation', 'karin', 'corrected_roll_phase',
'timing', 'wet_troposphere']`` ``corrected_roll_phase`` error will generate
already cross-calibrated roll-phase error whereas ``roll_phase`` will generate
errors before cross-calibration.

A :py:const:`repeat length <settings.len_repeat>` ``len_repeat = 2000`` key
defines the wavelength in km to repeat the noise and is used for all noises. The
path to the file that contains instrumental errors is mentioned in
:py:const:`error_spectrum <settings.error_spectrum>` key. It is possible to
specify the seed for Randomstate in :py:const:`nseed <settings.nseed>`. The
following parameters are specific to each noise component:

* KaRIN: The file noise that contains the spectrum is available in
  :py:const:`karin_noise <settings.karin_noise>`, a constant :py:const:`swh
  <settings.swh>` can be set ``swh=2`` or interpolated from a model SWH
  specified in the plugin :py:const:`swh_plugin <settings.swh_plugin>`.

* Roll-phase: To use the already cross-calibrated roll phase, download the file
  from the FTP (https://ftp.odl.bzh/swot) and specify the path in the
  :py:const:`corrected_roll_phase_dataset
  <settings.corrected_roll_phase_dataset>` key. So far, the following files are
  available:

     * ``data_sim_slope_v0.nc``: One year of cross-calibrated roll-phase
     * ``data_sim_slope_2cycles_v0.nc``: Two cycles of cross-calibrated
       roll-phase

* Wet troposphere: The number of beams used to correct the wet troposphere
  is set in :py:const:`nbeam <settings.nbeam>` variable. The beam position of
  each beam can be set as a list in :py:const:`beam_position
  <settings.beam_position>`, and the gaussian footprint can be changed using
  :py:const:`sigma <settings.sigma>` variable.

Nadir
-----

The SWOT simulator can be used to generate any nadir like observation. The
mission :py:const:`ephemeris <settings.ephemeris>` should be provided as
explained in Generation of the SWOT grid section. The corresponding
:py:const:`cycle_duration <settings.cycle_duration>` and :py:const:`height
<settings.height>` value should be specified. Default values are for the SWOT
mission: 

* ``cycle_duration=20.86455`` days

* ``height=891000`` m

Getting started 
---------------

To run the simulator:

.. code-block:: bash

   swot_simulator --first-date '<date>' --last-date \
     '<date>' --threads-per-worker <nproc> \
     <parameter_file> 

.. note ::

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

.. note ::

    It's possible to use Dask on a cluster to distribute the calculation.
    For example:

    .. code-block:: bash

        swot_simulator --first-date '2019-08-02T12:00:00' --last-date \
            '2019-08-06T00:00:00' params_aviso.py \
            --scheduler-file scheduler_info.json


Sample configuration file
-------------------------

.. _params-file:

.. literalinclude:: settings.py
