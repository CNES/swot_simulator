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

The plugin for the SSH is specified in the
:py:const:`ssh_plugin <settings.ssh_plugin>` key, and the one for the swh in
the :py:const:`swh_plugin <settings.swh_plugin>` key.

To disable interpolation of SWH or SSH models, the related parameters must be
set to ``None``.

Plugin versioning
=================

The simulator keeps track of the code version and parameters used to generate a
dataset. It generates a version.py file containing the simulator version, the
SSH plugin class and version, and the SWH plugin class and version. The
:func:`~swot_simulator.plugins.Interface.version()` method should be overriden
in the plugins to keep track of the changes in the plugins.
