The software
=============

The software is written in Python3, and is mainly tested on python 3.7 and 3.8.
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

.. toctree::
   :maxdepth: 1

   software/input.rst
   software/swot_grid.rst
   software/ssh.rst
   software/nadir.rst
   settings
   software/settings.rst
