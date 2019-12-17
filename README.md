# TODO

## Running the software.

This program is executed using the ``swot_simulator.py`` script. The syntax of
the command is as follows:

    swot_simulator.py settings.py

> This script is not present in the source code, because it is a Python entry
> point generated during the installation of the package.

The parameter of this script is a configuration file, defined in Python,
describing the simulation parameters. A sample file is available
[here](docs/source/settings.py).

The script options allow you to modify the duration of the simulation, the
level of verbiage of the script and to execute this program on a cluster using
[dask](https://docs.dask.org/en/latest/). The configuration file can define an
SSH interpolation plug-in. The software provides two plug-ins to interpolate
the SSH from the
[AVISO](http://marine.copernicus.eu/services-portfolio/access-to-products/?option=com_csw&view=details&product_id=SEALEVEL_GLO_PHY_CLIMATE_L4_REP_OBSERVATIONS_008_057)
and [MIT/GCM](http://online.kitp.ucsb.edu/online/blayers18/menemenlis/) grids.
