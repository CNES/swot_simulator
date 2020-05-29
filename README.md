[![Build Status](https://dev.azure.com/fbriol/swot_simulator/_apis/build/status/CNES.swot-simulator?branchName=master)](https://dev.azure.com/fbriol/swot_simulator/_build/latest?definitionId=2&branchName=master)
[![conda](https://anaconda.org/conda-forge/swot_simulator/badges/installer/conda.svg?service=github)](https://www.anaconda.com/distribution/)
[![downloads](https://anaconda.org/conda-forge/swot_simulator/badges/downloads.svg?service=github)](https://www.anaconda.com/distribution/)
[![platforms](https://anaconda.org/conda-forge/swot_simulator/badges/platforms.svg?service=github)](https://anaconda.org/conda-forge/swot_simulator)
[![latest-release-date](https://anaconda.org/conda-forge/swot_simulator/badges/latest_release_date.svg?service=github)](https://github.com/CNES/swot-simulator/commits/master)
[![license](https://anaconda.org/conda-forge/swot_simulator/badges/license.svg?service=github)](https://opensource.org/licenses/BSD-3-Clause)
[![Binder](https://binder.pangeo.io/badge_logo.svg)](https://binder.pangeo.io/v2/gh/CNES/swot_simulator/master?filepath=notebooks)

# SWOT Simulator for Ocean Science
## Description

This software simulates SWOT observations of sea surface height (SSH) that can
be applied to an ocean general circulation model (OGCM), allowing the
exploration of ideas and methods to optimize information retrieval from the SWOT
Mission in the future. From OGCM SSH inputs, the software generates SWOT-like
outputs on a swath along the orbit ground track and adds measurement error and
noise, which are generated according to technical characteristics published by
the SWOT project team. Not designed to directly simulate the payload instrument
performance, this SWOT simulator aims at providing statistically realistic
outputs for the science community with a simple software package released as an
open source in Python. The software is scalable and designed to support future
evolution of orbital parameters, error budget estimates from the project team
and suggestions from the science community.

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
