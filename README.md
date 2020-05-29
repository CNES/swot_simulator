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

## Documentation

Learn more about xarray in its official documentation at https://swot-simulator.readthedocs.io/
