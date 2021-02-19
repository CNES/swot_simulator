Installation
============

Required dependencies
---------------------

- Python (3.6 or later)
- setuptools
- `python-dateutil <https://github.com/dateutil/dateutil>`_
- `dask.distributed <https://github.com/dask/distributed>`_
- `netCDF4 <https://github.com/Unidata/netcdf4-python>`_
- `numba <https://github.com/numba/numba>`_
- `numpy <http://www.numpy.org/>`__ (1.15 or later)
- `pyinterp <https://github.com/CNES/pangeo-pyinterp>`_
- `scipy <https://github.com/scipy/scipy>`__
- `xarray <https://github.com/pydata/xarray>`__


Instructions
------------
This software is a pure Python package, but its dependencies are not. The
easiest way to get everything installed is to use conda_. To install
swot_simulator with its recommended dependencies using the conda command line
tool::

    $ conda install -c conda-forge swot_simulator

.. _conda: http://conda.io/

If you don't use conda, be sure you have the required dependencies installed
first. Then, install swot_simulator with pip::

    $ pip install swot_simulator

Testing
-------

To run the test suite after installing swot_simulator, install (via pypi or
conda) `pytest <https://pytest.org>`__ and run ``pytest`` in the root
directory of the swot_simulator repository.
