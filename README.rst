################################
SWOT Simulator for Ocean Science
################################

|Latest Documentation Status| |Conda| |Downloads statistics| |Platforms|
|Latest Release Date| |License| |Binder|

This software simulates SWOT observations of the sea surface height (SSH) that
can be applied to an ocean general circulation model (OGCM), allowing the
exploration of ideas and methods to optimize information retrieval from the SWOT
Mission in the future. From OGCM SSH inputs, the software generates SWOT-like
outputs on a swath along the orbit ground track and adds measurement errors and
noise, which are generated according to technical characteristics published by
the SWOT project team. Not designed to directly simulate the payload instrument
performance, this SWOT simulator aims at providing statistically realistic
outputs for the science community with a simple software package released as an
open source in Python. The software is scalable and designed to support future
evolution of orbital parameters, error budget estimates from the project team
and suggestions from the science community.

Tutorial and reference documentation is provided at
`swot-simulator.readthedocs.io <https://swot-simulator.readthedocs.io>`_.

About
-----

This project was created by Clement Ubelmann, Lucile Gaultier and Lee-Lueng Fu.

Jet Propulsion Laboratory, California Institute of Technology, CNES

License
-------

``swot_simulator`` is provided under a BSD-style license that can be found in
the `LICENSE <https://github.com/CNES/swot_simulator/blob/master/LICENSE>`_
file. By using, distributing, or contributing to this project, you agree to the
terms and conditions of this license.

.. |Latest Documentation Status| image:: https://dev.azure.com/fbriol/swot_simulator/_apis/build/status/CNES.swot_simulator?branchName=master
    :target: https://dev.azure.com/fbriol/swot_simulator/_build/latest?definitionId=2&branchName=master
.. |Conda| image:: https://anaconda.org/conda-forge/swot_simulator/badges/installer/conda.svg?service=github
    :target: https://www.anaconda.com/distribution/
.. |Downloads statistics| image:: https://anaconda.org/conda-forge/swot_simulator/badges/downloads.svg?service=github
    :target: https://www.anaconda.com/distribution/
.. |Platforms| image:: https://anaconda.org/conda-forge/swot_simulator/badges/platforms.svg?service=github
    :target: https://anaconda.org/conda-forge/swot_simulator
.. |Latest Release Date| image:: https://anaconda.org/conda-forge/swot_simulator/badges/latest_release_date.svg?service=github
    :target: https://github.com/CNES/swot_simulator/commits/master
.. |License| image:: https://anaconda.org/conda-forge/swot_simulator/badges/license.svg?service=github
    :target: https://opensource.org/licenses/BSD-3-Clause
.. |Binder| image:: https://binder.pangeo.io/badge_logo.svg
    :target: https://binder.pangeo.io/v2/gh/CNES/swot_simulator/master?filepath=notebooks
