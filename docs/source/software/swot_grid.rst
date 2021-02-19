Generation of the SWOT grid
----------------------------

:mod:`swot_simulator.orbit_propagator` module generates the SWOT grid. The
:py:const:`ephemeris <settings.ephemeris>` keyword sets the path to the orbit
file to process. This file defines three different columns (longitude,
latitude, and the corresponding time) the position of the satellite at each
point of the orbit. You can specify column order with the parameters
:py:const:`ephemeris_col <settings.ephemeris_col>` (default is ``[1, 2, 0]``).
The orbit is interpolated
at the along-track resolution specified by the user (parameter
:py:const:`delta_al <settings.delta_al>`). At each point of the orbit, the
swath is calculated with an across-track resolution defined by the
:py:const:`delta_ac <settings.delta_ac>` parameter.

The width of the swath (:py:const:`half_swath <settings.half_swath>`) and the
gap between the nadir and the swath (:py:const:`half_gap <settings.half_gap>`)
can also be defined according to :ref:`Fig. 2 <Fig2>`

By default, all measurements of the swath are valid. In reality, measurements
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
    :ref:`Proposed SWOT orbits <ProposedSWOTOrbits>`
