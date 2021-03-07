.. _ProposedSWOTOrbits:

Proposed SWOT orbits
--------------------

The software uses as an input the ground-tracks of the satellite orbit. The user
can choose between different orbits such as the fast sampling orbit (1-day
repeat), the science orbit (21-day repeat with a 10-day subcycle), and also the
contingency orbit (21-day repeat with 1-day subcycle). The table below shows the
characteristics of these three orbits:

+---------------------+--------------+--------------+------------+-------------+-----------+
|                     | Repeat Cycle | Repeat Cycle | Sub-cycles | Inclination | Elevation |
|                     | (days)       | (Orbits)     | (days)     |             | (km)      |
+=====================+==============+==============+============+=============+===========+
| Fast Sampling orbit | 0.99349      | 14           | N.A.       | 77.6        | 857       |
+---------------------+--------------+--------------+------------+-------------+-----------+
| Science Orbit       | 20.8646      | 292          | 1, 10      | 77.6        | 891       |
+---------------------+--------------+--------------+------------+-------------+-----------+
| Contingency orbit   | 20.8639      | 293          | 1          | 77.6        | 874       |
+---------------------+--------------+--------------+------------+-------------+-----------+

The ground-track coordinates corresponding to these orbits are given as input
ASCII files of 3 columns (longitude, latitude, time) for one complete cycle
sampled at every  ~5~km. The ascending node has been arbitrarily set to zero
degree of longitude, but the user can shift the orbit of any value in longitude.

Orbit files have been updated with the one provided by AVISO_ in September 2015
(https://www.aviso.altimetry.fr/en/missions/future-missions/swot/orbit.html).
There are two additional orbit files available in the last version of the
simulator. Input files are also ASCII with 3 columns (time, longitude,
latitude). Orbits are provided at low resolution and are interpolated
automatically by the simulator. `ephemeris_calval_june2015_ell.txt
<https://github.com/CNES/swot_simulator/blob/master/data/ephem_calval_june2015_ell.txt>`_
contains the updated fast sampling orbit and `ephemeris_science_sept2015_ell.txt
<https://github.com/CNES/swot_simulator/blob/master/data/ephem_science_sept2015_ell.txt>`_
the updated science orbit.

Other orbit files of the same format (time, longitude, latitude) can also be
used as an input. To avoid distortions in the SWOT grid, we recommend a minimum
of 10 km sampling between the ground-track points of the orbit.

.. _AVISO: https://www.aviso.altimetry.fr/en/missions/future-missions/swot/orbit.html
