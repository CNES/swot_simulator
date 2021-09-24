Changelog
=========

1.2.2 - 2021/09/24
------------------
* Logging settings
* Plugins to interpolate Schism maps
* Fixes issues on Nadir products

1.2.1 - 2021/09/15
------------------
* DATA folder is not available after pip install.

1.2 - 2021/09/13
----------------
* Bug fix: swaths are stored from right to left instead of left to right.
* Writing of Nadir products in proper format (official GDR reduced format).
* Addition of the orbital error.
* The simulator version is written in the generated products.

1.1 - 2021/07/21
----------------
* Support for unsmoothed products.
* Retains the simulation parameters.
* Bug fix: cross track distances are written in km.
* Set the satellite coordinates at the equator.
* Refactor the plugin architecture.

1.0 - 2021/03/12
----------------
* Restructuring of the random number generation.
* The long wavelengths are kept to simulate the low frequency signal.
* Added a method to get default settings.
* The signal is insufficient randomness at large wavelengths.
* Bug fix: The gyroscope effect is not considered.
* Adding the documentation.
* Added generation tests under Windows

0.3 - 2021/02/10
----------------

* Simulated errors are renamed at the request of JPL.
* Add SWH interpolation plugin and compute corresponding KaRIN
* Fix bug in KarIN
* Fix bug in naming of variables

0.2 - 2020/10/16
----------------

* Update product descriptions
* Add HYCOM model plugin

0.1 - 2020/05/29
----------------

* First beta release
