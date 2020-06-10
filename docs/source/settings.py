"""
Simulation settings
-------------------
"""
import os
import swot_simulator.plugins.ssh

# Geographical area to simulate defined by the minimum and maximum corner
# point :lon_min, lat_min, lon_max, lat_max
#
# Default: None equivalent to the area covering the Earth: -180, -90, 180, 90
area = None

# Distance, in km, between two points along track direction.
delta_al = 2.0

# Distance, in km, between two points across track direction.
delta_ac = 2.0

# Ephemeris file to read containing the satellite's orbit.
ephemeris = os.path.join("..", "..", "data", "ephem_science_sept2015_ell.txt")

# Index of columns to read in the ephemeris file containing, respectively,
# longitude in degrees, latitude in degrees and the number of seconds elapsed
# since the start of the orbit.
# Default: 1, 2, 0
ephemeris_cols = [1, 2, 0]

# If true, the swath, in the final dataset, will contain a center pixel
# divided in half by the reference ground track.
central_pixel = True

# If true, the generated netCDF file will be the complete product compliant
# with SWOT's Product Description Document (PDD), otherwise only the calculated
# variables will be written to the netCDF file.
complete_product = False

# Distance, in km, between the nadir and the center of the first pixel of the
# swath
half_gap = 2.0

# Distance, in km, between the nadir and the center of the last pixel of the
# swath
half_swath = 70.0

# Limits of SWOT swath requirements. Measurements outside the span will be set
# with fill values.
requirement_bounds = [10, 60]

# The next two parameters (cycle_duration and height) can be read from the
# ephemeris file if it includes these values in comments. The ephemeris
# delivered with this software contain this type of declaration

# Duration of a cycle.
# #cycle_duration=

# Satellite altitude (m)
# #height=

# True to generate Nadir products
nadir = False

# True to generate swath products
swath = True

# Type of SWOT product to be generated. Possible products are "expert",
# "basic" and "wind_wave". Default to expert
product_type = "expert"

# The plug-in handling the SSH interpolation under the satellite swath.
ssh_plugin = swot_simulator.plugins.ssh.AVISO("PATH to AVISO files")

# Orbit shift in longitude (degrees)
# #shift_lon=

# Orbit shift in time (seconds)
# #shift_time=

# The working directory. By default, files are generated in the user's root
# directory.
# #working_directory=

# Generation of measurement noise.

# The calculation of roll errors can be simulated, option "roll_phase", or
# interpolated, option "corrected_roll_phase", from the dataset specified by
# the option "roll_phase_dataset". Therefore, these two options are
# mutually exclusive. In other words, if the "roll_phase" option is present,
# the "corrected_roll_phase" option must be omitted, and vice versa.
noise = [
    'altimeter',
    'baseline_dilation',
    'karin',
    # 'corrected_roll_phase',
    'roll_phase',
    'timing',
    'wet_troposphere',
]

# repeat length
len_repeat = 20000

# File containing spectrum of instrument error
error_spectrum = os.path.join("..", "..", "data", "error_spectrum.nc")

# KaRIN file containing spectrum for several SWH
karin_noise = os.path.join("..", "..", "data", "karin_noise_v2.nc")

# Estimated roll phase dataset
# #corrected_roll_phase_dataset =

# SWH for the region
swh = 2.0

#  Number of km of random coefficients for KaRIN noise (recommended
# nrand_karin=1000)
nrand_karin = 1000

# Number of beam used to correct wet troposphere signal (1, 2 or 'both')
nbeam = 2

# Gaussian footprint of sigma km
sigma = 6.0

# Beam position if there are 2 beams (in km from nadir):
beam_position = [-20, 20],

# Seed for RandomState. Must be convertible to 32 bit unsigned integers.
nseed = 0
