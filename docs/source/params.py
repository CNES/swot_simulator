import os
import swot_simulator.plugins.ssh
import swot_simulator.plugins.swh

# Geographical area to simulate defined by the minimum and maximum corner
# point :lon_min, lat_min, lon_max, lat_max
#
# Default: None equivalent to the area covering the Earth: -180, -90, 180, 90
area = [0, -60, 20, -50]

# Distance, in km, between two points along track direction.
delta_al = 2.0

# Distance, in km, between two points across track direction.
delta_ac = 2.0

# Ephemeris file to read containing the satellite's orbit.
ephemeris = os.path.join("data_swot", "ephem_science_sept2015_ell.txt")

# Index of columns to read in the ephemeris file containing, respectively,
# longitude in degrees, latitude in degrees and the number of seconds elapsed
# since the start of the orbit.
# Default: 1, 2, 0
ephemeris_cols = [1, 2, 0]

# If true, the generated netCDF file will be the complete product compliant
# with SWOT's Product Description Document (PDD), otherwise only the calculated
# variables will be written to the netCDF file.
complete_product = True

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
nadir = True

# True to generate swath products
swath = True

# The plug-in handling the SSH interpolation under the satellite swath.
ssh_plugin = swot_simulator.plugins.ssh.AVISO(
    "/mnt/data/data_model/cmems_nrt_008_046")

# Orbit shift in longitude (degrees)
# #shift_lon=

# Orbit shift in time (seconds)
# #shift_time=

# The working directory. By default, files are generated in the user's root
# directory.
working_directory = '/mnt/data/test_swot'

# Generation of measurement noise.

#############################
##   ERROR PARAMETERS      ##
#############################
# The calculation of roll errors can be simulated, option "roll_phase", or
# interpolated, option "corrected_roll_phase", from the dataset specified by
# the option "roll_phase_dataset". Therefore, these two options are
# mutually exclusive. In other words, if the "roll_phase" option is present,
# the "corrected_roll_phase" option must be omitted, and vice versa.
noise = [
    'altimeter',
    'baseline_dilation',
    'karin',
    'corrected_roll_phase',
    # 'roll_phase',
    'timing',
    'wet_troposphere',
]

# repeat length
len_repeat = 20000

# ---- Karin noise
# KaRIN file containing spectrum for several SWH:
karin_noise = os.path.join('data_swot', 'karin_noise_v2.nc')

# SWH for the region:
#        if swh greater than 7 m, swh is set to 7
swh = 2

# SWH plugin to interpolate model SWH on the SWOT grid:
swh_plugin = swot_simulator.plugins.ssh.WW3("/mnt/data/data_model/ww3")

# Number of km of random coefficients for KaRIN noise (recommended nrandkarin=1000):
nrand_karin = 1000

# --- Other instrumental and atmospheric noise
# File containing spectrum of instrument error:
error_spectrum = os.path.join('data_swot', 'global_sim_instrument_error.nc')

# Seed for RandomState. Must be convertible to 32 bit unsigned integers.
nseed = 0

# Roll-phase simulation of correction file
corrected_roll_phase_dataset = os.path.join('data_swot',
                                            'data_sim_slope_2cycles_v0.nc')

# Beam print size (in km):
#        Gaussian footprint of sigma km
sigma = 8.

# Number of beam used to correct wet_tropo signal (1, 2 or 12 for both):
nbeam = 2

# ------ Beam position if there are 2 beams (in km from nadir):
beam_position = [-35, 35]
