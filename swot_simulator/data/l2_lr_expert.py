# Configuration file used to generate the SWOT_L2_LR_SSH_Expert files according
# to the SWOT PDD.
import pathlib
import swot_simulator

# Geographical area to simulate defined by the minimum and maximum corner
# point :lon_min, lat_min, lon_max, lat_max. Default: -180, -90, 180, 90
area = None

# Number of beam used to correct wet troposphere signal (1, 2 or 'both')
beam_position = [-20, 20]

# If true, the swath, in the final dataset, will contain a center pixel
# divided in half by the reference ground track
central_pixel = True

# If true, the generated netCDF file will be the complete product compliant
# with SWOT's Product Description Document (PDD), otherwise only the calculated
# variables will be written to the netCDF file
complete_product = True

# Duration of a cycle in number of fractional days. By default, this value
# is read from the ephemeris file)
# cycle_duration = None

# Distance, in km, between two points across track direction
delta_ac = 2.0

# Distance, in km, between two points along track direction
delta_al = 2.0

# Index of columns to read in the ephemeris file containing, respectively,
# longitude in degrees, latitude in degrees and the number of seconds elapsed
# since the start of the orbit. Default: [1, 2, 0]
ephemeris_cols = None

# Ephemeris file to read containing the satellite's orbit.
#ephemeris = os.path.join(DATA, 'ephemeris_calval_june2015_ell.txt')
ephemeris = swot_simulator.DATA / 'ephemeris_science_sept2015_ell.txt'

# File containing spectrum of instrument error
error_spectrum = swot_simulator.DATA / 'error_spectrum.nc'

# Estimated roll phase dataset. Default: None
corrected_roll_phase_dataset = None

# Distance, in km, between the nadir and the center of the first pixel
# of the swath
half_gap = 2.0

# Distance, in km, between the nadir and the center of the last pixel
# of the swath
half_swath = 70.0

# Satellite altitude (m). By default, this value is read from the ephemeris
# file.
# height = None

# KaRIN file containing spectrum for several SWH
karin_noise = swot_simulator.DATA / 'karin_noise_v2.nc'

# Repeat length
len_repeat = 20000.0

# True to generate Nadir products
nadir = True

# Number of beam used to correct wet troposphere signal (1, 2 or 'both')
nbeam = 2

# This option defined the error to simulate. Allowed values are "altimeter,"
# "baseline_dilation," "corrected_roll_phase," "karin," "orbital," "roll_phase,"
# "timing," and "wet_troposphere." The calculation of roll errors can be
# simulated, option "roll_phase," or interpolated option "corrected_roll_phase,"
# from the dataset specified by the value of the option "roll_phase_dataset."
# Therefore, these two options are mutually exclusive. In other words, if
# the "roll_phase" option is present, the "corrected_roll_phase" option must
# be omitted, and vice versa.
noise = [
    'altimeter',
    'baseline_dilation',
    'karin',
    'orbital',
    'roll_phase',
    'timing',
    # 'wet_troposphere',
]

# Seed for RandomState. Must be convertible to 32 bit unsigned integers
nseed = 0

# Type of SWOT product to be generated. Possible products are 'basic', 'expert',
# 'unsmoothed' and 'wind_wave'. Default to 'expert'
product_type = 'expert'

# Limits of SWOT swath requirements. Measurements outside the span will be
# set with fill values
requirement_bounds = [10, 60]

# Orbit shift in longitude (degrees). Default 0.
shift_lon = None

# Orbit shift in time (seconds). Default 0.
shift_time = None

# Gaussian footprint of sigma (km)
sigma = 6.0

# The plug-in handling the SSH interpolation under the satellite swath
ssh_plugin = None

# SWH plugin to interpolate model SWH on the SWOT grid. Use only 'expert'
# or 'wind_wave' products that have the swh in its output
swh_plugin = None

# True to generate swath products
swath = True

# SWH for the region
swh = 2.0

# The working directory. By default, files are generated in the user's root
# directory
working_directory = None
