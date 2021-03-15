Sampled SSH and error fields
----------------------------

At each pass, for each cycle, an output netCDF file containing the SSH
interpolated from the model (if :py:const:`ssh_plugin <settings.ssh_plugin>` is
not set to `None`), and the different errors. The naming and the file format are
compliant with SWOT's PDD (Product Description Document), set
:py:const:`complete_product <settings.complete_product>` key to ``True`` to
store all the variables and ``False`` to store only the one compute by the
simulator.

You can also select the output type
:py:const:`product_type <settings.product_type>` among ``expert``, ``ssh``, and
``windwave``. Default is ``expert`` and details are in the corresponding xml
file.

The output directory is defined in
:py:const:`working_directory <settings.working_directory>` key, default is the
user's root. The output file names are stored in the output directory in
``karin/<year>`` directory for SWOT and ``nadir/<year>`` directory for the
nadir. The naming follows the pattern :

* ``SWOT_L2_LR_Expert_[cycle]_[pass]_[start_time]_[stop_time]_DG10_01.nc``
  for the Karin products and
* ``SWOT_GPN_2P1P_[cycle]_[pass]_[start_time]_[stop_time].nc``
  for the nadir products.

The interpolation method is specified in your SSH plugin. `pangeo-pyinterp
<https://github.com/CNES/pangeo-pyinterp>`_ module is used in the examples.

The :py:const:`nadir <settings.nadir>` and :py:const:`swath <settings.swath>`
parameters enable (parameter set to ``True``) or disable (parameter set to
``False``) the generation of these products.

Computation of :py:const:`errors <settings.noise>` are specified as a list:

.. code-block:: python

  noise = ['altimeter',
           'baseline_dilation',
           'karin',
           'corrected_roll_phase',
           'timing',
           'wet_troposphere']

.. note::

    ``corrected_roll_phase`` error will generate already cross-calibrated
    roll-phase error whereas ``roll_phase`` will generate errors before
    cross-calibration.

A :py:const:`repeat length <settings.len_repeat>` ``len_repeat = 2000`` key
defines the wavelength in km to repeat the noise and is used for all noises. The
path to the file that contains instrumental errors is mentioned in
:py:const:`error_spectrum <settings.error_spectrum>` key. It is possible to
specify the seed for Randomstate in :py:const:`nseed <settings.nseed>`. The
following parameters are specific to each noise component:

* KaRIN: The file noise that contains the spectrum is available in
  :py:const:`karin_noise <settings.karin_noise>`, a constant
  :py:const:`swh <settings.swh>` can be set ``swh=2`` or interpolated from a
  model SWH specified in the plugin
  :py:const:`swh_plugin <settings.swh_plugin>`.
* Roll-phase: To use the already cross-calibrated roll phase, download the file
  from the FTP (https://ftp.odl.bzh/swot) and specify the path in the
  :py:const:`corrected_roll_phase_dataset <settings.corrected_roll_phase_dataset>`
  key. So far, the following files are available:

  * ``data_sim_slope_v0.nc``: One year of cross-calibrated roll-phase
  * ``data_sim_slope_2cycles_v0.nc``: Two cycles of cross-calibrated
    roll-phase

* Wet troposphere: The number of beams used to correct the wet troposphere
  is set in :py:const:`nbeam <settings.nbeam>` variable. The beam position of
  each beam can be set as a list in
  :py:const:`beam_position <settings.beam_position>`, and the gaussian
  footprint can be changed using :py:const:`sigma <settings.sigma>` variable.
