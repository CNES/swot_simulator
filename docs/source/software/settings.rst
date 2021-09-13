Sample configuration file
-------------------------

The simulator uses a configuration file using Python code. This paragraph
describes this file and the different parameters you can modify. You will also
find two configuration examples to generate the `L2_LR_SSH_Expert
<https://github.com/CNES/swot_simulator/blob/master/data/l2_lr_expert.py>`_ and
`L2_LR_SSH_Unsmoothed
<https://github.com/CNES/swot_simulator/blob/master/data/l2_lr_unsmoothed.py>`_
SWOT products in official format.

.. note::

    The configuration file can be generated automatically using the
    :ref:`main program <main_program>`, or using the library:

    .. code-block:: python

        import swot_simulator.settings

        settings = swot_simulator.settings.Parameters.load_default()

.. _params-file:

.. literalinclude:: ../settings.py
