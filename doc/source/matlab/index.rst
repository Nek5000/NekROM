.. _matlab_section_tag:

.. toctree::
   :maxdepth: 2
   :caption: MATLAB/Octave API:

.. _matlab_driver_section_tag:

MATLAB Driver
=============

The main driver script is `rom_online_solver.m.`, which loads the basis functions 
and runs the ROM in MATLAB/Octave using a BDF3/EXT3 time-stepper. Supporting
functions are defined in separate files and documented below.

.. _matlab_convection_operators_section_tag:

Convection operators
--------------------

.. mat:automodule:: matlab

.. mat:autofunction:: conv_fom

.. mat:autofunction:: conv_deim


.. _matlab_io_section_tag:

Input and Output
----------------

.. mat:automodule:: matlab.io

.. mat:autofunction:: load_full_ops

.. mat:autofunction:: get_r_dim_ops

.. mat:autofunction:: get_grid

.. mat:autofunction:: get_snaps


.. _matlab_point_generator_section_tag:

DEIM Point Generators
---------------------

.. mat:automodule:: matlab.point_generators

.. mat:autofunction:: s_opt

.. mat:autofunction:: gpode

.. mat:autofunction:: gnat

.. mat:autofunction:: gappy_pod

