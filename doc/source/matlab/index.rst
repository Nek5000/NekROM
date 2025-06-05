.. _matlab_section_tag:

.. toctree::
   :maxdepth: 2
   :caption: MATLAB/Octave API:

.. _matlab_driver_section_tag:

MATLAB/Octave Driver
====================

The main driver script is `drive/matlab/rom_online_solver.m`, which loads the basis functions 
and runs the ROM in MATLAB/Octave using a BDF3/EXT3 time-stepper. Supporting
functions are defined in separate files and documented below.

.. _matlab_convection_operators_section_tag:

Plotting
--------

The Matlab/Octave driver uses NekToolKit for plotting of 2D fields. To enable this,
clone the NekToolKit repository and either append it to the `MATLABPATH` environment
variable (for MATLAB) or call `addpath` within `octaverc` (for Octave).

.. code-block:: shell

   git clone https://github.com/kent0/NekToolKit
   export MATLABPATH=$(pwd)/NekToolKit/matlab


Operators
--------------------

.. mat:automodule:: matlab.operators

.. mat:autofunction:: conv_fom

.. mat:autofunction:: conv_deim

.. mat:autofunction:: lgrad

.. mat:autofunction:: lcurl

.. _matlab_operators_section_tag:

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
