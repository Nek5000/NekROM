.. _matlab_section_tag:

.. toctree::
   :maxdepth: 2
   :caption: MATLAB/Octave API:

.. _matlab_driver_section_tag:

MATLAB/Octave Driver
====================

The main driver script is `drive/matlab/driver.m`, which loads the basis functions and runs the ROM in MATLAB/Octave using a
BDF3/EXT3 time-stepper. The older `drive/matlab/old/rom_online_solver.m` file is kept as a reference implementation.

The driver expects `config.m` to select the case and runtime options, and then uses helper functions from the `io/`,
`operators/`, and `point_generators/` subdirectories.

The DEIM-family runtime exposes the `NEKROM_DEIM_FINEGRID`, `NEKROM_DEIM_DEALIAS`, and `NEKROM_DEIM_DEALIAS_QUAD`
environment variables through `config.m`. `deim` is the cheapest sampled path; `clsdeim` and `mclsdeim` improve robustness
with constrained or oversampled point selection; and `NEKROM_DEIM_DEALIAS_QUAD=1` forces strict 3/2-grid quadrature for
the DEIM-family convection evaluation. That strict-quadrature path is stable, but its cost is much closer to a
fully dealiased ROM evaluation than to sampled DEIM. The tensor operators remain the current low-cost dealiased option.

A natural future extension is compressed quadrature or ECSW-style sampling on the overintegrated grid. In that setting, the
offline stage would select both the active quadrature points and their weights, which could keep the online cost lower than
full strict quadrature while preserving dealiased integration behavior. That path is not implemented in the current code.

.. only:: has_matlab_ext

   .. _matlab_convection_operators_section_tag:

   Plotting
   --------

   The MATLAB/Octave driver uses NekToolKit for plotting 2D fields. To enable plotting,
   clone the NekToolKit repository and either append it to the `MATLABPATH` environment
   variable for MATLAB or call `addpath` from `octaverc` for Octave.

   .. code-block:: shell

      git clone https://github.com/kent0/NekToolKit
      export MATLABPATH=$(pwd)/NekToolKit/matlab

   Operators
   --------------------

   .. mat:automodule:: matlab.operators

   .. mat:autofunction:: conv_deim

   .. mat:autofunction:: conv_fom

   .. mat:autofunction:: conv_tensor_dense

   .. mat:autofunction:: conv_tensor

   .. mat:autofunction:: conv_tensor_sparse

   .. mat:autofunction:: gen_Au

   .. mat:autofunction:: get_Me

   .. mat:autofunction:: lgrad

   .. mat:autofunction:: lcurl

   .. _matlab_operators_section_tag:

   Input and Output
   ----------------

   .. mat:automodule:: matlab.io

   .. mat:autofunction:: get_grid

   .. mat:autofunction:: get_pod_basis_from_arrays

   .. mat:autofunction:: get_pod_basis

   .. mat:autofunction:: get_r_dim_ops

   .. mat:autofunction:: get_snaps

   .. mat:autofunction:: get_sort_order

   .. mat:autofunction:: load_full_ops

   .. mat:autofunction:: output_fields

   .. mat:autofunction:: write_field

   .. _matlab_point_generator_section_tag:

   DEIM Point Generators
   ---------------------

   .. mat:automodule:: matlab.point_generators

   .. mat:autofunction:: s_opt

   .. mat:autofunction:: gnat

   .. mat:autofunction:: gappy_pod

   The legacy `gpode` implementation lives under `drive/matlab/point_generators/old/gpode/` and is not part of the current
   auto-documented MATLAB API.

.. only:: not has_matlab_ext

   The MATLAB domain extension is not installed in this build, so the API reference is omitted.
