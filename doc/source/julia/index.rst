.. toctree::
  :maxdepth: 2
  :caption: Julia Driver:

Julia Driver
============

The Julia entry point is `drive/julia/drive.jl`. It reads a `mor.yaml` file from the working directory, loads the ROM
operators from the `ops/` directory, and writes field and coefficient output to `out/`.

Supporting code lives in:

* `drive/julia/functions.jl` for configuration parsing and solver helpers
* `drive/julia/setup.jl` for the package dependencies used by the driver

Minimal workflow
----------------

1. Generate the ROM operator files with the NekROM Fortran workflow.
2. Place a `mor.yaml` file in the Julia driver directory.
3. Run `julia setup.jl` to install the Julia dependencies.
4. Run `julia drive.jl` from the same directory.

The current Julia documentation is intentionally brief and focuses on the shipped driver entry points rather than full API
reference pages.
