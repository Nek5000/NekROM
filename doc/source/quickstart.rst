.. _quickstart_section_tag:
  
Quickstart
==========

Running the flow past a cylinder example

In its most straightforward workflow, using NekROM requires running a full-order model simulation to generate solution snapshots
and then executing NekROM targeting these snapshots to generate POD bases and run the reduced-order simulation. Here we illustrate
the process with the 2D flow past a cylinder example.

1. Clone NekROM
2. Add NekROM/bin and Nek5000/bin to your PATH and set the MOR_DIR environment variable to your NekROM directory.
3. cd to `NekROM/examples/cyl`
4. Run `./run_fom` to generate the FOM snapshots (and drag/lift data).
5. Run `./run_rom` to created the POD bases and simulate with the reduced-order model.
6. Run `visnek romcyl` and open `romcyl.nek5000` with Visit or Paraview to visualize the ROM output

The `run_fom` and `run_rom` scripts automate several tasks in the NekROM workflow. We explain these scripts in more detail below. TODO

Running the lid-driven cavity example (Might be the same process, worth adding?)
ADD ME

What are all of these files NekROM produces? (Rename)
ADD ME
