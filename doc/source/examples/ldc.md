# Lid-Driven Cavity

`examples/ldc` is a lid-driven cavity benchmark with prebuilt FOM and ROM scripts and archived outputs.

## FOM

The FOM workflow is driven by `run_fom`, which builds the Nek5000 case, runs the simulation, and stores the snapshots and
log file in `snaps/`.

## ROM

The `run_rom` script builds the ROM executable and runs it against the existing snapshot set. The `.mor` file selects
velocity-only ROM mode with 30 modes and $H^1_0$ POD.

## Notes

This directory is more of a regression/reference setup than a polished tutorial. It is still useful when you want to inspect
how a larger archived case is packaged.
