# 3D PERIODIC HILL (NekROM Example)

This example is similar to the one presented in Nek5000 documentation, however
the hill profile is consistent with Almeida et al. 1993 and described at
https://turbmodels.larc.nasa.gov/Other_LES_Data/2Dhill_periodic/hill-geometry.dat

## How To Run

This directory contains both a full-order model (FOM) run script to
generate snapshots, and a ROM run script that uses `phill.mor`.

1. Generate snapshots (FOM): `./run_fom <np>`
2. Build ROM operators / run ROM: `./run_rom <np>`

Notes:
- This is a 3D case sized for MPI. `SIZE` sets `lpmin=12`, and the provided
  scripts default to `np=lpmin`.
- `run_fom` uses `genbox` (`Nek5000/tools/bin/genbox`) to create `phill.re2`
  from `phill.box` when `phill.re2` is missing.
- `run_fom` writes snapshots to `snaps/` and records them in `file.list`.
- `phill.par` controls time stepping and output cadence. Ensure `phill.mor`
  `ns` is <= the number of snapshots in `file.list`.

REMARKS:
- phill.re2 is a rectangular channel
- the geometric parameters of the hill are defined in userdat2
- it is assumed that the Reynolds number is based on the average velocity
  over the inlet plane
