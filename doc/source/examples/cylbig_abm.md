# Cylinder with Augmented Basis Method

`examples/cylbig_abm` is a larger cylinder example configured for the augmented basis method.

## FOM

The case uses `run_fom` to generate the FOM snapshots. The runtime `.mor` file enables both velocity and temperature fields,
sets `augment = 2`, and uses `mode = all` so the offline and online paths are both exercised.

## ROM

The `run_rom` script adjusts `LMOR` for the larger basis size, builds the ROM, and uses `avg.list` alongside `file.list` to
select the snapshot and average inputs.

## Notes

This is the reference case to consult when you want to see NekROM's augmented-basis workflow rather than the minimal cylinder
setup.
