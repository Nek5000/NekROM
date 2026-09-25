# Annulus Convection

`examples/ann` is a buoyancy-driven convection case in an annular geometry. It is useful as a compact reference for coupled
velocity and temperature ROM setups.

## FOM

The case uses `ann.mor` with `mode = all`, `field = vt`, and a very small POD basis (`nb = 1`). It also enables buoyancy and
uses the `Grashof` file to parameterize the run-time forcing.

## ROM

The `ann_rom.usr` file mirrors the FOM user file but adds the ROM entry points and `rom_update` hook. The case is a useful
reference for coupled buoyancy and thermal ROM logic.

## Notes

This example is the shortest shipped case for looking at a temperature-coupled buoyancy problem rather than a pure fluid ROM.
