# Axisymmetric Rayleigh-Benard

`examples/rb_axi` demonstrates the parametric workflow used by NekROM for axisymmetric Rayleigh-Benard convection.

## FOM

The example runs five FOM snapshot sets over a parameter interval from `eps = 1.6` to `eps = 2.6`. Each case is launched from
`run_fom`, which writes the `eps` control file, runs the helper script, and stores results in `fom0` through `fom4`.

## ROM

The `run_rom` script first builds ROM bases from the endpoint snapshot sets, then combines them into an interpolation basis.
It then sweeps 21 parameter points across the same interval and runs the ROM for each point.

## Post-Processing

The bundled `plot_nus.jl` script compares the FOM and ROM Nusselt number histories and writes `nus.pdf`.

## Notes

This case is the best reference for parametric ROM setups in NekROM.
