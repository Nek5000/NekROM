# Double Shear-Layer Roll-Up

`examples/shear` documents a double shear-layer roll-up problem with passive scalar transport.

## Problem Setup

The case includes both thin and thick shear-layer initial conditions. The thin case is more demanding and is the better
stress test for stability. The example also demonstrates transport of multiple passive scalars.

## Numerical Notes

The case README emphasizes that aliasing in the nonlinear quadrature is a major stability issue. The supplied input uses
dealiasing and filtering choices that avoid blow-up at high Reynolds numbers.

This case is also a useful stress test for DEIM-family stabilization: sampled DEIM can still become unstable on long
runs, `clsdeim` and `mclsdeim` are more robust, and strict-quadrature DEIM is stable but much more expensive.

## Reference Use

This is a useful example when you need to see how NekROM handles passive scalars and more demanding stabilization settings.
