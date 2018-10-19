# Model Order Reduction (MOR)

[![Build Status](https://travis-ci.com/kent0/MOR.svg?token=nDCiae81x8NojggcMEcA&branch=master)](https://travis-ci.com/kent0/MOR)

To download the baffle case snapshots, goto MOR/bin and run `./gsnaps`.

For each case run `../../bin/linkc` in a case directory to link the source.

## Code

* rom.f - includes ROM subroutines
* pod.f - include POD subroutines
* aux.f - includes auxiliary subroutines

## Cases

### Baffle
### Cylinder

### Lid-Driven Cavity

A 2D case used to reproduce [1].

Domain: x,y \in [-1,1]
BCs: u(x,y) = \delta (y-1) * (1-x^2)^2

The top boundary condition is not the standard uniform velocity for consistency with [1]

[1]: Fick, Maday, Patera, Taddei; "A Reduced Basis Technique for Long-Time Unsteady Turbulent Flows"

