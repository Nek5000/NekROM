# Model Order Reduction (MOR)

[![Build Status](https://travis-ci.com/kent0/MOR.svg?token=nDCiae81x8NojggcMEcA&branch=master)](https://travis-ci.com/kent0/MOR)

To download the baffle case snapshots, go to MOR/bin and run `./gsnaps baf`.
To download the cyl case snapshots, go to MOR/bin and run `./gsnaps cyl`.

For each case run `../../bin/linkc` in a case directory to link the source.

## Code

* rom.f - includes ROM subroutines
* pod.f - include POD subroutines
* aux.f - includes auxiliary subroutines
* dump.f - includes routines that write out ROM states
* read.f - includes routines that read in ROM states
* part.f - includes routines used for partitioning scheme

## Cases

### Baffle

- [1000 snapshots](https://uofi.box.com/shared/static/ktkxit8tblr6mrngwkwb9rhai6rwnqrc.gz)
- [100 snapshots](https://uofi.box.com/shared/static/bnoflcpf6i6mmofiy4ufmm8gey6t83k1.gz)

### Cylinder

- [snapshots](https://uofi.box.com/shared/static/insxch1bjvzl4bnzr3lpr1c016i4g75p.gz)

### Lid-Driven Cavity

A 2D case used to reproduce [1].

Domain: x,y \in [-1,1]
BCs: u(x,y) = \delta (y-1) * (1-x^2)^2

The top boundary condition is not the standard uniform velocity for consistency with [1]

- [500 snapshots](https://uofi.box.com/shared/static/f7h48n1jnuc7qx4cbej009m50a5d3skj.gz)

[1]: Fick, Maday, Patera, Taddei; "A Reduced Basis Technique for Long-Time Unsteady Turbulent Flows"
