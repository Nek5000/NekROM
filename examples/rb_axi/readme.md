# Axisymmetric Rayleigh-Benard Problem

Studied by Tuckerman & Barkley [1,2], the Rayleigh-Benard problem solved in an axisymmetric configuration exhibits traveling wave solutions above a critical Rayleigh Number.

## FOM

The domain is unit-height and has radius of 5 with Prandtl number 10 (consistent with [1]). 5 FOM snapshots sets are produced equally spaced from $\epsilon=\{1.6,\dots,2.6\}$ by running the script `./run_fom` in directories fom\*.

## ROM

The `./run_rom` script generates a ROM by first performing PODs on the the first and last snapshot sets, then combine the POD basis to generate a ROM. The ROM is run using 21 equally spaced points in the same parameter span as the FOM.

## TDEIM note

This case includes a coupled velocity+temperature ROM. NekROM can optionally load thermal DEIM (“TDEIM”) operators from `ops/tdeim_*` to accelerate the temperature convection term. For `rb_axi`, enabling TDEIM can destabilize some parameter points unless additional stabilization is used (e.g. `decoupled = yes` in `rb.mor` to clamp coefficients to the training ranges).

For robustness, `./run_rom` defaults to disabling TDEIM online via `MOR_DISABLE_TDEIM=1`. To experiment with TDEIM, set `MOR_DISABLE_TDEIM=0` and consider enabling `decoupled = yes`.

## Post-Processing

The `./plot_nus.jl` Julia script produced a nus.pdf file that compares the FOM- and ROM-generated Nusselt number data which establishes that the ROM captures FOM behavior efficiently.


[1] L. Tuckerman and D. Barkley, “Global bifurcation to traveling waves in axisymmetric convection,” Phys. Rev. Lett., vol. 61, pp. 408–411, 4 1988.
[2] D. Barkley and L. Tuckerman, “Traveling waves in axisymmetric convection: The role of sidewall conductivity,” Physica D., vol. 37, pp. 288–294, 1989.
