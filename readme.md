# NekROM - Model Order Reduction Framework for Nek5000

This package includes tools for applying model-order reduction (MOR) to data produced by [Nek5000](https://github.com/Nek5000/Nek5000) and for generating reduced-order models (ROMs). The generated ROMs can run either in the Fortran driver embedded inside the Nek5000 `userchk` subroutine or through the provided MATLAB and Julia driver scripts. Users can also provide their own drivers that read the ROM operators and quantities of interest (QOI) factors from the `ops/` and `qoi/` directories.

# Documentation

Documentation for NekROM is available at [https://nekrom.readthedocs.io/en/latest/](https://nekrom.readthedocs.io/en/latest/) or can be built locally following the instructions in `doc/README.md`.

# Setup & Procedure

Set shell variables (in .bashrc for BASH users):

```
export MOR_DIR="/path/to/NekROM"
export PATH="$MOR_DIR/bin:$PATH"
```

Required files in NekROM case directory:

- Nek5000 case files e.g., .rea, .map, SIZE
- \$caserom.usr, .usr file specific for NekROM cases (see `$MOR_DIR/examples`)
- LMOR, specifies compile-time parameters
- $case.mor, specifies run-time parameters
- file.list, contains list of paths to the snapshots (relative path)

Optional file:

- avg.list, contains list of paths to the average files

After ensuring the required files are in the case directory, run `makerom $caserom` to make a Nek5000 executable for ROM.

# Parameters

Compile-time parameters (for setting memory allocation size) can be found in `LMOR`.

- `ls`, maximum number of snapshots
- `lb`, maximum number of total modes
- `lbnl`, maximum number of runtime nonlinear POD bases, including the zeroth mode

The authoritative compile-time template is `templates/LMOR.template`.

run-time parameters can be found in `$case.mor`.

- [GENERAL], header for general parameters
    - `mode`, off = offline, on = online, all = offline + online
    - `field`, v = velocity, t = temperature, vt = velocity + temperature
    - `nb`, number of POD modes (must be less than lb, default == lb)
- [POD], header for pod parameters
    - `type`, l2 = $L^2$ POD modes, h10, $H^1_0$ POD modes
    - `mode0`, avg = average 0th mode, state = user-defined in ub,vb,wb,tb
    - `augment`, 0 = no ABM, 1 = 0th interactions, 2 = diagonals, 3 = 1 + 2
- [QOI], header for qoi parameters
    - `freq`, frequency of QOI dump, if <1 freq=iostep
    - `drag`, drag based on OBJ data

Additional runtime options are documented in `templates/mpar.template`, including `avginit`, `rktol`, `nplay`, `combined`, `ratio`, `copt`, `leray`, `tneubc`, `gravity`, `forcing`, `filter`, `ei`, and `deim`.

The shipped example cases under `examples/` provide case-specific setup notes and run scripts. The cylinder example is the shortest end-to-end workflow; `examples/rb_axi` shows the parametric sweep pattern.

## DEIM Stability and Dealiasing

The `deim` runtime option selects sampled DEIM, which is the cheapest online path but can still become unstable on demanding cases. `clsdeim` and `mclsdeim` use constrained or oversampled point selection and are more robust in practice.

For stricter quadrature, the MATLAB driver also supports `NEKROM_DEIM_DEALIAS_QUAD=1`, which forces a 3/2-grid quadrature path for the DEIM-family convection evaluation. That path is stable, but its runtime cost is much closer to a fully dealiased ROM evaluation than to sampled DEIM. It is MATLAB-driver only and is kept in memory rather than being written back into the Fortran-loaded `ops/` bundle.

If you need a cheap dealiased online convection operator today, the tensor-based operators remain the supported option in the current codebase.

The main future direction for a cheaper stable DEIM path is a compressed quadrature layer, such as ECSW-style sampling. That would choose both a reduced set of points on the overintegrated grid and corresponding quadrature weights, so the online path stays closer to sampled DEIM while retaining dealiased integration behavior. This is not implemented yet.

# Contribution

Our procedure for updating the code is exclusively through pull requests (no pushing). Please submit issues and PR to [https://github.com/Nek5000/NekROM](https://github.com/Nek5000/NekROM). PRs should be the smallest coherent change to the code-base. Issue titles should describe the issue, for example `Error in x` or `Missing x`. PR titles should describe the modification made, for example `Fixed x` or `Improved x`. See the documentation for the coding style of this project when contributing.

# Parameter File Support

In addition to `.rea` support for setting internal parameters, `.mor` files are supported as a [par](https://nek5000.github.io/NekDoc/user_files.html)-like dictionary. The possible key/value pairs are described in `templates/mpar.template`.
