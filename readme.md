# NekROM - Model Order Reduction Framework for Nek5000

This package include tools to help apply model-order reduction (MOR) to data produced by [Nek5000](https://github.com/Nek5000/Nek5000) to generate reduced-order models (ROM). The generated ROMs can be run either in the Fortran driver embedded inside the Nek5000 userchk subroutine or can be run separately by provided driver scripts in Matlab, Python, and Julia. Users can also provide their own drivers which read the ROM operators and QOI factors from the `ops/` and `qoi/` directories.

# Setup & Procedure

Set shell variables (in .bashrc for BASH users):

```
export MOR_DIR='/path/to/NekROM'
export PATH="$MOR_DIR/bin:$PATH"
```

Required files in NekROM case directory:

- Nek5000 case files e.g., .rea, .map, SIZE
- $caserom.usr, .usr file specific for NekROM cases (see `$MOR_DIR/examples`)
- LMOR,      specifies compile-time parameters
- $case.mor, specifies run-time parameters
- file.list, contains list of paths to the snapshots (relative path)

Optional file:

- avg.list, contains list of paths to the average files

After ensuring the required files are in the case directory, run `makerom $caserom` to make a Nek5000 executable for ROM.

# Contribution

Our procedure for updating the code is exclusively through pull requests (no pushing). Please submit issues and PR to github.com/SEAL-UIUC/NekMOR. PRs should be the smallest coherent change to the code-base. Issue titles should describe the issue e.g., 'Error in x', 'Missing x', etc. PR titles should describe the modification made e.g., 'Fixed x', 'Improved x, etc. See conventions.txt for the coding style of this project when contributing.

# Parameter File Support

In addition to the .rea support for setting internal parameters, .mor files are supported for [par](https://nek5000.github.io/NekDoc/user_files.html)-like dictionary. The possible key/values are described in templates/mpar.template.
