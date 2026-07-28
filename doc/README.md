# Building Documentation for NekROM

## Instructions to Build Documentation

### 1. Install Dependencies

#### Set Up a Conda Environment (Optional)

You may use an existing Conda environment, but we recommend creating and activating a new one:
```sh
conda create -n nekrom-docs
conda activate nekrom-docs
conda install pip
```

#### Install Documentation Dependencies
```sh
pip install sphinx sphinx-fortran six sphinx-mathjax-offline sphinx-book-theme myst-parser sphinxcontrib-bibtex sphinxcontrib-matlabdomain "numpy>=2.2.5"
```

If `make html` fails with `Could not import extension sphinxfortran.fortran_domain`, the active environment
is missing the `sphinx-fortran` package. Re-run the install command above inside the same environment that
will execute `make html`.

### 2. Build the NekROM Documentation
Clone NekROM and build the HTML documentation:
```sh
git clone https://github.com/Nek5000/NekROM.git
cd NekROM/doc
make html
```

The generated HTML lives in `NekROM/doc/build/html`.

### 3. View the Documentation
To view the generated documentation, open `build/html/index.html` in your browser:
```sh
cd NekROM/doc/build/html
open index.html
```

## Useful Resources
- [Using Sphinx with Fortran](https://ostueker.github.io/Example_Fortran/UsingSphinxFortran.html)
- [Example Fortran Repository](https://github.com/ostueker/Example_Fortran)
- [Sphinx-Fortran Documentation](https://sphinx-fortran.readthedocs.io/en/latest/user.autodoc.html)
- [VACUMM Sphinx-Fortran](https://github.com/VACUMM/sphinx-fortran/tree/master)

## Refactoring Fortran Code for Documentation
Follow [this guide](https://ostueker.github.io/Example_Fortran/RefactoringFortranForSphinx.html) to refactor source code so that data appears correctly in the documentation.
