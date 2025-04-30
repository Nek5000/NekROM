# Building Documentation for NekROM

## Useful Resources
- [Using Sphinx with Fortran](https://ostueker.github.io/Example_Fortran/UsingSphinxFortran.html)
- [Example Fortran Repository](https://github.com/ostueker/Example_Fortran)
- [Sphinx-Fortran Documentation](https://sphinx-fortran.readthedocs.io/en/latest/user.autodoc.html)
- [VACUMM Sphinx-Fortran](https://github.com/VACUMM/sphinx-fortran/tree/master)

## Instructions to Build Documentation

### 1. Set Up a Conda Environment
If you don’t mind using a development version of NumPy, you can use an existing Conda environment. Otherwise, create and activate a new environment:
```sh
conda create -n nekrom-docs
conda activate nekrom-docs
conda install pip
```

### 2. Install Dependencies

#### Install Documentation Dependencies
```sh
pip install sphinx sphinx-fortran six sphinx-mathjax-offline sphinx-book-theme myst-parser "numpy>=2.2.5"
```

#### Install Development Version of NumPy
Note: This is no longer required after the release of NumPy version 2.2.5.

##### Install Dependencies for Building NumPy
```sh
conda install cython compilers openblas meson-python pkg-config
```

Clone and build development version of NumPy:
```sh
git clone https://github.com/numpy/numpy.git
cd numpy
git submodule update --init
pip install . --no-build-isolation
```
Refer to [NumPy's build instructions](https://numpy.org/doc/stable//building/index.html#building-from-source-to-use-numpy) for more details.

### 3. Build the NekROM Documentation
Clone NekROM and build the HTML documentation:
```sh
git clone https://github.com/Nek5000/NekROM.git
cd NekROM/doc
make html
```

### 4. View the Documentation
To view the generated documentation, navigate to the build directory and open `index.html` in your browser:
```sh
cd NekROM/doc/build
open index.html
```

## Refactoring Fortran Code for Documentation
Follow [this guide](https://ostueker.github.io/Example_Fortran/RefactoringFortranForSphinx.html) to refactor source code so that data appears correctly in the documentation.
