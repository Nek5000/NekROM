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
```

### 2. Install Dependencies

#### Install Dependencies for Building NumPy
```sh
conda install pip cython compilers openblas meson-python pkg-config
```

#### Install Documentation Dependencies
```sh
pip install sphinx sphinx-fortran six
```

### 3. Install NumPy from a Specific Branch
Clone and install NumPy from the `common-block-division` branch:
```sh
git clone https://github.com/nchristensen/numpy.git
cd numpy
git checkout common-block-division
git submodule update --init
pip install . --no-build-isolation
```
Refer to [NumPy's build instructions](https://numpy.org/doc/stable//building/index.html#building-from-source-to-use-numpy) for more details.

### 4. Build the NekROM Documentation
Clone the `documentation` branch of NekROM and build the HTML documentation:
```sh
git clone https://github.com/nchristensen/NekROM.git
cd NekROM/doc
make html
```

### 5. View the Documentation
To view the generated documentation, navigate to the build directory and open `index.html` in your browser:
```sh
cd NekROM/doc/build
open index.html
```

## Refactoring Fortran Code for Documentation
Follow [this guide](https://ostueker.github.io/Example_Fortran/RefactoringFortranForSphinx.html) to refactor source code so that data appears correctly in the documentation.

### Known Issue
- `sphinx-fortran` does not support multi-line comments for parameters.
