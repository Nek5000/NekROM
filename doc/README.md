# Building Documentation for NekROM

## Instructions to Build Documentation

### 1. Install Dependencies

#### Set Up a Conda Environment (Optional)

If you don’t mind using a development version of NumPy, you can use an existing Conda environment. Otherwise, create and activate a new environment:
```sh
conda create -n nekrom-docs
conda activate nekrom-docs
conda install pip
```

#### Install Documentation Dependencies
```sh
pip install sphinx sphinx-fortran six sphinx-mathjax-offline sphinx-book-theme myst-parser sphinxcontrib-bibtex "numpy>=2.2.5"
```

### 2. Build the NekROM Documentation
Clone NekROM and build the HTML documentation:
```sh
git clone https://github.com/Nek5000/NekROM.git
cd NekROM/doc
make html
```

### 3. View the Documentation
To view the generated documentation, navigate to the build directory and open `index.html` in your browser:
```sh
cd NekROM/doc/build
open index.html
```

## Useful Resources
- [Using Sphinx with Fortran](https://ostueker.github.io/Example_Fortran/UsingSphinxFortran.html)
- [Example Fortran Repository](https://github.com/ostueker/Example_Fortran)
- [Sphinx-Fortran Documentation](https://sphinx-fortran.readthedocs.io/en/latest/user.autodoc.html)
- [VACUMM Sphinx-Fortran](https://github.com/VACUMM/sphinx-fortran/tree/master)

## Refactoring Fortran Code for Documentation
Follow [this guide](https://ostueker.github.io/Example_Fortran/RefactoringFortranForSphinx.html) to refactor source code so that data appears correctly in the documentation.
