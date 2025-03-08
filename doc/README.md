Useful resources
https://ostueker.github.io/Example_Fortran/UsingSphinxFortran.html
https://github.com/ostueker/Example_Fortran
https://sphinx-fortran.readthedocs.io/en/latest/user.autodoc.html
https://github.com/VACUMM/sphinx-fortran/tree/master

Instructions to build documentation:
1. Set up a new conda environment and activate it or use an existing conda environment if you don't mind using a development version of numpy
2. Clone this branch of numpy: https://github.com/nchristensen/numpy/tree/common-block-division
3. cd into the numpy directory and install it into the conda environment using these instructions: https://numpy.org/doc/stable//building/index.html#building-from-source-to-use-numpy
4. run pip install sphinx and run pip install sphinx-fortran
5. Clone this branch of NekROM: https://github.com/nchristensen/NekROM/tree/documentation
6. cd into NekROM/doc
7. run make html
8. To view the documentation, cd into NekROM/doc/build and open index.html in your web browser
9. Refactor source code following this style to make data show up in the documentation. 
