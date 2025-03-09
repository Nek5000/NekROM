Useful resources
https://ostueker.github.io/Example_Fortran/UsingSphinxFortran.html
https://github.com/ostueker/Example_Fortran
https://sphinx-fortran.readthedocs.io/en/latest/user.autodoc.html
https://github.com/VACUMM/sphinx-fortran/tree/master

Instructions to build documentation:

# Set up a new conda environment and activate it or use an existing conda environment if you don't mind using a development version of numpy
conda create -n nekrom-docs
conda activate nekrom-docs
# Install dependencies to build numpy
conda install pip cython compilers openblas meson-python pkg-config

# Install documentation dependencies with pip
pip install sphinx 
pip install sphinx-fortran 
pip install six

# Clone this branch of numpy: https://github.com/nchristensen/numpy/tree/common-block-division
# cd into the numpy directory and install it into the conda environment as follow (uses these instructions: https://numpy.org/doc/stable//building/index.html#building-from-source-to-use-numpy)
cd numpy
git submodule update --init
pip install . --no-build-isolation

# Clone this branch of NekROM: https://github.com/nchristensen/NekROM/tree/documentation
cd NekROM/doc
make html
# To view the documentation, cd into NekROM/doc/build and open index.html in your web browser
# Refactor source code following this style to make data show up in the documentation. https://ostueker.github.io/Example_Fortran/RefactoringFortranForSphinx.html
