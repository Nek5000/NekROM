import os
import shutil

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'NekROM'
copyright = '2025, NekROM contributors'
author = 'NekROM contributors'
release = '1.0.0'

# -- Preprocessing ------------------------------------------------------

# Strip out the include statements for the docs
temp_code_dir = '/tmp/rom_code_for_docs/'
code_dir = '../../code/'
os.makedirs(temp_code_dir, exist_ok=True)
for filename in os.listdir(code_dir):
    f = open(code_dir + filename, mode='rt')
    tmpfile = open(temp_code_dir + filename, mode='wt')
    skip_if_continued = False
    for line in f:
        if 'common ' in line:
            skip_if_continued = True
        elif skip_if_continued == True and len(line) >= 6 and line[5] != ' ':
            pass # The skipped common block is continued onto the next line
        else:
            skip_if_continued = False
            tmpfile.write(line)

    f.close()
    tmpfile.close()

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              'sphinxfortran.fortran_domain',
              'sphinxfortran.fortran_autodoc', 
              'sphinx.ext.mathjax', 
              'sphinx-mathjax-offline']

#fortran_src = [os.path.abspath(code_dir)]
fortran_src = [os.path.abspath(temp_code_dir)]

fortran_ext = ["f"]

templates_path = ['_templates']

exclude_patterns = []

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
