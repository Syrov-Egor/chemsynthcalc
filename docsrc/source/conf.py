# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'chemsynthcalc'
copyright = '2022, Egor Syrov'
author = 'Egor Syrov'
release = '1.0.8'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
import sys
import os
import sphinx_theme
import sphinx_rtd_theme

#sys.path.insert(0, os.path.abspath("src"))
#

#print(chemsynthcalc.chem_errors.BadCoeffiecients)


#sys.path.append(path)
#import chemsynthcalc.chem_errors

    # If extensions (or modules to document with autodoc) are in another directory,
    # add these directories to sys.path here. If the directory is relative to the
    # documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
dirname=os.path.dirname
path = dirname(dirname(dirname(os.path.realpath(__file__))))+"\src"
sys.path.append(path)

extensions = ['sphinx.ext.autodoc', 
'sphinx.ext.autosummary', 
'sphinx.ext.napoleon', 
'sphinx_rtd_theme',
'sphinx_search.extension']
autodoc_member_order = 'bysource'
# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

master_doc = 'index'

templates_path = ['_templates']
exclude_patterns = []

language = 'English'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
#html_static_path = ['_static']