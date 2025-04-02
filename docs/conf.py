# Configuration file for the Sphinx documentation builder.
# See the full list of built-in configuration options at:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
project = "OsmoFold"
author = "Vincent Nicholson"
release = "0.4.1"  # The full version, including alpha/beta/rc tags

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",  # Automatically generate documentation from docstrings
    "sphinx.ext.napoleon",  # Support for Google-style and NumPy-style docstrings
    "sphinx.ext.viewcode",  # Add links to source code
]

master_doc = "index"
exclude_patterns = []  # Files and directories to ignore

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"  # Read the Docs theme

# -- Autodoc configuration ---------------------------------------------------
autodoc_mock_imports = []  # List modules to mock if they can't be imported
