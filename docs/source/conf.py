# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'scAtlas'
copyright = '2024, Gilbert Han'
author = 'Gilbert Han'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'myst_parser',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    "sphinx_markdown_tables",
    "sphinxcontrib.bibtex",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx_copybutton",
    "sphinx_autodoc_typehints",
    "sphinx_tippy",
    "sphinx_design",
    "nbsphinx",
    "IPython.sphinxext.ipython_console_highlighting"
]

exclude_patterns = ['_build', '**.ipynb_checkpoints']
nbsphinx_allow_errors = True
nbsphinx_execute = 'never'
master_doc = "index"
pygments_style = "tango"
pygments_dark_style = "monokai"

nitpicky = True

templates_path = ['_templates']
source_suffix = {
    ".rst": "restructuredtext"
}



# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# bibliography
bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"
bibtex_default_style = "alpha"


# hover
tippy_anchor_parent_selector = "div.content"
tippy_enable_mathjax = True
# no need because of sphinxcontrib-bibtex
tippy_enable_doitips = False
# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_static_path = ["_static"]
html_css_files = [
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.1.1/css/all.min.css",
    "css/override.css",
]

html_show_sphinx = False
html_show_sourcelink = False
html_theme_options = {
    "sidebar_hide_name": True,
    "light_logo": "tooth_log.png", 
    "dark_logo": "tooth_log.png",
    "light_css_variables": {
        "color-brand-primary": "#003262",
        "color-brand-content": "#003262",
        "admonition-font-size": "var(--font-size-normal)",
        "admonition-title-font-size": "var(--font-size-normal)",
    },
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/GilbertHan1011/toothAtlasManuscript",
            "html": "",
            "class": "fab fa-github",
        },
    ],
}