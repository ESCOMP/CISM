"""
Substitutions for Sphinx
"""

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.

# pylint: disable=invalid-name

#################################
### Standard Sphinx variables ###
#################################

# General information about the project.
project = "cism"
copyright = "2026, UCAR"  # pylint: disable=redefined-builtin
author = ""

# The short X.Y version.
version = "CISM3"

# The full version, including alpha/beta/rc tags.
release = "CISM main"

#####################################################
### Custom variables needed for doc-builder setup ###
#####################################################

# Version label used at the top of some pages.
version_label = "the latest development code"

#######################################################
### Custom variables optional for doc-builder setup ###
#######################################################

tex_category = "Miscellaneous"

# Used by HTML help builder
htmlhelp = {
    "basename": "cismdocdoc", # Output file base name
}

# Used for LaTeX output
latex = {
    "target_name": "cismdoc.tex",
    "title": "CISM Documentation",
    "documentclass": "manual", # howto, manual, or own class
    "category": tex_category,
}

# Used for man_pages and texinfo_documents
mantex = {
    "name": "cismdoc",
    "title": "cismdoc Documentation",
}

# Used for texinfo_documents
tex = {
    "dirmenu_entry": "cismdoc",
    "description": "One line description of project.",
    "category": tex_category,
}

###############################
### Purely custom variables ###
###############################

nonparamfile_disclaimer_md = (
    "**Note:** The values here should be up-to-date with those used in {{version_label}},"
    " but there may be mistakes."
)
nonparamfile_disclaimer_rst = (
    "**Note:** The values here should be up-to-date with those used in |version_label|,"
    " but there may be mistakes."
)