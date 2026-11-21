# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = source
DIRWITHCONFPY = $(if $(wildcard $(dir $(lastword $(MAKEFILE_LIST)))doc-builder),doc-builder,.)
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" -c "$(DIRWITHCONFPY)" $(SPHINXOPTS) $(O)

# 'make fetch-images' should be run before building the documentation. (If building via
# the build_docs command, this is run automatically for you.) This is needed because we
# have configured this repository (via an .lfsconfig file at the top level) to NOT
# automatically fetch any of the large files when cloning / fetching.
fetch-images:
	git lfs install --force
	git lfs pull --exclude="" --include=""

.PHONY: help fetch-images Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" -c "$(DIRWITHCONFPY)" $(SPHINXOPTS) $(O)
