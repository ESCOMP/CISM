$$$$$$$$$$$$$$$$$$$
This instance of doc-builder was copied from https://github.com/ESMCI/doc-builder/tree/v3.4
$$$$$$$$$$$$$$$$$$$

This tool wraps the build command to build sphinx-based documentation.

This tool assists with creating the correct documentation build commands
in cases including:
- Building the documentation from a container
- Building versioned documentation, where the documentation builds land
  in subdirectories named based on the source branch

This tool should be put somewhere in your path. Then it should be run
from the directory that contains the Makefile for building the
documentation.

Simple usage is:

    build_docs -b /path/to/doc/build/repo/some/subdirectory [-c] [-d]

    Common additional flags are:
    -c: Before building, run 'make clean'
    -d: Use a container to build the documentation

Usage for automatically determining the subdirectory in which to build,
based on the version indicated by the current branch, is:

    ./build_docs -r /path/to/doc/build/repo [-v DOC_VERSION]

    This will build the documentation in a subdirectory of the doc build
    repo, where the subdirectory is built from DOC_VERSION. If
    DOC_VERSION isn't given, it will be determined based on the git
    branch name in the doc source repository.

    In the above example, documentation will be built in:
    /path/to/doc/build/repo/versions/DOC_VERSION

    This usage also accepts the optional arguments described above.
