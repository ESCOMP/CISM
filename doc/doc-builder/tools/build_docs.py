#!/usr/bin/env python3

"""Main driver wrapper around the doc_builder/build_docs utility.

This tool wraps the build command to build sphinx-based documentation.

This should be run from the directory that contains the Makefile for
building the documentation.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from doc_builder import build_docs  # pylint: disable=import-error,wrong-import-position

if __name__ == "__main__":
    build_docs.main()
