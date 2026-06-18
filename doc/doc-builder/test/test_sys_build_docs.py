#!/usr/bin/env python3
"""
High-level system tests of the whole build_docs
"""

import unittest
import tempfile
import shutil
import os

# pylint: disable=import-error,no-name-in-module
from test.test_utils.git_helpers import (
    make_git_repo,
    add_git_commit,
    checkout_git_branch,
)
from doc_builder import build_docs


class TestBuildDocs(unittest.TestCase):
    """High-level system tests of build_docs"""

    # Allow long method names
    # pylint: disable=invalid-name

    # ------------------------------------------------------------------------
    # Helper methods
    # ------------------------------------------------------------------------

    def setUp(self):
        """Set up temporary source and repo-root directories, and chdir to
        the source directory"""
        self._return_dir = os.getcwd()
        self._sourcedir = tempfile.mkdtemp()
        self._build_reporoot = tempfile.mkdtemp()
        self._build_versions_dir = os.path.join(self._build_reporoot, "versions")
        os.makedirs(self._build_versions_dir)
        os.chdir(self._sourcedir)

    def tearDown(self):
        os.chdir(self._return_dir)
        shutil.rmtree(self._sourcedir, ignore_errors=True)
        shutil.rmtree(self._build_reporoot, ignore_errors=True)

    @staticmethod
    def write_makefile():
        """Write a fake makefile in the current directory

        The 'html' target of this Makefile just results in the text
        'hello world' being written to BUILDDIR/testfile.
        """

        makefile_contents = """
html:
\t@mkdir -p $(BUILDDIR)
\t@echo "hello world" > $(BUILDDIR)/testfile
"""

        with open("Makefile", "w", encoding="utf-8") as makefile:
            makefile.write(makefile_contents)

    def assert_file_contents_equal(self, expected, filepath, msg=None):
        """Asserts that the contents of the file given by 'filepath' are equal to
        the string given by 'expected'. 'msg' gives an optional message to be
        printed if the assertion fails."""

        with open(filepath, "r", encoding="utf-8") as myfile:
            contents = myfile.read()

        self.assertEqual(expected, contents, msg=msg)

    # ------------------------------------------------------------------------
    # Begin tests
    # ------------------------------------------------------------------------

    def test_with_repo_root(self):
        """Test with repo-root specified, but doc-version not specified (so doc version needs to
        be determined via git commands).
        """

        self.write_makefile()
        make_git_repo()
        add_git_commit()
        checkout_git_branch("foo_branch")
        build_path = os.path.join(self._build_versions_dir, "foo_branch")
        os.makedirs(build_path)

        args = ["--repo-root", self._build_reporoot]
        build_docs.main(args)

        self.assert_file_contents_equal(
            expected="hello world\n", filepath=os.path.join(build_path, "testfile")
        )

    def test_multiple_versions(self):
        """Test with multiple versions being specified at once"""

        self.write_makefile()
        build_path1 = os.path.join(self._build_versions_dir, "v1")
        build_path2 = os.path.join(self._build_versions_dir, "v2")
        # Note that, since we're specifying the versions explicitly, we
        # shouldn't need to make the v1 and v2 directories ahead of time.

        args = ["--repo-root", self._build_reporoot, "--doc-version", "v1", "v2"]
        build_docs.main(args)

        self.assert_file_contents_equal(
            expected="hello world\n", filepath=os.path.join(build_path1, "testfile")
        )
        self.assert_file_contents_equal(
            expected="hello world\n", filepath=os.path.join(build_path2, "testfile")
        )


if __name__ == "__main__":
    unittest.main()
