#!/usr/bin/env python3
"""Tests of check_permanent_file

These are integration tests, since they interact with the OS and git,
and so are slower than typical unit tests.
"""

import unittest
import tempfile
import shutil
import os

# pylint: disable=import-error,no-name-in-module
from test.test_utils.git_helpers import (
    make_git_repo,
    add_git_commit,
    add_file,
)
from doc_builder.sys_utils import check_permanent_file


class TestCheckPermanentFile(unittest.TestCase):
    """Test the check_permanent_file function"""

    # ------------------------------------------------------------------------
    # Helper methods
    # ------------------------------------------------------------------------

    def setUp(self):
        self._return_dir = os.getcwd()
        self._tempdir = os.path.realpath(tempfile.mkdtemp())
        os.chdir(self._tempdir)

    def tearDown(self):
        os.chdir(self._return_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    # ------------------------------------------------------------------------
    # Begin tests
    # ------------------------------------------------------------------------

    def test_file_doesnt_exist(self):
        """Test that correct error is thrown if file doesn't exist"""
        with self.assertRaises(FileNotFoundError):
            check_permanent_file("wniefuwegr972brt92u4br")

    def test_error_new_file(self):
        """Test that correct error is thrown if file exists but is untracked"""
        make_git_repo()
        filename = add_file()
        with self.assertRaisesRegex(
            RuntimeError, "Important file/submodule may contain uncommitted changes"
        ):
            check_permanent_file(filename)

    def test_error_file_uncommitted_changes(self):
        """
        Test that correct error is thrown if file exists and is tracked but has uncommitted changes
        """
        make_git_repo()
        filename = add_git_commit()
        with open(filename, "a", encoding="utf8") as myfile:
            myfile.write("Here are some new changes we won't commit")
        with self.assertRaisesRegex(
            RuntimeError, "Important file/submodule may contain uncommitted changes"
        ):
            check_permanent_file(filename)

    def test_clean_file_in_this_dir(self):
        """Test no error thrown if given a file in current dir with no uncommitted changes"""
        make_git_repo()
        filename = add_git_commit()
        check_permanent_file(filename)
        self.assertEqual(os.getcwd(), self._tempdir)  # Should still be in original directory

    def test_clean_file_in_subdir(self):
        """Test no error thrown if given a file in a different dir with no uncommitted changes"""
        make_git_repo()
        sub_dir = "some_subdir"
        os.makedirs(sub_dir)
        filename = os.path.join(sub_dir, "tmpfile")
        add_git_commit(filename)
        check_permanent_file(filename)
        self.assertEqual(os.getcwd(), self._tempdir)  # Should still be in original directory


if __name__ == "__main__":
    unittest.main()
