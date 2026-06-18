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
from doc_builder.sys_utils import is_under


class TestIsUnder(unittest.TestCase):
    """Test the is_under function"""

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

    def test_isunder_true(self):
        """Test if child is under parent"""
        my_dir = "some_dir"
        os.makedirs(my_dir)
        filename = os.path.join(my_dir, "some_file")
        with open(filename, "a", encoding="utf8") as myfile:
            myfile.write("ietbeconae0nr0einr")
        self.assertTrue(is_under(filename, my_dir))
        self.assertFalse(is_under(my_dir, filename))

    def test_isunder_true_extraslash(self):
        """Test if child is under parent and parent has a path separator at the end"""
        my_dir = "some_dir" + os.sep
        os.makedirs(my_dir)
        filename = os.path.join(my_dir, "some_file")
        with open(filename, "a", encoding="utf8") as myfile:
            myfile.write("ietbeconae0nr0einr")
        self.assertTrue(is_under(filename, my_dir))
        self.assertFalse(is_under(my_dir, filename))

    def test_isunder_false(self):
        """Test if true"""
        my_dir = "some_dir"
        os.makedirs(my_dir)
        filename = "some_file"
        with open(filename, "a", encoding="utf8") as myfile:
            myfile.write("ietbeconae0nr0einr")
        self.assertFalse(is_under(filename, my_dir))

    def test_isunder_true_neither_exist(self):
        """Test true even if neither exists"""
        my_dir = "some_dir"
        filename = "some_file"
        self.assertFalse(is_under(filename, my_dir))


if __name__ == "__main__":
    unittest.main()
