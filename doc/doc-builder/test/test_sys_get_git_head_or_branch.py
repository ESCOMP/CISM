#!/usr/bin/env python3
"""Tests of get_git_head_or_branch

These are integration tests, since they interact with the OS and git,
and so are slower than typical unit tests.
"""

import unittest
import tempfile
import shutil
import os
import re
import subprocess

# pylint: disable=import-error,no-name-in-module
from test.test_utils.git_helpers import (
    make_git_repo,
    add_git_commit,
    checkout_git_branch,
    make_git_tag,
    checkout_git_ref,
)
from doc_builder.sys_utils import get_git_head_or_branch

# Pattern for testing whether a string matches a git SHA (40 hexadecimal characters)
SHA_PATTERN = re.compile("^[0-9a-f]{40}")


class TestGetGitHeadOrBranch(unittest.TestCase):
    """Test the get_git_head_or_branch function"""

    # ------------------------------------------------------------------------
    # Helper methods
    # ------------------------------------------------------------------------

    def setUp(self):
        self._return_dir = os.getcwd()
        self._tempdir = tempfile.mkdtemp()
        os.chdir(self._tempdir)

    def tearDown(self):
        os.chdir(self._return_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    # ------------------------------------------------------------------------
    # Begin tests
    # ------------------------------------------------------------------------

    def test_error_not_git_repo(self):
        """If not a git repository, should error"""
        with self.assertRaises(subprocess.CalledProcessError):
            get_git_head_or_branch()

    def test_on_branch(self):
        """If on a git branch, should return branch name"""
        make_git_repo()
        add_git_commit()
        checkout_git_branch("foo")
        self.assertEqual(get_git_head_or_branch(), "foo")

    def test_not_on_branch(self):
        """If not on a git branch, should return a commit SHA"""
        make_git_repo()
        add_git_commit()
        checkout_git_branch("foo")
        add_git_commit("file2")
        checkout_git_ref("HEAD~1")
        match = SHA_PATTERN.match(get_git_head_or_branch())
        self.assertTrue(match)

    def test_at_tag(self):
        """If not on a git branch but on a tag, should still return SHA"""
        make_git_repo()
        add_git_commit()
        checkout_git_branch("foo")
        add_git_commit("file2")
        checkout_git_ref("HEAD~1")
        make_git_tag("new-tag")
        self.assertNotEqual(get_git_head_or_branch(), "new-tag")
        match = SHA_PATTERN.match(get_git_head_or_branch())
        self.assertTrue(match)


if __name__ == "__main__":
    unittest.main()
