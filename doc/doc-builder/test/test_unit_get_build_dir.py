#!/usr/bin/env python3

"""Unit test driver for get_build_dir function"""

import shutil
import unittest

try:
    # For python2; needs pip install mock
    import mock
except ImportError:
    # For python3
    from unittest import mock
import os

# pylint: disable=import-error,no-name-in-module
from test.test_utils.sys_utils_fake import (
    make_fake_isdir,
)
from doc_builder.build_commands import get_build_dir


class TestGetBuildDir(unittest.TestCase):
    """Test the get_build_dir function"""

    # Allow long method names
    # pylint: disable=invalid-name

    def setUp(self):
        """Run this before each test"""

        self.path_to_foo = os.path.join("path", "to", "foo")
        self.path_to_repo = os.path.join("path", "to", "repo")
        self.path_to_parent = os.path.join("path", "to", "parent")
        self.path_to_parent_repo = os.path.join(self.path_to_parent, "repo")

    def tearDown(self):
        """Run this after each test, whether it succeeded or failed"""
        test_dir_list = [
            self.path_to_foo,
            self.path_to_repo,
            self.path_to_parent,
            self.path_to_parent_repo,
        ]
        for test_dir in test_dir_list:
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)

    def test_with_builddir(self):
        """If given a build_dir, should return that"""
        build_path = self.path_to_foo
        build_dir = get_build_dir(build_dir=build_path, repo_root=None, version=None)
        expected = (build_path, "None")
        self.assertEqual(expected, build_dir)

    def test_builddir_and_reporoot(self):
        """If given both build_dir and repo_root, should raise an exception"""
        with self.assertRaises(RuntimeError):
            _ = get_build_dir(build_dir=self.path_to_foo, repo_root=self.path_to_repo, version=None)

    def test_builddir_and_version(self):
        """If given both build_dir and version, should raise an exception"""
        with self.assertRaises(RuntimeError):
            _ = get_build_dir(build_dir=self.path_to_foo, repo_root=None, version="v1.0")

    def test_no_builddir_or_reporoot(self):
        """If given neither build_dir nor repo_root, should raise an exception"""
        with self.assertRaises(RuntimeError):
            _ = get_build_dir(build_dir=None, repo_root=None, version=None)

    def test_reporoot_and_version(self):
        """If given both repo_root and version, should return correct build_dir"""
        ver = "v1.0"
        with mock.patch("os.path.isdir") as mock_isdir:
            path_to_repo_versions = os.path.join(self.path_to_repo, "versions")
            # /path/to/repo/versions exists; with version specified explicitly,
            # /path/to/repo/versions/v1.0 doesn't need to exist
            mock_isdir.side_effect = make_fake_isdir(dirs_exist=[path_to_repo_versions])
            build_dir = get_build_dir(build_dir=None, repo_root=self.path_to_repo, version=ver)
        expected = (os.path.join(path_to_repo_versions, ver), ver)
        self.assertEqual(expected, build_dir)

    def test_reporoot_and_version_dir_not_exist(self):
        """
        If given repo_root and version, with repo_root/versions not existing, should be okay.
        get_build_dir() should make any needed intermediate dirs.
        """
        get_build_dir(build_dir=None, repo_root=self.path_to_parent_repo, version="v1.0")

    def test_reporoot_no_version(self):
        """If given repo_root but no version, get version from git branch"""
        ver = "release-v2.0"
        with mock.patch("doc_builder.sys_utils.git_current_branch") as mock_git_current_branch:
            with mock.patch("os.path.isdir") as mock_isdir:
                mock_git_current_branch.return_value = (True, ver)
                path_to_repo = os.path.join("path", "to", "repo")
                path_to_versions = os.path.join(path_to_repo, "versions")
                expected_dir = os.path.join(path_to_versions, ver)
                expected = (expected_dir, ver)
                mock_isdir.side_effect = make_fake_isdir(
                    dirs_exist=[path_to_versions[0], expected_dir]
                )
                build_dir = get_build_dir(build_dir=None, repo_root=path_to_repo, version=None)

        self.assertEqual(expected, build_dir)

    def test_reporoot_no_version_git_branch_problem(self):
        """If given repo_root but no version, with a problem getting git
        branch, should raise an exception."""
        with mock.patch("doc_builder.sys_utils.git_current_branch") as mock_git_current_branch:
            mock_git_current_branch.return_value = (False, "")
            with self.assertRaises(RuntimeError):
                _ = get_build_dir(build_dir=None, repo_root=self.path_to_repo, version=None)

    def test_reporoot_no_version_dir_not_exist(self):
        """If given repo_root but no version, with the expected
        directory not existing, should raise an exception."""
        with mock.patch("doc_builder.sys_utils.git_current_branch") as mock_git_current_branch:
            with mock.patch("os.path.isdir") as mock_isdir:
                mock_git_current_branch.return_value = (True, "release-v2.0")
                path_to_repo = os.path.join("path", "to", "repo")
                path_to_repo_versions = os.path.join(path_to_repo, "versions")
                # /path/to/repo exists, but /path/to/repo/release-v2.0 does not
                mock_isdir.side_effect = make_fake_isdir(dirs_exist=[path_to_repo_versions])
                with self.assertRaises(RuntimeError):
                    _ = get_build_dir(build_dir=None, repo_root=path_to_repo, version=None)


if __name__ == "__main__":
    unittest.main()
