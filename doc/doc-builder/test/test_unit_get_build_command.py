#!/usr/bin/env python3

"""Unit test driver for get_build_command function"""

import os
import unittest
from unittest.mock import patch

# pylint: disable=import-error
from doc_builder.build_commands import (
    get_build_command,
    get_mount_arg,
    get_container_cli_tool,
    DEFAULT_IMAGE,
)

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

# pylint: disable=line-too-long


class TestGetBuildCommand(unittest.TestCase):
    """Test the get_build_command function"""

    def setUp(self):
        """Run this before each test"""

        # Get current user's UID and GID
        uid = os.getuid()
        gid = os.getgid()
        self.uid_gid = f"{uid}:{gid}"

    def test_basic(self):
        """Tests basic usage"""
        build_command = get_build_command(
            build_dir="/path/to/foo",
            run_from_dir="/irrelevant/path",
            build_target="html",
            num_make_jobs=4,
            container_name=None,
            version="None",
        )
        expected = [
            "make",
            "SPHINXOPTS=-W --keep-going",
            "BUILDDIR=/path/to/foo",
            "-j",
            "4",
            "html",
        ]
        self.assertEqual(expected, build_command)

    def test_custom_conf_py_path(self):
        """Tests usage with --conf-py-path as file"""
        conf_py_path = os.path.join(os.path.dirname(__file__), "conf.py")
        build_command = get_build_command(
            build_dir="/path/to/foo",
            run_from_dir="/irrelevant/path",
            build_target="html",
            num_make_jobs=4,
            container_name=None,
            version="None",
            conf_py_path=conf_py_path,
        )
        expected = [
            "make",
            f"SPHINXOPTS=-W --keep-going -c '{os.path.dirname(conf_py_path)}'",
            "BUILDDIR=/path/to/foo",
            "-j",
            "4",
            "html",
        ]
        self.assertEqual(expected, build_command)

    def test_nonexistent_conf_py_path(self):
        """Tests error with --conf-py-path as nonexistent file"""
        conf_py_path = "nwirefeirourboub"
        with self.assertRaisesRegex(FileNotFoundError, "--conf-py-path not found"):
            get_build_command(
                build_dir="/path/to/foo",
                run_from_dir="/irrelevant/path",
                build_target="html",
                num_make_jobs=4,
                container_name=None,
                version="None",
                conf_py_path=conf_py_path,
            )

    def test_custom_conf_py_path_dir(self):
        """Tests usage with --conf-py-path as directory"""
        conf_py_path = os.path.dirname(__file__)
        build_command = get_build_command(
            build_dir="/path/to/foo",
            run_from_dir="/irrelevant/path",
            build_target="html",
            num_make_jobs=4,
            container_name=None,
            version="None",
            conf_py_path=conf_py_path,
        )
        expected = [
            "make",
            f"SPHINXOPTS=-W --keep-going -c '{conf_py_path}'",
            "BUILDDIR=/path/to/foo",
            "-j",
            "4",
            "html",
        ]
        self.assertEqual(expected, build_command)

    @patch("doc_builder.sys_utils.get_toplevel_of_doc_builder_parent")
    def test_container(self, mock_get_toplevel_of_doc_builder_parent):
        """Tests usage with container"""
        mock_get_toplevel_of_doc_builder_parent.return_value = "/path/to/username"
        conf_py_path = os.path.join(os.path.dirname(__file__), "conf.py")
        build_command = get_build_command(
            build_dir="/path/to/username/foorepos/foodocs/versions/main",
            run_from_dir="/path/to/username/foorepos/foocode/doc",
            build_target="html",
            num_make_jobs=4,
            container_name="foo",
            version="None",
            conf_py_path=conf_py_path,
            container_cli_tool="abc123",
        )
        mount_option, mount_arg = get_mount_arg("/path/to/username")
        expected = [
            "abc123",
            "run",
            "--name",
            "foo",
            "--user",
            self.uid_gid,
            mount_option,
            mount_arg,
            "--workdir",
            "/home/user/mounted_home/foorepos/foocode/doc",
            "-t",
            "--rm",
            "-e",
            "current_version=None",
            DEFAULT_IMAGE,
            "make",
            f"SPHINXOPTS=-W --keep-going -c '{os.path.dirname(conf_py_path)}'",
            "BUILDDIR=/home/user/mounted_home/foorepos/foodocs/versions/main",
            "-j",
            "4",
            "html",
        ]
        print("build_command: +", " ".join(build_command))
        self.assertEqual(expected, build_command)

    @patch("doc_builder.sys_utils.get_toplevel_of_doc_builder_parent")
    def test_container_relpath(self, mock_get_toplevel_of_doc_builder_parent):
        """Tests usage with container, with a relative path to build_dir"""
        mock_get_toplevel_of_doc_builder_parent.return_value = "/path/to/username"
        build_command = get_build_command(
            build_dir="../../foodocs/versions/main",
            run_from_dir="/path/to/username/foorepos/foocode/doc",
            build_target="html",
            num_make_jobs=4,
            container_name="foo",
            version="None",
            container_cli_tool="abc123",
        )
        mount_option, mount_arg = get_mount_arg("/path/to/username")
        expected = [
            "abc123",
            "run",
            "--name",
            "foo",
            "--user",
            self.uid_gid,
            mount_option,
            mount_arg,
            "--workdir",
            "/home/user/mounted_home/foorepos/foocode/doc",
            "-t",
            "--rm",
            "-e",
            "current_version=None",
            DEFAULT_IMAGE,
            "make",
            "SPHINXOPTS=-W --keep-going",
            "BUILDDIR=/home/user/mounted_home/foorepos/foodocs/versions/main",
            "-j",
            "4",
            "html",
        ]
        self.assertEqual(expected, build_command)

    @patch("doc_builder.sys_utils.get_toplevel_of_doc_builder_parent")
    def test_container_no_clitool_given(self, mock_get_toplevel_of_doc_builder_parent):
        """Tests usage with container"""
        mock_get_toplevel_of_doc_builder_parent.return_value = "/path/to/username"
        conf_py_path = os.path.join(os.path.dirname(__file__), "conf.py")
        build_command = get_build_command(
            build_dir="/path/to/username/foorepos/foodocs/versions/main",
            run_from_dir="/path/to/username/foorepos/foocode/doc",
            build_target="html",
            num_make_jobs=4,
            container_name="foo",
            version="None",
            conf_py_path=conf_py_path,
            container_cli_tool=None,
        )
        mount_option, mount_arg = get_mount_arg("/path/to/username")
        expected = [
            "run",
            "--name",
            "foo",
            "--user",
            self.uid_gid,
            mount_option,
            mount_arg,
            "--workdir",
            "/home/user/mounted_home/foorepos/foocode/doc",
            "-t",
            "--rm",
            "-e",
            "current_version=None",
            DEFAULT_IMAGE,
            "make",
            f"SPHINXOPTS=-W --keep-going -c '{os.path.dirname(conf_py_path)}'",
            "BUILDDIR=/home/user/mounted_home/foorepos/foodocs/versions/main",
            "-j",
            "4",
            "html",
        ]
        print("build_command: +", " ".join(build_command))
        self.assertEqual(build_command, [get_container_cli_tool()] + expected)

    @patch("doc_builder.sys_utils.get_toplevel_of_doc_builder_parent")
    def test_container_builddir_not_in_db_checkout(self, mock_get_toplevel_of_doc_builder_parent):
        """If build_dir is not in parent of doc-builder, should raise an exception"""
        mock_get_toplevel_of_doc_builder_parent.return_value = "/path/to/username"
        with self.assertRaisesRegex(RuntimeError, "build directory must be somewhere in"):
            _ = get_build_command(
                build_dir="/path/to/other/foorepos/foodocs/versions/main",
                run_from_dir="/path/to/username/foorepos/foocode/doc",
                build_target="html",
                num_make_jobs=4,
                container_name="foo",
                version="None",
            )

    @patch("doc_builder.sys_utils.get_toplevel_of_doc_builder_parent")
    def test_container_runfromdir_not_in_db_parent(self, mock_get_toplevel_of_doc_builder_parent):
        """If run_from_dir is not in parent of doc-builder, should raise an exception"""
        mock_get_toplevel_of_doc_builder_parent.return_value = "/path/to/username"
        with self.assertRaisesRegex(RuntimeError, "build_docs must be run from somewhere in"):
            _ = get_build_command(
                build_dir="/path/to/username/foorepos/foodocs/versions/main",
                run_from_dir="/path/to/other/foorepos/foocode/doc",
                build_target="html",
                num_make_jobs=4,
                container_name="foo",
                version="None",
            )

    @patch("doc_builder.sys_utils.get_toplevel_of_doc_builder_parent")
    def test_container_runfromdir_relative(self, mock_get_toplevel_of_doc_builder_parent):
        """If run_from_dir is relative, should raise an exception"""
        mock_get_toplevel_of_doc_builder_parent.return_value = "/path/to/username"
        with self.assertRaisesRegex(
            RuntimeError,
            "Expect absolute path; got",
        ):
            _ = get_build_command(
                build_dir="/path/to/username/foorepos/foodocs/versions/main",
                run_from_dir="../doc",
                build_target="html",
                num_make_jobs=4,
                container_name="foo",
                version="None",
            )


if __name__ == "__main__":
    unittest.main()
