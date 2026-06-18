#!/usr/bin/env python3

"""Unit tests for run_build_command output behavior"""

import subprocess
import unittest
from types import SimpleNamespace
from unittest.mock import patch, MagicMock
from io import StringIO

from doc_builder.build_docs import (  # pylint: disable=import-error
    run_build_command,
    _maybe_start_container,
    _MSG_BUILD_FAILED,
    _MSG_BUILD_COMPLETED_WITH_PROBLEMS,
    _SPHINX_BUILD_FINISHED_WITH_PROBLEMS,
)


# A minimal options object for testing (non-container)
_BASE_OPTIONS = {
    "build_in_container": False,
    "version_display_name": None,
    "versions": False,
    "static_path": "_static",
    "templates_path": "_templates",
}

_FAKE_COMMAND = ["make", "SPHINXOPTS=", "BUILDDIR=/tmp/build", "-j", "4", "html"]
_FAKE_VERSION = "main"

# Simulated Sphinx output with warnings buried in noise
_SPHINX_STDOUT_WITH_WARNINGS = (
    "Running Sphinx v4.5.0\n"
    "building [html]: targets for 3 source files\n"
    "/path/to/file.rst:42: WARNING: toctree contains reference to nonexisting document 'missing'\n"
    "writing output... [100%] index\n"
    "WARNING: unknown config value 'bogus'\n"
)
_SPHINX_STDERR_WITH_ERROR = "/path/to/file.rst:42: ERROR: master file not found\n"
_SPHINX_STDERR_FINISHED_WITH_PROBLEMS = (
    "WARNING: unknown config value 'bogus'\n" + _SPHINX_BUILD_FINISHED_WITH_PROBLEMS + "\n"
)


def _make_options(verbose):
    """Create a test options namespace with the given verbose setting"""
    return SimpleNamespace(**_BASE_OPTIONS, verbose=verbose)


def _make_called_process_error(stdout_text, stderr_text):
    """Create a CalledProcessError with the given stdout/stderr"""
    err = subprocess.CalledProcessError(returncode=2, cmd=_FAKE_COMMAND)
    err.stdout = stdout_text.encode("utf-8")
    err.stderr = stderr_text.encode("utf-8")
    return err


class TestRunBuildCommandOutput(unittest.TestCase):
    """Tests for run_build_command output in verbose vs. non-verbose mode"""

    @patch("subprocess.run")
    def test_non_verbose_prints_building_message(self, mock_run):
        """In non-verbose mode, prints 'Building documentation...' before building"""
        mock_run.return_value = MagicMock(returncode=0)
        opts = _make_options(verbose=False)
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            run_build_command(_FAKE_COMMAND, _FAKE_VERSION, opts)
        self.assertIn("Building documentation...", mock_stdout.getvalue())

    @patch("subprocess.run")
    def test_success_non_verbose_prints_complete_message(self, mock_run):
        """On success in non-verbose mode, prints completion message"""
        mock_run.return_value = MagicMock(returncode=0)
        opts = _make_options(verbose=False)
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            run_build_command(_FAKE_COMMAND, _FAKE_VERSION, opts)
        self.assertIn("Done.", mock_stdout.getvalue())

    @patch("subprocess.run")
    def test_success_non_verbose_no_command_echo(self, mock_run):
        """On success in non-verbose mode, does not echo the build command"""
        mock_run.return_value = MagicMock(returncode=0)
        opts = _make_options(verbose=False)
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            run_build_command(_FAKE_COMMAND, _FAKE_VERSION, opts)
        self.assertNotIn("make", mock_stdout.getvalue())

    @patch("subprocess.run")
    def test_success_verbose_echoes_command(self, mock_run):
        """On success in verbose mode, echoes the build command"""
        mock_run.return_value = MagicMock(returncode=0)
        opts = _make_options(verbose=True)
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            run_build_command(_FAKE_COMMAND, _FAKE_VERSION, opts)
        self.assertIn("make", mock_stdout.getvalue())

    @patch("subprocess.run")
    def test_failure_non_verbose_shows_only_complaints(self, mock_run):
        """On failure in non-verbose mode, shows only WARNING/ERROR lines"""
        mock_run.side_effect = _make_called_process_error(
            _SPHINX_STDOUT_WITH_WARNINGS, _SPHINX_STDERR_WITH_ERROR
        )
        opts = _make_options(verbose=False)
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout, patch(
            "sys.stderr", new_callable=StringIO
        ) as mock_stderr:
            with self.assertRaises(SystemExit):
                run_build_command(_FAKE_COMMAND, _FAKE_VERSION, opts)
        combined = mock_stdout.getvalue() + mock_stderr.getvalue()
        self.assertIn("WARNING: toctree contains reference", combined)
        self.assertNotIn("WARNING: unknown config value", combined)  # Because WARNING at line start
        self.assertIn("ERROR: master file not found", combined)
        # Should NOT contain normal Sphinx noise
        self.assertNotIn("Running Sphinx", combined)
        self.assertNotIn("building [html]", combined)

    @patch("subprocess.run")
    def test_failure_non_verbose_shows_hint(self, mock_run):
        """On failure in non-verbose mode, shows hint to use --verbose"""
        mock_run.side_effect = _make_called_process_error(
            _SPHINX_STDOUT_WITH_WARNINGS, _SPHINX_STDERR_WITH_ERROR
        )
        opts = _make_options(verbose=False)
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout, patch(
            "sys.stderr", new_callable=StringIO
        ) as mock_stderr:
            with self.assertRaises(SystemExit):
                run_build_command(_FAKE_COMMAND, _FAKE_VERSION, opts)
        combined = mock_stdout.getvalue() + mock_stderr.getvalue()
        self.assertIn("--verbose", combined)

    @patch("subprocess.run")
    def test_failure_non_verbose_shows_failed_message(self, mock_run):
        """On hard failure (not 'finished with problems'), shows failed message"""
        mock_run.side_effect = _make_called_process_error(
            _SPHINX_STDOUT_WITH_WARNINGS, _SPHINX_STDERR_WITH_ERROR
        )
        opts = _make_options(verbose=False)
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout, patch(
            "sys.stderr", new_callable=StringIO
        ) as mock_stderr:
            with self.assertRaises(SystemExit):
                run_build_command(_FAKE_COMMAND, _FAKE_VERSION, opts)
        combined = mock_stdout.getvalue() + mock_stderr.getvalue()
        self.assertIn(_MSG_BUILD_FAILED, combined)
        self.assertNotIn(_MSG_BUILD_COMPLETED_WITH_PROBLEMS, combined)

    @patch("subprocess.run")
    def test_failure_non_verbose_finished_with_problems(self, mock_run):
        """When Sphinx says 'build finished with problems', shows softer message"""
        mock_run.side_effect = _make_called_process_error(
            _SPHINX_STDOUT_WITH_WARNINGS, _SPHINX_STDERR_FINISHED_WITH_PROBLEMS
        )
        opts = _make_options(verbose=False)
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout, patch(
            "sys.stderr", new_callable=StringIO
        ) as mock_stderr:
            with self.assertRaises(SystemExit):
                run_build_command(_FAKE_COMMAND, _FAKE_VERSION, opts)
        combined = mock_stdout.getvalue() + mock_stderr.getvalue()
        self.assertIn(_MSG_BUILD_COMPLETED_WITH_PROBLEMS, combined)
        self.assertNotIn(_MSG_BUILD_FAILED, combined)

    @patch("subprocess.run")
    def test_failure_verbose_shows_full_output(self, mock_run):
        """On failure in verbose mode, dumps full stdout and stderr"""
        mock_run.side_effect = _make_called_process_error(
            _SPHINX_STDOUT_WITH_WARNINGS, _SPHINX_STDERR_WITH_ERROR
        )
        opts = _make_options(verbose=True)
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout, patch(
            "sys.stderr", new_callable=StringIO
        ) as mock_stderr:
            with self.assertRaises(subprocess.CalledProcessError):
                run_build_command(_FAKE_COMMAND, _FAKE_VERSION, opts)
        combined = mock_stdout.getvalue() + mock_stderr.getvalue()
        # Verbose mode shows everything including normal lines
        self.assertIn("Running Sphinx", combined)
        self.assertIn("building [html]", combined)


class TestMaybeStartContainer(unittest.TestCase):
    """Tests for _maybe_start_container output"""

    @patch("doc_builder.build_docs.start_container_software")
    def test_no_output(self, _mock_start):
        """Container startup should not print anything"""
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            _maybe_start_container(["podman", "run", "image"])
        self.assertEqual("", mock_stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
