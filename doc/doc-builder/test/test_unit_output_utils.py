#!/usr/bin/env python3

"""Unit tests for output_utils"""

import unittest
from doc_builder.output_utils import extract_sphinx_complaints  # pylint: disable=import-error


# Sample lines that represent Sphinx complaints
_WARNING_WITH_LOCATION = (
    "/path/to/file.rst:42: WARNING: toctree contains reference to nonexisting document 'missing'"
)
_CRITICAL_WITH_LOCATION = "/path/to/file.rst:42: CRITICAL: Duplicate ID: 'equation-n-cost-fix'"
_WARNING_WITHOUT_LOCATION = "WARNING: unknown config value 'foo'"
_CRITICAL_WITHOUT_LOCATION = "CRITICAL: Duplicate ID: 'equation-n-cost-fix'"
_ERROR_WITH_LOCATION = "/path/to/file.rst:10: ERROR: unknown directive type 'bogus'"
_ERROR_WITHOUT_LOCATION = "ERROR: master file not found"

# Sample lines that are NOT complaints
_NORMAL_LINES = [
    "Running Sphinx v4.5.0",
    "building [html]: targets for 1 source files",
    "writing output... [100%] index",
    "build succeeded.",
]


class TestExtractSphinxComplaints(unittest.TestCase):
    """Tests for extract_sphinx_complaints"""

    def test_empty_input(self):
        """Empty string returns no complaints"""
        self.assertEqual(extract_sphinx_complaints(""), [])

    def test_no_complaints(self):
        """Normal Sphinx output returns no complaints"""
        output = "\n".join(_NORMAL_LINES)
        self.assertEqual(extract_sphinx_complaints(output), [])

    def test_warning_with_location(self):
        """Captures WARNING lines with file:line prefix"""
        output = "\n".join(_NORMAL_LINES + [_WARNING_WITH_LOCATION])
        self.assertEqual(extract_sphinx_complaints(output), [_WARNING_WITH_LOCATION])

    def test_warning_without_location(self):
        """Does *not* capture WARNING lines without file:line prefix"""
        output = "\n".join(_NORMAL_LINES + [_WARNING_WITHOUT_LOCATION])
        self.assertEqual(extract_sphinx_complaints(output), [])

    def test_error_with_location(self):
        """Captures ERROR lines with file:line prefix"""
        output = "\n".join(_NORMAL_LINES + [_ERROR_WITH_LOCATION])
        self.assertEqual(extract_sphinx_complaints(output), [_ERROR_WITH_LOCATION])

    def test_error_without_location(self):
        """Does *not* capture ERROR lines without file:line prefix"""
        output = "\n".join(_NORMAL_LINES + [_ERROR_WITHOUT_LOCATION])
        self.assertEqual(extract_sphinx_complaints(output), [])

    def test_critical_with_location(self):
        """Captures CRITICAL lines with file:line prefix"""
        output = "\n".join(_NORMAL_LINES + [_CRITICAL_WITH_LOCATION])
        self.assertEqual(extract_sphinx_complaints(output), [_CRITICAL_WITH_LOCATION])

    def test_critical_without_location(self):
        """Does *not* capture CRITICAL lines without file:line prefix"""
        output = "\n".join(_NORMAL_LINES + [_CRITICAL_WITHOUT_LOCATION])
        self.assertEqual(extract_sphinx_complaints(output), [])

    def test_multiple_complaints(self):
        """Extracts multiple complaint lines in order"""
        output = "\n".join(
            [
                _NORMAL_LINES[0],
                _WARNING_WITH_LOCATION,
                _NORMAL_LINES[1],
                _ERROR_WITH_LOCATION,
                _WARNING_WITHOUT_LOCATION,
            ]
        )
        expected = [
            _WARNING_WITH_LOCATION,
            _ERROR_WITH_LOCATION,
        ]
        self.assertEqual(extract_sphinx_complaints(output), expected)

    def test_combines_stdout_and_stderr(self):
        """When given two strings, extracts from both in order"""
        stdout = "\n".join([_NORMAL_LINES[0], _WARNING_WITH_LOCATION])
        stderr = "\n".join([_ERROR_WITH_LOCATION, _NORMAL_LINES[1]])
        result = extract_sphinx_complaints(stdout, stderr)
        self.assertEqual(result, [_WARNING_WITH_LOCATION, _ERROR_WITH_LOCATION])

    def test_zero_arguments(self):
        """Calling with zero arguments returns empty list"""
        self.assertEqual(extract_sphinx_complaints(), [])


if __name__ == "__main__":
    unittest.main()
