"""Utilities for processing and filtering build output."""

import re

# Requiring at least one character before the keyword prevents the "WARNING: image platform ...
# does not match the expected platform ..." message from matching.
_SPHINX_COMPLAINT_PATTERN = re.compile(r"^.+(?:WARNING|ERROR|CRITICAL).*$", re.MULTILINE)


def extract_sphinx_complaints(*outputs):
    """Extract Sphinx warning/error lines from one or more output strings.

    Args:
        *outputs: One or more strings of build output (e.g., stdout, stderr).

    Returns:
        List of lines containing WARNING or ERROR, in the order they appeared.
    """
    complaints = []
    for output in outputs:
        complaints.extend(_SPHINX_COMPLAINT_PATTERN.findall(output))
    return complaints
