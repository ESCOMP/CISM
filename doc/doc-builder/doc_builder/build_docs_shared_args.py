"""
build_docs and build_docs_to_publish share some args. This module adds them to a parser or parser
group.
"""

import os

# pylint: disable=import-error,no-name-in-module
from .build_commands import DEFAULT_IMAGE, COMPATIBLE_CLI_TOOLS


def bd_parser(parser, site_root_required=False):
    """
    Add arguments that build_docs has in its overall parser.

    # site_root_required: Should be True from build_docs_to_publish, False from build_docs
    """
    parser.add_argument(
        "--site-root",
        required=site_root_required,
        help="URL or absolute file path that should contain the top-level index.html",
    )
    parser.add_argument(
        "-d",
        "--build-in-container",
        action="store_true",
        help="Use a container to build the documentation,\n"
        "rather than relying on locally-installed versions of Sphinx, etc.\n"
        "This checks that a compatible container tool is installed and running on your system.\n"
        "\n"
        "NOTE: This mounts your repo checkout in the container.\n"
        "Therefore, both the current directory (containing the Makefile for\n"
        "building the documentation) and the documentation build directory\n"
        "must reside somewhere within your repo checkout."
        "\n"
        f"Default image: {DEFAULT_IMAGE}\n"
        "This can be changed with -i/--container-image.",
    )
    parser.add_argument(
        "--container-cli-tool",
        help="Container command-line tool to use. Options: " + ", ".join(COMPATIBLE_CLI_TOOLS),
        default=None,
    )
    parser.add_argument(
        "--conf-py-path",
        help="Path to conf.py",
        # For some reason Sphinx can't handle an absolute path here
        default=os.path.relpath(os.path.join(os.path.dirname(__file__), os.pardir, "conf.py")),
    )
    parser.add_argument(
        "--static-path",
        help="Path to _static/. If relative, must be relative to conf.py.",
        default="_static",
    )
    parser.add_argument(
        "--templates-path",
        help="Path to _templates/. If relative, must be relative to conf.py.",
        default="_templates",
    )
    parser.add_argument(
        "-V",
        "--verbose",
        action="store_true",
        help="Show full build output. Default shows only errors/warnings.",
    )
    return parser


def bd_dir_group(parser_or_group, repo_root_default=None):
    """
    Add arguments that build_docs has in its dir_group
    """
    parser_or_group.add_argument(
        "-r",
        "--repo-root",
        default=repo_root_default,
        help="Root directory of the repository holding documentation builds.\n"
        "(If there are other path elements between the true repo root and\n"
        "the 'versions' directory, those should be included in this path.)",
    )
    return parser_or_group


def main(parser):
    """
    Add all arguments to parser, even if build_docs has them in dir_group
    """

    # Settings for build_docs_to_publish, because main() should only ever be called from there
    site_root_required = True
    repo_root_default = "_build"

    parser = bd_parser(parser, site_root_required)
    parser = bd_dir_group(parser, repo_root_default)

    return parser
