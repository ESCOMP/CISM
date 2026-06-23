#!/usr/bin/env python3

"""
Loop through all versions of the documentation, building each and moving it to a directory for
publication.

Adapted from https://www.codingwiththomas.com/blog/my-sphinx-best-practice-for-a-multiversion-
documentation-in-different-languages
(last visited 2025-05-20)
"""

import sys
import os
import subprocess
import argparse

# pylint: disable=import-error,no-name-in-module
from doc_builder.build_docs import (
    main as build_docs,
)
from doc_builder.build_docs_shared_args import main as build_docs_shared_args
from doc_builder.sys_utils import (
    get_git_head_or_branch,
    check_permanent_file,
    get_toplevel_git_dir,
)
from doc_builder.build_commands import get_container_cli_tool

# Change to the parent director of doc-builder and add to Python path
os.chdir(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, os.getcwd())

# Import our definitions of each documentation version.
# pylint: disable=wrong-import-position
from version_list import (
    LATEST_REF,
    VERSION_LIST,
)


# Path to certain important files
SOURCE = "source"
VERSIONS_PY = os.path.join("version_list.py")
MAKEFILE = "Makefile"


def get_static_templates_path_relative_to_here(args, path):
    """
    _static and _templates paths are relative to conf.py, but for operations here we need them
    relative to our current directory.
    """
    full_path = os.path.join(os.path.dirname(args.conf_py_path), path)
    return os.path.relpath(full_path)


def setup_this_ref(args):
    """
    Check that the repo is ready for the checkouts that are about to happen
    """
    # Get the current branch, or SHA if detached HEAD
    orig_ref = get_git_head_or_branch()

    # Some files/directories/submodules must stay the same for all builds. We list these in
    # the permanent_files list.
    permanent_files = [
        VERSIONS_PY,
        "doc-builder",
        MAKEFILE,
        args.conf_py_path,
        get_static_templates_path_relative_to_here(args, args.templates_path),
        get_static_templates_path_relative_to_here(args, args.static_path),
    ]

    # Check some things about "permanent" files before checkout
    for filename in permanent_files:
        check_permanent_file(filename)

    return orig_ref, permanent_files


def checkout_and_build(version, permanent_files, args):
    """
    Check out docs for a version and build
    """

    # Check out the git reference of this version (branch name, tag, or commit SHA)
    subprocess.check_output("git checkout " + version.ref, shell=True)

    # Check out LATEST_REF version of permanent files
    our_toplevel = get_toplevel_git_dir(os.getcwd())
    for filepath in permanent_files:
        # Skip if not in our working tree (i.e., in a submodule or a totally different place)
        file_toplevel = get_toplevel_git_dir(filepath)
        if file_toplevel != our_toplevel:
            continue

        # Check out file
        filename = os.path.basename(filepath)
        subprocess.check_output(f"git checkout {LATEST_REF} -- {filename}", shell=True)

    # Build the docs for this version
    build_args = [
        "-r",
        args.repo_root,
        "-v",
        version.short_name,
        "--version-display-name",
        version.display_name,
        "--versions",
        "--site-root",
        args.site_root,
        "--clean",
    ]
    if args.build_in_container:
        build_args += ["-d"]
    if args.conf_py_path:
        build_args += ["--conf-py-path", args.conf_py_path]
    if args.static_path:
        build_args += ["--static-path", args.static_path]
    if args.templates_path:
        build_args += ["--templates-path", args.templates_path]
    if args.container_cli_tool:
        build_args += ["--container-cli-tool", args.container_cli_tool]
    print(" ".join(build_args))
    build_docs(build_args)


def reset_git(orig_ref):
    """
    Restore the repo to how it was before we did our checkouts
    """
    # 1. Get the current ref's version of doc-builder to avoid "would be overwritten by checkout"
    #    errors.
    cmd = "git submodule update --checkout doc-builder".split(" ")
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    # If doc-builder not tracked by git, that's fine. Otherwise error.
    if result.returncode and "did not match any file(s) known to git" not in result.stderr:
        print("PWD: " + os.getcwd())
        print(result.stdout)
        print(result.stderr)
        raise subprocess.CalledProcessError(result.returncode, cmd)

    # 2. Check out the original git ref (branch or commit SHA)
    subprocess.check_output("git checkout " + orig_ref, shell=True)

    # 3. Restore the current version's doc-builder
    subprocess.check_output("git submodule update --checkout doc-builder", shell=True)


def check_version_list():
    """
    Check version list for problems
    """
    has_default = False
    for version in VERSION_LIST:
        # Expect at most one version with landing_version True
        if version.landing_version:
            if has_default:
                raise RuntimeError("Expected at most one version with landing_version True")
            has_default = True


def main():
    """
    Loop through all versions of the documentation, building each and moving it to a directory for
    publication.
    """
    # Set up parser
    parser = argparse.ArgumentParser()

    # Arguments shared with build_docs
    parser = build_docs_shared_args(parser)

    # Custom arguments for build_docs_to_publish
    parser.add_argument(
        "--publish-dir",
        default="_publish",
        help="Where the docs should be moved after being built",
    )

    # Parse arguments
    args = parser.parse_args()
    if args.container_cli_tool is None:
        args.container_cli_tool = get_container_cli_tool()

    # Check version list for problems
    check_version_list()

    # Loop over all documentation versions
    for version in VERSION_LIST:
        # Check that repo is ready and get our current state
        orig_ref, permanent_files = setup_this_ref(args)

        # Build this version
        try:
            checkout_and_build(version, permanent_files, args)
        # Restore the repo to how it was before we did our checkouts, even if checkout_and_build()
        # errored
        finally:
            reset_git(orig_ref)

        # Copy this version to the publication directory
        src = os.path.join(args.repo_root, "versions", version.short_name, "html")
        if version.landing_version:
            dst = args.publish_dir
        else:
            dst = os.path.join(args.publish_dir, version.short_name)
        if not os.path.exists(dst):
            os.makedirs(dst)
        subprocess.check_output(f"mv '{src}'/* '{dst}'/", shell=True)


if __name__ == "__main__":
    main()
