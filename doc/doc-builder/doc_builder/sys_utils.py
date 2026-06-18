"""
Functions that wrap system calls, including calls to the OS, git, etc.
"""

import subprocess
import os
import platform


def check_permanent_file(filepath):
    """
    Check a "permanent" file (one that we don't want to change between doc version builds)
    """

    # Ensure file exists
    if not os.path.exists(filepath):
        raise FileNotFoundError(filepath)

    # Change to directory this file is in. Necessary in case this file is in a submodule.
    orig_dir = os.getcwd()
    parent_dir, filename = os.path.split(filepath)
    if parent_dir:
        os.chdir(parent_dir)

    # Error if file contains uncommitted changes
    cmd = f"git add . && git diff --quiet {filename} && git diff --cached --quiet {filename}"
    try:
        subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as exception:
        subprocess.check_output("git reset", shell=True)  # Unstage files staged by `git add`
        msg = f"Important file/submodule may contain uncommitted changes: '{filename}'"
        raise RuntimeError(msg) from exception

    # Change back to original directory
    if parent_dir:
        os.chdir(orig_dir)


def get_git_head_or_branch():
    """
    Get the name of the current branch. If detached HEAD, get current commit SHA.
    """
    try:
        result = subprocess.run(
            ["git", "symbolic-ref", "--short", "-q", "HEAD"],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=True,
        )
        output = result.stdout.strip()
    except subprocess.CalledProcessError:
        output = ""

    if not output:
        # Fallback to commit SHA
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=True,
        )
        output = result.stdout.strip()

    return output


def get_toplevel_git_dir(file_or_dir):
    """
    Given a file or directory, get the top-level path of its git working tree.
    Return None if file isn't in a working tree.
    """
    # Change to directory, or file's parent directory
    orig_dir = os.getcwd()
    if os.path.isdir(file_or_dir):
        new_dir = file_or_dir
    else:
        new_dir = os.path.abspath(os.path.dirname(file_or_dir))
    os.chdir(new_dir)

    # Get the top-level path of working tree
    cmd = ["git", "rev-parse", "--show-toplevel"]
    result = subprocess.run(
        cmd,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.stderr:
        if "not a git repository" in result.stderr:
            toplevel = None
        else:
            print(result.stdout)
            print(result.stderr)
            raise subprocess.CalledProcessError(result.returncode, cmd)
    else:
        toplevel = result.stdout.strip()

    # Change back to original directory
    os.chdir(orig_dir)

    return toplevel


def get_toplevel_of_doc_builder_parent():
    """Get the top level of the repo of which doc-builder is a submodule"""
    topdir_doc_builder = get_toplevel_git_dir(__file__)
    topdir_parent = get_toplevel_git_dir(os.path.join(topdir_doc_builder, os.pardir))
    return topdir_parent


def git_current_branch():
    """Determines the name of the current git branch

    Returns a tuple, (branch_found, branch_name), where branch_found is
    a logical specifying whether a branch name was found for HEAD. (If
    branch_found is False, then branch_name is ''.) (branch_found will
    also be false if we're not in a git repository.)
    """
    cmd = ["git", "symbolic-ref", "--short", "-q", "HEAD"]
    with open(os.devnull, "w", encoding="utf-8") as devnull:
        try:
            # Suppress stderr because we don't want to clutter output with
            # git's message, e.g., if we're not in a git repository.
            branch_name = subprocess.check_output(cmd, stderr=devnull, universal_newlines=True)
        except subprocess.CalledProcessError:
            branch_found = False
            branch_name = ""
        else:
            branch_found = True
            branch_name = branch_name.strip()

    return branch_found, branch_name


def is_mac():
    """Returns True if running on a Mac"""
    return platform.system() == "Darwin"


def is_under(child, parent_dir):
    """Check whether child is in parent_dir

    Args:
        child (str): Path of child directory or file
        parent_dir (str): Path of parent directory
    """
    if parent_dir[-1] != os.sep:
        parent_dir += os.sep
    return child.startswith(parent_dir)
