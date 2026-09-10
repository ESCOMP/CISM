"""
Implementation of the top-level logic for build_docs.
"""

import subprocess
import argparse
import os
import random
import string
import sys
from urllib.parse import urlparse
import signal

# pylint: disable=import-error,no-name-in-module
from doc_builder.build_commands import (
    get_build_dir,
    get_build_command,
    get_container_cli_tool,
    DEFAULT_IMAGE,
)
from doc_builder.build_docs_shared_args import bd_dir_group, bd_parser
from doc_builder.output_utils import extract_sphinx_complaints
from doc_builder.sys_utils import is_mac


def is_web_url(url_string):
    """
    Checks if a string is a valid web URL.

    Args:
        url_string: The string to check.

    Returns:
        True if the string is a valid web URL, False otherwise.
    """
    result = urlparse(url_string)
    return all([result.scheme, result.netloc])


def commandline_options(cmdline_args=None):
    """Process the command-line arguments.

    cmdline_args, if present, should be a string giving the command-line
    arguments. This is typically just used for testing.
    """

    description = """
This tool wraps the build command to build sphinx-based documentation.

This tool assists with creating the correct documentation build commands
in cases including:
- Building the documentation from a container
- Building versioned documentation, where the documentation builds land
  in subdirectories named based on the source branch

This tool should be put somewhere in your path. Then it should be run
from the directory that contains the Makefile for building the
documentation.

Simple usage is:

    build_docs -b /path/to/doc/build/repo/some/subdirectory [-c] [-d] [-i]

    Common additional flags are:
    -c: Before building, run 'make clean'
    -d: Use the {DEFAULT_IMAGE} version of the container to build the documentation

Usage for automatically determining the subdirectory in which to build,
based on the version indicated by the current branch, is:

    ./build_docs -r /path/to/doc/build/repo [-v DOC_VERSION]

    This will build the documentation in a subdirectory of the doc build
    repo, where the subdirectory is built from DOC_VERSION. If
    DOC_VERSION isn't given, it will be determined based on the git
    branch name in the doc source repository.

    In the above example, documentation will be built in:
    /path/to/doc/build/repo/versions/DOC_VERSION

    This usage also accepts the optional arguments described above.
"""

    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )

    dir_group = parser.add_mutually_exclusive_group(required=True)

    dir_group.add_argument(
        "-b",
        "--build-dir",
        default=None,
        help="Full path to the directory in which the doc build should go.",
    )

    # Add argument(s) to dir_group that are also in build_docs_to_publish's parser
    dir_group = bd_dir_group(dir_group)

    # Add argument(s) to parser that are also in build_docs_to_publish's parser
    parser = bd_parser(parser)

    parser.add_argument(
        "-v",
        "--doc-version",
        nargs="+",
        default=[None],
        help="Version name to build,\n"
        "corresponding to a directory name under repo root.\n"
        "Not applicable if --build-dir is specified.\n"
        "Multiple versions can be specified, in which case a build\n"
        "will be done for each version (with the same source).",
    )

    parser.add_argument(
        "--version-display-name",
        default=None,
        help="Version name for display in dropdown menu. If absent, uses -v/--version.",
    )

    parser.add_argument(
        "-c", "--clean", action="store_true", help="Before building, run 'make clean'."
    )

    parser.add_argument(
        "-i",
        "--container-image",
        default=None,
        help="Container image version to use. Implies -d.",
    )

    parser.add_argument(
        "-t",
        "--build-target",
        default="html",
        help="Target for the make command.\n" "Default is 'html'.",
    )

    parser.add_argument(
        "--num-make-jobs",
        default=4,
        help="Number of parallel jobs to use for the make process.\n" "Default is 4.",
    )

    parser.add_argument(
        "--versions",
        action="store_true",
        help="Build multiple versions of the docs, with drop-down switcher menu.",
    )

    parser.add_argument(
        "-w",
        "--warnings-as-warnings",
        action="store_true",
        help="Treat sphinx warnings as warnings, not errors.",
    )

    options = parser.parse_args(cmdline_args)

    if options.verbose:
        print(f"options: {options}")

    if options.versions:
        if not options.site_root:
            raise RuntimeError("--site-root must be provided when --versions is enabled")
        if not is_web_url(options.site_root) and not os.path.isabs(options.site_root):
            raise RuntimeError(
                f"--site-root is neither a web URL nor an absolute path: '{options.site_root}'"
            )

    if options.container_image:
        options.container_image = options.container_image.lower()
        options.build_in_container = True
    elif options.build_in_container:
        options.container_image = DEFAULT_IMAGE

    return options


def setup_env_var(build_command, env, env_var, value, container):
    """
    Set up an environment variable, depending on whether using container or not
    """
    if container:
        # Need to pass to container via the build command
        build_command.insert(-3, "-e")
        build_command.insert(-3, f"{env_var}={value}")
    else:
        env[env_var] = value
    return build_command, env


def start_container_software(cmd):
    """
    Starts our Podman machine or Docker Desktop, if needed.
    """
    subprocess.run(
        cmd,
        shell=True,
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )


_SPHINX_BUILD_FINISHED_WITH_PROBLEMS = "build finished with problems"

_MSG_BUILD_FAILED = "Documentation build failed."
_MSG_BUILD_COMPLETED_WITH_PROBLEMS = (
    "Documentation build completed, but with problems that must be resolved."
)


def _report_build_failure(err, verbose):
    """Report a build failure to the user, filtering output based on verbosity.

    In verbose mode, dumps full stdout and stderr. In non-verbose mode, shows
    only Sphinx WARNING/ERROR lines plus a hint to re-run with --verbose.
    """
    stdout_text = err.stdout.decode("utf-8", errors="replace")
    stderr_text = err.stderr.decode("utf-8", errors="replace")
    if verbose:
        sys.stdout.write(stdout_text)
        sys.stderr.write(stderr_text)
    else:
        complaints = extract_sphinx_complaints(stdout_text, stderr_text)
        if complaints:
            sys.stderr.write("\n".join(complaints) + "\n")
        if _SPHINX_BUILD_FINISHED_WITH_PROBLEMS in stderr_text + stdout_text:
            sys.stderr.write(_MSG_BUILD_COMPLETED_WITH_PROBLEMS + "\n")
        else:
            sys.stderr.write(_MSG_BUILD_FAILED + "\n")
        sys.stderr.write("Re-run with --verbose for full output.\n")


def _maybe_start_container(build_command):
    """Start container software/VM if the build command uses one."""
    if "podman" in build_command:
        start_container_software("podman machine start")
    elif "docker" in build_command and is_mac():
        start_container_software("docker desktop start")


def run_build_command(build_command, version, options):
    """Echo and then run the given build command"""
    verbose = options.verbose

    if not verbose:
        if build_command[-1] == "clean":
            print("Cleaning documentation build directory...")
        else:
            print("Building documentation...")

    env = os.environ.copy()

    # Set version display name (in drop-down menu)
    if options.version_display_name:
        value = options.version_display_name
    else:
        value = version
    build_command, env = setup_env_var(
        build_command, env, "version_display_name", value, options.build_in_container
    )

    # Set paths to certain directories
    build_command, env = setup_env_var(
        build_command, env, "html_static_path", options.static_path, options.build_in_container
    )
    build_command, env = setup_env_var(
        build_command, env, "templates_path", options.templates_path, options.build_in_container
    )

    # Things to do/set based on whether including version dropdown
    if options.versions:
        version_dropdown = "True"
        build_command, env = setup_env_var(
            build_command,
            env,
            "pages_root",
            options.site_root,
            options.build_in_container,
        )
    else:
        version_dropdown = ""
    build_command, env = setup_env_var(
        build_command,
        env,
        "version_dropdown",
        version_dropdown,
        options.build_in_container,
    )

    _maybe_start_container(build_command)

    if verbose:
        print(" ".join(build_command))

    result = _try_build_command(build_command, verbose, env)

    # Print "The HTML pages are in ..." bit, which if we're not verbose is part of the suppressed
    # stdout
    if not verbose and build_command[-1] == "html":
        result_lines = result.stdout.decode("utf-8").splitlines()
        for line in result_lines:
            if line.startswith("The HTML pages are in"):
                print(line)
                break

    print("Done.")


def _try_build_command(build_command, verbose, env, allow_recursion=True):
    """Try the docs build command, retrying via recursion if needed"""
    try:
        result = subprocess.run(build_command, env=env, check=True, capture_output=True)
    except subprocess.CalledProcessError as err:
        stderr_text = err.stderr.decode("utf-8", errors="replace")

        # Unless we got this specific error message, exit.
        if "failed to chown recursively host path" not in stderr_text or not allow_recursion:
            _report_build_failure(err, verbose)
            # We do sys.exit(1) instead of raise if we're not verbose so we don't see the whole
            # Python traceback, which can distract from the ultimate source of the problem.
            if verbose:
                raise
            sys.exit(1)

        # If we made it here, we're going to change user mount and retry via recursion
        if verbose:
            print("Container failed due to missing subuid/subgid mappings.")
            print("Retrying without :U mount flag and with --user 0:0...")
        build_command = _fix_command_for_missing_subids(build_command)
        if verbose:
            print(" ".join(build_command))
        result = _try_build_command(build_command, verbose, env, allow_recursion=False)
    return result


def _fix_command_for_missing_subids(build_command):
    """Modify a container build command to work without subuid/subgid mappings.

    This removes the :U suffix from volume mounts and changes --user to 0:0
    (in rootless podman, UID 0 in the container maps to the host user).
    """
    fixed = []
    for i, arg in enumerate(build_command):
        if arg.endswith(":U"):
            fixed.append(arg[:-2])
        elif i > 0 and build_command[i - 1] == "--user":
            fixed.append("0:0")
        else:
            fixed.append(arg)
    return fixed


def setup_for_container():
    """Do some setup for running with container

    Returns a name that should be used in the container run command
    """

    container_name = "build_docs_" + "".join(
        random.choice(string.ascii_lowercase) for _ in range(8)
    )

    # It seems that, if we kill the build_docs process with Ctrl-C, the container process
    # continues. Handle that by implementing a signal handler. There may be a better /
    # more pythonic way to handle this, but this should work.
    def sigint_kill_container(signum, frame):
        """Signal handler: kill container process before exiting"""
        # pylint: disable=unused-argument
        container_kill_cmd = [get_container_cli_tool(), "kill", container_name]
        subprocess.check_call(container_kill_cmd)
        sys.exit(1)

    signal.signal(signal.SIGINT, sigint_kill_container)

    return container_name


def main(cmdline_args=None):
    """Top-level function implementing build_docs.

    cmdline_args, if present, should be a string giving the command-line
    arguments. This is typically just used for testing.
    """
    opts = commandline_options(cmdline_args)

    if opts.build_in_container:
        # We potentially reuse the same name for multiple container processes: the
        # clean and the actual build. However, since a given process should end before the
        # next one begins, and because we use '--rm' in the container run command, this
        # should be okay.
        container_name = setup_for_container()
    else:
        container_name = None

    # Note that we do a separate build for each version. This is
    # inefficient (assuming that the desired end result is for the
    # different versions to be identical), but was an easy-to-implement
    # solution to add convenience for building multiple versions of
    # documentation with short build times (i.e., rather than requiring
    # you to rerun build_docs multiple times). If this
    # multiple-versions-at-once option starts to be used a lot, we could
    # reimplement it to build just one version then copy the builds to
    # the other versions (if that gives the correct end result).
    for version in opts.doc_version:
        build_dir, version = get_build_dir(
            build_dir=opts.build_dir, repo_root=opts.repo_root, version=version
        )

        if opts.clean:
            clean_command = get_build_command(
                build_dir=build_dir,
                run_from_dir=os.getcwd(),
                build_target="clean",
                num_make_jobs=opts.num_make_jobs,
                version=version,
                container_name=container_name,
                container_image=opts.container_image,
                container_cli_tool=opts.container_cli_tool,
            )
            run_build_command(build_command=clean_command, version=version, options=opts)

        build_command = get_build_command(
            build_dir=build_dir,
            run_from_dir=os.getcwd(),
            build_target=opts.build_target,
            num_make_jobs=opts.num_make_jobs,
            version=version,
            container_name=container_name,
            container_image=opts.container_image,
            warnings_as_warnings=opts.warnings_as_warnings,
            conf_py_path=opts.conf_py_path,
            container_cli_tool=opts.container_cli_tool,
        )
        run_build_command(build_command=build_command, version=version, options=opts)
