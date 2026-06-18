"""General-purpose functions to help with testing"""

import subprocess

def check_call_suppress_output(args):
    """Make a subprocess call with the given args, suppressing all output unless there's an error"""
    try:
        subprocess.run(
            args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
            shell=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}: {e.cmd}")
        print("stdout:")
        print(e.stdout)
        print("stderr:")
        print(e.stderr)
        raise
