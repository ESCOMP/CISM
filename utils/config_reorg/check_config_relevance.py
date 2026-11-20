#!/usr/bin/env python3
"""Compare the keys of CESM-generated CISM configs with the keys CISM logs as relevant.

CISM-wrapper's buildnml should write a key to cism.<icesheet>.config only if it is
relevant for the selected options, or if the user set it in user_nl_cism. Glissade logs
each config section in config-file style, listing only the keys that are relevant.
For each CESM run directory, this script compares the two and reports:

  not logged   keys written to the config that CISM does not log (should be empty,
               apart from the keys in NOT_LOGGED below and keys the user set)
  defaults     (with --defaults) keys CISM logs that the config does not set; CISM uses
               its own default for these

The log is the newest glc.log* in the run directory (gzipped or not). In a run with
several ice sheets, the config log of each ice sheet starts with its [grid] section, in
the order of icesheet_names. Keys set in the case's user_nl_cism* files (found through
the run directory's CASEROOT link or the --caseroot option) are reported separately.

Usage:
  check_config_relevance.py <run dir> [<run dir> ...] [--defaults] [--summary]

Exit status: 0 if no unexpected key is found, 1 otherwise.
"""

import argparse
import glob
import gzip
import os
import re
import sys

# Keys that the wrapper writes but CISM does not log: key -> reason
NOT_LOGGED = {
    'log_level': 'a local variable passed to glimmer_log, not logged',
}

# Sections that are not Glissade config sections (read by GLAD, or I/O sections)
SKIP_SECTIONS = re.compile(r'^(CF |GLAD climate$|esm_output$)')

SECTION_RE = re.compile(r'^\s*\[([^\]]+)\]\s*$')
CONFIG_KEY_RE = re.compile(r'^\s*([A-Za-z_][A-Za-z_0-9]*)\s*=')
LOG_KEY_RE = re.compile(r'^\s+([A-Za-z_][A-Za-z_0-9]*)\s+=\s')
NML_KEY_RE = re.compile(r'^\s*([A-Za-z_][A-Za-z_0-9]*)\s*=')


def config_keys(path):
    """Return {(section, key)} for the Glissade sections of a config file."""
    keys, section = set(), None
    for line in open(path):
        m = SECTION_RE.match(line)
        if m:
            section = m.group(1)
            continue
        m = CONFIG_KEY_RE.match(line)
        if m and section and not SKIP_SECTIONS.match(section):
            keys.add((section, m.group(1)))
    return keys


def read_log(run_dir, inst=''):
    logs = sorted(glob.glob(os.path.join(run_dir, f'glc{inst}.log*')), key=os.path.getmtime)
    if not logs:
        return None, None
    path = logs[-1]
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', errors='replace') as f:
        return path, f.read().split('\n')


def logged_keys(lines, nsheets):
    """Return a list (one per ice sheet) of {(section, key)} from the config-style log."""
    blocks, section = [], None
    for line in lines:
        m = SECTION_RE.match(line)
        if m and line.startswith(' '):          # log sections are indented; config lines are not
            section = m.group(1)
            if section == 'grid':
                blocks.append(set())
            continue
        if not blocks or section is None:
            continue
        if line.strip() == '':
            continue
        m = LOG_KEY_RE.match(line)
        if m:
            blocks[-1].add((section, m.group(1)))
        elif not line.startswith('     '):        # left the config log
            section = None
    return blocks[:nsheets]


def user_keys(caseroot, icesheet):
    keys = set()
    if not caseroot:
        return keys
    for name in ('user_nl_cism', f'user_nl_cism_{icesheet}'):
        path = os.path.join(caseroot, name)
        if os.path.isfile(path):
            for line in open(path):
                if line.lstrip().startswith('!'):
                    continue
                m = NML_KEY_RE.match(line)
                if m:
                    keys.add(m.group(1))
    return keys


def find_caseroot(run_dir, given):
    if given:
        return given
    link = os.path.join(run_dir, 'CASEROOT')
    if os.path.islink(link) or os.path.isdir(link):
        return os.path.realpath(link)
    if os.path.isfile(link):
        return open(link).read().strip()
    return None


def instances(run_dir):
    """Return the instance suffixes of a run: [''] or ['_0001', '_0002', ...] (multi-instance)."""
    insts = sorted({re.match(r'glc(_\d{4})\.log', os.path.basename(p)).group(1)
                    for p in glob.glob(os.path.join(run_dir, 'glc_[0-9][0-9][0-9][0-9].log*'))})
    return insts or ['']


def check_run(run_dir, args):
    nbad, unexpected = 0, {}
    for inst in instances(run_dir):
        n, u = check_instance(run_dir, inst, args)
        nbad += n
        for k, v in u.items():
            unexpected.setdefault(k, []).extend(v)
    return nbad, unexpected


def check_instance(run_dir, inst, args):
    configs = sorted(glob.glob(os.path.join(run_dir, f'cism{inst}.*.config')))
    log_path, lines = read_log(run_dir, inst)
    if not configs or lines is None:
        print(f'{run_dir}: no cism{inst}.*.config or glc{inst}.log; skipped')
        return 0, {}
    names = None
    for line in lines:
        m = re.match(r'\s*icesheet_names:\s*(.*)', line)
        if m:
            names = m.group(1).split()
            break
    if names is None:
        names = [os.path.basename(c).split('.')[-2] for c in configs]
    blocks = logged_keys(lines, len(names))
    caseroot = find_caseroot(run_dir, args.caseroot)
    nbad, unexpected = 0, {}
    print(f'{run_dir}\n  log: {os.path.basename(log_path)}')
    for sheet, logged in zip(names, blocks):
        cfg = os.path.join(run_dir, f'cism{inst}.{sheet}.config')
        if not os.path.isfile(cfg):
            print(f'  {sheet}: no {os.path.basename(cfg)}')
            continue
        written = config_keys(cfg)
        user = user_keys(caseroot, sheet)
        logged_names = {k for _, k in logged}
        extra = sorted((s, k) for s, k in written - logged if k not in NOT_LOGGED and k not in logged_names)
        moved = sorted((s, k) for s, k in written - logged if k in logged_names)
        by_user = [x for x in extra if x[1] in user]
        bad = [x for x in extra if x[1] not in user]
        print(f'  {sheet}: {len(written)} keys written, {len(logged)} logged; '
              f'{len(bad)} written but not logged' + (f', {len(by_user)} more set by the user' if by_user else ''))
        if not args.summary:
            for s, k in bad:
                print(f'    not logged  [{s}] {k}')
            for s, k in by_user:
                print(f'    user set    [{s}] {k} (not logged)')
            for s, k in moved:
                print(f'    other section: {k} written in [{s}], logged elsewhere')
            if args.defaults:
                for s, k in sorted(logged - written):
                    print(f'    default     [{s}] {k}')
        for x in bad:
            unexpected.setdefault(x, []).append(f'{os.path.basename(os.path.dirname(run_dir.rstrip("/")))}{inst}:{sheet}')
        nbad += len(bad)
    return nbad, unexpected


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('run_dirs', nargs='+')
    parser.add_argument('--caseroot', help='case directory with the user_nl_cism files (one run dir only)')
    parser.add_argument('--defaults', action='store_true', help='also list logged keys the config does not set')
    parser.add_argument('--summary', action='store_true', help='print counts per run, then the union of unexpected keys')
    args = parser.parse_args()
    total, union = 0, {}
    for d in args.run_dirs:
        n, u = check_run(d, args)
        total += n
        for k, v in u.items():
            union.setdefault(k, []).extend(v)
    if len(args.run_dirs) > 1 and union:
        print(f'\nWritten but not logged, over all runs ({len(union)} keys):')
        for (s, k), where in sorted(union.items()):
            print(f'  [{s}] {k}: {len(where)} run(s)')
    return 1 if total else 0


if __name__ == '__main__':
    sys.exit(main())
