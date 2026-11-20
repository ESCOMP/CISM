#!/usr/bin/env python3
"""Check that CISM-wrapper's namelist definition covers the Glissade config keys.

Every key that Glissade reads from a config section should have an entry in
CISM-wrapper's cime_config/namelist_definition_cism.xml, in the namelist group
of the same name as the section, so that CESM users can set it in user_nl_cism.
This script compares the two and reports:

  missing     keys that Glissade reads but the wrapper does not define
  misplaced   keys that the wrapper puts in a different group (section)
  unread      wrapper entries in a config group that Glissade does not read
              there (stale or misspelled names)
  values      options whose wrapper valid_values differ from Glissade's
              value names (the option_choice arrays in glissade_types.F90)
  category    entries whose <category> is not cism_config_<group>

The keys are taken from the GetValue, GetOption, GetLogical and GetDycore calls
in glissade/glissade_setup.F90: in the handle_<section> subroutines, and after a
GetSection(config, section, '<section>') call in other routines (the ocean and
atmosphere levels); the value
names from the option_choice arrays in glissade/glissade_types.F90.

Keys that are left out of the wrapper on purpose are listed in
ALLOWED_MISSING below, with the reason.

Usage:
  check_wrapper_keys.py <CISM-wrapper directory or namelist_definition_cism.xml>
                        [--cism <CISM directory>]   (default: this checkout)

Exit status: 0 if no problems are found, 1 otherwise.
"""

import argparse
import os
import re
import sys
import xml.etree.ElementTree as ET

# Keys that the wrapper leaves out on purpose: section -> {key: reason}
ALLOWED_MISSING = {
    'constants': {key: 'physical constants must stay consistent with the coupler'
                  for key in ('rhoi', 'rhoo', 'grav', 'shci', 'lhci', 'trpt')},
}

# Options whose wrapper valid_values are not compared with Glissade's names
SKIP_VALUES = {
    'dycore': 'GetDycore also accepts glide and the integers 0 and 2',
}

# Wrapper groups that are written to the CISM config file but read by GLAD,
# not by the Glissade config reader; they are not checked here
NOT_GLISSADE_GROUPS = {'glad_climate', 'esm_output'}

GET_CALL = re.compile(
    r"call\s+Get(Value|Option|Logical)\s*\(\s*section\s*,\s*'([^']+)'\s*,\s*(\w*)", re.I)


def strip_comment(line):
    """Remove a Fortran ! comment that is not inside a quoted string."""
    out, quote = [], None
    for c in line:
        if quote:
            if c == quote:
                quote = None
        elif c in "'\"":
            quote = c
        elif c == '!':
            break
        out.append(c)
    return ''.join(out)


def join_continuations(text):
    """Join Fortran continuation lines (ending in &) and drop comments."""
    lines, buf = [], ''
    for line in text.split('\n'):
        code = strip_comment(line).rstrip()
        if code.endswith('&'):
            buf += code[:-1] + ' '
        else:
            code = code.lstrip()
            if code.startswith('&'):
                code = code[1:]
            lines.append(buf + code)
            buf = ''
    return lines


def glissade_keys(cism_dir):
    """Return {section: {key: choices_name or None}} from glissade_setup.F90."""
    path = os.path.join(cism_dir, 'glissade', 'glissade_setup.F90')
    keys, section = {}, None
    for line in join_continuations(open(path).read()):
        m = re.match(r'\s*subroutine\s+handle_(\w+)\s*\(', line, re.I)
        if m:
            section = m.group(1).lower()
            keys.setdefault(section, {})
            continue
        if re.match(r'\s*end\s+subroutine', line, re.I):
            section = None
            continue
        # other routines (e.g. glissade_get_zocn) get a section by name
        m = re.search(r"call\s+GetSection\s*\(\s*config\s*,\s*section\s*,\s*'(\w+)'", line, re.I)
        if m:
            section = m.group(1).lower()
            keys.setdefault(section, {})
            continue
        if section is None:
            continue
        if re.search(r'call\s+GetDycore\s*\(', line, re.I):
            keys[section]['dycore'] = 'dycore_choices'
        for kind, key, third in GET_CALL.findall(line):
            choices = third.lower() if kind.lower() == 'option' else None
            keys[section][key] = choices
    return {s: k for s, k in keys.items() if k}


def glissade_choices(cism_dir):
    """Return {choices_name: [value names]} from the option_choice arrays."""
    path = os.path.join(cism_dir, 'glissade', 'glissade_types.F90')
    choices = {}
    for line in join_continuations(open(path).read()):
        m = re.match(r'\s*type\s*\(\s*option_choice\s*\).*::\s*(\w+)\s*=\s*\[(.*)\]', line, re.I)
        if m:
            names = re.findall(r"option_choice\s*\(\s*\w+\s*,\s*'([^']*)'", m.group(2), re.I)
            choices[m.group(1).lower()] = names
    return choices


def wrapper_entries(xml_path):
    """Return {key: (group, category, valid_values or None)} for the config entries."""
    root = ET.parse(xml_path).getroot()
    entries = {}
    for e in root.iter('entry'):
        category = (e.findtext('category') or '').strip()
        group = (e.findtext('group') or '').strip()
        if not category.startswith('cism_config'):
            continue
        valid = e.findtext('valid_values')
        valid = [v.strip() for v in valid.split(',')] if valid else None
        entries[e.get('id')] = (group, category, valid)
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('wrapper', help='CISM-wrapper directory or path to namelist_definition_cism.xml')
    parser.add_argument('--cism', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'),
                        help='CISM directory (default: the checkout that holds this script)')
    args = parser.parse_args()

    xml_path = args.wrapper
    if os.path.isdir(xml_path):
        xml_path = os.path.join(xml_path, 'cime_config', 'namelist_definition_cism.xml')
    if not os.path.isfile(xml_path):
        sys.exit(f'Error: {xml_path} not found')

    keys = glissade_keys(args.cism)
    choices = glissade_choices(args.cism)
    entries = wrapper_entries(xml_path)
    sections = set(keys)
    where = {key: sec for sec, ks in keys.items() for key in ks}

    problems = []

    # Glissade keys missing from the wrapper or in a different group
    for sec in keys:
        for key, choices_name in keys[sec].items():
            if key not in entries:
                if key in ALLOWED_MISSING.get(sec, {}):
                    continue
                problems.append(f'missing    [{sec}] {key}')
                continue
            group, category, valid = entries[key]
            if group != sec:
                problems.append(f'misplaced  {key}: wrapper group {group}, Glissade section [{sec}]')
            if choices_name and key not in SKIP_VALUES:
                names = choices.get(choices_name)
                if names is None:
                    problems.append(f'values     {key}: no option_choice array {choices_name} found')
                elif valid is None:
                    problems.append(f'values     {key}: wrapper has no valid_values; Glissade: {",".join(names)}')
                elif set(valid) != set(names):
                    extra = sorted(set(valid) - set(names))
                    lacking = sorted(set(names) - set(valid))
                    msg = f'values     {key}:'
                    if lacking:
                        msg += f' not in wrapper: {",".join(lacking)};'
                    if extra:
                        msg += f' not in Glissade: {",".join(extra)};'
                    problems.append(msg.rstrip(';'))

    # Wrapper entries in a Glissade config group that Glissade does not read there
    for key, (group, category, valid) in sorted(entries.items(), key=lambda kv: (kv[1][0], kv[0])):
        if group in NOT_GLISSADE_GROUPS or category == 'cism_config_history':
            continue
        if category != f'cism_config_{group}':
            problems.append(f'category   {key}: category {category}, group {group}')
        if key not in where:
            if group in sections or category.startswith('cism_config_'):
                problems.append(f'unread     {key} (wrapper group {group}): not read by Glissade')

    nkeys = sum(len(k) for k in keys.values())
    nallowed = sum(len(v) for v in ALLOWED_MISSING.values())
    print(f'Glissade reads {nkeys} keys in {len(keys)} sections; '
          f'{nallowed} are left out of the wrapper on purpose ({", ".join(ALLOWED_MISSING)}).')
    print(f'Wrapper: {len(entries)} config entries in {xml_path}')
    if problems:
        print(f'\n{len(problems)} problem(s):')
        for p in problems:
            print('  ' + p)
        return 1
    print('No problems found.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
