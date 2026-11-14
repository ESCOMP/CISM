#!/usr/bin/env python3
"""Convert CISM (Glissade) config files, or CESM user_nl_cism files, to the reorganized config format.

The config reorganization moves keys out of the long [options], [ho_options] and [parameters]
sections into short sections by topic (e.g. [calving], [climate], [basal_physics]), and renames
some keys (e.g. marine_margin -> calving_law).  Option values are not changed.

Usage:
    convert_cism_config.py FILE [FILE ...]              write FILE.converted next to each FILE
    convert_cism_config.py FILE -o NEWFILE              write one converted file to NEWFILE
    convert_cism_config.py FILE [FILE ...] --in-place   overwrite each FILE
    convert_cism_config.py FILE [FILE ...] --check      report what would change; exit 1 if anything would
    convert_cism_config.py --self-test                  run the built-in tests

A file with no [section] lines is treated as a CESM user_nl_cism file, in which only key names change
(the CESM namelist has no sections).  Use --user-nl to force this.

How a config file is converted:
  - Each key listed in MOVES is moved to its new section, under its new name, with its value and
    inline comment.  Comment lines directly above the key (with no blank line in between) move with it.
  - A new section is created if needed, before the first [CF ...] section (or at the end of the file).
  - An old section left with no keys is removed, together with any comment lines left in it.
  - Keys not listed in MOVES are left where they are; the model warns about keys it does not read.
  - Converting a file that is already converted changes nothing.
  - Glide configs ([options] dycore = 0) are left unchanged: the Glide config format is frozen.

Only the standard library is used.
"""

import argparse
import re
import sys

# (old section, old key) -> (new section, new key).
# Filled in one topic at a time, in the commits of the config reorganization.
MOVES = {
    # [calving]
    ('options', 'marine_margin'): ('calving', 'calving_law'),
    ('options', 'calving_init'): ('calving', 'calving_init'),
    ('options', 'calving_domain'): ('calving', 'calving_domain'),
    ('options', 'apply_calving_mask'): ('calving', 'apply_calving_mask'),
    ('options', 'remove_icebergs'): ('calving', 'remove_icebergs'),
    ('options', 'remove_isthmuses'): ('calving', 'remove_isthmuses'),
    ('options', 'limit_marine_cliffs'): ('calving', 'limit_marine_cliffs'),
    ('ho_options', 'which_ho_calving_front'): ('calving', 'subgrid_calving_front'),
    ('ho_options', 'which_ho_calvingmip_domain'): ('calving', 'calvingmip_domain'),
    ('parameters', 'marine_limit'): ('calving', 'marine_limit'),
    ('parameters', 'calving_fraction'): ('calving', 'calving_fraction'),
    ('parameters', 'calving_timescale'): ('calving', 'calving_timescale'),
    ('parameters', 'calving_minthck'): ('calving', 'calving_minthck'),
    ('parameters', 'dthck_dx_cf'): ('calving', 'dthck_dx_cf'),
    ('parameters', 'thck_effective_min'): ('calving', 'thck_effective_min'),
    ('parameters', 'eigenconstant'): ('calving', 'eigenconstant'),
    ('parameters', 'tau_eigenconstant1'): ('calving', 'tau_eigenconstant1'),
    ('parameters', 'tau_eigenconstant2'): ('calving', 'tau_eigenconstant2'),
    ('parameters', 'stress_threshold'): ('calving', 'stress_threshold'),
    ('parameters', 'lateral_rate_min'): ('calving', 'lateral_rate_min'),
    ('parameters', 'effec_stress_min'): ('calving', 'effec_stress_min'),
    ('parameters', 'length_scale'): ('calving', 'calving_length_scale'),
    ('parameters', 'damage_constant'): ('calving', 'damage_constant'),
    ('parameters', 'taumax_cliff'): ('calving', 'taumax_cliff'),
    ('parameters', 'cliff_timescale'): ('calving', 'cliff_timescale'),
    ('parameters', 'calving_front_x'): ('calving', 'calving_front_x'),
    ('parameters', 'calving_front_y'): ('calving', 'calving_front_y'),
    ('parameters', 'calving_front_radius'): ('calving', 'calving_front_radius'),
    ('parameters', 'cf_advance_retreat_amplitude'): ('calving', 'cf_advance_retreat_amplitude'),
    ('parameters', 'cf_advance_retreat_period'): ('calving', 'cf_advance_retreat_period'),
    ('parameters', 'f_ground_threshold'): ('calving', 'f_ground_threshold'),
}

# Order of the new sections when they are created.
SECTION_ORDER = [
    'grid', 'time', 'options', 'numerics', 'constants', 'initialization',
    'velocity', 'velocity_solver', 'advection', 'flow_law', 'temperature',
    'basal_physics', 'basal_hydro', 'grounding_line', 'basal_melt', 'ocean_data',
    'climate', 'calving', 'removal', 'lateral_melt', 'inversion', 'geothermal',
    'isostasy', 'glaciers', 'diagnostics',
]

# CESM user_nl_cism: old key -> new key.  Kept separate from MOVES because a few key names occur in
# more than one CISM section, while the CESM namelist has one entry per name.
USER_NL_RENAMES = {
    # [calving]
    'marine_margin': 'calving_law',
    'which_ho_calving_front': 'subgrid_calving_front',
    'which_ho_calvingmip_domain': 'calvingmip_domain',
    'length_scale': 'calving_length_scale',
}

# Match the CISM reader (glimmer_config.F90): a section header is a line starting with '[' in
# column 1, named by the text before the first ']' (anything after it is ignored); a key line
# starts in column 1 (lines starting with a space are ignored by CISM, so they are never moved).
SECTION_RE = re.compile(r'^\[([^\]]*)\]')
KEY_RE = re.compile(r'^()([A-Za-z_][A-Za-z0-9_]*)(\s*[=:].*)$')
COMMENT_RE = re.compile(r'^\s*[#;!]')
NAMELIST_KEY_RE = re.compile(r'^(\s*)([A-Za-z_][A-Za-z0-9_]*)(\s*=.*)$')


class Section:
    def __init__(self, name, header):
        self.name = name        # None for lines before the first section
        self.header = header    # the original header line
        self.lines = []


def parse(text):
    """Split a config file into sections, keeping every line."""
    sections = [Section(None, None)]
    for line in text.split('\n'):
        m = SECTION_RE.match(line)
        if m:
            sections.append(Section(m.group(1).strip(), line))
        else:
            sections[-1].lines.append(line)
    return sections


def render(sections):
    out = []
    for s in sections:
        if s.header is not None:
            out.append(s.header)
        out.extend(s.lines)
    return '\n'.join(out)


def has_keys(section):
    return any(KEY_RE.match(l) and not COMMENT_RE.match(l) for l in section.lines)


def is_glide(sections):
    """True if [options] sets dycore = 0 (Glide)."""
    for s in sections:
        if s.name == 'options':
            for line in s.lines:
                m = KEY_RE.match(line)
                if m and not COMMENT_RE.match(line) and m.group(2) == 'dycore':
                    value = re.split(r'[#;!]', m.group(3).lstrip(' \t=:'))[0].strip()
                    return value == '0'
    return False


def convert_config(text, moves=None):
    """Return (new_text, list of changes) for a CISM config file.  Glide configs are returned unchanged."""
    moves = MOVES if moves is None else moves
    sections = parse(text)
    if is_glide(sections):
        return text, []
    changes = []
    moved = {}   # new section -> list of lines to add
    emptied = set()   # sections that lost keys

    for s in sections:
        if s.name is None:
            continue
        kept = []
        for line in s.lines:
            m = KEY_RE.match(line)
            if m and not COMMENT_RE.match(line) and (s.name, m.group(2)) in moves:
                new_sec, new_key = moves[(s.name, m.group(2))]
                # take the comment lines directly above this key along with it
                carried = []
                while kept and kept[-1].strip() and COMMENT_RE.match(kept[-1]):
                    carried.insert(0, kept.pop())
                moved.setdefault(new_sec, []).extend(carried + [new_key + m.group(3)])
                emptied.add(s.name)
                changes.append(f'[{s.name}] {m.group(2)} -> [{new_sec}] {new_key}')
            else:
                kept.append(line)
        s.lines = kept

    for new_sec, lines in moved.items():
        target = next((s for s in sections if s.name == new_sec), None)
        if target is None:
            target = Section(new_sec, f'[{new_sec}]')
            target.lines = ['']
            # insert before the first [CF ...] section, keeping SECTION_ORDER among new sections
            pos = next((i for i, s in enumerate(sections) if s.name and s.name.startswith('CF ')), len(sections))
            sections.insert(pos, target)
            changes.append(f'new section [{new_sec}]')
        existing = {KEY_RE.match(l).group(2) for l in target.lines if KEY_RE.match(l) and not COMMENT_RE.match(l)}
        for l in lines:
            m = KEY_RE.match(l)
            if m and not COMMENT_RE.match(l) and m.group(2) in existing:
                raise ValueError(f'[{new_sec}] already has a key {m.group(2)}')
        # add before the trailing blank lines of the target section
        n = len(target.lines)
        while n > 0 and not target.lines[n - 1].strip():
            n -= 1
        target.lines[n:n] = lines

    # order newly created sections by SECTION_ORDER within the block where they were inserted
    def order_key(s):
        return SECTION_ORDER.index(s.name) if s.name in SECTION_ORDER else len(SECTION_ORDER)
    new_names = {c[len('new section ['):-1] for c in changes if c.startswith('new section')}
    i = 0
    while i < len(sections):
        if sections[i].name in new_names:
            j = i
            while j < len(sections) and sections[j].name in new_names:
                j += 1
            sections[i:j] = sorted(sections[i:j], key=order_key)
            i = j
        else:
            i += 1

    # remove old sections that have no keys left
    for s in list(sections):
        if s.name in emptied and s.name not in moved and not has_keys(s):
            sections.remove(s)
            changes.append(f'removed empty section [{s.name}]')

    return render(sections), changes


def convert_user_nl(text, renames=None):
    """Return (new_text, list of changes) for a CESM user_nl_cism file."""
    renames = USER_NL_RENAMES if renames is None else renames
    out, changes = [], []
    for line in text.split('\n'):
        m = NAMELIST_KEY_RE.match(line)
        if m and not COMMENT_RE.match(line) and m.group(2) in renames:
            out.append(m.group(1) + renames[m.group(2)] + m.group(3))
            changes.append(f'{m.group(2)} -> {renames[m.group(2)]}')
        else:
            out.append(line)
    return '\n'.join(out), changes


def is_config(text):
    return any(SECTION_RE.match(l) for l in text.split('\n'))


def self_test():
    moves = {('options', 'marine_margin'): ('calving', 'calving_law'),
             ('parameters', 'calving_timescale'): ('calving', 'calving_timescale'),
             ('ho_options', 'which_ho_babc'): ('basal_physics', 'sliding_law')}
    old = '\n'.join([
        '[grid]', 'ewn = 31', '',
        '[options]', 'dycore = 2', '# calving: 4 = thickness threshold', 'marine_margin = 4   # see docs', '',
        '[ho_options]', 'which_ho_babc = 9', '',
        '[parameters]', 'calving_timescale = 1.0', 'ice_limit = 1.', '',
        '[CF output]', 'name = out.nc', ''])
    new, changes = convert_config(old, moves)
    assert '[calving]\n# calving: 4 = thickness threshold\ncalving_law = 4   # see docs\ncalving_timescale = 1.0\n\n' in new, new
    assert '[basal_physics]\nsliding_law = 9\n\n' in new, new
    assert '[ho_options]' not in new, new                      # emptied section removed
    assert new.index('[basal_physics]') < new.index('[calving]') < new.index('[CF output]'), new
    assert 'dycore = 2' in new and 'ice_limit = 1.' in new
    again, changes2 = convert_config(new, moves)                # idempotent
    assert again == new and not changes2, changes2
    nl, c = convert_user_nl("which_ho_babc = 9\n! comment\nice_limit = 1.\n", {'which_ho_babc': 'sliding_law'})
    assert nl == "sliding_law = 9\n! comment\nice_limit = 1.\n" and c == ['which_ho_babc -> sliding_law'], nl
    assert convert_config(old, {})[0] == old                    # empty table: no change
    commented = old.replace('[ho_options]', '[ho_options]   # Higher order options').replace(
        'dycore = 2', 'dycore = 2\n marine_margin = 1')
    new2, _ = convert_config(commented, moves)                  # header with a trailing comment
    assert '[basal_physics]\nsliding_law = 9\n' in new2 and '[ho_options]' not in new2, new2
    assert ' marine_margin = 1' in new2                         # indented key (ignored by CISM) left alone
    glide = old.replace('dycore = 2', 'dycore = 0   # Glide')
    assert convert_config(glide, moves) == (glide, [])          # Glide config: no change
    print('self-test passed')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('files', nargs='*')
    ap.add_argument('-o', '--output')
    ap.add_argument('--in-place', action='store_true')
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--user-nl', action='store_true', help='treat the input as a CESM user_nl_cism file')
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args()
    if a.self_test:
        self_test()
        return 0
    if not a.files:
        ap.error('no input files')
    if a.output and len(a.files) != 1:
        ap.error('-o needs exactly one input file')
    any_change = False
    for f in a.files:
        text = open(f).read()
        if a.user_nl or not is_config(text):
            new, changes = convert_user_nl(text)
        else:
            new, changes = convert_config(text)
        any_change |= bool(changes)
        glide = not (a.user_nl or not is_config(text)) and is_glide(parse(text))
        if a.check:
            print(f'{f}: ' + ('Glide config, left unchanged' if glide else '; '.join(changes) if changes else 'no changes'))
            continue
        dest = a.output or (f if a.in_place else f + '.converted')
        if changes or dest != f:
            open(dest, 'w').write(new)
        print(f'{f}: {len(changes)} changes' + (f', written to {dest}' if dest != f else ''))
    return 1 if (a.check and any_change) else 0


if __name__ == '__main__':
    sys.exit(main())
