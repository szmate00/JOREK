#!/usr/bin/env python3
"""Catch `use <module>, only:` lists that are missing a name the file references.

mod_boundary_conditions.f90 cannot be compiled by the unit-test harness (it pulls in
MPI, vacuum, tr_module, ...), so an omission from its only-list is invisible to the
tests and only shows up on the cluster. That happened for mach1_weak. This is a cheap
static stand-in: for every file that imports phys_module with an only-list, any
phys_module name used in the body must be imported or declared locally.
"""
import io, re, sys, os

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC  = os.path.join(REPO, 'models', 'phys_module.f90')
TARGETS = [
    'models/model600/mod_boundary_conditions.f90',
    'models/model600/mod_boundary_matrix_open.f90',
    'models/model600/mod_elt_matrix_fft.f90',
]
DECL = re.compile(r'^\s*(real|integer|logical|character|type)\b.*::\s*(.*)$', re.I)

def read(p):
    return io.open(p, encoding='utf-8', errors='surrogateescape').read()

def declared_names(text):
    """Names declared in phys_module's specification part."""
    out = set()
    for line in text.split('\n'):
        s = line.split('!')[0]
        m = DECL.match(s)
        if not m:
            continue
        for item in m.group(2).split(','):
            item = item.strip()
            item = item.split('=')[0].split('(')[0].strip()
            if re.fullmatch(r'[A-Za-z]\w*', item or ''):
                out.add(item.lower())
    return out

def full_use_stmt(lines, i):
    buf = lines[i]
    while buf.rstrip().endswith('&'):
        i += 1
        buf = buf.rstrip().rstrip('&') + ' ' + lines[i]
    return buf, i

def main():
    phys = declared_names(read(SRC))
    bad = 0
    for rel in TARGETS:
        path = os.path.join(REPO, rel)
        if not os.path.exists(path):
            continue
        text = read(path)
        lines = text.split('\n')
        imported, has_only, local, body = set(), False, set(), []
        i = 0
        while i < len(lines):
            s = lines[i].split('!')[0]
            if re.match(r'\s*use\s+phys_module\s*,\s*only\s*:', s, re.I):
                stmt, i = full_use_stmt(lines, i)
                has_only = True
                for n in stmt.split(':', 1)[1].split(','):
                    n = n.split('=>')[-1].strip().lower()
                    if n:
                        imported.add(n)
            elif re.match(r'\s*use\s+phys_module\s*$', s, re.I):
                has_only = False
                imported |= phys          # whole-module import: everything is visible
            else:
                m = DECL.match(s)
                if m:
                    for item in m.group(2).split(','):
                        item = item.strip().split('=')[0].split('(')[0].strip()
                        if re.fullmatch(r'[A-Za-z]\w*', item or ''):
                            local.add(item.lower())
                body.append(s)
            i += 1
        if not has_only:
            print('  %-52s whole-module import, nothing to check' % rel)
            continue
        # Derived-type components share names with phys_module variables
        # (bcs(i)%dirichlet%u, %mach1, %rho ...). Strip every `%component`
        # before tokenising, or they all read as missing imports.
        KEYWORDS = {'type', 'if', 'then', 'else', 'endif', 'do', 'enddo', 'end',
                    'call', 'return', 'use', 'only', 'implicit', 'none', 'true',
                    'false', 'and', 'or', 'not'}
        used = set()
        for s in body:
            # identifiers inside string literals are not references
            s = re.sub(r"'[^']*'", ' ', s)
            s = re.sub(r'"[^"]*"', ' ', s)
            s = re.sub(r'%\s*[A-Za-z]\w*', ' ', s)
            for tok in re.findall(r'[A-Za-z]\w*', s):
                used.add(tok.lower())
        used -= KEYWORDS
        missing = sorted((used & phys) - imported - local)
        if missing:
            bad += 1
            print('  %-52s MISSING: %s' % (rel, ', '.join(missing)))
        else:
            print('  %-52s only-list complete (%d imported)' % (rel, len(imported)))
    bad += extent_check()
    if bad:
        print('IMPORT CHECK FAIL: %d problem(s)' % bad)
        return 1
    print('IMPORT CHECK PASS')
    return 0


def extent_check():
    """Whole-array assignments between arrays of different declared extent.

    mod_boundary_conditions declares fd_* with FD_NT = 12 while the diagnostics
    module declares fw_* with FW_NT = 9, so `fd_loc_rho = fw_mom_max` is a shape
    mismatch that ifort rejects. The harness cannot compile that file, so this is
    caught statically instead.
    """
    bad = 0
    files = ['models/model600/mod_boundary_conditions.f90',
             'models/model600/mod_floating_transport_diag.f90']
    # resolve named integer parameters, so an extent written as 9 and one written
    # as FW_NT (= 9) compare equal instead of reporting a false mismatch
    const = {}
    for rel in files:
        path = os.path.join(REPO, rel)
        if not os.path.exists(path):
            continue
        for m in re.finditer(r'integer\s*,\s*parameter[^:]*::\s*([A-Za-z]\w*)\s*=\s*(\d+)',
                             read(path), re.I):
            const[m.group(1).lower()] = m.group(2)

    def resolve(e):
        return const.get(e.lower(), e)

    extent = {}
    for rel in files:
        path = os.path.join(REPO, rel)
        if not os.path.exists(path):
            continue
        for line in read(path).split('\n'):
            s = line.split('!')[0]
            if not re.match(r'\s*(real|integer)', s, re.I) or '::' not in s:
                continue
            for item in s.split('::', 1)[1].split(','):
                m = re.match(r'\s*([A-Za-z]\w*)\s*\(\s*([A-Za-z]\w*|\d+)\s*\)', item)
                if m:
                    extent.setdefault(m.group(1).lower(), m.group(2))
    for rel in files:
        path = os.path.join(REPO, rel)
        if not os.path.exists(path):
            continue
        for n, line in enumerate(read(path).split('\n'), 1):
            s = line.split('!')[0].strip()
            m = re.fullmatch(r'([A-Za-z]\w*)\s*=\s*([A-Za-z]\w*)', s)
            if not m:
                continue
            l, r = m.group(1).lower(), m.group(2).lower()
            if l in extent and r in extent and resolve(extent[l]) != resolve(extent[r]):
                print('  %s:%d  EXTENT MISMATCH: %s(%s) = %s(%s)'
                      % (rel, n, l, extent[l], r, extent[r]))
                bad += 1
    if not bad:
        print('  %-52s no cross-extent whole-array assignments' % 'array extents')
    return bad

if __name__ == '__main__':
    sys.exit(main())
