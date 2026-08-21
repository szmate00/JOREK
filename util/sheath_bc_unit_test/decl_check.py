#!/usr/bin/env python3
"""Duplicate-declaration check for the files the stub compile cannot reach.

Fortran is case insensitive, so a new local called sh_C collides with an existing sh_c and the
build fails with "This name has already been assigned a data type" - but only after a full compile
on the cluster. This catches it in a second. Scopes are per subroutine/function, and initialiser
values are stripped so that 0.d0 is not mistaken for a name. Run from the repo root:

    python3 util/sheath_bc_unit_test/decl_check.py
"""
import re, sys, pathlib

FILES = ["models/model600/mod_boundary_conditions.f90",
         "models/model600/mod_boundary_matrix_open.f90",
         "models/model600/mod_sheath_bc.f90",
         "models/model600/mod_sheath_diag.f90"]

# names that legitimately appear in several declarations of one scope as array dimensions
IGNORE = {"n_var","n_order","n_degrees","n_degrees_1d","n_tor","n_plane","n_gauss",
          "n_vertex_max","max_bnd_types","kind","dimension","int_all","hsize_t","len"}

DECL = re.compile(r"\s*(real\s*\*\s*8|real|integer|logical|character|type\s*\()", re.I)
SCOPE = re.compile(r"^\s*(pure\s+|elemental\s+|recursive\s+)*(subroutine|function|program|module)\s+(\w+)", re.I)
ENDSCOPE = re.compile(r"^\s*end\s*(subroutine|function|program|module)\b", re.I)

def names_in(decl_rhs):
    """declared names only: strip initialisers and array specs, split at top-level commas"""
    out, depth, item = [], 0, ""
    for ch in decl_rhs.split("!")[0]:
        if ch in "([": depth += 1
        elif ch in ")]": depth -= 1
        if ch == "," and depth == 0:
            out.append(item); item = ""
        else:
            item += ch
    out.append(item)
    names = []
    for it in out:
        it = it.split("=")[0].strip()          # drop "= 0.d0"
        m = re.match(r"([A-Za-z_]\w*)", it)
        if m: names.append(m.group(1).lower())
    return names

status = 0
for f in FILES:
    p = pathlib.Path(f)
    if not p.exists():
        continue
    scope, decl, bad = "(file)", {}, {}
    for i, line in enumerate(p.read_text().split("\n")):
        if SCOPE.match(line) or ENDSCOPE.match(line):
            for k, v in decl.items():
                if len(v) > 1: bad[f"{scope}:{k}"] = v
            m = SCOPE.match(line)
            scope, decl = (m.group(3) if m else "(file)"), {}
            continue
        if not DECL.match(line) or "::" not in line:
            continue
        for nm in names_in(line.split("::", 1)[1]):
            if nm in IGNORE: continue
            decl.setdefault(nm, []).append(i + 1)
    for k, v in decl.items():
        if len(v) > 1: bad[f"{scope}:{k}"] = v
    if bad:
        status = 1
        print(f"FAIL  {f}")
        for k, v in sorted(bad.items()):
            print(f"        {k} declared at lines {v}")
    else:
        print(f"PASS  {f}")
sys.exit(status)
