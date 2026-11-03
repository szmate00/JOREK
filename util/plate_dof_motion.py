#!/usr/bin/env python3
"""Are the plate DOFs pinned or free?  Relative change of each plate quantity between the first and the
last step of a boundary_quantities file, on the plate interior (end elements excluded).

  python3 util/plate_dof_motion.py B='runB/postproc/boundary_quantities*.dat' [--F0 4.2313]

A quantity constant to < 1e-6 while T_e moves is a pinned DOF (Dirichlet row with zero right-hand side);
a quantity that moves with T_e is free or row-controlled.  zj is read back from Jn_par / Bnorm.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sheath_row_check import load                                                     # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+", help="LABEL=GLOB of boundary_quantities files")
ap.add_argument("--F0", type=float, default=4.23130791805)
a = ap.parse_args()

for spec in a.runs:
    label, pat = spec.split("=", 1)
    r = load(label, pat, "k", 1.0, a.F0)
    i0, i1 = 0, r.steps.size - 1
    for pn in r.plates:
        keep = r.row_mask(pn, i0)
        g = lambda n, i: np.asarray(r.avg(n, pn, i))[keep]
        rows = {
            "vpar": (g("vpar", i0), g("vpar", i1)),
            "zj (Jn_par/Bnorm)": (g("Jn_par", i0) / g("Bnorm", i0), g("Jn_par", i1) / g("Bnorm", i1)),
            "Phi": (g("Phi", i0), g("Phi", i1)),
            "Psi_N": (g("Psi_N", i0), g("Psi_N", i1)),
            "T_e": (g("T_e", i0), g("T_e", i1)),
            "ne": (g("ne", i0), g("ne", i1)),
        }
        print(f"\n{label} {pn}: steps {r.steps[i0]} -> {r.steps[i1]}, {keep.sum()} interior samples")
        print(f"  {'quantity':18s} {'max|dv|/max|v|':>15s} {'median|dv|/|v|':>15s}  verdict")
        for name, (v0, v1) in rows.items():
            scale = np.max(np.abs(v0)) + 1e-300
            dv = np.abs(v1 - v0)
            rel_max = float(np.max(dv) / scale)
            rel_med = float(np.median(dv / (np.abs(v0) + 1e-300)))
            verdict = "PINNED (constant)" if rel_max < 1e-6 else ("moves" if rel_max > 1e-2 else "almost constant")
            print(f"  {name:18s} {rel_max:15.3e} {rel_med:15.3e}  {verdict}")
