#!/usr/bin/env python3
"""T_e at the two ends of each plate at every step: a pinned corner node repeats the same digits at every step.
  python3 util/plate_end_pin.py NEW='run/postproc/boundary_quantities*.dat'
Samples sit at (m-1)/n_elm_pts of each element, so the first sample of a plate is a node and the last one is
inside the final element; the pinned end therefore repeats exactly, the other end stays nearly constant."""
import math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sheath_row_check import load, MU0, AMU

label, pat = sys.argv[1].split("=", 1)
r = load(label, pat, "k", math.sqrt(MU0 * 0.76464396e20 * 2.01410177811 * AMU), 4.23130791805)   # Cathey case units
for pn in ("inner", "outer"):
    s = r.plates[pn]["s"] * 100
    print(f"\n{pn}: s[cm]  {s[0]:7.2f} {s[1]:7.2f}   ...   {s[-2]:7.2f} {s[-1]:7.2f}")
    for i, st in enumerate(r.steps):
        Te = r.avg("T_e", pn, i)
        print(f"  step {st:5d}  {Te[0]:8.4f} {Te[1]:8.4f}   ...   {Te[-2]:8.4f} {Te[-1]:8.4f}  eV")
