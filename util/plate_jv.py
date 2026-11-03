#!/usr/bin/env python3
"""Sheath characteristic of the plates, in the style of verif_sheath.pdf (util/paper_verification_figs.py):
j_n/j_sat against e(Phi - Phi_fl)/T_e, inner plate blue circles, outer plate orange triangles, first run filled,
second run open, the law 1 - exp(-x) in black. End elements (two per plate end, as in sheath_row_check) are left
out unless --show-ends, which draws them in red.
  python3 util/plate_jv.py A='runA/postproc/boundary_quantities*.dat' [B='runB/...'] [--steps 1300 1300]
                           [--lambda 3.0] [--per-plane] [--show-ends] [-o plate_jv.pdf]
One step per run (--steps, default: the last step of each run); the legend shows the SI time of that step
after the label (an empty label, "=GLOB", gives the time alone; --no-time gives the label alone). Samples further than --tol from the law are
printed with plate, s, T_e and whether they sit in an end element. --per-plane (3D runs) uses one sample per
toroidal plane instead of the plane mean, so the ratio is never taken of averaged quantities."""
import argparse, math, os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sheath_row_check import load, MU0, AMU

C_PLATE = {"inner": "#00B1EA", "outer": "#EF7C00"}
MARK = {"inner": "o", "outer": "^"}
INK, MUTED, C_END = "k", "0.5", "#D03B3B"

ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+", help="LABEL=GLOB (one or two runs)")
ap.add_argument("--steps", nargs="*", type=int, help="one restart step per run (default: last)")
ap.add_argument("--lambda", dest="lam", type=float, default=3.0, help="sheath_Lambda of the run")
ap.add_argument("--tol", type=float, default=0.2, help="print samples with |j/j_sat - f(x)| above this")
ap.add_argument("--per-plane", action="store_true", help="3D runs: one sample per plane instead of the plane mean")
ap.add_argument("--show-ends", action="store_true", help="draw the end-element samples in red")
ap.add_argument("--central-density", type=float, default=0.76464396)
ap.add_argument("--F0", type=float, default=4.23130791805)
ap.add_argument("-o", default="plate_jv.pdf")
ap.add_argument("--no-time", action="store_true", help="legend without the SI time of the sample")
a = ap.parse_args()
t_unit = math.sqrt(MU0 * a.central_density * 1e20 * 2.01410177811 * AMU)
plt.rcdefaults()                       # default matplotlib style, as verif_sheath.pdf (the loader modules set their own)

fig, ax = plt.subplots(figsize=(4.6, 3.6), constrained_layout=True)
x_all = []
for k, item in enumerate(a.runs[:2]):
    label, pat = item.split("=", 1)
    r = load(label, pat, "k", t_unit, a.F0)
    i = int(np.argmin(np.abs(r.steps - a.steps[k]))) if a.steps and k < len(a.steps) else r.steps.size - 1
    print(f"{label}: step {r.steps[i]}, t = {r.t[i]:.3f} ms")
    tag = label if a.no_time else (f"{label}, " if label else "") + f"t = {r.t[i]:.2f} ms"
    for pn in ("inner", "outer"):
        s = r.plates[pn]["s"] * 100
        keep = r.row_mask(pn, i)                                         # False in the end elements
        if a.per_plane:
            j, js, Phi, Te = (r.q(n, pn, i).ravel() for n in ("Jn_par", "Jn_sat", "Phi", "T_e"))
            nphi = r.q("T_e", pn, i).shape[0]
            keep, s = np.tile(keep, nphi), np.tile(s, nphi)
        else:
            j, js, Phi, Te = (r.avg(n, pn, i) for n in ("Jn_par", "Jn_sat", "Phi", "T_e"))
        ok = (Te > 0.1) & (js > 0.02 * js.max())
        x, y = Phi / Te - a.lam, j / js
        f = 1.0 - np.exp(-np.maximum(x, -a.lam))
        col = C_PLATE[pn]
        kw = dict(ls="none", marker=MARK[pn], ms=6 if pn == "inner" else 7, markeredgewidth=1.2)
        ax.plot(x[ok & keep], y[ok & keep], color=col, markerfacecolor=col if k == 0 else "none",
                markeredgecolor=col, zorder=3 + k, label=f"{tag}, {pn} plate", **kw)
        if a.show_ends:
            ax.plot(x[ok & ~keep], y[ok & ~keep], color=C_END, markerfacecolor=C_END if k == 0 else "none",
                    markeredgecolor=C_END, zorder=2, label="end elements" if (k == 0 and pn == "inner") else None, **kw)
        x_all.append(x[ok & keep])
        dev = np.abs(y - f)
        print(f"  {label}, {pn} plate: {int((ok & keep).sum())} interior points, rms deviation "
              f"{np.sqrt(np.mean(dev[ok & keep] ** 2)):.3f}, largest {dev[ok & keep].max():.3f}")
        for m in np.flatnonzero(ok & (dev > a.tol)):
            print(f"    step {r.steps[i]:5d} s = {s[m]:6.2f} cm  x = {x[m]:6.2f}  j/j_sat = {y[m]:6.2f}  law = {f[m]:5.2f}"
                  f"  T_e = {Te[m]:5.1f} eV  {'END ELEMENT' if not keep[m] else 'interior'}")
x_all = np.concatenate(x_all)
xs = np.linspace(min(x_all.min(), -1.0) - 0.2, max(x_all.max(), 2.0) + 0.2, 400)
ax.plot(xs, 1.0 - np.exp(-np.maximum(xs, -a.lam)), color=INK, lw=1.2, zorder=0,
        label=r"$1-\exp[-e(\Phi-\Phi_\mathrm{fl})/T_e]$")
ax.axhline(0.0, color=MUTED, lw=0.6)
ax.axvline(0.0, color=MUTED, lw=0.6)
ax.set_xlabel(r"$e\,(\Phi-\Phi_\mathrm{fl})/T_e$ on the plate", fontsize=14)
ax.set_ylabel(r"$j_n/j_\mathrm{sat}$", fontsize=14)
ax.legend(loc="lower right", frameon=False)
fig.savefig(a.o)
print(f"wrote {a.o}")
plt.show()
