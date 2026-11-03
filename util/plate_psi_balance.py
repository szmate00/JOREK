#!/usr/bin/env python3
"""Plate current profiles of both targets on one psi axis, and their sum.
  python3 util/plate_psi_balance.py A='runA/postproc/boundary_quantities*.dat' B='runB/...' [--steps 1300 1300]
                                    [-o plate_psi_balance.pdf]
x = sign * (Psi_N - 1): positive on the SOL side of the strike point of that step (the Psi_N minimum on the plate,
s > 0 = SOL), negative on the private-flux side (the postproc mirrors Psi_N there, so |Psi_N - 1| is used).
Top: j_n [kA/m^2] of the inner and the outer plate (Jn_par, > 0 = ions collected).
Bottom: the current per unit Psi_N, dI/dPsi_N = 2 pi R j_n dl/dPsi_N [kA], which is what a flux tube carries, and
the sum inner + outer interpolated onto a common grid: zero where the tube's current returns through the other
plate, non-zero where it crosses flux surfaces or leaves through other boundaries. Third panel: the plate potential Phi and the local floating potential Lambda*T_e of both plates: the
sign of j_n is the sign of Phi - Lambda*T_e at that point. Totals are printed."""
import argparse, math, os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sheath_row_check import load, MU0, AMU

C = {"inner": "#00B1EA", "outer": "#EF7C00", "sum": "k"}
ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+", help="LABEL=GLOB")
ap.add_argument("--steps", nargs="*", type=int, help="one restart step per run (default: last)")
ap.add_argument("--central-density", type=float, default=0.76464396)
ap.add_argument("--F0", type=float, default=4.23130791805)
ap.add_argument("--lambda", dest="lam", type=float, default=3.0, help="sheath_Lambda of the run")
ap.add_argument("-o", default="plate_psi_balance.pdf")
a = ap.parse_args()
t_unit = math.sqrt(MU0 * a.central_density * 1e20 * 2.01410177811 * AMU)
plt.rcdefaults()

fig, ax = plt.subplots(3, 1, figsize=(6.5, 9), constrained_layout=True, sharex=True)
LS = ["-", "--", ":", "-."]
for k, item in enumerate(a.runs):
    label, pat = item.split("=", 1)
    r = load(label, pat, "k", t_unit, a.F0)
    i = int(np.argmin(np.abs(r.steps - a.steps[k]))) if a.steps and k < len(a.steps) else r.steps.size - 1
    ls = LS[k % 4]
    prof = {}
    for pn in ("inner", "outer"):
        p = r.plates[pn]
        o = np.argsort(p["s"])
        s, ring = p["s"][o], p["ring"][o]
        psi, j = r.avg("Psi_N", pn, i)[o], r.avg("Jn_par", pn, i)[o] * 1e-3           # kA/m^2
        k0 = int(np.argmin(psi))
        x = np.where(s >= s[k0], 1.0, -1.0) * np.abs(psi - 1.0)
        dpsi = np.abs(np.gradient(psi))
        if 0 < k0 < psi.size - 1:                                                    # kink of Psi_N at the strike point
            dpsi[k0] = 0.5 * (dpsi[k0 - 1] + dpsi[k0 + 1])
        dIdpsi = j * ring / np.maximum(dpsi, 1e-12)                                  # kA per unit Psi_N
        prof[pn] = (x, j, dIdpsi, np.sum(j * ring))
        print(f"{label} {pn}: step {r.steps[i]}, t = {r.t[i]:.2f} ms, strike point s = {s[k0] * 100:+.2f} cm, "
              f"I = {np.sum(j * ring):+.3f} kA")
        ax[0].plot(x, j, color=C[pn], ls=ls, label=f"{label}, {pn}")
        Phi, Te = r.avg("Phi", pn, i)[o], r.avg("T_e", pn, i)[o]
        ax[2].plot(x, Phi, color=C[pn], ls=ls, label=f"{label}, {pn}: $\\Phi$")
        ax[2].plot(x, a.lam * Te, color=C[pn], ls=ls, lw=0.8, alpha=0.6, label=f"{label}, {pn}: $\\Lambda T_e$")
        ax[1].plot(x, dIdpsi, color=C[pn], ls=ls, label=f"{label}, {pn}")
    xg = np.linspace(max(prof["inner"][0].min(), prof["outer"][0].min()),
                     min(prof["inner"][0].max(), prof["outer"][0].max()), 400)
    tot = sum(np.interp(xg, *sorted_xy) for sorted_xy in
              ((prof[pn][0][np.argsort(prof[pn][0])], prof[pn][2][np.argsort(prof[pn][0])]) for pn in prof))
    ax[1].plot(xg, tot, color=C["sum"], ls=ls, lw=1.8, label=f"{label}, inner + outer")
    print(f"{label}: inner + outer = {prof['inner'][3] + prof['outer'][3]:+.3f} kA")
for x_ in ax:
    x_.axhline(0, color="0.5", lw=0.6)
    x_.axvline(0, color="0.5", lw=0.6)
ax[0].set_ylabel(r"$j_n$ [kA/m$^2$]")
ax[1].set_ylabel(r"$dI/d\Psi_N$ [kA]")
ax[2].set_ylabel(r"$\Phi$ and $\Lambda T_e$ [V]")
ax[2].set_xlabel(r"$\pm(\Psi_N-1)$ at the plate  (left: private flux, right: SOL)")
ax[2].legend(frameon=False, fontsize=7, ncol=2)
ax[0].legend(frameon=False, fontsize=8)
ax[1].legend(frameon=False, fontsize=8)
fig.savefig(a.o)
print(f"wrote {a.o}")
plt.show()
