#!/usr/bin/env python3
"""The 1D thermoelectric prediction of jorek2_SOLcurrent against the plate current of a sheath_j run.
  python3 util/compare_solcurrent.py --sol Spitzer=runB/SOLcurrent_spitzer/step001300.dat eta=runB/SOLcurrent_eta/step001300.dat \\
      --boundary 'runB/postproc/boundary_quantities*.dat' --step 1300 [--plate-zmax -0.5] [-o solcurrent_vs_plates.pdf]
jorek2_SOLcurrent traces every open field line from a wall point to the other target and solves Stangeby 17.29 on it
(Staebler-Hinton circuit: both sheaths, 0.71 thermal force, electron-pressure integral, Spitzer-Haerm conductivity
averaged along the line). Its I_loc is the current the tube would carry if it closed plate to plate with no
cross-field current; its sign is positive for current OUT of the wall at the start point. Here it is converted to
the plate convention of the boundary output (Jn_par > 0 = ions collected = current INTO the wall): j_plate = -j_wall.
Every start point is assigned to the nearest plate sample of the boundary file to get its Psi_N and side.
Figure: per plate, the diagnostic's thermoelectric j against our Jn_par on the +-(Psi_N-1) axis, plus their sums.
Table: totals per plate, the open-circuit voltage pieces (Lambda dTe, pe term) and the Spitzer tube resistance per
tube band, next to our measured Phi_out - Phi_in."""
import argparse, math, os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sheath_row_check import load, MU0, AMU

COLS = ["R_b", "Z_b", "Phi_b", "R_e", "Z_e", "Phi_e", "ne_b", "ne_e", "Te_b", "Te_e", "L", "av_sigma_par",
        "dP", "pe_term", "j_hat_par", "j_par", "bdotn", "j_wall", "dl", "I_loc"]
C = {"inner": "#00B1EA", "outer": "#EF7C00"}

ap = argparse.ArgumentParser()
ap.add_argument("--sol", required=True, nargs="+", help="SOLcurrent/stepXXXXXX.dat, or LABEL=path (several: e.g. Spitzer=... eta=...)")
ap.add_argument("--boundary", required=True, help="boundary_quantities file(s) of the same run (LABEL not needed)")
ap.add_argument("--step", type=int, help="boundary step (default: last)")
ap.add_argument("--lambda", dest="lam", type=float, default=3.0)
ap.add_argument("--central-density", type=float, default=0.76464396)
ap.add_argument("--F0", type=float, default=4.23130791805)
ap.add_argument("--match", type=float, default=0.01, help="max distance [m] start point to plate sample")
ap.add_argument("-o", default="solcurrent_vs_plates.pdf")
a = ap.parse_args()
plt.rcdefaults()

sols = []
for item in a.sol:
    label, _, path = item.rpartition("=") if "=" in item else ("1D thermoelectric", "", item)
    d = np.genfromtxt(path)
    if d.ndim != 2 or d.shape[1] != len(COLS):
        sys.exit(f"{path}: expected {len(COLS)} columns, got {d.shape}")
    sols.append((label, d, {n: d[:, k] for k, n in enumerate(COLS)}))
LS = ["--", ":", "-."]
t_unit = math.sqrt(MU0 * a.central_density * 1e20 * 2.01410177811 * AMU)
r = load("run", a.boundary, "k", t_unit, a.F0)
i = int(np.argmin(np.abs(r.steps - a.step))) if a.step is not None else r.steps.size - 1
print(f"boundary step {r.steps[i]}, t = {r.t[i]:.3f} ms")
for label, d, s in sols:
    print(f"  {label}: {d.shape[0]} start points, sum I_loc = {s['I_loc'].sum() * 1e-3:+.3f} kA (~0 by construction)")

fig, ax = plt.subplots(2, 1, figsize=(6.5, 7), sharex=True, constrained_layout=True)
tot = {}
for pn in ("inner", "outer"):
    p = r.plates[pn]
    Rp, Zp = r.q("R", pn, i)[0], r.q("Z", pn, i)[0]
    psi, jn, Te, Phi, ring = (r.avg(n, pn, i) for n in ("Psi_N", "Jn_par", "T_e", "Phi", "ring")) if "ring" in r.c else \
        (*(r.avg(n, pn, i) for n in ("Psi_N", "Jn_par", "T_e", "Phi")), p["ring"])
    k0 = int(np.argmin(psi))
    x = np.where(p["s"] >= p["s"][k0], 1.0, -1.0) * np.abs(psi - 1.0)
    I_us = float(np.sum(jn * ring))
    ax[0].plot(x, jn * 1e-3, color=C[pn], lw=1.8, label=f"{pn}: sheath row (Jn_par)")
    tot[pn] = [I_us]
    for m, (label, d, s) in enumerate(sols):
        # assign SOLcurrent start points to this plate's samples
        dist = np.hypot(s["R_b"][:, None] - Rp[None, :], s["Z_b"][:, None] - Zp[None, :])
        near = np.argmin(dist, axis=1)
        ok = dist[np.arange(d.shape[0]), near] < a.match
        j_te = -s["j_wall"][ok]                                       # plate convention: + = into the wall
        x_te = x[near[ok]]
        I_te = -s["I_loc"][ok].sum()
        tot[pn].append(I_te)
        o = np.argsort(x_te)
        ax[0].plot(x_te[o], j_te[o] * 1e-3, color=C[pn], lw=1.2, ls=LS[m % 3], label=f"{pn}: {label} (jorek2_SOLcurrent)")
        # open-circuit pieces per start point, for the table
        Th, Tc = np.maximum(s["Te_b"][ok], s["Te_e"][ok]), np.minimum(s["Te_b"][ok], s["Te_e"][ok])
        V_lam = a.lam * (Th - Tc)                                      # V
        V_pe = s["pe_term"][ok] * Th                                   # pe_term = integral/(e Th) -> V
        R_sp = s["L"][ok] / np.maximum(s["av_sigma_par"][ok], 1e-30)  # Ohm m^2 (per unit area)
        print(f"\n{pn} plate, {label}: {ok.sum()} tubes matched; I_thermo(1D) = {I_te * 1e-3:+.3f} kA, I(sheath row) = {I_us * 1e-3:+.3f} kA")
        print(f"  {'x=+-(PsiN-1)':>13s} {'Th[eV]':>7s} {'Tc[eV]':>7s} {'L[m]':>6s} {'Lam dT [V]':>10s} {'pe term [V]':>11s} {'0.71dT [V]':>10s}"
              f" {'R_tube [Ohm m2]':>15s} {'j_1D':>8s} {'j_row':>8s}  [kA/m2]")
        for k in np.argsort(x_te)[:: max(1, ok.sum() // 12)]:
            print(f"  {x_te[k]:13.4f} {Th[k]:7.1f} {Tc[k]:7.1f} {s['L'][ok][k]:6.1f} {V_lam[k]:10.1f} {V_pe[k]:11.1f} {0.71 * (Th[k] - Tc[k]):10.1f}"
                  f" {R_sp[k]:15.2e} {j_te[k] * 1e-3:8.2f} {jn[near[ok]][k] * 1e-3:8.2f}")
print(f"\nboth plates: sheath row sum {sum(v[0] for v in tot.values()) * 1e-3:+.3f} kA; "
      + "; ".join(f"{label} sum {sum(v[1 + m] for v in tot.values()) * 1e-3:+.3f} kA" for m, (label, _, _) in enumerate(sols))
      + " (1D sums are zero by construction)")
ax[0].axhline(0, color="0.5", lw=0.6); ax[0].axvline(0, color="0.5", lw=0.6)
ax[0].set_ylabel(r"$j_n$ into the plate [kA/m$^2$]")
ax[0].legend(frameon=False, fontsize=8)
# bottom: the two plates' sums on a common grid for both models
xg = np.linspace(-0.03, 0.14, 300)
for name, ls, sol in [("sheath row", "-", None)] + [(label, LS[m % 3], (d, s)) for m, (label, d, s) in enumerate(sols)]:
    tot_curve = np.zeros_like(xg)
    for pn in ("inner", "outer"):
        p = r.plates[pn]
        psi, jn = r.avg("Psi_N", pn, i), r.avg("Jn_par", pn, i)
        k0 = int(np.argmin(psi)); x = np.where(p["s"] >= p["s"][k0], 1.0, -1.0) * np.abs(psi - 1.0)
        if sol is None:
            y, xx = jn, x
        else:
            d, s = sol
            Rp, Zp = r.q("R", pn, i)[0], r.q("Z", pn, i)[0]
            dist = np.hypot(s["R_b"][:, None] - Rp[None, :], s["Z_b"][:, None] - Zp[None, :])
            near = np.argmin(dist, axis=1); ok = dist[np.arange(d.shape[0]), near] < a.match
            y, xx = -s["j_wall"][ok], x[near[ok]]
        o = np.argsort(xx)
        tot_curve += np.interp(xg, xx[o], y[o], left=0, right=0)
    ax[1].plot(xg, tot_curve * 1e-3, color="k", ls=ls, label=f"inner + outer, {name}")
ax[1].axhline(0, color="0.5", lw=0.6); ax[1].axvline(0, color="0.5", lw=0.6)
ax[1].set_ylabel(r"$j_{n,\mathrm{in}} + j_{n,\mathrm{out}}$ [kA/m$^2$]")
ax[1].set_xlabel(r"$\pm(\Psi_N - 1)$ at the plate (left: private flux, right: SOL)")
ax[1].legend(frameon=False, fontsize=8)
fig.savefig(a.o)
print(f"wrote {a.o}")
plt.show()
