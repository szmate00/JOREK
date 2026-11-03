#!/usr/bin/env python3
"""Does the sheath current row bind on the divertor plates? Axisymmetric (or any) boundary_quantities files.

  python3 util/sheath_row_check.py BASE='runB0/postproc/boundary_quantities*.dat' [GROUNDED=... TWOT=...] [--steps 1500 3000]

For every run and restart step, per plate (inner / outer, found from Bnorm and R as in compare_wall_runs):
  mismatch   rms(Jn_par - Jn_sheath) / rms(Jn_sheath) on the plate interior (the two end elements excluded)
             -> 0.01 when the row binds (inxflow), ~1 when it does not
  Phi/Te     median e*Phi/T_e on the interior (Lambda ~ 3 when floating, 7-9 in the broken AUG run)
  I_int      current collected by the plate interior [kA] (+ = ions)
  I_ends     current of the two end elements [kA]
  jsat/Bohm  median Jn_sat / (e n c_s |b.n|)  (1.3 on AUG with gamma = 5/3 in c_s)
  Mach       median |Mach_par| on the interior
Window 1: these against time for every run. Window 2: profiles along the plate at --steps (default: first and
last): Jn_par with Jn_sheath, Phi with Lambda*T_e, T_e, Mach; end elements shaded. Window 3: Jn_par vs Jn_sheath.
Sign convention: Jn > 0 = ions collected. Needs numpy, matplotlib; imports Run from util/compare_wall_runs.py.
"""
import argparse
import glob
import math
import os
import sys
import tempfile

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_wall_runs import Run, MU0, AMU, E_CHG                                   # noqa: E402

COLORS = ["0.35", "#d6392a", "#2a6fd6", "#2aa05a"]
LS = {"outer": "-", "inner": "--"}


def load(label, spec, color, t_unit, F0):
    files = []
    for pat in spec.split("+"):
        hits = sorted(glob.glob(pat))
        if not hits:
            sys.exit(f"{label}: no file matches {pat}")
        files += hits
    path = files[0]
    if len(files) > 1:
        tmp = tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False)
        for f in files:
            tmp.write(open(f).read() + "\n")
        tmp.close()
        path = tmp.name
    r = Run(path, label, color, t_unit)
    _, keep = np.unique(r.steps, return_index=True)
    r.steps, r.t, r.t_jorek, r.d = r.steps[keep], r.t[keep], r.t_jorek[keep], r.d[keep]
    r.F0 = F0
    print(f"{label}: {r.steps.size} steps {r.steps[0]}..{r.steps[-1]}, t = {r.t[0]:.4f}..{r.t[-1]:.4f} ms, "
          f"{r.d.shape[1]} plane(s)")
    return r


def plate_state(r, pn, i, mass, lam):
    """Per-plate scalars at step index i and the along-plate profiles used for the plots."""
    p = r.plates[pn]
    keep = r.row_mask(pn, i)
    g = lambda n: r.avg(n, pn, i)                                    # toroidal mean (one plane: the plane itself)
    j, jsh, js, Phi, Te, Ti, ne, bn = (g(n) for n in ("Jn_par", "Jn_sheath", "Jn_sat", "Phi", "T_e", "T_i", "ne", "Bnorm"))
    B = g("B_abs") if "B_abs" in r.c else np.hypot(r.F0 / p["R"], bn)
    mach = np.abs(g("Mach_par")) if "Mach_par" in r.c else np.full_like(j, np.nan)
    ok = keep & (Te > 0.1) & (js > 1e-3 * np.abs(js).max())
    cs = np.sqrt(5.0 / 3.0 * E_CHG * (Te + Ti) / mass)
    j_bohm = E_CHG * ne * cs * np.abs(bn) / B
    out = {
        "mismatch": np.sqrt(np.mean((j - jsh)[ok] ** 2)) / (np.sqrt(np.mean(jsh[ok] ** 2)) + 1e-30),
        "Phi/Te": float(np.median((Phi / Te)[ok])),
        "I_int": float(np.sum((j * p["ring"])[keep])) * 1e-3,
        "I_ends": float(np.sum((j * p["ring"])[~keep])) * 1e-3,
        "jsat/Bohm": float(np.median((js / j_bohm)[ok])),
        "Mach": float(np.nanmedian(mach[ok])),
        "angle": float(np.degrees(np.median(np.arcsin(np.clip(np.abs(bn) / B, 0, 1))[keep]))),
    }
    prof = dict(s=p["s"], keep=keep, j=j, jsh=jsh, js=js, Phi=Phi, Phi_fl=lam * Te, Te=Te, mach=mach)
    return out, prof


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", metavar="LABEL=PATTERN")
    ap.add_argument("--steps", nargs="*", type=int, default=None, help="restart steps for the profile window")
    ap.add_argument("--sheath-lambda", type=float, default=2.8)
    ap.add_argument("--F0", type=float, default=4.23130791805)
    ap.add_argument("--central-density", type=float, default=0.76464396)
    ap.add_argument("--central-mass", type=float, default=2.01410177811)
    args = ap.parse_args()
    mass = args.central_mass * AMU
    t_unit = math.sqrt(MU0 * args.central_density * 1e20 * mass)
    runs = []
    for k, item in enumerate(args.runs):
        label, _, spec = item.partition("=")
        runs.append(load(label, spec, COLORS[k % len(COLORS)], t_unit, args.F0))

    # boundary types along each plate (bnd_type_lo/hi = lowest/highest type among the element's nodes; 0 = none flagged)
    for r in runs:
        for pn in ("inner", "outer"):
            if "bnd_type_lo" in r.c and "bnd_type_hi" in r.c:
                lo = r.q("bnd_type_lo", pn, 0)[0].astype(int)
                hi = r.q("bnd_type_hi", pn, 0)[0].astype(int)
                seq = lo[::4].tolist()
                ok = all(v == 1 for v in seq[1:-1]) and seq[0] in (1, 3) and seq[-1] in (1, 3)
                print(f"{r.label} {pn}: boundary types per element {'as expected (1 inside, 1/3 at the ends)' if ok else seq}")

    # is psi pinned on the plate? Dirichlet psi is applied by the same node visits as every other boundary row, so a
    # drift of Psi_N along the plate between the first and the last step means those nodes receive no rows at all
    for r in runs:
        if r.steps.size < 2:
            continue
        for pn in ("inner", "outer"):
            keep = r.row_mask(pn)
            d0, d1 = r.q("Psi_N", pn, 0)[0], r.q("Psi_N", pn, r.steps.size - 1)[0]
            dR = np.abs(r.q("R", pn, r.steps.size - 1)[0] - r.q("R", pn, 0)[0]).max()
            print(f"{r.label} {pn}: |dPsi_N| steps {r.steps[0]}->{r.steps[-1]}: interior {np.abs(d1 - d0)[keep].max():.1e}, "
                  f"ends {np.abs(d1 - d0)[~keep].max():.1e}")

    keys = ["mismatch", "Phi/Te", "I_int", "I_ends", "jsat/Bohm", "Mach"]
    fig1, ax1 = plt.subplots(3, 2, figsize=(12, 9), sharex=True, constrained_layout=True, num="sheath row vs time")
    ax1 = ax1.ravel()
    print(f"\n{'run':8s}{'step':>7s}{'t [ms]':>9s} plate  " + "".join(f"{k:>10s}" for k in keys) + f"{'angle':>7s}   [kA for I, deg]")
    series = {}
    for r in runs:
        for pn in ("inner", "outer"):
            S = {k: [] for k in keys}
            for i in range(r.steps.size):
                st, _ = plate_state(r, pn, i, mass, args.sheath_lambda)
                for k in keys:
                    S[k].append(st[k])
                if i in (0, r.steps.size - 1):
                    print(f"{r.label:8s}{r.steps[i]:7d}{r.t[i]:9.4f} {pn:6s}" + "".join(f"{st[k]:10.3f}" for k in keys) + f"{st['angle']:7.2f}")
            series[(r.label, pn)] = S
            for a, k in zip(ax1, keys):
                a.plot(r.t, S[k], color=r.color, ls=LS[pn], lw=1.8, label=f"{r.label} {pn}")

    titles = ["rms(Jn_par - Jn_sheath)/rms(Jn_sheath), interior  (0.01 = row binds)", "median e Phi / T_e, interior",
              "interior current [kA] (+ ions)", "end-element current [kA]", "median Jn_sat / Bohm", "median |Mach_par|"]
    for a, k, ttl in zip(ax1, keys, titles):
        a.set_title(ttl, loc="left", fontsize=10)
        a.axhline(0, color="0.6", lw=0.6)
    ax1[1].axhline(args.sheath_lambda, color="0.6", lw=0.8, ls="--")
    ax1[4].axhline(1.0, color="0.6", lw=0.8, ls="--")
    ax1[5].axhline(1.0, color="0.6", lw=0.8, ls="--")
    for a in ax1[4:]:
        a.set_xlabel("time [ms]")
    ax1[0].legend(fontsize=8, frameon=False, ncol=2)

    # profiles along the plates
    picks = args.steps or [int(runs[0].steps[0]), int(runs[0].steps[-1])]
    fig2, ax2 = plt.subplots(4, 2 * len(picks), figsize=(6 * len(picks), 11), constrained_layout=True, num="plate profiles")
    fig3, ax3 = plt.subplots(1, len(picks), figsize=(5.5 * len(picks), 4.8), constrained_layout=True, squeeze=False, num="Jn_par vs Jn_sheath")
    for c, step in enumerate(picks):
        for r in runs:
            i = int(np.argmin(np.abs(r.steps - step)))
            for q, pn in enumerate(("inner", "outer")):
                _, P = plate_state(r, pn, i, mass, args.sheath_lambda)
                col = 2 * c + q
                s = P["s"] * 100
                a = ax2[:, col]
                a[0].plot(s, P["j"] * 1e-3, color=r.color, lw=1.8, label=f"{r.label} Jn_par")
                a[0].plot(s, P["jsh"] * 1e-3, color=r.color, lw=1.0, ls=":", label=f"{r.label} Jn_sheath (target)")
                a[1].plot(s, P["Phi"], color=r.color, lw=1.8, label=f"{r.label} Phi")
                a[1].plot(s, P["Phi_fl"], color=r.color, lw=1.0, ls=":", label=f"{r.label} Lambda T_e")
                a[2].plot(s, P["Te"], color=r.color, lw=1.8, label=r.label)
                a[3].plot(s, P["mach"], color=r.color, lw=1.8, label=r.label)
                ok = P["keep"]
                ax3[0, c].scatter(P["jsh"][ok] * 1e-3, P["j"][ok] * 1e-3, s=10, color=r.color, marker="o" if pn == "outer" else "^",
                                  alpha=0.6, label=f"{r.label} {pn} (step {r.steps[i]})")
                if r is runs[0]:
                    for x in a:
                        for lo, hi in ((s[0], s[np.argmax(ok)]), (s[len(ok) - 1 - np.argmax(ok[::-1])], s[-1])):
                            x.axvspan(lo, hi, color="0.9", zorder=0)
                a[0].set_title(f"{pn}, step {r.steps[i]}: j_n [kA/m^2]", loc="left", fontsize=9)
                a[1].set_title("Phi [V], dotted Lambda T_e", loc="left", fontsize=9)
                a[2].set_title("T_e [eV]", loc="left", fontsize=9)
                a[3].set_title("|Mach_par|", loc="left", fontsize=9)
                a[3].set_xlabel("s along the plate [cm] (shaded: end elements)")
        lim = ax3[0, c].get_xlim() + ax3[0, c].get_ylim()
        lo, hi = min(lim), max(lim)
        ax3[0, c].plot([lo, hi], [lo, hi], color="k", lw=0.8)
        ax3[0, c].set_xlabel("Jn_sheath (row target) [kA/m^2]")
        ax3[0, c].set_ylabel("Jn_par [kA/m^2]")
        ax3[0, c].legend(fontsize=8, frameon=False)
    for x in ax2.ravel():
        x.axhline(0, color="0.6", lw=0.6)
    for x in ax2[0]:                                   # y range from the actual currents; targets may run off-scale
        ys = np.concatenate([l.get_ydata() for l in x.get_lines() if l.get_linestyle() == "-"] or [np.zeros(1)])
        span = max(np.nanmax(np.abs(ys)), 1e-3)
        x.set_ylim(-1.3 * span, 1.3 * span)
    ax2[0, 0].legend(fontsize=7, frameon=False)
    ax2[1, 0].legend(fontsize=7, frameon=False)
    plt.show()


if __name__ == "__main__":
    main()
