#!/usr/bin/env python3
"""
Which way does the ExB circulation run through the private flux region?

Reads a `pol_line` export (util/postproc_pfr_line.in) and projects the ExB
velocity onto the traverse direction. For a line drawn from the INNER leg to
the OUTER leg, the projection v_l is positive when plasma moves inner -> outer.

    v_l > 0  ->  transport INNER -> OUTER   (opposes HFSHD)
    v_l < 0  ->  transport OUTER -> INNER   (the HFSHD direction)

Reported alongside the particle flux n*v_l, because the density varies strongly
across the PFR and the flux is what matters, and alongside a coherence measure
|mean(v)|/mean(|v|): where that is small the direction rotates along the line
and the mean is not a direction.

Usage:
  ./analyse_pfr_line.py pol_line_s03800..03900.dat
  ./analyse_pfr_line.py base.dat puff.dat --labels baseline 10x-puff
"""

import argparse
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from analyse_targets import read_boundary_file, pick_block, _smooth, _trapz  # noqa: E402


def analyse(label, path, step, smooth, psin_min):
    step, t_now, d = pick_block(read_boundary_file(path), step)
    R, Z = d["R"], d["Z"]
    phi = -d["Phi"]                       # SIGN FIX, as in analyse_targets
    vr, vz = _smooth(d["V_ExB_R"], smooth), _smooth(d["V_ExB_Z"], smooth)
    ne, psin = d["ne"], d["Psi_N"]

    # Distance along the traverse, its unit tangent, and its poloidal NORMAL.
    l = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(R), np.diff(Z)))])
    t = np.array([R[-1] - R[0], Z[-1] - Z[0]])
    t = t / np.linalg.norm(t)
    # Rotate the tangent by -90 deg in the poloidal plane. For a vertical cut
    # (dZ only) this gives n_hat = +R_hat, so positive flux = outward in R.
    nrm = np.array([t[1], -t[0]])

    # FLUX is through the surface the curve sweeps toroidally, so it takes the
    # NORMAL component. v.t_hat is flow ALONG the curve and is not a flux
    # through anything - it is kept only as a secondary diagnostic.
    vn = vr * nrm[0] + vz * nrm[1]
    vl = vr * t[0] + vz * t[1]
    vmag = np.hypot(vr, vz)

    # Keep the part of the line genuinely inside the PFR.
    m = psin > psin_min
    if m.sum() < 5:
        sys.exit(f"{path}: only {m.sum()} points with psi_N > {psin_min}; "
                 "the line is probably not inside the PFR - check its endpoints")

    vn_m, vl_m = np.nanmean(vn[m]), np.nanmean(vl[m])
    # Gamma = 2*pi * int n (v.n_hat) R dl  [particles/s through the revolved surface]
    integ = (ne * vn * R)[m]
    ll = l[m]
    o = np.argsort(ll)
    gamma = 2.0 * np.pi * float(_trapz(integ[o], ll[o]))
    coh = abs(np.nanmean(vn[m])) / max(np.nanmean(vmag[m]), 1e-30)

    print(f"\n{'=' * 72}\n{label}   step {step}\n{'=' * 72}")
    print(f"  line      ({R[0]:.3f}, {Z[0]:.3f})  ->  ({R[-1]:.3f}, {Z[-1]:.3f})"
          f"   length {l[-1]:.3f} m")
    print(f"  inside PFR (psi_N > {psin_min}): {m.sum()} of {m.size} points, "
          f"psi_N {psin[m].min():.3f}-{psin[m].max():.3f}")
    print(f"\n  Phi            {np.nanmin(phi[m]):9.3f} .. {np.nanmax(phi[m]):9.3f} V"
          f"   (drop across PFR {phi[m][0] - phi[m][-1]:+.3f} V)")
    print(f"  n_e            {np.nanmean(ne[m]) / 1e20:9.4f} x1e20 m^-3 (mean)")
    print(f"  T_e            {np.nanmean(d['T_e'][m]):9.4f} eV (mean)")
    print(f"  v.n_hat        {vn_m:9.2f} m/s   (normal to the cut)")
    print(f"  v.t_hat        {vl_m:9.2f} m/s   (along the cut - NOT a flux)")
    print(f"  |v_ExB|        {np.nanmean(vmag[m]):9.2f} m/s (mean magnitude)")
    print(f"  n_hat          ({nrm[0]:+.3f}, {nrm[1]:+.3f})   "
          f"+Gamma = flow along +n_hat")
    print(f"  GAMMA_ExB      {gamma:9.3e} particles/s "
          f"  = 2*pi * int n (v.n) R dl")
    print(f"  coherence      {coh:9.3f}")

    if coh < 0.3:
        print("\n  ** DIRECTION NOT RESOLVED - the flow rotates along the line **")
    elif gamma > 0:
        print("\n  ** net ExB transport along +n_hat"
              "  (for a vertical cut: INNER -> OUTER, opposing HFSHD) **")
    else:
        print("\n  ** net ExB transport along -n_hat"
              "  (for a vertical cut: OUTER -> INNER, the HFSHD direction) **")
    return dict(label=label, l=l, phi=phi, vl=vn, ne=ne, psin=psin, m=m,
                flux=ne * vn * R, vl_m=vn_m, flux_m=gamma, coh=coh)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--labels", nargs="*")
    ap.add_argument("--step", type=int)
    ap.add_argument("--smooth", type=int, default=9)
    ap.add_argument("--psin-min", type=float, default=1.0,
                    help="keep points with psi_N above this (default 1.0)")
    ap.add_argument("-o", "--out", default="pfr_line.png")
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()

    labels = a.labels or [f.split("/")[-1] for f in a.files]
    res = [analyse(lb, f, a.step, a.smooth, a.psin_min)
           for f, lb in zip(a.files, labels)]

    if len(res) == 2:
        print(f"\n{'=' * 72}\n  {res[0]['label']}  ->  {res[1]['label']}")
        print(f"    v.n_hat   : {res[0]['vl_m']:+9.2f}  ->  {res[1]['vl_m']:+9.2f} m/s")
        print(f"    GAMMA_ExB : {res[0]['flux_m']:+9.3e}  ->  "
              f"{res[1]['flux_m']:+9.3e}  particles/s")
        print("=" * 72)

    if a.no_plot:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib not available - text only.")
        return
    fig, ax = plt.subplots(1, 4, figsize=(18, 4))
    for r, c in zip(res, ["#1f77b4", "#d62728", "#2ca02c"]):
        m = r["m"]
        for i, (k, lab, sc) in enumerate([("phi", r"$\Phi$ [V]", 1),
                                          ("ne", r"$n_e$ [$10^{20}$m$^{-3}$]", 1e20),
                                          ("vl", r"$v\cdot\hat n$ [m/s]", 1),
                                          ("flux", r"$n\,(v\cdot\hat n)\,R$ [m$^{-1}$s$^{-1}$]", 1)]):
            ax[i].plot(r["l"][m], r[k][m] / sc, color=c, lw=1.6, label=r["label"])
            ax[i].set_ylabel(lab)
            ax[i].set_xlabel("distance along traverse [m]")
            ax[i].grid(alpha=0.3)
            if i >= 2:
                ax[i].axhline(0, color="k", lw=0.8)
    ax[0].legend(fontsize=8)
    fig.suptitle("PFR traverse, inner leg -> outer leg")
    fig.tight_layout()
    fig.savefig(a.out, dpi=140)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
