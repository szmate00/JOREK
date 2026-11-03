#!/usr/bin/env python3
"""Parallel current along single flux tubes: is the current that leaves one target the current that arrives at the
other, or does it leak across the field on the way?

  python3 util/paper_jpar_tubes.py --levels 0.99 1.01 1.03 --exprs "R Z Psi_N Phi T_e T_i ne pres_e currdens eta_T dPsi_dt BR BZ Btor vpar" \\
      -o figs_bias -- -100V=run_m100/postproc/exprs_..._s000900.h5 0V=run_0/postproc/exprs_..._s000900.h5 +60V=...

  jpar_tubes.pdf  top row: j_par/B along each flux tube against the poloidal distance from its inner end, one column per
                  Psi_N level, one line per run; bottom row: the potential Phi along the same tubes (where the voltage
                  between the two targets drops).

In steady state div j = 0 gives B.grad(j_par/B) = -div(j_perp): without cross-field current j_par/B is constant along a
flux tube, so a change of j_par/B between the two ends is current that crosses the field on the way (diamagnetic,
viscous, polarisation). In reduced MHD the parallel current is carried by the toroidal one: j_par = j_phi B_phi/B, with
j_phi = `currdens`. The sign is chosen so that + means current flowing along the tube from its inner to its outer end
(inner/outer = the end at smaller/larger R). Levels above 1 are SOL tubes, below 1 private-region tubes (the rectangle
output does not mirror Psi_N). Rectangle files as for paper_axisym_flow.py; put "--" before labels starting with "-".
Printed per tube and run: j_par/B averaged over --end cm next to each end (but --gap cm away from the wall, outside the
last element), their ratio, and the largest |change| along the tube with its position.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_flux_tubes import read_rectangle, bilinear                           # noqa: E402

FS = 13
MUTED = "0.5"
COLORS = ["#4FAEE5", "0.35", "#E08230", "#C8D34D", "#00B1EA", "#EF7C00"]


def contour_pieces(R, Z, field, level):
    f = plt.figure()
    try:
        segs = f.gca().contour(R, Z, field.T, levels=[level], corner_mask=False).allsegs[0]
    finally:
        plt.close(f)
    return [np.asarray(s) for s in segs if len(s) > 10]


def tube_path(R, Z, fld, inside, level):
    """The flux tube Psi_N = level from its inner (small R) to its outer end: for SOL levels the longest open contour
    piece, for private-region levels the longest piece below the X-point."""
    psi = np.ma.masked_where(~inside, fld["Psi_N"])
    pieces = contour_pieces(R, Z, psi, level)
    if level < 1.0:                                                      # private region: below the X-point
        closed = max(contour_pieces(R, Z, psi, 0.995), key=len)
        z_x = closed[:, 1].min()
        pieces = [p for p in pieces if p[:, 1].mean() < z_x]
    pieces = [p for p in pieces if not np.allclose(p[0], p[-1])]
    if not pieces:
        return None
    path = max(pieces, key=lambda p: np.sum(np.hypot(*np.diff(p, axis=0).T)))
    return path if path[0, 0] <= path[-1, 0] else path[::-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", metavar="LABEL=FILE")
    ap.add_argument("--levels", type=float, nargs="+", default=[0.99, 1.01, 1.03])
    ap.add_argument("--exprs", default="R Z Psi_N Phi T_e T_i ne pres_e currdens eta_T dPsi_dt BR BZ Btor vpar")
    ap.add_argument("--end", type=float, default=3.0, help="length [cm] of the averaging window next to each end")
    ap.add_argument("--gap", type=float, default=1.0, help="distance [cm] from the wall left out (last element)")
    ap.add_argument("-o", "--outdir", default=".")
    args = ap.parse_args()
    plt.rcdefaults()
    os.makedirs(args.outdir, exist_ok=True)
    exprs = args.exprs.split()
    for need in ("currdens", "BR", "BZ", "Btor", "Phi", "Psi_N"):
        if need not in exprs:
            sys.exit(f"--exprs must contain {need}")

    nl = len(args.levels)
    fig, ax = plt.subplots(2, nl, figsize=(4.6 * nl, 7.5), sharex="col", constrained_layout=True, squeeze=False)
    print(f"{'run':>8s}{'Psi_N':>7s}{'L [m]':>7s}{'jB_in':>10s}{'jB_out':>10s}{'out/in':>8s}{'max |change|':>14s}{'at [m]':>8s}"
          f"   [j_par/B in kA m^-2 T^-1, + from the inner to the outer end]")
    for k, item in enumerate(args.runs):
        label, _, path = item.partition("=")
        col = COLORS[k % len(COLORS)]
        fld, R, Z, inside, t_now = read_rectangle(path, exprs)
        for j, lev in enumerate(args.levels):
            p = tube_path(R, Z, fld, inside, lev)
            if p is None:
                print(f"{label:>8s}{lev:7.3f}  no open contour")
                continue
            g = {n: bilinear(R, Z, fld[n], p) for n in ("currdens", "BR", "BZ", "Btor", "Phi")}
            B = np.sqrt(g["BR"] ** 2 + g["BZ"] ** 2 + g["Btor"] ** 2)
            jpar = g["currdens"] * g["Btor"] / B                          # reduced MHD: j_par = j_phi B_phi / B
            t = np.gradient(p, axis=0)
            t /= np.maximum(np.hypot(t[:, 0], t[:, 1]), 1e-30)[:, None]
            sign = np.sign(np.nanmean(g["BR"] * t[:, 0] + g["BZ"] * t[:, 1]) * np.nanmean(g["Btor"]))
            q = sign * jpar / B * 1e-3                                     # + along the tube from the inner to the outer end
            ell = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
            L = ell[-1]
            win_in = (ell >= args.gap / 100) & (ell <= (args.gap + args.end) / 100)
            win_out = (ell <= L - args.gap / 100) & (ell >= L - (args.gap + args.end) / 100)
            q_in, q_out = float(np.nanmean(q[win_in])), float(np.nanmean(q[win_out]))
            inner = (ell >= args.gap / 100) & (ell <= L - args.gap / 100)
            dev = np.abs(q - q_in)
            dev[~inner] = np.nan
            im = int(np.nanargmax(dev))
            print(f"{label:>8s}{lev:7.3f}{L:7.2f}{q_in:10.3f}{q_out:10.3f}{q_out / q_in if q_in else np.nan:8.2f}"
                  f"{dev[im]:14.3f}{ell[im]:8.2f}")
            ax[0, j].plot(ell, np.where(inner, q, np.nan), color=col, lw=1.8, label=label)
            ax[1, j].plot(ell, np.where(inner, g["Phi"], np.nan), color=col, lw=1.8, label=label)
            if k == 0:
                ax[0, j].set_title(rf"$\Psi_N$ = {lev:g} ({'private region' if lev < 1 else 'SOL'})", fontsize=12, loc="left")
    for j in range(nl):
        ax[0, j].axhline(0.0, color=MUTED, lw=0.6)
        ax[1, j].set_xlabel("poloidal distance from the inner end [m]", fontsize=FS - 1)
    ax[0, 0].set_ylabel(r"$j_\parallel / B$ [kA m$^{-2}$ T$^{-1}$]" + "\n(+ from the inner to the outer end)", fontsize=FS - 1)
    ax[1, 0].set_ylabel(r"$\Phi$ [V]", fontsize=FS)
    ax[0, -1].legend(frameon=False, fontsize=10)
    out = os.path.join(args.outdir, "jpar_tubes.pdf")
    fig.savefig(out)
    print(f"wrote {out}")
    import matplotlib
    if matplotlib.get_backend().lower() not in ("agg", "pdf", "ps", "svg", "cairo", "template"):
        plt.show()


if __name__ == "__main__":
    main()
