#!/usr/bin/env python3
"""Which term of JOREK's vorticity equation balances the parallel current along a flux tube: the discrete answer,
from JOREK's own per-term output (RHS_terms_vtk) integrated over the tube.
  jorek2_postproc < rhs_u_aug.pp      ->  postproc_rhs_u/RHS.01300.vtk   (u row, equation index 2)
  python3 util/rhs_tube_budget.py postproc_rhs_u/RHS.01300.vtk --rect postproc/exprs_..._s001300.h5 \\
      --exprs "R Z Psi Psi_N BR BZ Btor currdens JR JZ Phi T_e T_i ne omega" --F0 4.23130791805 \\
      --bands 1.02:1.04 1.05:1.10 [--gap 0.01] [-o rhs_tube_budget.pdf]
RHS_terms_vtk writes the finite-element projection of every term group of the u row (one factor(var_u,i) each,
names in model600/mod_model_settings.f90) on an element subdivision: exactly what JOREK assembles, all terms with the
same weights (time step, R factors of the weak form). The row is charge conservation, F0 R div J = 0, with
  u_Eq__JxB            divergence of the parallel current  [psi, j]        u_Eq__grad_p      diamagnetic (-[R^2, p])
  u_Eq__visco_term     perpendicular viscosity                              u_Eq__visco_num_term  hyperviscosity
  u_Eq__rho_v.grad_v   Reynolds stress of the ExB flow                      u_Eq__diamag_term  gyroviscous (tauIC)
  u_Eq__zeta_timevol_term  polarisation (time change of the vorticity)      u_Eq__tg_num_term  Taylor-Galerkin
  u_Eq__ext_dens_source    momentum change by particle sources              u_Eq__neoclassical_term, rep/epf_pressure
For each Psi_N band the script takes the tube of util/check_current_conservation.py (same contours, same end windows),
integrates every term over the tube area between the windows and prints: the raw integral, its ratio to the JxB term,
and the current it corresponds to, scaled so that the JxB integral equals the parallel current gained along the
tube from the rectangle (dI_par). The sum over all terms is the residual of the implicit step (0 in a steady state).
Figure: cumulative integral of each significant term along the tube, against the poloidal distance from the inner end,
next to the parallel current itself: where along the tube each term acts."""
import argparse, os, re, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_current_conservation as cc                                             # noqa: E402
from verify_flux_tubes import read_rectangle                                        # noqa: E402

MU0 = 4e-7 * np.pi


def read_vtk_cells(path):
    """Points (n, 2), quad cells (m, 4) and {name: point values} of a legacy binary (big-endian) VTK of JOREK."""
    b = open(path, "rb").read()

    def header(key, pos):
        i = b.index(key, pos)
        j = b.index(b"\n", i)
        return b[i:j].decode().split(), j + 1

    tok, pos = header(b"POINTS", 0)
    n = int(tok[1])
    xyz = np.frombuffer(b, dtype=">f4", count=3 * n, offset=pos).reshape(n, 3).astype(float)
    pos += 12 * n
    tok, pos = header(b"CELLS", pos)
    ncell, size = int(tok[1]), int(tok[2])
    raw = np.frombuffer(b, dtype=">i4", count=size, offset=pos)
    pos += 4 * size
    cells, k = [], 0
    while k < size:
        c = raw[k]
        cells.append(raw[k + 1:k + 1 + c])
        k += 1 + c
    tok, pos = header(b"CELL_TYPES", pos)
    pos += 4 * int(tok[1])
    _, pos = header(b"POINT_DATA", pos)
    fields = {}
    while True:
        i = b.find(b"SCALARS ", pos)
        if i < 0:
            break
        j = b.index(b"\n", i)
        name = re.sub(r"\s+float$", "", b[i + 8:j].decode()).strip()
        k = b.index(b"LOOKUP_TABLE default\n", j) + len(b"LOOKUP_TABLE default\n")
        fields[name] = np.frombuffer(b, dtype=">f4", count=n, offset=k).astype(float)
        pos = k + 4 * n
    quads = np.array([c for c in cells if c.size == 4])
    return xyz[:, :2], quads, fields


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vtk", help="RHS.<step>.vtk of RHS_terms_vtk for the u row")
    ap.add_argument("--rect", required=True, help="rectangle expressions file of the same step (tube geometry, I_par)")
    ap.add_argument("--exprs", default="R Z Psi Psi_N BR BZ Btor currdens JR JZ Phi T_e T_i ne omega")
    ap.add_argument("--bands", nargs="+", default=["1.02:1.04", "1.05:1.10"])
    ap.add_argument("--F0", type=float, default=4.23130791805)
    ap.add_argument("--points", type=int, default=400)
    ap.add_argument("--gap", type=float, default=0.01)
    ap.add_argument("--frac", type=float, default=0.02, help="terms below this fraction of the largest are grouped as 'other'")
    ap.add_argument("-o", default="rhs_tube_budget.pdf")
    args = ap.parse_args()
    plt.rcdefaults()

    pts, quads, fields = read_vtk_cells(args.vtk)
    fields = {k: v for k, v in fields.items() if np.any(v != 0)}
    if not fields:
        sys.exit(f"{args.vtk}: no non-zero terms")
    # cell areas and cell-mean values
    P = pts[quads]                                                                # (m, 4, 2)
    area = 0.5 * np.abs(np.sum(P[:, :, 0] * np.roll(P[:, :, 1], -1, axis=1) - np.roll(P[:, :, 0], -1, axis=1) * P[:, :, 1], axis=1))
    cen = P.mean(axis=1)
    cval = {k: v[quads].mean(axis=1) for k, v in fields.items()}
    print(f"{args.vtk}: {pts.shape[0]} points, {quads.shape[0]} cells, terms: {', '.join(k.split('__')[-1] for k in fields)}")

    exprs = args.exprs.split()
    fld, R, Z, inside, t_now = read_rectangle(args.rect, exprs)
    fld["_beta"] = cc.unmirror_psin(fld, inside)
    bands = [tuple(float(x) for x in b.split(":")) for b in args.bands]
    fig, ax = plt.subplots(len(bands), 1, figsize=(8, 3.6 * len(bands)), squeeze=False, constrained_layout=True)
    for r, (a, b) in enumerate(bands):
        tb = cc.tube_budget(R, Z, fld, inside, a, b, args.F0, args.points, args.gap, 1.0)
        if tb is None:
            print(f"{a}:{b}: no open contour")
            continue
        poly = Path(np.vstack([tb["paths"]["a"][0], tb["paths"]["b"][0][::-1]]))
        in_tube = poly.contains_points(cen)
        c_pts = tb["paths"]["c"][0]
        i_near = np.argmin(np.hypot(cen[in_tube, None, 0] - c_pts[None, :, 0], cen[in_tube, None, 1] - c_pts[None, :, 1]), axis=1)
        s_c = tb["s"][i_near]
        s0, s1 = tb["s"][np.flatnonzero(tb["win_in"])[0]], tb["s"][np.flatnonzero(tb["win_out"])[-1]]
        keep = (s_c >= s0) & (s_c <= s1)
        idx = np.flatnonzero(in_tube)[keep]
        w = area[idx]
        integ = {k: float(np.sum(cval[k][idx] * w)) for k in fields}
        jxb = [k for k in fields if k.endswith("JxB")]
        if not jxb:
            sys.exit("no u_Eq__JxB term in the file: is this the u row (RHS_terms_vtk 2)?")
        I_jxb = integ[jxb[0]]
        dIpar = tb["I_par_out"] - tb["I_par_in"]                                     # A, parallel current gained
        scale = dIpar / I_jxb if I_jxb != 0 else np.nan                               # A per unit of term integral
        top = max(abs(v) for v in integ.values())
        print(f"\nband {a}:{b}, L = {tb['L']:.2f} m, {idx.size} cells in the tube between the windows; "
              f"dI_par (rectangle) = {dIpar * 1e-3:+.3f} kA")
        print(f"  {'term':28s}{'integral':>14s}{'/ JxB':>10s}{'entering [kA]':>15s}")
        tot = 0.0
        for k, v in sorted(integ.items(), key=lambda kv: -abs(kv[1])):
            tot += v
            ent = -v * scale * 1e-3 if not k.endswith("JxB") else np.nan             # entering by carrier X = -int div J_X
            print(f"  {k.split('__')[-1]:28s}{v:14.4e}{v / I_jxb:10.3f}{ent:15.3f}")
        print(f"  {'sum (residual of the step)':28s}{tot:14.4e}{tot / I_jxb:10.3f}")
        print("  entering = current the term deposits in the tube, with the JxB term set equal to dI_par; the time term"
              " (zeta_timevol) is the polarisation current; a large tg_num/visco_num share means numerical closure.")
        # cumulative along the tube
        A = ax[r, 0]
        order = np.argsort(s_c[keep])
        ell = tb["s"] * tb["L"]
        keep_names = [k for k in fields if abs(integ[k]) >= args.frac * top]
        for k in keep_names:
            cum = np.cumsum((cval[k][idx] * w)[order])
            A.plot(s_c[keep][order] * tb["L"], -cum * scale * 1e-3 if not k.endswith("JxB") else cum * scale * 1e-3,
                   lw=1.6 if k.endswith("JxB") else 1.2, label=k.split("__")[-1] + (" (= dI_par)" if k.endswith("JxB") else ""))
        A.plot(ell, (tb["I_par"] - tb["I_par_in"]) * 1e-3, "k--", lw=1.0, label=r"$I_\parallel(s) - I_\parallel(0)$ (rectangle)")
        A.axhline(0, color="0.6", lw=0.6)
        A.set_title(rf"$\Psi_N$ = {a:g} .. {b:g}: cumulative current deposited by each term along the tube", loc="left", fontsize=10)
        A.set_ylabel("current [kA]")
        A.legend(frameon=False, fontsize=8)
        if r == len(bands) - 1:
            A.set_xlabel("poloidal distance from the inner end [m]")
    fig.savefig(args.o)
    print(f"wrote {args.o}")
    plt.show()


if __name__ == "__main__":
    main()
