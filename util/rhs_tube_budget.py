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


def check_line(pts, fields, fld, R, Z, inside, F0):
    """Pointwise calibration of the VTK values: along the outer midplane (Z of the magnetic axis, Psi_N 0.97..1.15)
    print c*R of the JxB and grad_p terms next to the strong forms [psi, j] and -[R^2, p] from the rectangle
    (normalised units: j = -mu0 R J_phi, p = mu0 e n (T_e+T_i)). Their ratio is the common weight of the VTK terms
    (time step etc.) and must be the same constant for both terms and along the line; where it is not, the projection
    does not represent the strong form there (grid-scale structure, leakage from a sharp layer)."""
    MU0 = 4e-7 * np.pi; E = 1.602176634e-19
    RR = np.meshgrid(R, Z, indexing="ij")[0]
    dR, dZ = R[1] - R[0], Z[1] - Z[0]
    psi = np.where(inside, fld["Psi"], np.nan)
    j = np.where(inside, -MU0 * RR * fld["currdens"], np.nan)
    Ti = fld["T_i"] if "T_i" in fld else fld["T_e"]
    pn = np.where(inside, MU0 * E * fld["ne"] * (fld["T_e"] + Ti), np.nan)
    br = lambda a, b: np.gradient(a, dR, axis=0) * np.gradient(b, dZ, axis=1) - np.gradient(a, dZ, axis=1) * np.gradient(b, dR, axis=0)
    jxb_strong = br(psi, j)
    gp_strong = -br(RR ** 2, pn)
    # magnetic axis: minimum of Psi_N inside
    pN = np.where(inside, fld["Psi_N"], np.inf); ia, ja = np.unravel_index(np.argmin(pN), pN.shape)
    Zax = Z[ja]; jz = ja
    sel = np.flatnonzero(inside[:, jz] & (fld["Psi_N"][:, jz] > 0.97) & (fld["Psi_N"][:, jz] < 1.15) & (R > R[ia]))
    kj = [k for k in fields if k.endswith("JxB")][0]; kp = [k for k in fields if k.endswith("grad_p")][0]
    print(f"\ncheck along the outer midplane, Z = {Zax:.3f} m: VTK c*R against the strong forms from the rectangle")
    print(f"  {'R [m]':>7s} {'Psi_N':>7s} {'JxB vtk*R':>12s} {'[psi,j]':>12s} {'ratio':>10s} {'grad_p vtk*R':>13s} {'-[R2,p]':>12s} {'ratio':>10s}")
    ratios = []
    for i in sel[:: max(1, sel.size // 14)]:
        Rl = R[i]
        d = np.hypot(pts[:, 0] - Rl, pts[:, 1] - Zax)
        near = np.argsort(d)[:8]                                                   # the subdivision points around (R, Z)
        vj = np.mean(fields[kj][near]) * Rl; vp = np.mean(fields[kp][near]) * Rl
        print(f"  {Rl:7.3f} {fld['Psi_N'][i, jz]:7.4f} {vj:12.3e} {jxb_strong[i, jz]:12.3e} {vj / jxb_strong[i, jz]:10.3e}"
              f" {vp:13.3e} {gp_strong[i, jz]:12.3e} {vp / gp_strong[i, jz]:10.3e}")
        if 1.05 <= fld["Psi_N"][i, jz] <= 1.12 and np.isfinite(gp_strong[i, jz]) and gp_strong[i, jz] != 0:
            ratios.append(vp / gp_strong[i, jz])
    W = float(np.median(ratios)) if ratios else np.nan
    print(f"  weight of the VTK terms (median grad_p ratio, Psi_N 1.05..1.12): W = {W:.3e}; absolute current of a term"
          f" = -(2 pi/(mu0 F0)) * integral / W")
    return W


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
    ap.add_argument("--check-line", action="store_true",
                    help="compare the VTK's JxB and grad_p densities with [psi, j] and -[R^2, p] from the rectangle along the outer midplane")
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
    W = check_line(pts, fields, fld, R, Z, inside, args.F0) if args.check_line else np.nan
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
        w = area[idx] * cen[idx, 0]                                                # the mass matrix carries R dA: tau = c R
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
        absW = (2 * np.pi / (MU0 * args.F0)) / W if np.isfinite(W) else np.nan          # A per unit integral, from the calibration
        print(f"  {'term':28s}{'integral':>14s}{'/ JxB':>10s}{'entering [kA]':>15s}{'[kA, W-calib.]':>16s}")
        tot = 0.0
        for k, v in sorted(integ.items(), key=lambda kv: -abs(kv[1])):
            tot += v
            ent = -v * scale * 1e-3 if not k.endswith("JxB") else np.nan             # entering by carrier X = -int div J_X
            entW = (-v if not k.endswith("JxB") else v) * absW * 1e-3
            print(f"  {k.split('__')[-1]:28s}{v:14.4e}{v / I_jxb:10.3f}{ent:15.3f}{entW:16.3f}")
        print(f"  {'sum (residual of the step)':28s}{tot:14.4e}{tot / I_jxb:10.3f}")
        if np.isfinite(W):
            print(f"  W-calibrated: JxB = parallel current gained = {I_jxb * absW * 1e-3:+.3f} kA against {dIpar * 1e-3:+.3f} kA from the rectangle")
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
