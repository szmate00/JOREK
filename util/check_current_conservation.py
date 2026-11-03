#!/usr/bin/env python3
"""Current conservation along flux tubes (inxflow with the sheath current BC): does the current that the sheath row
puts into one target come out of the other, and where does it cross the field on the way?

  python3 util/check_current_conservation.py --rect V0/postproc/exprs_..._s001000.h5 \\
      --boundary V0/postproc/boundary_quantities_....dat --bands 1.01:1.03 1.04:1.08 0.97:0.99 -o figs

For every band Psi_N in [a, b] (SOL for a > 1, private region for b < 1, physical normalisation: the postproc mirrors
Psi_N below the X-point, Psi_N -> 2 - Psi_N, which the script undoes with the unmirrored `Psi`) the script builds the
flux tube from the rectangle output, Psi_N contours from the inner (small R) to the outer end, and evaluates the
axisymmetric current budget of the tube, all currents in kA for the full torus (2 pi x per radian):

  I_par(s)   parallel current through the tube cross-section at the poloidal distance s from the inner end (straight
             sections perpendicular to the central contour, between the surfaces a and b), + from the inner to the
             outer end: I = 2 pi int J_par b_pol.t R dn with J_par/B = J_phi R / F0 (reduced MHD, the current the
             sheath row sets through Jn_par = -j B_n/(mu0 F0)).
  I_pol(s)   the same with the full poloidal current of the postproc (`JR`, `JZ` = parallel + diamagnetic part).
  I_side     cross-field current leaving the tube through its two sides (J.grad(psi)/|grad psi|, cumulative in s).
  residual   I_pol(s) - I_pol(0) + I_side(s): the current the postproc does not carry (polarisation/viscous/inertial
             cross-field current of the vorticity equation) that left the tube between its inner end and s. In a
             converged steady state with negligible polarisation current it is flat at zero.
  div check  the same residual from the area integral of div(J_pol) over the tube (independent bookkeeping).
  ends       the plate currents of the boundary_quantities output in the band, 2 pi sum Jn R dl over the target points
             with a <= Psi_N <= b, for Jn_par (the current carried) and Jn_sheath (the characteristic on the local state).
             Current conservation at the ends: I_inner,plate = -I_par(0), I_outer,plate = +I_par(1) (Jn_* are positive
             INTO the wall). The last --gap cm next to each plate are left out of the tube averages (last element).

Current composition (Rozhansky, Kaveeva, Senichenkov, Sorokina, Vekshina, Coster, McCarthy, Khromov, NME 25 (2020)
100840; Sorokina et al. 2020; Zhao et al. arXiv:2601.14749): the plate currents of a tube are split into
  through     (I_out,pl - I_in,pl)/2: the part flowing from one plate to the other along the tube = thermoelectric
              current (Staebler-Hinton) plus the Pfirsch-Schlueter current that reaches the plates; opposite signs at
              the two plates (+: inner -> outer in the plasma, i.e. the inner plate collects ions, the outer electrons)
  closing     I_in,pl + I_out,pl: the "plate closing current" (PCC) of Rozhansky et al., the same sign at both
              plates; it enters the tube through its sides and is split into
                dia, below X   diamagnetic (grad-B) current through the sides below the X-point height = PCC type 1
                               (the vertical grad-B current of the divertor legs closed through the plates; Zhao's
                               source d j_par/ds = (1/BR) dp/dZ), the postproc's JR/JZ carry it exactly
                dia, above X   the same above the X-point = PCC type 2 / Pfirsch-Schlueter closure in the SOL
                residual       current the postproc's J does not carry: polarisation, viscous and inertial currents
                               of the vorticity equation (and the discretisation error), below / above the X-point
  thermo est. the thermoelectric current the open-circuit voltage Lambda (T_e,out - T_e,in) would drive through
              the tube's parallel Spitzer resistance (sum eta_par dl_par / A_par; --eta-si overrides eta_par),
              Zhao's estimate j_th = 2.8 (T_e,ou - T_e,in) sigma_par / (e L_par), given in the through-current sign convention, with the measured through-current
              and the end-to-end potential difference Phi_out - Phi_in next to it.
Written as a table and current_composition.pdf (stacked bars per band). If JR JZ are not in the rectangle but ne T_e
T_i are, the poloidal current is rebuilt as J_par b_pol (currdens R/F0 B_pol) + B x grad(p_tot)/B^2 with
p_tot = e n_e (T_e + T_i), the same decomposition the postproc uses.

The rectangle needs --exprs with at least R Z Psi Psi_N BR BZ Btor currdens JR JZ (same order as in the postproc
command; JR JZ may be replaced by ne T_e T_i, see above); the boundary file needs R Z Psi_N bnd_dl Bnorm Jn_par (Jn_sheath, Jnorm, bnd_type_lo/hi optional). The
boundary step closest in time to the rectangle is used unless --step is given. Rectangle and boundary files as for
util/paper_jpar_tubes.py and util/compare_wall_runs.py. Also printed: the total current into each plate, into the whole
wall and per boundary type, for Jn_par (parallel, what the sheath row sets) and Jnorm (parallel + diamagnetic). The
model's wall current is Jn_par only: at the released plate nodes the vorticity row drops the boundary term of the
non-parallel currents, so Jnorm is a pointwise field value, not a discrete wall flux. Summing the vorticity rows of
all u value DOFs gives F0 x (whole-wall Jn_par total) = sum of the residuals at the PINNED u value DOFs (types with
Dirichlet/floating u: the type-2 flux-surface walls, corners, plate nodes below min_sheath_angle) minus the volume
integral of the strong-form tauIC gyroviscous terms. A non-zero Jn_par total is therefore the "ground current" taken
up by the pinned rows (review 2026-10-23), not a diagnostic or postproc error.
"""
import argparse
import math
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_flux_tubes import read_rectangle, bilinear                           # noqa: E402
from paper_jpar_tubes import contour_pieces                                      # noqa: E402
from compare_wall_runs import read_fast                                          # noqa: E402

FS = 11
COLORS = {"par": "#1f5fa8", "pol": "#E08230", "res": "0.25", "div": "#4FAEE5", "a": "#7a9f35", "b": "#b03a2e"}


def resample(path, n):
    """Path (m, 2) resampled to n points uniform in arclength; returns points, unit tangents, fraction s, length."""
    seg = np.hypot(*np.diff(path, axis=0).T)
    ell = np.concatenate([[0.0], np.cumsum(seg)])
    L = ell[-1]
    s = np.linspace(0.0, 1.0, n)
    pts = np.column_stack([np.interp(s * L, ell, path[:, 0]), np.interp(s * L, ell, path[:, 1])])
    t = np.gradient(pts, axis=0)
    t /= np.maximum(np.hypot(t[:, 0], t[:, 1]), 1e-30)[:, None]
    return pts, t, s, L


def unmirror_psin(fld, inside):
    """The postproc's Psi_N is mirrored below the X-point (private region: Psi_N -> 2 Psi_X - Psi_N, get_psi_n with Z
    in equil_info.f90). Psi is not: fit Psi = alpha + beta Psi_N on the closed surfaces (Psi_N < 0.9) and replace Psi_N
    by (Psi - alpha)/beta everywhere, so that private-region surfaces have their physical Psi_N < 1."""
    sel = inside & (fld["Psi_N"] < 0.9)
    beta, alpha = np.polyfit(fld["Psi_N"][sel], fld["Psi"][sel], 1)
    fld["Psi_N"] = np.where(inside, (fld["Psi"] - alpha) / beta, 0.0)
    return beta


def xpoint_z(R, Z, fld, inside):
    """Z of the lower X-point: the lowest point of the longest closed contour Psi_N = 0.995."""
    psi = np.ma.masked_where(~inside, fld["Psi_N"])
    closed = max(contour_pieces(R, Z, psi, 0.995), key=len)
    return closed[:, 1].min()


def tube_path(R, Z, fld, inside, level, z_x):
    """The open contour Psi_N = level from its inner (small R) to its outer end: the longest open piece (SOL), or for
    level < 1 the longest open piece below the X-point (private region; the closed core surface is excluded)."""
    psi = np.ma.masked_where(~inside, fld["Psi_N"])
    pieces = [p for p in contour_pieces(R, Z, psi, level) if not np.allclose(p[0], p[-1])]
    if level < 1.0:
        pieces = [p for p in pieces if p[:, 1].mean() < z_x]
    if not pieces:
        return None
    path = max(pieces, key=lambda p: np.sum(np.hypot(*np.diff(p, axis=0).T)))
    return path if path[0, 0] <= path[-1, 0] else path[::-1]


def plates(bq, names, bnorm_frac=0.05, z_max=None):
    """Indices of the inner and outer target points of the boundary output (|B.n| above bnorm_frac of its maximum,
    split at the largest gap in R)."""
    c = {n: i for i, n in enumerate(names)}
    R, Z = bq[0, :, c["R"]], bq[0, :, c["Z"]]
    bn = np.abs(bq[:, :, c["Bnorm"]]).mean(axis=0)
    sel = bn > bnorm_frac * bn.max()
    if z_max is not None:
        sel &= Z < z_max
    tgt = np.where(sel)[0]
    Rs = np.sort(R[tgt])
    k = int(np.argmax(np.diff(Rs)))
    R_split = 0.5 * (Rs[k] + Rs[k + 1])
    pl = {"inner": tgt[R[tgt] < R_split], "outer": tgt[R[tgt] > R_split]}
    psiN = bq[0, :, c["Psi_N"]]
    for pn, side in (("inner", +1), ("outer", -1)):                          # private side: larger R on the inner plate
        i = pl[pn]
        R_strike = R[i][np.argmin(np.abs(psiN[i] - 1.0))]
        pl[pn + "_pfr"] = i[side * (R[i] - R_strike) > 0]
        pl[pn + "_sol"] = i[side * (R[i] - R_strike) <= 0]
    return pl, c


def plate_current(bq, c, idx, name, lo=None, hi=None, pfr=False):
    """2 pi sum Jn R dl [A] over the points idx of the (phi-averaged) boundary block, optionally within lo <= Psi_N <= hi
    (physical Psi_N; the boundary output's Psi_N is mirrored in the private region like the rectangle's)."""
    if name not in c:
        return np.nan
    g = bq.mean(axis=0)
    sel = idx
    if lo is not None:
        psi = g[idx, c["Psi_N"]]
        if pfr:
            lo, hi = 2.0 - hi, 2.0 - lo
        sel = idx[(psi >= lo) & (psi <= hi)]
    if sel.size == 0:
        return np.nan
    return 2.0 * np.pi * float(np.sum(g[sel, c[name]] * g[sel, c["R"]] * g[sel, c["bnd_dl"]]))


def wall_totals(bq, c, pl, z_x=None):
    """Total current into the wall [kA] from the boundary output: the parallel current Jn_par (what the sheath row
    sets; the model's only discrete wall current) and, if present, the pointwise poloidal current Jnorm (parallel +
    diamagnetic, for information). Per plate, whole boundary, and per boundary type (bnd_type_lo, elements whose two
    wall nodes have the same type; corner elements are listed under both types they touch)."""
    names = [n for n in ("Jn_par", "Jnorm") if n in c]
    everything = np.arange(bq.shape[1])
    for name in names:
        tot = {pn: plate_current(bq, c, pl[pn], name) for pn in ("inner", "outer")}
        allp = plate_current(bq, c, everything, name)
        label = "parallel current (Jn_par)   " if name == "Jn_par" else "total pol. current (Jnorm)  "
        print(f"{label} into the wall [kA]: inner plate {tot['inner'] * 1e-3:9.3f}, outer plate {tot['outer'] * 1e-3:9.3f},"
              f" sum {(tot['inner'] + tot['outer']) * 1e-3:9.3f}, whole boundary {allp * 1e-3:9.3f}")
    if "Jnorm" in c and "Jn_par" in c:
        dia = plate_current(bq, c, everything, "Jnorm") - plate_current(bq, c, everything, "Jn_par")
        print(f"  diamagnetic current through the whole boundary (Jnorm - Jn_par) {dia * 1e-3:9.3f} kA")
    if "bnd_type_lo" in c and "bnd_type_hi" in c:
        g = bq.mean(axis=0)
        lo, hi = g[:, c["bnd_type_lo"]].round().astype(int), g[:, c["bnd_type_hi"]].round().astype(int)
        types = sorted(set(lo.tolist()) | set(hi.tolist()))
        row = "  per boundary type [kA]:"
        R, Z = g[:, c["R"]], g[:, c["Z"]]
        R_in = R[pl["inner"]].max() if pl["inner"].size else -np.inf          # the plates' inner edges in R
        R_out = R[pl["outer"]].min() if pl["outer"].size else np.inf
        for t in types:
            idx = everything[(lo == t) | (hi == t)]
            vals = ", ".join(f"{n} {plate_current(bq, c, idx, n) * 1e-3:8.3f}" for n in names)
            row += f"\n    type {t:2d} ({idx.size:4d} pts): {vals}"
            if z_x is not None and t not in (1, 3) and idx.size:
                floor = idx[(Z[idx] < z_x) & (R[idx] > R_in) & (R[idx] < R_out)]      # PFR floor: below the X-point, between the plates
                edge = np.setdiff1d(idx, floor)
                for name, sub in (("PFR floor", floor), ("SOL edge ", edge)):
                    if sub.size:
                        vals = ", ".join(f"{n} {plate_current(bq, c, sub, n) * 1e-3:8.3f}" for n in names)
                        row += f"\n      {name} ({sub.size:4d} pts): {vals}"
        print(row)
        if "Jnorm" in c:
            tot = plate_current(bq, c, everything, "Jnorm") * 1e-3
            print(f"  charge balance of the postproc current (parallel + diamagnetic) over the whole boundary: {tot:+.3f} kA"
                  " (0 if that current alone closes; the rest is the vorticity equation's polarisation/viscous flux)")


def contiguous(mask, j0):
    """Of a boolean array (n, m), the run of True values containing column j0 in every row."""
    left = np.cumprod(mask[:, j0::-1], axis=1)[:, ::-1].astype(bool)
    right = np.cumprod(mask[:, j0:], axis=1).astype(bool)
    out = np.zeros_like(mask)
    out[:, :j0 + 1] = left
    out[:, j0:] = right
    return out


def section_integral(F, psiN, ins, sec, a, b, h):
    """int F dn over each row's contiguous run sec, with the band limits a, b located by linear interpolation of
    psiN into the neighbouring cells (trapezoid inside, partial cells at the ends). Rows without a run give 0."""
    n_s, n_d = F.shape
    has = sec.any(axis=1)
    j1 = np.argmax(sec, axis=1)
    j2 = n_d - 1 - np.argmax(sec[:, ::-1], axis=1)
    rows = np.arange(n_s)
    Fm = np.where(sec, F, 0.0)
    inner = Fm.sum(axis=1) * h - 0.5 * h * (F[rows, j1] + F[rows, j2])
    inner = np.where(has, inner, 0.0)

    def end_piece(j, jn):
        okn = has & (jn >= 0) & (jn < n_d)
        jn_c = np.clip(jn, 0, n_d - 1)
        pn, pj = psiN[rows, jn_c], psiN[rows, j]
        okn &= ins[rows, jn_c] & np.isfinite(pn)
        lim = np.where(pn < pj, a, b)                                       # the band limit crossed towards jn
        x = np.where(okn & (np.abs(pj - pn) > 0), (pj - lim) / np.where(pj == pn, 1.0, pj - pn), 0.0)
        x = np.clip(x, 0.0, 1.0)
        Fc = F[rows, j] + x * (F[rows, jn_c] - F[rows, j])
        return np.where(okn, x * h * 0.5 * (F[rows, j] + Fc), 0.0)

    return inner + end_piece(j1, j1 - 1) + end_piece(j2, j2 + 1)


def tube_budget(R, Z, fld, inside, a, b, F0, n_s, gap, d_max=1.0, eta_override=None):
    z_x = xpoint_z(R, Z, fld, inside)
    """Currents along the flux tube a <= Psi_N <= b. Returns a dict or None (no open contour).

    Cross-sections are straight lines through the points of the central contour Psi_N = (a+b)/2, perpendicular to it,
    between the surfaces Psi_N = a and b; the current through them is 2 pi int J.t R dn with t the tangent of the
    central contour (+ from the inner to the outer end)."""
    paths = {}
    for key, lev in (("a", a), ("c", 0.5 * (a + b)), ("b", b)):
        p = tube_path(R, Z, fld, inside, lev, z_x)
        if p is None:
            return None
        paths[key] = resample(p, n_s)
    pts, t, s, L = paths["c"]
    nrm = np.column_stack([-t[:, 1], t[:, 0]])
    h = 0.25 * min(R[1] - R[0], Z[1] - Z[0])
    d = np.arange(-d_max, d_max + h / 2, h)
    j0 = int(np.argmin(np.abs(d)))
    q = pts[:, None, :] + d[None, :, None] * nrm[:, None, :]                  # (n_s, n_d, 2)
    qf = q.reshape(-1, 2)
    ok_dom = (qf[:, 0] >= R[0]) & (qf[:, 0] <= R[-1]) & (qf[:, 1] >= Z[0]) & (qf[:, 1] <= Z[-1])
    psiN = np.where(ok_dom, bilinear(R, Z, fld["Psi_N"], qf), np.nan).reshape(q.shape[:2])
    ins = (bilinear(R, Z, inside.astype(float), qf) > 0.999).reshape(q.shape[:2]) & ok_dom.reshape(q.shape[:2])
    sec = contiguous(ins & (psiN >= a) & (psiN <= b), j0)
    g = {n: np.nan_to_num(bilinear(R, Z, fld[n], qf)).reshape(q.shape[:2]) for n in ("currdens", "BR", "BZ", "JR", "JZ")}
    Rq = q[..., 0]
    tx, ty = t[:, 0][:, None], t[:, 1][:, None]
    fpar = g["currdens"] * Rq / F0 * (g["BR"] * tx + g["BZ"] * ty)           # J_par b_pol . t  (J_par/B = J_phi R/F0)
    fpol = g["JR"] * tx + g["JZ"] * ty
    I_par = 2 * np.pi * section_integral(fpar * Rq, psiN, ins, sec, a, b, h)   # [A]
    I_pol = 2 * np.pi * section_integral(fpol * Rq, psiN, ins, sec, a, b, h)
    width = section_integral(np.ones_like(Rq), psiN, ins, sec, a, b, h)
    # sides: current leaving through the surface Psi_N = a (towards smaller Psi_N) and Psi_N = b (towards larger)
    sgn_n = np.sign(fld["_beta"])                                             # grad(Psi_N) direction
    side = {}
    for key, out_sign in (("a", -1.0), ("b", +1.0)):
        pk, tk, sk, Lk = paths[key]
        gk = {n: bilinear(R, Z, fld[n], pk) for n in ("BR", "BZ", "JR", "JZ")}
        gpsi = np.column_stack([-pk[:, 0] * gk["BZ"], pk[:, 0] * gk["BR"]])    # (psi_R, psi_Z) = (-R B_Z, R B_R)
        nk = gpsi / np.maximum(np.hypot(gpsi[:, 0], gpsi[:, 1]), 1e-30)[:, None] * sgn_n
        jn = (gk["JR"] * nk[:, 0] + gk["JZ"] * nk[:, 1]) * out_sign            # A/m^2 leaving the tube
        dl = np.gradient(sk) * Lk
        cum_k = 2 * np.pi * np.cumsum(jn * pk[:, 0] * dl)
        # map onto the central path's s: the side point's nearest central point
        i_c = np.argmin(np.hypot(pk[:, None, 0] - pts[None, :, 0], pk[:, None, 1] - pts[None, :, 1]), axis=1)
        s_k = s[i_c]
        order = np.argsort(s_k, kind="stable")
        cum = np.interp(s, s_k[order], cum_k[order])
        side[key] = {"jn": jn, "ell": sk * Lk, "cum": cum, "L": Lk, "Z": pk[:, 1], "dI": 2 * np.pi * jn * pk[:, 0] * dl}
    win_in = (s * L >= gap) & (s * L <= gap + 0.03)
    win_out = (s * L <= L - gap) & (s * L >= L - gap - 0.03)
    i_ref = int(np.argmax(win_in))                                            # budget relative to the inner window's start
    I_side = side["a"]["cum"] + side["b"]["cum"]
    I_side = I_side - I_side[i_ref]
    residual = I_pol - I_pol[i_ref] + I_side
    # divergence check on the rectangle: area integral of div(J_pol) over the tube polygon up to each s (grid points
    # assigned to the s of their nearest central point, i.e. the same perpendicular sections)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    JR, JZ = np.where(inside, fld["JR"], np.nan), np.where(inside, fld["JZ"], np.nan)
    div = np.gradient(RR * JR, R, axis=0) / RR + np.gradient(JZ, Z, axis=1)
    poly = Path(np.vstack([paths["a"][0], paths["b"][0][::-1]]))
    in_tube = poly.contains_points(np.column_stack([RR.ravel(), ZZ.ravel()])).reshape(RR.shape) & inside
    dA = (R[1] - R[0]) * (Z[1] - Z[0])
    pts_in = np.column_stack([RR[in_tube], ZZ[in_tube]])
    i_near = np.argmin(np.hypot(pts_in[:, None, 0] - pts[None, :, 0], pts_in[:, None, 1] - pts[None, :, 1]), axis=1)
    contrib = 2 * np.pi * np.nan_to_num(div[in_tube]) * RR[in_tube] * dA
    div_cum = np.cumsum(np.bincount(i_near, weights=contrib, minlength=n_s))
    div_cum = div_cum - div_cum[i_ref]
    # --- composition: the side current split by its carrier (diamagnetic in JR/JZ vs the residual) and by the
    # --- X-point height (divertor legs below = PCC type 1, SOL above = type 2 / Pfirsch-Schlueter closure); only the
    # --- part between the two end windows counts (the last --gap m next to each plate are excluded, as in the budget)
    i_out = int(np.flatnonzero(win_out)[-1])
    comp = {}
    for key in ("a", "b"):
        sd = side[key]
        on = (sd["ell"] >= gap) & (sd["ell"] <= sd["L"] - gap)
        below = sd["Z"] < z_x
        comp["dia_below"] = comp.get("dia_below", 0.0) + float(np.sum(sd["dI"][on & below]))
        comp["dia_above"] = comp.get("dia_above", 0.0) + float(np.sum(sd["dI"][on & ~below]))
    d_res = np.diff(residual)                                                 # per step of s, placed at the central path
    z_c = 0.5 * (pts[1:, 1] + pts[:-1, 1])
    on_c = np.zeros(n_s - 1, bool); on_c[i_ref:i_out] = True
    comp["res_below"] = float(np.sum(d_res[on_c & (z_c < z_x)]))
    comp["res_above"] = float(np.sum(d_res[on_c & (z_c >= z_x)]))
    # --- parallel length and Spitzer resistance of the tube (for the thermoelectric estimate)
    gB = {n: np.nan_to_num(bilinear(R, Z, fld[n], pts)) for n in ("BR", "BZ", "Btor")}
    B_pol = np.hypot(gB["BR"], gB["BZ"]); B_tot = np.hypot(B_pol, gB["Btor"])
    dl_pol = np.gradient(s) * L
    dl_par = dl_pol * B_tot / np.maximum(B_pol, 1e-12)
    comp["L_par"] = float(np.sum(dl_par[i_ref:i_out]))
    comp["R_tube"] = np.nan
    if "T_e" in fld:
        Te_c = np.maximum(np.nan_to_num(bilinear(R, Z, fld["T_e"], pts)), 0.1)      # eV
        eta_par = 5.2e-5 * 12.0 / Te_c ** 1.5                                          # Spitzer parallel, Z = 1, ln Lambda 12
        if eta_override is not None:
            eta_par = np.full_like(Te_c, eta_override)
        A_par = 2 * np.pi * pts[:, 0] * np.maximum(width, 1e-6) * B_pol / np.maximum(B_tot, 1e-12)   # cross-section perp. to B
        comp["R_tube"] = float(np.sum((eta_par * dl_par / A_par)[i_ref:i_out]))    # Ohm
        comp["Te_tube_min"] = float(Te_c[i_ref:i_out].min())
    return {"s": s, "L": L, "I_par": I_par, "I_pol": I_pol, "I_side": I_side, "residual": residual, "div_cum": div_cum,
            "side": side, "win_in": win_in, "win_out": win_out, "paths": paths, "width": width,
            "in_tube": in_tube, "i_near": i_near, "comp": comp, "z_x": z_x,
            "I_par_in": float(I_par[win_in].mean()), "I_par_out": float(I_par[win_out].mean()),
            "I_pol_in": float(I_pol[win_in].mean()), "I_pol_out": float(I_pol[win_out].mean())}


def _fill(f, inside):
    """f with every point outside the domain replaced by the value of the nearest point inside (so that finite
    differences near the boundary see no jump)."""
    from scipy.ndimage import distance_transform_edt
    idx = distance_transform_edt(~inside, return_distances=False, return_indices=True)
    return f[tuple(idx)]


def _smooth(f, inside, sig):
    """Nearest-value extrapolation outside the domain, then Gaussian smoothing (sig in cells, 0 = none)."""
    g = _fill(np.asarray(f, float), inside)
    if sig <= 0:
        return g
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(g, sig)


def carrier_budget(R, Z, fld, inside, tb, F0, visco, visco_num, n0, t0, sig, fld_prev=None, dt_norm=None, margin=6):
    """Net current each carrier of the vorticity equation deposits in the tube, from the divergences of the model's
    currents integrated over the tube area (+ = current ENTERING the tube, i.e. collected by the plates), in A.
    Phi is Gaussian-smoothed by sig cells before differentiation; n and p are not.
    Normalised units of model600 (lengths m, u = Phi_SI*t0/F0, rho = n/n0, p = mu0 p_SI), F0 R div J = rhs term:
      diamagnetic   div J = 2 dp/dZ / F0                      (term [R^2, p])
      viscous       div J = visco * Lap(omega) / F0,  omega = Lap_pol u   (term R mu Lap omega)
      hyperviscous  div J = -visco_num * Lap(Lap(omega)) / F0  (sign as the code's stabilising term; indicative only)
      Reynolds      div J = ( 1/2 [R^2 |grad u|^2, R^2 rho] + [R^4 rho omega, u] ) / (F0 R)
      polarisation  div J = -div( rho R^2 grad d_t u ) / F0,  d_t u from a second rectangle (fld_prev, dt_norm)
    J_SI = J/mu0. The diamagnetic entry is a check against the side integral of the postproc's JR/JZ."""
    MU0 = 4e-7 * np.pi; E = 1.602176634e-19
    RR = np.meshgrid(R, Z, indexing="ij")[0]
    dR, dZ = R[1] - R[0], Z[1] - Z[0]
    def dr(f): return np.gradient(f, dR, axis=0)
    def dz(f): return np.gradient(f, dZ, axis=1)
    def lap(f): return dr(dr(f)) + dr(f) / RR + dz(dz(f))
    def br(a, b): return dr(a) * dz(b) - dz(a) * dr(b)
    def div(fr, fz): return dr(RR * fr) / RR + dz(fz)
    u = _smooth(fld["Phi"] * t0 / F0, inside, sig)                 # third and fourth derivatives: smoothed
    rho = _fill(fld["ne"] / n0, inside)                               # first derivatives only: never smoothed (a steep
    Ti = fld["T_i"] if "T_i" in fld else fld["T_e"]                   # pedestal smeared across the separatrix would
    p = _fill(MU0 * E * fld["ne"] * (fld["T_e"] + Ti), inside)       # leak core pressure into the SOL tubes)
    om = lap(u)
    terms = {}
    terms["diamagnetic"] = 2.0 * dz(p) / F0
    terms["viscous"] = visco * lap(om) / F0
    if visco_num:
        terms["hyperviscous"] = -visco_num * lap(lap(om)) / F0
    gu2 = dr(u) ** 2 + dz(u) ** 2
    terms["Reynolds"] = (0.5 * br(RR ** 2 * gu2, RR ** 2 * rho) + br(RR ** 4 * rho * om, u)) / (F0 * RR)
    if fld_prev is not None and dt_norm:
        u_prev = _smooth(fld_prev["Phi"] * t0 / F0, inside, sig)
        dtu = (u - u_prev) / dt_norm
        terms["polarisation d/dt"] = -div(rho * RR ** 2 * dr(dtu), rho * RR ** 2 * dz(dtu)) / F0
    dA = dR * dZ
    m = tb["in_tube"] & inside
    # exclude the two end windows like the side budget (last --gap m next to each plate): points nearest to s outside
    s_pt = tb["s"][tb["i_near"]]
    keep = np.zeros_like(m); keep[m] = (s_pt >= tb["s"][np.flatnonzero(tb["win_in"])[0]]) & (s_pt <= tb["s"][np.flatnonzero(tb["win_out"])[-1]])
    # the 3rd/4th-derivative stencils of the viscous, Reynolds and polarisation terms reach `margin` cells inward from
    # the domain boundary (the fields are only C0 there after the nearest-value fill): leave that strip out of them
    from scipy.ndimage import distance_transform_edt
    deep = distance_transform_edt(inside) > margin + 2.0 * sig
    out = {}
    for k, v in terms.items():
        sel = keep if k == "diamagnetic" else (keep & deep)
        out[k] = -2.0 * np.pi * float(np.nansum(np.where(sel, v, 0.0) * RR * dA)) / MU0     # entering = -int div J
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rect", required=True, help="rectangle expressions file (h5)")
    ap.add_argument("--boundary", help="boundary_quantities file of the same run")
    ap.add_argument("--exprs", default="R Z Psi Psi_N BR BZ Btor currdens JR JZ Phi T_e ne",
                    help="expression names of the rectangle file, in order")
    ap.add_argument("--bands", nargs="+", default=["1.01:1.03", "1.04:1.08"], metavar="A:B",
                    help="Psi_N bands of the flux tubes")
    ap.add_argument("--F0", type=float, default=3.0, help="F0 of the run [T m] (inxflow: 3)")
    ap.add_argument("--points", type=int, default=400, help="points (cross-sections) along the tube")
    ap.add_argument("--max-width", type=float, default=1.0, help="longest cross-section searched [m]")
    ap.add_argument("--gap", type=float, default=0.01, help="length [m] next to each plate left out of the end averages")
    ap.add_argument("--step", type=int, help="boundary step to use (default: closest in time to the rectangle)")
    ap.add_argument("--plate-zmax", type=float, help="plate points only below this Z [m] (AUG)")
    ap.add_argument("--zoom", type=float, default=0.0,
                    help="length [m]: also plot the last --zoom m at each end without the gap mask (current_conservation_ends.pdf);"
                         " use with --points 2000 or more so that the sections resolve the last element")
    ap.add_argument("--sheath-lambda", type=float, default=3.0, help="Lambda of the run (open-circuit voltage Lambda dT_e)")
    ap.add_argument("--visco", type=float, default=2.3e-8, help="namelist visco (normalised, constant) for the carrier budget")
    ap.add_argument("--visco-num", type=float, default=0.0, help="namelist visco_num (hyperviscosity); 0 = not evaluated")
    ap.add_argument("--central-density", type=float, default=0.76464396)
    ap.add_argument("--central-mass", type=float, default=2.01410177811)
    ap.add_argument("--smooth", type=float, default=1.5, help="Gaussian smoothing [cells] of Phi before its 3rd/4th derivatives (n, p never smoothed)")
    ap.add_argument("--rect-prev", help="rectangle file at an earlier time (same --exprs): enables the polarisation d/dt term")
    ap.add_argument("--eta-si", type=float, help="constant parallel resistivity [Ohm m] instead of Spitzer in the thermoelectric estimate")
    ap.add_argument("-o", "--outdir", default=".")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    exprs = args.exprs.split()
    for need in ("R", "Z", "Psi", "Psi_N", "BR", "BZ", "Btor", "currdens"):
        if need not in exprs:
            sys.exit(f"--exprs must contain {need}")
    if not ("JR" in exprs and "JZ" in exprs) and not all(n in exprs for n in ("ne", "T_e", "T_i")):
        sys.exit("--exprs must contain JR JZ, or ne T_e T_i to rebuild them")
    fld, R, Z, inside, t_now = read_rectangle(args.rect, exprs)
    if "JR" not in fld:                                                       # parallel + diamagnetic, as the postproc's JR/JZ
        RR = np.meshgrid(R, Z, indexing="ij")[0]
        p = 1.602176634e-19 * fld["ne"] * (fld["T_e"] + fld["T_i"])            # Pa (si-units output: ne m^-3, T eV)
        p = np.where(inside, p, np.nan)
        pR, pZ = np.gradient(p, R, axis=0), np.gradient(p, Z, axis=1)
        B2 = fld["BR"] ** 2 + fld["BZ"] ** 2 + fld["Btor"] ** 2
        jpar_over_B = fld["currdens"] * RR / args.F0                            # J_par/B = J_phi R/F0 (reduced MHD)
        fld["JR"] = np.nan_to_num(jpar_over_B * fld["BR"] - fld["Btor"] * pZ / np.maximum(B2, 1e-30))
        fld["JZ"] = np.nan_to_num(jpar_over_B * fld["BZ"] + fld["Btor"] * pR / np.maximum(B2, 1e-30))
        print("JR/JZ rebuilt from currdens (parallel) and B x grad(e n_e (T_e + T_i))/B^2 (diamagnetic); one grid cell"
              " next to the boundary carries no diamagnetic part")
    fld["_beta"] = unmirror_psin(fld, inside)
    print(f"rectangle {args.rect}: t = {t_now:.4g} (JOREK units), {R.size} x {Z.size} points, "
          f"psi_bnd - psi_axis = {fld['_beta']:.4f} Wb/rad")

    AMU = 1.66053906660e-27
    n0 = args.central_density * 1e20
    t0 = math.sqrt(4e-7 * np.pi * n0 * args.central_mass * AMU)
    fld_prev = dt_norm = None
    if args.rect_prev:
        fld_prev, Rp_, Zp_, inside_p, t_prev = read_rectangle(args.rect_prev, exprs)
        if Rp_.size != R.size or Zp_.size != Z.size:
            sys.exit("--rect-prev must be sampled on the same rectangle")
        dt_norm = t_now - t_prev
        print(f"previous rectangle {args.rect_prev}: t = {t_prev:.4g}, dt = {dt_norm:.4g} JOREK units = {dt_norm * t0 * 1e3:.3f} ms")
    bq = c = pl = None
    if args.boundary:
        names, steps, times, data = read_fast(args.boundary)
        i = int(np.argmin(np.abs(times - t_now))) if args.step is None else int(np.where(steps == args.step)[0][0])
        bq = data[i]
        pl, c = plates(bq, names, z_max=args.plate_zmax)
        print(f"boundary {args.boundary}: step {steps[i]}, t = {times[i]:.4g} (dt to rectangle {times[i] - t_now:+.3g}),"
              f" {pl['inner'].size} inner / {pl['outer'].size} outer target points")
        wall_totals(bq, c, pl, xpoint_z(R, Z, fld, inside))

    bands = [tuple(float(x) for x in b.split(":")) for b in args.bands]
    nb = len(bands)
    fig, ax = plt.subplots(nb, 2, figsize=(11, 3.4 * nb), squeeze=False, constrained_layout=True)
    hdr = (f"{'band':>12s}{'L [m]':>7s}{'Ipar_in':>9s}{'Ipar_out':>9s}{'Ipol_in':>9s}{'Ipol_out':>9s}{'side_a':>9s}{'side_b':>9s}"
           f"{'resid':>9s}{'div':>9s}")
    if bq is not None:
        hdr += f"{'-I_in,pl':>10s}{'I_out,pl':>10s}{'sheath_in':>10s}{'sheath_out':>11s}"
    print(hdr + "   [kA, + from the inner to the outer end; plate currents in the band, Jn into the wall]")
    ends = []                                                                 # (a, b, tb, Iin, Iout, Sin, Sout) for --zoom
    for r, (a, b) in enumerate(bands):
        tb = tube_budget(R, Z, fld, inside, a, b, args.F0, args.points, args.gap, args.max_width, args.eta_si)
        if tb is not None and all(n in fld for n in ("Phi", "ne", "T_e")):
            tb["carriers"] = carrier_budget(R, Z, fld, inside, tb, args.F0, args.visco, args.visco_num, n0, t0,
                                            args.smooth, fld_prev, dt_norm)
        if tb is None:
            print(f"{a:.3f}:{b:.3f}   no open contour (private-region bands need the contour to stay below the X-point;"
                  " try a band further from 1, e.g. 0.975:0.99)")
            ax[r, 0].set_visible(False); ax[r, 1].set_visible(False)
            continue
        k = 1e-3
        Iin = Iout = Sin = Sout = np.nan
        line = (f"{a:6.3f}:{b:5.3f}{tb['L']:7.2f}{tb['I_par_in'] * k:9.3f}{tb['I_par_out'] * k:9.3f}"
                f"{tb['I_pol_in'] * k:9.3f}{tb['I_pol_out'] * k:9.3f}"
                f"{tb['side']['a']['cum'][-1] * k:9.3f}{tb['side']['b']['cum'][-1] * k:9.3f}"
                f"{tb['residual'][tb['win_out']].mean() * k:9.3f}{tb['div_cum'][tb['win_out']].mean() * k:9.3f}")
        if bq is not None:
            reg = "_pfr" if b < 1.0 else "_sol"
            Iin = plate_current(bq, c, pl["inner" + reg], "Jn_par", a, b, b < 1.0)
            Iout = plate_current(bq, c, pl["outer" + reg], "Jn_par", a, b, b < 1.0)
            Sin = plate_current(bq, c, pl["inner" + reg], "Jn_sheath", a, b, b < 1.0)
            Sout = plate_current(bq, c, pl["outer" + reg], "Jn_sheath", a, b, b < 1.0)
            line += f"{-Iin * k:10.3f}{Iout * k:10.3f}{-Sin * k:10.3f}{Sout * k:11.3f}"
        print(line)
        ends.append((a, b, tb, Iin, Iout, Sin, Sout))
        s, ell = tb["s"], tb["s"] * tb["L"]
        keep = (ell >= args.gap) & (ell <= tb["L"] - args.gap)                 # the last element next to each plate is left out
        m = lambda y: np.where(keep, y, np.nan)
        A = ax[r, 0]
        A.plot(ell, m(tb["I_par"]) * k, color=COLORS["par"], lw=1.8, label=r"$I_\parallel$ (parallel)")
        A.plot(ell, m(tb["I_pol"]) * k, color=COLORS["pol"], lw=1.4, ls="--", label=r"$I_{pol}$ (parallel + diamagnetic)")
        A.plot(ell, m(tb["I_pol"] + tb["I_side"]) * k, color=COLORS["res"], lw=1.4,
               label=r"$I_{pol}$ + current left through the sides")
        A.plot(ell, m(tb["I_pol_in"] + tb["div_cum"]) * k, color=COLORS["div"], lw=1.0, ls=":", label=r"$I_{pol}$(inner end) + $\int$div $J_{pol}$")
        if bq is not None and np.isfinite(Iin) and np.isfinite(Iout):
            A.plot([0.0], [-Iin * k], "v", color="k", ms=7, label="plate current (boundary output)")
            A.plot([tb["L"]], [Iout * k], "v", color="k", ms=7)
            if np.isfinite(Sin):
                A.plot([0.0, tb["L"]], [-Sin * k, Sout * k], "o", mfc="none", color="k", ms=7, label="sheath characteristic")
        for w in (tb["win_in"], tb["win_out"]):
            A.axvspan(ell[w].min(), ell[w].max(), color="0.9", zorder=0)
        A.axhline(0.0, color="0.6", lw=0.6)
        A.set_title(rf"$\Psi_N$ = {a:g} .. {b:g} ({'SOL' if a >= 1 else 'private region'}), L = {tb['L']:.2f} m", loc="left", fontsize=FS)
        A.set_ylabel("current [kA], + inner → outer", fontsize=FS)
        if r == 0:
            A.legend(frameon=False, fontsize=8)
        B = ax[r, 1]
        for key, col, lab in (("a", COLORS["a"], rf"side $\Psi_N$ = {a:g}"), ("b", COLORS["b"], rf"side $\Psi_N$ = {b:g}")):
            sd = tb["side"][key]
            kp = (sd["ell"] >= args.gap) & (sd["ell"] <= sd["L"] - args.gap)
            B.plot(sd["ell"], np.where(kp, sd["jn"], np.nan) * k, color=col, lw=1.4, label=lab)
        B.axhline(0.0, color="0.6", lw=0.6)
        B.set_ylabel(r"$J_{pol}\cdot n$ leaving the tube [kA/m$^2$]", fontsize=FS)
        B.legend(frameon=False, fontsize=8)
        if r == nb - 1:
            A.set_xlabel("poloidal distance from the inner end [m]", fontsize=FS)
            B.set_xlabel("poloidal distance from the inner end [m]", fontsize=FS)
    out = os.path.join(args.outdir, "current_conservation.pdf")
    fig.savefig(out)
    print(f"wrote {out}")
    composition(ends, bq, c, pl, args.sheath_lambda, args.outdir)
    if args.zoom > 0 and ends:
        plot_ends(ends, args.zoom, args.outdir)
    div_budget(R, Z, fld, inside, args.outdir)
    import matplotlib
    if matplotlib.get_backend().lower() not in ("agg", "pdf", "ps", "svg", "cairo", "template"):
        plt.show()


def composition(ends, bq, c, pl, lam, outdir):
    """Current composition per band in the taxonomy of Rozhansky et al. (NME 2020), see the module docstring.
    Plate currents from the boundary output when present (Jn into the wall, + = ions), else the tube's end values."""
    if not ends:
        return
    k = 1e-3
    print("\ncurrent composition per band [kA] (Rozhansky et al. 2020 taxonomy; + = ions into the plate for I_in, I_out;"
          " through > 0 = inner -> outer in the plasma):")
    print(f"{'band':>12s}{'I_in':>8s}{'I_out':>8s}{'through':>9s}{'closing':>9s} | {'dia<X':>8s}{'dia>X':>8s}{'res<X':>8s}{'res>X':>8s}"
          f"{'sum':>8s} | {'L_par':>7s}{'R_tube':>9s}{'LdTe':>7s}{'dPhi':>7s}{'I_th':>8s}")
    rows = []
    for a, b, tb, Iin, Iout, Sin, Sout in ends:
        if not (np.isfinite(Iin) and np.isfinite(Iout)):
            Iin, Iout = -tb["I_par_in"], tb["I_par_out"]                        # from the tube when no boundary file
        thr = 0.5 * (Iout - Iin)
        clo = Iin + Iout
        cp = tb["comp"]
        sides = -(cp["dia_below"] + cp["dia_above"]) + (cp["res_below"] + cp["res_above"])   # entering the tube
        LdT = dPhi = Ith = np.nan
        if bq is not None and "T_e" in c:
            reg = "_pfr" if b < 1.0 else "_sol"
            def band_mean(name, pn):
                g = bq.mean(axis=0); idx = pl[pn + reg]; psi = g[idx, c["Psi_N"]]
                lo, hi = (2.0 - b, 2.0 - a) if b < 1.0 else (a, b)
                sel = idx[(psi >= lo) & (psi <= hi)]
                return float(np.mean(g[sel, c[name]])) if sel.size else np.nan
            LdT = lam * (band_mean("T_e", "outer") - band_mean("T_e", "inner"))
            if "Phi" in c:
                dPhi = band_mean("Phi", "outer") - band_mean("Phi", "inner")
            if np.isfinite(cp.get("R_tube", np.nan)) and cp["R_tube"] > 0:
                Ith = -LdT / cp["R_tube"]            # A, in the through-current convention (+ inner -> outer): a hotter outer plate drives it outer -> inner
        print(f"{a:6.3f}:{b:5.3f}{Iin * k:8.3f}{Iout * k:8.3f}{thr * k:9.3f}{clo * k:9.3f} | {-cp['dia_below'] * k:8.3f}"
              f"{-cp['dia_above'] * k:8.3f}{cp['res_below'] * k:8.3f}{cp['res_above'] * k:8.3f}{sides * k:8.3f} | "
              f"{cp['L_par']:7.1f}{cp.get('R_tube', np.nan):9.2e}{LdT:7.1f}{dPhi:7.1f}{Ith * k:8.3f}")
        rows.append((a, b, Iin, Iout, thr, clo, cp, sides))
    if any("carriers" in tb for _, _, tb, *_ in ends):
        names = [k for k in ("diamagnetic", "viscous", "hyperviscous", "Reynolds", "polarisation d/dt") if any(k in tb.get("carriers", {}) for _, _, tb, *_ in ends)]
        print("\ncarriers of the vorticity equation, current ENTERING the tube [kA] (area integral of the model's div J per"
              " term; 'measured' = closing from the plates and the side residual of the postproc current):")
        print(f"{'band':>12s}" + "".join(f"{n:>18s}" for n in names) + f"{'sum':>9s}{'dI_par':>9s}{'closing':>9s}{'dia chk':>9s}")
        for a, b, tb, Iin, Iout, *_ in ends:
            cr = tb.get("carriers", {}); cp = tb["comp"]
            vals = [cr.get(n, np.nan) * k for n in names]
            dIpar = (tb["I_par_out"] - tb["I_par_in"]) * k                       # parallel current gained along the tube
            dia_ends = ((tb["I_pol_out"] - tb["I_par_out"]) - (tb["I_pol_in"] - tb["I_par_in"]))   # diamagnetic poloidal current out through the ends
            dia_chk = -(cp["dia_below"] + cp["dia_above"] + dia_ends) * k         # entering: -(sides + ends)
            print(f"{a:6.3f}:{b:5.3f}" + "".join(f"{v:18.3f}" for v in vals) + f"{np.nansum(vals):9.3f}{dIpar:9.3f}{(Iin + Iout) * k:9.3f}{dia_chk:9.3f}")
        print("  Each column is -int div(J_X) over the tube area between the end windows = the current carrier X deposits in"
              " the tube (through sides and ends). Their sum must equal dI_par, the parallel current gained between the"
              " windows (closing = the plates' version of it). 'dia chk' is the same for the diamagnetic current from the"
              " postproc's JR/JZ fluxes (sides + ends) and should match the diamagnetic column. Fourth derivatives of Phi:"
              " repeat with --smooth 2 and 4; a carrier that moves a lot between the two is not resolved.")
    print("  closing = I_in + I_out should equal 'sum' (side current entering the tube: diamagnetic = grad-B closing current"
          " of JR/JZ, res = polarisation/viscous/inertial, below/above the X-point); LdTe = Lambda (T_e,out - T_e,in) [V],"
          " dPhi = Phi_out - Phi_in [V] at the plates, I_th = -LdTe / R_tube (Spitzer unless --eta-si), same sign convention as through.")
    fig, ax = plt.subplots(figsize=(1.6 * len(rows) + 3, 4), constrained_layout=True)
    x = np.arange(len(rows))
    parts = [("dia_below", "grad-B closing, below X-point (PCC 1)", "#E08230"), ("dia_above", "grad-B closing, above X-point (PCC 2)", "#f2c28b"),
             ("res_below", "polarisation/viscous, below X-point", "#1f5fa8"), ("res_above", "polarisation/viscous, above X-point", "#9ec2e6")]
    pos = np.zeros(len(rows)); neg = np.zeros(len(rows))
    for key, lab, col in parts:
        sign = -1.0 if key.startswith("dia") else 1.0                          # side arrays count current LEAVING; residual counts entering
        v = np.array([sign * r[6][key] * k for r in rows])
        base = np.where(v >= 0, pos, neg)
        ax.bar(x, v, 0.6, bottom=base, color=col, label=lab)
        pos += np.maximum(v, 0); neg += np.minimum(v, 0)
    ax.plot(x, [r[5] * k for r in rows], "kv", ms=8, label="closing = I_in + I_out (plates)")
    ax.plot(x, [r[4] * k for r in rows], "ko", mfc="none", ms=8, label="through = (I_out - I_in)/2")
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels([f"{r[0]:g}..{r[1]:g}" for r in rows])
    ax.set_xlabel(r"$\Psi_N$ band", fontsize=FS); ax.set_ylabel("current [kA]", fontsize=FS)
    ax.legend(frameon=False, fontsize=8)
    out = os.path.join(outdir, "current_composition.pdf")
    fig.savefig(out)
    print(f"wrote {out}")


def div_budget(R, Z, fld, inside, outdir):
    """Where the current the postproc does not carry is: 2 pi int div(J_pol) R dA [kA] per region of the whole domain
    (grid points more than one cell from the boundary; div from finite differences of JR, JZ, element-scale noise
    averages out). In a steady state with div(J_pol + J_polarisation) = 0 the sum over the domain equals minus the
    polarisation/viscous current leaving through the boundary, i.e. it should match the Jnorm total into the wall
    (boundary output) up to the current through the outer flux surface. Also writes current_conservation_div.pdf."""
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    JR, JZ = np.where(inside, fld["JR"], np.nan), np.where(inside, fld["JZ"], np.nan)
    div = np.gradient(RR * JR, R, axis=0) / RR + np.gradient(JZ, Z, axis=1)      # NaN within one cell of the boundary
    dA = (R[1] - R[0]) * (Z[1] - Z[0])
    w = 2.0 * np.pi * div * RR * dA
    ok = np.isfinite(w)
    psiN = fld["Psi_N"]
    z_x = xpoint_z(R, Z, fld, inside)
    pfr = (psiN < 1.0) & (ZZ < z_x) & ok
    regions = [("core  Psi_N < 0.9", (psiN < 0.9) & ~pfr), ("core  0.9 - 1", (psiN >= 0.9) & (psiN < 1.0) & ~pfr),
               ("private region", pfr), ("SOL   1 - 1.02", (psiN >= 1.0) & (psiN < 1.02)),
               ("SOL   1.02 - 1.05", (psiN >= 1.02) & (psiN < 1.05)), ("SOL   1.05 - 1.1", (psiN >= 1.05) & (psiN < 1.1)),
               ("SOL   > 1.1", psiN >= 1.1)]
    shell = inside.copy()                                                     # points within ~10 mm of the boundary
    for _ in range(max(1, int(round(0.01 / (R[1] - R[0]))))):
        shell[1:-1, 1:-1] &= shell[:-2, 1:-1] & shell[2:, 1:-1] & shell[1:-1, :-2] & shell[1:-1, 2:]
        shell[0, :] = shell[-1, :] = shell[:, 0] = shell[:, -1] = False
    shell = inside & ~shell & ok
    print("current the postproc current does not conserve, 2 pi int div(J_pol) R dA per region [kA] (+ = J_pol diverges"
          " there, i.e. polarisation/viscous current converges); second column: the part within 10 mm of the boundary")
    tot = tot_sh = 0.0
    for name, m in regions:
        m = m & ok
        val, val_sh = np.nansum(w[m]) * 1e-3, np.nansum(w[m & shell]) * 1e-3
        tot += val
        tot_sh += val_sh
        print(f"    {name:18s} {val:9.3f} {val_sh:9.3f}   ({m.sum():7d} points, {m.sum() * dA:7.4f} m^2)")
    print(f"    {'whole domain':18s} {tot:9.3f} {tot_sh:9.3f}   (compare with the Jnorm total into the wall; one cell next"
          " to the boundary is excluded)")
    # the same total with the domain shrunk by k cells: the flux of J_pol through a surface k cells inside the boundary.
    # A total that changes with k says div(J_pol) is concentrated within that distance of the boundary (last element).
    h = 0.5 * ((R[1] - R[0]) + (Z[1] - Z[0]))
    line = "    flux of J_pol out of the domain shrunk by d from its boundary [kA]:"
    m = inside.copy()
    for k in range(1, 17):
        m[1:-1, 1:-1] &= m[:-2, 1:-1] & m[2:, 1:-1] & m[1:-1, :-2] & m[1:-1, 2:]
        m[0, :] = m[-1, :] = m[:, 0] = m[:, -1] = False
        if k in (1, 2, 4, 8, 16):
            JRk, JZk = np.where(m, fld["JR"], np.nan), np.where(m, fld["JZ"], np.nan)
            divk = np.gradient(RR * JRk, R, axis=0) / RR + np.gradient(JZk, Z, axis=1)
            line += f"  d = {k * h * 1e3:4.1f} mm: {np.nansum(2.0 * np.pi * divk * RR * dA) * 1e-3:8.3f}"
    print(line)
    fig, ax = plt.subplots(figsize=(6, 9), constrained_layout=True)
    fdiv = div * 1e-3
    lim = np.nanpercentile(np.abs(fdiv[ok]), 99)
    pc = ax.pcolormesh(R, Z, fdiv.T, cmap="RdBu_r", vmin=-lim, vmax=lim, shading="auto", rasterized=True)
    ax.contour(R, Z, np.ma.masked_where(~inside, psiN).T, levels=[0.9, 1.0, 1.02, 1.05, 1.1], colors="k", linewidths=0.5)
    ax.set_aspect("equal")
    ax.set_xlabel("R [m]", fontsize=FS)
    ax.set_ylabel("Z [m]", fontsize=FS)
    ax.set_title(r"div $J_{pol}$ (parallel + diamagnetic) [kA/m$^3$]", loc="left", fontsize=FS)
    fig.colorbar(pc, ax=ax, shrink=0.6)
    out = os.path.join(outdir, "current_conservation_div.pdf")
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


def plot_ends(ends, zoom, outdir):
    """The last `zoom` m next to each plate for every band, without the gap mask: the parallel current against the
    plate current (does the sheath current arrive as parallel current, or is it fed across the field in the last
    element?) and the side current densities with the poloidal distance to the plate."""
    k = 1e-3
    nb = len(ends)
    fig, ax = plt.subplots(2 * nb, 2, figsize=(11, 4.6 * nb), squeeze=False, constrained_layout=True)
    for r, (a, b, tb, Iin, Iout, Sin, Sout) in enumerate(ends):
        ell = tb["s"] * tb["L"]
        for col, (name, d, Ipl, Spl) in enumerate((("inner", ell, -Iin, -Sin), ("outer", tb["L"] - ell, Iout, Sout))):
            A, B = ax[2 * r, col], ax[2 * r + 1, col]
            sel = d <= zoom
            A.plot(d[sel], tb["I_par"][sel] * k, color=COLORS["par"], lw=1.8, label=r"$I_\parallel$ (parallel)")
            A.plot(d[sel], tb["I_pol"][sel] * k, color=COLORS["pol"], lw=1.4, ls="--", label=r"$I_{pol}$ (parallel + diamagnetic)")
            A.plot(d[sel], (tb["I_pol"] + tb["I_side"])[sel] * k, color=COLORS["res"], lw=1.4,
                   label=r"$I_{pol}$ + current left through the sides")
            if np.isfinite(Ipl):
                A.plot([0.0], [Ipl * k], "v", color="k", ms=7, label="plate current (boundary output)")
            if np.isfinite(Spl):
                A.plot([0.0], [Spl * k], "o", mfc="none", color="k", ms=7, label="sheath characteristic")
            A.axhline(0.0, color="0.6", lw=0.6)
            A.set_title(rf"$\Psi_N$ = {a:g} .. {b:g}, {name} end", loc="left", fontsize=FS)
            A.set_ylabel("current [kA], + inner → outer", fontsize=FS)
            if r == 0 and col == 0:
                A.legend(frameon=False, fontsize=8)
            for key, colr, lab in (("a", COLORS["a"], rf"side $\Psi_N$ = {a:g}"), ("b", COLORS["b"], rf"side $\Psi_N$ = {b:g}")):
                sd = tb["side"][key]
                dk = sd["ell"] if name == "inner" else sd["L"] - sd["ell"]
                sk = dk <= zoom
                B.plot(dk[sk], sd["jn"][sk] * k, color=colr, lw=1.4, label=lab)
            B.axhline(0.0, color="0.6", lw=0.6)
            B.set_ylabel(r"$J_{pol}\cdot n$ leaving the tube [kA/m$^2$]", fontsize=FS)
            B.set_xlabel(f"poloidal distance from the {name} plate [m]", fontsize=FS)
            if r == 0:
                B.legend(frameon=False, fontsize=8)
    out = os.path.join(outdir, "current_conservation_ends.pdf")
    fig.savefig(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
