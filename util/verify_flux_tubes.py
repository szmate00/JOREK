#!/usr/bin/env python3
"""Ohm's law along the SOL flux tubes of an axisymmetric run, from jorek2_postproc `rectangle` output (SI units).

  python3 util/verify_flux_tubes.py exprs_Rmin..._s01500.h5 --thermo-coef 0.71 [--psi 1.01 1.20 40]

For n = 0 the parallel Ohm's law of model600 with the thermal force, integrated along a poloidal flux contour from
the inner to the outer end of a SOL flux tube, reads (Phi = F0 u in volts, (R, Z, phi) right-handed as in JOREK)

  Phi_out - Phi_in =   int (B_phi/B_p^2) (E_phi - eta j_phi) B_p.dl                     resistive + inductive
                     + int (B_phi/B_p^2) div(eta_num grad(R j_phi)) B_p.dl              hyper-resistivity (--eta-num)
                     + int (B_phi^2/B^2) dp_e/(e n_e)                                   electron pressure
                     + c int (B_phi^2/B^2) dT_e/e                                       thermal force, c = thermoelectric_coef

This is the induction equation as assembled in models/model600/mod_elt_matrix_fft.f90 for n = 0, without current
source and bootstrap current (keep_current_prof = .f., bootstrap = .f.) and without tgnum_psi:

  (1/R) dpsi/dt = eta zj/R + [psi, u] - div(eta_num grad zj)
                  - 2 tauIC F0^2/(rho B^2 R^2) [psi, p_e] - 2 tauIC c F0^2/(B^2 R^2) [psi, T_e],     zj = -R mu0 j_phi,

with [a, b] = a_R b_Z - a_Z b_R = -R B_p.grad b and 2 tauIC F0 = 1/e in normalised units (tauIC at its physical value).
The pressure term is divided by the corrected density corr_neg_dens(rho), which differs from rho below rho_1 (or
rho_min_neg): --dens-corr-scale. eta_num is constant unless eta_num_T_dependent or eta_num_psin_dependent is set (then --eta-num is not enough).

with E_phi = -(dPsi/dt)/(R t_0), j_phi = `currdens`, and B_p.dl the projection of the path on the poloidal field.
The script evaluates both sides on Psi_N contours and plots them: the right-hand side with the thermal force
(--thermo-coef 0.71: runs with thermoelectric_ohm) and without (0: runs without). A run verifies the
implementation if its left-hand side follows the curve with its own coefficient.

The contours end one grid cell before the plates (no extrapolation), so both sides refer to the same end points.
The HDF5 file keeps only the first letter of each expression name; the order is taken from --exprs, which must be
the `expressions` line used for the `rectangle` command (default: doc/elm_sheath_thermo/postproc/rectangle.pp).
"""
import argparse
import math
import sys

import numpy as np
import matplotlib.pyplot as plt

MU0 = 4.0e-7 * np.pi
AMU = 1.66053906660e-27
E_CHG = 1.602176634e-19
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
DEFAULT_EXPRS = "R Z Psi_N Phi T_e T_i ne pres_e currdens eta_T dPsi_dt BR BZ Btor"
trapezoid = getattr(np, "trapezoid", None) or np.trapz



def show_or_save(prefix):
    """plt.show() on an interactive backend; otherwise (no display) save every open figure as <prefix>_<n>.png."""
    import matplotlib
    if matplotlib.get_backend().lower() in ("agg", "pdf", "ps", "svg", "cairo", "template"):
        for n in plt.get_fignums():
            name = f"{prefix}_{n}.png"
            plt.figure(n).savefig(name, dpi=150)
            print(f"no display: saved {name}")
    else:
        plt.show()

def read_rectangle(path, names):
    """Fields on the rectangle as {name: array[iR, iZ]}, the axes R[iR], Z[iZ] and the mask of points in the domain."""
    import h5py
    with h5py.File(path, "r") as f:
        v = np.array(f["values"])
        t_now = float(np.array(f["t_now"]).ravel()[0]) if "t_now" in f else float("nan")
    ax_e = [k for k, n in enumerate(v.shape) if n == len(names)]
    if len(ax_e) != 1:
        sys.exit(f"{path}: cannot find the expression axis (shape {v.shape}, {len(names)} names); check --exprs")
    v = np.moveaxis(v, ax_e[0], -1)
    fld = {n: np.array(v[..., k], dtype=float) for k, n in enumerate(names)}
    inside = fld["Btor"] != 0.0

    def mean(a, ax):                                                                # mean over the points in the domain
        n = inside.sum(axis=ax)
        return np.where(n > 0, np.where(inside, a, 0.0).sum(axis=ax) / np.maximum(n, 1), np.nan)

    if np.nanstd(mean(fld["R"], 1)) < np.nanstd(mean(fld["R"], 0)):                 # R varies along axis 1: swap
        fld = {n: a.T for n, a in fld.items()}
        inside = inside.T

    # points inside the domain that the sampler missed (zeros surrounded by valid points) would cut the contours:
    # fill them from their valid neighbours
    n_filled = 0
    for _ in range(4):
        pad = np.pad(inside, 1)
        nb = sum(pad[1 + di:pad.shape[0] - 1 + di, 1 + dj:pad.shape[1] - 1 + dj]
                 for di in (-1, 0, 1) for dj in (-1, 0, 1) if (di, dj) != (0, 0))
        hole = (~inside) & (nb >= 6)
        if not hole.any():
            break
        for n, a in fld.items():
            pa = np.pad(np.where(inside, a, 0.0), 1)
            tot = sum(pa[1 + di:pa.shape[0] - 1 + di, 1 + dj:pa.shape[1] - 1 + dj]
                      for di in (-1, 0, 1) for dj in (-1, 0, 1) if (di, dj) != (0, 0))
            a[hole] = tot[hole] / nb[hole]
        inside = inside | hole
        n_filled += int(hole.sum())
    if n_filled:
        print(f"filled {n_filled} isolated points inside the domain that the sampler had left empty")

    def axis(a, ax):
        m = mean(a, ax)
        k = np.arange(m.size)
        ok = np.isfinite(m)
        return np.polyval(np.polyfit(k[ok], m[ok], 1), k)                           # uniform; fills empty rows

    return fld, axis(fld["R"], 1), axis(fld["Z"], 0), inside, t_now


def bilinear(R, Z, a, pts):
    """a[iR, iZ] on the uniform axes R, Z at the points pts (n, 2)."""
    x = (pts[:, 0] - R[0]) / (R[1] - R[0])
    y = (pts[:, 1] - Z[0]) / (Z[1] - Z[0])
    i = np.clip(np.floor(x).astype(int), 0, R.size - 2)
    j = np.clip(np.floor(y).astype(int), 0, Z.size - 2)
    fx, fy = x - i, y - j
    return (a[i, j] * (1 - fx) * (1 - fy) + a[i + 1, j] * fx * (1 - fy) + a[i, j + 1] * (1 - fx) * fy
            + a[i + 1, j + 1] * fx * fy)


def flux_tube(R, Z, fld, inside, level, t_unit, coef, log_density=True, n_corr_scale=0.0):
    """Terms of the integrated Ohm's law on the longest open contour Psi_N = level, from its inner to its outer end."""
    psi = np.ma.masked_where(~inside, fld["Psi_N"])
    fig = plt.figure()
    # corner_mask=False: contour only through cells whose four corners are inside the domain. With the default, a cell
    # with one corner outside is cut into a triangle and the contour enters it; the bilinear interpolation would then
    # mix in the empty corner (a spurious jump of every field at the two ends of the path)
    segs = fig.gca().contour(R, Z, psi.T, levels=[level], corner_mask=False).allsegs[0]
    plt.close(fig)
    segs = [s for s in segs if len(s) > 10]
    if not segs:
        return None
    path = max(segs, key=lambda s: np.sum(np.hypot(*np.diff(s, axis=0).T)))
    if np.allclose(path[0], path[-1]):
        return None                                                    # closed contour: not a SOL flux tube
    n_seg = len(segs)
    if path[0, 0] > path[-1, 0]:
        path = path[::-1]                                              # start at the inner (small R) end
    g = {n: bilinear(R, Z, fld[n], path) for n in ("Phi", "T_e", "ne", "pres_e", "currdens", "eta_T", "dPsi_dt",
                                                    "BR", "BZ", "Btor", "hyper", "ln_ne")}
    dl = np.diff(path, axis=0)
    mid = lambda a: 0.5 * (a[1:] + a[:-1])
    Rm = mid(path[:, 0])
    bp_dl = mid(g["BR"]) * dl[:, 0] + mid(g["BZ"]) * dl[:, 1]
    geom = mid(g["Btor"]) / (mid(g["BR"]) ** 2 + mid(g["BZ"]) ** 2) * bp_dl
    E_phi = -mid(g["dPsi_dt"]) / (Rm * t_unit)
    seg_res = -geom * mid(g["eta_T"]) * mid(g["currdens"])
    seg_ind = geom * E_phi
    seg_hyp = geom * mid(g["hyper"])
    # the pressure and thermal-force terms of the psi equation carry the factor B_phi^2/B^2 (mod_elt_matrix_fft)
    bphi2 = mid(g["Btor"]) ** 2 / (mid(g["Btor"]) ** 2 + mid(g["BR"]) ** 2 + mid(g["BZ"]) ** 2)
    if log_density:                                    # dp_e/(e n_e) = T_e dln(n_e) + dT_e  [V, T_e in eV]
        seg_pe = bphi2 * (mid(g["T_e"]) * np.diff(g["ln_ne"]) + np.diff(g["T_e"]))
    else:
        seg_pe = bphi2 * np.diff(g["pres_e"]) / (E_CHG * mid(g["ne"]))
    if n_corr_scale > 0.0:
        # JOREK divides the pressure gradient by the CORRECTED density corr_neg_dens(rho) (models/corr_neg_include.f90):
        # n_corr = n above the scale, L1 + L2 exp((n - L1 - L2)/L2) below it, L1 = L2 = scale/2 (default coefficients).
        # The scale is rho_min_neg if that is set (>= 0), otherwise rho_1 of the density profile.
        n, half = mid(g["ne"]), 0.5 * n_corr_scale
        n_corr = np.where(n < 2.0 * half, half + half * np.exp((n - 2.0 * half) / half), n)
        seg_pe = seg_pe * n / n_corr
    seg_T = bphi2 * np.diff(g["T_e"])
    cum = lambda a: np.concatenate([[0.0], np.cumsum(a)])
    bp = np.hypot(mid(g["BR"]), mid(g["BZ"]))
    turns = np.sum(np.abs(mid(g["Btor"])) / (Rm * bp) * np.hypot(dl[:, 0], dl[:, 1])) / (2.0 * np.pi)
    return {
        "turns": turns, "v_loop": float(np.mean(2.0 * np.pi * Rm * E_phi)),
        "n_seg": n_seg, "path": path, "l": cum(np.hypot(dl[:, 0], dl[:, 1])),
        "cum_lhs": g["Phi"] - g["Phi"][0], "cum_res": cum(seg_res), "cum_ind": cum(seg_ind), "cum_pe": cum(seg_pe),
        "cum_hyp": cum(seg_hyp),
        "cum_thermo": coef * cum(seg_T), "cum_T": cum(seg_T), "Te_path": g["T_e"],
        "lhs": g["Phi"][-1] - g["Phi"][0],
        "res": np.sum(seg_res),
        "ind": np.sum(geom * E_phi), "hyp": np.sum(seg_hyp),
        "pe": np.sum(seg_pe),
        "dTe": g["T_e"][-1] - g["T_e"][0],
        "thermo": coef * np.sum(seg_T),
        "Te_in": g["T_e"][0], "Te_out": g["T_e"][-1], "length": np.sum(np.hypot(dl[:, 0], dl[:, 1])),
        "ends": (path[0], path[-1]),
    }


def analyse(path, exprs, levels, coef, central_density=1.0, central_mass=2.01410177811, max_end_height=0.15,
            eta_num=0.0, log_density=True, dens_corr_scale=0.0):
    t_unit = math.sqrt(MU0 * central_density * 1e20 * central_mass * AMU)
    fld, R, Z, inside, t_now = read_rectangle(path, exprs)
    # hyper-resistivity: div(eta_num grad(R j_phi)) [V/m] with constant eta_num (JOREK units -> SI: mu0/t_unit)
    S = fld["R"] * fld["currdens"]
    ok = inside.copy()
    ok[1:-1, 1:-1] &= inside[2:, 1:-1] & inside[:-2, 1:-1] & inside[1:-1, 2:] & inside[1:-1, :-2]
    ok[0, :] = ok[-1, :] = ok[:, 0] = ok[:, -1] = False
    lap = np.zeros_like(S)
    lap[1:-1, 1:-1] = ((S[2:, 1:-1] - 2.0 * S[1:-1, 1:-1] + S[:-2, 1:-1]) / (R[1] - R[0]) ** 2
                       + (S[1:-1, 2:] - 2.0 * S[1:-1, 1:-1] + S[1:-1, :-2]) / (Z[1] - Z[0]) ** 2)
    fld["hyper"] = np.where(ok, eta_num * MU0 / t_unit * lap, 0.0)
    # ln(n_e) for the pressure term: the density varies by orders of magnitude towards the inner plate in the far
    # SOL, where interpolating n_e and p_e linearly between the sampling points is not accurate enough
    fld["ln_ne"] = np.log(np.maximum(fld["ne"], 1e10))
    z_low = Z[np.where(inside.any(axis=0))[0][0]]                      # lowest point of the domain: the plates are there
    out, skipped = [], []
    for lev in levels:
        tube = flux_tube(R, Z, fld, inside, lev, t_unit, coef, log_density, dens_corr_scale * central_density * 1e20)
        if tube is None:
            continue
        tube["psi"] = lev
        ends = tube["ends"]
        if max(ends[0][1], ends[1][1]) > z_low + max_end_height or abs(ends[0][0] - ends[1][0]) < 0.1:
            skipped.append(tube)                                       # does not run from plate to plate
        else:
            out.append(tube)
    geom = {"R": R, "Z": Z, "inside": inside, "skipped": skipped}
    if not out:
        print(f"{path}: no flux tube accepted of {len(levels)} levels (lowest domain point Z = {z_low:.3f} m, "
              f"ends allowed up to {max_end_height:.2f} m above it); the longest open contour pieces end at:")
        for t in skipped[:8]:
            e = t["ends"]
            print(f"  Psi_N = {t['psi']:.3f}: ({e[0][0]:.3f}, {e[0][1]:.3f}) -> ({e[1][0]:.3f}, {e[1][1]:.3f})")
        if not skipped:
            print("  (no open contour found at all: check --exprs and that Psi_N is in the file)")
        print("  -> raise --end-height if the ends are on the targets, or narrow the Psi_N range to tubes that reach them")
    return out, t_now * t_unit * 1e3, geom


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file")
    ap.add_argument("--thermo-coef", type=float, required=True,
                    help="thermoelectric_coef of the run if thermoelectric_ohm = .t. (0.71), else 0")
    ap.add_argument("--psi", type=float, nargs=3, default=[1.01, 1.20, 40], metavar=("MIN", "MAX", "N"))
    ap.add_argument("--exprs", default=DEFAULT_EXPRS)
    ap.add_argument("--reference", help="rectangle file of the run without the thermal force (same Psi_N levels): "
                    "the coefficient is also fitted on the difference of the two runs, which removes what they share")
    ap.add_argument("--eta-num", type=float, default=0.0,
                    help="eta_num of the namelist (JOREK units, constant): adds the hyper-resistive term")
    ap.add_argument("--dens-corr-scale", type=float, default=0.0,
                    help="scale of JOREK's density correction in the pressure term, JOREK units: rho_min_neg if it is "
                         "set (>= 0) in the namelist, otherwise rho_1 (0.01 for inxflow); 0 = no correction")
    ap.add_argument("--linear-density", action="store_true",
                    help="pressure term from linearly interpolated p_e and n_e (default: T_e dln(n_e) + dT_e)")
    ap.add_argument("--trim", type=float, default=0.10,
                    help="poloidal distance [m] from each end excluded from the comparison and the fit")
    ap.add_argument("--end-height", type=float, default=0.15,
                    help="a contour counts as a flux tube if both ends lie within this height [m] of the lowest point")
    ap.add_argument("--central-density", type=float, default=1.0)
    ap.add_argument("--central-mass", type=float, default=2.01410177811)
    args = ap.parse_args()

    levels = np.linspace(args.psi[0], args.psi[1], int(args.psi[2]))
    tubes, t_ms, geom = analyse(args.file, args.exprs.split(), levels, args.thermo_coef, args.central_density,
                                args.central_mass, args.end_height, args.eta_num, not args.linear_density, args.dens_corr_scale)
    for t in geom["skipped"]:
        print(f"skipped Psi_N = {t['psi']:.3f}: longest contour piece runs from ({t['ends'][0][0]:.3f}, {t['ends'][0][1]:.3f}) to "
              f"({t['ends'][1][0]:.3f}, {t['ends'][1][1]:.3f}), {t['n_seg']} pieces at this level")
    if not tubes:
        sys.exit("no contour running from plate to plate found in the given Psi_N range")
    # The potential changes by ~100 V within the last centimetres in front of the plates, and the psi equation is
    # replaced by its Dirichlet condition on the plate nodes, so end-to-end differences depend on where exactly the
    # contour stops. Both sides are therefore compared between the points --trim metres (poloidal) from each end, and
    # the thermal-force coefficient is fitted along the path, where T_e varies by tens of eV:
    #   Phi(l) - [resistive + inductive + pressure](l) = a + c (T_e(l) - T_e(0))
    psi = np.array([t["psi"] for t in tubes])
    lhs, base, dTe, c_fit, Te_in, Te_out, Te_up, parts = [], [], [], [], [], [], [], {"res": [], "ind": [], "pe": [], "hyp": []}
    c_diff, ref = [], None
    if args.reference:
        tubes0, _, _ = analyse(args.reference, args.exprs.split(), levels, 0.0, args.central_density, args.central_mass,
                               args.end_height, args.eta_num, not args.linear_density, args.dens_corr_scale)
        ref = {round(t["psi"], 6): t for t in tubes0}
    for t in tubes:
        m = (t["l"] >= args.trim) & (t["l"] <= t["l"][-1] - args.trim)
        i0, i1 = np.where(m)[0][[0, -1]]
        d = lambda key: t[key][i1] - t[key][i0]
        lhs.append(d("cum_lhs"))
        base.append(d("cum_res") + d("cum_ind") + d("cum_pe") + d("cum_hyp"))
        for key in parts:
            parts[key].append(d("cum_" + key))
        dTe.append(d("cum_T"))
        Te_in.append(t["Te_path"][i0]); Te_out.append(t["Te_path"][i1]); Te_up.append(t["Te_path"][m].max())
        resid = (t["cum_lhs"] - t["cum_res"] - t["cum_ind"] - t["cum_pe"] - t["cum_hyp"])[m]
        x = t["cum_T"][m]
        c_fit.append(np.polyfit(x, resid, 1)[0] if np.ptp(x) > 1.0 else np.nan)
        if ref is not None:                                            # same fit on the difference to the reference run
            t0 = ref.get(round(t["psi"], 6))
            if t0 is None:
                c_diff.append(np.nan)
            else:
                r0 = t0["cum_lhs"] - t0["cum_res"] - t0["cum_ind"] - t0["cum_pe"] - t0["cum_hyp"]
                r0 = np.interp(t["l"][m] / t["l"][-1], t0["l"] / t0["l"][-1], r0)
                c_diff.append(np.polyfit(x, resid - r0, 1)[0] if np.ptp(x) > 1.0 else np.nan)
    lhs, base, dTe, c_fit, c_diff = (np.array(a) for a in (lhs, base, dTe, c_fit, c_diff))
    rhs = base + args.thermo_coef * dTe

    print(f"{args.file}: t = {t_ms:.3f} ms, {len(tubes)} flux tubes, compared between the points {args.trim} m from each end")
    print(f"{'Psi_N':>7s}{'Te_in':>8s}{'Te_up':>8s}{'Te_out':>8s}{'LHS':>9s}{'resist.':>9s}{'induct.':>9s}{'hyper':>9s}{'pe':>9s}{'thermo':>9s}"
          f"{'RHS':>9s}{'LHS-RHS':>9s}{'c fit':>8s}   [V; eV]")
    for k, t in enumerate(tubes):
        print(f"{t['psi']:7.3f}{Te_in[k]:8.1f}{Te_up[k]:8.1f}{Te_out[k]:8.1f}{lhs[k]:9.2f}{parts['res'][k]:9.2f}{parts['ind'][k]:9.2f}"
              f"{parts['hyp'][k]:9.2f}{parts['pe'][k]:9.2f}{args.thermo_coef * dTe[k]:9.2f}{rhs[k]:9.2f}{lhs[k] - rhs[k]:9.2f}{c_fit[k]:8.3f}")
    print(f"\nrms(LHS - RHS) = {np.sqrt(np.mean((lhs - rhs) ** 2)):.2f} V, rms(LHS) = {np.sqrt(np.mean(lhs ** 2)):.2f} V")
    print(f"thermal-force coefficient fitted along the paths: median {np.nanmedian(c_fit):.3f}, "
          f"quartiles {np.nanpercentile(c_fit, 25):.3f} .. {np.nanpercentile(c_fit, 75):.3f} (run value: {args.thermo_coef})")
    if ref is not None:
        print(f"the same fit on the difference to {args.reference}: median {np.nanmedian(c_diff):.3f}, quartiles "
              f"{np.nanpercentile(c_diff, 25):.3f} .. {np.nanpercentile(c_diff, 75):.3f}, range {np.nanmin(c_diff):.3f} .. "
              f"{np.nanmax(c_diff):.3f}")
        print(f"{'Psi_N':>7s}{'c':>8s}{'c - c_ref':>11s}")
        for k, t in enumerate(tubes):
            print(f"{t['psi']:7.3f}{c_fit[k]:8.3f}{c_diff[k]:11.3f}")

    # --- the induced toroidal field (decaying plasma current, no current source) as a drive of the SOL current
    print("\ninduced field along the flux tubes (whole tube, plate to plate):")
    print(f"{'Psi_N':>7s}{'turns':>8s}{'V_loop':>9s}{'induct.':>9s}{'resist.':>9s}{'ind/res':>9s}   [V]  "
          "ind/res > 0: the induced field opposes the tube current, < 0: it drives it")
    for t in tubes:
        ratio = t["ind"] / t["res"] if t["res"] != 0 else float("nan")
        print(f"{t['psi']:7.3f}{t['turns']:8.2f}{t['v_loop']:9.2f}{t['ind']:9.2f}{t['res']:9.2f}{ratio:9.2f}")

    fig, ax = plt.subplots(3, 1, figsize=(8, 10), sharex=True, constrained_layout=True)
    ax[0].plot(psi, lhs, "o", ms=4, color=INK, label=r"$\Delta\Phi$ (simulation)")
    ax[0].plot(psi, base + 0.71 * dTe, color=C1, label="Ohm's law with the thermal force (c = 0.71)")
    ax[0].plot(psi, base, color=C2, label="Ohm's law without it (c = 0)")
    ax[0].set(ylabel="V", title=f"Potential difference between the points {args.trim} m from the plates, t = {t_ms:.2f} ms")
    ax[0].legend(loc="best")
    ax[1].plot(psi, c_fit, "o", ms=4, color=INK, label="fitted along each flux tube")
    if ref is not None:
        ax[1].plot(psi, c_diff, "s", ms=4, mfc="white", mec=C3, label="fitted on the difference to the reference run")
    ax[1].axhline(0.71, color=C1, lw=1.2, label="0.71")
    ax[1].axhline(0.0, color=C2, lw=1.2, label="0")
    ax[1].set(ylabel="c", title="Thermal-force coefficient from $\\Phi(l)$ against $T_e(l)$ along the path")
    ax[1].legend(loc="best")
    for key, c, lab in (("res", C1, "resistive"), ("ind", C3, "inductive"), ("pe", C2, "electron pressure")):
        ax[2].plot(psi, parts[key], color=c, label=lab)
    if args.eta_num:
        ax[2].plot(psi, parts["hyp"], color=INK, lw=1.0, label="hyper-resistivity")
    ax[2].plot(psi, 0.71 * dTe, color=INK2, ls=(0, (4, 2)), label=r"$0.71\,\Delta T_e/e$")
    ax[2].axhline(0.0, color=MUTED, lw=0.8)
    ax[2].set(ylabel="V", xlabel=r"$\Psi_N$", title="Terms of the right-hand side")
    ax[2].legend(loc="best")
    for a in ax:
        a.grid(True, color="#e4e3df", lw=0.6)

    # --- where the contours run, and both sides accumulated along three of them
    pick = [tubes[k] for k in sorted({0, len(tubes) // 2, len(tubes) - 1})]
    fig = plt.figure(figsize=(13, 8), constrained_layout=True)
    gs = fig.add_gridspec(len(pick), 2, width_ratios=[1.0, 1.4])
    am = fig.add_subplot(gs[:, 0])
    am.contourf(geom["R"], geom["Z"], geom["inside"].T.astype(float), levels=[0.5, 1.5], colors=["#e9e8e4"])
    for t in tubes:
        am.plot(t["path"][:, 0], t["path"][:, 1], color=MUTED, lw=0.5)
        am.plot(*t["ends"][0], "o", ms=3, color=C1)
        am.plot(*t["ends"][1], "o", ms=3, color=C2)
    for t in geom["skipped"]:
        am.plot(t["path"][:, 0], t["path"][:, 1], color="#e34948", lw=1.0)
    for t in pick:
        am.plot(t["path"][:, 0], t["path"][:, 1], color=INK, lw=1.2)
    am.set_aspect("equal")
    am.set(xlabel="R [m]", ylabel="Z [m]", title="contours used (ends: blue inner, orange outer; red: skipped)")
    for k, t in enumerate(pick):
        a = fig.add_subplot(gs[k, 1])
        tot = t["cum_res"] + t["cum_ind"] + t["cum_pe"] + t["cum_hyp"] + t["cum_thermo"]
        a.plot(t["l"], t["cum_lhs"], color=INK, lw=2.0, label=r"$\Phi(l)-\Phi_\mathrm{in}$ (simulation)")
        a.plot(t["l"], tot, color=C1, lw=1.5, label=f"Ohm's law, c = {args.thermo_coef}")
        a.plot(t["l"], t["cum_pe"], color=C2, lw=0.9, label="electron pressure")
        a.plot(t["l"], t["cum_res"] + t["cum_ind"], color=C3, lw=0.9, label="resistive + inductive")
        if args.eta_num:
            a.plot(t["l"], t["cum_hyp"], color=MUTED, lw=0.9, label="hyper-resistivity")
        a.axhline(0.0, color=MUTED, lw=0.8)
        a.grid(True, color="#e4e3df", lw=0.6)
        a.set(ylabel="V", title=rf"accumulated from the inner end, $\Psi_N$ = {t['psi']:.3f}")
        if k == 0:
            a.legend(loc="best", fontsize=8)
    a.set_xlabel("poloidal distance along the contour [m]")
    show_or_save("flux_tubes")


if __name__ == "__main__":
    main()
