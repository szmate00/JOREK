#!/usr/bin/env python3
"""
analyse_flux_function.py -- is Phi a flux function, and does the departure move density?

Reads the HDF5 written by `rectangle` in jorek2_postproc (see postproc_flux_function.in)
and forms

    v_psi = (v_ExB . grad psi)/|grad psi| = ( V_ExB_Z*BR - V_ExB_R*BZ ) / Btheta   [m/s]

which is IDENTICALLY ZERO if and only if Phi = Phi(psi). Nothing is differentiated
numerically: grad psi is reconstructed from the exported B components, since
B_pol = (psi_Z, -psi_R)/R gives grad psi = (-R*BZ, R*BR) and |grad psi| = R*Btheta, and
the R factors cancel.

Reported alongside it:

    sin_ang = v_psi / |v_ExB|      dimensionless, 0 = flux function, +-1 = purely radial
    Gamma_psi     = ne * v_psi                    cross-field particle flux density
    Gamma_par_pol = ne * vpar * Btheta/B_abs      parallel flux projected poloidally

Their ratio is the number the whole exercise is for: how much the ExB drags ACROSS flux
surfaces relative to what streams ALONG them, both as poloidal-plane flux densities.

Usage
-----
    analyse_flux_function.py FILE.h5                    one case
    analyse_flux_function.py REF.h5 TEST.h5             A/B, with the difference

    --pfr Rmin Rmax Zmin Zmax   PFR box (default: AUG 38773)
    --labels A B                names for the two cases
    --plot FILE.png             maps of v_psi, sin_ang and (for an A/B) the difference
    --vars "R Z Psi_N ..."      expression order, if you edited the .in file
"""

import argparse
import math
import os
import sys

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("ERROR: h5py is required.  pip install h5py")

# Must match the `expressions` line of postproc_flux_function.in. The names stored IN the
# file are useless: write_hdf5_2d (mod_diag_output.f90:552) copies each name into a
# character(len=1) array, so every one is truncated to its first letter. Order is the only
# reliable key, and a mismatch here silently mislabels every column -- hence the checks in
# load().
DEFAULT_VARS = ("R Z Psi_N Phi ne T_e T_i vpar V_ExB_R V_ExB_Z BR BZ Btheta B_abs").split()

# AUG 38773: strike points (1.256,-1.055) and (1.593,-1.150), X-point near (1.43,-0.93).
# The box sits inside both legs and below the X-point.
DEFAULT_PFR = (1.30, 1.52, -1.14, -1.00)


def load(path, varnames):
    with h5py.File(path, "r") as f:
        vals = np.array(f["values"])          # h5py returns C order: (n_var, nZ, nR)
        n_var = int(np.array(f["n_var"]))
        dim = np.array(f["dim"]).astype(int)  # Fortran (nR, nZ)
        t_now = float(np.array(f["t_now"]))
        time = float(np.array(f["time"])) if "time" in f else float("nan")
        idx = int(np.array(f["index_now"])) if "index_now" in f else -1

    if n_var != len(varnames):
        raise SystemExit(
            "ERROR: %s holds %d variables but %d names were given.\n"
            "       The file stores no usable names (they are truncated to one letter by\n"
            "       write_hdf5_2d), so the list must match the .in file exactly.\n"
            "       Pass the right one with --vars." % (path, n_var, len(varnames)))

    if vals.shape[0] != n_var:
        # Be forgiving about which axis h5py put first.
        if vals.shape[-1] == n_var:
            vals = np.moveaxis(vals, -1, 0)
        else:
            raise SystemExit("ERROR: %s has values shape %s, incompatible with n_var=%d."
                             % (path, vals.shape, n_var))

    d = {name: vals[i] for i, name in enumerate(varnames)}
    d["_t_now"], d["_time"], d["_index"] = t_now, time, idx
    d["_shape"], d["_dim"], d["_path"] = vals.shape[1:], tuple(dim), path

    # The R and Z columns are exported precisely so the layout can be verified rather than
    # assumed: R must vary along one axis and Z along the other.
    R, Z = d["R"], d["Z"]
    dR = (np.nanmax(np.abs(np.diff(R, axis=1))), np.nanmax(np.abs(np.diff(R, axis=0))))
    if dR[0] < dR[1]:
        # R varies down axis 0 -> the grid is transposed relative to what we assume.
        for k in varnames:
            d[k] = d[k].T
        d["_shape"] = d["R"].shape
    return d


def derive(d):
    """Everything the report needs, computed once."""
    vR, vZ = d["V_ExB_R"], d["V_ExB_Z"]
    BR, BZ, Bth, Babs = d["BR"], d["BZ"], d["Btheta"], d["B_abs"]

    with np.errstate(divide="ignore", invalid="ignore"):
        # v_ExB . grad(psi) / |grad(psi)|   -- the R factors cancel exactly.
        d["v_psi"] = np.where(Bth > 0, (vZ * BR - vR * BZ) / Bth, np.nan)
        d["vE_mag"] = np.hypot(vR, vZ)
        d["sin_ang"] = np.where(d["vE_mag"] > 0, d["v_psi"] / d["vE_mag"], np.nan)
        d["Gamma_psi"] = d["ne"] * d["v_psi"]
        d["Gamma_par_pol"] = np.where(Babs > 0, d["ne"] * d["vpar"] * Bth / Babs, np.nan)
        d["flux_ratio"] = np.where(np.abs(d["Gamma_par_pol"]) > 0,
                                   d["Gamma_psi"] / d["Gamma_par_pol"], np.nan)

    # Sign check. grad(psi) = (-R*BZ, R*BR) is asserted, not measured, so verify it against
    # the numerically differentiated Psi_N. Only the SIGN is taken from the numerical
    # gradient, which is robust even on a coarse grid.
    R, Z, pn = d["R"], d["Z"], d["Psi_N"]
    try:
        # `rectangle` lays out a uniform grid, so differentiate on the INDEX grid and
        # rescale by the constant spacing. Passing the R/Z columns to np.gradient as
        # coordinate arrays fails outright when they carry NaNs outside the plasma --
        # which they do -- and the check then silently never runs.
        rr = R[np.isfinite(R)]
        zz = Z[np.isfinite(Z)]
        nR, nZ = R.shape[1], R.shape[0]
        dR = (np.max(rr) - np.min(rr)) / max(nR - 1, 1)
        dZ = (np.max(zz) - np.min(zz)) / max(nZ - 1, 1)
        gz, gr = np.gradient(pn)                         # axis 0 = Z, axis 1 = R
        gr, gz = gr / dR, gz / dZ
        pr_a, pz_a = -R * BZ, R * BR                     # asserted grad psi
        m = np.isfinite(gr) & np.isfinite(gz) & np.isfinite(pr_a) & np.isfinite(pz_a)
        m &= (np.hypot(gr, gz) > 0) & (np.hypot(pr_a, pz_a) > 0)
        if m.sum() > 100:
            cos = ((gr[m] * pr_a[m] + gz[m] * pz_a[m])
                   / (np.hypot(gr[m], gz[m]) * np.hypot(pr_a[m], pz_a[m])))
            d["_sign_cos"] = float(np.nanmedian(cos))
            d["_sign_n"] = int(m.sum())
        else:
            d["_sign_cos"] = float("nan")
    except Exception as e:
        d["_sign_cos"] = float("nan")
        d["_sign_err"] = str(e)
    return d


def _ang(sin_val):
    """sin_ang -> degrees. It is a SINE, so this is arcsin, not a radian conversion:
    sin_ang = 1 is 90 degrees (flow purely across flux surfaces), not 57.3."""
    if sin_val is None or not np.isfinite(sin_val):
        return float("nan")
    return math.degrees(math.asin(min(1.0, max(-1.0, sin_val))))


def box_mask(d, pfr):
    Rmin, Rmax, Zmin, Zmax = pfr
    return ((d["R"] >= Rmin) & (d["R"] <= Rmax)
            & (d["Z"] >= Zmin) & (d["Z"] <= Zmax)
            & np.isfinite(d["v_psi"]))


def stat(a, m):
    v = a[m]
    v = v[np.isfinite(v)]
    if v.size == 0:
        return dict(n=0, mean=float("nan"), absmean=float("nan"),
                    p95=float("nan"), mx=float("nan"))
    return dict(n=v.size, mean=float(np.mean(v)), absmean=float(np.mean(np.abs(v))),
                p95=float(np.percentile(np.abs(v), 95)), mx=float(np.max(np.abs(v))))


def fmt(x, p=4):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/a"
    a = abs(x)
    if a >= 1e4 or (a < 1e-3 and a > 0):
        return "%.*e" % (p - 1, x)
    return "%.*f" % (p, x)


def report_one(d, pfr, label, out):
    m = box_mask(d, pfr)
    allm = np.isfinite(d["v_psi"])

    out.append("")
    out.append("  %s   t_now = %s,  time = %s s,  step %d,  grid %s"
               % (label, fmt(d["_t_now"]), fmt(d["_time"]), d["_index"], d["_shape"]))
    sc = d.get("_sign_cos", float("nan"))
    if not np.isfinite(sc):
        # A check that silently does not run is worse than no check at all.
        out.append("    SIGN CHECK DID NOT RUN (grad Psi_N could not be differentiated --")
        out.append("    NaNs in the R/Z columns outside the plasma?). Magnitudes below are")
        out.append("    unaffected; the SIGN of v_psi is unverified.")
    else:
        if sc > 0.9:
            note = "grad psi orientation confirmed against Psi_N (+v_psi = outward)"
        elif sc < -0.9:
            note = "grad psi is ANTI-parallel to grad Psi_N: +v_psi points INWARD"
        else:
            note = ("WEAK sign check (cos = %s) -- do not trust the SIGN of v_psi, only "
                    "its magnitude" % fmt(sc, 3))
        out.append("    %s" % note)

    if m.sum() == 0:
        out.append("    WARNING: the PFR box contains no valid points. Check --pfr against")
        out.append("             the actual strike points and X-point of this equilibrium.")
        return

    out.append("    PFR box R %.3f..%.3f  Z %.3f..%.3f   (%d of %d points)"
               % (pfr[0], pfr[1], pfr[2], pfr[3], m.sum(), allm.sum()))
    out.append("")
    out.append("      %-22s %12s %12s %12s"
               % ("quantity (PFR box)", "mean|.|", "95th pct", "max|.|"))
    out.append("      " + "-" * 62)
    for key, name, unit in (("v_psi", "v_psi", "m/s"),
                            ("vE_mag", "|v_ExB|", "m/s"),
                            ("sin_ang", "sin_ang", ""),
                            ("Gamma_psi", "Gamma_psi", "m^-2 s^-1"),
                            ("Gamma_par_pol", "Gamma_par,pol", "m^-2 s^-1"),
                            ("flux_ratio", "Gamma_psi/Gamma_par", "")):
        s = stat(d[key], m)
        out.append("      %-22s %12s %12s %12s"
                   % (name + (" [%s]" % unit if unit else ""),
                      fmt(s["absmean"]), fmt(s["p95"]), fmt(s["mx"])))
    out.append("      " + "-" * 62)

    sa = stat(d["sin_ang"], m)
    fr = stat(d["flux_ratio"], m)

    # RATIO OF MEANS, not the mean of the pointwise ratio. Gamma_par,pol carries a factor
    # Btheta, which vanishes at the X-point, and vpar stagnates elsewhere in the PFR -- so
    # the pointwise ratio is singular on a set of measure zero that the grid samples
    # anyway, and its mean is dominated by those points rather than by the transport. The
    # flux-weighted ratio below is the physical one; fr["absmean"] is kept only as a
    # diagnostic of how bad the pointwise version is.
    gp = stat(d["Gamma_psi"], m)
    gq = stat(d["Gamma_par_pol"], m)
    flux_w = gp["absmean"] / gq["absmean"] if gq["absmean"] else float("nan")

    ve = stat(d["vE_mag"], m)
    vp = stat(d["v_psi"], m)
    sin_w = vp["absmean"] / ve["absmean"] if ve["absmean"] else float("nan")

    # Circulation or transport? A closed ExB cell has a large |Gamma_psi| everywhere and
    # zero NET flux. Comparing the signed mean with the mean magnitude separates the two.
    gsig = d["Gamma_psi"][m]
    gsig = gsig[np.isfinite(gsig)]
    net = float(np.mean(gsig)) if gsig.size else float("nan")
    netfrac = abs(net) / gp["absmean"] if gp["absmean"] else float("nan")
    out.append("")
    out.append("      flux-weighted   <|v_psi|>/<|v_ExB|>       = %s  (%.1f deg)"
               % (fmt(sin_w, 3), _ang(sin_w)))
    out.append("      flux-weighted   <|G_psi|>/<|G_par,pol|>   = %s" % fmt(flux_w, 3))
    out.append("      net vs circulating  |<G_psi>|/<|G_psi>|   = %s" % fmt(netfrac, 3))
    out.append("")
    out.append("      IS PHI A FLUX FUNCTION IN THE PFR?")
    if not np.isfinite(sa["absmean"]):
        out.append("        cannot tell -- no valid points.")
    elif sa["absmean"] < 0.02:
        out.append("        Essentially yes: mean |sin_ang| = %s, i.e. the ExB is within"
                   % fmt(sa["absmean"], 3))
        out.append("        %.2f degrees of the flux surfaces." % _ang(sa["absmean"]))
    else:
        out.append("        NO: mean |sin_ang| = %s, i.e. the ExB makes a mean angle of"
                   % fmt(sa["absmean"], 3))
        out.append("        %.1f degrees to the flux surfaces (peak %.1f deg)."
                   % (_ang(sa["absmean"]), _ang(sa["mx"])))
    out.append("")
    out.append("      DOES IT MOVE DENSITY?")
    if not np.isfinite(flux_w):
        out.append("        cannot tell -- parallel flux is zero or missing.")
    else:
        out.append("        <|Gamma_psi|> / <|Gamma_par,pol|> = %s" % fmt(flux_w, 3))
        if flux_w < 0.01:
            out.append("        Below 1%: the cross-field ExB is negligible against parallel")
            out.append("        streaming, and this structure cannot drive an asymmetry.")
        elif flux_w < 0.1:
            out.append("        A few per cent -- real but subdominant.")
        else:
            out.append("        Above 10%: the cross-field ExB is a leading-order transport")
            out.append("        channel in the PFR.")
        if np.isfinite(netfrac):
            out.append("")
            if netfrac < 0.1:
                out.append("        BUT |<Gamma_psi>|/<|Gamma_psi|> = %s: the flux very nearly"
                           % fmt(netfrac, 3))
                out.append("        CANCELS. This is a closed circulation cell, not a net")
                out.append("        transfer between the legs. A large circulating flux moves")
                out.append("        no inventory. For the NET inter-leg transfer use the signed")
                out.append("        integral through a dividing surface -- postproc_pfr_line.in")
                out.append("        with analyse_pfr_line.py, whose cut is vertical at the")
                out.append("        X-point major radius precisely for this reason.")
            else:
                out.append("        |<Gamma_psi>|/<|Gamma_psi|> = %s, so a real fraction of this"
                           % fmt(netfrac, 3))
                out.append("        is NET, not circulating. Confirm the direction and the")
                out.append("        magnitude with analyse_pfr_line.py.")
        if fr["absmean"] > 3 * flux_w:
            out.append("")
            out.append("        (The pointwise mean ratio is %s -- inflated by points where"
                       % fmt(fr["absmean"], 3))
            out.append("        Btheta or vpar approach zero. Ignore it; use the flux-weighted")
            out.append("        number above.)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", metavar="FILE.h5", help="one file, or REF and TEST")
    ap.add_argument("--pfr", nargs=4, type=float, default=list(DEFAULT_PFR),
                    metavar=("Rmin", "Rmax", "Zmin", "Zmax"))
    ap.add_argument("--labels", nargs=2, default=None)
    ap.add_argument("--vars", default=None, help="expression order, space separated")
    ap.add_argument("--plot", default=None, metavar="FILE.png")
    args = ap.parse_args()

    if len(args.files) > 2:
        sys.exit("ERROR: give one file, or two for an A/B.")
    varnames = args.vars.split() if args.vars else list(DEFAULT_VARS)
    for need in ("R", "Z", "V_ExB_R", "V_ExB_Z", "BR", "BZ", "Btheta"):
        if need not in varnames:
            sys.exit("ERROR: '%s' must be in the expression list." % need)

    labels = args.labels or [os.path.basename(p) for p in args.files]
    ds = [derive(load(p, varnames)) for p in args.files]

    out = ["=" * 78, "FLUX-FUNCTION TEST   v_psi = (v_ExB . grad psi)/|grad psi|", "=" * 78,
           "  Zero everywhere iff Phi = Phi(psi). Derivatives are analytic (Bezier basis),",
           "  not differenced -- grad psi comes from the exported B components."]

    for d, lab in zip(ds, labels):
        report_one(d, tuple(args.pfr), lab, out)

    if len(ds) == 2:
        a, b = ds
        if a["_shape"] != b["_shape"]:
            out += ["", "ERROR: the two exports are on different grids %s vs %s -- rerun both"
                    % (a["_shape"], b["_shape"]), "       with the same `rectangle` line."]
        else:
            dt = abs(a["_time"] - b["_time"])
            rel = dt / max(abs(a["_time"]), abs(b["_time"]), 1e-30)
            out += ["", "=" * 78, "A/B   %s  minus  %s" % (labels[1], labels[0]), "=" * 78]
            if np.isfinite(rel) and rel > 1e-6:
                out += ["  WARNING: the two frames are at DIFFERENT times (%s vs %s s)."
                        % (fmt(a["_time"]), fmt(b["_time"])),
                        "           A difference between different times is not an A/B."]
            m = box_mask(a, tuple(args.pfr)) & box_mask(b, tuple(args.pfr))
            out.append("")
            out.append("  %-22s %12s %12s %12s %10s"
                       % ("quantity (PFR box)", labels[0][:12], labels[1][:12],
                          "delta", "rel"))
            out.append("  " + "-" * 74)
            for key, name in (("v_psi", "mean |v_psi| [m/s]"),
                              ("sin_ang", "mean |sin_ang|"),
                              ("flux_ratio", "mean |Gamma ratio|"),
                              ("Phi", "mean Phi [V]"),
                              ("ne", "mean ne [m^-3]")):
                if key not in a:
                    continue
                use = np.abs if key in ("v_psi", "sin_ang", "flux_ratio") else (lambda x: x)
                va, vb = use(a[key])[m], use(b[key])[m]
                va, vb = va[np.isfinite(va)], vb[np.isfinite(vb)]
                if va.size == 0 or vb.size == 0:
                    continue
                ma, mb = float(np.mean(va)), float(np.mean(vb))
                dd = mb - ma
                r = 100.0 * dd / abs(ma) if ma else float("nan")
                out.append("  %-22s %12s %12s %12s %9s%%"
                           % (name, fmt(ma), fmt(mb), fmt(dd),
                              fmt(r, 3) if np.isfinite(r) else "n/a"))
            out.append("  " + "-" * 74)

            # Where does the difference live? A change concentrated in the PFR means the
            # transport channel between the legs was modified; one spread flat over the
            # whole domain is a global offset with no structural consequence.
            dv = b["v_psi"] - a["v_psi"]
            fin = np.isfinite(dv)
            if fin.sum() and m.sum():
                inb = float(np.mean(np.abs(dv[m])))
                allb = float(np.mean(np.abs(dv[fin])))
                conc = inb / allb if allb else float("nan")
                out += ["",
                        "  mean |delta v_psi|:  %s m/s in the PFR box,  %s over the whole"
                        % (fmt(inb), fmt(allb)),
                        "  export.  Concentration ratio %s." % fmt(conc, 3)]
                if not np.isfinite(conc):
                    pass
                elif conc > 1.5:
                    out += ["  Above 1.5: the change is CONCENTRATED in the PFR, i.e. it",
                            "  modifies the inter-leg transport channel."]
                elif conc > 0.7:
                    out += ["  Near 1: a roughly uniform offset, with no particular",
                            "  structural preference for the PFR."]
                else:
                    out += ["  Below 0.7: the change AVOIDS the PFR -- it is %.1fx weaker"
                            % (1.0 / conc if conc else float("nan"),),
                            "  there than elsewhere in the export. Whatever moved, it was",
                            "  not the inter-leg channel; look at the legs and the targets."]

    out.append("")
    out.append("=" * 78)
    print("\n".join(out))

    if args.plot:
        make_plot(args.plot, ds, labels, tuple(args.pfr))
    return 0


def make_plot(path, ds, labels, pfr):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not available -- skipping --plot)")
        return
    ncol = len(ds) + (1 if len(ds) == 2 else 0)
    fig, axes = plt.subplots(2, ncol, figsize=(5.2 * ncol, 8.4), squeeze=False)

    def show(ax, d, key, title, lim=None, cmap="RdBu_r"):
        R, Z, V = d["R"], d["Z"], d[key]
        v = V[np.isfinite(V)]
        if lim is None:
            lim = float(np.percentile(np.abs(v), 99)) if v.size else 1.0
        im = ax.pcolormesh(R, Z, V, cmap=cmap, vmin=-lim, vmax=lim, shading="auto")
        ax.plot([pfr[0], pfr[1], pfr[1], pfr[0], pfr[0]],
                [pfr[2], pfr[2], pfr[3], pfr[3], pfr[2]], "k-", lw=1.0)
        if "Psi_N" in d:
            ax.contour(R, Z, d["Psi_N"], levels=[1.0], colors="k", linewidths=1.2)
        ax.set_title(title, fontsize=9)
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, fraction=0.046)
        return lim

    # Common limits across cases, so a 2% difference cannot hide in the colour scaling.
    lim_v = max(float(np.percentile(np.abs(d["v_psi"][np.isfinite(d["v_psi"])]), 99))
                for d in ds)
    for j, (d, lab) in enumerate(zip(ds, labels)):
        show(axes[0][j], d, "v_psi", "v_psi  [m/s]   %s" % lab, lim=lim_v)
        show(axes[1][j], d, "sin_ang", "sin_ang   %s" % lab, lim=1.0)

    if len(ds) == 2:
        a, b = ds
        diff = {k: b[k] - a[k] for k in ("v_psi", "sin_ang")}
        diff.update({"R": a["R"], "Z": a["Z"], "Psi_N": a["Psi_N"]})
        show(axes[0][2], diff, "v_psi", "delta v_psi  [m/s]")
        show(axes[1][2], diff, "sin_ang", "delta sin_ang")

    fig.suptitle("Flux-function test: v_psi = 0 everywhere iff Phi = Phi(psi)"
                 "   (black box = PFR, black line = separatrix)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=130)
    print("\nwrote %s" % path)


if __name__ == "__main__":
    sys.exit(main())
