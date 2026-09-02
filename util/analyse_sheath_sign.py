#!/usr/bin/env python3
"""Decide the sheath a_n sign question from a boundary_quantities export.

    jorek2_postproc < util/postproc_sheath_sign.in
    util/analyse_sheath_sign.py boundary_quantities_*.dat

THE TEST (doc/sheath_sign_review.md section 9). At a strike point in ion saturation the
conventional ion current and the ion parallel velocity point the same way along B, so with
F0 > 0, sign(Jpar) == sign(Jpar_ionsat). Jpar comes from Jtor = -zj0/R; Jpar_ionsat is
rho*vpar*Btot. Neither passes through mod_sheath_bc, so this is independent of the
coefficient under suspicion.

DISAGREE  -> a_n's sign is wrong; revert it (c_sat follows automatically) and add the
             missing minus in mod_sheath_diag.
AGREE     -> the proposed fix is WRONG. Do not apply it.

Sampling filters, all of which matter - a point that fails any of them cannot settle the
question and is excluded rather than counted as a disagreement.
"""
import sys, re, argparse, pathlib
import numpy as np

# --- copied verbatim from util/analyse_targets.py, which is known to parse real
# --- jorek2_postproc output. The marker line is "# time step #  N ... t_now = X".
BLOCK_RE = re.compile(r"time step #\s*(\d+).*?t_now\s*=\s*([-\dEe.+]+)")


def _isnum(t):
    try:
        float(t); return True
    except ValueError:
        return False


def read_boundary_file(path):
    """Same format as util/analyse_targets.py: '# name name ...' then rows."""
    names, blocks, rows, step, t_now = None, [], [], None, None

    def flush():
        if step is not None and rows:
            a = np.asarray(rows, float)
            blocks.append((step, t_now, {n: a[:, i] for i, n in enumerate(names)}))

    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                m = BLOCK_RE.search(s)
                if m:
                    flush(); rows = []
                    step, t_now = int(m.group(1)), float(m.group(2))
                else:
                    cand = s.lstrip("#").split()
                    if cand and not _isnum(cand[0]):
                        names = cand
                continue
            rows.append([float(v) for v in s.split()])
    flush()
    if names is None:
        sys.exit(f"{path}: no header line (need header=.true.)")
    if not blocks and rows:
        # --- No "time step #" marker anywhere, but the file does have data: it is a
        # --- single-step export. Take the step from the filename and use it as one block.
        m = re.search(r"s0*(\d+)", pathlib.Path(path).name)
        a = np.asarray(rows, float)
        blocks = [(int(m.group(1)) if m else 0, 0.0,
                   {n: a[:, i] for i, n in enumerate(names)})]
    if not blocks:
        sys.exit(f"{path}: no data blocks. First lines of the file were:\n  " +
                 "\n  ".join(open(path).read().splitlines()[:8]))
    return blocks


def analyse(d, a):
    need = ["Jpar", "Jpar_ionsat", "Bnorm", "B_abs", "ne", "T_e", "R", "Z"]
    miss = [n for n in need if n not in d]
    if miss:
        sys.exit("missing expressions in the export: " + " ".join(miss)
                 + "\n(add them to the `expressions` line in postproc_sheath_sign.in)")

    jp, js = d["Jpar"], d["Jpar_ionsat"]
    bn_frac = np.abs(d["Bnorm"]) / np.maximum(d["B_abs"], 1e-300)
    rx = np.hypot(d["R"] - a.xpt_r, d["Z"] - a.xpt_z)

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.abs(jp) / np.abs(js)

    cuts = [
        ("finite",              np.isfinite(jp) & np.isfinite(js)),
        (f"|B.n|/B > {a.bnmin}", bn_frac > a.bnmin),
        (f"ne > {a.nemin:g}",    d["ne"] > a.nemin),
        (f"T_e > {a.temin} eV",  d["T_e"] > a.temin),
        (f"> {a.xrad} m from X-point", rx > a.xrad),
        (f"{a.rmin} < |Jpar/Jsat| < {a.rmax} (ion branch)",
         np.isfinite(ratio) & (ratio > a.rmin) & (ratio < a.rmax)),
    ]

    print(f"  {len(jp)} boundary points in this block")
    keep = np.ones_like(jp, bool)
    for label, c in cuts:
        keep &= c
        if keep.sum():
            Rs, Zs = d["R"][keep], d["Z"][keep]
            print(f"    after cut  {label:<38s} : {keep.sum():5d}"
                  f"   R {Rs.min():.3f}..{Rs.max():.3f}  Z {Zs.min():+.3f}..{Zs.max():+.3f}")
        else:
            print(f"    after cut  {label:<38s} : {keep.sum():5d}")

    n = int(keep.sum())
    if n < a.minpts:
        print(f"\n  ONLY {n} POINTS SURVIVE (need >= {a.minpts}). Inconclusive - relax a cut,")
        print("  or export a step where the targets are more clearly in ion saturation.")
        return None

    agree = np.sign(jp[keep]) == np.sign(js[keep])
    frac = agree.mean()
    print(f"\n  sign(Jpar) == sign(Jpar_ionsat) on {agree.sum()}/{n} points  ({100*frac:.1f} %)")
    print(f"  median |Jpar|/|Jpar_ionsat| = {np.median(ratio[keep]):.3f}"
          "   (near 1 => the two are the same physical quantity)")

    # --- Split by target. A CONVENTION error is global: both targets disagree at
    # --- about the same rate. A target-asymmetric result means something else
    # --- (a local j_sat collapse, one leg detached) and does NOT implicate a_n.
    Rk = d["R"][keep]
    for label, m in (("inner (R < %.2f)" % a.rsplit, Rk < a.rsplit),
                     ("outer (R > %.2f)" % a.rsplit, Rk >= a.rsplit)):
        if m.sum():
            print(f"    {label:<20s}: {agree[m].sum():4d}/{m.sum():<4d} agree "
                  f"({100*agree[m].mean():5.1f} %)")
        else:
            print(f"    {label:<20s}: no qualifying points")

    edges = np.linspace(Rk.min(), Rk.max(), 9)
    print("    agreement vs R (look for TWO separated clusters = two targets):")
    for b in range(len(edges) - 1):
        m = (Rk >= edges[b]) & (Rk <= edges[b + 1] if b == len(edges) - 2 else Rk < edges[b + 1])
        bar = "#" * int(round(20 * m.sum() / max(1, len(Rk))))
        if m.sum():
            print(f"      R {edges[b]:.3f}-{edges[b+1]:.3f} : {agree[m].sum():4d}/{m.sum():<4d}"
                  f" ({100*agree[m].mean():5.1f} %) {bar}")
        else:
            print(f"      R {edges[b]:.3f}-{edges[b+1]:.3f} :    -        (empty)")

    if "Phi" in d:
        print(f"  Phi over the same points: median {np.median(d['Phi'][keep]):+.2f} V"
              f"   (expect ~ +Lambda*Te, i.e. POSITIVE, after the fix)")
    return frac


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+")
    p.add_argument("--bnmin", type=float, default=0.02,
                   help="minimum |B.n|/|B| - excludes grazing incidence (default 0.02, "
                        "4x sheath_min_bn)")
    p.add_argument("--nemin", type=float, default=1e18, help="minimum ne [m^-3]")
    p.add_argument("--temin", type=float, default=2.0,
                   help="minimum T_e [eV] - stays clear of the t_min floor at ~0.49 eV")
    p.add_argument("--xrad",  type=float, default=0.05, help="exclusion radius around the X-point [m]")
    p.add_argument("--xpt-r", type=float, default=1.43)
    p.add_argument("--xpt-z", type=float, default=-0.93)
    p.add_argument("--rmin",  type=float, default=0.3, help="lower bound on |Jpar/Jpar_ionsat|")
    p.add_argument("--rmax",  type=float, default=3.0, help="upper bound on |Jpar/Jpar_ionsat|")
    p.add_argument("--minpts", type=int,  default=20)
    p.add_argument("--rsplit", type=float, default=1.42,
                   help="R dividing inner from outer target [m], as sheath_diag_R_split")
    a = p.parse_args()

    results = []
    for f in a.files:
        for step, t_now, d in read_boundary_file(f):
            print(f"\n=== {f}   step {step}   t_now {t_now:.6g}")
            r = analyse(d, a)
            if r is not None:
                results.append(r)

    if not results:
        sys.exit("\nno block produced enough qualifying points - inconclusive")

    m = float(np.mean(results))
    print("\n" + "=" * 72)
    if m < 0.25:
        print(f"VERDICT: the two DISAGREE ({100*(1-m):.0f} % of qualifying points).")
        print("  => a_n's sign is wrong. Revert it (c_sat = -0.5*a_n follows automatically)")
        print("     AND add the missing minus in mod_sheath_diag.f90:159-160. Apply together.")
    elif m > 0.75:
        print("VERDICT: the two AGREE almost everywhere.")
        print("  => the proposed fix is WRONG. Do NOT apply it. The +zj convention in")
        print("     mod_sheath_diag.f90:160 is right and the derivation needs revisiting.")
    else:
        print(f"VERDICT: MIXED ({100*m:.0f} % agree). Not a clean answer - the sample is")
        print("  probably contaminated. Tighten --bnmin, --temin and --rmin/--rmax and re-run.")
    print("=" * 72)


if __name__ == "__main__":
    main()
