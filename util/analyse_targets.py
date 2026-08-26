#!/usr/bin/env python3
"""
Divertor target profile analysis for the model600 sheath boundary condition.

Reads a `boundary_quantities_*.dat` written by jorek2_postproc (see
util/postproc_targets.in) and answers the question the SHEATH log cannot:

  - What is Phi in VOLTS at each target, and what is the potential difference
    across the private flux region?  The log reports ePhi/kTe, which hides the
    Te contrast: a 6 % difference in ePhi/kTe was a factor ~4 in volts.
  - What is the tangential electric field along each target, i.e. the E that
    drives the ExB drift through the PFR?
  - Do n_e and T_e separate the way HFSHD requires?  j_sat ~ n*sqrt(T), so the
    currents alone cannot distinguish "inner detaching" (n up, T down) from
    "inner starved" (n down, T flat).

Two sign/naming traps are handled here; see postproc_targets.in for the detail:
  * the `Phi` expression is +F0*u, but model600 has Phi = -F0*u  -> negated
  * `Te` is the mean temperature; the real ones are `T_e` and `T_i`

Usage
-----
  ./analyse_targets.py run.dat
  ./analyse_targets.py unpuffed.dat puffed.dat --labels baseline 10x-puff
  ./analyse_targets.py run.dat --step 3900 --rsplit 1.42 --zmax -0.90
"""

import argparse
import re
import sys

import numpy as np

# --- Physical constants (SI)
E_CHG = 1.602176634e-19
M_U = 1.66053906660e-27
GAMMA = 5.0 / 3.0

BLOCK_RE = re.compile(r"time step #\s*(\d+).*?t_now\s*=\s*([-\dEe.+]+)")


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------
def read_boundary_file(path):
    """Parse a boundary_quantities .dat into a list of (step, t, {name: array}).

    Blocks are delimited by the '# time step #NNNNNN,  t_now = ...' comment
    that write_ascii_1d emits, so blank-line conventions do not matter.
    """
    names, blocks = None, []
    step = t_now = None
    rows = []

    def flush():
        if step is not None and rows:
            arr = np.asarray(rows, dtype=float)
            blocks.append((step, t_now, {n: arr[:, i] for i, n in enumerate(names)}))

    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                m = BLOCK_RE.search(s)
                if m:
                    flush()
                    rows = []
                    step, t_now = int(m.group(1)), float(m.group(2))
                else:
                    # header line: expression names, no spaces inside a name
                    cand = s.lstrip("#").split()
                    if cand and not _is_number(cand[0]):
                        names = cand
                continue
            rows.append([float(v) for v in s.split()])
    flush()

    if names is None:
        sys.exit(f"{path}: no header line found - was `header=.true.` written?")
    if not blocks:
        sys.exit(f"{path}: no data blocks found.")
    return blocks


def _is_number(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def pick_block(blocks, step):
    if step is None:
        return blocks[-1]
    for b in blocks:
        if b[0] == step:
            return b
    avail = ", ".join(str(b[0]) for b in blocks)
    sys.exit(f"step {step} not in file. Available: {avail}")


# ---------------------------------------------------------------------------
# Target extraction
# ---------------------------------------------------------------------------
class Target:
    """One divertor target, ordered by arc length from its strike point."""

    def __init__(self, name, d, mask, mass_amu):
        self.name = name
        self.n = int(mask.sum())
        if self.n < 3:
            self.valid = False
            return
        self.valid = True

        g = {k: v[mask] for k, v in d.items()}

        # Order along the surface: monotonic in Z on a near-vertical plate,
        # in R on a near-horizontal one. Pick whichever spans further.
        R, Z = g["R"], g["Z"]
        key = R if (R.max() - R.min()) >= (Z.max() - Z.min()) else Z
        o = np.argsort(key)
        g = {k: v[o] for k, v in g.items()}

        self.R, self.Z = g["R"], g["Z"]
        self.psin = g["Psi_N"]
        self.bnorm = g.get("Bnorm", np.zeros_like(self.R))

        # SIGN FIX: the postproc `Phi` expression returns +F0*u; model600 has
        # Phi = -F0*u. Without this the ExB direction is inverted.
        self.phi = -g["Phi"]

        self.te = g["T_e"]
        self.ti = g.get("T_i", g["T_e"])
        self.ne = g["ne"]
        self.vr = g.get("V_ExB_R", np.zeros_like(self.R))
        self.vz = g.get("V_ExB_Z", np.zeros_like(self.R))

        # Arc length along the target, zeroed at the strike point (Psi_N = 1).
        ds = np.hypot(np.diff(self.R), np.diff(self.Z))
        s = np.concatenate([[0.0], np.cumsum(ds)])
        self.s = s - self._strike_value(s)

        # Tangential electric field, E_t = -dPhi/ds. This is the component that
        # drives an ExB flow ACROSS field lines, i.e. through the PFR.
        self.et = -np.gradient(self.phi, self.s)

        # Ion saturation current density, j_sat = e * n * c_s * |b.n|.
        cs = np.sqrt(GAMMA * (self.ti + self.te) * E_CHG / (mass_amu * M_U))
        self.jsat = E_CHG * self.ne * cs * np.abs(self.bnorm) if np.any(self.bnorm) \
            else E_CHG * self.ne * cs
        self.vexb = np.hypot(self.vr, self.vz)

    def _strike_value(self, arr):
        """Value of `arr` where Psi_N crosses 1, by linear interpolation."""
        p = self.psin
        i = np.argmin(np.abs(p - 1.0))
        for j in range(len(p) - 1):
            if (p[j] - 1.0) * (p[j + 1] - 1.0) <= 0 and p[j] != p[j + 1]:
                w = (1.0 - p[j]) / (p[j + 1] - p[j])
                return arr[j] + w * (arr[j + 1] - arr[j])
        return arr[i]

    def at_strike(self, arr):
        return self._strike_value(arr)

    def peak(self, arr):
        i = int(np.argmax(np.abs(arr)))
        return arr[i], self.s[i], self.psin[i]


def split_targets(d, rsplit, zmax, mass_amu):
    """Divertor points only (Z < zmax), split inner/outer on major radius."""
    div = d["Z"] < zmax
    if not div.any():
        sys.exit(f"no boundary points with Z < {zmax}. Adjust --zmax.")
    inner = div & (d["R"] < rsplit)
    outer = div & (d["R"] >= rsplit)
    return (Target("INNER", d, inner, mass_amu),
            Target("OUTER", d, outer, mass_amu))


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def report(label, step, t_now, inner, outer):
    print(f"\n{'=' * 74}\n{label}   step {step}   t = {t_now:.6g} s\n{'=' * 74}")

    if not (inner.valid and outer.valid):
        print("  one target has too few points - check --rsplit / --zmax")
        return None

    print(f"{'':22s}{'INNER':>16s}{'OUTER':>16s}{'in/out':>12s}")

    def row(name, a, b, fmt="{:16.4g}", ratio=True):
        r = f"{a / b:12.3f}" if (ratio and b != 0) else " " * 12
        print(f"  {name:20s}{fmt.format(a)}{fmt.format(b)}{r}")

    pi, po = inner.at_strike(inner.phi), outer.at_strike(outer.phi)
    tei, teo = inner.at_strike(inner.te), outer.at_strike(outer.te)
    nei, neo = inner.at_strike(inner.ne), outer.at_strike(outer.ne)

    print("\n  -- at the strike point (Psi_N = 1) --")
    row("Phi [V]", pi, po)
    row("T_e [eV]", tei, teo)
    row("T_i [eV]", inner.at_strike(inner.ti), outer.at_strike(outer.ti))
    row("n_e [1e20 m^-3]", nei / 1e20, neo / 1e20)
    row("j_sat [kA/m^2]", inner.at_strike(inner.jsat) / 1e3,
        outer.at_strike(outer.jsat) / 1e3)

    print("\n  -- peak over the target --")
    # E_t is reported SIGNED at the point of largest magnitude: the sign is the
    # direction of the drive, which is the whole question for HFSHD.
    for nm, attr, sc, unit in [("n_e", "ne", 1e20, "1e20 m^-3"),
                               ("T_e", "te", 1.0, "eV"),
                               ("E_t at max|E_t|", "et", 1.0, "V/m"),
                               ("|v_ExB|", "vexb", 1.0, "m/s")]:
        a, sa, _ = inner.peak(getattr(inner, attr))
        b, sb, _ = outer.peak(getattr(outer, attr))
        row(f"{nm} [{unit}]", a / sc, b / sc)
        print(f"  {'':20s}{'at s=' + format(sa, '.3f') + ' m':>16s}"
              f"{'at s=' + format(sb, '.3f') + ' m':>16s}")

    dphi = pi - po
    print(f"\n  ** PFR potential difference  Phi_in - Phi_out = {dphi:+.2f} V **")
    print(f"     ePhi/kTe would report {pi / max(tei, 1e-12):.2f} vs "
          f"{po / max(teo, 1e-12):.2f} - a {abs(pi / max(tei,1e-12) - po / max(teo,1e-12)) / max(pi / max(tei,1e-12), 1e-12) * 100:.0f} % contrast "
          f"for a factor {abs(pi / po) if po else float('nan'):.2f} in volts.")
    print(f"     n_e(in)/n_e(out) at the strike point = {nei / neo:.3f}"
          "   (> 1 is the HFSHD direction)")
    return dphi


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def plot(cases, outfile):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib not available - text report only "
              "(pip install matplotlib for the figure).")
        return

    panels = [("phi", r"$\Phi$  [V]", 1.0),
              ("te", r"$T_e$  [eV]", 1.0),
              ("ne", r"$n_e$  [$10^{20}$ m$^{-3}$]", 1e20),
              ("jsat", r"$j_{sat}$  [kA m$^{-2}$]", 1e3),
              ("et", r"$E_t = -\partial\Phi/\partial s$  [V/m]", 1.0),
              ("vexb", r"$|v_{E\times B}|$  [m/s]", 1.0)]

    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]

    for ax, (attr, ylab, sc) in zip(axes.ravel(), panels):
        for ci, (label, _, _, inner, outer) in enumerate(cases):
            c = colors[ci % len(colors)]
            for tgt, ls in ((inner, "-"), (outer, "--")):
                if tgt.valid:
                    ax.plot(tgt.psin, getattr(tgt, attr) / sc, ls, color=c, lw=1.6,
                            label=f"{label} {tgt.name}")
        ax.set_ylabel(ylab)
        ax.axvline(1.0, color="k", lw=0.8, ls=":")
        ax.grid(alpha=0.3)

    for ax in axes[1]:
        ax.set_xlabel(r"$\psi_N$      (< 1 = PFR,  > 1 = SOL)")
    axes[0, 0].legend(fontsize=7, ncol=1)
    fig.suptitle("Divertor target profiles  -  solid = inner, dashed = outer", y=0.99)
    fig.tight_layout()
    fig.savefig(outfile, dpi=140)
    print(f"\nwrote {outfile}")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="boundary_quantities_*.dat")
    ap.add_argument("--labels", nargs="*", default=None)
    ap.add_argument("--step", type=int, default=None,
                    help="time step to analyse (default: last in file)")
    ap.add_argument("--rsplit", type=float, default=1.42,
                    help="major radius dividing inner/outer, matching "
                         "sheath_diag_R_split (default 1.42)")
    ap.add_argument("--zmax", type=float, default=-0.90,
                    help="only boundary points below this Z are divertor "
                         "targets (default -0.90)")
    ap.add_argument("--mass", type=float, default=2.0141,
                    help="main ion mass in amu (default deuterium)")
    ap.add_argument("-o", "--out", default="target_profiles.png")
    ap.add_argument("--no-plot", action="store_true", help="text report only")
    args = ap.parse_args()

    labels = args.labels or [f.split("/")[-1] for f in args.files]
    if len(labels) != len(args.files):
        sys.exit("--labels must have one entry per file")

    cases, dphis = [], []
    for path, label in zip(args.files, labels):
        step, t_now, d = pick_block(read_boundary_file(path), args.step)
        inner, outer = split_targets(d, args.rsplit, args.zmax, args.mass)
        dphi = report(label, step, t_now, inner, outer)
        cases.append((label, step, t_now, inner, outer))
        if dphi is not None:
            dphis.append((label, dphi))

    if len(dphis) == 2:
        (l0, d0), (l1, d1) = dphis
        print(f"\n{'=' * 74}\n  PFR drive, {l0} -> {l1}: "
              f"{d0:+.2f} V -> {d1:+.2f} V   (change {d1 - d0:+.2f} V)\n{'=' * 74}")

    if not args.no_plot:
        plot(cases, args.out)


if __name__ == "__main__":
    main()
