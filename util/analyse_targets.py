#!/usr/bin/env python3
"""
Divertor target profile analysis for the model600 sheath boundary condition.

Reads `boundary_quantities_*.dat` from jorek2_postproc (see
util/postproc_targets.in) and reports what the SHEATH log cannot:

  - Phi in VOLTS at each target, and the potential difference across the
    private flux region. The log reports ePhi/kTe, which divides out the Te
    contrast that carries most of the asymmetry.
  - The tangential electric field along the target, E_t = -dPhi/ds, which is
    what drives an ExB flow across field lines through the PFR.
  - n_e and T_e separately. j_sat ~ n*sqrt(T), so currents alone cannot tell
    "inner detaching" (n up, T down) from "inner starved" (n down, T flat).

GEOMETRY NOTES
--------------
* psi_N does NOT separate SOL from PFR: it exceeds 1 on BOTH sides of the
  separatrix leg. The abscissa here is signed arc length from the strike point,
  positive toward the PFR (i.e. toward the other target) on both targets.
* Boundary points arrive in boundary-element order from bnd_pos(), so file
  order is path order. Do not sort them - the target curves back on itself and
  any sort by R or Z scrambles it.
* Z < zmax alone also selects the extended wall (type 5 is ~50 m^2 here), so
  points are additionally required to carry a significant normal field.

UNIT / SIGN TRAPS
-----------------
* Use `set units 0` in the postproc script. With JOREK units n_e comes out
  ~1e-2 and T_e ~1e-3; the script detects this and refuses rather than
  reporting nonsense.
* The `Phi` expression returns +F0*u (mod_expression.f90:1757, flagged
  "### sign?"), but model600 has Phi = -F0*u. Negated here.
* `Te` is (Ti+Te)/2, not the electron temperature. Use `T_e` and `T_i`.
"""

import argparse
import re
import sys

import numpy as np

_trapz = getattr(np, "trapezoid", None) or np.trapz

E_CHG = 1.602176634e-19
M_U = 1.66053906660e-27
GAMMA = 5.0 / 3.0

BLOCK_RE = re.compile(r"time step #\s*(\d+).*?t_now\s*=\s*([-\dEe.+]+)")


# ---------------------------------------------------------------------------
def read_boundary_file(path):
    """Parse into [(step, t_now, {name: array})], preserving file order."""
    names, blocks, rows = None, [], []
    step = t_now = None

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
                    flush()
                    rows = []
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
    if not blocks:
        sys.exit(f"{path}: no data blocks")
    return blocks


def _isnum(t):
    try:
        float(t)
        return True
    except ValueError:
        return False


def pick_block(blocks, step):
    if step is None:
        return blocks[-1]
    for b in blocks:
        if b[0] == step:
            return b
    sys.exit(f"step {step} not found. Have: {', '.join(str(b[0]) for b in blocks)}")


def check_units(d, path):
    """JOREK-unit output is off by ~1e20 in n and ~5e4 in T. Catch it early."""
    n, t = np.nanmax(d["ne"]), np.nanmax(d["T_e"])
    if n < 1e12 or t < 0.05:
        sys.exit(
            f"\n{path}: values look like JOREK units, not SI\n"
            f"    max n_e = {n:.3g}  (expect ~1e18-1e21 m^-3)\n"
            f"    max T_e = {t:.3g}  (expect ~1-100 eV)\n"
            f"  Add `set units 0` to the postproc script and re-export.\n"
            f"  (units defaults to 0=SI, so something in your session set it to 1.)\n")


# ---------------------------------------------------------------------------
def contiguous_runs(mask):
    """[(start, stop)] index ranges of True runs."""
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return []
    brk = np.flatnonzero(np.diff(idx) > 1)
    starts = np.concatenate([[idx[0]], idx[brk + 1]])
    stops = np.concatenate([idx[brk], [idx[-1]]])
    return list(zip(starts, stops + 1))


def _smooth(a, w):
    w = int(w)
    if w <= 1 or a.size <= w:
        return a
    k = np.ones(w) / float(w)
    return np.convolve(np.pad(a, w // 2, mode="edge"), k, mode="valid")[:a.size]


class Target:
    """One divertor target, in boundary-path order, s=0 at the strike point."""

    def __init__(self, name, d, sl, mass_amu, pfr_toward_larger_R,
                 half_width=0.12, smooth=9):
        self.name, self.valid, self.edge_warn = name, False, False
        g = {k: v[sl] for k, v in d.items()}
        if len(g["R"]) < 5:
            return

        # STRIKE POINT = the minimum of psi_N along the surface. psi_N exceeds 1
        # on both sides of the leg and equals 1 at the separatrix, so its minimum
        # IS the strike point. This is geometric: unlike the j_sat peak it does
        # not wander when the target detaches and j_sat goes flat.
        i0 = int(np.nanargmin(g["Psi_N"]))
        if i0 < 3 or i0 > len(g["R"]) - 4:
            self.edge_warn = True

        # Trim to +-half_width of path length around it, so the target is the
        # target and not the baffle above it. Z-cuts alone cannot do this: the
        # surface runs continuously from the strike point up the vessel wall.
        RR, ZZ = g["R"], g["Z"]
        sfull = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(RR), np.diff(ZZ)))])
        keep = np.abs(sfull - sfull[i0]) <= half_width
        g = {k: v[keep] for k, v in g.items()}
        if len(g["R"]) < 5:
            return
        self.istrike = int(np.nanargmin(g["Psi_N"]))
        self.valid = True

        self.R, self.Z, self.psin = g["R"], g["Z"], g["Psi_N"]
        self.bnorm = g.get("Bnorm", np.zeros_like(self.R))
        self.phi = -g["Phi"]          # SIGN FIX: postproc returns +F0*u
        self.te, self.ti, self.ne = g["T_e"], g.get("T_i", g["T_e"]), g["ne"]
        self.vr = g.get("V_ExB_R", np.zeros_like(self.R))
        self.vz = g.get("V_ExB_Z", np.zeros_like(self.R))
        self.vexb = np.hypot(self.vr, self.vz)

        cs = np.sqrt(GAMMA * (self.ti + self.te) * E_CHG / (mass_amu * M_U))
        self.jsat = E_CHG * self.ne * cs * np.abs(self.bnorm) / \
            max(np.nanmax(np.abs(self.bnorm)), 1e-30)

        # Arc length in PATH order (never sorted), zeroed at the strike point.
        ds = np.hypot(np.diff(self.R), np.diff(self.Z))
        s = np.concatenate([[0.0], np.cumsum(ds)])
        s = s - s[self.istrike]

        # Orient so that positive s points into the PFR on BOTH targets, making
        # the two directly comparable. psi_N cannot do this: it is > 1 on both
        # sides of the leg.
        span = self.R[-1] - self.R[0]
        self.flip = 1.0 if ((span > 0) == pfr_toward_larger_R) else -1.0
        self.s = self.flip * s

        # --- ExB direction, decomposed on the surface.
        # t_hat is the unit tangent pointing the same way as +s, i.e. TOWARD THE
        # PFR on both targets, so v_t is directly comparable between them.
        # n_hat is the outward normal as exported by the boundary expressions.
        ds_path = np.gradient(s)
        with np.errstate(invalid="ignore", divide="ignore"):
            tR = self.flip * np.gradient(self.R) / ds_path
            tZ = self.flip * np.gradient(self.Z) / ds_path
        tn = np.hypot(tR, tZ)
        tn[tn == 0] = 1.0
        self.tR, self.tZ = tR / tn, tZ / tn

        # v_ExB comes from grad(u); on a C1 bicubic that derivative is only C0, so
        # with nsub_bnd sub-sampling it carries element-scale ripple. Smooth the
        # components on the same window as Phi before taking directions.
        self.vr = _smooth(self.vr, smooth)
        self.vz = _smooth(self.vz, smooth)
        self.vt = self.vr * self.tR + self.vz * self.tZ     # + = toward the PFR
        nR = g.get("bnd_normal_R", np.zeros_like(self.R))
        nZ = g.get("bnd_normal_Z", np.zeros_like(self.R))
        self.vn = self.vr * nR + self.vz * nZ               # + = out of the plasma
        self.gam_t = self.ne * self.vt                      # particle flux [m^-2 s^-1]

        # E_t = -dPhi/ds. Guard repeated abscissa values, which would give inf.
        # Phi is smoothed first: nsub_bnd sub-sampling puts several points inside
        # one element, and differentiating the raw trace gives 1e4 V/m spikes that
        # are element-scale noise, not a field.
        good = np.concatenate([[True], np.abs(np.diff(self.s)) > 1e-12])
        self.et = np.full_like(self.s, np.nan)
        if good.sum() > 2:
            o = np.argsort(self.s[good])
            ss, pp = self.s[good][o], self.phi[good][o]
            if smooth > 1 and pp.size > smooth:
                k = np.ones(int(smooth)) / float(int(smooth))
                pp = np.convolve(np.pad(pp, int(smooth) // 2, mode="edge"), k,
                                 mode="valid")[:pp.size]
            self.et[good] = np.interp(self.s[good], ss, -np.gradient(pp, ss))

    def at_strike(self, a):
        return a[self.istrike]

    def peak(self, a):
        m = np.isfinite(a)
        if not m.any():
            return np.nan, np.nan
        i = np.flatnonzero(m)[int(np.nanargmax(np.abs(a[m])))]
        return a[i], self.s[i]

    def integral(self, a, half=None):
        """int a ds over the target window - the extensive measure.

        The strike-point VALUE conflates two different things: more plasma on
        this leg, and the same plasma redistributed toward the strike point.
        HFSHD is a statement about the inner divertor's total content, so the
        integral is the metric that answers it.
        """
        m = np.isfinite(a) if half is None else (np.abs(self.s) <= half) & np.isfinite(a)
        if m.sum() < 2:
            return np.nan
        # s decreases with index on a flipped target, which would make trapz
        # return a negative line integral. Integrate along increasing s.
        o = np.argsort(self.s[m])
        return float(_trapz(a[m][o], self.s[m][o]))

    def reach(self):
        """Largest half-width over which this target is symmetric about s=0."""
        return min(abs(self.s.min()), abs(self.s.max()))

    def pfr_mean(self, a, width):
        m = (self.s > 0) & (self.s < width) & np.isfinite(a)
        return np.nanmean(a[m]) if m.any() else np.nan


def split_targets(d, rsplit, zmax, bnfrac, mass_amu, half_width, smooth,
                  verbose=False):
    div = d["Z"] < zmax
    if not div.any():
        sys.exit(f"no boundary points with Z < {zmax}; adjust --zmax")

    bn = np.abs(d.get("Bnorm", np.ones_like(d["R"])))
    thr = bnfrac * np.nanmax(bn[div])
    sel = div & (bn > thr)

    runs = contiguous_runs(sel)
    if verbose:
        print(f"  {len(runs)} contiguous boundary runs with Z<{zmax}, "
              f"|Bnorm|>{thr:.4g}:")
        for a, b in runs:
            print(f"    idx {a:5d}-{b:5d}  n={b-a:4d}  "
                  f"R {d['R'][a:b].min():.3f}-{d['R'][a:b].max():.3f}  "
                  f"Z {d['Z'][a:b].min():.3f}-{d['Z'][a:b].max():.3f}")

    def best(side):
        cand = [(a, b) for a, b in runs
                if (np.mean(d["R"][a:b]) < rsplit) == side]
        return max(cand, key=lambda ab: ab[1] - ab[0]) if cand else None

    out = []
    for side, name, pfr_larger_R in ((True, "INNER", True), (False, "OUTER", False)):
        ab = best(side)
        if ab is None:
            out.append(Target(name, d, slice(0, 0), mass_amu, pfr_larger_R,
                              half_width, smooth))
        else:
            out.append(Target(name, d, slice(*ab), mass_amu, pfr_larger_R,
                              half_width, smooth))
    return out


# ---------------------------------------------------------------------------
def report(label, step, t_now, inner, outer, pfr_width):
    print(f"\n{'=' * 76}\n{label}   step {step}   t_now = {t_now:.6g} (JOREK units)"
          f"\n{'=' * 76}")
    if not (inner.valid and outer.valid):
        print("  a target has too few points - rerun with --verbose and check "
              "--rsplit / --zmax / --bnfrac")
        return None

    for t in (inner, outer):
        if t.edge_warn:
            print(f"  WARNING: {t.name} strike point (min psi_N) sits at the edge of "
                  f"the selected run.\n           The real one is probably outside it - "
                  f"check --zmax / --bnfrac / --rsplit.")

    print(f"{'':24s}{'INNER':>15s}{'OUTER':>15s}{'in/out':>11s}")

    def row(nm, a, b):
        r = f"{a / b:11.3f}" if (b and np.isfinite(a / b)) else " " * 11
        print(f"  {nm:22s}{a:15.4g}{b:15.4g}{r}")

    for t in (inner, outer):
        t._sp = dict(phi=t.at_strike(t.phi), te=t.at_strike(t.te),
                     ti=t.at_strike(t.ti), ne=t.at_strike(t.ne),
                     R=t.R[t.istrike], Z=t.Z[t.istrike])

    print("\n  -- strike point (min psi_N) --")
    row("R [m]", inner._sp["R"], outer._sp["R"])
    row("Z [m]", inner._sp["Z"], outer._sp["Z"])
    row("psi_N", inner.psin[inner.istrike], outer.psin[outer.istrike])
    row("Phi [V]", inner._sp["phi"], outer._sp["phi"])
    row("T_e [eV]", inner._sp["te"], outer._sp["te"])
    row("T_i [eV]", inner._sp["ti"], outer._sp["ti"])
    row("n_e [1e20 m^-3]", inner._sp["ne"] / 1e20, outer._sp["ne"] / 1e20)

    print(f"\n  -- averaged over the PFR side, 0 < s < {pfr_width} m --")
    row("Phi [V]", inner.pfr_mean(inner.phi, pfr_width),
        outer.pfr_mean(outer.phi, pfr_width))
    row("T_e [eV]", inner.pfr_mean(inner.te, pfr_width),
        outer.pfr_mean(outer.te, pfr_width))
    row("n_e [1e20 m^-3]", inner.pfr_mean(inner.ne, pfr_width) / 1e20,
        outer.pfr_mean(outer.ne, pfr_width) / 1e20)
    row("E_t [V/m]", inner.pfr_mean(inner.et, pfr_width),
        outer.pfr_mean(outer.et, pfr_width))
    row("|v_ExB| [m/s]", inner.pfr_mean(inner.vexb, pfr_width),
        outer.pfr_mean(outer.vexb, pfr_width))

    print("\n  -- peak over the target --")
    for nm, at, sc in (("n_e [1e20 m^-3]", "ne", 1e20), ("T_e [eV]", "te", 1.0),
                       ("E_t [V/m]", "et", 1.0), ("|v_ExB| [m/s]", "vexb", 1.0)):
        a, sa = inner.peak(getattr(inner, at))
        b, sb = outer.peak(getattr(outer, at))
        row(nm, a / sc, b / sc)
        print(f"  {'':22s}{'s=' + format(sa, '+.3f'):>15s}"
              f"{'s=' + format(sb, '+.3f'):>15s}")

    # Extensive comparison. Both targets MUST be integrated over the same
    # window: the selected runs are not equally long about their strike points
    # (the outer one is routinely cut short), and integrating each over
    # whatever it happens to span turns a length ratio into a density ratio.
    hw = min(inner.reach(), outer.reach())
    li, lo = inner.integral(inner.ne, hw), outer.integral(outer.ne, hw)
    ji, jo = inner.integral(inner.jsat, hw), outer.integral(outer.jsat, hw)
    print(f"\n  -- integrated over a COMMON window, |s| < {hw:.3f} m --")
    print(f"  {'':22s}{'(inner reach ' + format(inner.reach(), '.3f'):>15s}"
          f"{'  outer ' + format(outer.reach(), '.3f') + ')':>15s}")
    row("int n_e ds [1e20/m^2]", li / 1e20, lo / 1e20)
    row("int j_sat ds [a.u.]", ji, jo)

    # The PFR-side average difference is the meaningful drive, not the value at
    # the strike point: the potential hill sits a few cm INTO the PFR, so the two
    # strike points can agree while a large gradient exists between them. The
    # mean E_t is a poor proxy too - it averages sign changes to near zero.
    dphi_pfr = (inner.pfr_mean(inner.phi, pfr_width)
                - outer.pfr_mean(outer.phi, pfr_width))
    print(f"\n  ** PFR-side potential difference (0 < s < {pfr_width} m) = "
          f"{dphi_pfr:+.2f} V **")

    dphi = inner._sp["phi"] - outer._sp["phi"]
    nr = inner._sp["ne"] / outer._sp["ne"]
    print(f"     (at the strike points themselves: {dphi:+.3f} V)")
    print(f"     n_e(in)/n_e(out) at the strike point = {nr:.3f}"
          "     (> 1 is the HFSHD direction)")
    lr = li / lo if lo else float("nan")
    print(f"     INTEGRATED  int n_in ds / int n_out ds = {lr:.3f}"
          "     <- the HFSHD metric")

    # --- Which way does the PFR ExB run?
    vt_i = inner.pfr_mean(inner.vt, pfr_width)
    vt_o = outer.pfr_mean(outer.vt, pfr_width)

    def coherence(t):
        m = np.hypot(t.pfr_mean(t.vt, pfr_width), t.pfr_mean(t.vn, pfr_width))
        return m / max(t.pfr_mean(t.vexb, pfr_width), 1e-30)

    coh_i, coh_o = coherence(inner), coherence(outer)
    print(f"\n  -- ExB direction at the targets (v_t > 0 = toward the PFR) --")
    row("v_t [m/s]", vt_i, vt_o)
    row("n*v_t [1e22 m^-2 s^-1]", inner.pfr_mean(inner.gam_t, pfr_width) / 1e22,
        outer.pfr_mean(outer.gam_t, pfr_width) / 1e22)
    row("v_n [m/s] (+ into wall)", inner.pfr_mean(inner.vn, pfr_width),
        outer.pfr_mean(outer.vn, pfr_width))
    # |mean(v)| / mean(|v|). 1 = the flow keeps one direction across the window;
    # near 0 = it rotates, and the mean components carry no information.
    row("direction coherence", coh_i, coh_o)

    COH_MIN = 0.5
    if min(coh_i, coh_o) < COH_MIN:
        bad = " and ".join(n for n, c in (("INNER", coh_i), ("OUTER", coh_o))
                           if c < COH_MIN)
        print(f"\n  ** DIRECTION NOT RESOLVED on {bad} (coherence < {COH_MIN}) **\n"
              f"     The ExB direction rotates within the averaging window, so the\n"
              f"     mean components are not a flow direction. The transfer between\n"
              f"     the targets happens in the PFR VOLUME anyway - a boundary trace\n"
              f"     only sees where field lines land. Use a pol_line traverse across\n"
              f"     the PFR instead (see util/postproc_pfr_line.in).")
    elif np.isfinite(vt_i) and np.isfinite(vt_o):
        if vt_i > 0 and vt_o < 0:
            verdict = ("enters the PFR at the INNER leg and leaves at the OUTER\n"
                       "     -> net PFR transport INNER -> OUTER, which OPPOSES HFSHD")
        elif vt_i < 0 and vt_o > 0:
            verdict = ("enters the PFR at the OUTER leg and leaves at the INNER\n"
                       "     -> net PFR transport OUTER -> INNER, the HFSHD DIRECTION")
        elif vt_i > 0 and vt_o > 0:
            verdict = ("converges into the PFR from BOTH legs\n"
                       "     -> no net inner/outer transfer; outflow must be elsewhere")
        else:
            verdict = ("diverges out of the PFR at BOTH legs\n"
                       "     -> no net inner/outer transfer; inflow must be elsewhere")
        print(f"\n  ** The ExB flow {verdict} **")

    return dphi_pfr, nr, lr


# ---------------------------------------------------------------------------
def plot(cases, outfile, smax):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib not available - text report only.")
        return

    panels = [("phi", r"$\Phi$  [V]", 1.0), ("te", r"$T_e$  [eV]", 1.0),
              ("ne", r"$n_e$  [$10^{20}$ m$^{-3}$]", 1e20),
              ("jsat", r"$j_{sat}$  [a.u.]", 1.0),
              ("et", r"$E_t=-\partial\Phi/\partial s$  [V/m]", 1.0),
              ("vexb", r"$|v_{E\times B}|$  [m/s]", 1.0),
              ("vt", r"$v_t$  [m/s]   (+ = toward PFR)", 1.0),
              ("gam_t", r"$n\,v_t$  [$10^{22}$ m$^{-2}$s$^{-1}$]", 1e22)]
    fig, axes = plt.subplots(2, 4, figsize=(19, 8), sharex=True)
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]

    for ax, (attr, ylab, sc) in zip(axes.ravel(), panels):
        for ci, (label, _, _, inner, outer) in enumerate(cases):
            c = colors[ci % len(colors)]
            for t, ls in ((inner, "-"), (outer, "--")):
                if t.valid:
                    m = np.abs(t.s) <= smax
                    ax.plot(t.s[m] * 100, getattr(t, attr)[m] / sc, ls, color=c,
                            lw=1.6, label=f"{label} {t.name}")
        ax.set_ylabel(ylab)
        ax.axvline(0.0, color="k", lw=0.8, ls=":")
        if attr in ("et", "vt", "gam_t"):
            ax.axhline(0.0, color="k", lw=0.8, ls="-", alpha=0.4)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("s from strike point [cm]      (s > 0 = toward the PFR)")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Divertor target profiles  -  solid = inner, dashed = outer", y=0.99)
    fig.tight_layout()
    fig.savefig(outfile, dpi=140)
    print(f"\nwrote {outfile}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--labels", nargs="*")
    ap.add_argument("--step", type=int)
    ap.add_argument("--rsplit", type=float, default=1.42)
    ap.add_argument("--zmax", type=float, default=-0.90)
    ap.add_argument("--bnfrac", type=float, default=0.10,
                    help="keep boundary points with |Bnorm| above this fraction "
                         "of its divertor maximum; rejects the wall (default 0.10)")
    ap.add_argument("--pfr-width", type=float, default=0.05,
                    help="PFR averaging window in m from the strike point")
    ap.add_argument("--half-width", type=float, default=0.12,
                    help="keep +-this many m of path around the strike point; the "
                         "boundary runs continuously from the target up the vessel "
                         "wall, so a Z cut alone cannot isolate the plate (default 0.12)")
    ap.add_argument("--smooth", type=int, default=9,
                    help="points to average Phi over before differentiating for "
                         "E_t (default 9); 1 disables")
    ap.add_argument("--smax", type=float, default=0.15, help="plot range in m")
    ap.add_argument("--mass", type=float, default=2.0141)
    ap.add_argument("-o", "--out", default="target_profiles.png")
    ap.add_argument("--verbose", action="store_true",
                    help="list the boundary runs found (use if selection looks wrong)")
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()

    labels = a.labels or [f.split("/")[-1] for f in a.files]
    if len(labels) != len(a.files):
        sys.exit("--labels needs one entry per file")

    cases, summ = [], []
    for path, label in zip(a.files, labels):
        step, t_now, d = pick_block(read_boundary_file(path), a.step)
        check_units(d, path)
        if a.verbose:
            print(f"\n{label}: {len(d['R'])} boundary points")
        inner, outer = split_targets(d, a.rsplit, a.zmax, a.bnfrac, a.mass,
                                     a.half_width, a.smooth, a.verbose)
        r = report(label, step, t_now, inner, outer, a.pfr_width)
        cases.append((label, step, t_now, inner, outer))
        if r:
            summ.append((label, *r))

    if len(summ) == 2:
        (l0, d0, n0, r0), (l1, d1, n1, r1) = summ
        print(f"\n{'=' * 76}\n  {l0}  ->  {l1}")
        print(f"    PFR drive (PFR-side mean)   : {d0:+.3f} V  ->  {d1:+.3f} V"
              f"   ({d1 - d0:+.3f} V)")
        print(f"    strike-point   n_in/n_out   : {n0:.3f}      ->  {n1:.3f}"
              f"       ({n1 - n0:+.3f})")
        print(f"    INTEGRATED     n_in/n_out   : {r0:.3f}      ->  {r1:.3f}"
              f"       ({r1 - r0:+.3f})")
        print(f"{'=' * 76}")

    if not a.no_plot:
        plot(cases, a.out, a.smax)


if __name__ == "__main__":
    main()
