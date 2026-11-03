#!/usr/bin/env python3
"""Compare two runs through their jorek2_postproc `boundary_quantities` output: a run with the sheath current
boundary condition and the thermoelectric terms against a reference run with the wall potential pinned.

  python3 util/compare_wall_runs.py util/thermo.dat util/nothermo.dat -o doc/elm_sheath_thermo

Writes the figures (PDF) into <outdir>/figs, the table of the pre-ELM state into <outdir>/table_pre.tex and the numbers
quoted in doc/elm_sheath_thermo/elm_sheath_thermo.tex into <outdir>/numbers.txt. Expressions and units as for util/plot_wall_current.py (SI, one toroidal period).
"""
import argparse
import math
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_wall_current import (MU0, AMU, C_INNER, C_OUTER, C_THIRD, INK, INK2, MUTED, DIVERGING, eich, eich_fit,  # noqa: E402
                               trapezoid)

E_CHG = 1.602176634e-19
GAMMA = 5.0 / 3.0
C_SH, C_REF = C_INNER, C_OUTER                     # one colour per run, fixed
L_SH, L_REF = "sheath + thermoelectric", "reference"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "legend.fontsize": 8, "lines.linewidth": 1.6,
                     "savefig.facecolor": "white", "figure.facecolor": "white", "axes.facecolor": "white"})


def read_fast(path):
    """names, steps, times [JOREK units], data(step, phi, point, column)."""
    names, steps, times, chunks, buf = None, [], [], [], []

    def flush():
        if buf:
            chunks.append(np.array(" ".join(buf).split(), dtype=float))
            buf.clear()

    with open(path) as f:
        for line in f:
            if line[0] == "#":
                if "time step" in line:
                    flush()
                    steps.append(int(line.split("#")[2].split(",")[0]))
                    times.append(float(line.split("=")[-1]))
                elif names is None:
                    names = line.lstrip("#").split()
            elif line.strip():
                buf.append(line)
    flush()
    if names is None:
        raise ValueError(f"{path}: no '# R Z ...' header line with the expression names")
    blocks = [c.reshape(-1, len(names)) for c in chunks]
    sizes = np.array([b.shape[0] for b in blocks])
    n_rows = int(np.bincount(sizes).argmax())                       # the usual block length
    keep = sizes == n_rows
    if not keep.all():
        bad = [steps[i] for i in np.where(~keep)[0]]
        print(f"{path}: {len(bad)} block(s) with a different number of points dropped (steps {bad[:6]}"
              f"{'...' if len(bad) > 6 else ''}); the file seems to hold output of more than one postproc run")
    blocks = [b for b, k in zip(blocks, keep) if k]
    steps, times = [x for x, k in zip(steps, keep) if k], [x for x, k in zip(times, keep) if k]
    # points per toroidal plane: the first index at which the opening sequence of points repeats and that divides the
    # block (a closed wall ends on its first point without starting a new plane)
    b0 = blocks[0]
    n_bnd = b0.shape[0]
    k = min(4, n_bnd)
    for i in np.where(np.all(np.isclose(b0[1:, :2], b0[0, :2], rtol=0, atol=1e-12), axis=1))[0] + 1:
        if n_bnd % i == 0 and i + k <= n_bnd and np.allclose(b0[i:i + k, :2], b0[:k, :2], rtol=0, atol=1e-12):
            n_bnd = int(i)
            break
    return names, np.array(steps), np.array(times), np.array([b.reshape(-1, n_bnd, len(names)) for b in blocks])


def phi_mean(a):
    """Toroidal average over axis 0 (one period, both end points included)."""
    return a[0] if a.shape[0] == 1 else trapezoid(a, axis=0) / (a.shape[0] - 1)


PLATE_SELECT = {"z_max": None, "max_s": None, "types": None}     # set by the scripts (--plate-zmax, --plate-length) for walls
                                                  # where more than the two targets see the field (AUG)


class Run:
    def __init__(self, path, label, color, t_unit, bnorm_frac=0.05, z_max=None, max_s=None):
        self.label, self.color = label, color
        z_max = PLATE_SELECT["z_max"] if z_max is None else z_max
        max_s = PLATE_SELECT["max_s"] if max_s is None else max_s
        names, self.steps, times, self.d = read_fast(path)
        self.c = {n: i for i, n in enumerate(names)}
        self.t = times * t_unit * 1e3                                    # ms
        self.t_jorek = times
        g = self.d[0, 0]
        self.R, self.Z, self.psi, self.dl = (g[:, self.c[k]] for k in ("R", "Z", "Psi_N", "bnd_dl"))
        bn = np.abs(self.d[0, :, :, self.c["Bnorm"]]).mean(axis=0)
        sel = bn > bnorm_frac * bn.max()
        if z_max is not None:
            sel &= self.Z < z_max
        tgt = np.where(sel)[0]
        Rs = np.sort(self.R[tgt])
        R_split = 0.5 * (Rs[np.argmax(np.diff(Rs))] + Rs[np.argmax(np.diff(Rs)) + 1])
        i_in, i_out = tgt[self.R[tgt] < R_split], tgt[self.R[tgt] > R_split]
        centre = lambda i: np.array([self.R[i].mean(), self.Z[i].mean()])
        self.plates = {"inner": self._plate(i_in, centre(i_out), max_s), "outer": self._plate(i_out, centre(i_in), max_s)}
        for pn in ("inner", "outer"):
            p = self.plates[pn]
            print(f"  {label}: {pn} plate {p['idx'].size} points, s = {p['s'].min() * 100:.1f}..{p['s'].max() * 100:.1f} cm, "
                  f"strike point R = {p['R'][np.argmin(np.abs(p['s']))]:.3f} m")

    def _plate(self, idx, other, max_s=None):
        x = np.column_stack([self.R[idx], self.Z[idx]])
        mean = x.mean(axis=0)
        direction = np.linalg.svd(x - mean)[2][0]
        if np.dot(mean - other, direction) < 0.0:
            direction = -direction
        idx = idx[np.argsort((x - mean) @ direction)]
        s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(self.R[idx]), np.diff(self.Z[idx])))])
        s -= s[np.argmin(self.psi[idx])]
        if max_s is not None:                                           # the target only, around its strike point
            keep = np.abs(s) <= max_s
            idx, s = idx[keep], s[keep]
        return {"idx": idx, "s": s, "R": self.R[idx], "ring": 2.0 * np.pi * self.R[idx] * self.dl[idx]}

    def transition_mask(self, plate, i=0, points_per_element=4):
        """Plate points of the second element from each plate end: the sheath row acts there, but the node it shares
        with the corner element is not released (see row_mask)."""
        n = self.plates[plate]["idx"].size
        first = last = points_per_element
        if "bnd_type_lo" in self.c and "bnd_type_hi" in self.c:
            corner = self.q("bnd_type_lo", plate, i)[0] != self.q("bnd_type_hi", plate, i)[0]
            first = int(np.argmin(corner)) if corner[0] else 0
            last = int(np.argmin(corner[::-1])) if corner[-1] else 0
        mask = np.zeros(n, bool)
        mask[first: 2 * first] = True
        if last:
            mask[n - 2 * last: n - last] = True
        return mask

    def row_mask(self, plate, i=0, points_per_element=4, n_elements=2):
        """Plate points where the sheath row acts with both of its end nodes released.

        A plate node is released only if every wall element it belongs to carries the row, so the element next to a
        corner (no row: one node is the corner) AND the element after it (row, but its outer node is pinned to the
        floating potential and the frozen current) are left out at each end of the plate. With the expressions
        bnd_type_lo / bnd_type_hi in the file the corner elements are found from the node types; otherwise each
        element is taken to have points_per_element points (nsub_bnd)."""
        n = self.plates[plate]["idx"].size
        keep = np.ones(n, bool)
        if "bnd_type_lo" in self.c and "bnd_type_hi" in self.c:
            corner = self.q("bnd_type_lo", plate, i)[0] != self.q("bnd_type_hi", plate, i)[0]
            first = int(np.argmin(corner)) if corner[0] else 0              # corner points at the start of the plate
            last = int(np.argmin(corner[::-1])) if corner[-1] else 0        # ... and at its end
            if corner.all():
                return ~keep
            first, last = first or points_per_element, last or points_per_element
        else:
            first = last = points_per_element
        if n_elements > 0:                                  # n_elements at each end of the plate are left out
            keep[: n_elements * first] = False
            keep[n - n_elements * last:] = False
        return keep

    def type_mask(self, plate, i=0, types=None, points_per_element=4, n_elements=1):
        """Plate points on wall elements of the given node types (both end nodes), minus n_elements on either side
        of every change of type (the junction nodes are pinned like the corners). types None: PLATE_SELECT['types'];
        still None, or no bnd_type columns in the file: every point."""
        n = self.plates[plate]["idx"].size
        keep = np.ones(n, bool)
        types = PLATE_SELECT["types"] if types is None else types
        if "bnd_type_lo" not in self.c or "bnd_type_hi" not in self.c:
            return keep
        lo, hi = self.q("bnd_type_lo", plate, i)[0], self.q("bnd_type_hi", plate, i)[0]
        if types is not None:
            keep &= np.isin(lo, types) & np.isin(hi, types)
        change = np.where((np.diff(lo) != 0) | (np.diff(hi) != 0) | (lo[:-1] != hi[:-1]))[0]
        for j in change:                                              # within n_elements of a junction
            keep[max(0, j - n_elements * points_per_element + 1): j + n_elements * points_per_element + 1] = False
        return keep

    def q(self, name, plate, i):
        """(n_phi, n_point) of an expression on a plate at step index i."""
        return self.d[i][:, self.plates[plate]["idx"], self.c[name]]

    def avg(self, name, plate, i):
        return phi_mean(self.q(name, plate, i))

    def integral(self, name, plate):
        p = self.plates[plate]
        return np.array([np.sum(self.avg(name, plate, i) * p["ring"]) for i in range(len(self.steps))])

    def trace(self, name, plate, fun=np.max):
        return np.array([fun(self.avg(name, plate, i)) for i in range(len(self.steps))])

    def fluct(self, name, plate):
        """rms over the plate and the toroidal angle of the deviation from the toroidal mean."""
        out = np.empty(len(self.steps))
        for i in range(len(self.steps)):
            a = self.q(name, plate, i)[:-1]
            out[i] = np.sqrt(np.mean((a - a.mean(axis=0)) ** 2))
        return out

    def index_at(self, t_ms):
        return int(np.argmin(np.abs(self.t - t_ms)))

    def exb_normal(self, plate, i):
        """ExB velocity through the plate [m/s], positive out of the plasma, from the toroidally averaged potential:
        v_E.n = -(dPhi/ds)/B_phi for a tangent t with outward normal (t_Z, -t_R); B_phi = F0/R (lower divertor plates)."""
        p = self.plates[plate]
        R, Z, s = self.R[p["idx"]], self.Z[p["idx"]], p["s"]
        tR, tZ = np.gradient(R, s), np.gradient(Z, s)
        sign = np.where(tR > 0.0, 1.0, -1.0)      # (t_Z, -t_R) is the outward normal of a plate at the bottom if t_R > 0
        return -np.gradient(self.avg("Phi", plate, i), s) / (self.F0 / R) * sign

    def bohm_normal(self, plate, i, mass):
        """cs |b.n| [m/s] with cs = sqrt(gamma (Te + Ti)/m_i), as in the sheath row."""
        p = self.plates[plate]
        cs = np.sqrt(GAMMA * E_CHG * (self.avg("T_e", plate, i) + self.avg("T_i", plate, i)) / mass)
        bn = self.avg("Bnorm", plate, i)
        return cs * np.abs(bn) / np.hypot(self.F0 / p["R"], bn)


def growth(run, plate, lo=1e-7, hi=5e-2):
    """Exponential fit of the T_e footprint amplitude [eV] in its linear phase; returns gamma [1/t_A], amplitude."""
    amp = run.fluct("T_e", plate)
    m = (amp > lo) & (amp < hi) & (run.t_jorek > run.t_jorek[0] + 300.0)
    g, b = np.polyfit(run.t_jorek[m], np.log(amp[m]), 1)
    return g, b, amp, m


def time_of(amp, t, level):
    j = int(np.argmax(amp > level))
    return float(np.interp(math.log(level), np.log(amp[j - 1:j + 1]), t[j - 1:j + 1]))


def panel_label(ax, text):
    ax.text(0.02, 0.96, text, transform=ax.transAxes, ha="left", va="top", fontsize=9, fontweight="bold", color=INK)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sheath_file")
    ap.add_argument("reference_file")
    ap.add_argument("-o", "--outdir", default="doc/elm_sheath_thermo")
    ap.add_argument("--F0", type=float, default=3.0)
    ap.add_argument("--central-density", type=float, default=1.0)
    ap.add_argument("--central-mass", type=float, default=2.01410177811)
    ap.add_argument("--sheath-lambda", type=float, default=3.0)
    ap.add_argument("--t-pre", type=float, default=1.388, help="time [ms] of the pre-ELM comparison")
    args = ap.parse_args()
    plt.switch_backend("Agg")                 # this script only writes files (not at import: other scripts use Run)
    figdir = os.path.join(args.outdir, "figs")
    os.makedirs(figdir, exist_ok=True)
    mass = args.central_mass * AMU
    t_unit = math.sqrt(MU0 * args.central_density * 1e20 * mass)

    sh = Run(args.sheath_file, L_SH, C_SH, t_unit)
    ref = Run(args.reference_file, L_REF, C_REF, t_unit)
    runs = (sh, ref)
    for r in runs:
        r.F0 = args.F0
    plates = ("inner", "outer")
    N = {}                                                              # numbers for the text

    def save(fig, name):
        fig.savefig(os.path.join(figdir, name + ".pdf"))
        plt.close(fig)

    # ------------------------------------------------------------------ geometry
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.4), constrained_layout=True, gridspec_kw={"width_ratios": [1, 1.5]})
    for a in ax:
        a.plot(sh.R, sh.Z, ".", ms=1.0, color=MUTED)
        for pn, col in zip(plates, (C_INNER, C_OUTER)):
            p = sh.plates[pn]
            a.plot(p["R"], sh.Z[p["idx"]], color=col, lw=3, label=f"{pn} plate")
            k = int(np.argmin(np.abs(p["s"])))
            a.plot(p["R"][k], sh.Z[p["idx"]][k], "o", ms=5, mfc="white", mec=INK, mew=1.2)
        a.set_aspect("equal")
        a.set_xlabel("R [m]")
    ax[0].set_ylabel("Z [m]")
    ax[0].legend(loc="center")
    ax[1].set_xlim(2.35, 3.6)
    ax[1].set_ylim(-1.72, -1.45)
    ax[1].set_title("divertor region (circles: strike points)")
    ax[0].set_title("boundary of the computational domain")
    save(fig, "geometry")
    for pn in plates:
        p = sh.plates[pn]
        N[f"len_pfr_{pn}"] = -p["s"].min() * 100
        N[f"len_sol_{pn}"] = p["s"].max() * 100
        k = int(np.argmin(np.abs(p["s"])))
        N[f"angle_{pn}"] = math.degrees(math.asin(abs(sh.avg("Bnorm", pn, 0)[k]) / math.hypot(args.F0 / p["R"][k],
                                                                                         sh.avg("Bnorm", pn, 0)[k])))
        N[f"Rstrike_{pn}"] = p["R"][k]

    # ------------------------------------------------------------------ growth of the footprint
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0), sharey=True, constrained_layout=True)
    for k, pn in enumerate(plates):
        for r in runs:
            g, b, amp, m = growth(r, pn)
            tag = "sh" if r is sh else "ref"
            N[f"gamma_{tag}_{pn}"] = g
            N[f"gamma_si_{tag}_{pn}"] = g / t_unit
            N[f"t1eV_{tag}_{pn}"] = time_of(amp, r.t, 1.0) if amp.max() > 1.0 else float("nan")
            N[f"t01eV_{tag}_{pn}"] = time_of(amp, r.t, 0.1)
            N[f"ampend_{tag}_{pn}"] = amp[-1]
            ok = amp > 1e-9
            ax[k].semilogy(r.t[ok], amp[ok], color=r.color, label=f"{r.label}: $\\gamma\\tau_A$ = {g:.4f}")
            ax[k].semilogy(r.t[m], np.exp(g * r.t_jorek[m] + b), color=INK, lw=0.8, ls=(0, (4, 2)))
        ax[k].set_title(f"{pn} plate")
        ax[k].set_xlabel("time [ms]")
        ax[k].legend(loc="upper left")
    ax[0].set_ylabel(r"rms $\tilde T_e$ on the plate [eV]")
    save(fig, "growth")
    t_on = {"sh": N["t1eV_sh_outer"], "ref": N["t1eV_ref_outer"]}       # onset: 1 eV rms on the outer plate
    N["lead_ms"] = t_on["ref"] - t_on["sh"]
    N["gamma_ratio"] = N["gamma_sh_outer"] / N["gamma_ref_outer"]

    # ------------------------------------------------------------------ time traces
    rows = [("heatF_total", 1e-6, "power to the plate [MW]", "integral"),
            ("Jn_sat", 1e-3, "ion saturation current [kA]", "integral"),
            ("T_e", 1.0, "peak $T_e$ [eV]", "max"),
            ("ne", 1e-19, "peak $n_e$ [$10^{19}$ m$^{-3}$]", "max"),
            ("Jn_par", 1e-3, "plate current [kA]", "integral")]
    fig, ax = plt.subplots(len(rows), 2, figsize=(7.0, 9.0), sharex=True, constrained_layout=True)
    for j, (name, fac, ylabel, kind) in enumerate(rows):
        for k, pn in enumerate(plates):
            for r in runs:
                y = (r.integral(name, pn) if kind == "integral" else r.trace(name, pn)) * fac
                ax[j, k].plot(r.t, y, color=r.color, label=r.label)
                ax[j, k].axvline(t_on["sh" if r is sh else "ref"], color=r.color, lw=0.7, ls=":")
            if name == "Jn_par":
                ax[j, k].axhline(0.0, color=MUTED, lw=0.6)
            if k == 0:
                ax[j, k].set_ylabel(ylabel)
            if j == 0:
                ax[j, k].set_title(f"{pn} plate")
    ax[0, 0].legend(loc="upper left")
    for a in ax[-1]:
        a.set_xlabel("time [ms]")
    save(fig, "traces")

    # ------------------------------------------------------------------ table of the pre-ELM state
    i_pre = {r: r.index_at(args.t_pre) for r in runs}
    N["t_pre"] = sh.t[i_pre[sh]]
    table = []
    for pn in plates:
        for r in runs:
            tag = "sh" if r is sh else "ref"
            i = i_pre[r]
            P = r.integral("heatF_total", pn)[i] * 1e-6
            Isat = r.integral("Jn_sat", pn)[i] * 1e-3
            Ipar = r.integral("Jn_par", pn)[i] * 1e-3
            te, ti, ne = r.avg("T_e", pn, i), r.avg("T_i", pn, i), r.avg("ne", pn, i)
            q = r.avg("heatF_total", pn, i)
            fit = eich_fit(r.plates[pn]["s"], q)
            vals = dict(P=P, Isat=Isat, Ipar=Ipar, Te=te.max(), Ti=ti.max(), ne=ne.max() * 1e-19, neTe=ne[np.argmax(te)] * 1e-19,
                        qmax=(eich(np.linspace(-0.05, 0.15, 800), *fit).max() if fit else np.nan) * 1e-6,
                        lam=fit[1] * 1e3 if fit else np.nan, S=fit[2] * 1e3 if fit else np.nan,
                        dphi_min=r.avg("dPhi_float", pn, i).min(), dphi_max=r.avg("dPhi_float", pn, i).max(),
                        jmin=r.avg("Jn_par", pn, i).min() * 1e-3, jmax=r.avg("Jn_par", pn, i).max() * 1e-3,
                        jsh_min=r.avg("Jn_sheath", pn, i).min() * 1e-3,
                        jsh_ratio=np.min(r.avg("Jn_sheath", pn, i) / np.maximum(r.avg("Jn_sat", pn, i), 1e-30)))
            for key, v in vals.items():
                N[f"pre_{key}_{tag}_{pn}"] = v
            table.append((pn, r.label, vals))
    for key in ("P", "Isat"):
        for pn in plates:
            N[f"pre_{key}_ratio_{pn}"] = N[f"pre_{key}_sh_{pn}"] / N[f"pre_{key}_ref_{pn}"]
        for tag in ("sh", "ref"):
            N[f"pre_{key}_tot_{tag}"] = N[f"pre_{key}_{tag}_inner"] + N[f"pre_{key}_{tag}_outer"]
            N[f"pre_{key}_asym_{tag}"] = N[f"pre_{key}_{tag}_outer"] / N[f"pre_{key}_{tag}_inner"]
    with open(os.path.join(args.outdir, "table_pre.tex"), "w") as f:
        f.write("\\begin{tabular}{llrrrrrrrr}\n\\toprule\nplate & run & $P$ & $I_\\mathrm{sat}$ & $I_\\parallel$ & $T_e^\\mathrm{max}$ & "
                "$T_i^\\mathrm{max}$ & $n_e^\\mathrm{max}$ & $q^\\mathrm{max}$ & $\\lambda_q^\\mathrm{t}$, $S^\\mathrm{t}$ \\\\\n"
                " & & MW & kA & kA & eV & eV & $10^{19}\\,$m$^{-3}$ & MW/m$^2$ & mm \\\\\n\\midrule\n")
        for pn, lab, v in table:
            f.write(f"{pn} & {lab} & {v['P']:.2f} & {v['Isat']:.1f} & {v['Ipar']:.2f} & {v['Te']:.1f} & {v['Ti']:.1f} & "
                    f"{v['ne']:.3f} & {v['qmax']:.2f} & {v['lam']:.0f}, {v['S']:.0f} \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")

    # ------------------------------------------------------------------ pre-ELM profiles of the plasma state
    rows = [("ne", 1e-19, "$n_e$ [$10^{19}$ m$^{-3}$]"), ("T_e", 1.0, "$T_e$ [eV]"), ("T_i", 1.0, "$T_i$ [eV]"),
            ("Phi", 1.0, r"$\Phi$ [V]"), ("heatF_total", 1e-6, "heat flux [MW/m$^2$]")]
    fig, ax = plt.subplots(len(rows), 2, figsize=(7.0, 9.0), sharex="col", constrained_layout=True)
    for j, (name, fac, ylabel) in enumerate(rows):
        for k, pn in enumerate(plates):
            a = ax[j, k]
            for r in runs:
                s = r.plates[pn]["s"] * 100
                y = r.avg(name, pn, i_pre[r]) * fac
                if name == "heatF_total":
                    a.plot(s, y, ls="none", marker="o", ms=1.6, color=r.color)
                    fit = eich_fit(r.plates[pn]["s"], y / fac)
                    if fit:
                        fine = np.linspace(r.plates[pn]["s"][0], r.plates[pn]["s"][-1], 400)
                        a.plot(fine * 100, eich(fine, *fit) * fac, color=r.color, label=r.label + " (Eich fit)")
                else:
                    a.plot(s, y, color=r.color, label=r.label)
                if name == "Phi" and r is sh:
                    a.plot(s, args.sheath_lambda * r.avg("T_e", pn, i_pre[r]), color=r.color, lw=0.9, ls=(0, (4, 2)),
                           label=r"floating, $\Lambda T_e/e$")
            a.axvline(0.0, color=MUTED, lw=0.6, ls=":")
            if k == 0:
                a.set_ylabel(ylabel)
            if j == 0:
                a.set_title(f"{pn} plate, t = {N['t_pre']:.2f} ms")
    ax[0, 0].legend(loc="upper right")
    ax[3, 1].legend(loc="upper right", fontsize=7)
    for a in ax[-1]:
        a.set_xlabel("distance from the strike point [cm]")
    save(fig, "profiles_state")

    # ------------------------------------------------------------------ current, sheath state and ExB flow
    fig, ax = plt.subplots(3, 2, figsize=(7.0, 6.6), sharex="col", constrained_layout=True)
    for k, pn in enumerate(plates):
        i = i_pre[sh]
        s = sh.plates[pn]["s"] * 100
        jsat = sh.avg("Jn_sat", pn, i) * 1e-3
        a = ax[0, k]
        a.fill_between(s, -jsat, jsat, color="#e4e3df", lw=0, label=r"$\pm j_\mathrm{sat}$")
        a.plot(s, sh.avg("Jnorm", pn, i) * 1e-3, color=C_THIRD, lw=1.0, label="incl. diamagnetic current")
        a.plot(s, sh.avg("Jn_sheath", pn, i) * 1e-3, color=INK, lw=0.9, ls=(0, (4, 2)), label="sheath characteristic")
        a.plot(s, sh.avg("Jn_par", pn, i) * 1e-3, color=C_SH, label="parallel current")
        a.axhline(0.0, color=MUTED, lw=0.6)
        a.set_title(f"{pn} plate, t = {N['t_pre']:.2f} ms")
        a = ax[1, k]
        for r in runs:
            sr = r.plates[pn]["s"] * 100
            a.plot(sr, r.avg("dPhi_float", pn, i_pre[r]) / r.avg("T_e", pn, i_pre[r]), color=r.color, label=r.label)
        a.axhline(0.0, color=MUTED, lw=0.6)
        a.axhline(-args.sheath_lambda, color=MUTED, lw=0.6, ls=(0, (1, 2)))
        a = ax[2, k]
        a.plot(s, sh.bohm_normal(pn, i, mass) * 1e-3, color=INK2, lw=1.0, label=r"Bohm, $c_s|\mathbf{b}\cdot\mathbf{n}|$")
        a.plot(s, sh.exb_normal(pn, i) * 1e-3, color=C_SH, label=r"$\mathbf{v}_E\cdot\mathbf{n}$ (sheath run)")
        a.axhline(0.0, color=MUTED, lw=0.6)
        a.set_xlabel("distance from the strike point [cm]")
        for a in ax[:, k]:
            a.axvline(0.0, color=MUTED, lw=0.6, ls=":")
        vE, vB = sh.exb_normal(pn, i), sh.bohm_normal(pn, i, mass)
        ne, ring = sh.avg("ne", pn, i), sh.plates[pn]["ring"]
        N[f"vE_max_{pn}"] = np.abs(vE).max() * 1e-3
        N[f"vE_ratio_{pn}"] = np.max(np.abs(vE) / vB)
        N[f"vB_strike_{pn}"] = vB[np.argmin(np.abs(sh.plates[pn]["s"]))] * 1e-3
        N[f"G_par_{pn}"] = np.sum(ne * vB * ring)
        N[f"G_E_out_{pn}"] = np.sum(ne * np.clip(vE, 0, None) * ring)
        N[f"G_E_in_{pn}"] = -np.sum(ne * np.clip(vE, None, 0) * ring)
        pE = 2.5 * ne * E_CHG * (sh.avg("T_e", pn, i) + sh.avg("T_i", pn, i)) * vE
        N[f"P_E_net_{pn}"] = np.sum(pE * ring) * 1e-6
        N[f"phi_max_{pn}"] = sh.avg("Phi", pn, i).max()
    ax[0, 0].set_ylabel("current into the plate [kA/m$^2$]")
    ax[1, 0].set_ylabel(r"$e(\Phi-\Phi_\mathrm{fl})/T_e$")
    ax[2, 0].set_ylabel("velocity through the plate [km/s]")
    ax[0, 0].legend(loc="lower right", fontsize=7)
    ax[1, 1].legend(loc="center right", fontsize=7)
    ax[2, 1].legend(loc="center right", fontsize=7)
    save(fig, "current_sheath")

    # ------------------------------------------------------------------ the two ends of the SOL flux tubes
    def on_flux(r, name, pn, i, psig, weight=None):
        """Expression on the SOL side of a plate, interpolated to Psi_N; with weight="flux": current per unit Psi_N."""
        p = r.plates[pn]
        m = p["s"] > 0.0
        psi, sp = r.psi[p["idx"]][m], p["s"][m]
        o = np.argsort(psi)
        y = r.avg(name, pn, i)[m]
        if weight == "flux":
            y = y * 2.0 * np.pi * p["R"][m] * np.gradient(sp[o], psi[o])[np.argsort(o)]
        return np.interp(psig, psi[o], y[o])

    psig = np.linspace(1.003, 1.21, 120)
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.9), constrained_layout=True)
    for r in runs:
        for pn, ls in (("inner", "-"), ("outer", (0, (4, 2)))):
            ax[0].plot(psig, on_flux(r, "T_e", pn, i_pre[r], psig), color=r.color, ls=ls)
            ax[1].semilogy(psig, on_flux(r, "ne", pn, i_pre[r], psig) * 1e-19, color=r.color, ls=ls,
                           label=f"{r.label}, {pn}")
    dI = {pn: on_flux(sh, "Jn_par", pn, i_pre[sh], psig, "flux") * 1e-3 for pn in plates}
    ax[2].plot(psig, dI["inner"], color=C_INNER, label="inner plate")
    ax[2].plot(psig, dI["outer"], color=C_OUTER, label="outer plate")
    ax[2].plot(psig, 0.5 * (dI["outer"] - dI["inner"]), color=INK, lw=1.0, ls=(0, (4, 2)), label="through-current")
    ax[2].axhline(0.0, color=MUTED, lw=0.6)
    ax[0].set(xlabel=r"$\psi_N$", ylabel="$T_e$ at the plate [eV]")
    ax[1].set(xlabel=r"$\psi_N$", ylabel="$n_e$ at the plate [$10^{19}$ m$^{-3}$]")
    ax[2].set(xlabel=r"$\psi_N$", ylabel=r"$\mathrm{d}I/\mathrm{d}\psi_N$ into the plate [kA]")
    ax[1].set_ylim(1e-3, 0.8)
    ax[1].legend(loc="lower left", fontsize=6)
    ax[2].legend(loc="lower right", fontsize=6.5)
    for a, lab in zip(ax, "abc"):
        a.set_title(f"({lab})", loc="left")
    save(fig, "flux_tubes")
    N["I_through_sol"] = trapezoid(0.5 * (dI["outer"] - dI["inner"]), psig)
    N["I_common_sol"] = trapezoid(dI["outer"] + dI["inner"], psig)
    for r in runs:
        tag = "sh" if r is sh else "ref"
        for pn in plates:
            p = r.plates[pn]
            m = p["s"] > 0.0
            for name, fac, key in (("Jn_par", 1e-3, "I"), ("Jn_sat", 1e-3, "Isat"), ("heatF_total", 1e-6, "P")):
                y = r.avg(name, pn, i_pre[r]) * p["ring"] * fac
                N[f"{key}_sol_{tag}_{pn}"], N[f"{key}_pfr_{tag}_{pn}"] = y[m].sum(), y[~m].sum()
        far = (psig > 1.06) & (psig < 1.16)
        N[f"dTe_far_{tag}"] = np.mean(on_flux(r, "T_e", "inner", i_pre[r], psig)[far] - on_flux(r, "T_e", "outer", i_pre[r], psig)[far])
        k110 = int(np.argmin(np.abs(psig - 1.10)))
        N[f"ne_ratio_110_{tag}"] = on_flux(r, "ne", "outer", i_pre[r], psig)[k110] / on_flux(r, "ne", "inner", i_pre[r], psig)[k110]
        N[f"ne_in_110_{tag}"] = on_flux(r, "ne", "inner", i_pre[r], psig)[k110] * 1e-19
        N[f"ne_out_110_{tag}"] = on_flux(r, "ne", "outer", i_pre[r], psig)[k110] * 1e-19

    # ------------------------------------------------------------------ net current and its closure
    I = {pn: sh.integral("Jn_par", pn) * 1e-3 for pn in plates}
    In = {pn: sh.integral("Jnorm", pn) * 1e-3 for pn in plates}
    N["Inet_start"], N["Inet_pre"], N["Inet_end"] = [(I["inner"] + I["outer"])[i] for i in (0, i_pre[sh], -1)]
    N["Inorm_net_pre"] = (In["inner"] + In["outer"])[i_pre[sh]]
    for pn in plates:
        N[f"I_start_{pn}"], N[f"I_pre_{pn}"], N[f"I_end_{pn}"] = I[pn][0], I[pn][i_pre[sh]], I[pn][-1]

    # ------------------------------------------------------------------ response of the currents to the mode
    fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.2), constrained_layout=True)
    a = ax[0, 0]
    a.plot(sh.t, I["inner"], color=C_INNER, label="inner plate")
    a.plot(sh.t, I["outer"], color=C_OUTER, label="outer plate")
    a.plot(sh.t, I["inner"] + I["outer"], color=INK2, lw=1.0, label="sum")
    a.axhline(0.0, color=MUTED, lw=0.6)
    a.axvline(t_on["sh"], color=MUTED, lw=0.7, ls=":")
    a.set(xlabel="time [ms]", ylabel="plate current [kA]", title="toroidally averaged plate current")
    a.legend(loc="center left")
    a = ax[0, 1]
    for pn, col in zip(plates, (C_INNER, C_OUTER)):
        fl = sh.fluct("Jn_par", pn) * 1e-3
        ok = fl > 1e-9
        a.semilogy(sh.t[ok], fl[ok], color=col, label=f"{pn} plate")
        N[f"dJ_end_{pn}"] = fl[-1]
        N[f"dPhi_end_{pn}"] = sh.fluct("Phi", pn)[-1]
    a.axvline(t_on["sh"], color=MUTED, lw=0.7, ls=":")
    a.set(xlabel="time [ms]", ylabel=r"rms $\tilde j_n$ [kA/m$^2$]", title="non-axisymmetric plate current")
    a.legend(loc="upper left")
    i_a, i_b = i_pre[sh], len(sh.steps) - 1
    for k, pn in enumerate(plates):
        a = ax[1, k]
        s = sh.plates[pn]["s"] * 100
        jb = sh.q("Jn_par", pn, i_b) * 1e-3
        a.fill_between(s, jb.min(axis=0), jb.max(axis=0), color=sh.color, alpha=0.25, lw=0, label="range over the toroidal angle")
        a.plot(s, sh.avg("Jn_par", pn, i_a) * 1e-3, color=INK2, lw=1.0, ls=(0, (4, 2)), label=f"t = {sh.t[i_a]:.2f} ms")
        a.plot(s, sh.avg("Jn_par", pn, i_b) * 1e-3, color=sh.color, label=f"t = {sh.t[i_b]:.2f} ms")
        a.axhline(0.0, color=MUTED, lw=0.6)
        a.axvline(0.0, color=MUTED, lw=0.6, ls=":")
        a.set(xlabel="distance from the strike point [cm]", title=f"{pn} plate")
        a.legend(loc="upper right" if pn == "inner" else "lower right", fontsize=7)
        N[f"jmin_end_{pn}"] = sh.avg("Jn_par", pn, i_b).min() * 1e-3
        N[f"jmax_end_{pn}"] = sh.avg("Jn_par", pn, i_b).max() * 1e-3
        N[f"jrange_end_{pn}"] = np.max(jb.max(axis=0) - jb.min(axis=0))
    ax[1, 0].set_ylabel("current into the plate [kA/m$^2$]")
    save(fig, "current_response")

    # ------------------------------------------------------------------ footprints at equal mode amplitude
    i_m = {sh: sh.index_at(t_on["sh"]), ref: ref.index_at(t_on["ref"])}
    N["t_match_sh"], N["t_match_ref"] = sh.t[i_m[sh]], ref.t[i_m[ref]]
    phi_deg = np.linspace(0.0, 60.0, sh.d.shape[1])
    fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.6), constrained_layout=True)
    panels = [(ax[0, 0], sh, "T_e", 1.0, r"$\tilde T_e$ [eV], " + L_SH), (ax[0, 1], ref, "T_e", 1.0, r"$\tilde T_e$ [eV], " + L_REF),
              (ax[1, 0], sh, "Jn_par", 1e-3, r"$\tilde j_n$ [kA/m$^2$], " + L_SH), (ax[1, 1], sh, "Phi", 1.0, r"$\tilde\Phi$ [V], " + L_SH)]
    vmax_T = max(np.abs(r.q("T_e", "outer", i_m[r]) - r.avg("T_e", "outer", i_m[r])).max() for r in runs)
    for a, r, name, fac, title in panels:
        dev = (r.q(name, "outer", i_m[r]) - r.avg(name, "outer", i_m[r])) * fac
        vmax = vmax_T if name == "T_e" else np.abs(dev).max()
        m = a.pcolormesh(r.plates["outer"]["s"] * 100, phi_deg, dev, cmap=DIVERGING, vmin=-vmax, vmax=vmax, shading="nearest",
                         rasterized=True)
        a.grid(False)
        a.set_title(title)
        a.set_xlim(-3, 10)
        fig.colorbar(m, ax=a)
    for a in ax[:, 0]:
        a.set_ylabel("toroidal angle [deg]")
    for a in ax[1]:
        a.set_xlabel("distance from the strike point [cm]")
    save(fig, "footprints")

    # n = 6 component (first harmonic of the period) of the footprints on the outer plate of the sheath run:
    # amplitude ratios and phases against the temperature footprint where that is largest
    def harmonic(r, name, pn, i):
        a = r.q(name, pn, i)[:-1]
        k = np.arange(a.shape[0])
        return 2.0 * np.mean(a * np.exp(-2j * np.pi * k / a.shape[0])[:, None], axis=0)

    hT, hJ, hP = (harmonic(sh, n_, "outer", i_m[sh]) for n_ in ("T_e", "Jn_par", "Phi"))
    hS = harmonic(sh, "Jn_sat", "outer", i_m[sh])
    kmax = int(np.argmax(np.abs(hT)))
    N["foot_s_Tmax"] = sh.plates["outer"]["s"][kmax] * 100
    N["foot_T_amp"] = abs(hT[kmax])
    N["foot_J_amp"] = abs(hJ[kmax]) * 1e-3
    N["foot_J_phase"] = math.degrees(np.angle(hJ[kmax] / hT[kmax]))
    N["foot_Phi_amp"] = abs(hP[kmax])
    N["foot_Phi_phase"] = math.degrees(np.angle(hP[kmax] / hT[kmax]))
    N["foot_J_max"] = np.abs(hJ).max() * 1e-3
    N["foot_Jsat_amp"] = abs(hS[kmax]) * 1e-3
    N["foot_Phi_max"] = np.abs(hP).max()
    N["foot_s_Phimax"] = sh.plates["outer"]["s"][int(np.argmax(np.abs(hP)))] * 100
    N["foot_jsat_local"] = sh.avg("Jn_sat", "outer", i_m[sh])[kmax] * 1e-3
    # linear response of the characteristic at that point: dj = jsat e^x (e dPhi/Te) + ...
    w = np.abs(hT) > 0.3 * np.abs(hT).max()
    N["foot_corr_JT"] = float(np.real(np.vdot(hT[w], hJ[w])) / (np.linalg.norm(hT[w]) * np.linalg.norm(hJ[w])))
    N["foot_corr_PhiT"] = float(np.real(np.vdot(hT[w], hP[w])) / (np.linalg.norm(hT[w]) * np.linalg.norm(hP[w])))

    # radial structure of the n = 6 perturbation at equal amplitude
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8), sharey=True, constrained_layout=True)
    for k, pn in enumerate(plates):
        for r in runs:
            a = r.q("T_e", pn, i_m[r])[:-1]
            ax[k].plot(r.plates[pn]["s"] * 100, np.sqrt(np.mean((a - a.mean(axis=0)) ** 2, axis=0)), color=r.color, label=r.label)
            tag = "sh" if r is sh else "ref"
            rms = np.sqrt(np.mean((a - a.mean(axis=0)) ** 2, axis=0))
            N[f"spk_{tag}_{pn}"] = r.plates[pn]["s"][np.argmax(rms)] * 100
        ax[k].axvline(0.0, color=MUTED, lw=0.6, ls=":")
        ax[k].set(xlabel="distance from the strike point [cm]", title=f"{pn} plate")
    ax[0].set_ylabel(r"toroidal rms of $\tilde T_e$ [eV]")
    ax[1].legend(loc="upper right")
    save(fig, "mode_structure")

    # ------------------------------------------------------------------ heat-flux footprint in time
    fig, ax = plt.subplots(3, 2, figsize=(7.0, 6.4), sharex=True, constrained_layout=True)
    for k, pn in enumerate(plates):
        for r in runs:
            tag = "sh" if r is sh else "ref"
            fits = [eich_fit(r.plates[pn]["s"], r.avg("heatF_total", pn, i)) for i in range(len(r.steps))]
            par = np.array([f if f else (np.nan,) * 5 for f in fits])
            fine = np.linspace(r.plates[pn]["s"][0], r.plates[pn]["s"][-1], 400)
            qpk = np.array([eich(fine, *f).max() if f else np.nan for f in fits]) * 1e-6
            lam_int = r.integral("heatF_total", pn) * 1e-6 / (2 * np.pi * N[f"Rstrike_{pn}"]) / qpk * 1e3   # integral width [mm]
            ax[0, k].plot(r.t, qpk, color=r.color, label=r.label)
            ax[1, k].plot(r.t, lam_int, color=r.color)
            ax[2, k].plot(r.t, par[:, 2] * 1e3, color=r.color)
            for a in ax[:, k]:
                a.axvline(t_on[tag], color=r.color, lw=0.7, ls=":")
            N[f"lamint_pre_{tag}_{pn}"] = lam_int[i_pre[r]]
            N[f"lamint_end_{tag}_{pn}"] = lam_int[-1]
            N[f"qpk_pre_{tag}_{pn}"], N[f"qpk_end_{tag}_{pn}"] = qpk[i_pre[r]], qpk[-1]
            N[f"S_pre_{tag}_{pn}"], N[f"S_end_{tag}_{pn}"] = par[i_pre[r], 2] * 1e3, par[-1, 2] * 1e3
            N[f"P_end_{tag}_{pn}"] = r.integral("heatF_total", pn)[-1] * 1e-6
            i1 = r.index_at(t_on[tag])
            N[f"P_on_{tag}_{pn}"] = r.integral("heatF_total", pn)[i1] * 1e-6
        ax[0, k].set_title(f"{pn} plate")
        ax[2, k].set_xlabel("time [ms]")
    ax[0, 0].set_ylabel("peak heat flux (fit) [MW/m$^2$]")
    ax[1, 0].set_ylabel(r"integral width $\lambda_\mathrm{int}^\mathrm{t}$ [mm]")
    ax[2, 0].set_ylabel(r"Eich spreading $S^\mathrm{t}$ [mm]")
    ax[0, 0].legend(loc="upper left")
    save(fig, "heat_footprint")

    # ------------------------------------------------------------------ numbers for the text
    with open(os.path.join(args.outdir, "numbers.txt"), "w") as g:
        for key in sorted(N):
            g.write(f"{key:28s} {N[key]:.6g}\n")
    print(open(os.path.join(args.outdir, "numbers.txt")).read())
    print(f"figures in {figdir}")


if __name__ == "__main__":
    main()
