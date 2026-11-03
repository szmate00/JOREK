#!/usr/bin/env python3
"""Plate currents, potentials and heat loads from jorek2_postproc `boundary_quantities` output.

Expects SI units and the expressions of doc/sheath_j.md:
  R Z Psi_N bnd_dl Bnorm Jn_par Jn_sat Jn_sheath Jnorm dPhi_float Phi T_e T_i ne heatF_total
(any order, extra columns are ignored), one or several toroidal angles per time step covering one toroidal period.

  python3 util/plot_wall_current.py boundary_quantities_s00700..01700.dat [--n-period 6]

The target plates are the boundary points where the field crosses the wall (|Bnorm| above a fraction of its maximum);
they are split into inner and outer plate at the largest gap in R. Currents are positive from the plasma into the wall.
Figures (plt.show()): plate currents in time, current profiles, potential/temperature/heat-flux profiles (heat flux
as dots with its Eich fit as the line), toroidal map, and the Eich fit parameters in time. The Eich widths are in
target coordinates (lambda_q and S times the flux expansion) unless --flux-expansion is given.
"""
import argparse
import math
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

MU0 = 4.0e-7 * np.pi
AMU = 1.66053906660e-27
trapezoid = getattr(np, "trapezoid", None) or np.trapz

# categorical slots (fixed order), neutral inks, diverging pair with a neutral midpoint
C_INNER, C_OUTER, C_THIRD = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, SURFACE = "#0b0b0b", "#52514e", "#898781", "#fcfcfb"
DIVERGING = LinearSegmentedColormap.from_list("blue_gray_orange", ["#104281", C_INNER, "#f0efec", C_OUTER, "#8f3410"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK2, "ytick.labelcolor": INK2, "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0, "font.size": 10,
    "axes.titlesize": 10, "axes.titleweight": "bold", "legend.frameon": False,
})



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

def read_boundary_quantities(path):
    """Return column names, step numbers, times (JOREK units) and a list of (n_rows, n_col) arrays, one per step."""
    names, steps, times, blocks, rows = None, [], [], [], []

    def flush():
        if rows:
            blocks.append(np.array(rows, dtype=float))
            rows.clear()

    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                if "time step" in s:
                    flush()
                    steps.append(int(s.split("#")[2].split(",")[0]))
                    times.append(float(s.split("=")[-1]))
                elif names is None:
                    names = s.lstrip("#").split()
                continue
            rows.append(s.split())
    flush()
    if names is None or not blocks:
        sys.exit(f"{path}: no header or no data found")
    if len(blocks) != len(steps):
        sys.exit(f"{path}: {len(steps)} time-step headers but {len(blocks)} data blocks")
    return names, np.array(steps), np.array(times), blocks


def split_toroidal(block):
    """(n_rows, n_col) -> (n_phi, n_bnd, n_col); the boundary points repeat for every toroidal angle."""
    same = np.where(np.all(np.isclose(block[1:, :2], block[0, :2], rtol=0, atol=1e-12), axis=1))[0]
    n_bnd = same[0] + 1 if same.size else block.shape[0]
    if block.shape[0] % n_bnd:
        sys.exit("boundary points do not repeat regularly over the toroidal angles")
    return block.reshape(block.shape[0] // n_bnd, n_bnd, block.shape[1])


class Plate:
    """One target plate: points ordered along the plate, s = 0 at the strike point, s > 0 in the SOL.

    Psi_N cannot order the points: JOREK mirrors it in the private region (Psi_N > 1 on both sides of the strike
    point), so the strike point is its minimum and the order is geometric (projection on the plate's direction).
    The SOL side is the one facing away from the other plate.
    """

    def __init__(self, name, color, idx, R, Z, psi, other_RZ):
        x = np.column_stack([R[idx], Z[idx]])
        mean = x.mean(axis=0)
        direction = np.linalg.svd(x - mean)[2][0]
        if np.dot(mean - other_RZ, direction) < 0.0:
            direction = -direction
        order = np.argsort((x - mean) @ direction)
        self.name, self.color, self.idx = name, color, idx[order]
        R, Z = R[self.idx], Z[self.idx]
        s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(R), np.diff(Z)))])
        self.s = s - s[np.argmin(psi[self.idx])]
        self.R = R


_erfc = np.vectorize(math.erfc, otypes=[float])


def eich(s, q0, lam, S, s0, qbg):
    """Eich profile, T. Eich et al., PRL 107 (2011) 215001, in the coordinate of s (lam and S in the same units)."""
    x = s - s0
    return 0.5 * q0 * np.exp(np.clip((0.5 * S / lam) ** 2 - x / lam, -700.0, 700.0)) * _erfc(0.5 * S / lam - x / S) + qbg


def eich_fit(s, q):
    """Least-squares Eich fit (Levenberg-Marquardt, lam and S through their logarithms). Returns
    (q0, lam, S, s0, qbg) or None if there is no usable peak or the fit does not converge."""
    ok = np.isfinite(q)
    s, q = s[ok], q[ok]
    if s.size < 8 or q.max() <= 0.0:
        return None
    model = lambda p: eich(s, p[0], math.exp(p[1]), math.exp(p[2]), p[3], p[4])
    width = max(trapezoid(np.clip(q, 0.0, None), s) / q.max(), 1e-3)
    p = np.array([2.0 * q.max(), math.log(width), math.log(0.3 * width), s[np.argmax(q)] - 0.3 * width, 0.0])
    r = model(p) - q
    cost, damp = r @ r, 1e-2
    for _ in range(200):
        J = np.empty((s.size, p.size))
        for k in range(p.size):
            dp = np.zeros(p.size)
            dp[k] = 1e-6 * max(abs(p[k]), 1e-3)
            J[:, k] = (model(p + dp) - model(p - dp)) / (2.0 * dp[k])
        A, g = J.T @ J, J.T @ r
        try:
            step = np.linalg.solve(A + damp * np.diag(np.diag(A) + 1e-30), -g)
        except np.linalg.LinAlgError:
            return None
        r_new = model(p + step) - q
        if np.all(np.isfinite(r_new)) and r_new @ r_new < cost:
            converged = cost - r_new @ r_new < 1e-12 * cost
            p, r, cost, damp = p + step, r_new, r_new @ r_new, damp * 0.3
            if converged:
                break
        else:
            damp *= 10.0
            if damp > 1e12:
                break
    q0, lam, S = p[0], math.exp(p[1]), math.exp(p[2])
    if not (q0 > 0.0 and lam < 10.0 * np.ptp(s) and S < 10.0 * np.ptp(s)):
        return None
    return q0, lam, S, p[3], p[4]


def phi_mean(a):
    """Toroidal average over the first axis; the angles span one period including both ends."""
    return a[0] if a.shape[0] == 1 else trapezoid(a, axis=0) / (a.shape[0] - 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file")
    ap.add_argument("--n-period", type=int, default=6, help="toroidal periodicity of the run (extent of the map)")
    ap.add_argument("--central-density", type=float, default=1.0, help="in 1e20 m^-3, for the time unit")
    ap.add_argument("--central-mass", type=float, default=2.01410177811)
    ap.add_argument("--bnorm-frac", type=float, default=0.05, help="target = |Bnorm| above this fraction of its maximum")
    ap.add_argument("--flux-expansion", type=float, nargs=2, metavar=("INNER", "OUTER"), default=[1.0, 1.0],
                    help="divide the Eich widths by these to map them to the outer midplane (default: target values)")
    ap.add_argument("--steps", type=int, nargs=3, metavar=("PRE", "PEAK", "LATE"),
                    help="time steps for the profile plots (default: first, maximum outer-plate power, last)")
    args = ap.parse_args()

    names, steps, times, blocks = read_boundary_quantities(args.file)
    col = {n: i for i, n in enumerate(names)}
    missing = [n for n in ("R", "Z", "Psi_N", "bnd_dl", "Bnorm", "Jn_par", "Jn_sat", "Jn_sheath", "dPhi_float", "T_e",
                           "heatF_total") if n not in col]
    if missing:
        sys.exit(f"missing expressions in {args.file}: {' '.join(missing)}")
    t_ms = times * np.sqrt(MU0 * args.central_density * 1e20 * args.central_mass * AMU) * 1e3

    # --- target plates from the first step
    first = split_toroidal(blocks[0])
    n_phi = first.shape[0]
    R, Z, psi = first[0, :, col["R"]], first[0, :, col["Z"]], first[0, :, col["Psi_N"]]
    bn = np.abs(first[:, :, col["Bnorm"]]).mean(axis=0)
    target = np.where(bn > args.bnorm_frac * bn.max())[0]
    Rs = np.sort(R[target])
    R_split = 0.5 * (Rs[np.argmax(np.diff(Rs))] + Rs[np.argmax(np.diff(Rs)) + 1])
    i_in, i_out = target[R[target] < R_split], target[R[target] > R_split]
    centre = lambda i: np.array([R[i].mean(), Z[i].mean()])
    plates = [Plate("inner plate", C_INNER, i_in, R, Z, psi, centre(i_out)),
              Plate("outer plate", C_OUTER, i_out, R, Z, psi, centre(i_in))]
    print(f"{len(steps)} steps ({steps[0]}..{steps[-1]}), {n_phi} toroidal angles, {first.shape[1]} boundary points")
    for p in plates:
        print(f"  {p.name}: {p.idx.size} points, R = {p.R.min():.3f}..{p.R.max():.3f} m, "
              f"private region {-p.s.min():.3f} m, SOL {p.s.max():.3f} m")

    # --- data on the plates only: q[name][plate] has shape (n_step, n_phi, n_point)
    keep = ["Jn_par", "Jn_sat", "Jn_sheath", "dPhi_float", "T_e", "heatF_total"]
    q = {n: [np.empty((len(steps), n_phi, p.idx.size)) for p in plates] for n in keep}
    for i, block in enumerate(blocks):
        b = split_toroidal(block)
        for n in keep:
            for k, p in enumerate(plates):
                q[n][k][i] = b[:, p.idx, col[n]]
    del blocks

    dl = [first[0, p.idx, col["bnd_dl"]] for p in plates]
    ring = [2.0 * np.pi * p.R * d for p, d in zip(plates, dl)]          # area per point, full torus [m^2]

    def integral(name, k):                                               # toroidally averaged plate integral (n_step)
        return np.array([np.sum(phi_mean(a) * ring[k]) for a in q[name][k]])

    I_par = [integral("Jn_par", k) * 1e-3 for k in range(2)]            # kA
    I_sat = [integral("Jn_sat", k) * 1e-3 for k in range(2)]
    P = [integral("heatF_total", k) * 1e-6 for k in range(2)]           # MW
    Te_max = [np.array([phi_mean(a).max() for a in q["T_e"][k]]) for k in range(2)]

    if args.steps:
        sel = [int(np.argmin(np.abs(steps - s))) for s in args.steps]
    else:
        sel = [0, int(np.argmax(P[1])), len(steps) - 1]
    sel_label = [f"step {steps[i]} ({t_ms[i]:.3f} ms)" for i in sel]
    sel_color = [C_THIRD, C_OUTER, C_INNER]

    # --- 1. plate current, target temperature and power against time
    fig, ax = plt.subplots(3, 1, figsize=(8, 8.5), sharex=True, constrained_layout=True)
    for k, p in enumerate(plates):
        ax[0].plot(t_ms, I_par[k], color=p.color, label=p.name)
        ax[0].plot(t_ms, I_sat[k], color=p.color, lw=1.0, ls=(0, (4, 2)), label=f"{p.name}, ion saturation")
        ax[1].plot(t_ms, Te_max[k], color=p.color, label=p.name)
        ax[2].plot(t_ms, P[k], color=p.color, label=p.name)
    ax[0].plot(t_ms, I_par[0] + I_par[1], color=INK2, lw=1.2, label="sum of both plates")
    ax[0].axhline(0.0, color=MUTED, lw=0.8)
    ax[0].set(title="Current from the plasma into the plate (parallel current, toroidal average)", ylabel="kA")
    ax[1].set(title="Peak electron temperature on the plate", ylabel="eV")
    ax[2].set(title="Power to the plate", ylabel="MW", xlabel="time [ms]")
    for a in ax:
        for i in sel:
            a.axvline(t_ms[i], color=MUTED, lw=0.8, ls=":")
        a.legend(loc="best", fontsize=9)

    # --- 2. current profiles along the plates at the three times
    fig, ax = plt.subplots(3, 2, figsize=(10, 9), sharex="col", constrained_layout=True)
    for r, i in enumerate(sel):
        for k, p in enumerate(plates):
            a, s = ax[r, k], p.s * 100.0
            jsat = phi_mean(q["Jn_sat"][k][i]) * 1e-3
            a.fill_between(s, -jsat, jsat, color="#e9e8e4", lw=0, label="within ion saturation")
            a.plot(s, phi_mean(q["Jn_sheath"][k][i]) * 1e-3, color=C_THIRD, lw=1.4, ls=(0, (4, 2)),
                   label="sheath characteristic")
            a.plot(s, phi_mean(q["Jn_par"][k][i]) * 1e-3, color=p.color, label="parallel current")
            a.axhline(0.0, color=MUTED, lw=0.8)
            a.axvline(0.0, color=MUTED, lw=0.8, ls=":")
            a.set_title(f"{p.name}, {sel_label[r]}")
            if k == 0:
                a.set_ylabel("kA/m$^2$ into the plate")
            if r == 2:
                a.set_xlabel("distance from the strike point [cm] (negative: private region)")
            a.legend(loc="best", fontsize=8)

    # --- Eich fits of the toroidally averaged heat flux, every step
    fx = args.flux_expansion
    where = "at the midplane" if args.flux_expansion != [1.0, 1.0] else "on the target"
    fits = [[eich_fit(p.s, phi_mean(a)) for a in q["heatF_total"][k]] for k, p in enumerate(plates)]

    # --- 3. potential above floating, electron temperature and heat flux along the plates
    rows = [("dPhi_float", 1.0, "V", "Potential above floating"), ("T_e", 1.0, "eV", "Electron temperature"),
            ("heatF_total", 1e-6, "MW/m$^2$", "Heat flux and Eich fit")]
    fig, ax = plt.subplots(3, 2, figsize=(10, 9), sharex="col", constrained_layout=True)
    for r, (name, fac, unit, title) in enumerate(rows):
        for k, p in enumerate(plates):
            a = ax[r, k]
            for i, c, lab in zip(sel, sel_color, sel_label):
                if name != "heatF_total":
                    a.plot(p.s * 100.0, phi_mean(q[name][k][i]) * fac, color=c, label=lab)
                    continue
                # heat flux: the data as dots, its Eich fit as the line
                a.plot(p.s * 100.0, phi_mean(q[name][k][i]) * fac, ls="none", marker="o", ms=2.5, color=c,
                       label=None if fits[k][i] else lab)
                if fits[k][i]:
                    fine = np.linspace(p.s[0], p.s[-1], 400)
                    a.plot(fine * 100.0, eich(fine, *fits[k][i]) * fac, color=c, label=f"{lab}: "
                           f"$\\lambda_q$ {fits[k][i][1] / fx[k] * 1e3:.1f} mm, S {fits[k][i][2] / fx[k] * 1e3:.1f} mm")
            if name == "dPhi_float":
                a.axhline(0.0, color=MUTED, lw=0.8)
            a.axvline(0.0, color=MUTED, lw=0.8, ls=":")
            a.set_title(f"{title}, {p.name}")
            if k == 0:
                a.set_ylabel(unit + (f" (widths {where})" if name == "heatF_total" else ""))
            if r == 2:
                a.set_xlabel("distance from the strike point [cm] (negative: private region)")
            a.legend(loc="best", fontsize=8)

    # --- 4. toroidal structure of the current on the plates at the peak
    if n_phi > 1:
        i = sel[1]
        phi_deg = np.linspace(0.0, 360.0 / args.n_period, n_phi)
        fig, ax = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
        vmax = max(np.abs(q["Jn_par"][k][i]).max() for k in range(2)) * 1e-3
        for k, p in enumerate(plates):
            m = ax[k].pcolormesh(p.s * 100.0, phi_deg, q["Jn_par"][k][i] * 1e-3, cmap=DIVERGING, vmin=-vmax, vmax=vmax,
                                 shading="nearest")
            ax[k].grid(False)
            ax[k].set(title=f"{p.name}, {sel_label[1]}", xlabel="distance from the strike point [cm]")
        ax[0].set_ylabel("toroidal angle [deg]")
        fig.colorbar(m, ax=ax, label="parallel current into the plate [kA/m$^2$]")
    else:
        print("one toroidal angle only: no toroidal map")

    # --- 5. Eich fit parameters in time
    par = [np.array([f if f else (np.nan,) * 5 for f in fits[k]]) for k in range(2)]
    fig, ax = plt.subplots(3, 1, figsize=(8, 8.5), sharex=True, constrained_layout=True)
    for k, p in enumerate(plates):
        ax[0].plot(t_ms, par[k][:, 1] / fx[k] * 1e3, color=p.color, label=p.name)
        ax[1].plot(t_ms, par[k][:, 2] / fx[k] * 1e3, color=p.color, label=p.name)
        ax[2].plot(t_ms, par[k][:, 0] * 1e-6, color=p.color, label=p.name)
    ax[0].set(title=f"Eich power fall-off length $\\lambda_q$ ({where})", ylabel="mm")
    ax[1].set(title=f"Eich spreading width S ({where})", ylabel="mm")
    ax[2].set(title="Eich peak parameter $q_0$", ylabel="MW/m$^2$", xlabel="time [ms]")
    for a in ax:
        a.legend(loc="best", fontsize=9)

    # --- numbers
    print(f"\n{'':24s}{'inner plate':>14s}{'outer plate':>14s}")
    for lab, v in (("I_par min [kA]", [x.min() for x in I_par]), ("I_par max [kA]", [x.max() for x in I_par]),
                   ("I_sat at peak [kA]", [x[sel[1]] for x in I_sat]), ("P max [MW]", [x.max() for x in P]),
                   ("energy [kJ]", [trapezoid(x, t_ms) for x in P]), ("Te max [eV]", [x.max() for x in Te_max]),
                   (f"lambda_q, {steps[sel[1]]} [mm]", [par[k][sel[1], 1] / fx[k] * 1e3 for k in range(2)]),
                   (f"S, {steps[sel[1]]} [mm]", [par[k][sel[1], 2] / fx[k] * 1e3 for k in range(2)])):
        print(f"{lab:24s}{v[0]:14.4g}{v[1]:14.4g}")
    show_or_save("wall_current")


if __name__ == "__main__":
    main()
