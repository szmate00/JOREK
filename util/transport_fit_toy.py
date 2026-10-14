#!/usr/bin/env python3
"""
transport_fit_toy.py -- a 1D stand-in for JOREK + jorek2_postproc used to test transport_fit.py
and the driver transport_fit.sh without a cluster.

Steady-state radial transport on psi_N in [0, 1.3] with a hidden D(psi) and chi(psi):
  d/dpsi [ D du/dpsi - v u ] + S(psi) - u/tau(psi) = 0        (density)
  d/dpsi [ chi n dT/dpsi ]   + Q(psi) - n T/tau_T(psi) = 0    (heat, conductivity = chi * n as with
                                                              use_zkperp_times_density)
with a core source, an edge (recycling-like) source, parallel losses for psi_N > 1, an optional
inward pinch v (the "drift" that is not diffusion) and a flat region in the reference.  The run
relaxes towards the steady state with a finite equilibration time so the steadiness gate is
exercised.  Files are written in the jorek2_postproc ASCII format that transport_fit.py reads,
plus empty jorek??????.h5 placeholders so the driver's step bookkeeping works.

Commands (all take --cfg FILE; the run dir is RUN from the cfg):
  init      create the run dir: fake namelist, the hidden-truth reference exports for step 0,
            probe file from the true target profiles (+ noise), constant-coefficient profile files
  run       --nstep N   advance one segment with the current profile files
  export    --step N    write the midplane and boundary exports for a stored step
  selftest  run the whole loop in a temporary directory and check convergence
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import transport_fit as tf  # noqa: E402

PSI_MAX = 1.3
N_CELL = 520
TOY_STATE = "toy_state.json"
TOY_PARAMS = {"pinch_frac": 0.3, "tau_eq": 120.0, "noise": 0.01, "temperatures": 2, "flat_core": True}


# ----------------------------------------------------------------------------- hidden truth
def psi_grid():
    h = PSI_MAX / N_CELL
    return (np.arange(N_CELL) + 0.5) * h, h


def D_true(psi):
    core = 1.0 - 0.8 * np.exp(-((psi - 0.95) / 0.05)**2)          # transport barrier at the pedestal
    sol = 1.0 + 1.5 * np.clip((psi - 1.0) / 0.3, 0.0, 1.0)         # rising into the far SOL
    return np.where(psi < 1.0, core, sol)


def chi_true(psi, which):
    f = 1.0 if which == "e" else 0.7
    return f * (2.0 - 1.6 * np.exp(-((psi - 0.93) / 0.06)**2) + 1.0 * np.clip((psi - 1.0) / 0.3, 0.0, 1.0))


def sources(psi, flat_core):
    S = 2.0 * np.exp(-((psi - 0.3) / 0.25)**2) + 6.0 * np.exp(-((psi - 1.03) / 0.03)**2)
    if flat_core:
        S = np.where(psi < 0.45, 0.0, S)                               # no source inside -> flat core density
    tau = np.where(psi > 1.0, 0.02, np.inf)
    Q = {"e": 3.0 * np.exp(-((psi - 0.4) / 0.3)**2), "i": 2.0 * np.exp(-((psi - 0.4) / 0.3)**2)}
    tau_T = np.where(psi > 1.0, 0.01, np.inf)
    return S, tau, Q, tau_T


def solve_steady(psi, h, coef, src, tau, n=None, v=None, u_edge=0.0):
    """Finite-volume steady state of d/dpsi[c (n) du/dpsi - v u] + src - u/tau = 0, Neumann at 0,
    Dirichlet u_edge at PSI_MAX.  coef is D (n None) or chi (times n)."""
    N = len(psi)
    cf = 0.5 * (coef[1:] + coef[:-1])
    if n is not None:
        cf = cf * 0.5 * (n[1:] + n[:-1])
    A = np.zeros((N, N))
    b = src.copy()
    for i in range(N):
        A[i, i] += 1.0 / tau[i] if np.isfinite(tau[i]) else 0.0
    for f in range(N - 1):                                             # face between f and f+1
        k = cf[f] / h**2
        A[f, f] += k
        A[f, f + 1] -= k
        A[f + 1, f + 1] += k
        A[f + 1, f] -= k
        if v is not None:                                              # upwind pinch, v < 0 = inward
            vf = 0.5 * (v[f] + v[f + 1]) / h
            if vf < 0:                                                 # flux -v u taken from cell f+1
                A[f, f + 1] += vf
                A[f + 1, f + 1] -= vf
            else:
                A[f, f] += vf
                A[f + 1, f] -= vf
    k = coef[-1] * (n[-1] if n is not None else 1.0) / h**2 * 2.0     # Dirichlet ghost at the edge
    A[-1, -1] += k
    b[-1] += k * u_edge
    return np.linalg.solve(A, b)


def steady_profiles(psi, h, D, chi_e, chi_i, params):
    S, tau, Q, tau_T = sources(psi, params["flat_core"])
    v = -params["pinch_frac"] * D * 20.0 * np.exp(-((psi - 0.97) / 0.06)**2)   # inward pinch near the separatrix
    n = solve_steady(psi, h, D, S, tau, v=v, u_edge=0.02)
    Te = solve_steady(psi, h, chi_e, Q["e"], tau_T, n=n, u_edge=0.0005)      # edge below the loss-dominated SOL
    if params["temperatures"] == 2:
        Ti = solve_steady(psi, h, chi_i, Q["i"], tau_T, n=n, u_edge=0.0005)
        return {"ne": n, "T_e": Te, "T_i": Ti}
    return {"ne": n, "T": Te}


# ----------------------------------------------------------------------------- fake geometry and exports
def fake_wall():
    """A closed wall: a circle (Bnorm ~ 0) plus two straight target plates (Bnorm = 1)."""
    th = np.linspace(0, 2 * np.pi, 240, endpoint=False)
    Rc, Zc = 1.6 + 0.6 * np.cos(th), 0.1 + 0.9 * np.sin(th)
    bn_c = 0.01 * np.ones_like(Rc)
    plates = {}
    for name, R0, sign in (("inner", 1.25, -1.0), ("outer", 1.95, 1.0)):
        s = np.linspace(-0.05, 0.12, 70)                               # m along the plate, strike point at 0
        R, Z = R0 + sign * 0.3 * s, -1.0 + 0.95 * s
        plates[name] = {"s": s, "R": R, "Z": Z}
    return (Rc, Zc, bn_c), plates


def plate_psi(s):
    """psi_N along a plate: SOL side compresses the flux (lambda_s = 1.5 cm per 0.03 psi), PFR mirrored."""
    return np.where(s > 0, 1.0 + s / 0.5, 1.0 + np.abs(s) / 0.35)


def target_profiles(prof, psi_s):
    psi, _ = psi_grid()
    out = {}
    for k, u in prof.items():
        out[k] = np.interp(psi_s, psi, u)
    out["ne_t"] = 0.6 * out["ne"]
    out["Te_t"] = 0.35 * out["T_e" if "T_e" in out else "T"]
    return out


def write_export(path, names, cols, step, boundary=False):
    with open(path, "w") as f:
        f.write("# " + "".join(f"{n:<23s}" for n in names) + "\n")
        f.write(f"# time step #{step:06d}" + (",  t_now = 0.0" if boundary else "") + "\n")
        for row in np.column_stack(cols):
            f.write("".join(f"{v:23.15e}" for v in row) + "\n")


def export_step(run, cfg, step, prof):
    psi, _ = psi_grid()
    pp = os.path.join(run, "postproc")
    os.makedirs(pp, exist_ok=True)
    R = 1.65 + 0.5 * psi
    two = "T_e" in prof
    names = tf.midplane_exprs(cfg)
    # exports in SI-like units: densities x 1e20 m^-3, temperatures x 100 eV
    cols = [R, psi, prof["ne"] / 3.0, prof["ne"] * 1e20] + ([prof["T_e"] * 100, prof["T_i"] * 100] if two else [prof["T"] * 100])
    write_export(os.path.join(pp, f"exprs_outer-midplane_s{step:06d}.dat"), names, cols, step)
    (Rc, Zc, bnc), plates = fake_wall()
    Rl, Zl, psil, bnl, nel, tel = [Rc], [Zc], [1.0 + 0.3 * np.ones_like(Rc)], [bnc], [np.zeros_like(Rc)], [np.zeros_like(Rc)]
    for name, p in plates.items():
        ps = plate_psi(p["s"])
        t = target_profiles(prof, ps)
        Rl.append(p["R"])
        Zl.append(p["Z"])
        psil.append(ps)
        bnl.append(np.ones_like(ps))
        nel.append(t["ne_t"])
        tel.append(t["Te_t"])
    names = tf.boundary_exprs(cfg)
    ne_all = np.concatenate(nel) * 1e20
    cols = [np.concatenate(Rl), np.concatenate(Zl), np.concatenate(psil), np.concatenate(bnl), ne_all,
            np.concatenate(tel) * 100, 0.1 * ne_all, ne_all]
    write_export(os.path.join(pp, f"boundary_quantities_s{step:06d}.dat"), names, cols, step, boundary=True)


# ----------------------------------------------------------------------------- run dir
def write_namelist(path, two):
    zk = ("  ZK_e_perp_file = 'ZK_e_toy.dat'\n  ZK_i_perp_file = 'ZK_i_toy.dat'\n" if two
          else "  ZK_perp_file = 'ZK_toy.dat'\n")
    with open(path, "w") as f:
        f.write(" toy namelist for transport_fit_toy.py\n &in1\n  restart = .f.\n  nstep_n = 0\n  tstep_n = 1.\n"
                "  nout = 1\n  central_density = 1.0\n  central_mass = 2.0141017781\n"
                "  use_zkperp_times_density = .true.\n  D_perp_file = 'D_toy.dat'\n" + zk + " /\n")


def coeff_files(run, cfg, two):
    nml = tf.parse_namelist(os.path.join(run, cfg.s("NAMELIST")))
    nf = tf.norm_factors(nml)
    psi, _ = psi_grid()
    D = tf.read_prof(os.path.join(run, nml["d_perp_file"]))
    D = tf.pchip(D[0], D[1], psi) * nf["D_to_SI"]
    if two:
        e = tf.read_prof(os.path.join(run, nml["zk_e_perp_file"]))
        i = tf.read_prof(os.path.join(run, nml["zk_i_perp_file"]))
        return D, tf.pchip(e[0], e[1], psi) * nf["ZK_to_SI"], tf.pchip(i[0], i[1], psi) * nf["ZK_to_SI"]
    z = tf.read_prof(os.path.join(run, nml["zk_perp_file"]))
    return D, tf.pchip(z[0], z[1], psi) * nf["ZK_to_SI"], None


def cmd_init(cfg, args):
    run = cfg.s("RUN")
    os.makedirs(run, exist_ok=True)
    params = dict(TOY_PARAMS)
    params["temperatures"] = cfg.i("TEMPERATURES")
    params["pinch_frac"] = args.pinch if args.pinch is not None else params["pinch_frac"]
    params["tau_eq"] = args.tau_eq if args.tau_eq is not None else params["tau_eq"]
    two = params["temperatures"] == 2
    write_namelist(os.path.join(run, cfg.s("NAMELIST")), two)
    nml = tf.parse_namelist(os.path.join(run, cfg.s("NAMELIST")))
    nf = tf.norm_factors(nml)
    psi, h = psi_grid()
    truth = steady_profiles(psi, h, D_true(psi), chi_true(psi, "e"), chi_true(psi, "i"), params)
    # constant first guess, in JOREK units
    x = np.linspace(0.0, 1.4, 200)
    tf.write_prof(os.path.join(run, "D_toy.dat"), x, np.full(200, args.d0 / nf["D_to_SI"]))
    for name in (("ZK_e_toy.dat", "ZK_i_toy.dat") if two else ("ZK_toy.dat",)):
        tf.write_prof(os.path.join(run, name), x, np.full(200, args.chi0 / nf["ZK_to_SI"]))
    rng = np.random.default_rng(1)
    # probes on both plates from the true target profiles, with noise
    _, plates = fake_wall()
    with open(os.path.join(run, "probes.dat"), "w") as f:
        f.write("# R Z ne[1e19] Te[eV]  (toy probes from the hidden truth, 8% noise)\n")
        for name, p in plates.items():
            sel = np.arange(4, len(p["s"]), 5)
            t = target_profiles(truth, plate_psi(p["s"][sel]))
            for j, k in enumerate(sel):
                ne = t["ne_t"][j] * (1 + 0.08 * rng.standard_normal()) * 10.0      # "1e19" units of a 1e20 value
                te = t["Te_t"][j] * (1 + 0.08 * rng.standard_normal()) * 100.0
                f.write(f"{p['R'][k]:.5f} {p['Z'][k]:.5f} {ne:.4f} {te:.3f}\n")
    state = {"step": 0, "params": params, "prof": {k: v.tolist() for k, v in truth.items()},
             "truth": {k: v.tolist() for k, v in truth.items()}}
    with open(os.path.join(run, TOY_STATE), "w") as f:
        json.dump(state, f)
    export_step(run, cfg, 0, truth)
    open(os.path.join(run, "jorek000000.h5"), "w").close()
    open(os.path.join(run, "jorek_restart.h5"), "w").close()
    print(f"toy init: run dir {run}, {'two' if two else 'one'} temperature(s), pinch {params['pinch_frac']}, "
          f"tau_eq {params['tau_eq']} steps, first guess D {args.d0} chi {args.chi0} m2/s")


def cmd_run(cfg, args):
    run = cfg.s("RUN")
    with open(os.path.join(run, TOY_STATE)) as f:
        st = json.load(f)
    params = st["params"]
    two = params["temperatures"] == 2
    if os.path.exists(os.path.join(run, "toy_fail_once")):              # test hook for the driver's crash path
        os.remove(os.path.join(run, "toy_fail_once"))
        raise SystemExit("toy run: forced failure (toy_fail_once)")
    nml = tf.parse_namelist(os.path.join(run, cfg.s("NAMELIST")))
    if not bool(tf.nml_get(nml, "restart", False)):
        raise SystemExit("toy run: restart is not .t. in the namelist (the driver should set it)")
    nstep = int(tf.nml_get(nml, "nstep_n", 0)) if args.nstep is None else args.nstep
    D, chi_e, chi_i = coeff_files(run, cfg, two)
    if np.any(D <= 0) or np.any(chi_e <= 0):
        raise SystemExit("toy run: non-positive coefficient -> 'crash'")
    psi, h = psi_grid()
    target = steady_profiles(psi, h, D, chi_e, chi_i if two else chi_e, params)
    nout = max(1, int(tf.nml_get(nml, "nout", 1)))
    first = st["step"]
    saved = st.setdefault("saved", {})
    for k_step in range(nout, nstep + 1, nout):                         # exponential approach, stored every nout
        alpha = 1.0 - np.exp(-k_step / params["tau_eq"])
        saved[str(first + k_step)] = {k: (np.array(u0) + alpha * (target[k] - np.array(u0))).tolist()
                                      for k, u0 in st["prof"].items()}
    st["prof"] = saved[str(first + nstep)]
    st["step"] = first + nstep
    for key in [k for k in saved if int(k) < st["step"] - 2 * nstep]:  # keep the recent segments only
        del saved[key]
    with open(os.path.join(run, TOY_STATE), "w") as f:
        json.dump(st, f)
    open(os.path.join(run, f"jorek{st['step']:06d}.h5"), "w").close()
    print(f"toy run: {nstep} steps -> step {st['step']} (relaxation {alpha:.2f} of the way to steady state)")


def cmd_export(cfg, args):
    run = cfg.s("RUN")
    with open(os.path.join(run, TOY_STATE)) as f:
        st = json.load(f)
    if args.step == 0:
        prof = {k: np.array(v) for k, v in st["truth"].items()}
    else:
        saved = st.get("saved", {})
        if str(args.step) not in saved:
            raise SystemExit(f"toy export: step {args.step} not stored")
        prof = {k: np.array(v) for k, v in saved[str(args.step)].items()}
        rng = np.random.default_rng(args.step + 7)                     # measurement noise, fixed per step
        psi, _ = psi_grid()
        prof = {k: v * (1.0 + st["params"]["noise"] * rng.standard_normal(len(v)) * np.exp(-((psi - 1.0) / 0.4)**2))
                for k, v in prof.items()}
    export_step(run, cfg, args.step, prof)


# ----------------------------------------------------------------------------- self test
def cmd_selftest(cfg_unused, args):
    here = os.path.dirname(os.path.abspath(__file__))
    tmp = tempfile.mkdtemp(prefix="tfit_toy_")
    run = os.path.join(tmp, "run")
    cfgp = os.path.join(tmp, "tfit.cfg")
    two = args.temperatures
    with open(cfgp, "w") as f:
        f.write(f"RUN={run}\nNAMELIST=input\nTEMPERATURES={two}\nPROBE_FILE=probes.dat\nREF_STEP=0\n"
                f"NSTEP_SEG={args.nstep or 200}\nNSTEP_MAX=1600\nMAX_ITER={args.max_iter}\nTOL=0.03\nRELAX={args.relax}\n"
                f"CLAMP=2.0\nPSI_MIN=0.7\nW_TARGET={args.w_target}\nLINEPOINTS=0\nSMOOTH=1\n")
    cfg = tf.Cfg(cfgp)
    py = [sys.executable, os.path.join(here, "transport_fit.py")]
    toy = [sys.executable, os.path.abspath(__file__)]
    a = argparse.Namespace(d0=args.d0, chi0=args.chi0, pinch=args.pinch, tau_eq=args.tau_eq)
    cmd_init(cfg, a)

    def sh(cmd, ok=(0,)):
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode not in ok:
            print(r.stdout)
            print(r.stderr)
            raise SystemExit(f"command failed ({r.returncode}): {' '.join(cmd)}")
        return r

    print(sh(py + ["check", "--cfg", cfgp]).stdout)
    print(sh(py + ["reference", "--cfg", cfgp, "--step", "0"]).stdout)
    it = 0
    while True:
        status = dict(kv.split("=") for kv in sh(py + ["status", "--cfg", cfgp]).stdout.split())
        if status["DONE"] == "1":
            break
        nstep = int(status["NSTEP"])
        steps = sorted(int(os.path.basename(p)[5:11]) for p in __import__("glob").glob(os.path.join(run, "jorek[0-9]*.h5")))
        first = steps[-1]
        sh(py + ["set-namelist", "--cfg", cfgp, "--set", "restart=.t.", "--set", f"nstep_n={nstep}"])
        sh(toy + ["run", "--cfg", cfgp, "--nstep", str(nstep)])
        last = first + nstep
        mid = first + nstep // 2
        sh(toy + ["export", "--cfg", cfgp, "--step", str(mid)])
        sh(toy + ["export", "--cfg", cfgp, "--step", str(last)])
        r = sh(py + ["update", "--cfg", cfgp, "--first", str(first), "--mid", str(mid), "--last", str(last)],
               ok=(tf.EXIT_CONTINUE, tf.EXIT_CONVERGED, tf.EXIT_RERUN, tf.EXIT_STOP))
        print(r.stdout)
        it += 1
        if r.returncode in (tf.EXIT_CONVERGED, tf.EXIT_STOP):
            break
    print(sh(py + ["report", "--cfg", cfgp]).stdout)
    # --- verdict: recovered coefficients vs the hidden truth over the fitted, non-flat range
    nml = tf.parse_namelist(os.path.join(run, "input"))
    nf = tf.norm_factors(nml)
    st = tf.load_state(cfg)
    knots = np.array(st["knots"])
    worst = 0.0
    for name, c in st["channels"].items():
        x, y = tf.channel_file_values(st, c)
        got = tf.pchip(x, y, knots) * (nf["D_to_SI"] if c["kind"] == "D" else nf["ZK_to_SI"])
        true = D_true(knots) if c["kind"] == "D" else chi_true(knots, "i" if name == "T_i" else "e")
        sel = np.array([fl in ("", "wrong-sign") for fl in c["flags"]]) & (knots >= 0.75) & (knots <= 1.2)
        err = np.abs(np.log(got[sel] / true[sel]))
        print(f"  {name:4s}: coefficient recovered within a factor {np.exp(err.max()):.2f} (max), "
              f"{np.exp(np.median(err)):.2f} (median) over the {sel.sum()} fitted knots in psi_N 0.75..1.2")
        worst = max(worst, float(np.median(err)))
    last_row = st["history"][-1]
    mism = max(v for k, v in last_row.items() if k.startswith("mism_"))
    print(f"  final profile mismatch {mism * 100:.1f}% after {st['iter']} iterations; run dir {run}")
    ok = mism < 0.05 and worst < np.log(1.3)
    print("SELFTEST " + ("PASSED" if ok else "FAILED"))
    if not args.keep:
        shutil.rmtree(tmp, ignore_errors=True)
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["init", "run", "export", "selftest"])
    ap.add_argument("--cfg")
    ap.add_argument("--nstep", type=int)
    ap.add_argument("--step", type=int)
    ap.add_argument("--d0", type=float, default=0.5, help="first-guess D [m2/s]")
    ap.add_argument("--chi0", type=float, default=1.0, help="first-guess chi [m2/s]")
    ap.add_argument("--pinch", type=float, default=None, help="pinch strength (0 = pure diffusion)")
    ap.add_argument("--tau-eq", type=float, default=None, help="equilibration time in steps")
    ap.add_argument("--temperatures", type=int, default=2)
    ap.add_argument("--max-iter", type=int, default=15)
    ap.add_argument("--relax", type=float, default=0.7)
    ap.add_argument("--w-target", type=float, default=0.5)
    ap.add_argument("--keep", action="store_true", help="keep the temporary run dir of the self test")
    args = ap.parse_args()
    if args.command == "selftest":
        cmd_selftest(None, args)
        return
    if not args.cfg:
        raise SystemExit("--cfg is required")
    cfg = tf.Cfg(args.cfg)
    {"init": cmd_init, "run": cmd_run, "export": cmd_export}[args.command](cfg, args)


if __name__ == "__main__":
    main()
