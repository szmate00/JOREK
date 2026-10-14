#!/usr/bin/env python3
"""
transport_fit.py -- automated fitting of the tabulated D_perp / ZK_perp profiles of a
JOREK model600 run so that the simulated upstream profiles hold the step-0 profiles and
the target profiles follow Langmuir-probe data.

Driven by util/transport_fit.sh (settings at the top of that script are written into a
small KEY=VALUE file, tfit.cfg, that every command here reads).  Only numpy is required;
matplotlib is used for the figures when it is importable.

Commands (all take --cfg FILE):
  check                       pre-flight checks, state initialisation, profile files written
  reference  --step N         register the exported reference (upstream) step and the plates
  status                      print ITER=.. NSTEP=.. DONE=.. for the bash driver
  update     --first N --mid M --last L
                              one iteration: steadiness gate (mid vs last), gradient ratios, new files
                              exit code 0 continue, 10 converged, 20 rerun longer, 30 stop
  crash                       segment failed: restore files, halve the clamp
  report                      summary of tfit_log.csv

Method: per knot in psi_N, ln D <- ln D + step, step = -relax * ln(g) / s, where
g = |du/dpsi|_sim / |du/dpsi|_ref (density weighted for the heat channels when
use_zkperp_times_density is on) and s is the measured sensitivity d ln|du/dpsi| / d ln D
(-1 for pure diffusion).  The written file is the user's original table times
exp(PCHIP(ln factor at knots)).  See doc/transport_fit/PLAN.md.
"""
import argparse
import glob
import json
import math
import os
import re
import shutil
import sys

import numpy as np

MU0 = 4.0e-7 * math.pi
M_AMU = 1.66053906660e-27
GAMMA = 5.0 / 3.0
EXIT_CONTINUE, EXIT_CONVERGED, EXIT_RERUN, EXIT_STOP = 0, 10, 20, 30

STATE_FILE = "tfit_state.json"
LOG_FILE = "tfit_log.csv"
REF_MID = "tfit_ref_midplane.dat"
REF_BND = "tfit_ref_boundary.dat"

CFG_DEFAULTS = {
    "RUN": ".", "NAMELIST": "input", "TEMPERATURES": "2", "PROBE_FILE": "",
    "REF_STEP": "0", "NSTEP_SEG": "200", "NSTEP_MAX": "800", "MAX_ITER": "12",
    "TOL": "0.05", "RELAX": "0.7", "CLAMP": "2.0", "PSI_MIN": "0.7", "W_TARGET": "0.5",
    "TE_TARGET_WEIGHT": "0.5", "CHI_EDGE_MAX": "10.0", "D_MIN": "0.01", "D_MAX": "50.0",
    "CHI_MIN": "0.01", "CHI_MAX": "100.0", "STEADY_RATIO": "0.5", "FLAT_REL": "0.3",
    "LINEPOINTS": "2000", "DENSITY_EXPR": "ne", "KNOTS": "", "PROBE_MAX_DIST": "0.02",
    "PLATE_ZMAX": "", "PROBE_NE_UNIT": "1e19", "PROBE_ERR_REL": "0.15",
    "PROBE_NE_ERR_FLOOR": "0.1", "PROBE_TE_ERR_FLOOR": "2.0", "SMOOTH": "1",
    "FIT_HEAT_FROM_ITER": "3", "EDGE_HOLD": "0.08", "SENS_MIN_STEP": "0.15", "SMALL_FRAC": "0.03",
    "TSTEP": "",
}


# ----------------------------------------------------------------------------- config
class Cfg:
    def __init__(self, path):
        d = dict(CFG_DEFAULTS)
        with open(path) as f:
            for line in f:
                line = line.split("#", 1)[0].strip()
                if not line or "=" not in line:
                    continue
                k, v = line.split("=", 1)
                d[k.strip().upper()] = v.strip().strip('"').strip("'")
        self.d = d
        self.cfg_path = path

    def s(self, k):
        return self.d[k]

    def f(self, k):
        return float(self.d[k])

    def i(self, k):
        return int(float(self.d[k]))

    def b(self, k):
        return self.d[k].lower() in ("1", "t", "true", "yes", ".t.", ".true.")

    def path(self, k):
        p = self.d[k]
        if not p:
            return ""
        return p if os.path.isabs(p) else os.path.join(self.d["RUN"], p)


def say(msg=""):
    print(msg, flush=True)


def fail(msg):
    say("ERROR: " + msg)
    sys.exit(2)


# ----------------------------------------------------------------------------- namelist
_ASSIGN = re.compile(r"([A-Za-z_][A-Za-z0-9_%]*)\s*(\(\s*([0-9:,\s]+)\s*\))?\s*=")


def _strip_comment(line):
    out, q = [], None
    for ch in line:
        if q:
            out.append(ch)
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
            out.append(ch)
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out)


def _parse_value(tok):
    t = tok.strip()
    if not t:
        return None
    if t[0] in "'\"":
        j = t.find(t[0], 1)
        return t[1:j] if j > 0 else t.strip("'\"")
    tl = t.lower()
    if tl in (".t.", ".true.", "t", "true"):
        return True
    if tl in (".f.", ".false.", "f", "false"):
        return False
    try:
        return float(tl.replace("d", "e"))
    except ValueError:
        return t


def _split_values(text):
    vals, buf, q = [], [], None
    for ch in text:
        if q:
            buf.append(ch)
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
            buf.append(ch)
        elif ch == ",":
            vals.append("".join(buf))
            buf = []
        elif ch in "/&":                      # end of the namelist group
            break
        else:
            buf.append(ch)
    vals.append("".join(buf))
    return [v for v in (x.strip() for x in vals) if v != ""]


def parse_namelist(path):
    """Tolerant Fortran-namelist reader.  Scalars -> value; arrays -> {index: value} (1-based)."""
    nml = {}
    with open(path, errors="replace") as f:
        text = "\n".join(_strip_comment(l) for l in f)
    matches = list(_ASSIGN.finditer(text))
    for i, m in enumerate(matches):
        key = m.group(1).lower()
        idx = m.group(3)
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        vals = [_parse_value(v) for v in _split_values(text[m.end():end])]
        vals = [v for v in vals if v is not None]
        if not vals:
            continue
        if idx is None and len(vals) == 1:
            if key in nml and isinstance(nml[key], dict):
                nml[key][1] = vals[0]
            else:
                nml[key] = vals[0]
            continue
        arr = nml.get(key)
        if not isinstance(arr, dict):
            arr = {} if arr is None else {1: arr}
        if idx is None:
            start = 1
        elif ":" in idx:
            start = int(idx.split(":")[0])
        else:
            start = int(idx.split(",")[0])
        for j, v in enumerate(vals):
            arr[start + j] = v
        nml[key] = arr
    return nml


def nml_get(nml, key, default=None):
    v = nml.get(key.lower(), default)
    if isinstance(v, dict):
        return v.get(1, default)
    return v


def nml_array(nml, key, n, default=0.0):
    v = nml.get(key.lower())
    out = np.full(n, default, dtype=float)
    if isinstance(v, dict):
        for i, x in v.items():
            if 1 <= i <= n and isinstance(x, (int, float)):
                out[i - 1] = x
    elif isinstance(v, (int, float)):
        out[0] = v
    return out


def set_namelist(path, assignments):
    """Replace the value of each key (first and only occurrence) or insert after &in1."""
    with open(path, errors="replace") as f:
        lines = f.readlines()
    for key, val in assignments.items():
        # value = everything up to a comment, the end of the line, or the next "key =" on the same line
        pat = re.compile(r"^(\s*" + re.escape(key) + r"\s*(\(\s*1\s*\))?\s*=\s*)([^!\n]*?)(?=\s*,\s*[A-Za-z_]\w*\s*(\([^)]*\))?\s*=|\s*!|\s*\n|$)", re.I)
        hits = [i for i, l in enumerate(lines) if pat.match(l)]
        if len(hits) > 1:
            fail(f"{path}: '{key}' appears {len(hits)} times; make it unique before running the loop")
        if hits:
            i = hits[0]
            m = pat.match(lines[i])
            rest = lines[i][m.end():]
            lines[i] = f"{m.group(1)}{val}{'  ' + rest if rest.strip() else chr(10)}"
            if not lines[i].endswith("\n"):
                lines[i] += "\n"
        else:
            j = next((i for i, l in enumerate(lines) if l.strip().lower().startswith("&in1")), None)
            if j is None:
                fail(f"{path}: no &in1 group found, cannot insert '{key}'")
            lines.insert(j + 1, f"  {key} = {val}\n")
    with open(path, "w") as f:
        f.writelines(lines)


def fortran_logical(b):
    return ".t." if b else ".f."


# ----------------------------------------------------------------------------- units
def norm_factors(nml):
    n0 = float(nml_get(nml, "central_density", 1.0)) * 1e20
    mass = float(nml_get(nml, "central_mass", 2.01410177811))
    rho0 = mass * M_AMU * n0
    s = math.sqrt(MU0 * rho0)
    return {"n0": n0, "rho0": rho0, "sqrt_mu0rho0": s,
            "D_to_SI": 1.0 / s,                       # D_SI = D_jorek / sqrt(mu0 rho0)
            "ZK_to_SI": 1.0 / ((GAMMA - 1.0) * s)}    # chi_SI = ZK_file / ((gamma-1) sqrt(mu0 rho0))  [times density]


# ----------------------------------------------------------------------------- profiles
def read_prof(path):
    x, y = [], []
    with open(path) as f:
        for line in f:
            t = line.split()
            if len(t) < 2:
                break
            try:
                x.append(float(t[0].replace("d", "e").replace("D", "e")))
                y.append(float(t[1].replace("d", "e").replace("D", "e")))
            except ValueError:
                break
    if len(x) < 2:
        fail(f"{path}: fewer than two numeric rows (readProf format: two columns, no header)")
    x, y = np.array(x), np.array(y)
    if np.any(np.diff(x) <= 0):
        fail(f"{path}: psi_N must be strictly increasing")
    return x, y


def write_prof(path, x, y):
    with open(path, "w") as f:
        for a, b in zip(x, y):
            f.write(f"  {a:.10e}   {b:.10e}\n")


def tanh_profile(c, psi, model_ge_300=True):
    """get_dperp1 / get_zkperp1 analytic form; c = 6 namelist coefficients."""
    c = np.asarray(c, dtype=float)
    c4 = c[3] if c[3] != 0.0 else 1e-12
    y = c[0] * ((1.0 - c[1]) + c[1] * (0.5 - 0.5 * np.tanh((psi - c[4]) / c4)))
    if model_ge_300:
        y = y + c[5] * c[1] * (0.5 - 0.5 * np.tanh((-psi + c[4] + c[2]) / c4))
    return y


def pchip(x, y, xq):
    """Monotone cubic Hermite interpolation (Fritsch-Carlson), constant outside [x0, xn]."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    if n == 1:
        return np.full_like(np.asarray(xq, float), y[0])
    h = np.diff(x)
    d = np.diff(y) / h
    m = np.zeros(n)
    if n == 2:
        m[:] = d[0]
    else:
        for k in range(1, n - 1):
            if d[k - 1] * d[k] <= 0.0:
                m[k] = 0.0
            else:
                w1, w2 = 2.0 * h[k] + h[k - 1], h[k] + 2.0 * h[k - 1]
                m[k] = (w1 + w2) / (w1 / d[k - 1] + w2 / d[k])
        m[0] = ((2.0 * h[0] + h[1]) * d[0] - h[0] * d[1]) / (h[0] + h[1])
        if np.sign(m[0]) != np.sign(d[0]):
            m[0] = 0.0
        elif np.sign(d[0]) != np.sign(d[1]) and abs(m[0]) > 3.0 * abs(d[0]):
            m[0] = 3.0 * d[0]
        m[-1] = ((2.0 * h[-1] + h[-2]) * d[-1] - h[-1] * d[-2]) / (h[-1] + h[-2])
        if np.sign(m[-1]) != np.sign(d[-1]):
            m[-1] = 0.0
        elif np.sign(d[-1]) != np.sign(d[-2]) and abs(m[-1]) > 3.0 * abs(d[-1]):
            m[-1] = 3.0 * d[-1]
    xq = np.asarray(xq, float)
    i = np.clip(np.searchsorted(x, xq, side="right") - 1, 0, n - 2)
    t = (xq - x[i]) / h[i]
    t = np.clip(t, 0.0, 1.0)
    h00 = 2 * t**3 - 3 * t**2 + 1
    h10 = t**3 - 2 * t**2 + t
    h01 = -2 * t**3 + 3 * t**2
    h11 = t**3 - t**2
    return h00 * y[i] + h10 * h[i] * m[i] + h01 * y[i + 1] + h11 * h[i] * m[i + 1]


# ----------------------------------------------------------------------------- postproc exports
def read_export(path):
    """jorek2_postproc ASCII: '# ' + 23-char padded names, '# time step' comments, blocks of rows."""
    names, blocks, buf = None, [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if names is None and "time step" not in line:
                    body = line[2:].rstrip("\n")
                    names = [body[i:i + 23].strip() for i in range(0, len(body), 23)]
                    names = [n for n in names if n]
                elif "time step" in line and buf:
                    blocks.append(buf)
                    buf = []
            elif line.strip():
                buf.append(line)
            elif buf:
                blocks.append(buf)
                buf = []
    if buf:
        blocks.append(buf)
    if names is None or not blocks:
        fail(f"{path}: no header / no data")
    arrs = []
    for b in blocks:
        a = np.array(" ".join(b).split(), dtype=float)
        if a.size % len(names):
            fail(f"{path}: row length does not match the {len(names)} header names")
        arrs.append(a.reshape(-1, len(names)))
    return names, arrs


def read_midplane(path):
    names, arrs = read_export(path)
    a = arrs[-1]
    d = {n: a[:, i] for i, n in enumerate(names)}
    if "Psi_N" not in d:
        fail(f"{path}: Psi_N column missing")
    ok = np.all(np.isfinite(a), axis=1)
    order = np.argsort(d["Psi_N"][ok])
    return {n: v[ok][order] for n, v in d.items()}


def read_boundary(path):
    names, arrs = read_export(path)
    a = arrs[0]
    d = {n: a[:, i] for i, n in enumerate(names)}
    for k in ("R", "Z", "Psi_N"):
        if k not in d:
            fail(f"{path}: {k} column missing")
    # one toroidal plane only (boundary_quantities without arguments); drop a repeated closing point
    R, Z = d["R"], d["Z"]
    n = len(R)
    for i in range(4, n):
        if abs(R[i] - R[0]) < 1e-12 and abs(Z[i] - Z[0]) < 1e-12 and n % i == 0:
            n = i
            break
    d = {k: v[:n] for k, v in d.items()}
    ok = np.all(np.isfinite(np.column_stack(list(d.values()))), axis=1)
    return {k: v[ok] for k, v in d.items()}


# ----------------------------------------------------------------------------- knots and gradients
def default_knots(psi_max):
    k = [0.0, 0.3, 0.5, 0.6, 0.7, 0.75, 0.8, 0.84, 0.87, 0.90, 0.92, 0.94, 0.96, 0.98, 1.00,
         1.02, 1.04, 1.06, 1.08, 1.10, 1.13, 1.16, 1.20, 1.25, 1.30]
    x = 1.4
    while x < psi_max + 0.1:
        k.append(round(x, 3))
        x += 0.1
    k = np.array([v for v in k if v <= psi_max + 0.1])
    if k[-1] < psi_max:
        k = np.append(k, psi_max)
    return k


def local_linear(x, u, x0, h):
    """Tricube-weighted local linear fit at x0 with bandwidth h: (value, slope, n_points)."""
    for _ in range(5):
        d = np.abs(x - x0) / h
        w = np.where(d < 1.0, (1.0 - d**3)**3, 0.0)
        n = int(np.count_nonzero(w > 0))
        if n >= 4:
            break
        h *= 2.0
    if n < 3:
        return np.nan, np.nan, n
    X = np.column_stack([np.ones_like(x), x - x0])
    W = w[:, None]
    A = X.T @ (W * X)
    b = X.T @ (w * u)
    try:
        a0, a1 = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, n
    return a0, a1, n


def knot_bandwidths(knots):
    h = np.empty_like(knots)
    for i in range(len(knots)):
        lo = knots[i] - knots[i - 1] if i > 0 else knots[i + 1] - knots[i]
        hi = knots[i + 1] - knots[i] if i < len(knots) - 1 else knots[i] - knots[i - 1]
        h[i] = 0.6 * max(lo, hi)
    return h


def profile_at_knots(mid, expr, knots):
    psi, u = mid["Psi_N"], mid[expr]
    h = knot_bandwidths(knots)
    val, slope = np.full(len(knots), np.nan), np.full(len(knots), np.nan)
    for i, k in enumerate(knots):
        if k < psi.min() - h[i] or k > psi.max() + h[i]:
            continue
        val[i], slope[i], _ = local_linear(psi, u, k, h[i])
    return val, slope


# ----------------------------------------------------------------------------- plates and probes
def read_probes(path, cfg):
    """Columns: R Z ne Te [ne_err Te_err]; ne in PROBE_NE_UNIT (default 1e19 m^-3), Te in eV."""
    rows = []
    with open(path) as f:
        for line in f:
            t = line.split("#", 1)[0].split()
            if len(t) >= 4:
                rows.append([float(v) for v in t[:6]] + [np.nan] * (6 - min(6, len(t))))
    if not rows:
        fail(f"{path}: no probe rows (R Z ne Te [ne_err Te_err])")
    a = np.array(rows, dtype=float)
    unit = cfg.f("PROBE_NE_UNIT")
    rel, fn, ft = cfg.f("PROBE_ERR_REL"), cfg.f("PROBE_NE_ERR_FLOOR"), cfg.f("PROBE_TE_ERR_FLOOR")
    ne, te = a[:, 2] * unit, a[:, 3]
    ne_err = np.where(np.isfinite(a[:, 4]), a[:, 4] * unit, np.maximum(rel * ne, fn * unit))
    te_err = np.where(np.isfinite(a[:, 5]), a[:, 5], np.maximum(rel * te, ft))
    return {"R": a[:, 0], "Z": a[:, 1], "ne": ne, "Te": te, "ne_err": ne_err, "Te_err": te_err}


def _components(R, Z, link):
    n = len(R)
    parent = np.arange(n)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(n):
        d = np.hypot(R - R[i], Z - Z[i])
        for j in np.where((d < link) & (np.arange(n) > i))[0]:
            a, b = find(i), find(j)
            if a != b:
                parent[a] = b
    roots = np.array([find(i) for i in range(n)])
    return [np.where(roots == r)[0] for r in np.unique(roots)]


def find_plates(bnd, probes, cfg):
    """Target plates = connected groups of wall points with |Bnorm| > 5% of max, ordered along
    the plate, s = 0 at the Psi_N minimum (strike point), s > 0 on the SOL side."""
    R, Z, psi = bnd["R"], bnd["Z"], bnd["Psi_N"]
    if "Bnorm" not in bnd:
        fail("boundary export has no Bnorm column")
    bn = np.abs(bnd["Bnorm"])
    sel = bn > 0.05 * bn.max()
    if cfg.s("PLATE_ZMAX"):
        sel &= Z < cfg.f("PLATE_ZMAX")
    idx = np.where(sel)[0]
    if idx.size < 10:
        fail("fewer than 10 wall points with |Bnorm| > 5% of the maximum; no target plates found")
    spacing = np.median(np.sort(np.hypot(R[idx][:, None] - R[idx][None, :],
                                         Z[idx][:, None] - Z[idx][None, :]) + np.eye(idx.size) * 1e9, axis=1)[:, 0])
    comps = [idx[c] for c in _components(R[idx], Z[idx], 3.5 * spacing)]
    comps = [c for c in comps if c.size >= 5]
    if probes is not None:
        keep = []
        for c in comps:
            d = np.min(np.hypot(R[c][:, None] - probes["R"][None, :], Z[c][:, None] - probes["Z"][None, :]))
            if d < cfg.f("PROBE_MAX_DIST"):
                keep.append(c)
        if keep:
            comps = keep
    comps = sorted(comps, key=lambda c: -c.size)[:2]
    if len(comps) < 1:
        fail("no target plate found")
    comps = sorted(comps, key=lambda c: R[c].mean())
    names = ["inner", "outer"] if len(comps) == 2 else ["plate"]
    plates = []
    for ic, (name, c) in enumerate(zip(names, comps)):
        other = comps[1 - ic] if len(comps) == 2 else None
        x = np.column_stack([R[c], Z[c]])
        mean = x.mean(axis=0)
        direction = np.linalg.svd(x - mean)[2][0]
        if other is not None:
            ref = np.array([R[other].mean(), Z[other].mean()])
            if np.dot(mean - ref, direction) < 0.0:
                direction = -direction
        order = np.argsort((x - mean) @ direction)
        c = c[order]
        s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(R[c]), np.diff(Z[c])))])
        i_sp = int(np.argmin(psi[c]))
        s = s - s[i_sp]
        plates.append({"name": name, "idx": c, "s": s, "R": R[c], "Z": Z[c], "psi": psi[c],
                       "R_sp": float(R[c][i_sp]), "Z_sp": float(Z[c][i_sp])})
    return plates


def project_probes(plates, probes, cfg):
    """Nearest plate point per probe -> (plate index, s, psi_N); far probes are dropped."""
    out = []
    dmax = cfg.f("PROBE_MAX_DIST")
    for i in range(len(probes["R"])):
        best = (np.inf, None, None)
        for ip, p in enumerate(plates):
            d = np.hypot(p["R"] - probes["R"][i], p["Z"] - probes["Z"][i])
            j = int(np.argmin(d))
            if d[j] < best[0]:
                best = (d[j], ip, j)
        d, ip, j = best
        if ip is None or d > dmax:
            out.append(None)
            continue
        p = plates[ip]
        # linear position between j and its nearer neighbour along the plate
        js = [k for k in (j - 1, j + 1) if 0 <= k < len(p["s"])]
        s, psi = p["s"][j], p["psi"][j]
        if js:
            k = min(js, key=lambda k: np.hypot(p["R"][k] - probes["R"][i], p["Z"][k] - probes["Z"][i]))
            seg = np.array([p["R"][k] - p["R"][j], p["Z"][k] - p["Z"][j]])
            L2 = float(seg @ seg)
            if L2 > 0:
                t = float(np.clip((np.array([probes["R"][i] - p["R"][j], probes["Z"][i] - p["Z"][j]]) @ seg) / L2, 0, 1))
                s = p["s"][j] + t * (p["s"][k] - p["s"][j])
                psi = p["psi"][j] + t * (p["psi"][k] - p["psi"][j])
        out.append({"plate": ip, "s": float(s), "psi": float(psi), "dist": float(d)})
    return out


def _weighted_slope(x, y, w):
    W = np.sum(w)
    xm, ym = np.sum(w * x) / W, np.sum(w * y) / W
    sxx = np.sum(w * (x - xm)**2)
    if sxx <= 0:
        return np.nan, np.nan
    b = np.sum(w * (x - xm) * (y - ym)) / sxx
    return ym - b * xm, b


def target_widths(plate, bnd, probes, proj, field_bnd, field_probe, cfg):
    """SOL-side e-folding length in psi_N of ln(field) along the plate: simulation vs probes.
    Returns dict or None when fewer than 3 SOL-side probes."""
    if field_bnd not in bnd:
        return None
    sel = [i for i, pr in enumerate(proj) if pr is not None and pr["plate"] == plate["id"] and pr["s"] > 0]
    if len(sel) < 3:
        return None
    xp = np.array([proj[i]["psi"] for i in sel])
    yp = probes[field_probe][sel]
    ep = probes[field_probe + "_err"][sel]
    ok = (yp > 0) & np.isfinite(yp)
    if np.count_nonzero(ok) < 3 or np.ptp(xp[ok]) < 0.005:
        return None
    xp, yp, ep = xp[ok], yp[ok], ep[ok]
    wp = 1.0 / np.maximum(ep / yp, 0.05)**2
    _, bp = _weighted_slope(xp, np.log(yp), wp)
    lo, hi = max(1.0, xp.min() - 0.005), xp.max() + 0.005
    ps, us = plate["psi"], bnd[field_bnd][plate["idx"]]
    m = (plate["s"] > 0) & (ps >= lo) & (ps <= hi) & (us > 0)
    if np.count_nonzero(m) < 3:
        return None
    order = np.argsort(ps[m])
    sim_at_probes = np.interp(xp, ps[m][order], us[m][order])
    if np.any(sim_at_probes <= 0):
        return None
    _, bs = _weighted_slope(xp, np.log(sim_at_probes), wp)          # same estimator and weights as the probes
    level = float(np.sqrt(np.mean(((sim_at_probes - yp) / yp)**2)))
    s_sorted = np.argsort(ps[m])
    dpsi_ds = np.polyfit(plate["s"][m][s_sorted], ps[m][s_sorted], 1)[0] if np.count_nonzero(m) > 3 else np.nan
    return {"lam_probe": (-1.0 / bp) if bp < 0 else np.nan, "lam_sim": (-1.0 / bs) if bs < 0 else np.nan,
            "n_probes": int(len(xp)), "psi_lo": float(lo), "psi_hi": float(hi), "level_rms": level,
            "dpsi_ds": float(dpsi_ds), "peak_sim": float(us[m].max()), "peak_probe": float(yp.max())}


# ----------------------------------------------------------------------------- channels
def channel_defs(cfg, nml):
    two = cfg.i("TEMPERATURES") == 2
    dens = cfg.s("DENSITY_EXPR")
    ch = [{"name": dens, "kind": "D", "key": "d_perp_file", "tanh": "d_perp", "expr": dens, "weighted": False,
           "bnd_expr": "ne", "probe": "ne", "target_w": cfg.f("W_TARGET")}]
    wt = cfg.f("W_TARGET") * cfg.f("TE_TARGET_WEIGHT")
    if two:
        ch.append({"name": "T_e", "kind": "ZK", "key": "zk_e_perp_file", "tanh": "zk_e_perp", "expr": "T_e",
                   "weighted": True, "bnd_expr": "T_e", "probe": "Te", "target_w": wt})
        ch.append({"name": "T_i", "kind": "ZK", "key": "zk_i_perp_file", "tanh": "zk_i_perp", "expr": "T_i",
                   "weighted": True, "bnd_expr": None, "probe": None, "target_w": 0.0})
    else:
        ch.append({"name": "T", "kind": "ZK", "key": "zk_perp_file", "tanh": "zk_perp", "expr": "T",
                   "weighted": True, "bnd_expr": "Te", "probe": "Te", "target_w": wt})
    return ch


def midplane_exprs(cfg):
    two = cfg.i("TEMPERATURES") == 2
    return ["R", "Psi_N", "rho", "ne", "T_e", "T_i"] if two else ["R", "Psi_N", "rho", "ne", "T"]


def boundary_exprs(cfg):
    two = cfg.i("TEMPERATURES") == 2
    return ["R", "Z", "Psi_N", "Bnorm", "ne", "T_e" if two else "Te", "heatF_total", "partF_total"]


# ----------------------------------------------------------------------------- state
def load_state(cfg, must_exist=True):
    p = os.path.join(cfg.s("RUN"), STATE_FILE)
    if not os.path.exists(p):
        if must_exist:
            fail(f"{p} missing: run 'transport_fit.py check' first")
        return None
    with open(p) as f:
        return json.load(f)


def save_state(cfg, st):
    p = os.path.join(cfg.s("RUN"), STATE_FILE)
    tmp = p + ".tmp"
    with open(tmp, "w") as f:
        json.dump(st, f, indent=1)
    os.replace(tmp, p)


def channel_file_values(st, ch):
    """Current table (x, y) for a channel = original table * exp(pchip(lnF))."""
    x, y0 = np.array(ch["orig_x"]), np.array(ch["orig_y"])
    lnF = pchip(np.array(st["knots"]), np.array(ch["lnF"]), x)
    return x, y0 * np.exp(lnF)


def write_channel_file(cfg, st, ch):
    x, y = channel_file_values(st, ch)
    y = np.maximum(y, 1e-30)
    write_prof(os.path.join(cfg.s("RUN"), ch["file"]), x, y)


# ----------------------------------------------------------------------------- check
def cmd_check(cfg, args):
    run = cfg.s("RUN")
    nmlp = os.path.join(run, cfg.s("NAMELIST"))
    if not os.path.isfile(nmlp):
        fail(f"namelist {nmlp} not found")
    nml = parse_namelist(nmlp)
    nf = norm_factors(nml)
    problems, warnings = [], []
    say(f"transport_fit check: run dir {os.path.abspath(run)}, namelist {cfg.s('NAMELIST')}")
    say(f"  central_density {nf['n0']:.3e} m^-3, central_mass {nml_get(nml, 'central_mass', 2.014):.4f}, "
        f"sqrt(mu0 rho0) = {nf['sqrt_mu0rho0']:.4e}")
    say(f"  D = 1 m2/s  <->  D_jorek = {1 / nf['D_to_SI']:.4e};  chi = 1 m2/s  <->  ZK_file = {1 / nf['ZK_to_SI']:.4e}")

    two = cfg.i("TEMPERATURES") == 2
    say(f"  temperatures: {'two (T_e, T_i; ZK_e_perp_file, ZK_i_perp_file)' if two else 'one (T; ZK_perp_file)'}")
    if not bool(nml_get(nml, "use_zkperp_times_density", False)):
        problems.append("use_zkperp_times_density is not .true. in the namelist; this tool assumes it is on "
                        "(the ZK file then holds a chi-like quantity multiplied by the density in the code)")
    if two and nml_get(nml, "zk_perp_file", "none") != "none" and nml_get(nml, "zk_e_perp_file", "none") == "none":
        warnings.append("TEMPERATURES=2 but only ZK_perp_file is set; with with_TiTe the code reads "
                        "ZK_e_perp_file / ZK_i_perp_file")
    if (not two) and nml_get(nml, "zk_e_perp_file", "none") != "none":
        warnings.append("TEMPERATURES=1 but ZK_e_perp_file is set; a one-temperature build reads ZK_perp_file")

    nout = int(nml_get(nml, "nout", 1))
    nstep = cfg.i("NSTEP_SEG")
    if nstep % nout:
        warnings.append(f"NSTEP_SEG={nstep} is not a multiple of nout={nout}; segments are rounded up")
    for key in ("nstep_n", "tstep_n"):
        v = nml.get(key)
        if isinstance(v, dict) and len([i for i in v if i > 1 and v[i]]) > 0:
            msg = (f"{key} has entries beyond slot 1 (a time-step ramp?): every segment runs with tstep_n(1) = "
                   f"{nml_get(nml, 'tstep_n', '?')} for nstep_n(1) steps and ignores the other slots")
            if key == "tstep_n" and not cfg.s("TSTEP"):
                problems.append(msg + "; set TSTEP in the driver to the production time step, or make tstep_n a single value")
            else:
                warnings.append(msg)
    tstep_eff = float(cfg.s("TSTEP")) if cfg.s("TSTEP") else float(nml_get(nml, "tstep_n", nml_get(nml, "tstep", 0.0)) or 0.0)
    if tstep_eff > 0:
        t_ms = cfg.i("NSTEP_SEG") * tstep_eff * nf["sqrt_mu0rho0"] * 1e3
        say(f"  segment: {cfg.i('NSTEP_SEG')} steps x tstep {tstep_eff:g} = {t_ms:.2f} ms of plasma time "
            f"(edge transport and SOL transit need of order 1-3 ms per segment)")
        if t_ms < 0.5:
            warnings.append(f"a segment covers only {t_ms:.2f} ms; the edge will not re-equilibrate and the steadiness "
                            f"gate will keep doubling it. Raise NSTEP_SEG or TSTEP")
    if "nstep_n" not in nml and "nstep" not in nml:
        warnings.append("neither nstep_n nor nstep in the namelist; 'nstep_n' will be inserted")

    h5 = sorted(glob.glob(os.path.join(run, "jorek[0-9]*.h5")))
    if not h5:
        problems.append("no jorek??????.h5 step files in the run dir (the loop restarts from the newest one)")
    rst = os.path.join(run, "jorek_restart.h5")
    if h5 and os.path.exists(rst) and os.path.getmtime(rst) > max(os.path.getmtime(p) for p in h5) + 1.0:
        warnings.append("jorek_restart.h5 is newer than the newest jorek??????.h5; the loop uses the newest step file")
    ref = os.path.join(run, f"jorek{cfg.i('REF_STEP'):06d}.h5")
    if not os.path.exists(ref):
        problems.append(f"reference step file {os.path.basename(ref)} missing (REF_STEP={cfg.i('REF_STEP')})")
    if not os.path.exists(os.path.join(run, "part_restart.h5")):
        warnings.append("part_restart.h5 not found: the first segment initialises the kinetic particles "
                        "(restart_particles=.f.); later segments continue them")
    pf = cfg.path("PROBE_FILE")
    if pf:
        if not os.path.isfile(pf):
            problems.append(f"PROBE_FILE {pf} not found")
        else:
            pr = read_probes(pf, cfg)
            say(f"  probe file: {len(pr['R'])} probes, ne {pr['ne'].min():.2e}..{pr['ne'].max():.2e} m^-3, "
                f"Te {pr['Te'].min():.1f}..{pr['Te'].max():.1f} eV")
    else:
        warnings.append("no PROBE_FILE: upstream-only fit, target profiles are reported only")

    # --- profile files and state
    st = load_state(cfg, must_exist=False)
    chans = channel_defs(cfg, nml)
    if st is None:
        st = {"iter": 0, "nstep": nstep, "clamp": cfg.f("CLAMP"), "crashes": 0, "done": 0, "converged_count": 0,
              "knots": None, "psi_max": None, "channels": {}, "history": [], "ref_step": None}
        assignments = {}
        for ch in chans:
            fname = nml_get(nml, ch["key"], "none")
            if fname in ("none", None, ""):
                c = nml_array(nml, ch["tanh"], 6)
                if ch["tanh"] == "zk_e_perp" and not np.any(c):
                    c = nml_array(nml, "zk_perp", 6)
                if ch["tanh"] == "zk_i_perp" and not np.any(c):
                    c = nml_array(nml, "zk_perp", 6)
                if c[0] == 0.0:
                    problems.append(f"{ch['key']} is 'none' and {ch['tanh']}(1) is 0: nothing to bootstrap from")
                    continue
                x = np.linspace(0.0, 1.5, 200)
                y = tanh_profile(c, x)
                fname = f"{ch['tanh']}_tfit.dat"
                write_prof(os.path.join(run, fname), x, y)
                assignments[ch["key"]] = f"'{fname}'"
                say(f"  {ch['key']} was 'none': bootstrapped {fname} from the tanh coefficients "
                    f"{', '.join(f'{v:g}' for v in c)} (psi_N 0..1.5)")
            else:
                src = os.path.join(run, fname)
                if not os.path.isfile(src):
                    problems.append(f"{ch['key']} = {fname} not found")
                    continue
            x, y = read_prof(os.path.join(run, fname))
            if np.any(y <= 0):
                problems.append(f"{fname}: non-positive values (JOREK stops on them)")
            user_copy = os.path.join(run, fname + ".user")
            if not os.path.exists(user_copy):
                shutil.copy2(os.path.join(run, fname), user_copy)
            st["channels"][ch["name"]] = {"file": fname, "key": ch["key"], "kind": ch["kind"], "expr": ch["expr"],
                                          "weighted": ch["weighted"], "orig_x": x.tolist(), "orig_y": y.tolist(),
                                          "lnF": None, "prev_lnF": None, "s_est": None, "last_dlnF": None,
                                          "last_lng": None, "flags": None, "active": ch["kind"] == "D"}
        if assignments:
            set_namelist(nmlp, assignments)
        if not problems:
            save_state(cfg, st)
            say(f"  state initialised in {STATE_FILE} (knots are set by 'reference' from the exported psi_N range)")
    else:
        say(f"  existing state: iteration {st['iter']}, nstep {st['nstep']}, clamp {st['clamp']:.2f}, "
            f"crashes {st['crashes']}, done {st['done']}")

    for ch in chans:
        c = st["channels"].get(ch["name"])
        if not c:
            continue
        x, y = np.array(c["orig_x"]), np.array(c["orig_y"])
        f = nf["D_to_SI"] if c["kind"] == "D" else nf["ZK_to_SI"]
        say(f"  {ch['name']:4s} <- {c['file']}: {len(x)} points, psi_N {x[0]:.3f}..{x[-1]:.3f}, "
            f"{'D' if c['kind'] == 'D' else 'chi'} = {y.min() * f:.3g}..{y.max() * f:.3g} m2/s")
        if len(x) < 50 or len(x) > 500:
            warnings.append(f"{c['file']}: {len(x)} points (JOREK recommends 50..500); the written file keeps this grid")

    for w in warnings:
        say("  WARNING: " + w)
    for p in problems:
        say("  PROBLEM: " + p)
    if problems:
        sys.exit(2)
    say("check: OK")


# ----------------------------------------------------------------------------- reference
def export_paths(cfg, step):
    pp = os.path.join(cfg.s("RUN"), "postproc")
    return (os.path.join(pp, f"exprs_outer-midplane_s{step:06d}.dat"),
            os.path.join(pp, f"boundary_quantities_s{step:06d}.dat"))


def cmd_reference(cfg, args):
    run = cfg.s("RUN")
    st = load_state(cfg)
    nml = parse_namelist(os.path.join(run, cfg.s("NAMELIST")))
    midp, bndp = export_paths(cfg, args.step)
    for p in (midp, bndp):
        if not os.path.isfile(p):
            fail(f"reference export {p} missing (jorek2_postproc failed?)")
    mid, bnd = read_midplane(midp), read_boundary(bndp)
    for e in midplane_exprs(cfg)[1:]:
        if e not in mid:
            fail(f"midplane export lacks '{e}' (available: {', '.join(mid)})")
    shutil.copy2(midp, os.path.join(run, REF_MID))
    shutil.copy2(bndp, os.path.join(run, REF_BND))
    psi_edge = float(mid["Psi_N"].max())
    psi_wall = float(bnd["Psi_N"].max())
    psi_max = max(psi_edge, psi_wall)
    say(f"reference step {args.step}: midplane {len(mid['Psi_N'])} points, psi_N up to {psi_edge:.3f} at the grid edge; "
        f"wall psi_N up to {psi_wall:.3f} (PFR mirrored)")
    if psi_edge < 1.02:
        fail("the outer-midplane cut does not reach psi_N = 1.02: no SOL on the midplane, nothing to fit there")

    knots = np.array([float(v) for v in cfg.s("KNOTS").split(",")]) if cfg.s("KNOTS") else default_knots(psi_edge)
    knots = np.unique(knots[knots <= psi_edge + 1e-9])
    st["knots"], st["psi_max"], st["psi_edge"], st["ref_step"] = knots.tolist(), psi_max, psi_edge, int(args.step)
    nf = norm_factors(nml)
    for name, c in st["channels"].items():
        x, y = np.array(c["orig_x"]), np.array(c["orig_y"])
        if x[-1] < psi_max + 0.02:
            say(f"  {c['file']}: table ends at psi_N {x[-1]:.3f} < grid maximum {psi_max:.3f}: extended flat "
                f"(JOREK extrapolates linearly outside the table, which can go negative)")
            x, y = np.append(x, psi_max + 0.05), np.append(y, y[-1])
            c["orig_x"], c["orig_y"] = x.tolist(), y.tolist()
        if c["lnF"] is None:
            c["lnF"] = [0.0] * len(knots)
            c["s_est"] = [None] * len(knots)
            c["flags"] = [""] * len(knots)
        write_channel_file(cfg, st, c)
        f = nf["D_to_SI"] if c["kind"] == "D" else nf["ZK_to_SI"]
        vk = pchip(x, y, knots) * np.exp(np.array(c["lnF"]))
        say(f"  {name:4s}: {len(knots)} knots, "
            f"{'D' if c['kind'] == 'D' else 'chi'} at psi_N 0.9/1.0/1.1 = "
            f"{np.interp(0.9, knots, vk) * f:.3g}/{np.interp(1.0, knots, vk) * f:.3g}/{np.interp(1.1, knots, vk) * f:.3g} m2/s")
        v, s = profile_at_knots(mid, c["expr"], knots)
        rel = np.abs(s) / np.maximum(np.abs(v), 1e-300)
        nflat = int(np.count_nonzero((rel < cfg.f("FLAT_REL")) & (knots >= cfg.f("PSI_MIN"))))
        say(f"        reference {c['expr']}: {nflat} flat knots in the fit range (held), "
            f"value at psi_N 1 = {np.interp(1.0, knots, v):.3g}")

    pf = cfg.path("PROBE_FILE")
    probes = read_probes(pf, cfg) if pf else None
    plates = find_plates(bnd, probes, cfg)
    for i, p in enumerate(plates):
        p["id"] = i
        say(f"  {p['name']} plate: {p['idx'].size} wall points, s = {p['s'].min() * 100:.1f}..{p['s'].max() * 100:.1f} cm, "
            f"strike point R = {p['R_sp']:.3f} Z = {p['Z_sp']:.3f} m, psi_N {p['psi'].min():.3f}..{p['psi'].max():.3f}")
    st["plates"] = [{"name": p["name"], "idx": p["idx"].tolist()} for p in plates]
    if probes is not None:
        proj = project_probes(plates, probes, cfg)
        n_ok = sum(1 for q in proj if q is not None)
        say(f"  probes: {n_ok}/{len(proj)} within {cfg.f('PROBE_MAX_DIST') * 100:.1f} cm of a plate")
        for i, q in enumerate(proj):
            if q is None:
                say(f"    probe {i + 1} at R={probes['R'][i]:.3f} Z={probes['Z'][i]:.3f}: no plate nearby, dropped")
        for p in plates:
            n_sol = sum(1 for q in proj if q and q["plate"] == p["id"] and q["s"] > 0)
            n_pfr = sum(1 for q in proj if q and q["plate"] == p["id"] and q["s"] <= 0)
            say(f"    {p['name']}: {n_sol} probes on the SOL side (fitted), {n_pfr} on the PFR side (reported only)")
            if n_sol < 3:
                say(f"    WARNING: fewer than 3 SOL-side probes on the {p['name']} plate: no target width from it")
    save_state(cfg, st)
    say("reference: OK")


# ----------------------------------------------------------------------------- status
def cmd_status(cfg, args):
    st = load_state(cfg)
    nml = parse_namelist(os.path.join(cfg.s("RUN"), cfg.s("NAMELIST")))
    nout = max(1, int(nml_get(nml, "nout", 1)))
    nstep = int(math.ceil(st["nstep"] / nout) * nout)
    done = 1 if (st["done"] or st["iter"] >= cfg.i("MAX_ITER")) else 0
    has_parts = os.path.exists(os.path.join(cfg.s("RUN"), "part_restart.h5"))
    print(f"ITER={st['iter']} NSTEP={nstep} DONE={done} CLAMP={st['clamp']:.3f} CRASHES={st['crashes']} "
          f"RESTART_PARTICLES={fortran_logical(has_parts)} NOUT={nout}")


def cmd_set_namelist(cfg, args):
    assignments = {}
    for kv in args.set:
        k, v = kv.split("=", 1)
        assignments[k] = v
    set_namelist(os.path.join(cfg.s("RUN"), cfg.s("NAMELIST")), assignments)


# ----------------------------------------------------------------------------- update
def smooth3(v, mask):
    out = v.copy()
    idx = np.where(mask)[0]
    for n, i in enumerate(idx):
        a = v[idx[n - 1]] if n > 0 else v[i]
        b = v[idx[n + 1]] if n + 1 < len(idx) else v[i]
        out[i] = 0.25 * a + 0.5 * v[i] + 0.25 * b
    return out


def rel_rms(u, ref, mask):
    scale = np.maximum(np.abs(ref[mask]), 0.05 * np.nanmax(np.abs(ref[mask])))
    return float(np.sqrt(np.nanmean(((u[mask] - ref[mask]) / scale)**2)))


def cmd_update(cfg, args):
    run = cfg.s("RUN")
    st = load_state(cfg)
    if st.get("knots") is None:
        fail("no knots in the state: run 'reference' first")
    nml = parse_namelist(os.path.join(run, cfg.s("NAMELIST")))
    nf = norm_factors(nml)
    knots = np.array(st["knots"])
    it = st["iter"] + 1
    relax, clamp, tol = cfg.f("RELAX"), float(st["clamp"]), cfg.f("TOL")
    psi_min = cfg.f("PSI_MIN")
    fit = (knots >= psi_min) & (knots <= st["psi_edge"] - cfg.f("EDGE_HOLD") + 1e-9)

    ref = read_midplane(os.path.join(run, REF_MID))
    midp1, _ = export_paths(cfg, args.mid)
    midp2, bndp2 = export_paths(cfg, args.last)
    for p in (midp1, midp2, bndp2):
        if not os.path.isfile(p):
            fail(f"export {p} missing")
    first, last, bnd = read_midplane(midp1), read_midplane(midp2), read_boundary(bndp2)   # 'first' = mid-segment export
    dens = cfg.s("DENSITY_EXPR")
    n_ref, _ = profile_at_knots(ref, dens, knots)
    n_cur, _ = profile_at_knots(last, dens, knots)

    # --- target data
    probes = read_probes(cfg.path("PROBE_FILE"), cfg) if cfg.path("PROBE_FILE") else None
    plates = find_plates(bnd, probes, cfg)
    for i, p in enumerate(plates):
        p["id"] = i
    proj = project_probes(plates, probes, cfg) if probes is not None else None
    tgt_info = {}
    if probes is not None:
        for p in plates:
            for fb, fp in (("ne", "ne"), (boundary_exprs(cfg)[5], "Te")):
                r = target_widths(p, bnd, probes, proj, fb, fp, cfg)
                if r is not None:
                    tgt_info[(p["name"], fp)] = r

    row = {"iter": it, "first": args.first, "last": args.last, "nstep": st["nstep"], "clamp": clamp}
    say(f"update: iteration {it}, segment steps {args.first}..{args.last} (steadiness from step {args.mid}), "
        f"clamp {clamp:.2f}, relax {relax:.2f}")

    # --- activate heat channels
    dens_ch = st["channels"][dens]
    mism_dens_prev = st["history"][-1].get(f"mism_{dens}", 1.0) if st["history"] else 1.0
    heat_on = it >= cfg.i("FIT_HEAT_FROM_ITER") or mism_dens_prev < 2.0 * tol
    for name, c in st["channels"].items():
        if c["kind"] == "ZK":
            c["active"] = bool(heat_on)

    # --- steadiness gate (all channels, fit range, on the midplane points)
    not_steady = []
    for name, c in st["channels"].items():
        e = c["expr"]
        psi_pts = last["Psi_N"]
        m = (psi_pts >= psi_min) & (psi_pts <= st["psi_edge"])
        u_last = last[e]
        u_first = np.interp(psi_pts, first["Psi_N"], first[e])
        u_ref = np.interp(psi_pts, ref["Psi_N"], ref[e])
        mism = rel_rms(u_last, u_ref, m)
        drift = rel_rms(u_last, u_first, m)
        row[f"mism_{name}"], row[f"drift_{name}"] = mism, drift
        bands = [(psi_min, 0.9), (0.9, 1.0), (1.0, 1.1), (1.1, st["psi_edge"])]
        btxt = "  ".join(f"{lo:.2f}-{hi:.2f}: {rel_rms(u_last, u_ref, m & (psi_pts >= lo) & (psi_pts < hi)) * 100:4.1f}%"
                         for lo, hi in bands if np.count_nonzero(m & (psi_pts >= lo) & (psi_pts < hi)) > 3)
        say(f"  {name:4s}: mismatch to reference {mism * 100:5.1f}%  drift over the second half of the segment "
            f"{drift * 100:5.1f}%{'' if c['active'] else '  (not fitted yet)'}")
        say(f"        by band  {btxt}")
        if c["active"] and mism > 0.5 * tol and drift > cfg.f("STEADY_RATIO") * mism and drift > 0.5 * tol:
            not_steady.append(name)
    if not_steady:
        if st["nstep"] < cfg.i("NSTEP_MAX"):
            st["nstep"] = min(2 * st["nstep"], cfg.i("NSTEP_MAX"))
            row["code"] = EXIT_RERUN
            st["history"].append(row)
            write_log(cfg, st)
            save_state(cfg, st)
            say(f"  not steady ({', '.join(not_steady)}): segment length doubled to {st['nstep']} steps, "
                f"no update, rerun from the same state")
            plot_iteration(cfg, st, it, ref, first, last, bnd, plates, probes, proj, tgt_info, nf, suffix="_unsteady")
            sys.exit(EXIT_RERUN)
        say(f"  not steady ({', '.join(not_steady)}) at NSTEP_MAX: updating with half the relaxation")
        relax *= 0.5

    # --- gradient ratios and steps
    maxstep, nflag_total = 0.0, 0
    for name, c in st["channels"].items():
        v_ref, s_ref = profile_at_knots(ref, c["expr"], knots)
        v_cur, s_cur = profile_at_knots(last, c["expr"], knots)
        if c["weighted"]:
            g = (n_cur * s_cur) / (n_ref * s_ref)
        else:
            g = s_cur / s_ref
        rel_ref = np.abs(s_ref) / np.maximum(np.abs(v_ref), 1e-300)
        flat = rel_ref < cfg.f("FLAT_REL")
        v_sep = float(np.interp(1.0, ref["Psi_N"], ref[c["expr"]]))
        flags = [""] * len(knots)
        lng = np.full(len(knots), np.nan)
        for k in range(len(knots)):
            if not c["active"] or knots[k] < psi_min:
                flags[k] = "held"
                continue
            if not fit[k]:
                flags[k] = "edge-held"
                continue
            if not np.isfinite(s_ref[k]) or not np.isfinite(s_cur[k]):
                flags[k] = "nodata"
                continue
            if flat[k]:
                flags[k] = "flat-ref"
                continue
            if knots[k] > 1.0 and np.isfinite(v_sep) and v_ref[k] < cfg.f("SMALL_FRAC") * v_sep:
                flags[k] = "small-ref"       # below a few % of the separatrix value: nothing to fit against
                continue
            if g[k] <= 0.05:          # simulated profile flat or hollow where the reference is not
                if v_cur[k] < 0.2 * v_ref[k]:
                    flags[k] = "low-sim"          # the whole profile is low here: the gradient says nothing about D
                    continue
                lng[k] = -math.log(2.0)           # bounded push towards less transport
                flags[k] = "flat-sim"
                continue
            lng[k] = math.log(g[k])
        # level correction on the held edge band next to the Dirichlet grid edge: the profile level inside is
        # set by the gradient across this band, so scale its coefficient by the level ratio at the last fitted knot
        if c["active"]:
            k_fit = np.where(fit & np.isfinite(v_cur) & np.isfinite(v_ref) & np.isfinite(lng))[0]
            band = np.array([fl == "edge-held" for fl in flags]) & np.isfinite(s_ref)
            adjacent = k_fit.size and not any(flags[k] in ("small-ref", "low-sim", "nodata") for k in range(k_fit[-1] + 1, len(knots)))
            if k_fit.size and band.any() and adjacent:
                kl = k_fit[-1]
                u_edge = float(ref[c["expr"]][-1])                     # Dirichlet value at the grid edge
                if v_ref[kl] > 0 and v_cur[kl] > 0:
                    # d ln u_k / d ln D_band = -(u_k - u_edge)/u_k for an outward flux (falling profile);
                    # the sign flips when the band is fed from the edge (rising reference profile)
                    sign = -1.0 if np.nanmean(s_ref[band]) > 0 else 1.0
                    gain = min(5.0, v_ref[kl] / max(abs(v_ref[kl] - u_edge), 0.2 * v_ref[kl]))
                    lng_level = sign * gain * math.log(v_cur[kl] / v_ref[kl])
                    for k in range(len(knots)):
                        if band[k]:
                            lng[k] = lng_level
                            flags[k] = "edge-level"
                    row[f"lng_edge_{name}"] = lng_level
        # blend the target widths into the SOL knots
        tw = 0.0
        if probes is not None and c.get("kind") and c["active"]:
            chdef = next(d for d in channel_defs(cfg, nml) if d["name"] == name)
            if chdef["probe"]:
                tw = chdef["target_w"]
                num, den, hi = 0.0, 0.0, 1.0
                for (pname, fp), r in tgt_info.items():
                    if fp == chdef["probe"] and np.isfinite(r["lam_probe"]) and np.isfinite(r["lam_sim"]):
                        num += r["n_probes"] * math.log(r["lam_probe"] / r["lam_sim"])
                        den += r["n_probes"]
                        hi = max(hi, r["psi_hi"])
                if den > 0 and tw > 0:
                    lng_t = num / den
                    for k in range(len(knots)):
                        if 1.0 <= knots[k] <= hi + 0.02 and np.isfinite(lng[k]):
                            lng[k] = (1.0 - tw) * lng[k] + tw * lng_t
                        elif 1.0 <= knots[k] <= hi + 0.02 and flags[k] in ("flat-ref", "nodata"):
                            lng[k] = tw * lng_t
                            flags[k] = ""
                    row[f"lng_target_{name}"] = lng_t
        # sensitivity from the previous step
        s_est = np.array([np.nan if v is None else v for v in (c["s_est"] or [None] * len(knots))], dtype=float)
        if c["last_dlnF"] is not None and c["last_lng"] is not None:
            d_prev = np.array(c["last_dlnF"], dtype=float)
            lng_prev = np.array([np.nan if v is None else v for v in c["last_lng"]], dtype=float)
            for k in range(len(knots)):
                if abs(d_prev[k]) > cfg.f("SENS_MIN_STEP") and np.isfinite(lng_prev[k]) and np.isfinite(lng[k]):
                    s_obs = (lng[k] - lng_prev[k]) / d_prev[k]
                    s_est[k] = s_obs if not np.isfinite(s_est[k]) else 0.5 * s_est[k] + 0.5 * s_obs
        # fixed-point step with the pure-diffusion slope (-1); the measured sensitivity only ever
        # reduces a step (a clearly wrong sign halves it and flags the knot), never enlarges it
        clamp_ch = min(clamp, 1.5) if c.get("first_active", True) else clamp
        step = np.zeros(len(knots))
        for k in range(len(knots)):
            if not np.isfinite(lng[k]):
                continue
            r = relax
            if np.isfinite(s_est[k]) and s_est[k] > 0.2:
                r = 0.5 * relax
                flags[k] = (flags[k] + " wrong-sign").strip()
            step[k] = float(np.clip(r * lng[k], -math.log(clamp_ch), math.log(clamp_ch)))
        active = np.isfinite(lng)
        if cfg.b("SMOOTH"):
            inner = active & np.array([fl != "edge-level" for fl in flags])
            step = np.where(inner, smooth3(np.where(inner, step, 0.0), inner), step)
        lnF = np.array(c["lnF"]) + step
        # bounds in SI
        x0, y0 = np.array(c["orig_x"]), np.array(c["orig_y"])
        yk = pchip(x0, y0, knots)
        f = nf["D_to_SI"] if c["kind"] == "D" else nf["ZK_to_SI"]
        lo = cfg.f("D_MIN") if c["kind"] == "D" else cfg.f("CHI_MIN")
        hi = cfg.f("D_MAX") if c["kind"] == "D" else cfg.f("CHI_MAX")
        if c["kind"] == "ZK" and it <= 3:
            edge = knots >= st["psi_edge"] - 0.1
            hi_k = np.where(edge, min(hi, cfg.f("CHI_EDGE_MAX")), hi)
        else:
            hi_k = np.full(len(knots), hi)
        lnF = np.minimum(lnF, np.log(hi_k / f / yk))
        lnF = np.maximum(lnF, np.log(lo / f / yk))
        step = lnF - np.array(c["lnF"])
        # commit
        backup = os.path.join(run, f"{c['file']}.iter{it - 1:02d}")
        shutil.copy2(os.path.join(run, c["file"]), backup)
        c["prev_lnF"], c["lnF"] = c["lnF"], lnF.tolist()
        c["last_dlnF"] = step.tolist()
        c["last_lng"] = [None if not np.isfinite(v) else float(v) for v in lng]
        c["s_est"] = [None if not np.isfinite(v) else float(v) for v in s_est]
        c["flags"] = flags
        c["first_active"] = not c["active"]
        c.setdefault("lnF_history", []).append([round(float(v), 4) for v in lnF])
        write_channel_file(cfg, st, c)
        nflag = sum(1 for fl in flags if "wrong-sign" in fl or fl in ("flat-sim", "low-sim"))
        nflag_total += nflag
        maxstep = max(maxstep, float(np.max(np.abs(step))) if len(step) else 0.0)
        row[f"maxfactor_{name}"] = float(np.exp(np.max(np.abs(step)))) if len(step) else 1.0
        row[f"nflag_{name}"] = nflag
        if c["active"]:
            vk = yk * np.exp(lnF) * f
            say(f"  {name:4s}: max factor this step {np.exp(np.max(np.abs(step))):.2f}, flagged knots {nflag}, "
                f"{'D' if c['kind'] == 'D' else 'chi'} at psi_N 0.9/1.0/1.1 = "
                f"{np.interp(0.9, knots, vk):.3g}/{np.interp(1.0, knots, vk):.3g}/{np.interp(1.1, knots, vk):.3g} m2/s"
                + (f", target weight {tw:.2f}" if tw > 0 else ""))
            for k in range(len(knots)):
                if flags[k] and flags[k] not in ("held", "edge-held", "edge-level", "small-ref", "flat-ref"):
                    say(f"        knot psi_N {knots[k]:.3f}: {flags[k]}"
                        + (f", sensitivity {s_est[k]:+.2f}" if np.isfinite(s_est[k]) else ""))

    # --- target report
    for (pname, fp), r in tgt_info.items():
        say(f"  target {pname} {fp}: e-folding in psi_N sim {r['lam_sim']:.3f} vs probes {r['lam_probe']:.3f}"
            f" ({r['n_probes']} probes), level rms {r['level_rms'] * 100:.0f}%, peak sim/probe "
            f"{r['peak_sim']:.3g}/{r['peak_probe']:.3g}")
        row[f"lam_sim_{pname}_{fp}"], row[f"lam_probe_{pname}_{fp}"] = r["lam_sim"], r["lam_probe"]
        row[f"level_{pname}_{fp}"] = r["level_rms"]
    if "partF_total" in bnd and "heatF_total" in bnd:
        row["wall_partF_max"] = float(np.nanmax(np.abs(bnd["partF_total"])))
        row["wall_heatF_max"] = float(np.nanmax(np.abs(bnd["heatF_total"])))

    # --- convergence
    active_names = [n for n, c in st["channels"].items() if c["active"]]
    all_ok = all(row[f"mism_{n}"] < tol for n in active_names) and heat_on
    small = maxstep < math.log(1.03) and heat_on
    st["converged_count"] = st["converged_count"] + 1 if all_ok else 0
    st["iter"] = it
    code = EXIT_CONTINUE
    if st["converged_count"] >= 2:
        code = EXIT_CONVERGED
        st["done"] = 1
        say(f"  CONVERGED: all channels within {tol * 100:.0f}% for two iterations")
    elif small:
        code = EXIT_CONVERGED
        st["done"] = 1
        say("  coefficient change below 3% everywhere: stopping (profiles within "
            + ", ".join(f"{row[f'mism_{n}'] * 100:.1f}%" for n in active_names) + ")")
    elif it >= cfg.i("MAX_ITER"):
        code = EXIT_STOP
        st["done"] = 1
        say(f"  MAX_ITER={cfg.i('MAX_ITER')} reached")
    row["code"] = code
    st["history"].append(row)
    write_log(cfg, st)
    save_state(cfg, st)
    plot_iteration(cfg, st, it, ref, first, last, bnd, plates, probes, proj, tgt_info, nf)
    sys.exit(code)


def write_log(cfg, st):
    keys = []
    for r in st["history"]:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(cfg.s("RUN"), LOG_FILE), "w") as f:
        f.write(",".join(keys) + "\n")
        for r in st["history"]:
            f.write(",".join("" if k not in r else (f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k])) for k in keys) + "\n")


# ----------------------------------------------------------------------------- crash
def cmd_crash(cfg, args):
    run = cfg.s("RUN")
    st = load_state(cfg)
    st["crashes"] += 1
    restored = []
    for name, c in st["channels"].items():
        if c.get("prev_lnF") is not None:
            c["lnF"] = c["prev_lnF"]
            c["last_dlnF"], c["last_lng"] = None, None
            write_channel_file(cfg, st, c)
            restored.append(c["file"])
    st["clamp"] = max(1.1, 0.5 * (st["clamp"] - 1.0) + 1.0)
    if st["iter"] > 0:
        st["iter"] -= 1
        if st["history"] and st["history"][-1].get("code") in (EXIT_CONTINUE, EXIT_CONVERGED):
            st["history"][-1]["code"] = "crashed"
    say(f"crash: segment failed (crash {st['crashes']}); restored {', '.join(restored) or 'nothing'}; "
        f"clamp now {st['clamp']:.2f}; iteration counter back to {st['iter']}")
    if st["crashes"] >= 3:
        st["done"] = 1
        say("crash: three failures, stopping. Look at the flagged knots in the last figure and the log.")
    save_state(cfg, st)
    write_log(cfg, st)


# ----------------------------------------------------------------------------- report / plot
def cmd_report(cfg, args):
    st = load_state(cfg)
    if not st["history"]:
        say("no iterations yet")
        return
    names = list(st["channels"])
    say(f"{'it':>3s} {'steps':>13s} " + " ".join(f"{'mism_' + n:>10s}" for n in names)
        + " " + " ".join(f"{'fac_' + n:>8s}" for n in names) + "  code")
    for r in st["history"]:
        say(f"{r['iter']:3d} {r['first']:6d}..{r['last']:<6d} "
            + " ".join(f"{r.get('mism_' + n, float('nan')) * 100:9.1f}%" for n in names) + " "
            + " ".join(f"{r.get('maxfactor_' + n, float('nan')):8.2f}" for n in names) + f"  {r.get('code', '')}")
    last = st["history"][-1]
    for k, v in last.items():
        if k.startswith("lam_sim_"):
            tag = k[len("lam_sim_"):]
            say(f"  target {tag}: width sim {v:.3f} vs probes {last.get('lam_probe_' + tag, float('nan')):.3f} (psi_N), "
                f"level rms {last.get('level_' + tag, float('nan')) * 100:.0f}%")
    say(f"state: iteration {st['iter']}, done {st['done']}, crashes {st['crashes']}, clamp {st['clamp']:.2f}")
    say(f"profile files: " + ", ".join(c["file"] for c in st["channels"].values()) + "  (originals kept as .user)")


def plot_iteration(cfg, st, it, ref, first, last, bnd, plates, probes, proj, tgt_info, nf, suffix=""):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    knots = np.array(st["knots"])
    names = list(st["channels"])
    ncol = max(3, len(names))
    fig, ax = plt.subplots(3, ncol, figsize=(4.2 * ncol, 10))
    for j, name in enumerate(names):
        c = st["channels"][name]
        e = c["expr"]
        a = ax[0, j]
        a.plot(ref["Psi_N"], ref[e], "k-", lw=2, label="reference (step 0)")
        a.plot(first["Psi_N"], first[e], "c--", lw=1, label="mid segment")
        a.plot(last["Psi_N"], last[e], "r-", lw=1.2, label="segment end")
        a.axvline(1.0, color="gray", lw=0.5)
        a.axvspan(0, cfg.f("PSI_MIN"), color="gray", alpha=0.1)
        a.set_xlim(max(0.0, cfg.f("PSI_MIN") - 0.1), st["psi_edge"])
        a.set_title(f"{e}  it {it}")
        a.set_xlabel("psi_N")
        if j == 0:
            a.legend(fontsize=8)
        b = ax[1, j]
        x0, y0 = np.array(c["orig_x"]), np.array(c["orig_y"])
        f = nf["D_to_SI"] if c["kind"] == "D" else nf["ZK_to_SI"]
        xx, yy = channel_file_values(st, c)
        b.plot(x0, y0 * f, "k:", label="original")
        b.plot(xx, yy * f, "b-", label="current")
        vk = pchip(x0, y0, knots) * np.exp(np.array(c["lnF"])) * f
        b.plot(knots, vk, "bo", ms=3)
        for k, fl in enumerate(c["flags"] or []):
            if fl and fl not in ("held", "edge-held", "edge-level"):
                b.plot(knots[k], vk[k], "rx" if ("wrong" in fl or "flat-sim" in fl) else "y^", ms=7)
        b.set_yscale("log")
        b.set_xlim(max(0.0, cfg.f("PSI_MIN") - 0.1), st["psi_edge"])
        b.set_ylabel("D [m2/s]" if c["kind"] == "D" else "chi [m2/s]")
        b.set_xlabel("psi_N")
        b.axvline(1.0, color="gray", lw=0.5)
        if j == 0:
            b.legend(fontsize=8)
    for j in range(len(names), ncol):
        ax[0, j].axis("off")
        ax[1, j].axis("off")
    # target panels
    te_b = boundary_exprs(cfg)[5]
    for j, (fb, fp, lab) in enumerate([("ne", "ne", "n_e [m^-3]"), (te_b, "Te", "T_e [eV]")]):
        a = ax[2, j]
        for p, col in zip(plates, ("tab:blue", "tab:red")):
            if fb in bnd:
                a.plot(p["s"] * 100, bnd[fb][p["idx"]], "-", color=col, label=f"{p['name']} sim")
            if probes is not None:
                ii = [i for i, q in enumerate(proj) if q is not None and q["plate"] == p["id"]]
                if ii:
                    a.errorbar([proj[i]["s"] * 100 for i in ii], probes[fp][ii], yerr=probes[fp + "_err"][ii],
                               fmt="o", color=col, ms=4, label=f"{p['name']} probes")
        a.axvline(0, color="gray", lw=0.5)
        a.set_xlabel("distance from strike point along the plate [cm]  (>0 SOL)")
        a.set_ylabel(lab)
        if fp == "ne":
            a.set_yscale("log")
        a.legend(fontsize=8)
    a = ax[2, 2]
    a.axis("off")
    txt = [f"iteration {it}" + (" (not steady)" if suffix else "")]
    for r in st["history"][-1:]:
        for n in names:
            if f"mism_{n}" in r:
                txt.append(f"{n}: mismatch {r[f'mism_{n}'] * 100:.1f}%, drift {r[f'drift_{n}'] * 100:.1f}%")
    for (pname, fp), r in tgt_info.items():
        txt.append(f"{pname} {fp}: lambda_psi sim {r['lam_sim']:.3f} probes {r['lam_probe']:.3f}")
    a.text(0.0, 1.0, "\n".join(txt), va="top", fontsize=9, family="monospace")
    for j in range(3, ncol):
        ax[2, j].axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(cfg.s("RUN"), f"tfit_iter{it:02d}{suffix}.png"), dpi=110)
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["check", "reference", "status", "update", "crash", "report", "set-namelist",
                                        "exprs"])
    ap.add_argument("--cfg", required=True)
    ap.add_argument("--step", type=int, default=None)
    ap.add_argument("--first", type=int, default=None)
    ap.add_argument("--mid", type=int, default=None)
    ap.add_argument("--last", type=int, default=None)
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="for set-namelist")
    args = ap.parse_args()
    cfg = Cfg(args.cfg)
    if args.command == "check":
        cmd_check(cfg, args)
    elif args.command == "reference":
        if args.step is None:
            args.step = cfg.i("REF_STEP")
        cmd_reference(cfg, args)
    elif args.command == "status":
        cmd_status(cfg, args)
    elif args.command == "update":
        if args.first is None or args.last is None:
            fail("update needs --first, --mid and --last")
        if args.mid is None:
            args.mid = args.first
        cmd_update(cfg, args)
    elif args.command == "crash":
        cmd_crash(cfg, args)
    elif args.command == "report":
        cmd_report(cfg, args)
    elif args.command == "set-namelist":
        cmd_set_namelist(cfg, args)
    elif args.command == "exprs":
        print("MIDPLANE=" + " ".join(midplane_exprs(cfg)))
        print("BOUNDARY=" + " ".join(boundary_exprs(cfg)))


if __name__ == "__main__":
    main()
