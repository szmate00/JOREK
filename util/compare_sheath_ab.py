#!/usr/bin/env python3
"""
compare_sheath_ab.py -- A/B analysis of two JOREK runs from their SHEATH log output.

Written for the thermoelectric_ohm A/B test: run the same case twice, once with
thermoelectric_ohm = .false. (reference) and once with .true. (test), and ask whether
the 0.71*grad_par(Te) thermal force in Ohm's law changes what the sheath does.

    ./compare_sheath_ab.py nothermo.log thermo.log

The thermoelectric term drives current along field lines from the hot target to the cold
one, so its signature is in the ANTISYMMETRIC part of the target currents -- the loop
current I_loop = (I_inner - I_outer)/2 -- and not in the net wall current, which is set by
quasi-neutrality. Both are reported, plus the sheath potential asymmetry that any such
current has to produce.

ALIGNMENT. The two runs are compared over the same interval of PHYSICAL TIME (t_now), not
the same step numbers. JOREK's tstep is adaptive, so the step index is not a clock: two
runs that have diverged will sit at different t_now at the same step, and a step-aligned
comparison then measures the difference between two different times. The runs need not
have the same number of steps, the same tstep, or the same output cadence -- averages are
time-weighted (trapezoidal in t_now), so a run that took smaller steps does not get more
weight for it. Step alignment is available with --align step if a log has no t_now.

GATES. Every comparison is gated on the two runs being converged, steady, and imposing the
same constraint. Those gates are not decoration: reading a transient as a result is the
most common way to get a wrong answer out of these logs, and the script refuses to quote a
delta it cannot support.

Usage
-----
    compare_sheath_ab.py REF.log TEST.log [options]

    --labels A B        names for the two runs        (default: file basenames)
    --window N          average over the last N steps of the SHORTER run, within the
                        common time overlap           (default: 200)
    --from-time T       explicit time window instead of --window
    --to-time T
    --align {auto,time,step}                          (default: auto -> time if available)
    --weak-tol X        convergence gate on the weak residual   (default: 1e-2)
    --plot FILE.png     time-series comparison figure
    --csv PREFIX        dump the parsed series to PREFIX_{a,b}.csv
    --all-blocks        keep every SHEATH block; default keeps the last one per step
                        (construct_matrix can be called more than once per step)
"""

import argparse
import math
import os
import re
import sys
from collections import OrderedDict

NAN = float("nan")


def _f(tok):
    """Fortran writes *** on field overflow, and a killed run leaves half-written lines."""
    try:
        return float(tok)
    except (TypeError, ValueError):
        return NAN


def isnan(x):
    return x is None or (isinstance(x, float) and math.isnan(x))


# ----------------------------------------------------------------------------------------
# Parsing
# ----------------------------------------------------------------------------------------

RE_STEP = re.compile(
    r"\*\s+time step\s*:\s*(\d+)\s+(\d+)\s+(\d+)\s+([-\d.EedD+]+)\s+([-\d.EedD+]+)")
RE_AFTER = re.compile(r"After step\s+(\d+)\s*\(t_now=\s*([-\d.EedD+]+)\s*\)")

RE_HEAD = re.compile(
    r"SHEATH:\s*I_wall=\s*(\S+)\s*A\s*\(Ampere\s+(\S+)\s*A\)\s*"
    r"ePhi/kTe min/mean/max=\s*(\S+)\s*/\s*(\S+)\s*/\s*(\S+)\s*"
    r"max\|j/jsat\|=\s*(\S+)\s*e-limited\s*(\S+)\s*%")
RE_ACTIVE = re.compile(
    r"ePhi/kTe max where the sheath is ACTIVE=\s*(\S+)\s+gated-off area\s+(\S+)\s*%")
RE_NODES = re.compile(r"\|zj-zj_sh\|/\|zj_sat\|\s+nodes=\s*(\S+)\s+gauss=\s*(\S+)")
RE_WEAK = re.compile(r"weak \|F_a/D_a\|/\|S_a/D_a\|=\s*(\S+)\s+over\s+(\d+)\s+trace samples")
RE_TARGET = re.compile(
    r"(INNER|OUTER) target:\s*I_sheath=\s*(\S+)\s*A\s*I_Ampere=\s*(\S+)\s*"
    r"ePhi/kTe=\s*(\S+)\s*max\|j/jsat\|=\s*(\S+)\s*area=\s*(\S+)\s*m\^2")
RE_BND = re.compile(
    r"bnd type\s+(\d+):\s*I_sheath=\s*(\S+)\s*A\s*I_Ampere=\s*(\S+)\s*A\s*area=\s*(\S+)\s*m\^2")
RE_ROWS = re.compile(
    r"sheath trace rows:\s*(\d+)\s+accumulated,\s*(\d+)\s+replaced,\s*(\d+)\s+below the floor")


def closure(a, b):
    """Relative disagreement between two independent measures of the same current."""
    d = max(abs(a), abs(b))
    return abs(a - b) / d if d > 0 else NAN


def derive(r):
    ii, io = r.get("I_inner", NAN), r.get("I_outer", NAN)
    # The thermoelectric signature: antisymmetric (loop) vs symmetric (net) current.
    r["I_loop"] = 0.5 * (ii - io)
    pi, po = r.get("ePhi_inner", NAN), r.get("ePhi_outer", NAN)
    r["dePhi"] = pi - po
    r["clo_inner"] = closure(ii, r.get("IA_inner", NAN))
    r["clo_outer"] = closure(io, r.get("IA_outer", NAN))
    r["clo_wall"] = closure(r.get("I_wall", NAN), r.get("I_wall_A", NAN))


def parse_log(path, keep_all=False):
    """Return (list of records in file order, n_blocks_seen).

    A SHEATH block is attributed to the step named by the most recent 'time step :' banner
    above it, and to the t_now printed by the 'After step' line below it.
    """
    recs = []
    by_step = {}
    cur = None
    step = None
    tstep = NAN
    seen_banner = False
    n_blocks = 0

    def flush():
        nonlocal cur
        if cur is None:
            return
        derive(cur)
        s = cur["step"]
        if (not keep_all) and s in by_step:
            recs[by_step[s]] = cur          # keep the LAST block of a repeated step
        else:
            by_step[s] = len(recs)
            recs.append(cur)
        cur = None

    with open(path, "r", errors="replace") as fh:
        for line in fh:
            m = RE_STEP.search(line)
            if m:
                flush()
                step = int(m.group(2))       # istep
                tstep = _f(m.group(4))
                seen_banner = True
                continue

            m = RE_AFTER.search(line)
            if m:
                # t_now is printed at the END of a step, so it belongs to the block above.
                s = int(m.group(1))
                tgt = None
                if cur is not None and (cur["step"] == s or not seen_banner):
                    tgt = cur
                elif s in by_step:
                    tgt = recs[by_step[s]]
                elif (not seen_banner) and recs:
                    tgt = recs[-1]
                if tgt is not None:
                    tgt["t_now"] = _f(m.group(2))
                    # With no banner the block's "step" is only a counter. The After-step
                    # line carries the real number, so adopt it -- otherwise a log from a
                    # RESTARTED run (counter 1,2,3... vs istep 1000,1001,...) never matches
                    # and t_now is silently dropped, taking time alignment with it.
                    if not seen_banner:
                        tgt["step"] = s
                continue

            m = RE_HEAD.search(line)
            if m:
                flush()
                n_blocks += 1
                cur = {
                    "step": step if seen_banner else n_blocks,
                    "step_is_guess": not seen_banner,
                    "tstep": tstep, "t_now": NAN,
                    "I_wall": _f(m.group(1)), "I_wall_A": _f(m.group(2)),
                    "ePhi_min": _f(m.group(3)), "ePhi_mean": _f(m.group(4)),
                    "ePhi_max": _f(m.group(5)), "jsat_max": _f(m.group(6)),
                    "elim_pct": _f(m.group(7)), "bnd": {},
                }
                continue

            if cur is None:
                continue

            m = RE_ACTIVE.search(line)
            if m:
                cur["ePhi_max_active"] = _f(m.group(1)); cur["gated_pct"] = _f(m.group(2)); continue
            m = RE_NODES.search(line)
            if m:
                cur["nodes_res"] = _f(m.group(1)); cur["gauss"] = _f(m.group(2)); continue
            m = RE_WEAK.search(line)
            if m:
                cur["weak"] = _f(m.group(1)); cur["n_trace_samples"] = int(m.group(2)); continue
            m = RE_TARGET.search(line)
            if m:
                t = m.group(1).lower()
                cur["I_%s" % t] = _f(m.group(2)); cur["IA_%s" % t] = _f(m.group(3))
                cur["ePhi_%s" % t] = _f(m.group(4)); cur["jsat_%s" % t] = _f(m.group(5))
                cur["area_%s" % t] = _f(m.group(6)); continue
            m = RE_BND.search(line)
            if m:
                cur["bnd"][int(m.group(1))] = (_f(m.group(2)), _f(m.group(3)), _f(m.group(4)))
                continue
            m = RE_ROWS.search(line)
            if m:
                cur["rows_acc"] = int(m.group(1)); cur["rows_repl"] = int(m.group(2))
                cur["rows_floor"] = int(m.group(3)); continue

    flush()
    return recs, n_blocks


# ----------------------------------------------------------------------------------------
# Time-weighted statistics
#
# The two runs may sample time at different rates (adaptive tstep, different output
# cadence). A plain mean over records would then weight the finer-stepping run's late
# behaviour more heavily. Every average below is trapezoidal in t_now, so both runs are
# weighted by physical time and are directly comparable.
# ----------------------------------------------------------------------------------------

class Window:
    """A run's records inside one comparison window, with its abscissa."""

    def __init__(self, recs, xkey):
        self.recs = recs
        self.xkey = xkey

    def __len__(self):
        return len(self.recs)

    def xy(self, key):
        out = []
        for r in self.recs:
            x, y = r.get(self.xkey, NAN), r.get(key, NAN)
            if not isnan(x) and not isnan(y):
                out.append((x, y))
        out.sort()
        return out

    def span(self):
        xs = [r.get(self.xkey, NAN) for r in self.recs]
        xs = [x for x in xs if not isnan(x)]
        return (min(xs), max(xs)) if xs else (NAN, NAN)


def _weights(pts):
    """Trapezoidal weights in x for a sorted list of (x, y)."""
    n = len(pts)
    if n == 0:
        return []
    if n == 1:
        return [1.0]
    w = [0.0] * n
    for i in range(n - 1):
        dx = pts[i + 1][0] - pts[i][0]
        if dx <= 0:
            continue
        w[i] += 0.5 * dx
        w[i + 1] += 0.5 * dx
    if sum(w) <= 0:                       # degenerate (all same x): fall back to uniform
        w = [1.0] * n
    return w


def wmean(win, key):
    pts = win.xy(key)
    if not pts:
        return NAN
    w = _weights(pts)
    tw = sum(w)
    return sum(wi * y for wi, (_, y) in zip(w, pts)) / tw if tw > 0 else NAN


def wstd(win, key):
    pts = win.xy(key)
    if len(pts) < 2:
        return NAN
    w = _weights(pts)
    tw = sum(w)
    if tw <= 0:
        return NAN
    m = sum(wi * y for wi, (_, y) in zip(w, pts)) / tw
    v = sum(wi * (y - m) ** 2 for wi, (_, y) in zip(w, pts)) / tw
    # small-sample correction so short windows do not look artificially quiet
    n = len(pts)
    return math.sqrt(v * n / (n - 1))


def wdrift(win, key):
    """Change across the window: mean(second half) - mean(first half), split at the
    MIDPOINT IN X, not by record count -- otherwise a run whose steps shrink late puts
    most of its records in the second half and the split is not a split in time.

    Returned in the metric's own units, deliberately NOT as a number of sigma: a trend
    inflates the very standard deviation you would divide by, so it would hide itself.
    It is judged against the A/B delta it competes with -- see verdict_for().
    """
    pts = win.xy(key)
    if len(pts) < 8:
        return NAN
    xlo, xhi = pts[0][0], pts[-1][0]
    if xhi <= xlo:
        return NAN
    mid = 0.5 * (xlo + xhi)
    a = [p for p in pts if p[0] <= mid]
    b = [p for p in pts if p[0] > mid]
    if len(a) < 2 or len(b) < 2:
        return NAN

    def m(sub):
        w = _weights(sub)
        tw = sum(w)
        return sum(wi * y for wi, (_, y) in zip(w, sub)) / tw if tw > 0 else NAN

    return m(b) - m(a)


def verdict_for(delta, pooled_sd, drift_a, drift_b):
    """Classify one metric's A/B difference. Returns (tag, worst_drift)."""
    dr = [abs(x) for x in (drift_a, drift_b) if not isnan(x)]
    worst = max(dr) if dr else NAN
    if isnan(delta):
        return "?", worst
    # Drift is checked BEFORE scatter. While a run is still relaxing, its trend inflates
    # its own standard deviation, so a "within noise" test would fire for the wrong reason
    # and report a transient as a null result.
    if not isnan(worst) and worst > abs(delta):
        return "DRIFT > delta", worst
    if not isnan(pooled_sd) and pooled_sd > 0 and abs(delta) < pooled_sd:
        return "within noise", worst
    if isnan(worst):
        return "?", worst
    if worst > 0.5 * abs(delta):
        return "drift ~ delta", worst
    return "ok", worst


# ----------------------------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------------------------

def fmt(x, w=11, p=4):
    if isnan(x):
        return " " * max(0, w - 3) + "n/a"
    if x == 0:
        return "%*.1f" % (w, 0.0)
    a = abs(x)
    if a >= 1e4 or a < 1e-3:
        return "%*.*e" % (w, p - 1, x)
    return "%*.*f" % (w, p, x)


def pct(x):
    return "     n/a" if isnan(x) else "%+7.2f%%" % x


class Row:
    def __init__(self, name, key, unit="", note=""):
        self.name, self.key, self.unit, self.note = name, key, unit, note


MAIN_ROWS = [
    Row("I_loop = (I_in-I_out)/2", "I_loop", "A", "THE thermoelectric observable"),
    Row("I_inner", "I_inner", "A"),
    Row("I_outer", "I_outer", "A"),
    Row("I_wall (net)", "I_wall", "A", "set by quasi-neutrality, expect ~no change"),
    Row("ePhi/kTe inner", "ePhi_inner", ""),
    Row("ePhi/kTe outer", "ePhi_outer", ""),
    Row("d(ePhi/kTe) in-out", "dePhi", "", "in units of Te, NOT volts"),
    Row("ePhi/kTe mean", "ePhi_mean", ""),
    Row("ePhi/kTe max", "ePhi_max", ""),
    Row("max|j/jsat|", "jsat_max", "", "saturation headroom"),
    Row("e-limited area", "elim_pct", "%"),
    Row("gated-off area", "gated_pct", "%"),
]

HEALTH_ROWS = [
    Row("weak residual", "weak", "", "BC convergence"),
    Row("gauss residual", "gauss", ""),
    Row("closure INNER", "clo_inner", "", "|I_sheath-I_Ampere|/max"),
    Row("closure OUTER", "clo_outer", ""),
    Row("closure wall", "clo_wall", ""),
]

W = 112


def compare_table(title, rows, wa, wb, la, lb, out):
    out.append("")
    out.append(title)
    out.append("-" * W)
    out.append("  %-24s %12s %12s %12s %9s %11s   %s"
               % ("metric", la[:12], lb[:12], "delta", "rel", "drift", "verdict"))
    out.append("-" * W)
    for row in rows:
        ma, mb = wmean(wa, row.key), wmean(wb, row.key)
        d = mb - ma
        rel = 100.0 * d / abs(ma) if (not isnan(ma) and ma != 0) else NAN
        sa, sb = wstd(wa, row.key), wstd(wb, row.key)
        pooled = math.sqrt((sa ** 2 + sb ** 2) / 2) if not (isnan(sa) or isnan(sb)) else NAN
        tag, worst = verdict_for(d, pooled, wdrift(wa, row.key), wdrift(wb, row.key))
        out.append("  %-24s %12s %12s %12s %9s %11s   %s"
                   % (row.name + (" [%s]" % row.unit if row.unit else ""),
                      fmt(ma, 12), fmt(mb, 12), fmt(d, 12), pct(rel), fmt(worst, 11), tag))
        if row.note:
            out.append("  %-24s ( %s )" % ("", row.note))
    out.append("-" * W)


def _wrap(s, w):
    words, lines, cur = s.split(), [], ""
    for word in words:
        if len(cur) + len(word) + 1 > w:
            lines.append(cur); cur = word
        else:
            cur = (cur + " " + word).strip()
    if cur:
        lines.append(cur)
    return lines


# ----------------------------------------------------------------------------------------
# Window selection
# ----------------------------------------------------------------------------------------

def choose_window(ra, rb, args, out):
    """Return (win_a, win_b, xkey, description) or (None, None, None, reason)."""
    ta = [r["t_now"] for r in ra if not isnan(r.get("t_now", NAN))]
    tb = [r["t_now"] for r in rb if not isnan(r.get("t_now", NAN))]
    have_time = (len(ta) >= 0.5 * len(ra)) and (len(tb) >= 0.5 * len(rb)) and ta and tb

    mode = args.align
    if mode == "auto":
        mode = "time" if have_time else "step"
    if mode == "time" and not have_time:
        out.append("  WARNING: --align time requested but t_now is missing from one of the logs")
        out.append("           ('After step NNNNNN (t_now= ...)'). Falling back to step alignment.")
        mode = "step"

    if mode == "step":
        out.append("  ALIGNMENT: by STEP NUMBER.")
        out.append("             tstep is adaptive in JOREK, so equal step numbers are NOT equal")
        out.append("             times. If the two runs diverged, this compares different times.")
        sa = {r["step"]: r for r in ra}
        sb = {r["step"]: r for r in rb}
        common = sorted(set(sa) & set(sb))
        if not common:
            return None, None, None, "the two runs share no common step number"
        if args.s_from is not None or args.s_to is not None:
            lo = args.s_from if args.s_from is not None else common[0]
            hi = args.s_to if args.s_to is not None else common[-1]
            sel = [s for s in common if lo <= s <= hi]
            how = "steps %g..%g" % (lo, hi)
        else:
            sel = common[-args.window:]
            how = "last %d common steps" % len(sel)
        if not sel:
            return None, None, None, "the requested step range is empty"
        out.append("  common steps : %d  (%s .. %s)" % (len(common), common[0], common[-1]))
        return (Window([sa[s] for s in sel], "step"),
                Window([sb[s] for s in sel], "step"), "step", how)

    # --- time alignment -----------------------------------------------------------------
    lo_ov, hi_ov = max(min(ta), min(tb)), min(max(ta), max(tb))
    out.append("  ALIGNMENT: by PHYSICAL TIME t_now (correct under adaptive tstep).")
    out.append("  time range : %-14s t_now %s .. %s   (%d steps)"
               % (args.labels[0][:14] if args.labels else "run A",
                  fmt(min(ta)).strip(), fmt(max(ta)).strip(), len(ta)))
    out.append("  time range : %-14s t_now %s .. %s   (%d steps)"
               % (args.labels[1][:14] if args.labels else "run B",
                  fmt(min(tb)).strip(), fmt(max(tb)).strip(), len(tb)))
    if hi_ov <= lo_ov:
        return None, None, None, ("the two runs do not overlap in time (A ends at %s, B starts "
                                  "at %s)" % (fmt(max(ta)).strip(), fmt(min(tb)).strip()))
    out.append("  overlap    : t_now %s .. %s" % (fmt(lo_ov).strip(), fmt(hi_ov).strip()))

    if args.t_from is not None or args.t_to is not None:
        lo = args.t_from if args.t_from is not None else lo_ov
        hi = args.t_to if args.t_to is not None else hi_ov
        if lo < lo_ov or hi > hi_ov:
            out.append("  NOTE: the requested time window extends outside the overlap; clipped.")
        lo, hi = max(lo, lo_ov), min(hi, hi_ov)
        how = "t_now %s..%s (requested)" % (fmt(lo).strip(), fmt(hi).strip())
    else:
        # Count back --window steps in whichever run has FEWER samples in the overlap, and
        # use that run's t_now as the window start. Both runs are then cut at the same time.
        ina = sorted(t for t in ta if lo_ov <= t <= hi_ov)
        inb = sorted(t for t in tb if lo_ov <= t <= hi_ov)
        scarcer = ina if len(ina) <= len(inb) else inb
        lo = scarcer[-args.window] if len(scarcer) > args.window else lo_ov
        hi = hi_ov
        how = "last %d steps of the sparser run -> t_now %s..%s" % (
            args.window, fmt(lo).strip(), fmt(hi).strip())

    sel_a = [r for r in ra if not isnan(r.get("t_now", NAN)) and lo <= r["t_now"] <= hi]
    sel_b = [r for r in rb if not isnan(r.get("t_now", NAN)) and lo <= r["t_now"] <= hi]
    if len(sel_a) < 2 or len(sel_b) < 2:
        return None, None, None, ("fewer than 2 samples in the window for one of the runs "
                                  "(A: %d, B: %d)" % (len(sel_a), len(sel_b)))
    return Window(sel_a, "t_now"), Window(sel_b, "t_now"), "t_now", how


# ----------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="A/B comparison of two JOREK runs from their SHEATH log output.")
    ap.add_argument("ref", help="reference logfile (e.g. thermoelectric_ohm = .false.)")
    ap.add_argument("test", help="test logfile (e.g. thermoelectric_ohm = .true.)")
    ap.add_argument("--labels", nargs=2, metavar=("A", "B"), default=None)
    ap.add_argument("--window", type=int, default=200,
                    help="average over the last N steps of the sparser run (default 200)")
    ap.add_argument("--from-time", dest="t_from", type=float, default=None)
    ap.add_argument("--to-time", dest="t_to", type=float, default=None)
    ap.add_argument("--from-step", dest="s_from", type=int, default=None)
    ap.add_argument("--to-step", dest="s_to", type=int, default=None)
    ap.add_argument("--align", choices=("auto", "time", "step"), default="auto")
    ap.add_argument("--weak-tol", type=float, default=1e-2)
    ap.add_argument("--plot", default=None, metavar="FILE.png")
    ap.add_argument("--csv", default=None, metavar="PREFIX")
    ap.add_argument("--all-blocks", action="store_true")
    args = ap.parse_args()

    if not args.labels:
        args.labels = [os.path.basename(args.ref), os.path.basename(args.test)]
    la, lb = args.labels

    ra, na = parse_log(args.ref, args.all_blocks)
    rb, nb = parse_log(args.test, args.all_blocks)

    out = ["=" * W, "SHEATH A/B  --  %s   vs   %s" % (la, lb), "=" * W]

    if not ra or not rb:
        print("\n".join(out))
        for lab, r, p in ((la, ra, args.ref), (lb, rb, args.test)):
            if not r:
                print("\nERROR: no SHEATH blocks found in %s (%s)." % (p, lab))
        print("\nThe SHEATH block is written by rank 0 only. If this is a per-rank log, use the")
        print("rank-0 file. Check also that the sheath BC was actually enabled in the namelist.")
        return 1

    for lab, r, n in ((la, ra, na), (lb, rb, nb)):
        out.append("  %-14s %6d blocks, %6d kept, steps %s .. %s"
                   % (lab, n, len(r), r[0]["step"], r[-1]["step"]))
    if na != len(ra) or nb != len(rb):
        out.append("  (construct_matrix printed more than once per step; kept the last."
                   " Use --all-blocks to keep all.)")
    if any(r.get("step_is_guess") for r in ra + rb):
        out.append("  WARNING: no 'time step :' banner found -- steps are block counters.")

    # Where does each run stand right now? The window mean can be dominated by a
    # transient the run has already left, so the last sample is reported separately.
    out.append("")
    out.append("  CURRENT STATE (last sample in each log)")
    for lab, r in ((la, ra), (lb, rb)):
        last = r[-1]
        out.append("    %-14s t_now %-10s weak %-11s I_wall %-11s ePhi/kTe mean %s"
                   % (lab, fmt(last.get("t_now", NAN)).strip(),
                      fmt(last.get("weak", NAN)).strip(),
                      fmt(last.get("I_wall", NAN)).strip(),
                      fmt(last.get("ePhi_mean", NAN)).strip()))
    out.append("")

    wa, wb, xkey, how = choose_window(ra, rb, args, out)
    if wa is None:
        out += ["", "ERROR: " + how + ".", "",
                "       An A/B needs a shared interval. Comparing non-overlapping stretches of",
                "       two runs measures when you looked, not what you changed."]
        print("\n".join(out))
        return 1

    xa, xb = wa.span(), wb.span()

    def sp(v):
        return ("%d" % v) if xkey == "step" and not isnan(v) else fmt(v).strip()

    out.append("  window used: %s" % how)
    out.append("               %-14s %4d samples, %s %s .. %s"
               % (la[:14], len(wa), xkey, sp(xa[0]), sp(xa[1])))
    out.append("               %-14s %4d samples, %s %s .. %s"
               % (lb[:14], len(wb), xkey, sp(xb[0]), sp(xb[1])))
    if xkey == "t_now":
        full_a = ra[-1].get("t_now", NAN) - ra[0].get("t_now", NAN)
        if not isnan(full_a) and full_a > 0 and (xa[1] - xa[0]) > 0.9 * full_a:
            out.append("  WARNING: the window spans essentially the WHOLE of %s, startup" % la)
            out.append("           included. Averages over a relaxation are not steady-state")
            out.append("           values. Use --from-time to cut the transient off.")
    if min(len(wa), len(wb)) < 8:
        out.append("  WARNING: fewer than 8 samples in a window -- no drift estimate is possible.")
    elif min(len(wa), len(wb)) < 20:
        out.append("  WARNING: fewer than 20 samples -- the scatter estimate is weak.")
    if len(wa) and len(wb):
        ratio = max(len(wa), len(wb)) / min(len(wa), len(wb))
        if ratio > 2.0:
            out.append("  NOTE: one run has %.1fx more samples in the same time window (different"
                       % ratio)
            out.append("        tstep or output cadence). Averages are time-weighted, so this is")
            out.append("        handled -- but the sparser run resolves the scatter less well.")

    # --- gate 1: convergence -------------------------------------------------------------
    out.append("")
    out.append("GATE 1 -- BC convergence (weak residual, tol %.1e)" % args.weak_tol)
    bad = []
    for lab, w in ((la, wa), (lb, wb)):
        v = wmean(w, "weak")
        ok = (not isnan(v)) and v < args.weak_tol
        out.append("  %-14s weak = %s   %s" % (lab, fmt(v), "ok" if ok else "NOT CONVERGED"))
        if not ok:
            bad.append(lab)
    if bad:
        out.append("  ==> %s has not converged the sheath BC in this window. The physics"
                   % " and ".join(bad))
        out.append("      difference below is not separable from BC error.")

    # --- gate 2: same constraint ---------------------------------------------------------
    out.append("")
    out.append("GATE 2 -- trace rows (both runs must impose the SAME constraint)")
    reps = []
    for lab, w in ((la, wa), (lb, wb)):
        acc = [r.get("rows_acc", NAN) for r in w.recs]
        acc = [a for a in acc if not isnan(a)]
        if acc:
            rep = wmean(w, "rows_repl")
            reps.append(rep)
            out.append("  %-14s accumulated %g, replaced %g, below floor %g   %s"
                       % (lab, wmean(w, "rows_acc"), rep, wmean(w, "rows_floor"),
                          "constant" if min(acc) == max(acc) else "VARYING across the window"))
        else:
            reps.append(NAN)
            out.append("  %-14s no trace-row line (nodal route, or an older build)" % lab)
    if len(reps) == 2 and not any(isnan(x) for x in reps) and abs(reps[0] - reps[1]) > 0.5:
        out.append("  ==> The runs replace DIFFERENT numbers of rows (%g vs %g); they are not"
                   % (reps[0], reps[1]))
        out.append("      imposing the same BC, so this is not a controlled comparison.")

    # --- gate 3: same geometry -----------------------------------------------------------
    for k, nm in (("area_inner", "INNER area"), ("area_outer", "OUTER area")):
        aa, ab = wmean(wa, k), wmean(wb, k)
        if not (isnan(aa) or isnan(ab)) and aa > 0 and abs(ab - aa) / aa > 1e-6:
            out.append("  ==> %s differs: %s vs %s -- different grids, not an A/B."
                       % (nm, fmt(aa).strip(), fmt(ab).strip()))

    compare_table("PHYSICS  (time-weighted window mean; 'delta' is %s minus %s)" % (lb, la),
                  MAIN_ROWS, wa, wb, la, lb, out)
    compare_table("HEALTH  (both columns should be small and similar)",
                  HEALTH_ROWS, wa, wb, la, lb, out)

    # --- per boundary type ---------------------------------------------------------------
    types = sorted(set().union(*[set(r["bnd"]) for r in wa.recs + wb.recs]) or set())
    if types:
        out += ["", "PER BOUNDARY TYPE  (I_sheath, time-weighted window mean)", "-" * W,
                "  %-24s %12s %12s %12s %9s" % ("", la[:12], lb[:12], "delta", "rel")]
        for t in types:
            for w in (wa, wb):
                for r in w.recs:
                    r["_bnd%d" % t] = r["bnd"].get(t, (NAN,))[0]
            ma, mb = wmean(wa, "_bnd%d" % t), wmean(wb, "_bnd%d" % t)
            d = mb - ma
            rel = 100.0 * d / abs(ma) if (not isnan(ma) and ma != 0) else NAN
            out.append("  %-24s %12s %12s %12s %9s"
                       % ("bnd type %d [A]" % t, fmt(ma, 12), fmt(mb, 12), fmt(d, 12), pct(rel)))
        out.append("-" * W)

    # --- verdict -------------------------------------------------------------------------
    out += ["", "READING THIS", "-" * W]
    ma, mb = wmean(wa, "I_loop"), wmean(wb, "I_loop")
    d = mb - ma
    sa, sb = wstd(wa, "I_loop"), wstd(wb, "I_loop")
    pooled = math.sqrt((sa ** 2 + sb ** 2) / 2) if not (isnan(sa) or isnan(sb)) else NAN
    dra, drb = wdrift(wa, "I_loop"), wdrift(wb, "I_loop")
    tag, worst = verdict_for(d, pooled, dra, drb)

    out += ["  I_loop is the antisymmetric target current -- what the 0.71*grad_par(Te) thermal",
            "  force acts on directly. I_wall is the symmetric part and is pinned by",
            "  quasi-neutrality, so an effect appearing ONLY in I_wall is more likely a transport",
            "  change than a thermoelectric one.",
            "",
            "    I_loop  %-14s = %s A" % (la, fmt(ma)),
            "    I_loop  %-14s = %s A" % (lb, fmt(mb)),
            "    delta                   = %s A   (%s)"
            % (fmt(d), pct(100 * d / abs(ma) if (not isnan(ma) and ma) else NAN)),
            "    run-to-run scatter      = %s A" % fmt(pooled),
            "    within-run drift        = %s A (%s) / %s A (%s)"
            % (fmt(dra), la[:12], fmt(drb), lb[:12])]

    if bad:
        v = "BLOCKED: %s not converged (gate 1). Nothing above is trustworthy." % " and ".join(bad)
    elif tag == "DRIFT > delta":
        v = ("BLOCKED: I_loop drifts by %s A across the window, MORE than the %s A difference "
             "between the runs. These are two transients caught at different phases, not a "
             "measured effect. Extend the runs, or move the window later and re-check."
             % (fmt(worst).strip(), fmt(d).strip()))
    elif tag == "within noise":
        v = ("NEGATIVE: the shift in I_loop (%s A) is smaller than the runs' own scatter (%s A). "
             "The thermoelectric term is not changing the target current balance at this level."
             % (fmt(d).strip(), fmt(pooled).strip()))
    elif tag == "drift ~ delta":
        v = ("MARGINAL: the difference is %s A but each run still drifts by up to %s A inside the "
             "window. Real, but do not quote the magnitude yet."
             % (fmt(d).strip(), fmt(worst).strip()))
    elif tag == "ok":
        v = ("POSITIVE: I_loop shifts by %s A (%.0fx the scatter), against a residual within-run "
             "drift of %s A. The thermoelectric term is moving the loop current."
             % (fmt(d).strip(), abs(d) / pooled if pooled else NAN, fmt(worst).strip()))
    else:
        v = "Cannot classify: too few samples in the window for a drift estimate (need 8)."
    for i, chunk in enumerate(_wrap(v, W - 6)):
        out.append("  " + ("  " if i else "") + chunk)

    out += ["",
            "  Note: d(ePhi/kTe) is in units of the LOCAL Te, so an unchanged d(ePhi/kTe) across",
            "  a run with a different Te profile is NOT an unchanged potential in volts. The",
            "  SHEATH line carries no Te -- take that from the target profiles instead.",
            "-" * W]

    print("\n".join(out))

    if args.csv:
        for lab, r, suf in ((la, ra, "a"), (lb, rb, "b")):
            write_csv("%s_%s.csv" % (args.csv, suf), r)
            print("wrote %s_%s.csv  (%s)" % (args.csv, suf, lab))
    if args.plot:
        make_plot(args.plot, ra, rb, la, lb, wa, wb, xkey)
    return 0


CSV_KEYS = ["step", "t_now", "tstep", "I_wall", "I_wall_A", "I_inner", "I_outer", "I_loop",
            "ePhi_inner", "ePhi_outer", "dePhi", "ePhi_mean", "ePhi_max", "jsat_max",
            "elim_pct", "gated_pct", "weak", "gauss", "clo_inner", "clo_outer", "clo_wall",
            "rows_acc", "rows_repl", "rows_floor"]


def write_csv(path, recs):
    with open(path, "w") as fh:
        fh.write(",".join(CSV_KEYS) + "\n")
        for r in recs:
            fh.write(",".join("" if isnan(r.get(k, NAN)) else repr(r[k]) for k in CSV_KEYS) + "\n")


def make_plot(path, ra, rb, la, lb, wa, wb, xkey):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not available -- skipping --plot)")
        return
    panels = [("I_loop", "I_loop = (I_in-I_out)/2  [A]"),
              ("I_wall", "net wall current  [A]"),
              ("dePhi", "ePhi/kTe  inner - outer"),
              ("ePhi_max", "max ePhi/kTe"),
              ("jsat_max", "max |j/jsat|"),
              ("weak", "weak residual  (log)")]
    fig, axes = plt.subplots(3, 2, figsize=(11, 9), sharex=True)
    for ax, (key, title) in zip(axes.ravel(), panels):
        for recs, lab, col in ((ra, la, "#3B6EA5"), (rb, lb, "#B5482A")):
            pts = sorted((r.get(xkey, NAN), r.get(key, NAN)) for r in recs)
            pts = [(x, y) for x, y in pts if not isnan(x) and not isnan(y)]
            if pts:
                ax.plot([p[0] for p in pts], [p[1] for p in pts], lw=1.2, color=col, label=lab)
        lo = min(wa.span()[0], wb.span()[0])
        hi = max(wa.span()[1], wb.span()[1])
        if not (isnan(lo) or isnan(hi)):
            ax.axvspan(lo, hi, color="0.85", zorder=0)
        ax.set_title(title, fontsize=9)
        ax.grid(alpha=0.3, lw=0.5)
        if key == "weak":
            ax.set_yscale("log")
    axes[0, 0].legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel(xkey)
    fig.suptitle("SHEATH A/B: %s vs %s   (shaded = averaging window, aligned on %s)"
                 % (la, lb, xkey), fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(path, dpi=130)
    print("\nwrote %s" % path)


if __name__ == "__main__":
    sys.exit(main())
