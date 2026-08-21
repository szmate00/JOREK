#!/usr/bin/env python3
"""Fit the exponential growth rate of a boundary mode from max|w| samples.

Why this rather than the crash cycle: if the mode grows as exp(gamma t) and the run ends when
it reaches some fixed amplitude, then

    t_crash = (1/gamma) * ln(A_crit/A_0),

so the crash time depends only LOGARITHMICALLY on the seed. Halving the seed moves it by
0.69/gamma - a few percent - while a 20% change in gamma moves it by 20%. Comparing crash
cycles therefore hides exactly the quantity a formulation change is meant to alter, and
comparing growth rates exposes it.

It is also much cheaper: gamma is measurable a few hundred cycles in, long before any crash.

Input: whitespace-separated "time max_w" pairs (a third column is ignored), from a file or
stdin. Time in whatever unit you like; gamma comes back in its inverse.

    ./growth_rate.py samples.txt [--tau_a 6.53e-7] [--label "gauge removal"]
"""
import sys
import numpy as np

def main():
    args = sys.argv[1:]
    tau_a, label = None, None
    if '--tau_a' in args:
        i = args.index('--tau_a'); tau_a = float(args[i+1]); args = args[:i]+args[i+2:]
    if '--label' in args:
        i = args.index('--label'); label = args[i+1]; args = args[:i]+args[i+2:]
    src = open(args[0]) if args else sys.stdin
    rows = [l.split() for l in src if l.strip() and not l.lstrip().startswith('#')]
    d = np.array([[float(r[0]), float(r[1])] for r in rows])
    if len(d) < 3:
        sys.exit("need at least 3 samples")
    t, w = d[:,0], np.abs(d[:,1])
    if (w <= 0).any():
        sys.exit("max|w| must be positive")

    # fit ln w = ln A0 + gamma t
    A = np.vstack([np.ones_like(t), t]).T
    coef, res, *_ = np.linalg.lstsq(A, np.log(w), rcond=None)
    lnA0, gam = coef
    pred = A @ coef
    ss = 1.0 - np.sum((np.log(w)-pred)**2)/max(np.sum((np.log(w)-np.mean(np.log(w)))**2), 1e-300)

    if label: print(f"  {label}")
    print(f"  samples          : {len(t)}  over t = {t.min():.4g} .. {t.max():.4g}")
    print(f"  growth rate gamma: {gam:.6g} per unit time")
    print(f"  e-folding time   : {1/gam:.6g}" + (f" = {1/gam/tau_a:.0f} tau_A" if tau_a else ""))
    print(f"  fitted seed A0   : {np.exp(lnA0):.4g}")
    print(f"  R^2 of ln-fit    : {ss:.4f}" + ("" if ss > 0.9 else "   <-- POOR: not a clean exponential"))
    print()
    print("  Compare gamma between configurations, not crash cycles: t_crash depends on the")
    print("  seed only as ln(A_crit/A0), so a 2x seed change moves it by only 0.69/gamma.")
    if tau_a:
        for r in (2.0, 5.0, 10.0):
            print(f"    a {r:.0f}x seed change would shift the crash by {np.log(r)/gam/tau_a:.0f} tau_A")

if __name__ == '__main__':
    main()
