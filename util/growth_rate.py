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

Two input modes:

  --log FILE    parse a JOREK logfile directly. JOREK already writes, every step,

                    After step 000705 (t_now= 1.421E+03):
                    W_mag,_kin      =  ... ...,  W_kin(1) ... W_kin(n_tor),

                W_kin(1) is the n=0 kinetic energy - the ExB energy the boundary structure
                lives in. It is a global integral, so far less noisy than a field maximum,
                and it needs no VTK post-processing. Energy scales as amplitude squared, so
                the amplitude growth rate is HALF the fitted rate; that factor is applied.

  positional    whitespace-separated "time value" pairs from a file or stdin, for a field
                maximum read off by hand. Add --energy if the values are energies.

    ./growth_rate.py --log jorek.log [--tau_a 6.53e-7] [--label "natural w"]
    ./growth_rate.py samples.txt --tau_a 6.53e-7
"""
import sys, re
import numpy as np

def parse_jorek_log(path):
    """Pull (t_now, W_kin(1)) from a JOREK logfile.

    Formats being matched, from mod_jorek_timestepping:
        130 format(1x,a,i6.6,a,es10.3,a)                     -> " After step NNNNNN (t_now= X):"
        131 format(1x,a,2(2(es10.2,' ...',es10.2,',')))      -> " W_mag,_kin = a ... b, c ... d,"
    so W_kin(1) is the THIRD number on the W_mag,_kin line.
    """
    num = r'[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?'
    out, t = [], None
    for line in open(path, errors='ignore'):
        m = re.search(r't_now\s*=\s*(' + num + r')', line)
        if m:
            t = float(m.group(1).replace('D', 'E').replace('d', 'e'))
            continue
        if 'W_mag' in line and '_kin' in line and t is not None:
            vals = [float(v.replace('D', 'E').replace('d', 'e'))
                    for v in re.findall(num, line.split('=', 1)[-1])]
            if len(vals) >= 3 and vals[2] > 0.0:
                out.append((t, vals[2]))
            t = None
    if not out:
        sys.exit("no 'After step ... t_now' / 'W_mag,_kin' pairs found - is this a JOREK logfile?")
    return out

def main():
    args = sys.argv[1:]
    tau_a, label = None, None
    is_energy = False
    logfile = None
    if '--label' in args:
        i = args.index('--label'); label = args[i+1]; args = args[:i]+args[i+2:]
    if '--log' in args:
        i = args.index('--log'); logfile = args[i+1]; args = args[:i]+args[i+2:]
        is_energy = True
        # JOREK's t_now is already in Alfven times, so no conversion is wanted here.
        # Passing --tau_a with --log would divide an Alfven time by seconds.
        tau_a = 1.0
    if '--energy' in args:
        i = args.index('--energy'); is_energy = True; args = args[:i]+args[i+1:]
    if '--tau_a' in args:
        i = args.index('--tau_a')
        if logfile is None:
            tau_a = float(args[i+1])
        else:
            print("  note: --tau_a ignored with --log; t_now is already in Alfven times")
        args = args[:i]+args[i+2:]
    if '--label' in args:
        i = args.index('--label'); label = args[i+1]; args = args[:i]+args[i+2:]
    if logfile:
        d = np.array(parse_jorek_log(logfile))
        print(f"  parsed {len(d)} steps from {logfile}")
        # drop the early transient: fit the last two thirds
        if len(d) > 30:
            d = d[len(d)//3:]
            print(f"  fitting the last {len(d)} (early transient dropped)")
    else:
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
    if is_energy:
        gam = 0.5 * gam        # energy ~ amplitude^2
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
