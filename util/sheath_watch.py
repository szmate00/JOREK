#!/usr/bin/env python3
"""Read the SHEATH: diagnostic out of a JOREK log and say whether the sheath
boundary condition is converging or running away.

    python3 util/sheath_watch.py run.log [run2.log ...]

The decisive quantity is |I_Ampere - I_wall|: I_wall is the current the sheath
characteristic asks for, I_Ampere is what the plasma actually delivers. With the
correct sheath_flux_sign the gap SHRINKS. With the wrong one it GROWS - that is an
anti-damped boundary condition, and no amount of tuning will save it.
"""
import re, sys

PAT = re.compile(
    r"SHEATH:\s*I_wall=\s*([-\d.Ee+]+)\s*A\s*\(Ampere\s*([-\d.Ee+]+)\s*A\)"
    r"\s*ePhi/kTe min/mean/max=\s*([-\d.]+)\s*/\s*([-\d.]+)\s*/\s*([-\d.]+)"
    r"\s*max\|j/jsat\|=\s*([-\d.Ee+]+)\s*e-limited\s*([-\d.]+)\s*%")


def read(path):
    rows = []
    with open(path, errors="replace") as fh:
        for line in fh:
            m = PAT.search(line)
            if m:
                iw, ia, xmin, xmean, xmax, jr, lim = (float(g) for g in m.groups())
                rows.append(dict(I_wall=iw, I_amp=ia, gap=abs(ia - iw),
                                 phi_min=xmin, phi_mean=xmean, phi_max=xmax,
                                 jratio=jr, lim=lim))
    return rows


def trend(vals, frac=0.25):
    """Compare the mean of the last quarter with the mean of the first quarter."""
    if len(vals) < 8:
        return None
    n = max(2, int(len(vals) * frac))
    a = sum(vals[:n]) / n
    b = sum(vals[-n:]) / n
    if a == 0:
        return None
    return b / a


def report(path):
    rows = read(path)
    print(f"\n=== {path} ===")
    if not rows:
        print("  no SHEATH: lines found - is the diagnostic active "
              "(bcs(:)%natural%u or %sheath_u)?")
        return
    print(f"  {len(rows)} records\n")
    hdr = f"  {'#':>5} {'I_wall[A]':>12} {'I_Ampere[A]':>12} {'|gap|':>12} " \
          f"{'ePhi/kTe mean':>14} {'max|j/jsat|':>12} {'e-lim%':>7}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    idx = list(range(len(rows)))
    show = idx if len(idx) <= 14 else idx[:5] + [None] + idx[-8:]
    for i in show:
        if i is None:
            print(f"  {'...':>5}")
            continue
        r = rows[i]
        print(f"  {i+1:>5} {r['I_wall']:>12.4e} {r['I_amp']:>12.4e} {r['gap']:>12.4e} "
              f"{r['phi_mean']:>14.2f} {r['jratio']:>12.3e} {r['lim']:>7.1f}")

    g = trend([r["gap"] for r in rows])
    print()
    if g is None:
        print("  too few records to judge a trend")
    elif g < 0.7:
        print(f"  GAP SHRINKING  (last/first = {g:.2f})  -> sign is RIGHT, BC is damping")
    elif g > 1.5:
        print(f"  GAP GROWING    (last/first = {g:.2f})  -> sign is likely WRONG "
              f"(flip sheath_flux_sign), or the BC is anti-damped")
    else:
        print(f"  gap roughly flat (last/first = {g:.2f}) - run longer, or the ramp "
              f"has not bitten yet")

    ratio = [r["I_amp"] / r["I_wall"] for r in rows[-8:] if r["I_wall"] != 0]
    if len(ratio) >= 4:
        m = sum(ratio) / len(ratio)
        spread = max(ratio) - min(ratio)
        if abs(spread) < 0.1 * max(abs(m), 1e-30) and abs(m - 1.0) > 0.15:
            print(f"  NOTE I_Ampere/I_wall settling at {m:.3f}, not 1 -> suspect a missing"
                  f" constant factor (R, 1/R or F0) in the surface term")

    lim = [r["lim"] for r in rows]
    t = trend(lim)
    if t is not None and t > 1.5 and lim[-1] > 5.0:
        print(f"  WARNING e-limited fraction climbing ({lim[0]:.1f}% -> {lim[-1]:.1f}%):"
              f" the electron-saturation cap is becoming load-bearing")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    for p in sys.argv[1:]:
        report(p)
