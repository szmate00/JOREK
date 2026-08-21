#!/usr/bin/env python3
"""Analyse the bc_geom_diag_*.dat files written by model600's boundary_conditions.

Each row is one visit to a boundary node:

    R  Z  bnd_type  iv_dir  bn  direction  factor

A node at a corner is visited once from each adjacent boundary element, and each
visit computes its own outward normal. direction = sign(B.n) and the Chodura
factor follow from that normal, so the two visits can disagree - most easily
where B.n passes through zero, which is what a corner between a target (B
near-normal) and a side wall (B grazing) produces. Because
boundary_conditions_add_one_entry overwrites rather than accumulates, a
disagreement means v_par = +c_s or -c_s is decided by element ordering.

The element loop is distributed, so the two visits may land on different ranks.
Merging every rank's file here is what makes those corners visible.

Usage:
    ./bc_geom_analyze.py [files ...] [--region Rmin Rmax Zmin Zmax] [--floor VALUE]
"""
import sys, glob, warnings
import numpy as np

def load(paths):
    # Ranks that own no boundary element write a header and nothing else, which is normal and
    # not worth a warning - only every file being empty means something went wrong.
    rows, empty = [], 0
    for p in paths:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            d = np.loadtxt(p, comments='#', ndmin=2)
        if d.size:
            rows.append(d)
        else:
            empty += 1
    if empty:
        print(f"note: {empty} of {len(paths)} rank files hold no rows "
              f"(ranks owning no boundary element)")
    if not rows:
        sys.exit("every file is header-only: the dump ran but wrote no rows. Check the log for\n"
                 "  'BCDIAG [boundary_conditions]: call N, mach1 boundary-node visits ...'\n"
                 "A visit count > 0 with no rows means the build predates the newunit fix.")
    return np.vstack(rows)

def main():
    args = sys.argv[1:]
    region, floor = None, None
    if '--floor' in args:
        i = args.index('--floor')
        floor = float(args[i+1])
        args = args[:i] + args[i+2:]
    if '--region' in args:
        i = args.index('--region')
        region = [float(x) for x in args[i+1:i+5]]
        args = args[:i] + args[i+5:]
    paths = args or sorted(glob.glob('bc_geom_diag_*.dat'))
    if not paths:
        sys.exit("no bc_geom_diag_*.dat files given or found in the current directory")

    d = load(paths)
    R, Z, bt, ivd, bn, dirn, fac = (d[:,0], d[:,1], d[:,2].astype(int),
                                    d[:,3].astype(int), d[:,4], d[:,5], d[:,6])
    psb = d[:,7] if d.shape[1] > 7 else None
    print(f"read {len(paths)} file(s), {len(d)} visits")

    # --- group visits by node position (coordinates are written with fixed precision)
    key = np.round(np.stack([R, Z], axis=1), 8)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    print(f"{len(uniq)} distinct boundary nodes, "
          f"{np.bincount(inv).max()} visits to the most-visited node")

    conflict, fspread = [], []
    for i in range(len(uniq)):
        m = inv == i
        ds = np.unique(np.sign(dirn[m]))
        if len(ds) > 1:
            conflict.append((uniq[i,0], uniq[i,1], bt[m][0],
                             np.abs(bn[m]).min(), fac[m].min(), fac[m].max()))
        elif fac[m].max() - fac[m].min() > 1e-3:
            fspread.append((uniq[i,0], uniq[i,1], bt[m][0],
                            np.abs(bn[m]).min(), fac[m].min(), fac[m].max()))

    def show(title, rows, note):
        print(f"\n=== {title}: {len(rows)} ===")
        if not rows:
            print("    none")
            return
        print(f"    {note}")
        print(f"    {'R':>12} {'Z':>12} {'type':>5} {'min|bn|':>12} {'factor lo':>12} {'factor hi':>12}")
        for r in sorted(rows, key=lambda r: r[3]):
            print(f"    {r[0]:12.6f} {r[1]:12.6f} {r[2]:5d} {r[3]:12.3e} {r[4]:12.4f} {r[5]:12.4f}")

    show("NODES WHERE direction DISAGREES BETWEEN VISITS", conflict,
         "v_par sign here is decided by element ordering, not physics.")
    show("NODES WHERE factor DISAGREES BETWEEN VISITS", fspread,
         "smoothing fires on one side of the corner only.")

    # --- direction flips between neighbouring nodes, independent of the per-node check
    un_dir = np.array([np.sign(dirn[inv == i]).mean() for i in range(len(uniq))])
    un_bn  = np.array([bn[inv == i].mean() for i in range(len(uniq))])
    # Check the two nearest neighbours, not one: along a boundary a node has a neighbour on
    # either side, and with even spacing argmin picks whichever comes first in the array. Dedupe
    # on the node pair rather than on i < j - the nearest-neighbour relation is not symmetric, so
    # an i < j test silently drops flips that only the higher-indexed node can see.
    seen, flips = set(), []
    for i in range(len(uniq)):
        dist = np.hypot(uniq[:,0]-uniq[i,0], uniq[:,1]-uniq[i,1])
        dist[i] = np.inf
        for j in np.argsort(dist)[:2]:
            j = int(j)
            if not np.isfinite(dist[j]):
                continue
            pair = (min(i, j), max(i, j))
            if pair in seen:
                continue
            if un_dir[i]*un_dir[j] < 0:
                seen.add(pair)
                flips.append((uniq[i,0], uniq[i,1], dist[j],
                              abs(un_bn[i]), abs(un_bn[j])))
    print(f"\n=== direction FLIPS BETWEEN NEAREST-NEIGHBOUR NODES: {len(flips)} ===")
    if flips:
        print("    A flip is expected where the field grazes the wall; it is dangerous when |bn|")
        print("    is NOT small there, i.e. a genuine sign change imposed at finite incidence.")
        print(f"    {'R':>12} {'Z':>12} {'spacing':>12} {'|bn| a':>12} {'|bn| b':>12}")
        for f in sorted(flips, key=lambda f: -max(f[3], f[4])):
            print(f"    {f[0]:12.6f} {f[1]:12.6f} {f[2]:12.3e} {f[3]:12.3e} {f[4]:12.3e}")
    else:
        print("    none")

    # --- node spacing along the boundary. A constant spacing comes from a uniform
    # --- float(j-1)/float(n-1) distribution; a smoothly varying one comes from meshac2 with
    # --- SIG_*. Where a constant run butts against a varying one, two different constructions
    # --- have been joined, and the cell size steps abruptly at the seam.
    print("\n=== ABRUPT SPACING CHANGES ALONG THE BOUNDARY (ratio > 2) ===")
    order = np.argsort(np.arctan2(uniq[:,1] - uniq[:,1].mean(),
                                  uniq[:,0] - uniq[:,0].mean()))
    ou = uniq[order]
    gap = np.hypot(np.diff(ou[:,0]), np.diff(ou[:,1]))
    steps = []
    for i in range(1, len(gap)):
        if gap[i] <= 0 or gap[i-1] <= 0:
            continue
        r = max(gap[i]/gap[i-1], gap[i-1]/gap[i])
        # skip the seam of the angular sort itself and any wrap-around artefact
        if r > 2.0 and max(gap[i], gap[i-1]) < 0.5:
            steps.append((ou[i,0], ou[i,1], gap[i-1], gap[i], r))
    if steps:
        print(f"    {'R':>12} {'Z':>12} {'before':>12} {'after':>12} {'ratio':>8}")
        for st in sorted(steps, key=lambda s: -s[4])[:25]:
            print(f"    {st[0]:12.6f} {st[1]:12.6f} {st[2]:12.6f} {st[3]:12.6f} {st[4]:8.1f}")
        print(f"    ({len(steps)} in total; ordering is by polar angle about the mesh centre,")
        print("     so a few entries near the sort seam may be spurious - check R,Z)")
    else:
        print("    none")

    # --- ps0_b, the quantity the Mach-1 ExB term divides by. mach1_psib_floor should sit below
    # --- the values on the targets (so the physical term is untouched there) and above the
    # --- near-tangential tail (so the amplification is bounded).
    if psb is not None:
        a = np.abs(psb[psb != 0.0])
        if a.size:
            print("\n=== |ps0_b| DISTRIBUTION (for choosing mach1_psib_floor) ===")
            for q in (0.1, 1, 5, 25, 50):
                print(f"    {q:5.1f}th percentile: {np.percentile(a, q):.4e}")
            print(f"    minimum          : {a.min():.4e}")
            i = int(np.argmin(np.abs(psb)))
            print(f"    smallest |ps0_b| at R={R[i]:.4f} Z={Z[i]:.4f}, type {bt[i]}, bn={bn[i]:.3e}")
            print("    A floor near the 1st-5th percentile bounds the tail while leaving the")
            print("    bulk of the boundary, including the targets, unchanged.")

    # --- impact of a candidate mach1_psib_floor, so the choice can be checked before running
    if psb is not None and floor is not None:
        a = np.abs(psb)
        raw  = np.where(a > 0, 1.0/np.where(a > 0, a, 1.0), 0.0)
        regd = a / (a**2 + floor**2)
        with np.errstate(divide='ignore', invalid='ignore'):
            frac = np.where(raw > 0, regd/raw, 1.0)     # 1 = untouched, <1 = damped
        print(f"\n=== IMPACT OF mach1_psib_floor = {floor:.3e} ===")
        print(f"    bound on 1/ps0_b        : {1.0/(2*floor):.3e}   (was up to {raw.max():.3e})")
        for tag, lo, hi in (("untouched  (>99%)", 0.99, 1.01),
                            ("mild  (90-99%)   ", 0.90, 0.99),
                            ("damped (50-90%)  ", 0.50, 0.90),
                            ("strong (<50%)    ", -1.0, 0.50)):
            m = (frac > lo) & (frac <= hi)
            print(f"    {tag}: {m.sum():5d} visits ({100.0*m.sum()/len(frac):5.1f}%)")
        for t in np.unique(bt):
            m = bt == t
            print(f"    type {t}: median retention {np.median(frac[m]):.4f}, "
                  f"worst {frac[m].min():.4f}")
        print("    Retention near 1 on the target types means the physical term is unchanged")
        print("    there; the damping should fall on the near-tangential wall only.")

    if region:
        Rmin, Rmax, Zmin, Zmax = region
        m = (R >= Rmin) & (R <= Rmax) & (Z >= Zmin) & (Z <= Zmax)
        print(f"\n=== ALL VISITS IN R=[{Rmin},{Rmax}] Z=[{Zmin},{Zmax}]: {m.sum()} ===")
        print(f"    {'R':>12} {'Z':>12} {'type':>5} {'iv_dir':>7} {'bn':>13} {'dir':>6} {'factor':>12}")
        for i in np.argsort(-Z[m]):
            print(f"    {R[m][i]:12.6f} {Z[m][i]:12.6f} {bt[m][i]:5d} {ivd[m][i]:7d} "
                  f"{bn[m][i]:13.4e} {dirn[m][i]:6.1f} {fac[m][i]:12.4f}")

        # distinct nodes in the box, ordered along the boundary, with the spacing to the next
        ru = np.unique(np.round(np.stack([R[m], Z[m]], axis=1), 8), axis=0)
        ru = ru[np.argsort(-ru[:,1])]
        print(f"\n    --- {len(ru)} distinct nodes, spacing to the next one ---")
        print(f"    {'R':>12} {'Z':>12} {'spacing':>12}   note")
        prev = None
        for i in range(len(ru)-1):
            sp = float(np.hypot(ru[i+1,0]-ru[i,0], ru[i+1,1]-ru[i,1]))
            note = ""
            if prev and prev > 0:
                r = max(sp/prev, prev/sp)
                if r > 2.0:
                    note = f"<-- spacing steps {r:.1f}x"
            print(f"    {ru[i,0]:12.6f} {ru[i,1]:12.6f} {sp:12.6f}   {note}")
            prev = sp

if __name__ == '__main__':
    main()
