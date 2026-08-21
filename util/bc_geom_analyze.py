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
    ./bc_geom_analyze.py [files ...] [--region Rmin Rmax Zmin Zmax]
"""
import sys, glob
import numpy as np

def load(paths):
    rows = []
    for p in paths:
        d = np.loadtxt(p, comments='#', ndmin=2)
        if d.size:
            rows.append(d)
    if not rows:
        sys.exit("no data found - check that bc_geom_diag_*.dat exist and are non-empty")
    return np.vstack(rows)

def main():
    args = sys.argv[1:]
    region = None
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

    if region:
        Rmin, Rmax, Zmin, Zmax = region
        m = (R >= Rmin) & (R <= Rmax) & (Z >= Zmin) & (Z <= Zmax)
        print(f"\n=== ALL VISITS IN R=[{Rmin},{Rmax}] Z=[{Zmin},{Zmax}]: {m.sum()} ===")
        print(f"    {'R':>12} {'Z':>12} {'type':>5} {'iv_dir':>7} {'bn':>13} {'dir':>6} {'factor':>12}")
        for i in np.argsort(-Z[m]):
            print(f"    {R[m][i]:12.6f} {Z[m][i]:12.6f} {bt[m][i]:5d} {ivd[m][i]:7d} "
                  f"{bn[m][i]:13.4e} {dirn[m][i]:6.1f} {fac[m][i]:12.4f}")

if __name__ == '__main__':
    main()
