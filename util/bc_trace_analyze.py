#!/usr/bin/env python3
"""Analyse bc_float_trace_*.dat: is the imposed floating trace smooth or noisy?

Columns: R Z bnd_type iv_dir target_u u_val u_dof Te_val Te_dof target_dof

The derivative row slaves du/dl to dTe/dl through the RAW nodal temperature derivative DOF.
Te at the wall carries a natural BC, so that DOF is unconstrained, and grid-scale noise in it
goes straight into the boundary electric field and into w. This measures whether the noise is
actually there, rather than assuming it.

The diagnostic quantity is the ROUGHNESS: how much of a quantity's variation sits at the
grid scale. A smooth physical profile has roughness << 1; a signal alternating node to node
has roughness ~ 1.

Usage: ./bc_trace_analyze.py [files ...] [--region Rmin Rmax Zmin Zmax]
"""
import sys, glob, warnings
import numpy as np

def roughness(y):
    """Fraction of variance at the Nyquist (node-to-node alternating) scale.

    Second difference picks out grid-scale content; normalising by the total variation
    makes it independent of amplitude, so smooth-but-large and small-but-noisy are
    distinguished rather than confused.
    """
    if len(y) < 5:
        return float('nan')
    d2 = y[2:] - 2*y[1:-1] + y[:-2]
    d1 = np.abs(np.diff(y))
    if d1.sum() == 0:
        return 0.0
    return float(np.sqrt(np.mean(d2**2)) / (2.0*np.mean(d1) + 1e-300))

def main():
    args = sys.argv[1:]
    region = None
    if '--region' in args:
        i = args.index('--region')
        region = [float(x) for x in args[i+1:i+5]]
        args = args[:i] + args[i+5:]
    paths = args or sorted(glob.glob('bc_float_trace_*.dat'))
    if not paths:
        sys.exit("no bc_float_trace_*.dat found")
    rows = []
    for p in paths:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            d = np.loadtxt(p, comments='#', ndmin=2)
        if d.size:
            rows.append(d)
    if not rows:
        sys.exit("all trace files are empty")
    d = np.vstack(rows)
    print(f"read {len(paths)} file(s), {len(d)} visits")

    R, Z, bt = d[:,0], d[:,1], d[:,2].astype(int)
    key = np.round(np.stack([R, Z], axis=1), 8)
    uniq, idx = np.unique(key, axis=0, return_index=True)
    d = d[idx]
    R, Z, bt = d[:,0], d[:,1], d[:,2].astype(int)
    print(f"{len(d)} distinct boundary nodes")

    if region:
        m = (R>=region[0]) & (R<=region[1]) & (Z>=region[2]) & (Z<=region[3])
        d, R, Z, bt = d[m], R[m], Z[m], bt[m]
        print(f"restricted to region: {len(d)} nodes")

    # order along the boundary by walking nearest-neighbour from an endpoint
    order = np.argsort(-Z)
    d = d[order]
    names = ['target_u','u_val','u_dof','Te_val','Te_dof','target_dof']
    print(f"\n=== ROUGHNESS (fraction of variation at the node-to-node scale) ===")
    print("    < 0.2  smooth      0.2-0.5  some grid-scale content      > 0.5  noise-dominated")
    print(f"    {'quantity':>12} {'roughness':>11} {'min':>13} {'max':>13}")
    for j, nm in enumerate(names):
        y = d[:, 4+j]
        r = roughness(y)
        flag = '' if r < 0.2 else ('  <-- grid-scale content' if r < 0.5 else '  <-- NOISE DOMINATED')
        print(f"    {nm:>12} {r:11.3f} {y.min():13.4e} {y.max():13.4e}{flag}")
    print("\n  target_dof is the slope actually imposed on u by the derivative row.")
    print("  If its roughness is much higher than target_u's, the derivative row is")
    print("  injecting grid-scale structure that the value row is not - which is the case")
    print("  for building slopes from a reconstructed target instead of raw Te DOFs.")

if __name__ == '__main__':
    main()
