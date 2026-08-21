#!/usr/bin/env python3
"""Round selected corners of a JOREK wall polygon.

Format: first line is the point count, then one "R Z" pair per line, the polygon
closed by repeating the first point last.

Why round rather than refine: find_wall_crossing projects grid rays onto the
polygon, so inserting COLLINEAR points along a straight segment changes nothing
at all - a ray crosses the same line in the same place. Only the corners change
where crossings land. A hard corner joining two long segments makes the crossing
spacing vary sharply as the ray sweeps past it, which shows up as a cluster of
boundary nodes far tighter than their neighbours.

Each selected vertex is replaced by a quadratic Bezier from a point `radius`
back along the incoming edge to a point `radius` along the outgoing edge, using
the original vertex as the control point. That is tangent-continuous with both
edges and stays inside the original corner, so the wall never moves outward.

Usage:
    ./round_wall_corners.py in.dat out.dat --corner R Z [--corner R Z ...]
                            [--radius 0.03] [--points 7] [--drop-collinear 0.5]
"""
import sys
import numpy as np

def load(path):
    with open(path) as f:
        tok = f.read().split()
    n = int(tok[0])
    pts = np.array([float(x) for x in tok[1:]]).reshape(-1, 2)
    if len(pts) != n:
        print(f"note: header says {n} points, file holds {len(pts)}", file=sys.stderr)
    return pts

def save(path, pts):
    with open(path, 'w') as f:
        f.write(f"{len(pts)}\n")
        for r, z in pts:
            f.write(f"{r:.6f}  {z:.6f}\n")

def round_corner(pts, idx, radius, npt):
    closed = np.allclose(pts[0], pts[-1])
    n = len(pts) - 1 if closed else len(pts)
    prev, here, nxt = pts[(idx-1) % n], pts[idx], pts[(idx+1) % n]

    din, dout = here - prev, nxt - here
    lin, lout = np.linalg.norm(din), np.linalg.norm(dout)
    r = min(radius, 0.45*lin, 0.45*lout)     # never eat more than half an edge
    a = here - din/lin*r
    b = here + dout/lout*r

    t = np.linspace(0.0, 1.0, npt)[:, None]
    arc = (1-t)**2 * a + 2*(1-t)*t*here + t**2 * b
    turn = np.degrees(np.arccos(np.clip(np.dot(din/lin, dout/lout), -1, 1)))
    print(f"  vertex {idx} at ({here[0]:.6f}, {here[1]:.6f}): turns {turn:.1f} deg, "
          f"edges {lin:.3f}/{lout:.3f} m, rounded with r = {r:.4f} m, {npt} points")
    return np.vstack([pts[:idx], arc, pts[idx+1:]])

def drop_collinear(pts, tol_deg):
    """Remove vertices that turn by less than tol_deg - they do not change the shape.

    Worth doing before rounding: a redundant vertex sitting a couple of mm from a real corner
    caps the rounding radius at a fraction of that 2 mm stub, so the corner cannot be rounded at
    all while it is there. Removing it is exactly shape-preserving.
    """
    closed = np.allclose(pts[0], pts[-1])
    core = pts[:-1] if closed else pts
    keep, dropped = [], []
    n = len(core)
    for i in range(n):
        a = core[i] - core[(i-1) % n]
        b = core[(i+1) % n] - core[i]
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        if na == 0 or nb == 0:
            dropped.append(core[i]); continue
        ang = np.degrees(np.arccos(np.clip(a@b/na/nb, -1, 1)))
        (keep if ang > tol_deg else dropped).append(core[i])
    for d in dropped:
        print(f"  dropped collinear vertex ({d[0]:.6f}, {d[1]:.6f})")
    out = np.array(keep)
    if closed:
        out = np.vstack([out, out[0]])
    print(f"  {len(dropped)} vertices dropped, {len(out)} remain")
    return out

def main():
    args = sys.argv[1:]
    if len(args) < 2:
        sys.exit(__doc__)
    src, dst = args[0], args[1]
    radius, npt, corners, collin = 0.03, 7, [], None
    i = 2
    while i < len(args):
        if args[i] == '--corner':
            corners.append((float(args[i+1]), float(args[i+2]))); i += 3
        elif args[i] == '--radius':
            radius = float(args[i+1]); i += 2
        elif args[i] == '--points':
            npt = int(args[i+1]); i += 2
        elif args[i] == '--drop-collinear':
            collin = float(args[i+1]); i += 2
        else:
            sys.exit(f"unrecognised argument: {args[i]}")
    if not corners:
        sys.exit("give at least one --corner R Z")

    pts = load(src)
    print(f"read {len(pts)} points from {src}")
    if collin is not None:
        print(f"removing vertices that turn by less than {collin} deg:")
        pts = drop_collinear(pts, collin)
    # round from the highest index down, so earlier indices stay valid
    idxs = []
    for R, Z in corners:
        d = np.hypot(pts[:, 0]-R, pts[:, 1]-Z)
        j = int(np.argmin(d))
        if d[j] > 1e-4:
            sys.exit(f"no vertex within 1e-4 of ({R}, {Z}); nearest is "
                     f"({pts[j,0]:.6f}, {pts[j,1]:.6f}) at {d[j]:.2e}")
        idxs.append(j)
    for j in sorted(idxs, reverse=True):
        pts = round_corner(pts, j, radius, npt)

    save(dst, pts)
    print(f"wrote {len(pts)} points to {dst}")

if __name__ == '__main__':
    main()
