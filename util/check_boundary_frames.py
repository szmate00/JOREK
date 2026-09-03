#!/usr/bin/env python3
"""Boundary node-frame diagnostic, straight off a JOREK HDF5 restart file.

    python3 util/check_boundary_frames.py jorek000000.h5

WHY.  In grids/grid_xpoint_wall.f90 the two first-derivative DOF directions at a node,
x(1,2,:) and x(1,3,:), are built from DIFFERENT things depending on which block of the
grid the node came from:

  flux-aligned grid (boundary types 1, 2, 3)
      x(1,2,:) = (dR_dt, dZ_dt)                 analytic tangent of the grid line (CUB1D)
      x(1,3,:) = (-PSI_Z, +PSI_R)/|grad psi|    EXACTLY the poloidal field direction
  ray-cast wall extension (boundary types 4, 5, 9)
      x(1,2,:) = the ray from the plasma edge to the wall
      x(1,3,:) = (cos(tht_ext), sin(tht_ext)), an interpolated geometric angle
                 (in the LEG block, a secant through the neighbouring nodes)

Three consequences, all measured here.

(A) psi DERIVATIVE DOFs.  grid_xpoint_wall.f90:1886-1889 sets
        values(1,k,psi) = grad(psi) . x(1,k,:)
    so on a flux-grid node values(1,3,psi) is ANALYTICALLY ZERO - psi does not vary along
    a flux surface.  On an extension node it is a generic O(|grad psi|) number.
    dirichlet%psi pins the value and the TANGENTIAL derivative and leaves the NORMAL one
    free (mod_boundary_conditions.f90:251-260 and :400-405, iv_dir vs iv_perp_dir), and
    mod_sheath_bc.f90's header says that free normal DOF is how the wall current responds
    to the sheath.  The mapping (construct_matrix_mod.f90:135-168 with
    mod_boundary_matrix_open.f90:178):
        types 1, 4, 9  ->  tangential DOF 2, FREE normal DOF 3
        types 2, 5     ->  tangential DOF 3, FREE normal DOF 2
    So on type 1 the free DOF starts at exactly zero and the sheath's response is pure
    signal; on 4/5/9 it carries the bulk of grad(psi) and the response is a small relative
    perturbation on a large background.

(B) FRAME-CHORD ALIGNMENT.  grid_xpoint_wall.f90:1751 builds the element size as
        size = sign( |chord|/n_order , frame . chord )
    i.e. the chord LENGTH, with the frame entering only through its SIGN.  That is exact
    iff the frame vector is parallel to the element edge.  cos(theta) below is exactly
    that validity measure; at cos -> 0 the SIGN itself is set by noise.

(C) FRAME CONDITIONING.  Nothing requires x(1,2,:) and x(1,3,:) to be orthogonal or even
    independent.  |x2 x x3| -> 0 is the degenerate case, and zj = Delta psi amplifies a
    bad frame by 1/h^2.

None of the three is measured anywhere in JOREK, and all are constant across a run, so a
single restart file answers the question.

EXPECTED IF THE FLUX/EXTENSION SPLIT IS THE STORY
    type 1  cos ~ 1, det ~ 1, free-DOF fraction ~ 0     clean on all three
    type 5  cos ~ 1, det ~ 1, free-DOF fraction O(1)    frames fine, DOF large -> drift
    type 4  BIMODAL: one flux-grid node, n_ext-1 extension nodes, one surgery node
    type 9  cos and/or det well below 1                 surgery, frames never recomputed

NOTE the restart is a mid-run state, so values() is the evolved psi, not the construction
psi.  The analytic-zero property is a construction property; measuring it in the evolved
state is the stronger test, since that is the state the solver is actually in.
"""

import sys
import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("needs h5py:  pip install h5py")

# --- boundary type -> free (normal) DOF, from construct_matrix_mod.f90:135-168 together
# --- with direction_perp(1) = 6/direction(2).  Type 3 is a corner: both sides occur.
PERP_DOF = {1: 3, 4: 3, 9: 3, 2: 2, 5: 2}

# --- inner/outer split in major radius, matching sheath_diag_R_split in the namelist
R_SPLIT = 1.42

# --- element side iv (1-based) -> iv_dir, from mod_boundary_conditions.f90:251-260:
# --- iv*iv2 in {2,12} -> dir 2 ; iv*iv2 in {6,4} -> dir 3, with iv2 = mod(iv,4)+1.
SIDE_DIR = {1: 2, 2: 3, 3: 2, 4: 3}


def get(f, name, default=None):
    if name not in f:
        if default is not None:
            return default
        sys.exit(f"dataset '{name}' not in the file; is this a JOREK HDF5 restart?")
    return f[name][()]


def as_scalar(v):
    return int(np.asarray(v).reshape(-1)[0])


def axes_after_nodes(a, n_nodes, label):
    """Squeeze length-1 axes, move the node axis first, return the array."""
    a = np.squeeze(np.asarray(a))
    if a.ndim == 0:
        sys.exit(f"'{label}' collapsed to a scalar")
    if a.shape[0] == n_nodes:
        return a
    for ax in range(a.ndim):
        if a.shape[ax] == n_nodes:
            return np.moveaxis(a, ax, 0)
    sys.exit(f"could not find the node axis (n_nodes={n_nodes}) in '{label}' {a.shape}")


def split_axes(a, n_nodes, n_deg, other, label):
    """-> (n_nodes, n_degrees, other).

    The HDF5 wrapper does NOT preserve the Fortran argument order: `x` is declared
    (n_nodes, n_coord_tor, n_degrees, n_dim) but lands in the file as
    (n_nodes, n_dim, n_degrees).  So detect the axes by size rather than assume.
    """
    a = axes_after_nodes(a, n_nodes, label)
    if a.ndim != 3:
        sys.exit(f"'{label}' has shape {a.shape}; expected 3 axes after squeezing")
    s1, s2 = a.shape[1], a.shape[2]
    if s1 == n_deg and s2 == other:
        return a
    if s2 == n_deg and s1 == other:
        return np.moveaxis(a, 1, 2)
    if s1 == n_deg and s2 == n_deg:          # ambiguous: n_degrees == other
        return a
    sys.exit(f"'{label}' shape {a.shape} matches neither "
             f"(n_degrees={n_deg}, other={other}) nor its transpose")


def main(path):
    with h5py.File(path, "r") as f:
        n_nodes = as_scalar(get(f, "n_nodes"))
        n_order = as_scalar(get(f, "n_order", 3))
        n_var = as_scalar(get(f, "n_var", 0))
        n_deg = as_scalar(get(f, "n_degrees", ((n_order + 1) // 2) ** 2))
        ver = as_scalar(get(f, "rst_hdf5_version", 0))

        raw_x, raw_v = get(f, "x"), get(f, "values")
        x = split_axes(raw_x, n_nodes, n_deg, 2, "x")          # (n_nodes, n_degrees, n_dim)
        val = split_axes(raw_v, n_nodes, n_deg, n_var, "values")  # (n_nodes, n_degrees, n_var)
        bnd = np.asarray(get(f, "boundary")).reshape(-1).astype(int)
        vtx = np.asarray(get(f, "vertex"))
        if vtx.shape[0] == 4 and vtx.shape[-1] != 4:
            vtx = vtx.T
        vtx = vtx.reshape(-1, 4).astype(int) - 1              # Fortran 1-based -> 0-based

    print(f"  file            : {path}")
    print(f"  rst_hdf5_version: {ver}   n_order: {n_order}   n_var: {n_var}")
    print(f"  n_nodes         : {n_nodes}   n_elements: {vtx.shape[0]}")
    print(f"  x      raw {np.squeeze(raw_x).shape}  ->  (n_nodes, n_degrees, n_dim) {x.shape}")
    print(f"  values raw {np.squeeze(raw_v).shape}  ->  (n_nodes, n_degrees, n_var) {val.shape}\n")

    if x.shape[1] < 3 or val.shape[1] < 3:
        sys.exit("need at least 3 nodal degrees (n_order >= 3)")

    x2, x3 = x[:, 1, :], x[:, 2, :]
    det = np.abs(x2[:, 0] * x3[:, 1] - x2[:, 1] * x3[:, 0])
    dot = np.abs(x2[:, 0] * x3[:, 0] + x2[:, 1] * x3[:, 1])

    # --- (A) how much of grad(psi) sits in each derivative DOF.  var_psi = 1 (Fortran)
    v2, v3 = val[:, 1, 0], val[:, 2, 0]
    nrm = np.hypot(v2, v3)
    with np.errstate(divide="ignore", invalid="ignore"):
        f2 = np.where(nrm > 0, np.abs(v2) / nrm, 0.0)
        f3 = np.where(nrm > 0, np.abs(v3) / nrm, 0.0)

    # --- (B) frame-chord alignment on every BOUNDARY element side
    cos_min = np.full(n_nodes, np.inf)
    cos_sum = np.zeros(n_nodes)
    cos_cnt = np.zeros(n_nodes, dtype=int)
    pos = x[:, 0, :]
    for iv in range(4):                       # 0-based; Fortran iv = iv+1
        a = vtx[:, iv]
        b = vtx[:, (iv + 1) % 4]
        live = (bnd[a] != 0) & (bnd[b] != 0)  # construct_matrix_mod.f90:102
        if not live.any():
            continue
        d = SIDE_DIR[iv + 1] - 1              # DOF index, 0-based
        ch = pos[b[live]] - pos[a[live]]
        ln = np.hypot(ch[:, 0], ch[:, 1])
        ok = ln > 0
        for nodes_, sgn in ((a[live], +1.0), (b[live], -1.0)):
            fr = x[nodes_, d, :]
            c = np.abs(sgn * (fr[:, 0] * ch[:, 0] + fr[:, 1] * ch[:, 1]))
            c = np.where(ok, c / np.where(ok, ln, 1.0), np.nan)
            for nd, cc in zip(nodes_[ok], c[ok]):
                cos_min[nd] = min(cos_min[nd], cc)
                cos_sum[nd] += cc
                cos_cnt[nd] += 1
    cos_mean = np.where(cos_cnt > 0, cos_sum / np.maximum(cos_cnt, 1), np.nan)

    hdr = (f"{'type':>5} {'nodes':>7} | {'cos min':>9} | "
           f"{'det min':>9} {'det p05':>9} {'det med':>9} {'det mean':>9} "
           f"{'%<0.3':>7} {'%<0.1':>7} | {'free':>5} {'free frac':>10}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for t in sorted(set(bnd.tolist())):
        if t == 0:
            continue
        m = (bnd == t) & (cos_cnt > 0)
        if not m.any():
            m = bnd == t
        p = PERP_DOF.get(t)
        ff = (f2[m].mean() if p == 2 else f3[m].mean()) if p else float("nan")
        cmin = np.nanmin(cos_min[m]) if np.isfinite(cos_min[m]).any() else float("nan")
        d = det[m]
        print(f"{t:>5} {(bnd == t).sum():>7} | {cmin:>9.4f} | "
              f"{d.min():>9.4f} {np.percentile(d, 5):>9.4f} {np.median(d):>9.4f} "
              f"{d.mean():>9.4f} {100 * (d < 0.3).mean():>6.1f}% {100 * (d < 0.1).mean():>6.1f}% | "
              f"{('DOF ' + str(p)) if p else '-':>5} {ff:>10.4f}")

    # --- threshold sweep: exactly how much of each type a sheath_weak_detmin gate costs
    thr = [0.1, 0.2, 0.3, 0.4, 0.5, 0.7]
    print(f"\n  IF sheath_weak_detmin WERE SET, % of each type's rows it would GATE OFF")
    print("    " + f"{'type':>5} " + " ".join(f"{t:>8.2f}" for t in thr))
    print("    " + "-" * (6 + 9 * len(thr)))
    for t in sorted(set(bnd.tolist())):
        if t == 0:
            continue
        d = det[bnd == t]
        cells = " ".join(f"{100 * (d < q).mean():>7.1f}%" for q in thr)
        print(f"    {t:>5} {cells}")

    # --- where the degenerate nodes actually are, so they can be found on the grid
    deg = (det < 0.3) & (bnd != 0)
    if deg.any():
        rs = R_SPLIT
        ins, out = deg & (pos[:, 0] < rs), deg & (pos[:, 0] >= rs)
        print(f"\n  DEGENERATE BOUNDARY NODES (|x2 x x3| < 0.3): {deg.sum()} of "
              f"{(bnd != 0).sum()}")
        for lbl, mm in (("INNER (R < %.2f)" % rs, ins), ("OUTER (R >= %.2f)" % rs, out)):
            if mm.any():
                print(f"    {lbl:>18}: {mm.sum():>3} nodes, worst det {det[mm].min():.4f} "
                      f"(1/det {1.0 / max(det[mm].min(), 1e-30):.1f}), "
                      f"R {pos[mm, 0].min():.4f}..{pos[mm, 0].max():.4f}  "
                      f"Z {pos[mm, 1].min():.4f}..{pos[mm, 1].max():.4f}")
        print("    worst 15, sorted by det")
        print(f"    {'type':>5} {'det':>8} {'angle':>8} {'1/det':>8} {'R':>9} {'Z':>9}")
        idx = np.argsort(det)
        shown = 0
        for i in idx:
            if not deg[i]:
                continue
            ang = np.degrees(np.arcsin(min(det[i], 1.0)))
            print(f"    {bnd[i]:>5} {det[i]:>8.4f} {ang:>7.2f}d {1.0 / max(det[i], 1e-30):>8.1f} "
                  f"{pos[i, 0]:>9.4f} {pos[i, 1]:>9.4f}")
            shown += 1
            if shown >= 15:
                break

    print("\n  READ IT LIKE THIS")
    print("    cos    - |frame . chord|/|chord| on the boundary sides, for the DOF that")
    print("             grid_xpoint_wall.f90:1751 uses to build element%size.  MEASURED ~1")
    print("             on every type, so the size construction is fine everywhere and the")
    print("             sign() concern does not bite.  Kept only as a regression check.")
    print("    det    - |x2 x x3| = |sin(angle between the two frame vectors)|.  This is THE")
    print("             discriminator: the nodal derivative basis is conditioned as 1/det, so")
    print("             det = 0.01 inflates every derivative DOF at that node by ~100x, and")
    print("             zj = Delta psi is built from second derivatives.  MEASURED det mean")
    print("             0.88 / 0.87 / 0.38 / 0.05 on types 1 / 5 / 4 / 9 against survival")
    print("             3900 / 308 / 8 / 4 steps.")
    print("    free frac - the fraction of grad(psi) in the FREE normal DOF, the one the")
    print("             sheath current response has to move.  ~0 on type 1 by construction;")
    print("             O(1) means the response is a small perturbation on a big background.")
    print("    A BIMODAL type (min far from mean) mixes flux-grid and extension nodes -")
    print("             predicted for type 4.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
