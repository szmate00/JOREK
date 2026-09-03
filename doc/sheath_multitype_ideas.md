# Getting the weak sheath j-V BC to run on multiple boundary types

Written 2026-09-03 against `sheath-jsat-vpar-38ab278`. Everything below is either read
off the source or arithmetic from numbers already in the logs; nothing here has been run.

## The fact any explanation has to fit

    type 1 alone          ~3900 steps, TIMEOUT, converged, weak residual still FALLING
    type 5 alone           308  crash
    type 1 + 5             305  crash   (weak 3.3e-3, I closure to 4 digits)
    type 1 + 4             220  crash   (weak flat 0.20, closure ~1%)
    type 4 alone             8  crash
    type 9 alone             4  crash

Type 1 alone stops being a stability problem. Every multi-type case dies **while still
evolving**. That is a difference of kind. All the gate/geometry/weight hypotheses were
about degree, which is why none of them discriminated.

And the 1+5 failure is not a spike. It is a *monotone ramp* on the inner target only:

    inner ePhi/kTe   3.96  4.18  4.37  4.73  5.03  5.50  6.89  8.68  10.02   -> blow-up
    outer ePhi/kTe   1.05 ................................................ 1.06

A one-sided monotone drift over 300 steps is an integrator, not an instability.

---

## Idea 1 (top pick): the row loses its authority over `u` exactly where it is needed

`mod_sheath_bc.f90:560` — the only thing coupling `u` into the replaced row is

    dzj_du = zj_sat * fp * dx_du ,    fp = exp(-X)*dxlim_dx + s*sigma(X)

With `sheath_sat_slope = 0.03` and `Lambda = 3`, using `X = ePhi/kTe - Lambda`:

| where | ePhi/kTe | X | fp | u-authority |
|---|---|---|---|---|
| outer target       | 1.05  | -1.95 | 7.03  | 1x |
| inner, step ~10    | 3.96  | +0.96 | 0.405 | 17x weaker |
| inner, step ~300   | 10.02 | +7.02 | 0.031 | **227x weaker** |
| inner, 1+4 at 219  | 15.3  | +12.3 | 0.030 | 234x weaker |

So there is a **positive feedback loop**: more ion saturation -> smaller `fp` -> weaker
`u` column -> `u` free to drift further positive -> more ion saturation. The loop gain is
set by `fp`, which falls exponentially. This is a self-reinforcing drift with a slow
timescale, which is the observed signature exactly.

**Why only multi-type.** With type 1 alone, types 2/3/4/5/9 keep `dirichlet%u = .true.`.
Type 5 is the big wall — 724 trace rows against type 1's 196. That Dirichlet holds `u`
over most of the SOL, and the elliptic problem propagates it into the leg, so `u` at the
strike point is pinned by its *neighbours* even when `fp -> 0`. Release type 5 and that
anchor is gone: `u` over the whole SOL is now held only by rows whose coupling is `0.03`.
Nothing else changes by a factor of 200 when a second type is enabled.

This also explains why `sheath_weak_wmin`, `sheath_min_bn`, ufade, geometry, `b_n` sign
and cell size all failed as discriminators — none of them touch `fp`.

### The cheap test (do this first)

`sheath_dfdx_min` is exactly the regulariser for this and, per
[[sheath-multitype-working-config]], it was only ever tested against the **poisoned
`sheath_init_u = .true.` baseline** (8 -> 9 steps). It has never been run on the 305-step
1+5 or the 220-step 1+4 baseline.

    sheath_dfdx_min = 0.3

At the floating potential `fp = 1`, so 0.3 is inactive in the regime of interest; on the
outer target `fp = 7`, inactive; on the drifting inner target it restores a factor of ~10
of the u-authority precisely where the loop is running. `mod_sheath_bc.f90:539` applies it
to `fp` only and never to `f`, so the fixed point does not move.

**Prediction if idea 1 is right:** the inner `ePhi/kTe` ramp stalls or slows markedly,
and the run passes 305 without the `max|j/jsat|` excursion. **If the ramp is unchanged,
idea 1 is dead** — that is a clean falsification, unlike most of this campaign.

### The proper fix, if the test is positive

`sheath_dfdx_min` is inconsistent: it changes `J` without changing `F`, so it damps but
also degrades the Newton rate. The clean version is to **write the row in potential space
instead of current space**. `F = 0` and `F/fp = 0` have the identical solution set, so
dividing residual and every Jacobian column by `max(fp, fp_min)` *at each Gauss point*
(it varies within a row, so it does not cancel) is a legal re-weighting of the Galerkin
test function that leaves the fixed point exactly where it was, and makes the `u` column
O(`zj_sat * dx_du`) — independent of how deep in saturation the surface sits.

Concretely: return `fp` from `sheath_current` (one extra optional argument), and in
`mod_boundary_matrix_open.f90` around :671-677 and :893-914 carry a
`wk_fpi = 1/max(fp, sheath_weak_fpmin)` factor on `wk_res` and on every `tr_J`/`tr_Jp`
entry. This is ~15 lines, a new namelist flag `sheath_weak_fpmin`, default 0 = current
behaviour.

**Note this is the opposite of every gate tried so far.** Gates *remove* rows where the
demand is high; measured, all of them were harmful ([[sheath-multitype-1plus5-works]]:
`wmin=0.5` took 305 steps to 2). This *strengthens* the row there instead. The whole
campaign has been pulling the wrong direction on this lever.

---

## Idea 2: nothing enforces global charge balance, and with two free surfaces that matters

Per-type totals in the healthy part of the 1+5 run:

    type 1  -1451 A    type 5  -1603 A    type 3  -120.0 A    type 9  -147.9 A
    sum ~ -3322 A

Charge conservation requires `oint j.n dS = 0` over the *whole* wall in steady state. The
diagnostic (`mod_sheath_diag.f90:399-401`) only sums types where the sheath is active, so
the return path through the Dirichlet-`u` wall is invisible — we literally cannot see
whether the books balance.

When one surface is free, the Dirichlet wall is an infinite-capacity current sink and the
sheath just picks its local potential. When two are free, the split between them is set by
a *global* compatibility condition that the purely local j-V rows never state. The
constant mode of `u` is then determined only by the residual Dirichlet types (2/3/4/9 —
corners and the leg-parallel surface), and that is a weak, badly placed anchor.

Two things follow, both cheap:

**2a. Measure it.** Accumulate `zj0 * B.n * dS` over *every* boundary element, not only
sheath-active ones, and print the whole-wall total each output. If it drifts monotonically
alongside the `ePhi/kTe` ramp, ideas 1 and 2 are the same phenomenon seen from two sides.
If it stays near zero, idea 2 is out and idea 1 stands alone. ~20 lines in
`mod_sheath_diag.f90`.

**2b. Close the circuit.** `sheath_V_wall_at(BigR)` already threads a wall potential
through the entire chain (`mod_sheath_bc.f90:227`). Make its *scalar* part adapt between
steps:

    V_wall(n+1) = V_wall(n) - gain * I_wall_total / (dI/dV)

with `dI/dV` estimated from the already-accumulated `sum(zj_sat*fp*dx_du*dS)`. This is a
zero-dimensional wall circuit model (Artola-style), needs **no new matrix structure at
all** — one MPI reduction and one scalar update per step — and it removes the free
constant mode by construction. It is also physically the right object: a real vessel is
one conductor, and the two targets are electrically connected through it.

---

## Idea 3: the demand is a cold-leg resistivity artefact, and that is fixable

The characteristic can pass `f ~ 1.45` at `s = 0.03`; the inner target was demanding
`f = 26`. Everyone has treated this as a BC problem. It is a `zj0` problem.

[[cold-leg-clamps]]: `T_min` freezes `eta` below ~5 eV, so at 1 eV the resistivity is
**11x too small**. Too-small `eta` in the leg means too-large parallel current arriving at
the target. Meanwhile `j_sat = c_sat*rho*g*cs/|B|` is small there because the leg is cold
and thin. Both errors push `|j0/j_sat|` up, and both are in the plasma, not the sheath.

`T_min_eta` exists on `keep-current-prof-cutoff` and is **not on this branch** — I checked;
only `T_min_sheath` is here, and its docstring references `T_min_eta` as its model. Porting
it is a small, self-contained cherry-pick, and it attacks the numerator of the one ratio
that every failure mode runs through.

This is worth doing *regardless* of ideas 1-2: `sheath_sat_slope 0.03 -> 0.3` (the lever
memory currently recommends) makes an unphysical demand *reachable*, which is a workaround.
Fixing `eta` makes the demand physical.

---

## Idea 4: staged continuation in space, not just in time

`sheath_ramp_time` is global and keyed on `t_start`
(`mod_boundary_matrix_open.f90:459`). Every multi-type run therefore switches both
surfaces on simultaneously, from a restart converged with neither.

Make it **per type**: `bcs(i)%sheath_ramp_time`, so one can start from the converged
type-1-alone state (which exists and is excellent — weak 6.5e-4, closure 0.00%) and ramp
type 5 in over hundreds of steps while type 1 stays fully on. Since `09c96aa39` the ramp
scales the residual only and leaves the `u` coupling intact, so a partially-ramped surface
is still constrained — the failure mode that made the ramp useless before is gone.

This directly targets "multi-type dies while still evolving": give each surface time to
find its own steady state rather than asking the solver to find both at once from a state
consistent with neither. `bcs%` already carries per-type sheath flags, so this is a
`type_bcs` field, a broadcast entry, and one line at :459-460.

---

## Idea 5: which boundary `zj` DOFs got no row at all — currently unknowable

`mod_sheath_trace.f90` header admits the MPI limitation, and `sheath_trace_apply` is
purely rank-local: a row is written only if `st_row(is)` lies in `[index_min, index_max]`.
A boundary trace DOF **owned** by a rank that has no local boundary element touching it
gets accumulated by the neighbour (which then discards it) and written by nobody. With
`dirichlet%zj = .false.` and `natural%zj = .false.` that DOF's equation is the volume weak
form *missing its surface term* — incomplete, exactly the runaway in
[[sheath-d-collapse-crash]].

Multi-type multiplies the number of partition-boundary trace DOFs, and the effect is
partition-dependent, which would explain why no parameter story ever accounted for it.

**Diagnostic (do this before any fix):** count boundary `zj` DOFs in `[index_min,
index_max]` on a sheath-enabled type that received no `sheath_trace_add`. If it is 0, this
idea is dead in one run. The `416 accumulated / 372 replaced` figure does **not** answer
this — the 44 are most likely halo duplicates the owner also wrote.

**If it is non-zero, the cheap safety net already exists:** `bcs(i)%natural%zj = .true.`.
The `mod_sheath_bc.f90` header says so explicitly ("the reason to turn natural%zj on if a
boundary zj row ever falls below the trace accumulator's degeneracy floor"), memory lists
it as an immediate lever, and I find no evidence it was ever actually run. It gives every
unreplaced row a complete equation and costs nothing on rows that *are* replaced, because
those are overwritten with `zbig`.

---

## Idea 6: a diagnostic that would have settled this months ago

There is no measurement anywhere of the row's **authority over `u`**. Add, per boundary
type, the area-integrated sheath conductance

    G = integral of  |zj_sat * fp * dx_du| dS      [the total dI/dV of that surface]

plus `fp` min / area-mean. This is the single number ideas 1 and 2 both turn on, it costs
one accumulator in `sheath_diag_add`, and it is *predictive* rather than post-hoc:

- idea 1 right => `G_inner` collapses ~200x over the 305 steps while `G_outer` is flat,
  and it starts falling **before** `max|j/jsat|` moves;
- idea 1 wrong => `G` is steady and the blow-up is genuinely a `zj0` transport event
  (-> idea 3).

Unlike `|Fd/S|`, which [[sheath-multitype-1plus5-works]] confirmed is useless as a
discriminator (type 1 read 1.77e3 while the dying type 4 read 2.32e2), `G` has a
mechanism attached to it.

---

## Two smaller code observations

**`wk_dfac` is computed and never used.** `mod_boundary_matrix_open.f90:672-677` sets it,
its declaration at :86 calls it "the factor every Jacobian column carries", and no column
carries it — the `tr_J` entries at :893-914 take `wk_wrx` only. As a trust region this is
arguably the *right* choice (bounded `F` with a full-strength `J` gives a smaller step,
whereas applying `dfac` would shrink `J` and grow it), but it means `J != dF/dx`, so the
Newton rate is degraded whenever the cap is active — with `rmax = 2` and a demand of 26,
`den = 14` and `dfac = 1/196`, i.e. the linearisation is off by two orders. Either use it
and accept the larger step, or fix the comment. Right now the code and its documentation
disagree about which is happening.

**`sheath_sat_slope_e` is the wrong tool and should probably be retired.** Its own comment
says it makes `f` unbounded below; memory records it as "an enabler, not a stabiliser".
Idea 1's `fp` floor achieves what it was reaching for — nonzero `u` authority deep on a
branch — *without* moving the solution, which `sat_slope_e` does by construction because it
alters `f`.

---

## Suggested order

1. `sheath_dfdx_min = 0.3` on the 1+5 305-step baseline. One line, one run, decisive
   either way. (idea 1)
2. Add the `G` / `fp` diagnostic and the whole-wall current sum. Two small patches, and
   re-run the same case. (ideas 6, 2a)
3. Whichever of `sheath_weak_fpmin` (row in potential space) or the `V_wall` circuit the
   step-2 measurement points at. (ideas 1, 2b)
4. Port `T_min_eta`. Independent of the above and worth doing anyway. (idea 3)
5. Per-type `sheath_ramp_time`, to start multi-type runs from the converged type-1 state.
   (idea 4)

Idea 5's diagnostic is cheap enough to fold into step 2.

---

# ADDENDUM 2026-09-03: what type 1 actually IS

Read off `grids/grid_xpoint_wall.f90` and `matrix/construct_matrix_mod.f90`. This supersedes
the run-outcome reasoning above as the primary lead.

## 1. Type 1 is a bitfield, and it is the only sheath type born in the flux-aligned grid

Types 1/2/3 are never assigned literally. They are accumulated from two flags
(`grid_xpoint_wall.f90:1116-1122`, and again at `:1174-1180` for the leg):

    boundary = 0
    if (k == 1 .or. k == n_open+n_private+1) boundary = boundary + 2   ! radial extreme
    if (j == 1)                              boundary = boundary + 1   ! poloidal end = TARGET

so  1 = target only,  2 = radial extreme only,  3 = both = the corner. That is the
"1: TARGET side 2 / 2: TANGENT side 3 / 3: CORNER between them" recap, generated.

Types **4, 5, 9 are assigned literally, and only inside `if (extend)`** — the ray-cast
extension from the plasma edge out to the wall polygon (`:1258, :1322, :1326, :1385, :1389`
and `:1633-1642, :1682-1692`). **Type 1 and type 5 come from two different mesh
generators.** No amount of namelist symmetry makes them the same kind of object.

## 2. The node frames are built from different things

Flux-grid node (types 1/2/3), `:1113-1114`:

    x(1,2,:) = (dR_dt, dZ_dt)   / |..|      ! grid-line tangent, from the CUB1D spline
    x(1,3,:) = (-PSI_Z, +PSI_R) / |grad psi|  ! EXACTLY the poloidal field direction

Extension node (types 4/5/9), `:1251-1252` and `:1310-1312`:

    x(1,2,:) = (dRtmp, dZtmp) / |..|              ! the ray-cast direction to the wall
    x(1,3,:) = (cos(tht_ext), sin(tht_ext))       ! an interpolated GEOMETRIC ANGLE
                                                  ! (for the leg: atan2 to the neighbours)

`x(1,2,:)` and `x(1,3,:)` are the directions of the two first-derivative DOFs. So on a
type-1 node the second derivative DOF is **aligned with B_pol**; on a type-4/5/9 node it is
a geometric angle with no relation to the field or to psi.

## 3. And that is exactly the DOF the sheath closure runs through

`mod_boundary_matrix_open.f90:178`:

    direction_perp(1) = 6 / direction(2)          ! =3 if direction(2)=2, =2 if direction(2)=3

`construct_matrix_mod.f90:135-168` gives `direction(2) = 2` for types 1, 4, 9-9 and
`= 3` for types 2, 5. Therefore:

| | along-edge DOF | NORMAL-derivative DOF | what that direction IS |
|---|---|---|---|
| type 1 | 2 | **3** | **B_pol / flux-surface tangent** |
| type 5 | 3 | **2** | the ray-cast direction to the wall |

`dirichlet%psi` pins the value and the TANGENTIAL derivative and leaves the normal one free
— and `mod_sheath_bc.f90`'s header says in as many words that this free normal DOF "is the
degree of freedom through which the wall current responds to the sheath".

**On type 1 that free DOF is `dpsi/dl` along B_pol, which is identically ZERO** — psi is
constant along a flux surface. The whole of `grad psi` sits in the PINNED tangential DOF.
The quantity the sheath has to move starts at exactly zero, with no background to hide in.

**On type 5 the free DOF is the ray-cast direction, which crosses flux surfaces**, so it
carries `~|grad psi|` — large. The sheath's response is a small relative perturbation on a
big number, while the pinned DOF is the geometric angle. Same code, same flags, inverted
conditioning.

This is a difference of KIND and it is specific to this boundary condition, because this
is the only BC whose closure goes through the free psi normal-derivative DOF. It also
explains, without any new physics, why `sheath_psi_jacobian` (`04b79e59f`, columns on
`direction_perp`) never moved the needle: on type 1 it is a column on a DOF that is zero.

**Prediction:** the failure severity should scale with `|grad psi|` at the boundary, i.e.
with how steeply the ray-cast extension crosses flux surfaces. Measurable directly.

## 4. Type 4 is type 1, relabelled, but made of extension nodes

`:1641-1642` and `:1691-1692`:

    newnode_list%node(...)%boundary = 4  ! 1   <- the trailing comment is the OLD value

Someone changed these from 1 to 4 (and the corners from 3 to 9). So
`update_boundary_types.f90:627` "4: TARGET, side 2 (Same as type-1)" is TRUE at the level of
the label and of `direction()`. What it is NOT is the same at the level of the node frame:
these are `index_ext2 + ...` nodes, i.e. extension nodes with the `(cos tht_ext, sin tht_ext)`
frame. **The "open contradiction" recorded in [[sheath-multitype-1plus5-works]] is resolved:
type 4 is type 1's label on type 5's geometry.** That is precisely the worst combination —
`direction(2) = 2` selects DOF 3 as the normal derivative, and on an extension node DOF 3 is
the interpolated angle, not B_pol.

## 5. The frame conditioning is never checked anywhere

The two frame vectors are not required to be orthogonal, or even independent. On flux nodes
they nearly are (a target cuts across the field). On extension nodes DOF 2 is the ray-cast
direction and DOF 3 is the direction to the neighbouring extension nodes along the wall —
and at a leg end those can become nearly PARALLEL. There is one guard, at `:1319`, and it
only fires for a ray of exactly zero length:

    if ((dRtmp**2 + dZtmp**2) .lt. 1d-8) x(1,2,:) = (/ -sin(tht_ext), cos(tht_ext) /)

Nothing checks near-degeneracy. A near-singular nodal frame makes every derivative DOF at
that node blow up, and `zj = Delta psi` amplifies it by `1/h^2` — which is exactly the
`D min` collapse of 8-10 orders seen on types 4 and 9 and never on 1 or 5.

### THE MEASUREMENT TO RUN FIRST

Per boundary type, print the frame determinant and the psi alignment:

    det_a = | x(1,2,:) CROSS x(1,3,:) |            ! 1 = orthonormal, 0 = degenerate
    alp_a = | grad(psi) . x(1,direction_perp(1),:) | / |grad psi|   ! 0 on an ideal type 1

Expected if this is the story: `det ~ 1, alp ~ 0` on type 1; `det` with a small minimum and
`alp ~ 1` on types 4 and 9; type 5 in between. This is ~30 lines in `mod_sheath_diag.f90`,
needs no physics change, and it discriminates in a single step of a single run — unlike
`|Fd/S|`, `W`, `b_n` and `dS`, all of which have already been shown not to.

If it holds, the fix is not a sheath parameter at all. It is either
 (a) rebuild the extension node frames so DOF 3 follows `(-PSI_Z, PSI_R)` like the flux
     grid does — a few lines in `grid_xpoint_wall.f90`, but it changes the grid, so it needs
     a fresh equilibrium rather than a restart; or
 (b) restrict the weak sheath to types whose frames are field-aligned (1, and 5 where it
     measures clean), and state the uncovered area honestly — which is what the 1+5 result
     is already doing by accident.

## 6. Bookkeeping note

`use_simple_bnd_types = .false.` by default (`preset_parameters.f90:826`) and
`max_bnd_types = 30`, so types 11/15/19/20/21 are addressable. The logs only ever show
1/3/4/5/9, which means either the namelist sets it `.true.` or `update_boundary_types` is
not on this grid's path. **Worth confirming from the namelist**: if 11/15 exist as distinct
node types they are the OUTWARD-field halves of the type-1 and type-5 targets, they carry
real area, and they have never had a `bcs()` block. Cheap to check, and expensive to be
wrong about.

---

# DEEP DIVE 2026-09-03: the grid construction

All line numbers `grids/grid_xpoint_wall.f90` unless stated. Confirmed path: with
`grid_to_wall = .true.`, `n_wall_blocks = 0`, `xcase < UPPER_XPOINT` and
`RZ_grid_inside_wall = .false.`, `core/mod_flux_grid.f90:56-63` calls `grid_xpoint_wall`,
which carries the comment `!!rks only for ITER wall for the moment`.

## A. PROVEN: the psi normal-derivative DOF is analytically zero on type 1 and only there

`:1886-1889` sets every node's psi DOFs as the projection of grad(psi) on the node frame:

    values(1,2,psi) = PSI_R*x(1,2,1) + PSI_Z*x(1,2,2)
    values(1,3,psi) = PSI_R*x(1,3,1) + PSI_Z*x(1,3,2)

and on a flux-grid node `x(1,3,:) = (-PSI_Z, +PSI_R)/|grad psi|` (`:1114`), so

    values(1,3,psi) = ( PSI_R*(-PSI_Z) + PSI_Z*(+PSI_R) ) / |grad psi| = 0    EXACTLY

Now chain it to the boundary condition. `mod_boundary_conditions.f90:251-260` sets
`iv_dir`/`iv_perp_dir` from the element side, and `:400-405` ("Fix derivatives in one
direction") pins the value and `iv_dir` only - **the normal derivative `iv_perp_dir` is
left free**. `mod_sheath_bc.f90`'s header states that this free normal DOF is "the degree
of freedom through which the wall current responds to the sheath". And
`mod_boundary_matrix_open.f90:178` with `construct_matrix_mod.f90:135-168` gives

    types 1, 4, 9  ->  tangential DOF 2, FREE normal DOF 3
    types 2, 5     ->  tangential DOF 3, FREE normal DOF 2

**On type 1 the free DOF is `values(1,3,psi)` = exactly zero at construction.** All of
grad(psi) sits in the PINNED DOF. The quantity the sheath has to move starts at zero with
no background to hide in.

**On types 4/5/9 the free DOF carries O(|grad psi|)** - `x(1,3,:)` there is
`(cos(tht_ext), sin(tht_ext))`, an interpolated geometric angle (`:1252`, `:1311`), and
`x(1,2,:)` is the ray-cast direction to the wall. Neither is a null direction of psi. The
sheath's response is then a small relative perturbation on a large number.

This is no longer an inference from run outcomes. It is two lines of the grid builder and
two lines of the boundary-condition loop.

Corroborating detail: `:1891` explicitly forces `values(1,3,psi) = 0` for boundary type 2
- which is already zero analytically, so it only cleans up interpolation roundoff. The
author knew the property mattered and enforced it. It is applied to type 2 and to nothing
else; types 4/5/9 never get it, and could not (their DOF 3 is not a null direction).

## B. PROVEN: the element size assumes the frame is parallel to the edge, and nothing checks it

`:1751-1752`:

    size_0 = sign( sqrt((R0-RP)**2+(Z0-ZP)**2)/n_order ,  dR0*(RP-R0) + dZ0*(ZP-Z0) )

The MAGNITUDE is the chord length; the frame vector enters through its SIGN only. This is
exact iff the frame vector is parallel to the element edge. At angle `theta` the
isoparametric map is distorted by `tan(theta)`, and as `theta -> 90 deg` the sign is set by
noise. `iv = 1,3` uses `x(1,2,:)`, `iv = 2,4` uses `x(1,3,:)`.

Where does each type sit?

* **type 1** - both frame vectors ARE the analytic tangents of the two grid-line families
  (the CUB1D tangent, and B_pol which is tangent to the flux surfaces the other family
  follows). `theta ~ 0` to O(h^2). The element is a proper isoparametric map.
* **type 5** - `x(1,2,:)` is the CUB1D tangent of the polar ray, and at `i = n_ext`
  `tht_ext = T_wall_par(j)` EXACTLY (`:1249`), i.e. the wall tangent, and the nodes lie on
  the wall. So `theta ~ 0` on the type-5 boundary too. **That is why type 5 works.**
* **type 4 / 9** - in the LEG block `tht_ext` is a secant through the neighbouring
  extension nodes, and at the ends it degrades to a ONE-SIDED secant (`:1288-1300`:
  `j > 2` and `j < n_leg` fall back to `Rtmp2 = Rtmp` / `Rtmp3 = Rtmp`). Type 4 and type 9
  are exactly the leg-end nodes.

**Nothing anywhere measures the frame-chord angle, the frame determinant, or the frame
orthogonality.** The only guard (`:1319`) fires solely for a zero-length ray.

## C. PROVEN: `tht_ext` is computed twice in the leg block and the first value is dead

`:1286` computes the blend from the flux-surface tangent to the wall tangent - the same
construction the main-wall block uses - preceded by four lines of careful 2*pi
branch-fixing (`:1282-1284`). `:1309` then overwrites it unconditionally with
`atan2(Ztmp2-Ztmp3, Rtmp2-Rtmp3)`. **The main-wall extension and the leg extension use two
different frame conventions, and the leg one is arrived at by silently discarding the
other.** Nobody can tell from the source which was intended.

## D. PROVEN: type 4 mixes two frame conventions inside one boundary type

`:1641-1642` relabels `vertex(3)` and `vertex(4)` to type 4 for `j = n_leg-1, i != n_ext`
(the trailing comment `! 1` is the OLD value, so these were type 1). At `i = 1` the element
vertices are (`:1611-1614`)

    vertex(1), vertex(4) = index_leg1 + ...    <- FLUX-GRID leg nodes  (frame: tangent, B_pol)
    vertex(2), vertex(3) = index_ext2 + ...    <- EXTENSION nodes      (frame: ray, secant)

so the relabel puts **one flux-grid node and n_ext-1 extension nodes into the same boundary
type**, meeting across a shared edge. Their DOF-3 directions mean different things and
their `values(1,3,psi)` differ by O(|grad psi|) across that edge.

## E. PROVEN: the type-9 surgery rewrites types without touching frames

`:1923-1949` removes the "small edge triangles":

```fortran
do i = 1, n_elements
  do j = 1, n_vertex_max
    iv = element(i)%vertex(j)
    if (node(iv)%boundary .eq. 9) then
      remove_elements(...) = i ;  remove_nodes(...) = iv
      node(iv2)%boundary = 99 ;  node(iv3)%boundary = 99 ;  node(iv4)%boundary = 99
    endif
  enddo
enddo
```

Four separate problems, in descending order of how much I would bet on them:

1. **The promoted nodes never get a frame suitable for the boundary they now sit on.**
   Removing the quad creates a NEW boundary edge, and no `x(1,·,:)` is recomputed to be
   tangent to it. The promotion itself is necessary - those nodes really do become
   boundary - but combined with B it means type 9's frame-chord angle is essentially
   arbitrary. That is the worst case of B, and it lands precisely on the type that dies in
   4 steps.
2. **`= 99` is unconditional, so a genuine type-5 wetted wall node adjacent to a removed
   triangle is DEMOTED to 9** and silently loses its type - and with it its `bcs()` block.
   Real wall area changes classification as a side effect of a mesh repair.
3. **No `exit` after a hit.** An element with two type-9 vertices is pushed onto
   `remove_elements` twice and removes two nodes.
4. **Order dependence.** The pass reads `.eq. 9` while writing `99` into other nodes, so a
   type-9 node that is overwritten before it is visited never has its own element removed.

## F. Type 4's chain is cut by that surgery

`:1636-1637` sets `vertex(3) -> 9` and `vertex(4) -> 4` on the element at
`(j = n_leg-1, i = n_ext)`. That element has a type-9 vertex, so E deletes it, and its
other three vertices - the type-4 one and two type-5 ones - all become `99 -> 9`. So the
type-4 run of nodes at `i = 1..n_ext-1` **terminates against a hole**, at nodes whose type
AND frame were both rewritten. Combined with D, type 4 is: one flux-grid node at one end,
extension nodes in the middle, and a surgery-rewritten node at the other end. Three
conventions in 45 nodes.

## G. Minor, conditional

`construct_matrix_mod.f90:137-144`: in the `(9,9)` branch `direction(3)` is never assigned,
while every other branch assigns it. `direction` is declared at `:48` in the enclosing
scope, so it would carry over from the previous boundary element. **Dead at the current
`n_order = 3`** (`mod_settings.f90:11`, so `n_degrees_1d = 2` and the `n_order >= 5` guards
never fire) - but it is a live bug for anyone building at higher order.

## Ranking, and what it says about the run results

| type | frames | free psi DOF | surgery | observed |
|---|---|---|---|---|
| 1 | analytic tangents, both | **exactly 0** | no | 3900, converged |
| 5 | wall tangent at i=n_ext, clean | O(\|grad psi\|) | edges only | 305-308, slow drift |
| 4 | THREE conventions in one type | O(\|grad psi\|) | cut at one end | 8 |
| 9 | never recomputed after surgery | O(\|grad psi\|) | is the surgery | 4 |

The ordering of the failure severity matches the ordering of the grid pathology exactly,
and it separates the two failure MODES that the run logs already distinguish: type 5's is a
slow monotone drift (a conditioning problem - B is fine, A is not), types 4 and 9 blow up
in single-digit steps (a correctness problem - B and E are both violated). That is the
"difference of kind" the run data demanded and no namelist hypothesis supplied.

## THE TEST - no rebuild, no cluster

The HDF5 restart already carries everything (`mod_export_restart.f90:615-703` writes the
`x`, `values`, `boundary`, `vertex` and `size` datasets), so:

    python3 util/check_boundary_frames.py jorek000000.h5

prints, per boundary type:

* **cos** - `|frame . chord| / |chord|` on every boundary element side, for the DOF that
  `:1751` actually uses to build `element%size`. This is the direct measure of finding B:
  1.0 means the size construction is exact. A low MIN with a high MEAN is a localised bad
  patch (surgery, or a one-sided secant); a low mean is a systematically misaligned type.
* **det** - `|x2 x x3|`, the frame determinant (finding C). Below ~0.3 the frame is near
  degenerate and every derivative DOF there is inflated, with `zj = Delta psi` amplifying
  by `1/h^2`.
* **free frac** - the fraction of `grad(psi)` in the FREE normal DOF (finding A), using the
  per-type tangential/normal mapping. ~0 on type 1 by construction; O(1) elsewhere.

Verified end to end against synthetic files with a deliberately misaligned, squashed node,
in both axis layouts: clean rows read `cos 1.0000 / det 1.0000 / free 0.0000`, the bad node
reads `cos 0.6600 / det 0.0639 / free 0.6247`.

**HDF5 axis-order gotcha**, measured on a real file: the wrapper does NOT preserve the
Fortran argument order. `x`, declared `(n_nodes, n_coord_tor, n_degrees, n_dim)`, lands as
`(n_nodes, n_dim=2, n_degrees=4)` - trailing dims reversed - and `values`, declared
`(n_nodes, n_tor, n_degrees, n_var)`, lands as `(n_nodes, n_var=8, n_degrees=4)`, i.e.
`values[:, var, deg]` with psi = var 0. The script detects the axes by size and prints the
raw -> resolved mapping, so a layout surprise shows up as a printed line, not a wrong table.

## What to do if it confirms

**Do not start by changing the grid.** Order by cost:

1. **`bcs(4)` and `bcs(9)` stay off.** Already the 1+5 configuration. The honest statement
   in a paper is now much stronger than "types 4 and 9 are unstable": they are a leg-end
   patch with three frame conventions and a mesh repair that rewrites boundary types
   without recomputing frames. ~1 % of the wetted area.
2. **Fix E cheaply and correctly** - add `exit` after the hit, and guard the `= 99` write
   so it cannot demote an existing wetted type. Both one-liners, neither changes node
   positions, so no new equilibrium is needed. Then re-measure with the script.
3. **Add the assertion that should always have been there**: in the size loop at `:1751`,
   warn when `|frame . chord| / |chord| < 0.5`. Five lines, permanent, and it would have
   caught all of B, D and E at grid-build time.
4. **Normalise the psi column of the weak sheath row by the local `|grad psi|`.** This is
   the only change that addresses A without touching the grid: `sheath_psi_jacobian`
   currently writes a column on `direction_perp`, whose scale differs by orders between a
   type-1 node (DOF ~ 0) and a type-5 node (DOF ~ |grad psi|). That is a per-type scale
   error in the Jacobian, and it is why the flag "8 -> 9 steps alone" never discriminated.
5. **Only then**, if a regrid is on the table: give extension nodes a field-aligned DOF 3
   like the flux grid. `PSI_R`/`PSI_Z` are already computed at every node at `:1880-1882`,
   so the data is there. But the frame plays TWO roles - geometry parametrisation (B) and
   field-DOF basis (A) - and the extension needs different things from each: wall-tangent
   for B, flux-tangent for A. **You cannot have both from one vector.** That trade-off, not
   a missing line, is the real design problem, and 4 is the way to sidestep it.


---

# MEASURED 2026-09-03: the frame determinant is the discriminator

Run on the real grid (`baseline_cutoff_github_jv_weak_1_sat_slope0.03_novpar_node15`,
`jorek000000.h5`, 48292 nodes / 47829 elements).

**Validation first:** the script reports 98 type-1 and 362 type-5 boundary nodes. Times two
trace DOFs that is **196 and 724** - exactly the per-type trace-row counts in the run logs
([[sheath-multitype-1plus5-works]]). It is reading the right objects.

    type  nodes | cos min  cos mean | det min  det mean  |x2.x3|max | |v3|/|v| | free frac
       1     98 |  0.9998   1.0000  |  0.7152   0.8815     0.6989   |  0.0000  |  0.0000
       2     65 |  0.9999   1.0000  |  0.9383   0.9674     0.3458   |  0.0000  |  1.0000
       3      2 |  1.0000   1.0000  |  0.9344   0.9502     0.3563   |  0.0000  |     -
       4     88 |  1.0000   1.0000  |  0.0187   0.3782     0.9998   |  0.4505  |  0.4505
       5    362 |  0.8584   0.9993  |  0.0810   0.8745     0.9967   |  0.2163  |  0.9412
       9      6 |  1.0000   1.0000  |  0.0092   0.0511     1.0000   |  0.6892  |  0.6892

## (A) The psi-DOF prediction is CONFIRMED, exactly

`|v3|/|v| = 0.0000` on types **1, 2 and 3** - every flux-grid type - and nonzero on **4, 5
and 9** - every extension type. The analytic argument (`values(1,3,psi) = grad(psi) .
(-PSI_Z,+PSI_R) = 0`) is confirmed to four decimals, and confirmed *in the evolved state*,
not just at construction: this is a mid-run restart, so the solver is maintaining it.

**But it is NOT the discriminator, and I over-weighted it.** Type 2's free-DOF fraction is
1.0000 and type 2 is perfectly healthy; type 5's is 0.9412 and it runs 308 steps; type 9's
is only 0.6892 and it dies at 4. The ordering does not follow. A is real, it explains why
`sheath_psi_jacobian` writes a column whose scale differs by orders between types, and it
is worth fixing - but it is not what kills types 4 and 9.

## (B) The frame-chord / element-size prediction is FALSIFIED

`cos` is ~1.0000 on every type, and the *worst* value in the whole grid (0.8584) is on type
**5**, the type that works. The `size = sign(|chord|, frame.chord)` construction at
`grid_xpoint_wall.f90:1751` is therefore exact everywhere in practice: the frames really are
parallel to the edges they parametrise, including on the leg extension and on the surgery
nodes. **Drop this hypothesis.** The one-sided-secant argument (finding B of the deep dive)
predicted types 4 and 9 would be worst here; they read 1.0000.

## (C) The frame DETERMINANT is CONFIRMED and is the discriminator

`det = |x2 x x3| = |sin(angle between the two frame vectors)|`, so the nodal derivative
basis is conditioned as `1/det`:

    type   det min   angle    1/det   det mean   mean angle   steps survived
       1   0.7152   45.66d      1.4     0.8815      61.82d    3900 (timeout)
       2   0.9383   69.77d      1.1     0.9674      75.33d    - (never sheathed)
       3   0.9344   69.13d      1.1     0.9502      71.84d    -
       5   0.0810    4.65d     12.3     0.8745      60.99d     308
       4   0.0187    1.07d     53.5     0.3782      22.22d       8
       9   0.0092    0.53d    108.7     0.0511       2.93d       4

**`det mean` orders the four sheath-capable types exactly as their survival**, across three
orders of magnitude in step count, and it separates them into two clean groups: 0.88 / 0.87
(types 1 and 5 - both work) against 0.38 / 0.05 (types 4 and 9 - both fail alone). On type 9
the two frame vectors are parallel to four decimals (`|x2.x3| = 1.0000`, mean angle 2.93
degrees) over the WHOLE type, not just a tail.

**Mechanism, and it is geometric.** On an extension node `x(1,2,:)` is the ray from the
plasma edge to the wall and `x(1,3,:)` is the along-wall/secant direction. `det -> 0` means
**the ray to the wall runs nearly ALONG the wall** - i.e. the extension meets the wall at
grazing incidence. That is exactly what happens at the leg-end corner, which is where types
9 and the bad half of type 4 live. The nodal derivative basis there is ill-conditioned by
50-110x, every derivative DOF at those nodes is inflated by that factor, and `zj = Delta
psi` is built from SECOND derivatives. The weak sheath row REPLACES the zj row at precisely
those DOFs, so it is writing a constraint on a numerically degraded quantity.

**Type 4 is bimodal, as predicted** - `det` min 0.0187 against mean 0.3782 over 88 nodes,
consistent with the deep dive's finding D/F (one flux-grid node, extension nodes, and a
surgery-rewritten node in a single boundary type).

**Type 5 has a bad tail, not a bad bulk** - min 0.0810, mean 0.8745. That is the natural
reading of its failure mode: it runs 308 steps with a slow drift rather than blowing up in
single digits, because only a few of its 362 nodes are degenerate.

## What to do

**1. Gate the weak sheath row on the frame determinant.** New flag `sheath_weak_detmin`
(default 0 = current behaviour), skipping the replaced row on a node whose `det` is below
it and falling through to the frozen-increment fallback `da8a4b555` already provides.

This is **not** another version of `sheath_weak_wmin`, and the difference matters: every
gate tried so far keyed on the SOLUTION (`|zj0/zj_sat|`, `D_a/D0_a`), so it could close in
response to the very divergence it was meant to prevent - which is why
[[sheath-multitype-1plus5-works]] measured `wmin = 0.5` taking a 305-step case down to 2.
A determinant gate keys on a **static geometric property of the mesh**, fixed for the whole
run and identical on every rank. It cannot participate in a feedback loop, and which nodes
it removes is knowable before the run starts - from the table this script prints.

Get the per-type `%` of nodes below the threshold from the updated script before choosing
one; `det < 0.3` is the natural first cut (type 1 loses nothing, type 9 loses everything).

**2. Then re-test `bcs(4)` and `bcs(9)`** with the gate on. This is the first proposal in
the whole campaign that has a measured, run-independent reason to expect a different answer.

**3. The grid fix, if a regrid is ever on the table:** orthogonalise the extension frames -
Gram-Schmidt `x(1,3,:)` against `x(1,2,:)` at `:1252` and `:1311`, or detect `det` below a
threshold and rotate. Note this changes `element%size` through `:1751` and therefore changes
the geometry, so it needs a fresh equilibrium, not a restart. Also note (B): the current
frames are chord-aligned, and that property is worth not breaking.

**4. Still worth doing, independently:** normalise the psi column of the weak sheath row by
local `|grad psi|` (finding A). Confirmed real, just not the killer.


---

# MEASURED, part 2: the degenerate nodes are TWO CORNERS, and the inner one is worse

Distribution and locations from the same file.

    type  nodes | det min  det p05  det med  det mean | %<0.3  %<0.1 | nodes <0.3
       1     98 |  0.7152   0.7952   0.8893    0.8815 |  0.0%   0.0% |    0 of 98
       2     65 |  0.9383   0.9497   0.9619    0.9674 |  0.0%   0.0% |    0 of 65
       3      2 |  0.9344   0.9359   0.9502    0.9502 |  0.0%   0.0% |    0 of  2
       4     88 |  0.0187   0.0475   0.3612    0.3782 | 42.0%  14.8% |   37 of 88
       5    362 |  0.0810   0.5769   0.9384    0.8745 |  2.2%   0.3% |    8 of 362
       9      6 |  0.0092   0.0107   0.0464    0.0511 |100.0%  83.3% |    6 of  6

51 of 621 boundary nodes are degenerate, and they sit at exactly **two points**:

    cluster A   R 1.2654 .. 1.2678   Z -1.0011 .. -0.9914     2.4 mm wide
    cluster B   R 1.6050 .. 1.6057   Z -1.1084 .. -1.1060     0.7 mm wide

Both are leg-end corners, and **types 4, 5 and 9 interleave inside each one** - which is the
direct fingerprint of the surgery at `grid_xpoint_wall.f90:1923-1949` rewriting boundary
types in a small neighbourhood (findings E and F of the deep dive). `det` rises monotonically
with distance from each corner, exactly as the "ray to the wall runs along the wall" picture
requires.

## The inner corner is 1.65x worse - and every failure in this campaign is on the inner target

With `sheath_diag_R_split = 1.42`, cluster A is the INNER leg end and cluster B the OUTER:

    cluster A  INNER   worst det 0.0092   1/det 108.7   10 of the worst 15 nodes
    cluster B  OUTER   worst det 0.0152   1/det  65.8    5 of the worst 15 nodes

Now compare against what the logs have always said and never explained:

* [[sheath-d-collapse-crash]]: "a **LOCAL collapse on the INNER target**, not a global
  drift"; INNER `max|j/jsat|` 19.2 -> 88.7 -> 9.6e4 while OUTER stays 11.6-12.5 throughout.
* [[sheath-multitype-working-config]]: the 220-step 1+4 failure is on the inner target;
  "the OUTER target is untouched throughout" (`ePhi/kTe` 1.05 -> 1.06).
* [[sheath-multitype-1plus5-works]]: the 305-step drift is the INNER `ePhi/kTe` ramping
  3.96 -> 10.02 with the outer flat.

**Every documented failure is on the inner target, and the inner corner is measurably the
more degenerate one.** That correlation was not designed into this test - the script has no
knowledge of which target fails - so it is genuine corroboration.

## Why type 1 alone survives despite the corners being in its grid

The corners are static: they are in the mesh for every run, including the 3900-step type-1
one. But type 1 has **0.0 %** of its nodes below 0.3, so with only `bcs(1)` enabled no
replaced row is ever WRITTEN on a degenerate node - those nodes keep their ordinary
Dirichlet. Enable 4, 5 or 9 and the weak route starts writing a constraint on
`zj = Delta psi` at nodes whose derivative basis is conditioned 50-110x. That is the whole
mechanism, and it explains the "difference of kind" without any parameter.

    type   max 1/det   degenerate nodes   steps survived
       1         1.4          0 of  98      3900 (timeout)
       5        12.3          8 of 362       308
       4        53.5         37 of  88         8
       9       108.7          6 of   6         4

Monotone in both columns. Four points, so do not fit a law to it - but the ranking is
unambiguous and it is a property of the MESH, fixed before the run starts.

## THE NEXT RUN

`sheath_weak_detmin`, gating the replaced row on the frame determinant. The updated script
prints the exact cost of any threshold; from the table above, **0.3** costs type 1 nothing,
removes all 6 type-9 nodes, 37 of 88 type-4 and 8 of 362 type-5.

**Test it on 1+5 first, not on 1+4.** That is the working 305-step configuration, the gate
touches only 8 of its 920 trace rows, and the prediction is sharp and directional: those 8
nodes include the single worst node in the grid, on the inner target, which is where the
drift that ends the run lives. If 1+5 with `detmin = 0.3` runs past 305, the mechanism is
confirmed on the best baseline available and only then is it worth re-enabling 4 and 9.

**One honest caveat, and it is the same one that sank `sheath_weak_wmin`.** Skipping the row
freezes that node's `zj` increment but does NOT restore a boundary condition on `u` - on a
weak type `dirichlet%u` and `natural%u` are both required `.false.`, so u at a gated node is
left to the interior equation and its neighbours. With 8 gated nodes among 354 constrained
ones on type 5 that should be harmless; with 37 of 88 on type 4 it may well not be, and a
fragmented type 4 could reproduce the wmin failure for a different reason.

What is genuinely different from wmin is *when* and *why* the gate closes: wmin keyed on the
SOLUTION and could close in response to the divergence it was meant to prevent, so it fired
mid-run and changed the boundary condition discontinuously during evolution. A determinant
gate keys on a static property of the mesh - fixed from step 1, identical on every rank, and
the exact node set is knowable before the run from the table this script prints. It cannot
feed back. That removes the feedback failure mode; it does not remove the u-constraint one.

If 1+5 improves but 4 and 9 still fail with the gate on, the answer is the grid fix
(Gram-Schmidt the extension frames at `:1252` and `:1311`), which needs a fresh equilibrium
because rotating `x(1,3,:)` changes `element%size` through `:1751` and therefore the mapped
geometry. Note also finding B: the current frames ARE chord-aligned everywhere, and any
frame fix must not break that.


---

# COMPLETION: the gate has to freeze u as well, or 1+4 is untestable

The first cut of `sheath_weak_detmin` only skipped the trace ROW. That leaves a gated node
with **zj frozen and u free** - which on a weak type means u has no boundary condition at
all, because `dirichlet%u` and `natural%u` are both required `.false.` there. Memory records
exactly that combination as the one that produced the boundary current filament
([[sheath-bc-tests-flux-dragging]]). On type 5 it was 8 rows among 354 and probably
harmless; on type 4 it would be 37 of 88, which is why the honest advice was "do not test
this on 1+4" - and a gate that cannot be tested on 1+4 does not address the actual goal.

**A degenerate node is now taken out of the sheath ENTIRELY**, which is the policy the code
already applies at a Dirichlet-u corner (`7f2e44427`), so all three mechanisms already exist:

1. `mod_boundary_conditions.f90` applies the **u Dirichlet** at that node regardless of type.
2. `mod_boundary_matrix_open.f90:1090` skips its trace row, so it keeps its Dirichlet zj row -
   the same `cycle` that already handled Dirichlet-u corners.
3. `sh_ufree` at `:218` sees it as frozen, so **`sheath_weak_ufade`** fades its Gauss points
   out of the free neighbour's row. That is precisely the residual leak `f6ef86fbd` was
   written to fix, and it is why `sheath_weak_ufade = .true.` must be set alongside
   `sheath_weak_detmin` - `initialise_parameters` now prints a NOTE if it is not.

The node is then simply not a sheath node: u frozen, zj frozen, no contradiction, and no
half-constrained DOF. `sheath_frame_det` / `sheath_frame_frozen` live in `mod_sheath_bc.f90`
so all three call sites share one definition and cannot drift apart.

## What this actually buys, stated plainly

* **1+5**: extends the best existing configuration. 8 of 920 rows change. Modest but real,
  and it is the clean single-variable test of the mechanism.
* **1+4**: now a legitimate test. Type 4 keeps 51 of its 88 rows as genuine sheath rows and
  gives up 37 to a Dirichlet corner - which is what the far end of type 4 has effectively
  been all along, except that until now it was being given a replaced row it could not
  support. If 1+4 then runs, the coverage question is answered.
* **1+4+5+9**: type 9 is 100 % gated, i.e. it reverts to its pre-sheath Dirichlet. That is
  not a loss - it is 6 nodes and 0.039 m^2 - and it means the full multi-type namelist can be
  used without type 9 killing it in 4 steps.

The honest limit: this removes 8.2 % of the boundary nodes from the sheath rather than fixing
them. The fix is to orthogonalise the extension frames in the grid builder, which needs a
fresh equilibrium. This makes the existing grid usable now and tells you exactly what it
costs.

## Run order

    sheath_weak_detmin = 0.3
    sheath_weak_ufade  = .true.

1. **1+5**, expect `8 frame-gated` in the log and `det min` 0.081 / 0.715 on types 5 / 1.
   Prediction: past 305 steps.
2. **1+4**, expect `37 frame-gated` on type 4. This is the one that matters.
3. **1+4+5+9** if 2 works.
