# Extending the weak sheath BC beyond boundary type 1 — test and implementation plan

Branch `sheath-jsat-vpar-38ab278` (= 38ab27888 + `sheath_jsat_from_vpar`). Written 2026-09-01
from three independent reviews plus direct verification. Every item is marked
**[V]** verified at source, **[I]** inferred, or **[?]** unresolved.

**The goal is not "cover the whole wall".** It is: (a) make sure the sheath is on BOTH targets,
because it may not be; (b) stop the uncovered wall from injecting artefacts into the PFR; (c)
keep the run alive when a target detaches. HFSHD is not reachable by boundary conditions alone —
the measured sensitivity is 0.0005 per volt on the integrated in-out ratio, so no rearrangement
of the wall reaches the factor ~2 the experiment shows. The value here is removing artefacts that
make the NEGATIVE result untrustworthy.

---

## Stage 0 — free checks, no runs. Do these first; two can invalidate prior results.

**0.1 RESOLVED 2026-09-01: type 1 carries BOTH targets (HFS and LFS), confirmed by the user.**
So the sheath is on both, the measured -23.10 V IS a genuine target-to-target difference, and no
prior number needs withdrawing. Stage 2 (type 5) therefore loses its urgency - it becomes a
completeness item, not a correctness fix. The rest of this section stands.

Original concern, kept for the reasoning: types 1 and 5 are the
SAME physical class - field-crossing target - separated only by which logical element side the
wall presents (`update_boundary_types.f90`; the side map at `construct_matrix_mod.f90:112-118`
sends 1,11,12,4 to side 2 and 5,15,2 to side 3). If your two targets are split between 1 and 5,
the measured -23.10 V "target-to-target" dPhi is (floating target) - (Dirichlet-frozen target)
and every in-out number in this campaign is uninterpretable.
CHECK: `sheath_diag_report` already prints per-type area and an INNER/OUTER split with areas.
Read them off any current log. Near-zero area on one of INNER/OUTER = sheath on one target only.

**0.2 Node count and wetted area per boundary type.** Decides (a) whether type 2 exists - it is
the only remaining `dirichlet%u` gauge if 1/4/5/9 are all covered, and losing the gauge leaves
u's constant mode undetermined; (b) whether type 5 in YOUR grid is target or grazing wall.

**0.3 Are elements already being dropped?** `grep -c "boundary element incoherent"` and
`grep -c "boundary element not included"` on current logs. Non-zero means boundary elements are
already losing ALL their physics (Mach1, sheath heat flux, density flux), silently.

**0.4 Mean and min |b_n| PER TYPE.** **[V]** `mach1_diag_add` bins only by major radius
(`mod_mach1_diag.f90:96-110`) and is called from every `mach1` type, so the quoted mean
|b_n| = 0.0080 is diluted by grazing wall and is NOT the target incidence. Small change: pass the
boundary type and bin on it. Until then, treat every |b_n|-derived number as unattributed.

**0.5 Fingerprint executables with `strings`, never the log.** **[V]** `mod_log_params.f90` gates
the whole sheath parameter dump behind a guard that omits `sheath_zj_weak`, so for this
configuration NO sheath parameters are printed. A log grep reads the namelist echo, i.e. what you
typed, not what was compiled.
    strings <exe> | grep -c sheath_jsat_from_vpar        # b01cb62d0
    strings <exe> | grep -c sheath_diag_weak_accumulate  # a3076fd86
    strings <exe> | grep -c bohm_drift_compatible        # 683298bec

---

## Stage 1 — source fixes that are prerequisites. All small, all independent.

**1.1 Make the validation actually run.** **[V]** `initialise_parameters.f90` guards every
`sheath_zj_weak` check behind `any(natural%u) .or. any(natural%w) .or. any(natural%zj) .or.
any(sheath_zj)` - omitting `sheath_zj_weak`. With this configuration none of them has EVER run:
not the `dirichlet%zj` stop, not the `dirichlet%u` stop, not the gauge check, not the n_tor
refusal. Add `.or. any(bcs(:)%sheath_zj_weak)`. The comment directly above says `sheath_zj` had
to be added for exactly this reason. **This is the highest-value line in the plan**: it converts
a forgotten `dirichlet%u/zj = .false.` on a newly covered type from a silent runtime blow-up into
a setup-time stop.

**1.2 Same guard in `mod_log_params.f90`** so the parameter dump prints (see 0.5), and add a
`sheath_zj_weak bnd types` line.

**1.3 `wk_rat` guard inversion.** **[V]** `mod_boundary_matrix_open.f90` ~`:561-563`: if
`|zj_sat|` underflows the 1e-30 guard, `wk_rat` stays 0 so `wk_wgt = 1` - FULL weight exactly
where the gate should close hardest. Add `else wk_wgt = 0.d0`.

**1.4 `factor` non-negativity.** **[V]** ~`:382-386`: nothing clips
`0.25*(1+tanh(...))**2 - c_3`. A slightly larger `vpar_smoothing_coef(3)` makes it negative,
flipping the sign of `g_bn` and swapping the ion and electron branches on the tangential wall.
`factor = max(factor, 0.d0)`.

**1.5 Headroom.** `st_max_col` 64 -> 96, `st_max_row` 4000 -> 8000
(`mod_sheath_trace.f90:43,52`). A genuine 4-edge junction needs 60 of 64 columns; overflow is a
fatal `stop` on rank 0 that hangs the other ranks in the next collective.

**1.6 Pin `sheath_weak_rmax = 2.0` in the namelist**, and add a setup ERROR if it is <= 0 with
`sheath_zj_weak`. **[V]** Everything in the trace row is proportional to `zj_sat`, and the cap
`wk_cap = sheath_weak_rmax*max(|zj_sat|,1e-30)` is what makes the row degrade to `dzj = 0` as
`b_n -> 0`. With rmax = 0 the row becomes `zj = 0` at zbig over the tangential wall. A parameter
documented as startup insurance silently becomes the tangential-wall boundary condition.

**1.7 Fix stale documentation that keeps regenerating errors.** **[V]** (a) `Phi = -F0*u` in
`mod_sheath_bc.f90:26-32`, `doc/sheath_bc_whitepaper.tex`, `doc/HANDOFF.md` - the correct
convention is `Phi = +F0*u`; the JOREK wiki defines `y = -R sin(phi)` ("phi goes clockwise from
above"), so JOREK's (R,Z,phi) is right-handed and its e_phi is minus the standard one. See
`doc/jorek_solps_parallel_ohm.md`. A reviewer reproduced the old error this week purely from
these comments. (b) `phys_module.f90:252` describes `sheath_zj_weak` as "a penalty integrated
over the boundary Gauss points" - the implementation is row REPLACEMENT, and
`mod_sheath_trace.f90`'s header documents the penalty approach as tried and abandoned.

**1.8 Remove dead `wk_dfac`** (declared, zeroed, assigned twice, never read). Keep the comment
that explains why the Jacobian is deliberately left unscaled.

---

## Stage 2 — type 5. DEMOTED: type 1 carries both targets, so this is completeness, not correctness.

Only worth doing if 0.2 shows type 5 carries significant wetted area that is genuinely
field-crossing. If type 5 is grazing wall in this grid, it belongs in Stage 4 instead.

**2.1 BLOCKER first.** **[V]** Type 1 -> side 2, type 5 -> side 3, so a boundary element spanning
one node of each hits the "should never happen" branch (`construct_matrix_mod.f90:129-134`),
prints a per-element WARNING and `cycle`s - dropping the ENTIRE element: Mach1, sheath heat flux,
density flux, everything. Under 32-rank OpenMP interleaving that warning is effectively
invisible. Fix before covering 5: replace the 2-vs-3 hard `cycle` with a per-node local
direction, or at minimum count the drops and report them at the end of assembly.

**2.2 Then cover type 5 exactly as type 1:**
    bcs(5)%sheath_zj_weak = .true.
    bcs(5)%dirichlet%zj   = .false.
    bcs(5)%dirichlet%u    = .false.
    bcs(5)%natural%u      = .false.     ! REQUIRED off - hard stop otherwise
    bcs(5)%natural%zj     = .true.      ! completes the row if any trace row is skipped
    bcs(5)%dirichlet%w    = .true.
    bcs(5)%mach1          = .true.      ! preset default for 3:5

**2.3 Re-measure the target-to-target dPhi and the integrated in-out ratio.** If 0.1 showed the
sheath was on one target only, this number replaces the -23.10 V and the bias-scan slope has to
be re-derived.

---

## Stage 3 — type 9 (corners). Release, do not flag.

**[V]** `direction(2)` is never assigned when both nodes are type 9 (`side1 = side2 = 0` matches
neither branch and the mismatch test only catches 2-vs-3), and
`mod_boundary_matrix_open.f90:158` then evaluates `direction_perp(1) = 6 / direction(2)` on a
stale or uninitialised value - integer divide-by-zero or silent garbage geometry into the trace
accumulator. Note the `grid_to_wall = .false.` branch DOES handle (9,9) explicitly; the
`grid_to_wall = .true.` branch, which this run uses, does not.

**3.1** Extend the side map to cover type 9, and add an explicit branch for `side1 = side2 = 0`.
**3.2** Then `bcs(9)%dirichlet%u = .false.`, `bcs(9)%dirichlet%zj = .false.`, and **leave
`sheath_zj_weak = .false.`**. Corners inherit the projection from adjacent covered edges, which
is what the trace accumulator already does by summing both adjacent edges into one row. This
closes the j-V loop at the seam, which is half-open today.

---

## Stage 4 — CORRECTED 2026-09-01. Type 4 is a TARGET; type 2 is the tangent wall.

**The earlier version of this stage was wrong and is retracted.** It claimed type 4 was the
tangent wall and that covering it with the weak route was a no-op. Both false. The RECAP comment
in `update_boundary_types.f90` is correct:

    1: TARGET, side 2        2: TANGENT, side 3      3: CORNER, between type-1 and type-2
    4: TARGET, side 2        5: TARGET, side 3       9: CORNER, between type-4 and type-5

**[V]** `grid_xpoint_wall.f90:1637-1642` (and `:1685`) assign boundary 4 at `j .eq. n_leg-1` with
the comment `! LEFT LEG` / `! RIGHT LEG` - the divertor leg END FACES. The `12 -> 4` remap in
`update_boundary_types.f90` is gated on `use_simple_bnd_types`, which is `.false.` by default and
not in the namelist, so it never fires and cannot be used to infer what type 4 is. An earlier
campaign note already recorded the answer: type 4 carries 269 A/m^2 against type 1's 625 A/m^2 at
a consistent wall potential - real current, not a vanishing j_sat.

**The adjacency graph the grid actually builds:**

    side 2:  1 (target)  4 (target, leg end faces)
    side 3:  2 (TANGENT wall)  5 (target, side 3)
    corners: 3 between 1 and 2      9 between 4 and 5

The corner types sit exactly on the side-2/side-3 transitions, so the grid does not create direct
(1,2) or (4,5) edges - the dropped-element hazard is mediated by 3 and 9, which are unmapped and
assemble with either side. **{1,4} is one contiguous side-2 target set; {2,5} is the side-3 set.**

**4.1 THE FIRST STEP IS TYPE 1 + TYPE 4.** Same physical class, same logical side, edges assemble,
and the leg end faces carry ~43% of type 1's current density. Same configuration as type 1:
    bcs(4)%sheath_zj_weak = .true.
    bcs(4)%dirichlet%zj   = .false.
    bcs(4)%dirichlet%u    = .false.
    bcs(4)%natural%u      = .false.
    bcs(4)%natural%zj     = .true.
    bcs(4)%dirichlet%w    = .true.
    bcs(4)%mach1          = .true.

**4.2 The previous 1+4 failure was the GAUGE, not the geometry.** **[V]** The code records
"max ePhi/kTe 10.6 on one sheath type, 35.8 on two, then blow-up. Set sheath_init_u = .true. AND
sheath_init_u_all = .true." Doubling the free-boundary area doubles the harmonic excursion from
the u = 0 pinned gauge nodes. `sheath_init_u_all = .true.` is the documented fix and is being
tested independently.

**4.3 Then the corners, 3 and 9** - release `dirichlet%u` and `dirichlet%zj`, do NOT set
`sheath_zj_weak`. They already receive sheath rows through the edge-level OR while their u stays
frozen, which is a half-open j-V loop today.

**4.4 Type 2 (TANGENT wall, includes the sub-PFR dome) is the Robin candidate**, not type 4.
`natural%u = .true.` with `sheath_wall_pen`, keeping `dirichlet%zj = .true.` This is where the
dome experiment lives: `sheath_V_wall_asym` was R-antisymmetric and target-only, so it never
moved the dome potential, and the dome is the surface most directly setting the PFR E_r. Blocked
by the per-edge route selection (`:182-190` ORs over both nodes, so a mixed weak/Robin
configuration silently picks Robin at every seam).
MEASURE: `int n*v_ExB . dl` across the PFR traverse, dome-pinned vs dome-floating.

**4.5 Type 5** stays as in Stage 2 - completeness, and it needs the side-map work because
(1,5)/(4,5) would be direct side-2/side-3 edges if the grid ever produces them.

## Stage 5 — the electron branch. Independent of coverage; keeps Stage 4 alive when a target detaches.

**[V, mechanism corrected 2026-09-01]** The row loses its grip on u through
`dxlim_dx = 1/(1+exp(-z))` (`mod_sheath_bc.f90:314`) with `z = (X - X_min)/sheath_smooth_dX`,
NOT through `x_frozen` (which needs `x_lim < -30` and is unreachable while the limiter bounds
x_lim near X_min = -3). At the observed crash value `ePhi/kTe = -28.9`, `z = -51.8` and
`dxlim_dx ~ 3e-23`. **E-folding width is `sheath_smooth_dX = 0.5` in X, i.e. 0.5*Te in VOLTS** -
at Te = 5 eV the row loses its grip within 2.5 V of the clamp.

**5.1 Add `sheath_sat_slope_e`**, the exact mirror of `sheath_sat_slope`: subtract
`s_e*softplus(-(X - X_min))` so `df/dX -> -s_e /= 0` as X -> -inf and the row ALWAYS retains a u
column. Same status as the ion-side slope: a regularisation with a convergence knob - report at
the smallest s_e that runs and show the answer does not move.

**5.2 Optional: `sheath_X_min < 0` as a sentinel for `X_min = -Lambda(Ti,Te)` locally.** NOTE
this run has `sheath_Lambda_local = .false.`, so Lambda = 3.0 and `X_min = -3.0 = -Lambda`
exactly - the clamp already sits at Phi = 0, electron saturation, correctly placed. This item
only matters if `Lambda_local` is enabled.

**5.3 Better watchdog.** `weak` is the WRONG indicator for this failure: it stayed flat at 8e-3
through the entire crash because the zj row WAS being satisfied - the divergence was in the u
equation. Report instead the sheath conductance relative to the polarisation diagonal,
`|dzj_du*B.n*dl*theta*tstep| / (rho*R^3)`, which `sheath_stiff_max` already computes
(`mod_boundary_matrix_open.f90:~454-461`). Its collapse is the leading indicator.

---

## Stage 6 — the skip criterion. Needed before coverage is trustworthy at scale.

**[V]** The current test is `st_D < 1e-14*d_max` with `d_max` RANK-LOCAL
(`mod_sheath_trace.f90:205-211`), so the same geometric DOF is skipped on a rank owning a strike
point and kept on a rank owning only tangential wall. It also treats value and derivative rows
completely differently: measured `D min 6.34e-14`, `max 1.62e-3`.

**6.1** Accumulate `D_a^raw` - the same integral as `wk_D` with `wk_wgt` forced to 1 - and use
`A_a = D_a/D_a^raw` in [0,1]: dimensionless, h-independent, identical on every rank, no MPI
reduction, same treatment for value and derivative rows. Skip when `A_a < sheath_weak_amin`.
**6.2** Hysteresis: `A_a` is state-dependent, so a single threshold can limit-cycle (rows
flipping replaced/assembled every Newton step). Two thresholds plus per-DOF state that
deliberately does NOT reset per matrix construction. Ship with the upper threshold at 0
(hysteresis off) until measured.
**6.3** Store the boundary type per trace row from `nodes(i)%boundary`, not the element's
`bnd_type1` - exact at a junction between two covered types.
**6.4** Observability, all currently missing: per-type row counts (accumulated / replaced /
skipped), an `A_a` histogram with a bin for the hysteresis band, a count of rows that flipped
activity this step, and an explicit WARNING when `ePhi/kTe min` goes negative with a streak
counter - it preceded the last crash by several outputs and nothing calls it out.

---

## Stage 7 — longer term: one physical floor instead of two tuning knobs

**[I]** `j_sat = e*n*(c_s*g(b_n)*|b_n| + v_perp)` with `v_perp ~ D_perp/lambda_n`. With
D_perp ~ 1 m^2/s and lambda_n ~ 1 cm, `v_perp ~ 100 m/s` against c_s ~ 2.2e4 m/s at 5 eV, i.e.
`v_perp/c_s ~ 4.5e-3` - within 10% of the empirically chosen `sheath_min_bn = 0.005`. The gate
has been reproducing the cross-field flux floor by hand. Making it explicit removes
`sheath_min_bn` AND `sheath_wall_pen`, makes the characteristic well-posed at exact tangency with
no gate, and lets types 1/2/4/5 share one uniform statement - which is what the preset comment
has been asking for.

NOTE **[V]**: `sheath_min_bn` currently does NOTHING inside `sheath_current` - it is imported and
never used (`g_eff = g_bn`). On the weak route it survives only in the diagnostic weight. The
setup warning claiming it puts a floor on the obliqueness factor is stale and misleading.

---

## What NOT to do

- **Do not port the obliqueness gate into the weak route.** **[V]** `wk_wgt` already multiplies
  `wk_D`, `wk_F`, `wk_S` and `tr_J` identically, and the row is written as `F/D`, so any weight
  uniform over a row's support cancels exactly. Gating a REPLACED row is not a well-posed
  operation - you can only choose whether to write it.
- **Do not port `sh_pen_c` to a weak type.** **[V]** Blocked by a hard stop, and its motivation
  (gating the sheath leaves u free) is specific to the Robin route, where the surface term IS u's
  boundary condition. On the weak route the sheath never touches the u row.
- **Do not use the weak route on type 4 or 2.** See Stage 4.
- **Do not tune `sheath_sat_slope` against the electron branch.** The softplus is exponentially
  small for X < 0; it cannot affect that failure in either direction.

---

## Monitoring, every run

    sheath trace rows: N accumulated, M replaced, K below the floor   ! M<N or K>0 = degrading
    max|j/jsat|                        ! <20 fine; ~40 gate at 6%; ~157 rows start crossing the floor
    ePhi/kTe min                       ! going NEGATIVE is the earliest precursor
    per-type I_sheath vs I_Ampere      ! divergence isolated to one type = that type misconfigured
    D min/max                          ! collapse = gate closing
    grep -c "boundary element incoherent|not included"   ! silent element drops
