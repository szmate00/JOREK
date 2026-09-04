# SUPERSEDED - see doc/HANDOFF_sheath_multitype.md

> This document describes branch `bc-tests` and is kept for history only. Several of its
> conclusions have since been falsified by measurement - in particular it reports the 1+4
> configuration as failed, which was true only before `sheath_weak_detmin` existed. Do not use its
> namelist recommendations. The authoritative account is
> **`doc/HANDOFF_sheath_multitype.md`, section 7**.

# Handoff: JOREK sheath boundary condition + HFSHD investigation

Repo: JOREK MHD. Branch `bc-tests`, commit `8f6534ad8`. model600, ASDEX Upgrade shot 38773,
`n_tor = 1` (axisymmetric). Cluster runs; laptop has no MPI toolchain, so nothing here is
compiled locally.

**Question:** does a self-consistently evolved sheath potential at the divertor targets produce
the in-out density asymmetry (HFSHD, high-field-side high-density)?

---

## Part 1 — the boundary condition: DONE, working

A weak Galerkin row replacement imposing the Langmuir j-V characteristic on the boundary trace
space. Replaces the previous Dirichlet `u = 0` at the wall, which forbade the physics
($\Phi$ = const means $E_r$ = 0 by construction).

    closure  I_sheath == I_Ampere to 4 significant figures, both targets   (0.00 %)
    weak residual  6.5e-4      Gauss residual  9.0e-3
    3900 steps to wall-clock timeout, no crash, NO tuning parameter

Files: `models/model600/mod_sheath_bc.f90` (all sheath physics, stateless),
`mod_sheath_trace.f90` (row accumulator), `mod_boundary_matrix_open.f90` (assembly),
`mod_boundary_conditions.f90` (apply), `mod_sheath_diag.f90` (diagnostics).

Full writeup with a from-scratch finite-element primer: `doc/sheath_bc_whitepaper.pdf` (19 pp).

Namelist:

    bcs(1)%sheath_zj_weak = .true.
    bcs(1)%dirichlet%zj   = .false.   ! REQUIRED - a pinned zj opens the j-V loop
    bcs(1)%natural%zj     = .false.   ! NOT required on the weak route, see below
    bcs(1)%dirichlet%u    = .false.
    bcs(1)%dirichlet%w    = .true.
    bcs(1)%mach1          = .true.
    bc_natural_open       = .t.

Keep `dirichlet%u = .true.` on at least one other boundary type as a gauge, or the constant
mode of `u` is undetermined.

`natural%zj` is required on the OLD natural%u route, where nothing replaces the boundary zj rows
and releasing their Dirichlet leaves an incomplete weak form. On the WEAK route those rows are
replaced by the characteristic itself, so the surface term would only feed rows that are then
overwritten, and the j-V loop closes through the constraint: zj at the wall is not pinned, it is
set to zj_sh(u,rho,Ti,Te), and the interior zj = Delta*psi equation makes psi's normal derivative
follow. Measured: left off, 3900 steps at 0.00 % closure. Not required here.

**This result stands regardless of how the physics question resolves.**

---

## Part 2 — the physics: NEGATIVE so far

`sheath_V_wall_asym` imposes an antisymmetric target bias, so a known potential difference can
be applied and the density response measured. This separates "does the transport chain work"
from "is the self-consistent potential large enough", which no self-consistent run can do.

At asym = 40 against an identical unbiased run from the same branch point, step 700:

    dPhi at strike points      -23.10 V  ->  -59.30 V    (-36.2 V)
    n_in/n_out at strike pt      0.531   ->    0.866     (+63 %)   POINTWISE - MISLEADING
    int n ds ratio, common w     1.171   ->    1.189     (+1.5 %)  EXTENSIVE - the real one

**The pointwise signal is a rearrangement, not a transfer.** On the inner leg the strike-point
value rose 22 % while the peak FELL 12 % and moved 1.2 cm toward the strike point; both legs
gained ~16 % together. The probe sits exactly where plasma is being pushed to.

Sensitivity on the extensive metric: **0.0005 per volt**. A factor-2 asymmetry from 1.17 would
need ~1660 V against the -23 V the plasma generates by itself, i.e. an amplification of ~70x.

Consistent with two independent estimates: the direct drift contribution
$\Gamma_{ExB}/\Gamma_\parallel \approx 12\,\%$, and the ExB direction being unresolvable on
the inner leg (coherence 0.023) because there is no coherent inter-leg flow to find.

**Interpretation:** in a sheath-limited divertor the sheath ExB drift reshapes profiles but
does not transfer plasma across the private flux region. It is NOT falsified as a mechanism,
because the amplifier has never been present - every run so far has been either hot and
tenuous (drive but no recycling) or cold and detached (recycling but no drive).

---

## Part 3 — what to do next

### Fix first, in EVERY run (config only, no rebuild)

    T_min_sheath = 1.d-6      ! unset -> falls back to T_min_neg, knee 1.47 eV, which
                              ! compresses a 19 % target Te contrast to 2 %

All four wall polygons are traversed in bowtie order, so `inside_polygon` ray-casting returns
two triangles rather than a rectangle. Measured coverage 50.5 / 30.1 / 50.3 / 50.5 % of the
intended boxes - both pumps at half area and the wrong shape.
`mod_particle_wall_interaction.f90:83` warns about ordering explicitly. Swap the last two Z
values in each:

    valves(1)%poly_Z = -1.11, -1.05, -1.05, -1.11
    valves(2)%poly_Z = -1.11, -1.05, -1.09, -1.15
    P1     %poly_Z = -1.14, -1.12, -1.12, -1.14
    P2     %poly_Z = -1.2,  -1.151, -1.151, -1.2

### The main experiment: a 2D scan

Compute is not a constraint. Run the full grid rather than sequentially:

    puff x {1, 2, 3, 5}  X  sheath_V_wall_asym {0, 20, 40}

(1/0 and 1/40 already exist.) This gives d(ratio)/d(dPhi) **at each recycling level**, which is
the actual question. Plot the slope against target density: rising steeply means recycling
amplifies the drift and HFSHD is reachable; flat across the scan means the mechanism is too
weak at any accessible drive and the answer lies in the missing physics below.

Target for the puff levels: `T_e` 5-15 eV and `n` 1e19-1e20 at the targets - the
high-recycling regime, which has never been sampled. **Each point must reach its own
equilibrium**; the 10x puff case was still visibly evolving at step 2000.

### In parallel, zero compute, on existing checkpoints

1. **Measure the PFR flux directly.** Everything so far infers the transfer from target
   densities. `util/postproc_pfr_line.in` + `util/analyse_pfr_line.py` project the ExB
   velocity onto a traverse across the PFR, so `n*v_l` is the transfer itself. Run on the
   baseline and asym=40 checkpoints. Confirms or overturns the negative by direct measurement.

2. **Fix the volume metric.** The common integration window is only +-5.5 cm, set by the outer
   target's reach, and it EXCLUDES the inner leg's density peak at +5.6 cm. Target traces
   cannot properly measure divertor content. Use `rectangle Rmin Rmax nR Zmin Zmax nZ 0.`
   over the divertor and integrate n per leg. Reader script not yet written (HDF5 output).

---

## Gotchas that WILL bite

1. **The postproc `Phi` expression has the wrong sign for model600.** It returns `+F0*u`
   (`mod_expression.f90:1757`, flagged `!### sign?` in the source) but model600 implements
   `v_pol = +R grad(u) x e_phi`, so `Phi = -F0*u`. **Negate it.** Unflipped, the ExB
   direction inverts and so does the conclusion. `V_ExB_R/Z` are consistent and safe as-is.
2. **`Te` is not the electron temperature.** It is `T0*fact_T/2` with `T0 = Ti + Te`, i.e.
   the mean. Use `T_e` and `T_i` (`mod_expression.f90:1999-2003`).
3. **Put `set units 0` in every postproc script.** The default is 0 = SI, but exports have
   come out in JOREK units (n ~ 1e-2, T ~ 1e-3); `analyse_targets.py` now refuses those.
4. **psi_N does NOT separate SOL from PFR** - it exceeds 1 on BOTH sides of a separatrix leg.
   Use signed arc length from the strike point. The strike point is `min(psi_N)`, which is
   geometric; the `j_sat` peak wanders when the target detaches.
5. **Never use the strike-point density ratio as the HFSHD metric.** Use the line-integrated
   one, over a COMMON window on both targets. Unequal windows turn a 37 % length difference
   into a density ratio - that is how a baseline with a pointwise ratio of 0.53 once reported
   an integrated ratio of 2.05.
6. **`I_sheath` and `I_Ampere` share the `zj*(B.n)/(F0*mu0)` conversion**, so a
   normalisation error cancels in their ratio. They verify the solver, not `sheath_current`.
   A factor of two in `a_n` passes every check the code prints. Real verification needs a
   circular-limiter case with uniform Te, where `ePhi/kTe` must equal Lambda exactly.
7. **The ExB direction cannot be measured from a boundary trace.** Coherence
   `|mean(v)|/mean(|v|)` is 0.02 on the inner leg - the direction rotates within the window.
   The transfer happens in the PFR volume. Use the pol_line traverse.

---

## Falsified - do not retry

- **`sheath_weak_relax` < 1** - worse on every metric at timeout. The period-2 alternation it
  targeted was a transient that cleared itself.
- **`sheath_X_min` below -3** - actively harmful. X_min is the CLAMP on the electron branch,
  so -6 raises the floor on f from -19 to -445. Made a failing case fail faster.
- **Sheath on boundary types 1 + 4 together** - four attempts, including a genuine bug fix
  (rows were being written for nodes whose `u` is Dirichlet-frozen) that changed nothing.
  Two sheath surfaces close a circuit through `u` with no damping in the loop. Enabling
  types 5 and 9 as well to remove the seam is worse: type 5 is ~50 m^2 of grazing wall.
- **`sheath_wall_vel`, `sheath_init_u_all`, a beta penalty at any value, the `t_min`
  resistivity clamp** - all falsified during the boundary condition work.

---

## Tools and documents

    util/analyse_targets.py     util/postproc_targets.in    target profiles, both traps handled
    util/analyse_pfr_line.py    util/postproc_pfr_line.in   PFR traverse, ExB direction
    doc/sheath_bc_whitepaper.pdf                            the BC, with an FEM primer
    doc/hfshd_test_plan.html                                the test plan, ordered
    doc/solps_vs_jorek_temperature_equations.md             what model600 is missing vs SOLPS

**Known gaps in model600**, independent of the boundary condition but relevant here: no
$\nabla B$ terms in the energy equations, and no thermoelectric heat-flux divergence
$\nabla\cdot(0.71 T_e j_\parallel/e)$. Both are in SOLPS and both act in the in-out channel.
The 0.71 grad-parallel-Te thermal force in Ohm's law exists on branch `u-bcs-tests-mate`
(`thermoelectric_ohm`) but is NOT on `bc-tests` - cherry-pick if needed.
