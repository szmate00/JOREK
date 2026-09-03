# Handoff: the weak-form sheath j-V boundary condition, and why multiple node types are hard

Repo: JOREK. Branch **`sheath-jsat-vpar-38ab278`**, head `ab2c8caa8`. `model600`, ASDEX Upgrade
shot 38773, `n_tor = 1` (axisymmetric), `n_order = 3`. Cluster runs; the laptop has no MPI
toolchain, so nothing here is compiled locally beyond single-module syntax checks.

**The question this exists to answer:** does a self-consistently evolved sheath potential at the
divertor targets produce the in-out density asymmetry (HFSHD, high-field-side high-density)?

This document supersedes `doc/HANDOFF.md` (branch `bc-tests`, older). Companion documents:
`doc/sheath_bc_whitepaper.pdf` (19 pp, from-scratch derivation plus a finite-element primer) and
`doc/sheath_multitype_ideas.md` (the multi-type investigation in full, with every falsified
hypothesis recorded).

---

# 1. The physical model

`model600` is JOREK's **reduced-MHD** tokamak model. "Reduced" means the magnetic field and the
velocity are represented by scalar potentials rather than full vectors, which removes the fast
magnetosonic wave and lets you take timesteps set by the physics rather than by the fastest wave
in the system.

## 1.1 The ansatz

    B    = F0/R e_phi  +  1/R grad(psi) x e_phi          (+ toroidal field, poloidal flux)
    v_pol = +R grad(u) x e_phi                            (ExB stream function)
    v     = v_pol + v_par B

`psi` is the poloidal flux, `u` the electric potential stream function, `v_par` the parallel
velocity. `F0 = R*B_phi` is a constant (the toroidal field function).

**Sign convention, and it has caused real errors here.** Hoelzl et al. 2021 eq. 26 defines
`u = Phi/F0` with `v_pol = -R grad(u) x e_phi`. `model600` implements `v_pol = +R grad(u) x e_phi`.
The code's `u` is therefore *minus* the paper's:

    Phi = -F0 * u

and the same inversion applies to `psi` and to the current variable: the code's `zj = Delta* psi`
is minus the paper's `j`. Memory note `phi sign convention` records that a postprocessing script
was "corrected" in the wrong direction on the basis of the paper before this was pinned down.

## 1.2 The variables

This build has `with_TiTe = .true.` and `with_vpar = .true.`, so **`n_var = 8`**:

| index | variable | meaning |
|---|---|---|
| 1 | `psi`  | poloidal flux |
| 2 | `u`    | electric potential stream function, `Phi = -F0 u` |
| 3 | `zj`   | toroidal current, defined as `Delta* psi` |
| 4 | `w`    | vorticity, defined as `Delta* u` |
| 5 | `rho`  | density |
| 6 | `Ti`   | ion temperature |
| 7 | `Vpar` | parallel velocity |
| 8 | `Te`   | electron temperature |

Note the repo default in `models/model600/mod_model_settings.f90` is `with_TiTe = .false.`;
`util/config.sh` overrides it at build time. Check `n_var` in the HDF5 restart if in doubt.

**`zj` and `w` are auxiliary variables**, not independent physics. They exist so the equations stay
second order: instead of a fourth-order operator on `psi`, you solve two second-order ones. Their
"equations" have exactly one term each (`mod_model_settings.f90:118-123`):

    zj_Eq__DeltaStar_Psi     zj  - Delta* psi = 0
    w_Eq__DeltaStar_u        w   - Delta* u   = 0

**This matters enormously for the boundary condition.** The sheath is a statement about the
*current* at the wall, so it is imposed on `zj`, and `zj`'s equation is a definition that has been
integrated by parts.

## 1.3 The equations, by their term names in the code

The RHS diagnostic names every term (`mod_model_settings.f90:78-130`), which is the most reliable
inventory of what is actually in this model:

* **`psi`** (7 terms): `eta_J`, `B.grad_u`, `eta_num_term`, `diamag_term`, `zeta_timevol_term`,
  `RE_coupling`, `thermal_force`. This is Ohm's law: `dpsi/dt = eta(j - j0) + B.grad(u) + ...`
* **`u`** (13 terms): `rho_v.grad_v`, `JxB`, `visco_term`, `grad_p`, `visco_num_term`,
  `tg_num_term`, `diamag_term`, `diamag_visco`, ... This is the vorticity equation.
* **`rho`** (14), **`Ti`** (18), **`Te`** (20), **`vpar`** (13): continuity, two temperature
  equations, parallel momentum.

**One measured fact about the `u` equation you should know before touching this BC.** `JxB` and
`grad_p` are each `2.85e-5` and agree to **0.03 %**, cancelling down to the `~1e-7` scale where the
diamagnetic term (`1.2e-7`), `rho_v.grad_v` (`7.3e-8`) and viscosity (`5e-7`) live. The sheath BC
steers `u` as a perturbation on a cancellation of two terms 200x larger. That is a structural
reason this boundary condition is numerically unforgiving, and it is **not** a pathology: `grad_p`
is `div j_dia` and `JxB` is the `div_par j_par` closure, so "they cancel" is the statement that the
diamagnetic current closes as parallel current to the plates. Which is exactly the circuit the
sheath is supposed to terminate.

---

# 2. Discretisation

## 2.1 Poloidal plane: 2D cubic Bezier finite elements

Quadrilateral elements, **bicubic Bezier (Bernstein) basis**, `n_order = 3`. Per node per
variable there are `n_degrees = ((n_order+1)/2)^2 = 4` degrees of freedom:

    DOF 1   the value
    DOF 2   the first derivative along the node's frame direction x(1,2,:)
    DOF 3   the first derivative along the node's frame direction x(1,3,:)
    DOF 4   the mixed second derivative

The basis is **C1 continuous** across element boundaries. `elements/basis_at_gaussian.f90`,
`elements/bezier_1d.f90`, `elements/mod_node_indices.f90`.

**The node "frame" is the key object for everything in section 6.** Each node carries two unit
vectors, `x(1,2,:)` and `x(1,3,:)`, which are the *physical directions* the two first-derivative
DOFs measure along. They are set by the grid generator and are **not required to be orthogonal or
even independent**. The element geometry is mapped with the same basis: `element%size(vertex, 2)`
and `(vertex, 3)` scale the frame vectors to reproduce the element edges.

## 2.2 Toroidal: Fourier

`n_tor` harmonics; this campaign runs `n_tor = 1`, i.e. axisymmetric. `n_tor_local > 1` is refused
at setup with `sheath_zj_weak`.

## 2.3 Time: linearised implicit, ONE solve per step

`core/mod_jorek_timestepping.f90:377` calls `construct_matrix` once and `:390` calls
`solve_sparse_system` once, **with no loop**. There is *no Newton iteration*.

Any comment in this codebase justifying an omitted Jacobian term because "the geometry is frozen
within a Newton iteration" is therefore describing a **systematic per-step error**, not an
iteration-count cost. This has been re-derived by reviewers several times; it is recorded in
memory as `sheath weak route code facts`.

---

# 3. How boundary conditions are applied

## 3.1 Boundary types

Every node carries `node%boundary`, an integer type. `bcs(t)` is a per-type namelist structure with
`%dirichlet%<var>`, `%natural%<var>`, `%mach1`, `%sheath_zj`, `%sheath_zj_weak`, ... The defaults
(`models/preset_parameters.f90:426-427`) are

    bcs(:)%dirichlet%u  = .true.
    bcs(:)%dirichlet%zj = .true.

**Remember that pair.** Section 6.6 is about a bug that came from forgetting the second one is a
default that a sheath-enabled type overrides.

## 3.2 Dirichlet rows

`models/model600/mod_boundary_conditions.f90`. Loops over boundary nodes; for each, the element
side determines `iv_dir` (the DOF tangential to the boundary) and `iv_perp_dir` (the normal one),
`:251-260`. The Dirichlet then pins **the value and the TANGENTIAL derivative only** (`:400-405`,
"Fix derivatives in one direction"), writing the row with a large number `zbig` that annihilates
everything else in it.

**The normal-derivative DOF is deliberately left free.** For `psi` that free DOF is how the wall
current responds to the sheath (`mod_sheath_bc.f90` header). Section 6 returns to this.

## 3.3 Natural (weak / surface-integral) rows

`models/model600/mod_boundary_matrix_open.f90`, reached when `bc_natural_open = .true.`. When a
variable's equation is integrated by parts, a surface term appears. `bcs%natural%<var> = .true.`
assembles it; `.false.` discards it.

Discarding is harmless **only while a Dirichlet freezes that variable's trace** - the code says so
explicitly. Release the Dirichlet without adding the surface term and the boundary equation is the
volume weak form *missing its surface term*, i.e. **incomplete**. That is the single most common
failure mode in this whole subsystem and it produces an unbounded runaway, not a small error.

## 3.4 Which side is which

`matrix/construct_matrix_mod.f90:135-168` maps the two nodes' boundary types to `direction(2)`,
the DOF running *along* the boundary edge; `mod_boundary_matrix_open.f90:178` then sets
`direction_perp(1) = 6/direction(2)`.

    types 1, 4, 9   ->  direction(2) = 2,  normal-derivative DOF = 3
    types 2, 5      ->  direction(2) = 3,  normal-derivative DOF = 2

Note also: the `grid_to_wall` branch at `:111` is guarded by `n_wall_blocks > 0`, and
`n_wall_blocks` defaults to **0**, so the **ELSE branch at `:135-168` is what runs here**. Several
reviews have analysed the wrong branch.

---

# 4. `keep_current_prof` and the fake HFSHD

Read this before believing any in-out asymmetry result from this setup.

`keep_current_prof` (`phys_module.f90:969`) adds an artificial source `eta*(j - j0)` to hold the
initial current profile. `j0` is rebuilt from the **initial** FF'/n/T profiles evaluated at the
**evolving** psi. In the SOL the profile tanh tails give near-core values, and in the private flux
region `psi_N < 1`, so the source is large and spurious on open field lines. Its stiffness goes as
`eta_T ~ T^-3/2`, which is enormous in a cold leg.

Chain, confirmed causally: spurious SOL/PFR current -> `{psi,j}` torque -> localised X-point
convection cell -> PFR `perp_convection` dipole -> an in-out density asymmetry that **scales with
eta** and looks exactly like HFSHD.

`keep_current_prof_confined` (`:980`) masks the source to the confined region with smooth tanh
masks in `psi_N` and in `Z` beyond the X-point.

**What we saw (bracket completed 2026-07-23):** the asymmetry dies under all three of
(a) `keep_current_prof = .false.`, (b) `keep_current_prof_confined` with the mask on, and
(c) mask off but low `eta`. Only *SOL-active source + high eta* produces it. **The observed
asymmetry was 100 % artifact.**

Two consequences that shaped everything afterwards, from the 50-term RHS comparison
(`rhs_comparison-2.pdf`, term changes normalised against the dominant term in each equation):

1. `psi eta_J` changed by **45 %** of its equation's dominant term - the artifact, localised. But
   `Te` and `Ti` changed by **< 0.15 %**. So the fake HFSHD was a full density asymmetry with
   **zero temperature asymmetry** - purely current-driven ExB. The sheath route works through
   `Phi = Lambda kTe/e` and therefore *needs* a `Te` asymmetry. **The fake HFSHD is not a template
   for the real one, and its disappearance is not evidence against the sheath route.**
2. `w = Delta* u` is a pure wall-hugging boundary layer in both runs, interior `~0`. So `u` is
   nearly **harmonic** inside, and the PFR ExB pattern is set almost entirely by `u`'s BOUNDARY
   VALUES. That is precisely the handle the sheath BC has - and it also explains why the BC is so
   unforgiving: no bulk term damps boundary noise, so grid-scale noise on the wall propagates
   harmonically straight into the private flux region.

---

# 5. The weak-form sheath j-V boundary condition

## 5.1 The physics

Langmuir probe characteristic (Stangeby eq. 2.68), referenced to the wall:

    j_par = j_sat (1 - exp(-X)) ,   X = e*Phi/(k Te) - Lambda ,   Phi = V_sheath - V_wall

`X = 0` is the floating potential `Phi_f = Lambda kTe/e` (zero net current); `X -> +inf` is ion
saturation; `X = -Lambda` (grounded wall) is electron saturation, drawing
`j_sat(1 - e^Lambda) ~ -19 j_sat`. The branch asymmetry is the whole difficulty.

It is evaluated **forward** - current as a function of potential. That direction is single-valued
and `dj/dPhi ~ exp(-X)` tends smoothly to zero at ion saturation. Inverting for `Phi` is singular
exactly at `j -> j_sat`, which is where divertor conditions actually sit. (The older nodal route
inverts; that is why it is the older route.)

Two regularisations, both in `mod_sheath_bc.f90`:

    f(X) = (1 - exp(-X)) + s ln(1 + exp(X))          s = sheath_sat_slope
    X_lim = X_min + dX ln(1 + exp((X-X_min)/dX))     electron-side limiter

The softplus is exponentially small for `X << 0`, so the electron branch and the floating potential
are untouched while `f` becomes unbounded above and every demanded current is reachable at finite
`X`. **Trap:** lowering `X_min` to "give the solver room" is backwards - `X_min` is what *stops*
the electron current growing, and below the clamp `dX_lim/dX -> 0`, so `d(zj_sh)/du -> 0` and the
sheath stops responding to the potential at all.

In code units (`mod_sheath_bc.f90` header):

    e*Phi/(k Te) = (a_n u / 2 - v_w)/Te      a_n = -2 e F0 sqrt(mu0 rho0)/m_i
    zj_sat = c_sat rho g(b_n) c_s / |B|      c_sat = -a_n/2,  c_s = sqrt(gamma (Ti+Te))

`g(b_n)` is the Chodura-Riemann function *including its sign*, the same one `bcs%mach1` imposes on
`Vpar` at the same wall - so the two are consistent by construction.

## 5.2 Why WEAK rather than nodal

The nodal route (`bcs%sheath_zj`) imposes `zj = zj_sh` pointwise at boundary nodes. Measured, it
does that essentially exactly - `|zj - zj_sh|/|zj_sat|` reaches `1.8e-4` at the nodes - **and the
integrated currents still differ by 33 % on the inner target**, because the cubic trace *between*
the nodes is uncontrolled. The weak residual

    F_a = integral over Gamma of  N_a (zj_sh - zj) dS

stays at O(1) while the nodal residual falls four decades. So a projection onto the trace space has
real work to do that the nodal constraint does not do.

## 5.3 Why REPLACE the row rather than penalise it

`zj = Delta* psi` is integrated by parts and its surface term is refused (its Jacobian needs
normal-derivative columns the trial loop cannot produce). Harmless while `dirichlet%zj` freezes the
trace - but the weak route must release it, so the equation is then incomplete. **Adding a penalty
to a wrong equation cannot fix it:** measured, `beta = 1e-2` left the incomplete equation dominant
and `beta = 1` drove a period-2 divergence. Writing the row with `zbig` annihilates it instead,
exactly as every Dirichlet in the code already does.

## 5.4 Why row normalisation is mandatory

A Galerkin trace block has internal scale structure a pointwise Dirichlet does not: its mass matrix
goes as `h` (value x value), `h^2` (value x derivative) and `h^3` (derivative x derivative), so at
`h ~ 1 mm` the block spans `~1e6` internally. A uniform coefficient therefore cannot work at any
magnitude - `beta = 1e9` died in one step. Dividing each row by its own diagonal
`D_a = int N_a N_a dS` makes every row O(1); a single `zbig` then makes them uniformly dominant and
no penalty parameter is needed.

**Corollary that has misled three parameter scans:** a weight `wk_wgt` that multiplies `D_a`, `F_a`
and every `J_ab` identically **cancels exactly** out of a row written as `J/D` and `F/D`. It cannot
fade a replaced row - it only shrinks `D_a`. Only a weight that VARIES within a row does anything.
Row replacement is binary, so any "fade" has to be the decision whether to write the row at all.

## 5.5 The files

    mod_sheath_bc.f90        all sheath physics, stateless; also the node-frame helpers
    mod_sheath_trace.f90     the row accumulator and the row writer
    mod_boundary_matrix_open.f90   assembly of D_a, F_a, S_a and the Jacobian block
    mod_boundary_conditions.f90    Dirichlet rows, and the frame freeze
    mod_sheath_diag.f90      I_sheath / I_Ampere / ePhi(kTe) / per-type diagnostics

## 5.6 The working namelist

    bcs(1)%sheath_zj_weak = .true.
    bcs(1)%dirichlet%zj   = .false.   ! REQUIRED - a pinned zj OPENS the j-V loop
    bcs(1)%natural%zj     = .false.   ! NOT required on the weak route
    bcs(1)%dirichlet%u    = .false.
    bcs(1)%natural%u      = .false.
    bcs(1)%dirichlet%w    = .true.
    bcs(1)%mach1          = .true.
    bc_natural_open       = .true.
    sheath_Lambda = 3.0 / sheath_sat_slope = 0.03 / sheath_zj_ratio_max = 20
    sheath_weak_rmax = 2.0 / sheath_weak_relax = 1.0
    sheath_init_u = .false.          ! see below
    sheath_weak_wmin = 0.0           ! see below

Keep `dirichlet%u = .true.` on at least one other boundary type as a gauge.

**`dirichlet%zj = .false.` is not optional.** Measured with it `.true.`: `I_Ampere` stayed constant
to four digits for a whole run while `I_wall` ran from -15540 A to -23 A. The sheath adapted all the
way to floating, the plasma current could not move at all, and `u` diverged trying to reconcile them.

**Two flags that are actively harmful and are documented as such:**

* `sheath_init_u = .true.` - setting `u` to the floating potential at startup. 104 steps -> 4 on one
  config, 219 -> 8 on another. **Never turn it on.** (It was nevertheless recommended throughout an
  earlier campaign, so several "falsified" hypotheses were falsified against a poisoned baseline.)
* `sheath_weak_wmin > 0` - gating the row on the *solution-dependent* validity weight. Clean A/B on
  1+5: `0` runs 305 steps, `0.5` dies at **2**, because at step 2 it skipped 706 of 920 rows.
  Skipping does not soften a row, it deletes that DOF's equation. **Leave it at 0.** This is the
  cautionary tale that section 6.5 is designed around.

---

# 6. Multiple boundary types: what is hard and why

## 6.1 The result that needs explaining

| config | outcome |
|---|---|
| type 1 alone | **~3900 steps, TIMEOUT**, converged, weak residual 6.5e-4 and still falling, `I_sheath = I_Ampere` to 4 s.f. on both targets |
| type 5 alone | 308, crash |
| type 1 + 5 | 305, crash (weak 3.3e-3, 4-digit closure) |
| type 4 alone | 8, crash |
| type 9 alone | 4, crash |
| type 1 + 4 | 8, crash |

Type 1 alone reaches a converged steady state and **stops being a stability problem**. Every
multi-type configuration dies *while still evolving*. That is a difference of kind, not degree, and
it is what defeated every namelist-level hypothesis - all of them explain a factor in step count.

## 6.2 The grid: where the boundary types come from

`grids/grid_xpoint_wall.f90`, reached from `core/mod_flux_grid.f90:56-63` with
`grid_to_wall = .true.`, `n_wall_blocks = 0`, `xcase < UPPER_XPOINT`. It carries the comment
`!!rks only for ITER wall for the moment`. **The grid is NOT rebuilt on restart**
(`jorek2_main.f90:374` opens `if (.not. restart)`), so every boundary type, node position and frame
comes from the restart file and is constant across every run being compared.

**Types 1, 2, 3 are a bitfield generated inside the FLUX-ALIGNED grid** (`:1116-1122`, `:1174-1180`):

    boundary = 0
    if (k == 1 .or. k == n_open+n_private+1) boundary = boundary + 2   ! radial extreme
    if (j == 1)                              boundary = boundary + 1   ! poloidal end = TARGET

so 1 = target, 2 = tangent wall, 3 = both = corner.

**Types 4, 5, 9 are assigned literally and only inside `if (extend)`** - the ray-cast extension out
to the wall polygon (`:1258, :1322, :1326, :1385, :1389, :1633-1642, :1682-1692`).

**Two different mesh generators.** No namelist symmetry makes them the same kind of object.

## 6.3 The node frames differ, and that is the whole story

Flux-grid node (`:1113-1114`):

    x(1,2,:) = (dR_dt, dZ_dt)/|..|             grid-line tangent, from the CUB1D spline
    x(1,3,:) = (-PSI_Z, +PSI_R)/|grad psi|     EXACTLY the poloidal field direction

Extension node (`:1251-1252`, `:1310-1312`):

    x(1,2,:) = the ray-cast direction to the wall
    x(1,3,:) = (cos(tht_ext), sin(tht_ext))    an interpolated GEOMETRIC angle
               (in the LEG block, a secant through the neighbouring nodes)

Two consequences, both measured on the real grid with `util/check_boundary_frames.py` (see 6.7).

**(a) The free psi DOF is analytically zero on type 1 and only there.** `:1886-1889` sets
`values(1,k,psi) = grad(psi) . x(1,k,:)`, and on a flux-grid node `x(1,3,:)` is the flux-surface
tangent, so `values(1,3,psi) = 0` **exactly**. Type 1's *normal-derivative* DOF is DOF 3
(section 3.4), and `dirichlet%psi` leaves the normal derivative free - so on type 1 the DOF through
which the wall current responds to the sheath starts at exactly zero, with no background to hide in.
On types 4/5/9 it carries `O(|grad psi|)` and the sheath's response is a small relative perturbation
on a large number. **Measured `|v3|/|v|`: 0.0000 on types 1, 2, 3; nonzero on 4, 5, 9** - confirmed
to four decimals, and in the *evolved* state, not just at construction.

This is real, it explains why `sheath_psi_jacobian` never discriminated (on type 1 it adds a column
on a DOF that is zero), **and it is NOT the discriminator** - type 2's free-DOF fraction is 1.0 and
type 2 is healthy, type 9's is 0.69 and it dies at 4.

**(b) The frame DETERMINANT is the discriminator.** `det = |x2 x x3| = |sin(angle between them)|`,
so the nodal derivative basis is conditioned as `1/det` - and `zj = Delta* psi`, the quantity the
weak row REPLACES, is built from second derivatives.

    type   nodes   det min   det mean   max 1/det   steps survived
       1      98    0.7152     0.8815         1.4   3900 (timeout)
       2      65    0.9383     0.9674         1.1   - (never sheathed)
       3       2    0.9344     0.9502         1.1   -
       5     362    0.0810     0.8745        12.3    308
       4      88    0.0187     0.3782        53.5      8
       9       6    0.0092     0.0511       108.7      4

`det mean` orders the four sheath-capable types **exactly** as their survival, across three orders
of magnitude in step count, and splits them cleanly: 0.88 / 0.87 (work) against 0.38 / 0.05 (fail
alone).

## 6.4 What `det -> 0` actually means

Do not read it as "ill-conditioned basis". Combine it with the other measurement: the frame-chord
alignment `|frame.chord|/|chord|` is `~1.0000` on *every* type, so the along-edge frame vector is
parallel to the boundary everywhere. `det -> 0` therefore means the **other** one - the
`direction_perp` DOF, the NORMAL-derivative DOF - is *also* nearly parallel to the boundary. So:

> **At these nodes the element has no degree of freedom pointing off the wall.**

Both derivative DOFs lie in the surface. And `dirichlet%psi` deliberately leaves the
normal-derivative DOF free precisely because that is how the wall current responds to the sheath.
The sheath is being asked to drive a DOF that does not point where it needs to.

Geometrically: on an extension node `x(1,2,:)` is the ray to the wall and `x(1,3,:)` the along-wall
direction, so `det -> 0` means **the ray runs nearly ALONG the wall** - grazing incidence, which is
what happens at a leg-end corner.

**Where they are.** All 51 degenerate boundary nodes (`det < 0.3`, out of 621) sit at exactly two
points:

    cluster A   R 1.2654 .. 1.2678   Z -1.0011 .. -0.9914     2.4 mm wide
    cluster B   R 1.6050 .. 1.6057   Z -1.1084 .. -1.1060     0.7 mm wide

both leg-end corners, with **types 4, 5 and 9 interleaved inside each one** - the fingerprint of the
type-9 grid surgery at `:1923-1949`, which removes "small edge triangles" and promotes their
neighbours to type 9 *without recomputing any frame for the new boundary edge* (and, unconditionally,
demotes any genuine type-5 wall node it touches). With `sheath_diag_R_split = 1.42`, cluster A is
**inner** (worst `1/det` = 108.7) and B **outer** (65.8).

**The inner corner is 1.65x worse - and every failure recorded for this BC is on the inner target**
(``sheath d collapse crash``: "LOCAL collapse on the INNER target", outer steady throughout;
``sheath multitype 1plus5 works``: inner `ePhi/kTe` ramps 3.96 -> 10.02, outer flat). The script
has no knowledge of which target fails, so that correlation is unforced. It is the first
explanation this campaign has had for the inner/outer asymmetry of the failures.

**And it explains why type 1 alone survives although the corners are in its mesh.** Type 1 has
**0.0 %** of its nodes below 0.3, so with only `bcs(1)` enabled *no replaced row is ever written on
a degenerate node* - they keep their ordinary Dirichlet. Enable 4, 5 or 9 and the weak route starts
writing a `zbig`-dominant constraint on `zj = Delta* psi` at nodes whose derivative basis is
conditioned 50-110x.

**Framing that should not be lost: this is a MESH defect, not a sheath defect.** The same degenerate
elements are used by every equation in the model. The sheath is simply the only one that writes a
dominant *replaced* row on a second-derivative quantity there, so an ill-posed representation
becomes a dominant wrong equation instead of a small error.

## 6.5 The fix: `sheath_weak_detmin`

    sheath_weak_detmin = 0.3
    sheath_weak_ufade  = .true.     ! REQUIRED alongside it, see below

A node whose frame determinant is below the threshold is taken out of the sheath **entirely** - it
becomes an ordinary non-sheath boundary node whatever its type says. Three mechanisms, all of which
already existed:

1. `mod_boundary_conditions.f90` applies **BOTH** its Dirichlets, `u` and `zj`.
2. `mod_boundary_matrix_open.f90:1090` skips its trace row, via the same `cycle` that already
   handled Dirichlet-u corners.
3. `sh_ufree` at `:218` sees it as frozen, so `sheath_weak_ufade` fades its Gauss points out of the
   free neighbour's row. That leak is real and was measured: `f6ef86fbd` cut type 4's `|Fd/S|max` at
   step 1 from **897 to 8.76**, a 102x reduction. Hence `ufade` is not optional here.

**Why this is not another `sheath_weak_wmin`.** Every previous gate keyed on the SOLUTION
(`|zj0/zj_sat|`, `D_a/D0_a`), so it could close in response to the very divergence it existed to
prevent. This one keys on a **static property of the mesh**: fixed from step 1, identical on every
rank, and the exact node set it removes is knowable *before the run* from
`util/check_boundary_frames.py`. It cannot feed back.

**Runtime confirmation that the gate is live:** the per-type `det min` in the sheath report must
come out `>= sheath_weak_detmin`, because every row below it is removed before accumulation. The
`N frame-gated` counter is a backstop and should read **0**.

## 6.6 A bug that was shipped and fixed - read this before editing the gate

The first version put the frame test into the `cycle` at the **accumulation** site only, reusing the
existing Dirichlet-u corner guard. That guard's comment says a skipped node "keeps its Dirichlet zj
row" - **true for the type-3/9 corners it was written for**, because `bcs(:)%dirichlet%zj = .true.`
is the DEFAULT. It is **false for a sheath-enabled type**, where the weak route requires
`dirichlet%zj = .false.` in the namelist. So a frame-gated type-4/5 node got **no zj row at all** -
the incomplete-equation runaway of section 3.3.

Fixed in `ab2c8caa8` by applying both Dirichlets in `mod_boundary_conditions`, which visits every
boundary node unconditionally. That matters rather than being tidier: a fallback inside the trace
accumulator *can* be bypassed, because `ufade` drives `wk_wgt` to zero on an edge with **both** ends
frozen and the accumulator skips a row with `wk_D <= 0` - certain at these corners, where the
degenerate nodes are contiguous (10 in a row on the inner one).

## 6.7 The diagnostic

    python3 util/check_boundary_frames.py jorek000000.h5

Straight off an HDF5 restart - no rebuild, no cluster. Prints per boundary type: the frame-chord
alignment `cos`, the determinant (min / p05 / median / mean, and % below 0.3 and 0.1), the fraction
of `grad psi` in each derivative DOF, a threshold sweep for choosing `detmin`, and the `(R,Z)` of
the worst nodes with an inner/outer split.

**HDF5 axis gotcha:** the wrapper does not preserve the Fortran argument order. `x`, declared
`(n_nodes, n_coord_tor, n_degrees, n_dim)`, lands as `(n_nodes, 2, 4)`; `values` lands as
`(n_nodes, n_var, n_degrees)`, i.e. `values[:, var, deg]` with psi = var 0. The script detects the
axes by size and prints the raw -> resolved mapping.

---

# 7. Current status and what to do next

## 7.1 NEW, 2026-09-03: 1+4 converges for the first time

With `sheath_weak_detmin = 0.3`, on **1+4**, which had never survived past 8 steps:

    output    weak     |Fd/S| type 1   |Fd/S| type 4
       0     1.934          673            8.76
       3     0.230           40.5          3.32
       6     0.136        1.15e5           1.35     <- transient spike, recovered
      10     0.106            1.02         0.254

Weak residual **18x down and monotone**; type 4's `|Fd/S|` **34x down**. Gate confirmed live:
`bnd type 4 ... det min = 0.304` (>= 0.30) with rows down from 220 to 152, and `0 frame-gated` as
designed. Type 1 unaffected: `det min = 0.715`, 196 rows, nothing removed.

**And this was run with `sheath_init_u = .true.`** - the flag documented in 5.6 as harmful (219 -> 8
on an earlier 1+4). It converged anyway. **Re-run with `sheath_init_u = .false.`; it should be
better still.** That is the single highest-value next run.

Caution before celebrating: `ePhi/kTe max` rises monotonically 3.00 -> 6.03 over these outputs, and
the inner target goes 2.97 -> 4.53 while the outer falls 2.97 -> 2.08. The BC is being satisfied
better and better, but the plasma state is still evolving away from the restart. Watch whether the
inner target settles or keeps climbing - a monotone inner ramp with a flat outer is the 1+5 signature
that ended at 305.

## 7.2 Run order

1. **1+4 with `sheath_init_u = .false.`** (in flight config otherwise unchanged).
2. **1+5 with `detmin = 0.3`** - queued. Expect ~8 rows gated, `det min >= 0.3` on type 5.
3. **1+4+5+9** if both hold. Type 9 is 100 % gated, i.e. reverts to its pre-sheath Dirichlet
   (6 nodes, 0.039 m^2), so the full multi-type namelist becomes usable.

## 7.3 Open items, ranked

* **`min |xjac|` diagnostic.** `det` at a node is a PROXY. The quantity the equations feel is the
  element mapping Jacobian: `mod_boundary_matrix_open` computes `xjac = x_s y_t - x_t y_s` and
  `grad_t = (-y_s, x_s)/xjac` at every boundary Gauss point, and `grad_t` carries the `1/xjac` into
  the normal-derivative terms. At the corner the Jacobian is `size(:,2) x2 wedge size(:,3) x3`, so
  it is proportional to `det(x2,x3)` **up to the two `element%size` factors**. A per-type
  `min|xjac|` is nearly free and should come before serious threshold tuning. NOT IMPLEMENTED.
* **The grid fix.** Gram-Schmidt `x(1,3,:)` against `x(1,2,:)` for extension nodes
  (`grid_xpoint_wall.f90:1252`, `:1311`). This changes `element%size` through `:1751` and therefore
  the mapped geometry, so it needs a **fresh equilibrium**, not a restart. Note the current frames
  are chord-aligned everywhere (`cos ~ 1`) and any fix must not break that. `detmin` is containment;
  this is the cure.
* **Repair the type-9 surgery** at `:1923-1949`: add `exit` after a hit (an element with two type-9
  vertices is currently removed twice) and guard the `= 99` write so it cannot demote a genuine
  wetted type-5 node. Neither moves a node, so no new equilibrium is needed.
* **Normalise the psi column of the weak row by local `|grad psi|`.** Addresses 6.3(a) without
  touching the grid, and is why `sheath_psi_jacobian` never discriminated.
* **`wk_dfac` is computed and never used** (`mod_boundary_matrix_open.f90:672-677`); its declaration
  at `:86` claims every Jacobian column carries it. Defensible as a trust region (bounded `F` with a
  full-strength `J` gives a *smaller* step), but code and comment disagree - fix one of them.
* **`direction(3)` is never assigned in the `(9,9)` branch** of `construct_matrix_mod.f90:137-144`
  while every other branch assigns it. **Dead at `n_order = 3`**; a live bug above it.
* **Is `use_simple_bnd_types` set?** It defaults `.false.` (`preset_parameters.f90:826`) and
  `max_bnd_types = 30`, yet the logs only ever show types 1/3/4/5/9. If types **11/15** exist as
  distinct node types they are the OUTWARD-field halves of the type-1 and type-5 targets, carry real
  area, and have never had a `bcs()` block. Cheap to check from the namelist; expensive to be wrong.

## 7.4 Things that are settled - do not re-litigate

Recorded so they are not re-derived a fourth time. Each has been claimed otherwise by at least one
review.

* There is **no Newton iteration** (section 2.3).
* The **`n_wall_blocks = 0` else-branch** is what runs in `construct_matrix_mod` (section 3.4).
* **The grid is not rebuilt on restart**, so no grid-side property can be the trigger for something
  that happens at step 219 - though it can explain why a surface is awkward.
* **`sheath_min_bn` is a no-op on the weak route** - `sh_wgt_bn` lives inside
  `if (apply_natural_bc(var_u))`, which the weak route requires `.false.`. It survives only in the
  diagnostic weight, so `gated-off area` in the output carries no information about the weak row.
* **A uniform weight cancels exactly** out of a replaced row (section 5.4).
* **`Bdotn` in the diagnostics is dimensional `B.n`, not `b_n = B.n/|B|`.** With `Btot ~ 2`, divide
  the printed range by ~2 before comparing against `sheath_min_bn`.
* **The frame-chord / `element%size` `sign()` hypothesis is FALSIFIED** - `cos ~ 1.0000` on every
  type, worst 0.8584 on type 5, which works.
* **Type 1 alone runs to TIMEOUT, not to a crash.** The "1350 steps" figure quoted through much of
  this campaign came from a live progress report, not a run limit.
