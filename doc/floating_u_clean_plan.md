# Floating potential boundary condition, clean implementation plan

Branch `floating-u-clean`, off `develop` at 46882606a. Written 2026-09-14.

## Goal and definition of done

A floating-potential boundary condition on the model600 electrostatic potential,

    Phi_plasma - V_wall = Lambda * k_B*Te / e        (Lambda = 3 by default)

that works out of the box: no parameter scans, no stabilisers, no clips, no thresholds, no floors
added for its sake. Any reasonable model600 divertor setup switches it on with the boundary-type flags
and nothing else.

Done means, on the reference case with the DEFAULT namelist values of every new switch:

1. Lambda = 0 reproduces the plain Dirichlet `u = 0` run to solver tolerance.
2. Lambda = 3 on every connected wall segment runs longer than every tuned run of the old
   `floating-u` branch (best 1100 steps), with `min_sheath_angle` and `D_perp_sc_num` at their
   develop defaults.
3. rho, Ti, Te stay positive everywhere, including between nodes, without a floor doing it.
4. Every new Jacobian column is finite-difference checked in `tests/floating_transport`, and the
   check fails when the column is removed.
5. Particle and energy balance at the wall: the fluid wall flux, the sheath energy flux and the
   kinetic recycling flux all use the same total normal flow.

## Why the potential row is not the problem

The old branch imposed the row to 1e-12 V and still died. The imposed potential follows Te, so it has
a tangential gradient of order Lambda*dTe/dl at every strike point. That gradient drives an ExB flow
normal to the wall of the same order as, and on one side of every Te peak opposite to, the parallel
sonic outflow. This is physics, and it stays under the full sheath current characteristic. Measured
on the old branch at step 500 (`doc/pp_out.txt`): |vE.n| up to 40 km/s, larger than the parallel
normal flow on 8% of wall points, net inflow on 10% of wall points.

JOREK's density and temperature equations advect with the undifferentiated test function, so they
carry no boundary flux at all, and the natural boundary rows only know the parallel flow. A wall
point with net inflow therefore has no inflow datum: the boundary is ill posed there. That, plus a
time update that never checks positivity, is what the campaign kept hitting. The plan below fixes
those two things and makes the four wall channels (momentum, particles, energy, kinetic recycling)
see one and the same total normal flow. Everything else is scaffolding.

Nothing from the old branch's stabilisation history is carried over: no drift saturation, no SOLPS
clip, no zero-sum stabiliser, no density diffusion sensor, no incidence floor, no grazing drop, no
T_min branch in the potential row, no current-profile mask.

## Conventions pinned once

- Phi = +F0*u. JOREK's (R,Z,phi) basis is right-handed; the poloidal ExB velocity the fluid advects
  with is v_E = (-R*u_Z, +R*u_R) (`mod_elt_matrix_fft.f90` fu_velocity, `particles/mod_fields.f90`
  581-582).
- Outward unit normal from `normal_direction` (node minus opposite vertex), never from the tangent
  alone. Everything below is outward-positive.
- On a wall Gauss point with tangent (x_s, y_s)/dl: vE.n = -orient*R*u_s/dl where
  orient = sign(y_s*n_R - x_s*n_Z). Parallel normal flow: Vpar*(B_pol.n) = Vpar*(psi_y*n_R - psi_x*n_Z)/R,
  a velocity because Vpar is v/|B|. Unit incidence |b.n| = |B_pol.n|/|B|.
- Total normal flow: vn = Vpar*(B_pol.n) + vE.n.
- Per-Gauss-point branch decisions use the n=0 state. One solve per step, no Newton loop, so every
  column of a residual must be assembled: a lagged column is a per-step error, not a convergence cost.

## Items, in dependency order

### 1. Potential row: `bcs(i)%floating_u`

New: `bcs%floating_u` (type_bcs), `sheath_Lambda` (default 3), `sheath_V_wall` (default 0, volts),
`models/model600/mod_floating_u.f90` with `floating_u_norm(a_n, C_T, C_V)`,
a_n = 2*e*F0*sqrt(mu0*rho0)/m_i, C_T = 2*Lambda/a_n (halved under single-T), C_V = sqrt(mu0*rho0)/F0,
and a self-test for both signs of F0.

Row: in the Dirichlet loop of `mod_boundary_conditions.f90`, when `k == var_u` and the type is
floating, every u trace DOF (value and tangential derivatives, the same `iv_dir` loop as Dirichlet)
gets the affine relation u = C_T*Te + C_V*V_wall: diagonal zbig, column zbig*(-C_T) on the same Te
DOF, RHS -zbig*(u - target). V_wall enters the value DOF of the n=0 harmonic only. No temperature
floor and no branch in this row; item 5 owns positivity.

Files: `phys_module.f90`, `preset_parameters.f90`, `mod_log_params.f90`, `initialise_parameters.f90`
(namelist read + broadcast), `mod_boundary_conditions.f90`, new `mod_floating_u.f90`.
Anchor: Lambda = 0 is Dirichlet u = 0 exactly.

### 2. Total-flow Bohm condition, weak: `mach1_weak`

The develop nodal Mach row carries `factor/Btot*R^2*u_b/psi_b` in its value row only, with no
counterpart in the slope row and a division by psi_b that diverges at grazing incidence. It cannot
be used with a wall potential that varies along the wall. It is left untouched for everyone else.

New flag `mach1_weak` (default false; setup aborts if any `floating_u` is set without it, so the
combination is validated rather than tuned). It replaces the nodal Mach rows on `mach1` types with
one Galerkin residual per boundary Gauss point, assembled in `mod_boundary_matrix_open.f90` next to
the existing `mach_one_bnd_integral` penalty:

    res = (B_pol.n)*Vpar - max( cs*|b.n| - vE.n , 0 )
    row = Zbig * (B_pol.n) * res        (weight = d res / d Vpar)

Columns: Vpar (Zbig*(B_pol.n)^2), Ti/Te through cs, and the u trace DOFs through vE.n on the active
branch. The Vpar column carries (B.n)^2, so a grazing point loses authority continuously and the
natural Vpar condition takes over; there is no threshold and no division. One-sidedness is the Bohm
inequality. With `mach1_weak` the nodal block in `mod_boundary_conditions.f90` is not entered at all
(no expression with 1/psi_b is evaluated), while `apply_cs` still suppresses the Dirichlet Vpar row.
Type 3 keeps no Mach row, as today.

Diagnostic: per boundary type, the normalised weighted moment and pointwise residual, min/mean/max,
with the location of the max.

### 3. One total normal flow in the particle, energy and recycling channels

- Sheath energy flux (`mod_boundary_matrix_open.f90`): the (gamma_sheath - 1)*rho*T*vn collection
  uses vn_out = max(Vpar*(B_pol.n) + vE.n, 0) instead of the parallel part alone, with columns on
  rho, T, Vpar and the u trace DOFs. Closed branch (vn <= 0) collects nothing. The c_angle term is
  unchanged and stays at its develop default.
- Density reflection term: same substitution where `density_reflection` is nonzero.
- Kinetic recycling (`particles/mod_particle_wall_interaction.f90`): Gamma_d from
  max(v_par.n + v_ExB.n, 0), with `calc_EBpsiU` returning the fluid v_ExB = R grad(u) x e_phi as an
  optional trailing argument. This keeps the fluid loss and the kinetic source equal on every face.
- The existing `calc_NeTevpar` single-T normalisation under WITH_TiTe is fixed in passing if it is
  confirmed on develop (it halves Te in the recycling projection).

### 4. Weak inflow closure for rho, Ti, Te (new)

Where the total normal flow is inward the boundary needs a datum. Boundary-only, no penalty
parameter: the standard upwind weak inflow term added to the natural rows,

    R_X(v) +=  - oint  v * min(vn, 0) * (X - X_in) dl

for X in {rho, Ti, Te}. Its coefficient is the physical inflow rate, it vanishes identically where
the flow is outward, and it is consistent (zero when X = X_in). Inflow states at a material wall:
rho_in = 0 (no plasma enters from a wall), and for the temperatures X_in = X (no term: the
incoming energy flux is zero once the incoming density is zero). Columns: X (|vn|), u trace DOFs and
Vpar (through vn, times X - X_in). Branch on the n=0 vn per Gauss point, exact derivative of min.

This term is what makes the boundary well posed at the grazing wall segments where item 2 has,
by design, no authority. It is expected to be active on a few percent of the wall.

Verification: FD every column; check that with u = 0 and no parallel inflow the term is
identically zero, so item 10a is unaffected.

### 5. Admissible update (new)

`core/mod_jorek_timestepping.f90`: after the solve and before `update_values`, evaluate the
candidate rho, Ti, Te at every node and at the Gauss points of every element (the trace is cubic,
so nodal positivity is not enough). If any is <= 0: reject the step, halve dt, re-assemble and
re-solve; on success recover dt geometrically toward the requested value. Flag
`admissible_update` (default true when any `floating_u` is set). Gears history must not advance on a
rejected step; check how `index_now` and the previous-step storage interact with a re-solve. Print
where the minimum was and how many rejections occurred. This is not a floor: a rejected state is
never written.

### 6. Exterior sides only (new)

`construct_matrix_mod.f90` gives the open-boundary integral to any element side whose two endpoints
are boundary nodes. Build the exterior-side table from mesh connectivity once per matrix
construction (`mod_floating_boundary_edges.f90`, conforming unrefined meshes, aborts otherwise) and
skip non-exterior sides whenever `floating_u` or `mach1_weak` is active. Report once how many sides
were skipped; zero is the expected answer and proves the point cheaply.

### 7. Boundary diagnostics and postproc

Per timestep table, per boundary type: |u - u_float| in volts, vE.n min/max (signed, outward),
inflow fraction of wall length, weak Mach moment, min rho and min T with (R,Z). Postproc
expressions: `Phi`, `Te_float`, `vparB_norm`, `vExB_norm`, `vflow_norm`, `bnd_type`, `bnd_dl` along
the boundary, so a run can be judged from a table rather than a lifetime.

### 8. Tests, `tests/floating_transport`

gfortran serial harness (no MPI): fixtures for a wall element, FD Jacobians of every column added in
items 2, 3, 4 (verified to fail when a column is dropped), the B.n sign crossing of the weak Mach row,
vE.n against `mod_fields`, fluid and kinetic wall flux agreement, the Lambda = 0 reduction, and the
`check_imports.py` guard for `only:` lists and cross-extent assignments in files the harness
cannot compile.

### 9. Documentation

`doc/floating_u.md`: what the condition is and is not (zero local current, no thermoelectric current),
the namelist, the conventions above, the diagnostics, and the ladder results. The derivation of a_n,
C_T, C_V with each step referenced to a file and line.

### 10. Validation ladder (cluster)

a. `mach1_weak` + items 3-6, Lambda = 0: must match the develop Dirichlet-u run to solver tolerance.
   Any difference is a bug in 2-6, not physics.
b. Lambda = 3, all connected wall types, default namelist otherwise. Done criteria 2 and 3.
c. Same with kinetic neutrals on. Done criterion 5 from the wall-flux table.
d. dt halved and mesh refined at the targets: the inflow fraction and the wall fluxes must converge.
   If they do not, item 4's inflow state is wrong, not under-resolved.

## Status (2026-09-14)

| item | state | commit / test |
|---|---|---|
| 1 potential row | done | `mod_floating_u.f90`, self-test both field signs |
| 2 weak Bohm row | done | FD every column, negative controls, B.n sweep |
| 3 one total flow | done | Ti/Te rows FD, closed wall, kinetic recycling on the same flow; calc_NeTevpar Te fix |
| 4 inflow closure | done | rho rows FD, absent for outward flow |
| 5 admissible update | REMOVED (user decision): the run must go through the production ramp untouched; positivity is a property of the equations, not of a step policy | commits 8466ca5fb..ab6153da5, removed after |
| 6 exterior sides | done | connectivity test |
| 7 diagnostics | done (wall table); boundary postproc expressions not ported | equations untouched |
| 8 tests | done | `tests/floating_transport/run.sh` |
| 9 documentation | done | `doc/floating_u.md` |
| 10 ladder | not started: needs the cluster build | |

Ladder step a is expected to differ from develop only where the inflow closure is active
(grazing points with inward parallel flow, which the weak row leaves natural); everywhere else
the rows are identical to develop with u = 0, as the harness checks.

## Revision after the first cluster runs (2026-09-14, evening)

Facts from the runs: the build works, the weak Bohm row holds the 12.5 km/s inner-strike-point
drift at mom ~1e-7, inflow is 2% of the type-5 wall and zero elsewhere, and the first negative
value (Ti, -3e-4 eV, location unknown) appeared at step 438. The step-rejection policy (old item 5)
was removed: the run must go through the production ramp
`tstep_n = 1e-3, 1e-2, 1e-1, 0.3, 1, 2, 10` untouched. Positivity therefore has to be a property of
the discretisation. Three mechanisms can break it, none of them a boundary-condition matter, and
each has a parameter-free treatment:

### R1. The ramp switches drive BDF2 outside its stability bound

JOREK's Gears scheme is the exact variable-step BDF2 (checked: zeta = dt/(dt+dt_prev) with the
history scaled by dt/dt_prev reproduces the (1+2r)/(1+r), -(1+r), r^2/(1+r) coefficients). A
variable-step BDF2 is zero-stable only for step ratios r < 1+sqrt(2) = 2.41. The ramp has ratios
10, 10, 10, 3, 3.3, 2, 5. At a switch with ratio r, any field that decayed by more than
(1+r)^2/r^2 - 1 over the previous step is extrapolated NEGATIVE regardless of the right-hand side:
at r = 5 a 31% drop per step suffices (numerator 6*y1 - 25/6*y0 < 0 for y1/y0 < 25/36). A cold
strike-point layer with a 13 us step and a 10 us parallel transit does that; a wall with u = 0
does not, which is why develop never saw it. Every old floating-u run died at steps 598-608, the
2 -> 10 switch.

Treatment (standard for variable-step multistep codes, no parameter): take the first step after
a ratio beyond 1+sqrt(2) as a one-step implicit Euler step (zeta = 0), then resume BDF2. One flag
set by the timestepper, read where zeta is formed (mod_elt_matrix_fft, mod_boundary_matrix_open).
Changes nothing for ramps with ratios below the bound.

### R2. Cell Peclet number of the ExB advection at the strike point

Continuous Galerkin without upwinding is non-monotone for advection once Pe_h = v*h/D > 2. The
imposed potential makes Pe_h large exactly at the strike point: with the measured 12.5 km/s,
h = 2.9 mm and the D_perp of the reference case (5.3 m^2/s), Pe_h = 6.9. This is the dispersive
density dipole the earlier branch saw (one-cell negative spot at the strike point one step before
Ti and w react), and the reason `use_sc` with `D_perp_sc_num = 10` was the only thing that ever
gave long runs. That coefficient is a dimensional tuning knob and is out.

Treatment, in order of preference:
  (a) Resolution: refine the target region until Pe_h < 2, i.e. h < 2*D_perp/v_E, about 0.8 mm
      here. A grid criterion the user controls; no code. Cost: elements.
  (b) Streamline-upwind stabilisation of the poloidal (ExB) advection of rho, Ti, Te with the
      standard optimal coefficient tau = h/(2|v|) * (coth Pe_h - 1/Pe_h), applied along the flow
      only. No free constant; it vanishes as Pe_h^2 where the mesh resolves the flow, so it is
      inactive everywhere a develop run is resolved. Streamline-diffusion form first (a few lines
      in the rho and T rows, FD-testable); consistent SUPG weighting of the full residual only if
      the diffusion form measurably smears the target profiles.
  (c) `use_sc` as it exists, with its coefficient scanned. Rejected by principle.

### R3. Locate every first violation

Add the volume minima of rho, Ti, Te with (R,Z) to the `[floating_u]` table (information only,
no control), so the step-438 Ti zero and any later one are attributed to a wall, an outer
boundary profile or the strike-point layer before anything is changed.

### What stays

Items 1-4 and 6-9 are unchanged and verified. The definition of done is unchanged: default
namelist, production ramp, no stabiliser coefficient, no floor, Lambda 3 on every wall type.

## Order of work

1 and 2 first (they define the potential and the momentum channel), then 6 (removes an unknown from
every later run), then 4, then 5, then 3, then 7-9. Ladder step a after 6, step b after 5.
