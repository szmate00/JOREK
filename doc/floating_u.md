# Floating-potential boundary condition (model600)

## What it is

`bcs(i)%floating_u = .true.` imposes on every u trace DOF of boundary type i

    Phi - V_wall = Lambda * k_B*Te / e        i.e.        u = C_T*Te + C_V*V_wall

with `sheath_Lambda` (default 3) and `sheath_V_wall` in volts (default 0): the zero-local-current limit of
the sheath characteristic. It carries no net wall current. It gives a wall potential that follows Te,
hence an ExB drift normal to the wall wherever Te varies along it, of the order of the sonic outflow at
the strike points. The row is affine and exact. `Phi = +F0*u` (right-handed (R,Z,phi)); `a_n` carries the
sign of F0 so the physical potential is independent of the field direction (self-test at setup).

## Namelist

```fortran
bcs(1)%floating_u = .t.      ! every connected wall segment: 1, 3, 4, 5, 9 on the usual grids
bcs(3)%floating_u = .t.
bcs(4)%floating_u = .t.
bcs(5)%floating_u = .t.
bcs(9)%floating_u = .t.
mach1_weak        = .t.      ! required: the drift-compatible Bohm condition
sheath_Lambda     = 3.d0
```

`dirichlet%u` stays `.true.` on those types. Nothing else changes: the production timestep ramp,
`min_sheath_angle`, `D_perp_sc_num` and every other parameter stay as in a develop run. The former switches
`mach1_weak_drift` and `mach1_weak_drift_cut` are gone (remove them from the namelist): the Bohm condition is
drift-compatible everywhere, with no marginal form and no angle below which the drift is ignored.

Status: the earlier 2000+ step reference ran with the kinetic recycling clipped to INFLOW (inward
`wall_normal_vector`), i.e. target recycling at the c_angle floor. With the recycling on the outward normal
the same namelist died at step 729; with the sheath-set flux but the grazing cut still on, at 808 (same
inner-strike-point depletion, positive density). The cut made the Vpar row and the wall flux disagree below
the angle; both now impose the same condition. Not yet run in this form.

## What is assembled

1. **Potential row** (`mod_boundary_conditions.f90`, `mod_floating_u.f90`): u = C_T*Te + C_V*V_wall on
   every u trace DOF, replacing the Dirichlet rows. C_T = 2*Lambda/a_n (halved in a single-T build),
   C_V = sqrt(mu0*rho0)/F0, a_n = 2*e*F0*sqrt(mu0*rho0)/m_i.
2. **Drift-compatible Bohm condition** (`mach1_weak`, `mod_boundary_matrix_open.f90`): the total normal flow
   into the sheath `vn = Vpar*(B_pol.n) + vE.n >= cs*|b.n|`, imposed in the same form at every wall Gauss
   point of edges whose both endpoints are `mach1` types, by three consistent parts:
   - **Vpar row**: `res = (B_pol.n)*Vpar - max(cs*|b.n| - vE.n, 0)`, weight `B_pol.n`, one residual per
     Gauss point on the Vpar trace (SOLPS non-marginal form): the parallel flow supplies the outward normal
     flow the ExB drift does not and is never asked to reverse (target 0 where the drift alone exceeds
     `cs*|b.n|`). The nodal Mach1 rows are not assembled; type 3 gets no Dirichlet Vpar row. Why weak: the
     nodal row's drift term sits in the value row only and divides by psi_b.
   - **Wall fluxes**: the volume advection of rho, rho*Ti and rho*Te (parallel and ExB) is not integrated by
     parts, so it carries an implicit wall flux `q*vn*R*dl` in both directions, with no inflow datum where
     `vn < 0`; central Galerkin gains energy `-1/2 oint q^2 vn` there. The same condition is imposed on the
     flux, pointwise: `Gamma = n*max(vn, cs*|b.n|)`, energy `gamma_sh*T*Gamma`, i.e. the sheath rows on
     `fx_n = max(vn, cs*|b.n|)*R*dl` plus the excess sink `-q*max(cs*|b.n| - vn, 0)*R*dl`. Zero wherever the
     row holds; where it does not pointwise (the row fixes moments), a sink proportional to q. Boundary
     energy term `oint q^2 (vn/2 - max(vn, cs|b.n|)) <= 0` for every vn.
   - **Kinetic recycling** (`mod_particle_wall_interaction.f90`): exactly what the fluid loses at that wall
     point, `n*max(vn, cs*|b.n|)` on edges whose both endpoints are `mach1` types (the same edge test as the
     fluid), `n*max(vn, 0)` elsewhere, plus the c_angle floor, all with the fluid's sound speed
     `sqrt(gamma*(Ti+Te))` (fluid Ti via `calc_NeTeTi`), on the outward normal (`wall_normal_vector` points
     inward). `calc_EBpsiU` returns the fluid ExB velocity. `calc_NeTevpar` reads Te itself in a
     two-temperature model600 build (it read Ti/2).
   Exact columns on rho, Ti, Te, Vpar, u and the psi trace; the Btot dependence on the free normal psi
   derivative is lagged (trace-DOF loop), as are the `corr_neg` derivatives of cs (pre-existing).
3. **Diagnostics** (`wall_diag`, default on under `mach1_weak`; `wall_diag_every`, `wall_diag_profile_every`;
   `mod_wall_diag.f90`): per step in the log, `[wall]` per boundary type (inflow and sub-Bohm fractions, sink
   share, Bohm-row moment, extreme vE.n with location, max Mach, potential range), `[wall flow]` (ranges of
   Vpar*B.n, vE.n, vn, cs|b.n|, min |b.n|, max |dPhi/ds|), `[wall min]` (min rho, Ti, Te with location),
   `[wall pt]` (full local state at the min-rho, min-Te, max-|vE.n| and max-sink points), `[volume]` (node
   minima of rho, Ti, Te and max Te with location and boundary type), and `[wall prof]`, the full wall
   profile, automatically on any step where the wall minimum of rho or Te is non-positive or has halved.
4. **Exterior sides only** (`mod_boundary_edges.f90`): the open-boundary integral is skipped on interior
   sides with two labelled endpoints (table from connectivity, once per matrix construction).

Conventions: outward-positive normal flows; `vE.n = -orient*R*u_s/dl`, `orient = sign(y_s*n_R - x_s*n_Z)`
for the edge tangent `(x_s, y_s)/dl`; `vpar0*ps0_s*normal_sign3 = Vpar*(B_pol.n)*R*dl`; branches decided on
the n=0 state per Gauss point and differentiated exactly on the branch taken.

## Not covered

`n_order >= 5` trace DOFs beyond value and first derivative; boundary postproc expressions along the wall.
