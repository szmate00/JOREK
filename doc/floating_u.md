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
mach1_weak           = .t.   ! required
mach1_weak_drift     = .t.   ! drift-compatible Bohm condition, drift bounded to +-2 cs|b.n| (SOLPS-ITER BCMOM=13) ...
mach1_weak_drift_cut = .t.   ! ... except where the field grazes the wall
sheath_Lambda        = 3.d0
```

`dirichlet%u` stays `.true.` on those types. Nothing else changes: the production timestep ramp,
`min_sheath_angle`, `D_perp_sc_num` and every other parameter stay as in a develop run.

Status: the earlier 2000+ step reference ran with the kinetic recycling clipped to INFLOW (inward
`wall_normal_vector`), i.e. target recycling at the c_angle floor. With the recycling on the outward normal
the same namelist died at step 729 (wall injecting density at ExB inflow points); with the sheath-set flux at
808 (density positive, same spot); with the flux consistent with the cut at 700: the unbounded drift
compensation asked for Mach 105 at the cold strike-point front (measured). The drift bound is the fix for that.
Not yet run in this form.

## What is assembled

1. **Potential row** (`mod_boundary_conditions.f90`, `mod_floating_u.f90`): u = C_T*Te + C_V*V_wall on
   every u trace DOF, replacing the Dirichlet rows. C_T = 2*Lambda/a_n (halved in a single-T build),
   C_V = sqrt(mu0*rho0)/F0, a_n = 2*e*F0*sqrt(mu0*rho0)/m_i.
2. **Weak Bohm condition and wall fluxes** (`mach1_weak`, `mod_boundary_matrix_open.f90`), on edges whose both
   endpoints are `mach1` types, the same condition in every part at every Gauss point:
   - **Vpar row**: `res = (B_pol.n)*Vpar - target`, weight `B_pol.n`, one residual per Gauss point on the
     Vpar trace. With `mach1_weak_drift` (SOLPS-ITER BCMOM=13 non-marginal, wide-grid form):
     `target = cs*|b.n| - vE_r` with the drift bounded, `vE_r = clamp(vE.n, -2 cs|b.n|, +2 cs|b.n|)`: the
     parallel flow supplies the outward normal flow the ExB drift does not, but is never asked for more than
     2 cs of correction (target within `[-cs|b.n|, 3 cs|b.n|]`, so the demand never diverges at grazing
     incidence or at a cold strike-point front, where the unbounded form asked for Mach 100). The u column of
     the target vanishes while the bound is active. Without `mach1_weak_drift`: `target = cs*|b.n|`. With
     `mach1_weak_drift_cut` the drift is not compensated where `|b.n| < sin(min_sheath_angle)`. The nodal
     Mach1 rows are not assembled; type 3 gets no Dirichlet Vpar row. Why weak: the nodal row's drift term
     sits in the value row only and divides by psi_b.
   - **Wall fluxes**: the volume advection of rho, rho*Ti and rho*Te is not integrated by parts, so it carries
     an implicit wall flux `q*vn*R*dl` in both directions (`vn = Vpar*(B_pol.n) + vE.n`), with no inflow
     datum where `vn < 0`. It is replaced by `Gamma = n*max(vn, v_fl)`, energy `gamma_sh*T*Gamma`, with
     `v_fl = cs*|b.n|` where the row compensates the drift (Gamma is the Bohm flux wherever the row holds) and
     `v_fl = 0` where it does not (below the cut, or the marginal row): there the wall removes what flows out
     and supplies nothing where the flow is inward. In the rows: `fx_n = max(vn, v_fl)*R*dl` plus the excess
     sink `-q*max(v_fl - vn, 0)*R*dl`, proportional to q; boundary energy term
     `oint q^2 (vn/2 - max(vn, v_fl)) <= 0` for every vn.
   - **Kinetic recycling** (`mod_particle_wall_interaction.f90`): exactly what the fluid loses at that wall
     point: `n*max(vn, cs*|b.n|)` where the fluid row compensates the drift (edge with both endpoints `mach1`
     types, `mach1_weak_drift`, not below the cut), `n*max(vn, 0)` otherwise, plus the c_angle floor, all
     with the fluid's sound speed `sqrt(gamma*(Ti+Te))` (fluid Ti via `calc_NeTeTi`), on the outward normal
     (`wall_normal_vector` points inward). `calc_EBpsiU` returns the fluid ExB velocity. `calc_NeTevpar`
     reads Te itself in a two-temperature model600 build (it read Ti/2).
   Exact columns on rho, Ti, Te, Vpar, u and the psi trace; the Btot dependence on the free normal psi
   derivative is lagged (trace-DOF loop), as are the `corr_neg` derivatives of cs (pre-existing).
3. **Diagnostics** (`wall_diag`, default on under `mach1_weak`; `wall_diag_every`; `mod_wall_diag.f90`): two
   tables per step in the log, one line per boundary type with (R,Z): `[floating_u]` (floating-row residual
   |u - C_T*Te - C_V*V_wall| in V, largest outward and inward vE.n, min rho, min Te, inflow fraction, share of the
   wall loss supplied by the excess sink) and `[mach1_weak]` (Bohm-row moment, pointwise |res| min/mean/max in
   m/s with location, max Mach with location, max target, max cs|b.n|, max unbounded drift demand
   |vE.n|/(cs|b.n|), fractions with drift compensation and with the bound active). With
   `wall_diag_profile_every > 0` the full wall profile `[wall prof]` every that many steps and on a collapse of
   the wall minimum of rho or Te.
4. **Exterior sides only** (`mod_boundary_edges.f90`): the open-boundary integral is skipped on interior
   sides with two labelled endpoints (table from connectivity, once per matrix construction).

Conventions: outward-positive normal flows; `vE.n = -orient*R*u_s/dl`, `orient = sign(y_s*n_R - x_s*n_Z)`
for the edge tangent `(x_s, y_s)/dl`; `vpar0*ps0_s*normal_sign3 = Vpar*(B_pol.n)*R*dl`; branches decided on
the n=0 state per Gauss point and differentiated exactly on the branch taken.

## Not covered

`n_order >= 5` trace DOFs beyond value and first derivative; boundary postproc expressions along the wall.
