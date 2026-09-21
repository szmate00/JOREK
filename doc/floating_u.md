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
mach1_weak_drift     = .t.   ! SOLPS non-marginal Bohm condition ...
mach1_weak_drift_cut = .t.   ! ... except where the field grazes the wall
sheath_Lambda        = 3.d0
```

`dirichlet%u` stays `.true.` on those types. Nothing else changes: the production timestep ramp,
`min_sheath_angle`, `D_perp_sc_num` and every other parameter stay as in a develop run. Reference case:
2000+ steps from the equilibrium through the ramp `1e-3 .. 10`.

## What is assembled

1. **Potential row** (`mod_boundary_conditions.f90`, `mod_floating_u.f90`): u = C_T*Te + C_V*V_wall on
   every u trace DOF, replacing the Dirichlet rows. C_T = 2*Lambda/a_n (halved in a single-T build),
   C_V = sqrt(mu0*rho0)/F0, a_n = 2*e*F0*sqrt(mu0*rho0)/m_i.
2. **Weak Bohm condition** (`mach1_weak`, `mod_boundary_matrix_open.f90`): one residual per wall Gauss
   point on edges whose both endpoints are `mach1` types,
   `res = (B_pol.n)*Vpar - target`, weight `d(res)/d(Vpar) = B_pol.n`, so it fades as (B.n)^2 at grazing
   incidence and the natural Vpar condition takes over; no threshold, no division. `target = cs*|b.n|`
   (marginal), or `max(cs*|b.n| - vE.n, 0)` with `mach1_weak_drift` (the parallel flow supplies the
   outward normal flow the drift does not, never asked to reverse). With `mach1_weak_drift_cut` the drift
   is not compensated where |b.n| < sin(min_sheath_angle): the angle below which the sheath fluxes already
   come from the c_angle floor model. The nodal Mach1 rows are not assembled under `mach1_weak`; type 3
   gets no Dirichlet Vpar row. Columns on Vpar, Ti, Te and the u trace are exact; the |B| dependence on
   the free normal psi derivative is lagged (trace-DOF loop).
   Why weak: the nodal row's drift term sits in the value row only and divides by psi_b; it cannot take a
   wall potential that varies along the wall. Why the cut: compensating the drift where the field grazes
   demands an unbounded parallel flow (Mach 400 measured); a smooth 2cs bound failed at the same wall;
   leaving Vpar free there runs away within tens of steps. All measured on this case.
3. **One total normal flow** `max(Vpar*(B_pol.n) + vE.n, 0)` in the sheath energy transmission and density
   reflection rows (exact u, Vpar, psi columns). Off the weak route the expressions are unchanged. The
   kinetic recycling flux is not changed here (to follow once upstream PR #32 is merged).
4. **Exterior sides only** (`mod_boundary_edges.f90`): the open-boundary integral is skipped on interior
   sides with two labelled endpoints (table from connectivity, once per matrix construction).

Conventions: outward-positive normal flows; `vE.n = -orient*R*u_s/dl`, `orient = sign(y_s*n_R - x_s*n_Z)`
for the edge tangent `(x_s, y_s)/dl`; `vpar0*ps0_s*normal_sign3 = Vpar*(B_pol.n)*R*dl`; branches decided on
the n=0 state per Gauss point and differentiated exactly on the branch taken.

## Not covered

`n_order >= 5` trace DOFs beyond value and first derivative; boundary postproc expressions along the wall.
