# Floating potential boundary condition (model600)

Branch `floating-u-clean`. Plan and status: `doc/floating_u_clean_plan.md`.

## What it is

`bcs(i)%floating_u = .true.` imposes, on every u trace DOF of boundary type i,

    Phi - V_wall = Lambda * k_B*Te / e        i.e.        u = C_T*Te + C_V*V_wall

with `sheath_Lambda` (default 3) and `sheath_V_wall` in volts (default 0). This is the zero
local-current limit of the sheath characteristic: an insulating wall, or a conducting one on
which no net current flows. It carries no thermoelectric current. What it does give is a wall
potential that follows Te, hence a tangential electric field and an ExB drift normal to the wall
wherever Te varies along it, of the same order as the parallel sonic outflow at a strike point.

## Switching it on

```fortran
bcs(1)%floating_u = .t.      ! every connected wall segment: 1, 3, 4, 5, 9 on the usual grids
bcs(3)%floating_u = .t.
bcs(4)%floating_u = .t.
bcs(5)%floating_u = .t.
bcs(9)%floating_u = .t.
mach1_weak        = .t.      ! required, checked at setup
sheath_Lambda     = 3.d0     ! default
```

`dirichlet%u` stays `.true.` on those types (the floating row replaces it). `admissible_update`
is switched on automatically. Nothing else needs setting; `min_sheath_angle`, `D_perp_sc_num`
and every other parameter stay at their defaults. `floating_u_diag = .t.` prints the wall table.

## What is assembled

1. **Potential row** (`mod_boundary_conditions.f90`, `mod_floating_u.f90`): affine, no floor, no
   branch. `a_n = 2*e*F0*sqrt(mu0*rho0)/m_i`, `C_T = 2*Lambda/a_n` (halved in a single-T build),
   `C_V = sqrt(mu0*rho0)/F0`. Phi = +F0*u in JOREK's right-handed (R,Z,phi) basis; a_n carries
   the sign of F0 so the physical potential does not depend on the field direction. Self-test
   at setup.
2. **Weak Bohm condition on the total normal flow** (`mach1_weak`, `mod_boundary_matrix_open.f90`):
   one Galerkin residual per wall Gauss point on edges whose both endpoints are `mach1` types,
   `res = (B_pol.n)*Vpar - max(cs*|b.n| - vE.n, 0)`, weighted by `d(res)/d(Vpar) = B_pol.n`.
   The parallel flow supplies the outward normal flow the drift does not, is never asked to
   reverse, and loses authority as (B.n)^2 at grazing incidence with no threshold. The nodal
   Mach rows are not assembled, and type 3 gets no Dirichlet Vpar row. Columns on Vpar, Ti, Te
   and the u trace are exact; the |B| dependence on the free normal psi derivative is lagged.
3. **Inflow closure** on the density row where the total normal flow is inward:
   `-oint v*min(vn,0)*(rho - 0) dl`, exact columns on rho, u, Vpar. Zero where the flow is
   outward. The temperatures get no term.
4. **One total normal flow** `max(Vpar*(B_pol.n) + vE.n, 0)` in the sheath energy transmission
   and density reflection rows (exact u, Vpar, psi columns) and in the kinetic recycling flux.
5. **Admissible update** (`core/mod_jorek_timestepping.f90`, `core/mod_state_check.f90`): a step
   whose rho, Ti or Te is non-positive at any Gauss point is undone, dt is halved and the step
   re-solved (up to 8 halvings, then abort). dt recovers geometrically; `tstep_prev` is the step
   actually taken.
6. **Exterior sides only** (`mod_boundary_edges.f90`): the open-boundary integral is skipped on
   interior sides with two labelled endpoints; their number is printed once.

## Conventions

- Outward-positive normal flows. `vE.n = -orient*R*u_s/dl` with `orient = sign(y_s*n_R - x_s*n_Z)`
  and the edge tangent `(x_s, y_s)/dl`; `Vpar*(B_pol.n)` is a velocity since v = Vpar*B.
- `vpar0*ps0_s*normal_sign3 = Vpar*(B_pol.n)*R*dl` exactly; the total-flow measure reduces to it
  when u is constant along the wall.
- Branches (`max`, `min`) are decided on the n=0 state per Gauss point and differentiated exactly
  on the branch taken.

## Reading `floating_u_diag`

```
 [floating_u] type  inflow   max vE.n[m/s]  at (R,Z)   mom   min rho  at (R,Z)   min Te[eV]  at (R,Z)
```

`inflow` is the fraction of that type's wall length with inward total flow (where the inflow
closure is active). `mom = |sum Bn*res*dl| / sum |Bn|*cs*dl` is what the weak row imposes; the
pointwise residual at grazing incidence is not controlled by design and is not reported. The
minima are at the wall Gauss points; the volume minima are what `admissible_update` checks.

## Tests

`bash tests/floating_transport/run.sh` (gfortran, no MPI) compiles the production assembler
against fixture modules and checks: the normalisation for both field signs and a wall bias;
between-node undershoot detection; exterior-side classification; every column of the weak row,
the inflow term and the energy rows against central finite differences (each verified to fail
when a column is dropped); the saturated and closed branches; evenness, monotonicity and
vanishing of the row through B.n = 0; and that the diagnostics change no equation.

## Not covered

`n_order >= 5` trace DOFs beyond value and first derivative; non-axisymmetric evaluation of the
branches (per plane, on the FFT grid as the rest of the boundary assembly); boundary postproc
expressions along the wall; any cluster run. See the ladder in the plan.
