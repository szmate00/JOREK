# Floating-potential boundary condition (model600), nodal Mach-1 variant

Branch `floating-u-nodal` off `floating-u-clean-pr`: develop plus two things, nothing else. The particle side,
the boundary integrals (`mod_boundary_matrix_open.f90`) and the matrix construction are develop's.

## 1. `bcs(i)%floating_u`

On every u trace DOF (value and tangential derivative) of boundary type i the Dirichlet u row is replaced by

    Phi - V_wall = Lambda * k_B*max(Te,T_min) / e        i.e.        u = C_T*max(Te,T_min) + C_V*V_wall

with `sheath_Lambda` (default 3) and `sheath_V_wall` in volts (default 0): the zero-local-current limit of the
sheath characteristic. `C_T = 2*Lambda/a_n` (halved in a single-T build), `C_V = sqrt(mu0*rho0)/F0`,
`a_n = 2*e*F0*sqrt(mu0*rho0)/m_i`; the Te column and the RHS make the row exact in one solve
(`mod_boundary_conditions.f90`, `mod_floating_u.f90`, self-test at setup). `dirichlet%u` stays `.true.` on those
types. `Phi = +F0*u` (right-handed (R,Z,phi)); `a_n` carries the sign of F0.

The temperature floor is the one the nodal Mach-1 rows use (`max(T, T_min)` on the axisymmetric node value).
Without it a wall node whose Te crosses zero imposes a negative potential, i.e. a reversed ExB drift, while cs at
the same node stays floored (measured 2026-09-06: type-1 min Te crossed zero and the Mach residual started
growing on the next output). The branch is taken on the node value; the Te column carries the exact derivative
of `max()` on that branch (1 or 0), so the relation is piecewise affine and one solve imposes it exactly while the
branch does not change. A clamped node has the constant potential `C_T*T_min + C_V*V_wall`: its derivative DOFs
and harmonics have target 0 and no Te column, so it drives no ExB flow. `T_min` is the namelist floor (0.1 eV in
the current runs); the `[floating_u]` residual column uses the same floor.

The wall potential then follows Te along the wall, so there is an ExB drift normal to the wall wherever Te varies
along it; the rest of the wall treatment (nodal Mach-1 row, sheath particle/energy fluxes, kinetic recycling) is
develop's, in which only the parallel flow enters the wall fluxes.

## 2. `mach1_omit_drift`

Develop's nodal Mach-1 row is

    Mach1BC = -Vpar + direction/Btot*factor*cs + factor/Btot*R^2*u_b/psi_b

The last term is the ExB drift correction (it divides by psi_b ~ B.n, and it exists in the value row only; the
slope row has no drift term at bicubic order). With `mach1_omit_drift = .true.` that term and its u column are
dropped (also the n_order >= 5 u_bb term), so the row is the marginal Bohm condition `Vpar = +-cs/|B|` everywhere.
Default `.false.` = develop.

## Namelist

```fortran
bcs(1)%floating_u = .t.      ! every connected wall segment: 1, 3, 4, 5, 9 on the usual grids
bcs(3)%floating_u = .t.
bcs(4)%floating_u = .t.
bcs(5)%floating_u = .t.
bcs(9)%floating_u = .t.
sheath_Lambda     = 3.d0
mach1_omit_drift  = .t.      ! .f. = develop's nodal row with its drift term
```

Everything else as in a develop run.

## 3. Wall diagnostics (`wall_diag`, default `.true.`, printed only when some `bcs%floating_u` is set)

`mod_wall_diag.f90`, sampled at the wall Gauss points of `mod_boundary_matrix_open` (first toroidal plane) on every
matrix construction (`wall_diag_every`, default 1). One line per boundary type, velocities in m/s, locations (R,Z):

    [floating_u] |u-uf|[V]: max floating-row residual at the wall nodes; vE.n out/in: largest outward and inward
                 ExB normal speed; min rho, min Te; inflow: fraction of the wall length with net inflow
                 (Vpar*B.n + vE.n < 0); exb>cs: fraction with |vE.n| > cs|b.n|
    [mach1]      |res| = |B.n*Vpar - cs|b.n|| min/mean/max (what the nodal row imposes under mach1_omit_drift);
                 max Mach |Vpar*B|/cs; max cs|b.n|; drift/cs = max |vE.n|/(cs|b.n|); Gauss points

`wall_diag_profile_every = N > 0` adds `[wall prof]`, the full wall profile (one line per Gauss point) every N
steps and whenever the wall minimum of rho or Te goes non-positive or halves. No equation is touched.

## Not covered

`n_order >= 5` trace DOFs beyond value and first derivative; boundary postproc expressions along the wall.
