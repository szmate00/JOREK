# Floating-potential boundary condition (model600), nodal Mach-1 variant

Branch `floating-u-nodal` off `floating-u-clean-pr`: develop plus two things, nothing else. The particle side,
the boundary integrals (`mod_boundary_matrix_open.f90`) and the matrix construction are develop's.

## 1. `bcs(i)%floating_u`

On every u trace DOF (value and tangential derivative) of boundary type i the Dirichlet u row is replaced by

    Phi - V_wall = Lambda * k_B*Te / e        i.e.        u = C_T*Te + C_V*V_wall

with `sheath_Lambda` (default 3) and `sheath_V_wall` in volts (default 0): the zero-local-current limit of the
sheath characteristic. `C_T = 2*Lambda/a_n` (halved in a single-T build), `C_V = sqrt(mu0*rho0)/F0`,
`a_n = 2*e*F0*sqrt(mu0*rho0)/m_i`; the Te column and the RHS make the row exact in one solve
(`mod_boundary_conditions.f90`, `mod_floating_u.f90`, self-test at setup). `dirichlet%u` stays `.true.` on those
types. `Phi = +F0*u` (right-handed (R,Z,phi)); `a_n` carries the sign of F0.

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

## Not covered

`n_order >= 5` trace DOFs beyond value and first derivative; boundary postproc expressions along the wall.
