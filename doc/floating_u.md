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

## 4. Sheath-set wall flux on the total normal flow (automatic with `floating_u`)

Develop's rho/Ti/Te wall rows carry only the parallel flow, while the volume advection (not integrated by parts)
carries the full velocity, so with a wall potential that varies along the wall an ExB wall flux `q*vE.n` flows
implicitly in both directions with no inflow datum where `vE.n < 0`. Measured (2026-09-24, no `use_sc`): the run dies
where that inflow is largest (outer target, vE.n = -6e4 m/s, 11x the Bohm normal speed), density negative there
first, the Mach row still satisfied. On edges whose both endpoints are `floating_u` + `mach1` types the sheath rows
therefore impose the flux on the TOTAL normal flow (`mod_boundary_matrix_open.f90`, `sf_on`):

    Gamma = n*max(vn, v_fl),   vn = Vpar*(B_pol.n) + vE.n,   v_fl = factor*cs*|b.n|,   energy gamma_sh*T*Gamma

`v_fl` is exactly the parallel normal flow the nodal Mach-1 row imposes (`factor` = the `vpar_smoothing` weight, 1
without it), so with the row holding `Gamma = n*(v_fl + max(vE.n,0))`: unchanged where the ExB points out of the
plasma, the implicit emission replaced by the Bohm flux where it points in (the wall absorbs, never emits; SOLEDGE's
Bohm-Chodura inequality on the wall flux, SOLPS's `U_out` floor). In the rows: `max(vn, v_fl)*R*dl` replaces the
parallel measure and the excess sink `-q*max(v_fl - vn, 0)*R*dl` is added; boundary energy term
`oint q^2 (vn/2 - max(vn, v_fl)) <= 0` at every incidence angle; nothing divides by `b.n`. Exact columns for the
trial u_s, Vpar, psi_s (trace) and T (through cs); the |B| dependence on the normal psi DOF is lagged, as in develop.
The `c_angle` floors are untouched. Kinetic recycling (`mod_particle_wall_interaction.f90`) uses the same `Gamma` on
those edges, with the fluid's n, Te, Ti and the outward normal (`wall_normal_vector` points inward); the fluid ExB
velocity comes from `calc_EBpsiU`'s new optional `v_ExB = (-R*u_Z, R*u_R, 0)`. Off those edges, and in every run
without `floating_u`, fluid and kinetic expressions are develop's.

Serial checks (scratchpad harness, production routine + stub modules): residuals equal develop's for constant u and
flow above the floor; every u/Vpar/Ti/Te/rho column matches central FD (worst 1.5e-7) in the sub-Bohm, ExB-inflow,
ExB-outflow, near-branch-point and `vpar_smoothing` states; psi trace columns to the lagged-|B| level (1e-3);
with `floating_u` off the assembled edge matrix equals develop's to 3e-14 (re-associated products).

## Not covered

`n_order >= 5` trace DOFs beyond value and first derivative; boundary postproc expressions along the wall.
