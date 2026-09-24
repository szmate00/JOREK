# Sheath current boundary condition (model600), nodal branch

Branch `floating-u-nodal-sheath-j` off `floating-u-nodal`: the current-row form that ran on `sheath-j-clean`
(types 1, 4, 5 with 3 and 9 floating, 2026-09-21), ported onto the nodal Mach-1 row (`mach1_omit_drift`), the
total-flow sheath fluxes and the T_min-floored floating row of `floating-u-nodal` (see `doc/floating_u.md`).

## The row

`bcs(i)%sheath_j = .true.` puts the sheath current-voltage characteristic in the **zj row** of the wall,

    zj = j_sat * ( 1 - exp(min(x, Lambda)) ),     x = Lambda - a_n*(u - C_V*V_wall) / (2*Te),

one residual per wall Gauss point, Zbig*dl weight, exact columns on zj, u, rho, Ti, Te (corr_neg Te and rho,
with their derivatives), on edges whose both endpoints are sheath types and where |b.n| >= sin(min_sheath_angle)
(`mod_boundary_matrix_open.f90`). The potential at those nodes then comes from the **vorticity equation**: the
Dirichlet u and zj rows are released per DOF (`mod_boundary_conditions.f90`): the value DOF if any incident wall
edge is a sheath edge, a derivative DOF only if the edge in its own direction is. Non-released u DOFs of a sheath
type (below the angle, or a derivative along a non-sheath edge) get the floating row `u = C_T*max(Te,T_min) +
C_V*V_wall`. psi and w stay Dirichlet; `dirichlet%u` and `dirichlet%zj` must stay `.true.` on sheath types. The
electron current is capped at x >= Lambda (potential not below the wall): there the u column vanishes.

`j_sat = c_sat * rho * (+-v_fl/|b.n|) / |B|`, `v_fl = factor*cs*|b.n|`: the **parallel Bohm flow the nodal Mach-1
row imposes** (`factor` = the `vpar_smoothing` weight, 1 without it). This is the form that ran; it differs from
the rho/Ti/Te wall flux `n*(v_fl + max(vE.n,0))` by the ExB part, which carries ions and electrons together and
adds no net current. The alternative j_sat on the total flow (SOLEDGE/SOLPS) divides the ExB part by b.n and is
left for a later test. `c_sat` carries the sign of F0 like zj, so the current into the wall, `-zj*(B_pol.n)/F0`,
is independent of the field sign. Normalisation in `mod_floating_u.f90` (`sheath_j_norm`, Artola eqs. 5, 7, 8).

`thermoelectric_ohm` (`thermoelectric_coef` = 0.71): the thermal force `-0.71*grad_par(Te)/e` in Ohm's law,
written as the tauIC electron-pressure term with Pe/rho -> 0.71*Te (needs tauIC /= 0). Cherry-picked from
`sheath-j-clean` unchanged.

## Ion-saturation slope and switch-on ramp (2026-09-24)

With the plain characteristic the row has no solution where the plasma delivers more than j_sat (f <= 1), and
the potential runs away there in one solve; every from-0 run died this way at the inner strike point.

`sheath_j_ion_slope` (s, default 0): f = 1 - e^y - s*min(y, 0), the finite slope of the ion-saturation branch
(sheath expansion, the ion branch of a Langmuir characteristic). A node pushed to j > j_sat then has a finite
potential, Phi = Te*(Lambda + (j/j_sat - 1)/s) above floating, and the row keeps its u column there.

`sheath_j_ramp_time`, `sheath_j_ramp_alpha0` (default 0 and 0.05): switch-on ramp of the characteristic,
y = x/alpha with alpha from alpha0 at t = 0 to 1 at the end of the timestep ramp (0), at the given time (> 0), or
no ramp (< 0). A small alpha is the sheath of a plasma at alpha*Te: a stiff I-V that pins the potential near
floating and passes any current, so the equilibrium's wall current, which is no sheath current, can be relaxed by
the induction equation while the sheath tightens, as JOREK ramps RMPs. The electron cap sits at y = Lambda. The end
state is the same with or without the ramp; on a restart past the ramp alpha = 1. `[sheath_j] ramp alpha` prints
the current value.

`sheath_flux_on_sheath_j` (default .t.): A/B switch, .f. gives develop's parallel-only wall rows and recycling on
the sheath types (the first sheath_j run's configuration).

## Namelist (the Sept-21 configuration on this branch)

```fortran
bcs(1)%sheath_j   = .t.        ! targets: sheath current
bcs(4)%sheath_j   = .t.
bcs(5)%sheath_j   = .t.
bcs(3)%floating_u = .t.        ! corners / tangential wall: floating, pinned current
bcs(9)%floating_u = .t.
sheath_j_current_row = .t.     ! the only form here (default .t.; kept so the old namelist reads)
sheath_Lambda     = 3.d0
mach1_omit_drift  = .t.
thermoelectric_ohm = .t.       ! optional, needs tauIC /= 0
sheath_j_ion_slope = 0.03d0    ! finite ion-saturation branch (0 = hard saturation, the Sept-21 form)
sheath_j_ramp_time = 0.d0      ! default: characteristic switched on over the timestep ramp
```

Not on this branch (the namelist read aborts on them): `mach1_weak*`, `floating_u_diag`, `sheath_j_float_u`,
`sheath_j_ion_slope`, `sheath_j_e_slope`, `sheath_j_ramp_time`. Diagnostics are on by default (`wall_diag`).

## Reading the log

    [sheath_j] type  e-sat  |j|>jsat   j/jsat min     max    max|j/jsat| at (R,Z)      Phi[V] min      max     Inet/Isat

per type carrying the row: fraction of its length at electron saturation; fraction with |j/j_sat| > 1 (the part
the characteristic cannot hold); min and max of j/j_sat (negative = electron current, +1 = ion saturation) and where
the largest |j/j_sat| sits; wall potential min/max in volts; net current into the wall over the saturation current
integrated over the type. The `[floating_u]` and `[mach1]` tables are printed as before.

## Tests

`tests/floating_u_nodal/run.sh`, `test_sheath_j`: no zj row with sheath_j off or below the angle, present above it
and no u row; every zj/u/rho/Ti/Te column of the zj rows vs central FD at floating, on the ion side, beyond
electron saturation (u column zero) and with `vpar_smoothing` (worst 1e-7); the ion saturation current flows into
the wall for both signs of F0; develop identity with everything off (3e-14). The nodal side (row release,
`node_incidence`) is not compiled by the harness.
