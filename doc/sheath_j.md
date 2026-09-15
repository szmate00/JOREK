# Sheath current boundary condition (model600), Option III

Branch `sheath-j-clean` off `floating-u-clean`. Plan and rationale: `doc/sheath_j_plan.md`.

## What it is

`bcs(i)%sheath_j = .true.` puts the sheath current-voltage characteristic

    zj = j_sat * ( 1 - exp(x) ) ,   x = Lambda - e*(Phi - V_wall)/(k_B*Te) ,   j_sat = c_sat*rho*(+-cs/|B|)

into the current-definition slot at the wall, and leaves the potential to the vorticity equation
(charge continuity). The wall then says: no perpendicular current into the wall, parallel current
equal to the sheath current, potential from continuity. The floating potential is its j -> 0 limit.

The wall is sorted by incidence with the same angle as the Bohm row:

| |b.n| >= sin(min_sheath_angle)             | |b.n| < sin(min_sheath_angle)              |
|---|---|
| u row: vorticity equation (retained)        | u row: floating potential (Dirichlet)      |
| zj row: weak sheath row (Gauss points)      | zj row: Dirichlet zj (as develop)          |

The electron current saturates where the plasma potential falls below the wall potential
(x >= Lambda): the exponent is capped there, f = 1 - e^Lambda. That is electron saturation, a physics
statement; it also bounds the u column. c_sat carries the sign of F0 like zj, so the current INTO the
wall, -zj*(B_pol.n)/F0, is independent of the field sign (checked in the harness for both signs).

## Switching it on

```fortran
bcs(1)%sheath_j = .t.        ! instead of bcs(i)%floating_u, never both on one type
bcs(3)%sheath_j = .t.
bcs(4)%sheath_j = .t.
bcs(5)%sheath_j = .t.
bcs(9)%sheath_j = .t.
mach1_weak            = .t.
mach1_weak_drift      = .t.  ! the D configuration of floating-u-clean
mach1_weak_drift_cut  = .t.
floating_u_diag       = .t.
```

`dirichlet%u` and `dirichlet%zj` stay `.true.` on those types (the condition takes their rows over).
`sheath_Lambda` (3) and `sheath_V_wall` (0) as for the floating potential.
`sheath_j_pin_current = .t.` is the stage-1 test: the u row is released but zj stays Dirichlet, so the
continuity-set wall potential is tested on its own against D.

## Reading the log

```
 [sheath_j]   type  e-sat    j/jsat min     max     Phi[V] min      max     Inet/Isat
```

per type, over the wall length carrying the row: fraction of that length on the electron-saturated
branch; min and max of j/j_sat (negative = electron current, +1 = ion saturation); min and max of
the wall potential in volts; and the net current into the wall as a fraction of the saturation current
integrated over the type (Inet = -sum zj*(B_pol.n)*R*dl, Isat = sum |j_sat*(B_pol.n)|*R*dl). The
`[floating_u]` table is printed as before.

## Tests

`tests/floating_transport/run.sh`, `test_sheath_j`: no zj row when off, pinned, or below the angle;
zj = 0 with u at the floating value is a root; every column (zj, u, rho, Ti, Te) matches finite
differences on the electron branch and on the saturated branch, where the u column vanishes; the
saturation current flows into the wall for both signs of F0; the diagnostics change no equation.
Dropping the u, Te or rho column makes the FD test fail.

## Not covered

The nodal side (`mod_boundary_conditions.f90`: node incidence, release of the u and zj rows) is not
compiled by the harness; `check_imports.py` covers its only-lists. No MPI build, no run.
