# Sheath current boundary condition (model600)

Branch `sheath-j-clean` off `floating-u-clean`. History and the runs that led here: `doc/sheath_j_plan.md`.

## What it is

`bcs(i)%sheath_j = .true.` imposes the sheath current-voltage characteristic, solved for the potential,

    u = C_V*V_wall + (2Te/a_n) * ( Lambda - ln X ) ,    X = 1 - zj/j_sat ,    j_sat = c_sat*rho*(+-cs/|B|) ,

as a weak row on the u trace (one residual per wall Gauss point, weight one, exact columns on u, zj, rho,
Ti, Te), and lets the wall current be a current: the zj row keeps its definition zj = Delta*psi, completed
at the wall by the surface term  oint v (dpsi/dn)/R dl  that the volume form drops (which is why JOREK
pins zj otherwise). psi stays Dirichlet, so the induction row at the wall stays dropped exactly as in every
fixed-boundary run, and w stays Dirichlet. The wall current is then the current the interior induction
equation drives into the last element, on its own time scale; the sheath sets the potential from it.

X is bounded to [e^-Lambda, e^Lambda]: the upper bound is electron saturation (potential not below the
wall), the lower one says the characteristic is trusted up to twice the floating drop, a model statement
that keeps the row finite on the ion-saturated branch. At a bound the zj/rho/cs columns vanish and the
u column stays, so the potential is anchored everywhere. c_sat carries the sign of F0 like zj, so the
current into the wall, -zj*(B_pol.n)/F0 = e*n*cs*|b.n|, is independent of the field sign.

The wall is sorted by incidence with the angle the Bohm row uses:

| |b.n| >= sin(min_sheath_angle)                       | |b.n| < sin(min_sheath_angle)        |
|---|---|
| u rows: weak sheath row                               | u rows: floating potential (Dirichlet) |
| zj rows: current definition + surface term            | zj rows: Dirichlet zj (as develop)   |

The decision is made once per node, over both of its wall edges. The surface term is assembled on every
Gauss point of every edge with a sheath node at either end, with no angle gate: a released node needs it
over its whole support (measured: gating it per Gauss point left the row as dpsi/dn = 0 weakly over part
of the support and the wall current absorbed a flux of order R*B_t/h, thousands of j_sat, at the edges
of the sheath region).

## Why this form (the five runs of 2026-09-15)

- A wall potential set by charge continuity, the vorticity row, is not anchored where the density is
  low (rung 1, three variants; the "current row in the zj slot" variant): the potential must have its own
  row with itself on the diagonal.
- Freeing the current by swapping the psi Dirichlet condition into the zj row and keeping the induction
  row (wall Ohm's law) makes the wall current a slack variable for the last element's parallel field: it
  reported +-1000 j_sat under the floating potential, steady, dipolar. The current must keep its
  definition.
- Corners must be released per node, not per visit; the sheath is never imposed where the field grazes.

## Switching it on

```fortran
bcs(1)%sheath_j   = .t.
bcs(4)%sheath_j   = .t.
bcs(5)%sheath_j   = .t.
bcs(9)%sheath_j   = .t.
bcs(3)%floating_u = .t.      ! the target/outer-boundary corner stays floating with a pinned current
mach1_weak            = .t.
mach1_weak_drift      = .t.  ! the D configuration of floating-u-clean, unchanged at a restart
mach1_weak_drift_cut  = .t.
floating_u_diag       = .t.
```

`dirichlet%u` and `dirichlet%zj` stay `.true.` on the sheath types. `sheath_Lambda` (3) and
`sheath_V_wall` (0) as for the floating potential.

`sheath_j_float_u = .t.` keeps the floating row on u and only frees the current (definition + surface
term): the measurement of the wall current the plasma delivers under the floating potential, in the
`[sheath_j]` table. Run it first.

## Reading the log

```
 [sheath_j]   type  e-sat    j/jsat min     max    max|j/jsat| at (R,Z)      Phi[V] min      max     Inet/Isat
```

per type, over the wall length carrying the rows: fraction of that length with an X bound active; fraction
with |j/j_sat| > 1 (the part of the target the characteristic cannot hold); min and
max of j/j_sat (negative = electron current, +1 = ion saturation) and where the largest |j/j_sat| sits;
min and max of the wall potential in volts; net current into the wall over the saturation current
integrated over the type (Inet = -sum zj*(B_pol.n)*R*dl, Isat = sum |j_sat*(B_pol.n)|*R*dl).

## Tests

`tests/floating_transport/run.sh`, `test_sheath_j`: rows present/absent (off, grazing, on, float_u); the
surface term integrates to -b*ln2 for psi = a*R + b*Z on the fixture edge; every psi column of the zj
rows, trace and normal-derivative DOFs, matches finite differences (dropping the normal-derivative loop
fails it); the potential row has zj = 0 at floating u as a root, matches FD inside the bounds, and keeps
only its u column at a bound; the saturation current flows into the wall for both signs of F0; the
diagnostics change no equation. The nodal side (node incidence, row release) is not compiled by the
harness; `check_imports.py` covers its only-lists.

## Start-up ramp

The sheath row's current dependence is ramped: u = C_V*V_wall + (2Te/a_n)*(Lambda - alpha(t)*ln X), alpha
from 0 (the floating row, current free but not yet acted on) to 1 (the full characteristic). By default
alpha reaches 1 when the timestep ramp ends, i.e. when the last `tstep_n` phase begins; `sheath_j_ramp_time`
sets the time explicitly, a negative value disables the ramp. This is what lets a run take step 0 with the
sheath flags on: the equilibrium's residual wall current relaxes under a nearly floating potential while the
gain grows, as JOREK ramps RMPs.

## Initialisation (alternative)

The equilibrium wall current is not a sheath current: it runs +-1.5 j_sat and changes sign along the
target, and the potential row turns that pattern into 80 V next to 220 V within one step (measured: the
outer target collapsed in ten steps from the equilibrium). Start the sheath run from a state in which the
wall current has relaxed under the floating potential: run `sheath_j_float_u = .t.` first (current free by
its definition, potential floating), then restart from it with `sheath_j_float_u = .f.`. The `|j|>jsat`
column of the measurement run says when the target interior is within the characteristic.

## Current-slot form (`sheath_j_current_row = .t.`)

The structure of the old `sheath-jsat-vpar-38ab278` weak-trace route, which ran ~3900 steps to timeout
on boundary type 1 alone (converged, I_sheath = I_Ampere to four figures): the characteristic
zj = j_sat*(1 - exp(x)) IS the zj row at the wall (weight one, electron saturation at x >= Lambda), and u
is left to the vorticity equation. No surface term, no potential row. Restored 2026-09-16 for the
type-1-only comparison; the earlier runs of this structure here all included types 4 and 9, which the
old branch found to fail within 4-8 steps on their own, so it was never tested fairly.
