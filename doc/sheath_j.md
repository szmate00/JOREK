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
sheath_j_cancel_flux  = .t.  ! default: magnetisation, viscous and kinetic wall fluxes out of the wall balance (see below)
sheath_j_cancel_vis   = .f.  ! 2026-09-26 review configuration: magnetisation only ...
sheath_j_cancel_kin   = .f.
sheath_j_ion_slope    = 1.d0 ! ... with the tangent-continued characteristic: f = 1 - e^x - x (x < 0),
sheath_j_e_slope      = 20.085537d0 ! f = 1 - e^Lambda - e^Lambda (x - Lambda) (x > Lambda), Lambda = 3
```

`dirichlet%u` and `dirichlet%zj` stay `.true.` on the sheath types. `sheath_Lambda` (3) and
`sheath_V_wall` (0) as for the floating potential.

`sheath_j_float_u = .t.` keeps the floating row on u and only frees the current (definition + surface
term): the measurement of the wall current the plasma delivers under the floating potential, in the
`[sheath_j]` table. Run it first.

## Zero current on floating wall segments (`bcs(i)%zj_zero`, 2026-09-24)

The Dirichlet zj row alone keeps the wall current at its t = 0 value, the equilibrium's current density at the wall,
for the whole run. On the floating types (3, 9) that is a permanent current source next to the sheath types and
contradicts the floating potential imposed on the same nodes (the zero-current point of the characteristic); the
independent review (`doc/review_sheath_j_2026-09-24.md`) identified the frozen type-9 current at the inner strike
point as the source that every late runaway starts from, its ratio to j_sat growing as the target cools. With
`bcs(i)%zj_zero = .t.` the pinned zj trace DOFs of that type are driven to zero (value and tangential derivative,
all harmonics, exactly in one solve), as SOLEDGE does on every non-sheath boundary. Recommended on 3 and 9 (and 2);
on a sheath type it acts only on the DOFs that are not released (below the angle, or a derivative along a
non-sheath edge).

## The [wall J] diagnostic (2026-09-25)

`[wall J]` (with `floating_u_diag`): per boundary type, the implicit boundary currents of the released vorticity row
against the sheath's capacity, all per unit edge parameter as they enter the u row with the same test function
(s along the wall, n outward); capacity j_sat*|psi_s| (the wall flux of v*[psi,zj]):

| column | current | boundary part of the assembled term |
|---|---|---|
| mag | R^2 \|p_s\| | pressure bracket R^2 [v,p] |
| exb | rho R^4 \|w\| \|u_s\| | ExB advection of vorticity rho R^4 w [v,u] |
| vis | visco R \|d_n(R^2 w)\| dl (old setup: visco R \|d_n w\| dl) | -visco R grad v . grad(R^2 w) |
| kin | (v_E^2/2) \|d_s(R^2 rho)\| | -(v_E^2/2) [v, R^2 rho] |
| dia | 2 \|tauIC\| R^3 \|Pi_Z\| \|d_n u\| dl | -2 tauIC R^3 Pi_Z grad v . grad u |

Integrated ratios per type and the largest local ratio of each with its (R,Z). Anything of order 0.3 or more at a
strike point is a current the sheath is asked to pass that no sheath can (independent review,
`doc/review_sheath_j_2026-09-24.md`, item F4). Not measurable at the boundary (element-local data): the
ionisation / kinetic particle-source term R^3 S grad v . grad u and the kinetic pressure coupling; zero in these
runs: tg_num, Wdia, the toroidal viscous part. With `sheath_j_cancel_flux` the mag (total-derivative part), vis
and kin terms are cancelled in the row; the print still shows their size. The exb column had R^2 instead of R^4 before
2026-09-25 (harmless: w is Dirichlet zero on the wall, so the term is zero anyway).

A SOLPS-style linearisation slope in the Jacobian (`sheath_j_patankar`, floor max(e^x, 1) on df/dx, residual
unchanged) was tried on 2026-09-25 and removed again: it left the row unsatisfied over 17-21% of the inner target
(j/j_sat 1.6-2.4 where the residual allows at most 1) and the run died at 608 against 672 without it.

## Cancelled wall fluxes (`sheath_j_cancel_flux`, 2026-09-25)

The vorticity equation is assembled integrated by parts. While u was Dirichlet the boundary terms of that
integration never mattered; at a released sheath node they are in the row, as currents through the wall that no
sheath sets. `[wall J]` measured two of them at O(1-100) j_sat where the runs died:

- the **viscous vorticity flux** `visco*R*d_n(R^2 w)` (`visco_old_setup`: `visco*R*d_n w`): w is Dirichlet at the
  wall, frozen at its t = 0 value, while the interior w follows the potential, so `d_n w` at the wall is a grid
  mismatch, not a flux; 4 j_sat steady at the cold end of the inner target with a matching 5 j_sat electron
  current, 18 j_sat at the point that ran away (run of 2026-09-25, step 608);
- the **magnetisation current** `d_s(R^2 p)`, the total-derivative part of the pressure-bracket flux `R^2 d_s p`
  (`R^2 d_s p = d_s(R^2 p) - 2 R R_s p`): a target does not collect the diamagnetic current (Rozhansky et al.
  2001; SOLPS closes the target balance without it), only the grad-B part `2 R p n_Z` is real.

- the **kinetic-energy flux** `(v_E^2/2) d_s(R^2 rho)`, the wall part of the convective polarisation term
  `-(v_E^2/2)[v, R^2 rho]`: with the first two cancelled it led the collapse at the inner 4/9 junction (run of
  2026-09-26, step 372: 45 -> 93 -> 678 -> 9e3 j_sat locally, two prints ahead of everything else).

`sheath_j_cancel_flux = .t.` (default) adds all three back with the opposite sign as surface terms of the u row
(`sheath_j_cancel_vis`, `sheath_j_cancel_kin`, default .t., switch the viscous and kinetic parts off separately;
the 2026-09-26 review's configuration is magnetisation only),
with exact columns on w (trace and normal-derivative DOFs), rho, Ti, Te (including the viscosity's temperature
dependence) and u (through v_E^2), on every edge with a sheath node at either end, no angle gate, current-row form
only. What is left in the wall balance is the polarisation flux (the capacitor), the grad-B pressure flux, the
small diamagnetic part and the sheath current: the closure of the drift-fluid codes, no perpendicular current
through the target. Weakly this is the `d_n w = 0` (zero viscous current) condition with w itself still
Dirichlet; the cleaner end state, w released with its definition row plus a Neumann penalty on its normal
derivative, is the next step if this one moves the crash. Harness: value-DOF rows integrate to the analytic
viscous + magnetisation + kinetic sum for w = c Z, u = c_u R, rho = rho0 + rho1 R (old setup: `int R dR` for the
viscous part), every column matches FD, no row without a sheath node or with the flag off.

`[wall B]` (with `floating_u_diag`): the **signed** wall contents of the u row integrated per type, in units of the
type's saturation current and in the sign they have in the row's right-hand side (s = n rotated by +90 degrees,
`d_s psi = R B_pol.n`): `Ish = -zj d_s psi` (the sheath current; compare with `Inet/Isat` of `[sheath_j]`, same
number up to sign convention and an R weighting), `gradB = +2 R p n_Z` (what the cancellation leaves of the
pressure flux), `pol = +rho R^3 d_n(delta u)/tstep` (the capacitor, LAGGED: last increment rescaled to this step; the mass term is
assembled with a minus sign, so this is its content on the F side; the sign was inverted until 2026-09-26),
`dia`, then `rest = -(Ish + gradB + pol + dia)`, which is what the boundary cannot evaluate: the half-cell volume
parts of every term, essentially the interior parallel current delivered to the wall cell. (The kinetic recycling
source is not in the u row at all, so it is not a candidate.) The cancelled fluxes mag, vis, kin are printed signed
as well. A type whose `rest` is large and of the sign of the electron current is being fed from the interior, not
from any wall flux.

## Angle gate of the current row and the wall profile print (2026-09-27)

`sheath_j_min_angle` (degrees, default -1 = use `min_sheath_angle`): the incidence angle above which the current row
is carried at a Gauss point and the u/zj DOFs are released at a node (same number in `mod_boundary_matrix_open` and
`mod_boundary_conditions`). Below it a sheath-type node keeps the pinned floating potential and zero current. Unlike
`min_sheath_angle` it leaves the particle/heat flux floors (`cs*sin(min_sheath_angle)`) and the kinetic recycling
untouched. Purpose: the 4/9 corner cell, where |b.n| ~ 0.02 (1.1 deg) and the ExB inflow Lambda*dTe/ds/B is ~20x
the Bohm outflow cs|b.n|, collapses under the current row (every run, last at step 672 with the continued
characteristic); with the gate above the corner's incidence the row stays on the plates only where the field lines
end steeply enough to refill the cell. The price is the current on the grazing corner cells.

`floating_u_prof_every` (steps, default 0): with `floating_u_diag`, prints `[wall prof]`, one line per wall Gauss
point (all ranks), every N steps and on the step where the wall minimum of rho or Te halves or goes non-positive:
type, R, Z, rho, Ti, Te, Phi, Vpar*Bn, vE.n, cs|b.n|, vn, b.n, dTe/ds, dPhi/ds, j/jsat, x (the last two only where
the current row is carried), r = dPhi/ds / (Lambda dTe/ds) (1 where the potential gradient is the floating one) and
delta = (Phi - Lambda Te)/Te (0 at floating). Read b.n and dTe/ds at the corner from it before choosing the angle.

## Tangentially filtered wall Te (`sheath_Te_smooth`, 2026-09-27)

The wall potential follows Te (u = C_T Te on floating segments; the current row anchors u near C_T Te where the
current is small), so its tangential gradient is Lambda dTe/ds and the ExB normal flow it drives is Lambda dTe/ds/B.
At the inner 4/9 corner that was 1-4e4 m/s, ~20x the Bohm outflow cs|b.n|, in every run, from a Te step of tens of
eV across one or two elements, and the corner cell emptied under it. b.n is smooth there (the step-1 `[wall prof]`:
1.15 deg at the corner, 1.1-1.4 on the neighbours), so the angle gate cannot single the corner out.

`sheath_Te_smooth = N` (default 0): the Te the wall rows see is filtered along the wall, N passes of the 3-point
[1/4 1/2 1/4] average over the two wall neighbours of every wall node (all four DOF components, n = 1 harmonic; one
pass = one element of smoothing, the length scale is the grid). `mod_wall_smooth.f90`: the wall chain from the
exterior labelled sides of the local elements, allreduced; built once per matrix construction. It enters (i) the
floating target u = C_T Te_sm (mod_boundary_conditions; lagged, so no Te column there) and (ii) the characteristic's
x and j_sat's cs at the sheath Gauss points (mod_boundary_matrix_open, `tesm_g`; lagged, Te columns of the zj row
zero). The plasma's own Te, the Bohm row, the sheath particle/heat fluxes and the recycling use the raw Te. Harness:
on a one-element loop with a Te step the filtered values are 0.005/0.007, and the zj rows equal the rows assembled
with the raw Te set to the filtered value, with a zero Te column. In `[wall prof]`, r = dPhi/ds/(Lambda dTe/ds)
drops below 1 where the filter acts.

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

`sheath_j_ion_slope` (default 0): finite slope of the ion-saturation branch in the current-slot form,
f = 1 - exp(x) - s*x for x < 0. With a hard saturation the characteristic has no voltage root wherever
the plasma delivers j >= j_sat, and the vorticity row drives Phi to infinity there (measured, type 1
alone, no use_sc: Phi max 234 -> 683 -> 2600 V at the outer target over 178 steps while the current
stayed inside [-7, +1] j_sat). The old route's 3900-step run had sat_slope = 0.03; physically it is the
sheath-expansion slope of a Langmuir-probe I-V curve in ion saturation.

`sheath_j_e_slope` (default 0): the same beyond electron saturation, f = 1 - e^Lambda - s_e*(x - Lambda)
for x > Lambda. Measured (type 1, no use_sc, no slopes): the run converged from its 2600 V transient to
Phi 73-410 V and Inet/Isat +0.17 at 280, then at the outer target the plasma pushed j/j_sat = -26 against
the cap's -19, Phi fell below the wall with no root, the ExB reached 6e5 m/s and rho went negative there;
crash at 438. Both rails of a hard characteristic leave the potential unanchored wherever the plasma
demands more than the sheath can pass.

## Corner nodes: release per DOF (2026-09-18)

A wall edge carries the sheath row if both endpoints are sheath types and the incidence along it is
above the angle. At a node the VALUE DOF is released if any incident edge carries the row (the 4/9-corner
lesson: a per-visit decision let one edge pin what the other released). A DERIVATIVE DOF is released
only if an edge in its own direction carries the row; two collinear wall edges share the derivative DOF
and either may release it. At a target/flux-surface corner (type 3 next to type 2) the flux-surface
derivative therefore stays pinned to its floating value like the type-2 neighbours. Runs I and J
(type 3 released wholesale) went to -19 kV at that node in one step; this is the fix.
