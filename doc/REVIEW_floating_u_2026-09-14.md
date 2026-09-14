# Floating potential and Mach boundary review, 2026-09-14

Reviewed `floating-u` at `6546f4e68`. Sources and tests were exported to
`/tmp/floating-u-review-6546f4e68`; the working checkout and its existing changes
were left alone. All source line references below refer to that exported revision.

## Assessment for the supplied configuration

The supplied configuration uses `mach1_weak=.true.`, `min_sheath_angle=0.3`,
Lambda=3, zero wall bias, and floating potential on types 1, 3, 4, 5 and 9.
Both floating diagnostics are enabled, with the probe at R=1.60, Z=-1.11.

**The floating-potential row can be satisfied very accurately while the coupled
temperature solution becomes negative.** It transfers the wall's Te gradients
directly into normal ExB transport; neither the weak Mach projection nor the
temperature time update preserves positivity. The heat boundary term responds
to this total normal transport and retains an additional minimum-angle loss.

This gives a credible mechanism for the inner-target failure: a steep Te profile
drives rapid outward wall transport and thermal depletion, a discrete update or
polynomial reconstruction undershoots, and the invalid temperature then changes
the imposed potential, velocities and transport coefficients. However, **the
first term causing the latest inner-strike-point undershoot is not established**:
the latest run's log, full input, restart, and first failing field were not
provided. Existing local logs are older runs with different Mach implementations.
They must not be represented as measurements of the supplied weak-Mach case.

The implementation also has specific defects and differences between formulations,
listed below. A claim in the old campaign document that the row's small residual
rules out implementation errors in the coupled system is not justified.

## What the supplied switches actually do

- `mach1_weak` replaces the nodal Mach rows with a boundary penalty integral.
- `mach1_omit_drift` does not remove drift from that weak residual. It controls
  the nodal path; changing it is not a weak-Mach drift-off experiment.
- `min_sheath_angle=0.3` sets an additive particle/heat collection speed
  `c_angle*cs`, with `c_angle=0.00523598776`. It does not regularize or cut off
  weak Mach enforcement. `mach1_drop_grazing` does not control the weak integral.
- The listed floating flags set the potential trace; they do not independently
  specify the Mach masks or natural rho/Ti/Te masks. The older full inputs set
  `bcs(2:3)%mach1` false and alter their thermal BCs; those settings need to be
  included when identifying which equations reach each corner.
- The fixed diagnostic probe is at the outer target. It does not provide an
  inner-target history, although the extrema diagnostics can locate minima.
- `floating_u_wall_flux` and `floating_u_mach_flux` were not included in the
  supplied fragment. They default false; analysis of the default heat/particle
  path below assumes they remain false in the complete input.

## Evidence for the inner-target transport mechanism

The supplied older postprocessor tables contain a strong outward jet near
R=1.256382, Z=-1.045352. At that same point:

| Existing file | Te [eV] | Parallel normal speed [m/s] | ExB normal speed [m/s] |
|---|---:|---:|---:|
| `doc/pp_out.txt` | 53.4174 | approximately 0 | 25,034.2 |
| `doc/pp_omit_out.txt` | 41.8102 | 2,564.9 | 23,289.4 |
| `doc/pp_omit_heatflux_out.txt` | 23.8842 | 2,271.2 | 15,539.5 |

At these points the tabulated `Phi-3*Te` error is only about 1e-10 V or smaller.
Thus an accurately imposed floating potential coexists with substantial wall
transport. File names alone do not certify the exact executable or all A/B inputs.
These snapshots establish an existing inner-target jet, not the location or cause
of the latest negative-temperature event.

In the production boundary assembler, with `a=|B.n|/|B|`, `C=cs*a`,
`vEn=v_ExB.n`, and `P=Vpar*(B.n)`, the additional electron heat term on a
floating edge is proportional to

```text
-(gamma_sheath_e - 1) * rho * Te * [max(P + vEn, 0) + cs*c_angle].
```

See `models/model600/mod_boundary_matrix_open.f90:433-445,548-550`.
The strong volume equation supplies advection and compression separately. On an
outgoing branch, large positive vEn increases both physical transport and the
sheath correction. Applying the weak Mach condition does not cap positive vEn.
For the electron transmission factor 5.5 in the older inputs, the converted
`gamma_sheath_e` is 3, so the explicit boundary correction coefficient is 2.
The ion coefficient has a different conversion and should not be interpreted as
the same cooling term; see `initialise_parameters.f90:286-301`.

At `a=1e-3`, the remaining minimum-angle speed is 5.24 times the projected sonic
speed C. At `a=1e-4`, it is 52.4 times C. These are illustrative local incidences,
not measured values at the latest failing point. The 0.3-degree input still
introduces a finite loss where the magnetic collection tends to zero.

## Findings in the implementation

### 1. P1 — The accepted update does not enforce temperature or density positivity

`core/mod_jorek_timestepping.f90:393-428` accepts `solver%step_success`, updates all
nodal values, runs a diagnostic, and advances time. There is no test that rho,
Ti or Te remains admissible, no retry based on that test, and no nonlinear
reconvergence of the new wall active sets on this path. A direct linear solve
can succeed even when its update gives a negative temperature.

The `corr_neg` functions make selected rate/transport evaluations usable at cold
or negative raw temperatures. They do not constrain the evolved temperature
field. The floating floor also does not serve as a temperature limiter.

The older run logs use `time_evol_scheme='Gears'` and increase dt from 2 to 10.
If those settings are retained, the history term creates an additional concrete
undershoot mechanism. From the production delta scaling at
`mod_elt_matrix_fft.f90:446-448`, zeta at `:286-289`, and mass/history terms at
`:1965-1966,3794`, a scalar cooling equation `dT/dt=-k*T` gives

```text
r = dt_new/dt_old; zeta = r/(1+r)
T_new = [(1+r)*T_now - r*r/(1+r)*T_previous]
        / [1+zeta+k*dt_new].
```

For r=5, the numerator is negative when `T_now/T_previous < 25/36`, even if
both history values are positive. Using exact cooling history 1 -> 0.5 and a
fivefold step increase gives **T_new=-0.22016445**, whereas the exact solution is
**+0.015625**. This is a reproduction of the update algebra, not a simulated
JOREK temperature history. It shows why an implicit scheme and a converged solve
do not establish positivity. At constant dt, BDF2 likewise lacks unconditional
positivity. The default Crank–Nicholson option has its own positivity restriction.

Repair requires an admissible nonlinear update and a temperature/density spatial
discretization with a controlled positivity property. A constant arbitrary
sheath damping coefficient is not a substitute. Timestep reduction or a first-order
implicit comparison can diagnose the mechanism but cannot certify the spatial
scheme or serve as a universal tuned fix.

### 2. P1 — The nodal and weak Mach formulations prescribe different flows

For no smoothing, define P, C and vEn as above. The implementations prescribe:

| Path | Prescribed P at finite incidence |
|---|---|
| Nodal Mach (`mod_boundary_conditions.f90:878-886`) | `C + max(-vEn,0)` |
| New `mach1_weak` (`mod_boundary_matrix_open.f90:459-467`) | `max(C-vEn,0)` |
| Legacy `mach_one_bnd_integral` (`:559`) | `C` |
| Experimental `floating_u_mach_flux` (`mod_floating_transport.f90:37-57`) | `C + min(max(-vEn,0),2*C)` |

The last two include smoothing factors when enabled. Only one path should be
interpreted as active for a given input.

For vEn=C/2, nodal Mach gives P=C while the new weak route gives P=C/2.
For vEn=2C, nodal Mach gives P=C while the new weak route gives **P=0**.
For vEn=-5C, both new nodal and weak routes demand P=6C, whereas the older
experimental clipped route gives P=3C.

The code comment that the weak branch “impose[s] nothing” once drift is sufficient
is incorrect: the RHS target becomes zero, but
`amat(var_vpar,var_vpar) ∝ (B.n)^2` remains. It penalizes Vpar toward zero.
It changes parallel momentum, shear and particle/energy fluxes on the outgoing
inner-target branch. This is a physical change as well as a discretization change.

Choose the intended Bohm/presheath condition, including the treatment of already
supersonic outgoing flow, then compare nodal and weak discretizations of that
same condition. The supplied SOLPS manual, pp. 82-85, explicitly distinguishes
marginal/non-marginal options and links momentum, density, energy and potential
BC choices; their names alone do not establish equivalence here.

### 3. P1 — Weak residuals do not ensure bounded or pointwise outward flow

The new penalty is `1e12 * Bn * (Bn*Vpar - target)`; its Vpar derivative contains
`1e12*Bn**2`. At exact tangency the weak term vanishes, but at small nonzero Bn
and finite inward drift its desired root still has

```text
Vpar = [cs*abs(Bn)/B - vEn]/Bn.
```

That root is unbounded as Bn -> 0 with vEn<0 fixed. Multiplying the equation by
Bn removes a division from assembly; it does not make the requested solution
bounded. Where the bulk equation wins, the requested normal flow is instead
not enforced. A small weighted residual can coexist with a physically large
local flow error in that transition.

There is a second limitation independent of tangency. I assembled the unchanged
production weak boundary block using the existing serial fixtures on a regular
element, positive Ti/Te, and an exactly floating linear potential trace. Solving
its four Vpar trace equations gave a relative algebraic residual 8.85e-17, but
the cubic trace had a negative parallel-normal flow between Gauss points, about
0.68% of the local projected sonic flux. The prescribed target was nonnegative.
This is the boundary-only projection limit; the full volume system was not solved.
It demonstrates that passing the moment test is not a positivity test.

The finite-incidence root and its domain of validity need to be distinguished
from a numerically finite assembled penalty. Tangential electric fields are
physically consequential at shallow incidence, as analyzed by
[Geraldini, Brunner and Parra (2024)](https://arxiv.org/html/2401.07385v2).

### 4. P1 when present in the mesh — Weak Mach can be applied to an interior edge

`mod_boundary_matrix_open.f90:180-187` computes whether a floating edge is exterior
but then sets `mw_on` from `mach1_weak` and endpoint Mach labels without using
that result. The caller's exterior-edge rejection and connectivity-based direction
selection at `matrix/construct_matrix_mod.f90:104-127` are enabled only for the
older `floating_u_mach_flux` / `floating_u_wall_flux` switches.

Consequently, under the supplied switches, an interior side whose two endpoints
carry boundary labels can receive the weak momentum BC. I reproduced this with
two conforming elements sharing a side: the production topology helper reports
`exterior=.false.`, yet the production assembler gives a nonzero weak Vpar block
(maximum absolute coefficient approximately 307009).

Restrict assembly to verified exterior material sides and use their actual side
orientation. Whether such a side exists near the user's inner strike point
requires the actual mesh; this review does not assume it does. The existing
single-element tests do not exercise the full caller's edge-selection logic.

### 5. P2 — Selecting weak Mach still evaluates the singular nodal expressions

`mod_boundary_conditions.f90:865` sets `m1_skip` for weak mode, but lines 867-868
divide by `ps0_b` before lines 873-874 overwrite the results. The cubic slope
block at lines 930-948 also evaluates coefficients containing `1/ps0_b**2`
without a skip guard. With `mach1_omit_drift=.false.`, this block is entered even
in weak mode. Legacy smoothing calculations are also performed before the skip.

Thus weak-mode production still executes code that can trigger floating-point
exceptions at exact tangency. The production weak-assembler tests do not compile
or execute this nodal routine. Skip the entire nodal evaluation before forming
these expressions. This is a separate crash route; it is not evidence that it
causes a temperature undershoot before a crash in the reported run.

### 6. P1 for the cold transition — The floating floor cannot stop the first undershoot

`mod_boundary_conditions.f90:557-574` chooses the floor branch from the old
axisymmetric nodal temperature. If Told>T_min but the accepted update produces
Tnew<T_min, that solve still imposes the unfloored relation on Tnew. A negative
potential is possible in that first invalid step; the next assembly changing
the branch is too late to prevent it. It is not an exactly imposed max relation
when an update crosses the branch.

Moreover, all endpoint values may remain positive while a cubic Hermite trace
is negative between them. For example, endpoints T=1,1 with derivatives -10,+10
give `T(s)=1-10*s+10*s*s` and T(0.5)=-1.5. The endpoint floor remains inactive,
and the floating relation transmits the interior temperature undershoot to u.
The derivative DOFs, not just the nodal values, must be considered in any
positivity treatment. This also explains why nodal minima alone cannot identify
the first invalid temperature in a high-order simulation.

### 7. P2 — The default closed-wall correction does not close incoming particle advection

With `floating_u_wall_flux` false, the natural density correction remains
reflection plus the minimum-angle loss at `mod_boundary_matrix_open.f90:538-539`.
For negative total vn the floating heat code cancels its additional sheath
collection, but the strong volume advection is still incoming. Kinetic recycling
now uses `max(vn,0)` in `particles/mod_particle_wall_interaction.f90:1783-1786`.
Those are not identical charged-particle balances on an incoming patch.

Exact pointwise weak-Mach roots would make total vn outgoing at finite incidence.
The finite penalty, polynomial projection, omitted corners and grazing transition
do not guarantee those pointwise roots. Therefore the incoming branch needs an
explicit coupled transport model. The opt-in absorbing-wall experiment addresses
part of this issue, but does not establish positivity or fix the other findings.

There is also a derivative bug on non-floating edges: `fu_ven_open` initializes
to zero, and the unconditional heat/Vpar Jacobian correction at `:734-739`
cancels the legacy parallel heat derivative even when `fu_edge` is false and
the RHS still contains that flux. It matters only on non-floating edges that
actually carry a natural thermal row. The older type-2/type-3 thermal overrides
disable those rows, so this is not asserted to explain the supplied case.
An additional production-assembler finite-difference test confirms the bug:
the Te/Vpar matrix entry is approximately zero (-3.45e-23) while the required
negative RHS derivative is 5.68727e-7. Reproduction: `nonfloating_heat.f90` in
the exported review directory.

## What has improved, and what the tests establish

- The positive-temperature affine floating normalization is consistent with the
  code's velocity convention. Its small nodal residual does not validate the
  surrounding transport, but there is no basis here to blame a simple sign error.
- The latest nodal cubic drift-slope u columns have been added. An old diagnosis
  that these columns are wholly absent is no longer true at this revision.
- The corrected-temperature slope is now included in the natural sound-speed
  derivative. The old cold-state Jacobian omission is fixed in that path.
- The current weak assembler's tested derivatives pass finite differences, and
  the boundary term itself vanishes continuously at zero B.n.
- All existing `tests/floating_transport/run.sh` tests passed with gfortran,
  bounds checking and invalid/division/overflow traps. They exercise local kernels,
  a serial boundary assembler and selected analytic consistency checks. They do
  not evolve the coupled temperature equations or test the actual nodal wrapper,
  full mesh selection, full MPI system, or global energy/particle conservation.
- Additional reproduction sources are in `/tmp/floating-u-review-6546f4e68/`:
  `weak_projection.f90`, `weak_projection_matrix.txt`, and `internal_edge.f90`.
  The suite output is `/tmp/floating-u-review-tests.log`.
- No full JOREK/MPI simulation was run. The working checkout has no Makefile.inc
  and no mpif90 was found. The user's current working source changes were preserved.

## What would identify the first cause in the actual run

Use the last positive state and the first negative update from the **same** run
and binary. Establish whether Ti, Te or rho becomes negative first, and whether
the first event is nodal or between nodes. Record both old and updated values at
that same location, not minima whose locations change.

At that element and its incident wall side, evaluate the actual assembled energy
budget: parallel and ExB advection, compression, parallel/perpendicular conduction,
sheath flux including the 0.3-degree term, radiation/ionization and kinetic sources,
and the mass/history terms. Also record dt, dt_prev, the timestep scheme, positive
field minima, the unweighted normal-flow residual, actual vn, B.n/|B|, and exterior
edge identity. The existing fixed probe at the outer strike point cannot provide
this diagnosis for the inner strike point; the older inner jet is near R=1.2564,
Z=-1.0454, but the failing run's own minimum location should select the probe.

Resolve the definite assembly defects, then compare equivalent nodal/weak physical
conditions. An admissible implicit update, positivity control for high-order
temperature/density fields, and measured particle/energy balance are the substantive
requirements for stable simulation without empirical sheath damping. The available
evidence does not support another isolated Lambda, angle or stabilizer scan as a
complete fix.
