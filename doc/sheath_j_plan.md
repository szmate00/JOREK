# Sheath current boundary condition on floating-u-clean: plan

Written 2026-09-15, after the floating-potential A/B/C/D/E campaign on `floating-u-clean`.
Predecessors: `doc/floating_u_clean_plan.md` (transport package, results), the old branches
`sheath-current-bc`, `u-sheath-jv-branch`, `sheath-potential-bc-artola` (derivation notes only; no code
is taken from them). Source of the characteristic: Artola, *Sheath boundary conditions for the electric
potential in JOREK* (2026-07-30), `potential_BC_JOREK_ionsat_v2.pdf`, eqs. (5), (6).

## Goal and definition of done

Replace the floating potential by the sheath current-voltage characteristic

    Phi - V_wall = (Te/e) * [ Lambda - ln(1 - j/j_sat) ] ,     j_sat = e n cs |b.n| ,

with the wall current free, so that a thermoelectric current can flow between targets of unequal
temperature. Out of the box in the sense established for floating-u: default namelist, the production
timestep ramp, no stabiliser, no floor, no Mach constant; the wall sorted by incidence, not by type
label.

Done means, on the reference case:

1. With the wall current held at zero the run reproduces D (floating) to solver tolerance.
2. With the current free and both targets at the same temperature, the net wall current per type is
   zero to the accuracy of the row, and the potential is the floating one.
3. With unequal target temperatures the thermoelectric current flows from the hotter to the colder
   target with the Harbour/Staebler-Hinton sign, and the wall potential on the two targets sits on
   opposite sides of floating, as the wall table shows.
4. Survives the hot phase (first 600 steps, 16 km/s drift) and the dt = 10 phase to 2000 steps.
5. Every new Jacobian column finite-difference checked; the harness fails when a column is dropped.

## What carries over unchanged (do not touch)

The transport package that made D run is independent of how the potential row is written:

- weak Bohm row, marginal below the grazing angle, SOLPS compensation above it (`mach1_weak`,
  `mach1_weak_drift`, `mach1_weak_drift_cut`);
- inflow closure on the density row where the total normal flow is inward (`mach1_weak_inflow`);
- one total normal flow in the sheath energy fluxes, density reflection and kinetic recycling;
- exterior sides only; wall diagnostics table.

Lessons from the campaign that bind the design below:

- one residual per Gauss point, weighted by its own sensitivity, no separate slope row (nodal
  value+slope pairs came apart twice; the weak Bohm row never did);
- a wall DOF is never left without a row (E: free Vpar on the grazing wall ran to Mach 2e4 in 23
  steps; the old sheath module's saturated branch left u without a row);
- one existing angle sorts the wall into "sheath model applies" and "floor model applies"
  (D: types 1, 3, 4, 5, 9 needed no individual treatment);
- the bound was what hurt in C, not the compensation: do not introduce a smooth cap whose scale
  varies along the wall inside a boundary row;
- the drift is 16 km/s in the hot phase and decays as the divertor cools; anything the current BC is
  for (asymmetry, HFSHD) lives in the first 600 steps.

## The wall system per boundary node

Four fields, four rows. Both definition rows are integrated by parts with the surface term dropped
(`mod_elt_matrix_fft.f90` current definition and vorticity definition), so at a boundary test
function they silently impose a zero normal derivative; that is why JOREK pins zj and w.

| row (test fn) | equation                          | D as it is                | Option III (FIRST)                        | Option I (alternative)                    | Option II (derive)                  |
|---|---|---|---|---|---|
| psi           | induction                          | replaced: Dirichlet psi   | replaced: Dirichlet psi (unchanged)       | RETAINED = wall Ohm's law, sets zj        | retained                            |
| zj            | zj = Delta* psi (no surface term)  | replaced: Dirichlet zj    | replaced: SHEATH row, zj = j_sat*f(u)     | replaced: Dirichlet psi (row swap)        | replaced: Dirichlet psi (row swap)  |
| u             | vorticity (charge continuity)      | replaced: floating row    | RETAINED: charge continuity sets u        | replaced: SHEATH row (potential form)     | RETAINED                            |
| w             | w = Delta* u (no surface term)     | replaced: Dirichlet w     | replaced: Dirichlet w (unchanged)         | replaced: Dirichlet w                     | replaced: SHEATH row (row swap)     |

**Option III (user's proposal, 2026-09-15, minus freeing the induction row):** the sheath in the
current-definition slot, the vorticity equation retained. This is the electrostatic edge codes' sheath
condition: nabla.j = 0 with the wall entering through the parallel current, whose wall value is the
sheath current. Properties:
- u never loses its row: on the ion-saturated branch f -> 1 and the u column of the sheath row
  vanishes, but u keeps the vorticity equation. The orphaned-DOF problem disappears structurally.
- the current form is the right form: res = zj - j_sat*f(u), d(res)/d(zj) = 1, no log, no X, no
  branch. The 220 kV objection applied to a row that solved for u.
- no grazing singularity in the row: the characteristic is for j.n = e*n*cs*|b.n|*f; in terms of
  the parallel current JOREK carries, |b.n| cancels, zj = c_sat*n*cs*f(u). Weight one.
- the retained vorticity row's dropped surface terms ARE the wanted physics: the polarisation term
  integrated by parts drops rho*dn(du/dt) (no perpendicular polarisation current into the wall), the
  viscosity drops dn(w). Wall statement: no perpendicular current, parallel current = sheath current,
  potential from continuity.
- psi stays Dirichlet and the induction row stays dropped at the wall exactly as today. Retaining
  it (psi free at the wall) would let B.n evolve: a leaky wall without a vacuum response. Not a BC.
- unlike Option I, the wall current does not depend on the last element's eta.
- electron branch: f = 1 - exp(Lambda - e*Phi/Te) is unbounded for Phi below the wall potential;
  physically the electron current saturates there, f(Phi_wall) = 1 - e^Lambda (about -19 for D).
  Capping the exponent at Lambda is electron saturation, a physics statement, not a clamp; it also
  removes the overflow.

**Option I:** JOREK's row swap frees zj: the psi Dirichlet condition goes into the zj row, the
induction equation is kept and, with psi frozen, reduces to stationary Ohm's law over the last
element, which sets the wall current; the sheath (potential form) sits in the u slot. Complete
without a surface integral (eta*zj enters undifferentiated; the one dropped term is the numerical
hyper-resistivity). Kept as the alternative if III shows the wall potential from continuity is ill
conditioned somewhere.

**Option II:** both swaps, sheath in the w slot, vorticity and induction retained. Derivation
exercise; not first.

## The sheath row

- **Weak, per Gauss point**, on the zj trace, in the edge loop the Bohm row uses, weight
  d(res)/d(zj) = 1. Corners accumulate instead of assign; the characteristic is evaluated per
  toroidal plane, i.e. on the real 3D wall state, which the nodal n=0-only version could not do.
- **Current form** (Option III): res = zj - c_sat*rho*cs*f(u), f = 1 - exp(min(x, Lambda)),
  x = Lambda - a_n*u/(2Te). Exact columns on zj, u, rho, Ti, Te (through cs and through x). The u
  column is -c_sat*rho*cs*exp(x)*a_n/(2Te) on the electron branch and vanishes on the saturated ion
  branch, where u is carried by the vorticity row alone.
- **Sign structure**: the characteristic is for the current INTO the wall on both targets and for
  both signs of F0. Anchor on Artola eqs. (2)-(4) and on mod_expression's Jpol; self-test as for
  floating_u. The sign of the outflow (direction of Vpar in the marginal Bohm row) is the sign of
  the ion saturation current.
- **Potential form** (Option I only): res = u - (2Te/a_n)*(Lambda - ln X), X = 1 - j/j_sat; needs a
  saturated-branch treatment that keeps a u row.
- **Below the grazing angle**: floating row and zj Dirichlet, exactly as D keeps the marginal Bohm
  row there. j_sat -> 0 at grazing incidence would otherwise put the whole grazing wall on the
  saturated branch, and no current flows into a wall the field does not reach.
- **The floating anchor.** With j_sat -> infinity the row is the floating row; that is the regression
  against D. With the current free, "floating" should mean zj = 0 at the wall (insulating), not zj at
  its restart value; decide and document.
- **Energy channel**: gamma_sheath_e stays fixed in the first version; with current the electron
  transmission depends on the drop (SOLPS B.52). Second order; document, measure, then decide.

## Stages

0. **Paper.** One page: the four rows per node under D, Option I, Option II, every retained row's
   dropped surface term marked. Sign of the field-aligned normal current for both targets and both F0
   signs, anchored on Artola eqs. (2)-(4) and on mod_expression's Jpol.
1. **Vorticity row retained alone.** Release the u Dirichlet/floating row above the grazing angle
   with zj still pinned, i.e. Option III with the sheath row replaced by Dirichlet zj. This tests the
   retained vorticity row as an equation for the wall potential on its own, against D from the same
   restart. Expected: a potential that is no longer Lambda*Te pointwise; report Phi min/max per type.
   If this is ill conditioned anywhere, Option I becomes first.
2. **Sheath row** in the zj slot (Option III), weak, current form with electron saturation, with the
   grazing-angle sort. Harness: FD every column, both branches (electron, ion-saturated), both F0
   signs, the reduction to the floating potential when the sheath row is made stiff (j_sat -> 0 with
   zj = 0: x -> 0 exactly, Phi = Lambda*Te), bit-for-bit against D where applicable.
3. **Ladder** on the cluster: (a) stiff sheath (zj = 0) = D; (b) current free, equal targets: zero net
   current per type; (c) unequal targets (the reference case has them): sign and magnitude of the
   thermoelectric current, potentials on opposite sides of floating; (d) hot phase and dt = 10 phase
   to 2000 steps.
4. **Option I** if stage 1 shows the continuity-set potential ill conditioned; **Option II** as a
   derivation exercise only.
5. Physics refinements, each measured before being built: drop-dependent gamma_e; Lambda(Ti/Te);
   secondary electron emission; the diamagnetic part of the model's own normal current.

## Diagnostics to add to the wall table

Per type: j/j_sat min and max, saturated fraction of the wall length, voltage defect max with
position, integral of the field-aligned normal current (the net current to the vessel per type),
Phi min and max in volts. Same Gauss-point accumulator as the existing table.

## Open physics, stated so it is not mistaken for numerics

- Whether the sheath collects the full sonic ion flux where the tangential field drives the drift
  away from the wall (marginal vs non-marginal bracket; D is the non-marginal form).
- The grazing wall carries no current in this model by construction (floor model). Real walls at
  grazing incidence do collect current through the magnetic presheath; out of scope.
- The vessel is grounded: net current to the wall is allowed. A per-tile floating wall (zero net
  current per segment) is a later scalar constraint, not part of this plan.

### Rung-1 result and Option I (2026-09-15)

Rung 1 (u released, zj pinned) died within 14-19 steps in three configurations (equilibrium start, D@200
start, drift and marginal rows), always at a type-4/9 corner, the potential kink growing x2-3 per step at
dt = 6 ns while the Bohm row held at Mach 1.00: the continuity-set potential has no anchor without the
sheath conductance, and that conductance vanishes on the saturated branch and at low density. III is
therefore only anchored where the plasma is dense and near the characteristic. Option I coded under
`sheath_j_ohm` (row swap + weak potential row on u, bounded X); III kept for the rung-2a comparison
(D@200 start, current free, marginal row) which has not been run.

### Final formulation (2026-09-15, after the Option I measurement)

The wall-Ohm's-law variant reported +-1000 j_sat, steady and dipolar, under the floating potential: with
the current definition swapped out, the wall zj was the only unknown in the induction row and absorbed the
last element's parallel-field mismatch at 1/eta. Kept instead: psi Dirichlet in its own row (induction
row dropped at the wall as always), zj = Delta*psi at the wall with the surface term oint v (dpsi/dn)/R dl
restored, the weak potential-form sheath row on u. Options III and the row swap removed from the branch.
Next: the measurement `sheath_j_float_u` (floating u, current by definition), then the sheath row.
