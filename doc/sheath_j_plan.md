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

1. With j_sat sent to infinity (X = 1 - j/j_sat -> 1) the run reproduces D to solver tolerance.
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

| row (test fn)  | equation                          | D as it is                      | Option I (first)                       | Option II (derive in parallel)        |
|---|---|---|---|---|
| psi            | induction                          | replaced: Dirichlet psi         | RETAINED = wall Ohm's law, determines zj | retained                              |
| zj             | zj = Delta* psi (no surface term)  | replaced: Dirichlet zj          | replaced: Dirichlet psi (row swap)     | replaced: Dirichlet psi (row swap)    |
| u              | vorticity (charge continuity)      | replaced: floating row          | replaced: SHEATH row                   | RETAINED: charge continuity at the wall |
| w              | w = Delta* u (no surface term)     | replaced: Dirichlet w           | replaced: Dirichlet w                  | replaced: SHEATH row (row swap)       |

The row swap is JOREK's own way of freeing zj (and w): the Dirichlet condition of the potential-like
field goes into the row of its Laplacian, and the evolution equation is kept. With psi frozen at the
node, the retained induction row is stationary Ohm's law over the last element, eta*j = -grad_par(Phi)
plus the remaining induction terms; eta*zj enters undifferentiated, so that row is a complete equation
for the boundary current. No surface integral is needed. Its one dropped surface term is the numerical
hyper-resistivity (grad v . grad zj), small.

Option I is the smallest step from D: one row changes content (floating -> sheath), one Dirichlet flag
changes meaning (zj swapped, not pinned). The boundary pair Ohm + sheath is a nonlinear algebraic
system per node, solved implicitly each step; its Jacobian is well conditioned for finite eta, and in a
cold SOL eta is large so the current stays small there of its own accord.

Option II enforces charge continuity at the wall test functions and sets the potential through the w
slot, as the electrostatic edge codes do with the sheath current as the vorticity boundary flux. In
JOREK the vorticity row carries no surface flux; the coupling is through the value of zj. The retained
vorticity row contains int grad v . grad(du/dt) and the viscosity, both integrated by parts with dropped
surface terms, so its boundary row carries an implicit condition on normal derivatives. Derive before
trusting. Keep for the case where Option I shows a wall charge-balance violation (measurable: the
integral of the field-aligned normal current per type in the wall table).

## The sheath row

- **Weak, per Gauss point**, on the u trace, in the edge loop the Bohm row uses, weight
  d(res)/du = 1 (no fading: the potential must be set everywhere the sheath model applies). Corners
  accumulate instead of assign; the characteristic is evaluated per toroidal plane, i.e. on the real
  3D wall state, which the nodal n=0-only version could not do.
- **Residual in the potential**, the voltage defect: res = u - (2Te/a_n)*(Lambda - ln X),
  X = 1 - j/j_sat. Exact columns on u, Te, rho, zj (and Ti, Te through cs in j_sat). The current form
  asks for 1 - exp(x) steps (220 kV at 100 V off the root); measured on the old branch.
- **j_sat = c_sat * n * cs * |b.n|** written explicitly, cs from Ti+Te. With the marginal Bohm row
  there is no drift term to be inconsistent with, and the wall Vpar noise stays out of the potential.
  |b.n| carries the sign structure; the field-aligned normal current is j*(B.n)/(mu0*F0)-type, so the
  characteristic must be written for the current INTO the wall, both targets, both signs of F0
  (self-test, as for floating_u).
- **Saturation, u never without a row.** Where X <= 0 the characteristic has no voltage root. The row
  must still constrain u there. Candidate: the potential residual with X replaced by max(X, X_min)
  where X_min is not a tuning constant but the value at which the potential form's slope equals the
  current form's, i.e. continuity of the row; to be derived, then the branch is per Gauss point and
  continuous, like the Bohm inequality. The current limit j <= j_sat is the circuit's business (Ohm's
  law with rising Phi), not the u row's.
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
1. **Swap alone.** Implement the psi -> zj row swap (flag), keep the floating row. Run against D from
   the same restart. Expected: potential identical to D, wall current relaxes to its Ohm's-law value
   in the first steps; report min/max/integral of j per type in the wall table. This tests the freed
   current on its own.
2. **Sheath row** in the u slot (Option I), weak, potential form, with the grazing-angle sort. Harness:
   FD every column, both branches, both F0 signs, the j_sat -> infinity reduction to the floating row
   bit-for-bit.
3. **Ladder** on the cluster: (a) j_sat -> infinity = D; (b) current free, equal targets: zero net
   current per type; (c) unequal targets (the reference case has them): sign and magnitude of the
   thermoelectric current, potentials on opposite sides of floating; (d) hot phase and dt = 10 phase
   to 2000 steps.
4. **Option II** only if (b) or (c) show a wall charge-balance violation the u-slot placement cannot
   fix.
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
