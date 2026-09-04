# Plan: robust prescribed floating-potential BC in model600

Status: design only; no implementation in this document.

## 0. Clean-room rule

Do **not** merge or cherry-pick `feature/floating-potential-bc`. It is an experimental branch with
useful measurements but also superseded sign conventions, abandoned gates, relaxation variants,
Mach-1 experiments, gauge removal and competing row choices. Use its history only to construct
regression tests and to avoid repeating falsified ideas.

The implementation should be written afresh on the current branch and kept deliberately small:

* one boundary-mode flag;
* one audited physical-to-JOREK conversion shared with diagnostics;
* one linear trace constraint, `u-C_T*Te-C_V*V_wall=0`;
* one setup validator;
* one residual/normalisation diagnostic;
* focused tests.

No code from the old branch should enter merely because it already exists. Any line copied in
spirit must be justified from the governing equation and covered by a test. Experimental controls
such as `floating_u_relax`, `floating_min_bn`, `floating_gauge_removal`, `floating_amp_ramp`,
`floating_keep_u_row`, `floating_u_value_only`, `mach1_psib_floor` and `mach1_exb_term` are outside
the clean implementation.

Scope: impose the deliberately reduced sheath model

\[
V_p-V_{wall}=\Lambda T_e/e,\qquad \Lambda=3,
\]

on the plasma potential, without imposing the nonlinear sheath current--voltage characteristic.
This is the local zero-current (floating) limit of that characteristic. It is not a model of a
conducting wall carrying net current and it cannot reproduce thermoelectric target currents.

## 1. What has already been learned

This is not a new idea. Branch `feature/floating-potential-bc` implemented
`Phi=Lambda*Te/e+V_wall` as a coupled nodal Dirichlet condition on `u`.

Measured history:

* Value plus tangential-derivative trace rows survived about 718 steps without transport
  stabilisation. Value-only died in about 20 steps because the unconstrained Hermite tangent
  derivative generated large curvature and `w=Delta*u` read it directly.
* Freeing/failing to define `w` died in about 29 steps. A varying `u` trace exposes the fact that
  freezing `w` at its restart value is not a mathematically complete treatment of the mixed
  `u`--`w` system.
* Changing the grid, removing the constant gauge component, changing the boundary relaxation and
  several Mach-1 variants did not materially change the roughly 700-step failure.
* The first field to fail was `rho`, as a one-cell dispersive undershoot. The measured cell Peclet
  number was above the central-Galerkin limit. With the already existing shock-capturing transport,
  `use_sc=.true.` and `D_perp_sc_num=10`, the full `Lambda=3` case ran for more than 3500 steps.
  The remaining required check--whether this setting measurably smears the SOL solution--was not
  completed.
* Applying the same condition on every connected plasma-wall boundary type removed boundary-type
  jumps. Applying it on only one part of a continuous wall creates a one-element potential jump;
  since `v_E.n` is proportional to the tangential derivative of `u`, that is an artificial wall jet.
* The later weak-current campaign exposed failures that this reduced BC does not contain:
  `j_sat`, `B.n`, grazing incidence, characteristic saturation, row authority and inversion of the
  current characteristic all disappear when the prescribed condition is linear in `Te`.
* The weak-current campaign also exposed practices that must not be copied: evolving-solution
  gates, deleting a row without supplying a complete alternative equation, endpoint-OR edge
  activation, rank-local trace ownership, fitted determinant cutoffs and regularising the physics
  until a crash disappears.

The historical floating branch predates the sign correction described in
`doc/sheath_sign_review.md`. Its stability result remains useful because the magnitude was right,
but the physical potential and ExB direction in those runs were reversed.

## 2. Normalisation and sign, derived independently

Let

\[
n_0={\tt central\_density}\,10^{20},\qquad
m_i={\tt central\_mass}\,m_u,\qquad
\rho_0=n_0m_i.
\]

The model600/particle-field conversion is

\[
V_p[\mathrm V]=\frac{F_0 u}{\sqrt{\mu_0\rho_0}},
\qquad
k_BT_e[\mathrm J]=\frac{T_e^{JOREK}}{\mu_0n_0}.
\]

Therefore the prescribed physical condition

\[
V_p=V_{wall}+\Lambda k_BT_e/e
\]

becomes

\[
u_b=\frac{2\left(\Lambda T_e^{JOREK}+v_w\right)}{a_n},
\quad
a_n=\frac{2eF_0\sqrt{\mu_0\rho_0}}{m_i},
\quad
v_w=eV_{wall}\mu_0n_0.
\]

There is no additional minus sign. JOREK uses the right-handed `(R,Z,phi)` basis with clockwise
`phi`; the components in model600 implement the reference velocity ansatz in that basis. Thus
`a_n` has the sign of `F0`, and `u` automatically changes sign when `F0` does while the physical
floating potential remains positive relative to the wall.

For the supplied AUG case (`central_density=1.011088`, `central_mass=2.01410174369812`,
`F0=2.972306`):

```
rho0                         = 3.381578385970e-7 kg/m^3
sqrt(mu0*rho0)               = 6.518755039086e-7
physical volts per unit u    = 4.559622170458e6 V
a_n                          = 1.856385066569e2
2*Lambda/a_n, Lambda=3       = 3.232088055465e-2
Te(JOREK) corresponding 1 eV = 2.035678612360e-5
u corresponding to 1 eV Te  = 6.579492527774e-7
converted physical potential = 3.000000000000 V
```

This numerical identity should become a unit test. The current
`sheath-jsat-vpar-38ab278` branch still has the pre-correction negative `a_n`; that helper must not
be reused until the sign-fix commit is incorporated or its formula is corrected.

With separate temperatures, use `Te` directly. In the single-temperature build `Te=T/2`, hence

\[
u_b=\Lambda T/a_n+2v_w/a_n.
\]

## 3. Proposed discrete boundary condition

Introduce one explicit per-boundary-type mode, provisionally `bcs(i)%floating_u`, which is mutually
exclusive with `sheath_u`, `sheath_zj`, `sheath_zj_weak` and `natural%u`.

Define the lifting variable

\[
g=u-C_TT_e-C_VV_{wall},\qquad C_T=2\Lambda/a_n,
\]

and impose the homogeneous essential condition `g=0` on the boundary trace. For fixed
`Lambda=3` this relation is exactly linear, so it requires neither Newton iteration nor clipping,
gates, saturation slopes, incidence-angle floors or relaxation parameters.

### Trace degrees of freedom

Use exactly the same trace space as ordinary JOREK Dirichlet data:

* node value;
* every pure tangential derivative supported at the chosen polynomial order;
* no normal derivative and no mixed derivative.

For fixed `Lambda`, differentiation is exact degree-of-freedom by degree-of-freedom:

```
u_value       - C_T*Te_value       = C_V*V_wall
u_tangent     - C_T*Te_tangent     = 0
u_tangent_2   - C_T*Te_tangent_2   = 0       (when present)
```

This is why the condition is insensitive to the type-4/type-9 frame degeneracy: it never converts
logical derivatives into `(R,Z)`, divides by `xjac`, or evaluates `B.n`. The `u` and `Te` trace
DOFs share the same stored frame, so the geometric scaling cancels identically.

Do not use value-only imposition. It was measured to be much worse and it does not impose the
Dirichlet trace between nodes.

Implement this in a dedicated topology-driven trace pass, not inside the existing Mach-1 geometry
block. The historical branch made `floating_u` depend on `mach1` merely because that block already
computed boundary geometry; the new linear constraint needs none of that geometry. It must be
independent of `mach1`, `bc_natural_open`, side number, `B.n`, `xjac` and element orientation.

Build the set of constrained trace DOFs from incident boundary edges, deduplicate it, and write each
row once. At a smooth edge this gives value plus its tangent DOFs. At a corner with two incident
wall edges it supplies the trace constraints required by both edges without depending on which
element visits the node first. This is the mechanism that makes the same implementation valid for
types 1, 4, 5, 9 and other legitimate wall labels.

### Algebraic enforcement

First implementation should use the ordinary Dirichlet row-replacement path, keeping exactly the
same `u` diagonal scaling as the baseline and adding only the `Te` cross-column and the wall-voltage
RHS. This makes the matrix conditioning differ minimally from standard `u=0`.

The production-quality option is exact constraint elimination/static condensation of boundary `u`
trace increments in favour of `Te` increments. That removes the arbitrary `zbig` scale entirely,
but should be attempted only after the row-replacement version reproduces the historical result.

The constant wall voltage belongs only to the axisymmetric value equation. The proportional
`u=C_T*Te` relation applies mode by mode and derivative by derivative to all toroidal harmonics.

### Temperature positivity

The physical equation uses the evolved electron temperature. During a nonlinear excursion it must
not turn a negative numerical `Te` into a reversed sheath potential. Use the same smooth positive
temperature mapping as the model equations and include its exact derivative in the `Te` column.
Do not add a separate fitted sheath-temperature cap. In the physical positive-temperature regime
the map must be the identity to numerical precision.

## 4. Boundary topology

Apply one identical closure over every connected material-wall edge in scope. For the current AUG
wall this is expected to include types 1, 4, 5 and 9, but the implementation must determine and
report the incident physical edges rather than infer an edge from an OR of its endpoint labels.

At a junction:

* identical closures from both incident wall edges are harmless and must produce identical rows;
* a floating/non-floating junction is rejected unless it is explicitly declared as a physical
  change of wall model;
* stable edge and node identities are printed so MPI decomposition and element ordering can be
  checked.

The implementation must not use `boundary_index` for storage. Wall-grid types can be assigned after
that index is constructed; an earlier floating implementation indexed garbage for such nodes.

## 5. The coupled `u`--`w` and magnetic-boundary issue

Prescribing a spatially varying `u=3Te` is numerically simple but dynamically consequential:

\[
v_E\cdot n \propto R\,\partial_\ell u
                 = R C_T\,\partial_\ell T_e.
\]

The earlier runs proved that this real drift can expose two independent discretisation problems:

1. `w=Delta*u` is a mixed formulation. Freezing boundary `w` at a restart value is generally
   inconsistent with an evolving nonconstant `u` trace; simply releasing `w` without a complete
   replacement equation is worse.
2. A nonzero tangential gradient of `u` can advect poloidal flux against a frozen-`psi` boundary,
   producing a thin resistive current layer. This cannot be fixed honestly by forcing
   `[u,psi]=0`, because that would contradict `u=3Te` whenever `Te` is not a flux function.

The minimal reproducibility implementation should retain the baseline `dirichlet%w` and `psi`
conditions, because that is the combination that previously ran 3500+ steps with transport
stabilisation. It must be labelled as the prescribed-potential experiment, not yet as a complete
mathematical wall closure.

In parallel, measure rather than guess:

* residual of `w-Delta*u` on the boundary and first interior layer;
* current-sheet thickness relative to element size;
* boundary magnetic-flux work and `v_E.n`;
* fixed-location histories of `rho`, `w`, `zj`, `u-CTe` and solver residuals.

If these fail convergence under mesh refinement, the production implementation needs a consistent
`u`--`w` boundary pair (weak/Nitsche imposition or a discrete boundary-vorticity/influence method)
and a magnetic-flux boundary condition compatible with the allowed normal drift. Do not revive the
mesh-dependent `sheath_wall_vel` as the final answer.

## 6. Transport stability is part of the acceptance test

The previous full-amplitude run did not fail first in the boundary row. The imposed physical drift
pushed the density advection above the stable cell-Peclet range of the central scheme, producing a
one-cell negative-density dipole. The current user input already contains

```
use_sc        = .true.
D_perp_sc_num = 10
```

so the first reproduction should use that unchanged. This is evidence, not yet a universal default.
For an out-of-the-box implementation, compute and report the local cell Peclet number associated
with the imposed boundary drift. A long-term robust solution should use a consistency-preserving,
Peclet-based transport stabilisation whose strength follows from local speed and cell size, rather
than requiring a scan of `D_perp_sc_num` for every equilibrium.

Stability with `D_perp_sc_num=10` is not sufficient: density and temperature profiles at the outer
midplane, both targets and across the PFR must be compared with a resolved reference to demonstrate
that shock capturing is inactive away from under-resolved gradients and does not set the answer.

## 7. Implementation sequence and gates

### Phase A: clean reimplementation with the old result as an external A/B reference

1. Re-derive and write the fixed-`Lambda` coupled trace constraint on the current branch. Consult
   `feature/floating-potential-bc` only after the design and tests are specified, and do not merge
   or cherry-pick it.
2. Correct the `a_n` sign and share one audited normalisation helper.
3. Add strict mutual-exclusion and topology validation.
4. Run the supplied AUG restart with the existing shock-capturing settings and all connected wall
   types enabled.

Gate A passes when the run reproduces the historical lifetime without any floating-BC tuning and
the imposed trace residual is at roundoff.

### Phase B: verification independent of plasma evolution

1. Unit-test the SI/JOREK conversion numerically, including both signs of `F0` and both temperature
   builds.
2. On synthetic straight and curved meshes, impose analytic `Te` traces and verify `Phi=3Te` at
   nodes and boundary Gauss points.
3. Rotate an edge between local side families; the trace and integrated `E_t` must be invariant.
4. Repeat with different MPI decompositions and element orderings; rows and results must be
   bitwise identical where feasible.
5. Test `n_order=3` and `n_order=5`, and at least one non-axisymmetric harmonic.

### Phase C: coupled-physics diagnostics

1. Begin with `Lambda=0`; this must reproduce ordinary homogeneous `u` Dirichlet exactly.
2. Run `Lambda=3` from the same restart and track boundary work, cell Peclet number, magnetic flux
   advection and the `w=Delta*u` residual.
3. Demonstrate that the first unstable quantity, if any, is identified at a stable edge/node rather
   than inferred from global extrema.
4. Verify the physical sign directly: for positive `F0` and positive `Te`, the reconstructed plasma
   potential must be `+3Te[eV]` volts above a grounded wall, and reversing `F0` must reverse `u`
   but not the reconstructed volts.

### Phase D: remove remaining numerical dependencies

1. Replace or validate `zbig` enforcement by exact constraint elimination.
2. If boundary `w` is not convergent, implement and verify a consistent weak `u`--`w` boundary
   pair; do not merely unfreeze `w`.
3. If frozen `psi` produces an unresolved current layer, derive a compatible magnetic boundary
   treatment from the induction equation and wall model.
4. Make advection stabilisation respond automatically to local Peclet number, then show convergence
   as mesh resolution is increased and artificial diffusion decreases.

## 8. Acceptance criteria

The prescribed BC is ready only when all of the following hold:

* no `detmin`, `min_bn`, `ratio_max`, saturation slope, current clip or hand-tuned relaxation is
  needed;
* the trace identity is satisfied at boundary Gauss points, not only at nodes;
* results are invariant to boundary side family, boundary-type subdivision, element ordering and
  MPI decomposition;
* `u`/potential sign and magnitude pass the independent SI conversion test;
* no one-cell negative-density dipole, boundary checkerboard or unresolved current sheet develops;
* target/PFR profiles converge with mesh refinement and are insensitive to reducing numerical
  stabilisation;
* the documentation states honestly that `Phi=3Te/e` is a prescribed local-floating,
  zero-current approximation, not the full current-carrying sheath solution.

## 9. Minimal eventual input surface

The intended user interface should need only a per-type switch and physical data already present:

```
bcs(1)%floating_u = .true.
bcs(4)%floating_u = .true.
bcs(5)%floating_u = .true.
bcs(9)%floating_u = .true.

sheath_Lambda       = 3.0d0
sheath_Lambda_local = .false.
sheath_V_wall       = 0.0d0
```

The ordinary trace marker `dirichlet%u` may remain true internally, as on the historical branch,
but the user should not have to coordinate contradictory flags manually. Setup should configure or
validate the required row ownership explicitly. All sheath-current flags must be false for this
experiment.
