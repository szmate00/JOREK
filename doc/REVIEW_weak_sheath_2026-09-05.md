# Weak sheath-current review, 2026-09-05

Reviewed branch: `sheath-jsat-vpar-38ab278`, commit
`2cedbf2f5f6e1803167edf76fbfb6a70fbf2fed9`.
All source line references below refer to this commit, not the checked-out `floating-u` branch.
Sources were inspected using `git show`; no production code or branch was changed for this review.
The existing grid is a constraint. This is a code review and implementation proposal, not a claim
that a complete simulation has been fixed or validated.

## Findings that change the diagnosis

### 1. The branch still has the wrong potential and current-conversion signs

`models/model600/mod_sheath_bc.f90:210` sets `a_n` negative for positive F0.
`models/model600/mod_sheath_diag.f90:180` converts `zj` to outward current with a plus sign.
Both disagree with the implemented right-handed (R,Z,phi) convention:

```
particles/mod_fields.f90:581: v_R = -R*U_Z + Vpar*psi_Z/R
particles/mod_fields.f90:582: v_Z = +R*U_R - Vpar*psi_R/R
diagnostics/new_diag/mod_expression.f90:1354: Jtor = -zj0/BigR
```

The corresponding relations, with rho_ref = n_ref*m_i, are

```
Phi_SI = F0*u / sqrt(mu0*rho_ref)
a_n    = +2*e*F0*sqrt(mu0*rho_ref)/m_i
c_sat  = -a_n/2
J_n,field-aligned = -zj*(B.n)/(mu0*F0)
```

Independent component/curl derivation: `doc/sheath_sign_review.md` already identifies this bug
on the reviewed branch. The executable still has the old signs. The same conclusion follows
from the reduced velocity ansatz in Artola et al., equations (6)-(7):
<https://arxiv.org/pdf/2101.01755>.

Compiled reproduction using the branch's unchanged `mod_sheath_bc.f90`, parameter stubs, and the
user's F0, mass and reference density: a_n = -185.6385066569. At nominal floating, Te=1 eV,
Lambda=3, Vwall=0, the physical potential reconstructed independently is **-3 V**, not +3 V.

Fix normalization, initialization, and current/potential diagnostics together. Validate both signs
of F0 and B.n. The old type-1 lifetime is useful numerical history but its signed potential is
not a physical reference that a corrected implementation must reproduce.

### 2. The sheath ion flux and the actual Mach-1 condition are different

`mod_boundary_matrix_open.f90:428` applies the clipped tanh smoothing whenever `vpar_smoothing`
is true. `mod_boundary_conditions.f90:561` applies nodal smoothing only when `bn_1*bn_2<0`;
otherwise its factor is 1. The nodal branch also does not apply the same lower clip.

For an edge with the same incidence sign at both ends, and the user's smoothing coefficients:

| absolute b_n | Gauss-point sheath factor | nodal Mach-1 factor |
|---|---:|---:|
| 1.2e-6 | 2.05907e-6 | 1 |
| 1.7e-4 | 2.30782e-4 | 1 |
| 0.03 | 0.598441 | 1 |
| 0.05 | 0.948819 | 1 |

These are formula comparisons at specified incidence, not measurements of a supplied restart.
They falsify the code comment that both routes use EXACTLY the same Bohm velocity. This mismatch
is especially relevant to type 5, where the measured problem is grazing incidence.
The rounded coefficients also leave factor(0)=4.63476e-7 rather than exactly zero; multiplying
by sign(b_n) leaves a small jump even without the v_perp floor. A shared smoothing law should
enforce its intended oddness and g(0)=0 algebraically rather than through rounded coefficients.

There are additional differences:

* `mod_boundary_conditions.f90:646` includes
  `factor/Btot * R**2 * u_b/psi_b` in the imposed Vpar.
  The default sheath-current estimate does not include it.
* This drift correction divides by the tangential flux derivative, proportional to B.n.
  It can grow without bound on grazing edges, and `0/0` is possible at exact tangency.
  The sign-change-only smoothing does not protect every grazing edge.
* For cubic elements, `dMach1BC` at :650 omits the tangential derivative of this drift term;
  the `u_bb` term at :684 is only added for n_order>=5, and is still not the complete derivative.
* The nodal sound speed uses `max(T,T_min)` (:587), natural particle/heat terms use
  `corr_neg_temp1`, and the sheath uses `sheath_temp_floor`. They are different cold-state fluxes.

Normalization of the drift correction also needs an explicit manufactured test. From the
implemented component velocities, exactly

```
-v_E.n/(B.n) = R**2*u_b/psi_b.
```

If the intended condition is cancellation of ExB normal flow by parallel flow, the added Vpar
is this quantity. The extra `factor/Btot` in :646 prevents exact cancellation. This identifies
a discrepancy with that interpretation; it is not permission to change the Bohm criterion
without deciding which presheath entrance and flow condition the model represents.

Implement one shared, explicitly derived wall-flow evaluator for the Mach condition, sheath,
heat/particle closure and diagnostics. Either use a field-aligned entrance condition or derive a
drift-aware entrance condition and its Jacobian. Never obtain a universally applicable grazing
closure by dividing a prescribed normal flow by B.n.

### 3. `sheath_jsat_from_vpar` has a zero-velocity discontinuity and an electron-flux error

`mod_sheath_bc.f90:459` treats exactly zero Vpar as absent data and returns to the full Bohm
branch. Arbitrarily small nonzero Vpar uses the fractional floor. Compiled test with
`sheath_jsat_vpar_min=0.1`: **jsat(0)/jsat(1e-100)=10**. Zero velocity is a valid solution state.
Use the optional argument/model flag to represent availability; never use a numerical sentinel.

More fundamentally, :611 evaluates `zj_sh=zj_sat*f`. Changing the ion estimate therefore also
rescales the electron-current amplitude, with Lambda otherwise unchanged. At fixed density,
temperature and geometry, this incorrectly preserves the same floating drop when only the ion
collection changes.

Use independently defined ion and electron normal fluxes. On the electron-repelling branch:

```
chi = e*(Phi_plasma-Phi_wall)/(k_B*Te)
J_n,sheath = e * [Gamma_i,n - Gamma_e0,n*exp(-chi)]
```

`Gamma_e0,n` must come from the adopted electron collection model. Its dependence on incidence,
cross-field transport and magnetic presheath physics must be explicit. If the expression is
factored, Lambda_eff=log(Gamma_e0,n/Gamma_i,n); that factoring becomes unsuitable when the ion
flux vanishes. Do not fix the problem by dividing by a floored ion flux again.

### 4. The supposed cross-field flux floor is not a normal collection floor

`mod_sheath_bc.f90:498` uses

```
zj_sat = c_sat*rho*sign(B.n)*(abs(g)*cs + sheath_v_perp)/Btot.
```

Conversion back to normal current multiplies this by B.n. Consequently the alleged normal
cross-field contribution is proportional to `sheath_v_perp*abs(b_n)` and vanishes at tangency.
It also introduces a sign discontinuity in the constrained `zj` as B.n changes sign.

A real normal ion collection contribution belongs in Gamma_i,n, without another incidence
projection. Whether ambipolar diffusion or ExB contributes net current must be derived for BOTH
species; a particle flux is not automatically an uncompensated electrical current.

### 5. The regularizations alter the physical and discrete fixed points

`mod_sheath_bc.f90:562` adds `s*softplus(X_lim)`. At X=0 its leading contribution is s*log(2),
so the claim that the floating point is unchanged is false. Compiled tests, default limiter
(-3,0.5), no electron tail:

| evaluation | result |
|---|---:|
| j/jsat at raw X=0, s=0 | 0.00123707676 |
| j/jsat at raw X=0, s=0.3 | 0.20936696477 |
| actual zero-current raw X, s=0.3 | -0.17039431749 |

The softplus also allows unlimited super-saturation ion current. The electron tail at :588
changes electron saturation. Subtracting s*log(2) would restore a zero but would not repair
the excessive ion-current capacity.

`mod_boundary_matrix_open.f90:687` clips the residual BEFORE weak integration. In general
`integral N*r=0` and `integral N*clip(r)=0` have different roots. For example, quadrature residuals
(1,-2), weights (2,1), and clip(r)=r/(1+abs(r)) give raw moment 0 but clipped moment 1/3.
Thus `sheath_weak_rmax` is not purely a globalization mechanism that preserves the Galerkin
solution. The claimed current-step bound also does not bound the coupled potential increment.

Assemble the intended residual unchanged. Stabilize the update using a consistent nonlinear
solve with an accepted-step line search/trust region or timestep rejection, evaluated against
that residual. This is a numerical requirement to test, not evidence that absent Newton
iterations caused the documented crashes. Do not automatically implement a general timestepper
rewrite before testing the coupled boundary block and the existing step acceptance mechanism.

### 6. `zj=zj_sheath` has the wrong mathematical limit at tangency

The weak trace block at `mod_boundary_matrix_open.f90:925` retains its zj mass matrix even as
the incidence-dependent potential coupling vanishes. With no added floor and g->0 it asks for
`zj->0`; zero normal parallel current at tangency does NOT require zero tangential parallel current.
The ratio gate and determinant gate do not change this fact.

Multiplying the entire normalized row by B.n or a small weight does not fix it. It cancels where
uniform and degenerates at zero. Adding a floating-u penalty while retaining zj=zj_sheath pins
both potential and, near floating, tangential current; that is a different wall model.

The proposed production direction is an edge-based normal-current weak closure in the reduced
charge/vorticity balance, with a complete electromagnetic and flow boundary pair. Retain the
current definition on the appropriate unconstrained test space. The model's own
`models/mod_poloidal_currents.f90:398` reconstructs a pressure-inclusive wall current from
`zj*psi_s + R**2*p_s`; the sheath diagnostic currently tests only its field-aligned part.
Derive the boundary flux from the actual reduced operator, including pressure/polarization
contributions where appropriate. Do not copy this diagnostic estimate into the PDE blindly.

At exact tangency with no normal collection, the parallel sheath contribution supplies no
potential condition. The remaining wall/flow/induction problem must still be well posed.
Grounding the metal specifies Phi_wall, not Phi_plasma. A global reference is only a gauge
when it does not alter physically specified plasma-wall voltage differences.

Magnetic presheath theory has explicit incidence/order assumptions; extending the ordinary
Boltzmann-electron sheath relation to arbitrarily shallow angles is not automatically physical.
See Geraldini, Parra and Militello, especially the validity condition in their abstract:
<https://arxiv.org/abs/1907.09421>.

### 7. Trace assembly and boundary selection have concrete defects

* **Endpoint activation:** `mod_boundary_matrix_open.f90:228` ORs node flags over an entire edge.
  The emission guard at :1109 checks frozen u/frame but not that node's sheath flag. A mixed
  corner can activate its other incident edge. Type 2 is the user's artificial PFR boundary;
  type 3 must not turn it into a material sheath merely by sharing a node with type 1.
* **Topology:** `matrix/construct_matrix_mod.f90:102` starts from endpoint labels, followed by
  type-pair rules. Reconstruct actual exterior edges from connectivity on restart, then assign
  physical boundary roles to edges. Preserve corner one-sided traces and shared scalar values.
  The existing `grids/mod_boundary` boundary-element metadata is a candidate to reuse and audit.
* **Approximate replacement:** `mod_sheath_trace.f90:313` overwrites selected entries with zbig
  through `models/mod_assembly.f90:58`, but never clears the complete row. Unlisted volume
  entries survive. Their smallness is assumed; skew geometry can make this assumption weak.
  Clear the row or eliminate the constraint algebraically, and test independence from zbig.
* **Scaling:** dividing by diagonal D makes a unit diagonal, not an O(1) matrix. A Hermite block
  with scales [[h,h^2],[h^2,h^3]] becomes [[1,h],[1/h,1]]. Use both row and unknown scaling,
  or a locally orthonormal trace basis/Gram factorization, with consistent mapping to all
  coupled variables. This can be done without moving nodes or changing the represented fields.
* **Dynamic fallback:** :280 can reject a row based on solution weights and then freeze only
  its zj increment. A complete replacement boundary problem is required if a row is removed;
  this is distinct from the history showing that the ratio gate was not necessary for failure.
* **Construction lifetime:** `construct_matrix_mod.f90:476` resets the module accumulator only
  for the main matrix. Boundary accumulation and application also run during direct harmonic
  construction (`matrix/mod_direct_construction.f90:28`). Contributions can therefore mix
  between main and preconditioner builds. Give each construction its own reset/context. This is
  a conditional code defect, not an established trigger in the user's solver configuration.
* **Capacity/coverage:** fixed 8000-row/96-column storage and the two-function trace loops limit
  scalability and basis coverage. Size from the trace graph. n_tor>1 is explicitly rejected;
  higher-order/refined cases need correct trace and constraint transforms or explicit rejection.

### 8. MPI diagnostics are wrong in a different way from the module's warning

The warning at `mod_sheath_trace.f90:30` claims the owning rank can lack an adjacent edge.
For ordinary unconstrained rows under the reviewed distribution this is false:
`communication/distribute_nodes_elements.f90:99` puts EVERY element touching an owned DOF into
that rank's local list. Both incident edges are available. Do not add a blind accumulator
all-reduce; that would count already replicated contributions again.

However, `sheath_diag_add` and `sheath_geom_add` run on every such local element and their global
reports use MPI_SUM without a unique edge owner. Areas and integrated currents can be inflated
by partition overlap, with nonuniform multiplicity. The ratios can look convincing because
both numerator and denominator are duplicated. Some trace-row diagnostic statistics also count
unowned rows (`mod_sheath_trace.f90:381`).

Use a unique diagnostic owner per exterior edge and an owned-row filter for row reports.
Verify 1/2/4-rank equality of physical totals AND normalized algebraic rows. Hanging-node
constraints/refinement must be audited separately: trace accumulation bypasses the ordinary
element-matrix constraint transforms.

## What the previous runs do and do not establish

Read: `doc/HANDOFF_sheath_multitype.md`, `doc/sheath_multitype_ideas.md`,
`doc/sheath_sign_review.md`, and the superseded `doc/HANDOFF.md`; history was inspected after
the source/formulation review.

* The ~3900-step type-1 run ended by timeout. It demonstrates a stable run of that implementation,
  not correct physical signs, complete current balance, or validation on other boundaries.
* The 8-to-517 comparison changes initialization and saturation slope as well as detmin. It
  cannot measure a 65-fold benefit attributable to detmin alone.
* Some type-5 detmin comparisons also change slope. The near-equal crash times are useful
  evidence, but should not be described as controlled single-parameter experiments.
* A global argmax that moves between outputs is not a fixed-location trajectory.
* qjac<0.9 does not mean an element is unusable. qjac measures angular skew only; unequal metric
  lengths matter too. Finite small determinants imply poor conditioning, not a missing DOF.
* Boundary sign samples do not establish that every volume element is unfolded.
* Missing current-definition surface terms cannot be inferred for normal derivative test DOFs
  with zero boundary trace. At corners, check the trace on BOTH incident edges explicitly.
* Strong-form volume ExB transport already exists. A new particle/enthalpy surface sink can
  double-count it unless the volume form is changed consistently.

Section 7.8 records that the natural-u flux correction, with freed zj and repaired natural-zj,
already failed after roughly four steps. Rewriting zj_sheath*B.n as q_sheath is algebraically the
same term. The genuinely new proposal must include corrected physical signs/fluxes, consistent
Mach and transport closures, verified boundary energy terms, and controlled nonlinear updates.
The history warrants a small coupled benchmark before another full simulation, not abandonment
of weak fluxes and not an unchanged retry of the old flags.

## Five implementation packages, in order

1. **Audited physical wall law.** Correct both signs. Separate ion/electron collection and remove
   the zero-Vpar sentinel. Unify normalization and flow-state evaluation. Remove saturation tails
   from the intended physical residual. Tests: 1 eV -> +3 V in the simple floating case; F0 and
   incidence reversals; independent electron amplitude; true ion/electron limits; finite-difference
   Jacobians including branch transitions. Diagnostic reference must be independently reconstructed.
2. **Shared sheath/Bohm boundary block.** Unify smoothing and thermal-state conventions, resolve
   the drift correction's normalization and remove division by vanishing psi_b. Formulate the
   entrance condition as a weak coupled constraint with the correct grazing limit. Test planar
   oblique surfaces through tangency with a tangential potential gradient; measure the actual
   normal particle flux and Vpar, not just the prescribed target value.
3. **Exterior-edge topology and exact trace algebra on this grid.** Rebuild metadata from existing
   connectivity; give type-3/9 junctions explicit edge roles; keep artificial type 2 separate.
   Implement full row elimination/clearing and two-sided block scaling; isolate assembly contexts;
   fix diagnostic ownership. Tests: rotated side-2/side-3 equivalence, mixed corners, manufactured
   fields on skewed existing mappings, MPI/OpenMP and preconditioner-build invariance.
4. **Complete weak normal-current closure.** In a small coupled psi/zj/u/w test, derive all boundary
   terms and the energy balance, retain current definitions, and couple the forward normal wall
   current without imposing tangential j=0 at grazing incidence. Establish the remaining tangent
   wall flow/induction closure explicitly. Advance to the real grid only after discrete charge
   and energy balances and boundary-rank tests pass. This is the substantial formulation work.
5. **Reliable evolution and acceptance.** Use the true residual with controlled coupled updates,
   positive density/temperature step acceptance and suitable conservative transport stabilization.
   Initialize by a compatible boundary projection/interior lift; check w and existing current.
   Compare 1, 1+4, 1+5 and all material edges from the same positive restart. Require convergence
   under timestep/resolution changes, no deleted material patches, independent wall-current and
   power balance, and physical current/potential profiles. Bound numerical tolerances by error
   estimates; avoid equilibrium-specific detmin, v_perp, or saturation-slope fitting.

Packages 1-3 contain actionable source fixes on the current grid. Package 4 is needed for a
physical all-incidence model; successful execution of 1-3 alone would not validate the existing
grazing constraint. A genuinely singular volume mapping cannot be repaired by a sheath BC, but
the available diagnostics do not establish that this grid has such a mapping.

## Checks actually performed

The exact branch `mod_sheath_bc.f90` was compiled with gfortran, bounds checks and invalid/zero/
overflow traps, using small parameter modules in `/tmp/jorek-sheath-review.zHJmdW`.
`review.f90` reproduces the normalization, shifted floating root, zero-Vpar jump and signed
v_perp-floor results above. `derivatives.f90` compares all seven partial derivatives (u, rho,
Ti, Te, g, B, Vpar) to centered finite differences in smooth Mach and actual-Vpar regimes.
Maximum relative errors were 2.80e-10 in each case. This validates those local derivatives away
from switches; it does not validate the assembled boundary Jacobian or physical model.

No MPI JOREK build or cluster simulation was run. No production source changes, commits,
branch switches or pushes were made by this review.
