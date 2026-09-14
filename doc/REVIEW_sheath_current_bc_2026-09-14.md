# Review of `sheath-current-bc`, 2026-09-14

Reviewed implementation commit: `fda378495`, relative to parent `46882606a` / `develop`.
During the review HEAD advanced to `4e0c518c3`; its only change adds the sheath flag
to the boundary-condition documentation table. It was inspected and does not
change the source findings or reproduction results.
The branch was already checked out. No production source, branch, input, or existing
untracked file was changed. This report concerns the current implementation, not the
older weak-sheath branches discussed in other local review documents.

**Assessment: this branch does not yet meet the requested physical correctness,
grazing-angle robustness, or out-of-box stability requirements.** The local
constant-geometry, positive-temperature, axisymmetric current-law linearization
is sound. Its integration with the boundary system and its broader physical claims
have significant problems. Some are introduced here; others are existing conditions
that the new closure relies on and therefore must address to meet the requested goal.

## Findings

### 1. P1 — The documented activation leaves the boundary current fixed

Locations: `models/model600/mod_boundary_conditions.f90:329-353`,
`models/preset_parameters.f90:376-378`, `doc/sheath_current_bc.md:11-21`.

The default is `bcs(:)%dirichlet%zj = .true.`. Enabling `sheath_j` suppresses only
the Dirichlet **u** rows. On fixed-boundary harmonics, the old current rows still
receive the large Dirichlet diagonal, effectively freezing the boundary current
and its tangential derivatives. The special exemption for STARWALL-supported
free-boundary harmonics does not fix the default fixed-boundary case.

Consequently, starting with zero target current yields an approximately floating
potential while the target current remains effectively zero. Unequal target
temperatures cannot freely establish the intended SOL current through those rows.
Starting with a prescribed current beyond the instantaneous saturation current
can instead leave the sheath equation with no finite potential solution.

The documented flag recipe is therefore insufficient. Define and implement the
complete compatible psi/zj/u/w boundary system, including current-definition
surface terms when freeing current. Simply telling users to unset `dirichlet%zj`
would not by itself demonstrate that the resulting system is correctly closed.

### 2. P1 — The existing Bohm coupling still has an exact grazing singularity

Locations: new sheath coupling at `mod_boundary_conditions.f90:583-595`;
existing Bohm expressions at `:675-678`; `mod_sheath_current.f90:134-137,169`.

The sheath module eliminates Vpar using `direction*cs/|B|`, claiming that this is
what the Mach condition enforces. The actual nodal condition also contains

```text
factor / Btot * R**2 * u_b / psi_b
```

and its u Jacobian divides by `psi_b` as well. This tangential flux derivative
vanishes at magnetic tangency. The default smoothing flag is false. At exact
tangency, even `u_b=0` evaluates `0/0`; nonzero `u_b` produces an unbounded flow
as incidence decreases. A floating potential proportional to a varying Te provides
exactly such a nonzero tangential u gradient.

Compiled reproduction of the unchanged Mach expression with invalid/division
traps terminates at `psi_b=0`, while the same expression at `psi_b=1` runs.
For a positive 10 eV state, R=2, |B|=2, F0=3, and `u_b=u_float`, decreasing
`psi_b` from 1e-2 to 1e-8 raises the imposed `Vpar*B/cs` from 1.10 to about
99,749. These are manufactured local states, not measurements from a restart.

There is also a normalization issue in the drift cancellation. From the actual
velocity components in `particles/mod_fields.f90:581-582`,
`-v_u.n/(B.n) = R**2*u_b/psi_b`. If cancellation of the normal u-driven flow is
the intention, the additional `factor/Btot` in the imposed Vpar prevents it in
general. The shear/current module and natural heat fluxes consequently cannot
assume they share one Bohm particle flux.

The required repair is a consistently derived normal-flux/sheath/Bohm closure
with a defined tangency limit. Bounding the coefficients of the new u row does
not repair divisions in the coupled Vpar row. A finite-incidence derivation does
not automatically extend to exact tangency by cancelling B.n.

### 3. P1 — The tangential derivative row rejects exact solutions when |B| varies

Location: `models/model600/mod_boundary_conditions.f90:587-595`.

Let `G = j - jsat*f`, with `jsat ∝ rho*sqrt(Ti+Te)/B`. Away from an incidence
sign change, differentiating along the boundary includes the geometric term

```text
G_b = G_u*u_b + G_Te*Te_b + G_Ti*Ti_b + G_rho*rho_b + G_j*j_b
      + jsat*f*B_b/B.
```

The assembled derivative residual includes only the first five terms. In its
raw scaled form the missing term is `-xi*f*B_b/B`. Holding geometry fixed during
a Newton update does not mean it is spatially constant along the wall.

Counterexample: keep u, rho, Ti and Te constant; vary B along the boundary; set
`j(b)=jsat(b)*f`. The value condition is satisfied everywhere, and the exact
derivative is zero. Nevertheless, the implemented derivative residual is
`+xi*f*B_b/B`. The compiled production row reproduces this error: for x=1,
10 eV and B_b/B=0.1, its raw residual is 3.700135409e-7 in u units, equivalent
to 1.71828 V per chosen boundary coordinate. It drives an artificial potential
gradient or distorts the current derivatives of a valid sheath profile.

Derive the derivative residual from the complete spatial characteristic, then
linearize that residual. The current derivative Jacobian also omits the
value-column contributions from differentiating its state-dependent coefficients;
that is an approximate iteration even at constant geometry, not the exact
Jacobian asserted by the comments.

### 4. P1 — Row rescaling does not provide nonlinear stability near saturation

Location: `models/model600/mod_sheath_current.f90:172-187,196-211`.

Multiplying an equation and its Jacobian by the same frozen factor leaves its
Newton update unchanged in exact arithmetic. It improves representation and may
affect linear conditioning; it is not damping or a bound on the potential update.

For fixed positive density/temperatures, fixed j=0, and
`x = e*Phi/Te - Lambda`, the returned row gives exactly

```text
delta_x = 1 - exp(x).
```

The compiled production routine at Te=10 eV and x=10 requests
`delta_Phi = -220254.6579 V`, despite being only 100 V above the floating root.
The example lies well inside the exponent guard. It is especially relevant when
other boundary rows hold j fixed or when a temperature drop moves a previously
reasonable potential into saturation. The finite-row self-test cannot detect it.

Use a nonlinear solution strategy that checks the actual sheath residual and
physical admissibility, with compatible initialization and a treatment of the
saturation limit. This can be designed without equilibrium-specific physical
stabilizers; the present normalization alone is insufficient. No full-simulation
crash is claimed from this local counterexample.

### 5. P1 — The current closure omits the pressure contribution in the model's current reconstruction

Location: `models/model600/mod_sheath_current.f90:166-169`.
Independent implementation reference: `diagnostics/new_diag/mod_expression.f90:1325-1326`.

The derivation supplied in `potential_BC_JOREK_ionsat_v2.pdf`, equations (2)-(4),
uses `Jpol = -j*Bpol/F0`. The model600 diagnostic reconstructs instead

```text
JpolR = (-zj*BR - R*P_Z)/F0
JpolZ = (-zj*BZ + R*P_R)/F0.
```

Thus its outward current contains
`R*(-P_Z*n_R + P_R*n_Z)/F0` in addition to the field-aligned contribution.
At a tangential pressure gradient, setting `zj=0` and `Phi=Lambda*Te/e` satisfies
the new zero-current characteristic, but generally does not give zero normal
current according to the existing model reconstruction. This contribution is
particularly consequential at small B.n, where the field-aligned contribution
shrinks. Temperature and density gradients are explicitly part of the requested
operating regime.

Establish which current and which presheath entrance the closure represents, and
derive its connection to the reduced momentum/current equations. Either include
the required pressure/drift contributions consistently or explicitly restrict
the model's physical domain. This finding concerns the model's reconstructed
current; it does not assert that the reduced magnetic ansatz resolves every
physical current component independently.

### 6. P1 for nonlinear 3D runs — Only the axisymmetric current law is evaluated

Locations: `models/model600/mod_boundary_conditions.f90:563-568,646-662`.

Every state sample comes from harmonic 1. The same background coefficients are
placed on each harmonic's diagonal block, and all non-axisymmetric residuals
are overwritten with zero. This constrains increments through a homogeneous
background linearization; it does not enforce the nonlinear characteristic on
the existing 3D state or provide its Fourier-mode coupling.

For constant density and temperature, take the entirely electron-repelling state
`e*Phi/Te = Lambda + cos(phi)` with Lambda=3. The true mean characteristic is
`<j/jsat> = 1 - I0(1) = -0.266065878`. The background-only closure evaluates it
as zero. The mean rectification and generated harmonics are missing even though
the potential remains positive throughout the toroidal period. A restart with
an existing nonzero-harmonic BC error also receives no residual correction.

Evaluate the full wall state in toroidal quadrature and assemble the transformed
residual/Jacobian, or clearly limit this implementation to axisymmetry and the
appropriate small-perturbation regime.

### 7. P2 — A fixed Lambda gives the wrong electron flux dependence in Ti/Te builds

Location: `models/model600/mod_sheath_current.f90:169-185`.

The implemented unretarded electron amplitude is `jsat*exp(Lambda)`, hence it
scales as `sqrt(Ti+Te)` at fixed density and geometry. Maxwellian electron
collection instead scales as `sqrt(Te)` at fixed Te and density. Changing Ti/Te
from 1 to 9 therefore incorrectly multiplies the electron amplitude by sqrt(5),
or 2.236, while leaving the floating potential fixed.

This is an approximation even though the code differentiates its chosen formula
correctly. Consistency with the ion sound speed used here would require
`Lambda_eff = log(sqrt(mi*Te/(2*pi*me*GAMMA*(Ti+Te))))` under the simple
field-aligned Maxwellian assumptions, including its derivatives. Equivalently,
calculate ion and electron fluxes independently. The stated 2.84 deuterium value
also assumes a sound-speed convention without the default GAMMA=5/3 multiplier.

The Maxwellian electron collection relation is derived in the primary research
tutorial [Myra, 2021](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/tutorial-on-radio-frequency-sheath-physics-for-magnetically-confined-fusion-devices/37BDBCF6273303327AD6EC2DD4E80533).
That simple relation is itself not a universal all-incidence closure.

### 8. P2 — Flooring either temperature removes both temperature responses

Locations: `models/model600/mod_boundary_conditions.f90:574-585`;
`models/model600/mod_sheath_current.f90:190-193`.

Ti and Te are floored independently when their values are sampled, but their
Jacobian/derivative response is disabled jointly through one OR flag. If Te is
below its floor while Ti remains above it, the residual still varies with Ti,
yet `c_Ti` is set to zero. The derivative boundary condition also loses its
physical Ti-gradient term. The converse problem occurs when only Ti is floored.

A production-row reproduction with Te fixed at a 1 eV floor and Ti=10 eV gives
a required raw Ti coefficient of 8.353619656e-4, but the active flag returns zero.
Use independent temperature active sets and match both spatial and Newton
derivatives to the actual clipped values.

### 9. P2 — The reported residual in volts is not a voltage error

Location: `models/model600/mod_sheath_current.f90:426-428`.

The diagnostic converts `res`, already divided by the maximum coefficient,
directly into volts. That normalization involves coefficients of different
variables and changes with the state and chosen normalization. Even the raw
scaled residual is a Newton residual, not the exact distance to the logarithmic
root far from the root.

In the 100 V mismatch example above the diagnostic reports approximately
3739.87 V, not 100 V. Changing `central_density` from 1 to 100, while keeping
Te=10 eV and the same 100 V distance from floating, changes it to 37398.74 V.
Thus it cannot support the promised interpretation of distance from the sheath
characteristic or comparisons between normalizations.

Report an independently defined physical normal-current residual and, where
`j/jsat<1`, the actual voltage defect relative to
`Phi_required = Te/e * [Lambda - log(1-j/jsat)]`. Explicitly distinguish states
with no finite voltage root. Also count owned unique rows: the current diagnostic
increments once per element traversal, before checking whether the row is locally
owned, so `rows` is not a count of independent boundary constraints.

## Additional limitations relevant to the requested outcome

- **The coupled transport BCs still use an incidence parameter.**
  `mod_boundary_matrix_open.f90:77,343-357` adds particle/heat losses proportional
  to `min_sheath_angle`, whose default is one degree. These remain finite at
  tangency although the new field-aligned normal current tends to zero. Heat
  transmission factors also remain independent of the new electron current.
  Current, particle and energy fluxes need a common derivation and balance test.
- **Parallel thermoelectric physics is incomplete for a Braginskii SOL target.**
  `mod_elt_matrix_fft.f90:1549-1563` contains resistivity and an electron-pressure
  term controlled by tauIC, but no independent parallel electron thermal-force
  term. The standard collisional law includes such a term; see the project's
  independent comparator [GRILLIX's Braginskii derivation](https://grillix-2a40ff.pages.mpcdf.de/page/equations/2_braginskii/1_derivation.html).
  A sheath change alone cannot validate the complete thermoelectric loop. This is
  a pre-existing model limitation, not a newly introduced algebraic error.
- **The no-Vpar model extension is not guarded.** Enabling the new flag with
  `with_vpar=.false.` now enters the shared geometry block and reads
  `values(...,var_vpar)` at `mod_boundary_conditions.f90:537-538`, where var_vpar
  is zero. Reject unsupported configurations or guard those accesses.
- **Higher-order elements are not fully covered.** The new ownership rule skips
  all old u Dirichlet trace derivatives, but replaces only value and first
  derivative rows. For n_order>=5, additional trace DOFs need an explicit boundary
  treatment; they cannot be assumed constrained by these two equations.
- **Exact tangency needs a physical model, not merely finite arithmetic.**
  The current module uses a discontinuous sign and nonzero jsat at B.n=0.
  Near tangency, tangential electric fields can compete with parallel collection;
  this is explicitly analyzed by [Geraldini, Brunner and Parra, 2024](https://arxiv.org/html/2401.07385v2).
  No scalar cancellation proves all-angle validity. The exponential electron
  expression should also not be presented as a general electron-attracting
  sheath model when the plasma potential falls below the wall potential.

## Checks performed and confidence limits

- Read the complete branch diff, the supplied Artola derivation, the nodal
  assembly and its matrix helpers, current/potential normalization, the natural
  particle/heat BCs, the induction equation, and the current/velocity diagnostics.
- Compiled the **unchanged production** `mod_sheath_current.f90` with gfortran,
  `-O0 -g -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow`. Used the repository's
  physical constants and minimal parameter/MPI stubs. MPI routines were not
  exercised. The production self-test passes for F0=+3 and -3; it checks local
  smooth-state derivatives and finite coefficients, not the assembled BC system.
- Reproduced the varying-B derivative defect, large fixed-current Newton update,
  temperature-floor response, diagnostic misinterpretation, and toroidal mean
  rectification example. Reproduced the existing Bohm division failure separately
  with the exact expression and floating-point traps.
- Reproduction sources, executable and numerical output are in
  `/tmp/sheath-current-review/`; `results.txt` contains the successful numerical
  checks. The expected grazing trap is a separate executable, `grazing 0`.
- **No full MPI build or JOREK simulation was run.** This checkout has no
  `Makefile.inc`, and no mpif90 was found. No restart data was used to infer an
  actual crash or a stability lifetime. No claim of long-time stability follows
  from these local tests.

## What would establish readiness

First resolve the current/flow/energy boundary formulation and row ownership,
then repair the derivative and nonlinear assembly. Validate a small coupled
psi/zj/u/w/rho/T/Vpar problem, including the actual assembled residual/Jacobian
and compatible initialization. Test positive states with both field signs,
unequal target temperatures, spatially varying |B|, independently floored Ti/Te,
approach to tangency from both signs, mixed boundary corners, and nonlinear
toroidal perturbations. Compare independent normal-current, particle and power
balances and the resulting potential gradients and flows.

Production evidence should then include serial/MPI agreement, timestep and mesh
convergence, and sustained runs from documented inputs without empirical
stabilizer adjustment. Numerical stability must be assessed against the chosen
physical model; eliminating numerical blow-ups cannot mean suppressing physical
instabilities or assuming arbitrary initial states are admissible.
