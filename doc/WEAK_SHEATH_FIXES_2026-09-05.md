# Weak sheath fixes: implementation and validation status

Branch: `sheath-jsat-vpar-38ab278`; starting point `2cedbf2f5`.

**This is a tested source-level repair, not a demonstrated all-boundary production
solution.** No full JOREK executable or AUG restart was run on the development
machine. The remaining tangent-wall formulation and evolution issues below must
not be concealed by saying that the unit tests establish simulation stability.
The grid coordinates and frames have not been changed.

## Implemented

1. **Normalization and bounded forward law.**
   `Phi_SI = F0*u/sqrt(mu0*rho_ref)` and the field-aligned current conversion is
   `J_n = -zj*(B.n)/(mu0*F0)`. All routes use the shared normalization routine.
   The current and voltage diagnostics have matching signs.

   Define `chi=e*(Phi_plasma-Phi_wall)/(kB*Te)`,
   `z_ref=c_sat*rho*g*cs/B`. The forward target is now

   ```
   z_sheath = z_ion - z_ref * exp(Lambda - max(chi,0))
   ```

   With prescribed Bohm collection, `z_ion=z_ref`. With
   `sheath_jsat_from_vpar`, only the outgoing ion term is changed. Exactly zero
   Vpar is a real stagnation state, not an absent-data sentinel. The default ion
   floor is zero. The electron reference does not get rescaled by solved Vpar.
   With fixed Lambda this remains a calibrated floating-drop model; with local
   Lambda the existing thermal correction is used. This is a **field-aligned
   collection approximation**, not a first-principles magnetic-presheath model.

   Electron collection saturates at chi <= 0. No soft-clipped floating root,
   unbounded ion tail, or artificial electron tail is present. The derivative at
   chi=0 is the repelling-side generalized derivative; away from that kink the
   supplied derivatives differentiate the actual law.

2. **Weak Mach condition on the same quadrature as the current.**
   On selected weak-sheath edges, Mach1 is projected as
   `integral N*(Vpar - g*cs/B) dl/R = 0`. It is accumulated across incident
   edges, including corners, instead of overwriting a nodal corner value for each
   edge visited. It uses the same incidence evaluator and temperature map as the
   current law, with thermal and magnetic Jacobian columns.

   This explicitly adopts a field-aligned entrance condition. It does **not**
   prescribe a total normal outflow and invert it using division by B.n. The
   singular `u_s/psi_s` drift correction is absent from this weak route.
   Normal ExB flow still exists in the bulk equations and must be measured;
   this change does not cancel it, remove it, or add a compensating particle sink.
   Non-weak legacy nodal Mach conditions are not rewritten by this change.

   At a selected corner next to an artificial edge, the artificial-edge Vpar
   derivative Dirichlet can remain, but may not overwrite the shared Mach value.
   Per-node `mach1=.false.` remains respected: no Mach trace row is written on
   that node.

3. **Current-grid exterior connectivity and explicit coverage.**
   The weak route rebuilds an exterior-side table from the conforming element
   connectivity before every main and directly built preconditioner matrix.
   Interior sides with two nonzero boundary labels are excluded. Canonical
   vertex value DOFs identify shared corner vertices; no coordinate tolerance,
   node displacement, or determinant cutoff is used for topology.

   A weak-current edge requires **both** endpoint types to enable
   `sheath_zj_weak`. Thus 1--3 can be selected while 2--3 remains artificial.
   Enable types 3/9 explicitly when their material incident sides are wanted;
   the old type-1-only input no longer silently covers its mixed corner edges.
   Existing side/direction metadata is still used; broader nonconforming and
   higher-order trace transforms are not implemented.

4. **Exact, independently keyed current and Mach row replacement.**
   The full owned scalar row is cleared across every adjacent node block,
   variable and harmonic before replacement. The row is divided by its largest
   absolute coefficient. The result is independent of the old zbig penalty
   magnitude; unlisted volume coefficients cannot survive.

   This is **row equilibration, not full column/basis equilibration**. It does
   not prove good conditioning of the coupled global operator.
   Current and Vpar equations at the same geometric DOF have distinct slots.
   Row storage is allocated from the boundary graph; the per-row column bound
   remains checked. Exhausted storage, nonfinite coefficients or unsupported
   zero rows cause an explicit error, not a silently truncated equation or a
   frozen-zj/free-u fallback.

   Every matrix construction resets its accumulator, including direct PC
   construction. Row owners retain contributions from all their incident local
   elements; no duplicate all-reduce is introduced.

5. **Raw weak moments and diagnostic ownership.**
   Current residuals are integrated without per-Gauss clipping or evolving
   demand weights. Legacy dynamic gates and residual-clipping parameters are
   rejected for the weak route.

   Only one owner per exterior edge contributes current/area/geometry
   diagnostics, and only owned algebraic rows enter row statistics. PC
   construction does not add to physical diagnostics. Electron-saturation area
   uses chi<=0, not the obsolete fitted exponent clamp.
   Trace row totals now include **current plus Mach** equations.

## Input migration — not a validated operating recipe

Old saturation tails and the alleged cross-field floor are rejected, rather than
silently ignored or reinterpreted. Start from the following parameter choices:

```fortran
sheath_sat_slope      = 0.d0
sheath_sat_slope_e    = 0.d0
sheath_v_perp         = 0.d0
sheath_dfdx_min       = 0.d0
sheath_zj_ratio_max   = 0.d0
sheath_weak_wmin      = 0.d0
sheath_weak_rmax      = 0.d0
sheath_weak_detmin    = 0.d0
sheath_weak_ufade     = .false.
sheath_weak_relax     = 1.d0
sheath_weak_beta      = 0.d0
sheath_psi_jacobian   = .true.
sheath_jsat_from_vpar = .true.
sheath_jsat_vpar_min  = 0.d0
T_min_sheath         = -1.d0
bc_natural_open      = .true.
mach_one_bnd_integral = .false.
```

Lambda and wall voltage remain the chosen physical model parameters. For the
specific Phi=3Te calibration, use `sheath_Lambda=3`,
`sheath_Lambda_local=.false.`, `sheath_V_wall=0`; do not confuse this calibration
with the local thermal-Lambda model.

For each selected material type (1, 3, 4, 5, 9 on the supplied grid), the current
constraint requires `sheath_zj_weak=.true.`, `sheath_zj=.false.`,
`dirichlet%zj=.false.`, `dirichlet%u=.false.`, `natural%u=.false.`.
Keep the existing vorticity-definition condition, typically
`dirichlet%w=.true.`, `natural%w=.false.`. Set `mach1=.true.` where the weak
field-aligned entrance constraint is wanted. Type 2 is not automatically selected.
Do not indiscriminately alter density/temperature BCs or recycling settings.

`sheath_X_min` and `sheath_smooth_dX` no longer change the forward law.
`sheath_min_bn` is a diagnostic incidence weight on the weak route, not a cure
for its tangent limit. Nonzero `sheath_jsat_vpar_min` still represents an
explicit artificial ion floor and is warned about.

A restart generated with the old sign or with super-saturation tails is not a
physical steady reference for the corrected law. `sheath_init_u` is only the
existing nominal floating **boundary guess**, not a consistent interior lift,
nor an exact floating projection with independently solved ion flux. Do not
claim it guarantees a gentle startup. Use a positive compatible checkpoint and
inspect the first small accepted steps before extending the timestep.

## What is NOT fixed yet

* The current trace equation still imposes `zj=z_sheath`. At exact tangency,
  zero normal parallel current does not imply zero tangential parallel current.
  The replacement row therefore still lacks a physically justified tangent-wall
  limit. The new finite Mach limit does not resolve that electromagnetic issue.
* A derived weak normal-current closure with the accompanying induction/flow
  conditions and discrete charge/energy balance remains to be implemented and
  validated in a coupled benchmark. Renaming the old natural-u mismatch is not
  a new formulation; the earlier four-step natural-route failures remain relevant.
* No global nonlinear line search, positive-density/temperature step rejection,
  or conservative positivity-preserving transport discretization was added.
  Strict physical saturation can expose an incompatible demanded current;
  artificial super-saturation is deliberately not used to hide it.
* Shared weak Mach/current temperature maps do not by themselves complete the
  particle and heat closure audit. In particular, a separate T_min_sheath or
  the single-temperature variant can differ from existing transport maps.
  No extra ExB particle/enthalpy boundary sink has been added: strong-form
  volume transport must not be double-counted.
* C1 corner regularity, full block/column scaling, a singular volume mapping,
  and simultaneous wall constraints still require coupled-operator tests.
  A nonzero legacy detmin freezes/excludes material nodes; it is not a solution.
* Weak trace support remains axisymmetric, cubic and unrefined, now explicitly
  checked at setup. General harmonic coupling/refinement is not claimed.

## Checks performed locally

Run `bash tests/sheath/run.sh`. It compiles the production wall-law, trace,
row-clearing, sparse-location and exterior-connectivity modules with dependency
fixtures and gfortran `-fcheck=all -ffpe-trap=invalid,zero,overflow`.

The tests cover:

* independently reconstructed 1 eV -> +3 V for either sign of F0;
* current signs for both incidence directions; exact floating and saturation
  limits; continuity at zero Vpar; electron amplitude independent of ion Vpar;
* the floating-drop shift when ion collection halves;
* seven analytic wall-law derivatives across local-Lambda, repelling,
  electron-saturated, incoming/outgoing and floored states (maximum measured
  relative error approximately 2.94e-10);
* odd incidence smoothing, exact zero at tangency, and both smoothing branches;
* cold-state entrance value/tangential thermal Jacobians;
* all-column/all-harmonic row clearing, unowned-row isolation, separate equations
  at one DOF, contribution merging, reset and zbig-independence;
* exterior/shared/reversed edges, corner aliases and topology rebuilding.

The modified integration files also passed Fortran-2008 syntax parsing after
preprocessing for model600 (fparser 0.2.5). This is **not** a linked production
build, a full assembled boundary-Jacobian test, or an MPI/OpenMP simulation test.
This machine has gfortran but no MPI Fortran compiler or configured Makefile.inc.

## Required next acceptance checks

1. Build the full model600 executable with production dependencies and checking
   enabled; finite-difference the assembled current/Mach boundary rows.
2. Verify physical diagnostic totals and assembled owned rows on 1/2/4 MPI ranks
   and multiple OpenMP thread counts, including a main -> PC -> main sequence.
3. Use identical positive restarts for the material-coverage ladder. Begin with
   short small-timestep runs, recording true weak residuals, current and power
   balance, minimum density/temperatures and normal ExB/parallel particle fluxes.
4. Do not certify type-5 tangency or long-run stability until the normal-current
   formulation and nonlinear/positivity acceptance gaps above have been closed.
