# Floating-u transport A/B experiments

Implemented on `floating-u`, based on `72217bcb2`. These are opt-in experiments,
not a claim that the 607-step failure has been reproduced or cured locally.
The prescribed affine potential relation and its normalisation are unchanged.
No timestep changes, nonlinear corrector, density clipping, grid-coordinate changes,
wall recombination, or kinetic-neutral reflection changes are included.

## Switches (all default false)

| Input | Change |
|---|---|
| `floating_u_mach_flux` | Quadrature normal-flow penalty replaces the nodal Mach value/derivative prescription on floating material edges. |
| `floating_u_wall_flux` | Experimental absorbing charged-wall particle/pressure surface corrections based on total normal velocity. |
| `floating_u_rho_stabilise` | Additional conservative, isotropic density-sensitive numerical diffusion in the volume. |
| `floating_u_transport_diag` | Actual assembly-state transport samples and immediate post-update density samples. No equation changes. |

The first supported configuration is model600 with **one toroidal harmonic,
bicubic elements, separate Ti/Te, density and parallel velocity**, an unrefined
conforming mesh, and no **fluid** neutral/impurity equations. Kinetic neutrals are
unaffected and supported in the sense that their coupling is not changed.
Unsupported combinations fail at setup/topology construction instead of silently
using a partial implementation. The checked-in default model settings may differ
from your configured production build: confirm `with_TiTe`, `with_rho`, and
`with_vpar` in that build.

## 1. Normal-flow experiment

Let `Bn=B.n`, `a=Bn/|B|`, `vEn=vE.n`, and `Vpar` be JOREK's coefficient multiplying B.
The added boundary residual is

```
F = a * (Bn*Vpar + vEn - abs(a)*f(a)*cs)
```

integrated against the trace test functions with the existing boundary penalty
scale `1e12`. The existing bulk parallel-momentum equation is retained.
`f` uses the existing smoothing coefficients, evaluated at every boundary Gauss
point rather than only on edges whose endpoint signs differ. Its tiny negative
round-off-offset tail is clamped to zero.

At finite incidence, F=0 is the intended normal sound-flux target including ExB.
There is **no division by Bn**, and no extra `factor/|B|` multiplying an explicit
drift-cancellation velocity. At exact tangency F and its Vpar/u/temperature
derivatives vanish: this boundary term does not force tangential parallel current
or velocity to zero. The bulk Vpar equation remains. It does **not** prove the
finite-incidence target can be achieved pointwise near tangency or resolve the
magnetic presheath. The condition is imposed weakly with a finite penalty.

Magnetic geometry, including smoothing weights, is **Picard-lagged**. The Vpar,
u and sound-speed derivatives are consistent with that lagged residual, including
the derivative of the temperature correction. This is not a full Newton Jacobian
with respect to psi; point 4 (nonlinear timestep corrections) is not implemented.

Exterior sides come from the existing mesh connectivity, not endpoint label OR
tests. Both endpoints must be floating. Test rows are supplied only for endpoint
types with `mach1=true`; deliberate corner exceptions are preserved. Neither the
legacy `mach_one_bnd_integral` nor `no_mach1_bc` may be combined with this flag.

## 2. Absorbing charged-wall experiment

This flag makes an explicit **modelling assumption**: the material wall supplies
no incoming charged-plasma reservoir. It does not modify neutral reflection.
Keep `density_reflection=0` (the charged-fluid coefficient, not the kinetic-neutral
wall action). The existing `min_sheath_angle` loss is retained unchanged.

With `vn=Bn*Vpar+vEn`, `c=cs*min_sheath_angle*pi/180`, define

```
collection speed = max(vn,0) + c
extra particle flux / rho = max(vn,0) + c - vn
extra pressure flux / (rho*T) = gamma_here*max(vn,0) - vn + (gamma_here-1)*c
```

`gamma_here` is the existing JOREK `gamma_sheath_i/e`, **not** the raw Stangeby input.
These are surface **flux differences**: strong volume advection already carries
`rho*vn` and `rho*T*vn`. Replacing the existing natural corrections with these
expressions avoids adding a duplicate ExB particle loss. At zero ExB and outward
parallel flow, the particle/pressure residuals reduce to the existing expressions
for zero charged reflection. With inward total advection, the diffusive boundary
correction cancels charged-particle injection. This is a total-flux boundary model,
not assignment of an arbitrary incoming density.

The separate viscous heat-flux term is retained; the experimental Jacobian includes
its independent normal/mixed Vpar columns. Magnetic geometry is lagged as above.
Finite diffusion/conduction must be able to support the prescribed corrections;
this is not a justified collisionless, zero-diffusion grazing-wall closure.

The artificial type-2 boundary is **not** floated or assigned this wall closure by
a neighbouring floating type-3 corner. Keep your type-2 and type-3 overrides.

## 3. Density-sensitive diffusion

The existing pressure-based shock capturing remains unchanged. The additional
coefficient uses total poloidal velocity and the existing element mapping:

```
rate = abs(v.grad(s)) + abs(v.grad(t))
h = |v| / rate
sensor = h*|grad(rho)| / (abs(rho) + h*|grad(rho)|)
Dadd = 0.5*h*|v|*sensor
```

Zero speed or zero gradient gives zero added diffusion. It is nonnegative and
rotation/coordinate-reversal invariant. `Dadd` is added equally to perpendicular
and parallel density diffusivities, so the increment is isotropic. The same
lagged coefficient enters the residual and Jacobian. It acts throughout the
volume, not only at the boundary: the density hole can start inside the mesh.

This is conservative numerical diffusion, **not a discrete positivity theorem**,
not a limiter, and not a density floor. It can alter resolved density/temperature
profiles. Retain the comparison with this flag off even if it prevents a crash.

## 4. Diagnostics

The old `floating_u_diag` now has the correct ExB inflow/outflow sign. Its old
Peclet column is explicitly labelled `Pe_core`: it still uses `D_perp(1)`.

`floating_u_transport_diag=true` additionally prints:

- Volume minimum rho, maximum advective CFL, and the volume point nearest the probe.
- Actual Dperp/Dpar (including existing SC and low-density replacement), Dadd, total
  poloidal velocity, metric-based CFL and directional Peclet using the actual
  diffusion tensor. Hyperdiffusion is not folded into this second-order Peclet.
- Kinetic density source, strong advective/compressive density contribution, and
  qjac at those same locations. These are selected terms, not a complete integrated
  particle budget.
- Wall minimum/maximum **total** vn and the nearest wall probe, with separate signed
  ExB and parallel components, endpoint types and the Gauss-point potential residual.
- Immediately after each accepted update: interior minimum density and a fixed
  nearest-probe sample with **new rho, reconstructed old rho, and delta at that same
  location**. This is diagnostic only; it neither clips nor rejects the update.

Records carry the four original vertex IDs and side/Gauss indices. Probe selection
depends on fixed geometry, not the evolving solution. Default probe coordinates:

```fortran
floating_u_probe_R = 1.60d0
floating_u_probe_Z = -1.11d0
```

Extrema tolerate duplicated halo contributions without inflated area sums.
MPI-decomposition invariance still requires a production check. Gauss sampling is
not a certificate of positivity everywhere in a high-order element. Diagnostics
cost additional interpolation and output; enable identically in all A/B cases.

## 5. Run matrix

Start every case from the **same original restart**, not a checkpoint produced by
another case. Preserve particle settings, solver settings, w BC, transport inputs,
and your complete original timestep ramp:

```fortran
tstep_n = 1.d-3, 1.d-2, 1.d-1, 3.d-1, 1.d0, 2.d0, 10.d0
nstep_n = 100, 100, 100, 100, 100, 100, 1000000

bcs(1)%floating_u = .true.
bcs(3)%floating_u = .true.
bcs(4)%floating_u = .true.
bcs(5)%floating_u = .true.
bcs(9)%floating_u = .true.
bcs(2)%floating_u = .false.
sheath_Lambda = 3.d0
sheath_V_wall = 0.d0

floating_u_diag = .true.
floating_u_transport_diag = .true.
floating_u_probe_R = 1.60d0
floating_u_probe_Z = -1.11d0

! Baseline A; replace these three values according to the table.
floating_u_mach_flux = .false.
floating_u_wall_flux = .false.
floating_u_rho_stabilise = .false.
```

Do not remove the rest of your existing namelist, especially corner/PFR settings.
No wall recombination is required or enabled. Do not simultaneously set the legacy
`mach_one_bnd_integral=true`. New flags default false even when omitted.

| Case | mach_flux | wall_flux | rho_stabilise | Purpose |
|---|---|---|---|---|
| A | F | F | F | Same-binary baseline, corrected diagnostics only |
| B | T | F | F | Isolate the normal-flow prescription |
| C | F | T | F | Isolate the absorbing total-flux model |
| D | F | F | T | Isolate density-sensitive diffusion |
| E | T | T | F | Consistent flow plus wall collection |
| F | T | T | T | Combined candidate |

First run A and B, then C/D; use E/F to test interactions. This is a small switch
comparison, not a threshold scan. Case C/F uses the explicitly stated absorbing-wall
assumption: do not interpret stabilisation alone as validation of that model.

Keep the same scheduled dt in every case. Compare at equal physical time, including
the 600->601 transition. A first checkpoint around step 1000 tests passage beyond
the old failure, not long-time validation. Record first negative density (including
post-update diagnostics), GMRES history, potential/temperature profiles, density
minimum location, total wall flow and the magnitude of Dadd. A surviving run with
substantially changed profiles is a changed solution, not automatically a fix.

## Validation completed locally

Run `bash tests/floating_transport/run.sh`.

- Compiles and runs the **production** pure transport kernels.
- Compiles the **production boundary assembler**, floating normalisation, exterior
  topology helper and diagnostic module against minimal serial dependencies.
- Finite-difference checks of u, Vpar, rho, Ti and Te trace/normal/mixed columns for
  the new flow/particle/heat contributions, including viscous heat derivatives,
  inward/outward flow and sub-eV temperature corrections.
- Kernel tests for exact finite-incidence drift compensation, the tangency limit,
  flux-difference balances, density sensor rotation/reversal and constant states.
- Identical covered-edge results for labels 1,3,4,5,9; mixed 3/2 edge exclusion.
- Optional `tests/floating_transport/check_syntax.py` checks integration syntax with
  fparser. Existing GNU consecutive-sign notation in the volume source is normalised
  only for that standards-based parser, not changed in the production source.

**Not completed:** full MPI build, cluster restart run, MPI-decomposition comparison,
full coupled nonlinear residual convergence, or a proof of long-time stability or
positivity. The local environment has gfortran but no MPI Fortran compiler/build
configuration. The post-update interpolation path is syntax/interface checked, not
exercised by the boundary fixture.
