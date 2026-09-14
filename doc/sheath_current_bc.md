# Sheath current boundary condition (model600)

Imposes the sheath current-voltage characteristic on the electrostatic potential at a material
wall, replacing the Dirichlet `u` rows on the boundary types you flag.

Characteristic and ion saturation current: Artola, *Sheath boundary conditions for the electric
potential in JOREK* (2026-07-30), `potential_BC_JOREK_ionsat_v2.pdf`, eqs. (5) and (6). The
linearisation is deliberately **not** his eq. (17) — see `models/model600/mod_sheath_current.f90`,
which carries the reasoning and is not duplicated here.

> ## Not runnable yet
>
> One blocker remains, and it is not a namelist issue: **boundary row ownership is undefined**.
> `bcs(:)%dirichlet%zj = .true.` by default, so the boundary current is pinned and cannot respond
> to the potential this condition sets. Freeing `zj` is not a one-line change — the psi/zj/u/w
> boundary system has to be defined together, including the current-definition surface terms.
> Until that is settled, enabling the flag constrains `u` against a frozen `j`. See
> `doc/REVIEW_sheath_current_bc_2026-09-14.md` finding 1.

## Switching it on (once ownership is settled)

```fortran
bcs(1)%sheath_j = .t.
bcs(3)%sheath_j = .t.
bcs(4)%sheath_j = .t.
bcs(5)%sheath_j = .t.
bcs(9)%sheath_j = .t.
sheath_Lambda   = 3.d0
```

Retype each index — the `bcs(:)` array form does not work in a namelist. You do not need to touch
`dirichlet%u`; the sheath row claims those rows on the flagged types. `dirichlet%w` is unaffected.

`sheath_Lambda` is the floating sheath drop in units of `Te/e`,
`-0.5*log(2*pi*(me/mi)*(1+Ti/Te))`. It is the only number the condition takes and it is a physical
constant. There is no gate, threshold, clip, slope limiter, exponent clamp, density or temperature
floor, or relaxation gain anywhere in this boundary condition.

**Requirements and stated restrictions**

- `with_vpar` is required, and is checked at the assembly site. `j_sat = c_sat*rho*Vpar` uses the
  actual parallel velocity rather than the Bohm value, which is what removes all geometry from the
  residual.
- The characteristic presumes **ion outflow** at the wall, `Vpar*(B.n) > 0`, which is what
  `bcs(i)%mach1` provides. Where the outflow vanishes so does `j_sat`, and the condition correctly
  stops constraining the potential there.
- At the type-3/type-2 corner the existing exception is preserved: the sheath row is not assembled
  and the Dirichlet `u` rows are left in place, so `u` is never left without an equation.
- `n_order >= 5` is **not** covered: the flag removes all `u` Dirichlet trace-derivative rows but
  only the value and first-derivative rows are replaced. Bicubic (`n_order = 3`) is exact.
- Only the n=0 background is evaluated, as for Mach1. This is a background linearisation and does
  not enforce the nonlinear characteristic on a non-axisymmetric wall state.

## Reading the log

Two lines per matrix construction, whenever at least one sheath row was assembled on an owned node:

```
 [sheath_j] rows=   1248 nsat=     3  ePhi/Te-Lam= -2.104E+00  4.551E+00  j/jsat= -8.12E-01  9.98E-01
 [sheath_j] max Vdef [V]=  3.418E+00  at R,Z=   1.6012  -1.1104  bnd type  4
```

| field | meaning | what is wrong if it moves |
|---|---|---|
| `rows` | sheath rows on locally owned nodes | a drop means boundary nodes stopped qualifying |
| `nsat` | rows on the saturated branch, where the condition becomes `j = j_sat` and stops setting the potential | a handful is expected; a whole target means the sheath is no longer setting `Phi` there |
| `ePhi/Te-Lam` | normalised potential measured from floating | the two targets on opposite signs is the thermoelectric current, and is physical |
| `j/jsat` | current relative to ion saturation | at or past 1 the sheath is being asked for at least the full ion current |
| `max Vdef` | largest **voltage defect**: `u` minus the potential the characteristic requires, from the unnormalised residual, so it is a genuine voltage — with that node's position and boundary type | **this is the number that grows before a crash**, and the boundary type says which node family loses the condition first |

Saturated nodes are excluded from `max Vdef` because there is no voltage root there to be distant
from. If every row saturates, the second line says so instead.

## What has been checked, and what has not

Verified locally with `gfortran -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow`, for both signs
of `F0`, `central_mass` 2 and 2.5, `sheath_Lambda` 3 and 2.84, and `central_density` 1 and 100:

- The module compiles with no diagnostics, and the assembly block compiles against stubs carrying
  the real `mod_assembly` signatures.
- The floating zero-current state is a root of eq. (6) coded independently of the production row.
- **The residual is the voltage defect**: displacing `u` by 137 V returns exactly 137 V. So a node
  100 V off the characteristic asks for a 100 V step, not `exp(100 e/Te)`.
- Every column of both branches matches a central finite difference of the raw residual to 1e-5.
- The two branches join: a hair either side of `X = 0` the normalised rows and residuals agree, and
  the `u` and `Te` columns fade out of the unsaturated row exactly as the saturated row drops them.

**Not done:** MPI build, cluster run, regression comparison, or any statement about stability. No
run of any kind has exercised this code, and the ownership blocker above must be resolved first.

## Review findings and their status

From `doc/REVIEW_sheath_current_bc_2026-09-14.md`:

| # | Finding | Status |
|---|---|---|
| 1 | Activation leaves the boundary current fixed | **Open — blocker.** Documented above. |
| 2 | Mach row's `u_b/psi_b` grazing singularity | Partly addressed: `j_sat` no longer asserts the Bohm value, so this module makes no claim about `Vpar`. The division itself is pre-existing in the Mach row and unfixed on this base. |
| 3 | Derivative row omitted `j_sat*f*B_b/B` | **Fixed structurally**: the residual contains no geometry, so the derivative row is exactly `d/dl` of the value row. Its Jacobian omits value-DOF columns, which is now stated rather than claimed away. |
| 4 | Row rescaling is not nonlinear stability | **Fixed**: the residual is the voltage defect and the update is bounded. Measured: at `x = 10` the old row asked for `1-exp(10)`; the new row asks for exactly `-100 V`. |
| 5 | Pressure term in the model's current reconstruction | Open. `JpolR = (-zj*BR - R*P_Z)/F0` does carry a diamagnetic piece Artola's eq. (2) does not. Needs settling, not assumed absent. |
| 6 | Only the axisymmetric law is evaluated | Open, stated as a restriction above. Not live for `n_tor = 1`. |
| 7 | Fixed Lambda gives `sqrt(Ti+Te)` electron scaling | Open. A regression against the older `mod_sheath_bc.f90`, which had `sheath_get_lambda(Ti,Te)` with derivatives. |
| 8 | Flooring either temperature removed both responses | **Fixed**: only `Te` enters the residual, so there is one independent flag and no `Ti` column at all. |
| 9 | Reported residual was not a voltage; `rows` overcounted | **Fixed**: the unnormalised residual is converted, and `sheath_diag_add` takes the caller's ownership test. |
| — | `with_vpar = .false.` unguarded | **Fixed**: `apply_sheath_j` requires `with_vpar`. |
| — | `n_order >= 5` trace DOFs | Open, stated as a restriction above. |
