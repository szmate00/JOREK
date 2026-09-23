# PR draft: floating-potential boundary condition for model600

Branch `floating-u-clean-pr` into `develop`. One commit, 13 files, +523/-33.

## Summary

Adds a floating-potential boundary condition on the electrostatic potential at material walls,
`Phi - V_wall = Lambda*Te/e`, together with the wall treatment it needs to run out of the box: a weak
(Galerkin) Bohm condition on the parallel velocity that survives a wall potential varying along the
target, one consistent total normal flow in the sheath particle and energy fluxes and in the kinetic
recycling, and an exterior-side filter for the open-boundary integral. Default namelist otherwise;
no stabiliser, no floor, no clip. Everything is off unless `bcs(i)%floating_u` is set; develop runs are
bit-for-bit unchanged.

Reference case (AUG, model600, production timestep ramp 1e-3 .. 10): 2000+ steps from the equilibrium
with the potential floating on all five wall types, drift-compensating Bohm condition, no tuning.

## New features

1. **`bcs(i)%floating_u`** (`mod_boundary_conditions.f90`, `mod_floating_u.f90`). Replaces the Dirichlet
   u rows of boundary type i by the affine relation `u = C_T*Te + C_V*V_wall` on every u trace DOF
   (value and tangential derivative). `C_T = 2*Lambda/a_n`, halved in a single-temperature build;
   `C_V = sqrt(mu0*rho0)/F0`; `a_n = 2*e*F0*sqrt(mu0*rho0)/m_i`, so that `e*Phi/(k_B*Te) = a_n*u/(2*Te)`.
   The Te column and the RHS make the row exact in one linear solve. A setup self-test checks the
   normalisation for the current F0, mass and density and stops the run if it fails.
   Parameters: `sheath_Lambda` (3), `sheath_V_wall` (0 V). `dirichlet%u` must stay `.true.` on those types.

2. **`mach1_weak`: weak Bohm condition** (`mod_boundary_matrix_open.f90`). Replaces the nodal Mach-1
   value and slope rows on `mach1` boundary types by one Galerkin residual per wall Gauss point,
   `res = (B_pol.n)*Vpar - target`, weighted by `d(res)/d(Vpar) = B_pol.n` and projected on the Vpar trace
   basis. The Vpar column carries (B.n)^2, so a grazing point loses authority continuously and the natural
   Vpar condition takes over; no threshold and no division by B.n anywhere. Assembled on edges whose both
   endpoints are `mach1` types. The nodal rows are not assembled under this flag; type 3 gets no Dirichlet
   Vpar row. Columns on Vpar, Ti, Te and the u trace are exact.
   - The target is `max(cs*|b.n| - vE.n, 0)` everywhere (SOLPS non-marginal form): the parallel flow
     supplies the outward normal flow the ExB drift does not and is never asked to reverse. There is no
     marginal variant and no grazing-angle cut: the Vpar row, the wall fluxes (item 3) and the kinetic
     recycling impose the same drift-compatible condition `Vpar*(B_pol.n) + vE.n >= cs*|b.n|` at every point.
   Why weak: the nodal Mach-1 row on develop carries its drift term `factor/Btot*R^2*u_b/psi_b` in the
   value row only, divided by psi_b; it is dormant while u is constant along the wall and cannot take a
   wall potential that varies along it. `floating_u` therefore requires `mach1_weak` (checked at setup).

3. **One total normal flow at the wall** (`mod_boundary_matrix_open.f90`, `mod_particle_wall_interaction.f90`,
   `mod_fields.f90`). Under `mach1_weak` the sheath energy transmission and the density reflection terms
   use the total outgoing normal flow `max(Vpar*(B_pol.n) + vE.n, 0)` instead of the parallel part alone,
   with exact u, Vpar and psi columns; a face the plasma flows away from collects nothing beyond the
   grazing-incidence floor. The kinetic recycling flux uses the same total flow on the outward normal
   (`wall_normal_vector` points inward), so the neutral source equals the fluid ion loss face by face;
   `calc_EBpsiU` returns the fluid ExB velocity `R grad(u) x e_phi` as an optional trailing argument
   (existing callers unchanged). Off the weak route all expressions are unchanged. Overlaps with PR #32,
   which rewrites this call site; to be rebased onto it once merged.

4. **Exterior sides only** (`mod_boundary_edges.f90`, `construct_matrix_mod.f90`). The open-boundary
   integral was applied to every element side whose two endpoints carry a boundary label; such a side can
   be interior. A connectivity table (a side is exterior iff exactly one element owns it, nodes identified
   by their value-DOF index) is built once per matrix construction when `mach1_weak` or a floating type is
   active, and interior sides are skipped. Conforming unrefined meshes only (aborts otherwise).

## Fixes

5. **`calc_NeTevpar` electron temperature** (`mod_fields.f90`). In a two-temperature model600 build the
   kinetic recycling projection read the ion temperature and applied the single-temperature halving;
   it now reads Te itself, T/2 in a single-temperature build. Superseded by PR #32 once merged.

## Files

| file | change |
|---|---|
| `models/model600/mod_floating_u.f90` | new: normalisation (`floating_u_norm`), volts conversion, setup self-test |
| `models/model600/mod_boundary_edges.f90` | new: exterior-side table from connectivity |
| `models/model600/mod_boundary_conditions.f90` | floating row in the Dirichlet loop; nodal Mach-1 block skipped under `mach1_weak`; Dirichlet Vpar skipped on `mach1` types under `mach1_weak` |
| `models/model600/mod_boundary_matrix_open.f90` | weak Bohm row (residual, columns); total-flow measure `fx_*` in the sheath particle/energy rows with u/Vpar/psi columns |
| `matrix/construct_matrix_mod.f90` | build the exterior-side table per construction; skip interior sides in the side loop (model600 only) |
| `models/model600/initialise_parameters.f90` | namelist entries; setup checks (`mach1_weak` needs `with_vpar`, excludes `mach_one_bnd_integral`; `floating_u` needs `mach1_weak` and `dirichlet%u`); self-test call |
| `models/phys_module.f90`, `models/preset_parameters.f90`, `models/mod_log_params.f90`, `communication/broadcast_phys.f90` | the five parameters and `bcs%floating_u`: declaration, defaults, log, MPI broadcast |
| `particles/mod_fields.f90` | optional `v_ExB` output of `calc_EBpsiU`; `calc_NeTevpar` Te fix |
| `particles/mod_particle_wall_interaction.f90` | recycling flux on the total outgoing normal flow under `mach1_weak` |
| `doc/floating_u.md` | one-page description, namelist, conventions |

## Namelist of the reference run

```fortran
bcs(1)%floating_u = .t. ; bcs(3)%floating_u = .t. ; bcs(4)%floating_u = .t.
bcs(5)%floating_u = .t. ; bcs(9)%floating_u = .t.
mach1_weak = .t.
sheath_Lambda = 3.d0
```

## Conventions (for reviewers)

`Phi = +F0*u` in JOREK's right-handed (R,Z,phi) basis; the fluid ExB velocity is `(-R*u_Z, +R*u_R)`.
Outward-positive normal flows; on a wall Gauss point `vE.n = -orient*R*u_s/dl` with
`orient = sign(y_s*n_R - x_s*n_Z)`; `Vpar*(B_pol.n)` is a velocity (v = Vpar*B);
`vpar0*ps0_s*normal_sign3 = Vpar*(B_pol.n)*R*dl` exactly. Branches (`max`) are decided on the n=0 state
per Gauss point and differentiated exactly on the branch taken.

## Known limits

- `n_order >= 5`: only the value and first tangential-derivative u DOFs are replaced.
- The |B| dependence of the weak row on the free normal psi derivative is lagged (the column loop covers
  trace DOFs; measured 9e-4 of the Vpar column on a test element).
- STALE: the reference-run claims above predate the recycling sign fix and the sheath-set flux; to be
  rewritten after the current runs (see doc/floating_u.md, Status).
- Every Jacobian column of the new rows was finite-difference checked in a serial harness during
  development (not part of this PR). No MPI regression case is added.
