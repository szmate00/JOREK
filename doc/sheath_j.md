# Sheath current boundary condition (model600, `bcs(i)%sheath_j`)

Branch `nodal-sheath-j-min`: develop + the floating-potential BC with develop's nodal Mach-1 row
(`doc/floating_u.md`) + this. Nothing else.

## The row

On wall edges whose both endpoints are sheath types, at the Gauss points where |b.n| >= sin(`min_sheath_angle`),
the **zj row** carries the sheath current-voltage characteristic (`mod_boundary_matrix_open.f90`):

    zj = j_sat * f(x),   x = Lambda - e (Phi - V_wall) / Te,
    f  = 1 - e^x                              0 <= x <= Lambda
    f  = 1 - e^x - s_i * x                    x < 0       (Phi above floating, ion side)
    f  = 1 - e^Lambda - s_e * (x - Lambda)    x > Lambda  (Phi below the wall, beyond electron saturation)

with a Zbig*dl weight and exact columns on zj, u, rho, Ti, Te (corr_neg-corrected Te and rho, as in every natural
row). j_sat = c_sat rho (+-v_fl/|b.n|)/|B|, v_fl = factor*cs*|b.n| the parallel Bohm flow the nodal Mach-1 row
imposes (factor = the `vpar_smoothing` weight). c_sat carries the sign of F0 like zj, so the current into the wall,
-zj (B_pol.n)/F0, is independent of the field sign. Normalisation: `sheath_j_norm` in `mod_floating_u.f90`.

`sheath_j_ion_slope` (s_i) and `sheath_j_e_slope` (s_e), default 0 (hard saturation on both sides). With s_i = 1 and
s_e = e^Lambda the characteristic is continued along its tangent beyond floating and beyond electron saturation:
a node asked for more than j_sat, or more than the thermal electron current, then sits at a potential a few Te off
floating (Phi ~ Te (Lambda - 1 + j/j_sat) on the ion side) instead of having no root. That is the form that ran
through detachment on `sheath-j-clean` (2026-09-28).

## Which rows the wall nodes keep (`mod_boundary_conditions.f90`)

- **zj**: the value DOF is released (its Dirichlet row dropped) if any incident wall edge carries the row, a
  tangential-derivative DOF only if an edge in its own direction does.
- **u**: only the **value** DOF is released; the potential's level there follows from the vorticity equation. The
  tangential-derivative DOFs keep the floating row u_s = C_T Te_s. Released slope DOFs have no anchor (the
  characteristic fixes the level of u at the Gauss points, not its slope), and on the finely resolved ray-cast target
  segments they formed a sub-element sawtooth of dPhi/ds whose alternating ExB emptied the wall cell next to a
  floating segment; that was the killer of every earlier sheath-current run. The wall potential's slope within an
  element is therefore the floating slope, its level follows the current.
- Below the angle a sheath-type node keeps the floating row on u and Dirichlet on zj (frozen at its initial value).
- psi and w stay Dirichlet everywhere.

## Namelist

```fortran
bcs(1)%sheath_j   = .t.     ! target plates
bcs(4)%sheath_j   = .t.
bcs(5)%sheath_j   = .t.
bcs(3)%floating_u = .t.     ! the rest of the wall floating (current frozen at its initial value)
bcs(9)%floating_u = .t.
mach1_omit_drift  = .t.     ! nodal Mach-1 row without its ExB term: no u column in the Vpar rows
sheath_Lambda     = 3.d0
sheath_j_ion_slope = 1.d0
sheath_j_e_slope   = 20.085537d0   ! exp(sheath_Lambda)
```

Setup checks (initialise_parameters): a type is not both sheath_j and floating_u; `dirichlet%u` and `dirichlet%zj`
stay .true. on sheath types (the release is per DOF, done here); `mach1` on sheath types.

Not on this branch: the weak Bohm row, total-flow wall fluxes, zero current on floating segments (`zj_zero`),
wall-flux cancellations, Te smoothing, separate angle gate, wall diagnostics.

## Tests

`tests/sheath_j/run.sh` (serial, gfortran): builds the production assembler against stub modules; no zj row with
sheath_j off or below the angle, never a u row; every zj/u/rho/Ti/Te column of the zj rows vs central FD at floating,
on the ion side, beyond electron saturation (hard cap: no u column), with vpar_smoothing, and with the continued
characteristic on both sides (worst 1.4e-7); the continued residual is continuous at floating and C1 at the cap; the
saturation current flows into the wall for both signs of F0; with sheath_j off the assembled edge is bitwise
identical to develop's. The DOF release (`mod_boundary_conditions.f90`) needs the full build.
