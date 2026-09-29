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
row). j_sat = c_sat rho (+-cs)/|B|, the parallel Bohm flux the nodal Mach-1 row imposes on a target (its
`vpar_smoothing` weight acts only across a tangency, where b.n changes sign between an edge's nodes; until
2026-09-29 j_sat carried that weight everywhere, i.e. 0.28-0.92 of the imposed flux along the inner target,
the only difference inside the sheath row from sheath-j-clean). c_sat carries the sign of F0 like zj, so the current into the wall,
-zj (B_pol.n)/F0, is independent of the field sign. Normalisation: `sheath_j_norm` in `mod_floating_u.f90`.

`sheath_j_ion_slope` (s_i) and `sheath_j_e_slope` (s_e), default 0 (hard saturation on both sides, the form that
runs, see below). With s_i = 1 and s_e = e^Lambda the characteristic is continued along its tangent beyond floating
and beyond electron saturation: a node asked for more than j_sat, or more than the thermal electron current, then
sits at a potential a few Te off floating (Phi ~ Te (Lambda - 1 + j/j_sat) on the ion side) instead of having no
root. Not needed for stability once only the value DOF of u is released.

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

## Namelist (the working setup)

```fortran
bcs(1)%sheath_j   = .t.     ! target plates
bcs(4)%sheath_j   = .t.
bcs(5)%sheath_j   = .t.
bcs(3)%floating_u = .t.     ! the rest of the wall floating (current frozen at its initial value)
bcs(9)%floating_u = .t.
mach1_omit_drift  = .t.     ! nodal Mach-1 row without its ExB term: no u column in the Vpar rows
```

Ran from t = 0 through the production timestep ramp and past step 1000 (2026-09-28), AUG-like divertor with
D puff and kinetic recycling, both targets detaching, with sheath_Lambda = 3 and hard saturation (both slopes 0).
The same minimal setup also ran with the weak Bohm row on `sheath-j-clean` (no zj_zero, no wall-flux
cancellation, no slopes, value-only release). Every earlier setup, on either Mach row and with any of those
additions, died between steps 14 and 680 at the inner-target corner; the one change that removed that failure is
releasing only the value DOF of u.

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

## Thermoelectric terms (physics options, independent of the BC)

- `thermoelectric_ohm` (default .f.): the thermal force in Ohm's law, E_par gains -c*grad_par(Te)/e with
  c = `thermoelectric_coef` (0.71, Braginskii Z = 1). Written as the tauIC electron-pressure term with Pe/rho ->
  c*Te, so it shares that term's normalisation; needs tauIC /= 0.
- `thermoelectric_heat` (default .f.): its Onsager partner, the heat flux carried by the current,
  q_e,par = -c*Te*j_par/e along b, in the electron energy equation as -(gamma-1) div(q b) (weak form
  +(gamma-1) int q (b.grad v) dV, the same (gamma-1) as the Ohmic heating). 1/e -> 2*tauIC*F0 as in the thermal
  force, j_par = -zj*F0/(R**2 B). No surface term: this flux carries no heat through the wall (the sheath heat flux
  is the wall closure, as for parallel conduction). Two-temperature model and tauIC /= 0 required. Exact columns
  on Te, zj and psi, poloidal and toroidal parts of b.grad v. Not compiled by the serial harness (volume assembler).

Both use the same coefficient, as the Onsager relation requires; they can be switched separately.

## Wall fluxes on the total outgoing flow (options; sheath-j-clean's form)

- `sheath_heat_total_flow` (default .f.): the Ti/Te sheath heat sinks use the total outgoing normal flow
  max(Vpar B_pol.n + vE.n, 0) R dl instead of develop's parallel measure Vpar psi_s sign, on edges whose both
  endpoints carry the Mach-1 row: an ExB-inflow face loses no sheath energy, an ExB-outflow face loses it at the
  total rate. Written as the difference to develop's terms (bitwise develop with the flag off), exact columns on
  psi, rho, Ti/Te, Vpar and u. The density row is unchanged (its wall sink is the strong-form volume advection plus
  the n cs sin(min_sheath_angle) floor, identical on both branches).
- `recycling_total_flow` (default .f.): the kinetic recycling flux on the total outgoing normal flow,
  n max(-(Vpar B.n_in + v_ExB.n_in), 0) + n cs sin(min_sheath_angle) with the inward wall normal of
  `wall_normal_vector` (particles/mod_particle_wall_interaction.f90), v_ExB the fluid's (-R u_Z, R u_R) from
  `calc_EBpsiU`; and the correct Te in `calc_NeTevpar` (develop's model600 path reads variable 6 = Ti in a
  two-temperature build and halves it, so the recycling cs was low by ~sqrt(2)). Both behind the flag so the
  sheath-j-clean particle side can be switched as one; with it off the particle code is develop's.
