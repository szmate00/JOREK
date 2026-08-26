# SOLPS vs JOREK model600: the temperature equations, term by term

Comparison of the SOLPS-ITER ion and electron internal-energy equations (Eqs. 3-4 of the
source, with the rearranged `q^∇T` forms) against what JOREK model600 actually assembles.

JOREK side read from `models/model600/mod_elt_matrix_fft.f90` (`rhs_ij(var_Ti)` and
`rhs_ij(var_Te)`) with term names from `models/model600/mod_model_settings.f90`
(`Ti_term_names`, `Te_term_names`). Branch `bc-tests`.

**Conventions.** `Phi = -F0*u`, `v_pol = R grad(u) x e_phi`, so `v_R = -R u_Z`, `v_Z = +R u_R`
and hence

    div(v_ExB) = (1/R) d_R(R v_R) + d_Z(v_Z) = -2 d_Z(u)

which is what lets us identify the E×B compression term below. `BigR**2*(a_s*b_t - a_t*b_s)`
is the Poisson bracket `R^2 [a,b]`, i.e. E×B advection.


## 1. Term-by-term correspondence

### Ion energy equation

| SOLPS term | Meaning | JOREK model600 | Status |
|---|---|---|---|
| `(3/2) n_i dT_i/dt` | time evolution | `Ti_Eq__zeta_timevol` (9) | present |
| `(3/2) Gamma_i . grad T_i` | advection | `Vperp.grad_Pi` (2) + `Vpar.grad_Pi` (4) | present, **E×B only** |
| `div q_ix^gradT` | parallel conduction | `parallel_conduct` (5) | present, B-projected anisotropic |
| `div q_iy^gradT` | perpendicular conduction | `perp_conduction` (6) | present |
| `-sum_a p_a div V_a,par` | parallel compression | `gamma_Pi_div_V` (3) | present |
| `-sum_a p_a div V^ExB` | E×B compression | `gamma_Pi_div_V` (3) | present (see §2) |
| `-sum_a p_a div V_a^gradB` | **magnetic-drift compression** | — | **ABSENT** |
| `-sum_a V_a^gradB . grad p_a` | **magnetic-drift advection** | — | **ABSENT** |
| `-(3/2)(T_i - T_0) S_i` | particle-source enthalpy | `neutral_friction` (10), `recomb_thermal_loss` (18) | present |
| `Q_CX` | charge exchange | only via kinetic coupling `aux_energy_source` (15) | **no fluid term** |
| `Q_vis` | viscous heating | `viscopar_heating` (12) + `visco_heating` (14) | present, par + perp |
| `Q_ei` | equipartition | `TiTe_energy_exch` (11) | present |
| `Q_i^fr` | friction heating | `neutral_friction` (10) | present |

### Electron energy equation

| SOLPS term | Meaning | JOREK model600 | Status |
|---|---|---|---|
| `(3/2) n_e dT_e/dt` | time evolution | `Te_Eq__zeta_timevol_ter` (10) | present |
| `(3/2) Gamma_e . grad T_e` | advection | `Vperp.grad_Pe` (2) + `Vpar.grad_Pe` (4) | present, **E×B only** |
| `div q_ex^gradT`, `div q_ey^gradT` | conduction | `parallel_conduct` (5), `perp_conduction` (6) | present |
| `div q_e,par^u` | **current-driven electron heat flux** | — | **ABSENT** |
| `-p_e div V_e,par` | parallel compression | `gamma_Pe_div_V` (3) | present |
| `-p_e div V^ExB` | E×B compression | `gamma_Pe_div_V` (3) | present |
| `-p_e div V_e^gradB`, `-V_e^gradB . grad p_e` | **magnetic drift** | — | **ABSENT** |
| `-(3/2) T_e S_e` | particle-source enthalpy | `ionization_sink` (12) | present |
| `Q_ion` | ionisation cost | `ionization_sink` (12), `imp_ionization` (17) | present |
| `Q_rad` | radiation | `line_radiation` (13), `Brems_radiation` (14), `backg_imp_radiat` (15), `main_imp_radiat` (16) | present, more resolved than SOLPS's single term |
| `-Q_ei` | equipartition | `TiTe_energy_exch` (11) | present |
| `Q_e^fr` | friction/ohmic | `ohmic_heating` (9) = `(gamma-1) eta (zj/R)^2` | present as ohmic only |


## 2. The E×B compression IS there (easy to miss)

`factor(var_Ti,3)` / `factor(var_Te,3)` contains

    + v * (r0 + rimp0*alpha) * T0 * 2*GAMMA * BigR * u0_y

`u0_y = d_Z(u)`, and from the convention above `div(v_ExB) = -2 d_Z(u)`. So this term is exactly
`-gamma * p * div(v_ExB)`, i.e. SOLPS's `-sum_a p_a div V^ExB`. The same `factor(...,3)` also
carries the two parallel-compression pieces. So compression by both the parallel flow and the
E×B flow is fully represented.


## 3. What is genuinely missing

### 3.1 Magnetic (gradB / curvature) drift terms — the largest gap

SOLPS carries four separate magnetic-drift contributions across the two equations:
`p_a div V_a^gradB` and `V_a^gradB . grad p_a` for ions and electrons. **JOREK model600 has
none of them in the temperature equations.** The only perpendicular advection is the Poisson
bracket with `u`, i.e. pure E×B.

Why it matters here: the gradB/curvature drift is the classical driver of in-out divertor
asymmetries, and crucially **its direction reverses with the toroidal field**, which is how
experiments separate drift-driven asymmetry from everything else. A model without it cannot
reproduce the forward/reversed-field asymmetry reversal, and any in-out asymmetry it does
produce must come from E×B and parallel physics alone.

Note this is a statement about the **temperature** equations specifically. The density and
vorticity equations do carry diamagnetic terms (`rho_Eq__diamag_term` and `u_Eq__diamag_term`
are both large - the former is one of the two dominant terms in the density equation), so the
model is not drift-free overall. The energy channel is where the magnetic drift is dropped.

### 3.2 Current-driven electron parallel heat flux `div q_e,par^u`

SOLPS transports electron heat with the parallel current:

    q_e,par^u = 0.71 * (T_e/e) * j_par

JOREK model600 has no such term. Its only current-related electron energy term is ohmic
heating `(gamma-1) * eta * (zj/R)^2`, which is a positive-definite **source**, not a flux -
it cannot move heat from one target to the other.

Why it matters here: the sheath BC diagnostics show a thermoelectric current loop of
**~2.5 kA** circulating between the inner target (ion branch, `ePhi/kTe = 9.96`, +620 A) and
the outer (electron branch, `ePhi/kTe = 0.77`, -2558 A). The heat flux that physically
accompanies that current, `0.71 (T_e/e) I ~ 0.71 * 13 V * 2500 A ~ 23 kW`, is simply not being
carried. In absolute terms that is a small fraction of P_SOL, but it is **antisymmetric between
the targets** - it cools one and heats the other - which is precisely the kind of term that
matters for an in-out asymmetry rather than for the global power balance.

This is the energy-equation partner of the `0.71 grad_par(Te)` thermal force in Ohm's law
(`thermoelectric_ohm`, branch `u-bcs-tests-mate` @ 485c35c64). The two come from the same
Braginskii closure and are normally implemented together; currently neither is on `bc-tests`.

### 3.3 Charge exchange has no fluid term

SOLPS has an explicit `Q_CX` in the ion equation. model600 has no fluid CX sink; CX enters only
through the kinetic neutral coupling (`aux_energy_source`, term 15, fed by `aux_E0_Ti`).

This is arguably *better* physics - kinetic neutrals with `use_kin_cx = .t.` beat a fluid CX
closure - but it means CX energy transfer exists only while the particle model is running, and
its fidelity is tied to the neutral particle statistics rather than to a rate coefficient.

### 3.4 Diamagnetic heat convection

`Vperp.grad_P` is built from the Poisson bracket with `u` alone, so the perpendicular
temperature advection is E×B with no diamagnetic contribution. A `diamag_heat_conv` term
(`+/- GAMMA * 2 tau' * d_Z(n T_s^2)`, the Poletaeva XPR mechanism) was written previously but
is **not on this branch**.


## 4. What JOREK has that SOLPS does not

These are numerical, not physical, and have no SOLPS counterpart:

- `tg_num_terms` (8) - Taylor-Galerkin stabilisation, `O(tstep^2)`
- `ZK_perp_num_term` (7) - fourth-order numerical perpendicular diffusion
- `implicit_heating` (13 / 19) - a smooth heating floor that prevents the temperature going
  negative. Note it is *not* small where the plasma is cold, which is exactly the divertor.
- `power_teleported` (Te 18)
- `aux_*` (15-17 / 20) - kinetic-particle coupling sources for energy, particles and parallel
  momentum

The radiation channel is also **more** resolved in model600 than in the SOLPS equations as
written: line, bremsstrahlung, background-impurity and main-impurity radiation are separate
terms rather than a single `Q_rad`.


## 5. Bearing on the HFSHD work

Ranked by how much they could plausibly matter for an in-out divertor density asymmetry:

1. **Magnetic-drift terms (§3.1)** - the classical asymmetry driver, absent from the energy
   channel. Their absence means the field-reversal test that experiments use to identify
   drift-driven asymmetry cannot be reproduced, and it removes one of the two mechanisms that
   would give an in-out Te difference in the first place.

2. **Current-driven electron heat flux (§3.2)** - antisymmetric between targets by
   construction, and a ~2.5 kA thermoelectric loop is already present in the solution, so the
   term would be active immediately rather than needing to be driven up. It is also the
   cheapest of these to add, and its Ohm's-law partner already exists on another branch.

3. **Diamagnetic heat convection (§3.4)** - written but not merged here.

4. **Fluid CX (§3.3)** - probably not worth adding given the kinetic neutrals.

**Caveat on ordering:** the sheath boundary condition is currently saturated
(`mean|j/jsat| ~ 3.4-5`), which means the wall passes `j_sat` largely regardless of `Phi`.
Until that is resolved, adding energy-equation terms will change the temperatures but the
potential will not respond to them in the intended way. Both §3.1 and §3.2 are worth having,
but neither is a substitute for getting the sheath below saturation first.
