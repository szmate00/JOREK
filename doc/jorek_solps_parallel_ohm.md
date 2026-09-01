# SOLPS vs JOREK model600: the parallel Ohm's law / induction equation, term by term

Companion to `doc/solps_vs_jorek_temperature_equations.md`, which does the same treatment for the
ion and electron energy equations. Same conventions, same normalisation statements.

**Scope.** The JOREK model600 poloidal-flux (induction) equation `rhs_ij(var_psi)` as assembled in
`models/model600/mod_elt_matrix_fft.f90:1557-1577`, branch `thermoelectric-ohm`, with
`WITH_TiTe`, `n_tor = 1`, bicubic Bezier elements. Compared against the SOLPS-ITER parallel
electron momentum balance, Eq. (B.12) of the SOLPS-ITER manual 3.0.9 (p. 355).

**Status of statements.** Everything with a `file:line` was *read*. Everything labelled
**[derived]** is algebra taking the code as ground truth. Section 9 flags things that look
wrong; nothing in sections 1-7 depends on section 9.

**One result up front, because it inverts a working assumption.** This document concludes
`Phi = +F0 * u`, not `-F0 * u`. The evidence is in section 8.1. The earlier `-F0*u` result verified the
*components* of `v_pol` (which are handedness-independent) and then converted to `Phi` assuming a
right-handed `(R, phi, Z)` basis. JOREK's basis is right-handed in the order `(R, Z, phi)`, which
its own curl and current diagnostics confirm. Independently re-derived and confirmed 2026-08-31:
`rot_tmp` is the right-handed curl for coordinates ordered (R,Z,phi) with scale factors (1,1,R),
and applying it to the code's own `A` reproduces the code's own `B` exactly; since `B = curl A`
is not optional, the basis is settled. Corroborated three ways: the diamagnetic term acquires the
Braginskii sign only under `+F0*u`; `Jtor = -zj0/BigR` matches the third curl component in that
basis; and the measured `min u ~ -1320` becomes a potential WELL, consistent with the observed
trapping of positive impurity ions.

---

## 1. Conventions and normalisations

### 1.1 Coordinates and handedness

JOREK stores 3-vectors in the order `(R, Z, phi)` and that basis is **right-handed**:

$$\mathbf{e}_R \times \mathbf{e}_Z = \mathbf{e}_\varphi,\qquad
\mathbf{e}_Z \times \mathbf{e}_\varphi = \mathbf{e}_R,\qquad
\mathbf{e}_\varphi \times \mathbf{e}_R = \mathbf{e}_Z .$$

Equivalently, JOREK's $\mathbf{e}_\varphi$ is **minus** the $\mathbf{e}_\varphi$ of the usual
right-handed cylindrical system $(R,\varphi,Z)$.

**PRIMARY SOURCE -- the JOREK wiki, "Basic Cylindrical Coordinate System"**, which defines
$(u^1,u^2,u^3) = (R,Z,\phi)$ by

$$x = R\cos\phi,\qquad y = -R\sin\phi,\qquad z = Z$$

and states explicitly: *"Thus, $\phi$ goes **clockwise** if looked at from above!"* The minus on
$y$ is the whole story. It gives

$$\mathbf e_R = (\cos\phi,\,-\sin\phi,\,0),\qquad
\mathbf e_\varphi = (-\sin\phi,\,-\cos\phi,\,0),\qquad
\mathbf e_Z = (0,0,1)$$

$$\Rightarrow\quad \mathbf e_R\times\mathbf e_Z
= (\cos\phi,-\sin\phi,0)\times(0,0,1) = (-\sin\phi,-\cos\phi,0) = \mathbf e_\varphi\quad\checkmark$$

Matching a physical point between the two conventions gives $\phi_J = -\phi_{\rm std}$, hence
$\mathbf e_\varphi^{J} = -\mathbf e_\varphi^{\rm std}$ pointwise. (The wiki's table of contents
also has a COCOS section -- consult it before comparing signs against any equilibrium code, SOLPS
included.)

The code confirms the same thing four more times, independently:

| Evidence | file:line |
|---|---|
| `vvector = [v_R, v_Z, v_phi]` -- component order is `(R,Z,phi)` | `particles/mod_fields.f90:585` |
| `rot_tmp` is exactly the right-handed curl in $(R,Z,\varphi)$: `rotA(1)=dA(3,2)-dA(2,3)/R`, `rotA(2)=dA(1,3)/R-dA(3,1)-A(3)/R`, `rotA(3)=dA(2,1)-dA(1,2)` | `particles/mod_fields.f90:590-596` |
| applied to `A = (-F0 Z/(2R), F0 ln R/2, psi/R)` it returns exactly `B = [psi_Z, -psi_R, F0]/R`, which the code also sets directly. In a $(R,\varphi,Z)$-handed reading the same $A$ gives $B_\varphi=-F_0/R$, contradicting the code. | `particles/mod_fields.f90:775, 791` |
| `Jtor = -zj0/BigR` and `currdens = -zj0/R/mu0` -- the physical toroidal current density is **minus** $\Delta^*\psi/(\mu_0 R)$, which is only true in a $(R,Z,\varphi)$-handed basis | `diagnostics/new_diag/mod_expression.f90:1354, 1901` |
| `JpolR = (R*BP_Z - BZ_p)/R` = $(\nabla\times\mathbf B)_R$ in $(R,Z,\varphi)$ | `diagnostics/new_diag/mod_expression.f90:1304` |

Everything below is written in this basis. Poloidal Poisson bracket:

$$[a,b] \equiv \partial_R a\,\partial_Z b - \partial_Z a\, \partial_R b .$$

In the code `(a_s*b_t - a_t*b_s) = xjac * [a,b]`, so terms written with `_s/_t` derivatives and
*no* explicit `xjac` carry the same integration weight as terms written with `_x/_y` and an
explicit `xjac`. `_x` = $\partial_R$, `_y` = $\partial_Z$, `_p` = $\partial_\varphi$.

### 1.2 Field, current, potential

$$\mathbf{B} = \frac{1}{R}\left(\psi_Z,\; -\psi_R,\; F_0\right)
\qquad\text{i.e.}\qquad
\mathbf{B} = \frac{F_0}{R}\mathbf{e}_\varphi + \frac{1}{R}\nabla\psi\times\mathbf{e}_\varphi$$

`particles/mod_fields.f90:791`. $F_0 = R B_\varphi$ is a constant (reduced MHD: no poloidal
current perturbation).

$$\mathbf{A} = \left(-\tfrac{F_0 Z}{2R},\; \tfrac{F_0}{2}\ln R,\; \tfrac{\psi}{R}\right)
\qquad\Rightarrow\qquad A_\varphi = +\psi/R$$

`particles/mod_fields.f90:775`.

$$zj \equiv \Delta^*\psi = \psi_{RR} - \frac{\psi_R}{R} + \psi_{ZZ},
\qquad \mu_0 j_\varphi = -\frac{zj}{R}$$

$zj=\Delta^*\psi$ is read from the weak `var_zj` equation
`rhs_ij(var_zj) = -(v_x*ps0_x + v_y*ps0_y + v*zj0)/BigR * xjac`
(`models/model600/mod_elt_matrix_fft.f90:1657`), whose strong form is
$zj/R = \nabla\!\cdot\!(R^{-1}\nabla\psi) = \Delta^*\psi/R$ **[derived]**. The minus sign in
$\mu_0 j_\varphi$ is `diagnostics/new_diag/mod_expression.f90:1901`.

$$\mathbf{v} = R\,\mathbf{e}_\varphi\times\nabla u + v_{\rm par}\,\mathbf{B}
= \left(-R u_Z,\;+R u_R,\; \tfrac{F_0}{R}v_{\rm par}\right)$$

`particles/mod_fields.f90:582-584`. Note $R\,\mathbf{e}_\varphi\times\nabla u
= -R\,\nabla u\times\mathbf{e}_\varphi$ in this basis, i.e. the JOREK reference paper's
$\mathbf{v}_{\rm pol}=-R\nabla u\times\mathbf{e}_\varphi$ is reproduced *exactly*.
`Vpar` is $v_\parallel/|B|$, confirmed at `models/model600/mod_boundary_conditions.f90:617`,
`diagnostics/new_diag/mod_expression.f90:1680`.

$$\boxed{\ \Phi = +F_0\,u\ }\qquad\text{[derived, see section 8.1]}$$

### 1.3 Parallel-gradient operator (as coded)

$$\mathbf{B}\cdot\nabla f
= \frac{1}{R}\left(\frac{F_0}{R}\partial_\varphi f + f_R\psi_Z - f_Z\psi_R\right)
= \frac{F_0}{R^2}\partial_\varphi f - \frac{1}{R}[\psi,f]$$

This is literally `Bgrad_T = (F0/BigR*T0_p + T0_x*ps0_y - T0_y*ps0_x)/BigR`,
`models/model600/mod_elt_matrix_fft.f90:1533`. And $\nabla_\parallel f = |B|^{-1}\mathbf B\cdot\nabla f$.

`BB2` $=|B|^2=(F_0^2+\psi_R^2+\psi_Z^2)/R^2$, `models/model600/mod_elt_matrix_fft.f90:1537`.

### 1.4 Normalisation (identical to the temperature document)

With $n_0 \equiv$ `central_density` $\times 10^{20}\,$m$^{-3}$,
$m_i\equiv$ `central_mass` $\times$ `ATOMIC_MASS_UNIT`,
$\rho_0 = n_0 m_i$ (`models/mod_plasma_functions.f90:27-29`):

| quantity | JOREK to SI | source |
|---|---|---|
| length | metres, unchanged | -- |
| $\psi$ | T m$^2$, unchanged | -- |
| time | $t_{\rm SI} = t_J \sqrt{\mu_0\rho_0}$, `sqrt_mu0_rho0` | `mod_plasma_functions.f90:28` |
| velocity | $v_{\rm SI} = v_J/\sqrt{\mu_0\rho_0}$ | -- |
| density | $n_{\rm SI} = \rho_J\, n_0$ | `r0_corr=1` $\equiv n_0$ |
| pressure | $p_J = \mu_0 p_{\rm SI}$ | -- |
| temperature | $T_J = 2.01\times10^{-5}\cdot$`central_density`$\cdot T[{\rm eV}]$ | `models/phys_module.f90:1018` (comment) |
| resistivity | $\eta_{\rm SI} = \eta_J\sqrt{\mu_0/\rho_0}$ = $\eta_J\cdot$`sqrt_mu0_over_rho0` | `mod_plasma_functions.f90:46` (`eta_Spitzer`) |
| potential | $\Phi_{\rm SI}[\mathrm V] = F_0 u_J/\sqrt{\mu_0\rho_0}$ | `diagnostics/new_diag/mod_expression.f90:1758` |

**Diamagnetic coefficient.** `tauIC` is a free namelist input. Its physical value is

$$\tau_{IC}^{\rm nom} = \frac{m_i}{2\,e\,F_0\sqrt{\mu_0\rho_0}}\qquad\text{(dimensionless)}$$

`models/mod_plasma_functions.f90:45`. model600 always uses `tauIC*2.`, i.e. the combination
$2\tau_{IC} = m_i/(e F_0\sqrt{\mu_0\rho_0})$ when set to nominal. The key identity used
repeatedly below **[derived]**:

$$\frac{2\tau_{IC}\,p_J}{\rho_J} \;\longleftrightarrow\; \frac{p_e^{\rm SI}}{e\,n_{\rm SI}\,F_0}
\quad\text{(a velocity stream function, same units as } u).$$

The $\tfrac12$ in $\tau_{IC}^{\rm nom}$ is the $T_i=T_e$ split: single-temperature models use
`tauIC` with the *total* $p$, two-temperature model600 uses `tauIC*2.` with $p_e$
(compare `mod_elt_matrix_fft.f90:2448-2455` against `:2435-2437`).

> **Caveat.** The header comment in `namelist/model333/injet_grid:588-591` states
> `tauIC = m_i/(2 e F0) = 1.482419e-09` -- it omits `sqrt_mu0_rho0` and is dimensionally wrong
> (it has units of s/m). Dividing by `sqrt(rho0*mu0) = 5.377472e-07` from the same header gives
> $\tau_{IC}\approx 2.76\times10^{-3}$, which is the correct dimensionless value and is
> consistent with the $10^{-3}$ used in the model600 regression tests
> (`reg_tests/testcases/inxflow_600_vpar_TiTe_all/input:14`). Trust
> `mod_plasma_functions.f90:45`, not the namelist comment.

### 1.5 Time discretisation and the `factor(...)` mechanism

The assembled system is JOREK's linearised Crank-Nicolson/Gears scheme

$$(1+\zeta)\,\delta\psi^{n} - \zeta\,\delta\psi^{n-1}
= \Delta t\left[\theta\,\frac{\partial F}{\partial \psi}\delta\psi^n + F(\psi^n)\right],$$

with $\theta=$ `time_evol_theta`, $\zeta =$ `time_evol_zeta`$\cdot 2\Delta t/(\Delta t+\Delta t_{\rm prev})$
(`mod_elt_matrix_fft.f90:292`). The $(1+\zeta)$ sits in
`amat(var_psi,var_psi) = v*psi/BigR*xjac*(1.d0+zeta)` (`:2402`), and the $-\zeta\delta\psi^{n-1}$
appears on the RHS as term 5.

`factor(var_psi,i)` is the **term-isolation diagnostic**: when `get_terms` is present, all
`factor` entries are zeroed except term `i`, which is set to `1/tstep`
(`mod_elt_matrix_fft.f90:1543-1551`). This is what produces the per-term RHS comparisons used in
this campaign. Term names: `models/model600/mod_model_settings.f90:94-100`.

---

## 2. The JOREK psi equation as implemented

### 2.1 Weak form, verbatim

`models/model600/mod_elt_matrix_fft.f90:1557-1577`. `v` is the test function, `xjac` the
element Jacobian, `_0` suffixes denote the current (explicit) solution.

```fortran
rhs_ij(var_psi) = v * eta_T  * (zj0 - current_source(ms,mt) - Jb)/ BigR           * xjac * tstep * factor(var_psi,1) &
          + v * (ps0_s * u0_t - ps0_t * u0_s)                                            * tstep * factor(var_psi,2) &
          - v * F0 / BigR  * u0_p                                                 * xjac * tstep * factor(var_psi,2) &
          + eta_num_T * (v_x * zj0_x + v_y * zj0_y)                               * xjac * tstep * factor(var_psi,3) &
          - v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * (ps0_s * Pe0_t - ps0_t * Pe0_s) * tstep * factor(var_psi,4) &
          + v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3 * Pe0_p                    * xjac * tstep * factor(var_psi,4) &
          - v * tec * tauIC*2./BB2 * F0**2/BigR**2 * (ps0_s * Te0_t - ps0_t * Te0_s)     * tstep * factor(var_psi,7) &
          + v * tec * tauIC*2./BB2 * F0**3/BigR**3 * Te0_p                        * xjac * tstep * factor(var_psi,7) &
          + zeta * v * delta_g(mp,var_psi,ms,mt) / BigR                           * xjac         * factor(var_psi,5) &
          - v * eta_T  * aux_jre_ind / BigR                                       * xjac * tstep * factor(var_psi,6)
```

with `amat(var_psi,var_psi) = v * psi / BigR * xjac * (1.d0 + zeta) - ...` (`:2402`).

### 2.2 Integration by parts -- term by term

**This must be established before any strong form can be written.** The test function appears as
bare `v` in every term except one:

| term | test function | integrated by parts? |
|---|---|---|
| 1 `eta_J` | `v` | **no** |
| 2 `B.grad_u` (both lines) | `v` | **no** |
| 3 `eta_num_term` | `v_x`, `v_y` | **yes, once** (in the poloidal plane only) |
| 4 `diamag_term` (both lines) | `v` | **no** |
| 5 `zeta_timevol_term` | `v` | **no** |
| 6 `RE_coupling` | `v` | **no** |
| 7 `thermal_force` (both lines) | `v` | **no** |
| mass matrix | `v` | **no** |

So only the hyper-resistivity is in weak form. Its surface term
$\oint \eta_{\rm num}\, v\, \partial_n(zj)\,\mathrm dS$ is dropped; `psi` and `zj` normally carry
Dirichlet conditions (`bcs(bnd_type)%dirichlet%psi`,
`models/model600/mod_boundary_conditions.f90:376`), which makes $v=0$ on the boundary and the
surface term vanish for the Galerkin rows that are actually solved.

Also note the hyper-resistivity term is the **only** term with no $1/R$ weight, and it uses only
`v_x, v_y` -- it is a **poloidal-plane** operator, with no $\partial_\varphi$ contribution.

### 2.3 Strong form **[derived]**

Divide out the common $\int v\,\mathrm{xjac}$, drop `tstep`/`factor`, use
`(a_s*b_t - a_t*b_s) = xjac*[a,b]`, and multiply through by $R$ (the mass matrix carries $1/R$):

$$
\begin{aligned}
\frac{\partial \psi}{\partial t}
=\;& \underbrace{\eta\left(zj - j_{\rm src} - J_b\right)}_{\textbf{1 }\texttt{eta\_J}}
\;\underbrace{-\;\eta\, j_{\rm RE}}_{\textbf{6 }\texttt{RE\_coupling}}\\[4pt]
&\underbrace{+\;R\,[\psi,u] \;-\; F_0\,\partial_\varphi u}_{\textbf{2 }\texttt{B.grad\_u}\;=\;-R^2\,\mathbf B\cdot\nabla u}\\[4pt]
&\underbrace{-\;R\,\nabla_{\rm pol}\!\cdot\!\left(\eta_{\rm num}\nabla_{\rm pol}\, zj\right)}_{\textbf{3 }\texttt{eta\_num\_term}}\\[4pt]
&\underbrace{+\;\frac{2\tau_{IC}\,F_0^{2}}{\rho_{\rm corr}\,|B|^{2}}\;\mathbf B\cdot\nabla p_e}_{\textbf{4 }\texttt{diamag\_term}}
\;\underbrace{+\;\frac{c_{\rm te}\,2\tau_{IC}\,F_0^{2}}{|B|^{2}}\;\mathbf B\cdot\nabla T_e}_{\textbf{7 }\texttt{thermal\_force}}\\[4pt]
&\underbrace{+\;\frac{\zeta}{\Delta t}\,\Delta\psi^{\,n-1}}_{\textbf{5 }\texttt{zeta\_timevol\_term}}
\end{aligned}
$$

The identification of terms 4 and 7 as $\mathbf B\cdot\nabla$ is the algebra
$-\tfrac1R[\psi,f] + \tfrac{F_0}{R^2}\partial_\varphi f = \mathbf B\cdot\nabla f$ applied to the
two lines of each block; the relative factor of $F_0/R$ between the two lines of each block is
exactly what makes this work, and it is a useful check that neither line has been mistyped.

### 2.4 What each symbol is

| symbol | meaning | file:line |
|---|---|---|
| `eta_T` | local resistivity, JOREK units, $T$- and $Z_{\rm eff}$-dependent | `mod_elt_matrix_fft.f90:1154`, `models/mod_plasma_functions.f90:135-233` |
| `zj0` | $\Delta^*\psi$ at the current iterate | `mod_elt_matrix_fft.f90:1657` |
| `current_source(ms,mt)` | `keep_current_prof` source, section 5 | `mod_elt_matrix_fft.f90:467-468`, `models/current.f90` |
| `Jb` | bootstrap current *minus* its initial-profile value; `0` unless `bootstrap` | `mod_elt_matrix_fft.f90:1195-1218` |
| `aux_jre_ind` | kinetic runaway-electron current, `= aux_jre` if `use_rep` **and** `.not. keep_current_prof`, else `0` | `mod_elt_matrix_fft.f90:768-773` |
| `eta_num_T` | hyper-resistivity, optionally $T^{-3}$ or $\psi_N$-profiled | `models/mod_plasma_functions.f90:255-268` |
| `r0_corr` | `corr_neg_dens(r0)`, positivity-corrected total density | `mod_elt_matrix_fft.f90:587` |
| `Pe0` | $(\rho + \rho_{\rm imp}\alpha_e)T_e$, i.e. $n_e T_e$ | `mod_elt_matrix_fft.f90:5687`, `alpha_e` at `:5458` |
| `tec` | `thermoelectric_coef` if `thermoelectric_ohm`, else `0.d0` | `mod_elt_matrix_fft.f90:185, 269-270` |
| `delta_g(mp,var_psi,..)` | previous time step's $\Delta\psi$ | `mod_elt_matrix_fft.f90:1574` |

### 2.5 Resistivity, in full

`models/mod_plasma_functions.f90:135-233`, called at `mod_elt_matrix_fft.f90:1154-1157`:

$$\eta = \eta_0\left(\frac{T_e^{\rm corr}}{T_{e,0}}\right)^{-3/2}\times
\underbrace{\frac{Z_{\rm eff}\frac{1+1.198Z_{\rm eff}+0.222Z_{\rm eff}^2}{1+2.966Z_{\rm eff}+0.753Z_{\rm eff}^2}}
{\frac{1+1.198+0.222}{1+2.966+0.753}}}_{\text{only if }\texttt{with\_impurities}}
\times \underbrace{\frac{\ln\Lambda}{\ln\Lambda_{\rm centre}}}_{\text{only if }\texttt{eta\_coul\_log\_dep}}$$

with clamps (`:172-186`):

- if $T^{\rm corr} >$ `T_max_eta`, freeze at $\eta_0(T_{\max}/T_0)^{-3/2}$ -- default `T_max_eta = 1.d99`, i.e. off (`models/preset_parameters.f90:42`);
- if the **raw** $T <$ `T_min`, freeze at $\eta_0(T_{\min}/T_0)^{-3/2}$ and zero the derivatives.

The floor uses the global `T_min` on this branch. `T_min` is in JOREK units,
`2.01e-5 * central_density * T_min[eV]` (`models/phys_module.f90:1018`). A separate `T_min_eta`
exists only on `keep-current-prof-cutoff`, not here.

`eta_T_dependent = .false.` gives $\eta=\eta_0$ everywhere. A *second, independent* resistivity
`eta_T_ohm` with its own `T_max_eta_ohm` is computed at `:1160-1163` and used only for Ohmic
heating in the $T_e$ equation -- the psi equation and the Ohmic heating term can therefore use
different resistivities.

---

## 3. Parallel Ohm's law in SI

### 3.1 Rearranging into $E_\parallel = \dots$ **[derived]**

With $A_\varphi = \psi/R$ and $b_\varphi = B_\varphi/|B| = F_0/(R|B|)$,

$$E_\parallel = \mathbf b\cdot\mathbf E = -\mathbf b\cdot\nabla\Phi - b_\varphi\,\partial_t A_\varphi
= -\mathbf b\cdot\nabla\Phi - \frac{F_0}{R^{2}|B|}\,\frac{\partial\psi}{\partial t}$$

$$\Longleftrightarrow\qquad
\frac{\partial\psi}{\partial t} = -\frac{R^{2}|B|}{F_0}\left[E_\parallel + \mathbf b\cdot\nabla\Phi\right].$$

Substituting the strong form of section 2.3 and using $\Phi=F_0u$ (which makes the `B.grad_u` block
cancel identically against $-\mathbf b\cdot\nabla\Phi$), $\mu_0 j_\varphi=-zj/R$, and the
reduced-MHD closure $j_\parallel = b_\varphi j_\varphi$:

$$\boxed{\;
E_\parallel \;=\; \eta_\parallel\left(j_\parallel - j_\parallel^{\rm src} - j_\parallel^{\rm boot} - j_\parallel^{\rm RE}\right)
\;-\;\frac{1}{e\,n_e}\nabla_\parallel p_e
\;-\;\frac{0.71}{e}\nabla_\parallel (k T_e)
\;+\;b_\varphi\,\nabla_{\rm pol}\!\cdot\!\left(\eta_{\rm num}\nabla_{\rm pol}\,zj\right)
\;+\;\mathcal T_\zeta \;}$$

Every coefficient in SI:

| JOREK term | SI form | how the coefficient converts |
|---|---|---|
| 1 `eta_J` | $\eta_\parallel j_\parallel$ | $\eta_{\rm SI} = \eta_J\sqrt{\mu_0/\rho_0}$; $j_\parallel = b_\varphi\, j_\varphi = -\dfrac{F_0\,zj}{\mu_0 R^2|B|}$ |
| 2 `B.grad_u` | cancels against $-\mathbf b\cdot\nabla\Phi$ -- it **is** $E_\parallel^{\rm ideal}$ | $\Phi = F_0 u$, $\Phi_{\rm SI}=F_0 u_J/\sqrt{\mu_0\rho_0}$ |
| 4 `diamag_term` | $-\dfrac{1}{e n_e}\nabla_\parallel p_e$ | $\dfrac{2\tau_{IC}F_0^{2}}{\rho_{\rm corr}|B|^{2}} \leftrightarrow \dfrac{R^{2}}{F_0 e n_e}$, using $F_0^2/|B|^2 = R^2 b_\varphi^2 \simeq R^2$ |
| 7 `thermal_force` | $-\dfrac{c_{\rm te}}{e}\nabla_\parallel (kT_e)$, $c_{\rm te}=0.71$ | same as above with $p_e/n \to c_{\rm te}T_e$; the $1/n$ cancels, which is why term 7 correctly has **no** `r0_corr` and no `amat(var_psi,var_rho)` entry (code comment `:1566-1570`) |
| 3 `eta_num_term` | $+b_\varphi\nabla_{\rm pol}\!\cdot(\eta_{\rm num}\nabla_{\rm pol} zj)$ | purely numerical; $\propto \nabla_{\rm pol}^2 j_\varphi$, poloidal only |
| 5 `zeta_timevol` | $\mathcal T_\zeta$, time-discretisation only | Gears/BDF2 residual, not physics |
| 6 `RE_coupling` | $-\eta_\parallel j_\parallel^{\rm RE}$ | folded into term 1's bracket |

The exact geometric factors, before the $b_\varphi^2\simeq1$ simplification, are
$F_0^2/|B|^2 = R^2 b_\varphi^2$ and $R^2|B|/F_0 = R/b_\varphi$. With $|b_n|\sim0.46^\circ$ field
line pitch at the targets (`doc/jorek_vpar_bc_drifts.tex`), $b_\varphi^2 = 1 - b_{\rm pol}^2$
differs from 1 by well under a per cent everywhere in the SOL -- the code's $F_0^2/|B|^2$ form is
the *exact* one and the $R^2$ form is the approximation, not the other way round.

### 3.2 Two structural approximations worth naming

1. **$j_\parallel \simeq b_\varphi j_\varphi$.** Only the toroidal current feels resistivity. The
   poloidal current's contribution $\mathbf j_{\rm pol}\cdot\mathbf B_{\rm pol}$ to $j_\parallel$
   is dropped. This is inherent to reduced MHD, not a model600 choice.
2. **No $\partial_\varphi$ in the hyper-resistivity.** Term 3 is a poloidal-plane operator.

---

## 4. The SOLPS form

SOLPS-ITER manual 3.0.9, Appendix B, p. 355, Eq. (B.12) -- "the parallel momentum balance for
electrons has the standard form":

$$j_\parallel = \sigma_\parallel\left[\frac{b_x}{e}\left(\frac{1}{h_x}\frac{\partial n T_e}{n\,\partial x}
+ 0.71\frac{\partial T_e}{\partial x}\right) - \frac{b_x}{h_x}\frac{\partial \Phi}{\partial x}\right]$$

$x$ is the poloidal coordinate, $h_x$ its metric coefficient, $b_x=B_x/B$ the poloidal pitch, so
$(b_x/h_x)\partial_x = \nabla_\parallel$ and $T$ is in energy units. Rearranged into the same
shape as section 3.1, with $\eta_\parallel = 1/\sigma_\parallel$:

$$\boxed{\;E_\parallel^{\rm SOLPS} \;=\; -\nabla_\parallel\Phi
\;=\;\eta_\parallel j_\parallel \;-\;\frac{1}{e n_e}\nabla_\parallel p_e\;-\;\frac{0.71}{e}\nabla_\parallel T_e\;}$$

**SOLPS's Ohm's law is electrostatic** -- there is no $\partial\mathbf A/\partial t$. $\Phi$ is
obtained not from this equation but from **current continuity**, manual p. 356, Eq. (B.13):

$$\frac{1}{\sqrt g}\frac{\partial}{\partial x}\!\left(\frac{\sqrt g}{h_x}\tilde\jmath_x\right)
+\frac{1}{\sqrt g}\frac{\partial}{\partial y}\!\left(\frac{\sqrt g}{h_y}\tilde\jmath_y\right)=0,
\qquad \tilde\jmath_x = b_z\tilde\jmath_\perp + b_x\tilde\jmath_\parallel,$$

with $\mathbf j = \tilde{\mathbf j}^{\rm (dia)} + \mathbf j^{\rm (in)} + \mathbf j^{\rm (vis)}
+ \mathbf j^{\rm (s)} + \mathbf j_\parallel$ (Eq. B.15), i.e. diamagnetic (B.16-17), inertia +
gyroviscosity (B.22-23), parallel/perpendicular/heat viscosity (B.18-21, B.24-25), and
ion-neutral friction (B.27-33). Boundary closure via `BCPOT` (manual p. 88), of which type 3
and 11 are the sheath conditions.

Two SOLPS knobs have no equation-level analogue in JOREK: `FLAG_SIG`/`PARM_SIG`, an **anomalous
current conductivity**, and `FLAG_ALF`/`PARM_ALF`, an **anomalous thermo-electric current**
(manual p. 91).

The energy-equation partner of the thermal force appears in SOLPS's electron heat flux
$\tilde q_{ex}$ as $-0.71\,b_x j_\parallel T_e/e$ (Eq. B.40, p. 359). JOREK model600 has no such
term -- see `doc/solps_vs_jorek_temperature_equations.md` section 3.2.

**Physical reference point for the HFSHD mechanism**, $\Delta\Phi = \int \eta j_\parallel\,\mathrm dl$
along a field line, is Senichenkov et al., *Contrib. Plasma Phys.* **62**(5-6), e202100177 (2022).

> **Note on the repo file.** `122307_1_online.pdf` in the repository root is **not** the
> Senichenkov paper. It is Loizu, Ricci, Halpern & Jolliet, *Boundary conditions for plasma fluid
> models at the magnetic presheath entrance*, Phys. Plasmas **19**, 122307 (2012). The Senichenkov
> statement above is cited from literature, not verified against a file in this repository.

---

## 5. `keep_current_prof`: the term with no SOLPS counterpart

$j_{\rm src}$ = `current_source(ms,mt)`, computed once per Gauss point at
`models/model600/mod_elt_matrix_fft.f90:467-468` when `keep_current_prof = .true.`
(default `.true.`, `models/preset_parameters.f90:713`). It converts term 1 from $\eta j$ into
$\eta(j-j_0)$ -- a relaxation of the current towards the *initial equilibrium* profile rather than
resistive decay.

### 5.1 The exact implemented form

`models/current.f90:56`:

```fortran
zjz = zFFprime - R*R * (zn * dT_dpsi + dn_dpsi * zT)
```

i.e. **[read]**

$$j_{\rm src}(R,Z) \;=\; \left(FF'\right)_{\rm prof} \;-\; R^{2}\left(n\,\frac{\partial T}{\partial\psi} + \frac{\partial n}{\partial\psi}\,T\right)
\;=\; \left(FF'\right)_{\rm prof} - R^{2}\,\frac{\partial p}{\partial\psi}$$

which is the Grad-Shafranov right-hand side evaluated on the *input analytic profiles* at the
*current* value of $\psi$. Under `WITH_TiTe`, $T = T_i + T_e$ and
$\partial T/\partial\psi = \partial T_i/\partial\psi + \partial T_e/\partial\psi$
(`models/current.f90:42-43`).

### 5.2 The $FF'$ profile, with its two cutoffs and the `FF_1` offset

`models/FFprime.f90`. With $\psi_N = (\psi-\psi_{\rm axis})/(\psi_{\rm bnd}-\psi_{\rm axis})$,
clamped to $[0,2]$ (`:64`):

**(a) base profile + tanh barrier** (analytic branch, `:76-108`)

$$\mathrm{prof}_0(\psi_N) = (FF_0 - FF_1)\left(1 + c_1\psi_N + c_2\psi_N^2 + c_3\psi_N^3\right)
+ \frac{c_6}{2c_8}\,\mathrm{sech}^2\!\left(\frac{\psi_N-c_7}{c_8}\right)\frac{1}{\Delta\psi}$$

$$\mathrm{prof}_1 = \mathrm{prof}_0 \cdot \underbrace{\tfrac12\left(1-\tanh\frac{\psi_N-c_5}{c_4}\right)}_{\text{tanh cutoff towards the SOL}}$$

$c_i = $ `FF_coef(i)`, $c_4=$ `sig_F`, $c_5=$ `psi_barrier`. If `num_ffprime` is set, $\mathrm{prof}_1$
is instead linearly interpolated from a file (`:110-136`).

**(b) $Z$ cutoff beyond the X-point(s)** (`:138-175`), **hard-coded width $\sigma_Z = 0.1$ m**:

$$\mathrm{prof}_1 \;\to\; \mathrm{prof}_1 \cdot
\underbrace{\tfrac12\left(1-\tanh\frac{Z_{X,\rm lower}-Z}{0.1}\right)}_{\text{unless }\texttt{UPPER\_XPOINT}}
\cdot
\underbrace{\tfrac12\left(1-\tanh\frac{Z-Z_{X,\rm upper}}{0.1}\right)}_{\text{unless }\texttt{LOWER\_XPOINT}}$$

**(c) free-boundary rescale** (`:187-196`), only if `freeboundary_equil .and. num_ffprime`:
multiply by `current_FB_fact`.

**(d) the `FF_1` offset, added AFTER everything** -- `models/FFprime.f90:198`:

```fortran
FFprime_profile = FFprime_profile + FF_1
```

$$\left(FF'\right)_{\rm prof} = \mathrm{prof}_1 \cdot (\text{Z cutoffs}) \cdot (\text{FB factor}) \;+\; FF_1$$

**This is the trap.** $FF_1$ is *outside* both tanh cutoffs. Deep in the SOL and in the private
flux region, where the $\psi_N$ tanh and the $Z$ tanh have driven $\mathrm{prof}_1$ to zero, the
current source does **not** go to zero -- it goes to $FF_1 - R^2\partial_\psi p$, i.e. a constant
$FF_1$ plus whatever the (also non-zero-tailed) $n$ and $T$ input profiles give. Since the base
polynomial is scaled by $(FF_0-FF_1)$, $FF_1$ is the intended *edge* value of $FF'$; but that
means a uniform $\eta\,FF_1$ drive survives on open field lines, exactly where $\eta$ is largest.

### 5.3 The confinement mask

`models/current.f90:69-88`, active if `keep_current_prof_confined .and. keep_current_confine_strength > 0`:

$$\text{mask} = \tfrac12\left(1-\tanh\frac{\psi_N - \psi_N^{\rm cut}}{\sigma_{\psi_N}}\right)
\cdot \tfrac12\left(1-\tanh\frac{Z_{X,\rm low}-Z}{\sigma_Z^{\rm kc}}\right)
\cdot \tfrac12\left(1-\tanh\frac{Z-Z_{X,\rm up}}{\sigma_Z^{\rm kc}}\right)$$

$$\text{mask}_{\rm eff} = 1 - s\,(1-\text{mask}),\qquad j_{\rm src}\to j_{\rm src}\cdot\text{mask}_{\rm eff}$$

with $\psi_N^{\rm cut}=$ `keep_current_psin_cutoff`, $\sigma_{\psi_N}=$ `keep_current_psin_sig`,
$\sigma_Z^{\rm kc}=$ `keep_current_z_sig`, $s=$ `keep_current_confine_strength` $\in[0,1]$.
The $Z$ factors are skipped for single-X-point cases per `xcase2`. $s=1$ removes the source
entirely outside the confined region; $s=0$ is bit-identical to the flag being off (and is
short-circuited before the tanh calls). Note this mask *does* catch the $FF_1$ offset, because it
multiplies the assembled `zjz`, whereas the `FFprime.f90` cutoffs do not.

### 5.4 Size and consequences

- Measured in this campaign at **~45% of the psi-equation RHS when uncut**
  (`rhs-term-comparison-cutoff`). It is the second-largest term in the equation.
- It is why $\eta$-scans of the HFSHD mechanism were misleading: the artificial
  $\eta(j-j_0)$ drive on open field lines, not $\eta j_\parallel$, produced the apparent
  $\eta$-dependence (`hfshd-eta-dependence`, bracket completed 2026-07-23: 100% artefact).
- **It also silently disables the runaway-electron coupling**: `aux_jre_ind` is forced to `0`
  whenever `keep_current_prof` is `.true.` (`mod_elt_matrix_fft.f90:768-772`).
- **SOLPS has nothing like it.** SOLPS's $\Phi$ comes from $\nabla\cdot\mathbf j=0$ with no
  reference to an initial current profile. Any comparison of $\eta j_\parallel$ between the two
  codes is invalid while `keep_current_prof` is on and unmasked.

---

## 6. Term-by-term correspondence

| # | JOREK term (`psi_term_names`) | strong form | SOLPS term (B.12) | status | comment |
|---|---|---|---|---|---|
| -- | mass matrix, `amat(var_psi,var_psi)` `:2402` | $\partial_t\psi/R$ | -- | **JOREK only** | SOLPS is electrostatic; there is no $\partial_t\mathbf A$. This is the single biggest structural difference. |
| 1 | `psi_Eq__eta_J` `:1557` | $\eta\,(zj - j_{\rm src} - J_b)$ | $j_\parallel/\sigma_\parallel$ | **present, different form** | Same physics. JOREK: $j_\parallel\simeq b_\varphi j_\varphi$, $\eta$ from Spitzer $T^{-3/2}$ + $Z_{\rm eff}$ + $\ln\Lambda$. SOLPS: $\sigma_\parallel$ Balescu-based. JOREK adds $-j_{\rm src}$ and $-J_b$, which SOLPS does not have. |
| 2 | `psi_Eq__B.grad_u` `:1559-1560` | $-R^2\,\mathbf B\cdot\nabla u$ | $-\sigma_\parallel\nabla_\parallel\Phi$ | **present, different form** | Identical physics: this term *is* $-\mathbf b\cdot\nabla\Phi$ under $\Phi=F_0u$. In JOREK it is an evolution term for $\psi$; in SOLPS it is the unknown $\Phi$ in a $\nabla\cdot\mathbf j = 0$ elliptic solve. |
| 3 | `psi_Eq__eta_num_term` `:1561` | $-R\,\nabla_{\rm pol}\!\cdot(\eta_{\rm num}\nabla_{\rm pol} zj)$ | -- | **JOREK only, numerical** | Hyper-resistivity, $O(k^4)$ damping. Poloidal-plane only. Not physics. |
| 4 | `psi_Eq__diamag_term` `:1563-1564` | $+\dfrac{2\tau_{IC}F_0^2}{\rho_{\rm corr}|B|^2}\mathbf B\cdot\nabla p_e \equiv -\dfrac{1}{en_e}\nabla_\parallel p_e$ | $\sigma_\parallel\dfrac{b_x}{e}\dfrac{\partial(nT_e)}{h_x n\partial x}$ | **present, same physics** | Identical Braginskii term. **Only active if `tauIC` $\neq 0$** -- it is 0 by default (`models/preset_parameters.f90:568`), and the model600 regression tests use `1.d-3` against a nominal $\approx 2.8\times10^{-3}$. Check the campaign namelist. |
| 5 | `psi_Eq__zeta_timevol_term` `:1574` | $+\zeta\Delta\psi^{n-1}/\Delta t$ | -- | **JOREK only, numerical** | Gears/BDF2 second-order time scheme. Active iff `time_evol_zeta` $\neq 0$; with `time_evol_zeta = 0` and $\theta=0.5$ the scheme is Crank-Nicolson and this term is identically zero. |
| 6 | `psi_Eq__RE_coupling` `:1577` | $-\eta\,j_{\rm RE}$ | -- | **JOREK only** | Kinetic runaway-electron current subtracted from the resistive drive. Active only if `use_rep = .true.` **and** `keep_current_prof = .false.` (`:768-772`). Off in this campaign. |
| 7 | `psi_Eq__thermal_force` `:1571-1572` | $+\dfrac{c_{\rm te}2\tau_{IC}F_0^2}{|B|^2}\mathbf B\cdot\nabla T_e \equiv -\dfrac{0.71}{e}\nabla_\parallel T_e$ | $\sigma_\parallel\dfrac{b_x}{e}0.71\dfrac{\partial T_e}{\partial x}$ | **present, same physics** | Braginskii thermal force. **Same coefficient 0.71**, same sign, same closure. Gated on `thermoelectric_ohm`; coefficient `thermoelectric_coef`. **Also inherits the `tauIC` gate** -- see section 8.2. |

---

## 7. What each code has that the other does not

### 7.1 JOREK: the inductive term $\partial_t\mathbf A$

SOLPS's Ohm's law (B.12) is purely electrostatic; JOREK evolves $\psi$. So JOREK carries
$-\partial_t A_\parallel = -(F_0/R^2|B|)\partial_t\psi$ and SOLPS does not. In a steady SOL this
is small, but every transient -- ELM, sheath-BC relaxation, timestep jump -- drives $E_\parallel$
in JOREK through a channel SOLPS structurally cannot have. Conversely, SOLPS enforces
$\nabla\cdot\mathbf j = 0$ *exactly* at every step; JOREK enforces it only implicitly, through
the vorticity (`var_u`) equation, which is the perpendicular momentum balance whose curl is the
charge-continuity statement. So the two codes determine $\Phi$ by structurally different routes:
**SOLPS solves an elliptic equation for $\Phi$; JOREK gets $\Phi=F_0u$ from a vorticity
evolution equation and gets $\psi$ from this induction equation.** Term-by-term agreement of the
Ohm's-law *terms* therefore does not imply the two codes will produce the same $\Phi$.

### 7.2 JOREK: `keep_current_prof`

Section 5. ~45% of the psi RHS uncut, no SOLPS counterpart, and it corrupts $\eta$-scans.

### 7.3 JOREK: numerical terms

`eta_num_term` (hyper-resistivity) and `zeta_timevol_term` (time scheme). Neither is physics.
`eta_num` values in the $10^{-13}$ range are typical, so the hyper-resistivity is normally
negligible against $\eta\sim10^{-6}$ *unless* $zj$ has grid-scale structure -- which it does near
the X-point and at the targets, exactly where the HFSHD signal lives. Worth measuring rather than
assuming.

### 7.4 SOLPS: everything in $\mathbf j_\perp$

SOLPS writes the *full* current vector: diamagnetic (B.16-17), inertial + gyroviscous
(B.22-23), parallel/perpendicular/heat-driven viscous (B.18-21, B.24-25), and ion-neutral
friction (B.27-33) currents. JOREK's psi equation has none of these -- it is a parallel Ohm's
law and nothing else. The perpendicular current physics lives in the `var_u` equation instead,
in a different (Boussinesq, reduced) form. This is the same gap identified in the temperature
document (section 3.1 there): **the magnetic-drift / $\nabla B$ channel is thin in model600**, and
$\tilde\jmath^{\rm(dia)}$ is precisely the current that closes the $\nabla B$-drift-driven
in-out asymmetry.

### 7.5 SOLPS: anomalous conductivity and anomalous thermoelectric coefficient

`FLAG_SIG`/`PARM_SIG` and `FLAG_ALF`/`PARM_ALF` (manual p. 91). Ad-hoc turbulent additions to
$\sigma_\parallel$ and to the thermoelectric coefficient, with the same Bohm/density scaling
options as the anomalous $D$ and $\chi$. JOREK model600 has no equation-level equivalent --
its only tunable is $\eta_0$ and the $T^{-3/2}$ law.

### 7.6 SOLPS: the current-driven electron heat flux

$q_{e\parallel}^{u} = 0.71\,(T_e/e)\,j_\parallel$ (in $\tilde q_{ex}$, B.40, p. 359). This is
the energy-equation partner of term 7. **JOREK has term 7 but not its partner** -- see the
temperature document section 3.2, where the missing flux was estimated at $\sim23$ kW across the
measured 2.5 kA thermoelectric loop, and, crucially, *antisymmetric between the targets*. The
two come from the same closure and are normally implemented together. On this branch the Ohm's
law half exists and the energy half does not.

### 7.7 Relative sizes of the three physical $E_\parallel$ terms

Order-of-magnitude, for a divertor leg with $L_\parallel\sim20$ m, $T_e\sim10$-30 eV,
$j_\parallel\sim10^5$ A/m$^2$, $\ln\Lambda\sim15$:

| term | estimate |
|---|---|
| $\eta_\parallel j_\parallel$ | $\eta_{\rm Sp}(10\,{\rm eV})\approx2.5\times10^{-5}\,\Omega$m, so $\sim2.5$ V/m |
| $-\frac{1}{en_e}\nabla_\parallel p_e$ | $\sim (T_e/e)\nabla_\parallel\ln p_e \sim 1$ V/m |
| $-\frac{0.71}{e}\nabla_\parallel T_e$ | $0.71\times20\,{\rm V}/20\,{\rm m}\sim 0.7$ V/m |

All three are the same order in the leg. This is the whole point of adding term 7: it is not a
correction, it is a co-equal driver of $\nabla_\parallel\Phi$. The relative importance shifts
strongly with $T_e$: $\eta\propto T_e^{-3/2}$ grows fast in a cold leg, while the thermal force
grows only linearly, so on a 1 eV target `eta_J` dominates -- but $T_e$ there may be set by the
`T_min` clamp rather than by physics, so check `T_min` against the leg temperature first.

---

## 8. Sign and gauge appendix -- the traps

### 8.1 The handedness trap, and the consequence for $\Phi$

**This is the trap that matters.** JOREK stores vectors as `(R, Z, phi)` and treats that basis
as right-handed (section 1.1, four independent confirmations). Doing any cross product in the standard
$(R,\varphi,Z)$ cylindrical convention while reading JOREK components gives the **opposite
sign** for every $\mathbf a\times\mathbf b$, hence for $\mathbf v_{E\times B}$, for
$\nabla\times\mathbf A$, and for the relation between $\Phi$ and $u$.

What is handedness-**independent** (and was previously verified correctly):

- $\mathbf v_{\rm pol} = (-Ru_Z,\,+Ru_R)$ -- these are *components*, `mod_fields.f90:582-583`;
- the density equation's ExB bracket `+v*BigR**2*(r0_s*u0_t - r0_t*u0_s)` = $R[\rho,u] = -\mathbf v_{\rm pol}\cdot\nabla\rho$;
- the psi equation's `+v*(ps0_s*u0_t - ps0_t*u0_s)` = $R[\psi,u] = -\mathbf v_{\rm pol}\cdot\nabla\psi$;
- $\mathrm{div}(\mathbf v_{E\times B}) = -2\partial_Z u$ (used in the temperature document).

What is handedness-**dependent**:

$$\mathbf v_{E\times B} = \frac{\mathbf E\times\mathbf B}{|B|^2}
\;\xrightarrow{\;\mathbf B\simeq \frac{F_0}{R}\mathbf e_\varphi,\;\mathbf E=-\nabla\Phi\;}\;
\left(-\frac{R}{F_0}\Phi_Z,\;+\frac{R}{F_0}\Phi_R\right)
\;\stackrel{!}{=}\;\left(-Ru_Z,\;+Ru_R\right)
\;\Longrightarrow\; \Phi = +F_0 u .$$

Cross-check from the $\partial_\varphi$ term. With $A_\varphi=+\psi/R$,
$E_\varphi = -\frac1R\partial_\varphi\Phi - \frac1R\partial_t\psi$, so
$\partial_t\psi = -R E_\varphi - \partial_\varphi\Phi$. The ideal part gives
$-R E_\varphi = -\mathbf v_{\rm pol}\cdot\nabla\psi = R[\psi,u]$, matching `:1559`, and then
`:1560`'s $-F_0\partial_\varphi u$ must equal $-\partial_\varphi\Phi$, i.e. $\Phi=+F_0u$. Same
answer by a second route.

**Therefore:**

- `diagnostics/new_diag/mod_expression.f90:1758`, `res = u0 * F0 / fact_time  !### sign?` -- the
  comment's doubt is unfounded. **The postproc expression is correct.**
- The convention `Phi = -F0*u` recorded from earlier work is **wrong**, and anything that used it
  to assign a *sign* to $\Phi$ -- the direction of the PFR potential hill, the sign of $E_r$ in
  the sheath-BC diagnostics, the polarity of an imposed bias -- should be re-examined. Anything
  that used it only to get $\mathbf v_{E\times B}$ *components* is unaffected.
- Corroboration 1: the diamagnetic term 4 acquires exactly the Braginskii sign
  under $\Phi=+F_0u$ and exactly the *wrong* sign under $\Phi=-F_0u$. JOREK's diamagnetic
  reduced MHD has been validated for a decade (diamagnetic ELM stabilisation); the reading that
  makes it correct is the right reading.
- Corroboration 2: the measured `min u ~ -1320` against `max u ~ 18` becomes a deep potential
  WELL under $\Phi=+F_0u$ (since $F_0 = 2.97 > 0$), which traps positive ions -- consistent with
  the observed accumulation of N+ in that structure. Under $\Phi=-F_0u$ it would be a +3900 V
  hill, which would expel them.

### 8.2 `tauIC` gates the thermal force

Term 7 is written as `tec * tauIC*2./BB2 * ...` (`:1571-1572`). **If `tauIC = 0` the
thermoelectric term is identically zero regardless of `thermoelectric_ohm` and
`thermoelectric_coef`.** `tauIC = 0` is the code default (`models/preset_parameters.f90:568`).
This is structurally reasonable -- $2\tau_{IC}$ *is* the $1/(eF_0\sqrt{\mu_0\rho_0})$ unit
conversion, not a physics switch -- but it means enabling the thermal force *also* requires
turning on the diamagnetic terms everywhere else in the model (density, vorticity, $T_i$), which
is a much larger change than the flag name suggests. It also means that setting
`tauIC = 1.d-3` against a nominal $2.8\times10^{-3}$ silently scales the thermal force by
$0.36$, so the effective coefficient is $0.71\times0.36\approx0.26$, not $0.71$. **Check
`tauIC` against `tauIC_nominal` in the logfile before quoting a 0.71 result.**

### 8.3 Sign conventions of the ideal term, anchored

The $R[\psi,u]$ sign is anchored, not asserted:
$-\mathbf v_{\rm pol}\cdot\nabla\psi$ with $v_R=-Ru_Z$, $v_Z=+Ru_R$ gives
$R(\psi_R u_Z - \psi_Z u_R) = R[\psi,u]$, which is `+v*(ps0_s*u0_t - ps0_t*u0_s)` exactly.
And the resistive sign is anchored by requiring $\partial_t\psi = \eta\Delta^*\psi$ to be
*diffusive*, which the code's `+` satisfies. Note that the resistive term's sign is
**insensitive** to handedness (both $A_\varphi$ and $j_\varphi$ flip), so it cannot be used to
settle section 8.1 -- only the curl diagnostics and the diamagnetic term can.

### 8.4 The parallel flow does not advect $\psi$

There is no $v_\parallel\partial_\varphi\psi$ term, and there should not be:
$(v_\parallel\mathbf b\times\mathbf B)_\varphi \equiv 0$, so the parallel flow makes no
contribution to $E_\varphi$. **[derived]** If you go looking for a missing advection term, this
is not one.

### 8.5 Gauge

$\Phi$ in this model is defined only up to a function of time; $u$ is fixed by the vorticity
equation and its boundary conditions. Nothing in the psi equation determines the absolute level
of $\Phi$ -- only $\nabla_\parallel\Phi$ enters. Any absolute potential (sheath BC targets,
"potential hill" measurements) is set by the `var_u` boundary conditions, not here. Note that
`sheath_init_u` was measured to be actively harmful (104 steps to 4), see `sheath-gauge-and-tail`.

---

## 9. Things that look wrong (kept separate from the reference above)

1. **`FF_1` is added after the cutoffs.** `models/FFprime.f90:198`. The $\psi_N$ tanh and the
   $Z$ tanh both multiply `prof1` only; `FF_1` is added afterwards, so $FF'$ does not go to zero
   in the SOL or PFR -- it goes to $FF_1$. This is a mechanism for the SOL/PFR
   `keep_current_prof` artefact. The `keep_current_prof_confined` mask in `models/current.f90:88`
   *does* catch it, because it multiplies the fully assembled `zjz`. **If you rely on
   `FFprime`'s own cutoffs to zero the source outside the separatrix, you are wrong by $FF_1$.**

2. **The diamagnetic term divides by `r0_corr`, not by $n_e$.** `:1563-1564` uses
   `Pe0 = (r0 + rimp0*alpha_e)*Te0` -- correctly $n_e T_e$ -- but the prefactor is
   `1/(r0_corr*BB2)` with `r0_corr = corr_neg_dens(r0)` (`:587`), which is the *total/main*
   density, not $n_e$. The exact Braginskii term is $-\nabla_\parallel p_e/(e n_e)$. With
   `with_impurities` and $\alpha_e = (m_i/m_{\rm imp})Z_{\rm imp}-1 \neq 0$ (`:5458`) these
   differ. Without fluid impurities $\alpha_e=0$ and the term is exact. Low priority, but note
   it if impurity runs are compared.

3. **The hyper-resistivity is poloidal-only and $R$-unweighted.** `:1561` has no `/BigR` and no
   `v_p*zj0_p/R**2` piece, unlike every other term. At $n_{\rm tor}=1$ this is nearly moot, but
   it means term 3 is not the toroidally-symmetric analogue of the others and its magnitude
   scales differently with $R$ across the machine. Probably deliberate (it is a stabiliser), but
   it should not be described as "$\eta_{\rm num}\nabla^2 j$" without qualification -- the term
   name in `mod_model_settings.f90:96` says exactly that, and it is misleading.

4. **`keep_current_prof` silently kills the RE coupling.** `:768-772`. Term 6 is unreachable in
   any run with the default `keep_current_prof = .true.`. There is no warning printed. If
   somebody enables `use_rep` expecting the runaway current to feed back into $\psi$, it will
   not, and nothing will say so.

5. **Two resistivities coexist.** `eta_T` (psi equation, clamped at `T_max_eta`) and `eta_T_ohm`
   (Ohmic heating in $T_e$, clamped at `T_max_eta_ohm`), `:1154-1163`. If only one of
   `T_max_eta`/`T_max_eta_ohm` is set, the induction equation and the energy equation see
   different resistivities and the Ohmic power will not equal $\eta j^2$ from the psi equation.
   Both default to `1.d99` (off), so this only bites if someone sets one of them.

6. **The `!### sign?` comment at `mod_expression.f90:1758` should be resolved.** Per section 8.1 the
   expression is correct. Leaving the doubt in place is what allowed the $-F_0u$ convention to
   propagate.

---

## 10. Bearing on the HFSHD campaign

1. **Re-check the sign of $\Phi$** (section 8.1). This is upstream of everything: the sign of the
   measured potential hill, the polarity of the sheath-BC bias scan, and the direction of the
   PFR ExB drift all depend on it. The bias scan was already found to be negative
   (`hfshd-bias-scan-negative`); if the imposed bias had the sign opposite to what was
   intended, that result needs re-reading.

2. **Verify `tauIC` is non-zero and at nominal** (section 8.2), otherwise the `thermoelectric_ohm` work
   is measuring either nothing or a 36%-strength version of the term.

3. **`keep_current_prof` must be masked before any $\eta j_\parallel$ comparison with SOLPS**
   (section 5). At 45% of the RHS it is not a perturbation, and it has no SOLPS counterpart at all.
   Note that the `FFprime` cutoffs alone do not do this -- only `keep_current_prof_confined` with
   `keep_current_confine_strength` $> 0$ does (section 9.1).

4. **The Ohm's-law side of the Braginskii thermoelectric closure now exists; the energy side
   does not** (section 7.6). The measured 2.5 kA loop will now push $\nabla_\parallel\Phi$, but the
   $0.71(T_e/e)j_\parallel$ heat flux that physically accompanies it still will not move heat
   between the targets. Since it is the antisymmetric part that matters for in-out asymmetry,
   having one half of the pair may be worse than having neither for interpreting the asymmetry --
   the current loop will be driven correctly, but its thermal consequence will be absent, so
   $T_e$ at the two targets will not respond.

5. **The whole term-by-term agreement in section 6 does not guarantee agreement in $\Phi$** (section 7.1).
   SOLPS gets $\Phi$ from $\nabla\cdot\mathbf j = 0$ including all the perpendicular currents;
   JOREK gets it from the vorticity equation. If the JOREK and SOLPS potentials differ, the
   parallel Ohm's law is probably not where the difference lives -- look at $\mathbf j_\perp$ and
   at the sheath boundary condition first.
