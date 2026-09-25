# Independent review: sheath current BC (current-row form), model600

Date: 2026-09-24. Reviewer: independent (not the author). Read-only review; no code was changed.

Code read: `floating-u-nodal-sheath-j` (f1f8eeb15, the working tree) and `sheath-j-clean` (15fed7377):
`models/model600/mod_boundary_matrix_open.f90`, `mod_boundary_conditions.f90`, `mod_floating_u.f90`,
`mod_elt_matrix_fft.f90` (induction, vorticity, current and vorticity definitions and their Jacobians),
`models/current.f90`, `models/preset_parameters.f90`, `doc/sheath_j.md`, `sheath-j-clean:doc/sheath_j_plan.md`.
References read: SOLPS-ITER `b2stbc_phys.F` (BCPOT=11, lines 10868-11064), `transport/b2tfch.F` (current
assembly, lines 300-420), SOLEDGE2D (Bufferand et al., NME 12 (2017) 852, eqs. 27-35), the Artola note, Rozhansky
et al., NF 41 (2001) 387. Measured run histories: the campaign notes and items 1-9 of the brief.

Numbers below are order-of-magnitude estimates at Spitzer resistivity and AUG-like parameters (B about 3.3 T at the
inner target, D, JOREK time unit t0 = sqrt(mu0*rho0) = 0.65 us for central_density = 1.01). They are there to show
which term wins. None of them was measured.

---

## Summary

1. **The current-row form has no bug; it has a missing anchor.** At a released wall node, `u` is the Lagrange
   multiplier of the wall half-cell's charge balance. When `f' = 0` (ion saturation, or the electron cap), that
   balance can hold `u` only through (a) the polarisation capacitance, which is tiny and proportional to rho, and
   (b) the response of the interior parallel current to `u` within ONE linear solve. In reduced MHD that response
   is inductive (the current comes from `Delta* psi`, and `psi` moves by induction) unless
   `dt >> mu0 h^2 / eta` at the element scale h. So a node asked for more than `j_sat` gets a potential of order
   `L dI/dt` over one step: kV early in the ramp, and without bound when the excess comes from a current source
   that does not depend on `u`. No choice of slope, cap, weight or ramp inside the row changes this. The slope
   turns "no root" into a root at `Te*(Lambda + excess/s)`, which is kV again (run 5). The alpha ramp makes the
   exponential steeper and produced the period-2 cycle (run 4).
2. **The excess current is mostly not sheath physics.** The code shows four sources that no sheath can carry and
   that the BC does not see: (i) the **frozen equilibrium current on the floating type-9 segment**, which sits
   exactly at the inner strike point where every late crash begins; (ii) `keep_current_prof`, which is
   `.true.` by default and **not masked on either sheath branch** (`keep_current_prof_confined` is not an ancestor
   of either); on open field lines it acts as a current source of strength `j_src(psi)`; (iii) the equilibrium's
   open-field-line current at t = 0; (iv) **implicit boundary fluxes in the retained vorticity row**. The largest of
   these is the divergence-free magnetisation current `R^2 d_s p` from the pressure bracket, which I estimate at
   one to several `j_sat` at hot grazing strike points.
3. **Two SOLPS features are missing and neither is a tunable stabiliser.** First, BCPOT=11 linearises the sheath
   current with the slope `max(j_i, j_e)*e/Te`. That slope is never zero: at ion saturation it is the floating
   conductance, and it stays finite at the electron cap. The residual stays exact. Second, SOLPS's potential
   equation contains the Ohmic parallel conductivity and the ion-neutral conductivity implicitly, so the flux tube
   anchors the wall potential inside the same matrix. JOREK has neither of these at the wall.
4. **What to implement first: zero current on the non-sheath wall types (2, 3, 9)** in place of the frozen
   equilibrium value. It is about 20 lines, adds no parameter, is what SOLEDGE does, and is the only form consistent
   with the floating potential that is already imposed on the same nodes. It targets the measured failure location
   directly. Before that, check `keep_current_prof` in the logs of the failed runs. Second: the SOLPS (Patankar)
   slope. Third: a consistent initial state (floating potential on open field lines, no open-field-line current).
   The Robin form, the total-flux `j_sat` and the ion-neutral Pedersen current are the longer-term structural
   pieces.

---

## 0. Code facts that change the diagnosis

**0.1 `keep_current_prof` is live on open field lines in the sheath binaries unless the namelist switches it off.**
`models/preset_parameters.f90:672` sets `keep_current_prof = .true.`. The induction row, `mod_elt_matrix_fft.f90:1558`,
carries `eta_T*(zj0 - current_source)`, and `current_source` comes from `models/current.f90`, which has no mask on
`floating-u-nodal-sheath-j` or on `sheath-j-clean` (I grepped both with `git show`: no `keep_current_prof_confined`).
The mask exists only on `keep-current-prof-cutoff-clean` (79e2b9f89). The campaign notes record that namelists
setting the mask variables abort on branches without them, so a run that started on these binaries **did not have
the mask**. It had the source unless it set `keep_current_prof = .false.`.

In circuit terms, `eta*(j - j_src)` on an open field line is a Norton source: a current `j_src` in parallel with a
conductance `1/(eta*L)`. As the leg cools, eta grows as `T^-3/2` and the shunt disappears, which leaves a pure current
source of `j_src` driving into the targets, while `j_sat ~ n*sqrt(T)` falls. The PFR matters here: psi_N < 1 again,
so the profiles evaluate to near-separatrix core values, and the PFR borders both strike points. The HFSHD-artefact
analysis (hfshd-eta-dependence) found this term is 45 % of the dominant psi-equation term in the divertor.
**Check this first.** `mod_log_params` prints `keep_current_prof`. If it was on and unmasked, runs 2-6 are not a
test of the sheath BC.

**0.2 Time integration is Gear's (theta = 1, zeta = 0.5)** in the production namelist (`run_with.txt:161`). The code
default is Crank-Nicolson (`preset_parameters.f90:14`), which would turn every stiff wall mode into a
period-2 oscillation. That is not what happened here. The period-2 cycle of run 4 is a Newton overshoot (A.4).

**0.3 The zj wall trace is an L2 projection, not a pointwise characteristic.** The row is
`sum_g w_g*Zbig*dl*v(g)*[zj_h(g) - j_sat(g) f(x(g))]` tested with the zj trace functions
(`mod_boundary_matrix_open.f90:412-437, 477, 559-570`): a Galerkin projection of `j_sat*f(u)` onto the Hermite
trace (value plus tangential derivative per node), evaluated at several Gauss points per edge. Where `j_sat` varies
by a factor 2 across one element (strike point, junction), the cubic projection overshoots. A pointwise
`zj/j_sat` of 1.1-1.4 (run 2, "1.2-1.4 at the junction for 20 steps") is therefore not in itself proof of
an unsatisfiable demand. It is a known artefact of this form.

---

## A. Why the current-row form runs away at j > j_sat

### A.1 The wall-node system as coded

At a sheath node whose `u` and `zj` Dirichlet rows are released (`mod_boundary_conditions.f90:240-273, 399-403`), with
psi and w Dirichlet:

- **zj row** (Zbig): `zj_h = P[j_sat f(x)]`, `x = Lambda - a_n (u - C_V V_w)/(2 Te)`, with f = 1 - exp(min(x, Lambda))
  (- s min(x,0)). Linearised column on u: `J_u = j_sat f'(x) dx/du = j_sat (a_n/2Te) e^x`. This is zero at the
  cap (x >= Lambda) and exponentially small on the ion side.
- **u row** = the vorticity equation tested with the wall basis function `phi_w`
  (`mod_elt_matrix_fft.f90:1579-1616`, Jacobian `:2486` polarisation, `:2534` parallel current). Integrating the
  parallel-current term by parts over the wall half-cell,
  `int phi_w [psi, zj] = - int zj [psi, phi_w] + oint_wall phi_w zj d_s psi ds`,
  gives the discrete charge balance of the half-cell:

      C_w (1+zeta) du/dt  =  I_in(zj_i)  +  I_impl  -  I_sh(zj_w)                                  (A1)

  Here `C_w = int R^3 rho_corr |grad phi_w|^2` is the polarisation capacitance (proportional to rho), `I_in` is the
  parallel current delivered from the interior, `I_sh = oint phi_w zj_h d_s psi` is the current leaving through the
  wall (= the sheath current by the zj row), and `I_impl` collects the implicit boundary fluxes of every
  by-parts/bracket term of the row (A.3).

One linear solve per step gives

    [ C_w(1+zeta)/dt + theta*G_p + theta*G_sh ] du  =  I_in + I_impl - I_sh(u0)                    (A2)

with `G_sh = oint phi_w^2 J_u |d_s psi|` (the sheath conductance, which goes to 0 at saturation) and `G_p` the change
of `I_in` per unit `du` within the same solve.

### A.2 What G_p is in reduced MHD, and why it is usually too small

`I_in` depends on `u_w` only through induction: the psi rows of the interior nodes contain `theta*dt*[psi, du]`
(`:2407`) and `theta*dt*eta*zj` (`:2411`), and Ampere (`:1652, 2739`) returns `dzj = Delta* dpsi`. For a current
structure of width h near the wall:

    dzj_i ~ - theta dt k^2 (grad_par du) / (1 + theta dt eta k^2),   k ~ 1/h                       (A3)

- Ohmic limit (`dt >> tau_h = mu0 h^2/eta`): `G_p ~ grad_par/eta`, a finite loop conductance. A root exists at
  `dPhi = excess/G_p`, and it is bounded if the excess comes from an EMF.
- Inductive limit (`dt << tau_h`): `G_p ~ theta dt grad_par^2/h^2`, small. The interior current is a stiff current
  source over one step.

At Spitzer with h = 2 mm: `tau_h ~ 2.8 us = 4 t0` at 50 eV, and `8 ns = 0.012 t0` at 1 eV. During the ramp
(dt = 1e-3 ... 2) hot targets are therefore in the inductive regime and cold targets in the Ohmic one.

In the inductive regime, the voltage the solve must produce to change an element-scale current by `dj` within one
step is the inductive EMF:

    V ~ mu0 h^2 dj l_par / dt  ~  7.7 kV * (1e-3/dt)      [h = 2 mm, dj = 1e5 A/m^2 ~ 2 j_sat at 50 eV, l_par = 10 m]

That is the order of runs 2-5: kV at t ~ 2, hundreds of volts at dt ~ 0.1-0.3. **The solver is computing
L dI/dt correctly. What is wrong is the demand.** A physical SOL plasma never carries more than `j_sat` into a
sheath, so a state that does is an inconsistent initial or forcing condition, not a transient the BC should
tolerate.

With `f' = 0`, (A2) leaves `u_w` to `C_w/dt + theta G_p`. Where the inflow sink has depleted rho, `C_w` is proportional to
`rho_corr`, so capacitance and `j_sat` vanish together. That is why the collapse finishes in about 3 steps.

### A.3 Where the excess comes from (sources the BC does not see)

1. **Frozen current on non-sheath walls.** Every non-released zj DOF is Dirichlet in increment form (matrix entry
   only, no RHS: `mod_boundary_conditions.f90:405-408`), so it stays at its t = 0 value for ever. On type 9 that
   value is the equilibrium current density at the inner strike point, noted as `1e5-1e6 A/m^2`, i.e. 2-20 j_sat
   in run 2. It reaches the junction node through three channels: the wall trace of the type-9-side edge in the
   junction's u row, the interior Ohm's law `eta*zj` (Galerkin overlap), and the zj-definition mass matrix.
   It is a constant current source. See section C.
2. **`keep_current_prof`** (0.1): a Norton current source of strength `j_src` on open field lines.
3. **The equilibrium's open-field-line current at t = 0** (FF'/p' tails, PFR): about 1 Isat into types 1/4 at
   step 1 (run 4 `[wall I]`), 3.4 j_sat locally. It relaxes by induction in about 180 steps. Until then it is an
   inductive current source (A.2).
4. **Implicit boundary fluxes of the retained vorticity row (I_impl).** With Dirichlet `u` they acted on overwritten
   rows and were harmless. With `u` released they are wall currents:
   - pressure bracket `+ R^2 (v_s p_t - v_t p_s)` (`:1588`). By parts,
     `int R^2[v,p] = - int v [R^2,p] + oint v R^2 d_s p`. The volume part is the curvature divergence of the
     diamagnetic current (correct). The boundary part is the normal component of the divergence-free
     **magnetisation current** `grad x (p b/B)`, about `d_s p / B` for b ~ e_phi. Its ratio to the sheath
     capacity is

         j_M.n / (j_sat |b.n|)  =  (1 + Ti/Te) rho_s / (L_s |b.n|),     L_s = p/|d_s p| along the wall

     At 50 eV, rho_s ~ 0.44 mm; with L_s = 5 mm and |b.n| = 0.05 this is about 3.5, and about 17 for L_s = 1 mm (the
     measured 50 eV/mm fronts). At 1 eV it is about 0.5-2.5. It is a dipole about the pressure peak: ion-side
     excess on one flank, electron-side on the other. B2.5 cancels exactly these divergence-free parts
     analytically (Rozhansky 2001, sec. 2.2: "effective" diamagnetic velocity and current); SOLEDGE and SOLPS never
     let them into the wall balance.
   - ExB vorticity advection `- r0_hat R^2 w0 (v_s u0_t - v_t u0_s)` (`:1580`), which gives
     `oint v rho R^2 w d_s u`: polarisation charge carried through the wall by vE.n (1e4-1e5 m/s at the failures).
   - viscosity `- visco_T R^2 grad v . grad w0` (`:1583`), which gives `oint v visco R^2 d_n w`, driven by the
     difference between the **frozen Dirichlet w at the wall** and the evolving interior w. It is another frozen
     value acting as a source.
   - smaller terms: `0.5 vv2 [v, rho]` (`:1579`), tauIC terms (`:1598-1602`), hyperviscosity.

   None of these is bounded by `j_sat`, and none depends on the sheath potential in a way that anchors it.

The sheath-j-clean plan (Option III, "the retained vorticity row's dropped surface terms ARE the wanted physics")
checked only the polarisation and viscous terms. The pressure bracket was listed under "physics refinements", but
with u released it is not a refinement: it is part of the wall charge balance.

### A.4 Can any choice within the form be well-posed in a one-solve code?

Each step solves a nonsingular linear system (`C_w > 0`), so it is always solvable. The problem is conditioning,
and whether the nonlinear system has a root at all. The nonlinear wall balance
`I_in + I_impl = I_sh(u) in [-19 j_sat, +j_sat]` has a root only if the delivered current lies inside the range
of the characteristic. The row cannot supply what the plasma side lacks:

| choice | effect | verdict |
|---|---|---|
| ion slope s | root at Phi = Te(Lambda + excess/s): 83 Te = 3.4 kV for excess 2.4, s = 0.03 (run 5) | an s-dependent fake root; a fit |
| e_slope past the cap | same on the electron side | a fit |
| alpha ramp y = x/alpha | makes the exponential steeper: a 3 V step, zero-slope cap 10 V from floating; the one-solve Newton overshoots onto the cap, then back (run 4, period 2) | wrong homotopy |
| Zbig weight | the row is a constraint; the weight changes nothing | irrelevant |
| potential row u = (2Te/a_n)(Lambda - ln(1 - j/j_sat)) | the degeneracy moves to the zj column (dPhi/dj -> inf at j -> j_sat); measured: rails at 2600 V | same ill-posedness |
| **Jacobian slope floored at the floating conductance (SOLPS)** | residual exact, steady state unchanged; step bounded by Te*excess/j_sat; the cap keeps slope e^Lambda | the only in-form change that is not a fit (F2) |
| **range from the total ion flux (SOLPS fchi = e*fna)** | more ion capacity where ExB flows out; consistent with the rho/T rows and recycling | consistency, not a crash cure (F5) |

Run 4's period-2 cycle, step by step: from x0 = -3 the coded slope is e^-3 = 0.05. For a target f = 0.5 the linear
step lands at x = -3 + 0.45/0.05 = +6, i.e. on the cap. At the cap the u column is exactly 0, so the vorticity row
alone places u, and the next solve throws it back. With SOLPS's slope `max(e^x, 1)` the step from the ion side is
conservative (slope 1 >= the true slope on (x0, 0)). From the electron side, Newton on the convex exponential is
monotone. The cap keeps the slope e^Lambda ~ 20. The cycle cannot form.

**Verdict on A.** The current-row form is the right continuum statement: charge balance with the sheath current as
the wall flux, the same as SOLEDGE. It is ill-posed in JOREK because the demanded current contains sources that are
not EMF-limited, and because the loop that would throttle the current is inductive at the element scale during the
ramp. No slope, cap, weight or ramp fixes that. Fix the demand (C, E, F4), and put the SOLPS slope in the Jacobian
so that each step stays bounded while the plasma responds.

---

## B. How SOLPS-ITER and SOLEDGE make it well-posed, and the Robin form in JOREK variables

### B.1 SOLPS-ITER, BCPOT = 11 (`b2stbc_phys.F:10868-11064`)

- The sheath enters the **potential equation** as a charge source in the boundary cell:
  `sch0 = [(1-gamma_see) fche(phi_old) - fchi + t0*phi_old] - t0*phi_new`.
  The residual is exact at phi_old. The linearisation slope is
  **`t0 = max(fchi, (1-gamma_see) fche) * e / Te`**. At ion saturation that is `fchi*e/Te`, the floating
  conductance, not the zero derivative. At the electron cap `fche` uses `expu2(max(-50, min(0, -e(phi-phi_w)/Te)))`,
  so the current is capped but `t0 = fche*e/Te` stays finite. This is Patankar's rule of a non-positive source
  slope, and it is the default definition of BCPOT = 11. It is not the optional `stab_coeff_sheath_*`
  (`:2111, 6835, 8942`), which the user rightly excludes.
- `fchi = e * sum Z * fna`: the **total** ion particle flux, the same flux as the BCCON = 14 particle BC. The
  consistency check at `:589-641` requires `istyle_fchi = 0` with the drift-compatible set. The electron current
  `fche = e n sqrt(Te/(2 pi m_e)) S exp(-e dphi/Te)` does not depend on the ion flux.
- The potential equation assembles the parallel current as the **Ohmic** `fch = csig*ehx` (`b2tfch.F:307-321`, ehx
  including the pressure and thermal-force terms) and adds the ion-neutral conductivity `csigin` to `conc`
  (`:333-343, 383-392`) plus viscous and inertial currents. The whole flux tube's conductance, parallel and
  perpendicular, is therefore implicit in the same matrix as the sheath.
- B2.5 iterates its equations nonlinearly within a step.

### B.2 SOLEDGE2D (Bufferand 2017, eqs. 27-34)

- `j_par = -sigma_par (grad_par phi - grad_par p_e/(n e) - 0.71 grad_par Te/e)` (eq. 29) is substituted into the
  vorticity equation, which is then solved **for phi** as an anisotropic Laplacian
  `(Delta_perp + sigma_par Delta_par) phi` (eq. 34). The sheath (eq. 31) is a **Robin** condition on that operator.
- The paper states the operator "is invertible as long as boundary conditions do not degenerate to Neumann". At
  saturation they do degenerate locally. The problem stays solvable because `sigma_par Delta_par` connects the
  saturated end to the other end of the tube, and `Delta_perp` plus the vorticity closure `zeta grad_perp Omega`
  (eq. 33) connect it to the neighbouring tubes.
- Lambda follows the actual ion flux (eq. 32). "On any other boundary surface ... or any tangential wall, currents
  are forced to zero."

### B.3 Why they do not run away, and what JOREK lacks

In both codes the current into a saturated sheath is throttled inside the same linear solve by the Ohmic
conductance of the flux tube (electrostatic `j_par`, no `dA_par/dt`), and SOLPS never linearises with a zero slope.
JOREK has an electromagnetic `j_par` (Ampere plus induction), so the throttle exists only when `dt >> mu0 h^2/eta`
(A.2). It has no ion-neutral conductivity in the vorticity row for kinetic neutrals: the only friction-like term is
the fluid-neutral `Sion_T*rn0` term at `:1611`, and `aux_rho0` enters rho, vpar and T but not u. It also uses the
zero derivative as its slope. Transplanting the Robin BC alone does not fix saturation.

### B.4 The Robin (conductance) form in JOREK variables

In the continuum this is the same condition as the current-row form. It changes which row carries the sheath and
removes the projected zj trace.

Rows at sheath edges (Gamma_s = edges with both endpoints sheath types and |b.n| >= sin(theta_min), as now):

- **u row** (vorticity, released per DOF as now): replace the implicit wall outflow of the parallel-current term,
  `+ dt oint v zj d_s psi ds`, by the sheath current:

      RHS_u(v) += dt * sigma_o * oint_Gs v * [ j_sh(u, rho, T) - zj ] * d_s psi0 ds                 (B1)
      j_sh = j_sat(rho, T) * f(x),   j_sat = c_sat rho cs sign(B.n)/|B|   (or the total-flux form, F5)

  Here sigma_o is the orientation sign that makes `oint v zj d_s psi` the outflow of `int v [psi, zj]` (the same
  orientation as `sf_orient`, `mod_boundary_matrix_open.f90:380`). Jacobian:

      A(u,u)   = - theta dt sigma_o oint v * j_sat * (df/du)_P * phi * d_s psi0 ds,
                   (df/du)_P = -(a_n/2Te) * max(e^min(x,Lambda), 1)        [SOLPS slope, residual exact]
      A(u,zj)  = + theta dt sigma_o oint v * phi * d_s psi0 ds             [removes the implicit zj outflow]
      A(u,rho), A(u,Te/Ti) = exact, through j_sat and x (as the present zj-row columns)
      psi column: none needed. psi is Dirichlet on the wall trace (value and along-wall derivative,
                  mod_boundary_conditions.f90:389-397), so d_s dpsi = 0 exactly and B_pol.n is fixed, not lagged.

- **zj row at the wall**: back to Ampere, *with* the surface term,
  `int v zj/R + int grad v . grad psi/R - oint v d_n psi/R dl = 0`. This existed on sheath-j-clean as part of the
  potential-row route and behaved as a current measurement. The zj Dirichlet row is not released at all, which
  removes half of the per-DOF release logic.
- psi: Dirichlet (unchanged). w: Dirichlet (unchanged), or its definition row with `+ oint v R d_n u` (F4).

What it fixes:
1. The sheath current enters at each Gauss point with `|j_sh| <= j_sat(g)`, so there is no projection overshoot
   (0.3).
2. zj stays `Delta* psi` in the last element, so the interior Ohm's law no longer sees a forced trace.
3. The frozen wall-trace problem disappears for sheath nodes.
4. The SOLPS slope is simply the Robin coefficient.

What it does not fix: with an excess from a non-EMF source and an inductive loop, (A2) is unchanged. It is the
right long-term form, not the cure.

---

## C. The type-4/type-9 junction at the inner strike point

**A frozen nonzero current next to released nodes is a defect.** It is also inconsistent with the potential imposed
on the same nodes. `u = C_T Te + C_V V_w` on type 9 is the zero-current point of the characteristic, and `zj = zj_eq`
contradicts it: the user's rule 2 fails at the segment level. Three consequences:

- It is a **constant current source** at exactly the point where every late crash starts (runs 2, 3, 6). The
  junction node's u row integrates the type-9-side edge's wall current, and that trace interpolates to the frozen
  value.
- Its ratio to the local `j_sat` **grows as the target cools**: `j_sat ~ n cs` falls while `zj_eq` is fixed. This
  matches run 6 exactly: j/j_sat 0.36 -> 0.8 -> 1.3 -> 2.1 -> 6.8 over 8 steps during detachment, starting at the
  junction.
- It is the only return path for the sheath types' net current (run 2 `[wall I]`: all of 1/4/5 collect net electron
  current, returned through the pinned 2/3/9). Global charge conservation is violated by construction.

The campaign notes dismissed `zj_wall_zero` because the start-up hypothesis was falsified. That argument concerns
the switch-on. It does not touch run 6, where the mechanism above applies.

**Zero current (a floating segment) is right**, as SOLEDGE does ("currents are forced to zero" on every non-sheath
boundary). In the floor regime `|b.n| < sin(theta_min)` the model has no sheath physics to give anything else, and
zero is the value consistent with the floating potential imposed there.

**Treatment of short, degenerate-frame segments without making them sheath types:**
- zj: Dirichlet with target 0 on the value DOF and the along-wall derivative DOF. This is the same pattern as the
  floating row's `fu_target`, i.e. `RHS = -Zbig*(zj_old - 0)` at `mod_boundary_conditions.f90:405`. The same applies
  to types 2 and 3.
- u: keep the floating Dirichlet row (u released on type 9 died in 5-19 steps: rung 1 and "type 9 as sheath type").
  Dirichlet rows do not care about frame quality; the volume vorticity row on those elements does.
- The junction node keeps its per-DOF release, as now.
- Accept, and measure, the one remaining inconsistency. Pinned-u nodes do not enforce charge balance, so the
  interior parallel current arriving at type 9 is absorbed. With zj = 0 on all non-sheath walls,
  `sum over sheath types of I_net` in `[wall I]` is exactly this leak, which gives a direct diagnostic.
- Expect a potential step of order `Te*ln(1 - j/j_sat)` over one element between a current-carrying junction and
  the floating segment. That is the physics of a current-carrying region's edge, not an error.

---

## D. The inflow sink under a sheath potential

**Is it a genuine positive feedback?** Half of the loop is physical and half is the A defect.
- Physical half: where vE.n < 0, the wall absorbs and never emits, and the volume ExB carries the wall cell's content
  inward. The cell depletes on `h/|vE.n| ~ 2 mm / 3e4 m/s ~ 0.1 t0`, i.e. within one step. So `j_sat ~ n cs`
  drops. SOLPS and SOLEDGE have the same drift-driven rarefaction.
- Numerical half: a depleted node with the same delivered current saturates. `u` then runs away (A.2), with
  `C_w ~ rho` falling at the same time. The resulting potential hill drives ExB circulation that strengthens the
  inflow on one flank, which gives more depletion. The loop gain is infinite only because `dPhi/dI` at saturation is
  unbounded.

So the sink is the trigger (run 2 died at t = 15 with it, run 1 at t = 109 without it). It is not the defect. Removing
it brings back the inflow ill-posedness (rho < 0 at the inflow point, runs 1 and 482), so it must stay.

**The consistent form (SOLPS):** one flux `Gamma_n = n max(v_n, v_fl)` (with the drift bounded to +-2 cs|b.n| as
on floating-u-bounded, if the Vpar row is bounded the same way) in:
- the rho, Ti and Te rows (now);
- the recycling (now);
- the **sheath ion current**, `j_i,n = e Gamma_n`, i.e. in zj units `j_sat = c_sat rho max(v_n, v_fl)/(|b.n| |B|)`.

The electron current stays separate: `j_e,n = e n sqrt(Te/2 pi m_e) |b.n| exp(-e dPhi/Te)`, capped at dPhi = 0. On
the Bohm branch this reduces exactly to the present form with `Lambda = 0.5 ln(m_i/(2 pi m_e (1+Ti/Te)))`. On the
ExB-outflow flank it raises the ion capacity by `1 + vE.n/v_fl`, which is up to 10x at the measured drift/cs. The
present option (a) violates rule 2: the rho row sees the ExB outflow but the sheath's ion current does not.
SOLPS does not run away on the same n because of B.3, not because of the flux form.

---

## E. Switch-on from t = 0

The alpha ramp failed because it steepened the characteristic. A blend "pinned zj -> sheath" starts in the
rung-1 configuration (u released, zj pinned), which died in 14-19 steps. **Do not ramp the characteristic. Make the
initial state satisfy the final equations**, and use the SOLPS slope so that any residual mismatch relaxes at a
bounded rate:

1. **Initial wall potential at floating.** The Sept-21 log starts at Phi ~ -1 V (Dirichlet u = 0 baseline), i.e.
   x = Lambda. Every sheath node starts on the cap, where the u column is zero. That is the -19 j_sat of step 2.
   Initialise `u = C_T Te + C_V V_w` on open field lines (smoothly masked to the core value across the separatrix)
   and set `w = Delta u` consistently. Every sheath node then starts at x ~ 0 (j ~ 0, slope `j_sat a_n/2Te`), the
   best-conditioned point of the characteristic. This changes only `initial_conditions.f90`.
2. **No open-field-line current at t = 0.** Confine `keep_current_prof` (merge the mask) and remove the equilibrium's
   SOL/PFR current, either by masking FF'/p' in the GS source outside the separatrix or by a relaxation baseline
   with zj = 0 walls (F1) and floating u. The Sept-21 chain was built on a cutoff baseline.
3. **The SOLPS slope (F2).** An excess `dj` then moves the potential by at most `(Te/e)*dj/j_sat` per step, while
   induction relaxes the current. The same equations apply at step 1 and step 10^4, with no ramp parameter.

If a homotopy is still wanted, ramp the source of the mismatch (the initial open-field-line current), not the
stiffness.

---

## F. Options, ranked

**Before anything (0.1 day):** check whether `keep_current_prof` was `.true.` in runs 2-6 (the logfile prints it)
and whether the binary had `keep_current_prof_confined`. If it was on and unmasked, re-run the run-6 case with it
off or masked before judging any BC.

| # | option | rows / columns touched | effort | P(fixes the observed crash class) | risk | diagnostic within ~50 steps |
|---|---|---|---|---|---|---|
| **F1** | **zj = 0 on non-sheath wall types (2, 3, 9)**, u floating as now | `mod_boundary_conditions.f90:399-437`: zj Dirichlet RHS target 0 on value and along-wall DOFs | 0.5 d | **high** for the late junction crashes (runs 3 and 6); low for early hot-target ones | step-1 jump where zj_eq is large (hot type 9); lower if combined with E | `[wall I]`: sum over sheath types ~ 0 (the leak); junction j/jsat max < 1 while the inner target cools past t ~ 1000; Inet/Isat of opposite signs on hot and cold targets (thermoelectric) instead of all negative |
| **F2** | **SOLPS (Patankar) slope** in the sheath row's u column: `J_u = j_sat (a_n/2Te) max(e^min(x,Lambda), 1)`, residual exact | `mod_boundary_matrix_open.f90:428-435` (`sc_dfdu`), harness FD test of u replaced by a slope test | 0.5-1 d | medium alone (turns 3-step blow-ups into a bounded drift of Te*excess/j_sat per step); high combined with F1/E | formally a lagged-Jacobian term, zero at steady state; coefficient = physical floating conductance, no free number. **The user must decide whether this falls under the stabiliser ban**; in my view it does not (it is the definition of BCPOT = 11, same class as the existing lagged Btot/sign(B.n) columns) | per-step max `|dPhi|` at sheath nodes <= Te*|j/jsat - 1| (+ margin); no alternation of the e-sat fraction; run-4 namelist (alpha ramp) no longer shows period 2 |
| **F3** | **Consistent initial state**: u = C_T Te on open field lines, w consistent, open-field-line current removed, `keep_current_prof` confined | `initial_conditions.f90`, merge 79e2b9f89 (`current.f90`), equilibrium FF'/p' mask | 1-2 d | high for the switch-on (runs 4/5 type) | touches the equilibrium path; needs regression on a run without a sheath | step-1 `[sheath_j]`: e-sat 0, `|Inet/Isat| << 1`, Phi within [0, 2 Lambda Te]; `[wall I]` flat through the ramp |
| **F4** | **Remove implicit non-sheath currents** from the wall charge balance: magnetisation current (pressure bracket), ExB vorticity advection, viscous flux from frozen w | first a diagnostic in `mod_wall_diag.f90` (per node `oint phi_w R^2 d_s p`, `oint phi_w rho R^2 w d_s u`, `oint phi_w visco R^2 d_n w` vs `oint phi_w j_sat |d_s psi|`); then compensating u-row boundary integrals on sheath edges in `mod_boundary_matrix_open.f90` (u rows only, like fx/ex) and the w definition with `+oint v R d_n u` | diag 1 d; fix 2-3 d | medium (magnitude predicted O(1-10) j_sat at hot grazing strike points; unmeasured) | sign conventions: FD harness mandatory | diag: any term > 0.3 I_sat at a strike point confirms; after the fix the electron and ion flanks of the strike-point dipole at t < 20 fall inside the characteristic |
| F5 | **j_sat on the total ion flux** (same Gamma as the rho/T rows and recycling) with the SOLPS electron term and Lambda(Ti/Te) | sheath block of `mod_boundary_matrix_open.f90:404-437` plus columns (vpar, u through vE.n, rho, T); floating row Lambda consistency | 1-2 d | low for the crash; required by rule 2 | the ExB column makes the row depend on u_s (new sparsity entry) | saturated fraction on the ExB-outflow flank drops; Inet/Isat recomputed with the total-flux Isat |
| F6 | **Robin form** (B.4): sheath current as the wall flux of the u row, zj back to Ampere with its surface term | u-row boundary integral in `mod_boundary_matrix_open.f90`; zj definition surface term (exists on sheath-j-clean's potential-row route); drop the zj release in `mod_boundary_conditions.f90` | 3-5 d | medium (removes projection overshoot and the frozen/forced trace; not the saturation cure) | the structural change is contained; the u release is unchanged | pointwise `|j_sh| <= j_sat` by construction; Tier-1 check: I_sheath = I_Ampere per type |
| F7 | **Ion-neutral (Pedersen) current** in the vorticity row for kinetic neutrals (ionisation from aux_rho0; CX needs a projected CX rate) | u row volume term like `:1611` with aux sources; aux projection for CX | 3-5 d | low-medium for the crash (gives a saturated node a cross-field path, estimated root at ~10^2 V instead of kV); **high for HFSHD fidelity** (SOLPS has csigin/fchin) | new physics, needs validation | Phi at a saturated node stays finite and depends on n_0 as predicted |
| - | slopes s / e_slope, alpha ramp, potential row, wall Ohm's law (psi swapped into zj), types 4/9 released, type 9 as sheath, SOLPS stab_coeff, use_sc | - | - | measured failures or fits | - | do not retry |

**On option (ii) of the brief**, bounding the current the interior pushes into a node. The only physical bound is
the loop impedance `R + L/dt` seen by the EMFs that drive the current. An algebraic bound on `zj_i` in the last
element would break Ampere and charge conservation, so it is a clip. The legitimate version is: remove the Norton
sources (0.1, F1, F4), start consistently (F3), and let the loop act (Ohmic when `dt >> mu0 h^2/eta`) with
bounded steps (F2).

**On option (iii), the potential-row form revisited.** It is the same degeneracy in the other column
(`dPhi/dj -> inf` at `j -> j_sat`). It is viable only where the delivered current is already EMF-limited, and there
the current-row form works too. Not recommended.

**From SOLPS/SOLEDGE, not yet tried:** the Patankar slope (F2); ion current from the total flux with a separate
electron term and flux-dependent Lambda (F5, SOLEDGE eq. 32); zero current on all non-sheath boundaries (F1,
SOLEDGE); ion-neutral conductivity in the potential equation (F7, SOLPS csigin/fchin); analytic cancellation of the
divergence-free magnetisation current in the charge balance (F4, Rozhansky 2001).

### The one to implement first: F1

- Every late crash (runs 2, 3 and 6) starts at the type-4/type-9 junction on the inner strike point. The frozen
  type-9 current is the one source there that is (a) fixed in time, (b) large (2-20 j_sat in the run-2 estimate),
  and (c) growing relative to the local `j_sat` as the target cools, which is run 6's signature.
- It is required by the user's own consistency rule (floating potential and frozen current contradict each other)
  and matches SOLEDGE.
- It is about 20 lines with no parameter, and it produces a clean discriminator: the run-6 setup (sheath-j-clean +
  recycling sign fix, Sept-21 namelist, `keep_current_prof` verified) with only zj = 0 on 2/3/9 changed.
  - If the inner target now cools through t ~ 1050 without junction saturation, the BC is usable for the HFSHD
    test, and F2/F3 become robustness work.
  - If it still saturates at the junction, the excess comes from the interior (0.1 or F4). The F4 diagnostic,
    which is cheap and can run in the same build, tells which.

Run the F4 diagnostic in the same build: it only prints.

---

## Verdict

The current-row form is the correct continuum statement of the sheath: charge balance with the sheath current as
the wall flux, the same as SOLEDGE. It fails in JOREK for two reasons that no parameter of the characteristic can
touch. First, the wall nodes are asked to carry currents that no sheath can carry and that the BC does not see:
the frozen equilibrium current of the floating type-9 segment at the inner strike point, a `keep_current_prof`
current source that is unmasked on both sheath branches unless the namelist disabled it, the equilibrium's
open-field-line current at t = 0, and the implicit magnetisation, vorticity-advection and viscous fluxes of the
released vorticity row. Second, in reduced MHD the flux tube's ability to throttle that current in response to the
wall potential within one linear solve is inductive at the element scale during the ramp, and the row's own
conductance is linearised to zero at both ends. The potential is then set by the polarisation capacitance and
L dI/dt, which gives kV in a few steps. The slope, the ramp and the potential-row variants move this singularity
around without removing it. The fixes are: remove the non-sheath current sources (zero current on non-sheath
walls first, then a masked current source, a consistent initial state, and the magnetisation current out of the
wall balance), take SOLPS's non-zero sheath linearisation slope, and in the longer term move to the Robin form with
`j_sat` from the total ion flux and add the ion-neutral Pedersen current that SOLPS has and model600 lacks for
kinetic neutrals. The physical HFSHD states need tens of volts, bounded by the thermoelectric EMF, not kV. A BC that
produces kV is being fed a current that the physics would never deliver.
