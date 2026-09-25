# Why a sheath current boundary condition is hard in JOREK

*A note for SOLPS-ITER and SOLEDGE users, 2026-09-25.*

**Scope.** This note explains why a sheath current-voltage (j-V) boundary condition that works routinely in
SOLPS-ITER (B2.5, `BCPOT = 11`) and SOLEDGE2D/3X is structurally harder in JOREK (model600: reduced MHD with a
kinetic particle coupling), and what would make it work. It stands alone, but it draws on an earlier independent
review of the same problem (`doc/review_sheath_j_2026-09-24.md`, cited below as "the review") and on the measured run
histories of the 2026 sheath campaign.

**Code references.** `file:line` references are to branch `sheath-j-clean` (commit 10131b184, the branch that
carries the sheath row with the weak Mach row and the zero-current option) unless another branch is named. The volume
equations in `mod_elt_matrix_fft.f90` are the same on `floating-u-nodal-sheath-j` (f1f8eeb15). SOLPS references are
to `SOLPS-ITER/modules/B2.5/src`.

**Numbers.** Estimates are order of magnitude. They assume Spitzer resistivity with ln Lambda_C = 10, deuterium,
B about 3 T in the divertor, and the JOREK time unit t0 = sqrt(mu0 rho0) = 0.65 us (central density 1.0e20 m^-3).
Statements marked *measured* come from the campaign's run logs.

---

## 1. What a sheath current boundary condition is

### 1.1 The characteristic

A plasma in contact with a conducting wall at potential V_w develops a thin (a few Debye lengths) non-neutral sheath.
Let Phi be the plasma potential at the sheath entrance (the magnetic presheath edge), n and Te the density and
electron temperature there, c_s = sqrt((Te + Ti)/m_i) the sound speed, and let j be the parallel current density
**into the wall**, counted positive for a net ion current. The Bohm criterion fixes the ion flux into the sheath at
n c_s. Electrons are Boltzmann-repelled by the sheath drop, so their flux is the one-sided thermal flux times
exp(-e(Phi - V_w)/Te). Together:

    j(Phi) = j_sat * [ 1 - exp( Lambda - e(Phi - V_w)/Te ) ],        j_sat = e n c_s               (1)
    Lambda = ln( sqrt(m_i / (2 pi m_e)) * sqrt(Te/(Te + Ti)) )  ~ 2.8 for D with Ti = Te (3 in the runs)

Everything is per unit area normal to B. The flux through a wall element carries the extra factor |b.n|, where b is
the field direction and n the wall normal. Three points on the curve matter:

- **Floating point**, j = 0: Phi_f = V_w + Lambda Te/e, about 3 Te/e. No net current flows.
- **Ion saturation branch**, Phi - V_w >> Lambda Te/e: all electrons are reflected and j -> +j_sat. Raising Phi
  further does not raise the current. The plasma can deliver at most j_sat of positive current, because the ion
  flux is fixed by the Bohm condition and not by the potential. (In a real Langmuir probe the ion branch keeps a
  small slope from sheath expansion. That slope is geometric, weak, and not a physics model for a divertor target.)
- **Electron saturation**, Phi <= V_w: the sheath no longer repels electrons and j -> -(e^Lambda - 1) j_sat,
  about -19 j_sat for Lambda = 3. This is the cap in the code.

The slope of the characteristic, the **sheath conductance** per unit area, is

    g_sh = dj/dPhi = (e/Te) * j_sat * exp(Lambda - e(Phi - V_w)/Te) = (e/Te) * (j_sat - j)            (2)

It equals j_sat e/Te at floating and falls exponentially on the ion branch. It is **zero at ion saturation**, and
zero again when the electron current is capped.

Worked numbers:

| | 50 eV, n = 1e19 m^-3 | 1 eV, n = 9e19 m^-3 |
|---|---|---|
| c_s | 6.9e4 m/s | 9.8e3 m/s |
| j_sat | 1.1e5 A/m^2 | 1.4e5 A/m^2 |
| g_sh at floating | 2.2e3 S/m^2 | 1.4e5 S/m^2 |

### 1.2 Why Phi is a free variable

Phi is not a boundary datum. Only the wall potential V_w is imposed. Phi is whatever potential makes the charge
balance of the plasma consistent, div j = 0, with the current leaving through the sheath given by (1). The sheath
condition is therefore a **relation between the boundary value of Phi and the normal current**: a nonlinear Robin
condition, not a Dirichlet or Neumann condition. Two limits are common simplifications:

- **Floating** (Dirichlet): Phi = V_w + Lambda Te/e. Every flux tube carries zero current at each end.
- **Insulating** (Neumann): j.n = 0.

Neither lets current close through the wall. The full characteristic lets a flux tube carry current from one target
to the other, bounded by what the two sheaths and the tube resistance allow.

"The plasma can only deliver up to j_sat" has a sharp consequence. In a steady state with no external circuit, a flux
tube can carry at most j_sat into the wall at the positive (ion) end, and each sheath sits at a potential of a few
Te. A state that pushes more than j_sat into a sheath is not a strongly biased state. It has no steady solution at
all, because the potential has no value at which the wall accepts the current. The only physical currents that come
near j_sat are the ones driven by EMFs of order Te/e (the thermoelectric effect, the parallel pressure gradient, and
the 0.71 grad_par Te thermal force), which are limited by the same sheaths. A well-posed model must therefore keep
every current that reaches a sheath inside the range [-19 j_sat, +j_sat] by construction.

### 1.3 Why we want it: the thermoelectric mechanism of the HFS high-density region

Senichenkov et al. (Contrib. Plasma Phys. 2022, `hfshd_currents.pdf`, building on Rozhansky) explain the AUG
high-field-side high-density front (HFSHD) through the thermoelectric current. The hot, attached outer target
floats at a higher sheath potential (about 3 Te,out) than the cold, detached inner one. A current therefore flows
in the SOL from the hot target to the cold one, through the highly resistive cold inner leg. At 1 eV,
eta_par ~ 5e-4 Ohm m and j ~ 1e4 A/m^2 (the measured AUG outer-target magnitude), so E_par ~ 5 V/m, or 5-25 V over
1-5 m of detached leg. The resulting potential hill near the X-point drives an E x B flow across the private flux
region from the outer to the inner divertor, which feeds the HFSHD. A floating BC forbids this current by
construction, so the mechanism needs a sheath j-V condition on both targets. The potentials involved are tens of
volts. A boundary condition that produces kilovolts is not modelling this physics.

---

## 2. How SOLPS-ITER and SOLEDGE make it well-posed

### 2.1 The potential equation is elliptic, and the current is algebraic in Phi

B2.5 solves for the potential `po` from charge continuity, div j = 0, in each cell (`transport/b2tfch.F`). The
parallel current is **Ohm's law**, an algebraic function of the potential gradient:

    j_par = sigma_par * ( -grad_par Phi + grad_par p_e/(e n) + 0.71 grad_par Te/e )                   (3)

In b2tfch this is `fch = csig*ehx + fchdia + fchin + fchvispar + fchvisper` (`b2tfch.F:307-321`). `ehx` contains the
potential, pressure and thermal-force terms of (3). `csig` is the Spitzer parallel conductivity. The perpendicular
currents are added as extra fluxes: diamagnetic (`fchdia`, in Rozhansky's "effective" form with the divergence-free
magnetisation part removed analytically), ion-neutral friction or Pedersen (`fchin`, with the conductivity `csigin`
added to the implicit coefficient `conc`, `b2tfch.F:333-343`), viscous and inertial. Every one of these currents
responds to Phi **inside the same linear system**. The current is not a state variable. It is recomputed from Phi
whenever Phi changes, and it has no memory.

SOLEDGE does the same in a different form (Bufferand et al., NME 12 (2017) 852, eqs. 27-34). It substitutes Ohm's law
(eq. 29) into the vorticity equation, which is charge continuity with the polarisation current as the time
derivative. It then solves that equation for phi as a strongly anisotropic 2D Laplacian:

    ( Delta_perp + sigma_par Delta_par ) phi^{t+dt} = Delta_perp phi^t - A(phi^t) + RHS                (SOLEDGE eq. 34)

Delta_perp carries the polarisation capacitance (density- and metric-weighted). sigma_par Delta_par is the implicit
Ohmic parallel conductance.

### 2.2 The sheath is a boundary conductance

In both codes the sheath is simply the boundary flux of this elliptic problem.

- **SOLPS `BCPOT = 11`** (`sources/b2stbc_phys.F:10868-11064`) adds the sheath current to the boundary cell as a
  charge source `sch0`. The ion current `fchi = e * sum Z Gamma_i` is the ion particle flux of the same boundary
  cell. With `istyle_fchi = 0` this is the **total** flux `fna`, the one the particle BC `BCCON = 14` sets
  (`:10884-10899`). A consistency check requires that choice with the drift-compatible set
  `BCCON = 14 / BCMOM = 13 / BCENE = 15` (`:589-641`). The electron current is
  `fche = e n v_te/sqrt(2 pi) * exp(max(-50, min(0, -e(po - V_w)/Te)))` (`:10929-10933`). The `min(0, ...)` is the
  electron cap. The source is linearised as

      sch0 = [ (1-gamma_see) fche(po_old) - fchi + t0 * po_old ] - t0 * po_new
      t0   = max( fchi, (1-gamma_see) fche ) * e / Te                                     (b2stbc_phys.F:10953-10961)

  The residual is exact at the old iterate: the t0 terms cancel at convergence. The linearisation slope t0 is
  **never zero**. On the ion side it is fchi e/Te, which is at least the true slope fche e/Te, so the step is
  conservative. At the electron cap fche is constant, so the true derivative is zero, but t0 = fche e/Te stays
  finite. This is Patankar's rule for source terms (a linearised source S = S_C + S_P Phi must have S_P <= 0,
  which keeps the matrix diagonally dominant), applied to the sheath. It is the definition of `BCPOT = 11`, not the
  optional `stab_coeff_sheath_*` stabiliser.
- **SOLEDGE** writes the sheath current (eq. 31, j = e Sum Z Gamma_par exp(Lambda - (phi - phi_w)/Te), with Lambda
  from the actual ion flux, eq. 32) as a **Robin** condition on the Laplacian. The paper states that the operator
  is invertible "as long as boundary conditions do not degenerate to Neumann", and that "on any other boundary
  surface ... or any tangential wall, currents are forced to zero". Only a sheath lets current in or out.
- **Bohm-Chodura on the total flux.** In both codes the particle, momentum and energy BCs act on the total normal
  flow, parallel plus E x B plus other drifts (SOLPS: `BCMOM = 13` imposes the Bohm inequality using the interior
  u_par, with the drift contribution bounded; SOLEDGE eq. 35). The ion current in the sheath BC is that same flux.
  The sheath current, particle BC and energy BC therefore all use one ion flux.
- **Nonlinear iteration.** B2.5 iterates its coupled equations within a time step and is normally run to a steady
  state, where any linearisation error, including the lagged t0, vanishes. SOLEDGE re-solves the potential after
  each transport update. Neither code relies on a single linearised solve to land on the right branch of an
  exponential.

### 2.3 The 1D picture: a resistor chain between two sheaths

Take one flux tube of length L between target A (s = 0) and target B (s = L), both walls grounded (V_w = 0). Let j
be the current density along +s, so j enters wall B and -j enters wall A. Integrating (3) along the tube:

    j = G (Phi_A - Phi_B + E),     G = sigma_par / L,     E = int ( grad_par p_e/(e n) + 0.71 grad_par Te/e ) ds   (4)

The ends satisfy j = I_B(Phi_B) and -j = I_A(Phi_A), with I from (1). Eliminating j gives the two-sheath (Harbour)
condition I_A(Phi_A) + I_B(Phi_B) = 0 and one Ohmic drop. The thermoelectric current follows because
I_A and I_B have different Te and j_sat.

Discretised in N cells this is a resistor chain. Node potentials Phi_i are joined by conductances
G_{i+1/2} = sigma_par/Delta s, and each end node has a conductance to ground g_A, g_B (the sheath slopes). One
Newton step solves a tridiagonal system:

    -G_{i-1/2} dPhi_{i-1} + (G_{i-1/2} + G_{i+1/2}) dPhi_i - G_{i+1/2} dPhi_{i+1} = r_i
    end B:  (G_{N-1/2} + g_B) dPhi_N - G_{N-1/2} dPhi_{N-1} = r_N                                      (5)

Now let end B saturate, g_B -> 0. The matrix stays nonsingular. Node B hangs on the chain, and the effective
conductance from B to ground is

    g_eff(B) = g_B + [ 1/G_tube + 1/g_A ]^-1,     G_tube = sigma_par/L                                  (6)

This stays positive as long as end A is not saturated too. **The potential at a saturated end is anchored through
the flux tube to the other end.** Numbers: at 10 eV and L = 10 m, G_tube = sigma_par/L ~ 6e3 S/m^2, comparable to
g_sh ~ 5e3 S/m^2 at floating for n = 1e19. A saturated node with an excess demand of 1 j_sat moves by about
j_sat/g_eff ~ 10-20 V, a few Te. Only if both ends saturate does the 1D problem degenerate to Neumann. In 2D the
perpendicular currents (polarisation, Pedersen, anomalous) still connect the node to its neighbours, and the
Patankar slope keeps even the local block regular. In steady state the demand never exceeds the range in the first
place, because every current in (3) is itself limited by the same sheaths.

Three ingredients make this well-posed:
1. an implicit, memoryless Ohmic conductance along the tube;
2. a sheath slope that is never zero in the Jacobian;
3. no current enters the charge balance that is not either Ohmic (so it responds to Phi) or a sheath current. Every
   other boundary is insulating.

JOREK lacks all three by construction.

---

## 3. What JOREK solves instead

### 3.1 Variables and equations (reduced MHD, model600)

JOREK does not have a potential equation in the B2.5 sense. It evolves, on a 2D grid of bicubic Hermite finite
elements (each node carries a value plus derivatives), times toroidal Fourier harmonics:

- **psi**: poloidal magnetic flux. It is proportional to R A_phi, the toroidal vector potential, which is
  A_par to leading order.
- **u**: the E x B stream function, with **Phi = F0 u** (F0 = R0 B0). The E x B velocity is v_E = -R^2 grad u x grad phi.
- **zj**: the toroidal current density variable, **zj = Delta* psi** (Ampere's law, a separate equation and unknown,
  `mod_elt_matrix_fft.f90:1652`). The parallel current is essentially zj/R (b is close to e_phi). The current into
  the wall is proportional to zj (B_pol.n).
- **w**: the vorticity, **w = Delta_pol u** (`:1658`).
- rho, T (or Ti, Te) and v_par as in any fluid code.

The two equations that replace B2.5's potential equation are:

- **Induction** (`:1558`, Jacobian `:2407-2411`), which is Ohm's law with inductance:

      d psi/dt = eta (zj - j_src)/R + R [psi, u] - F0 d u/d phi + (p_e and thermal-force terms)            (7)

  Here [a, b] is the poloidal Poisson bracket, so [psi, u] is B_pol.grad u. In physical terms this is
  E_par = -grad_par Phi - d A_par/dt = eta j_par - grad_par p_e/(e n) - ..., i.e. (3) **plus the inductive term
  dA_par/dt**. j_src is the optional `keep_current_prof` source, which must be masked to the confined region.
- **Vorticity** (`:1579-1616`, Jacobian `:2486`, `:2534`), which is charge continuity, div j = 0:

      d/dt (rho R^2 w) - (polarisation terms)  =  [psi, zj] - F0 d zj/d phi  +  R^2 [p, R^2-stuff] + viscosity + ...   (8)

  [psi, zj] is B.grad(j_par), the parallel current divergence. The pressure bracket (`:1588`) is the diamagnetic
  and curvature current divergence. The left side is the polarisation current. Its Jacobian is the capacitance
  term -R^3 rho grad v . grad u (1 + zeta) (`:2486`).

In B2.5 terms: where B2.5 has `fch = csig * ehx(po)`, an algebraic function of po, JOREK has `[psi, zj]`, where zj
is an **independent unknown**. zj is tied to psi by Ampere, and psi is tied to u only through the induction equation,
i.e. through time. The polarisation (inertial) current, a small correction in B2.5, is the time derivative that
makes JOREK's (8) an evolution equation for u.

### 3.2 Time stepping: one linear solve per step

JOREK advances all variables together with a theta scheme (Crank-Nicolson by default; Gear's BDF2, theta = 1 and
zeta = 0.5, in the production runs). The nonlinear system is linearised about the old state once, and **one** sparse
linear system is solved for the increment dX = X^{n+1} - X^n. There is no Newton loop and no inner iteration. Each
step is exactly one Newton step from the previous state. A strongly nonlinear boundary relation such as (1)
therefore gets exactly one linearisation per time step, evaluated at the old state.

### 3.3 The parallel current is a state variable with inductance

Consider one discrete "resistor" of the chain (5) in JOREK: a current channel of perpendicular width h (at least
the element size) and parallel length l. Per unit cross-section it has resistance R = eta l and inductance
L ~ mu0 h^2 l (the vector potential of a channel of width h is A ~ mu0 h^2 j, up to a logarithm). One implicit step
of (7) gives

    j^{n+1} = [ (L/dt) j^n + theta (dPhi + E) ] / ( L/dt + theta R )                                    (9)

Each resistor of SOLPS's chain becomes, over one step, two elements:

- a **conductance** G_eff = theta / (L/dt + theta R), and
- in parallel with it, a **current source** (L/dt) j^n / (L/dt + theta R), which carries the old current forward
  whatever the potential does.

The crossover is the magnetic diffusion time at the channel scale:

    tau_h = L/R = mu0 h^2 / eta                                                                         (10)

- **dt >> tau_h** (Ohmic): G_eff -> 1/(eta l) and the current source vanishes. This is SOLPS's chain.
- **dt << tau_h** (inductive): G_eff -> theta dt / L, small, and the current source is simply j^n. **The tube
  delivers its old current regardless of the wall potential.**

Numbers for h = 2 mm (the element scale at the strike points):

| Te | eta_par | tau_h = mu0 h^2/eta | in t0 units |
|---|---|---|---|
| 50 eV | 1.5e-6 Ohm m | ~3 us | ~4 t0 |
| 10 eV | 1.6e-5 Ohm m | ~0.3 us | ~0.5 t0 |
| 1 eV | 5e-4 Ohm m | ~10 ns | ~0.015 t0 |

The JOREK runs start with dt = 1e-3 t0 and ramp to 0.3-10 t0. **Hot targets are therefore in the inductive regime
throughout the start-up, and cold (detached) targets are Ohmic.** For the whole tube (l = 10 m, h = 2 mm,
L ~ 5e-11 H m^2 per unit area), the inductive conductance theta dt/L is about 13 S/m^2 at dt = 1e-3 t0 and about
4e3 S/m^2 at dt = 0.3 t0. The Ohmic value 1/(eta l) is about 7e4 S/m^2 at 50 eV. Moving an element-scale current by
dj = 1e5 A/m^2 (about 1 j_sat at 50 eV) within one step takes an EMF of

    V ~ L dj/dt ~ mu0 h^2 l dj / dt ~ 7.7 kV x (1e-3 t0 / dt)                                           (11)

This is the voltage of an inductive kick. An electrical engineer would recognise the situation: a current-carrying
inductor whose circuit is suddenly opened.

### 3.4 The wall-node balance

At a sheath node the Dirichlet rows for u and zj are released (section 4). u then comes from the vorticity equation
tested with the node's basis function phi_w over the wall half-cell. Integrating [psi, zj] by parts turns (8) into the
discrete charge balance of that half-cell:

    C_w (1 + zeta) du/dt = I_in(zj_interior) + I_impl - I_sh(u)                                         (12)

- C_w = int R^3 rho |grad phi_w|^2 is the **polarisation capacitance**. It is proportional to the mass density rho.
- I_in is the parallel current delivered from the interior.
- I_sh is the current leaving through the wall, which the zj row sets equal to the sheath current (1).
- I_impl collects the boundary fluxes implied by the other terms of (8) (section 4b).

One linear solve gives

    [ C_w (1+zeta)/dt + theta G_p + theta G_sh ] du = I_in + I_impl - I_sh(u^n)                          (13)

G_sh is the sheath conductance (2), evaluated at the old state. G_p is the response of I_in to du **within the
same solve**. By (9), G_p is the inductive G_eff in the start-up regime.

Now let the node saturate, as it does whenever the demand I_in + I_impl exceeds j_sat:

- **G_sh = 0.** The row's u column is j_sat f'(x) (df/du), which vanishes on the ion side as e^x and is exactly zero
  at the electron cap (`mod_boundary_matrix_open.f90:458-477`, column `:646`).
- **G_p is inductive**: small at small dt. The current it throttles is the tube current, not the excess.
- **C_w/dt** is the only anchor left. It is proportional to rho, and rho falls exactly where the target cools or E x B
  empties the wall cell. The polarisation capacitance of a 2 mm wall cell is small. The RC time of the cell against
  the sheath conductance is of order 1e-2 t0, so C_w/dt is small compared with g_sh at every production dt.
- The excess on the right does not depend on u at all when it comes from a frozen current, from the inductive memory
  j^n of the tube, or from I_impl.

The potential then moves by du ~ excess/(C_w/dt + theta G_p) per step. It moves again at the next step, because the
excess is still there. It stops only when the inductive current j^n has changed, which by (11) takes kV-scale EMFs
unless dt >> tau_h. If the excess is a true current source (does not depend on u) and exceeds the combined ion
capacity of every path it can reach (this sheath plus the far sheath through the tube plus cross-field paths), the
continuous problem has **no root**. The potential then rises without bound, like the open-switch inductor.

The feedback that finishes the run is physical. A potential of 1 kV across 2 mm at 3 T is an E x B drift of
5e5 V/m / 3 T ~ 1.7e5 m/s. The campaign logs show |vE.n| ~ 1.7e5 m/s and d Phi/ds ~ 5e5 V/m one step before the
crash. That drift empties the wall cell within a step (h/vE ~ 10 ns), so rho, j_sat and C_w all fall together, and
rho goes negative.

**Measured.** With correct recycling, every run with the sheath row died at the inner strike point, next to the short
floating segment (boundary type 9) whose Dirichlet current is frozen at the equilibrium value. As the target cooled
during detachment (j_sat ~ n sqrt(Te) falling), the demanded current walked through j/j_sat = 0.36, 0.8, 1.3, 2.1,
6.8, 390 in about 8 steps. The potential then ran to kV in 2-3 steps. The same recycling with the floating potential
alone (u = 3 Te/e, no current) ran stably into a detached divertor at 0.4 eV and 9e19 m^-3. The transport part of
the model can reach the detached state. The sheath current row cannot survive the approach to it.

### 3.5 What changes, in one sentence each

- **Electromagnetic instead of Ohmic parallel current.** The tube's conductance becomes an R-L branch. For
  dt << mu0 h^2/eta it is a weak conductance in parallel with a current source that carries the old current, so a
  saturated sheath is anchored through the tube only in the Ohmic limit. SOLPS and SOLEDGE are always in that limit,
  because their j_par has no dA_par/dt.
- **One linear solve instead of an inner iteration.** The linearisation of an exponential at the old state is the
  whole step. With a zero slope at saturation and at the cap, a single step can overshoot from one branch onto the
  other, or produce an unbounded du. B2.5 avoids both with the Patankar slope and with iteration to a steady state.

---

## 4. JOREK-specific complications

**(a) Dirichlet rows are penalties on the increment, so a "Dirichlet current" is a permanent current source.**
JOREK imposes a Dirichlet condition by overwriting the diagonal of that DOF's row with Zbig = 1e12
(`mod_assembly.f90:23-60`, which sets `a_mat%val = ZBIG`; `mod_boundary_conditions.f90:164, 410-412`). The matrix
acts on the increment dX, so the row reads Zbig dX = RHS. When no RHS is added (`boundary_conditions_add_RHS`,
`mod_assembly.f90:77-99`), dX = 0 and the variable stays **at its t = 0 value for ever**. For u on a floating wall a
target is supplied (u = C_T Te + C_V V_w, `mod_boundary_conditions.f90:430-440`). For zj on every non-sheath wall
type, until the zero-current option below, none is supplied: the wall current is frozen at the equilibrium current
density of the Grad-Shafranov solution. At the inner strike point (the floating segment of type 9) that is
1e5-1e6 A/m^2, about 2-20 j_sat, and it is fixed while j_sat falls as the target cools. In B2.5 language, it is as
if a boundary face carried a fixed `fch` imposed for ever, while its neighbouring faces obey `BCPOT = 11`. It also
contradicts the floating potential imposed on the same nodes, which is the zero-current point of the characteristic.
SOLEDGE forces these currents to zero. The option `bcs(i)%zj_zero` (`mod_boundary_conditions.f90:414-425`, commit
10131b184) supplies the RHS -Zbig zj^n, which drives every pinned zj trace DOF to zero in one solve.

**(b) Galerkin advection without integration by parts gives implicit wall fluxes.** JOREK tests brackets such as
[psi, zj], [p, R^2] and [u, w] with the basis function v and does not integrate them by parts. The boundary flux is
therefore not a prescribed quantity. It is whatever the discrete interior solution implies at the wall: a "natural"
outflow. For rho and T this gave an ill-posed inflow problem at points where vE.n points into the plasma. No inflow
datum exists there, and rho went negative at the maximum E x B inflow point. The fix was a sheath-set wall flux,
Gamma = n max(v_n, v_fl) (`fx_n`, `ex_n`, `mod_boundary_matrix_open.f90:364-398` on `floating-u-nodal-sheath-j`,
rows `:483-510`). For the vorticity row the consequence is new. As long as u was Dirichlet, the implicit boundary
fluxes of (8) sat in overwritten rows and did nothing. Once u is released they are wall currents in (12), and the
sheath BC does not see them:
- the **pressure bracket** R^2 [v, p] (`:1588`). By parts it gives the boundary integral of R^2 d_s p, the normal
  component of the divergence-free magnetisation current curl(p b/B). Its ratio to the sheath capacity is about
  (1 + Ti/Te) rho_s / (L_s |b.n|), where rho_s is the ion sound gyroradius and L_s the along-wall pressure scale.
  That ratio is about 3.5 at 50 eV with L_s = 5 mm and |b.n| = 0.05, and about 17 at the measured 1 mm fronts
  (estimate from the review, not measured). B2.5 removes exactly this divergence-free part analytically (Rozhansky
  et al., NF 41 (2001) 387, "effective" diamagnetic current).
- **E x B advection of vorticity** (`:1580`): polarisation charge carried through the wall by vE.n.
- **viscosity** (`:1583`): a flux driven by the difference between the **frozen Dirichlet w** at the wall and the
  evolving interior w. This is another frozen value acting as a source.

None of these is bounded by j_sat, and none depends on the sheath potential in a way that anchors it.

**(c) The wall current row is an L2 projection on cubic Hermite traces, not a pointwise condition.** The sheath row
is sum_g w_g Zbig dl v(g) [zj_h(g) - j_sat(g) f(x(g))], accumulated over the Gauss points of each wall edge
(`mod_boundary_matrix_open.f90:428-477`, columns `:641-651`). Each wall node carries a value DOF and an along-wall
derivative DOF per harmonic. The characteristic is therefore imposed as a Galerkin projection of j_sat f(u) onto a
cubic trace, not cell by cell as in SOLPS's finite volumes. Where j_sat changes by a factor 2 across one element (a
strike point, a type junction), the projection overshoots. A pointwise j/j_sat of 1.1-1.4 is then an artefact of the
form, not proof of an impossible demand. Which DOFs are released also takes care. The value DOF is released if any
incident wall edge is a sheath edge. A derivative DOF is released only if the sheath edge runs in its direction
(`mod_boundary_conditions.f90:260-273, 404-408`, `node_incidence` `:860-897`). Otherwise a corner node would have its
flux-surface derivative freed by a wall edge in the other direction.

**(d) Drift consistency, and grazing wall segments with degenerate element frames.** SOLPS's drift-compatible set
applies the Bohm inequality, the particle and energy fluxes and the sheath ion current to **one** total normal flux
that includes vE.n. JOREK's Mach row pins v_par (to c_s, or with a separate drift term), while the rho and T rows see
the total flow. On the nodal branch the sheath's j_sat uses only the parallel Bohm flow
(`j_sat = c_sat rho v_fl/|B|`, `:418` there). The ion current in the charge balance and the ion flux in the
particle balance therefore differ by the E x B part, which is large at grazing incidence, where
|b.n| ~ 0.01-0.05 and vE.n ~ c_s |b.n| is easy to reach. The geometry adds a JOREK-specific problem. B2.5's targets
are the poloidal ends of a flux-aligned grid. JOREK's grid is flux-aligned inside and extended to the real wall by
ray-cast elements. Parts of that wall (types 4, 5, 9) meet the field at grazing angles through elements with
near-degenerate frames (Jacobian determinant ~0.05 measured at the type-4/9 corner). The sheath row is therefore
applied only where |b.n| >= sin(theta_min) (`mod_boundary_matrix_open.f90:428`). Type 9 stays floating. Releasing
types 4/9 wholesale died in 5-19 steps, and type 9 as a sheath type died in 15 (measured). The consequence is that
short floating segments sit next to current-carrying ones, exactly where (a) bites.

**(e) The switch-on starts from a state that violates the final equations.** A B2.5 run has no current to "relax"
at t = 0. j_par is recomputed from po at the first iteration. A JOREK run starts from a Grad-Shafranov equilibrium
whose FF' and p' profiles carry current on open field lines (the SOL and the private flux region). That current is a
state variable (psi, zj) and can decay only inductively. *Measured*: about 1 I_sat of electron current flowed into
the targets at step 1 (3.4 j_sat locally) and took about 180 steps to relax. In addition the wall potential starts at
the baseline u = 0 (Phi ~ -1 V, below floating). Every sheath node therefore starts on the electron cap, where its u
column is exactly zero. *Measured*: step 1 had e-saturation 1.0 on all sheath types, step 2 had j/j_sat = -24 and
Phi up to 510 V. The obvious homotopy, a switch-on ramp of the characteristic y = x/alpha with alpha ramped from
0.05 to 1 (a sheath at alpha Te), makes things worse. A small alpha makes the characteristic **steeper**: the whole
transition from floating to the cap spans a few volts, and the cap, with its zero slope, sits about 10 V from
floating. One linear step from the ion side, where the slope is e^x, overshoots onto the cap. The next solve throws u
back. *Measured*: a period-2 flip-flop of the e-saturated fraction (0.43 / 0.005 / 0.38 / 0.000) and of Phi_min
(-290 / +90 / -175 / +126 V) every step, and a crash earlier than without the ramp. The homotopy that works in a
one-solve code relaxes the **source of the mismatch** (the initial open-field-line current and the initial wall
potential), not the stiffness of the boundary relation.

---

## 5. Why the naive fixes fail

Each item below changes the sheath row's shape or the solver's damping. None of them changes the three structural
facts of section 3: the demand contains currents that do not depend on u, the tube responds inductively, and the
one-step Jacobian has zero slope at saturation.

**A finite ion-saturation slope** (f = 1 - e^x - s min(x, 0), s = 0.03). This gives the row a root on the ion side
at Phi = Te(Lambda + excess/s), where excess is (j - j_sat)/j_sat. For excess = 2.4 that is 83 Te, about 3.4 kV at a
hot target. *Measured*: Phi_max = 3437 V, then an E x B drift of 1e5 m/s and death. At a cold target the same slope
is a few volts per j_sat. *Measured*: an identical walk at the inner strike point (j/j_sat 0.38, 1.7, 1.9, 8.9, 88,
then 8e4) and a crash at step 668 instead of 671. The slope turns "no root" into a root whose location is set by a
number with no physics behind it. It is a fit.

**An electron-side slope beyond the cap** (e_slope). This is the same construction on the other branch, with the
same verdict.

**Caps.** Capping the electron current at x = Lambda (Phi = V_w) is physical for the value, but the code also uses
the capped derivative: the u column is exactly 0. In a one-solve code that is the landing zone of every overshoot
(section 4e). SOLPS caps the value and keeps the slope (t0 = fche e/Te).

**Ramps of the characteristic.** A ramp in alpha steepens the exponential (4e). A ramp that blends "pinned zj" into
"sheath" starts in the configuration with u released and zj pinned, which died in 14-19 steps. Both are homotopies
in the wrong parameter.

**Potential-row variants.** One can invert (1) and impose u = (2Te/a_n)(Lambda - ln(1 - j/j_sat)) in the u row,
with zj free. That moves the singularity from the u column to the zj column (dPhi/dj -> infinity at j -> j_sat).
*Measured*: the potential ran to the rails at 2600 V. With zj pinned it is floating with extra steps, since no
current can flow. With psi swapped into the zj row (a local wall Ohm's law), it gave +-1000 j_sat at steady state,
because the wall zj was the only slack variable of the last element's parallel-field mismatch.

**The penalty weight Zbig.** The zj row is a constraint. Its weight does not change the solution.

**Stabilisers** (shock capturing `use_sc`, Taylor-Galerkin `tg_num`, SOLPS's `stab_coeff_sheath_*`). They add
numerical dissipation in the volume. They do not provide a current path that a saturated node lacks, and they do
not remove a u-independent current source. They can delay the runaway by smoothing the gradients that the E x B
feedback amplifies, which hides the ill-posedness rather than removing it. The campaign excludes them as a rule.

**Releasing more wall types.** Releasing types 4/9 wholesale, or making type 9 a sheath type, puts the vorticity
row on degenerate-frame elements (4d). *Measured*: death in 5-19 steps.

---

## 6. What would make it well-posed in JOREK

The review's ranked options, in plain language, together with their SOLPS or SOLEDGE counterparts:

**0. Mask the current source.** `keep_current_prof` is `.true.` by default (`models/preset_parameters.f90:672`) and
enters the induction equation as eta (zj - j_src) (`mod_elt_matrix_fft.f90:1558`). On an open field line this is a
current source j_src in parallel with a conductance that disappears as the leg cools. It must be confined to closed
field lines. The campaign confirmed it was masked in all measured runs, so it is not the cause there, but it has to
stay masked.

**F1. Zero current on non-sheath wall segments** (types 2, 3, 9), with u floating as now. This is SOLEDGE's rule,
and the only choice consistent with the floating potential already imposed on those nodes. It removes the frozen
equilibrium current at the inner strike point: the one source there that is fixed in time, large, and growing
relative to j_sat as the target cools, which is the measured signature. Implemented as `bcs(i)%zj_zero`
(`mod_boundary_conditions.f90:414-425`, 10131b184). *Measured (2026-09-25):* the discriminating run, the
configuration that died at step 671 with only zj_zero on 2/3/9 changed, died at step 672 with the identical junction
walk. The frozen current is therefore NOT the driver there. The same log showed min Te at the junction falling
2.4 -> 0.7 -> 0.3 eV in three steps while the integrated type-4 current stayed at 0.13-0.22 I_sat: j_sat collapsed
under a roughly constant demand. F1 stays as the consistent choice, but it is not the cure. The original plan: Success means the inner target cools through t ~ 1000 with the
junction j/j_sat < 1, and the hot and cold targets carry net currents of opposite sign (the thermoelectric
signature) instead of all collecting electrons. A remaining imbalance, the sum of the net currents over the sheath
types, then directly measures what the pinned-u segments absorb.

**F2. The Patankar slope in the Jacobian.** Keep the residual exact and replace the u column's j_sat f'(x) by
j_sat (a_n/2Te) max(e^min(x, Lambda), 1): at least the floating conductance on the ion side, and e^Lambda at the cap.
This is exactly `BCPOT = 11`'s t0. It changes no converged state. It bounds each step's potential change to about
Te x (excess/j_sat), and it removes the overshoot mechanism behind the period-2 cycle. It is a lagged-Jacobian term,
the same class as the lagged B and sign(B.n) columns already in the row, and its coefficient is a physical
conductance with no free number. Whether it counts as a "stabiliser" under the campaign's rules is the user's call.
The review argues that it does not. *Status:* implemented as `sheath_j_patankar` together with the `[wall J]`
print diagnostic of F4 (commit 276b91f40), run once and removed again (2026-09-25): the floored slope left the row
unsatisfied over 17-21% of the inner target (j/j_sat 1.6-2.4 where the residual allows at most 1) and the run died
at step 608 against 672 without it. The diagnostic stays.

**F3. A consistent initial state.** Start with u at the floating value on open field lines (smoothly masked across
the separatrix) and w = Delta u consistent with it, so that every sheath node starts at x ~ 0, the best-conditioned
point of the curve. Remove the equilibrium's open-field-line current, either by masking FF'/p' outside the separatrix
in the Grad-Shafranov source, or by a relaxation run with zero-current walls and floating u. If a homotopy is still
wanted, ramp this mismatch, not the stiffness.

**F4. Remove the implicit non-sheath currents from the wall balance.** First add a print diagnostic, per node, of the
boundary integrals of the magnetisation current (R^2 d_s p), the E x B vorticity advection and the viscous flux,
compared with the sheath capacity j_sat |B_pol.n|. If any of them is above ~0.3 I_sat at a strike point, subtract it
on the sheath edges with compensating boundary integrals in the u rows, and give w its definition row with the
surface term instead of a frozen Dirichlet value. This is the analogue of Rozhansky's cancellation of the
divergence-free diamagnetic current in B2.5. *Status (2026-09-25):* the diagnostic found the viscous flux at
4-20 j_sat and the magnetisation current at 0.6-12 j_sat at the strike points, both rising with the Te front; the
compensating u-row terms are implemented (`sheath_j_cancel_flux`, default on) with w still Dirichlet, i.e. the
weak `d_n w = 0` of the drift-fluid codes. Releasing w needs a second condition (a Neumann penalty on its normal
derivative), otherwise the fourth-order u-w problem is one boundary condition short; that is the next step.

**F5. j_sat from the total ion flux, with a separate electron term.** Use the same Gamma = n max(v_n, v_fl) that the
rho and T rows and the recycling use, and the electron current e n v_te/sqrt(2 pi) |b.n| exp(-e dPhi/Te) with
Lambda(Ti/Te) (SOLPS `istyle_fchi = 0`, SOLEDGE eqs. 31-32). This is a consistency requirement (one ion flux
everywhere), not a crash cure.

**F6. The Robin form.** Take the sheath current out of the zj row and put it into the u row as the wall flux of the
vorticity equation: replace the implicit outflow of [psi, zj] on sheath edges by the sheath current, with the F2
slope as the Robin coefficient, and restore zj = Delta* psi with its surface term at the wall. This is SOLEDGE's
structure in JOREK variables. It removes the projection overshoot (4c), since |j_sh| <= j_sat holds at every Gauss
point, and it removes the forced zj trace. It does not by itself cure saturation, because (13) is unchanged. 3-5 days
of work.

**F7. The ion-neutral (Pedersen) current.** SOLPS includes `csigin`/`fchin` in the potential equation. model600's
vorticity row has a friction-like term only for fluid neutrals (`:1611`), none for the kinetic neutrals used in these
runs. Adding it gives a saturated node a cross-field current path whose conductance grows with the neutral density,
exactly in the detached leg. This matters for HFSHD fidelity, and it gives section 2.3's "cross-field anchor" in 2D.

### 6.1 The deeper option: an Ohmic parallel conductance at the wall

Everything above removes the unphysical demand and bounds the steps. None of it restores the ingredient that makes
SOLPS robust even when the demand is momentarily wrong: an implicit, memoryless Ohmic conductance along the flux tube,
coupled to the sheath in the same matrix. There are three ways to get it in JOREK.

1. **Take dt >> mu0 h^2/eta.** Then (9) is Ohmic. This holds already at cold targets (tau_h ~ 0.01 t0 at 1 eV), which
   is why cold targets behave better. At hot targets it requires dt >> 4 t0 during the start-up, which the other
   physics does not allow. It is also not a formulation fix, since it holds only for a range of dt.
2. **An electrostatic wall layer.** In the last element row next to a sheath edge, replace the inductive Ohm's law
   by the electrostatic one: substitute j_par = sigma_par(-grad_par Phi + ...) for [psi, zj] in the wall half-cell's
   charge balance, as SOLEDGE does everywhere. This gives the sheath its anchor. The cost is structural. In that layer
   zj is no longer Delta* psi, so Ampere and magnetic energy conservation are broken there, and a matching condition
   is needed between the layer's electrostatic current and the interior electromagnetic one. At 50 eV and element
   scale the added conductance sigma_par/l ~ 1e7 S/m^2 also makes the u block of JOREK's coupled matrix far stiffer
   than anything in it now. That hits the preconditioner, which for n != 0 relies on per-harmonic blocks.
3. **A nonlinear inner iteration** (Newton or Picard) within a step, at least on the wall rows. This recovers
   B2.5's "iterate to consistency" but not its Ohmic anchor. Each iteration costs a matrix assembly and, for the
   axisymmetric problem, a new factorisation, which is JOREK's dominant cost. It also cures only the one-step
   overshoot, which F2 already fixes at almost no cost.

The honest summary is that JOREK's physics is right. A current in a flux tube does have inductance, and a sheath at
saturation is an open circuit. What SOLPS has and JOREK lacks is the separation of time scales that makes the
electrostatic model valid for transport. Rather than imposing that separation (option 2), the practical route is to
make sure the inductive current reaching a sheath is always a current that the sheath and the tube can carry (F0,
F1, F3, F4), and to linearise the sheath the way SOLPS does (F2). The Robin form (F6), the total ion flux (F5) and the
Pedersen current (F7) are the long-term structure.

---

## 7. One-page summary for a supervisor meeting

**The problem.** We want the sheath current-voltage condition on the divertor targets because the HFSHD mechanism
(Senichenkov 2022) is a thermoelectric current from the hot outer to the cold inner target. A floating-potential BC
forbids that current. In JOREK the sheath row runs when recycling is low, but with correct recycling every run dies
as the inner target cools. At the inner strike point the demanded current climbs past the ion saturation current
(j/j_sat 0.4 to 390 in about 8 steps), and the wall potential runs to kilovolts in 2-3 steps. The same plasma with a
floating wall potential detaches cleanly to 0.4 eV.

**Why it works in SOLPS and SOLEDGE.** Their parallel current is Ohm's law, an algebraic function of the potential,
solved implicitly in one elliptic potential equation. A saturated sheath is anchored through the flux tube's
conductance to the other target, SOLPS linearises the sheath with a slope that is never zero, it iterates each step,
and every boundary other than a sheath carries zero current. A sheath is never asked for more current than the
physics can deliver.

**Why JOREK is different.** Three things are specific to JOREK:
1. **The parallel current has inductance.** It is a state variable (Ampere plus induction), not a function of the
   potential. Over one time step shorter than mu0 h^2/eta (about 4 JOREK time units at 50 eV and 2 mm), a flux tube
   behaves as a current source that keeps delivering its old current whatever the wall potential does. When the
   sheath saturates, the wall is an open circuit on an inductor, and the potential jumps by L dI/dt, which is kV.
2. **One linear solve per step, and the sheath's slope is linearised to zero at saturation and at the electron cap.**
   Nothing anchors a saturated node except a tiny polarisation capacitance that is proportional to density and
   vanishes as the cell empties.
3. **Currents the sheath cannot see.** The boundary treatment adds currents the sheath BC does not know about: the
   equilibrium current frozen by a Dirichlet row on the short floating wall segment at the inner strike point
   (2-20 j_sat, constant while j_sat falls), the equilibrium's open-field-line current at t = 0, and the implicit
   boundary fluxes of the vorticity equation (magnetisation current, E x B and viscous vorticity fluxes). No sheath
   can pass these currents, and no choice of slope, cap or ramp inside the row can fix a demand that has no solution.
   That is why every such variant failed or merely moved the crash.

**What to do, in order.**
1. **Zero current on the non-sheath wall segments** (SOLEDGE's rule). Implemented (`bcs%zj_zero`) and run: it was
   not the killer, the crash repeated unchanged (step 672). Kept for consistency.
2. **SOLPS's sheath linearisation** (Patankar slope, residual exact). Tried and removed: with the slope floored the
   row is effectively unenforced on the ion-saturated part of the target, and the run ended earlier (608 vs 672).
3. **Start consistently**: floating wall potential and no open-field-line current at t = 0, instead of ramping the
   characteristic.
4. **Measure, then remove, the implicit vorticity-row wall currents.**
5. Longer term: the Robin form of the sheath row, j_sat from the same total ion flux as the particle balance, and
   the ion-neutral Pedersen current for kinetic neutrals.

**The deeper option.** Giving JOREK's wall an electrostatic, SOLPS-like parallel conductance is possible but
expensive and breaks Ampere's law locally. Do it only if items 1-4 do not make the potentials settle at the tens of
volts the thermoelectric physics predicts.
