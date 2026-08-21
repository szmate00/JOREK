# Sheath Boundary Condition for the Electric Potential in JOREK model600

**Branch:** `bc-tests`  ·  **Model:** 600 (reduced MHD + kinetic neutrals)  ·  **Case:** AUG #38773, lower single null, `grid_to_wall = .t.`
**Status:** implemented and running; not yet stable through the full timestep ramp.
**Last updated:** 2026-08-18

---

## 0. TL;DR for the impatient

We want the SOL electric potential to be *set by the sheath physics at the divertor targets* instead of being pinned to zero by a Dirichlet condition. The physics is a textbook Langmuir characteristic, and the implementation is faithful to it. The trouble is not the physics — it is that **in JOREK the magnetic field is a dynamic variable, and a potential that varies along the wall drives an E×B flow through the wall that drags poloidal flux into a resistive layer far thinner than the mesh.** That produces a grid-scale current filament at the divertor which grows and eventually kills the run.

We have identified that mechanism, confirmed it experimentally, and built a partial cure (`sheath_wall_vel`, a thin-resistive-wall relaxation) which reduces the filament amplitude by roughly three orders of magnitude. It is not yet a complete cure. The remaining growth is rate-limited by the boundary condition's per-step relaxation gain, which currently scales with `tstep` and therefore destabilises at every timestep increase.

**Next three runs:** a no-sheath control, a pinned per-step gain, and a `sheath_wall_vel` scan. Details in [§9](#9-next-steps-and-why).

---

## 1. What we are trying to achieve

In the scrape-off layer the plasma terminates on material surfaces. A Debye sheath forms there, and it sets the electrostatic potential of the plasma relative to the wall. That potential is not a free parameter — it is determined by the requirement that ion and electron fluxes to the surface balance (or, if a current flows, by how far they are out of balance).

This matters for us because:

- The **radial electric field in the SOL** and hence the E×B flow pattern is largely determined by the sheath potential and its variation along the target.
- **Parallel currents** (thermoelectric currents, in particular) flow between the two targets whenever the target electron temperatures differ. These currents are a genuine part of divertor physics and are entirely absent if the potential is pinned.
- Any study of **divertor detachment asymmetries, HFS high-density regions, or target current patterns** needs the potential to be free.

With the standard JOREK boundary condition, `u = 0` at the wall, none of this exists. The goal of this work is to replace that with the sheath's own current–voltage relation, so that the potential is *evolved self-consistently from the boundary*.

> [!note] What "u" is
> In model600, `u` is the electric potential stream function. The electrostatic potential is $\Phi = -F_0\,u$, and the E×B velocity is $\mathbf{v}_E = R\,\nabla u\times \mathbf{e}_\varphi$. The sign matters a great deal here — see [§3.3](#33-sign-conventions-the-thing-that-bites).

### 1.1 The actual goal: HFSHD

Everything in this document exists to answer one physics question: **can JOREK reproduce the high-field-side high-density (HFSHD) region observed in AUG?**

The accepted picture from SOLPS/EDGE2D is that the in–out density asymmetry is substantially **E×B drift driven**:

1. The sheath sets $\Phi \approx \Lambda T_e/e$ at each target.
2. $T_e$ varies along and across the target, so $\Phi$ does too — giving a radial electric field in the SOL, $E_r \approx \Lambda\,\nabla T_e/e$.
3. $E_r \times B$ drives a poloidal drift, notably **through the private flux region**, which carries plasma from the outer divertor toward the inner one.
4. Density piles up on the high-field side.

With the standard `u = 0` Dirichlet condition, step 1 does not happen, so steps 2–4 cannot. **JOREK is structurally unable to produce drift-driven HFSHD without this boundary condition.** That is the entire motivation.

An earlier apparent $\eta$-scaling of HFSHD turned out to be **100% an artefact** of `keep_current_prof` forcing the current profile in the SOL and PFR (bracketed and closed 2026-07-23). The remaining candidate knobs identified at that time were: the leg $\eta$ floor, the **sheath potential boundary condition**, and the viscosities. This work is the second of those.

> [!success] Good news on the leading order
> Even with `dirichlet%zj = .true.` — i.e. with the boundary current frozen — the **dominant driver survives**. The potential is $\Phi = (\Lambda + X)\,T_e/e$, and the $\Lambda T_e/e$ part responds to the evolving $T_e$. So $E_r \approx \Lambda\nabla T_e/e$ is captured. What is missing is the *thermoelectric* contribution, which needs the boundary current to be live — see [§8](#8-why-other-codes-do-not-have-this-problem) on free boundary.

---

## 2. The physics

### 2.1 The sheath current–voltage characteristic

The classical result (Stangeby, *The Plasma Boundary of Magnetic Fusion Devices*, eq. 2.68; equivalently the Langmuir probe characteristic) is:

$$
j = j_{\rm sat}\left(1 - e^{-X}\right),
\qquad
X = \frac{e\Phi}{k T_e} - \Lambda,
\qquad
\Phi = V_{\rm sheath\ entrance} - V_{\rm wall}
$$

Reading it physically:

| Regime | Condition | Meaning |
|---|---|---|
| **Floating** | $X = 0$, i.e. $\Phi = \Lambda T_e/e$ | Zero net current. The surface charges up until ion and electron fluxes balance. |
| **Ion saturation** | $X \to +\infty$ | Potential is very positive; electrons are fully repelled; $j \to j_{\rm sat}$ and stops responding. |
| **Electron saturation** | $\Phi \to 0$ ($X \to -\Lambda$) | Potential at wall level; electrons stream in freely; $|j|$ becomes large. |

The sheath factor is

$$
\Lambda = \Lambda_0 - \ln\sqrt{\gamma\left(1 + T_i/T_e\right)},
\qquad
\Lambda_0 = \ln\sqrt{\frac{m_i}{2\pi m_e}} \approx 3 \ \ \text{(deuterium)}
$$

so $\Lambda \approx 2.4$–$2.6$ at $T_i \approx T_e$. The floating potential at a 20 eV target is therefore about **50 V**, which is not a small number — it is comparable to the entire potential variation across the SOL.

The ion saturation current uses the same parallel outflow that JOREK's Mach-1 boundary condition already imposes:

$$
j_{\rm sat} = c_{\rm sat}\,\rho\;\frac{g(b_n)\,c_s}{|B|},
\qquad
c_s = \sqrt{\gamma\,(T_i+T_e)}
$$

where $g(b_n)$ is the Chodura–Riemann smoothing function — the same `direction * factor` that the Mach-1 condition uses, including its sign. That reuse is deliberate: if the sheath current and the Mach-1 outflow disagreed about how much plasma reaches the wall, the two boundary conditions would fight each other.

### 2.2 Forward versus inverted

There are two ways to impose this relation, and they are *not* numerically equivalent.

**Forward** — impose the current given the potential, $j = j(\Phi)$. This is single-valued, and its derivative

$$\frac{\partial j}{\partial \Phi} \propto e^{-X} \longrightarrow 0 \quad\text{as}\quad j\to j_{\rm sat}$$

goes smoothly to zero at ion saturation. The linearisation never becomes singular. No clipping is needed on the ion side.

**Inverted** — impose the potential given the current, $\Phi = \Phi(j)$:

$$\frac{\partial \Phi}{\partial j} \propto \frac{1}{1 - j/j_{\rm sat}}$$

which **diverges exactly at ion saturation** — which is precisely where a divertor target sits. This form must be clipped, and clipping has consequences we will come back to in [§7.3](#73-mechanism-c-the-clip-makes-the-row-explicit).

> [!important] Both forms are implemented
> The forward form is the natural-BC path (`bcs%natural%u`). The inverted form is the nodal path (`bcs%sheath_u`). We built the forward one *because* it is better conditioned — but as it turns out, the nodal one runs much further. See [§7.1](#71-mechanism-a-the-natural-w-and-zj-terms-were-mis-linearised).

### 2.3 What we should see if it works

- $e\Phi/kT_e \approx \Lambda \approx 2.4$ wherever the target current is small.
- The electric field along the target $\approx \Lambda \nabla T_e / e$ — i.e. the potential tracks the target temperature profile.
- Net current between inner and outer target when their temperatures differ (thermoelectric current).
- Total current to the wall integrating to approximately zero (charge conservation), which the diagnostic reports as `I_wall`.

---

## 3. How a boundary condition like this goes into JOREK

This section is for the reader who knows JOREK but has never opened the boundary-condition machinery. It is the part that makes this problem harder than it looks.

### 3.1 The finite elements

JOREK uses **bicubic Bézier elements with C¹ continuity**. Each node carries `n_order = 3`... more precisely, each node has **four degrees of freedom** per variable:

| DOF index | Meaning |
|---|---|
| `1` | the value |
| `iv_dir` (`2` or `3`) | first derivative along one logical direction |
| `iv_perp_dir` | first derivative along the other |
| `4` | the mixed second derivative $\partial^2/\partial s\,\partial t$ |

At a boundary node, one of the two first-derivative directions runs **along** the boundary (tangential) and the other runs **into** the domain (normal). This distinction is the source of most of the difficulty.

### 3.2 The two ways to impose a boundary condition

**(a) Penalty Dirichlet — `models/model600/mod_boundary_conditions.f90`**

JOREK imposes Dirichlet conditions by *overwriting* matrix rows with a large diagonal entry, `zbig = 1e12`. The routines are

```fortran
call boundary_conditions_add_one_entry( row_node, row_var, in, col_node, col_var, in, value, ... )
call boundary_conditions_add_RHS      ( row_node, row_var, in, ..., value, ... )
```

> [!warning] These **overwrite**, they do not accumulate
> If two pieces of code write the same row, the last one wins silently. This is why the wall-relaxation block has to explicitly re-add `loop_voltage` — otherwise it would destroy it.

The great advantage of this route: you can write **any linear relation** into the row, coupling any variables at that node. That is how the nodal sheath condition works — the row reads

$$
\delta u + \sum_k c_k\,\delta x_k = u_{\rm target} - u^{\,0},
\qquad c_k = -\frac{\partial u_{\rm target}}{\partial x_k}
$$

with $x_k \in \{\rho, zj, T_i, T_e\}$. It is a **collocation** condition: exact at the nodes, nothing said in between.

**(b) Natural / surface term — `models/model600/mod_boundary_matrix_open.f90`**

The weak form of every equation generates a surface integral when you integrate by parts. JOREK normally throws these away (equivalent to assuming the corresponding flux is zero) and imposes Dirichlet instead. A "natural" boundary condition puts the surface term *back*, with the physical flux you want:

```fortran
rhs_ij(var_u) = - v * BigR * ( zj_sh - zj0 ) * sh_Bn * dl * tstep * sheath_ramp
```

This is the charge-continuity equation with the sheath current as the boundary flux. It is a **weak** condition — satisfied in an integrated sense over the boundary, not pointwise — which is generally much better behaved numerically.

> [!danger] The restriction that shapes everything
> `mod_boundary_matrix_open` writes rows and columns **only at DOFs `direction(j)`**, that is: the value and *one* tangential derivative. It never touches the normal-derivative DOF. A surface integral can therefore only reach test functions that are non-zero on the boundary edge — which are exactly the rows a Dirichlet on the same variable overwrites.
>
> Consequence: **a natural BC on a variable that also has a Dirichlet does nothing at all.** And a natural BC whose residual depends on a *normal* derivative cannot be correctly linearised, because the required Jacobian columns are not in the available set. The code synthesises a fictitious normal derivative via `psi_t = H1(k,l,ms) * element_size_kl * element_size_perp`, which is not the real thing.

### 3.3 Sign conventions — the thing that bites

The JOREK reference paper (Hoelzl et al 2021, eq. 26) defines

$$u \equiv \Phi/F_0, \qquad \mathbf{v}_{\rm pol} = -R\,\nabla u\times\mathbf{e}_\varphi$$

but the **model600 element matrix implements** $\mathbf{v}_{\rm pol} = +R\,\nabla u\times\mathbf{e}_\varphi$ (check the E×B advection and compression terms of the density equation, and the $[u,\psi]$ term of the induction equation). Therefore

$$\boxed{u_{\rm code} = -\,u_{\rm paper}, \qquad \Phi = -F_0\,u_{\rm code}}$$

and likewise $zj_{\rm code} = \Delta^*\psi = -j_{\rm paper}$.

This propagates into a single leading minus sign in the normalisation constant:

```fortran
a_n   = - 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * rho0) / m_i   ! the MINUS is the convention above
c_sat = - 0.5d0 * a_n
```

Every other coefficient derives from `a_n`, so getting this one sign right fixes the whole chain. We got it wrong initially, and the symptom was a boundary condition that pushed the potential the wrong way — an anti-damped feedback loop.

### 3.4 Which boundary types

With `grid_to_wall = .t.`, the grid boundary follows the physical first wall. The boundary types in this case are `0, 1, 2, 3, 4, 5, 9`, with `0` being the interior (not a boundary). The sheath is applied on `1, 4, 5, 9` — the types that carry the strike points.

> [!caution] A single boundary type covers very different physics
> Type 1 covers both the divertor targets **and** stretches of the main chamber wall where the field is nearly tangential. At a tangential wall there is no parallel flux to the surface, hence no sheath at all — and worse, moving along such a wall keeps you on roughly the same flux surface, which is the geometry in which an imposed potential drags flux hardest. This is what `sheath_min_bn` is for ([§6.4](#64-grazing-incidence-and-the-obliqueness-gate)).

---

## 4. What was implemented

All of the following is committed on `bc-tests`, current head `2d649c4e2`.

### 4.1 New files

| File | Purpose |
|---|---|
| `models/model600/mod_sheath_bc.f90` | **The single home of the sheath physics.** Stateless; derives its constants from `phys_module` on the fly. Public: `sheath_norm`, `sheath_get_lambda`, `sheath_x_limited`, `sheath_current`. Used by all three consumers so they can never drift apart. |
| `models/model600/mod_sheath_diag.f90` | Wall-current diagnostic (`sheath_diag_reset/add/report`, array-only MPI collectives). Also holds `sheath_init_potential` (floating-potential initialisation) and `sheath_store_psi0` (records ψ's DOFs at first call, for the resistive wall). |
| `util/sheath_bc_unit_test/` | Standalone test harness — see [§8](#8-testing-without-a-cluster). |

### 4.2 Modified files

| File | What changed |
|---|---|
| `models/model600/mod_boundary_matrix_open.f90` | Three surface terms (`var_u`, `var_w`, `var_zj`), gated by `apply_natural_bc(...)`, with their full Jacobians. Stiffness cap, ramp, diagnostic hook. |
| `models/model600/mod_boundary_conditions.f90` | The nodal `sheath_u` path: characteristic inversion, obliqueness gate, relaxation, `sheath_u_align_psi`, `sheath_u_value_only`, and the thin-resistive-wall block. |
| `models/phys_module.f90` | `type_natural_bc` gained `u`, `w`, `zj`; `type_bcs` gained `sheath_u`; ~14 new scalar parameters. |
| `models/preset_parameters.f90` | Defaults — **all new switches off**, so the baseline is bit-for-bit unchanged. |
| `communication/broadcast_phys.f90` | MPI pack/unpack (verified: +14 reals/+14 logicals, divergence index shifted by exactly the number added). |
| `models/mod_log_params.f90`, `models/model600/initialise_parameters.f90` | Logging, namelist entries, validation. |
| `matrix/construct_matrix_mod.f90` | `#if JOREK_MODEL == 600` guarded diagnostic reset/report. |
| `jorek2_main.f90` | `#if JOREK_MODEL == 600` guarded `sheath_init_potential` call. |

### 4.3 The characteristic, in code

From `mod_sheath_bc.f90`:

```fortran
a_n    = - 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * rho0) / m_i
c_sat  = - 0.5d0 * a_n
vw     =   EL_CHG * sheath_V_wall * MU_ZERO * central_density * 1.d20
lam    =   lam0 - 0.5d0 * log( GAMMA * T_l / Te_l )        ! sheath_Lambda_local
x      = ( 0.5d0 * a_n * u - vw ) / Te_l - lam
x_lim  =   sheath_X_min + dx * log( 1.d0 + exp(z) )        ! smooth electron-side limiter
zj_sat =   c_sat * rho_l * g_eff * cs / Btot
zj_sh  =   zj_sat * ( 1 - exp(-x_lim) )
```

The limiter `sheath_x_limited` is a softplus: C¹, monotone, and finite over $X \in [-200, 200]$ (verified in the unit test). It replaces a hard clip on the electron-saturation side, where the exponential would otherwise overflow.

---

## 5. The two implementation routes, side by side

| | **Natural BC** (`bcs%natural%u`) | **Nodal** (`bcs%sheath_u`) |
|---|---|---|
| File | `mod_boundary_matrix_open.f90` | `mod_boundary_conditions.f90` |
| Form | forward, $j = j(\Phi)$ | inverted, $\Phi = \Phi(j)$ |
| Enforcement | weak (integrated over the edge) | collocation (exact at nodes) |
| Conditioning | benign; $\partial j/\partial\Phi \to 0$ at saturation | singular at saturation; **must clip** |
| `dirichlet%u` | must be `.false.` | must be `.true.` (the sheath overwrites the row) |
| Mechanism C (per-step gain) | **absent** — no clip, full Jacobian always | present — clip zeroes `sh_dr` |
| Best run | ~19 steps *(contaminated — see §9.6)* | **202 steps** |

The forward form is theoretically the better object, and it is the one the module header documents as the recommended path. Empirically the nodal form runs an order of magnitude further. The reason is not that the forward form is wrong — it is that its *companions* were wrong. See next section.

---

## 6. The knobs, and what each is for

All defaults are off / neutral, so nothing here changes a baseline run.

### 6.1 Core

| Parameter | Default | Meaning |
|---|---|---|
| `sheath_Lambda` | `3.d0` | $\Lambda_0$. Set $\le 0$ to compute it from `central_mass`. |
| `sheath_Lambda_local` | `.true.` | Use $\Lambda = \Lambda_0 - \ln\sqrt{\gamma(1+T_i/T_e)}$, consistent with $c_s$. |
| `sheath_V_wall` | `0.d0` | Wall potential in volts. 0 = grounded. |

### 6.2 Clipping (nodal path only)

| Parameter | Default | Meaning |
|---|---|---|
| `sheath_u_exp_max` | `2.d0` | Upper clip on $X$. Caps $\Phi$ at $(\Lambda + X_{\max})T_e/e$ and $\lvert j\rvert$ at $(1-e^{-X_{\max}})j_{\rm sat}$. |
| `sheath_u_exp_min` | `-3.d0` | Lower clip. Effective electron-saturation limit. |

With `exp_max = 2.0` the clip sits at $r = 1 - e^{-2} = 0.865$. Our diagnostic reports $\lvert j/j_{\rm sat}\rvert \approx 6.5$ at the targets, so **most target nodes are clipped**. This turns out to matter enormously — [§7.3](#73-mechanism-c-the-clip-makes-the-row-explicit).

### 6.3 Relaxation

| Parameter | Default | Meaning |
|---|---|---|
| `sheath_u_relax` | `1.d0` | Per-step fraction: `u` moves this far toward the characteristic each step. 1 = no relaxation. |
| `sheath_u_relax_time` | `-1.d0` | If > 0, replaces the above by `min(1, tstep/sheath_u_relax_time)`. |

`u` at the wall is slaved *algebraically* to $T_e$, with no inertia and no dissipation anywhere in the loop

$$T_e \longrightarrow u \longrightarrow \mathbf{v}_E \longrightarrow T_e$$

so the loop needs a response time of its own. Expressing it as a physical time seemed the principled choice, since the meaning of a per-step fraction changes at every `tstep_n` block boundary. **That reasoning turns out to be backwards for this failure mode** — see [§7.3](#73-mechanism-c-the-clip-makes-the-row-explicit).

### 6.4 Grazing incidence and the obliqueness gate

`sheath_min_bn` (default `0.d0`) does two related jobs.

Through a point where the field goes tangent to the wall, $g(b_n)$ passes through zero and **changes sign**, so the ratio $r = zj/j_{\rm sat}$ runs to $+\infty$ on one side and $-\infty$ on the other. Both ends clip, and the imposed potential jumps by $(X_{\max}-X_{\min})T_e$ between two neighbouring nodes — a positive blob right next to a negative one, anchored at the tangency point. Physically nothing is wrong there: no parallel flux reaches the wall, the sheath carries no current, and the surface simply floats.

The code writes the derivative as

```fortran
sh_dr_dzj = sh_g / ( sh_C * ( sh_g**2 + sheath_min_bn**2 ) )
```

which both weights the ratio by $g^2/(g^2+g_{\min}^2)$ *and* removes the division by zero. The same weight is folded into the relaxation

```fortran
sh_wgt_bn = sh_g**2 / ( sh_g**2 + sheath_min_bn**2 )
sh_relax  = sh_relax * sh_wgt_bn
```

so the row degenerates **smoothly to $\delta u = 0$** — i.e. exactly the standard Dirichlet — wherever the field is tangential. This is what protects the main-chamber portions of boundary type 1.

### 6.5 The thin resistive wall

| Parameter | Default | Meaning |
|---|---|---|
| `sheath_wall_vel` | `0.d0` | Speed at which the wall lets poloidal flux through: $\partial\psi/\partial t = -v_w\left(\partial\psi/\partial n - \left.\partial\psi/\partial n\right\rvert_{t_0}\right)$ |

This is the cure for the dominant failure mode. It is discussed in full in [§7.2](#72-mechanism-b-flux-dragging--the-main-event).

> [!warning] It is mesh-dependent
> `sheath_wall_vel` is expressed per degree-of-freedom unit of $\partial\psi/\partial n$, so its useful value depends on the mesh. **Scan it**; do not port a number between cases.

### 6.6 Experimental switches (all tried, all currently off)

| Parameter | What it does | Verdict |
|---|---|---|
| `sheath_u_align_psi` | Forces $[u,\psi]=0$ at the wall, i.e. the E×B flow runs along flux surfaces. Replaces the vorticity row at `u`'s normal-derivative DOF. | Divertor went clean; trouble **moved** to the main chamber. Worse overall (79/104 steps). |
| `sheath_u_value_only` | Imposes the characteristic on the node value only, leaving `du/dl` to the vorticity equation. | Much worse (34 steps): `w` runs free. |
| `sheath_stiff_max` | Caps the sheath Robin term relative to the row's own polarisation diagonal. | No measurable effect. |
| `sheath_ramp_time` | Linearly ramps the surface term in from `t_start`. | Was set far too large relative to `tstep` in early tests — the term was effectively never on. |
| `sheath_flux_sign` | Debug ±1 multiplier on the surface term. | Setting it to 0 **still crashed**, proving the sheath term itself was not the unstable object in the natural-BC runs. |
| `sheath_wall_diff` | Earlier attempt at wall flux relief, driven by `zj`. | **Inert.** `dirichlet%zj = .true.` freezes `zj`, so the term had nothing to act on: 3e-6 and 1e-5 gave *identical* crashes at step 208. Removed from consideration. |

---

## 7. What we ran, and what it taught us

### 7.0 The run ladder

All runs on the same AUG case with

```
tstep_n = 1.d-3, 1.d-2, 1.d-1, 3.d-1, 1.d+0, 2.d+0, 10.
nstep_n =   100,   100,   100,   100,   100,   100,   1000000
```

so the timestep blocks are: **steps 1–100 @ 1e-3, 101–200 @ 1e-2, 201–300 @ 1e-1, 301–400 @ 3e-1.**

| # | Configuration | Died at | Block | Character |
|---|---|---|---|---|
| V2 | `natural%zj` | 19 | 1e-3 | Boundary current structures, then inward |
| V3 | `natural%w` | 24 | 1e-3 | Same, localised at HFS divertor |
| V4 | `natural%u` + all variants (sign flip, stiffness cap, floating start, `flux_sign=0`) | ~19 | 1e-3 | Same |
| — | nodal `sheath_u` | 174 | 1e-2 | ± potential blob at HFS divertor, grows fast |
| — | + `sheath_min_bn = 0.2` | 199 | 1e-2 | Same, slower |
| — | + `exp_max = 0` (floating limit) | 173 | 1e-2 | Same |
| — | + `sheath_u_relax_time = 1.0` | 210 | 1e-1 | Same |
| — | **`eta = 2e-6` (×100)** | **305** | **3e-1** | Clean to 301, then fast growth |
| — | `sheath_wall_diff` = 3e-6 / 1e-5 | 208 / 208 | 1e-1 | Identical → term inert |
| — | `sheath_u_align_psi` (no gate) | 79 | 1e-2 | `w`, `zj` blew up on **main chamber** |
| — | `sheath_u_align_psi` + obliqueness gate | 104 | 1e-2 | Same |
| — | `sheath_u_value_only` | 34 | 1e-3 | `w` runs free |
| — | **`sheath_wall_vel = 4.d-3`** | **202** | **1e-1** | **Fields ~3000× smaller; filament still present** |

### 7.1 Mechanism A: the `natural%w` and `natural%zj` terms were mis-linearised

The first three variants all died at around 20 steps with the same signature. The cause is the restriction described in [§3.2](#32-the-two-ways-to-impose-a-boundary-condition):

The residuals of the `w` and `zj` surface terms depend on a **normal derivative** ($\nabla u\cdot\mathbf{n}$ and $\nabla\psi\cdot\mathbf{n}$ respectively), but the trial-function loop in `mod_boundary_matrix_open` runs over `l2 = direction(l)` only — the value and tangential DOFs. The normal derivative is *synthesised* from the value/tangential basis functions, which is not the true Hermite normal-derivative basis function.

The result is **a missing Jacobian entry plus a spurious one** — i.e. an effectively explicit boundary term of order $1/h$. Explicit $O(1/h)$ terms grow boundary-localised structures over ~20 steps and then propagate them inward. That is exactly what was observed.

> [!success] Action taken
> `bcs%natural%w` and `bcs%natural%zj` still exist in the code but are now **refused by `initialise_parameters`** with an explanatory message. Lifting the restriction properly needs an extra trial index carrying `direction_perp(1)` with the correct Hermite basis (zero on the edge, unit normal derivative) — a real but bounded piece of work in the element machinery.

A corollary worth recording, because we initially got it wrong: **there is no hidden $\nabla u\cdot\mathbf{n}=0$ condition** being imposed by dropping the `w` surface term, and no spurious $\nabla\psi\cdot\mathbf{n}/h$ current to repair. With `dirichlet%w` and `dirichlet%zj` on, the dropped surface terms impose *nothing at all*, because they can only reach rows that the Dirichlet overwrites anyway. So keeping those Dirichlets is correct and costs nothing.

### 7.2 Mechanism B: flux dragging — the main event

This is the important one, and it is specific to JOREK.

Once $u$ varies along the wall there is an E×B flow **through** the boundary:

$$v_E\cdot\mathbf{n} = R\,\frac{\partial u}{\partial \ell}$$

The parallel part of the flow is harmless, because $\mathbf{B}\cdot\nabla\psi = 0$. The *perpendicular* part is not: through the $R[u,\psi]$ term of the induction equation it **drags poloidal flux toward the boundary**. Against a Dirichlet $\psi$, that dragged flux has nowhere to go. It piles up until resistive diffusion can carry it away, which happens in a layer of width

$$\delta = \frac{\eta}{v_n}$$

At Spitzer resistivity ($\eta = 2\times10^{-8}$) this is of order **10 µm**, against a mesh of order **1 mm**. The layer is unresolvable by a factor of a hundred, and what the code produces instead is a grid-scale alternating current structure — the filament we have been fighting all along.

**The evidence is consistent and quantitative:**

1. **η is the only knob that ever bought a real factor.** Increasing $\eta$ by 100× thickened the layer by 100× and moved the crash from 199 to 305 — a whole `tstep` block further. Nothing else came close.
2. **The mode lives exactly where the drag is strongest.** It appears at the HFS divertor and along stretches of the main chamber where the wall runs nearly parallel to a flux surface — the geometry in which an imposed potential is maximally *not* a flux function, and therefore drags hardest.
3. **`sheath_u_align_psi` moved the problem rather than solving it.** Forcing $[u,\psi]=0$ cleaned the divertor completely and pushed the trouble onto the main chamber wall — where boundary type 1 has no business hosting a sheath at all.
4. **The cure works.** `sheath_wall_vel = 4.d-3` lets the wall pass flux at a finite rate,
   $$\frac{\partial\psi}{\partial t} = -v_w\left(\frac{\partial\psi}{\partial n} - \left.\frac{\partial\psi}{\partial n}\right|_{t_0}\right)$$
   and the amplitudes at cycle 200 dropped from `zj` max **4387 → 1.53** MA/m² and `w` max **2344 → 16.8**. That is a factor of ~3000 in `zj`, at a *later* cycle, at physical resistivity.

> [!important] Why the deviation, not the raw gradient
> The relaxation is driven by the **deviation of $\partial\psi/\partial n$ from its value at $t_{\rm start}$**. Using the raw gradient would bleed flux even in a perfectly quiet plasma and slowly destroy the equilibrium. `sheath_psi0` (in `mod_sheath_diag`) stores the reference state, and is filled **lazily inside `boundary_conditions`** from the very `node_list` that routine indexes — doing it in `jorek2_main` caused a SIGSEGV, because MPI node redistribution happens in between and the array was then indexed with mismatched node numbering.

**And yet: the filament is still there.** Zooming into the X-point region at cycle 200 with a ±1 MA/m² scale shows unmistakable mesh-scale alternating stripes at the HFS target ($R \approx 1.26$–$1.30$, $Z \approx -0.97$ to $-1.02$), a sharp line at the wall corner near $Z \approx -1.13$, and matching striping along the outer boundary at $R \approx 1.60$–$1.67$. Same location, same character, three orders of magnitude smaller. **`sheath_wall_vel` reduces the drive; it does not remove the mechanism.**

### 7.3 Mechanism C: the clip makes the row explicit

Why do so many crashes land within a few steps of a `tstep` increase?

The nodal row is

```
zbig*du + zbig*relax*sum_k coef_k*dx_k  =  zbig*relax*(u_target - u0)
```

The coefficients `coef_k` are the linearisation of $u_{\rm target}$ with respect to $\rho$, $zj$, $T_i$, $T_e$ — that is what makes the condition *implicit* in the plasma state. But look at what happens when the ratio is clipped:

```fortran
sh_dr = 1.d0
if ( sh_ratio /= sh_ratio_raw ) sh_dr = 0.d0
```

and every one of `coef_rho`, `coef_zj`, `coef_Ti` is proportional to `sh_dr`. **When clipped, they all vanish**, leaving only the direct $T_e$ term:

$$\text{clipped:}\qquad \delta u + \text{relax}\cdot c_{T_e}\,\delta T_e = \text{relax}\,\left(u_{\rm target} - u^0\right)$$

So $u$ is slaved to $T_e$ alone, with $\rho$ and $zj$ **lagged at the old timestep**. That is an explicit coupling, and explicit couplings have a per-step gain limit. The gain is exactly

$$\text{relax} = \frac{\Delta t}{\tau}$$

| Block | `tstep` | `relax` (τ = 1) | Outcome |
|---|---|---|---|
| 101–200 | 1e-2 | **0.01** | survived all 100 steps |
| 201–300 | 1e-1 | **0.1** | died at step 202 |

Since `sheath_u_exp_max = 2.0` clips at $r = 0.865$ while the targets sit at $\lvert j/j_{\rm sat}\rvert \approx 6.5$, **most target nodes are in the clipped, explicit branch**.

> [!important] The unified picture
> The stripe is a numerical instability of the lagged boundary coupling. **`sheath_wall_vel` reduces its drive; `relax` sets its growth rate.** At `relax = 0.01` it grows slowly — 100 steps to reach 1.5 MA/m². At `relax = 0.1` it crosses the stability threshold and blows up in two. Runs that die *mid-block* (34, 79, 104, 173, 174, 199) were killed by the filament reaching amplitude; runs that die *just after a jump* (202, 208, 210, 305) were killed by the growth rate suddenly increasing. Same mode, two ways to reach the end.

### 7.4 A note on what is *not* the problem

Several things were ruled out along the way, and are worth recording so nobody re-derives them:

- **`bc_natural_open` was already in the baseline** — it is not a confound.
- **The sheath term itself is not the unstable object** in the natural-BC runs: `sheath_flux_sign = 0` turns the term completely off and the run *still* crashed at ~19 steps.
- **The standard Dirichlet DOF set (value + tangential) is optimal.** Deviating in *either* direction is worse: freeing the tangential row (`sheath_u_value_only`) gives 34 steps, additionally pinning the normal derivative (`sheath_u_align_psi`) gives 79–104.
- **`sheath_wall_diff` was inert**, as established above.
- **The `find_flux_surfaces` "another solution" warnings** in the crash logs are unrelated to this work; buffered Fortran output merely truncated the log at that point.
- The retracted $(B/B_\varphi)^2$ "fix" from the first review: Artola's $J_\parallel = -jB/F_0$ assumes a force-free SOL current, which is the physically appropriate interpretation. Only the sign was wrong.

---

## 8. Why other codes do not have this problem

This question came up repeatedly and the answer is structural, not a matter of cleverness.

**SOLPS-ITER, SOLEDGE2D/3D, EDGE2D and GBS all hold the magnetic field fixed.** They solve transport on a prescribed equilibrium. There is no induction equation, no $[u,\psi]$ term, and therefore **no flux to drag**. The E×B flow at the wall is a transport flow and nothing more. Their entire difficulty with sheath BCs is elliptic-solver conditioning, which is a well-understood and much milder problem.

What they *do* struggle with is instructive anyway, and we have borrowed from it:

| Code | Practice | Our analogue |
|---|---|---|
| SOLPS | `BCPOT` options: 1 = value, 2 = gradient, 4 = weak/mixed, 5 = current density, 11 = sheath, 15 = feedback on total current | We have the value form (nodal) and the current-density/weak form (natural) |
| SOLPS | Artificial cross-field conductivity $\sigma_{\rm AN}$ to regularise the potential equation | No direct analogue yet; `eta_num` is the closest |
| SOLPS | Drift ramping: 10% → 50% → 80% → 100% over many steps | `sheath_ramp_time` |
| SOLPS | BC under-relaxation from 0.01 down to 0.001 | `sheath_u_relax` — **and note how small their numbers are** |
| SOLPS | Cells below ~1 mm produce potential oscillations | We are at ~1 mm and see exactly this |
| GBS | $v_{\parallel e} = \pm c_s\exp(\Lambda_0 - e\phi/T_e)$; Dirichlet $\phi = \Lambda T_e/e$ away from the strike points | Our floating-potential initialisation `sheath_init_u` |
| GBS | *"The Poisson equation is ill-defined if a Neumann BC is applied at the four walls"* | Why at least one boundary type must keep `dirichlet%u = .true.` |

The one published JOREK result with a Langmuir sheath BC is **Artola et al., COMPASS VDE modelling (arXiv:2101.01755)**. It is worth reading carefully because of what it says and what it uses:

- It describes the Langmuir BC as *"numerically very challenging"* — in a JOREK paper, from the person who derived the formulation.
- It uses a **one-sided limiter geometry**, not a divertor.
- It uses penalty $Z = 10^{12}$, $\gamma_{\rm sh} = 11$.
- **It uses free boundary with STARWALL.**

That last point is not incidental. A resistive wall lets $\partial\psi/\partial t \neq 0$ at the boundary, so the dragged flux is *absorbed* rather than piling into an unresolvable layer. Our `sheath_wall_vel` is a local, cheap emulation of exactly that.

> [!question] Does free boundary replace the sheath BC?
> **No.** The sheath BC is orthogonal to it and everything we have built stays exactly as it is. What changes is only the ψ treatment: `is_freebound(in,k)` starts returning true for the vacuum-coupled harmonics, so `boundary_conditions` skips the ψ *and* `zj` Dirichlets on its own and STARWALL's response matrix takes over. The wall-relaxation block already guards on `is_freebound`, so the two cannot fight.
>
> There is a bonus: with `zj` no longer frozen at the boundary, the current entering the characteristic becomes **live**, and the j–V loop finally closes — the self-consistency that was the original point and that `dirichlet%zj = .true.` has been blocking all along. The cost is that the inverted characteristic's saturation behaviour becomes live too, so `sheath_u_exp_max` and `sheath_min_bn` become physics choices rather than safety rails.

---

## 9. Next steps, and why

### 9.0 The timescale constraint — this sets the priority

Before choosing what to run, work out how long the run has to be.

From the code's own output at cycle 200: $t = 7.17715\times10^{-7}$ s after $100\times10^{-3} + 100\times10^{-2} = 1.1$ JOREK time units, so

$$\tau_A = 6.53\times 10^{-7}\ \text{s}$$

Now the timescale HFSHD actually needs. The E×B drift speed is $v_E \sim E/B \sim (50\,\text{V}/0.05\,\text{m})/2.5\,\text{T} \approx 400$ m/s; over a poloidal distance of order 0.5 m that is $\sim1.2$ ms. Parallel equilibration is $L_\parallel/c_s \sim 20\,\text{m}/3\times10^4\,\text{m s}^{-1} \approx 0.7$ ms. Either way:

$$t_{\rm HFSHD} \sim 1\ \text{ms} \approx 1500\,\tau_A$$

Against that, here is what each timestep block can deliver:

| `tstep` | steps needed for 1500 $\tau_A$ | verdict |
|---|---|---|
| 1e-2 | 150 000 | ✗ not feasible |
| 1e-1 | 15 000 | ✗ marginal at best |
| **1e0** | **1 500** | ✓ usable |
| **1e1** | **150** | ✓ comfortable |

> [!danger] The ramp is not optional
> The 202-step run reached $t = 1.1\,\tau_A$ — **three orders of magnitude short** of what HFSHD needs. Getting the boundary condition to run stably at `tstep = 1`–`10` is therefore not a convenience, it is the entire critical path. This is also why the base case's schedule ends at `tstep = 10` with `nstep = 1000000`.

This reprioritises the runs below. **Run B** is the one that matters among the nodal variants — and **[Run D](#96-run-d--retry-the-forward--natural-route-promoted-this-may-be-the-answer)** may matter more than any of them, because the forward/natural route has no per-step gain limit at all. Read §9.6 before committing machine time. Restated: **Run B is the one that matters**, because it is the only one that directly attacks "works at every timestep". Runs A and C are cheap and informative, but B is the one that decides whether this project can reach its physics goal.

### 9.1 Run A — the no-sheath control *(highest information, still not done)*

```fortran
 bcs(1)%sheath_u = .false.     ! and 4, 5, 9
```

Everything else identical, including `sheath_wall_vel`.

**Why.** Every crash in this campaign has been attributed to the boundary condition, but we have never checked whether *this case with this timestep ramp* survives without it. Note that the wall-relaxation block sits inside `if (apply_sheath_u)`, so switching `sheath_u` off also switches the wall relaxation off — this is a true baseline.

**What it decides.** If the control also dies around step 200, the ramp is the limit and we have been debugging something that was always there. If it sails through to `tstep = 1e0`, the boundary condition owns the limit and Runs B and C are the tests that matter. Either answer saves a lot of time.

### 9.2 Run B — pin the per-step gain

```fortran
 sheath_u_relax      = 0.01
 sheath_u_relax_time = 0.d0     ! disables the tstep/tau form
```

Everything else as-is, `sheath_wall_vel = 4.d-3` kept.

**Why.** [§7.3](#73-mechanism-c-the-clip-makes-the-row-explicit) shows the instability is governed by the *per-step* gain, not by a physical response time. `relax = 0.01` is the value that survived a full 100-step block; `relax = 0.1` died in two. Pinning it at 0.01 keeps the same gain at `tstep = 1e-1`, `3e-1`, `1e0` and beyond.

This reverses my earlier reasoning, and the reversal is the point: expressing the relaxation as a physical time is right when you are damping a *physical* loop, but this instability is numerical, and numerical instabilities are governed by per-step amplification factors. The knob has to match the thing being damped.

**Caveat to watch.** At `tstep = 10` and `relax = 0.01`, the effective response time is 1000 JOREK time units — the BC becomes slow to respond. In a steady state that is harmless, since `u` has already converged. If you see the potential lagging the temperature during a transient, that is the reason.

### 9.3 Run C — scan the drive

```fortran
 sheath_wall_vel = 1.2d-2      ! 3x the current value
```

**Why.** We went from 0 straight to `4.d-3` and got a factor of ~3000 in amplitude. We have no idea where the optimum is, or whether the response is still steep. If the cycle-200 amplitude drops again — say 1.5 → 0.5 — the mechanism is confirmed twice over and there is a value that holds the stripe below the noise floor.

**What to watch.** ψ at the boundary must not drift. Too fast a wall and the equilibrium stops being pinned at all: check that the boundary flux and the plasma current are not wandering. Remember the parameter is **mesh-dependent** ([§6.5](#65-the-thin-resistive-wall)).

**Expectation.** Runs B and C address different halves of the same instability — growth rate and drive — so they should compose. That combination is what I would take into a production run.

### 9.4 Before any long campaign — check the drift direction

As soon as a run survives into the `1e-1` / `3e-1` blocks, spend an hour on this. It is cheap and it is decisive.

1. **Magnitude.** Is $e\Phi/kT_e \approx \Lambda \approx 2.4$ wherever the target current is small? The diagnostic prints this directly.
2. **Sign along the target.** Is $\Phi$ higher where $T_e$ is higher? It must be, since $\Phi = \Lambda T_e/e$.
3. **Direction of the drift.** Does the poloidal E×B flow in the PFR and around the X-point run **from the outer divertor toward the inner one**?

> [!danger] If the drift runs the wrong way, no amount of runtime produces HFSHD
> It would produce the *opposite* asymmetry — density piling on the low-field side. And the place a sign error would live is exactly [§3.3](#33-sign-conventions--the-thing-that-bites): $\Phi = -F_0 u$, with the code's `u` being minus the reference paper's. We fixed that sign once already, in `a_n`. Verifying the drift direction on a short run de-risks a week of machine time.

Only once those three check out is it worth committing to the long ramp.

> [!note] Correction to earlier advice
> I previously suggested holding at `tstep = 1e-2` and "getting physics now". Given §9.0 that is wrong for this goal: 1000 steps at `1e-2` reaches $10\,\tau_A$, still 150× short of HFSHD. Short runs at `1e-2` are useful for the *checks above* and for nothing else.

### 9.5 If B and C do **not** help — free boundary

If the stripe is neither drive-limited nor rate-limited, the honest conclusion is that a boundary which cannot let flux move cannot host a potential that varies along it. Then:

1. Run STARWALL to produce the response matrices.
2. Switch to a free-boundary equilibrium — PF coil currents plus pressure/current profiles, instead of `R_Z_psi_bnd_file`.
3. Set the wall resistivity.

This is **setup work, not code work**: the sheath BC itself needs no modification. Artola's COMPASS paper is a working recipe, including a grid whose boundary follows the first wall — the same geometry choice `grid_to_wall = .t.` already gives us.

### 9.6 Run D — retry the forward / natural route *(promoted: this may be the answer)*

The forward characteristic imposed as a surface term (`bcs%natural%u`) is the better-conditioned object, and we stopped pursuing it on evidence that turns out to be contaminated.

**Why the ~19-step verdict does not stand.** In that campaign the run with `sheath_flux_sign = 0` — which multiplies the sheath surface term by zero, switching it off entirely — **also died at ~19 steps**. Whatever was unstable, it was not the sheath term. The other things active in those runs were `natural%w` and `natural%zj`, which are mis-linearised ([§7.1](#71-mechanism-a-the-naturalw-and-naturalzj-terms-were-mis-linearised)). They are now **refused by `initialise_parameters.f90:393`**, so a `natural%u`-only configuration has never actually been run. `sheath_wall_vel` did not exist then either.

**Why it may succeed where the nodal route stalls.** Mechanism C — the per-step gain limit that is blocking `tstep = 1e-1` and above — is a property of the *nodal* path alone. It arises because the clip sets `sh_dr = 0` and strands `rho` and `zj` at the old timestep. The natural path:

- **never clips on the ion side**, because $\partial j/\partial\Phi \propto e^{-X} \to 0$ smoothly at saturation;
- **carries a full Jacobian at every Gauss point** — `amat(var_u, var_u/var_zj/var_rho/var_Ti/var_Te)` are written unconditionally;
- enforces the condition **weakly**, integrated over the edge, rather than by collocation.

So there is no explicit branch and no per-step gain limit. Given [§9.0](#90-the-timescale-constraint--this-sets-the-priority) — that HFSHD needs `tstep = 1`–`10` and the nodal route currently fails at `1e-1` — **this route may be the one that actually reaches the physics.**

**It does not need free boundary.** The two are orthogonal: the natural BC changes *how the sheath condition is imposed on the `u` equation*; free boundary changes *how ψ behaves at the wall*. But it hits **the same flux dragging** ([§7.2](#72-mechanism-b-flux-dragging--the-main-event)) — that mechanism depends only on `u` varying along the wall, not on how `u` got there — so it needs the same relief.

> [!warning] Code change required, and now made
> `sheath_wall_vel` was nested inside `if ( apply_sheath_u )`, so it was **inert with `natural%u`** — and validation makes `natural%u` and `sheath_u` mutually exclusive (`initialise_parameters.f90:418`), correctly, since they are two implementations of one condition. The combination that should work was therefore unreachable.
>
> The wall block has been hoisted out of that guard and now runs for either route, gated on `(apply_sheath_u .or. apply_natural_u)`, with the enclosing geometry block widened to `apply_cs .or. apply_sheath_u .or. (apply_natural_u .and. sheath_wall_vel > 0)`. The Mach1 section inside is separately guarded by `apply_cs`, so nothing else changes. `index_node` is re-set inside the block, since the Mach1 section can leave it pointing at a derivative DOF.

**Namelist:**

```fortran
 bcs(1)%natural%u     = .true.     ! the sheath BC, forward form
 bcs(1)%dirichlet%u   = .false.    ! u is FREE and set by charge continuity
 bcs(1)%dirichlet%w   = .true.     ! KEEP
 bcs(1)%dirichlet%zj  = .true.     ! KEEP
 bcs(1)%mach1         = .true.
 bcs(1)%sheath_u      = .false.    ! mutually exclusive with natural%u
 ! ... identical for 4, 5, 9 ...

 bcs(2)%dirichlet%u   = .true.     ! MUST keep at least one type pinning u
 bcs(3)%dirichlet%u   = .true.

 bc_natural_open      = .true.
 sheath_Lambda_local  = .true.
 sheath_min_bn        = 0.2
 sheath_X_min         = -3.0
 sheath_smooth_dX     = 0.5
 sheath_wall_vel      = 4.d-3      ! now actually active on this route
 sheath_ramp_time     = 0.d0       ! NB: must be >> tstep to mean anything; 0 = off
 sheath_stiff_max     = 0.d0
```

> [!caution] Do **not** set `natural%w` or `natural%zj`
> Validation will refuse them. That refusal is the whole point of this retry.

**Longer term**, making them usable requires lifting the restriction in [§3.2](#32-the-two-ways-to-impose-a-boundary-condition): an extra trial index carrying `direction_perp(1)` with the correct Hermite basis function (zero on the edge, unit normal derivative) in `mod_boundary_matrix_open`. That would benefit any future natural BC whose flux involves a normal derivative, not just this one. It is not needed for Run D.

### 9.7 Small build improvement

Add to `Makefile.inc` (Intel):

```make
FLAGS += -g -traceback
```

The SIGSEGV backtrace we had to diagnose by inspection carried no line numbers. This costs nothing at runtime.

---

## 10. Testing without a cluster

`util/sheath_bc_unit_test/` exists because the laptop has no MPI toolchain and a cluster round-trip is expensive.

| Script | What it checks |
|---|---|
| `run_test.sh` | **Physics V0**: Λ against its analytic value; $j = 0$ at the floating potential; $\Phi_{\rm float} = 51.59$ V at $T_e = 20$ eV; electron saturation at $\Phi = 0$; **12 finite-difference derivative checks** agreeing to ~1e-10; limiter C¹, monotone and finite over $X\in[-200,200]$. |
| `syntax_check.sh` | Stub-compiles the sheath sources in **4 configurations** of `with_TiTe` × `with_vpar`. |
| `decl_check.py` | Per-scope duplicate-declaration check. Needed because `mod_boundary_conditions.f90` cannot be stub-compiled, and this class of bug (`sh_c` colliding with an existing `sh_C` — Fortran is case-insensitive) costs a full cluster build to discover. |

**Naming collisions caught this way, all real:** `Bdotn` vs the existing `bdotn`; the routine `sheath_lambda` vs the namelist variable `sheath_Lambda`; `sh_c`/`sh_s` vs the existing `sh_C`. In a codebase this size with case-insensitive identifiers, this check pays for itself.

---

## 11. Quick reference

### 11.1 Namelist — current best configuration

```fortran
 ! --- sheath BC on the strike-point boundary types
 bcs(1)%sheath_u      = .true.
 bcs(1)%dirichlet%u   = .true.     ! the sheath row overwrites it; KEEP
 bcs(1)%dirichlet%w   = .true.     ! KEEP
 bcs(1)%dirichlet%zj  = .true.     ! KEEP
 bcs(1)%mach1         = .true.     ! j_sat assumes it
 ! ... identical for 4, 5, 9 ...

 bcs(2)%dirichlet%u   = .true.     ! at least one type must pin u
 bcs(3)%dirichlet%u   = .true.

 sheath_Lambda        = 3.0
 sheath_Lambda_local  = .true.
 sheath_V_wall        = 0.0
 sheath_u_exp_max     = 2.0
 sheath_u_exp_min     = -3.0
 sheath_min_bn        = 0.2
 sheath_u_align_psi   = .false.
 sheath_u_value_only  = .false.
 sheath_wall_vel      = 4.d-3

 ! --- Run B changes these two:
 sheath_u_relax       = 0.01
 sheath_u_relax_time  = 0.d0
```

### 11.2 The three things to remember

1. **$\Phi = -F_0 u$** in model600 — the code's `u` is minus the reference paper's.
2. **A natural BC on a variable that also has a Dirichlet does nothing**, because surface integrals only reach rows the Dirichlet overwrites.
3. **A potential varying along the wall drags poloidal flux.** This is the whole difficulty, it is specific to codes that evolve **B**, and it is why a resistive wall — local (`sheath_wall_vel`) or proper (STARWALL) — is part of the answer rather than an optional extra.

### 11.3 File map

```
models/model600/mod_sheath_bc.f90          the characteristic (single source of truth)
models/model600/mod_sheath_diag.f90        wall-current diagnostic, psi0 storage, u init
models/model600/mod_boundary_conditions.f90   nodal path + resistive wall  <- the working one
models/model600/mod_boundary_matrix_open.f90  natural/surface-term path
models/phys_module.f90                     parameter declarations + doc comments
models/preset_parameters.f90               defaults (all off)
models/model600/initialise_parameters.f90  namelist + validation
communication/broadcast_phys.f90           MPI pack/unpack
util/sheath_bc_unit_test/                  offline physics + syntax + declaration checks
```

---

## 12. ⚠ Build / branch discrepancy — resolve this first

The namelist used on the cluster sets

```fortran
 keep_current_prof          = .true.
 keep_current_prof_confined = .true.
 keep_current_psin_cutoff   = 0.995
 keep_current_psin_sig      = 0.005
 keep_current_z_sig         = 0.02
```

but **none of the last four exist on `bc-tests`**. They live on `keep-current-prof-cutoff`, which is *not* an ancestor of `bc-tests`. Since `models/model600/initialise_parameters.f90:257` reads the namelist as

```fortran
read(42,in1)      ! no iostat
```

a build from plain `bc-tests` would **abort at namelist read** on an unrecognised group member. The runs therefore came from a build that merges `bc-tests` with `keep-current-prof-cutoff` — not from `bc-tests` as it stands.

**Why this matters, twice over:**

1. **Reproducibility and review.** The tree reviewed in this document is not exactly the tree that ran. The sheath sources are almost certainly identical, but that should be confirmed rather than assumed.
2. **HFSHD directly.** `keep_current_psin_cutoff = 0.995` is what masks `keep_current_prof` off outside the separatrix. Without it, the artificial $\eta(j - j_0)$ current source acts in the SOL and PFR — which is precisely what faked the $\eta$-scaling of HFSHD before, and which would also fight the sheath-driven boundary currents this whole boundary condition exists to create. `T_min_eta`, from the same branch, is likewise absent here.

**Action:** merge `keep-current-prof-cutoff` into `bc-tests` (or record explicitly which merge the cluster executable was built from) before the next round of runs.

---

## 13. Plan: making the natural (forward) sheath BC work robustly

This section takes the seven recommendations (11–17) received externally, audits each against what is actually in the code, and turns the gaps into an ordered plan. No nodal path, no new Dirichlet rows.

### 13.0 The headline finding: the Newton solver has never been switched on

`use_newton` defaults to `.false.` and the namelist does not set it. **Every run in this campaign solved a strongly nonlinear boundary condition with one linearised solve per timestep and no residual check whatsoever.** If the linearisation is poor — which it is whenever $X$ moves appreciably within a step — the step is simply wrong, and nothing detects it.

JOREK already has an inexact Newton solver, `solvers/mod_newton_solver.f90`, reached with `use_newton = .true.`:

| Feature | Where |
|---|---|
| Up to `maxNewton = 20` iterations | `newton_loop` |
| Convergence on $\lVert R_k\rVert/\lVert R_0\rVert \le$ `iter_tol` | `tol_Newton` |
| **Matrix re-assembled every iteration** — the sheath Jacobian is re-evaluated at the updated state | `construct_matrix` inside the loop |
| Eisenstat–Walker forcing, `tol_Gmres = gamma_Newton*(‖R_k‖/‖R_{k-1}‖)^alpha_Newton` | `tol_Gmres` |
| Step rejection via `solver%step_success = .false.` | on `inewton == maxNewton` |

That is recommendation 17's machinery, already written. Two caveats, both important:

> [!warning] What Newton does *not* give you here
> - **No line search / damping.** There is no $\lambda$. That half of recommendation 17 is genuinely absent, and adding it would mean modifying `solve_newton`.
> - **On failure JOREK aborts, it does not retry.** `jorek2_main.f90:789–797` prints `NO CONVERGENCE ... ABORTING` and exits the time loop; the adaptive-timestep block immediately below is commented out (*"in progress..."*). So treat a Newton abort as the signal to reduce `tstep` by hand.
> - **`gmres_tol` changes meaning.** `solver%iter_tol = gmres_tol`, and in the Newton path that is the *nonlinear* tolerance; the linear tolerance becomes adaptive. Your `gmres_tol = 1.d-7` therefore asks Newton for a 10⁻⁷ residual reduction within 20 iterations.

This may also account for a good part of Mechanism C. "Crashes a few steps after a timestep jump" is precisely the signature of one linearisation per step: bigger step → bigger state change → worse linearisation → nothing checks it.

### 13.1 Audit of recommendations 11–17

| # | Recommendation | Status in the code | Action |
|---|---|---|---|
| 11 | Smooth cap on the electron side, chain factor in the Jacobian | **Already done** — `sheath_x_limited` is a softplus, C¹ and monotone, and `fp = expmx * dxlim_dx` carries `dX_lim/dX` into residual *and* Jacobian. Asymptotic branches keep the debug traps alive. | none |
| 11 | `expm1` near $X = 0$ | **Was missing.** `f = 1 - exp(-X)` lost ~4 digits at $X = 10^{-4}$ — and $X \approx 0$ *is* the floating condition, i.e. the regime of interest. | **fixed** (§13.2) |
| 12 | Positive effective $\rho, T_i, T_e$ with a consistent Jacobian | **Already done, and stronger than asked** — `mod_boundary_matrix_open` applies `dcorr_neg_dens_drho1(r0)`, `dcorr_neg_temp_dT1(Ti0)`, `dcorr_neg_temp_dT1(Te0)` as chain factors. `corr_neg` is a smooth correction, not a hard `max`. The internal floors (`1e-14`, `1e-10`) zero their own sensitivities consistently and sit far below `corr_neg`'s floor. | monitor via the diagnostic |
| 13 | Reuse Mach-1 conventions exactly | **Already done** — `sh_Bn = bdotn*Btot` uses the very `bdotn` computed for the particle/heat fluxes, and `g_bn = normal_sign*factor` is the Chodura–Riemann function from `vpar_smoothing_coef`. Same normal, same $c_s$, same density normalisation. | none |
| 13 | Avoid raw `sign()` in a constitutive law near $b_n = 0$ | **A genuine defect.** `sheath_min_bn > 0` applies `g_eff = sgn(b_n)*max(|g|, g_min)` — which jumps from $-g_{\min}$ to $+g_{\min}$ across tangency. | **fixed** (§13.2) |
| 14 | No conflicting strong BC on the same boundary | **Already enforced** by `initialise_parameters`: `natural%u` requires `dirichlet%u = .false.`, refuses `sheath_u`, refuses `natural%w`/`%zj`, and requires ≥1 type retaining `dirichlet%u` as the gauge. `dirichlet%zj` is *not* a conflict — the frozen `zj` trace cancels between the volume term and the surface flux. | **one open question**, §13.6 |
| 15 | Ion saturation leaves the potential weakly constrained | Correct and unavoidable. The gauge exists; we are not going back to the inverted law. | none |
| 16 | Continuation ramp on residual **and** Jacobian | **Mechanism already correct** — `sheath_ramp` multiplies `rhs_ij(var_u)` and every `amat(var_u,*)` entry, at the same Gauss point. | **the value was always wrong**, §13.3 |
| 17 | Damped Newton, step rejection, convergence-driven timestep | Newton + rejection exist (`use_newton`); **no line search**; abort rather than retry. | §13.0, §13.4 |

**Why the grazing floor is not needed at all in the forward form.** The term added to the `u` row is $-\oint v\,R\,(zj_{\rm sh} - zj)\,(\mathbf{B}\cdot\mathbf{n})\,d\ell$. As the field goes tangent, $\mathbf{B}\cdot\mathbf{n}\to 0$ and the whole term vanishes *whatever* $j_{\rm sat}$ does; `dzj_du` vanishes with it, so the Robin diagonal goes smoothly to zero and `u` is left to the vorticity equation. That is also the physics — a tangential field delivers no parallel flux, so there is no sheath. The floor exists only for the **nodal** path, where `u` is slaved to the characteristic and a vanishing $j_{\rm sat}$ makes the row singular. And no continuous function can floor a magnitude while preserving a sign that changes, so the floor is *inherently* discontinuous: it must go, not be smoothed.

### 13.2 Phase 0 — code changes (done)

| Change | File | Why |
|---|---|---|
| Wall block hoisted out of `if (apply_sheath_u)`, gated on `(apply_sheath_u .or. apply_natural_u)`; enclosing geometry guard widened | `mod_boundary_conditions.f90` | `sheath_wall_vel` was inert with `natural%u`, so forward-BC + resistive wall was unreachable |
| Accurate $1-e^{-X}$ below $\lvert X\rvert = 10^{-4}$ (4-term series, rel. error ~10⁻¹⁸) | `mod_sheath_bc.f90` | rec 11; the floating condition is exactly where the cancellation bites |
| `sheath_min_bn > 0` now **refused** with `natural%u`, with the reasoning in the message; the existing `NOTE` for `= 0` rewritten to say it is correct rather than a compromise | `initialise_parameters.f90` | rec 13; removes a discontinuous constitutive law |

Verified offline: 18/18 physics and finite-difference derivative tests pass (worst error 1.2e-9), `if`/`do` balance 0, declaration check clean on all four sheath files, 4-configuration stub compile clean, no new lines over 132 characters. **Not compiled with a real toolchain** — that needs the cluster.

### 13.3 Phase 1 — the Run D namelist, line by line

```fortran
 ! --- route
 bcs(1)%natural%u     = .true.     ! forward characteristic as a surface term
 bcs(1)%dirichlet%u   = .false.    ! u is FREE, solved by the vorticity equation
 bcs(1)%dirichlet%w   = .true.     ! KEEP - its surface term reaches no row a Dirichlet doesn't
 bcs(1)%dirichlet%zj  = .true.     ! KEEP - the frozen trace cancels exactly
 bcs(1)%mach1         = .true.     ! j_sat reuses the Mach 1 outflow
 bcs(1)%sheath_u      = .false.
 ! ... identical for 4, 5, 9 ...
 bcs(2)%dirichlet%u   = .true.     ! gauge - see 13.6
 bcs(3)%dirichlet%u   = .true.
 bc_natural_open      = .true.

 ! --- characteristic
 sheath_Lambda        = 3.0        ! Lambda_0; <=0 computes it from central_mass
 sheath_Lambda_local  = .true.     ! Lambda(Ti/Te), consistent with c_s: 2.4 not 3.0
 sheath_V_wall        = 0.0        ! grounded vessel
 sheath_X_min         = -3.0       ! electron saturation at Phi = 0
 sheath_smooth_dX     = 0.5
 sheath_min_bn        = 0.0        ! REQUIRED now; see 13.1
 sheath_flux_sign     = 1.0

 ! --- robustness
 sheath_stiff_max     = 1.0        ! keep the Robin diagonal <= the polarisation diagonal
 sheath_ramp_time     = 1.0        ! full strength at ~step 200; see below
 sheath_init_u        = .true.     ! start u at the floating potential, not at 0
 sheath_wall_vel      = 4.d-3      ! flux relief - now active on this route

 ! --- nonlinear solve
 use_newton           = .true.
 maxNewton            = 20
 gamma_Newton         = 0.5
 alpha_Newton         = 2.0
 gmres_tol            = 1.d-5      ! NB: this is the NEWTON tolerance when use_newton = .true.
```

**`sheath_ramp_time = 1.0`, and why the old value did nothing.** The ramp is $(t_{\rm now}-t_{\rm start})/\texttt{sheath\_ramp\_time}$, clipped to $[0,1]$. With your schedule the elapsed time is 0.1 after 100 steps and 1.1 after 200. A `ramp_time` of 100 therefore gave a ramp of $\sim10^{-3}$ — **the sheath term was effectively never switched on in any V4 run.** At 1.0 the ramp is 10% around step 110, 50% around step 150, and full by step 200: a genuine continuation, completed before the timestep grows.

**`sheath_stiff_max = 1.0` matters more than it looks.** Uncapped, the Robin diagonal $R\lvert\partial zj/\partial u\rvert\lvert B_n\rvert\,d\ell\,\theta\,\Delta t$ can exceed the row's own polarisation diagonal $\rho R^3$ by three orders of magnitude at a cold dense target. At that point the natural BC has silently become a pointwise Dirichlet — and you are back to the nodal pathology, with node-to-node $T_e$ noise imprinted on `u` and thus on the flux-dragging velocity. Capping it scales residual and Jacobian by the same $\alpha$, so the fixed point is untouched.

**`gmres_tol = 1.d-5`.** With `use_newton = .true.` this becomes the nonlinear tolerance. 10⁻⁷ within 20 iterations is a demanding ask for a stiff boundary condition, and failure means an abort rather than a retry. Start at 10⁻⁵; tighten later if Newton converges in a handful of iterations.

### 13.4 Phase 2 — staged runs

**D0 — instrumented control.** Everything above, but `sheath_flux_sign = 0.0`. The characteristic is evaluated and the diagnostic reports it, but the term contributes nothing to residual or matrix.

*Purpose:* (i) confirm Newton converges on the base problem, so any later failure is attributable to the sheath; (ii) read off `ePhi/kTe`, `max|j/jsat|` and `e-limited %` for the *unperturbed* state, which tells you how far the plasma is from the characteristic before the BC pushes it; (iii) this is the control the natural route never had — the old V4 `flux_sign = 0` run crashed at ~19 steps, but with `natural%w`/`%zj` on. **Expect it to run clean now.** If it does not, stop: the problem is not the sheath.

**D1 — full run.** `sheath_flux_sign = 1.0`, ramp as above. Watch the `SHEATH:` line every step through the ramp.

**D2 — push the timestep.** Only after D1 is stable through the `1e-1` block. Extend `nstep_n` at `1e-1`, then let it ramp. §9.0 says you need `tstep = 1`–`10` for HFSHD; this is where you find out whether the forward route gets there.

### 13.5 Phase 3 — what to watch, and the decision rules

The diagnostic already prints everything needed:

```
SHEATH: I_wall=... A (Ampere ... A)  ePhi/kTe min/mean/max=... max|j/jsat|=...  e-limited ... %
```

| Signal | Healthy | What it means if not |
|---|---|---|
| `ePhi/kTe` mean | → $\Lambda \approx 2.4$ | The BC is not reaching its fixed point — check the ramp, then `stiff_max` |
| `I_wall` vs `I_Ampere` | converging toward each other | The plasma is not delivering the current the sheath asks for; the mismatch is what drives `u` |
| `I_wall` total | → 0 for a grounded wall in steady state | Net charge is leaving the domain |
| `e-limited %` | small and **not growing** | The electron-saturation cap is load-bearing rather than protective — a physical electron-saturation model would be needed |
| `max|j/jsat|` | O(1) | Values ≫1 mean the plasma is far off the characteristic |
| Newton iterations | few, stable | Rising count is the early warning that used to be invisible |

**Decision rules.** Newton abort → halve `tstep` at that block and restart from the last output; it is not a code failure, it is the rejection working. `e-limited %` climbing steadily → the wall is being driven to electron saturation, which is physics, and `sheath_X_min` becomes a modelling choice. Boundary structures returning in `zj`/`w` → scan `sheath_wall_vel` upward as in Run C; the forward route reduces the drive but does not remove flux dragging.

### 13.6 Open questions and what is still missing

1. **What are boundary types 2 and 3?** They currently carry `natural%rho = .t.`, `mach1 = .f.`, and would carry the `dirichlet%u` gauge. If they are far-field or open ends, `u = 0` there is fine. **If they are physical wall, it is wrong** — `u = 0` pins the *plasma* potential at the sheath entrance to zero, whereas the sheath says it should be $\approx\Lambda T_e/e$. In that case they need `natural%u` too, and the gauge has to shrink to a single reference node.
2. **The floating-conductor machinery of recommendation 14 is not needed here.** It applies to a wall that floats. The AUG vessel is grounded, so `sheath_V_wall = 0` with one Dirichlet reference is physics, not a gauge hack. Revisit only if you model a floating component.
3. **No line search.** The one part of recommendation 17 with no existing implementation. If Newton proves erratic, adding a damping factor $\lambda$ to `solve_newton` is a contained change — but try it unmodified first.
4. **No automatic timestep reduction on Newton failure.** The block exists in `jorek2_main.f90` but is commented out. Reviving it would turn every abort into a retry, which is the single most valuable robustness addition after Newton itself.
5. **`natural%w` / `natural%zj`** remain refused. Lifting that needs the `direction_perp(1)` trial index with the correct Hermite basis (§9.6). Not required for Run D.

---

## 14. Audit of the six specific concerns

Each checked against the code, not reasoned about abstractly. **Four are clean, one needed three fixes, and one is a real open problem.**

### 14.1 Orientation applied twice in `g_bn`, `zj_sh`, `sh_Bn`? — **No. Applied exactly once.**

```fortran
bdotn       = (+ ps0_y*normal(1) - ps0_x*normal(2)) / x_g(ms) / Btot   ! signed
normal_sign = sign(1.d0, bdotn)
factor      = 0.25*(1 + tanh((abs(bdotn) - c_1)/c_2))**2 - c_3         ! function of |bdotn| only
g_bn        = normal_sign * factor
sh_Bn       = bdotn * Btot
```

`factor` depends on `abs(bdotn)`, so it is a **positive magnitude** carrying no orientation; `normal_sign` supplies the sign once. Then in `sheath_current`, `zj_sat = c_sat*rho*g_bn*cs/Btot` with `c_sat = -a_n/2 > 0`. The product that enters the residual is therefore

$$zj_{\rm sh}\cdot \texttt{sh\_Bn} \;\propto\; \underbrace{\mathrm{sign}(b_n)\,|g|}_{g_{bn}}\cdot\frac{1}{|B|}\cdot \underbrace{b_n |B|}_{\texttt{sh\_Bn}} \;=\; |g|\,|b_n|\,f$$

— manifestly **sign-definite in the geometry**, with the sign carried by $f = 1-e^{-X}$ alone: $f>0$ is ion current out of the plasma, $f<0$ is electron current in. That is correct on both targets, and it is what you want: the orientation cancels rather than accumulating. ✅

### 14.2 Is `zj_sh - zj0` validated against the exact strong-form volume term? — **No, and my derivation disagrees with the implemented sign.**

> [!danger] This is the one real open problem
> The `u` equation's current term in `mod_elt_matrix_fft.f90:1571` is
> ```fortran
> + v * (ps0_s * zj0_t - ps0_t * zj0_s) * tstep * factor(var_u,2)
> ```
> The derivatives are on **`zj`**, not on the test function `v`. This is the **strong form**. It is assembled directly and **generates no surface term at all** — so nothing is being "put back". That contradicts the derivation in the `mod_sheath_bc.f90` header, which reasons from a weak form that integrates by parts and produces $-\oint v R\,zj\,(\mathbf{B}\!\cdot\!\mathbf{n})\,d\ell$. **That weak form is not what the code assembles.**

Re-deriving from what *is* assembled. With JOREK's $\mathbf{B}_{\rm pol} = \nabla\psi\times\nabla\varphi$, the identity is $[\psi,zj] = R\,\nabla\!\cdot\!(zj\,\mathbf{B}_{\rm pol})$, and `ps0_s*zj0_t - ps0_t*zj0_s` $= [\psi,zj]\cdot\texttt{xjac}$. So the volume term is $+\int v\,R\,\nabla\!\cdot\!(zj\,\mathbf{B}_{\rm pol})\,dV$. Applying the divergence theorem *to that*:

$$\int v R\,\nabla\!\cdot\!(zj\,\mathbf{B}_{\rm pol})\,dV \;=\; \oint v R\,zj\,B_n\,dS \;-\; \int \nabla(vR)\cdot(zj\,\mathbf{B}_{\rm pol})\,dV$$

To make the outflow $zj_{\rm sh}$ instead of $zj$, replace the first piece — which means **adding**

$$+\oint v\,R\,(zj_{\rm sh} - zj)\,B_n\,dS$$

to the same RHS. The code has

```fortran
rhs_ij(var_u) = - v * BigR * ( zj_sh - zj0 ) * sh_Bn * dl * tstep * sheath_ramp
```

i.e. the **opposite sign**, *provided* `sh_Bn = +B·n`. And that proviso is exactly what I cannot settle on paper: `bdotn` evaluates to $(\psi_Z n_R - \psi_R n_Z)/(R|B|)$, which is $-\mathbf{B}\!\cdot\!\mathbf{n}$ under the reference-paper convention but $+\mathbf{B}\!\cdot\!\mathbf{n}$ once you fold in that the code's $\psi$ is minus the paper's ([§3.3](#33-sign-conventions--the-thing-that-bites)). Two sign conventions compose here and I will not guess which wins.

> [!success] But there is a clean empirical test, and it costs one namelist line
> `sheath_flux_sign` was built for exactly this: it multiplies the surface term (residual **and** Jacobian) by $\pm1$ without a rebuild. Run D1 twice, `+1` and `-1`, and read the diagnostic:
> ```
> SHEATH: I_wall=... A (Ampere ... A)  ePhi/kTe min/mean/max=...
> ```
> `I_wall` is $\oint zj_{\rm sh}R B_n$, the current the sheath *wants*; `I_Ampere` is $\oint zj\,R B_n$, what the plasma actually delivers. **The correct sign is the one where `I_Ampere` moves toward `I_wall`.** The wrong sign drives them apart — an anti-damped boundary condition, which is precisely the failure signature this campaign has been chasing.
>
> Note the earlier "sign flip" test is *not* evidence: those V4 runs also had `natural%w`/`%zj` on, so both signs died at ~19 steps for an unrelated reason ([§9.6](#96-run-d--retry-the-forward--natural-route-promoted-this-may-be-the-answer)).

A second consequence: the header's claim that *"the frozen `zj` trace cancels exactly between the strong-form volume term and the added surface flux"* also rests on the weak-form picture. With the strong form, boundary `zj` enters the volume integral directly through `zj0_s`/`zj0_t` at boundary Gauss points while being Dirichlet-frozen. Nothing cancels; the surface term simply adds. **Treat that header paragraph as unverified.**

### 14.3 The local stiffness limiter — **fixed the kink; the omitted derivative is harmless**

Two separate criticisms, with different answers.

> [!danger] RETRACTED 2026-08-18 — this claim was wrong
> I argued that scaling the sheath term by $\alpha$ leaves the fixed point untouched, because the root of $-\alpha(zj_{\rm sh}-zj_0)B_n$ is $zj_{\rm sh}=zj_0$ for any $\alpha\neq0$. **That is not the fixed point of the coupled system.** The `u` row is $R_{\rm bulk} + \alpha R_{\rm sh} = 0$, whose solution satisfies $R_{\rm bulk} = -\alpha R_{\rm sh}$ — it does not force $R_{\rm sh}=0$. The imposed wall current is effectively $(1-\alpha)\,zj + \alpha\,zj_{\rm sh}$: a weighted average of what the plasma delivers and what the sheath asks for. **$\alpha<1$ therefore changes the physical solution**, not just the transient. My argument tacitly assumed the system independently drives the sheath mismatch to zero, which holds only at an exact steady state with vanishing perpendicular current.
>
> Consequence: `sheath_stiff_max` is a **physics modification, not a solver control**. Use `sheath_stiff_max = 0` (off) for any correctness work, and handle stiffness through timestep, line search and preconditioning instead. If it is ever needed, it must be treated as regularisation: frozen during each nonlinear solve or differentiated consistently, its value output per Gauss point, and the result shown to converge as $\alpha\to1$.

**"It has an omitted derivative."** True: the exact Jacobian is $-[\alpha'(zj_{\rm sh}-zj_0) + \alpha\,zj_{\rm sh}']B_n$ and the code has only the second piece. But the missing term is **proportional to the residual itself**, and $(zj_{\rm sh}-zj_0)\to0$ at the solution — so it vanishes there and Newton keeps its local convergence rate. This one is genuinely benign.

**What was *not* benign** was the hard switch `if (sh_d_robin > stiff_max*sh_d_pol)`, which put a kink in the residual as a function of the state — costing Newton iterations for nothing. Replaced with a smooth cap having the same two limits:

$$\alpha = \frac{1}{1+r},\qquad r = \frac{d_{\rm robin}}{\texttt{sheath\_stiff\_max}\cdot d_{\rm pol}}$$

$\alpha\to1$ for $r\ll1$, $\alpha\to \texttt{stiff\_max}\,d_{\rm pol}/d_{\rm robin}$ for $r\gg1$, C<sup>∞</sup> throughout. Note the recalibration: at $r=1$ the old form gave $\alpha=1$, the new gives $\alpha=0.5$, so damping starts a little earlier.

### 14.4 Safe exponential handling in `sheath_current` — **already sound; one inconsistency fixed**

Branch analysis with `sheath_exp_max = 30`:

| Branch | `expmx` | Reachable? | Consistent? |
|---|---|---|---|
| `x_lim > 30` | `0` | yes (ion saturation) | ✅ true derivative $e^{-X}\approx0$ too |
| `-30 ≤ x_lim ≤ 30` | `exp(-x_lim)` ≤ 1e13 | yes | ✅ exact |
| `x_lim < -30` | `exp(30)` frozen | only if `sheath_X_min < -30` | ❌ **was inconsistent** |

The softplus guarantees `x_lim ≥ sheath_X_min`, so with the default `-3.0` the third branch is unreachable — but if a user disables the limiter with a very negative `X_min`, the residual is frozen there while `fp = expmx*dxlim_dx` still reported a derivative of $10^{13}$. **Fixed:** an `x_frozen` flag now sets `fp = 0` to match the frozen residual. Overflow itself was never possible — `exp(30) ≈ 1.07\times10^{13}$, and `sheath_x_limited` guards its own `exp(z)` the same way. ✅

### 14.5 Reset `rhs_ij` and `amat` inside their loops — **no bug today, but done anyway**

`rhs_ij = 0.d0` / `amat = 0.d0` were at routine entry (line 82–83), outside every loop. Auditing the assignments:

- every condition guarding a write (`apply_natural_bc(...)`, `with_vpar`, `with_TiTe`, `with_neutrals`) is **loop-invariant**, so each component that is ever written is rewritten on every iteration before use;
- `rhs_ij` is consumed in the same `(i,j,im)` iteration it is built (line 473), and the RHS loop skips any variable without `apply_natural_bc`;
- `amat` is never accumulated onto itself — the only self-accumulation anywhere is `rhs_ij(var_rhon)`, which builds on an assignment made earlier in the *same* iteration.

So there is no stale-value bug. **But the invariant is undocumented and one future component added under an index-dependent condition would leak silently.** Resets moved inside the `im` loop (for `rhs_ij`) and the `in` loop (for `amat`), with a comment saying why. Cost: `n_var` and `n_var²` stores on boundary elements only. ✅

### 14.6 FD Jacobian and charge-conservation tests — **half exists; here is the rest**

**Constitutive level: exists and passes.** `util/sheath_bc_unit_test/run_test.sh` finite-differences `sheath_current` in three regimes × four variables, 12 checks, worst error 1.2e-9, plus Λ, floating potential, electron saturation, and limiter C¹/monotonicity over $X\in[-200,200]$. 18/18 pass after today's changes.

**Assembly level: verifiable by inspection, and it checks out.** With $\texttt{rhs}_u = -vR(zj_{\rm sh}-zj_0)B_n\,d\ell\,\Delta t\,\alpha$ and the convention $\texttt{amat} = -\theta\,\partial\,\texttt{rhs}/\partial x$:

| | expected | code (lines 515–522) |
|---|---|---|
| $\partial/\partial u$ | $+vR\,\texttt{dzj\_du}\,\psi B_n d\ell\,\theta\Delta t\,\alpha$ | ✅ matches |
| $\partial/\partial zj$ | $-vR\,\psi B_n d\ell\,\theta\Delta t\,\alpha$ | ✅ matches |
| $\partial/\partial\rho,T_i,T_e$ | $+vR\,\texttt{dzj\_d}\ast\,\psi B_n d\ell\,\theta\Delta t\,\alpha$ | ✅ matches |
| single-$T$: $\partial/\partial T$ | $\tfrac12(\texttt{dzj\_dTi}+\texttt{dzj\_dTe})$, since `Ti0 = 0.5*T0`, `Te0 = Ti0` | ✅ matches (line 522) |

**Charge conservation: already instrumented, and it is the key runtime test.** `I_wall` vs `I_Ampere` per boundary type plus totals. Two things to watch:
1. **Sign** — §14.2. The correct `sheath_flux_sign` is the one where they converge rather than diverge.
2. **Magnitude** — if a factor of $R$, $1/R$ or $F_0$ were wrong in the surface term, they would converge to a fixed *ratio* $\neq1$ rather than to each other. That ratio, if it appears, names the missing factor directly.

For a grounded wall in steady state, `I_wall` should also approach 0 overall.

### 14.7 Summary of changes made in this pass

| Change | File | Concern |
|---|---|---|
| Smooth stiffness cap $\alpha = 1/(1+r)$ replacing the hard switch | `mod_boundary_matrix_open.f90` | 14.3 |
| `rhs_ij` reset in the `im` loop, `amat` reset in the `in` loop | `mod_boundary_matrix_open.f90` | 14.5 |
| `x_frozen` flag → `fp = 0` in the capped branch | `mod_sheath_bc.f90` | 14.4 |

All offline checks pass: 18/18 unit tests, 4-configuration stub compile, declaration check, `if`/`do` balance.

**Not changed, deliberately:** the surface-term sign (§14.2). It should be settled by measurement with `sheath_flux_sign`, not by me picking a side in a two-convention argument — and the measurement is one namelist line.

---

## 15. Correction: the S0 control was invalid, and so was the ramp

**Observed 2026-08-18.** S0 (`sheath_flux_sign = 0.0`) crashed at step 205 with `zj` reaching 369 MA/m² and `w` 1751 in a thin filament along the **LFS divertor leg** — two orders of magnitude worse than the nodal route with `sheath_wall_vel` at a comparable cycle (`zj` 1.53, `w` 16.8).

### What went wrong with the test design

`sheath_flux_sign = 0` zeroes the surface term but leaves `bcs%dirichlet%u = .false.` in force. **`u` at the divertor then has no boundary condition at all.** It is free to develop whatever along-wall gradient the vorticity equation produces, and $v_E\!\cdot\!\mathbf{n} = R\,\partial u/\partial\ell$ is exactly the flux-dragging velocity of [§7.2](#72-mechanism-b-flux-dragging--the-main-event). This is the *maximum* of Mechanism B, not a control for it.

Compare the configurations honestly:

| | `u` at the sheath boundary | flux dragging |
|---|---|---|
| baseline | Dirichlet `u = 0` | none (no along-wall variation) |
| nodal `sheath_u` | slaved to the characteristic | grid-scale, from `Te` noise |
| **S0 as run** | **nothing** | **unbounded** |
| natural `%u`, term active | set by charge continuity | smoothed by the vorticity operator |

### The same flaw was in S+ and S−

With `sheath_ramp_time = 1.0` the term is near zero for the first ~100 steps, so **during the ramp S+ and S− are S0**. The filament grows before the boundary condition ever engages.

> [!danger] The ramp is a trap on this route
> Recommendation 16 assumes you are ramping *a term*. Here the surface term **is** the boundary condition for `u` — ramping it to zero removes the condition rather than softening it. This is the opposite of a continuation: it starts from an ill-posed problem and continues toward a well-posed one.
>
> `sheath_init_u = .true.` already does what the ramp was for: it starts `u` at the floating potential $\Lambda T_e/e$, the fixed point of the characteristic at $j = 0$, so the initial nonlinear mismatch is small. Newton handles the rest. **Use `sheath_ramp_time = 0.d0`.**

Both conditions are now warned about in `initialise_parameters.f90`.

### The corrected run set

| Run | Configuration | Purpose |
|---|---|---|
| **C** | `natural%u = .false.`, `dirichlet%u = .true.` on 1/4/5/9; `use_newton = .t.`; `sheath_wall_vel = 4.d-3` | The **real** control: standard BC. Does the case survive the ramp at all? This is what S0 should have been. |
| **S+** | `natural%u`, `sheath_flux_sign = +1`, **`sheath_ramp_time = 0`**, `sheath_init_u = .t.` | implemented sign |
| **S−** | as S+ with `sheath_flux_sign = -1` | the sign §14.2 derives |

With the ramp gone, the sign test is also *faster* to read: the term acts from step 1, so `sheath_watch.py` should show the `I_wall`/`I_Ampere` gap opening or closing within the first few tens of steps rather than after 150.

### What this does not change

The `zj`/`w` filament on the LFS leg in S0 is **not** evidence against the forward route — the route was not being exercised. It is, however, a clean demonstration of Mechanism B in isolation: remove every constraint on `u` at the wall and the flux dragging is immediate and violent. That is the strongest confirmation yet that the mechanism identified in §7.2 is the controlling one.

---

## 16. ~~Root cause: the characteristic has no solution at grazing incidence~~ (SUPERSEDED by §17)

> [!warning] Partly wrong — read §17
> The diagnosis below was built on `max|j/jsat| = 1.32e3`, which turned out to be measured at
> **gated-off** points and is therefore meaningless. With the gate-aware diagnostic the active-point
> ratio is **1.03**. Solvability was then falsified outright: `sheath_sat_slope = 0.05` makes the
> characteristic unbounded above, so every demanded current is reachable at finite potential — and
> the run was unchanged (crash 32 vs 33, max 7.00 vs 7.26). What survives from this section is the
> smooth gate itself, which is correct and necessary; what does not is the claim that unsolvability
> was driving the runaway.

### 16.0 (original text follows)

**Diagnostic, first 20 records (2026-08-18), sheath on types 1/4/5/9, `wall_vel = 0`, `min_bn = 0`:**

| quantity | behaviour |
|---|---|
| `max\|j/jsat\|` = `\|zj0/zj_sat\|` | **1.32e3**, constant |
| `ePhi/kTe` max | 0.06 → **12.10**, ~+0.53 per record, no saturation |
| `ePhi/kTe` mean | 0.00 → 2.4 (≈ Λ) → **3.23 and climbing** |
| `ePhi/kTe` min | 0.01, pinned (the gauge boundary) |
| `I_wall` | −15540 → −150 A, converging to a **non-zero** mismatch |
| `I_Ampere` | 868 A, constant — `dirichlet%zj` freezes it |

The forward characteristic $f = 1-e^{-X}$ has range $\left(-(e^{-X_{\min}}-1),\,1\right] = (-19,\,1]$. A required ratio of **1320 has no solution for any $u$**. At those points the residual can never vanish, so $u$ is driven monotonically — exactly the linear climb in the max while the mean passes through $\Lambda$ and keeps going.

**Where the ratio diverges:** at grazing incidence. $zj_{\rm sat}\propto g(b_n)\to0$ as the field goes tangent, while the frozen $zj_0$ does not, so $f = zj_0/zj_{\rm sat}\to\infty$.

> [!danger] My removal of `sheath_min_bn` for the natural route was wrong (§13.1, §13.2)
> I argued the forward form needs no grazing floor because the surface term vanishes with $\mathbf{B}\!\cdot\!\mathbf{n}$. That is true of the **forcing** and false of the **solvability**: as $g\to0$ both $zj_{\rm sat}$ and the Robin stiffness vanish, but the $u$ that satisfies the characteristic *diverges*. This is why the nodal path ran 199 steps with `min_bn = 0.2` and 174 with it off.

**The fix — a smooth gate, not a floor.** Flooring $|g|$ while keeping $\mathrm{sign}(b_n)$ is discontinuous through tangency (§14.1). Instead the whole term, residual **and** Jacobian, is multiplied by

$$w(b_n) = \frac{b_n^2}{b_n^2 + \texttt{sheath\_min\_bn}^2}$$

folded into `sheath_ramp` so residual and Jacobian can never disagree. Smooth, never changes sign, $\to1$ away from tangency, $\to0$ at it. It **removes** the term where it is unsolvable rather than distorting it. The discontinuous floor is deleted from `sheath_current`; validation now *recommends* `sheath_min_bn > 0` instead of refusing it.

**Suggested value:** `sheath_min_bn = 0.05`. `vpar_smoothing_coef(1) = 0.02` sets where the Chodura–Riemann function itself transitions, so gating a little above that removes the unsolvable band without touching the resolved target. The nodal path used 0.2, which also worked but gates off considerably more wall.

### Open question about the diagnostic cadence

`sheath_diag_report` is called from `construct_matrix`, which Newton calls **once per nonlinear iteration**, not once per step. So the records above may be Newton iterations rather than timesteps. `grep -c SHEATH logfile` against the step count settles it, and it matters: a monotone climb across *timesteps* is a physical runaway, whereas the same climb across *Newton iterations* within one step would mean the nonlinear solve is not converging — despite no abort message.

---

## 17. Root cause: the gate leaves 98.7% of the boundary with no condition on `u`

The gate-aware diagnostic settles it in two numbers:

```
gated-off area  98.7 %
ePhi/kTe max where the sheath is ACTIVE = 1.31 1.31 1.30 1.29 1.28 1.27 ...   (stable, declining)
ePhi/kTe max globally                   = 3.75 4.07 4.37 4.65 4.91 5.18 ...   (monotone climb)
```

**Where the sheath acts it is completely stable.** It settles at 1.31 and *declines* over thirteen records while the global max climbs without bound. The boundary condition works. The runaway is entirely on the 98.7% where the gate removed the term — and there, with `dirichlet%u = .false.`, `u` has **no boundary condition at all**. That is the same null space that killed `sheath_flux_sign = 0` at step 34, and $\partial u/\partial\ell$ is the flux-dragging velocity.

**It explains the whole ladder:**

| sheath on | freed boundary with no effective condition | crash |
|---|---|---|
| type 1 only | small | 205 |
| types 1, 4, 5, 9 | large | 30–34 |
| 1/4/5/9, `min_bn` 0.05 → 0.005 | smaller | 48 |

### 17.1 Two distinct causes

**(a) `sheath_min_bn = 0.05` was ~10× too large.** $b_n = \mathbf{B}\!\cdot\!\mathbf{n}/|B| = \sin\alpha$, and with $|B_{\rm pol}|/|B|\approx0.1$–$0.2$ and a 1–5° strike angle, $b_n\approx0.02$–$0.09$. At a 2° strike point $b_n\approx0.035$, giving gate weight $0.035^2/(0.035^2+0.05^2) = 0.33$ — **the gate was switching off the strike points themselves.** `sheath_min_bn = 0.005` gives 0.98 there. (32 → 48 steps confirms the direction.)

**(b) Most of a sheath-enabled boundary genuinely is tangential wall.** With `grid_to_wall`, types 1/4/5/9 cover most of the vessel perimeter; the strike zones are ~10–20 cm of a ~5 m poloidal perimeter. The sheath stiffness scales as $b_n^2$, so over that majority it is intrinsically too weak to hold `u` — correctly so, since no parallel flux reaches a tangential surface. **Gating the term off there is right; leaving nothing behind is not.**

> [!danger] The design error
> The **nodal** path has a fallback: `sh_relax = sh_relax * sh_wgt_bn` degenerates its row smoothly to $\delta u = 0$, i.e. the standard Dirichlet. I gated the **natural** path off without giving it anywhere to degenerate *to*.

### 17.2 The fallback: `sheath_wall_pen`

Where the gate removes the sheath, relax `u` toward the **local floating potential** with weight $(1-\text{gate})$ and the stiffness the sheath itself would have at $|b_n| = $ `sheath_wall_pen`:

$$\text{added to the } u \text{ row:}\qquad -\oint v\,R\,(1-w)\,\kappa\,\bigl(u - u_{\rm float}\bigr)\,d\ell\,\Delta t,
\qquad \kappa = c_{\rm sat}\rho\,c_s\,\frac{a_n}{2T_e}\,\texttt{sheath\_wall\_pen}^2$$

with $u_{\rm float} = 2(T_e\Lambda + v_w)/a_n$, i.e. $X = 0$.

**Why this and not a Dirichlet `u = 0`.** A surface carrying no net current *floats*, and the sheath solution itself tends to the floating potential as $j\to0$ — so the two regimes agree at the crossover and there is **no junction step**, which is what made every earlier mixed-boundary attempt fail. `u = 0` would reintroduce a $\Lambda T_e/e \approx 50$ V step across one element.

**Why it cannot be a weighted penalty row.** `boundary_conditions_add_one_entry` overwrites a single matrix entry with `zbig = 1e12`; even $10^{-3}\times$ that still dominates the assembled O(1) entries, so a weighted strong row is a strong row. The blend must live in the weak form, which is where the gate already is.

**Jacobian.** $\partial u_{\rm float}/\partial T_i$ and $\partial T_e$ are carried (including $\partial\Lambda/\partial T$). $\partial\kappa/\partial\text{state}$ is deliberately omitted: it multiplies $(u - u_{\rm float})$, which vanishes at this term's own fixed point, so Newton keeps its local rate — the same argument that fails for `sheath_stiff_max` (§14.3) but holds here.

### 17.3 Falsified along the way — do not retry

| Hypothesis | Killed by |
|---|---|
| Sheath **sign** wrong | `±1` behave identically (30 vs 33) |
| ψ **wall relaxation** | `wall_vel = 0` identical (33 vs 30) |
| **Solvability** at saturation | `sat_slope = 0.05` makes $f$ unbounded — zero effect (32 vs 33) |
| **Stiffness** degeneration | `stiff_max = 0` halved the growth rate, same crash step |
| Boundary-type **junctions** | adding type 3 changed nothing |
| Jacobian errors | Newton converges every step, no aborts |

**`dirichlet%zj = .false.` is catastrophic** — 1e32 within four steps. The dropped `zj` surface term weakly imposes $\nabla\psi\!\cdot\!\mathbf{n} = 0$, which fights `dirichlet%psi`. That route needs the `direction_perp(1)` Hermite trial column (§9.6) first.

---

## 18. `natural%zj` implemented properly — freeing the wall current

### 18.1 What was wrong

`mod_boundary_matrix_open` evaluates a field's normal derivative correctly:

```fortran
eq_t(mp,k,ms) += nodes(i)%values(in,j3,k) * element_size_ij * H1(i,j,ms) * HZ * element_size_perp
```

with `j3 = direction_perp(j)` — the DOF that genuinely carries $\partial/\partial n$. But the **trial** loop built

```fortran
psi_t = H1(k,l,ms) * element_size_kl * HZ(in,mp) * element_size_perp
```

and filed it under column `l2 = direction(l)` — the *value* DOF. So the true Jacobian entry was missing and a spurious one took its place: an effectively explicit $O(1/h)$ boundary term, which is what grew boundary structures over ~20 steps in the V2/V3 runs.

> [!important] The key realisation
> `psi_t` was numerically the **correct** normal-derivative coefficient all along. It was simply filed under the wrong column. Only the column assignment (and one sign) needed to change — no new basis function had to be derived.

### 18.2 The fix

$\nabla f\!\cdot\!\mathbf{n}$ splits exactly into the half carried by the tangential derivative and the half carried by the normal derivative:

$$\nabla f\!\cdot\!\mathbf{n} = \underbrace{\frac{y_t n_1 - x_t n_2}{J}}_{\texttt{gpn\_s}} f_s \;+\; \underbrace{\frac{-y_s n_1 + x_s n_2}{J}}_{\texttt{gpn\_t}} f_t$$

Verified against the original expression to $6\times10^{-13}$ over 20 000 random configurations. The `f_s` half stays in `amat` at column `direction(l)`; the `f_t` half goes into a new `amat_p`, assembled at `index_kl_p` with `l3 = direction_perp(l)`.

The pairing is the same one the field evaluation uses: `l=1` (value) ↔ normal first derivative, `l=2` (tangential derivative) ↔ mixed second derivative. With `n_order = 3`, `n_degrees = 4`, so `l3 = 4` is in range — the maximum `index_kl_p` lands exactly on the `ELM` bound. Note the **rows** still only reach value/tangential DOFs, which is correct: a surface integral cannot reach a test function that vanishes on the edge. The asymmetry is intended.

### 18.3 Two real bugs found on the way

| Bug | Effect |
|---|---|
| `l3 = direction_perp(j)` used **`j`**, the test-function index, instead of `l` | dead code until now, would have picked the wrong DOF |
| The trial block omitted the `vertex(1)*vertex(2) == 2` orientation flip that the field evaluation applies to `element_size_perp` | the normal-derivative Jacobian had the **wrong sign on half the boundary elements** |

The second is the same class of defect as the `sheath_wall_vel` orientation bug (§ earlier): a logical-coordinate derivative used without its orientation correction.

### 18.4 Why this frees the current

`dirichlet%psi = .true.` pins ψ's **value and tangential derivative** at the boundary — but **not its normal derivative.** The `zj` surface term couples `zj` at the wall to $\nabla\psi\!\cdot\!\mathbf{n}$, which is therefore free. That is the degree of freedom through which the wall current can respond to the sheath, and it is why `dirichlet%zj = .false.` without `natural%zj` blew up: the dropped surface term weakly imposed $\nabla\psi\!\cdot\!\mathbf{n} = 0$, fighting the ψ Dirichlet.

### 18.5 Validation changes

- `natural%zj` **requires** `dirichlet%zj = .false.` (a surface term only reaches rows the Dirichlet overwrites).
- `natural%u` **without** `natural%zj` now warns, quoting the measurement: `I_Ampere` constant to four digits while `I_wall` ran to zero and `u` diverged.
- `natural%zj` with `dirichlet%psi = .true.` prints a NOTE saying that is the intended combination, with the reason.
- The blanket refusal of `natural%w` / `natural%zj` is lifted; `natural%w` still needs `dirichlet%w = .false.`, and is not needed by the sheath.

### 18.6 Namelist

```fortran
 bcs(1)%natural%u     = .true.
 bcs(1)%dirichlet%u   = .false.
 bcs(1)%natural%zj    = .true.      ! NEW
 bcs(1)%dirichlet%zj  = .false.     ! NEW
 bcs(1)%dirichlet%w   = .true.
 bcs(1)%mach1         = .true.
 ! ... identical for 4, 5, 9; types 2 and 3 keep their defaults
```

**What to watch:** `I_Ampere` must stop being constant. If it moves toward `I_wall`, the j–V loop is closed for the first time and the boundary condition is doing what it was built to do.
