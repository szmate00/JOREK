# Prescribed floating-potential BC: campaign record

Branch `floating-u`, based on `keep-current-prof-cutoff-clean`.
Status as of 2026-09-05: the boundary condition is **exactly imposed** and **destabilises
the run**. Λ=0 reaches timeout; Λ=3 crashes at step 607. The failure is a transport
problem at the target, not an implementation error.

Companion documents: `doc/floating_u_derivation.tex` (normalisation and linearisation),
`doc/sheath_sign_review.md` (the Φ sign), `doc/HANDOFF_phi_3Te_bc_plan.md` (the plan this
implements).

---

## 1. What is implemented

The floating condition `V_p − V_wall = Λ k_B T_e / e` written as an **affine nodal
constraint** on the `u` trace,

```
u  =  C_T · T_e  +  C_V · V_wall
```

with (`models/model600/mod_floating_u.f90`)

```
a_n = + 2 e F0 sqrt(mu0 rho0) / m_i        rho0 = central_density*1e20 * central_mass * u
C_T = 2 Λ / a_n                            (Λ halved when .not. with_TiTe, since Te = T/2)
C_V = sqrt(mu0 rho0) / F0                  V_p[V] = F0 u / sqrt(mu0 rho0)
```

The relation is **exactly affine**, so it is imposed exactly in the single linearised
solve per step — there is no Newton loop to converge (`mod_jorek_timestepping.f90:377/390`).
The row is written in `mod_boundary_conditions.f90` inside the existing Dirichlet trace
loop, on exactly the DOFs that loop already selects (the value and the pure **tangential**
derivatives), as the ordinary `zbig` diagonal plus one `T_e` cross-column and the residual
RHS.

### Sign

`a_n` is **positive**, i.e. `Φ = +F0·u`. `mod_sheath_bc.f90:207` on the sheath branch and
`feature/floating-potential-bc` both use the opposite sign; `mod_boundary_conditions.f90:876`
uses this one. Settled in `doc/sheath_sign_review.md` from three non-comment anchors:
JOREK's `(e_R, e_Z, e_phi)` is right-handed. A first derivation that gave `Φ = −F0 u`
assumed a right-handed `(R, φ, Z)` and was wrong.

### Deliberate omissions

- **No positivity map.** An earlier version applied `corr_neg_temp1(T_e)` inside the row.
  Removed in `1780e841c` for two independent reasons: the value row's Jacobian used `−C_T`
  instead of `−C_T·f'(T_e)` (4.92× too large at 0.3 eV), and `corr_neg` was being applied to
  **Fourier coefficients** — a zero non-axisymmetric `T_e` coefficient maps to 0.84 eV and
  manufactures a spurious ~2.5 V harmonic. Harmless at `n_tor=1`, fatal in 3D. Deleting it
  fixes both and restores exact affineness. Negative temperatures are the global positivity
  scheme's job.
- **No amplitude ramp, no local Λ(T_i/T_e).** `feature/floating-potential-bc` had both.

### Commits

```
8e33551dc  prescribed floating-potential BC
e2448edb7  drop the natural%u check (component absent on this base)
1780e841c  drop the positivity map
245b41dba  per-type boundary diagnostic
09b58d899  stale comment
cdd385126  Pe column reported rank 0's local sentinel
0d4079e3f  signed vE.n, split into outflow and inflow
```

---

## 2. The diagnostic

`floating_u_diag = .true.` prints, once per step and per boundary type, at the
**axisymmetric harmonic only**:

| column | meaning |
|---|---|
| `\|u-uf\|[V]` | trace residual `\|u − C_T T_e − C_V V_wall\|`, in volts — the acceptance gate |
| `\|vE.n\|[m/s]` | normal ExB speed the BC imposes |
| `Pe` | `\|v_E·n\|·h_perp / D_perp(1)` |
| `Pe at (R,Z)` | location of the worst Pe |
| `min rho`, `min T[eV]` | boundary-trace minima |
| `vE.n out`, `vE.n IN` | signed split (added 0d4079e3f) |

**`v_E·n = R · ∂u/∂ℓ` is exact, not an estimate.** For `v_E = R ∇u × e_phi` in the
right-handed basis, with edge tangent `t = (R_b,Z_b)/dl` and normal `n = (Z_b,−R_b)/dl`,

```
v_E·n = R (∂_R u · R_b + ∂_Z u · Z_b) / dl = R · u_b / dl
```

So the normal ExB flow is driven by the **tangential** derivative of `u` — precisely what
`u = C_T T_e` manufactures wherever `T_e` varies along the wall.

### Two traps this diagnostic is built to avoid

- **Only MAX and MIN are accumulated.** Both are idempotent, so halo elements visited by
  several ranks cannot inflate them. Area integrals in the earlier weak-sheath campaign
  were inflated ~11 % by exactly this.
- **Orientation.** The tangent-derived normal `(Z_b,−R_b)` is *not* outward by
  construction — element vertex ordering and the `element_size` sign flips both enter it.
  Magnitudes are immune (the flips cancel in `|u_b|/|dl|`), but the signed columns are not,
  so they are oriented against `normal_direction`, which the loop already builds as
  node-minus-opposite-vertex. An uncorrected version of this produced a bogus "inverted =
  1448" count earlier in the campaign.

### One bug found in the diagnostic itself

`cdd385126`: `fd_pe_max` was carried into the `MPI_MAXLOC` through `fd_loc` so the worst
location could be broadcast, but the print used `fd_pe_max` itself. On any rank owning no
floating boundary node — rank 0 in practice — that is still the `−1` initialisation
sentinel, so every run before the fix printed `Pe = -1.00E+00`. The `(R,Z)` column was
unaffected because it is broadcast from the owning rank.

---

## 3. The four runs

Floating types 1, 3, 4, 5, 9. Same restart, same equilibrium, one variable changed each.

| run | Λ | `dirichlet%w` | outcome | failure site |
|---|---|---|---|---|
| **A** | +3 | `.true.` (frozen) | **crash 607** | OUTER target, R ≈ 1.60, Z ≈ −1.11 |
| **B** | +3 | `.false.` (free) | crash 510 | INNER target, R ≈ 1.26, Z ≈ −1.00 |
| **C** | −3 | `.true.` | crash 427 | INNER target, R ≈ 1.26 |
| **D** | 0 | `.true.` | **TIMEOUT, no crash** | — (control) |

### The row is exact

Residual at step 142, every type, both Λ=0 and Λ=3:

```
type 1  4.96e-12 V      type 3  4.64e-13 V      type 4  7.85e-12 V
type 5  8.68e-12 V      type 9  2.75e-12 V
```

The implementation is not in question. Reviewer priority 10 (solver / row failure) is dead.

### The null is perfect

`|v_E·n|` at step ~140:

```
Λ = 0    ~1e-12 m/s      (noise; worst-point location wanders randomly)
Λ = 3    ~3e3  m/s       (worst-point location FIXED across every step)
```

Fifteen orders of magnitude. Nothing else in the code produces this flow.

### The control implicates the BC

Λ=0 ran to **timeout with no crash**, at the *same* 0.646 eV type-1 divertor floor as A.
So the 607-step failure is not the underlying cold-divertor state — it is the imposed drift.

---

## 4. What the crash looks like (run A)

**No degradation trend.** Three outputs before the end, A was as healthy as the Λ=0
control (type 4 min ρ 5.7e-2 against D's 4.9e-2; residual 1.9e-8 V). Then four steps to
blow-up. This is a trigger, not an accumulation.

**The collapse oscillates, it does not drain.** Type 9 `min rho`:

```
4.086e-02  →  1.513e-04  →  1.859e-03  →  -1.538e-01
```

A factor of 270 down in one step, partial recovery, then negative. Monotone
over-evacuation cannot do that; a dispersive under/overshoot can.

**Thermal collapse is a consequence.** Type 1 `min T` goes 0.648 → 0.0358 eV one output
*after* ρ goes negative — the same causality as in the weak-sheath campaign.

**Types 4 and 9 go negative first in all three failing runs.** Those are exactly the two
types that failed *alone* on the weak route (crash at 8 and 4 steps) while types 1 and 5
ran. Two unrelated BC formulations, the same two types leading. That is independent
evidence for the grid defect, not for either BC.

---

## 5. Hypotheses falsified by measurement

Including three of my own.

1. **"The drift grows exponentially and runs away."** Read as +1.5 %/step off seven points
   and projected a doubling every 47 steps. **Wrong** — it saturated at ~3.2e4 m/s and run
   A was healthy at step 582 with the drift slightly *decreasing*.

2. **"Frozen `w` is the defect, free it."** `bcs(:)%dirichlet%w = .true.` writes only the
   `zbig` diagonal and no RHS; JOREK solves for increments, so boundary `w` is **frozen at
   its restart value**, not zeroed — a genuine inconsistency once `u` varies, since `Δu ≠ 0`.
   But freeing it crashed **earlier** (510 vs 607) and degraded the `u` residual by five
   orders (type 4: 1.8e-9 → 1.2e-4), three outputs *before* ρ went negative. Frozen `w` is
   helping. Same pattern as releasing `dirichlet%u` on the corners. The inconsistency is
   real; this is not the remedy.

3. **"The historical 3500-step run survived because of its sign."**
   `feature/floating-potential-bc` used `flt_a_n = −2 e F0 sqrt(mu0 rho0)/m_i`, opposite to
   HEAD. But Λ=−3 crashed at **427**, the earliest of the three. Something else carried
   that branch (different equilibrium, its C2 amplitude ramp, local Λ(T_i/T_e), or
   `corr_neg`). Weak corroboration that +3 is the better-behaved sign.

4. **"A cold divertor kills it."** A ran at 0.65 eV on type 1 and reached 607; C died at
   4.2 eV.

5. **"Mixed corner topology imposes incompatible conditions"** (reviewer P3) and
   **"solver failure"** (P10). Both dead by the Λ=0 control: the same corners are
   constrained on the same DOFs with the same structure at Λ=0, and that run times out.

6. **"ExB boundary flux is double-counted"** (reviewer P9). Dead — see §6; nothing was
   integrated by parts, so there is no surface term to double-count.

7. **"ρ is unconstrained on the floating types."** Wrong: `natural%rho = .true.` on
   exactly types 1, 4, 5, 9 (`preset_parameters.f90:461-465`), so ρ has the proper weak
   surface-flux treatment. Type 3 is the anomaly — `dirichlet%rho = .true.`,
   `natural%rho = .false.`, so its ρ is *frozen*. It stayed positive during A's collapse
   because it is pinned, which proves nothing.

---

## 6. Leading hypothesis: unspecified ExB **inflow** on the density

The ExB convection of ρ is **advective, not conservative**
(`mod_elt_matrix_fft.f90:1655`):

```fortran
+ v * BigR**2 * ( r0_s * u0_t - r0_t * u0_s) * tstep          ! R^2 [rho,u]
+ v * 2.d0 * BigR * r0 * u0_y * xjac * tstep                  ! compressibility
```

Both multiply `v`, the **test function**, never `v_x`/`v_y`. Two consequences:

1. **No ExB surface term is missing.** Nothing was integrated by parts, so there is no
   boundary term to be absent. Reviewer P9 is dead — for this reason, not the one first
   given.

2. **But an advective operator needs data wherever characteristics enter.** The natural ρ
   BC (`mod_boundary_matrix_open.f90:343`) carries only the **parallel** flux:

```fortran
rhs_ij(var_rho) = + v * density_reflection * r0 * vpar0 * ps0_s * normal_sign3 * tstep &
                  - v * r0 * cs0 * BigR * dl * c_angle * tstep    ! min-angle particle flux
```

   Nothing ExB, nothing `u`-dependent.

With `u = 0` this was harmless: `v_E·n = 0`, so the boundary is characteristically neutral
and no inflow condition is needed. With `u = C_T T_e`, `v_E·n = R ∂ℓu` **reverses sign
across the T_e peak** at the strike point, so part of the target has ExB characteristics
entering the domain with **no incoming density specified at all**.

This covers every observation:

| observation | explained |
|---|---|
| oscillatory, not monotone | yes — underdetermined inflow, not over-evacuation |
| localised at the strike point | yes — that is where `v_E·n` changes sign |
| Λ=0 immune | yes — `v_E·n = 0` |
| **Λ=−3 also fails, worse** | yes — flipping the sign moves inflow to the other side of the peak |
| residual stays exact until the end | yes — the `u` row is not what is failing |

### The tension underneath it

`v_E·n = R ∂ℓ u` and `u = C_T T_e`, so **wherever T_e varies along the target the BC
necessarily drives ExB flow through a solid wall.** That is intrinsic to prescribing
Φ = Λ k_B T_e/e pointwise on a surface with temperature structure — not a numerical
artefact. Whether that flow is physical, and what the correct boundary treatment of it is,
is the question this campaign has actually reached.

---

## 7. Open items and next steps

**No grid changes.** The BC must work out of the box on any reasonable setup; refinement
would be a workaround, and the missing inflow condition is a defect on *any* grid — an
advective operator with undefined incoming characteristics is underdetermined, not
under-resolved. Refining changes how the ambiguity manifests, not whether it exists.

1. **Rebuild and read the `vE.n IN` column** (`0d4079e3f`). If inflow is negligible on the
   failing types the hypothesis is dead and the cell-Péclet story returns; if it is
   comparable to `vE.n out`, it is live. This needs no grid change.
2. **Free exclusion test: halve `tstep`.** Neither candidate mechanism is fixed by a
   smaller timestep — cell Péclet and inflow well-posedness are both spatial. If halving it
   materially extends the run, both are wrong and it is a CFL problem.
3. **Free namelist test:** `bcs(4)%floating_u = .false.`, `bcs(9)%floating_u = .false.`
   (floating on 1, 3, 5). Discriminates the 4/9 grid defect from the drift in general.
   Leaves type 4's 0.539 m² uncovered, ~1 % of wetted area, and puts a potential jump at
   the junction — a diagnostic, not a candidate configuration.
4. **Work out the inflow closure**: what incoming ρ and T should be where `v_E·n < 0` at a
   material surface, and how it is expressed in this weak form. This is the piece that has
   to be right for the BC to work anywhere.

### Not to be done

- **No `D_perp_sc_num` scan, no fitted thresholds.** The BC must work with defaults.
- **Do not free `w`** — measured worse on both step count and residual.
- **Do not pin ρ** on the floating types — that removes real physics and only reproduces
  type 3's artificial stability.

### Known-weak points in the current implementation

- The diagnostic samples **boundary nodes only**. The interior density hole near the strike
  point is not sampled; `min rho` must not be read as exonerating anything there.
- The `u` row omits `∂(u_target)/∂ψ` — but the target has no ψ dependence, so unlike the
  weak-sheath row this is not an omission, merely worth stating.
- Boundary `w` remains frozen and inconsistent with `Δu ≠ 0`. Measured harmless-to-helpful,
  but unresolved in principle. A Nitsche or influence-matrix treatment of the `u`–`w` pair
  is the proper answer if it ever becomes load-bearing.
