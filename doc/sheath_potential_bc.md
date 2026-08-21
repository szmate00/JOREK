---
tags: [jorek, model600, boundary-conditions, sheath, sol, hfshd]
branch: bc-tests
commit: bb5d7fc31
status: implemented, not yet run
---

# Sheath $j$–$V$ boundary condition for the electric potential (JOREK model600)

Implementation of the Bohm/Stangeby sheath current–voltage characteristic as a boundary
condition on the electrostatic potential, following J. Artola, *"Sheath boundary conditions
for the electric potential in JOREK"* (revised 2026).

> [!abstract] One-line summary
> The sheath forces a relation between the current reaching a material surface and the
> potential drop in front of it. Imposing it makes the plasma potential at the target follow
> $\Phi \approx \Lambda T_e/e$, which gives a radial electric field, which drives poloidal
> $E\times B$ flows in the SOL — the leading candidate mechanism for in–out asymmetries such
> as HFSHD.

---

## 1. The physics

### 1.1 Why a sheath forms

Electrons are far more mobile than ions ($v_{th,e}/c_s \sim \sqrt{m_i/m_e} \sim 60$ for
deuterium). A surface immersed in plasma is therefore hit by electrons much faster than by
ions and charges **negative**. The resulting electric field, concentrated in a few Debye
lengths at the surface, repels electrons and accelerates ions until the two fluxes balance.

That thin layer is the **sheath**. It is far below MHD resolution, so it enters a fluid code
only as a boundary condition.

### 1.2 The current–voltage characteristic

Let $\Phi = V_\text{sheath entrance} - V_\text{wall}$. Only electrons with enough energy to
climb the potential hill reach the wall, so their flux is Boltzmann-suppressed by
$e^{-e\Phi/k T_e}$, while the ion flux is fixed at the Bohm value. The net current density is

$$
j \;=\; j_\text{sat}\left(1 - e^{-\left(\frac{e\Phi}{k T_e} - \Lambda\right)}\right)
\;\equiv\; j_\text{sat}\, f
$$

with the **sheath factor**

$$
\Lambda \;=\; \ln\sqrt{\frac{m_i}{2\pi m_e}} \;\approx\; 3 \quad\text{(deuterium)}
$$

> [!important] $\Phi$ is referenced to the **wall**, not to the floating potential
> This is the single most consequential detail, and the one the first version of Artola's note
> got wrong (missing $e^{+\Lambda}$).
>
> | condition | $X \equiv e\Phi/kT_e - \Lambda$ | $f$ | meaning |
> |---|---|---|---|
> | $\Phi = \Lambda T_e/e$ | $0$ | $0$ | **floating** — zero net current |
> | $\Phi \to \infty$ | $\to\infty$ | $\to 1$ | ion saturation |
> | $\Phi = 0$ | $-\Lambda$ | $1-e^{\Lambda} \approx -19$ | **electron saturation** |
>
> So $\Phi = 0$ is *not* zero current. It is the largest current in the problem. Getting this
> backwards puts the boundary condition at completely the wrong operating point.

The floating potential is positive because the field must point *at* the wall to hold the
electrons back:

| $T_e$ | 3 eV | 5 eV | 10 eV | 25 eV | 37 eV |
|---|---|---|---|---|---|
| $\Phi_f = \Lambda T_e/e$ | 9 V | 15 V | 30 V | 75 V | 111 V |

### 1.3 Ion saturation current

The maximum current the sheath can pass on the ion side. In SI,
$\mathbf J_\text{sat} = e n v_\parallel \mathbf B$ with $v_\parallel$ the JOREK parallel-velocity
variable (velocity divided by $|B|$), and matching the *normal* component to the poloidal
current $\mathbf J_\text{pol} = -j\mathbf B_\text{pol}/F_0$:

$$
-\frac{j_\text{sat}}{F_0}B_n = e n v_\parallel B_n
\quad\Longrightarrow\quad
j_\text{sat} = -e\,n\,v_\parallel F_0
$$

> [!note] $B_n$ cancels
> $j_\text{sat}$ carries **no** field-incidence factor. This is correct — but it also means the
> relation degenerates ($0/0$) where the field is exactly tangent to the wall. See
> [[#6 Limitations|Limitations]].

### 1.4 What this buys you for SOL physics

At low current the characteristic pins $\Phi \approx \Lambda T_e/e$ at the target. Since $T_e$
varies **radially** across the target, so does $\Phi$:

$$
E_r = -\partial_r \Phi \approx -\frac{\Lambda}{e}\,\partial_r T_e
$$

Parallel dynamics equalise potential along field lines quickly, so the target sets the potential
of the whole SOL flux tube. The resulting $\mathbf E\times\mathbf B$ drift is **poloidal** and
reverses with $B_\phi$ — the experimental fingerprint of drift-driven asymmetries.

| $|\nabla T_e|$ | $|E| = \Lambda|\nabla T_e|$ | $v_{E\times B} = E/B$ |
|---|---|---|---|
| 500 eV/m (10 eV over 2 cm) | 1500 V/m | ~700 m/s |
| 1250 eV/m (25 eV over 2 cm) | 3750 V/m | ~1700 m/s |

That is comparable to or faster than cross-field diffusive transport at typical $D_\perp$, so it
is not a small correction.

> [!warning] With the default `u = 0` boundary condition this drive is **absent, not approximated**
> $u=0 \Rightarrow \Phi = 0 \Rightarrow E_r = 0$. Any SOL $E_r$ in a standard run comes only from
> interior dynamics.

---

## 2. JOREK units

### 2.1 Normalisations

With $\rho_0 = n_0 m_c m_u$ the central mass density and $n_0 = $ `central_density` $\times 10^{20}$:

$$
\Phi = F_0\, u_\text{SI},\qquad
u_\text{SI} = \frac{u}{\sqrt{\mu_0\rho_0}},\qquad
k T_\text{SI} = \frac{T}{\mu_0 n_0},\qquad
j_\text{JOREK} = \mu_0\, j_\text{SI}
$$

Hence

$$
\frac{e\Phi}{kT_e} = \frac{a_n u}{2T_e},
\qquad
\boxed{\;a_n \equiv \frac{2 e F_0 \sqrt{\mu_0\rho_0}}{m_c m_u}\;}
$$

(the factor 2 is $T_e = T/2$ in the single-temperature model; the code uses $T_e$ directly so it
generalises to `with_TiTe`), and

$$
j_\text{sat} = c_\text{sat}\,\rho\, v_\parallel,
\qquad
c_\text{sat} \equiv -e F_0 n_0\sqrt{\mu_0/\rho_0} \;=\; -\tfrac{1}{2}a_n
$$

> [!tip] Useful numbers for AUG ($F_0 = 2.972$, `central_density` $=1.011$, D)
> - $a_n = 186.29$, $c_\text{sat} = -93.15$
> - **1 unit of $u$ = $4.576\times10^{6}$ V** — use this to convert output
> - $j_\text{sat} \approx 0.06$–$0.12$ MA/m² at 37 eV, $n\sim(0.5$–$1)\times10^{19}$
> - toroidal current: $j_\phi \,[\text{A/m}^2] = -\,zj/(\mu_0 R)$

### 2.2 With a wall potential

$\Phi = F_0 u_\text{SI} - V_\text{wall}$ gives

$$
X = \frac{\frac{1}{2}a_n u - e V_\text{wall}\mu_0 n_0}{T_e} - \Lambda
$$

`sheath_V_wall` is in volts; 0 means a grounded wall.

---

## 3. Discretisation

### 3.1 Invert, don't iterate

The naive route is Newton on $N = j - j_\text{sat}f(u) = 0$, which is Artola's eq. 17. **Don't.**
Near saturation $\partial f/\partial u \to 0$: the relation carries no information about $u$ at
all, so any form solved for $u$ is singular exactly where a divertor target usually sits.

Instead invert exactly. With $r = j/j_\text{sat}$:

$$
X = -\ln(1-r), \qquad
u_\text{target} = \frac{2\left(T_e\left(X+\Lambda\right) + v_w\right)}{a_n}
$$

### 3.2 Clip $r$, never the exponent

Clipping $r$ to $[f_\text{min}, f_\text{max}]$ bounds $X$ to
$[\texttt{exp\_min}, \texttt{exp\_max}]$ **by construction**, so $u$ can never run away, and

$$
\frac{\partial r}{\partial(\text{state})} = 0 \quad\text{outside the range}
$$

makes the Jacobian *honest*: where the plasma demands more current than the sheath can pass,
$u$ is simply the cap and has no $\rho$ or $j$ dependence at all.

> [!danger] Clipping the exponent instead is a trap
> It keeps the row invertible but leaves every coefficient at its exponentially amplified value
> $\xi \sim e^{X_\text{max}}$. The potential then becomes slaved to grid-scale noise in $T$,
> $\rho$ and $j$ — which produced element-scale checkerboards on the target plates.

`sheath_u_exp_min` should equal $-\Lambda$: that is where $\Phi = 0$, i.e. where electron
saturation physically is, giving $|j| \le (e^{\Lambda}-1)j_\text{sat} \approx 19\,j_\text{sat}$
— close to the deuterium value $\tfrac12\sqrt{m_i/\pi m_e}\approx 17$.

### 3.3 Coefficients

The row is $\delta u + \sum_k C_k \,\delta x_k = u_\text{target} - u_0$ with
$C_k = -\partial u_\text{target}/\partial x_k$ and $\xi = 2T_e/[a_n(1-r)]$:

$$
C_\rho = \frac{\xi r}{\rho},\quad
C_j = -\frac{\xi}{j_\text{sat}},\quad
C_{T_i} = \frac{\xi r}{2T},\quad
C_{T_e} = \frac{\xi r}{2T} - \frac{2(X+\Lambda)}{a_n}
$$

all multiplied by $\partial r$ (0 when clipped) except the explicit $T_e$ term.
Verified against an exact numerical Jacobian in the unsaturated, ion-clip and electron branches.

### 3.4 Derivative rows

The constraint differentiated along the boundary is applied to the tangential-derivative DOFs,
**with the $\partial j/\partial\ell$ term removed**.

> [!bug] Why $\partial j/\partial\ell$ must go
> $zj$ at a boundary node is a *second* derivative of $\psi$, so its derivative DOF is a
> *third* derivative — discontinuous across $C^1$ cubic Bézier elements. It is grid noise, not
> a gradient, and it dominates that row by 4–10×. Since $\partial u/\partial\ell$ is the
> $E\times B$ flow through the wall, that noise gets injected straight back into the plasma.
>
> The value row keeps the full coupling, so the characteristic still holds **exactly at the
> nodes**; only interpolation between them is affected.

Derivative rows are written in nodal-DOF units (no `element_size` scaling — it is common to all
variables and cancels), which keeps the diagonal at `zbig`.

---

## 4. Where it goes in the matrix

### 4.1 Rows at a boundary node

| row | carries |
|---|---|
| `psi` | $\delta\psi = 0$ (standard Dirichlet) |
| `u` | **the sheath characteristic** |
| `w` | $\delta w = 0$ (standard) |
| `zj` | $\delta j = 0$ (standard — frozen) |

Only the vorticity evolution equation is discarded, which a standard run discards anyway.

### 4.2 The surface-term constraint

> [!danger] The governing constraint on any boundary-condition work in model600
> The `psi`, `u`, `w` and `zj` weak forms are **all assembled without their surface terms**.
> Only `rho, T, Ti, Te, rhon, vpar` have natural-BC support in `mod_boundary_matrix_open`.
>
> Leaving *any* of those four rows to its own equation at a boundary node makes that node
> absorb the missing surface integral.

For example the current definition is assembled as
$\int(\nabla v\cdot\nabla\psi + v\,j)/R = 0$, whereas the weak form of $j = \Delta^*\psi$ is

$$
\int \frac{\nabla v\cdot\nabla\psi + v j}{R} \;-\; \oint v\,\frac{\nabla\psi\cdot\hat n}{R} = 0
$$

A free boundary $zj$ therefore picks up $\approx \nabla\psi\cdot\hat n / h$ — tens of MA/m² in a
one-element layer. The same holds for $w$ with $\nabla u\cdot\hat n/h$.

**Adding the term in `mod_boundary_matrix_open` does not fix it**: the surface and volume terms
must cancel to a relative $\sim zj\,h/|\nabla\psi| \approx 3\times10^{-3}$, and the two routines
use different representations of the normal derivative. It would have to be assembled inside
`mod_elt_matrix_fft`.

### 4.3 Consequence

The current stays Dirichlet, so there is no self-consistent sheath current and no thermoelectric
current. `initialise_parameters` **stops** if `dirichlet%zj` or `dirichlet%w` is `.false.` on a
type with `sheath_u`, and explains why.

---

## 5. Using it

### 5.1 Parameters

| name | default | meaning |
|---|---|---|
| `bcs(N)%sheath_u` | `.false.` | enable per boundary type |
| `sheath_Lambda` | `3.d0` | $\ln\sqrt{m_i/2\pi m_e}$ |
| `sheath_V_wall` | `0.d0` | wall potential [V] |
| `sheath_u_exp_max` | `2.d0` | ion-side clip; caps $\Phi$ at $(\Lambda+2)T_e/e$ |
| `sheath_u_exp_min` | `-3.d0` | electron-side clip; **keep at $-\Lambda$** |
| `sheath_u_relax` | `1.d0` | under-relaxation; fixed point unchanged |

### 5.2 Input

```fortran
 ! enable on EVERY boundary type bounding the plasma
 bcs(1)%sheath_u = .true.
 bcs(2)%sheath_u = .true.
 bcs(3)%sheath_u = .true.
 bcs(4)%sheath_u = .true.
 bcs(5)%sheath_u = .true.
 bcs(9)%sheath_u = .true.

 sheath_Lambda   =  3.d0
 sheath_V_wall   =  0.d0
 ! dirichlet%zj and dirichlet%w stay .true. (defaults) - enforced
```

> [!warning] Enable it on **all** boundary types, not just the targets
> $u$ is continuous along the boundary. A type with the BC next to one without it puts a step of
> order $\Lambda T_e/e$ across a single element, and since $v\cdot\hat n = R\,\partial_\ell u$
> that is a large artificial $E\times B$ jet at the junction. This produced "current structures
> at the edge of boundary type 1" in testing.
>
> Find your types with `Number of nodes per boundary type` in the grid log. Wall-aligned grids
> use 1,2,3,4,5,9 (simple) or 11,12,15,19 (final numbering) — **not necessarily type 1**.

### 5.3 Checks, in order

1. **BC is live** — plot $\Phi = 4.576\times10^6 \cdot u$ along the target against
   $\Lambda\,T_e[\text{eV}]$. They should coincide where the current is small.
2. **Drive exists** — $E_r$ and poloidal $v_{E\times B}$ in the SOL/PFR. Expect ~1 km/s.
3. **Current bounded** — $|j|/j_\text{sat} \le 0.865$ on the ion side.
4. **Transport response** — density redistribution. Slowest; drift transit ~0.3 ms, so allow
   several ms.
5. **Reverse $B_\phi$** — a drift-driven asymmetry must flip. This is the decisive test.

---

## 6. Limitations

- **No self-consistent sheath current.** $\delta j = 0$ at the boundary, so no thermoelectric
  current. Only the potential route is active. Requires the surface-term fix.
- **Near-tangential incidence.** $j_\text{sat}$ has no $B_n$ factor (correct, since $B_n$
  cancels), but the relation degenerates as $B_n\to 0$. `direction` $=\text{sign}(B_\text{pol}\cdot\hat n)$
  flips there. Consider a $|b\cdot\hat n|$ gate if trouble appears at grazing angles.
- **Corner nodes.** Visited twice with different `direction`/$B_\text{tot}$; the two rows sum.
  Pre-existing weakness of Mach1 too.
- **`T_min`.** Floors $T_e$ in the BC. At `T_min = 1.d-4` that is 4.9 eV, which inflates both
  $\Phi$ and $j_\text{sat}$ in a cold divertor. Consider lowering toward `t_min_neg`.
- **$E\times B$ through the wall.** A varying $u$ means $v\cdot\hat n \neq 0$. Tolerable only
  because both `psi` and `zj` rows are Dirichlet, so the inconsistency is discarded rather than
  compensated by a current.

---

## 7. Failure log

What was tried before arriving here. Kept because each entry is a trap that looks reasonable.

| attempt | outcome | cause |
|---|---|---|
| free `zj`, no other change | crash @6 | missing surface term → spurious boundary current |
| free `zj` + PR-803 row swap | crash @208 | `zj` becomes the multiplier enforcing $\partial_t\psi=0$, scaling as $\{\psi,u\}/\eta$ ($5\times10^7$ at $\eta=2\times10^{-8}$) |
| free `zj` + surface term added in `mod_boundary_matrix_open` | crash @105 | volume/surface cancellation to $3\times10^{-3}$ across two routines with different normal-derivative representations |
| free `u`, `w = 0` in the `w` row | crash @270 | `w` row carries the *definition*; `w` froze while `u` drifted |
| free `u`, `w = 0` moved to the `u` row | crash @50 | `w` row's own missing surface term, $\nabla u\cdot\hat n/h$ |
| clip on the exponent | checkerboard on plates | coefficients left at $\xi\sim e^{X_\text{max}}$ |
| $\partial j/\partial\ell$ in derivative rows | element-scale noise | third derivative of $\psi$ |
| missing $\Lambda$ | no in–out asymmetry possible | $f$ has one sign along the whole boundary |
| BC on type 1 only | structures at the type-1 edge | step in $u$ at the junction |

> [!note] Related upstream regression
> Commit `dcc0dfd36` (Artola, 2023, ITER PR 803 / IMAS-4585) added the `i_var1`/`i_var2` swap
> letting `zj` or `w` be free. It is an ancestor of `develop` but the code is **absent** — a
> later merge took the competing refactor of that block. Still missing on GitHub. Worth
> reporting.

---

## 8. References

- J. Artola, *Sheath boundary conditions for the electric potential in JOREK*, note, rev. 2026
  — eqs. 1–18; the corrected version carries the $+\Lambda$ in the exponent.
- P. Stangeby, *The Plasma Boundary of Magnetic Fusion Devices*, eq. 2.68.
- Artola et al., halo-current paper — uses a **linear penalty** form
  $\Phi = \Phi_w + \alpha(\mathbf J_\parallel - \mathbf J_\text{sat})\cdot\mathbf n$ instead,
  explicitly noting the Langmuir formulation is *"numerically very challenging"*. Worth
  revisiting if this proves fragile.

## 9. Code map

| file | what |
|---|---|
| `models/model600/mod_boundary_conditions.f90` | the BC itself, in the shared Mach1 geometry block |
| `models/phys_module.f90` | parameters, `bcs(:)%sheath_u` |
| `models/preset_parameters.f90` | defaults |
| `models/model600/initialise_parameters.f90` | namelist + consistency checks |
| `communication/broadcast_phys.f90` | MPI pack/unpack |
| `models/mod_log_params.f90` | logging |
