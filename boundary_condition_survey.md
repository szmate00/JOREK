# The `u`/`w` boundary-condition problem: a literature survey

## Why this document exists

JOREK's `u` (velocity stream function) and `w` (vorticity) satisfy `w = Δu`, with `w` advected by
a transport equation. Stripped of plasma physics this is the **streamfunction-vorticity (ψ-ω)
formulation**, equivalently the **mixed (Ciarlet-Raviart) formulation of a fourth-order problem**.
The pairing is identical in plate bending (deflection / bending moment) and Cahn-Hilliard
(concentration / chemical potential).

That family has one structural difficulty, and it is the oldest and most-studied boundary-condition
question in computational fluid dynamics:

> **The physical problem supplies two conditions on `u` and none on `w`. The mixed system needs one
> on each. So `w` at the wall is not data - it is an unknown.**

Every stable method below is a way of *computing* `w|_Γ`. JOREK's `dirichlet%w` *invents* it (freezes
it at its restart value). That is invisible while `u = 0` on the wall, because the terms that would
expose it vanish identically. It stops being invisible the moment `u` varies - which is exactly what
a floating-potential boundary condition does.

---

## Family A - local pointwise wall-vorticity formulas

**Thom (1933), Woods, Jensen, Briley.** Compute `w|_Γ` from the wall value of `u`, the first
interior value, and the mesh spacing.

- Thom is only first-order accurate locally but preserves second-order global accuracy.
- Briley is `O(h^3)` and preserves `L^2` stability.
- E & Liu: Thom's formula is *the correct implicit coupling*; treating `w|_Γ` explicitly or
  inconsistently with `u` **produces instability**. Existing wall-vorticity formulas are properly
  read as "the discrete counterpart of the Neumann condition for the stream function", not as
  Dirichlet data for vorticity.
- Review verdict: local formulas "had limited success"; vorticity boundary conditions "should be
  global in the sense that computing boundary vorticity values should involve nodes in both the
  interior and boundary".

**Relevance to JOREK.** Frozen `w` is worse than the explicit treatment E & Liu warn about - it is
*stale* data, never updated. The literature's minimum standard is a formula that reads the current
interior solution.

---

## Family B - influence / capacitance matrix

**Glowinski & Pironneau (1979); Kleiser & Schumann (spectral).** Treat `w|_Γ` as a vector of
unknowns. Split the solution into a particular part plus one **discrete harmonic** per boundary
node. Solve a small dense system - the influence matrix - for the boundary vorticities that make
the *second* condition on `u` hold.

- The influence matrix is **precomputed once**; each step costs only a small reduced solve.
- Spectral versions add a "tau correction" for the a priori lacking vorticity boundary conditions.
- Still actively developed (discrete-harmonics formulations, 2024).
- Cost note: the influence matrix is about twice the order one would expect from the analogous
  finite-difference or finite-element method.

**Relevance to JOREK.** This is the mature answer. JOREK has the ingredients (direct solver,
per-boundary-node DOFs, the definition row). The missing piece is the one-off harmonic
precomputation. Real project, known cost, known payoff.

---

## Family C - integral / global conditions

**Quartapelle & Valz-Gris (1981).** The velocity boundary conditions imply conditions of *integral*
type on the vorticity: the vorticity field must be orthogonal to the harmonic vector fields. This
determines a projection of `w` onto that manifold and **requires no pointwise boundary condition on
the vorticity at all**.

**Relevance to JOREK.** Conceptually the cleanest resolution of "there is no `w|_Γ` data": it says
there should not be. Harder to retrofit into an existing nodal-BC framework than Family B, but it
explains *why* every pointwise recipe feels arbitrary.

---

## Family D - weak (Nitsche) imposition of the Dirichlet data

**Modern; isogeometric analysis, Kirchhoff plates, mixed biharmonic.**

The diagnosis is stated directly for our configuration: with `u` enforced **strongly**, the mixed
formulation has *"no flexibility to properly represent σ at the boundary"* (σ = our `w`), giving
suboptimal convergence and reduced stability. Imposing `u = g` **weakly** restores optimal
estimates.

Three boundary terms on the `w`-definition row, test function τ:

```
(w,τ) + (∇u,∇τ)  -  ∮ τ ∂ₙu      -  ∮ (u-g) ∂ₙτ     +  (λ/h) ∮ (u-g) τ
                     consistency      adjoint-          penalty
                     [= natural%w]    consistency
```

- Skew-symmetric / non-symmetric variants are **parameter-free** for linear boundary conditions,
  at the cost of Galerkin orthogonality; reported to give *better* accuracy in derived quantities
  (bending moments, fluxes) - and here the derived quantity IS `w`.
- IGA uses C1 NURBS/Bezier elements on fourth-order problems: the same element technology as JOREK.

**Relevance to JOREK.** `natural%w` is the first of the three terms. `mod_boundary_matrix_open`
already has surface quadrature, the gradient split `∇f·n = gpn_s f_s + gpn_t f_t`, and the split
Jacobian (`amat`/`amat_p`) that routes tangential and normal-derivative DOFs to the right columns.
Missing: test-function derivatives `v_s`/`v_t` for `∂ₙτ`. Cheaper than Family B and keeps the `w`
equation untouched.

---

## Family E - energy-stable open boundaries / backflow stabilization

**Dong; Bazilevs; Esmaily Moghadam et al.** A *different* concern from all of the above: not
accuracy, but an energy argument.

At an open boundary the energy balance carries a term `∮ (1/2)|v|^2 (v·n)`. Where `v·n < 0`
(inflow) this **injects** energy. If nothing controls it, kinetic energy can grow without bound;
the literature reports exponential growth and "instant blow-up" at moderate Reynolds number.

Remedy: add a stabilization active **only** on the inflow portion,

```
-β ∮ (1/2)|v|^2 (v·n)₋ ,     (v·n)₋ = max(0, -v·n)
```

leaving outflow regions untouched.

**Relevance to JOREK.** On any boundary, `v_E·n = -R ∂u/∂τ` exactly (verified numerically to
1e-12). A constant `u` gives zero; a floating potential does not, and `∂u/∂τ` **changes sign along
the wall**, so parts of the boundary are inflow. Nothing in the model controls the energy influx
there. This is the only family that predicts *growth* rather than *error*, and it is untested.

---

## Family F - lifting / harmonic extension

Standard technique for nonhomogeneous essential data: construct a smooth extension `G` of the
boundary data into the domain, solve for `v = u - G` with homogeneous conditions. A *harmonic*
lift is the natural choice.

**Relevance to JOREK.** Avoids starting the nonzero trace as a one-cell boundary layer, and lets
the well-tested homogeneous machinery carry the rest. Attractive because the zero-target
configuration is known to be stable. Needs `G` in every place gradients of `u` appear.

---

## Family G - curvature consistency (Babuska-Sapondzhyan paradox)

For a simply supported plate, `w = 0` is correct only on a **straight** boundary with a
**constant** trace. On a curved boundary the true condition carries a curvature term
(`~ κ ∂ₙu`, plus the trace's second tangential derivative). Impose `w = 0` there and the scheme
converges - stably - **to the wrong solution**. Approximating a circle by polygons makes the error
grow with the number of sides.

**Relevance to JOREK.** The divertor wall is curved, and since `floating_u` the trace varies. Frozen
`w` is wrong on both counts, and was right before precisely because `u ≡ 0` killed both terms.
Caveat: this family explains *wrongness*, not blow-up.

---

## Family H - do not introduce the auxiliary variable

**M3D-C1** uses C1 triangular elements and applies Galerkin directly to the fourth-order operator
"without introduction of new auxiliary variables (such as the vorticity or the current density)".
Fewer auxiliary fields means fewer boundary conditions to invent. JOREK's elements are C1 too, so
the mixed formulation is a *choice* - but not one that can be revisited from a namelist.

---

## The plasma-side answer that already exists

**Loizu, Ricci, Halpern & Jolliet, Phys. Plasmas 19, 122307 (2012), "Boundary conditions for plasma
fluid models at the magnetic presheath entrance".** Analytically derives, from first principles and
in agreement with kinetic simulations, the boundary conditions at the magnetic presheath entrance
for `v_∥i`, `v_∥e`, `n`, `φ`, `T_e` **and for the vorticity `ω = ∇⊥²φ`**. GBS uses these.

This is the closest published answer to the exact question "what is the vorticity boundary condition
at a sheath, consistent with the potential boundary condition" - the `φ`/`ω` pair derived *together*
rather than one prescribed and the other invented. Automated extraction of the formulas from the PDF
failed; the paper should be read directly.

---

## Why this has not surfaced in the JOREK community

1. `u ≡ 0` on the wall hides the entire problem: every term that would expose the inconsistency
   vanishes identically for a constant trace on any boundary.
2. The nearest neighbour code (M3D-C1) has no auxiliary variables and therefore cannot have hit it.
3. The vocabulary does not overlap. Nothing in plasma-boundary language surfaces *Ciarlet-Raviart*,
   *Glowinski-Pironneau*, *influence matrix*, *tau correction*, *discrete harmonics*,
   *Babuska-Sapondzhyan*, or *backflow stabilization*.

---

## Ranked shortlist for this problem

| # | approach | family | cost | status |
|---|---|---|---|---|
| 1 | `natural%w` - definition supplies `w|_Γ` | D (partial) | done | queued |
| 2 | Read Loizu 2012 and adopt its `φ`/`ω` pair | plasma | reading | **not done** |
| 3 | Energy-influx diagnostic on the inflow part of the wall | E | small | **not done** |
| 4 | Full Nitsche: drop the `zbig` rows, add adjoint + penalty | D | ~50 lines | if 1 fails |
| 5 | Harmonic lifting of the trace | F | moderate | alternative to 4 |
| 6 | Glowinski-Pironneau influence matrix | B | project | last resort |

Items 2 and 3 are the cheapest untried things and come from the two families that were missing from
the earlier analysis.

## Sources

- E & Liu, *Vorticity boundary condition and related issues for finite difference schemes*
- *A review of vorticity conditions in the numerical solution of the ζ-ψ equations*, Comput. Fluids
- Glowinski & Pironneau, *Numerical methods for the first biharmonic equation and for the
  two-dimensional Stokes problem*, SIAM Review 21 (1979)
- Quartapelle & Valz-Gris, *Projection conditions on the vorticity in viscous incompressible flows*,
  IJNMF (1981)
- *Discrete harmonics for stream function-vorticity Stokes problem*, arXiv:2412.09996
- *An optimal error estimate for a mixed FEM for a biharmonic problem with clamped boundary
  conditions*, arXiv:2305.00407
- *Nitsche's method for Kirchhoff plates*, arXiv:2007.00403
- *Skew-symmetric Nitsche's formulation in isogeometric analysis*, CMAME
- *Energy-stable boundary conditions based on a quadratic form: outflow/open-boundary problems*,
  arXiv:1807.07056
- *On a new mixed formulation of Kirchhoff plates on curvilinear polygonal domains*, arXiv:1711.10260
- *Necessary and sufficient conditions for avoiding Babuska's paradox on simplicial meshes*,
  arXiv:2401.05897
- Jardin, *A high-order implicit finite element method*, JCP (2007) - M3D-C1
- Loizu, Ricci, Halpern & Jolliet, Phys. Plasmas 19, 122307 (2012)
