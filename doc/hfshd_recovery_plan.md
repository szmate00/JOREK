# Getting a pronounced, elongated HFSHD in JOREK

**Goal.** Reproduce the HFS high-density region observed in the experiment (AUG 38773),
pronounced and poloidally elongated, well enough for a quantitative experimental
comparison. This document is a recovery plan, not a code comparison. "JOREK cannot do
what SOLPS does" is not an acceptable output of this work; every branch below ends in
something to fix, measure or implement.

Companion documents: `hfshd_findings.tex` (what has been established and falsified),
`drift_bc_audit.tex` (boundary-condition audit), `jorek_vpar_bc_drifts.tex` (the vpar BC
derivation), `hfshd_test_plan.html` (the sheath-drive test sequence).

---

## 0. THE BINDING CONSTRAINT (added 2026-08-31, verified in code)

**Kinetic ionisation and charge exchange are switched off below 1 eV by a numerical guard, in
a divertor that sits at ~1 eV.**

    particles/mod_particle_evolution.f90:396
      limits = (n_e_raw .le. 1e14) .or. (T_e_raw*K_BOLTZ/EL_CHG .le. 1.d0) &
                                   .or. (T_i_raw*K_BOLTZ/EL_CHG .le. 1.d0)   !ADAS limits

`limits` gates ionisation (`:434`, `:558`), charge exchange (`:464`) and radiation (`:428`,
`:572`). The case's inner leg is flat ~1 eV with 0.66 eV at the strike point, so across the
whole leg there is **no ionisation source and no CX momentum sink**, while volume
recombination — on the fluid side, ungated — keeps draining it. The leg stays cold, dilute,
sonic and transparent to neutrals, and there is no path out of that state from inside the run.

This kills both halves of the loop at once: the parallel sink is never throttled (§1: nothing
accumulates, so nothing can be elongated) and Te never falls in a way that raises eta (so the
`int eta j_par dl` hill never gets its amplification). **The "ionisation front at the X-point"
is the 1 eV isotherm — it is the gate, not physics.**

The fluid channel does the same job correctly: `models/mod_atomic_coeff_deuterium.f90:77`
CLAMPS the ADAS argument to `max(Te, 0.2 eV)` and only zeroes ionisation/radiation below
**0.2 eV**, with an explicit comment that recombination stays allowed. The kinetic path uses a
floor 5x higher, applies it as an off-switch rather than a clamp, and additionally trips on
`T_i` — which protects no table lookup, since the CX rate is indexed on `T_e` (`:465`). The
code's own comment at `:464` says "CX uses adas as well. Te limit could be lower."

**Fix, highest value in this document:** clamp instead of gate, mirroring the fluid path; drop
`T_i` from the condition; put the floor behind a namelist parameter (`T_min_adas`, default 1.0
so the build stays bit-identical) so the A/B is one flag. Keep the `n_e <= 1e14` guard — that
one is a real table limit and never binds at 1e19-1e21.

**Decisive measurement first, zero compute:** on an existing checkpoint, export `T_e`/`T_i`
(not `Te` — that is `(Ti+Te)/2`, `mod_expression.f90:1999`) plus the kinetic coupling
projections `aux_rho0` (ionisation) and `aux_mom_par0` (parallel momentum, carrying CX
friction) over the divertor. Check whether both source terms are identically zero where
`min(Te,Ti) <= 1 eV`, and whether the 1 eV contour coincides with the peak of `aux_rho0`. If
they coincide, the gate is proven to be the constraint. If the sources are healthy and spread
along the leg, the gate is not binding and the next candidate is `use_sc = .false.`

---

## 1. The elongation is the most informative thing in the data

A poloidally *elongated* high-density region — extending from the inner strike point up
the inner leg toward and past the X-point — is not the same observable as an in-out
target density ratio, and it carries a much stronger constraint.

Plasma delivered into a divertor leg leaves it along the field at the parallel loss time

    tau_par ~ L_leg / v_par     ~  10 m / 3e4 m/s  ~  0.3 ms   (attached, sonic exit)

and is resupplied across the field by the drift at

    tau_fill ~ w_leg / v_ExB    ~  0.05 m / 100 m/s  ~  0.5 ms

These are comparable. An attached leg is a leaky bucket: whatever the drift delivers is
flushed to the target within a connection time, so the density rises only in a thin layer
near where it is delivered. **No drive of any magnitude produces an elongated structure in
an attached leg.** To pile plasma along the leg you need `tau_par >> tau_fill`, i.e. the
parallel sink throttled by one to two orders of magnitude. That is what detachment /
strong recycling does: momentum loss to neutrals drives `v_par -> 0` near the target and
the leg stops draining.

So the experimental elongation is direct evidence that the measured inner leg was
detached or close to it. The present JOREK case is not: flat ~1 eV divertor, 0.66 eV at
the strike point, sonic exit, and all four pump/valve polygons specified in bowtie order
so their effective areas are roughly halved (`hfshd_findings.tex`, next steps 5).

**Consequence for priorities.** The divertor regime is not an amplifier to be added later.
It is a precondition for the *shape* of the observable, not merely its amplitude. It is
also the explanation for the bias scan: 0.0005 per volt on the integrated in-out ratio is
what a leaky bucket looks like. Pushing harder on the drive in an attached leg was always
going to read as insensitive.

---

## 2. The chain, and the state of each link

| # | Link | Status |
|---|------|--------|
| 1 | Charge separation: curvature / diamagnetic drive in the u equation | **Present, correct, dominant.** `u_Eq__grad_p`, `mod_elt_matrix_fft.f90:1593`, measured 2.85e-5 = 57x `visco`. The 1.2e-7 figure was a misidentification: that is `u_Eq__diamag_term` (`:1603-1607`), the tauIC gyroviscous stress, not the charge-separation drive. |
| 2 | Sheath closure sets Phi at the targets | **Done.** Weak-trace Galerkin BC, `weak` 6.5e-4, closure 0.00% both targets, 3900 steps. |
| 3 | ExB transport across the PFR | Present, and larger than expected: 9% of the parallel flux crosses flux surfaces in the PFR and does not cancel. |
| 4 | Deposition/removal at the target knows about the drift | **Broken.** `var_rho`, `var_Ti/Te` boundary fluxes carry no ExB at all; `var_vpar`'s drift term is divided by Btot twice. |
| 5 | The receiving leg can accumulate (throttled parallel sink) | **Missing, and structurally prevented** — see §0 below. Attached, flat 1 eV, no detachment, pumps at half area. |
| 6 | Neutral source sustains and localises the structure | **Broken/untested.** Same polygon bug; recycling not established. |

Links 2 and 3 — the ones the campaign spent its effort on — are the two that work.

---

## 3. Stage 0: the cheap decisive checks (hours, no new runs)

**0a. Field direction.** HFSHD is drift-driven and its sign flips with the toroidal field
direction — that is the experimental control. Confirm the case was built with the shot's
field direction, and confirm the PFR ExB actually points outer-leg -> inner-leg in the
run. The measured pointwise `n_in/n_out = 0.53` (outer-dominant) against an integrated
1.17 is not by itself a sign error, but it has never been checked directly. If the drift
runs the wrong way, everything downstream is wasted effort. Cheapest possible test with
the largest possible consequence.

**0b. Replace the observable.** Target traces cannot measure divertor content — the common
integration window is only +-5.5 cm and excludes the inner leg's density peak at +5.6 cm.
Two changes:
  - volume metric per leg (`rectangle`), not target traces;
  - an **elongation metric**: n_e along a `pol_line` following the inner leg from the strike
    point through the X-point into the HFS SOL. Report the poloidal extent over which
    n_e exceeds a fixed fraction of its peak, plus the peak location. That is the quantity
    the experiment shows as "pronounced and elongated"; nothing currently measures it.

**0c. Synthetic diagnostic, not a code-internal number.** For the experimental comparison,
build the line integral the shot actually measured (interferometer chords / whatever 2D
emission diagnostic is available for 38773) and compare that, rather than comparing a
JOREK volume integral to a quoted experimental "density ratio". This removes an entire
class of arguments about what was integrated over what.

---

## 4. Stage 1: make a leg that can accumulate (highest payoff)

In order, all no-code:

1. **Fix the two PUMP polygons only** (CORRECTED 2026-08-31). Valves and pumps use
   DIFFERENT vertex conventions: valves are a bilinear map
   (`particles/mod_particle_puffing.f90:427`), which wants raster order and is therefore
   CORRECT as supplied - do not touch them. Pumps go through even-odd `inside_polygon`
   (`particles/mod_particle_wall_interaction.f90:1309`), where raster order IS a bowtie:
   P1/P2 keep the R-extreme triangles and drop the Z-extreme ones, losing half the area.
   Reorder cyclically, e.g. P1 R = 1.29,1.31,1.31,1.29 / Z = -1.14,-1.14,-1.12,-1.12.
   Estimated current throughput ~2e20/s against a 3e22/s puff = 0.6%, so the D inventory
   ratchets and the run has no steady state.
2. **`T_min_sheath = 1.d-6`**, then the puff/recycling scan, to drive the inner leg into
   high recycling and toward rollover. Acceptance: inner-target `v_par` falling well below
   c_s over an extended region; a pressure drop along the inner leg; an ionisation front
   that detaches from the plate.
3. **`use_sc = .false.`** (or step `D_perp_sc_num` 10 -> 3 -> 1). Shock-capturing diffusion
   acts exactly where gradients are steep, and an HFSHD *is* a steep gradient. This is the
   untested suspect in `hfshd_findings.tex` §"Numerical suppression" and it has never been
   run alone.
4. ~~Viscosity scan in the u equation.~~ **Dropped** — the premise was a misidentified term.
   The charge-separation drive is `u_Eq__grad_p` at 2.85e-5, 57x larger than `visco`, and it
   is balanced against `u_Eq__JxB` to 0.03%, which is the correct `div j_dia = -div_par j_par`
   closure. `visco_num_term` on the legs is 1.3% of the dominant term in its own equation — a
   co-symptom of steep leg gradients, not a short circuit. Do not spend a run on this.

**Decision point after Stage 1.** Repeat the drift A/B (or the bias scan) in the
high-recycling state. If the sensitivity per volt jumps by more than ~10x, the drive was
never the limiting factor and the remaining work is Stages 2-3. If it does not move,
escalate to §7.

---

## 5. Stage 2: boundary conditions that know about the drift

These pay off *in the detached regime specifically*: as `v_par -> 0` at the target the ExB
term stops being a correction and becomes a leading contribution to the boundary flux. The
drift correction goes as `v_E.n / b_n` and `|b_n|` averages 0.0080 (0.46 degrees), already
reaching 150% of c_s locally in the attached case. Doing these before Stage 1 would test
them in the regime where they matter least.

1. `Mach1BC`: remove the double 1/Btot conversion on the drift term and add SOLPS's
   `+-2 c_s |b_x|` clip (the measured would-reverse fraction is 0.89%, so the clip is
   needed). Behind `bohm_drift_compatible` for a one-build A/B. Already implemented on
   this branch — validate it.
2. ExB in the `var_rho` and `var_Ti/Te` target fluxes (SOLPS BCCON=14, BCENE/BCENI=15).
   The manual is explicit that the five drift-compatible BCs must be used in conjunction;
   we currently have one and a half of them.
3. Fluid and kinetic recycling **together** — `mod_boundary_matrix_open.f90:661` and
   `mod_particle_wall_interaction:1770` mirror each other and changing one alone breaks
   particle conservation across the interface. Watch `c_angle`: it is a grazing-incidence
   floor and the drift correction is largest at grazing incidence, so it may already be
   standing in for the drift.
4. Collision background velocity: pass the full `vvector` (ExB + parallel) instead of
   `P(1)*B/t_norm`. One line; currently drops 100% of the perpendicular background flow.

---

## 6. Stage 3: match the experimental profiles

With a leg that accumulates and BCs that transport, close the loop on the plasma state:
iterate `D_perp_file` / `ZK_perp_file` with `util/flux_match_diffusivities.py` against the
measured upstream and target profiles for the shot. An HFSHD sitting on the wrong
background n,T is not a match even if its shape is right.

---

## 7. Escalation: if it is still short after Stages 1-3

Only reachable after the above, and the output is an implementation task, not a null
result.

1. **Measure the potential drop the mechanism needs.** `Delta Phi` along the field line
   from the X-point to each strike point, against `integral eta*j_par dl` over the same
   path. SOLPS gets ~100 V. Tooling exists (`util/postproc_hill.in`,
   `util/analyse_hill.py`) and has never been run. This is the measurement that separates
   "mechanism present but subcritical" from "mechanism absent".
2. **Term budget of the current balance vs SOLPS, with matched n and T.** Comparing
   `div j_dia` between a flat 1 eV JOREK leg and a structured SOLPS one measures the case
   setup, not the code. Constrain the JOREK divertor state first (Stage 3), then compare
   term by term. Candidates for a genuine shortfall, in the order they are worth checking:
   the anomalous perpendicular current (SOLPS carries one; JOREK's counterpart is
   numerical viscosity), viscous currents, and the parts of the curvature drive dropped by
   the reduced-MHD `B_phi = F0/R` ansatz.
3. **Implement what the budget says is missing**, behind a flag, and A/B it. If a current
   term is genuinely absent from the reduced-MHD vorticity equation and the budget shows
   it carries the drive, that is the thesis contribution — a term added to JOREK that lets
   it reproduce an experimental divertor asymmetry it previously could not.
4. **Prescribed-Phi run as a bounding case, used as a diagnostic and never as the final
   answer.** `sheath_V_wall_asym` with the SOLPS or probe-measured target potential
   profile answers "given the correct Phi, does the rest of the model respond?" It
   localises the failure to drive vs response. It is a step in the debugging, not a result.

---

## 8. Standing rules from this campaign

- Integrated metrics only. A pointwise strike-point ratio moved 63% while the extensive
  metric moved 1.5%; the probe sits where plasma is being pushed to.
- Match A/B runs on `t_now`, never on step number — runs diverge in dt by 2.6x.
- `psi_N` does not separate SOL from PFR and is reflected below the X-point; use signed
  arc length, and `A_3` or Z-restricted `Psi_N` for gradient orientation.
- `si-units` / `jorek-units` are bare postproc words; `set units 0` is silently ignored.
- Keep `keep_current_prof_confined` on. The uncut term is a 45% artifact in the psi
  equation and produced a fake HFSHD once already.
