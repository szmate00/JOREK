# Questions for the authors of "Mechanisms of XPR formation from SOLPS-ITER modeling"

Poletaeva, Kaveeva, Senichenkov, Rozhansky, Mulyukov, Shtyrkhunov,
Phys. Plasmas 33, 012507 (2026), doi:10.1063/5.0289248

Framed to be directly useful for pushing JOREK (model600, reduced MHD + kinetic neutrals and
impurities) toward XPR simulation. Each question notes why it matters on the JOREK side.

---

## 1. Are the gradB terms in the ENERGY equations essential to the trigger?

Eqs. (3)-(4) carry four magnetic-drift contributions: `p_a div(V_a^gradB)` and
`V_a^gradB . grad(p_a)`, for both species. JOREK model600 has **none of them in the temperature
equations** - its only perpendicular energy advection is E×B. (The density and vorticity
equations do carry diamagnetic terms.)

You show that the far-core balance is neoclassical, i.e. parallel heat flux against gradB
convection (Fig. 6a), while at the LCFS the X-point cooling is dominated by `Q_conv` and radial
conduction. **Could you run the dynamic case with the gradB terms removed from the ENERGY
equations only** - keeping the drifts in the particle and momentum equations, and hence keeping
HFSHD - and report whether XPR still forms?

*Why:* this is the single most decisive test of whether model600's known gap is fatal for XPR or
merely changes the background poloidal Ti modulation.

## 2. How much did the 5/2 -> 3/2 convective-flux change move the threshold?

You replaced `div(5/2 n V_turb T)` with `div(3/2 n V_turb T)` for consistency with dropping the
turbulent-diffusion RHS terms, and note it matters for X-point cooling dynamics but not for
attached or semi-detached regimes. **By how much did that change the seeding rate or the timing
of XPR onset?**

*Why:* JOREK has no separate "anomalous convective" energy flux - perpendicular energy transport
is a diffusive `ZK_perp grad(T)` plus E×B advection of pressure. If the 3/2-vs-5/2 distinction
shifts the threshold significantly, there may be no clean JOREK counterpart and we need to know
what the equivalent choice is.

## 3. Is HFSHD necessary, or just the dominant channel in AUG? What amplitude?

You write that a strong HFSHD region provides the cold plasma for convective cooling, and that
"drifts determine the HFSHD formation." **What in-out density asymmetry do you actually get in
the benchmark** - `n_HFS/n_LFS` at the targets and near the X-point - **and is there a threshold
below which the convective channel is too weak to trigger XPR?**

*Why:* we are working hard to obtain a physically-driven HFSHD in JOREK. A target number turns
that from an open-ended effort into a pass/fail criterion. If XPR needs (say) a factor 3
asymmetry, and we can only produce 1.2, that tells us where we stand.

## 4. Amplitude and origin of the potential hill

Your chain is: density hill -> p_e hill -> potential hill via `en grad_par(phi) ~ grad_par(p_e)`
-> E×B vortex. **What is the hill amplitude in volts and in units of Te, and what E×B velocity
does the vortex reach?** And is the amplitude set by the parallel electron force balance alone,
or does the target sheath/current closure materially change it?

*Why:* we are implementing a self-consistent sheath current boundary condition on the potential
and need to know (a) what amplitude to expect, and (b) whether an incorrect wall potential can
suppress the vortex entirely.

## 5. What sheath/current boundary condition do these runs use?

**Specifically:** the sheath heat transmission coefficients, the potential BC (floating vs
current-carrying), whether the wall is treated as a global equipotential, and **whether there is
net current to the targets in the benchmark** (i.e. is the wall globally floating, or grounded
with a return path?).

*Why:* we find a ~2.5 kA thermoelectric loop between the targets and a net wall current that is
not zero, and we are unsure which of those is physical versus an artefact of our gauge choice.
Knowing SOLPS's convention would settle it.

## 6. Sensitivity to the parallel heat conductivities

You emphasise that Ti drops first because `chi_i,par << chi_e,par`. **How sensitive is XPR onset
to the ratio and to the absolute values?** If both were reduced by, say, a factor 10 (as is
sometimes done for numerical reasons), would the trigger survive with only a stretched timescale,
or would the mechanism change qualitatively?

*Why:* reducing `ZK_par` is a common numerical expedient in JOREK. If it breaks the Ti-first
ordering, it breaks the mechanism, and we need to know that before relying on such runs.

## 7. Timescale separation and the dynamic scheme

Transition takes ~30 ms at `U_N = 1e21` atoms/s, 2-2.5x above the stationary collapse limit.
**What timestep did the dynamic runs use, and what fraction of the 30 ms is the "initial cooling
stage" before radiation takes over?** Also: does onset correlate better with accumulated impurity
content than with seeding rate?

*Why:* JOREK's timestep is tied to MHD timescales. Knowing the separation between the cooling
time and the fastest dynamics tells us directly whether a JOREK XPR run is computationally
feasible, and whether we could shortcut by prescribing an impurity inventory rather than a rate.

## 8. How large is the non-coronal correction?

You stress that impurity radiation is essentially non-coronal during XPR formation, so a coronal
`L_z(Te)` is a rough estimate. **How large is the discrepancy at the critical Te ~ 40-60 eV -
factor 2, more?** And would a coronal-equilibrium radiation model reproduce the onset at all, or
only with a rescaled impurity concentration?

*Why:* our impurity radiation comes either from ADAS tables or from kinetic impurity particles
with their own charge-state evolution. Knowing the size of the non-coronal effect tells us
whether the kinetic treatment is essential or a refinement.

## 9. The reciprocal of the baffle test

Your baffle test (neutrals excluded, drifts on) still forms XPR, in contrast to Ref. 17 (neutrals
excluded, drifts off). **Have you run the reciprocal case in the dynamic setup - drifts off,
neutrals allowed?** And with drifts off, does HFSHD disappear entirely or merely weaken?

*Why:* your two tests differ in two variables at once, so "drifts are needed" is inferred rather
than isolated. The reciprocal test would close that, and it directly informs how much drift
physics a JOREK XPR attempt actually needs.

## 10. What single quantity would you use for a JOREK/SOLPS benchmark?

JOREK is time-dependent MHD with self-consistent perpendicular transport (no imposed
`D_perp`/`chi_perp`) but reduced drift physics; SOLPS has the full drift set with prescribed
anomalous coefficients. **At the XPR-formation stage, which single quantity would you consider
most discriminating for a cross-code comparison?**

Candidates we would propose: the X-point `Ti/Te` ratio versus time; the poloidal `Te` profile at
the LCFS; or the temperature at which `Q_conv` and `Q_rad` cross (~40 eV in your Fig. 28b).

*Why:* we would like any JOREK XPR result to be comparable to yours rather than merely
qualitatively similar, and it is much cheaper to agree the metric first.

---

## Bonus, if there is appetite for collaboration

- Ref. 31 reports turbulence modification around the XPR. Since JOREK can generate filamentary
  transport self-consistently rather than through prescribed coefficients, **what would you most
  want to learn from such a simulation** that SOLPS cannot provide?
- You leave the feedback-controlled stable XPR to future work. **What sensor and actuator pairing
  do you expect to work** (X-point Te on seeding rate?), and what is the characteristic response
  time the controller has to beat?
