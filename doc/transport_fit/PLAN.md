# Automated transport-profile fitting for JOREK model600 — plan

Status: implemented 2026-10-14 on branch `transport-fit` (util/transport_fit.py, util/transport_fit.sh,
util/transport_fit_toy.py, doc/transport_fit/README.md); validated on the toy, not yet run with JOREK.
See section 9 for what the toy validation changed.

Goal: a new person sets up a model600 run (neutrals, impurities, sheath BCs, drifts
all switched on), drops a Langmuir-probe file next to it, starts one script, and comes
back to tabulated D_perp / ZK_perp profile files that hold the step-0 upstream profiles
steady and reproduce the measured target profiles as well as the SOL transport allows,
with a convergence log and figures that say how good the match is and where it is not.

## 1. What the research says (one paragraph each)

**SOLPS practice (Canik 2011, Wilcox 2026).** The de-facto standard is a fixed-point
rule applied cell by cell on the outer-midplane slice: D_new = -Gamma_code / grad n_exp,
chi_new = -(q_code - 5/2 Gamma T) / (n grad T_exp). Converges in fewer than ten
iterations from a flat guess. Documented failure modes: runaway target T_e when the
target density gets too low, flat or hollow profile regions (a diffusive fit needs a
finite negative gradient), gradient/flux degeneracy (profiles matched with 5x different
throughput), and far-SOL chi driven to extreme values when energy leaks through the
radial grid edge, which then gets "baked in" through the code flux. Their hygiene:
update D first until the density is close, then chi; cap the boundary-cell chi early on.

**Gradient-based calibration (Carli 2022/2026, adjoint in SOLPS-ITER).** BFGS on a
relative L2 misfit of midplane and target data, drifts kept on throughout. 99% of the
cost reduction within ~10 iterations. Large parameter jumps crash drift cases, so they
ramp drifts back in after each change. Six knots per coefficient already over-fitted;
per-cell coefficients oscillated with a spike at the separatrix. Upstream-only data was
ill-conditioned; adding target data was what made the solution consistent at both ends.

**Bayesian / surrogate (Fu 2026 UEDGE, Luo 2025 EMC3).** Hundreds of forward runs for
five parameters. Robust to noise and crashes (failed runs are dropped or continued from
the nearest converged point), but far too expensive per evaluation for a kinetic JOREK
run. Useful ideas: log-space parameters, log-uniform bounds, error-bar floors
(sigma = 5% + absolute floor), and treating a crash as a data point.

**Core flux matching (TGYRO, PORTALS).** Relaxed Newton on log gradients with a step
cap and automatic halving when the residual does not decrease. The "simple relaxation"
mode caps the change at 5% per step and is slow but never diverges.

**JOREK papers.** Nobody fits D/chi automatically. Cathey 2023, Hoelzl 2021, Liang 2026
all use ad-hoc profiles "designed to maintain the initial profiles" and Liang admits
only one of the two goals (stationary XPR, matched midplane profiles) was reached.

**SOLPS-ITER manual (local copy).** No fitting recipe. Transport profiles are
(position, value) tables with monotone cubic interpolation, first and last points well
outside the grid, PFR and divertor regions excluded or rescaled by flags, and a legacy
undocumented Levenberg-Marquardt mode on constant coefficients. Separatrix density is
held by puff feedback while transport is fitted.

**JOREK plumbing (measured in this tree).** Profile files are headerless two-column
psi_N tables, re-read on every restart, 50-500 points, strictly increasing psi_N,
positive values, no clamping outside the table (linear extrapolation can go negative).
The private flux region reads the SOL table at the mirrored psi_N, so it cannot be set
separately. `midplane outer` exports a cut from the axis to the grid edge including the
SOL. `boundary_quantities` exports R, Z, Psi_N and the sheath quantities along the wall.
Which heat-conductivity file is used depends on the compiled `with_TiTe`. With
`use_zkperp_times_density` the file holds a chi-like quantity. Particles need
`part_restart.h5` carried between segments.

## 2. Design decisions

1. **Update rule: gradient-ratio fixed point in log space, per channel, per knot.**
   ln D_{k+1} = ln D_k + r * ln(|du/dpsi|_sim / |du/dpsi|_ref), clamped per
   iteration to [1/f, f]. Same for chi_e, chi_i; with `use_zkperp_times_density` the
   ratio is density-weighted, (rho |dT/dpsi|)_sim / (rho |dT/dpsi|)_ref. This is
   Canik's rule when the flux is diffusive, needs no source or flux export, and is a
   fixed point whatever the drift flux is. Default r = 0.7 without drifts, 0.5 with.
2. **Adaptive relaxation from iteration 3 on.** The loop measures its own sensitivity
   per knot, s = d ln|du/dpsi| / d ln D from the last two iterations (expected -1 for
   pure diffusion, smaller magnitude where drifts carry the flux). The step becomes
   -ln(ratio)/s, clamped to the same [1/f, f]. Knots with s > 0 (wrong sign: changing D
   made the mismatch worse) are frozen and flagged. This replaces hand tuning of r and
   is the TGYRO/Broyden idea in its diagonal form.
3. **Coarse knots, fine file.** Updates live on ~25 knots in psi_N (dense 0.85-1.10,
   coarse elsewhere); the written file is 200 points from a monotone (PCHIP)
   interpolation. Ratios are smoothed over neighbouring knots before applying. This is
   the regularisation Carli found necessary and it keeps the separatrix spike out.
4. **Upstream reference = step 0 of the run itself.** Exported with the same `midplane
   outer` command as every iteration, so no mapping of experimental data is needed
   upstream. Fit range psi_N in [PSI_MIN, psi_N(grid edge)], default PSI_MIN = 0.7:
   the core re-equilibrates too slowly to be fitted per segment and barely moves
   anyway. Core knots are held.
5. **Downstream reference = probe file (R, Z, n_e, T_e, optional errors).** Each probe
   is projected onto the wall polyline from the boundary export, giving a plate, a
   distance along the plate from the strike point (the Psi_N minimum on that plate),
   and a psi_N from the simulation's own equilibrium. The SOL side of each plate is
   compared in psi_N, the PFR side in distance only (check, not fit, because of the
   mirrored-table coupling).
6. **Downstream enters the fit through the same rule.** For psi_N > 1 the gradient
   ratio is a weighted blend of the upstream SOL mismatch and the target mismatch,
   both taken against psi_N: w_up * ln(ratio_midplane) + w_tgt * ln(ratio_target),
   defaults w_up = w_tgt = 0.5, user-settable (w_tgt = 0 gives a pure upstream fit).
   The target n_e gradient in psi_N is the SOL density width seen at the plate, so this
   is the width fit without a separate secant. Absolute target levels (peak n_e, T_e)
   are reported, not fitted: they are set by fuelling and pumping, which stay fixed
   during the fit (Horvath/Wilcox degeneracy; Rule 3 below).
7. **Staging.** Iterations 1-2 update D only. chi_e, chi_i join once the density
   residual has dropped below 2x its final tolerance or after iteration 3, whichever is
   first. Boundary-knot (grid edge) chi is capped at CHI_EDGE_MAX for the first three
   iterations (Wilcox's leak guard).
8. **Steadiness gate.** Every segment exports its first and last step. If the profile
   drift over the segment exceeds half of the mismatch to the reference in any fitted
   channel, the state is not steady: no update, segment length doubled (up to
   NSTEP_MAX), rerun. Global checks logged alongside: particle inventory and thermal
   energy from the JOREK macroscopic output, and the wall particle and heat fluxes from
   the boundary export.
9. **Crash handling.** A crashed segment restores the previous coefficient files,
   halves the clamp factor, and restarts from the last good restart file. Three
   consecutive crashes stop the loop with a report naming the knots that changed most.
   With drifts on, the clamp default is 1.5 instead of 2.
10. **Stopping.** Converged when every fitted channel's weighted RMS relative mismatch
    is below TOL (default 5%) for two consecutive iterations, or the coefficient change
    is below 3% everywhere, or MAX_ITER (default 12) is reached. Each iteration writes
    a CSV line (mismatch per channel, max factor applied, steadiness metric, wall
    fluxes) and a PNG (reference vs simulated profiles, current D and chi in m2/s,
    target comparison per plate, flagged knots).
11. **Fuelling and pumping are frozen during the fit.** The loop never touches the
    puff rate or pump weight. Target density level is a reported diagnostic; if it is
    off by more than its error bar at convergence the report says so and points at
    fuelling, not transport. The XPR puff controller (xpr-controller branch) is the
    place for that, as a separate step.
12. **Everything resumable.** A `state.json` in the run dir records the iteration
    counter, last good restart, clamp factor and sensitivities, so a killed driver
    resumes where it stopped. Coefficient files get `.iterNN` backups.
13. **Units.** The files stay in JOREK units. The tool converts to m2/s for all
    reporting using central_density and central_mass from the namelist
    (D_SI = D / sqrt(mu0 rho0); chi_SI = ZK / ((gamma-1) rho_norm sqrt(mu0 rho0)),
    or ZK / ((gamma-1) sqrt(mu0 rho0)) with `use_zkperp_times_density`).

## 3. Deliverables

| File | Role |
|---|---|
| `util/transport_fit.py` | Python core, numpy only (matplotlib optional for figures). Subcommands: `check`, `reference`, `update`, `report`, `plot`. |
| `util/transport_fit.sh` | Bash driver, settings as variables at the top, flat copy into `$BASE`. Owns the loop: segment, postproc, update, gate, resume. |
| `util/transport_fit_toy.py` | 1D steady-state solver with a hidden D/chi, a pinch term, noise and a flat region; closed-loop test of the update rule. |
| `doc/transport_fit/README.md` | The guide for the new person: inputs, one command, how to read the figures, what the flags mean. |
| `doc/transport_fit/PLAN.md` | This file. |

Existing pieces reused: `util/compare_wall_runs.py` (boundary-export parser, plate
split, strike-point ordering), `util/setinput.sh` (namelist edits),
`util/continue_run.sh` pattern (restart file handling).

## 4. Interfaces

### 4.1 Driver settings (top of `transport_fit.sh`)

```
BASE=/tokp/work/<user>/<case>        # run dir; everything flat in here
NAMELIST=input
JOREK=$BASE/jorek_model600
POSTPROC=$BASE/jorek2_postproc
SUBMIT="sbatch --wait job.sh"        # or "mpirun -np 8 $JOREK < $NAMELIST"
PROBE_FILE=$BASE/probes.dat          # R Z ne Te [ne_err Te_err], '#' comments
REF_STEP=0                           # upstream reference = this step's export
NSTEP_SEG=200   NSTEP_MAX=800        # steps per segment, cap after doubling
MAX_ITER=12     TOL=0.05
RELAX=0.7       CLAMP=2.0            # 0.5 / 1.5 recommended with drifts
PSI_MIN=0.7     W_TARGET=0.5
CHANNELS="rho T_e T_i"               # which of D, ZK_e, ZK_i to fit
CHI_EDGE_MAX=10                      # m2/s, first 3 iterations
```

### 4.2 Probe file

```
# R[m]  Z[m]  ne[1e19 m-3]  Te[eV]  [ne_err]  [Te_err]
1.602  -1.183  2.1  18.5  0.3  4.0
...
```
Inner and outer plates are found automatically from the nearest wall segment. Missing
error columns default to 15% plus floors (0.1e19, 2 eV). T_e from probes is known to
be high in cold divertors; the default target weighting is n_e 1.0, T_e 0.5.

### 4.3 Postproc commands (written by the driver)

```
namelist input
si-units
set linepoints 400
expressions R Psi_N rho ne T_e T_i dprof zkeprof zkiprof
mark_coords 1
for step A to B by (B-A) do
  midplane outer
done
expressions R Z Psi_N bnd_dl Bnorm ne T_e Jn_sat heatF_total partF_total bnd_type_lo bnd_type_hi
boundary_quantities
```
(The `bnd_type_*` expressions are the uncommitted working-tree additions; the tool
works without them, the plate masks just fall back to the Bnorm threshold.)

### 4.4 `transport_fit.py check` (pre-flight, run first, run always)

Prints a checklist and exits non-zero on a hard failure:
- namelist readable; `restart`, `nstep_n`, `nout`, `nout_particles`, `restart_particles`
  consistent with segment running; `D_perp_file` and the right ZK file(s) for the
  compiled `with_TiTe` present (bootstraps them from the tanh coefficients if absent);
  `use_zkperp_times_density` read and echoed; table psi_N range covers the grid's
  psi_N range (extrapolation would go negative);
- `jorek_restart.h5` and `part_restart.h5` present; `jorek2_postproc` runs on
  REF_STEP and the midplane cut reaches past psi_N = 1;
- probe file parses, every probe lands within 2 cm of the wall polyline, both plates
  found, strike points found;
- normalisation factors echoed in m2/s so the user sees what the current D, chi are.

### 4.5 Loop (one iteration)

1. `setinput.sh` sets `restart=.t.`, `nstep_n(1)=NSTEP`; copy the newest `jorek*.h5`
   to `jorek_restart.h5` (particles: `part_restart.h5` is already in place).
2. `$SUBMIT`; on non-zero exit go to crash handling.
3. Postproc the first and last step of the segment (midplane + boundary).
4. `transport_fit.py update`: steadiness gate; blend upstream and target ratios;
   sensitivity update; clamp, smooth, interpolate; write files; append CSV; PNG.
5. Exit code tells the driver: continue / converged / not-steady (rerun longer) /
   stop (3 crashes or flagged everywhere).

## 5. Validation plan

1. **Toy closed loop** (`transport_fit_toy.py`): hidden D with an edge barrier, a
   pinch of 30% of the diffusive flux, 3% noise, one flat region. Pass criteria:
   profile error below 1% within 10 iterations, flat region untouched, pinch region
   converges with the adaptive sensitivity, wrong-sign knots flagged when the pinch is
   raised to 150%.
2. **Driver dry run**: a `SUBMIT="true"` mode that replays an existing run's h5 files
   as if they were segments, so the loop mechanics (restart copying, postproc, resume
   from `state.json`, crash branch) are tested on a laptop without JOREK.
3. **Digital twin on JOREK**: take an existing attached AUG case, run it to a steady
   state with a known tabulated D/chi, export its midplane and target profiles, build a
   fake probe file from the target export. Start the loop from constant D/chi and
   check it recovers the known profiles. First without drifts and kinetic particles,
   then with both on.
4. **Real case**: the HFSHD AUG setup with the real probe file. Document what matched,
   what did not, and which knots were flagged.

## 6. Known limits and candidate code changes (not part of v1; ask first)

- **PFR shares the SOL table.** The target PFR side cannot be fitted. A small
  model600 change (a `D_perp_pfr_factor` / `ZK_perp_pfr_factor` multiplier applied
  where `get_psi_n` reports a mirrored value) would mirror SOLPS's
  `b2tqna_pfr_rescale`. Flagged as a v2 code change.
- **No in-run ramp of coefficients.** JOREK reads the files at start; a large change
  is applied as a jump. v1 mitigates with the clamp; v2 could ramp over the first N
  steps of a segment (`D_perp_file_prev` + ramp length), which is what Carli needed for
  drift cases.
- **Impurity transport** follows D (`D_perp_imp` is not file-driven in model600).
  Neutral diffusivity is physics, not fitted.
- **Non-uniqueness** is real: with fuelling frozen the fit is unique up to the
  smoothing, but a different puff would give a different D. The report states the
  frozen fuelling and the resulting wall flux so the degeneracy direction is visible.
- **Compute.** One iteration is one restart segment plus two postproc calls. For the
  AUG kinetic case that is of order one to three hours per iteration on the cluster,
  so a full fit is one to two days unattended.

## 7. Order of work

1. Spec freeze: this document, plus a 30-minute check on an existing run that the
   step-0 midplane export reproduces rho_file / T_file (it should, within interpolation).
2. `transport_fit.py` core and the toy test.
3. `transport_fit.sh` and the dry-run mode.
4. Digital-twin test on JOREK, attached, no drifts.
5. Digital-twin with drifts and kinetic particles; tune the defaults (RELAX, CLAMP,
   NSTEP_SEG) from what that shows.
6. README for the new person; commit everything.

## 8. Sources

Canik et al., JNM 415 (2011) S409, doi:10.1016/j.jnucmat.2010.11.084.
Wilcox et al., NF 66 (2026) 016023, doi:10.1088/1741-4326/ae1bb9.
Carli et al., CPP 62 (2022) e202100184; Carli et al., arXiv:2609.27435 (2026).
Fu et al., PPCF 68 (2026) 015010, doi:10.1088/1361-6587/ae2b0b.
Horvath et al., PPCF 65 (2023) 044003.  Boeyaert et al., arXiv:2311.14607.
Kaveeva et al., NF 60 (2020) 046019.  Zito et al., PPCF 63 (2021) 075003.
Candy et al., PoP 16 (2009) 060704 (TGYRO).  Rodriguez-Fernandez et al., NF 64 (2024) (PORTALS).
Luo et al., NF 65 (2025) 096016.  Bold et al., arXiv:2407.12072.
Cathey et al., NF 2023 (arXiv:2301.09066).  Liang et al., PoP 33 (2026) 052508.
SOLPS-ITER manual 3.0.9, sections 3.5.3, 3.5.5, 5.x (local copy).

## 9. What the toy validation changed (2026-10-14)

- **Steadiness gate** compares the mid-segment export with the last one, not the segment start
  with its end: start-vs-end measures the response to the coefficient change, not non-steadiness,
  and fired on fully equilibrated segments.
- **Edge level rule** (Wilcox eqs 4-5): the knots held next to the Dirichlet grid edge keep a
  wrong gradient there, which offsets the whole profile inside and capped the fit at ~5%. The band
  now carries a level correction, sign-aware (flux direction from the reference gradient) and with
  the gain u_k / (u_k - u_edge).
- **Measured sensitivity is diagnostic only.** Using it as a Newton slope amplified steps where
  neighbouring knots moved against each other (Carli's oscillating-coefficient failure). The step
  is now relax * ln(g) with the pure-diffusion slope; a clearly wrong measured sign halves it.
- **Small-reference hold**: SOL knots where the reference value is below SMALL_FRAC (3%) of its
  separatrix value are held. Ratios of meaningless tiny temperatures drove the far-SOL chi to the
  floor otherwise.
- **Target widths** use the same estimator for simulation and probes (simulation sampled at the
  probe psi_N with the probe weights); a uniform-weight simulation fit over the far SOL, where the
  probes sit on their error floor, gave inconsistent widths and pushed chi up against the upstream
  signal.
- **Freshly activated channels** start with clamp 1.5; the flat-simulation push is bounded to a
  factor 2 and skipped where the whole simulated profile is low.
- Result: on the toy (hidden D with a pedestal barrier, chi_e, chi_i, 30% inward pinch, 1% noise,
  equilibration time of 120 steps, 200-step segments) the mismatch goes from 20-30% to ~2% in
  12-14 iterations and the fitted knots recover the hidden coefficients within 2-7% (median), for
  one and two temperatures.

## 10. Fuelling feedback (2026-10-14, after the user's question about the rising inventory)

The design above froze the puff. That only works when the run is already stationary at the
right fuelling, which the user's kinetic runs are not (continuous puff, inventory ramps). The
loop now has `FUEL_CONTROL`: between segments the rates of one valve are scaled by
(n_sim/n_ref at psi_N = FUEL_PSI)^-FUEL_GAIN, clamped, so the puff owns the separatrix
density while D owns the gradients (SOLPS's nesepm feedback at segment granularity).

The toy showed the expected degeneracy immediately: with only the separatrix level held, the
loop settled at 81% of the true puff with a correspondingly lower SOL D. Profiles alone do not
separate throughput from transport (Wilcox's 5x scan, Horvath's 10x D). The independent
constraint is the target density level from the probes: for a matched separatrix density, too
little density at the plates means too little SOL transport, so the SOL D knots get
-W_TLEVEL * <ln(n_sim/n_probe)> added. With that, the toy (whose target density scales with
the radial flux feeding the plate) recovers the puff from 60% to within 2% of the truth, the
profiles within 2%, the fitted knots within 3-7% median. The edge-level rule uses the level
relative to the separatrix in this mode, so the inventory does not leak into the edge band.
With `FUEL_CONTROL=0` (puff known and fixed) the absolute edge-level rule is kept and the
separatrix level is reported as an inventory diagnostic.
