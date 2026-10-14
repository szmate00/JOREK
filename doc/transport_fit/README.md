# Automated transport-profile fitting for model600 (`util/transport_fit.sh`)

The loop adjusts the tabulated `D_perp_file` and `ZK_perp_file` (or `ZK_e_perp_file` /
`ZK_i_perp_file`) of a kinetic model600 run until the simulated outer-midplane profiles stay on
the step-0 profiles and the simulated target profiles follow Langmuir-probe data.  You set up
the run as usual, put a probe file next to it, start one script, and come back to fitted
profile files, a log and one figure per iteration.  It runs unattended on the cluster and
resumes if it is interrupted.

The method and its background are in `PLAN.md` in this directory.

## 1. What you need in the run directory

Everything flat in one directory (`$BASE` in the script):

| file | what |
|---|---|
| `input` (or your namelist name) | the run's namelist with everything switched on that you want on during the fit: neutrals, impurities, sheath BCs, drifts. `use_zkperp_times_density = .true.` is required. |
| `kinetic_main`, `jorek2_postproc` | the executables, built with the `with_TiTe` you want (see TEMPERATURES below) |
| `jorek000000.h5` + the newest `jorek??????.h5` | the reference step (the initial profiles) and the state the loop restarts from (`nout` must have written them) |
| `part_restart.h5` | optional; if missing the first segment initialises the particles |
| the profile files named in the namelist | or tanh coefficients, from which the loop bootstraps files |
| `probes.dat` | optional; `R Z ne Te [ne_err Te_err]`, one probe per line, `#` comments, R and Z in m, ne in 1e19 m-3, Te in eV. Without it the fit is upstream-only and the targets are only reported. |
| `transport_fit.py`, `transport_fit.sh` | copied from `util/` |

**Fuelling.** The fit targets a stationary phase: constant inventory, puff balanced by pump
and wall. Set the pump from the experiment. For the puff there are two modes:

- `FUEL_CONTROL=1` (default): the loop scales the rates of one valve (`PUFF_KEY`) between
  segments so that the outer-midplane density at `FUEL_PSI` (the separatrix) stays on its
  step-0 value, SOLPS-style. The probe density level then fixes the throughput: too little
  density at the plates for the matched separatrix density means more SOL transport, and
  the SOL D knots are moved accordingly (`W_TLEVEL`). The converged puff rate is a result
  of the fit, printed in the report, and the final state is stationary at that fuelling,
  which is the state to start a production study (e.g. a puff ramp) from.
- `FUEL_CONTROL=0`: puff fixed at what you set. The separatrix density level is reported
  every iteration; a systematic drift there is inventory, not transport, and the fit cannot
  fix it. Use this only when the puff is known and the run is already stationary.

Without a probe file the throughput is unconstrained either way: the same profiles can be
held with different (puff, D) pairs (Horvath 2023, Wilcox 2026).

## 2. Settings

All settings are variables at the top of `transport_fit.sh`. The ones you normally touch:

| variable | meaning | default |
|---|---|---|
| `BASE` | the run directory | |
| `RUN_MODE` | `slurm` (sbatch --wait), `local` (mpirun), `toy` (test without JOREK) | slurm |
| `TEMPERATURES` | 1 = one-temperature build (`T`, `ZK_perp_file`); 2 = `with_TiTe` (`T_e`, `T_i`, `ZK_e_perp_file`, `ZK_i_perp_file`) | 2 |
| `PROBE_FILE` | probe file, empty for upstream-only | probes.dat |
| `NSTEP_SEG`, `NSTEP_MAX` | steps per segment; doubled automatically while the state is not steady | 200 / 800 |
| `MAX_ITER`, `TOL` | iteration cap; converged when all channels are within TOL (relative rms) for two iterations | 12 / 0.05 |
| `RELAX`, `CLAMP` | relaxation and maximum factor per knot and iteration; use 0.5 / 1.5 with drifts | 0.7 / 2.0 |
| `PSI_MIN` | knots below this psi_N are held (the core does not equilibrate within a segment) | 0.7 |
| `W_TARGET`, `TE_TARGET_WEIGHT` | weight of the measured target width in the SOL knots; T_e from probes gets the extra factor | 0.5 / 0.5 |
| `CHI_EDGE_MAX` | cap on chi next to the grid edge in the first three iterations (m2/s) | 10 |
| `TSTEP` | time step for every segment (production value, no ramp); empty keeps `tstep_n(1)` | |
| `FUEL_CONTROL`, `PUFF_KEY`, `FUEL_PSI`, `FUEL_GAIN`, `W_TLEVEL` | puff feedback on the separatrix density and the throughput term from the probe density level (section 1) | 1 / valve 1 / 1.0 / 0.5 / 0.5 |
| the `JOB_*` block | the slurm header of the job script the loop writes (`tfit_job.sh`) | |
| `EXTRA_CFG` | further KEY=VALUE pairs, e.g. `PLATE_ZMAX=-0.5` when more than the two targets see the field, `KNOTS=0,0.5,...` for your own knot set, `EDGE_HOLD=0.08`, `SMALL_FRAC=0.03` | |

The segment length is the main cost knob. A segment has to be long enough for the edge to
re-equilibrate after the coefficients changed; the loop measures that (section 4) and doubles
the segment when needed, so start short.

## 3. Running

```
cp util/transport_fit.py util/transport_fit.sh $BASE/
# edit the settings block at the top of $BASE/transport_fit.sh
cd $BASE && nohup ./transport_fit.sh > transport_fit.out 2>&1 &
```

The script first runs the pre-flight `check` (namelist, files, units, probe file), exports the
reference step, finds the target plates and projects the probes onto them, and prints all of it.
Read that block once: it tells you the D and chi in m2/s that the namelist currently means, how
many probes landed on which plate, and where the strike points are. Then it loops:

1. copy the newest `jorek??????.h5` to `jorek_restart.h5`, set `restart`, `nstep_n`,
   `restart_particles` in the namelist, run one segment;
2. export the outer-midplane cut and the boundary quantities for the middle and the last step
   of the segment with `jorek2_postproc`;
3. update the profile files (`transport_fit.py update`), append a line to `tfit_log.csv`,
   write `tfit_iterNN.png`;
4. stop when converged, when the coefficients stop changing, at `MAX_ITER`, or after three
   failed segments.

Re-running the script resumes from `tfit_state.json`. `tfit_driver.log` has everything the
loop printed. The original profile files are kept as `<file>.user`, every iteration's input as
`<file>.iterNN`.

## 4. Reading the output

Per iteration the log prints, per channel (`ne`, `T_e`, `T_i` or `T`):

- **mismatch to reference**: relative rms difference between the simulated and the step-0
  profile over the fitted psi_N range, also split into bands (0.7-0.9, 0.9-1.0, 1.0-1.1, SOL).
- **drift over the second half of the segment**: how much the profile still moved. If it is
  larger than half the mismatch the state is not steady, no update is made and the segment is
  rerun twice as long (exit code 20 in the log).
- **max factor this step** and **D/chi at psi_N 0.9/1.0/1.1** in m2/s.
- **flagged knots** with a reason:
  - `flat-ref`: the reference profile has no gradient here, nothing to fit, held;
  - `small-ref`: reference value below a few percent of its separatrix value (far SOL T),
    held;
  - `flat-sim`: the simulation is flat where the reference is not, transport reduced by a
    bounded step;
  - `low-sim`: the simulated profile is far below the reference here, the gradient says
    nothing about the coefficient, held;
  - `wrong-sign`: the last change moved the gradient the wrong way (drift-dominated
    region or neighbour interaction), the step is halved. A knot flagged `wrong-sign` for many
    iterations is a region that D or chi cannot control; say so in the write-up instead of
    fighting it;
  - `edge-level`: the band next to the grid edge, which only carries the level correction.
- **density level at psi_N 1.00**: simulation over reference. With `FUEL_CONTROL=1` the puff
  factor applied this iteration and the cumulative factor since the start follow.
- **target inner/outer ne/Te**: SOL-side e-folding length in psi_N of the simulated profile
  along the plate versus the probes, the number of probes used, the rms level difference and
  the peaks. The width enters the SOL knots with weight `W_TARGET`. With `FUEL_CONTROL=1` the
  density level at the probes also moves the SOL D knots (throughput); otherwise it is only
  reported.

The figure has the same content: top row reference / mid-segment / end-of-segment profiles,
middle row the coefficients in m2/s (dotted = your original file, dots = knots, triangles =
held knots, crosses = flagged), bottom row the target profiles against distance from the strike
point with the probes.

`transport_fit.py report --cfg tfit.cfg` prints the iteration table at any time.

## 5. When it stops early

- **three failed segments**: the coefficients of the last successful iteration are back in
  place, the clamp was halved twice. Look at the flagged knots in the last figure and at
  `tfit_crash_N/logfile`. If the job hit the wall time rather than crashing, raise `WALLTIME`
  or lower `NSTEP_SEG` and rerun.
- **not steady at NSTEP_MAX**: the loop updates anyway with half the relaxation and says so.
  Raise `NSTEP_MAX` if the drift stays large.
- **check fails**: the message says what to fix (missing file, `use_zkperp_times_density`
  off, non-positive values in a profile file, profile table ending before the grid's psi_N).

## 6. Limits you should know

- The private flux region reads the SOL table at the mirrored psi_N, so it cannot be fitted
  separately; the PFR side of the target plots is a check only.
- Knots within `EDGE_HOLD` of the grid edge carry only the level correction, because the
  Dirichlet boundary there dominates the gradient.
- The core below `PSI_MIN` is held; its profiles barely move within a segment anyway.
- Impurities use the main-ion D (no file for `D_perp_imp` in model600); neutral diffusion is
  physics and is not fitted.
- Several (D, chi) sets hold the same profiles when fuelling is free; keep it fixed.

## 7. Testing

Without JOREK:

```
python3 util/transport_fit_toy.py selftest                  # two temperatures, pinch, noise, slow equilibration
python3 util/transport_fit_toy.py selftest --temperatures 1
python3 util/transport_fit_toy.py selftest --pinch 0 --tau-eq 10
```

Each runs the full loop (the real `transport_fit.py` commands) against a 1D model with a hidden
D and chi, from constant first guesses, and passes when the profiles match within 5% and the
fitted knots recover the hidden coefficients within a factor 1.3 (typically 2-7%).

The driver itself can be exercised with `RUN_MODE=toy`: `transport_fit_toy.py init --cfg` a
small cfg (`RUN=...`, `NAMELIST=input`, `TEMPERATURES=2`) creates a fake run directory with a
namelist, step files, exports and a probe file; point `BASE` at it and run the script. A file
`toy_fail_once` in that directory makes the next segment fail, which exercises the crash path.

With JOREK, the first real test is a digital twin: run an existing case to a steady state with
known tabulated profiles, use its step files as "experiment" (its last step as the reference
step and a probe file made from its boundary export), start from constant coefficients and
check that the loop recovers the known files. Then the real case.
