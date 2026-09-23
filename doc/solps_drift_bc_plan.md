# Plan: SOLPS-ITER drift-compatible wall BCs in model600 (floating u, later sheath j)

Branch `floating-u-sheath-flux`. Status: PHASE 1 IMPLEMENTED (`mach1_weak_drift_style = 1`, `mach1_weak_qalf_min`,
per-node active set in the weak row, `mod_bohm_active.f90`), serially tested, not yet run. Phase 2 (item 7) open. Written 2026-09-23 after reading the
SOLPS-ITER sources in `SOLPS-ITER/` (standard B2.5 and the wide-grid branch `feature/wg-release` of the B2.5
submodule) and an independent review of this branch.

## Campaign rules (apply to every item below)

- No grid changes.
- Every wall BC drift-compatible and consistent with every other one: the Vpar row, the rho/Ti/Te sheath fluxes,
  the energy transmission, the kinetic (ncs) recycling and the potential BC impose the same condition at the same
  points, with the same drift.
- Out of the box: no reliance on tg_num, shock capturing, extra diffusion, or thresholds fitted to crashing runs.
  The factor 2 of the drift bound is SOLPS's model constant (accepted). The field-alignment threshold is an input
  parameter to be tested.

## The SOLPS reference (wide-grid branch, `src/sources/b2stbc_phys.F`)

Per boundary face, with `b_n = b.n` (signed incidence), `vE_n = vE.n` (outward ExB normal speed), `cs`:

| SOLPS | definition |
|---|---|
| drift bound (`vbc_is`) | `vE_r = clamp(vE_n, -2 cs |b_n|, +2 cs |b_n|)` |
| parallel normal speed (BCMOM=13, non-marginal) | `vbc = max(cs |b_n| - vE_r, b_n u_par,interior)`: an INEQUALITY, the interior flow passes if it is faster |
| particle flux (BCCON=14) | `Gamma = n * max(cs |b_n|, vbc + vE_n)`, full `vE_n` in the total; BCCON=14 requires BCMOM=13 |
| electron energy (BCENE=15, style 1) | parallel electron thermal flux `n vte/sqrt(2 pi) |b_n| exp(-e dphi/Te)` times `(1+gamma_e+ENEPAR) Te + (1-gamma_e) e dphi` |
| potential (BCPOT=11) | current continuity, ion current `e Z Gamma` with the SAME `Gamma` as BCCON=14, electron current `(1-seec) e n vte/sqrt(2 pi) |b_n| exp(-e dphi/Te)`; BCPOT=11 requires BCCON=14 |
| field-aligned faces, `|b_n| < Qalfmin` (default 1e-3) | no Bohm: u_par extrapolated from the interior; flux = leakage `n cs Qalfmax` + `max(0, outward ExB flux)`; potential extrapolated |

The code checks at start-up that the set is used together (`BCCON=14, BCMOM=13, BCENE/I=15 style 1,
BCPOT=3|11, MOMPAR(:,:,2) >= 0.5`) and stops or warns otherwise.

## Mapping onto model600 (per wall Gauss point, JOREK units)

Notation as in `mod_boundary_matrix_open.f90`: `Bn = B_pol.n`, `b_n = bdotn`, `cs = cs0`, `vE_n = mw_vEn`,
`vn = Bn*Vpar + vE_n`.

1. **Drift bound**: `vE_r = max(-2 cs|b_n|, min(2 cs|b_n|, vE_n))`. Used only to build the Vpar target, as in SOLPS.
   Jacobian: the u column of `vE_r` exists only while the bound is inactive.
2. **Vpar row (inequality, implemented with a PER-NODE active set in the weak row; the full nodal rewrite is the
   fallback)**: target `T = cs|b_n| - vE_r` (parallel normal speed, may be negative down to
   `-cs|b_n|`, as in the wide-grid branch). Residual `res = Bn*Vpar - T`, imposed only where `res < 0` (the flow is
   slower than drift-compatible Bohm); where `res >= 0` the row is inactive and Vpar follows its own equation. Weight
   and assembly as the present weak row. The active set is decided per node from the node's own state (OR over its
   wall edges), collected during one matrix construction and used in the next (a one-step lag; all nodes active
   before the first collection), and differentiated exactly on the branch taken. Replaces both the equality and the
   `Vpar = 0` pin of the present code.
3. **Wall fluxes**: `Gamma = n * max(vn, cs|b_n|)`, energy `gamma_sh * T * Gamma`: the present sheath-set flux with
   floor `cs|b_n|` at every point (full `vE_n` in `vn`, as in SOLPS BCCON=14). No change of form, only the cut
   dependence goes.
4. **Field-aligned points**, `|b_n| < mach1_weak_qalf_min` (new input parameter): no Vpar row, flux floor 0, i.e.
   `Gamma = n * max(vn, 0)`, plus the existing c_angle floor (the analogue of SOLPS's leakage `n cs Qalfmax`).
   Default proposal: `1.d-3` as in SOLPS; A/B against `0.d0` (no threshold: the (B.n)^2 weighting of the row and the
   bound alone).
5. **Kinetic recycling** (`mod_particle_wall_interaction.f90`): the same `Gamma` point by point:
   `n_e * max(v_n,tot, c_s cos_alpha)` on mach1 edges above the threshold, `n_e * max(v_n,tot, 0)` below it, plus
   the c_angle floor, fluid sound speed. Same edge test as the fluid (already in place).
6. **min_sheath_angle** keeps only its original role (the c_angle floor / leakage); it no longer switches the drift.
7. **Floating potential, made consistent with the flux (phase 2, decision needed)**: with zero local current SOLPS's
   BCPOT=11 gives `e*dphi/Te = ln[(1-seec) vte |b_n| / (sqrt(2 pi) Gamma/n)]`, i.e. with `Gamma = n max(vn, cs|b_n|)`
   `Lambda_eff = Lambda - ln( max(vn, cs|b_n|) / (cs|b_n|) )`. The present row uses a constant `Lambda` everywhere, so
   the potential ignores the drift-enhanced outflow. Consistent version: `Lambda_eff` in the floating row (nonlinear,
   needs Vpar/u/T columns at the nodes), and at field-aligned points no potential condition (SOLPS extrapolates it).

## Switches

- New: `mach1_weak_drift_style` (0 = present equality, unbounded drift; 1 = SOLPS bounded inequality), mirrors
  SOLPS's `drift_style`. Needs `mach1_weak_drift = .t.`.
- New: `mach1_weak_qalf_min` (field-alignment threshold in `|b.n|`), used only with style 1.
- `mach1_weak_drift_cut`: redundant with style 1; refused together with it at setup. Keep or remove: user decision.
- Setup consistency check like SOLPS: style 1 requires `mach1_weak`, `mach1_weak_drift`; warn if `floating_u` is on
  without the phase-2 potential (item 7).

## Implementation tasks

| # | file | task |
|---|---|---|
| 1 | `models/phys_module.f90`, `preset_parameters.f90`, `mod_log_params.f90`, `broadcast_phys.f90`, `model600/initialise_parameters.f90` | the two switches, defaults, log, broadcast, namelist, setup checks |
| 2 | `models/model600/mod_boundary_matrix_open.f90` | bound `vE_r`, inequality row (active set, columns), threshold branch, flux floor per point |
| 3 | `particles/mod_particle_wall_interaction.f90` | same point test (bound is not needed there: the flux uses the full `vn`), threshold instead of the cut |
| 4 | `models/model600/mod_wall_diag.f90` | columns: fraction of wall with the inequality active, with the bound active, below the threshold |
| 5 | `doc/floating_u.md` | the model, the switches, the SOLPS correspondence |
| 6 (phase 2) | `models/model600/mod_boundary_conditions.f90`, `mod_floating_u.f90` | `Lambda_eff` floating row, potential free at field-aligned points |

## Tests

- Serial harness (production assembler against stubs, as for the sheath-set flux): FD of every column in each branch:
  bound inactive / active with either sign, inequality active / inactive, below the threshold; the drift bound
  continuous across `|vE_n| = 2 cs|b_n|`; style 0 bit-identical to the present code; develop path bit-identical.
- Consistency check: the fluid `Gamma` and the kinetic formula evaluated on the same state agree point by point.
- Runs (same equilibrium, namelist of the 808 run with `mach1_weak_drift_cut` removed):
  1. result of the run in progress (980cf2177) first: `[wall pt] minrho / maxMach`, `[wall] sink`;
  2. style 1, `mach1_weak_qalf_min = 1.d-3`;
  3. style 1, `mach1_weak_qalf_min = 0`;
  4. the better of 2/3 with `keep_current_prof_confined`;
  5. phase 2 (item 7) if 2-4 run.

## Later: sheath j

SOLPS BCPOT=11 is the sheath-current BC of `sheath-j-clean` with two changes to port: the ion current from the same
drift-compatible `Gamma` (not `c_sat rho cs/|B|`), and the electron current from the parallel electron thermal flux
with `|b_n|`. Field-aligned points: potential extrapolated. To be planned after the floating results.

## Open decisions

1. Default of `mach1_weak_qalf_min`: SOLPS's `1e-3`, or `0`.
2. Keep or remove `mach1_weak_drift_cut` once style 1 works.
3. Phase 2 (drift-consistent `Lambda_eff`, potential free at field-aligned points): yes / no / after the runs.
4. Negative Vpar target (parallel reversal down to `-cs|b_n|` when the drift is strongly outward, as in the wide
   grid) or never reversed (`max(0, .)`, as in standard SOLPS BCMOM=13).
