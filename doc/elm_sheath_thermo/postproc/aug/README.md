# AUG #33616 ELM comparison (run A: reference, run B: sheath_j + thermoelectric_ohm), n_tor = 13

Both runs are sampled at the SAME SI times with `for-si-units` + `for time`, so steps never have to match.
Run in each run directory (adjust `namelist input` to the namelist file name):

    ./jorek2_postproc < boundary_elm_aug.pp      # -> boundary_quantities_sAAAAA..sBBBBB.dat (plate state, 32 planes)
    ./jorek2_postproc < boundary_crash_aug.pp    # after the onset is known: T_ON / T_END in seconds, identical in A and B
    ./jorek2_postproc < rectangle_aug.pp         # volume state at a few matched times (T1, T2, ... in seconds)
    echo timesteps | ./jorek2_postproc           # step <-> time table of the restart files present

`macroscopic_vars.dat` is written by the running job (produce_live_data default): per-harmonic energies,
thermal energy inside the separatrix, boundary fluxes. No postproc needed for it.

Analysis (the script merges several boundary files per run and drops duplicate steps):

    python3 util/compare_elm_aug.py \
        A=runA/postproc/boundary_quantities_*.dat B=runB/postproc/boundary_quantities_*.dat \
        --macro A=runA/macroscopic_vars.dat B=runB/macroscopic_vars.dat \
        --times 14.0 15.0 15.5 --fit-window 13.0 14.2 -o doc/elm_sheath_thermo/figs_aug

Quote the file name patterns; the script expands them. `--align onset` shifts every trace to its own onset.

## Physics checks against the ELM-current literature (2026-11-02)

Command files (same `T_ON`, `T_END`, `T_PRE` in seconds in BOTH run directories; comments with `#` only):

    ./jorek2_postproc < plates_circuit_aug.pp     # plates every 1e-5 s: + Psi, dPsi_dt, Btor, B_abs, vpar, Mach_par, heatF_sheath, partF_total
    ./jorek2_postproc < midplane_aug.pp           # toroidally averaged outer/inner midplane profiles, q profile, equil_params

`plates_circuit_aug.pp` writes `boundary_quantities<range>.dat`, `midplane_aug.pp` writes `exprs_outer-midplane<range>.dat`,
`exprs_inner-midplane<range>.dat`, `qprofile<range>.dat`, `equil_params<range>.dat`. The midplane comment lines carry only
the step number, so use the same (or a subset of the) time list as the plate file: the script maps step -> time with
the plate data of the same run.

    python3 util/elm_physics_checks_aug.py \
        A='runA/postproc/boundary_quantities*.dat' B='runB/postproc/boundary_quantities*.dat' \
        --macro A=runA/macroscopic_vars.dat B=runB/macroscopic_vars.dat \
        --midplane A='runA/postproc/exprs_outer-midplane*.dat' B='runB/postproc/exprs_outer-midplane*.dat' \
        --window 0.4 --pre 0.1 --L-par 60 -o doc/elm_sheath_thermo/figs_aug

Output: `checks_charge_energy.pdf` (Q and E per plate, E_out/E_in vs charge; AUG forward field: Q_out < 0, E_in > E_out),
`checks_timing.pdf` (10 %/50 % crossings of |I|, Phi, T_e, q, footprint, E_mag, -dW/dt; Pitts 2003, Takahashi 2004),
`checks_circuit_<run>.pdf` (per Psi_N bin: currents at both plates and closure, T_e, Phi_in - Phi_out, 0.71 dT_e,
inductive EMF from dPsi_dt x L_par, remainder = resistive drop; Staebler-Hinton, Kallenbach), `checks_jv_law.pdf`,
`checks_filaments.pdf` (I+/I-, closure index, largest n>0 lobe vs Vianello's 1.5 kA, |j|/j_sat, corr(j, T_e), footprint,
toroidal spectrum), `checks_sheath_gamma.pdf` (P / sum j_sat T_e vs gamma_sheath; PIC 5-6 x rise), `checks_pedestal.pdf`
(pedestal flattening, V_loop = -2 pi dPsi/dt at the midplane and from the plate dPsi_dt), `checks_numbers.txt`.
Conventions: Jn > 0 = ions collected by the plate (the j-V figure confirms it: the sheath run follows 1 - exp(-x));
dPsi_dt is per JOREK time unit in the postproc output, the script divides by sqrt(mu0 rho0); the inductive EMF is an
estimate (plate E_par taken as the tube average, --L-par).

## Current composition along flux tubes (2026-11-03)

Why the two plates carry the same sign: with a self-consistent sheath the plate current of the axisymmetric B run is a
"plate closing current" in the sense of Rozhansky et al. (NME 25 (2020) 100840): the vertical grad-B (diamagnetic)
current of the divertor legs is closed by parallel current to BOTH plates (Zhao et al. arXiv:2601.14749, source
d j_par/ds = (1/BR) dp/dZ), while the thermoelectric plate-to-plate current is the small antisymmetric part. The
per-tube rule I_in = -I_out holds only without cross-field currents. util/plate_psi_balance.py shows the profiles and
their sum from the boundary files; the decomposition needs the rectangle output:

    ./jorek2_postproc < rectangle_currents_aug.pp        # T1 = time of the boundary step, writes exprs_..._s<step>.h5

    python3 util/check_current_conservation.py --rect runB/postproc/exprs_*_s01300.h5 \
        --boundary 'runB/postproc/boundary_quantities*.dat' --F0 4.23130791805 --plate-zmax -0.5 \
        --exprs "R Z Psi Psi_N BR BZ Btor currdens JR JZ Phi T_e T_i ne" \
        --bands 1.005:1.015 1.02:1.04 1.05:1.10 0.97:0.995 --sheath-lambda 3 -o doc/elm_sheath_thermo/figs_aug

Output: the budget table and current_conservation.pdf as before, plus the composition table and
current_composition.pdf: per band I_in, I_out, through = (I_out - I_in)/2 (thermoelectric + PS reaching the plates),
closing = I_in + I_out (PCC), the side current split into diamagnetic below / above the X-point (PCC type 1 / 2) and
the polarisation-viscous residual below / above, the parallel length, the Spitzer tube resistance, Lambda dT_e,
Phi_out - Phi_in and the thermoelectric estimate I_th = Lambda dT_e / R_tube (--eta-si for the run's eta). The
boundary file must carry T_e and Phi for the last three columns (plates_n0_aug.pp / plates_circuit_aug.pp do).
If the rectangle has no JR JZ, ne T_e T_i are enough: the script rebuilds J_par b_pol + B x grad p / B^2.
