namelist input
si-units
for-si-units
set nsub_bnd 4
# Plate data for the physics checks (util/elm_physics_checks_aug.py): charge/energy asymmetry, current -> potential
# -> T_e timing, flux-tube circuit (sheath drops, thermoelectric EMF, inductive EMF from dPsi_dt), filament current
# closure and the effective sheath heat transmission. Same T_ON / T_END (seconds) in run A and run B.
# Cadence 1e-5 s = every second restart file at nout 10 / tstep 1; use 5.7e-6 for every file.
# dPsi_dt is written in Wb/rad per JOREK time unit even with si-units (the script divides by sqrt(mu0 rho0)).
for time T_ON to T_END by 1.0e-5 do
  expressions R Z Psi_N Psi dPsi_dt bnd_dl bnd_type_lo bnd_type_hi Bnorm Btor B_abs Jn_par Jn_sat Jn_sheath Jnorm dPhi_float Phi T_e T_i ne vpar Mach_par heatF_total heatF_sheath heatF_par_cd heatF_par_cv heatF_prp_cv partF_total
  boundary_quantities 0.d0 3.14159265 32
done
