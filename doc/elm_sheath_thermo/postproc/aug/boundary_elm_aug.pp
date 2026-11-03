namelist input
si-units
for-si-units
set nsub_bnd 4
# AUG #33616 (n_tor = 13, n_period = 2): one toroidal period is pi; 32 planes resolve n <= 12.
# Times in SECONDS. The loop stops at the last restart file present, so rerun it as the job advances;
# steps already written are skipped ("Step number already selected") only within one loop, so give
# each rerun a start time after the end of the previous file.
# Loop 1: linear phase, every 0.05 ms.
for time 0.0118 to 0.0300 by 5.0e-5 do
  expressions R Z Psi_N bnd_dl bnd_type_lo bnd_type_hi Bnorm Jn_par Jn_sat Jn_sheath Jnorm dPhi_float Phi T_e T_i ne heatF_total heatF_par_cv heatF_prp_cv heatF_tot_cd
  boundary_quantities 0.d0 3.14159265 32
done
