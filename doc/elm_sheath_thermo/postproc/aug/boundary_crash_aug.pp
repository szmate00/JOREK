namelist input
si-units
for-si-units
set nsub_bnd 4
# Crash window at the restart-file cadence (nout = 10 steps of tstep 1 = 5.7 us): every 0.01 ms.
# Replace T_ON by the onset time in seconds (same value for BOTH runs, e.g. 0.0152), end = T_ON + 0.6 ms.
for time T_ON to T_END by 1.0e-5 do
  expressions R Z Psi_N bnd_dl bnd_type_lo bnd_type_hi Bnorm Jn_par Jn_sat Jn_sheath Jnorm dPhi_float Phi T_e T_i ne heatF_total heatF_par_cv heatF_prp_cv heatF_tot_cd
  boundary_quantities 0.d0 3.14159265 32
done
