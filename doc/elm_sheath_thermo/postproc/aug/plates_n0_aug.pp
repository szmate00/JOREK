namelist input
si-units
set nsub_bnd 4
# Axisymmetric (n_tor = 1) sheath-row check: plate quantities at every available restart file, one toroidal plane.
# Set the step range from `echo timesteps | ./jorek2_postproc`; the loop takes the closest restart for every
# requested step and skips duplicates, so a `by` equal to the restart interval (nout) is enough.
# Analyse with: python3 util/sheath_row_check.py LABEL='postproc/boundary_quantities*.dat'
for step STEP0 to STEP1 by NOUT do
  expressions R Z Psi_N bnd_dl bnd_type_lo bnd_type_hi Bnorm B_abs Jn_par Jnorm Jn_sat Jn_sheath dPhi_float Phi T_e T_i ne Mach_par vpar heatF_total
  boundary_quantities 0.d0 0.d0 1
done
