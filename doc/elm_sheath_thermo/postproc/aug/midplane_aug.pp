namelist input
si-units
for-si-units
set linepoints 400
set surfaces 150
# Toroidally averaged outer- and inner-midplane profiles through the ELM: pedestal flattening, the pedestal current,
# E_r, and Psi(R,t) for the toroidal electric field E_phi = -(1/R) dPsi/dt in the SOL (Kallenbach's loop-voltage
# EMF). Use the SAME time list as plates_circuit_aug.pp (or a subset of it): the comment lines of these files carry
# only the step number, and the script takes the step -> time map from the plate files.
# Start 0.3 ms before T_ON; cadence 2e-5 s.
for time T_PRE to T_END by 2.0e-5 do
  expressions R Z Psi_N Psi T_e T_i ne pres currdens Er Phi
  mark_coords 1
  midplane outer
  midplane inner
  qprofile
done
# Axis and X-point positions (equil_params<range>.dat) at one time, for the plots:
for time T_ON do
  equil_params
done
