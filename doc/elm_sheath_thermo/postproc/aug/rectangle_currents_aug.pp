namelist input
si-units
for-si-units
# Current composition (util/check_current_conservation.py): the tube budget needs Psi (unmirrored), the poloidal field,
# the parallel current and the full poloidal current JR JZ (parallel + diamagnetic), plus ne T_e T_i for the pressure,
# Phi and T_e for the thermoelectric estimate. One toroidal plane (axisymmetric phase) or the n=0 average later.
# Set T1 in seconds to the time of the boundary_quantities step to compare with (the script matches by t_now).
for time T1 do
  expressions R Z Psi Psi_N BR BZ Btor currdens JR JZ Phi T_e T_i ne
  rectangle 0.95 2.35 561 -1.40 1.35 1101 0.
done
