namelist input
si-units
for-si-units
# Volume state at matched SI times (seconds): one loop per time, phi = 0 and half a period (pi/2).
# 2.5 mm over the AUG domain. Replace T1.. by the same list in both runs (pre-onset, onset, power peak, +0.2 ms).
for time T1 do
  expressions R Z Psi_N ne T_e T_i Phi pres_e currdens vpar Mach_par Er Btor
  rectangle 0.95 2.35 561 -1.40 1.35 1101 0.
  rectangle 0.95 2.35 561 -1.40 1.35 1101 1.5707963
done
for time T2 do
  expressions R Z Psi_N ne T_e T_i Phi pres_e currdens vpar Mach_par Er Btor
  rectangle 0.95 2.35 561 -1.40 1.35 1101 0.
  rectangle 0.95 2.35 561 -1.40 1.35 1101 1.5707963
done
