namelist input
# Per-term output of the u row (vorticity equation = charge conservation) on the finite elements, for
# util/rhs_tube_budget.py: which term balances the parallel current along a flux tube. Equation index 2 = u.
set_postproc_dir postproc_rhs_u
for step 1300 do
  RHS_terms_vtk 2
done
