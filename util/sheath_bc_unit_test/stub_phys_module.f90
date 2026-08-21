!> Minimal stand-in for phys_module, containing only the variables that
!! models/model600/mod_sheath_bc.f90 reads. It exists so that the sheath
!! characteristic can be unit tested on a laptop without building JOREK
!! (no MPI, no solver, ~2 s). Not used by the JOREK build itself: the util/
!! directory is not in DIRS in the Makefile.
module phys_module
  implicit none
  real*8  :: F0, GAMMA, central_density, central_mass
  real*8  :: sheath_V_wall, sheath_Lambda, sheath_X_min, sheath_smooth_dX, sheath_min_bn
  real*8  :: sheath_sat_slope
  logical :: sheath_Lambda_local
end module phys_module
