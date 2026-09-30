!> Prescribed local sink for MARFE studies: the (R,Z) shape and time ramp shared by the model equations
!> (models/model600/mod_elt_matrix_fft.f90) and the integrated diagnostics (diagnostics/mod_integrals3D.f90).
!> The parameters marfe_sink_* live in phys_module.
module mod_marfe_sink

  use constants,   only: PI
  use phys_module, only: marfe_sink_R, marfe_sink_Z, marfe_sink_dR, marfe_sink_dZ, marfe_sink_t_start, marfe_sink_t_ramp

  implicit none

contains

  !> Shape of the prescribed local sink (MARFE studies) at a given position and time: a Gaussian in (R,Z) around
  !> (marfe_sink_R, marfe_sink_Z) times a smooth ramp in time, between 0 and 1. The rates marfe_sink_nu_* multiply it.
  pure real*8 function marfe_sink_shape(R, Z, time)

  ! --- Routine parameters.
  real*8, intent(in) :: R
  real*8, intent(in) :: Z
  real*8, intent(in) :: time

  ! --- Local variables
  real*8 :: ramp

  marfe_sink_shape = 0.d0

  if ( (marfe_sink_dR .le. 0.d0) .or. (marfe_sink_dZ .le. 0.d0) .or. (time .lt. marfe_sink_t_start) ) return

  if ( time .lt. marfe_sink_t_start + marfe_sink_t_ramp ) then
    ramp = 0.5d0 - 0.5d0*cos(PI * (time - marfe_sink_t_start) / marfe_sink_t_ramp)
  else
    ramp = 1.d0
  endif

  marfe_sink_shape = ramp * exp( - ((R - marfe_sink_R)/marfe_sink_dR)**2 - ((Z - marfe_sink_Z)/marfe_sink_dZ)**2 )

  end function marfe_sink_shape

end module mod_marfe_sink
