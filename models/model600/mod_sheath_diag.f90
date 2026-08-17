!> Wall current and sheath potential diagnostic for the charge-conserving sheath boundary
!! condition (bcs(:)%natural%u, model600).
!!
!! The quantities are accumulated at the boundary Gauss points inside
!! mod_boundary_matrix_open.f90, where the geometry and the plasma state are available anyway,
!! and reduced and printed once per matrix construction from construct_matrix_mod.f90. Doing it
!! that way avoids duplicating the boundary element identification, at the price of a critical
!! section per boundary Gauss point, which is negligible: boundary elements are a small fraction
!! of the mesh.
!!
!! Printed per time step:
!!   I_sheath  total current through the wall using the current the sheath actually passes [A]
!!   I_Ampere  the same integral using the interior current zj = Delta*psi [A]
!! The two agree once the potential has adjusted, so their difference is a direct convergence
!! measure of the boundary condition, and I_sheath tells you whether a floating wall model would
!! change anything (it should be small compared to the current flowing in each direction).
module mod_sheath_diag

  use phys_module, only: max_bnd_types

  implicit none

  private

  public :: sheath_diag_reset, sheath_diag_add, sheath_diag_report

  real*8,  save :: sd_I_sheath(max_bnd_types) = 0.d0  !< current through the wall, sheath value [A]
  real*8,  save :: sd_I_amp(max_bnd_types)    = 0.d0  !< the same with the interior Ampere current [A]
  real*8,  save :: sd_area(max_bnd_types)     = 0.d0  !< wetted area [m^2], for averages
  real*8,  save :: sd_phi_min                 =  1.d30
  real*8,  save :: sd_phi_max                 = -1.d30
  real*8,  save :: sd_phi_sum                 = 0.d0  !< area weighted, for the mean
  real*8,  save :: sd_ratio_max               = 0.d0  !< max |j/j_sat| demanded by the interior
  real*8,  save :: sd_lim_area                = 0.d0  !< area sitting on the electron side limiter

contains

!> Zero the accumulators. Called once per matrix construction, before the element loop.
subroutine sheath_diag_reset()
  implicit none
  sd_I_sheath = 0.d0
  sd_I_amp    = 0.d0
  sd_area     = 0.d0
  sd_phi_min  =  1.d30
  sd_phi_max  = -1.d30
  sd_phi_sum  = 0.d0
  sd_ratio_max= 0.d0
  sd_lim_area = 0.d0
end subroutine sheath_diag_reset


!> Add the contribution of one boundary Gauss point on one toroidal plane.
!!
!! @param bnd_type  boundary type of the node
!! @param zj_sh     sheath current (zj units)
!! @param zj0       interior current at the same point (zj units)
!! @param zj_sat    ion saturation current (zj units)
!! @param x_lim     sheath exponent actually used, after limiting
!! @param u0        potential variable
!! @param Te0       electron temperature (JOREK units)
!! @param Bdotn     B.n
!! @param dS        surface element of this sample: ws*dl*R*(2*pi/n_plane) [m^2]
subroutine sheath_diag_add(bnd_type, zj_sh, zj0, zj_sat, x_lim, u0, Te0, Bdotn, dS)

  use constants,     only: MU_ZERO
  use phys_module,   only: F0, sheath_X_min, sheath_smooth_dX
  use mod_sheath_bc, only: sheath_norm

  implicit none
  integer, intent(in) :: bnd_type
  real*8,  intent(in) :: zj_sh, zj0, zj_sat, x_lim, u0, Te0, Bdotn, dS

  real*8 :: jn_sheath, jn_amp, phi_over_te, ratio, a_n, c_sat, vw

  if ( bnd_type .lt. 1 .or. bnd_type .gt. max_bnd_types ) return

  ! --- J.n = J_par*(b.n) = zj*(B.n)/(F0*mu0)   [A/m^2]
  jn_sheath = zj_sh * Bdotn / (F0 * MU_ZERO)
  jn_amp    = zj0   * Bdotn / (F0 * MU_ZERO)

  call sheath_norm(a_n, c_sat, vw)
  phi_over_te = ( 0.5d0*a_n*u0 - vw ) / max(Te0, 1.d-14)     ! e*Phi/(k*Te)

  ratio = 0.d0
  if ( abs(zj_sat) .gt. 0.d0 ) ratio = abs(zj0 / zj_sat)

  !$omp critical (sheath_diag_accumulate)
  sd_I_sheath(bnd_type)  = sd_I_sheath(bnd_type) + jn_sheath * dS
  sd_I_amp(bnd_type)     = sd_I_amp(bnd_type)    + jn_amp    * dS
  sd_area(bnd_type)      = sd_area(bnd_type)     + dS
  sd_phi_sum             = sd_phi_sum            + phi_over_te * dS
  sd_phi_min             = min(sd_phi_min, phi_over_te)
  sd_phi_max             = max(sd_phi_max, phi_over_te)
  sd_ratio_max           = max(sd_ratio_max, ratio)
  ! --- area where the electron side limiter is biting, i.e. where the wall is close to
  ! --- electron saturation and the characteristic is being held back
  if ( x_lim .lt. sheath_X_min + 2.d0*max(sheath_smooth_dX,1.d-3) ) sd_lim_area = sd_lim_area + dS
  !$omp end critical (sheath_diag_accumulate)

end subroutine sheath_diag_add


!> Reduce over MPI ranks and print one line. Called after the element loop; rank 0 prints.
!!
!! Every rank calls all three collectives unconditionally (no early return), and only arrays are
!! passed, so this is safe with any MPI interface and cannot deadlock.
subroutine sheath_diag_report(my_id)

  use mpi_mod

  implicit none
  integer, intent(in) :: my_id

  integer, parameter :: ns = 3*max_bnd_types + 2
  real*8  :: loc_sum(ns), glo_sum(ns), loc_max(2), glo_max(2), loc_min(1), glo_min(1)
  real*8  :: area_tot, I_sh_tot, I_am_tot, phi_mean, lim_frac
  integer :: ierr, i, i0, i1, i2

  i0 = 0                      ! offsets into the packed reduction buffer
  i1 =   max_bnd_types
  i2 = 2*max_bnd_types

  loc_sum(i0+1:i0+max_bnd_types) = sd_I_sheath
  loc_sum(i1+1:i1+max_bnd_types) = sd_I_amp
  loc_sum(i2+1:i2+max_bnd_types) = sd_area
  loc_sum(ns-1)                  = sd_phi_sum
  loc_sum(ns)                    = sd_lim_area

  loc_max(1) = sd_phi_max
  loc_max(2) = sd_ratio_max
  loc_min(1) = sd_phi_min

  call MPI_Reduce(loc_sum, glo_sum, ns, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(loc_max, glo_max,  2, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(loc_min, glo_min,  1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)

  if ( my_id .ne. 0 ) return

  area_tot = sum( glo_sum(i2+1:i2+max_bnd_types) )
  if ( area_tot .le. 0.d0 ) return                     ! the sheath BC is not active anywhere

  I_sh_tot = sum( glo_sum(i0+1:i0+max_bnd_types) )
  I_am_tot = sum( glo_sum(i1+1:i1+max_bnd_types) )
  phi_mean = glo_sum(ns-1) / area_tot
  lim_frac = 1.d2 * glo_sum(ns) / area_tot

  write(*,'(A,es11.3,A,es11.3,A,f7.2,A,f7.2,A,f7.2,A,es9.2,A,f5.1,A)')             &
    ' SHEATH: I_wall=', I_sh_tot, ' A (Ampere ', I_am_tot,                         &
    ' A)  ePhi/kTe min/mean/max=', glo_min(1), ' /', phi_mean, ' /', glo_max(1),   &
    '  max|j/jsat|=', glo_max(2), '  e-limited ', lim_frac, ' %'

  do i = 1, max_bnd_types
    if ( glo_sum(i2+i) .le. 0.d0 ) cycle
    write(*,'(A,i3,A,es11.3,A,es11.3,A,es10.3,A)')                                &
      '         bnd type', i, ': I_sheath=', glo_sum(i0+i), ' A  I_Ampere=',      &
      glo_sum(i1+i), ' A  area=', glo_sum(i2+i), ' m^2'
  enddo

end subroutine sheath_diag_report

end module mod_sheath_diag
