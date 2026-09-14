!> Wall diagnostics for the floating potential and the weak Bohm condition (floating_u_diag).
!!
!! Accumulated at the wall Gauss points during matrix construction, reduced over MPI ranks and
!! printed once per matrix construction on rank 0, one line per boundary type:
!!
!!   inflow   fraction of the wall length where the total normal flow is inward
!!   vE.n     largest outward ExB normal speed [m/s] and where
!!   mom      normalised weighted moment of the Bohm residual, |sum Bn*res*dl| / sum |Bn|*cs*dl,
!!            which is what the Galerkin row imposes (pointwise |res| is not controlled at grazing
!!            incidence and is deliberately not reported)
!!   max M    largest wall Mach number |Vpar*B|/cs: the parallel flow the momentum row demands
!!   min rho, min Ti, min Te   at the wall Gauss points, and where (rho, Te)
module mod_floating_diag

  implicit none
  private

  public :: floating_diag_reset, floating_diag_add, floating_diag_report

  integer, parameter :: nt = 30            !< max_bnd_types
  real*8, save :: s_len(nt), s_in(nt), s_mom(nt), s_den(nt)
  real*8, save :: x_ven(nt), x_ven_R(nt), x_ven_Z(nt)
  real*8, save :: n_rho(nt), n_rho_R(nt), n_rho_Z(nt)
  real*8, save :: n_Te(nt),  n_Te_R(nt),  n_Te_Z(nt)
  real*8, save :: x_mach(nt), n_Ti(nt)

contains


subroutine floating_diag_reset()
  implicit none
  s_len = 0.d0 ; s_in = 0.d0 ; s_mom = 0.d0 ; s_den = 0.d0
  x_ven = -huge(1.d0) ; x_ven_R = 0.d0 ; x_ven_Z = 0.d0
  n_rho =  huge(1.d0) ; n_rho_R = 0.d0 ; n_rho_Z = 0.d0
  n_Te  =  huge(1.d0) ; n_Te_R  = 0.d0 ; n_Te_Z  = 0.d0
  x_mach = 0.d0 ; n_Ti = huge(1.d0)
end subroutine floating_diag_reset


!> One wall Gauss point. Quantities in JOREK units except where noted; dl includes the weight.
subroutine floating_diag_add(bnd_type, dl, vn, vEn, Bn, res, cs_bn, mach, rho, Ti, Te, R, Z)
  implicit none
  integer, intent(in) :: bnd_type
  real*8,  intent(in) :: dl, vn, vEn, Bn, res, cs_bn, mach, rho, Ti, Te, R, Z
  if ( bnd_type .lt. 1 .or. bnd_type .gt. nt ) return
  !$omp critical (floating_diag)
  s_len(bnd_type) = s_len(bnd_type) + dl
  if ( vn .lt. 0.d0 ) s_in(bnd_type) = s_in(bnd_type) + dl
  s_mom(bnd_type) = s_mom(bnd_type) + Bn * res * dl
  s_den(bnd_type) = s_den(bnd_type) + abs(Bn) * cs_bn * dl
  if ( vEn .gt. x_ven(bnd_type) ) then
    x_ven(bnd_type) = vEn ; x_ven_R(bnd_type) = R ; x_ven_Z(bnd_type) = Z
  endif
  x_mach(bnd_type) = max(x_mach(bnd_type), mach)
  n_Ti(bnd_type)   = min(n_Ti(bnd_type), Ti)
  if ( rho .lt. n_rho(bnd_type) ) then
    n_rho(bnd_type) = rho ; n_rho_R(bnd_type) = R ; n_rho_Z(bnd_type) = Z
  endif
  if ( Te .lt. n_Te(bnd_type) ) then
    n_Te(bnd_type) = Te ; n_Te_R(bnd_type) = R ; n_Te_Z(bnd_type) = Z
  endif
  !$omp end critical (floating_diag)
end subroutine floating_diag_add


subroutine floating_diag_report(my_id)

  use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT, EL_CHG
  use phys_module, only: central_density, central_mass
  use mpi_mod

  implicit none
  integer, intent(in) :: my_id

  real*8  :: len(nt), inflow(nt), mom(nt), den(nt), ven(nt), rmin(nt), tmin(nt), mach(nt), timin(nt)
  real*8  :: loc(3,2,nt), loc_g(3,2,nt), v_norm, T_eV
  integer :: it, ierr

  call MPI_ALLREDUCE(s_len, len,    nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_in,  inflow, nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_mom, mom,    nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_den, den,    nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(x_ven, ven,    nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(n_rho, rmin,   nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(n_Te,  tmin,   nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(x_mach, mach,  nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(n_Ti,  timin,  nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

  ! --- Locations: the rank holding each extremum sends its (R,Z), the others send -huge
  loc = -huge(1.d0)
  do it = 1, nt
    if ( x_ven(it) .eq. ven(it)  ) loc(1,:,it) = (/ x_ven_R(it), x_ven_Z(it) /)
    if ( n_rho(it) .eq. rmin(it) ) loc(2,:,it) = (/ n_rho_R(it), n_rho_Z(it) /)
    if ( n_Te(it)  .eq. tmin(it) ) loc(3,:,it) = (/ n_Te_R(it),  n_Te_Z(it)  /)
  enddo
  call MPI_ALLREDUCE(loc, loc_g, 6*nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

  if ( my_id .ne. 0 ) return
  if ( all(len .le. 0.d0) ) return

  v_norm = 1.d0 / sqrt(MU_ZERO * central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT)   ! m/s per unit
  T_eV   = 1.d0 / (EL_CHG * MU_ZERO * central_density * 1.d20)                                ! eV per unit

  write(*,'(A)') ' [floating_u] type  inflow   max vE.n[m/s]  at (R,Z)             mom      max M    min rho    at (R,Z)           min Ti[eV] min Te[eV]  at (R,Z)'
  do it = 1, nt
    if ( len(it) .le. 0.d0 ) cycle
    write(*,'(A,I4,F9.3,ES14.3,2F9.4,ES11.2,F8.2,ES11.2,2F9.4,2ES11.2,2F9.4)')       &
      ' [floating_u] ', it, inflow(it)/len(it), ven(it)*v_norm, loc_g(1,:,it), &
      abs(mom(it))/max(den(it), tiny(1.d0)), mach(it), rmin(it), loc_g(2,:,it), timin(it)*T_eV, tmin(it)*T_eV, loc_g(3,:,it)
  enddo

end subroutine floating_diag_report


end module mod_floating_diag
