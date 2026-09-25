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

  public :: floating_diag_reset, floating_diag_add, floating_diag_report, sheath_diag_add, wallj_diag_add, wallb_diag_add

  integer, parameter :: nt = 30            !< max_bnd_types
  real*8, save :: s_len(nt), s_in(nt), s_mom(nt), s_den(nt)
  real*8, save :: x_ven(nt), x_ven_R(nt), x_ven_Z(nt)
  real*8, save :: n_rho(nt), n_rho_R(nt), n_rho_Z(nt)
  real*8, save :: n_Te(nt),  n_Te_R(nt),  n_Te_Z(nt)
  real*8, save :: x_mach(nt), n_Ti(nt)
  ! --- sheath current row: wall length carrying it, electron-saturated length, min/max j/j_sat, min/max u,
  ! --- net current into the wall and the saturation current, both as sum j*(B_pol.n)*dl (JOREK units)
  real*8, save :: s_slen(nt), s_esat(nt), j_min(nt), j_max(nt), u_min(nt), u_max(nt), s_inet(nt), s_isat(nt)
  real*8, save :: j_abs(nt), j_abs_R(nt), j_abs_Z(nt)   !< largest |j/j_sat| and where
  real*8, save :: s_over(nt)                              !< wall length with |j/j_sat| > 1
  ! --- [wall J]: implicit boundary currents of the released vorticity row vs the sheath capacity, integrated per
  ! --- type (per unit edge parameter) and the largest local ratio of each with its location
  integer, parameter, public :: nwj = 5                 !< mag, exb, vis, kin, dia (see wallj_diag_add)
  real*8, save :: wj_cap(nt), wj_int(nwj,nt), wj_rmax(nwj,nt), wj_loc(2,nwj,nt)
  integer, parameter, public :: nwb = 8                 !< [wall B]: Ish, gradB, pol, dia, mag, vis, kin, Isat (signed, see wallb_diag_add)
  real*8, save :: wb_int(nwb,nt)

contains


subroutine floating_diag_reset()
  implicit none
  s_len = 0.d0 ; s_in = 0.d0 ; s_mom = 0.d0 ; s_den = 0.d0
  x_ven = -huge(1.d0) ; x_ven_R = 0.d0 ; x_ven_Z = 0.d0
  n_rho =  huge(1.d0) ; n_rho_R = 0.d0 ; n_rho_Z = 0.d0
  n_Te  =  huge(1.d0) ; n_Te_R  = 0.d0 ; n_Te_Z  = 0.d0
  x_mach = 0.d0 ; n_Ti = huge(1.d0)
  s_slen = 0.d0 ; s_esat = 0.d0 ; j_min = huge(1.d0) ; j_max = -huge(1.d0)
  u_min = huge(1.d0) ; u_max = -huge(1.d0) ; s_inet = 0.d0 ; s_isat = 0.d0
  j_abs = 0.d0 ; j_abs_R = 0.d0 ; j_abs_Z = 0.d0 ; s_over = 0.d0
  wj_cap = 0.d0 ; wj_int = 0.d0 ; wj_rmax = 0.d0 ; wj_loc = 0.d0 ; wb_int = 0.d0
end subroutine floating_diag_reset


!> One wall Gauss point, every type: sheath capacity and the implicit vorticity-row wall currents, per unit edge
!! parameter (w = Gauss weight). cur(1:nwj) = magnetisation R^2|p_s|, ExB advection of vorticity rho R^4|w||u_s|,
!! viscous flux, kinetic-energy flux (v_E^2/2)|d_s(R^2 rho)|, diamagnetic 2|tauIC| R^3 |Pi_Z| |du/dn| dl.
subroutine wallj_diag_add(bnd_type, w, R, Z, cap, cur)
  implicit none
  integer, intent(in) :: bnd_type
  real*8,  intent(in) :: w, R, Z, cap, cur(nwj)
  real*8  :: rr
  integer :: i
  if ( bnd_type .lt. 1 .or. bnd_type .gt. nt ) return
  !$omp critical (wallj_diag)
  wj_cap(bnd_type) = wj_cap(bnd_type) + cap * w
  do i = 1, nwj
    wj_int(i,bnd_type) = wj_int(i,bnd_type) + cur(i) * w
    rr = cur(i) / max(cap, tiny(1.d0))
    if ( rr .gt. wj_rmax(i,bnd_type) ) then ; wj_rmax(i,bnd_type) = rr ; wj_loc(:,i,bnd_type) = (/ R, Z /) ; endif
  enddo
  !$omp end critical (wallj_diag)
end subroutine wallj_diag_add


!> One wall Gauss point, every type: the SIGNED wall contents of the u row per unit edge parameter (w = Gauss weight),
!! cur(1:nwb) = sheath current, grad-B pressure flux, lagged polarisation flux, diamagnetic flux, magnetisation flux,
!! viscous flux, kinetic-energy flux, saturation current (normalisation).
subroutine wallb_diag_add(bnd_type, w, cur)
  implicit none
  integer, intent(in) :: bnd_type
  real*8,  intent(in) :: w, cur(nwb)
  if ( bnd_type .lt. 1 .or. bnd_type .gt. nt ) return
  !$omp critical (wallb_diag)
  wb_int(:,bnd_type) = wb_int(:,bnd_type) + cur(:) * w
  !$omp end critical (wallb_diag)
end subroutine wallb_diag_add


!> One wall Gauss point carrying the sheath current row.
subroutine sheath_diag_add(bnd_type, dl, zj, jsat, esat, u, Bn, R, Z, capped)
  implicit none
  integer, intent(in) :: bnd_type
  real*8,  intent(in) :: dl, zj, jsat, u, Bn, R, Z
  logical, intent(in) :: esat, capped   !< electron-saturated (III) / X bound active (I)
  if ( bnd_type .lt. 1 .or. bnd_type .gt. nt ) return
  !$omp critical (sheath_diag)
  s_slen(bnd_type) = s_slen(bnd_type) + dl
  if ( esat .or. capped ) s_esat(bnd_type) = s_esat(bnd_type) + dl
  if ( jsat .ne. 0.d0 ) then
    j_min(bnd_type) = min(j_min(bnd_type), zj/jsat) ; j_max(bnd_type) = max(j_max(bnd_type), zj/jsat)
    if ( abs(zj/jsat) .gt. 1.d0 ) s_over(bnd_type) = s_over(bnd_type) + dl
    if ( abs(zj/jsat) .gt. j_abs(bnd_type) ) then
      j_abs(bnd_type) = abs(zj/jsat) ; j_abs_R(bnd_type) = R ; j_abs_Z(bnd_type) = Z
    endif
  endif
  u_min(bnd_type) = min(u_min(bnd_type), u) ; u_max(bnd_type) = max(u_max(bnd_type), u)
  s_inet(bnd_type) = s_inet(bnd_type) - zj   * Bn * R * dl     ! current into the wall ~ -zj*(B_pol.n)/F0
  s_isat(bnd_type) = s_isat(bnd_type) + abs(jsat * Bn) * R * dl
  !$omp end critical (sheath_diag)
end subroutine sheath_diag_add


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
  use phys_module, only: central_density, central_mass, F0
  use mpi_mod

  implicit none
  integer, intent(in) :: my_id

  real*8  :: len(nt), inflow(nt), mom(nt), den(nt), ven(nt), rmin(nt), tmin(nt), mach(nt), timin(nt)
  real*8  :: slen(nt), esat(nt), jmn(nt), jmx(nt), umn(nt), umx(nt), inet(nt), isat(nt), u_volt
  real*8  :: jab(nt), jloc(2,nt), jloc_g(2,nt), over(nt)
  real*8  :: gcap(nt), gint(nwj,nt), grmx(nwj,nt), wl(2,nwj,nt), wl_g(2,nwj,nt), gb(nwb,nt), bn
  integer :: iw
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
  call MPI_ALLREDUCE(s_slen, slen,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_esat, esat,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(j_min,  jmn,   nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(j_max,  jmx,   nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(u_min,  umn,   nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(u_max,  umx,   nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_inet, inet,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_isat, isat,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(j_abs,  jab,   nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_over, over,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(wj_cap,  gcap,  nt,     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(wj_int,  gint,  nwj*nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(wj_rmax, grmx,  nwj*nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(wb_int,  gb,    nwb*nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  wl = -huge(1.d0)
  do it = 1, nt
    do iw = 1, nwj
      if ( wj_rmax(iw,it) .eq. grmx(iw,it) ) wl(:,iw,it) = wj_loc(:,iw,it)
    enddo
  enddo
  call MPI_ALLREDUCE(wl, wl_g, 2*nwj*nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  jloc = -huge(1.d0)
  do it = 1, nt
    if ( j_abs(it) .eq. jab(it) ) jloc(:,it) = (/ j_abs_R(it), j_abs_Z(it) /)
  enddo
  call MPI_ALLREDUCE(jloc, jloc_g, 2*nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

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
  u_volt = F0 / sqrt(MU_ZERO * central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT)     ! volts per unit u

  write(*,'(A)') ' [floating_u] type  inflow   max vE.n[m/s]  at (R,Z)             mom      max M    min rho    at (R,Z)           min Ti[eV] min Te[eV]  at (R,Z)'
  do it = 1, nt
    if ( len(it) .le. 0.d0 ) cycle
    write(*,'(A,I4,F9.3,ES14.3,2F9.4,ES11.2,F8.2,ES11.2,2F9.4,2ES11.2,2F9.4)')       &
      ' [floating_u] ', it, inflow(it)/len(it), ven(it)*v_norm, loc_g(1,:,it), &
      abs(mom(it))/max(den(it), tiny(1.d0)), mach(it), rmin(it), loc_g(2,:,it), timin(it)*T_eV, tmin(it)*T_eV, loc_g(3,:,it)
  enddo

  ! --- sheath current row, where it is carried
  if ( any(slen .gt. 0.d0) ) then
    write(*,'(A)') ' [sheath_j]   type  e-sat  |j|>jsat   j/jsat min     max    max|j/jsat| at (R,Z)      Phi[V] min      max     Inet/Isat'
    do it = 1, nt
      if ( slen(it) .le. 0.d0 ) cycle
      write(*,'(A,I4,2F9.3,2ES11.2,2F9.4,2F12.3,ES12.3)') ' [sheath_j] ', it, esat(it)/slen(it), over(it)/slen(it), jmn(it), jmx(it), &
        jloc_g(:,it), umn(it)*u_volt, umx(it)*u_volt, inet(it)/max(isat(it), tiny(1.d0))
    enddo
  endif

  ! --- implicit vorticity-row wall currents vs the sheath capacity, integrated per type and the largest local ratios
  if ( any(gcap .gt. 0.d0) ) then
    write(*,'(A)') ' [wall J]   type  mag/Icap  exb/Icap  vis/Icap  kin/Icap  dia/Icap | max mag/jcap at (R,Z)        max exb/jcap at (R,Z)        max vis/jcap at (R,Z)        max kin/jcap at (R,Z)        max dia/jcap at (R,Z)'
    do it = 1, nt
      if ( gcap(it) .le. 0.d0 ) cycle
      write(*,'(A,I4,5ES10.2,A,5(ES10.2,2F8.4,1X))') ' [wall J] ', it, ( gint(iw,it)/gcap(it), iw = 1, nwj ), ' |', &
        ( grmx(iw,it), wl_g(:,iw,it), iw = 1, nwj )
    enddo
    ! --- signed wall balance of the u row, integrated per type, in units of the type's saturation current, in the
    ! --- sign of the row's right-hand side: rest = -(Ish + gradB + pol + dia) is what the boundary cannot see (the
    ! --- half-cell volume parts, essentially the interior parallel current); mag/vis/kin are the cancelled fluxes
    write(*,'(A)') ' [wall B]   type     Ish/Isat   gradB/Isat  pol/Isat(lag) dia/Isat |    rest/Isat | cancelled:   mag/Isat    vis/Isat    kin/Isat'
    do it = 1, nt
      if ( gcap(it) .le. 0.d0 ) cycle
      bn = max( gb(8,it), tiny(1.d0) )
      write(*,'(A,I4,4ES12.3,A,ES12.3,A,3ES12.3)') ' [wall B] ', it, gb(1,it)/bn, gb(2,it)/bn, gb(3,it)/bn, gb(4,it)/bn, ' |', &
        - ( gb(1,it) + gb(2,it) + gb(3,it) + gb(4,it) ) / bn, ' |            ', gb(5,it)/bn, gb(6,it)/bn, gb(7,it)/bn
    enddo
  endif

end subroutine floating_diag_report


end module mod_floating_diag
