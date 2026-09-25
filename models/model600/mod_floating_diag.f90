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

  public :: floating_diag_reset, floating_diag_add, floating_diag_report, sheath_diag_add, wallj_diag_add

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
  real*8, save :: wj_cap(nt), wj_mag(nt), wj_exb(nt), wj_vis(nt), wj_rmag(nt), wj_rexb(nt), wj_rvis(nt), wj_loc(2,3,nt)

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
  wj_cap = 0.d0 ; wj_mag = 0.d0 ; wj_exb = 0.d0 ; wj_vis = 0.d0 ; wj_rmag = 0.d0 ; wj_rexb = 0.d0 ; wj_rvis = 0.d0 ; wj_loc = 0.d0
end subroutine floating_diag_reset


!> One wall Gauss point, every type: sheath capacity and the three implicit vorticity-row wall currents, per unit edge
!! parameter (w = Gauss weight).
subroutine wallj_diag_add(bnd_type, w, R, Z, cap, mag, exb, vis)
  implicit none
  integer, intent(in) :: bnd_type
  real*8,  intent(in) :: w, R, Z, cap, mag, exb, vis
  real*8 :: rm, re, rv
  if ( bnd_type .lt. 1 .or. bnd_type .gt. nt ) return
  rm = mag / max(cap, tiny(1.d0)) ; re = exb / max(cap, tiny(1.d0)) ; rv = vis / max(cap, tiny(1.d0))
  !$omp critical (wallj_diag)
  wj_cap(bnd_type) = wj_cap(bnd_type) + cap * w
  wj_mag(bnd_type) = wj_mag(bnd_type) + mag * w
  wj_exb(bnd_type) = wj_exb(bnd_type) + exb * w
  wj_vis(bnd_type) = wj_vis(bnd_type) + vis * w
  if ( rm .gt. wj_rmag(bnd_type) ) then ; wj_rmag(bnd_type) = rm ; wj_loc(:,1,bnd_type) = (/ R, Z /) ; endif
  if ( re .gt. wj_rexb(bnd_type) ) then ; wj_rexb(bnd_type) = re ; wj_loc(:,2,bnd_type) = (/ R, Z /) ; endif
  if ( rv .gt. wj_rvis(bnd_type) ) then ; wj_rvis(bnd_type) = rv ; wj_loc(:,3,bnd_type) = (/ R, Z /) ; endif
  !$omp end critical (wallj_diag)
end subroutine wallj_diag_add


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
  real*8  :: gcap(nt), gmag(nt), gexb(nt), gvis(nt), grm(nt), gre(nt), grv(nt), wl(2,3,nt), wl_g(2,3,nt)
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
  call MPI_ALLREDUCE(wj_cap, gcap,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(wj_mag, gmag,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(wj_exb, gexb,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(wj_vis, gvis,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(wj_rmag, grm,  nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(wj_rexb, gre,  nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(wj_rvis, grv,  nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  wl = -huge(1.d0)
  do it = 1, nt
    if ( wj_rmag(it) .eq. grm(it) ) wl(:,1,it) = wj_loc(:,1,it)
    if ( wj_rexb(it) .eq. gre(it) ) wl(:,2,it) = wj_loc(:,2,it)
    if ( wj_rvis(it) .eq. grv(it) ) wl(:,3,it) = wj_loc(:,3,it)
  enddo
  call MPI_ALLREDUCE(wl, wl_g, 6*nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
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
    write(*,'(A)') ' [wall J]   type  mag/Icap  exb/Icap  vis/Icap | max mag/jcap at (R,Z)        max exb/jcap at (R,Z)        max vis/jcap at (R,Z)'
    do it = 1, nt
      if ( gcap(it) .le. 0.d0 ) cycle
      write(*,'(A,I4,3ES10.2,A,3(ES10.2,2F8.4,1X))') ' [wall J] ', it, gmag(it)/gcap(it), gexb(it)/gcap(it), gvis(it)/gcap(it), ' |', &
        grm(it), wl_g(:,1,it), gre(it), wl_g(:,2,it), grv(it), wl_g(:,3,it)
    enddo
  endif

end subroutine floating_diag_report


end module mod_floating_diag
