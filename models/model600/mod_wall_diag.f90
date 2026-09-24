!> Wall diagnostics for the floating-potential boundary (bcs%floating_u): potential, parallel flow, ExB drift
!! and the plasma state at the wall, written to the log once per matrix construction when wall_diag is set
!! (every wall_diag_every steps). Sampled at the wall Gauss points of the open-boundary integral
!! (mod_boundary_matrix_open), first toroidal plane; Bn = B_pol.n, b_n = Bn/|B|, vn = Vpar*Bn + vE.n the total
!! outward normal flow, cs|b.n| the Bohm normal speed. Velocities in m/s, one line per boundary type.
!!
!! [floating_u] max floating-row residual |u - C_T*max(Te,T_min) - C_V*V_wall| [V] at the wall nodes, largest outward and
!!              largest inward ExB normal speed, min rho and min Te at the wall Gauss points, each with (R,Z);
!!              fraction of the wall length with net inflow (vn < 0) and with the ExB normal speed above the
!!              Bohm normal speed (|vE.n| > cs|b.n|); sink/Bohm: excess sink of the sheath-set wall flux integrated
!!              over the type, relative to its Bohm floor (0 without the flux)
!! [mach1]      residual of the drift-free Mach-1 condition in normal-speed units, |res| = |Bn*Vpar - cs|b.n||
!!              min/mean/max with the location of the max (with mach1_omit_drift this is what the nodal row
!!              imposes; with the drift term the row's target differs by its ExB term), max Mach |Vpar*B|/cs
!!              with location, max cs|b.n|, max drift demand |vE.n|/(cs|b.n|), number of Gauss points
!! [sheath_j]   per type carrying the sheath current row: fraction of its length at electron saturation, fraction
!!              with |j/j_sat| > 1, min/max j/j_sat (negative = electron current, +1 = ion saturation), largest
!!              |j/j_sat| with (R,Z), wall potential min/max [V], net current into the wall over the saturation current
!! [wall prof]  only if wall_diag_profile_every > 0: the full wall profile, one line per Gauss point, every that
!!              many steps and on the step where the wall minimum of rho or Te is non-positive or has halved
module mod_wall_diag

  implicit none
  private

  public :: wall_diag_reset, wall_diag_add, wall_diag_sheath_add, wall_diag_report, wall_diag_now

  integer, parameter :: nt   = 30                !< max_bnd_types
  integer, parameter :: ns   = 15                !< state vector length
  integer, parameter :: ncap = 50000             !< wall Gauss points kept per rank for [wall prof]
  real*8, save :: s_len(nt), s_in(nt), s_exb(nt), s_npt(nt), s_sink(nt), s_bohm(nt)
  real*8, save :: r_min(nt), r_sum(nt), r_max(nt), l_res(2,nt), c_max(nt), d_max(nt)
  real*8, save :: x_veo(nt), x_vei(nt), x_mach(nt)
  real*8, save :: n_rho(nt), n_Te(nt)
  real*8, save :: l_rho(2,nt), l_Te(2,nt), l_veo(2,nt), l_vei(2,nt), l_mach(2,nt)
  ! --- sheath current row: wall length carrying it, electron-saturated length, length with |j/j_sat| > 1, min/max
  ! --- j/j_sat, largest |j/j_sat| and where, min/max u, net current into the wall and the saturation current
  real*8, save :: j_len(nt), j_esat(nt), j_over(nt), j_min(nt), j_max(nt), j_abs(nt), l_jabs(2,nt), ju_min(nt), ju_max(nt), j_inet(nt), j_isat(nt)
  real*8, save, allocatable :: buf(:,:)
  integer, save :: nbuf = 0
  real*8, save :: prev_rho = huge(1.d0), prev_Te = huge(1.d0)
  logical, save :: active = .false.

contains


!> True on the matrix constructions whose tables are printed.
logical function wall_diag_now()
  use phys_module, only: wall_diag, wall_diag_every, bcs, index_now
  implicit none
  wall_diag_now = wall_diag .and. ( any(bcs(:)%floating_u) .or. any(bcs(:)%sheath_j) ) .and. ( mod(index_now, max(wall_diag_every, 1)) .eq. 0 )
end function wall_diag_now


subroutine wall_diag_reset()
  implicit none
  active = wall_diag_now()
  s_len = 0.d0 ; s_in = 0.d0 ; s_exb = 0.d0 ; s_npt = 0.d0 ; s_sink = 0.d0 ; s_bohm = 0.d0
  r_min = huge(1.d0) ; r_sum = 0.d0 ; r_max = -huge(1.d0) ; l_res = 0.d0
  c_max = -huge(1.d0) ; d_max = -huge(1.d0)
  x_veo = -huge(1.d0) ; x_vei = huge(1.d0) ; x_mach = 0.d0
  n_rho = huge(1.d0) ; n_Te = huge(1.d0)
  l_rho = 0.d0 ; l_Te = 0.d0 ; l_veo = 0.d0 ; l_vei = 0.d0 ; l_mach = 0.d0
  j_len = 0.d0 ; j_esat = 0.d0 ; j_over = 0.d0 ; j_min = huge(1.d0) ; j_max = -huge(1.d0) ; j_abs = 0.d0 ; l_jabs = 0.d0
  ju_min = huge(1.d0) ; ju_max = -huge(1.d0) ; j_inet = 0.d0 ; j_isat = 0.d0
  if ( active .and. .not. allocated(buf) ) allocate(buf(ns,ncap))
  nbuf = 0
end subroutine wall_diag_reset


!> One wall Gauss point, JOREK units; w = Gauss weight * dl. u is the potential variable, Te_s and u_s the
!! tangential derivatives per unit length.
!! sink = rho*max(v_fl - vn, 0)*R*dl the excess sink of the sheath-set flux at this point, bohm = rho*v_fl*R*dl its floor.
subroutine wall_diag_add(bnd_type, w, R, Z, rho, Ti, Te, Vpar, Btot, Bn, b_n, vEn, cs, u, Te_s, u_s, sink, bohm)
  implicit none
  integer, intent(in) :: bnd_type
  real*8,  intent(in) :: w, R, Z, rho, Ti, Te, Vpar, Btot, Bn, b_n, vEn, cs, u, Te_s, u_s, sink, bohm
  real*8 :: vn, cb, res, mach, st(ns)
  if ( .not. active ) return
  if ( bnd_type .lt. 1 .or. bnd_type .gt. nt ) return
  vn   = Bn * Vpar + vEn
  cb   = cs * abs(b_n)
  res  = Bn * Vpar - cb
  mach = abs(Vpar) * Btot / max(cs, tiny(1.d0))
  st   = (/ dble(bnd_type), R, Z, rho, Ti, Te, u, Bn*Vpar, vEn, cb, vn, b_n, mach, Te_s, u_s /)
  !$omp critical (wall_diag)
  s_len(bnd_type) = s_len(bnd_type) + w
  if ( vn .lt. 0.d0 )      s_in(bnd_type)  = s_in(bnd_type)  + w
  if ( abs(vEn) .gt. cb )  s_exb(bnd_type) = s_exb(bnd_type) + w
  s_npt(bnd_type) = s_npt(bnd_type) + 1.d0
  s_sink(bnd_type) = s_sink(bnd_type) + sink * w ; s_bohm(bnd_type) = s_bohm(bnd_type) + bohm * w
  r_min(bnd_type) = min(r_min(bnd_type), abs(res)) ; r_sum(bnd_type) = r_sum(bnd_type) + abs(res)
  if ( abs(res) .gt. r_max(bnd_type) ) then ; r_max(bnd_type) = abs(res) ; l_res(:,bnd_type) = (/ R, Z /) ; endif
  c_max(bnd_type) = max(c_max(bnd_type), cb)
  d_max(bnd_type) = max(d_max(bnd_type), abs(vEn) / max(cb, tiny(1.d0)))
  if ( vEn .gt. x_veo(bnd_type) ) then ; x_veo(bnd_type) = vEn ; l_veo(:,bnd_type) = (/ R, Z /) ; endif
  if ( vEn .lt. x_vei(bnd_type) ) then ; x_vei(bnd_type) = vEn ; l_vei(:,bnd_type) = (/ R, Z /) ; endif
  if ( mach .gt. x_mach(bnd_type) ) then ; x_mach(bnd_type) = mach ; l_mach(:,bnd_type) = (/ R, Z /) ; endif
  if ( rho .lt. n_rho(bnd_type) ) then ; n_rho(bnd_type) = rho ; l_rho(:,bnd_type) = (/ R, Z /) ; endif
  if ( Te  .lt. n_Te(bnd_type)  ) then ; n_Te(bnd_type)  = Te  ; l_Te(:,bnd_type)  = (/ R, Z /) ; endif
  if ( nbuf .lt. ncap ) then
    nbuf = nbuf + 1
    buf(:,nbuf) = st
  endif
  !$omp end critical (wall_diag)
end subroutine wall_diag_add


!> One wall Gauss point carrying the sheath current row: zj and j_sat (JOREK units), Bn = B_pol.n, u, capped =
!! electron saturation active (x >= Lambda).
subroutine wall_diag_sheath_add(bnd_type, w, R, Z, zj, jsat, Bn, u, capped)
  implicit none
  integer, intent(in) :: bnd_type
  real*8,  intent(in) :: w, R, Z, zj, jsat, Bn, u
  logical, intent(in) :: capped
  if ( .not. active ) return
  if ( bnd_type .lt. 1 .or. bnd_type .gt. nt ) return
  !$omp critical (wall_diag_sheath)
  j_len(bnd_type) = j_len(bnd_type) + w
  if ( capped ) j_esat(bnd_type) = j_esat(bnd_type) + w
  if ( jsat .ne. 0.d0 ) then
    j_min(bnd_type) = min(j_min(bnd_type), zj/jsat) ; j_max(bnd_type) = max(j_max(bnd_type), zj/jsat)
    if ( abs(zj/jsat) .gt. 1.d0 ) j_over(bnd_type) = j_over(bnd_type) + w
    if ( abs(zj/jsat) .gt. j_abs(bnd_type) ) then ; j_abs(bnd_type) = abs(zj/jsat) ; l_jabs(:,bnd_type) = (/ R, Z /) ; endif
  endif
  ju_min(bnd_type) = min(ju_min(bnd_type), u) ; ju_max(bnd_type) = max(ju_max(bnd_type), u)
  j_inet(bnd_type) = j_inet(bnd_type) - zj * Bn * R * w          ! current into the wall ~ -zj*(B_pol.n)/F0
  j_isat(bnd_type) = j_isat(bnd_type) + abs(jsat * Bn) * R * w
  !$omp end critical (wall_diag_sheath)
end subroutine wall_diag_sheath_add


subroutine wall_diag_report(my_id, node_list)

  use constants,          only: MU_ZERO, ATOMIC_MASS_UNIT, EL_CHG
  use phys_module,        only: central_density, central_mass, F0, t_now, index_now, wall_diag_profile_every, bcs, sheath_V_wall, T_min
  use mod_floating_u,     only: floating_u_norm
  use mod_parameters,     only: var_rho, var_Te, var_Ti, var_T, var_u, with_TiTe
  use data_structure,     only: type_node_list
  use mpi_mod

  implicit none
  integer,              intent(in) :: my_id
  type(type_node_list), intent(in) :: node_list

  real*8  :: g_len(nt), g_in(nt), g_exb(nt), g_npt(nt), g_sink(nt), g_bohm(nt), sink_bohm(nt), gr_min(nt), gr_sum(nt), gr_max(nt), gc_max(nt), gd_max(nt)
  real*8  :: g_veo(nt), g_vei(nt), g_mach(nt), g_rho(nt), g_Te(nt)
  real*8  :: loc(2,5,nt), loc_g(2,5,nt), lm(2,nt), lm_g(2,nt)
  real*8  :: f_res(nt), f_loc(2,nt), fu_a_n, fu_C_T, fu_C_V, fres, ur
  real*8  :: v_norm, T_eV, n_20, phi_V, wall_rho, wall_Te
  real*8  :: gj_len(nt), gj_esat(nt), gj_over(nt), gj_min(nt), gj_max(nt), gj_abs(nt), gju_min(nt), gju_max(nt), gj_inet(nt), gj_isat(nt)
  real*8  :: lj(2,nt), lj_g(2,nt)
  integer :: it, ib, ierr, inode, k
  logical :: prof

  if ( .not. active ) return

  call MPI_ALLREDUCE(s_len,  g_len,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_in,   g_in,   nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_exb,  g_exb,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_npt,  g_npt,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_sink, g_sink, nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_bohm, g_bohm, nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(r_min,  gr_min, nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(r_sum,  gr_sum, nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(r_max,  gr_max, nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(c_max,  gc_max, nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(d_max,  gd_max, nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(x_veo,  g_veo,  nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(x_vei,  g_vei,  nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(x_mach, g_mach, nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(n_rho,  g_rho,  nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(n_Te,   g_Te,   nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

  ! --- Locations: the rank holding each extremum sends its coordinates, the others -huge
  loc = -huge(1.d0) ; lm = -huge(1.d0)
  do it = 1, nt
    if ( x_veo(it)  .eq. g_veo(it)  ) loc(:,1,it) = l_veo(:,it)
    if ( x_vei(it)  .eq. g_vei(it)  ) loc(:,2,it) = l_vei(:,it)
    if ( n_rho(it)  .eq. g_rho(it)  ) loc(:,3,it) = l_rho(:,it)
    if ( n_Te(it)   .eq. g_Te(it)   ) loc(:,4,it) = l_Te(:,it)
    if ( r_max(it)  .eq. gr_max(it) ) loc(:,5,it) = l_res(:,it)
    if ( x_mach(it) .eq. g_mach(it) ) lm(:,it)    = l_mach(:,it)
  enddo
  call MPI_ALLREDUCE(loc, loc_g, 10*nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(lm,  lm_g,   2*nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(j_len,  gj_len,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(j_esat, gj_esat, nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(j_over, gj_over, nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(j_min,  gj_min,  nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(j_max,  gj_max,  nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(j_abs,  gj_abs,  nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(ju_min, gju_min, nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(ju_max, gju_max, nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(j_inet, gj_inet, nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(j_isat, gj_isat, nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  lj = -huge(1.d0)
  do it = 1, nt
    if ( j_abs(it) .eq. gj_abs(it) ) lj(:,it) = l_jabs(:,it)
  enddo
  call MPI_ALLREDUCE(lj, lj_g, 2*nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

  if ( sum(g_len) .le. 0.d0 ) return
  sink_bohm = 0.d0
  do it = 1, nt
    if ( g_bohm(it) .gt. 0.d0 ) sink_bohm(it) = g_sink(it) / g_bohm(it)
  enddo

  ! --- Full profile (only if wall_diag_profile_every > 0): every that many steps, and on the step where the wall
  ! --- minimum of rho or Te is non-positive or has halved
  wall_rho = minval(g_rho) ; wall_Te = minval(g_Te)
  prof = .false.
  if ( wall_diag_profile_every .gt. 0 ) prof = ( mod(index_now, wall_diag_profile_every) .eq. 0 )  &
         .or. wall_rho .le. 0.d0 .or. wall_Te .le. 0.d0                                          &
         .or. wall_rho .lt. 0.5d0*prev_rho .or. wall_Te .lt. 0.5d0*prev_Te
  prev_rho = wall_rho ; prev_Te = wall_Te

  v_norm = 1.d0 / sqrt(MU_ZERO * central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT)   ! m/s per unit
  T_eV   = 1.d0 / (EL_CHG * MU_ZERO * central_density * 1.d20)                                ! eV per unit
  n_20   = central_density                                                                    ! 1e20 m^-3 per unit
  phi_V  = F0 * v_norm                                                                        ! V per unit of u

  ! --- Floating-row residual |u - C_T*max(Te,T_min) - C_V*V_wall| at the wall nodes of floating types (n=0
  ! --- value DOF), with the row's own temperature floor
  f_res = -1.d0 ; f_loc = 0.d0
  call floating_u_norm(fu_a_n, fu_C_T, fu_C_V)
  do inode = 1, node_list%n_nodes
    ib = node_list%node(inode)%boundary
    if ( ib .lt. 1 .or. ib .gt. nt ) cycle
    if ( .not. bcs(ib)%floating_u ) cycle
    if ( with_TiTe ) then
      ur = node_list%node(inode)%values(1,1,var_Te)
    else
      ur = node_list%node(inode)%values(1,1,var_T)
    endif
    fres = abs( node_list%node(inode)%values(1,1,var_u) - fu_C_T*max(ur, T_min) - fu_C_V*sheath_V_wall )
    if ( fres .gt. f_res(ib) ) then ; f_res(ib) = fres ; f_loc(:,ib) = node_list%node(inode)%x(1,1,1:2) ; endif
  enddo

  if ( my_id .eq. 0 ) then
    write(*,'(A,I8,A,ES13.5)') ' [floating_u] step', index_now, '  t_now', t_now
    write(*,'(A)') ' [floating_u] type  |u-uf|[V]  at (R,Z)           vE.n out[m/s] at (R,Z)          vE.n in[m/s]  at (R,Z)          min rho[1e20] at (R,Z)         min Te[eV]  at (R,Z)           inflow  exb>cs  sink/Bohm'
    do it = 1, nt
      if ( g_len(it) .le. 0.d0 ) cycle
      write(*,'(A,I4,ES11.3,2F8.4,4(ES11.3,2F8.4),2F8.3,ES10.2)') ' [floating_u] ', it, max(f_res(it),0.d0)*phi_V, f_loc(:,it), &
        g_veo(it)*v_norm, loc_g(:,1,it), g_vei(it)*v_norm, loc_g(:,2,it), g_rho(it)*n_20, loc_g(:,3,it), &
        g_Te(it)*T_eV, loc_g(:,4,it), g_in(it)/g_len(it), g_exb(it)/g_len(it), sink_bohm(it)
    enddo
    write(*,'(A)') ' [mach1] type  |res| min / mean / max [m/s]        at (R,Z)         max M    at (R,Z)          max cs|b.n|  drift/cs   npts'
    do it = 1, nt
      if ( g_len(it) .le. 0.d0 ) cycle
      write(*,'(A,I4,3ES11.3,2F8.4,F8.2,2F8.4,2ES12.3,F7.0)') ' [mach1] ', it, &
        gr_min(it)*v_norm, gr_sum(it)/max(g_npt(it),1.d0)*v_norm, gr_max(it)*v_norm, loc_g(:,5,it), &
        g_mach(it), lm_g(:,it), gc_max(it)*v_norm, gd_max(it), g_npt(it)
    enddo
    if ( any(gj_len .gt. 0.d0) ) then
      write(*,'(A)') ' [sheath_j] type  e-sat  |j|>jsat   j/jsat min     max    max|j/jsat| at (R,Z)      Phi[V] min      max     Inet/Isat'
      do it = 1, nt
        if ( gj_len(it) .le. 0.d0 ) cycle
        write(*,'(A,I4,2F9.3,2ES11.2,2F9.4,2F12.3,ES12.3)') ' [sheath_j] ', it, gj_esat(it)/gj_len(it), gj_over(it)/gj_len(it), &
          gj_min(it), gj_max(it), lj_g(:,it), gju_min(it)*phi_V, gju_max(it)*phi_V, gj_inet(it)/max(gj_isat(it), tiny(1.d0))
      enddo
    endif
    if ( prof ) write(*,'(A,I8,A)') ' [wall prof] step', index_now, &
      '  full wall profile follows (one line per Gauss point, all ranks): type R Z rho[1e20] Ti Te[eV] Phi[V] Vpar*Bn vE.n cs|b.n| vn[m/s] b.n M dTe/ds[eV/m] dPhi/ds[V/m]'
  endif
  if ( prof ) then
    do k = 1, nbuf
      call write_state(' [wall prof]', buf(:,k))
    enddo
  endif

contains

  subroutine write_state(label, st)
    character(len=*), intent(in) :: label
    real*8,           intent(in) :: st(ns)
    write(*,'(A,I5,2F8.4,12ES11.3)') label, nint(st(1)), st(2:3), st(4)*n_20, st(5)*T_eV, st(6)*T_eV, &
      st(7)*phi_V, st(8:11)*v_norm, st(12), st(13), st(14)*T_eV, st(15)*phi_V
  end subroutine write_state

end subroutine wall_diag_report


end module mod_wall_diag
