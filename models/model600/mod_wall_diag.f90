!> Boundary diagnostics for the drift-compatible sheath boundary (mach1_weak, floating_u): potential, parallel
!! flow, ExB and the sheath fluxes at the wall, written to the log once per matrix construction when wall_diag is
!! set (every wall_diag_every steps). Wall Gauss points of the first toroidal plane; Bn = B_pol.n, b_n = Bn/|B|,
!! vn = Vpar*Bn + vE.n the total outward normal flow, cs|b.n| the Bohm normal speed. Velocities in m/s.
!!
!! [wall]       per boundary type
!!              len     wall length carrying the sheath rows [m]
!!              in      fraction of that length with vn < 0 (net inflow)
!!              sub     fraction with vn < cs*|b.n| (sheath-set excess sink active)
!!              sink    int rho*(cs|b.n| - vn)^+ R dl / int rho*max(vn, cs|b.n|) R dl
!!              mom     |sum Bn*res*dl| / sum |Bn|*cs|b.n|*dl, the moment of the Bohm row
!!              vE.n    largest outward and largest inward ExB normal speed, with (R,Z)
!!              M       largest |Vpar*B|/cs
!!              Phi     wall potential range [V]
!! [wall flow]  per type: min/max of Vpar*Bn, vE.n, vn and cs|b.n|, min |b.n|, max |dPhi/ds| [V/m]
!! [wall min]   per type: min rho [1e20 m^-3], min Ti, min Te [eV], each with (R,Z)
!! [wall pt]    full local state at the wall point with the lowest rho, the lowest Te, the largest |vE.n| and the
!!              largest excess-sink density, over all types and ranks
!! [wall prof]  the full wall profile, one line per wall Gauss point (printed by every rank for its own points),
!!              every wall_diag_profile_every steps (0: never) and automatically on the step where the wall minimum
!!              of rho or Te is non-positive or has fallen below half its value at the previous report
!! [volume]     node values (n=0 part): min rho, min Ti, min Te, max Te with (R,Z) and boundary type (0 = interior)
!! State columns of [wall pt] and [wall prof]:
!!   type R Z rho[1e20] Ti[eV] Te[eV] Phi[V] Vpar*Bn vE.n cs|b.n| vn b.n M dTe/ds[eV/m] dPhi/ds[V/m]
module mod_wall_diag

  implicit none
  private

  public :: wall_diag_reset, wall_diag_add, wall_diag_report, wall_diag_now

  integer, parameter :: nt   = 30                !< max_bnd_types
  integer, parameter :: ns   = 15                !< state vector length
  integer, parameter :: npt  = 4                 !< [wall pt] entries
  integer, parameter :: ncap = 50000             !< wall Gauss points kept per rank for [wall prof]
  real*8, save :: s_len(nt), s_in(nt), s_sub(nt), s_flux(nt), s_sink(nt), s_mom(nt), s_den(nt)
  real*8, save :: x_veo(nt), x_vei(nt), x_mach(nt), x_phi(nt), n_phi(nt)
  real*8, save :: f_min(5,nt), f_max(5,nt)       !< Vpar*Bn, vE.n, vn, cs|b.n|, |b.n| (min) / |dPhi/ds| (max)
  real*8, save :: n_rho(nt), n_Ti(nt), n_Te(nt)
  real*8, save :: l_rho(2,nt), l_Ti(2,nt), l_Te(2,nt), l_veo(2,nt), l_vei(2,nt)
  real*8, save :: p_st(ns,npt), p_key(npt)       !< states and keys (min rho, min Te, max |vE.n|, max sink)
  real*8, save, allocatable :: buf(:,:)
  integer, save :: nbuf = 0
  real*8, save :: prev_rho = huge(1.d0), prev_Te = huge(1.d0)
  logical, save :: active = .false.

contains


!> True on the matrix constructions whose tables are printed.
logical function wall_diag_now()
  use phys_module, only: wall_diag, wall_diag_every, mach1_weak, index_now
  implicit none
  wall_diag_now = wall_diag .and. mach1_weak .and. ( mod(index_now, max(wall_diag_every, 1)) .eq. 0 )
end function wall_diag_now


subroutine wall_diag_reset()
  implicit none
  active = wall_diag_now()
  s_len = 0.d0 ; s_in = 0.d0 ; s_sub = 0.d0 ; s_flux = 0.d0 ; s_sink = 0.d0 ; s_mom = 0.d0 ; s_den = 0.d0
  x_veo = -huge(1.d0) ; x_vei = huge(1.d0) ; x_mach = 0.d0 ; x_phi = -huge(1.d0) ; n_phi = huge(1.d0)
  f_min = huge(1.d0) ; f_max = -huge(1.d0)
  n_rho = huge(1.d0) ; n_Ti = huge(1.d0) ; n_Te = huge(1.d0)
  l_rho = 0.d0 ; l_Ti = 0.d0 ; l_Te = 0.d0 ; l_veo = 0.d0 ; l_vei = 0.d0
  p_st = 0.d0 ; p_key = (/ huge(1.d0), huge(1.d0), -huge(1.d0), -huge(1.d0) /)
  if ( active .and. .not. allocated(buf) ) allocate(buf(ns,ncap))
  nbuf = 0
end subroutine wall_diag_reset


!> One wall Gauss point, JOREK units; w = Gauss weight * dl. u is the potential variable, Te_s and u_s the
!! tangential derivatives per unit length.
subroutine wall_diag_add(bnd_type, w, R, Z, rho, Ti, Te, Vpar, Btot, Bn, b_n, vEn, cs, res, u, Te_s, u_s)
  implicit none
  integer, intent(in) :: bnd_type
  real*8,  intent(in) :: w, R, Z, rho, Ti, Te, Vpar, Btot, Bn, b_n, vEn, cs, res, u, Te_s, u_s
  real*8 :: vn, cb, gam, sink, mach, st(ns)
  if ( .not. active ) return
  if ( bnd_type .lt. 1 .or. bnd_type .gt. nt ) return
  vn   = Bn * Vpar + vEn
  cb   = cs * abs(b_n)
  gam  = rho * max(vn, cb) * R
  sink = rho * max(cb - vn, 0.d0)
  mach = abs(Vpar) * Btot / max(cs, tiny(1.d0))
  st   = (/ dble(bnd_type), R, Z, rho, Ti, Te, u, Bn*Vpar, vEn, cb, vn, b_n, mach, Te_s, u_s /)
  !$omp critical (wall_diag)
  s_len(bnd_type)  = s_len(bnd_type)  + w
  if ( vn .lt. 0.d0 ) s_in(bnd_type)  = s_in(bnd_type)  + w
  if ( vn .lt. cb )   s_sub(bnd_type) = s_sub(bnd_type) + w
  s_flux(bnd_type) = s_flux(bnd_type) + gam * w
  s_sink(bnd_type) = s_sink(bnd_type) + sink * R * w
  s_mom(bnd_type)  = s_mom(bnd_type)  + Bn * res * w
  s_den(bnd_type)  = s_den(bnd_type)  + abs(Bn) * cb * w
  if ( vEn .gt. x_veo(bnd_type) ) then ; x_veo(bnd_type) = vEn ; l_veo(:,bnd_type) = (/ R, Z /) ; endif
  if ( vEn .lt. x_vei(bnd_type) ) then ; x_vei(bnd_type) = vEn ; l_vei(:,bnd_type) = (/ R, Z /) ; endif
  x_mach(bnd_type) = max(x_mach(bnd_type), mach)
  x_phi(bnd_type)  = max(x_phi(bnd_type), u)
  n_phi(bnd_type)  = min(n_phi(bnd_type), u)
  f_min(1:4,bnd_type) = min(f_min(1:4,bnd_type), (/ Bn*Vpar, vEn, vn, cb /))
  f_max(1:4,bnd_type) = max(f_max(1:4,bnd_type), (/ Bn*Vpar, vEn, vn, cb /))
  f_min(5,bnd_type)   = min(f_min(5,bnd_type), abs(b_n))
  f_max(5,bnd_type)   = max(f_max(5,bnd_type), abs(u_s))
  if ( rho .lt. n_rho(bnd_type) ) then ; n_rho(bnd_type) = rho ; l_rho(:,bnd_type) = (/ R, Z /) ; endif
  if ( Ti  .lt. n_Ti(bnd_type)  ) then ; n_Ti(bnd_type)  = Ti  ; l_Ti(:,bnd_type)  = (/ R, Z /) ; endif
  if ( Te  .lt. n_Te(bnd_type)  ) then ; n_Te(bnd_type)  = Te  ; l_Te(:,bnd_type)  = (/ R, Z /) ; endif
  if ( rho       .lt. p_key(1) ) then ; p_key(1) = rho       ; p_st(:,1) = st ; endif
  if ( Te        .lt. p_key(2) ) then ; p_key(2) = Te        ; p_st(:,2) = st ; endif
  if ( abs(vEn)  .gt. p_key(3) ) then ; p_key(3) = abs(vEn)  ; p_st(:,3) = st ; endif
  if ( sink      .gt. p_key(4) ) then ; p_key(4) = sink      ; p_st(:,4) = st ; endif
  if ( nbuf .lt. ncap ) then
    nbuf = nbuf + 1
    buf(:,nbuf) = st
  endif
  !$omp end critical (wall_diag)
end subroutine wall_diag_add


subroutine wall_diag_report(my_id, node_list)

  use constants,          only: MU_ZERO, ATOMIC_MASS_UNIT, EL_CHG
  use phys_module,        only: central_density, central_mass, F0, t_now, index_now, wall_diag_profile_every
  use mod_parameters,     only: var_rho, var_Te, var_Ti, var_T, with_TiTe
  use data_structure,     only: type_node_list
  use mpi_mod

  implicit none
  integer,              intent(in) :: my_id
  type(type_node_list), intent(in) :: node_list

  real*8  :: g_len(nt), g_in(nt), g_sub(nt), g_flux(nt), g_sink(nt), g_mom(nt), g_den(nt)
  real*8  :: g_veo(nt), g_vei(nt), g_mach(nt), g_phi(nt), m_phi(nt), g_rho(nt), g_Ti(nt), g_Te(nt)
  real*8  :: gf_min(5,nt), gf_max(5,nt)
  real*8  :: loc(2,5,nt), loc_g(2,5,nt), pt(ns,npt), pt_g(ns,npt), kmin(2), kmax(2)
  real*8  :: v_norm, T_eV, n_20, phi_V, wall_rho, wall_Te
  real*8  :: nr, nti, nte, xte, val
  real*8  :: lr(3), lti(3), lte(3), lxte(3)
  integer :: it, ierr, inode, k
  logical :: prof
  character(len=8) :: tag(npt) = (/ 'minrho  ', 'minTe   ', 'max|vE| ', 'maxsink ' /)

  if ( .not. active ) return

  call MPI_ALLREDUCE(s_len,  g_len,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_in,   g_in,   nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_sub,  g_sub,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_flux, g_flux, nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_sink, g_sink, nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_mom,  g_mom,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_den,  g_den,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(x_veo,  g_veo,  nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(x_vei,  g_vei,  nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(x_mach, g_mach, nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(x_phi,  g_phi,  nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(n_phi,  m_phi,  nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(f_min,  gf_min, 5*nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(f_max,  gf_max, 5*nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(n_rho,  g_rho,  nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(n_Ti,   g_Ti,   nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(n_Te,   g_Te,   nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

  ! --- Locations and point states: the rank holding each extremum sends its values, the others -huge
  loc = -huge(1.d0)
  do it = 1, nt
    if ( x_veo(it) .eq. g_veo(it) ) loc(:,1,it) = l_veo(:,it)
    if ( x_vei(it) .eq. g_vei(it) ) loc(:,2,it) = l_vei(:,it)
    if ( n_rho(it) .eq. g_rho(it) ) loc(:,3,it) = l_rho(:,it)
    if ( n_Ti(it)  .eq. g_Ti(it)  ) loc(:,4,it) = l_Ti(:,it)
    if ( n_Te(it)  .eq. g_Te(it)  ) loc(:,5,it) = l_Te(:,it)
  enddo
  call MPI_ALLREDUCE(loc, loc_g, 10*nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(p_key(1:2), kmin, 2, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(p_key(3:4), kmax, 2, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  pt = -huge(1.d0)
  do k = 1, 2
    if ( p_key(k)   .eq. kmin(k) ) pt(:,k)   = p_st(:,k)
    if ( p_key(k+2) .eq. kmax(k) ) pt(:,k+2) = p_st(:,k+2)
  enddo
  call MPI_ALLREDUCE(pt, pt_g, ns*npt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

  if ( sum(g_len) .le. 0.d0 ) return

  ! --- Full profile: on request, or when the wall minimum of rho or Te is non-positive or has halved
  wall_rho = minval(g_rho) ; wall_Te = minval(g_Te)
  prof = ( wall_diag_profile_every .gt. 0 .and. mod(index_now, max(wall_diag_profile_every, 1)) .eq. 0 ) &
         .or. wall_rho .le. 0.d0 .or. wall_Te .le. 0.d0                                                  &
         .or. wall_rho .lt. 0.5d0*prev_rho .or. wall_Te .lt. 0.5d0*prev_Te
  prev_rho = wall_rho ; prev_Te = wall_Te

  v_norm = 1.d0 / sqrt(MU_ZERO * central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT)   ! m/s per unit
  T_eV   = 1.d0 / (EL_CHG * MU_ZERO * central_density * 1.d20)                                ! eV per unit
  n_20   = central_density                                                                    ! 1e20 m^-3 per unit
  phi_V  = F0 * v_norm                                                                        ! V per unit of u

  if ( my_id .eq. 0 ) then
    write(*,'(A,I8,A,ES13.5)') ' [wall] step', index_now, '  t_now', t_now
    write(*,'(A)') ' [wall] type  len[m]    in     sub    sink      mom       max vE.n out at (R,Z)             max vE.n in at (R,Z)              max M    Phi min / max [V]'
    do it = 1, nt
      if ( g_len(it) .le. 0.d0 ) cycle
      write(*,'(A,I4,ES10.3,2F7.3,2ES10.2,2(ES12.3,2F8.4),F8.2,2ES11.3)') ' [wall] ', it, g_len(it), &
        g_in(it)/g_len(it), g_sub(it)/g_len(it), g_sink(it)/max(g_flux(it), tiny(1.d0)),         &
        abs(g_mom(it))/max(g_den(it), tiny(1.d0)),                                                &
        g_veo(it)*v_norm, loc_g(:,1,it), g_vei(it)*v_norm, loc_g(:,2,it), g_mach(it),             &
        m_phi(it)*phi_V, g_phi(it)*phi_V
    enddo
    write(*,'(A)') ' [wall flow] type  Vpar*Bn min / max    vE.n min / max       vn min / max         cs|b.n| min / max    min|b.n|   max|dPhi/ds|[V/m]'
    do it = 1, nt
      if ( g_len(it) .le. 0.d0 ) cycle
      write(*,'(A,I4,8ES11.3,2ES11.3)') ' [wall flow] ', it, (gf_min(k,it)*v_norm, gf_max(k,it)*v_norm, k = 1, 4), &
        gf_min(5,it), gf_max(5,it)*phi_V
    enddo
    write(*,'(A)') ' [wall min] type  rho[1e20]  at (R,Z)          Ti[eV]     at (R,Z)          Te[eV]     at (R,Z)'
    do it = 1, nt
      if ( g_len(it) .le. 0.d0 ) cycle
      write(*,'(A,I4,3(ES11.3,2F8.4))') ' [wall min] ', it, g_rho(it)*n_20, loc_g(:,3,it), &
        g_Ti(it)*T_eV, loc_g(:,4,it), g_Te(it)*T_eV, loc_g(:,5,it)
    enddo
    write(*,'(A)') ' [wall pt] what     type  R       Z       rho[1e20]  Ti[eV]     Te[eV]     Phi[V]     Vpar*Bn    vE.n       cs|b.n|    vn         b.n        M          dTe/ds     dPhi/ds'
    do k = 1, npt
      if ( pt_g(1,k) .le. 0.d0 ) cycle
      call write_state(' [wall pt] '//tag(k), pt_g(:,k))
    enddo
    if ( prof ) write(*,'(A,I8,A)') ' [wall prof] step', index_now, &
      '  full wall profile follows (one line per Gauss point, all ranks); columns as [wall pt]'
  endif
  if ( prof ) then
    do k = 1, nbuf
      call write_state(' [wall prof]', buf(:,k))
    enddo
  endif

  if ( my_id .ne. 0 ) return

  ! --- Volume: node values (value DOF, n=0 harmonic)
  nr = huge(1.d0) ; nti = huge(1.d0) ; nte = huge(1.d0) ; xte = -huge(1.d0)
  lr = 0.d0 ; lti = 0.d0 ; lte = 0.d0 ; lxte = 0.d0
  do inode = 1, node_list%n_nodes
    val = node_list%node(inode)%values(1,1,var_rho)
    if ( val .lt. nr ) then ; nr = val ; lr = (/ node_list%node(inode)%x(1,1,1:2), dble(node_list%node(inode)%boundary) /) ; endif
    if ( with_TiTe ) then
      val = node_list%node(inode)%values(1,1,var_Ti)
      if ( val .lt. nti ) then ; nti = val ; lti = (/ node_list%node(inode)%x(1,1,1:2), dble(node_list%node(inode)%boundary) /) ; endif
      val = node_list%node(inode)%values(1,1,var_Te)
    else
      val = 0.5d0 * node_list%node(inode)%values(1,1,var_T)
      if ( val .lt. nti ) then ; nti = val ; lti = (/ node_list%node(inode)%x(1,1,1:2), dble(node_list%node(inode)%boundary) /) ; endif
    endif
    if ( val .lt. nte ) then ; nte = val ; lte  = (/ node_list%node(inode)%x(1,1,1:2), dble(node_list%node(inode)%boundary) /) ; endif
    if ( val .gt. xte ) then ; xte = val ; lxte = (/ node_list%node(inode)%x(1,1,1:2), dble(node_list%node(inode)%boundary) /) ; endif
  enddo
  write(*,'(A)') ' [volume] min rho[1e20] at (R,Z) type | min Ti[eV] at (R,Z) type | min Te[eV] at (R,Z) type | max Te[eV] at (R,Z) type'
  write(*,'(A,4(ES11.3,2F8.4,I4))') ' [volume] ', nr*n_20, lr(1:2), nint(lr(3)), nti*T_eV, lti(1:2), nint(lti(3)), &
    nte*T_eV, lte(1:2), nint(lte(3)), xte*T_eV, lxte(1:2), nint(lxte(3))

contains

  subroutine write_state(label, st)
    character(len=*), intent(in) :: label
    real*8,           intent(in) :: st(ns)
    write(*,'(A,I5,2F8.4,12ES11.3)') label, nint(st(1)), st(2:3), st(4)*n_20, st(5)*T_eV, st(6)*T_eV, &
      st(7)*phi_V, st(8:11)*v_norm, st(12), st(13), st(14)*T_eV, st(15)*phi_V
  end subroutine write_state

end subroutine wall_diag_report


end module mod_wall_diag
