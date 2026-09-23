!> Boundary diagnostics for the drift-compatible sheath boundary (mach1_weak, floating_u): potential, parallel
!! flow, ExB and the sheath fluxes at the wall, written to the log once per matrix construction when wall_diag is
!! set (every wall_diag_every steps). Wall Gauss points of the first toroidal plane; Bn = B_pol.n, b_n = Bn/|B|,
!! vn = Vpar*Bn + vE.n the total outward normal flow, cs|b.n| the Bohm normal speed. Velocities in m/s.
!!
!! [floating_u] one line per boundary type: max floating-row residual |u - C_T*Te - C_V*V_wall| [V] at the wall
!!              nodes, largest outward and largest inward ExB normal speed, min rho and min Te at the wall Gauss
!!              points, each with (R,Z); fraction of the wall length with net inflow (vn < 0); share of the wall
!!              particle loss supplied by the sheath-set excess sink
!! [mach1_weak] one line per boundary type: Bohm-row moment |sum Bn*res*dl| / sum |Bn|*cs|b.n|*dl (what the
!!              Galerkin row imposes), pointwise |res| = |Bn*Vpar - target| min/mean/max [m/s] with the location
!!              of the max, max Mach |Vpar*B|/cs with location, max target, max cs|b.n|, max unbounded drift
!!              demand |vE.n|/(cs|b.n|), fractions of the wall with drift compensation and with the bound active,
!!              number of Gauss points
!! [wall prof]  only if wall_diag_profile_every > 0: the full wall profile, one line per Gauss point, every that
!!              many steps and on the step where the wall minimum of rho or Te is non-positive or has halved
module mod_wall_diag

  implicit none
  private

  public :: wall_diag_reset, wall_diag_add, wall_diag_report, wall_diag_now

  integer, parameter :: nt   = 30                !< max_bnd_types
  integer, parameter :: ns   = 15                !< state vector length
  integer, parameter :: npt  = 5                 !< [wall pt] entries
  integer, parameter :: ncap = 50000             !< wall Gauss points kept per rank for [wall prof]
  real*8, save :: s_len(nt), s_in(nt), s_sub(nt), s_flux(nt), s_sink(nt), s_mom(nt), s_den(nt), s_bnd(nt)
  real*8, save :: s_dc(nt), s_npt(nt), r_min(nt), r_sum(nt), r_max(nt), l_res(2,nt), t_max(nt), c_max(nt), d_max(nt)  ! [mach1_weak]
  real*8, save :: x_veo(nt), x_vei(nt), x_mach(nt), x_phi(nt), n_phi(nt)
  real*8, save :: f_min(5,nt), f_max(5,nt)       !< Vpar*Bn, vE.n, vn, cs|b.n|, |b.n| (min) / |dPhi/ds| (max)
  real*8, save :: n_rho(nt), n_Ti(nt), n_Te(nt)
  real*8, save :: l_rho(2,nt), l_Ti(2,nt), l_Te(2,nt), l_veo(2,nt), l_vei(2,nt), l_mach(2,nt)
  real*8, save :: p_st(ns,npt), p_key(npt)       !< states and keys (min rho, min Te, max |vE.n|, max sink, max M)
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
  s_len = 0.d0 ; s_in = 0.d0 ; s_sub = 0.d0 ; s_flux = 0.d0 ; s_sink = 0.d0 ; s_mom = 0.d0 ; s_den = 0.d0 ; s_bnd = 0.d0
  s_dc = 0.d0 ; s_npt = 0.d0 ; r_min = huge(1.d0) ; r_sum = 0.d0 ; r_max = -huge(1.d0) ; l_res = 0.d0
  t_max = -huge(1.d0) ; c_max = -huge(1.d0) ; d_max = -huge(1.d0)
  x_veo = -huge(1.d0) ; x_vei = huge(1.d0) ; x_mach = 0.d0 ; x_phi = -huge(1.d0) ; n_phi = huge(1.d0)
  f_min = huge(1.d0) ; f_max = -huge(1.d0)
  n_rho = huge(1.d0) ; n_Ti = huge(1.d0) ; n_Te = huge(1.d0)
  l_rho = 0.d0 ; l_Ti = 0.d0 ; l_Te = 0.d0 ; l_veo = 0.d0 ; l_vei = 0.d0 ; l_mach = 0.d0
  p_st = 0.d0 ; p_key = (/ huge(1.d0), huge(1.d0), -huge(1.d0), -huge(1.d0), -huge(1.d0) /)
  if ( active .and. .not. allocated(buf) ) allocate(buf(ns,ncap))
  nbuf = 0
end subroutine wall_diag_reset


!> One wall Gauss point, JOREK units; w = Gauss weight * dl. u is the potential variable, Te_s and u_s the
!! tangential derivatives per unit length, vfl the flux floor the sheath rows apply here (cs*|b.n| where the Bohm
!! row compensates the drift, 0 below the grazing cut or with the marginal row).
subroutine wall_diag_add(bnd_type, w, R, Z, rho, Ti, Te, Vpar, Btot, Bn, b_n, vEn, cs, res, u, Te_s, u_s, vfl, bnd, dc, tgt)
  implicit none
  integer, intent(in) :: bnd_type
  real*8,  intent(in) :: w, R, Z, rho, Ti, Te, Vpar, Btot, Bn, b_n, vEn, cs, res, u, Te_s, u_s, vfl
  logical, intent(in) :: bnd   !< the drift bound of the Vpar row is active here
  logical, intent(in) :: dc    !< the row compensates the drift here
  real*8,  intent(in) :: tgt   !< the row target (parallel normal speed)
  real*8 :: vn, cb, gam, sink, mach, st(ns)
  if ( .not. active ) return
  if ( bnd_type .lt. 1 .or. bnd_type .gt. nt ) return
  vn   = Bn * Vpar + vEn
  cb   = cs * abs(b_n)
  gam  = rho * max(vn, vfl) * R
  sink = rho * max(vfl - vn, 0.d0)
  mach = abs(Vpar) * Btot / max(cs, tiny(1.d0))
  st   = (/ dble(bnd_type), R, Z, rho, Ti, Te, u, Bn*Vpar, vEn, cb, vn, b_n, mach, Te_s, u_s /)
  !$omp critical (wall_diag)
  s_len(bnd_type)  = s_len(bnd_type)  + w
  if ( vn .lt. 0.d0 ) s_in(bnd_type)  = s_in(bnd_type)  + w
  if ( vn .lt. vfl )  s_sub(bnd_type) = s_sub(bnd_type) + w
  if ( bnd ) s_bnd(bnd_type) = s_bnd(bnd_type) + w
  if ( dc  ) s_dc(bnd_type)  = s_dc(bnd_type)  + w
  s_npt(bnd_type) = s_npt(bnd_type) + 1.d0
  r_min(bnd_type) = min(r_min(bnd_type), abs(res)) ; r_sum(bnd_type) = r_sum(bnd_type) + abs(res)
  if ( abs(res) .gt. r_max(bnd_type) ) then ; r_max(bnd_type) = abs(res) ; l_res(:,bnd_type) = (/ R, Z /) ; endif
  t_max(bnd_type) = max(t_max(bnd_type), abs(tgt)) ; c_max(bnd_type) = max(c_max(bnd_type), cb)
  d_max(bnd_type) = max(d_max(bnd_type), abs(vEn) / max(cb, tiny(1.d0)))
  s_flux(bnd_type) = s_flux(bnd_type) + gam * w
  s_sink(bnd_type) = s_sink(bnd_type) + sink * R * w
  s_mom(bnd_type)  = s_mom(bnd_type)  + Bn * res * w
  s_den(bnd_type)  = s_den(bnd_type)  + abs(Bn) * cb * w
  if ( vEn .gt. x_veo(bnd_type) ) then ; x_veo(bnd_type) = vEn ; l_veo(:,bnd_type) = (/ R, Z /) ; endif
  if ( vEn .lt. x_vei(bnd_type) ) then ; x_vei(bnd_type) = vEn ; l_vei(:,bnd_type) = (/ R, Z /) ; endif
  if ( mach .gt. x_mach(bnd_type) ) then ; x_mach(bnd_type) = mach ; l_mach(:,bnd_type) = (/ R, Z /) ; endif
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
  if ( mach      .gt. p_key(5) ) then ; p_key(5) = mach      ; p_st(:,5) = st ; endif
  if ( nbuf .lt. ncap ) then
    nbuf = nbuf + 1
    buf(:,nbuf) = st
  endif
  !$omp end critical (wall_diag)
end subroutine wall_diag_add


subroutine wall_diag_report(my_id, node_list)

  use constants,          only: MU_ZERO, ATOMIC_MASS_UNIT, EL_CHG
  use phys_module,        only: central_density, central_mass, F0, t_now, index_now, wall_diag_profile_every, bcs, sheath_V_wall
  use mod_floating_u,     only: floating_u_norm
  use mod_parameters,     only: var_rho, var_Te, var_Ti, var_T, var_u, with_TiTe
  use data_structure,     only: type_node_list
  use mpi_mod

  implicit none
  integer,              intent(in) :: my_id
  type(type_node_list), intent(in) :: node_list

  real*8  :: g_len(nt), g_in(nt), g_sub(nt), g_flux(nt), g_sink(nt), g_mom(nt), g_den(nt), g_bnd(nt)
  real*8  :: g_dc(nt), g_npt(nt), gr_min(nt), gr_sum(nt), gr_max(nt), gt_max(nt), gc_max(nt), gd_max(nt), lres(2,nt), lres_g(2,nt)
  real*8  :: f_res(nt), f_loc(2,nt), f_loc_g(2,nt), fu_a_n, fu_C_T, fu_C_V, fres, ur, sink_share(nt), lm(2,nt)
  integer :: ib, it2
  real*8  :: g_veo(nt), g_vei(nt), g_mach(nt), g_phi(nt), m_phi(nt), g_rho(nt), g_Ti(nt), g_Te(nt)
  real*8  :: gf_min(5,nt), gf_max(5,nt)
  real*8  :: loc(2,5,nt), loc_g(2,5,nt), pt(ns,npt), pt_g(ns,npt), kmin(2), kmax(3)
  real*8  :: v_norm, T_eV, n_20, phi_V, wall_rho, wall_Te
  real*8  :: nr, nti, nte, xte, val
  real*8  :: lr(3), lti(3), lte(3), lxte(3)
  integer :: it, ierr, inode, k
  logical :: prof
  character(len=8) :: tag(npt) = (/ 'minrho  ', 'minTe   ', 'max|vE| ', 'maxsink ', 'maxMach ' /)

  if ( .not. active ) return

  call MPI_ALLREDUCE(s_len,  g_len,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_in,   g_in,   nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_sub,  g_sub,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_flux, g_flux, nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_sink, g_sink, nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_mom,  g_mom,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_den,  g_den,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_bnd,  g_bnd,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_dc,   g_dc,   nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(s_npt,  g_npt,  nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(r_min,  gr_min, nt, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(r_sum,  gr_sum, nt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(r_max,  gr_max, nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(t_max,  gt_max, nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(c_max,  gc_max, nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(d_max,  gd_max, nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  lres = -huge(1.d0)
  do it = 1, nt
    if ( r_max(it) .eq. gr_max(it) ) lres(:,it) = l_res(:,it)
  enddo
  call MPI_ALLREDUCE(lres, lres_g, 2*nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
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
  call MPI_ALLREDUCE(p_key(3:5), kmax, 3, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  pt = -huge(1.d0)
  do k = 1, 2
    if ( p_key(k)   .eq. kmin(k) ) pt(:,k)   = p_st(:,k)
  enddo
  do k = 1, 3
    if ( p_key(k+2) .eq. kmax(k) ) pt(:,k+2) = p_st(:,k+2)
  enddo
  call MPI_ALLREDUCE(pt, pt_g, ns*npt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

  if ( sum(g_len) .le. 0.d0 ) return

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

  ! --- location of max M
  lm = -huge(1.d0)
  do it = 1, nt
    if ( x_mach(it) .eq. g_mach(it) ) lm(:,it) = l_mach(:,it)
  enddo
  call MPI_ALLREDUCE(lm, lres, 2*nt, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  sink_share = 0.d0
  do it = 1, nt
    if ( g_flux(it) .gt. 0.d0 ) sink_share(it) = g_sink(it) / g_flux(it)
  enddo

  ! --- Floating-row residual |u - C_T*Te - C_V*V_wall| at the wall nodes of floating types (n=0 value DOF)
  f_res = -1.d0 ; f_loc = 0.d0
  if ( any(bcs(:)%floating_u) ) then
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
      fres = abs( node_list%node(inode)%values(1,1,var_u) - fu_C_T*ur - fu_C_V*sheath_V_wall )
      if ( fres .gt. f_res(ib) ) then ; f_res(ib) = fres ; f_loc(:,ib) = node_list%node(inode)%x(1,1,1:2) ; endif
    enddo
  endif

  if ( my_id .eq. 0 ) then
    write(*,'(A,I8,A,ES13.5)') ' [floating_u] step', index_now, '  t_now', t_now
    write(*,'(A)') ' [floating_u] type  |u-uf|[V]  at (R,Z)           vE.n out[m/s] at (R,Z)          vE.n in[m/s]  at (R,Z)          min rho[1e20] at (R,Z)         min Te[eV]  at (R,Z)           inflow  sink'
    do it = 1, nt
      if ( g_len(it) .le. 0.d0 ) cycle
      write(*,'(A,I4,ES11.3,2F8.4,4(ES11.3,2F8.4),F8.3,ES10.2)') ' [floating_u] ', it, max(f_res(it),0.d0)*phi_V, f_loc(:,it), &
        g_veo(it)*v_norm, loc_g(:,1,it), g_vei(it)*v_norm, loc_g(:,2,it), g_rho(it)*n_20, loc_g(:,3,it), &
        g_Te(it)*T_eV, loc_g(:,5,it), g_in(it)/g_len(it), sink_share(it)
    enddo
    write(*,'(A)') ' [mach1_weak] type  mom      |res| min / mean / max [m/s]        at (R,Z)         max M    at (R,Z)          max tgt[m/s] max cs|b.n|  drift/cs  comp   bnd    npts'
    do it = 1, nt
      if ( g_len(it) .le. 0.d0 ) cycle
      write(*,'(A,I4,ES10.2,3ES11.3,2F8.4,F8.2,2F8.4,2ES12.3,ES10.2,2F7.3,F7.0)') ' [mach1_weak] ', it, &
        abs(g_mom(it))/max(g_den(it), tiny(1.d0)), gr_min(it)*v_norm, gr_sum(it)/max(g_npt(it),1.d0)*v_norm, gr_max(it)*v_norm, &
        lres_g(:,it), g_mach(it), lres(:,it), gt_max(it)*v_norm, gc_max(it)*v_norm, gd_max(it), g_dc(it)/g_len(it), &
        g_bnd(it)/g_len(it), g_npt(it)
    enddo
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
