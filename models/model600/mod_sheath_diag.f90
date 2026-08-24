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
  public :: sheath_store_psi0, sheath_psi0
  public :: sheath_init_potential

  !> psi's degrees of freedom at t_start, per node. The wall relaxation needs the DEVIATION of
  !! dpsi/dn from its equilibrium value; driving it with the raw value would make the flux drift
  !! even in a quiet plasma. Indexed (dof, node).
  real*8, allocatable, save :: sheath_psi0(:,:)

  real*8,  save :: sd_I_sheath(max_bnd_types) = 0.d0  !< current through the wall, sheath value [A]
  real*8,  save :: sd_I_amp(max_bnd_types)    = 0.d0  !< the same with the interior Ampere current [A]
  real*8,  save :: sd_area(max_bnd_types)     = 0.d0  !< wetted area [m^2], for averages
  real*8,  save :: sd_phi_min                 =  1.d30
  real*8,  save :: sd_phi_max                 = -1.d30
  real*8,  save :: sd_phi_sum                 = 0.d0  !< area weighted, for the mean
  !> Same, but only where the obliqueness gate leaves the sheath term ACTIVE. Where the gate has
  !! removed it, u has no boundary condition at all (dirichlet%u is .false. on these types), so a
  !! runaway there means an unconstrained null space rather than a sheath the plasma is overdriving.
  !! Comparing the two maxima separates those two completely different failures.
  real*8,  save :: sd_phi_max_act             = -1.d30
  real*8,  save :: sd_gate_off_area           = 0.d0
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
  sd_phi_max_act = -1.d30
  sd_gate_off_area = 0.d0
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
subroutine sheath_diag_add(bnd_type, zj_sh, zj0, zj_sat, x_lim, u0, Te0, Bdotn, dS, gate)

  use constants,     only: MU_ZERO
  use phys_module,   only: F0, sheath_X_min, sheath_smooth_dX
  use mod_sheath_bc, only: sheath_norm, sheath_V_wall_at

  implicit none
  integer, intent(in) :: bnd_type
  real*8,  intent(in) :: zj_sh, zj0, zj_sat, x_lim, u0, Te0, Bdotn, dS
  real*8,  intent(in), optional :: gate   !< obliqueness weight; the ratio below is only
                                          !< meaningful where the term is actually active
  real*8,  intent(in), optional :: V_wall_loc !< local wall bias, so e*Phi/kTe is measured against
                                          !< the SAME reference the constraint imposes. Omitting it
                                          !< silently reports against the global sheath_V_wall,
                                          !< which differs wherever sheath_V_wall_asym /= 0.

  real*8 :: jn_sheath, jn_amp, phi_over_te, ratio, a_n, c_sat, vw

  if ( bnd_type .lt. 1 .or. bnd_type .gt. max_bnd_types ) return

  ! --- J.n = J_par*(b.n) = zj*(B.n)/(F0*mu0)   [A/m^2]
  jn_sheath = zj_sh * Bdotn / (F0 * MU_ZERO)
  jn_amp    = zj0   * Bdotn / (F0 * MU_ZERO)

  if ( present(V_wall_loc) ) then
    call sheath_norm(a_n, c_sat, vw, V_wall_loc)
  else
    call sheath_norm(a_n, c_sat, vw)
  endif
  phi_over_te = ( 0.5d0*a_n*u0 - vw ) / max(Te0, 1.d-14)     ! e*Phi/(k*Te)

  ! --- Solvability ratio. The characteristic can only deliver f in (-(exp(-X_min)-1), 1], so
  ! --- abs(zj0/zj_sat) > 1 means NO u satisfies it at this point and u is driven without bound.
  ! --- Only report it where the obliqueness gate leaves the term active: at a gated-off point
  ! --- the ratio diverges harmlessly because the term is not there.
  ratio = 0.d0
  if ( abs(zj_sat) .gt. 0.d0 ) ratio = abs(zj0 / zj_sat)
  if ( present(gate) ) then
    if ( gate .lt. 0.5d0 ) ratio = 0.d0
  endif

  !$omp critical (sheath_diag_accumulate)
  sd_I_sheath(bnd_type)  = sd_I_sheath(bnd_type) + jn_sheath * dS
  sd_I_amp(bnd_type)     = sd_I_amp(bnd_type)    + jn_amp    * dS
  sd_area(bnd_type)      = sd_area(bnd_type)     + dS
  sd_phi_sum             = sd_phi_sum            + phi_over_te * dS
  sd_phi_min             = min(sd_phi_min, phi_over_te)
  sd_phi_max             = max(sd_phi_max, phi_over_te)
  if ( present(gate) ) then
    if ( gate .ge. 0.5d0 ) then
      sd_phi_max_act = max(sd_phi_max_act, phi_over_te)
    else
      sd_gate_off_area = sd_gate_off_area + dS
    endif
  else
    sd_phi_max_act = max(sd_phi_max_act, phi_over_te)
  endif
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
  real*8  :: loc_sum(ns+1), glo_sum(ns+1), loc_max(3), glo_max(3), loc_min(1), glo_min(1)
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
  loc_max(3) = sd_phi_max_act
  loc_sum(ns+1) = sd_gate_off_area
  loc_min(1) = sd_phi_min

  call MPI_Reduce(loc_sum, glo_sum, ns+1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(loc_max, glo_max,  3, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
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

  ! --- The two maxima separate the two failures: if the ACTIVE max settles while the global one
  ! --- runs away, the runaway is at gated-off points where u has no boundary condition at all.
  write(*,'(A,f8.2,A,f5.1,A)')                                                     &
    '         ePhi/kTe max where the sheath is ACTIVE=', glo_max(3),               &
    '   gated-off area ', 1.d2*glo_sum(ns+1)/max(area_tot,1.d-30), ' %'

  do i = 1, max_bnd_types
    if ( glo_sum(i2+i) .le. 0.d0 ) cycle
    write(*,'(A,i3,A,es11.3,A,es11.3,A,es10.3,A)')                                &
      '         bnd type', i, ': I_sheath=', glo_sum(i0+i), ' A  I_Ampere=',      &
      glo_sum(i1+i), ' A  area=', glo_sum(i2+i), ' m^2'
  enddo

end subroutine sheath_diag_report


!> Initialise u on the natural%u boundary types to the floating potential, Lambda*Te/e.
!!
!! Without this the run starts at u = 0, i.e. Phi = 0, which is electron saturation: the sheath
!! immediately demands the full electron saturation current everywhere and the boundary condition
!! has to travel ~Lambda*Te away from where the plasma is, through an exponential, in one implicit
!! step. Starting at X = 0 puts the state where the linearisation is valid.
!!
!! Only the axisymmetric component is touched. Lambda is treated as locally constant so that the
!! derivative degrees of freedom can carry the same relation, keeping grad(u) consistent with
!! grad(Te) rather than leaving it at the old field's value.
subroutine sheath_init_potential(node_list, my_id)

  use mod_parameters
  use data_structure
  use phys_module,   only: bcs
  use mod_sheath_bc, only: sheath_norm, sheath_V_wall_at, sheath_get_lambda
  use mpi_mod

  implicit none

  type (type_node_list), intent(inout) :: node_list
  integer,               intent(in)    :: my_id

  integer :: i, id, ib, ierr
  real*8  :: a_n, c_sat, vw, Ti0, Te0, T0, lam, dlam_dTi, dlam_dTe, cfac
  real*8  :: n_loc(1), n_glo(1)

  n_loc(1) = 0.d0

  do i = 1, node_list%n_nodes

    ib = node_list%node(i)%boundary
    if ( ib .lt. 1 .or. ib .gt. max_bnd_types ) cycle
    ! --- both sheath routes want u to start at its own fixed point rather than at 0, which is
    ! --- deep electron saturation (X = -Lambda) and demands ~19*j_sat of electron current
    if ( .not. (bcs(ib)%natural%u .or. bcs(ib)%sheath_zj) ) cycle

    if ( with_TiTe ) then
      Ti0 = node_list%node(i)%values(1,1,var_Ti)
      Te0 = node_list%node(i)%values(1,1,var_Te)
    else
      T0  = node_list%node(i)%values(1,1,var_T)
      Ti0 = 0.5d0 * T0
      Te0 = 0.5d0 * T0
    endif
    if ( Te0 .le. 0.d0 ) cycle

    ! --- per node: with a differentially biased wall the floating potential is a function of R
    call sheath_norm(a_n, c_sat, vw, sheath_V_wall_at(node_list%node(i)%x(1,1,1)))

    call sheath_get_lambda(Ti0, Te0, lam, dlam_dTi, dlam_dTe)

    ! --- zero net current (X = 0) sits at e*Phi/(k*Te) = Lambda, i.e. a_n*u/2 - vw = Lambda*Te
    cfac = 2.d0 * lam / a_n
    node_list%node(i)%values(1,1,var_u) = cfac * Te0 + 2.d0 * vw / a_n

    do id = 2, n_degrees
      if ( with_TiTe ) then
        node_list%node(i)%values(1,id,var_u) = cfac * node_list%node(i)%values(1,id,var_Te)
      else
        node_list%node(i)%values(1,id,var_u) = cfac * 0.5d0 * node_list%node(i)%values(1,id,var_T)
      endif
    enddo

    n_loc(1) = n_loc(1) + 1.d0

  enddo

  call MPI_Reduce(n_loc, n_glo, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if ( my_id .eq. 0 ) write(*,'(A,i0,A)')                                            &
    ' SHEATH: sheath_init_u set u to the floating potential on ', nint(n_glo(1)),    &
    ' boundary nodes'

end subroutine sheath_init_potential


!> Record psi's degrees of freedom at the start of the run, for the resistive wall relaxation.
!! Called once from jorek2_main after the restart is in place.
subroutine sheath_store_psi0(node_list)

  use mod_parameters
  use data_structure

  implicit none
  type (type_node_list), intent(in) :: node_list

  integer :: i, id

  if ( allocated(sheath_psi0) ) deallocate(sheath_psi0)
  allocate( sheath_psi0(n_degrees, node_list%n_nodes) )
  sheath_psi0 = 0.d0

  do i = 1, node_list%n_nodes
    do id = 1, n_degrees
      sheath_psi0(id,i) = node_list%node(i)%values(1,id,var_psi)
    enddo
  enddo

end subroutine sheath_store_psi0

end module mod_sheath_diag
