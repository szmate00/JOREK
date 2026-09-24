!> Serial checks of the sheath current row (bcs%sheath_j, current-row form) in the production boundary assembler:
!!  1. no zj row with sheath_j off, none below the grazing angle, present (and no u row) above it
!!  2. every zj/u/rho/Ti/Te column of the zj rows vs central FD at floating potential, on the ion side (Phi above
!!     floating), beyond electron saturation (u column vanishes there) and with vpar_smoothing (factor < 1)
!!  3. the ion saturation current flows INTO the wall for both signs of F0
program test_sheath_j
  use mod_parameters
  use phys_module
  use data_structure
  use mod_boundary_matrix_open
  use mod_floating_u, only: floating_u_norm, sheath_j_ramp
  use basis_at_gaussian, only: set_basis
  implicit none
  integer, parameter :: nd = 4*4*n_var
  integer, parameter :: fd_vars(5) = [var_zj, var_u, var_rho, var_Ti, var_Te]
  type(type_element) :: e
  type(type_node)    :: nodes(4), base(4)
  real*8  :: a(nd,nd), r(nd), ap(nd,nd), rp(nd), am(nd,nd), rm(nd), eps, worst, scale, a_n, C_T, C_V, ufl, rs
  integer :: i, row, col, icase
  call set_basis()
  base(1)%x(1,1,:) = [1.d0, 0.d0]; base(2)%x(1,1,:) = [2.d0, 0.d0]
  base(3)%x(1,1,:) = [2.d0, 1.d0]; base(4)%x(1,1,:) = [1.d0, 1.d0]
  do i = 1, 4
    base(i)%x(1,2,:) = [1.d0, 0.d0]; base(i)%x(1,3,:) = [0.d0, 1.d0/3.d0]
    base(i)%values(1,1,var_rho) = 0.2d0 ; base(i)%values(1,2,var_rho) = 0.03d0
    base(i)%values(1,1,var_Ti)  = 0.003d0 ; base(i)%values(1,2,var_Ti) = 0.0004d0
    base(i)%values(1,1,var_Te)  = 0.004d0 ; base(i)%values(1,2,var_Te) = 0.0005d0
    base(i)%values(1,1,var_vpar)= 0.05d0
    base(i)%values(1,1,var_zj)  = 0.01d0 ; base(i)%values(1,2,var_zj) = 0.002d0
  enddo
  base%boundary = 1
  eps = 1.d-8
  call floating_u_norm(a_n, C_T, C_V)
  ufl = C_T * 0.004d0            ! floating u at the wall Te

  ! ---------------------------------------------------------------- 1. presence
  call set_psi(0.08d0); call set_u(ufl)
  bcs(1)%sheath_j = .false.; nodes = base; call assemble(a, r)
  if ( anyrow(var_zj) ) error stop 'FAIL 1: zj row assembled with sheath_j off'
  bcs(1)%sheath_j = .true.
  call set_psi(0.02d0); nodes = base; call assemble(a, r)        ! |b.n| ~ 0.007 < sin(1 deg)
  if ( anyrow(var_zj) ) error stop 'FAIL 1: zj row assembled below the grazing angle'
  call set_psi(0.08d0); nodes = base; call assemble(a, r)        ! |b.n| ~ 0.027
  if ( .not. anyrow(var_zj) ) error stop 'FAIL 1: no zj row above the grazing angle'
  if ( anyrow(var_u) ) error stop 'FAIL 1: a u row was assembled (u belongs to the vorticity equation)'
  write(*,'(a)') ' PASS 1: zj row: none with sheath_j off, none below the angle, present above it; no u row'
  ! the total-flow wall flux must be on for a sheath type too: with u varying along the wall the rho rows get a u column
  call set_u(ufl); base(:)%values(1,2,var_u) = 0.05d0; nodes = base; call assemble(a, r)
  worst = 0.d0
  do row = 1, nd
    if ( varof(row) /= var_rho ) cycle
    do col = 1, nd
      if ( varof(col) == var_u ) worst = max(worst, abs(a(row,col)))
    enddo
  enddo
  if ( worst .eq. 0.d0 ) error stop 'FAIL 1: total-flow wall flux not active on a sheath_j type (no u column in the rho rows)'
  write(*,'(a)') ' PASS 1b: total-flow wall flux active on the sheath type (rho rows carry the u column)'

  ! ---------------------------------------------------------------- 2. FD
  do icase = 1, 4
    vpar_smoothing = .false.
    select case (icase)
    case (1); call set_u( 1.0d0*ufl)                                   ! floating: x = 0, zj should relax to 0
    case (2); call set_u( 1.5d0*ufl)                                   ! ion side, x < 0
    case (3); call set_u(-0.5d0*ufl)                                   ! beyond electron saturation, x > Lambda
    case (4); call set_u( 1.5d0*ufl); vpar_smoothing = .true.; vpar_smoothing_coef = [0.05d0, 0.02d0, 0.d0]
    end select
    nodes = base; call assemble(a, r)
    worst = 0.d0
    eps = 1.d-8
    if ( icase == 3 ) eps = 1.d-6      ! capped: the residual carries (e^Lambda - 1)*j_sat ~ 20x more, 1e-8 steps drown in roundoff
    do col = 1, nd
      if ( .not. any(varof(col) == fd_vars) ) cycle
      nodes = base; call bump(col, +eps); call assemble(ap, rp)
      nodes = base; call bump(col, -eps); call assemble(am, rm)
      scale = max(1.d-30, maxval(abs(a(:,col))), maxval(abs(rp-rm))/(2*eps))
      do row = 1, nd
        if ( varof(row) /= var_zj ) cycle
        worst = max(worst, abs(a(row,col) + (rp(row)-rm(row))/(2*eps)) / scale)
        if ( abs(a(row,col) + (rp(row)-rm(row))/(2*eps)) / scale > 1.d-6 ) then
          write(*,'(a,4i5,2es14.5)') ' FAIL 2: case,row,col,var, amat, -fd', icase, row, col, varof(col), a(row,col), -(rp(row)-rm(row))/(2*eps)
          error stop 1
        endif
        if ( icase == 3 .and. varof(col) == var_u .and. a(row,col) /= 0.d0 ) error stop 'FAIL 2: u column kept beyond electron saturation'
      enddo
    enddo
    write(*,'(a,i2,a,es9.2)') ' PASS 2 case', icase, ': zj/u/rho/Ti/Te columns of the zj rows vs FD, worst ', worst
  enddo

  ! ---------------------------------------------------------------- 2c. ion-saturation slope and switch-on ramp
  ! FD of every column of the zj rows with s > 0 on the ion side, and with alpha < 1 in all three states; the ramp
  ! function: alpha0 at t = 0, 1 at the end of the timestep ramp, linear between, 1 with the ramp disabled.
  do icase = 1, 5
    sheath_j_ion_slope = 0.d0 ; sheath_j_ramp_time = -1.d0 ; t_now = 0.d0
    select case (icase)
    case (1); sheath_j_ion_slope = 0.03d0 ; call set_u( 1.5d0*ufl)                 ! ion side with slope
    case (2); sheath_j_ion_slope = 0.03d0 ; call set_u( 0.6d0*ufl)                 ! electron side: slope inactive
    case (3); sheath_j_ramp_time = 1.d0 ; t_now = 0.d0 ; call set_u( 1.05d0*ufl)   ! alpha0 = 0.05, ion side
    case (4); sheath_j_ramp_time = 1.d0 ; t_now = 0.d0 ; call set_u( 0.98d0*ufl)   ! alpha0, electron side below the cap
    case (5); sheath_j_ramp_time = 1.d0 ; t_now = 0.5d0 ; sheath_j_ion_slope = 0.03d0 ; call set_u( 1.3d0*ufl)  ! mid-ramp with slope
    end select
    nodes = base; call assemble(a, r)
    worst = 0.d0
    do col = 1, nd
      if ( .not. any(varof(col) == fd_vars) ) cycle
      eps = 1.d-8
      ! alpha = 0.05: the exponent varies 20x faster in u and Te, central FD needs a smaller step on those columns
      if ( icase >= 3 .and. (varof(col) == var_u .or. varof(col) == var_Te) ) eps = 3.d-10
      nodes = base; call bump(col, +eps); call assemble(ap, rp)
      nodes = base; call bump(col, -eps); call assemble(am, rm)
      scale = max(1.d-30, maxval(abs(a(:,col))), maxval(abs(rp-rm))/(2*eps))
      do row = 1, nd
        if ( varof(row) /= var_zj ) cycle
        worst = max(worst, abs(a(row,col) + (rp(row)-rm(row))/(2*eps)) / scale)
        if ( abs(a(row,col) + (rp(row)-rm(row))/(2*eps)) / scale > 1.d-6 ) then
          write(*,'(a,4i5,2es14.5)') ' FAIL 2c: case,row,col,var, amat, -fd', icase, row, col, varof(col), a(row,col), -(rp(row)-rm(row))/(2*eps)
          error stop 1
        endif
      enddo
    enddo
    write(*,'(a,i2,a,es9.2)') ' PASS 2c case', icase, ': slope/ramp: zj-row columns vs FD, worst ', worst
  enddo
  sheath_j_ion_slope = 0.d0 ; sheath_j_ramp_time = -1.d0 ; t_now = 0.d0
  ! ramp function
  tstep_n(1:3) = [1.d-2, 1.d-1, 1.d0] ; nstep_n(1:3) = [10, 10, 1000]      ! ramp phases end at t = 1.1
  sheath_j_ramp_time = 0.d0
  if ( abs(sheath_j_ramp(0.d0) - 0.05d0) > 1.d-14 )  error stop 'FAIL 2c: alpha(0) /= alpha0'
  if ( abs(sheath_j_ramp(0.55d0) - 0.525d0) > 1.d-14 ) error stop 'FAIL 2c: alpha not linear over the timestep ramp'
  if ( sheath_j_ramp(1.1d0) /= 1.d0 .or. sheath_j_ramp(5.d0) /= 1.d0 ) error stop 'FAIL 2c: alpha /= 1 after the ramp'
  sheath_j_ramp_time = -1.d0
  if ( sheath_j_ramp(0.d0) /= 1.d0 ) error stop 'FAIL 2c: alpha /= 1 with the ramp disabled'
  tstep_n = 0.d0 ; nstep_n = 0
  write(*,'(a)') ' PASS 2c: ramp alpha: alpha0 at t = 0, linear to 1 at the end of the timestep ramp, 1 when disabled'

  ! ---------------------------------------------------------------- 3. sign of the saturation current
  ! u far above floating: f -> 1, residual of the zj row with zj = 0 is +Zbig*dl*j_sat per unit test function.
  ! On this edge B_pol.n > 0, and the current into the wall is -zj*(B_pol.n)/F0, so j_sat/F0 must be negative.
  do icase = 1, 2
    F0 = merge(2.97d0, -2.97d0, icase == 1)
    call floating_u_norm(a_n, C_T, C_V); ufl = C_T * 0.004d0
    call set_u(50.d0*ufl); nodes = base; nodes(:)%values(1,1,var_zj) = 0.d0; nodes(:)%values(1,2,var_zj) = 0.d0
    call assemble(a, r)
    rs = r(var_zj) + r(4*n_var + var_zj)              ! value DOFs of the two wall nodes
    if ( rs / F0 .ge. 0.d0 ) then
      write(*,*) 'FAIL 3: j_sat/F0 not negative: F0, sum residual', F0, rs; error stop 1
    endif
  enddo
  F0 = 2.97d0
  write(*,'(a)') ' PASS 3: ion saturation current flows into the wall for both signs of F0'

contains
  integer function varof(idx)
    integer, intent(in) :: idx
    varof = mod(idx-1, n_var) + 1
  end function
  integer function dofof(idx)
    integer, intent(in) :: idx
    dofof = mod((idx-1)/n_var, 4) + 1
  end function
  integer function nodeof(idx)
    integer, intent(in) :: idx
    nodeof = (idx-1)/(4*n_var) + 1
  end function
  logical function anyrow(var)
    integer, intent(in) :: var
    integer :: ir
    anyrow = .false.
    do ir = 1, nd
      if ( varof(ir) == var .and. ( any(a(ir,:) /= 0.d0) .or. r(ir) /= 0.d0 ) ) anyrow = .true.
    enddo
  end function
  subroutine bump(idx, d)
    integer, intent(in) :: idx
    real*8,  intent(in) :: d
    nodes(nodeof(idx))%values(1,dofof(idx),varof(idx)) = nodes(nodeof(idx))%values(1,dofof(idx),varof(idx)) + d
  end subroutine
  subroutine set_psi(pslope)
    real*8, intent(in) :: pslope
    integer :: ii
    do ii = 1, 4
      base(ii)%values(1,1,var_psi) = pslope*base(ii)%x(1,1,1) ; base(ii)%values(1,2,var_psi) = pslope ; base(ii)%values(1,3,var_psi) = 0.25d0*pslope
    enddo
  end subroutine
  subroutine set_u(u0)
    real*8, intent(in) :: u0
    integer :: ii
    do ii = 1, 4
      base(ii)%values(1,1,var_u) = u0 ; base(ii)%values(1,2,var_u) = 0.1d0*u0
    enddo
  end subroutine
  subroutine assemble(mat, rhs)
    real*8, intent(out) :: mat(nd,nd), rhs(nd)
    integer :: vertices(2), directions(2)
    vertices = [1,2]; directions = [1,2]
    mat = 0.d0; rhs = 0.d0
    call boundary_matrix_open(vertices, directions, e, nodes, .true., 1, 1.d0, 0.d0, 0.d0, 1.d0, &
                              [1.d0,1.d0], [0.d0,0.d0], mat, rhs, 1, 1)
  end subroutine
end program test_sheath_j
