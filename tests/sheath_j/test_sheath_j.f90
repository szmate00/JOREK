!> Serial checks of the sheath current row (bcs%sheath_j) in the production boundary assembler:
!!  1. no zj row with sheath_j off, none below the grazing angle, present above it; never a u row
!!  2. every zj/u/rho/Ti/Te column of the zj rows vs central FD: at floating, on the ion side, beyond electron
!!     saturation, with vpar_smoothing, and with the tangent-continued characteristic (ion slope 1, electron slope
!!     e^Lambda) on the ion side beyond floating and beyond electron saturation
!!  3. the continued characteristic is continuous at floating (x = 0, where 1 - e^x - x has a kink by construction)
!!     and continuously differentiable at the electron cap (x = Lambda, with s_e = e^Lambda)
!!  4. the ion saturation current flows INTO the wall for both signs of F0
!!  (2b. j_sat carries no vpar_smoothing weight: the zj rows are identical with it on and off)
!!  5. sheath_heat_total_flow: without ExB the Ti/Te rows equal develop's; with ExB every psi/rho/Ti/Te/Vpar/u
!!     column of the Ti/Te rows matches FD at ExB outflow and at ExB inflow
!!  6. outer-plate bias: all edge nodes beyond R_xpoint = a uniform wall potential of the same value; none beyond
!!     it, or another boundary type = no effect; one node beyond it = in between (value-basis ramp)
program test_sheath_j
  use mod_parameters
  use phys_module
  use data_structure
  use mod_boundary_matrix_open
  use mod_floating_u, only: floating_u_norm
  use basis_at_gaussian, only: set_basis
  implicit none
  integer, parameter :: nd = 4*4*n_var
  integer, parameter :: fd_vars(5) = [var_zj, var_u, var_rho, var_Ti, var_Te]
  type(type_element) :: e
  type(type_node)    :: nodes(4), base(4)
  real*8  :: a(nd,nd), r(nd), ap(nd,nd), rp(nd), am(nd,nd), rm(nd), eps, worst, scale, a_n, C_T, C_V, ufl, rs
  real*8  :: rlo(nd), rhi(nd), ucap, us
  integer, parameter :: tf_vars(6) = [var_psi, var_rho, var_Ti, var_Te, var_vpar, var_u]
  integer :: ic2
  real*8  :: rxp = 1.d0                   ! R of the lower X-point passed to the assembler
  real*8  :: rb(nd), r0b(nd), vb
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

  ! ---------------------------------------------------------------- 2. FD of every column of the zj rows
  ! x = Lambda*(1 - u/ufl) on this state: u = ufl is floating, u > ufl the ion side, u < 0 beyond the cap
  do icase = 1, 7
    vpar_smoothing = .false. ; sheath_j_ion_slope = 0.d0 ; sheath_j_e_slope = 0.d0
    select case (icase)
    case (1); call set_u( 1.0d0*ufl)                                   ! floating
    case (2); call set_u( 1.5d0*ufl)                                   ! ion side, hard saturation
    case (3); call set_u(-0.5d0*ufl)                                   ! beyond electron saturation, hard cap: u column 0
    case (4); call set_u( 1.5d0*ufl); vpar_smoothing = .true.; vpar_smoothing_coef = [0.05d0, 0.02d0, 0.d0]
    case (5); call set_u( 1.5d0*ufl); sheath_j_ion_slope = 1.d0        ! continued ion branch
    case (6); call set_u(-0.5d0*ufl); sheath_j_e_slope = exp(sheath_Lambda)   ! continued electron branch
    case (7); call set_u( 0.5d0*ufl); sheath_j_ion_slope = 1.d0; sheath_j_e_slope = exp(sheath_Lambda)  ! between: slopes inactive
    end select
    nodes = base; call assemble(a, r)
    worst = 0.d0
    eps = 1.d-8
    if ( icase == 3 .or. icase == 6 ) eps = 1.d-6   ! beyond the cap the residual carries ~20 j_sat: 1e-8 steps drown in roundoff
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
        if ( icase == 3 .and. varof(col) == var_u .and. a(row,col) /= 0.d0 ) error stop 'FAIL 2: u column kept beyond a hard electron cap'
      enddo
      if ( icase == 6 .and. varof(col) == var_u .and. all(a(:,col) == 0.d0) .and. dofof(col) == 1 .and. nodeof(col) <= 2 ) &
        error stop 'FAIL 2: continued electron branch has no u column'
    enddo
    write(*,'(a,i2,a,es9.2)') ' PASS 2 case', icase, ': zj/u/rho/Ti/Te columns of the zj rows vs FD, worst ', worst
  enddo

  ! ---------------------------------------------------------------- 2b. j_sat does not carry the vpar_smoothing weight
  ! the nodal row applies the weight only across a tangency, so on a target the zj rows must be the same with it on or off
  call set_u( 1.5d0*ufl ); sheath_j_ion_slope = 0.d0 ; sheath_j_e_slope = 0.d0
  vpar_smoothing = .false. ; nodes = base; call assemble(a, r)
  vpar_smoothing = .true. ; vpar_smoothing_coef = [0.05d0, 0.02d0, 0.d0] ; nodes = base; call assemble(ap, rp)
  do row = 1, nd
    if ( varof(row) /= var_zj ) cycle
    if ( rp(row) /= r(row) .or. any(ap(row,:) /= a(row,:)) ) error stop 'FAIL 2b: zj rows depend on vpar_smoothing'
  enddo
  vpar_smoothing = .false.
  write(*,'(a)') ' PASS 2b: zj rows (j_sat) independent of the vpar_smoothing weight'

  ! ---------------------------------------------------------------- 3. continuity of the continued characteristic
  ! the residual is continuous in u across x = 0 and x = Lambda with the slopes on (a kink in df/dx at x = 0, which
  ! is the stated form 1 - e^x - x; C1 at x = Lambda with s_e = e^Lambda)
  sheath_j_ion_slope = 1.d0 ; sheath_j_e_slope = exp(sheath_Lambda)
  ucap = 0.d0                                   ! u = 0 gives x = Lambda on this state
  ! a jump shows as a difference across the point that does not shrink with the step: the ratio of the differences for
  ! steps 1e-6 and 1e-8 (relative to ufl) is ~100 for a continuous residual and ~1 for a discontinuous one
  if ( jump_ratio(ufl) < 30.d0 )  error stop 'FAIL 3: residual jumps at floating'
  if ( jump_ratio(ucap) < 30.d0 ) error stop 'FAIL 3: residual jumps at the electron cap'
  call set_u(ucap - 1.d-8*ufl); nodes = base; call assemble(a, rlo)
  call set_u(ucap + 1.d-8*ufl); nodes = base; call assemble(ap, rhi)
  do col = 1, nd
    if ( varof(col) /= var_u ) cycle
    if ( maxval(abs(a(:,col) - ap(:,col))) > 1.d-5 * max(1.d-30, maxval(abs(ap(:,col)))) ) error stop 'FAIL 3: u column not continuous at the cap (s_e = e^Lambda)'
  enddo
  sheath_j_ion_slope = 0.d0 ; sheath_j_e_slope = 0.d0
  write(*,'(a)') ' PASS 3: continued characteristic: residual continuous at floating and at the cap, slope continuous at the cap'

  ! ---------------------------------------------------------------- 5. sheath heat sinks on the total outgoing flow
  bcs(1)%sheath_j = .false.
  call set_u(ufl)
  ! (a) no ExB (u constant along the wall): the total and the parallel measure agree
  sheath_heat_total_flow = .false.; nodes = base; call assemble(a, r)
  sheath_heat_total_flow = .true. ; nodes = base; call assemble(ap, rp)
  worst = 0.d0
  do row = 1, nd
    if ( varof(row) /= var_Ti .and. varof(row) /= var_Te ) cycle
    worst = max(worst, abs(rp(row) - r(row)) / max(1.d-30, maxval(abs(r))))
  enddo
  if ( worst > 1.d-12 ) error stop 'FAIL 5a: total-flow sink differs from the parallel one without ExB'
  write(*,'(a,es9.2)') ' PASS 5a: without ExB the total-flow sheath heat sink equals the parallel one, rel diff ', worst
  ! (b) FD of every column the flag adds, at ExB outflow (us > 0 here) and at ExB inflow (us < 0): the difference
  !     flag on minus flag off of the Jacobian against the FD of the difference of the residual (develop's own
  !     columns are not complete: e.g. the c_angle floor's cs has no cross-temperature column, not tested here)
  do ic2 = 1, 2
    us = merge(0.01d0, -0.01d0, ic2 == 1)
    do i = 1, 4
      base(i)%values(1,1,var_u) = ufl + us*base(i)%x(1,1,1) ; base(i)%values(1,2,var_u) = us
    enddo
    sheath_heat_total_flow = .true. ; nodes = base; call assemble(a, r)
    sheath_heat_total_flow = .false.; nodes = base; call assemble(ap, rp)
    a = a - ap ; r = r - rp
    worst = 0.d0
    do col = 1, nd
      if ( .not. any(varof(col) == tf_vars) ) cycle
      eps = 1.d-8
      sheath_heat_total_flow = .true. ; nodes = base; call bump(col, +eps); call assemble(ap, rp)
      sheath_heat_total_flow = .false.; nodes = base; call bump(col, +eps); call assemble(am, rm)
      rlo = rp - rm
      sheath_heat_total_flow = .true. ; nodes = base; call bump(col, -eps); call assemble(ap, rp)
      sheath_heat_total_flow = .false.; nodes = base; call bump(col, -eps); call assemble(am, rm)
      rhi = rp - rm                                          ! rlo: difference at +eps, rhi: at -eps
      scale = max(1.d-30, maxval(abs(a(:,col))), maxval(abs(rlo-rhi))/(2*eps), 1.d-9*maxval(abs(rp))/eps)
      do row = 1, nd
        if ( varof(row) /= var_Ti .and. varof(row) /= var_Te ) cycle
        worst = max(worst, abs(a(row,col) + (rlo(row)-rhi(row))/(2*eps)) / scale)
        if ( abs(a(row,col) + (rlo(row)-rhi(row))/(2*eps)) / scale > 1.d-6 ) then
          write(*,'(a,4i5,2es14.5)') ' FAIL 5b: case,row,col,var, damat, -dfd', ic2, row, col, varof(col), a(row,col), -(rlo(row)-rhi(row))/(2*eps)
          error stop 1
        endif
      enddo
    enddo
    write(*,'(a,a,a,es9.2)') ' PASS 5b ', merge('outflow','inflow ', ic2 == 1), ': added Ti/Te-row columns vs FD of the added residual, worst ', worst
  enddo
  sheath_heat_total_flow = .false. ; bcs(1)%sheath_j = .true. ; call set_u(ufl)
  do i = 1, 4
    base(i)%values(1,2,var_u) = 0.d0
  enddo

  ! ---------------------------------------------------------------- 6. outer-plate bias
  call set_u(1.2d0*ufl) ; vb = 40.d0
  sheath_bias_V = 0.d0 ; sheath_V_wall = vb ; nodes = base; call assemble(a, rb)          ! uniform wall potential vb
  sheath_V_wall = 0.d0 ; nodes = base; call assemble(a, r0b)                              ! no bias
  sheath_bias_V = vb ; rxp = 0.5d0 ; nodes = base; call assemble(a, r)                    ! both nodes (R = 1, 2) beyond R_xpoint
  worst = 0.d0
  do row = 1, nd
    if ( varof(row) == var_zj ) worst = max(worst, abs(r(row) - rb(row)) / max(1.d-30, maxval(abs(rb))))
  enddo
  if ( worst > 1.d-12 ) error stop 'FAIL 6: bias on all nodes differs from the uniform wall potential'
  rxp = 3.d0 ; nodes = base; call assemble(a, r)                                          ! no node beyond R_xpoint
  if ( any(r /= r0b) ) error stop 'FAIL 6: bias acts on nodes inside R_xpoint'
  rxp = 0.5d0 ; sheath_bias_type = 2 ; nodes = base; call assemble(a, r)                  ! other boundary type
  if ( any(r /= r0b) ) error stop 'FAIL 6: bias acts on another boundary type'
  sheath_bias_type = 1 ; rxp = 1.5d0 ; nodes = base; call assemble(a, r)                   ! node 2 (R = 2) only
  do row = 1, nd
    if ( varof(row) /= var_zj ) cycle
    if ( (r(row) - r0b(row)) * (rb(row) - r(row)) < -1.d-12*maxval(abs(rb)) ) error stop 'FAIL 6: one-node bias not between none and all'
  enddo
  if ( all(r == r0b) .or. all(r == rb) ) error stop 'FAIL 6: one-node bias equals none or all'
  sheath_bias_V = 0.d0 ; rxp = 1.d0
  write(*,'(a,es9.2)') ' PASS 6: outer-plate bias: all nodes = uniform V_wall, none/other type = no effect, one node between; rel ', worst

  ! ---------------------------------------------------------------- 4. sign of the saturation current
  ! u far above floating: f -> 1, residual of the zj row with zj = 0 is +Zbig*dl*j_sat per unit test function.
  ! On this edge B_pol.n > 0, and the current into the wall is -zj*(B_pol.n)/F0, so j_sat/F0 must be negative.
  do icase = 1, 2
    F0 = merge(2.97d0, -2.97d0, icase == 1)
    call floating_u_norm(a_n, C_T, C_V); ufl = C_T * 0.004d0
    call set_u(50.d0*ufl); nodes = base; nodes(:)%values(1,1,var_zj) = 0.d0; nodes(:)%values(1,2,var_zj) = 0.d0
    call assemble(a, r)
    rs = r(var_zj) + r(4*n_var + var_zj)              ! value DOFs of the two wall nodes
    if ( rs / F0 .ge. 0.d0 ) then
      write(*,*) 'FAIL 4: j_sat/F0 not negative: F0, sum residual', F0, rs; error stop 1
    endif
  enddo
  F0 = 2.97d0
  write(*,'(a)') ' PASS 4: ion saturation current flows into the wall for both signs of F0'

contains
  real*8 function jump_ratio(u0)
    real*8, intent(in) :: u0
    real*8 :: d(2), h
    integer :: kk
    do kk = 1, 2
      h = merge(1.d-6, 1.d-8, kk == 1) * ufl
      call set_u(u0 - h); nodes = base; call assemble(a, rlo)
      call set_u(u0 + h); nodes = base; call assemble(a, rhi)
      d(kk) = maxval(abs(rhi - rlo))
    enddo
    jump_ratio = d(1) / max(d(2), tiny(1.d0))
  end function
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
      base(ii)%values(1,1,var_u) = u0 ; base(ii)%values(1,2,var_u) = 0.d0
    enddo
  end subroutine
  subroutine assemble(mat, rhs)
    real*8, intent(out) :: mat(nd,nd), rhs(nd)
    integer :: vertices(2), directions(2)
    vertices = [1,2]; directions = [1,2]
    mat = 0.d0; rhs = 0.d0
    call boundary_matrix_open(vertices, directions, e, nodes, .true., 1, 1.d0, 0.d0, 0.d0, 1.d0, &
                              [rxp,rxp], [0.d0,0.d0], mat, rhs, 1, 1)
  end subroutine
end program test_sheath_j
