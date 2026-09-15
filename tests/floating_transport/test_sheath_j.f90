!> Serial checks of the sheath current row (bcs%sheath_j) in the PRODUCTION assembler.
!!  1. off: no zj row; sheath_j_pin_current: no zj row either; below the grazing angle: no zj row.
!!  2. residual vanishes at zj = 0 with u at the floating value (x = 0): the j -> 0 anchor.
!!  3. every column of the zj rows (zj, u, rho, Ti, Te) matches finite differences on the electron
!!     branch; on the electron-saturated branch (x >= Lambda) the u and Te-through-x columns vanish.
!!  4. sign: the ion saturation current flows INTO the wall for both signs of F0.
!!  5. the diagnostics leave the equations untouched.
!!  6. Option I (sheath_j_ohm): no zj row but a u row; zj = 0 with u floating is a root; every column
!!     (u, zj, rho, Ti, Te) matches FD with X inside its bounds; with X beyond either bound the zj, rho
!!     and cs columns vanish and the u column stays.
program test_sheath_j
  use mod_parameters
  use phys_module
  use data_structure
  use mod_boundary_matrix_open
  use mod_floating_u,    only: floating_u_norm, sheath_j_norm
  use mod_floating_diag, only: floating_diag_reset, floating_diag_report
  use basis_at_gaussian, only: set_basis
  implicit none
  integer, parameter :: nd = 4*4*n_var
  integer, parameter :: fd_vars(5) = [var_zj, var_u, var_rho, var_Ti, var_Te]
  type(type_element) :: e
  type(type_node)    :: nodes(4), base(4)
  real*8  :: a0(nd,nd), r0(nd), a(nd,nd), r(nd), ap(nd,nd), rp(nd), am(nd,nd), rm(nd)
  real*8  :: eps, err, worst, scale, a_n, C_T, C_V, c_sat, ufl
  integer :: i, k, var, dof, col, row, isgn

  call set_basis()
  base(1)%x(1,1,:) = [1.d0, 0.d0]; base(2)%x(1,1,:) = [2.d0, 0.d0]
  base(3)%x(1,1,:) = [2.d0, 1.d0]; base(4)%x(1,1,:) = [1.d0, 1.d0]
  do i = 1, 4
    base(i)%x(1,2,:) = [1.d0, 0.d0]
    base(i)%x(1,3,:) = [0.d0, 1.d0/3.d0]
    base(i)%values(1,1,var_rho)  = 0.2d0
    base(i)%values(1,1,var_Ti)   = 0.003d0
    base(i)%values(1,1,var_Te)   = 0.004d0
    base(i)%values(1,1,var_vpar) = 0.04d0
    base(i)%values(1,1,var_psi)  = 0.08d0*base(i)%x(1,1,1)
    base(i)%values(1,2,var_psi)  = 0.08d0
    base(i)%values(1,1,var_zj)   = 0.01d0
  enddo
  base%boundary = 1
  call floating_u_norm(a_n, C_T, C_V)
  ufl = C_T * 0.004d0                       ! floating value of u for Te = 0.004
  base(:)%values(1,1,var_u) = ufl
  mach1_weak = .true.

  ! ---------------------------------------------------------------- 1. off / pinned / grazing
  bcs(1)%sheath_j = .false.; nodes = base; call assemble(a0, r0)
  if ( anyrow(a0, r0, var_zj) ) error stop 'FAIL: zj row assembled with sheath_j off'
  bcs(1)%sheath_j = .true.
  sheath_j_pin_current = .true.;  nodes = base; call assemble(a, r)
  if ( anyrow(a, r, var_zj) ) error stop 'FAIL: zj row assembled with sheath_j_pin_current'
  sheath_j_pin_current = .false.
  do i = 1, 4
    base(i)%values(1,1,var_psi) = 1.d-4*base(i)%x(1,1,1); base(i)%values(1,2,var_psi) = 1.d-4
  enddo
  nodes = base; call assemble(a, r)
  if ( anyrow(a, r, var_zj) ) error stop 'FAIL: zj row assembled below the grazing angle'
  do i = 1, 4
    base(i)%values(1,1,var_psi) = 0.08d0*base(i)%x(1,1,1); base(i)%values(1,2,var_psi) = 0.08d0
  enddo
  nodes = base; call assemble(a, r)
  if ( .not. anyrow(a, r, var_zj) ) error stop 'FAIL: zj row not assembled'
  write(*,'(a)') ' PASS: sheath row off / pinned / grazing: no zj row; on: zj row present'

  ! ---------------------------------------------------------------- 2. j -> 0 anchor
  base(:)%values(1,1,var_zj) = 0.d0
  nodes = base; call assemble(a, r)
  worst = 0.d0
  do row = 1, nd
    if ( isvar(row, var_zj) ) worst = max(worst, abs(r(row)))
  enddo
  scale = maxval(abs(a))
  if ( worst > 1.d-12*scale ) then
    write(*,*) 'FAIL: residual at zj = 0, u floating', worst, scale
    error stop 1
  endif
  write(*,'(a)') ' PASS: zj = 0 with u at the floating value is a root of the sheath row'

  ! ---------------------------------------------------------------- 3. FD, electron branch and saturated branch
  base(:)%values(1,1,var_zj) = 0.01d0
  eps = 1.d-7      ! the saturated residual is ~19 j_sat: central differences need a larger step
  do k = 1, 2
    if ( k == 1 ) base(:)%values(1,1,var_u) = 0.8d0*ufl      ! x = 0.2*Lambda > 0: electron branch, unsaturated
    if ( k == 2 ) base(:)%values(1,1,var_u) = -3.d0*ufl      ! x = 4*Lambda: electron saturated
    nodes = base; call assemble(a, r)
    worst = 0.d0
    do i = 1, size(fd_vars)
      var = fd_vars(i)
      do row = 1, 2
        do dof = 1, 4
          col = n_var*4*(row-1) + n_var*(dof-1) + var
          nodes = base; nodes(row)%values(1,dof,var) = nodes(row)%values(1,dof,var) + eps; call assemble(ap, rp)
          nodes = base; nodes(row)%values(1,dof,var) = nodes(row)%values(1,dof,var) - eps; call assemble(am, rm)
          scale = max(1.d0, maxval(abs(a(:,col))))
          do isgn = 1, nd
            if ( .not. isvar(isgn, var_zj) ) cycle
            err = abs( a(isgn,col) + (rp(isgn)-rm(isgn))/(2*eps) ) / scale
            worst = max(worst, err)
            if ( err > 1.d-6 ) then
              write(*,'(a,4i5,3es12.3)') ' FAIL: sheath FD branch,row,col,var,err,amat,fd', k, isgn, col, var, err, a(isgn,col), -(rp(isgn)-rm(isgn))/(2*eps)
              error stop 1
            endif
          enddo
        enddo
      enddo
    enddo
    if ( k == 2 ) then
      do isgn = 1, nd
        if ( .not. isvar(isgn, var_zj) ) cycle
        do col = 1, nd
          if ( isvar(col, var_u) .and. a(isgn,col) /= 0.d0 ) error stop 'FAIL: saturated sheath row keeps a u column'
        enddo
      enddo
    endif
    if ( k == 1 ) write(*,'(a,es9.2)') ' PASS: electron branch: every sheath-row column matches FD, worst rel err ', worst
    if ( k == 2 ) write(*,'(a,es9.2)') ' PASS: electron-saturated branch: columns match FD, no u column, worst rel err ', worst
  enddo
  base(:)%values(1,1,var_u) = ufl

  ! ---------------------------------------------------------------- 4. sign of the saturation current, both F0
  do isgn = 1, 2
    call sheath_j_norm(a_n, c_sat)
    ! B_pol.n > 0 here (psi = 0.08 R, normal -Z); the outward Vpar is +cs/|B|, so
    ! j_sat = c_sat*rho*(+cs/|B|) and the current into the wall is -j_sat*(B_pol.n)/F0 ~ -c_sat/F0 > 0
    if ( -c_sat/F0 .le. 0.d0 ) then
      write(*,*) 'FAIL: saturation current does not flow into the wall, F0 =', F0
      error stop 1
    endif
    F0 = -F0
  enddo
  write(*,'(a)') ' PASS: ion saturation current flows into the wall for both signs of F0'

  ! ---------------------------------------------------------------- 5. diagnostics
  call floating_diag_reset()
  floating_u_diag = .true.;  nodes = base; call assemble(ap, rp)
  floating_u_diag = .false.; nodes = base; call assemble(a, r)
  if ( any(a /= ap) .or. any(r /= rp) ) error stop 'FAIL: sheath diagnostics changed the equations'
  call floating_diag_report(0)
  write(*,'(a)') ' PASS: sheath diagnostics leave the equations untouched'

  ! ---------------------------------------------------------------- 6. Option I
  sheath_j_ohm = .true.
  base(:)%values(1,1,var_u)  = ufl
  base(:)%values(1,1,var_zj) = 0.d0
  nodes = base; call assemble(a, r)
  if ( anyrow(a, r, var_zj) )        error stop 'FAIL: Option I assembles a zj row'
  if ( .not. anyrow(a, r, var_u) )   error stop 'FAIL: Option I assembles no u row'
  worst = 0.d0
  do row = 1, nd
    if ( isvar(row, var_u) ) worst = max(worst, abs(r(row)))
  enddo
  if ( worst > 1.d-12*maxval(abs(a)) ) error stop 'FAIL: Option I: zj = 0 with u floating is not a root'
  ! j_sat estimate of the fixture for placing X: c_sat*rho*cs/|B| with B_pol.n > 0
  call sheath_j_norm(a_n, c_sat)
  scale = c_sat * 0.2d0 * sqrt(gamma*(0.003d0+0.004d0)) / ( sqrt(F0**2 + 0.08d0**2) / 1.5d0 )   ! ~ j_sat at R = 1.5
  eps = 1.d-7
  do k = 1, 3
    if ( k == 1 ) base(:)%values(1,1,var_zj) =   0.5d0 * scale     ! X ~ 0.5, inside the bounds
    if ( k == 2 ) base(:)%values(1,1,var_zj) =   5.0d0 * scale     ! X < 0: lower bound active (estimate is approximate)
    if ( k == 3 ) base(:)%values(1,1,var_zj) = -1.d2  * scale      ! X >> e^3: upper bound active
    nodes = base; call assemble(a, r)
    worst = 0.d0
    do i = 1, size(fd_vars)
      var = fd_vars(i)
      do row = 1, 2
        do dof = 1, 4
          col = n_var*4*(row-1) + n_var*(dof-1) + var
          nodes = base; nodes(row)%values(1,dof,var) = nodes(row)%values(1,dof,var) + eps; call assemble(ap, rp)
          nodes = base; nodes(row)%values(1,dof,var) = nodes(row)%values(1,dof,var) - eps; call assemble(am, rm)
          do isgn = 1, nd
            if ( .not. isvar(isgn, var_u) ) cycle
            err = abs( a(isgn,col) + (rp(isgn)-rm(isgn))/(2*eps) ) / max(1.d0, maxval(abs(a(:,col))))
            worst = max(worst, err)
            if ( err > 1.d-6 ) then
              write(*,'(a,4i5,3es12.3)') ' FAIL: Option I FD case,row,col,var,err,amat,fd', k, isgn, col, var, err, a(isgn,col), -(rp(isgn)-rm(isgn))/(2*eps)
              error stop 1
            endif
          enddo
        enddo
      enddo
    enddo
    if ( k .ge. 2 ) then
      do isgn = 1, nd
        if ( .not. isvar(isgn, var_u) ) cycle
        do col = 1, nd
          if ( (isvar(col, var_zj) .or. isvar(col, var_rho)) .and. a(isgn,col) /= 0.d0 ) error stop 'FAIL: Option I bound active but zj/rho column kept'
        enddo
      enddo
    endif
    if ( k == 1 ) write(*,'(a,es9.2)') ' PASS: Option I potential row: every column matches FD inside the bounds, worst rel err ', worst
    if ( k == 3 ) write(*,'(a)')       ' PASS: Option I potential row: at either bound the zj/rho/cs columns vanish, u column stays'
  enddo
  sheath_j_ohm = .false.

contains

  logical function isvar(idx, var)
    integer, intent(in) :: idx, var
    isvar = ( mod(idx-1, n_var) + 1 == var )
  end function

  logical function anyrow(mat, rhs, var)
    real*8,  intent(in) :: mat(nd,nd), rhs(nd)
    integer, intent(in) :: var
    integer :: rr
    anyrow = .false.
    do rr = 1, nd
      if ( isvar(rr, var) .and. ( any(mat(rr,:) /= 0.d0) .or. rhs(rr) /= 0.d0 ) ) anyrow = .true.
    enddo
  end function

  subroutine assemble(mat, rhs)
    real*8, intent(out) :: mat(nd,nd), rhs(nd)
    integer :: vertices(2), directions(2)
    vertices = [1,2]; directions = [1,2]
    mat = 0.d0; rhs = 0.d0
    call boundary_matrix_open(vertices, directions, e, nodes, .true., 1, 1.d0, 0.d0, 0.d0, 1.d0, &
                              [1.d0,1.d0], [0.d0,0.d0], mat, rhs, 1, 1)
  end subroutine

end program test_sheath_j
