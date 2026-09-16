!> Serial checks of the sheath current BC (bcs%sheath_j) in the PRODUCTION assembler:
!!  1. off: no u or zj row; below the grazing angle: the surface term but no potential row; on: both;
!!     sheath_j_float_u: zj row (surface term) present, no u row.
!!  2. surface term of the current definition: for psi = a*R + b*Z on the R-edge with outward normal -Z the
!!     value-DOF rows of zj must sum to  - b * int dl/R = - b * ln 2  (partition of unity of the value basis).
!!  3. every psi column of the zj rows (value, s-, t- and st-derivative DOFs of both edge nodes) matches a
!!     central finite difference.
!!  4. potential row: zj = 0 with u at the floating value is a root; every column (u, zj, rho, Ti, Te) matches
!!     FD with X inside its bounds; at either bound the zj/rho/cs columns vanish and the u column stays.
!!  5. the ion saturation current flows into the wall for both signs of F0.
!!  6. the diagnostics leave the equations untouched.
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
  integer, parameter :: fd_vars(5) = [var_u, var_zj, var_rho, var_Ti, var_Te]
  type(type_element) :: e
  type(type_node)    :: nodes(4), base(4)
  real*8  :: a(nd,nd), r(nd), ap(nd,nd), rp(nd), am(nd,nd), rm(nd)
  real*8  :: eps, err, worst, scale, a_n, C_T, C_V, c_sat, ufl, bslope, tot
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
    base(i)%values(1,3,var_psi)  = 0.02d0
    base(i)%values(1,1,var_zj)   = 0.d0
  enddo
  base%boundary = 1
  call floating_u_norm(a_n, C_T, C_V)
  ufl = C_T * 0.004d0
  base(:)%values(1,1,var_u) = ufl
  mach1_weak = .true.

  ! ---------------------------------------------------------------- 1. rows present or not
  bcs(1)%sheath_j = .false.; nodes = base; call assemble(a, r)
  if ( anyrow(a, r, var_zj) .or. anyrow(a, r, var_u) ) error stop 'FAIL: sheath rows assembled with sheath_j off'
  bcs(1)%sheath_j = .true.
  do i = 1, 4
    base(i)%values(1,1,var_psi) = 1.d-4*base(i)%x(1,1,1); base(i)%values(1,2,var_psi) = 1.d-4
  enddo
  nodes = base; call assemble(a, r)
  if ( anyrow(a, r, var_u) )        error stop 'FAIL: potential row assembled below the grazing angle'
  if ( .not. anyrow(a, r, var_zj) ) error stop 'FAIL: surface term missing below the grazing angle (it follows the node, not the angle)'
  do i = 1, 4
    base(i)%values(1,1,var_psi) = 0.08d0*base(i)%x(1,1,1); base(i)%values(1,2,var_psi) = 0.08d0
  enddo
  nodes = base; call assemble(a, r)
  if ( .not. (anyrow(a, r, var_zj) .and. anyrow(a, r, var_u)) ) error stop 'FAIL: sheath rows not assembled'
  sheath_j_float_u = .true.;  nodes = base; call assemble(ap, rp)
  if ( anyrow(ap, rp, var_u) )        error stop 'FAIL: sheath_j_float_u assembles a u row'
  if ( .not. anyrow(ap, rp, var_zj) ) error stop 'FAIL: sheath_j_float_u drops the surface term'
  sheath_j_float_u = .false.
  write(*,'(a)') ' PASS: sheath rows: off none, grazing surface term only, on both, float_u surface term only'

  ! ---------------------------------------------------------------- 2. surface term value
  ! psi = a*R + b*Z: dpsi/dn = -b on the edge (outward normal -Z); nodal t-derivative DOF = b*dZ/dt = b/3
  bslope = 0.06d0
  do i = 1, 4
    base(i)%values(1,1,var_psi) = 0.08d0*base(i)%x(1,1,1) + bslope*base(i)%x(1,1,2)
    base(i)%values(1,2,var_psi) = 0.08d0
    base(i)%values(1,3,var_psi) = bslope/3.d0
    base(i)%values(1,4,var_psi) = 0.d0
  enddo
  nodes = base; call assemble(a, r)
  tot = r(n_var*4*0 + var_zj) + r(n_var*4*1 + var_zj)      ! value-DOF rows of nodes 1 and 2
  if ( abs(tot - (-bslope*log(2.d0))) > 1.d-6*bslope ) then    ! 4-point Gauss quadrature of 1/R
    write(*,*) 'FAIL: surface term integral', tot, -bslope*log(2.d0)
    error stop 1
  endif
  write(*,'(a)') ' PASS: surface term integrates to -b*ln2 for psi = a*R + b*Z'

  ! ---------------------------------------------------------------- 3. FD of the zj rows in every psi DOF
  eps = 1.d-7
  worst = 0.d0
  do i = 1, 2
    do dof = 1, 4
      col = n_var*4*(i-1) + n_var*(dof-1) + var_psi
      nodes = base; nodes(i)%values(1,dof,var_psi) = nodes(i)%values(1,dof,var_psi) + eps; call assemble(ap, rp)
      nodes = base; nodes(i)%values(1,dof,var_psi) = nodes(i)%values(1,dof,var_psi) - eps; call assemble(am, rm)
      scale = max(1.d0, maxval(abs(a(:,col))))
      do row = 1, nd
        if ( .not. isvar(row, var_zj) ) cycle
        err = abs( a(row,col) + (rp(row)-rm(row))/(2*eps) ) / scale
        worst = max(worst, err)
        if ( err > 1.d-6 ) then
          write(*,'(a,3i5,3es12.3)') ' FAIL: surface term FD node,dof,row,err,amat,fd', i, dof, row, err, a(row,col), -(rp(row)-rm(row))/(2*eps)
          error stop 1
        endif
      enddo
    enddo
  enddo
  write(*,'(a,es9.2)') ' PASS: surface term: every psi column (trace and normal-derivative DOFs) matches FD, worst ', worst
  do i = 1, 4
    base(i)%values(1,1,var_psi) = 0.08d0*base(i)%x(1,1,1); base(i)%values(1,3,var_psi) = 0.02d0
  enddo

  ! ---------------------------------------------------------------- 4. potential row
  nodes = base; call assemble(a, r)
  worst = 0.d0
  do row = 1, nd
    if ( isvar(row, var_u) ) worst = max(worst, abs(r(row)))
  enddo
  if ( worst > 1.d-12*maxval(abs(a)) ) error stop 'FAIL: zj = 0 with u floating is not a root of the potential row'
  call sheath_j_norm(a_n, c_sat)
  scale = c_sat * 0.2d0 * sqrt(gamma*(0.003d0+0.004d0)) / ( sqrt(F0**2 + 0.08d0**2) / 1.5d0 )   ! ~ j_sat at R = 1.5
  do k = 1, 3
    if ( k == 1 ) base(:)%values(1,1,var_zj) =   0.5d0 * scale     ! X ~ 0.5, inside the bounds
    if ( k == 2 ) base(:)%values(1,1,var_zj) =   5.0d0 * scale     ! X < 0: lower bound active
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
              write(*,'(a,4i5,3es12.3)') ' FAIL: potential row FD case,row,col,var,err,amat,fd', k, isgn, col, var, err, a(isgn,col), -(rp(isgn)-rm(isgn))/(2*eps)
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
          if ( (isvar(col, var_zj) .or. isvar(col, var_rho)) .and. a(isgn,col) /= 0.d0 ) error stop 'FAIL: bound active but zj/rho column kept'
        enddo
      enddo
    endif
    if ( k == 1 ) write(*,'(a,es9.2)') ' PASS: potential row: every column matches FD inside the bounds, worst rel err ', worst
    if ( k == 3 ) write(*,'(a)')       ' PASS: potential row: at either bound the zj/rho/cs columns vanish, u column stays'
  enddo
  base(:)%values(1,1,var_zj) = 0.d0

  ! --- ramp: alpha = 0 must reduce the potential row to the floating row (u column only, residual u - ufl),
  ! --- alpha from the timestep ramp must reach 1 when the last tstep_n phase begins
  base(:)%values(1,1,var_zj) = 0.5d0 * scale
  sheath_j_ramp_time = 0.d0 ; t_now = 0.d0
  tstep_n(1:3) = [1.d-3, 1.d-2, 1.d0] ; nstep_n(1:3) = [100, 100, 1000]
  nodes = base; call assemble(a, r)
  do isgn = 1, nd
    if ( .not. isvar(isgn, var_u) ) cycle
    do col = 1, nd
      if ( (isvar(col, var_zj) .or. isvar(col, var_rho)) .and. a(isgn,col) /= 0.d0 ) error stop 'FAIL: alpha = 0 keeps a current column'
    enddo
  enddo
  t_now = 0.55d0 ; nodes = base; call assemble(ap, rp)      ! half way through the ramp (t_end = 1.1)
  t_now = 1.1d0  ; nodes = base; call assemble(am, rm)      ! ramp end: full characteristic
  err = 0.d0 ; worst = 0.d0
  do isgn = 1, nd
    if ( .not. isvar(isgn, var_u) ) cycle
    do col = 1, nd
      if ( .not. isvar(col, var_zj) ) cycle
      err = max(err, abs(ap(isgn,col) - 0.5d0*am(isgn,col))) ; worst = max(worst, abs(am(isgn,col)))
    enddo
  enddo
  if ( worst <= 0.d0 .or. err > 1.d-12*worst ) error stop 'FAIL: ramp factor is not linear in time / does not reach 1'
  sheath_j_ramp_time = -1.d0 ; t_now = 0.d0 ; base(:)%values(1,1,var_zj) = 0.d0
  write(*,'(a)') ' PASS: sheath ramp: alpha = 0 is the floating row, alpha follows the timestep ramp to 1'

  ! --- current-slot form: zj row is the characteristic, no u row, no surface term; FD every column
  sheath_j_current_row = .true.
  base(:)%values(1,1,var_zj) = 0.01d0 ; base(:)%values(1,1,var_u) = 0.8d0*ufl
  nodes = base; call assemble(a, r)
  if ( anyrow(a, r, var_u) )        error stop 'FAIL: current-slot form assembles a u row'
  if ( .not. anyrow(a, r, var_zj) ) error stop 'FAIL: current-slot form assembles no zj row'
  worst = 0.d0
  do i = 1, size(fd_vars)
    var = fd_vars(i)
    do row = 1, 2
      do dof = 1, 4
        col = n_var*4*(row-1) + n_var*(dof-1) + var
        nodes = base; nodes(row)%values(1,dof,var) = nodes(row)%values(1,dof,var) + eps; call assemble(ap, rp)
        nodes = base; nodes(row)%values(1,dof,var) = nodes(row)%values(1,dof,var) - eps; call assemble(am, rm)
        do isgn = 1, nd
          if ( .not. isvar(isgn, var_zj) ) cycle
          err = abs( a(isgn,col) + (rp(isgn)-rm(isgn))/(2*eps) ) / max(1.d0, maxval(abs(a(:,col))))
          worst = max(worst, err)
          if ( err > 1.d-6 ) then
            write(*,'(a,3i5,3es12.3)') ' FAIL: current-slot FD row,col,var,err,amat,fd', isgn, col, var, err, a(isgn,col), -(rp(isgn)-rm(isgn))/(2*eps)
            error stop 1
          endif
        enddo
      enddo
    enddo
  enddo
  sheath_j_current_row = .false. ; base(:)%values(1,1,var_zj) = 0.d0 ; base(:)%values(1,1,var_u) = ufl
  write(*,'(a,es9.2)') ' PASS: current-slot form: zj row only, every column matches FD, worst rel err ', worst

  ! ---------------------------------------------------------------- 5. sign of the saturation current
  do isgn = 1, 2
    call sheath_j_norm(a_n, c_sat)
    if ( -c_sat/F0 .le. 0.d0 ) then
      write(*,*) 'FAIL: saturation current does not flow into the wall, F0 =', F0
      error stop 1
    endif
    F0 = -F0
  enddo
  write(*,'(a)') ' PASS: ion saturation current flows into the wall for both signs of F0'

  ! ---------------------------------------------------------------- 6. diagnostics
  call floating_diag_reset()
  floating_u_diag = .true.;  nodes = base; call assemble(ap, rp)
  floating_u_diag = .false.; nodes = base; call assemble(a, r)
  if ( any(a /= ap) .or. any(r /= rp) ) error stop 'FAIL: sheath diagnostics changed the equations'
  call floating_diag_report(0)
  write(*,'(a)') ' PASS: sheath diagnostics leave the equations untouched'

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
