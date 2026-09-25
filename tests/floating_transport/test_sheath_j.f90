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
!!  7. cancelled wall fluxes (sheath_j_cancel_flux, current-row form): value-DOF rows of u integrate to
!!     -visco*c*int R^3 dR - [R^2 p] + (cu^2/2) int R^2 d(R^2 rho) for w = c*Z, u = cu*R, rho = rho0 + rho1*R,
!!     uniform T (old setup: int R dR for the viscous part); every w/rho/Ti/Te/u column of the u rows (trace and
!!     normal-derivative DOFs) matches FD with a temperature-dependent viscosity; no u row with the flag off or
!!     without a sheath node.
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
  real*8  :: eps, err, worst, scale, a_n, C_T, C_V, c_sat, ufl, bslope, tot, cslope, pval, fdv, want
  integer, parameter :: cf_vars(5) = [var_w, var_rho, var_Ti, var_Te, var_u]
  real*8  :: rho0, rho1, cu
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
  write(*,'(a,es9.2)') ' PASS: current-slot form: zj row only, every column matches FD, worst rel err ', worst
  ! --- 7. cancelled wall fluxes (sheath_j_cancel_flux), current-row form still on
  ! --- (a) values on the bottom edge (outward normal -Z, s = +R): w = c*Z gives dw/dn = -c; u = cu*R gives
  ! ---     v_E^2 = R^2 cu^2; rho = rho0 + rho1*R, uniform T. Sum of the value-DOF rows of u (partition of unity,
  ! ---     tstep = 1):
  ! ---       viscous   -visco*c*int_1^2 R^3 dR = -15/4 visco c        (visco_old_setup: int R dR = 3/2)
  ! ---       magnetis. -[R^2 p]_1^2 = -(3 rho0 + 7 rho1) T
  ! ---       kinetic   +(cu^2/2) int_1^2 R^2 (2 R rho0 + 3 R^2 rho1) dR = cu^2 (15/4 rho0 + 93/10 rho1)
  sheath_j_cancel_flux = .true.
  cslope = 0.05d0 ; rho0 = 0.2d0 ; rho1 = 0.05d0 ; cu = 0.3d0*ufl ; pval = 0.003d0+0.004d0
  do i = 1, 4
    base(i)%values(1,1,var_w)   = cslope*base(i)%x(1,1,2)
    base(i)%values(1,3,var_w)   = cslope/3.d0
    base(i)%values(1,1,var_rho) = rho0 + rho1*base(i)%x(1,1,1)
    base(i)%values(1,2,var_rho) = rho1
    base(i)%values(1,1,var_u)   = cu*base(i)%x(1,1,1)
    base(i)%values(1,2,var_u)   = cu
  enddo
  nodes = base; call assemble(a, r)
  if ( .not. anyrow(a, r, var_u) ) error stop 'FAIL: cancelled fluxes assemble no u row'
  tot  = r(n_var*4*0 + var_u) + r(n_var*4*1 + var_u)
  want = -15.d0/4.d0*visco*cslope - (3.d0*rho0 + 7.d0*rho1)*pval + cu**2*(15.d0/4.d0*rho0 + 9.3d0*rho1)
  if ( abs(tot - want) > 1.d-9*abs(want) ) then
    write(*,*) 'FAIL: cancelled fluxes, default setup: value rows sum', tot, want
    error stop 1
  endif
  visco_old_setup = .true. ; nodes = base; call assemble(a, r)
  tot  = r(n_var*4*0 + var_u) + r(n_var*4*1 + var_u)
  want = -1.5d0*visco*cslope - (3.d0*rho0 + 7.d0*rho1)*pval + cu**2*(15.d0/4.d0*rho0 + 9.3d0*rho1)
  if ( abs(tot - want) > 1.d-9*abs(want) ) then
    write(*,*) 'FAIL: cancelled fluxes, old setup: value rows sum', tot, want
    error stop 1
  endif
  visco_old_setup = .false.
  do i = 1, 4
    base(i)%values(1,:,var_u) = 0.d0 ; base(i)%values(1,1,var_u) = 0.8d0*ufl
  enddo
  write(*,'(a)') ' PASS: cancelled wall fluxes: value rows integrate to the analytic viscous + magnetisation + kinetic sum'
  ! --- (b) FD of every w/rho/Ti/Te column of the u rows at a generic state, viscosity temperature dependent
  visco_T_dependent = .true. ; visco = 1.d-3
  do i = 1, 4
    base(i)%values(1,1,var_u)   = ufl*(0.8d0 + 0.3d0*base(i)%x(1,1,1) - 0.2d0*base(i)%x(1,1,2))
    base(i)%values(1,2,var_u)   = 0.3d0*ufl ; base(i)%values(1,3,var_u) = -0.2d0*ufl/3.d0 ; base(i)%values(1,4,var_u) = 0.05d0*ufl
    base(i)%values(1,1,var_w)   = 0.03d0 + 0.02d0*base(i)%x(1,1,1) + 0.05d0*base(i)%x(1,1,2)
    base(i)%values(1,2,var_w)   = 0.02d0 ; base(i)%values(1,3,var_w) = 0.05d0/3.d0 ; base(i)%values(1,4,var_w) = 0.01d0
    base(i)%values(1,1,var_rho) = 0.2d0 + 0.05d0*base(i)%x(1,1,1)
    base(i)%values(1,2,var_rho) = 0.05d0 ; base(i)%values(1,3,var_rho) = 0.01d0
    base(i)%values(1,1,var_Ti)  = 0.003d0 + 0.0005d0*base(i)%x(1,1,1)
    base(i)%values(1,2,var_Ti)  = 0.0005d0 ; base(i)%values(1,3,var_Ti) = 0.0002d0
    base(i)%values(1,1,var_Te)  = 0.004d0 + 0.0007d0*base(i)%x(1,1,1)
    base(i)%values(1,2,var_Te)  = 0.0007d0 ; base(i)%values(1,3,var_Te) = -0.0003d0
  enddo
  nodes = base; call assemble(a, r)
  worst = 0.d0
  do i = 1, size(cf_vars)
    var = cf_vars(i)
    tot = 0.d0                                ! non-empty columns of this variable
    do row = 1, 2
      do dof = 1, 4
        col = n_var*4*(row-1) + n_var*(dof-1) + var
        nodes = base; nodes(row)%values(1,dof,var) = nodes(row)%values(1,dof,var) + eps; call assemble(ap, rp)
        nodes = base; nodes(row)%values(1,dof,var) = nodes(row)%values(1,dof,var) - eps; call assemble(am, rm)
        scale = 0.d0
        do isgn = 1, nd
          if ( isvar(isgn, var_u) ) scale = max( scale, abs(a(isgn,col)), abs((rp(isgn)-rm(isgn))/(2*eps)) )
        enddo
        if ( scale == 0.d0 ) cycle            ! e.g. a value DOF of w: no normal derivative on an edge with n = -Z
        tot = tot + 1.d0
        do isgn = 1, nd
          if ( .not. isvar(isgn, var_u) ) cycle
          fdv = -(rp(isgn)-rm(isgn))/(2*eps)
          err = abs( a(isgn,col) - fdv ) / scale
          worst = max(worst, err)
          if ( err > 1.d-5 ) then
            write(*,'(a,4i5,3es12.3)') ' FAIL: cancelled fluxes FD var,node,dof,row,err,amat,fd', var, row, dof, isgn, err, a(isgn,col), fdv
            error stop 1
          endif
        enddo
      enddo
    enddo
    if ( tot < 3.d0 ) then
      write(*,'(a,i5,f5.0)') ' FAIL: cancelled fluxes: too few non-empty columns for var', var, tot
      error stop 1
    endif
  enddo
  write(*,'(a,es9.2)') ' PASS: cancelled wall fluxes: every w/rho/Ti/Te/u column of the u rows matches FD, worst rel err ', worst
  ! --- (c) gating: no u row without a sheath node or with the flag off
  bcs(1)%sheath_j = .false. ; nodes = base; call assemble(ap, rp)
  if ( anyrow(ap, rp, var_u) ) error stop 'FAIL: cancelled fluxes assembled without a sheath node'
  bcs(1)%sheath_j = .true.
  sheath_j_cancel_flux = .false. ; nodes = base; call assemble(ap, rp)
  if ( anyrow(ap, rp, var_u) ) error stop 'FAIL: cancelled fluxes assembled with the flag off'
  write(*,'(a)') ' PASS: cancelled wall fluxes: only on sheath edges and only with sheath_j_cancel_flux'
  visco_T_dependent = .false. ; visco = 1.d-5
  do i = 1, 4
    base(i)%values(1,:,var_u) = 0.d0 ; base(i)%values(1,1,var_u) = 0.8d0*ufl
    base(i)%values(1,:,var_w) = 0.d0
    base(i)%values(1,:,var_rho) = 0.d0 ; base(i)%values(1,1,var_rho) = 0.2d0
    base(i)%values(1,:,var_Ti)  = 0.d0 ; base(i)%values(1,1,var_Ti)  = 0.003d0
    base(i)%values(1,:,var_Te)  = 0.d0 ; base(i)%values(1,1,var_Te)  = 0.004d0
  enddo
  ! --- ion branch with a finite slope: x < 0 (u above floating), FD the u and Te columns
  sheath_j_ion_slope = 0.03d0 ; base(:)%values(1,1,var_u) = 2.d0*ufl
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
          if ( .not. isvar(isgn, var_zj) ) cycle
          err = abs( a(isgn,col) + (rp(isgn)-rm(isgn))/(2*eps) ) / max(1.d0, maxval(abs(a(:,col))))
          worst = max(worst, err)
          if ( err > 1.d-6 ) then
            write(*,'(a,3i5,3es12.3)') ' FAIL: ion-slope FD row,col,var,err,amat,fd', isgn, col, var, err, a(isgn,col), -(rp(isgn)-rm(isgn))/(2*eps)
            error stop 1
          endif
        enddo
      enddo
    enddo
  enddo
  sheath_j_ion_slope = 0.d0
  write(*,'(a,es9.2)') ' PASS: ion-branch slope: every column matches FD on the ion side, worst rel err ', worst
  ! --- beyond electron saturation with a finite slope: x > Lambda (u below the wall), FD the u and Te columns
  sheath_j_e_slope = 0.03d0 ; base(:)%values(1,1,var_u) = -3.d0*ufl
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
          if ( .not. isvar(isgn, var_zj) ) cycle
          err = abs( a(isgn,col) + (rp(isgn)-rm(isgn))/(2*eps) ) / max(1.d0, maxval(abs(a(:,col))))
          worst = max(worst, err)
          if ( err > 1.d-6 ) then
            write(*,'(a,3i5,3es12.3)') ' FAIL: e-slope FD row,col,var,err,amat,fd', isgn, col, var, err, a(isgn,col), -(rp(isgn)-rm(isgn))/(2*eps)
            error stop 1
          endif
        enddo
      enddo
    enddo
  enddo
  sheath_j_e_slope = 0.d0
  write(*,'(a,es9.2)') ' PASS: electron-saturation slope: every column matches FD beyond the cap, worst rel err ', worst
  sheath_j_current_row = .false. ; base(:)%values(1,1,var_zj) = 0.d0 ; base(:)%values(1,1,var_u) = ufl

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
