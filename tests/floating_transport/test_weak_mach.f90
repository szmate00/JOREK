!> Serial checks of the weak Bohm row (mach1_weak) in the PRODUCTION assembler
!! mod_boundary_matrix_open.f90, compiled against tests/floating_transport/fixtures.f90.
!!
!!  1. mach1_weak off: no Vpar row is assembled. Switched on with u constant along the wall, every
!!     other row keeps its residual and its columns, except for the new u columns of the sheath
!!     fluxes (a sensitivity that exists even where vE.n = 0).
!!  2. Drift-compensating form (mach1_weak_drift): every column of the Vpar rows (u, Vpar, Ti, Te,
!!     rho) matches finite differences; the default marginal form has no u column and matches too.
!!     Every column of the Vpar rows matches a
!!     central finite difference of the residual. Removing any column makes this fail.
!!  3. Saturated branch (drift alone beyond sonic outflow): the u and temperature columns are
!!     exactly zero; only the Vpar column survives.
!!  4. B.n crossing zero along the wall: the row fades continuously with the incidence, with no
!!     jump between the two signs.
!!  5. Inflow closure on the density row: absent while the total normal flow is outward; with
!!     inward flow every column (rho, u, Vpar, Ti, Te) of the rho rows matches finite differences.
!!  6. Sheath energy fluxes on the total outgoing flow: with a tangential u slope every column of
!!     the Ti and Te rows matches finite differences; with inward total flow they collect nothing
!!     beyond the grazing-incidence floor.
program test_weak_mach
  use mod_parameters
  use phys_module
  use data_structure
  use mod_boundary_matrix_open
  use basis_at_gaussian, only: set_basis
  use mod_floating_diag, only: floating_diag_reset, floating_diag_report
  implicit none
  integer, parameter :: nd = 4*4*n_var
  integer, parameter :: fd_vars(5) = [var_u, var_vpar, var_Ti, var_Te, var_rho]
  type(type_element) :: e
  type(type_node)    :: nodes(4), base(4)
  real*8  :: a0(nd,nd), r0(nd), a(nd,nd), r(nd), ap(nd,nd), rp(nd), am(nd,nd), rm(nd)
  real*8  :: eps, err, worst, scale, s, rowmax(7)
  integer :: i, k, var, dof, col, row, isl

  call set_basis()

  ! --- Side 1 of a unit-ish element: R from 1 to 2 at Z = 0, exterior normal pointing to -Z.
  base(1)%x(1,1,:) = [1.d0, 0.d0]; base(2)%x(1,1,:) = [2.d0, 0.d0]
  base(3)%x(1,1,:) = [2.d0, 1.d0]; base(4)%x(1,1,:) = [1.d0, 1.d0]
  do i = 1, 4
    base(i)%x(1,2,:) = [1.d0, 0.d0]
    base(i)%x(1,3,:) = [0.d0, 1.d0/3.d0]
    base(i)%values(1,1,var_rho)  = 0.2d0
    base(i)%values(1,1,var_Ti)   = 0.003d0
    base(i)%values(1,1,var_Te)   = 0.004d0
    base(i)%values(1,1,var_vpar) = 0.04d0
    base(i)%values(1,1,var_psi)  = 0.08d0*base(i)%x(1,1,1)   ! B_pol = -0.08/R e_Z: finite incidence
    base(i)%values(1,2,var_psi)  = 0.08d0
    base(i)%values(1,3,var_psi)  = 0.02d0                    ! free normal derivative, enters |B| only
  enddo
  base%boundary = 1

  ! ---------------------------------------------------------------- 1. off / on
  ! u constant along the wall: vE.n = 0, so with outward flow every row but Vpar's is unchanged
  call set_u(0.d0)
  mach1_weak = .false.
  nodes = base; call assemble(a0, r0)
  do row = 1, nd
    if ( isvar(row, var_vpar) .and. any(a0(row,:) /= 0.d0) ) error stop 'FAIL: Vpar row assembled with mach1_weak off'
  enddo
  mach1_weak = .true.; mach1_weak_drift = .true.   ! drift-compensating form first, its u column is the hard one
  nodes = base; call assemble(a, r)
  ! (to roundoff: the total-flow measure max(vn,0)*R*dl equals vpar0*ps0_s*normal_sign3 analytically;
  !  the u columns are new and are checked against finite differences in 6)
  do row = 1, nd
    if ( isvar(row, var_vpar) ) cycle
    scale = max(1.d0, maxval(abs(a0(row,:))), abs(r0(row)))
    if ( abs(r(row) - r0(row)) > 1.d-12*scale ) error stop 'FAIL: mach1_weak changed a residual'
    do col = 1, nd
      if ( isvar(col, var_u) ) cycle
      if ( abs(a(row,col) - a0(row,col)) > 1.d-12*scale ) error stop 'FAIL: mach1_weak changed a non-u column'
    enddo
  enddo
  worst = 0.d0
  do row = 1, nd
    if ( isvar(row, var_vpar) ) worst = max(worst, maxval(abs(a(row,:))))
  enddo
  if ( worst <= 0.d0 ) error stop 'FAIL: weak row not assembled'
  write(*,'(a)') ' PASS: mach1_weak off leaves no Vpar row; on keeps every other residual and column but u'
  call set_u(1.d-4)      ! small tangential u slope: inside the compensating branch, away from the kink
  nodes = base; call assemble(a, r)

  ! ---------------------------------------------------------------- 2. FD every column
  ! rhs = -w*res and amat = +w*d(res)/dx, so amat + d(rhs)/dx must vanish. psi carries no column:
  ! B.n depends only on its Dirichlet trace DOFs, and the |B| dependence on the free normal
  ! derivative is lagged (the assembler's column loop covers trace DOFs only). Its size relative
  ! to the Vpar column is reported below, not asserted.
  eps = 1.d-8
  worst = 0.d0
  do k = 1, size(fd_vars)
    var = fd_vars(k)
    do i = 1, 2
      do dof = 1, 4
        col = n_var*4*(i-1) + n_var*(dof-1) + var
        nodes = base; nodes(i)%values(1,dof,var) = nodes(i)%values(1,dof,var) + eps; call assemble(ap, rp)
        nodes = base; nodes(i)%values(1,dof,var) = nodes(i)%values(1,dof,var) - eps; call assemble(am, rm)
        scale = max(1.d0, maxval(abs(a(:,col))))
        do row = 1, nd
          if ( .not. isvar(row, var_vpar) ) cycle
          err = abs( a(row,col) + (rp(row)-rm(row))/(2*eps) ) / scale
          worst = max(worst, err)
          if ( err > 1.d-6 ) then
            write(*,'(a,3i5,3es12.3)') ' FAIL: weak row FD row,col,var,err,amat,fd', row, col, var, err, a(row,col), -(rp(row)-rm(row))/(2*eps)
            error stop 1
          endif
        enddo
      enddo
    enddo
  enddo
  write(*,'(a,es9.2)') ' PASS: every column of the weak row matches finite differences, worst rel err ', worst
  ! --- lagged |B| dependence: FD of the residual in the free normal psi derivative vs the Vpar column
  nodes = base; nodes(1)%values(1,3,var_psi) = nodes(1)%values(1,3,var_psi) + eps; call assemble(ap, rp)
  nodes = base; nodes(1)%values(1,3,var_psi) = nodes(1)%values(1,3,var_psi) - eps; call assemble(am, rm)
  scale = 0.d0; err = 0.d0
  do row = 1, nd
    if ( .not. isvar(row, var_vpar) ) cycle
    err   = max(err,   abs((rp(row)-rm(row))/(2*eps)))
    do i = 1, 2
      do dof = 1, 4
        scale = max(scale, abs(a(row, n_var*4*(i-1)+n_var*(dof-1)+var_vpar)))
      enddo
    enddo
  enddo
  write(*,'(a,es9.2)') ' INFO: lagged psi (normal derivative) sensitivity relative to the Vpar column: ', err/scale

  ! --- marginal form (default): no u column at all, every other column still matches FD
  mach1_weak_drift = .false.
  nodes = base; call assemble(a, r)
  do row = 1, nd
    if ( .not. isvar(row, var_vpar) ) cycle
    do col = 1, nd
      if ( isvar(col, var_u) .and. a(row,col) /= 0.d0 ) error stop 'FAIL: marginal weak row carries a u column'
    enddo
  enddo
  do k = 1, size(fd_vars)
    var = fd_vars(k)
    do i = 1, 2
      do dof = 1, 4
        col = n_var*4*(i-1) + n_var*(dof-1) + var
        nodes = base; nodes(i)%values(1,dof,var) = nodes(i)%values(1,dof,var) + eps; call assemble(ap, rp)
        nodes = base; nodes(i)%values(1,dof,var) = nodes(i)%values(1,dof,var) - eps; call assemble(am, rm)
        scale = max(1.d0, maxval(abs(a(:,col))))
        do row = 1, nd
          if ( .not. isvar(row, var_vpar) ) cycle
          if ( abs( a(row,col) + (rp(row)-rm(row))/(2*eps) ) / scale > 1.d-6 ) error stop 'FAIL: marginal weak row FD'
        enddo
      enddo
    enddo
  enddo
  write(*,'(a)') ' PASS: marginal weak row (default) has no u column and matches FD'
  mach1_weak_drift = .true.

  ! ---------------------------------------------------------------- 3. saturated branch
  call set_u(-1.d0)      ! vE.n far beyond sonic outflow: target pinned at zero
  nodes = base; call assemble(a, r)
  do row = 1, nd
    if ( .not. isvar(row, var_vpar) ) cycle
    do col = 1, nd
      if ( isvar(col, var_vpar) ) cycle
      if ( a(row,col) /= 0.d0 ) then
        write(*,'(a,2i5,es12.3)') ' FAIL: saturated weak row keeps a column', row, col, a(row,col)
        error stop 1
      endif
    enddo
  enddo
  write(*,'(a)') ' PASS: saturated branch keeps only the Vpar column'

  ! ---------------------------------------------------------------- 4. B.n through zero
  call set_u(1.d-4)
  do isl = -3, 3
    s = 0.02d0*isl
    do i = 1, 4
      base(i)%values(1,1,var_psi) = s*base(i)%x(1,1,1)
      base(i)%values(1,2,var_psi) = s
    enddo
    nodes = base; call assemble(a, r)
    rowmax(isl+4) = 0.d0
    do row = 1, nd
      if ( isvar(row, var_vpar) ) rowmax(isl+4) = max(rowmax(isl+4), maxval(abs(a(row,:))), abs(r(row)))
    enddo
  enddo
  ! symmetric in the sign of B.n, monotone in |B.n|, and vanishing at tangency
  do isl = 1, 3
    if ( abs(rowmax(4+isl) - rowmax(4-isl)) > 1.d-9*rowmax(4+isl) ) error stop 'FAIL: row not symmetric in sign(B.n)'
    if ( rowmax(4+isl) <= rowmax(4+isl-1) ) error stop 'FAIL: row not monotone in |B.n|'
  enddo
  if ( rowmax(4) > 1.d-12*rowmax(7) ) error stop 'FAIL: row does not vanish at tangency'
  write(*,'(a)') ' PASS: weak row is even in B.n, monotone in |B.n| and vanishes at tangency'

  ! ---------------------------------------------------------------- 5. inflow closure
  do i = 1, 4
    base(i)%values(1,1,var_psi) = 0.08d0*base(i)%x(1,1,1)
    base(i)%values(1,2,var_psi) = 0.08d0
  enddo
  ! outward total flow: the rho rows must be exactly those of the mach1_weak-off assembly
  mach1_weak = .false.; nodes = base; call assemble(a0, r0)
  mach1_weak = .true.;  nodes = base; call assemble(a, r)
  do row = 1, nd
    if ( .not. isvar(row, var_rho) ) cycle
    if ( any(a(row,:) /= a0(row,:)) .or. r(row) /= r0(row) ) error stop 'FAIL: inflow term active with outward flow'
  enddo
  ! inward total flow: reverse Vpar so Vpar*(B.n) < 0 with a small drift
  base(:)%values(1,1,var_vpar) = -0.04d0
  nodes = base; call assemble(a, r)
  worst = 0.d0
  do k = 1, size(fd_vars)
    var = fd_vars(k)
    do i = 1, 2
      do dof = 1, 4
        col = n_var*4*(i-1) + n_var*(dof-1) + var
        nodes = base; nodes(i)%values(1,dof,var) = nodes(i)%values(1,dof,var) + eps; call assemble(ap, rp)
        nodes = base; nodes(i)%values(1,dof,var) = nodes(i)%values(1,dof,var) - eps; call assemble(am, rm)
        scale = max(1.d0, maxval(abs(a(:,col))))
        do row = 1, nd
          if ( .not. isvar(row, var_rho) ) cycle
          err = abs( a(row,col) + (rp(row)-rm(row))/(2*eps) ) / scale
          worst = max(worst, err)
          if ( err > 1.d-6 ) then
            write(*,'(a,3i5,3es12.3)') ' FAIL: inflow FD row,col,var,err,amat,fd', row, col, var, err, a(row,col), -(rp(row)-rm(row))/(2*eps)
            error stop 1
          endif
        enddo
      enddo
    enddo
  enddo
  nodes = base; call assemble(a0, r0)
  if ( all(r0(pack([(row, row=1,nd)], [(isvar(row,var_rho), row=1,nd)])) == 0.d0) ) error stop 'FAIL: inflow term not assembled'
  write(*,'(a,es9.2)') ' PASS: inflow closure: inactive for outward flow, every rho-row column matches FD, worst rel err ', worst

  ! ---------------------------------------------------------------- 6. energy fluxes on the total flow
  base(:)%values(1,1,var_vpar) = 0.04d0
  call set_u(1.d-4)
  nodes = base; call assemble(a, r)
  worst = 0.d0
  do k = 1, size(fd_vars)
    var = fd_vars(k)
    do i = 1, 2
      do dof = 1, 4
        col = n_var*4*(i-1) + n_var*(dof-1) + var
        nodes = base; nodes(i)%values(1,dof,var) = nodes(i)%values(1,dof,var) + eps; call assemble(ap, rp)
        nodes = base; nodes(i)%values(1,dof,var) = nodes(i)%values(1,dof,var) - eps; call assemble(am, rm)
        scale = max(1.d0, maxval(abs(a(:,col))))
        do row = 1, nd
          if ( .not. (isvar(row, var_Ti) .or. isvar(row, var_Te)) ) cycle
          err = abs( a(row,col) + (rp(row)-rm(row))/(2*eps) ) / scale
          worst = max(worst, err)
          if ( err > 1.d-6 ) then
            write(*,'(a,3i5,3es12.3)') ' FAIL: energy-row FD row,col,var,err,amat,fd', row, col, var, err, a(row,col), -(rp(row)-rm(row))/(2*eps)
            error stop 1
          endif
        enddo
      enddo
    enddo
  enddo
  ! the u column must actually be there
  err = 0.d0
  do row = 1, nd
    if ( .not. isvar(row, var_Te) ) cycle
    do i = 1, 2
      do dof = 1, 4
        err = max(err, abs(a(row, n_var*4*(i-1)+n_var*(dof-1)+var_u)))
      enddo
    enddo
  enddo
  if ( err <= 0.d0 ) error stop 'FAIL: energy row has no u column'
  ! inward total flow: the sheath collection is gone, only the c_angle floor remains, and the
  ! Te row must then be independent of u and Vpar
  base(:)%values(1,1,var_vpar) = -0.04d0
  nodes = base; call assemble(a, r)
  do row = 1, nd
    if ( .not. isvar(row, var_Te) ) cycle
    do i = 1, 2
      do dof = 1, 4
        if ( a(row, n_var*4*(i-1)+n_var*(dof-1)+var_u)    /= 0.d0 ) error stop 'FAIL: closed wall keeps a u column'
        if ( a(row, n_var*4*(i-1)+n_var*(dof-1)+var_vpar) /= 0.d0 ) error stop 'FAIL: closed wall keeps a Vpar column'
      enddo
    enddo
  enddo
  write(*,'(a,es9.2)') ' PASS: sheath energy fluxes on the total flow: Ti/Te rows match FD, closed wall collects nothing, worst rel err ', worst

  ! ---------------------------------------------------------------- 7. diagnostics change nothing
  call floating_diag_reset()
  floating_u_diag = .true.;  nodes = base; call assemble(ap, rp)
  floating_u_diag = .false.; nodes = base; call assemble(a, r)
  if ( any(a /= ap) .or. any(r /= rp) ) error stop 'FAIL: floating_u_diag changed the equations'
  call floating_diag_report(0)
  write(*,'(a)') ' PASS: wall diagnostics leave the equations untouched'

contains

  logical function isvar(idx, var)
    integer, intent(in) :: idx, var
    isvar = ( mod(idx-1, n_var) + 1 == var )
  end function

  !> u exactly linear along the wall with tangential slope g: u = 0.01 + g*(R-1).
  subroutine set_u(g)
    real*8, intent(in) :: g
    integer :: ii
    do ii = 1, 4
      base(ii)%values(1,1,var_u) = 0.01d0 + g*(base(ii)%x(1,1,1) - 1.d0)
      base(ii)%values(1,2,var_u) = g
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

end program test_weak_mach
