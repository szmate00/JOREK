!> Serial checks of the sheath-set wall flux on the total normal flow (bcs%floating_u + mach1 types):
!!  A. u constant along the wall, flow above the floor: rho/Ti/Te residuals equal develop (parallel-measure identity),
!!     columns equal except the new u columns and the sf_on-gated cross-species floor columns.
!!  B. every u/Vpar/Ti/Te/rho column of the rho, Ti, Te rows vs central FD in five states (sub-Bohm parallel flow,
!!     ExB inflow, ExB outflow, near the branch point, vpar_smoothing); psi trace columns to the lagged-|B| level.
program test_nodal_flux
  use mod_parameters
  use phys_module
  use data_structure
  use mod_boundary_matrix_open
  use basis_at_gaussian, only: set_basis
  implicit none
  integer, parameter :: nd = 4*4*n_var
  integer, parameter :: fd_vars(5) = [var_u, var_vpar, var_Ti, var_Te, var_rho]
  integer, parameter :: rows(3) = [var_rho, var_Ti, var_Te]
  type(type_element) :: e
  type(type_node)    :: nodes(4), base(4)
  real*8  :: a0(nd,nd), r0(nd), a(nd,nd), r(nd), ap(nd,nd), rp(nd), am(nd,nd), rm(nd)
  real*8  :: eps, worst, wpsi, scale
  integer :: i, row, col, icase
  logical :: cross
  call set_basis()
  base(1)%x(1,1,:) = [1.d0, 0.d0]; base(2)%x(1,1,:) = [2.d0, 0.d0]
  base(3)%x(1,1,:) = [2.d0, 1.d0]; base(4)%x(1,1,:) = [1.d0, 1.d0]
  do i = 1, 4
    base(i)%x(1,2,:) = [1.d0, 0.d0]; base(i)%x(1,3,:) = [0.d0, 1.d0/3.d0]
    base(i)%values(1,1,var_rho) = 0.2d0 ; base(i)%values(1,2,var_rho) = 0.03d0
    base(i)%values(1,1,var_Ti)  = 0.003d0 ; base(i)%values(1,2,var_Ti) = 0.0004d0
    base(i)%values(1,1,var_Te)  = 0.004d0
    base(i)%values(1,1,var_psi) = 0.08d0*base(i)%x(1,1,1) ; base(i)%values(1,2,var_psi) = 0.08d0 ; base(i)%values(1,3,var_psi) = 0.02d0
  enddo
  base%boundary = 1
  eps = 1.d-8
  density_reflection = 0.3d0
  call set_state(0.2d0, 0.d0)
  bcs(1)%floating_u = .false.; nodes = base; call assemble(a0, r0)
  bcs(1)%floating_u = .true.;  nodes = base; call assemble(a, r)
  do row = 1, nd
    if ( .not. isrow(row) ) cycle
    if ( abs(r(row) - r0(row)) > 1.d-12*max(abs(r0(row)),1.d-30) ) then
      write(*,*) 'FAIL A: residual differs', row, r(row), r0(row); error stop 1
    endif
    do col = 1, nd
      cross = (varof(row) == var_Ti .and. varof(col) == var_Te) .or. (varof(row) == var_Te .and. varof(col) == var_Ti)
      if ( varof(col) == var_u .or. cross ) cycle
      if ( abs(a(row,col) - a0(row,col)) > 1.d-12*max(abs(a0(row,col)),1.d-30) ) then
        write(*,*) 'FAIL A: column differs', row, col, a(row,col), a0(row,col); error stop 1
      endif
    enddo
  enddo
  write(*,'(a)') ' PASS A: u constant, flow above the floor: residuals and columns equal develop (new u / cross columns aside)'
  do icase = 1, 5
    vpar_smoothing = .false.
    select case (icase)
    case (1); call set_state( 0.01d0,  0.d0)
    case (2); call set_state( 0.2d0,  +0.05d0)
    case (3); call set_state( 0.2d0,  -0.05d0)
    case (4); call set_state( 0.01d0, +0.002d0)
    case (5); call set_state( 0.01d0,  0.d0); vpar_smoothing = .true.; vpar_smoothing_coef = [0.05d0, 0.02d0, 0.d0]
    end select
    nodes = base; call assemble(a, r)
    worst = 0.d0; wpsi = 0.d0
    do col = 1, nd
      if ( .not. any(varof(col) == [fd_vars, var_psi]) ) cycle
      if ( varof(col) == var_psi .and. dofof(col) > 2 ) cycle
      nodes = base; call bump(col, +eps); call assemble(ap, rp)
      nodes = base; call bump(col, -eps); call assemble(am, rm)
      scale = max(1.d-30, maxval(abs(a(:,col))), maxval(abs(rp-rm))/(2*eps))
      do row = 1, nd
        if ( .not. isrow(row) ) cycle
        if ( varof(col) == var_psi ) then
          wpsi = max(wpsi, abs(a(row,col) + (rp(row)-rm(row))/(2*eps)) / scale)
        else
          worst = max(worst, abs(a(row,col) + (rp(row)-rm(row))/(2*eps)) / scale)
          if ( abs(a(row,col) + (rp(row)-rm(row))/(2*eps)) / scale > 1.d-6 ) then
            write(*,'(a,4i5,2es14.5)') ' FAIL B: case,row,col,var, amat, -fd', icase, row, col, varof(col), a(row,col), -(rp(row)-rm(row))/(2*eps)
            error stop 1
          endif
        endif
      enddo
    enddo
    if ( wpsi > 2.d-3 ) error stop 'FAIL B: psi trace column beyond the lagged-Btot level'
    write(*,'(a,i2,a,es9.2,a,es9.2)') ' PASS B case', icase, ': u/Vpar/Ti/Te/rho columns vs FD worst ', worst, ', psi trace ', wpsi
  enddo
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
  logical function isrow(idx)
    integer, intent(in) :: idx
    isrow = any(varof(idx) == rows)
  end function
  subroutine bump(idx, d)
    integer, intent(in) :: idx
    real*8,  intent(in) :: d
    nodes(nodeof(idx))%values(1,dofof(idx),varof(idx)) = nodes(nodeof(idx))%values(1,dofof(idx),varof(idx)) + d
  end subroutine
  subroutine set_state(vp, g)
    real*8, intent(in) :: vp, g
    integer :: ii
    do ii = 1, 4
      base(ii)%values(1,1,var_vpar) = vp
      base(ii)%values(1,1,var_u)    = 0.01d0 + g*(base(ii)%x(1,1,1) - 1.d0)
      base(ii)%values(1,2,var_u)    = g
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
end program test_nodal_flux
