!> Gauss-point minima (mod_state_check): a cubic trace can be negative between nodes while every
!! nodal value is positive. rho = 1 at all four vertices with s-slopes -10 / +10 along the bottom
!! edge gives rho(0.5, 0) = -1.5; the check must see it, and must not with zero slopes.
program test_state_check
  use mod_parameters
  use data_structure
  use mod_state_check
  use basis_at_gaussian, only: set_basis
  implicit none
  type(type_element_list) :: el
  type(type_node_list)    :: nl
  real*8  :: minima(3), at(2,3)
  integer :: i

  call set_basis()
  nl%n_nodes = 4; nl%n_dof = 16
  allocate( nl%node(4) )
  el%n_elements = 1
  el%element(1)%vertex = [1,2,3,4]
  ! unit square: 1:(0,0) 2:(1,0) 3:(1,1) 4:(0,1), with the s and t derivative DOFs of the geometry
  nl%node(1)%x(1,1,:) = [0.d0,0.d0]; nl%node(2)%x(1,1,:) = [1.d0,0.d0]
  nl%node(3)%x(1,1,:) = [1.d0,1.d0]; nl%node(4)%x(1,1,:) = [0.d0,1.d0]
  do i = 1, 4
    nl%node(i)%x(1,2,:) = [1.d0,0.d0]; nl%node(i)%x(1,3,:) = [0.d0,1.d0]
  enddo
  do i = 1, 4
    nl%node(i)%values(1,1,var_rho) = 1.d0
    nl%node(i)%values(1,1,var_Ti)  = 2.d0
    nl%node(i)%values(1,1,var_Te)  = 3.d0
  enddo

  call state_gauss_minima(0, nl, el, [var_rho, var_Ti, var_Te], minima)
  if ( any(abs(minima - [1.d0, 2.d0, 3.d0]) > 1.d-12) ) then
    write(*,*) 'FAIL: constant fields', minima
    error stop 1
  endif

  nl%node(1)%values(1,2,var_rho) = -10.d0
  nl%node(2)%values(1,2,var_rho) = +10.d0
  call state_gauss_minima(0, nl, el, [var_rho, var_Ti, var_Te], minima, at)
  ! the undershoot sits on the bottom edge (t = 0) near its middle
  if ( abs(at(2,1)) > 0.1d0 .or. abs(at(1,1)-0.5d0) > 0.2d0 ) then
    write(*,*) 'FAIL: undershoot location', at(:,1)
    error stop 1
  endif
  if ( minima(1) >= 0.d0 .or. abs(minima(2)-2.d0) > 1.d-12 ) then
    write(*,*) 'FAIL: between-node undershoot not detected', minima
    error stop 1
  endif
  write(*,'(a,es10.2)') ' PASS: Gauss-point minima catch a between-node undershoot, min rho = ', minima(1)
end program
