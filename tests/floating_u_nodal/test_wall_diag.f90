!> Smoke test of the [floating_u]/[mach1] wall diagnostics tables on one fixture edge (vE.n = -R*g on this edge).
program test_wd
  use mod_parameters
  use phys_module
  use data_structure
  use mod_boundary_matrix_open
  use mod_wall_diag
  use basis_at_gaussian, only: set_basis
  implicit none
  integer, parameter :: nd = 4*4*n_var
  type(type_element) :: e
  type(type_node)    :: nodes(4)
  type(type_node_list) :: nl
  real*8  :: a(nd,nd), r(nd), g
  integer :: i, vertices(2), directions(2)
  call set_basis()
  nodes(1)%x(1,1,:) = [1.d0, 0.d0]; nodes(2)%x(1,1,:) = [2.d0, 0.d0]
  nodes(3)%x(1,1,:) = [2.d0, 1.d0]; nodes(4)%x(1,1,:) = [1.d0, 1.d0]
  g = 5.d-2
  do i = 1, 4
    nodes(i)%x(1,2,:) = [1.d0, 0.d0]; nodes(i)%x(1,3,:) = [0.d0, 1.d0/3.d0]
    nodes(i)%values(1,1,var_rho) = 0.2d0
    nodes(i)%values(1,1,var_Ti)  = 0.003d0
    nodes(i)%values(1,1,var_Te)  = 0.004d0
    nodes(i)%values(1,1,var_psi) = 0.08d0*nodes(i)%x(1,1,1) ; nodes(i)%values(1,2,var_psi) = 0.08d0
    nodes(i)%values(1,1,var_vpar)= 0.05d0
    nodes(i)%values(1,1,var_u)   = 0.01d0 + g*(nodes(i)%x(1,1,1) - 1.d0) ; nodes(i)%values(1,2,var_u) = g
    nodes(i)%index(1) = i
  enddo
  nodes%boundary = 1
  bcs(1)%floating_u = .true.
  nl%n_nodes = 4 ; allocate(nl%node(4)) ; nl%node = nodes
  index_now = 3 ; wall_diag_every = 1 ; wall_diag_profile_every = 0
  call wall_diag_reset()
  vertices = [1,2]; directions = [1,2]; a = 0.d0; r = 0.d0
  call boundary_matrix_open(vertices, directions, e, nodes, .true., 1, 1.d0, 0.d0, 0.d0, 1.d0, [1.d0,1.d0], [0.d0,0.d0], a, r, 1, 1)
  call wall_diag_report(0, nl)
end program
