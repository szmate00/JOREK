!> Dumps the assembled edge matrix and residual for two states with floating_u OFF; run.sh builds it against the
!! current routine and against develop's and compares (expected: same sparsity, <1e-13 relative, re-associated products).
program dump_dev
  use mod_parameters
  use phys_module
  use data_structure
  use mod_boundary_matrix_open
  use basis_at_gaussian, only: set_basis
  implicit none
  integer, parameter :: nd = 4*4*n_var
  type(type_element) :: e
  type(type_node)    :: nodes(4)
  real*8  :: a(nd,nd), r(nd)
  integer :: i, k, vertices(2), directions(2)
  call set_basis()
  nodes(1)%x(1,1,:) = [1.d0, 0.d0]; nodes(2)%x(1,1,:) = [2.d0, 0.d0]
  nodes(3)%x(1,1,:) = [2.d0, 1.d0]; nodes(4)%x(1,1,:) = [1.d0, 1.d0]
  do i = 1, 4
    nodes(i)%x(1,2,:) = [1.d0, 0.d0]; nodes(i)%x(1,3,:) = [0.d0, 1.d0/3.d0]
    nodes(i)%values(1,1,var_rho) = 0.2d0 ; nodes(i)%values(1,2,var_rho) = 0.03d0
    nodes(i)%values(1,1,var_Ti)  = 0.003d0 ; nodes(i)%values(1,2,var_Ti) = 0.0004d0
    nodes(i)%values(1,1,var_Te)  = 0.004d0
    nodes(i)%values(1,1,var_psi) = 0.08d0*nodes(i)%x(1,1,1) ; nodes(i)%values(1,2,var_psi) = 0.08d0 ; nodes(i)%values(1,3,var_psi) = 0.02d0
    nodes(i)%values(1,1,var_u)   = 0.01d0 + 0.05d0*(nodes(i)%x(1,1,1)-1.d0) ; nodes(i)%values(1,2,var_u) = 0.05d0
  enddo
  nodes%boundary = 1
  density_reflection = 0.3d0
  open(10, file='dump.bin', form='unformatted', access='stream', status='replace')
  do k = 1, 2
    nodes(:)%values(1,1,var_vpar) = merge(0.2d0, 0.01d0, k == 1)
    vertices = [1,2]; directions = [1,2]; a = 0.d0; r = 0.d0
    call boundary_matrix_open(vertices, directions, e, nodes, .true., 1, 1.d0, 0.d0, 0.d0, 1.d0, [1.d0,1.d0], [0.d0,0.d0], a, r, 1, 1)
    write(10) a, r
  enddo
  close(10)
end program
