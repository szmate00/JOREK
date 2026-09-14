!> Gauss-point minima of selected variables over the whole domain, for admissible_update.
!! Real-space values on every toroidal plane, reduced over MPI ranks. The elements are split
!! across ranks the same way as in Integrals_3D.
module mod_state_check

  implicit none
  private

  public :: state_gauss_minima

contains

subroutine state_gauss_minima(my_id, node_list, element_list, vars, minima, where)

  use mod_parameters,     only: n_var, n_tor, n_plane, n_vertex_max, n_degrees
  use data_structure,     only: type_node_list, type_element_list, type_element, type_node, make_deep_copy_node
  use gauss,              only: n_gauss
  use basis_at_gaussian,  only: H, HZ
  use mpi_mod

  implicit none
  integer,                  intent(in)  :: my_id
  type(type_node_list),     intent(in)  :: node_list
  type(type_element_list),  intent(in)  :: element_list
  integer,                  intent(in)  :: vars(:)        !< variable indices to check (0 entries are skipped)
  real*8,                   intent(out) :: minima(size(vars))
  real*8, optional,         intent(out) :: where(2, size(vars))   !< (R,Z) of each minimum

  type(type_element) :: element
  type(type_node)    :: nodes(n_vertex_max)
  real*8  :: eq_g(n_plane, n_gauss, n_gauss), local(size(vars)), x_g(n_gauss, n_gauss), y_g(n_gauss, n_gauss)
  real*8  :: loc(2, size(vars)), loc_g(2, size(vars))
  integer :: imin(3)
  integer :: n_mpi, ife_delta, ife_min, ife_max, ife, iv, i, j, k, kv, in, mp, ms, mt, ierr

  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_mpi, ierr)
  n_mpi     = max(n_mpi, 1)
  ife_delta = ceiling(dble(element_list%n_elements) / n_mpi)
  ife_min   =      my_id     * ife_delta + 1
  ife_max   = min((my_id +1) * ife_delta, element_list%n_elements)

  local = huge(1.d0)
  loc   = -huge(1.d0)

  do ife = ife_min, ife_max
    element = element_list%element(ife)
    do iv = 1, n_vertex_max
      call make_deep_copy_node(node_list%node(element%vertex(iv)), nodes(iv))
    enddo

    x_g = 0.d0 ; y_g = 0.d0
    do i = 1, n_vertex_max
      do j = 1, n_degrees
        x_g = x_g + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,:,:)
        y_g = y_g + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,:,:)
      enddo
    enddo

    do kv = 1, size(vars)
      k = vars(kv)
      if ( k .le. 0 .or. k .gt. n_var ) cycle
      eq_g = 0.d0
      do i = 1, n_vertex_max
        do j = 1, n_degrees
          do mp = 1, n_plane
            do ms = 1, n_gauss
              do mt = 1, n_gauss
                do in = 1, n_tor
                  eq_g(mp,ms,mt) = eq_g(mp,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt) * HZ(in,mp)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      if ( minval(eq_g) .lt. local(kv) ) then
        local(kv) = minval(eq_g)
        imin      = minloc(eq_g)
        loc(:,kv) = (/ x_g(imin(2),imin(3)), y_g(imin(2),imin(3)) /)
      endif
    enddo
  enddo

  call MPI_ALLREDUCE(local, minima, size(vars), MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

  ! --- location: the rank holding the global minimum sends its (R,Z), the others -huge
  if ( present(where) ) then
    do kv = 1, size(vars)
      if ( local(kv) .ne. minima(kv) ) loc(:,kv) = -huge(1.d0)
    enddo
    call MPI_ALLREDUCE(loc, loc_g, 2*size(vars), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    where = loc_g
  endif

end subroutine state_gauss_minima

end module mod_state_check
