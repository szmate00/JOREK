!> Exterior element sides from mesh connectivity (model600).
!!
!! The open-boundary integral (mod_boundary_matrix_open) is applied by construct_matrix to every
!! element side whose two endpoint nodes carry a boundary label. A side can satisfy that and still
!! be interior, e.g. across a region one element wide between two wall segments. A boundary
!! integral on such a side is wrong. This module identifies the genuinely exterior sides: a side is
!! exterior if and only if exactly one element owns it.
!!
!! Nodes are identified by the global index of their value DOF, so coincident node records (two
!! records at one geometric vertex) count as the same vertex. A side whose two endpoints coincide
!! (a collapsed side on the axis) is never exterior.
module mod_boundary_edges

  implicit none
  private

  public :: boundary_edges_build, boundary_edges_active, boundary_edge_is_exterior

  logical, allocatable, save :: exterior(:,:)   !< exterior(side, element)

contains


!> Build the exterior-side table. Called once per matrix construction; cheap (one pass over the
!! local element list).
subroutine boundary_edges_build(element_list, node_list, my_id)

  use data_structure, only: type_element_list, type_node_list
  use mpi_mod

  implicit none
  type(type_element_list), intent(in) :: element_list
  type(type_node_list),    intent(in) :: node_list
  integer,                 intent(in) :: my_id

  integer, allocatable :: first(:), hi(:), elm(:), sid(:)
  integer :: ne, nid, e, s, lo, up, i, j, ierr, n_dup

  ne  = element_list%n_elements
  nid = node_list%n_dof
  if ( allocated(exterior) ) deallocate(exterior)
  allocate( exterior(4, ne) )
  exterior = .false.

  ! --- Sides keyed by their lower endpoint id, CSR layout: count, prefix sum, fill
  allocate( first(nid+2), hi(4*ne), elm(4*ne), sid(4*ne) )
  first = 0
  do e = 1, ne
    if ( element_list%element(e)%n_sons .gt. 0 ) then
      write(*,*) 'ERROR: boundary_edges_build needs an unrefined conforming mesh, EXITING!'
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    endif
    do s = 1, 4
      call side_ids(e, s, lo, up)
      if ( lo .eq. up ) cycle
      first(lo+2) = first(lo+2) + 1
    enddo
  enddo
  first(1) = 1
  do i = 2, nid+2
    first(i) = first(i) + first(i-1)
  enddo
  do e = 1, ne
    do s = 1, 4
      call side_ids(e, s, lo, up)
      if ( lo .eq. up ) cycle
      j = first(lo+1)
      hi(j) = up ; elm(j) = e ; sid(j) = s
      first(lo+1) = j + 1
    enddo
  enddo

  ! --- A side owned by exactly one element is exterior
  do lo = 1, nid
    do i = first(lo), first(lo+1)-1
      n_dup = 0
      do j = first(lo), first(lo+1)-1
        if ( hi(j) .eq. hi(i) ) n_dup = n_dup + 1
      enddo
      if ( n_dup .eq. 1 ) exterior(sid(i), elm(i)) = .true.
      if ( n_dup .gt. 2 ) then
        write(*,*) 'ERROR: boundary_edges_build found a side shared by more than two elements, EXITING!'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      endif
    enddo
  enddo

  deallocate( first, hi, elm, sid )

contains

  subroutine side_ids(e, s, lo, up)
    integer, intent(in)  :: e, s
    integer, intent(out) :: lo, up
    integer :: ia, ib
    ia = node_list%node(element_list%element(e)%vertex(s))%index(1)
    ib = node_list%node(element_list%element(e)%vertex(mod(s,4)+1))%index(1)
    lo = min(ia, ib) ; up = max(ia, ib)
  end subroutine side_ids

end subroutine boundary_edges_build


logical function boundary_edges_active()
  implicit none
  boundary_edges_active = allocated(exterior)
end function boundary_edges_active


logical function boundary_edge_is_exterior(ielm, side)
  implicit none
  integer, intent(in) :: ielm, side
  boundary_edge_is_exterior = exterior(side, ielm)
end function boundary_edge_is_exterior


end module mod_boundary_edges
