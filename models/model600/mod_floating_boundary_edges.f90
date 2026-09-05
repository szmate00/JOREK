!> Exterior sides from the existing conforming mesh connectivity, not node labels.
!! Rebuilt before each floating-transport matrix construction, including restart and PC
!! paths. Coordinates/frames are untouched. Shared corner value indices identify
!! coincident topological vertices even when a grid uses separate node records.
module mod_floating_boundary_edges
  implicit none
  private
  public :: floating_edges_build, floating_edge_is_exterior
  logical, allocatable, save :: exterior(:,:)
contains
  subroutine floating_edges_build(element_list, node_list)
    use data_structure, only: type_element_list, type_node_list
    use, intrinsic :: iso_fortran_env, only: int64
    use mpi_mod
    type(type_element_list), intent(in) :: element_list
    type(type_node_list), intent(in) :: node_list
    integer, allocatable :: key1(:), key2(:), counts(:), first_element(:), first_side(:)
    integer :: ne, capacity, e, s, a, b, slot, i, ierr
    ne = element_list%n_elements
    capacity = max(16, 8*ne+1)
    if (allocated(exterior)) deallocate(exterior)
    allocate(exterior(4,ne))
    exterior = .false.
    allocate(key1(capacity),key2(capacity),counts(capacity),first_element(capacity),first_side(capacity))
    key1 = 0
    key2 = 0
    counts = 0
    do e=1,ne
      if (element_list%element(e)%n_sons > 0) then
        write(*,*) 'ERROR: floating transport exterior-edge lookup currently requires an unrefined conforming mesh.'
        call MPI_Abort(MPI_COMM_WORLD,1,ierr)
        return
      endif
      do s=1,4
        a = node_list%node(element_list%element(e)%vertex(s))%index(1)
        b = node_list%node(element_list%element(e)%vertex(mod(s,4)+1))%index(1)
        if (a == b) cycle ! collapsed axis side, not a material wall
        if (a <= 0 .or. b <= 0) then
          write(*,*) 'ERROR: floating transport requires assigned global vertex DOFs.'
          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
          return
        endif
        if (a > b) then
          i=a; a=b; b=i
        endif
        slot = 1+int(modulo(int(a,int64)*104729_int64+int(b,int64),int(capacity,int64)))
        do
          if (counts(slot) == 0) then
            key1(slot)=a; key2(slot)=b
            first_element(slot)=e; first_side(slot)=s
            exit
          endif
          if (key1(slot) == a .and. key2(slot) == b) exit
          slot = mod(slot,capacity)+1
        enddo
        counts(slot)=counts(slot)+1
        if (counts(slot) > 2) then
          write(*,*) 'ERROR: non-manifold edge in floating transport connectivity: ',a,b
          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
          return
        endif
      enddo
    enddo
    do slot=1,capacity
      if (counts(slot) == 1) exterior(first_side(slot),first_element(slot))=.true.
    enddo
  end subroutine floating_edges_build

  logical function floating_edge_is_exterior(element_id, side)
    integer, intent(in) :: element_id, side
    floating_edge_is_exterior = exterior(side,element_id)
  end function floating_edge_is_exterior
end module mod_floating_boundary_edges
