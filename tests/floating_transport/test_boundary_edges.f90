!> Exterior-side table (mod_boundary_edges): two conforming elements sharing one side, with both
!! endpoints of that side carrying boundary labels, a coincident node record, and a collapsed side.
program test_boundary_edges
  use data_structure
  use mod_boundary_edges
  implicit none
  type(type_element_list) :: el
  type(type_node_list)    :: nl
  integer :: i
  logical :: ext(4,3)

  ! nodes 1-6 with distinct value DOFs; node 7 is a second record of vertex 3 (same index(1))
  nl%n_nodes = 7; nl%n_dof = 6
  allocate( nl%node(7) )
  do i = 1, 6
    nl%node(i)%index(1) = i
  enddo
  nl%node(7)%index(1) = 3
  nl%node(:)%boundary = 0
  nl%node(2)%boundary = 1; nl%node(3)%boundary = 4; nl%node(7)%boundary = 4

  ! element 1: 1-2-3-4, element 2: 2-5-6-7 (side 2-3 shared through the coincident record),
  ! element 3: a degenerate element whose side 1 is collapsed (5-5)
  el%n_elements = 3
  el%element(1)%vertex = [1,2,3,4]
  el%element(2)%vertex = [2,5,6,7]
  el%element(3)%vertex = [5,5,6,4]

  call boundary_edges_build(el, nl, 0)
  if ( .not. boundary_edges_active() ) error stop 'FAIL: table not built'
  do i = 1, 4
    ext(i,1) = boundary_edge_is_exterior(1,i)
    ext(i,2) = boundary_edge_is_exterior(2,i)
    ext(i,3) = boundary_edge_is_exterior(3,i)
  enddo
  ! element 1: sides 1 (1-2), 3 (3-4) exterior, side 2 (2-3) shared, side 4 (4-1) exterior
  if ( any(ext(:,1) .neqv. [.true., .false., .true., .true.]) ) error stop 'FAIL: element 1 exterior sides'
  ! element 2: side 4 (7-2 = 3-2) shared, side 3 (6-7 = 6-3) exterior ... side 2 (5-6) is shared with element 3
  if ( any(ext(:,2) .neqv. [.true., .false., .true., .false.]) ) error stop 'FAIL: element 2 exterior sides'
  ! element 3: side 1 collapsed (never exterior), side 2 (5-6) shared, sides 3 (6-4), 4 (4-5) exterior
  if ( any(ext(:,3) .neqv. [.false., .false., .true., .true.]) ) error stop 'FAIL: element 3 exterior sides'
  write(*,'(a)') ' PASS: exterior sides from connectivity, shared side with two labelled endpoints is interior'
end program
