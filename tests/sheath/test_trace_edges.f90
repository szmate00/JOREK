program test_trace_edges
  use mod_parameters
  use data_structure
  use mod_assembly
  use mod_sheath_trace
  use mod_sheath_boundary_edges
  implicit none
  type(type_SP_MATRIX) :: matrix
  type(type_element_list) :: elements
  type(type_node_list) :: nodes
  real*8, allocatable :: before(:), first_result(:)
  real*8 :: rhs(3*n_var), vals(2), first_rhs
  integer :: i,nt,width,start,j,k,cols(2),vars(2),phase
  logical :: target
  ! Dense node blocks in the sparse structure; only global DOF 2 is owned.
  nt=2; width=n_var*nt
  allocate(matrix%val(3*width**2),matrix%irn(3*width**2),matrix%jcn(3*width**2))
  allocate(matrix%ijA_size(1),matrix%ijA_index(1,3),matrix%irn_jcn(1,3))
  matrix%i_tor_max=nt
  matrix%ijA_size=3
  do i=1,3
    matrix%ijA_index(1,i)=1+(i-1)*width**2
    matrix%irn_jcn(1,i)=i
  enddo
  matrix%val=17.d0
  matrix%irn=42
  matrix%jcn=43
  before=matrix%val
  call boundary_conditions_clear_row(2,var_zj,2,2,2,matrix)
  do i=1,size(matrix%val)
    j=mod(i-1,width**2)/width
    target=j==(var_zj-1)*nt+1
    if (target) then
      call require(matrix%val(i)==0.d0,'all columns/harmonics cleared')
    else
      call require(matrix%val(i)==before(i),'other scalar equations untouched')
    endif
  enddo
  call require(all(matrix%irn==42).and.all(matrix%jcn==43),'clear preserves sparse metadata')
  before=matrix%val
  call boundary_conditions_clear_row(1,var_zj,2,2,2,matrix)
  call require(all(matrix%val==before),'unowned clear is a no-op')
  deallocate(matrix%val,matrix%irn,matrix%jcn,before)

  nt=1; width=n_var
  matrix%i_tor_max=nt
  allocate(matrix%val(3*width**2),matrix%irn(3*width**2),matrix%jcn(3*width**2))
  do i=1,3
    matrix%ijA_index(1,i)=1+(i-1)*width**2
  enddo
  do phase=1,2
    ! A repeated construction must not retain any previous trace contributions.
    call sheath_trace_reset(4)
    matrix%val=17.d0
    rhs=19.d0
    cols=[2,3]; vars=[var_zj,1]; vals=[2.d-14,-4.d-10]
    call sheath_trace_add(2,5,1.d-14,1.d-14,1.d-14,3.d-9,0.d0,1.d0,0.04d0,2,cols,vars,vals)
    call sheath_trace_add(2,5,1.d-14,1.d-14,1.d-14,3.d-9,0.d0,1.d0,0.04d0,2,cols,vars,vals)
    ! A replicated row owned elsewhere must not be applied or reported as owned.
    call sheath_trace_add(1,1,1.d0,1.d0,1.d0,1.d0,0.d0,1.d0,1.d0,2,cols,vars,vals)
    ! A different equation on the SAME geometric DOF must have its own slot.
    cols=[2,2]; vars=[3,4]; vals=[2.d0,6.d0]
    call sheath_trace_add(2,5,2.d0,2.d0,2.d0,10.d0,0.d0,1.d0,0.04d0,2,cols,vars,vals,equation=3)
    call sheath_trace_apply(1,10.d0**(phase*10),2,2,matrix,rhs)
    start=matrix%ijA_index(1,2)+2*width
    call require(abs(matrix%val(start+2)-1.d0/3.d0)<1.d-15,'same-DOF flow equation isolation')
    call require(matrix%val(start+3)==1.d0,'flow thermal column')
    call require(abs(rhs(n_var+3)-5.d0/3.d0)<1.d-14,'flow RHS isolation')
    do i=1,3
      start=matrix%ijA_index(1,i)+(var_zj-1)*width
      do k=1,n_var
        if (i==2.and.k==var_zj) then
          call require(abs(matrix%val(start+k-1)-5.d-5)<1.d-18,'merged/scaled current coefficient')
        elseif (i==3.and.k==1) then
          call require(matrix%val(start+k-1)==-1.d0,'max coefficient normalization')
        else
          call require(matrix%val(start+k-1)==0.d0,'no volume contamination survives')
        endif
      enddo
    enddo
    call require(abs(rhs(n_var+var_zj)-7.5d0)<1.d-13,'raw merged RHS')
    if (phase==1) then
      first_result=matrix%val
      first_rhs=rhs(n_var+var_zj)
    else
      call require(all(first_result==matrix%val),'reset and penalty-scale independence')
      call require(first_rhs==rhs(n_var+var_zj),'RHS penalty-scale independence')
    endif
    call sheath_trace_report(0)
  enddo

  ! Two quads sharing a labelled edge: labels cannot make it a wall.
  do i=1,12
    nodes%node(i)%index(1)=i
  enddo
  elements%n_elements=2
  elements%element(1)%vertex=[1,2,3,4]
  elements%element(2)%vertex=[2,5,6,3]
  call sheath_edges_build(elements,nodes)
  call require(.not.sheath_edge_is_exterior(1,2),'shared side excluded')
  call require(.not.sheath_edge_is_exterior(2,4),'reversed shared side excluded')
  call require(sheath_edge_is_exterior(1,1),'ordinary exterior')
  call require(sheath_edge_is_exterior(2,3),'rotated exterior')
  ! Corner node aliases share value DOFs: identity is topological, not a label.
  nodes%node(7)%index(1)=2
  elements%element(2)%vertex=[7,5,6,3]
  call sheath_edges_build(elements,nodes)
  call require(.not.sheath_edge_is_exterior(1,2),'corner alias')
  ! Rebuilding after a topology change must discard old incidence counts.
  elements%n_elements=1
  call sheath_edges_build(elements,nodes)
  call require(sheath_edge_is_exterior(1,2),'topology rebuild')
  print *, 'PASS: exact row replacement, owner filtering, reset and exterior connectivity'
contains
  subroutine require(condition,name)
    logical,intent(in)::condition
    character(*),intent(in)::name
    if (.not.condition) then
      print *, 'FAIL: ',name
      stop 1
    endif
  end subroutine
end program
