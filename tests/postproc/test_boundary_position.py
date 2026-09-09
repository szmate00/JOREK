"""Compile the production boundary sampler against small analytic geometry fixtures.

Run: python3 tests/postproc/test_boundary_position.py
This exercises bnd_pos itself, without requiring the MPI postprocessor build.
"""
from pathlib import Path
import os
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
source = (root / "diagnostics/new_diag/mod_position.f90").read_text()
start = source.index("  function bnd_pos(")
end = source.index("  end function bnd_pos", start) + len("  end function bnd_pos")
production = source[start:end]
fixture = """
module fixture
implicit none
type node
  integer :: boundary=0
end type
type type_node_list
  type(node) :: node(8)
end type
type element
  integer :: vertex(4)
end type
type type_element_list
  type(element) :: element(2)
end type
type type_bnd_node_list
  integer :: unused
end type
type bnd_element
  integer :: side,element,vertex(2)
end type
type type_bnd_element_list
  integer :: n_bnd_elements
  type(bnd_element) :: bnd_element(8)
end type
type t_pol_pos
  integer :: ielm,bnd_contour,bnd_type,bnd_type1,bnd_type2
  real*8 :: s,t,R,Z,R_s,R_t,Z_s,Z_t,length,dl,bnd_normal(2)
end type
type t_pol_pos_list
  type(t_pol_pos),allocatable :: pos(:,:)
end type
real*8 :: curvature=0d0
contains
subroutine alloc_pol_pos(list,shape)
type(t_pol_pos_list) :: list
integer :: shape(2)
allocate(list%pos(shape(1),shape(2)))
end subroutine
subroutine interp_RZ(nodes,elements,ielm,s,t,R,Rs,Rt,Z,Zs,Zt)
type(type_node_list) :: nodes
type(type_element_list) :: elements
integer :: ielm
real*8 :: s,t,R,Rs,Rt,Z,Zs,Zt
R=10d0*(ielm-1)+s; Z=t+curvature*s*s
Rs=1d0; Rt=0d0; Zs=2d0*curvature*s; Zt=1d0
end subroutine
subroutine fill_pol_pos(pos,nodes,elements)
type(t_pol_pos) :: pos
type(type_node_list) :: nodes
type(type_element_list) :: elements
call interp_RZ(nodes,elements,pos%ielm,pos%s,pos%t,pos%R,pos%R_s,pos%R_t,pos%Z,pos%Z_s,pos%Z_t)
end subroutine
"""
checks = """
end module
program test
use fixture
implicit none
type(type_node_list) :: nodes
type(type_element_list) :: elements
type(type_bnd_node_list) :: bn
type(type_bnd_element_list) :: edges
type(t_pol_pos_list) :: points
integer :: i,j,n,offset
real*8 :: exact
elements%element(1)%vertex=(/1,2,3,4/)
elements%element(2)%vertex=(/5,6,7,8/)
nodes%node%boundary=4
nodes%node(2)%boundary=9
edges%n_bnd_elements=8
do j=1,2
  offset=4*(j-1)
  do i=1,4
    edges%bnd_element(offset+i)%element=j
    edges%bnd_element(offset+i)%side=i
    edges%bnd_element(offset+i)%vertex=(/offset+i,offset+mod(i,4)+1/)
  enddo
enddo
do n=1,4,3
  points=bnd_pos(nodes,elements,bn,edges,n)
  do j=1,2
    do i=1,4*n
      offset=(j-1)*4*n+i
      call close(points%pos(1,offset)%length,dble(i-1)/n,'square arc length')
      if(points%pos(1,offset)%bnd_contour/=j) error stop 'disconnected contour'
    enddo
  enddo
  ! Sides 3 and 4 must traverse backwards in local s/t.
  call close(points%pos(1,2*n+1)%R,1d0,'side 3 orientation R')
  call close(points%pos(1,2*n+1)%Z,1d0,'side 3 orientation Z')
  call close(points%pos(1,3*n+1)%Z,1d0,'side 4 orientation')
  call close(points%pos(1,1)%bnd_normal(2),-1d0,'bottom normal')
  call close(points%pos(1,n+1)%bnd_normal(1),1d0,'right normal')
  call close(points%pos(1,2*n+1)%bnd_normal(2),1d0,'top normal')
  call close(points%pos(1,3*n+1)%bnd_normal(1),-1d0,'left normal')
  if(points%pos(1,1)%bnd_type/=4.or.points%pos(1,n+1)%bnd_type/=9) error stop 'node types'
  if(n>1) then
    if(points%pos(1,2)%bnd_type/=-1) error stop 'mixed edge type'
    if(points%pos(1,2)%bnd_type1/=4.or.points%pos(1,2)%bnd_type2/=9) error stop 'endpoint types'
    if(points%pos(1,2*n+2)%bnd_type/=4) error stop 'uniform edge type'
  endif
enddo
! Reversed traversal of a curved edge must use arc length, not endpoint chords.
edges%n_bnd_elements=2
edges%bnd_element(1)%vertex=(/2,1/)
edges%bnd_element(2)%side=4
edges%bnd_element(2)%vertex=(/1,4/)
curvature=0.5d0
points=bnd_pos(nodes,elements,bn,edges,4)
exact=0.5d0*(sqrt(2d0)+asinh(1d0))
call close(points%pos(1,5)%length,exact,'curved complete edge')
call close(points%pos(1,1)%R,1d0,'reversed edge starts at vertex 2')
if(points%pos(1,1)%bnd_type/=9) error stop 'reversed node type'
call close(points%pos(1,8)%length,exact+0.75d0,'second edge accumulation')
print *, 'PASS: production bnd_pos arc length, orientation, contours, normals and types'
contains
subroutine close(a,b,label)
real*8 :: a,b
character(*) :: label
if(abs(a-b)>1d-11) then
  print *,label,a,b
  error stop 'boundary sampler mismatch'
endif
end subroutine
end program
"""
with tempfile.TemporaryDirectory(prefix="jorek-bnd-position-") as directory:
    path = Path(directory)
    (path / "test.f90").write_text(fixture + production + checks)
    subprocess.run(
        [os.environ.get("FC", "gfortran"), "-ffree-line-length-none", "-fcheck=all",
         "-ffpe-trap=invalid,zero,overflow", "-finit-real=snan", "test.f90", "-o", "test"],
        cwd=path, check=True,
    )
    subprocess.run([str(path / "test")], cwd=path, check=True)
