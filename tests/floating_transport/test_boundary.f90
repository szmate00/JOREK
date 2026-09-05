program test_boundary
  use mod_parameters
  use phys_module
  use data_structure
  use mod_boundary_matrix_open
  use mod_floating_boundary_edges
  use mod_floating_transport_diag, only: transport_diag_reset,transport_diag_report
  use basis_at_gaussian, only: set_basis
  implicit none
  integer,parameter :: nd=4*4*n_var
  type(type_element) :: e
  type(type_node) :: nodes(4),base(4)
  type(type_element_list) :: el
  type(type_node_list) :: nl
  real*8 :: a(nd,nd),ap(nd,nd),am(nd,nd),r(nd),rp(nd),rm(nd),err,scale,eps
  integer :: i,dof,var,col,mode,k,row,icase,ib
  integer,parameter :: variables(5)=[var_u,var_vpar,var_Ti,var_Te,var_rho]
  integer,parameter :: boundary_types(5)=[1,3,4,5,9]
  call set_basis()
  do i=1,4
    nl%node(i)%index(1)=i
  enddo
  call floating_edges_build(el,nl)
  ! Affine side 1: R runs from 1 to 2, Z=0, exterior normal points down.
  base(1)%x(1,1,:)=[1.d0,0.d0]; base(2)%x(1,1,:)=[2.d0,0.d0]
  base(3)%x(1,1,:)=[2.d0,1.d0]; base(4)%x(1,1,:)=[1.d0,1.d0]
  do i=1,4
    base(i)%x(1,2,:)=[1.d0,0.d0]
    base(i)%x(1,3,:)=[0.d0,1.d0/3.d0]
    base(i)%values(1,1,var_rho)=0.2d0
    base(i)%values(1,1,var_Ti)=0.003d0
    base(i)%values(1,1,var_Te)=0.004d0
    base(i)%values(1,1,var_vpar)=0.04d0
    base(i)%values(1,1,var_psi)=0.08d0*base(i)%x(1,1,1)
    base(i)%values(1,2,var_psi)=0.08d0
    base(i)%values(1,1,var_u)=0.01d0*base(i)%x(1,1,1)
    base(i)%values(1,2,var_u)=0.01d0
  enddo
  bcs(1)%floating_u=.true.
  eps=1.d-7
  visco_par_heating=0.02d0 ! exercises the independent normal-derivative heat columns
  do icase=1,3
    if (icase==2) then
      base%boundary=4; bcs(4)%floating_u=.true.
      do i=1,4
        base(i)%values(:,:,var_u)=-base(i)%values(:,:,var_u)
      enddo
    elseif (icase==3) then
      base%boundary=9; bcs(9)%floating_u=.true.
      do i=1,4
        base(i)%values(1,1,var_Ti)=1.d-5
        base(i)%values(1,1,var_Te)=1.d-5
      enddo
      eps=1.d-9
    endif
  do mode=1,3
    floating_u_mach_flux=mode/=2
    floating_u_wall_flux=mode/=1
    nodes=base; call assemble(a,r)
    do k=1,5
      var=variables(k)
      do i=1,2
        do dof=1,4
          col=n_var*4*(i-1)+n_var*(dof-1)+var
          nodes=base; nodes(i)%values(1,dof,var)=nodes(i)%values(1,dof,var)+eps
          call assemble(ap,rp)
          nodes=base; nodes(i)%values(1,dof,var)=nodes(i)%values(1,dof,var)-eps
          call assemble(am,rm)
          ! ELM is minus d(RHS)/dx; theta=1. Only overridden rows are tested.
          do row=1,nd
            if (mod(row-1,n_var)+1==var_vpar .and. .not.floating_u_mach_flux) cycle
            if (mod(row-1,n_var)+1/=var_vpar .and. .not.floating_u_wall_flux) cycle
            scale=max(1.d0,maxval(abs(a(:,col))))
            err=abs(a(row,col)+(rp(row)-rm(row))/(2*eps))/scale
            if(err>2.d-7) then
              write(*,*) 'FAIL boundary FD mode,row,col,err',mode,row,col,err
              error stop 1
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
  enddo
  ! Relabelling a fully covered physical edge must not change its equations.
  floating_u_mach_flux=.true.; floating_u_wall_flux=.true.
  nodes=base; call assemble(a,r)
  do ib=1,size(boundary_types)
    bcs(boundary_types(ib))%floating_u=.true.
    nodes=base; nodes%boundary=boundary_types(ib)
    call assemble(ap,rp)
    if(any(a/=ap).or.any(r/=rp)) error stop 'boundary type changed an identical physical closure'
  enddo
  ! Coverage: a type-3 corner must not activate a type-2 artificial edge.
  nodes=base; nodes(1)%boundary=3; nodes(2)%boundary=2
  bcs(3)%floating_u=.true.; bcs(2)%floating_u=.false.
  floating_u_mach_flux=.false.; floating_u_wall_flux=.false.; call assemble(a,r)
  floating_u_mach_flux=.true.; floating_u_wall_flux=.true.; call assemble(ap,rp)
  if(any(a/=ap).or.any(r/=rp)) error stop 'type 2 edge changed'
  ! Diagnostic-enabled and disabled matrices/residuals are identical.
  nodes=base
  floating_u_mach_flux=.false.; floating_u_wall_flux=.false.
  floating_u_transport_diag=.false.; call assemble(a,r)
  call transport_diag_reset()
  floating_u_transport_diag=.true.; call assemble(ap,rp)
  if(any(a/=ap).or.any(r/=rp)) error stop 'diagnostic changed equations'
  call transport_diag_report(0,'serial fixture')
  write(*,*) 'PASS: production boundary assembler finite-difference Jacobians and type-2 exclusion'
contains
  subroutine assemble(mat,rhs)
    real*8,intent(out)::mat(nd,nd),rhs(nd)
    integer::vertices(2),directions(2)
    vertices=[1,2]; directions=[1,2]
    mat=0.d0; rhs=0.d0
    call boundary_matrix_open(vertices,directions,e,nodes,.true.,1,1.d0,0.d0,0.d0,1.d0, &
        [1.d0,1.d0],[0.d0,0.d0],mat,rhs,1,1,1)
  end subroutine
end program
