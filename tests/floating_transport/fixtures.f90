! Small serial dependency fixtures for compiling the PRODUCTION boundary assembler.
! They do not stand in for an MPI production simulation.
module mod_parameters
  implicit none
  integer, parameter :: n_var=8,n_tor=1,n_plane=1,n_vertex_max=4,n_degrees=4,n_degrees_1d=2
  integer, parameter :: var_psi=1,var_u=2,var_zj=3,var_w=4,var_rho=5,var_Ti=6,var_vpar=7,var_Te=8
  integer, parameter :: var_T=6,var_rhon=5
  logical, parameter :: with_TiTe=.true.,with_vpar=.true.,with_neutrals=.false.
end module
module mod_model_settings
  use mod_parameters, only: with_TiTe
end module
module constants
  implicit none
  real*8, parameter :: PI=acos(-1.d0),MU_ZERO=4.d-7*PI,ATOMIC_MASS_UNIT=1.660538921d-27,EL_CHG=1.602176565d-19
end module
module phys_module
  use mod_parameters
  implicit none
  integer, parameter :: max_bnd_types=32
  type natural_type
    logical :: rho=.true.,Ti=.true.,Te=.true.,T=.true.,vpar=.true.,rhon=.false.
  end type
  type bc_type
    type(natural_type) :: natural
    logical :: floating_u=.false.,mach1=.true.
  end type
  type(bc_type) :: bcs(0:max_bnd_types)
  real*8 :: time_evol_theta=1.d0,time_evol_zeta=0.5d0,tstep=1.d0,tstep_prev=1.d0
  real*8 :: min_sheath_angle=1.d0,F0=2.97d0,gamma=5.d0/3.d0
  real*8 :: T_min_neg=3.d-5,T_1=0.01d0,corr_neg_temp_coef(2)=[0.5d0,0.5d0]
  real*8 :: central_density=1.d0,central_mass=2.014d0,sheath_Lambda=3.d0,sheath_V_wall=0.d0
  logical :: floating_u_mach_flux=.false.,floating_u_wall_flux=.false.,floating_u_transport_diag=.false.
  real*8 :: floating_u_probe_R=1.6d0,floating_u_probe_Z=-1.11d0
  logical :: vpar_smoothing=.false.,mach_one_bnd_integral=.false.
  real*8 :: vpar_smoothing_coef(3)=[0.02d0,0.016d0,0.005754d0]
  real*8 :: density_reflection=0.d0,neutral_reflection=0.d0,visco_par_heating=0.d0
  real*8 :: gamma_sheath_i=0.6d0,gamma_sheath_e=3.d0,gamma_sheath=3.d0
  real*8 :: neutral_line_R_start(10)=0.d0,neutral_line_R_end(10)=0.d0
  real*8 :: neutral_line_Z_start(10)=0.d0,neutral_line_Z_end(10)=0.d0,neutral_line_source(10)=0.d0
end module
module data_structure
  use mod_parameters
  implicit none
  type type_node
    real*8 :: x(1,4,2)=0.d0,values(1,4,n_var)=0.d0,deltas(1,4,n_var)=0.d0
    integer :: boundary=1,index(4)=0
  end type
  type type_element
    integer :: vertex(4)=[1,2,3,4],n_sons=0
    real*8 :: size(4,4)=1.d0
  end type
  type type_node_list
    type(type_node) :: node(12)
  end type
  type type_element_list
    integer :: n_elements=1
    type(type_element) :: element(3)
  end type
end module
module gauss
  implicit none
  integer, parameter :: n_gauss=4
  real*8 :: xgauss(4)=[0.0694318442029737d0,0.330009478207572d0,0.669990521792428d0,0.930568155797026d0]
  real*8 :: wgauss(4)=[0.173927422568727d0,0.326072577431273d0,0.326072577431273d0,0.173927422568727d0]
end module
module basis_at_gaussian
  use gauss
  implicit none
  real*8 :: H1(2,2,4),H1_s(2,2,4),H1_ss(2,2,4),HZ(1,1)=1.d0,HZ_p(1,1)=0.d0
contains
  subroutine set_basis()
    integer :: i
    real*8 :: s
    do i=1,4
      s=xgauss(i)
      H1(:,1,i)=[2*s**3-3*s**2+1,-2*s**3+3*s**2]
      H1(:,2,i)=[s**3-2*s**2+s,s**3-s**2]
      H1_s(:,1,i)=[6*s*s-6*s,-6*s*s+6*s]
      H1_s(:,2,i)=[3*s*s-4*s+1,3*s*s-2*s]
      H1_ss(:,1,i)=[12*s-6,-12*s+6]
      H1_ss(:,2,i)=[6*s-4,6*s-2]
    enddo
  end subroutine
end module
module corr_neg
contains
  real*8 function corr_neg_temp1(t) result(c)
    use phys_module, only: T_min_neg,corr_neg_temp_coef
    real*8,intent(in)::t
    real*8::knee
    knee=T_min_neg*sum(corr_neg_temp_coef)
    c=t
    if (t<knee) c=T_min_neg*corr_neg_temp_coef(1)+T_min_neg*corr_neg_temp_coef(2)*exp((t-knee)/(T_min_neg*corr_neg_temp_coef(2)))
  end function
  real*8 function corr_neg_dens(t) result(c)
    real*8,intent(in)::t
    c=max(t,1.d-10)
  end function
end module
module diffusivities
contains
  real*8 function get_dperp(t)
    real*8::t
    get_dperp=t
  end function
  real*8 function get_zkperp(t)
    real*8::t
    get_zkperp=t
  end function
end module
module mpi_mod
  implicit none
  integer,parameter :: MPI_COMM_WORLD=0,MPI_DOUBLE_PRECISION=1,MPI_2DOUBLE_PRECISION=2,MPI_MINLOC=3
contains
  subroutine MPI_Allreduce(send,recv,n,datatype,op,comm,ierr)
    real*8 :: send(*),recv(*)
    integer :: n,datatype,op,comm,ierr,k
    k=n
    if(datatype==MPI_2DOUBLE_PRECISION) k=2*n
    recv(1:k)=send(1:k); ierr=0
  end subroutine
  subroutine MPI_Bcast(buffer,n,datatype,root,comm,ierr)
    real*8::buffer(*)
    integer::n,datatype,root,comm,ierr
    ierr=0
  end subroutine
  subroutine MPI_Abort(comm,code,ierr)
    integer::comm,code,ierr
    ierr=code
    error stop 'MPI_Abort fixture'
  end subroutine
end module
module mod_interp
contains
  subroutine interp_PRZ_1(nodes,elements,e,vars,n,s,t,phi,p,ps,pt,pp,R,Rs,Rt,Z,Zs,Zt,deltas)
    use data_structure
    type(type_node_list),intent(in)::nodes
    type(type_element_list),intent(in)::elements
    integer,intent(in)::e,n,vars(n)
    real*8,intent(in)::s,t,phi
    real*8,intent(out)::p(n),ps(n),pt(n),pp(n),R,Rs,Rt,Z,Zs,Zt
    logical,optional,intent(in)::deltas
    error stop 'post-update interpolation is not exercised by the boundary fixture'
  end subroutine
end module
