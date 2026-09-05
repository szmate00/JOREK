! Minimal dependency fixtures. Tests compile the PRODUCTION sheath/assembly modules.
module constants
  implicit none
  real*8, parameter :: MU_ZERO=1.2566370614359173d-6, ATOMIC_MASS_UNIT=1.66053906660d-27
  real*8, parameter :: EL_CHG=1.602176634d-19, PI=3.141592653589793d0, MASS_ELECTRON=9.1093837015d-31
end module
module phys_module
  implicit none
  real*8 :: F0=2.972306d0, central_density=1.011088d0, central_mass=2.01410174369812d0
  real*8 :: GAMMA=5.d0/3.d0, sheath_Lambda=3.d0, sheath_V_wall=0.d0
  logical :: sheath_Lambda_local=.false., sheath_jsat_from_vpar=.false., vpar_smoothing=.true.
  real*8 :: vpar_smoothing_coef(3)=[0.02d0,0.016d0,0.005754d0]
  real*8 :: sheath_X_min=-3.d0, sheath_smooth_dX=0.5d0
  real*8 :: sheath_jsat_vpar_min=0.d0, sheath_weak_detmin=0.d0
  real*8 :: T_min_neg=3.d-5, T_min_sheath=-1.d0, corr_neg_temp_coef(2)=[0.1d0,0.1d0]
  real*8 :: sheath_V_wall_asym=0.d0, sheath_V_wall_R0=1.42d0, sheath_V_wall_dR=0.01d0
  integer, parameter :: max_bnd_types=32
end module
module mod_parameters
  implicit none
  integer, parameter :: n_var=7, var_zj=2
end module
module mod_integer_types
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  integer, parameter :: int_all=int64
end module
module data_structure
  use mod_integer_types
  implicit none
  type type_SP_MATRIX
    integer :: i_tor_min=1, i_tor_max=1
    integer(kind=int_all), allocatable :: irn(:), jcn(:), ijA_size(:), ijA_index(:,:), irn_jcn(:,:)
    real*8, allocatable :: val(:)
  end type
  type type_node
    integer :: index(4)=0
  end type
  type type_element
    integer :: vertex(4)=0, n_sons=0
  end type
  type type_node_list
    type(type_node) :: node(12)
  end type
  type type_element_list
    integer :: n_elements=0
    type(type_element) :: element(3)
  end type
end module
module mpi_mod
  implicit none
  integer, parameter :: MPI_COMM_WORLD=0, MPI_REAL8=1, MPI_SUM=2, MPI_MAX=3, MPI_MIN=4
contains
  subroutine MPI_Abort(comm,code,ierr)
    integer :: comm,code,ierr
    ierr=code
    error stop 1
  end subroutine
  subroutine MPI_Reduce(send,recv,n,datatype,op,root,comm,ierr)
    real*8 :: send(*),recv(*)
    integer :: n,datatype,op,root,comm,ierr
    recv(1:n)=send(1:n)
    ierr=0
  end subroutine
end module
