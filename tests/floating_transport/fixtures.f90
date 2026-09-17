! Serial dependency fixtures for compiling the PRODUCTION model600 boundary assembler
! (mod_boundary_matrix_open.f90) without MPI. They stand in for the modules it imports,
! not for a simulation. WITH_TiTe build, one toroidal plane, bicubic elements.
module mod_parameters
  implicit none
  integer, parameter :: n_var=8, n_tor=1, n_plane=1, n_vertex_max=4, n_degrees=4, n_degrees_1d=2
  integer, parameter :: var_psi=1, var_u=2, var_zj=3, var_w=4, var_rho=5, var_Ti=6, var_vpar=7, var_Te=8
  integer, parameter :: var_T=6, var_rhon=5
  logical, parameter :: with_TiTe=.true., with_vpar=.true., with_neutrals=.false.
end module
module mod_model_settings
  use mod_parameters, only: with_TiTe, with_vpar
end module
module constants
  implicit none
  real*8, parameter :: PI=acos(-1.d0), MU_ZERO=4.d-7*PI, ATOMIC_MASS_UNIT=1.660538921d-27, EL_CHG=1.602176565d-19
end module
module phys_module
  use mod_parameters
  implicit none
  integer, parameter :: max_bnd_types=32
  type natural_type
    logical :: rho=.true., Ti=.true., Te=.true., T=.true., vpar=.false., rhon=.false.
  end type
  type bc_type
    type(natural_type) :: natural
    logical :: floating_u=.false., mach1=.true., sheath_j=.false.
  end type
  type(bc_type) :: bcs(0:max_bnd_types)
  real*8  :: time_evol_theta=1.d0, time_evol_zeta=0.5d0, tstep=1.d0, tstep_prev=1.d0
  real*8  :: min_sheath_angle=1.d0, F0=2.97d0, gamma=5.d0/3.d0
  real*8  :: T_min_neg=3.d-5, T_1=0.01d0, corr_neg_temp_coef(2)=[0.5d0,0.5d0]
  real*8  :: central_density=1.d0, central_mass=2.014d0, sheath_Lambda=3.d0, sheath_V_wall=0.d0
  logical :: vpar_smoothing=.false., mach_one_bnd_integral=.false., mach1_weak=.false., mach1_weak_drift=.false., floating_u_diag=.false.
  real*8  :: mach1_weak_drift_bound=0.d0
  logical :: mach1_weak_drift_cut=.false., mach1_weak_cut=.false., mach1_weak_inflow=.true., sheath_j_float_u=.false.
  real*8  :: sheath_j_ramp_time=-1.d0, t_now=0.d0, tstep_n(7)=0.d0
  logical :: sheath_j_current_row=.false.
  real*8  :: sheath_j_ion_slope=0.d0, sheath_j_e_slope=0.d0
  integer :: nstep_n(7)=0
  real*8  :: vpar_smoothing_coef(3)=[0.02d0,0.016d0,0.005754d0]
  real*8  :: density_reflection=0.d0, neutral_reflection=0.d0, visco_par_heating=0.d0
  real*8  :: gamma_sheath_i=0.6d0, gamma_sheath_e=3.d0, gamma_sheath=3.d0
  real*8  :: neutral_line_R_start(10)=0.d0, neutral_line_R_end(10)=0.d0
  real*8  :: neutral_line_Z_start(10)=0.d0, neutral_line_Z_end(10)=0.d0, neutral_line_source(10)=0.d0
end module
module data_structure
  use mod_parameters
  implicit none
  type type_node
    real*8  :: x(1,4,2)=0.d0, values(1,4,n_var)=0.d0, deltas(1,4,n_var)=0.d0
    integer :: boundary=1, index(4)=0
  end type
  type type_element
    integer :: vertex(4)=[1,2,3,4], n_sons=0
    real*8  :: size(4,4)=1.d0
  end type
  type type_node_list
    integer :: n_nodes=0, n_dof=0
    type(type_node), allocatable :: node(:)
  end type
  type type_element_list
    integer :: n_elements=0
    type(type_element) :: element(8)
  end type
contains
  subroutine make_deep_copy_node(a, b)
    type(type_node), intent(in)  :: a
    type(type_node), intent(out) :: b
    b = a
  end subroutine
end module
module mpi_mod
  implicit none
  integer, parameter :: MPI_COMM_WORLD=0, MPI_DOUBLE_PRECISION=1, MPI_MIN=2, MPI_MAX=3, MPI_SUM=4
contains
  subroutine MPI_COMM_SIZE(comm, n, ierr)
    integer :: comm, n, ierr
    n = 1; ierr = comm
  end subroutine
  subroutine MPI_ALLREDUCE(send, recv, n, datatype, op, comm, ierr)
    real*8  :: send(*), recv(*)
    integer :: n, datatype, op, comm, ierr
    recv(1:n) = send(1:n); ierr = datatype + op + comm
  end subroutine
  subroutine MPI_ABORT(comm, code, ierr)
    integer :: comm, code, ierr
    ierr = code + comm
    error stop 'MPI_ABORT (fixture)'
  end subroutine
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
  real*8 :: H1(2,2,4), H1_s(2,2,4), H1_ss(2,2,4), HZ(1,1)=1.d0, HZ_p(1,1)=0.d0
  real*8 :: H(4,4,4,4)   ! 2D bicubic Hermite basis (vertex, dof, ms, mt) from the 1D products
contains
  subroutine set_basis()
    integer :: i
    real*8  :: s
    do i=1,4
      s=xgauss(i)
      H1(:,1,i)   =[ 2*s**3-3*s**2+1, -2*s**3+3*s**2 ]
      H1(:,2,i)   =[ s**3-2*s**2+s,    s**3-s**2      ]
      H1_s(:,1,i) =[ 6*s*s-6*s,       -6*s*s+6*s     ]
      H1_s(:,2,i) =[ 3*s*s-4*s+1,      3*s*s-2*s     ]
      H1_ss(:,1,i)=[ 12*s-6,          -12*s+6        ]
      H1_ss(:,2,i)=[ 6*s-4,            6*s-2         ]
    enddo
    call set_basis_2d()
  end subroutine
  !> vertices 1:(0,0) 2:(1,0) 3:(1,1) 4:(0,1); dofs 1:value 2:d/ds 3:d/dt 4:d2/dsdt
  subroutine set_basis_2d()
    integer :: iv, jd, ms, mt, ia, ib, ja, jb
    do iv = 1, 4
      ia = merge(1, 2, iv==1 .or. iv==4) ; ib = merge(1, 2, iv==1 .or. iv==2)
      do jd = 1, 4
        ja = merge(1, 2, jd==1 .or. jd==3) ; jb = merge(1, 2, jd==1 .or. jd==2)
        do ms = 1, 4
          do mt = 1, 4
            H(iv,jd,ms,mt) = H1(ia,ja,ms) * H1(ib,jb,mt)
          enddo
        enddo
      enddo
    enddo
  end subroutine
end module
module corr_neg
contains
  real*8 function corr_neg_temp1(t) result(c)
    use phys_module, only: T_min_neg, corr_neg_temp_coef
    real*8, intent(in) :: t
    real*8 :: knee
    knee=T_min_neg*sum(corr_neg_temp_coef)
    c=t
    if (t<knee) c=T_min_neg*corr_neg_temp_coef(1)+T_min_neg*corr_neg_temp_coef(2)*exp((t-knee)/(T_min_neg*corr_neg_temp_coef(2)))
  end function
  real*8 function corr_neg_dens(t) result(c)
    real*8, intent(in) :: t
    c=max(t,1.d-10)
  end function
end module
module diffusivities
contains
  real*8 function get_dperp(t)
    real*8 :: t
    get_dperp=t
  end function
  real*8 function get_zkperp(t)
    real*8 :: t
    get_zkperp=t
  end function
end module
module mod_interp
end module
