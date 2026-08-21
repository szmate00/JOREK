! Stub modules: just enough to let gfortran syntax-check models/model600/mod_boundary_matrix_open.f90
! outside a full JOREK build. Scaffolding only, never compiled into JOREK.
module mod_parameters
  implicit none
  integer, parameter :: n_var=8, n_plane=4, n_tor=1, n_degrees=4, n_vertex_max=4, n_degrees_1d=2, n_order=3
  integer, parameter :: var_psi=1, var_u=2, var_w=3, var_zj=4, var_rho=5, var_T=6, var_Ti=0, var_Te=0
  integer, parameter :: var_vpar=7, var_rhon=8, var_rhoimp=0, var_nre=0
  logical, parameter :: with_TiTe=.false., with_vpar=.true., with_neutrals=.true.
end module mod_parameters

module data_structure
  use mod_parameters
  implicit none
  type type_node
    integer :: boundary
    real*8  :: x(n_tor,n_degrees,2), values(n_tor,n_degrees,n_var), deltas(n_tor,n_degrees,n_var)
  end type
  type type_element
    real*8 :: size(n_vertex_max,n_degrees)
  end type
  type type_SP_MATRIX
    integer :: i_tor_min, i_tor_max
  end type
  type type_node_list
    integer :: n_nodes
    type(type_node), allocatable :: node(:)
  end type
end module data_structure

module gauss
  use mod_parameters
  implicit none
  integer, parameter :: n_gauss=4
  real*8 :: wgauss(n_gauss), xgauss(n_gauss)
end module gauss

module basis_at_gaussian
  use mod_parameters
  use gauss
  implicit none
  real*8 :: H1(n_vertex_max,n_degrees,n_gauss), H1_s(n_vertex_max,n_degrees,n_gauss)
  real*8 :: H1_ss(n_vertex_max,n_degrees,n_gauss)
  real*8 :: HZ(n_tor,n_plane), HZ_p(n_tor,n_plane)
end module basis_at_gaussian

module phys_module
  use mod_parameters
  implicit none
  integer, parameter :: max_bnd_types = 30
  type type_dirichlet_bc
    logical :: u, w, zj, psi
  end type
  type type_natural_bc
    logical :: u, w, zj, rho, T, Ti, Te, Vpar, rhon
  end type
  type type_bcs
    type(type_dirichlet_bc) :: dirichlet
    type(type_natural_bc)   :: natural
    logical :: mach1, sheath_u, sheath_zj
  end type
  type(type_bcs) :: bcs(max_bnd_types)
  real*8  :: F0, GAMMA, central_density, central_mass, time_evol_theta, time_evol_zeta
  real*8  :: tstep, tstep_prev, t_now, t_start, min_sheath_angle
  real*8  :: sheath_V_wall, sheath_Lambda, sheath_X_min, sheath_smooth_dX, sheath_min_bn
  real*8  :: sheath_sat_slope, sheath_wall_pen, sheath_zj_ratio_max, sheath_zj_relax
  real*8  :: sheath_ramp_time, vpar_smoothing_coef(3)
  real*8  :: sheath_stiff_max, sheath_flux_sign
  logical :: sheath_init_u
  real*8  :: gamma_sheath, gamma_sheath_i, gamma_sheath_e, density_reflection, neutral_reflection
  real*8  :: visco_par_heating
  real*8  :: neutral_line_R_start(10), neutral_line_R_end(10), neutral_line_Z_start(10)
  real*8  :: neutral_line_Z_end(10), neutral_line_source(10)
  logical :: sheath_Lambda_local, vpar_smoothing, mach_one_bnd_integral
end module phys_module

module corr_neg
  implicit none
contains
  real*8 function corr_neg_temp1(v);  real*8 :: v; corr_neg_temp1 = v; end function
  real*8 function corr_neg_dens(v);   real*8 :: v; corr_neg_dens  = v; end function
  real*8 function dcorr_neg_temp_dT1(v);   real*8 :: v; dcorr_neg_temp_dT1 = 1.d0; end function
  real*8 function dcorr_neg_dens_drho1(v); real*8 :: v; dcorr_neg_dens_drho1 = 1.d0; end function
end module corr_neg

module mod_interp
  implicit none
end module mod_interp

module diffusivities
  implicit none
contains
  subroutine get_dperp();  end subroutine
  subroutine get_zkperp(); end subroutine
end module diffusivities
