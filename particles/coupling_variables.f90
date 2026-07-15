!> storage of the variables required in each of the coupling schemes
!> used to construct the indices of the feedback from the particle evolution
!> and the indices of the aux_node_list
!> PLEASE CAREFULLY OBSERVE EXISTING PARAMETERS TO AVOID OVERLAP WHEN ADDING A NEW COUPLING SCHEME
module coupling_variables
  implicit none

  integer, parameter :: n_aux_var_max = 100

  ! ====== variables names and number of coupling schemes ====== !

  integer, parameter :: var_name_len  = 15
  character(len=var_name_len), dimension(n_aux_var_max) :: coupling_vars

  ! NCS (Coupling scheme for neutral particles)
#ifdef WITH_TiTe
    character(len=var_name_len), dimension(4) :: ncs_var_names = [character(len=var_name_len) :: &
      "rho",      & !> density
      "mom_par",  & !> parallel momentum
      "E_Te",     & !> electron energy
      "E_Ti"      & !> ion energy
    ]
#else
    character(len=var_name_len), dimension(3) :: ncs_var_names = [character(len=var_name_len) :: &
      "rho",      & !> density
      "mom_par",  & !> parallel momentum
      "E"         & !> total energy
    ]
#endif

  ! ICS (Coupling scheme for impurity particles)
  !> These are only the base variables, there is also the impurity charge density, unique to each impurity group
#ifdef WITH_TiTe
    character(len=var_name_len), dimension(3) :: ics_var_names = [character(len=var_name_len) :: &
      "mom_par",  & !> parallel momentum
      "E_Te",     & !> electron energy
      "E_Ti"      & !> ion energy
    ]
#else
    character(len=var_name_len), dimension(2) :: ics_var_names = [character(len=var_name_len) :: &
      "mom_par",  & !> parallel momentum
      "E"         & !> total energy
    ]
#endif

  ! ICS full force-density coupling channels (Strien 2022), only registered when use_ics_full_force_coupling
  character(len=var_name_len), dimension(6) :: ics_force_var_names = [character(len=var_name_len) :: &
    "fk_par",   & !> B.(total momentum to plasma per unit time): Lorentz reaction + collisions
    "fk_R",     & !> R component of total force density on the plasma
    "fk_Z",     & !> Z component of total force density on the plasma
    "Rk_par",   & !> B.(collisional momentum to plasma per unit time)
    "Rk_R",     & !> R component of collisional force density on the plasma
    "Rk_Z"      & !> Z component of collisional force density on the plasma
  ]

  ! REP (Pressure coupling scheme for runaway electrons)
  character(len=var_name_len), dimension(3) :: rep_var_names = [character(len=var_name_len) :: &
    "P_par",    & !> parallel component of dynamic pressure tensor
    "P_perp",   & !> perpendicular component of dynamic pressure tensor
    "j_Phi"     & !> Phi compoment of current  
  ]

  ! EPF (Full pressure coupling for energetic particles)
  character(len=var_name_len), dimension(7) :: epf_var_names = [character(len=var_name_len) :: &
    "PI_RR",     & !> (R,R)     component of full pressure tensor
    "PI_ZZ",     & !> (Z,Z)     component of full pressure tensor
    "PI_PHIPHI", & !> (Phi,Phi) component of full pressure tensor
    "PI_RZ",     & !> (R,Z)     component of full pressure tensor
    "PI_RPHI",   & !> (Z,Phi)   component of full pressure tensor
    "PI_ZPHI",   & !> (R,Phi)   component of full pressure tensor
    "rho_ep"     & !>           density
  ]


  ! =========== Storage variables for kinetic coupling indices ======== !

  !> variables indices
  integer :: rho_idx_kin      = 0
  integer :: mom_par_idx_kin  = 0
#ifdef WITH_TiTe
    integer :: E_Te_idx_kin   = 0
    integer :: E_Ti_idx_kin   = 0 
#else
    integer :: E_idx_kin      = 0
#endif
  integer :: P_par_idx_kin    = 0
  integer :: P_perp_idx_kin   = 0
  integer :: j_Phi_idx_kin    = 0
  integer :: rho_ep_idx_kin    = 0
  integer :: PI_RR_idx_kin     = 0
  integer :: PI_ZZ_idx_kin     = 0
  integer :: PI_PHIPHI_idx_kin = 0
  integer :: PI_RZ_idx_kin     = 0
  integer :: PI_RPHI_idx_kin   = 0
  integer :: PI_ZPHI_idx_kin   = 0

  !> index of coupling variables specific to each impurity group
  integer :: ics_indices_kin(n_aux_var_max) = -1

  !> full force-density coupling channels (use_ics_full_force_coupling, Strien 2022):
  !> fk_* = TOTAL momentum given to the plasma (Lorentz-force reaction + collisions) per unit time,
  !> Rk_* = COLLISIONAL-only momentum given to the plasma (for the V.R_k work terms in the energy equation)
  integer :: fk_par_idx_kin   = 0  !< B.(dp/dt) of total force density (enters parallel momentum equation)
  integer :: fk_R_idx_kin     = 0  !< R component of total force density (enters vorticity equation)
  integer :: fk_Z_idx_kin     = 0  !< Z component of total force density (enters vorticity equation)
  integer :: Rk_par_idx_kin   = 0  !< B.(dp/dt) of collisional force density (parallel flow work term)
  integer :: Rk_R_idx_kin     = 0  !< R component of collisional force density (ExB flow work term)
  integer :: Rk_Z_idx_kin     = 0  !< Z component of collisional force density (ExB flow work term)

  !> named diagnostic projection slots (previously hard-coded literals 6/7/8 in evolve_ncs_ics,
  !> which collide with registered variables for some scheme combinations)
  integer :: ncs_dens_diag_idx_kin = 0 !< NCS neutral density diagnostic projection
  integer :: ics_prad_diag_idx_kin = 0 !< ICS impurity radiated power diagnostic projection
  integer :: ics_dens_diag_idx_kin = 0 !< ICS impurity density diagnostic projection

  !> hybrid Crank-Nicolson coupling (use_kin_cn_coupling, Strien 2022 Eq. 3.55):
  !> prev_aux_idx(k) holds the aux slot containing the PREVIOUS fluid interval's value of coupling slot k
  !> (0 for slots that are not time-interval integrals, e.g. imp_q and the diagnostic slots)
  integer :: prev_aux_idx(n_aux_var_max) = 0

end module coupling_variables
