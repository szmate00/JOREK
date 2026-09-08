module mod_boundary_conditions
contains

!*******************************************************************************
!* Subroutine: boundary_condition                                              *
!*******************************************************************************
!*                                                                             *
!* Add boundary condition on the matrix.                                       *
!*                                                                             *
!* Parameters:                                                                 *
!*   my_id        - Identifier of the node in MPI_COMM_WORLD                   *
!*   node_list    - List of nodes                                              *
!*   element_list - List of all elements                                       *
!*   local_elms   - List of local elements                                     *
!*   n_local_elms - Number of local elements                                   *
!*   index_min    - Minimal index of local elements                            *
!*   index_max    - Maximal index of local elements (                          *
!*   xpoint2      -                                                            *
!*   xcase2       -                                                            *
!*   psi_axis     -                                                            *
!*   psi_bnd      -                                                            *
!*   Z_xpoint     -                                                            *
!*                                                                             *
!*******************************************************************************

subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,& 
                                n_local_elms, index_min, index_max, rhs_loc, xpoint2,     &
                                xcase2, R_axis, Z_axis, psi_axis, psi_bnd,                &
                                R_xpoint, Z_xpoint, psi_xpoint, a_mat)

use constants, only : PI, MU_ZERO, ATOMIC_MASS_UNIT, EL_CHG
use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
use data_structure
use vacuum, ONLY: is_freebound
use phys_module, only: F0, GAMMA, freeboundary, RMP_on, psi_RMP_cos, dpsi_RMP_cos_dR, dpsi_RMP_cos_dZ, &
       psi_RMP_sin, dpsi_RMP_sin_dR, dpsi_RMP_sin_dZ, t_now, RMP_growth_rate, RMP_ramp_up_time,            &
       RMP_start_time, tstep, RMP_har_cos, RMP_har_sin, T_min,                                             &
       mach_one_bnd_integral, Vpar_smoothing, vpar_smoothing_coef, no_mach1_bc,                            &
       Number_RMP_harmonics, RMP_har_cos_spectrum,RMP_har_sin_spectrum, grid_to_wall, n_wall_blocks, keep_n0_const, &
       bcs, loop_voltage, central_density, central_mass,                                                   &
       sheath_V_wall, floating_u_diag, D_perp, mach1_omit_drift, floating_u_mach_flux, min_sheath_angle
use mod_floating_u, only: floating_u_norm, mach1_uout_supplement
use tr_module
use mpi_mod
use mod_basisfunctions
use mod_interp
use mod_integer_types
use mod_node_indices

implicit none

! --- Routine parameters
integer,                            intent(in)    :: my_id
type (type_node_list),              intent(in)    :: node_list
type (type_element_list),           intent(in)    :: element_list
type (type_bnd_node_list),          intent(in)    :: bnd_node_list
integer,                            intent(in)    :: local_elms(*)
integer,                            intent(in)    :: n_local_elms
integer,                            intent(in)    :: index_min
integer,                            intent(in)    :: index_max
logical,                            intent(in)    :: xpoint2
integer,                            intent(in)    :: xcase2
real*8,                             intent(in)    :: R_axis
real*8,                             intent(in)    :: Z_axis
real*8,                             intent(in)    :: psi_axis
real*8,                             intent(in)    :: psi_bnd
real*8,                             intent(in)    :: R_xpoint(2)
real*8,                             intent(in)    :: Z_xpoint(2)
real*8,                             intent(in)    :: psi_xpoint(2)
real*8,                             intent(inout) :: rhs_loc(*)
type(type_SP_MATRIX)                              :: a_mat

! Internal parameters
real*8  :: zbig, zbig_backup,  T0, Ti0, Te0, Vpar0, bigR
real*8  :: R_s, R_t, Z, Z_s, Z_t, R_tt, Z_tt, ps0, ps0_s, ps0_t, ps0_tt, ps0_x, ps0_y, direction, xjac
real*8  :: ps0_b, T0_b, Ti0_b, Te0_b, u0_b, Vpar0_b, Vpar0_bb, T0_bb, Ti0_bb, Te0_bb, u0_bb, R_b, Z_b, R_bb, Z_bb, ps0_bb, grad_b(2)
real*8  :: Btot, grad_psi, u0_s, u0_t, u0_x, u0_y
real*8  :: element_size_s, element_size_t, element_size_0, element_size_3
real*8  :: H1(2,n_degrees_1d), H1_s(2,n_degrees_1d), H1_ss(2,n_degrees_1d)
integer :: i, in, iv, iv2, iv3, inode, inode2, inode3, k
integer :: index_large_i, index_node, index_node2, index_node3, ielm
integer(kind=int_all) :: ijA_position,ijA_position2
integer :: ilarge2, kv, kT, kTi, kTe, ku, kn, ilarge_vv, ilarge_vT, ilarge_vus, ilarge_vn
integer :: ilarge_vsvs, ilarge_vsTs, ilarge_vsT, ilarge_vut, ilarge_vtvt, ilarge_vtTt, ilarge_vtT
integer :: ierr
logical :: apply_psi_BC, apply_current_BC, s_constant_boundary, t_constant_boundary, apply_cs, apply_dirichlet_1234, apply_dirichlet_all

real*8, allocatable :: psi_RMP_cos1(:),dpsi_RMP_cos_dR1(:),dpsi_RMP_cos_dZ1(:)
real*8, allocatable :: psi_RMP_sin1(:),dpsi_RMP_sin_dR1(:),dpsi_RMP_sin_dZ1(:)
real*8  :: Rnode, dRnode_ds, Znode, dZnode_ds, dRnode_dt, dZnode_dt, establish_RMP
real*8  :: delta_psi_rmp, delta_psi_rmp_dR, delta_psi_rmp_dZ, delta_psi_rmp_ds, delta_psi_rmp_dt, psi_test, sigmo_fonc
real*8  :: R_mid, Z_mid, R_center, Z_center, direction2, normal(2), normal_direction(2), grad_s(2), grad_t(2)
real*8  :: factor, factor_b, factor_bb, c_1, c_2, c_3, bn, dl, dl_b
real*8  :: cs0, cs0_T, cs0_TT, cs0_TTT
real*8  :: bn_b, bn_b_abs, hfact_b, hfact_bb, bn_1, bn_2, ps2_b, element_size_2
integer :: ilarge_vp, ilarge_vp2, bnd_type
integer :: kp, j, err, itest, i_mid, i_bnd, idir, iv_dir, iv_perp_dir, k_max
!> Prescribed floating-potential BC (bcs%floating_u). fu_var_T is the temperature trace variable
!! this build evolves: Te under WITH_TiTe, otherwise the single T - and floating_u_norm has already
!! halved Lambda for that case, so ONE coefficient covers both builds.
real*8  :: fu_a_n, fu_C_T, fu_C_V, fu_u, fu_T, fu_targ, fu_act
integer :: fu_var_T
!> Floating-u boundary diagnostic (floating_u_diag). Per boundary type, MAX/MIN only.
integer, parameter :: FD_NT = 12
real*8  :: fd_res_max(FD_NT), fd_vn_max(FD_NT), fd_pe_max(FD_NT)
real*8  :: fd_vin_max(FD_NT), fd_vout_max(FD_NT)
!> Magnitude of the two competing pieces of the Mach1 VALUE row, per boundary type:
!! the sound-speed term and the ExB drift compensation whose tangential derivative is
!! missing from the slope row on bicubic elements. Their ratio IS the inconsistency.
real*8  :: fd_m1cs_max(FD_NT), fd_m1dr_max(FD_NT)
real*8  :: m1_dr
!> One-sided, floored, clipped drift-compatible Bohm supplement (SOLPS BCMOM=13,
!! non-marginal branch, without extrapolation - see mach1_uout_supplement).
!! m1_D = R^2*u_b/psi_b (exact -vE.n/(Bn*|B|), Vpar units); m1_bfl = min(1,|bn|/s0)
!! with s0 = min_sheath_angle in radians (the c_angle scale) floors the incidence;
!! m1_S = 2*cs/Btot is the SOLPS bound; m1_sup is the applied supplement and
!! m1_act/m1_clw select its branch. u0_bb_r reconstructs the second tangential
!! derivative of u from BOTH endpoints' value/slope DOFs (same stencil as ps0_bb)
!! for the bicubic slope-row residual.
real*8  :: m1_D, m1_S, m1_bfl, m1_smin, m1_Dfl, m1_sup, m1_act, m1_clw
real*8  :: m1_dDdb, m1_dSdb, m1_dslope, m1_dfdb, m1_clamp, m1_raw, u0_bb_r
real*8  :: fd_rho_min(FD_NT), fd_T_min(FD_NT), fd_pe_R(FD_NT), fd_pe_Z(FD_NT)
real*8  :: fd_loc(2,FD_NT)
real*8  :: fd_es, fd_ep, fd_dl, fd_h, fd_res, fd_vn, fd_pe, fd_sq, fd_sgn
integer :: fd_owner, fd_t
integer :: n_rmp_harm, N_rmp_har_block_size

real*8  :: R_out, Z_out, s_elm, t_elm, QR,QR_s,QR_t,QR_st,QR_ss,QR_tt,QZ,QZ_s,QZ_t,QZ_st,QZ_ss,QZ_tt
real*8  :: QPs0,QPs0_s,QPs0_t,QPs0_st,QPs0_ss,QPs0_tt
integer :: ifail, i_elm

real*8  ::   Mach1BC,   Mach1BC_v,   Mach1BC_T,   Mach1BC_u
real*8  ::  dMach1BC,  dMach1BC_v,  dMach1BC_T,  dMach1BC_Ti, dMach1BC_Te,  dMach1BC_Tb, dMach1BC_ubb
real*8  :: d2Mach1BC, d2Mach1BC_v, d2Mach1BC_T, d2Mach1BC_Tb, d2Mach1BC_Tbb

integer :: node_indices( (n_order+1)/2, (n_order+1)/2 ), index_tmp, kk, ll
logical, parameter :: include_2nd_derivatives = .false.

RMPspectrum: if (RMP_on .and. (n_tor .ge. 3)) then !*****
  
! for the moment it's done in a way that all RMP harmonics follow each other,i.e. n=2,n=3,n=4... 
! if you want for example n=2 and n=4 RMP you should consider n=2,3,4, but put zeros at the boundary in the input file for n=3 RMP
! example: ntor=13 and nperiod=1(so taking into account, toroidal numbers n=0,1,2....6) and  n=2 and n=3 are toroidal numbers of RMPs, 
! so Number_RMP_harmonics=2, RMP_har_cos_spectrum(1)=4,RMP_har_sin_spectrum(1)=5,RMP_har_cos_spectrum(2)=6,RMP_har_sin_spectrum(2)=7.  
  
  call tr_allocate(psi_RMP_cos1,1, bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"psi_RMP_cos1",CAT_UNKNOWN)
  call tr_allocate(dpsi_RMP_cos_dR1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_cos_dR1",CAT_UNKNOWN)
  call tr_allocate(dpsi_RMP_cos_dZ1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_cos_dZ1",CAT_UNKNOWN)
  call tr_allocate(psi_RMP_sin1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"psi_RMP_sin1",CAT_UNKNOWN)
  call tr_allocate(dpsi_RMP_sin_dR1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_sin_dR1",CAT_UNKNOWN)
  call tr_allocate(dpsi_RMP_sin_dZ1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_sin_dZ1",CAT_UNKNOWN)
  
  N_rmp_har_block_size=bnd_node_list%n_bnd_nodes
    
  psi_test =  node_list%node(bnd_node_list%bnd_node(1)%index_jorek)%values(RMP_har_cos_spectrum(1),1,1)
  
  ! if necessary, replace by:
  ! psi_test =  node_list%node(bnd_node_list%bnd_node(1)%index_jorek)%values(min(RMP_har_cos_spectrum(1), n_tor),1,1)
  write (*,*) 'psi_bnd at previous time step', psi_test
    
  if (abs(psi_test) .le. abs(psi_RMP_cos(1))) then
    sigmo_fonc = ( 1.d0 + exp(-RMP_growth_rate*( t_now - RMP_start_time - RMP_ramp_up_time/2.d0 )))**(-1) &
               - ( 1.d0 + exp(-RMP_growth_rate*( 0.d0 - RMP_ramp_up_time/2.d0 )))**(-1) 
    establish_RMP = (RMP_growth_rate*sigmo_fonc*(1-sigmo_fonc)+1.e-6)*tstep 
  else
    establish_RMP = 0.d0
  endif
  ! Other possibility (simpler) : if ( (t_now - RMP_start_time) .ge. 2.2*RMP_ramp_up_time/2.d0 ) then establish_RMP =0.0
  
  do j=1, bnd_node_list%n_bnd_nodes*Number_RMP_harmonics  
    psi_RMP_cos1(j)     = psi_RMP_cos(j)     * establish_RMP
    dpsi_RMP_cos_dR1(j) = dpsi_RMP_cos_dR(j) * establish_RMP
    dpsi_RMP_cos_dZ1(j) = dpsi_RMP_cos_dZ(j) * establish_RMP
    psi_RMP_sin1(j)     = psi_RMP_sin(j)     * establish_RMP
    dpsi_RMP_sin_dR1(j) = dpsi_RMP_sin_dR(j) * establish_RMP
    dpsi_RMP_sin_dZ1(j) = dpsi_RMP_sin_dZ(j) * establish_RMP
  end do

  if (my_id == 0) then
    write (*,*) 'psi_RMP_cos1(1) and derivatives after multiplication in boundary conditions'
    write (*,*) psi_RMP_cos1(1), dpsi_RMP_cos_dR1(1), dpsi_RMP_cos_dZ1(1)
    write (*,*) 'establish_RMP', establish_RMP
  endif

end if RMPspectrum

zbig        = 1.d12
zbig_backup = zbig

! --- calculate node_indices
call calculate_node_indices(node_indices)

if ( floating_u_diag ) then
  fd_res_max = -1.d0 ; fd_vn_max = -1.d0 ; fd_pe_max = -1.d0
  fd_vin_max = 0.d0 ; fd_vout_max = 0.d0
  fd_m1cs_max = -1.d0 ; fd_m1dr_max = -1.d0
  fd_rho_min = huge(1.d0) ; fd_T_min = huge(1.d0)
  fd_pe_R = 0.d0 ; fd_pe_Z = 0.d0
endif

do i=1, n_local_elms !=== do elements

  ielm = local_elms(i)

  i_bnd = 0

  do iv=1, n_vertex_max 
    inode = element_list%element(ielm)%vertex(iv)
    if (node_list%node(inode)%boundary .ne. 0) i_bnd = i_bnd + 1
  enddo
  
  if (i_bnd .lt. 2) cycle           

  R_mid = 0.d0; Z_mid = 0.d0; R_center = 0.d0; Z_center = 0.d0
  
  iv2 = 0; iv3 = 0

  do iv=1, n_vertex_max ! check vertices for being a boundary point

    inode = element_list%element(ielm)%vertex(iv)

    if (node_list%node(inode)%boundary .eq. 0) cycle 

    do idir=1, 2        ! check the two directions

      R_mid = node_list%node(inode)%x(1,1,1)
      Z_mid = node_list%node(inode)%x(1,1,2)

      if (idir .eq. 1) then
        iv2 = mod(iv  ,4) + 1
        iv3 = mod(iv+2,4) + 1
      else
        iv2 = mod(iv+2,4) + 1
        iv3 = mod(iv  ,4) + 1
      endif

      inode2 = element_list%element(ielm)%vertex(iv2)
      inode3 = element_list%element(ielm)%vertex(iv3)

      if (node_list%node(inode2)%boundary .eq. 0) cycle

      if ((iv*iv2 .eq. 2) .or. (iv*iv2 .eq. 12)) then
        s_constant_boundary = .false.
        t_constant_boundary = .true.
        iv_dir      = 2
        iv_perp_dir = 3
      elseif ((iv*iv2 .eq. 6) .or. (iv*iv2 .eq. 4)) then
        s_constant_boundary = .true.
        t_constant_boundary = .false.
        iv_dir      = 3
        iv_perp_dir = 2
      else
        write(*,*) 'THIS SHOULD NOT BE POSSIBLE'
      endif


      R_center = node_list%node(inode3)%x(1,1,1)
      Z_center = node_list%node(inode3)%x(1,1,2)

      normal_direction = (/R_mid - R_center, Z_mid - Z_center /) / norm2((/R_mid - R_center, Z_mid - Z_center /))

      bnd_type = node_list%node(inode)%boundary
      fu_var_T = var_T
      if ( with_TiTe ) fu_var_T = var_Te

      ! --- FLOATING-U BOUNDARY DIAGNOSTIC ------------------------------------------
      ! --- Per boundary type, at the AXISYMMETRIC harmonic only:
      ! ---   fd_res  = |u - C_T*T - C_V*V_wall|   trace residual (the acceptance gate)
      ! ---   fd_vn   = |v_E.n| = |R * du/dl|      normal ExB speed the BC imposes
      ! ---   fd_pe   = |v_E.n| * h_perp / D_perp  cell Peclet number of that flow
      ! ---   min rho, min T on the boundary trace
      ! ---
      ! --- v_E.n = -R*du/dl for the candidate normal below. For v_E = -R grad(u) x e_phi in
      ! --- the right-handed (e_R,e_Z,e_phi) basis, with edge tangent t = (R_b,Z_b)/dl
      ! --- and normal n = (Z_b,-R_b)/dl,
      ! ---     v_E.n = -R*(du/dR * R_b + du/dZ * Z_b)/dl = -R * u_b/dl.
      ! --- So the normal ExB flow is driven by the TANGENTIAL derivative of u - which
      ! --- is exactly what u = C_T*Te manufactures wherever Te varies along the wall.
      ! ---
      ! --- Only MAX and MIN are accumulated. Both are idempotent, so a halo element
      ! --- visited by several ranks cannot inflate them the way an area integral would.
      ! --- The element_size sign flips applied further down cancel in |u_b|/|dl|, so no
      ! --- orientation convention enters any magnitude reported here.
      if ( floating_u_diag .and. (bnd_type .ge. 1) .and. (bnd_type .le. FD_NT) ) then
      if ( bcs(bnd_type)%floating_u .and.                                              &
           (a_mat%i_tor_min .le. 1) .and. (a_mat%i_tor_max .ge. 1) ) then
        call floating_u_norm(fu_a_n, fu_C_T, fu_C_V)
        fd_es = element_list%element(ielm)%size(iv, iv_dir)      * H1_s(1,2)
        fd_ep = element_list%element(ielm)%size(iv, iv_perp_dir) * H1_s(1,2)
        fd_dl = sqrt( (node_list%node(inode)%x(1,iv_dir,1)      * fd_es)**2            &
                    + (node_list%node(inode)%x(1,iv_dir,2)      * fd_es)**2 )
        fd_h  = sqrt( (node_list%node(inode)%x(1,iv_perp_dir,1) * fd_ep)**2            &
                    + (node_list%node(inode)%x(1,iv_perp_dir,2) * fd_ep)**2 )
        fd_res = abs( node_list%node(inode)%values(1,1,var_u)                           &
                      - fu_C_T * node_list%node(inode)%values(1,1,fu_var_T)             &
                      - fu_C_V * sheath_V_wall )
        ! --- SIGNED v_E.n. The identity v_E.n = -R*du/dl pairs the edge tangent
        ! --- t = (R_b,Z_b)/dl with the normal n = (Z_b,-R_b)/dl. Which way that n
        ! --- points depends on the element's vertex ordering and on the element_size
        ! --- sign flips, so it is NOT an outward normal by construction. Orient it
        ! --- against normal_direction, which this loop already builds as node-minus-
        ! --- opposite-vertex and therefore genuinely points out of the domain.
        ! --- fd_vn > 0 is ExB outflow, fd_vn < 0 is ExB inflow. Only the TOTAL
        ! --- normal plasma velocity determines advective particle inflow.
        fd_vn = 0.d0
        if ( fd_dl .gt. 0.d0 ) then
          fd_vn = -node_list%node(inode)%x(1,1,1)                                       &
                  * node_list%node(inode)%values(1,iv_dir,var_u) * fd_es / fd_dl
          fd_sgn =   ( node_list%node(inode)%x(1,iv_dir,2)*fd_es) * normal_direction(1)  &
                   - ( node_list%node(inode)%x(1,iv_dir,1)*fd_es) * normal_direction(2)
          if ( fd_sgn .lt. 0.d0 ) fd_vn = -fd_vn
        endif
        fd_vout_max(bnd_type) = max( fd_vout_max(bnd_type),  fd_vn )
        fd_vin_max (bnd_type) = max( fd_vin_max (bnd_type), -fd_vn )
        fd_vn = abs(fd_vn)
        fd_pe = 0.d0
        if ( D_perp(1) .gt. 0.d0 ) fd_pe = fd_vn * fd_h / D_perp(1)
        fd_res_max(bnd_type) = max( fd_res_max(bnd_type), fd_res )
        fd_vn_max (bnd_type) = max( fd_vn_max (bnd_type), fd_vn  )
        fd_rho_min(bnd_type) = min( fd_rho_min(bnd_type),                               &
                                    node_list%node(inode)%values(1,1,var_rho) )
        fd_T_min  (bnd_type) = min( fd_T_min  (bnd_type),                               &
                                    node_list%node(inode)%values(1,1,fu_var_T) )
        if ( fd_pe .gt. fd_pe_max(bnd_type) ) then
          fd_pe_max(bnd_type) = fd_pe
          fd_pe_R  (bnd_type) = node_list%node(inode)%x(1,1,1)
          fd_pe_Z  (bnd_type) = node_list%node(inode)%x(1,1,2)
        endif
      endif
      endif

      do in=a_mat%i_tor_min, a_mat%i_tor_max  ! === do n_tor
      
        if (keep_n0_const  .and.  in .eq. 1 ) then
          zbig = 1.d15
        else
          zbig = zbig_backup
        endif

!================================= start RMPs (for both directions) ==================================
        if (RMP_on ) then

          do n_rmp_harm=1, Number_RMP_harmonics !=== do RMP harmonics

            if (((in.eq.RMP_har_cos_spectrum(n_rmp_harm)) .or. (in.eq.RMP_har_sin_spectrum(n_rmp_harm))) &
               .and. (.not. freeboundary)) then
                   ! in .eq. RMP_har_cos corresponds to cos(n_perturbation)
                   ! in .eq. RMP_har_sin corresponds to sin(n_perturbation)
                                 
              index_node = node_list%node(inode)%index(1)  !=== index in RHS (or matrix A not compressed)
                                          
              Rnode     = node_list%node(inode)%x(1,1,1) 
              dRnode_ds = node_list%node(inode)%x(1,iv_dir,1) 
              Znode     = node_list%node(inode)%x(1,1,2) 
              dZnode_ds = node_list%node(inode)%x(1,iv_dir,2) 
                  
              if (in.eq.RMP_har_cos_spectrum(n_rmp_harm)) then
                delta_psi_rmp = psi_RMP_cos1(node_list%node(inode)%boundary_index +N_rmp_har_block_size*(n_rmp_harm-1))
                delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
              else 
                delta_psi_rmp = psi_RMP_sin1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
              endif
                  
              delta_psi_rmp_ds = delta_psi_rmp_dR * dRnode_ds + delta_psi_rmp_dZ * dZnode_ds

              call boundary_conditions_add_one_entry(                &
                     index_node, var_psi, in, index_node, var_psi, in,         &
                     zbig, index_min, index_max, a_mat)

              call boundary_conditions_add_RHS(                      &
                     index_node, var_psi, in, index_min, index_max,       &
                     RHS_loc, ZBIG * delta_psi_rmp, a_mat%i_tor_min, a_mat%i_tor_max)
                  
              index_node2 = node_list%node(inode)%index(iv_dir)

              call boundary_conditions_add_one_entry(                 &
                     index_node2, var_psi, in, index_node2, var_psi, in,        &
                     zbig, index_min, index_max, a_mat)

              call boundary_conditions_add_RHS(                       &
                     index_node2, var_psi, in, index_min, index_max,       &
                     RHS_loc, ZBIG * delta_psi_rmp_ds, a_mat%i_tor_min, a_mat%i_tor_max)

            endif !=== endif selection RMP harmonics
        
          enddo   !=== enddo RMP harmonics   
        
        endif     !=== endif RMP

 
        do k=1, n_var ! === do variables
                                                                                                 
          !------------ Decide when Psi or Current need BCs --------------------------------------------------                      
          !----Psi
          apply_psi_BC = .false.
          if (k == var_psi) then                        
            if ( (RMP_on) .and. (in .lt. RMP_har_cos_spectrum(1))                    )   apply_psi_BC = .true.
            if ( (RMP_on) .and. (in .gt. RMP_har_sin_spectrum(Number_RMP_harmonics)) )   apply_psi_BC = .true.
            if ( (.not. RMP_on) .and. (in .ge. 2)              )                         apply_psi_BC = .true.
            if (in .eq. 1)                                                               apply_psi_BC = .true.
            if (is_freebound(in,k))                                                      apply_psi_BC = .false.                     
          endif
                
          !----Current
          apply_current_BC = .false.
          if (k == var_zj) then
            if ( .not. is_freebound(in,k) )   apply_current_BC = .true.
          endif
          !---------------------------------------------------------------------------------------------------                      

          
          !------------ Decide when to apply vpar=cs ---------------------------------------------------------                      
          apply_cs = .false.          
          if ( (.not. mach_one_bnd_integral) .and. bcs(bnd_type)%mach1 .and. with_vpar) then
            apply_cs = .true.
          endif
          if (floating_u_mach_flux .and. bcs(bnd_type)%floating_u) apply_cs = .false.
          !---------------------------------------------------------------------------------------------------                      

          if (  ( (k == var_psi     ) .and. bcs(bnd_type)%dirichlet%psi     )  .or.  &
                ( (k == var_u       ) .and. bcs(bnd_type)%dirichlet%u       )  .or.  &
                ( (k == var_zj      ) .and. bcs(bnd_type)%dirichlet%zj      )  .or.  &
                ( (k == var_w       ) .and. bcs(bnd_type)%dirichlet%w       )  .or.  &
                ( (k == var_rho     ) .and. bcs(bnd_type)%dirichlet%rho     )  .or.  &
                ( (k == var_T       ) .and. bcs(bnd_type)%dirichlet%T       )  .or.  &
                ( (k == var_Ti      ) .and. bcs(bnd_type)%dirichlet%Ti      )  .or.  &
                ( (k == var_Te      ) .and. bcs(bnd_type)%dirichlet%Te      )  .or.  &
                ( (k == var_Vpar    ) .and. bcs(bnd_type)%dirichlet%Vpar    )  .or.  &
                ( (k == var_rhon    ) .and. bcs(bnd_type)%dirichlet%rhon    )  .or.  &
                ( (k == var_rhoimp  ) .and. bcs(bnd_type)%dirichlet%rho_imp )  .or.  &
                ( (k == var_nre     ) .and. bcs(bnd_type)%dirichlet%nre     )        &
             ) then

            ! --- If special conditions apply (e.g. freeboundary, mach1), do not apply Dirichlet even if specified in the namelist
            if ( (k==var_psi  ) .and. (.not. apply_psi_BC    ) )       cycle
            if ( (k==var_zj   ) .and. (.not. apply_current_BC) )       cycle
            if ( (k==var_vpar ) .and.  apply_cs .and. (bnd_type/=3)  ) cycle  ! vpar=cs is a special case (this is done below)
                                                                              ! however bnd_type=3 needs both BCs for different directions
            if (k==var_vpar .and. floating_u_mach_flux .and. bcs(bnd_type)%floating_u &
                .and. bcs(bnd_type)%mach1) cycle ! quadrature flow constraint supplies this trace

!            if ((k.eq.7) .and. (node_list%node(inode)%boundary .eq. 3)) cycle  !=== better included for ITER extended wall

            ! --- Fix derivatives in one direction
            do kk = 1,(n_order+1)/2
              if ( (iv_dir .eq. 3) .and. (kk .gt. 1) ) cycle ! do only t-derivatives and node value
              do ll = 1,(n_order+1)/2
                if ( (iv_dir .eq. 2) .and. (ll .gt. 1) ) cycle ! do only s-derivatives and node value
                index_tmp = node_indices(kk,ll)
                index_node = node_list%node(inode)%index(index_tmp)
                call boundary_conditions_add_one_entry(                 &
                       index_node, k, in, index_node, k, in,            &
                       zbig, index_min, index_max, a_mat)

                ! --- PRESCRIBED FLOATING POTENTIAL, V_p - V_wall = Lambda*k_B*Te/e.
                ! --- The entry above is the ordinary Dirichlet diagonal; here it is completed
                ! --- into the affine constraint
                ! ---
                ! ---     u - C_T*Te - C_V*V_wall = 0
                ! ---
                ! --- by adding the Te cross-column and the residual RHS, so the matrix
                ! --- conditioning differs from plain u = 0 only by that one extra column.
                ! ---
                ! --- Written on exactly the trace DOFs the enclosing loop already selects - the
                ! --- value and the pure TANGENTIAL derivatives - because differentiating a
                ! --- constant-coefficient relation along the boundary is exact DOF by DOF:
                ! ---     value:      u_val = C_T*corr(Te_val) + C_V*V_wall
                ! ---     tangential: u_d   = C_T*corr'(Te_val)*Te_d
                ! --- Both are EXACT, not linearised: the Jacobian entries are constants and no
                ! --- term is neglected. V_wall enters the value equation only, and only the
                ! --- axisymmetric harmonic.
                ! ---
                ! --- Corners are safe by construction, not by ordering: the row is idempotent, a
                ! --- node shared by two floating edges is visited once per incident edge, each
                ! --- visit constrains that edge's own tangential DOF, and the value row is
                ! --- written identically both times.
                ! ---
                ! --- Nothing here forms B.n, divides by the element Jacobian, or converts
                ! --- logical derivatives into (R,Z). u and Te share the same nodal frame, so its
                ! --- scaling cancels identically between the two terms - which is why one closure
                ! --- serves every wall boundary type.
                ! --- RAW Te, NOT a positivity-mapped Te. Two reasons, and the second is fatal
                ! --- to the alternative:
                ! ---
                ! ---  1. A NONLINEAR MAP CANNOT BE APPLIED COEFFICIENT BY COEFFICIENT. `in` is a
                ! ---     toroidal HARMONIC index, not a physical location, so corr_neg(Te_n) is
                ! ---     meaningless for n /= 0: a ZERO non-axisymmetric temperature coefficient
                ! ---     maps to corr_neg(0) = 0.84 eV and would inject a spurious ~2.5 V
                ! ---     potential harmonic out of nothing. Harmless at n_tor = 1, where in = 1
                ! ---     IS the physical value, and fundamentally wrong in 3D.
                ! ---  2. With the map, g = u - C_T*f(Te) - C_V*Vw has the temperature column
                ! ---     -C_T*f'(Te), the tangential relation acquires an f''*Te_l column on the
                ! ---     VALUE DOF, and neither is constant. The relation stops being affine and
                ! ---     the "exact, not linearised" claim stops being true.
                ! ---
                ! --- Using the evolved Te makes u = C_T*Te + C_V*Vw exactly affine, exact mode by
                ! --- mode, with the constant Jacobian (1, -C_T) and no cross-DOF terms. Negative
                ! --- temperatures are the global positivity scheme's job; repairing them
                ! --- nonlinearly inside a spectral boundary condition is worse than the disease.
                if ( (k == var_u) .and. bcs(bnd_type)%floating_u ) then
                  call floating_u_norm(fu_a_n, fu_C_T, fu_C_V)
                  ! --- TEMPERATURE FLOOR, matching the rest of the boundary path.
                  ! --- mod_boundary_conditions:667/670 (legacy Mach) and mod_mach1_trace
                  ! --- both use max(T, T_min); this row used the RAW value and was the
                  ! --- only thing in the boundary path that did not.
                  ! ---
                  ! --- Measured consequence: with raw Te, a wall node whose Te crosses
                  ! --- zero imposes a NEGATIVE potential and therefore a REVERSED ExB
                  ! --- drift, while cs at the same node stays floored - the two halves of
                  ! --- the boundary condition disagreeing about the temperature. On the
                  ! --- weighted run, type-1 min Te sat at +0.64 eV for 51 outputs, crossed
                  ! --- zero, and the physical Mach residual |G|/cs began growing on the
                  ! --- very next output while it had been flat to 2 % before.
                  ! ---
                  ! --- With the floor, a clamped node has u = C_T*T_min, which is CONSTANT,
                  ! --- so its tangential derivative is zero and it drives no ExB flow at
                  ! --- all - rather than driving one the wrong way.
                  ! ---
                  ! --- The branch is decided by the AXISYMMETRIC VALUE, as at :667/670, and
                  ! --- fu_act is the exact derivative of max() on that branch - not a
                  ! --- smooth positivity map. The relation therefore stays piecewise
                  ! --- affine and is still imposed exactly by the single linear solve
                  ! --- whenever the branch does not change, which is the same caveat the
                  ! --- legacy Mach rows already carry. Applying a floor to a non
                  ! --- axisymmetric COEFFICIENT would be meaningless, so the harmonics and
                  ! --- the derivative DOFs carry fu_act linearly instead.
                  fu_act = 1.d0
                  if ( node_list%node(inode)%values(1,1,fu_var_T) .le. T_min ) fu_act = 0.d0
                  fu_u    = node_list%node(inode)%values(in, index_tmp, var_u)
                  fu_T    = node_list%node(inode)%values(in, index_tmp, fu_var_T)
                  fu_targ = fu_C_T * fu_act * fu_T
                  if ( (index_tmp .eq. 1) .and. (in .eq. 1) .and. (fu_act .eq. 0.d0) ) &
                    fu_targ = fu_C_T * T_min
                  ! --- V_wall is a constant, so it enters the VALUE equation only, and only the
                  ! --- axisymmetric harmonic. Every derivative equation is homogeneous.
                  if ( (index_tmp .eq. 1) .and. (in .eq. 1) ) &
                    fu_targ = fu_targ + fu_C_V * sheath_V_wall
                  call boundary_conditions_add_one_entry(                     &
                         index_node, var_u, in, index_node, fu_var_T, in,     &
                         - zbig * fu_C_T * fu_act, index_min, index_max, a_mat)
                  call boundary_conditions_add_RHS(                           &
                         index_node, var_u, in, index_min, index_max, RHS_loc,&
                         - zbig * ( fu_u - fu_targ ),                         &
                         a_mat%i_tor_min, a_mat%i_tor_max)
                endif
              enddo
            enddo
            

          endif
          

          if ( (.not. is_freebound(in,k)) ) then
            if ( ( loop_voltage .ne. 0.d0 ) ) then
              if ( (k == var_psi) .and. (in == 1) ) then
                index_node = node_list%node(inode)%index(1)
                call boundary_conditions_add_RHS(       &
                          index_node, k, in,     &
                          index_min, index_max,  &
                          RHS_loc, ZBIG*(loop_voltage*sqrt(MU_ZERO*central_density*central_mass*ATOMIC_MASS_UNIT*1.d20))* tstep, &
                          a_mat%i_tor_min, a_mat%i_tor_max)
              endif
            endif
          endif

        enddo !=== variables

        if ((node_list%node(inode)%boundary .eq.  3) .and. (node_list%node(inode2)%boundary .eq.  2)) cycle

        ! --- Mach1 Boundary Conditions
        if ( apply_cs ) then

          call basisfunctions1(0.d0, H1, H1_s, H1_ss)

          element_size_s = element_list%element(ielm)%size(iv,2) * H1_s(1,2)
          element_size_t = element_list%element(ielm)%size(iv,3) * H1_s(1,2)

          if ((iv .eq. 2) .or. (iv .eq. 3))  element_size_s = - element_size_s 
          if ((iv .eq. 3) .or. (iv .eq. 4))  element_size_t = - element_size_t 

          element_size_0 =   element_list%element(ielm)%size(iv, iv_dir) * H1_s(1,2) 
          element_size_2 =   element_list%element(ielm)%size(iv2,iv_dir) * H1_s(1,2) 

          if (t_constant_boundary .and. ((iv  .eq. 2) .or. (iv  .eq. 3))) element_size_0 = - element_size_0
          if (s_constant_boundary .and. ((iv  .eq. 3) .or. (iv  .eq. 4))) element_size_0 = - element_size_0

          if (t_constant_boundary .and. ((iv2 .eq. 2) .or. (iv2 .eq. 3))) element_size_2 = - element_size_2
          if (s_constant_boundary .and. ((iv2 .eq. 3) .or. (iv2 .eq. 4))) element_size_2 = - element_size_2
          
          if (n_order .ge. 5) then
            element_size_3 =  element_list%element(ielm)%size(iv, iv_dir+3) * H1_ss(1,2) 
            if (t_constant_boundary .and. ((iv  .eq. 2) .or. (iv  .eq. 3))) element_size_3 = - element_size_3
            if (s_constant_boundary .and. ((iv  .eq. 3) .or. (iv  .eq. 4))) element_size_3 = - element_size_3
          endif

          index_node    = node_list%node(inode)%index(1)             ! position of value
          index_node2   = node_list%node(inode)%index(iv_dir)        ! position of first deriative
          if (n_order .ge. 5) &
            index_node3 = node_list%node(inode)%index(iv_dir+3)      ! position of 2nd deriative

          ! --- Determine the direction of the BCs and apply smoothing factors if requested
          ps0       = node_list%node(inode)%values(1,1,var_psi)
          ps0_b     = node_list%node(inode)%values(1,iv_dir,var_psi)  * element_size_0 
          ps0_s     = node_list%node(inode)%values(1,2,var_psi)       * element_size_s
          ps0_t     = node_list%node(inode)%values(1,3,var_psi)       * element_size_t

          BigR      = node_list%node(inode)%x(1,1,1)
          R_b       = node_list%node(inode)%x(1,iv_dir,1) * element_size_0
          Z_b       = node_list%node(inode)%x(1,iv_dir,2) * element_size_0

          R_s       = node_list%node(inode)%x(1,2,1)      * element_size_s
          R_t       = node_list%node(inode)%x(1,3,1)      * element_size_t    
          Z_s       = node_list%node(inode)%x(1,2,2)      * element_size_s
          Z_t       = node_list%node(inode)%x(1,3,2)      * element_size_t    
          Z         = node_list%node(inode)%x(1,1,2)
          
          ps0_bb = element_list%element(ielm)%size(iv ,1)      * node_list%node(inode )%values(1,1,var_psi)      * H1_ss(1,1) &
                 + element_list%element(ielm)%size(iv ,iv_dir) * node_list%node(inode )%values(1,iv_dir,var_psi) * H1_ss(1,2) &
                 + element_list%element(ielm)%size(iv2,1)      * node_list%node(inode2)%values(1,1,var_psi)      * H1_ss(2,1) &
                 + element_list%element(ielm)%size(iv2,iv_dir) * node_list%node(inode2)%values(1,iv_dir,var_psi) * H1_ss(2,2)

          ! --- Second tangential derivative of u along the edge, reconstructed from
          ! --- both endpoints' value/slope DOFs exactly as ps0_bb above. Needed by the
          ! --- slope row's drift residual; no second-derivative nodal DOF exists at
          ! --- bicubic order.
          u0_bb_r= element_list%element(ielm)%size(iv ,1)      * node_list%node(inode )%values(1,1,var_u)        * H1_ss(1,1) &
                 + element_list%element(ielm)%size(iv ,iv_dir) * node_list%node(inode )%values(1,iv_dir,var_u)   * H1_ss(1,2) &
                 + element_list%element(ielm)%size(iv2,1)      * node_list%node(inode2)%values(1,1,var_u)        * H1_ss(2,1) &
                 + element_list%element(ielm)%size(iv2,iv_dir) * node_list%node(inode2)%values(1,iv_dir,var_u)   * H1_ss(2,2)

          R_bb = + element_list%element(ielm)%size(iv ,1)      * node_list%node(inode )%x(1,1,1)      * H1_ss(1,1)  &
                 + element_list%element(ielm)%size(iv ,iv_dir) * node_list%node(inode )%x(1,iv_dir,1) * H1_ss(1,2)  &
                 + element_list%element(ielm)%size(iv2,1)      * node_list%node(inode2)%x(1,1,1)      * H1_ss(2,1)  &
                 + element_list%element(ielm)%size(iv2,iv_dir) * node_list%node(inode2)%x(1,iv_dir,1) * H1_ss(2,2)  

          Z_bb = + element_list%element(ielm)%size(iv ,1)      * node_list%node(inode )%x(1,1,     2) * H1_ss(1,1)  &
                 + element_list%element(ielm)%size(iv ,iv_dir) * node_list%node(inode )%x(1,iv_dir,2) * H1_ss(1,2)  &
                 + element_list%element(ielm)%size(iv2,1)      * node_list%node(inode2)%x(1,1,     2) * H1_ss(2,1)  &
                 + element_list%element(ielm)%size(iv2,iv_dir) * node_list%node(inode2)%x(1,iv_dir,2) * H1_ss(2,2)  

          ps2_b     = node_list%node(inode2)%values(1,iv_dir,var_psi) * element_size_2 
          
          xjac  =  R_s*Z_t - R_t*Z_s
          ps0_x = (   Z_t * ps0_s - Z_s * ps0_t ) / xjac
          ps0_y = ( - R_t * ps0_s + R_s * ps0_t ) / xjac

          grad_s = (/  Z_t,  -R_t /) / xjac
          grad_t = (/ -Z_s,   R_s /) / xjac

          if (s_constant_boundary) then
            grad_b = grad_s
          elseif (t_constant_boundary) then
            grad_b = grad_t
          endif

          normal     = dot_product(grad_b,normal_direction) * grad_b      ! outward pointing normal
          normal     = normal / norm2(normal)
          direction  = sign(1.d0,dot_product((/ps0_y,-ps0_x/),normal))
          
          grad_psi = sqrt(ps0_x**2 + ps0_y**2)
          Btot     = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR
          dl       = sqrt(R_b**2 + Z_b**2)
          dl_b     = (R_b*R_bb + Z_b*Z_bb) / dl

          bn     = dot_product( (/ps0_y,-ps0_x/), normal ) /  (BigR*Btot)  ! B�n/Btot
          bn_b   = 1.d0 / (Btot*dl*BigR) * (ps0_bb - ps0_b * dl_b /dl )

          bn_1 = + ps0_b/(BigR*Btot*dl)
          bn_2 = + ps2_b/(BigR*Btot*dl)

          if ((s_constant_boundary) .and. ((iv  .eq. 1) .or. (iv  .eq. 4))) bn_1 = - bn_1
          if ((t_constant_boundary) .and. ((iv  .eq. 3) .or. (iv  .eq. 4))) bn_1 = - bn_1
          if ((s_constant_boundary) .and. ((iv2 .eq. 1) .or. (iv2 .eq. 4))) bn_2 = - bn_2
          if ((t_constant_boundary) .and. ((iv2 .eq. 3) .or. (iv2 .eq. 4))) bn_2 = - bn_2

          ! --- Apply Smoothing?
          c_1 = vpar_smoothing_coef(1); c_2 = vpar_smoothing_coef(2); c_3 = vpar_smoothing_coef(3)
          if ((vpar_smoothing) .and. (bn_1*bn_2 .lt. 0.d0)) then
            if (c_2 .gt. 0d0) then
              factor    = 0.25d0 * ( 1.d0 + tanh( (abs(bn) - c_1) / c_2 ) )**2 - c_3
              factor_b  = 0.5d0  * ( 1.d0 + tanh( (abs(bn) - c_1) / c_2 ) )           & 
                        * (bn_b * bn/abs(bn) /c_2) /(cosh( (abs(bn) - c_1) / c_2 ) )**2
              factor_bb = 0.5d0 * (bn_b * bn/abs(bn) /c_2)**2 / (cosh( (abs(bn) - c_1) / c_2 ) )**4 &
                         - 1.0d0  * ( 1.d0 + tanh( (abs(bn) - c_1) / c_2 ) ) &
                        * (bn_b * bn/abs(bn) /c_2)**2 /(cosh( (abs(bn) - c_1) / c_2 ) )**3 * sinh( (abs(bn) - c_1) / c_2 )
             else
               factor    = tanh(bn/c_1)
               factor_b  = bn_b /c_1 / cosh(bn/c_1)**2
               factor_bb = -2.d0 * (bn_b/c_1)**2 / cosh(bn/c_1)**3 * sinh(bn/c_1)
               direction = 1.d0                            
            endif                       
          else
            factor    = 1.d0
            factor_b  = 0.d0
            factor_bb = 0.d0
          endif
          Hfact_b   = factor * R_b / BigR  + factor_b
          Hfact_bb  = factor * R_bb/ BigR - factor * R_b**2 / BigR**2  + factor_bb

          ! --- For the BC's the magnetic field is assumed constant, and we do
          ! --- not take derivatives of psi into account. Vpar, T and U are the
          ! --- only variables that will be used for the linearisation
          if ( with_TiTe ) then
            Ti0        = max(node_list%node(inode)%values(1,1,var_Ti), T_min)
            Ti0_b      = node_list%node(inode)%values(1,iv_dir,var_Ti)    * element_size_0 

            Te0        = max(node_list%node(inode)%values(1,1,var_Te), T_min)
            Te0_b      = node_list%node(inode)%values(1,iv_dir,var_Te)    * element_size_0

            T0   = Ti0 + Te0
            T0_b = Ti0_b + Te0_b
          else
            T0        = max(node_list%node(inode)%values(1,1,var_T), T_min)
            T0_b      = node_list%node(inode)%values(1,iv_dir,var_T)    * element_size_0

            Ti0   = T0   / 2.d0
            Ti0_b = T0_b / 2.d0
            Te0   = T0   / 2.d0
            Te0_b = T0_b / 2.d0
          end if

          Vpar0     = node_list%node(inode)%values(1,1,var_vpar)
          Vpar0_b   = node_list%node(inode)%values(1,iv_dir,var_Vpar) * element_size_0 

          u0_b      = node_list%node(inode)%values(1,iv_dir,var_u)    * element_size_0 

          if (n_order .ge. 5) then
            if ( with_TiTe ) then
              Ti0_bb    = node_list%node(inode)%values(1,iv_dir+3,var_Ti)   * element_size_3 
              Te0_bb    = node_list%node(inode)%values(1,iv_dir+3,var_Te)   * element_size_3 
              T0_bb     = Ti0_bb + Te0_bb
            else
              T0_bb     = node_list%node(inode)%values(1,iv_dir+3,var_T)    * element_size_3 
              Ti0_bb = T0_bb / 2.d0
              Te0_bb = T0_bb / 2.d0
            endif
            Vpar0_bb  = node_list%node(inode)%values(1,iv_dir+3,var_Vpar) * element_size_3 
            u0_bb     = node_list%node(inode)%values(1,iv_dir+3,var_u)    * element_size_3 
          endif

          ! --- Mach1 BC's and derivatives
          cs0      =   sqrt(gamma*(Ti0+Te0))
          cs0_T    =   0.5d0  * gamma    / cs0
          cs0_TT   = - 0.25d0 * gamma**2 / cs0**3 
          cs0_TTT  = 3.d0/8.d0* gamma**3 / cs0**5 

          ! --- DRIFT-COMPATIBLE BOHM SUPPLEMENT (SOLPS BCMOM=13, NON-MARGINAL branch -
          ! --- the one the manual marks "recommended for cases with drifts" - without
          ! --- the interior extrapolation a nodal row cannot have):
          ! ---
          ! ---     Vpar = direction*cs/Btot + supplement,   supplement >= 0 outward
          ! ---
          ! --- The supplement is the component of the kinematic cancellation
          ! --- D = R^2*u_b/psi_b = -vE.n/(Bn*|B|) that INCREASES the outward parallel
          ! --- flow, i.e. it compensates INWARD ExB drift so the target stays a net
          ! --- outflow boundary. Outward drift leaves Vpar at exactly sonic: the Bohm
          ! --- condition is an inequality, extra outward flux is allowed, and the
          ! --- recommended SOLPS branch never demands subsonic or reversed parallel
          ! --- flow. (The MARGINAL branch - Vpar = cs - vE.n/bn down to full reversal -
          ! --- was implemented first and crashed at 466 steps; outward-drift limit of
          ! --- the recommended branch is exactly the plain sonic row, which is the
          ! --- configuration measured stable.)
          ! ---
          ! --- Three safeguards, each with a measured failure behind it:
          ! ---  1. INCIDENCE FLOOR m1_bfl = min(1,|bn|/s0), s0 = min_sheath_angle in
          ! ---     radians - EXACTLY the c_angle scale that already floors the sheath
          ! ---     particle and heat fluxes (mod_boundary_matrix_open:93), applied to
          ! ---     the momentum channel with the same meaning: below s0 the sheath-
          ! ---     entrance model is floor-dominated and the effective collection
          ! ---     angle saturates (magnetic presheath / finite Larmor radius). This
          ! ---     bounds every Jacobian column by ~1/s0 instead of 1/bn - the
          ! ---     unbounded intra-solve column was how the marginal form died - and
          ! ---     widens the grazing response band from 2*cs*bn (a step function in
          ! ---     u) to 2*cs*s0 (a resolvable ramp). NOT a new threshold: the code's
          ! ---     existing definition of grazing, third use. bn is frozen in time
          ! ---     (psi is Dirichlet on the wall), so the floor is a static map.
          ! ---  2. SOLPS CLIP at 2*cs/Btot (manual 3.0.9 p.407/411, b2stbc_cbc=1.0),
          ! ---     stated in cs units. Hard, not tanh: piecewise-linear branches are
          ! ---     exact within one frozen solve; the smooth version degenerated to a
          ! ---     fixed-point iteration (period-2, crash at 319).
          ! ---  3. ONE-SIDED and anchored at zero, so wall points whose potential
          ! ---     gradient wobbles around zero produce nothing.
          ! ---
          ! --- With this row active the sheath ExB heat flux (mod_boundary_matrix_open)
          ! --- stays consistent: total collection remains >= cs*|bn| and bounded.
          ! ---
          ! --- m1_dr = 0 (mach1_omit_drift) remains the A/B control: plain Vpar = cs.
          m1_dr = 1.d0
          if ( mach1_omit_drift ) m1_dr = 0.d0

          ! --- Same convention as c_angle (radians of min_sheath_angle); bn is the
          ! --- incidence SINE - identical to within 5e-5 at 1 degree. A non-positive
          ! --- angle disables the floor (m1_bfl = 1), never the supplement.
          m1_smin = min_sheath_angle * PI / 180.d0
          m1_bfl  = 1.d0
          if ( m1_smin .gt. 0.d0 ) m1_bfl = min( 1.d0, abs(bn)/m1_smin )

          m1_D   = BigR**2 * U0_b / ps0_b
          m1_Dfl = m1_D * m1_bfl
          m1_S   = 2.d0 * cs0 / Btot
          call mach1_uout_supplement(direction*m1_Dfl, m1_S, m1_sup, m1_act, m1_clw)

          Mach1BC     = - Vpar0   + direction / Btot * factor  * cs0     + m1_dr * direction * m1_sup
          Mach1BC_v   = - 1.0
          Mach1BC_T   =           + direction / Btot * factor  * cs0_T                  &
                                  + m1_dr * m1_clw * direction * 2.d0 * cs0_T / Btot
          Mach1BC_u   =             m1_dr * m1_act * m1_bfl * BigR**2 * element_size_0/ps0_b

          ! --- Report the APPLIED supplement and the sonic term. Their ratio is
          ! --- bounded by 2/factor by construction - the acceptance check.
          if ( floating_u_diag .and. (bnd_type .ge. 1) .and. (bnd_type .le. FD_NT) ) then
          if ( bcs(bnd_type)%floating_u ) then
            fd_m1cs_max(bnd_type) = max( fd_m1cs_max(bnd_type),                        &
                                         abs( direction / Btot * factor * cs0 ) )
            fd_m1dr_max(bnd_type) = max( fd_m1dr_max(bnd_type), abs(m1_sup) )
          endif
          endif
          dMach1BC    = - Vpar0_b + direction / Btot * factor  * cs0_T * (Ti0_b+Te0_b)  &
                                  + direction / Btot * Hfact_b * cs0         
          dMach1BC_v  = - element_size_0
          dMach1BC_T  =           + direction / Btot * factor  * cs0_TT* T0_b   &
                                  + direction / Btot * Hfact_b * cs0_T
          dMach1BC_Ti =           + direction / Btot * factor  * cs0_TT* T0_b   &
                                  + direction / Btot * Hfact_b * cs0_T
          dMach1BC_Te =           + direction / Btot * factor  * cs0_TT* T0_b   &
                                  + direction / Btot * Hfact_b * cs0_T
          dMach1BC_Tb =           + direction / Btot * factor  * cs0_T * element_size_0

          ! --- SLOPE ROW OF THE SUPPLEMENT (bicubic). The tangential derivative of
          ! --- the same expression the value row imposes, on the branch frozen for
          ! --- this solve:
          ! ---   active:  d/db [m1_bfl * R^2*u_b/psi_b], with u_bb reconstructed from
          ! ---            both endpoints (u0_bb_r), psi_bb = ps0_bb, and the floor
          ! ---            m1_bfl treated as frozen geometry (bn is Dirichlet-static);
          ! ---   clipped: d/db [direction*2*cs/Btot], temperature part (field lagged,
          ! ---            as the whole block already lags magnetic geometry);
          ! ---   inactive: zero.
          ! --- RESIDUAL-ONLY for the u dependence: the active-branch derivative
          ! --- involves BOTH endpoints' u DOFs, and cross-node columns written into a
          ! --- row that several edges visit are exactly the stale-column trap fixed
          ! --- for the value row below. With the incidence floor the term is smooth
          ! --- and bounded, and the row keeps its strong Vpar_b diagonal, so lagging
          ! --- it is plain Picard on a bounded term. The clipped branch's exact
          ! --- temperature columns ARE carried (dMach1BC_Ti/Te/Tb below).
          ! ---
          ! --- SAFETY CLAMP at 2*m1_S: the imposed function b -> sup(b) lives in
          ! --- [0,S], so its mean slope across one edge cannot exceed the band per
          ! --- unit parameter; steeper pointwise structure is sub-element and not
          ! --- representable in the Hermite slope DOF anyway. Derived from the clip,
          ! --- not a new threshold; with the floor active it is rarely reached.
          m1_dslope = 0.d0
          m1_clamp  = 0.d0
          if ( n_order .eq. 3 .and. m1_dr .ne. 0.d0 ) then
            ! --- d/db of the FLOORED cancellation is f*D' + f'*D, not f*D' alone.
            ! --- The floor f = min(1,|bn|/s0) is frozen in TIME (psi is Dirichlet on
            ! --- the wall) but it is not constant along the boundary: bn varies, and
            ! --- its tangential derivative bn_b is already formed above. On the
            ! --- unfloored branch f = 1 identically, so f' = 0 there and the term is
            ! --- present only where the floor is actually doing something.
            m1_dfdb = 0.d0
            if ( m1_smin .gt. 0.d0 .and. abs(bn) .lt. m1_smin ) &
              m1_dfdb = sign(1.d0,bn) * bn_b / m1_smin
            m1_dDdb = ( 2.d0*BigR*R_b*U0_b + BigR**2*u0_bb_r                           &
                        - BigR**2*U0_b*ps0_bb/ps0_b ) / ps0_b
            m1_dSdb = 2.d0 * cs0_T * (Ti0_b+Te0_b) / Btot
            m1_raw  = m1_act * ( m1_bfl*m1_dDdb + m1_dfdb*m1_D )                       &
                    + m1_clw * direction * m1_dSdb
            ! --- The safety clamp is a branch like any other, so differentiate the
            ! --- branch that is actually taken: on the clamp the imposed slope is
            ! --- +-2*m1_S = +-4*cs/Btot, whose temperature derivative is +-4*cs_T/Btot
            ! --- and which does not depend on T_b at all. Adding the unclamped clip
            ! --- columns there would describe a term the residual is not using.
            if ( m1_raw .gt. 2.d0*m1_S ) then
              m1_clamp = +1.d0
            elseif ( m1_raw .lt. -2.d0*m1_S ) then
              m1_clamp = -1.d0
            endif
            m1_dslope = max( -2.d0*m1_S, min( 2.d0*m1_S, m1_raw ) )
            dMach1BC    = dMach1BC + m1_dr * m1_dslope
            if ( m1_clamp .ne. 0.d0 ) then
              dMach1BC_T  = dMach1BC_T  + m1_dr * m1_clamp * 4.d0 * cs0_T / Btot
              dMach1BC_Ti = dMach1BC_Ti + m1_dr * m1_clamp * 4.d0 * cs0_T / Btot
              dMach1BC_Te = dMach1BC_Te + m1_dr * m1_clamp * 4.d0 * cs0_T / Btot
            else
              dMach1BC_T  = dMach1BC_T  + m1_dr * m1_clw * direction * 2.d0 * cs0_TT * T0_b / Btot
              dMach1BC_Ti = dMach1BC_Ti + m1_dr * m1_clw * direction * 2.d0 * cs0_TT * T0_b / Btot
              dMach1BC_Te = dMach1BC_Te + m1_dr * m1_clw * direction * 2.d0 * cs0_TT * T0_b / Btot
              dMach1BC_Tb = dMach1BC_Tb + m1_dr * m1_clw * direction * 2.d0 * cs0_T * element_size_0 / Btot
            endif
          endif

          if (n_order .ge. 5) then
            ! --- Same convention as the value row: floored, gated off outside the
            ! --- active branch where the supplement has no u dependence.
            dMach1BC     = dMach1BC + m1_dr * m1_act * m1_bfl * BigR**2 * U0_bb/ps0_b
            dMach1BC_ubb = + m1_dr * m1_act * m1_bfl * BigR**2 * element_size_3/ps0_b
            d2Mach1BC    = - Vpar0_bb + direction / Btot * factor   * cs0_TT * (Ti0_b+Te0_b)**2   &
                                      + direction / Btot * factor   * cs0_T  * (Ti0_bb+Te0_bb)   !&
                                      !+ direction / Btot * Hfact_b  * cs0_T  * T0_b *2.0 !&
                                      !+ direction / Btot * Hfact_bb * cs0         
            d2Mach1BC_v  = - element_size_3
            d2Mach1BC_T  =            + direction / Btot * factor   * cs0_TTT * (Ti0_b+Te0_b)**2   &
                                      + direction / Btot * factor   * cs0_TT  * (Ti0_bb+Te0_bb)    &
                                      + direction / Btot * Hfact_b  * cs0_TT  * (Ti0_b+Te0_b) *2.0 &
                                      + direction / Btot * Hfact_bb * cs0_T 
            d2Mach1BC_Tb =            + direction / Btot * factor   * cs0_TT  * (Ti0_b+Te0_b) * 2.0 * element_size_0 &
                                      + direction / Btot * Hfact_b  * cs0_T                   * 2.0 * element_size_0 
            d2Mach1BC_Tbb=            + direction / Btot * factor   * cs0_T   * element_size_3 
          endif

          ! --- Apply Mach1
          ku = var_u
          kv = var_Vpar
          kT  = var_T
          kTi = var_Ti
          kTe = var_Te

          ! --- Impose Mach1 on node values
          call boundary_conditions_add_one_entry(             &
               index_node, kv, in, index_node, kv, in,        &
               - zbig * Mach1BC_v,                            &
               index_min, index_max, a_mat)

          if ( with_TiTe ) then
            call boundary_conditions_add_one_entry(             &
                 index_node, kv, in, index_node, kTi, in,       &
                 - zbig * Mach1BC_T,                            &
                 index_min, index_max, a_mat)
  
            call boundary_conditions_add_one_entry(             &
                 index_node, kv, in, index_node, kTe, in,       &
                 - zbig * Mach1BC_T,                            &
                 index_min, index_max, a_mat)
          else
            call boundary_conditions_add_one_entry(             &
                 index_node, kv, in, index_node, kT, in,        &
                 - zbig * Mach1BC_T,                            &
                 index_min, index_max, a_mat)
          endif

          ! --- CORNER DOUBLE-WRITE FIX. add_one_entry ASSIGNS, and at a geometric
          ! --- corner the two incident edges write this row's u column at DIFFERENT
          ! --- columns (node index(2) vs index(3)), so both used to survive while
          ! --- only the last RHS did - a spurious extra ExB column with its own
          ! --- 1/ps0_b amplification. Each visit now zeroes the OTHER tangential
          ! --- column first, so the last-visiting edge's relation is the only one
          ! --- left in the row. On a straight boundary the other column is written
          ! --- to the zero it already holds. Dormant when Mach1BC_u = 0 - which is
          ! --- exactly why it went unnoticed under dirichlet u.
          call boundary_conditions_add_one_entry(                            &
               index_node,  kv, in, node_list%node(inode)%index(5-iv_dir), ku, in, &
               0.d0,                                                         &
               index_min, index_max, a_mat)
          call boundary_conditions_add_one_entry(             &
               index_node,  kv, in, index_node2, ku, in,      &
               - zbig * Mach1BC_u,                            &
               index_min, index_max, a_mat)

          if (in .eq. 1) then
            call boundary_conditions_add_RHS(                        &
                   index_node, kv, in,index_min, index_max, RHS_loc, &
                   Zbig * Mach1BC,                                   &
                   a_mat%i_tor_min, a_mat%i_tor_max)
          else
            call boundary_conditions_add_RHS(                         &
                   index_node, kv, in, index_min, index_max, RHS_loc, &
                   0.d0,                                              &
                   a_mat%i_tor_min, a_mat%i_tor_max)
          endif
 
          ! --- Impose Mach1 on node derivatives
          call boundary_conditions_add_one_entry(               &
                 index_node2, kv, in, index_node2, kv, in,      &
                 - zbig * dMach1BC_v,                           &
                 index_min, index_max, a_mat)

          if ( with_TiTe ) then
            call boundary_conditions_add_one_entry(               &
                   index_node2, kv, in, index_node2, kTi, in,     &
                   - zbig * dMach1BC_Tb,                          &
                   index_min, index_max, a_mat)
 
            call boundary_conditions_add_one_entry(               &
                   index_node2, kv, in, index_node2, kTe, in,     &
                   - zbig * dMach1BC_Tb,                          &
                   index_min, index_max, a_mat)
 
            call boundary_conditions_add_one_entry(               &
                   index_node2, kv, in, index_node,  kTi, in,     &
                   - zbig * dMach1BC_Ti,                          & 
                   index_min, index_max, a_mat)
 
            call boundary_conditions_add_one_entry(               &
                   index_node2, kv, in, index_node,  kTe, in,     &
                   - zbig * dMach1BC_Te,                          & 
                   index_min, index_max, a_mat)
          else
            call boundary_conditions_add_one_entry(               &
                   index_node2, kv, in, index_node2, kT, in,      &
                   - zbig * dMach1BC_Tb,                          &
                   index_min, index_max, a_mat)
 
            call boundary_conditions_add_one_entry(               &
                   index_node2, kv, in, index_node,  kT, in,      &
                   - zbig * dMach1BC_T,                           & 
                   index_min, index_max, a_mat)
          endif

          if ( include_2nd_derivatives .and. (n_order .ge. 5) ) then
            call boundary_conditions_add_one_entry(               &
                   index_node2, kv, in, index_node3, ku, in,      &
                   - zbig * dMach1BC_ubb,                         &
                   index_min, index_max, a_mat)
          endif

          if (in .eq. 1) then
            call boundary_conditions_add_RHS(                          &
                   index_node2, kv, in, index_min, index_max, RHS_loc, &
                   Zbig * dMach1BC,                                    &
                   a_mat%i_tor_min, a_mat%i_tor_max)
          else
             call boundary_conditions_add_RHS(                         &
                   index_node2, kv, in, index_min, index_max, RHS_loc, &
                   0.d0,                                               &
                   a_mat%i_tor_min, a_mat%i_tor_max) 
          endif

          ! --- Impose Mach1 on node 2nd derivatives
!DOESNT WORK. I DONT KNOW WHY...
          if ( include_2nd_derivatives .and. (n_order .ge. 5) ) then
            call boundary_conditions_add_one_entry(               &
                   index_node3, kv, in, index_node3, kv,  in,     &
                   - zbig * d2Mach1BC_v,                          &
                   index_min, index_max, a_mat)
            call boundary_conditions_add_one_entry(               &
                   index_node3, kv, in, index_node , kTi, in,     &
                   - zbig * d2Mach1BC_T,                          &
                   index_min, index_max, a_mat)
            call boundary_conditions_add_one_entry(               &
                   index_node3, kv, in, index_node , kTe, in,     &
                   - zbig * d2Mach1BC_T,                          &
                   index_min, index_max, a_mat)
            call boundary_conditions_add_one_entry(               &
                   index_node3, kv, in, index_node2, kTi, in,     &
                   - zbig * d2Mach1BC_Tb,                         &
                   index_min, index_max, a_mat)
            call boundary_conditions_add_one_entry(               &
                   index_node3, kv, in, index_node2, kTe, in,     &
                   - zbig * d2Mach1BC_Tb,                         &
                   index_min, index_max, a_mat)
            call boundary_conditions_add_one_entry(               &
                   index_node3, kv, in, index_node3, kTi, in,     &
                   - zbig * d2Mach1BC_Tbb,                        &
                   index_min, index_max, a_mat)
            call boundary_conditions_add_one_entry(               &
                   index_node3, kv, in, index_node3, kTe, in,     &
                   - zbig * d2Mach1BC_Tbb,                        &
                   index_min, index_max, a_mat)
            if (in .eq. 1) then
              call boundary_conditions_add_RHS(                           &
                     index_node3, kv, in, index_min, index_max, RHS_loc,  &
                     Zbig * d2Mach1BC,                                    &
                     a_mat%i_tor_min, a_mat%i_tor_max)
            else
               call boundary_conditions_add_RHS(                         &
                     index_node3, kv, in, index_min, index_max, RHS_loc, &
                     0.d0,                                               &
                     a_mat%i_tor_min, a_mat%i_tor_max) 
            endif
          endif

          ! --- Fix derivatives in one direction
          k = var_Vpar
          do kk = 1,(n_order+1)/2
            do ll = 1,(n_order+1)/2
              if ( (iv_dir .eq. 2) .and. (ll .gt. 1) ) cycle ! do only pure s derivatives, not cross _st
              if ( (iv_dir .eq. 3) .and. (kk .gt. 1) ) cycle ! do only pure t derivatives, not cross _st
              if ( (iv_dir .eq. 2) .and. (kk .lt. 3) ) cycle ! do only node value, 1st and 2nd derivatives, fix the rest
              if ( (iv_dir .eq. 3) .and. (ll .lt. 3) ) cycle ! do only node value, 1st and 2nd derivatives, fix the rest
              index_tmp = node_indices(kk,ll)
              index_node = node_list%node(inode)%index(index_tmp)
              call boundary_conditions_add_one_entry(                 &
                     index_node, k, in, index_node, k, in,            &
                     zbig, index_min, index_max, a_mat)
            enddo
          enddo

        endif   !=== apply_cs
        
      enddo     !=== enddo loop n_tor
    
    enddo       !=== enddo directions (i_dir)

  enddo         !=== enddo vertex
 
enddo           !=== do elements

! --- FLOATING-U BOUNDARY DIAGNOSTIC: reduce and report -------------------------
! --- MAX/MIN reductions are safe under halo duplication. Ranks that do not own the
! --- axisymmetric harmonic contribute the sentinels and drop out of the extrema.
if ( floating_u_diag ) then
  do fd_t = 1, FD_NT
    fd_loc(1,fd_t) = fd_pe_max(fd_t)
    fd_loc(2,fd_t) = real(my_id,8)
  enddo
  call MPI_AllReduce(MPI_IN_PLACE, fd_res_max, FD_NT, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, err)
  call MPI_AllReduce(MPI_IN_PLACE, fd_vn_max,  FD_NT, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, err)
  call MPI_AllReduce(MPI_IN_PLACE, fd_rho_min, FD_NT, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, err)
  call MPI_AllReduce(MPI_IN_PLACE, fd_T_min,   FD_NT, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, err)
  call MPI_AllReduce(MPI_IN_PLACE, fd_vin_max, FD_NT, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, err)
  call MPI_AllReduce(MPI_IN_PLACE, fd_vout_max,FD_NT, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, err)
  call MPI_AllReduce(MPI_IN_PLACE, fd_m1cs_max,FD_NT, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, err)
  call MPI_AllReduce(MPI_IN_PLACE, fd_m1dr_max,FD_NT, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, err)
  call MPI_AllReduce(MPI_IN_PLACE, fd_loc,     FD_NT, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD, err)
  ! --- fd_pe_max is reduced through fd_loc (MAXLOC carries the owning rank so the
  ! --- location can be broadcast). Copy the reduced value back, otherwise rank 0
  ! --- reports its own local value - which is the -1 sentinel whenever rank 0 owns
  ! --- no floating boundary node.
  fd_pe_max(:) = fd_loc(1,:)
  call floating_u_norm(fu_a_n, fu_C_T, fu_C_V)
  fd_sq = fu_C_V * F0        ! = sqrt(mu0*rho0); v[m/s] = v[JOREK]/fd_sq
  do fd_t = 1, FD_NT
    if ( fd_loc(1,fd_t) .lt. 0.d0 ) cycle
    fd_owner = nint(fd_loc(2,fd_t))
    call MPI_Bcast(fd_pe_R(fd_t), 1, MPI_DOUBLE_PRECISION, fd_owner, MPI_COMM_WORLD, err)
    call MPI_Bcast(fd_pe_Z(fd_t), 1, MPI_DOUBLE_PRECISION, fd_owner, MPI_COMM_WORLD, err)
  enddo
  if ( my_id .eq. 0 ) then
    write(*,'(A)') ' [floating_u] type   |u-uf|[V]    |vE.n|[m/s]    Pe_core  Pe at (R,Z)          min rho    min T[eV]   vE.n out    vE.n IN     m1 cs      m1 drift   drift/cs'
    do fd_t = 1, FD_NT
      if ( fd_res_max(fd_t) .lt. 0.d0 ) cycle
      write(*,'(A,I3,4X,ES10.3,4X,ES10.3,4X,ES9.2,2X,A,F6.3,A,F7.3,A,4X,ES10.3,3X,ES10.3,3X,ES10.3,2X,ES10.3,2X,ES10.3,2X,ES10.3,2X,ES9.2)') &
        ' [floating_u] ', fd_t,                                    &
        fd_res_max(fd_t) / abs(fu_C_V),                            &
        fd_vn_max(fd_t)  / fd_sq,                                  &
        fd_pe_max(fd_t),                                           &
        '(', fd_pe_R(fd_t), ',', fd_pe_Z(fd_t), ')',               &
        fd_rho_min(fd_t),                                          &
        fd_T_min(fd_t) / ( MU_ZERO * central_density * 1.d20 * EL_CHG ),                &
        fd_vout_max(fd_t) / fd_sq, fd_vin_max(fd_t) / fd_sq,                          &
        fd_m1cs_max(fd_t), fd_m1dr_max(fd_t),                                          &
        ! --- A type never visited by a Mach row keeps its -1 sentinels; dividing the
        ! --- sentinel by tiny() printed -4.5e307. Report 0 for such types instead.
        merge( fd_m1dr_max(fd_t) / max(fd_m1cs_max(fd_t), tiny(1.d0)),                 &
               0.d0, fd_m1cs_max(fd_t) .gt. 0.d0 )
    enddo
  endif
endif

if (RMP_on) then
  if (allocated(psi_RMP_cos1))         call tr_deallocate(psi_RMP_cos1,"psi_RMP_cos1",CAT_UNKNOWN)
  if (allocated(dpsi_RMP_cos_dR1))     call tr_deallocate(dpsi_RMP_cos_dR1,"dpsi_RMP_cos_dR1",CAT_UNKNOWN)
  if (allocated(dpsi_RMP_cos_dZ1))     call tr_deallocate(dpsi_RMP_cos_dZ1,"dpsi_RMP_cos_dZ1",CAT_UNKNOWN)
  if (allocated(psi_RMP_sin1))         call tr_deallocate(psi_RMP_sin1,"psi_RMP_sin1",CAT_UNKNOWN)
  if (allocated(dpsi_RMP_sin_dR1))     call tr_deallocate(dpsi_RMP_sin_dR1,"dpsi_RMP_sin_dR1",CAT_UNKNOWN)
  if (allocated(dpsi_RMP_sin_dZ1))     call tr_deallocate(dpsi_RMP_sin_dZ1,"dpsi_RMP_sin_dZ1",CAT_UNKNOWN)
endif

return
end subroutine boundary_conditions 

end module mod_boundary_conditions
