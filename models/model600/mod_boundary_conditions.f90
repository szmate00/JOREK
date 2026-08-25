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
use corr_neg, only: corr_neg_temp, corr_neg_dens, dcorr_neg_dens_drho1
use mod_sheath_bc, only: sheath_get_lambda, sheath_current, sheath_V_wall_at, sheath_temp_floor
use mod_sheath_diag, only: sheath_psi0, sheath_store_psi0, sheath_diag_add_nodal
use phys_module, only: F0, GAMMA, freeboundary, RMP_on, psi_RMP_cos, dpsi_RMP_cos_dR, dpsi_RMP_cos_dZ, &
       psi_RMP_sin, dpsi_RMP_sin_dR, dpsi_RMP_sin_dZ, t_now, RMP_growth_rate, RMP_ramp_up_time,            &
       RMP_start_time, tstep, RMP_har_cos, RMP_har_sin, T_min,                                             &
       mach_one_bnd_integral, Vpar_smoothing, vpar_smoothing_coef, no_mach1_bc,                            &
       Number_RMP_harmonics, RMP_har_cos_spectrum,RMP_har_sin_spectrum, grid_to_wall, n_wall_blocks, keep_n0_const, &
       bcs, loop_voltage, central_density, central_mass,                                                    &
       sheath_Lambda, sheath_V_wall, sheath_u_exp_max, sheath_u_exp_min, sheath_u_relax, sheath_min_bn, &
       sheath_u_relax_time, sheath_wall_vel, sheath_u_align_psi, sheath_u_value_only, &
       sheath_zj_relax, sheath_ramp_time, t_start, sheath_zj_ratio_max
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
logical :: apply_sheath_u
logical :: apply_natural_u

! --- Sheath j-V boundary condition for the electric potential
real*8  :: u0, sh_a_n, sh_c_sat, sh_vw, sh_rho, sh_zj, sh_Ti, sh_Te, sh_T, sh_cs, sh_jsat
real*8  :: sh_lam, sh_dlam_dTi, sh_dlam_dTe
real*8  :: sh_g, sh_C, sh_dr_dzj
real*8  :: sh_relax
real*8  :: sh_wall_rhs
real*8  :: sh_perp_sign
real*8  :: sh_pl, sh_pn, sh_ul, sh_un, sh_nrm, sh_ca, sh_sa, sh_wgt
! --- nodal sheath current constraint on the zj row (bcs%sheath_zj)
logical :: apply_sheath_zj
real*8  :: szj_sh, szj_du, szj_drho, szj_dTi, szj_dTe, szj_sat, szj_x
real*8  :: szj_rel, szj_g, szj_coef(n_var), szj_R, szj_Rb, szj_ratio
integer :: k_szj
real*8  :: sh_wgt_bn
integer :: index_node_p
real*8  :: sh_ratio, sh_ratio_raw, sh_f_min, sh_f_max, sh_x, sh_xi, sh_dr, sh_u_targ
real*8  :: sh_R, sh_R_b, sh_R_bb
real*8  :: sh_coef(0:n_var)     ! index 0 is a scratch slot: var_T / var_Ti / var_Te are 0 in
real*8  :: sh_coef_d(0:n_var)   ! the model variants that lack them, so never index-checked
integer :: k_sh

!> Density floor (JOREK units) keeping 1/j_sat finite. Far below any physical SOL density.
!! corr_neg_dens is deliberately not used: its floor scales with rho_1, which is itself of the
!! order of the SOL density and would bias j_sat.
real*8, parameter :: sheath_rho_floor = 1.d-6

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

! --- Record psi at the first call, from the very node_list this routine indexes. Doing it here
! --- rather than in jorek2_main avoids any question of node redistribution or restart ordering
! --- changing the size or the contents between the two.
if ( (sheath_wall_vel .gt. 0.d0) .and. (.not. allocated(sheath_psi0)) ) &
  call sheath_store_psi0(node_list)

zbig        = 1.d12
zbig_backup = zbig

! --- calculate node_indices
call calculate_node_indices(node_indices)

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
            ! --- with the surface term of the current definition restored, zj is the honest
            ! --- Delta*psi and must not be frozen (see mod_boundary_matrix_open)
            if ( bcs(bnd_type)%natural%zj )   apply_current_BC = .false.
          endif
          !---------------------------------------------------------------------------------------------------                      

          
          !------------ Decide when to apply vpar=cs ---------------------------------------------------------                      
          apply_cs = .false.          
          if ( (.not. mach_one_bnd_integral) .and. bcs(bnd_type)%mach1 .and. with_vpar) then
            apply_cs = .true.
          endif

          !------------ Decide when to replace the Dirichlet BC on u by the sheath j-V BC --------
          apply_sheath_u = bcs(bnd_type)%sheath_u
          apply_sheath_zj = bcs(bnd_type)%sheath_zj
          ! --- The forward/natural route imposes the same sheath condition as a surface term,
          ! --- so it needs the same thin resistive wall. It does not use the nodal rows below.
          apply_natural_u = bcs(bnd_type)%natural%u
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
            if ( (k==var_u    ) .and.  apply_sheath_u                )  cycle  ! u is set by the sheath BC below
            if ( (k==var_zj   ) .and.  apply_sheath_zj               )  cycle  ! zj is set by the sheath current
            if ( (k==var_u    ) .and.  bcs(bnd_type)%natural%u        )  cycle  ! u is free: the charge-continuity
            if ( (k==var_w    ) .and.  bcs(bnd_type)%natural%w        )  cycle  ! surface terms take over (see
                                                                                ! mod_boundary_matrix_open)
            if ( (k==var_vpar ) .and.  apply_cs .and. (bnd_type/=3)  ) cycle  ! vpar=cs is a special case (this is done below)
                                                                              ! however bnd_type=3 needs both BCs for different directions

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

        ! --- Boundary conditions needing the local geometry (Mach1, sheath j-V BC for u)
        if ( apply_cs .or. apply_sheath_u .or. apply_sheath_zj .or.           &
             ( apply_natural_u .and. (sheath_wall_vel .gt. 0.d0) ) ) then

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

          ! --- the sheath BC does not need vpar, so this block can be reached without it
          ! --- (var_Vpar is 0 in that case, which would index out of bounds)
          if ( with_vpar ) then
            Vpar0   = node_list%node(inode)%values(1,1,var_vpar)
            Vpar0_b = node_list%node(inode)%values(1,iv_dir,var_Vpar) * element_size_0
          else
            Vpar0   = 0.d0
            Vpar0_b = 0.d0
          endif

          u0        = node_list%node(inode)%values(1,1,var_u)
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
            if ( with_vpar ) then
              Vpar0_bb = node_list%node(inode)%values(1,iv_dir+3,var_Vpar) * element_size_3
            else
              Vpar0_bb = 0.d0
            endif
            u0_bb     = node_list%node(inode)%values(1,iv_dir+3,var_u)    * element_size_3 
          endif

          Mach1: if ( apply_cs ) then

          ! --- Mach1 BC's and derivatives
          cs0      =   sqrt(gamma*(Ti0+Te0))
          cs0_T    =   0.5d0  * gamma    / cs0
          cs0_TT   = - 0.25d0 * gamma**2 / cs0**3 
          cs0_TTT  = 3.d0/8.d0* gamma**3 / cs0**5 

          Mach1BC     = - Vpar0   + direction / Btot * factor  * cs0               + factor / Btot * BigR**2 * U0_b/ps0_b 
          Mach1BC_v   = - 1.0
          Mach1BC_T   =           + direction / Btot * factor  * cs0_T 
          Mach1BC_u   =                                                            + factor / Btot * BigR**2 * element_size_0/ps0_b 
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


          if (n_order .ge. 5) then
            dMach1BC     = dMach1BC + factor / Btot * BigR**2 * U0_bb/ps0_b
            dMach1BC_ubb = + factor / Btot * BigR**2 * element_size_3/ps0_b
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

          endif Mach1

          !--------------------------------------------------------------------------------------
          ! --- Sheath j-V characteristic as boundary condition for the electric potential
          !
          ! --- J. Artola, "Sheath boundary conditions for the electric potential in JOREK".
          ! --- Stangeby 2.68 with the sheath factor:
          !
          ! ---     j = j_sat * f ,   f = 1 - exp(-X) ,   X = e*Phi/(k*Te) - Lambda
          !
          ! --- Phi = V_sheath_entrance - V_wall is referenced to the WALL, so f = 0 (zero current)
          ! --- sits at Phi = Lambda*Te/e, the floating potential, and Phi = 0 is electron
          ! --- saturation. In JOREK units, with Phi = F0*u,
          !
          ! ---     e*Phi/(k*Te) = ( a_n*u/2 - e*V_wall*mu0*n0 ) / Te
          ! ---     a_n   = 2*e*F0*sqrt(mu0*rho0)/(m_c*m_amu)
          ! ---     j_sat = c_sat*rho*v_par ,  c_sat = -e*F0*n0*sqrt(mu0/rho0) = -a_n/2
          ! ---     v_par = direction * c_s / |B|          (Mach 1 at the sheath entrance)
          !
          ! --- The characteristic is INVERTED for u rather than iterated on. With r = j/j_sat the
          ! --- exact solution of j = j_sat*(1-exp(-X)) is X = -ln(1-r), u = 2*(Te*(X+Lambda)+vw)/a_n.
          ! --- Clipping r - not the exponent - to the reachable range bounds X to
          ! --- [exp_min, exp_max] by construction, so u can never run away, and dr = 0 outside
          ! --- that range makes the Jacobian honest: where the requested current is beyond what
          ! --- the sheath can pass, u is simply the cap and has no rho or j dependence at all.
          ! --- Clipping the exponent instead keeps the row invertible but leaves every coefficient
          ! --- at its exponentially amplified value, which slaves u to the boundary noise in T,
          ! --- rho and j. Away from saturation this is Artola eq. 17 with f0 replaced by r.
          !
          ! --- The BC takes the u row. The current stays Dirichlet: psi, u, w and zj weak forms
          ! --- are all assembled without their surface terms (only rho, T, Ti, Te, rhon and vpar
          ! --- have natural BC support), so leaving any of them to its own equation at a boundary
          ! --- node makes that node absorb the missing surface integral.
          !--------------------------------------------------------------------------------------
          if ( apply_sheath_u ) then

            ! --- Effective under-relaxation of this boundary condition.
            ! --- u at the wall is slaved algebraically to Te, with no inertia and no dissipation
            ! --- anywhere in the Te -> u -> ExB flow -> Te loop, so that loop needs a response
            ! --- time of its own. Given as a per-step fraction (sheath_u_relax) its meaning
            ! --- changes every time tstep changes, which for a ramped tstep_n schedule means the
            ! --- damping silently weakens or strengthens at every block boundary. Given as a time
            ! --- (sheath_u_relax_time) the response is the same physical rate at any step size.
            sh_relax = sheath_u_relax
            if ( sheath_u_relax_time .gt. 0.d0 ) &
              sh_relax = min( 1.d0, tstep / sheath_u_relax_time )

            ! --- Normalisation constants (Artola eqs. 5 and 8; c_sat = -a_n/2 identically).
            ! --- NOTE the leading minus: the electrostatic potential is Phi = -F0*u in the code's
            ! --- variables. The JOREK reference paper (Hoelzl et al 2021 eq. 26) defines u = Phi/F0
            ! --- for v_pol = -R grad(u) x e_phi, while model600 implements v_pol = +R grad(u) x e_phi,
            ! --- so the code's u is minus the paper's. Flipping a_n propagates the correct sign to
            ! --- c_sat, to u_target and to every coefficient below, since all of them derive from it.
            sh_a_n   = - 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT) &
                     / (central_mass * ATOMIC_MASS_UNIT)
            sh_c_sat = - 0.5d0 * sh_a_n
            ! --- local wall potential, so an antisymmetric bias between the targets is felt
            sh_vw    = EL_CHG * sheath_V_wall_at(BigR) * MU_ZERO * central_density * 1.d20

            ! --- State at the boundary node (axisymmetric component)
            sh_rho   = max( node_list%node(inode)%values(1,1,var_rho), sheath_rho_floor )
            sh_zj    =      node_list%node(inode)%values(1,1,var_zj )
            sh_Ti    = corr_neg_temp(Ti0)
            sh_Te    = corr_neg_temp(Te0)
            sh_T     = sh_Ti + sh_Te

            ! --- Shared sheath factor, so this legacy nodal path and the charge-conserving
            ! --- surface-term form always use the same Lambda. NOTE: the dLambda/dT contribution
            ! --- to the coefficients below is neglected here (as are the corr_neg derivatives);
            ! --- this path is kept for A/B comparison only.
            call sheath_get_lambda(sh_Ti, sh_Te, sh_lam, sh_dlam_dTi, sh_dlam_dTe)

            ! --- Ion saturation current in JOREK units (Artola eqs. 5 and 16)
            sh_cs    = sqrt(GAMMA * sh_T)
            ! --- direction*factor is Artola's g(b_n): with vpar_smoothing and
            ! --- vpar_smoothing_coef(2) <= 0, direction is forced to +1 above and the sign of the
            ! --- field projection is carried by factor alone, so both are needed here.
            sh_g     = direction * factor
            sh_C     = sh_c_sat * sh_rho * sh_cs / Btot
            sh_jsat  = sh_C * sh_g

            ! --- Grazing incidence. Through a point where the field goes tangent to the wall g(b_n)
            ! --- passes through zero and changes sign, so the ratio r = zj/j_sat runs to +infinity
            ! --- on one side and -infinity on the other. Both ends clip, and the imposed potential
            ! --- jumps by (exp_max - exp_min)*Te between two neighbouring nodes - a + blob next to
            ! --- a - blob, anchored at the tangency point. Physically nothing is wrong there: no
            ! --- parallel flux reaches the wall, so the sheath carries no current and the surface
            ! --- simply floats. Weighting the ratio by g^2/(g^2 + g_min^2) says exactly that, and
            ! --- doing it as g/(g^2 + g_min^2) rather than as a separate factor removes the
            ! --- division by zero as well. sheath_min_bn = 0 reproduces the unweighted behaviour.
            sh_dr_dzj = sh_g / ( sh_C * ( sh_g**2 + sheath_min_bn**2 ) )

            ! --- How obliquely does the field actually meet the wall here? With grid_to_wall the
            ! --- same boundary type covers the divertor targets AND the main chamber, where the
            ! --- field is nearly tangential: no parallel flux reaches such a surface, so there is
            ! --- no sheath to impose, and worse, moving along that wall stays on roughly the same
            ! --- flux surface, so a potential varying along it is maximally NOT a flux function
            ! --- ([u,psi] = -psi_n*u_l with u_l imposed) and drags flux as hard as possible.
            ! --- This weight folds into the relaxation, so the row degenerates smoothly to du = 0,
            ! --- i.e. exactly the standard Dirichlet, where the field is tangential.
            sh_wgt_bn = sh_g**2 / ( sh_g**2 + sheath_min_bn**2 )
            sh_relax  = sh_relax * sh_wgt_bn

            ! --- Invert the characteristic, clipping the current ratio
            sh_f_min     = 1.d0 - exp(-sheath_u_exp_min)
            sh_f_max     = 1.d0 - exp(-sheath_u_exp_max)
            sh_ratio_raw = sh_zj * sh_dr_dzj
            sh_ratio     = min( max( sh_ratio_raw, sh_f_min ), sh_f_max )

            sh_x      = - log( 1.d0 - sh_ratio )              ! inside the clip range by construction
            sh_u_targ = 2.d0 * ( sh_Te * (sh_x + sh_lam) + sh_vw ) / sh_a_n
            sh_xi     = 2.d0 * sh_Te / ( sh_a_n * (1.d0 - sh_ratio) )

            sh_dr     = 1.d0
            if ( sh_ratio /= sh_ratio_raw ) sh_dr = 0.d0

            ! --- Row:  du + sum_k coef(k)*dx_k = u_target - u0 ,  coef = -d(u_target)/dx
            sh_coef             = 0.d0
            sh_coef(var_u  )    =   1.d0
            sh_coef(var_rho)    =   sh_dr * sh_xi * sh_ratio / sh_rho
            sh_coef(var_zj )    = - sh_dr * sh_xi * sh_dr_dzj
            if ( with_TiTe ) then
              sh_coef(var_Ti)   =   sh_dr * sh_xi * sh_ratio / (2.d0 * sh_T)
              sh_coef(var_Te)   =   sh_dr * sh_xi * sh_ratio / (2.d0 * sh_T) &
                                  - 2.d0 * (sh_x + sh_lam) / sh_a_n
            else
              sh_coef(var_T )   =   sh_dr * sh_xi * sh_ratio / (2.d0 * sh_T) &
                                  -        (sh_x + sh_lam) / sh_a_n
            endif
            sh_R                =   sh_u_targ - u0

            ! --- Same constraint differentiated along the boundary, for the derivative degrees of
            ! --- freedom. |B| and R are treated as constant along the boundary, as for Mach1.
            ! --- The dj/dl term is dropped: zj at a boundary node is a second derivative of psi,
            ! --- so its derivative degree of freedom is a third derivative, discontinuous across
            ! --- C1 cubic Bezier elements. It is grid noise rather than a gradient, and it
            ! --- dominates this row - and since du/dl is the ExB flow through the wall,
            ! --- v.n = R*du/dl, that noise would be fed straight back into the plasma. The value
            ! --- row keeps the full coupling, so the characteristic still holds exactly at the
            ! --- nodes; only the interpolation between them is affected.
            sh_coef_d           = sh_coef
            sh_coef_d(var_zj)   = 0.d0

            ! --- The derivative rows are written in nodal degree of freedom units, without the
            ! --- element_size scaling: it is common to all variables and cancels out, which keeps
            ! --- the diagonal entry at zbig.
            sh_R_b = 0.d0
            do k_sh = 1, n_var
              if ( sh_coef_d(k_sh) .ne. 0.d0 ) &
                sh_R_b = sh_R_b - sh_coef_d(k_sh) * node_list%node(inode)%values(1,iv_dir,k_sh)
            enddo

            index_node    = node_list%node(inode)%index(1)          ! position of value
            index_node2   = node_list%node(inode)%index(iv_dir)     ! position of first derivative

            ! --- Impose the characteristic on the node value of u. The diagonal is never relaxed,
            ! --- so the row stays well scaled and the zero-relaxation limit is du = 0 rather than
            ! --- a singular row.
            do k_sh = 1, n_var
              if ( sh_coef(k_sh) .eq. 0.d0 ) cycle
              if ( k_sh .eq. var_u ) then
                call boundary_conditions_add_one_entry(                      &
                       index_node, var_u, in, index_node, k_sh, in,          &
                       zbig * sh_coef(k_sh), index_min, index_max, a_mat)
              else
                call boundary_conditions_add_one_entry(                      &
                       index_node, var_u, in, index_node, k_sh, in,          &
                       zbig * sh_relax * sh_coef(k_sh), index_min, index_max, a_mat)
              endif
            enddo

            if (in .eq. 1) then
              call boundary_conditions_add_RHS(                              &
                     index_node, var_u, in, index_min, index_max, RHS_loc,   &
                     zbig * sh_relax * sh_R, a_mat%i_tor_min, a_mat%i_tor_max)
            else
              call boundary_conditions_add_RHS(                              &
                     index_node, var_u, in, index_min, index_max, RHS_loc,   &
                     0.d0, a_mat%i_tor_min, a_mat%i_tor_max)
            endif

            ! --- ... and on the first derivative of u along the boundary. Skipped when
            ! --- sheath_u_value_only: that row slaves du/dl to the tangential-derivative DOF of
            ! --- Te, and du/dl is v_E.n, the flux-dragging velocity, so the row hands grid-scale
            ! --- Te noise straight to the induction equation. Without it, u = u_target still holds
            ! --- at every node and the vorticity equation - which has dissipation - decides the
            ! --- interpolation between them.
            if ( .not. sheath_u_value_only ) then

            ! --- ... and on the first derivative of u along the boundary
            do k_sh = 1, n_var
              if ( sh_coef_d(k_sh) .eq. 0.d0 ) cycle
              if ( k_sh .eq. var_u ) then
                call boundary_conditions_add_one_entry(                      &
                       index_node2, var_u, in, index_node2, k_sh, in,        &
                       zbig * sh_coef_d(k_sh), index_min, index_max, a_mat)
              else
                call boundary_conditions_add_one_entry(                      &
                       index_node2, var_u, in, index_node2, k_sh, in,        &
                       zbig * sh_relax * sh_coef_d(k_sh), index_min, index_max, a_mat)
              endif
            enddo

            if (in .eq. 1) then
              call boundary_conditions_add_RHS(                              &
                     index_node2, var_u, in, index_min, index_max, RHS_loc,  &
                     zbig * sh_relax * sh_R_b, a_mat%i_tor_min, a_mat%i_tor_max)
            else
              call boundary_conditions_add_RHS(                              &
                     index_node2, var_u, in, index_min, index_max, RHS_loc,  &
                     0.d0, a_mat%i_tor_min, a_mat%i_tor_max)
            endif

            ! --- Second derivative along the boundary, higher order elements only. The variation
            ! --- of the coefficients along the boundary is neglected, so this row is consistent
            ! --- to first order only.
            if (n_order .ge. 5) then
              index_node3 = node_list%node(inode)%index(iv_dir+3)

              sh_R_bb = 0.d0
              do k_sh = 1, n_var
                if ( sh_coef_d(k_sh) .ne. 0.d0 ) &
                  sh_R_bb = sh_R_bb - sh_coef_d(k_sh) * node_list%node(inode)%values(1,iv_dir+3,k_sh)
              enddo

              do k_sh = 1, n_var
                if ( sh_coef_d(k_sh) .eq. 0.d0 ) cycle
                if ( k_sh .eq. var_u ) then
                  call boundary_conditions_add_one_entry(                    &
                         index_node3, var_u, in, index_node3, k_sh, in,      &
                         zbig * sh_coef_d(k_sh), index_min, index_max, a_mat)
                else
                  call boundary_conditions_add_one_entry(                    &
                         index_node3, var_u, in, index_node3, k_sh, in,      &
                         zbig * sh_relax * sh_coef_d(k_sh), index_min, index_max, a_mat)
                endif
              enddo

              if (in .eq. 1) then
                call boundary_conditions_add_RHS(                            &
                       index_node3, var_u, in, index_min, index_max, RHS_loc, &
                       zbig * sh_relax * sh_R_bb, a_mat%i_tor_min, a_mat%i_tor_max)
              else
                call boundary_conditions_add_RHS(                            &
                       index_node3, var_u, in, index_min, index_max, RHS_loc, &
                       0.d0, a_mat%i_tor_min, a_mat%i_tor_max)
              endif
            endif
            endif   !=== .not. sheath_u_value_only


            !------------------------------------------------------------------------------------
            ! --- Make the ExB flow at the wall run along flux surfaces:  [u,psi] = 0.
            !
            ! --- v_E.n = R du/dl is non-zero as soon as the sheath potential varies along the
            ! --- wall, and that is fine by itself: a flow along flux surfaces drags no poloidal
            ! --- flux. What drags flux is grad(u) not being parallel to grad(psi), and the piece
            ! --- that decides this - du/dn - is otherwise free, governed by the vorticity equation
            ! --- with nothing in it that knows about psi. The dragged flux then has to be balanced
            ! --- by a resistive current of order [u,psi]/eta, which at Spitzer resistivity forms a
            ! --- boundary layer of width eta/v_n, some 10 microns against a millimetre mesh. That
            ! --- is the field-aligned current filament, and it is why this boundary condition is
            ! --- hard in an MHD code and free in SOLPS or GBS, which hold B fixed.
            ! --- The small parallel resistivity of the SOL demands Phi ~ Phi(psi) anyway, so
            ! --- imposing it is physics rather than a numerical patch.
            !
            ! --- Written on raw degrees of freedom, so the element size scalings cancel exactly
            ! --- and no division by psi_l is needed - it vanishes where the field is tangent to
            ! --- the wall. There the constraint is meaningless (no flux is dragged either way) and
            ! --- the row blends smoothly back to leaving du/dn alone.
            !------------------------------------------------------------------------------------
            if ( sheath_u_align_psi .and. (sh_wgt_bn .gt. 0.25d0) ) then

              index_node_p = node_list%node(inode)%index(iv_perp_dir)

              sh_pl = node_list%node(inode)%values(1,iv_dir     ,var_psi)
              sh_pn = node_list%node(inode)%values(1,iv_perp_dir,var_psi)
              sh_ul = node_list%node(inode)%values(1,iv_dir     ,var_u  )
              sh_un = node_list%node(inode)%values(1,iv_perp_dir,var_u  )

              sh_nrm = sqrt( sh_pl**2 + sh_pn**2 )
              if ( sh_nrm .gt. 0.d0 ) then
                sh_ca = sh_pl / sh_nrm
                sh_sa = - sh_pn / sh_nrm
              else
                sh_ca = 1.d0
                sh_sa = 0.d0
              endif

              ! --- blend to du/dn = 0 where the field is tangent to the wall (sh_ca -> 0), so the
              ! --- row can never lose its diagonal
              sh_wgt = sh_ca**2 / ( sh_ca**2 + 2.5d-3 )      ! 2.5e-3 = (0.05)^2

              call boundary_conditions_add_one_entry(                             &
                     index_node_p, var_u, in, index_node_p, var_u, in,            &
                     zbig * ( sh_wgt * sh_ca + (1.d0 - sh_wgt) ),                  &
                     index_min, index_max, a_mat)
              call boundary_conditions_add_one_entry(                             &
                     index_node_p, var_u, in, index_node2, var_u, in,             &
                     zbig * sh_wgt * sh_sa, index_min, index_max, a_mat)

              if (in .eq. 1) then
                call boundary_conditions_add_RHS(                                 &
                       index_node_p, var_u, in, index_min, index_max, RHS_loc,    &
                       - zbig * sh_wgt * ( sh_ca * sh_un + sh_sa * sh_ul ),         &
                       a_mat%i_tor_min, a_mat%i_tor_max)
              else
                call boundary_conditions_add_RHS(                                 &
                       index_node_p, var_u, in, index_min, index_max, RHS_loc,    &
                       0.d0, a_mat%i_tor_min, a_mat%i_tor_max)
              endif

            endif

          endif   !=== apply_sheath_u

          !--------------------------------------------------------------------------------------
          ! --- SHEATH CURRENT ON THE zj ROW (bcs%sheath_zj)
          !
          ! --- The sheath sets the CURRENT that crosses the wall, so that is a statement about
          ! --- zj, and this imposes it directly:
          ! ---     zj = zj_sheath(u, rho, Ti, Te) = zj_sat * ( 1 - exp(-X) )
          ! --- replacing the Dirichlet that otherwise freezes zj at its initial value. The
          ! --- characteristic is used in the FORWARD direction, current as a function of
          ! --- potential, which is single valued and whose derivative tends smoothly to zero at
          ! --- ion saturation - unlike the inverted form used by the sheath_u path above, which
          ! --- is singular exactly at divertor conditions.
          !
          ! --- Why not a surface term on the u equation: that equation is assembled in STRONG
          ! --- form (mod_elt_matrix_fft, the +v*(ps0_s*zj0_t - ps0_t*zj0_s) term has its
          ! --- derivatives on zj, not on the test function), so integration by parts was never
          ! --- performed and there is no boundary flux for a natural BC to replace. Adding one
          ! --- injects a spurious source at the wall instead of redirecting a flux.
          !
          ! --- u stays free here (dirichlet%u = .false.): charge continuity determines it, and it
          ! --- settles where the current the plasma delivers matches what the sheath can pass.
          ! --- Keep a Dirichlet on u on at least one other boundary type as a gauge.
          !--------------------------------------------------------------------------------------
          if ( apply_sheath_zj ) then

            szj_g = direction * factor            ! Chodura-Riemann g(b_n), signed, as for Mach1

            ! --- corr_neg_dens, NOT a bare max(): the Gauss-point path uses r0_corr, and the
            ! --- node/Gauss residual comparison is only meaningful if both evaluate the SAME
            ! --- characteristic from the same corrected state. Numerically almost identical here
            ! --- (target densities sit well above rho_min_neg) but the inconsistency made the two
            ! --- diagnostics measure subtly different things.
            call sheath_current( u0,                                                  &
                                 corr_neg_dens(node_list%node(inode)%values(1,1,var_rho)), &
                                 sheath_temp_floor(Ti0), sheath_temp_floor(Te0),      &
                                 szj_g, sign(1.d0, bn), Btot,                         &
                                 szj_sh, szj_du, szj_drho, szj_dTi, szj_dTe,          &
                                 szj_sat, szj_x, sheath_V_wall_at(BigR) )

            ! --- Effective strength. The diagonal below is never relaxed, so this factor going to
            ! --- zero leaves the row as d(zj) = 0, i.e. exactly the Dirichlet freeze this
            ! --- replaces. That makes both ends of the continuation well posed.
            szj_rel = sheath_zj_relax
            if ( sheath_ramp_time .gt. 0.d0 ) &
              szj_rel = szj_rel * max(0.d0, min(1.d0, (t_now - t_start) / sheath_ramp_time))

            ! --- Obliqueness. Gate on bn = B.n/|B| itself, NOT on szj_g: the Chodura factor is
            ! --- 1.0 everywhere except on the few edges that straddle a sign change in b_n (see
            ! --- the vpar_smoothing branch above), so szj_g = direction*factor is +-1 over almost
            ! --- the whole boundary and a gate built from it would be identically 1.
            ! --- Two things need suppressing where the field grazes the wall: the sheath carries
            ! --- no current there, and direction = sign(b_n) flips discontinuously through the
            ! --- tangency point, which would flip the imposed zj with it. bn^2/(bn^2+min_bn^2)
            ! --- is smooth, never negative, and hands those nodes back to the frozen value.
            if ( sheath_min_bn .gt. 0.d0 ) &
              szj_rel = szj_rel * bn**2 / ( bn**2 + sheath_min_bn**2 )

            ! --- Validity gate. Where the plasma demands |zj/zj_sat| far above 1 the sheath
            ! --- cannot pass that current at ANY potential: on the ion side f -> 1, so the
            ! --- residual settles at a constant and u is driven without bound. Worse, if j_sat
            ! --- collapses locally (a cooling or rarefying spot on the target) the demand runs
            ! --- away and u chases it - measured 5 -> 265 while u reached 500*Te. Outside that
            ! --- range the characteristic carries no useful information, so fade the constraint
            ! --- out and let the node keep its frozen zj, which is the baseline and is stable.
            ! --- The weight is smooth, monotone, and 1 to within 6% for ratio < 0.5*max.
            szj_ratio = 0.d0
            if ( abs(szj_sat) .gt. 1.d-30 ) &
              szj_ratio = abs( node_list%node(inode)%values(1,1,var_zj) / szj_sat )
            if ( sheath_zj_ratio_max .gt. 0.d0 ) &
              szj_rel = szj_rel / ( 1.d0 + (szj_ratio/sheath_zj_ratio_max)**4 )

            ! --- Report the strength the row is ACTUALLY assembled with. Without this the only
            ! --- gate visible in the SHEATH output is the obliqueness one, and a small szj_rel
            ! --- silently degrades the constraint to the Dirichlet freeze it is meant to replace.
            call sheath_diag_add_nodal(szj_rel, szj_ratio, BigR,                &
                                       node_list%node(inode)%x(1,1,2),          &
                                       szj_sh - node_list%node(inode)%values(1,1,var_zj), &
                                       szj_sat)

            ! --- Row:  d(zj) - sum_k (d zj_sh / d x_k) d x_k = zj_sh - zj0
            szj_coef            = 0.d0
            szj_coef(var_zj )   =   1.d0
            szj_coef(var_u  )   = - szj_du
            ! --- chain the density correction through, now that sheath_current is given the
            ! --- corrected density: szj_drho is d/d(rho_corr), the row needs d/d(rho)
            szj_drho = szj_drho * dcorr_neg_dens_drho1(node_list%node(inode)%values(1,1,var_rho))
            szj_coef(var_rho)   = - szj_drho
            if ( with_TiTe ) then
              szj_coef(var_Ti)  = - szj_dTi
              szj_coef(var_Te)  = - szj_dTe
            else
              szj_coef(var_T )  = - 0.5d0 * ( szj_dTi + szj_dTe )
            endif
            szj_R = szj_sh - node_list%node(inode)%values(1,1,var_zj)

            index_node  = node_list%node(inode)%index(1)
            index_node2 = node_list%node(inode)%index(iv_dir)

            do k_szj = 1, n_var
              if ( szj_coef(k_szj) .eq. 0.d0 ) cycle
              if ( k_szj .eq. var_zj ) then
                call boundary_conditions_add_one_entry(                        &
                       index_node, var_zj, in, index_node, k_szj, in,          &
                       zbig * szj_coef(k_szj), index_min, index_max, a_mat)
              else
                call boundary_conditions_add_one_entry(                        &
                       index_node, var_zj, in, index_node, k_szj, in,          &
                       zbig * szj_rel * szj_coef(k_szj), index_min, index_max, a_mat)
              endif
            enddo

            if (in .eq. 1) then
              call boundary_conditions_add_RHS(                                &
                     index_node, var_zj, in, index_min, index_max, RHS_loc,    &
                     zbig * szj_rel * szj_R, a_mat%i_tor_min, a_mat%i_tor_max)
            else
              call boundary_conditions_add_RHS(                                &
                     index_node, var_zj, in, index_min, index_max, RHS_loc,    &
                     0.d0, a_mat%i_tor_min, a_mat%i_tor_max)
            endif

            ! --- The same constraint differentiated along the boundary, for the tangential
            ! --- derivative degree of freedom. |B| and R are treated as constant along the
            ! --- boundary, as for Mach1 and for the sheath_u path.
            szj_Rb = 0.d0
            do k_szj = 1, n_var
              if ( szj_coef(k_szj) .ne. 0.d0 ) &
                szj_Rb = szj_Rb - szj_coef(k_szj) * node_list%node(inode)%values(1,iv_dir,k_szj)
            enddo

            do k_szj = 1, n_var
              if ( szj_coef(k_szj) .eq. 0.d0 ) cycle
              if ( k_szj .eq. var_zj ) then
                call boundary_conditions_add_one_entry(                        &
                       index_node2, var_zj, in, index_node2, k_szj, in,        &
                       zbig * szj_coef(k_szj), index_min, index_max, a_mat)
              else
                call boundary_conditions_add_one_entry(                        &
                       index_node2, var_zj, in, index_node2, k_szj, in,        &
                       zbig * szj_rel * szj_coef(k_szj), index_min, index_max, a_mat)
              endif
            enddo

            if (in .eq. 1) then
              call boundary_conditions_add_RHS(                                &
                     index_node2, var_zj, in, index_min, index_max, RHS_loc,   &
                     zbig * szj_rel * szj_Rb, a_mat%i_tor_min, a_mat%i_tor_max)
            else
              call boundary_conditions_add_RHS(                                &
                     index_node2, var_zj, in, index_min, index_max, RHS_loc,   &
                     0.d0, a_mat%i_tor_min, a_mat%i_tor_max)
            endif

          endif   !=== apply_sheath_zj

          !------------------------------------------------------------------------------------
          ! --- Thin resistive wall for the poloidal flux.
          ! --- Applies to BOTH sheath routes (bcs%sheath_u and bcs%natural%u): the flux dragging
          ! --- below is caused by u varying along the wall, not by how that u was imposed.
          !
          ! --- Once u varies along the wall there is an ExB flow through it, v_E.n = R du/dl,
          ! --- and that flow drags poloidal flux, through the R[u,psi] term of the induction
          ! --- equation. Parallel flow is harmless (B.grad psi = 0); only the perpendicular flow
          ! --- drags. Against a Dirichlet psi the dragged flux has nowhere to go and piles up in
          ! --- a resistive layer of width eta/v_n - of order 10 microns at Spitzer resistivity
          ! --- against a millimetre mesh. Unresolvable, and it appears as the grid-scale current
          ! --- filament. This is the whole reason the boundary condition is hard here and free
          ! --- in SOLPS, SOLEDGE or GBS, which hold B fixed: no flux to drag, no layer.
          !
          ! --- Let the wall pass flux at a finite speed instead of freezing it:
          ! ---     dpsi/dt = - sheath_wall_vel * ( dpsi/dn - dpsi/dn at t_start )
          ! --- The DEVIATION from the initial state is what drives it; using the raw gradient
          ! --- would bleed flux even in a quiet plasma. The sign is a relaxation: one sided
          ! --- diffusion, dpsi/dt = eta*grad^2(psi) ~ -(2*eta/h)*dpsi/dn, so an excess outward
          ! --- gradient moves psi to reduce it. A useful magnitude is v_n itself - the wall has
          ! --- to pass flux about as fast as the flow delivers it.
          !------------------------------------------------------------------------------------
          if ( (sheath_wall_vel .gt. 0.d0)                                                 &
               .and. (apply_sheath_u .or. apply_natural_u .or. apply_sheath_zj)            &
               .and. (.not. is_freebound(in,var_psi))                                      &
               .and. allocated(sheath_psi0) ) then
            if ( inode .gt. size(sheath_psi0,2) ) then
              if (my_id .eq. 0) write(*,*) 'WARNING: sheath_psi0 too small, wall relaxation off'
            else

            ! --- Mach1 above may have left index_node pointing at a derivative DOF
            index_node   = node_list%node(inode)%index(1)
            index_node_p = node_list%node(inode)%index(iv_perp_dir)

            ! --- Orientation. values(1,iv_perp_dir,...) is a derivative with respect to a LOGICAL
            ! --- coordinate whose sense depends on which vertex the boundary edge hangs off; the
            ! --- Mach1 block above applies exactly this correction to element_size_s/t. Without
            ! --- it the relaxation dpsi/dt = -v_w*(dpsi/dn - dpsi/dn|_0) is dissipative on some
            ! --- boundary elements and ANTI-dissipative on the rest, growing a mesh-scale
            ! --- alternating structure in psi - visible first in zj = Delta*psi, only later in u.
            sh_perp_sign = 1.d0
            if ( (iv_perp_dir .eq. 2) .and. ((iv .eq. 2) .or. (iv .eq. 3)) ) sh_perp_sign = -1.d0
            if ( (iv_perp_dir .eq. 3) .and. ((iv .eq. 3) .or. (iv .eq. 4)) ) sh_perp_sign = -1.d0

            sh_wall_rhs = - tstep * sheath_wall_vel * sh_perp_sign                             &
                          * ( node_list%node(inode)%values(1,iv_perp_dir,var_psi)              &
                            - sheath_psi0(iv_perp_dir,inode) )
            ! --- loop_voltage drives the same row and would otherwise be overwritten here
            if ( loop_voltage .ne. 0.d0 ) sh_wall_rhs = sh_wall_rhs                            &
              + loop_voltage * sqrt(MU_ZERO*central_density*central_mass*ATOMIC_MASS_UNIT*1.d20) * tstep

            call boundary_conditions_add_one_entry(                         &
                   index_node, var_psi, in, index_node, var_psi, in,        &
                   zbig, index_min, index_max, a_mat)
            call boundary_conditions_add_one_entry(                         &
                   index_node, var_psi, in, index_node_p, var_psi, in,      &
                   zbig * tstep * sheath_wall_vel * sh_perp_sign, index_min, index_max, a_mat)
            if (in .eq. 1) then
              call boundary_conditions_add_RHS(                             &
                     index_node, var_psi, in, index_min, index_max, RHS_loc,&
                     zbig * sh_wall_rhs, a_mat%i_tor_min, a_mat%i_tor_max)
            else
              call boundary_conditions_add_RHS(                             &
                     index_node, var_psi, in, index_min, index_max, RHS_loc,&
                     0.d0, a_mat%i_tor_min, a_mat%i_tor_max)
            endif
            endif

          endif


        endif   !=== apply_cs .or. sheath_u .or. sheath_zj .or. natural_u+wall
        
      enddo     !=== enddo loop n_tor
    
    enddo       !=== enddo directions (i_dir)

  enddo         !=== enddo vertex
 
enddo           !=== do elements

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
