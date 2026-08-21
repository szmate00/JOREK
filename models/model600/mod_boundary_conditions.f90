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

use constants, only : PI, MU_ZERO, ATOMIC_MASS_UNIT, EL_CHG, MASS_ELECTRON
use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
use data_structure
use vacuum, ONLY: is_freebound
use phys_module, only: F0, GAMMA, freeboundary, RMP_on, psi_RMP_cos, dpsi_RMP_cos_dR, dpsi_RMP_cos_dZ, &
       psi_RMP_sin, dpsi_RMP_sin_dR, dpsi_RMP_sin_dZ, t_now, RMP_growth_rate, RMP_ramp_up_time,            &
       RMP_start_time, tstep, RMP_har_cos, RMP_har_sin, T_min,                                             &
       mach_one_bnd_integral, Vpar_smoothing, vpar_smoothing_coef, no_mach1_bc,                            &
       Number_RMP_harmonics, RMP_har_cos_spectrum,RMP_har_sin_spectrum, grid_to_wall, n_wall_blocks, keep_n0_const, &
       bcs, loop_voltage, central_density, central_mass,                                        &
       floating_Lambda, floating_Lambda_local, floating_V_wall, floating_u_relax,              &
       floating_ramp_time, t_start, floating_u_value_only, floating_min_bn,            &
       floating_start_time, floating_gauge_removal, floating_amp_ramp, mach1_psib_floor, &
       mach1_exb_term
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
! --- floating-potential boundary condition on u (bcs%floating_u)
logical :: apply_floating_u
real*8  :: flt_Ti, flt_Te, flt_T, flt_lam, flt_dlam_dTi, flt_dlam_dTe
real*8  :: flt_a_n, flt_vw, flt_u, flt_du_dTi, flt_du_dTe, flt_rel
real*8  :: flt_coef(n_var), flt_R, flt_Rb, flt_lam0
integer :: k_flt, flt_row
logical :: flt_clip_Ti, flt_clip_Te
real*8  :: flt_amp, flt_x, flt_t0
real*8  :: psib_inv
integer, save     :: flt_tr_unit = -1
logical, save     :: flt_tr_open = .false., flt_tr_ok = .false.
character(len=64) :: flt_tr_file
real*8            :: flt_tr_Te, flt_tr_Td, flt_tr_tgt
! --- Gauge reference: one constant subtracted from the floating target over the whole connected
! --- wall. Accumulated during a sweep and used in the NEXT one, so no extra pass over the mesh is
! --- needed. The lag is harmless - the wall temperature moves slowly, and more fundamentally ANY
! --- spatially uniform value is a valid gauge, so this constant does not need to be accurate. It
! --- only has to be the same everywhere on the wall and roughly centred, which is also why an
! --- unweighted nodal mean is enough and arclength weighting buys nothing here.
real*8, save :: flt_gauge = 0.d0
real*8, save :: flt_gsum  = 0.d0, flt_gcnt = 0.d0
real*8       :: flt_gs_glob, flt_gc_glob
integer :: flt_trt, flt_trt2, flt_bnd2, flt_njump
logical, save :: flt_scanned = .false.
! --- corner-consistency diagnostic for the mach1 direction/factor (see below)
real*8,  allocatable, save :: bcdiag_dir(:), bcdiag_fac(:)
integer, allocatable, save :: bcdiag_vis(:)
integer, save              :: bcdiag_unit = -1
logical, save              :: bcdiag_open = .false.   ! newunit= returns a NEGATIVE unit, so the
                                                      ! sign of bcdiag_unit says nothing about
                                                      ! whether the file is open. Track it here.
logical, save              :: bcdiag_ok   = .false.   ! open succeeded (bcdiag_ios is a local and
                                                      ! is undefined on visits that skip the open)
integer, save              :: bcdiag_calls = 0, bcdiag_hits = 0
integer                    :: bcdiag_ios
character(len=64)          :: bcdiag_file
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
          endif
          !---------------------------------------------------------------------------------------------------                      

          
          !------------ Decide when to hold u at the floating potential --------------------------------------
          apply_floating_u = bcs(bnd_type)%floating_u

          !------------ Decide when to apply vpar=cs ---------------------------------------------------------                      
          apply_cs = .false.          
          if ( (.not. mach_one_bnd_integral) .and. bcs(bnd_type)%mach1 .and. with_vpar) then
            apply_cs = .true.
          endif
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
            if ( (k==var_u    ) .and.  apply_floating_u             )  cycle  ! u is set to Phi_float below
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
          !--------------------------------------------------------------------------------------------
          ! --- FLOATING-POTENTIAL BOUNDARY CONDITION ON u  (bcs%floating_u)
          !
          ! --- Hold the plasma potential at the wall at the value a surface takes up when it draws no
          ! --- net current - the floating potential:
          ! ---     Phi_float = Lambda*Te/e + V_wall ,   Lambda = Lambda_0 - ln sqrt(gamma*(1+Ti/Te))
          ! --- The standard Dirichlet pins u = 0, i.e. Phi = 0, which is not a neutral choice: it is
          ! --- the electron-saturation point, and it suppresses the SOL electric field entirely.
          !
          ! --- Because Te varies along the target, so does Phi, giving E ~ Lambda*grad(Te)/e along the
          ! --- wall and the associated ExB drift. That drift - notably through the private flux region
          ! --- - is the leading-order driver of the in-out density asymmetry, so this alone is enough
          ! --- to let a SOL potential structure develop. What it does NOT capture is the deviation of
          ! --- Phi from floating that a net current would produce (thermoelectric currents); that
          ! --- needs the full j-V characteristic and a free potential.
          !
          ! --- NOTE the sign convention: the electrostatic potential is Phi = -F0*u in the code's
          ! --- variables. model600 implements v_pol = +R grad(u) x e_phi, whereas the JOREK reference
          ! --- paper (Hoelzl et al 2021 eq. 26) defines u = Phi/F0 with v_pol = -R grad(u) x e_phi, so
          ! --- the code's u is minus the paper's. That is where the minus in a_n comes from, and every
          ! --- coefficient below follows from it.
          !
          ! --- The row is the Dirichlet it replaces, with a state-dependent target instead of zero:
          ! ---     du - (dPhi_float/dTi) dTi - (dPhi_float/dTe) dTe = u_float - u0
          ! --- The diagonal is never relaxed, so floating_u_relax -> 0 reproduces the plain Dirichlet
          ! --- exactly. Together with floating_ramp_time that gives a continuation whose two ends are
          ! --- both well posed - the baseline at one end, the floating potential at the other.
          !--------------------------------------------------------------------------------------------
          ! --- CONSISTENCY OF direction/factor AT CORNER NODES
          ! ---
          ! --- A boundary node at a corner is visited once from each adjacent boundary element,
          ! --- and each visit computes its own outward normal. direction = sign(B.n) and the
          ! --- Chodura factor follow from that normal, so the two visits can disagree - most
          ! --- easily where B.n passes through zero, which is precisely what happens at a corner
          ! --- between a target (B near-normal) and a side wall (B grazing).
          ! ---
          ! --- boundary_conditions_add_one_entry OVERWRITES rather than accumulates, so when they
          ! --- disagree the value imposed is whichever element the loop reaches last: v_par = +c_s
          ! --- or -c_s decided by element ordering, which is not physics. vpar_smoothing only
          ! --- covers this when the edge straddles a sign change (bn_1*bn_2 < 0); if bn is small
          ! --- but keeps its sign across the corner, factor = 1 and the full +-c_s is imposed at
          ! --- grazing incidence.
          ! ---
          ! --- Report it once per run, with coordinates, so a corner artefact is distinguishable
          ! --- from a physical structure at the same place.
          if ( (.not. flt_scanned) .and. bcs(bnd_type)%mach1 ) then
            if ( .not. allocated(bcdiag_vis) ) then
              allocate( bcdiag_vis(node_list%n_nodes) ) ; bcdiag_vis = 0
              allocate( bcdiag_dir(node_list%n_nodes) ) ; bcdiag_dir = 0.d0
              allocate( bcdiag_fac(node_list%n_nodes) ) ; bcdiag_fac = 0.d0
            endif
            ! --- The comparison below needs both of a corner's elements on THIS rank, and the
            ! --- element loop is distributed, so a corner split across ranks is missed. Dump the
            ! --- raw geometry too: adjacent nodes with opposite direction are visible in the file
            ! --- regardless of how the mesh was partitioned.
            bcdiag_hits = bcdiag_hits + 1
            if ( .not. bcdiag_open ) then
              bcdiag_open = .true.
              write(bcdiag_file,'(A,I0,A)') 'bc_geom_diag_', my_id, '.dat'
              open(newunit=bcdiag_unit, file=trim(bcdiag_file), status='replace',                 &
                   action='write', iostat=bcdiag_ios)
              if ( bcdiag_ios .ne. 0 ) then
                write(*,'(A,A,A,I0)') ' WARNING [boundary_conditions]: could not open ',          &
                      trim(bcdiag_file), ', iostat = ', bcdiag_ios
                bcdiag_unit = -1
              else
                bcdiag_ok = .true.
                write(bcdiag_unit,'(A)') '#            R               Z  bnd_type    iv_dir'//   &
                                         '              bn       direction          factor'//   &
                                         '          ps0_b'
              endif
            endif
            if ( bcdiag_ok ) then
              ! --- ps0_b is dumped so mach1_psib_floor can be chosen from the actual distribution
              ! --- rather than guessed: it is the quantity the Mach-1 ExB term divides by.
              write(bcdiag_unit,'(2f16.8,2i10,4es16.6)')                                          &
                node_list%node(inode)%x(1,1,1), node_list%node(inode)%x(1,1,2),                   &
                bnd_type, iv_dir, bn, direction, factor, ps0_b
            endif

            if ( bcdiag_vis(inode) .eq. 0 ) then
              bcdiag_vis(inode) = 1
              bcdiag_dir(inode) = direction
              bcdiag_fac(inode) = factor
            else
              if ( (bcdiag_dir(inode)*direction .lt. 0.d0) .or.                                    &
                   (abs(bcdiag_fac(inode)-factor) .gt. 1.d-3) ) then
                write(*,'(A,2f9.4,A,2f7.2,A,2es11.3)')                                             &
                  ' WARNING [boundary_conditions]: mach1 corner disagreement at R,Z = ',           &
                  node_list%node(inode)%x(1,1,1), node_list%node(inode)%x(1,1,2),                  &
                  '  direction ', bcdiag_dir(inode), direction, '  factor ', bcdiag_fac(inode), factor
              endif
            endif
          endif

          ! --- CONSISTENCY OF THE u BC ACROSS BOUNDARY-TYPE JUNCTIONS
          ! ---
          ! --- bnd_type is a property of the NODE, so two nodes at opposite ends of the same
          ! --- boundary edge may carry different types and therefore different conditions on u.
          ! --- If one is held at the floating potential and its neighbour is pinned by the plain
          ! --- Dirichlet, u jumps by the full Phi_float across a single element - tens of volts at
          ! --- ordinary target temperatures. The resulting grad(u) is a spurious ExB velocity
          ! --- localised at that junction, with no physical content whatsoever: it is an artefact
          ! --- of the two types having been configured differently.
          ! ---
          ! --- This is silent otherwise, and it is worst exactly where the cells are smallest,
          ! --- so it is reported once per run, with coordinates.
          if ( .not. flt_scanned ) then
            flt_bnd2 = node_list%node(inode2)%boundary
            flt_trt  = 0
            flt_trt2 = 0
            if ( bnd_type .gt. 0 ) then
              if ( bcs(bnd_type)%floating_u ) then
                flt_trt = 2
              elseif ( bcs(bnd_type)%dirichlet%u ) then
                flt_trt = 1
              endif
            endif
            if ( flt_bnd2 .gt. 0 ) then
              if ( bcs(flt_bnd2)%floating_u ) then
                flt_trt2 = 2
              elseif ( bcs(flt_bnd2)%dirichlet%u ) then
                flt_trt2 = 1
              endif
            endif
            if ( (flt_trt .ne. flt_trt2) .and. (max(flt_trt,flt_trt2) .eq. 2) ) then
              write(*,'(A,I3,A,I3,A,2f9.4)')                                                       &
                ' WARNING [boundary_conditions]: u BC changes across a boundary edge, type ',      &
                bnd_type, ' -> ', flt_bnd2, '  at R,Z = ',                                         &
                node_list%node(inode)%x(1,1,1), node_list%node(inode)%x(1,1,2)
              write(*,'(A)')                                                                       &
                '          u jumps by the full floating potential over one element here. Set'//    &
                ' floating_u on'
              write(*,'(A)')                                                                       &
                '          BOTH types (normally every type where field lines strike the wall,'//   &
                ' i.e. every'
              write(*,'(A)')                                                                       &
                '          type carrying mach1) or the jump drives a spurious ExB flow at this'//  &
                ' junction.'
            endif
          endif

          if ( apply_floating_u ) then

            ! --- T_min is a HARD clip, so where it bites u_float no longer depends on the
            ! --- temperature and neither may the Jacobian. Carrying the unclipped derivative there
            ! --- would linearise a constant, which is exactly the sort of inconsistency that shows
            ! --- up as a slow boundary instability rather than an obvious error.
            ! --- The clip is applied to Ti and Te SEPARATELY. A single flag covering both zeroed
            ! --- du/dTe whenever Ti clipped, which discards a derivative that is perfectly valid:
            ! --- Te is what sets Phi_float, and Ti only enters through Lambda.
            flt_clip_Ti = .false.
            flt_clip_Te = .false.
            if ( with_TiTe ) then
              if ( node_list%node(inode)%values(1,1,var_Ti) .lt. T_min ) flt_clip_Ti = .true.
              if ( node_list%node(inode)%values(1,1,var_Te) .lt. T_min ) flt_clip_Te = .true.
            else
              if ( node_list%node(inode)%values(1,1,var_T ) .lt. T_min ) then
                flt_clip_Ti = .true.
                flt_clip_Te = .true.
              endif
            endif

            if ( with_TiTe ) then
              flt_Ti = max(node_list%node(inode)%values(1,1,var_Ti), T_min)
              flt_Te = max(node_list%node(inode)%values(1,1,var_Te), T_min)
            else
              flt_T  = max(node_list%node(inode)%values(1,1,var_T ), T_min)
              flt_Ti = 0.5d0 * flt_T
              flt_Te = 0.5d0 * flt_T
            endif
            flt_T = flt_Ti + flt_Te

            ! --- Lambda_0 = ln sqrt(m_i/(2*pi*m_e)), about 3 for deuterium
            if ( floating_Lambda .gt. 0.d0 ) then
              flt_lam0 = floating_Lambda
            else
              flt_lam0 = log( sqrt( central_mass * ATOMIC_MASS_UNIT / (2.d0*PI*MASS_ELECTRON) ) )
            endif

            if ( floating_Lambda_local ) then
              flt_lam      = flt_lam0 - 0.5d0 * log( GAMMA * flt_T / flt_Te )
              flt_dlam_dTi = - 0.5d0 / flt_T
              flt_dlam_dTe =   0.5d0 * flt_Ti / (flt_Te * flt_T)
            else
              flt_lam      = flt_lam0
              flt_dlam_dTi = 0.d0
              flt_dlam_dTe = 0.d0
            endif

            ! --- e*Phi/(k*Te) = (a_n*u/2 - vw)/Te, so Phi = Phi_float means u = 2*(Te*Lambda + vw)/a_n
            flt_a_n = - 2.d0 * EL_CHG * F0                                                          &
                      * sqrt( MU_ZERO * central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT ) &
                      / (central_mass * ATOMIC_MASS_UNIT)
            flt_vw  = EL_CHG * floating_V_wall * MU_ZERO * central_density * 1.d20

            flt_u      = 2.d0 * ( flt_Te * flt_lam + flt_vw ) / flt_a_n
            flt_du_dTi = 2.d0 *   flt_Te * flt_dlam_dTi / flt_a_n
            flt_du_dTe = 2.d0 * ( flt_lam + flt_Te * flt_dlam_dTe ) / flt_a_n

            ! --- AMPLITUDE RAMP (floating_amp_ramp) versus the old relaxation ramp.
            ! ---
            ! --- The old form multiplied flt_rel, the under-relaxation factor. Because the row's
            ! --- diagonal is never relaxed, the fixed point of that row is u = flt_u whatever
            ! --- flt_rel is - so scaling flt_rel never reduced the target, it only approached the
            ! --- FULL target more slowly. It was never a ramp on the boundary condition.
            ! ---
            ! --- It also measured time as t_now - t_start, i.e. from the original start of the
            ! --- simulation rather than from when the BC was switched on. After a restart at
            ! --- t >> floating_ramp_time that expression saturates at 1 on the very first step,
            ! --- so the ramp did nothing whatsoever.
            ! ---
            ! --- floating_start_time is the time at which the condition begins, and the amplitude
            ! --- multiplies the TARGET with a C2 smoothstep so the boundary data has no kink at
            ! --- either end. Leave floating_start_time < 0 to take t_start, the old behaviour.
            flt_rel = floating_u_relax
            flt_amp = 1.d0
            if ( floating_amp_ramp .and. (floating_ramp_time .gt. 0.d0) ) then
              flt_t0 = floating_start_time
              if ( flt_t0 .lt. 0.d0 ) flt_t0 = t_start
              flt_x  = max(0.d0, min(1.d0, (t_now - flt_t0) / floating_ramp_time))
              flt_amp = flt_x**3 * ( 10.d0 - 15.d0*flt_x + 6.d0*flt_x*flt_x )
            elseif ( floating_ramp_time .gt. 0.d0 ) then
              flt_rel = flt_rel * max(0.d0, min(1.d0, (t_now - t_start) / floating_ramp_time))
            endif

            ! --- OPTIONAL obliqueness gate, off by default and NOT required by the physics.
            ! --- Phi_float = Lambda*Te/e is a property of the sheath in front of a material
            ! --- surface and does not project with the incidence angle: a surface at 2 degrees
            ! --- floats at the same potential as one at 45. What projects is the current DENSITY
            ! --- through the wall, which is why a j-V sheath condition needs a gate like this and
            ! --- the floating potential does not. The Mach 1 condition just below follows the same
            ! --- convention - vpar = +-c_s over the whole boundary type, with the Chodura factor
            ! --- smoothing only the edges where b_n changes sign.
            ! --- Kept as a numerical experiment: if boundary structures show up on near-tangential
            ! --- stretches, gating them out isolates whether they originate there. The gated limit
            ! --- is exactly du = 0, i.e. the plain Dirichlet, because the diagonal is never relaxed.
            if ( floating_min_bn .gt. 0.d0 ) &
              flt_rel = flt_rel * bn**2 / ( bn**2 + floating_min_bn**2 )

            if ( flt_clip_Ti ) flt_du_dTi = 0.d0
            if ( flt_clip_Te ) flt_du_dTe = 0.d0

            ! --- GAUGE REMOVAL (floating_gauge_removal)
            ! ---
            ! --- The model uses u only through its derivatives - v_pol = R grad(u) x e_phi and
            ! --- w = Delta*u - so a spatially uniform shift of u is a gauge and drives nothing.
            ! --- It is NOT harmless as a boundary condition, though. Imposing u = C on the wall
            ! --- while the interior still sits near zero forces the offset to diffuse inwards
            ! --- through a boundary layer, and while that layer is thin
            ! ---     du/dn ~ C/h ,   w = Delta*u ~ C/h^2 ,
            ! --- which is a large transient vorticity source with no physical content. This is the
            ! --- mechanism behind the observation that a zero target is stable while a CONSTANT
            ! --- nonzero target is not - a constant target has no tangential electric field at all,
            ! --- so nothing physical distinguishes it from zero.
            ! ---
            ! --- Subtracting the wall-average removes exactly that component and leaves the physics
            ! --- untouched, since d/ds (q - qbar) = d/ds q. flt_gauge is one arclength-weighted mean
            ! --- over the WHOLE connected floating wall, computed in the previous sweep and reduced
            ! --- across ranks; taking a separate mean per boundary type would reintroduce a jump at
            ! --- every type transition.
            ! ---
            ! --- NOTE this is only legitimate while the wall draws no net current. A j-V sheath
            ! --- condition depends on ePhi/kTe absolutely, not just on its gradient, and must NOT
            ! --- use this.
            flt_gsum = flt_gsum + flt_u
            flt_gcnt = flt_gcnt + 1.d0
            if ( floating_gauge_removal ) flt_u = flt_u - flt_gauge

            ! --- Ramp the TARGET, not the relaxation (see above).
            flt_u      = flt_amp * flt_u
            flt_du_dTi = flt_amp * flt_du_dTi
            flt_du_dTe = flt_amp * flt_du_dTe

            flt_coef          = 0.d0
            flt_coef(var_u )  =   1.d0
            if ( with_TiTe ) then
              flt_coef(var_Ti) = - flt_du_dTi
              flt_coef(var_Te) = - flt_du_dTe
            else
              flt_coef(var_T ) = - 0.5d0 * ( flt_du_dTi + flt_du_dTe )
            endif
            flt_R = flt_u - node_list%node(inode)%values(1,1,var_u)

            ! --- TRACE DUMP: the imposed boundary data, node by node.
            ! ---
            ! --- The derivative row slaves du/dl to dTe/dl through the RAW nodal temperature
            ! --- derivative DOF. Te at the wall carries a natural BC, so that DOF is not
            ! --- constrained by anything, and any grid-scale noise in it is handed straight to the
            ! --- boundary electric field, to w = grad^2 u, and to the Mach-1 u_b/ps0_b term where
            ! --- it is amplified by 1/bn. Whether that noise actually exists is a question about
            ! --- this run, not about the formulation, so measure it rather than assume it: dump
            ! --- the target and the DOFs it is built from and look at them along the wall.
            if ( .not. flt_scanned ) then
              if ( .not. flt_tr_open ) then
                flt_tr_open = .true.
                write(flt_tr_file,'(A,I0,A)') 'bc_float_trace_', my_id, '.dat'
                open(newunit=flt_tr_unit, file=trim(flt_tr_file), status='replace',               &
                     action='write', iostat=bcdiag_ios)
                if ( bcdiag_ios .ne. 0 ) then
                  flt_tr_unit = -1
                else
                  flt_tr_ok = .true.
                  write(flt_tr_unit,'(A)') '#            R               Z  bnd_type    iv_dir'// &
                        '        target_u           u_val           u_dof          Te_val'//     &
                        '          Te_dof     target_dof'
                endif
              endif
              if ( flt_tr_ok ) then
                flt_tr_Te = 0.d0
                flt_tr_Td = 0.d0
                if ( with_TiTe ) then
                  flt_tr_Te = node_list%node(inode)%values(1,1,var_Te)
                  flt_tr_Td = node_list%node(inode)%values(1,iv_dir,var_Te)
                else
                  flt_tr_Te = node_list%node(inode)%values(1,1,var_T)
                  flt_tr_Td = node_list%node(inode)%values(1,iv_dir,var_T)
                endif
                flt_tr_tgt = 0.d0
                do k_flt = 1, n_var
                  if ( flt_coef(k_flt) .ne. 0.d0 .and. k_flt .ne. var_u )                          &
                    flt_tr_tgt = flt_tr_tgt - flt_coef(k_flt)                                      &
                               * node_list%node(inode)%values(1,iv_dir,k_flt)
                enddo
                write(flt_tr_unit,'(2f16.8,2i10,6es16.6)')                                        &
                  node_list%node(inode)%x(1,1,1), node_list%node(inode)%x(1,1,2),                 &
                  bnd_type, iv_dir, flt_u,                                                        &
                  node_list%node(inode)%values(1,1,var_u),                                        &
                  node_list%node(inode)%values(1,iv_dir,var_u),                                   &
                  flt_tr_Te, flt_tr_Td, flt_tr_tgt
              endif
            endif

            ! --- Which ROW carries the constraint. With a Dirichlet on w the natural place is the u
            ! --- row. With w free, JOREK's idiom is to put the condition on u into the w EQUATION
            ! --- instead and leave the u row to the vorticity equation - the same swap the code uses
            ! --- for "fixed psi but free zj". The constraint count is preserved either way: one
            ! --- condition, one row. That matters here because u is no longer constant along the
            ! --- wall, so Delta*u is not the frozen w and pinning both would be inconsistent.
            flt_row = var_u
            if ( .not. bcs(bnd_type)%dirichlet%w ) flt_row = var_w

            index_node  = node_list%node(inode)%index(1)
            index_node2 = node_list%node(inode)%index(iv_dir)

            do k_flt = 1, n_var
              if ( flt_coef(k_flt) .eq. 0.d0 ) cycle
              if ( k_flt .eq. var_u ) then
                call boundary_conditions_add_one_entry(                      &
                       index_node, flt_row, in, index_node, k_flt, in,       &
                       zbig * flt_coef(k_flt), index_min, index_max, a_mat)
              else
                call boundary_conditions_add_one_entry(                      &
                       index_node, flt_row, in, index_node, k_flt, in,       &
                       zbig * flt_rel * flt_coef(k_flt), index_min, index_max, a_mat)
              endif
            enddo

            if (in .eq. 1) then
              call boundary_conditions_add_RHS(                              &
                     index_node, flt_row, in, index_min, index_max, RHS_loc, &
                     zbig * flt_rel * flt_R, a_mat%i_tor_min, a_mat%i_tor_max)
            else
              call boundary_conditions_add_RHS(                              &
                     index_node, flt_row, in, index_min, index_max, RHS_loc, &
                     0.d0, a_mat%i_tor_min, a_mat%i_tor_max)
            endif

            ! --- The same constraint differentiated along the boundary, for the tangential
            ! --- derivative degree of freedom. Skipped when floating_u_value_only: that row slaves
            ! --- du/dl to dTe/dl, and du/dl IS the along-wall electric field, so grid-scale noise in
            ! --- the boundary temperature is handed straight to E, to the ExB flow, and to
            ! --- w = Delta*u, which differentiates it again. Without the row, u still equals the
            ! --- floating potential at every node - so Phi still varies with Te along the target,
            ! --- which is the physics - and only the interpolation between nodes is decided by the
            ! --- vorticity equation, which has dissipation.
            if ( .not. floating_u_value_only ) then

            flt_Rb = 0.d0
            do k_flt = 1, n_var
              if ( flt_coef(k_flt) .ne. 0.d0 ) &
                flt_Rb = flt_Rb - flt_coef(k_flt) * node_list%node(inode)%values(1,iv_dir,k_flt)
            enddo

            do k_flt = 1, n_var
              if ( flt_coef(k_flt) .eq. 0.d0 ) cycle
              if ( k_flt .eq. var_u ) then
                call boundary_conditions_add_one_entry(                      &
                       index_node2, flt_row, in, index_node2, k_flt, in,     &
                       zbig * flt_coef(k_flt), index_min, index_max, a_mat)
              else
                call boundary_conditions_add_one_entry(                      &
                       index_node2, flt_row, in, index_node2, k_flt, in,     &
                       zbig * flt_rel * flt_coef(k_flt), index_min, index_max, a_mat)
              endif
            enddo

            if (in .eq. 1) then
              call boundary_conditions_add_RHS(                              &
                     index_node2, flt_row, in, index_min, index_max, RHS_loc,&
                     zbig * flt_rel * flt_Rb, a_mat%i_tor_min, a_mat%i_tor_max)
            else
              call boundary_conditions_add_RHS(                              &
                     index_node2, flt_row, in, index_min, index_max, RHS_loc,&
                     0.d0, a_mat%i_tor_min, a_mat%i_tor_max)
            endif

            endif   !=== .not. floating_u_value_only

          endif   !=== apply_floating_u

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

          ! --- REGULARISED 1/ps0_b (mach1_psib_floor)
          ! ---
          ! --- This ExB term is identically zero under a Dirichlet u, because U0_b is then zero
          ! --- everywhere on the boundary. Any condition that gives u structure along the wall -
          ! --- floating_u among them - switches it on for the first time.
          ! ---
          ! --- ps0_b is the tangential derivative of psi along the wall, and from the definition
          ! --- of bn a few hundred lines above, ps0_b = bn * BigR * Btot * dl. So 1/ps0_b is 1/bn
          ! --- up to smooth factors: an amplification of about 30 where the field meets a target
          ! --- at 2 degrees, and of order 1e4 on the near-tangential stretches of the main chamber
          ! --- wall where bn ~ 1e-4. The existing Chodura smoothing does not protect this - it only
          ! --- engages on edges where bn CHANGES SIGN (bn_1*bn_2 < 0), not where bn is merely small.
          ! ---
          ! --- Replace 1/x by x/(x^2 + eps^2): identical for |x| >> eps, bounded by 1/(2 eps), and
          ! --- odd, so the sign of the term is preserved. mach1_psib_floor <= 0 keeps the raw form.
          if ( mach1_psib_floor .gt. 0.d0 ) then
            psib_inv = ps0_b / ( ps0_b**2 + mach1_psib_floor**2 )
          else
            psib_inv = 1.d0 / ps0_b
          endif
          ! --- mach1_exb_term = .false. removes this coupling outright rather than bounding it.
          ! --- That is a diagnostic, not a physical option: the term is the ExB contribution to the
          ! --- flow crossing the boundary and belongs in the Bohm condition. But it is identically
          ! --- zero under a Dirichlet u and becomes active only once u varies along the wall, so
          ! --- switching it off is the one clean way to attribute boundary structure to it rather
          ! --- than to the floating value itself. If structures persist with it off, the coupling
          ! --- is exonerated and the cause is in the imposed trace.
          if ( .not. mach1_exb_term ) psib_inv = 0.d0
          Mach1BC     = - Vpar0   + direction / Btot * factor  * cs0               + factor / Btot * BigR**2 * U0_b*psib_inv 
          Mach1BC_v   = - 1.0
          Mach1BC_T   =           + direction / Btot * factor  * cs0_T 
          Mach1BC_u   =                                                            + factor / Btot * BigR**2 * element_size_0*psib_inv 
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
            dMach1BC     = dMach1BC + factor / Btot * BigR**2 * U0_bb*psib_inv
            dMach1BC_ubb = + factor / Btot * BigR**2 * element_size_3*psib_inv
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

        endif   !=== apply_cs
        
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

! --- Report the floating amplitude actually in force. Two separate no-ops have already shipped
! --- in this boundary condition - a ramp that scaled the relaxation instead of the target, and a
! --- ramp measured from the original t_start so that it saturated instantly after a restart - and
! --- neither was visible from the output. The amplitude is what decides whether the condition is
! --- doing anything at all, so print it rather than leaving it to be inferred.
if ( any(bcs(:)%floating_u) .and. (my_id .eq. 0) ) then
  flt_amp = 1.d0
  if ( floating_amp_ramp .and. (floating_ramp_time .gt. 0.d0) ) then
    flt_t0 = floating_start_time
    if ( flt_t0 .lt. 0.d0 ) flt_t0 = t_start
    flt_x   = max(0.d0, min(1.d0, (t_now - flt_t0) / floating_ramp_time))
    flt_amp = flt_x**3 * ( 10.d0 - 15.d0*flt_x + 6.d0*flt_x*flt_x )
  endif
  write(*,'(A,es12.5,A,f8.5,A,es12.5)') ' FLOATING: t = ', t_now,                                 &
        '   amplitude = ', flt_amp, '   gauge ref = ', flt_gauge
  if ( flt_amp .lt. 1.d-3 ) write(*,'(A)')                                                        &
        '          NOTE amplitude is essentially zero - the floating BC is not acting. Check'//   &
        ' floating_start_time against the current time.'
endif

! --- Reduce the floating-target gauge reference over all ranks for use in the next sweep. Any
! --- uniform constant is a valid gauge, so the only thing that matters is that every rank ends up
! --- with the SAME value - hence the global reduction rather than a per-rank mean.
if ( floating_gauge_removal ) then
  call MPI_Allreduce(flt_gsum, flt_gs_glob, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(flt_gcnt, flt_gc_glob, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  if ( flt_gc_glob .gt. 0.d0 ) flt_gauge = flt_gs_glob / flt_gc_glob
endif
flt_gsum = 0.d0
flt_gcnt = 0.d0

! --- The scans above are one-off reports. Retire the one-shot only once the mach1 geometry
! --- block has actually been reached at least once - setting it unconditionally means an early
! --- call in which no local boundary element qualifies consumes the report and nothing is ever
! --- written. Announce the outcome so a silent run is diagnosable: no BCDIAG line at all means
! --- the binary predates this code, a line with 0 visits means the block is never reached.
if ( .not. flt_scanned ) then
  bcdiag_calls = bcdiag_calls + 1
  if ( my_id .eq. 0 ) then
    write(*,'(A,I0,A,I0)') ' BCDIAG [boundary_conditions]: call ', bcdiag_calls,                  &
          ', mach1 boundary-node visits on rank 0 = ', bcdiag_hits
  endif
  if ( (bcdiag_hits .gt. 0) .or. (bcdiag_calls .ge. 10) ) flt_scanned = .true.
endif
if ( bcdiag_ok ) then
  close(bcdiag_unit)
  bcdiag_ok = .false.
endif
if ( flt_tr_ok ) then
  close(flt_tr_unit)
  flt_tr_ok = .false.
endif
if ( allocated(bcdiag_vis) ) deallocate(bcdiag_vis)
if ( allocated(bcdiag_dir) ) deallocate(bcdiag_dir)
if ( allocated(bcdiag_fac) ) deallocate(bcdiag_fac)

return
end subroutine boundary_conditions 

end module mod_boundary_conditions
