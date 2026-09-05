!> Initialise input parameters and read the input namelist
subroutine initialise_parameters(my_id, filename)

use tr_module
use phys_module
use mod_plasma_functions, only: initialise_reference_parameters
use vacuum
use pellet_module
use live_data

implicit none

! --- Routine parameters
integer,                      intent(in) :: my_id
character(len=*),             intent(in) :: filename
real*8 :: vacuum_fraction, b_over_a, a_over_b

! --- Local variables
integer :: ierr,err,i
!> Effective resistivity floor and the corr_neg asymptote that usually sets it (see the NOTE
!! printed below): resistivity() evaluates eta at T_corr, which never goes below the asymptote.
real*8 :: eta_floor_cn, eta_floor_eff

! --- Namelist with input parameters.
namelist /in1/  tstep, nstep, tstep_n, nstep_n,                     &
                rst_hdf5, rst_hdf5_version, keep_current_prof,      &
                thermoelectric_ohm, thermoelectric_coef,            &
                keep_current_prof_confined, keep_current_psin_cutoff, &
                keep_current_psin_sig, keep_current_z_sig,          &
                keep_current_confine_strength,                      &
                eta, visco, visco_par, visco_par_par,               &
                restart, rst_format, regrid, bootstrap, write_ps,   &
                bootstrap_psin_cutoff,                              &
                regrid_from_rz, force_horizontal_Xline,             &
                n_R, n_Z, n_radial, n_pol, n_tht, n_flux,           &
                n_open, n_private, n_leg, n_leg_out, n_ext,         &
                n_outer, n_inner, n_up_priv, n_up_leg, n_up_leg_out,&
                n_tht_equidistant,                                  &
                psi_axis_init, XR_r, SIG_r, XR_tht, SIG_tht,        &
                XR_z, SIG_z, bgf_r, bgf_z, bgf_rpolar, bgf_tht,     &
                xr_closed,                                          &
                SIG_closed, SIG_open, SIG_private, SIG_theta,       &
                SIG_leg_0, SIG_leg_1, dPSI_open, dPSI_private,      &
                SIG_up_leg_0, SIG_up_leg_1, SIG_up_priv,            &
                SIG_outer, SIG_inner, SIG_theta_up,                 &
                dPSI_outer, dPSI_inner, dPSI_up_priv,               &
                nout, nout_projection, nout_particles,              &
                xr1, sig1, xr2, sig2,                               &
                R_begin, R_end, Z_begin, Z_end,                     &
                rect_grid_vac_psi,                                  &
                R_geo, Z_geo, amin, mf, fbnd, fpsi, mode,           &
                R_boundary, Z_boundary, psi_boundary, n_boundary,   &
                R_Z_psi_bnd_file,                                   &
                n_pfc, manipulate_psi_map,                          &
                Rmin_pfc, Rmax_pfc, Zmin_pfc, Zmax_pfc, current_pfc,&
                n_jropes,                                           &
                R_jropes, Z_jropes, w_jropes, current_jropes,       &
                extend_existing_grid, no_mach1_bc,                  &
                grid_to_wall, RZ_grid_inside_wall, eqdsk_psi_fact,  &
                RZ_grid_jump_thres,                                 &
                n_wall_blocks, n_ext_block, corner_block,           &
                n_ext_equidistant,                                  &
                n_block_points_left,  n_block_points_right,         &
                R_block_points_left,  R_block_points_right,         &
                Z_block_points_left,  Z_block_points_right,         &
                use_simple_bnd_types,                               &
                tokamak_device, thermalization,                     &
                F0,                                                 &
                gamma_stangeby,gamma_i_stangeby,gamma_e_stangeby,   &
                gamma_sheath, gamma_sheath_i, gamma_sheath_e,       &
                sheath_Lambda, sheath_V_wall, sheath_V_wall_asym, sheath_V_wall_R0, sheath_V_wall_dR, sheath_u_relax,       &
                sheath_u_relax_time,                                &
                sheath_wall_vel,                                    &
                sheath_u_align_psi,                                 &
                sheath_u_value_only,                                &
                sheath_u_exp_max, sheath_u_exp_min,                 &
                sheath_Lambda_local, sheath_X_min, sheath_smooth_dX,&
                sheath_sat_slope, sheath_sat_slope_e, sheath_weak_wmin, &
                sheath_psi_jacobian, sheath_v_perp, sheath_weak_ufade, &
                sheath_weak_detmin,                                 &
                sheath_dfdx_min, &
                sheath_wall_pen, sheath_zj_relax, &
                sheath_jsat_from_vpar, sheath_jsat_vpar_min,        &
                sheath_zj_ratio_max,                                &
                sheath_min_bn, sheath_ramp_time,                    &
                sheath_stiff_max, sheath_init_u, sheath_flux_sign,  &
                sheath_init_u_all, sheath_diag_R_split,             &
                T_min_sheath, sheath_weak_beta, sheath_weak_rmax,   &
                sheath_weak_relax,                                  &
                deuterium_adas, deuterium_adas_1e20,                &
                old_deuterium_atomic,                               &
                density_reflection,                                 &
                diag_mach1, mach1_diag_R_split,                     &
                mach_one_bnd_integral, Vpar_smoothing,              &
                Vpar_smoothing_coef,                                &
                zjz_0, zjz_1, zj_coef,                              &
                rho_0, rho_1, rho_coef,                             &
                rhon_0, rhon_1, rhon_coef,                          &
                T_0,   T_1,   T_coef,                               &
                Ti_0,  Ti_1,  Ti_coef,                              &
                Te_0,  Te_1,  Te_coef,                              &
                FF_0,  FF_1,  FF_coef,                              &
                ZK_par, ZK_i_par, ZK_e_par, ZK_par_max,             &
                ZK_perp, ZK_i_perp, ZK_e_perp, D_par, D_perp,       &
                V_pinch_gauss, V_pinch_psin, V_pinch_sig, v_pinch_file, &
                heatsource_e, heatsource_i, heatsource,             &
                particlesource, tauIC, Wdia,                        &
                eta_num, visco_num, visco_par_num,                  &
                D_perp_num,  D_perp_num_tanh,                       &
                D_perp_num_tanh_psin, D_perp_num_tanh_sig,          &
                ZK_perp_num, ZK_perp_num_tanh,                      &
                ZK_perp_num_tanh_psin, ZK_perp_num_tanh_sig,        &
                eta_num_T_dependent, visco_num_T_dependent,         &
                Dn_perp_num, time_evol_scheme,                      &
                ZK_i_perp_num, ZK_i_perp_num_tanh,                  &
                ZK_i_perp_num_tanh_psin, ZK_i_perp_num_tanh_sig,    &
                ZK_e_perp_num, ZK_e_perp_num_tanh,                  &
                ZK_e_perp_num_tanh_psin, ZK_e_perp_num_tanh_sig,    &
                pellet_amplitude, pellet_R, pellet_Z, pellet_phi,   &
                pellet_radius, pellet_sig, pellet_length,           &
                pellet_psi, pellet_delta_psi, pellet_density,       &
                pellet_velocity_R, pellet_velocity_Z, pellet_theta, &
                pellet_ellipse,                                     &
                central_density, central_mass,                      &
                pellet_particles, use_pellet,                       &
                ellip,tria_u,tria_l,quad_u,quad_l,                  &
                xampl,xwidth,xsig,xtheta,xshift,xleft, xpoint,      &
                forceSDN,                                           &
                xcase, SDN_threshold, D_perp_file, ZK_perp_file,    &
                ZK_e_perp_file, ZK_i_perp_file,                     &
                rho_file, T_file, Ti_file, Te_file, ffprime_file,   &
                rot_file, normalized_velocity_profile,              &
                freeboundary_equil, freeboundary,  freeb_change_indices, &
                resistive_wall,                                     &
                wall_resistivity, wall_resistivity_fact,            &
                bc_natural_open,                                    &
                use_mumps_eq, use_pastix_eq, use_strumpack_eq,      &
                use_mumps_prj, use_pastix_prj, use_strumpack_prj,   &
                use_mumps, mumps_ordering,                          &
                use_BLR_compression, epsilon_BLR, just_in_time_BLR, &
                use_pastix, use_murge, use_murge_element, use_wsmp, &
                n_tor_fft_thresh, use_strumpack,                    &
                refinement, force_central_node,                     &
                fix_axis_nodes,                                     &
                adaptive_time, equil, bench_without_plot,           &
                eta_T_dependent, visco_T_dependent, T_max_visco,    &
                zkpar_T_dependent, T_max_eta, T_max_eta_ohm,        & 
                heatsource_psin, heatsource_sig,                    &
                heatsource_e_psin, heatsource_e_sig,                &
                heatsource_i_psin, heatsource_i_sig,                &
                particlesource_psin, particlesource_sig,            &
                edgeparticlesource, edgeparticlesource_psin,        &
                edgeparticlesource_sig,                             &
                particlesource_gauss,    heatsource_gauss,          &
                heatsource_gauss_i,      heatsource_gauss_e,        &
                heatsource_gauss_psin,   heatsource_gauss_sig,      &
                heatsource_gauss_i_psin, heatsource_gauss_i_sig,    &
                heatsource_gauss_e_psin, heatsource_gauss_e_sig,    &
                particlesource_gauss_psin, particlesource_gauss_sig,&
                neutral_line_source,                                &
                neutral_line_R_start, neutral_line_Z_start,         &
                neutral_line_R_end,   neutral_line_Z_end,           &
                produce_live_data, gmres, gmres_max_iter,           &
                gmres_m, gmres_4, gmres_tol, iter_precon,           &
                pastix_pivot, max_steps_noUpdate,                   &
                keep_n0_const, linear_run, export_for_nemec,        &
                RMP_on, RMP_har_cos,RMP_har_sin,                    &
                RMP_growth_rate, RMP_ramp_up_time,                  &
                RMP_psi_cos_file, RMP_psi_sin_file,                 &
                V_0,V_1,V_coef, output_bnd_elements,                &
                n_limiter, R_limiter, Z_limiter,                    &
                first_target_point, last_target_point,              &
                R_Z_psi_bnd_file, wall_file,time_evol_scheme,       &
                spi_tor_rot, tor_frequency, spi_num_vol,            &
                NEO, neo_file, aki_neo_const, amu_neo_const,        &
                D_prof_neg_thresh, ZK_prof_neg_thresh, T_min,       &
                D_prof_imp_neg_thresh, D_prof_tot_neg_thresh,       &
                ZK_par_neg_thresh, D_imp_extra_neg_thresh,          &
                T_min_neg,rho_min_neg,implicit_heat_source,         &
                ne_SI_min, Te_eV_min, rn0_min,                      &
                D_neutral_x, D_neutral_y, D_neutral_p,              &
                neutral_reflection, rho_min,                        &
                corr_neg_temp_coef,                                 &
                corr_neg_dens_coef, D_prof_neg, ZK_prof_neg,        &  
                ZK_par_neg, ZK_par_neg_thresh,                      & 
                ZK_e_par_neg, ZK_i_par_neg, ZK_e_prof_neg, ZK_i_prof_neg,   &
                ZK_e_prof_neg_thresh, ZK_i_prof_neg_thresh,         &
                ZK_e_par_neg_thresh, ZK_i_par_neg_thresh,           &
                ns_deltaphi, ns_delta_minor_rad, ksi_ion, spi_rnd_seed, &
                ns_amplitude, ns_R, ns_Z, ns_phi, ns_radius,        &
                spi_Vel_Rref,spi_Vel_Zref, using_spi, n_spi, n_inj, &
                spi_Vel_RxZref, spi_quantity, spi_abl_model,        &
                ns_radius_ratio, ns_radius_min, spi_angle,          &
                spi_L_inj, spi_L_inj_diff, spi_abl_mag_reduction,   &
                spi_abl_history_old,                                &
                drift_distance, energy_teleported,                  &
                K_Dmv, A_Dmv, L_tube, V_Dmv, P_Dmv,                 &
                spi_Vel_diff, t_ns, JET_MGI, ASDEX_MGI,             &
                delta_n_convection, nimp_bg, output_prad_phi,       &
                RMP_on, RMP_har_cos,RMP_har_sin, spi_shard_file,    &
                spi_plume_file, spi_plume_hdf5,                     &
                RMP_growth_rate, RMP_ramp_up_time,                  &
                RMP_psi_cos_file, RMP_psi_sin_file,                 &
                Number_RMP_harmonics,RMP_har_cos_spectrum,          &
                RMP_har_sin_spectrum, imp_type, adas_dir, n_adas,   &
                index_main_imp,                                     &
                amix, amix_freeb, equil_accuracy, use_imp_adas,     &
                equil_accuracy_freeb, current_ref, FB_Ip_position,  &
                FB_Ip_integral, Z_axis_ref, FB_Zaxis_position,      &
                FB_Zaxis_derivative,FB_Zaxis_integral, start_VFB,   &
                n_feedback_current, n_feedback_vertical,            &
                n_iter_freeb, n_pf_coils, pf_coils, R_axis_ref,     &
                axis_srch_radius,                                   &
                starwall_equil_coils, freeb_equil_iterate_area,     &
                psi_offset_freeb, diag_coils, rmp_coils,            &
                voltage_coils, vert_FB_amp, find_pf_coil_currents,  &
                delta_psi_GS, newton_GS_fixbnd, newton_GS_freebnd,  &
                pastix_maxthrd, eta_ohmic, centralize_harm_mat,     &
                vert_FB_amp_ts, vert_FB_gain, vert_pos_file,        & 
                vert_FB_tact, start_VFB_ts, I_coils_max, rad_FB_amp,&
                autodistribute_modes, modes_per_family,             &
                mode_families_modes, n_mode_families,               &
                weights_per_family, autodistribute_ranks,           &
                ranks_per_family, treat_axis, Z_xpoint_limit, tgnum_rhoimp,     &
                tgnum_psi, tgnum_u, tgnum_zj, tgnum_w, tgnum_rho,   &
                tgnum_T, tgnum_Ti, tgnum_Te, tgnum_vpar, tgnum_rhon,&
                tgnum_nre, tgnum_AR, tgnum_AZ, tgnum_A3,            &
                tstep_particles, nstep_particles,                   &
                nsubstep_particles,                                 &
                filter_perp,    filter_hyper,    filter_par,        &
                filter_perp_n0, filter_hyper_n0, filter_par_n0,     &
                apply_dirichlet_proj, restart_particles,            &
                proj_collection_period,                             &
                part_group_configs, part_groups_in_use, valves,     &
                fluid_configs, init_particles_only,                 &
                find_RZ_nearby_iter, find_RZ_nearby_tol,            &
                min_sheath_angle, bcs, part_kill_ratio,             &
                use_sc, add_sources_in_sc, visco_sc_num,            &
                D_perp_sc_num, D_par_sc_num, ZK_perp_sc_num,        &
                ZK_par_sc_num, ZK_i_perp_sc_num, ZK_i_par_sc_num,   &
                ZK_e_perp_sc_num, ZK_e_par_sc_num, visco_par_sc_num,&
                Dn_pol_sc_num, Dn_p_sc_num, cte_current_FB_fact,    &
                D_perp_imp_sc_num, D_par_imp_sc_num,                &
                eta_num_prof, eta_num_psin_dependent, D_par_imp,    &
                D_perp_imp, spi_quantity_bg, pellet_density_bg,     &
                visco_par_heating, constant_imp_source,             &
                T_min_ZKpar,Ti_min_ZKpar,Te_min_ZKpar,T_min_eta,    &
                CARIDDI_mode, use_newton, maxNewton, gamma_Newton,  &
                alpha_Newton, vacuum_min, strumpack_matching,       &
                visco_old_setup, visco_heating, eta_coul_log_dep,   &
                export_polar_boundary, xpoint_search_tries,         &
                use_manual_random_seed, manual_seed,                &
                use_fixed_rng_value, fixed_rng_value,               &
                loop_voltage, export_aux_node_list,                 &
                use_zkperp_times_density, zkperp_density_floor


if (my_id .eq. 0) then

  ! --- Preset input parameters to reasonable default values.
  call preset_parameters()

  call vacuum_preset(my_id, freeboundary_equil, freeboundary, resistive_wall)
  
  ! --- Model-specific presets
  particlesource_psin = 100.d0
  
  ! --- Read input parameters from namelist.
  if (trim(filename) .ne. "__NO_FILENAME__" ) then
    open(42, file=filename, status='old', action='read', iostat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'ERROR: COULD NOT OPEN NAMELIST FILE "', trim(filename), '".'
      stop
    end if
    read(42,in1)
    close(42)
  else
    read(5,in1)
  endif

  if ( thermoelectric_ohm ) then
    write(*,*) 'NOTE: thermoelectric_ohm is ON. This adds the Braginskii thermal force'
    write(*,*) '      -(0.71/e)*grad_par(Te) to Ohm''s law, alongside the electron pressure'
    write(*,*) '      term model600 already carries. Two things it does: it makes Phi vary'
    write(*,*) '      ALONG field lines, so Phi is no longer a flux function and the ExB'
    write(*,*) '      drift acquires a component ACROSS flux surfaces; and it drives a'
    write(*,*) '      thermoelectric current between targets at different Te, which acts'
    write(*,*) '      again through eta*j_par.'
    write(*,*) '      INCOMPLETE PAIR: the reciprocal electron heat flux'
    write(*,*) '      -div(0.71*Te*j_par/e) is NOT implemented. Current and potential will'
    write(*,*) '      respond without the electron energy transport that accompanies them,'
    write(*,*) '      so total energy is not conserved to Braginskii order. Adequate for'
    write(*,*) '      establishing whether the potential structure appears; NOT adequate'
    write(*,*) '      for a quantitative energy or detachment result.'
    if ( (.not. any(bcs(:)%natural%zj)) .and. (.not. any(bcs(:)%sheath_zj_weak)) .and.      &
         any(bcs(:)%dirichlet%zj) ) then
      write(*,*) '      WARNING: zj is pinned at the wall (dirichlet%zj, no sheath_zj_weak and'
      write(*,*) '      no natural%zj), so the wall current cannot respond and the'
      write(*,*) '      thermoelectric drive will move the potential only.'
    endif
  endif

  ! --- keep_current_confine_strength multiplies the DEPTH of the mask, so values outside
  ! --- [0,1] are not a stronger cut: mask_eff = 1 - s*(1 - mask) goes NEGATIVE for s > 1
  ! --- wherever mask < 1 - 1/s, which reverses the sign of the current source in the SOL
  ! --- and private flux region rather than removing it. s < 0 amplifies it there instead.
  if ( (keep_current_confine_strength .lt. 0.d0) .or.                                     &
       (keep_current_confine_strength .gt. 1.d0) ) then
    write(*,*) 'ERROR: keep_current_confine_strength must lie in [0,1]. It is the fraction'
    write(*,*) '       of the keep_current_prof source that the confinement mask removes'
    write(*,*) '       outside the confined region: 1 = mask in full (the default and the'
    write(*,*) '       previous behaviour), 0 = no masking, identical to'
    write(*,*) '       keep_current_prof_confined = .false. Outside that range the applied'
    write(*,*) '       mask leaves the range [0,1] and the source is reversed or amplified'
    write(*,*) '       on open field lines instead of being suppressed.'
    stop
  endif
  if ( keep_current_prof_confined .and. (.not. keep_current_prof) ) then
    write(*,*) 'WARNING: keep_current_prof_confined is set but keep_current_prof is .false.,'
    write(*,*) '         so there is no current source for the mask to act on.'
  endif

  if ( ( n_tor .eq. 1 ) .and. freeboundary .and. (.not. freeboundary_equil) ) then
    write(*,*) 'WARNING: The parameter freeboundary is automatically changed to .false. since n_tor==1 and freeboundary_equil is .false.'
    freeboundary= .false.
  end if

  ! --- Calculate normalisation factor for MGI source (related to its toroidal shape)
  ns_tor_norm = ns_deltaphi * PI**0.5 * ERF(PI/ns_deltaphi)

  if (trim(R_Z_psi_bnd_file) .ne. 'none') then
    ! --- Open the file.
    OPEN(UNIT=243, FILE=R_Z_psi_bnd_file, FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=err)
    if ( err /= 0 ) then
      write(*,*) 'ERROR in initialise_parameters: Cannot open file '//TRIM(R_Z_psi_bnd_file)//'.'
      stop
    endif
    write(*,'(A)') ' boundary info from R_Z_psi_bnd_file: R_boundary, Z_boundary, psi_boundary '
    do i=1,n_boundary
      read(243,*) R_boundary(i),Z_boundary(i),psi_boundary(i)
      write(*,*) R_boundary(i),Z_boundary(i),psi_boundary(i)
    enddo
  endif

  ! --- Calculate JOREK gamma_sheath from gamma_stangeby if provided (otherwise the other way around)
  if ( with_TiTe ) then
    if (gamma_e_stangeby > -1.d89) then
      gamma_sheath_e = (gamma-1.d0) * (gamma_e_stangeby - 1.d0)
    else
      gamma_e_stangeby = gamma_sheath_e / (gamma-1.d0) + 1.d0
    end if
    if (gamma_i_stangeby > -1.d89) then
      gamma_sheath_i = (gamma-1.d0) * (gamma_i_stangeby - 1.d0 - gamma)
    else
      gamma_i_stangeby = gamma_sheath_i / (gamma-1.d0) + 1.d0 + gamma
    end if
  else
    if (gamma_stangeby > -1.d89) then
      gamma_sheath = (gamma-1.d0) * (0.5d0*gamma_stangeby - 1.d0 - 0.5d0*gamma)
    else
      gamma_stangeby = 2.d0 * ( gamma_sheath / (gamma-1.d0) + 1.d0 + 0.5d0 * gamma )
    end if
  end if

  ! --- Consistency checks for the sheath j-V boundary condition on the electric potential
  if ( any(bcs(:)%sheath_u) ) then

    if ( sheath_u_exp_min .ge. sheath_u_exp_max ) then
      write(*,*) 'ERROR: sheath_u_exp_min must be smaller than sheath_u_exp_max.'
      stop
    endif

    ! --- The lower clip is the electron saturation limit. At X = -Lambda the potential is zero
    ! --- and the current is the full electron thermal flux, which is where that limit belongs.
    if ( abs(sheath_u_exp_min + sheath_Lambda) .gt. 1.d-8 ) then
      write(*,*) 'NOTE: sheath_u_exp_min is not -sheath_Lambda. The lower clip is the electron'
      write(*,*) '      saturation limit and it sits at Phi = 0 when it equals -Lambda, which is'
      write(*,*) '      where electron saturation physically is. Present setting caps |j| at'
      write(*,*) '      (exp(-sheath_u_exp_min)-1) =', exp(-sheath_u_exp_min)-1.d0, ' times j_sat.'
    endif

    if ( sheath_wall_vel .gt. 0.d0 ) then
      write(*,*) 'NOTE: sheath_wall_vel > 0: psi is no longer frozen on the sheath boundary types.'
      write(*,*) '      It relaxes as dpsi/dt = -sheath_wall_vel*(dpsi/dn - its t_start value).'
      write(*,*) '      This exists because a potential varying along the wall drives an ExB flow'
      write(*,*) '      through it, which drags poloidal flux into a layer of width eta/v_n - some'
      write(*,*) '      10 microns at Spitzer resistivity, far below any divertor cell. Useful'
      write(*,*) '      values are around v_n = R du/dl at the target; too large lets the'
      write(*,*) '      equilibrium field diffuse away, so scan it and watch psi_bnd.'
    endif

    if ( sheath_u_relax_time .gt. 0.d0 ) then
      write(*,*) 'NOTE: sheath_u_relax_time > 0, so the boundary condition is under-relaxed by'
      write(*,*) '      min(1, tstep/sheath_u_relax_time) and sheath_u_relax is ignored. This'
      write(*,*) '      keeps the response time of the Te -> u -> ExB -> Te loop fixed as tstep'
      write(*,*) '      is ramped up through tstep_n.'
    endif

    do i = 1, max_bnd_types
      if ( .not. bcs(i)%sheath_u ) cycle

      ! --- The characteristic occupies the u row, so the plain Dirichlet on u is skipped there.
      ! --- The current stays Dirichlet: every attempt to leave psi, u, w or zj to its own weak
      ! --- form at a boundary node fails, because those forms are assembled without their
      ! --- surface terms (only rho, T, Ti, Te, rhon and vpar have natural BC support in
      ! --- mod_boundary_matrix_open). A free zj picks up the missing grad(psi).n integral as a
      ! --- spurious current of order grad(psi).n / h.
      if ( .not. bcs(i)%dirichlet%zj ) then
        write(*,*) 'ERROR: bcs(', i, ')%sheath_u needs bcs(', i, ')%dirichlet%zj = .true.'
        write(*,*) '       The current definition weak form has no surface term, so a free zj at'
        write(*,*) '       the boundary absorbs it as a spurious current.'
        stop
      endif
      if ( .not. bcs(i)%dirichlet%w ) then
        write(*,*) 'ERROR: bcs(', i, ')%sheath_u needs bcs(', i, ')%dirichlet%w = .true.'
        write(*,*) '       The w row carries the definition w = grad^2 u, which has no surface'
        write(*,*) '       term either.'
        stop
      endif
      if ( .not. bcs(i)%mach1 ) then
        write(*,*) 'WARNING: bcs(', i, ')%sheath_u is active but bcs(', i, ')%mach1 is .false.;'
        write(*,*) '         j_sat is evaluated with v_par = +-c_s/|B|, which assumes the Mach 1'
        write(*,*) '         condition holds at the same boundary.'
      endif
    enddo

    ! --- u is continuous along the boundary. A type with the BC adjacent to one without it puts
    ! --- a step of order Lambda*Te/e across a single element, and since v.n = R*du/dl that is a
    ! --- large artificial ExB flow through the wall right at the junction.
    !do i = 1, max_bnd_types
    !  if ( bcs(i)%sheath_u .or. .not. bcs(i)%dirichlet%u ) cycle
    !  write(*,*) 'NOTE: bcs(', i, ')%sheath_u is .false. while other types have it. If this'
    !  write(*,*) '      boundary type touches one that has it, u steps from ~Lambda*Te/e to 0'
    !  write(*,*) '      across one element there. Enable it on every type bounding the plasma.'
    !enddo

  endif

  ! --- Consistency checks for the sheath boundary conditions: the surface terms of the u, w and
  ! --- zj weak forms (mod_boundary_matrix_open) and the nodal sheath current on the zj row
  ! --- (mod_boundary_conditions). sheath_zj MUST be in this guard - it is a nodal row, so it needs
  ! --- none of the natural flags, and leaving it out meant none of its checks ever ran.
  if ( any(bcs(:)%natural%u) .or. any(bcs(:)%natural%w) .or. any(bcs(:)%natural%zj) &
       .or. any(bcs(:)%sheath_zj) .or. any(bcs(:)%sheath_zj_weak) ) then

    ! --- only the SURFACE terms live in that branch; the nodal sheath_zj rows do not
    if ( ( any(bcs(:)%natural%u) .or. any(bcs(:)%natural%w) .or. any(bcs(:)%natural%zj) &
           .or. any(bcs(:)%sheath_zj_weak) ) &
         .and. (.not. bc_natural_open) ) then
      write(*,*) 'ERROR: bcs(:)%natural%u / %w / %zj require bc_natural_open = .true.,'
      write(*,*) '       because the boundary integrals are only assembled in that branch of'
      write(*,*) '       construct_matrix.'
      stop
    endif

    if (any(bcs(:)%natural%u) .or. any(bcs(:)%sheath_zj) .or. any(bcs(:)%sheath_zj_weak)) then
      if (F0 == 0.d0 .or. central_mass <= 0.d0 .or. central_density <= 0.d0) then
        write(*,*) 'ERROR: sheath normalization needs nonzero F0 and positive reference mass/density.'
        stop 1
      endif
      if (sheath_Lambda > 30.d0) then
        write(*,*) 'ERROR: sheath_Lambda > 30 is outside the supported collection model.'
        stop 1
      endif
      if (sheath_sat_slope /= 0.d0 .or. sheath_sat_slope_e /= 0.d0 .or. &
          sheath_v_perp /= 0.d0 .or. sheath_dfdx_min /= 0.d0) then
        write(*,*) 'ERROR: remove sheath_sat_slope, sheath_sat_slope_e, sheath_v_perp and sheath_dfdx_min.'
        write(*,*) '       Artificial saturation tails/floors are not the bounded physical wall law.'
        write(*,*) '       See doc/WEAK_SHEATH_FIXES_2026-09-05.md before restarting an old run.'
        stop 1
      endif
      if (sheath_jsat_from_vpar .and. .not. with_vpar) then
        write(*,*) 'ERROR: sheath_jsat_from_vpar requires a model with Vpar; zero is not missing data.'
        stop 1
      endif
      if (sheath_jsat_vpar_min < 0.d0 .or. sheath_jsat_vpar_min > 1.d0) then
        write(*,*) 'ERROR: sheath_jsat_vpar_min must be in [0,1]. Default zero adds no ion source.'
        stop 1
      endif
      if (sheath_jsat_from_vpar .and. sheath_jsat_vpar_min > 0.d0) &
        write(*,*) 'WARNING: sheath_jsat_vpar_min adds an artificial ion-collection floor.'
      if (vpar_smoothing .and. vpar_smoothing_coef(2) <= 0.d0 .and. vpar_smoothing_coef(1) <= 0.d0) then
        write(*,*) 'ERROR: tanh incidence smoothing requires a positive width.'
        stop 1
      endif
      write(*,*) 'SHEATH: Phi=+F0*u; separate ion/electron collection; bounded electron saturation at Phi=Vwall.'
      write(*,*) '        sheath_X_min and sheath_smooth_dX no longer modify the forward current law.'
    endif
    if (any(bcs(:)%sheath_zj_weak)) then
      if (mach_one_bnd_integral) then
        write(*,*) 'ERROR: weak sheath uses its own Galerkin Mach1 projection; disable mach_one_bnd_integral.'
        stop 1
      endif
      if (n_order /= 3 .or. refinement) then
        write(*,*) 'ERROR: weak sheath trace currently supports cubic, unrefined elements only.'
        write(*,*) '       Higher-order and hanging-node trace transforms are not implemented.'
        stop 1
      endif
      if (sheath_zj_ratio_max > 0.d0 .or. sheath_weak_wmin > 0.d0 .or. sheath_weak_rmax > 0.d0) then
        write(*,*) 'ERROR: weak sheath requires sheath_zj_ratio_max=0, sheath_weak_wmin=0, sheath_weak_rmax=0.'
        write(*,*) '       Solution-dependent row gates and clipped quadrature residuals are not supported.'
        stop 1
      endif
      write(*,*) 'WARNING: zj trace replacement is NOT a validated tangent-wall closure.'
      write(*,*) '         Weak-sheath Mach1 uses field-aligned flow, without division by B.n.'
    endif

    do i = 1, max_bnd_types


      ! --- natural%w / natural%zj were refused while their Jacobians were wrong: the residuals
      ! --- depend on a normal derivative, and the trial loop produced columns only at the
      ! --- value/tangential DOFs, so the true entry was missing and a spurious one was added
      ! --- through psi_t. mod_boundary_matrix_open now splits grad(.).n into its tangential and
      ! --- normal halves and assembles the second at direction_perp(l), so both are correct.
      ! --- What they DO need is for the same variable's Dirichlet to be off: a surface integral
      ! --- only reaches rows at the value and tangential DOFs, which a Dirichlet overwrites.
      if ( bcs(i)%natural%zj .and. bcs(i)%dirichlet%zj ) then
        write(*,*) 'ERROR: bcs(', i, ')%natural%zj needs dirichlet%zj = .false. With the Dirichlet'
        write(*,*) '       on, the surface term only reaches rows the Dirichlet overwrites, so it'
        write(*,*) '       imposes nothing and zj at the wall stays frozen - which is exactly what'
        write(*,*) '       keeps the sheath j-V loop open.'
        stop
      endif
      if ( bcs(i)%natural%w .and. bcs(i)%dirichlet%w ) then
        write(*,*) 'ERROR: bcs(', i, ')%natural%w needs dirichlet%w = .false.'
        stop
      endif
      if ( bcs(i)%natural%zj .and. bcs(i)%dirichlet%psi ) then
        write(*,*) 'NOTE: bcs(', i, ')%natural%zj with dirichlet%psi = .true. That is the intended'
        write(*,*) '      combination: the Dirichlet pins psi''s value and tangential derivative,'
        write(*,*) '      while the surface term couples zj to psi''s NORMAL derivative, which is'
        write(*,*) '      free. That is the link that lets the wall current respond to the sheath.'
      endif

      ! --- Nodal sheath current on the zj row. This is the route that puts the sheath where it
      ! --- belongs: the u equation is assembled in STRONG form, so it has no boundary flux a
      ! --- surface term could replace, and adding one injects a spurious source at the wall.
      if ( bcs(i)%sheath_zj_weak ) then
        if ( sheath_weak_beta .gt. 0.d0 ) then
          write(*,*) 'ERROR: sheath_weak_beta > 0, but the penalty route it controlled has been'
          write(*,*) '       replaced by a Galerkin ROW REPLACEMENT and no longer exists. A'
          write(*,*) '       penalty could not work here at any coefficient: the zj equation is'
          write(*,*) '       missing its surface term once dirichlet%zj is released, so adding to'
          write(*,*) '       it cannot repair it (beta = 1e-2 left the wrong equation dominant,'
          write(*,*) '       beta = 1 gave period-2 divergence), and a Galerkin trace block spans'
          write(*,*) '       ~1e6 internally between value and derivative rows so a uniform large'
          write(*,*) '       coefficient wrecks the conditioning (beta = 1e9 died in one step).'
          write(*,*) '       Row normalisation by the trace mass diagonal replaces it. Remove the'
          write(*,*) '       parameter from the namelist.'
          stop
        endif
        if ( sheath_weak_detmin .lt. 0.d0 .or. sheath_weak_detmin .ge. 1.d0 ) then
          write(*,*) 'ERROR: sheath_weak_detmin must lie in [0,1). It is |sin(angle)| between'
          write(*,*) '       the two first-derivative DOF directions of a boundary node, so 1'
          write(*,*) '       would gate every row including the perfectly orthogonal ones.'
          write(*,*) '       0 disables the gate. Run util/check_boundary_frames.py on a'
          write(*,*) '       restart file to see what any threshold costs per boundary type.'
          stop
        endif
        if ( sheath_weak_detmin .gt. 0.d0 .and. .not. sheath_weak_ufade ) then
          write(*,*) 'NOTE: sheath_weak_detmin freezes u on the nodes it removes, so those'
          write(*,*) '      Gauss points sit at X = -Lambda and demand ~-19*j_sat, which leaks'
          write(*,*) '      into the FREE neighbour row unless sheath_weak_ufade fades them'
          write(*,*) '      out. Set sheath_weak_ufade = .true. alongside it.'
        endif
        if ( sheath_weak_relax .le. 0.d0 .or. sheath_weak_relax .gt. 1.d0 ) then
          write(*,*) 'ERROR: sheath_weak_relax must lie in (0,1]. It is the fraction of the'
          write(*,*) '       sheath mismatch the replacement row closes per step. <= 0 would'
          write(*,*) '       leave the row as M*d(zj) = 0, i.e. zj FROZEN on the boundary, and'
          write(*,*) '       > 1 over-corrects, which is the instability it exists to damp.'
          stop
        endif
        if (.not. bcs(i)%mach1 .and. .not. sheath_jsat_from_vpar) then
          write(*,*) 'WARNING: weak sheath without Mach1 uses a prescribed Bohm ion flux, not the solved Vpar.'
        endif
        if ( n_tor .gt. 1 ) then
          write(*,*) 'ERROR: bcs(', i, ')%sheath_zj_weak with n_tor > 1. The trace accumulator'
          write(*,*) '       stores one toroidal index per row. The sheath characteristic is'
          write(*,*) '       nonlinear, so a correct 3D form has to be evaluated in real toroidal'
          write(*,*) '       planes and transformed back, coupling harmonics through the row -'
          write(*,*) '       which this does not yet do. Verify axisymmetric first.'
          stop
        endif
        if ( bcs(i)%sheath_zj ) then
          write(*,*) 'ERROR: bcs(', i, ')%sheath_zj and %sheath_zj_weak are both set. They are'
          write(*,*) '       two ways to impose the SAME characteristic - one as a nodal row'
          write(*,*) '       replacement, one as a boundary penalty - and together the nodal row'
          write(*,*) '       overwrites zj after the penalty has been assembled, so the penalty'
          write(*,*) '       is silently discarded. Choose one.'
          stop
        endif
        if ( bcs(i)%dirichlet%zj ) then
          write(*,*) 'ERROR: bcs(', i, ')%sheath_zj_weak with dirichlet%zj = .true. The Dirichlet'
          write(*,*) '       freezes zj at these nodes, which overrides the penalty completely.'
          write(*,*) '       Set dirichlet%zj = .false. so the zj equation is assembled and the'
          write(*,*) '       penalty can act on it.'
          stop
        endif
        if ( bcs(i)%dirichlet%u ) then
          write(*,*) 'ERROR: bcs(', i, ')%sheath_zj_weak with dirichlet%u = .true. (the DEFAULT).'
          write(*,*) '       u is then frozen at this boundary, so the sheath characteristic is'
          write(*,*) '       evaluated at a fixed potential and there is no j-V feedback at all -'
          write(*,*) '       the whole point of the boundary condition. Set dirichlet%u = .false.'
          write(*,*) '       and keep a Dirichlet on u on at least one other type as the gauge.'
          stop
        endif
        if ( sheath_v_perp .le. 0.d0 ) then
          write(*,*) 'WARNING: bcs(', i, ')%sheath_zj_weak with sheath_v_perp = 0 (the DEFAULT).'
          write(*,*) '         zj_sat can then vanish where the field grazes the wall. The weak'
          write(*,*) '         replacement row cannot be switched off continuously by a uniform'
          write(*,*) '         weight, because that weight cancels during row normalisation.'
        endif
        if ( sheath_min_bn .gt. 0.d0 ) then
          write(*,*) 'NOTE: sheath_min_bn does not gate a weak replacement row. On this route it'
          write(*,*) '      only sets the active-area weight in the sheath diagnostics; use it to'
          write(*,*) '      interpret the output, not as a weak-row stabilisation parameter.'
        endif
        if ((.not. sheath_init_u) .and. sheath_ramp_time <= 0.d0) &
          write(*,*) 'NOTE: initialize from a compatible plasma-wall voltage; abrupt changes can launch flows.'
        if ( sheath_diag_R_split .le. 0.d0 ) then
          write(*,*) 'NOTE: bcs(', i, ')%sheath_zj_weak with sheath_diag_R_split = 0. The SHEATH'
          write(*,*) '      output then reports both targets summed, and an inner-outer'
          write(*,*) '      difference - which is what this boundary condition exists to produce'
          write(*,*) '      - cancels in that sum. Set it to about the X-point major radius.'
        endif
        if ( bcs(i)%natural%u ) then
          write(*,*) 'ERROR: bcs(', i, ')%sheath_zj_weak with natural%u = .true. The sheath'
          write(*,*) '       current is then already imposed through the u surface term, and the'
          write(*,*) '       Gauss-point characteristic the penalty needs is not evaluated on'
          write(*,*) '       that path. Use one route or the other.'
          stop
        endif
        write(*,*) 'NOTE: weak sheath replaces current trace rows with raw Galerkin moments, not a penalty.'
        write(*,*) '      An edge is selected only when BOTH endpoint types enable the weak sheath.'

      endif

      if ( bcs(i)%sheath_zj ) then
        if ( bcs(i)%natural%u ) then
          write(*,*) 'ERROR: bcs(', i, ')%sheath_zj and %natural%u are two different routes for the'
          write(*,*) '       same boundary condition. Enable only one.'
          stop
        endif
        if ( bcs(i)%sheath_u ) then
          write(*,*) 'ERROR: bcs(', i, ')%sheath_zj and %sheath_u both write the sheath; sheath_u'
          write(*,*) '       inverts the characteristic onto the u row and is singular at ion'
          write(*,*) '       saturation. Enable only one.'
          stop
        endif
        if ( bcs(i)%natural%zj ) then
          write(*,*) 'ERROR: bcs(', i, ')%sheath_zj and %natural%zj both target the zj rows.'
          write(*,*) '       Use natural%zj for a genuinely free wall current, sheath_zj to set it'
          write(*,*) '       from the sheath. Not both.'
          stop
        endif
        if ( .not. bcs(i)%dirichlet%zj ) then
          write(*,*) 'ERROR: bcs(', i, ')%sheath_zj needs dirichlet%zj = .true. The sheath rows'
          write(*,*) '       REPLACE the Dirichlet rows on zj (value and tangential derivative);'
          write(*,*) '       the flag marks which rows are taken over.'
          stop
        endif
        if ( .not. bcs(i)%mach1 ) then
          write(*,*) 'WARNING: bcs(', i, ')%sheath_zj without mach1; the saturation current uses'
          write(*,*) '         v_par = g(b_n)*c_s/|B|, which assumes the Mach 1 condition here.'
        endif
        if ( bcs(i)%dirichlet%u .and. (.not. sheath_init_u) ) then
          write(*,*) 'ERROR: bcs(', i, ')%sheath_zj with dirichlet%u = .true. but sheath_init_u ='
          write(*,*) '       .false. u is then FROZEN at whatever the restart left it at, which is'
          write(*,*) '       normally ~0 - and Phi = 0 is X = -Lambda, DEEP ELECTRON SATURATION.'
          write(*,*) '       The constraint would drive the wall current hard toward'
          write(*,*) '       (1-exp(Lambda))*j_sat, of order -19*j_sat, which is not a state you'
          write(*,*) '       want to relax towards. Observed: I_Ampere ran 868 -> -24560 A in 13'
          write(*,*) '       steps with max|j/jsat| reaching 200.'
          write(*,*) '       Set sheath_init_u = .true. so u starts at the local floating'
          write(*,*) '       potential, where the sheath carries no net current.'
          stop
        endif
        if ( bcs(i)%dirichlet%u ) then
          write(*,*) 'NOTE: bcs(', i, ')%sheath_zj with dirichlet%u = .true. The sheath current is'
          write(*,*) '      evaluated at a FROZEN potential, so there is no j-V feedback - useful'
          write(*,*) '      for relaxing the wall current toward a sheath-consistent state, not as'
          write(*,*) '      physics. Set dirichlet%u = .false. to let the potential respond,'
          write(*,*) '      keeping a Dirichlet on u on at least one other type as the gauge.'
        endif
        if ( (.not. bcs(i)%dirichlet%u) .and. (.not. sheath_init_u_all) ) then
          write(*,*) 'NOTE: bcs(', i, ')%sheath_zj with a Dirichlet u gauge on other boundary'
          write(*,*) '      types that is NOT initialised to the floating potential. Those nodes'
          write(*,*) '      are frozen at whatever u the restart holds - usually 0, i.e. Phi = 0,'
          write(*,*) '      X = -Lambda, deep electron saturation - while this boundary floats at'
          write(*,*) '      Lambda*Te/e. The step between them appears as grad(Phi) along the'
          write(*,*) '      wall, which IS ExB flow through it. Delta*u is ~0 in the interior, so'
          write(*,*) '      u is harmonic and the excursion grows with the free boundary area:'
          write(*,*) '      observed max ePhi/kTe 10.6 on one sheath type, 35.8 on two, then'
          write(*,*) '      blow-up. Set sheath_init_u = .true. AND sheath_init_u_all = .true. to'
          write(*,*) '      put the whole wall at the floating potential instead.'
        endif
        if ( .not. bc_natural_open ) then
          write(*,*) 'NOTE: bcs(', i, ')%sheath_zj without bc_natural_open. The boundary condition'
          write(*,*) '      itself is unaffected - it is a nodal row in mod_boundary_conditions -'
          write(*,*) '      but the SHEATH: wall-current diagnostic is evaluated in'
          write(*,*) '      boundary_matrix_open, which is only called in that branch, so it will'
          write(*,*) '      not be printed. Set bc_natural_open = .true. to get it.'
        endif
        if ( sheath_min_bn .le. 0.d0 ) then
          write(*,*) 'WARNING: bcs(', i, ')%sheath_zj with sheath_min_bn = 0. Where the field grazes'
          write(*,*) '         the wall the saturation current vanishes and the constraint would'
          write(*,*) '         drive zj to zero, which is wrong - zj = Delta*psi does not vanish'
          write(*,*) '         there. Use sheath_min_bn ~ 0.005 to hand those nodes back to the'
          write(*,*) '         frozen value.'
        endif
      endif

      if ( .not. bcs(i)%natural%u ) cycle

      if ( bcs(i)%dirichlet%u ) then
        write(*,*) 'ERROR: bcs(', i, ')%natural%u needs dirichlet%u = .false.; the whole point'
        write(*,*) '       is that u is free and set by charge continuity at the wall.'
        stop
      endif
      if ( bcs(i)%sheath_u ) then
        write(*,*) 'ERROR: bcs(', i, ')%natural%u and %sheath_u are two different implementations'
        write(*,*) '       of the same boundary condition. Enable only one of them.'
        stop
      endif
      if ( .not. bcs(i)%natural%zj ) then
        write(*,*) 'WARNING: bcs(', i, ')%natural%u without %natural%zj. zj at the wall is then'
        write(*,*) '         frozen by its Dirichlet, so the current the plasma delivers cannot'
        write(*,*) '         respond and the j-V loop is OPEN: the sheath can only adapt u, and'
        write(*,*) '         if the delivered current exceeds j_sat no u satisfies the'
        write(*,*) '         characteristic. Measured: I_Ampere constant to 4 digits while I_wall'
        write(*,*) '         ran to zero and u diverged. Set natural%zj = .true. and'
        write(*,*) '         dirichlet%zj = .false. on the same boundary types.'
      endif

      if ( .not. bcs(i)%mach1 ) then
        write(*,*) 'WARNING: bcs(', i, ')%natural%u is active but mach1 is .false.; the sheath'
        write(*,*) '         current uses v_par = g(b_n)*c_s/|B|, which assumes the Mach 1'
        write(*,*) '         condition holds at the same boundary.'
      endif

      ! --- The grazing-incidence floor is inherently DISCONTINUOUS: g_eff = sgn(b_n)*max(|g|,
      ! --- g_min) jumps from -g_min to +g_min as the field crosses tangency, so j_sat flips sign
      ! --- by a finite amount between neighbouring Gauss points. No continuous function can floor
      ! --- a magnitude while preserving a sign that changes.
      ! --- The natural form does not need it. The added term is -(zj_sh - zj)*(B.n), which
      ! --- vanishes with B.n whatever j_sat does, and dzj/du -> 0 with it, so the Robin diagonal
      ! --- vanishes smoothly too and u stays governed by the vorticity equation there. That is
      ! --- also the physics: a tangential field delivers no parallel flux, so there is no sheath.
      ! --- The floor exists for the NODAL path, where u is slaved to the characteristic and a
      ! --- vanishing j_sat makes the row singular.
      ! --- The ramp is a trap on this route. For a normal term, ramping it in from zero is a
      ! --- continuation. Here the surface term IS the boundary condition for u: with
      ! --- dirichlet%u = .false. there is nothing else holding u at this boundary, so a ramp
      ! --- factor near zero leaves u free to develop any along-wall gradient it likes - and
      ! --- du/dl is v_E.n, the velocity that drags poloidal flux. Observed: a ramp over the
      ! --- first ~100 steps grew a current filament on the divertor leg before the term ever
      ! --- engaged. sheath_init_u serves the purpose the ramp was meant to serve, by starting u
      ! --- at the floating potential, which is the fixed point of the characteristic at j = 0.
      ! --- ...but only for the ROBIN route. On the weak route the ramp multiplies the RELAX
      ! --- factor, not a surface term, and relax = 0 is exactly the Dirichlet freeze on zj, so
      ! --- every stage of the continuation is well posed and no fallback is needed.
      if ( (sheath_ramp_time .gt. 0.d0) .and. (sheath_wall_pen .le. 0.d0)                    &
           .and. any(bcs(:)%natural%u) ) then
        write(*,*) 'ERROR: sheath_ramp_time > 0 needs sheath_wall_pen > 0. The surface term is the'
        write(*,*) '       ONLY condition on u here, so ramping it from zero without the fallback'
        write(*,*) '       leaves u unconstrained at the wall and it drags flux freely. With the'
        write(*,*) '       fallback the two are complements - floating-potential relaxation at'
        write(*,*) '       ramp = 0, sheath at ramp = 1 - so u is always well posed.'
        stop
      endif
      if ( .not. sheath_init_u ) then
        write(*,*) 'WARNING: bcs(', i, ')%natural%u without sheath_init_u. u then starts ~Lambda*Te'
        write(*,*) '         away from its own fixed point, i.e. ~50 V at a 20 eV target.'
      endif

      if ( sheath_wall_pen .le. 0.d0 ) then
        write(*,*) 'WARNING: sheath_wall_pen = 0 with bcs(', i, ')%natural%u. Where the obliqueness'
        write(*,*) '         gate removes the sheath term, u then has NO boundary condition at all'
        write(*,*) '         and drifts freely - and du/dl is the flux-dragging velocity. With'
        write(*,*) '         grid_to_wall most of a sheath boundary is near-tangential wall, so'
        write(*,*) '         that is most of it. Use sheath_wall_pen ~ 0.05 (3 degree incidence).'
      endif

      if ( sheath_min_bn .le. 0.d0 ) then
        write(*,*) 'WARNING: sheath_min_bn = 0 with bcs(', i, ')%natural%u. Where the field grazes'
        write(*,*) '         the wall, zj_sat -> 0 while the frozen zj0 does not, so the'
        write(*,*) '         characteristic is asked for a current ratio outside its range and NO'
        write(*,*) '         u satisfies it - the residual cannot vanish and u runs away. Use'
        write(*,*) '         sheath_min_bn ~ 0.05 (the Chodura transition sits near |b_n| = 0.02).'
      endif

    enddo

    ! --- u enters the vorticity equation only through its gradient, so its constant mode is
    ! --- pinned by the sheath term alone, and that term loses its grip on u in ion saturation
    ! --- (dj/du -> 0). Keep at least one boundary type with a Dirichlet on u, as GBS does by
    ! --- using phi = Lambda*Te/e on the walls without strike points.
    ! --- Mixing the sheath and a Dirichlet on u between NEIGHBOURING target boundary types puts
    ! --- a step of Lambda*Te/e - some 50 V at a 20 eV target - across a single element wherever
    ! --- they meet. du/dl is v_E.n, so that step is the flux-dragging velocity at the largest
    ! --- value the mesh can represent, concentrated at a point. Observed: applying natural%u to
    ! --- one target type only drove w to 8692, five times worse than applying no boundary
    ! --- condition at all. mach1 marks the types that carry parallel outflow, i.e. the ones that
    ! --- have a sheath; if some of those have natural%u and others do not, say so.
    do i = 1, max_bnd_types
      if ( bcs(i)%mach1 .and. (.not. bcs(i)%natural%u) .and. (.not. bcs(i)%sheath_u) ) then
        write(*,*) 'WARNING: bcs(', i, ')%mach1 is set but it has no sheath BC, while another'
        write(*,*) '         boundary type does. Where the two meet, u steps by ~Lambda*Te/e'
        write(*,*) '         across one element and drags flux hard. Apply the sheath to every'
        write(*,*) '         boundary type carrying strike points, or to none.'
      endif
    enddo

    if ( .not. any( bcs(:)%dirichlet%u .and. .not. bcs(:)%natural%u .and. .not. bcs(:)%sheath_zj ) ) then
      write(*,*) 'ERROR: no boundary type pins u any more (every type has natural%u or sheath_zj,'
      write(*,*) '       or no Dirichlet). Setting the wall CURRENT leaves u a Neumann problem,'
      write(*,*) '       determined only up to a constant, so a reference is still required.'
      write(*,*) '       Keep dirichlet%u = .true. on at least one type (typically'
      write(*,*) '       the main chamber wall, where the sheath current is negligible anyway).'
      stop
    endif

    ! --- mod_boundary_matrix_open divides by vpar_smoothing_coef(2) without guarding it, and the
    ! --- sheath current uses that same smoothing function g(b_n) = normal_sign*factor
    if ( vpar_smoothing .and. (vpar_smoothing_coef(2) .le. 0.d0) ) then
      write(*,*) 'ERROR: vpar_smoothing = .true. with vpar_smoothing_coef(2) <= 0 divides by zero'
      write(*,*) '       in the boundary integrals, which the sheath current depends on.'
      write(*,*) '       Use the Chodura values (0.02, 0.016, 0.005754) or vpar_smoothing = .false.'
      stop
    endif

    if ( sheath_X_min .gt. 0.d0 ) then
      write(*,*) 'WARNING: sheath_X_min > 0 also limits the ion side of the characteristic,'
      write(*,*) '         which the forward form does not need. Use a negative value.'
    endif
    if ( sheath_min_bn .le. 0.d0 ) then
      write(*,*) 'NOTE: sheath_min_bn = 0, so the sheath current vanishes where the field is'
      write(*,*) '      tangent to the wall. That is correct for the forward form: no parallel'
      write(*,*) '      flux reaches such a surface, so there is no sheath, and the surface term'
      write(*,*) '      switches itself off smoothly through B.n. The floor is for the NODAL'
      write(*,*) '      path only, where it keeps the slaved row from going singular.'
    endif

  endif

  if (sum(nstep_n) .gt. 0) then
    nstep = sum(nstep_n)
    tstep = tstep_n(1)
  else
    tstep_n    = 0.d0
    tstep_n(1) = tstep
    nstep_n    = 0
    nstep_n(1) = nstep
  endif

  ! --- Fill the same ablation model to others if not specified to keep the old behavior
  do i = 2,n_inj
    if (spi_abl_model(i) < 0) then
      spi_abl_model(i) = spi_abl_model(1)
    end if
  end do

  call allocate_live_data()

endif


keep_n0_const  = ( keep_n0_const .or. linear_run )
! --- Read numerical profiles for rho, T, and ff'.
call read_num_profiles(my_id)

! --- Determine the derivatives of the numerical input profiles.
call derive_num_profiles(my_id)

! --- Initialize the shattered pellet position

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
spi_quantity_bg   = spi_quantity
pellet_density_bg = pellet_density
#endif

if ( my_id == 0 ) then
  if (2*PI/(n_tor*n_period) >= ns_deltaphi) then
    write(*,*) "WARNING! ns_deltaphi too small for the n_tor, BEWARE!"
    if (t_now > minval(t_ns)) then
      write(*,*) "EXITING NOW!!!"
      stop
    end if
  end if


  if (n_inj > n_inj_max .or. n_inj < 1) then
    write(*,*) "ERROR! Do not support n_inj larger than n_inj_max or smaller than 1, EXITING!"
    stop
  end if  

  do i = 1, n_inj_max
    if (n_spi(i)/=0 .and. i > n_inj) then
      write(*,*) "ERROR! Something wrong with n_inj, double check, EXITING!", n_spi, n_inj
      stop
    end if
  end do 
  
  if (n_adas > n_imp_max) then 
    write(*,*) "ERROR: n_adas should be no larger than n_imp_max, EXITING!"
    stop
  end if

  if (n_adas > 1 .and. (.not. use_imp_adas)) then
    write(*,*) "ERROR: Only support ADAS data for more than one impurities, through setting use_imp_adas to true, EXITING!"
    stop
  end if

  if (index_main_imp < 0 .or. index_main_imp > n_adas) then 
    write(*,*) "ERROR: Illegal value of index_main_imp, EXITING!"
    write(*,*) "ERROR: index_main_imp:", index_main_imp
    stop
  end if

  do i = 1,n_inj
    if (drift_distance(i) < 0.d0 .or. energy_teleported(i) < 0.d0) then 
      write(*,*) "ERROR: drift_distance and energy_teleported should be 0 or positive as signs already handled in codes, EXITING!"
      stop
    end if
  end do

  if (using_spi) call init_spi_all()
  ! ---- This is an ad hoc way to stop duplication of neutral source and ionized main species source, should find a permanenty solution later

  if (with_impurities .and. with_neutrals .and. using_spi .and. (maxval(spi_quantity_bg)>0.0)) then
    write(*,*) "WARNNING! Currently do not support both with_neutrals and with impurities enabled at the same time while injecting background species! Should be fixed sonn. EXITING!"
    stop
  endif

  call initialise_reference_parameters()

  ! --- Reported here, AFTER initialise_reference_parameters(), because that is where the
  ! --- T_min_eta sentinel is resolved to T_min. Printed earlier it would show -1.d12.
  ! --- THE EFFECTIVE RESISTIVITY FLOOR, and why T_min_eta often cannot move it.
  ! ---
  ! --- resistivity() evaluates eta_T at T_corr, which the CALLER has already passed through
  ! --- corr_neg_temp1 (mod_elt_matrix_fft passes T_or_Te_corr). corr_neg_temp1 asymptotes to
  ! --- L1 = T_min_neg*corr_neg_temp_coef(1) from above, so T_corr can NEVER go below L1 no
  ! --- matter how low T falls. The explicit clamp inside resistivity() keys on the RAW
  ! --- temperature, `T_raw < T_min_eta`, and only zeroes the derivative there - so whenever
  ! --- T_min_eta < L1 it cannot change the eta VALUE at all, and the floor is set by
  ! --- T_min_neg through corr_neg rather than by T_min or T_min_eta.
  ! ---
  ! --- This is easy to get wrong in both directions (a hand-set T_min_eta well below L1 has
  ! --- been tried and did nothing but leave deta_dT live and steep across [T_min_eta, 2*L1]),
  ! --- so print the number that actually applies instead of leaving it to be re-derived.
  if ( eta_T_dependent ) then
    eta_floor_cn  = T_min_neg * corr_neg_temp_coef(1)
    eta_floor_eff = max( T_min_eta, eta_floor_cn )
    write(*,'(A,es10.3,A,es10.3,A)')                                                  &
      ' NOTE: effective resistivity floor = ', eta_floor_eff,                          &
      '  (T_min_eta = ', T_min_eta, ')'
    write(*,'(A,es10.3,A)')                                                            &
      '       corr_neg asymptote T_min_neg*corr_neg_temp_coef(1) = ', eta_floor_cn, ''
    if ( T_min_eta .lt. eta_floor_cn ) then
      write(*,*) '      T_min_eta lies BELOW that asymptote, so it cannot lower the eta VALUE:'
      write(*,*) '      resistivity() sees T_corr, which never goes below the asymptote. Lower'
      write(*,*) '      T_min_neg if you want a genuinely more resistive cold leg - but note it'
      write(*,*) '      is the GLOBAL positivity floor and moves every T-dependent term with it.'
    endif
  endif

end if

return
end subroutine initialise_parameters
