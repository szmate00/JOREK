module mod_boundary_matrix_open
  implicit none
contains

subroutine boundary_matrix_open(vertex, direction, element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, &
                                psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, i_tor_min, i_tor_max, diagnostic_owner)
!---------------------------------------------------------------------
! calculates the matrix contribution of the boundaries of one element
! implements the natural boundary conditions
!---------------------------------------------------------------------
use constants
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use corr_neg
use mod_interp
use diffusivities, only: get_dperp, get_zkperp
use mod_sheath_bc, only: sheath_frame_det, sheath_frame_frozen, sheath_incidence, sheath_bohm_state
use mod_sheath_bc, only: sheath_current, sheath_norm, sheath_get_lambda, sheath_V_wall_at, &
                         sheath_temp_floor, dsheath_temp_floor_dT
use mod_sheath_diag, only: sheath_diag_add, sheath_diag_add_weak
use mod_sheath_trace, only: sheath_trace_add
use mod_sheath_geom_diag, only: sheath_geom_add

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)        ! the two nodes containing the boundary nodes
integer, intent(in)   :: i_tor_min   
integer, intent(in)   :: i_tor_max   
logical, intent(in), optional :: diagnostic_owner
logical :: collect_diagnostics

real*8     :: x_g(n_gauss), x_s(n_gauss), x_t(n_gauss), x_ss(n_gauss)
real*8     :: y_g(n_gauss), y_s(n_gauss), y_t(n_gauss), y_ss(n_gauss)

real*8     :: eq_g(n_plane,n_var,n_gauss), eq_s(n_plane,n_var,n_gauss), eq_p(n_plane,n_var,n_gauss)
real*8     :: eq_t(n_plane,n_var,n_gauss), eq_ss(n_plane,n_var,n_gauss)
real*8     :: delta_g(n_plane,n_var,n_gauss), delta_s(n_plane,n_var,n_gauss)

real*8     :: ELM(n_vertex_max*n_var*n_degrees*n_tor,n_vertex_max*n_var*n_degrees*n_tor)
real*8     :: RHS(n_vertex_max*n_var*n_degrees*n_tor)
real*8     :: rhs_ij(n_var), amat(n_var,n_var)
! --- Jacobian columns at the NORMAL-derivative degrees of freedom, see the block that
! --- assembles them below. Kept separate from amat because they land in different columns.
real*8     :: amat_p(n_var,n_var)
real*8     :: gpn_s, gpn_t, es_perp_sign
integer    :: index_kl_p

integer    :: vertex(2), direction(2), direction_perp(2), bnd_type1, bnd_type2
integer    :: i, j, j2, j3, ms, mt, mp, k, l, l2, l3, index_ij, index_kl, index, xcase2, is
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, ij8, kl1, kl2, kl3, kl4, kl5, kl6, kl7, kl8, i_var, k_var
real*8     :: ws, xjac,  dl, BigR, phi, eps_cyl, Btot
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
real*8     :: rhs_ij_5, rhs_ij_6, rhs_ij_7, rhs_ij_8
real*8     :: theta, zeta, Zbig, BB2, bdotn, gradvpar0dotn, gradvpardotn, factor, psi_ss, vpar_ss
real*8     :: R_inside, Z_inside, R_mid, Z_mid, R_cnt, Z_cnt, normal(2), normal_direction(2)
real*8     :: normal_sign, normal_sign3

real*8     :: v, v_x, v_y, v_s, v_p, v_ss, v_xx, v_yy, v_xs, v_ys
real*8     :: ps0, ps0_s, ps0_t, ps0_x, ps0_y, Vpar0, r0_corr, T0_corr, Ti0_corr, Te0_corr, cs0  
real*8     :: Ti0_sh, Te0_sh   ! temperatures for the SHEATH characteristic only, own floor
real*8     :: vpar0_s, vpar0_t, vpar0_x, vpar0_y 
real*8     :: vpar_s, vpar_t, vpar_x, vpar_y 
real*8     :: psi, psi_s, psi_t, vpar, T, Ti, Te, cs_T, cs_Ti, cs_Te
real*8     :: T0,   T0_s,  T0_t, T0_p
real*8     :: Ti0, Ti0_s, Ti0_t, Ti0_x, Ti0_y, Ti0_p
real*8     :: Te0, Te0_s, Te0_t, Te0_x, Te0_y, Te0_p
real*8     :: r0, r0_s, r0_t, r0_p, r0_x, r0_y, rho, rho_s, rho_t, rho_x, rho_y
real*8     :: c_1, c_2, c_3, c_angle, neutral_source
real*8     :: element_size_ij, element_size_kl, element_size_perp
real*8     :: grad_t(2), B0_R, B0_Z, factor_cs_bnd_integral
logical    :: xpoint2
integer    :: n_tor_local 
logical    :: apply_natural_bc(0:n_var)

! --- Charge-conserving sheath boundary condition and the two supporting surface terms
real*8     :: u0, u0_s, u0_t, u0_x, u0_y, zj0, sh_Bn, g_bn, sheath_ramp, sh_wgt_bn
real*8     :: sh_pen_c, sh_u_float, sh_an, sh_csat, sh_vw, sh_lam, sh_dlTi, sh_dlTe
real*8     :: sh_ramp_t, sh_act
! --- diagnostic-only evaluation for the nodal bcs%sheath_zj route. That route writes its rows in
! --- mod_boundary_conditions, which has no surface quadrature, so the wall-current diagnostic is
! --- evaluated here instead where dS is available and correctly weighted. Contributes nothing to
! --- rhs_ij or amat.
logical    :: diag_sheath_zj
logical    :: weak_sheath_zj   ! current trace projection at Gauss points
real*8     :: wk_res            ! raw current mismatch (with optional continuation)
real*8     :: wk_wgt  ! validity weight and the ratio it is built from
real*8     :: wk_gate         !< the VALIDITY weight alone, before the u-fade is applied
real*8     :: wk_Dv(2,2,n_tor) !< D_a carrying the validity weight ONLY. W = D/D0 conflates two
                              !! different questions once the fade multiplies wk_wgt: "is the
                              !! characteristic solvable here" and "is u free here". A row beside
                              !! a Dirichlet-u seam is faded ON PURPOSE and must not be gated as
                              !! though its physics had failed - type 4 is short with seams at
                              !! both ends, so its healthy W of 0.83 is mostly fade. sheath_weak_wmin
                              !! gates on Dv/D0 so it means only what its name says.
real*8     :: sh_ufree(2)     ! 1 where this node's u can respond, 0 where it is frozen
real*8     :: wk_wrx          ! wk_wgt * sheath_weak_relax, for the under-relaxed columns
! --- Galerkin trace diagnostics: D_a = int N_a N_a dS, F_a = int N_a (zj - zj_sh) dS,
! --- S_a = int N_a zj_sat dS. Purely diagnostic - see sheath_diag_add_weak.
real*8     :: wk_D(2,2,n_tor), wk_F(2,2,n_tor), wk_S(2,2,n_tor)
real*8     :: wk_D0(2,2,n_tor)   !< D_a with wk_wgt forced to 1: the row's UNWEIGHTED
                                 !! support, so W_a = D_a/D0_a is the mean validity
                                 !! weight over that support - scale free, and the
                                 !! only thing wk_wgt can still act through, since it
                                 !! cancels out of J_ab/D_a and F_a/D_a exactly.
integer    :: wk_i, wk_j, wk_m
! --- Galerkin trace block for the WEAK sheath row: J(a,b,var) = int N_a * dres/dvar * N_b dS,
! --- with a = (i,j) the test trace DOF and b = (k,l) the trial one. Axisymmetric only, so no
! --- toroidal index - initialise_parameters refuses n_tor_local > 1 with sheath_zj_weak.
real*8     :: tr_J(2,2,2,2,n_var), tr_F(2,2)
! Mach1 uses the SAME edge test space, instead of last-writer nodal corner rows.
logical :: weak_mach
real*8 :: ma_J(2,2,2,2,n_var), ma_Jp(2,2,2,2), ma_F(2,2), ma_S(2,2)
real*8 :: ma_v, ma_vb, ma_dT(2), ma_dbT(2), ma_dbTb(2), ma_dps, ma_dpt, ma_cs
real*8     :: tr_Jp(2,2,2,2)      !< psi column of the weak sheath row, on the
                                  !! direction_perp DOFs (the ones dirichlet%psi
                                  !! leaves free). Separate array because the
                                  !! existing tr_J columns all sit on direction().
real*8     :: dzj_dpt, dzj_dg, dzj_dB, dpx_dpt, dpy_dpt, dB_dpt, dbn_dpt
real*8     :: dzj_dps, dpx_dps, dpy_dps, dB_dps, dbn_dps
real*8     :: dg_dpt, dfac_dbn, zt_arg
integer    :: tr_col(4), tr_var(4), tr_k, tr_l, tr_nc
real*8     :: tr_vals(4)
real*8     :: dzj_sh, dzj_sat, dzj_x, dzj_d1, dzj_d2, dzj_d3, dzj_d4, dzj_d5, dzj_wgt
real*8     :: sh_duf_dTi, sh_duf_dTe
real*8     :: gradu0dotn, gradps0dotn, gradudotn, gradpsidotn
real*8     :: zj_sh, dzj_du, dzj_drho, dzj_dTi, dzj_dTe, dzj_dvpar, zj_sat_g, x_sheath
!> Node-frame determinant |x(1,2,:) x x(1,3,:)| of the row's OWN node, normalised to unit
!! vectors so it is |sin(angle between the two first-derivative DOF directions)|. The nodal
!! derivative basis is conditioned as 1/det, and zj = Delta*psi is built from SECOND
!! derivatives, so a near-degenerate frame inflates exactly the quantity the weak row
!! replaces. Nothing in the grid builder requires the two frame vectors to be independent.
real*8     :: sh_det
real*8     :: sheath_alpha, sh_d_pol, sh_d_robin

type (type_node)         :: tmp_node

theta = time_evol_theta
!zeta  = time_evol_zeta
! change zeta for variable dt
zeta  = time_evol_zeta * 2.0d0 * tstep / (tstep + tstep_prev)

Zbig = 1.d12

rhs_ij = 0.d0
amat   = 0.d0

c_angle = min_sheath_angle     * PI / 180.d0 ! --- angle factor for minimum heat and particle fluxes (in radians here)

!--------------------- reorder the nodes to have the same direction as full element (maybe not necesary)
if ((vertex(1) .eq. 3) .and. (vertex(2) .eq. 4)) then
  tmp_node  = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 4
  vertex(2) = 3
endif
if ((vertex(1) .eq. 4) .and. (vertex(2) .eq. 1)) then
  tmp_node  = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 1
  vertex(2) = 4
endif
if ((vertex(1) .eq. 3) .and. (vertex(2) .eq. 2)) then
  tmp_node  = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 2
  vertex(2) = 3
endif
if ((vertex(1) .eq. 2) .and. (vertex(2) .eq. 1)) then
  tmp_node  = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 1
  vertex(2) = 2
endif


!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; x_ss  = 0.d0; 
y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0; y_ss  = 0.d0; 
eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0; eq_ss = 0.d0; eq_p = 0.d0;

delta_g = 0.d0; delta_s = 0.d0;

direction_perp(1) = 6 / direction(2)     ! =3 if direction(2)=2, =2 if direction(2)=3
direction_perp(2) = 4

! --- Orientation of the transverse logical coordinate, exactly as applied to the FIELDS when
! --- eq_t is accumulated below. The trial-function block used to omit this flip, so on half the
! --- boundary elements its normal-derivative Jacobian had the wrong sign.
es_perp_sign = -1.d0
if ((vertex(1)*vertex(2) .eq. 2)) es_perp_sign = +1.d0

R_mid = sum(nodes(1:2)%x(1,1,1)) / 2.d0     ! mid point on boundary (approx.)
Z_mid = sum(nodes(1:2)%x(1,1,2)) / 2.d0
R_cnt = sum(nodes(1:4)%x(1,1,1)) / 4.d0     ! center point within element (approx.)
Z_cnt = sum(nodes(1:4)%x(1,1,2)) / 4.d0

normal_direction = (/R_mid - R_cnt, Z_mid - Z_cnt /) / norm2((/R_mid - R_cnt, Z_mid - Z_cnt /))

apply_natural_bc(:) = .false.
diag_sheath_zj      = .false.
weak_sheath_zj      = .false.

bnd_type1 = nodes(1)%boundary 
bnd_type2 = nodes(2)%boundary 

! --- A j-V relation is a statement about a surface FREE to respond. At a node whose u is
! --- Dirichlet-frozen with u = 0 and V_wall = 0, X = -Lambda exactly, so f = 1 - exp(Lambda)
! --- = -19.1: the Gauss points near it demand -19*j_sat of electron current and no potential
! --- change can answer. Skipping the ROW on that node (further down) does not remove those
! --- Gauss points from its FREE neighbour's row, because wk_F/wk_D/wk_S/tr_J are accumulated
! --- over the whole edge - so the neighbour gets dragged onto the electron branch.
! ---
! --- sheath_weak_detmin freezes u on a node whose frame is degenerate (mod_boundary_conditions
! --- applies the Dirichlet), so such a node is frozen for this purpose too and its Gauss points
! --- must be faded out of the neighbour's row by exactly the same mechanism.
sh_ufree(1) = 1.d0
if ( bcs(bnd_type1)%dirichlet%u .or.                                                     &
     sheath_frame_frozen(nodes(1)%x(1,2,1:2), nodes(1)%x(1,3,1:2)) ) sh_ufree(1) = 0.d0
sh_ufree(2) = 1.d0
if ( bcs(bnd_type2)%dirichlet%u .or.                                                     &
     sheath_frame_frozen(nodes(2)%x(1,2,1:2), nodes(2)%x(1,3,1:2)) ) sh_ufree(2) = 0.d0

! --- If one of the nodes has a boundary type where natural BCs are applied, apply boundary integral for the full bnd element
diag_sheath_zj = bcs(bnd_type1)%sheath_zj .or. bcs(bnd_type2)%sheath_zj
! Explicit coverage at BOTH ends: a type-3 corner must not sheath a 2--3 PFR edge.
weak_sheath_zj = bcs(bnd_type1)%sheath_zj_weak .and. bcs(bnd_type2)%sheath_zj_weak
weak_mach = weak_sheath_zj .and. with_vpar .and. (bcs(bnd_type1)%mach1 .or. bcs(bnd_type2)%mach1)
diag_sheath_zj = diag_sheath_zj .or. weak_sheath_zj
collect_diagnostics = .false.
if (present(diagnostic_owner)) collect_diagnostics = diagnostic_owner

do i_var=1, n_var
  if ( (i_var==var_rho ) .and. (bcs(bnd_type1)%natural%rho  .or. bcs(bnd_type2)%natural%rho ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_T   ) .and. (bcs(bnd_type1)%natural%T    .or. bcs(bnd_type2)%natural%T   ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_Ti  ) .and. (bcs(bnd_type1)%natural%Ti   .or. bcs(bnd_type2)%natural%Ti  ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_Te  ) .and. (bcs(bnd_type1)%natural%Te   .or. bcs(bnd_type2)%natural%Te  ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_rhon) .and. (bcs(bnd_type1)%natural%rhon .or. bcs(bnd_type2)%natural%rhon))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_vpar) .and. (bcs(bnd_type1)%natural%vpar .or. bcs(bnd_type2)%natural%vpar))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_u   ) .and. (bcs(bnd_type1)%natural%u    .or. bcs(bnd_type2)%natural%u   ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_w   ) .and. (bcs(bnd_type1)%natural%w    .or. bcs(bnd_type2)%natural%w   ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_zj  ) .and. (bcs(bnd_type1)%natural%zj   .or. bcs(bnd_type2)%natural%zj  ))  apply_natural_bc(i_var)=.true.
enddo

do i=1,2    ! sum over 2 verices
  
  do j=1,2  ! sum over two basis functions

    j2 = direction(j)
    element_size_ij = element%size(vertex(i),j2)

    j3 = direction_perp(j)
    element_size_perp = - element%size(vertex(i),direction_perp(1)) * 3.d0

    if ((vertex(1)*vertex(2) .eq. 2)) then
      element_size_perp = + element%size(vertex(i),direction_perp(1)) * 3.d0
    endif

    do ms=1, n_gauss

      x_g(ms)  = x_g(ms)  + nodes(i)%x(1,j2,1) * element_size_ij * H1(i,j,ms)
      x_s(ms)  = x_s(ms)  + nodes(i)%x(1,j2,1) * element_size_ij * H1_s(i,j,ms)
      x_t(ms)  = x_t(ms)  + nodes(i)%x(1,j3,1) * element_size_ij * H1(i,j,ms)   * element_size_perp

      y_g(ms)  = y_g(ms)  + nodes(i)%x(1,j2,2) * element_size_ij * H1(i,j,ms)
      y_s(ms)  = y_s(ms)  + nodes(i)%x(1,j2,2) * element_size_ij * H1_s(i,j,ms)
      y_t(ms)  = y_t(ms)  + nodes(i)%x(1,j3,2) * element_size_ij * H1(i,j,ms)   * element_size_perp

      do mp=1,n_plane

        do k=1,n_var

          do in=1,n_tor

            eq_g(mp,k,ms)  = eq_g(mp,k,ms)  + nodes(i)%values(in,j2,k) * element_size_ij * H1(i,j,ms)   * HZ(in,mp)
            eq_s(mp,k,ms)  = eq_s(mp,k,ms)  + nodes(i)%values(in,j2,k) * element_size_ij * H1_s(i,j,ms) * HZ(in,mp)
            eq_t(mp,k,ms)  = eq_t(mp,k,ms)  + nodes(i)%values(in,j3,k) * element_size_ij * H1(i,j,ms)   * HZ(in,mp) * element_size_perp
            eq_p(mp,k,ms)  = eq_p(mp,k,ms)  + nodes(i)%values(in,j2,k) * element_size_ij * H1(i,j,ms)   * HZ_p(in,mp)
            eq_ss(mp,k,ms) = eq_ss(mp,k,ms) + nodes(i)%values(in,j2,k) * element_size_ij * H1_ss(i,j,ms)* HZ(in,mp)

            delta_g(mp,k,ms) = delta_g(mp,k,ms) + nodes(i)%deltas(in,j2,k) * element_size_ij * H1(i,j,ms)   * HZ(in,mp)
            delta_s(mp,k,ms) = delta_s(mp,k,ms) + nodes(i)%deltas(in,j2,k) * element_size_ij * H1_s(i,j,ms) * HZ(in,mp)

          enddo
        enddo
      enddo

    enddo
  enddo
enddo

! changes deltas for variable time steps
delta_g = delta_g * tstep / tstep_prev
delta_s = delta_s * tstep / tstep_prev

n_tor_local = i_tor_max - i_tor_min +1

wk_D = 0.d0; wk_D0 = 0.d0; wk_Dv = 0.d0; wk_F = 0.d0; wk_S = 0.d0; tr_J = 0.d0; tr_F = 0.d0
tr_Jp = 0.d0
ma_J=0.d0; ma_Jp=0.d0; ma_F=0.d0; ma_S=0.d0

!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

  ws = wgauss(ms)

  dl   = sqrt(x_s(ms)**2 + y_s(ms)**2) 
  xjac = x_s(ms)*y_t(ms) - x_t(ms)*y_s(ms)
  BigR = x_g(ms)
  if (dl <= 0.d0 .or. xjac == 0.d0 .or. BigR <= 0.d0) then
    write(*,*) 'ERROR: singular boundary map (side, Gauss point, dl, xjac, R): ',vertex(1),ms,dl,xjac,BigR
    error stop 1
  endif

  grad_t = (/ - y_s(ms),   x_s(ms) /) / xjac

!  normal_direction = (/R_mid - R_cnt, Z_mid - Z_cnt /) / norm2((/R_mid - R_cnt, Z_mid - Z_cnt /))
  normal_direction = (/x_g(ms) - R_cnt, y_g(ms) - Z_cnt /) / norm2((/x_g(ms) - R_cnt, y_g(ms) - Z_cnt /))

  normal = dot_product(grad_t,normal_direction) * grad_t      ! outward pointing normal
  normal = normal / norm2(normal)

  ! --- grad(f).n split into the part carried by the TANGENTIAL derivative of f and the part
  ! --- carried by its NORMAL derivative:  grad(f).n = gpn_s * f_s + gpn_t * f_t . The two halves
  ! --- belong to different degrees of freedom, so the Jacobian has to keep them apart.
  gpn_s = (   y_t(ms) * normal(1) - x_t(ms) * normal(2) ) / xjac
  gpn_t = ( - y_s(ms) * normal(1) + x_s(ms) * normal(2) ) / xjac

  neutral_source = 0.d0

  ! --- Neutral sources at the boundary
  do is = 1, 10
    if     ( ((x_g(ms) - neutral_line_R_start(is))*(x_g(ms) - neutral_line_R_end(is)) .lt. 0.d0) &
       .and. ((y_g(ms) - neutral_line_Z_start(is))*(y_g(ms) - neutral_line_Z_end(is)) .lt. 0.d0) ) then
       neutral_source = neutral_source + neutral_line_source(is)
    endif
  enddo

  do mp = 1, n_plane

    ps0   = eq_g(mp,var_psi,ms)
    ps0_s = eq_s(mp,var_psi,ms) 
    ps0_t = eq_t(mp,var_psi,ms)   
    ps0_x = (   y_t(ms) * ps0_s - y_s(ms) * ps0_t ) / xjac
    ps0_y = ( - x_t(ms) * ps0_s + x_s(ms) * ps0_t ) / xjac

    B0_R =   ps0_y / x_g(ms)
    B0_Z = - ps0_x / x_g(ms)

    r0    = eq_g(mp,var_rho,ms)
    r0_s  = eq_s(mp,var_rho,ms)
    r0_t  = eq_t(mp,var_rho,ms)
    r0_p  = eq_p(mp,var_rho,ms)
    r0_x = (   y_t(ms) * r0_s - y_s(ms) * r0_t ) / xjac
    r0_y = ( - x_t(ms) * r0_s + x_s(ms) * r0_t ) / xjac

    if (with_TiTe) then
      Ti0    = eq_g(mp,var_Ti,ms)
      Ti0_s  = eq_s(mp,var_Ti,ms)
      Ti0_t  = eq_t(mp,var_Ti,ms)
      Ti0_p  = eq_p(mp,var_Ti,ms)
     
      Te0    = eq_g(mp,var_Te,ms)
      Te0_s  = eq_s(mp,var_Te,ms)
      Te0_t  = eq_t(mp,var_Te,ms)
      Te0_p  = eq_p(mp,var_Te,ms)

      T0     = Te0   + Ti0
      T0_s   = Te0_s + Ti0_s
      T0_t   = Te0_t + Ti0_t
      T0_p   = Te0_p + Ti0_p
    else
      T0     = eq_g(mp,var_T,ms)
      T0_s   = eq_s(mp,var_T,ms)
      T0_t   = eq_t(mp,var_T,ms)
      T0_p   = eq_p(mp,var_T,ms)

      Ti0    = T0    * 0.5d0  
      Ti0_s  = T0_s  * 0.5d0 
      Ti0_t  = T0_t  * 0.5d0 
      Ti0_p  = T0_p  * 0.5d0 

      Te0    = Ti0
      Te0_s  = Ti0_s
      Te0_t  = Ti0_t
      Te0_p  = Ti0_p
    endif


    Ti0_x = (   y_t(ms) * Ti0_s - y_s(ms) * Ti0_t ) / xjac
    Ti0_y = ( - x_t(ms) * Ti0_s + x_s(ms) * Ti0_t ) / xjac
    
    Te0_x = (   y_t(ms) * Te0_s - y_s(ms) * Te0_t ) / xjac
    Te0_y = ( - x_t(ms) * Te0_s + x_s(ms) * Te0_t ) / xjac

    if (with_vpar) then
      Vpar0   = eq_g(mp,var_vpar,ms)
      vpar0_s = eq_s(mp,var_vpar,ms) 
      vpar0_t = eq_t(mp,var_vpar,ms)   
      vpar0_x = (   y_t(ms) * vpar0_s - y_s(ms) * vpar0_t ) / xjac
      vpar0_y = ( - x_t(ms) * vpar0_s + x_s(ms) * vpar0_t ) / xjac
    else
      Vpar0   = 0.d0
      vpar0_s = 0.d0 
      vpar0_t = 0.d0 
      vpar0_x = 0.d0 
      vpar0_y = 0.d0 
    endif

    T0_corr  = corr_neg_temp1(T0)
    Ti0_corr = corr_neg_temp1(Ti0)
    Te0_corr = corr_neg_temp1(Te0)
    ! --- the sheath characteristic gets its own, milder floor: see sheath_temp_floor
    Ti0_sh   = sheath_temp_floor(Ti0)
    Te0_sh   = sheath_temp_floor(Te0)
    r0_corr  = corr_neg_dens(r0)

    if (with_TiTe) then
      cs0    = sqrt(gamma*(Ti0_corr+Te0_corr))
    else
      cs0    = sqrt(gamma*T0_corr)
    endif

    Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR

    BB2 = Btot**2

    bdotn = (+ ps0_y * normal(1) - ps0_x * normal(2)) / x_g(ms) / Btot
    gradvpar0dotn = (+ vpar0_x * normal(1) + vpar0_y * normal(2)) 

    normal_sign  = sign(1.d0,bdotn)
    normal_sign3 = sign(1.d0,ps0_s) * normal_sign

    c_1 = vpar_smoothing_coef(1); c_2 = vpar_smoothing_coef(2); c_3 = vpar_smoothing_coef(3)
    call sheath_incidence(bdotn, g_bn, dfac_dbn)
    factor = abs(g_bn)

    factor_cs_bnd_integral = 0.d0
    if (mach_one_bnd_integral) factor_cs_bnd_integral = 1.d0

    !-------------------------------------------------------------------------------------------
    ! --- Charge-conserving sheath boundary condition (see models/model600/mod_sheath_bc.f90)
    !
    ! --- The current term of the vorticity equation is assembled in strong form,
    ! ---     - v R (B.grad zj) xjac  ==  - v R div(zj B) xjac ,
    ! --- which is identically equal to its conservative form INCLUDING the surface integral
    ! ---     - oint v R zj (B.n) dl .
    ! --- Replacing the parallel current that leaves the domain by the current the sheath can
    ! --- actually pass therefore only requires adding the difference,
    ! ---     - oint v R ( zj_sheath - zj ) (B.n) dl ,
    ! --- to the u row. The mismatch is taken up by the perpendicular (polarisation and viscous)
    ! --- current of the same equation: charge piles up, the potential moves, the ExB response
    ! --- redistributes it. The term vanishes identically when the two currents agree.
    !
    ! --- The characteristic is evaluated in the forward direction, j = j(Phi), whose derivative
    ! --- tends smoothly to zero at ion saturation. Its contribution to amat has the same sign as
    ! --- the (negative definite) polarisation operator, i.e. it is dissipative.
    !-------------------------------------------------------------------------------------------
    u0    = eq_g(mp,var_u ,ms)
    u0_s  = eq_s(mp,var_u ,ms)
    u0_t  = eq_t(mp,var_u ,ms)
    u0_x  = (   y_t(ms) * u0_s - y_s(ms) * u0_t ) / xjac
    u0_y  = ( - x_t(ms) * u0_s + x_s(ms) * u0_t ) / xjac
    zj0   = eq_g(mp,var_zj,ms)

    gradu0dotn  = u0_x  * normal(1) + u0_y  * normal(2)
    gradps0dotn = ps0_x * normal(1) + ps0_y * normal(2)

    sh_Bn = bdotn * Btot                                          ! B.n, from the value already
                                                                  ! computed above for the fluxes
    g_bn  = normal_sign * factor                                  ! Chodura-Riemann g(b_n), signed

    ! --- Continuation factor, kept separately from the rest of sheath_ramp because the
    ! --- tangential-wall fallback below has to be its complement.
    sh_ramp_t = 1.d0
    if ( sheath_ramp_time .gt. 0.d0 ) &
      sh_ramp_t = max(0.d0, min(1.d0, (t_now - t_start) / sheath_ramp_time))
    sheath_ramp = sh_ramp_t

    zj_sh = 0.d0; dzj_du = 0.d0; dzj_drho = 0.d0; dzj_dTi = 0.d0; dzj_dTe = 0.d0
    dzj_dvpar = 0.d0
    if ( apply_natural_bc(var_u) ) then
      call sheath_current(u0, r0_corr, Ti0_sh, Te0_sh, g_bn, normal_sign, Btot, &
                          zj_sh, dzj_du, dzj_drho, dzj_dTi, dzj_dTe, zj_sat_g, x_sheath, &
                          sheath_V_wall_at(BigR), Vpar0, dzj_dvpar)
      ! --- chain rule through the corr_neg corrections, so the Jacobian stays exact where the
      ! --- density and temperature floors are active (which is exactly the cold divertor)
      dzj_drho = dzj_drho * dcorr_neg_dens_drho1(r0)
      dzj_dTi  = dzj_dTi  * dsheath_temp_floor_dT(Ti0)
      dzj_dTe  = dzj_dTe  * dsheath_temp_floor_dT(Te0)

      ! --- Cap the stiffness of the Robin term. Its diagonal scales as
      ! --- R*|dzj/du|*|B.n|*dl*theta*tstep while the row's own polarisation diagonal scales as
      ! --- rho*R^3; at a cold dense target the first can exceed the second by three orders of
      ! --- magnitude, which turns this boundary condition into a pointwise Dirichlet. u is then
      ! --- slaved to the local Te, rho and zj, so node-to-node variation in those is imprinted
      ! --- straight onto u - and grad(u) along the wall is ExB flow through it. Scaling the
      ! --- residual AND the Jacobian by the same alpha leaves the fixed point untouched while
      ! --- handing control of u back to the vorticity equation.
      sheath_alpha = 1.d0
      if ( sheath_stiff_max .gt. 0.d0 ) then
        sh_d_pol   = r0_corr * BigR**3
        sh_d_robin = BigR * abs(dzj_du) * abs(sh_Bn) * dl * max(theta,1.d-2) * tstep
        ! --- Smooth cap, alpha = 1/(1+r) with r = d_robin/(stiff_max*d_pol): alpha -> 1 for
        ! --- r << 1 and -> stiff_max*d_pol/d_robin for r >> 1, i.e. the same two limits as the
        ! --- hard switch it replaces, but C-infinity in between. The switch put a kink in the
        ! --- residual as a function of the state, which costs Newton iterations for nothing.
        ! --- NOTE the omitted d(alpha)/d(state) in the Jacobian is harmless: the missing term is
        ! --- proportional to alpha' * (zj_sh - zj0), and (zj_sh - zj0) -> 0 at the solution, so
        ! --- it vanishes there and the local convergence rate is preserved.
        sheath_alpha = 1.d0 / ( 1.d0 + sh_d_robin / max(sheath_stiff_max * sh_d_pol, 1.d-300) )
      endif
      ! --- Obliqueness gate. Where the field grazes the wall, zj_sat -> 0 while the frozen zj0
      ! --- does not, so the characteristic is asked for f = zj0/zj_sat far outside its range
      ! --- (-(exp(-X_min)-1), 1] and NO u satisfies it: the residual cannot vanish and u is driven
      ! --- without bound. Observed: max|zj0/zj_sat| = 1.3e3 with the gate off, and ePhi/kTe
      ! --- climbing linearly past 12 while the mean sat near Lambda. The weight is smooth, never
      ! --- changes sign, and multiplies residual AND Jacobian, so it removes the term where it is
      ! --- unsolvable instead of distorting it. sheath_min_bn = 0 recovers the ungated behaviour.
      sh_wgt_bn = 1.d0
      if ( sheath_min_bn .gt. 0.d0 ) &
        sh_wgt_bn = bdotn**2 / ( bdotn**2 + sheath_min_bn**2 )
      sh_act = sh_wgt_bn * sh_ramp_t          ! how much of the sheath is actually switched on here
      sheath_ramp = sheath_ramp * sheath_alpha * sheath_flux_sign * sh_wgt_bn  ! every Gauss point

      ! --- Tangential-wall fallback. The gate above removes the sheath term where the field
      ! --- grazes the wall - but with dirichlet%u = .false. that leaves u with NO boundary
      ! --- condition there at all, and du/dl is v_E.n, the flux-dragging velocity. With
      ! --- grid_to_wall the great majority of a sheath-enabled boundary IS near-tangential wall,
      ! --- so this is most of the boundary, not an edge case (measured: 98.7% gated off).
      ! --- Relax u toward the local FLOATING potential instead, with weight (1 - sh_act) and the
      ! --- stiffness the sheath itself would have at |b_n| = sheath_wall_pen. That is continuous
      ! --- with the sheath solution - a surface carrying no net current floats - so there is no
      ! --- step between the two regimes, and it is the physically right statement for a
      ! --- conducting surface that no parallel flux reaches.
      ! --- Only the u_float sensitivity is carried into the Jacobian below: d(stiffness)/dstate
      ! --- multiplies (u0 - u_float), which vanishes at this term's own fixed point.
      sh_pen_c = 0.d0; sh_u_float = 0.d0; sh_duf_dTi = 0.d0; sh_duf_dTe = 0.d0
      if ( sheath_wall_pen .gt. 0.d0 ) then
        call sheath_norm(sh_an, sh_csat, sh_vw, sheath_V_wall_at(BigR))
        call sheath_get_lambda(Ti0_corr, Te0_corr, sh_lam, sh_dlTi, sh_dlTe)
        ! --- dzj_du*B.n with g(b_n)*b_n replaced by the reference sheath_wall_pen^2 and f' = 1
        sh_pen_c   = sh_csat * r0_corr * sqrt(GAMMA*(Ti0_corr+Te0_corr))                     &
                   * 0.5d0 * sh_an / Te0_corr * sheath_wall_pen**2
        ! --- Complement of the ACTIVE sheath strength, not just of the gate. u then has a
        ! --- boundary condition at every point and at every stage of the continuation: pure
        ! --- floating-potential relaxation at ramp = 0, the sheath where it is switched on at
        ! --- ramp = 1, a smooth blend in between. The two agree when the sheath carries no net
        ! --- current, so handing over between them introduces no step. This is what makes
        ! --- sheath_ramp_time usable at all: ramping the sheath term alone would remove the only
        ! --- condition on u (dirichlet%u is .false. here) and let it drift freely.
        sh_pen_c   = sh_pen_c * ( 1.d0 - sh_act )
        sh_u_float = 2.d0 * ( Te0_corr * sh_lam + sh_vw ) / sh_an
        sh_duf_dTi = 2.d0 *   Te0_corr * sh_dlTi / sh_an
        sh_duf_dTe = 2.d0 * ( sh_lam + Te0_corr * sh_dlTe ) / sh_an
      endif

      ! --- wall current / potential diagnostic; dS is the toroidally integrated surface element
      if (collect_diagnostics) call sheath_diag_add(bnd_type1, zj_sh, zj0, zj_sat_g, x_sheath, u0, Te0_sh, sh_Bn, &
                           ws * dl * BigR * TWOPI / dble(n_plane), sh_wgt_bn,             &
                           sheath_V_wall_at(BigR), BigR)
    endif

    ! --- Same diagnostic for the nodal sheath_zj route. Evaluated at the Gauss point purely to be
    ! --- reported: the constraint itself is imposed nodally in mod_boundary_conditions, and the
    ! --- state used here (u0, r0_corr, Ti0_corr, Te0_corr, g_bn, Btot) is the same, so the numbers
    ! --- describe the same characteristic. The obliqueness weight matches the one the nodal block
    ! --- applies, so max|j/jsat| is reported only where the constraint is actually active.
    ! --- Zeroed unconditionally: the weak sheath term below reads these, and they are only
    ! --- filled when the nodal/diagnostic branch runs. A hard zero makes the term vanish rather
    ! --- than pick up whatever the previous Gauss point left behind.
    dzj_sh = 0.d0; dzj_d1 = 0.d0; dzj_d2 = 0.d0; dzj_d3 = 0.d0; dzj_d4 = 0.d0; dzj_d5 = 0.d0
    dzj_sat = 0.d0; wk_res = 0.d0; wk_wgt = 0.d0; wk_wrx = 0.d0
    dzj_dpt = 0.d0; dzj_dps = 0.d0; dzj_dg = 0.d0; dzj_dB = 0.d0

    if ( diag_sheath_zj .and. (.not. apply_natural_bc(var_u)) ) then
      call sheath_current(u0, r0_corr, Ti0_sh, Te0_sh, g_bn, normal_sign, Btot, &
                          dzj_sh, dzj_d1, dzj_d2, dzj_d3, dzj_d4, dzj_sat, dzj_x,      &
                          sheath_V_wall_at(BigR), Vpar0, dzj_d5, dzj_dg, dzj_dB)
      dzj_wgt = 1.d0
      if ( sheath_min_bn .gt. 0.d0 ) &
        dzj_wgt = bdotn**2 / ( bdotn**2 + sheath_min_bn**2 )
      if (collect_diagnostics) call sheath_diag_add(bnd_type1, dzj_sh, zj0, dzj_sat, dzj_x, u0, Te0_sh, sh_Bn, &
                           ws * dl * BigR * TWOPI / dble(n_plane), dzj_wgt,             &
                           sheath_V_wall_at(BigR), BigR)

      ! --- sheath_current differentiates w.r.t. the FLOORED state, so chain the floors through
      ! --- before these are used as Jacobian entries against the solution variables.
      dzj_d2 = dzj_d2 * dcorr_neg_dens_drho1(r0)
      dzj_d3 = dzj_d3 * dsheath_temp_floor_dT(Ti0)
      dzj_d4 = dzj_d4 * dsheath_temp_floor_dT(Te0)

      ! --- d(zj_sh)/d(psi). The whole psi dependence is through the field geometry:
      ! ---     Btot  = sqrt(F0^2 + ps0_x^2 + ps0_y^2)/BigR
      ! ---     bdotn = (ps0_y*n_R - ps0_x*n_Z)/(R*Btot)
      ! --- and dirichlet%psi pins psi's VALUE and TANGENTIAL derivative while leaving the NORMAL
      ! --- derivative free, so the only live DOFs are direction_perp, reached through ps0_t.
      ! --- Omitting this column was justified as "the geometry is frozen within a Newton
      ! --- iteration", but there is no Newton iteration - one construct_matrix and one solve per
      ! --- step - so it is a systematic per-step error that compounds with the sheath area.
      if ( sheath_psi_jacobian .and. Btot .gt. 1.d-30 ) then
        dpx_dpt = - y_s(ms) / xjac
        dpy_dpt = + x_s(ms) / xjac
        dB_dpt  = ( ps0_x*dpx_dpt + ps0_y*dpy_dpt ) / ( BigR**2 * Btot )
        dbn_dpt = ( dpy_dpt*normal(1) - dpx_dpt*normal(2) ) / ( x_g(ms) * Btot )               &
                  - bdotn * dB_dpt / Btot
        ! --- d(g_bn)/d(bdotn). g_bn = normal_sign*factor(|bdotn|) and normal_sign = sign(bdotn),
        ! --- so the two sign factors cancel and this is just d(factor)/d|bdotn|. Zero where the
        ! --- factor was clipped at 0, and zero without vpar_smoothing (factor is then constant).
        ! dfac_dbn is the derivative returned by the SAME incidence evaluator.
        dg_dpt  = dfac_dbn * dbn_dpt
        dzj_dpt = dzj_dg * dg_dpt + dzj_dB * dB_dpt
        dpx_dps = +y_t(ms)/xjac
        dpy_dps = -x_t(ms)/xjac
        dB_dps = (ps0_x*dpx_dps+ps0_y*dpy_dps)/(BigR**2*Btot)
        dbn_dps = (dpy_dps*normal(1)-dpx_dps*normal(2))/(BigR*Btot)-bdotn*dB_dps/Btot
        dzj_dps = dzj_dg*dfac_dbn*dbn_dps+dzj_dB*dB_dps
      endif
      if (weak_mach) then
        call sheath_bohm_state([Ti0,Te0], [0.d0,0.d0], g_bn/Btot, 0.d0, &
                               ma_v,ma_vb,ma_dT,ma_dbT,ma_dbTb)
        ma_cs=sqrt(GAMMA*(Ti0_sh+Te0_sh))
        ma_dps=0.d0; ma_dpt=0.d0
        if (sheath_psi_jacobian) then
          ma_dps=ma_cs*(dfac_dbn*dbn_dps/Btot-g_bn*dB_dps/Btot**2)
          ma_dpt=ma_cs*(dfac_dbn*dbn_dpt/Btot-g_bn*dB_dpt/Btot**2)
        endif
      endif

      ! Fixed test space. Dynamic gates change the equation and can delete live rows.
      ! Unsupported legacy gate parameters are rejected during initialization.
      wk_wgt = 1.d0

      wk_gate = wk_wgt          ! validity only; the fade below must not enter the gate

      ! --- Fade toward a frozen-u node. H1(i,1,ms) are the two VALUE basis functions along the
      ! --- edge and form a partition of unity, so this is 1 at a free node, 0 at a frozen one and
      ! --- linear-in-basis between. It varies WITHIN the row, which is the point: a weight that is
      ! --- uniform over a row's support cancels exactly out of F_a/D_a and J_ab/D_a and can fade
      ! --- nothing (see mod_sheath_trace). wk_D0 is deliberately left unweighted so W = D_a/D0_a
      ! --- reports the fade instead of hiding it.
      if ( sheath_weak_ufade )                                                             &
        wk_wgt = wk_wgt * ( sh_ufree(1)*H1(1,1,ms) + sh_ufree(2)*H1(2,1,ms) )

      ! --- Under-relaxation. The replacement row is M*d(zj) - sum_x C_x*d(x) = F, which asks for
      ! --- the FULL step onto the linearised characteristic every time. Scaling C_x and F - but
      ! --- NOT M, and not the D_a that normalises the row - turns it into
      ! --- zj_new = zj_old + relax*(zj_sh_lin - zj_old), the same structure as the nodal szj_rel.
      ! --- Needed because the relation is stiff and the linearisation is deliberately incomplete
      ! --- (no psi column, so |B| and g_bn are frozen within an iteration): the effective gain can
      ! --- exceed 1 and the state alternates between two branches on successive steps. relax = 1
      ! --- reproduces the un-relaxed row exactly.
      ! --- Under-relaxation only. sh_ramp_t is NOT applied here - see the residual below.
      wk_wrx = wk_wgt * sheath_weak_relax

      ! --- Continuation ramps the TARGET, not the COUPLING. wk_wrx scales the Jacobian columns
      ! --- and the residual together, so putting sh_ramp_t there drives the u COLUMN to zero as
      ! --- well - and on a weak type u has no other boundary condition (dirichlet%u = .false.,
      ! --- natural%u = .false.), so ramping the coupling down leaves u completely unconstrained.
      ! --- That is precisely the failure the Robin route needed sheath_wall_pen to avoid, and it
      ! --- was MEASURED here: at step 2 with sh_ramp_t = 0.02 the potential still collapsed and
      ! --- e-limited reached 100 % of the area by step 3 with max|j/jsat| only 2.25 - a tiny
      ! --- demand with the whole surface on the clamp, i.e. u adrift, not over-driven.
      ! ---
      ! --- Scaling the RESIDUAL alone gives J*dx = ramp*F: the full linearised characteristic
      ! --- still couples u to zj at every stage, so u is constrained throughout, and only the
      ! --- distance travelled toward the target per step is eased in.
      wk_res  = ( dzj_sh - zj0 ) * sh_ramp_t
      ! Integrate the RAW moment. Clipping individual Gauss residuals changes its roots.

      ! --- DIAGNOSTIC ONLY (mod_sheath_geom_diag). Nothing here feeds back into any row.
      ! --- Placed at the END of this block on purpose: wk_wgt is only final after the validity
      ! --- gate and the u-fade above, and it is what separates the Gauss points that actually
      ! --- contribute to a replaced row from those the endpoint OR merely reaches.
      ! ---
      ! --- bdotn is passed UNMODIFIED. It is already b_n = B.n/|B| (:421 divides by Btot; :468
      ! --- recovers dimensional B.n as bdotn*Btot; :521 compares it directly with sheath_min_bn).
      ! --- An earlier version divided by Btot a second time and understated the grazing angle by
      ! --- ~2.5-2.8x. g_bn and Btot go too, so zj_sat = c_sat*rho*(|g|cs + v_perp)/|B| decomposes
      ! --- exactly - g is a tanh of |b_n| (:429), not proportional to it, so it cannot be
      ! --- reconstructed from the angle. dS is the same measure sheath_diag_add is given, so
      ! --- areas are comparable between the two diagnostics.
      if (collect_diagnostics) call sheath_geom_add(bnd_type1, vertex(1), xjac, x_s(ms), y_s(ms), x_t(ms), y_t(ms), &
                           x_g(ms), y_g(ms), r0_corr, Ti0_corr, Te0_corr,                  &
                           abs(bdotn), cs0, g_bn, Btot, zj0, dzj_sat,                      &
                           ws * dl * BigR * TWOPI / dble(n_plane), wk_wgt)
    endif

    do i=1,2                ! loop over nodes

      do j=1,2              ! loop over basis functions

        j2 = direction(j)
        element_size_ij = element%size(vertex(i),j2)

        do im=i_tor_min, i_tor_max

          ! --- Reset here, not once at routine entry. Every component written below sits under a
          ! --- loop-invariant condition today, so nothing is stale - but that invariant is
          ! --- undocumented and one component added under an index-dependent condition would leak
          ! --- silently into the next iteration. Zeroing per iteration removes the hazard; the
          ! --- cost is n_var stores on boundary elements only.
          rhs_ij = 0.d0

          v   =  H1(i,j,ms) * element_size_ij * HZ(im,mp)         ! test function

          ! --- Sheath current: charge continuity closed at the wall (see the derivation above)
          if ( apply_natural_bc(var_u) ) &
            rhs_ij(var_u)  = - v * BigR * ( ( zj_sh - zj0 ) * sh_Bn * sheath_ramp        &
                                          + sh_pen_c * ( u0 - sh_u_float ) )              &
                                        * dl * tstep

          ! --- Surface term of the vorticity definition (w = Delta_pol u). REFUSED by
          ! --- initialise_parameters: its Jacobian needs columns for the normal-derivative DOFs,
          ! --- which the trial function loop below cannot produce (l2 = direction(l) only), so the
          ! --- term ends up effectively explicit and grows boundary structures. It is also
          ! --- unnecessary while dirichlet%w = .true., because these rows are overwritten anyway.
          ! --- Kept as the starting point for a correct implementation, see mod_sheath_bc.f90.
          if ( apply_natural_bc(var_w) ) &
            rhs_ij(var_w)  = + v * BigR * gradu0dotn * dl

          ! Raw weak moments, with the same trace basis and measure as row assembly.
          if ( diag_sheath_zj ) then
            wk_D(i,j,im-i_tor_min+1) = wk_D(i,j,im-i_tor_min+1)                        &
                                     + v * v * wk_wgt       * ws * dl / BigR
            wk_D0(i,j,im-i_tor_min+1) = wk_D0(i,j,im-i_tor_min+1)                      &
                                     + v * v                * ws * dl / BigR
            wk_Dv(i,j,im-i_tor_min+1) = wk_Dv(i,j,im-i_tor_min+1)                      &
                                     + v * v * wk_gate      * ws * dl / BigR
            wk_F(i,j,im-i_tor_min+1) = wk_F(i,j,im-i_tor_min+1)                        &
                                     + v * ( zj0 - dzj_sh ) * wk_wgt * ws * dl / BigR
            wk_S(i,j,im-i_tor_min+1) = wk_S(i,j,im-i_tor_min+1)                        &
                                     + v * abs(dzj_sat)     * wk_wgt * ws * dl / BigR
          endif

          ! Current and Mach equations share quadrature, not a penalty.
          if (weak_sheath_zj .and. .not. apply_natural_bc(var_u)) then
            tr_F(i,j)=tr_F(i,j)+v*wk_res*wk_wrx*ws*dl/BigR
            if (weak_mach) then
              ma_F(i,j)=ma_F(i,j)+v*(ma_v-Vpar0)*ws*dl/BigR
              ma_S(i,j)=ma_S(i,j)+v*abs(ma_v)*ws*dl/BigR
            endif
          endif

          ! --- Surface term of the current definition (zj = Delta*psi). REFUSED for the same
          ! --- reason as the w term above, and equally unnecessary: the frozen zj trace cancels
          ! --- exactly between the strong-form volume term and the sheath surface flux.
          if ( apply_natural_bc(var_zj) ) &
            rhs_ij(var_zj) = + v * gradps0dotn / BigR * dl

          ! --- Neutral sources
          if (with_neutrals) then
            rhs_ij(var_rhon) =  v * neutral_source * BigR * dl * tstep     
          endif

          ! --- Most B.C.s need vpar
          if (with_vpar) then

            ! --- Density reflection and minimum particle flux
            rhs_ij(var_rho)   = + v * density_reflection * r0      * vpar0 * ps0_s * normal_sign3 * tstep     &
                                - v * r0      * cs0 * BigR * dl * c_angle * tstep     ! particle flux at 1 degree angle  

            ! --- Sheath heat flux (c_angle for mininum heat fluxes at grazing angles)
            if (with_TiTe) then
              rhs_ij(var_Ti)  = - v * (gamma_sheath_i-1.d0) * r0 * Ti0 * vpar0 * ps0_s * normal_sign3 * tstep &
                                - v * (gamma_sheath_i-1.d0) * r0 * Ti0 * cs0    * BigR * dl * c_angle * tstep & 
                                - v * (GAMMA - 1.d0) * vpar0 * visco_par_heating * gradvpar0dotn * BigR * dl  * tstep  

              rhs_ij(var_Te)  = - v * (gamma_sheath_e-1.d0) * r0 * Te0 * vpar0 * ps0_s * normal_sign3 * tstep &
                                - v * (gamma_sheath_e-1.d0) * r0 * Te0 * cs0  * BigR * dl * c_angle   * tstep  
            else
              rhs_ij(var_T)   = - v * (gamma_sheath  -1.d0) * r0 * T0  * vpar0 * ps0_s * normal_sign3 * tstep &
                                - v * (gamma_sheath  -1.d0) * r0 * T0  * cs0    * BigR * dl * c_angle * tstep & 
                                - v * (GAMMA - 1.d0) * vpar0 * visco_par_heating * gradvpar0dotn * BigR * dl  * tstep  
            endif

            ! --- Mach=1 through boundary integral penalization method
            rhs_ij(var_vpar) = - v * (vpar0 * Btot * normal_sign - cs0 * factor) * dl * Zbig  * factor_cs_bnd_integral 

            ! --- Fluid neutral reflection
            if (with_neutrals) then 
              rhs_ij(var_rhon) = rhs_ij(var_rhon)                                                   &
                               + v * neutral_reflection * r0 * vpar0 * ps0_s * normal_sign3 * tstep &
                               + v * neutral_reflection * r0 * cs0 * BigR * dl * c_angle    * tstep ! particle flux at 1 degree angle  
            endif ! with_neutrals

          endif ! with_vpar
          index_ij = n_tor_local*n_var*n_degrees*(vertex(i)-1) + n_tor_local * n_var * (j2-1) + im - i_tor_min +1  ! index in the ELM matrix

          do i_var = 1, n_var
            ! --- var_zj is carried when the weak sheath term is on even though natural%zj is
            ! --- .false. - that flag selects the zj = Delta*psi surface term, a different thing
            ! --- which is refused because its Jacobian needs normal-derivative columns. The
            ! --- penalty below needs only value columns, so it is assemblable here.
            if ( .not. apply_natural_bc(i_var) ) cycle
            RHS(index_ij+(i_var-1)*(n_tor_local)) = RHS(index_ij+(i_var-1)*(n_tor_local)) + rhs_ij(i_var) * ws
          enddo


          do k=1,2                                                          ! loop over nodes

            do l=1,2                                                        ! loop over basis functions

              l2 = direction(l)
              element_size_kl = element%size(vertex(k),l2)

              ! --- direction_perp(l), NOT (j): the normal-derivative DOF that pairs with this
              ! --- trial basis function. l=1 (value) pairs with the normal first derivative,
              ! --- l=2 (tangential derivative) pairs with the mixed second derivative - exactly
              ! --- the pairing the field evaluation uses when it builds eq_t.
              l3 = direction_perp(l)
              element_size_perp = es_perp_sign * element%size(vertex(k),direction_perp(1)) * 3.d0

              do in = i_tor_min, i_tor_max                                              ! loop over toroidal harmonics

                amat   = 0.d0      ! see the note on rhs_ij above; amat is never accumulated onto
                amat_p = 0.d0      ! itself, so a per-iteration reset is always safe here

                psi    = H1(k,l,ms)    * element_size_kl * HZ(in,mp)
                psi_s  = H1_s(k,l,ms)  * element_size_kl * HZ(in,mp)
                psi_ss = H1_ss(k,l,ms) * element_size_kl * HZ(in,mp)
                psi_t  = H1(k,l,ms)    * element_size_kl * HZ(in,mp) * element_size_perp

                rho   = psi
                rho_s = psi_s
                rho_t = psi_t
                rho_x = (   y_t(ms) * rho_s - y_s(ms) * rho_t ) / xjac
                rho_y = ( - x_t(ms) * rho_s + x_s(ms) * rho_t ) / xjac

                T = psi; Ti = psi; Te = psi; vpar = psi; vpar_ss = psi_ss

                vpar_s = psi_s   
                vpar_t = psi_t
                vpar_x = (   y_t(ms) * vpar_s - y_s(ms) * vpar_t ) / xjac
                vpar_y = ( - x_t(ms) * vpar_s + x_s(ms) * vpar_t ) / xjac

                gradvpardotn  = (+ vpar_x * normal(1) + vpar_y * normal(2)) 

                gradudotn   = rho_x * normal(1) + rho_y * normal(2)   ! grad(trial function).n
                gradpsidotn = gradudotn                               ! (all trial functions equal)

                ! --- Jacobian of the sheath current term. amat = -d(rhs)/dx, and the diagonal
                ! --- entry adds to the negative definite polarisation operator, i.e. it damps.
                if ( apply_natural_bc(var_u) ) then
                  amat(var_u,var_u  ) = + v * BigR * ( dzj_du * sh_Bn * sheath_ramp + sh_pen_c ) &
                                            * psi * dl * theta * tstep
                  amat(var_u,var_zj ) = - v * BigR             * psi * sh_Bn * dl * theta * tstep * sheath_ramp
                  amat(var_u,var_rho) = + v * BigR * dzj_drho * psi * sh_Bn * dl * theta * tstep * sheath_ramp
                  ! --- j_sat built from Vpar (sheath_jsat_from_vpar) makes the row depend on the
                  ! --- parallel velocity too; dzj_dvpar is identically zero otherwise, so this
                  ! --- entry is inert in the default configuration.
                  if ( var_vpar .gt. 0 ) &
                    amat(var_u,var_vpar) = + v * BigR * dzj_dvpar * psi * sh_Bn * dl * theta * tstep * sheath_ramp
                  if ( with_TiTe ) then
                    amat(var_u,var_Ti) = + v * BigR * ( dzj_dTi * sh_Bn * sheath_ramp            &
                                                      - sh_pen_c * sh_duf_dTi )                &
                                             * psi * dl * theta * tstep
                    amat(var_u,var_Te) = + v * BigR * ( dzj_dTe * sh_Bn * sheath_ramp            &
                                                      - sh_pen_c * sh_duf_dTe )                &
                                             * psi * dl * theta * tstep
                  else
                    amat(var_u,var_T ) = + v * BigR * ( 0.5d0*(dzj_dTi+dzj_dTe)*sh_Bn*sheath_ramp &
                                                      - 0.5d0*sh_pen_c*(sh_duf_dTi+sh_duf_dTe) ) &
                                             * psi * dl * theta * tstep
                  endif
                endif

                ! --- Vorticity and current definition surface terms (algebraic equations: no
                ! --- theta, no tstep). Their residuals depend on grad(.).n, so their Jacobians
                ! --- need columns at BOTH the value/tangential DOFs (through f_s) and the
                ! --- normal-derivative DOFs (through f_t). The first go into amat as usual; the
                ! --- second into amat_p, which is assembled at l3 = direction_perp(l) below.
                ! --- Before this split, the whole grad(.).n was written into the value/tangential
                ! --- column using psi_t as if the value DOF carried a normal derivative: the true
                ! --- entry was missing and a spurious one was added in its place, leaving the
                ! --- terms effectively explicit with an O(1/h) coefficient.
                if ( apply_natural_bc(var_w) ) then
                  amat  (var_w,var_u)   = - v * BigR * gpn_s * psi_s * dl
                  amat_p(var_w,var_u)   = - v * BigR * gpn_t * psi_t * dl
                endif

                if ( apply_natural_bc(var_zj) ) then
                  amat  (var_zj,var_psi) = - v * gpn_s * psi_s / BigR * dl
                  amat_p(var_zj,var_psi) = - v * gpn_t * psi_t / BigR * dl
                endif

                ! --- Galerkin trace block. Same measure and test/trial functions as the rows
                ! --- the matrix uses, so the projection is consistent with the discretisation it
                ! --- will replace. Sign convention matches the nodal route: +1 on the zj column,
                ! --- -d(zj_sh)/dx on the state columns, residual (zj_sh - zj) on the RHS.
                if ( weak_sheath_zj .and. (.not. apply_natural_bc(var_u)) ) then
                  tr_J(i,j,k,l,var_zj ) = tr_J(i,j,k,l,var_zj )                        &
                                        + v * psi * wk_wgt * ws * dl / BigR
                  tr_J(i,j,k,l,var_u  ) = tr_J(i,j,k,l,var_u  )                        &
                                        - v * dzj_d1 * psi * wk_wrx * ws * dl / BigR
                  ! --- psi enters only through ps0_t, whose trial factor is psi_t. The DOF is
                  ! --- direction_perp, not direction, so it needs its own array and its own
                  ! --- column list at emission.
                  tr_Jp(i,j,k,l)        = tr_Jp(i,j,k,l)                                   &
                                        - v * dzj_dpt * psi_t * wk_wrx * ws * dl / BigR
                  tr_J(i,j,k,l,var_psi) = tr_J(i,j,k,l,var_psi) &
                                        - v * dzj_dps * psi_s * wk_wrx * ws * dl / BigR
                  if (weak_mach) then
                    ma_J(i,j,k,l,var_vpar)=ma_J(i,j,k,l,var_vpar)+v*psi*ws*dl/BigR
                    ma_J(i,j,k,l,var_psi)=ma_J(i,j,k,l,var_psi)-v*ma_dps*psi_s*ws*dl/BigR
                    ma_Jp(i,j,k,l)=ma_Jp(i,j,k,l)-v*ma_dpt*psi_t*ws*dl/BigR
                    if (with_TiTe) then
                      ma_J(i,j,k,l,var_Ti)=ma_J(i,j,k,l,var_Ti)-v*ma_dT(1)*psi*ws*dl/BigR
                      ma_J(i,j,k,l,var_Te)=ma_J(i,j,k,l,var_Te)-v*ma_dT(2)*psi*ws*dl/BigR
                    else
                      ma_J(i,j,k,l,var_T)=ma_J(i,j,k,l,var_T)-v*0.5d0*sum(ma_dT)*psi*ws*dl/BigR
                    endif
                  endif
                  tr_J(i,j,k,l,var_rho) = tr_J(i,j,k,l,var_rho)                        &
                                        - v * dzj_d2 * psi * wk_wrx * ws * dl / BigR
                  if ( var_vpar .gt. 0 ) &
                    tr_J(i,j,k,l,var_vpar) = tr_J(i,j,k,l,var_vpar)                    &
                                        - v * dzj_d5 * psi * wk_wrx * ws * dl / BigR
                  if ( with_TiTe ) then
                    tr_J(i,j,k,l,var_Ti) = tr_J(i,j,k,l,var_Ti)                        &
                                         - v * dzj_d3 * psi * wk_wrx * ws * dl / BigR
                    tr_J(i,j,k,l,var_Te) = tr_J(i,j,k,l,var_Te)                        &
                                         - v * dzj_d4 * psi * wk_wrx * ws * dl / BigR
                  else
                    tr_J(i,j,k,l,var_T ) = tr_J(i,j,k,l,var_T )                        &
                                         - v * 0.5d0*(dzj_d3+dzj_d4) * psi * wk_wrx * ws * dl / BigR
                  endif
                endif

                cs_T   = gamma * T  / (2.d0 * cs0)
                cs_Ti  = gamma * Ti / (2.d0 * cs0)
                cs_Te  = gamma * Te / (2.d0 * cs0)

                ! --- Most of natural BCs need vpar
                if (with_vpar) then

                  ! --- Density reflection and minimum particle flux (c_angle)
                  amat(var_rho,var_psi)   = - v * density_reflection * r0  * vpar0 * psi_s * normal_sign3 * theta * tstep 
                  amat(var_rho,var_rho)   = - v * density_reflection * rho * vpar0 * ps0_s * normal_sign3 * theta * tstep &
                                            + v                      * rho * cs0   * BigR * dl * c_angle  * theta * tstep 
                  amat(var_rho,var_vpar)  = - v * density_reflection * r0  * vpar  * ps0_s * normal_sign3 * theta * tstep 

                  ! --- Sheath heat flux
                  if (with_TiTe) then                
                    amat(var_rho,var_Ti)  = + v * r0 * cs_Ti * BigR * dl * c_angle  * theta * tstep
                    amat(var_rho,var_Te)  = + v * r0 * cs_Te * BigR * dl * c_angle  * theta * tstep
                  else
                    amat(var_rho,var_T)   = + v * r0 * cs_T  * BigR * dl * c_angle  * theta * tstep
                  endif

                  ! --- Sheath heat flux
                  if (with_TiTe) then                
                    amat(var_Ti,var_psi)  = + v * (gamma_sheath_i-1.d0) * r0  * Ti0 * vpar0 * psi_s * normal_sign3 * theta * tstep 
                    amat(var_Ti,var_rho)  = + v * (gamma_sheath_i-1.d0) * rho * Ti0 * vpar0 * ps0_s * normal_sign3 * theta * tstep & 
                                            + v * (gamma_sheath_i-1.d0) * rho * Ti0 * cs0   * BigR  * dl * c_angle * theta * tstep 
                    amat(var_Ti,var_Ti)   = + v * (gamma_sheath_i-1.d0) * r0  * Ti  * vpar0 * ps0_s * normal_sign3 * theta * tstep & 
                                            + v * (gamma_sheath_i-1.d0) * r0  * Ti  * cs0   * BigR  * dl * c_angle * theta * tstep &
                                            + v * (gamma_sheath_i-1.d0) * r0  * Ti0 * cs_Ti * BigR  * dl * c_angle * theta * tstep

                    amat(var_Te,var_psi)  = + v * (gamma_sheath_e-1.d0) * r0  * Te0 * vpar0 * psi_s * normal_sign3 * theta * tstep 
                    amat(var_Te,var_rho)  = + v * (gamma_sheath_e-1.d0) * rho * Te0 * vpar0 * ps0_s * normal_sign3 * theta * tstep & 
                                            + v * (gamma_sheath_e-1.d0) * rho * Te0 * cs0   * BigR  * dl * c_angle * theta * tstep 
                    amat(var_Te,var_Te)   = + v * (gamma_sheath_e-1.d0) * r0  * Te  * vpar0 * ps0_s * normal_sign3 * theta * tstep &
                                            + v * (gamma_sheath_e-1.d0) * r0  * Te  * cs0   * BigR  * dl * c_angle * theta * tstep &
                                            + v * (gamma_sheath_e-1.d0) * r0  * Te0 * cs_Te * BigR  * dl * c_angle * theta * tstep

                    amat(var_Ti,var_vpar) = + v * (gamma_sheath_i-1.d0) * r0  * Ti0 * vpar  * ps0_s * normal_sign3 * theta * tstep &
                                            + v * (GAMMA - 1.d0) * vpar * visco_par_heating * gradvpar0dotn * BigR * dl    * theta * tstep &
                                            + v * (GAMMA - 1.d0) * vpar0 * visco_par_heating * gradvpardotn * BigR * dl    * theta * tstep
                    amat(var_Te,var_vpar) = + v * (gamma_sheath_e-1.d0) * r0  * Te0 * vpar  * ps0_s * normal_sign3 * theta * tstep 
                  else
                    amat(var_T,var_psi)   = + v * (gamma_sheath  -1.d0) * r0  *  T0 * vpar0 * psi_s * normal_sign3 * theta * tstep 
                    amat(var_T,var_rho)   = + v * (gamma_sheath  -1.d0) * rho *  T0 * vpar0 * ps0_s * normal_sign3 * theta * tstep &
                                            + v * (gamma_sheath  -1.d0) * rho *  T0 * cs0   * BigR  * dl * c_angle * theta * tstep 
                    amat(var_T,var_T)     = + v * (gamma_sheath  -1.d0) * r0  *  T  * vpar0 * ps0_s * normal_sign3 * theta * tstep &
                                            + v * (gamma_sheath  -1.d0) * r0  *  T  * cs0   * BigR  * dl * c_angle * theta * tstep &
                                            + v * (gamma_sheath  -1.d0) * r0  *  T0 * cs_T  * BigR  * dl * c_angle * theta * tstep

                    amat(var_T,var_vpar)  = + v * (gamma_sheath  -1.d0) * r0  * T0  * vpar  * ps0_s * normal_sign3 * theta * tstep & 
                                            + v * (GAMMA - 1.d0) * vpar * visco_par_heating * gradvpar0dotn * BigR * dl    * theta * tstep &
                                            + v * (GAMMA - 1.d0) * vpar0 * visco_par_heating * gradvpardotn * BigR * dl    * theta * tstep
                  endif ! with_TiTe

                  ! --- Mach 1 condition through penalization boundary integral method
                  amat(var_vpar,var_vpar) =   v * (vpar * Btot * normal_sign) * dl * Zbig * factor_cs_bnd_integral

                  if (with_TiTe) then
                    amat(var_vpar,var_Ti) =   v * ( - cs_Ti) * factor         * dl * Zbig * factor_cs_bnd_integral
                    amat(var_vpar,var_Te) =   v * ( - cs_Te) * factor         * dl * Zbig * factor_cs_bnd_integral
                  else
                    amat(var_vpar,var_T)  =   v * ( - cs_T)  * factor         * dl * Zbig * factor_cs_bnd_integral
                  endif

                  ! --- Fluid neutral sources and reflection
                  if (with_neutrals) then
                    amat(var_rhon,var_psi) = - v * neutral_reflection * r0  * vpar0 * psi_s * normal_sign3      * theta * tstep 
  
                    amat(var_rhon,var_rho) = - v * neutral_reflection * rho     * vpar0 * ps0_s * normal_sign3      * theta * tstep &
                                             - v * neutral_reflection * rho     * cs0 * BigR * dl * c_angle * theta * tstep 
 
                    if (with_TiTe) then 
                      amat(var_rhon,var_Ti) = - v * neutral_reflection * r0 * cs_Ti * BigR * dl * c_angle * theta * tstep 
                      amat(var_rhon,var_Te) = - v * neutral_reflection * r0 * cs_Te * BigR * dl * c_angle * theta * tstep 
                    else
                      amat(var_rhon,var_T)  = - v * neutral_reflection * r0 * cs_T  * BigR * dl * c_angle * theta * tstep 
                    endif
  
                    amat(var_rhon,var_vpar) = - v * neutral_reflection * r0 * vpar  * ps0_s * normal_sign3     * theta * tstep 
                  endif ! with neutrals

                endif   ! with_vpar
                index_kl = n_tor_local*n_var*n_degrees*(vertex(k)-1) + n_tor_local * n_var * (l2-1) + in - i_tor_min +1  ! index in the ELM matrix

                ! --- same, at the normal-derivative degree of freedom l3 = direction_perp(l)
                index_kl_p = n_tor_local*n_var*n_degrees*(vertex(k)-1) + n_tor_local * n_var * (l3-1) + in - i_tor_min +1

                ! --- Add contributions to ELM matrix                 
                do k_var = 1, n_var
                  do i_var = 1, n_var

                    if ( .not. apply_natural_bc(i_var) ) cycle

                    ELM(index_ij+(i_var-1)*(n_tor_local),index_kl+(k_var-1)*(n_tor_local)) = &
                    ELM(index_ij+(i_var-1)*(n_tor_local),index_kl+(k_var-1)*(n_tor_local))   &
                      + amat(i_var,k_var) * ws

                    if ( amat_p(i_var,k_var) .ne. 0.d0 )                                       &
                      ELM(index_ij+(i_var-1)*(n_tor_local),index_kl_p+(k_var-1)*(n_tor_local)) = &
                      ELM(index_ij+(i_var-1)*(n_tor_local),index_kl_p+(k_var-1)*(n_tor_local))   &
                        + amat_p(i_var,k_var) * ws

                  enddo
                enddo


              enddo
            enddo
          enddo

        enddo
      enddo
    enddo

  enddo
enddo

! --- Hand the element's trace rows to the diagnostic. Normalising by D_a is what makes value
! --- rows and derivative rows comparable: the trace mass block scales as h, h^2, h^3 across the
! --- value/derivative combinations, so an unnormalised Galerkin row block spans orders of
! --- magnitude internally - which is why a large penalty on it fails where a pointwise Dirichlet,
! --- assigning the same zbig to both row types, does not.
if ( diag_sheath_zj .and. collect_diagnostics ) then
  do wk_m = 1, n_tor_local
    do wk_i = 1, 2
      do wk_j = 1, 2
        if ( wk_D(wk_i,wk_j,wk_m) .gt. 0.d0 )                                          &
          call sheath_diag_add_weak( wk_F(wk_i,wk_j,wk_m) / wk_D(wk_i,wk_j,wk_m),      &
                                     wk_S(wk_i,wk_j,wk_m) / wk_D(wk_i,wk_j,wk_m),      &
                                     wk_D(wk_i,wk_j,wk_m) )
      enddo
    enddo
  enddo
endif

! --- Hand this element's trace rows to the weak-sheath accumulator, which sums them with the
! --- adjacent edges' contributions, normalises by D_a and writes them with zbig in
! --- boundary_conditions. Global DOF indices come from nodes(i)%index(dof) - the same map the
! --- nodal route uses. Axisymmetric only: one toroidal index, checked at setup.
if ( weak_sheath_zj .and. (.not. apply_natural_bc(var_u)) ) then
  do wk_i = 1, 2

    ! Only explicitly selected, free-u nodes receive a current constraint.
    ! The legacy static frame exclusion is not a mesh repair.
    if (.not. bcs(nodes(wk_i)%boundary)%sheath_zj_weak) cycle
    if ( bcs(nodes(wk_i)%boundary)%dirichlet%u .or.                                    &
         sheath_frame_frozen(nodes(wk_i)%x(1,2,1:2), nodes(wk_i)%x(1,3,1:2)) ) cycle

    sh_det = sheath_frame_det( nodes(wk_i)%x(1,2,1:2), nodes(wk_i)%x(1,3,1:2) )

    do wk_j = 1, 2
      if ( wk_D(wk_i,wk_j,1) .le. 0.d0 ) cycle
      tr_nc = 0
      do tr_k = 1, 2
        do tr_l = 1, 2
          tr_nc = tr_nc + 1
          tr_col(tr_nc) = nodes(tr_k)%index(direction(tr_l))
        enddo
      enddo
      ! --- one variable at a time: the column list is the same four trace DOFs for each
      do tr_k = 1, n_var
        if ( (tr_k .ne. var_zj) .and. (tr_k .ne. var_u) .and. (tr_k .ne. var_psi) .and. (tr_k .ne. var_rho) .and.  &
             (tr_k .ne. var_Ti) .and. (tr_k .ne. var_Te) .and. (tr_k .ne. var_T)  .and.  &
             ( (tr_k .ne. var_vpar) .or. (.not. sheath_jsat_from_vpar) ) ) cycle
        if ( tr_k .eq. 0 ) cycle
        tr_nc = 0
        do tr_l = 1, 2
          do wk_m = 1, 2
            tr_nc = tr_nc + 1
            tr_col(tr_nc)  = nodes(tr_l)%index(direction(wk_m))
            tr_var(tr_nc)  = tr_k
            tr_vals(tr_nc) = tr_J(wk_i,wk_j,tr_l,wk_m,tr_k)
          enddo
        enddo
        if ( tr_k .eq. var_zj ) then
          call sheath_trace_add( nodes(wk_i)%index(direction(wk_j)), nodes(wk_i)%boundary, &
                                 wk_D(wk_i,wk_j,1), wk_D0(wk_i,wk_j,1),                 &
                                 wk_Dv(wk_i,wk_j,1),                                    &
                                 tr_F(wk_i,wk_j),                                         &
                                 wk_F(wk_i,wk_j,1), wk_S(wk_i,wk_j,1), sh_det,            &
                                 tr_nc, tr_col, tr_var, tr_vals )
        else
          call sheath_trace_add( nodes(wk_i)%index(direction(wk_j)), nodes(wk_i)%boundary, &
                                 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, sh_det,              &
                                 tr_nc, tr_col, tr_var, tr_vals )
        endif
      enddo

      ! --- The psi column, on the direction_perp DOFs (the ones dirichlet%psi leaves free).
      ! --- D/D0/F/Fd/S are already accumulated by the loop above, so this call adds columns only.
      if ( sheath_psi_jacobian ) then
        tr_nc = 0
        do tr_l = 1, 2
          do wk_m = 1, 2
            tr_nc = tr_nc + 1
            tr_col(tr_nc)  = nodes(tr_l)%index(direction_perp(wk_m))
            tr_var(tr_nc)  = var_psi
            tr_vals(tr_nc) = tr_Jp(wk_i,wk_j,tr_l,wk_m)
          enddo
        enddo
        call sheath_trace_add( nodes(wk_i)%index(direction(wk_j)), nodes(wk_i)%boundary, &
                               0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, sh_det,              &
                               tr_nc, tr_col, tr_var, tr_vals )
      endif
    enddo
  enddo
endif

! Do not gate the flow projection with a current-validity or frame threshold.
! A corner value receives moments from every selected incident material edge.
if (weak_mach) then
  do wk_i=1,2
    if (.not. bcs(nodes(wk_i)%boundary)%mach1) cycle
    sh_det=sheath_frame_det(nodes(wk_i)%x(1,2,1:2),nodes(wk_i)%x(1,3,1:2))
    do wk_j=1,2
      if (wk_D0(wk_i,wk_j,1) <= 0.d0) cycle
      do tr_k=1,n_var
        if (tr_k /= var_vpar .and. tr_k /= var_psi .and. tr_k /= var_T .and. &
            tr_k /= var_Ti .and. tr_k /= var_Te) cycle
        tr_nc=0
        do tr_l=1,2
          do wk_m=1,2
            tr_nc=tr_nc+1
            tr_col(tr_nc)=nodes(tr_l)%index(direction(wk_m))
            tr_var(tr_nc)=tr_k
            tr_vals(tr_nc)=ma_J(wk_i,wk_j,tr_l,wk_m,tr_k)
          enddo
        enddo
        if (tr_k == var_vpar) then
          call sheath_trace_add(nodes(wk_i)%index(direction(wk_j)),nodes(wk_i)%boundary, &
            wk_D0(wk_i,wk_j,1),wk_D0(wk_i,wk_j,1),wk_D0(wk_i,wk_j,1), &
            ma_F(wk_i,wk_j),-ma_F(wk_i,wk_j),ma_S(wk_i,wk_j),sh_det, &
            tr_nc,tr_col,tr_var,tr_vals,equation=var_vpar)
        else
          call sheath_trace_add(nodes(wk_i)%index(direction(wk_j)),nodes(wk_i)%boundary, &
            0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,sh_det,tr_nc,tr_col,tr_var,tr_vals,equation=var_vpar)
        endif
      enddo
      if (sheath_psi_jacobian) then
        tr_nc=0
        do tr_l=1,2
          do wk_m=1,2
            tr_nc=tr_nc+1
            tr_col(tr_nc)=nodes(tr_l)%index(direction_perp(wk_m))
            tr_var(tr_nc)=var_psi
            tr_vals(tr_nc)=ma_Jp(wk_i,wk_j,tr_l,wk_m)
          enddo
        enddo
        call sheath_trace_add(nodes(wk_i)%index(direction(wk_j)),nodes(wk_i)%boundary, &
          0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,sh_det,tr_nc,tr_col,tr_var,tr_vals,equation=var_vpar)
      endif
    enddo
  enddo
endif

return

end subroutine

end module mod_boundary_matrix_open
