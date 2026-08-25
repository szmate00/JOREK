module mod_boundary_matrix_open
  implicit none
contains

subroutine boundary_matrix_open(vertex, direction, element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, &
                                psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, i_tor_min, i_tor_max)
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
use mod_sheath_bc, only: sheath_current, sheath_norm, sheath_get_lambda, sheath_V_wall_at, &
                         sheath_temp_floor, dsheath_temp_floor_dT
use mod_sheath_diag, only: sheath_diag_add

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)        ! the two nodes containing the boundary nodes
integer, intent(in)   :: i_tor_min   
integer, intent(in)   :: i_tor_max   

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
logical    :: weak_sheath_zj   ! penalty enforcement of the characteristic at the Gauss points
logical    :: wk_here          ! ... and whether it is active at THIS Gauss point
real*8     :: wk_res            ! bounded sheath residual driving the weak term
real*8     :: wk_dfac           ! d(bounded)/d(raw), the factor every Jacobian column carries
real*8     :: wk_cap, wk_den   ! the cap itself and the saturating denominator
real*8     :: dzj_sh, dzj_sat, dzj_x, dzj_d1, dzj_d2, dzj_d3, dzj_d4, dzj_wgt
real*8     :: sh_duf_dTi, sh_duf_dTe
real*8     :: gradu0dotn, gradps0dotn, gradudotn, gradpsidotn
real*8     :: zj_sh, dzj_du, dzj_drho, dzj_dTi, dzj_dTe, zj_sat_g, x_sheath
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

! --- If one of the nodes has a boundary type where natural BCs are applied, apply boundary integral for the full bnd element
diag_sheath_zj = bcs(bnd_type1)%sheath_zj .or. bcs(bnd_type2)%sheath_zj
weak_sheath_zj = ( bcs(bnd_type1)%sheath_zj_weak .or. bcs(bnd_type2)%sheath_zj_weak )   &
                 .and. ( sheath_weak_beta .gt. 0.d0 )
diag_sheath_zj = diag_sheath_zj .or. weak_sheath_zj

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
!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

  ws = wgauss(ms)

  dl   = sqrt(x_s(ms)**2 + y_s(ms)**2) 
  xjac = x_s(ms)*y_t(ms) - x_t(ms)*y_s(ms)
  BigR = x_g(ms)

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
    if (vpar_smoothing) then
      factor = 0.25d0 * ( 1.d0 + tanh( (abs(bdotn) - c_1) / c_2 ) )**2 - c_3
    else
      factor = 1.d0
    endif

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
    if ( apply_natural_bc(var_u) ) then
      call sheath_current(u0, r0_corr, Ti0_sh, Te0_sh, g_bn, normal_sign, Btot, &
                          zj_sh, dzj_du, dzj_drho, dzj_dTi, dzj_dTe, zj_sat_g, x_sheath, &
                          sheath_V_wall_at(BigR))
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
      call sheath_diag_add(bnd_type1, zj_sh, zj0, zj_sat_g, x_sheath, u0, Te0_corr, sh_Bn, &
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
    dzj_sh = 0.d0; dzj_d1 = 0.d0; dzj_d2 = 0.d0; dzj_d3 = 0.d0; dzj_d4 = 0.d0
    dzj_sat = 0.d0; wk_res = 0.d0; wk_dfac = 0.d0
    wk_here = weak_sheath_zj .and. (.not. apply_natural_bc(var_u))

    if ( diag_sheath_zj .and. (.not. apply_natural_bc(var_u)) ) then
      call sheath_current(u0, r0_corr, Ti0_sh, Te0_sh, g_bn, normal_sign, Btot, &
                          dzj_sh, dzj_d1, dzj_d2, dzj_d3, dzj_d4, dzj_sat, dzj_x,      &
                          sheath_V_wall_at(BigR))
      dzj_wgt = 1.d0
      if ( sheath_min_bn .gt. 0.d0 ) &
        dzj_wgt = bdotn**2 / ( bdotn**2 + sheath_min_bn**2 )
      call sheath_diag_add(bnd_type1, dzj_sh, zj0, dzj_sat, dzj_x, u0, Te0_corr, sh_Bn, &
                           ws * dl * BigR * TWOPI / dble(n_plane), dzj_wgt,             &
                           sheath_V_wall_at(BigR), BigR)

      ! --- sheath_current differentiates w.r.t. the FLOORED state, so chain the floors through
      ! --- before these are used as Jacobian entries against the solution variables.
      dzj_d2 = dzj_d2 * dcorr_neg_dens_drho1(r0)
      dzj_d3 = dzj_d3 * dsheath_temp_floor_dT(Ti0)
      dzj_d4 = dzj_d4 * dsheath_temp_floor_dT(Te0)

      ! --- Trust region on the residual. zj_sh - zj is unbounded on the electron branch - at
      ! --- u = 0 the characteristic sits at X = -Lambda, f = 1-exp(Lambda) ~ -19 - so from a
      ! --- restart whose wall is not already near the floating potential the penalty is asked for
      ! --- ~20*j_sat of driving force on the first step. Measured doing exactly that: ePhi/kTe
      ! --- 0.00/0.01/0.03, e-limited 100%, I_sheath -12 kA against I_Ampere +1 kA, blow-up in 4
      ! --- steps. r/(1+|r|/cap) is the IDENTITY for |r| << cap, so the fixed point and the local
      ! --- convergence rate are untouched, and saturates at cap however far off the state starts.
      ! --- Its derivative 1/(1+|r|/cap)^2 decays only algebraically, unlike a tanh clip whose
      ! --- exponential decay would leave the row effectively explicit during the approach.
      wk_res  = dzj_sh - zj0
      wk_dfac = 1.d0
      if ( sheath_weak_rmax .gt. 0.d0 ) then
        wk_cap  = sheath_weak_rmax * max(abs(dzj_sat), 1.d-30)
        wk_den  = 1.d0 + abs(wk_res) / wk_cap
        wk_res  = wk_res  / wk_den
        wk_dfac = 1.d0 / wk_den**2
      endif
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

          ! --- WEAK sheath characteristic: integral(beta * v * (zj_sh - zj)) over the boundary,
          ! --- added to the zj equation instead of replacing its rows. The nodal route satisfies
          ! --- the characteristic at the nodes (measured |zj-zj_sh|/|zj_sat| = 1.8e-2) and misses
          ! --- it by 3.9 at the Gauss points where the currents are integrated - so I_Ampere and
          ! --- I_sheath cannot close no matter how the characteristic or the gates are tuned.
          ! --- Enforcing it at the Gauss points closes that gap by construction, and zj_sh is
          ! --- evaluated there with the LOCAL g_bn, |B| and R, so the incomplete chain rule of
          ! --- the nodal tangential-derivative row (which treats |B| and R as constant along the
          ! --- boundary) does not arise either. Algebraic equation: no theta, no tstep.
          ! --- sheath_ramp_time ramps the penalty in, exactly as it ramps the nodal route. This
          ! --- matters because the characteristic is unbounded on the electron side: starting
          ! --- from u = 0 gives X = -Lambda and f = 1 - exp(Lambda) ~ -19, so the penalty is
          ! --- asked to drag zj to 19*j_sat of electron current on the first step. Observed:
          ! --- e-limited 100%, I_sheath -12 kA against I_Ampere +1 kA, blow-up in four steps.
          ! --- /BigR, matching the zj equation this is added to:
          ! ---   rhs_ij(var_zj) = -( v_x*ps0_x + v_y*ps0_y + v*zj0 ) / BigR * xjac
          ! --- so the penalty carries the same geometric weight as the term it competes with.
          if ( wk_here ) &
            rhs_ij(var_zj) = rhs_ij(var_zj)                                             &
                           + v * sheath_weak_beta * sh_ramp_t * wk_res / BigR * dl

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
            if ( .not. ( apply_natural_bc(i_var) .or.                                  &
                         ( wk_here .and. (i_var .eq. var_zj) ) ) ) cycle
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

                ! --- Jacobian of the weak sheath term. amat = -d(rhs)/dx, and rhs carries
                ! --- +beta*(zj_sh - zj), so the zj column is POSITIVE (it damps) and the state
                ! --- columns carry -d(zj_sh)/dx. Value columns only - psi, never psi_t - which is
                ! --- why this assembles where the zj = Delta*psi surface term above cannot.
                if ( wk_here ) then
                  amat(var_zj,var_zj ) = amat(var_zj,var_zj )                          &
                                       + v * sheath_weak_beta * sh_ramp_t * wk_dfac * psi / BigR * dl
                  amat(var_zj,var_u  ) = amat(var_zj,var_u  )                          &
                                       - v * sheath_weak_beta * sh_ramp_t * wk_dfac * dzj_d1 * psi / BigR * dl
                  amat(var_zj,var_rho) = amat(var_zj,var_rho)                          &
                                       - v * sheath_weak_beta * sh_ramp_t * wk_dfac * dzj_d2 * psi / BigR * dl
                  if ( with_TiTe ) then
                    amat(var_zj,var_Ti) = amat(var_zj,var_Ti)                          &
                                        - v * sheath_weak_beta * sh_ramp_t * wk_dfac * dzj_d3 * psi / BigR * dl
                    amat(var_zj,var_Te) = amat(var_zj,var_Te)                          &
                                        - v * sheath_weak_beta * sh_ramp_t * wk_dfac * dzj_d4 * psi / BigR * dl
                  else
                    amat(var_zj,var_T ) = amat(var_zj,var_T )                          &
                                        - v * sheath_weak_beta * sh_ramp_t * wk_dfac      &
                                            * 0.5d0*(dzj_d3+dzj_d4) * psi / BigR * dl
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

return
end subroutine

end module mod_boundary_matrix_open
