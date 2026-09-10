module mod_boundary_matrix_open
  implicit none
contains

subroutine boundary_matrix_open(vertex, direction, element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, &
                                psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, i_tor_min, i_tor_max, element_id)
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
use mod_floating_transport, only: floating_mach_flux, floating_wall_flux, floating_temperature_slope
use mod_floating_transport_diag, only: transport_diag_wall, weak_mach_diag_sample
use mod_floating_boundary_edges, only: floating_edge_is_exterior
use mod_floating_u, only: floating_u_norm

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)        ! the two nodes containing the boundary nodes
integer, intent(in)   :: i_tor_min   
integer, intent(in)   :: i_tor_max   
integer, intent(in), optional :: element_id

real*8     :: x_g(n_gauss), x_s(n_gauss), x_t(n_gauss), x_ss(n_gauss)
real*8     :: y_g(n_gauss), y_s(n_gauss), y_t(n_gauss), y_ss(n_gauss)

real*8     :: eq_g(n_plane,n_var,n_gauss), eq_s(n_plane,n_var,n_gauss), eq_p(n_plane,n_var,n_gauss)
real*8     :: eq_t(n_plane,n_var,n_gauss), eq_ss(n_plane,n_var,n_gauss)
real*8     :: delta_g(n_plane,n_var,n_gauss), delta_s(n_plane,n_var,n_gauss)

real*8     :: ELM(n_vertex_max*n_var*n_degrees*n_tor,n_vertex_max*n_var*n_degrees*n_tor)
real*8     :: RHS(n_vertex_max*n_var*n_degrees*n_tor)
real*8     :: rhs_ij(n_var), amat(n_var,n_var)

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
logical :: fu_edge, fu_mach, fu_wall
!! WEAK (Galerkin) DRIFT-INCLUSIVE BOHM CONDITION on Vpar. One residual per boundary
!! quadrature point, projected onto the trace test functions, replacing the nodal
!! value+slope rows entirely:
!!
!!    res = (B.n)*Vpar - max( cs*|b.n| - vE.n , 0 )
!!
!! i.e. the parallel normal flow must supply whatever the drift does not, and is never
!! asked to reverse (Bohm is an inequality: if the drift alone already exceeds sonic
!! outflow there is nothing to impose). The row is weighted by d(res)/d(Vpar) = B.n,
!! which is the Galerkin projection onto the space Vpar can control - and which makes
!! the constraint fade smoothly as the field grazes the wall, with NO threshold: the
!! Vpar column carries (B.n)^2, so a tangential point simply loses authority and the
!! natural condition grad(Vpar).n = 0 takes over continuously.
logical :: mw_on
real*8  :: mw_bnu, mw_tgt, mw_res, mw_act, mw_w, mw_bnj, mw_btj, mw_bnuj
real*8 :: fu_ven, fu_bn, fu_vn, fu_orient, fu_mres, fu_mjac(3), fu_ven_trial
!> Outward ExB normal speed entering the SHEATH TRANSMISSION term, and the exact
!! derivative flag of the clip that bounds it. See the block where they are set.
real*8 :: fu_ven_sh, fu_ven_act, fu_ven_open, fu_vtot
real*8 :: fu_particle, fu_heat_i, fu_heat_e, fu_dp(2), fu_dhi(2), fu_dhe(2)
real*8 :: fu_slope_i, fu_slope_e, fu_knee, fu_area
real*8 :: fu_a,fu_ct,fu_cv,fu_qjac,fu_res
real*8 :: fu_normal_grad_trial,fu_transverse_trial,fu_visc_col
real*8 :: fu_cs(3)
integer :: fu_side,fu_normal_dof,fu_normal_col

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

R_mid = sum(nodes(1:2)%x(1,1,1)) / 2.d0     ! mid point on boundary (approx.)
Z_mid = sum(nodes(1:2)%x(1,1,2)) / 2.d0
R_cnt = sum(nodes(1:4)%x(1,1,1)) / 4.d0     ! center point within element (approx.)
Z_cnt = sum(nodes(1:4)%x(1,1,2)) / 4.d0

normal_direction = (/R_mid - R_cnt, Z_mid - Z_cnt /) / norm2((/R_mid - R_cnt, Z_mid - Z_cnt /))

apply_natural_bc(:) = .false.

bnd_type1 = nodes(1)%boundary 
bnd_type2 = nodes(2)%boundary 
! Both endpoints must be covered: a floating type-3 corner does not turn type 2
! into a material wall. Exterior topology is checked by the shared caller.
fu_edge = bcs(bnd_type1)%floating_u .and. bcs(bnd_type2)%floating_u
fu_side=0
select case(vertex(1)*vertex(2))
case(2)
  fu_side=1
case(6)
  fu_side=2
case(12)
  fu_side=3
case(4)
  fu_side=4
end select
! The sheath ExB energy flux below needs this too, and construct_matrix_mod now
! builds the topology whenever any boundary type carries the floating potential,
! so the refinement is no longer tied to the opt-in transport experiments.
if (fu_edge) then
  if (.not. present(element_id)) error stop 'floating transport requires the exterior edge identity'
  fu_edge=floating_edge_is_exterior(element_id,fu_side)
endif
if (floating_u_transport_diag) call floating_u_norm(fu_a,fu_ct,fu_cv)
! --- mach1_weak REPLACES every other Mach route. Precedence is enforced here rather
! --- than trusted to the namelist so the three can never be assembled together.
mw_on = mach1_weak .and. with_vpar .and. (bcs(bnd_type1)%mach1 .or. bcs(bnd_type2)%mach1)
fu_mach = floating_u_mach_flux .and. fu_edge .and. (.not. mach1_weak) .and. &
          (bcs(bnd_type1)%mach1 .or. bcs(bnd_type2)%mach1)
fu_wall = floating_u_wall_flux .and. fu_edge

! --- If one of the nodes has a boundary type where natural BCs are applied, apply boundary integral for the full bnd element
do i_var=1, n_var
  if ( (i_var==var_rho ) .and. (bcs(bnd_type1)%natural%rho  .or. bcs(bnd_type2)%natural%rho ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_T   ) .and. (bcs(bnd_type1)%natural%T    .or. bcs(bnd_type2)%natural%T   ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_Ti  ) .and. (bcs(bnd_type1)%natural%Ti   .or. bcs(bnd_type2)%natural%Ti  ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_Te  ) .and. (bcs(bnd_type1)%natural%Te   .or. bcs(bnd_type2)%natural%Te  ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_rhon) .and. (bcs(bnd_type1)%natural%rhon .or. bcs(bnd_type2)%natural%rhon))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_vpar) .and. (bcs(bnd_type1)%natural%vpar .or. bcs(bnd_type2)%natural%vpar))  apply_natural_bc(i_var)=.true.
enddo

if (fu_mach .or. mw_on) apply_natural_bc(var_vpar)=.true.

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
    r0_corr  = corr_neg_dens(r0)

    if (with_TiTe) then
      cs0    = sqrt(gamma*(Ti0_corr+Te0_corr))
    else
      cs0    = sqrt(gamma*T0_corr)
    endif

    Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR

    BB2 = Btot**2

    ! The sign is obtained from the actual outward normal, not a type/side label.
    fu_orient = sign(1.d0,y_s(ms)*normal(1)-x_s(ms)*normal(2))
    fu_ven = -fu_orient*BigR*eq_s(mp,var_u,ms)/dl
    fu_bn = (ps0_y*normal(1)-ps0_x*normal(2))/BigR
    fu_vn = fu_bn*Vpar0+fu_ven
    ! --- SHEATH TRANSMISSION MUST BE CHARGED ON THE FLUX THAT ACTUALLY ARRIVES.
    ! ---
    ! --- The volume energy equation convects with the FULL velocity (the ExB terms
    ! --- Ti0_s*u0_t-Ti0_t*u0_s and 2*GAMMA*R*u0_y sit beside the parallel ones in
    ! --- mod_elt_matrix_fft) and its strong form integrates to the surface flux. The
    ! --- boundary term supplies the sheath transmission in EXCESS of that convection,
    ! --- hence (gamma_sheath-1). It was computed from vpar0*ps0_s*normal_sign3, which
    ! --- is identically (fu_bn*Vpar0)*BigR*dl - the PARALLEL normal flux alone - so
    ! --- the two halves of one sheath transmission used two different flows. With
    ! --- dirichlet u the trace of u along the wall is constant, fu_ven vanishes and
    ! --- the inconsistency is invisible; bcs%floating_u ties u to Te and it becomes
    ! --- the dominant term on a grazing wall.
    ! ---
    ! --- The collection is therefore the TOTAL normal flow,
    ! ---
    ! ---     fu_vtot = fu_bn*Vpar0 + clip(fu_ven, +-2*cs*|bn|),
    ! ---
    ! --- floored at zero because a material wall absorbs and never emits, and the
    ! --- term added below is the difference from the parallel-only expression that
    ! --- is already there, so nothing without a wall-tangential potential gradient
    ! --- changes at all.
    ! ---
    ! --- WHY BOTH SIGNS OF fu_ven, not max(fu_ven,0) as first written: that form
    ! --- silently assumed Vpar*Bn = cs*|bn|, which holds only with the drift term
    ! --- disabled. Once the Mach row carries the drift supplement, Vpar is
    ! --- supersonic by exactly the amount needed to cancel an INWARD ExB, and
    ! --- charging "parallel + outward-ExB-only" bills the sheath for the supplement
    ! --- while the drift it compensates carries the flux back out. Measured on the
    ! --- production case: supplement 1.45*cs on boundary type 1 gave a 2.45x
    ! --- over-charge, driving that type alone onto the temperature floor while
    ! --- types 4 and 9 (supplement identically zero) stayed at 30 eV. Using the
    ! --- total flow is what SOLPS's linked BCCON=14/BCENE,I=15 set does by sharing
    ! --- one U_out with BCMOM=13.
    ! ---
    ! --- It gives exactly the Bohm collection cs*|bn| only where the Mach row's
    ! --- incidence floor is INACTIVE. That floor multiplies the supplement by
    ! --- f = min(1,|bn|/s0), so with an inward drift the total is cs*|bn| +
    ! --- (1-f)*vE.n: at grazing incidence part of the drift is deliberately left
    ! --- uncompensated, and this term reports that honestly rather than pretending
    ! --- it was cancelled. The c_angle minimum-flux term is additional to both.
    ! ---
    ! --- The sheath energy flux uses THE SAME total normal flow the Mach1 row does:
    ! --- fu_bn*Vpar0 + fu_ven, unclipped. There used to be a SOLPS bound at
    ! --- 2*cs*|b.n| here, whose only justification was to match the Mach row's
    ! --- 2*cs/Btot clip on its own supplement. That clip is gone, so keeping this one
    ! --- would leave the momentum and energy conditions bounding the same drift
    ! --- differently - the exact inconsistency the bound was introduced to avoid.
    ! ---
    ! --- Branches, each with an exact derivative (fu_ven_open selects an open wall,
    ! --- fu_ven_act the ExB column):
    ! ---   fu_vtot <= 0 : collection 0    d/du 0      d/dVpar -fu_bn
    ! ---   fu_vtot >  0 : sh = fu_ven     d/du trial  d/dVpar 0
    fu_ven_sh   = 0.d0
    fu_ven_act  = 0.d0
    fu_ven_open = 0.d0
    if (fu_edge) then
      fu_vtot = fu_bn*Vpar0 + fu_ven
      if (fu_vtot > 0.d0) then
        fu_ven_open = 1.d0
        fu_ven_sh   = fu_ven
        fu_ven_act  = 1.d0
      else
        ! Wall closed: total inflow, collect nothing rather than emit.
        fu_ven_sh = -fu_bn*Vpar0
      endif
    endif
    if (floating_u_transport_diag .and. fu_edge .and. mp==1) then
      fu_qjac=abs(xjac)/(dl*norm2((/x_t(ms),y_t(ms)/)))
      fu_res=eq_g(mp,var_u,ms)-fu_ct*Te0-fu_cv*sheath_V_wall
      call transport_diag_wall(element%vertex(1:4),fu_side,ms,BigR,y_g(ms),r0,Ti0,Te0, &
          fu_ven,fu_bn*Vpar0,cs0,fu_bn/Btot,fu_qjac,fu_res,(/bnd_type1,bnd_type2/))
    endif
    if (fu_mach) call floating_mach_flux(Vpar0,fu_ven,fu_bn,Btot,cs0,vpar_smoothing, &
                                       vpar_smoothing_coef,fu_mres,fu_mjac)

    ! --- WEAK BOHM RESIDUAL. fu_bn is B_pol.n and carries a factor |B|, which is why
    ! --- fu_bn*Vpar0 is already a velocity (Vpar is v/|B|); |b.n| therefore needs the
    ! --- explicit /Btot. No clip, no floor, no angle threshold enters.
    mw_bnu = abs(fu_bn) / Btot
    mw_tgt = cs0*mw_bnu - fu_ven
    mw_act = 1.d0
    if ( mw_tgt .le. 0.d0 ) then
      mw_tgt = 0.d0                  ! drift alone already sonic or more: impose nothing
      mw_act = 0.d0
    endif
    mw_res = fu_bn*Vpar0 - mw_tgt
    mw_w   = fu_bn                   ! = d(res)/d(Vpar)
    ! --- Sampled only under floating_u_diag, which is also the flag the reset and the
    ! --- print are gated on. Otherwise the running max and the point counts would
    ! --- accumulate across every timestep with nothing ever clearing them.
    if (mw_on .and. floating_u_diag) &
      call weak_mach_diag_sample(bnd_type1,mw_res,cs0*mw_bnu,mw_tgt,mw_bnu,mw_act)
    if (fu_wall) then
      call floating_wall_flux(fu_vn,cs0,c_angle,gamma_sheath_i,fu_particle,fu_heat_i,fu_dp,fu_dhi)
      call floating_wall_flux(fu_vn,cs0,c_angle,gamma_sheath_e,fu_particle,fu_heat_e,fu_dp,fu_dhe)
      if (.not. with_TiTe) &
        call floating_wall_flux(fu_vn,cs0,c_angle,gamma_sheath,fu_particle,fu_heat_i,fu_dp,fu_dhi)
    endif
    ! --- UNCONDITIONAL. cs0 above is built from corr_neg_temp1(T), so d(cs)/dT carries
    ! --- corr_neg_temp1'(T) - which is exactly what floating_temperature_slope returns.
    ! --- It used to be formed only for the opt-in wall/mach experiments, so every
    ! --- production run differentiated a corrected sound speed as if it were
    ! --- uncorrected. Below the corr_neg knee that overstates d(cs)/dT badly (a factor
    ! --- 3 at 0.65 eV with T_min_neg = 3e-5, and exponentially worse below), and the
    ! --- divertor wall in this campaign sits exactly there - so the implicit solve
    ! --- believed the sheath sink relaxes with T far more strongly than it does.
    fu_knee = T_min_neg
    if (fu_knee < 0.d0) fu_knee=T_1
    fu_slope_i = floating_temperature_slope(Ti0,fu_knee,corr_neg_temp_coef)
    fu_slope_e = floating_temperature_slope(Te0,fu_knee,corr_neg_temp_coef)
    if (.not. with_TiTe) fu_slope_i=floating_temperature_slope(T0,fu_knee,corr_neg_temp_coef)

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
    if (mach_one_bnd_integral .and. .not. mach1_weak) factor_cs_bnd_integral = 1.d0

    do i=1,2                ! loop over nodes

      do j=1,2              ! loop over basis functions

        j2 = direction(j)
        element_size_ij = element%size(vertex(i),j2)

        do im=i_tor_min, i_tor_max

          v   =  H1(i,j,ms) * element_size_ij * HZ(im,mp)         ! test function

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
                                - v * (gamma_sheath_i-1.d0) * r0 * Ti0 * fu_ven_sh * BigR * dl        * tstep &
                                - v * (gamma_sheath_i-1.d0) * r0 * Ti0 * cs0    * BigR * dl * c_angle * tstep & 
                                - v * (GAMMA - 1.d0) * vpar0 * visco_par_heating * gradvpar0dotn * BigR * dl  * tstep  

              rhs_ij(var_Te)  = - v * (gamma_sheath_e-1.d0) * r0 * Te0 * vpar0 * ps0_s * normal_sign3 * tstep &
                                - v * (gamma_sheath_e-1.d0) * r0 * Te0 * fu_ven_sh * BigR * dl        * tstep &
                                - v * (gamma_sheath_e-1.d0) * r0 * Te0 * cs0  * BigR * dl * c_angle   * tstep  
            else
              rhs_ij(var_T)   = - v * (gamma_sheath  -1.d0) * r0 * T0  * vpar0 * ps0_s * normal_sign3 * tstep &
                                - v * (gamma_sheath  -1.d0) * r0 * T0  * fu_ven_sh * BigR * dl        * tstep &
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
          if (fu_mach) rhs_ij(var_vpar) = -v*dl*Zbig*fu_mres
          if (mw_on)   rhs_ij(var_vpar) = -v*BigR*dl*Zbig*mw_w*mw_res
          if (fu_wall) then
            fu_area = v*BigR*dl*tstep
            ! REPLACE the legacy diffusive boundary corrections, not the volume
            ! advection. Charged-particle reflection is zero for this experiment.
            rhs_ij(var_rho) = -fu_area*r0*fu_particle
            if (with_TiTe) then
              rhs_ij(var_Ti) = -fu_area*r0*Ti0*fu_heat_i &
                  -fu_area*(GAMMA-1.d0)*Vpar0*visco_par_heating*gradvpar0dotn
              rhs_ij(var_Te) = -fu_area*r0*Te0*fu_heat_e
            else
              rhs_ij(var_T) = -fu_area*r0*T0*fu_heat_i &
                  -fu_area*(GAMMA-1.d0)*Vpar0*visco_par_heating*gradvpar0dotn
            endif
          endif
          index_ij = n_tor_local*n_var*n_degrees*(vertex(i)-1) + n_tor_local * n_var * (j2-1) + im - i_tor_min +1  ! index in the ELM matrix

          do i_var = 1, n_var
            if ( .not. apply_natural_bc(i_var) ) cycle
            if ((fu_mach .or. mw_on) .and. i_var==var_vpar) then
              if (.not. bcs(nodes(i)%boundary)%mach1) cycle
            endif
            RHS(index_ij+(i_var-1)*(n_tor_local)) = RHS(index_ij+(i_var-1)*(n_tor_local)) + rhs_ij(i_var) * ws
          enddo


          do k=1,2                                                          ! loop over nodes

            do l=1,2                                                        ! loop over basis functions

              l2 = direction(l)
              element_size_kl = element%size(vertex(k),l2)

              l3 = direction_perp(j)
              element_size_perp = - element%size(vertex(k),direction_perp(1)) * 3.d0

              do in = i_tor_min, i_tor_max                                              ! loop over toroidal harmonics

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

                ! --- d(cs)/dT of the CORRECTED sound speed: the corr_neg_temp1 slope
                ! --- belongs here, not only in the opt-in fu_cs below.
                cs_T   = gamma * T  / (2.d0 * cs0) * fu_slope_i
                cs_Ti  = gamma * Ti / (2.d0 * cs0) * fu_slope_i
                cs_Te  = gamma * Te / (2.d0 * cs0) * fu_slope_e
                if (fu_mach .or. fu_wall) then
                  fu_cs=(/cs_T,cs_Ti,cs_Te/)   ! slope now already in cs_T/cs_Ti/cs_Te
                endif
                fu_ven_trial = -fu_orient*BigR*psi_s/dl
                ! The true trace basis has zero t derivative here. Independent
                ! normal/mixed DOFs supply eq_t; do not substitute a scaled trace.
                fu_normal_grad_trial=(y_t(ms)*normal(1)-x_t(ms)*normal(2))*psi_s/xjac

                ! --- Most of natural BCs need vpar
                if (with_vpar) then

                  ! --- Density reflection and minimum particle flux (c_angle)
                  amat(var_rho,var_psi)   = - v * density_reflection * r0  * vpar0 * psi_s * normal_sign3 * theta * tstep 
                  amat(var_rho,var_rho)   = - v * density_reflection * rho * vpar0 * ps0_s * normal_sign3 * theta * tstep &
                                            + v                      * rho * cs0   * BigR * dl * c_angle  * theta * tstep 
                  ! --- ZERO-SUM SHEATH STABILISER on the particle BC (SOLPS
                  ! --- b2stbc_stab_coeff_sheath_ni). Jacobian only, as for Ti/Te.
                  ! --- NOTE the particle sink here is only the c_angle floor and the
                  ! --- reflected fraction: the bulk of the wall particle loss is carried
                  ! --- by the strong-form volume advection, which this cannot damp. So
                  ! --- this is the weakest of the three by construction.
                  amat(var_rho,var_rho)   = amat(var_rho,var_rho)                                             &
                                          + v * stab_coeff_sheath_ni * rho * ( vpar0 * ps0_s * normal_sign3   &
                                              + fu_ven_sh * BigR * dl + cs0 * BigR * dl * c_angle ) * theta * tstep
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
                                            + v * (gamma_sheath_i-1.d0) * rho * Ti0 * fu_ven_sh * BigR * dl        * theta * tstep &
                                            + v * (gamma_sheath_i-1.d0) * rho * Ti0 * cs0   * BigR  * dl * c_angle * theta * tstep 
                    amat(var_Ti,var_Ti)   = + v * (gamma_sheath_i-1.d0) * r0  * Ti  * vpar0 * ps0_s * normal_sign3 * theta * tstep & 
                                            + v * (gamma_sheath_i-1.d0) * r0  * Ti  * fu_ven_sh * BigR * dl        * theta * tstep &
                                            + v * (gamma_sheath_i-1.d0) * r0  * Ti  * cs0   * BigR  * dl * c_angle * theta * tstep &
                                            + v * (gamma_sheath_i-1.d0) * r0  * Ti0 * cs_Ti * BigR  * dl * c_angle * theta * tstep
                    ! --- Exact derivative of the outward ExB flux. fu_ven_act is 1 on the
                    ! --- open branch and 0 on the closed one, so the column is exact in both.
                    ! --- With no clip there is no bound left to differentiate, hence no
                    ! --- temperature column from this term.
                    amat(var_Ti,var_u)    = + v * (gamma_sheath_i-1.d0) * r0  * Ti0 * fu_ven_act * fu_ven_trial * BigR * dl * theta * tstep
                    ! --- ZERO-SUM SHEATH STABILISER (SOLPS b2stbc_stab_coeff_sheath_*).
                    ! --- alpha multiplies the SAME sheath flux prefactor as the transmission
                    ! --- coefficient and enters the JACOBIAN ONLY; the residual keeps the
                    ! --- physical gamma. Net source therefore gains alpha*prefactor*(X_old -
                    ! --- X_new), which vanishes identically at X_new = X_old - the steady
                    ! --- state is unchanged and only the approach to it is damped.
                    amat(var_Ti,var_Ti)   = amat(var_Ti,var_Ti)                                       &
                                          + v * stab_coeff_sheath_ti * r0 * Ti * ( vpar0 * ps0_s * normal_sign3         &
                                              + fu_ven_sh * BigR * dl + cs0 * BigR * dl * c_angle ) * theta * tstep
                    ! --- was the clip bound's cross-temperature column; the bound is gone.
                    ! --- Assigned, not deleted: this is a plain "=" slot, and an unwritten
                    ! --- entry would keep whatever value was already in it.
                    amat(var_Ti,var_Te)   = 0.d0

                    amat(var_Te,var_psi)  = + v * (gamma_sheath_e-1.d0) * r0  * Te0 * vpar0 * psi_s * normal_sign3 * theta * tstep 
                    amat(var_Te,var_rho)  = + v * (gamma_sheath_e-1.d0) * rho * Te0 * vpar0 * ps0_s * normal_sign3 * theta * tstep & 
                                            + v * (gamma_sheath_e-1.d0) * rho * Te0 * fu_ven_sh * BigR * dl        * theta * tstep &
                                            + v * (gamma_sheath_e-1.d0) * rho * Te0 * cs0   * BigR  * dl * c_angle * theta * tstep 
                    amat(var_Te,var_Te)   = + v * (gamma_sheath_e-1.d0) * r0  * Te  * vpar0 * ps0_s * normal_sign3 * theta * tstep &
                                            + v * (gamma_sheath_e-1.d0) * r0  * Te  * fu_ven_sh * BigR * dl        * theta * tstep &
                                            + v * (gamma_sheath_e-1.d0) * r0  * Te  * cs0   * BigR  * dl * c_angle * theta * tstep &
                                            + v * (gamma_sheath_e-1.d0) * r0  * Te0 * cs_Te * BigR  * dl * c_angle * theta * tstep
                    amat(var_Te,var_u)    = + v * (gamma_sheath_e-1.d0) * r0  * Te0 * fu_ven_act * fu_ven_trial * BigR * dl * theta * tstep
                    ! --- ZERO-SUM SHEATH STABILISER (SOLPS b2stbc_stab_coeff_sheath_*).
                    ! --- alpha multiplies the SAME sheath flux prefactor as the transmission
                    ! --- coefficient and enters the JACOBIAN ONLY; the residual keeps the
                    ! --- physical gamma. Net source therefore gains alpha*prefactor*(X_old -
                    ! --- X_new), which vanishes identically at X_new = X_old - the steady
                    ! --- state is unchanged and only the approach to it is damped.
                    amat(var_Te,var_Te)   = amat(var_Te,var_Te)                                       &
                                          + v * stab_coeff_sheath_te * r0 * Te * ( vpar0 * ps0_s * normal_sign3         &
                                              + fu_ven_sh * BigR * dl + cs0 * BigR * dl * c_angle ) * theta * tstep
                    amat(var_Te,var_Ti)   = 0.d0

                    ! --- The closed-wall term below must be folded into THESE assignments,
                    ! --- not accumulated earlier: these are plain "=" and would overwrite it,
                    ! --- leaving a residual that says the collection is zero and a Jacobian
                    ! --- that still responds as if the parallel collection were present.
                    amat(var_Ti,var_vpar) = + v * (gamma_sheath_i-1.d0) * r0  * Ti0 * vpar  * ps0_s * normal_sign3 * theta * tstep &
                                            - v * (gamma_sheath_i-1.d0) * r0  * Ti0 * (1.d0-fu_ven_open) * fu_bn * vpar * BigR * dl * theta * tstep &
                                            + v * (GAMMA - 1.d0) * vpar * visco_par_heating * gradvpar0dotn * BigR * dl    * theta * tstep &
                                            + v * (GAMMA - 1.d0) * vpar0 * visco_par_heating * gradvpardotn * BigR * dl    * theta * tstep
                    amat(var_Te,var_vpar) = + v * (gamma_sheath_e-1.d0) * r0  * Te0 * vpar  * ps0_s * normal_sign3 * theta * tstep &
                                            - v * (gamma_sheath_e-1.d0) * r0  * Te0 * (1.d0-fu_ven_open) * fu_bn * vpar * BigR * dl * theta * tstep 
                  else
                    amat(var_T,var_psi)   = + v * (gamma_sheath  -1.d0) * r0  *  T0 * vpar0 * psi_s * normal_sign3 * theta * tstep 
                    amat(var_T,var_rho)   = + v * (gamma_sheath  -1.d0) * rho *  T0 * vpar0 * ps0_s * normal_sign3 * theta * tstep &
                                            + v * (gamma_sheath  -1.d0) * rho *  T0 * fu_ven_sh * BigR * dl        * theta * tstep &
                                            + v * (gamma_sheath  -1.d0) * rho *  T0 * cs0   * BigR  * dl * c_angle * theta * tstep 
                    amat(var_T,var_T)     = + v * (gamma_sheath  -1.d0) * r0  *  T  * vpar0 * ps0_s * normal_sign3 * theta * tstep &
                                            + v * (gamma_sheath  -1.d0) * r0  *  T  * fu_ven_sh * BigR * dl        * theta * tstep &
                                            + v * (gamma_sheath  -1.d0) * r0  *  T  * cs0   * BigR  * dl * c_angle * theta * tstep &
                                            + v * (gamma_sheath  -1.d0) * r0  *  T0 * cs_T  * BigR  * dl * c_angle * theta * tstep
                    amat(var_T,var_u)     = + v * (gamma_sheath  -1.d0) * r0  *  T0 * fu_ven_act * fu_ven_trial * BigR * dl * theta * tstep
                    ! --- ZERO-SUM SHEATH STABILISER (SOLPS b2stbc_stab_coeff_sheath_*).
                    ! --- alpha multiplies the SAME sheath flux prefactor as the transmission
                    ! --- coefficient and enters the JACOBIAN ONLY; the residual keeps the
                    ! --- physical gamma. Net source therefore gains alpha*prefactor*(X_old -
                    ! --- X_new), which vanishes identically at X_new = X_old - the steady
                    ! --- state is unchanged and only the approach to it is damped.
                    amat(var_T,var_T)   = amat(var_T,var_T)                                       &
                                          + v * stab_coeff_sheath_ti * r0 * T * ( vpar0 * ps0_s * normal_sign3         &
                                              + fu_ven_sh * BigR * dl + cs0 * BigR * dl * c_angle ) * theta * tstep

                    amat(var_T,var_vpar)  = + v * (gamma_sheath  -1.d0) * r0  * T0  * vpar  * ps0_s * normal_sign3 * theta * tstep & 
                                            - v * (gamma_sheath  -1.d0) * r0  *  T0 * (1.d0-fu_ven_open) * fu_bn * vpar * BigR * dl * theta * tstep &
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
                ! Magnetic geometry is lagged for these experiments. All vpar,
                ! u, rho and temperature dependencies are differentiated.
                if (mw_on) then
                  ! --- WEAK BOHM JACOBIAN. Every column of the residual is here; the
                  ! --- row replaces whatever the natural vpar terms put in it.
                  ! ---   res  = (B.n)*Vpar - max(cs*|b.n| - vE.n, 0)
                  ! ---   d/dVpar = B.n            d/du = +act*d(vE.n)/du
                  ! ---   d/dT    = -act*|b.n|*d(cs)/dT
                  ! ---   d/dpsi  via B.n, using ps0_s*normal_sign3 == fu_bn*BigR*dl so
                  ! ---           the trial contribution to B.n is psi_s*normal_sign3/(BigR*dl)
                  ! --- The Galerkin WEIGHT mw_w is left lagged. That omission is safe in a
                  ! --- way the nodal row's missing u column was NOT: d(w)/dx multiplies the
                  ! --- residual itself, so it vanishes as the constraint is met, whereas a
                  ! --- missing d(res)/dx term does not vanish and accumulates every step.
                  amat(var_vpar,:)=0.d0
                  amat(var_vpar,var_vpar)= v*BigR*dl*Zbig*mw_w * mw_w * vpar
                  amat(var_vpar,var_u)   = v*BigR*dl*Zbig*mw_w * mw_act * fu_ven_trial
                  ! --- psi column. The row that reaches the RHS is w*res with w = B.n,
                  ! --- and BOTH factors depend on psi, so the derivative is
                  ! ---     w*d(res)/dpsi + res*d(w)/dpsi
                  ! --- The second term is NOT negligible: it is multiplied by the
                  ! --- residual, which is only small once the constraint is met. During
                  ! --- a transient it dominated the column by a factor ~2.7 in the unit
                  ! --- test. Omitting a term because "it vanishes at the solution" is
                  ! --- precisely the reasoning that left the nodal slope row without its
                  ! --- u columns, so it is carried here.
                  ! --- d(B.n)/ddof uses ps0_s*normal_sign3 == fu_bn*BigR*dl, and
                  ! --- d(Btot)/ddof the trial R/Z derivatives, so |b.n| is exact too.
                  mw_bnj  = psi_s*normal_sign3/(BigR*dl)
                  mw_btj  = ( ps0_x*rho_x + ps0_y*rho_y ) / ( BigR**2 * Btot )
                  mw_bnuj = sign(1.d0,fu_bn)*mw_bnj/Btot - abs(fu_bn)*mw_btj/Btot**2
                  amat(var_vpar,var_psi) = v*BigR*dl*Zbig                                   &
                        * (   mw_w * ( Vpar0*mw_bnj - mw_act*cs0*mw_bnuj )                  &
                            + mw_res * mw_bnj )
                  if (with_TiTe) then
                    amat(var_vpar,var_Ti)= -v*BigR*dl*Zbig*mw_w * mw_act * cs_Ti * mw_bnu
                    amat(var_vpar,var_Te)= -v*BigR*dl*Zbig*mw_w * mw_act * cs_Te * mw_bnu
                  else
                    amat(var_vpar,var_T) = -v*BigR*dl*Zbig*mw_w * mw_act * cs_T  * mw_bnu
                  endif
                endif
                if (fu_mach) then
                  amat(var_vpar,:)=0.d0
                  amat(var_vpar,var_vpar)=v*dl*Zbig*fu_mjac(1)*vpar
                  amat(var_vpar,var_u)=v*dl*Zbig*fu_mjac(2)*fu_ven_trial
                  if (with_TiTe) then
                    amat(var_vpar,var_Ti)=v*dl*Zbig*fu_mjac(3)*fu_cs(2)
                    amat(var_vpar,var_Te)=v*dl*Zbig*fu_mjac(3)*fu_cs(3)
                  else
                    amat(var_vpar,var_T)=v*dl*Zbig*fu_mjac(3)*fu_cs(1)
                  endif
                endif
                if (fu_wall) then
                  fu_area=v*BigR*dl*tstep*theta
                  amat(var_rho,:)=0.d0
                  amat(var_rho,var_rho)=fu_area*rho*fu_particle
                  amat(var_rho,var_u)=fu_area*r0*fu_dp(1)*fu_ven_trial
                  amat(var_rho,var_vpar)=fu_area*r0*fu_dp(1)*fu_bn*vpar
                  if (with_TiTe) then
                    amat(var_rho,var_Ti)=fu_area*r0*fu_dp(2)*fu_cs(2)
                    amat(var_rho,var_Te)=fu_area*r0*fu_dp(2)*fu_cs(3)
                    amat(var_Ti,:)=0.d0
                    amat(var_Te,:)=0.d0
                    amat(var_Ti,var_rho)=fu_area*rho*Ti0*fu_heat_i
                    amat(var_Te,var_rho)=fu_area*rho*Te0*fu_heat_e
                    amat(var_Ti,var_u)=fu_area*r0*Ti0*fu_dhi(1)*fu_ven_trial
                    amat(var_Te,var_u)=fu_area*r0*Te0*fu_dhe(1)*fu_ven_trial
                    amat(var_Ti,var_vpar)=fu_area*r0*Ti0*fu_dhi(1)*fu_bn*vpar &
                        +fu_area*(GAMMA-1.d0)*visco_par_heating*(vpar*gradvpar0dotn+Vpar0*fu_normal_grad_trial)
                    amat(var_Te,var_vpar)=fu_area*r0*Te0*fu_dhe(1)*fu_bn*vpar
                    amat(var_Ti,var_Ti)=fu_area*r0*(Ti*fu_heat_i+Ti0*fu_dhi(2)*fu_cs(2))
                    amat(var_Ti,var_Te)=fu_area*r0*Ti0*fu_dhi(2)*fu_cs(3)
                    amat(var_Te,var_Te)=fu_area*r0*(Te*fu_heat_e+Te0*fu_dhe(2)*fu_cs(3))
                    amat(var_Te,var_Ti)=fu_area*r0*Te0*fu_dhe(2)*fu_cs(2)
                  else
                    amat(var_rho,var_T)=fu_area*r0*fu_dp(2)*fu_cs(1)
                    amat(var_T,:)=0.d0
                    amat(var_T,var_rho)=fu_area*rho*T0*fu_heat_i
                    amat(var_T,var_u)=fu_area*r0*T0*fu_dhi(1)*fu_ven_trial
                    amat(var_T,var_vpar)=fu_area*r0*T0*fu_dhi(1)*fu_bn*vpar &
                        +fu_area*(GAMMA-1.d0)*visco_par_heating*(vpar*gradvpar0dotn+Vpar0*fu_normal_grad_trial)
                    amat(var_T,var_T)=fu_area*r0*(T*fu_heat_i+T0*fu_dhi(2)*fu_cs(1))
                  endif
                endif
                index_kl = n_tor_local*n_var*n_degrees*(vertex(k)-1) + n_tor_local * n_var * (l2-1) + in - i_tor_min +1  ! index in the ELM matrix
                if (fu_wall) then
                  ! Complete the viscous heat-flux derivative on independent
                  ! normal/mixed Vpar DOFs. These are not trace trial functions.
                  fu_normal_dof=direction_perp(l)
                  fu_transverse_trial=-element%size(vertex(k),direction_perp(1))*3.d0
                  if (vertex(1)*vertex(2)==2) fu_transverse_trial=-fu_transverse_trial
                  fu_transverse_trial=fu_transverse_trial*H1(k,l,ms)*element_size_kl*HZ(in,mp)
                  fu_visc_col=fu_area*(GAMMA-1.d0)*visco_par_heating*Vpar0*fu_transverse_trial &
                       *(-y_s(ms)*normal(1)+x_s(ms)*normal(2))/xjac
                  fu_normal_col=n_tor_local*n_var*n_degrees*(vertex(k)-1) &
                       +n_tor_local*n_var*(fu_normal_dof-1)+in-i_tor_min+1+(var_vpar-1)*n_tor_local
                  i_var=var_T
                  if (with_TiTe) i_var=var_Ti
                  if (apply_natural_bc(i_var)) ELM(index_ij+(i_var-1)*n_tor_local,fu_normal_col) = &
                       ELM(index_ij+(i_var-1)*n_tor_local,fu_normal_col)+fu_visc_col*ws
                endif

                ! --- Add contributions to ELM matrix                 
                do k_var = 1, n_var
                  do i_var = 1, n_var

                    if ( .not. apply_natural_bc(i_var) ) cycle
                    if ((fu_mach .or. mw_on) .and. i_var==var_vpar) then
                      if (.not. bcs(nodes(i)%boundary)%mach1) cycle
                    endif

                    ELM(index_ij+(i_var-1)*(n_tor_local),index_kl+(k_var-1)*(n_tor_local)) = &
                    ELM(index_ij+(i_var-1)*(n_tor_local),index_kl+(k_var-1)*(n_tor_local))   &
                      + amat(i_var,k_var) * ws

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
