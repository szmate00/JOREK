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
use mod_wall_diag, only: wall_diag_add, wall_diag_sheath_add
use mod_floating_u, only: floating_u_norm, sheath_j_norm

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
real*8     :: wd_orient, wd_vEn, wd_jr                           ! wall diagnostics: outward ExB normal speed, zj/j_sat
logical    :: sf_on                                              ! sheath-set wall flux on this edge (floating_u + mach1 types)
real*8     :: sf_orient, sf_vEn, sf_Bn, sf_vn, sf_fl            ! outward ExB / parallel / total normal speed, its floor
real*8     :: fx_n, fx_v, fx_p, fx_u, fx_T                      ! normal-flow measure of the sheath fluxes and its columns
real*8     :: ex_n, ex_v, ex_p, ex_u, ex_T                      ! excess sink of the sheath-set wall flux and its columns
logical    :: sj_on, sj_here                                     ! sheath current row on this edge / at this Gauss point
real*8     :: sj_an, sj_csat, sj_CT, sj_CV, sj_jsat, sj_vfl      ! its normalisation, j_sat and the Bohm parallel flow it scales with
real*8     :: sc_x, sc_ex, sc_f, sc_dfdu, sc_dfdTe, sc_res, sc_w, sc_Tc, sc_dTc, sc_drc   ! the characteristic and its columns
logical    :: xpoint2
integer    :: n_tor_local 
logical    :: apply_natural_bc(0:n_var)

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

! --- Sheath-set wall flux on the total normal flow: on edges whose both endpoints carry a wall potential that varies
! --- along the wall (bcs%floating_u or bcs%sheath_j) and the Mach-1 row, so that the flux and the row describe the
! --- same wall
sf_on = with_vpar .and. ( bcs(bnd_type1)%floating_u .or. ( bcs(bnd_type1)%sheath_j .and. sheath_flux_on_sheath_j ) ) &
                  .and. ( bcs(bnd_type2)%floating_u .or. ( bcs(bnd_type2)%sheath_j .and. sheath_flux_on_sheath_j ) ) &
                  .and. bcs(bnd_type1)%mach1 .and. bcs(bnd_type2)%mach1

! --- If one of the nodes has a boundary type where natural BCs are applied, apply boundary integral for the full bnd element
do i_var=1, n_var
  if ( (i_var==var_rho ) .and. (bcs(bnd_type1)%natural%rho  .or. bcs(bnd_type2)%natural%rho ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_T   ) .and. (bcs(bnd_type1)%natural%T    .or. bcs(bnd_type2)%natural%T   ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_Ti  ) .and. (bcs(bnd_type1)%natural%Ti   .or. bcs(bnd_type2)%natural%Ti  ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_Te  ) .and. (bcs(bnd_type1)%natural%Te   .or. bcs(bnd_type2)%natural%Te  ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_rhon) .and. (bcs(bnd_type1)%natural%rhon .or. bcs(bnd_type2)%natural%rhon))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_vpar) .and. (bcs(bnd_type1)%natural%vpar .or. bcs(bnd_type2)%natural%vpar))  apply_natural_bc(i_var)=.true.
enddo

! --- Sheath current BC (bcs%sheath_j): on edges whose both endpoints are sheath types the zj rows carry the
! --- characteristic zj = j_sat*(1 - exp(x)), x = Lambda - a_n*(u - C_V*V_wall)/(2Te), at the Gauss points where
! --- |b.n| >= sin(min_sheath_angle) (the nodes there have their Dirichlet u and zj rows released,
! --- mod_boundary_conditions; u follows from the vorticity equation). j_sat = c_sat*rho*(+-v_fl/|b.n|)/|B| with
! --- v_fl = factor*cs*|b.n| the parallel Bohm flow the nodal Mach-1 row imposes: the sheath-accelerated flux; the
! --- ExB part of the wall flux carries ions and electrons together and adds no net current.
sj_on = sheath_j_current_row .and. bcs(bnd_type1)%sheath_j .and. bcs(bnd_type2)%sheath_j
sj_an = 0.d0 ; sj_csat = 0.d0 ; sj_CT = 0.d0 ; sj_CV = 0.d0
if ( sj_on ) then
  apply_natural_bc(var_zj) = .true.
  call sheath_j_norm(sj_an, sj_csat)
  call floating_u_norm(sj_an, sj_CT, sj_CV)
endif

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

    ! --- Sheath-set wall flux (sf_on). The volume advection of rho and p is not integrated by parts, so with a wall
    ! --- potential that varies along the wall it carries an implicit ExB wall flux q*vE.n*R*dl in both directions,
    ! --- with no inflow datum where vE.n < 0. The sheath rows then impose the flux on the TOTAL normal flow,
    ! ---   Gamma = n*max(vn, v_fl),   vn = Vpar*(B_pol.n) + vE.n,   v_fl = factor*cs*|b.n|,
    ! --- the wall absorbs and never emits (SOLEDGE's Bohm-Chodura inequality on the wall flux; SOLPS's U_out floor).
    ! --- v_fl is exactly the parallel normal flow the nodal Mach-1 row imposes (Vpar = +-factor*cs/|B|, factor the
    ! --- vpar_smoothing weight, 1 without it), so with the row holding Gamma = n*(v_fl + max(vE.n,0)): unchanged
    ! --- where the ExB points out, the implicit emission replaced by the Bohm flux where it points in. In the rows,
    ! --- fx_n = max(vn, v_fl)*R*dl replaces develop's parallel measure vpar0*ps0_s*normal_sign3 = Vpar*(B_pol.n)*R*dl
    ! --- and the excess sink -q*ex_n, ex_n = max(v_fl - vn, 0)*R*dl, is added (q = rho, rho*Ti, rho*Te); boundary
    ! --- energy term oint q^2 (vn/2 - max(vn, v_fl)) <= 0 at every incidence angle. Nothing divides by b.n.
    ! --- fx_v, fx_p, fx_u, fx_T: coefficients of the trial Vpar, psi_s, u_s and cs(T) columns; on the floor branch
    ! --- |b.n|*R*dl = |ps0_s|/Btot with the Btot dependence on the normal psi derivative lagged.
    ! --- vE.n = -orient*R*u_s/dl outward (v_E = (-R*u_Z, +R*u_R)). The kinetic recycling uses the same Gamma
    ! --- (mod_particle_wall_interaction). Off the route: develop's expressions (same sparsity, products re-associated:
    ! --- measured 3e-14 relative against develop's routine in the serial harness).
    fx_n = vpar0 * ps0_s * normal_sign3
    fx_v =         ps0_s * normal_sign3
    fx_p = vpar0         * normal_sign3
    fx_u = 0.d0
    fx_T = 0.d0
    ex_n = 0.d0; ex_v = 0.d0; ex_p = 0.d0; ex_u = 0.d0; ex_T = 0.d0
    sf_fl = max(factor, 0.d0) * cs0 * abs(bdotn)
    if ( sf_on ) then
      sf_orient = sign(1.d0, y_s(ms)*normal(1) - x_s(ms)*normal(2))
      sf_vEn = - sf_orient * BigR * eq_s(mp,var_u,ms) / dl
      sf_Bn  = bdotn * Btot
      sf_vn  = sf_Bn * Vpar0 + sf_vEn
      if ( sf_vn .ge. sf_fl ) then          ! flow satisfies the condition: Gamma = n*vn
        fx_n = sf_vn * BigR * dl
        fx_v = sf_Bn * BigR * dl
        fx_p = vpar0 * normal_sign3
        fx_u = - sf_orient * BigR**2
      else                                  ! Gamma = n*v_fl, excess sink n*(v_fl - vn)
        fx_n = sf_fl * BigR * dl
        fx_v = 0.d0
        fx_p = max(factor, 0.d0) * cs0 * sign(1.d0, ps0_s) / Btot
        fx_u = 0.d0
        fx_T = max(factor, 0.d0) * abs(bdotn) * BigR * dl
        ex_n = ( sf_fl - sf_vn ) * BigR * dl
        ex_v = - sf_Bn * BigR * dl
        ex_p = fx_p - vpar0 * normal_sign3
        ex_u = + sf_orient * BigR**2
        ex_T = fx_T
      endif
    endif

    ! --- Sheath current row (sj_on) at this Gauss point: residual zj - j_sat*f(x), f = 1 - exp(min(x, Lambda))
    ! --- (electron current saturated at x >= Lambda, potential not below the wall), Zbig*dl weight, exact columns
    ! --- on zj, u, rho and T (through cs and through x). corr_neg-corrected Te and rho as in every natural row: a
    ! --- raw Te <= 0 would flip the sign of x and a raw rho <= 0 the sign of j_sat. The Btot and sign(B.n)
    ! --- dependence on psi is lagged, as everywhere in this routine.
    sj_here = sj_on .and. ( abs(bdotn) .ge. sin(c_angle) )
    sj_vfl  = max(factor, 0.d0) * cs0
    sj_jsat = 0.d0
    sc_w = 0.d0 ; sc_x = 0.d0 ; sc_ex = 1.d0 ; sc_f = 0.d0 ; sc_dfdu = 0.d0 ; sc_dfdTe = 0.d0 ; sc_res = 0.d0
    sc_Tc = 1.d0 ; sc_dTc = 1.d0 ; sc_drc = 1.d0
    if ( sj_here ) then
      sj_jsat = sj_csat * r0_corr * normal_sign * sj_vfl / Btot
      sc_Tc   = Te0_corr
      sc_dTc  = dcorr_neg_temp_dT(Te0)
      sc_drc  = dcorr_neg_dens_drho(r0)
      sc_x    = sheath_Lambda - sj_an * ( eq_g(mp,var_u,ms) - sj_CV*sheath_V_wall ) / (2.d0*sc_Tc)
      sc_ex   = exp( min(sc_x, sheath_Lambda) )
      sc_f    = 1.d0 - sc_ex
      if ( sc_x .lt. sheath_Lambda ) then
        sc_dfdu  =   sc_ex * sj_an / (2.d0*sc_Tc)
        sc_dfdTe = - sc_ex * sj_an * ( eq_g(mp,var_u,ms) - sj_CV*sheath_V_wall ) / (2.d0*sc_Tc**2) * sc_dTc
      endif
      sc_res  = eq_g(mp,var_zj,ms) - sj_jsat * sc_f
      sc_w    = Zbig * dl
      if ( mp .eq. 1 ) call wall_diag_sheath_add(bnd_type1, ws*dl, BigR, y_g(ms), eq_g(mp,var_zj,ms), sj_jsat, &
                                                 bdotn*Btot, eq_g(mp,var_u,ms), sc_x .ge. sheath_Lambda)
      if ( mp .eq. 1 .and. bnd_type2 .ne. bnd_type1 ) &
        call wall_diag_sheath_add(bnd_type2, ws*dl, BigR, y_g(ms), eq_g(mp,var_zj,ms), sj_jsat, &
                                  bdotn*Btot, eq_g(mp,var_u,ms), sc_x .ge. sheath_Lambda)
    endif

    ! --- Wall diagnostics (wall_diag): first toroidal plane, attributed to the types of both edge endpoints.
    ! --- vE.n = -orient*R*u_s/dl is the outward ExB normal speed (v_E = (-R*u_Z, +R*u_R)); B_pol.n = bdotn*Btot.
    if ( with_vpar .and. mp .eq. 1 ) then
      wd_orient = sign(1.d0, y_s(ms)*normal(1) - x_s(ms)*normal(2))
      wd_jr = 0.d0
      if ( sj_jsat .ne. 0.d0 ) wd_jr = eq_g(mp,var_zj,ms) / sj_jsat
      wd_vEn    = - wd_orient * BigR * eq_s(mp,var_u,ms) / dl
      call wall_diag_add(bnd_type1, ws*dl, BigR, y_g(ms), r0, Ti0, Te0, Vpar0, Btot, bdotn*Btot, bdotn, wd_vEn, cs0, &
                         eq_g(mp,var_u,ms), Te0_s/dl, eq_s(mp,var_u,ms)/dl, r0*ex_n, r0*sf_fl*BigR*dl, wd_jr, sc_x, eq_g(mp,var_zj,ms))
      if ( bnd_type2 .ne. bnd_type1 ) &
        call wall_diag_add(bnd_type2, ws*dl, BigR, y_g(ms), r0, Ti0, Te0, Vpar0, Btot, bdotn*Btot, bdotn, wd_vEn, cs0, &
                           eq_g(mp,var_u,ms), Te0_s/dl, eq_s(mp,var_u,ms)/dl, r0*ex_n, r0*sf_fl*BigR*dl, wd_jr, sc_x, eq_g(mp,var_zj,ms))
    endif


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

          ! --- Sheath current row (bcs%sheath_j)
          rhs_ij(var_zj) = - v * sc_w * sc_res

          ! --- Most B.C.s need vpar
          if (with_vpar) then

            ! --- Density reflection and minimum particle flux
            rhs_ij(var_rho)   = + v * density_reflection * r0      * fx_n * tstep     &
                                - v * r0      * ex_n * tstep                          & ! sheath-set wall flux: excess sink
                                - v * r0      * cs0 * BigR * dl * c_angle * tstep     ! particle flux at 1 degree angle  

            ! --- Sheath heat flux (c_angle for mininum heat fluxes at grazing angles)
            if (with_TiTe) then
              rhs_ij(var_Ti)  = - v * (gamma_sheath_i-1.d0) * r0 * Ti0 * fx_n * tstep &
                                - v *                         r0 * Ti0 * ex_n * tstep &
                                - v * (gamma_sheath_i-1.d0) * r0 * Ti0 * cs0    * BigR * dl * c_angle * tstep & 
                                - v * (GAMMA - 1.d0) * vpar0 * visco_par_heating * gradvpar0dotn * BigR * dl  * tstep  

              rhs_ij(var_Te)  = - v * (gamma_sheath_e-1.d0) * r0 * Te0 * fx_n * tstep &
                                - v *                         r0 * Te0 * ex_n * tstep &
                                - v * (gamma_sheath_e-1.d0) * r0 * Te0 * cs0  * BigR * dl * c_angle   * tstep  
            else
              rhs_ij(var_T)   = - v * (gamma_sheath  -1.d0) * r0 * T0  * fx_n * tstep &
                                - v *                         r0 * T0  * ex_n * tstep &
                                - v * (gamma_sheath  -1.d0) * r0 * T0  * cs0    * BigR * dl * c_angle * tstep & 
                                - v * (GAMMA - 1.d0) * vpar0 * visco_par_heating * gradvpar0dotn * BigR * dl  * tstep  
            endif

            ! --- Mach=1 through boundary integral penalization method
            rhs_ij(var_vpar) = - v * (vpar0 * Btot * normal_sign - cs0 * factor) * dl * Zbig  * factor_cs_bnd_integral 

            ! --- Fluid neutral reflection
            if (with_neutrals) then 
              rhs_ij(var_rhon) = rhs_ij(var_rhon)                                                   &
                               + v * neutral_reflection * r0 * fx_n * tstep &
                               + v * neutral_reflection * r0 * cs0 * BigR * dl * c_angle    * tstep ! particle flux at 1 degree angle  
            endif ! with_neutrals

          endif ! with_vpar
          index_ij = n_tor_local*n_var*n_degrees*(vertex(i)-1) + n_tor_local * n_var * (j2-1) + im - i_tor_min +1  ! index in the ELM matrix

          do i_var = 1, n_var
            if ( .not. apply_natural_bc(i_var) ) cycle
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

                cs_T   = gamma * T  / (2.d0 * cs0)
                cs_Ti  = gamma * Ti / (2.d0 * cs0)
                cs_Te  = gamma * Te / (2.d0 * cs0)

                ! --- Sheath current row: exact columns of zj - j_sat*f (amat is assigned per (k,l,in), never accumulated)
                amat(var_zj,var_zj)  =   v * sc_w * psi
                amat(var_zj,var_rho) = - v * sc_w * sj_csat * normal_sign * sj_vfl / Btot * sc_f * sc_drc * rho
                amat(var_zj,var_u)   = - v * sc_w * sj_jsat * sc_dfdu * psi
                if (with_TiTe) then
                  amat(var_zj,var_Ti)  = - v * sc_w * sj_csat * r0_corr * normal_sign * max(factor, 0.d0) / Btot * sc_f * cs_Ti
                  amat(var_zj,var_Te)  = - v * sc_w * ( sj_csat * r0_corr * normal_sign * max(factor, 0.d0) / Btot * sc_f * cs_Te &
                                                        + sj_jsat * sc_dfdTe * Te )
                else
                  amat(var_zj,var_T)   = - v * sc_w * ( sj_csat * r0_corr * normal_sign * max(factor, 0.d0) / Btot * sc_f * cs_T &
                                                        + sj_jsat * sc_dfdTe * 0.5d0 * T )
                endif

                ! --- Most of natural BCs need vpar
                if (with_vpar) then

                  ! --- Density reflection and minimum particle flux (c_angle)
                  amat(var_rho,var_psi)   = - v * density_reflection * r0  * psi_s * fx_p * theta * tstep &
                                            + v                      * r0  * psi_s * ex_p * theta * tstep
                  amat(var_rho,var_rho)   = - v * density_reflection * rho * fx_n * theta * tstep &
                                            + v                      * rho * ex_n * theta * tstep &
                                            + v                      * rho * cs0   * BigR * dl * c_angle  * theta * tstep 
                  amat(var_rho,var_vpar)  = - v * density_reflection * r0  * vpar * fx_v * theta * tstep &
                                            + v                      * r0  * vpar * ex_v * theta * tstep
                  amat(var_rho,var_u)     = - v * density_reflection * r0  * fx_u * psi_s * theta * tstep &
                                            + v                      * r0  * ex_u * psi_s * theta * tstep

                  ! --- Sheath heat flux
                  if (with_TiTe) then                
                    amat(var_rho,var_Ti)  = + v * r0 * cs_Ti * BigR * dl * c_angle  * theta * tstep &
                                            + v * r0 * cs_Ti * ( ex_T - density_reflection * fx_T ) * theta * tstep
                    amat(var_rho,var_Te)  = + v * r0 * cs_Te * BigR * dl * c_angle  * theta * tstep &
                                            + v * r0 * cs_Te * ( ex_T - density_reflection * fx_T ) * theta * tstep
                  else
                    amat(var_rho,var_T)   = + v * r0 * cs_T  * BigR * dl * c_angle  * theta * tstep &
                                            + v * r0 * cs_T  * ( ex_T - density_reflection * fx_T ) * theta * tstep
                  endif

                  ! --- Sheath heat flux
                  if (with_TiTe) then                
                    amat(var_Ti,var_psi)  = + v * (gamma_sheath_i-1.d0) * r0  * Ti0 * psi_s * fx_p * theta * tstep &
                                            + v *                         r0  * Ti0 * psi_s * ex_p * theta * tstep
                    amat(var_Ti,var_rho)  = + v * (gamma_sheath_i-1.d0) * rho * Ti0 * fx_n * theta * tstep & 
                                            + v *                         rho * Ti0 * ex_n * theta * tstep &
                                            + v * (gamma_sheath_i-1.d0) * rho * Ti0 * cs0   * BigR  * dl * c_angle * theta * tstep 
                    amat(var_Ti,var_Ti)   = + v * (gamma_sheath_i-1.d0) * r0  * Ti  * fx_n * theta * tstep & 
                                            + v *                         r0  * Ti  * ex_n * theta * tstep &
                                            + v * (gamma_sheath_i-1.d0) * r0  * Ti  * cs0   * BigR  * dl * c_angle * theta * tstep &
                                            + v * (gamma_sheath_i-1.d0) * r0  * Ti0 * cs_Ti * BigR  * dl * c_angle * theta * tstep

                    amat(var_Te,var_psi)  = + v * (gamma_sheath_e-1.d0) * r0  * Te0 * psi_s * fx_p * theta * tstep &
                                            + v *                         r0  * Te0 * psi_s * ex_p * theta * tstep
                    amat(var_Te,var_rho)  = + v * (gamma_sheath_e-1.d0) * rho * Te0 * fx_n * theta * tstep & 
                                            + v *                         rho * Te0 * ex_n * theta * tstep &
                                            + v * (gamma_sheath_e-1.d0) * rho * Te0 * cs0   * BigR  * dl * c_angle * theta * tstep 
                    amat(var_Te,var_Te)   = + v * (gamma_sheath_e-1.d0) * r0  * Te  * fx_n * theta * tstep & 
                                            + v *                         r0  * Te  * ex_n * theta * tstep &
                                            + v * (gamma_sheath_e-1.d0) * r0  * Te  * cs0   * BigR  * dl * c_angle * theta * tstep &
                                            + v * (gamma_sheath_e-1.d0) * r0  * Te0 * cs_Te * BigR  * dl * c_angle * theta * tstep

                    amat(var_Ti,var_vpar) = + v * (gamma_sheath_i-1.d0) * r0  * Ti0 * vpar * fx_v * theta * tstep &
                                            + v *                         r0  * Ti0 * vpar * ex_v * theta * tstep &
                                            + v * (GAMMA - 1.d0) * vpar * visco_par_heating * gradvpar0dotn * BigR * dl    * theta * tstep &
                                            + v * (GAMMA - 1.d0) * vpar0 * visco_par_heating * gradvpardotn * BigR * dl    * theta * tstep
                    amat(var_Te,var_vpar) = + v * (gamma_sheath_e-1.d0) * r0  * Te0 * vpar * fx_v * theta * tstep &
                                            + v *                         r0  * Te0 * vpar * ex_v * theta * tstep
                    amat(var_Ti,var_u)    = + v * (gamma_sheath_i-1.d0) * r0  * Ti0 * fx_u * psi_s * theta * tstep &
                                            + v *                         r0  * Ti0 * ex_u * psi_s * theta * tstep
                    amat(var_Te,var_u)    = + v * (gamma_sheath_e-1.d0) * r0  * Te0 * fx_u * psi_s * theta * tstep &
                                            + v *                         r0  * Te0 * ex_u * psi_s * theta * tstep
                    ! --- Sheath-set wall flux, floor branch: fx_n and ex_n depend on T through cs (both species). The
                    ! --- cross-species columns also carry the c_angle floor's cs dependence (sf_on only, so that runs
                    ! --- without the route keep develop's Jacobian pattern).
                    amat(var_Ti,var_Ti)   = amat(var_Ti,var_Ti) &
                                            + v * r0 * Ti0 * cs_Ti * ( (gamma_sheath_i-1.d0) * fx_T + ex_T ) * theta * tstep
                    amat(var_Ti,var_Te)   = + v * r0 * Ti0 * cs_Te * ( (gamma_sheath_i-1.d0) * fx_T + ex_T ) * theta * tstep &
                                            + v * (gamma_sheath_i-1.d0) * r0 * Ti0 * cs_Te * BigR * dl * c_angle * theta * tstep &
                                              * merge(1.d0, 0.d0, sf_on)
                    amat(var_Te,var_Te)   = amat(var_Te,var_Te) &
                                            + v * r0 * Te0 * cs_Te * ( (gamma_sheath_e-1.d0) * fx_T + ex_T ) * theta * tstep
                    amat(var_Te,var_Ti)   = + v * r0 * Te0 * cs_Ti * ( (gamma_sheath_e-1.d0) * fx_T + ex_T ) * theta * tstep &
                                            + v * (gamma_sheath_e-1.d0) * r0 * Te0 * cs_Ti * BigR * dl * c_angle * theta * tstep &
                                              * merge(1.d0, 0.d0, sf_on)
                  else
                    amat(var_T,var_psi)   = + v * (gamma_sheath  -1.d0) * r0  *  T0 * psi_s * fx_p * theta * tstep &
                                            + v *                         r0  *  T0 * psi_s * ex_p * theta * tstep
                    amat(var_T,var_rho)   = + v * (gamma_sheath  -1.d0) * rho *  T0 * fx_n * theta * tstep &
                                            + v *                         rho *  T0 * ex_n * theta * tstep &
                                            + v * (gamma_sheath  -1.d0) * rho *  T0 * cs0   * BigR  * dl * c_angle * theta * tstep 
                    amat(var_T,var_T)     = + v * (gamma_sheath  -1.d0) * r0  *  T  * fx_n * theta * tstep &
                                            + v *                         r0  *  T  * ex_n * theta * tstep &
                                            + v * (gamma_sheath  -1.d0) * r0  *  T  * cs0   * BigR  * dl * c_angle * theta * tstep &
                                            + v * (gamma_sheath  -1.d0) * r0  *  T0 * cs_T  * BigR  * dl * c_angle * theta * tstep &
                                            + v * r0 * T0 * cs_T * ( (gamma_sheath-1.d0) * fx_T + ex_T ) * theta * tstep

                    amat(var_T,var_u)     = + v * (gamma_sheath  -1.d0) * r0  * T0  * fx_u * psi_s * theta * tstep &
                                            + v *                         r0  * T0  * ex_u * psi_s * theta * tstep
                    amat(var_T,var_vpar)  = + v * (gamma_sheath  -1.d0) * r0  * T0  * vpar * fx_v * theta * tstep & 
                                            + v *                         r0  * T0  * vpar * ex_v * theta * tstep &
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
                    amat(var_rhon,var_psi) = - v * neutral_reflection * r0  * psi_s * fx_p      * theta * tstep 
  
                    amat(var_rhon,var_rho) = - v * neutral_reflection * rho     * fx_n      * theta * tstep &
                                             - v * neutral_reflection * rho     * cs0 * BigR * dl * c_angle * theta * tstep 
 
                    if (with_TiTe) then 
                      amat(var_rhon,var_Ti) = - v * neutral_reflection * r0 * cs_Ti * ( BigR * dl * c_angle + fx_T ) * theta * tstep 
                      amat(var_rhon,var_Te) = - v * neutral_reflection * r0 * cs_Te * ( BigR * dl * c_angle + fx_T ) * theta * tstep 
                    else
                      amat(var_rhon,var_T)  = - v * neutral_reflection * r0 * cs_T  * ( BigR * dl * c_angle + fx_T ) * theta * tstep 
                    endif
  
                    amat(var_rhon,var_vpar) = - v * neutral_reflection * r0 * vpar * fx_v     * theta * tstep 
                    amat(var_rhon,var_u)    = - v * neutral_reflection * r0 * fx_u * psi_s    * theta * tstep 
                  endif ! with neutrals

                endif   ! with_vpar
                index_kl = n_tor_local*n_var*n_degrees*(vertex(k)-1) + n_tor_local * n_var * (l2-1) + in - i_tor_min +1  ! index in the ELM matrix

                ! --- Add contributions to ELM matrix                 
                do k_var = 1, n_var
                  do i_var = 1, n_var

                    if ( .not. apply_natural_bc(i_var) ) cycle

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
