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
use mod_floating_diag, only: floating_diag_add, sheath_diag_add, wallj_diag_add
use mod_floating_u,    only: sheath_j_norm, floating_u_norm, sheath_j_ramp
use mod_plasma_functions, only: viscosity

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
logical    :: mw_on                                              ! weak Bohm condition (mach1_weak) on this edge
real*8     :: mw_orient, mw_vEn, mw_Bn, mw_vn, mw_tgt, mw_act, mw_cs, mw_res, mw_w, mw_in, mw_s, mw_x
real*8     :: fx_n, fx_v, fx_p, fx_u, fx_out                    ! normal-flow measure of the sheath fluxes and its columns
logical    :: sj_on, sj_here, sj_surf                            ! sheath rows on this edge / at this Gauss point / surface term on this edge
real*8     :: sj_an, sj_csat, sj_CT, sj_CV, sj_jsat, sj_psin, sj_w, sj_esp, sj_Tt
real*8     :: sc_x, sc_ex, sc_f, sc_dfdu, sc_dfdTe, sc_res, sc_w, sc_Tc, sc_dTc, sc_drc   ! current-slot form (sheath_j_current_row)
real*8     :: so_X, so_Xc, so_g, so_res, so_cu, so_czj, so_crho, so_cTe, so_ccs, so_w, so_al   ! sheath potential row
logical    :: so_capped
logical    :: sj_cf                                              ! cancel the implicit viscous and magnetisation wall currents of the u row on this edge (sheath_j_cancel_flux)
real*8     :: cf_fo, cf_fn, cf_nu, cf_dnu, cf_w0, cf_wx, cf_wy, cf_dwdn, cf_G, cf_vis, cf_mag, cf_p0, cf_p0x, cf_p0y
real*8     :: cf_dsx, cf_dsy, cf_dtx, cf_dty, cf_dp, cf_dpx, cf_dpy, cf_tx, cf_ty, cf_gTx, cf_gTy
real*8     :: wj_an, wj_csat, wj_cap, wj_mag, wj_exb, wj_vis, wj_p0s, wj_w0, wj_w0s, wj_w0t, wj_w0x, wj_w0y   ! [wall J] diagnostic
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

! --- If one of the nodes has a boundary type where natural BCs are applied, apply boundary integral for the full bnd element
do i_var=1, n_var
  if ( (i_var==var_rho ) .and. (bcs(bnd_type1)%natural%rho  .or. bcs(bnd_type2)%natural%rho ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_T   ) .and. (bcs(bnd_type1)%natural%T    .or. bcs(bnd_type2)%natural%T   ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_Ti  ) .and. (bcs(bnd_type1)%natural%Ti   .or. bcs(bnd_type2)%natural%Ti  ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_Te  ) .and. (bcs(bnd_type1)%natural%Te   .or. bcs(bnd_type2)%natural%Te  ))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_rhon) .and. (bcs(bnd_type1)%natural%rhon .or. bcs(bnd_type2)%natural%rhon))  apply_natural_bc(i_var)=.true.
  if ( (i_var==var_vpar) .and. (bcs(bnd_type1)%natural%vpar .or. bcs(bnd_type2)%natural%vpar))  apply_natural_bc(i_var)=.true.
enddo

! --- Weak Bohm condition on the total normal flow (mach1_weak). Carried by edges whose BOTH endpoints
! --- are mach1 types, i.e. a target edge and never a flux-surface edge; scattered through the Vpar rows.
mw_on = mach1_weak .and. with_vpar .and. bcs(bnd_type1)%mach1 .and. bcs(bnd_type2)%mach1
if ( mw_on ) apply_natural_bc(var_vpar) = .true.

! --- Sheath current BC (bcs%sheath_j) on edges whose both endpoints are sheath types, where
! --- |b.n| >= sin(min_sheath_angle): the u rows carry the weak sheath row (potential form) and the zj rows
! --- receive the surface term of the current definition, oint v*(dpsi/dn)/R dl, that the volume form drops
! --- (which is why JOREK pins zj); with it zj = Delta*psi holds at the wall and the wall current is the
! --- current the interior induction equation drives into the last element.
! --- The surface term follows the NODE: a released node's definition row needs it over the node's whole
! --- support, so it is assembled at every Gauss point of every edge with a sheath node at either end, with
! --- no angle gate (pinned nodes have their rows overwritten by Dirichlet anyway). Without this the row is
! --- the volume form alone over part of the support, i.e. dpsi/dn = 0 weakly, and the wall current absorbs
! --- a missing flux of order R*B_t/h: thousands of j_sat at the edge of the sheath region (measured).
sj_surf = ( bcs(bnd_type1)%sheath_j .or. bcs(bnd_type2)%sheath_j ) .and. (.not. sheath_j_current_row)
sj_on   = bcs(bnd_type1)%sheath_j .and. bcs(bnd_type2)%sheath_j
sj_an = 0.d0 ; sj_csat = 0.d0 ; sj_CT = 0.d0 ; sj_CV = 0.d0
if ( sj_surf ) apply_natural_bc(var_zj) = .true.
if ( sj_on ) then
  if ( sheath_j_current_row ) then
    apply_natural_bc(var_zj) = .true.        ! the characteristic as the zj row; u keeps the vorticity equation
  else if ( .not. sheath_j_float_u ) then
    apply_natural_bc(var_u)  = .true.
  endif
  call sheath_j_norm(sj_an, sj_csat)
  call floating_u_norm(sj_an, sj_CT, sj_CV)
endif
so_al = sheath_j_ramp(t_now)     ! ramp factor of the current dependence, 0 (floating) -> 1 (full characteristic)

! --- Cancelled wall fluxes (sheath_j_cancel_flux). The vorticity equation is assembled integrated by parts, so at
! --- a released wall node its row contains the boundary fluxes of those terms, currents no sheath sets: the
! --- viscous vorticity flux from the Dirichlet w (frozen at its t=0 wall value while the interior w follows the
! --- potential: measured 4-20 j_sat at the strike points, 2026-09-25) and the magnetisation current, the
! --- total-derivative part of the pressure-bracket flux R^2*dp/ds, which a target does not collect (Rozhansky
! --- 2001; SOLPS closes the target balance without it). Both are added back with the opposite sign on every
! --- edge with a sheath node at either end (a released value DOF's support spans both incident edges; on pinned
! --- DOFs the Dirichlet weight makes the term invisible). What remains in the balance is the polarisation flux,
! --- the grad-B part of the pressure flux, 2*R*p*n_Z, and the sheath current. Weakly this is the zero viscous
! --- current condition dw/dn = 0 of the drift-fluid codes, with w itself still Dirichlet.
sj_cf = sheath_j_cancel_flux .and. sheath_j_current_row .and. ( bcs(bnd_type1)%sheath_j .or. bcs(bnd_type2)%sheath_j )
if ( sj_cf ) apply_natural_bc(var_u) = .true.

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

    ! --- Weak Bohm condition (marginal form):  res = (B_pol.n)*Vpar - cs*|b.n|,  i.e. Vpar = +-cs/|B|,
    ! --- weighted by d(res)/d(Vpar) = B_pol.n so that a grazing point loses authority as (B.n)^2 and the
    ! --- natural Vpar condition takes over continuously. The ExB drift is NOT compensated by the parallel
    ! --- flow: compensating an inward drift of 10 km/s at a few percent incidence would demand a
    ! --- several-times-sonic parallel flow in the last element, whose divergence must then cancel the
    ! --- ExB inflow term to O(1) within that element - the dispersive density dipole at the strike
    ! --- point. Where the total flow is inward the inflow closure below supplies the density datum.
    ! --- mach1_weak_drift restores the SOLPS non-marginal form, target max(cs*|b.n| - vE.n, 0). With
    ! --- mach1_weak_drift_bound = c > 0 the compensation of an INWARD drift d = -vE.n > 0 is bounded
    ! --- smoothly at s = c*cs*|b.n| (SOLPS b2stbc_cbc, c = 2): target cs*|b.n| + s*tanh(d/s), whose u
    ! --- column carries sech^2(d/s) and so fades out where the bound takes over, with no branch. The
    ! --- flow beyond the bound is inward and is handled by the inflow closure below. With
    ! --- mach1_weak_drift_cut the compensation is dropped where |b.n| < sin(min_sheath_angle), the angle
    ! --- below which the sheath particle and heat fluxes are already taken from the c_angle floor model;
    ! --- the row there is the marginal one. Static map: psi is Dirichlet on the wall.
    ! --- vE.n = -orient*R*u_s/dl is the outward ExB normal flow for v_E = (-R*u_Z, +R*u_R) and the edge
    ! --- tangent (x_s, y_s)/dl. Vpar*(B_pol.n) is the parallel normal flow, a velocity since v = Vpar*B.
    mw_orient = sign(1.d0, y_s(ms)*normal(1) - x_s(ms)*normal(2))
    mw_vEn    = - mw_orient * BigR * eq_s(mp,var_u,ms) / dl
    mw_Bn     = bdotn * Btot
    mw_tgt    = cs0 * abs(bdotn)
    mw_act    = 0.d0                      ! coefficient of the u column: drift not in the row
    mw_cs     = 1.d0                      ! coefficient of the temperature columns (cs in the target)
    if ( mach1_weak_drift .and. .not. ( mach1_weak_drift_cut .and. abs(bdotn) .lt. sin(c_angle) ) ) then
      mw_tgt = cs0 * abs(bdotn) - mw_vEn
      mw_act = 1.d0
      if ( mw_tgt .le. 0.d0 ) then        ! the drift alone gives sonic outflow or more: nothing to impose
        mw_tgt = 0.d0
        mw_act = 0.d0
        mw_cs  = 0.d0
      elseif ( mach1_weak_drift_bound .gt. 0.d0 .and. mw_vEn .lt. 0.d0 ) then
        mw_s   = mach1_weak_drift_bound * cs0 * abs(bdotn)
        mw_x   = - mw_vEn / mw_s
        mw_tgt = cs0 * abs(bdotn) + mw_s * tanh(mw_x)
        mw_act = 1.d0 / cosh(mw_x)**2
        mw_cs  = 1.d0 + mach1_weak_drift_bound * ( tanh(mw_x) - mw_x / cosh(mw_x)**2 )   ! d(target)/d(cs) / |b.n|
      endif
    endif
    mw_res    = mw_Bn * Vpar0 - mw_tgt
    ! --- mach1_weak_cut: below sin(min_sheath_angle) no Vpar row at all - the parallel viscosity's natural
    ! --- condition grad(Vpar).n = 0 holds there - while the inflow closure and the total-flow fluxes stay.
    mw_w      = 0.d0
    if ( mw_on ) mw_w = Zbig * mw_Bn * dl
    if ( mach1_weak_cut .and. abs(bdotn) .lt. sin(c_angle) ) mw_w = 0.d0

    ! --- Inflow closure. Where the TOTAL normal flow vn = Vpar*(B_pol.n) + vE.n is inward the density
    ! --- equation, advected with the undifferentiated test function, has no boundary datum. The weak
    ! --- inflow term  -oint v*min(vn,0)*(rho - rho_in) dl  with rho_in = 0 (no plasma enters from a wall)
    ! --- supplies it: its coefficient is the physical inflow rate, it vanishes identically where the flow
    ! --- is outward, and it is what takes over where the weak Bohm row above loses authority at grazing
    ! --- incidence. Same edges as that row. The temperatures get no term: with rho_in = 0 nothing enters.
    mw_vn = mw_Bn * Vpar0 + mw_vEn
    mw_in = 0.d0
    if ( mw_on .and. mach1_weak_inflow .and. (mw_vn .lt. 0.d0) ) mw_in = 1.d0

    ! --- Normal flow carried by the sheath particle and energy fluxes below: Vpar*(B_pol.n)*R*dl, which is
    ! --- what vpar0*ps0_s*normal_sign3 is. Under mach1_weak it is the TOTAL outgoing flow max(vn,0)*R*dl,
    ! --- so the momentum, particle, energy and kinetic recycling channels all see one and the same flow.
    ! --- fx_v, fx_p, fx_u are the coefficients of the trial Vpar, psi_s and psi_s(u) in its columns.
    fx_n = vpar0 * ps0_s * normal_sign3
    fx_v =         ps0_s * normal_sign3
    fx_p = vpar0         * normal_sign3
    fx_u = 0.d0
    if ( mw_on ) then
      fx_out = 0.d0
      if ( mw_vn .gt. 0.d0 ) fx_out = 1.d0
      fx_n = max(mw_vn, 0.d0) * BigR * dl
      fx_v = fx_out * mw_Bn * BigR * dl
      fx_p = fx_out * vpar0 * normal_sign3
      fx_u = - fx_out * mw_orient * BigR**2
    endif

    ! --- Sheath BC at this Gauss point. j_sat = c_sat*rho*Vpar_Bohm, Vpar_Bohm = sign(B.n)*cs/|B| (Artola eq. 5
    ! --- with the marginal Bohm row). Written for the current INTO the wall, -zj*(B_pol.n)/F0 = e*n*cs*|b.n|,
    ! --- independent of the sign of F0.
    sj_here = sj_on .and. ( abs(bdotn) .ge. sin(c_angle) )
    sj_jsat = 0.d0 ; sj_psin = 0.d0 ; sj_w = 0.d0
    if ( sj_here ) sj_jsat = sj_csat * r0_corr * normal_sign * cs0 / Btot     ! corr_neg rho: j_sat keeps its sign
    if ( sj_surf ) then
      sj_psin = ps0_x * normal(1) + ps0_y * normal(2)     ! dpsi/dn, outward
      sj_w    = dl / BigR                                 ! weight of the surface term
    endif

    ! --- The sheath row: the characteristic solved for the potential,
    ! ---   u = C_V*V_wall + (2Te/a_n)*(Lambda - ln X),   X = 1 - zj/j_sat,
    ! --- as a weak row on the u trace with weight one: d(res)/du = 1, so the potential is anchored
    ! --- everywhere, including where the sheath conductance vanishes. X is bounded to [e^-Lambda, e^Lambda]:
    ! --- the upper bound is electron saturation (potential not below the wall), the lower one says the
    ! --- characteristic is trusted up to twice the floating drop (a model statement) and keeps the row
    ! --- finite on the ion-saturated branch. Where a bound is active the zj/rho/cs columns vanish.
    so_w = 0.d0 ; so_res = 0.d0 ; so_cu = 0.d0 ; so_czj = 0.d0 ; so_crho = 0.d0 ; so_cTe = 0.d0 ; so_ccs = 0.d0
    so_capped = .false.
    ! --- Current-slot form (sheath_j_current_row): zj = j_sat*(1 - exp(x)), x = Lambda - a_n*(u - C_V*V_wall)/(2Te),
    ! --- electron current saturated at x >= Lambda, weight one, in the zj row; u is then set by the vorticity
    ! --- equation. The structure of the old weak-trace route that ran ~3900 steps on type 1 alone.
    sc_w = 0.d0 ; sc_x = 0.d0 ; sc_ex = 1.d0 ; sc_f = 0.d0 ; sc_dfdu = 0.d0 ; sc_dfdTe = 0.d0 ; sc_res = 0.d0 ; sc_Tc = 1.d0 ; sc_dTc = 1.d0 ; sc_drc = 1.d0
    if ( sj_here .and. sheath_j_current_row ) then
      ! --- corr_neg-corrected Te and rho, as in every other natural row: a raw Te <= 0 would flip the sign
      ! --- of x and a raw rho <= 0 the sign of j_sat. sc_dTc = d(Te_corr)/dTe carries into the Te column.
      sc_Tc  = Te0_corr
      sc_dTc = dcorr_neg_temp_dT(Te0)
      ! --- Ion branch (x < 0, Phi above floating): f = 1 - exp(x) - s*x with s = sheath_j_ion_slope, the
      ! --- finite slope of ion saturation (sheath expansion). With s = 0 the characteristic has no voltage
      ! --- root wherever the plasma delivers j >= j_sat, and the vorticity row then drives Phi to infinity
      ! --- there (measured: 234 -> 683 -> 2600 V at the outer target in 178 steps).
      sc_x     = sheath_Lambda - sj_an * ( eq_g(mp,var_u,ms) - sj_CV*sheath_V_wall ) / (2.d0*sc_Tc)
      sc_ex    = exp( min(sc_x, sheath_Lambda) )
      sc_f     = 1.d0 - sc_ex - sheath_j_ion_slope * min(sc_x, 0.d0)
      if ( sc_x .lt. sheath_Lambda ) then
        sc_dfdu  =   sc_ex * sj_an / (2.d0*sc_Tc)
        sc_dfdTe = - sc_ex * sj_an * ( eq_g(mp,var_u,ms) - sj_CV*sheath_V_wall ) / (2.d0*sc_Tc**2) * sc_dTc
      endif
      if ( sc_x .lt. 0.d0 ) then                                      ! d(-s*x)/du and /dTe
        sc_dfdu  = sc_dfdu  + sheath_j_ion_slope * sj_an / (2.d0*sc_Tc)
        sc_dfdTe = sc_dfdTe - sheath_j_ion_slope * sj_an * ( eq_g(mp,var_u,ms) - sj_CV*sheath_V_wall ) / (2.d0*sc_Tc**2) * sc_dTc
      endif
      ! --- Beyond electron saturation (x > Lambda, Phi below the wall): f = 1 - e^Lambda - s_e*(x - Lambda),
      ! --- so the row keeps a conductance where the plasma pushes more than the thermal electron flux
      ! --- (measured: j/j_sat = -26 against the cap's -19 at the outer target, Phi then unanchored below 0).
      if ( sc_x .gt. sheath_Lambda ) then
        sc_f     = sc_f - sheath_j_e_slope * ( sc_x - sheath_Lambda )
        sc_dfdu  = sc_dfdu  + sheath_j_e_slope * sj_an / (2.d0*sc_Tc)
        sc_dfdTe = sc_dfdTe - sheath_j_e_slope * sj_an * ( eq_g(mp,var_u,ms) - sj_CV*sheath_V_wall ) / (2.d0*sc_Tc**2) * sc_dTc
      endif
      sc_res   = eq_g(mp,var_zj,ms) - sj_jsat * sc_f
      sc_drc   = dcorr_neg_dens_drho(r0)                       ! d(rho_corr)/d(rho) for the rho column
      sc_w     = Zbig * dl
      so_capped = ( sc_x .ge. sheath_Lambda )
    endif

    if ( sj_here .and. .not. sheath_j_float_u .and. .not. sheath_j_current_row ) then
      so_X  = 1.d0
      if ( sj_jsat .ne. 0.d0 ) so_X = 1.d0 - eq_g(mp,var_zj,ms) / sj_jsat
      so_Xc = min( max( so_X, exp(-sheath_Lambda) ), exp(sheath_Lambda) )
      so_capped = ( so_Xc .ne. so_X ) .or. ( sj_jsat .eq. 0.d0 )
      so_g   = 2.d0*Te0/sj_an * ( sheath_Lambda - so_al*log(so_Xc) )
      so_res = eq_g(mp,var_u,ms) - sj_CV*sheath_V_wall - so_g
      so_cu  = 1.d0
      so_cTe = - 2.d0/sj_an * ( sheath_Lambda - so_al*log(so_Xc) )           ! d(res)/dTe through the prefactor
      if ( .not. so_capped ) then
        so_czj  = - so_al * 2.d0*Te0/sj_an / (so_Xc * sj_jsat)                        ! d(res)/dzj
        so_crho =   so_al * 2.d0*Te0/sj_an * eq_g(mp,var_zj,ms) / (so_Xc * sj_jsat * r0_corr) * dcorr_neg_dens_drho(r0)   ! through j_sat ~ rho_corr
        so_ccs  =   so_al * 2.d0*Te0/sj_an * eq_g(mp,var_zj,ms) / (so_Xc * sj_jsat * cs0)  ! through j_sat ~ cs, times cs_T
      endif
      so_w = Zbig * dl
    endif

    if ( floating_u_diag .and. sj_here ) &
      call sheath_diag_add(bnd_type1, ws*dl, eq_g(mp,var_zj,ms), sj_jsat, so_capped, &
                           eq_g(mp,var_u,ms), mw_Bn, BigR, y_g(ms), so_capped)

    ! --- [wall J] diagnostic (print only): the implicit boundary currents of the released vorticity row at this
    ! --- wall point against the sheath's capacity, all per unit edge parameter as they enter the u row with the
    ! --- same test function: capacity j_sat*|psi_s| (the wall flux of v*[psi,zj]); magnetisation current
    ! --- R^2*|p_s| (wall flux of the pressure bracket R^2*[v,p]); ExB advection of vorticity rho*R^2*|w|*|u_s|
    ! --- (wall flux of rho*R^2*w*[v,u]); viscous flux visco*R^3*|dw/dn|*dl (wall flux of -visco*R^3*grad v.grad w,
    ! --- with the constant visco as an estimate of visco_T). Every type, first plane.
    if ( floating_u_diag .and. mp .eq. 1 ) then
      call sheath_j_norm(wj_an, wj_csat)
      wj_cap = abs( wj_csat * r0_corr * cs0 / Btot ) * abs(ps0_s)
      if (with_TiTe) then
        wj_p0s = r0_s * (Ti0 + Te0) + r0 * (Ti0_s + Te0_s)
      else
        wj_p0s = r0_s * T0 + r0 * T0_s
      endif
      wj_mag = BigR**2 * abs(wj_p0s)
      wj_w0  = eq_g(mp,var_w,ms) ; wj_w0s = eq_s(mp,var_w,ms) ; wj_w0t = eq_t(mp,var_w,ms)
      wj_w0x = (   y_t(ms) * wj_w0s - y_s(ms) * wj_w0t ) / xjac
      wj_w0y = ( - x_t(ms) * wj_w0s + x_s(ms) * wj_w0t ) / xjac
      wj_exb = r0 * BigR**2 * abs(wj_w0) * abs(eq_s(mp,var_u,ms))
      wj_vis = visco * BigR**3 * abs( wj_w0x*normal(1) + wj_w0y*normal(2) ) * dl
      call wallj_diag_add(bnd_type1, ws, BigR, y_g(ms), wj_cap, wj_mag, wj_exb, wj_vis)
      if ( bnd_type2 .ne. bnd_type1 ) call wallj_diag_add(bnd_type2, ws, BigR, y_g(ms), wj_cap, wj_mag, wj_exb, wj_vis)
    endif

    ! --- Cancelled wall fluxes at this Gauss point (per unit edge parameter, added to the u row as +v*(...)*tstep):
    ! ---   viscous:       + visco_T*R*( fo*R^2*dw/dn + fn*2*R*w*n_R )*dl,  the boundary flux of the assembled
    ! ---                  -visco_T*R*grad v.grad(R^2 w) (fo = 1, fn = 1; visco_old_setup: fo = 1/R^2, fn = 0);
    ! ---   magnetisation: - d/ds(R^2 p)*dl = -( -n_Z*(2Rp + R^2 p_R) + n_R*R^2 p_Z )*dl,  the total-derivative
    ! ---                  part of the flux R^2*dp/ds of the assembled R^2*[v,p]; p = rho*(Ti+Te).
    ! --- n is the outward normal, s runs along the wall with n rotated by +90 degrees. cf_ds*/cf_dt*: (R,Z)
    ! --- derivatives per unit s- and t-derivative of a basis function, for the columns.
    cf_vis = 0.d0 ; cf_mag = 0.d0 ; cf_G = 0.d0 ; cf_nu = 0.d0 ; cf_dnu = 0.d0 ; cf_fo = 0.d0 ; cf_fn = 0.d0
    cf_dsx = 0.d0 ; cf_dsy = 0.d0 ; cf_dtx = 0.d0 ; cf_dty = 0.d0
    if ( sj_cf ) then
      cf_fo = 1.d0 ; cf_fn = 1.d0
      if ( visco_old_setup ) then
        cf_fo = 1.d0 / BigR**2 ; cf_fn = 0.d0
      endif
      if (with_TiTe) then
        call viscosity(visco, Te0, Te0_corr, Te_0, cf_nu, cf_dnu)
        cf_dnu = cf_dnu * dcorr_neg_temp_dT(Te0)
      else
        call viscosity(visco, T0, T0_corr, T_0, cf_nu, cf_dnu)
        cf_dnu = cf_dnu * dcorr_neg_temp_dT(T0)
      endif
      cf_dsx =   y_t(ms) / xjac ; cf_dsy = - x_t(ms) / xjac
      cf_dtx = - y_s(ms) / xjac ; cf_dty =   x_s(ms) / xjac
      cf_w0   = eq_g(mp,var_w,ms)
      cf_wx   = cf_dsx * eq_s(mp,var_w,ms) + cf_dtx * eq_t(mp,var_w,ms)
      cf_wy   = cf_dsy * eq_s(mp,var_w,ms) + cf_dty * eq_t(mp,var_w,ms)
      cf_dwdn = cf_wx * normal(1) + cf_wy * normal(2)
      cf_G    = BigR * ( cf_fo * BigR**2 * cf_dwdn + cf_fn * 2.d0 * BigR * cf_w0 * normal(1) ) * dl
      cf_vis  = cf_nu * cf_G
      cf_p0   = r0 * T0
      cf_p0x  = r0_x * T0 + r0 * ( Ti0_x + Te0_x )
      cf_p0y  = r0_y * T0 + r0 * ( Ti0_y + Te0_y )
      cf_mag  = - ( - normal(2) * ( 2.d0 * BigR * cf_p0 + BigR**2 * cf_p0x ) + normal(1) * BigR**2 * cf_p0y ) * dl
    endif

    if ( floating_u_diag .and. mw_on ) then
      call floating_diag_add(bnd_type1, ws*dl, mw_vn, mw_vEn, mw_Bn, mw_res, cs0*abs(bdotn), abs(Vpar0)*Btot/cs0, r0, Ti0, Te0, BigR, y_g(ms))
      if ( bnd_type2 .ne. bnd_type1 ) &
        call floating_diag_add(bnd_type2, ws*dl, mw_vn, mw_vEn, mw_Bn, mw_res, cs0*abs(bdotn), abs(Vpar0)*Btot/cs0, r0, Ti0, Te0, BigR, y_g(ms))
    endif

    c_1 = vpar_smoothing_coef(1); c_2 = vpar_smoothing_coef(2); c_3 = vpar_smoothing_coef(3)
    if (vpar_smoothing) then
      factor = 0.25d0 * ( 1.d0 + tanh( (abs(bdotn) - c_1) / c_2 ) )**2 - c_3
    else
      factor = 1.d0
    endif

    factor_cs_bnd_integral = 0.d0
    if (mach_one_bnd_integral) factor_cs_bnd_integral = 1.d0

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
            rhs_ij(var_rho)   = + v * density_reflection * r0      * fx_n * tstep     &
                                - v * r0      * cs0 * BigR * dl * c_angle * tstep     ! particle flux at 1 degree angle  

            ! --- Weak inflow term (mach1_weak): rho relaxes to rho_in = 0 at the inflow rate |vn|
            rhs_ij(var_rho)   = rhs_ij(var_rho) + v * mw_in * mw_vn * r0 * BigR * dl * tstep

            ! --- Sheath BC (bcs%sheath_j): surface term of the current definition on the zj row, and the
            ! --- weak potential row on the u row. The volume row is -(grad v . grad psi + v*zj)/R with the
            ! --- boundary integral + oint v*(dpsi/dn)/R dl dropped; here it is put back.
            rhs_ij(var_zj)    = + v * sj_w * sj_psin - v * sc_w * sc_res
            rhs_ij(var_u)     = - v * so_w * so_res

            ! --- Cancelled wall fluxes of the vorticity row (sheath_j_cancel_flux)
            rhs_ij(var_u)     = rhs_ij(var_u) + v * ( cf_vis + cf_mag ) * tstep

            ! --- Sheath heat flux (c_angle for mininum heat fluxes at grazing angles)
            if (with_TiTe) then
              rhs_ij(var_Ti)  = - v * (gamma_sheath_i-1.d0) * r0 * Ti0 * fx_n * tstep &
                                - v * (gamma_sheath_i-1.d0) * r0 * Ti0 * cs0    * BigR * dl * c_angle * tstep & 
                                - v * (GAMMA - 1.d0) * vpar0 * visco_par_heating * gradvpar0dotn * BigR * dl  * tstep  

              rhs_ij(var_Te)  = - v * (gamma_sheath_e-1.d0) * r0 * Te0 * fx_n * tstep &
                                - v * (gamma_sheath_e-1.d0) * r0 * Te0 * cs0  * BigR * dl * c_angle   * tstep  
            else
              rhs_ij(var_T)   = - v * (gamma_sheath  -1.d0) * r0 * T0  * fx_n * tstep &
                                - v * (gamma_sheath  -1.d0) * r0 * T0  * cs0    * BigR * dl * c_angle * tstep & 
                                - v * (GAMMA - 1.d0) * vpar0 * visco_par_heating * gradvpar0dotn * BigR * dl  * tstep  
            endif

            ! --- Mach=1 through boundary integral penalization method
            rhs_ij(var_vpar) = - v * (vpar0 * Btot * normal_sign - cs0 * factor) * dl * Zbig  * factor_cs_bnd_integral 

            ! --- Weak Bohm condition on the total normal flow (mach1_weak)
            rhs_ij(var_vpar) = rhs_ij(var_vpar) - v * mw_w * mw_res

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

                ! --- Most of natural BCs need vpar
                if (with_vpar) then

                  ! --- Density reflection and minimum particle flux (c_angle)
                  amat(var_rho,var_psi)   = - v * density_reflection * r0  * psi_s * fx_p * theta * tstep 
                  amat(var_rho,var_rho)   = - v * density_reflection * rho * fx_n * theta * tstep &
                                            + v                      * rho * cs0   * BigR * dl * c_angle  * theta * tstep 
                  amat(var_rho,var_vpar)  = - v * density_reflection * r0  * vpar * fx_v * theta * tstep 

                  ! --- Weak inflow term (mach1_weak): exact columns of vn*rho, vn = Vpar*(B_pol.n) - orient*R*u_s/dl
                  amat(var_rho,var_rho)   = amat(var_rho,var_rho)  - v * mw_in * mw_vn * rho  * BigR * dl * theta * tstep
                  amat(var_rho,var_vpar)  = amat(var_rho,var_vpar) - v * mw_in * mw_Bn * vpar * r0 * BigR * dl * theta * tstep
                  amat(var_rho,var_u)     = - v * density_reflection * r0 * fx_u * psi_s * theta * tstep &
                                            + v * mw_in * mw_orient * BigR**2 * psi_s * r0 * theta * tstep

                  ! --- Surface term of the current definition: dpsi/dn = [ps0_s*(y_t*n1 - x_t*n2) + ps0_t*(-y_s*n1 + x_s*n2)]/xjac.
                  ! --- Trace DOFs (psi_s) here; the normal-derivative DOFs (psi_t) in the extra loop below.
                  amat(var_zj,var_psi) = - v * sj_w * ( y_t(ms)*normal(1) - x_t(ms)*normal(2) ) / xjac * psi_s

                  ! --- Current-slot form: exact columns of zj - j_sat*f
                  amat(var_zj,var_zj)  =   v * sc_w * psi          ! amat is assigned per (k,l,in), never accumulated
                  amat(var_zj,var_rho) = - v * sc_w * sj_csat * normal_sign * cs0 / Btot * sc_f * sc_drc * rho
                  amat(var_zj,var_u)   = - v * sc_w * sj_jsat * sc_dfdu * psi
                  if (with_TiTe) then
                    amat(var_zj,var_Ti)  = - v * sc_w * sj_csat * r0 * normal_sign / Btot * sc_f * cs_Ti
                    amat(var_zj,var_Te)  = - v * sc_w * ( sj_csat * r0 * normal_sign / Btot * sc_f * cs_Te + sj_jsat * sc_dfdTe * Te )
                  else
                    amat(var_zj,var_T)   = - v * sc_w * ( sj_csat * r0 * normal_sign / Btot * sc_f * cs_T + sj_jsat * sc_dfdTe * 0.5d0 * T )
                  endif

                  ! --- Sheath potential row: exact columns of u - C_V*V_wall - (2Te/a_n)*(Lambda - ln X)
                  amat(var_u,var_u)    =   v * so_w * so_cu  * psi
                  amat(var_u,var_zj)   =   v * so_w * so_czj * psi
                  amat(var_u,var_rho)  =   v * so_w * so_crho * rho
                  if (with_TiTe) then
                    amat(var_u,var_Ti) =   v * so_w * so_ccs * cs_Ti
                    amat(var_u,var_Te) =   v * so_w * ( so_cTe * Te + so_ccs * cs_Te )
                  else
                    amat(var_u,var_T)  =   v * so_w * ( so_cTe * 0.5d0 * T + so_ccs * cs_T )
                  endif

                  ! --- Cancelled wall fluxes: columns of the trace DOFs (value psi, edge derivative psi_s); the
                  ! --- normal-derivative DOFs in the extra loop below. amat = -theta*tstep*d(rhs/tstep)/d(dof).
                  amat(var_u,var_w) = - v * cf_nu * BigR * ( cf_fo * BigR**2 * ( cf_dsx*normal(1) + cf_dsy*normal(2) ) * psi_s &
                                                           + cf_fn * 2.d0 * BigR * normal(1) * psi ) * dl * theta * tstep
                  if ( sj_cf ) then
                    cf_gTx = Ti0_x + Te0_x ; cf_gTy = Ti0_y + Te0_y
                    ! rho: dp = psi*T0, dp_R = psi_R*T0 + psi*T0_R, dp_Z likewise
                    cf_dp  = psi * T0 ; cf_dpx = cf_dsx*psi_s * T0 + psi * cf_gTx ; cf_dpy = cf_dsy*psi_s * T0 + psi * cf_gTy
                    amat(var_u,var_rho) = amat(var_u,var_rho) + v * ( - normal(2) * ( 2.d0*BigR*cf_dp + BigR**2*cf_dpx ) &
                                                                      + normal(1) * BigR**2 * cf_dpy ) * dl * theta * tstep
                    ! temperatures: dp = r0*psi, dp_R = r0_R*psi + r0*psi_R; the viscosity's own temperature dependence
                    cf_dp  = r0 * psi ; cf_dpx = r0_x * psi + r0 * cf_dsx*psi_s ; cf_dpy = r0_y * psi + r0 * cf_dsy*psi_s
                    if (with_TiTe) then
                      amat(var_u,var_Ti) = amat(var_u,var_Ti) + v * ( - normal(2) * ( 2.d0*BigR*cf_dp + BigR**2*cf_dpx ) &
                                                                      + normal(1) * BigR**2 * cf_dpy ) * dl * theta * tstep
                      amat(var_u,var_Te) = amat(var_u,var_Te) + v * ( - normal(2) * ( 2.d0*BigR*cf_dp + BigR**2*cf_dpx ) &
                                                                      + normal(1) * BigR**2 * cf_dpy ) * dl * theta * tstep &
                                                              - v * cf_dnu * cf_G * psi * theta * tstep
                    else
                      amat(var_u,var_T)  = amat(var_u,var_T)  + v * ( - normal(2) * ( 2.d0*BigR*cf_dp + BigR**2*cf_dpx ) &
                                                                      + normal(1) * BigR**2 * cf_dpy ) * dl * theta * tstep &
                                                              - v * cf_dnu * cf_G * psi * theta * tstep
                    endif
                  endif

                  ! --- Sheath heat flux
                  if (with_TiTe) then                
                    amat(var_rho,var_Ti)  = + v * r0 * cs_Ti * BigR * dl * c_angle  * theta * tstep
                    amat(var_rho,var_Te)  = + v * r0 * cs_Te * BigR * dl * c_angle  * theta * tstep
                  else
                    amat(var_rho,var_T)   = + v * r0 * cs_T  * BigR * dl * c_angle  * theta * tstep
                  endif

                  ! --- Sheath heat flux
                  if (with_TiTe) then                
                    amat(var_Ti,var_psi)  = + v * (gamma_sheath_i-1.d0) * r0  * Ti0 * psi_s * fx_p * theta * tstep 
                    amat(var_Ti,var_rho)  = + v * (gamma_sheath_i-1.d0) * rho * Ti0 * fx_n * theta * tstep & 
                                            + v * (gamma_sheath_i-1.d0) * rho * Ti0 * cs0   * BigR  * dl * c_angle * theta * tstep 
                    amat(var_Ti,var_Ti)   = + v * (gamma_sheath_i-1.d0) * r0  * Ti  * fx_n * theta * tstep & 
                                            + v * (gamma_sheath_i-1.d0) * r0  * Ti  * cs0   * BigR  * dl * c_angle * theta * tstep &
                                            + v * (gamma_sheath_i-1.d0) * r0  * Ti0 * cs_Ti * BigR  * dl * c_angle * theta * tstep

                    amat(var_Te,var_psi)  = + v * (gamma_sheath_e-1.d0) * r0  * Te0 * psi_s * fx_p * theta * tstep 
                    amat(var_Te,var_rho)  = + v * (gamma_sheath_e-1.d0) * rho * Te0 * fx_n * theta * tstep & 
                                            + v * (gamma_sheath_e-1.d0) * rho * Te0 * cs0   * BigR  * dl * c_angle * theta * tstep 
                    amat(var_Te,var_Te)   = + v * (gamma_sheath_e-1.d0) * r0  * Te  * fx_n * theta * tstep &
                                            + v * (gamma_sheath_e-1.d0) * r0  * Te  * cs0   * BigR  * dl * c_angle * theta * tstep &
                                            + v * (gamma_sheath_e-1.d0) * r0  * Te0 * cs_Te * BigR  * dl * c_angle * theta * tstep

                    amat(var_Ti,var_vpar) = + v * (gamma_sheath_i-1.d0) * r0  * Ti0 * vpar * fx_v * theta * tstep &
                                            + v * (GAMMA - 1.d0) * vpar * visco_par_heating * gradvpar0dotn * BigR * dl    * theta * tstep &
                                            + v * (GAMMA - 1.d0) * vpar0 * visco_par_heating * gradvpardotn * BigR * dl    * theta * tstep
                    amat(var_Te,var_vpar) = + v * (gamma_sheath_e-1.d0) * r0  * Te0 * vpar * fx_v * theta * tstep 
                    amat(var_Ti,var_u)    = + v * (gamma_sheath_i-1.d0) * r0  * Ti0 * fx_u * psi_s * theta * tstep
                    amat(var_Te,var_u)    = + v * (gamma_sheath_e-1.d0) * r0  * Te0 * fx_u * psi_s * theta * tstep
                  else
                    amat(var_T,var_psi)   = + v * (gamma_sheath  -1.d0) * r0  *  T0 * psi_s * fx_p * theta * tstep 
                    amat(var_T,var_rho)   = + v * (gamma_sheath  -1.d0) * rho *  T0 * fx_n * theta * tstep &
                                            + v * (gamma_sheath  -1.d0) * rho *  T0 * cs0   * BigR  * dl * c_angle * theta * tstep 
                    amat(var_T,var_T)     = + v * (gamma_sheath  -1.d0) * r0  *  T  * fx_n * theta * tstep &
                                            + v * (gamma_sheath  -1.d0) * r0  *  T  * cs0   * BigR  * dl * c_angle * theta * tstep &
                                            + v * (gamma_sheath  -1.d0) * r0  *  T0 * cs_T  * BigR  * dl * c_angle * theta * tstep

                    amat(var_T,var_u)     = + v * (gamma_sheath  -1.d0) * r0  * T0  * fx_u * psi_s * theta * tstep
                    amat(var_T,var_vpar)  = + v * (gamma_sheath  -1.d0) * r0  * T0  * vpar * fx_v * theta * tstep & 
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

                  ! --- Weak Bohm condition on the total normal flow (mach1_weak): exact columns of mw_res
                  amat(var_vpar,var_vpar) = amat(var_vpar,var_vpar) + v * mw_w * mw_Bn * vpar
                  if (with_TiTe) then
                    amat(var_vpar,var_Ti) = amat(var_vpar,var_Ti) - v * mw_w * mw_cs * abs(bdotn) * cs_Ti
                    amat(var_vpar,var_Te) = amat(var_vpar,var_Te) - v * mw_w * mw_cs * abs(bdotn) * cs_Te
                  else
                    amat(var_vpar,var_T)  = amat(var_vpar,var_T)  - v * mw_w * mw_cs * abs(bdotn) * cs_T
                  endif
                  amat(var_vpar,var_u)    = - v * mw_w * mw_act * mw_orient * BigR * psi_s / dl
                  ! --- No psi column: B_pol.n depends on the tangential psi derivative only, which is Dirichlet on
                  ! --- the wall. |B| also sees the free normal derivative, but this loop only carries the trace
                  ! --- DOFs, and |B| is F0/R to O(B_pol^2/F0^2) at the wall, so that dependence is lagged.

                  ! --- Fluid neutral sources and reflection
                  if (with_neutrals) then
                    amat(var_rhon,var_psi) = - v * neutral_reflection * r0  * psi_s * fx_p      * theta * tstep 
  
                    amat(var_rhon,var_rho) = - v * neutral_reflection * rho     * fx_n      * theta * tstep &
                                             - v * neutral_reflection * rho     * cs0 * BigR * dl * c_angle * theta * tstep 
 
                    if (with_TiTe) then 
                      amat(var_rhon,var_Ti) = - v * neutral_reflection * r0 * cs_Ti * BigR * dl * c_angle * theta * tstep 
                      amat(var_rhon,var_Te) = - v * neutral_reflection * r0 * cs_Te * BigR * dl * c_angle * theta * tstep 
                    else
                      amat(var_rhon,var_T)  = - v * neutral_reflection * r0 * cs_T  * BigR * dl * c_angle * theta * tstep 
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

          ! --- Surface term of the current definition, columns of the psi normal-derivative DOFs (3 and 4 of
          ! --- both edge nodes). Their trace vanishes, their t-derivative on the edge is what the state loop
          ! --- above uses: size(vertex,direction(l))*H1*HZ*element_size_perp with the same sign rule.
          if ( sj_surf ) then
            do k=1,2
              sj_esp = - element%size(vertex(k),direction_perp(1)) * 3.d0
              if ((vertex(1)*vertex(2) .eq. 2)) sj_esp = + element%size(vertex(k),direction_perp(1)) * 3.d0
              do l=1,2
                l3 = direction_perp(l)
                do in = i_tor_min, i_tor_max
                  sj_Tt = H1(k,l,ms) * element%size(vertex(k),direction(l)) * HZ(in,mp) * sj_esp
                  index_kl = n_tor_local*n_var*n_degrees*(vertex(k)-1) + n_tor_local * n_var * (l3-1) + in - i_tor_min +1
                  ELM(index_ij+(var_zj-1)*(n_tor_local),index_kl+(var_psi-1)*(n_tor_local)) = &
                  ELM(index_ij+(var_zj-1)*(n_tor_local),index_kl+(var_psi-1)*(n_tor_local))   &
                    - v * sj_w * ( - y_s(ms)*normal(1) + x_s(ms)*normal(2) ) / xjac * sj_Tt * ws
                enddo
              enddo
            enddo
          endif

          ! --- Cancelled wall fluxes: columns of the normal-derivative DOFs of w (viscous flux through dw/dn) and
          ! --- of rho, Ti, Te (magnetisation flux through p_R, p_Z); same DOF rule as the surface term above.
          if ( sj_cf ) then
            do k=1,2
              sj_esp = - element%size(vertex(k),direction_perp(1)) * 3.d0
              if ((vertex(1)*vertex(2) .eq. 2)) sj_esp = + element%size(vertex(k),direction_perp(1)) * 3.d0
              do l=1,2
                l3 = direction_perp(l)
                do in = i_tor_min, i_tor_max
                  sj_Tt = H1(k,l,ms) * element%size(vertex(k),direction(l)) * HZ(in,mp) * sj_esp
                  cf_tx = cf_dtx * sj_Tt ; cf_ty = cf_dty * sj_Tt          ! (R,Z) derivatives of this basis function on the edge
                  index_kl = n_tor_local*n_var*n_degrees*(vertex(k)-1) + n_tor_local * n_var * (l3-1) + in - i_tor_min +1
                  ELM(index_ij+(var_u-1)*(n_tor_local),index_kl+(var_w-1)*(n_tor_local)) = &
                  ELM(index_ij+(var_u-1)*(n_tor_local),index_kl+(var_w-1)*(n_tor_local))   &
                    - v * cf_nu * BigR * cf_fo * BigR**2 * ( cf_tx*normal(1) + cf_ty*normal(2) ) * dl * theta * tstep * ws
                  ELM(index_ij+(var_u-1)*(n_tor_local),index_kl+(var_rho-1)*(n_tor_local)) = &
                  ELM(index_ij+(var_u-1)*(n_tor_local),index_kl+(var_rho-1)*(n_tor_local))   &
                    + v * BigR**2 * ( - normal(2) * cf_tx + normal(1) * cf_ty ) * T0 * dl * theta * tstep * ws
                  if (with_TiTe) then
                    ELM(index_ij+(var_u-1)*(n_tor_local),index_kl+(var_Ti-1)*(n_tor_local)) = &
                    ELM(index_ij+(var_u-1)*(n_tor_local),index_kl+(var_Ti-1)*(n_tor_local))   &
                      + v * BigR**2 * ( - normal(2) * cf_tx + normal(1) * cf_ty ) * r0 * dl * theta * tstep * ws
                    ELM(index_ij+(var_u-1)*(n_tor_local),index_kl+(var_Te-1)*(n_tor_local)) = &
                    ELM(index_ij+(var_u-1)*(n_tor_local),index_kl+(var_Te-1)*(n_tor_local))   &
                      + v * BigR**2 * ( - normal(2) * cf_tx + normal(1) * cf_ty ) * r0 * dl * theta * tstep * ws
                  else
                    ELM(index_ij+(var_u-1)*(n_tor_local),index_kl+(var_T-1)*(n_tor_local)) = &
                    ELM(index_ij+(var_u-1)*(n_tor_local),index_kl+(var_T-1)*(n_tor_local))   &
                      + v * BigR**2 * ( - normal(2) * cf_tx + normal(1) * cf_ty ) * r0 * dl * theta * tstep * ws
                  endif
                enddo
              enddo
            enddo
          endif

        enddo
      enddo
    enddo

  enddo
enddo

return
end subroutine

end module mod_boundary_matrix_open
