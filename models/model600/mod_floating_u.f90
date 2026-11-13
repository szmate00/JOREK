!> Sheath wall conditions for model600: the floating-potential row (bcs%floating_u) and the sheath current as the
!! wall current of the vorticity row (bcs%sheath_u).
!!
!! Floating row: the zero-local-current limit of the sheath characteristic on the potential at a material wall,
!!
!!     Phi - V_wall = Lambda * k_B*Te / e   <=>   u = C_T * Te + C_V * V_wall ,
!!
!! with Lambda = sheath_Lambda, V_wall = sheath_V_wall. Affine in JOREK variables, imposed exactly by the linear solve.
!!
!! Sheath current (bcs%sheath_u): on the value DOF of u the vorticity row is kept with its parallel term integrated
!! by parts, and the wall integral of the normal current is set to the sheath current (mod_boundary_matrix_open):
!!
!!     J_par b_n + J_perp.n = J_sh(Phi) = e n c_s |b_n| [1 - exp(Lambda - e(Phi - V_wall)/Te)] .
!!
!! zj at the wall is then determined by the induction row (Ohm's law), which needs dirichlet%zj = .f. on the type
!! (the dpsi = 0 penalty moves to the zj row, mod_boundary_conditions). The slope DOFs of u keep the floating row.
!!
!! Sign convention: Phi = +F0*u. JOREK's (R,Z,phi) basis is right-handed and the poloidal ExB velocity the fluid
!! advects with is v = (-R*u_Z, +R*u_R) = -R grad(u) x e_phi, Hoelzl et al. 2021 eq. 26 with u = Phi/F0. So a_n
!! below carries the sign of F0 and the physical potential does not depend on the field direction.
module mod_floating_u

  implicit none
  private

  public :: floating_u_norm, floating_u_volts, floating_u_selftest, sheath_j_norm, wall_node_bn, sheath_u_edge

contains


!> Normalisation of the floating condition in JOREK units: u = C_T*Te + C_V*V_wall.
!!
!! From Phi = F0*u/sqrt(mu0*rho0) and k_B*Te [J] = Te_JOREK/(mu0*n0):
!!   a_n = 2*e*F0*sqrt(mu0*rho0)/m_i   so that   e*Phi/(k_B*Te) = a_n*u/(2*Te_JOREK)
!!   C_T = 2*Lambda/a_n                (halved in a single-temperature build, where Te = T/2)
!!   C_V = sqrt(mu0*rho0)/F0           (u per volt)
pure subroutine floating_u_norm(a_n, C_T, C_V)

  use constants,          only: MU_ZERO, ATOMIC_MASS_UNIT, EL_CHG
  use phys_module,        only: F0, central_density, central_mass, sheath_Lambda
  use mod_model_settings, only: with_TiTe

  implicit none
  real*8, intent(out) :: a_n, C_T, C_V

  real*8 :: m_i, rho0, lam

  m_i  = central_mass * ATOMIC_MASS_UNIT
  rho0 = central_density * 1.d20 * m_i

  a_n  = 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * rho0) / m_i

  lam  = sheath_Lambda
  if ( .not. with_TiTe ) lam = 0.5d0 * lam   ! single-T build evolves T = Ti + Te, Te = T/2

  C_T  = 2.d0 * lam / a_n
  C_V  = sqrt(MU_ZERO * rho0) / F0

end subroutine floating_u_norm


!> Physical potential in volts from u.
pure real*8 function floating_u_volts(u)

  use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT
  use phys_module, only: F0, central_density, central_mass

  implicit none
  real*8, intent(in) :: u

  floating_u_volts = F0 * u / sqrt(MU_ZERO * central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT)

end function floating_u_volts


!> Self-test of the normalisation at the current namelist values: a_n carries the sign of F0 and Te = 1 eV
!! reconstructs to exactly Lambda volts above the wall.
logical function floating_u_selftest(my_id)

  use constants,          only: MU_ZERO, EL_CHG
  use phys_module,        only: F0, central_density, sheath_Lambda
  use mod_model_settings, only: with_TiTe

  implicit none
  integer, intent(in) :: my_id

  real*8, parameter :: tol = 1.d-12
  real*8 :: a_n, C_T, C_V, Te_1eV, u_1eV, volts

  call floating_u_norm(a_n, C_T, C_V)

  Te_1eV = EL_CHG * MU_ZERO * central_density * 1.d20      ! 1 eV in JOREK temperature units
  u_1eV  = C_T * Te_1eV
  if ( .not. with_TiTe ) u_1eV = C_T * 2.d0 * Te_1eV       ! trace variable is T = 2*Te
  volts  = floating_u_volts(u_1eV)

  floating_u_selftest = ( a_n * F0 .gt. 0.d0 ) .and. &
                        ( abs(volts - sheath_Lambda) .le. tol * max(abs(sheath_Lambda), 1.d0) )

  if ( my_id .eq. 0 .and. .not. floating_u_selftest ) &
    write(*,'(A,2es14.6)') ' ERROR: floating_u normalisation selftest failed (a_n, volts at 1 eV):', a_n, volts

end function floating_u_selftest


!> Normalisation of the sheath current (Artola eqs. 5, 7, 8): the ion saturation current in the toroidal-current
!! variable is j_sat = c_sat*rho*Vpar with c_sat = -e*F0*n0*sqrt(mu0/rho0), and the exponent of the characteristic
!! is e*Phi/(k_B*Te) = a_n*u/(2*Te_JOREK). c_sat carries the sign of F0 like zj itself, so the current INTO the
!! wall, -zj*(B_pol.n)/F0, is independent of the field sign.
pure subroutine sheath_j_norm(a_n, c_sat)

  use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT, EL_CHG
  use phys_module, only: F0, central_density, central_mass

  implicit none
  real*8, intent(out) :: a_n, c_sat

  real*8 :: m_i, n_0, rho0

  m_i   = central_mass * ATOMIC_MASS_UNIT
  n_0   = central_density * 1.d20
  rho0  = n_0 * m_i

  a_n   =   2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * rho0) / m_i
  c_sat = - EL_CHG * F0 * n_0 * sqrt(MU_ZERO / rho0)

end subroutine sheath_j_norm


!> |b.n| = |B_pol.n|/|B| at a wall node along one of its wall edges, from the node's own DOFs (lagged psi, which is
!! Dirichlet on the wall, so this is a static map). Shared by mod_boundary_conditions (release of the u value DOF)
!! and mod_boundary_matrix_open (assembly of the sheath edge term) so that both decide from the same number.
!!   xn(1:3,1:2)  node value, s- and t-derivative DOFs of R (:,1) and Z (:,2)
!!   psi_s, psi_t the s- and t-derivative DOFs of psi at the node
!!   size_s, size_t  element%size(iv,2), element%size(iv,3) of the element the edge belongs to
!!   iv           vertex index of the node in that element
!!   s_const      the edge runs along t (s = const), else along s
pure real*8 function wall_node_bn(xn, psi_s, psi_t, size_s, size_t, iv, s_const)

  use phys_module,        only: F0
  use mod_basisfunctions, only: basisfunctions1
  use mod_settings,       only: n_degrees_1d

  implicit none
  real*8,  intent(in) :: xn(3,2), psi_s, psi_t, size_s, size_t
  integer, intent(in) :: iv
  logical, intent(in) :: s_const

  real*8 :: Hb(2,n_degrees_1d), Hb_s(2,n_degrees_1d), Hb_ss(2,n_degrees_1d)
  real*8 :: es_s, es_t, p_s, p_t, r_s, r_t, z_s, z_t, xj, p_x, p_y, g_b(2), bt, rr

  call basisfunctions1(0.d0, Hb, Hb_s, Hb_ss)
  es_s = size_s * Hb_s(1,2)
  es_t = size_t * Hb_s(1,2)
  if ((iv .eq. 2) .or. (iv .eq. 3)) es_s = - es_s
  if ((iv .eq. 3) .or. (iv .eq. 4)) es_t = - es_t
  p_s = psi_s * es_s
  p_t = psi_t * es_t
  rr  = xn(1,1)
  r_s = xn(2,1) * es_s ; r_t = xn(3,1) * es_t
  z_s = xn(2,2) * es_s ; z_t = xn(3,2) * es_t
  xj  = r_s*z_t - r_t*z_s
  p_x = (   z_t * p_s - z_s * p_t ) / xj
  p_y = ( - r_t * p_s + r_s * p_t ) / xj
  if ( s_const ) then
    g_b = (/  z_t, -r_t /) / xj          ! grad(s): normal to the edge s = const
  else
    g_b = (/ -z_s,  r_s /) / xj          ! grad(t)
  endif
  g_b = g_b / norm2(g_b)
  bt  = sqrt(F0**2 + p_x**2 + p_y**2) / rr
  wall_node_bn = abs( p_y*g_b(1) - p_x*g_b(2) ) / (rr * bt)

end function wall_node_bn


!> Does a wall edge carry the sheath edge term? Both end nodes of sheath_u types and the incidence above
!! sin(min_sheath_angle) at both ends. The release of the u value DOF (mod_boundary_conditions) requires this on
!! EVERY wall edge of the node: an edge of a released node that kept the strong form would add an indefinite
!! operator to the potential row (doc/sheath_option_a, Sec. 3).
pure logical function sheath_u_edge(btype1, btype2, bn1, bn2)

  use constants,   only: PI
  use phys_module, only: bcs, min_sheath_angle, max_bnd_types

  implicit none
  integer, intent(in) :: btype1, btype2
  real*8,  intent(in) :: bn1, bn2

  sheath_u_edge = .false.
  if ( btype1 .lt. 1 .or. btype1 .gt. max_bnd_types ) return
  if ( btype2 .lt. 1 .or. btype2 .gt. max_bnd_types ) return
  sheath_u_edge = bcs(btype1)%sheath_u .and. bcs(btype2)%sheath_u .and. &
                  ( min(bn1, bn2) .ge. sin(min_sheath_angle * PI / 180.d0) )

end function sheath_u_edge


end module mod_floating_u
