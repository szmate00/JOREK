!> Sheath current-voltage characteristic for the electric potential boundary condition (model600).
!!
!! This module is the single place where the sheath physics lives. It is used by
!!  - mod_boundary_matrix_open.f90 : the charge-continuity surface term (the boundary condition)
!!  - mod_sheath_diag.f90          : the wall current / potential diagnostic
!!  - mod_boundary_conditions.f90  : the older nodal variant, kept for comparison
!! so that the three can never drift apart.
!!
!! It is deliberately stateless: every routine derives its constants from phys_module on the fly.
!! That costs a square root and a logarithm per Gauss point, which is nothing next to the exp()
!! and the surrounding assembly, and it removes any dependence on an initialisation order and any
!! risk of a race when called from an OpenMP region.
!!
!! PHYSICS (Stangeby 2.68 / Langmuir probe characteristic, referenced to the wall):
!!
!!     j = j_sat * ( 1 - exp(-X) ) ,   X = e*Phi/(k*Te) - Lambda ,   Phi = V_sheath - V_wall
!!
!! j = 0 at X = 0, i.e. Phi = Lambda*Te/e (the floating potential); Phi = 0 (X = -Lambda) is
!! electron saturation. The characteristic is evaluated in the FORWARD direction, current as a
!! function of potential. That direction is single valued and its derivative dj/dPhi ~ exp(-X)
!! tends smoothly to zero once the ion saturation current is reached, so no clipping is needed on
!! the ion side and the linearisation never becomes singular. Inverting it for Phi (as the older
!! nodal implementation does) is singular exactly at j -> j_sat, which is where divertor
!! conditions actually sit.
!!
!! JOREK UNITS: the electrostatic potential is Phi = -F0*u in the code's variables. The JOREK
!! reference paper (Hoelzl et al 2021, eq. 26) defines u = Phi/F0 together with a velocity ansatz
!! v_pol = -R grad(u) x e_phi, whereas the model600 element matrix implements
!! v_pol = +R grad(u) x e_phi (see the ExB advection and compression terms of the density
!! equation, and the [u,psi] term of the induction equation). The code's u is therefore minus the
!! paper's, which is where the minus sign in a_n comes from. The same holds for psi and for the
!! current variable: the code's zj = Delta*psi is minus the paper's j.
!!
!!     e*Phi/(k*Te) = ( a_n*u/2 - vw ) / Te
!!     a_n          = -2*e*F0*sqrt(mu0*rho0)/m_i
!!     vw           = e*V_wall*mu0*n0
!!     c_sat        = -a_n/2
!!     zj_sat       = c_sat * rho * g_bn * cs / |B|      (Artola eqs. 5 and 16, in the code's sign)
!!
!! zj_sat uses the same parallel velocity g(b_n)*cs/|B| that the Mach 1 boundary condition
!! imposes, i.e. the Chodura-Riemann smoothing function including its sign. It corresponds to a
!! force-free (field aligned) SOL current, J_par = -j*B/F0 in the paper's variables, which is the
!! physically appropriate interpretation in the scrape-off layer.
!! HOW TO ENABLE IT (namelist in1), on the boundary types that carry the strike points:
!!
!!   bcs(1)%natural%u    = .true.    ! the sheath boundary condition itself
!!   bcs(1)%natural%w    = .true.    ! required: keeps grad(u).n free, i.e. keeps E_r alive
!!   bcs(1)%natural%zj   = .true.    ! required: lets the boundary current follow Ampere's law
!!   bcs(1)%dirichlet%u  = .false.
!!   bcs(1)%dirichlet%w  = .false.
!!   bcs(1)%dirichlet%zj = .false.
!!   bcs(1)%mach1        = .true.    ! keep it: j_sat assumes the Mach 1 condition at the same wall
!!   bc_natural_open     = .true.    ! the boundary integrals live in that branch of construct_matrix
!!
!! Leave at least one boundary type with dirichlet%u = .true. (typically the main chamber wall):
!! u enters the vorticity equation only through its gradient, and the sheath term loses its grip
!! on u in ion saturation, so its constant mode would otherwise be undetermined. Optional knobs:
!! sheath_Lambda (Lambda_0, <=0 computes it from central_mass), sheath_Lambda_local,
!! sheath_V_wall, sheath_X_min, sheath_smooth_dX, sheath_min_bn, sheath_ramp_time.
!!
module mod_sheath_bc

  implicit none

  private

  public :: sheath_norm, sheath_get_lambda, sheath_x_limited, sheath_current

  !> Below these values the characteristic is evaluated with the floor instead. Callers are
  !! expected to pass corr_neg corrected quantities; this is only a last defence against a
  !! division by zero in a debug build with -ffpe-trap=zero.
  real*8, parameter :: sheath_t_floor   = 1.d-14
  real*8, parameter :: sheath_rho_floor = 1.d-10

  !> exp() is only evaluated inside this range; outside it the asymptotic branch is used, which
  !! keeps the debug build (-ffpe-trap=overflow,denormal) alive for any state the solver may pass
  !! through on its way to the solution.
  real*8, parameter :: sheath_exp_max = 30.d0

contains

!> Normalisation constants of the characteristic in JOREK units.
!! @param a_n    -2*e*F0*sqrt(mu0*rho0)/m_i , so that e*Phi/(k*Te) = (a_n*u/2 - vw)/Te
!! @param c_sat  -a_n/2, the ion saturation current coefficient
!! @param vw     e*V_wall*mu0*n0, the wall bias in the same units
pure subroutine sheath_norm(a_n, c_sat, vw)

  use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT, EL_CHG
  use phys_module, only: F0, central_density, central_mass, sheath_V_wall

  implicit none
  real*8, intent(out) :: a_n, c_sat, vw

  real*8 :: m_i, rho0

  m_i  = central_mass * ATOMIC_MASS_UNIT
  rho0 = central_density * 1.d20 * m_i

  ! --- NOTE the minus sign: Phi = -F0*u in the code's variables (see the module header)
  a_n   = - 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * rho0) / m_i
  c_sat = - 0.5d0 * a_n
  vw    =   EL_CHG * sheath_V_wall * MU_ZERO * central_density * 1.d20

end subroutine sheath_norm


!> Sheath factor Lambda and its derivatives.
!!
!!   Lambda = ln( (1/4)*v_th,e / c_s ) = Lambda_0 - ln sqrt( gamma*(1 + Ti/Te) )
!!
!! Lambda_0 = ln sqrt(m_i/(2*pi*m_e)) is taken from sheath_Lambda when that is positive, and
!! computed from central_mass otherwise. The Ti/Te correction is not cosmetic: with
!! c_s = sqrt(gamma*(Ti+Te)) and Ti = Te it is -ln sqrt(2*gamma) = -0.6, i.e. Lambda = 2.4 rather
!! than 3.0, which is 12 V at Te = 20 eV. Set sheath_Lambda_local = .false. to freeze Lambda_0.
pure subroutine sheath_get_lambda(Ti, Te, lam, dlam_dTi, dlam_dTe)

  use constants,   only: PI, ATOMIC_MASS_UNIT, MASS_ELECTRON
  use phys_module, only: GAMMA, central_mass, sheath_Lambda, sheath_Lambda_local

  implicit none
  real*8, intent(in)  :: Ti, Te
  real*8, intent(out) :: lam, dlam_dTi, dlam_dTe

  real*8 :: Ti_l, Te_l, T_l, lam0

  if ( sheath_Lambda .gt. 0.d0 ) then
    lam0 = sheath_Lambda
  else
    lam0 = log( sqrt( central_mass * ATOMIC_MASS_UNIT / (2.d0*PI*MASS_ELECTRON) ) )
  endif

  if ( .not. sheath_Lambda_local ) then
    lam      = lam0
    dlam_dTi = 0.d0
    dlam_dTe = 0.d0
    return
  endif

  Ti_l = max(Ti, sheath_t_floor)
  Te_l = max(Te, sheath_t_floor)
  T_l  = Ti_l + Te_l

  lam      = lam0 - 0.5d0 * log( GAMMA * T_l / Te_l )
  dlam_dTi = - 0.5d0 / T_l
  dlam_dTe =   0.5d0 * Ti_l / (Te_l * T_l)

end subroutine sheath_get_lambda


!> Smooth lower limit on the sheath exponent X (the electron saturation side).
!!
!!   X_lim = X_min + dx*ln(1 + exp((X-X_min)/dx))
!!
!! C1 continuous, monotone, dX_lim/dX in (0,1), and X_lim -> X for X >> X_min. The ion side needs
!! no limit at all: 1-exp(-X) <= 1 by construction, which is the whole point of using the forward
!! form of the characteristic.
pure subroutine sheath_x_limited(x, x_lim, dxlim_dx)

  use phys_module, only: sheath_X_min, sheath_smooth_dX

  implicit none
  real*8, intent(in)  :: x
  real*8, intent(out) :: x_lim, dxlim_dx

  real*8 :: z, dx

  dx = max(sheath_smooth_dX, 1.d-3)
  z  = (x - sheath_X_min) / dx

  if ( z .gt. sheath_exp_max ) then          ! exp(z) would overflow; X_lim = X to 1e-13
    x_lim    = x
    dxlim_dx = 1.d0
  elseif ( z .lt. -sheath_exp_max ) then     ! exp(z) would underflow towards a denormal
    x_lim    = sheath_X_min
    dxlim_dx = 0.d0
  else
    x_lim    = sheath_X_min + dx * log( 1.d0 + exp(z) )
    dxlim_dx = 1.d0 / ( 1.d0 + exp(-z) )
  endif

end subroutine sheath_x_limited


!> The sheath current in JOREK current-variable (zj) units, with its exact derivatives.
!!
!! @param u      electric potential variable at the wall (JOREK units, Phi = -F0*u)
!! @param rho    density (JOREK units, corr_neg corrected by the caller)
!! @param Ti,Te  ion/electron temperature (JOREK units, corr_neg corrected by the caller)
!! @param g_bn   Chodura-Riemann function g(b_n) INCLUDING its sign (= normal_sign*factor)
!! @param sgn_bn sign of B.n, giving the grazing incidence floor a well defined direction
!! @param Btot   |B|
!! @param zj_sh  sheath current, expressed as the zj value that carries it
!! @param zj_sat ion saturation current in the same units (diagnostic)
!! @param x_out  the exponent X actually used, after limiting (diagnostic)
subroutine sheath_current(u, rho, Ti, Te, g_bn, sgn_bn, Btot,        &
                          zj_sh, dzj_du, dzj_drho, dzj_dTi, dzj_dTe, &
                          zj_sat, x_out)

  use phys_module, only: GAMMA, sheath_min_bn

  implicit none

  real*8, intent(in)  :: u, rho, Ti, Te, g_bn, sgn_bn, Btot
  real*8, intent(out) :: zj_sh, dzj_du, dzj_drho, dzj_dTi, dzj_dTe
  real*8, intent(out) :: zj_sat, x_out

  real*8 :: a_n, c_sat, vw
  real*8 :: rho_l, Ti_l, Te_l, T_l, cs, g_eff
  real*8 :: lam, dlam_dTi, dlam_dTe
  real*8 :: x, x_lim, dxlim_dx, expmx, f, fp
  real*8 :: dx_du, dx_dTi, dx_dTe

  call sheath_norm(a_n, c_sat, vw)

  rho_l = max(rho, sheath_rho_floor)
  Ti_l  = max(Ti,  sheath_t_floor)
  Te_l  = max(Te,  sheath_t_floor)
  T_l   = Ti_l + Te_l

  cs = sqrt( GAMMA * T_l )

  ! --- Grazing incidence: g(b_n) vanishes where the field is tangent to the wall, which would
  ! --- switch the sheath off there, and with it the only term that ties u to the current. Floor
  ! --- the magnitude while keeping the direction of B.n, in the spirit of the minimum particle
  ! --- and heat fluxes already applied in mod_boundary_matrix_open.
  g_eff = g_bn
  if ( sheath_min_bn .gt. 0.d0 ) g_eff = sgn_bn * max( abs(g_bn), sheath_min_bn )

  zj_sat = c_sat * rho_l * g_eff * cs / Btot

  call sheath_get_lambda(Ti_l, Te_l, lam, dlam_dTi, dlam_dTe)

  ! --- X = e*Phi/(k*Te) - Lambda
  x = ( 0.5d0 * a_n * u - vw ) / Te_l - lam

  call sheath_x_limited(x, x_lim, dxlim_dx)
  x_out = x_lim

  if ( x_lim .lt. -sheath_exp_max ) then
    expmx = exp(sheath_exp_max)     ! electron branch, capped. Only reachable when the user has
                                    ! effectively disabled the limiter with a very negative X_min
  elseif ( x_lim .gt. sheath_exp_max ) then
    expmx = 0.d0
  else
    expmx = exp(-x_lim)
  endif

  f  = 1.d0 - expmx                 ! the characteristic
  fp = expmx * dxlim_dx             ! d f / d X

  ! --- derivatives of X (Lambda contributes through Ti and Te when sheath_Lambda_local)
  dx_du  =   0.5d0 * a_n / Te_l
  dx_dTi = - dlam_dTi
  dx_dTe = - ( x + lam ) / Te_l - dlam_dTe

  zj_sh    = zj_sat * f
  dzj_du   = zj_sat * fp * dx_du
  dzj_drho = zj_sat * f / rho_l
  dzj_dTi  = zj_sat * f / (2.d0 * T_l) + zj_sat * fp * dx_dTi
  dzj_dTe  = zj_sat * f / (2.d0 * T_l) + zj_sat * fp * dx_dTe

  ! --- the floors are hard clips, so drop the corresponding sensitivity to stay consistent
  if ( rho .lt. sheath_rho_floor ) dzj_drho = 0.d0
  if ( Ti  .lt. sheath_t_floor   ) dzj_dTi  = 0.d0
  if ( Te  .lt. sheath_t_floor   ) dzj_dTe  = 0.d0

end subroutine sheath_current

end module mod_sheath_bc
