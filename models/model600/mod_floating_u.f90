!> Floating-potential boundary condition for model600 (bcs(i)%floating_u).
!!
!! Imposes the zero-local-current limit of the sheath characteristic on the plasma potential at a
!! material wall,
!!
!!     Phi - V_wall = Lambda * k_B*Te / e ,
!!
!! with Lambda = sheath_Lambda and V_wall = sheath_V_wall. In JOREK variables this is affine,
!!
!!     u = C_T * Te + C_V * V_wall ,
!!
!! so the boundary row has a constant Jacobian and is imposed exactly by the linear solve. It carries
!! no net wall current, hence no thermoelectric current; what it does give is a wall potential that
!! follows Te, i.e. a tangential electric field and an ExB drift wherever Te varies along the wall.
!!
!! Sign convention: Phi = +F0*u. JOREK's (R,Z,phi) basis is right-handed and the poloidal ExB
!! velocity the fluid advects with is v = (-R*u_Z, +R*u_R) = -R grad(u) x e_phi (mod_elt_matrix_fft,
!! particles/mod_fields), which is Hoelzl et al. 2021 eq. 26 with u = Phi/F0. So a_n below carries
!! the sign of F0 and the physical potential does not depend on the field direction.
module mod_floating_u

  implicit none
  private

  public :: floating_u_norm, floating_u_volts, floating_u_selftest, sheath_j_norm, sheath_j_ramp

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


!> Physical potential in volts from u. One place, shared by the row and the diagnostics.
pure real*8 function floating_u_volts(u)

  use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT
  use phys_module, only: F0, central_density, central_mass

  implicit none
  real*8, intent(in) :: u

  floating_u_volts = F0 * u / sqrt(MU_ZERO * central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT)

end function floating_u_volts


!> Self-test of the normalisation at the current namelist values. Checks that a_n carries the sign
!! of F0 and that Te = 1 eV reconstructs to exactly Lambda volts above the wall.
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


!> Normalisation of the sheath current row (Artola eqs. 5, 7, 8): the ion saturation current in the
!! toroidal-current variable is j_sat = c_sat*rho*Vpar with c_sat = -e*F0*n0*sqrt(mu0/rho0), and the
!! exponent of the characteristic is e*Phi/(k_B*Te) = a_n*u/(2*Te_JOREK). c_sat carries the sign of F0
!! like zj itself, so the current INTO the wall, -zj*(B_pol.n)/F0, is independent of the field sign.
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


!> Switch-on ramp alpha(t) of the sheath characteristic zj = j_sat*(1 - exp(x/alpha)): alpha = 1 is the
!! characteristic itself; a small alpha is the sheath of a plasma at alpha*Te, a stiff I-V that pins the potential
!! near floating and lets any current pass, so the equilibrium's wall current (which is no sheath current) can be
!! relaxed by the induction equation while the sheath tightens, as JOREK ramps RMPs. Linear from
!! sheath_j_ramp_alpha0 at t = 0 to 1 at t_end: sheath_j_ramp_time > 0 sets t_end, 0 takes the end of the
!! timestep ramp (the last tstep_n phase begins), < 0 disables the ramp. The end state is the same either way.
pure real*8 function sheath_j_ramp(t)

  use phys_module, only: sheath_j_ramp_time, sheath_j_ramp_alpha0, tstep_n, nstep_n

  implicit none
  real*8, intent(in) :: t

  real*8  :: t_end
  integer :: i, n_last

  sheath_j_ramp = 1.d0
  if ( sheath_j_ramp_time .lt. 0.d0 ) return
  t_end = sheath_j_ramp_time
  if ( t_end .eq. 0.d0 ) then
    n_last = 1
    do i = 1, size(nstep_n)
      if ( nstep_n(i) .gt. 0 ) n_last = i
    enddo
    do i = 1, n_last-1
      t_end = t_end + tstep_n(i) * nstep_n(i)
    enddo
  endif
  if ( t_end .le. 0.d0 ) return
  sheath_j_ramp = sheath_j_ramp_alpha0 + ( 1.d0 - sheath_j_ramp_alpha0 ) * min( 1.d0, max( 0.d0, t / t_end ) )

end function sheath_j_ramp


end module mod_floating_u
