!> Prescribed floating-potential boundary condition for model600.
!!
!! Imposes the LOCAL ZERO-CURRENT limit of the sheath characteristic,
!!
!!     V_p - V_wall = Lambda * k_B*Te / e ,     Lambda = sheath_Lambda (3 by default)
!!
!! on the plasma potential at a material wall. Full derivation, with every step referenced to a
!! file and line of the source: doc/floating_u_derivation.tex (and .pdf).
!!
!! WHAT THIS IS NOT. It cannot carry a net wall current, so it produces no thermoelectric target
!! current and it is not a model of a conducting vessel closing current through the wall. What it
!! does give - and what the full characteristic also gives - is a wall potential that follows Te,
!! hence a tangential electric field and an ExB drift wherever Te varies along the wall.
!!
!! WHY IT IS ROBUST WHERE THE CURRENT-CARRYING ROUTE IS NOT. In JOREK variables the condition is
!! AFFINE in the evolved variables:
!!
!!     u = C_T*Te + C_V*V_wall
!!
!! so the boundary row has the constant Jacobian (1, -C_T), needs no Newton iteration, and has no
!! saturation, no inversion and no vanishing derivative. In particular it never forms B.n, never
!! divides by the element Jacobian, and never converts logical derivatives into (R,Z): the u and Te
!! trace DOFs are stored in the SAME nodal frame, so whatever that frame is - including a
!! near-degenerate one - its scaling cancels identically between the two terms. That is why the
!! same closure is valid on boundary types 1, 4, 5 and 9 alike, and why none of detmin, min_bn,
!! saturation slopes, ratio gates, current clips, ramps or relaxation parameters appears here.
!!
!! THE SIGN. sheath_norm in mod_sheath_bc.f90 returns a_n < 0 and MUST NOT BE REUSED - it encodes
!! Phi = -F0*u, which is wrong. Two comments in the current source disagree
!! (mod_sheath_bc.f90:207 says -F0*u, mod_boundary_conditions.f90:876 says +F0*u), so neither can
!! be trusted; the convention is settled from the coordinate basis instead. JOREK uses
!! (e_R, e_Z, e_phi) RIGHT-HANDED with phi clockwise from above, established three ways that are
!! not comments:
!!   * diagnostics/jorek2_povray.f90:102-103 and vacuum/mod_vacuum_fields.f90:1172-1173 map
!!     x = R cos(phi), y = -R sin(phi);
!!   * particles/mod_fields.f90:590-596 rot_tmp is the right-handed curl for (R,Z,phi) and applied
!!     to the code's own A returns the code's own B exactly (in the standard basis it returns -B);
!!   * diagnostics/new_diag/mod_expression.f90:1354 has Jtor = -zj/R.
!! In that basis grad(u) x e_phi = (u_Z, -u_R), and the implemented poloidal velocity
!! (particles/mod_fields.f90:581-582) is v_pol = (-R*u_Z, +R*u_R) = -R grad(u) x e_phi, which is
!! exactly the ansatz of Hoelzl et al 2021 eq. 26. Hence u_code = u_paper = Phi/F0, i.e.
!!
!!     Phi = +F0*u    and    a_n = +2*e*F0*sqrt(mu0*rho0)/m_i  >  0 for F0 > 0.
!!
!! a_n carries the sign of F0, so reversing the toroidal field reverses u while leaving the
!! reconstructed physical potential unchanged - as it must, since the floating condition contains
!! no field direction. floating_u_selftest() checks exactly that.
module mod_floating_u

  implicit none
  private

  public :: floating_u_norm, floating_u_target, floating_u_volts, floating_u_selftest

contains

!> Normalisation of the prescribed floating condition, in JOREK units.
!!
!! Returns the three quantities that define  u = C_T*Te + C_V*V_wall  and nothing else. Derived
!! from V_p = F0*u/sqrt(mu0*rho0) and k_B*Te[J] = Te_JOREK/(mu0*n0); see the module header for the
!! sign and doc/floating_u_derivation.tex for the algebra.
!!
!! @param a_n  2*e*F0*sqrt(mu0*rho0)/m_i, POSITIVE for F0 > 0
!! @param C_T  2*Lambda/a_n, the u-per-unit-Te coefficient. Halved automatically in a
!!             single-temperature build, where the model evolves T = Ti+Te and Te = T/2.
!! @param C_V  sqrt(mu0*rho0)/F0, the u-per-volt-of-V_wall coefficient. Identically 2*e*mu0*n0/a_n;
!!             both forms are evaluated and compared by floating_u_selftest.
pure subroutine floating_u_norm(a_n, C_T, C_V)

  use constants,        only: MU_ZERO, ATOMIC_MASS_UNIT, EL_CHG
  use phys_module,      only: F0, central_density, central_mass, sheath_Lambda
  use mod_model_settings, only: with_TiTe

  implicit none
  real*8, intent(out) :: a_n, C_T, C_V

  real*8 :: m_i, n_0, rho0, lam

  m_i  = central_mass * ATOMIC_MASS_UNIT
  n_0  = central_density * 1.d20
  rho0 = n_0 * m_i

  ! --- NO leading minus. See the module header: Phi = +F0*u in JOREK's (R,Z,phi) basis.
  a_n  = 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * rho0) / m_i

  lam  = sheath_Lambda
  ! --- Single-temperature build evolves T = Ti + Te, so the trace variable is T and Te = T/2.
  ! --- Folding the factor in here keeps ONE coefficient for the caller and one place to be wrong.
  if ( .not. with_TiTe ) lam = 0.5d0 * lam

  C_T  = 2.d0 * lam / a_n
  C_V  = sqrt(MU_ZERO * rho0) / F0

end subroutine floating_u_norm


!> The prescribed boundary value of u and its exact derivative with respect to the temperature
!! trace variable.
!!
!! RAW temperature, deliberately. A positivity map here would make the relation nonlinear (so the
!! "exact, not linearised" property would be lost) and, worse, would be applied to Fourier
!! COEFFICIENTS at the assembly site, which is meaningless for n /= 0. Negative temperatures are
!! the global positivity scheme's responsibility.
!!
!! @param T_raw   the evolved temperature trace variable (Te with WITH_TiTe, otherwise T)
!! @param V_wall  wall bias in volts
!! @param u_b     prescribed u
!! @param du_dT   d(u_b)/d(T_raw), exact, including the positivity map's derivative
subroutine floating_u_target(T_raw, V_wall, u_b, du_dT)

  implicit none
  real*8, intent(in)  :: T_raw, V_wall
  real*8, intent(out) :: u_b, du_dT

  real*8 :: a_n, C_T, C_V

  call floating_u_norm(a_n, C_T, C_V)

  ! --- RAW T. No positivity map: see the note at the assembly site. The relation must stay
  ! --- affine and mode-decoupled, so du/dT is the constant C_T.
  u_b   = C_T * T_raw + C_V * V_wall
  du_dT = C_T

end subroutine floating_u_target


!> Reconstruct the physical potential in VOLTS from u. One place, shared by the diagnostics, so the
!! conversion used to impose the condition and the conversion used to check it cannot drift apart.
pure real*8 function floating_u_volts(u)

  use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT
  use phys_module, only: F0, central_density, central_mass

  implicit none
  real*8, intent(in) :: u
  real*8 :: rho0

  rho0 = central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT
  floating_u_volts = F0 * u / sqrt(MU_ZERO * rho0)

end function floating_u_volts


!> Self-test of the normalisation, run once at setup on rank 0. Cheap, and it is the only thing
!! standing between a sign error and a run whose ExB drift points the wrong way.
!!
!! Checks, at the CURRENT namelist's normalisation:
!!   1. the two independent forms of C_V agree;
!!   2. a_n carries the sign of F0;
!!   3. the round trip closes: Te of 1 eV must reconstruct to exactly Lambda volts.
!! It does NOT check the reversed-F0 case or a nonzero V_wall - those are covered by the
!! standalone driver, not here, and this docstring previously overclaimed them.
!! Returns .false. and prints on failure. It does not stop; the caller decides.
logical function floating_u_selftest(my_id)

  use constants,        only: MU_ZERO, ATOMIC_MASS_UNIT, EL_CHG
  use phys_module,      only: F0, central_density, central_mass, sheath_Lambda
  use mod_model_settings, only: with_TiTe

  implicit none
  integer, intent(in) :: my_id

  real*8, parameter :: tol = 1.d-12
  real*8 :: a_n, C_T, C_V, C_V_alt, n_0, rho0, Te_1eV, u_1eV, volts, lam_eff

  floating_u_selftest = .true.

  call floating_u_norm(a_n, C_T, C_V)
  n_0  = central_density * 1.d20
  rho0 = n_0 * central_mass * ATOMIC_MASS_UNIT

  ! --- 1. the two closed forms of C_V must agree
  C_V_alt = 2.d0 * EL_CHG * MU_ZERO * n_0 / a_n
  if ( abs(C_V - C_V_alt) .gt. tol * abs(C_V) ) then
    floating_u_selftest = .false.
    if (my_id .eq. 0) write(*,'(A,2es22.14)') &
      ' FLOATING_U SELFTEST FAIL: C_V forms disagree:', C_V, C_V_alt
  endif

  ! --- 2. a_n must carry the sign of F0
  if ( a_n * F0 .le. 0.d0 ) then
    floating_u_selftest = .false.
    if (my_id .eq. 0) write(*,'(A,2es22.14)') &
      ' FLOATING_U SELFTEST FAIL: sign(a_n) /= sign(F0):', a_n, F0
  endif

  ! --- 3. round trip. Te = 1 eV must reconstruct to exactly Lambda volts above the wall.
  ! --- In a single-T build the trace variable is T = 2*Te, and floating_u_norm has already
  ! --- halved Lambda, so feeding T = 2*Te_1eV is the consistent test.
  Te_1eV = EL_CHG * MU_ZERO * n_0
  if ( with_TiTe ) then
    u_1eV = C_T * Te_1eV
  else
    u_1eV = C_T * 2.d0 * Te_1eV
  endif
  volts = floating_u_volts(u_1eV)
  lam_eff = sheath_Lambda
  if ( abs(volts - lam_eff) .gt. tol * max(abs(lam_eff), 1.d0) ) then
    floating_u_selftest = .false.
    if (my_id .eq. 0) write(*,'(A,2es22.14)') &
      ' FLOATING_U SELFTEST FAIL: 1 eV must give Lambda volts:', volts, lam_eff
  endif

  if ( my_id .eq. 0 ) then
    write(*,'(A)')            ' --- floating_u normalisation ---'
    write(*,'(A,es22.14)')    '   a_n                     = ', a_n
    write(*,'(A,es22.14)')    '   C_T  (u per unit Te)    = ', C_T
    write(*,'(A,es22.14)')    '   C_V  (u per volt V_wall)= ', C_V
    write(*,'(A,es22.14)')    '   volts per unit u        = ', floating_u_volts(1.d0)
    write(*,'(A,f18.12,A)')   '   Te = 1 eV reconstructs to ', volts, ' V above the wall'
    if ( floating_u_selftest ) then
      write(*,'(A)')          '   selftest PASSED'
    else
      write(*,'(A)')          '   selftest FAILED - do not run'
    endif
  endif

end function floating_u_selftest


end module mod_floating_u
