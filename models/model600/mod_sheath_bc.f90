!> Model600 field-aligned sheath current law and normalization.
!!
!! Right-handed (R,Z,phi), with phi clockwise viewed from above:
!!   Phi_SI = F0*u/sqrt(mu0*rho_ref), Jtor = -zj/(mu0*R).
!! Hence J_n = -zj*(B.n)/(mu0*F0) for a field-aligned current.
!!
!! Ion and electron collection are evaluated separately. The electron reference
!! uses the prescribed Bohm/Chodura state and Lambda; changing the solution Vpar
!! changes ONLY the ion term. This is still a field-aligned collection model,
!! not a magnetic-presheath or cross-field wall-collection model.
!!
!! The repelling branch is Boltzmann; at Phi_plasma <= Phi_wall the electron
!! collection is saturated. There are no super-saturation tails or displaced
!! floating roots. At grazing incidence this evaluator does NOT supply the
!! missing tangent-wall closure of a zj trace constraint. See
!! doc/WEAK_SHEATH_FIXES_2026-09-05.md for scope and validation requirements.
module mod_sheath_bc

  implicit none

  private

  public :: sheath_frame_det, sheath_frame_frozen
  public :: sheath_norm, sheath_get_lambda, sheath_x_limited, sheath_current
  public :: sheath_V_wall_at, sheath_incidence
  public :: sheath_temp_floor, dsheath_temp_floor_dT
  public :: sheath_bohm_state

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

!> Determinant of a boundary node's frame: |sin(angle)| between the two first-derivative DOF
!! directions, normalised so the result is in [0,1] whatever the stored vector lengths.
!!
!! The grid builder never requires these two to be orthogonal or even independent. In the
!! ray-cast wall extension x(1,2,:) is the ray to the wall and x(1,3,:) an interpolated
!! along-wall angle (grid_xpoint_wall.f90:1252, :1311), and where the extension meets the wall
!! at grazing incidence the two become PARALLEL. The nodal derivative basis is then conditioned
!! as 1/det, and zj = Delta*psi - the very quantity the weak sheath row replaces - is built
!! from SECOND derivatives.
pure real*8 function sheath_frame_det(x2, x3)

  implicit none
  real*8, intent(in) :: x2(2), x3(2)
  real*8 :: n2, n3

  n2 = sqrt( x2(1)**2 + x2(2)**2 )
  n3 = sqrt( x3(1)**2 + x3(2)**2 )

  if ( n2 .gt. 1.d-30 .and. n3 .gt. 1.d-30 ) then
    sheath_frame_det = abs( x2(1)*x3(2) - x2(2)*x3(1) ) / (n2 * n3)
  else
    sheath_frame_det = 0.d0
  endif

end function sheath_frame_det


!> Is this node's frame too degenerate for the sheath to be imposed on it?
!!
!! A node that fails this test is treated as a FULLY frozen node - u Dirichlet as well as zj -
!! exactly as the code already treats a Dirichlet-u corner. Freezing zj alone would leave u
!! with no boundary condition at all on a weak type (dirichlet%u and natural%u are both
!! required .false. there), which is the "u free, zj pinned" configuration that produced the
!! boundary current filament. Freezing both is self-consistent: the node simply is not a
!! sheath node.
!!
!! The corresponding residual leak into the free neighbour's row is what sheath_weak_ufade
!! exists to remove, so that flag belongs on whenever this one is used.
pure logical function sheath_frame_frozen(x2, x3)

  use phys_module, only: sheath_weak_detmin

  implicit none
  real*8, intent(in) :: x2(2), x3(2)

  sheath_frame_frozen = .false.
  if ( sheath_weak_detmin .gt. 0.d0 ) &
    sheath_frame_frozen = sheath_frame_det(x2, x3) .lt. sheath_weak_detmin

end function sheath_frame_frozen


!> Normalisation constants of the characteristic in JOREK units.
!! @param a_n    +2*e*F0*sqrt(mu0*rho0)/m_i , so that e*Phi/(k*Te) = (a_n*u/2 - vw)/Te
!! @param c_sat  -a_n/2, the ion saturation current coefficient
!! @param vw     e*V_wall*mu0*n0, the wall bias in the same units
pure subroutine sheath_norm(a_n, c_sat, vw, V_wall_loc)

  use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT, EL_CHG
  use phys_module, only: F0, central_density, central_mass, sheath_V_wall

  implicit none
  real*8, intent(out) :: a_n, c_sat, vw
  ! Local metal potential in volts. A wall bias alone is not a gauge change
  ! when other plasma/wall voltage references remain fixed.
  real*8, intent(in), optional :: V_wall_loc
  ! --- NOTE: not "Vw" - Fortran is case insensitive, so that collides with the vw dummy above
  real*8 :: V_wall_use

  real*8 :: m_i, rho0

  V_wall_use = sheath_V_wall
  if ( present(V_wall_loc) ) V_wall_use = V_wall_loc

  m_i  = central_mass * ATOMIC_MASS_UNIT
  rho0 = central_density * 1.d20 * m_i

  ! --- The component velocity and curl conventions give Phi = +F0*u.
  a_n   = + 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * rho0) / m_i
  c_sat = - 0.5d0 * a_n
  vw    =   EL_CHG * V_wall_use * MU_ZERO * central_density * 1.d20

end subroutine sheath_norm


!> Positivity floor for the temperatures the SHEATH characteristic sees.
!!
!! Same functional form as corr_neg_temp1, but with its own floor T_min_sheath. The global
!! T_min_neg exists to keep every other term well behaved and is usually set well ABOVE divertor
!! temperatures: with T_min_neg = 3e-5 the knee sits at 1.93 eV, so a real target profile of
!! 0.3 -> 1.5 eV is delivered to the sheath as 1.14 -> 1.58 eV. Phi = Lambda*Te/e then inherits
!! that compression and the target's RADIAL potential gradient - which is what drives the E_r,
!! hence the PFR ExB drift, hence the in-out density asymmetry - is flattened by ~3.5x before the
!! boundary condition ever sees it.
!!
!! T_min_sheath < 0 (default) falls back to corr_neg_temp1, i.e. no change.
pure real*8 function sheath_temp_floor(T)

  use phys_module, only: T_min_neg, T_min_sheath, corr_neg_temp_coef

  implicit none
  real*8, intent(in) :: T
  real*8 :: L1, L2, Tf

  Tf = T_min_sheath
  if ( Tf .lt. 0.d0 ) Tf = T_min_neg
  L1 = Tf * corr_neg_temp_coef(1)
  L2 = Tf * corr_neg_temp_coef(2)

  if (L2 <= 0.d0) then
    sheath_temp_floor = max(T, L1, sheath_t_floor)
    return
  endif
  sheath_temp_floor = T
  if ( T .lt. L1 + L2 ) sheath_temp_floor = L1 + L2 * exp( (T - (L1+L2)) / L2 )

end function sheath_temp_floor


!> d(sheath_temp_floor)/dT, for the Jacobian chain rule.
pure real*8 function dsheath_temp_floor_dT(T)

  use phys_module, only: T_min_neg, T_min_sheath, corr_neg_temp_coef

  implicit none
  real*8, intent(in) :: T
  real*8 :: L1, L2, Tf

  Tf = T_min_sheath
  if ( Tf .lt. 0.d0 ) Tf = T_min_neg
  L1 = Tf * corr_neg_temp_coef(1)
  L2 = Tf * corr_neg_temp_coef(2)

  if (L2 <= 0.d0) then
    dsheath_temp_floor_dT = 0.d0
    if (T > max(L1,sheath_t_floor)) dsheath_temp_floor_dT = 1.d0
    return
  endif
  dsheath_temp_floor_dT = 1.d0
  if ( T .lt. L1 + L2 ) dsheath_temp_floor_dT = exp( (T - (L1+L2)) / L2 )

end function dsheath_temp_floor_dT


!> Wall potential at major radius R, in volts.
!>
!>   V_wall(R) = sheath_V_wall + sheath_V_wall_asym * tanh( (R - R0) / dR )
!>
!> An ANTISYMMETRIC bias about R0: the inner target sees about
!> sheath_V_wall - sheath_V_wall_asym and the outer about + it, so the
!> target-to-target difference is roughly 2*sheath_V_wall_asym while the mean is
!> unchanged. tanh rather than a step, because a discontinuity in the imposed
!> trace is turned straight into boundary vorticity by w = grad^2 u.
!>
!> This is a DIAGNOSTIC, not physics. It answers: if a potential difference of
!> this size existed between the targets, would the PFR ExB drift carry enough
!> to produce the in-out asymmetry? A positive answer means the transport chain
!> works and the self-consistent Phi is merely too small; a negative answer means
!> the chain is broken downstream of the boundary condition.
pure real*8 function sheath_V_wall_at(BigR)

  use phys_module, only: sheath_V_wall, sheath_V_wall_asym, &
                         sheath_V_wall_R0, sheath_V_wall_dR

  implicit none
  real*8, intent(in) :: BigR

  sheath_V_wall_at = sheath_V_wall
  if ( (sheath_V_wall_asym .ne. 0.d0) .and. (sheath_V_wall_dR .gt. 0.d0) )      &
    sheath_V_wall_at = sheath_V_wall + sheath_V_wall_asym                       &
                     * tanh( (BigR - sheath_V_wall_R0) / sheath_V_wall_dR )

end function sheath_V_wall_at


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


!> Signed Chodura factor and its derivative with respect to b_n=B.n/|B|.
!! Shared by sheath quadrature and the weak-sheath nodal flow condition.
!! Subtract at least the analytic zero-incidence value, avoiding a sign jump
!! caused by rounded c3. No cosh overflow and no division by b_n at tangency.
pure subroutine sheath_incidence(bn, g, dg)
  use phys_module, only: vpar_smoothing, vpar_smoothing_coef
  real*8, intent(in) :: bn
  real*8, intent(out) :: g, dg
  real*8 :: c1, c2, c3, t, baseline, magnitude
  c1 = vpar_smoothing_coef(1)
  c2 = vpar_smoothing_coef(2)
  c3 = vpar_smoothing_coef(3)
  if (.not. vpar_smoothing) then
    g = sign(1.d0, bn)
    if (bn == 0.d0) g = 0.d0
    dg = 0.d0
  elseif (c2 > 0.d0) then
    t = tanh((abs(bn)-c1)/c2)
    baseline = max(c3, 0.25d0*(1.d0+tanh(-c1/c2))**2)
    magnitude = max(0.d0, 0.25d0*(1.d0+t)**2-baseline)
    g = sign(magnitude, bn)
    dg = 0.d0
    if (magnitude > 0.d0 .or. (bn == 0.d0 .and. baseline <= &
        0.25d0*(1.d0+tanh(-c1/c2))**2)) &
      dg = 0.5d0*(1.d0+t)*(1.d0-t*t)/c2
  else
    ! c1 > 0 is validated at namelist initialization.
    t = tanh(bn/c1)
    g = t
    dg = (1.d0-t*t)/c1
  endif
end subroutine sheath_incidence

!> Field-aligned entrance state and thermal Jacobian, with frozen geometry A=g/B.
!! Ab is its tangential derivative supplied by the geometric discretization.
!! This condition does not cancel ExB flux by dividing it by B.n.
pure subroutine sheath_bohm_state(T, Tb, A, Ab, velocity, velocity_b, dT, dbT, dbTb)
  use phys_module, only: GAMMA, T_min_neg, T_min_sheath, corr_neg_temp_coef
  real*8, intent(in) :: T(2), Tb(2), A, Ab
  real*8, intent(out) :: velocity, velocity_b, dT(2), dbT(2), dbTb(2)
  real*8 :: Tf(2), fp(2), fpp(2), total_b, cs, csT, csTT, floor_T, L1, L2
  integer :: k
  floor_T = T_min_sheath
  if (floor_T < 0.d0) floor_T = T_min_neg
  L1 = floor_T*corr_neg_temp_coef(1)
  L2 = floor_T*corr_neg_temp_coef(2)
  do k=1,2
    Tf(k) = sheath_temp_floor(T(k))
    fp(k) = dsheath_temp_floor_dT(T(k))
    fpp(k) = 0.d0
    if (L2 > 0.d0 .and. T(k) < L1+L2) fpp(k) = fp(k)/L2
  enddo
  cs = sqrt(GAMMA*sum(Tf))
  csT = GAMMA/(2.d0*cs)
  csTT = -GAMMA**2/(4.d0*cs**3)
  total_b = sum(fp*Tb)
  velocity = A*cs
  velocity_b = Ab*cs+A*csT*total_b
  dT = A*csT*fp
  dbT = Ab*csT*fp+A*(csTT*fp*total_b+csT*fpp*Tb)
  dbTb = A*csT*fp
end subroutine sheath_bohm_state

!> Field-aligned current-variable target, with partial derivatives.
!! g_bn is signed. Vpar is the coefficient v_parallel/|B|, not a speed.
!! A present zero Vpar is a real stagnation state, not an absent argument.
!! Lambda defines the electron reference relative to the prescribed ion state.
!! x_out is the RAW chi-Lambda, never a soft-clipped exponent.
subroutine sheath_current(u, rho, Ti, Te, g_bn, sgn_bn, Btot,        &
                          zj_sh, dzj_du, dzj_drho, dzj_dTi, dzj_dTe, &
                          zj_sat, x_out, V_wall_loc, vpar, dzj_dvpar, &
                          dzj_dg, dzj_dB)
  use phys_module, only: GAMMA, sheath_jsat_from_vpar, sheath_jsat_vpar_min
  implicit none
  real*8, intent(in) :: u, rho, Ti, Te, g_bn, sgn_bn, Btot
  real*8, intent(out) :: zj_sh, dzj_du, dzj_drho, dzj_dTi, dzj_dTe, zj_sat, x_out
  real*8, intent(in), optional :: V_wall_loc, vpar
  real*8, intent(out), optional :: dzj_dvpar, dzj_dg, dzj_dB
  real*8 :: an, c, vw, rl, til, tel, total, cs, lam, lti, lte, chi
  real*8 :: zref, zi, ze, vi, vf, dvi, ziT, zig, ziB, eref, dchi, f
  logical :: actual_ion

  call sheath_norm(an, c, vw, V_wall_loc)
  rl = max(rho, sheath_rho_floor)
  til = max(Ti, sheath_t_floor)
  tel = max(Te, sheath_t_floor)
  total = til + tel
  cs = sqrt(GAMMA*total)
  call sheath_get_lambda(til, tel, lam, lti, lte)
  chi = (0.5d0*an*u-vw)/tel
  x_out = chi-lam

  ! Btot>0 follows from nonzero F0 and R>0, checked at setup/geometry.
  zref = c*rl*g_bn*cs/Btot
  zi = zref
  ziT = zref/(2.d0*total)
  zig = c*rl*cs/Btot
  ziB = -zref/Btot
  dvi = 0.d0
  actual_ion = sheath_jsat_from_vpar .and. present(vpar)
  if (actual_ion) then
    ! Outgoing ions only; no hidden full-Bohm reset at vpar=0.
    ! A nonzero user floor is an explicit legacy ion-collection approximation.
    vf = sheath_jsat_vpar_min*abs(g_bn)*cs/Btot
    vi = max(sgn_bn*vpar, vf)
    zi = c*rl*sgn_bn*vi
    if (sgn_bn*vpar > vf) then
      dvi = 1.d0
      ziT = 0.d0
      zig = 0.d0
      ziB = 0.d0
    else
      ziT = zi/(2.d0*total)
      zig = c*rl*sheath_jsat_vpar_min*cs/Btot
      ziB = -zi/Btot
    endif
  endif
  zj_sat = zi

  ! exp(lam-max(chi,0)) is bounded on the electron-attracting side.
  ! A generalized derivative from the repelling side is used at chi=0.
  eref = exp(lam-max(chi,0.d0))
  ze = zref*eref
  dchi = 0.d0
  if (chi >= 0.d0) dchi = 1.d0
  zj_sh = zi-ze
  if (.not. actual_ion .and. abs(x_out) < 1.d-4 .and. chi >= 0.d0) then
    ! Cancellation-free floating limit: 1-exp(-X).
    f = x_out*(1.d0-x_out*(0.5d0-x_out*(1.d0/6.d0-x_out/24.d0)))
    zj_sh = zref*f
  endif
  dzj_du = ze*dchi*0.5d0*an/tel
  dzj_drho = zj_sh/rl
  dzj_dTi = ziT-ze*(0.5d0/total+lti)
  dzj_dTe = ziT-ze*(0.5d0/total+lte+dchi*chi/tel)
  if (present(dzj_dvpar)) dzj_dvpar = c*rl*dvi
  if (present(dzj_dg)) dzj_dg = zig-c*rl*cs/Btot*eref
  if (present(dzj_dB)) dzj_dB = ziB+ze/Btot

  if (rho < sheath_rho_floor) dzj_drho = 0.d0
  if (Ti < sheath_t_floor) dzj_dTi = 0.d0
  if (Te < sheath_t_floor) dzj_dTe = 0.d0
end subroutine sheath_current

end module mod_sheath_bc
