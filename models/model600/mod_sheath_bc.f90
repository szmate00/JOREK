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
!!   bcs(1)%natural%u    = .true.    ! the sheath boundary condition
!!   bcs(1)%dirichlet%u  = .false.   ! u is free and set by charge continuity at the wall
!!   bcs(1)%natural%zj   = .true.    ! REQUIRED - see below; without it the j-V loop is open
!!   bcs(1)%dirichlet%zj = .false.   ! ... and its Dirichlet has to come off
!!   bcs(1)%dirichlet%w  = .true.    ! KEEP, and leave natural%w = .false.
!!   bcs(1)%mach1        = .true.    ! j_sat assumes the Mach 1 condition at the same wall
!!   bc_natural_open     = .true.    ! the boundary integrals live in that branch of construct_matrix
!!
!! Leave at least one boundary type with dirichlet%u = .true. (typically the main chamber wall):
!! u enters the vorticity equation only through its gradient, and the sheath term loses its grip
!! on u in ion saturation, so its constant mode would otherwise be undetermined. Optional knobs:
!! sheath_Lambda (Lambda_0, <=0 computes it from central_mass), sheath_Lambda_local,
!! sheath_V_wall, sheath_X_min, sheath_smooth_dX, sheath_min_bn, sheath_ramp_time.
!!
!! WHY natural%zj IS REQUIRED (corrected 2026-08-19). This header used to say that dirichlet%zj
!! could stay on because "the frozen zj trace cancels exactly between the strong-form volume term
!! and the added surface flux". That is FALSE, and the diagnostic falsified it directly: with
!! dirichlet%zj = .true., I_Ampere = oint zj (B.n) stayed constant to four digits for a whole run
!! while I_wall ran from -15540 A to -23 A. The sheath adapted all the way to floating, the plasma
!! current could not move at all, and u diverged trying to reconcile them. The j-V loop is OPEN
!! whenever zj at the wall is pinned: the boundary condition can then only move the potential, and
!! if the delivered current exceeds j_sat no potential satisfies the characteristic.
!!
!! What that measurement shows is that dirichlet%zj must be OFF. It does not by itself show that
!! natural%zj must be on, and the two routes differ here:
!!
!!  - NATURAL%U route: nothing replaces the boundary zj rows, so releasing the Dirichlet leaves
!!    them with an incomplete weak form. natural%zj = .true. is then genuinely required; it
!!    couples zj at the boundary to grad(psi).n - the NORMAL derivative of psi, which
!!    dirichlet%psi does NOT pin (that fixes the value and the tangential derivative only), and
!!    that is the degree of freedom through which the wall current responds to the sheath.
!!
!!  - WEAK (sheath_zj_weak) route: the boundary zj trace rows are REPLACED by the characteristic
!!    itself, so the surface term would only contribute to rows that are then overwritten. The
!!    j-V loop closes through the constraint: zj at the wall is not pinned, it is set to
!!    zj_sh(u,rho,Ti,Te), and the interior zj = Delta*psi equation makes psi's normal derivative
!!    follow. MEASURED: dirichlet%zj = .false. with natural%zj left OFF ran 3900 steps with
!!    I_sheath = I_Ampere to four significant figures on both targets. natural%zj is therefore
!!    NOT required on this route.
!!
!!    The one thing it still touches there: the trace space spans only the value and tangential
!!    DOFs, so the zj rows at the NORMAL-derivative DOFs of boundary nodes are not replaced and
!!    keep an incomplete weak form without it. That is second order and was not measurable at
!!    0.00 % closure, but it is the reason to turn natural%zj on if a boundary zj row ever falls
!!    below the trace accumulator's degeneracy floor.
!!
!! Its Jacobian was wrong until 2026-08-19: the trial loop produced columns only at
!! l2 = direction(l), and psi_t synthesised a fictitious normal derivative on the value DOF, so the
!! true entry was missing and a spurious one took its place - an effectively explicit O(1/h)
!! boundary term that grew boundary structures over ~20 steps. mod_boundary_matrix_open now splits
!!     grad(f).n = gpn_s * f_s + gpn_t * f_t
!! and assembles the second half at direction_perp(l), the DOF that actually carries the normal
!! derivative (l=1 value pairs with the normal first derivative, l=2 tangential with the mixed
!! second derivative - the same pairing the field evaluation uses when it builds eq_t). The split
!! is algebraically exact; only the column assignment and an orientation sign changed.
!!
!! dirichlet%w still stays on: the w surface term is not needed by the sheath, and with the
!! Dirichlet present it reaches only rows that Dirichlet overwrites, so it imposes nothing.
module mod_sheath_bc

  implicit none

  private

  public :: sheath_norm, sheath_get_lambda, sheath_x_limited, sheath_current
  public :: sheath_V_wall_at
  public :: sheath_temp_floor, dsheath_temp_floor_dT

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
pure subroutine sheath_norm(a_n, c_sat, vw, V_wall_loc)

  use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT, EL_CHG
  use phys_module, only: F0, central_density, central_mass, sheath_V_wall

  implicit none
  real*8, intent(out) :: a_n, c_sat, vw
  ! --- Local wall potential in VOLTS, overriding the uniform sheath_V_wall. A uniform bias is a
  ! --- gauge: Phi shifts everywhere, E = -grad(Phi) does not, and no drift responds to it. Only a
  ! --- bias that DIFFERS between the two targets makes a potential difference across the private
  ! --- flux region, which is what drives a PFR ExB flow.
  real*8, intent(in), optional :: V_wall_loc
  ! --- NOTE: not "Vw" - Fortran is case insensitive, so that collides with the vw dummy above
  real*8 :: V_wall_use

  real*8 :: m_i, rho0

  V_wall_use = sheath_V_wall
  if ( present(V_wall_loc) ) V_wall_use = V_wall_loc

  m_i  = central_mass * ATOMIC_MASS_UNIT
  rho0 = central_density * 1.d20 * m_i

  ! --- NOTE the minus sign: Phi = -F0*u in the code's variables (see the module header)
  a_n   = - 2.d0 * EL_CHG * F0 * sqrt(MU_ZERO * rho0) / m_i
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
!! @param dzj_dg derivative of zj_sh with respect to g_bn, holding the state and Btot fixed
!! @param dzj_dB derivative of zj_sh with respect to Btot, holding the state and g_bn fixed
subroutine sheath_current(u, rho, Ti, Te, g_bn, sgn_bn, Btot,        &
                          zj_sh, dzj_du, dzj_drho, dzj_dTi, dzj_dTe, &
                          zj_sat, x_out, V_wall_loc, vpar, dzj_dvpar, &
                          dzj_dg, dzj_dB)

  use phys_module, only: GAMMA, sheath_sat_slope, sheath_sat_slope_e, &
                         sheath_v_perp, sheath_dfdx_min, &
                         sheath_jsat_from_vpar, sheath_jsat_vpar_min

  implicit none

  real*8, intent(in)  :: u, rho, Ti, Te, g_bn, sgn_bn, Btot
  real*8, intent(out) :: zj_sh, dzj_du, dzj_drho, dzj_dTi, dzj_dTe
  real*8, intent(out) :: zj_sat, x_out
  real*8, intent(in), optional :: V_wall_loc
  !> Parallel velocity at the wall in Vpar units (v_par/|B|), i.e. the model's own var_vpar.
  !> Only consulted when sheath_jsat_from_vpar is set; absent restores the Mach 1 behaviour.
  real*8, intent(in),  optional :: vpar
  real*8, intent(out), optional :: dzj_dvpar
  real*8, intent(out), optional :: dzj_dg, dzj_dB

  real*8 :: a_n, c_sat, vw
  real*8 :: rho_l, Ti_l, Te_l, T_l, cs, g_eff
  real*8  :: v_mach, v_eff, v_flr, dveff_dvpar, dlnv_dT, v_par_sat
  real*8  :: dzsat_dg, dzsat_dB
  logical :: use_vpar
  real*8 :: lam, dlam_dTi, dlam_dTe
  real*8  :: x, x_lim, dxlim_dx, expmx, f, fp, sp, dsp
  real*8  :: dx_du, dx_dTi, dx_dTe
  logical :: x_frozen

  if ( present(V_wall_loc) ) then
    call sheath_norm(a_n, c_sat, vw, V_wall_loc)
  else
    call sheath_norm(a_n, c_sat, vw)
  endif

  rho_l = max(rho, sheath_rho_floor)
  Ti_l  = max(Ti,  sheath_t_floor)
  Te_l  = max(Te,  sheath_t_floor)
  T_l   = Ti_l + Te_l

  cs = sqrt( GAMMA * T_l )

  ! --- No floor on g(b_n) here. Flooring the magnitude while keeping the sign of B.n is
  ! --- discontinuous through tangency, and the real problem it was aimed at is SOLVABILITY, not
  ! --- magnitude: as g -> 0 the saturation current zj_sat -> 0, so the characteristic is asked
  ! --- for f = zj0/zj_sat -> infinity, which has no solution (f <= 1 on the ion side, and
  ! --- >= -(exp(-X_min)-1) on the electron side) and drives u monotonically. That is handled in
  ! --- mod_boundary_matrix_open by gating the whole term - residual and Jacobian together - with
  ! --- the smooth obliqueness weight b_n^2/(b_n^2 + sheath_min_bn^2).
  g_eff = g_bn

  ! --- Parallel velocity that carries the ion saturation current.
  ! ---
  ! --- v_mach = g(b_n)*cs/|B| is EXACTLY what bcs(...)%mach1 imposes on Vpar at this same wall
  ! --- (Mach1BC = -Vpar0 + direction/Btot*factor*cs0 + drift, and direction*factor is g(b_n)),
  ! --- so the default below is consistent with the Mach 1 condition by construction - which is
  ! --- what the module header asserts.
  ! ---
  ! --- It is also an ASSUMPTION: it fixes the ion flux at the Bohm value. A detaching target is
  ! --- sub-sonic, and there the Mach 1 value overestimates the flux and with it j_sat, which is
  ! --- the wrong direction for a leg one is trying to detach. sheath_jsat_from_vpar swaps in the
  ! --- Vpar the solution actually has - same units, so it is a one-symbol change - and it
  ! --- inherits whatever drift correction Mach1BC carries into Vpar.
  ! --- Whether to consult the solution's Vpar at all. vpar == 0 exactly means the caller has
  ! --- none to offer (a model built without var_vpar hands over a hard zero), NOT a stagnation
  ! --- point, which is measure-zero in a continuous field.
  dveff_dvpar = 0.d0
  dzsat_dg    = 0.d0
  dzsat_dB    = 0.d0
  use_vpar    = sheath_jsat_from_vpar .and. present(vpar)
  if ( use_vpar ) then
    if ( vpar .eq. 0.d0 ) use_vpar = .false.
  endif

  if ( use_vpar ) then
    ! --- v_mach = g(b_n)*cs/|B| is EXACTLY what bcs(...)%mach1 imposes on Vpar at this same wall
    ! --- (Mach1BC = -Vpar0 + direction/Btot*factor*cs0 + drift, and direction*factor is g(b_n)),
    ! --- so the default below is consistent with the Mach 1 condition by construction. It is also
    ! --- an ASSUMPTION: it fixes the ion flux at the Bohm value, and a detaching target is
    ! --- SUB-SONIC, where the Mach 1 value overestimates the flux and with it j_sat - the wrong
    ! --- direction for a leg one is trying to detach. Swapping in the Vpar the solution actually
    ! --- has is a one-symbol change (same units) and inherits Mach1BC's drift correction.
    v_mach = g_eff * cs / Btot
    v_flr  = sheath_jsat_vpar_min * abs(v_mach)
    ! --- Take Vpar only where it flows INTO the wall and clears the floor. The floor is not
    ! --- optional: zj_sat -> 0 asks the characteristic for f = zj0/zj_sat -> inf, which has no
    ! --- solution. Below it, and for flow directed away from the wall, fall back to the floor -
    ! --- a hard clip, so like the rho and T floors it drops its own sensitivity.
    if ( abs(v_mach) .gt. 1.d-300 .and. vpar * v_mach .gt. 0.d0 &
         .and. abs(vpar) .gt. v_flr ) then
      v_eff       = vpar
      dveff_dvpar = 1.d0
      dlnv_dT     = 0.d0        ! v_eff no longer carries the cs temperature dependence
      ! --- The solution variable Vpar is independent of the magnetic geometry in this branch.
      ! --- In particular, do not reuse the Mach-1 g/B derivative for it.
      dzsat_dg    = 0.d0
      dzsat_dB    = 0.d0
    else
      v_eff       = sign(v_flr, v_mach)
      dlnv_dT     = 0.5d0 / T_l
      ! --- v_eff = sheath_jsat_vpar_min*g_bn*cs/Btot on the fallback branch.
      dzsat_dg    = c_sat * rho_l * sheath_jsat_vpar_min * cs / Btot
    endif
    zj_sat = c_sat * rho_l * v_eff
    if ( dveff_dvpar .eq. 0.d0 ) dzsat_dB = - zj_sat / Btot
  else
    ! --- DEFAULT PATH, written out literally rather than as a special case of the above, so that
    ! --- sheath_jsat_from_vpar = .false. is BIT-IDENTICAL to the pre-flag code and not merely
    ! --- algebraically equal. A 1e-16 re-association is enough to change which side of a marginal
    ! --- configuration a long run falls on, which makes an A/B uninterpretable.
    if ( sheath_v_perp .gt. 0.d0 ) then
      ! --- CROSS-FIELD FLUX FLOOR. j_sat = e*n*(c_s*g(b_n) + v_perp). The parallel flux vanishes
      ! --- where the field grazes the wall, but perpendicular transport still delivers particles,
      ! --- so the saturation current does not go to zero at tangency and the characteristic stays
      ! --- solvable there. sgn_bn rather than sign(g_eff) so the direction is still defined at
      ! --- exact tangency, where g_eff = 0. Only the magnitude gets the floor.
      v_par_sat = abs(g_eff) * cs
      zj_sat    = c_sat * rho_l * sgn_bn * ( v_par_sat + sheath_v_perp ) / Btot
      ! --- d(ln zj_sat)/dT: only the parallel part carries the c_s ~ sqrt(T) dependence, so the
      ! --- 1/(2T) is diluted by the floor and -> 0 as the floor dominates.
      dlnv_dT   = 0.5d0 / T_l * v_par_sat / max( v_par_sat + sheath_v_perp, 1.d-300 )
    else
      zj_sat  = c_sat * rho_l * g_eff * cs / Btot
      dlnv_dT = 0.d0            ! unused on this path; the derivative below is written out too
    endif
    ! --- On either Mach-1 branch, d[sgn(g)*(|g|*cs+v_perp)]/dg = cs away
    ! --- from the sign-change point. This direct derivative remains finite as g -> 0; forming it
    ! --- as zj_sat/g would spuriously create a v_perp/g singularity when the cross-field floor is
    ! --- active. The caller sets dg/dpsi=0 where the clipped Chodura factor is exactly zero.
    dzsat_dg = c_sat * rho_l * cs / Btot
    dzsat_dB = - zj_sat / Btot
  endif

  call sheath_get_lambda(Ti_l, Te_l, lam, dlam_dTi, dlam_dTe)

  ! --- X = e*Phi/(k*Te) - Lambda
  x = ( 0.5d0 * a_n * u - vw ) / Te_l - lam

  call sheath_x_limited(x, x_lim, dxlim_dx)
  x_out = x_lim

  x_frozen = .false.
  if ( x_lim .lt. -sheath_exp_max ) then
    expmx    = exp(sheath_exp_max)  ! electron branch, capped. Only reachable when the user has
                                    ! effectively disabled the limiter with a very negative X_min
    x_frozen = .true.               ! the residual no longer depends on X here, so neither may the
                                    ! Jacobian: fp must be 0, not expmx (which would be 1e13 off)
  elseif ( x_lim .gt. sheath_exp_max ) then
    expmx = 0.d0                    ! ion saturation; the true derivative is exp(-X) ~ 0 too
  else
    expmx = exp(-x_lim)
  endif

  ! --- 1 - exp(-X) cancels catastrophically near X = 0 - which is exactly the floating
  ! --- condition a grounded wall sits at, so it is the regime of interest rather than a corner
  ! --- case. At X = 1e-4 the naive form keeps only ~12 significant digits. Fortran has no expm1
  ! --- intrinsic, so use the Taylor series of 1 - exp(-X) below the crossover; four terms give
  ! --- a relative error of ~1e-18 there, well inside double precision.
  if ( abs(x_lim) .lt. 1.d-4 ) then
    f = x_lim * ( 1.d0 - x_lim * ( 0.5d0 - x_lim * ( 1.d0/6.d0 - x_lim/24.d0 ) ) )
  else
    f = 1.d0 - expmx                ! the characteristic
  endif
  fp = expmx * dxlim_dx             ! d f / d X (no cancellation: expmx -> 1 as X -> 0)
  if ( x_frozen ) fp = 0.d0         ! stay consistent with the frozen residual above

  ! --- Finite conductance at ion saturation. The forward characteristic saturates EXACTLY at
  ! --- j_sat, so a plasma delivering even 3% more current than the sheath can pass - which a
  ! --- restart equilibrium built without this boundary condition has no reason to respect - has
  ! --- NO solution: as X -> infinity the residual settles at a small constant and drives u
  ! --- linearly for ever. Adding s*ln(1+exp(X)) makes f unbounded above, so every demanded
  ! --- current is reachable at a finite X, while leaving the electron branch and the floating
  ! --- potential essentially untouched (the softplus is exponentially small for X << 0).
  ! --- C-infinity, and applied to f and its derivative together.
  if ( sheath_sat_slope .ne. 0.d0 ) then
    if ( x_lim .gt. sheath_exp_max ) then
      sp = x_lim; dsp = 1.d0
    elseif ( x_lim .lt. -sheath_exp_max ) then
      sp = 0.d0;  dsp = 0.d0
    else
      sp = log( 1.d0 + exp(x_lim) ); dsp = 1.d0 / ( 1.d0 + exp(-x_lim) )
    endif
    f = f + sheath_sat_slope * sp
    if ( .not. x_frozen ) fp = fp + sheath_sat_slope * dsp * dxlim_dx
  endif

  ! --- Finite conductance on the ELECTRON side, the exact mirror of the block above and the
  ! --- prerequisite for TWO sheath surfaces. As X -> -inf the limiter drives dxlim_dx -> 0, so
  ! --- fp -> 0 and the replaced zj row loses its u column: the sheath keeps no authority over u
  ! --- there, and with two floating patches nothing damps the exchange between them.
  ! ---
  ! --- (x - x_lim) is exactly the amount the limiter clipped: identically 0 on the ion branch,
  ! --- and -> x - X_min (linear, NEGATIVE) deep on the electron branch. Its derivative is
  ! --- (1 - dxlim_dx). So f acquires a LINEAR tail and fp -> s_e /= 0, while the ion branch and
  ! --- the fixed point are untouched to O(s_e).
  ! ---
  ! --- NOT to be done by blending x_lim itself: f and fp are EXPONENTIAL in x_lim (see expmx
  ! --- above), so blending inside the exponent multiplies the demanded electron current by
  ! --- exp(s_e*|X|) - at X = -19.6, s_e = 0.1 that is -104*j_sat instead of -19.1*j_sat. That
  ! --- was tried, and it made the two-surface blow-up worse rather than better.
  if ( sheath_sat_slope_e .gt. 0.d0 ) then
    f  = f  + sheath_sat_slope_e * ( x - x_lim )
    fp = fp + sheath_sat_slope_e * ( 1.d0 - dxlim_dx )
  endif

  ! --- Quasi-Newton slope floor. sheath_weak_rmax bounds the RESIDUAL, while a small df/dX weakens
  ! --- the potential column on both plateaus. Flooring fp regularises that column. This is not a
  ! --- Levenberg-Marquardt update and, because the replacement row also contains a free zj column
  ! --- and other coupled variables, it is not by itself a strict bound on the global delta-X.
  ! ---
  ! --- Applied to fp ONLY, never to f, so the fixed point is exactly where it was: at convergence
  ! --- the residual vanishes and the Jacobian's magnitude no longer matters. That is the advantage
  ! --- over sheath_sat_slope_e, which alters f itself and therefore moves the solution. df/dX = 1
  ! --- at the floating potential, so any floor below 1 is inactive in the regime of interest.
  ! --- fp is non-negative by construction (exp(-x_lim)*dxlim_dx plus non-negative terms), so a
  ! --- plain max is the right form.
  if ( sheath_dfdx_min .gt. 0.d0 ) fp = max( fp, sheath_dfdx_min )

  ! --- derivatives of X (Lambda contributes through Ti and Te when sheath_Lambda_local)
  dx_du  =   0.5d0 * a_n / Te_l
  dx_dTi = - dlam_dTi
  dx_dTe = - ( x + lam ) / Te_l - dlam_dTe

  zj_sh    = zj_sat * f
  dzj_du   = zj_sat * fp * dx_du
  dzj_drho = zj_sat * f / rho_l
  if ( use_vpar .or. sheath_v_perp .gt. 0.d0 ) then
    dzj_dTi = zj_sat * f * dlnv_dT + zj_sat * fp * dx_dTi
    dzj_dTe = zj_sat * f * dlnv_dT + zj_sat * fp * dx_dTe
  else
    ! --- original form, byte for byte
    dzj_dTi = zj_sat * f / (2.d0 * T_l) + zj_sat * fp * dx_dTi
    dzj_dTe = zj_sat * f / (2.d0 * T_l) + zj_sat * fp * dx_dTe
  endif
  ! --- new column: zj_sat is linear in v_eff, so this is exact where the floor is not active
  if ( present(dzj_dvpar) ) dzj_dvpar = c_sat * rho_l * dveff_dvpar * f
  if ( present(dzj_dg) )    dzj_dg    = dzsat_dg * f
  if ( present(dzj_dB) )    dzj_dB    = dzsat_dB * f

  ! --- the floors are hard clips, so drop the corresponding sensitivity to stay consistent
  if ( rho .lt. sheath_rho_floor ) dzj_drho = 0.d0
  if ( Ti  .lt. sheath_t_floor   ) dzj_dTi  = 0.d0
  if ( Te  .lt. sheath_t_floor   ) dzj_dTe  = 0.d0

end subroutine sheath_current

end module mod_sheath_bc
