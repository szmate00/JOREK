!> Sheath current boundary condition for the electrostatic potential (model600).
!!
!! Imposes the full sheath current-voltage characteristic (Stangeby 2.68) at a material wall,
!!
!!     j = j_sat * ( 1 - exp( Lambda - e*Phi/(k_B*Te) ) )  ==  j_sat * f ,
!!
!! as ONE scalar constraint per boundary node, replacing the u (potential) node rows. Here j is the
!! JOREK toroidal current variable zj, j_sat is the ion saturation current expressed in the same
!! toroidal variable, and Lambda = sheath_Lambda is the floating sheath potential drop in units of
!! Te/e (2.84 for deuterium with Ti=Te; 3 by convention).
!!
!! The derivation and the linearisation are Artola, "Sheath boundary conditions for the electric
!! potential in JOREK" (2026-07-30), eqs. (1)-(18). This module implements his eq. (17) verbatim,
!! generalised to a separate Ti/Te build. Everything below is stated in JOREK units.
!!
!! WHY THIS FORM AND NOT j - j_sat*f = 0.
!! The residual written directly in the current is unusable: at u = 0 the exponent is +Lambda, so
!! f0 = 1 - e^3 = -19.1 and the row demands 19 ion saturation currents out of a node that has not
!! yet moved. Artola's scaling by
!!
!!     xi = (2*Te/a_n) * exp( e*Phi/(k_B*Te) - Lambda )
!!
!! removes exactly that: xi*(j0/j_sat0 - f0) is a POTENTIAL, bounded and of order Te/e, and the
!! coefficient of delta_u is identically 1. Two consequences matter.
!!   * Unit diagonal in u. The u trace DOF and the Te/rho/zj trace DOFs live in the same nodal
!!     frame, so no element Jacobian, no B.n and no logical-to-(R,Z) conversion ever enters this
!!     row. Its conditioning is therefore independent of the boundary node type: the row is the
!!     same object on a flux-aligned type 1 node and on a ray-cast type 4/9 leg-end node, and it
!!     does not weaken at grazing incidence the way a residual proportional to j_sat does.
!!   * Graceful branch limits. On the electron branch (Phi far below floating) xi -> 0 and the row
!!     degenerates to delta_u = 0. Deep in ion saturation xi -> infinity and, after the exact row
!!     normalisation below, the row degenerates to a condition on delta_j, which is the correct
!!     physics: a saturated sheath constrains the current, not the potential.
!!
!! WHAT IS LAGGED, EXACTLY. Artola linearises the CURRENT residual
!!
!!     G = j - j_sat*f
!!
!! and only then multiplies the resulting linear equation by the scalar -xi/j_sat0 (his eqs. 9-14).
!! So the row is  -(xi/j_sat0) * [ dG/dx . delta_x + G0 ]  with xi/j_sat0 FROZEN. This is not a
!! shortcut, it is the whole construction: differentiating xi as well would add a term
!! res*a_n/(2*Te) to the u column and put the exponential stiffness straight back into the matrix,
!! whereas a frozen positive multiplier changes only the step length and cannot move the fixed
!! point, since G0 -> 0 there. The self-test therefore checks each coefficient against
!! -(xi/j_sat0) * dG/dx with G coded independently from eq. (6), NOT against a derivative of the
!! scaled residual - those are different objects and the difference is the whole point.
!!
!! Two further quantities are lagged in the ordinary way, both evaluated at the current state: the
!! magnetic geometry (|B| and sign(B.n)), as everywhere else in the boundary path, and, inside the
!! Jacobian coefficients only, Artola's substitution j0/j_sat0 -> f0 (equal when the BC holds). The
!! latter keeps those coefficients bounded by 1 instead of proportional to a possibly large current
!! mismatch, and differs from the exact derivative at second order in the step.
!!
!! WHAT IS NOT HERE. No ratio gate, no minimum-|B.n|, no saturation-slope limiter, no current clip,
!! no relaxation gain, no frame-determinant threshold and no incidence floor. The row needs none:
!! see the two consequences above. The only non-physics numbers in this file are an IEEE range guard
!! on the exponent and a division guard on the density, both marked as such.
module mod_sheath_current

  implicit none
  private

  public :: sheath_current_norm, sheath_current_row, sheath_current_selftest
  public :: sheath_current_volts
  public :: sheath_diag_reset, sheath_diag_add, sheath_diag_print

  !> IEEE range guard on the exponent e*Phi/(k_B*Te) - Lambda. exp() of +-354 is representable and
  !! its square is not; a boundary node 354 e-foldings away from the characteristic carries no
  !! information beyond "saturated", which is what the clamped value already expresses. This is a
  !! representability bound, not a physics threshold, and the solution never approaches it.
  real*8, parameter :: exp_arg_max = 354.d0

  !> Division guard for the density appearing in j_sat. Only stops a hard zero-divide; a boundary
  !! node with rho <= 0 is a failed solution that the positivity scheme owns, not this row.
  real*8, parameter :: rho_floor   = 1.d-10

  ! --- Monitoring state, accumulated over one matrix construction and printed once per timestep.
  ! --- Read-only as far as the physics is concerned: nothing here feeds back into the row.
  integer :: nd_rows  = 0        !< number of sheath rows assembled
  integer :: nd_sat   = 0        !< of those, how many are at or beyond ion saturation, j/j_sat >= 1
  real*8  :: xd_min   =  1.d99   !< min/max of e*Phi/(k_B*Te) - Lambda over those rows
  real*8  :: xd_max   = -1.d99
  real*8  :: jd_min   =  1.d99   !< min/max of j/j_sat
  real*8  :: jd_max   = -1.d99
  real*8  :: sd_min   =  1.d99   !< min/max of the row scale, i.e. how far into a branch we are
  real*8  :: sd_max   = -1.d99
  real*8  :: rd_max   = -1.d99   !< max |residual| in volts, and where it sits
  real*8  :: rd_R     = 0.d0
  real*8  :: rd_Z     = 0.d0
  integer :: rd_type  = 0
  logical :: nd_first = .true.   !< run the normalisation self-test once, on the first active step

contains


!> Normalisation constants of the characteristic, in JOREK units. Artola eqs. (5) and (8).
!!
!! @param a_n    2*e*F0*sqrt(mu0*rho0)/m_i. The normalised potential is e*Phi/(k_B*Te) = a_n*u/(2*Te),
!!               with Te the JOREK electron temperature variable. Carries the sign of F0, so
!!               reversing the toroidal field reverses u and leaves Phi unchanged.
!! @param c_sat  -e*F0*n_0*sqrt(mu0/rho0), so that j_sat = c_sat*rho*v_par.
pure subroutine sheath_current_norm(a_n, c_sat)

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

end subroutine sheath_current_norm


!> The linearised sheath current row at one boundary node: Artola eq. (17).
!!
!! Returns the residual and the coefficients of the single scalar constraint
!!
!!     c_u*delta_u + c_T*delta_T + c_Te*delta_Te + c_Ti*delta_Ti + c_rho*delta_rho
!!                 + c_zj*delta_zj  =  - res ,
!!
!! already normalised so that max|c| = 1 (an exact rescaling of a constraint row, which is what
!! makes the ion-saturation limit finite). `scal` is the factor divided out, so `res*scal` and
!! `c_X*scal` are the raw row; the self-test needs this because the normalisation is itself a
!! function of the state. Use c_T in a single-temperature build and c_Te/c_Ti in a separate-
!! temperature build; all of them are always set, and the caller consumes the ones its model has.
!!
!! v_par is eliminated in favour of the Bohm value sign(B.n)*sqrt(GAMMA*(Ti+Te))/|B| (Artola eq. 16),
!! which is what the Mach BC enforces anyway. j_sat therefore depends on (rho, Ti, Te) and geometry
!! only, and this row never couples to the parallel velocity column.
!!
!! @param u0        potential variable u at the node
!! @param Ti0, Te0  ion and electron temperature at the node, already floored by the caller
!! @param rho0      density at the node
!! @param zj0       toroidal current variable zj at the node
!! @param Btot      |B| at the node
!! @param bn_sign   sign of B.n, outward positive (the boundary path's `direction`)
!! @param T_frozen  .true. if the caller's temperature floor is active, so that dTe = dTi = 0
subroutine sheath_current_row(u0, Ti0, Te0, rho0, zj0, Btot, bn_sign, T_frozen, &
                              res, c_u, c_T, c_Te, c_Ti, c_rho, c_zj, scal, xnorm, jratio)

  use phys_module, only: GAMMA, sheath_Lambda

  implicit none
  real*8,  intent(in)  :: u0, Ti0, Te0, rho0, zj0, Btot, bn_sign
  logical, intent(in)  :: T_frozen
  real*8,  intent(out) :: res, c_u, c_T, c_Te, c_Ti, c_rho, c_zj
  real*8,  intent(out) :: scal    !< factor the row was divided by; multiply back for the raw row
  real*8,  intent(out) :: xnorm   !< e*Phi/(k_B*Te) - Lambda, the normalised potential above floating
  real*8,  intent(out) :: jratio  !< j/j_sat at this node

  real*8 :: a_n, c_sat, T0, rho_s, jsat0, x0, expmx, f0, xi

  call sheath_current_norm(a_n, c_sat)

  T0    = Ti0 + Te0
  rho_s = sign(max(abs(rho0), rho_floor), rho0)

  ! --- Ion saturation current in the toroidal variable. The B.n of Artola eq. (3) cancels between
  ! --- the two sides, so j_sat carries no incidence factor - only the sign of B.n, which differs
  ! --- between the two divertor targets exactly as the physical current direction does.
  jsat0 = c_sat * rho_s * bn_sign * sqrt(GAMMA * T0) / Btot

  ! --- Normalised potential measured from the floating value, and the characteristic.
  x0    = a_n * u0 / (2.d0 * Te0) - sheath_Lambda
  x0    = max(-exp_arg_max, min(exp_arg_max, x0))    ! IEEE range guard, see module header
  expmx = exp(-x0)
  f0    = 1.d0 - expmx
  xi    = (2.d0 * Te0 / a_n) / expmx

  ! --- Artola eq. (17), with his single T split into Te (which sets the exponent) and Ti+Te (which
  ! --- sets c_s in j_sat). Setting Ti=Te=T/2 recovers c_T = (xi*f0/2 - u0)/T0 exactly.
  jratio = zj0 / jsat0
  xnorm  = x0
  res   = xi * (f0 - jratio)
  c_u   = 1.d0
  c_Te  = - u0 / Te0 + xi * f0 / (2.d0 * T0)
  c_Ti  =              xi * f0 / (2.d0 * T0)
  c_rho = xi * f0 / rho_s
  c_zj  = - xi / jsat0
  c_T   = 0.5d0 * (c_Te + c_Ti)

  ! --- If the caller's temperature floor is active, Te and Ti are constants at this node and the
  ! --- consistent derivative is zero. Anything else makes the value and slope rows disagree.
  if ( T_frozen ) then
    c_T = 0.d0 ; c_Te = 0.d0 ; c_Ti = 0.d0
  endif

  ! --- Exact row normalisation. A constraint row may be scaled by any positive number, and doing
  ! --- it here is load-bearing rather than cosmetic. The assembler ASSIGNS each entry as
  ! --- zbig*coefficient and leaves the columns nobody writes at their volume-assembled values, so
  ! --- the row is a replacement only as long as its largest entry dwarfs those leftovers. Fixing
  ! --- max|c| = 1 guarantees that at every state: without it the coefficients vanish with xi on
  ! --- the electron branch and the leftover volume terms would silently take the row back over,
  ! --- and they blow up with xi in deep ion saturation. It also gives the branch limits their
  ! --- correct form - delta_u = 0 on the electron side, a condition on delta_j when saturated.
  scal  = max(abs(c_u), abs(c_T), abs(c_Te), abs(c_Ti), abs(c_rho), abs(c_zj))
  res   = res   / scal
  c_u   = c_u   / scal
  c_T   = c_T   / scal
  c_Te  = c_Te  / scal
  c_Ti  = c_Ti  / scal
  c_rho = c_rho / scal
  c_zj  = c_zj  / scal

end subroutine sheath_current_row


!> Reconstruct the physical potential in VOLTS from u. One place, so that the conversion used to
!! impose the condition and the conversion used to report it cannot drift apart.
pure real*8 function sheath_current_volts(u)

  use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT
  use phys_module, only: F0, central_density, central_mass

  implicit none
  real*8, intent(in) :: u
  real*8 :: rho0

  rho0 = central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT
  sheath_current_volts = F0 * u / sqrt(MU_ZERO * rho0)

end function sheath_current_volts


!> Self-test of the characteristic and of its linearisation, run once at setup on rank 0.
!!
!! Checks, at the current namelist's normalisation:
!!   1. a_n carries the sign of F0, and Te of 1 eV floats at exactly sheath_Lambda volts;
!!   2. the zero-current state at the floating potential is a root of eq. (6);
!!   3. every coefficient equals a central finite difference of the residual, evaluated AWAY from
!!      the root so that no column can be accidentally right only where the residual vanishes;
!!   4. the single-temperature coefficient reduces to Artola eq. (17);
!!   5. the row stays finite 300 e-foldings into both branches.
!! Returns .false. and prints on failure. It does not stop; the caller decides.
logical function sheath_current_selftest(my_id)

  use constants,   only: MU_ZERO, EL_CHG
  use phys_module, only: F0, GAMMA, central_density, sheath_Lambda

  implicit none
  integer, intent(in) :: my_id

  real*8, parameter :: tol = 1.d-6
  real*8 :: a_n, c_sat, Te_1eV, Te0, Ti0, T0, rho0, Btot, bn, u_float, jsat0, u0, zj0
  real*8 :: res, c_u, c_T, c_Te, c_Ti, c_rho, c_zj, scal, xnorm, jratio, pref

  sheath_current_selftest = .true.
  call sheath_current_norm(a_n, c_sat)

  ! --- A representative boundary state: 10 eV, unit density, unit field, normal incidence sign +1.
  Te_1eV = EL_CHG * MU_ZERO * central_density * 1.d20
  Te0    = 10.d0 * Te_1eV
  Ti0    = Te0
  T0     = Ti0 + Te0
  rho0   = 1.d0
  Btot   = 1.d0
  bn     = 1.d0

  ! --- 1. sign, and the floating level in volts per eV
  if ( a_n * F0 .le. 0.d0 ) then
    sheath_current_selftest = .false.
    if (my_id .eq. 0) write(*,'(A,2es22.14)') ' SHEATH_J SELFTEST FAIL: sign(a_n) /= sign(F0):', a_n, F0
  endif
  u_float = 2.d0 * sheath_Lambda * Te0 / a_n
  if ( abs(sheath_current_volts(u_float) - 10.d0*sheath_Lambda) .gt. tol * 10.d0*abs(sheath_Lambda) ) then
    sheath_current_selftest = .false.
    if (my_id .eq. 0) write(*,'(A,2es22.14)') ' SHEATH_J SELFTEST FAIL: floating level [V]:', &
      sheath_current_volts(u_float), 10.d0*sheath_Lambda
  endif

  ! --- 2. the floating, zero-current state must be a root
  if ( abs(gcur(u_float, Ti0, Te0, rho0, 0.d0)) .gt. tol * abs(c_sat) ) then
    sheath_current_selftest = .false.
    if (my_id .eq. 0) write(*,'(A,es22.14)') ' SHEATH_J SELFTEST FAIL: floating state not a root:', &
      gcur(u_float, Ti0, Te0, rho0, 0.d0)
  endif

  ! --- 3. finite-difference check of every column, at a state OFF the characteristic
  jsat0 = c_sat * rho0 * bn * sqrt(GAMMA * T0) / Btot
  u0    = 0.7d0 * u_float
  zj0   = 0.3d0 * jsat0
  ! ---    What is differentiated is Artola's CURRENT residual G = j - j_sat*f, coded here
  ! ---    independently from eq. (6), and each column is compared against pref*dG/dx where
  ! ---    pref = -xi/j_sat0 is taken from the module itself as c_zj*scal (exact, because
  ! ---    dG/dj = 1). This checks eqs. (12)-(17) against eq. (6) with one scalar in common.
  call sheath_current_row(u0, Ti0, Te0, rho0, zj0, Btot, bn, .false., &
                          res, c_u, c_T, c_Te, c_Ti, c_rho, c_zj, scal, xnorm, jratio)
  pref = c_zj * scal

  call check_col('u  ', c_u  *scal, pref*fd_u  (u0, Ti0, Te0, rho0, zj0))
  call check_col('Te ', c_Te *scal, pref*fd_Te (u0, Ti0, Te0, rho0, zj0))
  call check_col('Ti ', c_Ti *scal, pref*fd_Ti (u0, Ti0, Te0, rho0, zj0))
  call check_col('rho', c_rho*scal, pref*fd_rho(u0, Ti0, Te0, rho0, zj0))

  ! ---    and the residual itself must be pref*G0
  if ( abs(res*scal - pref*gcur(u0, Ti0, Te0, rho0, zj0)) .gt. 1.d-10 * abs(res*scal) ) then
    sheath_current_selftest = .false.
    if (my_id .eq. 0) write(*,'(A,2es22.14)') ' SHEATH_J SELFTEST FAIL: residual /= pref*G0:', &
      res*scal, pref*gcur(u0, Ti0, Te0, rho0, zj0)
  endif

  ! --- 4. single-temperature reduction, Artola eq. (17)
  if ( abs(c_T - 0.5d0*(c_Te + c_Ti)) .gt. tol * max(abs(c_T), 1.d0) ) then
    sheath_current_selftest = .false.
    if (my_id .eq. 0) write(*,'(A,2es22.14)') ' SHEATH_J SELFTEST FAIL: c_T /= (c_Te+c_Ti)/2:', &
      c_T, 0.5d0*(c_Te + c_Ti)
  endif

  ! --- 5. both branches must stay finite far from the characteristic
  call check_finite('ion saturation branch',  300.d0*u_float, Ti0, Te0, rho0, zj0)
  call check_finite('electron branch      ', -300.d0*u_float, Ti0, Te0, rho0, zj0)

  if ( my_id .eq. 0 ) then
    write(*,'(A)')         ' --- sheath current BC normalisation ---'
    write(*,'(A,es22.14)') '   a_n                        = ', a_n
    write(*,'(A,es22.14)') '   c_sat                      = ', c_sat
    write(*,'(A,f10.4)')   '   sheath_Lambda              = ', sheath_Lambda
    write(*,'(A,es22.14)') '   volts per unit u           = ', sheath_current_volts(1.d0)
    write(*,'(A,es22.14)') '   floating u at Te = 10 eV   = ', u_float
    if ( sheath_current_selftest ) then
      write(*,'(A)')       '   selftest PASSED'
    else
      write(*,'(A)')       '   selftest FAILED - do not run'
    endif
  endif

contains

  !> Artola eq. (6) written out independently of the production row: the current residual
  !! G = j - j_sat*(1 - exp(Lambda - e*Phi/(k_B*Te))). This is what the row linearises.
  real*8 function gcur(u, Ti, Te, rho, zj)
    real*8, intent(in) :: u, Ti, Te, rho, zj
    real*8 :: js, ff
    js   = c_sat * rho * bn * sqrt(GAMMA*(Ti+Te)) / Btot
    ff   = 1.d0 - exp(sheath_Lambda - a_n*u/(2.d0*Te))
    gcur = zj - js*ff
  end function gcur

  real*8 function fd_u(u, Ti, Te, rho, zj)
    real*8, intent(in) :: u, Ti, Te, rho, zj
    real*8 :: h
    h = 1.d-7 * abs(u)
    fd_u = (gcur(u+h, Ti, Te, rho, zj) - gcur(u-h, Ti, Te, rho, zj)) / (2.d0*h)
  end function fd_u

  real*8 function fd_Te(u, Ti, Te, rho, zj)
    real*8, intent(in) :: u, Ti, Te, rho, zj
    real*8 :: h
    h = 1.d-7 * Te
    fd_Te = (gcur(u, Ti, Te+h, rho, zj) - gcur(u, Ti, Te-h, rho, zj)) / (2.d0*h)
  end function fd_Te

  real*8 function fd_Ti(u, Ti, Te, rho, zj)
    real*8, intent(in) :: u, Ti, Te, rho, zj
    real*8 :: h
    h = 1.d-7 * Ti
    fd_Ti = (gcur(u, Ti+h, Te, rho, zj) - gcur(u, Ti-h, Te, rho, zj)) / (2.d0*h)
  end function fd_Ti

  real*8 function fd_rho(u, Ti, Te, rho, zj)
    real*8, intent(in) :: u, Ti, Te, rho, zj
    real*8 :: h
    h = 1.d-7 * rho
    fd_rho = (gcur(u, Ti, Te, rho+h, zj) - gcur(u, Ti, Te, rho-h, zj)) / (2.d0*h)
  end function fd_rho

  subroutine check_col(name, c, d)
    character(*), intent(in) :: name
    real*8,       intent(in) :: c, d
    if ( abs(c - d) .gt. 1.d-5 * max(abs(c), abs(d), 1.d-30) ) then
      sheath_current_selftest = .false.
      if (my_id .eq. 0) write(*,'(A,A,A,2es22.14)') &
        ' SHEATH_J SELFTEST FAIL: column ', name, ' /= finite difference:', c, d
    endif
  end subroutine check_col

  subroutine check_finite(name, u, Ti, Te, rho, zj)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    character(*), intent(in) :: name
    real*8,       intent(in) :: u, Ti, Te, rho, zj
    real*8 :: r, b1, b2, b3, b4, b5, b6, sc, xn, jr
    call sheath_current_row(u, Ti, Te, rho, zj, Btot, bn, .false., r, b1, b2, b3, b4, b5, b6, sc, xn, jr)
    if ( .not. (ieee_is_finite(r)  .and. ieee_is_finite(b1) .and. ieee_is_finite(b2) .and. &
                ieee_is_finite(b3) .and. ieee_is_finite(b4) .and. ieee_is_finite(b5) .and. &
                ieee_is_finite(b6) .and. ieee_is_finite(sc)) ) then
      sheath_current_selftest = .false.
      if (my_id .eq. 0) write(*,'(A,A)') ' SHEATH_J SELFTEST FAIL: non-finite row on the ', name
    endif
  end subroutine check_finite

end function sheath_current_selftest


!> Clear the monitoring accumulators. Called once at the start of each matrix construction.
subroutine sheath_diag_reset()
  implicit none
  nd_rows = 0 ; nd_sat = 0
  xd_min  =  1.d99 ; xd_max = -1.d99
  jd_min  =  1.d99 ; jd_max = -1.d99
  sd_min  =  1.d99 ; sd_max = -1.d99
  rd_max  = -1.d99 ; rd_R = 0.d0 ; rd_Z = 0.d0 ; rd_type = 0
end subroutine sheath_diag_reset


!> Record one assembled sheath row.
subroutine sheath_diag_add(res, scal, xnorm, jratio, BigR, Z, bnd_type)
  implicit none
  real*8,  intent(in) :: res, scal, xnorm, jratio, BigR, Z
  integer, intent(in) :: bnd_type
  real*8 :: rv

  nd_rows = nd_rows + 1
  if ( jratio .ge. 1.d0 ) nd_sat = nd_sat + 1
  xd_min  = min(xd_min, xnorm)  ; xd_max = max(xd_max, xnorm)
  jd_min  = min(jd_min, jratio) ; jd_max = max(jd_max, jratio)
  sd_min  = min(sd_min, scal)   ; sd_max = max(sd_max, scal)

  ! --- The residual is a potential, so report it in volts: that is the number to watch, since it
  ! --- is how far this node's potential is from the value the characteristic demands.
  rv = abs(sheath_current_volts(res))
  if ( rv .gt. rd_max ) then
    rd_max = rv ; rd_R = BigR ; rd_Z = Z ; rd_type = bnd_type
  endif
end subroutine sheath_diag_add


!> Reduce the monitoring accumulators over all ranks and print two lines on rank 0.
!!
!! Reading the output when a run is going wrong:
!!   * `ePhi/Te-Lam` is the normalised potential measured from floating. Near zero everywhere means
!!     the wall sits close to floating, i.e. little net current. Strongly positive at a node means
!!     it is heading into ion saturation, strongly negative means the electron branch. The two
!!     divertor targets ending up on opposite sides is the thermoelectric current and is physical.
!!   * `j/jsat` above 1 is the sheath being asked to pass more ion current than it can. `nsat` counts
!!     those nodes. A handful is the expected saturation; a whole target is not.
!!   * `rowscale` is the factor divided out of each row. It grows like exp(ePhi/Te), so its spread
!!     measures how far apart the branches are. A scale running away on its own is not physical.
!!   * `max|res|` is how far the worst node's potential is from the characteristic, in volts, with
!!     that node's position and boundary type. This is the number that grows before a crash, and the
!!     boundary type it reports says which node family is losing the condition first.
subroutine sheath_diag_print(my_id)

  use mpi_mod

  implicit none
  integer, intent(in) :: my_id

  integer :: ierr, cnt_in(2), cnt_out(2)
  real*8  :: mn_in(3), mn_out(3), mx_in(3), mx_out(3)
  real*8  :: loc_in(2), loc_out(2), payload(3)

  cnt_in = (/ nd_rows, nd_sat /)
  call MPI_ALLREDUCE(cnt_in, cnt_out, 2, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  if ( cnt_out(1) .eq. 0 ) return      ! collective: every rank sees the same total and returns

  ! --- Check the normalisation once, the first time the condition is actually assembled: by now
  ! --- F0 and the central density/mass are the values the run will use, which is not true at
  ! --- namelist time. A failure here is a wrong potential scale, so it must not be survivable.
  if ( nd_first ) then
    nd_first = .false.
    if ( .not. sheath_current_selftest(my_id) ) then
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    endif
  endif

  mn_in = (/ xd_min, jd_min, sd_min /)
  mx_in = (/ xd_max, jd_max, sd_max /)
  call MPI_ALLREDUCE(mn_in, mn_out, 3, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(mx_in, mx_out, 3, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)

  ! --- Find the rank owning the worst residual, then take its location from that rank.
  loc_in = (/ rd_max, dble(my_id) /)
  call MPI_ALLREDUCE(loc_in, loc_out, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD, ierr)
  payload = (/ rd_R, rd_Z, dble(rd_type) /)
  call MPI_BCAST(payload, 3, MPI_REAL8, nint(loc_out(2)), MPI_COMM_WORLD, ierr)

  if ( my_id .eq. 0 ) then
    write(*,'(A,I7,A,I6,A,2es11.3,A,2es11.3)')                            &
      ' [sheath_j] rows=', cnt_out(1), ' nsat=', cnt_out(2),              &
      '  ePhi/Te-Lam=', mn_out(1), mx_out(1), '  j/jsat=', mn_out(2), mx_out(2)
    write(*,'(A,2es11.3,A,es11.3,A,2f9.4,A,I3)')                          &
      ' [sheath_j] rowscale=', mn_out(3), mx_out(3),                      &
      '  max|res|[V]=', loc_out(1),                                       &
      '  at R,Z=', payload(1), payload(2), '  bnd type', nint(payload(3))
  endif

end subroutine sheath_diag_print


end module mod_sheath_current
