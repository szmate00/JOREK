!> Sheath current boundary condition for the electrostatic potential (model600).
!!
!! Imposes the sheath current-voltage characteristic (Stangeby 2.68; Artola eq. 6)
!!
!!     j = j_sat * ( 1 - exp( Lambda - e*Phi/(k_B*Te) ) )
!!
!! at a material wall, as ONE scalar constraint per boundary node replacing the u (potential) node
!! value row and its boundary-derivative row. j is the JOREK toroidal current variable zj and
!! j_sat = c_sat*rho*Vpar is the ion saturation current in the same variable (Artola eqs. 1-5).
!! Lambda = sheath_Lambda is the floating sheath drop in units of Te/e.
!!
!! Source: Artola, "Sheath boundary conditions for the electric potential in JOREK" (2026-07-30),
!! `potential_BC_JOREK_ionsat_v2.pdf`. The characteristic and j_sat are his eqs. (5) and (6)
!! unchanged. The LINEARISATION here is deliberately NOT his eq. (17) - see "why not eq. 17" below.
!!
!! ------------------------------------------------------------------------------------------------
!! THE TWO BRANCHES
!!
!! Solving the characteristic for the potential gives, wherever a voltage root exists,
!!
!!     u = (2*Te/a_n) * [ Lambda - ln X ] ,        X = 1 - j/j_sat > 0 ,
!!
!! with e*Phi/(k_B*Te) = a_n*u/(2*Te). Residual POT: res = u - (2*Te/a_n)*(Lambda - ln X).
!! Every column of this is exact, `d(res)/du = 1` identically, and - the property that matters -
!! the Newton update IS the voltage defect: a node 100 V from the characteristic asks for 100 V,
!! not for exp(100 e/Te).
!!
!! X <= 0 means the plasma is demanding at least the full ion saturation current, and the
!! characteristic then has no voltage root at all: the sheath is saturated. There the physical
!! statement is a condition on the CURRENT, residual SAT: res = j - j_sat, with no u column.
!!
!! The branch test is written as `(j_sat - j)*j_sat > 0`, which is X > 0 without dividing, so it
!! also selects the saturated branch when j_sat is exactly zero (no ion outflow, no sheath current
!! capacity, hence j = 0). This is the only branch in the module and it is a physical boundary, not
!! a tuned switch.
!!
!! THE BRANCHES JOIN. Each row is normalised by its largest coefficient - an exact rescaling of a
!! constraint row, needed anyway for the assembler reason below. As X -> 0+ the POT coefficients
!! c_rho, c_Vpar and c_j all diverge like 1/X while c_u and c_Te stay O(1), so after normalisation
!! c_u, c_Te -> 0 and the surviving three tend to the ratios (1/rho : 1/Vpar : -1/j_sat), which are
!! exactly the normalised SAT coefficients. Both normalised residuals tend to zero. So the switch
!! is continuous, and a node oscillating about saturation does not see two unrelated rows.
!!
!! ------------------------------------------------------------------------------------------------
!! WHY NOT ARTOLA EQ. (17)
!!
!! Eq. (17) linearises the CURRENT residual G = j - j_sat*f and then multiplies the linear equation
!! by a frozen -xi/j_sat, xi = (2*Te/a_n)*exp(e*Phi/kTe - Lambda). That makes the u column exactly 1
!! and looks well conditioned, but multiplying a residual and its Jacobian by the same frozen
!! scalar cannot change the Newton update. Solving that row with j held fixed gives
!!
!!     delta_x = 1 - exp(x) ,      x = e*Phi/(k_B*Te) - Lambda ,
!!
!! so a state 100 V above the floating root at Te = 10 eV requests a 220 kV step. That was measured
!! with the previous version of this module, not argued. The POT residual above is the same
!! characteristic written so that the update is bounded by construction.
!!
!! WHY j_sat USES THE ACTUAL Vpar AND NOT THE BOHM VALUE. Substituting Vpar = sign(B.n)*cs/|B|
!! (Artola eq. 16) makes j_sat depend on Ti+Te and on |B|. Two things go wrong. |B| varies ALONG
!! the wall, so the boundary-derivative row then needs a d|B|/dl term, and omitting it drives a
!! spurious potential gradient of order a volt per unit boundary coordinate. And the substitution
!! asserts what the Mach BC imposes, which is not what it imposes: the nodal Mach row also carries
!! a drift term factor/Btot*R^2*u_b/psi_b. Using j_sat = c_sat*rho*Vpar removes the geometry from
!! the residual entirely - so the derivative row below is exactly d/dl of the value row with no
!! missing term - and stops asserting anything about Vpar.
!!
!! The price is a stated domain restriction rather than a hidden one: the characteristic presumes
!! ion OUTFLOW at the wall, Vpar*(B.n) > 0, which is what `bcs(i)%mach1` provides. The condition
!! therefore requires with_vpar, checked at the assembly site. Where the outflow vanishes, j_sat
!! vanishes, the saturated branch is selected and the row correctly stops constraining u.
!!
!! ------------------------------------------------------------------------------------------------
!! WHAT IS LAGGED, AND WHAT IS NOT
!!
!! The value row is the exact Newton row of the residual: every column is the analytic derivative,
!! finite-difference checked against the characteristic coded independently from eq. (6).
!!
!! The boundary-derivative row applies the value row's coefficient vector to the derivative DOFs.
!! Because the residual contains no geometry, that IS d/dl of the value condition exactly. Its
!! JACOBIAN, however, omits the value-DOF columns that come from differentiating the state-
!! dependent coefficients: it is a Picard lag of a term which is itself a gradient, not an exact
!! Jacobian, and it is not claimed to be one.
!!
!! WHAT IS NOT HERE. No ratio gate, no minimum-|B.n|, no saturation-slope limiter, no current clip,
!! no exponent clamp, no relaxation gain, no frame-determinant threshold, no density or temperature
!! floor of its own, and no exponential anywhere - the logarithm is the only transcendental, and its
!! argument is positive by the branch test. sheath_Lambda is the only input.
module mod_sheath_current

  implicit none
  private

  public :: sheath_current_norm, sheath_current_row, sheath_current_volts, sheath_current_selftest
  public :: sheath_diag_reset, sheath_diag_add, sheath_diag_print

  ! --- Monitoring state, accumulated over one matrix construction and printed once per build.
  ! --- Nothing here feeds back into the row.
  integer :: nd_rows  = 0        !< sheath rows assembled on LOCALLY OWNED nodes
  integer :: nd_sat   = 0        !< of those, how many took the saturated branch
  real*8  :: xd_min   =  1.d99   !< min/max of e*Phi/(k_B*Te) - Lambda
  real*8  :: xd_max   = -1.d99
  real*8  :: jd_min   =  1.d99   !< min/max of j/j_sat
  real*8  :: jd_max   = -1.d99
  real*8  :: vd_max   = -1.d99   !< max voltage defect on the unsaturated branch, and where
  real*8  :: vd_R     = 0.d0
  real*8  :: vd_Z     = 0.d0
  integer :: vd_type  = 0
  logical :: nd_first = .true.   !< run the normalisation self-test once, on the first active build

contains


!> Normalisation constants of the characteristic, in JOREK units. Artola eqs. (5) and (8).
!!
!! @param a_n    2*e*F0*sqrt(mu0*rho0)/m_i, so that e*Phi/(k_B*Te) = a_n*u/(2*Te) with Te the JOREK
!!               electron temperature variable. Carries the sign of F0, so reversing the toroidal
!!               field reverses u and leaves the physical potential unchanged.
!! @param c_sat  -e*F0*n_0*sqrt(mu0/rho0), so that j_sat = c_sat*rho*Vpar.
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


!> The sheath row at one boundary node.
!!
!! Returns the residual and the coefficients of the single scalar constraint
!!
!!     c_u*du + c_T*dT + c_Te*dTe + c_rho*drho + c_vpar*dVpar + c_zj*dzj  =  - res ,
!!
!! normalised so that max|c| = 1. Use c_T in a single-temperature build and c_Te in a separate-
!! temperature build; both are always set. There is no Ti column: with j_sat taken from Vpar, the
!! residual does not contain Ti at all.
!!
!! `scal` is the factor divided out, so `res_raw = res*scal` and `c_X*scal` is the raw coefficient;
!! the self-test needs it because the normalisation is itself a function of the state.
!! `res_raw` is the residual before normalisation. On the unsaturated branch it is the VOLTAGE
!! DEFECT in u units, i.e. u minus the potential the characteristic requires, which is the only
!! quantity in here that may honestly be reported in volts. On the saturated branch it is a current
!! and `saturated` is .true., so the caller must not convert it.
!!
!! @param u0         potential variable u at the node
!! @param Te0        electron temperature at the node (T/2 in a single-temperature build)
!! @param rho0       density at the node
!! @param vpar0      parallel velocity coefficient at the node (Vpar, not a speed)
!! @param zj0        toroidal current variable zj at the node
!! @param Te_frozen  .true. if the caller's temperature floor is active, so that dTe = 0
subroutine sheath_current_row(u0, Te0, rho0, vpar0, zj0, Te_frozen,           &
                              res, c_u, c_T, c_Te, c_rho, c_vpar, c_zj,       &
                              res_raw, scal, saturated, xnorm, jratio)

  use phys_module, only: sheath_Lambda

  implicit none
  real*8,  intent(in)  :: u0, Te0, rho0, vpar0, zj0
  logical, intent(in)  :: Te_frozen
  real*8,  intent(out) :: res, c_u, c_T, c_Te, c_rho, c_vpar, c_zj
  real*8,  intent(out) :: res_raw    !< residual before normalisation; a potential if .not. saturated
  real*8,  intent(out) :: scal       !< factor the row was divided by: res_raw = res*scal, c_X_raw = c_X*scal
  logical, intent(out) :: saturated  !< .true. if the characteristic has no voltage root here
  real*8,  intent(out) :: xnorm      !< e*Phi/(k_B*Te) - Lambda, the normalised potential above floating
  real*8,  intent(out) :: jratio     !< j/j_sat, or zero where j_sat vanishes

  real*8 :: a_n, c_sat, jsat0, X, w

  call sheath_current_norm(a_n, c_sat)

  jsat0 = c_sat * rho0 * vpar0
  xnorm = a_n * u0 / (2.d0 * Te0) - sheath_Lambda

  ! --- X > 0 written without a division, so j_sat = 0 falls to the saturated branch.
  saturated = .not. ( (jsat0 - zj0) * jsat0 .gt. 0.d0 )

  if ( saturated ) then

    ! --- No voltage root: impose the current limit j = j_sat. No u and no Te column.
    jratio  = 0.d0
    if ( jsat0 .ne. 0.d0 ) jratio = zj0 / jsat0
    res_raw = zj0 - jsat0
    c_u     = 0.d0
    c_Te    = 0.d0
    c_rho   = - c_sat * vpar0
    c_vpar  = - c_sat * rho0
    c_zj    =   1.d0

  else

    ! --- Unsaturated: the characteristic solved for the potential. rho0 and vpar0 are both nonzero
    ! --- here, because j_sat is, so none of the divisions below need a guard.
    jratio  = zj0 / jsat0
    X       = 1.d0 - jratio
    w       = 2.d0 * Te0 / a_n
    res_raw = u0 - w * (sheath_Lambda - log(X))
    c_u     =   1.d0
    c_Te    = - (2.d0 / a_n) * (sheath_Lambda - log(X))
    c_rho   =   w * jratio / (rho0  * X)
    c_vpar  =   w * jratio / (vpar0 * X)
    c_zj    = - w / (X * jsat0)

  endif

  ! --- If the caller's temperature floor is active, Te is a constant at this node and the
  ! --- consistent column is zero. Anything else makes the value and derivative rows disagree.
  if ( Te_frozen ) c_Te = 0.d0

  ! --- Single-temperature build: the model evolves T = Ti + Te with Te = T/2, so dTe = dT/2.
  c_T = 0.5d0 * c_Te

  ! --- Exact row normalisation, load-bearing for two separate reasons. First, the assembler
  ! --- ASSIGNS each entry as zbig*coefficient and leaves the columns nobody writes at their
  ! --- volume-assembled values, so the row replaces the equation only while its largest entry
  ! --- dwarfs those leftovers; fixing max|c| = 1 guarantees that at every state. Second, it is
  ! --- what makes the two branches join continuously at X = 0, as derived in the module header.
  scal  = max(abs(c_u), abs(c_T), abs(c_Te), abs(c_rho), abs(c_vpar), abs(c_zj))
  res   = res_raw / scal
  c_u   = c_u   / scal
  c_T   = c_T   / scal
  c_Te  = c_Te  / scal
  c_rho = c_rho / scal
  c_vpar= c_vpar/ scal
  c_zj  = c_zj  / scal

end subroutine sheath_current_row


!> Reconstruct a physical potential in VOLTS from a quantity in u units. One place, so that the
!! conversion used to impose the condition and the conversion used to report it cannot drift apart.
pure real*8 function sheath_current_volts(u)

  use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT
  use phys_module, only: F0, central_density, central_mass

  implicit none
  real*8, intent(in) :: u
  real*8 :: rho0

  rho0 = central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT
  sheath_current_volts = F0 * u / sqrt(MU_ZERO * rho0)

end function sheath_current_volts


!> Clear the monitoring accumulators. Called once at the start of each matrix construction.
subroutine sheath_diag_reset()
  implicit none
  nd_rows = 0 ; nd_sat = 0
  xd_min  =  1.d99 ; xd_max = -1.d99
  jd_min  =  1.d99 ; jd_max = -1.d99
  vd_max  = -1.d99 ; vd_R = 0.d0 ; vd_Z = 0.d0 ; vd_type = 0
end subroutine sheath_diag_reset


!> Record one assembled sheath row. `owned` must be the caller's local-ownership test, so that a
!! node visited from several elements or from a halo is counted once and only where it is solved.
subroutine sheath_diag_add(res_raw, saturated, xnorm, jratio, BigR, Z, bnd_type, owned)
  implicit none
  real*8,  intent(in) :: res_raw, xnorm, jratio, BigR, Z
  logical, intent(in) :: saturated, owned
  integer, intent(in) :: bnd_type
  real*8 :: vdef

  if ( .not. owned ) return

  nd_rows = nd_rows + 1
  xd_min  = min(xd_min, xnorm)  ; xd_max = max(xd_max, xnorm)
  jd_min  = min(jd_min, jratio) ; jd_max = max(jd_max, jratio)

  if ( saturated ) then
    ! --- res_raw is a current here and there is no voltage root, so it is NOT converted.
    nd_sat = nd_sat + 1
  else
    vdef = abs(sheath_current_volts(res_raw))
    if ( vdef .gt. vd_max ) then
      vd_max = vdef ; vd_R = BigR ; vd_Z = Z ; vd_type = bnd_type
    endif
  endif
end subroutine sheath_diag_add


!> Reduce the monitoring accumulators over all ranks and print two lines on rank 0.
!!
!! Reading the output when a run is going wrong:
!!   * `ePhi/Te-Lam` is the normalised potential measured from floating. Near zero everywhere means
!!     the wall sits close to floating, i.e. little net current. The two divertor targets ending up
!!     on opposite signs is the thermoelectric current and is physical.
!!   * `j/jsat` at or above 1 is the sheath being asked to pass at least the full ion saturation
!!     current. `nsat` counts the nodes that took the saturated branch, where the condition becomes
!!     j = j_sat and stops constraining the potential. A handful is expected; a whole target means
!!     the potential is no longer being set by the sheath there.
!!   * `max Vdef` is the largest voltage defect - u minus the potential the characteristic requires,
!!     converted from the UNNORMALISED residual, so it is a genuine voltage - together with the
!!     position and boundary type of that node. This is the number that grows before a crash, and
!!     the boundary type says which node family is losing the condition first. Nodes on the
!!     saturated branch are excluded, because there is no voltage root there to be distant from.
subroutine sheath_diag_print(my_id)

  use mpi_mod

  implicit none
  integer, intent(in) :: my_id

  integer :: ierr, cnt_in(2), cnt_out(2)
  real*8  :: mn_in(2), mn_out(2), mx_in(2), mx_out(2)
  real*8  :: loc_in(2), loc_out(2), payload(3)

  cnt_in = (/ nd_rows, nd_sat /)
  call MPI_ALLREDUCE(cnt_in, cnt_out, 2, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  if ( cnt_out(1) .eq. 0 ) return      ! collective: every rank sees the same total and returns

  ! --- Check the normalisation once, the first time the condition is actually assembled: by now
  ! --- F0 and the central density/mass are the values the run will use, which is not true at
  ! --- namelist time. A failure here is a wrong potential scale, so it must not be survivable.
  if ( nd_first ) then
    nd_first = .false.
    if ( .not. sheath_current_selftest(my_id) ) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  endif

  mn_in = (/ xd_min, jd_min /)
  mx_in = (/ xd_max, jd_max /)
  call MPI_ALLREDUCE(mn_in, mn_out, 2, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(mx_in, mx_out, 2, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)

  ! --- Find the rank owning the worst voltage defect, then take its location from that rank.
  loc_in = (/ vd_max, dble(my_id) /)
  call MPI_ALLREDUCE(loc_in, loc_out, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD, ierr)
  payload = (/ vd_R, vd_Z, dble(vd_type) /)
  call MPI_BCAST(payload, 3, MPI_REAL8, nint(loc_out(2)), MPI_COMM_WORLD, ierr)

  if ( my_id .eq. 0 ) then
    write(*,'(A,I7,A,I6,A,2es11.3,A,2es11.3)')                            &
      ' [sheath_j] rows=', cnt_out(1), ' nsat=', cnt_out(2),              &
      '  ePhi/Te-Lam=', mn_out(1), mx_out(1), '  j/jsat=', mn_out(2), mx_out(2)
    if ( cnt_out(1) .gt. cnt_out(2) ) then
      write(*,'(A,es11.3,A,2f9.4,A,I3)')                                  &
        ' [sheath_j] max Vdef [V]=', loc_out(1),                          &
        '  at R,Z=', payload(1), payload(2), '  bnd type', nint(payload(3))
    else
      write(*,'(A)') ' [sheath_j] every row saturated: no voltage root anywhere on the wall'
    endif
  endif

end subroutine sheath_diag_print


!> Self-test of the characteristic and of its linearisation, run once on the first active build.
!!
!! Checks, at the current namelist's normalisation:
!!   1. a_n carries the sign of F0, and Te of 1 eV floats at exactly sheath_Lambda volts;
!!   2. the zero-current state at the floating potential is a root of eq. (6), and the unsaturated
!!      residual there is zero;
!!   3. the unsaturated residual IS the voltage defect: displacing u by a known number of volts
!!      must return exactly that number;
!!   4. every unsaturated column equals a central finite difference of the residual;
!!   5. every saturated column equals a central finite difference of j - j_sat;
!!   6. the two branches join: approaching X = 0 from below, the normalised rows and residuals of
!!      the two branches converge.
!! Returns .false. and prints on failure.
logical function sheath_current_selftest(my_id)

  use constants,   only: MU_ZERO, EL_CHG
  use phys_module, only: F0, central_density, sheath_Lambda

  implicit none
  integer, intent(in) :: my_id

  real*8, parameter :: tol = 1.d-8
  real*8 :: a_n, c_sat, Te_1eV, Te0, rho0, vpar0, u_float, jsat0, u0, zj0, dV
  real*8 :: res, c_u, c_T, c_Te, c_rho, c_vpar, c_zj, res_raw, scal, xnorm, jratio
  logical :: sat

  sheath_current_selftest = .true.
  call sheath_current_norm(a_n, c_sat)

  ! --- A representative boundary state: 10 eV, unit density, outgoing sonic-ish Vpar.
  Te_1eV  = EL_CHG * MU_ZERO * central_density * 1.d20
  Te0     = 10.d0 * Te_1eV
  rho0    = 1.d0
  vpar0   = 3.d-2
  u_float = 2.d0 * sheath_Lambda * Te0 / a_n
  jsat0   = c_sat * rho0 * vpar0

  ! --- 1. sign, and the floating level in volts
  if ( a_n * F0 .le. 0.d0 ) call fail_2('sign(a_n) /= sign(F0)', a_n, F0)
  if ( abs(sheath_current_volts(u_float) - 10.d0*sheath_Lambda) .gt. tol * 10.d0*abs(sheath_Lambda) ) &
    call fail_2('floating level [V]', sheath_current_volts(u_float), 10.d0*sheath_Lambda)

  ! --- 2. the floating, zero-current state is a root of eq. (6) and of the residual
  if ( abs(gcur(u_float, Te0, rho0, vpar0, 0.d0)) .gt. tol * abs(jsat0) ) &
    call fail_1('floating state not a root of eq. 6', gcur(u_float, Te0, rho0, vpar0, 0.d0))
  call sheath_current_row(u_float, Te0, rho0, vpar0, 0.d0, .false., &
                          res, c_u, c_T, c_Te, c_rho, c_vpar, c_zj, res_raw, scal, sat, xnorm, jratio)
  if ( sat ) call fail_1('floating state reported saturated', 0.d0)
  if ( abs(res_raw) .gt. tol * abs(u_float) ) call fail_1('residual nonzero at the root', res_raw)

  ! --- 3. the residual IS the voltage defect. Displace u by 137 V and demand 137 V back.
  dV = 137.d0
  u0 = u_float + dV / sheath_current_volts(1.d0)
  call sheath_current_row(u0, Te0, rho0, vpar0, 0.d0, .false., &
                          res, c_u, c_T, c_Te, c_rho, c_vpar, c_zj, res_raw, scal, sat, xnorm, jratio)
  if ( abs(sheath_current_volts(res_raw) - dV) .gt. tol * dV ) &
    call fail_2('residual is not the voltage defect', sheath_current_volts(res_raw), dV)

  ! --- 4. finite-difference check of every unsaturated column, at a state OFF the characteristic
  u0  = 0.7d0 * u_float
  zj0 = 0.3d0 * jsat0
  call sheath_current_row(u0, Te0, rho0, vpar0, zj0, .false., &
                          res, c_u, c_T, c_Te, c_rho, c_vpar, c_zj, res_raw, scal, sat, xnorm, jratio)
  if ( sat ) call fail_1('test state unexpectedly saturated', zj0)
  ! ---    The raw row is compared against the raw residual: the normalisation factor is itself a
  ! ---    function of the state, so normalised coefficients and a differenced residual would be
  ! ---    two different scalings.
  call check('POT u   ', c_u   *scal, fd(1, u0, Te0, rho0, vpar0, zj0))
  call check('POT Te  ', c_Te  *scal, fd(2, u0, Te0, rho0, vpar0, zj0))
  call check('POT rho ', c_rho *scal, fd(3, u0, Te0, rho0, vpar0, zj0))
  call check('POT vpar', c_vpar*scal, fd(4, u0, Te0, rho0, vpar0, zj0))
  call check('POT zj  ', c_zj  *scal, fd(5, u0, Te0, rho0, vpar0, zj0))

  ! --- 5. the same on the saturated branch, where the residual is j - j_sat
  zj0 = 1.4d0 * jsat0
  call sheath_current_row(u0, Te0, rho0, vpar0, zj0, .false., &
                          res, c_u, c_T, c_Te, c_rho, c_vpar, c_zj, res_raw, scal, sat, xnorm, jratio)
  if ( .not. sat ) call fail_1('state beyond saturation not reported saturated', zj0)
  call check('SAT rho ', c_rho *scal, fd(3, u0, Te0, rho0, vpar0, zj0))
  call check('SAT vpar', c_vpar*scal, fd(4, u0, Te0, rho0, vpar0, zj0))
  call check('SAT zj  ', c_zj  *scal, fd(5, u0, Te0, rho0, vpar0, zj0))
  if ( abs(c_u) + abs(c_Te) .ne. 0.d0 ) call fail_1('saturated row still constrains u or Te', c_u)

  ! --- 6. the branches join. Approach X = 0 from below and compare the normalised rows.
  call check_join(1.d-12)

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
  !! G = j - j_sat*(1 - exp(Lambda - e*Phi/(k_B*Te))), with j_sat = c_sat*rho*Vpar from eq. (5).
  real*8 function gcur(u, Te, rho, vp, zj)
    real*8, intent(in) :: u, Te, rho, vp, zj
    gcur = zj - c_sat*rho*vp * (1.d0 - exp(sheath_Lambda - a_n*u/(2.d0*Te)))
  end function gcur

  !> The raw residual of whichever branch the caller wants to differentiate.
  real*8 function raw(u, Te, rho, vp, zj)
    real*8, intent(in) :: u, Te, rho, vp, zj
    real*8 :: r, a1,a2,a3,a4,a5,a6, rr, sc, xn, jr
    logical :: st
    call sheath_current_row(u, Te, rho, vp, zj, .false., r, a1,a2,a3,a4,a5,a6, rr, sc, st, xn, jr)
    raw = rr
  end function raw

  !> Central difference of the raw residual with respect to variable `iv`.
  real*8 function fd(iv, u, Te, rho, vp, zj)
    integer, intent(in) :: iv
    real*8,  intent(in) :: u, Te, rho, vp, zj
    real*8 :: h
    fd = 0.d0
    select case (iv)
    case (1) ; h = 1.d-7*abs(u)
               fd = (raw(u+h,Te,rho,vp,zj) - raw(u-h,Te,rho,vp,zj))/(2.d0*h)
    case (2) ; h = 1.d-7*Te
               fd = (raw(u,Te+h,rho,vp,zj) - raw(u,Te-h,rho,vp,zj))/(2.d0*h)
    case (3) ; h = 1.d-7*rho
               fd = (raw(u,Te,rho+h,vp,zj) - raw(u,Te,rho-h,vp,zj))/(2.d0*h)
    case (4) ; h = 1.d-7*abs(vp)
               fd = (raw(u,Te,rho,vp+h,zj) - raw(u,Te,rho,vp-h,zj))/(2.d0*h)
    case (5) ; h = 1.d-7*abs(zj)
               fd = (raw(u,Te,rho,vp,zj+h) - raw(u,Te,rho,vp,zj-h))/(2.d0*h)
    end select
  end function fd

  !> Both branches at X = +eps and X = -eps must give the same normalised row and residual.
  !> Both branches evaluated a hair either side of X = 0 must give the same NORMALISED row and
  !! residual. This is the continuity property derived in the module header; without it a node
  !! oscillating about saturation would see two unrelated equations.
  subroutine check_join(eps)
    real*8, intent(in) :: eps
    real*8 :: zp, zm
    real*8 :: r1,u1,t1,e1,d1,v1,j1,rr1,sc1,x1,q1
    real*8 :: r2,u2,t2,e2,d2,v2,j2,rr2,sc2,x2,q2
    logical :: s1, s2
    zp = jsat0 * (1.d0 - eps)        ! X = +eps, unsaturated
    zm = jsat0 * (1.d0 + eps)        ! X = -eps, saturated
    call sheath_current_row(u0,Te0,rho0,vpar0,zp,.false., r1,u1,t1,e1,d1,v1,j1, rr1,sc1,s1,x1,q1)
    call sheath_current_row(u0,Te0,rho0,vpar0,zm,.false., r2,u2,t2,e2,d2,v2,j2, rr2,sc2,s2,x2,q2)
    if ( s1 .or. (.not. s2) ) call fail_1('branch test misclassified X = +-eps', eps)
    ! --- A constraint row is defined up to an overall sign; align them before comparing.
    if ( j1*j2 .lt. 0.d0 ) then
      d2 = -d2 ; v2 = -v2 ; j2 = -j2 ; r2 = -r2 ; u2 = -u2 ; e2 = -e2
    endif
    call check('join rho ', d1, d2)
    call check('join vpar', v1, v2)
    call check('join zj  ', j1, j2)
    if ( abs(r1) + abs(r2) .gt. 1.d-6 ) call fail_2('branches do not join in residual', r1, r2)
    ! --- The u and Te columns must be vanishing on the unsaturated side, and exactly zero on the
    ! --- saturated side. They fade like eps, so the bound is set by eps, not by a fixed number.
    if ( abs(u1) + abs(e1) .gt. 1.d-3 ) call fail_2('u/Te columns do not fade at saturation', u1, e1)
    if ( abs(u2) + abs(e2) .ne. 0.d0  ) call fail_2('saturated row constrains u or Te', u2, e2)
  end subroutine check_join

  subroutine check(name, c, d)
    character(*), intent(in) :: name
    real*8,       intent(in) :: c, d
    if ( abs(c - d) .gt. 1.d-5 * max(abs(c), abs(d), 1.d-30) ) then
      sheath_current_selftest = .false.
      if (my_id .eq. 0) write(*,'(A,A,A,2es22.14)') &
        ' SHEATH_J SELFTEST FAIL: ', name, ' /= finite difference:', c, d
    endif
  end subroutine check

  subroutine fail_1(name, a)
    character(*), intent(in) :: name
    real*8,       intent(in) :: a
    sheath_current_selftest = .false.
    if (my_id .eq. 0) write(*,'(A,A,A,es22.14)') ' SHEATH_J SELFTEST FAIL: ', name, ':', a
  end subroutine fail_1

  subroutine fail_2(name, a, b)
    character(*), intent(in) :: name
    real*8,       intent(in) :: a, b
    sheath_current_selftest = .false.
    if (my_id .eq. 0) write(*,'(A,A,A,2es22.14)') ' SHEATH_J SELFTEST FAIL: ', name, ':', a, b
  end subroutine fail_2

end function sheath_current_selftest


end module mod_sheath_current
