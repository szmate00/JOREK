!> Wall current and sheath potential diagnostic for the charge-conserving sheath boundary
!! condition (bcs(:)%natural%u, model600).
!!
!! The quantities are accumulated at the boundary Gauss points inside
!! mod_boundary_matrix_open.f90, where the geometry and the plasma state are available anyway,
!! and reduced and printed once per matrix construction from construct_matrix_mod.f90. Doing it
!! that way avoids duplicating the boundary element identification, at the price of a critical
!! section per boundary Gauss point, which is negligible: boundary elements are a small fraction
!! of the mesh.
!!
!! Printed per time step:
!!   I_sheath  total current through the wall using the current the sheath actually passes [A]
!!   I_Ampere  the same integral using the interior current zj = Delta*psi [A]
!! The two agree once the potential has adjusted, so their difference is a direct convergence
!! measure of the boundary condition, and I_sheath tells you whether a floating wall model would
!! change anything (it should be small compared to the current flowing in each direction).
module mod_sheath_diag

  use phys_module, only: max_bnd_types

  implicit none

  private

  public :: sheath_diag_reset, sheath_diag_add, sheath_diag_report
  public :: sheath_diag_add_nodal, sheath_diag_add_weak
  public :: sheath_store_psi0, sheath_psi0
  public :: sheath_init_potential

  !> psi's degrees of freedom at t_start, per node. The wall relaxation needs the DEVIATION of
  !! dpsi/dn from its equilibrium value; driving it with the raw value would make the flux drift
  !! even in a quiet plasma. Indexed (dof, node).
  real*8, allocatable, save :: sheath_psi0(:,:)

  real*8,  save :: sd_I_sheath(max_bnd_types) = 0.d0  !< current through the wall, sheath value [A]
  real*8,  save :: sd_I_amp(max_bnd_types)    = 0.d0  !< the same with the interior Ampere current [A]
  real*8,  save :: sd_area(max_bnd_types)     = 0.d0  !< wetted area [m^2], for averages
  ! --- Geometry of each boundary type, to separate two failure modes that no magnitude gate can
  ! --- tell apart. b_n SIGN REVERSAL within one type flips zj_sat mid-face, so the constraint
  ! --- drives zj opposite ways on the two halves and the integrated I_sheath partially cancels.
  ! --- Cell size shows whether a type is resolved on a very different scale from the others -
  ! --- zj = Delta*psi amplifies an interpolated psi by 1/h^2.
  real*8,  save :: sd_bn_min(max_bnd_types)   =  1.d30 !< min SIGNED B.n
  real*8,  save :: sd_bn_max(max_bnd_types)   = -1.d30 !< max SIGNED b_n
  real*8,  save :: sd_area_bn_neg(max_bnd_types) = 0.d0 !< area with b_n < 0
  real*8,  save :: sd_I_bn_neg(max_bnd_types)    = 0.d0 !< I_sheath over that area
  ! --- |zj_sat| per type. Nothing else distinguishes the two ways max|j/jsat| can explode:
  ! --- zj0 rising (the plasma delivering more current) or zj_sat collapsing (the sheath losing
  ! --- its capacity). They need opposite fixes and the ratio alone cannot tell them apart.
  real*8,  save :: sd_sat_min(max_bnd_types)  =  1.d30 !< min |zj_sat| on this type
  real*8,  save :: sd_sat_sum(max_bnd_types)  =  0.d0  !< area-weighted sum of |zj_sat|
  real*8,  save :: sd_ds_min(max_bnd_types)   =  1.d30 !< smallest area element on this type
  real*8,  save :: sd_ds_max(max_bnd_types)   = -1.d30
  real*8,  save :: sd_phi_min                 =  1.d30
  real*8,  save :: sd_phi_max                 = -1.d30
  real*8,  save :: sd_phi_sum                 = 0.d0  !< area weighted, for the mean
  !> Same, but only where the obliqueness gate leaves the sheath term ACTIVE. Where the gate has
  !! removed it, u has no boundary condition at all (dirichlet%u is .false. on these types), so a
  !! runaway there means an unconstrained null space rather than a sheath the plasma is overdriving.
  !! Comparing the two maxima separates those two completely different failures.
  real*8,  save :: sd_phi_max_act             = -1.d30
  real*8,  save :: sd_gate_off_area           = 0.d0
  real*8,  save :: sd_ratio_max               = 0.d0  !< max |j/j_sat| demanded by the interior
  real*8,  save :: sd_lim_area                = 0.d0  !< area sitting on the electron side limiter

  !> Effective strength the NODAL sheath_zj constraint is actually applied with, accumulated in
  !! mod_boundary_conditions. szj_rel is the product of sheath_zj_relax, the ramp and the two
  !! gates; where it is small the row degenerates towards d(zj) = 0, i.e. the Dirichlet freeze the
  !! constraint replaces, and zj is then free to sit far off the sheath characteristic. That is
  !! invisible in every other number here - `gated-off area` only sees the obliqueness gate - and
  !! it is the natural explanation for I_Ampere and I_sheath differing by a constant factor.
  !! Zeroed in sheath_diag_report AFTER packing, not in sheath_diag_reset: boundary_conditions is
  !! called after the report within one matrix construction, so these describe the PREVIOUS pass.
  real*8,  save :: sd_zjrel_sum               = 0.d0  !< sum of szj_rel over sheath_zj nodes
  real*8,  save :: sd_zjrel_min               = 1.d30
  real*8,  save :: sd_zjrel_n                 = 0.d0  !< node count, as a real for one reduction
  real*8,  save :: sd_zjratio_sum             = 0.d0  !< sum of |zj/zj_sat| demanded, same nodes

  !> Inner/outer target split, keyed on major radius against sheath_diag_R_split. Boundary type 1
  !! carries BOTH targets, so every total above averages the two together - and an in-out
  !! difference is exactly what HFSHD is. Opposite errors on the two targets cancel in the sum.
  !! Index 1 = inner (R < split), 2 = outer. Second index: 1 I_sheath, 2 I_Ampere, 3 area, 4 phi*dS
  real*8,  save :: sd_io_sum(2,4)             = 0.d0
  real*8,  save :: sd_io_max(2)               = 0.d0  !< max |j/j_sat| on each target

  !> Where the worst NODE is. A single node whose |zj/zj_sat| grows without bound is the observed
  !! failure on the target plate: the validity gate shuts it off progressively (strength min decays
  !! monotonically) and it diverges anyway. Knowing whether it sits at a corner, at the strike
  !! point, or on a leg edge decides what to fix, and no other number here localises it.
  real*8,  save :: sd_worst_ratio             = -1.d0
  real*8,  save :: sd_worst_RZ(2)             = 0.d0

  !> Acceptance metric for the constraint itself: how far zj is from the characteristic, as a
  !! fraction of the current the sheath can actually pass. Reported separately at NODES (where the
  !! constraint is imposed) and at GAUSS points (where the currents are integrated), because the
  !! two failing differently means different things - node small / Gauss large is an interpolation
  !! or derivative-DOF problem, both large means the row is not being satisfied at all.
  !!     eps = sqrt(sum (zj - zj_sh)^2 / sum zj_sat^2)
  real*8,  save :: sd_res_node(2)             = 0.d0  !< 1: sum (zj-zj_sh)^2   2: sum zj_sat^2
  real*8,  save :: sd_res_gauss(2)            = 0.d0  !< same, weighted by dS

  !> WEAK (Galerkin) residual of the sheath characteristic on the boundary trace space:
  !!     F_a = integral over Gamma of  N_a * (zj - zj_sh) dS,   normalised by the trace mass
  !!     diagonal D_a = integral of N_a*N_a dS.
  !! This is what a projection would drive to zero, and it is NOT the pointwise residual: a
  !! Galerkin condition only needs (zj - zj_sh) ORTHOGONAL to the trace space, not zero. Reporting
  !! both separates "the trace cannot represent zj_sh" from "the trace is not being projected onto
  !! it", which the pointwise number alone cannot do. Purely diagnostic - nothing is imposed.
  real*8,  save :: sd_wk_res2                 = 0.d0  !< sum of (F_a/D_a)^2
  real*8,  save :: sd_wk_ref2                 = 0.d0  !< sum of (S_a/D_a)^2, S_a = int N_a zj_sat
  real*8,  save :: sd_wk_n                    = 0.d0  !< number of trace rows sampled
  real*8,  save :: sd_wk_dmin                 =  1.d30
  real*8,  save :: sd_wk_dmax                 = -1.d30

contains

!> Zero the accumulators. Called once per matrix construction, before the element loop.
subroutine sheath_diag_reset()
  implicit none
  sd_I_sheath = 0.d0
  sd_I_amp    = 0.d0
  sd_area     = 0.d0
  sd_bn_min      =  1.d30 ; sd_bn_max      = -1.d30
  sd_area_bn_neg =  0.d0  ; sd_I_bn_neg    =  0.d0
  sd_ds_min      =  1.d30 ; sd_ds_max      = -1.d30
  sd_sat_min     =  1.d30 ; sd_sat_sum     =  0.d0
  sd_phi_min  =  1.d30
  sd_phi_max  = -1.d30
  sd_phi_sum  = 0.d0
  sd_phi_max_act = -1.d30
  sd_gate_off_area = 0.d0
  sd_ratio_max= 0.d0
  sd_io_sum   = 0.d0
  sd_io_max   = 0.d0
  sd_worst_ratio = -1.d0
  sd_worst_RZ    = 0.d0
  sd_res_gauss   = 0.d0
  sd_wk_res2 = 0.d0; sd_wk_ref2 = 0.d0; sd_wk_n = 0.d0
  sd_wk_dmin = 1.d30; sd_wk_dmax = -1.d30

  sd_lim_area = 0.d0
end subroutine sheath_diag_reset


!> Add the contribution of one boundary Gauss point on one toroidal plane.
!!
!! @param bnd_type  boundary type of the node
!! @param zj_sh     sheath current (zj units)
!! @param zj0       interior current at the same point (zj units)
!! @param zj_sat    ion saturation current (zj units)
!! @param x_lim     sheath exponent actually used, after limiting
!! @param u0        potential variable
!! @param Te0       electron temperature (JOREK units)
!! @param Bdotn     B.n
!! @param dS        surface element of this sample: ws*dl*R*(2*pi/n_plane) [m^2]
subroutine sheath_diag_add(bnd_type, zj_sh, zj0, zj_sat, x_lim, u0, Te0, Bdotn, dS, gate, &
                           V_wall_loc, R_loc)

  use constants,     only: MU_ZERO
  use phys_module,   only: F0, sheath_X_min, sheath_smooth_dX, sheath_diag_R_split
  use mod_sheath_bc, only: sheath_norm

  implicit none
  integer, intent(in) :: bnd_type
  real*8,  intent(in) :: zj_sh, zj0, zj_sat, x_lim, u0, Te0, Bdotn, dS
  real*8,  intent(in), optional :: gate   !< obliqueness weight; the ratio below is only
                                          !< meaningful where the term is actually active
  real*8,  intent(in), optional :: V_wall_loc !< local wall bias, so e*Phi/kTe is measured against
                                          !< the SAME reference the constraint imposes. Omitting it
                                          !< silently reports against the global sheath_V_wall,
                                          !< which differs wherever sheath_V_wall_asym /= 0.
  real*8,  intent(in), optional :: R_loc   !< major radius, for the inner/outer target split

  real*8 :: jn_sheath, jn_amp, phi_over_te, ratio, a_n, c_sat, vw
  integer :: io

  if ( bnd_type .lt. 1 .or. bnd_type .gt. max_bnd_types ) return

  ! --- J.n = J_par*(b.n) = zj*(B.n)/(F0*mu0)   [A/m^2]
  jn_sheath = zj_sh * Bdotn / (F0 * MU_ZERO)
  jn_amp    = zj0   * Bdotn / (F0 * MU_ZERO)

  if ( present(V_wall_loc) ) then
    call sheath_norm(a_n, c_sat, vw, V_wall_loc)
  else
    call sheath_norm(a_n, c_sat, vw)
  endif
  phi_over_te = ( 0.5d0*a_n*u0 - vw ) / max(Te0, 1.d-14)     ! e*Phi/(k*Te)

  ! --- Solvability ratio. The characteristic can only deliver f in (-(exp(-X_min)-1), 1], so
  ! --- abs(zj0/zj_sat) > 1 means NO u satisfies it at this point and u is driven without bound.
  ! --- Only report it where the obliqueness gate leaves the term active: at a gated-off point
  ! --- the ratio diverges harmlessly because the term is not there.
  ratio = 0.d0
  if ( abs(zj_sat) .gt. 0.d0 ) ratio = abs(zj0 / zj_sat)
  if ( present(gate) ) then
    if ( gate .lt. 0.5d0 ) ratio = 0.d0
  endif

  !$omp critical (sheath_diag_accumulate)
  sd_I_sheath(bnd_type)  = sd_I_sheath(bnd_type) + jn_sheath * dS
  sd_I_amp(bnd_type)     = sd_I_amp(bnd_type)    + jn_amp    * dS
  sd_area(bnd_type)      = sd_area(bnd_type)     + dS

  ! --- Bdotn is the signed, dimensional B.n. Its sign is the discriminator here; the magnitude
  ! --- printed below is in JOREK's magnetic-field units, not the dimensionless b_n=B.n/|B|.
  sd_bn_min(bnd_type) = min(sd_bn_min(bnd_type), Bdotn)
  sd_bn_max(bnd_type) = max(sd_bn_max(bnd_type), Bdotn)
  sd_sat_min(bnd_type) = min(sd_sat_min(bnd_type), abs(zj_sat))
  sd_sat_sum(bnd_type) = sd_sat_sum(bnd_type) + abs(zj_sat) * dS
  sd_ds_min(bnd_type) = min(sd_ds_min(bnd_type), dS)
  sd_ds_max(bnd_type) = max(sd_ds_max(bnd_type), dS)
  if ( Bdotn .lt. 0.d0 ) then
    sd_area_bn_neg(bnd_type) = sd_area_bn_neg(bnd_type) + dS
    sd_I_bn_neg(bnd_type)    = sd_I_bn_neg(bnd_type)    + jn_sheath * dS
  endif
  sd_phi_sum             = sd_phi_sum            + phi_over_te * dS
  sd_phi_min             = min(sd_phi_min, phi_over_te)
  sd_phi_max             = max(sd_phi_max, phi_over_te)
  if ( present(gate) ) then
    if ( gate .ge. 0.5d0 ) then
      sd_phi_max_act = max(sd_phi_max_act, phi_over_te)
    else
      sd_gate_off_area = sd_gate_off_area + dS
    endif
  else
    sd_phi_max_act = max(sd_phi_max_act, phi_over_te)
  endif
  sd_ratio_max           = max(sd_ratio_max, ratio)
  sd_res_gauss(1)        = sd_res_gauss(1) + (zj0 - zj_sh)**2 * dS
  sd_res_gauss(2)        = sd_res_gauss(2) + zj_sat**2        * dS
  ! --- inner/outer split; ratio here is already gate-masked, so it reports where the term is live
  if ( present(R_loc) .and. (sheath_diag_R_split .gt. 0.d0) ) then
    io = 1
    if ( R_loc .ge. sheath_diag_R_split ) io = 2
    sd_io_sum(io,1) = sd_io_sum(io,1) + jn_sheath * dS
    sd_io_sum(io,2) = sd_io_sum(io,2) + jn_amp    * dS
    sd_io_sum(io,3) = sd_io_sum(io,3) + dS
    sd_io_sum(io,4) = sd_io_sum(io,4) + phi_over_te * dS
    sd_io_max(io)   = max(sd_io_max(io), ratio)
  endif
  ! --- area where the electron side limiter is biting, i.e. where the wall is close to
  ! --- electron saturation and the characteristic is being held back
  if ( x_lim .lt. sheath_X_min + 2.d0*max(sheath_smooth_dX,1.d-3) ) sd_lim_area = sd_lim_area + dS
  !$omp end critical (sheath_diag_accumulate)

end subroutine sheath_diag_add


!> Record the strength the nodal sheath_zj constraint is applied with at one boundary node.
!!
!! Called from mod_boundary_conditions once szj_rel is final, i.e. after the relaxation, the ramp,
!! the obliqueness gate and the validity gate have all been folded in.
!!
!! @param szj_rel   effective row strength, 1 = constraint imposed fully, 0 = zj left frozen
!! @param zj_ratio  |zj/zj_sat| the interior is demanding at this node
subroutine sheath_diag_add_nodal(szj_rel, zj_ratio, R_node, Z_node, zj_resid, zj_sat)

  implicit none
  real*8, intent(in) :: szj_rel, zj_ratio
  real*8, intent(in), optional :: R_node, Z_node   !< position, to locate the worst node
  real*8, intent(in), optional :: zj_resid         !< zj_sheath - zj at this node
  real*8, intent(in), optional :: zj_sat           !< ion saturation current, to normalise it

  sd_zjrel_sum   = sd_zjrel_sum   + szj_rel
  sd_zjrel_min   = min(sd_zjrel_min, szj_rel)
  sd_zjrel_n     = sd_zjrel_n     + 1.d0
  sd_zjratio_sum = sd_zjratio_sum + zj_ratio

  if ( present(zj_resid) ) sd_res_node(1) = sd_res_node(1) + zj_resid**2
  if ( present(zj_sat)   ) sd_res_node(2) = sd_res_node(2) + zj_sat**2

  if ( zj_ratio .gt. sd_worst_ratio ) then
    sd_worst_ratio = zj_ratio
    if ( present(R_node) ) sd_worst_RZ(1) = R_node
    if ( present(Z_node) ) sd_worst_RZ(2) = Z_node
  endif

end subroutine sheath_diag_add_nodal


!> One boundary trace row's weak sheath residual, already divided by its trace mass diagonal.
!!
!! @param res_over_D  F_a / D_a, in zj units
!! @param ref_over_D  S_a / D_a, the same projection applied to zj_sat, to normalise with
!! @param diag        D_a itself, to expose the value/derivative row scaling spread
subroutine sheath_diag_add_weak(res_over_D, ref_over_D, diag)

  implicit none
  real*8, intent(in) :: res_over_D, ref_over_D, diag

  sd_wk_res2 = sd_wk_res2 + res_over_D**2
  sd_wk_ref2 = sd_wk_ref2 + ref_over_D**2
  sd_wk_n    = sd_wk_n    + 1.d0
  sd_wk_dmin = min(sd_wk_dmin, diag)
  sd_wk_dmax = max(sd_wk_dmax, diag)

end subroutine sheath_diag_add_weak


!> Reduce over MPI ranks and print one line. Called after the element loop; rank 0 prints.
!!
!! Every rank calls all three collectives unconditionally (no early return), and only arrays are
!! passed, so this is safe with any MPI interface and cannot deadlock.
subroutine sheath_diag_report(my_id)

  use mpi_mod

  implicit none
  integer, intent(in) :: my_id

  integer, parameter :: ns = 3*max_bnd_types + 2
  integer, parameter :: nx = ns + 20                 ! +1 gate-off, +2..5 nodal, +6..13 in/out, +14..17 resid, +18..20 weak
  real*8  :: loc_sum(nx), glo_sum(nx), loc_max(6), glo_max(6), loc_min(3), glo_min(3)
  real*8  :: area_tot, I_sh_tot, I_am_tot, phi_mean, lim_frac, zjrel_mean, zjratio_mean
  real*8  :: wr_loc(1), wr_glo(1), wr_own(3), wr_all(3), eps_node, eps_gauss, wk_eps
  integer :: io, ib0
  character(len=5), parameter :: io_name(2) = (/ 'INNER', 'OUTER' /)
  integer :: ierr, i, i0, i1, i2
  real*8  :: gbn_min(max_bnd_types), gbn_max(max_bnd_types)
  real*8  :: gds_min(max_bnd_types), gds_max(max_bnd_types)
  real*8  :: gan_neg(max_bnd_types), gin_neg(max_bnd_types)
  real*8  :: gst_min(max_bnd_types), gst_sum(max_bnd_types)

  i0 = 0                      ! offsets into the packed reduction buffer
  i1 =   max_bnd_types
  i2 = 2*max_bnd_types

  loc_sum(i0+1:i0+max_bnd_types) = sd_I_sheath
  loc_sum(i1+1:i1+max_bnd_types) = sd_I_amp
  loc_sum(i2+1:i2+max_bnd_types) = sd_area
  loc_sum(ns-1)                  = sd_phi_sum
  loc_sum(ns)                    = sd_lim_area

  loc_max(1) = sd_phi_max
  loc_max(2) = sd_ratio_max
  loc_max(3) = sd_phi_max_act
  loc_sum(ns+1) = sd_gate_off_area
  loc_sum(ns+2) = sd_zjrel_sum
  loc_sum(ns+3) = sd_zjrel_n
  loc_sum(ns+4) = sd_zjratio_sum
  loc_sum(ns+5) = 0.d0
  loc_sum(ns+6 :ns+9 ) = sd_io_sum(1,:)
  loc_sum(ns+10:ns+13) = sd_io_sum(2,:)
  loc_sum(ns+14:ns+15) = sd_res_node
  loc_sum(ns+16:ns+17) = sd_res_gauss
  loc_sum(ns+18) = sd_wk_res2
  loc_sum(ns+19) = sd_wk_ref2
  loc_sum(ns+20) = sd_wk_n
  loc_max(6) = sd_wk_dmax
  ! --- loc_min goes through MPI_MIN, so the raw value is already the global minimum. Negating it
  ! --- (as a "MIN via MAX" trick) was wrong twice over: the reduction is MIN, not MAX, and a rank
  ! --- with no sheath boundary carries the 1.d30 sentinel whose negation then WINS the MIN,
  ! --- reporting 1.d30. Ranks with no data lose a plain MIN, which is the behaviour wanted.
  loc_min(3) = sd_wk_dmin
  loc_max(4) = sd_io_max(1)
  loc_max(5) = sd_io_max(2)
  loc_min(1) = sd_phi_min
  loc_min(2) = sd_zjrel_min

  call MPI_Reduce(loc_sum, glo_sum, nx, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(loc_max, glo_max,  6, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(loc_min, glo_min,  3, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)

  ! --- Locate the worst node. Take the global max, then have ONLY the owning rank contribute its
  ! --- position to a SUM reduction, so the position arrives on rank 0 and the line is printed
  ! --- there with everything else. Printing from the owning rank directly would be simpler but
  ! --- many MPI launchers only capture rank 0's stdout, which would silently lose the line.
  wr_loc(1) = sd_worst_ratio
  call MPI_Allreduce(wr_loc, wr_glo, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
  wr_own = 0.d0
  if ( (wr_glo(1) .gt. 0.d0) .and. (sd_worst_ratio .ge. wr_glo(1)) ) then
    wr_own(1) = sd_worst_RZ(1)
    wr_own(2) = sd_worst_RZ(2)
    wr_own(3) = 1.d0                      ! owner count, so a tie can be divided out
  endif
  call MPI_Reduce(wr_own, wr_all, 3, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  ! --- zero the nodal block here, not in sheath_diag_reset: boundary_conditions runs AFTER this
  ! --- report within one matrix construction, so zeroing at reset would wipe the pass we want.
  sd_zjrel_sum = 0.d0; sd_zjrel_n = 0.d0; sd_zjratio_sum = 0.d0; sd_zjrel_min = 1.d30
  sd_worst_ratio = -1.d0; sd_worst_RZ = 0.d0; sd_res_node = 0.d0

  ! --- These are COLLECTIVE and must be called on every rank, so they belong above the
  ! --- my_id/=0 return. Placing them below it would have rank 0 enter a reduction alone and
  ! --- hang every other rank at the next collective.
  call MPI_Reduce(sd_bn_min,      gbn_min, max_bnd_types, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sd_bn_max,      gbn_max, max_bnd_types, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sd_ds_min,      gds_min, max_bnd_types, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sd_ds_max,      gds_max, max_bnd_types, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sd_area_bn_neg, gan_neg, max_bnd_types, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sd_I_bn_neg,    gin_neg, max_bnd_types, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sd_sat_min,     gst_min, max_bnd_types, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sd_sat_sum,     gst_sum, max_bnd_types, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if ( my_id .ne. 0 ) return

  area_tot = sum( glo_sum(i2+1:i2+max_bnd_types) )
  if ( area_tot .le. 0.d0 ) return                     ! the sheath BC is not active anywhere

  I_sh_tot = sum( glo_sum(i0+1:i0+max_bnd_types) )
  I_am_tot = sum( glo_sum(i1+1:i1+max_bnd_types) )
  phi_mean = glo_sum(ns-1) / area_tot
  lim_frac = 1.d2 * glo_sum(ns) / area_tot

  write(*,'(A,es11.3,A,es11.3,A,f7.2,A,f7.2,A,f7.2,A,es9.2,A,f5.1,A)')             &
    ' SHEATH: I_wall=', I_sh_tot, ' A (Ampere ', I_am_tot,                         &
    ' A)  ePhi/kTe min/mean/max=', glo_min(1), ' /', phi_mean, ' /', glo_max(1),   &
    '  max|j/jsat|=', glo_max(2), '  e-limited ', lim_frac, ' %'

  ! --- The two maxima separate the two failures: if the ACTIVE max settles while the global one
  ! --- runs away, the runaway is at gated-off points where u has no boundary condition at all.
  write(*,'(A,f8.2,A,f5.1,A)')                                                     &
    '         ePhi/kTe max where the sheath is ACTIVE=', glo_max(3),               &
    '   gated-off area ', 1.d2*glo_sum(ns+1)/max(area_tot,1.d-30), ' %'

  if ( (wr_glo(1) .gt. 0.d0) .and. (wr_all(3) .gt. 0.d0) )                          &
    write(*,'(A,es10.3,A,f7.4,A,f8.4,A)')                                           &
      '         worst node |zj/zj_sat|=', wr_glo(1),                                &
      ' at R=', wr_all(1)/wr_all(3), ' Z=', wr_all(2)/wr_all(3), ' m'

  ! --- Acceptance metric: is zj actually ON the characteristic? eps = ||zj - zj_sh|| / ||zj_sat||,
  ! --- at nodes (where the row is imposed) and at Gauss points (where the currents are
  ! --- integrated). Small at nodes but large at Gauss points would mean the constraint holds
  ! --- pointwise and the interpolation between nodes does not; both large means the row itself is
  ! --- not being satisfied, which is what a flat characteristic (small df/dX) produces.
  eps_node  = -1.d0
  eps_gauss = -1.d0
  if ( glo_sum(ns+15) .gt. 0.d0 ) eps_node  = sqrt( glo_sum(ns+14) / glo_sum(ns+15) )
  if ( glo_sum(ns+17) .gt. 0.d0 ) eps_gauss = sqrt( glo_sum(ns+16) / glo_sum(ns+17) )
  if ( (eps_node .ge. 0.d0) .or. (eps_gauss .ge. 0.d0) )                             &
    write(*,'(A,es10.3,A,es10.3)')                                                   &
      '         |zj-zj_sh|/|zj_sat|  nodes=', eps_node, '  gauss=', eps_gauss

  ! --- WEAK residual: what a Galerkin projection onto the boundary trace space would drive to
  ! --- zero. A projection only needs (zj - zj_sh) ORTHOGONAL to that space, so this can be far
  ! --- smaller than the pointwise `gauss` figure above without anything being wrong. If it is
  ! --- already small, a weak formulation has little left to gain and the closure gap lies
  ! --- elsewhere. D_max/D_min is the value-row vs derivative-row scaling spread, which is what
  ! --- makes a Galerkin row block hard to condition.
  if ( glo_sum(ns+20) .gt. 0.d0 ) then
    wk_eps = 0.d0
    if ( glo_sum(ns+19) .gt. 0.d0 ) wk_eps = sqrt( glo_sum(ns+18) / glo_sum(ns+19) )
    write(*,'(A,es10.3,A,i0,A,es9.2,A,es9.2)')                                       &
      '         weak |F_a/D_a|/|S_a/D_a|=', wk_eps, '  over ', nint(glo_sum(ns+20)), &
      ' trace samples;  D min=', glo_min(3), ' max=', glo_max(6)
  endif

  ! --- How hard the nodal constraint is actually being pushed. mean(szj_rel) well below 1 means
  ! --- the row is mostly the Dirichlet freeze, so zj is under no obligation to match zj_sheath and
  ! --- I_Ampere/I_sheath cannot close. Compare mean|j/jsat| with max|j/jsat| above: a mean near 1
  ! --- with a large max is a hot spot, a large mean is a wall-wide violation.
  if ( glo_sum(ns+3) .gt. 0.d0 ) then
    zjrel_mean   = glo_sum(ns+2) / glo_sum(ns+3)
    zjratio_mean = glo_sum(ns+4) / glo_sum(ns+3)
    write(*,'(A,f7.4,A,f7.4,A,f8.2,A,i0,A)')                                       &
      '         sheath_zj strength mean=', zjrel_mean, ' min=', glo_min(2),        &
      '   mean|j/jsat|=', zjratio_mean, '  over ', nint(glo_sum(ns+3)), ' nodes'
  endif

  ! --- Inner vs outer target. HFSHD IS an in-out difference, so the totals above - which average
  ! --- the two targets together and let opposite errors cancel - cannot answer the question.
  do io = 1, 2
    ib0 = ns + 2 + 4*io          ! inner block starts at ns+6, outer at ns+10
    if ( glo_sum(ib0+2) .le. 0.d0 ) cycle                      ! +0 I_sheath  +1 I_Ampere
    write(*,'(A,A,A,es11.3,A,es11.3,A,f7.2,A,es9.2,A,es10.3,A)')                   &
      '         ', io_name(io), ' target: I_sheath=', glo_sum(ib0),                &
      ' A  I_Ampere=', glo_sum(ib0+1),                                             &
      '  ePhi/kTe=', glo_sum(ib0+3) / glo_sum(ib0+2),                              &
      '  max|j/jsat|=', glo_max(3+io), '  area=', glo_sum(ib0+2), ' m^2'
  enddo

  do i = 1, max_bnd_types
    if ( glo_sum(i2+i) .le. 0.d0 ) cycle
    write(*,'(A,i3,A,es11.3,A,es11.3,A,es10.3,A)')                                &
      '         bnd type', i, ': I_sheath=', glo_sum(i0+i), ' A  I_Ampere=',      &
      glo_sum(i1+i), ' A  area=', glo_sum(i2+i), ' m^2'
  enddo

  ! --- Geometry per boundary type. Two things this separates that nothing else does:
  ! ---  * B.n range straddling zero => the field REVERSES through the face, so zj_sat flips sign
  ! ---    mid-type, the constraint drives zj opposite ways on the two halves and the integrated
  ! ---    I_sheath partially cancels. A magnitude gate (sheath_min_bn) cannot see this.
  ! ---  * dS spread => how finely this type is resolved. zj = Delta*psi amplifies an interpolated
  ! ---    psi by 1/h^2, so a type packed into a short chord is far more sensitive than a type
  ! ---    spread over the aligned grid.

  do i = 1, max_bnd_types
      if ( glo_sum(i2+i) .le. 0.d0 ) cycle
      write(*,'(A,i3,A,es10.2,A,es10.2,A,f6.1,A,es9.2,A,es9.2,A,es10.2,A)')        &
        '         bnd type', i, ' geom: B.n ', gbn_min(i), ' ..', gbn_max(i),      &
        '   area(b_n<0) ', 1.d2*gan_neg(i)/max(glo_sum(i2+i),1.d-30), ' %   dS ', &
        gds_min(i), ' ..', gds_max(i),                                             &
        '   I(b_n<0)=', gin_neg(i), ' A'
    ! --- 6 descriptors, 6 items. |zj_sat| min vs mean says WHICH way max|j/jsat| blew up:
    ! --- a collapsing min with a steady mean is a local loss of sheath capacity; a steady min
    ! --- means the numerator zj0 rose and the plasma is delivering the current.
    write(*,'(A,i3,A,es10.2,A,es10.2)')                                            &
      '         bnd type', i, ' |zj_sat| min ', gst_min(i),                        &
      '  area-mean ', gst_sum(i)/max(glo_sum(i2+i),1.d-30)
  enddo

end subroutine sheath_diag_report


!> Initialise u on the natural%u boundary types to the floating potential, Lambda*Te/e.
!!
!! Without this the run starts at u = 0, i.e. Phi = 0, which is electron saturation: the sheath
!! immediately demands the full electron saturation current everywhere and the boundary condition
!! has to travel ~Lambda*Te away from where the plasma is, through an exponential, in one implicit
!! step. Starting at X = 0 puts the state where the linearisation is valid.
!!
!! Only the axisymmetric component is touched. Lambda is treated as locally constant so that the
!! derivative degrees of freedom can carry the same relation, keeping grad(u) consistent with
!! grad(Te) rather than leaving it at the old field's value.
subroutine sheath_init_potential(node_list, my_id)

  use mod_parameters
  use data_structure
  use phys_module,   only: bcs, sheath_init_u_all
  use mod_sheath_bc, only: sheath_norm, sheath_V_wall_at, sheath_get_lambda
  use mpi_mod

  implicit none

  type (type_node_list), intent(inout) :: node_list
  integer,               intent(in)    :: my_id

  integer :: i, id, ib, ierr
  real*8  :: a_n, c_sat, vw, Ti0, Te0, T0, lam, dlam_dTi, dlam_dTe, cfac
  real*8  :: n_loc(2), n_glo(2), u_loc(2), u_glo(2), u_wall
  logical :: is_sheath_type

  ! --- Pass 1: the mean floating potential over the sheath boundary. The GAUGE types are then set
  ! --- to that single CONSTANT rather than to their own local Lambda*Te/e.
  ! ---
  ! --- Why a constant: a metal wall is an EQUIPOTENTIAL. Te varies by orders of magnitude along
  ! --- it (0.66 eV at a cold strike point against tens of eV in the outer SOL), so most of the
  ! --- wall is NOT at its local floating potential - net current flows and closes through the
  ! --- wall, which is the physics the sheath BC exists to capture. Imposing Lambda*Te/e pointwise
  ! --- is neither an equipotential nor a self-consistent sheath: it makes Phi vary as strongly as
  ! --- Te, and since Delta*u ~ 0 in the interior u is harmonic and immediately smooths that
  ! --- variation, leaving cold nodes at a Phi set by their hot neighbours. Measured: ePhi/kTe
  ! --- starts at exactly 3.00 everywhere, then reaches 23.69 after ONE step - about the Te
  ! --- contrast between adjacent wall regions - and that excursion is completely insensitive to
  ! --- the sheath row (20x weaker constraint gave 23.69 to three digits), so it is the imposed
  ! --- profile, not the boundary condition.
  ! ---
  ! --- The derivative DOFs are left at ZERO on the gauge types for the same reason: an
  ! --- equipotential has no tangential gradient, and grad(u) along the wall IS ExB flow through it.
  u_loc = 0.d0

  if ( sheath_init_u_all ) then
    do i = 1, node_list%n_nodes
      ib = node_list%node(i)%boundary
      if ( ib .lt. 1 .or. ib .gt. max_bnd_types ) cycle
      if ( .not. (bcs(ib)%natural%u .or. bcs(ib)%sheath_zj .or. bcs(ib)%sheath_zj_weak) ) cycle
      if ( with_TiTe ) then
        Ti0 = node_list%node(i)%values(1,1,var_Ti)
        Te0 = node_list%node(i)%values(1,1,var_Te)
      else
        T0  = node_list%node(i)%values(1,1,var_T)
        Ti0 = 0.5d0 * T0
        Te0 = 0.5d0 * T0
      endif
      if ( Te0 .le. 0.d0 ) cycle
      call sheath_norm(a_n, c_sat, vw, sheath_V_wall_at(node_list%node(i)%x(1,1,1)))
      call sheath_get_lambda(Ti0, Te0, lam, dlam_dTi, dlam_dTe)
      u_loc(1) = u_loc(1) + 2.d0 * ( lam * Te0 + vw ) / a_n
      u_loc(2) = u_loc(2) + 1.d0
    enddo
    call MPI_Allreduce(u_loc, u_glo, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    u_wall = 0.d0
    if ( u_glo(2) .gt. 0.d0 ) u_wall = u_glo(1) / u_glo(2)
  endif

  n_loc = 0.d0

  do i = 1, node_list%n_nodes

    ib = node_list%node(i)%boundary
    if ( ib .lt. 1 .or. ib .gt. max_bnd_types ) cycle
    ! --- both sheath routes want u to start at its own fixed point rather than at 0, which is
    ! --- deep electron saturation (X = -Lambda) and demands ~19*j_sat of electron current
    !
    ! --- With sheath_init_u_all the Dirichlet GAUGE types are set as well. They are the same
    ! --- physical wall: freezing them at u = 0 asserts Phi = 0 there while an adjacent sheath
    ! --- boundary floats at Lambda*Te/e, and the difference appears as a node-to-node potential
    ! --- step whose along-wall gradient is ExB flow through the boundary. For pure ExB the
    ! --- offset would be arbitrary, but the sheath characteristic depends on ABSOLUTE Phi
    ! --- (X = e*Phi/kTe - Lambda), so once a sheath BC exists the gauge is fixed by physics.
    ! --- Dirichlet still freezes these nodes - it freezes them at a physical value instead of 0.
    is_sheath_type = bcs(ib)%natural%u .or. bcs(ib)%sheath_zj .or. bcs(ib)%sheath_zj_weak
    if ( .not. is_sheath_type ) then
      if ( .not. sheath_init_u_all ) cycle
    endif

    if ( with_TiTe ) then
      Ti0 = node_list%node(i)%values(1,1,var_Ti)
      Te0 = node_list%node(i)%values(1,1,var_Te)
    else
      T0  = node_list%node(i)%values(1,1,var_T)
      Ti0 = 0.5d0 * T0
      Te0 = 0.5d0 * T0
    endif
    if ( Te0 .le. 0.d0 ) cycle

    ! --- per node: with a differentially biased wall the floating potential is a function of R
    call sheath_norm(a_n, c_sat, vw, sheath_V_wall_at(node_list%node(i)%x(1,1,1)))

    call sheath_get_lambda(Ti0, Te0, lam, dlam_dTi, dlam_dTe)

    if ( is_sheath_type ) then
      ! --- zero net current (X = 0) sits at e*Phi/(k*Te) = Lambda, i.e. a_n*u/2 - vw = Lambda*Te.
      ! --- u is free on these types, so this is only a starting guess and the local value is the
      ! --- right one: it is where the sheath itself carries no net current.
      cfac = 2.d0 * lam / a_n
      node_list%node(i)%values(1,1,var_u) = cfac * Te0 + 2.d0 * vw / a_n

      do id = 2, n_degrees
        if ( with_TiTe ) then
          node_list%node(i)%values(1,id,var_u) = cfac * node_list%node(i)%values(1,id,var_Te)
        else
          node_list%node(i)%values(1,id,var_u) = cfac * 0.5d0 * node_list%node(i)%values(1,id,var_T)
        endif
      enddo
    else
      ! --- gauge: one constant over the whole wall, no tangential gradient. Dirichlet freezes
      ! --- these, so this is the value they hold for the run.
      node_list%node(i)%values(1,1,var_u) = u_wall
      do id = 2, n_degrees
        node_list%node(i)%values(1,id,var_u) = 0.d0
      enddo
    endif

    if ( is_sheath_type ) then
      n_loc(1) = n_loc(1) + 1.d0
    else
      n_loc(2) = n_loc(2) + 1.d0
    endif

  enddo

  call MPI_Reduce(n_loc, n_glo, 2, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if ( my_id .eq. 0 ) then
    write(*,'(A,i0,A,i0,A)')                                                        &
      ' SHEATH: sheath_init_u set u to the floating potential on ', nint(n_glo(1)),  &
      ' sheath nodes and ', nint(n_glo(2)), ' gauge nodes'
    if ( sheath_init_u_all ) write(*,'(A,es12.4,A)')                                &
      '         gauge nodes held at the CONSTANT wall potential u = ', u_wall,       &
      ' (equipotential, zero tangential gradient)'
  endif

end subroutine sheath_init_potential


!> Record psi's degrees of freedom at the start of the run, for the resistive wall relaxation.
!! Called once from jorek2_main after the restart is in place.
subroutine sheath_store_psi0(node_list)

  use mod_parameters
  use data_structure

  implicit none
  type (type_node_list), intent(in) :: node_list

  integer :: i, id

  if ( allocated(sheath_psi0) ) deallocate(sheath_psi0)
  allocate( sheath_psi0(n_degrees, node_list%n_nodes) )
  sheath_psi0 = 0.d0

  do i = 1, node_list%n_nodes
    do id = 1, n_degrees
      sheath_psi0(id,i) = node_list%node(i)%values(1,id,var_psi)
    enddo
  enddo

end subroutine sheath_store_psi0

end module mod_sheath_diag
