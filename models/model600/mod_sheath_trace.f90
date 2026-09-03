!> Weak (Galerkin) sheath characteristic on the boundary trace space.
!!
!! The nodal route (bcs%sheath_zj) imposes zj = zj_sh POINTWISE at boundary nodes. Measured, it
!! does that essentially exactly - |zj-zj_sh|/|zj_sat| reaches 1.8e-4 at the nodes - and yet the
!! integrated currents still differ by 33% on the inner target, because the cubic trace BETWEEN
!! the nodes is not controlled. The weak residual
!!
!!     F_a = integral over Gamma of  N_a * (zj_sh - zj) dS
!!
!! stays at O(1) while the nodal residual falls by four decades, so a projection onto the trace
!! space has real work to do that the nodal constraint is not doing.
!!
!! This module accumulates that projection and replaces the boundary zj rows with it.
!!
!! WHY REPLACE RATHER THAN PENALISE. zj = Delta*psi is integrated by parts and its surface term
!! is refused (its Jacobian needs normal-derivative columns the trial loop cannot produce). That
!! is harmless while dirichlet%zj freezes the trace - the code says so - but the weak route has to
!! release it, and the equation is then incomplete at the boundary. Adding a penalty to a wrong
!! equation cannot fix it: measured, beta = 1e-2 left the incomplete equation dominant and beta = 1
!! drove a period-2 divergence. Writing the row with zbig annihilates it instead, exactly as every
!! Dirichlet in this code already does.
!!
!! WHY ROW NORMALISATION IS MANDATORY. A Galerkin trace block has internal scale structure that a
!! pointwise Dirichlet does not: its mass matrix goes as h (value x value), h^2 (value x
!! derivative) and h^3 (derivative x derivative), so at h ~ 1 mm the block spans ~1e6 internally.
!! A uniform coefficient on it therefore cannot work at any magnitude - beta = 1e9 died in one
!! step. Dividing each row by its own diagonal D_a = int N_a N_a dS makes every row O(1), after
!! which a single zbig makes them all uniformly dominant and no penalty parameter is needed.
!!
!! KNOWN LIMITATION (MPI). Rows are written only by the rank that owns them, so a trace DOF whose
!! two adjacent boundary edges live on different ranks is assembled from the local edge only. Both
!! F_a and D_a lose the same contribution, so the NORMALISED row is still a valid weighted-residual
!! statement - just tested against a function supported on one edge rather than two. A full fix
!! needs a halo exchange of the accumulator. Affected DOFs are those on an MPI partition boundary
!! that also lie on the sheath boundary.
module mod_sheath_trace

  implicit none
  private

  public :: sheath_trace_reset, sheath_trace_add, sheath_trace_apply, sheath_trace_report

  integer, parameter :: st_max_row = 8000    !< trace DOFs per rank; 424 observed on two types
  !> Columns per row. Only the variables the characteristic actually depends on appear: zj (the
  !! unit diagonal), u, rho, Ti/Te (or T), optionally vpar, and the normal-derivative psi DOFs that
  !! change g(B.n) and |B|. w is absent because the characteristic does not depend on it.
  !! A shared trace DOF collects columns from BOTH adjacent edges: 3 distinct nodes x 2 trace DOFs
  !! x 5 variables = 30, and a junction where three boundary segments meet takes 4 nodes = 40.
  !! 64 leaves headroom; the overflow is fatal rather than silent, so being generous costs memory.
  integer, parameter :: st_max_col = 96

  integer, save :: st_n = 0
  integer, save :: st_row(st_max_row)
  real*8,  save :: st_D(st_max_row)
  real*8,  save :: st_F(st_max_row)
  real*8,  save :: st_D0(st_max_row)    !< UNWEIGHTED diagonal; W_a = st_D/st_D0 in [0,1]
  real*8,  save :: st_Dv(st_max_row)    !< diagonal with the VALIDITY weight only (no u-fade), so
                                        !! W_valid = st_Dv/st_D0 answers "is the characteristic
                                        !! solvable here" without being contaminated by "is u free
                                        !! here". sheath_weak_wmin gates on this one.
  real*8,  save :: st_Fd(st_max_row)    !< DIAGNOSTIC residual: raw (zj0-zj_sh), weight wk_wgt.
                                        !! NOT st_F, which is the trust-region-bounded residual
                                        !! carrying a different weight (wk_wrx) and is what the
                                        !! row actually asks for. Only st_Fd/st_S is comparable
                                        !! to the printed global |F_a/D_a|/|S_a/D_a|.
  real*8,  save :: st_S(st_max_row)     !< scale S_a, so |Fd_a/S_a| is the row residual in j_sat
  real*8,  save :: st_det(st_max_row)   !< node-frame determinant |x2 x x3| of the row's OWN
                                        !! node: |sin(angle)| between the two first-derivative
                                        !! DOF directions. A STATIC property of the mesh, so
                                        !! unlike every other gate here it is fixed from step 1
                                        !! and cannot respond to the solution.
  integer, save :: st_bnd(st_max_row)   !< boundary type of the row's OWN node (not the element's)
  integer, save :: st_nc(st_max_row)
  integer, save :: st_col(st_max_col, st_max_row)
  integer, save :: st_var(st_max_col, st_max_row)
  real*8,  save :: st_val(st_max_col, st_max_row)

  logical, save :: st_over_row = .false.     !< ran out of rows
  logical, save :: st_over_col = .false.     !< ran out of columns in some row
  integer, save :: st_n_applied = 0
  integer, save :: st_n_skipped = 0
  integer, save :: st_n_detgate = 0          !< of those, how many the determinant gate took

contains

!> Clear the accumulator. Once per matrix construction, before the element loop.
subroutine sheath_trace_reset()
  implicit none
  st_n = 0
  st_over_row = .false.
  st_over_col = .false.
  st_n_applied = 0
end subroutine sheath_trace_reset


!> Slot for a global row index, creating it on first use. Linear search backwards: adjacent
!! boundary edges share DOFs and are visited consecutively, so the hit is normally immediate.
integer function st_slot(irow)
  implicit none
  integer, intent(in) :: irow
  integer :: i

  do i = st_n, 1, -1
    if ( st_row(i) .eq. irow ) then
      st_slot = i
      return
    endif
  enddo

  if ( st_n .ge. st_max_row ) then
    st_over_row = .true.
    st_slot = 0
    return
  endif

  st_n            = st_n + 1
  st_row(st_n)    = irow
  st_D(st_n)      = 0.d0
  st_F(st_n)      = 0.d0
  st_D0(st_n)     = 0.d0
  st_Dv(st_n)     = 0.d0
  st_Fd(st_n)     = 0.d0
  st_S(st_n)      = 0.d0
  ! --- Large, so a row that is never told its determinant is never gated: an unknown frame
  ! --- must recover the previous behaviour, not silently delete the row.
  st_det(st_n)    = 1.d30
  st_bnd(st_n)    = 0
  st_nc(st_n)     = 0
  st_slot         = st_n

end function st_slot


!> Accumulate one element's contribution to one trace row.
!!
!! @param irow  global DOF index of the test function
!! @param ibnd  boundary type of the node owning this row (nodes(i)%boundary, NOT bnd_type1)
!! @param dD    contribution to D_a = int N_a N_a dS   (carries wk_wgt)
!! @param dD0   the same integral with wk_wgt forced to 1, for the validity weight W_a
!! @param dDv   the same integral with the VALIDITY weight only, excluding the u-fade
!! @param dFd   DIAGNOSTIC residual (raw, weight wk_wgt) - report only, never written to the row
!! @param dS    contribution to the scale S_a, so |Fd_a/S_a| is the residual in units of j_sat
!! @param dnod  node-frame determinant of the row's own node (see st_det). A per-node constant,
!!              so every call for a given row carries the same value.
!! @param dF    contribution to F_a = int N_a (zj_sh - zj) dS   (note the sign: this is the RHS)
!! @param nc    number of columns in this contribution
!! @param icol  global DOF indices of the trial functions
!! @param ivar  variable index of each column
!! @param vals  contribution to int N_a (d(residual)/d x) N_b dS for each column
subroutine sheath_trace_add(irow, ibnd, dD, dD0, dDv, dF, dFd, dS, dnod, nc, icol, ivar, vals)
  implicit none
  integer, intent(in) :: irow, ibnd, nc
  real*8,  intent(in) :: dD, dD0, dDv, dF, dFd, dS, dnod
  integer, intent(in) :: icol(nc), ivar(nc)
  real*8,  intent(in) :: vals(nc)

  integer :: is, ic, jc
  logical :: found

  ! --- The element loop that calls this is OpenMP-threaded (construct_matrix_mod passes omp_tid),
  ! --- and every st_* array here is module-level shared state that this routine both searches and
  ! --- extends. Without the critical section two threads race on st_n and write past each other,
  ! --- producing corrupted row indices - which is a crash, not a wrong answer. sheath_diag_add
  ! --- protects its own accumulators the same way. Boundary elements are a small fraction of the
  ! --- mesh, so serialising here costs little.
  !$omp critical (sheath_trace_accumulate)

  is = st_slot(irow)

  if ( is .gt. 0 ) then

    st_D(is)  = st_D(is)  + dD
    st_D0(is) = st_D0(is) + dD0
    st_Dv(is) = st_Dv(is) + dDv
    st_F(is) = st_F(is) + dF
    st_Fd(is) = st_Fd(is) + dFd
    st_S(is)  = st_S(is)  + dS
    ! --- A per-node constant, so every contribution carries the same value and this is an
    ! --- assignment in all but name; min() is the conservative reading if that ever changes.
    st_det(is) = min(st_det(is), dnod)
    ! --- The row's OWN node type, not the element's bnd_type1: at a seam between two covered
    ! --- types the element label is ambiguous and every per-type number built from it is wrong.
    st_bnd(is) = max(st_bnd(is), ibnd)

  do ic = 1, nc
    found = .false.
    do jc = 1, st_nc(is)
      if ( (st_col(jc,is) .eq. icol(ic)) .and. (st_var(jc,is) .eq. ivar(ic)) ) then
        st_val(jc,is) = st_val(jc,is) + vals(ic)
        found = .true.
        exit
      endif
    enddo
    if ( .not. found ) then
      if ( st_nc(is) .ge. st_max_col ) then
        st_over_col = .true.
        cycle
      endif
      st_nc(is)              = st_nc(is) + 1
      st_col(st_nc(is), is)  = icol(ic)
      st_var(st_nc(is), is)  = ivar(ic)
      st_val(st_nc(is), is)  = vals(ic)
    endif
  enddo

  endif

  !$omp end critical (sheath_trace_accumulate)

end subroutine sheath_trace_add


!> Normalise every accumulated row by its own diagonal and write it into the matrix with zbig,
!! replacing the incomplete boundary zj equation. Called after the element loop, from
!! boundary_conditions, which owns a_mat and RHS_loc.
subroutine sheath_trace_apply(in, zbig, index_min, index_max, a_mat, RHS_loc)
  use phys_module, only: sheath_weak_wmin, sheath_weak_detmin

  use mod_parameters
  use data_structure, only: type_SP_MATRIX
  use mod_assembly,   only: boundary_conditions_add_one_entry, boundary_conditions_add_RHS

  implicit none
  integer,              intent(in)    :: in, index_min, index_max
  real*8,               intent(in)    :: zbig
  type(type_SP_MATRIX), intent(inout) :: a_mat
  real*8,               intent(inout) :: RHS_loc(*)

  integer :: is, ic
  real*8  :: sc

  st_n_applied = 0
  st_n_skipped = 0
  st_n_detgate = 0

  ! --- Degeneracy floor, relative to the largest diagonal on this rank. The spread between value
  ! --- and derivative rows is genuine element-size scaling (measured 2.6e10 on this mesh) and
  ! --- normalising it is exactly the point, so the floor must be far below that - it is here only
  ! --- to catch a row that received essentially no boundary support, where 1/D_a would manufacture
  ! --- an enormous row out of numerical noise.
  do is = 1, st_n

    ! --- THE VALIDITY GATE. wk_wgt multiplies D_a, F_a and every J_ab identically, and the row is
    ! --- written as J/D and F/D, so a wk_wgt uniform over a row's support cancels EXACTLY - it
    ! --- cannot fade a row that is written with zbig, it only shrinks D_a. Measured: the divertor
    ! --- leg-end faces, whose whole support sits outside the characteristic's validity range,
    ! --- diverge in 4-8 steps with D collapsing 8-10 orders, while the near-tangential wall at
    ! --- 42 % gated area is fine, because there the weight VARIES within each row and does not
    ! --- cancel. Row replacement is binary, so the fade has to be the decision whether to write
    ! --- the row at all. W_a = D_a/D0_a is the mean validity weight over this row's own support:
    ! --- dimensionless, mesh independent, and identical on every rank - unlike the 1e-14*d_max
    ! --- floor it replaces, which was RANK-LOCAL, so whether a given DOF was enforced depended on
    ! --- which other rows happened to land on that rank and on which types were enabled.
    ! --- DEGENERACY FLOOR, RESTORED. This gate replaced a `st_D < 1e-14*d_max` floor, and at the
    ! --- default sheath_weak_wmin = 0 the test reads `st_D < 0`, which the check above already
    ! --- excludes - so the floor was silently unreachable. Its purpose, per the original comment,
    ! --- is to catch a row that received essentially no boundary support, where 1/D_a manufactures
    ! --- an enormous row out of numerical noise. Keyed on the row's OWN unweighted support D0_a it
    ! --- is rank-independent and treats value and derivative rows alike, which the old d_max form
    ! --- did not. wmin is an optional STRONGER gate on top; it can never disable this.
    ! --- Two separate tests. The 1e-12 degeneracy floor is on st_D, because it protects the
    ! --- 1/D_a normalisation and that uses st_D. The user-facing validity gate is on st_Dv, which
    ! --- excludes the u-fade: a row beside a Dirichlet-u seam is faded deliberately and gating it
    ! --- as though its physics had failed would freeze the seam rows first, for the wrong reason.
    ! --- THE FRAME GATE. Every other test here keys on the SOLUTION, so it can close in
    ! --- response to the very divergence it exists to prevent - measured, sheath_weak_wmin =
    ! --- 0.5 took a 305-step case to 2 by deleting 706 of 920 rows at step 2. This one keys on
    ! --- a STATIC property of the mesh: the node-frame determinant is fixed before the run
    ! --- starts, identical on every rank, and the exact set of rows it removes is knowable in
    ! --- advance (util/check_boundary_frames.py prints it from a restart file). It cannot feed
    ! --- back.
    ! ---
    ! --- What it does NOT do is restore a boundary condition on u at the nodes it gates: on a
    ! --- weak type dirichlet%u and natural%u are both required .false., so u there is left to
    ! --- the interior equation and its neighbours. On boundary type 5 that is 8 rows among
    ! --- 354; on type 4 it would be 37 of 88, and a type fragmented that far may fail for that
    ! --- reason instead. Read the gated count in the report before reading a null result.
    if ( st_D(is) .le. 0.d0 .or. st_D0(is) .le. 0.d0 .or.                          &
         st_D(is)  .lt. 1.d-12 * st_D0(is) .or.                                    &
         st_Dv(is) .lt. sheath_weak_wmin * st_D0(is) .or.                          &
         ( sheath_weak_detmin .gt. 0.d0 .and.                                      &
           st_det(is) .lt. sheath_weak_detmin ) ) then
      st_n_skipped = st_n_skipped + 1
      if ( sheath_weak_detmin .gt. 0.d0 .and. st_det(is) .lt. sheath_weak_detmin ) &
        st_n_detgate = st_n_detgate + 1

      ! --- A rejected sheath projection still needs a controlled boundary row. Falling through
      ! --- to the assembled current-definition row changes the boundary condition abruptly as W
      ! --- crosses the floor and was the last step in every observed D-collapse runaway. Freeze
      ! --- this increment instead: the unit self-column is nonsingular, bounds zj while the local
      ! --- characteristic is unusable, and recovers the pre-sheath Dirichlet behaviour for that
      ! --- row. Other assembled entries are O(1) beside zbig and therefore inactive, exactly as
      ! --- for the ordinary Dirichlet rows written by boundary_conditions.
      call boundary_conditions_add_one_entry(                                    &
             st_row(is), var_zj, in, st_row(is), var_zj, in, zbig,               &
             index_min, index_max, a_mat)
      call boundary_conditions_add_RHS(                                           &
             st_row(is), var_zj, in, index_min, index_max, RHS_loc, 0.d0,         &
             a_mat%i_tor_min, a_mat%i_tor_max)
      if ( st_row(is) .ge. index_min .and. st_row(is) .le. index_max )            &
        st_n_applied = st_n_applied + 1
      cycle
    endif

    ! --- Belt and braces on the same failure. sc multiplies st_val AND st_F uniformly, so clamping
    ! --- it rescales the row in the matrix without changing the equation it states. The fallback
    ! --- above makes this unreachable today; it is here so that a future change to the floor
    ! --- cannot resurrect the 1/D_a overflow.
    sc = 1.d0 / max( st_D(is), 1.d-12 * st_D0(is) )

    do ic = 1, st_nc(is)
      call boundary_conditions_add_one_entry(                                   &
             st_row(is), var_zj, in, st_col(ic,is), st_var(ic,is), in,           &
             zbig * st_val(ic,is) * sc, index_min, index_max, a_mat)
    enddo

    if ( in .eq. 1 ) then
      call boundary_conditions_add_RHS(                                          &
             st_row(is), var_zj, in, index_min, index_max, RHS_loc,              &
             zbig * st_F(is) * sc, a_mat%i_tor_min, a_mat%i_tor_max)
    else
      call boundary_conditions_add_RHS(                                          &
             st_row(is), var_zj, in, index_min, index_max, RHS_loc,              &
             0.d0, a_mat%i_tor_min, a_mat%i_tor_max)
    endif

    ! --- Count only rows this rank actually WROTE. boundary_conditions_add_one_entry
    ! --- silently returns for a row outside [index_min,index_max], so a halo DOF that
    ! --- local elements touched was being accumulated, counted as replaced, and written
    ! --- nowhere - inflating 'N accumulated, N replaced' into false reassurance.
    if ( st_row(is) .ge. index_min .and. st_row(is) .le. index_max ) &
      st_n_applied = st_n_applied + 1

  enddo

end subroutine sheath_trace_apply


!> One line on what was written, and a hard stop on overflow - a silently truncated row would be
!! a wrong equation rather than a missing one.
subroutine sheath_trace_report(my_id)

  use mpi_mod
  use phys_module, only: max_bnd_types

  implicit none
  integer, intent(in) :: my_id

  integer :: ierr, is, ib
  real*8  :: loc(5), glo(5), dloc(2), dglo(2)
  integer, parameter :: nbt = max_bnd_types
  ! --- cols 1-2 reduce with SUM, cols 3-6 with MAX, so each group is contiguous in memory
  ! --- (Fortran column-major) and travels in one call: 1 count, 2 sum(W), 3 max|Fd/S|,
  ! --- 4 -min(D), 5 max(D), 6 -min(W).
  real*8  :: tloc(nbt,8), tglo(nbt,8), r, w

  loc(1) = dble(st_n)
  loc(2) = dble(st_n_applied)
  loc(3) = 0.d0
  if ( st_over_row .or. st_over_col ) loc(3) = 1.d0
  loc(4) = dble(st_n_skipped)
  loc(5) = dble(st_n_detgate)

  dloc(1) = 1.d30; dloc(2) = -1.d30
  if ( st_n .gt. 0 ) then
    dloc(1) = minval(st_D(1:st_n), mask = st_D(1:st_n) .gt. 0.d0)
    dloc(2) = maxval(st_D(1:st_n))
  endif

  ! --- per-type: count (SUM), max|F/S| (MAX), D min (as -D, so one MAX reduce), D max (MAX)
  tloc = 0.d0
  tglo = 0.d0          ! MPI_Reduce fills this on rank 0 only; the debug build inits reals to
                       ! signalling NaN, so an untouched tglo on any other rank would trap
  tloc(:,4) = -1.d30
  tloc(:,5) = -1.d30
  tloc(:,6) = -1.d30
  tloc(:,7) = -1.d30
  tloc(:,8) = -1.d30
  do is = 1, st_n
    ib = st_bnd(is)
    if ( ib .lt. 1 .or. ib .gt. nbt ) cycle
    tloc(ib,1) = tloc(ib,1) + 1.d0
    if ( abs(st_S(is)) .gt. 1.d-300 ) then
      r = abs( st_Fd(is) / st_S(is) )
      tloc(ib,3) = max(tloc(ib,3), r)
    endif
    if ( st_D(is) .gt. 0.d0 ) then
      tloc(ib,4) = max(tloc(ib,4), -st_D(is))
      tloc(ib,5) = max(tloc(ib,5),  st_D(is))
    endif
    if ( st_D0(is) .gt. 0.d0 ) then
      w = st_D(is) / st_D0(is)
      tloc(ib,2) = tloc(ib,2) + w
      tloc(ib,6) = max(tloc(ib,6), -w)
      tloc(ib,7) = max(tloc(ib,7), -st_Dv(is) / st_D0(is))   ! min W_valid, as a negative
    endif
    ! --- min node-frame determinant, as a negative so it rides the same MAX reduce. 1.d30
    ! --- means no row on this type was ever told its frame, which cannot happen on the weak
    ! --- route but would otherwise print as a spurious 1e30.
    if ( st_det(is) .lt. 1.d29 ) tloc(ib,8) = max(tloc(ib,8), -st_det(is))
  enddo

  call MPI_Reduce(tloc(1,1), tglo(1,1), nbt*2, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(tloc(1,3), tglo(1,3), nbt*6, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

  call MPI_Reduce(loc,  glo,  5, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(dloc(1), dglo(1), 1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(dloc(2), dglo(2), 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

  if ( my_id .ne. 0 ) return

  tglo(:,4) = -tglo(:,4)      ! D min, W min and det min were reduced as negatives through MAX
  tglo(:,6) = -tglo(:,6)
  tglo(:,7) = -tglo(:,7)
  tglo(:,8) = -tglo(:,8)

  ! --- 12 edit descriptors, 12 output items. Counted, because getting this wrong is a run-time
  ! --- severe(61) that kills rank 0 mid-write and hangs every other rank in the next collective.
  write(*,'(A,i0,A,i0,A,i0,A,i0,A,es9.2,A,es9.2)')                              &
    '         sheath trace rows: ', nint(glo(1)),                                &
    ' accumulated, ',              nint(glo(2)),                                 &
    ' replaced, ',                 nint(glo(4)),                                 &
    ' below the floor (',          nint(glo(5)),                                 &
    ' frame-gated);  D min=',      dglo(1),                                      &
    ' max=',                       dglo(2)

  ! --- PER BOUNDARY TYPE. Every aggregate number in this campaign has been un-attributed:
  ! --- a single global |F/D| cannot say whether the residual lives on the target, on a corner,
  ! --- or on one type of a two-type configuration, and three parameter scans were run without
  ! --- that being answerable. |F_a/S_a| is the row residual in units of j_sat - near 0 is a
  ! --- converged row, order 19 is a row pinned at electron saturation with no voltage feedback.
  do ib = 1, nbt
    if ( nint(tglo(ib,1)) .le. 0 ) cycle
    ! --- 12 edit descriptors, 12 output items. Counted: a mismatch is a runtime severe(61)
    ! --- that kills rank 0 mid-write and hangs every other rank in the next collective.
    ! --- W is the validity weight D_a/D0_a. W near 1 = the characteristic is valid over the whole
    ! --- row; W -> 0 = the row is being written where the characteristic has no solution, which
    ! --- the wk_wgt weighting CANNOT fade because it cancels out of J/D and F/D.
    ! --- det min is the node-frame determinant, a STATIC mesh property: it is the same at
    ! --- every step of every run on this grid, so it belongs here as the label that says
    ! --- WHY a type is or is not solvable, not as something to watch evolve.
    write(*,'(A,i2,A,i0,A,f6.3,A,f6.3,A,f6.3,A,es9.2,A,es9.2,A,f6.3)')                   &
      '           bnd type ', ib,                                                  &
      ': rows=',              nint(tglo(ib,1)),                                     &
      '  W min=',             tglo(ib,6),                                           &
      ' mean=',               tglo(ib,2)/max(tglo(ib,1),1.d0),                      &
      '  Wv min=',            tglo(ib,7),                                           &
      '  |Fd/S|max=',         tglo(ib,3),                                           &
      '  D min=',             tglo(ib,4),                                           &
      '  det min=',           tglo(ib,8)
  enddo

  if ( glo(3) .gt. 0.d0 ) then
    write(*,*) 'ERROR: the sheath trace accumulator overflowed. A truncated row is a WRONG'
    write(*,*) '       equation, not a missing one, so this cannot be allowed to continue.'
    write(*,*) '       Raise st_max_row / st_max_col in mod_sheath_trace.f90.'
    stop
  endif

end subroutine sheath_trace_report

end module mod_sheath_trace
