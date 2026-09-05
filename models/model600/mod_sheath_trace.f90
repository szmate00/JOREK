!> Weak sheath trace moments, accumulated from every incident local element.
!! Row owners have every incident element under distribute_nodes_elements.
!! Do not all-reduce the replicated contributions. Only owned rows are applied
!! and included in the global row statistics.
!!
!! Replacement clears the full scalar row, then divides by its largest absolute
!! coefficient. This is row equilibration, NOT column/basis equilibration.
module mod_sheath_trace

  implicit none
  private

  public :: sheath_trace_reset, sheath_trace_add, sheath_trace_apply, sheath_trace_report

  integer, save :: st_capacity = 0 ! allocated from the boundary-node graph at construction
  !> Columns per row. Only the variables the characteristic actually depends on appear: zj (the
  !! unit diagonal), u, rho, Ti/Te (or T), optionally vpar, and the normal-derivative psi DOFs that
  !! change g(B.n) and |B|. w is absent because the characteristic does not depend on it.
  !! A shared trace DOF collects columns from BOTH adjacent edges: 3 distinct nodes x 2 trace DOFs
  !! x 5 variables = 30, and a junction where three boundary segments meet takes 4 nodes = 40.
  !! 64 leaves headroom; the overflow is fatal rather than silent, so being generous costs memory.
  integer, parameter :: st_max_col = 96

  integer, save :: st_owner_min = 1, st_owner_max = 0
  integer, save :: st_n = 0
  integer, allocatable, save :: st_row(:), st_equation(:)
  real*8,  allocatable, save :: st_D(:)
  real*8,  allocatable, save :: st_F(:)
  real*8,  allocatable, save :: st_D0(:)    !< UNWEIGHTED diagonal; W_a = st_D/st_D0 in [0,1]
  real*8,  allocatable, save :: st_Dv(:)    !< diagonal with the VALIDITY weight only (no u-fade), so
                                        !! W_valid = st_Dv/st_D0 answers "is the characteristic
                                        !! solvable here" without being contaminated by "is u free
                                        !! here". sheath_weak_wmin gates on this one.
  real*8,  allocatable, save :: st_Fd(:)    !< DIAGNOSTIC residual: raw (zj0-zj_sh), weight wk_wgt.
                                        !! NOT st_F, which is the trust-region-bounded residual
                                        !! carrying a different weight (wk_wrx) and is what the
                                        !! row actually asks for. Only st_Fd/st_S is comparable
                                        !! to the printed global |F_a/D_a|/|S_a/D_a|.
  real*8,  allocatable, save :: st_S(:)     !< scale S_a, so |Fd_a/S_a| is the row residual in j_sat
  real*8,  allocatable, save :: st_det(:)   !< node-frame determinant |x2 x x3| of the row's OWN
                                        !! node: |sin(angle)| between the two first-derivative
                                        !! DOF directions. A STATIC property of the mesh, so
                                        !! unlike every other gate here it is fixed from step 1
                                        !! and cannot respond to the solution.
  integer, allocatable, save :: st_bnd(:)   !< boundary type of the row's OWN node (not the element's)
  integer, allocatable, save :: st_nc(:)
  integer, allocatable, save :: st_col(:,:)
  integer, allocatable, save :: st_var(:,:)
  real*8,  allocatable, save :: st_val(:,:)

  logical, save :: st_over_row = .false.     !< ran out of rows
  logical, save :: st_over_col = .false.     !< ran out of columns in some row
  integer, save :: st_n_applied = 0
  integer, save :: st_n_skipped = 0
  integer, save :: st_n_detgate = 0          !< of those, how many the determinant gate took

contains

!> Clear the accumulator. Once per matrix construction, before the element loop.
subroutine sheath_trace_reset(capacity)
  implicit none
  integer, optional, intent(in) :: capacity
  integer :: requested
  requested = 8000 ! standalone/legacy caller fallback, not the production mesh limit
  if (present(capacity)) requested = max(1,capacity)
  if (requested /= st_capacity) then
    if (allocated(st_row)) deallocate(st_equation,st_row,st_D,st_F,st_D0,st_Dv,st_Fd,st_S,st_det,st_bnd,st_nc,st_col,st_var,st_val)
    allocate(st_equation(requested),st_row(requested),st_D(requested),st_F(requested),st_D0(requested),st_Dv(requested), &
             st_Fd(requested),st_S(requested),st_det(requested),st_bnd(requested),st_nc(requested), &
             st_col(st_max_col,requested),st_var(st_max_col,requested),st_val(st_max_col,requested))
    st_capacity=requested
  endif
  st_n = 0
  st_over_row = .false.
  st_over_col = .false.
  st_n_applied = 0
  st_n_skipped = 0
  st_n_detgate = 0
  st_owner_min = 1
  st_owner_max = 0
end subroutine sheath_trace_reset


!> Slot for a global row index, creating it on first use. Linear search backwards: adjacent
!! boundary edges share DOFs and are visited consecutively, so the hit is normally immediate.
integer function st_slot(irow, equation)
  implicit none
  integer, intent(in) :: irow, equation
  integer :: i

  do i = st_n, 1, -1
    if (st_row(i) == irow .and. st_equation(i) == equation) then
      st_slot = i
      return
    endif
  enddo

  if ( st_n .ge. st_capacity ) then
    st_over_row = .true.
    st_slot = 0
    return
  endif

  st_n            = st_n + 1
  st_row(st_n)    = irow
  st_equation(st_n) = equation
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
!! @param equation scalar equation being replaced (default var_zj; also used for weak Mach1)
subroutine sheath_trace_add(irow, ibnd, dD, dD0, dDv, dF, dFd, dS, dnod, nc, icol, ivar, vals, equation)
  use mod_parameters, only: var_zj
  implicit none
  integer, optional, intent(in) :: equation
  integer, intent(in) :: irow, ibnd, nc
  real*8,  intent(in) :: dD, dD0, dDv, dF, dFd, dS, dnod
  integer, intent(in) :: icol(nc), ivar(nc)
  real*8,  intent(in) :: vals(nc)

  integer :: is, ic, jc, row_variable
  logical :: found

  ! --- The element loop that calls this is OpenMP-threaded (construct_matrix_mod passes omp_tid),
  ! --- and every st_* array here is module-level shared state that this routine both searches and
  ! --- extends. Without the critical section two threads race on st_n and write past each other,
  ! --- producing corrupted row indices - which is a crash, not a wrong answer. sheath_diag_add
  ! --- protects its own accumulators the same way. Boundary elements are a small fraction of the
  ! --- mesh, so serialising here costs little.
  !$omp critical (sheath_trace_accumulate)

  row_variable = var_zj
  if (present(equation)) row_variable = equation
  is = st_slot(irow, row_variable)

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


!> Replace owned current equations exactly, independently of a penalty magnitude.
subroutine sheath_trace_apply(in, zbig, index_min, index_max, a_mat, RHS_loc)
  use mod_parameters
  use data_structure, only: type_SP_MATRIX
  use mod_assembly, only: boundary_conditions_clear_row, boundary_conditions_add_one_entry, &
                         boundary_conditions_add_RHS
  use mpi_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, intent(in) :: in, index_min, index_max
  real*8, intent(in) :: zbig ! retained for call compatibility; no penalty is used
  type(type_SP_MATRIX), intent(inout) :: a_mat
  real*8, intent(inout) :: RHS_loc(*)
  integer :: is, ic, ierr
  real*8 :: row_norm

  st_owner_min = index_min
  st_owner_max = index_max
  st_n_applied = 0
  st_n_skipped = 0
  st_n_detgate = 0
  if (st_over_row .or. st_over_col) then
    write(*,*) 'ERROR: sheath trace storage exhausted; refusing a truncated boundary system.'
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    return
  endif
  do is=1,st_n
    if (st_row(is) < index_min .or. st_row(is) > index_max) cycle
    if (st_nc(is) <= 0 .or. st_D(is) <= 0.d0 .or. st_D0(is) <= 0.d0) then
      write(*,*) 'ERROR: unsupported sheath trace row (DOF,type,D,D0): ', &
                  st_row(is), st_bnd(is), st_D(is), st_D0(is)
      ! Freezing only zj strands the free u. Do not invent that fallback.
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      return
    endif
    if (.not. all(ieee_is_finite(st_val(1:st_nc(is),is))) .or. .not. ieee_is_finite(st_F(is))) then
      write(*,*) 'ERROR: nonfinite sheath trace row: ', st_row(is)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      return
    endif
    row_norm = maxval(abs(st_val(1:st_nc(is),is)))
    if (row_norm <= tiny(row_norm)) then
      write(*,*) 'ERROR: zero sheath equation at DOF ', st_row(is)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      return
    endif
    call boundary_conditions_clear_row(st_row(is), st_equation(is), in, index_min, index_max, a_mat)
    do ic=1,st_nc(is)
      call boundary_conditions_add_one_entry(st_row(is), st_equation(is), in, st_col(ic,is), st_var(ic,is), in, &
          st_val(ic,is)/row_norm, index_min, index_max, a_mat)
    enddo
    call boundary_conditions_add_RHS(st_row(is), st_equation(is), in, index_min, index_max, RHS_loc, &
                                    st_F(is)/row_norm, a_mat%i_tor_min, a_mat%i_tor_max)
    st_n_applied = st_n_applied+1
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

  loc(1) = dble(count(st_row(1:st_n) >= st_owner_min .and. st_row(1:st_n) <= st_owner_max))
  loc(2) = dble(st_n_applied)
  loc(3) = 0.d0
  if ( st_over_row .or. st_over_col ) loc(3) = 1.d0
  loc(4) = dble(st_n_skipped)
  loc(5) = dble(st_n_detgate)

  dloc(1) = 1.d30; dloc(2) = -1.d30
  do is=1,st_n
    if (st_row(is) < st_owner_min .or. st_row(is) > st_owner_max) cycle
    if (st_D(is) > 0.d0) dloc(1) = min(dloc(1),st_D(is))
    dloc(2) = max(dloc(2),st_D(is))
  enddo

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
    if (st_row(is) < st_owner_min .or. st_row(is) > st_owner_max) cycle
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
    '         sheath trace equations (current + Mach1): ', nint(glo(1)),                                &
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
    write(*,*) '       Check graph capacity / st_max_col in mod_sheath_trace.f90.'
    stop
  endif

end subroutine sheath_trace_report

end module mod_sheath_trace
