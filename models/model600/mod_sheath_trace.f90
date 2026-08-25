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

  integer, parameter :: st_max_row = 4000    !< trace DOFs per rank; 204 observed, so ample
  !> Columns per row. Only the variables the characteristic actually depends on appear: zj (the
  !! unit diagonal), u, rho and Ti/Te - or T without WITH_TiTe. vpar and w are absent because
  !! sheath_current takes neither, so their derivatives are identically zero; psi is absent by
  !! CHOICE, since zj_sat ~ 1/Btot and g_bn do depend on the field but sheath_current returns no
  !! d(zj_sh)/d(psi) - the geometry is frozen within a Newton iteration (quasi-Newton).
  !! A shared trace DOF collects columns from BOTH adjacent edges: 3 distinct nodes x 2 trace DOFs
  !! x 5 variables = 30, and a junction where three boundary segments meet takes 4 nodes = 40.
  !! 64 leaves headroom; the overflow is fatal rather than silent, so being generous costs memory.
  integer, parameter :: st_max_col = 64

  integer, save :: st_n = 0
  integer, save :: st_row(st_max_row)
  real*8,  save :: st_D(st_max_row)
  real*8,  save :: st_F(st_max_row)
  integer, save :: st_nc(st_max_row)
  integer, save :: st_col(st_max_col, st_max_row)
  integer, save :: st_var(st_max_col, st_max_row)
  real*8,  save :: st_val(st_max_col, st_max_row)

  logical, save :: st_over_row = .false.     !< ran out of rows
  logical, save :: st_over_col = .false.     !< ran out of columns in some row
  integer, save :: st_n_applied = 0
  integer, save :: st_n_skipped = 0

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
  st_nc(st_n)     = 0
  st_slot         = st_n

end function st_slot


!> Accumulate one element's contribution to one trace row.
!!
!! @param irow  global DOF index of the test function
!! @param dD    contribution to D_a = int N_a N_a dS
!! @param dF    contribution to F_a = int N_a (zj_sh - zj) dS   (note the sign: this is the RHS)
!! @param nc    number of columns in this contribution
!! @param icol  global DOF indices of the trial functions
!! @param ivar  variable index of each column
!! @param vals  contribution to int N_a (d(residual)/d x) N_b dS for each column
subroutine sheath_trace_add(irow, dD, dF, nc, icol, ivar, vals)
  implicit none
  integer, intent(in) :: irow, nc
  real*8,  intent(in) :: dD, dF
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

    st_D(is) = st_D(is) + dD
    st_F(is) = st_F(is) + dF

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

  use mod_parameters
  use data_structure, only: type_SP_MATRIX
  use mod_assembly,   only: boundary_conditions_add_one_entry, boundary_conditions_add_RHS

  implicit none
  integer,              intent(in)    :: in, index_min, index_max
  real*8,               intent(in)    :: zbig
  type(type_SP_MATRIX), intent(inout) :: a_mat
  real*8,               intent(inout) :: RHS_loc(*)

  integer :: is, ic
  real*8  :: sc, d_max

  st_n_applied = 0
  st_n_skipped = 0

  ! --- Degeneracy floor, relative to the largest diagonal on this rank. The spread between value
  ! --- and derivative rows is genuine element-size scaling (measured 2.6e10 on this mesh) and
  ! --- normalising it is exactly the point, so the floor must be far below that - it is here only
  ! --- to catch a row that received essentially no boundary support, where 1/D_a would manufacture
  ! --- an enormous row out of numerical noise.
  d_max = 0.d0
  do is = 1, st_n
    d_max = max(d_max, st_D(is))
  enddo

  do is = 1, st_n

    ! --- A row with no diagonal received no boundary-edge contribution at all; writing zbig into
    ! --- it would create a singular equation, so leave the assembled row alone.
    if ( st_D(is) .le. 0.d0 ) cycle
    if ( st_D(is) .lt. 1.d-14 * d_max ) then
      st_n_skipped = st_n_skipped + 1
      cycle
    endif

    sc = 1.d0 / st_D(is)

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

    st_n_applied = st_n_applied + 1

  enddo

end subroutine sheath_trace_apply


!> One line on what was written, and a hard stop on overflow - a silently truncated row would be
!! a wrong equation rather than a missing one.
subroutine sheath_trace_report(my_id)

  use mpi_mod

  implicit none
  integer, intent(in) :: my_id

  integer :: ierr
  real*8  :: loc(4), glo(4), dloc(2), dglo(2)

  loc(1) = dble(st_n)
  loc(2) = dble(st_n_applied)
  loc(3) = 0.d0
  if ( st_over_row .or. st_over_col ) loc(3) = 1.d0
  loc(4) = dble(st_n_skipped)

  dloc(1) = 1.d30; dloc(2) = -1.d30
  if ( st_n .gt. 0 ) then
    dloc(1) = minval(st_D(1:st_n), mask = st_D(1:st_n) .gt. 0.d0)
    dloc(2) = maxval(st_D(1:st_n))
  endif

  call MPI_Reduce(loc,  glo,  4, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(dloc(1), dglo(1), 1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(dloc(2), dglo(2), 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

  if ( my_id .ne. 0 ) return

  write(*,'(A,i0,A,i0,A,i0,A,es9.2,A,es9.2)')                                   &
    '         sheath trace rows: ', nint(glo(1)), ' accumulated, ',              &
    nint(glo(2)), ' replaced, ', nint(glo(4)), ' below the degeneracy floor;',   &
    '  D min=', dglo(1), ' max=', dglo(2)

  if ( glo(3) .gt. 0.d0 ) then
    write(*,*) 'ERROR: the sheath trace accumulator overflowed. A truncated row is a WRONG'
    write(*,*) '       equation, not a missing one, so this cannot be allowed to continue.'
    write(*,*) '       Raise st_max_row / st_max_col in mod_sheath_trace.f90.'
    stop
  endif

end subroutine sheath_trace_report

end module mod_sheath_trace
