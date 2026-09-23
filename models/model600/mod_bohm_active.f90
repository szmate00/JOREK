!> Per-node active set of the drift-compatible Bohm inequality (mach1_weak_drift_style = 1).
!!
!! The Vpar row imposes (B_pol.n)*Vpar >= cs*|b.n| - vE_r only where the flow is below it (SOLPS-ITER BCMOM=13,
!! non-marginal). The decision is made per wall node, from the node's own state, so that the two wall edges sharing
!! a node agree: a node is active if on any of its wall edges the condition is violated at the node. The flags are
!! collected during one matrix construction and used in the next one (the active set lags one time step); before
!! the first collection every node is active, i.e. the equality form. Nodes are identified by the global index of
!! their value DOF.
module mod_bohm_active

  implicit none
  private

  public :: bohm_active_begin, bohm_active_mark, bohm_active_node, bohm_active_finish

  real*8, allocatable, save :: act_used(:), act_new(:)
  logical, save :: have_flags = .false.

contains


!> Start collecting flags for a matrix construction.
subroutine bohm_active_begin(n_dof)
  implicit none
  integer, intent(in) :: n_dof
  if ( allocated(act_new) ) then
    if ( size(act_new) .ne. n_dof ) then
      deallocate(act_new, act_used)
      have_flags = .false.
    endif
  endif
  if ( .not. allocated(act_new) ) then
    allocate(act_new(n_dof), act_used(n_dof))
    act_used = 1.d0
  endif
  act_new = 0.d0
end subroutine bohm_active_begin


!> The Bohm inequality is violated at this node on one of its wall edges.
subroutine bohm_active_mark(idx)
  implicit none
  integer, intent(in) :: idx
  if ( .not. allocated(act_new) ) return
  if ( idx .lt. 1 .or. idx .gt. size(act_new) ) return
  !$omp critical (bohm_active)
  act_new(idx) = 1.d0
  !$omp end critical (bohm_active)
end subroutine bohm_active_mark


!> 1 if the Vpar row is imposed on this node's test functions, 0 if the node is free.
real*8 function bohm_active_node(idx)
  implicit none
  integer, intent(in) :: idx
  bohm_active_node = 1.d0
  if ( .not. have_flags ) return
  if ( idx .lt. 1 .or. idx .gt. size(act_used) ) return
  bohm_active_node = act_used(idx)
end function bohm_active_node


!> Combine the flags of all ranks; they are used from the next matrix construction on.
subroutine bohm_active_finish()
  use mpi_mod
  implicit none
  integer :: ierr
  if ( .not. allocated(act_new) ) return
  call MPI_ALLREDUCE(act_new, act_used, size(act_new), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  have_flags = .true.
end subroutine bohm_active_finish


end module mod_bohm_active
