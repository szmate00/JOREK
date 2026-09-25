!> Tangentially smoothed electron temperature for the wall rows (model600, sheath_Te_smooth > 0).
!!
!! The wall potential follows Te: on a floating segment u = C_T*Te exactly, and the sheath current row anchors u
!! near C_T*Te wherever the current is small. Its tangential gradient is then Lambda*dTe/ds, and the ExB normal flow
!! it drives, Lambda*dTe/ds/B, exceeded the Bohm outflow cs|b.n| by ~20x at the inner target corner in every run
!! (measured 1-4e4 m/s with a Te step of tens of eV across one or two elements, 2026-09-27). Here the Te the wall
!! rows see is filtered along the wall before it enters them: a 3-point [1/4 1/2 1/4] average over the two wall
!! neighbours of each wall node, applied npass times to every DOF component (value and derivatives) of the n = 1
!! harmonic. One pass is one element of smoothing: the length scale is the grid, not a free number. The plasma's
!! own Te is untouched; only the floating target and the characteristic's Te (in x and in j_sat's cs) use the
!! filtered field, lagged (evaluated at matrix construction), so their Te columns vanish.
!!
!! Wall neighbours: from the exterior element sides with both endpoints labelled (mod_boundary_edges when active,
!! else every such side), collected over the local elements and allreduced, so every rank has the whole chain.
!! Coincident node records (two records at one vertex) break the chain there and stay unfiltered.
module mod_wall_smooth

  implicit none
  private

  public :: wall_smooth_build, wall_smooth_active, wall_smooth_te

  real*8,  allocatable, save :: te_sm(:,:)     !< (4, n_nodes): filtered Te DOFs (n = 1 harmonic), raw off the wall
  logical, save :: active = .false.

contains


subroutine wall_smooth_build(element_list, node_list, npass)

  use data_structure,     only: type_element_list, type_node_list
  use mod_parameters,     only: var_Te, var_T
  use mod_model_settings, only: with_TiTe
  use mod_boundary_edges, only: boundary_edges_active, boundary_edge_is_exterior
  use mpi_mod

  implicit none
  type(type_element_list), intent(in) :: element_list
  type(type_node_list),    intent(in) :: node_list
  integer,                 intent(in) :: npass

  integer, allocatable :: nbr(:,:), nbr_g(:,:)
  real*8,  allocatable :: tmp(:,:), rn(:,:), rn_g(:,:)     ! the neighbour table is reduced as reals (one allreduce type everywhere)
  integer :: nn, e, s, a, b, i, kT, ipass, ierr

  nn = node_list%n_nodes
  kT = var_T
  if ( with_TiTe ) kT = var_Te
  active = .false.
  if ( nn .le. 0 .or. npass .le. 0 ) return

  ! --- wall neighbours of every wall node, from the exterior labelled sides of the local elements
  allocate( nbr(2,nn), nbr_g(2,nn) )
  nbr = 0
  do e = 1, element_list%n_elements
    do s = 1, 4
      a = element_list%element(e)%vertex(s)
      b = element_list%element(e)%vertex(mod(s,4)+1)
      if ( a .eq. b ) cycle
      if ( node_list%node(a)%boundary .eq. 0 .or. node_list%node(b)%boundary .eq. 0 ) cycle
      if ( boundary_edges_active() ) then
        if ( .not. boundary_edge_is_exterior(e, s) ) cycle
      endif
      call add_nbr(a, b)
      call add_nbr(b, a)
    enddo
  enddo
  allocate( rn(2,nn), rn_g(2,nn) )
  rn = dble(nbr)
  call MPI_ALLREDUCE(rn, rn_g, 2*nn, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  nbr_g = nint(rn_g)
  deallocate( rn, rn_g )

  ! --- the filter, npass times, on nodes with two wall neighbours
  if ( allocated(te_sm) ) deallocate(te_sm)
  allocate( te_sm(4,nn), tmp(4,nn) )
  do i = 1, nn
    te_sm(:,i) = node_list%node(i)%values(1,1:4,kT)
  enddo
  do ipass = 1, npass
    tmp = te_sm
    do i = 1, nn
      if ( nbr_g(1,i) .le. 0 .or. nbr_g(2,i) .le. 0 ) cycle
      te_sm(:,i) = 0.25d0 * tmp(:,nbr_g(1,i)) + 0.5d0 * tmp(:,i) + 0.25d0 * tmp(:,nbr_g(2,i))
    enddo
  enddo
  deallocate( nbr, nbr_g, tmp )
  active = .true.

contains

  subroutine add_nbr(p, q)
    integer, intent(in) :: p, q
    if ( nbr(1,p) .eq. q .or. nbr(2,p) .eq. q ) return
    if ( nbr(1,p) .eq. 0 ) then
      nbr(1,p) = q
    else if ( nbr(2,p) .eq. 0 ) then
      nbr(2,p) = q
    endif
  end subroutine add_nbr

end subroutine wall_smooth_build


logical function wall_smooth_active()
  implicit none
  wall_smooth_active = active
end function wall_smooth_active


!> The filtered Te DOF idof (1..4, n = 1 harmonic) of node inode, or raw when the filter is not built
real*8 function wall_smooth_te(inode, idof, raw)
  implicit none
  integer, intent(in) :: inode, idof
  real*8,  intent(in) :: raw
  wall_smooth_te = raw
  if ( .not. active ) return
  if ( inode .lt. 1 .or. inode .gt. size(te_sm,2) ) return
  wall_smooth_te = te_sm(idof, inode)
end function wall_smooth_te

end module mod_wall_smooth
