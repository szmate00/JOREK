!> A/B diagnostic for the Mach-1 (Bohm) boundary condition.
!!
!! PURPOSE. JOREK imposes the sheath entrance condition through Mach1BC in
!! mod_boundary_conditions, which sets
!!
!!     v_par = +-factor*c_s  +  factor * R^2 * (du/dn) / (dpsi/dn)
!!
!! The second term is a drift correction and is UNDOCUMENTED - it appears
!! identically in model401, model502 and model600 with no comment, and the only
!! commit touching it is "model600 bcs". SOLPS-ITER's drift-compatible
!! Bohm-Chodura condition (manual 3.0.9, p.411) instead requires
!!
!!     poloidal(v_par) + poloidal(v_ExB)  >=  poloidal(c_s)
!!
!! i.e. an added v_ExB.n / b_n on the parallel velocity.
!!
!! Working the geometry: iv_dir indexes the TANGENTIAL derivative (t_constant
!! boundary gives iv_dir = 2 = d/ds, and s runs along the boundary there), so
!! u0_b and ps0_b are tangential. With v_ExB = R grad(u) x e_phi one gets
!! v_ExB.n = -R du/dt and dpsi/dt = -R B_pol.n, hence
!!
!!     R^2 * (du/dt) / (dpsi/dt)  =  v_ExB.n / (B_pol.n)  =  v_ExB.n / (b_n |B|)
!!
!! which is the SOLPS correction up to normalisation. So JOREK very likely
!! ALREADY has the drift-compatible condition, and this module is a verification
!! rather than a discovery. Note it also means the term vanishes identically
!! wherever dirichlet%u holds - that pins the value and the TANGENTIAL
!! derivative, so du/dt = 0 - and therefore the whole question only arises on a
!! boundary where u is free, i.e. under the sheath boundary condition.
!!
!! This module answers that by measuring BOTH at the same nodes, imposing
!! NEITHER. It changes no behaviour: with diag_mach1 = .false. (default) nothing
!! is accumulated or printed.
!!
!! WHAT TO READ. The ratio's SPREAD matters more than its value, because the
!! Vpar normalisation may leave a constant factor the derivation above does not
!! pin down:
!!
!!   min ~ max, any value   the two are the SAME quantity differing by a constant
!!                          normalisation. Nothing to implement.
!!   ratio varies           genuinely different quantities, and the difference is
!!                          what an implementation would have to supply.
!!   ratio ~ 0              the correction is absent (expected wherever
!!                          dirichlet%u holds, since du/dt is then pinned to 0).
module mod_mach1_diag

  implicit none
  private
  public :: mach1_diag_reset, mach1_diag_add, mach1_diag_report

  ! --- Index 1 = inner (R < mach1_diag_R_split), 2 = outer. Weighted by |dl| so
  ! --- that long boundary elements are not under-represented; this is a nodal
  ! --- route, so there is no wetted area to weight by.
  real*8, save :: md_w(2)      = 0.d0   !< sum of weights
  real*8, save :: md_cs(2)     = 0.d0   !< sum w*c_s
  real*8, save :: md_dj(2)     = 0.d0   !< sum w*|JOREK drift term|
  real*8, save :: md_ds(2)     = 0.d0   !< sum w*|SOLPS drift term|
  real*8, save :: md_ratio(2)  = 0.d0   !< sum w*(JOREK/SOLPS), over resolved nodes only
  real*8, save :: md_rw(2)     = 0.d0   !< sum of w over those nodes, so the mean is not
                                        !< diluted by the ones where no ratio was formed
  real*8, save :: md_rmax(2)   = -1.d30 !< max of that ratio
  real*8, save :: md_rmin(2)   =  1.d30 !< min of that ratio
  real*8, save :: md_vn(2)     = 0.d0   !< sum w*(v_ExB . n), + = into the wall
  real*8, save :: md_vt(2)     = 0.d0   !< sum w*(v_ExB . t)
  real*8, save :: md_bn(2)     = 0.d0   !< sum w*|b_n|
  real*8, save :: md_djmax(2)  = 0.d0   !< max |JOREK drift| / c_s
  real*8, save :: md_dsmax(2)  = 0.d0   !< max |SOLPS drift| / c_s
  real*8, save :: md_bnmin(2)  = 1.d30  !< min |b_n|
  real*8, save :: md_rev(2)    = 0.d0   !< weight where the SOLPS form would reverse the outflow
  real*8, save :: md_n(2)      = 0.d0   !< node count

contains

!> Clear the accumulators. Once per boundary_conditions call, before the loop.
subroutine mach1_diag_reset()
  implicit none
  md_w=0.d0; md_cs=0.d0; md_dj=0.d0; md_ds=0.d0; md_ratio=0.d0
  md_vn=0.d0; md_vt=0.d0; md_bn=0.d0; md_djmax=0.d0; md_dsmax=0.d0
  md_bnmin=1.d30; md_rev=0.d0; md_n=0.d0
  md_rmax=-1.d30; md_rmin=1.d30; md_rw=0.d0
end subroutine mach1_diag_reset


!> Accumulate one boundary node.
!!
!! All velocities in the SAME units, i.e. multiplied through by Btot relative to
!! the Vpar variable, so that they are directly comparable with c_s.
!!
!! @param BigR   major radius of the node
!! @param cs     sound speed sqrt(gamma*(Ti+Te))
!! @param d_jrk  the JOREK drift term as it enters v_par: factor*R^2*u_n/psi_n
!! @param vE_n   E x B velocity component along the OUTWARD normal
!! @param vE_t   E x B velocity component along the boundary tangent
!! @param bn     B.n/|B|, the obliqueness (signed as in mod_boundary_conditions)
!! @param wgt    weight, the local boundary length element
subroutine mach1_diag_add(BigR, cs, d_jrk, vE_n, vE_t, bn, wgt)

  use phys_module, only: mach1_diag_R_split

  implicit none
  real*8, intent(in) :: BigR, cs, d_jrk, vE_n, vE_t, bn, wgt

  integer :: k
  real*8  :: d_slp, rat, absbn

  if ( cs .le. 0.d0 ) return

  k = 2
  if ( (mach1_diag_R_split .gt. 0.d0) .and. (BigR .lt. mach1_diag_R_split) ) k = 1

  ! --- SOLPS form: the E x B contribution reprojected onto the parallel
  ! --- velocity is v_ExB.n / b_n. Guard the grazing-incidence divergence; a
  ! --- vanishing b_n is exactly where this correction is unbounded, which is
  ! --- itself worth knowing and is why SOLPS clips it at 2*c_s.
  absbn = abs(bn)
  d_slp = 0.d0
  if ( absbn .gt. 1.d-8 ) d_slp = vE_n / absbn

  ! --- Only form the ratio where the denominator is a meaningful fraction of c_s;
  ! --- below that it is dominated by whatever noise is in du/dt and its spread
  ! --- would swamp the signal this diagnostic exists to show.
  rat = 0.d0
  if ( abs(d_slp) .gt. 1.d-3 * cs ) rat = d_jrk / d_slp

  !$omp critical (mach1_diag_accumulate)
  md_n(k)     = md_n(k)     + 1.d0
  md_w(k)     = md_w(k)     + wgt
  md_cs(k)    = md_cs(k)    + wgt * cs
  md_dj(k)    = md_dj(k)    + wgt * abs(d_jrk)
  md_ds(k)    = md_ds(k)    + wgt * abs(d_slp)
  if ( abs(d_slp) .gt. 1.d-3 * cs ) then
    md_ratio(k) = md_ratio(k) + wgt * rat
    md_rw(k)    = md_rw(k)    + wgt
    md_rmax(k)  = max(md_rmax(k), rat)
    md_rmin(k)  = min(md_rmin(k), rat)
  endif
  md_vn(k)    = md_vn(k)    + wgt * vE_n
  md_vt(k)    = md_vt(k)    + wgt * vE_t
  md_bn(k)    = md_bn(k)    + wgt * absbn
  md_djmax(k) = max(md_djmax(k), abs(d_jrk)/cs)
  md_dsmax(k) = max(md_dsmax(k), abs(d_slp)/cs)
  md_bnmin(k) = min(md_bnmin(k), absbn)
  ! --- Where |v_ExB.n / b_n| exceeds c_s the drift-compatible condition would
  ! --- demand a reversed parallel outflow. SOLPS clips its E x B contribution
  ! --- at 2*c_s*|b_x|, so a large fraction here means the clip is mandatory.
  if ( abs(d_slp) .gt. cs ) md_rev(k) = md_rev(k) + wgt
  !$omp end critical (mach1_diag_accumulate)

end subroutine mach1_diag_add


!> Reduce and print. Called after the element loop, from boundary_conditions.
subroutine mach1_diag_report(my_id)

  use mpi_mod
  use phys_module, only: mach1_diag_R_split

  implicit none
  integer, intent(in) :: my_id

  integer :: ierr, k
  real*8  :: loc(2,11), glo(2,11), lmax(2,3), gmax(2,3), lmin(2,2), gmin(2,2)
  real*8  :: w, cs, dj, ds, rat, vn, vt, bn, rev
  character(len=5) :: lbl(2)

  loc(:,1)=md_n;  loc(:,2)=md_w;  loc(:,3)=md_cs; loc(:,4)=md_dj; loc(:,5)=md_ds
  loc(:,6)=md_ratio; loc(:,7)=md_vn; loc(:,8)=md_vt; loc(:,9)=md_bn; loc(:,10)=md_rev
  loc(:,11)=md_rw
  lmax(:,1)=md_djmax; lmax(:,2)=md_dsmax; lmax(:,3)=md_rmax
  lmin(:,1)=md_bnmin; lmin(:,2)=md_rmin

  call MPI_Reduce(loc,  glo,  22, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(lmax, gmax,  6, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(lmin, gmin,  4, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)

  if ( my_id .ne. 0 ) return
  if ( sum(glo(:,1)) .lt. 0.5d0 ) return

  lbl(1) = 'INNER'; lbl(2) = 'OUTER'

  write(*,'(A)') ' MACH1 A/B: v_par = +-factor*c_s + JOREK_drift ;  SOLPS form adds v_ExB.n/b_n'
  do k = 1, 2
    if ( glo(k,1) .lt. 0.5d0 ) cycle
    if ( (mach1_diag_R_split .le. 0.d0) .and. (k .eq. 1) ) cycle
    w = glo(k,2)
    if ( w .le. 0.d0 ) cycle
    cs=glo(k,3)/w; dj=glo(k,4)/w; ds=glo(k,5)/w
    rat = 0.d0
    if ( glo(k,11) .gt. 0.d0 ) rat = glo(k,6)/glo(k,11)
    vn=glo(k,7)/w; vt=glo(k,8)/w; bn=glo(k,9)/w; rev=glo(k,10)/w

    write(*,'(A,A,A,i0,A,es10.3,A)')                                              &
      '   ', lbl(k), ' target: ', nint(glo(k,1)), ' nodes,  c_s=', cs, ' (mean)'
    write(*,'(A,es10.3,A,f7.2,A,f7.2,A)')                                         &
      '     JOREK drift term  = ', dj, '   mean/c_s=', 100.d0*dj/max(cs,1.d-30),  &
      ' %   max/c_s=', 100.d0*gmax(k,1), ' %'
    write(*,'(A,es10.3,A,f7.2,A,f7.2,A)')                                         &
      '     SOLPS drift term  = ', ds, '   mean/c_s=', 100.d0*ds/max(cs,1.d-30),  &
      ' %   max/c_s=', 100.d0*gmax(k,2), ' %'
    write(*,'(A,f9.3,A,f9.3,A,f9.3,A)')                                           &
      '     ratio JOREK/SOLPS = ', rat, '   min=', gmin(k,2), '  max=', gmax(k,3), &
      '   (constant => same quantity)'
    write(*,'(A,f6.2,A)') '     ratio resolved on ', 100.d0*glo(k,11)/w,                &
      ' % of the boundary weight (|SOLPS drift| > 1e-3 c_s)'
    write(*,'(A,es10.3,A,es10.3,A)')                                              &
      '     v_ExB.n (into wall)=', vn, '   v_ExB.t=', vt, ' [same units as c_s]'
    write(*,'(A,f8.4,A,f8.4,A,f6.2,A)')                                           &
      '     |b_n| mean=', bn, '  min=', gmin(k,1), '   would-reverse fraction=',   &
      100.d0*rev, ' %'
  enddo

end subroutine mach1_diag_report

end module mod_mach1_diag
