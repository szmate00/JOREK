!> Geometry and demand diagnostics at the boundary Gauss points, for the weak sheath route.
!!
!! DIAGNOSTIC ONLY. Nothing here touches a matrix entry, a row or the RHS. Every number below is
!! read from quantities the assembly has already computed, so a build with this module active is
!! bit-identical to one without it and can be run on an existing restart for a single step.
!!
!! WHY IT EXISTS. Two questions ended the last four parameter runs and neither is answerable from
!! the existing output.
!!
!! 1. WHERE. The per-type report gives `|zj_sat| min` and an area-mean, but not the LOCATION of the
!!    offending point, nor whether it is the same point from run to run, nor whether the point that
!!    carries `max|zj/zj_sat|` is the same one that carries the minimum `|zj_sat|`. Measured on
!!    boundary type 5: min 2.47e-6 against an area-mean of 1.16e-3, a factor 470. The row's grip on
!!    u is `dzj_du = zj_sat*fp*dx_du`, i.e. LINEAR in zj_sat, so at that point the sheath has
!!    essentially no authority over u - which is why the run dies at ~310 steps whether
!!    `sheath_zj_ratio_max` deletes the row (D collapses, 317) or leaves it in place and inert
!!    (D constant, 307). Fixing that needs to know if it is one Gauss point or half the surface.
!!
!! 2. WHY. `zj_sat = c_sat*rho*(|g(b_n)|*cs + v_perp)/|B|` vanishes for TWO independent reasons -
!!    grazing incidence (`g -> 0`) and low density (`rho -> 0`). A weight keyed on `b_n` catches
!!    only the first. Recording rho, Ti, Te and `|B.n|/|B|` AT the offending point decides which,
!!    and therefore what a grazing blend must key on.
!!
!! It also records the element mapping Jacobian, which is the quantity the equations actually
!! experience:
!!
!!     qjac = |xjac| / (|x_s| * |x_t|)
!!
!! a dimensionless local conditioning measure in [0,1]. The nodal frame determinant used by
!! `sheath_weak_detmin` is only the CORNER LIMIT of this, up to the two `element%size` factors, so
!! `qjac` both validates that proxy and supersedes it. `sign(xjac) < 0` is an INVERTED element -
!! a hard mesh error, not a threshold to tune - and nothing in JOREK currently checks for it.
!!
!! Histograms rather than thresholds, deliberately: no threshold can be chosen before the
!! distribution is known, and choosing one here would repeat the mistake that produced
!! `sheath_weak_wmin`.
module mod_sheath_geom_diag

  use phys_module, only: max_bnd_types

  implicit none
  private

  public :: sheath_geom_reset, sheath_geom_add, sheath_geom_report

  integer, parameter :: nbt  = max_bnd_types
  integer, parameter :: nbin = 12          !< histogram bins, see sg_bin()
  integer, parameter :: npay = 12          !< payload doubles per located point, see sg_pack()

  ! --- Default-initialised at declaration, matching mod_sheath_diag. sheath_geom_reset is called
  ! --- under `.not. harmonic_matrix` while sheath_geom_add is not obviously guarded the same way,
  ! --- and the debug build initialises reals to signalling NaN - so an accumulator that were only
  ! --- ever set by reset() would trap on the first add() of a preconditioner construction.
  real*8, save :: sg_n   (nbt)      = 0.d0    !< Gauss points seen
  real*8, save :: sg_qsum(nbt)      = 0.d0    !< sum of qjac, for the mean
  real*8, save :: sg_qmin(nbt)      = -1.d30  !< min qjac, held NEGATED so one MAX reduce serves
  real*8, save :: sg_zsum(nbt)      = 0.d0    !< sum of |zj_sat|
  real*8, save :: sg_zmin(nbt)      = -1.d30  !< min |zj_sat|, negated
  real*8, save :: sg_neg (nbt)      = 0.d0    !< count of xjac < 0  -- inverted elements
  real*8, save :: sg_qh(nbin,nbt)   = 0.d0    !< histogram of qjac, linear bins over [0,1]
  real*8, save :: sg_zh(nbin,nbt)   = 0.d0    !< histogram of log10|zj_sat|, one decade per bin

  real*8, save :: sg_wd(npay) = 0.d0, sg_ws(npay) = 0.d0

contains

!> Bin index for a value in [lo,hi], clamped into [1,nbin].
integer function sg_bin(val, lo, hi)
  implicit none
  real*8, intent(in) :: val, lo, hi
  sg_bin = 1 + int( (val - lo) / max(hi - lo, 1.d-300) * dble(nbin) )
  sg_bin = max(1, min(nbin, sg_bin))
end function sg_bin


!> Clear the accumulators. Called once per matrix construction, before the element loop.
subroutine sheath_geom_reset()
  implicit none
  sg_n = 0.d0 ; sg_qsum = 0.d0 ; sg_zsum = 0.d0 ; sg_neg = 0.d0
  sg_qh = 0.d0 ; sg_zh = 0.d0
  ! --- minima are carried negated so that a single MPI_MAX reduce serves; -1.d30 is "nothing yet"
  sg_qmin = -1.d30 ; sg_zmin = -1.d30
  ! --- payload slot 1 is the metric and is compared, so it must start lower than any real value
  sg_wd = 0.d0 ; sg_wd(1) = -1.d30
  sg_ws = 0.d0 ; sg_ws(1) = -1.d30
end subroutine sheath_geom_reset


!> Pack one located point. Slot 1 is the comparison metric; the rest is context.
subroutine sg_pack(p, metric, R, Z, qjac, axjac, rho, Ti, Te, bn, cs, zj, zjsat)
  implicit none
  real*8, intent(out) :: p(npay)
  real*8, intent(in)  :: metric, R, Z, qjac, axjac, rho, Ti, Te, bn, cs, zj, zjsat
  p(1)  = metric ; p(2)  = R     ; p(3)  = Z    ; p(4)  = qjac
  p(5)  = axjac  ; p(6)  = rho   ; p(7)  = Ti   ; p(8)  = Te
  p(9)  = bn     ; p(10) = cs    ; p(11) = zj   ; p(12) = zjsat
end subroutine sg_pack


!> Accumulate one boundary Gauss point.
!!
!! @param ibnd   boundary type of the element's first node (bnd_type1, as sheath_diag_add uses)
!! @param xjac   the element mapping Jacobian AT this Gauss point, signed
!! @param xs,ys  d(x,y)/ds at this point
!! @param xt,yt  d(x,y)/dt at this point
!! @param R,Z    position
!! @param rho    density, corr_neg corrected
!! @param Ti,Te  temperatures, corr_neg corrected
!! @param bn     |B.n|/|B|, the DIMENSIONLESS grazing angle (note the per-type `geom:` line in the
!!               other diagnostic prints DIMENSIONAL B.n; this one is normalised)
!! @param cs     sound speed
!! @param zj     the plasma's current, zj0
!! @param zjsat  the saturation current the sheath can pass
subroutine sheath_geom_add(ibnd, xjac, xs, ys, xt, yt, R, Z, rho, Ti, Te, bn, cs, zj, zjsat)

  implicit none
  integer, intent(in) :: ibnd
  real*8,  intent(in) :: xjac, xs, ys, xt, yt, R, Z, rho, Ti, Te, bn, cs, zj, zjsat

  real*8  :: ns, nt, qjac, az, dem
  integer :: ib, k

  if ( ibnd .lt. 1 .or. ibnd .gt. nbt ) return

  ns = sqrt(xs*xs + ys*ys)
  nt = sqrt(xt*xt + yt*yt)
  if ( ns .gt. 1.d-300 .and. nt .gt. 1.d-300 ) then
    qjac = abs(xjac) / (ns * nt)
  else
    qjac = 0.d0
  endif
  qjac = min(qjac, 1.d0)          ! |sin| cannot exceed 1; guard against rounding

  az  = abs(zjsat)
  dem = 0.d0
  if ( az .gt. 1.d-300 ) dem = abs(zj) / az

  ib = ibnd

  ! --- The element loop that reaches here is OpenMP-threaded, and every sg_* array is shared
  ! --- module state that this routine both reads and updates. Same reasoning as
  ! --- sheath_trace_add: without the critical section two threads race and the counts are wrong.
  !$omp critical (sheath_geom_accumulate)

  sg_n(ib)    = sg_n(ib)    + 1.d0
  sg_qsum(ib) = sg_qsum(ib) + qjac
  sg_zsum(ib) = sg_zsum(ib) + az
  sg_qmin(ib) = max(sg_qmin(ib), -qjac)
  if ( az .gt. 0.d0 ) sg_zmin(ib) = max(sg_zmin(ib), -az)
  if ( xjac .lt. 0.d0 ) sg_neg(ib) = sg_neg(ib) + 1.d0

  k = sg_bin(qjac, 0.d0, 1.d0)
  sg_qh(k,ib) = sg_qh(k,ib) + 1.d0

  ! --- log10 histogram spanning 1e-12 .. 1e0, one decade per bin. Clamped at both ends by
  ! --- sg_bin, so an exact zero lands in bin 1 rather than producing -Infinity.
  if ( az .gt. 1.d-300 ) then
    k = sg_bin(log10(az), -12.d0, 0.d0)
  else
    k = 1
  endif
  sg_zh(k,ib) = sg_zh(k,ib) + 1.d0

  if ( dem .gt. sg_wd(1) ) &
    call sg_pack(sg_wd,  dem, R, Z, qjac, abs(xjac), rho, Ti, Te, bn, cs, zj, zjsat)
  if ( az .gt. 0.d0 .and. -az .gt. sg_ws(1) ) &
    call sg_pack(sg_ws, -az, R, Z, qjac, abs(xjac), rho, Ti, Te, bn, cs, zj, zjsat)

  !$omp end critical (sheath_geom_accumulate)

end subroutine sheath_geom_add


!> Reduce across ranks and print. Rank 0 only.
subroutine sheath_geom_report(my_id)

  use mpi_mod

  implicit none
  integer, intent(in) :: my_id

  integer :: ierr, ib, k
  real*8  :: gn(nbt), gqs(nbt), gqm(nbt), gzs(nbt), gzm(nbt), gng(nbt)
  real*8  :: gqh(nbin,nbt), gzh(nbin,nbt)
  ! --- scalars are passed to MPI as 1-element arrays throughout this file, matching the
  ! --- rest of the code base: mpi_mod uses mpif.h, whose implicit interfaces would accept
  ! --- a bare scalar, but an explicit-interface MPI module would reject it.
  real*8  :: mx(2), gmx(2), pay(npay), gpay(npay), nhit(1), gnhit(1)
  integer :: ip

  ! --- MPI_Reduce fills the receive buffers on rank 0 only; the debug build initialises reals to
  ! --- signalling NaN, so an untouched buffer on another rank would trap when read below.
  gn = 0.d0 ; gqs = 0.d0 ; gqm = 0.d0 ; gzs = 0.d0 ; gzm = 0.d0 ; gng = 0.d0
  gqh = 0.d0 ; gzh = 0.d0 ; gmx = 0.d0 ; gpay = 0.d0 ; gnhit = 0.d0 ; nhit = 0.d0

  call MPI_Reduce(sg_n,    gn,  nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_qsum, gqs, nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_zsum, gzs, nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_neg,  gng, nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_qmin, gqm, nbt, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_zmin, gzm, nbt, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_qh, gqh, nbin*nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_zh, gzh, nbin*nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  gqm = -gqm ; gzm = -gzm       ! the minima travelled negated

  if ( my_id .eq. 0 ) then
    if ( sum(gn) .le. 0.d0 ) return
    write(*,'(A)') '         --- sheath geometry / demand diagnostic (DIAGNOSTIC ONLY) ---'
    ! --- 12 edit descriptors, 12 output items. Counted: a mismatch is a runtime severe(61) that
    ! --- kills rank 0 mid-write and hangs every other rank in the next collective.
    do ib = 1, nbt
      if ( gn(ib) .le. 0.d0 ) cycle
      write(*,'(A,i3,A,i0,A,f6.4,A,f6.4,A,es10.3,A,es10.3,A,i0)')                       &
        '         bnd type', ib,                                                        &
        ' gauss=',      nint(gn(ib)),                                                    &
        '  qjac min=',  gqm(ib),                                                         &
        ' mean=',       gqs(ib)/gn(ib),                                                  &
        '  |zj_sat| min=', gzm(ib),                                                      &
        ' mean=',       gzs(ib)/gn(ib),                                                  &
        '  inverted=',  nint(gng(ib))
    enddo
  endif

  ! --- The two located points. Each rank keeps its own worst; the global max of the metric is
  ! --- shared, and the single rank holding it contributes its payload while every other rank
  ! --- contributes zeros, so a SUM reduce delivers the winner without a point-to-point exchange.
  ! --- nhit counts how many ranks matched, so an exact tie - which would corrupt the sum - is
  ! --- reported rather than printed as a plausible-looking wrong point.
  do ip = 1, 2
    if ( ip .eq. 1 ) then
      mx(1) = sg_wd(1) ; pay = sg_wd
    else
      mx(1) = sg_ws(1) ; pay = sg_ws
    endif
    call MPI_Allreduce(mx(1), gmx(1), 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    nhit(1) = 0.d0
    if ( mx(1) .ne. gmx(1) ) then
      pay = 0.d0
    else
      nhit(1) = 1.d0
    endif
    call MPI_Reduce(pay,  gpay,  npay, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(nhit(1), gnhit(1), 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if ( my_id .ne. 0 ) cycle
    if ( gmx(1) .le. -1.d29 ) cycle
    if ( nint(gnhit(1)) .ne. 1 ) then
      write(*,'(A,i0,A)') '         (', nint(gnhit(1)),                                     &
        ' ranks tied on this metric; the located point below is a SUM and is not meaningful)'
    endif
    if ( ip .eq. 1 ) then
      write(*,'(A,es10.3)') '         WORST DEMAND  max|zj/zj_sat| = ', gpay(1)
    else
      write(*,'(A,es10.3)') '         WEAKEST GRIP  min|zj_sat|    = ', -gpay(1)
    endif
    ! --- 16 edit descriptors, 16 output items. Counted, same reason as above.
    write(*,'(A,f8.4,A,f9.4,A,f6.4,A,es9.2,A,es9.2,A,es9.2,A,es9.2)')                    &
      '           R=',        gpay(2),  ' Z=',      gpay(3),                             &
      '  qjac=',              gpay(4),  '  |xjac|=', gpay(5),                            &
      '  rho=',               gpay(6),  '  Ti=',     gpay(7),                            &
      '  Te=',                gpay(8)
    write(*,'(A,es9.2,A,es9.2,A,es9.2,A,es9.2)')                                         &
      '           |B.n|/|B|=', gpay(9),  '  cs=',    gpay(10),                            &
      '  zj=',                 gpay(11), '  zj_sat=', gpay(12)
  enddo

  if ( my_id .ne. 0 ) return

  ! --- Histograms. No threshold is applied anywhere in this module on purpose: none can be
  ! --- justified before the distribution is known, and inventing one here would repeat exactly
  ! --- the mistake that produced sheath_weak_wmin.
  write(*,'(A)') '         qjac histogram, 10 % bins from 0 to 1 (low = badly conditioned):'
  do ib = 1, nbt
    if ( gn(ib) .le. 0.d0 ) cycle
    write(*,'(A,i3,A,12(1x,i6))') '           bnd type', ib, ':', (nint(gqh(k,ib)), k=1,nbin)
  enddo
  write(*,'(A)') '         |zj_sat| histogram, one decade per bin, 1e-12 .. 1e0:'
  do ib = 1, nbt
    if ( gn(ib) .le. 0.d0 ) cycle
    write(*,'(A,i3,A,12(1x,i6))') '           bnd type', ib, ':', (nint(gzh(k,ib)), k=1,nbin)
  enddo

end subroutine sheath_geom_report

end module mod_sheath_geom_diag
