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
!!
!! FIVE CORRECTIONS after the first measurement, all of which changed a number that was quoted:
!!
!! * `bdotn` is ALREADY `B.n/|B|` (`mod_boundary_matrix_open.f90:421` divides by Btot; `:468`
!!   recovers the dimensional value as `sh_Bn = bdotn*Btot`, and `:521` compares it directly with
!!   `sheath_min_bn`). The first version divided by Btot a SECOND time, understating the grazing
!!   angle by ~2.5-2.8x. The memory note claiming "Bdotn is dimensional" refers to `sh_Bn` in the
!!   other diagnostic's `geom:` line, a different variable.
!! * The Chodura factor `g` is NOT proportional to `b_n` - `:429` makes it a tanh of `|bdotn|` -
!!   so reconstructing `|g|cs` from `b_n` was a second approximation on top of the first. `g_bn`
!!   and `Btot` are now recorded directly, which closes
!!   `zj_sat = c_sat*rho*(|g|cs + v_perp)/|B|` exactly from recorded quantities.
!! * Point counts are not area. Type 5's `dS` spans 4.0e-3 .. 1.6e-1, a factor of FORTY, so a
!!   fraction of Gauss points is not interchangeable with a fraction of wetted area - and area is
!!   what a coverage claim needs. Everything is now area-weighted, with point counts kept
!!   alongside for reference.
!! * REQUESTED is not ACTIVE. This routine is called from the `diag_sheath_zj` block, which is an
!!   OR over the edge's two endpoint types, so it sees Gauss points that are later faded out by
!!   `sheath_weak_ufade` or that sit on `sheath_weak_detmin`-gated nodes - which is why the first
!!   run reported type 9 at all, although type 9 contributes no trace rows whatsoever. Both
!!   distributions are now reported, the active one weighted by the effective trace weight.
!! * The raw SIGN of `xjac` is an orientation convention, not a mesh defect: it flips with
!!   `element_size_perp`, so it tracks the element's side family. The first version reported
!!   "inverted=1448" for all of type 5 and 0 for type 1, which is the signature of a convention
!!   and not of a broken mesh. Signs are now binned by (type, side family) and only a MIXTURE
!!   within one cell - a genuine fold - is flagged.
module mod_sheath_geom_diag

  use phys_module, only: max_bnd_types

  implicit none
  private

  public :: sheath_geom_reset, sheath_geom_add, sheath_geom_report

  integer, parameter :: nbt  = max_bnd_types
  integer, parameter :: nbin = 12          !< histogram bins, see sg_bin()
  integer, parameter :: npay = 16          !< payload doubles per located point, see sg_pack()

  ! --- Default-initialised at declaration, matching mod_sheath_diag. sheath_geom_reset is called
  ! --- under `.not. harmonic_matrix` while sheath_geom_add is not obviously guarded the same way,
  ! --- and the debug build initialises reals to signalling NaN - so an accumulator that were only
  ! --- ever set by reset() would trap on the first add() of a preconditioner construction.
  real*8, save :: sg_n   (nbt)      = 0.d0    !< Gauss points seen
  real*8, save :: sg_a   (nbt)      = 0.d0    !< REQUESTED area, sum of dS
  real*8, save :: sg_aact(nbt)      = 0.d0    !< ACTIVE area, sum of w_eff*dS
  real*8, save :: sg_qsum(nbt)      = 0.d0    !< area-weighted sum of qjac
  real*8, save :: sg_qmin(nbt)      = -1.d30  !< min qjac, held NEGATED so one MAX reduce serves
  real*8, save :: sg_qmna(nbt)      = -1.d30  !< min qjac over ACTIVE points only, negated
  real*8, save :: sg_zsum(nbt)      = 0.d0    !< area-weighted sum of |zj_sat|
  real*8, save :: sg_zmin(nbt)      = -1.d30  !< min |zj_sat|, negated
  !> Signs of xjac binned by (side family, type). The family is mod(side,2), which is exactly what
  !! selects direction(2) = 2 or 3, and the sign is constant within a family by construction. A
  !! MIXTURE inside one cell is a genuine fold; a difference BETWEEN cells is only the convention.
  real*8, save :: sg_sgn(2,2,nbt)   = 0.d0    !< (1=positive|2=negative, family, type)
  real*8, save :: sg_qh (nbin,nbt)  = 0.d0    !< area-weighted qjac histogram, REQUESTED
  real*8, save :: sg_qha(nbin,nbt)  = 0.d0    !< the same, ACTIVE (weighted by w_eff)
  real*8, save :: sg_zh (nbin,nbt)  = 0.d0    !< area-weighted log10|zj_sat| histogram, REQUESTED
  real*8, save :: sg_zha(nbin,nbt)  = 0.d0    !< the same, ACTIVE

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
  sg_n = 0.d0 ; sg_a = 0.d0 ; sg_aact = 0.d0
  sg_qsum = 0.d0 ; sg_zsum = 0.d0 ; sg_sgn = 0.d0
  sg_qh = 0.d0 ; sg_qha = 0.d0 ; sg_zh = 0.d0 ; sg_zha = 0.d0
  ! --- minima are carried negated so that a single MPI_MAX reduce serves; -1.d30 is "nothing yet"
  sg_qmin = -1.d30 ; sg_qmna = -1.d30 ; sg_zmin = -1.d30
  ! --- payload slot 1 is the metric and is compared, so it must start lower than any real value
  sg_wd = 0.d0 ; sg_wd(1) = -1.d30
  sg_ws = 0.d0 ; sg_ws(1) = -1.d30
end subroutine sheath_geom_reset


!> Pack one located point. Slot 1 is the comparison metric; the rest is context.
subroutine sg_pack(p, metric, R, Z, qjac, axjac, rho, Ti, Te, bn, cs, zj, zjsat, gbn, Btot, dS, weff)
  implicit none
  real*8, intent(out) :: p(npay)
  real*8, intent(in)  :: metric, R, Z, qjac, axjac, rho, Ti, Te, bn, cs, zj, zjsat, gbn, Btot, dS, weff
  p(1)  = metric ; p(2)  = R     ; p(3)  = Z    ; p(4)  = qjac
  p(5)  = axjac  ; p(6)  = rho   ; p(7)  = Ti   ; p(8)  = Te
  p(9)  = bn     ; p(10) = cs    ; p(11) = zj   ; p(12) = zjsat
  p(13) = gbn    ; p(14) = Btot  ; p(15) = dS   ; p(16) = weff
end subroutine sg_pack


!> Accumulate one boundary Gauss point.
!!
!! @param ibnd   boundary type of the element's first node (bnd_type1, as sheath_diag_add uses)
!! @param iside  the element's local side index; only its PARITY is used, and that is exactly what
!!               selects direction(2) = 2 or 3, i.e. the side family the xjac sign convention
!!               follows
!! @param xjac   the element mapping Jacobian AT this Gauss point, signed
!! @param xs,ys  d(x,y)/ds at this point
!! @param xt,yt  d(x,y)/dt at this point
!! @param R,Z    position
!! @param rho    density, corr_neg corrected
!! @param Ti,Te  temperatures, corr_neg corrected
!! @param bn     b_n = B.n/|B|, ALREADY DIMENSIONLESS - pass bdotn unmodified, do not divide by
!!               Btot again (mod_boundary_matrix_open.f90:421 has already done it)
!! @param cs     sound speed
!! @param gbn    the Chodura-Riemann factor g(b_n), signed. Recorded rather than reconstructed
!!               from bn, because it is a tanh of |bn| and not proportional to it.
!! @param Btot   |B|, so that zj_sat = c_sat*rho*(|g|cs + v_perp)/|B| closes exactly
!! @param zj     the plasma's current, zj0
!! @param zjsat  the saturation current the sheath can pass
!! @param dS     the surface measure at this point, ws*dl*R*2pi/n_plane - the same one
!!               sheath_diag_add is given, so areas here are comparable with that diagnostic
!! @param weff   the EFFECTIVE trace weight wk_wgt, after the validity gate and the u-fade. Points
!!               with weff = 0 are requested by the endpoint OR but contribute nothing to any row.
subroutine sheath_geom_add(ibnd, iside, xjac, xs, ys, xt, yt, R, Z, rho, Ti, Te, bn, cs, &
                           gbn, Btot, zj, zjsat, dS, weff)

  implicit none
  integer, intent(in) :: ibnd, iside
  real*8,  intent(in) :: xjac, xs, ys, xt, yt, R, Z, rho, Ti, Te, bn, cs
  real*8,  intent(in) :: gbn, Btot, zj, zjsat, dS, weff

  real*8  :: ns, nt, qjac, az, dem, wa, wact
  integer :: ib, k, fam

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

  ib   = ibnd
  fam  = 2 - mod(abs(iside), 2)             ! 1 for odd sides, 2 for even
  wa   = max(dS, 0.d0)                      ! requested area weight
  wact = wa * max(min(weff, 1.d0), 0.d0)    ! active area weight

  ! --- The element loop that reaches here is OpenMP-threaded, and every sg_* array is shared
  ! --- module state that this routine both reads and updates. Same reasoning as
  ! --- sheath_trace_add: without the critical section two threads race and the counts are wrong.
  !$omp critical (sheath_geom_accumulate)

  sg_n(ib)    = sg_n(ib)    + 1.d0
  sg_a(ib)    = sg_a(ib)    + wa
  sg_aact(ib) = sg_aact(ib) + wact
  sg_qsum(ib) = sg_qsum(ib) + qjac * wa
  sg_zsum(ib) = sg_zsum(ib) + az   * wa
  sg_qmin(ib) = max(sg_qmin(ib), -qjac)
  if ( wact .gt. 0.d0 ) sg_qmna(ib) = max(sg_qmna(ib), -qjac)
  if ( az .gt. 0.d0 )   sg_zmin(ib) = max(sg_zmin(ib), -az)

  if ( xjac .ge. 0.d0 ) then
    sg_sgn(1,fam,ib) = sg_sgn(1,fam,ib) + 1.d0
  else
    sg_sgn(2,fam,ib) = sg_sgn(2,fam,ib) + 1.d0
  endif

  k = sg_bin(qjac, 0.d0, 1.d0)
  sg_qh (k,ib) = sg_qh (k,ib) + wa
  sg_qha(k,ib) = sg_qha(k,ib) + wact

  ! --- log10 histogram spanning 1e-12 .. 1e0, one decade per bin. Clamped at both ends by
  ! --- sg_bin, so an exact zero lands in bin 1 rather than producing -Infinity.
  if ( az .gt. 1.d-300 ) then
    k = sg_bin(log10(az), -12.d0, 0.d0)
  else
    k = 1
  endif
  sg_zh (k,ib) = sg_zh (k,ib) + wa
  sg_zha(k,ib) = sg_zha(k,ib) + wact

  if ( dem .gt. sg_wd(1) ) &
    call sg_pack(sg_wd,  dem, R, Z, qjac, abs(xjac), rho, Ti, Te, bn, cs, zj, zjsat, &
                 gbn, Btot, dS, weff)
  if ( az .gt. 0.d0 .and. -az .gt. sg_ws(1) ) &
    call sg_pack(sg_ws, -az, R, Z, qjac, abs(xjac), rho, Ti, Te, bn, cs, zj, zjsat, &
                 gbn, Btot, dS, weff)

  !$omp end critical (sheath_geom_accumulate)

end subroutine sheath_geom_add


!> Reduce across ranks and print. Rank 0 only.
subroutine sheath_geom_report(my_id)

  use mpi_mod

  implicit none
  integer, intent(in) :: my_id

  integer :: ierr, ib, k, ip, f
  real*8  :: gn(nbt), ga(nbt), gaa(nbt), gqs(nbt), gqm(nbt), gqa(nbt), gzs(nbt), gzm(nbt)
  real*8  :: gsg(2,2,nbt)
  real*8  :: gqh(nbin,nbt), gqha(nbin,nbt), gzh(nbin,nbt), gzha(nbin,nbt)
  real*8  :: mx(2), gmx(2), pay(npay), gpay(npay), rnk(1), grnk(1)

  ! --- MPI_Reduce fills the receive buffers on rank 0 only; the debug build initialises reals to
  ! --- signalling NaN, so an untouched buffer on another rank would trap when read below.
  gn = 0.d0 ; ga = 0.d0 ; gaa = 0.d0 ; gqs = 0.d0 ; gqm = 0.d0 ; gqa = 0.d0
  gzs = 0.d0 ; gzm = 0.d0 ; gsg = 0.d0
  gqh = 0.d0 ; gqha = 0.d0 ; gzh = 0.d0 ; gzha = 0.d0
  gmx = 0.d0 ; gpay = 0.d0 ; grnk = 0.d0 ; rnk = 0.d0

  call MPI_Reduce(sg_n,    gn,  nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_a,    ga,  nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_aact, gaa, nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_qsum, gqs, nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_zsum, gzs, nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_qmin, gqm, nbt, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_qmna, gqa, nbt, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_zmin, gzm, nbt, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_sgn,  gsg, 4*nbt,    MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_qh,   gqh, nbin*nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_qha, gqha, nbin*nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_zh,   gzh, nbin*nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(sg_zha, gzha, nbin*nbt, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  gqm = -gqm ; gqa = -gqa ; gzm = -gzm       ! the minima travelled negated

  if ( my_id .eq. 0 ) then
    if ( sum(gn) .le. 0.d0 ) return
    write(*,'(A)') '         --- sheath geometry / demand diagnostic (DIAGNOSTIC ONLY) ---'
    ! --- REQUESTED area is every Gauss point the endpoint OR reaches; ACTIVE is weighted by the
    ! --- effective trace weight, so it excludes what sheath_weak_ufade fades and what
    ! --- sheath_weak_detmin gates. The two differ a lot and only ACTIVE describes the rows.
    ! --- 14 edit descriptors, 14 output items. Counted AND EXERCISED: a mismatch is a runtime
    ! --- severe(61) that kills rank 0 mid-write and hangs every other rank in the next
    ! --- collective. An earlier draft of this very line had 15 items against 14 descriptors.
    do ib = 1, nbt
      if ( gn(ib) .le. 0.d0 ) cycle
      ! --- gqa is the -1.d30 sentinel when a type has NO active points (measured: type 9 under
      ! --- sheath_weak_detmin), and f6.4 prints that as ******. Show a plain 0 and let the
      ! --- 'active=  0.0 %' beside it carry the meaning.
      if ( gqa(ib) .lt. 0.d0 .or. gqa(ib) .gt. 1.d0 ) gqa(ib) = 0.d0
      write(*,'(A,i3,A,i0,A,es10.3,A,f5.1,A,f6.4,A,f6.4,A,es10.3)')                     &
        '         bnd type', ib,                                                        &
        ' pts=',        nint(gn(ib)),                                                    &
        '  area=',      ga(ib),                                                          &
        ' active=',     1.d2*gaa(ib)/max(ga(ib),1.d-300),                                &
        ' %  qjac min=', gqm(ib),                                                        &
        ' (active ',    gqa(ib),                                                         &
        ')  |zj_sat| min=', gzm(ib)
    enddo

    ! --- The sign of xjac is an ORIENTATION CONVENTION: it flips with element_size_perp and so
    ! --- follows the side family. A difference between families is expected and meaningless. Only
    ! --- a MIXTURE within one (type, family) cell is a genuine fold, and only that is flagged.
    do ib = 1, nbt
      if ( gn(ib) .le. 0.d0 ) cycle
      do f = 1, 2
        if ( gsg(1,f,ib) .gt. 0.d0 .and. gsg(2,f,ib) .gt. 0.d0 ) then
          write(*,'(A,i3,A,i1,A,i0,A,i0,A)')                                             &
            '         *** WARNING bnd type', ib, ' side family ', f,                     &
            ': xjac changes sign within the family (', nint(gsg(1,f,ib)), ' pos, ',      &
            nint(gsg(2,f,ib)), ' neg) - INVERTED/FOLDED ELEMENT, not a convention'
        endif
      enddo
    enddo
  endif

  ! --- The two located points. Each rank keeps its own worst; the global max of the metric is
  ! --- shared, and exactly ONE rank contributes its payload while every other contributes zeros,
  ! --- so a SUM reduce delivers the winner without a point-to-point exchange.
  ! ---
  ! --- TIES ARE SYSTEMATIC, NOT ACCIDENTAL, and an earlier version assumed otherwise. Elements on
  ! --- an MPI partition boundary are processed by more than one rank, so the SAME Gauss point is
  ! --- recorded twice with bit-identical values - a guaranteed exact tie, not a 1-in-2^52 fluke.
  ! --- MEASURED on 1+4: two ranks tied and the summed payload printed R = 3.21, Z = -2.23 (outside
  ! --- the vessel), g_bn = 1.84 (the Chodura factor cannot exceed 1) and w_eff = 1.63 (a weight
  ! --- cannot exceed 1); every field halved to a valid value. So the winner is now chosen by
  ! --- LOWEST RANK among those holding the maximum, which is unique by construction.
  do ip = 1, 2
    if ( ip .eq. 1 ) then
      mx(1) = sg_wd(1) ; pay = sg_wd
    else
      mx(1) = sg_ws(1) ; pay = sg_ws
    endif
    call MPI_Allreduce(mx(1), gmx(1), 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    ! --- lowest rank holding the maximum wins; 1.d30 for every rank that does not hold it
    rnk(1) = 1.d30
    if ( mx(1) .eq. gmx(1) ) rnk(1) = dble(my_id)
    call MPI_Allreduce(rnk(1), grnk(1), 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
    if ( dble(my_id) .ne. grnk(1) ) pay = 0.d0
    call MPI_Reduce(pay, gpay, npay, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if ( my_id .ne. 0 ) cycle
    if ( gmx(1) .le. -1.d29 ) cycle
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
    ! --- b_n here is DIMENSIONLESS (bdotn as the assembly computes it). g_bn and Btot are
    ! --- recorded so zj_sat = c_sat*rho*(|g|cs + v_perp)/|B| can be decomposed with no assumption
    ! --- about how g depends on b_n - it is a tanh of |b_n|, not proportional to it.
    ! --- 16 edit descriptors, 16 output items.
    write(*,'(A,es9.2,A,es9.2,A,es9.2,A,es9.2,A,es9.2,A,es9.2,A,f6.4,A,es9.2)')          &
      '           b_n=',      gpay(9),  '  g_bn=',   gpay(13),                            &
      '  |B|=',               gpay(14), '  cs=',     gpay(10),                            &
      '  zj=',                gpay(11), '  zj_sat=', gpay(12),                            &
      '  w_eff=',             gpay(16), '  dS=',     gpay(15)
  enddo

  if ( my_id .ne. 0 ) return

  ! --- Histograms as a PERCENTAGE OF AREA, requested and active. No threshold is applied anywhere
  ! --- in this module on purpose: none can be justified before the distribution is known, and
  ! --- inventing one here would repeat the mistake that produced sheath_weak_wmin.
  write(*,'(A)') '         qjac histogram, % of area per 10 % bin (low = badly conditioned):'
  call sg_hist(gqh,  ga,  gn, 'req')
  call sg_hist(gqha, gaa, gn, 'act')
  write(*,'(A)') '         |zj_sat| histogram, % of area per decade, 1e-12 .. 1e0:'
  call sg_hist(gzh,  ga,  gn, 'req')
  call sg_hist(gzha, gaa, gn, 'act')

end subroutine sheath_geom_report


!> One histogram block, printed as a percentage of the relevant area so that the rows of different
!! boundary types are directly comparable and no reader has to renormalise by hand.
subroutine sg_hist(h, norm, cnt, tag)
  implicit none
  real*8,       intent(in) :: h(nbin,nbt), norm(nbt), cnt(nbt)
  character(3), intent(in) :: tag
  integer :: ib, k
  do ib = 1, nbt
    if ( cnt(ib) .le. 0.d0 .or. norm(ib) .le. 0.d0 ) cycle
    ! --- 3 edit descriptors before the loop, 3 items; then nbin pairs. Counted.
    write(*,'(A,i3,1x,A,A,12(1x,f6.2))') '           bnd type', ib, tag, ':',            &
      (1.d2*h(k,ib)/norm(ib), k=1,nbin)
  enddo
end subroutine sg_hist


end module mod_sheath_geom_diag
