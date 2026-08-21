program test_sheath_bc
  use constants,     only: PI, MU_ZERO, ATOMIC_MASS_UNIT, EL_CHG, MASS_ELECTRON
  use phys_module
  use mod_sheath_bc
  implicit none

  real*8 :: rho0, m_i, tnorm, phifac, Te, Ti, u, zj, du, dr, dti, dte, zjs, x
  real*8 :: sheath_a_n, sheath_c_sat, sheath_vw
  real*8 :: zj2, d2(5), zjp, zjm, fd, lam, dl1, dl2, xl, dxl, xl2, dxl2
  integer :: i, nfail
  real*8, parameter :: TeV = 20.d0

  nfail = 0

  ! --- reference case: AUG-like, D plasma, n0 = 1e20
  F0 = 4.1d0; GAMMA = 5.d0/3.d0; central_density = 1.d0; central_mass = 2.d0
  sheath_V_wall = 0.d0; sheath_Lambda = -1.d0; sheath_Lambda_local = .true.
  sheath_X_min = -1.d3; sheath_smooth_dX = 0.5d0; sheath_min_bn = 0.d0   ! limiter effectively off
  sheath_sat_slope = 0.d0                                               ! exact saturation

  m_i    = central_mass * ATOMIC_MASS_UNIT
  rho0   = central_density * 1.d20 * m_i
  tnorm  = EL_CHG * MU_ZERO * central_density * 1.d20      ! T_jorek = tnorm * T_eV
  phifac = -F0 / sqrt(MU_ZERO*rho0)                        ! Phi[V] = phifac * u
  Te = tnorm*TeV; Ti = tnorm*TeV

  call sheath_norm(sheath_a_n, sheath_c_sat, sheath_vw)
  write(*,'(A,es12.4,A,es12.4)') ' a_n =', sheath_a_n, '   c_sat =', sheath_c_sat
  call sheath_get_lambda(Ti, Te, lam, dl1, dl2)
  write(*,'(A,f8.4,A,f8.4)') ' Lambda =', lam, '   Lambda0 - 0.5*ln(2*gamma) =', &
        log(sqrt(m_i/(2.d0*PI*MASS_ELECTRON))) - 0.5d0*log(2.d0*GAMMA)
  call check('Lambda(Ti=Te)', lam, log(sqrt(m_i/(2.d0*PI*MASS_ELECTRON)))-0.5d0*log(2.d0*GAMMA), 1.d-12, nfail)

  ! --- A: zero current at Phi = Lambda*Te/e
  u = 2.d0*(Te*lam + sheath_vw)/sheath_a_n
  call sheath_current(u, 0.05d0, Ti, Te, 1.d0, 1.d0, 2.5d0, zj, du, dr, dti, dte, zjs, x)
  call check('A: j at floating potential', zj, 0.d0, 1.d-20, nfail)
  call check('A: Phi_float [V]', phifac*u, lam*TeV, 1.d-10, nfail)
  write(*,'(A,f10.4,A)') ' floating potential =', phifac*u, ' V'

  ! --- F/B: strongly ion-saturated (Phi >> Lambda Te/e) -> zj -> zj_sat > 0 for g_bn = +1
  u = 2.d0*(Te*(lam+20.d0))/sheath_a_n
  call sheath_current(u, 0.05d0, Ti, Te, 1.d0, 1.d0, 2.5d0, zj, du, dr, dti, dte, zjs, x)
  call check('B: j -> j_sat', zj/zjs, 1.d0, 1.d-8, nfail)
  zj2 = du                                              ! dj/du deep in ion saturation
  u = 2.d0*(Te*lam)/sheath_a_n                          ! ... versus at the floating potential
  call sheath_current(u, 0.05d0, Ti, Te, 1.d0, 1.d0, 2.5d0, zjp, du, dr, dti, dte, zjs, x)
  call check('B: dj/du suppressed by exp(-X)', abs(zj2/du), 0.d0, 1.d-8, nfail)
  u = 2.d0*(Te*(lam+20.d0))/sheath_a_n
  if (zjs .le. 0.d0) then
    write(*,*) 'FAIL  F: zj_sat must be > 0 for g_bn=+1 and F0>0, got', zjs; nfail = nfail+1
  else
    write(*,'(A,es12.4)') ' PASS  F: zj_sat > 0 for outward field :', zjs
  endif

  ! --- C: electron saturation at Phi = 0
  call sheath_current(0.d0, 0.05d0, Ti, Te, 1.d0, 1.d0, 2.5d0, zj, du, dr, dti, dte, zjs, x)
  call check('C: j(Phi=0)/j_sat', zj/zjs, 1.d0-exp(lam), 1.d-10, nfail)

  ! --- D: analytic vs finite-difference derivatives, three regimes
  do i = 1, 3
    if (i==1) u = 2.d0*(Te*(lam+0.3d0))/sheath_a_n      ! unsaturated
    if (i==2) u = 2.d0*(Te*(lam+8.d0 ))/sheath_a_n      ! ion saturated
    if (i==3) u = 2.d0*(Te*(lam-2.d0 ))/sheath_a_n      ! electron side
    call sheath_current(u, 0.05d0, Ti, Te, 1.d0, 1.d0, 2.5d0, zj, du, dr, dti, dte, zjs, x)
    d2 = (/ du, dr, dti, dte, 0.d0 /)
    call fdcheck(i, 1, u,      0.05d0, Ti, Te, d2(1), nfail)
    call fdcheck(i, 2, u,      0.05d0, Ti, Te, d2(2), nfail)
    call fdcheck(i, 3, u,      0.05d0, Ti, Te, d2(3), nfail)
    call fdcheck(i, 4, u,      0.05d0, Ti, Te, d2(4), nfail)
  enddo

  ! --- D2: same, with the finite saturation conductance switched on. The slope term is added to
  ! --- BOTH f and df/dX, so an inconsistency there would show up here and nowhere else.
  sheath_sat_slope = 0.05d0
  do i = 1, 3
    if (i==1) u = 2.d0*(Te*(lam+0.3d0))/sheath_a_n
    if (i==2) u = 2.d0*(Te*(lam+8.d0 ))/sheath_a_n
    if (i==3) u = 2.d0*(Te*(lam-2.d0 ))/sheath_a_n
    call sheath_current(u, 0.05d0, Ti, Te, 1.d0, 1.d0, 2.5d0, zj, du, dr, dti, dte, zjs, x)
    d2 = (/ du, dr, dti, dte, 0.d0 /)
    call fdcheck(10+i, 1, u,   0.05d0, Ti, Te, d2(1), nfail)
    call fdcheck(10+i, 2, u,   0.05d0, Ti, Te, d2(2), nfail)
    call fdcheck(10+i, 3, u,   0.05d0, Ti, Te, d2(3), nfail)
    call fdcheck(10+i, 4, u,   0.05d0, Ti, Te, d2(4), nfail)
  enddo

  ! --- and that the characteristic is now unbounded above, so any demanded current is reachable
  u = 2.d0*(Te*(lam+20.d0))/sheath_a_n
  call sheath_current(u, 0.05d0, Ti, Te, 1.d0, 1.d0, 2.5d0, zj, du, dr, dti, dte, zjs, x)
  if ( zj/zjs .gt. 1.5d0 ) then
    write(*,'(A,es10.2)') ' PASS  D2: f exceeds 1 at large X, ratio ', zj/zjs
  else
    write(*,'(A,es10.2)') ' FAIL  D2: f still capped at 1, ratio ', zj/zjs
    nfail = nfail + 1
  endif
  sheath_sat_slope = 0.d0

  ! --- E: limiter is C1, monotone, and safe far outside the range
  sheath_X_min = -3.d0; sheath_smooth_dX = 0.5d0
  call sheath_norm(sheath_a_n, sheath_c_sat, sheath_vw)
  do i = -400, 400
    x = dble(i)*0.5d0
    call sheath_x_limited(x, xl, dxl)
    if (xl .ne. xl .or. abs(xl) .gt. 1.d30) then
      write(*,*) 'FAIL  E: limiter not finite at X =', x; nfail = nfail+1; exit
    endif
    if (dxl .lt. -1.d-14 .or. dxl .gt. 1.d0+1.d-14) then
      write(*,*) 'FAIL  E: dXlim/dX out of [0,1] at X =', x, dxl; nfail = nfail+1; exit
    endif
    call sheath_x_limited(x+1.d-6, xl2, dxl2)
    if (xl2 .lt. xl - 1.d-14) then
      write(*,*) 'FAIL  E: limiter not monotone at X =', x; nfail = nfail+1; exit
    endif
    fd = (xl2-xl)/1.d-6
    if (abs(fd-dxl) .gt. 1.d-4) then
      write(*,*) 'FAIL  E: limiter derivative wrong at X =', x, fd, dxl; nfail = nfail+1; exit
    endif
  enddo
  if (nfail .eq. 0) write(*,*) 'PASS  E: limiter C1, monotone and finite over X in [-200,200]'

  ! --- E2: full characteristic must stay finite at extreme potentials
  do i = -20, 20
    u = 2.d0*(Te*(lam+dble(i)*10.d0))/sheath_a_n
    call sheath_current(u, 0.05d0, Ti, Te, 1.d0, 1.d0, 2.5d0, zj, du, dr, dti, dte, zjs, x)
    if (zj .ne. zj .or. abs(zj) .gt. 1.d30 .or. du .ne. du .or. abs(du) .gt. 1.d30) then
      write(*,*) 'FAIL  E2: not finite at X =', dble(i)*10.d0, zj, du; nfail = nfail+1; exit
    endif
  enddo
  if (nfail .eq. 0) write(*,*) 'PASS  E2: characteristic finite for X in [-200,200]'

  write(*,*)
  if (nfail .eq. 0) then
    write(*,*) '=== ALL SHEATH UNIT TESTS PASSED ==='
  else
    write(*,*) '=== ', nfail, ' FAILURES ==='
    stop 1
  endif

contains

  subroutine check(name, got, want, tol, nf)
    character(len=*), intent(in) :: name
    real*8, intent(in) :: got, want, tol
    integer, intent(inout) :: nf
    real*8 :: err
    err = abs(got-want) / max(abs(want), 1.d-30)
    if (abs(want) .lt. 1.d-25) err = abs(got-want)
    if (err .gt. tol) then
      write(*,'(A,A,A,es14.6,A,es14.6,A,es10.2)') ' FAIL  ', name, ' got', got, ' want', want, ' err', err
      nf = nf + 1
    else
      write(*,'(A,A,A,es10.2)') ' PASS  ', name, '   err', err
    endif
  end subroutine

  subroutine fdcheck(ireg, ivar, u0, r0, Ti0, Te0, dana, nf)
    integer, intent(in) :: ireg, ivar
    real*8,  intent(in) :: u0, r0, Ti0, Te0, dana
    integer, intent(inout) :: nf
    real*8 :: h, a(4), ap(4), am(4), zp, zm, dum(5), fd, err
    character(len=4), parameter :: nm(4) = (/ 'u   ','rho ','Ti  ','Te  ' /)
    a = (/ u0, r0, Ti0, Te0 /)
    h = max(abs(a(ivar)), 1.d-30) * 1.d-6
    ap = a; ap(ivar) = a(ivar) + h
    am = a; am(ivar) = a(ivar) - h
    call sheath_current(ap(1), ap(2), ap(3), ap(4), 1.d0, 1.d0, 2.5d0, zp, dum(1),dum(2),dum(3),dum(4), dum(5), dum(5))
    call sheath_current(am(1), am(2), am(3), am(4), 1.d0, 1.d0, 2.5d0, zm, dum(1),dum(2),dum(3),dum(4), dum(5), dum(5))
    fd  = (zp-zm)/(2.d0*h)
    err = abs(fd-dana)/max(abs(fd),abs(dana),1.d-30)
    if (err .gt. 1.d-6) then
      write(*,'(A,I2,A,A,A,es14.6,A,es14.6,A,es10.2)') ' FAIL  D regime',ireg,' d/d',nm(ivar), &
            ' analytic', dana, ' fd', fd, ' err', err
      nf = nf + 1
    else
      write(*,'(A,I2,A,A,A,es10.2)') ' PASS  D regime',ireg,' d/d',nm(ivar),'   err', err
    endif
  end subroutine

end program test_sheath_bc
