!> The Mach1 slope row's u columns must BE the derivatives of its residual.
!>
!> mod_boundary_conditions.f90 imposes, on the bicubic tangential-derivative row,
!>     dMach1BC += m1_dr * m1_dsdx * m1_dDdb,   m1_dDdb = C1*U0_b + C2*u0_bb_r
!> where U0_b is this node's own u tangential-derivative DOF and u0_bb_r is the
!> reconstructed second tangential derivative, built from the value AND derivative
!> DOFs of BOTH edge endpoints. So the residual depends on FOUR u DOFs.
!>
!> For a long time it depended on them with no Jacobian columns at all: the only
!> u column ever written for that row was gated behind
!>     include_2nd_derivatives .and. n_order >= 5
!> with include_2nd_derivatives a parameter = .false. and n_order = 3 in production.
!> With one linearised solve per timestep and no Newton iteration, a missing column
!> is a systematic per-step error, and these columns scale as 1/psi_b - worst
!> exactly at grazing incidence. This test finite-differences the residual with
!> respect to each of the four DOFs and asserts the coefficients match.
program test_mach_slope_jacobian
  implicit none
  integer :: nfail
  nfail = 0
  !              BigR    R_b     ps0_b    ps0_bb   e0      label
  call check( 1.60d0, 0.31d0,  4.0d-2, -1.1d-2, 0.55d0, 'typical wall node    ')
  call check( 1.26d0,-0.72d0,  3.0d-3,  2.4d-2, 0.48d0, 'near-grazing, 1/psi_b')
  call check( 2.05d0, 0.04d0, -6.0d-2,  0.0d0 , 0.61d0, 'sign-flipped psi_b   ')
  call check( 1.45d0, 1.10d0,  1.2d-1, -8.0d-2, 0.39d0, 'strong psi curvature ')
  if ( nfail == 0 ) then
    print *, 'PASS: Mach1 slope-row u columns are the exact residual derivatives'
  else
    print *, 'FAIL: ', nfail, ' slope-row Jacobian check(s)'
    stop 1
  endif

contains

  !> m1_dDdb exactly as the production line forms it.
  pure function dDdb(BigR, R_b, ps0_b, ps0_bb, e0, s1v, s1d, s2v, s2d, Hss, u) result(r)
    real*8, intent(in) :: BigR, R_b, ps0_b, ps0_bb, e0, s1v, s1d, s2v, s2d
    real*8, intent(in) :: Hss(2,2), u(4)
    real*8 :: r, U0_b, u0_bb_r, C1, C2
    U0_b    = e0 * u(2)                                   ! this node's deriv DOF
    u0_bb_r = s1v*u(1)*Hss(1,1) + s1d*u(2)*Hss(1,2)   &   ! this node value/deriv
            + s2v*u(3)*Hss(2,1) + s2d*u(4)*Hss(2,2)       ! neighbour value/deriv
    C1 = ( 2.d0*BigR*R_b - BigR**2 * ps0_bb / ps0_b ) / ps0_b
    C2 =   BigR**2 / ps0_b
    r  = C1 * U0_b + C2 * u0_bb_r
  end function dDdb

  subroutine check(BigR, R_b, ps0_b, ps0_bb, e0, label)
    real*8, intent(in) :: BigR, R_b, ps0_b, ps0_bb, e0
    character(*), intent(in) :: label
    real*8 :: s1v, s1d, s2v, s2d, Hss(2,2), u(4), up(4), um(4)
    real*8 :: C1, C2, an(4), fd(4), h, tol
    integer :: k

    ! --- element sizes and Hermite second-derivative weights: arbitrary but fixed
    s1v = 0.83d0 ; s1d = 0.47d0 ; s2v = 1.21d0 ; s2d = 0.66d0
    Hss(1,1) = -1.7d0 ; Hss(1,2) = -0.9d0
    Hss(2,1) =  1.7d0 ; Hss(2,2) = -0.5d0
    u = (/ 0.21d0, -0.34d0, 0.08d0, 0.52d0 /)
    h = 1.d-6 ; tol = 1.d-7

    ! --- the four analytic columns, exactly as the production lines form them
    C1 = ( 2.d0*BigR*R_b - BigR**2 * ps0_bb / ps0_b ) / ps0_b
    C2 =   BigR**2 / ps0_b
    an(1) = C2 * s1v * Hss(1,1)                        ! dMach1BC_uv
    an(2) = C1 * e0  + C2 * s1d * Hss(1,2)             ! dMach1BC_ud
    an(3) = C2 * s2v * Hss(2,1)                        ! dMach1BC_unv
    an(4) = C2 * s2d * Hss(2,2)                        ! dMach1BC_und

    do k = 1, 4
      up = u ; um = u
      up(k) = u(k) + h ; um(k) = u(k) - h
      fd(k) = ( dDdb(BigR,R_b,ps0_b,ps0_bb,e0,s1v,s1d,s2v,s2d,Hss,up)    &
              - dDdb(BigR,R_b,ps0_b,ps0_bb,e0,s1v,s1d,s2v,s2d,Hss,um) ) / (2.d0*h)
      if ( abs(fd(k)-an(k)) > tol*max(1.d0,abs(an(k))) ) then
        nfail = nfail + 1
        print *, 'MISMATCH ', label, ' dof ', k, ' fd=', fd(k), ' analytic=', an(k)
      endif
    enddo

    ! --- the C1 term must actually be present in the deriv-DOF column: dropping it
    ! --- (the historical bug was omitting the whole set) must be detectable
    if ( abs(an(2) - C2*s1d*Hss(1,2)) < 1.d-12 ) then
      nfail = nfail + 1
      print *, 'DEGENERATE ', label, ' C1 term vanished from the deriv column'
    endif

    print '(a,a,a,es10.2,a,es10.2,a,es10.2)', '   ', label,                &
      ' C1=', C1, ' C2=', C2, ' max|fd-an|=', maxval(abs(fd-an))
  end subroutine check

end program test_mach_slope_jacobian
