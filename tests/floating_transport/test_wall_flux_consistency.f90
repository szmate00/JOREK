!> Locks the fluid-kinetic wall particle-flux convention.
!>
!> The fluid loses n_e*(v_par.n + v_ExB.n) through a boundary face. The kinetic
!> recycling source in mod_particle_wall_interaction must return the SAME flux,
!> or particles are destroyed where the drift pushes into the wall and created
!> where it pulls away. The two sides build that flux from different variables
!> in different bases, so this test recomputes both from one state and asserts
!> they agree:
!>
!>   fluid   (mod_expression):   vExB_norm = (-R*u_Z*nR + R*u_R*nZ)/fact_time
!>                               vparB_norm = vpar*(B.n)/fact_time
!>   kinetic (mod_fields):       v_ExB   = (/ -R*U_Z, R*U_R, 0 /)/t_norm
!>           (wall interaction): v_par_n = vpar*dot_product(B, n)
!>
!> The printed old/new column is the previous parallel-only formula divided by
!> the corrected one, i.e. the factor by which recycling was wrong.
program test_wall_flux_consistency
  implicit none
  integer :: nfail
  nfail = 0

  !                 R      |b.n|     v_par.n     vE.n        label
  call check_case( 1.60d0,  0.026d0,  3.12d3,   1.00d4, 'inner: drift INTO wall  ')
  call check_case( 1.59d0, -0.044d0,  3.17d4,  -2.90d4, 'outer: drift OFF wall   ')
  call check_case( 1.70d0,  0.030d0,  5.00d4,   0.00d0, 'no drift (unchanged)    ')
  call check_case( 1.50d0, -0.005d0,  0.00d0,  -2.00d3, 'grazing, inflow only    ')
  call check_case( 1.55d0,  0.020d0, -4.00d3,   1.00d3, 'parallel flow off wall  ')

  ! --- the normal/tangential split of both flow components
  call check_decomp( 0.9d0,  1.60d0,  0.70d0, -0.30d0,  1.9d0, -0.4d0, 4.8d4 )
  call check_decomp(-2.1d0,  1.26d0, -0.05d0,  0.12d0, -0.6d0,  2.2d0,-1.1d4 )
  call check_decomp( 0.0d0,  2.00d0,  0.00d0,  0.00d0,  1.0d0,  1.0d0, 1.0d0 )
  call check_decomp( 3.3d0,  1.05d0,  1.00d0,  0.00d0,  0.0d0,  1.0d0, 0.0d0 )

  if ( nfail == 0 ) then
    print *, 'PASS: fluid and kinetic wall normal flux agree (drift included, sign-correct)'
  else
    print *, 'FAIL: ', nfail, ' wall flux consistency check(s)'
    stop 1
  endif

contains

  !> vpn_target and vEn_target are the OUTWARD-POSITIVE normal components [m/s]
  subroutine check_case(R, bn_unit, vpn_target, vEn_target, label)
    real*8, intent(in) :: R, bn_unit, vpn_target, vEn_target
    character(*), intent(in) :: label
    real*8 :: t_norm, Btot, nR, nZ, Bnorm, B(3), nvec(3), v_ExB(3)
    real*8 :: U_R, U_Z, vpar, c, theta, tol, n_e
    real*8 :: fluid_vExB_n, fluid_vpar_n, fluid_total
    real*8 :: kin_vE_n, kin_vpar_n, kin_total, gamma_old, gamma_new

    t_norm = 3.1d-7 ; Btot = 2.5d0 ; n_e = 1.0d19 ; tol = 1.0d-9

    ! --- wall normal at an arbitrary orientation
    theta = 0.7d0 ; nR = cos(theta) ; nZ = sin(theta)
    nvec = (/ nR, nZ, 0.d0 /)

    ! --- B with the prescribed incidence b.n and |B| = Btot
    Bnorm = bn_unit * Btot
    B(1)  = Bnorm*nR - 0.3d0*Btot*nZ
    B(2)  = Bnorm*nZ + 0.3d0*Btot*nR
    B(3)  = sqrt( max( Btot**2 - Bnorm**2 - (0.3d0*Btot)**2, 0.d0 ) )

    ! --- pick vpar [m/s/T] so that vpar*(B.n) hits the requested v_par.n
    vpar = vpn_target / Bnorm

    ! --- pick u derivatives so that the drift hits the requested vE.n:
    ! --- with U_R = c*nZ, U_Z = -c*nR the normal drift is exactly c*R/t_norm
    c   = vEn_target * t_norm / R
    U_R =  c*nZ
    U_Z = -c*nR

    ! --- FLUID side, exactly as the boundary diagnostics compute it [m/s]
    fluid_vExB_n = ( -R*U_Z*nR + R*U_R*nZ ) / t_norm
    fluid_vpar_n = vpar * dot_product(B, nvec)
    fluid_total  = fluid_vpar_n + fluid_vExB_n

    ! --- KINETIC side, as mod_fields + mod_particle_wall_interaction now do it
    v_ExB      = (/ -R*U_Z, R*U_R, 0.d0 /) / t_norm
    kin_vE_n   = dot_product(v_ExB, nvec)
    kin_vpar_n = vpar * dot_product(B, nvec)
    kin_total  = kin_vpar_n + kin_vE_n

    ! --- exclude the c_angle floor: it is common to both formulas
    gamma_new = n_e * max(kin_total, 0.d0)
    gamma_old = n_e * abs(vpar) * Btot * abs(bn_unit)

    if ( abs(kin_vE_n   - vEn_target ) > tol*max(1.d0,abs(vEn_target )) .or.  &
         abs(kin_vpar_n - vpn_target ) > tol*max(1.d0,abs(vpn_target )) .or.  &
         abs(kin_vE_n   - fluid_vExB_n) > tol*max(1.d0,abs(fluid_vExB_n)) .or. &
         abs(kin_vpar_n - fluid_vpar_n) > tol*max(1.d0,abs(fluid_vpar_n)) .or. &
         abs(kin_total  - fluid_total ) > tol*max(1.d0,abs(fluid_total )) ) then
      nfail = nfail + 1
      print *, 'MISMATCH ', label
      print *, '   fluid   vE.n/vpar.n/total = ', fluid_vExB_n, fluid_vpar_n, fluid_total
      print *, '   kinetic vE.n/vpar.n/total = ', kin_vE_n, kin_vpar_n, kin_total
    endif

    print '(a,a,a,es10.2,a,es10.2,a,es10.2,a,f8.2)',                          &
      '   ', label, ' vpar.n=', kin_vpar_n, ' vE.n=', kin_vE_n,               &
      ' Gamma_new/n_e=', gamma_new/n_e, '   old/new=',                        &
      gamma_old / max(gamma_new, 1.d-30)
  end subroutine check_case

  !> vExB_norm/vExB_tan and vparB_norm/vpar_tan must be an ORTHOGONAL split of the
  !> same poloidal vector, so each pair has to satisfy
  !>     (v.n)^2 + (v.t)^2 = |v_pol|^2
  !> with t = (n_Z,-n_R). If either expression had a swapped or sign-flipped
  !> component the identity fails, which a single-component test cannot see.
  subroutine check_decomp(theta, R, U_R, U_Z, BR, BZ, vpar)
    real*8, intent(in) :: theta, R, U_R, U_Z, BR, BZ, vpar
    real*8 :: nR, nZ, t_norm, tol
    real*8 :: vE_n, vE_t, vE_pol2, vp_n, vp_t, vp_pol2

    t_norm = 3.1d-7 ; tol = 1.0d-12
    nR = cos(theta) ; nZ = sin(theta)

    ! --- exactly the expressions in mod_expression.f90
    vE_n = ( -R*U_Z*nR + R*U_R*nZ ) / t_norm            ! vExB_norm
    vE_t = -R * ( U_R*nR + U_Z*nZ ) / t_norm            ! vExB_tan
    vp_n = vpar * ( BR*nR + BZ*nZ ) / t_norm            ! vparB_norm (poloidal part)
    vp_t = vpar * ( BR*nZ - BZ*nR ) / t_norm            ! vpar_tan

    vE_pol2 = ( R/t_norm )**2 * ( U_R*U_R + U_Z*U_Z )
    vp_pol2 = ( vpar/t_norm )**2 * ( BR*BR + BZ*BZ )

    if ( abs(vE_n**2 + vE_t**2 - vE_pol2) > tol*max(1.d0,vE_pol2) .or. &
         abs(vp_n**2 + vp_t**2 - vp_pol2) > tol*max(1.d0,vp_pol2) ) then
      nfail = nfail + 1
      print *, 'MISMATCH decomposition at theta=', theta
      print *, '   ExB : ', vE_n**2 + vE_t**2, ' vs ', vE_pol2
      print *, '   par : ', vp_n**2 + vp_t**2, ' vs ', vp_pol2
    endif
  end subroutine check_decomp

end program test_wall_flux_consistency
