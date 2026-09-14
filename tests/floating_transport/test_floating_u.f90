!> Normalisation of the floating-potential row: self-test for both signs of F0 and a nonzero wall
!! bias reconstructed to the volt.
program test_floating_u
  use phys_module,   only: F0, sheath_Lambda, sheath_V_wall
  use constants,     only: EL_CHG, MU_ZERO
  use mod_floating_u
  implicit none
  real*8 :: a_n, C_T, C_V, Te_1eV, u, volts

  if ( .not. floating_u_selftest(0) ) error stop 'FAIL: floating_u selftest, F0 > 0'
  F0 = -F0
  if ( .not. floating_u_selftest(0) ) error stop 'FAIL: floating_u selftest, F0 < 0'
  F0 = -F0

  ! --- Te = 2 eV with a 5 V wall must reconstruct to 2*Lambda + 5 volts (WITH_TiTe fixture)
  sheath_V_wall = 5.d0
  call floating_u_norm(a_n, C_T, C_V)
  Te_1eV = EL_CHG * MU_ZERO * 1.d20
  u      = C_T * 2.d0*Te_1eV + C_V * sheath_V_wall
  volts  = floating_u_volts(u)
  if ( abs(volts - (2.d0*sheath_Lambda + 5.d0)) > 1.d-10 ) then
    write(*,*) 'FAIL: wall bias reconstruction', volts
    error stop 1
  endif
  write(*,'(a)') ' PASS: floating_u normalisation, both field signs and wall bias'
end program
