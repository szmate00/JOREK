! Experimental, independently switchable transport closures for prescribed floating u.
! No state, MPI, floors on evolved fields, or division by the normal magnetic field.
module mod_floating_transport
  implicit none
  private
  public :: floating_mach_flux, floating_wall_flux, density_transport_diffusion
  public :: floating_temperature_slope
contains
  pure real*8 function floating_temperature_slope(raw, knee, coef) result(slope)
    real*8, intent(in) :: raw, knee, coef(2)
    slope = 1.d0
    if (raw < knee*sum(coef)) slope = exp((raw-knee*sum(coef))/(knee*coef(2)))
  end function
  ! Weighted normal-flow residual. Magnetic geometry is Picard-lagged by the caller.
  ! At finite incidence its zero is vn = |b.n| f cs. At exact tangency it contributes
  ! no parallel-momentum constraint: the bulk equation remains, not Vpar=0.
  pure subroutine floating_mach_flux(vpar, ven, bn, bmag, cs, smoothing, coef, residual, jac)
    real*8, intent(in) :: vpar, ven, bn, bmag, cs, coef(3)
    logical, intent(in) :: smoothing
    real*8, intent(out) :: residual, jac(3) ! derivatives wrt vpar, ven, cs
    real*8 :: a, f, target
    a = bn / bmag
    f = 1.d0
    if (smoothing) f = max(0.d0, 0.25d0*(1.d0+tanh((abs(a)-coef(1))/coef(2)))**2-coef(3))
    target = abs(a)*f
    residual = a*(bn*vpar + ven - target*cs)
    jac = (/ a*bn, a, -a*target /)
  end subroutine

  ! Absorbing charged-particle wall: no net plasma injection. Returned quantities
  ! are additional outward fluxes relative to the strong volume advection.
  ! gamma_here is the existing JOREK temperature-equation transmission coefficient,
  ! NOT the Stangeby gamma. Geometry/velocity and cs derivatives are supplied below.
  pure subroutine floating_wall_flux(vn, cs, angle, gamma_here, particle, heat, dp, dh)
    real*8, intent(in) :: vn, cs, angle, gamma_here
    real*8, intent(out) :: particle, heat, dp(2), dh(2) ! wrt vn, cs
    real*8 :: collect, outward
    outward = 0.d0
    if (vn > 0.d0) outward = 1.d0
    collect = max(vn,0.d0) + angle*cs
    particle = collect-vn
    ! Strong pressure advection already transports rho*T*vn. Cancel incoming
    ! advection as for particles, retaining the legacy minimum-angle heat term.
    heat = (gamma_here-1.d0)*collect + max(vn,0.d0)-vn
    dp = (/ outward-1.d0, angle /)
    dh = (/ gamma_here*outward-1.d0, (gamma_here-1.d0)*angle /)
  end subroutine

  ! Density-gradient sensor independent of pressure. Isotropic conservative
  ! diffusion, frozen within a linear solve. No assertion of discrete positivity.
  ! s,t span [0,1]; gradients of these coordinates give a flow-direction cell size.
  pure real*8 function density_transport_diffusion(rho, grad_rho, velocity, grad_s, grad_t) result(diff)
    real*8, intent(in) :: rho, grad_rho(2), velocity(2), grad_s(2), grad_t(2)
    real*8 :: speed, rate, h, variation, sensor
    diff = 0.d0
    speed = norm2(velocity)
    rate = abs(dot_product(velocity,grad_s)) + abs(dot_product(velocity,grad_t))
    if (speed <= tiny(speed) .or. rate <= tiny(rate)) return
    h = speed/rate
    variation = h*norm2(grad_rho)
    if (variation <= tiny(variation)) return
    sensor = variation/(abs(rho)+variation)
    diff = 0.5d0*h*speed*sensor
  end function
end module
