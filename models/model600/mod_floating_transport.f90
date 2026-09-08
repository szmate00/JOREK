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
  ! Weighted normal-flow residual for the drift-compatible Bohm condition, in the
  ! NON-MARGINAL form (SOLPS BCMOM=13 with MOMPAR(,,2) >= 0.5, the branch the manual
  ! marks "recommended for cases with drifts"). The imposed target is on the PARALLEL
  ! normal flux:
  !
  !     Bn*Vpar  =  |a|*f*cs  +  min( max(-ven,0), 2*cs*|a| )
  !
  ! i.e. sonic outflow, plus a one-sided supplement that compensates INWARD ExB only,
  ! bounded by SOLPS's +-2*cs*|b_x| (b2stbc_cbc = 1.0). Outward drift leaves the target
  ! at plain sonic - the Bohm condition is an inequality and extra outward flux is
  ! allowed - so the outward limit is exactly the configuration measured stable.
  !
  ! WHY THIS FORM AND NOT THE PREVIOUS ONE. The earlier residual was
  ! a*(bn*vpar + ven - |a|*f*cs), whose zero is total normal flow = sonic exactly:
  ! the MARGINAL branch, which permits and at strong outward drift demands REVERSED
  ! parallel flow. Imposed nodally that crashed at 466 steps. It also carried the RAW
  ! ven with no bound at all, so at the measured wall values it demanded ~10*cs -
  ! worse than the nodal route it was meant to replace.
  !
  ! Everything here is bounded and division-free: the supplement is formed as a FLUX,
  ! so unlike the nodal row there is no 1/bn inversion and hence NO incidence floor is
  ! needed. The overall weight `a` makes the constraint degenerate as a^2 at tangency,
  ! where it contributes no parallel-momentum constraint at all and the bulk equation
  ! remains - rather than dividing by a quantity that is going to zero.
  !
  ! Bn*Vpar is bounded into [|a|*f*cs, 3*cs*|a|], i.e. v_par in [f*cs, 3*cs], the same
  ! envelope as the nodal row. Branches are piecewise linear, so a branch frozen for
  ! one linear solve is exact.
  pure subroutine floating_mach_flux(vpar, ven, bn, bmag, cs, smoothing, coef, residual, jac)
    real*8, intent(in) :: vpar, ven, bn, bmag, cs, coef(3)
    logical, intent(in) :: smoothing
    real*8, intent(out) :: residual, jac(3) ! derivatives wrt vpar, ven, cs
    real*8 :: a, f, bound, sup, w_act, w_clip
    a = bn / bmag
    f = 1.d0
    if (smoothing) f = max(0.d0, 0.25d0*(1.d0+tanh((abs(a)-coef(1))/coef(2)))**2-coef(3))
    ! One-sided, clipped supplement. w_act selects the unclipped inward branch,
    ! w_clip the saturated one; both zero for outward drift.
    bound = 2.d0*cs*abs(a)
    sup   = 0.d0 ; w_act = 0.d0 ; w_clip = 0.d0
    if (-ven .ge. bound) then
      sup = bound ; w_clip = 1.d0
    elseif (-ven .gt. 0.d0) then
      sup = -ven  ; w_act  = 1.d0
    endif
    residual = a*( bn*vpar - abs(a)*f*cs - sup )
    jac = (/ a*bn, a*w_act, -a*abs(a)*( f + 2.d0*w_clip ) /)
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
