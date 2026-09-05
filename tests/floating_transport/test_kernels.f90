program test_kernels
  use mod_floating_transport
  implicit none
  real*8 :: x(3),xp(3),xm(3),coef(3),r,j(3),rp,rm,junk(3),eps,bn,bmag,cs
  real*8 :: p,h,dp(2),dh(2),pp,hp,pm,hm,vn,d1,d2,vel(2),gs(2),gt(2),grad(2)
  real*8 :: junk2(2)
  real*8 :: Rgeo,ur,uz,tangent(2),normal(2),ven,psir,psiz,ub,psib,vp
  integer :: k,sgn,i
  coef=(/0.02d0,0.016d0,0.005754d0/)
  eps=1.d-6
  do sgn=-1,1,2
    bn=sgn*0.13d0; bmag=1.9d0
    x=(/0.22d0,-0.03d0,0.08d0/)
    call floating_mach_flux(x(1),x(2),bn,bmag,x(3),.true.,coef,r,j)
    do k=1,3
      xp=x; xm=x; xp(k)=xp(k)+eps; xm(k)=xm(k)-eps
      call floating_mach_flux(xp(1),xp(2),bn,bmag,xp(3),.true.,coef,rp,junk)
      call floating_mach_flux(xm(1),xm(2),bn,bmag,xm(3),.true.,coef,rm,junk)
      call close(j(k),(rp-rm)/(2*eps),'Mach Jacobian')
    enddo
    ! With f=1, exact drift cancellation has NO extra factor/Bmag.
    vp=(abs(bn/bmag)*x(3)-x(2))/bn
    call floating_mach_flux(vp,x(2),bn,bmag,x(3),.false.,coef,r,j)
    call close(r,0.d0,'finite incidence normal Mach target')
  enddo
  call floating_mach_flux(100.d0,-0.1d0,0.d0,2.d0,0.1d0,.true.,coef,r,j)
  call close(r,0.d0,'tangent residual')
  if (any(j/=0.d0)) error stop 'tangent parallel row must have no constraint'
  do i=-2,2
    vn=i*0.15d0; cs=0.08d0
    call floating_wall_flux(vn,cs,0.01d0,3.d0,p,h,dp,dh)
    call close(vn+p,max(vn,0.d0)+0.01d0*cs,'particle flux difference: no duplicate advection')
    call close(vn+h,3.d0*max(vn,0.d0)+2.d0*0.01d0*cs,'pressure flux difference')
    if (vn+p<0.d0) error stop 'absorbing wall injected particles'
    if (i==0) cycle ! semismooth switch; one-sided derivative selected at zero
    call floating_wall_flux(vn+eps,cs,0.01d0,3.d0,pp,hp,junk(1:2),junk2)
    call floating_wall_flux(vn-eps,cs,0.01d0,3.d0,pm,hm,junk(1:2),junk2)
    call close(dp(1),(pp-pm)/(2*eps),'particle velocity Jacobian')
    call close(dh(1),(hp-hm)/(2*eps),'pressure velocity Jacobian')
    call floating_wall_flux(vn,cs+eps,0.01d0,3.d0,pp,hp,junk(1:2),junk2)
    call floating_wall_flux(vn,cs-eps,0.01d0,3.d0,pm,hm,junk(1:2),junk2)
    call close(dp(2),(pp-pm)/(2*eps),'particle sound-speed Jacobian')
    call close(dh(2),(hp-hm)/(2*eps),'pressure sound-speed Jacobian')
  enddo
  ! Sign/orientation manufactured check directly from JOREK velocity components.
  Rgeo=1.6d0; ur=0.03d0; uz=-0.02d0; psir=0.06d0; psiz=0.07d0
  tangent=(/0.6d0,0.8d0/); normal=(/tangent(2),-tangent(1)/)
  ven=dot_product((/-Rgeo*uz,Rgeo*ur/),normal)
  ub=ur*tangent(1)+uz*tangent(2); psib=psir*tangent(1)+psiz*tangent(2)
  bn=dot_product((/psiz,-psir/),normal)/Rgeo
  call close(ven,-Rgeo*ub,'signed ExB identity')
  call close(-ven/bn,Rgeo**2*ub/psib,'drift cancellation identity')
  ! Density sensor responds without any pressure input, vanishes for uniform rho,
  ! is nonnegative and invariant under rotation and coordinate reversal.
  gs=(/2.d0,0.3d0/); gt=(/0.2d0,3.d0/); vel=(/0.4d0,-0.2d0/); grad=(/1.d0,2.d0/)
  d1=density_transport_diffusion(0.2d0,grad,vel,gs,gt)
  if (d1<=0.d0) error stop 'density sensor inactive'
  d2=density_transport_diffusion(0.2d0,rot(grad),rot(vel),rot(gs),rot(gt))
  call close(d1,d2,'rotated density diffusion')
  d2=density_transport_diffusion(0.2d0,grad,vel,-gs,gt)
  call close(d1,d2,'reversed coordinate diffusion')
  call close(density_transport_diffusion(0.2d0,(/0.d0,0.d0/),vel,gs,gt),0.d0,'constant state')
  call close(density_transport_diffusion(0.d0,grad,(/0.d0,0.d0/),gs,gt),0.d0,'zero velocity')
  d2=density_transport_diffusion(-0.01d0,grad,vel,gs,gt)
  if (d2<0.d0) error stop 'negative density produced antidiffusion'
  call close(floating_temperature_slope(0.1d0,0.03d0,(/0.5d0,0.5d0/)),1.d0,'warm temperature slope')
  write(*,*) 'PASS: production floating transport kernels, Jacobians, flux balances and orientation'
contains
  subroutine close(a,b,name)
    real*8, intent(in) :: a,b
    character(*), intent(in) :: name
    if (abs(a-b)>2.d-9*max(1.d0,abs(a),abs(b))) then
      write(*,*) 'FAIL ',name,a,b
      error stop 1
    endif
  end subroutine
  pure function rot(a) result(b)
    real*8,intent(in)::a(2)
    real*8::b(2)
    b=(/-a(2),a(1)/)
  end function
end program
