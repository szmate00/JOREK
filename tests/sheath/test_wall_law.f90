program test_wall_law
  use constants
  use phys_module
  use mod_sheath_bc
  implicit none
  real*8 :: an,c,vw,te,u,z,sat,x,j(7),q(7),qp(7),qm(7),jp(7),jm(7),zp,zm,sp,sm,h,fd
  real*8 :: g,dg,gn,dgn,bn,err,phi,reference,zero_current,nonzero_current
  real*8 :: T(2),Tb(2),Tp(2),Tm(2),vel,velb,dT(2),dbT(2),dbTb(2),vp,vm,vbp,vbm,dummy(2)
  integer :: i,sgn,mode,k
  do sgn=-1,1,2
    F0=sgn*2.972306d0
    call sheath_norm(an,c,vw)
    te=MU_ZERO*central_density*1.d20*EL_CHG ! independently: 1 eV
    u=2.d0*sheath_Lambda*te/an
    phi=F0*u/sqrt(MU_ZERO*central_density*1.d20*central_mass*ATOMIC_MASS_UNIT)
    call close(phi,3.d0,1.d-13,'1 eV -> +3 V for either F0')
    do k=-1,1,2
      q=[u,0.01d0,te,te,k*0.4d0,2.d0,k*0.01d0]
      call evaluate(q,z,j,sat)
      call close(z,0.d0,1.d-13*abs(sat),'floating root')
      q(1)=2.d0*100.d0*te/an
      call evaluate(q,z,j,sat)
      call close(z,sat,1.d-13*abs(sat),'ion saturation')
      call require(-z*(k*0.1d0)/(F0*MU_ZERO)>0.d0,'outward ion current sign')
      q(1)=0.d0
      call evaluate(q,z,j,sat)
      reference=z
      q(1)=-2.d0*100.d0*te/an
      call evaluate(q,z,j,sat)
      call close(z,reference,1.d-13*abs(sat),'electron saturation')
    enddo
  enddo
  F0=2.972306d0
  call sheath_norm(an,c,vw)
  sheath_jsat_from_vpar=.true.
  sheath_jsat_vpar_min=0.1d0
  q=[2.d0*4.d0*2.d-4/an,0.01d0,3.d-4,2.d-4,0.4d0,2.d0,0.d0]
  call evaluate(q,zero_current,j,sat)
  q(7)=1.d-100
  call evaluate(q,nonzero_current,j,sp)
  call close(sp,sat,1.d-14,'zero Vpar continuity with legacy floor')
  call close(zero_current,nonzero_current,1.d-14,'current continuity at zero Vpar')
  sheath_jsat_vpar_min=0.d0
  q(7)=0.d0
  call evaluate(q,z,j,sat)
  call close(sat,0.d0,0.d0,'no hidden ion source at stagnation')
  reference=z-sat
  q(7)=0.02d0
  call evaluate(q,z,j,sat)
  call close(z-sat,reference,1.d-13,'electron flux independent of Vpar')
  ! Decreasing outgoing ion flux shifts the floating drop by log(1/fraction).
  q(7)=0.5d0*q(5)*sqrt(GAMMA*(q(3)+q(4)))/q(6)
  q(1)=2.d0*(sheath_Lambda+log(2.d0))*q(4)/an
  call evaluate(q,z,j,sat)
  call close(z,0.d0,1.d-13,'floating drop responds to ion collection')

  ! Smooth regions on repelling and electron-saturated branches, both signs,
  ! with local Lambda, stagnant/incoming/outgoing/floored ion collection.
  err=0.d0
  do mode=1,6
    sheath_jsat_from_vpar=(mode/=1)
    sheath_Lambda_local=mod(mode,2)==0
    sheath_jsat_vpar_min=0.d0
    if (mode==4) sheath_jsat_vpar_min=0.1d0
    do k=-1,1,2
      q=[2.d0*4.d0*2.d-4/an,0.01d0,3.d-4,2.d-4,k*0.4d0,2.d0,k*0.01d0]
      if (mode==3 .or. mode==4) q(7)=-q(7)
      if (mode==5) q(1)=-q(1)
      if (mode==6) q(3)=1.d-15 ! hard floor: sensitivity must vanish
      call evaluate(q,z,j,sat)
      do i=1,7
        h=1.d-5*abs(q(i))
        qp=q; qm=q
        qp(i)=q(i)+h; qm(i)=q(i)-h
        call evaluate(qp,zp,jp,sp)
        call evaluate(qm,zm,jm,sm)
        fd=(zp-zm)/(2.d0*h)
        err=max(err,abs(fd-j(i))/max(abs(fd),abs(j(i)),1.d-9))
      enddo
    enddo
  enddo
  call require(err<2.d-6,'seven wall-law Jacobian columns')
  write(*,'(A,es12.4)') 'wall-law Jacobian maximum relative error: ',err

  do mode=1,2
    if (mode==2) vpar_smoothing_coef=[0.02d0,0.d0,0.d0]
    call sheath_incidence(0.d0,g,dg)
    call close(g,0.d0,0.d0,'exact tangency')
    do i=-10,10
      bn=dble(i)*0.006d0
      call sheath_incidence(bn,g,dg)
      call sheath_incidence(-bn,gn,dgn)
      call close(g,-gn,1.d-14,'odd incidence')
      h=1.d-7
      call sheath_incidence(bn+h,zp,sp)
      call sheath_incidence(bn-h,zm,sm)
      call close((zp-zm)/(2*h),dg,2.d-4,'incidence derivative through zero')
    enddo
  enddo

  ! Cold-state chain rule for the shared field-aligned entrance condition.
  T=[1.d-6,4.d-6]; Tb=[3.d-5,-2.d-5]
  call sheath_bohm_state(T,Tb,0.3d0,0.2d0,vel,velb,dT,dbT,dbTb)
  do i=1,2
    h=1.d-10; Tp=T; Tm=T; Tp(i)=Tp(i)+h; Tm(i)=Tm(i)-h
    call sheath_bohm_state(Tp,Tb,0.3d0,0.2d0,vp,vbp,dummy,jp(1:2),jp(3:4))
    call sheath_bohm_state(Tm,Tb,0.3d0,0.2d0,vm,vbm,dummy,jm(1:2),jm(3:4))
    call close((vp-vm)/(2*h),dT(i),1.d-5,'Bohm temperature derivative')
    call close((vbp-vbm)/(2*h),dbT(i),1.d-4,'Bohm tangent-temperature derivative')
    Tp=Tb; Tm=Tb; Tp(i)=Tp(i)+h; Tm(i)=Tm(i)-h
    call sheath_bohm_state(T,Tp,0.3d0,0.2d0,vp,vbp,dummy,jp(1:2),jp(3:4))
    call sheath_bohm_state(T,Tm,0.3d0,0.2d0,vm,vbm,dummy,jm(1:2),jm(3:4))
    call close((vbp-vbm)/(2*h),dbTb(i),1.d-5,'Bohm tangent-gradient derivative')
  enddo
  call sheath_bohm_state(T,Tb,0.d0,0.d0,vel,velb,dT,dbT,dbTb)
  call close(vel,0.d0,0.d0,'finite tangent entrance limit')
  T_min_sheath=0.d0
  call require(sheath_temp_floor(-1.d-5)>0.d0,'zero-width floor stays finite')
  call close(dsheath_temp_floor_dT(-1.d-5),0.d0,0.d0,'zero-width floor derivative')
  print *, 'PASS: wall law, normalization, incidence and cold entrance-state tests'
contains
  subroutine evaluate(q,z,j,sat)
    real*8,intent(in)::q(7)
    real*8,intent(out)::z,j(7),sat
    real*8::x
    call sheath_current(q(1),q(2),q(3),q(4),q(5),sign(1.d0,q(5)),q(6), &
      z,j(1),j(2),j(3),j(4),sat,x,vpar=q(7),dzj_dvpar=j(7),dzj_dg=j(5),dzj_dB=j(6))
  end subroutine
  subroutine require(condition,name)
    logical,intent(in)::condition
    character(*),intent(in)::name
    if (.not.condition) then
      print *, 'FAIL: ',name
      stop 1
    endif
  end subroutine
  subroutine close(a,b,tol,name)
    real*8,intent(in)::a,b,tol
    character(*),intent(in)::name
    if (abs(a-b)>tol) then
      print *, 'FAIL: ',name,a,b,tol
      stop 1
    endif
  end subroutine
end program
