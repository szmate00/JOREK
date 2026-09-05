! Diagnostics only. Extrema are idempotent under halo replication; no area sums.
! Each record carries its location and mesh identity, not an unrelated global min.
module mod_floating_transport_diag
  implicit none
  private
  integer, parameter :: np=24, ns=6
  real*8, save :: score(ns)=huge(1.d0), payload(np,ns)=0.d0
  public :: transport_diag_reset, transport_diag_volume, transport_diag_wall
  public :: transport_diag_report, transport_diag_updated
contains
  subroutine transport_diag_reset()
    score=huge(1.d0)
    payload=0.d0
  end subroutine

  subroutine save_sample(slot,value,data)
    integer, intent(in) :: slot
    real*8, intent(in) :: value,data(np)
    integer :: k
    logical :: replace
    !$omp critical (floating_transport_samples)
    replace=value < score(slot)
    if (value==score(slot)) then
      do k=18,np ! deterministic location identity when extrema tie
        if (data(k)==payload(k,slot)) cycle
        replace=data(k)<payload(k,slot)
        exit
      enddo
    endif
    if (replace) then
      score(slot)=value
      payload(:,slot)=data
    endif
    !$omp end critical (floating_transport_samples)
  end subroutine

  subroutine transport_diag_volume(vertices,ms,mt,R,Z,rho,Ti,Te,dperp,dpar,dadd,vel,rate,source,divv,adv,jac,qjac,pe)
    use phys_module, only: floating_u_probe_R,floating_u_probe_Z,tstep
    integer, intent(in) :: vertices(4),ms,mt
    real*8, intent(in) :: R,Z,rho,Ti,Te,dperp,dpar,dadd,vel(2),rate,source,divv,adv,jac,qjac,pe
    real*8 :: p(np)
    p=0.d0
    p(1:17)=(/R,Z,rho,Ti,Te,dperp,dpar,dadd,vel,rate*tstep,source,divv,adv,jac,qjac,pe/)
    p(18:24)=real((/vertices,ms,mt,0/),8)
    call save_sample(1,rho,p)
    call save_sample(2,-rate*tstep,p)
    call save_sample(3,(R-floating_u_probe_R)**2+(Z-floating_u_probe_Z)**2,p)
  end subroutine

  subroutine transport_diag_wall(vertices,side,ms,R,Z,rho,Ti,Te,ven,vpn,cs,bn,qjac,res,btypes)
    use phys_module, only: floating_u_probe_R,floating_u_probe_Z
    integer, intent(in) :: vertices(4),side,ms,btypes(2)
    real*8, intent(in) :: R,Z,rho,Ti,Te,ven,vpn,cs,bn,qjac,res
    real*8 :: p(np)
    p=0.d0
    p(1:12)=(/R,Z,rho,Ti,Te,ven,vpn,ven+vpn,cs,bn,qjac,res/)
    p(13:14)=real(btypes,8)
    p(18:24)=real((/vertices,side,ms,0/),8)
    call save_sample(4,ven+vpn,p)
    call save_sample(5,-ven-vpn,p)
    call save_sample(6,(R-floating_u_probe_R)**2+(Z-floating_u_probe_Z)**2,p)
  end subroutine

  subroutine transport_diag_report(my_id,label)
    use mpi_mod
    use phys_module, only: central_density,central_mass,tstep
    use constants, only: MU_ZERO,ATOMIC_MASS_UNIT,EL_CHG
    integer, intent(in) :: my_id
    character(*), intent(in) :: label
    real*8 :: loc(2,ns),glob(2,ns),p(np),tn,teunit
    integer :: s,owner,ierr
    character(16), parameter :: names(ns)=[character(16):: 'volume min rho','volume max CFL', &
        'volume probe','wall min vn','wall max vn','wall probe']
    loc(1,:)=score
    loc(2,:)=real(my_id,8)
    call MPI_Allreduce(loc,glob,ns,MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,ierr)
    tn=sqrt(MU_ZERO*central_density*1.d20*central_mass*ATOMIC_MASS_UNIT)
    teunit=MU_ZERO*central_density*1.d20*EL_CHG
    if (my_id==0) write(*,'(A,A,A,ES12.4)') ' [floating transport] ',label,' dt=',tstep
    do s=1,ns
      if (glob(1,s)==huge(1.d0)) cycle
      owner=nint(glob(2,s))
      p=payload(:,s)
      call MPI_Bcast(p,np,MPI_DOUBLE_PRECISION,owner,MPI_COMM_WORLD,ierr)
      if (my_id/=0) cycle
      write(*,'(A,A,A,7I10)') ' [floating transport] ',trim(names(s)),' vertices/indices=',nint(p(18:24))
      write(*,'(A,2F11.6,A,ES12.4,A,2ES12.4)') '   R,Z=',p(1:2),' rho=',p(3),' Ti,Te[eV]=',p(4:5)/teunit
      if (s<=3) then
        write(*,'(A,3ES12.4,A,2ES12.4,A,ES12.4)') &
          '   Dperp,Dpar,Dadd[m2/s]=',p(6:8)/tn,' vR,vZ[m/s]=',p(9:10)/tn,' CFL=',p(11)
        write(*,'(A,3ES12.4,A,2ES12.4)') '   aux_rho,divv,adv_rho[code]=',p(12:14),' xjac,qjac=',p(15:16)
        write(*,'(A,ES12.4)') '   Pe along total poloidal velocity (actual diffusion tensor)=',p(17)
      else
        write(*,'(A,3ES12.4,A,ES12.4,A,2I4)') '   vEn,vparBn,total[m/s]=',p(6:8)/tn, &
             ' Bn/B=',p(10),' endpoint types=',nint(p(13:14))
        write(*,'(A,ES12.4,A,ES12.4)') '   qjac=',p(11),' floating trace residual[code u]=',p(12)
      endif
    enddo
  end subroutine

  ! Immediately AFTER update_values/update_deltas: sample actual new fields, not
  ! a linear predictor. The matching old density is new minus this accepted delta.
  ! Sampling Gauss points is detection, not a proof of polynomial positivity.
  subroutine transport_diag_updated(my_id,node_list,element_list,local_elms,n_local_elms)
    use data_structure, only: type_node_list,type_element_list
    use mod_parameters, only: var_rho,var_Ti,var_Te
    use gauss, only: n_gauss,xgauss
    use mod_interp, only: interp_PRZ
    use phys_module, only: floating_u_probe_R,floating_u_probe_Z
    use mpi_mod
    type(type_node_list), intent(in) :: node_list
    type(type_element_list), intent(in) :: element_list
    integer, intent(in) :: my_id,n_local_elms,local_elms(*)
    integer :: i,e,ms,mt,s,owner,ierr,vars(3)
    real*8 :: best(2),data(14,2),p(3),ps(3),pt(3),pp(3),d(3),rec(14),v(2),loc(2,2),glob(2,2)
    real*8 :: R,Rs,Rt,Z,Zs,Zt
    vars=(/var_rho,var_Ti,var_Te/)
    best=huge(1.d0); data=0.d0
    do i=1,n_local_elms
      e=local_elms(i)
      if (element_list%element(e)%n_sons>0) cycle
      do mt=1,n_gauss
        do ms=1,n_gauss
          call interp_PRZ(node_list,element_list,e,vars,3,xgauss(ms),xgauss(mt),0.d0, &
                            p,ps,pt,pp,R,Rs,Rt,Z,Zs,Zt)
          call interp_PRZ(node_list,element_list,e,vars,3,xgauss(ms),xgauss(mt),0.d0, &
                            d,ps,pt,pp,R,Rs,Rt,Z,Zs,Zt,deltas=.true.)
          rec=(/R,Z,p(1),p(1)-d(1),d(1),p(2:3),real(element_list%element(e)%vertex(1:4),8), &
                real(ms,8),real(mt,8),0.d0/)
          v=(/p(1),(R-floating_u_probe_R)**2+(Z-floating_u_probe_Z)**2/)
          do s=1,2
            if (v(s)<best(s)) then
              best(s)=v(s); data(:,s)=rec
            endif
          enddo
        enddo
      enddo
    enddo
    loc(1,:)=best; loc(2,:)=real(my_id,8)
    call MPI_Allreduce(loc,glob,2,MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,ierr)
    do s=1,2
      if (glob(1,s)==huge(1.d0)) cycle
      owner=nint(glob(2,s)); rec=data(:,s)
      call MPI_Bcast(rec,14,MPI_DOUBLE_PRECISION,owner,MPI_COMM_WORLD,ierr)
      if (my_id==0) then
        write(*,'(A,I1,A,6I10)') ' [floating transport POST] min/probe=',s,' vertices/GP=',nint(rec(8:13))
        write(*,'(A,2F11.6,A,3ES13.5,A,2ES13.5)') '   R,Z=',rec(1:2), &
          ' rho new,old,delta=',rec(3:5),' Ti,Te[code]=',rec(6:7)
      endif
    enddo
  end subroutine
end module
