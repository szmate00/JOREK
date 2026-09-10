program test_boundary
  use mod_parameters
  use phys_module
  use data_structure
  use mod_boundary_matrix_open
  use mod_floating_boundary_edges
  use mod_floating_transport_diag, only: transport_diag_reset,transport_diag_report,  &
                                         weak_mach_diag_reset,fw_tot_n,fw_cs_max,fw_res_max
  use basis_at_gaussian, only: set_basis
  implicit none
  integer,parameter :: nd=4*4*n_var
  type(type_element) :: e
  type(type_node) :: nodes(4),base(4)
  type(type_element_list) :: el
  type(type_node_list) :: nl
  real*8 :: a(nd,nd),ap(nd,nd),am(nd,nd),r(nd),rp(nd),rm(nd),err,scale,eps
  real*8 :: rc(nd),ac(nd,nd),slope,uc
  real*8 :: fdc(nd,nd),fdt(nd,nd)
  logical :: changed_plus,changed_minus
  integer :: i,dof,var,col,mode,k,row,icase,ib
  integer,parameter :: variables(5)=[var_u,var_vpar,var_Ti,var_Te,var_rho]
  integer,parameter :: boundary_types(5)=[1,3,4,5,9]
  call set_basis()
  do i=1,4
    nl%node(i)%index(1)=i
  enddo
  call floating_edges_build(el,nl)
  ! Affine side 1: R runs from 1 to 2, Z=0, exterior normal points down.
  base(1)%x(1,1,:)=[1.d0,0.d0]; base(2)%x(1,1,:)=[2.d0,0.d0]
  base(3)%x(1,1,:)=[2.d0,1.d0]; base(4)%x(1,1,:)=[1.d0,1.d0]
  do i=1,4
    base(i)%x(1,2,:)=[1.d0,0.d0]
    base(i)%x(1,3,:)=[0.d0,1.d0/3.d0]
    base(i)%values(1,1,var_rho)=0.2d0
    base(i)%values(1,1,var_Ti)=0.003d0
    base(i)%values(1,1,var_Te)=0.004d0
    base(i)%values(1,1,var_vpar)=0.04d0
    base(i)%values(1,1,var_psi)=0.08d0*base(i)%x(1,1,1)
    base(i)%values(1,2,var_psi)=0.08d0
    base(i)%values(1,1,var_u)=0.01d0*base(i)%x(1,1,1)
    base(i)%values(1,2,var_u)=0.01d0
  enddo
  bcs(1)%floating_u=.true.
  eps=1.d-7
  visco_par_heating=0.02d0 ! exercises the independent normal-derivative heat columns
  do icase=1,3
    if (icase==2) then
      base%boundary=4; bcs(4)%floating_u=.true.
      do i=1,4
        base(i)%values(:,:,var_u)=-base(i)%values(:,:,var_u)
      enddo
    elseif (icase==3) then
      base%boundary=9; bcs(9)%floating_u=.true.
      do i=1,4
        base(i)%values(1,1,var_Ti)=1.d-5
        base(i)%values(1,1,var_Te)=1.d-5
      enddo
      eps=1.d-9
    endif
  do mode=1,3
    floating_u_mach_flux=mode/=2
    floating_u_wall_flux=mode/=1
    nodes=base; call assemble(a,r)
    do k=1,5
      var=variables(k)
      do i=1,2
        do dof=1,4
          col=n_var*4*(i-1)+n_var*(dof-1)+var
          nodes=base; nodes(i)%values(1,dof,var)=nodes(i)%values(1,dof,var)+eps
          call assemble(ap,rp)
          nodes=base; nodes(i)%values(1,dof,var)=nodes(i)%values(1,dof,var)-eps
          call assemble(am,rm)
          ! ELM is minus d(RHS)/dx; theta=1. Only overridden rows are tested.
          do row=1,nd
            if (mod(row-1,n_var)+1==var_vpar .and. .not.floating_u_mach_flux) cycle
            if (mod(row-1,n_var)+1/=var_vpar .and. .not.floating_u_wall_flux) cycle
            scale=max(1.d0,maxval(abs(a(:,col))))
            err=abs(a(row,col)+(rp(row)-rm(row))/(2*eps))/scale
            if(err>2.d-7) then
              write(*,*) 'FAIL boundary FD mode,row,col,err',mode,row,col,err
              error stop 1
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
  enddo
  ! Relabelling a fully covered physical edge must not change its equations.
  floating_u_mach_flux=.true.; floating_u_wall_flux=.true.
  nodes=base; call assemble(a,r)
  do ib=1,size(boundary_types)
    bcs(boundary_types(ib))%floating_u=.true.
    nodes=base; nodes%boundary=boundary_types(ib)
    call assemble(ap,rp)
    if(any(a/=ap).or.any(r/=rp)) error stop 'boundary type changed an identical physical closure'
  enddo
  ! Coverage: a type-3 corner must not activate a type-2 artificial edge.
  nodes=base; nodes(1)%boundary=3; nodes(2)%boundary=2
  bcs(3)%floating_u=.true.; bcs(2)%floating_u=.false.
  floating_u_mach_flux=.false.; floating_u_wall_flux=.false.; call assemble(a,r)
  floating_u_mach_flux=.true.; floating_u_wall_flux=.true.; call assemble(ap,rp)
  if(any(a/=ap).or.any(r/=rp)) error stop 'type 2 edge changed'
  ! Diagnostic-enabled and disabled matrices/residuals are identical.
  nodes=base
  floating_u_mach_flux=.false.; floating_u_wall_flux=.false.
  floating_u_transport_diag=.false.; call assemble(a,r)
  call transport_diag_reset()
  floating_u_transport_diag=.true.; call assemble(ap,rp)
  if(any(a/=ap).or.any(r/=rp)) error stop 'diagnostic changed equations'
  call transport_diag_report(0,'serial fixture')
  ! --------------------------------------------------------------------------
  ! DEFAULT PATH (both opt-in flags false): the sheath ExB energy flux.
  ! The loop above never exercises this path - it only ran mach/wall combinations.
  ! u must be set as a genuinely LINEAR trace along the wall. Equal endpoint values
  ! with equal nonzero endpoint slopes is a cubic whose du/ds changes sign mid-edge,
  ! which would drive an outward ExB on part of the edge for BOTH signs.
  ! --------------------------------------------------------------------------
  floating_u_mach_flux=.false.; floating_u_wall_flux=.false.; floating_u_transport_diag=.false.
  ! PRE-EXISTING GAP, deliberately excluded here: on the default path the viscous
  ! heating column amat(var_Ti,var_vpar) differentiates gradvpardotn against the
  ! TRACE trial functions, whose normal derivative is zero. The fu_wall block adds
  ! the independent normal/mixed DOFs (fu_visc_col) precisely because of that; the
  ! default path never did, and the mode loop above never reached those rows. It is
  ! unrelated to the flux added below, so switch the heating off rather than let it
  ! mask the check.
  visco_par_heating=0.d0
  nodes=base; nodes%boundary=1; bcs(1)%floating_u=.true.; bcs(2)%floating_u=.false.
  uc=0.01d0
  slope=0.002d0
  ! 1. u CONSTANT along the wall: no tangential potential gradient, hence no ExB,
  !    so the term must be exactly absent. This is every non-floating production run.
  call set_u(0.d0,0.d0)
  call assemble(ac,rc)
  ! 2. BOTH signs act, with OPPOSITE sense: the collection is the total normal flow,
  !    so outward ExB adds to it and inward ExB subtracts. (The earlier outward-only
  !    form silently assumed Vpar*Bn = cs*|bn| and over-charged the sheath once the
  !    Mach row carried a drift supplement.)
  call set_u(+slope,0.d0); call assemble(ap,rp)
  changed_plus = any(rp/=rc) .or. any(ap/=ac)
  call set_u(-slope,0.d0); call assemble(am,rm)
  changed_minus = any(rm/=rc) .or. any(am/=ac)
  if (.not.changed_plus .or. .not.changed_minus) &
    error stop 'sheath ExB flux must respond to both signs of the wall potential gradient'
  do row=1,nd
    if (mod(row-1,n_var)+1/=var_Ti .and. mod(row-1,n_var)+1/=var_Te) cycle
    if ((rp(row)-rc(row))*(rm(row)-rc(row)) > 0.d0) &
      error stop 'outward and inward ExB must move the sheath flux in opposite senses'
  enddo
  ! Pick the outward sign for the remaining branch tests.
  if (sum(abs(rp-rc)) >= sum(abs(rm-rc))) then
    slope=+slope
  else
    slope=-slope
  endif
  ! 3. Only the SLOPE of u can drive an ExB flow; a constant offset cannot. To
  !    round-off, not exactly: fu_ven comes from eq_s, a floating-point sum over the
  !    basis, so a constant field's derivative cancels to ~1e-16 rather than to 0.
  call set_u(slope,0.d0);   call assemble(a,r)
  call set_u(slope,7.d0);   call assemble(ap,rp)
  if ( maxval(abs(a-ap)) > 1.d-12*max(1.d0,maxval(abs(a))) .or.   &
       maxval(abs(r-rp)) > 1.d-12*max(1.d0,maxval(abs(r))) )      &
    error stop 'sheath ExB flux must depend on du/dl only'
  ! 4. NO BOUND. The flux uses the same total flow as the Mach1 row, unclipped, so a
  !    further increase of the slope must keep changing it and the u column must stay
  !    non-zero. This is the inverse of the old assertion: it FAILS if a clip is
  !    reintroduced here while the Mach row has none.
  call set_u(slope*1.d3,0.d0); call assemble(a,r)
  call set_u(slope*2.d3,0.d0); call assemble(ap,rp)
  if ( maxval(abs(r-rp)) <= 0.d0 ) &
    error stop 'sheath ExB flux must not saturate - it carries the unclipped drift'
  do row=1,nd
    if (mod(row-1,n_var)+1/=var_Ti .and. mod(row-1,n_var)+1/=var_Te) cycle
    if ( all( a(row,var_u:nd:n_var) == 0.d0 ) .and. abs(r(row)) > 0.d0 ) &
      error stop 'unclipped sheath ExB flux must keep its u column at large slope'
  enddo
  ! 5. Verify the ADDED term's Jacobian by DIFFERENCING against a control with the
  !    term absent. The default-path columns have a pre-existing finite-difference
  !    discrepancy of ~2.5e-6 relative in amat(var_Te,var_Te) which is present with
  !    or without this flux and is not addressed here; differencing cancels it, so
  !    what is tested is exactly the contribution of the new term and nothing else.
  eps=1.d-9
  ! Sit well away from the one remaining branch boundary, the sign of fu_vtot.
  slope=slope*0.1d0
  call set_u(0.d0,0.d0);  base=nodes; call assemble(ac,rc); call sweep(fdc)
  call set_u(slope,0.d0); base=nodes; call assemble(a,r)
  if (all(a(:,var_u:nd:n_var)==0.d0)) error stop 'unclipped sheath ExB flux has no u column'
  call sweep(fdt)
  do col=1,nd
    do row=1,nd
      if (mod(row-1,n_var)+1/=var_Ti .and. mod(row-1,n_var)+1/=var_Te) cycle
      if (mod(col-1,n_var)+1==var_u) then
        ! The control has NO u column at all (the term is absent), and its own
        ! finite difference there straddles the kink of max(fu_ven,0), so it is not
        ! a derivative of anything. Test this column against the test-state FD alone.
        scale=max(1.d0,maxval(abs(a(:,col))))
        err=abs( a(row,col) + fdt(row,col) )/scale
      else
        scale=max(1.d0,maxval(abs(a(:,col)-ac(:,col))))
        err=abs( (a(row,col)-ac(row,col)) + (fdt(row,col)-fdc(row,col)) )/scale
      endif
      if (err>2.d-6) then
        write(*,*) 'FAIL sheath ExB FD row,col,err',row,col,err
        error stop 1
      endif
    enddo
  enddo
  ! 6. The wall never emits: drive a huge INWARD ExB and check the collection is
  !    floored at zero, i.e. the added term exactly cancels the parallel expression.
  call set_u(-slope*1.d3,0.d0); call assemble(a,r)
  call set_u(-slope*2.d3,0.d0); call assemble(ap,rp)
  if (any(a/=ap) .or. any(r/=rp)) error stop 'closed-wall branch is not saturated'
  ! 7. FD the CLOSED branch specifically. Its Vpar column cancels the parallel
  !    collection, and it lives in a block whose amat(*,var_vpar) entries are plain
  !    assignments further down - accumulating into them earlier is silently
  !    overwritten, leaving a residual that says the collection is zero and a
  !    Jacobian that still responds as if the parallel collection were present.
  !    Nothing in the open-branch sweep above can see that, so test it here.
  ! Enter it reliably: the branch is fu_bn*Vpar0 + clip(vE.n) <= 0, so weaken the
  ! parallel flow and take whichever drift sign actually closes the wall (test 2's
  ! `slope` marks neither sense now that both act).
  changed_plus=.false.
  do k=1,2
    if (k==1) then
      call set_u(+abs(slope)*1.d3,0.d0)
    else
      call set_u(-abs(slope)*1.d3,0.d0)
    endif
    do i=1,4
      nodes(i)%values(1,1,var_vpar)=1.d-4
    enddo
    base=nodes
    nodes=base; call assemble(a,r)
    nodes=base; nodes(1)%values(1,1,var_vpar)=nodes(1)%values(1,1,var_vpar)+1.d-3
    call assemble(ap,rp)
    changed_minus=.true.
    do row=1,nd
      if (mod(row-1,n_var)+1/=var_Ti .and. mod(row-1,n_var)+1/=var_Te) cycle
      if (abs(rp(row)-r(row))>1.d-14) changed_minus=.false.
    enddo
    if (changed_minus) then
      changed_plus=.true.
      exit
    endif
  enddo
  if (.not.changed_plus) error stop 'closed-wall branch never entered by either drift sign'
  ! Defining property of the branch: the collection is zero, so the Ti/Te residual
  ! does not depend on Vpar - verified just above - and therefore neither may its
  ! Vpar column. That column lives in a block whose amat(*,var_vpar) entries are
  ! plain assignments further down, so accumulating into them earlier is silently
  ! overwritten, leaving a residual that says the collection is zero and a Jacobian
  ! that still responds as if the parallel collection were present.
  nodes=base; call assemble(a,r)
  do row=1,nd
    if (mod(row-1,n_var)+1/=var_Ti .and. mod(row-1,n_var)+1/=var_Te) cycle
    do i=1,2
      do dof=1,4
        col=n_var*4*(i-1)+n_var*(dof-1)+var_vpar
        if (abs(a(row,col))>1.d-14) error stop 'closed-wall branch keeps a Vpar column the residual does not have'
      enddo
    enddo
  enddo
  ! 8. floating_temperature_slope must BE corr_neg_temp1', since cs0 is built from
  !    corr_neg_temp1(T) and d(cs)/dT is now formed as gamma*T/(2*cs0)*slope on the
  !    production path (it used to omit the slope entirely, which overstates d(cs)/dT
  !    by 1/corr' - a factor 3 at 0.65 eV and exponentially worse below, exactly where
  !    the divertor wall sits). Assembly-level sensitivity is too weak to test in this
  !    fixture: on the default path cs_T enters only the small c_angle term. Test the
  !    identity directly instead.
  call test_corr_slope()
  call test_sheath_stabiliser()
  call test_weak_mach()
  call test_weak_mach_sign()
  write(*,*) 'PASS: sheath ExB energy flux - absent without a gradient, both senses, unbounded, closed wall, FD Jacobian'
  write(*,*) 'PASS: production boundary assembler finite-difference Jacobians and type-2 exclusion'
contains
  !> The SOLPS-style zero-sum sheath stabiliser: RESIDUAL bit-identical, Jacobian
  !! diagonal linear in alpha, and nothing outside the Ti/Te/rho diagonal touched.
  subroutine test_sheath_stabiliser()
    real*8 :: a0(nd,nd),r0(nd),a1(nd,nd),r1(nd),a2(nd,nd),r2(nd)
    integer :: rr,cc,vr
    floating_u_mach_flux=.false.; floating_u_wall_flux=.false.
    nodes=base
    stab_coeff_sheath_te=0.d0; stab_coeff_sheath_ti=0.d0; stab_coeff_sheath_ni=0.d0
    call assemble(a0,r0)
    stab_coeff_sheath_te=10.d0; stab_coeff_sheath_ti=10.d0; stab_coeff_sheath_ni=10.d0
    call assemble(a1,r1)
    stab_coeff_sheath_te=20.d0; stab_coeff_sheath_ti=20.d0; stab_coeff_sheath_ni=20.d0
    call assemble(a2,r2)
    ! 1. ZERO-SUM: the residual must not move at all. This is the whole claim - the
    !    converged solution is untouched because the term vanishes at X_new = X_old.
    if (any(r1/=r0) .or. any(r2/=r0)) error stop 'sheath stabiliser changed the residual'
    ! 2. It must actually do something.
    if (all(a1==a0)) error stop 'sheath stabiliser did nothing to the Jacobian'
    ! 3. Exactly linear in alpha.
    do cc=1,nd
      do rr=1,nd
        if (abs((a2(rr,cc)-a0(rr,cc)) - 2.d0*(a1(rr,cc)-a0(rr,cc))) > &
            1.d-12*max(1.d-30,abs(a2(rr,cc)-a0(rr,cc)))) &
          error stop 'sheath stabiliser is not linear in alpha'
        ! 4. Confined to the Ti/Te/rho DIAGONAL blocks.
        if (a1(rr,cc)/=a0(rr,cc)) then
          vr=mod(rr-1,n_var)+1
          if (vr/=var_Ti .and. vr/=var_Te .and. vr/=var_rho) &
            error stop 'sheath stabiliser touched a row it should not'
          if (mod(cc-1,n_var)+1/=vr) &
            error stop 'sheath stabiliser touched an off-diagonal column'
        endif
      enddo
    enddo
    stab_coeff_sheath_te=0.d0; stab_coeff_sheath_ti=0.d0; stab_coeff_sheath_ni=0.d0
    write(*,*) 'PASS: zero-sum sheath stabiliser - residual untouched, Jacobian linear in alpha, diagonal only'
  end subroutine
  !> floating_temperature_slope == d/dT corr_neg_temp1, by central differences.
  subroutine test_corr_slope()
    use corr_neg, only: corr_neg_temp1
    use mod_floating_transport, only: floating_temperature_slope
    real*8 :: t,h,fd,an,knee
    integer :: n
    knee=T_min_neg
    do n=1,40
      t=0.05d0*dble(n)*knee*sum(corr_neg_temp_coef)*2.d0
      h=1.d-9*max(t,knee)
      fd=(corr_neg_temp1(t+h)-corr_neg_temp1(t-h))/(2.d0*h)
      an=floating_temperature_slope(t,knee,corr_neg_temp_coef)
      if (abs(fd-an)>1.d-5*max(1.d-3,abs(an))) then
        write(*,*) 'FAIL corr slope t,fd,analytic',t,fd,an
        error stop 1
      endif
    enddo
    write(*,*) 'PASS: floating_temperature_slope is the corr_neg_temp1 derivative'
  end subroutine
  !> d(RHS)/dx by central differences, for every DOF, into fd(row,col).
  !> WEAK (Galerkin) drift-inclusive Bohm condition: every column of the residual
  !! must be in the matrix. This is the check the nodal route never had - its slope
  !! row's residual depended on four u DOFs with no column at all, which is what made
  !! it a systematic per-step error with no Newton iteration to absorb it.
  !! Sweeps psi as well as u/Vpar/Ti/Te, because the weak row has a psi column
  !! (through B.n) that the shared `sweep` does not perturb.
  subroutine test_weak_mach()
    real*8 :: aw(nd,nd),rw(nd),fw(nd,nd),e2,sc,er,worst
    integer,parameter :: wvars(5)=[var_u,var_vpar,var_Ti,var_Te,var_psi]
    integer :: kk,vv,ii,dd,cc,rr
    mach1_weak=.true.
    nodes=base; nodes%boundary=1
    bcs(1)%floating_u=.true.; bcs(1)%mach1=.true.; bcs(2)%floating_u=.false.
    ! --- Sit WELL inside the compensating branch. In this fixture cs*|b.n| ~ 2.9e-3,
    ! --- so a slope of 2e-3 puts vE.n right on the max(.,0) kink and the finite
    ! --- difference straddles it - that is a bad test state, not a bad Jacobian.
    call set_u(1.d-4,0.d0)
    base=nodes
    call weak_mach_diag_reset()
    call assemble(aw,rw)
    ! --- The diagnostic must actually receive samples. It first shipped printing its
    ! --- header with no rows, because the reset lived at the top of
    ! --- boundary_conditions while the accumulation happens earlier, during element
    ! --- assembly - so every sample was wiped just before being reported.
    if (fw_tot_n(1) <= 0.d0) error stop 'weak Mach diagnostic received no samples'
    if (fw_cs_max(1) <= 0.d0) error stop 'weak Mach diagnostic has no sonic target'
    if (fw_res_max(1) < 0.d0) error stop 'weak Mach diagnostic residual not recorded'
    ! --- the row has to exist at all
    worst=0.d0
    do rr=1,nd
      if (mod(rr-1,n_var)+1/=var_vpar) cycle
      worst=max(worst,maxval(abs(aw(rr,:))))
    enddo
    if (worst<=0.d0) error stop 'weak Mach row was not assembled'
    ! --- finite-difference every column of the Vpar rows
    e2=1.d-9
    fw=0.d0
    do kk=1,5
      vv=wvars(kk)
      do ii=1,2
        do dd=1,4
          cc=n_var*4*(ii-1)+n_var*(dd-1)+vv
          nodes=base; nodes(ii)%values(1,dd,vv)=nodes(ii)%values(1,dd,vv)+e2
          call assemble(ap,rp)
          nodes=base; nodes(ii)%values(1,dd,vv)=nodes(ii)%values(1,dd,vv)-e2
          call assemble(am,rm)
          fw(:,cc)=(rp-rm)/(2*e2)
        enddo
      enddo
    enddo
    nodes=base
    ! --- rhs = -w*res and amat = +w*d(res)/dx, so amat + d(rhs)/dx must vanish
    worst=0.d0
    do kk=1,5
      vv=wvars(kk)
      do ii=1,2
        do dd=1,4
          cc=n_var*4*(ii-1)+n_var*(dd-1)+vv
          do rr=1,nd
            if (mod(rr-1,n_var)+1/=var_vpar) cycle
            sc=max(1.d0,maxval(abs(aw(:,cc))))
            er=abs(aw(rr,cc)+fw(rr,cc))/sc
            if (er>worst) worst=er
            if (er>1.d-5) then
              write(*,*) 'FAIL weak Mach FD row,col,var,err',rr,cc,vv,er,aw(rr,cc),fw(rr,cc)
              error stop 1
            endif
          enddo
        enddo
      enddo
    enddo
    ! --- SATURATED branch: drive the drift far past sonic outflow. Bohm is an
    ! --- inequality, so nothing is imposed - the target pins at zero and the u and
    ! --- temperature columns must be exactly zero, with the residual carrying only
    ! --- the parallel term. This is the branch the nodal one-sided form also took.
    ! --- vE.n scales as MINUS the slope here, so a negative slope is the one that
    ! --- drives the drift past sonic outflow and saturates the target at zero.
    call set_u(-1.d0,0.d0)
    base=nodes
    call assemble(aw,rw)
    do rr=1,nd
      if (mod(rr-1,n_var)+1/=var_vpar) cycle
      do kk=1,4
        vv=wvars(kk)
        if (vv==var_vpar) cycle
        do ii=1,2
          do dd=1,4
            cc=n_var*4*(ii-1)+n_var*(dd-1)+vv
            if (aw(rr,cc)/=0.d0) then
              write(*,*) 'FAIL saturated weak Mach row keeps a column',rr,cc,vv,aw(rr,cc)
              error stop 1
            endif
          enddo
        enddo
      enddo
    enddo
    nodes=base
    mach1_weak=.false.
    bcs(1)%mach1=.false.
    write(*,'(a,es9.2)') ' PASS: weak Mach Galerkin row - all Vpar columns match FD (worst rel err ',worst
    write(*,'(a)')       '       ), saturated branch drops its u and T columns'
  end subroutine
  !> B.n CHANGING SIGN ON THE WALL. The legacy Mach row imposes direction*cs with
  !! direction = sign(B.n), which JUMPS by 2*cs where the field crosses tangency -
  !! that discontinuity is the whole reason vpar_smoothing's tanh factor exists.
  !!
  !! The weak row cannot have that defect: the target is cs*|b.n|, continuous through
  !! zero, there is no `direction` at all (the sign is carried by fu_bn in
  !! fu_bn*Vpar0), and the assembled row is
  !!     fu_bn*(fu_bn*Vpar0 - target) = fu_bn^2*Vpar0 - fu_bn*target
  !! whose first term vanishes as (B.n)^2 and second linearly. The only
  !! non-smoothness left is the kink of |b.n| at exactly B.n = 0, and it is
  !! multiplied by fu_bn, which vanishes there too.
  !!
  !! So this sweeps the psi slope through zero and requires the row to fade away
  !! continuously instead of jumping. vpar_smoothing is left OFF throughout: if the
  !! weak row needed it, this test would fail.
  subroutine test_weak_mach_sign()
    real*8 :: g(9),mag(9),rowmax,big,small,asym
    integer :: kk,ii,rr
    g = (/ -4.d-2,-2.d-2,-1.d-3,-1.d-6, 0.d0, 1.d-6, 1.d-3, 2.d-2, 4.d-2 /)
    mach1_weak=.true.
    vpar_smoothing=.false.
    bcs(1)%floating_u=.true.; bcs(1)%mach1=.true.; bcs(2)%floating_u=.false.
    do kk=1,9
      nodes=base; nodes%boundary=1
      do ii=1,4
        nodes(ii)%values(1,1,var_psi)=g(kk)*nodes(ii)%x(1,1,1)
        nodes(ii)%values(1,2,var_psi)=g(kk)
      enddo
      call set_u(1.d-4,0.d0)
      base=nodes
      call assemble(a,r)
      ! --- EXCLUDE the psi column. The assembled row is w*res with w = B.n, so its
      ! --- exact psi derivative is w*d(res)/dpsi + res*d(w)/dpsi, and the second term
      ! --- does NOT vanish at tangency - there the row would read as a constraint on
      ! --- dpsi rather than on Vpar. That is inert in production because psi's trace
      ! --- AND its tangential derivative are Dirichlet on the wall (which is exactly
      ! --- why the code treats b.n as a static map), so those DOFs cannot move. This
      ! --- fixture imposes no psi BC, so the column is live here and would mask the
      ! --- property being tested: that the CONSTRAINT ON VPAR fades away.
      rowmax=0.d0
      do rr=1,nd
        if (mod(rr-1,n_var)+1/=var_vpar) cycle
        do ii=1,nd
          if (mod(ii-1,n_var)+1==var_psi) cycle
          rowmax=max(rowmax,abs(a(rr,ii)))
        enddo
        rowmax=max(rowmax,abs(r(rr)))
      enddo
      mag(kk)=rowmax
    enddo
    ! --- the row must FADE as the field goes tangential, not jump
    big   = max(mag(1),mag(9))
    small = max(mag(4),mag(6))
    if (small > 1.d-3*big) then
      write(*,*) 'FAIL weak Mach row does not fade at tangency',small,big
      error stop 1
    endif
    ! --- and it must be continuous ACROSS the sign change: the two samples either
    ! --- side of B.n = 0 must agree, where the legacy form would differ by ~2*cs
    asym = abs(mag(4)-mag(6)) / max(big,tiny(1.d0))
    if (asym > 1.d-4) then
      write(*,*) 'FAIL weak Mach row jumps across B.n = 0',mag(4),mag(6),asym
      error stop 1
    endif
    if (mag(5) /= 0.d0) then
      write(*,*) 'FAIL weak Mach row nonzero at exact tangency',mag(5)
      error stop 1
    endif
    mach1_weak=.false.
    bcs(1)%mach1=.false.
    nodes=base
    write(*,'(a,es9.2,a,es9.2)') ' PASS: weak Mach row fades through B.n = 0 without vpar_smoothing'// &
      ' (tangency/bulk ',small/big,', asymmetry ',asym
  end subroutine
  subroutine sweep(fd)
    real*8,intent(out)::fd(nd,nd)
    integer::kk,vv,ii,dd,cc
    fd=0.d0
    do kk=1,5
      vv=variables(kk)
      do ii=1,2
        do dd=1,4
          cc=n_var*4*(ii-1)+n_var*(dd-1)+vv
          nodes=base; nodes(ii)%values(1,dd,vv)=nodes(ii)%values(1,dd,vv)+eps
          call assemble(ap,rp)
          nodes=base; nodes(ii)%values(1,dd,vv)=nodes(ii)%values(1,dd,vv)-eps
          call assemble(am,rm)
          fd(:,cc)=(rp-rm)/(2*eps)
        enddo
      enddo
    enddo
    nodes=base
  end subroutine
  !> Set u as an exactly linear trace along the wall: u = uc + off + g*(R-1).
  subroutine set_u(g,off)
    real*8,intent(in)::g,off
    integer::ii
    do ii=1,4
      nodes(ii)%values(1,1,var_u)=uc+off+g*(nodes(ii)%x(1,1,1)-1.d0)
      nodes(ii)%values(1,2,var_u)=g
    enddo
  end subroutine
  subroutine assemble(mat,rhs)
    real*8,intent(out)::mat(nd,nd),rhs(nd)
    integer::vertices(2),directions(2)
    vertices=[1,2]; directions=[1,2]
    mat=0.d0; rhs=0.d0
    call boundary_matrix_open(vertices,directions,e,nodes,.true.,1,1.d0,0.d0,0.d0,1.d0, &
        [1.d0,1.d0],[0.d0,0.d0],mat,rhs,1,1,1)
  end subroutine
end program
