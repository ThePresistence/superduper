!
! created by Jie Liu on 09/16
!

module cis_solve

use es_global

public

contains

!======================================================== 
!     G : Inital guess vector
!     H : Orbital energy differences
!     I : Rank of orbital energy differences
!======================================================== 
subroutine eguess(G,Ha,Hb,Ia,Ib)
  implicit none
  integer Ia(*), Ib(*)
  real*8  Ha(*), Hb(*), G(*)
  integer iroot, ioff
  real*8  d,c
  logical done

  done = .false.
  if(imult.eq.0 .or. imult.eq.3) then
    iroot = 1
    c = sqrt(0.5d0)
    if(imult.eq.3 .or. sasfcis) c = 1.d0
    do while(.not.done)
      ioff = Ia(iroot) + (iroot-1)*NOV
      call vecinit(G((iroot-1)*NOV+1),NOV,0d0)
      G(ioff) = c
      if(imult.eq.0) G(ioff+NOVa) = c
      if(I1n.eq.0.and.imult.eq.0) G(ioff+NOVa) = -1*G(ioff+NOVa)
      iroot = iroot + 1
      if(iroot.gt.nroots) then
         ! The number of required states should not be more than NOVa
         if(iroot.gt.NOVa) then
           done = .true.
         else
           ! In case of degenerate states
           d = dabs(Ha(Ia(iroot))-Ha(Ia(iroot-1)))
           if(d.ge.1d-5) done = .true.
         endif
      endif
    end do
    if(iroot.gt.nroots+1) then
       write(nb6,*) ' The number of excited states is changed from'  &
                              ,nroots,' to ',iroot-1
       nroots = iroot-1
    endif
  else if(imult.eq.1) then
      write(nb6,*) "Guess X for open shell singlet"
      stop
  else if(imult.eq.2) then
      write(nb6,*) "Guess X for doublet"
  endif

  return
end subroutine eguess

!======================================================== 
! Sort up the orbital energy difference matrix E_ai
! the index is stored in I
!======================================================== 
subroutine esort(H,I,n)
  implicit none
  integer j,k,m,n,lsup,ncount
  integer I(*)
  real*8  H(*)

  do k=1,n
    I(k) = k
  enddo
  
  lsup = n
  do while (lsup > 1)
    ! ncout in the greatest element out of order
    ncount = 0  
    do j = 1, (lsup-1)
      if (H(I(j)) > H(I(j+1))) then
        m = I(j)
        I(j) = I(j+1)
        I(j+1) = m
        ncount = j
      endif 
    enddo
    lsup = ncount   
  enddo

  return
end subroutine esort

!======================================================== 
!     Calculate MO orbital energy difference
!     E_ai = (E_a - E_i)
!======================================================== 
subroutine eorbdif(V,Ha,Hb)
  implicit none
  integer i,j,k
  real*8  Ha(*),Hb(*),V(*)

  if(imult.le.2) then
    k = 1
    do i = 1,NOa
       do j = 1,NVa
          Ha(k) = V(jEa+NOa+j-1) - V(jEa+i-1) + cfac1 
          k = k + 1
       enddo
    enddo
  
    if(nden.eq.2) then
      k = 1
      do i = 1,NOb
        do j = 1,NVb
          Hb(k) = V(jEb+NOb+j-1) - V(jEb+i-1) + cfac1 
          k = k + 1
        enddo
      enddo
    else 
      call veccopy(Hb,Ha,NOVa)
    endif
  else if(imult.eq.3) then
    k = 1
    do i = 1,NOa
       do j = 1,NVb
          Ha(k) = V(jEb+NOb+j-1) - V(jEa+i-1) + cfac1
          Hb(k) = V(jEb+NOb+j-1) - V(jEa+i-1) + cfac1
          k = k + 1
       enddo
    enddo
    call matprnt(Ha,NVb,NOa,6)
  endif

  return
end subroutine eorbdif

!======================================================== 
subroutine emake_rao(Pa,Pb,R,Ca,Cb,n)
  implicit none
  integer i,n
  real*8 Pa(NBas*NBas,*),Pb(NBas*NBas,*),R(NOV,*),Ca(LM2,*),Cb(LM2,*)
  real*8, dimension(:), allocatable :: tmp

  allocate(tmp(NBas*NBas))

  if(imult.le.2) then
    ! Pao = Cv * X * Co^T
    do i=1,n
      call matmult(Ca(1,1+NOa),R(1,i),tmp,NBas,NOa,NVa,LM2,NVa,NBas,1)
      call matmult(tmp,Ca(1,1),Pa(1,i),NBas,NBas,NOa,NBas,LM2,NBas,3)
      if(nden.eq.2) then
        call matmult(Cb(1,1+NOb),R(1+NOVa,i),tmp,NBas,NOb,NVb,LM2,NVb,NBas,1)
        call matmult(tmp,Cb(1,1),Pb(1,i),NBas,NBas,NOb,NBas,LM2,NBas,3)
      endif
    enddo
  else if(imult.eq.3) then
    ! Pao = Cbv * X * Cao^T
    do i=1,n
      call matmult(Cb(1,1+NOb),R(1,i),tmp,NBas,NOa,NVb,LM2,NVb,NBas,1)
      call matmult(tmp,Ca(1,1),Pa(1,i),NBas,NBas,NOa,NBas,LM2,NBas,3)
      if(nden.eq.2) call vecinit(Pb(1,i),NBas*NBas,0.d0)
    enddo
  endif

  deallocate(tmp)
  return
end subroutine emake_rao

!======================================================== 
subroutine emake_amo(Az,Fa,Fb,R,Ha,Hb,Ca,Cb,n,rpa,jspin,ndim)
  implicit none
  integer i,m,n,jspin,ndim
  logical rpa
  real*8 Az(NOV,*),Fa(NBas*NBas,*),Fb(NBas*NBas,*),R(NOV,*)
  real*8 Ha(*),Hb(*),Ca(LM2,*),Cb(LM2,*)
  real*8, dimension(:), allocatable :: tmp,tmp2

  allocate(tmp(NBas*NBas),tmp2(NOV))

  do m=1,n
    if(jspin.eq.0) then
      ! A = Cv^T * F * Co
      call matmult(Ca(1,NOa+1),Fa(1,m),tmp,NVa,NBas,NBas,LM2,NBas,NVa,2)
      call matmult(tmp,Ca(1,1),Az(1,m),NVa,NOa,NBas,NVa,LM2,NVa,1)
      if(ndim.eq.2) then
        call matmult(Cb(1,NOb+1),Fb(1,m),tmp,NVb,NBas,NBas,LM2,NBas,NVb,2)
        call matmult(tmp,Cb(1,1),Az(1+NOVa,m),NVb,NOb,NBas,NVb,LM2,NVb,1)
      endif
    else if(jspin.eq.3) then
      ! A = Cbv^T * F * Cao
      call matmult(Cb(1,NOb+1),Fa(1,m),tmp,NVb,NBas,NBas,LM2,NBas,NVb,2)
      call matmult(tmp,Ca(1,1),Az(1,m),NVb,NOa,NBas,NVb,LM2,NVb,1)
    endif
    if(rpa) then
      call matmult(Ca(1,NOa+1),Fa(1,m),tmp,NVa,NBas,NBas,LM2,NBas,NVa,4)
      call matmult(tmp,Ca(1,1),tmp2,NVa,NOa,NBas,NVa,LM2,NVa,1)
      call vecadd(Az(1,m),Az(1,m),tmp2,NOa*NVa)
      if(ndim.eq.2) then
        call matmult(Cb(1,NOb+1),Fb(1,m),tmp,NVb,NBas,NBas,LM2,NBas,NVb,4)
        call matmult(tmp,Cb(1,1),tmp2,NVb,NOb,NBas,NVb,LM2,NVb,1)
        call vecadd(Az(1+NOVa,m),Az(1+NOVa,m),tmp2,NOb*NVb)
      endif
    endif

    ! A = (epsilon_a-epsilon_i) + Cv^T * F * Co
    do i=1,NOVa
      Az(i,m) = Az(i,m) + R(i,m)*Ha(i)
    enddo
    if(jspin.eq.0) then
      if(ndim.eq.2) then
        do i=1,NOVb
          Az(i+NOVa,m) = Az(i+NOVa,m) + R(i+NOVa,m)*Hb(i)
        enddo
      else
        call veccopy(Az(1+NOVa,m),Az(1,m),NOVa)
        if(I1n.eq.0) call vecscale(Az(1+NOVa,m),NOVa,-1d0)
      endif
    else if(jspin.eq.3) then
      call vecinit(Az(1+NOVa,m),NOVa,0.d0)
    endif
  enddo
  if(iprtci.ge.5) then
    write(nb6,*) "Az"
    call matprnt(Az,NOV,n,6)
  endif

  deallocate(tmp)
  deallocate(tmp2)
  return
end subroutine emake_amo


!======================================================== 
subroutine cisdiag(H)
  implicit none
  integer ierror,i
  real*8 H(*)
  real*8, dimension(:), allocatable :: sAs,E

  allocate(sAs(ns*ns),E(ns))

  call matmult(Tv,At,sAs,ns,ns,NOV,NOV,NOV,ns,2)
  if(iprtci.ge.5) then
    write(nb6,*) "A in subspace"
    call matprnt(sAs,ns,ns,6)
  endif

  ! diagonalization in subspace
  ierror = 0
  call matdiagsquare(sAs,E,ns,ierror)
  if(iprtci.ge.5) then
    write(nb6,*) "excitation energy:"
    call matprnt(E,ns,1,6)
  endif
  if(ierror.ne.0) then
    write(nb6,*) "CIS DIAGON fails"
    stop
  endif

  ! check for convergence
  call chkcnv(sAs,E,H)

  deallocate(sAs)
  deallocate(E)
  return
end subroutine cisdiag

!======================================================== 
subroutine chkcnv(Xs,E,H)
  implicit real*8 (a-h,o-z)
  integer m, isub, ioff
  real*8  tol, maxd, d
  real*8  Xs(*), E(*), H(*)
  real*8, dimension(:), allocatable :: Res
  allocate(Res(NOV))

  tol   = 1d1**(-nrmdav)
  maxd  = 0d0
  isub  = ns
  !ibeta = NOVa+1
  !if(sasfcis) ibeta = 1

  !call eorbdif(H,H(ibeta),V)

  ! Check convergence and expand subspace
  do m = 1,nroots
     call chkres(Res,Xs(ns*(m-1)+1),E(m),d,H)
     if(d.gt.maxd) maxd = d
     if(d.gt.tol)  then
       call ssexpand(Res,isub,NOV)
     endif
  enddo

  ! Save the information for subspace
  nl = isub - ns
  if(iprtci.gt.0) write(nb6,*) niter,nl,maxD
  ns = isub

  if(nl.eq.0) then
    allocate(Xv(NOV,nroots*(I13+1)),es(nroots,I13+1))
    ioff = 0
    if(I3.and.I1.and.I1n) ioff = nroots
    call matmult(Tv,Xs,Xv(1,ioff+1),NOV,nroots,ns,NOV,ns,NOV,1)
    es(1:nroots,I13+1) = E(1:nroots)
    if(.not.is_sf .and. E(1).lt.0) then
      write(6,*) "negative 1st excitation energy"
      stop
    endif
    if(iprtci.ge.5) then
      write(nb6,*) "cis amplitude"
      call matprnt(Xv,NOV,nroots,6)
    endif
    if(iprtci.gt.0) write(nb6,*) "CIS Converged"
  else if(ns.gt.maxdav) then
    ! If maximum dimension of Davidson subspace is reached,
    ! the subspace information is saved to restart cis calculation
    write(nb6,*) "Max Subspace is reached. Restart CIS calculation"
    deallocate(Res)
    allocate(Res(NOV*nroots))
    call matmult(Tv,Xs,Res,NOV,nroots,ns,NOV,ns,NOV,1)
    call veccopy(Tv,Res,NOV*nroots)
    call veccopy(Tv(1,nroots+1),Tv(1,ns-nl+1),nl*NOV)
    call matmult(At,Xs,Res,NOV,nroots,ns,NOV,ns,NOV,1)
    call veccopy(At,Res,NOV*nroots)
    ns = nroots + nl
  endif

  deallocate(Res)
  return
end subroutine chkcnv

!======================================================== 
subroutine chkres(Res,Xs,E,D,H)
  implicit none
  integer i, j
  real*8  D, E, Res(*), Xs(*), H(*)
  real*8, dimension(:), allocatable :: tmp

  allocate(tmp(NOV))

  ! error vector A*X-\omega*X
  call vecinit(tmp,NOV,0d0)
  do i = 1,ns
     do j = 1,NOV
        tmp(j) = tmp(j) + (At(j,i)-E*Tv(j,i))*Xs(i)
     enddo
  enddo

  call vecdot(d,tmp,tmp,NOV)
  d = sqrt(d/NOV)

  do i=1,NOV
     Res(i) = tmp(i)/(E-H(i))
  enddo

  deallocate(tmp)
  return
end subroutine chkres

!======================================================== 
subroutine ssexpand(Res,isub,ndim)
  implicit real*8 (a-h,o-z)
  integer i, isub, ndim
  real*8 tol, R, RR, RG
  real*8 Res(*)

  tol  = (2d-2)*(1d1**(-nrmdav))

  ! schmit orthogonization
  do i=1,isub
    call vecdot(RG,Res,Tv(1,i),ndim)
    call vecadd2(Res,-RG,Tv(1,i),Res,ndim)
  enddo

  ! add new trial vector to subspace
  call vecdot(RR,Res,Res,ndim)
  R = SQRT(RR/ndim)
  if(R.ge.tol) then
    if(isub.eq.maxdav) then
       write(nb6,*) "Maximum subspace in SSExpand"
       STOP
    endif
    call vecscale(Res,ndim,SQRT(1/RR))
    call veccopy(Tv(1,isub+1),Res,ndim)
    isub = isub + 1
  endif

  return
end subroutine ssexpand

!======================================================== 
subroutine get_cis_transition_density2(P,Ca,Cb)
  implicit none
  integer i,n2
  real*8  P(nbas*nbas,*),Ca(LM2,*),Cb(LM2,*)
  real*8, dimension(:), allocatable :: tmp,tmp2,tmp3

  n2 = nbas*nbas
  allocate(tmp(N2),tmp2(N2),tmp3(N2))
  call vecinit(P,N2*nroots,0.d0)

  do i=2,nroots
    if(imult.le.2) then
      call veccopy(tmp3,Xv(1,i),NOVa)
      if(irpa) call vecadd(tmp3,tmp3,Yv(1,i),NOVa)
      call matmult(Ca(1,1+NOa),tmp3,tmp,NBas,NOa,NVa,LM2,NVa,NBas,1)
      call matmult(tmp,Ca(1,1),P(1,i),NBas,NBas,NOa,NBas,LM2,NBas,3)
      if(nden.eq.2) then
        call veccopy(tmp3,Xv(1+NOVa,i),NOVb)
        if(irpa) call vecadd(tmp3,Yv(1+NOVa,i),NOVb)
        call matmult(Cb(1,1+NOb),Xv(1+NOVa,i),tmp,NBas,NOb,NVb,LM2,NVb,NBas,1)
        call matmult(tmp,Cb(1,1),tmp2,NBas,NBas,NOb,NBas,LM2,NBas,3)
        call vecadd(P(1,i),P(1,i),tmp2,n2)
      else
        call vecscale(P(1,i),n2,2.d0)
      endif
    else
      call matmult(Cb(1,1+NOb),Xv(1,1),tmp,NBas,NOa,NVb,LM2,NVb,NBas,1)
      call matmult(Cb(1,1+NOb),Xv(1,i),tmp2,NBas,NOa,NVb,LM2,NVb,NBas,1)
      call matmult(tmp,tmp2,P(1,i),NBas,NBas,NOa,NBas,NBas,NBas,3)
      call matmult(Ca(1,1),Xv(1,1),tmp,NBas,NVb,NOa,LM2,NVb,NBas,3)
      call matmult(Ca(1,1),Xv(1,i),tmp2,NBas,NVb,NOa,LM2,NVb,NBas,3)
      call matmult(tmp,tmp2,tmp3,NBas,NBas,NVb,NBas,NBas,NBas,3)
      call vecsub(P(1,i),P(1,i),tmp3,n2)
    endif
  enddo

  deallocate(tmp)
  deallocate(tmp2)
  deallocate(tmp3)
  return
end subroutine get_cis_transition_density2

!======================================================== 
subroutine get_cis_density(V,P,Ca,Cb)
  implicit none
  integer i,n2
  real*8  V(*),P(nbas*nbas,*),Ca(LM2,*),Cb(LM2,*)
  real*8, dimension(:), allocatable :: tmp,tmp2,tmp3

  n2 = nbas*nbas
  allocate(tmp(N2),tmp2(N2),tmp3(N2))

  do i=1,nroots
    if(imult.le.2) then
      if(i.gt.1) then
        call matmult(Ca(1,1+NOa),Xv(1,i),tmp,NBas,NOa,NVa,LM2,NVa,NBas,1)
        call matmult(tmp,tmp,P(1,i),NBas,NBas,NOa,NBas,NBas,NBas,3)
        call matmult(Ca(1,1),Xv(1,i),tmp,NBas,NVa,NOa,LM2,NVa,NBas,3)
        call matmult(tmp,tmp,tmp2,NBas,NBas,NVa,NBas,NBas,NBas,3)
        call vecsub(P(1,i),P(1,i),tmp2,NBas*NBas)
        if(irpa) then
          call matmult(Ca(1,1+NOa),Yv(1,i),tmp,NBas,NOa,NVa,LM2,NVa,NBas,1)
          call matmult(tmp,tmp,tmp2,NBas,NBas,NOa,NBas,NBas,NBas,3)
          call vecadd(P(1,i),P(1,i),tmp2,n2)
          call matmult(Ca(1,1),Yv(1,i),tmp,NBas,NVa,NOa,LM2,NVa,NBas,3)
          call matmult(tmp,tmp,tmp2,NBas,NBas,NVa,NBas,NBas,NBas,3)
          call vecsub(P(1,i),P(1,i),tmp2,n2)
        endif
        if(nden.eq.2) then
          call matmult(Cb(1,1+NOb),Xv(1+NOVa,i),tmp,NBas,NOb,NVb,LM2,NVb,NBas,1)
          call matmult(tmp,tmp,tmp2,NBas,NBas,NOb,NBas,NBas,NBas,3)
          call vecadd(P(1,i),P(1,i),tmp2,n2)
          call matmult(Cb(1,1),Xv(1+NOVa,i),tmp,NBas,NVb,NOb,LM2,NVb,NBas,3)
          call matmult(tmp,tmp,tmp2,NBas,NBas,NVb,NBas,NBas,NBas,3)
          call vecsub(P(1,i),P(1,i),tmp2,NBas*NBas)
          if(irpa) then
            call matmult(Cb(1,1+NOb),Yv(1+NOVa,i),tmp,NBas,NOb,NVb,LM2,NVb,NBas,1)
            call matmult(tmp,tmp,tmp2,NBas,NBas,NOb,NBas,NBas,NBas,3)
            call vecadd(P(1,i),P(1,i),tmp2,n2)
            call matmult(Cb(1,1),Yv(1+NOVa,i),tmp,NBas,NVb,NOb,LM2,NVb,NBas,3)
            call matmult(tmp,tmp,tmp2,NBas,NBas,NVb,NBas,NBas,NBas,3)
            call vecsub(P(1,i),P(1,i),tmp2,NBas*NBas)
          endif
        else
          call vecscale(P(1,i),n2,2.d0)
        endif
      endif
    else
      call matmult(Cb(1,1+NOb),Xv(1,i),tmp,NBas,NOa,NVb,LM2,NVb,NBas,1)
      call matmult(tmp,tmp,P(1,i),NBas,NBas,NOa,NBas,NBas,NBas,3)
      call matmult(Ca(1,1),Xv(1,i),tmp,NBas,NVb,NOa,LM2,NVb,NBas,3)
      call matmult(tmp,tmp,tmp3,NBas,NBas,NVb,NBas,NBas,NBas,3)
      call vecsub(P(1,i),P(1,i),tmp3,n2)
    endif
    call square(V(jPa),tmp,NBas,NBas,LM4)
    call vecadd(P(1,i),P(1,i),tmp,n2)
    call square(V(jPb),tmp,NBas,NBas,LM4)
    call vecadd(P(1,i),P(1,i),tmp,n2)
  enddo

  deallocate(tmp)
  deallocate(tmp2)
  deallocate(tmp3)
  return
end subroutine get_cis_density

end module cis_solve
