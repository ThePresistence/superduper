!======================================================== 
!
! Created by Jie Liu on 04/16
! 
! CISEneRGY finds the lowest CIS states  
!
!======================================================== 
module rpa_energy

use es_global
use cis_solve

public  rpaenergy

contains

!======================================================== 
subroutine rpaenergy(V)
  implicit none
  real*8 V(*)
 
  call erpa_init()

  call erpa_solve(V)

  call erpa_exit()

  return
end subroutine rpaenergy

!======================================================== 
subroutine erpa_init()
  if(allocated(At)) deallocate(At)
  if(allocated(Bt)) deallocate(Bt)
  if(allocated(Tv)) deallocate(Tv)
  allocate(At(NOV,maxdav))
  allocate(Bt(NOV,maxdav))
  allocate(Tv(NOV,maxdav))
end subroutine erpa_init

!======================================================== 
subroutine erpa_exit()
  if(allocated(At)) deallocate(At)
  if(allocated(Bt)) deallocate(Bt)
  if(allocated(Tv)) deallocate(Tv)
end subroutine erpa_exit

!======================================================== 
subroutine erpa_solve(V)
  implicit none
  real*8 V(*)
  integer i,m,iroots
  integer, dimension(:), allocatable :: map

  ! save nroots for print out
  call eiter2(V)    
  iroots = nroots         
  allocate(map(nroots*(I13+1)))

  ! running singlet calculation if both singlet and triplet states are required
  if(I13.eq.1) then
    I1n = 1 
    call eiter2(V)

    ! rank of excitation energies
    iroots = iroots + nroots
    call esort(es,map,iroots)
  else
    do m = 1,iroots
      map(m) = m
    enddo
  endif

  call erpasave

  if(iprtci.gt.0) call rpaprt(map,nroots,V)

  deallocate(map)

  return
end subroutine erpa_solve

!======================================================== 
subroutine eiter2(V)
  implicit none
  logical success
  real*8 V(*)
  real*8, dimension(:), allocatable :: H

  nl = nroots
  ns = nroots
  niter = 0
  success = .false.
  allocate(H(NOV))
  call eorbdif(V,H,H(NOVa+1))
  call eiter2init(V,H,H(NOVa+1))

  do while ((.not.success).and.(niter.le.kitdav))
     niter = niter + 1 
     call makea2(V,H) 
     call rpadiag(H)
     if (nl.eq.0) success=.true.
  enddo
  if (.not.success) then
     write(6,*) "Check max Davidson iteration"
     stop
  endif

  deallocate(H)
  return
end subroutine eiter2

subroutine eiter2init(V,Ha,Hb)
  implicit none
  integer ibeta
  real*8  V(*),Ha(*),Hb(*)
  integer, dimension(:), allocatable :: Ia,Ib

  allocate(Ia(NOVa),Ib(NOVb))

  ! Sort up the MO orbitals difference
  call esort(Ha,Ia,NOVa)
  call esort(Hb,Ib,NOVb)

  ! Initionalize the trial vectors
  call eguess(Tv,Ha,Hb,Ia,Ib)
  if(iprtci.ge.5) then
    write(nb6,*) "initial guess"
    call matprnt(Tv,NOV,nroots,6)
  endif

  deallocate(Ia)
  deallocate(Ib)
  return

end subroutine eiter2init
!======================================================== 
subroutine makea2(V,H)
  implicit none
  integer lb
  real*8, dimension(:), allocatable :: Pa,Fa,Ha,Hb
  real*8 V(*),H(*)

  allocate(Pa(NBas*NBas*nl*nden),Fa(NBas*NBas*nl*nden))
  allocate(Ha(NOVa),Hb(NOVb))

  lb = NBas*NBas*nl*(nden-1) + 1     

  ! Transition density transformed from MO to AO
  call emake_rao(Pa,Pa(lb),Tv(1,ns-nl+1),V(jCa),V(jCb),nl)

  ! Fock-like matrix in AO
  call makefao(Fa,Fa(lb),Pa,Pa(lb),V(jW),nl,NBas,LM6,nden,I1n,cfacJ,cfacK)
  deallocate(Pa)

  ! Building A matrix in MO
  !call eorbdif(V,Ha,Hb)
  call emake_amo(At(1,ns-nl+1),Fa,Fa(lb),Tv(1,ns-nl+1),H,H(NOVa+1),V(jCa),V(jCb),nl,.false.,imult,nden)
  call emake_bmo(Bt(1,ns-nl+1),Fa,Fa(lb),Tv(1,ns-nl+1),V(jCa),V(jCb),nl)

  deallocate(Fa)
  deallocate(Ha)
  deallocate(Hb)
  return
end subroutine makea2


!======================================================== 
subroutine emake_bmo(Az,Fa,Fb,R,Ca,Cb,n)
  implicit none
  integer i,m,n
  real*8 Az(NOV,*),Fa(NBas*NBas,*),Fb(NBas*NBas,*),R(NOV,*)
  real*8 Ca(LM2,*),Cb(LM2,*)
  real*8, dimension(:), allocatable :: tmp

  allocate(tmp(NBas*NBas))

  do m=1,n
    if(imult.eq.0 .or. imult.eq.2) then
      ! A = Cv^T * F * Co
      call matmult(Ca(1,NOa+1),Fa(1,m),tmp,NVa,NBas,NBas,LM2,NBas,NVa,4)
      call matmult(tmp,Ca(1,1),Az(1,m),NVa,NOa,NBas,NVa,LM2,NVa,1)
      if(nden.eq.2) then
        call matmult(Cb(1,NOb+1),Fb(1,m),tmp,NVb,NBas,NBas,LM2,NBas,NVb,4)
        call matmult(tmp,Cb(1,1),Az(1+NOVa,m),NVb,NOb,NBas,NVb,LM2,NVb,1)
      endif
    else
      write(6,*) "wrong multiplicity"
      stop
    !else if(imult.eq.3) then
      ! A = Cbv^T * F * Cao
      !call matmult(Cb(1,NOb+1),Fa(1,m),tmp,NVb,NBas,NBas,LM2,NBas,NVb,4)
      !call matmult(tmp,Ca(1,1),Az(1,m),NVb,NOa,NBas,NVb,LM2,NVb,1)
    endif

    if(imult.eq.0) then
      if(nden.eq.1) then
        call veccopy(Az(1+NOVa,m),Az(1,m),NOVa)
        if(I1n.eq.0) call vecscale(Az(1+NOVa,m),NOVa,-1d0)
      endif
    endif
  enddo
  if(iprtci.ge.5) then
    write(nb6,*) "Bz"
    call matprnt(Az,NOV,n,6)
  endif

  deallocate(tmp)
  return
end subroutine emake_bmo

!======================================================== 
subroutine rpadiag(H)
  implicit none
  integer ierror,i
  real*8  H(*), dot1, dot2
  real*8, dimension(:), allocatable :: sAs,sBs,E
  real*8, dimension(:,:), allocatable :: tmp,tmp2,tmp3,tmp4

  allocate(sAs(ns*ns),sBs(ns*ns),E(ns))
  allocate(tmp(ns,ns),tmp2(ns,ns),tmp3(ns,ns),tmp4(ns,ns))

  call matmult(Tv,At,sAs,ns,ns,NOV,NOV,NOV,ns,2)
  call matmult(Tv,Bt,sBs,ns,ns,NOV,NOV,NOV,ns,2)
  if(iprtci.ge.5) then
    write(nb6,*) "A in subspace"
    call matprnt(sAs,ns,ns,6)
    write(nb6,*) "B in subspace"
    call matprnt(sBs,ns,ns,6)
  endif
  call vecsub(tmp,sAs,sBs,ns*ns)
  call vecadd(sAs,sAs,sBs,ns*ns)
  call veccopy(sBs,tmp,ns*ns)
  

  ! tmp2 = (A-B)^1/2, tmp4 = (A-B)^-1
  call matdiagsquare(tmp,E,ns,ierror)
  if(ierror.ne.0) then
    write(nb6,*) "A-B diagonalization fails"
    stop
  endif
  call vecinit(tmp2,ns*ns,0.d0)
  call vecinit(tmp4,ns*ns,0.d0)
  do i=1,ns
    if(E(i).lt.0.d0) then
      write(nb6,*) "eigenvalue of A-B is negative"
      stop
    endif
    tmp2(i,i) = sqrt(E(i))
    tmp4(i,i) = 1/E(i)
  enddo
  call matmult(tmp2,tmp,tmp3,ns,ns,ns,ns,ns,ns,3)
  call matmult(tmp,tmp3,tmp2,ns,ns,ns,ns,ns,ns,1)
  call matmult(tmp4,tmp,tmp3,ns,ns,ns,ns,ns,ns,3)
  call matmult(tmp,tmp3,tmp4,ns,ns,ns,ns,ns,ns,1)
  
  ! tmp3 = (A-B)^1/2 (A+B) (A-B)^1/2
  call matmult(tmp2,sAs,tmp,ns,ns,ns,ns,ns,ns,1)
  call matmult(tmp,tmp2,tmp3,ns,ns,ns,ns,ns,ns,1)
  if(iprtci.ge.5) then
    write(nb6,*) "(A-B)^1/2 (A+B) (A-B)^1/2"
    call matprnt(tmp3,ns,ns,6)
  endif
  call matdiagsquare(tmp3,E,ns,ierror)
  if(ierror.ne.0) then
    write(nb6,*) "(A-B)^1/2 (A+B) (A-B)^1/2 diagonalization fails"
    stop
  endif
  do i=1,ns
     if(E(i).lt.0.d0) then
      write(nb6,*) "imaginary omega^2"
      stop
    endif
    E(i) = sqrt(E(i))
  enddo
  if(iprtci.ge.5) then
    write(nb6,*) "excitation energy:" 
    call matprnt(E,ns,1,6)
  endif

  ! tmp = (X+Y)
  call matmult(tmp2,tmp3,tmp,ns,ns,ns,ns,ns,ns,1)
  
  ! tmp4 = (X-Y)
  call matmult(tmp4,tmp,tmp2,ns,ns,ns,ns,ns,ns,1)
  call vecinit(tmp3,ns*ns,0.d0)
  do i=1,ns
    tmp3(i,i) = E(i)
  enddo
  call matmult(tmp2,tmp3,tmp4,ns,ns,ns,ns,ns,ns,1)

  
  ! sAs = X, sBs = Y
  call vecadd(sAs,tmp,tmp4,ns*ns)
  call vecsub(sBs,tmp,tmp4,ns*ns)
  do i=1,ns
    call vecdot(dot1,sAs(1+(i-1)*ns),sAs(1+(i-1)*ns),ns)
    call vecdot(dot2,sBs(1+(i-1)*ns),sBs(1+(i-1)*ns),ns)
    dot1 = 1.d0 / sqrt(dot1 - dot2)
    call vecscale(sAs(1+(i-1)*ns),ns,dot1)
    call vecscale(sBs(1+(i-1)*ns),ns,dot1)
  enddo

  ! check for convergence
  call chkcnv2(sAs,sBs,E,H)

  deallocate(sAs)
  deallocate(E)
  deallocate(tmp)
  deallocate(tmp2)
  deallocate(tmp3)
  deallocate(tmp4)
  return
end subroutine rpadiag

!======================================================== 
subroutine chkcnv2(Xs,Ys,E,H)
  implicit real*8 (a-h,o-z)
  integer m, isub, ioff
  real*8  tol, maxd, d
  real*8  Xs(*), Ys(*), E(*), H(*)
  real*8, dimension(:), allocatable :: Res

  allocate(Res(NOV))

  tol   = 1d1**(1-nrmdav)
  maxd  = 0d0
  isub  = ns

  !call eorbdif(V,H,H(NOVa+1))

  ! Check convergence and expand subspace
  do m = 1,nroots
     call chkres2(Res,Xs(ns*(m-1)+1),Ys(ns*(m-1)+1),E(m),d,H)
     if(d.gt.maxd) maxd = d
     if(d.gt.tol)  call ssexpand(Res,isub,NOV)
  enddo
  if(iprtci.ge.5) then
    write(nb6,*) "expanded subspace",ns,isub
    call matprnt(Tv(1,ns+1),NOV,isub-ns,6)
  endif

  ! Save the information for subspace
  nl = isub - ns
  write(nb6,*) niter,nl,maxD
  ns = isub

  if(nl.eq.0) then
    if(.not.allocated(Xv)) allocate(Xv(NOV,nroots*(I13+1)))
    if(.not.allocated(Yv)) allocate(Yv(NOV,nroots*(I13+1)))
    if(.not.allocated(es)) allocate(es(nroots,I13+1))
    ioff = 0
    if(I3.and.I1.and.I1n) ioff = nroots
    call matmult(Tv,Xs,Xv(1,ioff+1),NOV,nroots,ns,NOV,ns,NOV,1)
    call matmult(Tv,Ys,Yv(1,ioff+1),NOV,nroots,ns,NOV,ns,NOV,1)
    es(1:nroots,I13+1) = E(1:nroots)
    if(iprtci.gt.5) then
      write(nb6,*) "rpa X amplitude"
      call matprnt(Xv,NOV,nroots,6)
      write(nb6,*) "rpa Y amplitude"
      call matprnt(Yv,NOV,nroots,6)
    endif
    write(NB6,*) "RPA Converged"
  else if(ns.gt.maxdav) then
    ! If maximum dimension of Davidson subspace is reached,
    ! the subspace information is saved to restart rpa calculation
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
end subroutine chkcnv2

!======================================================== 
subroutine chkres2(Res,Xs,Ys,E,D,H)
  implicit none
  integer i, j
  real*8  D, D2, E, Res(*), Xs(*), Ys(*), H(*)
  real*8, dimension(:), allocatable :: tmp,tmp2

  allocate(tmp(NOV),tmp2(NOV))

  call vecinit(tmp,NOV,0d0)
  call vecinit(tmp2,NOV,0d0)
  do i = 1,ns
     do j = 1,NOV
        tmp(j)  = tmp(j)  + (At(j,i)-E*Tv(j,i))*Xs(i) + Bt(j,i)*Ys(i)
        tmp2(j) = tmp2(j) + (At(j,i)+E*Tv(j,i))*Ys(i) + Bt(j,i)*Xs(i)
     enddo
  enddo

  call vecdot(d,tmp,tmp,NOV)
  call vecdot(d2,tmp2,tmp2,NOV)
  d = sqrt(d/NOV)
  d2= sqrt(d2/NOV)

  if(d.ge.d2) then
    do i=1,NOV
       Res(i) = tmp(i)/(E-H(i))
    enddo
  else
    d = d2
    do i=1,NOV
       Res(i) = tmp2(i)/(E-H(i))
    enddo

  endif

  deallocate(tmp)
  deallocate(tmp2)
  return
end subroutine chkres2

subroutine rpaprt(map,iroots,V)
  use limit, only : lmprop, lmstat
  use cissym
  implicit none
  integer i, j, k, m, iroots, mult, ioff, isym
  real*8  V(*), dmax, escf, enuc, au2dby, ciprop
  real*8, dimension(:,:), allocatable :: tdip,P,R
  integer map(*)
  common  /energt/ escf,enuc
  common  /ciprp / ciprop(lmprop,lmstat)

  ! determine the symmetry labels for all active orbitals
  call ovsym(V)

  ! calculate the oscillator strengths and transition moments
  allocate(tdip(4,iroots))
  call vecinit(tdip,4*iroots,0.d0)
  if(iuvcd.gt.0) then
    allocate(P(nbas*nbas,iroots),R(nbas*nbas,iroots))
    call get_cis_density(V,P,V(jCa),V(jCb))
    call get_cis_transition_density2(R,V(jCa),V(jCb))
    call cisprop(P,R,Es,ciprop,jsym,labels,jsym,iroots,LM3)
    do i=1,iroots
      tdip(1,i) = ciprop(5,i)
    enddo
    deallocate(P)
    deallocate(R)
  endif

  au2dby = 0.529167D0*4.803242D0
  ! print out
  do m = 1,iroots
    ! Spin multiplicity
    If(I1.eq.1.and.I3.eq.1) Then
       mult = 1
       if(map(m).le.(iroots/2)) mult = 3
    else if(I3.eq.1) Then
       mult = 3
    else if(I1.eq.1) Then
       mult = 1
    else if(IMult.eq.2) Then
       mult = 2
    endif

    if(is_sym) then
      write(nb6,'(A/)') '-------------------------------------------------------------------------'
      write(nb6,'(" State",I3,",  Mult.",I2,",  ",A4,"(",I1,"), E-E(1)=",F10.6," eV,  E=",F14.6," eV")') &
            m, mult, labels(jsym(m)), jsym(m), es(map(m),1), es(map(m),1)+escf+enuc
      write(nb6,'(A/)') '-------------------------------------------------------------------------'
    else
      write(nb6,'(A/)') '-------------------------------------------------------------------------'
      write(nb6,'(" State",I3,",  Mult.",I2,", E-E(1)=",F10.6," eV,  E=",F14.6," eV")') &
            m, mult, es(map(m),1), es(map(m),1)+escf+enuc
      write(nb6,'(A/)') '-------------------------------------------------------------------------'
    endif

    write(*,'(" Px: ",F6.2,"     Py: ",F6.2,"     Pz: ",F6.2)') &
      tdip(2,map(m))*au2dby, tdip(3,map(m))*au2dby, tdip(4,map(m))*au2dby
    write(*,'(" Strength:      ",F6.2)') tdip(1,map(m))

    ioff = 0
    do j = 1,NOa
       do k = 1,NVa
          ioff = ioff + 1 
          if(abs(Xv(ioff,map(m))).ge.0.1) then
             If(nden.eq.1) Then
               write(*,'("X: ",I3," -->",I3,F10.3)') j,k+NOa, &
                     sqrt(2.0d0)*Xv(ioff,map(m))
             else
               write(*,'("X: ",I3," -->",I3,F10.3,"  Alpha")') j,k+NOa, &
                     sqrt(2.0d0)*Xv(ioff,map(m))
              endif
           endif
           if(abs(Yv(ioff,map(m))).ge.0.1) then
             If(nden.eq.1) Then
               write(*,'("Y: ",I3," -->",I3,F10.3)') j,k+NOa, &
                     sqrt(2.0d0)*Yv(ioff,map(m))
             else
               write(*,'("Y: ",I3," -->",I3,F10.3,"  Alpha")') j,k+NOa, &
                     sqrt(2.0d0)*Yv(ioff,map(m))
              endif
           endif
        enddo
     enddo
  enddo

  deallocate(tdip)
  return
end subroutine rpaprt

!======================================================== 
subroutine erpasave
  use limit, only : lmprop, lmstat
  implicit real*8 (a-h,o-z)
  common  /mmcom1/ emm
  common  /mmdp  / emmdp
  common  /ciprp / ciprop(lmprop,lmstat)
  common  /constf/ a0,afact,ev,evcal
  common  /inopt2/ in2(300)
  common  /energt/ e0,enuc,eat,atheat

  nroots = in2(139)

  if(.not.is_sf) then
    do i=nroots,2,-1
      es(i,1) = es(i-1,1)
    enddo
    es(1,1) = 0.d0
  endif

  ! Save data for dynamics. Energy expression depends on iaterg.
  if(in2(119).eq.-1) then
     ! New approach for ODM2 and ODM3.
     ciprop(1,1:nroots) = (es(1:nroots,1)+e0+enuc)*evcal+emm+emmdp
  else
     ! Traditional approach for all other methods.
     ciprop(1,1:nroots) = (es(1:nroots,1)+e0+enuc-eat)*evcal+atheat+emm+emmdp
  endif

  ! save the excited state energy specified by lroot
  etarget = es(lstate,1)

  ! save cis amplitude for state following
  if(stat_fol .and. allocated(Xn)) then
    call veccopy(Xn,Xv,NOV*in2(139))
  endif

  if(.not.is_sf) then
    do i=nroots,2,-1
      call veccopy(Xv(1,i),Xv(1,i-1),NOV)
    enddo
    call vecinit(Xv(1,1),NOV,0.d0)
  endif

  return
end subroutine erpasave

end module rpa_energy
