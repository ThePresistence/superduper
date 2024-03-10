!======================================================== 
!
! Created by Jie Liu on 04/16
! 
! cisenergy finds the lowest CIS states  
!
!======================================================== 
module cis_energy

use es_global
use cis_solve

public  cisenergy

contains

!======================================================== 
subroutine cisenergy(V)
  real*8 V(*)
 
  call ecisinit()

  call ecissolve(V)

  call ecisexit()

  return
end subroutine cisenergy

!======================================================== 
subroutine ecisinit()
  if(allocated(At)) deallocate(At)
  if(allocated(Tv)) deallocate(Tv)
  allocate(At(NOV,maxdav))
  allocate(Tv(NOV,maxdav))
end subroutine ecisinit

!======================================================== 
subroutine ecissave
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
end subroutine ecissave

!======================================================== 
subroutine ecisexit()
  if(allocated(At)) deallocate(At)
  if(allocated(Tv)) deallocate(Tv)
end subroutine ecisexit


!======================================================== 
subroutine ecissolve(V)
  implicit none
  integer, dimension(:), allocatable :: map
  integer i,m,iroots
  real*8 V(*)

  call eiter(V)    
  iroots = nroots         
  allocate(map(nroots*(I13+1)))

  ! running singlet calculation if both singlet and triplet states are required
  if(I13.eq.1) then
    I1n = 1 
    call eiter(V)
    ! rank of excitation energies
    iroots = iroots + nroots
    call esort(es,map,iroots)
  else
    do m = 1,iroots
      map(m) = m
    enddo
  endif

  call ecissave

  if(iprtci.ge.0) call cisprt(map,nroots,V)

!  if(iprtci.ge.2) call cisana(V(jCa),V(jCb))

  deallocate(map)
  return
end subroutine ecissolve

!======================================================== 
subroutine eiterinit(V,Ha,Hb)
  implicit none
  integer i
  real*8  V(*),Ha(*),Hb(*)
  integer, dimension(:), allocatable :: Ia,Ib

  allocate(Ia(NOV),Ib(NOVb))

  ! Sort up the MO orbitals difference
  call esort(Ha,Ia,NOVa)
  call esort(Hb,Ib,NOVb)

  ! Initionalize the trial vectors
  if(gread) then
    call veccopy(Tv,Xn,NOV*nroots)
  else
    call eguess(Tv,Ha,Hb,Ia,Ib)
  endif
  if(iprtci.ge.5) then
    write(nb6,*) "initial guess"
    call matprnt(Tv,NOV,nroots,6)
  endif

  deallocate(Ia)
  deallocate(Ib)
  return
end subroutine eiterinit

!======================================================== 
subroutine eiter(V)
  implicit none
  logical success
  real*8  V(*)
  real*8, dimension(:), allocatable :: H

  nl = nroots
  ns = nroots
  niter = 0
  success = .false.
  allocate(H(NOV))
  call eorbdif(V,H,H(NOVa+1))
  call eiterinit(V,H,H(NOVa+1))

  do while ((.not.success).and.(niter.le.kitdav))
     niter = niter + 1 
     call makea(V,H) 
     call cisdiag(H)
     if (nl.eq.0) success = .true.
  enddo
  if (.not.success) then
     write(6,*) "Check Maximum Davidson Iteration"
     stop
  endif

  deallocate(H)
  return
end subroutine eiter

!======================================================== 
subroutine makea(V,H)
  implicit none
  integer lb
  real*8  V(*),H(*)
  real*8, dimension(:), allocatable :: Pa,Fa

  allocate(Pa(NBas*NBas*nl*nden),Fa(NBas*NBas*nl*nden))

  lb = NBas*NBas*nl*(nden-1) + 1     

  ! Transition density transformed from MO to AO
  call emake_rao(Pa,Pa(lb),Tv(1,ns-nl+1),V(jCa),V(jCb),nl)

  ! Fock-like matrix in AO
  call makefao(Fa,Fa(lb),Pa,Pa(lb),V(jW),nl,NBas,LM6,nden,I1n,cfacJ,cfacK)
  deallocate(Pa)

  ! Building A matrix in MO
  call emake_amo(At(1,ns-nl+1),Fa,Fa(lb),Tv(1,ns-nl+1),H,H(NOVa+1),&
                 V(jCa),V(jCb),nl,.false.,imult,nden)

  deallocate(Fa)
  return
end subroutine makea

!======================================================== 
subroutine cisprt(map,iroots,V)
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

    ! print out the primary single excitations
    ioff = 0
    if(is_sf) then
    do j=1,NOa
      do k=1,NVb
        ioff = ioff + 1
        if(abs(Xv(ioff,map(m))).ge.0.1) then
          write(*,'(I3,"(alpha) -->",I3,"(beta)",F8.3)') j+NAlpha-NOa,k+NOb, &
                  Xv(ioff,map(m))
        endif
      enddo !j
    enddo !k
    else
    do j=1,NOa
      do k=1,NVa
        ioff = ioff + 1 
        if(abs(Xv(ioff,map(m))).ge.0.1) then
          if(nden.eq.1) Then
            write(*,'(I3," -->",I3,F8.3)') j+NAlpha-NOa,k+NAlpha, &
                  sqrt(2.0d0)*Xv(ioff,map(m))
          else
            write(*,'(I3," -->",I3,F8.3,"  Alpha")') j+NAlpha-NOa,k+NAlpha,Xv(ioff,map(m))
          endif
        endif
      enddo !j
    enddo !k
    if(nden.eq.2) then
    do j=1,NOb
      do k=1,NVb
        ioff = ioff + 1
        if(abs(Xv(ioff,map(m))).ge.0.1) then
          write(*,'(I3," -->",I3,F8.3,"  Beta")') j+NBeta-NOb,k+NBeta,Xv(ioff,map(m))
        endif
      enddo !j
    enddo !k    
    endif
    endif
  enddo

  deallocate(tdip)
  return
end subroutine cisprt

end module cis_energy
