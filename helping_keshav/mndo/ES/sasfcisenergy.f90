!==================================================================== 
!
! Created by Jie Liu on 04/16
! 
! sasfcisenergy finds the lowest spin-adapted spin-flip CIS states  
!
!====================================================================
module sasfcis_energy

use es_global
use cis_solve
use sasfcis_util

implicit none

real*8  et3
real*8, dimension(:), allocatable :: Fa, Fb

public  sasfcisenergy

contains

!======================================================== 
subroutine sasfcisenergy(V)
  implicit none
  real*8 V(*)

  call sasfecisinit

  call sasfecissolve(V)

  call sasfecisexit

  return
end subroutine sasfcisenergy

!======================================================== 
subroutine sasfecisinit

  if(iprtci.gt.1) then
    if(sf_xcis) then
      write(nb6,*) "**********  spin adapted spin-flipped xcis  *************"
    else
      write(nb6,*) "**********  spin adapted spin-flipped cis  *************"
    endif
    if(iprtci.gt.2) then
      write(nb6,*) "alpha occupied orbital: ",NOa
      write(nb6,*) "alpha virtual  orbital: ",NVa
      write(nb6,*) "alpha occupied orbital: ",NOb
      write(nb6,*) "beta  virtual  orbital: ",NVb 
    endif
  endif

  
  call sasfcis_setup

  allocate(At(NOV,maxdav))
  allocate(Tv(NOV,maxdav))

end subroutine sasfecisinit

!======================================================== 
subroutine sasfecisexit()
  use limit, only : lmprop, lmstat
  implicit real*8 (a-h,o-z)
  integer in2
  common  /mmcom1/ emm
  common  /mmdp  / emmdp
  common  /ciprp / ciprop(lmprop,lmstat)
  common  /constf/ a0,afact,ev,evcal
  common  /inopt2/ in2(300)
  common  /energt/ e0,enuc,eat,atheat

  if(stat_fol .and. allocated(Xn)) then
    call veccopy(Xn,Xv,NOV*in2(139))
  endif

  ! Save data for dynamics. Energy expression depends on iaterg.
  if(in2(119).eq.-1) then
     ! New approach for ODM2 and ODM3.
     ciprop(1,1:nroots) = (es(1:nroots,1)+e0+et3+enuc)*evcal+emm+emmdp
  else
     ! Traditional approach for all other methods.
     ciprop(1,1:nroots) = (es(1:nroots,1)+e0+et3+enuc-eat)*evcal+atheat+emm+emmdp
  endif 

  es(1:nroots,1) = es(1:nroots,1) + et3

  ! save the excited state energy specified by lroot
  etarget = es(lstate,1)
  
  if(allocated(At)) deallocate(At)
  if(allocated(Tv)) deallocate(Tv)

  return
end subroutine sasfecisexit

!======================================================== 
subroutine sasfecissolve(V)
  implicit none
  integer m
  real*8  V(*)

  if(I3.eq.1) then
    write(nb6,*) "triplet state not available in sasfcis"
    stop
  endif

  call sasfeiter(V)    

  call sasfcisprt(V,nroots)

  return
end subroutine sasfecissolve

!======================================================== 
subroutine sasfeiter(V)
  implicit none
  integer i,j
  logical success
  real*8  V(*)
  real*8, dimension(:), allocatable :: H

  nl = nroots
  ns = nroots
  niter = 0
  success = .false.
  allocate(H(NOV))
  allocate(Fa(n2),Fb(n2))
  allocate(Js(n2,3),Ks(n2,3))
  call sasfeiterinit(V,H)
  call sasfeguess(H)

  do while ((.not.success).and.(niter.le.kitdav))
     niter = niter + 1
     call sasfemakea(V)
     call cisdiag(H)
     if (nl.eq.0) success=.true.
  enddo

  if (.not.success) then
     write(6,*) "Check max Davidson iteration"
     stop
  endif
  deallocate(H)
  deallocate(Js)
  deallocate(Ks)
  deallocate(Fa)
  deallocate(Fb)

  return
end subroutine sasfeiter

!======================================================== 
subroutine sasfeguess(H)
  implicit none
  real*8 H(*)
  integer, dimension(:), allocatable :: I

  allocate(I(NOV))

  ! Sort up the MO orbitals difference
  call esort(H,I,NOV)

  ! Initionalize the trial vectors
  if(.not.gread) then
    call eguess(Tv,H,H,I,I)
  else
    call veccopy(Tv,Xn,NOV*nroots)
  endif
  if(iprtci.eq.5) then
    write(nb6,*) "initial guess"
    call matprnt(Tv,NOV,nroots,6)
  endif

  deallocate(I)
  return
end subroutine sasfeguess

!======================================================== 
subroutine sasfeiterinit(V,H)
  implicit none
  integer i,j
  real*8  V(*),H(*)
  real*8, dimension(:), allocatable :: Pa,tmp,tmp2

  allocate(Pa(n2*3))
  ! Transition density transformed from MO to AO
  call sasfemake_sao(Pa,V(jCa))
  ! Fock-like matrix in AO
  call makejk(Js,Js,Ks,Ks,Pa,Pa,V(jW),3,nbas,LM6,nden,I1n,cfacJ,cfacK)
  deallocate(Pa)
  ! AO->MO
  call sasfemake_smo(Js,Ks,V(jCa),3)
  
  ! Fock matrix in MO basis
  allocate(tmp(n2),tmp2(n2))
  call square(V(jFa),tmp,nbas,nbas,LM4)
  call matmult(V(jCa),tmp,tmp2,nbas,nbas,nbas,LM2,nbas,nbas,2)
  call matmult(tmp2,V(jCa),Fa,nbas,nbas,nbas,nbas,LM2,nbas,1)
  call veccopy(Fb,Fa,nbas*nbas)
  et3 = 0.d0
  if(sref) then
    et3 = Fa(s2s2+1) - Fa(s1s1+1) - Js(s1s1+1,3)
    call vecadd3(Fa,-one,Js(1,1),nbas*nbas)
    call vecadd3(Fa,one,Js(1,3),nbas*nbas)
    call vecadd3(Fa,one,Ks(1,3),nbas*nbas)
    call vecadd3(Fb,-one,Js(1,1),nbas*nbas)
    call vecadd3(Fb,-one,Ks(1,1),nbas*nbas)
    call vecadd3(Fb,one,Js(1,3),nbas*nbas)
  else
    call vecadd3(Fa,f12,Ks(1,1),nbas*nbas)
    call vecadd3(Fa,f12,Ks(1,3),nbas*nbas)
    call vecadd3(Fb,-f12,Ks(1,1),nbas*nbas)
    call vecadd3(Fb,-f12,Ks(1,3),nbas*nbas)
  endif
  if(iprtci.ge.5) then
    write(nb6,*) "alpha fock matrix in mo basis"
    call matprnt(Fa,nbas,nbas,6)
    write(nb6,*) "beta fock matrix in mo basis"
    call matprnt(Fb,nbas,nbas,6)
    write(nb6,*) "Js(s1,s1)"
    call matprnt(Ks,nbas,nbas,6)
    write(nb6,*) "Js(s2,s2)"
    call matprnt(Ks(1,3),nbas,nbas,6)
  endif
  deallocate(tmp)
  deallocate(tmp2)

  ! diagnol elements
  H(1) = Fb(s2s2+1) - Fa(s1s1+1) !- Js(s2s2+1,1) !V(jEa+NOb+1) - V(jEa+NOb)
  H(2) = Fb(s1s1+1) - Fa(s2s2+1) !- Js(s2s2+1,1) !V(jEa+NOb)   - V(jEa+NOb+1)
  H(3) = Fb(s1s1+1) - Fa(s1s1+1) !- Js(s1s1+1,1) + Js(s1s2+1,2) !0.d0
  do i=0,NOb-1
    H(os1+i) = Fb(s1s1+1) - Fa(i+1+i*nbas) !- Js(i+1+i*nbas,1) - Ks(i+1+i*nbas,3) !V(jEa+NOb)  -V(jEa+i)
    H(os2+i) = Fb(s2s2+1) - Fa(i+1+i*nbas) !- Ks(i+1+i*nbas,1) - Js(i+1+i*nbas,3) !V(jEa+NOb+1)-V(jEa+i)
  enddo
  do i=0,NVa-1
    H(vs1+i) = Fb(jv+i+(jv+i-1)*nbas) - Fa(s1s1+1) !- Js(jv+i+(jv+i-1)*nbas,1) - Ks(jv+i+(jv+i-1)*nbas,3) !(V(jEa+NOa+i)-V(jEa+NOb)
    H(vs2+i) = Fb(jv+i+(jv+i-1)*nbas) - Fa(s2s2+1) !- Ks(jv+i+(jv+i-1)*nbas,1) - Js(jv+i+(jv+i-1)*nbas,3) !(V(jEa+NOa+i)-V(jEa+NOb+1)
  enddo
  do i=0,NOb-1
    do j=0,NVa-1
      H(ov1+i*NVa+j) = Fb(jv+j+(jv+j-1)*nbas) - Fa(i+1+i*nbas) 
                    ! + Ks(jv+i+(jv+i-1)*nbas,3) - 0.5d0*Ks(jv+i+(jv+i-1)*nbas,1) &
                    ! + Ks(i+1+i*nbas,1) - 0.5d0*Ks(i+1+i*nbas,3) &
                    ! - 1.5d0*Ks(s2s2+1,1) !V(jEa+NOa+j) - V(jEa+i)
      H(ov2+i*NVa+j) = Fb(jv+j+(jv+j-1)*nbas) - Fa(i+1+i*nbas) 
                    ! + 0.5d0*Ks(jv+j+(jv+j-1)*nbas,1) - Ks(jv+j+(jv+j-1)*nbas,3) &
                    ! - Ks(i+1+i*nbas,1) + 0.5d0*Ks(i+1+i*nbas,3) &
                    ! - 0.5d0*Ks(s2s2+1,1)!V(jEa+NOa+j) - V(jEa+i)
      if(sf_xcis) then
        H(ov3+i*NVa+j) = Fa(jv+j+(jv+j-1)*nbas) - Fa(i+1+i*nbas) + Fb(s2s2+1) - Fa(s1s1+1) &
                       + Js(jv+j+(jv+j-1)*nbas,3) - Js(jv+j+(jv+j-1)*nbas,1) - Ks(jv+j+(jv+j-1)*nbas,1) &
                       - Js(i+1+i*nbas,3) + Js(i+1+i*nbas,1) + Js(i+1+i*nbas,1) - Js(s2s2+1,1)
        H(ov4+i*NVa+j) = Fa(jv+j+(jv+j-1)*nbas) - Fa(i+1+i*nbas) + Fb(s1s1+1) - Fa(s2s2+1) &
                       + Js(jv+j+(jv+j-1)*nbas,1) - Js(jv+j+(jv+j-1)*nbas,3) - Ks(jv+j+(jv+j-1)*nbas,3) &
                       - Js(i+1+i*nbas,1) + Js(i+1+i*nbas,3) + Js(i+1+i*nbas,3) - Js(s2s2+1,1)
      endif
    enddo
  enddo

  if(iprtci.ge.5) then
    write(nb6,*) "energy difference"
    do i=1,NOV
      write(nb6,*) i,H(i)
    enddo 
  endif

  return
end subroutine sasfeiterinit

!======================================================== 
subroutine sasfemakea(V)
  implicit none
  integer i,j
  real*8  V(*),time1,time2
  real*8, dimension(:), allocatable :: Pa,Ja,Ka

  allocate(Pa(nbas*nbas*npr*nl),Ja(nbas*nbas*npr*nl),Ka(nbas*nbas*npr*nl))
  call sasfemake_rao(Pa,Tv(1,ns-nl+1),V(jCa),nl)
  call makejk(Ja,Ja,Ka,Ka,Pa,Pa,V(jW),npr*nl,nbas,LM6,nden,I1n,cfacJ,cfacK)
  deallocate(Pa)
  ! Building A matrix in MO
  call sasfemake_amo(At(1,ns-nl+1),Ja,Ka,Tv(1,ns-nl+1),V(jCa),nl)
  deallocate(Ja)
  deallocate(Ka)

  return
end subroutine sasfemakea

!======================================================== 
subroutine sasfemake_amo(Az,Ja,Ka,R,Ca,n)
  implicit none
  integer i,m,n,ioff1,ioff2
  real*8  dot1,dot2,dot3,dot4
  real*8  Az(NOV,*),Ja(nbas*nbas*npr,*),Ka(nbas*nbas*npr,*)
  real*8  R(NOV,*),Ca(LM2,*)
  real*8, dimension(:), allocatable :: tmp, tmp2

  allocate(tmp(nbas*nbas),tmp2(nbas*nbas))
  
  call vecinit(Az,NOV*n,0.d0)

  do i=1,n
    call ao2mo2(Ja(1,i),Ca,npr)
    call ao2mo2(Ka(1,i),Ca,npr)

    if(I1.eq.1) then

    if(sf_xcis) then
      call vecadd3(Ka(jJov3,i),two,Ja(jJov3,i),NVb*NOa)
      call vecadd3(Ka(jJov4,i),two,Ja(jJov4,i),NVb*NOa)
    endif

    ! s1->s2
    Az(1,i) = -Js(s2s2+1,1)*R(1,i)-Js(s1s2+1,2)*R(2,i)-sq2*Js(s1s2+1,1)*R(3,i)   &
              +sq2*Ka(jJos1+s21,i)+sq2*Ka(jJos2+s21,i)                           &
              +sq2*Ka(jJvs1+s21,i)+sq2*Ka(jJvs2+s21,i)                           &
              +sq3*Ja(jJov1+s21,i)+Ja(jJov2+s21,i)+two*Ka(jJov2+s21,i) 
    ! Fock
    Az(1,i) = Az(1,i) + (Fb(s2s2+1)-Fa(s1s1+1))*R(1,i)
    Az(1,i) = Az(1,i) + sq2*Fb(s2s1+1)*R(3,i)
    call vecdot(dot1,R(os2,i),Fa(i1s1+1),NOb)
    call vecdot(dot2,R(vs1,i),Fb(a1s2+1),NVa)
    Az(1,i) = Az(1,i) - sq2*(dot1-dot2)
    if(sf_xcis) then
      Az(1,i) = Az(1,i) + ssq2*(Ka(jJov3+s22,i)-Ka(jJov3+s11,i))
      do m=0,NOb-1
        call vecdot(dot1,R(ov3+m*NVa,i),Fa(jv+m*nbas),NVa)
        call vecdot(dot2,R(ov3+m*NVa,i),Fb(jv+m*nbas),NVa)
        Az(1,i) = Az(1,i) + ssq2*(dot1+dot2)
      enddo
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    ! s2->s1
    Az(2,i) = -Js(s1s2+1,2)*R(1,i)-Js(s2s2+1,1)*R(2,i)-sq2*Js(s1s2+1,1)*R(3,i)     &
              +sq2*Ka(jJos1+s12,i)+sq2*Ka(jJos2+s12,i)                             &
              +sq2*Ka(jJvs1+s12,i)+sq2*Ka(jJvs2+s12,i)                             &
              -sq3*Ja(jJov1+s12,i)-sq3*Ka(jJov1+s12,i)                             &
              -Ja(jJov2+s12,i)+Ka(jJov2+s12,i)
    ! Fock
    Az(2,i) = Az(2,i) + (Fb(s1s1+1)-Fa(s2s2+1))*R(2,i)
    Az(2,i) = Az(2,i) - sq2*Fa(s2s1+1)*R(3,i)
    call vecdot(dot1,R(os1,i),Fa(i1s2+1),NOb)
    call vecdot(dot2,R(vs2,i),Fb(a1s1+1),NVa)
    Az(2,i) = Az(2,i) - sq2*(dot1-dot2)                  
    if(sf_xcis) then
    Az(2,i) = Az(2,i) + ssq2*(Ka(jJov4+s11,i)-Ka(jJov4+s22,i))
    do m=0,NOb-1
      call vecdot(dot1,R(ov4+m*NVa,i),Fa(jv+m*nbas),NVa)
      call vecdot(dot2,R(ov4+m*NVa,i),Fb(jv+m*nbas),NVa)
      Az(2,i) = Az(2,i) + ssq2*(dot1+dot2)
    enddo
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    endif ! I1
  
    ! s1->s1
    Az(3,i) = -sq2*Js(s1s2+1,1)*R(1,i)-sq2*Js(s1s2+1,1)*R(2,i)                     &
              -(Js(s1s1+1,1)-Js(s1s2+1,2))*R(3,i)                                  &
              +Ka(jJos1+s11,i)-Ka(jJos1+s22,i)+Ka(jJos2+s11,i)-Ka(jJos2+s22,i)     &
              +Ka(jJvs1+s11,i)-Ka(jJvs1+s22,i)+Ka(jJvs2+s11,i)-Ka(jJvs2+s22,i)     &
              +ssq6*Ka(jJov1+s22,i)+sq2*Ka(jJov2+s11,i)-ssq2*Ka(jJov2+s22,i)
    ! Fock
    Az(3,i) = Az(3,i) + sq2*Fb(s2s1+1)*R(1,i) - sq2*Fa(s2s1+1)*R(2,i)
    Az(3,i) = Az(3,i) + (Fb(s1s1+1)-Fa(s1s1+1))*R(3,i)
    call vecdot(dot1,R(os1,i),Fa(i1s1+1),NOb)
    call vecdot(dot2,R(os2,i),Fa(i1s2+1),NOb)
    call vecdot(dot3,R(vs1,i),Fb(a1s1+1),NVa)
    call vecdot(dot4,R(vs2,i),Fb(a1s2+1),NVa)
    Az(3,i) =  Az(3,i) - dot1 + dot2 + dot3 - dot4
    do m=0,NOb-1
      call vecdot(dot1,R(ov1+m*NVa,i),Fb(m*nbas+NOa+1),NVa)
      call vecdot(dot2,R(ov2+m*NVa,i),Fb(m*nbas+NOa+1),NVa)
      Az(3,i) = Az(3,i) + ssq6*dot1 + ssq2*dot2
    enddo
    if(sf_xcis) then
    Az(3,i) = Az(3,i) + Ka(jJov3+s12,i) - Ka(jJov4+s21,i)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! occ->s1
    ioff1 = 0
    ioff2 = 1
    do m=0,NOb-1
      Az(os1+m,i) =  sq2*Ks(i1s1+m+1,2)*R(1,i)+sq2*Ks(NOb+1+m*nbas,2)*R(2,i)      &
                    +(Ks(NOb+1+m*nbas,1)-Ks(NOb+1+m*nbas,3))*R(3,i)               &
                    +Ka(jJos1+ioff1,i)+Ja(jJos12+ioff2,i)                         &
                    +Ka(jJos2+ioff1,i)-Ja(jJos2+ioff1,i)                          &
                    +Ka(jJvs1+ioff1,i)                                            &
                    +Ja(jJvs2+ioff1,i)+two*Ka(jJvs2+ioff1,i)                      &
                    -ssq6*Ja(jJov1+ioff1,i)-ssq6*Ka(jJov1+ioff1,i)                &
                    -ssq2*Ja(jJov2+ioff1,i)+ssq2*Ka(jJov2+ioff1,i)
      ioff1 = ioff1 + NVb
      ioff2 = ioff2 + NVb
    enddo
    call matmult(R(ov1,i),Ks(a1s1+1,3),tmp,NOb,1,NVa,NVa,nbas,NOb,2)
    call vecadd3(Az(os1,i),-ssq6,tmp,NOb)
    call matmult(R(ov2,i),Ks(a1s1+1,3),tmp,NOb,1,NVa,NVa,nbas,NOb,2)
    call vecadd3(Az(os1,i),-ssq2,tmp,NOb)
    ! Fock contribution
    call vecadd3(Az(os1,i),-sq2*R(2,i),Fa(i1s2+1),NOb)
    call vecadd3(Az(os1,i),-R(3,i),Fa(i1s1+1),NOb)
    call vecadd3(Az(os1,i),Fb(s1s1+1),R(os1,i),NOb)
    call matmult(Fa,R(os1,i),tmp,NOb,1,NOb,nbas,NOb,NOb,1)
    call vecsub(Az(os1,i),Az(os1,i),tmp,NOb)
    call vecadd3(Az(os1,i),Fb(s1s2+1),R(os2,i),NOb)
    call matmult(R(ov1,i),Fb(a1s1+1),tmp,NOb,1,NVa,NVa,nbas,NOb,2)
    call vecadd3(Az(os1,i),-ssq6,tmp,NOb)
    call matmult(R(ov2,i),Fb(a1s1+1),tmp,NOb,1,NVa,NVa,nbas,NOb,2)
    call vecadd3(Az(os1,i),ssq2,tmp,NOb)
    if(sf_xcis) then
    ioff2 = 1
    do m=0,NOb-1
      Az(os1+m,i) = Az(os1+m,i) - Ka(jJov4+ioff2,i)
      ioff2 = ioff2 + NVb
    enddo
    call matmult(R(ov3,i),Js(a1s1+1,2),tmp,NOb,1,NVa,NVa,nbas,NOb,2)
    call vecsub(Az(os1,i),Az(os1,i),tmp,NOb)
    call matmult(R(ov4,i),Js(a1s2+1,1),tmp,NOb,1,NVa,NVa,nbas,NOb,2)
    call vecsub(Az(os1,i),Az(os1,i),tmp,NOb)
    call matmult(R(ov4,i),Fa(a1s2+1),tmp,NOb,1,NVa,NVa,nbas,NOb,2)
    call vecsub(Az(os1,i),Az(os1,i),tmp,NOb)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! occ->s2
    ioff1 = 0
    ioff2 = 1
    do m=0,NOb-1
      Az(os2+m,i) =  sq2*Ks(i1s2+m+1,2)*R(1,i)+sq2*Ks(js2+m*nbas,2)*R(2,i)        &
                    +(Ks(js2+m*nbas,1)-Ks(js2+m*nbas,3))*R(3,i)                   &
                    -Ja(jJos1+ioff2,i)+Ka(jJos1+ioff2,i)                          &
                    +Ja(jJos21+ioff1,i)+Ka(jJos2+ioff2,i)                         &
                    +two*Ka(jJvs1+ioff2,i)+Ja(jJvs1+ioff2,i)                      &
                    +Ka(jJvs2+ioff2,i)                                            &
                    +ssq6*Ja(jJov1+ioff2,i)                                       &
                    +ssq2*Ja(jJov2+ioff2,i)+sq2*Ka(jJov2+ioff2,i)
      ioff1 = ioff1 + NVb
      ioff2 = ioff2 + NVb
    enddo
    call matmult(R(ov1,i),Ks(a1s2+1,1),tmp,NOb,1,NVa,NVa,nbas,NOb,2)
    call vecadd3(Az(os2,i),ssq6,tmp,NOb)
    call matmult(R(ov2,i),Ks(a1s2+1,1),tmp,NOb,1,NVa,NVa,nbas,NOb,2)
    call vecadd3(Az(os2,i),ssq2,tmp,NOb)
    ! Fock contribution
    call vecadd3(Az(os2,i),-sq2*R(1,i),Fa(i1s1+1),NOb)
    call vecadd3(Az(os2,i),R(3,i),Fa(i1s2+1),NOb)
    call vecadd3(Az(os2,i),Fb(s1s2+1),R(os1,i),NOb)
    call vecadd3(Az(os2,i),Fb(s2s2+1),R(os2,i),NOb) !(s2i,s2j)
    call matmult(Fa,R(os2,i),tmp,NOb,1,NOb,nbas,NOb,NOb,1)
    call vecsub(Az(os2,i),Az(os2,i),tmp,NOb)
    call matmult(R(ov2,i),Fb(a1s2+1),tmp,NOb,1,NVa,NVa,nbas,NOb,2) !(s2i,aj2)
    call vecadd3(Az(os2,i),sq2,tmp,NOb)
    if(sf_xcis) then
    ioff1 = 0
    do m=0,NOb-1
      Az(os2+m,i) = Az(os2+m,i) - Ka(jJov3+ioff1,i)
      ioff1 = ioff1 + NVb
    enddo
    call matmult(R(ov3,i),Js(a1s1+1,3),tmp,NOb,1,NVa,NVa,nbas,NOb,2)
    call vecsub(Az(os2,i),Az(os2,i),tmp,NOb)
    call matmult(R(ov3,i),Fa(a1s1+1),tmp,NOb,1,NVa,NVa,nbas,NOb,2)
    call vecsub(Az(os2,i),Az(os2,i),tmp,NOb)
    call matmult(R(ov4,i),Js(a1s2+1,2),tmp,NOb,1,NVa,NVa,nbas,NOb,2)
    call vecsub(Az(os2,i),Az(os2,i),tmp,NOb)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! s1->vir
    do m=0,NVa-1
      Az(vs1+m,i) =  sq2*Ks(js1+(NOa+m)*nbas,2)*R(1,i)+sq2*Ks(a1s1+m+1,2)*R(2,i)  &
                    +(Ks(a1s1+m+1,1)-Ks(a1s1+m+1,3))*R(3,i)                       &
                    +Ka(jJos1+as1+m,i)+two*Ka(jJos2+as1+m,i)+Ja(jJos2+as1+m,i)    &
                    +Ka(jJvs1+as1+m,i)+Ja(jJvs12+as2+m,i)                         &
                    +Ka(jJvs2+as1+m,i)-Ja(jJvs2+as1+m,i)                          &
                    +ssq6*Ja(jJov1+as1+m,i)+ssq2*Ja(jJov2+as1+m,i)                &
                    +sq2*Ka(jJov2+as1+m,i)
    enddo
    call matmult(R(ov1,i),Ks(i1s1+1,3),tmp,NVa,1,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(vs1,i),ssq6,tmp,NVa)
    call matmult(R(ov2,i),Ks(i1s1+1,3),tmp,NVa,1,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(vs1,i),ssq2,tmp,NVa)
    ! Fock contribution
    call vecadd3(Az(vs1,i),sq2*R(1,i),Fb(a1s2+1),NVa)
    call vecadd3(Az(vs1,i),R(3,i),Fb(a1s1+1),NVa)
    call matmult(Fb(a1b1+1),R(vs1,i),tmp,NVa,1,NVa,nbas,NVa,NVa,1)
    call vecadd(Az(vs1,i),Az(vs1,i),tmp,NVa)
    call vecadd3(Az(vs1,i),-Fa(s1s1+1),R(vs1,i),NVa)
    call vecadd3(Az(vs1,i),-Fa(s1s2+1),R(vs2,i),NVa)
    call matmult(R(ov2,i),Fa(i1s1+1),tmp,NVa,1,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(vs1,i),-sq2,tmp,NVa)
    if(sf_xcis) then
    do m=0,NVa-1
      Az(vs1+m,i) = Az(vs1+m,i) + Ka(jJov3+as2+m,i)
    enddo
    call matmult(R(ov3,i),Js(i1s2+1,1),tmp,NVa,1,NOb,NVa,nbas,NVa,1)
    call vecadd(Az(vs1,i),Az(vs1,i),tmp,NVa)
    call matmult(R(ov3,i),Fb(i1s2+1),tmp,NVa,1,NOb,NVa,nbas,NVa,1)
    call vecsub(Az(vs1,i),Az(vs1,i),tmp,NVa)
    call matmult(R(ov4,i),Js(i1s1+1,2),tmp,NVa,1,NOb,NVa,nbas,NVa,1)
    call vecadd(Az(vs1,i),Az(vs1,i),tmp,NVa)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! s2->vir
    do m=0,NVa-1
      Az(vs2+m,i) =  sq2*Ks(js2+(NOa+m)*nbas,2)*R(1,i)+sq2*Ks(a1s2+m+1,2)*R(2,i)  &
                    +(Ks(a1s2+m+1,1)-Ks(a1s2+m+1,3))*R(3,i)                       &
                    +two*Ka(jJos1+as2+m,i)+Ja(jJos1+as2+m,i)+Ka(jJos2+as2+m,i)    &
                    +Ka(jJvs1+as2+m,i)-Ja(jJvs1+as2+m,i)                          &
                    +Ka(jJvs2+as2+m,i)+Ja(jJvs21+as1+m,i)                         &
                    -ssq6*Ja(jJov1+as2+m,i)-ssq6*Ka(jJov1+as2+m,i)                &
                    +ssq2*Ka(jJov2+as2+m,i)-ssq2*Ja(jJov2+as2+m,i)
    enddo
    call matmult(R(ov1,i),Ks(i1s2+1,1),tmp,NVa,1,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(vs2,i),-ssq6,tmp,NVa)
    call matmult(R(ov2,i),Ks(i1s2+1,1),tmp,NVa,1,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(vs2,i),-ssq2,tmp,NVa)
    ! Fock contribution
    call vecadd3(Az(vs2,i),sq2*R(2,i),Fb(a1s1+1),NVa)
    call vecadd3(Az(vs2,i),-R(3,i),Fb(a1s2+1),NVa)
    call vecadd3(Az(vs2,i),-Fa(s1s2+1),R(vs1,i),NVa)
    call matmult(Fb(a1b1+1),R(vs2,i),tmp,NVa,1,NVa,nbas,NVa,NVa,1)
    call vecadd(Az(vs2,i),Az(vs2,i),tmp,NVa)
    call vecadd3(Az(vs2,i),-Fa(s2s2+1),R(vs2,i),NVa)
    call matmult(R(ov1,i),Fa(i1s2+1),tmp,NVa,1,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(vs2,i),ssq6,tmp,NVa)
    call matmult(R(ov2,i),Fa(i1s2+1),tmp,NVa,1,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(vs2,i),-ssq2,tmp,NVa)
    if(sf_xcis) then
    do m=0,NVa-1
      Az(vs2+m,i) = Az(vs2+m,i) + Ka(jJov4+as1+m,i)
    enddo
    call matmult(R(ov3,i),Js(i1s2+1,2),tmp,NVa,1,NOb,NVa,nbas,NVa,1)
    call vecadd(Az(vs2,i),Az(vs2,i),tmp,NVa)
    call matmult(R(ov4,i),Js(i1s1+1,3),tmp,NVa,1,NOb,NVa,nbas,NVa,1)
    call vecadd(Az(vs2,i),Az(vs2,i),tmp,NVa)
    call matmult(R(ov4,i),Fb(i1s1+1),tmp,NVa,1,NOb,NVa,nbas,NVa,1)
    call vecsub(Az(vs2,i),Az(vs2,i),tmp,NVa)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! occ->vir A
    ioff1 = ov1
    do m=0,NOb-1
      call vecadd3(Az(ioff1,i),sq3*R(1,i),Js(jv+m*nbas,2),NVa)
      call vecadd(tmp,Js(jv+m*nbas,2),Ks(jv+m*nbas,2),NVa)
      call vecadd3(Az(ioff1,i),-sq3*R(2,i),tmp,NVa)
      call vecadd3(Az(ioff1,i),ssq6*R(3,i),Ks(jv+m*nbas,3),NVa)
      call vecadd(tmp,Ja(jJos1+2+m*NVb,i),Ka(jJos1+2+m*NVb,i),NVa)
      call vecadd3(Az(ioff1,i),-ssq6,tmp,NVa)
      call vecadd3(Az(ioff1,i),ssq6,Ja(jJos2+2+m*NVb,i),NVa)
      call vecadd3(Az(ioff1,i),ssq6,Ja(jJvs1+2+m*NVb,i),NVa)
      call vecadd(tmp,Ja(jJvs2+2+m*NVb,i),Ka(jJvs2+2+m*NVb,i),NVa)
      call vecadd3(Az(ioff1,i),-ssq6,tmp,NVa)
      call vecadd2(tmp,f32,Ja(jJov1+2+m*NVb,i),Ka(jJov1+2+m*NVb,i),NVa)
      call vecadd(Az(ioff1,i),Az(ioff1,i),tmp,NVa)
      call vecadd3(Az(ioff1,i),ssq3,Ja(jJov2+2+m*NVb,i),NVa)
      ! Fock contribution
      call vecadd3(Az(ioff1,i),ssq6*R(3,i),Fb(NOa+1+m*nbas),NVa)
      ioff1 = ioff1 + NVa
    enddo
    call matmult(Ks(a1s1+1,3),R(os1,i),tmp,NVa,NOb,1,nbas,NOb,NVa,3)
    call vecadd3(Az(ov1,i),-ssq6,tmp,NVa*NOb)
    call matmult(Ks(a1s2+1,1),R(os2,i),tmp,NVa,NOb,1,nbas,NOb,NVa,3)
    call vecadd3(Az(ov1,i),ssq6,tmp,NVa*NOb)
    call matmult(R(vs1,i),Ks(i1s1+1,3),tmp,NVa,NOb,1,NVa,nbas,NVa,3)
    call vecadd3(Az(ov1,i),ssq6,tmp,NVa*NOb)
    call matmult(R(vs2,i),Ks(i1s2+1,1),tmp,NVa,NOb,1,NVa,nbas,NVa,3)
    call vecadd3(Az(ov1,i),-ssq6,tmp,NVa*NOb)
    call vecadd2(tmp2,-f12,Ks(1,1),Ks(1,3),nbas*nbas)
    call matmult(tmp2(a1b1+1),R(ov1,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
    call vecadd(Az(ov1,i),Az(ov1,i),tmp,NVa*NOb)
    call vecadd2(tmp2,-f12,Ks(1,3),Ks(1,1),nbas*nbas)
    call matmult(R(ov1,i),tmp2,tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
    call vecadd(Az(ov1,i),Az(ov1,i),tmp,NVa*NOb)
    call matmult(Ks(a1b1+1,1),R(ov2,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
    call vecadd3(Az(ov1,i),ssq3,tmp,NVa*NOb)
    call matmult(R(ov2,i),Ks(1,3),tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(ov1,i),ssq3,tmp,NVa*NOb)
    ! diagnol
    call vecadd3(Az(ov1,i),f32*Js(s1s2+1,2),R(ov1,i),NVa*NOb)
    call vecadd3(Az(ov1,i),ssq3*Js(s1s2+1,2),R(ov2,i),NVa*NOb)
    ! Fock
    call matmult(Fb(a1s1+1),R(os1,i),tmp,NVa,NOb,1,nbas,NOb,NVa,3)
    call vecadd3(Az(ov1,i),-ssq6,tmp,NVa*NOb)
    call matmult(R(vs2,i),Fa(i1s2+1),tmp,NVa,NOb,1,NVa,nbas,NVa,3)
    call vecadd3(Az(ov1,i),ssq6,tmp,NVa*NOb)
    call matmult(Fb(a1b1+1),R(ov1,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
    call vecadd(Az(ov1,i),Az(ov1,i),tmp,NVa*NOb)
    call matmult(R(ov1,i),Fa,tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
    call vecsub(Az(ov1,i),Az(ov1,i),tmp,NVa*NOb)
    if(sf_xcis) then
    dot1 = ssq6*(Fa(s1s2+1)+Js(s1s2+1,3))
    call vecadd3(Az(ov1,i),dot1,R(ov3,i),NVa*NOb)
    call matmult(R(ov3,i),Js(1,2),tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(ov1,i),-ssq6,tmp,NVa*NOb)
    call matmult(Js(a1b1+1,2),R(ov3,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
    call vecadd3(Az(ov1,i),ssq6,tmp,NVa*NOb)

    dot1 = -ssq6*(Fa(s1s2+1)+Js(s1s2+1,1))
    call vecadd3(Az(ov1,i),dot1,R(ov4,i),NVa*NOb)
    call vecadd(tmp2,Js(1,2),Ks(1,2),nbas*nbas)
    call matmult(R(ov4,i),tmp2,tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(ov1,i),ssq6,tmp,NVa*NOb)
    call matmult(tmp2(a1b1+1),R(ov4,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
    call vecadd3(Az(ov1,i),-ssq6,tmp,NVa*NOb)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! occ->vir B
    ioff2 = ov2
    do m=0,NOb-1
      call mattrans(tmp2,Ks(1,2),nbas,nbas)
      call vecadd2(tmp,two,tmp2(jv+m*nbas),Js(jv+m*nbas,2),NVa)
      call vecadd3(Az(ioff2,i),R(1,i),tmp,NVa)
      call vecsub(tmp,Js(jv+m*nbas,2),Ks(jv+m*nbas,2),NVa)
      call vecadd3(Az(ioff2,i),-R(2,i),tmp,NVa)
      call vecadd2(tmp,-two,Ks(jv+m*nbas,1),Ks(jv+m*nbas,3),NVa)
      call vecadd3(Az(ioff2,i),-ssq2*R(3,i),tmp,NVa)
      call vecsub(tmp,Ja(jJos1+2+m*NVb,i),Ka(jJos1+2+m*NVb,i),NVa)
      call vecadd3(Az(ioff2,i),-ssq2,tmp,NVa)
      call vecadd2(tmp,two,Ka(jJos2+2+m*NVb,i),Ja(jJos2+2+m*NVb,i),NVa)
      call vecadd3(Az(ioff2,i),ssq2,tmp,NVa)
      call vecadd2(tmp,two,Ka(jJvs1+2+m*NVb,i),Ja(jJvs1+2+m*NVb,i),NVa)
      call vecadd3(Az(ioff2,i),ssq2,tmp,NVa)
      call vecsub(tmp,Ja(jJvs2+2+m*NVb,i),Ka(jJvs2+2+m*NVb,i),NVa)
      call vecadd3(Az(ioff2,i),-ssq2,tmp,NVa)
      call vecadd3(Az(ioff2,i),ssq3,Ja(jJov1+2+m*NVb,i),NVa)
      call vecadd2(tmp,f12,Ja(jJov2+2+m*NVb,i),Ka(jJov2+2+m*NVb,i),NVa)
      call vecadd(Az(ioff2,i),Az(ioff2,i),tmp,NVa)

      call vecadd3(Az(ioff2,i),ssq2*R(3,i),Fb(NOa+1+m*nbas),NVa)
      ioff2 = ioff2 + NVa
    enddo
    call matmult(Ks(a1s1+1,3),R(os1,i),tmp,NVa,NOb,1,nbas,NOb,NVa,3)
    call vecadd3(Az(ov2,i),-ssq2,tmp,NVa*NOb)
    call matmult(Ks(a1s2+1,1),R(os2,i),tmp,NVa,NOb,1,nbas,NOb,NVa,3)
    call vecadd3(Az(ov2,i),ssq2,tmp,NVa*NOb)
    call matmult(R(vs1,i),Ks(i1s1+1,3),tmp,NVa,NOb,1,NVa,nbas,NVa,3)
    call vecadd3(Az(ov2,i),ssq2,tmp,NVa*NOb)
    call matmult(R(vs2,i),Ks(i1s2+1,1),tmp,NVa,NOb,1,NVa,nbas,NVa,3)
    call vecadd3(Az(ov2,i),-ssq2,tmp,NVa*NOb)
    call vecadd2(tmp2,-f12,Ks(1,1),Ks(1,3),nbas*nbas)
    call matmult(tmp2(a1b1+1),R(ov2,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
    call vecsub(Az(ov2,i),Az(ov2,i),tmp,NVa*NOb)
    call vecadd2(tmp2,-f12,Ks(1,3),Ks(1,1),nbas*nbas)
    call matmult(R(ov2,i),tmp2,tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
    call vecsub(Az(ov2,i),Az(ov2,i),tmp,NVa*NOb)
    call matmult(Ks(a1b1+1,1),R(ov1,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
    call vecadd3(Az(ov2,i),ssq3,tmp,NVa*NOb)
    call matmult(R(ov1,i),Ks(1,3),tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(ov2,i),ssq3,tmp,NVa*NOb)
    ! diagnol
    call vecadd3(Az(ov2,i),f12*Js(s1s2+1,2),R(ov2,i),NVa*NOb)
    call vecadd3(Az(ov2,i),ssq3*Js(s1s2+1,2),R(ov1,i),NVa*NOb)
    ! Fock
    call matmult(Fb(a1s1+1),R(os1,i),tmp,NVa,NOb,1,nbas,NOb,NVa,3)
    call vecadd3(Az(ov2,i),ssq2,tmp,NVa*NOb)
    call matmult(Fb(a1s2+1),R(os2,i),tmp,NVa,NOb,1,nbas,NOb,NVa,3)
    call vecadd3(Az(ov2,i),sq2,tmp,NVa*NOb)
    call matmult(R(vs1,i),Fa(i1s1+1),tmp,NVa,NOb,1,NVa,nbas,NVa,3)
    call vecadd3(Az(ov2,i),-sq2,tmp,NVa*NOb)
    call matmult(R(vs2,i),Fa(i1s2+1),tmp,NVa,NOb,1,NVa,nbas,NVa,3)
    call vecadd3(Az(ov2,i),-ssq2,tmp,NVa*NOb)
    call matmult(Fb(a1b1+1),R(ov2,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
    call vecadd(Az(ov2,i),Az(ov2,i),tmp,NVa*NOb)
    call matmult(R(ov2,i),Fa,tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
    call vecsub(Az(ov2,i),Az(ov2,i),tmp,NVa*NOb)
    if(sf_xcis) then
      dot1 = ssq2*(Fa(s1s2+1)+Js(s1s2+1,3))
      call vecadd3(Az(ov2,i),dot1,R(ov3,i),NVa*NOb)
      call mattrans(tmp2,Ks(1,2),nbas,nbas)
      call vecadd2(tmp2,two,tmp2,Js(1,2),nbas*nbas)
      call matmult(R(ov3,i),tmp2,tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
      call vecadd3(Az(ov2,i),-ssq2,tmp,NVa*NOb)
      call matmult(tmp2(a1b1+1),R(ov3,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
      call vecadd3(Az(ov2,i),ssq2,tmp,NVa*NOb)
     
      dot1 = -ssq2*(Fa(s1s2+1)+Js(s1s2+1,1))
      call vecadd3(Az(ov2,i),dot1,R(ov4,i),NVa*NOb)
      call vecsub(tmp2,Js(1,2),Ks(1,2),nbas*nbas)
      call matmult(R(ov4,i),tmp2,tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
      call vecadd3(Az(ov2,i),ssq2,tmp,NVa*NOb)
      call matmult(tmp2(a1b1+1),R(ov4,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
      call vecadd3(Az(ov2,i),-ssq2,tmp,NVa*NOb)
    endif

    if(sf_xcis) then
    ! occ -> vir C
    ioff2 = ov3
    do m=0,NOb-1
      call vecadd2(tmp,two,Js(jv+m*nbas,3),Ks(jv+m*nbas,3),NVa)
      call vecadd3(Az(ioff2,i),ssq2*R(1,i),tmp,NVa)
      call vecadd2(tmp,two,Js(jv+m*nbas,1),Ks(jv+m*nbas,1),NVa)
      call vecadd3(Az(ioff2,i),-ssq2*R(1,i),tmp,NVa)
      call vecadd2(tmp,two,Js(jv+m*nbas,2),Ks(jv+m*nbas,2),NVa)
      call vecadd3(Az(ioff2,i),R(3,i),tmp,NVa)
      call vecadd3(Az(ioff2,i),-two,Ja(jJos21+2+m*NVb,i),NVa)
      call vecsub(Az(ioff2,i),Az(ioff2,i),Ka(jJos21+2+m*NVb,i),NVa)
      call vecadd3(Az(ioff2,i),two,Ja(jJvs12+2+m*NVb,i),NVa)
      call vecadd(Az(ioff2,i),Az(ioff2,i),Ka(jJvs12+2+m*NVb,i),NVa)
      call vecadd(Az(ioff2,i),Az(ioff2,i),Ka(jJov3+2+m*NVb,i),NVa)

      call vecadd3(Az(ioff2,i),ssq2*R(1,i),Fa(jv+m*nbas),NVa)
      call vecadd3(Az(ioff2,i),ssq2*R(1,i),Fb(jv+m*nbas),NVa)
      call vecadd3(Az(ioff2,i),-R(os1+m,i),Js(a1s1+1,2),NVa)
      call vecadd(tmp,Fa(a1s1+1),Js(a1s1+1,3),NVa)
      call vecadd3(Az(ioff2,i),-R(os2+m,i),tmp,NVa)
      dot1 = -Fb(i1s2+1+m) + Js(i1s2+1+m,1)
      call vecadd3(Az(ioff2,i),dot1,R(vs1,i),NVa)
      call vecadd3(Az(ioff2,i),Js(i1s2+1+m,2),R(vs2,i),NVa)
      ioff2 = ioff2 + NVa
    enddo
    dot1 = ssq6*(Fa(s1s2+1)+Js(s1s2+1,3))
    call vecadd3(Az(ov3,i),dot1,R(ov1,i),NVa*NOb)
    dot1 = ssq2*(Fa(s1s2+1)+Js(s1s2+1,3))
    call vecadd3(Az(ov3,i),dot1,R(ov2,i),NVa*NOb)
    dot1 = Fb(s2s2+1)-Fa(s1s1+1)-Js(s1s1+1,3)
    call vecadd3(Az(ov3,i),dot1,R(ov3,i),NVa*NOb)
    dot1 = -Js(s1s2+1,2)
    call vecadd3(Az(ov3,i),dot1,R(ov4,i),NVa*NOb)

    call matmult(Js(a1b1+1,2),R(ov1,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
    call vecadd3(Az(ov3,i),ssq6,tmp,NVa*NOb)
    call matmult(R(ov1,i),Js(1,2),tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(ov3,i),-ssq6,tmp,NVa*NOb)

    call vecadd2(tmp2,two,Ks(1,2),Js(1,2),nbas*nbas)
    call matmult(tmp2(a1b1+1),R(ov2,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
    call vecadd3(Az(ov3,i),ssq2,tmp,NVa*NOb)
    call matmult(R(ov2,i),tmp2,tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(ov3,i),-ssq2,tmp,NVa*NOb)

    call vecadd(tmp2,Fa,Js(1,3),nbas*nbas)
    call vecsub(tmp2,tmp2,Js(1,1),nbas*nbas)
    call vecsub(tmp2,tmp2,Ks(1,1),nbas*nbas)
    call matmult(tmp2(a1b1+1),R(ov3,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
    call vecadd(Az(ov3,i),Az(ov3,i),tmp,NVa*NOb)
    call matmult(R(ov3,i),tmp2,tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
    call vecsub(Az(ov3,i),Az(ov3,i),tmp,NVa*NOb)

    ! occ -> vir D
    ioff2 = ov4
    do m=0,NOb-1
      call vecadd2(tmp,two,Js(jv+m*nbas,3),Ks(jv+m*nbas,3),NVa)
      call vecadd3(Az(ioff2,i),-ssq2*R(2,i),tmp,NVa)
      call vecadd2(tmp,two,Js(jv+m*nbas,1),Ks(jv+m*nbas,1),NVa)
      call vecadd3(Az(ioff2,i),ssq2*R(2,i),tmp,NVa)
      call mattrans(tmp2,Ks(1,2),nbas,nbas)
      call vecadd2(tmp,two,Js(jv+m*nbas,2),tmp2(jv+m*nbas),NVa)
      call vecadd3(Az(ioff2,i),-R(3,i),tmp,NVa)
      call vecadd3(Az(ioff2,i),-two,Ja(jJos12+2+m*NVb,i),NVa)
      call vecsub(Az(ioff2,i),Az(ioff2,i),Ka(jJos12+2+m*NVb,i),NVa)
      call vecadd3(Az(ioff2,i),two,Ja(jJvs21+2+m*NVb,i),NVa)
      call vecadd(Az(ioff2,i),Az(ioff2,i),Ka(jJvs21+2+m*NVb,i),NVa)
      call vecadd(Az(ioff2,i),Az(ioff2,i),Ka(jJov4+2+m*NVb,i),NVa)

      call vecadd3(Az(ioff2,i),ssq2*R(2,i),Fa(jv+m*nbas),NVa)
      call vecadd3(Az(ioff2,i),ssq2*R(2,i),Fb(jv+m*nbas),NVa)
      call vecadd(tmp,Fa(a1s2+1),Js(a1s2+1,1),NVa)
      call vecadd3(Az(ioff2,i),-R(os1+m,i),tmp,NVa)
      call vecadd3(Az(ioff2,i),-R(os2+m,i),Js(a1s2+1,2),NVa)
      call vecadd3(Az(ioff2,i),Js(i1s1+1+m,2),R(vs1,i),NVa)
      dot1 = -Fb(i1s1+1+m) + Js(i1s1+1+m,3)
      call vecadd3(Az(ioff2,i),dot1,R(vs2,i),NVa)
      ioff2 = ioff2 + NVa
    enddo
    dot1 = -ssq6*(Fa(s1s2+1)+Js(s1s2+1,1))
    call vecadd3(Az(ov4,i),dot1,R(ov1,i),NVa*NOb)
    dot1 = -ssq2*(Fa(s1s2+1)+Js(s1s2+1,1))
    call vecadd3(Az(ov4,i),dot1,R(ov2,i),NVa*NOb)
    dot1 = -Js(s1s2+1,2)
    call vecadd3(Az(ov4,i),dot1,R(ov3,i),NVa*NOb)
    dot1 = Fb(s1s1+1)-Fa(s2s2+1)-Js(s1s1+1,3)
    call vecadd3(Az(ov4,i),dot1,R(ov4,i),NVa*NOb)

    call mattrans(tmp2,Ks(1,2),nbas,nbas)
    call vecadd(tmp2,Js(1,2),tmp2,nbas*nbas)
    call matmult(tmp2(a1b1+1),R(ov1,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
    call vecadd3(Az(ov4,i),-ssq6,tmp,NVa*NOb)
    call matmult(R(ov1,i),tmp2,tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(ov4,i),ssq6,tmp,NVa*NOb)

    call mattrans(tmp2,Ks(1,2),nbas,nbas)
    call vecsub(tmp2,Js(1,2),tmp2,nbas*nbas)
    call matmult(tmp2(a1b1+1),R(ov2,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
    call vecadd3(Az(ov4,i),-ssq2,tmp,NVa*NOb)
    call matmult(R(ov2,i),tmp2,tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
    call vecadd3(Az(ov4,i),ssq2,tmp,NVa*NOb)

    call vecadd(tmp2,Fa,Js(1,1),nbas*nbas)
    call vecsub(tmp2,tmp2,Js(1,3),nbas*nbas)
    call vecsub(tmp2,tmp2,Ks(1,3),nbas*nbas)
    call matmult(tmp2(a1b1+1),R(ov4,i),tmp,NVa,NOb,NVa,nbas,NVa,NVa,1)
    call vecadd(Az(ov4,i),Az(ov4,i),tmp,NVa*NOb)
    call matmult(R(ov4,i),tmp2,tmp,NVa,NOb,NOb,NVa,nbas,NVa,1)
    call vecsub(Az(ov4,i),Az(ov4,i),tmp,NVa*NOb)

    endif

  enddo

  deallocate(tmp)
  deallocate(tmp2)
  return
end subroutine sasfemake_amo

!======================================================== 
subroutine sasfcisprt(V,iroots)
  use limit, only : lmprop, lmstat
  use cissym
  implicit none
  integer i, j, k, m, iroots, mult, ioff, isym 
  real*8  V(*), dmax, escf, enuc, au2dby, ciprop
  real*8, dimension(:,:), allocatable :: tdip,P,R 
  common  /energt/ escf,enuc
  common  /ciprp / ciprop(lmprop,lmstat)

  if(iprtci.le.0) return

  ! determine the symmetry labels for all active orbitals
  call ovsym(V)

  ! calculate the oscillator strengths and transition moments
  allocate(tdip(4,iroots))
  call vecinit(tdip,4*iroots,0.d0)
  if(iuvcd.gt.0) then
    allocate(P(n2,nroots),R(n2,nroots))
    !call vecinit(P,n2*nroots,0.d0)
    call get_sasfcis_density(P,Xv,V(jCa),V(jCb),V(jPa),V(jPb))
    call get_sasfcis_transition_density(R,Xv,V(jCa),V(jCb))
    call cisprop(P,R,Es(1,1),ciprop,jsym,labels,jsym,iroots,LM3)
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
       if(m.le.(iroots/2)) mult = 3
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
            m, mult, labels(jsym(m)), jsym(m), es(m,1)-es(1,1), es(m,1)+escf+enuc+et3
      write(nb6,'(A/)') '-------------------------------------------------------------------------'
    else 
      write(nb6,'(A/)') '-------------------------------------------------------------------------'
      write(nb6,'(" State",I3,",  Mult.",I2,", E-E(1)=",F10.6," eV,  E=",F14.6," eV")') &
            m, mult, es(m,1)-es(1,1), es(m,1)+escf+enuc+et3
      write(nb6,'(A/)') '-------------------------------------------------------------------------'
    endif
 
    write(*,'(" Px: ",F6.2,"     Py: ",F6.2,"     Pz: ",F6.2)') &
      tdip(2,m)*au2dby, tdip(3,m)*au2dby, tdip(4,m)*au2dby
    write(*,'(" Strength:      ",F6.2)') tdip(1,m)

    ! print out the primary single excitations
    if(abs(Xv(1,m)).ge.1.d-1) then
       write(*,'("S  1 --> S  2",F15.3)') Xv(1,m)
    endif
    if(abs(Xv(2,m)).ge.1.d-1) then
       write(*,'("S  2 --> S  1",F15.3)') Xv(2,m)
    endif
    if(abs(Xv(3,m)).ge.1.d-1) then
       write(*,'("S  1 --> S  1",F15.3)') Xv(3,m)
    endif
    ioff = 4
    do i=1,2
      do j=1,NOb
        if(abs(Xv(ioff,m)).ge.1.d-1) then
          write(*,'("D",I3," --> S",I3,F15.3)') j,i,Xv(ioff,m)
        endif
        ioff = ioff + 1
      enddo
    enddo
    do i=1,2
      do j=1,NVa
        if(abs(Xv(ioff,m)).ge.1.d-1) then
          write(*,'("S",I3," --> V",I3,F15.3)') i,j,Xv(ioff,m)
        endif
        ioff = ioff + 1
      enddo
    enddo
    do i=1,NOb
      do j=1,NVa
        if(abs(Xv(ioff,m)).ge.1.d-1) then
          write(*,'("D",I3," -->aV",I3,F15.3)') i,j,Xv(ioff,m)
        endif
        if(abs(Xv(ioff+NVa*NOb,m)).ge.1.d-1) then
          write(*,'("D",I3," -->bV",I3,F15.3)') i,j,Xv(ioff+NVa*NOb,m)
        endif
        if(sf_xcis) then
        if(abs(Xv(ioff+2*NVa*NOb,m)).ge.1.d-1) then
          write(*,'("D",I3," -->cV",I3,F15.3)') i,j,Xv(ioff+2*NVa*NOb,m)
        endif
        if(abs(Xv(ioff+3*NVa*NOb,m)).ge.1.d-1) then
          write(*,'("D",I3," -->dV",I3,F15.3)') i,j,Xv(ioff+3*NVa*NOb,m)
        endif
        endif
        ioff = ioff + 1
      enddo
    enddo
  enddo

  deallocate(tdip)
  return
end subroutine sasfcisprt



end module sasfcis_energy
