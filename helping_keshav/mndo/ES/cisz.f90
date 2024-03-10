!
!  created by jie liu on 10/16
!
module cis_z

use es_global
use es_gradient
use cis_solve
use cisz_solve

public cis_nac_z, get_cis_density_mn

contains

!======================================================== 
subroutine cis_nac_z(V)
  implicit none
  real*8 V(*)

  call z_nac_init 

  call z_nac_rhs(V)

  call z_iter(V)

  call z_nac_exit

end subroutine cis_nac_z

!======================================================== 
subroutine z_nac_init

  if(iprtz.gt.0) then
    write(NB6,*)
    write(NB6,*) "**************************************"
    write(NB6,*) "*           ITER_NAC_Z               *"
    write(NB6,*) "**************************************"
    write(NB6,*)
  endif

  if(allocated(Tv)) deallocate(Tv)
  if(allocated(At)) deallocate(At)
  allocate(Tv(NOV,maxdav))
  allocate(At(NOV,maxdav))
  allocate(Zr(NOV,nz))
  call vecinit(Tv,NOV*maxdav,0.d0)
  call vecinit(At,NOV*maxdav,0.d0)
  call vecinit(Zr,NOV*nz,0.d0)

end subroutine z_nac_init

!======================================================== 
subroutine z_nac_exit
  if(allocated(Tv)) deallocate(Tv)
  if(allocated(At)) deallocate(At)
  if(allocated(Zr)) deallocate(Zr)
end subroutine z_nac_exit

!======================================================== 
subroutine z_nac_rhs(V)
  implicit real*8 (a-h,o-z)
  real*8 V(*)
  real*8, dimension(:,:,:), allocatable :: P,R,F,K
  
  allocate(P(n2,nz,nden),  F(n2,nz,nden))
  allocate(R(n2,nex,nden), K(n2,nex,nden))

  call get_nac_density(P,P(1,1,nden),Xv,V(jCa),V(jCb),nz) 
  if(iprtz.ge.5) then
    write(nb6,*) "nac_density"
    call matprnt(P,n2,nz*nden,8)
  endif

  call emake_rao(R,R(1,1,nden),Xv,V(jCa),V(jCb),nex)
  if(iprtz.ge.5) then
    write(nb6,*) "transition density"
    call matprnt(R,n2,nex*nden,8)
  endif

  call makefao(F,F(1,1,nden),P,P(1,1,nden),V(jW),nz,NBas,LM6,nden,1,cfacJ,cfacK)
  if(iprtz.ge.5) then
    write(nb6,*) "F[p]"
    call matprnt(F,n2,nz*nden,8)
  endif

  call makefao(K,K(1,1,nden),R,R(1,1,nden),V(jW),nex,NBas,LM6,nden,I1n,cfacJ,cfacK)
  if(iprtz.ge.5) then
    write(nb6,*) "K[R]"
    call matprnt(K,n2,nex*nden,8)
  endif
 
  m = 1
  do i=1,nex
    do j=i,nex
      call z_nac_rhs_mn(V(jCa),V(jCb),F(1,m,1),F(1,m,nden),K(1,i,1),K(1,i,nden),   &
                 K(1,j,1),K(1,j,nden),Xv(1,i),Xv(1+NOVa,i),        &
                 Xv(1,j), Xv(1+NOVa,j),Zr(1,m),Zr(1+NOVa,m))
      if(do_nac) then
        ! in case of ground state 
        if(imult.le.2 .and. i.ne.j .and. (Es(i,1).eq.0.d0 .or. Es(j,1).eq.0.d0)) then
          if(Es(i,1) .gt. Es(j,1)) istate = i
          if(Es(i,1) .lt. Es(j,1)) istate = j
          de = 0.5d0*(Es(j,1)-Es(i,1))
          call vecadd3(Zr(1,m),de,Xv(1,istate),NOV)
        endif
      endif
      m = m + 1
    enddo
  enddo

  if(iprtz.ge.5) then
    write(nb6,*) "r.h.s in cis_nac"
    call matprnt(Zr,NOV,nz,6)
  endif


  deallocate(P)
  deallocate(F)
  deallocate(R)
  deallocate(K)
  return
end subroutine z_nac_rhs

!======================================================== 
subroutine get_nac_density(Pa,Pb,R,Ca,Cb,n)
  implicit none
  integer i,j,n,m,istate,jstate
  real*8 Pa(n2,*),Pb(n2,*),R(NOV,*),Ca(LM2,*),Cb(LM2,*) 

  call vecinit(Pa,n2*n,0.d0)
  if(nden.eq.2) call vecinit(Pb,n2*n,0d0)

  !Pao = 0.5 * [Cv * (XI * XJ^T + XJ * XI^T) * Cv - Co * (XI^T * XJ + XJ^T * XI) * Co]
  m = 1
  do i=1,nex
    do j=i,nex
      call get_cis_density_mn(Pa(1,m),Pb(1,m),R(1,i),R(1+NOVa,i),  &
             R(1,j),R(1+NOVa,j),xx,xx,xx,xx,Ca,Cb,1,0,0)
      m = m + 1
    enddo
  enddo
  
  return
end subroutine get_nac_density

!======================================================== 
subroutine get_cis_density_mn(Pa,Pb,RIa,RIb,RJa,RJb,Za,Zb,Pa0,Pb0,Ca,Cb,iunrelax,irelax,igs)
  implicit none
  integer iunrelax,irelax,igs,n2
  real*8 Pa(*),Pb(*),RIa(*),RIb(*),RJa(*),RJb(*),Ca(LM2,*),Cb(LM2,*)
  real*8 Za(*),Zb(*),Pa0(*),Pb0(*)
  real*8, dimension(:), allocatable :: tmp,tmp2,tmp3

  n2 = nbas*nbas
  allocate(tmp(n2),tmp2(n2),tmp3(n2))
  call vecinit(Pa,n2,0.d0)

  if(iunrelax.eq.1) then
    ! RI * RJ^t
    if(imult.le.2) then
      call matmult(Ca(1,1+NOa),RIa,tmp,NBas,NOa,NVa,LM2,NVa,NBas,1)
      call matmult(Ca(1,1+NOa),RJa,tmp2,NBas,NOa,NVa,LM2,NVa,NBas,1)
      call matmult(tmp,tmp2,Pa,NBas,NBas,NOa,NBas,NBas,NBas,3)

      ! RI^t * RJ
      call matmult(Ca(1,1),RIa,tmp,NBas,NVa,NOa,LM2,NVa,NBas,3)
      call matmult(Ca(1,1),RJa,tmp2,NBas,NVa,NOa,LM2,NVa,NBas,3)
      call matmult(tmp,tmp2,tmp3,NBas,NBas,NVa,NBas,NBas,NBas,3)
      call vecsub(Pa,Pa,tmp3,n2)

    else
      call matmult(Ca(1,1),RIa,tmp,NBas,NVb,NOa,LM2,NVb,NBas,3)
      call matmult(Ca(1,1),RJa,tmp2,NBas,NVb,NOa,LM2,NVb,NBas,3)
      call matmult(tmp,tmp2,tmp3,NBas,NBas,NVb,NBas,NBas,NBas,3)
      call vecsub(Pa,Pa,tmp3,n2)
    endif
    call addmatsym(Pa,NBas)
    call vecscale(Pa,n2,0.5d0)
  endif

  if(irelax.eq.1) then
    call matmult(Ca(1,1+NOa),Za,tmp,NBas,NOa,NVa,LM2,NVa,NBas,1)
    call matmult(tmp,Ca(1,1),tmp2,NBas,NBas,NOa,NBas,LM2,NBas,3)
    call addmatsym(tmp2,NBas)
    call vecadd(Pa,Pa,tmp2,NBas*NBas)
  endif

  if(igs.eq.1) then
    call square(Pa0,tmp,NBas,NBas,LM4)
    call vecadd2(Pa,0.5d0,tmp,Pa,NBas*NBas)
  endif

  if(nden.eq.2) then
    call vecinit(Pb,n2,0.d0)

    if(iunrelax.eq.1) then
      if(imult.le.2) then
        call matmult(Cb(1,1+NOb),RIb,tmp, NBas,NOb,NVb,LM2,NVb,NBas,1)
        call matmult(Cb(1,1+NOb),RJb,tmp2,NBas,NOb,NVb,LM2,NVb,NBas,1)
        call matmult(tmp,tmp2,Pb,NBas,NBas,NOb,NBas,NBas,NBas,3)

        call matmult(Cb(1,1),RIb,tmp, NBas,NVb,NOb,LM2,NVb,NBas,3)
        call matmult(Cb(1,1),RJb,tmp2,NBas,NVb,NOb,LM2,NVb,NBas,3)
        call matmult(tmp,tmp2,tmp3,NBas,NBas,NVb,NBas,NBas,NBas,3)
        call vecsub(Pb,Pb,tmp3,n2)
      else
        call matmult(Cb(1,1+NOb),RIa,tmp,NBas,NOa,NVb,LM2,NVb,NBas,1)
        call matmult(Cb(1,1+NOb),RJa,tmp2,NBas,NOa,NVb,LM2,NVb,NBas,1)
        call matmult(tmp,tmp2,Pb,NBas,NBas,NOa,NBas,NBas,NBas,3)
      endif
      call addmatsym(Pb,NBas)
      call vecscale(Pb,n2,0.5d0)
    endif

    if(irelax.eq.1) then
      call matmult(Cb(1,1+NOb),Zb,tmp,NBas,NOb,NVb,LM2,NVb,NBas,1)
      call matmult(tmp,Cb(1,1),tmp2,NBas,NBas,NOb,NBas,LM2,NBas,3)
      call addmatsym(tmp2,NBas)
      call vecadd(Pb,Pb,tmp2,NBas*NBas)
    endif

    if(igs.eq.1) then
      call square(Pb0,tmp,NBas,NBas,LM4)
      call vecadd2(Pb,0.5d0,tmp,Pb,NBas*NBas)
    endif

  endif

  deallocate(tmp)
  deallocate(tmp2)
  deallocate(tmp3)
  return
end subroutine get_cis_density_mn

!======================================================== 
subroutine z_nac_rhs_0n(istate,n)
  implicit none
  integer istate,n
  real*8 de

  de = 0.5d0*Es(istate,1) 
  call veccopy(Zr(1,n),Xv(1,istate),NOV)
  call vecscale(Zr(1,n),NOV,de)

  return
end subroutine z_nac_rhs_0n

!======================================================== 
subroutine z_nac_rhs_mn(Ca,Cb,Fa,Fb,KIa,KIb,KJa,KJb,RIa,RIb,RJa,RJb,Za,Zb)
  implicit none
  real*8 Ca(LM2,*),Cb(LM2,*),Fa(*),Fb(*),KIa(*),KIb(*),KJa(*),KJb(*)
  real*8 RIa(*),RIb(*),RJa(*),RJb(*),Za(*),Zb(*)
  real*8, dimension(:), allocatable :: tmp,tmp2

  allocate(tmp(n2),tmp2(n2))

  call vecinit(Za,NOVa,0.d0)  

  ! Cv^t * G[P] * Co
  call matmult(Ca(1,1+NOa),Fa,tmp,NVa,NBas,NBas,LM2,NBas,NVa,2)
  call matmult(tmp,Ca,tmp2,NVa,NOa,NBas,NVa,LM2,NVa,1)
  call vecsub(Za,Za,tmp2,NVa*NOa)

  if(imult.le.2) then
    ! Cv^t * G[R_I^t] * Cv * X_J
    call matmult(Ca(1,1+NOa),KIa,tmp,NVa,NBas,NBas,LM2,NBas,NVa,4)
    call matmult(tmp,Ca(1,1+NOa),tmp2,NVa,NVa,NBas,NVa,LM2,NVa,1)
    call matmult(tmp2,RJa,tmp,NVa,NOa,NVa,NVa,NVa,NVa,1)
    !call vecadd2(Za,-1.d0,tmp,Za,NOVa) 
    call vecadd2(Za,-0.5d0,tmp,Za,NVa*NOa) 
 
    ! X_J * Co^t * G[R_I^t] * Co
    call matmult(Ca(1,1),KIa,tmp,NOa,NBas,NBas,LM2,NBas,NOa,4)
    call matmult(tmp,Ca(1,1),tmp2,NOa,NOa,NBas,NOa,LM2,NOa,1)
    call matmult(RJa,tmp2,tmp,NVa,NOa,NOa,NVa,NOa,NVa,1)
    !call vecadd2(Za,1.d0,tmp,Za,NOVa)
    call vecadd2(Za,0.5d0,tmp,Za,NVa*NOa)
 
    ! Cv^t * G[R_J^t] * Cv * X_I
    call matmult(Ca(1,1+NOa),KJa,tmp,NVa,NBas,NBas,LM2,NBas,NVa,4)
    call matmult(tmp,Ca(1,1+NOa),tmp2,NVa,NVa,NBas,NVa,LM2,NVa,1)
    call matmult(tmp2,RIa,tmp,NVa,NOa,NVa,NVa,NVa,NVa,1)
    call vecadd2(Za,-0.5d0,tmp,Za,NVa*NOa) 
 
    ! X_I * Co^t * G[R_J^t] * Co
    call matmult(Ca(1,1),KJa,tmp,NOa,NBas,NBas,LM2,NBas,NOa,4)
    call matmult(tmp,Ca(1,1),tmp2,NOa,NOa,NBas,NOa,LM2,NOa,1)
    call matmult(RIa,tmp2,tmp,NVa,NOa,NOa,NVa,NOa,NVa,1)
    call vecadd2(Za,0.5d0,tmp,Za,NVa*NOa)
  else
    call matmult(Ca(1,1+NOa),KIa,tmp,NVa,NBas,NBas,LM2,NBas,NVa,4)
    call matmult(tmp,Cb(1,1+NOb),tmp2,NVa,NVb,NBas,NVa,LM2,NVa,1)
    call matmult(tmp2,RJa,tmp,NVa,NOa,NVb,NVa,NVb,NVa,1)
    call vecadd2(Za,-0.5d0,tmp,Za,NVa*NOa) 

    call matmult(Ca(1,1+NOa),KJa,tmp,NVa,NBas,NBas,LM2,NBas,NVa,4)
    call matmult(tmp,Cb(1,1+NOb),tmp2,NVa,NVb,NBas,NVa,LM2,NVa,1)
    call matmult(tmp2,RIa,tmp,NVa,NOa,NVb,NVa,NVb,NVa,1)
    call vecadd2(Za,-0.5d0,tmp,Za,NVa*NOa) 
  endif

  if(nden.eq.2) then
    call vecinit(Zb,NOVb,0d0)

    ! Cv^t * G[P] * Co
    call matmult(Cb(1,1+NOb),Fb,tmp,NVb,NBas,NBas,LM2,NBas,NVb,2)
    call matmult(tmp,Cb,tmp2,NVb,NOb,NBas,NVb,LM2,NVb,1)
    call vecsub(Zb,Zb,tmp2,NVb*NOb)
 
    if(imult.le.2) then
      ! Cv^t * G[R_I^t] * Cv * X_J
      call matmult(Cb(1,1+NOb),KIb,tmp,NVb,NBas,NBas,LM2,NBas,NVb,4)
      call matmult(tmp,Cb(1,1+NOb),tmp2,NVb,NVb,NBas,NVb,LM2,NVb,1)
      call matmult(tmp2,RJb,tmp,NVb,NOb,NVb,NVb,NVb,NVb,1)
      call vecadd2(Zb,-0.5d0,tmp,Zb,NVb*NOb) 
     
      ! X_J * Co^t * G[R_I^t] * Co
      call matmult(Cb(1,1),KIb,tmp,NOb,NBas,NBas,LM2,NBas,NOb,4)
      call matmult(tmp,Cb(1,1),tmp2,NOb,NOb,NBas,NOb,LM2,NOb,1)
      call matmult(RJb,tmp2,tmp,NVb,NOb,NOb,NVb,NOb,NVb,1)
      call vecadd2(Zb,0.5d0,tmp,Zb,NVb*NOb)
     
      ! Cv^t * G[R_J^t] * Cv * X_I
      call matmult(Cb(1,1+NOb),KJb,tmp,NVb,NBas,NBas,LM2,NBas,NVb,4)
      call matmult(tmp,Cb(1,1+NOb),tmp2,NVb,NVb,NBas,NVb,LM2,NVb,1)
      call matmult(tmp2,RIb,tmp,NVb,NOb,NVb,NVb,NVb,NVb,1)
      call vecadd2(Zb,-0.5d0,tmp,Zb,NVb*NOb) 
     
      ! X_I * Co^t * G[R_J^t] * Co
      call matmult(Cb(1,1),KJb,tmp,NOb,NBas,NBas,LM2,NBas,NOb,4)
      call matmult(tmp,Cb(1,1),tmp2,NOb,NOb,NBas,NOb,LM2,NOb,1)
      call matmult(RIb,tmp2,tmp,NVb,NOb,NOb,NVb,NOb,NVb,1)
      call vecadd2(Zb,0.5d0,tmp,Zb,NVb*NOb)
    else
      ! X_I * Co^t * G[R_I^t] * Co
      call matmult(Ca(1,1),KIa,tmp,NOa,NBas,NBas,LM2,NBas,NOa,4)
      call matmult(tmp,Cb(1,1),tmp2,NOa,NOb,NBas,NOa,LM2,NOa,1)
      call matmult(RJa,tmp2,tmp,NVb,NOb,NOa,NVb,NOa,NVb,1)
      call vecadd2(Zb,0.5d0,tmp,Zb,NVb*NOb)

      ! X_I * Co^t * G[R_I^t] * Co
      call matmult(Ca(1,1),KJa,tmp,NOa,NBas,NBas,LM2,NBas,NOa,4)
      call matmult(tmp,Cb(1,1),tmp2,NOa,NOb,NBas,NOa,LM2,NOa,1)
      call matmult(RIa,tmp2,tmp,NVb,NOb,NOa,NVb,NOa,NVb,1)
      call vecadd2(Zb,0.5d0,tmp,Zb,NVb*NOb)
    endif
  else 
    call veccopy(Zb,Za,NVa*NOa)
  endif

  deallocate(tmp)
  deallocate(tmp2)
  return
end subroutine z_nac_rhs_mn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! iterative solving Z-Vector equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_iter(V)
  implicit none
  logical success
  real*8 V(*)
  real*8, dimension(:), allocatable :: H

  niter = 0
  success = .false.

  ! initial guess
  allocate(H(NOV))
  call z_iter_init(V,Tv,Zr,H,H(1+NOVa),nz)

  ! Enter the iterative diagonalizer
  do while ((.not.success).and.(niter.le.kitdav))
    niter = niter + 1
    call z_make_a(V,H)
    call z_solve(H,NOV)
    if (nl.eq.0) success=.true.
  enddo
  if (.not.success) then
     write(NB6,*) "Max Davidson Iteration: ", kitdav
     STOP
  endif

  !call vecinit(Zv,NOV*nz,0d0)
  if(iprtz.ge.6) then
    ns = nz
    nl = nz
    write(6,*) "Z vector"
    call matprnt(Zv,NOV,nz,6)
    call veccopy(Tv,Zv,NOV*nz)
    call z_make_a(V,H)
    write(6,*) "(A+B)*Z"
    call vecsub(At,At,Zr,NOV*nz)
    call matprnt(At,NOV,nz,6)
  endif

  deallocate(H)
  return
end subroutine z_iter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! trial vectors for Z vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_iter_init(V,G,B,Ha,Hb,n)
  implicit none
  integer i,j,k,m,n
  real*8 dot,tol
  real*8 V(*),G(NOV,*),Ha(*),Hb(*),B(NOV,*)

  call vecinit(Ha,NOVa,0.d0)
  call vecinit(Hb,NOVb,0.d0)
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

  tol = NOV*4d-4*1d1**(2*(1-NrmDav))
  call vecinit(G,NOV*n,0.d0)
  ns = 0
  do m=1,n
    do i=1,NOa*NVa
      G(i,ns+1) = B(i,m) / Ha(i)
    enddo
    do i=1,NOb*NVb
      G(i+NOVa,ns+1) = B(i+NOVa,m) / Hb(i)
    enddo
    do i=1,ns
      call vecdot(dot,G(1,ns+1),G(1,i),NOV)
      call vecadd2(G(1,ns+1),-dot,G(1,i),G(1,ns+1),NOV)
    enddo
    call vecdot(dot,G(1,ns+1),G(1,ns+1),NOV)
    if(dot.ge.tol) then
      call vecscale(G(1,ns+1),NOV,1/sqrt(dot))
      ns = ns + 1
    endif
  enddo
  if(ns.eq.0) then
    ns = 1
    if(dot.lt.1.d-8) G(1,1) = 1.d0
  endif
  nl = ns
  if(iprtz.ge.5) then
    write(NB6,*) "inital guess in iterz"
    call matprnt(G,NOV,ns,6)
  endif

  return
end subroutine z_iter_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! build A matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_make_a(V,H)
  implicit none
  integer lb
  real*8  V(*),H(*)
  real*8, dimension(:,:), allocatable :: P, F

  allocate(P(NBas*NBas*nl,nden))
  allocate(F(NBas*NBas*nl,nden))

  call z_make_r(P,P(1,nden),Tv(1,ns-nl+1),V(jCa),V(jCb),nl)
  if(iprtz.ge.6) then
    write(NB6,*) "trial vectors in iterz"
    call matprnt(P,NBas*NBas,nl,6)
  endif

  call makefao(F,F(1,nden),P,P(1,nden),V(jW),nl,NBas,LM6,nden,1,cfacJ,cfacK)
  !call makejk(F,F(1,nden),F2,F2(1,nden),P,P(1,nden),V(jW),nl,NBas,LM6,nden,1,cfacJ,cfacK)
  !call vecadd(F,F,F2,nbas*nbas*nl*nden)
  if(iprtz.ge.6) then
    write(NB6,*) "fock-like matrix in iterz"
    call matprnt(F,NBas*NBas,nl,6)
  endif

  call z_make_amo(At(1,ns-nl+1),F,F(1,nden),Tv(1,ns-nl+1),H,H(1+NOVa),V(jCa),V(jCb),nl)
  if(iprtz.ge.6) then
    write(NB6,*) "A in iterz"
    call matprnt(At(1,ns-nl+1),NOV,nl,6)
  endif

  deallocate(P)
  deallocate(F)
  return
end subroutine z_make_a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A matrix in mo basis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_make_amo(Az,Fa,Fb,R,Ha,Hb,Ca,Cb,n)
  implicit none
  integer n,m,i
  real*8  Az(NOV,*),Fa(NBas*NBas,*),Fb(NBas*NBas,*),R(NOV,*)
  real*8  Ha(*),Hb(*),Ca(LM2,*),Cb(LM2,*)
  real*8, dimension(:), allocatable :: tmp,tmp2

  allocate(tmp(n2),tmp2(n2))
  do m=1,n
    ! A = Cv^T * F * Co
    call matmult(Ca(1,NOa+1),Fa(1,m),tmp,NVa,NBas,NBas,LM2,NBas,NVa,2)
    call matmult(tmp,Ca(1,1),Az(1,m),NVa,NOa,NBas,NVa,LM2,NVa,1)
    call matmult(Ca(1,NOa+1),Fa(1,m),tmp,NVa,NBas,NBas,LM2,NBas,NVa,4)
    call matmult(tmp,Ca(1,1),tmp2,NVa,NOa,NBas,NVa,LM2,NVa,1)
    call vecadd(Az(1,m),Az(1,m),tmp2,NOa*NVa)
    if(nden.eq.2) then
      call matmult(Cb(1,NOb+1),Fb(1,m),tmp,NVb,NBas,NBas,LM2,NBas,NVb,2)
      call matmult(tmp,Cb(1,1),Az(1+NOVa,m),NVb,NOb,NBas,NVb,LM2,NVb,1)
      call matmult(Cb(1,NOb+1),Fb(1,m),tmp,NVb,NBas,NBas,LM2,NBas,NVb,4)
      call matmult(tmp,Cb(1,1),tmp2,NVb,NOb,NBas,NVb,LM2,NVb,1)
      call vecadd(Az(1+NOVa,m),Az(1+NOVa,m),tmp2,NOb*NVb)
    endif 

    do i=1,NOVa
      Az(i,m) = Az(i,m) + R(i,m)*Ha(i)
    enddo
    if(nden.eq.2) then
      do i=1,NOVb
        Az(i+NOVa,m) = Az(i+NOVa,m) + R(i+NOVa,m)*Hb(i)
      enddo
    else
      call veccopy(Az(1+NOVa,m),Az(1,m),NOVa)
    endif
  enddo
  deallocate(tmp)
  deallocate(tmp2)
   
  return
end subroutine z_make_amo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate AO trial density 
! R_I = Cv * X_I * Co^t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_make_r(Pa,Pb,R,Ca,Cb,n)
  implicit none
  integer n,i
  real*8  Pa(NBas*NBas,*),Pb(NBas*NBas,*),R(NOV,*),Ca(LM2,*),Cb(LM2,*)
  real*8, dimension(:), allocatable :: tmp

  allocate(tmp(NBas*NBas))

  do i=1,n
    call matmult(Ca(1,1+NOa),R(1,i),tmp,NBas,NOa,NVa,LM2,NVa,NBas,1)
    call matmult(tmp,Ca(1,1),Pa(1,i),NBas,NBas,NOa,NBas,LM2,NBas,3)
    if(nden.eq.2) then
      call matmult(Cb(1,1+NOb),R(1+NOVa,i),tmp,NBas,NOb,NVb,LM2,NVb,NBas,1)
      call matmult(tmp,Cb(1,1),Pb(1,i),NBas,NBas,NOb,NBas,LM2,NBas,3)
    endif
  enddo

  deallocate(tmp)
  return
end subroutine z_make_r

end module cis_z
