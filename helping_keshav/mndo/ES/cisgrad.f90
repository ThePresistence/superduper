!
! Created by Jie Liu on 05/16
! 
! cisgrad calculate the gradient for the specified state
!
!======================================================== 
! Xv	        transition density
! Tv	        trial vectors
! At            A matrix
! Zr            Lagrangian (r.h.s) in Z-vector equation
! Zv            Z vectors
!======================================================== 

module cis_gradient

use es_global
use cis_solve
use cisz_solve

implicit none

private

integer iprint
real*8  ecis,enuc

public  cisgrad, z_iter, z_make_f
public  get_cis_density_nn

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cisgrad(V)
  implicit none
  real*8 V(*)
  
  write(NB6,*) 
  write(NB6,*) "******************************************"
  write(NB6,*) "*****    cis gradient calculation   ******"
  write(NB6,*) "******************************************"
  write(NB6,*) 
  call cis_grad(V)

  return
end subroutine cisgrad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main cis gradient calculation subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cis_grad(V)
  implicit none
  real*8 V(*)

  call cis_grad_init()

  call cisz(V)

  call cis_n_grad(V)

  call cis_grad_exit()
  return
end subroutine cis_grad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initionalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cis_grad_init()
  use LIMIT, only : LEN
  implicit none
  integer IN2
  real*8  E
  common /INOPT2/ IN2(300)
  common /ENERGT/ E(2)

  nz = nroots
  if(imult.eq.3) nden = 2
   
  ecis = E(1)
  enuc = E(2)

  iprint = IN2(41)

  allocate(Zr(NOV,nz))
  call vecinit(Zr,NOV*nz,0d0)
  if(allocated(At)) deallocate(At)
  if(allocated(Tv)) deallocate(Tv)
  allocate(At(NOV,maxdav))
  allocate(Tv(NOV,maxdav))

end subroutine cis_grad_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! clean up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cis_grad_exit()
  implicit none
  real*8 E
  common /ENERGT/ E(2)
  E(1) = ecis
  E(2) = enuc
  if(allocated(Zr)) deallocate(Zr)
  if(allocated(Zv)) deallocate(Zv)
  if(allocated(At)) deallocate(At)
  if(allocated(Tv)) deallocate(Tv)
end subroutine cis_grad_exit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! build cis gradient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cis_n_grad(V)
  use LIMIT
  implicit none
  integer irelax,igs,i,j,N2,io,imode,natom,NUMAT
  real*8 V(*), CG
  real*8, dimension(:,:), allocatable :: P
  common /CGRAD / CG(3,LM1+LM1M)
  common /ATOMS / NUMAT

  if(iprint.gt.0) then
    write(NB6,*)
    write(NB6,*) "     Let's Build CIS Gradient      "
    write(NB6,*)
  endif

  allocate(P(NBas*NBas,3*nden))
  call vecinit(P,NBas*NBas*3*nden,0d0)

  irelax = 1
  igs = 1 
  ! P_I * H^x + W_I * S^x
  call cis_density_n(P,P(1,nden),V,lstate,irelax,igs)
  if(iprint.eq.6) then
    write(NB6,*) "cis density"
    call matprnt(P,NBas,NBas,6)
  endif

  call emake_rao(P(1,nden+1),P(1,2*nden),Xv(1,lstate),V(jCa),V(jCb),1) 
  if(iprint.eq.6) then
    write(NB6,*) "cis transition density"
    call matprnt(P(1,nden+1),NBas,NBas,6)
  endif

  call square(V(jPa),P(1,2*nden+1),NBas,NBas,LM4)
  if(nden.eq.2) call square(V(jPb),P(1,3*nden),NBas,NBas,LM4)

  imode = 1
  natom = NUMAT
  do i=1,natom
    do j=1,3
      call ee_cis_gradient(imode,i,j,V,LM5,CG(j,i),P,nden,iprint,NB6,I1n)
    enddo
  enddo
  if(iprint.gt.0) then
    write(NB6,*) "total gradient"
    call matprnt(CG,3,natom,6)
  endif
  deallocate(P)
end subroutine cis_n_grad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! density of cis state
! P = Cv X X^t Cv^t - Co X^t X Co^t + 1/2 Cv Z Co^t + 1/2 Co Z^t Cv^t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_cis_density_nn(Pa,Pb,Ra,Rb,Za,Zb,Pa0,Pb0,Ca,Cb,irelax,igs)
  implicit none
  integer irelax,igs
  real*8 Pa(*),Pb(*),Ra(*),Rb(*),Ca(LM2,*),Cb(LM2,*),Za(*),Zb(*),Pa0(*),Pb0(*)
  real*8, dimension(:), allocatable:: tmp, tmp2

  allocate(tmp(NBas*NBas),tmp2(NBas*NBas))

  if(imult.le.2) then
    call matmult(Ca(1,1+NOa),Ra,tmp,NBas,NOa,NVa,LM2,NVa,NBas,1)
    call matmult(tmp,tmp,Pa,NBas,NBas,NOa,NBas,NBas,NBas,3)

    call matmult(Ca(1,1),Ra,tmp,NBas,NVa,NOa,LM2,NVa,NBas,3)
    call matmult(tmp,tmp,tmp2,NBas,NBas,NVa,NBas,NBas,NBas,3)

    call vecsub(Pa,Pa,tmp2,NBas*NBas)
  else if(imult.eq.3) then
    call matmult(Ca(1,1),Ra,tmp,NBas,NVb,NOa,LM2,NVa,NBas,3)
    call matmult(tmp,tmp,Pa,NBas,NBas,NVb,NBas,NBas,NBas,3)
    call vecscale(Pa,NBas*NBas,-1.d0)
  endif

  if(irelax.eq.1) then
    call matmult(Ca(1,1+NOa),Za,tmp,NBas,NOa,NVa,LM2,NVa,NBas,1)
    call matmult(tmp,Ca(1,1),tmp2,NBas,NBas,NOa,NBas,LM2,NBas,3)
    call addmatsym(tmp2,NBas)
    call vecadd(Pa,Pa,tmp2,NBas*NBas)
  endif

  if(igs.eq.1) then
    ! 1/2 P0 II P0 (ground state)
    call square(Pa0,tmp,NBas,NBas,LM4)
    call vecadd2(Pa,0.5d0,tmp,Pa,NBas*NBas)
  endif

  if(nden.eq.2) then
    if(imult.le.2) then
      call matmult(Cb(1,1+NOb),Rb,tmp,NBas,NOb,NVb,LM2,NVb,NBas,1)
      call matmult(tmp,tmp,Pb,NBas,NBas,NOb,NBas,NBas,NBas,3)
     
      call matmult(Cb(1,1),Rb,tmp,NBas,NVb,NOb,LM2,NVb,NBas,3)
      call matmult(tmp,tmp,tmp2,NBas,NBas,NVb,NBas,NBas,NBas,3)
     
      call vecsub(Pb,Pb,tmp2,NBas*NBas)
    else if(imult.eq.3) then
      call matmult(Cb(1,1+NOb),Ra,tmp,NBas,NOa,NVb,LM2,NVb,NBas,1)
      call matmult(tmp,tmp,Pb,NBas,NBas,NOa,NBas,NBas,NBas,3)
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
  return
end subroutine get_cis_density_nn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! n-th cis density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cis_density_n(Pa,Pb,V,istate,irelax,igs)
  implicit none
  integer istate,irelax,igs,io
  real*8 Pa(*),Pb(*),V(*)

  call get_cis_density_nn(Pa,Pb,Xv(1,istate),Xv(1+NOVa,istate),Zv(1,istate),Zv(1+NOVa,istate),V(jPa),V(jPb),V(jCa),V(jCb),irelax,igs)

  return
end subroutine cis_density_n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! cis transition density 
! P = Cv X Co^t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_cis_transition_density(Pa,Pb,R,Ca,Cb)
  implicit none
  real*8 Pa(NBas*NBas),Pb(NBas*NBas),R(*),Ca(LM2,*),Cb(LM2,*)

  call emake_rao(Pa,Pb,R,Ca,Cb,1)

  return
end subroutine get_cis_transition_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! n-th cis transition density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cis_transition_density_n(Pa,Pb,V,istate)
  implicit none
  integer istate
  real*8  Pa(*),Pb(*),V(*)

  call get_cis_transition_density(Pa,Pb,Xv(1,istate),V(jCa),V(jCb))

  return
end subroutine cis_transition_density_n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate Z vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cisz(V)
  implicit none
  real*8 V(*)

  if(iprint.gt.0) then
    write(NB6,*) 
    write(NB6,*) "**************************************"
    write(NB6,*) "*               Iterz                *"  
    write(NB6,*) "**************************************"
    write(NB6,*) 
  endif

  call z_rhs(V)

  call z_iter(V)

  return  
end subroutine cisz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate the Lagrangian L = - Cv^t * G[P_I] * Co - Cv^t * G[R_I^t] * Cv * X_I
!                              + X_I * Co^t * G[R_I^t] * Co
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rhs(V)
  implicit none
  real*8  V(*)
 
  call z_rhs_p(V(jCa),V(jCb),V(jW),Zr)

  call z_rhs_r(V(jCa),V(jCb),V(jW),Zr)

  if(iprint.eq.6) then 
    write(NB6,*) "Lagrangian in Z-vector equation"
    call matprnt(Zr,NOV,nz,6)
  endif

  return
end subroutine z_rhs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! L_p = - Cv^t * G[P_I] * Co
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rhs_p(Ca,Cb,W,B)
  implicit none
  real*8 Ca(*),Cb(*),W(*),B(*)
  real*8, dimension(:), allocatable :: P, F 
  integer lb,ndim

  allocate(P(NBas*NBas*nz*nden))
  allocate(F(NBas*NBas*nz*nden))

  lb = NBas*NBas*nz*(nden-1) + 1 

  call z_make_p(P,P(lb),Xv,Ca,Cb,nz)
  if(iprint.eq.6) then
    write(NB6,*) "unrelaxed difference density matrix"
    call matprnt(P,NBas*NBas,nz,6)
  endif

  call z_make_f(F,F(lb),P,P(lb),W,nz)
  if(iprint.eq.6) then
    write(NB6,*) "fock-like matrix"
    call matprnt(F,NBas*NBas,nz,6)
  endif

  call z_make_lp(F,F(lb),P,Ca,Cb,nz)
  if(iprint.eq.6) then 
    write(NB6,*) "Lagrangian P in Z-vector equation"
    call matprnt(P,NOV,nz,6)
  endif

  call vecsub(B,B,P,NOV*nz)

  deallocate(P)
  deallocate(F)
  return
end subroutine z_rhs_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! L_r = - Cv^t * G[R_I^t] * Cv * X_I + X_I * Co^t * G[R_I^t] * Co
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rhs_r(Ca,Cb,W,B)
  implicit none
  real*8 Ca(*),Cb(*),W(*),B(*)
  real*8, dimension(:), allocatable :: P, F 
  integer lb

  allocate(P(NBas*NBas*nz*nden))
  allocate(F(NBas*NBas*nz*nden))

  lb = NBas*NBas*nz*(nden-1) + 1 

  call emake_rao(P,P(lb),Xv,Ca,Cb,nz)
  if(iprint.eq.6) then
    write(NB6,*) "transition density matrix"
    call matprnt(P,NBas*NBas,nz,6)
  endif

  call z_make_f(F,F(lb),P,P(lb),W,nz)
  if(iprint.eq.6) then
    write(NB6,*) "fock-like matrix"
    call matprnt(F,NBas*NBas,nz,6)
  endif

  call z_make_lr(F,F(lb),Xv,P,Ca,Cb,nz)
  if(iprint.eq.6) then 
    write(NB6,*) "Lagrangian R in Z-vector equation"
    call matprnt(P,NOV,nz,6)
  endif

  call vecsub(B,B,P,NOV*nz)

  deallocate(P)
  deallocate(F)
  return
end subroutine z_rhs_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate unrelaxed difference density matrix 
! P_I = Cv * X_I * X_I^t * Cv^t - Co * X_I^t * X_I * Co^t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_make_p(Pa,Pb,R,Ca,Cb,n)
  implicit none
  integer i,n
  real*8 Pa(NBas*NBas,*),Pb(NBas*NBas,*),R(NOV,*),Ca(LM2,*),Cb(LM2,*)
  real*8, dimension(:), allocatable :: tmp,tmp2

  allocate(tmp(NBas*NBas),tmp2(NBas*NBas))
  
  ! Pao = Cv * X * X^T * Cv - Co * X^T * X * Co
  if(imult.le.2) then
    do i=1,n
      call matmult(Ca(1,1+NOa),R(1,i),tmp,NBas,NOa,NVa,LM2,NVa,NBas,1)
      call matmult(tmp,tmp,Pa(1,i),NBas,NBas,NOa,NBas,NBas,NBas,3)
      
      call matmult(Ca(1,1),R(1,i),tmp,NBas,NVa,NOa,LM2,NVa,NBas,3)
      call matmult(tmp,tmp,tmp2,NBas,NBas,NVa,NBas,NBas,NBas,3)
  
      call vecsub(Pa(1,i),Pa(1,i),tmp2,NBas*NBas)
  
      if(nden.eq.2) then
        call matmult(Cb(1,1+NOb),R(1+NOVa,i),tmp,NBas,NOb,NVb,LM2,NVb,NBas,1)
        call matmult(tmp,tmp,Pb(1,i),NBas,NBas,NOb,NBas,NBas,NBas,3)
        
        call matmult(Cb(1,1),R(1+NOVa,i),tmp,NBas,NVb,NOb,LM2,NVb,NBas,3)
        call matmult(tmp,tmp,tmp2,NBas,NBas,NVb,NBas,NBas,NBas,3)
    
        call vecsub(Pb(1,i),Pb(1,i),tmp2,NBas*NBas)
      endif
    enddo
  else if(imult.eq.3) then
    do i=1,n
      call matmult(Ca(1,1),R(1,i),tmp,NBas,NVb,NOa,LM2,NVa,NBas,3)
      call matmult(tmp,tmp,Pa(1,i),NBas,NBas,NVb,NBas,NBas,NBas,3)
      call vecscale(Pa(1,i),NBas*NBas,-1.d0)
 
      call matmult(Cb(1,1+NOb),R(1,i),tmp,NBas,NOa,NVb,LM2,NVb,NBas,1)
      call matmult(tmp,tmp,Pb(1,i),NBas,NBas,NOa,NBas,NBas,NBas,3)
    enddo
  endif

  deallocate(tmp)
  deallocate(tmp2)
  return
end subroutine z_make_p

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate Fock-like matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_make_f(Fa,Fb,Pa,Pb,W,n)
  implicit none
  integer n
  real*8 Fa(*),Fb(*),Pa(*),Pb(*),W(*)

  call makefao(Fa,Fb,Pa,Pb,W,n,NBas,LM6,nden,I1n,1.d0,1.d0)
  
  return
end subroutine z_make_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! L_p in mo basis set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_make_lp(Fa,Fb,Lp,Ca,Cb,n)
  implicit none
  integer i,n
  real*8 Fa(NBas*NBas,*),Fb(NBas*NBas,*),Lp(NOV,*),Ca(LM2,*),Cb(LM2,*)
  real*8, dimension(:), allocatable :: tmp
  
  allocate(tmp(NBas*NBas))
  call vecinit(Lp,NOV*n,0.d0)

  do i=1,n
    call matmult(Ca(1,1+NOa),Fa(1,i),tmp,NVa,NBas,NBas,LM2,NBas,NVa,2)
    call matmult(tmp,Ca,Lp(1,i),NVa,NOa,NBas,NVa,LM2,NVa,1)
    if(nden.eq.2 .or. imult.eq.3) then
      call matmult(Cb(1,1+NOb),Fb(1,i),tmp,NVb,NBas,NBas,LM2,NBas,NVb,2)
      call matmult(tmp,Cb,Lp(1+NOVa,i),NVb,NOb,NBas,NVb,LM2,NVb,1)
    else 
      call veccopy(Lp(1+NOVa,i),Lp(1,i),NOVa)
    endif
  enddo

  deallocate(tmp)
  return
end subroutine z_make_lp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! L_r in mo basis set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_make_lr(Fa,Fb,R,Lr,Ca,Cb,n)
  implicit none
  integer i,n
  real*8 Fa(NBas*NBas,*),Fb(NBas*NBas,*),R(NOV,*)
  real*8 Lr(NOV,*),Ca(LM2,*),Cb(LM2,*)
  real*8, dimension(:), allocatable :: tmp,tmp2
  
  allocate(tmp(NBas*NBas),tmp2(NBas*NBas))
  call vecinit(Lr,NOV*n,0.d0)

  ! Cv^t * G[R_I^t] * Cv * X_I - X_I * Co^t * G[R_I^t] * Co
  if(imult.le.2) then
    do i=1,n
      ! Cv^t * G[R_I^t] * Cv * X_I
      call matmult(Ca(1,1+NOa),Fa(1,i),tmp,NVa,NBas,NBas,LM2,NBas,NVa,4)
      call matmult(tmp,Ca(1,1+NOa),tmp2,NVa,NVa,NBas,NVa,LM2,NVa,1)
      call matmult(tmp2,R(1,i),Lr(1,i),NVa,NOa,NVa,NVa,NVa,NVa,1)
 
      ! X_I * Co^t * G[R_I^t] * Co
      call matmult(Ca(1,1),Fa(1,i),tmp,NOa,NBas,NBas,LM2,NBas,NOa,4)
      call matmult(tmp,Ca(1,1),tmp2,NOa,NOa,NBas,NOa,LM2,NOa,1)
      call matmult(R(1,i),tmp2,tmp,NVa,NOa,NOa,NVa,NOa,NVa,1)
  
      call vecsub(Lr(1,i),Lr(1,i),tmp,NOVa)
 
      if(nden.eq.2) then
        call matmult(Cb(1,1+NOb),Fb(1,i),tmp,NVb,NBas,NBas,LM2,NBas,NVb,4)
        call matmult(tmp,Cb(1,1+NOb),tmp2,NVb,NVb,NBas,NVb,LM2,NVb,1)
        call matmult(tmp2,R(1+NOVa,i),Lr(1+NOVa,i),NVb,NOb,NVb,NVb,NVb,NVb,1)
       
        call matmult(Cb(1,1),Fb(1,i),tmp,NOb,NBas,NBas,LM2,NBas,NOb,4)
        call matmult(tmp,Cb(1,1),tmp2,NOb,NOb,NBas,NOb,LM2,NOb,1)
        call matmult(R(1+NOVa,i),tmp2,tmp,NVb,NOb,NOb,NVb,NOb,NVb,1)
       
        call vecsub(Lr(1+NOVa,i),Lr(1+NOVa,i),tmp,NOVb)
      else
        call veccopy(Lr(1+NOVa,i),Lr(1,i),NOVa)
      endif
    enddo  
  else
    do i=1,n
      ! Cv^t * G[R_I^t] * Cv * X_I
      call matmult(Ca(1,1+NOa),Fa(1,i),tmp,NVa,NBas,NBas,LM2,NBas,NVa,4)
      call matmult(tmp,Cb(1,1+NOb),tmp2,NVa,NVb,NBas,NVa,LM2,NVa,1)
      call matmult(tmp2,R(1,i),Lr(1,i),NVa,NOa,NVb,NVa,NVb,NVa,1)
      call vecscale(Lr(1,i),NVa*NOa,-1.d0)
 
      ! X_I * Co^t * G[R_I^t] * Co
      call matmult(Ca(1,1),Fa(1,i),tmp,NOa,NBas,NBas,LM2,NBas,NOa,4)
      call matmult(tmp,Cb(1,1),tmp2,NOa,NOb,NBas,NOa,LM2,NOa,1)
      call matmult(R(1,i),tmp2,Lr(1+NOVa,i),NVb,NOb,NOa,NVb,NOa,NVb,1)
    enddo
  endif

  deallocate(tmp)
  deallocate(tmp2)
  return
end subroutine z_make_lr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! trial vectors for Z vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_guess(G,B,Ha,Hb,n)
  implicit none 
  integer i,m,n
  real*8 dot,tol
  real*8 G(NOV,*),Ha(*),Hb(*),B(NOV,*)

  tol = NOV*4d-4*1d1**(2*(1-NrmDav))
  ns = 0
  do m=1,n
    do i=1,NOVa
      G(i,m) = B(i,m) / Ha(i)
    enddo
    do i=1,NOVb
      G(i+NOVa,m) = B(i+NOVa,m) / Hb(i)
    enddo
    do n=1,ns
      call vecdot(dot,G(1,m),G(1,n),NOV)
      call vecadd2(G(1,m),-dot,G(1,n),G(1,m),NOV)
    enddo
    call vecdot(dot,G(1,m),G(1,m),NOV)
    if(dot.ge.tol) then
      call vecscale(G(1,m),NOV,1/sqrt(dot))
      ns = ns + 1
    endif
  enddo
  nl = ns
  if(iprint.eq.6) then
    write(NB6,*) "inital guess in iterz"
    call matprnt(G,NOV,ns,6)
  endif

  return
end subroutine z_guess

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
  call eorbdif(V,H,H(1+NOVa))
  call z_guess(Tv,Zr,H,H(1+NOVa),nz)

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
  if(iprint.eq.6) then
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
! build A matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_make_a(V,H)
  implicit none
  real*8 V(*),H(*)
  real*8, dimension(:), allocatable :: P, F 
  integer lb,lz

  allocate(P(NBas*NBas*nl*nden))
  allocate(F(NBas*NBas*nl*nden))

  lb = NBas*NBas*nl*(nden-1) + 1 
  lz = (ns-nl)*NOV + 1

  call z_make_r(P,P(lb),Tv(1,ns-nl+1),V(jCa),V(jCb),nl)
  if(iprint.eq.6) then
    write(NB6,*) "trial vectors in iterz"
    call matprnt(P,NBas*NBas,nl,6)
  endif

  call z_make_f(F,F(lb),P,P(lb),V(jW),nl)
  if(iprint.eq.6) then
    write(NB6,*) "fock-like matrix in iterz"
    call matprnt(F,NBas*NBas,nl,6)
  endif

  call z_make_amo(At(1,ns-nl+1),F,F(lb),Tv(1,ns-nl+1),H,H(1+NOVa),V(jCa),V(jCb),nl)
  if(iprint.eq.6) then
    write(NB6,*) "inital guess in iterz"
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
  integer n
  real*8 Az(NOV,*),Fa(NBas*NBas,*),Fb(NBas*NBas,*),R(NOV,*)
  real*8 Ha(*),Hb(*),Ca(LM2,*),Cb(LM2,*)

  call emake_amo(Az,Fa,Fb,R,Ha,Hb,Ca,Cb,n,.true.,0,2)

  return
end subroutine z_make_amo



end module cis_gradient
