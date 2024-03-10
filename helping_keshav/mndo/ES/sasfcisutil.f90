!
! created by jie liu on 09/16
!
module sasfcis_util

use es_global

public

integer nex,n2,npr
integer js1,js2,jv
integer os1,os2,vs1,vs2,ov1,ov2,ov3,ov4
integer s11,s12,s21,s22,as1,as2
integer s1s1,s1s2,s2s1,s2s2,a1s1,a1s2,i1s1,i1s2,a1b1
integer jJos1,jJos2,jJos12,jJos21,jJvs1,jJvs2,jJvs12,jJvs21,jJov1,jJov2,jJov3,jJov4
real*8  one,two,sq2,sq3,sq6,ssq2,ssq3,ssq6,f12,f32,f14
logical do_nac,do_mecp,sref,tref
real*8, dimension(:,:), allocatable :: Js, Ks

contains

!======================================================== 
subroutine sasfcis_setup
  n2 = nbas*nbas
  
  one  = 1.d0
  two  = 2.d0
  sq2  = sqrt(two)
  sq3  = sqrt(3.d0)
  sq6  = sqrt(6.d0)
  ssq2 = sq2/two
  ssq3 = sq3/two
  ssq6 = sq6/two
  f12  = 1.d0/two
  f32  = 3.d0/two
  f14  = f12/two

  js1 = NOb + 1
  js2 = NOb + 2
  jv  = NOa + 1

  os1  = 4               ! i->s1
  os2  = os1 + NOb       ! i->s2
  vs1  = os2 + NOb       ! s1->a
  vs2  = vs1 + NVa       ! s2->a
  ov1  = vs2 + NVa       ! i->a A
  ov2  = ov1 + NVa*NOb   ! i->a B
  ov3  = ov2 + NVa*NOb   ! i->a B
  ov4  = ov3 + NVa*NOb   ! i->a B

  ! position in (NVb,NOa) matrix
  s11 = NOb*NVb           ! (s1,s1) 
  s12 = (NOb+1)*NVb       ! (s1,s2) 
  s21 = s11 + 1           ! (s2,s1)
  s22 = s12 + 1           ! (s2,s2)
  as1 = s11+2             ! (a1,s1) 
  as2 = s12+2             ! (a1,s2) 
  ! position in (NBas,NBas) matrix
  s1s1 = NOb*(NBas+1)     ! (s1,s1) 
  s2s1 = s1s1 + 1         ! (s2,s1)
  s2s2 = (NOb+1)*(NBas+1) ! (s2,s2) 
  s1s2 = s2s2 - 1         ! (s1,s2) 
  i1s1 = NOb*NBas         ! (i1,s1)
  i1s2 = (NOb+1)*NBas     ! (i1,s2) 
  a1s1 = s1s1+2           ! (a1,s1) 
  a1s2 = s1s2+2           ! (a1,s2) 
  a1b1 = NOa*(NBas+1)     ! (a1,b1)

  jJos1  = 1
  jJos2  = jJos1 + nbas*nbas
  jJos12 = jJos1 + nbas*nbas*2
  jJos21 = jJos1 + nbas*nbas*3
  jJvs1  = jJos1 + nbas*nbas*4
  jJvs2  = jJos1 + nbas*nbas*5
  jJvs12 = jJos1 + nbas*nbas*6
  jJvs21 = jJos1 + nbas*nbas*7
  jJov1  = jJos1 + nbas*nbas*8
  jJov2  = jJos1 + nbas*nbas*9
  jJov3  = jJos1 + nbas*nbas*10
  jJov4  = jJos1 + nbas*nbas*11

  sref = .false.
  tref = .false.
  if(imult.eq.0) then
    sref = .true.
  else if(imult.eq.3) then
    tref = .true.
  else
    write(nb6,*) "wrong imult for sasfcis"
    stop
  endif
  
  npr = 10
  if(sf_xcis) npr = 12

  return
end subroutine sasfcis_setup

!======================================================== 
subroutine sasfemake_sao(Pa,Ca)
  implicit none
  real*8  Pa(NBas*NBas,*),Ca(LM2,*)

  call matmult(Ca(1,js1),Ca(1,js1),Pa(1,1),NBas,NBas,1,LM2,LM2,NBas,3) !s1s1
  call matmult(Ca(1,js1),Ca(1,js2),Pa(1,2),NBas,NBas,1,LM2,LM2,NBas,3) !s1s2
  call matmult(Ca(1,js2),Ca(1,js2),Pa(1,3),NBas,NBas,1,LM2,LM2,NBas,3) !s2s2
  return
end subroutine sasfemake_sao

!======================================================== 
subroutine sasfemake_smo(Ja,Ka,Ca,njs)
  implicit none
  integer njs,i
  real*8  Ca(LM2,*),Ja(NBas*NBas,*),Ka(NBas*NBas,*)
  real*8, dimension(:), allocatable :: tmp

  allocate(tmp(NBas*NBas))
  do i=1,njs
    call matmult(Ca,Ja(1,i),tmp,NBas,NBas,NBas,LM2,NBas,NBas,2)
    call matmult(tmp,Ca,Ja(1,i),NBas,NBas,NBas,NBas,LM2,NBas,1)
    call matmult(Ca,Ka(1,i),tmp,NBas,NBas,NBas,LM2,NBas,NBas,2)
    call matmult(tmp,Ca,Ka(1,i),NBas,NBas,NBas,NBas,LM2,NBas,1)
  enddo
  deallocate(tmp)
  return
end subroutine sasfemake_smo

!======================================================== 
subroutine sasfemake_rao(Pa,R,Ca,n)
  implicit none
  integer i,n
  real*8  Pa(NBas*NBas*npr,*),R(NOV,*),Ca(LM2,*)
  real*8, dimension(:), allocatable :: tmp

  allocate(tmp(NBas*NBas))

  do i=1,n
    ! occ->s1
    call matmult(Ca,R(os1,i),tmp,NBas,1,NOb,LM2,NOb,NBas,1)
    call matmult(Ca(1,js1),tmp,Pa(jJos1 ,i),NBas,NBas,1,LM2,NBas,NBas,3)
    call matmult(Ca(1,js2),tmp,Pa(jJos12,i),NBas,NBas,1,LM2,NBas,NBas,3)

    ! occ->s2
    call matmult(Ca,R(os2,i),tmp,NBas,1,NOb,LM2,NOb,NBas,1)
    call matmult(Ca(1,js1),tmp,Pa(jJos21,i),NBas,NBas,1,LM2,NBas,NBas,3)
    call matmult(Ca(1,js2),tmp,Pa(jJos2 ,i),NBas,NBas,1,LM2,NBas,NBas,3)

    ! s1->vir
    call matmult(Ca(1,jv),R(vs1,i),tmp,NBas,1,NVa,LM2,NVa,NBas,1)
    call matmult(tmp,Ca(1,js1),Pa(jJvs1 ,i),NBas,NBas,1,NBas,LM2,NBas,3)
    call matmult(tmp,Ca(1,js2),Pa(jJvs12,i),NBas,NBas,1,NBas,LM2,NBas,3)

    ! s2-Vir
    call matmult(Ca(1,jv),R(vs2,i),tmp,NBas,1,NVa,LM2,NVa,NBas,1)
    call matmult(tmp,Ca(1,js1),Pa(jJvs21,i),NBas,NBas,1,NBas,LM2,NBas,3)
    call matmult(tmp,Ca(1,js2),Pa(jJvs2 ,i),NBas,NBas,1,NBas,LM2,NBas,3)

    ! occ-Vir
    call matmult(R(ov1,i),Ca,tmp,NVa,NBas,NOb,NVa,LM2,NVa,3)
    call matmult(Ca(1,jv),tmp,Pa(jJov1,i),NBas,NBas,NVa,LM2,NVa,NBas,1)
    call matmult(R(ov2,i),Ca,tmp,NVa,NBas,NOb,NVa,LM2,NVa,3)
    call matmult(Ca(1,jv),tmp,Pa(jJov2,i),NBas,NBas,NVa,LM2,NVa,NBas,1)
    
    ! occ-vir (s1s2,s2s1)
    if(sf_xcis) then
      call matmult(R(ov3,i),Ca,tmp,NVa,NBas,NOb,NVa,LM2,NVa,3)
      call matmult(Ca(1,jv),tmp,Pa(jJov3,i),NBas,NBas,NVa,LM2,NVa,NBas,1)
      call matmult(R(ov4,i),Ca,tmp,NVa,NBas,NOb,NVa,LM2,NVa,3)
      call matmult(Ca(1,jv),tmp,Pa(jJov4,i),NBas,NBas,NVa,LM2,NVa,NBas,1)
    endif
  enddo

  deallocate(tmp)
  return
end subroutine sasfemake_rao

!======================================================== 
subroutine ao2mo(Ja,Ca,n)
  implicit none
  integer n,i
  real*8  Ja(n2,*),Ca(LM2,*)
  real*8, dimension(:), allocatable :: tmp

  allocate(tmp(n2))
  do i=1,n
    call matmult(Ca,Ja(1,i),tmp,nbas,nbas,nbas,LM2,nbas,nbas,2)
    call matmult(tmp,Ca,Ja(1,i),nbas,nbas,nbas,nbas,LM2,nbas,1)
  enddo
  deallocate(tmp)
  return
end subroutine ao2mo

!======================================================== 
subroutine ao2mo2(Ja,Ca,n)
  implicit none
  integer n,i
  real*8  Ja(n2,*),Ca(LM2,*)
  real*8, dimension(:), allocatable :: tmp

  allocate(tmp(n2))
  do i=1,n
    call matmult(Ja(1,i),Ca,tmp,NBas,NOa,NBas,NBas,LM2,NBas,1)
    call matmult(Ca(1,js1),tmp,Ja(1,i),NVb,NOa,NBas,LM2,NBas,NVb,2)
  enddo
  deallocate(tmp)
  return
end subroutine ao2mo2

!======================================================== 
subroutine mo2ao(Pa,Ca,n)
  implicit none
  integer n,i
  real*8  Pa(n2,*),Ca(LM2,*)
  real*8, dimension(:), allocatable :: tmp

  allocate(tmp(n2))
  do i=1,n
    call matmult(Ca,Pa(1,i),tmp,nbas,nbas,nbas,LM2,nbas,nbas,1)
    call matmult(tmp,Ca,Pa(1,i),nbas,nbas,nbas,nbas,LM2,nbas,3)
  enddo
  deallocate(tmp)
  return
end subroutine mo2ao

!======================================================== 
subroutine sasfemake_pao(P,R,Ca,np,n,irelax,igs)
  implicit none 
  integer np,n,i,j,igs,irelax
  real*8  P(n2*n,*),R(NOV,*),Ca(LM2,*)
   
  call sasfemake_pmo(P,R,Ca,np,n,igs)
  if(irelax.eq.1) then
    call vecadd(P,P,Zv,n2*n)
  endif
  call mo2ao(P,Ca,np*n)

  return
end subroutine sasfemake_pao

!======================================================== 
subroutine sasfemake_pmo(P,R,Ca,np,n,igs)
  implicit none 
  integer np,n,i,j,m,igs
  real*8  P(n2*n,*),R(NOV,*),Ca(LM2,*)
  real*8, dimension(:), allocatable :: tmp

  allocate(tmp(n2*n))
  call sasfemake_pmo1_nonsym(P,tmp,R,Ca,n)
  call sasfemake_pmo2(P(1,2),P(1,3),P(1,6),P(1,4),P(1,5),P(1,7),R,Ca,n,igs)
 
  if(tref) then
    call vecadd3(P(1,2),-f12,P,n2*n)
    call vecadd3(P(1,2),f12,tmp,n2*n)
  
    call vecadd3(P(1,3),-f12,P,n2*n)
    call vecadd3(P(1,3),f12,tmp,n2*n)

    call vecadd(P,P,tmp,n2*n)
  else if(sref) then 
    call vecadd(P(1,2),P(1,2),tmp,n2*n)
    call vecsub(P(1,3),P(1,3),P,n2*n)
    call vecadd(P,P,tmp,n2*n)
    if(sf_xcis) then
      call vecsub(P(1,4),P(1,4),P,n2*n)
      call vecadd(P(1,5),P(1,5),P,n2*n)
    else
      call veccopy(P(1,4),P,n2*n)
    endif

    m = 0
    do i=1,nex
      do j=i,nex
        if(i.eq.j) then
          P(s2s2+1+m*n2,1) = P(s2s2+1+m*n2,1) + one
          P(s1s1+1+m*n2,1) = P(s1s1+1+m*n2,1) - one
       
          P(s1s1+1+m*n2,2) = P(s1s1+1+m*n2,2) - f12
          P(s2s2+1+m*n2,3) = P(s2s2+1+m*n2,3) - f12

          if(sf_xcis) then
            P(s2s2+1+m*n2,4) = P(s2s2+1+m*n2,4) - f12
            P(s1s1+1+m*n2,4) = P(s1s1+1+m*n2,4) + f12
            P(s2s2+1+m*n2,5) = P(s2s2+1+m*n2,5) + f12
            P(s1s1+1+m*n2,5) = P(s1s1+1+m*n2,5) - f12
          else 
            P(s2s2+1+m*n2,4) = P(s2s2+1+m*n2,4) + f12
            P(s1s1+1+m*n2,4) = P(s1s1+1+m*n2,4) - f12
          endif
        endif
        m = m + 1
      enddo
    enddo

  endif

  deallocate(tmp)
  return
end subroutine sasfemake_pmo

!======================================================== 
subroutine sasfemake_pao_nonsym(P,R,Ca,n)
  implicit none 
  integer n
  real*8  P(n2,n),R(NOV,*),Ca(LM2,*)
  real*8, dimension(:), allocatable :: tmp

  allocate(tmp(n2*n))
  call sasfemake_pmo1_nonsym(P,tmp,R,Ca,n)
  call vecadd(P,P,tmp,n2*n)
  call mo2ao(P,Ca,n)
  deallocate(tmp)
  return
end subroutine sasfemake_pao_nonsym

!======================================================== 
subroutine get_sasfcis_transition_density(P,R,Ca,Cb)
  implicit none 
  integer i,j,k
  real*8  P(n2,*),R(NOV,*),Ca(LM2,*),Cb(LM2,*)
  real*8, dimension(:), allocatable :: Pa,Pb

  allocate(Pa(n2),Pb(n2))
  call vecinit(P,n2,0.d0)
  k = 2
  do j=2,nroots
    call sasfemake_1e_density(Pa,Pb,R(1,1),R(1,j),Ca) 
    call vecadd(P(1,k),Pa,Pb,n2)
    k = k + 1
  enddo
  call mo2ao(P,Ca,nroots)
  
  deallocate(Pa)
  deallocate(Pb)
  return
end subroutine get_sasfcis_transition_density

!======================================================== 
subroutine get_sasfcis_density(P,R,Ca,Cb,Pa0,Pb0)
  implicit none 
  integer i,j,k
  real*8  P(n2,*),R(NOV,*),Ca(LM2,*),Cb(LM2,*),Pa0(*),Pb0(*)
  real*8, dimension(:), allocatable :: Pa,Pb

  allocate(Pa(n2),Pb(n2))

  k = 1
  do i=1,nroots
    call sasfemake_1e_density(Pa,Pb,R(1,i),R(1,i),Ca)
    call vecadd(P(1,k),Pa,Pb,n2)
    call mo2ao(P(1,k),Ca,1)
    call square(Pa0,Pa,NBas,NBas,LM4)
    call vecadd(P(1,k),P(1,k),Pa,n2)
    call square(Pb0,Pb,NBas,NBas,LM4)
    call vecadd(P(1,k),P(1,k),Pb,n2)
    k = k + 1
  enddo
  
  deallocate(Pa)
  deallocate(Pb)
  return
end subroutine get_sasfcis_density

!======================================================== 
subroutine sasfemake_pmo1_nonsym(Pa,Pb,R,Ca,n)
  implicit none 
  integer i,j,k,n
  real*8  Pa(n2,*),Pb(n2,*),R(NOV,*),Ca(LM2,*)

  k = 1
  do i=1,nex
    do j=i,nex
      call sasfemake_1e_density(Pa(1,k),Pb(1,k),R(1,i),R(1,j),Ca)
      k = k + 1
    enddo
  enddo
  return
end subroutine sasfemake_pmo1_nonsym

!======================================================== 
subroutine sasfemake_1e_density(Pa,Pb,Ri,Rj,Ca)
  implicit none 
  integer n,i,j,k,l,m
  real*8  dot,dot2,dot3,dot4
  real*8  Pa(n2),Pb(n2),Ri(NOV),Rj(NOV),Ca(LM2,*)
  real*8, dimension(:), allocatable :: tmp,tmp2

  allocate(tmp(n2),tmp2(n2))
  call vecinit(Pa,n2,0.d0)
  call vecinit(Pb,n2,0.d0)

  ! s2s2 alpha, beta
  call vecdot(dot,Ri(os2),Rj(os2),NOb)
  Pb(s2s2+1) = Pb(s2s2+1) + Ri(1)*Rj(1) + dot
  call vecdot(dot,Ri(vs2),Rj(vs2),NVa)
  Pa(s2s2+1) = Pa(s2s2+1) - Ri(2)*Rj(2) - dot
  
  ! s1s1 alpha, beta
  call vecdot(dot,Ri(os1),Rj(os1),NOb)
  Pb(s1s1+1) = Pb(s1s1+1) + Ri(2)*Rj(2) + Ri(3)*Rj(3) + dot
  call vecdot(dot,Ri(vs1),Rj(vs1),NVa)
  Pa(s1s1+1) = Pa(s1s1+1) - Ri(1)*Rj(1) - Ri(3)*Rj(3) - dot

  ! s2s1 alpha, beta
  call vecdot(dot ,Ri(os1),Rj(os2),NOb)
  call vecdot(dot2,Rj(os1),Ri(os2),NOb)
  Pb(s2s1+1) = Pb(s2s1+1) + sq2*Ri(1)*Rj(3) + dot2 
  Pb(s1s2+1) = Pb(s1s2+1) + sq2*Rj(1)*Ri(3) + dot
  call vecdot(dot ,Ri(vs1),Rj(vs2),NVa)
  call vecdot(dot2,Rj(vs1),Ri(vs2),NVa)
  Pa(s1s2+1) = Pa(s1s2+1) - sq2*Ri(2)*Rj(3) - dot2
  Pa(s2s1+1) = Pa(s2s1+1) - sq2*Rj(2)*Ri(3) - dot
 
  ! is1 is2 alpha
  call matmult(Ri(vs1),Rj(ov2),tmp ,1,NOb,NVa,NVa,NVa,1,2) 
  call matmult(Rj(vs1),Ri(ov2),tmp2,1,NOb,NVa,NVa,NVa,1,2) 
  do k=1,NOb
    Pa(k+(js1-1)*nbas) = Pa(k+(js1-1)*nbas) &
             - sq2*Ri(1)*Rj(os2+k-1) - Ri(3)*Rj(os1+k-1) - sq2*tmp(k)
    Pa(js1+(k-1)*nbas) = Pa(js1+(k-1)*nbas) &
             - sq2*Rj(1)*Ri(os2+k-1) - Rj(3)*Ri(os1+k-1) - sq2*tmp2(k)
  enddo

  call vecadd3(Pa(i1s2+1),-sq2*Ri(2),Rj(os1),NOb)
  call vecadd3(Pa(i1s2+1),     Ri(3),Rj(os2),NOb)
  call matmult(Ri(vs2),Rj(ov1),tmp,1,NOb,NVa,NVa,NVa,1,2) 
  call vecadd3(Pa(i1s2+1),ssq6,tmp,NOb)
  call matmult(Ri(vs2),Rj(ov2),tmp,1,NOb,NVa,NVa,NVa,1,2) 
  call vecadd3(Pa(i1s2+1),-ssq2,tmp,NOb)

  call matmult(Rj(vs2),Ri(ov1),tmp,1,NOb,NVa,NVa,NVa,1,2)
  call matmult(Rj(vs2),Ri(ov2),tmp2,1,NOb,NVa,NVa,NVa,1,2)      
  do k=1,NOb
    Pa(js2+(k-1)*nbas) = Pa(js2+(k-1)*nbas) &
             - sq2*Rj(2)*Ri(os1+k-1) + Rj(3)*Ri(os2+k-1) &
             + ssq6*tmp(k) - ssq2*tmp2(k)
  enddo

  ! as1 as2 beta
  call vecadd3(Pb(a1s1+1),sq2*Rj(2),Ri(vs2),NVa)
  call vecadd3(Pb(a1s1+1),Rj(3),Ri(vs1),NVa)
  call matmult(Ri(ov1),Rj(os1),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
  call vecadd3(Pb(a1s1+1),-ssq6,tmp,NVa)
  call matmult(Ri(ov2),Rj(os1),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
  call vecadd3(Pb(a1s1+1),ssq2,tmp,NVa)

  call matmult(Rj(ov1),Ri(os1),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
  call matmult(Rj(ov2),Ri(os1),tmp2,NVa,1,NOb,NVa,NOb,NVa,1)
  do k=1,NVa
    Pb(js1+(jv+k-2)*nbas) = Pb(js1+(jv+k-2)*nbas) &
              + sq2*Ri(2)*Rj(vs2+k-1) + Ri(3)*Rj(vs1+k-1) &
              - ssq6*tmp(k) + ssq2*tmp2(k)
  enddo
  
  call vecadd3(Pb(a1s2+1),sq2*Rj(1),Ri(vs1),NVa)
  call vecadd3(Pb(a1s2+1),-Rj(3),Ri(vs2),NVa)
  call matmult(Ri(ov2),Rj(os2),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
  call vecadd3(Pb(a1s2+1),sq2,tmp,NVa)

  call matmult(Rj(ov2),Ri(os2),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
  do k=1,NVa
    Pb(js2+(jv+k-2)*nbas) = Pb(js2+(jv+k-2)*nbas) &
              + sq2*Ri(1)*Rj(vs1+k-1) - Ri(3)*Rj(vs2+k-1) &
              + sq2*tmp(k)
  enddo

  ! ij alpha
  call vecinit(tmp,n2,0.d0)
  call matmult(Rj(os1),Ri(os1),tmp,NOb,NOb,1,NOb,NOb,nbas,3)
  call vecsub(Pa(1),Pa(1),tmp,n2)
  call vecinit(tmp,n2,0.d0)
  call matmult(Rj(os2),Ri(os2),tmp,NOb,NOb,1,NOb,NOb,nbas,3)
  call vecsub(Pa(1),Pa(1),tmp,n2)
  call vecinit(tmp,n2,0.d0)
  call matmult(Rj(ov1),Ri(ov1),tmp,NOb,NOb,NVa,NVa,NVa,nbas,2)
  call vecsub(Pa(1),Pa(1),tmp,n2)
  call vecinit(tmp,n2,0.d0)
  call matmult(Rj(ov2),Ri(ov2),tmp,NOb,NOb,NVa,NVa,NVa,nbas,2)
  call vecsub(Pa(1),Pa(1),tmp,n2)

  ! ab beta
  call vecinit(tmp,n2,0.d0)
  call matmult(Ri(vs1),Rj(vs1),tmp(a1b1+1),NVa,NVa,1,NVa,NVa,nbas,3)
  call vecadd(Pb(1),Pb(1),tmp,n2)
  call vecinit(tmp,n2,0.d0)
  call matmult(Ri(vs2),Rj(vs2),tmp(a1b1+1),NVa,NVa,1,NVa,NVa,nbas,3)
  call vecadd(Pb(1),Pb(1),tmp,n2)
  call vecinit(tmp,n2,0.d0)
  call matmult(Ri(ov1),Rj(ov1),tmp(a1b1+1),NVa,NVa,NOb,NVa,NVa,nbas,3)
  call vecadd(Pb(1),Pb(1),tmp,n2)
  call vecinit(tmp,n2,0.d0)
  call matmult(Ri(ov2),Rj(ov2),tmp(a1b1+1),NVa,NVa,NOb,NVa,NVa,nbas,3)
  call vecadd(Pb(1),Pb(1),tmp,n2)

  ! ai beta
  do l=0,NVa-1
  do k=0,NOb-1
    Pb(jv+l+k*nbas) = Pb(jv+l+k*nbas) + ssq6*Rj(3)*Ri(ov1+l+k*NVa)
    Pb(k+1+(jv+l-1)*nbas) = Pb(k+1+(jv+l-1)*nbas) + ssq6*Ri(3)*Rj(ov1+l+k*NVa)
    Pb(jv+l+k*nbas) = Pb(jv+l+k*nbas) + ssq2*Rj(3)*Ri(ov2+l+k*NVa)
    Pb(k+1+(jv+l-1)*nbas) = Pb(k+1+(jv+l-1)*nbas) + ssq2*Ri(3)*Rj(ov2+l+k*NVa)
  enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(sf_xcis) then
    do l=0,NVa-1
    do k=0,NOb-1
      dot = ssq2*( Rj(1)*Ri(ov3+l+k*NVa) + Rj(2)*Ri(ov4+l+k*NVa) )
      Pa(jv+l+k*nbas) = Pa(jv+l+k*nbas) + dot
      Pb(jv+l+k*nbas) = Pb(jv+l+k*nbas) + dot
      dot = ssq2*( Ri(1)*Rj(ov3+l+k*NVa) + Ri(2)*Rj(ov4+l+k*NVa) )
      Pa(k+1+(jv+l-1)*nbas) = Pa(k+1+(jv+l-1)*nbas) + dot
      Pb(k+1+(jv+l-1)*nbas) = Pb(k+1+(jv+l-1)*nbas) + dot 
    enddo
    enddo
     
    call matmult(Ri(ov3),Rj(os2),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
    call vecsub(Pa(a1s1+1),Pa(a1s1+1),tmp,NVa)
    call matmult(Ri(ov4),Rj(os1),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
    call vecsub(Pa(a1s2+1),Pa(a1s2+1),tmp,NVa)
    call matmult(Rj(ov3),Ri(os2),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
    call matmult(Rj(ov4),Ri(os1),tmp2,NVa,1,NOb,NVa,NOb,NVa,1)
    do k=1,NVa
      Pa(js1+(jv+k-2)*nbas) = Pa(js1+(jv+k-2)*nbas) - tmp(k) 
      Pa(js2+(jv+k-2)*nbas) = Pa(js2+(jv+k-2)*nbas) - tmp2(k)
    enddo
 
    call matmult(Ri(vs1),Rj(ov3),tmp,1,NOb,NVa,NVa,NVa,1,2)
    call vecsub(Pb(i1s2+1),Pb(i1s2+1),tmp,NOb)
    call matmult(Ri(vs2),Rj(ov4),tmp,1,NOb,NVa,NVa,NVa,1,2)
    call vecsub(Pb(i1s1+1),Pb(i1s1+1),tmp,NOb)
 
    call matmult(Rj(vs1),Ri(ov3),tmp ,1,NOb,NVa,NVa,NVa,1,2)
    call matmult(Rj(vs2),Ri(ov4),tmp2,1,NOb,NVa,NVa,NVa,1,2)
    do k=1,NOb
      Pb(js1+(k-1)*nbas) = Pb(js1+(k-1)*nbas) - tmp2(k)
      Pb(js2+(k-1)*nbas) = Pb(js2+(k-1)*nbas) - tmp(k)
    enddo
 
    call vecdot(dot ,Ri(ov1),Rj(ov3),NVa*NOb)
    call vecdot(dot2,Ri(ov2),Rj(ov3),NVa*NOb)
    call vecdot(dot3,Rj(ov1),Ri(ov4),NVa*NOb)
    call vecdot(dot4,Rj(ov2),Ri(ov4),NVa*NOb)
    Pa(s1s2+1) = Pa(s1s2+1) + ssq6*dot + ssq2*dot2 - ssq6*dot3 - ssq2*dot4      
    call vecdot(dot ,Rj(ov1),Ri(ov3),NVa*NOb)
    call vecdot(dot2,Rj(ov2),Ri(ov3),NVa*NOb)
    call vecdot(dot3,Ri(ov1),Rj(ov4),NVa*NOb)
    call vecdot(dot4,Ri(ov2),Rj(ov4),NVa*NOb)
    Pa(s2s1+1) = Pa(s2s1+1) + ssq6*dot + ssq2*dot2 - ssq6*dot3 - ssq2*dot4      
    call vecdot(dot ,Rj(ov3),Ri(ov3),NVa*NOb)
    call vecdot(dot2,Rj(ov4),Ri(ov4),NVa*NOb)
    Pa(s1s1+1) = Pa(s1s1+1) - dot
    Pb(s2s2+1) = Pb(s2s2+1) + dot
    Pb(s1s1+1) = Pb(s1s1+1) + dot2
    Pa(s2s2+1) = Pa(s2s2+1) - dot2
 
    call vecinit(tmp,n2,0.d0)
    call matmult(Ri(ov3),Rj(ov3),tmp(a1b1+1),NVa,NVa,NOb,NVa,NVa,nbas,3)
    call vecadd(Pa(1),Pa(1),tmp,n2)
    call vecinit(tmp,n2,0.d0)
    call matmult(Ri(ov4),Rj(ov4),tmp(a1b1+1),NVa,NVa,NOb,NVa,NVa,nbas,3)
    call vecadd(Pa(1),Pa(1),tmp,n2)
 
    call vecinit(tmp,n2,0.d0)
    call matmult(Rj(ov3),Ri(ov3),tmp,NOb,NOb,NVa,NVa,NVa,nbas,2)
    call vecsub(Pa(1),Pa(1),tmp,n2)
    call vecinit(tmp,n2,0.d0)
    call matmult(Rj(ov4),Ri(ov4),tmp,NOb,NOb,NVa,NVa,NVa,nbas,2)
    call vecsub(Pa(1),Pa(1),tmp,n2)
  
  endif !sf_xcis

  deallocate(tmp)
  deallocate(tmp2)
  return
end subroutine sasfemake_1e_density

!======================================================== 
subroutine sasfemake_pmo2(P11,P22,P12,P11J,P22J,P12J,R,Ca,n,igs)
  implicit none 
  integer n,i,j,k,m,igs
  real*8  dot,dot2
  real*8  P11(n2,*),P22(n2,*),P12(n2,*),P11J(n2,*),P22J(n2,*),P12J(n2,*),R(NOV,*),Ca(LM2,*)
  real*8, dimension(:), allocatable :: tmp

  allocate(tmp(n2))
  call vecinit(P11,n2*n,0.d0)
  call vecinit(P22,n2*n,0.d0)
  if(sf_xcis) then
  call vecinit(P12 ,n2*n,0.d0)
  call vecinit(P11J,n2*n,0.d0)
  call vecinit(P22J,n2*n,0.d0)
  call vecinit(P12J,n2*n,0.d0)
  endif
  m = 1
  do i=1,nex
    do j=i,nex
      ! ab 
      call vecinit(tmp,n2,0.d0)
      call matmult(R(ov1,i),R(ov1,j),tmp(a1b1+1),NVa,NVa,NOb,NVa,NVa,nbas,3)
      call vecsub(P22(1,m),P22(1,m),tmp,n2)
      call vecadd3(P11(1,m),f12,tmp,n2)
      call vecinit(tmp,n2,0.d0)
      call matmult(R(ov2,i),R(ov2,j),tmp(a1b1+1),NVa,NVa,NOb,NVa,NVa,nbas,3)
      call vecadd(P22(1,m),P22(1,m),tmp,n2)
      call vecadd3(P11(1,m),-f12,tmp,n2)
      call vecinit(tmp,n2,0.d0)
      call matmult(R(ov1,i),R(ov2,j),tmp(a1b1+1),NVa,NVa,NOb,NVa,NVa,nbas,3)
      call vecadd3(P11(1,m),-ssq3,tmp,n2)
      call vecinit(tmp,n2,0.d0)
      !call matmult(R(ov2,i),R(ov1,j),tmp(a1b1+1),NVa,NVa,NOb,NVa,NVa,nbas,3)
      call matmult(R(ov1,j),R(ov2,i),tmp(a1b1+1),NVa,NVa,NOb,NVa,NVa,nbas,3)
      call vecadd3(P11(1,m),-ssq3,tmp,n2)

      ! ij
      call vecinit(tmp,n2,0.d0)
      call matmult(R(ov1,i),R(ov1,j),tmp,NOb,NOb,NVa,NVa,NVa,nbas,2)
      call vecsub(P11(1,m),P11(1,m),tmp,n2)
      call vecadd3(P22(1,m),f12,tmp,n2)
      call vecinit(tmp,n2,0.d0)
      call matmult(R(ov2,i),R(ov2,j),tmp,NOb,NOb,NVa,NVa,NVa,nbas,2)
      call vecadd(P11(1,m),P11(1,m),tmp,n2)
      call vecadd3(P22(1,m),-f12,tmp,n2)
      call vecinit(tmp,n2,0.d0)
      call matmult(R(ov1,i),R(ov2,j),tmp,NOb,NOb,NVa,NVa,NVa,nbas,2)
      call vecadd3(P22(1,m),-ssq3,tmp,n2)
      call vecinit(tmp,n2,0.d0)
      !call matmult(R(ov2,i),R(ov1,j),tmp,NOb,NOb,NVa,NVa,NVa,nbas,2)
      call matmult(R(ov1,j),R(ov2,i),tmp,NOb,NOb,NVa,NVa,NVa,nbas,2)
      call vecadd3(P22(1,m),-ssq3,tmp,n2)

      ! is1 is2
      call matmult(R(vs1,i),R(ov1,j),tmp,1,NOb,NVa,NOb,NVa,1,2)
      !do k=0,NOb-1
      !  P22(js1+k*nbas,m) = P22(js1+k*nbas,m) - ssq6*tmp(k+1)
      !enddo
      call vecadd3(P22(i1s1+1,m),-ssq6,tmp,NOb)
      call matmult(R(vs1,j),R(ov1,i),tmp,1,NOb,NVa,NOb,NVa,1,2)
      !do k=0,NOb-1
      !  P22(js1+k*nbas,m) = P22(js1+k*nbas,m) - ssq6*tmp(k+1)
      !enddo
      call vecadd3(P22(i1s1+1,m),-ssq6,tmp,NOb)
      call matmult(R(vs1,i),R(ov2,j),tmp,1,NOb,NVa,NOb,NVa,1,2)
      !do k=0,NOb-1
      !  P22(js1+k*nbas,m) = P22(js1+k*nbas,m) - ssq2*tmp(k+1)
      !enddo
      call vecadd3(P22(i1s1+1,m),-ssq2,tmp,NOb)
      call matmult(R(vs1,j),R(ov2,i),tmp,1,NOb,NVa,NOb,NVa,1,2)
      !do k=0,NOb-1
      !  P22(js1+k*nbas,m) = P22(js1+k*nbas,m) - ssq2*tmp(k+1)
      !enddo
      call vecadd3(P22(i1s1+1,m),-ssq2,tmp,NOb)

      call matmult(R(vs2,i),R(ov1,j),tmp,1,NOb,NVa,NOb,NVa,1,2)
      !do k=0,NOb-1
      !  P11(js2+k*nbas,m) = P11(js2+k*nbas,m) + ssq6*tmp(k+1)
      !enddo
      call vecadd3(P11(i1s2+1,m),ssq6,tmp,NOb)
      call matmult(R(vs2,j),R(ov1,i),tmp,1,NOb,NVa,NOb,NVa,1,2)
      !do k=0,NOb-1
      !  P11(js2+k*nbas,m) = P11(js2+k*nbas,m) + ssq6*tmp(k+1)
      !enddo
      call vecadd3(P11(i1s2+1,m),ssq6,tmp,NOb)
      call matmult(R(vs2,i),R(ov2,j),tmp,1,NOb,NVa,NOb,NVa,1,2)
      !do k=0,NOb-1
      !  P11(js2+k*nbas,m) = P11(js2+k*nbas,m) + ssq2*tmp(k+1)
      !enddo
      call vecadd3(P11(i1s2+1,m),ssq2,tmp,NOb)
      call matmult(R(vs2,j),R(ov2,i),tmp,1,NOb,NVa,NOb,NVa,1,2)
      !do k=0,NOb-1
      !  P11(js2+k*nbas,m) = P11(js2+k*nbas,m) + ssq2*tmp(k+1)
      !enddo
      call vecadd3(P11(i1s2+1,m),ssq2,tmp,NOb)

      ! as1 as2
      call matmult(R(ov1,j),R(os1,i),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecadd3(P22(a1s1+1,m),ssq6,tmp,NVa)
      call matmult(R(ov1,i),R(os1,j),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecadd3(P22(a1s1+1,m),ssq6,tmp,NVa)
      call matmult(R(ov2,j),R(os1,i),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecadd3(P22(a1s1+1,m),ssq2,tmp,NVa)
      call matmult(R(ov2,i),R(os1,j),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecadd3(P22(a1s1+1,m),ssq2,tmp,NVa)

      call matmult(R(ov1,j),R(os2,i),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecadd3(P11(a1s2+1,m),-ssq6,tmp,NVa)
      call matmult(R(ov1,i),R(os2,j),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecadd3(P11(a1s2+1,m),-ssq6,tmp,NVa)
      call matmult(R(ov2,j),R(os2,i),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecadd3(P11(a1s2+1,m),-ssq2,tmp,NVa)
      call matmult(R(ov2,i),R(os2,j),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecadd3(P11(a1s2+1,m),-ssq2,tmp,NVa)

      call vecdot(dot,R(ov1,i),R(ov1,j),NVa*NOb)
      P22(s1s1+1,m) = P22(s1s1+1,m) + f32*dot
      call vecdot(dot,R(ov2,i),R(ov2,j),NVa*NOb)
      P22(s1s1+1,m) = P22(s1s1+1,m) + f12*dot
      call vecdot(dot,R(ov1,i),R(ov2,j),NVa*NOb)
      P22(s1s1+1,m) = P22(s1s1+1,m) + ssq3*dot
      call vecdot(dot,R(ov1,j),R(ov2,i),NVa*NOb)
      P22(s1s1+1,m) = P22(s1s1+1,m) + ssq3*dot
      if((igs.eq.1).and.(i.eq.j)) then
        if(tref) then
          P22(s1s1+1,m) = P22(s1s1+1,m) - f14
          P22(s2s2+1,m) = P22(s2s2+1,m) - f14
          P11(s1s1+1,m) = P11(s1s1+1,m) - f14
          P11(s2s2+1,m) = P11(s2s2+1,m) - f14 
        !else if(sref) then
        !  P11(s1s1+1,m) = P11(s1s1+1,m) + f12
        !  P22(s2s2+1,m) = P22(s2s2+1,m) + f12
        endif
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(sf_xcis) then
      call matmult(R(ov3,j),R(os1,i),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecsub(P12(a1s1+1,m),P12(a1s1+1,m),tmp,NVa)
      call matmult(R(ov3,i),R(os1,j),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecsub(P12(a1s1+1,m),P12(a1s1+1,m),tmp,NVa)
      call matmult(R(ov3,j),R(os2,i),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecsub(P12(a1s2+1,m),P12(a1s2+1,m),tmp,NVa)
      call matmult(R(ov3,i),R(os2,j),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecsub(P12(a1s2+1,m),P12(a1s2+1,m),tmp,NVa)

      call matmult(R(ov4,j),R(os1,i),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecsub(P11J(a1s2+1,m),P11J(a1s2+1,m),tmp,NVa)
      call matmult(R(ov4,i),R(os1,j),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecsub(P11J(a1s2+1,m),P11J(a1s2+1,m),tmp,NVa)
      call matmult(R(ov4,j),R(os2,i),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecsub(P12J(a1s2+1,m),P12J(a1s2+1,m),tmp,NVa)
      call matmult(R(ov4,i),R(os2,j),tmp,NVa,1,NOb,NVa,NOb,NVa,1)
      call vecsub(P12J(a1s2+1,m),P12J(a1s2+1,m),tmp,NVa)

      call matmult(R(vs1,i),R(ov3,j),tmp,1,NOb,NVa,NOb,NVa,1,2)
      call vecadd(P11J(i1s2+1,m),P11J(i1s2+1,m),tmp,NOb)
      call matmult(R(vs1,j),R(ov3,i),tmp,1,NOb,NVa,NOb,NVa,1,2)
      call vecadd(P11J(i1s2+1,m),P11J(i1s2+1,m),tmp,NOb)
      call matmult(R(vs2,i),R(ov3,j),tmp,1,NOb,NVa,NOb,NVa,1,2)
      call vecadd(P12J(i1s2+1,m),P12J(i1s2+1,m),tmp,NOb)
      call matmult(R(vs2,j),R(ov3,i),tmp,1,NOb,NVa,NOb,NVa,1,2)
      call vecadd(P12J(i1s2+1,m),P12J(i1s2+1,m),tmp,NOb)

      call matmult(R(vs1,i),R(ov4,j),tmp,1,NOb,NVa,NOb,NVa,1,2)
      call vecadd(P12J(i1s1+1,m),P12J(i1s1+1,m),tmp,NOb)
      call matmult(R(vs1,j),R(ov4,i),tmp,1,NOb,NVa,NOb,NVa,1,2)
      call vecadd(P12J(i1s1+1,m),P12J(i1s1+1,m),tmp,NOb)
      call matmult(R(vs2,i),R(ov4,j),tmp,1,NOb,NVa,NOb,NVa,1,2)
      call vecadd(P22J(i1s1+1,m),P22J(i1s1+1,m),tmp,NOb)
      call matmult(R(vs2,j),R(ov4,i),tmp,1,NOb,NVa,NOb,NVa,1,2)
      call vecadd(P22J(i1s1+1,m),P22J(i1s1+1,m),tmp,NOb)

      call vecdot(dot ,R(ov1,i),R(ov3,j),NVa*NOb)
      call vecdot(dot2,R(ov1,j),R(ov3,i),NVa*NOb)
      P22(s1s2+1,m) = P22(s1s2+1,m) + ssq6*dot + ssq6*dot2 ! AC
      call vecdot(dot ,R(ov2,i),R(ov3,j),NVa*NOb)
      call vecdot(dot2,R(ov2,j),R(ov3,i),NVa*NOb)
      P22(s1s2+1,m) = P22(s1s2+1,m) + ssq2*dot + ssq2*dot2 ! BC
      call vecdot(dot,R(ov3,i),R(ov3,j),NVa*NOb)
      P12(s1s2+1,m) = P12(s1s2+1,m) - dot                  ! CC
      call vecdot(dot ,R(ov4,i),R(ov3,j),NVa*NOb)
      call vecdot(dot2,R(ov4,j),R(ov3,i),NVa*NOb)
      P22(s1s1+1,m) = P22(s1s1+1,m) - dot - dot2           ! DC
      call vecdot(dot ,R(ov1,i),R(ov4,j),NVa*NOb)
      call vecdot(dot2,R(ov1,j),R(ov4,i),NVa*NOb)
      P11(s1s2+1,m) = P11(s1s2+1,m) - ssq6*dot - ssq6*dot2 ! AD
      call vecdot(dot ,R(ov2,i),R(ov4,j),NVa*NOb)
      call vecdot(dot2,R(ov2,j),R(ov4,i),NVa*NOb)
      P11(s1s2+1,m) = P11(s1s2+1,m) - ssq2*dot - ssq2*dot2 ! BD
      call vecdot(dot ,R(ov4,i),R(ov4,j),NVa*NOb)
      P12(s1s2+1,m) = P12(s1s2+1,m) - dot                  ! DD
      
      ! ij
      call matmult(R(ov1,i),R(ov3,j),tmp,NOb,NOb,NVa,NVa,NVa,NOb,2)
      call matadd2(P12J(1,m),-ssq6,tmp,nbas,NOb,NOb)
      call matmult(R(ov1,j),R(ov3,i),tmp,NOb,NOb,NVa,NVa,NVa,NOb,2)
      call matadd2(P12J(1,m),-ssq6,tmp,nbas,NOb,NOb)
      call matmult(R(ov2,i),R(ov3,j),tmp,NOb,NOb,NVa,NVa,NVa,NOb,2)
      call matadd2(P12J(1,m),-ssq2,tmp,nbas,NOb,NOb)
      call matadd2(P12(1,m),sq2,tmp,nbas,NOb,NOb)
      call matmult(R(ov2,j),R(ov3,i),tmp,NOb,NOb,NVa,NVa,NVa,NOb,2)
      call matadd2(P12J(1,m),-ssq2,tmp,nbas,NOb,NOb)
      call matadd2(P12(1,m),sq2,tmp,nbas,NOb,NOb)
      call matmult(R(ov3,i),R(ov3,j),tmp,NOb,NOb,NVa,NVa,NVa,NOb,2)
      call matadd(P11J(1,m),P11J(1,m),tmp,nbas,NOb,NOb)
      !write(6,*) "P11"
      !call matprnt(P11(1,m),nbas,nbas,6)
      !call matprnt(tmp,NOb,NOb,6)
      call matsub(P11(1,m),P11(1,m),tmp,nbas,NOb,NOb)
      !call matprnt(P11(1,m),nbas,nbas,6)
      call matsub(P22J(1,m),P22J(1,m),tmp,nbas,NOb,NOb)

      call matmult(R(ov4,i),R(ov1,j),tmp,NOb,NOb,NVa,NVa,NVa,NOb,2)
      call matadd2(P12J(1,m),ssq6,tmp,nbas,NOb,NOb)
      call matadd2(P12(1,m),-ssq6,tmp,nbas,NOb,NOb)
      call matmult(R(ov4,j),R(ov1,i),tmp,NOb,NOb,NVa,NVa,NVa,NOb,2)
      call matadd2(P12J(1,m),ssq6,tmp,nbas,NOb,NOb)
      call matadd2(P12(1,m),-ssq6,tmp,nbas,NOb,NOb)
      call matmult(R(ov4,i),R(ov2,j),tmp,NOb,NOb,NVa,NVa,NVa,NOb,2)
      call matadd2(P12J(1,m),ssq2,tmp,nbas,NOb,NOb)
      call matadd2(P12(1,m),ssq2,tmp,nbas,NOb,NOb)
      call matmult(R(ov4,j),R(ov2,i),tmp,NOb,NOb,NVa,NVa,NVa,NOb,2)
      call matadd2(P12J(1,m),ssq2,tmp,nbas,NOb,NOb)
      call matadd2(P12(1,m),ssq2,tmp,nbas,NOb,NOb)
      call matmult(R(ov4,i),R(ov4,j),tmp,NOb,NOb,NVa,NVa,NVa,NOb,2)
      call matsub(P11J(1,m),P11J(1,m),tmp,nbas,NOb,NOb)
      call matadd(P22J(1,m),P22J(1,m),tmp,nbas,NOb,NOb)
      call matsub(P22(1,m),P22(1,m),tmp,nbas,NOb,NOb)

      ! ab
      call matmult(R(ov1,i),R(ov3,j),tmp,NVa,NVa,NOb,NVa,NVa,NVa,3)
      call matadd2(P12J(a1b1+1,m),ssq6,tmp,nbas,NVa,NVa)
      call matmult(R(ov1,j),R(ov3,i),tmp,NVa,NVa,NOb,NVa,NVa,NVa,3)
      call matadd2(P12J(a1b1+1,m),ssq6,tmp,nbas,NVa,NVa)
      call matmult(R(ov3,i),R(ov2,j),tmp,NVa,NVa,NOb,NVa,NVa,NVa,3)
      call matadd2(P12J(a1b1+1,m),ssq2,tmp,nbas,NVa,NVa)
      call matadd2(P12(a1b1+1,m),-sq2,tmp,nbas,NVa,NVa)
      call matmult(R(ov3,j),R(ov2,i),tmp,NVa,NVa,NOb,NVa,NVa,NVa,3)
      call matadd2(P12J(a1b1+1,m),ssq2,tmp,nbas,NVa,NVa)
      call matadd2(P12(a1b1+1,m),-sq2,tmp,nbas,NVa,NVa)
      call matmult(R(ov3,i),R(ov3,j),tmp,NVa,NVa,NOb,NVa,NVa,NVa,3)
      call matsub(P11J(a1b1+1,m),P11J(a1b1+1,m),tmp,nbas,NVa,NVa)
      call matadd(P11(a1b1+1,m),P11(a1b1+1,m),tmp,nbas,NVa,NVa)
      call matadd(P22J(a1b1+1,m),P22J(a1b1+1,m),tmp,nbas,NVa,NVa)


      call matmult(R(ov1,i),R(ov4,j),tmp,NVa,NVa,NOb,NVa,NVa,NVa,3)
      call matadd2(P12J(a1b1+1,m),-ssq6,tmp,nbas,NVa,NVa)
      call matadd2(P12(a1b1+1,m),ssq6,tmp,nbas,NVa,NVa)
      call matmult(R(ov1,j),R(ov4,i),tmp,NVa,NVa,NOb,NVa,NVa,NVa,3)
      call matadd2(P12J(a1b1+1,m),-ssq6,tmp,nbas,NVa,NVa)
      call matadd2(P12(a1b1+1,m),ssq6,tmp,nbas,NVa,NVa)
      call matmult(R(ov2,i),R(ov4,j),tmp,NVa,NVa,NOb,NVa,NVa,NVa,3)
      call matadd2(P12J(a1b1+1,m),-ssq2,tmp,nbas,NVa,NVa)
      call matadd2(P12(a1b1+1,m),-ssq2,tmp,nbas,NVa,NVa)
      call matmult(R(ov2,j),R(ov4,i),tmp,NVa,NVa,NOb,NVa,NVa,NVa,3)
      call matadd2(P12J(a1b1+1,m),-ssq2,tmp,nbas,NVa,NVa)
      call matadd2(P12(a1b1+1,m),-ssq2,tmp,nbas,NVa,NVa)
      call matmult(R(ov4,i),R(ov4,j),tmp,NVa,NVa,NOb,NVa,NVa,NVa,3)
      call matadd(P11J(a1b1+1,m),P11J(a1b1+1,m),tmp,nbas,NVa,NVa)
      call matsub(P22J(a1b1+1,m),P22J(a1b1+1,m),tmp,nbas,NVa,NVa)
      call matadd(P22(a1b1+1,m),P22(a1b1+1,m),tmp,nbas,NVa,NVa)
      
      endif

      m = m + 1
    enddo
  enddo
  
  deallocate(tmp)
  return
end subroutine sasfemake_pmo2

subroutine sasfemake_a2e2(Za,Ja,Ka,Jss,Kss,R,Ca,nr)
  implicit none
  integer nr
  real*8  Za(NOV,*),Ja(n2*npr,*),Ka(n2*npr,*),Jss(n2,4),Kss(n2,4),R(NOV,*),Ca(LM2,*)
  real*8, dimension(:), allocatable :: tmp,tmp2

  allocate(tmp(n2*4),tmp2(n2*4))
  call veccopy(tmp,Kss,n2*4)
  call veccopy(tmp2,Jss,n2*4)
  call ao2mo(tmp,Ca,4)
  call ao2mo(tmp2,Ca,4)
  call sasfemake_a2e(Za,Ja,Ka,tmp2,tmp,R,Ca,nr)
  deallocate(tmp)
  deallocate(tmp2)
  return
end subroutine sasfemake_a2e2

subroutine sasfemake_a2e(Za,Ja,Ka,Jss,Kss,R,Ca,nr)
  implicit none
  integer nr,i,m,ioff1,ioff2
  real*8  Za(NOV,*),Ja(n2*npr,*),Ka(n2*npr,*),Jss(n2,4),Kss(n2,4),R(NOV,*),Ca(LM2,*)
  real*8, dimension(:), allocatable :: tmp,tmp2

  allocate(tmp(n2),tmp2(n2))
  call vecinit(Za,NOV*nr,0.d0)
  do i=1,nr
    call ao2mo2(Ja(1,i),Ca,npr)
    call ao2mo2(Ka(1,i),Ca,npr)

    if(sf_xcis) then
      call vecadd3(Ka(jJov3,i),two,Ja(jJov3,i),NVb*NOa)
      call vecadd3(Ka(jJov4,i),two,Ja(jJov4,i),NVb*NOa)
    endif

    ! s1->s2
    Za(1,i) = -Jss(s2s2+1,1)*R(1,i)-Jss(s1s2+1,2)*R(2,i)-sq2*Jss(s1s2+1,1)*R(3,i)   &
              +sq2*Ka(jJos1+s21,i)+sq2*Ka(jJos2+s21,i)                              &
              +sq2*Ka(jJvs1+s21,i)+sq2*Ka(jJvs2+s21,i)                              &
              +sq3*Ja(jJov1+s21,i)+Ja(jJov2+s21,i)+two*Ka(jJov2+s21,i)

    ! s2->s1
    Za(2,i) = -Jss(s1s2+1,2)*R(1,i)-Jss(s2s2+1,1)*R(2,i)-sq2*Jss(s1s2+1,1)*R(3,i)   &
              +sq2*Ka(jJos1+s12,i)+sq2*Ka(jJos2+s12,i)                              &
              +sq2*Ka(jJvs1+s12,i)+sq2*Ka(jJvs2+s12,i)                              &
              -sq3*Ja(jJov1+s12,i)-sq3*Ka(jJov1+s12,i)                              &
              -Ja(jJov2+s12,i)+Ka(jJov2+s12,i)

    ! s1->s1
    Za(3,i) = -sq2*Jss(s1s2+1,1)*R(1,i)-sq2*Jss(s1s2+1,1)*R(2,i)                    &
              -(Jss(s1s1+1,1)-Jss(s1s2+1,2))*R(3,i)                                 &
              +Ka(jJos1+s11,i)-Ka(jJos1+s22,i)+Ka(jJos2+s11,i)-Ka(jJos2+s22,i)      &
              +Ka(jJvs1+s11,i)-Ka(jJvs1+s22,i)+Ka(jJvs2+s11,i)-Ka(jJvs2+s22,i)      &
              +ssq6*Ka(jJov1+s22,i)+sq2*Ka(jJov2+s11,i)-ssq2*Ka(jJov2+s22,i)

    ! i ->s1
    ioff1 = 0
    ioff2 = 1
    do m=0,NOb-1
      Za(os1+m,i) =  sq2*Kss(i1s1+m+1,2)*R(1,i)+sq2*Kss(NOb+1+m*NBas,2)*R(2,i)      &
                    +(Kss(NOb+1+m*NBas,1)-Kss(NOb+1+m*NBas,3))*R(3,i)               &
                    +Ka(jJos1+ioff1,i)+Ja(jJos12+ioff2,i)                           &
                    +Ka(jJos2+ioff1,i)-Ja(jJos2+ioff1,i)                            &
                    +Ka(jJvs1+ioff1,i)                                              &
                    +Ja(jJvs2+ioff1,i)+two*Ka(jJvs2+ioff1,i)                        &
                    -ssq6*Ja(jJov1+ioff1,i)-ssq6*Ka(jJov1+ioff1,i)                  &
                    -ssq2*Ja(jJov2+ioff1,i)+ssq2*Ka(jJov2+ioff1,i)
      ioff1 = ioff1 + NVb
      ioff2 = ioff2 + NVb
    enddo

    ! i ->s2
    ioff1 = 0
    ioff2 = 1
    do m=0,NOb-1
      Za(os2+m,i) =  sq2*Kss(i1s2+m+1,2)*R(1,i)+sq2*Kss(js2+m*NBas,2)*R(2,i)        &
                    +(Kss(js2+m*NBas,1)-Kss(js2+m*NBas,3))*R(3,i)                   &
                    -Ja(jJos1+ioff2,i)+Ka(jJos1+ioff2,i)                            &
                    +Ja(jJos21+ioff1,i)+Ka(jJos2+ioff2,i)                           &
                    +two*Ka(jJvs1+ioff2,i)+Ja(jJvs1+ioff2,i)                        &
                    +Ka(jJvs2+ioff2,i)                                              &
                    +ssq6*Ja(jJov1+ioff2,i)                                         &
                    +ssq2*Ja(jJov2+ioff2,i)+sq2*Ka(jJov2+ioff2,i)
      ioff1 = ioff1 + NVb
      ioff2 = ioff2 + NVb
    enddo

    ! s1-> a
    do m=0,NVa-1
      Za(vs1+m,i) =  sq2*Kss(js1+(NOa+m)*NBas,2)*R(1,i)+sq2*Kss(a1s1+m+1,2)*R(2,i)  &
                    +(Kss(a1s1+m+1,1)-Kss(a1s1+m+1,3))*R(3,i)                       &
                    +Ka(jJos1+as1+m,i)+two*Ka(jJos2+as1+m,i)+Ja(jJos2+as1+m,i)      &
                    +Ka(jJvs1+as1+m,i)+Ja(jJvs12+as2+m,i)                           &
                    +Ka(jJvs2+as1+m,i)-Ja(jJvs2+as1+m,i)                            &
                    +ssq6*Ja(jJov1+as1+m,i)+ssq2*Ja(jJov2+as1+m,i)                  &
                    +sq2*Ka(jJov2+as1+m,i)                                          
    enddo
    
    ! s2-> a
    do m=0,NVa-1
      Za(vs2+m,i) =  sq2*Kss(js2+(NOa+m)*NBas,2)*R(1,i)+sq2*Kss(a1s2+m+1,2)*R(2,i)  &
                    +(Kss(a1s2+m+1,1)-Kss(a1s2+m+1,3))*R(3,i)                       &
                    +two*Ka(jJos1+as2+m,i)+Ja(jJos1+as2+m,i)+Ka(jJos2+as2+m,i)      &
                    +Ka(jJvs1+as2+m,i)-Ja(jJvs1+as2+m,i)                            &
                    +Ka(jJvs2+as2+m,i)+Ja(jJvs21+as1+m,i)                           &
                    -ssq6*Ja(jJov1+as2+m,i)-ssq6*Ka(jJov1+as2+m,i)                  &
                    +ssq2*Ka(jJov2+as2+m,i)-ssq2*Ja(jJov2+as2+m,i)                  
    enddo

    ! i -> a A
    ioff1 = ov1
    do m=0,NOb-1
      call vecadd3(Za(ioff1,i),sq3*R(1,i),Jss(jv+m*NBas,2),NVa)
      call vecadd(tmp,Jss(jv+m*NBas,2),Kss(jv+m*NBas,2),NVa)
      call vecadd3(Za(ioff1,i),-sq3*R(2,i),tmp,NVa)
      call vecadd3(Za(ioff1,i),ssq6*R(3,i),Kss(jv+m*NBas,3),NVa)
      call vecadd(tmp,Ja(jJos1+2+m*NVb,i),Ka(jJos1+2+m*NVb,i),NVa)
      call vecadd3(Za(ioff1,i),-ssq6,tmp,NVa)
      call vecadd3(Za(ioff1,i),ssq6,Ja(jJos2+2+m*NVb,i),NVa)
      call vecadd3(Za(ioff1,i),ssq6,Ja(jJvs1+2+m*NVb,i),NVa)
      call vecadd(tmp,Ja(jJvs2+2+m*NVb,i),Ka(jJvs2+2+m*NVb,i),NVa)
      call vecadd3(Za(ioff1,i),-ssq6,tmp,NVa)
      call vecadd2(tmp,f32,Ja(jJov1+2+m*NVb,i),Ka(jJov1+2+m*NVb,i),NVa)
      call vecadd(Za(ioff1,i),Za(ioff1,i),tmp,NVa)
      call vecadd3(Za(ioff1,i),ssq3,Ja(jJov2+2+m*NVb,i),NVa)
      ioff1 = ioff1 + NVa
    enddo    

    ! i -> a B 
    ioff2 = ov2
    do m=0,NOb-1
      call vecadd2(tmp,two,Kss(jv+m*NBas,4),Jss(jv+m*NBas,2),NVa)
      call vecadd3(Za(ioff2,i),R(1,i),tmp,NVa)
      call vecsub(tmp,Jss(jv+m*NBas,2),Kss(jv+m*NBas,2),NVa)
      call vecadd3(Za(ioff2,i),-R(2,i),tmp,NVa)
      call vecadd2(tmp,-two,Kss(jv+m*NBas,1),Kss(jv+m*NBas,3),NVa)
      call vecadd3(Za(ioff2,i),-ssq2*R(3,i),tmp,NVa)
      call vecsub(tmp,Ja(jJos1+2+m*NVb,i),Ka(jJos1+2+m*NVb,i),NVa)
      call vecadd3(Za(ioff2,i),-ssq2,tmp,NVa)
      call vecadd2(tmp,two,Ka(jJos2+2+m*NVb,i),Ja(jJos2+2+m*NVb,i),NVa)
      call vecadd3(Za(ioff2,i),ssq2,tmp,NVa)
      call vecadd2(tmp,two,Ka(jJvs1+2+m*NVb,i),Ja(jJvs1+2+m*NVb,i),NVa)
      call vecadd3(Za(ioff2,i),ssq2,tmp,NVa)
      call vecsub(tmp,Ja(jJvs2+2+m*NVb,i),Ka(jJvs2+2+m*NVb,i),NVa)
      call vecadd3(Za(ioff2,i),-ssq2,tmp,NVa)
      call vecadd3(Za(ioff2,i),ssq3,Ja(jJov1+2+m*NVb,i),NVa)
      call vecadd2(tmp,f12,Ja(jJov2+2+m*NVb,i),Ka(jJov2+2+m*NVb,i),NVa)
      call vecadd(Za(ioff2,i),Za(ioff2,i),tmp,NVa)
      ioff2 = ioff2 + NVa
    enddo

    if(sf_xcis) then
      Za(1,i) = Za(1,i) + ssq2*(Ka(jJov3+s22,i)-Ka(jJov3+s11,i))
 
      Za(2,i) = Za(2,i) + ssq2*(Ka(jJov4+s11,i)-Ka(jJov4+s22,i))
 
      Za(3,i) = Za(3,i) + Ka(jJov3+s12,i) - Ka(jJov4+s21,i)
 
      ioff2 = 1
      do m=0,NOb-1
        Za(os1+m,i) = Za(os1+m,i) - Ka(jJov4+ioff2,i)
        ioff2 = ioff2 + NVb
      enddo
 
      ioff1 = 0
      do m=0,NOb-1
        Za(os2+m,i) = Za(os2+m,i) - Ka(jJov3+ioff1,i)
        ioff1 = ioff1 + NVb
      enddo
 
      do m=0,NVa-1
        Za(vs1+m,i) = Za(vs1+m,i) + Ka(jJov3+as2+m,i)
      enddo
 
      do m=0,NVa-1
        Za(vs2+m,i) = Za(vs2+m,i) + Ka(jJov4+as1+m,i)
      enddo
 
      ioff2 = ov3
      do m=0,NOb-1
        call vecadd2(tmp,two,Jss(jv+m*nbas,3),Kss(jv+m*nbas,3),NVa)
        call vecadd3(Za(ioff2,i),ssq2*R(1,i),tmp,NVa)
        call vecadd2(tmp,two,Jss(jv+m*nbas,1),Kss(jv+m*nbas,1),NVa)
        call vecadd3(Za(ioff2,i),-ssq2*R(1,i),tmp,NVa)
        call vecadd2(tmp,two,Jss(jv+m*nbas,2),Kss(jv+m*nbas,2),NVa)
        call vecadd3(Za(ioff2,i),R(3,i),tmp,NVa)
        call vecadd3(Za(ioff2,i),-two,Ja(jJos21+2+m*NVb,i),NVa)
        call vecsub(Za(ioff2,i),Za(ioff2,i),Ka(jJos21+2+m*NVb,i),NVa)
        call vecadd3(Za(ioff2,i),two,Ja(jJvs12+2+m*NVb,i),NVa)
        call vecadd(Za(ioff2,i),Za(ioff2,i),Ka(jJvs12+2+m*NVb,i),NVa)
        call vecadd(Za(ioff2,i),Za(ioff2,i),Ka(jJov3+2+m*NVb,i),NVa)
        ioff2 = ioff2 + NVa
      enddo
 
      ioff2 = ov4
      do m=0,NOb-1
        call vecadd2(tmp,two,Jss(jv+m*nbas,3),Kss(jv+m*nbas,3),NVa)
        call vecadd3(Za(ioff2,i),-ssq2*R(2,i),tmp,NVa)
        call vecadd2(tmp,two,Jss(jv+m*nbas,1),Kss(jv+m*nbas,1),NVa)
        call vecadd3(Za(ioff2,i),ssq2*R(2,i),tmp,NVa)
        call mattrans(tmp2,Kss(1,2),nbas,nbas)
        call vecadd2(tmp,two,Jss(jv+m*nbas,2),tmp2(jv+m*nbas),NVa)
        call vecadd3(Za(ioff2,i),-R(3,i),tmp,NVa)
        call vecadd3(Za(ioff2,i),-two,Ja(jJos12+2+m*NVb,i),NVa)
        call vecsub(Za(ioff2,i),Za(ioff2,i),Ka(jJos12+2+m*NVb,i),NVa)
        call vecadd3(Za(ioff2,i),two,Ja(jJvs21+2+m*NVb,i),NVa)
        call vecadd(Za(ioff2,i),Za(ioff2,i),Ka(jJvs21+2+m*NVb,i),NVa)
        call vecadd(Za(ioff2,i),Za(ioff2,i),Ka(jJov4+2+m*NVb,i),NVa)
        ioff2 = ioff2 + NVa
      enddo
    endif ! sf_xcis

  enddo
  deallocate(tmp)
  deallocate(tmp2)
  return
end subroutine sasfemake_a2e

end module sasfcis_util
