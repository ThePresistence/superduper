!
! created by jie liu on 09/16
!
module sasfcis_z

use es_global
use cis_solve
use sasfcis_util
use cisz_solve

private

integer npair,npair1
integer, dimension(:), allocatable :: occ,imo,jmo

public sasfcis_nac_z 

contains

!======================================================== 
subroutine sasfcis_nac_z(V)
  implicit none
  real*8 V(*)

  if(iprtz.ge.5) then
    write(nb6,*)
    write(nb6,*) "**************************************"
    write(nb6,*) "*           ITER_NAC_Z                *"
    write(nb6,*) "**************************************"
    write(nb6,*)
  endif
  
  call sasfcis_nac_z_init(V)

  call sasfcisz_nac_rhs(V)

  call sasfcisz_iter(V)

  call sasfcis_nac_z_exit

  return
end subroutine sasfcis_nac_z

!======================================================== 
subroutine sasfcis_nac_z_init(V)
  implicit none
  integer i
  real*8  V(*)

  allocate(occ(nbas))
  if(sref) then
    do i=1,js1
      occ(i) = 2
    enddo
    do i=js2,nbas
      occ(i) = 0
    enddo
  else
    do i=1,NOb
      occ(i) = 2
    enddo
    occ(js1) = 1
    occ(js2) = 1
    do i=jv,nbas
      occ(i) = 0
    enddo
  endif
  allocate(imo(LM4))
  allocate(jmo(LM4))

  call mopair(V(jEa))
  if(allocated(Tv)) deallocate(Tv)
  if(allocated(At)) deallocate(At)
  allocate(Tv(npair1,maxdav))
  allocate(At(npair1,maxdav))
  allocate(Zr(npair1,nz))
  call vecinit(Tv,npair1*maxdav,0.d0)
  call vecinit(At,npair1*maxdav,0.d0)
  call vecinit(Zr,npair1*nz,0.d0)

  return
end subroutine sasfcis_nac_z_init

!======================================================== 
subroutine sasfcis_nac_z_exit
  if(allocated(Tv)) deallocate(Tv)
  if(allocated(At)) deallocate(At)
  if(allocated(Zr)) deallocate(Zr)
  if(allocated(occ)) deallocate(occ)
  if(allocated(imo)) deallocate(imo)
  if(allocated(jmo)) deallocate(jmo)
end subroutine sasfcis_nac_z_exit

!======================================================== 
subroutine sasfcisz_nac_rhs(V)
  implicit none
  integer i,j,np,nr,m,njs
  real*8  V(*),dot1,dot2,dot3,dot4
  real*8, dimension(:), allocatable :: Fa
  real*8, dimension(:,:), allocatable :: Pa,Ja,Ka,Zw
  
  np  = 3
  if(sref) np = 4
  if(sf_xcis) np = 7
  nr  = npr
  njs = 4
  allocate(Zw(n2,nz))
  call vecinit(Zw,n2*nz,0.d0)

  ! (pq|sisj),(psi|qsj) (ij=11,12,22,21)
  allocate(Pa(n2,njs),Js(n2,njs),Ks(n2,njs))
  call sasfemake_sao(Pa,V(jCa))
  call makejk(Js,Js,Ks,Ks,Pa,Pa,V(jW),njs-1,nbas,LM6,1,I1n,cfacJ,cfacK)
  call sasfemake_smo(Js,Ks,V(jCa),njs-1)
  call mattrans(Js(1,4),Js(1,2),nbas,nbas)
  call mattrans(Ks(1,4),Ks(1,2),nbas,nbas)
  deallocate(Pa)

  ! excited state densities
  allocate(Pa(n2*nz,np),Ja(n2*nz,np),Ka(n2*nz,np))

! 1. CSF part (P_ij*U_ij^x)
  call sasfemake_pmo1_nonsym(Pa,Pa(1,2),Xv,V(jCa),nz)
  call vecadd(Pa,Pa,Pa(1,2),n2*nz)
  if(iprtz.ge.5) then
    do i=1,nz
    write(nb6,*) "nonsymmetry one particle density",i
    call matprnt(Pa(1+(i-1)*n2,1),nbas,nbas,6)  
    enddo
  endif
  if(do_nac) then
    m = 1
    do i=1,nex
      do j=i,nex
        call vecadd3(Zw(1,m),Es(j,1)-Es(i,1),Pa(1+(m-1)*n2,1),n2)
        m = m + 1
      enddo
    enddo
  endif

  call sasfemake_pmo(Pa,Xv,V(jCa),np,nz,1) 
  if(iprtz.ge.5) then
    do i=1,nz
    write(nb6,*) "one particle density",i 
    call matprnt(Pa(1+(i-1)*n2,1),nbas,nbas,6)
    enddo
  endif

! 2. one electron part P_pq*F_rq*U_rp^x+P_pq*Fqr*U_rq^x
  ! Fock matrix in mo basis
  allocate(Fa(n2))
  call square(V(jFa),Fa,nbas,nbas,LM4)
  call matmult(V(jCa),Fa,Ja,nbas,nbas,nbas,LM2,nbas,nbas,2)
  call matmult(Ja,V(jCa),Fa,nbas,nbas,nbas,nbas,LM2,nbas,1)

    !call vecinit(Pa(1,1),nz*n2,0.d0)
  call sasfcisz_nac_rhs_diag1(Zw,Fa,Pa(1,1),nz,one)       ! P_pq * F_rq
  call sasfcisz_nac_rhs_diag1(Zw,Ks,Pa(1,2),nz,-one)      ! P_pq * (rs1|qs1)
  call sasfcisz_nac_rhs_diag1(Zw,Ks(1,3),Pa(1,3),nz,-one) ! P_pq * (rs2|qs2)
  if(sf_xcis) then
    call sasfcisz_nac_rhs_diag1(Zw,Ks(1,2),Pa(1,6),nz,-one) ! P_pq * (rs1|qs2)
    call sasfcisz_nac_rhs_diag1(Zw,Js(1,1),Pa(1,4),nz, one) ! P_pq * (rq|s1s1)
    call sasfcisz_nac_rhs_diag1(Zw,Js(1,2),Pa(1,7),nz, one) ! P_pq * (rq|s1s2)
    call sasfcisz_nac_rhs_diag1(Zw,Js(1,3),Pa(1,5),nz, one) ! P_pq * (rq|s2s2)
  else if(sref) then
    call sasfcisz_nac_rhs_diag1(Zw,Js,Pa(1,4),nz,-one)      ! P_pq * (rq|s1s1)
    call sasfcisz_nac_rhs_diag1(Zw,Js(1,3),Pa(1,4),nz,one)  ! P_pq * (rq|s2s2)
  endif
  deallocate(Fa)

! 3. one electron part F_rq[P_pq]*U_rq^x
  call mo2ao(Pa,V(jCa),nz*np)
  call makejk(Ja,Ja,Ka,Ka,Pa,Pa,V(jW),nz*np,nbas,LM6,1,I1n,cfacJ,cfacK)
  call ao2mo(Ja,V(jCa),nz*np)
  call ao2mo(Ka,V(jCa),nz*np)
  call sasfcisz_nac_rhs_diag2(Zw,Ja,Ka,nz)
  deallocate(Pa)
  deallocate(Ja)
  deallocate(Ka)

! 4. two electron part A_rq,uv \times R_rq \times R_uv
  allocate(Pa(n2,nex*nr),Ja(n2,nex*nr),Ka(n2,nex*nr))
  call sasfemake_rao(Pa,Xv,V(jCa),nex)
  call makejk(Ja,Ja,Ka,Ka,Pa,Pa,V(jW),nex*nr,nbas,lm6,1,I1n,cfacJ,cfacK)
  call ao2mo(Ja,V(jCa),nex*nr)
  call ao2mo(Ka,V(jCa),nex*nr)
  call sasfcisz_nac_rhs_offdiag(Zw,Ja,Ka,Js,Ks,Xv,nz)
  deallocate(Pa)
  deallocate(Ja)
  deallocate(Ka)
  deallocate(Js)
  deallocate(Ks)

! 5. r.h.s for ni > nj
  call sasfcisz_rhs(V,Zv,Zw,V(jCa),V(jEa),nz)

  if(iprtz.ge.5) then
    do i=1,nz
      write(nb6,*) "r.h.s in sasfcis_nac",i
      call matprnt(Zw(1,i),nbas,nbas,6)
    enddo
  endif
  deallocate(Zw)
  return
end subroutine sasfcisz_nac_rhs

!======================================================== 
subroutine sasfcisz_iter(V)
  implicit none
  logical success
  integer i,j
  real*8  V(*)
  real*8, dimension(:), allocatable :: H

  niter = 0
  success = .false.

  ! initial guess
  allocate(H(npair1))
  call sasfcisz_guess(Tv,Zr,V(jEa),H,nz)

  ! Enter the iterative diagonalizer
  do while ((.not.success).and.(niter.le.kitdav))
    niter = niter + 1
    call sasfcisz_make_a(At(1,ns-nl+1),Tv(1,ns-nl+1),V(jCa),V(jEa),V)
    call z_solve(H,npair1)
    if (nl.eq.0) success=.true.
  enddo
  if (.not.success) then
    write(nb6,*) "Max Davidson Iteration: ", kitdav
    stop
  endif
  deallocate(H)

  call pair2matrix(Zv,Zr,nz)
  if(iprtz.gt.5) then
    do i=1,nz
      write(nb6,*) "Z matrix",i 
      call matprnt(Zv(1,i),nbas,nbas,6)
    enddo
  endif

  return
end subroutine sasfcisz_iter

!======================================================== 
subroutine sasfcisz_guess(G,B,Ea,H,n)
  implicit none
  integer i,k,l,m,n
  real*8 dot,tol
  real*8 G(npair1,*),Ea(*),B(nbas,nbas,*),H(*)

  tol = npair1*4d-4*1d1**(2*(1-NrmDav))
  ! mo orbital energy difference
  do i=1,npair1
    k = imo(i)
    l = jmo(i)
    H(i) = Ea(l)-Ea(k)
  enddo

  ns = 0
  do m=1,n
    do i=1,npair1
      k = imo(i)
      l = jmo(i)
      G(i,m) = Zr(i,m) / H(i)
    enddo
    do n=1,ns
      call vecdot(dot,G(1,m),G(1,n),npair1)
      call vecadd2(G(1,m),-dot,G(1,n),G(1,m),npair1)
    enddo
    call vecdot(dot,G(1,m),G(1,m),npair1)
    if(dot.ge.tol) then
      call vecscale(G(1,m),npair1,1/sqrt(dot))
      ns = ns + 1
    endif
  enddo
  nl = ns
end subroutine sasfcisz_guess

!======================================================== 
subroutine sasfcisz_rhs(V,Za,B,Ca,Ea,n)
  implicit none
  integer n,i,k,l,m
  real*8  V(*),Za(n2,*),B(nbas,nbas,*),Ca(LM2,*),Ea(*)
  real*8, dimension(:,:), allocatable :: Pa,Ja,Ka

  call vecinit(Za,n2*n,0.d0)
  do m=1,n
    do i=npair1+1,npair
      k = imo(i)
      l = jmo(i)
      Za(k+(l-1)*nbas,m) = (B(k,l,m)-B(l,k,m)) / (Ea(l)-Ea(k))
    enddo
  enddo

  allocate(Pa(n2,n),Ja(n2,n),Ka(n2,n))
  call veccopy(Pa,Za,n2*n)  
  call mo2ao(Pa,Ca,n)
  call makejk(Ja,Ja,Ka,Ka,Pa,Pa,V(jW),n,nbas,LM6,1,I1n,cfacJ,cfacK)
  call ao2mo(Ja,Ca,n)
  call ao2mo(Ka,Ca,n)
  do m=1,n
  do i=1,npair1
    k = imo(i)
    l = jmo(i)
    Zr(i,m) = B(k,l,m) - B(l,k,m) &
            + (occ(l)-occ(k))*(two*Ja(k+(l-1)*nbas,m) &
            + f12*Ka(k+(l-1)*nbas,m)+f12*Ka(l+(k-1)*nbas,m))
  enddo
  if(iprtz.ge.5) then
    write(nb6,*) "r.h.s.",m
    write(nb6,'(5F10.6)') Zr(1:NOb,m)
    write(nb6,'(5F10.6)') Zr(1+NOb:2*NOb,m)
    do i=1,NVa
    write(nb6,'(7F10.6)') Zr(2*NOb+1+(i-1)*NOa:2*NOb+i*NOa,m)
    enddo
  endif
  enddo

  deallocate(Pa)
  deallocate(Ja)
  deallocate(Ka)
  return
end subroutine sasfcisz_rhs

!======================================================== 
subroutine sasfcisz_make_a(Az,R,Ca,Ea,V)
  implicit none
  integer i,j,m,n
  real*8  Az(npair1,*),R(npair1,*),Ca(LM2,*),Ea(*),V(*)
  real*8,dimension(:,:,:),allocatable :: Ja,Ka
  real*8,dimension(:,:),allocatable :: Pa
  
  allocate(Pa(n2,nl))
  allocate(Ja(nbas,nbas,nl),Ka(nbas,nbas,nl))

  call vecinit(Pa,n2*nl,0.d0)
  call pair2matrix(Pa,R,nl)
  call mo2ao(Pa,Ca,nl)
  call makejk(Ja,Ja,Ka,Ka,Pa,Pa,V(jW),nl,nbas,lm6,1,I1n,cfacJ,cfacK)
  deallocate(Pa)
  call ao2mo(Ja,Ca,nl)
  call ao2mo(Ka,Ca,nl)  

  do i=1,nl
    do j=1,npair1
      m = imo(j)
      n = jmo(j)
      Az(j,i) = Az(j,i) + (Ea(n)-Ea(m))*R(j,i) &
              + (Occ(n)-Occ(m))*(-f12*Ka(m,n,i)-f12*Ka(n,m,i)-two*Ja(m,n,i)) 
    enddo
  enddo 
  deallocate(Ja)
  deallocate(Ka)
  return
end subroutine sasfcisz_make_a

!======================================================== 
subroutine sasfcisz_nac_rhs_diag1(PF,Fa,Pa,n,scale)
  implicit none
  integer n,i
  real*8  scale
  real*8  PF(n2,*),Fa(*),Pa(n2,*)
  real*8, dimension(:), allocatable :: tmp
  
  allocate(tmp(n2))
  do i=1,n
    call matmult(Fa,Pa(1,i),tmp,nbas,nbas,nbas,nbas,nbas,nbas,3)
    call vecadd3(PF(1,i),scale,tmp,n2)
    call matmult(Fa,Pa(1,i),tmp,nbas,nbas,nbas,nbas,nbas,nbas,2)
    call vecadd3(PF(1,i),scale,tmp,n2)
  enddo
  deallocate(tmp)
  return
end subroutine sasfcisz_nac_rhs_diag1

!======================================================== 
subroutine sasfcisz_nac_rhs_diag2(Za,Ja,Ka,n)
  implicit none
  integer n,i,j,m
  real*8  Za(n2,*),Ja(n2,n,*),Ka(n2,n,*)

  do m=1,n
    do i=1,nbas
      do j=1,NOa
        Za(i+(j-1)*nbas,m) = Za(i+(j-1)*nbas,m)  &
                           + occ(j)*(Ja(i+(j-1)*nbas,m,1)+Ja(j+(i-1)*nbas,m,1)) &
                           + occ(j)*f12*(Ka(i+(j-1)*nbas,m,1)+Ka(j+(i-1)*nbas,m,1))
        if(j.eq.js1) then
          Za(i+(j-1)*nbas,m) = Za(i+(j-1)*nbas,m) &
                   - Ka(i+(j-1)*nbas,m,2) - Ka(j+(i-1)*nbas,m,2)
        endif
        if(j.eq.js2) then
          Za(i+(j-1)*nbas,m) = Za(i+(j-1)*nbas,m) &
                   - Ka(i+(j-1)*nbas,m,3) - Ka(j+(i-1)*nbas,m,3)
        endif
        if(sf_xcis) then
          if(j.eq.js1) then
            Za(i+(j-1)*nbas,m) = Za(i+(j-1)*nbas,m) &
                     - Ka(i+(js2-1)*nbas,m,6) + Ja(i+(js2-1)*nbas,m,7) &
                     + Ja(i+(j-1)*nbas,m,4) + Ja(j+(i-1)*nbas,m,4)
          endif
          if(j.eq.js2) then
            Za(i+(j-1)*nbas,m) = Za(i+(j-1)*nbas,m) &
                     - Ka(js1+(i-1)*nbas,m,6) + Ja(js1+(i-1)*nbas,m,7) &
                     + Ja(i+(j-1)*nbas,m,5) + Ja(j+(i-1)*nbas,m,5)
          endif
        else if(sref) then
          if(j.eq.js1) then
            Za(i+(j-1)*nbas,m) = Za(i+(j-1)*nbas,m) &
                     - Ja(i+(j-1)*nbas,m,4) - Ja(j+(i-1)*nbas,m,4)
          endif
          if(j.eq.js2) then
            Za(i+(j-1)*nbas,m) = Za(i+(j-1)*nbas,m) &
                     + Ja(i+(j-1)*nbas,m,4) + Ja(j+(i-1)*nbas,m,4)
          endif
        endif
      enddo
    enddo
  enddo
  return 
end subroutine sasfcisz_nac_rhs_diag2

!======================================================== 
subroutine sasfcisz_nac_rhs_offdiag(Za,Ja,Ka,Jss,Kss,R,n)
  implicit none
  integer n,i,j,k,m,ia,ib,ic,ir
  real*8  Za(n2,*),Ja(n2*npr,*),Ka(n2*npr,*),Jss(n2,*),Kss(n2,*),R(NOV,*)
  real*8, dimension(:), allocatable :: tmp,tmp2

  allocate(tmp(n2),tmp2(n2))
  if(sf_xcis) then
    do i=1,nex
      call vecadd3(Ka(jJov3,i),two,Ja(jJov3,i),n2)
      call vecadd3(Ka(jJov4,i),two,Ja(jJov4,i),n2)
    enddo
  endif    

  m = 1
  do i=1,nex
    do j=i,nex
    ! s1->s2
    ia = js2
    ib = js1
    ir = 1
    call vecinit(tmp,n2,0.d0)
    call vecadd3(tmp,one*R(1,i),Kss(1,4),n2)
    call vecadd3(tmp,one*R(2,i),Kss(1,2),n2)
    call vecadd3(tmp,sq2*R(3,i),Kss(1,1),n2)
    call vecadd3(tmp,sq2,Ka(jJos1,i),n2)
    call vecadd3(tmp,sq2,Ka(jJos2,i),n2)
    call vecadd3(tmp,sq2,Ka(jJvs1,i),n2)
    call vecadd3(tmp,sq2,Ka(jJvs2,i),n2)
    call vecadd3(tmp,sq3,Ja(jJov1,i),n2)
    call vecadd3(tmp,one,Ja(jJov2,i),n2)
    call vecadd3(tmp,two,Ka(jJov2,i),n2)
    !if(sf_xcis) then
    !  call vecadd3(tmp,-R(3,j)/R(1,j),Ka(jJov4,i),n2)
    !endif
    if(i.eq.j) then
      call sasfciszmat(Za(1,m),tmp,R(1,j),1,1,ia,ib,ir,two)
      if(sf_xcis) then
        call sasfciszmat(Za(1,m),Ka(jJov4,i),R(1,j),1,1,ia,ib,3,-two)
      endif
    else
      call sasfciszmat(Za(1,m),tmp,R(1,j),1,1,ia,ib,ir,one)
      call vecinit(tmp,n2,0.d0)
      call vecadd3(tmp,one*R(1,j),Kss(1,4),n2)
      call vecadd3(tmp,one*R(2,j),Kss(1,2),n2)
      call vecadd3(tmp,sq2*R(3,j),Kss(1,1),n2)
      call vecadd3(tmp,sq2,Ka(jJos1,j),n2)
      call vecadd3(tmp,sq2,Ka(jJos2,j),n2)
      call vecadd3(tmp,sq2,Ka(jJvs1,j),n2)
      call vecadd3(tmp,sq2,Ka(jJvs2,j),n2)
      call vecadd3(tmp,sq3,Ja(jJov1,j),n2)
      call vecadd3(tmp,one,Ja(jJov2,j),n2)
      call vecadd3(tmp,two,Ka(jJov2,j),n2)
      !if(sf_xcis) then
      !  call vecadd3(tmp,-R(3,i)/R(1,i),Ka(jJov4,j),n2)
      !endif
      call sasfciszmat(Za(1,m),tmp,R(1,i),1,1,ia,ib,ir,one)
      if(sf_xcis) then
        call sasfciszmat(Za(1,m),Ka(jJov4,j),R(1,i),1,1,ia,ib,3,-one)
        call sasfciszmat(Za(1,m),Ka(jJov4,i),R(1,j),1,1,ia,ib,3,-one)
      endif
    endif

    ! s2->s1
    ia = js1 
    ib = js2
    ir = 2
    call vecinit(tmp,n2,0.d0)
    call vecadd3(tmp,one*R(1,i),Kss(1,4),n2)
    call vecadd3(tmp,one*R(2,i),Kss(1,2),n2)
    call vecadd3(tmp,sq2*R(3,i),Kss(1,1),n2)
    call vecadd3(tmp,sq2,Ka(jJos1,i),n2)
    call vecadd3(tmp,sq2,Ka(jJos2,i),n2)
    call vecadd3(tmp,sq2,Ka(jJvs1,i),n2)
    call vecadd3(tmp,sq2,Ka(jJvs2,i),n2)
    call vecadd3(tmp,-sq3,Ja(jJov1,i),n2)
    call vecadd3(tmp,-sq3,Ka(jJov1,i),n2)
    call vecadd3(tmp,-one,Ja(jJov2,i),n2)
    call vecadd3(tmp,one,Ka(jJov2,i),n2)
    !if(sf_xcis) then
    !  call vecadd3(tmp,R(3,j)/R(2,j),Ka(jJov3,i),n2)
    !endif
    if(i.eq.j) then
      call sasfciszmat(Za(1,m),tmp,R(1,j),1,1,ia,ib,ir,two)
      if(sf_xcis) then
        call sasfciszmat(Za(1,m),Ka(jJov3,i),R(1,j),1,1,ia,ib,3,two)
      endif
    else
      call sasfciszmat(Za(1,m),tmp,R(1,j),1,1,ia,ib,ir,one)
      call vecinit(tmp,n2,0.d0)
      call vecadd3(tmp,one*R(1,j),Kss(1,4),n2)
      call vecadd3(tmp,one*R(2,j),Kss(1,2),n2)
      call vecadd3(tmp,sq2*R(3,j),Kss(1,1),n2)
      call vecadd3(tmp,sq2,Ka(jJos1,j),n2)
      call vecadd3(tmp,sq2,Ka(jJos2,j),n2)
      call vecadd3(tmp,sq2,Ka(jJvs1,j),n2)
      call vecadd3(tmp,sq2,Ka(jJvs2,j),n2)
      call vecadd3(tmp,-sq3,Ja(jJov1,j),n2)
      call vecadd3(tmp,-sq3,Ka(jJov1,j),n2)
      call vecadd3(tmp,-one,Ja(jJov2,j),n2)
      call vecadd3(tmp,one,Ka(jJov2,j),n2)
      !if(sf_xcis) then
      !  call vecadd3(tmp,R(3,i)/R(2,i),Ka(jJov3,j),n2)
      !endif
      call sasfciszmat(Za(1,m),tmp,R(1,i),1,1,ia,ib,ir,one)
      if(sf_xcis) then
        call sasfciszmat(Za(1,m),Ka(jJov3,j),R(1,i),1,1,ia,ib,3,one)
        call sasfciszmat(Za(1,m),Ka(jJov3,i),R(1,j),1,1,ia,ib,3,one)
      endif
    endif

    ! s1->s1
    ia = js1
    ib = js2
    ir = 3
    call vecinit(tmp,n2,0.d0)
    call vecadd3(tmp,sq2*R(1,i),Kss(1,4),n2)
    call vecadd3(tmp,sq2*R(2,i),Kss(1,2),n2)
    call vecadd3(tmp,    R(3,i),Kss(1,1),n2)
    call vecadd3(tmp,   -R(3,i),Kss(1,3),n2)
    call vecadd3(tmp,one,Ka(jJos1,i),n2)
    call vecadd3(tmp,one,Ka(jJos2,i),n2)
    call vecadd3(tmp,one,Ka(jJvs1,i),n2)
    call vecadd3(tmp,one,Ka(jJvs2,i),n2)
    call vecadd3(tmp,sq2,Ka(jJov2,i),n2)
    call vecinit(tmp2,n2,0.d0)
    call vecadd3(tmp2,  -R(3,i),Kss(1,1),n2)
    call vecadd3(tmp2,-one,Ka(jJos1,i),n2)
    call vecadd3(tmp2,-one,Ka(jJos2,i),n2)
    call vecadd3(tmp2,-one,Ka(jJvs1,i),n2)
    call vecadd3(tmp2,-one,Ka(jJvs2,i),n2)
    call vecadd3(tmp2,ssq6,Ka(jJov1,i),n2)
    call vecadd3(tmp2,-ssq2,Ka(jJov2,i),n2)
    !if(sf_xcis) then
    !  call vecadd3(tmp, -ssq2*R(1,j)/R(3,j),Ka(jJov3,i),n2)
    !  call vecadd3(tmp2, ssq2*R(1,j)/R(3,j),Ka(jJov3,i),n2)
    !  call vecadd3(tmp,  ssq2*R(2,j)/R(3,j),Ka(jJov4,i),n2)
    !  call vecadd3(tmp2,-ssq2*R(2,j)/R(3,j),Ka(jJov4,i),n2)
    !endif
    if(i.eq.j) then
      call vecadd3(tmp,   f12*R(3,i),Kss(1,3),n2)
      call vecadd3(tmp2,  f12*R(3,i),Kss(1,1),n2)
      call sasfciszmat(Za(1,m),tmp ,R(1,j),1,1,ia,ia,ir,two)
      call sasfciszmat(Za(1,m),tmp2,R(1,j),1,1,ib,ib,ir,two)
      if(sf_xcis) then
        call sasfciszmat(Za(1,m),Ka(jJov3,i),R(1,j),1,1,ia,ia,1,-two*ssq2)
        call sasfciszmat(Za(1,m),Ka(jJov3,i),R(1,j),1,1,ib,ib,1,two*ssq2)
        call sasfciszmat(Za(1,m),Ka(jJov4,i),R(1,j),1,1,ia,ia,2,two*ssq2)
        call sasfciszmat(Za(1,m),Ka(jJov4,i),R(1,j),1,1,ib,ib,2,-two*ssq2)
      endif
    else
      call sasfciszmat(Za(1,m),tmp ,R(1,j),1,1,ia,ia,ir,one)
      call sasfciszmat(Za(1,m),tmp2,R(1,j),1,1,ib,ib,ir,one)
      call vecinit(tmp,n2,0.d0)
      call vecadd3(tmp,sq2*R(1,j),Kss(1,4),n2)
      call vecadd3(tmp,sq2*R(2,j),Kss(1,2),n2)
      call vecadd3(tmp,    R(3,j),Kss(1,1),n2)
      call vecadd3(tmp,one,Ka(jJos1,j),n2)
      call vecadd3(tmp,one,Ka(jJos2,j),n2)
      call vecadd3(tmp,one,Ka(jJvs1,j),n2)
      call vecadd3(tmp,one,Ka(jJvs2,j),n2)
      call vecadd3(tmp,sq2,Ka(jJov2,j),n2)
      call vecinit(tmp2,n2,0.d0)
      call vecadd3(tmp2,-one,Ka(jJos1,j),n2)
      call vecadd3(tmp2,-one,Ka(jJos2,j),n2)
      call vecadd3(tmp2,-one,Ka(jJvs1,j),n2)
      call vecadd3(tmp2,-one,Ka(jJvs2,j),n2)
      call vecadd3(tmp2,ssq6,Ka(jJov1,j),n2)
      call vecadd3(tmp2,-ssq2,Ka(jJov2,j),n2)
      !if(sf_xcis) then
      !  call vecadd3(tmp, -ssq2*R(1,i)/R(3,i),Ka(jJov3,j),n2)
      !  call vecadd3(tmp2, ssq2*R(1,i)/R(3,i),Ka(jJov3,j),n2)
      !  call vecadd3(tmp,  ssq2*R(2,i)/R(3,i),Ka(jJov4,j),n2)
      !  call vecadd3(tmp2,-ssq2*R(2,i)/R(3,i),Ka(jJov4,j),n2)
      !endif
      call sasfciszmat(Za(1,m),tmp ,R(1,i),1,1,ia,ia,ir,one)
      call sasfciszmat(Za(1,m),tmp2,R(1,i),1,1,ib,ib,ir,one)
      if(sf_xcis) then
        call sasfciszmat(Za(1,m),Ka(jJov3,i),R(1,j),1,1,ia,ia,1,-ssq2)
        call sasfciszmat(Za(1,m),Ka(jJov3,i),R(1,j),1,1,ib,ib,1,ssq2)
        call sasfciszmat(Za(1,m),Ka(jJov4,i),R(1,j),1,1,ia,ia,2,ssq2)
        call sasfciszmat(Za(1,m),Ka(jJov4,i),R(1,j),1,1,ib,ib,2,-ssq2)
        call sasfciszmat(Za(1,m),Ka(jJov3,j),R(1,i),1,1,ia,ia,1,-ssq2)
        call sasfciszmat(Za(1,m),Ka(jJov3,j),R(1,i),1,1,ib,ib,1,ssq2)
        call sasfciszmat(Za(1,m),Ka(jJov4,j),R(1,i),1,1,ia,ia,2,ssq2)
        call sasfciszmat(Za(1,m),Ka(jJov4,j),R(1,i),1,1,ib,ib,2,-ssq2)
      endif
    endif


    ! occ->s1
    ia = js1
    ib = 1
    ic = js2
    ir = os1 
    call vecinit(tmp,n2,0.d0)
    call vecadd3(tmp,sq2*R(1,i),Kss(1,4),n2)
    call vecadd3(tmp,sq2*R(2,i),Kss(1,2),n2)
    call vecadd3(tmp,    R(3,i),Kss(1,1),n2)
    call vecadd3(tmp,   -R(3,i),Kss(1,3),n2)
    call vecadd3(tmp,one,Ka(jJos1,i),n2)
    call vecadd3(tmp,one,Ka(jJos2,i),n2)
    call vecadd3(tmp,-one,Ja(jJos2,i),n2)
    call vecadd3(tmp,one,Ka(jJvs1,i),n2)
    call vecadd3(tmp,two,Ka(jJvs2,i),n2)
    call vecadd3(tmp,one,Ja(jJvs2,i),n2)
    call vecadd3(tmp,-ssq6,Ka(jJov1,i),n2)
    call vecadd3(tmp,-ssq6,Ja(jJov1,i),n2)
    call vecadd3(tmp,ssq2,Ka(jJov2,i),n2)
    call vecadd3(tmp,-ssq2,Ja(jJov2,i),n2)
    if(i.eq.j) then
      call sasfciszmat(Za(1,m),tmp, R(1,j),1,NOb,ia,ib,ir,two)
      call sasfciszmat(Za(1,m),Ja(jJos12,i),R(1,j),1,NOb,ic,ib,ir,two)
      if(sf_xcis) then
        call sasfciszmat(Za(1,m),Ka(jJov4,i),R(1,j),1,NOb,ic,ib,ir,-two)
      endif
    else
      call sasfciszmat(Za(1,m),tmp, R(1,j),1,NOb,ia,ib,ir,one)
      call sasfciszmat(Za(1,m),Ja(jJos12,i),R(1,j),1,NOb,ic,ib,ir,one)
      call vecinit(tmp,n2,0.d0)
      call vecadd3(tmp,sq2*R(1,j),Kss(1,4),n2)
      call vecadd3(tmp,sq2*R(2,j),Kss(1,2),n2)
      call vecadd3(tmp,    R(3,j),Kss(1,1),n2)
      call vecadd3(tmp,   -R(3,j),Kss(1,3),n2)
      call vecadd3(tmp,one,Ka(jJos1,j),n2)
      call vecadd3(tmp,one,Ka(jJos2,j),n2)
      call vecadd3(tmp,-one,Ja(jJos2,j),n2)
      call vecadd3(tmp,one,Ka(jJvs1,j),n2)
      call vecadd3(tmp,two,Ka(jJvs2,j),n2)
      call vecadd3(tmp,one,Ja(jJvs2,j),n2)
      call vecadd3(tmp,-ssq6,Ka(jJov1,j),n2)
      call vecadd3(tmp,-ssq6,Ja(jJov1,j),n2)
      call vecadd3(tmp,ssq2,Ka(jJov2,j),n2)
      call vecadd3(tmp,-ssq2,Ja(jJov2,j),n2)
      call sasfciszmat(Za(1,m),tmp, R(1,i),1,NOb,ia,ib,ir,one)
      call sasfciszmat(Za(1,m),Ja(jJos12,j),R(1,i),1,NOb,ic,ib,ir,one)
      if(sf_xcis) then
        call sasfciszmat(Za(1,m),Ka(jJov4,i),R(1,j),1,NOb,ic,ib,ir,-one)
        call sasfciszmat(Za(1,m),Ka(jJov4,j),R(1,i),1,NOb,ic,ib,ir,-one)
      endif
    endif

    ! occ->s2
    ia = js2
    ib = 1
    ic = js1
    ir = os2 
    call vecinit(tmp,n2,0.d0)
    call vecadd3(tmp,sq2*R(1,i),Kss(1,4),n2)
    call vecadd3(tmp,sq2*R(2,i),Kss(1,2),n2)
    call vecadd3(tmp,    R(3,i),Kss(1,1),n2)
    call vecadd3(tmp,   -R(3,i),Kss(1,3),n2)
    call vecadd3(tmp,one,Ka(jJos2,i),n2)
    call vecadd3(tmp,one,Ka(jJos1,i),n2)
    call vecadd3(tmp,-one,Ja(jJos1,i),n2)
    call vecadd3(tmp,one,Ka(jJvs2,i),n2)
    call vecadd3(tmp,two,Ka(jJvs1,i),n2)
    call vecadd3(tmp,one,Ja(jJvs1,i),n2)
    call vecadd3(tmp,ssq6,Ja(jJov1,i),n2)
    call vecadd3(tmp,sq2,Ka(jJov2,i),n2)
    call vecadd3(tmp,ssq2,Ja(jJov2,i),n2)
    if(i.eq.j) then
      call sasfciszmat(Za(1,m),tmp, R(1,j),1,NOb,ia,ib,ir,two)
      call sasfciszmat(Za(1,m),Ja(jJos21,i),R(1,j),1,NOb,ic,ib,ir,two)
      if(sf_xcis) then
        call sasfciszmat(Za(1,m),Ka(jJov3,i),R(1,j),1,NOb,ic,ib,ir,-two)
      endif
    else
      call sasfciszmat(Za(1,m),tmp, R(1,j),1,NOb,ia,ib,ir,one)
      call sasfciszmat(Za(1,m),Ja(jJos21,i),R(1,j),1,NOb,ic,ib,ir,one)
      call vecinit(tmp,n2,0.d0)
      call vecadd3(tmp,sq2*R(1,j),Kss(1,4),n2)
      call vecadd3(tmp,sq2*R(2,j),Kss(1,2),n2)
      call vecadd3(tmp,    R(3,j),Kss(1,1),n2)
      call vecadd3(tmp,   -R(3,j),Kss(1,3),n2)
      call vecadd3(tmp,one,Ka(jJos2,j),n2)
      call vecadd3(tmp,one,Ka(jJos1,j),n2)
      call vecadd3(tmp,-one,Ja(jJos1,j),n2)
      call vecadd3(tmp,one,Ka(jJvs2,j),n2)
      call vecadd3(tmp,two,Ka(jJvs1,j),n2)
      call vecadd3(tmp,one,Ja(jJvs1,j),n2)
      call vecadd3(tmp,ssq6,Ja(jJov1,j),n2)
      call vecadd3(tmp,sq2,Ka(jJov2,j),n2)
      call vecadd3(tmp,ssq2,Ja(jJov2,j),n2)
      call sasfciszmat(Za(1,m),tmp, R(1,i),1,NOb,ia,ib,ir,one)
      call sasfciszmat(Za(1,m),Ja(jJos21,j),R(1,i),1,NOb,ic,ib,ir,one)
      if(sf_xcis) then
        call sasfciszmat(Za(1,m),Ka(jJov3,i),R(1,j),1,NOb,ic,ib,ir,-one)
        call sasfciszmat(Za(1,m),Ka(jJov3,j),R(1,i),1,NOb,ic,ib,ir,-one)
      endif
    endif


    ! s1->vir
    ia = jv
    ib = js1
    ic = js2
    ir = vs1 
    call vecinit(tmp,n2,0.d0)
    call vecadd3(tmp,sq2*R(1,i),Kss(1,4),n2)
    call vecadd3(tmp,sq2*R(2,i),Kss(1,2),n2)
    call vecadd3(tmp,    R(3,i),Kss(1,1),n2)
    call vecadd3(tmp,   -R(3,i),Kss(1,3),n2)
    call vecadd3(tmp,one,Ka(jJos1,i),n2)
    call vecadd3(tmp,two,Ka(jJos2,i),n2)
    call vecadd3(tmp,one,Ja(jJos2,i),n2)
    call vecadd3(tmp,one,Ka(jJvs1,i),n2)
    call vecadd3(tmp,one,Ka(jJvs2,i),n2)
    call vecadd3(tmp,-one,Ja(jJvs2,i),n2)
    call vecadd3(tmp,ssq6,Ja(jJov1,i),n2)
    call vecadd3(tmp,sq2,Ka(jJov2,i),n2)
    call vecadd3(tmp,ssq2,Ja(jJov2,i),n2)
    if(i.eq.j) then
      call sasfciszmat(Za(1,m),tmp, R(1,j),NVa,1,ia,ib,ir,two)
      call sasfciszmat(Za(1,m),Ja(jJvs12,i),R(1,j),NVa,1,ia,ic,ir,two)
      if(sf_xcis) then
        call sasfciszmat(Za(1,m),Ka(jJov3,i),R(1,j),NVa,1,ia,ic,ir,two)
      endif
    else
      call sasfciszmat(Za(1,m),tmp, R(1,j),NVa,1,ia,ib,ir,one)
      call sasfciszmat(Za(1,m),Ja(jJvs12,i),R(1,j),NVa,1,ia,ic,ir,one)
      call vecinit(tmp,n2,0.d0)
      call vecadd3(tmp,sq2*R(1,j),Kss(1,4),n2)
      call vecadd3(tmp,sq2*R(2,j),Kss(1,2),n2)
      call vecadd3(tmp,    R(3,j),Kss(1,1),n2)
      call vecadd3(tmp,   -R(3,j),Kss(1,3),n2)
      call vecadd3(tmp,one,Ka(jJos1,j),n2)
      call vecadd3(tmp,two,Ka(jJos2,j),n2)
      call vecadd3(tmp,one,Ja(jJos2,j),n2)
      call vecadd3(tmp,one,Ka(jJvs1,j),n2)
      call vecadd3(tmp,one,Ka(jJvs2,j),n2)
      call vecadd3(tmp,-one,Ja(jJvs2,j),n2)
      call vecadd3(tmp,ssq6,Ja(jJov1,j),n2)
      call vecadd3(tmp,sq2,Ka(jJov2,j),n2)
      call vecadd3(tmp,ssq2,Ja(jJov2,j),n2)
      call sasfciszmat(Za(1,m),tmp, R(1,i),NVa,1,ia,ib,ir,one)
      call sasfciszmat(Za(1,m),Ja(jJvs12,j),R(1,i),NVa,1,ia,ic,ir,one)
      if(sf_xcis) then
        call sasfciszmat(Za(1,m),Ka(jJov3,j),R(1,i),NVa,1,ia,ic,ir,one)
        call sasfciszmat(Za(1,m),Ka(jJov3,i),R(1,j),NVa,1,ia,ic,ir,one)
      endif
    endif

    ! s2->vir                                
    ia = jv
    ib = js2
    ic = js1
    ir = vs2 
    call vecinit(tmp,n2,0.d0)
    call vecadd3(tmp,sq2*R(1,i),Kss(1,4),n2)
    call vecadd3(tmp,sq2*R(2,i),Kss(1,2),n2)
    call vecadd3(tmp,    R(3,i),Kss(1,1),n2)
    call vecadd3(tmp,   -R(3,i),Kss(1,3),n2)
    call vecadd3(tmp,one,Ka(jJos2,i),n2)
    call vecadd3(tmp,two,Ka(jJos1,i),n2)
    call vecadd3(tmp,one,Ja(jJos1,i),n2)
    call vecadd3(tmp,one,Ka(jJvs2,i),n2)
    call vecadd3(tmp,one,Ka(jJvs1,i),n2)
    call vecadd3(tmp,-one,Ja(jJvs1,i),n2)
    call vecadd3(tmp,-ssq6,Ja(jJov1,i),n2)
    call vecadd3(tmp,-ssq6,Ka(jJov1,i),n2)
    call vecadd3(tmp,ssq2,Ka(jJov2,i),n2)
    call vecadd3(tmp,-ssq2,Ja(jJov2,i),n2)
    if(i.eq.j) then
      call sasfciszmat(Za(1,m),tmp, R(1,j),NVa,1,ia,ib,ir,two)
      call sasfciszmat(Za(1,m),Ja(jJvs21,i),R(1,j),NVa,1,ia,ic,ir,two)
      if(sf_xcis) then
        call sasfciszmat(Za(1,m),Ka(jJov4,i),R(1,j),NVa,1,ia,ic,ir,two)
      endif
    else
      call sasfciszmat(Za(1,m),tmp, R(1,j),NVa,1,ia,ib,ir,one)
      call sasfciszmat(Za(1,m),Ja(jJvs21,i),R(1,j),NVa,1,ia,ic,ir,one)
      call vecinit(tmp,n2,0.d0)
      call vecadd3(tmp,sq2*R(1,j),Kss(1,4),n2)
      call vecadd3(tmp,sq2*R(2,j),Kss(1,2),n2)
      call vecadd3(tmp,    R(3,j),Kss(1,1),n2)
      call vecadd3(tmp,   -R(3,j),Kss(1,3),n2)
      call vecadd3(tmp,one,Ka(jJos2,j),n2)
      call vecadd3(tmp,two,Ka(jJos1,j),n2)
      call vecadd3(tmp,one,Ja(jJos1,j),n2)
      call vecadd3(tmp,one,Ka(jJvs2,j),n2)
      call vecadd3(tmp,one,Ka(jJvs1,j),n2)
      call vecadd3(tmp,-one,Ja(jJvs1,j),n2)
      call vecadd3(tmp,-ssq6,Ja(jJov1,j),n2)
      call vecadd3(tmp,-ssq6,Ka(jJov1,j),n2)
      call vecadd3(tmp,ssq2,Ka(jJov2,j),n2)
      call vecadd3(tmp,-ssq2,Ja(jJov2,j),n2)
      call sasfciszmat(Za(1,m),tmp, R(1,i),NVa,1,ia,ib,ir,one)
      call sasfciszmat(Za(1,m),Ja(jJvs21,j),R(1,i),NVa,1,ia,ic,ir,one)
      if(sf_xcis) then
        call sasfciszmat(Za(1,m),Ka(jJov4,i),R(1,j),NVa,1,ia,ic,ir,one)
        call sasfciszmat(Za(1,m),Ka(jJov4,j),R(1,i),NVa,1,ia,ic,ir,one)
      endif
    endif

    ! occ->vir A
    ia = jv
    ib = 1
    ir = ov1 
    call vecinit(tmp,n2,0.d0)
    call vecadd3(tmp,sq3*R(1,i),Jss(1,4),n2)
    call vecadd3(tmp,-sq3*R(2,i),Kss(1,2),n2)
    call vecadd3(tmp,-sq3*R(2,i),Jss(1,2),n2)
    call vecadd3(tmp,ssq6*R(3,i),Kss(1,3),n2)
    call vecadd3(tmp,-ssq6,Ka(jJos1,i),n2)
    call vecadd3(tmp,-ssq6,Ja(jJos1,i),n2)
    call vecadd3(tmp,ssq6,Ja(jJos2,i),n2)
    call vecadd3(tmp,ssq6,Ja(jJvs1,i),n2)
    call vecadd3(tmp,-ssq6,Ka(jJvs2,i),n2)
    call vecadd3(tmp,-ssq6,Ja(jJvs2,i),n2)
    call vecadd3(tmp,f32,Ja(jJov1,i),n2)
    call vecadd3(tmp,one,Ka(jJov1,i),n2)
    call vecadd3(tmp,ssq3,Ja(jJov2,i),n2)
    if(i.eq.j) then
      call sasfciszmat(Za(1,m),tmp, R(1,j),NVa,NOb,ia,ib,ir,two)
    else
      call sasfciszmat(Za(1,m),tmp, R(1,j),NVa,NOb,ia,ib,ir,one)
      call vecinit(tmp,n2,0.d0)
      call vecadd3(tmp,sq3*R(1,j),Jss(1,4),n2)
      call vecadd3(tmp,-sq3*R(2,j),Kss(1,2),n2)
      call vecadd3(tmp,-sq3*R(2,j),Jss(1,2),n2)
      call vecadd3(tmp,ssq6*R(3,j),Kss(1,3),n2)
      call vecadd3(tmp,-ssq6,Ka(jJos1,j),n2)
      call vecadd3(tmp,-ssq6,Ja(jJos1,j),n2)
      call vecadd3(tmp,ssq6,Ja(jJos2,j),n2)
      call vecadd3(tmp,ssq6,Ja(jJvs1,j),n2)
      call vecadd3(tmp,-ssq6,Ka(jJvs2,j),n2)
      call vecadd3(tmp,-ssq6,Ja(jJvs2,j),n2)
      call vecadd3(tmp,f32,Ja(jJov1,j),n2)
      call vecadd3(tmp,one,Ka(jJov1,j),n2)
      call vecadd3(tmp,ssq3,Ja(jJov2,j),n2)
      call sasfciszmat(Za(1,m),tmp, R(1,i),NVa,NOb,ia,ib,ir,one)
    endif

    ! occ->vir B
    ia = jv
    ib = 1
    ir = ov2 
    call vecinit(tmp,n2,0.d0)
    call vecadd3(tmp,one*R(1,i),Jss(1,4),n2)
    call vecadd3(tmp,two*R(1,i),Kss(1,4),n2)
    call vecadd3(tmp,one*R(2,i),Kss(1,2),n2)
    call vecadd3(tmp,-one*R(2,i),Jss(1,2),n2)
    call vecadd3(tmp,sq2*R(3,i),Kss(1,1),n2)
    call vecadd3(tmp,-ssq2*R(3,i),Kss(1,3),n2)
    call vecadd3(tmp,ssq2,Ka(jJos1,i),n2)
    call vecadd3(tmp,-ssq2,Ja(jJos1,i),n2)
    call vecadd3(tmp,ssq2,Ja(jJos2,i),n2)
    call vecadd3(tmp,sq2,Ka(jJos2,i),n2)
    call vecadd3(tmp,ssq2,Ja(jJvs1,i),n2)
    call vecadd3(tmp,sq2,Ka(jJvs1,i),n2)
    call vecadd3(tmp,ssq2,Ka(jJvs2,i),n2)
    call vecadd3(tmp,-ssq2,Ja(jJvs2,i),n2)
    call vecadd3(tmp,f12,Ja(jJov2,i),n2)
    call vecadd3(tmp,one,Ka(jJov2,i),n2)
    call vecadd3(tmp,ssq3,Ja(jJov1,i),n2)
    if(i.eq.j) then
      call sasfciszmat(Za(1,m),tmp, R(1,j),NVa,NOb,ia,ib,ir,two)
    else
      call sasfciszmat(Za(1,m),tmp, R(1,j),NVa,NOb,ia,ib,ir,one)
      call vecinit(tmp,n2,0.d0)
      call vecadd3(tmp,one*R(1,j),Jss(1,4),n2)
      call vecadd3(tmp,two*R(1,j),Kss(1,4),n2)
      call vecadd3(tmp,one*R(2,j),Kss(1,2),n2)
      call vecadd3(tmp,-one*R(2,j),Jss(1,2),n2)
      call vecadd3(tmp,sq2*R(3,j),Kss(1,1),n2)
      call vecadd3(tmp,-ssq2*R(3,j),Kss(1,3),n2)
      call vecadd3(tmp,ssq2,Ka(jJos1,j),n2)
      call vecadd3(tmp,-ssq2,Ja(jJos1,j),n2)
      call vecadd3(tmp,ssq2,Ja(jJos2,j),n2)
      call vecadd3(tmp,sq2,Ka(jJos2,j),n2)
      call vecadd3(tmp,ssq2,Ja(jJvs1,j),n2)
      call vecadd3(tmp,sq2,Ka(jJvs1,j),n2)
      call vecadd3(tmp,ssq2,Ka(jJvs2,j),n2)
      call vecadd3(tmp,-ssq2,Ja(jJvs2,j),n2)
      call vecadd3(tmp,f12,Ja(jJov2,j),n2)
      call vecadd3(tmp,one,Ka(jJov2,j),n2)
      call vecadd3(tmp,ssq3,Ja(jJov1,j),n2)
      call sasfciszmat(Za(1,m),tmp, R(1,i),NVa,NOb,ia,ib,ir,one)
    endif

    if(sf_xcis) then
    ! occ->vir C
    ia = jv
    ib = 1
    ir = ov3 
    call vecinit(tmp,n2,0.d0)
    call vecadd3(tmp,sq2*R(1,i),Jss(1,3),n2)
    call vecadd3(tmp,ssq2*R(1,i),Kss(1,3),n2)
    call vecadd3(tmp,-sq2*R(1,i),Jss(1,1),n2)
    call vecadd3(tmp,-ssq2*R(1,i),Kss(1,1),n2)
    call vecadd3(tmp,two*R(3,i),Jss(1,2),n2)
    call vecadd3(tmp,R(3,i),Kss(1,2),n2)
    call vecadd3(tmp,-two,Ja(jJos21,i),n2)
    call vecadd3(tmp,-one,Ka(jJos21,i),n2)
    call vecadd3(tmp,two,Ja(jJvs12,i),n2)
    call vecadd3(tmp,one,Ka(jJvs12,i),n2)
    call vecadd3(tmp,one,Ka(jJov3,i),n2)
    if(i.eq.j) then
      call sasfciszmat(Za(1,m),tmp, R(1,j),NVa,NOb,ia,ib,ir,two)
    else
      call sasfciszmat(Za(1,m),tmp, R(1,j),NVa,NOb,ia,ib,ir,one)
      call vecinit(tmp,n2,0.d0)
      call vecadd3(tmp,sq2*R(1,j),Jss(1,3),n2)
      call vecadd3(tmp,ssq2*R(1,j),Kss(1,3),n2)
      call vecadd3(tmp,-sq2*R(1,j),Jss(1,1),n2)
      call vecadd3(tmp,-ssq2*R(1,j),Kss(1,1),n2)
      call vecadd3(tmp,two*R(3,j),Jss(1,2),n2)
      call vecadd3(tmp,R(3,j),Kss(1,2),n2)
      call vecadd3(tmp,-two,Ja(jJos21,j),n2)
      call vecadd3(tmp,-one,Ka(jJos21,j),n2)
      call vecadd3(tmp,two,Ja(jJvs12,j),n2)
      call vecadd3(tmp,one,Ka(jJvs12,j),n2)
      call vecadd3(tmp,one,Ka(jJov3,j),n2)
      call sasfciszmat(Za(1,m),tmp, R(1,i),NVa,NOb,ia,ib,ir,one)
    endif

    ! occ->vir D
    ia = jv
    ib = 1
    ir = ov4 
    call vecinit(tmp,n2,0.d0)
    call vecadd3(tmp,sq2*R(2,i),Jss(1,1),n2)
    call vecadd3(tmp,ssq2*R(2,i),Kss(1,1),n2)
    call vecadd3(tmp,-sq2*R(2,i),Jss(1,3),n2)
    call vecadd3(tmp,-ssq2*R(2,i),Kss(1,3),n2)
    call vecadd3(tmp,-two*R(3,i),Jss(1,4),n2)
    call vecadd3(tmp,-R(3,i),Kss(1,4),n2)
    call vecadd3(tmp,-two,Ja(jJos12,i),n2)
    call vecadd3(tmp,-one,Ka(jJos12,i),n2)
    call vecadd3(tmp,two,Ja(jJvs21,i),n2)
    call vecadd3(tmp,one,Ka(jJvs21,i),n2)
    call vecadd3(tmp,one,Ka(jJov4,i),n2)
    if(i.eq.j) then
      call sasfciszmat(Za(1,m),tmp, R(1,j),NVa,NOb,ia,ib,ir,two)
    else
      call sasfciszmat(Za(1,m),tmp, R(1,j),NVa,NOb,ia,ib,ir,one)
      call vecinit(tmp,n2,0.d0)
      call vecadd3(tmp,sq2*R(2,j),Jss(1,1),n2)
      call vecadd3(tmp,ssq2*R(2,j),Kss(1,1),n2)
      call vecadd3(tmp,-sq2*R(2,j),Jss(1,3),n2)
      call vecadd3(tmp,-ssq2*R(2,j),Kss(1,3),n2)
      call vecadd3(tmp,-two*R(3,j),Jss(1,4),n2)
      call vecadd3(tmp,-R(3,j),Kss(1,4),n2)
      call vecadd3(tmp,-two,Ja(jJos12,j),n2)
      call vecadd3(tmp,-one,Ka(jJos12,j),n2)
      call vecadd3(tmp,two,Ja(jJvs21,j),n2)
      call vecadd3(tmp,one,Ka(jJvs21,j),n2)
      call vecadd3(tmp,one,Ka(jJov4,j),n2)
      call sasfciszmat(Za(1,m),tmp, R(1,i),NVa,NOb,ia,ib,ir,one)
    endif
    endif ! sf_xcis

    m = m + 1
    enddo
  enddo

  return
end subroutine sasfcisz_nac_rhs_offdiag

!======================================================== 
subroutine sasfcisz_nac_rhs_offdiag2(Za,Ja,Ka,Jss,Kss,R,n)
  implicit none
  integer n,i,j,k,m,ia,ib,ic,ir
  real*8  Za(n2,*),Ja(n2*npr,*),Ka(n2*npr,*),Jss(n2,*),Kss(n2,*),R(NOV,*)

  m = 1
  do i=1,nex
    do j=i,nex
    ! s1->s2
    ia = js2
    ib = js1
    ir = 1
    call sasfciszmat(Za(1,m),Kss(1,4),   R(1,i),1,1,ia,ib,ir,one*R(1,j)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,4),   R(1,j),1,1,ia,ib,ir,one*R(1,i)) ! s1->s2 
    call sasfciszmat(Za(1,m),Kss(1,2),   R(1,i),1,1,ia,ib,ir,one*R(2,j)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,2),   R(1,j),1,1,ia,ib,ir,one*R(2,i)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,1),   R(1,i),1,1,ia,ib,ir,sq2*R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,1),   R(1,j),1,1,ia,ib,ir,sq2*R(3,i)) ! s1->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,j),R(1,i),1,1,ia,ib,ir,sq2) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,i),R(1,j),1,1,ia,ib,ir,sq2) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos2,j),R(1,i),1,1,ia,ib,ir,sq2) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJos2,i),R(1,j),1,1,ia,ib,ir,sq2) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJvs1,j),R(1,i),1,1,ia,ib,ir,sq2) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs1,i),R(1,j),1,1,ia,ib,ir,sq2) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,j),R(1,i),1,1,ia,ib,ir,sq2) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,i),R(1,j),1,1,ia,ib,ir,sq2) ! s1->vir
    call sasfciszmat(Za(1,m),Ja(jJov1,j),R(1,i),1,1,ia,ib,ir,sq3) ! occ->vir A
    call sasfciszmat(Za(1,m),Ja(jJov1,i),R(1,j),1,1,ia,ib,ir,sq3) ! occ->vir A
    call sasfciszmat(Za(1,m),Ja(jJov2,i),R(1,j),1,1,ia,ib,ir,one) ! occ->vir B
    call sasfciszmat(Za(1,m),Ja(jJov2,j),R(1,i),1,1,ia,ib,ir,one) ! occ->vir B
    call sasfciszmat(Za(1,m),Ka(jJov2,i),R(1,j),1,1,ia,ib,ir,two) ! occ->vir B
    call sasfciszmat(Za(1,m),Ka(jJov2,j),R(1,i),1,1,ia,ib,ir,two) ! occ->vir B

    ! s2->s1
    ia = js1 
    ib = js2
    ir = 2
    call sasfciszmat(Za(1,m),Kss(1,4),   R(1,i),1,1,ia,ib,ir,one*R(1,j)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,4),   R(1,j),1,1,ia,ib,ir,one*R(1,i)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,2),   R(1,i),1,1,ia,ib,ir,one*R(2,j)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,2),   R(1,j),1,1,ia,ib,ir,one*R(2,i)) ! s1->s2 
    call sasfciszmat(Za(1,m),Kss(1,1),   R(1,i),1,1,ia,ib,ir,sq2*R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,1),   R(1,j),1,1,ia,ib,ir,sq2*R(3,i)) ! s1->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,j),R(1,i),1,1,ia,ib,ir,sq2) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,i),R(1,j),1,1,ia,ib,ir,sq2) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos2,j),R(1,i),1,1,ia,ib,ir,sq2) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJos2,i),R(1,j),1,1,ia,ib,ir,sq2) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJvs1,j),R(1,i),1,1,ia,ib,ir,sq2) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs1,i),R(1,j),1,1,ia,ib,ir,sq2) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,j),R(1,i),1,1,ia,ib,ir,sq2) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,i),R(1,j),1,1,ia,ib,ir,sq2) ! s1->vir
    call sasfciszmat(Za(1,m),Ja(jJov1,j),R(1,i),1,1,ia,ib,ir,-sq3) ! occ->vir A
    call sasfciszmat(Za(1,m),Ja(jJov1,i),R(1,j),1,1,ia,ib,ir,-sq3) ! occ->vir A
    call sasfciszmat(Za(1,m),Ka(jJov1,j),R(1,i),1,1,ia,ib,ir,-sq3) ! occ->vir A
    call sasfciszmat(Za(1,m),Ka(jJov1,i),R(1,j),1,1,ia,ib,ir,-sq3) ! occ->vir A
    call sasfciszmat(Za(1,m),Ja(jJov2,i),R(1,j),1,1,ia,ib,ir,-one) ! occ->vir B
    call sasfciszmat(Za(1,m),Ja(jJov2,j),R(1,i),1,1,ia,ib,ir,-one) ! occ->vir B
    call sasfciszmat(Za(1,m),Ka(jJov2,i),R(1,j),1,1,ia,ib,ir,one) ! occ->vir B
    call sasfciszmat(Za(1,m),Ka(jJov2,j),R(1,i),1,1,ia,ib,ir,one) ! occ->vir B

    ! s1->s1
    ia = js1
    ib = js2
    ir = 3
    call sasfciszmat(Za(1,m),Kss(1,4),   R(1,i),1,1,ia,ia,ir,sq2*R(1,j)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,4),   R(1,j),1,1,ia,ia,ir,sq2*R(1,i)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,2),   R(1,i),1,1,ia,ia,ir,sq2*R(2,j)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,2),   R(1,j),1,1,ia,ia,ir,sq2*R(2,i)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,1),   R(1,i),1,1,ia,ia,ir,R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,1),   R(1,j),1,1,ia,ia,ir,R(3,i)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,3),   R(1,i),1,1,ia,ia,ir,-R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,1),   R(1,i),1,1,ib,ib,ir,-R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,j),R(1,i),1,1,ia,ia,ir,one) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,i),R(1,j),1,1,ia,ia,ir,one) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,j),R(1,i),1,1,ib,ib,ir,-one) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,i),R(1,j),1,1,ib,ib,ir,-one) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos2,j),R(1,i),1,1,ia,ia,ir,one) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJos2,i),R(1,j),1,1,ia,ia,ir,one) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJos2,j),R(1,i),1,1,ib,ib,ir,-one) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJos2,i),R(1,j),1,1,ib,ib,ir,-one) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJvs1,j),R(1,i),1,1,ia,ia,ir,one) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs1,i),R(1,j),1,1,ia,ia,ir,one) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs1,j),R(1,i),1,1,ib,ib,ir,-one) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs1,i),R(1,j),1,1,ib,ib,ir,-one) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,j),R(1,i),1,1,ia,ia,ir,one) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,i),R(1,j),1,1,ia,ia,ir,one) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,j),R(1,i),1,1,ib,ib,ir,-one) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,i),R(1,j),1,1,ib,ib,ir,-one) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJov1,j),R(1,i),1,1,ib,ib,ir,ssq6) ! occ->vir A
    call sasfciszmat(Za(1,m),Ka(jJov1,i),R(1,j),1,1,ib,ib,ir,ssq6) ! occ->vir A
    call sasfciszmat(Za(1,m),Ka(jJov2,j),R(1,i),1,1,ia,ia,ir,sq2) ! occ->vir B
    call sasfciszmat(Za(1,m),Ka(jJov2,i),R(1,j),1,1,ia,ia,ir,sq2) ! occ->vir B
    call sasfciszmat(Za(1,m),Ka(jJov2,j),R(1,i),1,1,ib,ib,ir,-ssq2) ! occ->vir B
    call sasfciszmat(Za(1,m),Ka(jJov2,i),R(1,j),1,1,ib,ib,ir,-ssq2) ! occ->vir B

    ! occ->s1
    ia = js1
    ib = 1
    ic = js2
    ir = os1 
    call sasfciszmat(Za(1,m),Kss(1,4),    R(1,i),1,NOb,ia,ib,ir,sq2*R(1,j)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,4),    R(1,j),1,NOb,ia,ib,ir,sq2*R(1,i)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,2),    R(1,i),1,NOb,ia,ib,ir,sq2*R(2,j)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,2),    R(1,j),1,NOb,ia,ib,ir,sq2*R(2,i)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,1),    R(1,i),1,NOb,ia,ib,ir,one*R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,1),    R(1,j),1,NOb,ia,ib,ir,one*R(3,i)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,3),    R(1,i),1,NOb,ia,ib,ir,-one*R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,3),    R(1,j),1,NOb,ia,ib,ir,-one*R(3,i)) ! s1->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,j), R(1,i),1,NOb,ia,ib,ir,one) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,i), R(1,j),1,NOb,ia,ib,ir,one) ! occ->s1
    call sasfciszmat(Za(1,m),Ja(jJos12,j),R(1,i),1,NOb,ic,ib,ir,one) ! occ->s1
    call sasfciszmat(Za(1,m),Ja(jJos12,i),R(1,j),1,NOb,ic,ib,ir,one) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos2,j), R(1,i),1,NOb,ia,ib,ir,one) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJos2,i), R(1,j),1,NOb,ia,ib,ir,one) ! occ->s2
    call sasfciszmat(Za(1,m),Ja(jJos2,j), R(1,i),1,NOb,ia,ib,ir,-one) ! occ->s2
    call sasfciszmat(Za(1,m),Ja(jJos2,i), R(1,j),1,NOb,ia,ib,ir,-one) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJvs1,j), R(1,i),1,NOb,ia,ib,ir,one) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs1,i), R(1,j),1,NOb,ia,ib,ir,one) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,j), R(1,i),1,NOb,ia,ib,ir,two) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,i), R(1,j),1,NOb,ia,ib,ir,two) ! s2->vir
    call sasfciszmat(Za(1,m),Ja(jJvs2,j), R(1,i),1,NOb,ia,ib,ir,one) ! s2->vir
    call sasfciszmat(Za(1,m),Ja(jJvs2,i), R(1,j),1,NOb,ia,ib,ir,one) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJov1,j), R(1,i),1,NOb,ia,ib,ir,-ssq6) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov1,i), R(1,j),1,NOb,ia,ib,ir,-ssq6) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov1,j), R(1,i),1,NOb,ia,ib,ir,-ssq6) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov1,i), R(1,j),1,NOb,ia,ib,ir,-ssq6) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov2,j), R(1,i),1,NOb,ia,ib,ir,ssq2) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov2,i), R(1,j),1,NOb,ia,ib,ir,ssq2) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov2,j), R(1,i),1,NOb,ia,ib,ir,-ssq2) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov2,i), R(1,j),1,NOb,ia,ib,ir,-ssq2) ! occ->vir

    ! occ->s2
    ia = js2
    ib = 1
    ic = js1
    ir = os2 
    call sasfciszmat(Za(1,m),Kss(1,4),    R(1,i),1,NOb,ia,ib,ir, sq2*R(1,j)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,4),    R(1,j),1,NOb,ia,ib,ir, sq2*R(1,i)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,2),    R(1,i),1,NOb,ia,ib,ir, sq2*R(2,j)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,2),    R(1,j),1,NOb,ia,ib,ir, sq2*R(2,i)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,1),    R(1,i),1,NOb,ia,ib,ir, one*R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,1),    R(1,j),1,NOb,ia,ib,ir, one*R(3,i)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,3),    R(1,i),1,NOb,ia,ib,ir,-one*R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,3),    R(1,j),1,NOb,ia,ib,ir,-one*R(3,i)) ! s1->s1
    call sasfciszmat(Za(1,m),Ka(jJos2,j), R(1,i),1,NOb,ia,ib,ir,one) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJos2,i), R(1,j),1,NOb,ia,ib,ir,one) ! occ->s2
    call sasfciszmat(Za(1,m),Ja(jJos21,j),R(1,i),1,NOb,ic,ib,ir,one) ! occ->s2
    call sasfciszmat(Za(1,m),Ja(jJos21,i),R(1,j),1,NOb,ic,ib,ir,one) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJos1,j), R(1,i),1,NOb,ia,ib,ir,one) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,i), R(1,j),1,NOb,ia,ib,ir,one) ! occ->s1
    call sasfciszmat(Za(1,m),Ja(jJos1,j), R(1,i),1,NOb,ia,ib,ir,-one) ! occ->s1
    call sasfciszmat(Za(1,m),Ja(jJos1,i), R(1,j),1,NOb,ia,ib,ir,-one) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJvs2,j), R(1,i),1,NOb,ia,ib,ir,one) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,i), R(1,j),1,NOb,ia,ib,ir,one) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJvs1,j), R(1,i),1,NOb,ia,ib,ir,two) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs1,i), R(1,j),1,NOb,ia,ib,ir,two) ! s1->vir
    call sasfciszmat(Za(1,m),Ja(jJvs1,j), R(1,i),1,NOb,ia,ib,ir,one) ! s1->vir
    call sasfciszmat(Za(1,m),Ja(jJvs1,i), R(1,j),1,NOb,ia,ib,ir,one) ! s1->vir
    call sasfciszmat(Za(1,m),Ja(jJov1,j), R(1,i),1,NOb,ia,ib,ir,ssq6) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov1,i), R(1,j),1,NOb,ia,ib,ir,ssq6) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov2,j), R(1,i),1,NOb,ia,ib,ir,sq2) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov2,i), R(1,j),1,NOb,ia,ib,ir,sq2) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov2,j), R(1,i),1,NOb,ia,ib,ir,ssq2) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov2,i), R(1,j),1,NOb,ia,ib,ir,ssq2) ! occ->vir

    ! s1->vir
    ia = jv
    ib = js1
    ic = js2
    ir = vs1 
    call sasfciszmat(Za(1,m),Kss(1,4),    R(1,i),NVa,1,ia,ib,ir, sq2*R(1,j)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,4),    R(1,j),NVa,1,ia,ib,ir, sq2*R(1,i)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,2),    R(1,i),NVa,1,ia,ib,ir, sq2*R(2,j)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,2),    R(1,j),NVa,1,ia,ib,ir, sq2*R(2,i)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,1),    R(1,i),NVa,1,ia,ib,ir, one*R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,1),    R(1,j),NVa,1,ia,ib,ir, one*R(3,i)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,3),    R(1,i),NVa,1,ia,ib,ir,-one*R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,3),    R(1,j),NVa,1,ia,ib,ir,-one*R(3,i)) ! s1->s1
    call sasfciszmat(Za(1,m),Ka(jJos2,j), R(1,i),NVa,1,ia,ib,ir,two) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJos2,i), R(1,j),NVa,1,ia,ib,ir,two) ! occ->s2
    call sasfciszmat(Za(1,m),Ja(jJos2,j), R(1,i),NVa,1,ia,ib,ir,one) ! occ->s2
    call sasfciszmat(Za(1,m),Ja(jJos2,i), R(1,j),NVa,1,ia,ib,ir,one) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJos1,j), R(1,i),NVa,1,ia,ib,ir,one) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,i), R(1,j),NVa,1,ia,ib,ir,one) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJvs2,j), R(1,i),NVa,1,ia,ib,ir,one) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,i), R(1,j),NVa,1,ia,ib,ir,one) ! s2->vir
    call sasfciszmat(Za(1,m),Ja(jJvs2,j), R(1,i),NVa,1,ia,ib,ir,-one) ! s2->vir
    call sasfciszmat(Za(1,m),Ja(jJvs2,i), R(1,j),NVa,1,ia,ib,ir,-one) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJvs1,j), R(1,i),NVa,1,ia,ib,ir,one) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs1,i), R(1,j),NVa,1,ia,ib,ir,one) ! s1->vir
    call sasfciszmat(Za(1,m),Ja(jJvs12,j),R(1,i),NVa,1,ia,ic,ir,one) ! s1->vir
    call sasfciszmat(Za(1,m),Ja(jJvs12,i),R(1,j),NVa,1,ia,ic,ir,one) ! s1->vir
    call sasfciszmat(Za(1,m),Ja(jJov1,j), R(1,i),NVa,1,ia,ib,ir,ssq6) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov1,i), R(1,j),NVa,1,ia,ib,ir,ssq6) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov2,j), R(1,i),NVa,1,ia,ib,ir,sq2) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov2,i), R(1,j),NVa,1,ia,ib,ir,sq2) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov2,j), R(1,i),NVa,1,ia,ib,ir,ssq2) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov2,i), R(1,j),NVa,1,ia,ib,ir,ssq2) ! occ->vir
                                              
    ! s2->vir                                
    ia = jv
    ib = js2
    ic = js1
    ir = vs2 
    call sasfciszmat(Za(1,m),Kss(1,4),    R(1,i),NVa,1,ia,ib,ir, sq2*R(1,j)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,4),    R(1,j),NVa,1,ia,ib,ir, sq2*R(1,i)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,2),    R(1,i),NVa,1,ia,ib,ir, sq2*R(2,j)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,2),    R(1,j),NVa,1,ia,ib,ir, sq2*R(2,i)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,1),    R(1,i),NVa,1,ia,ib,ir, one*R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,1),    R(1,j),NVa,1,ia,ib,ir, one*R(3,i)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,3),    R(1,i),NVa,1,ia,ib,ir,-one*R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,3),    R(1,j),NVa,1,ia,ib,ir,-one*R(3,i)) ! s1->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,j), R(1,i),NVa,1,ia,ib,ir,two) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,i), R(1,j),NVa,1,ia,ib,ir,two) ! occ->s1
    call sasfciszmat(Za(1,m),Ja(jJos1,j), R(1,i),NVa,1,ia,ib,ir,one) ! occ->s1
    call sasfciszmat(Za(1,m),Ja(jJos1,i), R(1,j),NVa,1,ia,ib,ir,one) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos2,j), R(1,i),NVa,1,ia,ib,ir,one) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJos2,i), R(1,j),NVa,1,ia,ib,ir,one) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJvs1,j), R(1,i),NVa,1,ia,ib,ir,one) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs1,i), R(1,j),NVa,1,ia,ib,ir,one) ! s1->vir
    call sasfciszmat(Za(1,m),Ja(jJvs1,j), R(1,i),NVa,1,ia,ib,ir,-one) ! s1->vir
    call sasfciszmat(Za(1,m),Ja(jJvs1,i), R(1,j),NVa,1,ia,ib,ir,-one) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,j), R(1,i),NVa,1,ia,ib,ir,one) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,i), R(1,j),NVa,1,ia,ib,ir,one) ! s2->vir
    call sasfciszmat(Za(1,m),Ja(jJvs21,j),R(1,i),NVa,1,ia,ic,ir,one) ! s2->vir
    call sasfciszmat(Za(1,m),Ja(jJvs21,i),R(1,j),NVa,1,ia,ic,ir,one) ! s2->vir
    call sasfciszmat(Za(1,m),Ja(jJov1,j), R(1,i),NVa,1,ia,ib,ir,-ssq6) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov1,i), R(1,j),NVa,1,ia,ib,ir,-ssq6) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov1,j), R(1,i),NVa,1,ia,ib,ir,-ssq6) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov1,i), R(1,j),NVa,1,ia,ib,ir,-ssq6) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov2,j), R(1,i),NVa,1,ia,ib,ir,ssq2) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov2,i), R(1,j),NVa,1,ia,ib,ir,ssq2) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov2,j), R(1,i),NVa,1,ia,ib,ir,-ssq2) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov2,i), R(1,j),NVa,1,ia,ib,ir,-ssq2) ! occ->vir

    ! occ->vir A
    ia = jv
    ib = 1
    ir = ov1 
    call sasfciszmat(Za(1,m),Jss(1,4),    R(1,i),NVa,NOb,ia,ib,ir, sq3*R(1,j)) ! s1->s2
    call sasfciszmat(Za(1,m),Jss(1,4),    R(1,j),NVa,NOb,ia,ib,ir, sq3*R(1,i)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,2),    R(1,i),NVa,NOb,ia,ib,ir,-sq3*R(2,j)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,2),    R(1,j),NVa,NOb,ia,ib,ir,-sq3*R(2,i)) ! s2->s1
    call sasfciszmat(Za(1,m),Jss(1,2),    R(1,i),NVa,NOb,ia,ib,ir,-sq3*R(2,j)) ! s2->s1
    call sasfciszmat(Za(1,m),Jss(1,2),    R(1,j),NVa,NOb,ia,ib,ir,-sq3*R(2,i)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,3),    R(1,i),NVa,NOb,ia,ib,ir,ssq6*R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,3),    R(1,j),NVa,NOb,ia,ib,ir,ssq6*R(3,i)) ! s1->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,j), R(1,i),NVa,NOb,ia,ib,ir,-ssq6) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,i), R(1,j),NVa,NOb,ia,ib,ir,-ssq6) ! occ->s1
    call sasfciszmat(Za(1,m),Ja(jJos1,j), R(1,i),NVa,NOb,ia,ib,ir,-ssq6) ! occ->s1
    call sasfciszmat(Za(1,m),Ja(jJos1,i), R(1,j),NVa,NOb,ia,ib,ir,-ssq6) ! occ->s1
    call sasfciszmat(Za(1,m),Ja(jJos2,j), R(1,i),NVa,NOb,ia,ib,ir,ssq6) ! occ->s2
    call sasfciszmat(Za(1,m),Ja(jJos2,i), R(1,j),NVa,NOb,ia,ib,ir,ssq6) ! occ->s2
    call sasfciszmat(Za(1,m),Ja(jJvs1,j), R(1,i),NVa,NOb,ia,ib,ir,ssq6) ! s1->vir
    call sasfciszmat(Za(1,m),Ja(jJvs1,i), R(1,j),NVa,NOb,ia,ib,ir,ssq6) ! s1->vir
    call sasfciszmat(Za(1,m),Ja(jJvs2,j), R(1,i),NVa,NOb,ia,ib,ir,-ssq6) ! s2->vir
    call sasfciszmat(Za(1,m),Ja(jJvs2,i), R(1,j),NVa,NOb,ia,ib,ir,-ssq6) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,j), R(1,i),NVa,NOb,ia,ib,ir,-ssq6) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,i), R(1,j),NVa,NOb,ia,ib,ir,-ssq6) ! s2->vir
    call sasfciszmat(Za(1,m),Ja(jJov1,j), R(1,i),NVa,NOb,ia,ib,ir,f32) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov1,i), R(1,j),NVa,NOb,ia,ib,ir,f32) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov1,j), R(1,i),NVa,NOb,ia,ib,ir,one) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov1,i), R(1,j),NVa,NOb,ia,ib,ir,one) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov2,j), R(1,i),NVa,NOb,ia,ib,ir,ssq3) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov2,i), R(1,j),NVa,NOb,ia,ib,ir,ssq3) ! occ->vir

    ! occ->vir B
    ia = jv
    ib = 1
    ir = ov2 
    call sasfciszmat(Za(1,m),Jss(1,4),    R(1,i),NVa,NOb,ia,ib,ir, one*R(1,j)) ! s1->s2
    call sasfciszmat(Za(1,m),Jss(1,4),    R(1,j),NVa,NOb,ia,ib,ir, one*R(1,i)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,4),    R(1,i),NVa,NOb,ia,ib,ir, two*R(1,j)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,4),    R(1,j),NVa,NOb,ia,ib,ir, two*R(1,i)) ! s1->s2
    call sasfciszmat(Za(1,m),Kss(1,2),    R(1,i),NVa,NOb,ia,ib,ir, one*R(2,j)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,2),    R(1,j),NVa,NOb,ia,ib,ir, one*R(2,i)) ! s2->s1
    call sasfciszmat(Za(1,m),Jss(1,2),    R(1,i),NVa,NOb,ia,ib,ir,-one*R(2,j)) ! s2->s1
    call sasfciszmat(Za(1,m),Jss(1,2),    R(1,j),NVa,NOb,ia,ib,ir,-one*R(2,i)) ! s2->s1
    call sasfciszmat(Za(1,m),Kss(1,1),    R(1,i),NVa,NOb,ia,ib,ir, sq2*R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,1),    R(1,j),NVa,NOb,ia,ib,ir, sq2*R(3,i)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,3),    R(1,i),NVa,NOb,ia,ib,ir,-ssq2*R(3,j)) ! s1->s1
    call sasfciszmat(Za(1,m),Kss(1,3),    R(1,j),NVa,NOb,ia,ib,ir,-ssq2*R(3,i)) ! s1->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,j), R(1,i),NVa,NOb,ia,ib,ir,ssq2) ! occ->s1
    call sasfciszmat(Za(1,m),Ka(jJos1,i), R(1,j),NVa,NOb,ia,ib,ir,ssq2) ! occ->s1
    call sasfciszmat(Za(1,m),Ja(jJos1,j), R(1,i),NVa,NOb,ia,ib,ir,-ssq2) ! occ->s1
    call sasfciszmat(Za(1,m),Ja(jJos1,i), R(1,j),NVa,NOb,ia,ib,ir,-ssq2) ! occ->s1
    call sasfciszmat(Za(1,m),Ja(jJos2,j), R(1,i),NVa,NOb,ia,ib,ir,ssq2) ! occ->s2
    call sasfciszmat(Za(1,m),Ja(jJos2,i), R(1,j),NVa,NOb,ia,ib,ir,ssq2) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJos2,j), R(1,i),NVa,NOb,ia,ib,ir,sq2) ! occ->s2
    call sasfciszmat(Za(1,m),Ka(jJos2,i), R(1,j),NVa,NOb,ia,ib,ir,sq2) ! occ->s2
    call sasfciszmat(Za(1,m),Ja(jJvs1,j), R(1,i),NVa,NOb,ia,ib,ir,ssq2) ! s1->vir
    call sasfciszmat(Za(1,m),Ja(jJvs1,i), R(1,j),NVa,NOb,ia,ib,ir,ssq2) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs1,j), R(1,i),NVa,NOb,ia,ib,ir,sq2) ! s1->vir
    call sasfciszmat(Za(1,m),Ka(jJvs1,i), R(1,j),NVa,NOb,ia,ib,ir,sq2) ! s1->vir
    call sasfciszmat(Za(1,m),Ja(jJvs2,j), R(1,i),NVa,NOb,ia,ib,ir,-ssq2) ! s2->vir
    call sasfciszmat(Za(1,m),Ja(jJvs2,i), R(1,j),NVa,NOb,ia,ib,ir,-ssq2) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,j), R(1,i),NVa,NOb,ia,ib,ir,ssq2) ! s2->vir
    call sasfciszmat(Za(1,m),Ka(jJvs2,i), R(1,j),NVa,NOb,ia,ib,ir,ssq2) ! s2->vir
    call sasfciszmat(Za(1,m),Ja(jJov2,j), R(1,i),NVa,NOb,ia,ib,ir,f12) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov2,i), R(1,j),NVa,NOb,ia,ib,ir,f12) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov2,j), R(1,i),NVa,NOb,ia,ib,ir,one) ! occ->vir
    call sasfciszmat(Za(1,m),Ka(jJov2,i), R(1,j),NVa,NOb,ia,ib,ir,one) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov1,j), R(1,i),NVa,NOb,ia,ib,ir,ssq3) ! occ->vir
    call sasfciszmat(Za(1,m),Ja(jJov1,i), R(1,j),NVa,NOb,ia,ib,ir,ssq3) ! occ->vir
    m = m + 1
    enddo
  enddo

  return
end subroutine sasfcisz_nac_rhs_offdiag2

!======================================================== 
subroutine sasfciszmat(Za,Fa,Ra,ll,lr,ia,ib,ir,scale)
  implicit none
  integer ll,lr,ia,ib,ir,pr,ps,rp
  real*8  scale
  real*8  Za(*),Fa(*),Ra(*)
  real*8, dimension(:), allocatable :: tmp

  pr = nbas*(ia-1)+1
  ps = nbas*(ib-1)+1
  rp = ia
  allocate(tmp(n2))
  call matmult(Fa(ps),Ra(ir),tmp,nbas,ll,lr,nbas,ll,nbas,3)
  call vecadd3(Za(pr),scale,tmp,nbas*ll)
  call matmult(Fa(rp),Ra(ir),tmp,nbas,lr,ll,nbas,ll,nbas,2)
  call vecadd3(Za(ps),scale,tmp,nbas*lr)
  deallocate(tmp)

  return
end subroutine sasfciszmat

!======================================================== 
subroutine pair2matrix(Z1,Z2,n)
  implicit none
  integer i,j,k,l,n
  real*8  Z1(n2,*),Z2(npair1,*)
  
  do l=1,n
  do i=1,npair1
    j = imo(i)
    k = jmo(i)
    Z1(j+(k-1)*nbas,l) = Z2(i,l)
  enddo
  enddo
  return
end subroutine

!======================================================== 
subroutine mopair(Ea)
  implicit none
  integer i,j
  real*8  tol
  real*8  Ea(*)

  tol = 1.d-4
  npair = 0
  do i=js1,js2
    do j=1,NOa
      if(occ(i).lt.occ(j)) then
        npair = npair + 1
        imo(npair) = i
        jmo(npair) = j
      endif
    enddo
  enddo
  do i=jv,nbas
    do j=1,NOa
      if(occ(i).lt.occ(j)) then
        npair = npair + 1
        imo(npair) = i
        jmo(npair) = j
      endif
    enddo
  enddo
  npair1 = npair
  do i=2,NOa
    do j=1,i-1
      if((Ea(i)-Ea(j)).gt.tol .and. occ(i).eq.occ(j)) then
        npair = npair + 1
        imo(npair) = i
        jmo(npair) = j
      endif
    enddo
  enddo
  !if((Ea(js2)-Ea(js1)).gt.tol) then
  !  npair = npair + 1
  !  imo(npair) = js2
  !  jmo(npair) = js1
  !endif
  do i=NOb+2,nbas
    do j=NOb+1,i-1
      if((Ea(i)-Ea(j)).gt.tol .and. occ(i).eq.occ(j)) then
        npair = npair + 1
        imo(npair) = i
        jmo(npair) = j
      endif
    enddo
  enddo

  if(iprtz.ge.5) then
    write(nb6,*) "MO pairs:" 
    do i=1,npair
      write(nb6,*) i,imo(i),jmo(i)
    enddo
  endif
  return
end subroutine mopair

end module sasfcis_z
