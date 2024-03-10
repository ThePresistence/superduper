!
! Created by Jie Liu on 05/16
! sasfcisnac calculates the nonadiabatic coupling for spin adapted spin-flipped cis
!
module sasfcis_nonadiabatic_coupling

use es_global
use cis_solve
use sasfcis_util
use sasfcis_z

private

integer ncigrd,iprint,ncigrdsave
real*8  xx(1)
integer, dimension(:), allocatable :: list

public sasfcisnac, sasfcisgrad

contains

subroutine sasfcisgrad(V)
  use LIMIT, only: LM1, LM1M
  implicit none
  real*8 V(*),CG
  common /CGRAD / CG(3,LM1+LM1M)

  call sasfcisnac(V,CG)

  return
end subroutine sasfcisgrad

!======================================================== 
subroutine sasfcisnac(V,Grad)
  implicit none
  real*8 V(*), Grad(*)

  call sasfcis_nac_init()

  call sasfcis_nac_z(V)

  call sasfcis_nac_coupling(V,Grad)

  call sasfcis_nac_exit()

  return
end subroutine sasfcisnac

!======================================================== 
subroutine sasfcis_nac_init()
  USE LIMIT, ONLY: LMGRD, LEN
  implicit none
  integer IGRST, in2, i, j, k, m
  real*8  dot, maxs
  common  /GRDORG/ IGRST(LMGRD)
  common  /INOPT2/ in2(300)

  call sasfcis_setup  

  iprint = in2(41)
  ncigrd = in2(159)
  if(iprint.gt.2) then
    if(ncigrd.gt.0) then
      write(nb6,*)
      write(nb6,*) "========================================================"
      write(nb6,*) "*        sasfcis nonadiabatic coupling calculation     *"
      write(nb6,*) "========================================================"
      write(nb6,*)
    else
      write(nb6,*)
      write(nb6,*) "========================================================"
      write(nb6,*) "*               sasfcis gradient calculation           *"
      write(nb6,*) "========================================================"
      write(nb6,*)
    endif
  endif

  do_nac = .false.
  do_mecp = .false.
  if(in2(160).ge.2) do_nac = .true.
  if(in2(160).ge.3 .and. in2(160).le.5) do_mecp = .true.

  if(.not.do_nac) then
  if(ncigrd.le.1) then
    maxs = 0.d0
    if(allocated(Xn)) then
      do i=1,nroots
        !call vecdot(dot,Xn(1,1),Xv(1,i),NOV)
        dot = 0.d0
        do j=1,NOV
          if(abs(Xn(j,lstate)*Xv(j,i)).gt.dot) dot = abs(Xn(j,lstate)*Xv(j,i))
        enddo
        if(abs(dot).gt.maxs) then
          maxs = abs(dot)
          m = i
        endif
      enddo
      if(m.ne.lstate) then
        write(nb6,*) "the target state is changed from ",lstate," to ",m
      endif
      in2(140) = m
      lstate = m
    endif
  else
    if(allocated(Xn)) then
      do i=1,ncigrd
        maxs = 0.d0
        do j=1,nroots
          !call vecdot(dot,Xn(1,i),Xv(1,j),NOV)
          dot = 0.d0
          do k=1,NOV
            if(abs(Xn(k,i)*Xv(k,j)).gt.dot) dot = abs(Xn(k,i)*Xv(k,j))
          enddo
          if(abs(dot).gt.maxs) then
            maxs = abs(dot)
            m = i
          endif
        enddo
        if(m.ne.IGRST(i)) then
          write(nb6,*) "the target state is changed from ",IGRST(i)," to ",m
        endif
        IGRST(i) = m
      enddo
      lstate = IGRST(1)
    endif
  endif
  endif

  ncigrdsave = ncigrd
  if(ncigrd.gt.1) then
    allocate(list(ncigrd))
    do i=1,ncigrd
      list(i) = IGRST(i)
    enddo
  else if(ncigrd.eq.1) then
    allocate(list(ncigrd))
    list(1) = in2(140)
  else  
    write(nb6,*) "wrong ncigrd: ", ncigrd
    stop
  endif

  nz  = (ncigrd+1)*ncigrd/2
  nex = ncigrd
  n2  = nbas*nbas
  if(iprint.ge.5) then
    write(nb6,*) "excited state gradient:",nex
    write(nb6,*) "nonadiabatic coupling derivatives:",nz-nex
  endif

  allocate(Zv(n2,nz))
  ! rearrage the excitation energies and transition coefficients
  do i=1,nex
    call veccopy(Zv(1+(i-1)*NOV,1),Xv(1,list(i)),NOV)
    call veccopy(Xv(1,list(i)),Xv(1,i),NOV)
  enddo
  call veccopy(Xv,Zv,NOV*nex)
  if(iprint.ge.5) then
    write(nb6,*) "rearraged amplitude"
    call matprnt(Xv,NOV,nex,6)
  endif
  do i=1,nex
    Zv(i,1) = Es(list(i),1)
  enddo
  call veccopy(Es,Zv,nex)
  call vecinit(Zv,n2*nz,0d0)

  return
end subroutine sasfcis_nac_init

!======================================================== 
subroutine sasfcis_nac_exit()
  ncigrd = ncigrdsave
  ! clean up the memory
  if(allocated(Js))   deallocate(Js)
  if(allocated(Ks))   deallocate(Ks)
  if(allocated(list)) deallocate(list)
  if(allocated(Zv))   deallocate(Zv)
  return
end subroutine sasfcis_nac_exit

!======================================================== 
subroutine sasfcis_nac_coupling(V,Grad)
  use limit
  use naccsf
  implicit none
  integer i,j,m,k,l,np,nd
  integer jP0,jPcc,jR,jP,jPs1,jPs2,jPs12
  real*8  au2kcal, ev2kcal, CG
  real*8  V(*), Grad(ncigrd,ncigrd,3,*)
  real*8, dimension(:), allocatable :: P,de
  real*8, dimension(:,:,:), allocatable :: G
  common  /CGRAD / CG(3,LM1+LM1M)

  if(iprint.ge.5) then
    write(nb6,*)
    write(nb6,*) "====================Let's Build CIS Derivatives======================="
    write(nb6,*)
  endif

  au2kcal = -1185.80676799574D0
  ev2kcal = 23.06055D0

  nd = 3
  if(sref) nd = 4
  if(sf_xcis) nd = 7
  jP0  = 1                     ! ground state density (1)
  jR   = jP0 + n2              ! transition densities (npr*ncigrd)
  jPcc = jR  + n2*npr*ncigrd   ! C_si*C_sj^t          (3)
  jP   = jPcc+ n2*3            ! excited state densities (difference) (nd)
  jPs1 = jP  + n2*nz           ! s1s1
  jPs2 = jP  + n2*nz*2         ! s2s2
  np   = 4+npr*nex+nd*nz

  allocate(de(nz-nex))
  allocate(P(n2*np))
  allocate(G(3,natom,nz-nex))
  call vecinit(P,n2*np,0d0)

  ! prepare density matrices
  ! 1. ground state density
  call square(V(jPa),P(jP0),nbas,nbas,LM4)
  call vecscale(P(jP0),n2,two)

  ! 2. C_s1*C_s1^t C_s1*C_s2^t and C_s2*C_s2^t
  call sasfemake_sao(P(jPcc),V(jCa))

  ! 3. transition densities
  call sasfemake_rao(P(jR),Xv,V(jCa),nex)
  if(iprint.ge.5) then
    write(nb6,*) "************sasfcis transition density**********"
    do i=1,nex*npr
      write(nb6,*) "sasfcis transition density",i
      call matprnt(P(jR+(i-1)*n2),nbas,nbas,6)
    enddo
  endif

  ! 4. excited state densities
  call sasfemake_pao(P(jP),Xv,V(jCa),nd,nz,1,1)
  if(iprint.ge.5) then
    write(nb6,*) "****************sasfcis density*****************"
    do i=1,nz
      write(nb6,*) "sasfcis density",i
      call matprnt(P(jP+(i-1)*n2),nbas,nbas,6)
    enddo
  endif

  ! CI and MO part
  do i=1,natom
    do j=1,3
      call nac_coupling(V,Grad(1,1,j,i),P,i,j)
      !if(do_mecp) CG(j,i) = -Grad(list(2),list(1),j,i)
      CG(j,i) = Grad(1,1,j,i)
    enddo
  enddo

  ! AO part (only existing in nonadiabatic couplings)
  if(do_nac) then
    call sasfemake_pao_nonsym(P(jP),Xv,V(jCa),nz)
  endif
  m = 0
  l = 0
  do i=1,ncigrd
    do j=i,ncigrd
      l = l + 1
      if(i.ne.j) then
        m = m + 1
        do k=1,natom
          G(1,k,m) = Grad(i,j,1,k)*au2kcal
          G(2,k,m) = Grad(i,j,2,k)*au2kcal
          G(3,k,m) = Grad(i,j,3,k)*au2kcal
        enddo 
        if(do_nac) then
          de(m) = (Es(j,1) - Es(i,1)) * ev2kcal
          call veccopy(P(jP+(m-1)*n2),P(jP+(l-1)*n2),n2)
          if(iprint.ge.6) then
            write(nb6,*) "gradient x"
            call matprnt(G,3,natom,6)
          endif
        endif
      endif
    enddo
  enddo

  if(m.gt.0) then
    if(do_nac) call nac_ovlp(P(jP),G,G,de,nbas,natom,m)
    m = 0
    do i=1,ncigrd
      do j=i,ncigrd
        if(i.ne.j) then
          m = m + 1
          do k=1,natom
           Grad(i,j,1,k) =  G(1,k,m)/au2kcal
           Grad(i,j,2,k) =  G(2,k,m)/au2kcal
           Grad(i,j,3,k) =  G(3,k,m)/au2kcal
           Grad(j,i,1,k) = - Grad(i,j,1,k)
           Grad(j,i,2,k) = - Grad(i,j,2,k)
           Grad(j,i,3,k) = - Grad(i,j,3,k)
          enddo
        endif
      enddo
    enddo
  endif


  if(iprint.ge.5) then
    write(nb6,*) "total gradient"
    call matprnt(Grad,ncigrd*ncigrd,3*natom,8)
  endif
 
  deallocate(de)
  deallocate(P)
  deallocate(G)
  return
end subroutine sasfcis_nac_coupling

!======================================================== 
subroutine nac_coupling(V,CG,P,icntr,ic)
  implicit none
  integer icntr,ic,i,j,m,np
  integer jJ0,jJ11,jJ12,jJ22,jJ21,jJr,jJ
  integer jP0,jP11,jP12,jP22,jR,jP,jPs1,jPs2,jPs3,jPs4,jPs5,jPs6
  real*8  EV,BOHR,dstep,scale,dscal,escalp
  real*8  dot,dot2,dot3,dot4,dot5,dot6,dot7,dot8,dot9,dot10,dot11
  real*8  V(*),CG(ncigrd,*),P(*)
  real*8, dimension(:,:), allocatable :: Za
  real*8, dimension(:), allocatable :: Ja,Ka,H1,Ps,tmp

  EV = 27.2100d0
  BOHR = 0.529167d0

  jP0  = 1                     ! ground state density
  jR   = jP0 + n2              ! transition densities
  jP11 = jR  + n2*npr*ncigrd   ! C_si*C_sj^t
  jP12 = jP11+ n2
  jP22 = jP11+ n2*2
  jP   = jP11+ n2*3            ! excited state densities (difference)
  jPs1 = jP  + n2*nz           ! Ks1s1
  jPs2 = jP  + n2*nz*2         ! Ks2s2
  jPs3 = jP  + n2*nz*3         ! Js1s2
  jPs4 = jP  + n2*nz*4         ! Js2s2
  jPs5 = jP  + n2*nz*5         ! Ks1s2
  jPs6 = jP  + n2*nz*6         ! Js1s2
  jJ0  = 1
  jJr  = jJ0 + n2
  jJ11 = jJr + n2*npr*ncigrd
  jJ12 = jJ11+ n2
  jJ22 = jJ11+ n2*2
  jJ21 = jJ11+ n2*3
  np   = 5 + npr*ncigrd

  allocate(Ja(n2*np),Ka(n2*np),Za(NOV,ncigrd))
  allocate(H1(n2),Ps(n2))
  allocate(tmp(n2))

  if(dstep.eq.0.d0) dstep = 2.0d-4
  scale = 1/(2*dstep*EV)
  
  call rstep(dstep,icntr,ic)
  call makejkx(V,LM5,P,P,Ja,Ja,Ka,Ka,np-1,H1,icntr,escalp,1)
  call mattrans(Ja(jJ21),Ja(jJ12),nbas,nbas)
  call mattrans(Ka(jJ21),Ka(jJ12),nbas,nbas)
  dscal = scale*escalp
  call sasfemake_a2e2(Za,Ja(jJr),Ka(jJr),Ja(jJ11),Ka(jJ11),Xv,V(jCa),ncigrd)
  call vecadd3(Ja(jJ0),f12,Ka(jJ0),n2)
  m = 0
  do i=1,ncigrd
  do j=i,ncigrd
    if(i.eq.j) then
    ! (P+P0)*h^x
      call vecadd(Ps,P(jP+m*n2),P(jP0),n2)
      call vecdot(dot,Ps,H1,n2)
    ! (P+1/2*P0)*II^x
      call vecadd2(Ps,f12,P(jP0),P(jP+m*n2),n2)
      call vecdot(dot2,Ps,Ja(jJ0),n2)
    else 
    ! P*h^x
      call vecdot(dot,P(jP+m*n2),H1,n2)
    ! P*II^x
      call vecdot(dot2,P(jP+m*n2),Ja(jJ0),n2)
    endif
    ! R*A[R]^[x]
    call vecdot(dot3,Za(1,i),Xv(1,j),NOV)
    ! Ps1*(\mu s1|\nu s1)^[x]
    call vecdot(dot4,Ka(jJ11),P(jPs1+m*n2),n2)
    ! Ps2*(\mu s2|\nu s2)^[x]
    call vecdot(dot5,Ka(jJ22),P(jPs2+m*n2),n2)
    if(sf_xcis) then
      call vecdot(dot6,Ja(jJ11),P(jPs3+m*n2),n2)
      call vecdot(dot7,Ja(jJ22),P(jPs4+m*n2),n2)
      call vecdot(dot8,Ka(jJ12),P(jPs5+m*n2),n2)
      call vecdot(dot9,Ja(jJ12),P(jPs6+m*n2),n2)
    else if(sref) then
      ! P2*(\mu \nu|s1 s1)^[x]
      call vecdot(dot6,Ja(jJ11),P(jPs3+m*n2),n2)
      call vecdot(dot7,Ja(jJ22),P(jPs3+m*n2),n2)
      dot6 = -dot6
      dot8 = 0.d0
      dot9 = 0.d0
    else
      dot6 = 0.d0
      dot7 = 0.d0
      dot8 = 0.d0
      dot9 = 0.d0
    endif

    CG(i,j) = (-dot-dot2-dot3+dot4+dot5-dot6-dot7+dot8-dot9)*scale*BOHR
    if(i.eq.j) CG(i,i) = CG(i,i) - dscal*BOHR
    !if(list(i).eq.2 .and. list(j).eq.2) write(nb6,*) "step1:",icntr,ic,CG(i,j),dscal*BOHR,dot,dot2,dot3,dot4,dot5
    m = m + 1
  enddo ! j
  enddo ! i

  call rstep(-2.0*dstep,icntr,ic)
  call makejkx(V,LM5,P,P,Ja,Ja,Ka,Ka,np-1,H1,icntr,escalp,1)
  call mattrans(Ja(jJ21),Ja(jJ12),nbas,nbas)
  call mattrans(Ka(jJ21),Ka(jJ12),nbas,nbas)
  dscal = scale*escalp
  call sasfemake_a2e2(Za,Ja(jJr),Ka(jJr),Ja(jJ11),Ka(jJ11),Xv,V(jCa),ncigrd)
  call vecadd3(Ja(jJ0),f12,Ka(jJ0),n2)
  m = 0
  do i=1,ncigrd
  do j=i,ncigrd
    if(i.eq.j) then
    ! (P+P0)*h^x
      call vecadd(Ps,P(jP+m*n2),P(jP0),n2)
      call vecdot(dot,Ps,H1,n2)
    ! (P+1/2*P0)*II^x
      call vecadd2(Ps,f12,P(jP0),P(jP+m*n2),n2)
      call vecdot(dot2,Ps,Ja(jJ0),n2)
    else 
    ! P*h^x
      call vecdot(dot,P(jP+m*n2),H1,n2)
    ! P*II^x
      call vecdot(dot2,P(jP+m*n2),Ja(jJ0),n2)
    endif
    ! R*A[R]^[x]
    call vecdot(dot3,Za(1,i),Xv(1,j),NOV)
    ! Ps1*(\mu s1|\nu s1)^[x]
    call vecdot(dot4,Ka(jJ11),P(jPs1+m*n2),n2)
    ! Ps2*(\mu s2|\nu s2)^[x]
    call vecdot(dot5,Ka(jJ22),P(jPs2+m*n2),n2)
    if(sf_xcis) then
      call vecdot(dot6,Ja(jJ11),P(jPs3+m*n2),n2)
      call vecdot(dot7,Ja(jJ22),P(jPs4+m*n2),n2)
      call vecdot(dot8,Ka(jJ12),P(jPs5+m*n2),n2)
      call vecdot(dot9,Ja(jJ12),P(jPs6+m*n2),n2)
    else if(sref) then
      ! P2*(\mu \nu|s1 s1)^[x]
      call vecdot(dot6,Ja(jJ11),P(jPs3+m*n2),n2)
      call vecdot(dot7,Ja(jJ22),P(jPs3+m*n2),n2)
      dot6 = -dot6
      dot8 = 0.d0
      dot9 = 0.d0
    else
      dot6 = 0.d0
      dot7 = 0.d0
      dot8 = 0.d0
      dot9 = 0.d0
    endif

    CG(i,j) = CG(i,j) - (-dot-dot2-dot3+dot4+dot5-dot6-dot7+dot8-dot9)*scale*BOHR
    if(i.eq.j) CG(i,i) = CG(i,j) + dscal*BOHR
    !if(list(i).eq.2 .and. list(j).eq.2) write(nb6,*) "step2:",icntr,ic,CG(i,j),dscal*BOHR,dot,dot2,dot3,dot4,dot5
    !CG(j,i) = CG(i,j)
    m = m + 1
  enddo
  enddo

  call rstep(dstep,icntr,ic)

  deallocate(Ja)
  deallocate(Ka)
  deallocate(Za)
  deallocate(H1)
  deallocate(Ps)
  return
end subroutine nac_coupling

!======================================================== 
subroutine rstep(dx,icntr,ic)
  use limit, only: LM1
  implicit none
  real*8 dx, coord
  integer icntr,ic
  common /atomc/ coord(3,LM1)

  coord(ic,icntr) = coord(ic,icntr) + dx

end subroutine rstep
!======================================================== 

end module sasfcis_nonadiabatic_coupling
