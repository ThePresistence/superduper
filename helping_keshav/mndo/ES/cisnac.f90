!
! Created by Jie Liu on 05/16
! cisnac calculates the nonadiabatic coupling
!
module cis_nonadiabatic_coupling

use es_global
use es_gradient
use cis_gradient
use cis_z
use cis_solve

public cisnac,cisgrad2

contains

!======================================================== 
subroutine cisgrad2(V)
  use LIMIT, only: LM1, LM1M
  implicit none
  real*8 V(*),CG
  common /CGRAD / CG(3,LM1+LM1M)

  call cisnac(V,CG)

  return
end subroutine cisgrad2

!======================================================== 
subroutine cisnac(V,Grad)
  implicit none
  real*8 V(*), Grad(*), time1, time2

  call cis_nac_init()

  call cpusec(time1)
  call cis_nac_z(V)
  call cpusec(time2)
  if(iprtz.gt.1) then
    write(6,'(//A,F10.4//)') "cisz time: ",time2-time1
  endif

  call cpusec(time1)
  call cis_nac_coupling(V,Grad)
  call cpusec(time2)
  if(iprtz.gt.1) then
    write(6,'(//A,F10.4//)') "cis coupling time: ",time2-time1
  endif

  call cis_nac_exit()

  return
end subroutine cisnac

!======================================================== 
subroutine cis_nac_init()

  call es_gradient_setup

  call cis_nac_print_title

  if(is_sf) nden = 2
  n2  = nbas*nbas

  allocate(Zv(NOV,nz))
  call vecinit(Zv,NOV*nz,0d0)
 
  return
end subroutine cis_nac_init

!======================================================== 
subroutine cis_nac_print_title
  if(iprtz.gt.0) then
    write(nb6,*)
    write(nb6,*)     "========================================================"
    if(is_sf) then
      if(do_nac) then
        write(nb6,*) "*       sfcis nonadiabatic coupling calculation        *"
      else
        write(nb6,*) "*            sfcis gradient calculation                *"
      endif
    else
      if(do_nac) then
        write(nb6,*) "*         cis nonadiabatic coupling calculation        *"
      else
        write(nb6,*) "*              cis gradient calculation                *"
      endif
    endif
    write(nb6,*)     "========================================================"
    write(nb6,*)
  endif
end subroutine cis_nac_print_title 

!======================================================== 
subroutine cis_nac_exit()
  if(is_sf) nden = 1
  deallocate(Zv)
end subroutine cis_nac_exit

!======================================================== 
subroutine cis_nac_coupling(V,Grad)
  use limit
  use naccsf
  implicit real*8 (a-h,o-z)
  real*8  V(*), Grad(ncigrd,ncigrd,3,*)
  real*8, dimension(:), allocatable :: de
  real*8, dimension(:,:,:), allocatable :: P,R,P0,G
  common  /CGRAD / CG(3,LM1+LM1M)
  common  /GRDORG/ igrst(lmgrd)

  if(iprtz.gt.0) then
    write(nb6,*)
    write(nb6,*) "     Let's Build CIS Derivatives      "
    write(nb6,*)
  endif

  au2kcal = -1185.80676799574D0
  ev2kcal = 23.06055D0

  call cpusec(time1)
  allocate(P(n2,nz,nden),R(n2,nex,nden),P0(n2,1,nden))
  allocate(G(3,natom,nz-nex),de(nz-nex))
  call vecinit(P, n2*nz*nden,0d0)
  call vecinit(R, n2*nex*nden,0d0)
  call vecinit(P0,n2*nden,0d0)

  call cis_nac_density(P,P(1,1,nden),V)
  if(iprtz.ge.5) then
    write(nb6,*) "cis density"
    call matprnt(P,n2,nz*nden,6)
  endif

  call emake_rao(R,R(1,1,nden),Xv,V(jCa),V(jCb),nex)
  if(iprtz.ge.5) then
    write(nb6,*) "cis transition density"
    call matprnt(R,n2,nex*nden,6)
  endif

  call square(V(jPa),P0,NBas,NBas,LM4)
  if(nden.eq.2) call square(V(jPb),P0(1,1,nden),NBas,NBas,LM4)

  do i=1,natom
    do j=1,3
      call nac_coupling(V,Grad(1,1,j,i),P,R,P0,i,j)
      !if(do_mecp) CG(j,i) = -Grad(igrst(2),igrst(1),j,i)
      CG(j,i) = Grad(1,1,j,i)
    enddo
  enddo
  if(iprtz.gt.0) then
    write(nb6,*) "total gradient"
    call matprnt(Grad,ncigrd*ncigrd,3*natom,6)
  endif
  call cpusec(time2)
  if(iprtz.gt.0) then
    write(6,'(//A,F10.4//)') "nac coupling part1 time: ",time2-time1
  endif

  call cpusec(time1)
  m = 0
  do i=1,ncigrd
    do j=i,ncigrd
      if(i.ne.j) then
        m = m + 1
        do k=1,natom
          G(1,k,m) = -Grad(i,j,1,k)*au2kcal
          G(2,k,m) = -Grad(i,j,2,k)*au2kcal
          G(3,k,m) = -Grad(i,j,3,k)*au2kcal
        enddo 
        if(do_nac) then
          call cpusec(time1)
          de(m) = (Es(j,1) - Es(i,1)) * ev2kcal
          if(imult.le.2 .and. (igrst(i).eq.1.or.igrst(j).eq.1)) then
            if(igrst(i).eq.1) k=j
            if(igrst(j).eq.1) k=i
            call vecadd(P(1,m,1),R(1,k,1),R(1,k,nden),n2)
          else
            call get_cis_density_mn2(P(1,m,1),P(1,m,nden),Xv(1,i),Xv(1+NOVa,i),Xv(1,j),Xv(1+NOVa,j), &
                     V(jCa),V(jCb))
            call vecadd(P(1,m,1),P(1,m,1),P(1,m,nden),n2)
          endif
          call cpusec(time2)
          if(iprtz.gt.0) then
            write(6,'(//A,F10.4//)') "naccsf time: ",time2-time1
          endif
        endif
      endif
    enddo
  enddo
  if(m.gt.0) then
    if(do_nac) call nac_ovlp(P,G,G,de,nbas,natom,m)
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
  call cpusec(time2)
  if(iprtz.gt.0) then
    write(6,'(//A,F10.4//)') "nac coupling part2 time: ",time2-time1
  endif

  if(iprtz.gt.0) then
    write(nb6,*) "total gradient"
    call matprnt(Grad,ncigrd*ncigrd,3*natom,6)
  endif
 
  deallocate(P)
  deallocate(R)
  deallocate(P0)
  deallocate(G)
  return
end subroutine cis_nac_coupling

!======================================================== 
subroutine cis_nac_density(Pa,Pb,V)
  implicit none
  integer i,j,m
  real*8 Pa(n2,*),Pb(n2,*),V(*)

  m = 1
  do i=1,ncigrd
    do j=i,ncigrd
      call get_cis_density_mn(Pa(1,m),Pb(1,m),Xv(1,i),Xv(1+NOVa,i),Xv(1,j),Xv(1+NOVa,j), &
                 Zv(1,m),Zv(1+NOVa,m),V,V,V(jCa),V(jCb),1,1,0)
      m = m + 1
    enddo
  enddo

  return
end subroutine cis_nac_density

!======================================================== 
subroutine nac_coupling(V,CG,P,R,P0,icntr,ic)
  implicit none
  integer icntr,ic,iprtz,i,j,m,ip,ir,istate,jstate
  real*8  EV,BOHR,dstep,scale,dscal,escalp,escalm
  real*8  dot,dot2,dot3,dot4,dot5,dot6
  real*8  V(*),CG(ncigrd,*),P(n2,nz,*),R(n2,nex,*),P0(n2,*)
  real*8, dimension(:,:,:), allocatable :: F21,F22
  real*8, dimension(:,:), allocatable :: F11,F12
  real*8, dimension(:), allocatable :: Pt,H1,H2

  EV   = 27.2100d0
  BOHR = 0.529167d0

  allocate(F11(n2,nden),F12(n2,nden))
  allocate(F21(n2,nex,nden),F22(n2,nex,nden))
  allocate(H1(n2),H2(n2))
  allocate(Pt(n2))

  if(dstep.eq.0.d0) dstep = 2.0d-4
  
  call rstep(dstep,icntr,ic)
  call ee132(V,LM5,P0,F11,1,R,F21,nex,H1,icntr,escalp,nden,i1n)

  call rstep(-2.0*dstep,icntr,ic)
  call ee132(V,LM5,P0,F12,1,R,F22,nex,H2,icntr,escalm,nden,i1n)

  scale = 1/(2*dstep*EV)
  call vecsub(F11,F11,F12,n2*nden)
  call vecsub(F21,F21,F22,nex*n2*nden)
  call vecsub(H1,H1,H2,n2)
  dscal = scale*(escalm-escalp)

  m = 0
  do i=1,ncigrd
    do j=i,ncigrd
      m = m + 1
      if(i.eq.j) then
        call vecadd2(Pt,0.5d0,P0,P(1,m,1),n2)
        call vecdot(dot,F11,Pt,n2)
        call vecadd2(Pt,1.0d0,P0,P(1,m,1),n2)
        call vecdot(dot2,H1,Pt,n2)
      else
        call vecdot(dot,F11,P(1,m,1),n2)
        call vecdot(dot2,H1,P(1,m,1),n2)
      endif
      call vecdot(dot3,F21(1,i,1),R(1,j,1),n2)
      CG(i,j) = - ( dot + dot2 + dot3) *BOHR *scale
      !write(6,*) i,j,CG(i,j),dscal*bohr,dot*bohr*scale,dot2*bohr*scale,dot3*bohr*scale
      if(i.eq.j) CG(i,j) = CG(i,j) + dscal*BOHR
      if(nden.eq.2) then
        if(i.eq.j) then
          call vecadd2(Pt,0.5d0,P0(1,nden),P(1,m,nden),n2)
          call vecdot(dot,F11(1,nden),Pt,n2)
          call vecadd2(Pt,1.0d0,P0(1,nden),P(1,m,nden),n2)
          call vecdot(dot2,H1,Pt,n2)
        else
          call vecdot(dot,F11(1,nden),P(1,m,nden),n2)
          call vecdot(dot2,H1,P(1,m,nden),n2)
        endif
        call vecdot(dot3,F21(1,i,nden),R(1,j,nden),n2)
      endif
      CG(i,j) = CG(i,j) - ( dot + dot2 + dot3) *BOHR *scale
      !write(6,*) i,j,CG(i,j),dot*bohr*scale,dot2*bohr*scale,dot3*bohr*scale
    enddo
  enddo

  call rstep(dstep,icntr,ic)

  deallocate(Pt)
  deallocate(F11)
  deallocate(F12)
  deallocate(F21)
  deallocate(F22)
  deallocate(H1)
  deallocate(H2)

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
subroutine get_cis_density_mn2(Pa,Pb,RIa,RIb,RJa,RJb,Ca,Cb)
  implicit none
  real*8 Pa(*),Pb(*),RIa(*),RIb(*),RJa(*),RJb(*),Ca(LM2,*),Cb(LM2,*)
  real*8, dimension(:), allocatable :: tmp,tmp2,tmp3

  allocate(tmp(n2),tmp2(n2),tmp3(n2))
  call vecinit(Pa,n2,0.d0)

  if(imult.le.2) then
    ! RI * RJ^t
    call matmult(Ca(1,1+NOa),RJa,tmp,NBas,NOa,NVa,LM2,NVa,NBas,1)
    call matmult(Ca(1,1+NOa),RIa,tmp2,NBas,NOa,NVa,LM2,NVa,NBas,1)
    call matmult(tmp,tmp2,Pa,NBas,NBas,NOa,NBas,NBas,NBas,3)
 
    ! RI^t * RJ
    call matmult(Ca(1,1),RJa,tmp,NBas,NVa,NOa,LM2,NVa,NBas,3)
    call matmult(Ca(1,1),RIa,tmp2,NBas,NVa,NOa,LM2,NVa,NBas,3)
    call matmult(tmp,tmp2,tmp3,NBas,NBas,NVa,NBas,NBas,NBas,3)
    call vecadd(Pa,Pa,tmp3,n2)
  else
    call matmult(Ca(1,1),RJa,tmp,NBas,NVb,NOa,LM2,NVb,NBas,3)
    call matmult(Ca(1,1),RIa,tmp2,NBas,NVb,NOa,LM2,NVb,NBas,3)
    call matmult(tmp,tmp2,tmp3,NBas,NBas,NVb,NBas,NBas,NBas,3)
    call vecsub(Pa,Pa,tmp3,n2)
  endif

  if(nden.eq.2) then
    if(imult.le.2) then
      call matmult(Cb(1,1+NOb),RJb,tmp, NBas,NOb,NVb,LM2,NVb,NBas,1)
      call matmult(Cb(1,1+NOb),RIb,tmp2,NBas,NOb,NVb,LM2,NVb,NBas,1)
      call matmult(tmp,tmp2,Pb,NBas,NBas,NOb,NBas,NBas,NBas,3)
     
      call matmult(Cb(1,1),RJb,tmp, NBas,NVb,NOb,LM2,NVb,NBas,3)
      call matmult(Cb(1,1),RIb,tmp2,NBas,NVb,NOb,LM2,NVb,NBas,3)
      call matmult(tmp,tmp2,tmp3,NBas,NBas,NVb,NBas,NBas,NBas,3)
      call vecadd(Pb,Pb,tmp3,n2)
    else
      call matmult(Cb(1,1+NOb),RJa,tmp,NBas,NOa,NVb,LM2,NVb,NBas,1)
      call matmult(Cb(1,1+NOb),RIa,tmp2,NBas,NOa,NVb,LM2,NVb,NBas,1)
      call matmult(tmp,tmp2,Pb,NBas,NBas,NOa,NBas,NBas,NBas,3)
    endif
  endif

  deallocate(tmp)
  deallocate(tmp2)
  deallocate(tmp3)
  return
end subroutine get_cis_density_mn2

!======================================================== 
end module cis_nonadiabatic_coupling
