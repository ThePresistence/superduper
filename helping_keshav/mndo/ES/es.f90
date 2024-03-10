! 
! Created by Jie Liu on 04/16
!

module excited_state

use es_global

integer method

public es0

contains

!======================================================== 
subroutine es0(V,icall)
  implicit none
  integer icall
  real*8  V(*)

  call esinit(V,icall)

  call esrun(V)

  call esexit()

  return 
end subroutine es0

!======================================================== 
subroutine esrun(V)
  implicit none
  real*8  V(*)
 
  if(method.eq.1) then
    call cismain(V)
  else if(method.eq.2) then
    call rpamain(V)
  else
    write(nb6,*) "wrong excited state method!!! please check $kci$"
    stop
  endif

  return
end subroutine esrun 

!======================================================== 
subroutine esinit(V,icall)
  implicit real*8 (a-h,o-z)
  integer icall
  real*8  V(*)
  common  /INOPT2/ in2(300)
  common  /ENERGT/ e0(4)

  ! excited state method
  method = 0
  irpa = .false.
  if(in2(77).eq.6 .or. in2(77).eq.7) method = 1
  if(in2(77).eq.8) then
    method = 2
    irpa = .true.
  endif

  ! job type of excited state methods
  icall1 = icall/10
  icall2 = icall - 10*icall1
  if(icall2.eq.0) jobtype = 0
  if(icall2.ge.1) jobtype = 1
  if(icall2.eq.3) jobtype = 3

  ! save ground state energies
  escf = e0

  ! initialize global setup
  call es_global_setup

  ! symmetrize the MOs
  call es_mo_sym(V)

  ! molecular orbital mapping 
  if(in2(160).gt.0 .and. jobtype.ne.3) then
    gread = .false.
    if(.not.allocated(Ca0)) then
      allocate(Ca0(LM2,LM3))
      call veccopy(Ca0,V(jCa),LM2*LM3)
      if(nden.eq.2) then
        allocate(Cb0(LM2,LM3))
        call veccopy(Cb0,V(jCb),LM2*LM3)
      endif
    else
      call momap(Ca0,Cb0,V(jCa),V(jCb),icall)
    endif
  endif

  ! save X(Y) amplitudes if the state following is requried
  if(stat_fol .and. (.not.allocated(Xn))) then
    allocate(Xn(NOV,nroots))
  endif

  return
end subroutine esinit

!======================================================== 
subroutine esexit()
  use limit, only: lmgrd, lmconf
  implicit real*8 (a-h,o-z)
  common  /INOPT2/ in2(300)
  common  /DYNVA / cicomp(lmconf,6),nciconf
  common  /energt/ e0(4)
  common  /GRDORG/ igrst(lmgrd)

  if(in2(10).ge.1) call esref()

  ! save energies and cis amplitudes for dynamics
  if(jobtype.ne.3) then
    e0    = escf
    e0(1) = escf(1) + es(lstate,1)
    if((in2(160).eq.6).or.((in2(160).eq.1 .or. in2(160).eq.2) &
       .and.(in2(18).eq.6 .or. in2(18).eq.7))) then
      cicomp=0.D0
      do i=1,ngrad
        nciconf = NOV
        cicomp(1:NOV,i) = Xv(1:NOV,i)
      enddo
      ! fixed for cis
      if(.not.is_sf) then
        nciconf = NOV + 1
        do i=1,ngrad
          if(IGRST(i).eq.1) cicomp(NOV+1,IGRST(i)) = 1.d0
        enddo
      endif
    endif
  endif

  if(allocated(Tv)) deallocate(Tv)
  if(allocated(At)) deallocate(At)
  if(allocated(Bt)) deallocate(Bt)
  if(allocated(Es)) deallocate(Es)
  if(allocated(Xv)) deallocate(Xv)
  if(allocated(Yv)) deallocate(Yv)
  if(allocated(jsym)) deallocate(jsym)

  return
end subroutine esexit

end module excited_state
