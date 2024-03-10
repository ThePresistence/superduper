!
!  created by Jie Liu 02/17
!
module es_gradient

use es_global

public

integer :: n2
integer :: nex
integer :: ncigrd
logical :: do_nac
logical :: do_mecp
real*8  :: xx(1)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine es_gradient_setup
  use limit, only : lmgrd
  implicit real*8 (a-h,o-z)
  real*8, dimension(:), allocatable :: tmp, et
  common /INOPT2/ in2(300)
  common /GRDORG/ igrst(lmgrd)

  nrmz    = nrmdav
  ncigrd  = in2(159)
  if(ncigrd.le.0 .and. in2(3).ge.0) ncigrd=1
  nz      = (ncigrd+1)*ncigrd/2
  nex     = ncigrd
  n2      = nbas*nbas
  if(ncigrd.eq.1) igrst(1) = lstate

  do_nac  = .false.
  do_mecp = .false.
  if(in2(160).ge.2 .and. in2(160).le.5 ) do_mecp = .true.
  if(in2(160).ge.6) do_nac = .true.

  if(allocated(Xn)) then
    do i=1,ncigrd
      dmax = 0.d0
      mstate = igrst(i)
      if(.not.is_sf) mstate = igrst(i)-1
      if(mstate.gt.0) then
        do j=1,nroots
          call vecdot(dot,Xn(1,mstate),Xv(1,j),NOV)
          if(abs(dot).gt.dmax) then
            dmax = abs(dot)
            m = j
          endif
        enddo
      else
        m = 1
      endif
      if(m.ne.igrst(i)) then
        write(nb6,*) "the target state is changed from ",igrst(i)," to ",m
      endif
      igrst(i) = m
    enddo
    lstate = igrst(1)
  endif

  allocate(tmp(NOV*nex),et(nex))
  do i=1,nex
    et(i) = Es(igrst(i),1)
    call veccopy(tmp(1+(i-1)*NOV),Xv(1,igrst(i)),NOV)
  enddo
  call veccopy(Es,et,nex)
  call veccopy(Xv,tmp,NOV*nex)
  
  deallocate(tmp)
  deallocate(et)

  return
end subroutine es_gradient_setup

end module es_gradient
