!
! created by jie liu 12/16
! the main call for cis
!
subroutine cismain(V)
  use limit, only : LM1, LM1M
  use es_global
  use cis_energy
  use sasfcis_energy
  use cis_gradient
  use cis_hessian
  use cis_nonadiabatic_coupling
  use sasfcis_nonadiabatic_coupling
  implicit real*8 (a-h,o-z)
  real*8  V(*)
  common  /CGRAD / CG(3,LM1+LM1M) 
  real*8, dimension(:,:,:,:), allocatable, save :: Grad
 
  if(iprtci.gt.0) then
    write(nb6,'(//5X,A)')  "**************************************"
    write(nb6,'(5X,A)')    "*            CIS CALCULATION         *"
    write(nb6,'(5X,A,//)') "**************************************"
  endif

  ! calculate vertical excitation energies
  if(jobtype.ge.0 .and. jobtype.le.2) then
    call cpusec(time1)
    if(iprtci.gt.1) then
      write(nb6,'(//5X,A,//)') "         VERTICAL RPA EXCITATION ENERGY"
    endif
    if(.not.sasfcis) then
      call cisenergy(V)
    else 
      call sasfcisenergy(V)
    endif
    call cpusec(time2)
    if(iprtci.gt.0) then
      write(6,'(//A,F10.4//)') "excitation energy time: ",time2-time1
    endif
  endif
  
  ! analytical cis gradient of k-th state
  if(jobtype.eq.1 .and. ngrad.eq.1 .and. .not.ifd) then
    call cpusec(time1)
    if(iprtci.gt.1) then
      write(nb6,'(//5X,A,//)') "             ANALYTICAL CIS GRADIENT"
    endif
    if(.not.sasfcis) then
      call cisgrad2(V)
    else
      call sasfcisgrad(V)
    endif
    call cpusec(time2)
    if(iprtci.gt.0) then
      write(6,'(//A,F10.4//)') "analytic gradient time: ",time2-time1
    endif
  endif

  ! analytical cis gradients or nonadiabatic couplings
  if(jobtype.eq.1 .and. ngrad .gt.1 .and. .not.ifd) then
    call cpusec(time1)
    if(iprtci.gt.1) then
      write(nb6,'(//5X,A,//)') "             ANALYTICAL CIS NONADIABATIC COUPLING"
    endif
    if(.not.allocated(Grad)) allocate(Grad(ngrad,ngrad,3,natom))
    if(.not.sasfcis) then
      call cisnac(V,Grad)
    else
      call sasfcisnac(V,Grad)
    endif
    call cpusec(time2)
    if(iprtci.gt.0) then
      write(6,'(//A,F10.4//)') "analytic nonadiabatic coupling time: ",time2-time1
    endif
  endif

  ! analytical cis hessian
  if(jobtype.eq.2) then
    if(iprtci.gt.1) then
      write(nb6,'(//5X,A,//)') "             ANALYTICAL CIS HESSIAN IS NOT AVAILABLE"
      stop
    endif
    call cishess(V)
  endif

  ! read cis gradients or nonadiabatic couplings
  if(jobtype.eq.3) then
    if(nstate.eq.0) then
      ! read cis gradients
      if(iprtci.gt.1) then
        write(nb6,'(//5X,A,//)') "             READ ANALYTICAL CIS GRADIENT"
      endif
      if(iprtci.gt.0) write(nb6,*) "READ GRADIENT ",mstate
      do i=1,natom
        do j=1,3
          CG(j,i) = Grad(mstate,mstate,j,i)
        enddo
      enddo
    else
      ! read cis nonadiabatic couplings
      if(iprtci.gt.1) then
        write(nb6,'(//5X,A,//)') "             READ ANALYTICAL CIS NONADIABATIC COUPLING"
      endif
      if(iprtci.gt.0) write(nb6,*) "READ COUPLING ", mstate, nstate
      do i=1,3*natom
        do j=1,3
          CG(j,i) = Grad(mstate,nstate,j,i)
        enddo
      enddo
    endif
  endif

  return
end subroutine cismain
