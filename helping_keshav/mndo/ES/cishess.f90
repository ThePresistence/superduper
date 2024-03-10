!
! Created by Jie Liu on 05/16
! 
! cishess perform frequency analysis on the CIS states  
!
module cis_hessian

public :: cishess

contains

subroutine cishess(V)
  ! EE         Excitation energies
  ! X          CIS amplitude
  ! TrialVec   Trial vetors used in Davidson iteration
  ! A          A matrix projected onto the trial vector space 
  Implicit None
  Real*8 V(*)
 
  stop 'cis hessian is not available'

end subroutine cishess

end module cis_hessian
