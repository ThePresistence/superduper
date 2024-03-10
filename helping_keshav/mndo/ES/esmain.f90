!
! the excited state interface to the main program
!

!======================================================== 
! esmain is used by scfcal, ...
!======================================================== 
subroutine esmain(V,icall)
  use excited_state
  implicit none
  integer icall
  real*8  V(*)

  call es0(V,icall)

  return
end subroutine esmain

!========================================================
!  initialization for excited state calculations
!========================================================
subroutine ines(iform,jprint)
  implicit real*8 (a-h,o-z)
  integer i, iform, jprint
  integer NBF,IN1,IN2,IMOPAC
  logical UHF
  common  /NBFILE/ NBF(20)
  common  /INOPT1/ IN1(300)
  common  /INOPT2/ IN2(300)  
  common  /MOPAC / IMOPAC
  common  /ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
  common  /UHF   / UHF

  ! read excited state arguments 
  if(imopac.eq.0) then
     if(iform.le.0) then
        read(NBF(5),'(20I4)') (IN1(i),i=131,150)
        read(NBF(5),'(20I4)') (IN1(i),i=151,170)
     else
        read(NBF(5),*)   (IN1(i),i=131,150)
        read(NBF(5),*)   (IN1(i),i=151,170)
     endif
  endif

  ! initialize the excited state method
  IN2(131:170) = IN1(131:170)

  ! the number of roots and the target state
  if(IN2(139).eq.0) IN2(139) = 6
  if(IN2(140).eq.0) IN2(140) = 1

  ! the number of molecular orbitals
  NMOS = NORBS

  ! active space
  if(IN2(131).le.0) then
    IN2(131) = NALPHA
  endif

  if(IN2(132).le.0) then
    IN2(132) = NORBS - NALPHA
  endif

  ! excitation level
  IN2(138) = 1
  
  ! multiplicity of the excited state
  if(IN2(142).eq.0) then
    IN2(142) = IN2(66)
    if(IN2(142).eq.0) IN2(142) = 1
  endif

  ! SF-CIS uses the open-shell HF reference
  UHF = .false.
  if(IN2(77).eq.6 .and. IN2(66).eq.3) UHF = .true.

  ! we will use triplet reference since there is no difference between
  ! the singlet and triplet reference
  if(IN2(66).eq.1) IN2(66) = 3

  return
end subroutine
