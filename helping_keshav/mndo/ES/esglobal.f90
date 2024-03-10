!======================================================== 
!
!  created by Jie Liu on 03/16
!
!  definition of excited state global variables
!
!======================================================== 

module es_global

public

! excited state job type
integer :: jobtype

! es control
integer :: i1        ! singlet
integer :: i3        ! triplet
integer :: i13       ! singlet and triplet
integer :: i1n       ! do singlet es calculation
integer :: imult     ! multiplicity of the ground state
integer :: nden      ! 1: closed shell, 2: open shell
integer :: nroots    ! number of exited states
integer :: iprtci    ! print out level
integer :: lstate    ! target excited state      
integer :: iprtz     ! print out details of gradient or NAC
logical :: sasfcis   ! spin adapted spin-flipped es
integer :: iuvcd     ! control the output of properties
logical :: gread     ! if reading the initial guess
integer :: nciref    ! number of references
logical :: sf_xcis   ! spin-flip XCIS
logical :: stat_fol  ! tracking the excited state
logical :: is_sf     ! spin-flip
logical :: ifd       ! the finite-difference method for derivatives
logical :: irpa      ! TDHF

! es dimension
integer :: NOa       ! number of alpha occupied orbital
integer :: NOb       ! number of beta occupied orbital
integer :: NVa       ! number of alpha virtual orbital
integer :: NVb       ! number of beta virtual orbital
integer :: NOVa      ! alpha OV space dimension
integer :: NOVb      ! beta OV space dimension
integer :: NOV       ! OV space dimension
integer :: NOrb      ! number of molecular orbitals
integer :: NAlpha    ! number of alpha electrons
integer :: NBeta     ! number of beta electrons
integer :: NBas      ! number of AO orbitals

! es energy (gradient, nac, ...) iteration
integer :: NL        ! number of unconverged states
integer :: NS        ! number of trial vectors
integer :: NZ        ! number of Z vectors
integer :: KitDav    ! maximum number of davidson iteration
integer :: NIter     ! current iteration
integer :: MaxDav    ! maximum dimension of davidson subspace
integer :: NrmDav    ! davidson convergence tolerance 
integer :: NrmZ      ! davidson convergence tolerance for Z vector

! es dynamic interface
integer :: ngrad
integer :: mstate
integer :: nstate

! ground state matirx pointer 
integer :: jEa       ! alpha orbital energies
integer :: jEb       ! beta orbital energies
integer :: jCa       ! alpha molecular coefficient
integer :: jCb       ! beta molecular coefficient
integer :: jFa       ! alpha fock matrix
integer :: jFb       ! beta fock matrix
integer :: jPa       ! alpha density matrix
integer :: jPb       ! beta density matrix
integer :: jW        ! 2e integral
integer :: jEa0      ! alpha orbital energies
integer :: jEb0      ! beta orbital energies
integer :: jCa0      ! alpha molecular coefficient
integer :: jCb0      ! beta molecular coefficient


! matrix dimension (see definition in DYNDEF)
integer :: LM2        
integer :: LM3
integer :: LM4
integer :: LM5
integer :: LM6

! parametrization 
real*8  :: cFac1     ! scaling of J and K matrices
real*8  :: cFacJ     ! scaling of J matrix
real*8  :: cFacK     ! scaling of K matrix
integer :: inrefd    ! reference data (see definition of input)

! output file
integer :: NB6       ! output file 

! number of atoms
integer :: natom     

! save energy
real*8  :: etarget   ! saved tagert state energy                               
real*8  :: escf(4)   ! saved ground state energy                               

! symmetry control
integer                             :: ncisym       ! number of symmetry  
logical                             :: is_sym       ! symmetry
character*4                         :: labels(8)    ! symmetry lables
integer, dimension(:), allocatable  :: jsym         ! state symmetry
real*8, dimension(:), allocatable   :: fv, dv       ! save oscillation strength and transition dipole moments

! allocate memory
real*8, dimension(:,:), allocatable :: Tv           ! trial vectors
real*8, dimension(:,:), allocatable :: At           ! A matrix projected on davidson subspace
real*8, dimension(:,:), allocatable :: Bt           ! excitation energies
real*8, dimension(:,:), allocatable :: Es           ! excitation energies 
real*8, dimension(:,:), allocatable :: Xv           ! eigenvector of excited states
real*8, dimension(:,:), allocatable :: Yv           ! eigenvector of excited states
real*8, dimension(:,:), allocatable :: Zr           ! r.h.s. of Z vector equation
real*8, dimension(:,:), allocatable :: Zv           ! Z vectors
real*8, dimension(:,:), allocatable, save :: Xn     ! saved old eigenvector of excited states
real*8, dimension(:,:), allocatable, save :: Ca0,Cb0

public es_global_setup

contains

!======================================================== 
! es global initialization
!======================================================== 
subroutine es_global_setup
  use limit, only : len,lmgrd
  implicit real*8 (a-h,o-z)
  real*8  c(3)
  common  /LMSCF / LS(9)
  common  /LMUHF / LU(8)
  common  /LIMITS/ LM(7)
  common  /NBFILE/ NBF(20)
  common  /INOPT2/ IN2(300)
  common  /ORBITS/ NUMB,NORBS,NMOS,NA,NB
  common  /ATOMS / NUMAT
  common  /GRDORG/ IGRST(LMGRD),ISTATE,JSTATE
  common  /SKIPA / ISKPA
 
  ! the output file
  nb6 = NBF(6)

  ! the number of atoms
  natom = NUMAT

  ! symmetry of excited states
  ncisym = IN2(143)

  ! es arguments
  call es_control

  ! es dimension
  call ovdim 

  ! dynamics interface
  ngrad = IN2(159)
  if(ngrad.gt.0) then
    do i=1,ngrad
      if(IGRST(i).eq.ISTATE) mstate = i
      if(IGRST(i).eq.JSTATE) nstate = i
    enddo
    if(JSTATE.eq.0) nstate = 0
  endif

  ! tracking the excited state when running optimization or surface hopping
  stat_fol = .false.
  if(jobtype.ge.1) stat_fol = .true.

  ! if the finite difference method is required 
  ifd = .false.
  if(iskpa.eq.1) ifd = .true.

  ! ground state matrix pointers
  jCa = LS(1) + (NAlpha-NOa)*LM(1)
  jEa = LS(2) + NAlpha - NOa
  jFa = LS(6)
  jPa = LS(8)
  jW  = LS(9)
  if(nden.eq.2) then
    jCb = LU(1) + (NBeta -NOb)*LM(1)
    jEb = LU(2) + NBeta - NOb
    jFb = LU(6)
    jPb = LU(8)
  else
    jCb = jCa
    jEb = jEa
    jFb = jFa
    jPb = jPa
  endif

  jCa0 = LS(1) 
  jEa0 = LS(2)
  if(nden.eq.2) then
    jCb0 = LU(1) + (NBeta -NOb)*LM(1)
    jEb0 = LU(2) + NBeta - NOb
  else
    jCb0 = jCa0
    jEb0 = jEa0
  endif

  ! matrix dimension
  LM2 = LM(1)
  LM3 = LM(2)
  LM4 = LM(3)
  LM5 = LEN
  LM6 = LM(4)

  ! Davidson iteration
  nl     = 0
  ns     = 0
  niter  = 0
  kitdav = IN2(163)
  maxdav = IN2(162)
  nrmdav = IN2(164)
  nrmz   = nrmdav 
  ! if gradients are required, increasing the cis convergence criterion
  if(jobtype.ge.1) nrmdav = nrmdav + 2


  ! es parameters
  cfac1 = 0.d0
  cfacJ = 1.d0
  cfacK = 1.d0

  ! According to JL, the following code was relevant
  ! only during program development (changed by AK).
  ! if(IN2(10).ne.0) then
  if(.false.) then
      open(unit = 930, file = "fort.14", status = "old")
      read(unit = 930, fmt = "(2I3,3X,F15.8)") i,j, cfac1
      c(i) = cfac1
      read(unit = 930, fmt = "(2I3,3X,F15.8)") i,j, cfacJ
      c(i) = cfacJ
      read(unit = 930, fmt = "(2I3,3X,F15.8)") i,j, cfacK
      c(i) = cfacK
      close(unit = 930)
      cfac1 = c(1)
      cfacJ = c(2)
      cfacK = c(3)
      if(iprtci.gt.1) then
        write(nb6,*) "cfac1",cfac1,cfacJ,cfacK
      endif
  endif

  ! print out the excited state variables
  call es_input_print(iprtci)

  return
end subroutine es_global_setup

!======================================================== 
subroutine es_control
  implicit real*8 (a-h,o-z)
  common  /INOPT2/ IN2(300)

  ! initialization
  I1     = 0
  I3     = 0
  I13    = 0
  I1n    = 1
  nden   = 1
  imult  = IN2(66)
  ispin  = IN2(142)
  nroots = IN2(139)
  iprtz  = IN2(41)
  iprtci = IN2(133)
  lstate = IN2(140)
  iuvcd  = IN2(146)
  nciref = IN2(136)

  ! number of ci derivatives
  if(jobtype.gt.0) then
    if(in2(159).le.0) in2(159) = 1
  endif

  ! davidson dimension
  if(in2(163).eq.0) in2(163) = 50
  if(in2(162).eq.0) in2(162) = in2(139)*in2(163)

  ! convergence criterion
  if(in2(164).eq.0) in2(164) = 7

  ! spin-adapted spin-flip cis
  sasfcis = .false.
  if(in2(77).eq.7) sasfcis = .true.
  if(sasfcis) then
    sf_xcis = .true.
    if(nciref.eq.1) sf_xcis = .false.
  endif

  ! singlet or triplet
  if(imult.eq.0 .or. imult.eq.3) then
    if(ispin.eq.1) then
      I1 = 1
    else if(ispin.eq.3) then
      I3 = 1
    else
      write(nb6,*) "incorrect multiplicity for closed-shell CIS",ispin
      stop
    endif
  else if(imult.eq.1 .or. imult.eq.2) then
    write(nb6,*) "let's check the validation of UCIS later"
    stop
  else
    write(nb6,*) "wrong multiplicity of the ground state"
    stop
  endif
  if(I3.eq.1) I1n = 0 
  if(I3.eq.1 .and. I1.eq.1) I13 = 1 
  if(imult.eq.3 .and. ispin.eq.1) I1n = 0

  ! open shell
  if(imult.eq.1 .or. imult.eq.2 .or. (imult.eq.3 .and. in2(77).eq.6)) nden = 2

  ! spin-flip
  is_sf = .false.
  if(imult.eq.3) is_sf = .true.

  return
end subroutine es_control

!======================================================== 
subroutine ovdim
  implicit real*8 (a-h,o-z)
  common  /LIMITS/ LM(7)
  common  /INOPT2/ IN2(300)
  common  /ORBITS/ Numb,Norbs,Nmos,Na,Nb

  NAlpha = Na
  NBeta  = Nb
  if(sasfcis .and. imult.eq.0) then
    NAlpha = Na + 1
    NBeta  = Nb - 1
  endif

  !-AK if(in2(131).eq.NA ) then
  if(in2(131).ge.NA) then
    if(Numb.ne.Na) then
      write(nb6,*) "check NUMB and NA",Num,Na
    endif
    NOa = NAlpha
    NVa = Norbs - NAlpha
    NOb = NBeta
    NVb = Norbs - NBeta
  else
    NOa = IN2(131)
    NVa = IN2(132)
    NOb = IN2(131) + NBeta - NAlpha
    NVb = IN2(132) + NAlpha - NBeta
  endif

  ! dimension
  NOrb   = NOa + NVa
  NOVa   = NOa*NVa
  NOVb   = NOb*NVb
  NOV    = NOa*NVa + NOb*NVb
  NBas   = LM(2)
  if(sasfcis) then
    NOV = 3+2*NOb+2*NVa+2*NVa*NOb
    if(sf_xcis) NOV = NOV + 2*NVa*NOb
    NOVa = NOV
    NOVb = NOV 
  else if(IN2(66).eq.3) then
    NOV  = 2*NOa*NVb
    NOVa = NOV/2
    NOVb = NOV/2
  endif

  if(nroots.gt.NOVa) then
    in2(139) = NOVa
    nroots   = NOVa
  endif

  return
endsubroutine ovdim

!======================================================== 
! save es properties
!======================================================== 
subroutine esref()
  use limit, only : lmprop,lmstat
  implicit real*8 (a-h,o-z)
  common  /CIPRP / ciprop(lmprop,lmstat),icisym(lmstat)

  ciprop(1,1:nroots-1) = Es(2:nroots,1)
  if(.not.allocated(jsym)) then
    write(nb6,*) "please allocate jsym"
    stop
  endif
  icisym(1:nroots-1)   = jsym(2:nroots)
  if(iuvcd.gt.0) then
    if(.not.(allocated(fv) .and. allocated(dv))) then
      write(nb6,*) "please allocate fv and dv"
      stop
    endif
    ciprop(3,1:nroots) = fv(1:nroots)
    ciprop(5,1:nroots) = dv(1:nroots)
  endif

  return
end subroutine esref

!======================================================
! symmetrize the molecular orbitals
!======================================================
subroutine es_mo_sym(V)
  implicit real*8 (a-h,o-z)
  real*8 V(*)
  common  /atoms / numat
  integer, dimension(:), allocatable :: nsym, isym

  is_sym = .false.
  allocate(nsym(lm3*nden),isym(3*numat))

  call mosym(V(jCa),V(jEa),lm2,lm3,norb,nsym,isym,isub,0)
  if(nden.eq.2) then
    call mosym(V(jCb),V(jEb),lm2,lm3,norb,nsym,isym,isub,0)
  endif

  if(isub.gt.0) is_sym = .true.

  deallocate(nsym)
  deallocate(isym)
  return
end subroutine es_mo_sym

!======================================================
subroutine es_input_print(iprint)
  implicit real*8 (a-h,o-z)
  integer iprint

  if(iprint.gt.1) then
     write(nb6,*)
     write(nb6,*) "NRoots: ", NRoots, " NOrb: ", NOrb, " NAlpha: ", NAlpha
     write(nb6,*) "NBeta:  ", NBeta,  " NOa:  ", NOa,  " NOb:    ", NOb
     write(nb6,*) "NVa:    ", NVa,    " NVb:  ", NVb,  " NBas:   ", NBas
     write(nb6,*) "kitdav: ", kitdav, " NOVa: ", NOVa, " NOVb:   ", NOVb
     write(nb6,*) "NOV:    ", NOV,    " I3:   ", I3,   " I1:     ", I1
     write(nb6,*) "LState: ", LState, " lm2   ", lm2,  " lm3:    ", lm3
     write(nb6,*) "I1n     ", I1n,    " NrmDav", nrmdav
  endif

  return
end subroutine es_input_print

!======================================================
! mapping the MOs with the old MOs
!======================================================
subroutine momap(Ca,Cb,Ca1,Cb1,icall)
  use limit, only : lmx
  implicit real*8 (a-h,o-z)
  integer icall
  real*8  Ca(lm2,*),Cb(lm2,*),Ca1(lm2,*),Cb1(lm2,*)
  integer, dimension(:), allocatable :: lmomap,llmomap
  real*8, dimension(:,:), allocatable :: pdot,pdot2
  common /CIMAP / IFMAP
  common /INOPT2/ IN2(300)
  common /CIMOS / imoci(lmx)

  ! flag for mapping failure
  ifmap = 0

  ! set up mapping thresh
  tol = 0.8d0
  if(in2(166).gt.0) tol = real(in2(166))/100d0

  allocate(pdot(norb,norb))
  if(nden.eq.2) allocate(pdot2(norb,norb))
  do i=1,norb
    do j=1,norb
      call vecdot(pdot(i,j),Ca(1,i),Ca1(1,j),lm3)
      if(nden.eq.2) call vecdot(pdot2(i,j),Cb(1,i),Cb1(1,j),lm3)
    enddo
  enddo

  iupdate = 0
  allocate(lmomap(norb))
  allocate(llmomap(norb))
  lmomap = 0
  llmomap = 0
  ikeepa = 0
  ikeepb = 0
  do i=1,norb
    ovlpmax = 0.d0
    do j=1,norb
      if(lmomap(j).eq.1) cycle
      do k=1,norb
        if(llmomap(k).eq.1) cycle
        if(dabs(pdot(k,j)).gt.dabs(ovlpmax)) then
          ovlpmax = pdot(k,j)
          jmax = j
          kmax = k
        endif
      enddo
    enddo
    imoci(kmax) = jmax
    lmomap(jmax) = 1
    llmomap(kmax)= 1
    if(jmax.ne.kmax) write(6,*) "orbital changed from:",jmax," to",kmax
    if(dabs(ovlpmax).lt.tol) then
      write(6,*) "mapping fails:",jmax,kmax,ovlpmax
      if((kmax.ge.nbeta .and. kmax.le.nalpha+1).and.(jmax.ge.nbeta .and. jmax.le.nalpha+1)) iupdate=1
      ifmap = 1
      if(in2(156).gt.0) icall = -1
    endif
    if(ovlpmax.lt.0) Ca1(:,jmax) = -Ca1(:,jmax)
    if(nden.eq.2) then
      write(nb6,*) "open shell orbital mapping NYI"
      stop
    endif
  enddo
  if(ifmap.eq.1) then
    write(nb6,*) "warning: orbital mapping failure"
    if(iprtz.ge.2) then
      do i=1,norb
        write(6,*) "mapped:",i,imoci(i),pdot(i,imoci(i))
      enddo
    endif
  endif

  if(iupdate.eq.0) then
  call veccopy(Ca,Ca1,lm2*norb)
  if(nden.eq.2) call veccopy(Cb,Cb1,lm2*norb)
  endif

  if(ifmap.eq.0) gread=.true.

  deallocate(lmomap)
  deallocate(llmomap)
  deallocate(pdot)
  if(nden.eq.2) deallocate(pdot2)
  return
end subroutine momap

end module es_global
