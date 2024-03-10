!
! Direct CI driver routines.
!
! Written by Axel Koslowski in 2000-2005 at MPI Muelheim.
!
! Shape-driven algorithm added by Axel Koslowski in 2012-2013.
!

Module gugadirect

  Use davidson
  Use gugaglobal
  Use gugaloops
  Use gugashape
  Use gugashape7
  Use gugashape8
  Use gugashape9
  Use gugashape10
  Use gugashape11
  Use gugautils
  Use matutils
  Use property

  Private

  Public :: DirectCI    ! Driver routine for solving the CI problem blockwise.
  Public :: DirectGrad  ! Calculate one- and two-particle (transition) density matrices for the analytic
                        ! CI gradient and non-adiabatic coupling vectors.
  Public :: DirectProp  ! Driver routine for population analysis and spectroscopic property calculation.

Contains

  ! Direct CI public driver routine.
  Subroutine DirectCI(Flags, iciref, SymLabel, CMO, Occ, F, W, NFirst, NLast, &
                      EE, Enuc, ciselt, CIProp, icisym, icall)
    Implicit None
    Type(ShavittControl),             Intent(InOut) :: Flags
    Integer,          Dimension(:,:), Intent(InOut) :: iciref
    Character(Len=4), Dimension(:),   Intent(In)    :: SymLabel
    Double Precision, Dimension(:,:), Intent(In)    :: CMO
    Double Precision, Dimension(:),   Intent(In)    :: Occ
    Double Precision, Dimension(:),   Intent(In)    :: F
    Double Precision, Dimension(:,:), Intent(In)    :: W
    Integer,          Dimension(:),   Intent(In)    :: NFirst, NLast
    Double Precision,                 Intent(InOut) :: EE
    Double Precision,                 Intent(In)    :: Enuc
    Double Precision,                 Intent(In)    :: ciselt
    Double Precision, Dimension(:,:), Intent(Out)   :: CIProp
    Integer,          Dimension(:),   Intent(Out)   :: icisym
    Integer,                          Intent(InOut) :: icall
    ! End of dummy parameters.
    Double Precision, Dimension(:,:), Allocatable   :: OneInt, TwoInt
    Double Precision, Dimension(:),   Allocatable   :: Hii
    Double Precision, Dimension(:),   Allocatable   :: Etmp
    Double Precision, Dimension(:,:), Allocatable   :: Ctmp
    Integer,          Dimension(:),   Allocatable   :: RootSort
    Integer                                         :: ncio, ncione
    Integer                                         :: isym, nciref, i
    Integer,          Dimension(1)                  :: imax
    Double Precision                                :: tuser, tsys, twall

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(//A)') ' Entering direct CI section.'
       Call GetTime(tuser, tsys, twall)
    End If

    ncio   = UBound(Flags%WhoIsWho,1)
    ncione = ncio * (ncio+1) / 2

    Allocate(OneInt(ncio,   ncio))
    Allocate(TwoInt(ncione, ncione))

    Call ActiveIntegralTransformation(CMO(:,Flags%WhoIsWho), F, W, OneInt, TwoInt, NFirst, NLast, Flags%PrintLevel)
    Allocate(Hii(Flags%TotalCSF))

    If (Flags%Algorithm > 0) Then
       Call ShapeDrivenHii(Flags, Hii, OneInt, TwoInt)
    Else
       Call CalcHii(Flags, Hii, OneInt, TwoInt, Occ(Flags%WhoIsWho), 0)
    End If

    If (Flags%PrintLevel >= 4) Then
       Write(Stdout,'(//1X,A/)') 'Diagonal elements of CI Hamiltonian:'
       Call PrintVector(Hii, Stdout)
    End If

    If (Flags%NumberOfCIRoots > 0) Then
       Call DoDirectCI(Flags, 0, Hii, OneInt, TwoInt, Occ, 0, SymLabel, Flags%CIE, Flags%CIC, icall)
       Call Symmetrize(Flags, Flags%CIE, Flags%CIC, SymLabel, Flags%PrintLevel)
    Else
       Allocate(Etmp(Flags%TotalRoots))
       Allocate(Ctmp(Flags%TotalCSF, Flags%TotalRoots))
       Allocate(RootSort(Flags%TotalRoots))
       RootSort = (/(i,i=1,Flags%TotalRoots)/)
       Ctmp     = 0.D0
       Do isym=1, Flags%NumSym
          If (Flags%NRoots(isym) > 0) Then
             Call DoDirectCI(Flags, Flags%FirstCSF(isym)-1, Hii(Flags%FirstCSF(isym):Flags%LastCSF(isym)),               &
                             OneInt, TwoInt, Occ, isym, SymLabel, Etmp(Flags%FirstRoot(isym):Flags%LastRoot(isym)),      &
                             Ctmp(Flags%FirstCSF(isym):Flags%LastCSF(isym), Flags%FirstRoot(isym):Flags%LastRoot(isym)), &
                             icall)
          End If
       End Do
       Call SortCIE(Etmp, RootSort)
       Flags%CIE = Etmp  (RootSort)
       Flags%CIC = Ctmp(:,RootSort)
       Deallocate(RootSort)
       Deallocate(Ctmp)
       Deallocate(Etmp)
    End If

    Deallocate(OneInt)
    Deallocate(TwoInt)

    If (Flags%mciref >= 3) Then
       nciref = Flags%nciref
       Call RefCheck(Flags, iciref, FLags%DRT, ciselt, icall)
       If (nciref < Flags%nciref .Or. icall < 0)  Return
    End If

    Do i=1, Flags%TotalRoots
       imax = MaxLoc(Abs(Flags%CIC(:,i)))
       Flags%RootSym(i) = Flags%SymVec(imax(1))
    End Do

    ! Find CI root of interest for negative lroot:
    If (Flags%LRoot < 0) Then
       Call RootOfInterest(Flags)
    End If

    Call PrintCIResults(Flags, Hii, EE, Enuc, SymLabel)

    Deallocate(Hii)

    EE = EE + Flags%CIE(Flags%iState)

    CIProp        = 0.D0
    i             = Min(Flags%TotalRoots, UBound(CIProp,2))
    CIProp(1,1:i) = Flags%CIE(1:i)
    icisym  (1:i) = Flags%RootSym(1:i)

    If (Flags%PrintLevel >= 2) then
       Write(Stdout,'(//A)') ' Direct CI section done.'
       Call PrintTime(tuser, tsys, twall)
    End If

    Return
  End Subroutine DirectCI



  ! Find eigenvalues and eigenvectors of one block of the CI Hamiltonian.
  Subroutine DoDirectCI(Flags, i0, Hii, OneInt, TwoInt, Occ, isym, SymLabel, CIE, CIC, icall)
    Implicit None
    Type(ShavittControl),             Intent(InOut) :: Flags
    Integer,                          Intent(In)    :: i0
    Double Precision, Dimension(:),   Intent(In)    :: Hii
    Double Precision, Dimension(:,:), Intent(In)    :: OneInt, TwoInt
    Double Precision, Dimension(:),   Intent(In)    :: Occ
    Integer,                          Intent(In)    :: isym
    Character(Len=4), Dimension(:),   Intent(In)    :: SymLabel
    Double Precision, Dimension(:),   Intent(Out)   :: CIE
    Double Precision, Dimension(:,:), Intent(Out)   :: CIC
    Integer,                          Intent(InOut) :: icall
    ! End of dummy parameters.
    Double Precision, Dimension(:),   Pointer       :: dummy
    Double Precision, Dimension(1),   Target        :: dummyTarget
    Integer,          Dimension(:),   Pointer       :: idummy
    Integer,          Dimension(1),   Target        :: idummyTarget
    Integer                                         :: nroots, cidiag
    Double Precision                                :: tuser, tsys, twall

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(//1X,A)') 'Solving CI eigenvalue problem (direct algorithm)...'
       Call GetTime(tuser, tsys, twall)
    End If

    dummy  => dummyTarget
    idummy => idummyTarget
    nroots =  UBound(CIE,1)
    cidiag =  Flags%WhichDiagonalizer
    If (cidiag == 0) Then
       If (nroots == 1) Then
          cidiag = 10
       Else
          cidiag = 13
       End If
    End If

    Select Case (cidiag)
    Case (10)
       Call Davidson1(Flags, Hii, dummy, idummy, idummy, OneInt, TwoInt, Occ(Flags%WhoIsWho), isym, Flags%iRefConf(1:Flags%nRefCSF), &
                      i0, SymLabel, CIE, CIC, nroots, Flags%MinDav, Flags%MaxDav, Flags%KitDav, Flags%QNorm, Flags%Algorithm,        &
                      Stdout, Flags%PrintLevel, icall)
    Case (11)
       Call Davidson2(Flags, Hii, dummy, idummy, idummy, OneInt, TwoInt, Occ(Flags%WhoIsWho), isym, Flags%iFirstRef(1:Flags%nciref), &
                      Flags%iLastRef(1:Flags%nciref), Flags%iRefConf(1:Flags%nRefCSF), i0, SymLabel, CIE(1), CIC(:,1), Flags%LRoot,  &
                      Flags%MinDav, Flags%MaxDav, Flags%KitDav, Flags%QNorm, Flags%Algorithm, Stdout, Flags%PrintLevel, icall)
       ! One state has been calculated.
       Flags%iState = 1
       Flags%jState = 1
    Case (12,13)
       Call Davidson3(Flags, Hii, dummy, idummy, idummy, OneInt, TwoInt, Occ(Flags%WhoIsWho), isym, Flags%iRefConf(1:Flags%nRefCSF), &
                      i0, SymLabel, CIE, CIC, nroots, Flags%MinDav, Flags%MaxDav, Flags%KitDav, Flags%QNorm, Flags%Algorithm,        &
                      Stdout, Flags%PrintLevel, icall)
    Case Default
       Write(Stdout,'(/1X,A,I2,A)') 'Diagonalizer ', Flags%WhichDiagonalizer, ' cannot be used with direct CI.'
       Stop 'DirectCI'
    End Select

    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)

    Return
  End Subroutine DoDirectCI



  ! Public driver routine to calculate one- and two-particle
  ! (transition ) density matrices for the analytic CI gradient
  ! and non-adiabatic coupling vectors.
  !
  Subroutine DirectGrad(Flags, OneBodyDen, TwoBodyDen, Occ)
    Implicit None
    Type(ShavittControl),             Intent(In)  :: Flags
    Double Precision, Dimension(:,:), Intent(Out) :: OneBodyDen
    Double Precision, Dimension(:),   Intent(Out) :: TwoBodyDen
    Double Precision, Dimension(:),   Intent(In)  :: Occ
    ! End of dummy parameters.
    Integer                                       :: nactive, ncione, i, j

    i = Flags%iState
    j = Flags%jState

    If (i == j) Then
       If (Flags%Algorithm > 0) Then
          Call ShapeDrivenDirectGradDen(Flags, OneBodyDen, TwoBodyDen, Flags%CIC(:,i), Flags%RootSym(i))
       Else
          Call CalcGradDen(Flags, OneBodyDen, TwoBodyDen, Flags%CIC(:,i), Occ(Flags%WhoIsWho), Flags%RootSym(i))
       End If
    Else
       If (Flags%Algorithm > 0) Then
          If (Flags%NumSym == 1) Then
             Call ShapeDrivenDirectNACDenC1(Flags, OneBodyDen, TwoBodyDen, Flags%CIC(:,i), Flags%CIC(:,j))
          Else
             If (.Not. Associated(Flags%CombiListAll13ijkl))  Call InitAllCombinations(Flags)
             Call ShapeDrivenDirectNACDen(Flags, OneBodyDen, TwoBodyDen, Flags%CIC(:,i), Flags%CIC(:,j))
          End If
       Else
          Call CalcNACDen(Flags, OneBodyDen, TwoBodyDen, Flags%CIC(:,i), Flags%CIC(:,j), Occ(Flags%WhoIsWho))
       End If
    End If

    If (Flags%PrintLevel >= 4) Then
       Write(Stdout,'(//1X,"One-particle density matrix:"/)')
       Call PrintMatrix(OneBodyDen, Stdout)
       Write(Stdout,'(//1X,"Two-particle density matrix:"/)')
       nactive = Flags%DRT(1)%k
       ncione  = (nactive * (nactive+1)) / 2
       Call PrintPackedLower(TwoBodyDen, ncione, Stdout)
    End If

    Return
  End Subroutine DirectGrad



  ! Public driver routine for natural orbital generation, population analysis
  ! and spectroscopic property calculation (UV and CD).
  !
  Subroutine DirectProp(Flags, NAT, Coord, NFirst, NLast, Core, zs, zp, zd, nsp, nd, CMO, Occ, &
                        SymLabel, CIProp, icisym, PAO, iuvcd, imcd, ipop, inatur, icall)
    Implicit None
    Type(ShavittControl),               Intent(In)    :: Flags
    Integer,          Dimension(:),     Intent(In)    :: NAT
    Double Precision, Dimension(:,:),   Intent(In)    :: Coord
    Integer,          Dimension(:),     Intent(In)    :: NFirst, NLast
    Double Precision, Dimension(:),     Intent(In)    :: Core, zs, zp, zd
    Integer,          Dimension(:),     Intent(In)    :: nsp, nd
    Double Precision, Dimension(:,:),   Intent(In)    :: CMO
    Double Precision, Dimension(:),     Intent(In)    :: Occ
    Character(Len=4), Dimension(:),     Intent(In)    :: SymLabel
    Double Precision, Dimension(:,:),   Intent(Out)   :: CIProp
    Integer,          Dimension(:),     Intent(Out)   :: icisym
    Double Precision, Dimension(:,:),   Intent(Out)   :: PAO
    Integer,                            Intent(In)    :: iuvcd, imcd, ipop, inatur
    Integer,                            Intent(InOut) :: icall
    ! End of dummy parameters.
    Integer                                           :: nPerm, mTrans, nTrans
    Integer                                           :: ncio
    Double Precision, Dimension(:,:,:), Allocatable   :: Perm, Trans
    Double Precision                                  :: tuser, tsys, twall

    If (inatur <=0 .And. ipop <= 0 .And. iuvcd <= 0)  Return

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(//1X,A)') 'Natural orbitals, population analysis and properties (direct algorithm)...'
       Call GetTime(tuser, tsys, twall)
    End If

    ncio   = UBound(Flags%WhoIsWho,1)
    nPerm  = Flags%TotalRoots
    mTrans = 0
    nTrans = 0

    Allocate(Perm(ncio, ncio, nPerm))

    If (iuvcd >= 3 .And. Flags%TotalRoots > 1) Then
       mTrans = Flags%TotalRoots-1
       nTrans = Flags%TotalRoots
       Allocate(Trans(ncio, ncio, Flags%TotalRoots*(Flags%TotalRoots-1)/2))
    Else If (iuvcd >= 2 .And. Flags%TotalRoots > 1) Then
       mTrans = 1
       nTrans = Flags%TotalRoots
       Allocate(Trans(ncio, ncio, Flags%TotalRoots-1))
    Else
       Allocate(Trans(1,1,1))
    End If

    If (Flags%Algorithm > 0) Then
       Call ShapeDrivenDirectPropDen(Flags, Perm, Trans, nPerm, mTrans, nTrans, Flags%CIC)
    Else
       Call CalcPropDen(Flags, Perm, Trans, nPerm, mTrans, nTrans, Flags%CIC, Occ(Flags%WhoIsWho), 0)
    End If

    Call NaturalOrb(Perm, Occ, Flags%WhoIsWho, nPerm, Stdout, Flags%iState, inatur, icall)

    Call Population(Perm, PAO, CMO, NAT, NFirst, NLast, Core, Occ(1:Flags%HOMOindex), &
                    Flags%WhoIsWho, nPerm, Stdout, Flags%iState, ipop, icall)

    Call Properties(Flags%CIE, Flags%RootSym, SymLabel, Perm, Trans, CMO, NAT, Coord, &
                    NFirst, NLast, Core, zs, zp, zd, nsp, nd, Occ(1:Flags%HOMOindex), &
                    Flags%WhoIsWho, CIProp, icisym, nPerm, mTrans, nTrans, Stdout,    &
                    Flags%iop, iuvcd, imcd, icall)

    Deallocate(Perm)
    Deallocate(Trans)

    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)

    Return
  End Subroutine DirectProp


End Module gugadirect
