!
! Subroutines useful in all types of CI calculations.
!
! First version by Michael E. Beck in 1998 at Zurich University.
!
! Final version by Axel Koslowski in 2000-2005 at MPI Muelheim.
!

Module gugautils

  Use gugaglobal
  Use matutils

  Private

  Public :: ActiveIntegralTransformation  ! Integral transformation from AO into active MO basis.
  Public :: SortCIE                       ! Sort CI roots according to increasing energy.
  Public :: Symmetrize                    ! For instance in D2d, CI eigenvectors may not be C2v-adapted.
  Public :: RootOfInterest                ! Find CI root of interest for negative lroot.
  Public :: RefCheck                      ! Add references so that the fraction in each CI root is sufficiently large.
  Public :: FindGelfandStates
  Public :: FindStepVector
  Public :: CmpBraKet
  Public :: PrintCIResults                ! Energy, leading configurations, and complete eigenvectors.
  Public :: PrintTime                     ! Print user and system CPU and wall clock time.
  Public :: PrintGenerators
  Public :: CanonicalIndex                ! Calculate canonical index from index pair.
  Public :: DoubleCanonicalIndex
  Public :: Split1
  Public :: Split2

Contains

  !
  ! Integral transformation from AO into active orbital basis.
  !

  ! Public driver routine.
  Subroutine ActiveIntegralTransformation(ActiveCMO, FockMatrix, TwoElectronAOIntegrals,  &
             &                            OneElectronMOIntegrals, TwoElectronMOIntegrals, &
             &                            NFirst, NLast, PrintLevel)
    Implicit None
    Double Precision, Dimension(:,:), Intent(In)    :: ActiveCMO
    Double Precision, Dimension(:),   Intent(In)    :: FockMatrix
    Double Precision, Dimension(:,:), Intent(In)    :: TwoElectronAOIntegrals
    Double Precision, Dimension(:,:), Intent(Out)   :: OneElectronMOIntegrals
    Double Precision, Dimension(:,:), Intent(Out)   :: TwoElectronMOIntegrals
    Integer, Dimension(:),            Intent(In)    :: NFirst, NLast
    Integer,                          Intent(In)    :: PrintLevel
    ! End of dummy parameters.
    Double Precision, Dimension(:,:), Allocatable   :: fock, CC
    Integer                                         :: nao, naoao, nact, ncione
    Double Precision                                :: tuser, tsys, twall

    If (PrintLevel >= 2) Then
       Write(Stdout,'(//1X,A)') 'Integral transformation (active MOs)...'
       Call GetTime(tuser, tsys, twall)
    End If

    nao    = NLast(UBound(NLast,1))
    naoao  = UBound(TwoElectronAOIntegrals,1)
    nact   = UBound(ActiveCMO,2)
    ncione = nact * (nact+1) / 2

    ! Form active MO basis one-electron integrals:
    Allocate(fock(nao,nao))
    Call FormSquare(fock, FockMatrix)
    Call SimilarityTransform(OneElectronMOIntegrals, fock, ActiveCMO)
    Deallocate(fock)

    ! Form active MO basis two-electron integrals:
    Allocate(CC(naoao,ncione))
    Call FormMOProducts(ActiveCMO, CC, NFirst, NLast)
    Call SimilarityTransform(TwoElectronMOIntegrals, TwoElectronAOIntegrals, CC)
    Deallocate(CC)

    If (PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)

    If (PrintLevel > 5) Then
       Write(Stdout,'(//A/)') ' Active MO basis one-electron integrals:'
       Call PrintMatrix(OneElectronMOIntegrals, Stdout)
       Write(Stdout,'(//A/)') ' Active MO basis two-electron integrals:'
       Call PrintMatrix(TwoElectronMOIntegrals, Stdout)
    End If

    Return
  End Subroutine ActiveIntegralTransformation



  Subroutine SimilarityTransform(C, A, B)  ! C = BT * A * B
    Implicit None
    Double Precision, Dimension(:,:)              :: C
    Double Precision, Dimension(:,:)              :: A, B
    ! End of dummy parameters.
    Double Precision, Dimension(:,:), Allocatable :: P
    Integer                                       :: k, l, ldb

    k   = UBound(A,1)
    l   = UBound(B,2)
    ldb = UBound(B,1)

    Allocate(P(k,l))
    Call DGEMM('N', 'N', k, l, k, 1.D0, A(1,1), k,   B(1,1), ldb, 0.D0, P,      k)
    Call DGEMM('T', 'N', l, l, k, 1.D0, B(1,1), ldb, P(1,1), k,   0.D0, C(1,1), l)
    Deallocate(P)

    Return
  End Subroutine SimilarityTransform



  Subroutine FormSquare(square, packed)
    Implicit None
    Double Precision, Dimension(:,:) :: square
    Double Precision, Dimension(:)   :: packed
    ! End of dummy parameters.
    Integer                          :: i, j, ij

    ij = 0
    Do j=1, UBound(square,2)
       Do i=1, j
          ij = ij + 1
          square(i,j) = packed(ij)
          square(j,i) = packed(ij)
       End Do
    End Do

    Return
  End Subroutine FormSquare



  Subroutine FormMOProducts(CMO, CC, NFirst, NLast)
    Implicit None
    Double Precision, Dimension(:,:) :: CMO, CC
    Integer,          Dimension(:)   :: NFirst, NLast
    ! End of dummy parameters.
    Integer                          :: i, j, ij

    ij = 0
    Do i=1, UBound(CMO,2)
       Do j=1, i
          ij = ij + 1
          Call FormMOProduct(CMO, CC, NFirst, NLast, i, j, ij)
       End Do
    End Do

    Return
  End Subroutine FormMOProducts



  Subroutine FormMOProduct(CMO, CC, NFirst, NLast, i, j, ij)
    Implicit None
    Double Precision, Dimension(:,:) :: CMO, CC
    Integer,          Dimension(:)   :: NFirst, NLast
    Integer                          :: i, j, ij
    ! End of dummy parameters.
    Integer                          :: ia, mu, nu, munu

    munu = 0
    Do ia=1, UBound(NFirst,1)
       Do mu=NFirst(ia), NLast(ia)
          Do nu=NFirst(ia), mu-1
             munu = munu + 1
             CC(munu,ij) = CMO(mu,i) * CMO(nu,j) + CMO(nu,i) * CMO(mu,j)
          End Do
          munu = munu + 1
          CC(munu,ij) = CMO(mu,i) * CMO(mu,j)
       End Do
    End Do

    Return
  End Subroutine FormMOProduct



  !
  ! Public driver routines for post-processing CI results.
  !

  ! Sort according to increasing energy.
  Subroutine SortCIE(CIE, isort)
    Double Precision, Dimension(:) :: CIE
    Integer,          Dimension(:) :: isort
    ! End of dummy parameters.
    Integer                        :: i, j, k, l, n

    n = UBound(CIE,1)

    Do i=1, n-1
       k = i
       Do j=i+1, n
          If (CIE(isort(j)) < CIE(isort(k)))  k = j
       End Do
       If (k /= i) Then
          l        = isort(i)
          isort(i) = isort(k)
          isort(k) = l
       End If
    End Do

    Return
  End Subroutine SortCIE



  ! For instance in D2d, CI eigenvectors may not be C2v-adapted.
  Subroutine Symmetrize(Flags, CIE, CIC, symLabel, prtLev)
    Implicit None
    Type(ShavittControl)                            :: Flags
    Double Precision, Dimension(:)                  :: CIE
    Double Precision, Dimension(:,:)                :: CIC
    Character(Len=4), Dimension(:)                  :: symLabel
    Integer                                         :: prtLev
    ! End of dummy parameters.
    Integer, Parameter                              :: nrmTol = 10
    Integer, Parameter                              :: cieTol =  4
    Double Precision, Dimension(UBound(SymLabel,1)) :: symNorm
    Logical, Dimension(UBound(CIC,2))               :: symOkay
    Integer, Dimension(2)                           :: imax
    Double Precision                                :: ca, cb, r, sina, cosa
    Integer                                         :: ni, nj, i, j, k, ja, jb

    ni = UBound(CIC,1)
    nj = UBound(CIC,2)

    If (prtLev >= 3) Then
       Write(Stdout,'(/1X,A/)')  'Symmetry of CI roots:'
       Write(Stdout,'(2X,8A15)') symLabel
    End If

    Do j=1, nj
       symNorm = 0.D0
       Do i=1, ni
          k = Flags%SymVec(i)
          symNorm(k) = symNorm(k) + CIC(i,j) * CIC(i,j)
       End Do
       If (prtLev >= 3) &
          Write(Stdout,'(1X,I4,8F15.10)') j, symNorm
       imax(1:1)  = MaxLoc(symNorm)
       symOkay(j) = symNorm(imax(1)) > 1.D0 - 1.D1**(-nrmTol)
    End Do

    If (prtLev >= 3)  Write(Stdout,*)

    j = 1
    Do While (j < nj)
       If (symOkay(j)) Then
          j = j + 1
       Else
          If (.Not. symOkay(j+1) .And. CIE(j+1)-CIE(j) < 1.D1**(-cieTol)) Then
             If (prtLev >= 3) &
                Write(Stdout,'(1X,A,I4,A,I4,A)') 'Symmetrizing CI roots ', j, ' and ', j+1, '.'
             imax = MaxLoc(Abs(CIC(:,j:j+1)))
             i    = imax(1)
             ja   = j - 1 + imax(2)
             jb   = j + 2 - imax(2)
             ca   = CIC(i,ja)
             cb   = CIC(i,jb)
             r    = Sqrt(ca*ca + cb*cb)
             cosa = ca / r
             sina = cb / r
             Do i=1, ni
                ca        = CIC(i,ja)
                cb        = CIC(i,jb)
                CIC(i,ja) = cosa * ca + sina * cb
                CIC(i,jb) = cosa * cb - sina * ca
             End Do
             j = j + 2
          Else
             If (prtLev > -9) &
                Write(Stdout,'(1X,A,I4,A)') 'WARNING: CI root ', j, ' not symmetrized.'
             j = j + 1
          End If
       End If
    End Do

    If (prtLev >= 3)  Write(Stdout,'(/)')

    Return
  End Subroutine Symmetrize



  ! Find CI root of interest for negative lroot.
  Subroutine RootOfInterest(Flags)
    Implicit None
    Type(ShavittControl), Intent(InOut) :: Flags
    ! End of dummy parameters.
    Double Precision                    :: flarge, frac
    Integer                             :: iref, i1, i2, ilarge, i

    iref = -Flags%LRoot

    If (iref > Flags%nciref) Then
       Write(Stdout,'(1X,A)') 'Negative lroot inconsistent with number of reference configurations.'
       Stop 'gugadriver::DirectCI'
    End If

    i1 = Flags%iFirstRef(iref)
    i2 = Flags%iLastRef(iref)

    If (i1 > i2) Then
       Write(Stdout,'(1X,A)') 'Negative lroot points to reference not represented by the DRT.'
       Stop 'gugadriver::DirectCI'
    End If

    If (Flags%PrintLevel >= 1) Then
       Write(Stdout,'(//1X,"List of CI basis functions for reference configuration ",I3,":")') iref
       Write(Stdout,'( /1X,10I12)') Flags%iRefConf(i1:i2)
    End If

    ilarge = 1
    flarge = Dot_Product(Flags%CIC(Flags%iRefConf(i1:i2), 1), Flags%CIC(Flags%iRefConf(i1:i2), 1))
    Write(Stdout,'(//1X,"Fraction of these CI basis functions in root ",I4,":",F10.2," %")') ilarge, 1.D2*flarge

    Do i=2, Flags%TotalRoots
       frac = Dot_Product(Flags%CIC(Flags%iRefConf(i1:i2), i), Flags%CIC(Flags%iRefConf(i1:i2), i))
       Write(Stdout,'(1X,"Fraction of these CI basis functions in root ",I4,":",F10.2," %")') i, 1.D2*frac
       If (frac > flarge) Then
          ilarge = i
          flarge = frac
       End If
    End Do

    Flags%iState = ilarge
    Flags%jState = ilarge

    If (Flags%PrintLevel > -9) Then
       Write(Stdout,'(//1X,A,I4)') 'CI root of interest from negative lroot: ', ilarge
    End If

    Return
  End Subroutine RootOfInterest



  ! Add references so that the fraction in each CI root is above ciselt.
  Subroutine RefCheck(Flags, iciref, DRT, ciselt, icall)
    Implicit None
    Type(ShavittControl),             Intent(InOut) :: Flags
    Integer,          Dimension(:,:), Intent(InOut) :: iciref
    Type(PaldusRow),  Dimension(:),   Intent(In)    :: DRT
    Double Precision,                 Intent(In)    :: ciselt
    Integer,                          Intent(InOut) :: icall
    ! End of dummy parameters.
    Double Precision, Dimension(:),   Allocatable   :: fraction
    Logical,          Dimension(:),   Allocatable   :: IsRefConf
    Integer,          Dimension(:),   Pointer       :: itmp
    Integer                                         :: nciref, ismall, ilarge, i, j
    Double Precision                                :: c

    Allocate(fraction(Flags%TotalRoots))
    Allocate(IsRefConf(Flags%TotalCSF))
    IsRefConf                                     = .False.
    Where (Flags%iRefConf(1:Flags%nRefCSF) > 0) &
       IsRefConf(Flags%iRefConf(1:Flags%nRefCSF)) = .True.

    nciref = Flags%nciref

    Do
       ! Determine fraction of references in all roots:
       fraction = 0.D0
       Do i=1, Flags%TotalRoots
          Do j=1, Flags%nRefCSF
             If (Flags%iRefConf(j) > 0) Then
                c           = Flags%CIC(Flags%iRefConf(j), i)
                fraction(i) = fraction(i) + c * c
             End If
          End Do
       End Do

       ! Print fractions:
       If (Flags%PrintLevel >= 2) Then
          If (Flags%nciref == nciref) Then
             Write(Stdout,'(//1X,A/)') 'Initial fractions:'
          Else
             Write(Stdout,'(//1X,A,I2,A/)') 'Fractions after adding ', Flags%nciref-nciref, ' reference(s):'
          End If
          Write(Stdout,'(1X,A)') 'CI root   fraction'
          Do i=1, Flags%TotalRoots
             Write(Stdout,'(1X,I5,F11.2," %")') i, 100.D0*fraction(i)
          End Do
       End If

       ! Find root with smallest fraction:
       ismall = 1
       Do i=2, Flags%TotalRoots
          If (fraction(i) < fraction(ismall))  ismall = i
       End Do
       If (Flags%PrintLevel >= 2)  Write(Stdout,'(/1X,"Root with smallest fraction:",I9)') ismall

       ! Exit loop if fraction is okay:
       If (fraction(ismall) >= ciselt)  Exit

       ! Find non-reference CI basis functions with largest coefficient:
       Do i=1, Flags%TotalCSF
          If (.Not. IsRefConf(i)) Then
             ilarge = i
             Exit
          End If
       End Do
       Do i=ilarge+1, Flags%TotalCSF
          If (.Not. IsRefConf(i) .And. Abs(Flags%CIC(i,ismall)) > Abs(Flags%CIC(ilarge,ismall)))  ilarge = i
       End Do
       If (Flags%PrintLevel >= 2)  Write(Stdout,'(1X,"CSF with largest coefficient:",I8)') ilarge

       ! Add corresponding reference configuration:
       If (Flags%nciref == UBound(iciref,2)) Then
          Write(Stdout,'(1X,A)') 'Cannot add another reference configuration.'
          icall = -1
          Return
       End If
       Flags%nciref = Flags%nciref + 1
       Call FindOccupation(DRT, iciref(:,Flags%nciref), 1, Flags%WalkIndex(ilarge))
       If (Flags%PrintLevel >= 2) Then
          Write(Stdout,'(/1X,A)')  'New reference configuration:'
          Write(Stdout,'( 30I4)')  Flags%WhoIsWho
          Write(Stdout,'( 30I4)')  iciref(:,Flags%nciref)
       End If

       ! Find all corresponding CI basis functions:
       Allocate(itmp(Flags%nciref))
       itmp(1:Flags%nciref-1) = Flags%iFirstRef
       Deallocate(Flags%iFirstRef)
       Flags%iFirstRef => itmp

       Allocate(itmp(Flags%nciref))
       itmp(1:Flags%nciref-1) = Flags%iLastRef
       Deallocate(Flags%iLastRef)
       Flags%iLastRef => itmp

       Flags%iFirstRef(Flags%nciref) = Flags%nRefCSF + 1
       Call FindGelfandStates(DRT, iciref(:,Flags%nciref), Flags%iRefWalk, Flags%nRefCSF, 1, 1)
       Flags%iLastRef (Flags%nciref) = Flags%nRefCSF
       Deallocate(Flags%iRefConf)
       Allocate(Flags%iRefConf(Flags%nRefCSF))
       Flags%iRefConf = Flags%IndVec(Flags%iRefWalk(1:Flags%nRefCSF))
       Where (Flags%iRefConf(1:Flags%nRefCSF) > 0) &
          IsRefConf(Flags%iRefConf(1:Flags%nRefCSF)) = .True.
       If (Flags%PrintLevel >= 2) Then
          Write(Stdout,'(/1X,A)')  'Corresponding CSFs:'
          Write(Stdout,'( 12I10)') Flags%iRefConf(Flags%iFirstRef(Flags%nciref) : Flags%iLastRef(Flags%nciref))
       End If
    End Do

    Deallocate(IsRefConf)
    Deallocate(fraction)
    Return
  End Subroutine RefCheck



  ! Find lexical indices of all Gelfand states corresponding
  ! to a given occupation pattern.
  Recursive Subroutine FindGelfandStates(DRT, iOcc, iCSF, nCSF, irow, ilex)
    Implicit None
    Type(PaldusRow), Dimension(:), Intent(In)    :: DRT
    Integer,         Dimension(:), Intent(In)    :: iOcc
    Integer,         Dimension(:), Pointer       :: iCSF
    Integer,                       Intent(InOut) :: nCSF
    Integer,                       Intent(In)    :: irow, ilex
    ! End of dummy parameters.
    Integer,         Dimension(:), Pointer       :: jCSF
    Integer                                      :: k

    If (irow == 0) Return

    k = DRT(irow)%k

    If (k == 0) Then
       If (nCSF == UBound(iCSF,1)) Then
          Allocate(jCSF(nCSF+10))
          jCSF(1:nCSF) = iCSF
          Deallocate(iCSF)
          iCSF => jCSF
       End If
       nCSF = nCSF + 1
       iCSF(nCSF) = ilex
       Return
    End If

    Select Case (iOcc(k))
    Case (0)
       Call FindGelfandStates(DRT, iOcc, iCSF, nCSF, DRT(irow)%jd(0), ilex+DRT(irow)%yd(0))
    Case (1)
       Call FindGelfandStates(DRT, iOcc, iCSF, nCSF, DRT(irow)%jd(1), ilex+DRT(irow)%yd(1))
       Call FindGelfandStates(DRT, iOcc, iCSF, nCSF, DRT(irow)%jd(2), ilex+DRT(irow)%yd(2))
    Case (2)
       Call FindGelfandStates(DRT, iOcc, iCSF, nCSF, DRT(irow)%jd(3), ilex+DRT(irow)%yd(3))
    End Select

    Return
  End Subroutine FindGelfandStates



  ! Find orbital occupation pattern corresponding to a walk
  ! with a given lexical index.
  Recursive Subroutine FindOccupation(DRT, iOcc, irow, ilex)
    Implicit None
    Type(PaldusRow), Dimension(:), Intent(In)    :: DRT
    Integer,         Dimension(:), Intent(Out)   :: iOcc
    Integer,                       Intent(In)    :: irow, ilex
    ! End of dummy parameters.
    Integer                                      :: k

    If (irow == 0) Then
       Write(Stdout,'(1X,A)') 'ERROR: Inconsistent DRT.'
       Stop 'gugadrt::FindOccupation'
    End If

    k = DRT(irow)%k

    If (k == 0) Return

    If (DRT(irow)%yd(3) < ilex) Then
       iOcc(k) = 2
       Call FindOccupation(DRT, iOcc, DRT(irow)%jd(3), ilex - DRT(irow)%yd(3))
    Else If (DRT(irow)%yd(2) < ilex) Then
       iOcc(k) = 1
       Call FindOccupation(DRT, iOcc, DRT(irow)%jd(2), ilex - DRT(irow)%yd(2))
    Else If (DRT(irow)%yd(1) < ilex) Then
       iOcc(k) = 1
       Call FindOccupation(DRT, iOcc, DRT(irow)%jd(1), ilex - DRT(irow)%yd(1))
    Else
       iOcc(k) = 0
       Call FindOccupation(DRT, iOcc, DRT(irow)%jd(0), ilex - DRT(irow)%yd(0))
    End If

    Return
  End Subroutine FindOccupation



  ! Find the step vector corresponding to a walk with a given lexical index.
  Recursive Subroutine FindStepVector(DRT, dvec, irow, ilex)
    Implicit None
    Type(PaldusRow), Dimension(:), Intent(In)    :: DRT
    Integer,         Dimension(:), Intent(InOut) :: dvec
    Integer,                       Intent(In)    :: irow, ilex
    ! End of dummy parameters.
    Integer                                      :: k

    If (irow == 0) Then
       Write(Stdout,'(1X,A)') 'ERROR: Inconsistent DRT.'
       Stop 'gugadrt::FindStepVector'
    End If

    k = DRT(irow)%k

    If (k == 0) Return

    If (DRT(irow)%yd(3) < ilex) Then
       dvec(k) = 3
       Call FindStepVector(DRT, dvec, DRT(irow)%jd(3), ilex - DRT(irow)%yd(3))
    Else If (DRT(irow)%yd(2) < ilex) Then
       dvec(k) = 2
       Call FindStepVector(DRT, dvec, DRT(irow)%jd(2), ilex - DRT(irow)%yd(2))
    Else If (DRT(irow)%yd(1) < ilex) Then
       dvec(k) = 1
       Call FindStepVector(DRT, dvec, DRT(irow)%jd(1), ilex - DRT(irow)%yd(1))
    Else
       dvec(k) = 0
       Call FindStepVector(DRT, dvec, DRT(irow)%jd(0), ilex - DRT(irow)%yd(0))
    End If

    Return
  End Subroutine FindStepVector



  ! Establish an ordering of partial loops with
  ! different (ibra, iket) characteristics.
  ! This function returns a value
  !   <  0 if part <  (ibra, iket)
  !   == 0 if part == (ibra, iket)
  !   >  0 if part >  (ibra, iket).
  Integer Function CmpBraKet(part, ibra, iket)
    Implicit None
    Type(PartialLoop) :: part
    Integer           :: ibra, iket
    ! End of dummy parameters.
    Integer           :: idiff

    idiff = part%ibra - ibra

    If (idiff /= 0) Then
       CmpBraKet = idiff
       Return
    End If

    CmpBraKet = part%iket - iket
  End Function CmpBraKet



  !
  ! Printing routines for CI results and CPU time.
  !

  ! Public driver routine for printing CI energies, leading configurations,
  ! or complete eigenvectors depending on the requested output level.
  Subroutine PrintCIResults(Flags, Hii, EE, Enuc, SymLabel)
    Implicit None
    Type(ShavittControl),             Intent(In) :: Flags
    Double Precision, Dimension(:),   Intent(In) :: Hii
    Double Precision,                 Intent(In) :: EE, Enuc
    Character(Len=4), Dimension(:),   Intent(In) :: SymLabel
    ! End of dummy parameters.
    Double Precision                             :: Eref, Ecorr

    Eref  = Enuc + EE
    Ecorr = Flags%CIE(Flags%iState)

    If (Flags%PrintLevel >= 0) Then
       Write(Stdout,*)
       Write(Stdout,*) '=========================================================='
       Write(Stdout,*)
       Write(Stdout,*) 'Results of CI calculation:'
       Write(Stdout,*)
       Write(Stdout,fmt=('(" Correlation energy for CI root",I4,":",F20.8," eV")')) Flags%iState, Ecorr
       Write(Stdout,fmt=('(" Total energy:                      ",F20.8," eV")'))   Eref + Ecorr
       Write(Stdout,fmt=('(" Electronic energy:                 ",F20.8," eV")'))   EE   + Ecorr
       Write(Stdout,fmt=('(" Nuclear energy:                    ",F20.8," eV")'))   Enuc
       Write(Stdout,*)
    End If

! MD INTERFACE    If (Flags%PrintLevel >= 1) Then
! MD INTERFACE: [Flags%PrintLevel added]
       Call PrintLeadingConfigurations(Flags, Hii, Flags%CIE, Flags%CIC, Flags%WhoIsWho, SymLabel, Eref, &
                                       Flags%Leading, Flags%TotalRoots, Flags%Multiplicity, Flags%PrintLevel)
       If (Flags%PrintLevel >= 3) Then
          Call PrintEigenSystem(Flags, Flags%CIE, Flags%CIC, SymLabel)
       End If
! MD INTERFACE   End If

    Return
  End Subroutine PrintCIResults



  Subroutine PrintLeadingConfigurations(Flags, Hii, CIE, CIC, iactive, SymLabel, Eshift, climit, maxstate, mult, iout)
    Use Limit, Only: LMCONF
    Implicit None
    Type(ShavittControl)                         :: Flags
    Double Precision, Dimension(:)               :: Hii, CIE
    Double Precision, Dimension(:,:)             :: CIC
    Integer,          Dimension(:)               :: iactive
    Character(Len=4), Dimension(:)               :: SymLabel
    Double Precision                             :: Eshift, climit
    Integer                                      :: maxstate, mult
    ! End of dummy parameters.
    Integer,          Dimension(:),  Allocatable :: iCSF, d
    Character(Len=3), Dimension(0:3)             :: dstring
    Integer                                      :: nactive, nblock, nCSF, mCSF, i, j, k, l, help, irow, dd, n1, n2
    Data dstring / '  -', '  a', '  b', ' ab' /


! MD INTERFACE

    INTEGER :: iout
    INTEGER :: de
    INTEGER :: naco
    INTEGER :: IN2
    INTEGER :: icross

    ! Note warning about bounds-checking below
    COMMON /DYNVA2/ de(40,LMCONF),naco
    COMMON /INOPT2/ IN2(300)

! END MD INTERFACE

    nactive = UBound(iactive,1)
    nCSF    = UBound(CIC,1)
    Allocate(iCSF(nCSF))
    Allocate(d(nactive))

! MD INTERFACE

    icross = IN2(160)
    IF (icross .EQ. 6) THEN

       naco = nactive

       ! Presumably DYNVA2 will be removed when this interface is formalised.
       ! In the meantime these tests need to be kept in sync with the hard-coded 
       ! limits of de.
       If (naco > 40) Then
          Write(Stdout,'(/A)') ' FATAL INTERNAL ERROR: naco exceeds bounds of de in DYNVA2'
          Stop 'PrintLeadingConfigurations'
       End If
       If (nCSF > LMCONF) Then
          Write(Stdout,'(/A)') ' FATAL INTERNAL ERROR: nCSF exceeds bounds of de in DYNVA2'
          Stop 'PrintLeadingConfigurations'
       End If

       k = 1
       Do i=1, UBound(Flags%IndVec,1)
          If (Flags%IndVec(i) > 0) Then
             iCSF(k) = i
             k = k + 1
          End If
       End Do
       Do i=1,nCSF
          help = iCSF(i)-1
          irow = 1
          Do k=nactive, 1, -1
             Do dd=3,0,-1
                If (Flags%DRT(irow)%jd(dd)>0 .And. Flags%DRT(irow)%yd(dd)<=help) Then
                   help = help - Flags%DRT(irow)%yd(dd)
                   irow = Flags%DRT(irow)%jd(dd)
                   d(k) = dd
                   Exit
                End If
             End Do
          End Do

          Do k = 1, nactive
             de(k,i) = d(k)
          End Do

       End Do

    END IF

    IF (iout .LT. 1) RETURN

! END MD INTERFACE


    Do j=1, maxstate
       k = 1
       Do i=1, UBound(Flags%IndVec,1)
          If (Flags%IndVec(i) > 0) Then
             iCSF(k) = i
             k = k + 1
          End If
       End Do

       mCSF = 0
       Do While (mCSF < nCSF)
          k = mCSF+1
          Do i=mCSF+2, nCSF
             If (Abs(CIC(Flags%IndVec(iCSF(i)),j)) > Abs(CIC(Flags%IndVec(iCSF(k)),j)))  k = i
          End Do
          If (mCSF > 0 .And. Abs(CIC(Flags%IndVec(iCSF(k)),j)) < climit)  Exit
          mCSF = mCSF + 1
          If (k /= mCSF) Then
             help       = iCSF(k)
             iCSF(k)    = iCSF(mCSF)
             iCSF(mCSF) = help
          End If
       End Do

       If (UBound(SymLabel,1) > 1) Then
          Write(Stdout,'(/A)') ' -------------------------------------------------------------------------'
          Write(Stdout,'(" State",I3,",  Mult.",I2,",  ",A4,"(",I1,"),  E-E(1)=",F10.6," eV,  E=",F14.6," eV")') &
                j, mult, SymLabel(Flags%SymVec(Flags%IndVec(iCSF(1)))), Flags%SymVec(Flags%IndVec(iCSF(1))),     &
                CIE(j) - CIE(1), CIE(j) + Eshift
          Write(Stdout,'(A/)') ' -------------------------------------------------------------------------'
       Else
          Write(Stdout,'(/A)') ' ---------------------------------------------------------------'
          Write(Stdout,'(" State",I3,",  Mult.",I2,",  E-E(1)=",F10.6," eV,  E=",F14.6," eV")') &
            & j, mult, CIE(j) - CIE(1), CIE(j) + Eshift
          Write(Stdout,'(A/)') ' ---------------------------------------------------------------'
       End If

       Do i=1,mCSF
          help = iCSF(i)-1
          irow = 1
          Do k=nactive, 1, -1
             Do dd=3,0,-1
                If (Flags%DRT(irow)%jd(dd)>0 .And. Flags%DRT(irow)%yd(dd)<=help) Then
                   help = help - Flags%DRT(irow)%yd(dd)
                   irow = Flags%DRT(irow)%jd(dd)
                   d(k) = dd
                   Exit
                End If
             End Do
          End Do
          If (help /= 0) Then
             Write(Stdout,'(A)') ' FATAL INTERNAL ERROR: Inconsistent lexical index.'
             Stop 'PrintLeadingConfigurations'
          End If
          n2 = Min(nactive, 20)
          Write(Stdout,'(F10.2,F14.6,I8,2X,20I4)') &
             Hii(Flags%IndVec(iCSF(i))) + Eshift, CIC(Flags%IndVec(iCSF(i)),j), Flags%IndVec(iCSF(i)), iactive(1:n2)
          Write(Stdout,'(34X,20A4)') dstring(d(1:n2))
          nblock = (nactive - 1) / 20
          Do l=1,nblock
             n1 = 20*l+1
             n2 = Min(20*l+20, nactive)
             Write(Stdout,*)
             Write(Stdout,'(34X,20I4)') iactive(n1:n2)
             Write(Stdout,'(34X,20A4)') dstring(d(n1:n2))
          End Do
          Write(Stdout,*)
       End Do
    End Do

    Deallocate(d)
    Deallocate(iCSF)

    Return
  End Subroutine PrintLeadingConfigurations



  Subroutine PrintEigenSystem(Flags, CIEigenValues, CIEigenVectors, SymLabel)
    Implicit None
    Type(ShavittControl)                         :: Flags
    Double Precision, Dimension(:),   Intent(In) :: CIEigenValues
    Double Precision, Dimension(:,:), Intent(In) :: CIEigenVectors
    Character(Len=4), Dimension(:),   Intent(In) :: SymLabel
    ! End of dummy parameters.
    Integer,          Parameter                  :: nColsPerBlock = 6
    Integer,          Dimension(nColsPerBlock)   :: MaxInd
    Integer                                      :: nci, nroots, nblocks
    Integer                                      :: i, j, k, l, m

    nci     = Ubound(CIEigenVectors,1)
    nroots  = Ubound(CIEigenVectors,2)
    nblocks = (nroots+nColsPerBlock-1) / nColsPerBlock

    Write(Stdout,'(//1X,A/)') 'GUGA-CI eigensystem:'

    Do k=1,nblocks
       l = nColsPerBlock*(k-1)+1
       m = Min(nroots, l+nColsPerBlock-1)
       Write(Stdout,'(" Root",20(I16,4X))') (j,j=l,m)
       If (UBound(SymLabel,1) > 1) Then
          Do j=l,m
             MaxInd(j-l+1:j-l+1) = MaxLoc(Abs(CIEigenVectors(:,j)))
          End Do
          Write(Stdout,'(" Symmetry",20(8X,A,"(",I1,")",5X))') &
            &   (SymLabel(Flags%SymVec(MaxInd(j-l+1))),Flags%SymVec(MaxInd(j-l+1)), j=l,m)
       End If
       Write(Stdout,'(" Energy ",10F20.14)') (CIEigenValues(j), j=l,m)
       Write(Stdout,'(1X,A)') Repeat('-', 20*(m-l+1)+7)
       Do i=1,nci
          Write(Stdout,'(1X,I7,20F20.14)') i, (CIEigenVectors(i,j), j=l,m)
       End Do
       Write(Stdout,'(/)')
    End Do

    Return
  End Subroutine PrintEigenSystem



  Subroutine PrintTime(tuser, tsys, twall)
    Double Precision :: tuser, tsys, twall
    ! End of dummy parameters.
    Double Precision :: tuser1, tsys1, twall1, tcpu

    tuser1 = tuser
    tsys1  = tsys
    twall1 = twall

    Call GetTime(tuser, tsys, twall)

    tuser1 = tuser  - tuser1
    tsys1  = tsys   - tsys1
    twall1 = twall  - twall1
    tcpu   = tuser1 + tsys1

    Write(Stdout,'(/1X,"User   CPU time:",F12.3," s")') tuser1
    Write(Stdout,'( 1X,"System CPU time:",F12.3," s")') tsys1
    Write(Stdout,'( 1X,"Total  CPU time:",F12.3," s")') tcpu
    Write(Stdout,'( 1X,"Wall clock time:",F12.3," s")') twall1

    Return
  End Subroutine PrintTime



  ! Public printing routine.
  Subroutine PrintGenerators(Gen, Ind, NB6, onetwo)
    Implicit None
    Type(Generator), Dimension( :) :: Gen
    Integer,         Dimension(0:) :: Ind
    Integer                        :: NB6, onetwo
    ! End of dummy parameters.
    Integer                        :: h, i, j, k, l, m, ij, kl

    Write(NB6,'(/10X,A/)') "#       VALUE      MU      NU     I   J   K   L"
    Do m=1,UBound(Ind,1)
       Do h=Ind(m-1)+1, Ind(m)
          If (onetwo == 1) Then
             i = 1
             j = 1
             k = Gen(h)%p
             l = Gen(h)%q
          Else
             ij = Gen(h)%p
             kl = Gen(h)%q
             Call Split1(ij, i, j)
             Call Split1(kl, k, l)
          End If
          Write(NB6,'(1X,I10,F12.6,2I8,2X,4I4)') &
            &   h, Gen(h)%value, m, Gen(h)%icol, i, j, k, l
       End Do
    End Do
    Return
  End Subroutine PrintGenerators



  ! Calculate canonical index [i,j] = i*(i-1)/2+j for i>=j.
  Integer Function CanonicalIndex(i, j)
    Implicit None
    Integer, Intent(In) :: i, j
    ! End of dummy parameters.
    Integer :: a, b
    a = Max(i,j)
    b = Min(i,j)
    CanonicalIndex = a * (a-1) / 2 + b
    Return
  End Function CanonicalIndex



  ! Calculate double canonical index [[i,j],[k,l]]
  ! (see comments in CanonicalIndex).
  Integer Function DoubleCanonicalIndex(i, j, k, l)
    Implicit None
    Integer, Intent(in) :: i, j, k, l
    ! End of dummy parameters.
    DoubleCanonicalIndex = CanonicalIndex(CanonicalIndex(i,j), CanonicalIndex(k,l))
    Return
  End Function DoubleCanonicalIndex



  ! Split canonical index into components.
  Subroutine Split1(ij, i, j)
    Implicit None
    Integer, Intent(In)  :: ij
    Integer, Intent(Out) :: i, j
    ! End of dummy parameters.
    j  = NInt(Sqrt(2.D0*ij-0.75D0))
    i  = ij - (j*(j-1))/2
    Return
  End Subroutine Split1



  ! Split double canonical index into components.
  ! (Currently not used.)
  Subroutine Split2(ijkl, i, j, k, l)
    Implicit None
    Integer :: ijkl, i, j, k, l
    ! End of dummy parameters.
    Integer :: ij, kl
    Call Split1(ijkl, ij, kl)
    Call Split1(ij, i, j)
    Call Split1(kl, k, l)
    Return
  End Subroutine Split2


End Module gugautils
