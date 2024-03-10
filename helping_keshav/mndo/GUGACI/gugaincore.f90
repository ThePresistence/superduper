!
! Subroutines for in-core CI calculations.
!
! First version by Michael E. Beck in 1998 at Zurich University.
!
! Final version by Axel Koslowski in 2000-2005 at MPI Muelheim.
!
! Shape-driven algorithm added by Axel Koslowski in 2012-2013.
!

Module gugaincore

  Use davidson
  Use gugaglobal
  Use gugaloops
  Use gugashape
  Use gugasrtgen
  Use gugautils
  Use matutils
  Use property

  Private

  Public :: CalculateGenerators  ! Calculate and store coupling coefficients in memory.
  Public :: InCoreCI             ! Driver routine for solving the CI problem blockwise.
  Public :: InCoreGrad           ! Calculate one- and two-particle (transition) density matrices for the analytic
                                 ! CI gradient and the non-adiabatic coupling vectors.
  Public :: InCoreProp           ! Driver routine for population analysis and spectroscopic property calculation.

Contains

  !
  ! Driver routines to calculate, store and print coupling coefficients.
  !

  ! Public driver routine to calculate and store coupling coefficients in memory.
  Subroutine CalculateGenerators(Flags, Occ, irrep)
    Implicit None
    Type(ShavittControl),           Intent(InOut) :: Flags               ! iRefConf will be modified.
    Double Precision, Dimension(:), Intent(In)    :: Occ                 ! Orbital occupation numbers.
    Integer,                        Intent(In)    :: irrep               ! Requested symmetry or 0 for all generators.
    ! End of dummy parameters.
    Type(GenBufList), Dimension(:), Pointer       :: OneList, TwoList
    Integer                                       :: i
    Double Precision                              :: tuser, tsys, twall

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(//1X,"Calculating generators...")')
       Call GetTime(tuser, tsys, twall)
    End If

    Call AllocGenList(OneList, Flags%TotalCSF)
    Call AllocGenList(TwoList, Flags%TotalCSF)
    Call CalcGenerators(Flags, OneList, TwoList, Occ(Flags%WhoIsWho), irrep)

    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)

    Allocate(Flags%OneIndex(0:Flags%TotalCSF))
    Allocate(Flags%TwoIndex(0:Flags%TotalCSF))
    Flags%OneIndex(0) = 0
    Flags%TwoIndex(0) = 0

    Do i=1, Flags%TotalCSF
       Flags%OneIndex(i) = Flags%OneIndex(i-1) + (OneList(i)%count-1)*genBufSize + OneList(i)%first%count
       Flags%TwoIndex(i) = Flags%TwoIndex(i-1) + (TwoList(i)%count-1)*genBufSize + TwoList(i)%first%count
    End Do

    If (Flags%PrintLevel > -9) Then
       Write(Stdout,'(//1X,"Number of one-electron generators:",I10)') Flags%OneIndex(Flags%TotalCSF)
       Write(Stdout,'(  1X,"Number of two-electron generators:",I10)') Flags%TwoIndex(Flags%TotalCSF)
    End If

    Allocate(Flags%OneGen(Flags%OneIndex(Flags%TotalCSF)))
    Allocate(Flags%TwoGen(Flags%TwoIndex(Flags%TotalCSF)))

    If (Flags%PrintLevel >= 2)  Write(Stdout,'(//1X,"Sorting generators...")')

    Call SortGenList(OneList, Flags%OneGen, Flags%OneIndex, Flags%PrintLevel)
    Call SortGenList(TwoList, Flags%TwoGen, Flags%TwoIndex, Flags%PrintLevel)

    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)

    Deallocate(OneList)
    Deallocate(TwoList)

    If (Flags%PrintLevel >= 6) Then
       Write(Stdout,'(//1X,"One-electron generators:")')
       Call PrintGenerators(Flags%OneGen, Flags%OneIndex, Stdout, 1)
       Write(Stdout,'(//1X,"Two-electron generators:")')
       Call PrintGenerators(Flags%TwoGen, Flags%TwoIndex, Stdout, 2)
    End If

    Return
  End Subroutine CalculateGenerators



  ! During calculation, the coupling coefficients are stored
  ! in linked lists of buffers, one list for each row of the
  ! CI Hamiltonian. The following subroutine allocates the
  ! first buffer of each list.
  !
  Subroutine AllocGenList(list, nGelf)
    Implicit None
    Type(GenBufList), Dimension(:), Pointer :: list
    Integer                                 :: nGelf
    ! End of dummy parameters.
    Integer                                 :: i

    Allocate(list(nGelf))
    Do i=1, UBound(list,1)
       Allocate(list(i)%first)
       Nullify(list(i)%first%next)
       list(i)%first%count = 0
       list(i)%count = 1
    End Do

    Return
  End Subroutine AllocGenList



  ! The following subroutine converts the linked lists of
  ! of generator matrix element blocks that were produced
  ! by calling CalcGenerators() into the OneGen/OneIndex,
  ! TwoGen/TwoIndex structures. OneGen and TwoGen are
  ! one-dimensional arrays containing the generator matrix
  ! elements sorted by row and column of the CI Hamiltonian.
  ! OneIndex(i) and TwoIndex(i) point to the last element
  ! of row i. OneIndex(0) and TwoIndex(0) are equal to 0.
  ! Thus OneIndex(i-1)+1 points to the first one-electron
  ! generator matrix element of row i and the number of
  ! one-electron generator matrix elements in row i is
  ! OneIndex(i)-OneIndex(i-1).
  !
  Subroutine SortGenList(list, gen, ind, prtLevel)
    Implicit None
    Type(GenBufList), Pointer :: list(:)
    Type(Generator),  Pointer :: gen(:)
    Integer,          Pointer :: ind(:)
    Integer                   :: prtLevel
    ! End of dummy parameters.
    Type(GenBuffer),  Pointer :: buffer
    Integer                   :: i, j

    j = 1
    Do i=1, UBound(list,1)
       Do While (Associated(list(i)%first))
          buffer                  => list(i)%first
          gen(j:j+buffer%count-1) =  buffer%gen(1:buffer%count)
          j                       =  j + buffer%count
          list(i)%first           => buffer%next
          Deallocate(buffer)
       End Do
       If (prtLevel < 6) Then
          ! Original sorting routine.
          Call OldQuickSortGen(gen, ind(i-1)+1, ind(i))
          ! Somewhat faster routine in module srtgen.
          ! Call QuickSortGen(gen, ind(i-1)+1, ind(i))
       Else
          ! Unique order for debugging purposes.
          Call DebugSortGen(gen, ind(i-1)+1, ind(i))
       End If
    End Do

    Return
  End Subroutine SortGenList



  ! Sort coupling coefficients of one row according to
  ! increasing column index using the quicksort algorithm.
  !
  Recursive Subroutine OldQuickSortGen(Gen, l, r)
    Implicit None
    Type(Generator), Pointer :: Gen(:)
    Integer                  :: l, r
    ! End of dummy parameters.
    Type(Generator)          :: g
    Integer                  :: i, j, icol

    If (l < r) Then
       i  = l-1
       j  = r
       icol = Gen(r)%icol
       Do
          Do
             i = i + 1
             If (Gen(i)%icol >= icol) Exit
          End Do
          Do
             j = j - 1
             If (Gen(j)%icol <= icol .Or. j <= i) Exit
          End Do
          g      = Gen(i)
          Gen(i) = Gen(j)
          Gen(j) = g
          If (j <= i) Exit
       End Do
       Gen(j) = Gen(i)
       Gen(i) = Gen(r)
       Gen(r) = g
       Call OldQuickSortGen(Gen, l,   i-1)
       Call OldQuickSortGen(Gen, i+1, r)
    End If

    Return
  End Subroutine OldQuickSortGen



  !
  ! Driver routines for the solution of the CI eigenvalue problem.
  !

  ! In-core CI public driver routine.
  Subroutine InCoreCI(Flags, EE, Enuc, ciselt, iciref, CMO, F, TwoIntAO,   &
                      NFirst, NLast, Occ, SymLabel, CIProp, icisym, icall)
    Implicit None
    Type(ShavittControl),             Intent(InOut) :: Flags
    Double Precision,                 Intent(InOut) :: EE
    Double Precision,                 Intent(In)    :: Enuc
    Double Precision,                 Intent(In)    :: ciselt
    Integer,          Dimension(:,:), Intent(InOut) :: iciref
    Double Precision, Dimension(:,:), Intent(In)    :: CMO
    Double Precision, Dimension(:),   Intent(In)    :: F
    Double Precision, Dimension(:,:), Intent(In)    :: TwoIntAO
    Integer,          Dimension(:),   Intent(In)    :: NFirst, NLast
    Double Precision, Dimension(:),   Intent(In)    :: Occ
    Character(Len=4), Dimension(:),   Intent(In)    :: SymLabel
    Double Precision, Dimension(:,:), Intent(Out)   :: CIProp
    Integer,          Dimension(:),   Intent(Out)   :: icisym
    Integer,                          Intent(InOut) :: icall
    ! End of dummy parameters.
    Double Precision, Dimension(:,:), Allocatable   :: OneIntMO
    Double Precision, Dimension(:,:), Allocatable   :: TwoIntMO
    Double Precision, Dimension(:),   Allocatable   :: Hii
    Double Precision, Dimension(:),   Allocatable   :: Etmp
    Double Precision, Dimension(:,:), Allocatable   :: Ctmp
    Integer,          Dimension(:),   Allocatable   :: RootSort
    Integer,          Dimension(1)                  :: imax
    Integer                                         :: nci, ncio, ncione
    Integer                                         :: nciref, isym, i
    Double Precision                                :: tuser, tsys, twall

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(//A)') ' Entering in-core CI section.'
       Call GetTime(tuser, tsys, twall)
    End If

    nci    = UBound(Flags%CIC,1)
    ncio   = UBound(Flags%WhoIsWho,1)
    ncione = ncio * (ncio+1) / 2

    If (Flags%PrintLevel > 5) Then
       Write(Stdout,'(//A/)')' AO-basis two-electron integrals:'
       Call PrintMatrix(TwoIntAO, Stdout)
    End If

    Allocate(OneIntMO(ncio,   ncio))
    Allocate(TwoIntMO(ncione, ncione))
    Allocate(Hii(nci))

    Call ActiveIntegralTransformation(CMO(:,Flags%WhoIsWho), F, TwoIntAO, OneIntMO, TwoIntMO, &
                                      NFirst, NLast, Flags%PrintLevel)

    If (Flags%NumberOfCIRoots > 0) Then
       Call DoInCoreCI(Flags, 0, Flags%TotalCSF, Hii, Flags%CIE, Flags%CIC, &
                       EE, Enuc, OneIntMO, TwoIntMO, SymLabel, icall)
       Call Symmetrize(Flags, Flags%CIE, Flags%CIC, SymLabel, Flags%PrintLevel)
    Else
       Allocate(Etmp(Flags%TotalRoots))
       Allocate(Ctmp(Flags%TotalCSF, Flags%TotalRoots))
       Allocate(RootSort(Flags%TotalRoots))
       RootSort = (/(i,i=1,Flags%TotalRoots)/)
       Ctmp     = 0.D0
       Do isym=1, Flags%NumSym
          If (Flags%NRoots(isym) > 0) Then
             If (Flags%Algorithm == 2) Then
                Call ShapeDrivenGenerators(Flags, isym)
             Else If (Flags%Algorithm == -2) Then
                Call CalculateGenerators(Flags, Occ, isym)
             End If
             Call DoInCoreCI(Flags, Flags%FirstCSF(isym)-1, Flags%LastCSF(isym),  &
                             Hii (Flags%FirstCSF (isym) : Flags%LastCSF (isym)),  &
                             Etmp(Flags%FirstRoot(isym) : Flags%LastRoot(isym)),  &
                             Ctmp(Flags%FirstCSF (isym) : Flags%LastCSF (isym),   &
                                  Flags%FirstRoot(isym) : Flags%LastRoot(isym)),  &
                             EE, Enuc, OneIntMO, TwoIntMO, SymLabel, icall)
             If (Flags%Algorithm == 2) Then
                Call FreeGenList(Flags%OneList)
                Call FreeGenList(Flags%TwoList)
             Else If (Flags%Algorithm == -2) Then
                Deallocate(Flags%OneGen)
                Deallocate(Flags%TwoGen)
                Deallocate(Flags%OneIndex)
                Deallocate(Flags%TwoIndex)
             End If
          End If
       End Do
       Call SortCIE(Etmp, RootSort)
       Flags%CIE = Etmp  (RootSort)
       Flags%CIC = Ctmp(:,RootSort)
       Deallocate(RootSort)
       Deallocate(Ctmp)
       Deallocate(Etmp)
    End If

    Deallocate(OneIntMO)
    Deallocate(TwoIntMO)

    If (Flags%mciref >= 3) Then
       nciref = Flags%nciref
       Call RefCheck(Flags, iciref, Flags%DRT, ciselt, icall)
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

    CIProp        = 0.D0
    i             = Min(Flags%TotalRoots, UBound(CIProp,2))
    CIProp(1,1:i) = Flags%CIE(1:i)
    icisym  (1:i) = Flags%RootSym(1:i)

    EE = EE + Flags%CIE(Flags%iState)

    If (Flags%PrintLevel >= 2) then
       Write(Stdout,'(//A)') ' In-core CI section done.'
       Call PrintTime(tuser, tsys, twall)
    End If

    Return
  End Subroutine InCoreCI



  ! Find eigenvalues and eigenvectors of one block of the CI Hamiltonian.
  Subroutine DoInCoreCI(Flags, icol0, icol1, Hii, CIE, CIC, EE, Enuc, &
                        OneIntMO, TwoIntMO, SymLabel, icall)
    Implicit None
    Type(ShavittControl),             Intent(InOut) :: Flags
    Integer,                          Intent(In)    :: icol0, icol1
    Double Precision, Dimension(:),   Intent(Out)   :: Hii
    Double Precision, Dimension(:),   Intent(Out)   :: CIE
    Double Precision, Dimension(:,:), Intent(Out)   :: CIC
    Double Precision,                 Intent(In)    :: EE
    Double Precision,                 Intent(In)    :: Enuc
    Double Precision, Dimension(:,:), Intent(In)    :: OneIntMO
    Double Precision, Dimension(:,:), Intent(In)    :: TwoIntMO
    Character(Len=4), Dimension(:),   Intent(In)    :: SymLabel
    Integer,                          Intent(InOut) :: icall
    ! End of dummy parameters.
    Double Precision, Dimension(:),   Allocatable   :: Hij
    Integer,          Dimension(:),   Allocatable   :: ifirst
    Integer,          Dimension(:),   Allocatable   :: icol
    Double Precision                                :: dntotal
    Integer                                         :: ndim, n2gen, nnz

    ndim = icol1 - icol0

    If (Flags%Algorithm > 0) Then
       Call ShapeDrivenSparseHamiltonian(Flags, icol0, icol1, OneIntMO, TwoIntMO, Hii, Hij, ifirst, icol)
    Else
       ! Use number of two-electron generator matrix elements as an upper
       ! bound to the number of non-zero CI Hamiltonian matrix elements:
       n2gen = Flags%TwoIndex(icol1) - Flags%TwoIndex(icol0)

       Allocate(Hij(n2gen))
       Allocate(ifirst(ndim+1))
       Allocate(icol(n2gen))

       Call SetupSparseHamiltonian(Flags, Flags%OneIndex(icol0:icol1), Flags%TwoIndex(icol0:icol1), &
                                   icol0, OneIntMO, TwoIntMO, Hii, Hij, ifirst, icol)
    End If

    If (Flags%PrintLevel >= 2) Then
       nnz     = ifirst(ndim+1) - ifirst(1)
       dntotal = 0.5D0 * Dble(ndim) * (Dble(ndim) + 1.0D0)
       Write(Stdout, '(/1X,"CI Hamiltonian statistics (one triangle):")')
       Write(Stdout, '( 1X,"Dimension:                  ",I12)')       ndim
       Write(Stdout, '( 1X,"Number of non-zero elements:",I12)')       nnz
       Write(Stdout, '( 1X,"Sparseness:                 ",F11.2,"%")') 100.D0 * (dntotal - Dble(nnz)) / dntotal
    End If

    If (Flags%PrintLevel >= 4) Then
       Write(Stdout,'(//A)')       ' Hamiltonian matrix (eV) in the Gelfand-Tsetlin basis:'
       Write(Stdout,'(A,F14.8,A)') ' (reference energy:', EE+Enuc, ' eV)'
       Call PrintSparseHamiltonian(Hij, ifirst, icol, ndim)
    End If

    Call DiagonalizerInterface(Flags, Hii, Hij, ifirst, icol, CIE, CIC, SymLabel, icol0, icall)

    Deallocate(icol)
    Deallocate(ifirst)
    Deallocate(Hij)

    Return
  End Subroutine DoInCoreCI



  ! The following subroutine forms the current block of the CI Hamiltonian matrix
  ! (of dimension nCSF) and stores it in sparse form. The structure of the sparse
  ! Hamiltonian is defined as follows (corresponding to CSR sparse matrix format
  ! with some additional requirements concerning the diagonal elements, see below):
  !
  ! The vector Hij contains the matrix elements of the upper triangle of the CI
  ! Hamiltonian stored row-wise. The vector icol contains the column indices of the
  ! corresponding matrix elements. These count the elements of current block only,
  ! in contrast to the column indices stored with the generator matrix elements,
  ! which count the CSFs of all relevant blocks. icol0 is the index of the last
  ! column of the previous block or zero if the current block is the first one.
  ! The vector ifirst (of dimension nCSF+1) contains the index of the first element
  ! of each row in the vector Hij, i.e. Hij(ifirst(k)) is the value of the first
  ! matrix element in row k. Its column index is icol(ifirst(k)). It is required
  ! that the first matrix element stored in each row is the diagonal element,
  ! and it must always be stored even if it is zero. The vector Hii (of dimension
  ! nCSF) contains a copy of the diagonal elements.
  !
  Subroutine SetupSparseHamiltonian(Flags, OneIndex, TwoIndex, icol0, &
                                    OneMO, GamMO, Hii, Hij, ifirst, icol)
    Implicit None
    Type(ShavittControl)             :: Flags
    Integer,          Dimension(0:)  :: OneIndex, TwoIndex
    Integer                          :: icol0
    Double Precision, Dimension(:,:) :: OneMO, GamMO
    Double Precision, Dimension(:)   :: Hii, Hij
    Integer,          Dimension(:)   :: ifirst, icol
    ! End of dummy parameters.
    Integer                          :: nCSF, irow, ic, lastcol, igen, iham, ij, kl, k, l
    Double Precision                 :: tuser, tsys, twall
    !Integer                          :: irowhs, icolhs
    !Double Precision                 :: Hsmall

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(//A)') ' Forming sparse CI Hamiltonian...'
       Call GetTime(tuser, tsys, twall)
    End If

    nCSF = UBound(Hii,1)
    iham = 1

    ! Two-electron contributions:
    ! Loop over all rows.
    Do irow=1, nCSF
       ifirst(irow) = iham
       Hij(iham)    = 0.D0
       icol(iham)   = irow
       lastcol      = irow
       ! The generator matrix elements of each row are sorted according to
       ! increasing column index, so the diagonal element comes first.
       Do igen=TwoIndex(irow-1)+1, TwoIndex(irow)
          ic = Flags%TwoGen(igen)%icol - icol0
          ij = Flags%TwoGen(igen)%p
          kl = Flags%TwoGen(igen)%q
          If (ic == lastcol) Then
             ! Same matrix element
             Hij(iham)  = Hij(iham) + Flags%TwoGen(igen)%value * GamMO(ij,kl)
          Else
             ! Next matrix element.
             iham       = iham + 1
             Hij(iham)  = Flags%TwoGen(igen)%value * GamMO(ij,kl)
             icol(iham) = ic
             lastcol    = ic
          End If
       End Do
       iham = iham + 1
    End Do
    ifirst(nCSF+1) = iham

    ! One-electron contributions:
    Do irow=1, nCSF
       iham = ifirst(irow)
       Do igen=OneIndex(irow-1)+1, OneIndex(irow)
          ic = Flags%OneGen(igen)%icol - icol0
          If (ic > nCSF) Exit
          k  = Flags%OneGen(igen)%p
          l  = Flags%OneGen(igen)%q
          ! Find CI matrix element corresponding to
          ! the current generator matrix element.
          Do While (icol(iham) < ic .And. iham < ifirst(irow+1)-1)
             iham = iham + 1
          End Do
          If (ic == icol(iham)) Then
             Hij(iham) = Hij(iham) + Flags%OneGen(igen)%value * OneMO(k,l)
          End If
       End Do
    End Do

    ! Copy diagonal elements to Hii:
    Do irow=1, nCSF
       Hii(irow) = Hij(ifirst(irow))
    End Do

    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)

    Return
  End Subroutine SetupSparseHamiltonian



  ! Print CI Hamiltonian that is stored in sparse form.
  Subroutine PrintSparseHamiltonian(Hij, ifirst, icol, n)
    Implicit None
    Double Precision, Dimension(:) :: Hij
    Integer,          Dimension(:) :: ifirst, icol
    Integer                        :: n
    ! End of dummy parameters.
    Integer                        :: i, j, k, l

    ! Loop over all rows.
    Do i=1, n
       ! Loop over all blocks of seven columns.
       Do k=ifirst(i), ifirst(i+1)-1, 8
          l = Min(k+7, ifirst(i+1)-1)
          Write(Stdout, '(8X,A,8(I11,4X))')  '|', (icol(j), j=k,l)
          Write(Stdout, '(9A)')  ' -------+-', ('---------------', j=k,l)
          Write(Stdout, '(I6,A3,8F15.8)')  i, '  |', (Hij(j), j=k,l)
          Write(Stdout,*)
       End Do
    End Do

    Return
  End Subroutine PrintSparseHamiltonian



  ! Call the requested lapack or Davidson diagonalizer.
  Subroutine DiagonalizerInterface(Flags, Hii, Hij, ifirst, icol, CIE, CIC, SymLabel, i0, icall)
    Implicit None
    Type(ShavittControl),             Intent(InOut)    :: Flags
    Double Precision, Dimension(:),   Intent(In)       :: Hii, Hij
    Integer,          Dimension(:),   Intent(In)       :: ifirst, icol
    Double Precision, Dimension(:),   Intent(Out)      :: CIE
    Double Precision, Dimension(:,:), Intent(Out)      :: CIC
    Character(Len=4), Dimension(:),   Intent(In)       :: SymLabel
    Integer,                          Intent(In)       :: i0
    Integer,                          Intent(InOut)    :: icall
    ! End of dummy parameters.
    Integer                                            :: cidiag, info, m, nv
    Integer                                            :: lwork, liwork, ndim
    Integer                                            :: nroots, ldcic, ldget
    Integer                                            :: ILAENV, NB
    Double Precision                                   :: abstol, DLAMCH
    Double Precision, Dimension(:,:), Allocatable      :: A, B, Z
    Double Precision, Dimension(:),   Allocatable      :: ASP, work, fv1, e, w
    Integer,          Dimension(:),   Allocatable      :: iwork
    Integer,          Dimension(:),   Allocatable      :: ifail
    Double Precision, Dimension(:),   Pointer          :: dummy1
    Double Precision, Dimension(:,:), Pointer          :: dummy2
    Double Precision, Dimension(1),   Target           :: dummyTarget1
    Double Precision, Dimension(1,1), Target           :: dummyTarget2
    Double Precision                                   :: tuser, tsys, twall
    External                                           :: cpuDavLiu, gpuDavLiu

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(//1X,A)') 'Solving CI eigenvalue problem...'
       Call GetTime(tuser, tsys, twall)
    End If

    dummy1 => dummyTarget1
    dummy2 => dummyTarget2
    ndim   =  UBound(Hii,1)
    nroots =  UBound(CIE,1)
    cidiag =  Flags%WhichDiagonalizer

    ! By default, choose LAPACK or Davidson diagonalizer
    ! depending on the dimension of the CI problem.
    If (cidiag == 0) Then
       If (ndim < 500) Then
          cidiag = 6
       Else
          If (nroots == 1) Then
             cidiag = 10
          Else
             cidiag = 12
          End If
       End If
    End If

    ! Choose CPU or GPU variant of the Davidson diagonalizer.
    If (cidiag == 12) Then
       ! If the Davidson diagonalizer variant is not specified and
       ! a GPU is present, use the GPU variant for in-core CI.
       If (IAbs(Flags%Algorithm) < 3  .And.  Flags%NVCapa >= 20) Then
          cidiag = 14
       Else
          cidiag = 13
       End If
    End If

    ! Leading dimension of the CI coefficient matrix.
    If (nroots > 1) Then
       ldcic = ldget(CIC(1,1), CIC(1,2))
    Else
       ldcic = UBound(CIC,1)
    End If

    ! Miscellaneous diagonalizers introduced (27.6.2000)
    ! Reworked for sparse CI Hamiltonian (11.7.2000)
    Select Case (cidiag)
      Case (1)  ! TDIAG  (EISPACK-BASED), PACKED MATRIX.
        nv = ndim * (ndim+1) / 2
        Allocate(fv1(2*ndim))
        Allocate(ASP(nv))
        Allocate(w(ndim))
        Allocate(Z(ndim,ndim))
        Call SPHamiltonian(ASP, Hij, ifirst, icol, ndim)
        Call TDIAG(ASP, Z, w, fv1, nv, ndim, ndim, nroots)
        CIE = w  (1:nroots)
        CIC = Z(:,1:nroots)
        Deallocate(Z)
        Deallocate(w)
        Deallocate(ASP)
        Deallocate(fv1)
      Case (2)  ! TQL2/TRED2 (EISPACK), SQUARE MATRIX.
        Allocate(w(ndim))
        Allocate(e(ndim))
        Allocate(A(ndim, ndim))
        Allocate(Z(ndim, ndim))
        Call SquareHamiltonian(A, Hij, ifirst, icol, ndim)
        Call TRED2(ndim, ndim, A, w, e, Z)
        Call TQL2 (ndim, ndim, w, e, Z, info)
        CIE = w  (1:nroots)
        CIC = Z(:,1:nroots)
        Deallocate(Z)
        Deallocate(A)
        Deallocate(e)
        Deallocate(w)
        If (info /= 0) Then
          Write(Stdout, '(/A,I5,A)') ' TRED2/TQL2: ', info, '-th eigenvalue failed to converge.'
          icall = -1
          Return
        End If
      Case (3)  ! DSPEV  (LAPACK), PACKED MATRIX.
        Allocate(work(3*ndim))
        Allocate(ASP(ndim*(ndim+1)/2))
        Allocate(w(ndim))
        Allocate(Z(ndim,ndim))
        Call SPHamiltonian(ASP, Hij, ifirst, icol, ndim)
        Call DSPEV('V', 'U', ndim, ASP, w, Z(1,1), ndim, work, info)
        CIE = w  (1:nroots)
        CIC = Z(:,1:nroots)
        Deallocate(Z)
        Deallocate(w)
        Deallocate(ASP)
        Deallocate(work)
        If (info < 0) Then
          Write(Stdout, '(/A,I1,A)') ' DSPEV: Argument #', -info, ' has an illegal value.'
          Stop 'DiagonalizerInterface'
        Else If (info > 0) Then
          Write(Stdout, '(/A)') ' DSPEV: The algorithm failed to converge.'
          icall = -1
          Return
        End If
      Case (4)  ! DSPEVX (LAPACK), PACKED MATRIX.
        Allocate(work(8*ndim))
        Allocate(iwork(5*ndim))
        Allocate(ifail(ndim))
        Allocate(ASP(ndim*(ndim+1)/2))
        Call SPHamiltonian(ASP, Hij, ifirst, icol, ndim)
        abstol = 2.D0 * DLAMCH('S')
        Call DSPEVX('V', 'I', 'U', ndim, ASP, 0.D0, 0.D0, 1, nroots,  &
                    abstol, m, CIE, CIC(1,1), ldcic, work, iwork, ifail, info)
        If (info < 0) Then
          Write(Stdout, '(/A,I2,A)') ' DSPEVX: Argument #', -info, ' has an illegal value.'
          Stop 'DiagonalizerInterface'
        Else If (info > 0) Then
          Write(Stdout, '(/A,I5,A)') ' DSPEVX: ', info, ' eigenvectors failed to converge:'
          Write(Stdout, *) ifail(1:info)
          icall = -1
          Return
        End If
        Deallocate(ASP)
        Deallocate(ifail)
        Deallocate(iwork)
        Deallocate(work)
        If (m < nroots) Then
          Write(Stdout,'(/A)') ' DSPEVX: Less CI roots found than requested.'
          icall = -1
          Return
        End If
      Case (5)  ! DSYEV  (LAPACK), SQUARE MATRIX.
        lwork = 3 * ndim
        If (Flags%PrintLevel >= 2) &
          Write(Stdout, '(/A,I10)') ' DSYEV: Minimal lwork = ', lwork
        NB = ILAENV(1, 'DSYEV', 'VU', ndim, -1, -1, -1)
        If (NB > 1)  lwork = (NB + 2) * ndim
        If (Flags%PrintLevel >= 2) &
          Write(Stdout, '(A,I10)')  ' DSYEV: Actual  lwork = ', lwork
        Allocate(work(lwork))
        Allocate(w(ndim))
        Allocate(Z(ndim,ndim))
        Call SquareHamiltonian(Z, Hij, ifirst, icol, ndim)
        Call DSYEV('V', 'U', ndim, Z(1,1), ndim, w, work, lwork, info)
        If (info == 0 .And. Flags%PrintLevel >= 2) Then
          Write(Stdout, '(A,I10)')   ' DSYEV: Optimal lwork = ', NInt(work(1))
        Else If (info < 0) Then
          Write(Stdout, '(/A,I1,A)') ' DSYEV: Argument #', -info, ' has an illegal value.'
          Stop 'DiagonalizerInterface'
        Else If (info > 0) Then
          Write(Stdout, '(/A)') ' DSYEV: The algorithm failed to converge.'
          icall = -1
          Return
        End If
        CIE =   w(1:nroots)
        CIC = Z(:,1:nroots)
        Deallocate(Z)
        Deallocate(w)
        Deallocate(work)
      Case (6)  ! DSYEVX (LAPACK), SQUARE MATRIX.
        lwork = 8 * ndim
        If (Flags%PrintLevel >= 2) &
          Write(Stdout, '(/A,I10)') ' DSYEVX: Minimal lwork = ', lwork
        NB = ILAENV(1, 'DSYEVX', 'VIU', ndim, -1, -1, -1)
        If (NB > 5)  lwork = (NB + 3) * ndim
        If (Flags%PrintLevel >= 2) &
          Write(Stdout, '(A,I10)')  ' DSYEVX: Actual  lwork = ', lwork
        Allocate(A(ndim, ndim))
        Allocate(work(lwork))
        Allocate(iwork(5*ndim))
        Allocate(ifail(ndim))
        Call SquareHamiltonian(A, Hij, ifirst, icol, ndim)
        abstol = 2.D0 * DLAMCH('S')
        Call DSYEVX('V', 'I', 'U', ndim, A, ndim, 0.D0, 0.D0, 1, nroots, abstol, &
                     m, CIE, CIC(1,1), ldcic, work, lwork, iwork, ifail, info)
        If (info == 0 .And. Flags%PrintLevel >= 2) Then
          Write(Stdout, '(A,I10)') ' DSYEVX: Optimal lwork = ', NInt(work(1))
        Else If (info < 0) Then
          Write(Stdout, '(/A,I2,A)') ' DSYEVX: Argument #', -info, ' has an illegal value.'
          Stop 'DiagonalizerInterface'
        Else If (info > 0) Then
          Write(Stdout, '(/A,I5,A)') ' DSYEVX: ', info, ' eigenvectors failed to converge:'
          Write(Stdout, *) ifail(1:info)
          icall = -1
          Return
        End If
        Deallocate(ifail)
        Deallocate(iwork)
        Deallocate(work)
        Deallocate(A)
        If (m < nroots) Then
          Write(Stdout,'(A)') ' DSYEVX: Less CI roots found than requested.'
          icall = -1
          Return
        End If
      Case (7)  ! DSPEVD (LAPACK), PACKED MATRIX, DIVIDE-CONQUER.
        ! With 64-bit integers lwork and liwork need to be twice as large as expected to avoid a bus error on an SGI Octane.
        lwork  = 2 * (69 + 2 * ndim) * ndim + 2  !!! lg N has been estimated not to exceed 32
        liwork = 2 * (2 + 5 * ndim)
        If (Flags%PrintLevel >= 2) Then
           Write(Stdout, '(/A,I10)') ' DSPEVD:         lwork  = ', lwork
           Write(Stdout, '( A,I10)') ' DSPEVD:         liwork = ', liwork
        End If
        Allocate(work(lwork))
        Allocate(iwork(liwork))
        Allocate(ASP(ndim*(ndim+1)/2))
        Allocate(w(ndim))
        Allocate(Z(ndim,ndim))
        Call SPHamiltonian(ASP, Hij, ifirst, icol, ndim)
        Call DSPEVD('V', 'U', ndim, ASP, w, Z(1,1), ndim, work, lwork, iwork, liwork, info)
        If (info == 0 .And. Flags%PrintLevel >= 2) Then
          Write(Stdout, '(A,I10)')   ' DSPEVD: Optimal lwork  = ', NInt(work(1))
          Write(Stdout, '(A,I10)')   ' DSPEVD: Optimal liwork = ', iwork(1)
        Else If (info < 0) Then
          Write(Stdout, '(/A,I1,A)') ' DSPEVD: Argument #', -info, ' has an illegal value.'
          Stop 'DiagonalizerInterface'
        Else If (info > 0) Then
          Write(Stdout, '(/A)') ' DSPEVD: The algorithm failed to converge.'
          icall = -1
          Return
        End If
        CIE =   w(1:nroots)
        CIC = Z(:,1:nroots)
        Deallocate(Z)
        Deallocate(w)
        Deallocate(ASP)
        Deallocate(iwork)
        Deallocate(work)
      Case (8)  ! DSYEVD (LAPACK), SQUARE MATRIX, DIVIDE-CONQUER.
        ! With 64-bit integers lwork and liwork need to be twice as large as expected to avoid a bus error on an SGI Octane.
        lwork  = 2 * (69 + 3 * ndim) * ndim + 2  !!! lg N has been estimated not to exceed 32
        liwork = 2 * (2 + 5 * ndim)
        If (Flags%PrintLevel >= 2) Then
           Write(Stdout, '(/A,I10)') ' DSYEVD:         lwork  = ', lwork
           Write(Stdout, '( A,I10)') ' DSYEVD:         liwork = ', liwork
        End If
        Allocate(work(lwork))
        Allocate(iwork(liwork))
        Allocate(w(ndim))
        Allocate(Z(ndim,ndim))
        Call SquareHamiltonian(Z, Hij, ifirst, icol, ndim)
        Call DSYEVD('V', 'U', ndim, Z(1,1), ndim, w, work, lwork, iwork, liwork, info)
        If (info == 0 .And. Flags%PrintLevel >= 2) Then
          Write(Stdout, '(A,I10)')   ' DSYEVD: Optimal lwork  = ', NInt(work(1))
          Write(Stdout, '(A,I10)')   ' DSYEVD: Optimal liwork = ', iwork(1)
        Else If (info < 0) Then
          Write(Stdout, '(/A,I1,A)') ' DSYEVD: Argument #', -info, ' has an illegal value.'
          Stop 'DiagonalizerInterface'
        Else If (info > 0) Then
          Write(Stdout, '(/A)') ' DSYEVD: The algorithm failed to converge.'
          icall = -1
          Return
        End If
        CIE =   w(1:nroots)
        CIC = Z(:,1:nroots)
        Deallocate(Z)
        Deallocate(w)
        Deallocate(iwork)
        Deallocate(work)
      Case (9)  ! EWWRSP (EISPACK-BASED), SQUARE MATRIX.
        Allocate(iwork(ndim))
        Allocate(B(ndim, 9))
        Allocate(A(ndim, ndim))
        Allocate(w(ndim))
        Allocate(Z(ndim, nroots))
        Call SquareHamiltonian(A, Hij, ifirst, icol, ndim)
        Call EWWRSP(ndim, nroots, ndim, A, B, iwork, w, Z, info)
        CIE =   w(1:nroots)
        CIC = Z(:,1:nroots)
        Deallocate(Z)
        Deallocate(w)
        Deallocate(A)
        Deallocate(B)
        Deallocate(iwork)
        If (info > 0) Then
           Write(Stdout, '(/A,I5,A)') ' EWWRSP: eigenvalue #', info, ' failed to converge.'
           icall = -1
           Return
        Else If (info < 0) Then
           Write(Stdout, '(/A,I5,A)') ' EWWRSP: eigenvector #', -info, ' failed to converge.'
           icall = -1
           Return
        End If
      Case (10)  ! Davidson's original algorithm, sparse matrix;
        ! The iroot lowest energy eigenvectors are calculated.
        ! E. R. Davidson, J. Comp. Phys. 17 (1975), 87-94.
        Call Davidson1(Flags, Hii, Hij, ifirst, icol, dummy2, dummy2, dummy1, 0, Flags%iRefConf(1:Flags%nRefCSF), &
                       i0, SymLabel, CIE, CIC, nroots, Flags%MinDav, Flags%MaxDav, Flags%KitDav, Flags%QNorm,     &
                       Flags%Algorithm, Stdout, Flags%Printlevel, icall)
        If (icall < 0) Return
      Case (11)  ! Davidson's algorithm with Butscher's and Kammer's modification, sparse matrix;
        ! Only the eigenvector of interest (lroot) is calculated.
        ! W. Butscher, W. E. Kammer, J. Comp. Phys. 20 (1976), 313-325.
        Call Davidson2(Flags, Hii, Hij, ifirst, icol, dummy2, dummy2, dummy1, 0, Flags%iFirstRef(1:Flags%nciref),       &
                       Flags%iLastRef(1:Flags%nciref), Flags%iRefConf(1:Flags%nRefCSF), i0, SymLabel, CIE(1), CIC(:,1), &
                       Flags%LRoot, Flags%MinDav, Flags%MaxDav, Flags%KitDav, Flags%QNorm, Flags%Algorithm, Stdout,     &
                       Flags%Printlevel, icall)
        If (icall < 0) Return
        ! Only one root has been calculated:
        Flags%iState = 1
        Flags%jState = 1
      Case (13)  ! Davidson's algorithm with Liu's modification, sparse matrix, Fortran;
        ! The iroot lowest energy eigenvectors are calculated.
        ! B. Liu, in C. Moler, I. Shavitt, Numerical Algorithms in Chemistry:
        ! Algebraic Methods, LBL-8158 Lawrence Berkeley Laboratory (1978), 49-53.
        Call Davidson3(Flags, Hii, Hij, ifirst, icol, dummy2, dummy2, dummy1, 0, Flags%iRefConf(1:Flags%nRefCSF), &
                       i0, SymLabel, CIE, CIC, nroots, Flags%MinDav, Flags%MaxDav, Flags%KitDav, Flags%QNorm,     &
                       Flags%Algorithm, Stdout, Flags%Printlevel, icall)
        If (icall < 0) Return
      Case (14)  ! Davidson's algorithm with Liu's modification, sparse matrix, C++, GPU;
        Call Davidson4(Flags, Hii, Hij, ifirst, icol, Flags%iRefConf(1:Flags%nRefCSF), i0, SymLabel, &
                       CIE, CIC, nroots, Flags%MinDav, Flags%MaxDav, Flags%KitDav, Flags%QNorm,      &
                       Stdout, Flags%Printlevel, icall, gpuDavLiu)
        If (icall < 0) Return
      Case (15)  ! Davidson's algorithm with Liu's modification, sparse matrix, C++, CPU;
        Call Davidson4(Flags, Hii, Hij, ifirst, icol, Flags%iRefConf(1:Flags%nRefCSF), i0, SymLabel, &
                       CIE, CIC, nroots, Flags%MinDav, Flags%MaxDav, Flags%KitDav, Flags%QNorm,      &
                       Stdout, Flags%Printlevel, icall, cpuDavLiu)
        If (icall < 0) Return
      Case Default
        Write(Stdout, '(/A,I10)') ' Invalid diagonalization method: ', cidiag
        Stop 'gugadriver::DiagonalizerInterface'
    End Select

    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)

    Return
  End Subroutine DiagonalizerInterface



  ! Copy sparse Hamiltonian into square matrix.
  Subroutine SquareHamiltonian(A, Hij, ifirst, icol, n)
    Implicit None
    Double Precision, Dimension(:,:) :: A
    Double Precision, Dimension(:)   :: Hij
    Integer,          Dimension(:)   :: ifirst, icol
    Integer                          :: n
    ! End of dummy parameters.
    Double Precision                 :: h
    Integer                          :: i, j, k

    A = 0.D0
    Do i=1, n
       Do k=ifirst(i), ifirst(i+1)-1
          h = Hij(k)
          j = icol(k)
          A(i,j) = h
          A(j,i) = h
       End Do
    End Do

    Return
  End Subroutine SquareHamiltonian



  ! Copy sparse Hamiltonian into symmetric packed form.
  Subroutine SPHamiltonian(ASP, Hij, ifirst, icol, n)
    Implicit None
    Double Precision, Dimension(:) :: ASP, Hij
    Integer,          Dimension(:) :: ifirst, icol
    Integer                        :: n
    ! End of dummy parameters.
    Integer                        :: i, j, k

    ASP = 0.D0
    Do i=1, n
       Do k=ifirst(i), ifirst(i+1)-1
          j = icol(k)
          ASP(j*(j-1)/2+i) = Hij(k)
       End Do
    End Do

    Return
  End Subroutine SPHamiltonian



  !
  ! Public driver routine to calculate one- and two-particle
  ! density matrices for the analytic CI gradient and NAC
  ! vectors.
  !
  ! Usage of the following subroutine is reasonable
  ! only if the coupling coefficients exist already.
  ! Otherwise it is more economic to use DirectGrad()
  ! since this saves the step of storing and sorting
  ! the generator matrix elements.
  !
  Subroutine InCoreGrad(Flags, OneBodyDen, TwoBodyDen)
    Implicit None
    Type(ShavittControl),             Intent(In)    :: Flags
    Double Precision, Dimension(:,:), Intent(Out)   :: OneBodyDen
    Double Precision, Dimension(:),   Intent(Out)   :: TwoBodyDen
    ! End of dummy parameters.
    Integer                                         :: nactive, nnactive, i, j
    Double Precision                                :: tuser, tsys, twall

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(//1X,A)') 'Density matrices for analytic CI gradient and NAC vectors (in-core algorithm)...'
       Call GetTime(tuser, tsys, twall)
    End If

    i = Flags%iState
    j = Flags%jState

    If (Flags%Algorithm == 1) Then
       Call ShapeDrivenInCoreOneBodyDen(Flags, Flags%CIC(:,i), Flags%CIC(:,j), OneBodyDen)
       Call ShapeDrivenInCoreTwoBodyDen(Flags, Flags%CIC(:,i), Flags%CIC(:,j), TwoBodyDen)
    Else If (Flags%Algorithm == -1) Then
       Call FormOneBodyDen(Flags%OneGen, Flags%OneIndex, Flags%CIC(:,i), Flags%CIC(:,j), OneBodyDen)
       Call FormTwoBodyDen(Flags%TwoGen, Flags%TwoIndex, Flags%CIC(:,i), Flags%CIC(:,j), TwoBodyDen)
    Else
       Write(Stdout,'(/1X,"Invalid usage of subroutine InCoreGrad.")')
       Stop 'InCoreGrad'
    End If

    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)

    If (Flags%PrintLevel >= 5) Then
       nactive  = UBound(Flags%WhoIsWho,1)
       nnactive = nactive * (nactive + 1) / 2
       Write(Stdout,'(/1X,A/)') 'One-particle density matrix:'
       Call PrintMatrix(OneBodyDen, Stdout)
       Write(Stdout,'(/1X,A/)') 'Two-particle density matrix:'
       Call PrintPackedLower(TwoBodyDen, nnactive, Stdout)
    End If

    Return
  End Subroutine InCoreGrad



  !
  ! Public driver routine for natural orbitals,
  ! population analysis and spectroscopic property
  ! calculation (UV, CD).
  !
  ! Usage of the following subroutine is reasonable
  ! only if the coupling coefficients exist already.
  ! Otherwise it is more economic to use DirectProp()
  ! since this saves the step of storing and sorting
  ! the generator matrix elements.
  !
  Subroutine InCoreProp(Flags, NAT, Coord, NFirst, NLast, Core, zs, zp, zd, nsp, nd, CMO, Occ, &
                        SymLabel, CIProp, icisym, PAO, iuvcd, imcd, ipop, inatur, icall)
    Implicit None
    Type(ShavittControl),             Intent(In)    :: Flags
    Integer,          Dimension(:),   Intent(In)    :: NAT
    Double Precision, Dimension(:,:), Intent(In)    :: Coord
    Integer,          Dimension(:),   Intent(In)    :: NFirst, NLast
    Double Precision, Dimension(:),   Intent(In)    :: Core, zs, zp, zd
    Integer,          Dimension(:),   Intent(In)    :: nsp, nd
    Double Precision, Dimension(:,:), Intent(In)    :: CMO
    Double Precision, Dimension(:),   Intent(In)    :: Occ
    Character(Len=4), Dimension(:),   Intent(In)    :: SymLabel
    Double Precision, Dimension(:,:), Intent(Out)   :: CIProp
    Integer,          Dimension(:),   Intent(Out)   :: icisym
    Double Precision, Dimension(:,:), Intent(Out)   :: PAO
    Integer,                          Intent(In)    :: iuvcd, imcd, ipop, inatur
    Integer,                          Intent(InOut) :: icall
    ! End of dummy parameters.
    Integer                                         :: nactive
    Integer                                         :: nPerm, mTrans, nTrans
    Integer                                         :: i, j, k
    Double Precision, Dimension(:,:,:), Allocatable :: Perm, Trans
    Double Precision                                :: tuser, tsys, twall

    If (iuvcd <= 0 .And. ipop <= 0 .And. inatur <=0)  Return

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(///1X,"Natural orbitals, population analysis and properties (in-core algorithm)...")')
       Call GetTime(tuser, tsys, twall)
    End If

    nactive  = UBound(Flags%WhoIsWho,1)
    nPerm    = Flags%TotalRoots
    mTrans   = 0
    nTrans   = 0

    Allocate(Perm(nactive, nactive, nPerm))
    Do i=1, Flags%TotalRoots
       If (Flags%Algorithm == 1) Then
          Call ShapeDrivenInCoreOneBodyDen(Flags, Flags%CIC(:,i), Flags%CIC(:,i), Perm(:,:,i))
       Else If (Flags%Algorithm == -1) Then
          Call FormOneBodyDen(Flags%OneGen, Flags%OneIndex, Flags%CIC(:,i), Flags%CIC(:,i), Perm(:,:,i))
       Else
          Write(Stdout,'(/1X,"Invalid usage of subroutine InCoreProp.")')
          Stop 'InCoreProp'
       End If
    End Do

    If (iuvcd >= 3 .And. Flags%TotalRoots > 1) Then
       mTrans = Flags%TotalRoots-1
       nTrans = Flags%TotalRoots
       Allocate(Trans(nactive, nactive, Flags%TotalRoots*(Flags%TotalRoots-1)/2))
    Else If (iuvcd >= 2 .And. Flags%TotalRoots > 1) Then
       mTrans = 1
       nTrans = Flags%TotalRoots
       Allocate(Trans(nactive, nactive, Flags%TotalRoots-1))
    Else
       Allocate(Trans(1,1,1))
    End If

    If (nTrans > 1) Then
       k = 0
       Do i=1, mTrans
          Do j=i+1, nTrans
             k = k + 1
             If (Flags%Algorithm > 0) Then
                Call ShapeDrivenInCoreOneBodyDen(Flags, Flags%CIC(:,i), Flags%CIC(:,j), Trans(:,:,k))
             Else
                Call FormOneBodyDen(Flags%OneGen, Flags%OneIndex, Flags%CIC(:,i), Flags%CIC(:,j), Trans(:,:,k))
             End If
          End Do
       End Do
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
  End Subroutine InCoreProp



  !
  ! Form one- and two-particle density matrices.
  !
  Subroutine FormOneBodyDen(OneGen, OneIndex, CIV1, CIV2, OneBodyDen)
    Implicit None
    Type(Generator),  Dimension( :)   :: OneGen
    Integer,          Dimension(0:)   :: OneIndex
    Double Precision, Dimension( :)   :: CIV1, CIV2
    Double Precision, Dimension( :,:) :: OneBodyDen
    ! End of dummy parameters.
    Integer                           :: i, j, k, l, m

    OneBodyDen = 0.D0

    Do i=1, UBound(CIV1,1)
       Do m=OneIndex(i-1)+1, OneIndex(i)
          j = OneGen(m)%icol
          k = OneGen(m)%p
          If (i == j) Then
             OneBodyDen(k,k) = OneBodyDen(k,k) + OneGen(m)%value * CIV1(i) * CIV2(i)
          Else
             l = OneGen(m)%q
             OneBodyDen(k,l) = OneBodyDen(k,l) + OneGen(m)%value * CIV1(i) * CIV2(j)
             OneBodyDen(l,k) = OneBodyDen(l,k) + OneGen(m)%value * CIV1(j) * CIV2(i)
          End If
       End Do
    End Do

    Return
  End Subroutine FormOneBodyDen



  Subroutine FormOneBodyDenDiag(OneGen, OneIndex, CIV1, CIV2, OneBodyDenDiag)
    Implicit None
    Type(Generator),  Dimension( :) :: OneGen
    Integer,          Dimension(0:) :: OneIndex
    Double Precision, Dimension( :) :: CIV1, CIV2
    Double Precision, Dimension( :) :: OneBodyDenDiag
    ! End of dummy parameters.
    Double Precision                :: temp
    Integer                         :: i, j, k, l, m

    OneBodyDenDiag = 0.D0

    Do i=1, UBound(CIV1,1)
       Do m=OneIndex(i-1)+1, OneIndex(i)
          k = OneGen(m)%p
          l = OneGen(m)%q
          If (k == l) Then
             j    = OneGen(m)%icol
             temp = CIV1(i) * CIV2(j)
             If (i /= j)  temp = temp + CIV1(j) * CIV2(i)
             OneBodyDenDiag(k) = OneBodyDenDiag(k) + OneGen(m)%value * temp
          End If
       End Do
    End Do

    Return
  End Subroutine FormOneBodyDenDiag



  !
  ! Form the two-body transition density between two CI vectors.
  ! m is the index of the current generator matrix element in
  ! linear storage. i and j are the corresponding configuration
  ! state functions. ij, kl and ijkl are the canonical and doubly
  ! canonical MO indices, respectively.
  !
  Subroutine FormTwoBodyDen(TwoGen, TwoIndex, CIV1, CIV2, TwoBodyDen)
    Implicit None
    Type(Generator),  Dimension( :) :: TwoGen
    Integer,          Dimension(0:) :: TwoIndex
    Double Precision, Dimension( :) :: CIV1, CIV2
    Double Precision, Dimension( :) :: TwoBodyDen
    ! End of dummy parameters.
    Double Precision                :: temp
    Integer                         :: i, j, m, ij, kl, ijkl

    TwoBodyDen = 0.D0

    Do i=1, UBound(CIV1,1)
       Do m=TwoIndex(i-1)+1, TwoIndex(i)
          j    = TwoGen(m)%icol
          ij   = TwoGen(m)%p
          kl   = TwoGen(m)%q
          ijkl = CanonicalIndex(ij,kl)
          temp = CIV1(i) * CIV2(j)
          If (i /= j)  temp = temp + CIV1(j) * CIV2(i)
          TwoBodyDen(ijkl) = TwoBodyDen(ijkl) + temp * TwoGen(m)%value
       End Do
    End Do

    Return
  End Subroutine FormTwoBodyDen


End Module gugaincore
