!
! This module provides driver routines for the shape-driven GUGA-CI approach.
! These call the subroutines that calculate upper and lower partial loops (in
! modules gugaupper and gugalower), that combine upper and lower partial loops
! (in modules gugashape1 through gugashape2) to calculate the generator matrix
! elements which are either stored (in-core CI) or used "on the fly" (direct CI),
! and that use the generator matrix elements stored in-core (this module).
!
! Written by Axel Koslowski in 2012-2013 at MPI Muelheim.
!

Module gugashape

Use gugaglobal
Use gugalower
Use gugashape1
Use gugashape2
Use gugasrtgen
Use gugaupper
Use gugautils

Implicit None

Private

Public :: CalculatePartialLoops
Public :: InitSymCombinations
Public :: InitAllCombinations
Public :: CleanupCombinations
Public :: ShapeDrivenGenerators
Public :: ShapeDrivenSparseHamiltonian
Public :: ShapeDrivenInCoreOneBodyDen
Public :: ShapeDrivenInCoreTwoBodyDen
Public :: FreeGenList

Contains

  ! This is the driver routine for partial loop calculation.
  !
  Subroutine CalculatePartialLoops(Flags, Occ)
    Implicit None
    Type(ShavittControl), Intent(InOut) :: Flags
    Double Precision,     Intent(In)    :: Occ(:)  ! Orbital occupation numbers.
    ! End of dummy parameters.
    Type(LoopBufLists)                  :: bufLists

    Call FindBorderLevel(Flags, Flags%kBorder)
    Call InitLoopBufLists(bufLists)
    Call CalcUpperPartialLoops(Flags, bufLists, Occ(Flags%WhoIsWho))
    Call CalcLowerPartialLoops(Flags, bufLists, Occ(Flags%WhoIsWho))
    If (Flags%PrintLevel >= 1) Then
       Write(Stdout,'(/1X,"Total number of upper partial loops:    ",I10)') bufLists%ncountU
       Write(Stdout,'( 1X,"Total number of lower partial loops:    ",I10)') bufLists%ncountL
       Write(Stdout,'( 1X,"Total number of complete one-body loops:",I10)') bufLists%ncount1
       Write(Stdout,'( 1X,"Total number of complete two-body loops:",I10)') bufLists%ncount2
    End If
    Allocate(Flags%UpperPart(bufLists%ncountU))
    Allocate(Flags%LowerPart(bufLists%ncountL))
    Allocate(Flags%OneLoop(bufLists%ncount1))
    Allocate(Flags%TwoLoop(bufLists%ncount2))
    If (Flags%PrintLevel >= 2) Write(Stdout,'(/1X,"Sorting upper partial loops...")')
    Call SortPartialLoops(Flags, bufLists%firstU, Flags%UpperPart, Flags%iFirstUp, Flags%iLastUp, numUpTypes)
    If (Flags%PrintLevel >= 2) Write(Stdout,'(/1X,"Sorting lower partial loops...")')
    Call SortPartialLoops(Flags, bufLists%firstL, Flags%LowerPart, Flags%iFirstLo, Flags%iLastLo, numLoTypes)
    If (Flags%PrintLevel >= 2) Write(Stdout,'(/1X,"Sorting complete one-body loops...")')
    Call SortCompleteLoops(Flags, bufLists%first1, Flags%OneLoop, Flags%iFirst1, Flags%iDiag1, Flags%iOff1, Flags%iLast1, Flags%NumSym)
    If (Flags%PrintLevel >= 2) Write(Stdout,'(/1X,"Sorting complete two-body loops...")')
    Call SortCompleteLoops(Flags, bufLists%first2, Flags%TwoLoop, Flags%iFirst2, Flags%iDiag2, Flags%iOff2, Flags%iLast2, Flags%NumSym)
    Call AssembleUpperWalks(Flags)
    If (Flags%PrintLevel >= 5) Then
       Write(Stdout,'(///1X,"Upper partial loops:")')
       Call PrintPartialLoops(Flags%UpperPart, Flags%iFirstUp, Flags%iLastUp, numUpTypes)
       Write(Stdout,'(///1X,"Lower partial loops:")')
       Call PrintPartialLoops(Flags%LowerPart, Flags%iFirstLo, Flags%iLastLo, numLoTypes)
       Call PrintOneBodyLoops(Flags)
       Call PrintTwoBodyLoops(Flags)
    End If
  End Subroutine CalculatePartialLoops



  ! This is the analogon of driver routine CalculateGenerators in module gugaincore.
  ! The generator matrix elements will be identical, but with a different storage format.
  ! This subroutine will create linked lists of buffers containing the generator matrix
  ! elements, one list for each diagonal of the CI Hamiltonian (in Flags%OneList and
  ! Flags%TwoList).
  !
  Subroutine ShapeDrivenGenerators(Flags, Irrep)
    Implicit None
    Type(ShavittControl)         :: Flags
    Integer                      :: Irrep
    ! End of dummy parameters.
    Type(Generator),     Pointer :: TmpGen(:)
    Integer,             Pointer :: TmpIndex(:)
    Integer                      :: i
    Double Precision             :: tuser, tsys, twall

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(/1X,"Storing generator matrix elements...")')
       Call GetTime(tuser, tsys, twall)
    End If

    ! Allocate and initialize Flags%OneList and Flags%TwoList.
    Allocate(Flags%OneList(0:Flags%TotalCSF-1))
    Allocate(Flags%TwoList(0:Flags%TotalCSF-1))
    Do i=0, Flags%TotalCSF-1
       Nullify(Flags%OneList(i)%first)
       Nullify(Flags%TwoList(i)%first)
       Flags%OneList(i)%count = 0
       Flags%TwoList(i)%count = 0
    End Do

    If (Irrep == 0) Then
       Call StoreGenerators(Flags)
    Else
       Call StoreGeneratorsIrrep(Flags, Irrep)
    End If

    ! Print timing if requested.
    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)

    ! Debug print of generator matrix elements if requested.
    If (Flags%PrintLevel >= 6) Then
       Write(Stdout,'(///1X,"One-electron generator matrix elements:")')
       Call ReorgGenerators(Flags%OneList, TmpGen, TmpIndex)
       Call SortGenerators(TmpGen, TmpIndex)
       Call PrintGenerators(TmpGen, TmpIndex, Stdout, 1)
       Deallocate(TmpGen, TmpIndex)
       Write(Stdout,'(///1X,"Two-electron generator matrix elements:")')
       Call ReorgGenerators(Flags%TwoList, TmpGen, TmpIndex)
       Call ArrangeIJKL(TmpGen)
       Call SortGenerators(TmpGen, TmpIndex)
       Call PrintGenerators(TmpGen, TmpIndex, Stdout, 2)
       Deallocate(TmpGen, TmpIndex)
    End If
  End Subroutine ShapeDrivenGenerators



  ! The following subroutine forms the current block or several adjacent blocks
  ! (of total dimension icol1-icol0) of the CI Hamiltonian matrix (of dimension
  ! nci) and stores it in sparse form. The structure of the sparse Hamiltonian
  ! is defined as follows (corresponding to CSR sparse matrix format with some
  ! additional requirements concerning the diagonal elements, see below):
  !
  ! The vector Hij contains the matrix elements of the upper triangle of the CI
  ! Hamiltonian stored row-wise. The vector icol contains the column indices of
  ! the corresponding matrix elements. These count the elements of the current
  ! block only, in contrast to the column indices stored with the generator
  ! matrix elements, which count the CSFs of all relevant blocks. icol0 is the
  ! index of the last column of the previous block or zero if the current block
  ! is the first one. The vector ifirst (of dimension nci+1) contains the index
  ! of the first element of each row in the vector Hij, i.e. Hij(ifirst(k)) is
  ! the value of the first matrix element in row k. Its column index is
  ! icol(ifirst(k)). It is required that the first matrix element stored in each
  ! row is the diagonal element, and it must always be stored even if it is zero.
  ! The vector Hii (of dimension nci) contains a copy of the diagonal elements.
  !
  Subroutine ShapeDrivenSparseHamiltonian(Flags, icol0, icol1, OneMO, GamMO, &
                                          Hii, Hij, ifirst, icol)
    Implicit None
    Type(ShavittControl)             :: Flags
    Integer                          :: icol0, icol1
    Double Precision, Dimension(:,:) :: OneMO, GamMO
    Double Precision, Dimension(:)   :: Hii
    Double Precision, Allocatable    :: Hij(:)
    Integer,          Allocatable    :: ifirst(:), icol(:)
    ! End of dummy parameters.
    Integer,          Allocatable    :: iwork(:)
    Type(GenBuffer),  Pointer        :: bufptr
    Integer                          :: nci, nnz, ioff, ic, ir, i, j, k, l
    Double Precision                 :: tuser, tsys, twall
    Double Precision, Parameter      :: Zero = 0.0D0

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(/1X,"Forming sparse CI Hamiltonian...")')
       Call GetTime(tuser, tsys, twall)
    End If

    ! Dimension of the CI problem.
    nci = UBound(Hii,1)
    ! Total number of non-zero matrix elements.
    nnz = nci

    ! Index of first matrix element in each row (1:nci).
    ! Temporarily, number of non-zero matrix elements in each row.
    Allocate(ifirst(nci+1))
    ! Temporarily, last column (1:nci) in each row (1:nci).
    Allocate(icol(nci))

    ! Initialize.
    Do i=1, nci
       ifirst(i) = 1
       icol(i)   = i
    End Do
    ifirst(nci+1) = 0

    ! Count non-zero matrix elements in each row and in total.
    Do ioff=0, Flags%TotalCSF-1
       bufptr => Flags%OneList(ioff)%first
       Do While (Associated(bufptr))
          Do i=1, bufptr%count
             ic = bufptr%gen(i)%icol
             ir = ic - ioff
             If (Flags%SymVec(ir) == Flags%SymVec(ic) .And. ic > icol0 .And. ic <= icol1) Then
                ic = ic - icol0
                ir = ir - icol0
                If (ic > icol(ir)) Then
                   icol(ir)   = ic
                   ifirst(ir) = ifirst(ir) + 1
                   nnz        = nnz        + 1
                End If
             End If
          End Do
          bufptr => bufptr%next
       End Do
       bufptr => Flags%TwoList(ioff)%first
       Do While (Associated(bufptr))
          Do i=1, bufptr%count
             ic = bufptr%gen(i)%icol
             If (ic > icol0 .And. ic <= icol1) Then
                ic = ic - icol0
                ir = ic - ioff
                If (ic > icol(ir)) Then
                   icol(ir)   = ic
                   ifirst(ir) = ifirst(ir) + 1
                   nnz        = nnz        + 1
                End If
             End If
          End Do
          bufptr => bufptr%next
       End Do
    End Do

    ! Deallocate temporary vector.
    Deallocate(icol)

    ! Index of first element in each row.
    i = 1
    Do ir=1, nci+1
       j          = ifirst(ir)
       ifirst(ir) = i
       i          = i + j
    End Do

    ! Working copy.
    Allocate(iwork(nci))
    iwork(:) = ifirst(1:nci)

    ! Allocate and initialize final buffers.
    Allocate(Hij(nnz))
    Allocate(icol(nnz))
    Do ir=1, nci
       i       = ifirst(ir)
       Hij(i)  = Zero
       icol(i) = ir
    End Do

    ! Assemble matrix elements of the CI Hamiltonian.
    Do ioff=0, Flags%TotalCSF-1
       bufptr => Flags%OneList(ioff)%first
       Do While (Associated(bufptr))
          Do i=1, bufptr%count
             ic = bufptr%gen(i)%icol
             ir = ic - ioff
             If (Flags%SymVec(ir) == Flags%SymVec(ic) .And. ic > icol0 .And. ic <= icol1) Then
                ic = ic - icol0
                ir = ir - icol0
                j  = iwork(ir)
                If (ic > icol(j)) Then
                   j         = j + 1
                   iwork(ir) = j
                   Hij(j)    = Zero
                   icol(j)   = ic
                End If
                k      = bufptr%gen(i)%p
                l      = bufptr%gen(i)%q
                Hij(j) = Hij(j) + bufptr%gen(i)%value * OneMO(k,l)
             End If
          End Do
          bufptr => bufptr%next
       End Do
       bufptr => Flags%TwoList(ioff)%first
       Do While (Associated(bufptr))
          Do i=1, bufptr%count
             ic = bufptr%gen(i)%icol
             If (ic > icol0 .And. ic <= icol1) Then
                ic = ic - icol0
                ir = ic - ioff
                j  = iwork(ir)
                If (ic > icol(j)) Then
                   j         = j + 1
                   iwork(ir) = j
                   Hij(j)    = Zero
                   icol(j)   = ic
                End If
                k      = bufptr%gen(i)%p
                l      = bufptr%gen(i)%q
                Hij(j) = Hij(j) + bufptr%gen(i)%value * GamMO(k,l)
             End If
          End Do
          bufptr => bufptr%next
       End Do
    End Do

    ! Copy diagonal elements of the CI Hamiltonian.
    Do ir=1, nci
       i       = ifirst(ir)
       Hii(ir) = Hij(i)
    End Do

    ! Print timing of this subroutine if requested.
    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)
  End Subroutine ShapeDrivenSparseHamiltonian



  ! This subroutine calculates the one-particle density matrix from
  ! the generator matrix elements formed by ShapeDrivenGenerators.
  !
  Subroutine ShapeDrivenInCoreOneBodyDen(Flags, CIV1, CIV2, OneBodyDen)
    Implicit None
    Type(ShavittControl)     :: Flags
    Double Precision         :: CIV1(:)
    Double Precision         :: CIV2(:)
    Double Precision         :: OneBodyDen(:,:)
    ! End of dummy parameters.
    Integer                  :: ioff, i, j, k, l, m
    Type(GenBuffer), Pointer :: ptr

    ! Initialize.
    OneBodyDen = 0.D0

    ! Diagonal elements.
    ptr => Flags%OneList(0)%first
    Do While (Associated(ptr))
       Do m=1, ptr%count
          i               = ptr%gen(m)%icol
          k               = ptr%gen(m)%p
          OneBodyDen(k,k) = OneBodyDen(k,k) + ptr%gen(m)%value * CIV1(i) * CIV2(i)
       End Do
       ptr => ptr%next
    End Do

    ! Off-diagonal elements.
    Do ioff=1, UBound(Flags%OneList,1)
       ptr => Flags%OneList(ioff)%first
       Do While (Associated(ptr))
          Do m=1, ptr%count
             j               = ptr%gen(m)%icol
             i               = j - ioff
             k               = ptr%gen(m)%p
             l               = ptr%gen(m)%q
             OneBodyDen(k,l) = OneBodyDen(k,l) + ptr%gen(m)%value * CIV1(i) * CIV2(j)
             OneBodyDen(l,k) = OneBodyDen(l,k) + ptr%gen(m)%value * CIV1(j) * CIV2(i)
          End Do
          ptr => ptr%next
       End Do
    End Do
  End Subroutine ShapeDrivenInCoreOneBodyDen



  ! This subroutine calculates the two-particle density matrix from
  ! the generator matrix elements formed by ShapeDrivenGenerators.
  !
  Subroutine ShapeDrivenInCoreTwoBodyDen(Flags, CIV1, CIV2, TwoBodyDen)
    Implicit None
    Type(ShavittControl)     :: Flags
    Double Precision         :: CIV1(:)
    Double Precision         :: CIV2(:)
    Double Precision         :: TwoBodyDen(:)
    ! End of dummy parameters.
    Integer                  :: ioff, ijkl, ij, kl, i, j, m
    Type(GenBuffer), Pointer :: ptr

    ! Initialize.
    TwoBodyDen = 0.D0

    ! Diagonal elements.
    ptr => Flags%TwoList(0)%first
    Do While (Associated(ptr))
       Do m=1, ptr%count
          i                = ptr%gen(m)%icol
          ij               = ptr%gen(m)%p
          kl               = ptr%gen(m)%q
          ijkl             = CanonicalIndex(ij,kl)
          TwoBodyDen(ijkl) = TwoBodyDen(ijkl) + ptr%gen(m)%value * CIV1(i) * CIV2(i)
       End Do
       ptr => ptr%next
    End Do

    ! Off-diagonal elements.
    Do ioff=1, UBound(Flags%TwoList,1)
       ptr => Flags%TwoList(ioff)%first
       Do While (Associated(ptr))
          Do m=1, ptr%count
             j                = ptr%gen(m)%icol
             i                = j - ioff
             ij               = ptr%gen(m)%p
             kl               = ptr%gen(m)%q
             ijkl             = CanonicalIndex(ij,kl)
             TwoBodyDen(ijkl) = TwoBodyDen(ijkl) + ptr%gen(m)%value * (CIV1(i) * CIV2(j) + CIV1(j) * CIV2(i))
          End Do
          ptr => ptr%next
       End Do
    End Do
  End Subroutine ShapeDrivenInCoreTwoBodyDen



  ! This subroutine frees the memory needed
  ! to store the generator matrix elements.
  !
  Subroutine FreeGenList(List)
    Implicit None
    Type(GenBufList), Pointer :: List(:)
    ! End of dummy parameters.
    Type(GenBuffer),  Pointer :: ptr
    Integer                   :: i

    Do i=LBound(List,1), UBound(List,1)
       Do While (Associated(List(i)%first))
          ptr        => List(i)%first
          List(i)%first => ptr%next
          Deallocate(ptr)
       End Do
    End Do

    Deallocate(List)
  End Subroutine FreeGenList



  ! This subroutine finds the orbital level at which
  ! to separate upper and lower part of the DRT.
  !
  Subroutine FindBorderLevel(Flags, kBorder)
    Implicit None
    Type(ShavittControl), Intent(In) :: Flags
    Integer                          :: kBorder
    ! End of dummy parameters.
    Integer :: numUpperWalks(0:Flags%DRT(1)%k)
    Integer :: numLowerWalks(0:Flags%DRT(1)%k)
    Integer :: i, k

    numUpperWalks(:) = 0
    numLowerWalks(:) = 0

    ! Calculate number of upper and lower walks from each orbital level.
    Do i=1, UBound(Flags%DRT,1)
       k                = Flags%DRT(i)%k
       numUpperWalks(k) = numUpperWalks(k) + Flags%DRT(i)%xu
       numLowerWalks(k) = numLowerWalks(k) + Flags%DRT(i)%yd(4)
    End Do

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(/1X,"Number of upper and lower walks from each orbital level:")')
       Write(Stdout,'(/1X,"Level   Upper walks   Lower walks")')
       Do k=Flags%DRT(1)%k, 0, -1
          Write(Stdout,'(1X,I4,2I12)') k, numUpperWalks(k), numLowerWalks(k)
       End Do
    End If

    ! Find orbital level with the smallest number of upper plus lower walks. ???
    ! kBorder = 1
    ! Do k=2, Flags%DRT(1)%k-1
    !    If (numUpperWalks(k) + numLowerWalks(k) < numUpperWalks(kBorder) + numLowerWalks(kBorder))  kBorder = k
    ! End Do

    ! Find orbital level with the smallest difference between upper and lower walk numbers.
    kBorder = 1
    Do k=2, Flags%DRT(1)%k-1
       If (IAbs(numUpperWalks(k) - numLowerWalks(k)) < IAbs(numUpperWalks(kBorder) - numLowerWalks(kBorder)))  kBorder = k
    End Do

    If (Flags%PrintLevel > 0) Then
       Write(Stdout,'(/1X,"Orbital level separating upper and lower part of the DRT:",I4)') kBorder
    End If
  End Subroutine FindBorderLevel



  ! This subroutine nullifies the pointers to the first buffers
  ! in the linked lists of buffers and resets the counters.
  !
  Subroutine InitLoopBufLists(bufLists)
    Implicit None
    Type(LoopBufLists) :: bufLists
    ! End of dummy parameters.

    Nullify(bufLists%firstU)
    Nullify(bufLists%firstL)
    Nullify(bufLists%first1)
    Nullify(bufLists%first2)
    bufLists%ncountU = 0
    bufLists%ncountL = 0
    bufLists%ncount1 = 0
    bufLists%ncount2 = 0
  End Subroutine InitLoopBufLists



  ! Input to this subroutine is the linked list of temporary partial loop
  ! buffers (first). Final buffer Part to which all partial loops will be
  ! copied must be allocated before calling this subroutine. Partial loops
  ! are arranged by type and partial loops of each type are sorted by
  ! ascending ibra, iket, symmetry, i.e. iProd(sybra, syket), and itip,
  ! with priority decreasing in this order. On return, elements of iFirst
  ! and iLast point to the first and last partial loop with the
  ! corresponding type.
  !
  Subroutine SortPartialLoops(Flags, first, Part, iFirst, iLast, numTypes)
    Implicit None
    Type(ShavittControl)             :: Flags
    Type(PartialLoopBuffer), Pointer :: first
    Type(PartialLoop),       Pointer :: Part(:)
    Integer                          :: numTypes
    Integer                          :: iFirst(numTypes)
    Integer                          :: iLast(numTypes)
    ! End of dummy parameters.
    Type(PartialLoopBuffer), Pointer :: bufptr
    Integer                          :: itype, isym, idx, i
    Double Precision                 :: tuser, tsys, twall

    ! Get initial timings if requested.
    If (Flags%PrintLevel >= 2)  Call GetTime(tuser, tsys, twall)

    ! Count partial loops of each type.
    ! The numbers are stored in iFirst(:).
    iFirst(:) =  0
    bufptr    => first
    Do While (Associated(bufptr))
       Do i=1, bufptr%ncount
          itype         = bufptr%part(i)%itype
          iFirst(itype) = iFirst(itype) + 1
       End Do
       bufptr => bufptr%next
    End Do

    ! Initialize iLast(itype) so that it points to the last
    ! element of type itype currently stored in buffer Part.
    idx = 0
    Do itype=1, numTypes
       iLast(itype)  = idx
       idx           = idx + iFirst(itype)
       iFirst(itype) = iLast(itype) + 1
    End Do

    ! Copy partial loops from linked list of temporary buffers
    ! to final buffer Part and deallocate all temporary buffers.
    Do While (Associated(first))
       bufptr => first
       Do i=1, bufptr%ncount
          itype        = bufptr%part(i)%itype
          idx          = iLast(itype) + 1
          iLast(itype) = idx
          Part(idx)    = bufptr%part(i)
       End Do
       first => bufptr%next
       Deallocate(bufptr)
    End Do

    ! Sort partial loops of each type according to the criteria
    ! described above.
    Do itype=1, numTypes
       Call QuickSortPartialLoops(Part, iFirst(itype), iLast(itype))
    End Do

    ! Print timing of this subroutine if requested.
    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)
  End Subroutine SortPartialLoops



  ! This subroutine sorts the partial loops in the buffer between (and
  ! including) the indices iFirst and iLast using the Quicksort algorithm
  ! by ascending difference lexket-lexbra, ascending ibra, ascending iket,
  ! and ascending itip, in this order.
  !
  Recursive Subroutine QuickSortPartialLoops(Buffer, iFirst, iLast)
    Implicit None
    Type(PartialLoop), Pointer :: Buffer(:)
    Integer                    :: iFirst, iLast
    ! End of dummy parameters.
    Type(PartialLoop)          :: Temp
    Integer                    :: iLeft, iRight, iCenter, i, j

    ! Quicksort is a recursive algorithm that works as follows.
    ! In every call to sort a group of elements, one of the elements
    ! is picked as the so-called partitioning element, and all elements
    ! are rearranged in two subgroups, one in which all elements are
    ! less than or equal to the partitioning element and one in which
    ! all elements are greater than or equal to the partitioning element.
    ! The strategy concerning elements equal to the partitioning element
    ! is subtle. According to Numerical Recipes, Sedgewick recommends
    ! treating them as incorrectly positioned.

    ! As the partitioning element we choose the median of the leftmost,
    ! the rightmost, and the center element. This way, an already sorted
    ! group of elements is no longer the worst case. Moreover, this gives
    ! us the chance to put so-called sentinels at the ends of the array
    ! at which the innermost loops are guaranteed to stop. This simplifies
    ! the innermost loops and may speed up the algorithm. However, due
    ! to the complicated comparison of two partial loops the benefit
    ! of this technique may be small. Exchanging two partial loops may
    ! also be relatively expensive due to the size of the corresponding
    ! data structures (64 bytes). Therefore, the overhead of recursive
    ! subroutine calls will be comparatively small so we will use this
    ! technique for the benefit of simpler code.

    ! Quicksort will be applied if there are at least eight elements,
    ! otherwise a straight insertion sort will be done.
    If (iLast - iFirst > 6) Then
       iCenter = (iFirst + iLast) / 2
       iLeft   = iFirst + 1
       iRight  = iLast

       ! Move the iCenter element to position iFirst+1.
       Temp             = Buffer(iFirst+1)
       Buffer(iFirst+1) = Buffer(iCenter)
       Buffer(iCenter)  = Temp

       ! Rearrange so that Buffer(iFirst) <= Buffer(iFirst+1) <= Buffer(iLast).
       If (ComparePartialLoops(Buffer(iFirst), Buffer(iLast)) > 0) Then
          Temp             = Buffer(iFirst)
          Buffer(iFirst)   = Buffer(iLast)
          Buffer(iLast)    = Temp
       End If
       If (ComparePartialLoops(Buffer(iFirst+1), Buffer(iLast)) > 0) Then
          Temp             = Buffer(iFirst+1)
          Buffer(iFirst+1) = Buffer(iLast)
          Buffer(iLast)    = Temp
       End If
       If (ComparePartialLoops(Buffer(iFirst), Buffer(iFirst+1)) > 0) Then
          Temp             = Buffer(iFirst)
          Buffer(iFirst)   = Buffer(iFirst+1)
          Buffer(iFirst+1) = Temp
       End If

       ! Buffer(iFirst+1) is the partitioning element.
       Do
          ! Find the leftmost element that is not less than Buffer(iFirst+1).
          Do
             iLeft = iLeft + 1
             If (ComparePartialLoops(Buffer(iLeft), Buffer(iFirst+1)) >= 0) Exit
          End Do

          ! Find the rightmost element that is not greater than Buffer(iFirst+1).
          Do
             iRight = iRight - 1
             If (ComparePartialLoops(Buffer(iRight), Buffer(iFirst+1)) <= 0) Exit
          End Do

          ! Stop if pointers crossed.
          If (iLeft > iRight) Exit

          ! We have found two incorrectly positioned elements.
          ! Switch them to get them into the correct order.
          Temp           = Buffer(iLeft)
          Buffer(iLeft)  = Buffer(iRight)
          Buffer(iRight) = Temp
       End Do

       ! At this point, the situation may look like this:
       ! iFirst   iFirst+1   iFirst+2   iFirst+3   iRight   iLeft   iLast-1   iLast
       !    3         5          1          2         4       6        7        8
       !
       ! Or, iLeft and iRight previously stopped at the same element (at iLast-3):
       ! iFirst   iFirst+1   iFirst+2    iRight   iLast-3   iLeft   iLast-1   iLast
       !    3         5          1          2         5       6        7        8

       Temp             = Buffer(iFirst+1)
       Buffer(iFirst+1) = Buffer(iRight)
       Buffer(iRight)   = Temp

       Call QuickSortPartialLoops(Buffer, iFirst, iRight-1)
       Call QuickSortPartialLoops(Buffer, iLeft,  iLast)
    Else
       ! Straight insertion sort.
       ! At every iteration of the outermost loop, insert
       ! the current element (at position i) correctly into
       ! the previously sorted sequence at iFirst ... i-1.
       Do i=iFirst+1, iLast
          Temp = Buffer(i)
          j    = i
          ! Workaround avoiding ifort compiler bug...
          ! Do While (j > iFirst .And. ComparePartialLoops(Temp, Buffer(j-1)) < 0)
          Do While (j > iFirst)
             If (ComparePartialLoops(Temp, Buffer(j-1)) >= 0) Exit
             Buffer(j) = Buffer(j-1)
             j         = j - 1
          End Do
          Buffer(j) = Temp
       End Do
    End If
  End Subroutine QuickSortPartialLoops



  ! This function compares two partial loops P and Q considering
  ! ibra, iket, iProd(sybra, syket), and itip, and returns an integer
  !    < 0 if P <  Q,
  !   == 0 if P == Q, and
  !    > 0 if P >  Q.
  !
  Integer Function ComparePartialLoops(P, Q)
    Type(PartialLoop) :: P, Q
    ! End of dummy parameters.
    Integer           :: idiff

    idiff = P%ibra - Q%ibra

    If (idiff /= 0) Then
       ComparePartialLoops = idiff
       Return
    End If

    idiff = P%iket - Q%iket

    If (idiff /= 0) Then
       ComparePartialLoops = idiff
       Return
    End If

    idiff = iProd(P%sybra, P%syket) - iProd(Q%sybra, Q%syket)

    If (idiff /= 0) Then
       ComparePartialLoops = idiff
       Return
    End If

    ComparePartialLoops = P%itip - Q%itip
    Return
  End Function ComparePartialLoops



  ! Input to this subroutine is the linked list of temporary complete loop
  ! buffers (first). Final buffer Loop to which all complete loops will be
  ! copied must be allocated before calling this subroutine. Loops are
  ! sorted by ascending symmetry, i.e. iProd(sybra, syket), ascending
  ! difference lexket-lexbra, ascending itop, and ascending ibot, in this
  ! order. On return, elements of iFirst, iDiag, iOff, and iLast point to
  ! the corresponding complete loops of the appropriate symmetry.
  !
  Subroutine SortCompleteLoops(Flags, first, Loop, iFirst, iDiag, iOff, iLast, numSym)
    Implicit None
    Type(ShavittControl)              :: Flags
    Type(CompleteLoopBuffer), Pointer :: first
    Type(CompleteLoop),       Pointer :: Loop(:)
    Integer                           :: numSym
    Integer                           :: iFirst(numSym)
    Integer                           :: iDiag(numSym)
    Integer                           :: iOff(numSym)
    Integer                           :: iLast(numSym)
    ! End of dummy parameters.
    Type(CompleteLoopBuffer), Pointer :: bufptr
    Integer                           :: isym, idx, i
    Double Precision                  :: tuser, tsys, twall

    ! Get initial timings if requested.
    If (Flags%PrintLevel >= 2)  Call GetTime(tuser, tsys, twall)

    ! Count complete loops by symmetry, i.e. iProd(sybra, syket).
    ! The total numbers are stored in iFirst(:), the numbers of
    ! diagonal loops in iDiag(:).
    iFirst(:) =  0
    iDiag(:)  =  0
    bufptr    => first
    Do While (Associated(bufptr))
       Do i=1, bufptr%ncount
          isym         = iProd(bufptr%loop(i)%sybra, bufptr%loop(i)%syket)
          iFirst(isym) = iFirst(isym) + 1

          If (bufptr%loop(i)%lexbra == bufptr%loop(i)%lexket) Then
             iDiag(isym) = iDiag(isym) + 1
          End If
       End Do
       bufptr => bufptr%next
    End Do

    ! Initialize iFirst(:), iDiag(:), and iOff(:) such that
    ! their elements point to the correct loop in buffer Loop.
    ! Initialize iLast(:) such that its elements may be used
    ! as pointers to the last loop of the corresponding symmetry
    ! currently stored in buffer Loop during copying (next step).
    idx = 0
    Do isym=1, numSym
       iLast(isym)  = idx
       idx          = idx + iFirst(isym)
       iFirst(isym) = iLast(isym)  + 1
       iOff(isym)   = iFirst(isym) + iDiag(isym)
       iDiag(isym)  = iOff(isym)   - 1
    End Do

    ! Copy complete loops from linked list of buffers to
    ! final buffer Loop and deallocate all temporary buffers.
    Do While (Associated(first))
       bufptr => first
       Do i=1, bufptr%ncount
          isym        = iProd(bufptr%loop(i)%sybra, bufptr%loop(i)%syket)
          idx         = iLast(isym) + 1
          iLast(isym) = idx
          Loop(idx)   = bufptr%loop(i)
       End Do
       first => bufptr%next
       Deallocate(bufptr)
    End Do

    ! Sort complete loops of each symmetry
    ! according to the criteria described above.
    Do isym=1, numSym
       Call QuickSortCompleteLoops(Loop, iFirst(isym), iLast(isym))
    End Do

    ! Print timing of this subroutine if requested.
    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)
  End Subroutine SortCompleteLoops



  ! This subroutine sorts the complete loops in the buffer between (and
  ! including) the indices iFirst and iLast using the Quicksort algorithm
  ! by ascending difference lexket-lexbra, ascending itop, and ascending
  ! ibot, in this order.
  !
  Recursive Subroutine QuickSortCompleteLoops(Buffer, iFirst, iLast)
    Implicit None
    Type(CompleteLoop), Pointer :: Buffer(:)
    Integer                     :: iFirst, iLast
    ! End of dummy parameters.
    Type(CompleteLoop)          :: Temp
    Integer                     :: iLeft, iRight, iCenter, i, j

    ! Quicksort is a recursive algorithm that works as follows.
    ! In every call to sort a group of elements, one of the elements
    ! is picked as the so-called partitioning element, and all elements
    ! are rearranged in two subgroups, one in which all elements are
    ! less than or equal to the partitioning element and one in which
    ! all elements are greater than or equal to the partitioning element.
    ! The strategy concerning elements equal to the partitioning element
    ! is subtle. According to Numerical Recipes, Sedgewick recommends
    ! treating them as incorrectly positioned.

    ! As the partitioning element we choose the median of the leftmost,
    ! the rightmost, and the center element. This way, an already sorted
    ! group of elements is no longer the worst case. Moreover, this gives
    ! us the chance to put so-called sentinels at the ends of the array
    ! at which the innermost loops are guaranteed to stop. This simplifies
    ! the innermost loops and may speed up the algorithm. However, due
    ! to the complicated comparison of two complete loops the benefit
    ! of this technique may be small. Exchanging two complete loops may
    ! also be relatively expensive due to the size of the corresponding
    ! data structures (40 bytes). Therefore, the overhead of recursive
    ! subroutine calls will be comparatively small so we will use this
    ! technique for the benefit of simpler code.

    ! Quicksort will be applied if there are at least eight elements,
    ! otherwise a straight insertion sort will be done.
    If (iLast - iFirst > 6) Then
       iCenter = (iFirst + iLast) / 2
       iLeft   = iFirst + 1
       iRight  = iLast

       ! Move the iCenter element to position iFirst+1.
       Temp             = Buffer(iFirst+1)
       Buffer(iFirst+1) = Buffer(iCenter)
       Buffer(iCenter)  = Temp

       ! Rearrange so that Buffer(iFirst) <= Buffer(iFirst+1) <= Buffer(iLast).
       If (CompareCompleteLoops(Buffer(iFirst), Buffer(iLast)) > 0) Then
          Temp             = Buffer(iFirst)
          Buffer(iFirst)   = Buffer(iLast)
          Buffer(iLast)    = Temp
       End If
       If (CompareCompleteLoops(Buffer(iFirst+1), Buffer(iLast)) > 0) Then
          Temp             = Buffer(iFirst+1)
          Buffer(iFirst+1) = Buffer(iLast)
          Buffer(iLast)    = Temp
       End If
       If (CompareCompleteLoops(Buffer(iFirst), Buffer(iFirst+1)) > 0) Then
          Temp             = Buffer(iFirst)
          Buffer(iFirst)   = Buffer(iFirst+1)
          Buffer(iFirst+1) = Temp
       End If

       ! Buffer(iFirst+1) is the partitioning element.
       Do
          ! Find the leftmost element that is not less than Buffer(iFirst+1).
          Do
             iLeft = iLeft + 1
             If (CompareCompleteLoops(Buffer(iLeft), Buffer(iFirst+1)) >= 0) Exit
          End Do

          ! Find the rightmost element that is not greater than Buffer(iFirst+1).
          Do
             iRight = iRight - 1
             If (CompareCompleteLoops(Buffer(iRight), Buffer(iFirst+1)) <= 0) Exit
          End Do

          ! Stop if pointers crossed.
          If (iLeft > iRight) Exit

          ! We have found two incorrectly positioned elements.
          ! Switch them to get them into the correct order.
          Temp           = Buffer(iLeft)
          Buffer(iLeft)  = Buffer(iRight)
          Buffer(iRight) = Temp
       End Do

       ! At this point, the situation may look like this:
       ! iFirst   iFirst+1   iFirst+2   iFirst+3   iRight   iLeft   iLast-1   iLast
       !    3         5          1          2         4       6        7        8
       !
       ! Or, iLeft and iRight previously stopped at the same element (at iLast-3):
       ! iFirst   iFirst+1   iFirst+2    iRight   iLast-3   iLeft   iLast-1   iLast
       !    3         5          1          2         5       6        7        8

       Temp             = Buffer(iFirst+1)
       Buffer(iFirst+1) = Buffer(iRight)
       Buffer(iRight)   = Temp

       Call QuickSortCompleteLoops(Buffer, iFirst, iRight-1)
       Call QuickSortCompleteLoops(Buffer, iLeft,  iLast)
    Else
       ! Straight insertion sort.
       ! At every iteration of the outermost loop, insert
       ! the current element (at position i) correctly into
       ! the previously sorted sequence at iFirst ... i-1.
       Do i=iFirst+1, iLast
          Temp = Buffer(i)
          j    = i
          ! Workaround avoiding ifort compiler bug...
          ! Do While (j > iFirst .And. CompareCompleteLoops(Temp, Buffer(j-1)) < 0)
          Do While (j > iFirst)
             If (CompareCompleteLoops(Temp, Buffer(j-1)) >= 0) Exit
             Buffer(j) = Buffer(j-1)
             j         = j - 1
          End Do
          Buffer(j) = Temp
       End Do
    End If
  End Subroutine QuickSortCompleteLoops



  ! This function compares two complete loops P and Q considering
  ! lexket-lexbra, itop, and ibot, and returns an integer
  !    < 0 if P <  Q,
  !   == 0 if P == Q, and
  !    > 0 if P >  Q.
  ! Symmetry is not considered because complete loops are already
  ! arranged with respect to symmetry.
  !
  Integer Function CompareCompleteLoops(P, Q)
    Type(CompleteLoop) :: P, Q
    ! End of dummy parameters.
    Integer            :: idiff

    ! Uncomment to consider symmetry as major sorting criterion.
    !
    ! idiff = iProd(P%sybra, P%syket) - iProd(Q%sybra, Q%syket)
    !
    ! If (idiff /= 0) Then
    !    CompareCompleteLoops = idiff
    !    Return
    ! End If

    idiff = (P%lexket - P%lexbra) - (Q%lexket - Q%lexbra)

    If (idiff /= 0) Then
       CompareCompleteLoops = idiff
       Return
    End If

    idiff = P%itop - Q%itop

    If (idiff /= 0) Then
       CompareCompleteLoops = idiff
       Return
    End If

    CompareCompleteLoops = P%ibot - Q%ibot
    Return
  End Function CompareCompleteLoops



  ! Driver routine storing upper walk lexical indices of all rows in the DRT
  ! but the last (the tail of the Shavitt graph) into array Flags%lexUpper.
  !
  Subroutine AssembleUpperWalks(Flags)
    Implicit None
    Type(ShavittControl) :: Flags
    ! End of dummy parameters.
    Integer              :: nrows, n, i, j, k, m, d, y

    ! Initialize all DRT(i)%iu such that they point to the first
    ! upper walk lexical index of the _previous_ row (row i-1).
    ! After calling CollectUpperWalks all DRT(i)%iu will
    ! point to the first upper walk lexical index of the
    ! corresponding row (row i).
    nrows           = UBound(Flags%DRT,1)
    Flags%DRT(1)%iu = 1
    n               = 1
    Do i=1, nrows-1
       Flags%DRT(i+1)%iu = n
       n                 = n + Flags%DRT(i)%xu
    End Do

    If (Flags%PrintLevel >= 2) &
       Write(Stdout,'(/1X,"Total number of upper walks of all rows but the last:",I10)') n-1

    Allocate(Flags%lexUpper(n-1))
    Call CollectUpperWalks(Flags, 1, 0)

    If (Flags%PrintLevel >= 6) Then
       Write(Stdout,'(/1X,"Debug print of upper walks:")')
       Write(Stdout,'(/1X,"Row   Upper walk lexical indices")')
       Do i=1, nrows-1
          Do j=Flags%DRT(i)%iu, Flags%DRT(i+1)%iu-1, 12
             m = Min(j+11, Flags%DRT(i+1)%iu-1)
             Write(Stdout,'(1X,I3,12I10)') i, (Flags%lexUpper(k), k=j,m)
          End Do
       End Do
    End If
  End Subroutine AssembleUpperWalks



  ! Worker routine called by AssembleUpperWalks.
  !
  Recursive Subroutine CollectUpperWalks(Flags, irow, y)
    Implicit None
    Type(ShavittControl) :: Flags
    Integer              :: irow, y
    ! End of dummy parameters.
    Integer              :: d, j

    If (Flags%DRT(irow)%k == 0) Return
    Flags%lexUpper(Flags%DRT(irow+1)%iu) = y
    Flags%DRT(irow+1)%iu = Flags%DRT(irow+1)%iu + 1

    Do d=0,3
       j = Flags%DRT(irow)%jd(d)
       If (j /= 0) Then
          Call CollectUpperWalks(Flags, j, y+Flags%DRT(irow)%yd(d))
       End If
    End Do
  End Subroutine CollectUpperWalks



  ! Print table of all (upper or lower) partial loops for debugging purposes.
  !
  Subroutine PrintPartialLoops(Part, iFirst, iLast, numTypes)
    Implicit None
    Type(ShavittControl)       :: Flags
    Type(PartialLoop), Pointer :: Part(:)
    Integer                    :: numTypes
    Integer                    :: iFirst(numTypes)
    Integer                    :: iLast(numTypes)
    ! End of dummy parameters.
    Integer                    :: itype, isym, i

    Write(Stdout,'(/1X,"index   type   ibra   iket   sym   itip    singlet     triplet    lexbra   lexket   sybra   syket   i   j   k   l")')
    Do itype=1, numTypes
       If (iFirst(itype) <= iLast(itype)) Then
          Do i=iFirst(itype), iLast(itype)
             isym = iProd(Part(i)%sybra, Part(i)%syket)
             Write(Stdout,'(1X,I5,I6,2I7,I6,I7,1X,2F12.6,I8,I9,I7,I8,I6,3I4)') &
               i, Part(i)%itype, Part(i)%ibra, Part(i)%iket, isym, Part(i)%itip, Part(i)%singlet, Part(i)%triplet, &
               Part(i)%lexbra, Part(i)%lexket, Part(i)%sybra, Part(i)%syket, Part(i)%i, Part(i)%j, Part(i)%k, Part(i)%l
          End Do
       End If
    End Do
  End Subroutine PrintPartialLoops



  ! Print table of all complete one-particle loops in the
  ! upper or lower part of the DRT for debugging purposes.
  !
  Subroutine PrintOneBodyLoops(Flags)
    Implicit None
    Type(ShavittControl) :: Flags
    ! End of dummy parameters.
    Integer              :: isym, ioff, i

    Write(Stdout,'(///1X,"Complete one-body loops:")')
    Write(Stdout,'(/1X,"index   sym   ioff   itop   ibot     value     lexbra   lexket   sybra   syket   i   j")')
    Do isym=1, Flags%NumSym
       If (Flags%iFirst1(isym) <= Flags%iLast1(isym)) Then
          Do i=Flags%iFirst1(isym), Flags%iLast1(isym)
             ioff = Flags%OneLoop(i)%lexket - Flags%OneLoop(i)%lexbra
             Write(Stdout,'(1X,2I5,3I7,F13.6,I8,I9,I7,I8,I6,I4)') i, isym, ioff, &
               Flags%OneLoop(i)%itop,   Flags%OneLoop(i)%ibot,   Flags%OneLoop(i)%value, &
               Flags%OneLoop(i)%lexbra, Flags%OneLoop(i)%lexket, &
               Flags%OneLoop(i)%sybra,  Flags%OneLoop(i)%syket,  &
               Flags%OneLoop(i)%p,      Flags%OneLoop(i)%q
          End Do
       End If
    End Do
  End Subroutine PrintOneBodyLoops



  ! Print table of all complete two-particle loops in the
  ! upper or lower part of the DRT for debugging purposes.
  !
  Subroutine PrintTwoBodyLoops(Flags)
    Implicit None
    Type(ShavittControl) :: Flags
    ! End of dummy parameters.
    Integer              :: isym, ioff, i, j, k, l, m

    Write(Stdout,'(///1X,"Complete two-body loops:")')
    Write(Stdout,'(/1X,"index   sym   ioff   itop   ibot     value     lexbra   lexket   sybra   syket   i   j   k   l")')
    Do isym=1, Flags%NumSym
       If (Flags%iFirst2(isym) <= Flags%iLast2(isym)) Then
          Do m=Flags%iFirst2(isym), Flags%iLast2(isym)
             Call Split1(Flags%TwoLoop(m)%p, i, j)
             Call Split1(Flags%TwoLoop(m)%q, k, l)
             ioff = Flags%TwoLoop(m)%lexket - Flags%TwoLoop(m)%lexbra
             Write(Stdout,'(1X,2I5,3I7,F13.6,I8,I9,I7,I8,I6,3I4)') m, isym, ioff, &
               Flags%TwoLoop(m)%itop,   Flags%TwoLoop(m)%ibot,   Flags%TwoLoop(m)%value, &
               Flags%TwoLoop(m)%lexbra, Flags%TwoLoop(m)%lexket, &
               Flags%TwoLoop(m)%sybra,  Flags%TwoLoop(m)%syket,  i, j, k, l
          End Do
       End If
    End Do
  End Subroutine PrintTwoBodyLoops



  ! Print matrix with numbers of upper partial loops for each combination of
  ! rows at the orbital level separating upper and lower part of the DRT.
  ! Currently not used, may be called when needed.
  !
  Subroutine PrintStatistics(Flags, bufLists, kBorder)
    Implicit None
    Type(ShavittControl) :: Flags
    Type(LoopBufLists)   :: bufLists
    Integer              :: kBorder
    ! End of dummy parameters.
    Integer, Allocatable :: matrix(:,:)
    Integer              :: nrows, last, ioffset, ibra, iket, i, j
    Type(PartialLoopBuffer), Pointer :: ptr

    nrows = 0
    Do i=1, Size(Flags%DRT)
       If (Flags%DRT(i)%k == kBorder) Then
          nrows = nrows + 1
          last  = i
       End If
    End Do
    ioffset = last - nrows

    Allocate(matrix(nrows, nrows))
    Do i=1, numUpTypes
       matrix(:,:) = 0
       ptr => bufLists%firstU
       Do While (Associated(ptr))
          Do j=1, ptr%ncount
             If (ptr%part(j)%itype == i) Then
                ibra = ptr%part(j)%ibra - ioffset
                iket = ptr%part(j)%iket - ioffset
                matrix(ibra, iket) = matrix(ibra, iket) + 1
             End If
          End Do
          ptr => ptr%next
       End Do
       ! Print matrix.
       Write(Stdout,'(1X,/"Upper partial loop connection points of type ",I2,":"/)') i
       Write(Stdout,'(9X,15I8)') (j+ioffset,j=1,nrows)
       Do j=1,nrows
          Write(Stdout,'(1X,16I8)') j+ioffset, matrix(j, 1:nrows)
       End Do
    End Do
    Deallocate(matrix)
  End Subroutine PrintStatistics



  ! Driver routine assembling lists of combinations of upper
  ! and lower partial loops with the same symmetry.
  !
  Subroutine InitSymCombinations(Flags)
    Implicit None
    Type(ShavittControl) :: Flags
    ! End of dummy parameters.
    Double Precision     :: tuser, tsys, twall

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(/1X,"Finding partial loop combinations of equal symmetry...")')
       Call GetTime(tuser, tsys, twall)
    End If

    Call CountSymCombinations(Flags)
    Allocate(Flags%CombiListOne      (Flags%NumOne))
    Allocate(Flags%CombiListSym13ijkl(Flags%NumSym13ijkl))
    Allocate(Flags%CombiListSym22ijkl(Flags%NumSym22ijkl))
    Allocate(Flags%CombiListSym22ikjl(Flags%NumSym22ikjl))
    Allocate(Flags%CombiListSym22iill(Flags%NumSym22iill))
    Allocate(Flags%CombiListSym22ilil(Flags%NumSym22ilil))
    Allocate(Flags%CombiListSym22half(Flags%NumSym22half))
    Allocate(Flags%CombiListSym31ijkl(Flags%NumSym31ijkl))
    Call StoreSymCombinations(Flags)

    ! Debug print of combinations.
    If (Flags%PrintLevel >= 5) Then
       Write(Stdout,'(///1X,"Flags%CombiListOne:")')
       Call PrintCombinations(Flags%CombiListOne,       Flags%NumOne)
       Write(Stdout,'(///1X,"Flags%CombiListSym13ijkl:")')
       Call PrintCombinations(Flags%CombiListSym13ijkl, Flags%NumSym13ijkl)
       Write(Stdout,'(///1X,"Flags%CombiListSym22ijkl:")')
       Call PrintCombinations(Flags%CombiListSym22ijkl, Flags%NumSym22ijkl)
       Write(Stdout,'(///1X,"Flags%CombiListSym22ikjl:")')
       Call PrintCombinations(Flags%CombiListSym22ikjl, Flags%NumSym22ikjl)
       Write(Stdout,'(///1X,"Flags%CombiListSym22iill:")')
       Call PrintCombinations(Flags%CombiListSym22iill, Flags%NumSym22iill)
       Write(Stdout,'(///1X,"Flags%CombiListSym22ilil:")')
       Call PrintCombinations(Flags%CombiListSym22ilil, Flags%NumSym22ilil)
       Write(Stdout,'(///1X,"Flags%CombiListSym22half:")')
       Call PrintCombinations(Flags%CombiListSym22half, Flags%NumSym22half)
       Write(Stdout,'(///1X,"Flags%CombiListSym31ijkl:")')
       Call PrintCombinations(Flags%CombiListSym31ijkl, Flags%NumSym31ijkl)
    End If

    ! Print timing if requested.
    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)
  End Subroutine InitSymCombinations



  ! Driver routine assembling lists of all combinations of upper
  ! and lower partial loops regardless of their symmetry.
  !
  Subroutine InitAllCombinations(Flags)
    Implicit None
    Type(ShavittControl) :: Flags
    ! End of dummy parameters.
    Double Precision     :: tuser, tsys, twall

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(/1X,"Finding all valid partial loop combinations...")')
       Call GetTime(tuser, tsys, twall)
    End If

    Call CountAllCombinations(Flags)
    Allocate(Flags%CombiListAll13ijkl(Flags%NumAll13ijkl))
    Allocate(Flags%CombiListAll22ijkl(Flags%NumAll22ijkl))
    Allocate(Flags%CombiListAll22ikjl(Flags%NumAll22ikjl))
    Allocate(Flags%CombiListAll22iill(Flags%NumAll22iill))
    Allocate(Flags%CombiListAll22ilil(Flags%NumAll22ilil))
    Allocate(Flags%CombiListAll22half(Flags%NumAll22half))
    Allocate(Flags%CombiListAll31ijkl(Flags%NumAll31ijkl))
    Call StoreAllCombinations(Flags)

    ! Debug print of combinations.
    If (Flags%PrintLevel >= 5) Then
       Write(Stdout,'(///1X,"Flags%CombiListAll13ijkl:")')
       Call PrintCombinations(Flags%CombiListAll13ijkl, Flags%NumAll13ijkl)
       Write(Stdout,'(///1X,"Flags%CombiListAll22ijkl:")')
       Call PrintCombinations(Flags%CombiListAll22ijkl, Flags%NumAll22ijkl)
       Write(Stdout,'(///1X,"Flags%CombiListAll22ikjl:")')
       Call PrintCombinations(Flags%CombiListAll22ikjl, Flags%NumAll22ikjl)
       Write(Stdout,'(///1X,"Flags%CombiListAll22iill:")')
       Call PrintCombinations(Flags%CombiListAll22iill, Flags%NumAll22iill)
       Write(Stdout,'(///1X,"Flags%CombiListAll22ilil:")')
       Call PrintCombinations(Flags%CombiListAll22ilil, Flags%NumAll22ilil)
       Write(Stdout,'(///1X,"Flags%CombiListAll22half:")')
       Call PrintCombinations(Flags%CombiListAll22half, Flags%NumAll22half)
       Write(Stdout,'(///1X,"Flags%CombiListAll31ijkl:")')
       Call PrintCombinations(Flags%CombiListAll31ijkl, Flags%NumAll31ijkl)
    End If

    ! Print timing if requested.
    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)
  End Subroutine InitAllCombinations



  ! Free memory needed for lists of combinations
  ! of upper and lower partial loops.
  !
  Subroutine CleanupCombinations(Flags)
    Implicit None
    Type(ShavittControl) :: Flags
    ! End of dummy parameters.

    If (Associated(Flags%CombiListOne))        Deallocate(Flags%CombiListOne)
    If (Associated(Flags%CombiListSym13ijkl))  Deallocate(Flags%CombiListSym13ijkl)
    If (Associated(Flags%CombiListSym22ijkl))  Deallocate(Flags%CombiListSym22ijkl)
    If (Associated(Flags%CombiListSym22ikjl))  Deallocate(Flags%CombiListSym22ikjl)
    If (Associated(Flags%CombiListSym22iill))  Deallocate(Flags%CombiListSym22iill)
    If (Associated(Flags%CombiListSym22ilil))  Deallocate(Flags%CombiListSym22ilil)
    If (Associated(Flags%CombiListSym22half))  Deallocate(Flags%CombiListSym22half)
    If (Associated(Flags%CombiListSym31ijkl))  Deallocate(Flags%CombiListSym31ijkl)
    If (Associated(Flags%CombiListAll13ijkl))  Deallocate(Flags%CombiListAll13ijkl)
    If (Associated(Flags%CombiListAll22ijkl))  Deallocate(Flags%CombiListAll22ijkl)
    If (Associated(Flags%CombiListAll22ikjl))  Deallocate(Flags%CombiListAll22ikjl)
    If (Associated(Flags%CombiListAll22iill))  Deallocate(Flags%CombiListAll22iill)
    If (Associated(Flags%CombiListAll22ilil))  Deallocate(Flags%CombiListAll22ilil)
    If (Associated(Flags%CombiListAll22half))  Deallocate(Flags%CombiListAll22half)
    If (Associated(Flags%CombiListAll31ijkl))  Deallocate(Flags%CombiListAll31ijkl)
  End Subroutine CleanupCombinations



  ! Count combinations of upper and lower
  ! partial loops of equal symmetry.
  !
  Subroutine CountSymCombinations(Flags)
    Implicit None
    Type(ShavittControl) :: Flags
    ! End of dummy parameters.

    Flags%NumOne       = 0
    Flags%NumSym13ijkl = 0
    Flags%NumSym22ijkl = 0
    Flags%NumSym22ikjl = 0
    Flags%NumSym22iill = 0
    Flags%NumSym22ilil = 0
    Flags%NumSym22half = 0
    Flags%NumSym31ijkl = 0

    Call CountAllCombis(Flags, Flags%NumOne,       iUpperRt,   iLowerRb)
    Call CountSymCombis(Flags, Flags%NumSym13ijkl, iUpperRt,   iLowerR)
    Call CountSymCombis(Flags, Flags%NumSym22iill, iUpperW,    iLowerW)
    Call CountSymCombis(Flags, Flags%NumSym22ijkl, iUpperW,    iLowerOneR)
    Call CountSymCombis(Flags, Flags%NumSym22ijkl, iUpperOneR, iLowerW)
    Call CountSymCombis(Flags, Flags%NumSym22ijkl, iUpperOneR, iLowerOneR)
    Call CountSymCombis(Flags, Flags%NumSym22ijkl, iUpperOneR, iLowerOneL)
    Call CountSymCombis(Flags, Flags%NumSym22ilil, iUpperRLdg, iLowerRLdg)
    Call CountSymCombis(Flags, Flags%NumSym22ikjl, iUpperRLdg, iLowerRLcp)
    Call CountSymCombis(Flags, Flags%NumSym22ikjl, iUpperRLco, iLowerRLdg)
    Call CountSymCombis(Flags, Flags%NumSym22ikjl, iUpperRLco, iLowerRLcp)
    Call CountSymCombis(Flags, Flags%NumSym22ikjl, iUpperRLco, iLowerRLcm)
    Call CountSymCombis(Flags, Flags%NumSym22ikjl, iUpperRL,   iLowerRL)
    Call CountSymCombis(Flags, Flags%NumSym22ikjl, iUpperRtRt, iLowerRbR)
    Call CountSymCombis(Flags, Flags%NumSym22ikjl, iUpperRRt,  iLowerRbRb)
    Call CountSymCombis(Flags, Flags%NumSym22ikjl, iUpperRRt,  iLowerRbR)
    Call CountSymCombis(Flags, Flags%NumSym22ikjl, iUpperRRt,  iLowerRRb)
    Call CountSymCombis(Flags, Flags%NumSym22half, iUpperRtRt, iLowerRbRb)
    Call CountSymCombis(Flags, Flags%NumSym31ijkl, iUpperR,    iLowerRb)
    Call CountSymCombis(Flags, Flags%NumSym31ijkl, iUpperL,    iLowerLb)
  End Subroutine CountSymCombinations



  ! Count combinations of upper and lower
  ! partial loops of any symmetry.
  !
  Subroutine CountAllCombinations(Flags)
    Implicit None
    Type(ShavittControl) :: Flags
    ! End of dummy parameters.

    Flags%NumAll13ijkl = 0
    Flags%NumAll22ijkl = 0
    Flags%NumAll22ikjl = 0
    Flags%NumAll22iill = 0
    Flags%NumAll22ilil = 0
    Flags%NumAll22half = 0
    Flags%NumAll31ijkl = 0

    Call CountAllCombis(Flags, Flags%NumAll13ijkl, iUpperRt,   iLowerR)
    Call CountAllCombis(Flags, Flags%NumAll22iill, iUpperW,    iLowerW)
    Call CountAllCombis(Flags, Flags%NumAll22ijkl, iUpperW,    iLowerOneR)
    Call CountAllCombis(Flags, Flags%NumAll22ijkl, iUpperOneR, iLowerW)
    Call CountAllCombis(Flags, Flags%NumAll22ijkl, iUpperOneR, iLowerOneR)
    Call CountAllCombis(Flags, Flags%NumAll22ijkl, iUpperOneR, iLowerOneL)
    Call CountAllCombis(Flags, Flags%NumAll22ilil, iUpperRLdg, iLowerRLdg)
    Call CountAllCombis(Flags, Flags%NumAll22ikjl, iUpperRLdg, iLowerRLcp)
    Call CountAllCombis(Flags, Flags%NumAll22ikjl, iUpperRLco, iLowerRLdg)
    Call CountAllCombis(Flags, Flags%NumAll22ikjl, iUpperRLco, iLowerRLcp)
    Call CountAllCombis(Flags, Flags%NumAll22ikjl, iUpperRLco, iLowerRLcm)
    Call CountAllCombis(Flags, Flags%NumAll22ikjl, iUpperRL,   iLowerRL)
    Call CountAllCombis(Flags, Flags%NumAll22ikjl, iUpperRtRt, iLowerRbR)
    Call CountAllCombis(Flags, Flags%NumAll22ikjl, iUpperRRt,  iLowerRbRb)
    Call CountAllCombis(Flags, Flags%NumAll22ikjl, iUpperRRt,  iLowerRbR)
    Call CountAllCombis(Flags, Flags%NumAll22ikjl, iUpperRRt,  iLowerRRb)
    Call CountAllCombis(Flags, Flags%NumAll22half, iUpperRtRt, iLowerRbRb)
    Call CountAllCombis(Flags, Flags%NumAll31ijkl, iUpperR,    iLowerRb)
    Call CountAllCombis(Flags, Flags%NumAll31ijkl, iUpperL,    iLowerLb)
  End Subroutine CountAllCombinations



  ! Store combinations of upper and lower
  ! partial loops of equal symmetry.
  !
  Subroutine StoreSymCombinations(Flags)
    Implicit None
    Type(ShavittControl) :: Flags
    ! End of dummy parameters.

    Flags%NumOne       = 0
    Flags%NumSym13ijkl = 0
    Flags%NumSym22ijkl = 0
    Flags%NumSym22ikjl = 0
    Flags%NumSym22iill = 0
    Flags%NumSym22ilil = 0
    Flags%NumSym22half = 0
    Flags%NumSym31ijkl = 0

    Call StoreAllCombis(Flags, Flags%CombiListOne,       Flags%NumOne,       iUpperRt,   iLowerRb)
    Call StoreSymCombis(Flags, Flags%CombiListSym13ijkl, Flags%NumSym13ijkl, iUpperRt,   iLowerR)
    Call StoreSymCombis(Flags, Flags%CombiListSym22iill, Flags%NumSym22iill, iUpperW,    iLowerW)
    Call StoreSymCombis(Flags, Flags%CombiListSym22ijkl, Flags%NumSym22ijkl, iUpperW,    iLowerOneR)
    Call StoreSymCombis(Flags, Flags%CombiListSym22ijkl, Flags%NumSym22ijkl, iUpperOneR, iLowerW)
    Call StoreSymCombis(Flags, Flags%CombiListSym22ijkl, Flags%NumSym22ijkl, iUpperOneR, iLowerOneR)
    Call StoreSymCombis(Flags, Flags%CombiListSym22ijkl, Flags%NumSym22ijkl, iUpperOneR, iLowerOneL)
    Call StoreSymCombis(Flags, Flags%CombiListSym22ilil, Flags%NumSym22ilil, iUpperRLdg, iLowerRLdg)
    Call StoreSymCombis(Flags, Flags%CombiListSym22ikjl, Flags%NumSym22ikjl, iUpperRLdg, iLowerRLcp)
    Call StoreSymCombis(Flags, Flags%CombiListSym22ikjl, Flags%NumSym22ikjl, iUpperRLco, iLowerRLdg)
    Call StoreSymCombis(Flags, Flags%CombiListSym22ikjl, Flags%NumSym22ikjl, iUpperRLco, iLowerRLcp)
    Call StoreSymCombis(Flags, Flags%CombiListSym22ikjl, Flags%NumSym22ikjl, iUpperRLco, iLowerRLcm)
    Call StoreSymCombis(Flags, Flags%CombiListSym22ikjl, Flags%NumSym22ikjl, iUpperRL,   iLowerRL)
    Call StoreSymCombis(Flags, Flags%CombiListSym22ikjl, Flags%NumSym22ikjl, iUpperRtRt, iLowerRbR)
    Call StoreSymCombis(Flags, Flags%CombiListSym22ikjl, Flags%NumSym22ikjl, iUpperRRt,  iLowerRbRb)
    Call StoreSymCombis(Flags, Flags%CombiListSym22ikjl, Flags%NumSym22ikjl, iUpperRRt,  iLowerRbR)
    Call StoreSymCombis(Flags, Flags%CombiListSym22ikjl, Flags%NumSym22ikjl, iUpperRRt,  iLowerRRb)
    Call StoreSymCombis(Flags, Flags%CombiListSym22half, Flags%NumSym22half, iUpperRtRt, iLowerRbRb)
    Call StoreSymCombis(Flags, Flags%CombiListSym31ijkl, Flags%NumSym31ijkl, iUpperR,    iLowerRb)
    Call StoreSymCombis(Flags, Flags%CombiListSym31ijkl, Flags%NumSym31ijkl, iUpperL,    iLowerLb)
  End Subroutine StoreSymCombinations



  ! Store combinations of upper and lower
  ! partial loops of any symmetry.
  !
  Subroutine StoreAllCombinations(Flags)
    Implicit None
    Type(ShavittControl) :: Flags
    ! End of dummy parameters.

    Flags%NumAll13ijkl = 0
    Flags%NumAll22ijkl = 0
    Flags%NumAll22ikjl = 0
    Flags%NumAll22iill = 0
    Flags%NumAll22ilil = 0
    Flags%NumAll22half = 0
    Flags%NumAll31ijkl = 0

    Call StoreAllCombis(Flags, Flags%CombiListAll13ijkl, Flags%NumAll13ijkl, iUpperRt,   iLowerR)
    Call StoreAllCombis(Flags, Flags%CombiListAll22iill, Flags%NumAll22iill, iUpperW,    iLowerW)
    Call StoreAllCombis(Flags, Flags%CombiListAll22ijkl, Flags%NumAll22ijkl, iUpperW,    iLowerOneR)
    Call StoreAllCombis(Flags, Flags%CombiListAll22ijkl, Flags%NumAll22ijkl, iUpperOneR, iLowerW)
    Call StoreAllCombis(Flags, Flags%CombiListAll22ijkl, Flags%NumAll22ijkl, iUpperOneR, iLowerOneR)
    Call StoreAllCombis(Flags, Flags%CombiListAll22ijkl, Flags%NumAll22ijkl, iUpperOneR, iLowerOneL)
    Call StoreAllCombis(Flags, Flags%CombiListAll22ilil, Flags%NumAll22ilil, iUpperRLdg, iLowerRLdg)
    Call StoreAllCombis(Flags, Flags%CombiListAll22ikjl, Flags%NumAll22ikjl, iUpperRLdg, iLowerRLcp)
    Call StoreAllCombis(Flags, Flags%CombiListAll22ikjl, Flags%NumAll22ikjl, iUpperRLco, iLowerRLdg)
    Call StoreAllCombis(Flags, Flags%CombiListAll22ikjl, Flags%NumAll22ikjl, iUpperRLco, iLowerRLcp)
    Call StoreAllCombis(Flags, Flags%CombiListAll22ikjl, Flags%NumAll22ikjl, iUpperRLco, iLowerRLcm)
    Call StoreAllCombis(Flags, Flags%CombiListAll22ikjl, Flags%NumAll22ikjl, iUpperRL,   iLowerRL)
    Call StoreAllCombis(Flags, Flags%CombiListAll22ikjl, Flags%NumAll22ikjl, iUpperRtRt, iLowerRbR)
    Call StoreAllCombis(Flags, Flags%CombiListAll22ikjl, Flags%NumAll22ikjl, iUpperRRt,  iLowerRbRb)
    Call StoreAllCombis(Flags, Flags%CombiListAll22ikjl, Flags%NumAll22ikjl, iUpperRRt,  iLowerRbR)
    Call StoreAllCombis(Flags, Flags%CombiListAll22ikjl, Flags%NumAll22ikjl, iUpperRRt,  iLowerRRb)
    Call StoreAllCombis(Flags, Flags%CombiListAll22half, Flags%NumAll22half, iUpperRtRt, iLowerRbRb)
    Call StoreAllCombis(Flags, Flags%CombiListAll31ijkl, Flags%NumAll31ijkl, iUpperR,    iLowerRb)
    Call StoreAllCombis(Flags, Flags%CombiListAll31ijkl, Flags%NumAll31ijkl, iUpperL,    iLowerLb)
  End Subroutine StoreAllCombinations



  ! Count valid combinations of upper partial loops of type UpType
  ! with lower partial loops of type LoType for which the symmetry
  ! of upper and lower partial loop agrees.
  !
  ! A combination is valid if the row indices of the partial bra walks
  ! and the row indices of the partial ket walks, respectively, agree
  ! at the border level.
  !
  ! The partial loops are sorted by increasing bra row index as the
  ! major sorting criterion, then by increasing ket row index, and
  ! finally by their symmetry.
  !
  ! The algorithm works by finding the last partial loop with the
  ! same characteristics as the current partial loop (separately for
  ! upper and lower partial loops). If the characteristics of upper
  ! and lower partial loops match, the combination is counted, and
  ! both pointers to the current upper and lower partial loops are
  ! advanced. If the characteristics of upper and lower partial
  ! loops do not match, only the pointer to the current partial
  ! loop is advanced which comes first in the sorting order.
  !
  Subroutine CountSymCombis(Flags, Num, UpType, LoType)
    Implicit None
    Type(ShavittControl) :: Flags
    Integer              :: Num, UpType, LoType
    ! End of dummy parameters.
    Integer              :: iUp, iUpE, iUpX
    Integer              :: iLo, iLoE, iLoX

    ! Current upper and lower partial loops.
    iUp  = Flags%iFirstUp(UpType)
    iLo  = Flags%iFirstLo(LoType)
    ! Last upper and lower partial loops with the current characteristics.
    iUpE = 0
    iLoE = 0
    ! Last upper and lower partial loops with the requested types.
    iUpX = Flags%iLastUp(UpType)
    iLoX = Flags%iLastLo(LoType)

    Do While (iUp <= iUpX .And. iLo <= iLoX)
       If (iUpE < iUp) Then
          ! Find the last upper partial loop with the same
          ! characteristics as the current upper partial loop.
          iUpE = iUp
          ! Workaround avoiding apparent ifort compiler bug...
          ! Do While (iUpE < iUpX .And. CompareSym(Flags%UpperPart(iUp), Flags%UpperPart(iUpE+1)) == 0)
          Do While (iUpE < iUpX)
             If (CompareSym(Flags%UpperPart(iUp), Flags%UpperPart(iUpE+1)) /= 0) Exit
             iUpE = iUpE + 1
          End Do
       End If
       If (iLoE < iLo) Then
          ! Find the last lower partial loop with the same
          ! characteristics as the current lower partial loop.
          iLoE = iLo
          ! Workaround avoiding apparent ifort compiler bug...
          ! Do While (iLoE < iLoX .And. CompareSym(Flags%LowerPart(iLo), Flags%LowerPart(iLoE+1)) == 0)
          Do While (iLoE < iLoX)
             If (CompareSym(Flags%LowerPart(iLo), Flags%LowerPart(iLoE+1)) /= 0) Exit
             iLoE = iLoE + 1
          End Do
       End If
       If (CompareSym(Flags%UpperPart(iUp), Flags%LowerPart(iLo)) == 0) Then
          ! Characteristics of current upper and lower partial loops match.
          ! Count combination.
          Num = Num  + 1
          ! Advance current upper and lower partial loops.
          iUp = iUpE + 1
          iLo = iLoE + 1
       Else If (CompareSym(Flags%UpperPart(iUp), Flags%LowerPart(iLo)) < 0) Then
          ! Advance current upper partial loop.
          iUp = iUpE + 1
       Else
          ! Advance current lower partial loop.
          iLo = iLoE + 1
       End If
    End Do
  End Subroutine CountSymCombis



  ! All combinations of upper and lower partial loops that have
  ! been counted in the previous subroutine are now stored.
  ! iLast indicates the index of the last combination.
  !
  Subroutine StoreSymCombis(Flags, Combi, iLast, UpType, LoType)
    Implicit None
    Type(ShavittControl)       :: Flags
    Type(Combination), Pointer :: Combi(:)
    Integer                    :: iLast, UpType, LoType
    ! End of dummy parameters.
    Integer                    :: iUp, iUpE, iUpX
    Integer                    :: iLo, iLoE, iLoX

    ! Current upper and lower partial loops.
    iUp  = Flags%iFirstUp(UpType)
    iLo  = Flags%iFirstLo(LoType)
    ! Last upper and lower partial loops with the current characteristics.
    iUpE = 0
    iLoE = 0
    ! Last upper and lower partial loops with the requested types.
    iUpX = Flags%iLastUp(UpType)
    iLoX = Flags%iLastLo(LoType)

    Do While (iUp <= iUpX .And. iLo <= iLoX)
       If (iUpE < iUp) Then
          ! Find the last upper partial loop with the same
          ! characteristics as the current upper partial loop.
          iUpE = iUp
          ! Workaround avoiding apparent ifort compiler bug...
          ! Do While (iUpE < iUpX .And. CompareSym(Flags%UpperPart(iUp), Flags%UpperPart(iUpE+1)) == 0)
          Do While (iUpE < iUpX)
             If (CompareSym(Flags%UpperPart(iUp), Flags%UpperPart(iUpE+1)) /= 0) Exit
             iUpE = iUpE + 1
          End Do
       End If
       If (iLoE < iLo) Then
          ! Find the last lower partial loop with the same
          ! characteristics as the current lower partial loop.
          iLoE = iLo
          ! Workaround avoiding apparent ifort compiler bug...
          ! Do While (iLoE < iLoX .And. CompareSym(Flags%LowerPart(iLo), Flags%LowerPart(iLoE+1)) == 0)
          Do While (iLoE < iLoX)
             If (CompareSym(Flags%LowerPart(iLo), Flags%LowerPart(iLoE+1)) /= 0) Exit
             iLoE = iLoE + 1
          End Do
       End If
       If (CompareSym(Flags%UpperPart(iUp), Flags%LowerPart(iLo)) == 0) Then
          ! Characteristics of current upper and lower partial loops match.
          ! Store combination.
          iLast             = iLast + 1
          Combi(iLast)%iUpB = iUp
          Combi(iLast)%iUpE = iUpE
          Combi(iLast)%iLoB = iLo
          Combi(iLast)%iLoE = iLoE
          ! Advance current upper and lower partial loops.
          iUp = iUpE + 1
          iLo = iLoE + 1
       Else If (CompareSym(Flags%UpperPart(iUp), Flags%LowerPart(iLo)) < 0) Then
          ! Advance current upper partial loop.
          iUp = iUpE + 1
       Else
          ! Advance current lower partial loop.
          iLo = iLoE + 1
       End If
    End Do
  End Subroutine StoreSymCombis



  ! Count valid combinations of upper partial loops of type UpType
  ! with lower partial loops of type LoType regardless of their
  ! symmetry.
  !
  Subroutine CountAllCombis(Flags, Num, UpType, LoType)
    Implicit None
    Type(ShavittControl) :: Flags
    Integer              :: Num, UpType, LoType
    ! End of dummy parameters.
    Integer              :: iUp, iUpE, iUpX
    Integer              :: iLo, iLoE, iLoX

    ! Current upper and lower partial loops.
    iUp  = Flags%iFirstUp(UpType)
    iLo  = Flags%iFirstLo(LoType)
    ! Last upper and lower partial loops with the current characteristics.
    iUpE = 0
    iLoE = 0
    ! Last upper and lower partial loops with the requested types.
    iUpX = Flags%iLastUp(UpType)
    iLoX = Flags%iLastLo(LoType)

    Do While (iUp <= iUpX .And. iLo <= iLoX)
       If (iUpE < iUp) Then
          ! Find the last upper partial loop with the same
          ! characteristics as the current upper partial loop.
          iUpE = iUp
          ! Workaround avoiding apparent ifort compiler bug...
          ! Do While (iUpE < iUpX .And. CompareAll(Flags%UpperPart(iUp), Flags%UpperPart(iUpE+1)) == 0)
          Do While (iUpE < iUpX)
             If (CompareAll(Flags%UpperPart(iUp), Flags%UpperPart(iUpE+1)) /= 0) Exit
             iUpE = iUpE + 1
          End Do
       End If
       If (iLoE < iLo) Then
          ! Find the last lower partial loop with the same
          ! characteristics as the current lower partial loop.
          iLoE = iLo
          ! Workaround avoiding apparent ifort compiler bug...
          ! Do While (iLoE < iLoX .And. CompareAll(Flags%LowerPart(iLo), Flags%LowerPart(iLoE+1)) == 0)
          Do While (iLoE < iLoX)
             If (CompareAll(Flags%LowerPart(iLo), Flags%LowerPart(iLoE+1)) /= 0) Exit
             iLoE = iLoE + 1
          End Do
       End If
       If (CompareAll(Flags%UpperPart(iUp), Flags%LowerPart(iLo)) == 0) Then
          ! Characteristics of current upper and lower partial loops match.
          ! Count combination.
          Num = Num  + 1
          ! Advance current upper and lower partial loops.
          iUp = iUpE + 1
          iLo = iLoE + 1
       Else If (CompareAll(Flags%UpperPart(iUp), Flags%LowerPart(iLo)) < 0) Then
          ! Advance current upper partial loop.
          iUp = iUpE + 1
       Else
          ! Advance current lower partial loop.
          iLo = iLoE + 1
       End If
    End Do
  End Subroutine CountAllCombis



  ! All combinations of upper and lower partial loops that have
  ! been counted in the previous subroutine are now stored.
  ! iLast indicates the index of the last combination.
  !
  Subroutine StoreAllCombis(Flags, Combi, iLast, UpType, LoType)
    Implicit None
    Type(ShavittControl)       :: Flags
    Type(Combination), Pointer :: Combi(:)
    Integer                    :: iLast, UpType, LoType
    ! End of dummy parameters.
    Integer                    :: iUp, iUpE, iUpX
    Integer                    :: iLo, iLoE, iLoX

    ! Current upper and lower partial loops.
    iUp  = Flags%iFirstUp(UpType)
    iLo  = Flags%iFirstLo(LoType)
    ! Last upper and lower partial loops with the current characteristics.
    iUpE = 0
    iLoE = 0
    ! Last upper and lower partial loops with the requested types.
    iUpX = Flags%iLastUp(UpType)
    iLoX = Flags%iLastLo(LoType)

    Do While (iUp <= iUpX .And. iLo <= iLoX)
       If (iUpE < iUp) Then
          ! Find the last upper partial loop with the same
          ! characteristics as the current upper partial loop.
          iUpE = iUp
          ! Workaround avoiding apparent ifort compiler bug...
          ! Do While (iUpE < iUpX .And. CompareAll(Flags%UpperPart(iUp), Flags%UpperPart(iUpE+1)) == 0)
          Do While (iUpE < iUpX)
             If (CompareAll(Flags%UpperPart(iUp), Flags%UpperPart(iUpE+1)) /= 0) Exit
             iUpE = iUpE + 1
          End Do
       End If
       If (iLoE < iLo) Then
          ! Find the last lower partial loop with the same
          ! characteristics as the current lower partial loop.
          iLoE = iLo
          ! Workaround avoiding apparent ifort compiler bug...
          ! Do While (iLoE < iLoX .And. CompareAll(Flags%LowerPart(iLo), Flags%LowerPart(iLoE+1)) == 0)
          Do While (iLoE < iLoX)
             If (CompareAll(Flags%LowerPart(iLo), Flags%LowerPart(iLoE+1)) /= 0) Exit
             iLoE = iLoE + 1
          End Do
       End If
       If (CompareAll(Flags%UpperPart(iUp), Flags%LowerPart(iLo)) == 0) Then
          ! Characteristics of current upper and lower partial loops match.
          ! Store combination.
          iLast             = iLast + 1
          Combi(iLast)%iUpB = iUp
          Combi(iLast)%iUpE = iUpE
          Combi(iLast)%iLoB = iLo
          Combi(iLast)%iLoE = iLoE
          ! Advance current upper and lower partial loops.
          iUp = iUpE + 1
          iLo = iLoE + 1
       Else If (CompareAll(Flags%UpperPart(iUp), Flags%LowerPart(iLo)) < 0) Then
          ! Advance current upper partial loop.
          iUp = iUpE + 1
       Else
          ! Advance current lower partial loop.
          iLo = iLoE + 1
       End If
    End Do
  End Subroutine StoreAllCombis



  ! This function compares two partial loops P and Q considering
  ! ibra, iket, and iProd(sybra, syket), and returns an integer
  !    < 0 if P <  Q,
  !   == 0 if P == Q, and
  !    > 0 if P >  Q.
  !
  Integer Function CompareSym(P, Q)
    Implicit None
    Type(PartialLoop) :: P, Q
    ! End of dummy parameters.
    Integer           :: idiff

    idiff = P%ibra - Q%ibra

    If (idiff /= 0) Then
       CompareSym = idiff
       Return
    End If

    idiff = P%iket - Q%iket

    If (idiff /= 0) Then
       CompareSym = idiff
       Return
    End If

    CompareSym = iProd(P%sybra, P%syket) - iProd(Q%sybra, Q%syket)
  End Function CompareSym



  ! This function compares two partial loops P and Q
  ! considering ibra and iket, and returns an integer
  !    < 0 if P <  Q,
  !   == 0 if P == Q, and
  !    > 0 if P >  Q.
  !
  Integer Function CompareAll(P, Q)
    Implicit None
    Type(PartialLoop) :: P, Q
    ! End of dummy parameters.
    Integer           :: idiff

    idiff = P%ibra - Q%ibra

    If (idiff /= 0) Then
       CompareAll = idiff
       Return
    End If

    CompareAll = P%iket - Q%iket
  End Function CompareAll



  ! Print list of combinations of upper and lower
  ! partial loops for debugging purposes.
  !
  Subroutine PrintCombinations(Combi, Num)
    Implicit None
    Type(Combination), Pointer :: Combi(:)
    Integer                    :: Num
    ! End of dummy parameters.
    Integer                    :: i

    Write(Stdout,'(/1X,"index   iUpB   iUpE   iLoB   iLoE")')
    Do i=1, Num
       Write(Stdout,'(1X,I5,4I7)') i, Combi(i)%iUpB, Combi(i)%iUpE, Combi(i)%iLoB, Combi(i)%iLoE
    End Do
  End Subroutine PrintCombinations

End Module gugashape
