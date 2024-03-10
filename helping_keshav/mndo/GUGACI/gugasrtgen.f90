!
! This module provides subroutines for sorting the generator matrix elements.
! These routines are useful mostly for debugging purposes where a unique order
! of the generator matrix elements is helpful.
!
! Written by Axel Koslowski in 2012-2013 at MPI Muelheim.
!

Module gugasrtgen

Use gugaglobal
Use gugautils, Only: PrintTime, Split1

Implicit None

Private

Public :: ReorgGenerators
Public :: ArrangeIJKL
Public :: SortGenerators
Public :: QuickSortGen
Public :: DebugSortGen

Contains

  ! Allocate and initialize temporary buffers TmpGen and TmpIdx
  ! with the conventions used in the loop driven algorithm.
  ! This routine is needed only for printing the list of
  ! generator matrix elements for debugging purposes.
  !
  Subroutine ReorgGenerators(List, TmpGen, TmpIdx)
    Implicit None
    Type(GenBufList), Pointer :: List(:)
    Type(Generator),  Pointer :: TmpGen(:)
    Integer,          Pointer :: TmpIdx(:)
    ! End of dummy parameters.
    Type(GenBuffer),  Pointer :: bufptr
    Type(Generator)           :: Tmp
    Integer                   :: maxgen, maxidx, ioff
    Integer                   :: irow, ncount, igen, i

    ! Obtain vector dimensions.
    maxidx = UBound(List,1) + 1
    maxgen = 0

    ! Allocate and initialize.
    Allocate(TmpIdx(0:maxidx))
    TmpIdx(:) = 0

    ! Count generator matrix elements in
    ! each row of the CI Hamiltonian.
    Do ioff=0, maxidx-1
       bufptr => List(ioff)%first
       Do While (Associated(bufptr))
          Do i=1, bufptr%count
             irow         = bufptr%gen(i)%icol - ioff
             TmpIdx(irow) = TmpIdx(irow) + 1
             maxgen       = maxgen       + 1
          End Do
          bufptr => bufptr%next
       End Do
    End Do

    ! Make TmpIdx(irow) point to the last element of the previous row
    ! so it may be used to index the current element of the current row
    ! during copying of the generator matrix elements.
    i = 0
    Do irow=1, maxidx
       ncount       = TmpIdx(irow)
       TmpIdx(irow) = i
       i            = i + ncount
    End Do

    ! Allocate temporary generator matrix element buffer.
    Allocate(TmpGen(maxgen))

    ! Copy generator matrix elements to the temporary buffer.
    ! TmpIdx(irow) will finally point to the last element
    ! of the current row.
    Do ioff=0, maxidx-1
       bufptr => List(ioff)%first
       Do While (Associated(bufptr))
          Do i=1, bufptr%count
             irow         = bufptr%gen(i)%icol - ioff
             igen         = TmpIdx(irow) + 1
             TmpGen(igen) = bufptr%gen(i)
             TmpIdx(irow) = igen
          End Do
          bufptr => bufptr%next
       End Do
    End Do
  End Subroutine ReorgGenerators



  ! Uniquely arrange compound indices p (ij) and q (kl) of the
  ! two-electron generator matrix elements before sorting and
  ! printing the latter for debugging purposes.
  !
  ! A disagreement of the order of p and q between the loop-driven
  ! and the shape-driven algorithm may arise from the mode of loop
  ! construction.
  !
  ! For the loop-driven algorithm, i (at the bottom of the loop)
  ! will always be the smallest orbital index and j, k, and l
  ! will be rearranged so that the corresponding two-electron
  ! integral is (ij|kl).
  !
  ! For the shape-driven algorithm, l (at the top of the loop)
  ! will always be the largest orbital index, and i, j, and k
  ! will be rearranged so that the corresponding two-electron
  ! integral is (ij|kl).
  !
  Subroutine ArrangeIJKL(TwoGen)
    Implicit None
    Type(Generator), Pointer :: TwoGen(:)
    ! End of dummy parameters.
    Integer                  :: ij, kl, h, i, j, k, l, m

    Do m=LBound(TwoGen,1), UBound(TwoGen,1)
       ij = TwoGen(m)%p
       kl = TwoGen(m)%q
       Call Split1(ij, i, j)
       Call Split1(kl, k, l)
       If (i > k) Then
          h           = TwoGen(m)%p
          TwoGen(m)%p = TwoGen(m)%q
          TwoGen(m)%q = h
       End If
    End Do
  End Subroutine ArrangeIJKL



  ! Sort generator matrix elements using straight insertion
  ! before printing for debugging purposes. Straight insertion
  ! is efficient if the elements are already in order with
  ! respect to icol (which is the case if computed by the
  ! shape-driven algorithm) and need to be moved only to meet
  ! the sorting criteria concerning the MO indices. On average,
  ! there will be on the order of two contributions to each
  ! matrix element of the CI Hamiltonian.
  !
  Subroutine SortGenerators(TmpGen, TmpIdx)
    Implicit None
    Type(Generator), Pointer :: TmpGen(:)
    Integer,         Pointer :: TmpIdx(:)
    ! End of dummy parameters.
    Type(Generator)          :: Tmp
    Integer                  :: maxidx, irow, ifirst, i, j

    maxidx = UBound(TmpIdx,1)
    Do irow=1, maxidx
       ifirst = TmpIdx(irow-1)+1
       Do i=ifirst+1, TmpIdx(irow)
          Tmp = TmpGen(i)
          j   = i
          ! Workaround avoiding apparent ifort compiler bug...
          ! Do While (j > iFirst .And. CompareGen(Tmp, TmpGen(j-1)) < 0)
          Do While (j > iFirst)
             If (CompareGen(Tmp, TmpGen(j-1)) >= 0) Exit
             TmpGen(j) = TmpGen(j-1)
             j         = j - 1
          End Do
          TmpGen(j) = Tmp
       End Do
    End Do
  End Subroutine SortGenerators



  ! Sort generator matrix elements of each row of the CI Hamiltonian.
  ! This routine has become obsolete since construction of the sparse
  ! CI Hamiltonian for in-core CI using the shape-driven algorithm
  ! no longer needs sorting.
  !
  Subroutine QuickSortGenerators(Flags)
    Implicit None
    Type(ShavittControl) :: Flags
    ! End of dummy parameters.
    Integer              :: i
    Double Precision     :: tuser, tsys, twall

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(/1X,"Sorting generator matrix elements...")')
       Call GetTime(tuser, tsys, twall)
    End If

    Do i=1, Flags%TotalCSF
       Call QuickSortGen(Flags%OneGen, Flags%OneIndex(i-1)+1, Flags%OneIndex(i))
       Call QuickSortGen(Flags%TwoGen, Flags%TwoIndex(i-1)+1, Flags%TwoIndex(i))
       ! Use instead for debugging purposes.
       ! Call DebugSortGen(Flags%OneGen, Flags%OneIndex(i-1)+1, Flags%OneIndex(i))
       ! Call DebugSortGen(Flags%TwoGen, Flags%TwoIndex(i-1)+1, Flags%TwoIndex(i))
    End Do

    ! Print timing of this subroutine if requested.
    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)
  End Subroutine QuickSortGenerators



  ! This subroutine sorts the generator matrix elements in the buffer Gen
  ! between (and including) the indices iFirst and iLast using the Quicksort
  ! algorithm by ascending icol. It may replace subroutine OldQuickSortGen
  ! in module gugaincore (in gugaincore.f90) as it is slightly more efficient.
  !
  Recursive Subroutine QuickSortGen(Gen, iFirst, iLast)
    Implicit None
    Type(Generator), Pointer :: Gen(:)
    Integer                  :: iFirst, iLast
    ! End of dummy parameters.
    Type(Generator)          :: Temp
    Integer                  :: iLeft, iRight, iCenter, i, j

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
    ! the innermost loops and may speed up the algorithm.

    ! Quicksort will be applied if there are at least eight elements,
    ! otherwise a straight insertion sort will be done.
    If (iLast - iFirst > 6) Then
       iCenter = (iFirst + iLast) / 2
       iLeft   = iFirst + 1
       iRight  = iLast

       ! Move the iCenter element to position iFirst+1.
       Temp          = Gen(iFirst+1)
       Gen(iFirst+1) = Gen(iCenter)
       Gen(iCenter)  = Temp

       ! Rearrange so that Gen(iFirst) <= Gen(iFirst+1) <= Gen(iLast).
       If (Gen(iFirst)%icol > Gen(iLast)%icol) Then
          Temp          = Gen(iFirst)
          Gen(iFirst)   = Gen(iLast)
          Gen(iLast)    = Temp
       End If
       If (Gen(iFirst+1)%icol > Gen(iLast)%icol) Then
          Temp          = Gen(iFirst+1)
          Gen(iFirst+1) = Gen(iLast)
          Gen(iLast)    = Temp
       End If
       If (Gen(iFirst)%icol > Gen(iFirst+1)%icol) Then
          Temp          = Gen(iFirst)
          Gen(iFirst)   = Gen(iFirst+1)
          Gen(iFirst+1) = Temp
       End If

       ! Gen(iFirst+1) is the partitioning element.
       Do
          ! Find the leftmost element that is not less than Gen(iFirst+1).
          Do
             iLeft = iLeft + 1
             If (Gen(iLeft)%icol >= Gen(iFirst+1)%icol) Exit
          End Do

          ! Find the rightmost element that is not greater than Gen(iFirst+1).
          Do
             iRight = iRight - 1
             If (Gen(iRight)%icol <= Gen(iFirst+1)%icol) Exit
          End Do

          ! Stop if pointers crossed.
          If (iLeft > iRight) Exit

          ! We have found two incorrectly positioned elements.
          ! Switch them to get them into the correct order.
          Temp        = Gen(iLeft)
          Gen(iLeft)  = Gen(iRight)
          Gen(iRight) = Temp
       End Do

       ! At this point, the situation may look like this:
       ! iFirst   iFirst+1   iFirst+2   iFirst+3   iRight   iLeft   iLast-1   iLast
       !    3         5          1          2         4       6        7        8
       !
       ! Or, iLeft and iRight previously stopped at the same element (at iLast-3):
       ! iFirst   iFirst+1   iFirst+2    iRight   iLast-3   iLeft   iLast-1   iLast
       !    3         5          1          2         5       6        7        8

       Temp          = Gen(iFirst+1)
       Gen(iFirst+1) = Gen(iRight)
       Gen(iRight)   = Temp

       Call QuickSortGen(Gen, iFirst, iRight-1)
       Call QuickSortGen(Gen, iLeft,  iLast)
    Else
       ! Straight insertion sort.
       ! At every iteration of the outermost loop, insert
       ! the current element (at position i) correctly into
       ! the previously sorted sequence at iFirst ... i-1.
       Do i=iFirst+1, iLast
          Temp = Gen(i)
          j    = i
          ! Workaround avoiding possible ifort compiler bug...
          ! Do While (j > iFirst .And. Temp%icol < Gen(j-1)%icol)
          Do While (j > iFirst)
             If (Temp%icol >= Gen(j-1)%icol) Exit
             Gen(j) = Gen(j-1)
             j      = j - 1
          End Do
          Gen(j) = Temp
       End Do
    End If
  End Subroutine QuickSortGen



  ! This subroutine sorts the generator matrix elements in the buffer Gen
  ! between (and including) the indices iFirst and iLast using the Quicksort
  ! algorithm by ascending icol, ascending p, and ascending q, in this order.
  ! It may replace subroutine OldQuickSortGen in module gugaincore (in
  ! gugaincore.f90) if a unique order of the generator matrix elements
  ! is needed for debugging purposes.
  !
  Recursive Subroutine DebugSortGen(Gen, iFirst, iLast)
    Implicit None
    Type(Generator), Pointer :: Gen(:)
    Integer                  :: iFirst, iLast
    ! End of dummy parameters.
    Type(Generator)          :: Temp
    Integer                  :: iLeft, iRight, iCenter, i, j

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
    ! to the complicated comparison of two generator matrix elements in
    ! this routine which is intended for debugging purposes, the benefit
    ! of this technique may be small. Therefore, we also use recursive
    ! subroutine calls for the benefit of simpler code.

    ! Quicksort will be applied if there are at least eight elements,
    ! otherwise a straight insertion sort will be done.
    If (iLast - iFirst > 6) Then
       iCenter = (iFirst + iLast) / 2
       iLeft   = iFirst + 1
       iRight  = iLast

       ! Move the iCenter element to position iFirst+1.
       Temp          = Gen(iFirst+1)
       Gen(iFirst+1) = Gen(iCenter)
       Gen(iCenter)  = Temp

       ! Rearrange so that Gen(iFirst) <= Gen(iFirst+1) <= Gen(iLast).
       If (CompareGen(Gen(iFirst), Gen(iLast)) > 0) Then
          Temp          = Gen(iFirst)
          Gen(iFirst)   = Gen(iLast)
          Gen(iLast)    = Temp
       End If
       If (CompareGen(Gen(iFirst+1), Gen(iLast)) > 0) Then
          Temp          = Gen(iFirst+1)
          Gen(iFirst+1) = Gen(iLast)
          Gen(iLast)    = Temp
       End If
       If (CompareGen(Gen(iFirst), Gen(iFirst+1)) > 0) Then
          Temp          = Gen(iFirst)
          Gen(iFirst)   = Gen(iFirst+1)
          Gen(iFirst+1) = Temp
       End If

       ! Gen(iFirst+1) is the partitioning element.
       Do
          ! Find the leftmost element that is not less than Gen(iFirst+1).
          Do
             iLeft = iLeft + 1
             If (CompareGen(Gen(iLeft), Gen(iFirst+1)) >= 0) Exit
          End Do

          ! Find the rightmost element that is not greater than Gen(iFirst+1).
          Do
             iRight = iRight - 1
             If (CompareGen(Gen(iRight), Gen(iFirst+1)) <= 0) Exit
          End Do

          ! Stop if pointers crossed.
          If (iLeft > iRight) Exit

          ! We have found two incorrectly positioned elements.
          ! Switch them to get them into the correct order.
          Temp        = Gen(iLeft)
          Gen(iLeft)  = Gen(iRight)
          Gen(iRight) = Temp
       End Do

       ! At this point, the situation may look like this:
       ! iFirst   iFirst+1   iFirst+2   iFirst+3   iRight   iLeft   iLast-1   iLast
       !    3         5          1          2         4       6        7        8
       !
       ! Or, iLeft and iRight previously stopped at the same element (at iLast-3):
       ! iFirst   iFirst+1   iFirst+2    iRight   iLast-3   iLeft   iLast-1   iLast
       !    3         5          1          2         5       6        7        8

       Temp          = Gen(iFirst+1)
       Gen(iFirst+1) = Gen(iRight)
       Gen(iRight)   = Temp

       Call DebugSortGen(Gen, iFirst, iRight-1)
       Call DebugSortGen(Gen, iLeft,  iLast)
    Else
       ! Straight insertion sort.
       ! At every iteration of the outermost loop, insert
       ! the current element (at position i) correctly into
       ! the previously sorted sequence at iFirst ... i-1.
       Do i=iFirst+1, iLast
          Temp = Gen(i)
          j    = i
          ! Workaround avoiding ifort compiler bug...
          ! Do While (j > iFirst .And. CompareGen(Temp, Gen(j-1)) < 0)
          Do While (j > iFirst)
             If (CompareGen(Temp, Gen(j-1)) >= 0) Exit
             Gen(j) = Gen(j-1)
             j      = j - 1
          End Do
          Gen(j) = Temp
       End Do
    End If
  End Subroutine DebugSortGen



  ! This function compares two generator matrix elements considering
  ! icol, p, and q, and returns an integer
  !    < 0 if gen1 <  gen2,
  !   == 0 if gen1 == gen2, and
  !    > 0 if gen1 >  gen2.
  !
  Integer Function CompareGen(gen1, gen2)
    Type(Generator) :: gen1, gen2
    ! End of dummy parameters.
    Integer         :: idiff

    idiff = gen1%icol - gen2%icol

    If (idiff /= 0) Then
       CompareGen = idiff
       Return
    End If

    idiff = gen1%p - gen2%p

    If (idiff /= 0) Then
       CompareGen = idiff
       Return
    End If

    CompareGen = gen1%q - gen2%q
  End Function CompareGen

End Module gugasrtgen
