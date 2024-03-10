!
! This module provides worker routines for the shape-driven GUGA-CI approach
! calculating and storing generator matrix elements into linked lists of
! buffers for each diagonal of the CI Hamiltonian. States of any symmetry
! are considered. For the two-electron generator matrix elements, the symmetry
! of both states involved must agree.
!
! Written by Axel Koslowski in 2012-2013 at MPI Muelheim.
!

Module gugashape1

Use gugaglobal
Use gugaloops, Only: Tiny

Implicit None

Private

Public :: StoreGenerators

Contains

  ! Assemble upper and lower partial loops and form generator matrix elements,
  ! also from complete loops of the upper and lower part of the DRT. Store the
  ! generator matrix elements into linked lists of buffers (Flags%OneList and
  ! Flags%TwoList) for each diagonal of the CI Hamiltonian. The generator
  ! matrix elements will be the same as formed by calling CalcGenerators
  ! (in module gugaloops) with Irrep=0.
  !
  Subroutine StoreGenerators(Flags)
    Implicit None
    Type(ShavittControl) :: Flags
    ! End of dummy parameters.
    Integer              :: isym

    ! Process combinations of upper and lower partial loops.
    ! A subroutine like CombineSym22iljk does not exist because
    ! indices i and j of the lower partial loop have been switched
    ! during partial loop construction such that the corresponding
    ! two-electron integral is (ik|jl) and this case is covered by
    ! CombineSym22ikjl.
    !
    Call CombineOne(Flags)
    Call CombineSym13ijkl(Flags)
    Call CombineSym22ijkl(Flags)
    Call CombineSym22iill(Flags)
    Call CombineSym22ikjl(Flags)
    Call CombineSym22ilil(Flags)
    Call CombineSym22half(Flags)
    Call CombineSym31ijkl(Flags)

    ! Process complete one-body loops of any symmetry in the upper and lower part of the DRT.
    Do isym=1, Flags%NumSym
       Call Complete(Flags, Flags%OneList, Flags%OneLoop, Flags%iFirst1(isym), Flags%iLast1(isym))
    End Do

    ! Process complete two-body loops in the upper and lower part of the DRT.
    ! Since the symmetry of both states involved in the matrix elements must
    ! agree, only totally symmetric loops need to be considered.
    Call Complete(Flags, Flags%TwoList, Flags%TwoLoop, Flags%iFirst2(1), Flags%iLast2(1))
  End Subroutine StoreGenerators



  ! Process all combinations of upper partial loops (with index l)
  ! and lower partial loops (with index i) forming an off-diagonal
  ! one-electron loop. Matrix elements between states of ANY SYMMETRY
  ! are considered. Since the states in the index vector are ordered
  ! by symmetry as the major sorting criterion, states of DIFFERENT
  ! SYMMETRY will in general NOT be SORTED according to increasing
  ! LEXICAL INDEX, so the order of the state indices mbra and mket
  ! inferred from the lexical indices lexbra and lexket (with
  ! lexbra < lexket) must be considered explicitly.
  !
  Subroutine CombineOne(Flags)
    Implicit None
    Type(ShavittControl)     :: Flags
    ! End of dummy parameters.
    Integer                  :: iCombi, iUp, iLo
    Integer                  :: lexbra, lexket
    Integer                  :: ybra, yket, mbra, mket
    Integer                  :: igen, ioff, i, j
    Type(GenBuffer), Pointer :: bufptr
    Double Precision         :: value

    ! Loop over all combinations.
    Do iCombi=1, Flags%NumOne
       ! Loop over all upper partial loops.
       Do iUp=Flags%CombiListOne(iCombi)%iUpB, Flags%CombiListOne(iCombi)%iUpE
          ! Loop over all lower partial loops.
          Do iLo=Flags%CombiListOne(iCombi)%iLoB, Flags%CombiListOne(iCombi)%iLoE
             value = Flags%UpperPart(iUp)%singlet * Flags%LowerPart(iLo)%singlet
             If (Dabs(value) > Tiny) Then
                lexbra = Flags%UpperPart(iUp)%lexbra + Flags%LowerPart(iLo)%lexbra
                lexket = Flags%UpperPart(iUp)%lexket + Flags%LowerPart(iLo)%lexket
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   yket = lexket + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      mket = Flags%IndVec(yket + j)
                      ! Store generator matrix element if both CSFs have been selected.
                      If (mbra /= 0 .And. mket /= 0) Then
                         If (mbra <= mket) Then
                            ioff = mket - mbra
                            ! Allocate new buffer if the current one is full.
                            If (.Not. Associated(Flags%OneList(ioff)%first) .Or. &
                                Flags%OneList(ioff)%first%count == genBufSize) Then
                               Allocate(bufptr)
                               bufptr%count              =  0
                               bufptr%next               => Flags%OneList(ioff)%first
                               Flags%OneList(ioff)%first => bufptr
                               Flags%OneList(ioff)%count =  Flags%OneList(ioff)%count + 1
                            End If
                            bufptr                 => Flags%OneList(ioff)%first
                            igen                   =  bufptr%count + 1
                            bufptr%gen(igen)%value =  value
                            bufptr%gen(igen)%icol  =  mket
                            bufptr%gen(igen)%p     =  Flags%LowerPart(iLo)%i
                            bufptr%gen(igen)%q     =  Flags%UpperPart(iUp)%l
                            bufptr%count           =  igen
                         Else
                            ioff = mbra - mket
                            ! Allocate new buffer if the current one is full.
                            If (.Not. Associated(Flags%OneList(ioff)%first) .Or. &
                                Flags%OneList(ioff)%first%count == genBufSize) Then
                               Allocate(bufptr)
                               bufptr%count              =  0
                               bufptr%next               => Flags%OneList(ioff)%first
                               Flags%OneList(ioff)%first => bufptr
                               Flags%OneList(ioff)%count =  Flags%OneList(ioff)%count + 1
                            End If
                            bufptr                 => Flags%OneList(ioff)%first
                            igen                   =  bufptr%count + 1
                            bufptr%gen(igen)%value =  value
                            bufptr%gen(igen)%icol  =  mbra
                            bufptr%gen(igen)%p     =  Flags%UpperPart(iUp)%l
                            bufptr%gen(igen)%q     =  Flags%LowerPart(iLo)%i
                            bufptr%count           =  igen
                         End If
                      End If
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineOne



  ! Process all combinations of upper partial loops with one index (l)
  ! and lower partial loops with three indices (i, j, and k) forming an
  ! off-diagonal two-electron loop with partial bra and ket walk of the
  ! same symmetry. In this case state index mbra (inferred from lexbra)
  ! will be less than mket (inferred from lexket) because states of equal
  ! symmetry are sorted according to increasing lexical index. During
  ! lower partial loop construction, indices i, j, and k have been
  ! rearranged such that the corresponding two-electron integral is
  ! always (ij|kl).
  !
  Subroutine CombineSym13ijkl(Flags)
    Implicit None
    Type(ShavittControl)     :: Flags
    ! End of dummy parameters.
    Integer                  :: iCombi, iUp, iLo
    Integer                  :: lexbra, lexket
    Integer                  :: ybra, yket, mbra, mket
    Integer                  :: igen, ioff, ij, kl, i, j, k, l
    Type(GenBuffer), Pointer :: bufptr
    Double Precision         :: value

    ! Loop over all combinations.
    Do iCombi=1, Flags%NumSym13ijkl
       ! Loop over all upper partial loops.
       Do iUp=Flags%CombiListSym13ijkl(iCombi)%iUpB, Flags%CombiListSym13ijkl(iCombi)%iUpE
          ! Loop over all lower partial loops.
          Do iLo=Flags%CombiListSym13ijkl(iCombi)%iLoB, Flags%CombiListSym13ijkl(iCombi)%iLoE
             value = Flags%UpperPart(iUp)%singlet * Flags%LowerPart(iLo)%singlet
             If (Dabs(value) > Tiny) Then
                lexbra = Flags%UpperPart(iUp)%lexbra + Flags%LowerPart(iLo)%lexbra
                lexket = Flags%UpperPart(iUp)%lexket + Flags%LowerPart(iLo)%lexket
                k      = Flags%LowerPart(iLo)%k
                l      = Flags%UpperPart(iUp)%l
                ij     = Flags%LowerPart(iLo)%l
                ! l as the only index from the upper partial loop is always largest here.
                kl     = (l*(l-1))/2 + k
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   yket = lexket + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      mket = Flags%IndVec(yket + j)
                      ! Store generator matrix element if both CSFs have been selected.
                      If (mbra /= 0 .And. mket /= 0) Then
                         ioff = mket - mbra
                         ! Allocate new buffer if the current one is full.
                         If (.Not. Associated(Flags%TwoList(ioff)%first) .Or. &
                             Flags%TwoList(ioff)%first%count == genBufSize) Then
                            Allocate(bufptr)
                            bufptr%count              =  0
                            bufptr%next               => Flags%TwoList(ioff)%first
                            Flags%TwoList(ioff)%first => bufptr
                            Flags%TwoList(ioff)%count =  Flags%TwoList(ioff)%count + 1
                         End If
                         bufptr                 => Flags%TwoList(ioff)%first
                         igen                   =  bufptr%count + 1
                         bufptr%gen(igen)%value =  value
                         bufptr%gen(igen)%icol  =  mket
                         bufptr%gen(igen)%p     =  ij
                         bufptr%gen(igen)%q     =  kl
                         bufptr%count           =  igen
                      End If
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineSym13ijkl



  ! Process all combinations of upper partial loops with indices
  ! k and l and lower partial loops with indices i and j forming
  ! an off-diagonal two-electron loop with partial bra and ket
  ! walk of the same symmetry. In this case state index mbra
  ! (inferred from lexbra) will be less than mket (inferred
  ! from lexket) because states of equal symmetry are sorted
  ! according to increasing lexical index. The corresponding
  ! two-electron integral is (ij|kl).
  !
  Subroutine CombineSym22ijkl(Flags)
    Implicit None
    Type(ShavittControl)     :: Flags
    ! End of dummy parameters.
    Integer                  :: iCombi, iUp, iLo
    Integer                  :: lexbra, lexket
    Integer                  :: ybra, yket, mbra, mket
    Integer                  :: igen, ioff, i, j
    Type(GenBuffer), Pointer :: bufptr
    Double Precision         :: value

    ! Loop over all combinations.
    Do iCombi=1, Flags%NumSym22ijkl
       ! Loop over all upper partial loops.
       Do iUp=Flags%CombiListSym22ijkl(iCombi)%iUpB, Flags%CombiListSym22ijkl(iCombi)%iUpE
          ! Loop over all lower partial loops.
          Do iLo=Flags%CombiListSym22ijkl(iCombi)%iLoB, Flags%CombiListSym22ijkl(iCombi)%iLoE
             value = Flags%UpperPart(iUp)%singlet * Flags%LowerPart(iLo)%singlet
             If (Dabs(value) > Tiny) Then
                lexbra = Flags%UpperPart(iUp)%lexbra + Flags%LowerPart(iLo)%lexbra
                lexket = Flags%UpperPart(iUp)%lexket + Flags%LowerPart(iLo)%lexket
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   yket = lexket + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      mket = Flags%IndVec(yket + j)
                      ! Store generator matrix element if both CSFs have been selected.
                      If (mbra /= 0 .And. mket /= 0) Then
                         ioff = mket - mbra
                         ! Allocate new buffer if the current one is full.
                         If (.Not. Associated(Flags%TwoList(ioff)%first) .Or. &
                             Flags%TwoList(ioff)%first%count == genBufSize) Then
                            Allocate(bufptr)
                            bufptr%count              =  0
                            bufptr%next               => Flags%TwoList(ioff)%first
                            Flags%TwoList(ioff)%first => bufptr
                            Flags%TwoList(ioff)%count =  Flags%TwoList(ioff)%count + 1
                         End If
                         bufptr                 => Flags%TwoList(ioff)%first
                         igen                   =  bufptr%count + 1
                         bufptr%gen(igen)%value =  value
                         bufptr%gen(igen)%icol  =  mket
                         bufptr%gen(igen)%p     =  Flags%LowerPart(iLo)%l  ! ij
                         bufptr%gen(igen)%q     =  Flags%UpperPart(iUp)%i  ! kl
                         bufptr%count           =  igen
                      End If
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineSym22ijkl



  ! Process all combinations of upper partial loops with identical
  ! indices k and l and lower partial loops with identical indices
  ! i and j forming a diagonal two-electron loop. The corresponding
  ! two-electron integral is (ii|ll).
  !
  Subroutine CombineSym22iill(Flags)
    Implicit None
    Type(ShavittControl)     :: Flags
    ! End of dummy parameters.
    Integer                  :: iCombi, iUp, iLo
    Integer                  :: lexbra, ybra, mbra
    Integer                  :: igen, i, j
    Integer, Parameter       :: ioff = 0
    Type(GenBuffer), Pointer :: bufptr
    Double Precision         :: value

    ! Loop over all combinations.
    Do iCombi=1, Flags%NumSym22iill
       ! Loop over all upper partial loops.
       Do iUp=Flags%CombiListSym22iill(iCombi)%iUpB, Flags%CombiListSym22iill(iCombi)%iUpE
          ! Loop over all lower partial loops.
          Do iLo=Flags%CombiListSym22iill(iCombi)%iLoB, Flags%CombiListSym22iill(iCombi)%iLoE
             value = Flags%UpperPart(iUp)%singlet * Flags%LowerPart(iLo)%singlet
             If (Dabs(value) > Tiny) Then
                lexbra = Flags%UpperPart(iUp)%lexbra + Flags%LowerPart(iLo)%lexbra
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      ! Store generator matrix element if CSF has been selected.
                      If (mbra /= 0) Then
                         ! Allocate new buffer if the current one is full.
                         If (.Not. Associated(Flags%TwoList(ioff)%first) .Or. &
                             Flags%TwoList(ioff)%first%count == genBufSize) Then
                            Allocate(bufptr)
                            bufptr%count              =  0
                            bufptr%next               => Flags%TwoList(ioff)%first
                            Flags%TwoList(ioff)%first => bufptr
                            Flags%TwoList(ioff)%count =  Flags%TwoList(ioff)%count + 1
                         End If
                         bufptr                 => Flags%TwoList(ioff)%first
                         igen                   =  bufptr%count + 1
                         bufptr%gen(igen)%value =  value
                         bufptr%gen(igen)%icol  =  mbra
                         bufptr%gen(igen)%p     =  Flags%LowerPart(iLo)%l  ! ii
                         bufptr%gen(igen)%q     =  Flags%UpperPart(iUp)%i  ! ll
                         bufptr%count           =  igen
                      End If
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineSym22iill



  ! Process all combinations of upper partial loops with indices
  ! k and l and lower partial loops with indices i and j forming
  ! an off-diagonal two-electron loop with partial bra and ket
  ! walk of the same symmetry. In this case state index mbra
  ! (inferred from lexbra) will be less than mket (inferred
  ! from lexket) because states of equal symmetry are sorted
  ! according to increasing lexical index. The corresponding
  ! two-electron integral is (ik|jl). By switching indices i and j
  ! during lower partial loop construction this routine covers
  ! also those cases in which the integral would be (il|jk).
  !
  Subroutine CombineSym22ikjl(Flags)
    Implicit None
    Type(ShavittControl)     :: Flags
    ! End of dummy parameters.
    Integer                  :: iCombi, iUp, iLo
    Integer                  :: lexbra, lexket
    Integer                  :: ybra, yket, mbra, mket
    Integer                  :: igen, ioff, ik, jl, i, j, k, l
    Type(GenBuffer), Pointer :: bufptr
    Double Precision         :: value

    ! Loop over all combinations.
    Do iCombi=1, Flags%NumSym22ikjl
       ! Loop over all upper partial loops.
       Do iUp=Flags%CombiListSym22ikjl(iCombi)%iUpB, Flags%CombiListSym22ikjl(iCombi)%iUpE
          ! Loop over all lower partial loops.
          Do iLo=Flags%CombiListSym22ikjl(iCombi)%iLoB, Flags%CombiListSym22ikjl(iCombi)%iLoE
             value = Flags%UpperPart(iUp)%singlet * Flags%LowerPart(iLo)%singlet + &
                     Flags%UpperPart(iUp)%triplet * Flags%LowerPart(iLo)%triplet
             If (Dabs(value) > Tiny) Then
                lexbra = Flags%UpperPart(iUp)%lexbra + Flags%LowerPart(iLo)%lexbra
                lexket = Flags%UpperPart(iUp)%lexket + Flags%LowerPart(iLo)%lexket
                i      = Flags%LowerPart(iLo)%i
                j      = Flags%LowerPart(iLo)%j
                k      = Flags%UpperPart(iUp)%k
                l      = Flags%UpperPart(iUp)%l
                ! k and l from the upper partial loop are both greater
                ! than i and j from the lower partial loop here.
                ik     = (k*(k-1))/2 + i
                jl     = (l*(l-1))/2 + j
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   yket = lexket + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      mket = Flags%IndVec(yket + j)
                      ! Store generator matrix element if both CSFs have been selected.
                      If (mbra /= 0 .And. mket /= 0) Then
                         ioff = mket - mbra
                         ! Allocate new buffer if the current one is full.
                         If (.Not. Associated(Flags%TwoList(ioff)%first) .Or. &
                             Flags%TwoList(ioff)%first%count == genBufSize) Then
                            Allocate(bufptr)
                            bufptr%count              =  0
                            bufptr%next               => Flags%TwoList(ioff)%first
                            Flags%TwoList(ioff)%first => bufptr
                            Flags%TwoList(ioff)%count =  Flags%TwoList(ioff)%count + 1
                         End If
                         bufptr                 => Flags%TwoList(ioff)%first
                         igen                   =  bufptr%count + 1
                         bufptr%gen(igen)%value =  value
                         bufptr%gen(igen)%icol  =  mket
                         bufptr%gen(igen)%p     =  ik
                         bufptr%gen(igen)%q     =  jl
                         bufptr%count           =  igen
                      End If
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineSym22ikjl



  ! Process all combinations of upper partial loops with identical
  ! indices k and l and lower partial loops with identical indices
  ! i and j forming a diagonal two-electron loop. The corresponding
  ! two-electron integral is (il|il).
  !
  Subroutine CombineSym22ilil(Flags)
    Implicit None
    Type(ShavittControl)     :: Flags
    ! End of dummy parameters.
    Integer                  :: iCombi, iUp, iLo
    Integer                  :: lexbra, ybra, mbra
    Integer                  :: igen, il, i, j, l
    Integer, Parameter       :: ioff = 0
    Type(GenBuffer), Pointer :: bufptr
    Double Precision         :: value

    ! Loop over all combinations.
    Do iCombi=1, Flags%NumSym22ilil
       ! Loop over all upper partial loops.
       Do iUp=Flags%CombiListSym22ilil(iCombi)%iUpB, Flags%CombiListSym22ilil(iCombi)%iUpE
          ! Loop over all lower partial loops.
          Do iLo=Flags%CombiListSym22ilil(iCombi)%iLoB, Flags%CombiListSym22ilil(iCombi)%iLoE
             value = Flags%UpperPart(iUp)%singlet * Flags%LowerPart(iLo)%singlet + &
                     Flags%UpperPart(iUp)%triplet * Flags%LowerPart(iLo)%triplet
             If (Dabs(value) > Tiny) Then
                lexbra = Flags%UpperPart(iUp)%lexbra + Flags%LowerPart(iLo)%lexbra
                i      = Flags%LowerPart(iLo)%i
                l      = Flags%UpperPart(iUp)%l
                ! l from the upper partial loop is always greater
                ! than i from the lower partial loop here.
                il     = (l*(l-1))/2 + i
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      ! Store generator matrix element if CSF has been selected.
                      If (mbra /= 0) Then
                         ! Allocate new buffer if the current one is full.
                         If (.Not. Associated(Flags%TwoList(ioff)%first) .Or. &
                             Flags%TwoList(ioff)%first%count == genBufSize) Then
                            Allocate(bufptr)
                            bufptr%count              =  0
                            bufptr%next               => Flags%TwoList(ioff)%first
                            Flags%TwoList(ioff)%first => bufptr
                            Flags%TwoList(ioff)%count =  Flags%TwoList(ioff)%count + 1
                         End If
                         bufptr                 => Flags%TwoList(ioff)%first
                         igen                   =  bufptr%count + 1
                         bufptr%gen(igen)%value =  value
                         bufptr%gen(igen)%icol  =  mbra
                         bufptr%gen(igen)%p     =  il
                         bufptr%gen(igen)%q     =  il
                         bufptr%count           =  igen
                      End If
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineSym22ilil



  ! Process all combinations of upper partial loops with identical
  ! indices k and l and lower partial loops with identical indices
  ! i and j forming RtRt-RbRb loops for which the corresponding
  ! two-electron integral (il|il) must be divided by two.
  !
  Subroutine CombineSym22half(Flags)
    Implicit None
    Type(ShavittControl)     :: Flags
    ! End of dummy parameters.
    Integer                  :: iCombi, iUp, iLo
    Integer                  :: lexbra, lexket
    Integer                  :: ybra, yket, mbra, mket
    Integer                  :: igen, ioff, il, i, j, l
    Type(GenBuffer), Pointer :: bufptr
    Double Precision         :: value

    ! Loop over all combinations.
    Do iCombi=1, Flags%NumSym22half
       ! Loop over all upper partial loops.
       Do iUp=Flags%CombiListSym22half(iCombi)%iUpB, Flags%CombiListSym22half(iCombi)%iUpE
          ! Loop over all lower partial loops.
          Do iLo=Flags%CombiListSym22half(iCombi)%iLoB, Flags%CombiListSym22half(iCombi)%iLoE
             value = 0.5D0 * (Flags%UpperPart(iUp)%singlet * Flags%LowerPart(iLo)%singlet + &
                              Flags%UpperPart(iUp)%triplet * Flags%LowerPart(iLo)%triplet)
             If (Dabs(value) > Tiny) Then
                lexbra = Flags%UpperPart(iUp)%lexbra + Flags%LowerPart(iLo)%lexbra
                lexket = Flags%UpperPart(iUp)%lexket + Flags%LowerPart(iLo)%lexket
                i      = Flags%LowerPart(iLo)%i
                l      = Flags%UpperPart(iUp)%l
                ! l from the upper partial loop is always greater
                ! than i from the lower partial loop here.
                il     = (l*(l-1))/2 + i
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   yket = lexket + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      mket = Flags%IndVec(yket + j)
                      ! Store generator matrix element if both CSFs have been selected.
                      If (mbra /= 0 .And. mket /= 0) Then
                         ioff = mket - mbra
                         ! Allocate new buffer if the current one is full.
                         If (.Not. Associated(Flags%TwoList(ioff)%first) .Or. &
                             Flags%TwoList(ioff)%first%count == genBufSize) Then
                            Allocate(bufptr)
                            bufptr%count              =  0
                            bufptr%next               => Flags%TwoList(ioff)%first
                            Flags%TwoList(ioff)%first => bufptr
                            Flags%TwoList(ioff)%count =  Flags%TwoList(ioff)%count + 1
                         End If
                         bufptr                 => Flags%TwoList(ioff)%first
                         igen                   =  bufptr%count + 1
                         bufptr%gen(igen)%value =  value
                         bufptr%gen(igen)%icol  =  mket
                         bufptr%gen(igen)%p     =  il
                         bufptr%gen(igen)%q     =  il
                         bufptr%count           =  igen
                      End If
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineSym22half



  ! Process all combinations of upper partial loops with three indices
  ! (j, k, and l) and lower partial loops with one index (i) forming
  ! an off-diagonal two-electron loop with partial bra and ket walk of
  ! the same symmetry. In this case state index mbra (inferred from
  ! lexbra) will be less than mket (inferred from lexket) because states
  ! of equal symmetry are sorted according to increasing lexical index.
  ! During upper partial loop construction, indices j, k, and l have
  ! been rearranged such that the corresponding two-electron integral
  ! is always (ij|kl).
  !
  Subroutine CombineSym31ijkl(Flags)
    Implicit None
    Type(ShavittControl)     :: Flags
    ! End of dummy parameters.
    Integer                  :: iCombi, iUp, iLo
    Integer                  :: lexbra, lexket
    Integer                  :: ybra, yket, mbra, mket
    Integer                  :: igen, ioff, ij, kl, i, j
    Type(GenBuffer), Pointer :: bufptr
    Double Precision         :: value

    ! Loop over all combinations.
    Do iCombi=1, Flags%NumSym31ijkl
       ! Loop over all upper partial loops.
       Do iUp=Flags%CombiListSym31ijkl(iCombi)%iUpB, Flags%CombiListSym31ijkl(iCombi)%iUpE
          ! Loop over all lower partial loops.
          Do iLo=Flags%CombiListSym31ijkl(iCombi)%iLoB, Flags%CombiListSym31ijkl(iCombi)%iLoE
             value = Flags%UpperPart(iUp)%singlet * Flags%LowerPart(iLo)%singlet
             If (Dabs(value) > Tiny) Then
                lexbra = Flags%UpperPart(iUp)%lexbra + Flags%LowerPart(iLo)%lexbra
                lexket = Flags%UpperPart(iUp)%lexket + Flags%LowerPart(iLo)%lexket
                i      = Flags%LowerPart(iLo)%i
                j      = Flags%UpperPart(iUp)%j
                kl     = Flags%UpperPart(iUp)%i
                ! i as the only index from the lower partial loop is always smallest here.
                ij     = (j*(j-1))/2 + i
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   yket = lexket + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      mket = Flags%IndVec(yket + j)
                      ! Store generator matrix element if both CSFs have been selected.
                      If (mbra /= 0 .And. mket /= 0) Then
                         ioff = mket - mbra
                         ! Allocate new buffer if the current one is full.
                         If (.Not. Associated(Flags%TwoList(ioff)%first) .Or. &
                             Flags%TwoList(ioff)%first%count == genBufSize) Then
                            Allocate(bufptr)
                            bufptr%count              =  0
                            bufptr%next               => Flags%TwoList(ioff)%first
                            Flags%TwoList(ioff)%first => bufptr
                            Flags%TwoList(ioff)%count =  Flags%TwoList(ioff)%count + 1
                         End If
                         bufptr                 => Flags%TwoList(ioff)%first
                         igen                   =  bufptr%count + 1
                         bufptr%gen(igen)%value =  value
                         bufptr%gen(igen)%icol  =  mket
                         bufptr%gen(igen)%p     =  ij
                         bufptr%gen(igen)%q     =  kl
                         bufptr%count           =  igen
                      End If
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineSym31ijkl



  ! Process complete one- or two-electron loops with partial bra and ket walks
  ! of any symmetry. Since the states in the index vector are ordered by symmetry
  ! as the major sorting criterion, states of different symmetry will in general
  ! not be sorted according to increasing lexical index, so the order of the state
  ! indices mbra and mket inferred from the lexical indices lexbra and lexket
  ! (with lexbra < lexket) must be considered explicitly.
  !
  Subroutine Complete(Flags, List, Loop, iFirst, iLast)
    Implicit None
    Type(ShavittControl)         :: Flags
    Type(GenBufList),    Pointer :: List(:)
    Type(CompleteLoop),  Pointer :: Loop(:)
    Integer                      :: iFirst
    Integer                      :: iLast
    ! End of dummy parameters.
    Integer                      :: ybra, yket, mbra, mket
    Integer                      :: iLoop, igen, ioff, i, j
    Type(GenBuffer),     Pointer :: bufptr

    ! Loop over all complete one-body or two-body loops.
    Do iLoop=iFirst, iLast
       ! Loop over all upper partial walks.
       Do i=Flags%DRT(Loop(iLoop)%itop)%iu, Flags%DRT(Loop(iLoop)%itop+1)%iu-1
          ybra = Loop(iLoop)%lexbra + Flags%lexUpper(i)
          yket = Loop(iLoop)%lexket + Flags%lexUpper(i)
          ! Loop over all lower partial walks.
          Do j=1, Flags%DRT(Loop(iLoop)%ibot)%yd(4)
             mbra = Flags%IndVec(ybra + j)
             mket = Flags%IndVec(yket + j)
             ! Store generator matrix element if both CSFs have been selected.
             If (mbra /= 0 .And. mket /= 0) Then
                If (mbra <= mket) Then
                   ioff = mket - mbra
                   ! Allocate new buffer if the current one is full.
                   If (.Not. Associated(List(ioff)%first) .Or. &
                       List(ioff)%first%count == genBufSize) Then
                      Allocate(bufptr)
                      bufptr%count     =  0
                      bufptr%next      => List(ioff)%first
                      List(ioff)%first => bufptr
                      List(ioff)%count =  List(ioff)%count + 1
                   End If
                   bufptr                 => List(ioff)%first
                   igen                   =  bufptr%count + 1
                   bufptr%gen(igen)%value =  Loop(iLoop)%value
                   bufptr%gen(igen)%icol  =  mket
                   bufptr%gen(igen)%p     =  Loop(iLoop)%p
                   bufptr%gen(igen)%q     =  Loop(iLoop)%q
                   bufptr%count           =  igen
                Else
                   ioff = mbra - mket
                   ! Allocate new buffer if the current one is full.
                   If (.Not. Associated(List(ioff)%first) .Or. &
                       List(ioff)%first%count == genBufSize) Then
                      Allocate(bufptr)
                      bufptr%count     =  0
                      bufptr%next      => List(ioff)%first
                      List(ioff)%first => bufptr
                      List(ioff)%count =  List(ioff)%count + 1
                   End If
                   bufptr                 => List(ioff)%first
                   igen                   =  bufptr%count + 1
                   bufptr%gen(igen)%value =  Loop(iLoop)%value
                   bufptr%gen(igen)%icol  =  mbra
                   bufptr%gen(igen)%p     =  Loop(iLoop)%q
                   bufptr%gen(igen)%q     =  Loop(iLoop)%p
                   bufptr%count           =  igen
                End If
             End If
          End Do
       End Do
    End Do
  End Subroutine Complete

End Module gugashape1
