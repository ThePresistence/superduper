!
! This module provides worker routines for the shape-driven GUGA-CI approach
! calculating the product of one block of the CI Hamiltonian and a test vector
! of specific symmetry for use with direct CI.
!
! Written by Axel Koslowski in 2012-2013 at MPI Muelheim.
!

Module gugashape4

Use gugaglobal
Use gugaloops, Only: Tiny

Implicit None

Private

Public :: ShapeDrivenProdVecIrrep

Contains

  ! Assemble upper and lower partial loops and form generator matrix elements,
  ! also from complete loops of the upper and lower part of the DRT. Use the
  ! generator matrix elements to calculate contributions to the CI Hamiltonian
  ! and form the corresponding products with the test vector. This subroutine
  ! is analogous to CalcProdVec in module gugaloops with Irrep>0. Arguments
  ! Prod, CIV, and Hii correspond to one block of the CI matrix.
  !
  Subroutine ShapeDrivenProdVecIrrep(Flags, Prod, CIV, Hii, OneInt, TwoInt, Irrep, m0)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Prod(:), CIV(:), Hii(:)
    Double Precision     :: OneInt(:,:), TwoInt(:,:)
    Integer              :: Irrep, m0
    ! End of dummy parameters.
    Integer              :: i

    ! Initialize the product vector with the products of the diagonal elements
    ! of the CI Hamiltonian and the elements of the test vector.
    Do i=1, UBound(Hii,1)
       Prod(i) = Hii(i) * CIV(i)
    End Do

    ! Process combinations of upper and lower partial loops
    ! considering only the off-diagonal part of the CI Hamiltonian.
    ! A subroutine like CombineSym22iljk does not exist because
    ! indices i and j of the lower partial loop have been switched
    ! during partial loop construction such that the corresponding
    ! two-electron integral is (ik|jl) and this case is covered by
    ! CombineSym22ikjl.
    !
    Call CombineOne      (Flags, Prod, CIV, OneInt, Irrep, m0)
    Call CombineSym13ijkl(Flags, Prod, CIV, TwoInt, Irrep, m0)
    Call CombineSym22ijkl(Flags, Prod, CIV, TwoInt, Irrep, m0)
    Call CombineSym22ikjl(Flags, Prod, CIV, TwoInt, Irrep, m0)
    Call CombineSym22half(Flags, Prod, CIV, TwoInt, Irrep, m0)
    Call CombineSym31ijkl(Flags, Prod, CIV, TwoInt, Irrep, m0)

    ! Process complete off-diagonal one- and two-body loops in the upper and lower
    ! part of the DRT. Since the symmetry of both states involved in the matrix
    ! elements must agree, only totally symmetric loops need to be considered.
    Call Complete(Flags, Prod, CIV, OneInt, Flags%OneLoop, Flags%iOff1(1), Flags%iLast1(1), Irrep, m0)
    Call Complete(Flags, Prod, CIV, TwoInt, Flags%TwoLoop, Flags%iOff2(1), Flags%iLast2(1), Irrep, m0)
  End Subroutine ShapeDrivenProdVecIrrep



  ! Process all combinations of upper partial loops (with index l)
  ! and lower partial loops (with index i) forming an off-diagonal
  ! one-electron loop. Only matrix elements between states of
  ! symmetry Irrep are considered. States of identical symmetry are
  ! ordered according to increasing lexical index, so mbra < mket,
  ! in contrast to the analogous subroutine in module gugashape1.
  !
  Subroutine CombineOne(Flags, Prod, CIV, OneInt, Irrep, m0)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Prod(:), CIV(:)
    Double Precision     :: OneInt(:,:)
    Integer              :: Irrep, m0
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo, lexbra, lexket
    Integer              :: ybra, yket, mbra, mket, i, j
    Double Precision     :: value, help

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
                i      = Flags%LowerPart(iLo)%i
                j      = Flags%UpperPart(iUp)%l
                help   = OneInt(i,j) * value
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   yket = lexket + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      ! The dummy element of Flags%SymVec(0) with value -1 will
                      ! make the comparisons fail if mbra == 0 or mket == 0.
                      mbra = Flags%IndVec(ybra + j)
                      If (Flags%SymVec(mbra) == Irrep) Then
                         mket = Flags%IndVec(yket + j)
                         If (Flags%SymVec(mket) == Irrep) Then
                            ! Calculate product if both CSFs have the correct symmetry.
                            mbra       = mbra       - m0
                            mket       = mket       - m0
                            Prod(mbra) = Prod(mbra) + CIV(mket) * help
                            Prod(mket) = Prod(mket) + CIV(mbra) * help
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
  ! same symmetry. Only states of symmetry Irrep are considered. In this
  ! case state index mbra (inferred from lexbra) will be less than mket
  ! (inferred from lexket) because states of equal symmetry are sorted
  ! according to increasing lexical index. During lower partial loop
  ! construction, indices i, j, and k have been rearranged such that
  ! the corresponding two-electron integral is always (ij|kl).
  !
  Subroutine CombineSym13ijkl(Flags, Prod, CIV, TwoInt, Irrep, m0)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Prod(:), CIV(:)
    Double Precision     :: TwoInt(:,:)
    Integer              :: Irrep, m0
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo
    Integer              :: lexbra, lexket
    Integer              :: ybra, yket, mbra, mket
    Integer              :: ij, kl, i, j, k, l
    Double Precision     :: value, help

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
                help   = TwoInt(ij,kl) * value
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   yket = lexket + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      If (Flags%SymVec(mbra) == Irrep) Then
                         mket = Flags%IndVec(yket + j)
                         If (mket /= 0) Then
                            ! Calculate product if both CSFs have the correct symmetry.
                            mbra       = mbra       - m0
                            mket       = mket       - m0
                            Prod(mbra) = Prod(mbra) + CIV(mket) * help
                            Prod(mket) = Prod(mket) + CIV(mbra) * help
                         End If
                      End If
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineSym13ijkl



  ! Process all combinations of upper partial loops with indices k and l
  ! and lower partial loops with indices i and j forming an off-diagonal
  ! two-electron loop with partial bra and ket walk of the same symmetry.
  ! Only states of symmetry Irrep are considered. In this case state index
  ! mbra (inferred from lexbra) will be less than mket (inferred from lexket)
  ! because states of equal symmetry are sorted according to increasing
  ! lexical index. The corresponding two-electron integral is (ij|kl).
  !
  Subroutine CombineSym22ijkl(Flags, Prod, CIV, TwoInt, Irrep, m0)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Prod(:), CIV(:)
    Double Precision     :: TwoInt(:,:)
    Integer              :: Irrep, m0
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo
    Integer              :: lexbra, lexket
    Integer              :: ybra, yket, mbra, mket
    Integer              :: ij, kl, i, j
    Double Precision     :: value, help

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
                ij     = Flags%LowerPart(iLo)%l
                kl     = Flags%UpperPart(iUp)%i
                help   = TwoInt(ij,kl) * value
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   yket = lexket + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      If (Flags%SymVec(mbra) == Irrep) Then
                         mket = Flags%IndVec(yket + j)
                         If (mket /= 0) Then
                            ! Calculate product if both CSFs have the correct symmetry.
                            mbra       = mbra       - m0
                            mket       = mket       - m0
                            Prod(mbra) = Prod(mbra) + CIV(mket) * help
                            Prod(mket) = Prod(mket) + CIV(mbra) * help
                         End If
                      End If
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineSym22ijkl



  ! Process all combinations of upper partial loops with indices k and l
  ! and lower partial loops with indices i and j forming an off-diagonal
  ! two-electron loop with partial bra and ket walk of the same symmetry.
  ! Only states of symmetry Irrep are considered. In this case state index
  ! mbra (inferred from lexbra) will be less than mket (inferred from lexket)
  ! because states of equal symmetry are sorted according to increasing
  ! lexical index. The corresponding two-electron integral is (ik|jl).
  ! By switching indices i and j during lower partial loop construction this
  ! routine covers also those cases in which the integral would be (il|jk).
  !
  Subroutine CombineSym22ikjl(Flags, Prod, CIV, TwoInt, Irrep, m0)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Prod(:), CIV(:)
    Double Precision     :: TwoInt(:,:)
    Integer              :: Irrep, m0
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo
    Integer              :: lexbra, lexket
    Integer              :: ybra, yket, mbra, mket
    Integer              :: ik, jl, i, j, k, l
    Double Precision     :: value, help

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
                help   = TwoInt(ik,jl) * value
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   yket = lexket + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      If (Flags%SymVec(mbra) == Irrep) Then
                         mket = Flags%IndVec(yket + j)
                         If (mket /= 0) Then
                            ! Calculate product if both CSFs have the correct symmetry.
                            mbra       = mbra       - m0
                            mket       = mket       - m0
                            Prod(mbra) = Prod(mbra) + CIV(mket) * help
                            Prod(mket) = Prod(mket) + CIV(mbra) * help
                         End If
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
  ! i and j forming RtRt-RbRb loops for which the corresponding
  ! two-electron integral (il|il) must be divided by two.
  !
  Subroutine CombineSym22half(Flags, Prod, CIV, TwoInt, Irrep, m0)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Prod(:), CIV(:)
    Double Precision     :: TwoInt(:,:)
    Integer              :: Irrep, m0
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo
    Integer              :: lexbra, lexket
    Integer              :: ybra, yket, mbra, mket
    Integer              :: il, i, j, l
    Double Precision     :: value, help

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
                help   = TwoInt(il,il) * value
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   yket = lexket + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      If (Flags%SymVec(mbra) == Irrep) Then
                         mket = Flags%IndVec(yket + j)
                         If (mket /= 0) Then
                            ! Calculate product if both CSFs have the correct symmetry.
                            mbra       = mbra       - m0
                            mket       = mket       - m0
                            Prod(mbra) = Prod(mbra) + CIV(mket) * help
                            Prod(mket) = Prod(mket) + CIV(mbra) * help
                         End If
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
  ! the same symmetry. Only states of symmetry Irrep are considered.
  ! In this case state index mbra (inferred from lexbra) will be less
  ! than mket (inferred from lexket) because states of equal symmetry
  ! are sorted according to increasing lexical index. During upper
  ! partial loop construction, indices j, k, and l have been rearranged
  ! such that the corresponding two-electron integral is always (ij|kl).
  !
  Subroutine CombineSym31ijkl(Flags, Prod, CIV, TwoInt, Irrep, m0)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Prod(:), CIV(:)
    Double Precision     :: TwoInt(:,:)
    Integer              :: Irrep, m0
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo
    Integer              :: lexbra, lexket
    Integer              :: ybra, yket, mbra, mket
    Integer              :: ij, kl, i, j
    Double Precision     :: value, help

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
                help   = TwoInt(ij,kl) * value
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   yket = lexket + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      If (Flags%SymVec(mbra) == Irrep) Then
                         mket = Flags%IndVec(yket + j)
                         If (mket /= 0) Then
                            ! Calculate product if both CSFs have the correct symmetry.
                            mbra       = mbra       - m0
                            mket       = mket       - m0
                            Prod(mbra) = Prod(mbra) + CIV(mket) * help
                            Prod(mket) = Prod(mket) + CIV(mbra) * help
                         End If
                      End If
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineSym31ijkl



  ! Process complete off-diagonal one- or two-electron loops
  ! with partial bra and ket walks of the same symmetry.
  ! Only states of symmetry Irrep are considered.
  !
  Subroutine Complete(Flags, Prod, CIV, MOInt, Loop, iOff, iLast, Irrep, m0)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Prod(:), CIV(:)
    Double Precision     :: MOInt(:,:)
    Type(CompleteLoop)   :: Loop(:)
    Integer              :: iOff, iLast
    Integer              :: Irrep, m0
    ! End of dummy parameters.
    Integer              :: iLoop, ybra, yket
    Integer              :: mbra, mket, i, j
    Double Precision     :: help

    ! Loop over all complete off-diagonal one-body loops
    ! with partial bra and ket walk of the same symmetry.
    Do iLoop=iOff, iLast
       i    = Loop(iLoop)%p
       j    = Loop(iLoop)%q
       help = MOInt(i,j) * Loop(iLoop)%value
       ! Loop over all upper partial walks.
       Do i=Flags%DRT(Loop(iLoop)%itop)%iu, Flags%DRT(Loop(iLoop)%itop+1)%iu-1
          ybra = Loop(iLoop)%lexbra + Flags%lexUpper(i)
          yket = Loop(iLoop)%lexket + Flags%lexUpper(i)
          ! Loop over all lower partial walks.
          Do j=1, Flags%DRT(Loop(iLoop)%ibot)%yd(4)
             mbra = Flags%IndVec(ybra + j)
             If (Flags%SymVec(mbra) == Irrep) Then
                mket = Flags%IndVec(yket + j)
                If (mket /= 0) Then
                   ! Calculate product if both CSFs have the correct symmetry.
                   mbra       = mbra       - m0
                   mket       = mket       - m0
                   Prod(mbra) = Prod(mbra) + CIV(mket) * help
                   Prod(mket) = Prod(mket) + CIV(mbra) * help
                End If
             End If
          End Do
       End Do
    End Do
  End Subroutine Complete

End Module gugashape4
