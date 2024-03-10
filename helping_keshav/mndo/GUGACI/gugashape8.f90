!
! This module provides worker routines for the shape-driven GUGA-CI approach
! calculating in a direct manner the one- and two-particle density matrices of
! a state with symmetry Irrep for the computation of the analytical gradient.
!
! Written by Axel Koslowski in 2012-2013 at MPI Muelheim.
!

Module gugashape8

Use gugaglobal
Use gugaloops, Only: Tiny
Use gugautils, Only: CanonicalIndex

Implicit None

Private

Public :: ShapeDrivenDirectGradDen

Contains

  ! Assemble upper and lower partial loops and form generator matrix elements,
  ! also from complete loops of the upper and lower part of the DRT. Use the
  ! generator matrix elements to calculate contributions to the one- and
  ! two-particle density matrices. This subroutine is analogous to CalcGradDen
  ! in module gugaloops.
  !
  Subroutine ShapeDrivenDirectGradDen(Flags, OneDen, TwoDen, CIV, Irrep)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: OneDen(:,:)
    Double Precision     :: TwoDen(:)
    Double Precision     :: CIV(:)
    Integer              :: Irrep
    ! End of dummy parameters.

    ! Initialize density matrices.
    OneDen = 0.0D0
    TwoDen = 0.0D0

    ! Process combinations of upper and lower partial loops.
    ! A subroutine like CombineSym22iljk does not exist because
    ! indices i and j of the lower partial loop have been switched
    ! during partial loop construction such that the corresponding
    ! two-electron integral is (ik|jl) and this case is covered by
    ! CombineSym22ikjl.
    !
    Call CombineOne      (Flags, OneDen, CIV, Irrep)
    Call CombineSym13ijkl(Flags, TwoDen, CIV, Irrep)
    Call CombineSym22ijkl(Flags, TwoDen, CIV, Irrep)
    Call CombineSym22iill(Flags, TwoDen, CIV, Irrep)
    Call CombineSym22ikjl(Flags, TwoDen, CIV, Irrep)
    Call CombineSym22ilil(Flags, TwoDen, CIV, Irrep)
    Call CombineSym22half(Flags, TwoDen, CIV, Irrep)
    Call CombineSym31ijkl(Flags, TwoDen, CIV, Irrep)

    ! Process complete one- and two-body loops in the upper and lower part of the DRT.
    ! Since the symmetry of both states involved in the matrix elements must be equal
    ! to Irrep, only totally symmetric loops need to be considered.
    Call CompleteOneDiag(Flags, OneDen, CIV, Irrep)
    Call CompleteOneOff (Flags, OneDen, CIV, Irrep)
    Call CompleteTwoDiag(Flags, TwoDen, CIV, Irrep)
    Call CompleteTwoOff (Flags, TwoDen, CIV, Irrep)
  End Subroutine ShapeDrivenDirectGradDen



  ! Process all combinations of upper partial loops (with index l)
  ! and lower partial loops (with index i) forming an off-diagonal
  ! one-electron loop. Only matrix elements between states of
  ! symmetry Irrep are considered. States of identical symmetry are
  ! ordered according to increasing lexical index, so mbra < mket,
  ! in contrast to the analogous subroutine in module gugashape1.
  !
  Subroutine CombineOne(Flags, OneDen, CIV, Irrep)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: OneDen(:,:)
    Double Precision     :: CIV(:)
    Integer              :: Irrep
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo, lexbra, lexket
    Integer              :: ybra, yket, mbra, mket, i, j, k, l
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
                k      = Flags%LowerPart(iLo)%i
                l      = Flags%UpperPart(iUp)%l
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
                            ! Add contribution to density matrix if
                            ! both CSFs have the correct symmetry.
                            help        = CIV(mbra) * CIV(mket) * value
                            OneDen(k,l) = OneDen(k,l) + help
                            OneDen(l,k) = OneDen(l,k) + help
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
  Subroutine CombineSym13ijkl(Flags, TwoDen, CIV, Irrep)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: TwoDen(:)
    Double Precision     :: CIV(:)
    Integer              :: Irrep
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo
    Integer              :: lexbra, lexket
    Integer              :: ybra, yket, mbra, mket
    Integer              :: ijkl, ij, kl, i, j, k, l
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
                ij     = Flags%LowerPart(iLo)%l
                k      = Flags%LowerPart(iLo)%k
                l      = Flags%UpperPart(iUp)%l
                ! l as the only index from the upper partial loop is always largest here.
                kl     = (l*(l-1))/2 + k
                ijkl   = (kl*(kl-1))/2 + ij
                help   = 2.0D0 * value
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
                            ! Add contribution to density matrix if
                            ! both CSFs have the correct symmetry.
                            TwoDen(ijkl) = TwoDen(ijkl) + CIV(mbra) * CIV(mket) * help
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
  Subroutine CombineSym22ijkl(Flags, TwoDen, CIV, Irrep)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: TwoDen(:)
    Double Precision     :: CIV(:)
    Integer              :: Irrep
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo
    Integer              :: lexbra, lexket
    Integer              :: ybra, yket, mbra, mket
    Integer              :: ijkl, ij, kl, i, j
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
                ! kl from the upper partial loop is larger than ij.
                ijkl   = (kl*(kl-1))/2 + ij
                help   = 2.0D0 * value
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
                            ! Add contribution to density matrix if
                            ! both CSFs have the correct symmetry.
                            TwoDen(ijkl) = TwoDen(ijkl) + CIV(mbra) * CIV(mket) * help
                         End If
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
  Subroutine CombineSym22iill(Flags, TwoDen, CIV, Irrep)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: TwoDen(:)
    Double Precision     :: CIV(:)
    Integer              :: Irrep
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo
    Integer              :: lexbra, ybra, mbra
    Integer              :: iill, ii, ll, i, j
    Double Precision     :: value, help

    ! Loop over all combinations.
    Do iCombi=1, Flags%NumSym22iill
       ! Loop over all upper partial loops.
       Do iUp=Flags%CombiListSym22iill(iCombi)%iUpB, Flags%CombiListSym22iill(iCombi)%iUpE
          ! Loop over all lower partial loops.
          Do iLo=Flags%CombiListSym22iill(iCombi)%iLoB, Flags%CombiListSym22iill(iCombi)%iLoE
             value = Flags%UpperPart(iUp)%singlet * Flags%LowerPart(iLo)%singlet
             If (Dabs(value) > Tiny) Then
                lexbra = Flags%UpperPart(iUp)%lexbra + Flags%LowerPart(iLo)%lexbra
                ii     = Flags%LowerPart(iLo)%l
                ll     = Flags%UpperPart(iUp)%i
                ! ll from the upper partial loop is larger then ii.
                iill   = (ll*(ll-1))/2 + ii
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      If (Flags%SymVec(mbra) == Irrep) Then
                         ! Add contribution to density matrix
                         ! if CSF has the correct symmetry.
                         help         = CIV(mbra)
                         TwoDen(iill) = TwoDen(iill) + help * help * value
                      End If
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineSym22iill



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
  Subroutine CombineSym22ikjl(Flags, TwoDen, CIV, Irrep)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: TwoDen(:)
    Double Precision     :: CIV(:)
    Integer              :: Irrep
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo
    Integer              :: lexbra, lexket
    Integer              :: ybra, yket, mbra, mket
    Integer              :: ikjl, ik, jl, i, j, k, l
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
                ikjl   = CanonicalIndex(ik,jl)
                help   = 2.0D0 * value
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
                            ! Add contribution to density matrix if
                            ! both CSFs have the correct symmetry.
                            TwoDen(ikjl) = TwoDen(ikjl) + CIV(mbra) * CIV(mket) * help
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
  ! i and j forming a diagonal two-electron loop. The corresponding
  ! two-electron integral is (il|il).
  !
  Subroutine CombineSym22ilil(Flags, TwoDen, CIV, Irrep)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: TwoDen(:)
    Double Precision     :: CIV(:)
    Integer              :: Irrep
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo
    Integer              :: lexbra, ybra, mbra
    Integer              :: ilil, il, i, j, l
    Double Precision     :: value, help

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
                ilil   = (il*(il+1))/2
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      If (Flags%SymVec(mbra) == Irrep) Then
                         ! Add contribution to density matrix
                         ! if CSF has the correct symmetry.
                         help         = CIV(mbra)
                         TwoDen(ilil) = TwoDen(ilil) + help * help * value
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
  Subroutine CombineSym22half(Flags, TwoDen, CIV, Irrep)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: TwoDen(:)
    Double Precision     :: CIV(:)
    Integer              :: Irrep
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo
    Integer              :: lexbra, lexket
    Integer              :: ybra, yket, mbra, mket
    Integer              :: ilil, il, i, j, l
    Double Precision     :: value

    ! Loop over all combinations.
    Do iCombi=1, Flags%NumSym22half
       ! Loop over all upper partial loops.
       Do iUp=Flags%CombiListSym22half(iCombi)%iUpB, Flags%CombiListSym22half(iCombi)%iUpE
          ! Loop over all lower partial loops.
          Do iLo=Flags%CombiListSym22half(iCombi)%iLoB, Flags%CombiListSym22half(iCombi)%iLoE
             ! The factor of 0.5 cancels with the factor 2.0 below, so both factors are left out.
             value = Flags%UpperPart(iUp)%singlet * Flags%LowerPart(iLo)%singlet + &
                     Flags%UpperPart(iUp)%triplet * Flags%LowerPart(iLo)%triplet
             If (Dabs(value) > Tiny) Then
                lexbra = Flags%UpperPart(iUp)%lexbra + Flags%LowerPart(iLo)%lexbra
                lexket = Flags%UpperPart(iUp)%lexket + Flags%LowerPart(iLo)%lexket
                i      = Flags%LowerPart(iLo)%i
                l      = Flags%UpperPart(iUp)%l
                ! l from the upper partial loop is always greater
                ! than i from the lower partial loop here.
                il     = (l*(l-1))/2 + i
                ilil   = (il*(il+1))/2
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
                            ! Add contribution to density matrix if
                            ! both CSFs have the correct symmetry.
                            TwoDen(ilil) = TwoDen(ilil) + CIV(mbra) * CIV(mket) * value
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
  Subroutine CombineSym31ijkl(Flags, TwoDen, CIV, Irrep)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: TwoDen(:)
    Double Precision     :: CIV(:)
    Integer              :: Irrep
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo
    Integer              :: lexbra, lexket
    Integer              :: ybra, yket, mbra, mket
    Integer              :: ijkl, ij, kl, i, j
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
                ijkl   = CanonicalIndex(ij,kl)
                help   = 2.0D0 * value
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
                            ! Add contribution to density matrix if
                            ! both CSFs have the correct symmetry.
                            TwoDen(ijkl) = TwoDen(ijkl) + CIV(mbra) * CIV(mket) * help
                         End If
                      End If
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineSym31ijkl



  ! Process complete diagonal one-electron loops.
  ! Only states of symmetry Irrep are considered.
  !
  Subroutine CompleteOneDiag(Flags, OneDen, CIV, Irrep)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: OneDen(:,:)
    Double Precision     :: CIV(:)
    Integer              :: Irrep
    ! End of dummy parameters.
    Integer              :: iLoop, ybra
    Integer              :: mbra, i, j, k
    Double Precision     :: value, help

    ! Loop over all complete diagonal one-body loops.
    Do iLoop=Flags%iFirst1(1), Flags%iDiag1(1)
       k     = Flags%OneLoop(iLoop)%p
       value = Flags%OneLoop(iLoop)%value
       ! Loop over all upper partial walks.
       Do i=Flags%DRT(Flags%OneLoop(iLoop)%itop)%iu, Flags%DRT(Flags%OneLoop(iLoop)%itop+1)%iu-1
          ybra = Flags%OneLoop(iLoop)%lexbra + Flags%lexUpper(i)
          ! Loop over all lower partial walks.
          Do j=1, Flags%DRT(Flags%OneLoop(iLoop)%ibot)%yd(4)
             mbra = Flags%IndVec(ybra + j)
             If (Flags%SymVec(mbra) == Irrep) Then
                ! Add contribution to density matrix
                ! if CSF has the correct symmetry.
                help        = CIV(mbra)
                OneDen(k,k) = OneDen(k,k) + help * help * value
             End If
          End Do
       End Do
    End Do
  End Subroutine CompleteOneDiag



  ! Process complete off-diagonal one-electron loops
  ! with partial bra and ket walks of the same symmetry.
  ! Only states of symmetry Irrep are considered.
  !
  Subroutine CompleteOneOff(Flags, OneDen, CIV, Irrep)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: OneDen(:,:)
    Double Precision     :: CIV(:)
    Integer              :: Irrep
    ! End of dummy parameters.
    Integer              :: iLoop, ybra, yket
    Integer              :: mbra, mket, i, j, k, l
    Double Precision     :: value, help

    ! Loop over all complete off-diagonal one-body loops
    ! with partial bra and ket walk of the same symmetry.
    Do iLoop=Flags%iOff1(1), Flags%iLast1(1)
       k     = Flags%OneLoop(iLoop)%p
       l     = Flags%OneLoop(iLoop)%q
       value = Flags%OneLoop(iLoop)%value
       ! Loop over all upper partial walks.
       Do i=Flags%DRT(Flags%OneLoop(iLoop)%itop)%iu, Flags%DRT(Flags%OneLoop(iLoop)%itop+1)%iu-1
          ybra = Flags%OneLoop(iLoop)%lexbra + Flags%lexUpper(i)
          yket = Flags%OneLoop(iLoop)%lexket + Flags%lexUpper(i)
          ! Loop over all lower partial walks.
          Do j=1, Flags%DRT(Flags%OneLoop(iLoop)%ibot)%yd(4)
             mbra = Flags%IndVec(ybra + j)
             If (Flags%SymVec(mbra) == Irrep) Then
                mket = Flags%IndVec(yket + j)
                If (mket /= 0) Then
                   ! Add contribution to density matrix if
                   ! both CSFs have the correct symmetry.
                   help        = CIV(mbra) * CIV(mket) * value
                   OneDen(k,l) = OneDen(k,l) + help
                   OneDen(l,k) = OneDen(l,k) + help
                End If
             End If
          End Do
       End Do
    End Do
  End Subroutine CompleteOneOff



  ! Process complete diagonal two-electron loops.
  ! Only states of symmetry Irrep are considered.
  !
  Subroutine CompleteTwoDiag(Flags, TwoDen, CIV, Irrep)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: TwoDen(:)
    Double Precision     :: CIV(:)
    Integer              :: Irrep
    ! End of dummy parameters.
    Integer              :: iLoop, ybra, mbra
    Integer              :: ijkl, ij, kl, i, j
    Double Precision     :: value, help

    ! Loop over all complete diagonal two-body loops.
    Do iLoop=Flags%iFirst2(1), Flags%iDiag2(1)
       ij    = Flags%TwoLoop(iLoop)%p
       kl    = Flags%TwoLoop(iLoop)%q
       ijkl  = CanonicalIndex(ij,kl)
       value = Flags%TwoLoop(iLoop)%value
       ! Loop over all upper partial walks.
       Do i=Flags%DRT(Flags%TwoLoop(iLoop)%itop)%iu, Flags%DRT(Flags%TwoLoop(iLoop)%itop+1)%iu-1
          ybra = Flags%TwoLoop(iLoop)%lexbra + Flags%lexUpper(i)
          ! Loop over all lower partial walks.
          Do j=1, Flags%DRT(Flags%TwoLoop(iLoop)%ibot)%yd(4)
             mbra = Flags%IndVec(ybra + j)
             If (Flags%SymVec(mbra) == Irrep) Then
                ! Add contribution to density matrix
                ! if CSF has the correct symmetry.
                help         = CIV(mbra)
                TwoDen(ijkl) = TwoDen(ijkl) + help * help * value
             End If
          End Do
       End Do
    End Do
  End Subroutine CompleteTwoDiag



  ! Process complete off-diagonal two-electron loops
  ! with partial bra and ket walks of the same symmetry.
  ! Only states of symmetry Irrep are considered.
  !
  Subroutine CompleteTwoOff(Flags, TwoDen, CIV, Irrep)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: TwoDen(:)
    Double Precision     :: CIV(:)
    Integer              :: Irrep
    ! End of dummy parameters.
    Integer              :: iLoop, ybra, yket, mbra, mket
    Integer              :: ijkl, ij, kl, i, j
    Double Precision     :: help

    ! Loop over all complete off-diagonal two-body loops
    ! with partial bra and ket walk of the same symmetry.
    Do iLoop=Flags%iOff2(1), Flags%iLast2(1)
       ij   = Flags%TwoLoop(iLoop)%p
       kl   = Flags%TwoLoop(iLoop)%q
       ijkl = CanonicalIndex(ij,kl)
       help = 2.0D0 * Flags%TwoLoop(iLoop)%value
       ! Loop over all upper partial walks.
       Do i=Flags%DRT(Flags%TwoLoop(iLoop)%itop)%iu, Flags%DRT(Flags%TwoLoop(iLoop)%itop+1)%iu-1
          ybra = Flags%TwoLoop(iLoop)%lexbra + Flags%lexUpper(i)
          yket = Flags%TwoLoop(iLoop)%lexket + Flags%lexUpper(i)
          ! Loop over all lower partial walks.
          Do j=1, Flags%DRT(Flags%TwoLoop(iLoop)%ibot)%yd(4)
             mbra = Flags%IndVec(ybra + j)
             If (Flags%SymVec(mbra) == Irrep) Then
                mket = Flags%IndVec(yket + j)
                If (mket /= 0) Then
                   ! Add contribution to density matrix if
                   ! both CSFs have the correct symmetry.
                   TwoDen(ijkl) = TwoDen(ijkl) + CIV(mbra) * CIV(mket) * help
                End If
             End If
          End Do
       End Do
    End Do
  End Subroutine CompleteTwoOff

End Module gugashape8
