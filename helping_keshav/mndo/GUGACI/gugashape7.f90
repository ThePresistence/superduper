!
! This module provides worker routines for the shape-driven GUGA-CI approach
! calculating and storing the diagonal elements of the CI Hamiltionian for
! states of any relevant symmetry.
!
! Written by Axel Koslowski in 2012-2013 at MPI Muelheim.
!

Module gugashape7

Use gugaglobal
Use gugaloops, Only: Tiny

Implicit None

Private

Public :: ShapeDrivenHii

Contains

  ! Assemble upper and lower partial loops and form generator matrix elements,
  ! also from complete loops of the upper and lower part of the DRT. Use the
  ! generator matrix elements to calculate the diagonal elements of the CI
  ! Hamiltonian storing them into Hii. This routine is analogous to subroutine
  ! CalcHii in module gugaloops which is called with Irrep=0 only.
  !
  Subroutine ShapeDrivenHii(Flags, Hii, OneInt, TwoInt)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Hii(:)
    Double Precision     :: OneInt(:,:)
    Double Precision     :: TwoInt(:,:)
    ! End of dummy parameters.

    ! Initialize diagonal element vector.
    Hii = 0.0D0

    ! Process combinations of diagonal
    ! upper and lower partial loops.
    Call CombineSym22iill(Flags, Hii, TwoInt)
    Call CombineSym22ilil(Flags, Hii, TwoInt)

    ! Process diagonal complete one- and two-body loops in the upper and lower part of the DRT.
    Call Complete(Flags, Hii, OneInt, Flags%OneLoop, Flags%iFirst1(1), Flags%iDiag1(1))
    Call Complete(Flags, Hii, TwoInt, Flags%TwoLoop, Flags%iFirst2(1), Flags%iDiag2(1))
  End Subroutine ShapeDrivenHii



  ! Process all combinations of upper partial loops with identical
  ! indices k and l and lower partial loops with identical indices
  ! i and j forming a diagonal two-electron loop. The corresponding
  ! two-electron integral is (ii|ll).
  !
  Subroutine CombineSym22iill(Flags, Hii, TwoInt)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Hii(:)
    Double Precision     :: TwoInt(:,:)
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo
    Integer              :: lexbra, ybra, mbra
    Integer              :: ii, ll, i, j
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
                help   = TwoInt(ii,ll) * value
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      ! Add contribution to CI matrix diagonal element if CSF has been selected.
                      If (mbra /= 0)  Hii(mbra) = Hii(mbra) + help
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineSym22iill



  ! Process all combinations of upper partial loops with identical
  ! indices k and l and lower partial loops with identical indices
  ! i and j forming a diagonal two-electron loop. The corresponding
  ! two-electron integral is (il|il).
  !
  Subroutine CombineSym22ilil(Flags, Hii, TwoInt)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Hii(:)
    Double Precision     :: TwoInt(:,:)
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo
    Integer              :: lexbra, ybra, mbra
    Integer              :: il, i, j, l
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
                help   = TwoInt(il,il) * value
                ! Loop over all upper partial walks.
                Do i=Flags%DRT(Flags%UpperPart(iUp)%itip)%iu, Flags%DRT(Flags%UpperPart(iUp)%itip+1)%iu-1
                   ybra = lexbra + Flags%lexUpper(i)
                   ! Loop over all lower partial walks.
                   Do j=1, Flags%DRT(Flags%LowerPart(iLo)%itip)%yd(4)
                      mbra = Flags%IndVec(ybra + j)
                      ! Add contribution to CI matrix diagonal element if CSF has been selected.
                      If (mbra /= 0)  Hii(mbra) = Hii(mbra) + help
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineSym22ilil



  ! Process complete diagonal one- or two-electron loops.
  !
  Subroutine Complete(Flags, Hii, MOInt, Loop, iFirst, iDiag)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Hii(:)
    Double Precision     :: MOInt(:,:)
    Type(CompleteLoop)   :: Loop(:)
    Integer              :: iFirst
    Integer              :: iDiag
    ! End of dummy parameters.
    Integer              :: ybra, mbra
    Integer              :: iLoop, i, j
    Double Precision     :: help

    ! Loop over all complete diagonal one-body or two-body loops.
    Do iLoop=iFirst, iDiag
       i    = Loop(iLoop)%p
       j    = Loop(iLoop)%q
       help = MOInt(i,j) * Loop(iLoop)%value
       ! Loop over all upper partial walks.
       Do i=Flags%DRT(Loop(iLoop)%itop)%iu, Flags%DRT(Loop(iLoop)%itop+1)%iu-1
          ybra = Loop(iLoop)%lexbra + Flags%lexUpper(i)
          ! Loop over all lower partial walks.
          Do j=1, Flags%DRT(Loop(iLoop)%ibot)%yd(4)
             mbra = Flags%IndVec(ybra + j)
             ! Add contribution to CI matrix diagonal element if CSF has been selected.
             If (mbra /= 0)  Hii(mbra) = Hii(mbra) + help
          End Do
       End Do
    End Do
  End Subroutine Complete

End Module gugashape7
