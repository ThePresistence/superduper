!
! This module provides worker routines for the shape-driven GUGA-CI approach
! calculating in a direct manner the one-particle density matrices between
! the specified states for the computation of spectroscopic properties.
!
! Written by Axel Koslowski in 2012-2013 at MPI Muelheim.
!

Module gugashape11

Use gugaglobal
Use gugaloops, Only: Tiny

Implicit None

Private

Public :: ShapeDrivenDirectPropDen

Contains

  ! Assemble upper and lower partial loops and form generator matrix elements,
  ! also from complete loops of the upper and lower part of the DRT. Use the
  ! generator matrix elements to calculate contributions to the one-particle
  ! density matrices for computing spectroscopic properties. This subroutine
  ! is analogous to CalcPropDen in module gugaloops.
  !
  Subroutine ShapeDrivenDirectPropDen(Flags, Perm, Trans, nPerm, mTrans, nTrans, CIC)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Perm(:,:,:), Trans(:,:,:)
    Integer              :: nPerm, mTrans, nTrans
    Double Precision     :: CIC(:,:)
    ! End of dummy parameters.

    ! Initialize density matrices.
    If (nPerm  > 0)  Perm  = 0.0D0
    If (nTrans > 1)  Trans = 0.0D0

    ! Process combinations of upper and lower
    ! partial loops forming one-body loops.
    Call CombineOne(Flags, Perm, Trans, nPerm, mTrans, nTrans, CIC)

    ! Process complete one-body loops in the
    ! upper and lower part of the DRT.
    Call CompleteOneDiag(Flags, Perm, Trans, nPerm, mTrans, nTrans, CIC)
    Call CompleteOneOff (Flags, Perm, Trans, nPerm, mTrans, nTrans, CIC)
  End Subroutine ShapeDrivenDirectPropDen



  ! Process all combinations of upper partial loops (with index l)
  ! and lower partial loops (with index i) forming an off-diagonal
  ! one-electron loop.
  !
  Subroutine CombineOne(Flags, Perm, Trans, nPerm, mTrans, nTrans, CIC)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Perm(:,:,:), Trans(:,:,:)
    Integer              :: nPerm, mTrans, nTrans
    Double Precision     :: CIC(:,:)
    ! End of dummy parameters.
    Integer              :: iCombi, iUp, iLo, lexbra, lexket
    Integer              :: ybra, yket, mbra, mket, i, j, k, l
    Integer              :: ist, jst, ij
    Double Precision     :: value

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
                      If (mbra /= 0) Then
                         mket = Flags%IndVec(yket + j)
                         If (mket /= 0) Then
                            ! Add contributions to density matrices if both CSFs have been selected.
                            Do ist=1, nPerm
                               Perm(k,l,ist) = Perm(k,l,ist) + CIC(mbra,ist) * CIC(mket,ist) * value
                               Perm(l,k,ist) = Perm(l,k,ist) + CIC(mbra,ist) * CIC(mket,ist) * value
                            End Do
                            ij = 1
                            Do ist=1, mTrans
                               Do jst=ist+1, nTrans
                                  Trans(k,l,ij) = Trans(k,l,ij) + CIC(mbra,ist) * CIC(mket,jst) * value
                                  Trans(l,k,ij) = Trans(l,k,ij) + CIC(mket,ist) * CIC(mbra,jst) * value
                                  ij = ij + 1
                               End Do
                            End Do
                         End If
                      End If
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CombineOne



  ! Process complete diagonal one-electron loops.
  !
  Subroutine CompleteOneDiag(Flags, Perm, Trans, nPerm, mTrans, nTrans, CIC)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Perm(:,:,:), Trans(:,:,:)
    Integer              :: nPerm, mTrans, nTrans
    Double Precision     :: CIC(:,:)
    ! End of dummy parameters.
    Integer              :: isym, iLoop, ybra
    Integer              :: mbra, i, j, k
    Integer              :: ist, jst, ij
    Double Precision     :: value, cmi

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
             If (mbra /= 0) Then
                ! Add contributions to density matrices if CSF has been selected.
                Do ist=1, nPerm
                   cmi           = CIC(mbra,ist)
                   Perm(k,k,ist) = Perm(k,k,ist) + cmi * cmi * value
                End Do
                ij = 1
                Do ist=1, mTrans
                   Do jst=ist+1, nTrans
                      Trans(k,k,ij) = Trans(k,k,ij) + CIC(mbra,ist) * CIC(mbra,jst) * value
                      ij = ij + 1
                   End Do
                End Do
             End If
          End Do
       End Do
    End Do
  End Subroutine CompleteOneDiag



  ! Process complete off-diagonal one-electron loops.
  !
  Subroutine CompleteOneOff(Flags, Perm, Trans, nPerm, mTrans, nTrans, CIC)
    Implicit None
    Type(ShavittControl) :: Flags
    Double Precision     :: Perm(:,:,:), Trans(:,:,:)
    Integer              :: nPerm, mTrans, nTrans
    Double Precision     :: CIC(:,:)
    ! End of dummy parameters.
    Integer              :: iLoop, ybra, yket
    Integer              :: isym, mbra, mket, i, j, k, l
    Integer              :: ist, jst, ij
    Double Precision     :: value

    ! Loop over all complete off-diagonal one-body loops.
    Do isym=1, Flags%NumSym
       Do iLoop=Flags%iOff1(isym), Flags%iLast1(isym)
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
                If (mbra /= 0) Then
                   mket = Flags%IndVec(yket + j)
                   If (mket /= 0) Then
                      ! Add contributions to density matrices if both CSFs have been selected.
                      Do ist=1, nPerm
                         Perm(k,l,ist) = Perm(k,l,ist) + CIC(mbra,ist) * CIC(mket,ist) * value
                         Perm(l,k,ist) = Perm(l,k,ist) + CIC(mbra,ist) * CIC(mket,ist) * value
                      End Do
                      ij = 1
                      Do ist=1, mTrans
                         Do jst=ist+1, nTrans
                            Trans(k,l,ij) = Trans(k,l,ij) + CIC(mbra,ist) * CIC(mket,jst) * value
                            Trans(l,k,ij) = Trans(l,k,ij) + CIC(mket,ist) * CIC(mbra,jst) * value
                            ij = ij + 1
                         End Do
                      End Do
                   End If
                End If
             End Do
          End Do
       End Do
    End Do
  End Subroutine CompleteOneOff

End Module gugashape11
