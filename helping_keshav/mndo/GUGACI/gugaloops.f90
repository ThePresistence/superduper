!
! This module provides subroutines for the loop construction in GUGA-CI.
! It uses the Hamiltonian repartitioning scheme (I. Shavitt, Lecture
! Notes in Chemistry 22 (1981), 51-99) with modification of the WW
! segment values to be consistent with the conventions in the analytical
! derivative section by Serguei Patchkovskii.
!
! Written by Axel Koslowski in 2000-2005 at MPI Muelheim.
!

Module gugaloops

Use gugaglobal
Use gugautils, Only: CanonicalIndex

Implicit None

Private

Public :: CalcGenerators       ! Calculate and store coupling coefficients.
Public :: CalcProdVec          ! Calculate CI Hamiltonian matrix times CI vector.
Public :: CalcProdMat          ! Calculate CI Hamiltonian matrix times CI coefficient matrix.
Public :: CalcHii              ! Calculate and store diagonal elements of CI Hamiltonian.
Public :: CalcGradDen          ! Calculate densities needed for the analytical CI gradient.
Public :: CalcNACDen           ! Calculate transition densities needed for the non-adiabatic
                               ! coupling vectors.
Public :: CalcPropDen          ! Calculate densities needed for the spectroscopic properties.
Public :: FuncA
Public :: FuncB
Public :: FuncC
Public :: FuncD
Public :: Sqrt2
Public :: t
Public :: Tiny

!------------------------------------------------------------------------------
!------  Global parameters  ---------------------------------------------------
!------------------------------------------------------------------------------

  Double Precision, Parameter :: Sqrt2 = 1.4142135623730950488D0
  Double Precision, Parameter :: t     = 0.7071067811865475244D0
  Double Precision, Parameter :: Tiny  = 1.0D-10


!------------------------------------------------------------------------------
!------  Global variables  ----------------------------------------------------
!------------------------------------------------------------------------------

  ! DRT, IndVec, SymVec:  Pointers referencing the corresponding elements in
  !                       ShavittFlags.
  !
  ! UpperVec:             Weights of the upper walks leading to the top vertex
  !                       of the loop.
  !
  ! OneList, TwoList:     Buffer lists for the one and two-electron generator
  !                       matrix elements.
  !
  ! OneInt, TwoInt:       One- and two-electron MO integrals.
  !
  ! Occ:                  Occupation numbers of active orbitals.
  !
  ! CIV1, CIV2, CIP:      CI vectors and H*CIV1 product vector.
  !
  ! CIC:                  CI coefficient matrix.
  !
  ! CICtr, CIPtr:         Transposed CI coefficient and H*C product matrices.
  !
  ! Hdiag:                Diagonal elements of CI Hamiltonian.
  !
  ! OneDen, TwoDen:       One- and two-particle density.
  !
  ! Perm:                 Permanent one-electron densities.
  !
  ! Trans:                One-electron transition density
  !                       (square matrix packed linearly).
  !
  ! nPerm:                Number of states in Perm.
  !
  ! mTrans, nTrans:       Highest states for transitions (mTrans -> nTrans).
  !
  ! nupper:               Number of entries in the upper walk vector.
  !
  ! job:                  What to do with the coupling coefficients.
  !                       (See OneBodyLoopComplete and TwoBodyLoopComplete.)
  !
  ! oneOnly:              Generate one-electron loops only.
  !
  ! diagOnly:             Generate loops corresponding to diagonal elements only.
  !
  ! allOneGen:            Store all one-electron generator matrix elements.
  !
  ! symFlag:              Store one- and two-electron generator matrix elements
  !                       involving Gelfand states of this symmetry.
  !
  ! index0:               Index to subtract from the indices of the global
  !                       CI problem for conversion to the local CI problem.
  !

  Type(PaldusRow),  Dimension(:),     Pointer           :: DRT
  Integer,          Dimension(:),     Pointer           :: IndVec, SymVec
  Integer,          Dimension(:),     Allocatable       :: UpperVec
  Type(GenBufList), Dimension(:),     Pointer           :: OneList, TwoList
  Double Precision, Dimension(:,:),   Pointer           :: OneInt, TwoInt
  Double Precision, Dimension(:),     Allocatable       :: Occ
  Double Precision, Dimension(:),     Pointer,     Save :: CIV1, CIV2, CIP, Hdiag
  Double Precision, Dimension(:,:),   Pointer,     Save :: CIC, CICtr, CIPtr
  Double Precision, Dimension(:,:),   Pointer           :: OneDen
  Double Precision, Dimension(:),     Pointer           :: TwoDen
  Double Precision, Dimension(:,:,:), Pointer           :: Perm, Trans
  Integer                                               :: nPerm, mTrans, nTrans
  Integer                                               :: nupper, job
  Logical                                               :: oneOnly, diagOnly
  Logical                                               :: allOneGen, allTwoGen
  Logical,          Dimension(8)                        :: symFlag
  Integer                                               :: index0

Contains

!------------------------------------------------------------------------------
!------  Driver routines  -----------------------------------------------------
!------------------------------------------------------------------------------

  ! The following subroutines are the public driver routines for the generation
  ! of the unitary group generator matrix elements in a Gel'fand-Zetlin basis
  ! for electronic structure calculations. The basic procedure is described in:
  !
  ! F. A. Matsen, R. Pauncz, The Unitary Group in Quantum Chemistry, Elsevier
  ! Amsterdam (1986). [Obsoleted by the next reference.]
  !
  ! R. Pauncz, The Symmetric Group in Quantum Chemistry, CRC Boca Raton FL (1995).
  !
  ! I. Shavitt, The graphical unitary group approach and its application to
  ! direct configuration interaction calculations,
  ! Lecture Notes in Chemistry 22 (1981), 51-99.
  !
  ! Our implementation works by constructing all possible one- and two-electron
  ! loops starting at an arbitary vertex in the Shavitt graph, simply by adding
  ! segments in all possible ways. After a loop has been completed, all upper
  ! and lower walks are determined and the value of the coupling coefficient
  ! is used for as many CI matrix elements as possible. This procedure is done
  ! for all vertices in the graph. For all loops starting at the same vertex,
  ! the determination of the upper walks is done in advance, and the corresponding
  ! partial walk weights are stored in the array UpperVec.
  !


  ! List1 and List2 have to be properly initialized
  ! before calling subroutine CalcGenerators!
  !
  Subroutine CalcGenerators(Flags, List1, List2, Occup, Irrep)
    Implicit None
    Type(ShavittControl)                   :: Flags
    Type(GenBufList), Dimension(:), Target :: List1, List2
    Double Precision, Dimension(:)         :: Occup
    Integer                                :: Irrep
    ! End of dummy parameters.

    job       =  0
    oneOnly   =  .False.
    diagOnly  =  .False.
    allTwoGen =  .False.
    OneList   => List1
    TwoList   => List2

    If (Irrep > 0) Then
       allOneGen      = .False.
       symFlag        = .False.
       symFlag(Irrep) = .True.
    Else
       allOneGen      = .True.
       symFlag        = .True.
    End If

    Call StartLoops(Flags, Occup)

    Return
  End Subroutine CalcGenerators



  Subroutine CalcProdVec(Flags, Prod, CIVec, MOInt1, MOInt2, Occup, Irrep, i0)
    Implicit None
    Type(ShavittControl)                     :: Flags
    Double Precision, Dimension(:),   Target :: Prod, CIVec
    Double Precision, Dimension(:,:), Target :: MOInt1, MOInt2
    Double Precision, Dimension(:)           :: Occup
    Integer                                  :: Irrep, i0

    job       =  1
    oneOnly   =  .False.
    diagOnly  =  .False.
    allOneGen =  .False.
    allTwoGen =  .False.
    CIP       => Prod
    CIV1      => CIVec
    OneInt    => MOInt1
    TwoInt    => MOInt2
    index0    =  i0
    CIP       =  0.D0

    If (Irrep > 0) Then
       symFlag        = .False.
       symFlag(Irrep) = .True.
    Else
       symFlag        = .True.
    End If

    Call StartLoops(Flags, Occup)

    Return
  End Subroutine CalcProdVec



  Subroutine CalcProdMat(Flags, ProdTrans, CICTrans, MOInt1, MOInt2, Occup, Irrep, i0)
    Implicit None
    Type(ShavittControl)                     :: Flags
    Double Precision, Dimension(:,:), Target :: ProdTrans, CICTrans
    Double Precision, Dimension(:,:), Target :: MOInt1, MOInt2
    Double Precision, Dimension(:)           :: Occup
    Integer                                  :: Irrep, i0

    job       =  2
    oneOnly   =  .False.
    diagOnly  =  .False.
    allOneGen =  .False.
    allTwoGen =  .False.
    CIPtr     => ProdTrans
    CICtr     => CICTrans
    OneInt    => MOInt1
    TwoInt    => MOInt2
    index0    =  i0
    CIPtr     =  0.D0

    If (Irrep > 0) Then
       symFlag        = .False.
       symFlag(Irrep) = .True.
    Else
       symFlag        = .True.
    End If

    Call StartLoops(Flags, Occup)

    Return
  End Subroutine CalcProdMat



  Subroutine CalcHii(Flags, Hii, MOInt1, MOInt2, Occup, Irrep)
    Implicit None
    Type(ShavittControl)                     :: Flags
    Double Precision, Dimension(:),   Target :: Hii
    Double Precision, Dimension(:,:), Target :: MOInt1, MOInt2
    Double Precision, Dimension(:)           :: Occup
    Integer                                  :: Irrep

    job       =  3
    oneOnly   =  .False.
    diagOnly  =  .True.
    allOneGen =  .False.
    allTwoGen =  .False.
    Hdiag     => Hii
    OneInt    => MOInt1
    TwoInt    => MOInt2
    Hii       =  0.D0

    If (Irrep > 0) Then
       symFlag        = .False.
       symFlag(Irrep) = .True.
    Else
       symFlag        = .True.
    End If

    Call StartLoops(Flags, Occup)

    Return
  End Subroutine CalcHii



  Subroutine CalcGradDen(Flags, OneBodyDen, TwoBodyDen, CIVec, Occup, Irrep)
    Implicit None
    Type(ShavittControl)                     :: Flags
    Double Precision, Dimension(:,:), Target :: OneBodyDen
    Double Precision, Dimension(:),   Target :: TwoBodyDen
    Double Precision, Dimension(:),   Target :: CIVec
    Double Precision, Dimension(:)           :: Occup
    Integer                                  :: Irrep

    job       =  4
    oneOnly   =  .False.
    diagOnly  =  .False.
    allOneGen =  .False.
    allTwoGen =  .False.
    OneDen    => OneBodyDen
    TwoDen    => TwoBodyDen
    CIV1      => CIVec
    CIV2      => CIVec
    OneDen    =  0.D0
    TwoDen    =  0.D0

    If (Irrep > 0) Then
       symFlag        = .False.
       symFlag(Irrep) = .True.
    Else
       symFlag        = .True.
    End If

    Call StartLoops(Flags, Occup)

    Return
  End Subroutine CalcGradDen



  Subroutine CalcNACDen(Flags, OneBodyDen, TwoBodyDen, CIVec1, CIVec2, Occup)
    Implicit None
    Type(ShavittControl)                     :: Flags
    Double Precision, Dimension(:,:), Target :: OneBodyDen
    Double Precision, Dimension(:),   Target :: TwoBodyDen
    Double Precision, Dimension(:),   Target :: CIVec1, CIVec2
    Double Precision, Dimension(:)           :: Occup

    job       =  4
    oneOnly   =  .False.
    diagOnly  =  .False.
    allOneGen =  .True.
    allTwoGen =  .True.
    OneDen    => OneBodyDen
    TwoDen    => TwoBodyDen
    CIV1      => CIVec1
    CIV2      => CIVec2
    OneDen    =  0.D0
    TwoDen    =  0.D0
    symFlag   = .True.

    Call StartLoops(Flags, Occup)

    Return
  End Subroutine CalcNACDen



  Subroutine CalcPropDen(Flags, PermDen, TransDen, nPermDen, mTransDen, nTransDen, CICoeff, Occup, Irrep)
    Implicit None
    Type(ShavittControl)                       :: Flags
    Double Precision, Dimension(:,:,:), Target :: PermDen, TransDen
    Integer                                    :: nPermDen, mTransDen, nTransDen
    Double Precision, Dimension(:,:),   Target :: CICoeff
    Double Precision, Dimension(:)             :: Occup
    Integer                                    :: Irrep

    job       =  5
    oneOnly   =  .True.
    diagOnly  =  .False.
    allTwoGen =  .False.
    Perm      => PermDen
    Trans     => TransDen
    nPerm     =  nPermDen
    mTrans    =  mTransDen
    nTrans    =  nTransDen
    CIC       => CICoeff

    If (Irrep > 0) Then
       allOneGen      = .False.
       symFlag        = .False.
       symFlag(Irrep) = .True.
    Else
       allOneGen      = .True.
       symFlag        = .True.
    End If

    If (nPerm > 0) Then
       Perm = 0.D0
    End If

    If (nTrans > 1) Then
       Trans = 0.D0
    End If

    Call StartLoops(Flags, Occup)

    Return
  End Subroutine CalcPropDen


!------------------------------------------------------------------------------
!------  This is the entry point for the loop construction  -------------------
!------------------------------------------------------------------------------

  Subroutine StartLoops(Flags, Occup)
    Implicit None
    Type(ShavittControl)           :: Flags
    Double Precision, Dimension(:) :: Occup
    ! End of dummy parameters.
    Integer                        :: itop

    DRT    => Flags%DRT
    IndVec => Flags%IndVec
    SymVec => Flags%SymVec

    Allocate(UpperVec(DRT(1)%yd(4)))
    Allocate(Occ(Size(Occup)))
    Occ = Occup
    Do itop=1, Size(DRT)-1
       nupper = 0
       Call UpperWalks(itop, 0)
       If (.Not. oneOnly)  Call StartWWLoop(itop)
       Call StartWLoop(itop)
       If (DRT(itop)%k <= 1)  Cycle
       If (.Not. oneOnly)  Call StartRtLtLoop(itop)
       If (diagOnly)  Cycle
       Call StartRtLoop(itop)
       If (oneOnly)  Cycle
       Call StartWRtLoop(itop)
       Call StartRtRtLoop(itop)
    End Do
    Deallocate(Occ)
    Deallocate(UpperVec)

    Return
  End Subroutine StartLoops


!------------------------------------------------------------------------------
!------  Six different ways to start a loop  ----------------------------------
!------------------------------------------------------------------------------

  Subroutine StartWWLoop(itop)
    Implicit None
    Integer          :: itop
    ! End of dummy parameters.
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3, h

    j0  = DRT(itop)%jd(0)
    j1  = DRT(itop)%jd(1)
    j2  = DRT(itop)%jd(2)
    j3  = DRT(itop)%jd(3)

    h   = DRT(itop)%k
    uu  = Occ(h)

    If (j0 /= 0)  Call AppendWW(j0, DRT(itop)%yd(0), 0.5D0* uu      * uu,       h)
    If (j1 /= 0)  Call AppendWW(j1, DRT(itop)%yd(1), 0.5D0* uu      *(uu-2.D0), h)
    If (j2 /= 0)  Call AppendWW(j2, DRT(itop)%yd(2), 0.5D0* uu      *(uu-2.D0), h)
    If (j3 /= 0)  Call AppendWW(j3, DRT(itop)%yd(3), 0.5D0*(uu-2.D0)*(uu-2.D0), h)

    Return
  End Subroutine StartWWLoop



  Subroutine StartWLoop(itop)
    Implicit None
    Integer          :: itop
    ! End of dummy parameters.
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3, h

    j0  = DRT(itop)%jd(0)
    j1  = DRT(itop)%jd(1)
    j2  = DRT(itop)%jd(2)
    j3  = DRT(itop)%jd(3)

    h   = DRT(itop)%k
    uu  = Occ(h)

    If (j0 /= 0)  Call AppendUpperW(j0, DRT(itop)%yd(0),     -uu, h)
    If (j1 /= 0)  Call AppendUpperW(j1, DRT(itop)%yd(1), 1.D0-uu, h)
    If (j2 /= 0)  Call AppendUpperW(j2, DRT(itop)%yd(2), 1.D0-uu, h)
    If (j3 /= 0)  Call AppendUpperW(j3, DRT(itop)%yd(3), 2.D0-uu, h)

    Return
  End Subroutine StartWLoop



  Subroutine StartWRtLoop(itop)
    Implicit None
    Integer          :: itop
    ! End of dummy parameters.
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3
    Integer          :: y0, y1, y2, y3
    Integer          :: b, h

    j0 = DRT(itop)%jd(0)
    j1 = DRT(itop)%jd(1)
    j2 = DRT(itop)%jd(2)
    j3 = DRT(itop)%jd(3)

    y0 = DRT(itop)%yd(0)
    y1 = DRT(itop)%yd(1)
    y2 = DRT(itop)%yd(2)
    y3 = DRT(itop)%yd(3)

    b  = DRT(itop)%b
    h  = DRT(itop)%k

    uu = Occ(h)

    If (j0 /= 0 .And. j1 /= 0)  Call AppendLowerR(j0, j1, y0, y1,      -0.5D0*uu,               h, h, h)  ! AppendWRt
    If (j0 /= 0 .And. j2 /= 0)  Call AppendLowerR(j0, j2, y0, y2,      -0.5D0*uu,               h, h, h)  ! AppendWRt
    If (j1 /= 0 .And. j3 /= 0)  Call AppendLowerR(j1, j3, y1, y3, (1.D0-0.5D0*uu)*FuncA(b,0,1), h, h, h)  ! AppendWRt
    If (j2 /= 0 .And. j3 /= 0)  Call AppendLowerR(j2, j3, y2, y3, (1.D0-0.5D0*uu)*FuncA(b,2,1), h, h, h)  ! AppendWRt

    Return
  End Subroutine StartWRtLoop



  Subroutine StartRtRtLoop(itop)
    Implicit None
    Integer :: itop
    ! End of dummy parameters.
    Integer :: j0, j3, y0, y3, h

    j0 = DRT(itop)%jd(0)
    j3 = DRT(itop)%jd(3)

    y0 = DRT(itop)%yd(0)
    y3 = DRT(itop)%yd(3)

    h  = DRT(itop)%k

    If (j0 /= 0 .And. j3 /= 0)  Call AppendRR(j0, j3, y0, y3, 1.D0, Sqrt2, 0.D0, h, h, .true.)  ! AppendRtRt

    Return
  End Subroutine StartRtRtLoop



  Subroutine StartRtLtLoop(itop)
    Implicit None
    Integer          :: itop
    ! End of dummy parameters.
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3
    Integer          :: y0, y1, y2, y3
    Integer          :: b, h

    j0 = DRT(itop)%jd(0)
    j1 = DRT(itop)%jd(1)
    j2 = DRT(itop)%jd(2)
    j3 = DRT(itop)%jd(3)

    y0 = DRT(itop)%yd(0)
    y1 = DRT(itop)%yd(1)
    y2 = DRT(itop)%yd(2)
    y3 = DRT(itop)%yd(3)

    b  = DRT(itop)%b
    h  = DRT(itop)%k

    uu = Occ(h)

    If (j0 /= 0)  Call AppendRL(j0, j0, y0, y0, 1.D0, -t*      uu,     0.D0,          h, h, .false.)  ! AppendRtLt
    If (j1 /= 0)  Call AppendRL(j1, j1, y1, y1, 1.D0,  t*(1.D0-uu), -t*FuncA(b,-1,1), h, h, .false.)  ! AppendRtLt
    If (j1 /= 0 .And. j2 /= 0 .And. .Not. diagOnly) &
    &             Call AppendRL(j1, j2, y1, y2, 1.D0,     0.D0,        1.D0,          h, h, .true. )  ! AppendRtLt
    If (j2 /= 0)  Call AppendRL(j2, j2, y2, y2, 1.D0,  t*(1.D0-uu),  t*FuncA(b, 3,1), h, h, .false.)  ! AppendRtLt
    If (j3 /= 0)  Call AppendRL(j3, j3, y3, y3, 1.D0,  t*(2.D0-uu),    0.D0,          h, h, .false.)  ! AppendRtLt

    Return
  End Subroutine StartRtLtLoop



  Subroutine StartRtLoop(itop)
    Implicit None
    Integer :: itop
    ! End of dummy parameters.
    Integer :: j0, j1, j2, j3
    Integer :: y0, y1, y2, y3
    Integer :: b, h

    j0 = DRT(itop)%jd(0)
    j1 = DRT(itop)%jd(1)
    j2 = DRT(itop)%jd(2)
    j3 = DRT(itop)%jd(3)

    y0 = DRT(itop)%yd(0)
    y1 = DRT(itop)%yd(1)
    y2 = DRT(itop)%yd(2)
    y3 = DRT(itop)%yd(3)

    b  = DRT(itop)%b
    h  = DRT(itop)%k

    If (j0 /= 0 .And. j1 /= 0)  Call AppendUpperR(j0, j1, y0, y1, 1.D0,         h)  ! AppendRt
    If (j0 /= 0 .And. j2 /= 0)  Call AppendUpperR(j0, j2, y0, y2, 1.D0,         h)  ! AppendRt
    If (j1 /= 0 .And. j3 /= 0)  Call AppendUpperR(j1, j3, y1, y3, FuncA(b,0,1), h)  ! AppendRt
    If (j2 /= 0 .And. j3 /= 0)  Call AppendUpperR(j2, j3, y2, y3, FuncA(b,2,1), h)  ! AppendRt

    Return
  End Subroutine StartRtLoop


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the top segments  ----------------
!------------------------------------------------------------------------------

  Subroutine AppendWW(irow, y, value, l)
    Implicit None
    Integer          :: irow, y
    Double Precision :: value
    Integer          :: l
    ! End of dummy parameters.

    If (Dabs(value) < Tiny)  Return

    Call TwoBodyLoopComplete(irow, y, y, 0.5D0*value, l, l, l, l)

    Return
  End Subroutine AppendWW



  Subroutine AppendUpperW(irow, y, value, l)
    Implicit None
    Integer          :: irow, y
    Double Precision :: value
    Integer          :: l
    ! End of dummy parameters.

    If (Dabs(value) < Tiny)  Return

    Call OneBodyLoopComplete(irow, y, y, value, l, l)

    If (DRT(irow)%k == 0 .Or. oneOnly)  Return

    Call WalkDown(irow, y, y, value, l, l, .true.)  ! dummy call

    Return
  End Subroutine AppendUpperW


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the upper joining segments  ------
!------------------------------------------------------------------------------

  Recursive Subroutine AppendUpperR(ibra, iket, ybra, yket, value, l)
    Implicit None
    Integer          :: ibra, iket, ybra, yket
    Double Precision :: value
    Integer          :: l
    ! End of dummy parameters.
    Double Precision :: uu
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab

    If (Dabs(value) < Tiny)  Return

    jbra0  =     DRT(ibra)%jd(0)
    jbra1  =     DRT(ibra)%jd(1)
    jbra2  =     DRT(ibra)%jd(2)
    jbra3  =     DRT(ibra)%jd(3)

    jket0  =     DRT(iket)%jd(0)
    jket1  =     DRT(iket)%jd(1)
    jket2  =     DRT(iket)%jd(2)
    jket3  =     DRT(iket)%jd(3)

    ybra0  =     DRT(ibra)%yd(0)
    ybra1  =     DRT(ibra)%yd(1)
    ybra2  =     DRT(ibra)%yd(2)
    ybra3  =     DRT(ibra)%yd(3)

    yket0  =     DRT(iket)%yd(0)
    yket1  =     DRT(iket)%yd(1)
    yket2  =     DRT(iket)%yd(2)
    yket3  =     DRT(iket)%yd(3)

    h      =     DRT(iket)%k
    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b

    uu     =     Occ(h)

    If (deltab == -1) Then

       If (jbra1 == jket0 .And. jket0 /= 0) Then
          Call AppendUpperRb(jket0, ybra+ybra1, yket+yket0,                              value, h, l)
          Call AppendBottom (jket0, ybra+ybra1, yket+yket0,                    -0.5D0*uu*value, h, h, h, l)  ! AppendRbW
       End If

       If (jbra3 == jket2 .And. jket2 /= 0) Then
          Call AppendUpperRb(jket2, ybra+ybra3, yket+yket2,                 FuncA(b,1,2)*value, h, l)
          Call AppendBottom (jket2, ybra+ybra3, yket+yket2, (1.D0-0.5D0*uu)*FuncA(b,1,2)*value, h, h, h, l)  ! AppendRbW
       End If

       If (h <= 1)  Return

       If (jbra0 /= 0 .And. jket0 /= 0) Then
          Call AppendUpperR(jbra0, jket0, ybra+ybra0, yket+yket0,                           value,           l)
          Call AppendLowerR(jbra0, jket0, ybra+ybra0, yket+yket0,                       -uu*value,           l, h, h)  ! AppendRW
          Call AppendLowerR(jbra0, jket0, ybra+ybra0, yket+yket0,                  0.5D0*uu*value,           h, h, l)  ! AppendRbRt
       End If

       If (jbra1 /= 0 .And. jket1 /= 0) Then
          Call AppendUpperR(jbra1, jket1, ybra+ybra1, yket+yket1,                          -value,           l)
          Call AppendLowerR(jbra1, jket1, ybra+ybra1, yket+yket1,                 (uu-1.D0)*value,           l, h, h)  ! AppendRW
          Call AppendLowerR(jbra1, jket1, ybra+ybra1, yket+yket1,           (1.D0-0.5D0*uu)*value,           h, h, l)  ! AppendRbRt
       End If

       If (jbra1 /= 0 .And. jket2 /= 0) Then
          Call AppendUpperR(jbra1, jket2, ybra+ybra1, yket+yket2,                          -value/Dble(b+2), l)
          Call AppendLowerR(jbra1, jket2, ybra+ybra1, yket+yket2,                 (uu-1.D0)*value/Dble(b+2), l, h, h)  ! AppendRW
          Call AppendLowerR(jbra1, jket2, ybra+ybra1, yket+yket2, (1.D0-0.5D0*uu/Dble(b+2))*value,           h, h, l)  ! AppendRbRt
       End If

       If (jbra2 /= 0 .And. jket2 /= 0) Then
          Call AppendUpperR(jbra2, jket2, ybra+ybra2, yket+yket2,                FuncC(b,2)*value,           l)
          Call AppendLowerR(jbra2, jket2, ybra+ybra2, yket+yket2,      (1.D0-uu)*FuncC(b,2)*value,           l, h, h)  ! AppendRW
          Call AppendLowerR(jbra2, jket2, ybra+ybra2, yket+yket2,       0.5D0*uu*FuncC(b,2)*value,           h, h, l)  ! AppendRbRt
       End If

       If (jbra3 /= 0 .And. jket3 /= 0) Then
          Call AppendUpperR(jbra3, jket3, ybra+ybra3, yket+yket3,                          -value,           l)
          Call AppendLowerR(jbra3, jket3, ybra+ybra3, yket+yket3,                 (uu-2.D0)*value,           l, h, h)  ! AppendRW
          Call AppendLowerR(jbra3, jket3, ybra+ybra3, yket+yket3,           (1.D0-0.5D0*uu)*value,           h, h, l)  ! AppendRbRt
       End If

       If (oneOnly)  Return

       If (jbra3 /= 0 .And. jket0 /= 0) &
     & Call AppendL(jbra3, jket0, ybra+ybra3, yket+yket0, FuncA(b,1,2)*value, h, h, l)  ! AppendRbLt

       If (jbra1 /= 0 .And. jket0 /= 0) &
     & Call AppendRL(jbra1, jket0, ybra+ybra1, yket+yket0, value,  t,              -t*FuncA(b,0,2), h, l, .true.)  ! AppendRLt
       If (jbra2 /= 0 .And. jket0 /= 0) &
     & Call AppendRL(jbra2, jket0, ybra+ybra2, yket+yket0, value,    0.D0,            1.D0,         h, l, .true.)  ! AppendRLt
       If (jbra3 /= 0 .And. jket1 /= 0) &
     & Call AppendRL(jbra3, jket1, ybra+ybra3, yket+yket1, value,    0.D0,           -FuncA(b,1,2), h, l, .true.)  ! AppendRLt
       If (jbra3 /= 0 .And. jket2 /= 0) &
     & Call AppendRL(jbra3, jket2, ybra+ybra3, yket+yket2, value,  t*FuncA(b,1,2), -t*FuncA(b,3,2), h, l, .true.)  ! AppendRLt

       If (jbra0 /= 0 .And. jket1 /= 0) &
     & Call AppendRR(jbra0, jket1, ybra+ybra0, yket+yket1, value,    0.D0,            1.D0,         h, l, .false.)  ! AppendRRt
       If (jbra0 /= 0 .And. jket2 /= 0) &
     & Call AppendRR(jbra0, jket2, ybra+ybra0, yket+yket2, value,  t*FuncA(b,1,2),  t*FuncA(b,3,2), h, l, .false.)  ! AppendRRt
       If (jbra1 /= 0 .And. jket3 /= 0) &
     & Call AppendRR(jbra1, jket3, ybra+ybra1, yket+yket3, value, -t,              -t*FuncA(b,0,2), h, l, .false.)  ! AppendRRt
       If (jbra2 /= 0 .And. jket3 /= 0) &
     & Call AppendRR(jbra2, jket3, ybra+ybra2, yket+yket3, value,    0.D0,            FuncA(b,3,2), h, l, .false.)  ! AppendRRt

    Else If (deltab == 1) Then

       If (jbra2 == jket0 .And. jket0 /= 0) Then
          Call AppendUpperRb(jket0, ybra+ybra2, yket+yket0,                              value, h, l)
          Call AppendBottom (jket0, ybra+ybra2, yket+yket0,                    -0.5D0*uu*value, h, h, h, l)  ! AppendRbW
       End If

       If (jbra3 == jket1 .And. jket1 /= 0) Then
          Call AppendUpperRb(jket1, ybra+ybra3, yket+yket1,                 FuncA(b,1,0)*value, h, l)
          Call AppendBottom (jket1, ybra+ybra3, yket+yket1, (1.D0-0.5D0*uu)*FuncA(b,1,0)*value, h, h, h, l)  ! AppendRbW
       End If

       If (h <= 1)  Return

       If (jbra0 /= 0 .And. jket0 /= 0) Then
          Call AppendUpperR(jbra0, jket0, ybra+ybra0, yket+yket0,                         value,         l)
          Call AppendLowerR(jbra0, jket0, ybra+ybra0, yket+yket0,                     -uu*value,         l, h, h)  ! AppendRW
          Call AppendLowerR(jbra0, jket0, ybra+ybra0, yket+yket0,                0.5D0*uu*value,         h, h, l)  ! AppendRbRt
       End If

       If (jbra1 /= 0 .And. jket1 /= 0) Then
          Call AppendUpperR(jbra1, jket1, ybra+ybra1, yket+yket1,              FuncC(b,0)*value,         l)
          Call AppendLowerR(jbra1, jket1, ybra+ybra1, yket+yket1,    (1.D0-uu)*FuncC(b,0)*value,         l, h, h)  ! AppendRW
          Call AppendLowerR(jbra1, jket1, ybra+ybra1, yket+yket1,     0.5D0*uu*FuncC(b,0)*value,         h, h, l)  ! AppendRbRt
       End If

       If (jbra2 /= 0 .And. jket1 /= 0) Then
          Call AppendUpperR(jbra2, jket1, ybra+ybra2, yket+yket1,                         value/Dble(b), l)
          Call AppendLowerR(jbra2, jket1, ybra+ybra2, yket+yket1,               (1.D0-uu)*value/Dble(b), l, h, h)  ! AppendRW
          Call AppendLowerR(jbra2, jket1, ybra+ybra2, yket+yket1, (1.D0+0.5D0*uu/Dble(b))*value,         h, h, l)  ! AppendRbRt
       End If

       If (jbra2 /= 0 .And. jket2 /= 0) Then
          Call AppendUpperR(jbra2, jket2, ybra+ybra2, yket+yket2,                        -value,         l)
          Call AppendLowerR(jbra2, jket2, ybra+ybra2, yket+yket2,               (uu-1.D0)*value,         l, h, h)  ! AppendRW
          Call AppendLowerR(jbra2, jket2, ybra+ybra2, yket+yket2,         (1.D0-0.5D0*uu)*value,         h, h, l)  ! AppendRbRt
       End If

       If (jbra3 /= 0 .And. jket3 /= 0) Then
          Call AppendUpperR(jbra3, jket3, ybra+ybra3, yket+yket3,                        -value,         l)
          Call AppendLowerR(jbra3, jket3, ybra+ybra3, yket+yket3,               (uu-2.D0)*value,         l, h, h)  ! AppendRW
          Call AppendLowerR(jbra3, jket3, ybra+ybra3, yket+yket3,         (1.D0-0.5D0*uu)*value,         h, h, l)  ! AppendRbRt
       End If

       If (oneOnly)  Return

       If (jbra3 /= 0 .And. jket0 /= 0) &
     & Call AppendL(jbra3, jket0, ybra+ybra3, yket+yket0, FuncA(b,1,0)*value, h, h, l)  ! AppendRbLt

       If (jbra1 /= 0 .And. jket0 /= 0) &
     & Call AppendRL(jbra1, jket0, ybra+ybra1, yket+yket0, value,    0.D0,            1.D0,          h, l, .true.)  ! AppendRLt
       If (jbra2 /= 0 .And. jket0 /= 0) &
     & Call AppendRL(jbra2, jket0, ybra+ybra2, yket+yket0, value,  t,               t*FuncA(b, 2,0), h, l, .true.)  ! AppendRLt
       If (jbra3 /= 0 .And. jket1 /= 0) &
     & Call AppendRL(jbra3, jket1, ybra+ybra3, yket+yket1, value,  t*FuncA(b,1,0),  t*FuncA(b,-1,0), h, l, .true.)  ! AppendRLt
       If (jbra3 /= 0 .And. jket2 /= 0) &
     & Call AppendRL(jbra3, jket2, ybra+ybra3, yket+yket2, value,    0.D0,           -FuncA(b, 1,0), h, l, .true.)  ! AppendRLt

       If (jbra0 /= 0 .And. jket1 /= 0) &
     & Call AppendRR(jbra0, jket1, ybra+ybra0, yket+yket1, value,  t*FuncA(b,1,0), -t*FuncA(b,-1,0), h, l, .false.)  ! AppendRRt
       If (jbra0 /= 0 .And. jket2 /= 0) &
     & Call AppendRR(jbra0, jket2, ybra+ybra0, yket+yket2, value,    0.D0,            1.D0,          h, l, .false.)  ! AppendRRt
       If (jbra1 /= 0 .And. jket3 /= 0) &
     & Call AppendRR(jbra1, jket3, ybra+ybra1, yket+yket3, value,    0.D0,            FuncA(b,-1,0), h, l, .false.)  ! AppendRRt
       If (jbra2 /= 0 .And. jket3 /= 0) &
     & Call AppendRR(jbra2, jket3, ybra+ybra2, yket+yket3, value, -t,               t*FuncA(b, 2,0), h, l, .false.)  ! AppendRRt

    End If

    Return
  End Subroutine AppendUpperR


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the upper middle segments  -------
!------------------------------------------------------------------------------

  Subroutine AppendUpperRb(irow, ybra, yket, value, k, l)
    Implicit None
    Integer          :: irow, ybra, yket
    Double Precision :: value
    Integer          :: k, l
    ! End of dummy parameters.

    If (Dabs(value) < Tiny) Return

    Call OneBodyLoopComplete(irow, ybra, yket, value, k, l)

    If (DRT(irow)%k == 0 .Or. oneOnly)  Return

    Call WalkDown(irow, ybra, yket, value, k, l, .false.)  ! dummy call

    Return
  End Subroutine AppendUpperRb


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the middle joining segments  -----
!------------------------------------------------------------------------------

  Recursive Subroutine WalkDown(irow, ybra, yket, value, k, l, fromW)
    Implicit None
    Integer          :: irow, ybra, yket
    Double Precision :: value
    Integer          :: k, l
    Logical          :: fromW
    ! End of dummy parameters.
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3, y0, y1, y2, y3, h, b

    j0 = DRT(irow)%jd(0)
    j1 = DRT(irow)%jd(1)
    j2 = DRT(irow)%jd(2)
    j3 = DRT(irow)%jd(3)

    y0 = DRT(irow)%yd(0)
    y1 = DRT(irow)%yd(1)
    y2 = DRT(irow)%yd(2)
    y3 = DRT(irow)%yd(3)

    h  = DRT(irow)%k
    b  = DRT(irow)%b

    uu = Occ(h)

    If (j0 /= 0)  Call AppendLowerW(j0, ybra+y0, yket+y0,      -uu *value, h, k, l)
    If (j1 /= 0)  Call AppendLowerW(j1, ybra+y1, yket+y1, (1.D0-uu)*value, h, k, l)
    If (j2 /= 0)  Call AppendLowerW(j2, ybra+y2, yket+y2, (1.D0-uu)*value, h, k, l)
    If (j3 /= 0)  Call AppendLowerW(j3, ybra+y3, yket+y3, (2.D0-uu)*value, h, k, l)

    If (h <= 1)  Return

    If (j0 /= 0)  Call WalkDown(j0, ybra+y0, yket+y0, value, k, l, fromW)
    If (j1 /= 0)  Call WalkDown(j1, ybra+y1, yket+y1, value, k, l, fromW)
    If (j2 /= 0)  Call WalkDown(j2, ybra+y2, yket+y2, value, k, l, fromW)
    If (j3 /= 0)  Call WalkDown(j3, ybra+y3, yket+y3, value, k, l, fromW)

    If (diagOnly)  Return

    If (j0 /= 0 .And. j1 /= 0)  Call AppendLowerR(j0, j1, ybra+y0, yket+y1,              value, h, k, l)  ! AppendLowerRt
    If (j0 /= 0 .And. j2 /= 0)  Call AppendLowerR(j0, j2, ybra+y0, yket+y2,              value, h, k, l)  ! AppendLowerRt
    If (j1 /= 0 .And. j3 /= 0)  Call AppendLowerR(j1, j3, ybra+y1, yket+y3, FuncA(b,0,1)*value, h, k, l)  ! AppendLowerRt
    If (j2 /= 0 .And. j3 /= 0)  Call AppendLowerR(j2, j3, ybra+y2, yket+y3, FuncA(b,2,1)*value, h, k, l)  ! AppendLowerRt

    If (fromW)  Return

    If (j1 /= 0 .And. j0 /= 0)  Call AppendL(j1, j0, ybra+y1, yket+y0,              value, h, k, l)  ! AppendLt
    If (j2 /= 0 .And. j0 /= 0)  Call AppendL(j2, j0, ybra+y2, yket+y0,              value, h, k, l)  ! AppendLt
    If (j3 /= 0 .And. j1 /= 0)  Call AppendL(j3, j1, ybra+y3, yket+y1, FuncA(b,0,1)*value, h, k, l)  ! AppendLt
    If (j3 /= 0 .And. j2 /= 0)  Call AppendL(j3, j2, ybra+y3, yket+y2, FuncA(b,2,1)*value, h, k, l)  ! AppendLt

    Return
  End Subroutine WalkDown



  Recursive Subroutine AppendRL(ibra, iket, ybra, yket, value, sval, tval, k, l, hasBeenOpen)
    Implicit None
    Integer          :: ibra, iket, ybra, yket
    Double Precision :: value, sval, tval
    Integer          :: k, l
    Logical          :: hasBeenOpen
    ! End of dummy parameters.
    Double Precision :: uu
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab

    If (Dabs(sval) < Tiny .And. Dabs(tval) < Tiny)  Return

    jbra0  =     DRT(ibra)%jd(0)
    jbra1  =     DRT(ibra)%jd(1)
    jbra2  =     DRT(ibra)%jd(2)
    jbra3  =     DRT(ibra)%jd(3)

    jket0  =     DRT(iket)%jd(0)
    jket1  =     DRT(iket)%jd(1)
    jket2  =     DRT(iket)%jd(2)
    jket3  =     DRT(iket)%jd(3)

    ybra0  =     DRT(ibra)%yd(0)
    ybra1  =     DRT(ibra)%yd(1)
    ybra2  =     DRT(ibra)%yd(2)
    ybra3  =     DRT(ibra)%yd(3)

    yket0  =     DRT(iket)%yd(0)
    yket1  =     DRT(iket)%yd(1)
    yket2  =     DRT(iket)%yd(2)
    yket3  =     DRT(iket)%yd(3)

    h      =     DRT(iket)%k
    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b

    uu     =     Occ(h)

    If (deltab == -2) Then

       If (jbra1 == jket2 .And. jket2 /= 0) &
     & Call AppendBottom(jket2, ybra+ybra1, yket+yket2, FuncC(b,2)*tval*value, h, k, h, l)  ! AppendRbLb

       If (h <= 1)  Return

       If (jbra0 /= 0 .And. jket0 /= 0) &
     & Call AppendRL(jbra0, jket0, ybra+ybra0, yket+yket0, value, 0.D0,             tval,           k, l, .true.)
       If (jbra1 /= 0 .And. jket1 /= 0) &
     & Call AppendRL(jbra1, jket1, ybra+ybra1, yket+yket1, value, 0.D0, -FuncC(b,2)*tval,           k, l, .true.)
       If (jbra1 /= 0 .And. jket2 /= 0 .And. .Not. diagOnly) &
     & Call AppendRL(jbra1, jket2, ybra+ybra1, yket+yket2, value, 0.D0,      -Sqrt2*tval/Dble(b+2), k, l, .true.)
       If (jbra2 /= 0 .And. jket2 /= 0) &
     & Call AppendRL(jbra2, jket2, ybra+ybra2, yket+yket2, value, 0.D0, -FuncC(b,2)*tval,           k, l, .true.)
       If (jbra3 /= 0 .And. jket3 /= 0) &
     & Call AppendRL(jbra3, jket3, ybra+ybra3, yket+yket3, value, 0.D0,             tval,           k, l, .true.)

       If (diagOnly)  Return

       If (jbra0 /= 0 .And. jket2 /= 0) &
     & Call AppendLowerR(jbra0, jket2, ybra+ybra0, yket+yket2,    FuncC(b,2)*tval*value, l, h, k)  ! AppendRLb
       If (jbra1 /= 0 .And. jket3 /= 0) &
     & Call AppendLowerR(jbra1, jket3, ybra+ybra1, yket+yket3, -FuncA(b,3,2)*tval*value, l, h, k)  ! AppendRLb

       If (jbra1 /= 0 .And. jket0 /= 0) &
     & Call AppendL(jbra1, jket0, ybra+ybra1, yket+yket0,    FuncC(b,2)*tval*value, k, h, l)  ! AppendRbL
       If (jbra3 /= 0 .And. jket2 /= 0) &
     & Call AppendL(jbra3, jket2, ybra+ybra3, yket+yket2, -FuncA(b,1,2)*tval*value, k, h, l)  ! AppendRbL

    Else If (deltab == 0) Then

       If (jbra0 == jket0 .And. jket0 /= 0) &
     & Call AppendBottom(jket0, ybra+ybra0, yket+yket0, t*value*( uu      *sval),                   h, k, h, l)  ! AppendRbLb
       If (jbra1 == jket1 .And. jket1 /= 0) &
     & Call AppendBottom(jket1, ybra+ybra1, yket+yket1, t*value*((uu-1.D0)*sval+FuncA(b,2,0)*tval), h, k, h, l)  ! AppendRbLb
       If (jbra2 == jket2 .And. jket2 /= 0) &
     & Call AppendBottom(jket2, ybra+ybra2, yket+yket2, t*value*((uu-1.D0)*sval-FuncA(b,0,2)*tval), h, k, h, l)  ! AppendRbLb
       If (jbra3 == jket3 .And. jket3 /= 0) &
     & Call AppendBottom(jket3, ybra+ybra3, yket+yket3, t*value*((uu-2.D0)*sval),                   h, k, h, l)  ! AppendRbLb

       If (h <= 1)  Return

       If (jbra0 /= 0 .And. jket0 /= 0) &
     & Call AppendRL(jbra0, jket0, ybra+ybra0, yket+yket0, value, sval,               tval, k, l, hasBeenOpen)
       If (jbra1 /= 0 .And. jket1 /= 0) &
     & Call AppendRL(jbra1, jket1, ybra+ybra1, yket+yket1, value, sval,    FuncD(b,0)*tval, k, l, hasBeenOpen)
       If (jbra1 /= 0 .And. jket2 /= 0 .And. .Not. diagOnly) &
     & Call AppendRL(jbra1, jket2, ybra+ybra1, yket+yket2, value, 0.D0, -FuncB(b,0,2)*tval, k, l, .true.)

       If (jbra2 /= 0 .And. jket1 /= 0 .And. hasBeenOpen) &
     & Call AppendRL(jbra2, jket1, ybra+ybra2, yket+yket1, value, 0.D0, -FuncB(b,0,2)*tval, k, l, .true.)
       If (jbra2 /= 0 .And. jket2 /= 0) &
     & Call AppendRL(jbra2, jket2, ybra+ybra2, yket+yket2, value, sval,    FuncD(b,1)*tval, k, l, hasBeenOpen)
       If (jbra3 /= 0 .And. jket3 /= 0) &
     & Call AppendRL(jbra3, jket3, ybra+ybra3, yket+yket3, value, sval,               tval, k, l, hasBeenOpen)

       If (diagOnly)  Return

       If (jbra0 /= 0 .And. jket1 /= 0) &
     & Call AppendLowerR(jbra0, jket1, ybra+ybra0, yket+yket1, -t*value*(             sval-FuncA(b,2,0)*tval), l, h, k)  ! AppendRLb
       If (jbra0 /= 0 .And. jket2 /= 0) &
     & Call AppendLowerR(jbra0, jket2, ybra+ybra0, yket+yket2, -t*value*(             sval+FuncA(b,0,2)*tval), l, h, k)  ! AppendRLb
       If (jbra1 /= 0 .And. jket3 /= 0) &
     & Call AppendLowerR(jbra1, jket3, ybra+ybra1, yket+yket3, -t*value*(FuncA(b,0,1)*sval+FuncA(b,2,1)*tval), l, h, k)  ! AppendRLb
       If (jbra2 /= 0 .And. jket3 /= 0) &
     & Call AppendLowerR(jbra2, jket3, ybra+ybra2, yket+yket3, -t*value*(FuncA(b,2,1)*sval-FuncA(b,0,1)*tval), l, h, k)  ! AppendRLb

       If (hasBeenOpen) Then
          If (jbra1 /= 0 .And. jket0 /= 0) &
        & Call AppendL(jbra1, jket0, ybra+ybra1, yket+yket0, -t*value*(             sval-FuncA(b,2,0)*tval), k, h, l)  ! AppendRbL
          If (jbra2 /= 0 .And. jket0 /= 0) &
        & Call AppendL(jbra2, jket0, ybra+ybra2, yket+yket0, -t*value*(             sval+FuncA(b,0,2)*tval), k, h, l)  ! AppendRbL
          If (jbra3 /= 0 .And. jket1 /= 0) &
        & Call AppendL(jbra3, jket1, ybra+ybra3, yket+yket1, -t*value*(FuncA(b,0,1)*sval+FuncA(b,2,1)*tval), k, h, l)  ! AppendRbL
          If (jbra3 /= 0 .And. jket2 /= 0) &
        & Call AppendL(jbra3, jket2, ybra+ybra3, yket+yket2, -t*value*(FuncA(b,2,1)*sval-FuncA(b,0,1)*tval), k, h, l)  ! AppendRbL
       End If

    Else If (deltab == 2) Then

       If (jbra2 == jket1 .And. jket1 /= 0) &
     & Call AppendBottom(jket1, ybra+ybra2, yket+yket1, FuncC(b,0)*tval*value, h, k, h, l)  ! AppendRbLb

       If (h <= 1)  Return

       If (jbra0 /= 0 .And. jket0 /= 0) &
     & Call AppendRL(jbra0, jket0, ybra+ybra0, yket+yket0, value, 0.D0,             tval,         k, l, .true.)
       If (jbra1 /= 0 .And. jket1 /= 0) &
     & Call AppendRL(jbra1, jket1, ybra+ybra1, yket+yket1, value, 0.D0, -FuncC(b,0)*tval,         k, l, .true.)
       If (jbra2 /= 0 .And. jket1 /= 0 .And. .Not. diagOnly) &
     & Call AppendRL(jbra2, jket1, ybra+ybra2, yket+yket1, value, 0.D0,      -Sqrt2*tval/Dble(b), k, l, .true.)
       If (jbra2 /= 0 .And. jket2 /= 0) &
     & Call AppendRL(jbra2, jket2, ybra+ybra2, yket+yket2, value, 0.D0, -FuncC(b,0)*tval,         k, l, .true.)
       If (jbra3 /= 0 .And. jket3 /= 0) &
     & Call AppendRL(jbra3, jket3, ybra+ybra3, yket+yket3, value, 0.D0,             tval,         k, l, .true.)

       If (diagOnly)  Return

       If (jbra0 /= 0 .And. jket1 /= 0) &
     & Call AppendLowerR(jbra0, jket1, ybra+ybra0, yket+yket1,     FuncC(b,0)*tval*value, l, h, k)  ! AppendRLb
       If (jbra2 /= 0 .And. jket3 /= 0) &
     & Call AppendLowerR(jbra2, jket3, ybra+ybra2, yket+yket3, -FuncA(b,-1,0)*tval*value, l, h, k)  ! AppendRLb

       If (jbra2 /= 0 .And. jket0 /= 0) &
     & Call AppendL(jbra2, jket0, ybra+ybra2, yket+yket0,    FuncC(b,0)*tval*value, k, h, l)  ! AppendRbL
       If (jbra3 /= 0 .And. jket1 /= 0) &
     & Call AppendL(jbra3, jket1, ybra+ybra3, yket+yket1, -FuncA(b,1,0)*tval*value, k, h, l)  ! AppendRbL

    End If

    Return
  End Subroutine AppendRL



  Recursive Subroutine AppendRR(ibra, iket, ybra, yket, value, sval, tval, k, l, fromRtRt)
    Implicit None
    Integer          :: ibra, iket, ybra, yket
    Double Precision :: value, sval, tval
    Integer          :: k, l
    Logical          :: fromRtRt
    ! End of dummy parameters.
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab
    Double Precision :: f

    If (Dabs(sval) < Tiny .And. Dabs(tval) < Tiny)  Return

    jbra0  =     DRT(ibra)%jd(0)
    jbra1  =     DRT(ibra)%jd(1)
    jbra2  =     DRT(ibra)%jd(2)
    jbra3  =     DRT(ibra)%jd(3)

    jket0  =     DRT(iket)%jd(0)
    jket1  =     DRT(iket)%jd(1)
    jket2  =     DRT(iket)%jd(2)
    jket3  =     DRT(iket)%jd(3)

    ybra0  =     DRT(ibra)%yd(0)
    ybra1  =     DRT(ibra)%yd(1)
    ybra2  =     DRT(ibra)%yd(2)
    ybra3  =     DRT(ibra)%yd(3)

    yket0  =     DRT(iket)%yd(0)
    yket1  =     DRT(iket)%yd(1)
    yket2  =     DRT(iket)%yd(2)
    yket3  =     DRT(iket)%yd(3)

    h      =     DRT(iket)%k
    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b

    If (deltab == -2) Then

       If (jbra0 /= 0 .And. jket0 /= 0) &
     & Call AppendRR(jbra0, jket0, ybra+ybra0, yket+yket0, value, 0.D0,              tval, k, l, fromRtRt)
       If (jbra1 /= 0 .And. jket1 /= 0) &
     & Call AppendRR(jbra1, jket1, ybra+ybra1, yket+yket1, value, 0.D0,              tval, k, l, fromRtRt)
       If (jbra1 /= 0 .And. jket2 /= 0) &
     & Call AppendRR(jbra1, jket2, ybra+ybra1, yket+yket2, value, 0.D0, FuncB(b,2,3)*tval, k, l, fromRtRt)
       If (jbra2 /= 0 .And. jket2 /= 0) &
     & Call AppendRR(jbra2, jket2, ybra+ybra2, yket+yket2, value, 0.D0,   FuncD(b,2)*tval, k, l, fromRtRt)
       If (jbra3 /= 0 .And. jket3 /= 0) &
     & Call AppendRR(jbra3, jket3, ybra+ybra3, yket+yket3, value, 0.D0,              tval, k, l, fromRtRt)

       If (jbra1 /= 0 .And. jket0 /= 0) &
     & Call AppendLowerR(jbra1, jket0, ybra+ybra1, yket+yket0,               tval*value, k, h, l)  ! AppendRbR
       If (jbra3 /= 0 .And. jket2 /= 0) &
     & Call AppendLowerR(jbra3, jket2, ybra+ybra3, yket+yket2,  FuncA(b,1,2)*tval*value, k, h, l)  ! AppendRbR

       If (fromRtRt) Return

       If (jbra1 /= 0 .And. jket0 /= 0) &
     & Call AppendLowerR(jbra1, jket0, ybra+ybra1, yket+yket0,              -tval*value, l, h, k)  ! AppendRRb
       If (jbra3 /= 0 .And. jket2 /= 0) &
     & Call AppendLowerR(jbra3, jket2, ybra+ybra3, yket+yket2, -FuncA(b,1,2)*tval*value, l, h, k)  ! AppendRRb

    Else If (deltab == 0) Then

       If (jbra3 == jket0 .And. jket0 /= 0) Then
          f = 1.0D0
          If (fromRtRt) f = 0.5D0
          Call AppendBottom(jket0, ybra+ybra3, yket+yket0, f*Sqrt2*sval*value, h, k, h, l)  ! AppendRbRb
       End If

       If (jbra0 /= 0 .And. jket0 /= 0) &
     & Call AppendRR(jbra0, jket0, ybra+ybra0, yket+yket0, value,  sval,              tval, k, l, fromRtRt)
       If (jbra1 /= 0 .And. jket1 /= 0) &
     & Call AppendRR(jbra1, jket1, ybra+ybra1, yket+yket1, value, -sval,  -FuncD(b,0)*tval, k, l, fromRtRt)
       If (jbra1 /= 0 .And. jket2 /= 0) &
     & Call AppendRR(jbra1, jket2, ybra+ybra1, yket+yket2, value,  0.D0, FuncB(b,1,2)*tval, k, l, fromRtRt)
       If (jbra2 /= 0 .And. jket1 /= 0) &
     & Call AppendRR(jbra2, jket1, ybra+ybra2, yket+yket1, value,  0.D0, FuncB(b,0,1)*tval, k, l, fromRtRt)
       If (jbra2 /= 0 .And. jket2 /= 0) &
     & Call AppendRR(jbra2, jket2, ybra+ybra2, yket+yket2, value, -sval,  -FuncD(b,1)*tval, k, l, fromRtRt)
       If (jbra3 /= 0 .And. jket3 /= 0) &
     & Call AppendRR(jbra3, jket3, ybra+ybra3, yket+yket3, value,  sval,              tval, k, l, fromRtRt)

       If (jbra1 /= 0 .And. jket0 /= 0) &
     & Call AppendLowerR(jbra1, jket0, ybra+ybra1, yket+yket0, t*value*(FuncA(b,0,1)*sval+FuncA(b,2,1)*tval), k, h, l)  ! AppendRbR
       If (jbra2 /= 0 .And. jket0 /= 0) &
     & Call AppendLowerR(jbra2, jket0, ybra+ybra2, yket+yket0, t*value*(FuncA(b,2,1)*sval-FuncA(b,0,1)*tval), k, h, l)  ! AppendRbR
       If (jbra3 /= 0 .And. jket1 /= 0) &
     & Call AppendLowerR(jbra3, jket1, ybra+ybra3, yket+yket1, t*value*(            -sval+FuncA(b,2,0)*tval), k, h, l)  ! AppendRbR
       If (jbra3 /= 0 .And. jket2 /= 0) &
     & Call AppendLowerR(jbra3, jket2, ybra+ybra3, yket+yket2, t*value*(            -sval-FuncA(b,0,2)*tval), k, h, l)  ! AppendRbR

       If (fromRtRt) Return

       If (jbra1 /= 0 .And. jket0 /= 0) &
     & Call AppendLowerR(jbra1, jket0, ybra+ybra1, yket+yket0, t*value*(FuncA(b,0,1)*sval-FuncA(b,2,1)*tval), l, h, k)  ! AppendRRb
       If (jbra2 /= 0 .And. jket0 /= 0) &
     & Call AppendLowerR(jbra2, jket0, ybra+ybra2, yket+yket0, t*value*(FuncA(b,2,1)*sval+FuncA(b,0,1)*tval), l, h, k)  ! AppendRRb
       If (jbra3 /= 0 .And. jket1 /= 0) &
     & Call AppendLowerR(jbra3, jket1, ybra+ybra3, yket+yket1, t*value*(            -sval-FuncA(b,2,0)*tval), l, h, k)  ! AppendRRb
       If (jbra3 /= 0 .And. jket2 /= 0) &
     & Call AppendLowerR(jbra3, jket2, ybra+ybra3, yket+yket2, t*value*(            -sval+FuncA(b,0,2)*tval), l, h, k)  ! AppendRRb

    Else If (deltab == 2) Then

       If (jbra0 /= 0 .And. jket0 /= 0) &
     & Call AppendRR(jbra0, jket0, ybra+ybra0, yket+yket0, value, 0.D0,               tval, k, l, fromRtRt)
       If (jbra1 /= 0 .And. jket1 /= 0) &
     & Call AppendRR(jbra1, jket1, ybra+ybra1, yket+yket1, value, 0.D0,   FuncD(b,-1)*tval, k, l, fromRtRt)
       If (jbra2 /= 0 .And. jket1 /= 0) &
     & Call AppendRR(jbra2, jket1, ybra+ybra2, yket+yket1, value, 0.D0, FuncB(b,-1,0)*tval, k, l, fromRtRt)
       If (jbra2 /= 0 .And. jket2 /= 0) &
     & Call AppendRR(jbra2, jket2, ybra+ybra2, yket+yket2, value, 0.D0,               tval, k, l, fromRtRt)
       If (jbra3 /= 0 .And. jket3 /= 0) &
     & Call AppendRR(jbra3, jket3, ybra+ybra3, yket+yket3, value, 0.D0,               tval, k, l, fromRtRt)

       If (jbra2 /= 0 .And. jket0 /= 0) &
     & Call AppendLowerR(jbra2, jket0, ybra+ybra2, yket+yket0,               tval*value, k, h, l)  ! AppendRbR
       If (jbra3 /= 0 .And. jket1 /= 0) &
     & Call AppendLowerR(jbra3, jket1, ybra+ybra3, yket+yket1,  FuncA(b,1,0)*tval*value, k, h, l)  ! AppendRbR

       If (fromRtRt) Return

       If (jbra2 /= 0 .And. jket0 /= 0) &
     & Call AppendLowerR(jbra2, jket0, ybra+ybra2, yket+yket0,              -tval*value, l, h, k)  ! AppendRRb
       If (jbra3 /= 0 .And. jket1 /= 0) &
     & Call AppendLowerR(jbra3, jket1, ybra+ybra3, yket+yket1, -FuncA(b,1,0)*tval*value, l, h, k)  ! AppendRRb
    End If

    Return
  End Subroutine AppendRR


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the lower middle segments  -------
!------------------------------------------------------------------------------

  ! These subroutines would contain dummy calls only and have been removed.


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the lower joining segments  ------
!------------------------------------------------------------------------------

  Recursive Subroutine AppendLowerR(ibra, iket, ybra, yket, value, j, k, l)
    Implicit None
    Integer          :: ibra, iket, ybra, yket
    Double Precision :: value
    Integer          :: j, k, l
    ! End of dummy parameters.
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab

    If (Dabs(value) < Tiny .Or. oneOnly) Return

    jbra0  =     DRT(ibra)%jd(0)
    jbra1  =     DRT(ibra)%jd(1)
    jbra2  =     DRT(ibra)%jd(2)
    jbra3  =     DRT(ibra)%jd(3)

    jket0  =     DRT(iket)%jd(0)
    jket1  =     DRT(iket)%jd(1)
    jket2  =     DRT(iket)%jd(2)
    jket3  =     DRT(iket)%jd(3)

    ybra0  =     DRT(ibra)%yd(0)
    ybra1  =     DRT(ibra)%yd(1)
    ybra2  =     DRT(ibra)%yd(2)
    ybra3  =     DRT(ibra)%yd(3)

    yket0  =     DRT(iket)%yd(0)
    yket1  =     DRT(iket)%yd(1)
    yket2  =     DRT(iket)%yd(2)
    yket3  =     DRT(iket)%yd(3)

    h      =     DRT(iket)%k
    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b

    If (deltab == -1) Then

       If (jbra1 == jket0 .And. jket0 /= 0) &
     & Call AppendBottom(jket0, ybra+ybra1, yket+yket0,              value, h, j, k, l)  ! AppendLowerRb
       If (jbra3 == jket2 .And. jket2 /= 0) &
     & Call AppendBottom(jket2, ybra+ybra3, yket+yket2, FuncA(b,1,2)*value, h, j, k, l)  ! AppendLowerRb

       If (h <= 1)  Return

       If (jbra0 /= 0 .And. jket0 /= 0) &
     & Call AppendLowerR(jbra0, jket0, ybra+ybra0, yket+yket0,            value,           j, k, l)
       If (jbra1 /= 0 .And. jket1 /= 0) &
     & Call AppendLowerR(jbra1, jket1, ybra+ybra1, yket+yket1,           -value,           j, k, l)
       If (jbra1 /= 0 .And. jket2 /= 0) &
     & Call AppendLowerR(jbra1, jket2, ybra+ybra1, yket+yket2,           -value/Dble(b+2), j, k, l)
       If (jbra2 /= 0 .And. jket2 /= 0) &
     & Call AppendLowerR(jbra2, jket2, ybra+ybra2, yket+yket2, FuncC(b,2)*value,           j, k, l)
       If (jbra3 /= 0 .And. jket3 /= 0) &
     & Call AppendLowerR(jbra3, jket3, ybra+ybra3, yket+yket3,           -value,           j, k, l)

    Else If (deltab == 1) Then

       If (jbra2 == jket0 .And. jket0 /= 0) &
     & Call AppendBottom(jket0, ybra+ybra2, yket+yket0,              value, h, j, k, l)  ! AppendLowerRb
       If (jbra3 == jket1 .And. jket1 /= 0) &
     & Call AppendBottom(jket1, ybra+ybra3, yket+yket1, FuncA(b,1,0)*value, h, j, k, l)  ! AppendLowerRb

       If (h <= 1)  Return

       If (jbra0 /= 0 .And. jket0 /= 0) &
     & Call AppendLowerR(jbra0, jket0, ybra+ybra0, yket+yket0,            value,         j, k, l)
       If (jbra1 /= 0 .And. jket1 /= 0) &
     & Call AppendLowerR(jbra1, jket1, ybra+ybra1, yket+yket1, FuncC(b,0)*value,         j, k, l)
       If (jbra2 /= 0 .And. jket1 /= 0) &
     & Call AppendLowerR(jbra2, jket1, ybra+ybra2, yket+yket1,            value/Dble(b), j, k, l)
       If (jbra2 /= 0 .And. jket2 /= 0) &
     & Call AppendLowerR(jbra2, jket2, ybra+ybra2, yket+yket2,           -value,         j, k, l)
       If (jbra3 /= 0 .And. jket3 /= 0) &
     & Call AppendLowerR(jbra3, jket3, ybra+ybra3, yket+yket3,           -value,         j, k, l)

    End If

    Return
  End Subroutine AppendLowerR



  Recursive Subroutine AppendL(ibra, iket, ybra, yket, value, j, k, l)
    Implicit None
    Integer          :: ibra, iket, ybra, yket
    Double Precision :: value
    Integer          :: j, k, l
    ! End of dummy parameters.
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab

    If (Dabs(value) < Tiny)  Return

    jbra0  =     DRT(ibra)%jd(0)
    jbra1  =     DRT(ibra)%jd(1)
    jbra2  =     DRT(ibra)%jd(2)
    jbra3  =     DRT(ibra)%jd(3)

    jket0  =     DRT(iket)%jd(0)
    jket1  =     DRT(iket)%jd(1)
    jket2  =     DRT(iket)%jd(2)
    jket3  =     DRT(iket)%jd(3)

    ybra0  =     DRT(ibra)%yd(0)
    ybra1  =     DRT(ibra)%yd(1)
    ybra2  =     DRT(ibra)%yd(2)
    ybra3  =     DRT(ibra)%yd(3)

    yket0  =     DRT(iket)%yd(0)
    yket1  =     DRT(iket)%yd(1)
    yket2  =     DRT(iket)%yd(2)
    yket3  =     DRT(iket)%yd(3)

    h      =     DRT(iket)%k
    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b

    If (deltab == -1) Then
       If (jbra0 == jket2 .And. jket2 /= 0) &
     & Call AppendBottom(jket2, ybra+ybra0, yket+yket2,              value, h, j, k, l)  ! AppendLb
       If (jbra1 == jket3 .And. jket3 /= 0) &
     & Call AppendBottom(jket3, ybra+ybra1, yket+yket3, FuncA(b,2,1)*value, h, j, k, l)  ! AppendLb

       If (h <= 1)  Return

       If (jbra0 /= 0 .And. jket0 /= 0) &
     & Call AppendL(jbra0, jket0, ybra+ybra0, yket+yket0,            value,           j, k, l)
       If (jbra1 /= 0 .And. jket1 /= 0) &
     & Call AppendL(jbra1, jket1, ybra+ybra1, yket+yket1, FuncC(b,1)*value,           j, k, l)
       If (jbra1 /= 0 .And. jket2 /= 0) &
     & Call AppendL(jbra1, jket2, ybra+ybra1, yket+yket2,            value/Dble(b+1), j, k, l)
       If (jbra2 /= 0 .And. jket2 /= 0) &
     & Call AppendL(jbra2, jket2, ybra+ybra2, yket+yket2,           -value,           j, k, l)
       If (jbra3 /= 0 .And. jket3 /= 0) &
     & Call AppendL(jbra3, jket3, ybra+ybra3, yket+yket3,           -value,           j, k, l)

    Else If (deltab == 1) Then
       If (jbra0 == jket1 .And. jket1 /= 0) &
     & Call AppendBottom(jket1, ybra+ybra0, yket+yket1,              value, h, j, k, l)  ! AppendLb
       If (jbra2 == jket3 .And. jket3 /= 0) &
     & Call AppendBottom(jket3, ybra+ybra2, yket+yket3, FuncA(b,0,1)*value, h, j, k, l)  ! AppendLb

       If (h <= 1)  Return

       If (jbra0 /= 0 .And. jket0 /= 0) &
     & Call AppendL(jbra0, jket0, ybra+ybra0, yket+yket0,            value,           j, k, l)
       If (jbra1 /= 0 .And. jket1 /= 0) &
     & Call AppendL(jbra1, jket1, ybra+ybra1, yket+yket1,           -value,           j, k, l)
       If (jbra2 /= 0 .And. jket1 /= 0) &
     & Call AppendL(jbra2, jket1, ybra+ybra2, yket+yket1,           -value/Dble(b+1), j, k, l)
       If (jbra2 /= 0 .And. jket2 /= 0) &
     & Call AppendL(jbra2, jket2, ybra+ybra2, yket+yket2, FuncC(b,1)*value,           j, k, l)
       If (jbra3 /= 0 .And. jket3 /= 0) &
     & Call AppendL(jbra3, jket3, ybra+ybra3, yket+yket3,           -value,           j, k, l)

    End If

    Return
  End Subroutine AppendL


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the bottom segments  -------------
!------------------------------------------------------------------------------

  Subroutine AppendLowerW(irow, ybra, yket, value, i, k, l)
    Implicit None
    Integer          :: irow, ybra, yket
    Double Precision :: value
    Integer          :: i, k, l
    ! End of dummy parameters.

    If (Dabs(value) < Tiny)  Return

    Call TwoBodyLoopComplete(irow, ybra, yket, value, i, i, k, l)

    Return
  End Subroutine AppendLowerW



  Subroutine AppendBottom(irow, ybra, yket, value, i, j, k, l)
    Implicit None
    Integer          :: irow, ybra, yket
    Double Precision :: value
    Integer          :: i, j, k, l
    ! End of dummy parameters.

    If (Dabs(value) < Tiny .Or. oneOnly)  Return

    Call TwoBodyLoopComplete(irow, ybra, yket, value, i, j, k, l)

    Return
  End Subroutine AppendBottom


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the upper and lower walks  -------
!------------------------------------------------------------------------------

  Recursive Subroutine UpperWalks(irow, y)
    Implicit None
    Integer :: irow, y
    ! End of dummy parameters.
    Integer :: d

    If (irow > 1) Then
       Do d=3, 0, -1
          Call UpperWalks(DRT(irow)%ju(d), y+DRT(irow)%yu(d))
       End Do
    Else If (irow == 1) Then
       nupper = nupper + 1
       UpperVec(nupper) = y
    End If

    Return
  End Subroutine UpperWalks



  Subroutine OneBodyLoopComplete(ibot, ybra, yket, value, k, l)
    Implicit None
    Integer          :: ibot, ybra, yket
    Double Precision :: value
    Integer          :: k, l
    ! End of dummy parameters.
    Integer          :: i, j, ij, iupper, ilower
    Integer          :: ybra1, yket1, m, n, sybra, syket
    Double Precision :: cmi, temp

    ! WARNING:
    ! In job 1 - 4 the one-electron matrix
    ! is presumed to be symmetrical.

    Do iupper=1, nupper
       ybra1 = ybra + UpperVec(iupper)
       yket1 = yket + UpperVec(iupper)
       Do ilower=1, DRT(ibot)%yd(4)
          m = IndVec(ybra1 + ilower)
          n = IndVec(yket1 + ilower)
          If (m /= 0 .And. n /= 0) Then
             sybra = SymVec(m)
             syket = SymVec(n)
             If (allOneGen .Or. (sybra == syket .And. symFlag(sybra))) Then
                !
                ! For a totally symmetric operator like the Hamiltonian,
                ! matrix elements connecting functions of different symmetry
                ! will vanish and the corresponding generator matrix elements,
                ! although generally non-zero, will not contribute.
                !
                ! Thus for the electric and magnetic dipole moment operators,
                ! all contributions have to be taken into account.
                !
                Select Case (job)
                Case (0)
                   If (m <= n) Then
                      Call AddGen(OneList, m, n, value, k, l)
                   Else
                      Call AddGen(OneList, n, m, value, l, k)
                   End If
                Case (1)
                   m = m - index0
                   n = n - index0
                   CIP(m) = CIP(m) + CIV1(n) * OneInt(k,l) * value
                   If (m /= n) Then
                      CIP(n) = CIP(n) + CIV1(m) * OneInt(k,l) * value
                   End If
                Case (2)
                   m = m - index0
                   n = n - index0
                   CIPtr(:,m) = CIPtr(:,m) + CICtr(:,n) * OneInt(k,l) * value
                   If (m /= n) Then
                      CIPtr(:,n) = CIPtr(:,n) + CICtr(:,m) * OneInt(k,l) * value
                   End If
                Case (3)
                   If (m /= n) Then
                      Write(Stdout,'(A)') ' Unexpected off-diagonal element.'
                      Stop 'OneBodyLoopComplete'
                   End If
                   Hdiag(m) = Hdiag(m) + OneInt(k,l) * value
                Case (4)
                   OneDen(k,l)    = OneDen(k,l) + CIV1(m) * CIV2(n) * value
                   If (k /= l) Then
                      OneDen(l,k) = OneDen(l,k) + CIV1(n) * CIV2(m) * value
                   End If
                Case (5)
                   If (m /= n) Then
                      Do i=1, nPerm
                         Perm(k,l,i) = Perm(k,l,i) + CIC(m,i) * CIC(n,i) * value
                         Perm(l,k,i) = Perm(l,k,i) + CIC(m,i) * CIC(n,i) * value
                      End Do
                      ij = 1
                      Do i=1, mTrans
                         Do j=i+1, nTrans
                            Trans(k,l,ij) = Trans(k,l,ij) + CIC(m,i) * CIC(n,j) * value
                            Trans(l,k,ij) = Trans(l,k,ij) + CIC(n,i) * CIC(m,j) * value
                            ij = ij + 1
                         End Do
                      End Do
                   Else
                      Do i=1, nPerm
                         cmi = CIC(m,i)
                         Perm(k,k,i) = Perm(k,k,i) + cmi * cmi * value
                      End Do
                      ij = 1
                      Do i=1, mTrans
                         Do j=i+1, nTrans
                            Trans(k,k,ij) = Trans(k,k,ij) + CIC(m,i) * CIC(m,j) * value
                            ij = ij + 1
                         End Do
                      End Do
                   End If
                Case Default
                   Write(Stdout,'(A,I1,A)') ' Unexpected job number (', job, ').'
                   Stop 'OneBodyLoopComplete'
                End Select
             End If
          End If
       End Do
    End Do

    Return
  End Subroutine OneBodyLoopComplete



  Subroutine TwoBodyLoopComplete(ibot, ybra, yket, value, i, j, k, l)
    Implicit None
    Integer          :: ibot, ybra, yket
    Double Precision :: value
    Integer          :: i, j, k, l
    ! End of dummy parameters.
    Integer          :: ij, kl, ijkl, iupper, ilower
    Integer          :: ybra1, yket1, m, n, sybra, syket
    Double Precision :: temp

    ij = CanonicalIndex(i,j)
    kl = CanonicalIndex(k,l)

    Do iupper=1, nupper
       ybra1 = ybra + UpperVec(iupper)
       yket1 = yket + UpperVec(iupper)
       Do ilower=1, DRT(ibot)%yd(4)
          m  = IndVec(ybra1 + ilower)
          n  = IndVec(yket1 + ilower)
          If (m /= 0 .And. n /= 0) Then
             sybra = SymVec(m)
             syket = SymVec(n)
             If (allTwoGen .Or. (sybra == syket .And. symFlag(sybra))) Then
                !
                ! For a totally symmetric operator like the Hamiltonian,
                ! matrix elements connecting functions of different symmetry
                ! will vanish and the corresponding generator matrix elements,
                ! although generally non-zero, will not contribute.
                !
                ! For the calculation of non-adiabatic coupling elements however,
                ! all generator matrix elements need to be taken into account.
                !
                Select Case (job)
                Case (0)
                   Call AddGen(TwoList, m, n, value, ij, kl)
                Case (1)
                   m = m - index0
                   n = n - index0
                   CIP(m) = CIP(m) + CIV1(n) * TwoInt(ij,kl) * value
                   If (m /= n) Then
                      CIP(n) = CIP(n) + CIV1(m) * TwoInt(ij,kl) * value
                   End If
                Case (2)
                   m = m - index0
                   n = n - index0
                   CIPtr(:,m) = CIPtr(:,m) + CICtr(:,n) * TwoInt(ij,kl) * value
                   If (m /= n) Then
                      CIPtr(:,n) = CIPtr(:,n) + CICtr(:,m) * TwoInt(ij,kl) * value
                   End If
                Case (3)
                   If (m /= n) Then
                      Write(Stdout,'(A)') ' Unexpected off-diagonal element.'
                      Stop 'TwoBodyLoopComplete'
                   End If
                   Hdiag(m) = Hdiag(m) + TwoInt(ij,kl) * value
                Case (4)
                   temp = CIV1(m) * CIV2(n)
                   ijkl = CanonicalIndex(ij,kl)
                   If (m /= n)  temp = temp + CIV1(n) * CIV2(m)
                   TwoDen(ijkl) = TwoDen(ijkl) + temp * value
                Case Default
                   Write(Stdout,'(A,I1,A)') ' Unexpected job number (', job, ').'
                   Stop 'TwoBodyLoopComplete'
                End Select
             End If
          End If
       End Do
    End Do

    Return
  End Subroutine TwoBodyLoopComplete


!------------------------------------------------------------------------------
!------  Subroutine to add a generator matrix element to a buffer  ------------
!------------------------------------------------------------------------------

  Subroutine AddGen(list, m, n, value, ij, kl)
    Implicit None
    Type(GenBufList), Dimension(:) :: list
    Integer                        :: m, n  ! m <= n (upper triangle!)
    Double Precision               :: value
    Integer                        :: ij, kl
    ! End of dummy parameters.
    Type(GenBuffer), Pointer       :: buffer
    Integer                        :: ibuf

    ibuf = list(m)%first%count + 1
    If (ibuf <= genBufSize) Then
       list(m)%first%gen(ibuf) = Generator(value, n, ij, kl)
       list(m)%first%count     = ibuf
    Else
       Allocate(buffer)
       buffer%gen(1) =  Generator(value, n, ij, kl)
       buffer%count  =  1
       buffer%next   => list(m)%first
       list(m)%first => buffer
       list(m)%count =  list(m)%count + 1
    End If

    Return
  End Subroutine AddGen


!------------------------------------------------------------------------------
!------  Auxiliary functions, eqs. (49), (50), (62), and (63)  ----------------
!------------------------------------------------------------------------------

  ! The following functions are designed according to these guidelines:
  !
  ! 1) Prefer integer addition over floating point addition, because
  !    integer addition has a latency of 1 clock cycle on Alpha 21264
  !    and AMD Athlon processors (assuming register or direct operands)
  !    compared to a latency of 4 for a floating point addition.
  !
  ! 2) Avoid integer multiplications, because they have a higher latency (7)
  !    than floating point multiplications (4) on the Alpha 21264, and are
  !    VectorPath operations on the AMD Athlon. Integer to floating point
  !    conversion takes 4 clock cycles on both processors, and several
  !    conversions can be done in parallel due to the superscalar
  !    architectures.

  Function FuncA(b,p,q)
    Implicit None
    Double Precision :: FuncA
    Integer          :: b, p, q
    FuncA = DSqrt(Dble(b+p) / Dble(b+q))
    Return
  End Function FuncA


  Function FuncB(b,p,q)
    Implicit None
    Double Precision :: FuncB
    Integer          :: b, p, q
    FuncB = DSqrt(2.D0 / (Dble(b+p) * Dble(b+q)))
    Return
  End Function FuncB


  Function FuncC(b,p)
    Implicit None
    Double Precision :: FuncC
    Integer          :: b, p, s
    s     = b + p
    FuncC = DSqrt(Dble(s-1) * Dble(s+1)) / Dble(s)
    Return
  End Function FuncC


  Function FuncD(b,p)
    Implicit None
    Double Precision :: FuncD
    Integer          :: b, p, s
    s     = b + p
    FuncD = DSqrt((Dble(s-1) * Dble(s+2)) / (Dble(s) * Dble(s+1)))
    Return
  End Function FuncD

End Module gugaloops
