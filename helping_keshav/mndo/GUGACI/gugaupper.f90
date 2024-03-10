!
! This module provides subroutines for the top-down construction of
! partial and complete one- and two-body loops in the upper part of
! the DRT using the Hamiltonian repartitioning scheme
! (I. Shavitt, Lecture Notes in Chemistry 22 (1981), 51-99).
!
! Adapted from gugaloops.f90 by Axel Koslowski in 2012-2013 at MPI Muelheim.
!

Module gugaupper

Use gugaglobal
Use gugaloops, Only: Sqrt2, t, Tiny, FuncA, FuncB, FuncC, FuncD
Use gugautils, Only: PrintTime, CanonicalIndex

Implicit None

Private

Public :: CalcUpperPartialLoops  ! Calculate and store upper partial loops and complete loops.


!------------------------------------------------------------------------------
!------  Global variables  ----------------------------------------------------
!------------------------------------------------------------------------------

  Double Precision,   Parameter   :: Zero = 0.0D0
  Double Precision,   Parameter   :: Pt5  = 0.5D0
  Double Precision,   Parameter   :: One  = 1.0D0
  Double Precision,   Parameter   :: Two  = 2.0D0
  Type(PaldusRow),    Pointer     :: DRT(:)      ! Pointer to DRT.
  Double Precision,   Allocatable :: Occ(:)      ! Occupation numbers of active orbitals.
  Integer,            Allocatable :: Irrep(:)    ! Irrep of each active MO.
  Type(LoopBufLists), Pointer     :: lists       ! Provides global access to all lists.
  Integer                         :: itop        ! Top row of the current loop.
  Integer                         :: kBorder     ! Orbital level that separates upper and
                                                 ! lower part of the DRT.

Contains

!------------------------------------------------------------------------------
!------  Driver routine  ------------------------------------------------------
!------------------------------------------------------------------------------

  ! The following subroutine is the public driver routine for the generation
  ! of the partial and complete one- and two-body loops in the upper part of
  ! the DRT. The basic procedure of loop construction is described in:
  !
  ! F. A. Matsen, R. Pauncz, The Unitary Group in Quantum Chemistry, Elsevier
  ! Amsterdam (1986).
  !
  ! I. Shavitt, The graphical unitary group approach and its application to
  ! direct configuration interaction calculations,
  ! Lecture Notes in Chemistry 22 (1981), 51-99.
  !
  ! Our implementation works by constructing all possible one- and two-electron
  ! loops starting at an arbitary vertex in the upper part of the Shavitt graph,
  ! simply by adding segments in all possible ways. After the border between upper
  ! and lower part of the DRT has been reached, the partial loop is stored into the
  ! corresponding array. If a loop is completed before the border has been reached,
  ! the complete loop is stored. This procedure is done for all vertices in the
  ! upper part of the Shavitt graph.
  !


  ! bufLists has to be properly initialized
  ! before calling subroutine CalcUpperPartialLoops!
  !
  Subroutine CalcUpperPartialLoops(Flags, bufLists, Occup)
    Implicit None
    Type(ShavittControl)             :: Flags
    Type(LoopBufLists), Target       :: bufLists
    Double Precision,   Dimension(:) :: Occup
    ! End of dummy parameters.
    Double Precision     :: tuser, tsys, twall

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(/1X,"Calculating upper partial loops...")')
       Call GetTime(tuser, tsys, twall)
    End If

    DRT   => Flags%DRT
    lists => bufLists
    Allocate(Occ(DRT(1)%k), Irrep(DRT(1)%k))
    Occ(:)   = Occup(:)
    Irrep(:) = Flags%MOIrreps(:)
    kBorder  = Flags%kBorder
    Do itop=1, Size(DRT)-1
       If (DRT(itop)%k == kBorder) Exit
       Call StartRtLoop
       Call StartWLoop
       Call StartRtRtLoop
       Call StartRtLtLoop
       Call StartWRtLoop
       Call StartWWLoop
    End Do
    Deallocate(Occ, Irrep)

    ! Print timing of this subroutine if requested.
    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)
  End Subroutine CalcUpperPartialLoops



!------------------------------------------------------------------------------
!------  Six different ways to start a loop  ----------------------------------
!------------------------------------------------------------------------------

  Subroutine StartWWLoop
    Implicit None
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3, h, sy

    j0 = DRT(itop)%jd(0)
    j1 = DRT(itop)%jd(1)
    j2 = DRT(itop)%jd(2)
    j3 = DRT(itop)%jd(3)

    h  = DRT(itop)%k
    uu = Occ(h)
    sy = Irrep(h)

    ! These are the modified WW segment values published in
    ! A. Koslowski, M. E. Beck, W. Thiel, JCC 24, 714-726 (2003).
    If (j0 /= 0)  Call AppendWW(j0, DRT(itop)%yd(0), 1,  Pt5* uu     * uu,      h)
    If (j1 /= 0)  Call AppendWW(j1, DRT(itop)%yd(1), sy, Pt5* uu     *(uu-Two), h)
    If (j2 /= 0)  Call AppendWW(j2, DRT(itop)%yd(2), sy, Pt5* uu     *(uu-Two), h)
    If (j3 /= 0)  Call AppendWW(j3, DRT(itop)%yd(3), 1,  Pt5*(uu-Two)*(uu-Two), h)

    Return
  End Subroutine StartWWLoop



  Subroutine StartWLoop
    Implicit None
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3, h, sy

    j0 = DRT(itop)%jd(0)
    j1 = DRT(itop)%jd(1)
    j2 = DRT(itop)%jd(2)
    j3 = DRT(itop)%jd(3)

    h  = DRT(itop)%k
    uu = Occ(h)
    sy = Irrep(h)

    If (j0 /= 0)  Call AppendUpperW(j0, DRT(itop)%yd(0), 1,     -uu, h)
    If (j1 /= 0)  Call AppendUpperW(j1, DRT(itop)%yd(1), sy, One-uu, h)
    If (j2 /= 0)  Call AppendUpperW(j2, DRT(itop)%yd(2), sy, One-uu, h)
    If (j3 /= 0)  Call AppendUpperW(j3, DRT(itop)%yd(3), 1,  Two-uu, h)

    Return
  End Subroutine StartWLoop



  Subroutine StartWRtLoop
    Implicit None
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3
    Integer          :: y0, y1, y2, y3
    Integer          :: b, h, sy

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
    sy = Irrep(h)

    If (j0 /= 0 .And. j1 /= 0)  Call AppendLowerR(j0, j1, y0, y1, 1,  sy,     -Pt5*uu,               h, h, h)  ! AppendWRt
    If (j0 /= 0 .And. j2 /= 0)  Call AppendLowerR(j0, j2, y0, y2, 1,  sy,     -Pt5*uu,               h, h, h)  ! AppendWRt
    If (j1 /= 0 .And. j3 /= 0)  Call AppendLowerR(j1, j3, y1, y3, sy, 1,  (One-Pt5*uu)*FuncA(b,0,1), h, h, h)  ! AppendWRt
    If (j2 /= 0 .And. j3 /= 0)  Call AppendLowerR(j2, j3, y2, y3, sy, 1,  (One-Pt5*uu)*FuncA(b,2,1), h, h, h)  ! AppendWRt

    Return
  End Subroutine StartWRtLoop



  Subroutine StartRtRtLoop
    Implicit None
    Integer :: j0, j3, y0, y3, h

    j0 = DRT(itop)%jd(0)
    j3 = DRT(itop)%jd(3)

    y0 = DRT(itop)%yd(0)
    y3 = DRT(itop)%yd(3)

    h  = DRT(itop)%k

    If (j0 /= 0 .And. j3 /= 0)  Call AppendRR(j0, j3, y0, y3, 1, 1, Sqrt2, Zero, h, h, iUpperRtRt)  ! AppendRtRt

    Return
  End Subroutine StartRtRtLoop



  Subroutine StartRtLtLoop
    Implicit None
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3
    Integer          :: y0, y1, y2, y3
    Integer          :: b, h, sy

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
    sy = Irrep(h)

    If (j0 /= 0)  Call AppendRL(j0, j0, y0, y0, 1,  1,  -t*     uu,     Zero,          h, h, iUpperRLdg)  ! AppendRtLt
    If (j1 /= 0)  Call AppendRL(j1, j1, y1, y1, sy, sy,  t*(One-uu), -t*FuncA(b,-1,1), h, h, iUpperRLdg)  ! AppendRtLt
    If (j1 /= 0 .And. j2 /= 0) &
                  Call AppendRL(j1, j2, y1, y2, sy, sy,     Zero,       One,           h, h, iUpperRL)    ! AppendRtLt
    If (j2 /= 0)  Call AppendRL(j2, j2, y2, y2, sy, sy,  t*(One-uu),  t*FuncA(b, 3,1), h, h, iUpperRLdg)  ! AppendRtLt
    If (j3 /= 0)  Call AppendRL(j3, j3, y3, y3, 1,  1,   t*(Two-uu),    Zero,          h, h, iUpperRLdg)  ! AppendRtLt

    Return
  End Subroutine StartRtLtLoop



  Subroutine StartRtLoop
    Implicit None
    Integer :: j0, j1, j2, j3
    Integer :: y0, y1, y2, y3
    Integer :: b, h, sy

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
    sy = Irrep(h)

    If (j0 /= 0 .And. j1 /= 0)  Call AppendUpperR(j0, j1, y0, y1, 1,  sy, One,          h)  ! AppendRt
    If (j0 /= 0 .And. j2 /= 0)  Call AppendUpperR(j0, j2, y0, y2, 1,  sy, One,          h)  ! AppendRt
    If (j1 /= 0 .And. j3 /= 0)  Call AppendUpperR(j1, j3, y1, y3, sy, 1,  FuncA(b,0,1), h)  ! AppendRt
    If (j2 /= 0 .And. j3 /= 0)  Call AppendUpperR(j2, j3, y2, y3, sy, 1,  FuncA(b,2,1), h)  ! AppendRt

    Return
  End Subroutine StartRtLoop


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the top segments  ----------------
!------------------------------------------------------------------------------

  Subroutine AppendWW(ibot, y, sy, value, l)
    Implicit None
    Integer          :: ibot, y, sy
    Double Precision :: value
    Integer          :: l
    ! End of dummy parameters.

    If (Dabs(value) < Tiny)  Return

    Call TwoBodyLoopComplete(ibot, y, y, sy, sy, Pt5*value, l, l, l, l)

    Return
  End Subroutine AppendWW



  Subroutine AppendUpperW(ibot, y, sy, value, l)
    Implicit None
    Integer          :: ibot, y, sy
    Double Precision :: value
    Integer          :: l
    ! End of dummy parameters.

    If (Dabs(value) < Tiny)  Return

    Call OneBodyLoopComplete(ibot, y, y, sy, sy, value, l, l)
    Call WalkDown           (ibot, y, y, sy, sy, value, l, l, iUpperW)

    Return
  End Subroutine AppendUpperW


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the upper joining segments  ------
!------------------------------------------------------------------------------

  Recursive Subroutine AppendUpperR(ibra, iket, ybra, yket, sybra, syket, value, l)
    Implicit None
    Integer          :: ibra, iket, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: l
    ! End of dummy parameters.
    Double Precision :: uu
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab, prbra, prket
    Double Precision :: x

    If (Dabs(value) < Tiny)  Return

    h      =     DRT(iket)%k

    If (h == kBorder) Then
       Call AddPartialLoop(ibra, iket, ybra, yket, sybra, syket, value, Zero, 0, 0, l, iUpperRt)
       Return
    End If

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

    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b
    uu     =     Occ(h)
    prbra  =     iProd(sybra, Irrep(h))
    prket  =     iProd(syket, Irrep(h))

    If (deltab == -1) Then

       If (jbra1 == jket0 .And. jket0 /= 0) Then
          Call AppendUpperRb(jket0, ybra+ybra1, yket+yket0, prbra, syket,                           value, h, l)
          Call AppendBottom (jket0, ybra+ybra1, yket+yket0, prbra, syket,                   -Pt5*uu*value, h, h, h, l)  ! AppendRbW
       End If

       If (jbra3 == jket2 .And. jket2 /= 0) Then
          Call AppendUpperRb(jket2, ybra+ybra3, yket+yket2, sybra, prket,              FuncA(b,1,2)*value, h, l)
          Call AppendBottom (jket2, ybra+ybra3, yket+yket2, sybra, prket, (One-Pt5*uu)*FuncA(b,1,2)*value, h, h, h, l)  ! AppendRbW
       End If

       If (jbra0 /= 0 .And. jket0 /= 0) Then
          Call AppendUpperR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,                     value, l)
          Call AppendLowerR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,                 -uu*value, l, h, h)  ! AppendRW
          Call AppendLowerR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,              Pt5*uu*value, h, h, l)  ! AppendRbRt
       End If

       If (jbra1 /= 0 .And. jket1 /= 0) Then
          Call AppendUpperR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket,                    -value, l)
          Call AppendLowerR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket,            (uu-One)*value, l, h, h)  ! AppendRW
          Call AppendLowerR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket,        (One-Pt5*uu)*value, h, h, l)  ! AppendRbRt
       End If

       If (jbra1 /= 0 .And. jket2 /= 0) Then
          x = One / Dble(b+2)
          Call AppendUpperR(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket,                  -x*value, l)
          Call AppendLowerR(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket,          (uu-One)*x*value, l, h, h)  ! AppendRW
          Call AppendLowerR(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket,      (One-Pt5*uu*x)*value, h, h, l)  ! AppendRbRt
       End If

       If (jbra2 /= 0 .And. jket2 /= 0) Then
          Call AppendUpperR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket,          FuncC(b,2)*value, l)
          Call AppendLowerR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, (One-uu)*FuncC(b,2)*value, l, h, h)  ! AppendRW
          Call AppendLowerR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket,   Pt5*uu*FuncC(b,2)*value, h, h, l)  ! AppendRbRt
       End If

       If (jbra3 /= 0 .And. jket3 /= 0) Then
          Call AppendUpperR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,                    -value, l)
          Call AppendLowerR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,            (uu-Two)*value, l, h, h)  ! AppendRW
          Call AppendLowerR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,        (One-Pt5*uu)*value, h, h, l)  ! AppendRbRt
       End If

       If (jbra3 /= 0 .And. jket0 /= 0) &
          Call AppendL(jbra3, jket0, ybra+ybra3, yket+yket0, sybra, syket, FuncA(b,1,2)*value, h, h, l)  ! AppendRbLt

       If (jbra1 /= 0 .And. jket0 /= 0) &
          Call AppendRL(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket,              t*value, -t*FuncA(b,0,2)*value, h, l, iUpperRL)  ! AppendRLt
       If (jbra2 /= 0 .And. jket0 /= 0) &
          Call AppendRL(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket,                 Zero,                 value, h, l, iUpperRL)  ! AppendRLt
       If (jbra3 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket,                 Zero,   -FuncA(b,1,2)*value, h, l, iUpperRL)  ! AppendRLt
       If (jbra3 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket, t*FuncA(b,1,2)*value, -t*FuncA(b,3,2)*value, h, l, iUpperRL)  ! AppendRLt

       If (jbra0 /= 0 .And. jket1 /= 0) &
          Call AppendRR(jbra0, jket1, ybra+ybra0, yket+yket1, sybra, prket,                 Zero,                 value, h, l, iUpperRRt)  ! AppendRRt
       If (jbra0 /= 0 .And. jket2 /= 0) &
          Call AppendRR(jbra0, jket2, ybra+ybra0, yket+yket2, sybra, prket, t*FuncA(b,1,2)*value,  t*FuncA(b,3,2)*value, h, l, iUpperRRt)  ! AppendRRt
       If (jbra1 /= 0 .And. jket3 /= 0) &
          Call AppendRR(jbra1, jket3, ybra+ybra1, yket+yket3, prbra, syket,             -t*value, -t*FuncA(b,0,2)*value, h, l, iUpperRRt)  ! AppendRRt
       If (jbra2 /= 0 .And. jket3 /= 0) &
          Call AppendRR(jbra2, jket3, ybra+ybra2, yket+yket3, prbra, syket,                 Zero,    FuncA(b,3,2)*value, h, l, iUpperRRt)  ! AppendRRt

    Else If (deltab == 1) Then

       If (jbra2 == jket0 .And. jket0 /= 0) Then
          Call AppendUpperRb(jket0, ybra+ybra2, yket+yket0, prbra, syket,                           value, h, l)
          Call AppendBottom (jket0, ybra+ybra2, yket+yket0, prbra, syket,                   -Pt5*uu*value, h, h, h, l)  ! AppendRbW
       End If

       If (jbra3 == jket1 .And. jket1 /= 0) Then
          Call AppendUpperRb(jket1, ybra+ybra3, yket+yket1, sybra, prket,              FuncA(b,1,0)*value, h, l)
          Call AppendBottom (jket1, ybra+ybra3, yket+yket1, sybra, prket, (One-Pt5*uu)*FuncA(b,1,0)*value, h, h, h, l)  ! AppendRbW
       End If

       If (jbra0 /= 0 .And. jket0 /= 0) Then
          Call AppendUpperR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,                     value, l)
          Call AppendLowerR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,                 -uu*value, l, h, h)  ! AppendRW
          Call AppendLowerR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,              Pt5*uu*value, h, h, l)  ! AppendRbRt
       End If

       If (jbra1 /= 0 .And. jket1 /= 0) Then
          Call AppendUpperR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket,          FuncC(b,0)*value, l)
          Call AppendLowerR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, (One-uu)*FuncC(b,0)*value, l, h, h)  ! AppendRW
          Call AppendLowerR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket,   Pt5*uu*FuncC(b,0)*value, h, h, l)  ! AppendRbRt
       End If

       If (jbra2 /= 0 .And. jket1 /= 0) Then
          x = One / Dble(b)
          Call AppendUpperR(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket,                   x*value, l)
          Call AppendLowerR(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket,          (One-uu)*x*value, l, h, h)  ! AppendRW
          Call AppendLowerR(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket,      (One+Pt5*uu*x)*value, h, h, l)  ! AppendRbRt
       End If

       If (jbra2 /= 0 .And. jket2 /= 0) Then
          Call AppendUpperR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket,                    -value, l)
          Call AppendLowerR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket,            (uu-One)*value, l, h, h)  ! AppendRW
          Call AppendLowerR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket,        (One-Pt5*uu)*value, h, h, l)  ! AppendRbRt
       End If

       If (jbra3 /= 0 .And. jket3 /= 0) Then
          Call AppendUpperR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,                    -value, l)
          Call AppendLowerR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,            (uu-Two)*value, l, h, h)  ! AppendRW
          Call AppendLowerR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,        (One-Pt5*uu)*value, h, h, l)  ! AppendRbRt
       End If

       If (jbra3 /= 0 .And. jket0 /= 0) &
          Call AppendL(jbra3, jket0, ybra+ybra3, yket+yket0, sybra, syket, FuncA(b,1,0)*value, h, h, l)  ! AppendRbLt

       If (jbra1 /= 0 .And. jket0 /= 0) &
          Call AppendRL(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket,                 Zero,                  value, h, l, iUpperRL)  ! AppendRLt
       If (jbra2 /= 0 .And. jket0 /= 0) &
          Call AppendRL(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket,              t*value,  t*FuncA(b, 2,0)*value, h, l, iUpperRL)  ! AppendRLt
       If (jbra3 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket, t*FuncA(b,1,0)*value,  t*FuncA(b,-1,0)*value, h, l, iUpperRL)  ! AppendRLt
       If (jbra3 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket,                 Zero,   -FuncA(b, 1,0)*value, h, l, iUpperRL)  ! AppendRLt

       If (jbra0 /= 0 .And. jket1 /= 0) &
          Call AppendRR(jbra0, jket1, ybra+ybra0, yket+yket1, sybra, prket, t*FuncA(b,1,0)*value, -t*FuncA(b,-1,0)*value, h, l, iUpperRRt)  ! AppendRRt
       If (jbra0 /= 0 .And. jket2 /= 0) &
          Call AppendRR(jbra0, jket2, ybra+ybra0, yket+yket2, sybra, prket,                 Zero,                  value, h, l, iUpperRRt)  ! AppendRRt
       If (jbra1 /= 0 .And. jket3 /= 0) &
          Call AppendRR(jbra1, jket3, ybra+ybra1, yket+yket3, prbra, syket,                 Zero,    FuncA(b,-1,0)*value, h, l, iUpperRRt)  ! AppendRRt
       If (jbra2 /= 0 .And. jket3 /= 0) &
          Call AppendRR(jbra2, jket3, ybra+ybra2, yket+yket3, prbra, syket,             -t*value,  t*FuncA(b, 2,0)*value, h, l, iUpperRRt)  ! AppendRRt

    End If

    Return
  End Subroutine AppendUpperR


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the upper middle segments  -------
!------------------------------------------------------------------------------

  Subroutine AppendUpperRb(ibot, ybra, yket, sybra, syket, value, k, l)
    Implicit None
    Integer          :: ibot, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: k, l
    ! End of dummy parameters.

    If (Dabs(value) < Tiny) Return

    Call OneBodyLoopComplete(ibot, ybra, yket, sybra, syket, value, k, l)
    Call WalkDown           (ibot, ybra, yket, sybra, syket, value, k, l, iUpperOneR)

    Return
  End Subroutine AppendUpperRb


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the middle joining segments  -----
!------------------------------------------------------------------------------

  Recursive Subroutine WalkDown(irow, ybra, yket, sybra, syket, value, k, l, itype)
    Implicit None
    Integer          :: irow, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: k, l, itype
    ! End of dummy parameters.
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3, y0, y1, y2, y3, h, b, prbra, prket

    h     = DRT(irow)%k

    If (h == kBorder) Then
       Call AddPartialLoop(irow, irow, ybra, yket, sybra, syket, value, Zero, 0, k, l, itype)
       Return
    End If

    j0    = DRT(irow)%jd(0)
    j1    = DRT(irow)%jd(1)
    j2    = DRT(irow)%jd(2)
    j3    = DRT(irow)%jd(3)

    y0    = DRT(irow)%yd(0)
    y1    = DRT(irow)%yd(1)
    y2    = DRT(irow)%yd(2)
    y3    = DRT(irow)%yd(3)

    b     = DRT(irow)%b
    uu    = Occ(h)
    prbra = iProd(sybra, Irrep(h))
    prket = iProd(syket, Irrep(h))

    If (j0 /= 0)  Call AppendBottom(j0, ybra+y0, yket+y0, sybra, syket,     -uu *value, h, h, k, l)  ! AppendLowerW
    If (j1 /= 0)  Call AppendBottom(j1, ybra+y1, yket+y1, prbra, prket, (One-uu)*value, h, h, k, l)  ! AppendLowerW
    If (j2 /= 0)  Call AppendBottom(j2, ybra+y2, yket+y2, prbra, prket, (One-uu)*value, h, h, k, l)  ! AppendLowerW
    If (j3 /= 0)  Call AppendBottom(j3, ybra+y3, yket+y3, sybra, syket, (Two-uu)*value, h, h, k, l)  ! AppendLowerW

    If (j0 /= 0)  Call WalkDown(j0, ybra+y0, yket+y0, sybra, syket, value, k, l, itype)
    If (j1 /= 0)  Call WalkDown(j1, ybra+y1, yket+y1, prbra, prket, value, k, l, itype)
    If (j2 /= 0)  Call WalkDown(j2, ybra+y2, yket+y2, prbra, prket, value, k, l, itype)
    If (j3 /= 0)  Call WalkDown(j3, ybra+y3, yket+y3, sybra, syket, value, k, l, itype)

    If (j0 /= 0 .And. j1 /= 0)  Call AppendLowerR(j0, j1, ybra+y0, yket+y1, sybra, prket,              value, h, k, l)  ! AppendLowerRt
    If (j0 /= 0 .And. j2 /= 0)  Call AppendLowerR(j0, j2, ybra+y0, yket+y2, sybra, prket,              value, h, k, l)  ! AppendLowerRt
    If (j1 /= 0 .And. j3 /= 0)  Call AppendLowerR(j1, j3, ybra+y1, yket+y3, prbra, syket, FuncA(b,0,1)*value, h, k, l)  ! AppendLowerRt
    If (j2 /= 0 .And. j3 /= 0)  Call AppendLowerR(j2, j3, ybra+y2, yket+y3, prbra, syket, FuncA(b,2,1)*value, h, k, l)  ! AppendLowerRt

    If (itype == iUpperW)  Return

    If (j1 /= 0 .And. j0 /= 0)  Call AppendL(j1, j0, ybra+y1, yket+y0, prbra, syket,              value, h, k, l)  ! AppendLt
    If (j2 /= 0 .And. j0 /= 0)  Call AppendL(j2, j0, ybra+y2, yket+y0, prbra, syket,              value, h, k, l)  ! AppendLt
    If (j3 /= 0 .And. j1 /= 0)  Call AppendL(j3, j1, ybra+y3, yket+y1, sybra, prket, FuncA(b,0,1)*value, h, k, l)  ! AppendLt
    If (j3 /= 0 .And. j2 /= 0)  Call AppendL(j3, j2, ybra+y3, yket+y2, sybra, prket, FuncA(b,2,1)*value, h, k, l)  ! AppendLt

    Return
  End Subroutine WalkDown



  Recursive Subroutine AppendRL(ibra, iket, ybra, yket, sybra, syket, sval, tval, k, l, itype)
    Implicit None
    Integer          :: ibra, iket, ybra, yket, sybra, syket
    Double Precision :: sval, tval
    Integer          :: k, l, itype
    ! End of dummy parameters.
    Double Precision :: uu
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab, prbra, prket, ityp

    If (Dabs(sval) < Tiny .And. Dabs(tval) < Tiny)  Return

    h      =     DRT(iket)%k

    If (h == kBorder) Then
       ityp = itype
       If (itype == iUpperRL .And. ibra == iket)  ityp = iUpperRLco
       Call AddPartialLoop(ibra, iket, ybra, yket, sybra, syket, sval, tval, 0, k, l, ityp)
       Return
    End If

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

    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b
    uu     =     Occ(h)
    prbra  =     iProd(sybra, Irrep(h))
    prket  =     iProd(syket, Irrep(h))

    If (deltab == -2) Then

       If (jbra1 == jket2 .And. jket2 /= 0) &
          Call AppendBottom(jket2, ybra+ybra1, yket+yket2, prbra, prket, FuncC(b,2)*tval, h, k, h, l)  ! AppendRbLb

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendRL(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket, Zero,             tval,           k, l, itype)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, Zero, -FuncC(b,2)*tval,           k, l, itype)
       If (jbra1 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket, Zero,      -Sqrt2*tval/Dble(b+2), k, l, iUpperRL)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, Zero, -FuncC(b,2)*tval,           k, l, itype)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendRL(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket, Zero,             tval,           k, l, itype)

       If (jbra0 /= 0 .And. jket2 /= 0) &
          Call AppendLowerR(jbra0, jket2, ybra+ybra0, yket+yket2, sybra, prket,    FuncC(b,2)*(sval+tval), l, h, k)  ! AppendRLb
       If (jbra1 /= 0 .And. jket3 /= 0) &
          Call AppendLowerR(jbra1, jket3, ybra+ybra1, yket+yket3, prbra, syket, -FuncA(b,3,2)*(sval+tval), l, h, k)  ! AppendRLb

       If (jbra1 /= 0 .And. jket0 /= 0) &
          Call AppendL(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket,    FuncC(b,2)*(sval+tval), k, h, l)  ! AppendRbL
       If (jbra3 /= 0 .And. jket2 /= 0) &
          Call AppendL(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket, -FuncA(b,1,2)*(sval+tval), k, h, l)  ! AppendRbL

    Else If (deltab == 0) Then

       If (jbra0 == jket0 .And. jket0 /= 0) &
          Call AppendBottom(jket0, ybra+ybra0, yket+yket0, sybra, syket, t*  uu     *sval,                    h, k, h, l)  ! AppendRbLb
       If (jbra1 == jket1 .And. jket1 /= 0) &
          Call AppendBottom(jket1, ybra+ybra1, yket+yket1, prbra, prket, t*((uu-One)*sval+FuncA(b,2,0)*tval), h, k, h, l)  ! AppendRbLb
       If (jbra2 == jket2 .And. jket2 /= 0) &
          Call AppendBottom(jket2, ybra+ybra2, yket+yket2, prbra, prket, t*((uu-One)*sval-FuncA(b,0,2)*tval), h, k, h, l)  ! AppendRbLb
       If (jbra3 == jket3 .And. jket3 /= 0) &
          Call AppendBottom(jket3, ybra+ybra3, yket+yket3, sybra, syket, t* (uu-Two)*sval,                    h, k, h, l)  ! AppendRbLb

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendRL(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket, sval,               tval, k, l, itype)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, sval,    FuncD(b,0)*tval, k, l, itype)
       If (jbra1 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket, Zero, -FuncB(b,0,2)*tval, k, l, iUpperRL)
       If (jbra2 /= 0 .And. jket1 /= 0 .And. itype /= iUpperRLdg) &
          Call AppendRL(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket, Zero, -FuncB(b,0,2)*tval, k, l, iUpperRL)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, sval,    FuncD(b,1)*tval, k, l, itype)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendRL(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket, sval,               tval, k, l, itype)

       If (jbra0 /= 0 .And. jket1 /= 0) &
          Call AppendLowerR(jbra0, jket1, ybra+ybra0, yket+yket1, sybra, prket,              -t*(sval-FuncA(b,2,0)*tval), l, h, k)  ! AppendRLb
       If (jbra0 /= 0 .And. jket2 /= 0) &                                          
          Call AppendLowerR(jbra0, jket2, ybra+ybra0, yket+yket2, sybra, prket,              -t*(sval+FuncA(b,0,2)*tval), l, h, k)  ! AppendRLb
       If (jbra1 /= 0 .And. jket3 /= 0) &
          Call AppendLowerR(jbra1, jket3, ybra+ybra1, yket+yket3, prbra, syket, -t*(FuncA(b,0,1)*sval+FuncA(b,2,1)*tval), l, h, k)  ! AppendRLb
       If (jbra2 /= 0 .And. jket3 /= 0) &
          Call AppendLowerR(jbra2, jket3, ybra+ybra2, yket+yket3, prbra, syket, -t*(FuncA(b,2,1)*sval-FuncA(b,0,1)*tval), l, h, k)  ! AppendRLb

       If (ybra < yket) Then
          If (jbra1 /= 0 .And. jket0 /= 0) &
             Call AppendL(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket,              -t*(sval-FuncA(b,2,0)*tval), k, h, l)  ! AppendRbL
          If (jbra2 /= 0 .And. jket0 /= 0) &                                     
             Call AppendL(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket,              -t*(sval+FuncA(b,0,2)*tval), k, h, l)  ! AppendRbL
          If (jbra3 /= 0 .And. jket1 /= 0) &
             Call AppendL(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket, -t*(FuncA(b,0,1)*sval+FuncA(b,2,1)*tval), k, h, l)  ! AppendRbL
          If (jbra3 /= 0 .And. jket2 /= 0) &
             Call AppendL(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket, -t*(FuncA(b,2,1)*sval-FuncA(b,0,1)*tval), k, h, l)  ! AppendRbL
       End If

    Else If (deltab == 2) Then

       If (jbra2 == jket1 .And. jket1 /= 0) &
          Call AppendBottom(jket1, ybra+ybra2, yket+yket1, prbra, prket, FuncC(b,0)*tval, h, k, h, l)  ! AppendRbLb

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendRL(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket, Zero,             tval,         k, l, itype)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, Zero, -FuncC(b,0)*tval,         k, l, itype)
       If (jbra2 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket, Zero,      -Sqrt2*tval/Dble(b), k, l, iUpperRL)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, Zero, -FuncC(b,0)*tval,         k, l, itype)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendRL(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket, Zero,             tval,         k, l, itype)

       If (jbra0 /= 0 .And. jket1 /= 0) &
          Call AppendLowerR(jbra0, jket1, ybra+ybra0, yket+yket1, sybra, prket,     FuncC(b,0)*tval, l, h, k)  ! AppendRLb
       If (jbra2 /= 0 .And. jket3 /= 0) &
          Call AppendLowerR(jbra2, jket3, ybra+ybra2, yket+yket3, prbra, syket, -FuncA(b,-1,0)*tval, l, h, k)  ! AppendRLb

       If (jbra2 /= 0 .And. jket0 /= 0) &
          Call AppendL(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket,    FuncC(b,0)*tval, k, h, l)  ! AppendRbL
       If (jbra3 /= 0 .And. jket1 /= 0) &
          Call AppendL(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket, -FuncA(b,1,0)*tval, k, h, l)  ! AppendRbL

    End If

    Return
  End Subroutine AppendRL



  Recursive Subroutine AppendRR(ibra, iket, ybra, yket, sybra, syket, sval, tval, k, l, itype)
    Implicit None
    Integer          :: ibra, iket, ybra, yket, sybra, syket
    Double Precision :: sval, tval
    Integer          :: k, l, itype
    ! End of dummy parameters.
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab, prbra, prket
    Double Precision :: x

    If (Dabs(sval) < Tiny .And. Dabs(tval) < Tiny)  Return

    h      =     DRT(iket)%k

    If (h == kBorder) Then
       Call AddPartialLoop(ibra, iket, ybra, yket, sybra, syket, sval, tval, 0, k, l, itype)
       Return
    End If

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

    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b
    prbra  =     iProd(sybra, Irrep(h))
    prket  =     iProd(syket, Irrep(h))

    If (deltab == -2) Then

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendRR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket, Zero,              tval, k, l, itype)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendRR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, Zero,              tval, k, l, itype)
       If (jbra1 /= 0 .And. jket2 /= 0) &
          Call AppendRR(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket, Zero, FuncB(b,2,3)*tval, k, l, itype)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendRR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, Zero,   FuncD(b,2)*tval, k, l, itype)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendRR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket, Zero,              tval, k, l, itype)

       If (jbra1 /= 0 .And. jket0 /= 0) &
          Call AppendLowerR(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket,               tval, k, h, l)  ! AppendRbR
       If (jbra3 /= 0 .And. jket2 /= 0) &
          Call AppendLowerR(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket,  FuncA(b,1,2)*tval, k, h, l)  ! AppendRbR

       If (itype == iUpperRtRt) Return

       If (jbra1 /= 0 .And. jket0 /= 0) &
          Call AppendLowerR(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket,              -tval, l, h, k)  ! AppendRRb
       If (jbra3 /= 0 .And. jket2 /= 0) &
          Call AppendLowerR(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket, -FuncA(b,1,2)*tval, l, h, k)  ! AppendRRb

    Else If (deltab == 0) Then

       If (jbra3 == jket0 .And. jket0 /= 0) Then
          x = One
          If (itype == iUpperRtRt) x = Pt5
          Call AppendBottom(jket0, ybra+ybra3, yket+yket0, sybra, syket, x*Sqrt2*sval, h, k, h, l)  ! AppendRbRb
       End If

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendRR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,  sval,              tval, k, l, itype)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendRR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, -sval,  -FuncD(b,0)*tval, k, l, itype)
       If (jbra1 /= 0 .And. jket2 /= 0) &
          Call AppendRR(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket,  Zero, FuncB(b,1,2)*tval, k, l, itype)
       If (jbra2 /= 0 .And. jket1 /= 0) &
          Call AppendRR(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket,  Zero, FuncB(b,0,1)*tval, k, l, itype)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendRR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, -sval,  -FuncD(b,1)*tval, k, l, itype)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendRR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,  sval,              tval, k, l, itype)

       If (jbra1 /= 0 .And. jket0 /= 0) &
          Call AppendLowerR(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket, t*(FuncA(b,0,1)*sval+FuncA(b,2,1)*tval), k, h, l)  ! AppendRbR
       If (jbra2 /= 0 .And. jket0 /= 0) &
          Call AppendLowerR(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket, t*(FuncA(b,2,1)*sval-FuncA(b,0,1)*tval), k, h, l)  ! AppendRbR
       If (jbra3 /= 0 .And. jket1 /= 0) &
          Call AppendLowerR(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket,             -t*(sval-FuncA(b,2,0)*tval), k, h, l)  ! AppendRbR
       If (jbra3 /= 0 .And. jket2 /= 0) &                                        
          Call AppendLowerR(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket,             -t*(sval+FuncA(b,0,2)*tval), k, h, l)  ! AppendRbR

       If (itype == iUpperRtRt) Return

       If (jbra1 /= 0 .And. jket0 /= 0) &
          Call AppendLowerR(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket, t*(FuncA(b,0,1)*sval-FuncA(b,2,1)*tval), l, h, k)  ! AppendRRb
       If (jbra2 /= 0 .And. jket0 /= 0) &
          Call AppendLowerR(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket, t*(FuncA(b,2,1)*sval+FuncA(b,0,1)*tval), l, h, k)  ! AppendRRb
       If (jbra3 /= 0 .And. jket1 /= 0) &
          Call AppendLowerR(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket,             -t*(sval+FuncA(b,2,0)*tval), l, h, k)  ! AppendRRb
       If (jbra3 /= 0 .And. jket2 /= 0) &                                                       
          Call AppendLowerR(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket,             -t*(sval-FuncA(b,0,2)*tval), l, h, k)  ! AppendRRb

    Else If (deltab == 2) Then

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendRR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket, Zero,               tval, k, l, itype)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendRR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, Zero,   FuncD(b,-1)*tval, k, l, itype)
       If (jbra2 /= 0 .And. jket1 /= 0) &
          Call AppendRR(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket, Zero, FuncB(b,-1,0)*tval, k, l, itype)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendRR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, Zero,               tval, k, l, itype)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendRR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket, Zero,               tval, k, l, itype)

       If (jbra2 /= 0 .And. jket0 /= 0) &
          Call AppendLowerR(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket,               tval, k, h, l)  ! AppendRbR
       If (jbra3 /= 0 .And. jket1 /= 0) &
          Call AppendLowerR(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket,  FuncA(b,1,0)*tval, k, h, l)  ! AppendRbR

       If (itype == iUpperRtRt) Return

       If (jbra2 /= 0 .And. jket0 /= 0) &
          Call AppendLowerR(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket,              -tval, l, h, k)  ! AppendRRb
       If (jbra3 /= 0 .And. jket1 /= 0) &
          Call AppendLowerR(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket, -FuncA(b,1,0)*tval, l, h, k)  ! AppendRRb
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

  Recursive Subroutine AppendLowerR(ibra, iket, ybra, yket, sybra, syket, value, j, k, l)
    Implicit None
    Integer          :: ibra, iket, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: j, k, l
    ! End of dummy parameters.
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab, prbra, prket

    If (Dabs(value) < Tiny) Return

    h      =     DRT(iket)%k

    If (h == kBorder) Then
       Call AddPartialLoop(ibra, iket, ybra, yket, sybra, syket, value, Zero, j, k, l, iUpperR)
       Return
    End If

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

    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b
    prbra  =     iProd(sybra, Irrep(h))
    prket  =     iProd(syket, Irrep(h))

    If (deltab == -1) Then

       If (jbra1 == jket0 .And. jket0 /= 0) &
          Call AppendBottom(jket0, ybra+ybra1, yket+yket0, prbra, syket,              value, h, j, k, l)  ! AppendLowerRb
       If (jbra3 == jket2 .And. jket2 /= 0) &
          Call AppendBottom(jket2, ybra+ybra3, yket+yket2, sybra, prket, FuncA(b,1,2)*value, h, j, k, l)  ! AppendLowerRb

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendLowerR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,            value,           j, k, l)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendLowerR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket,           -value,           j, k, l)
       If (jbra1 /= 0 .And. jket2 /= 0) &
          Call AppendLowerR(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket,           -value/Dble(b+2), j, k, l)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendLowerR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, FuncC(b,2)*value,           j, k, l)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendLowerR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,           -value,           j, k, l)

    Else If (deltab == 1) Then

       If (jbra2 == jket0 .And. jket0 /= 0) &
          Call AppendBottom(jket0, ybra+ybra2, yket+yket0, prbra, syket,              value, h, j, k, l)  ! AppendLowerRb
       If (jbra3 == jket1 .And. jket1 /= 0) &
          Call AppendBottom(jket1, ybra+ybra3, yket+yket1, sybra, prket, FuncA(b,1,0)*value, h, j, k, l)  ! AppendLowerRb

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendLowerR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,            value,         j, k, l)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendLowerR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, FuncC(b,0)*value,         j, k, l)
       If (jbra2 /= 0 .And. jket1 /= 0) &
          Call AppendLowerR(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket,            value/Dble(b), j, k, l)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendLowerR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket,           -value,         j, k, l)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendLowerR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,           -value,         j, k, l)

    End If

    Return
  End Subroutine AppendLowerR



  Recursive Subroutine AppendL(ibra, iket, ybra, yket, sybra, syket, value, j, k, l)
    Implicit None
    Integer          :: ibra, iket, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: j, k, l
    ! End of dummy parameters.
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab, prbra, prket

    If (Dabs(value) < Tiny)  Return

    h      =     DRT(iket)%k

    If (h == kBorder) Then
       Call AddPartialLoop(ibra, iket, ybra, yket, sybra, syket, value, Zero, j, k, l, iUpperL)
       Return
    End If

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

    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b
    prbra  =     iProd(sybra, Irrep(h))
    prket  =     iProd(syket, Irrep(h))

    If (deltab == -1) Then
       If (jbra0 == jket2 .And. jket2 /= 0) &
          Call AppendBottom(jket2, ybra+ybra0, yket+yket2, sybra, prket,              value, h, j, k, l)  ! AppendLb
       If (jbra1 == jket3 .And. jket3 /= 0) &
          Call AppendBottom(jket3, ybra+ybra1, yket+yket3, prbra, syket, FuncA(b,2,1)*value, h, j, k, l)  ! AppendLb

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendL(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,            value,           j, k, l)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendL(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, FuncC(b,1)*value,           j, k, l)
       If (jbra1 /= 0 .And. jket2 /= 0) &
          Call AppendL(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket,            value/Dble(b+1), j, k, l)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendL(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket,           -value,           j, k, l)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendL(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,           -value,           j, k, l)

    Else If (deltab == 1) Then
       If (jbra0 == jket1 .And. jket1 /= 0) &
          Call AppendBottom(jket1, ybra+ybra0, yket+yket1, sybra, prket,              value, h, j, k, l)  ! AppendLb
       If (jbra2 == jket3 .And. jket3 /= 0) &
          Call AppendBottom(jket3, ybra+ybra2, yket+yket3, prbra, syket, FuncA(b,0,1)*value, h, j, k, l)  ! AppendLb

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendL(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,            value,           j, k, l)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendL(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket,           -value,           j, k, l)
       If (jbra2 /= 0 .And. jket1 /= 0) &
          Call AppendL(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket,           -value/Dble(b+1), j, k, l)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendL(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, FuncC(b,1)*value,           j, k, l)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendL(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,           -value,           j, k, l)

    End If

    Return
  End Subroutine AppendL


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the bottom segments  -------------
!------------------------------------------------------------------------------

  Subroutine AppendBottom(ibot, ybra, yket, sybra, syket, value, i, j, k, l)
    Implicit None
    Integer          :: ibot, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: i, j, k, l
    ! End of dummy parameters.

    If (Dabs(value) < Tiny)  Return

    Call TwoBodyLoopComplete(ibot, ybra, yket, sybra, syket, value, i, j, k, l)

    Return
  End Subroutine AppendBottom


!------------------------------------------------------------------------------
!------  Subroutines to store partial and complete loops into a buffer  -------
!------------------------------------------------------------------------------

  Subroutine AddPartialLoop(ibra, iket, ybra, yket, sybra, syket, sval, tval, j, k, l, itype)
    Implicit None
    Integer          :: ibra, iket, ybra, yket, sybra, syket
    Double Precision :: sval, tval
    Integer          :: j, k, l, itype
    ! End of dummy parameters.
    Type(PartialLoopBuffer), Pointer :: newbuf

    If (.Not. Associated(lists%firstU)) Then
       Allocate(lists%firstU)
       Nullify(lists%firstU%next)
       lists%firstU%ncount = 1
    Else If (lists%firstU%ncount >= loopBufSize) Then
       Allocate(newbuf)
       newbuf%next   => lists%firstU
       lists%firstU  => newbuf
       newbuf%ncount =  1
    Else
       lists%firstU%ncount = lists%firstU%ncount + 1
    End If

    lists%firstU%part(lists%firstU%ncount)%singlet = sval
    lists%firstU%part(lists%firstU%ncount)%triplet = tval
    lists%firstU%part(lists%firstU%ncount)%itype   = itype
    lists%firstU%part(lists%firstU%ncount)%itip    = itop
    lists%firstU%part(lists%firstU%ncount)%ibra    = ibra
    lists%firstU%part(lists%firstU%ncount)%iket    = iket
    lists%firstU%part(lists%firstU%ncount)%lexbra  = ybra
    lists%firstU%part(lists%firstU%ncount)%lexket  = yket
    lists%firstU%part(lists%firstU%ncount)%sybra   = sybra
    lists%firstU%part(lists%firstU%ncount)%syket   = syket
    lists%firstU%part(lists%firstU%ncount)%i       = 0
    lists%firstU%part(lists%firstU%ncount)%j       = j
    lists%firstU%part(lists%firstU%ncount)%k       = k
    lists%firstU%part(lists%firstU%ncount)%l       = l
    lists%ncountU = lists%ncountU + 1

    If (itype == iUpperW    .Or. &
        itype == iUpperOneR .Or. &
        itype == iUpperR    .Or. &
        itype == iUpperL) Then
       lists%firstU%part(lists%firstU%ncount)%i = CanonicalIndex(k,l)
    End If
    
    Return
  End Subroutine AddPartialLoop



  Subroutine OneBodyLoopComplete(ibot, ybra, yket, sybra, syket, value, k, l)
    Implicit None
    Integer          :: ibot, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: k, l
    ! End of dummy parameters.
    Type(CompleteLoopBuffer), Pointer :: newbuf

    If (.Not. Associated(lists%first1)) Then
       Allocate(lists%first1)
       Nullify(lists%first1%next)
       lists%first1%ncount = 1
    Else If (lists%first1%ncount >= loopBufSize) Then
       Allocate(newbuf)
       newbuf%next   => lists%first1
       lists%first1  => newbuf
       newbuf%ncount =  1
    Else
       lists%first1%ncount = lists%first1%ncount + 1
    End If

    lists%first1%loop(lists%first1%ncount)%value  = value
    lists%first1%loop(lists%first1%ncount)%itop   = itop
    lists%first1%loop(lists%first1%ncount)%ibot   = ibot
    lists%first1%loop(lists%first1%ncount)%lexbra = ybra
    lists%first1%loop(lists%first1%ncount)%lexket = yket
    lists%first1%loop(lists%first1%ncount)%sybra  = sybra
    lists%first1%loop(lists%first1%ncount)%syket  = syket
    lists%first1%loop(lists%first1%ncount)%p      = k
    lists%first1%loop(lists%first1%ncount)%q      = l
    lists%ncount1 = lists%ncount1 + 1

    Return
  End Subroutine OneBodyLoopComplete



  Subroutine TwoBodyLoopComplete(ibot, ybra, yket, sybra, syket, value, i, j, k, l)
    Implicit None
    Integer          :: ibot, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: i, j, k, l
    ! End of dummy parameters.
    Type(CompleteLoopBuffer), Pointer :: newbuf

    If (.Not. Associated(lists%first2)) Then
       Allocate(lists%first2)
       Nullify(lists%first2%next)
       lists%first2%ncount = 1
    Else If (lists%first2%ncount >= loopBufSize) Then
       Allocate(newbuf)
       newbuf%next   => lists%first2
       lists%first2  => newbuf
       newbuf%ncount =  1
    Else
       lists%first2%ncount = lists%first2%ncount + 1
    End If

    lists%first2%loop(lists%first2%ncount)%value  = value
    lists%first2%loop(lists%first2%ncount)%itop   = itop
    lists%first2%loop(lists%first2%ncount)%ibot   = ibot
    lists%first2%loop(lists%first2%ncount)%lexbra = ybra
    lists%first2%loop(lists%first2%ncount)%lexket = yket
    lists%first2%loop(lists%first2%ncount)%sybra  = sybra
    lists%first2%loop(lists%first2%ncount)%syket  = syket
    lists%first2%loop(lists%first2%ncount)%p      = CanonicalIndex(i,j)
    lists%first2%loop(lists%first2%ncount)%q      = CanonicalIndex(k,l)
    lists%ncount2 = lists%ncount2 + 1

    Return
  End Subroutine TwoBodyLoopComplete


End Module gugaupper
