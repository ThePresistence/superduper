!
! This module provides subroutines for the bottom-up construction of
! partial and complete one- and two-body loops in the lower part of
! the DRT using the Hamiltonian repartitioning scheme
! (I. Shavitt, Lecture Notes in Chemistry 22 (1981), 51-99).
!
! Written by Axel Koslowski in 2012-2013 at MPI Muelheim.
!

Module gugalower

Use gugaglobal
Use gugaloops, Only: Sqrt2, t, Tiny, FuncA, FuncB, FuncC, FuncD
Use gugautils, Only: PrintTime, CanonicalIndex

Implicit None

Private

Public :: CalcLowerPartialLoops  ! Calculate and store upper partial loops and complete loops.


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
  Integer                         :: ibot        ! Bottom row of the current loop.
  Integer                         :: kBorder     ! Orbital level that separates upper and
                                                 ! lower part of the DRT.
Contains

!------------------------------------------------------------------------------
!------  Driver routine  ------------------------------------------------------
!------------------------------------------------------------------------------

  ! The following subroutine is the public driver routine for the generation
  ! of the partial and complete one- and two-body loops in the lower part of
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
  ! loops starting at an arbitary vertex in the lower part of the Shavitt graph,
  ! simply by adding segments in all possible ways above the previous segment.
  ! To keep the code as simple as possible, the tables containing the segment
  ! values have been converted such that they depend on the parameters of the
  ! orbital level at the bottom of the segment (in contrast to the original
  ! tables in which the segment values depend on the parameters of the orbital
  ! level at the top of the segment). This makes the new tables optimally
  ! suited for bottom-up loop construction. After the border between upper and
  ! lower part of the DRT has been reached, the partial loop is stored into
  ! the corresponding array. If a loop is completed before the border has been
  ! reached, the complete loop is stored. This procedure is done for all
  ! vertices in the lower part of the Shavitt graph.
  !


  ! bufLists has to be properly initialized
  ! before calling subroutine CalcLowerPartialLoops!
  !
  Subroutine CalcLowerPartialLoops(Flags, bufLists, Occup)
    Implicit None
    Type(ShavittControl)             :: Flags
    Type(LoopBufLists), Target       :: bufLists
    Double Precision,   Dimension(:) :: Occup
    ! End of dummy parameters.
    Double Precision     :: tuser, tsys, twall

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(/1X,"Calculating lower partial loops...")')
       Call GetTime(tuser, tsys, twall)
    End If

    DRT   => Flags%DRT
    lists => bufLists
    Allocate(Occ(DRT(1)%k), Irrep(DRT(1)%k))
    Occ(:)   = Occup(:)
    Irrep(:) = Flags%MOIrreps(:)
    kBorder  = Flags%kBorder
    Do ibot=Size(DRT), 1, -1
       If (DRT(ibot)%k == kBorder) Exit
       Call StartWWLoop
       Call StartWLoop
       Call StartWRbLoop
       Call StartRbRbLoop
       Call StartRbLbLoop
       Call StartRbLoop
       Call StartLbLoop
    End Do
    Deallocate(Occ, Irrep)

    ! Print timing of this subroutine if requested.
    If (Flags%PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall)
  End Subroutine CalcLowerPartialLoops



!------------------------------------------------------------------------------
!------  Seven different ways to start a loop  --------------------------------
!------------------------------------------------------------------------------

  Subroutine StartWWLoop
    Implicit None
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3, h, sy

    j0 = DRT(ibot)%ju(0)
    j1 = DRT(ibot)%ju(1)
    j2 = DRT(ibot)%ju(2)
    j3 = DRT(ibot)%ju(3)

    h  = DRT(ibot)%k + 1
    uu = Occ(h)
    sy = Irrep(h)

    ! These are the modified WW segment values published in
    ! A. Koslowski, M. E. Beck, W. Thiel, JCC 24, 714-726 (2003).
    If (j0 /= 0)  Call AppendWW(j0, DRT(ibot)%yu(0), 1,  Pt5* uu     * uu,      h)
    If (j1 /= 0)  Call AppendWW(j1, DRT(ibot)%yu(1), sy, Pt5* uu     *(uu-Two), h)
    If (j2 /= 0)  Call AppendWW(j2, DRT(ibot)%yu(2), sy, Pt5* uu     *(uu-Two), h)
    If (j3 /= 0)  Call AppendWW(j3, DRT(ibot)%yu(3), 1,  Pt5*(uu-Two)*(uu-Two), h)

    Return
  End Subroutine StartWWLoop



  Subroutine StartWLoop
    Implicit None
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3, h, sy

    j0 = DRT(ibot)%ju(0)
    j1 = DRT(ibot)%ju(1)
    j2 = DRT(ibot)%ju(2)
    j3 = DRT(ibot)%ju(3)

    h  = DRT(ibot)%k + 1
    uu = Occ(h)
    sy = Irrep(h)

    If (j0 /= 0)  Call AppendLowerW(j0, DRT(ibot)%yu(0), 1,     -uu, h)
    If (j1 /= 0)  Call AppendLowerW(j1, DRT(ibot)%yu(1), sy, One-uu, h)
    If (j2 /= 0)  Call AppendLowerW(j2, DRT(ibot)%yu(2), sy, One-uu, h)
    If (j3 /= 0)  Call AppendLowerW(j3, DRT(ibot)%yu(3), 1,  Two-uu, h)

    Return
  End Subroutine StartWLoop



  Subroutine StartWRbLoop
    Implicit None
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3
    Integer          :: y0, y1, y2, y3
    Integer          :: b, h, sy

    j0 = DRT(ibot)%ju(0)
    j1 = DRT(ibot)%ju(1)
    j2 = DRT(ibot)%ju(2)
    j3 = DRT(ibot)%ju(3)

    y0 = DRT(ibot)%yu(0)
    y1 = DRT(ibot)%yu(1)
    y2 = DRT(ibot)%yu(2)
    y3 = DRT(ibot)%yu(3)

    b  = DRT(ibot)%b
    h  = DRT(ibot)%k + 1
    uu = Occ(h)
    sy = Irrep(h)

    If (j1 /= 0 .And. j0 /= 0)  Call AppendUpperR(j1, j0, y1, y0, sy, 1,      -Pt5*uu,               h, h, h)  ! AppendWRb
    If (j2 /= 0 .And. j0 /= 0)  Call AppendUpperR(j2, j0, y2, y0, sy, 1,      -Pt5*uu,               h, h, h)  ! AppendWRb
    If (j3 /= 0 .And. j1 /= 0)  Call AppendUpperR(j3, j1, y3, y1, 1,  sy, (One-Pt5*uu)*FuncA(b,2,1), h, h, h)  ! AppendWRb
    If (j3 /= 0 .And. j2 /= 0)  Call AppendUpperR(j3, j2, y3, y2, 1,  sy, (One-Pt5*uu)*FuncA(b,0,1), h, h, h)  ! AppendWRb

    Return
  End Subroutine StartWRbLoop



  Subroutine StartRbRbLoop
    Implicit None
    Integer :: j0, j3, y0, y3, h

    j3 = DRT(ibot)%ju(3)
    j0 = DRT(ibot)%ju(0)

    y3 = DRT(ibot)%yu(3)
    y0 = DRT(ibot)%yu(0)

    h  = DRT(ibot)%k + 1

    If (j3 /= 0 .And. j0 /= 0)  Call AppendRR(j3, j0, y3, y0, 1, 1, Sqrt2, Zero, h, h, iLowerRbRb)  ! AppendRbRb

    Return
  End Subroutine StartRbRbLoop



  Subroutine StartRbLbLoop
    Implicit None
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3
    Integer          :: y0, y1, y2, y3
    Integer          :: b, h, sy

    j0 = DRT(ibot)%ju(0)
    j1 = DRT(ibot)%ju(1)
    j2 = DRT(ibot)%ju(2)
    j3 = DRT(ibot)%ju(3)

    y0 = DRT(ibot)%yu(0)
    y1 = DRT(ibot)%yu(1)
    y2 = DRT(ibot)%yu(2)
    y3 = DRT(ibot)%yu(3)

    b  = DRT(ibot)%b
    h  = DRT(ibot)%k + 1
    uu = Occ(h)
    sy = Irrep(h)

    If (j0 /= 0)  Call AppendRL(j0, j0, y0, y0, 1,  1,  t*      uu,             Zero, h, h, iLowerRLdg)  ! AppendRbLb
    If (j1 /= 0)  Call AppendRL(j1, j1, y1, y1, sy, sy, t*(uu-One),  t*FuncA(b, 3,1), h, h, iLowerRLdg)  ! AppendRbLb
    If (j1 /= 0 .And. j2 /= 0) &
                  Call AppendRL(j1, j2, y1, y2, sy, sy,       Zero,       FuncC(b,1), h, h, iLowerRL)    ! AppendRbLb
    If (j2 /= 0 .And. j1 /= 0) &
                  Call AppendRL(j2, j1, y2, y1, sy, sy,       Zero,       FuncC(b,1), h, h, iLowerRL)    ! AppendRbLb
    If (j2 /= 0)  Call AppendRL(j2, j2, y2, y2, sy, sy, t*(uu-One), -t*FuncA(b,-1,1), h, h, iLowerRLdg)  ! AppendRbLb
    If (j3 /= 0)  Call AppendRL(j3, j3, y3, y3, 1,  1,  t*(uu-Two),             Zero, h, h, iLowerRLdg)  ! AppendRbLb

    Return
  End Subroutine StartRbLbLoop



  Subroutine StartRbLoop
    Implicit None
    Integer :: j0, j1, j2, j3
    Integer :: y0, y1, y2, y3
    Integer :: b, h, sy

    j0 = DRT(ibot)%ju(0)
    j1 = DRT(ibot)%ju(1)
    j2 = DRT(ibot)%ju(2)
    j3 = DRT(ibot)%ju(3)

    y0 = DRT(ibot)%yu(0)
    y1 = DRT(ibot)%yu(1)
    y2 = DRT(ibot)%yu(2)
    y3 = DRT(ibot)%yu(3)

    b  = DRT(ibot)%b
    h  = DRT(ibot)%k + 1
    sy = Irrep(h)

    If (j1 /= 0 .And. j0 /= 0)  Call AppendLowerR(j1, j0, y1, y0, sy, 1,  One,          h)  ! AppendRb
    If (j2 /= 0 .And. j0 /= 0)  Call AppendLowerR(j2, j0, y2, y0, sy, 1,  One,          h)  ! AppendRb
    If (j3 /= 0 .And. j1 /= 0)  Call AppendLowerR(j3, j1, y3, y1, 1,  sy, FuncA(b,2,1), h)  ! AppendRb
    If (j3 /= 0 .And. j2 /= 0)  Call AppendLowerR(j3, j2, y3, y2, 1,  sy, FuncA(b,0,1), h)  ! AppendRb

    Return
  End Subroutine StartRbLoop



  Subroutine StartLbLoop
    Implicit None
    Integer :: j0, j1, j2, j3
    Integer :: y0, y1, y2, y3
    Integer :: b, h, sy

    j0 = DRT(ibot)%ju(0)
    j1 = DRT(ibot)%ju(1)
    j2 = DRT(ibot)%ju(2)
    j3 = DRT(ibot)%ju(3)

    y0 = DRT(ibot)%yu(0)
    y1 = DRT(ibot)%yu(1)
    y2 = DRT(ibot)%yu(2)
    y3 = DRT(ibot)%yu(3)

    b  = DRT(ibot)%b
    h  = DRT(ibot)%k + 1
    sy = Irrep(h)

    If (j0 /= 0 .And. j1 /= 0)  Call AppendL(j0, j1, y0, y1, 1,  sy, One,          h)  ! AppendLb
    If (j0 /= 0 .And. j2 /= 0)  Call AppendL(j0, j2, y0, y2, 1,  sy, One,          h)  ! AppendLb
    If (j1 /= 0 .And. j3 /= 0)  Call AppendL(j1, j3, y1, y3, sy, 1,  FuncA(b,2,1), h)  ! AppendLb
    If (j2 /= 0 .And. j3 /= 0)  Call AppendL(j2, j3, y2, y3, sy, 1,  FuncA(b,0,1), h)  ! AppendLb

    Return
  End Subroutine StartLbLoop


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the top segments  ----------------
!------------------------------------------------------------------------------

  Subroutine AppendWW(itop, y, sy, value, i)
    Implicit None
    Integer          :: itop, y, sy
    Double Precision :: value
    Integer          :: i
    ! End of dummy parameters.

    If (Dabs(value) < Tiny)  Return

    Call TwoBodyLoopComplete(itop, y, y, sy, sy, Pt5*value, i, i, i, i)

    Return
  End Subroutine AppendWW



  Subroutine AppendLowerW(itop, y, sy, value, i)
    Implicit None
    Integer          :: itop, y, sy
    Double Precision :: value
    Integer          :: i
    ! End of dummy parameters.

    If (Dabs(value) < Tiny)  Return

    Call OneBodyLoopComplete(itop, y, y, sy, sy, value, i, i)
    Call WalkUp             (itop, y, y, sy, sy, value, i, i, iLowerW)

    Return
  End Subroutine AppendLowerW


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the lower joining segments  ------
!------------------------------------------------------------------------------

  Recursive Subroutine AppendLowerR(ibra, iket, ybra, yket, sybra, syket, value, i)
    Implicit None
    Integer          :: ibra, iket, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: i
    ! End of dummy parameters.
    Double Precision :: uu
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab, prbra, prket
    Double Precision :: x

    If (Dabs(value) < Tiny)  Return

    If (DRT(iket)%k == kBorder) Then
       Call AddPartialLoop(ibra, iket, ybra, yket, sybra, syket, value, Zero, i, 0, 0, iLowerRb)
       Return
    End If

    jbra0  =     DRT(ibra)%ju(0)
    jbra1  =     DRT(ibra)%ju(1)
    jbra2  =     DRT(ibra)%ju(2)
    jbra3  =     DRT(ibra)%ju(3)

    jket0  =     DRT(iket)%ju(0)
    jket1  =     DRT(iket)%ju(1)
    jket2  =     DRT(iket)%ju(2)
    jket3  =     DRT(iket)%ju(3)

    ybra0  =     DRT(ibra)%yu(0)
    ybra1  =     DRT(ibra)%yu(1)
    ybra2  =     DRT(ibra)%yu(2)
    ybra3  =     DRT(ibra)%yu(3)

    yket0  =     DRT(iket)%yu(0)
    yket1  =     DRT(iket)%yu(1)
    yket2  =     DRT(iket)%yu(2)
    yket3  =     DRT(iket)%yu(3)

    h      =     DRT(iket)%k + 1
    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b
    uu     =     Occ(h)
    prbra  =     iProd(sybra, Irrep(h))
    prket  =     iProd(syket, Irrep(h))

    If (deltab == -1) Then

       If (jbra0 == jket1 .And. jket1 /= 0) Then
          Call AppendLowerRt(jket1, ybra+ybra0, yket+yket1, sybra, prket,                           value, i, h)
          Call AppendTop    (jket1, ybra+ybra0, yket+yket1, sybra, prket,                   -Pt5*uu*value, i, h, h, h)  ! AppendWRt
       End If

       If (jbra2 == jket3 .And. jket3 /= 0) Then
          Call AppendLowerRt(jket3, ybra+ybra2, yket+yket3, prbra, syket,              FuncA(b,2,1)*value, i, h)
          Call AppendTop    (jket3, ybra+ybra2, yket+yket3, prbra, syket, (One-Pt5*uu)*FuncA(b,2,1)*value, i, h, h, h)  ! AppendWRt
       End If

       If (jbra0 /= 0 .And. jket0 /= 0) Then
          Call AppendLowerR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,                     value, i)
          Call AppendUpperR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,                 -uu*value, h, h, i)  ! AppendRW
          Call AppendUpperR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,              Pt5*uu*value, i, h, h)  ! AppendRbRt
       End If

       If (jbra1 /= 0 .And. jket1 /= 0) Then
          Call AppendLowerR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket,                    -value, i)
          Call AppendUpperR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket,            (uu-One)*value, h, h, i)  ! AppendRW
          Call AppendUpperR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket,        (One-Pt5*uu)*value, i, h, h)  ! AppendRbRt
       End If

       If (jbra2 /= 0 .And. jket1 /= 0) Then
          x = One / Dble(b+1)
          Call AppendLowerR(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket,                   x*value, i)
          Call AppendUpperR(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket,          (One-uu)*x*value, h, h, i)  ! AppendRW
          Call AppendUpperR(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket,      (One+Pt5*uu*x)*value, i, h, h)  ! AppendRbRt
       End If

       If (jbra2 /= 0 .And. jket2 /= 0) Then
          Call AppendLowerR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket,          FuncC(b,1)*value, i)
          Call AppendUpperR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, (One-uu)*FuncC(b,1)*value, h, h, i)  ! AppendRW
          Call AppendUpperR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket,   Pt5*uu*FuncC(b,1)*value, i, h, h)  ! AppendRbRt
       End If

       If (jbra3 /= 0 .And. jket3 /= 0) Then
          Call AppendLowerR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,                    -value, i)
          Call AppendUpperR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,            (uu-Two)*value, h, h, i)  ! AppendRW
          Call AppendUpperR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,        (One-Pt5*uu)*value, i, h, h)  ! AppendRbRt
       End If

       If (jbra0 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra0, jket1, ybra+ybra0, yket+yket1, sybra, prket,              -t*value,  t*FuncA(b,3,1)*value, h, i, iLowerRL)  ! AppendRLb
       If (jbra0 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra0, jket2, ybra+ybra0, yket+yket2, sybra, prket,                  Zero,      FuncC(b,1)*value, h, i, iLowerRL)  ! AppendRLb
       If (jbra1 /= 0 .And. jket3 /= 0) &
          Call AppendRL(jbra1, jket3, ybra+ybra1, yket+yket3, prbra, syket,                  Zero,   -FuncA(b,3,2)*value, h, i, iLowerRL)  ! AppendRLb
       If (jbra2 /= 0 .And. jket3 /= 0) &
          Call AppendRL(jbra2, jket3, ybra+ybra2, yket+yket3, prbra, syket, -t*FuncA(b,2,1)*value,  t*FuncA(b,0,1)*value, h, i, iLowerRL)  ! AppendRLb

       If (jbra1 /= 0 .And. jket0 /= 0) Then
          Call AppendRR(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket,                  Zero,                 value, i, h, iLowerRbR)  ! AppendRbR
          Call AppendRR(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket,                  Zero,                -value, h, i, iLowerRRb)  ! AppendRRb
       End If

       If (jbra2 /= 0 .And. jket0 /= 0) Then
          Call AppendRR(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket,  t*FuncA(b,2,1)*value, -t*FuncA(b,0,1)*value, i, h, iLowerRbR)  ! AppendRbR
          Call AppendRR(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket,  t*FuncA(b,2,1)*value,  t*FuncA(b,0,1)*value, h, i, iLowerRRb)  ! AppendRRb
       End If

       If (jbra3 /= 0 .And. jket1 /= 0) Then
          Call AppendRR(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket,              -t*value,  t*FuncA(b,3,1)*value, i, h, iLowerRbR)  ! AppendRbR
          Call AppendRR(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket,              -t*value, -t*FuncA(b,3,1)*value, h, i, iLowerRRb)  ! AppendRRb
       End If

       If (jbra3 /= 0 .And. jket2 /= 0) Then
          Call AppendRR(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket,                  Zero,    FuncA(b,0,1)*value, i, h, iLowerRbR)  ! AppendRbR
          Call AppendRR(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket,                  Zero,   -FuncA(b,0,1)*value, h, i, iLowerRRb)  ! AppendRRb
       End If

    Else If (deltab == 1) Then

       If (jbra0 == jket2 .And. jket2 /= 0) Then
          Call AppendLowerRt(jket2, ybra+ybra0, yket+yket2, sybra, prket,                           value, i, h)
          Call AppendTop    (jket2, ybra+ybra0, yket+yket2, sybra, prket,                   -Pt5*uu*value, i, h, h, h)  ! AppendWRt
       End If

       If (jbra1 == jket3 .And. jket3 /= 0) Then
          Call AppendLowerRt(jket3, ybra+ybra1, yket+yket3, prbra, syket,              FuncA(b,0,1)*value, i, h)
          Call AppendTop    (jket3, ybra+ybra1, yket+yket3, prbra, syket, (One-Pt5*uu)*FuncA(b,0,1)*value, i, h, h, h)  ! AppendWRt
       End If

       If (jbra0 /= 0 .And. jket0 /= 0) Then
          Call AppendLowerR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,                     value, i)
          Call AppendUpperR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,                 -uu*value, h, h, i)  ! AppendRW
          Call AppendUpperR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,              Pt5*uu*value, i, h, h)  ! AppendRbRt
       End If

       If (jbra1 /= 0 .And. jket1 /= 0) Then
          Call AppendLowerR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket,          FuncC(b,1)*value, i)
          Call AppendUpperR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, (One-uu)*FuncC(b,1)*value, h, h, i)  ! AppendRW
          Call AppendUpperR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket,   Pt5*uu*FuncC(b,1)*value, i, h, h)  ! AppendRbRt
       End If

       If (jbra1 /= 0 .And. jket2 /= 0) Then
          x = One / Dble(b+1)
          Call AppendLowerR(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket,                  -x*value, i)
          Call AppendUpperR(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket,          (uu-One)*x*value, h, h, i)  ! AppendRW
          Call AppendUpperR(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket,      (One-Pt5*uu*x)*value, i, h, h)  ! AppendRbRt
       End If

       If (jbra2 /= 0 .And. jket2 /= 0) Then
          Call AppendLowerR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket,                    -value, i)
          Call AppendUpperR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket,            (uu-One)*value, h, h, i)  ! AppendRW
          Call AppendUpperR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket,        (One-Pt5*uu)*value, i, h, h)  ! AppendRbRt
       End If

       If (jbra3 /= 0 .And. jket3 /= 0) Then
          Call AppendLowerR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,                    -value, i)
          Call AppendUpperR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,            (uu-Two)*value, h, h, i)  ! AppendRW
          Call AppendUpperR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,        (One-Pt5*uu)*value, i, h, h)  ! AppendRbRt
       End If

       If (jbra0 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra0, jket1, ybra+ybra0, yket+yket1, sybra, prket,                  Zero,       FuncC(b,1)*value, h, i, iLowerRL)  ! AppendRLb
       If (jbra0 /= 0 .And. jket2 /= 0) &                                                                                    
          Call AppendRL(jbra0, jket2, ybra+ybra0, yket+yket2, sybra, prket,              -t*value, -t*FuncA(b,-1,1)*value, h, i, iLowerRL)  ! AppendRLb
       If (jbra1 /= 0 .And. jket3 /= 0) &                                                                                    
          Call AppendRL(jbra1, jket3, ybra+ybra1, yket+yket3, prbra, syket, -t*FuncA(b,0,1)*value, -t*FuncA(b, 2,1)*value, h, i, iLowerRL)  ! AppendRLb
       If (jbra2 /= 0 .And. jket3 /= 0) &                                                                                    
          Call AppendRL(jbra2, jket3, ybra+ybra2, yket+yket3, prbra, syket,                  Zero,   -FuncA(b,-1,0)*value, h, i, iLowerRL)  ! AppendRLb

       If (jbra1 /= 0 .And. jket0 /= 0) Then
          Call AppendRR(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket, t*FuncA(b,0,1)*value,  t*FuncA(b, 2,1)*value, i, h, iLowerRbR)  ! AppendRbR
          Call AppendRR(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket, t*FuncA(b,0,1)*value, -t*FuncA(b, 2,1)*value, h, i, iLowerRRb)  ! AppendRRb
       End If

       If (jbra2 /= 0 .And. jket0 /= 0) Then
          Call AppendRR(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket,                 Zero,                  value, i, h, iLowerRbR)  ! AppendRbR
          Call AppendRR(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket,                 Zero,                 -value, h, i, iLowerRRb)  ! AppendRRb
       End If

       If (jbra3 /= 0 .And. jket1 /= 0) Then
          Call AppendRR(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket,                 Zero,    FuncA(b, 2,1)*value, i, h, iLowerRbR)  ! AppendRbR
          Call AppendRR(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket,                 Zero,   -FuncA(b, 2,1)*value, h, i, iLowerRRb)  ! AppendRRb
       End If

       If (jbra3 /= 0 .And. jket2 /= 0) Then
          Call AppendRR(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket,             -t*value, -t*FuncA(b,-1,1)*value, i, h, iLowerRbR)  ! AppendRbR
          Call AppendRR(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket,             -t*value,  t*FuncA(b,-1,1)*value, h, i, iLowerRRb)  ! AppendRRb
       End If

    End If

    Return
  End Subroutine AppendLowerR



  Recursive Subroutine AppendL(ibra, iket, ybra, yket, sybra, syket, value, i)
    Implicit None
    Integer          :: ibra, iket, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: i, j, k
    ! End of dummy parameters.
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab, prbra, prket

    If (Dabs(value) < Tiny)  Return

    If (DRT(iket)%k == kBorder) Then
       Call AddPartialLoop(ibra, iket, ybra, yket, sybra, syket, value, Zero, i, 0, 0, iLowerLb)
       Return
    End If

    jbra0  =     DRT(ibra)%ju(0)
    jbra1  =     DRT(ibra)%ju(1)
    jbra2  =     DRT(ibra)%ju(2)
    jbra3  =     DRT(ibra)%ju(3)

    jket0  =     DRT(iket)%ju(0)
    jket1  =     DRT(iket)%ju(1)
    jket2  =     DRT(iket)%ju(2)
    jket3  =     DRT(iket)%ju(3)

    ybra0  =     DRT(ibra)%yu(0)
    ybra1  =     DRT(ibra)%yu(1)
    ybra2  =     DRT(ibra)%yu(2)
    ybra3  =     DRT(ibra)%yu(3)

    yket0  =     DRT(iket)%yu(0)
    yket1  =     DRT(iket)%yu(1)
    yket2  =     DRT(iket)%yu(2)
    yket3  =     DRT(iket)%yu(3)

    h      =     DRT(iket)%k + 1
    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b
    prbra  =     iProd(sybra, Irrep(h))
    prket  =     iProd(syket, Irrep(h))

    If (deltab == -1) Then

       If (jbra2 == jket0 .And. jket0 /= 0) &
          Call AppendLowerLt(jket0, ybra+ybra2, yket+yket0, prbra, syket,              value, i, h)
       If (jbra3 == jket1 .And. jket1 /= 0) &
          Call AppendLowerLt(jket1, ybra+ybra3, yket+yket1, sybra, prket, FuncA(b,1,2)*value, i, h)

       If (jbra3 /= 0 .And. jket0 /= 0) &
          Call AppendUpperR(jbra3, jket0, ybra+ybra3, yket+yket0, sybra, syket, FuncA(b,1,2)*value, i, h, h)  ! AppendRbLt

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendL(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,            value,           i)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendL(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, FuncC(b,2)*value,           i)
       If (jbra2 /= 0 .And. jket1 /= 0) &
          Call AppendL(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket,           -value/Dble(b+2), i)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendL(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket,           -value,           i)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendL(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,           -value,           i)

       If (jbra1 /= 0 .And. jket0 /= 0) &
          Call AppendRL(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket,                  Zero,      FuncC(b,2)*value, i, h, iLowerRL)  ! AppendRbL
       If (jbra2 /= 0 .And. jket0 /= 0) &
          Call AppendRL(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket,              -t*value, -t*FuncA(b,0,2)*value, i, h, iLowerRL)  ! AppendRbL
       If (jbra3 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket, -t*FuncA(b,1,2)*value, -t*FuncA(b,3,2)*value, i, h, iLowerRL)  ! AppendRbL
       If (jbra3 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket,                  Zero,   -FuncA(b,0,1)*value, i, h, iLowerRL)  ! AppendRbL

    Else If (deltab == 1) Then

       If (jbra1 == jket0 .And. jket0 /= 0) &
          Call AppendLowerLt(jket0, ybra+ybra1, yket+yket0, prbra, syket,              value, i, h)
       If (jbra3 == jket2 .And. jket2 /= 0) &
          Call AppendLowerLt(jket2, ybra+ybra3, yket+yket2, sybra, prket, FuncA(b,1,0)*value, i, h)

       If (jbra3 /= 0 .And. jket0 /= 0) &
          Call AppendUpperR(jbra3, jket0, ybra+ybra3, yket+yket0, sybra, syket, FuncA(b,1,0)*value, i, h, h)  ! AppendRbLt

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendL(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,            value,         i)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendL(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket,           -value,         i)
       If (jbra1 /= 0 .And. jket2 /= 0) &
          Call AppendL(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket,            value/Dble(b), i)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendL(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, FuncC(b,0)*value,         i)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendL(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,           -value,         i)

       If (jbra1 /= 0 .And. jket0 /= 0) &
          Call AppendRL(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket,              -t*value, t*FuncA(b, 2,0)*value, i, h, iLowerRL)  ! AppendRbL
       If (jbra2 /= 0 .And. jket0 /= 0) &
          Call AppendRL(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket,                  Zero,      FuncC(b,0)*value, i, h, iLowerRL)  ! AppendRbL
       If (jbra3 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket,                  Zero,  -FuncA(b, 2,1)*value, i, h, iLowerRL)  ! AppendRbL
       If (jbra3 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket, -t*FuncA(b,1,0)*value, t*FuncA(b,-1,0)*value, i, h, iLowerRL)  ! AppendRbL

    End If

    Return
  End Subroutine AppendL


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the lower middle segments  -------
!------------------------------------------------------------------------------

  Subroutine AppendLowerRt(itop, ybra, yket, sybra, syket, value, i, j)
    Implicit None
    Integer          :: itop, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: i, j
    ! End of dummy parameters.

    If (Dabs(value) < Tiny) Return

    Call OneBodyLoopComplete(itop, ybra, yket, sybra, syket, value, i, j)
    Call WalkUp             (itop, ybra, yket, sybra, syket, value, i, j, iLowerOneR)

    Return
  End Subroutine AppendLowerRt



  Subroutine AppendLowerLt(itop, ybra, yket, sybra, syket, value, i, j)
    Implicit None
    Integer          :: itop, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: i, j
    ! End of dummy parameters.

    If (Dabs(value) < Tiny) Return

    Call WalkUp(itop, ybra, yket, sybra, syket, value, i, j, iLowerOneL)

    Return
  End Subroutine AppendLowerLt


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the middle joining segments  -----
!------------------------------------------------------------------------------

  Recursive Subroutine WalkUp(irow, ybra, yket, sybra, syket, value, i, j, itype)
    Implicit None
    Integer          :: irow, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: i, j, itype
    ! End of dummy parameters.
    Double Precision :: uu
    Integer          :: j0, j1, j2, j3, y0, y1, y2, y3, h, b, prbra, prket

    If (DRT(irow)%k == kBorder) Then
       Call AddPartialLoop(irow, irow, ybra, yket, sybra, syket, value, Zero, i, j, 0, itype)
       Return
    End If

    j0    = DRT(irow)%ju(0)
    j1    = DRT(irow)%ju(1)
    j2    = DRT(irow)%ju(2)
    j3    = DRT(irow)%ju(3)

    y0    = DRT(irow)%yu(0)
    y1    = DRT(irow)%yu(1)
    y2    = DRT(irow)%yu(2)
    y3    = DRT(irow)%yu(3)

    h     = DRT(irow)%k + 1
    b     = DRT(irow)%b
    uu    = Occ(h)
    prbra = iProd(sybra, Irrep(h))
    prket = iProd(syket, Irrep(h))

    ! A W segment must not be attached above a LtLb one-body loop, otherwise the resulting
    ! two-body loop will correspond to the lower triangle of the CI Hamiltonian.
    If (itype /= iLowerOneL) Then
       If (j0 /= 0)  Call AppendTop(j0, ybra+y0, yket+y0, sybra, syket,     -uu *value, i, j, h, h)  ! AppendUpperW
       If (j1 /= 0)  Call AppendTop(j1, ybra+y1, yket+y1, prbra, prket, (One-uu)*value, i, j, h, h)  ! AppendUpperW
       If (j2 /= 0)  Call AppendTop(j2, ybra+y2, yket+y2, prbra, prket, (One-uu)*value, i, j, h, h)  ! AppendUpperW
       If (j3 /= 0)  Call AppendTop(j3, ybra+y3, yket+y3, sybra, syket, (Two-uu)*value, i, j, h, h)  ! AppendUpperW
    End If

    If (j0 /= 0)  Call WalkUp(j0, ybra+y0, yket+y0, sybra, syket, value, i, j, itype)
    If (j1 /= 0)  Call WalkUp(j1, ybra+y1, yket+y1, prbra, prket, value, i, j, itype)
    If (j2 /= 0)  Call WalkUp(j2, ybra+y2, yket+y2, prbra, prket, value, i, j, itype)
    If (j3 /= 0)  Call WalkUp(j3, ybra+y3, yket+y3, sybra, syket, value, i, j, itype)

    If (j1 /= 0 .And. j0 /= 0)  Call AppendUpperR(j1, j0, ybra+y1, yket+y0, prbra, syket,              value, i, j, h)  ! AppendUpperRb
    If (j2 /= 0 .And. j0 /= 0)  Call AppendUpperR(j2, j0, ybra+y2, yket+y0, prbra, syket,              value, i, j, h)  ! AppendUpperRb
    If (j3 /= 0 .And. j1 /= 0)  Call AppendUpperR(j3, j1, ybra+y3, yket+y1, sybra, prket, FuncA(b,2,1)*value, i, j, h)  ! AppendUpperRb
    If (j3 /= 0 .And. j2 /= 0)  Call AppendUpperR(j3, j2, ybra+y3, yket+y2, sybra, prket, FuncA(b,0,1)*value, i, j, h)  ! AppendUpperRb

    Return
  End Subroutine WalkUp



  Recursive Subroutine AppendRL(ibra, iket, ybra, yket, sybra, syket, sval, tval, i, j, itype)
    Implicit None
    Integer          :: ibra, iket, ybra, yket, sybra, syket
    Double Precision :: sval, tval
    Integer          :: i, j, itype
    ! End of dummy parameters.
    Double Precision :: uu
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab, prbra, prket, ityp

    If (Dabs(sval) < Tiny .And. Dabs(tval) < Tiny)  Return

    If (DRT(iket)%k == kBorder) Then
       ityp = itype
       If (itype == iLowerRL .And. ibra == iket) Then
          If (ybra < yket) Then
             ityp = iLowerRLcp
          Else
             ityp = iLowerRLcm
          End If
       End If
       Call AddPartialLoop(ibra, iket, ybra, yket, sybra, syket, sval, tval, i, j, 0, ityp)
       Return
    End If

    jbra0  =     DRT(ibra)%ju(0)
    jbra1  =     DRT(ibra)%ju(1)
    jbra2  =     DRT(ibra)%ju(2)
    jbra3  =     DRT(ibra)%ju(3)

    jket0  =     DRT(iket)%ju(0)
    jket1  =     DRT(iket)%ju(1)
    jket2  =     DRT(iket)%ju(2)
    jket3  =     DRT(iket)%ju(3)

    ybra0  =     DRT(ibra)%yu(0)
    ybra1  =     DRT(ibra)%yu(1)
    ybra2  =     DRT(ibra)%yu(2)
    ybra3  =     DRT(ibra)%yu(3)

    yket0  =     DRT(iket)%yu(0)
    yket1  =     DRT(iket)%yu(1)
    yket2  =     DRT(iket)%yu(2)
    yket3  =     DRT(iket)%yu(3)

    h      =     DRT(iket)%k + 1
    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b
    uu     =     Occ(h)
    prbra  =     iProd(sybra, Irrep(h))
    prket  =     iProd(syket, Irrep(h))

    If (deltab == -2) Then

       ! An RtLt segment with d'd=21 will always lead to a generator
       ! matrix element corresponding to the lower triangle of the
       ! CI Hamiltonian and must never be attached.

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendRL(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket, Zero,               tval, i, j, itype)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, Zero,   -FuncC(b,3)*tval, i, j, itype)
       If (jbra2 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket, Zero, -FuncB(b,1,3)*tval, i, j, iLowerRL)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, Zero,   -FuncC(b,1)*tval, i, j, itype)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendRL(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket, Zero,               tval, i, j, itype)

       If (jbra2 /= 0 .And. jket0 /= 0) &
          Call AppendUpperR(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket,               tval, i, h, j)  ! AppendRLt
       If (jbra3 /= 0 .And. jket1 /= 0) &
          Call AppendUpperR(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket, -FuncA(b,2,3)*tval, i, h, j)  ! AppendRLt

    Else If (deltab == 0) Then

       ! An RtLt segment with d'=d must only be attached if the resulting
       ! loop corresponds to the upper triangle of the CI Hamiltonian.
       If (jbra0 == jket0 .And. jket0 /= 0 .And. ybra <= yket) &
          Call AppendTop(jket0, ybra+ybra0, yket+yket0, sybra, syket, -t*      uu *sval,                    i, h, j, h)  ! AppendRtLt
       If (jbra1 == jket1 .And. jket1 /= 0 .And. ybra <= yket) &
          Call AppendTop(jket1, ybra+ybra1, yket+yket1, prbra, prket,  t*((One-uu)*sval-FuncA(b,0,2)*tval), i, h, j, h)  ! AppendRtLt
       If (jbra2 == jket2 .And. jket2 /= 0 .And. ybra <= yket) &
          Call AppendTop(jket2, ybra+ybra2, yket+yket2, prbra, prket,  t*((One-uu)*sval+FuncA(b,2,0)*tval), i, h, j, h)  ! AppendRtLt
       If (jbra3 == jket3 .And. jket3 /= 0 .And. ybra <= yket) &
          Call AppendTop(jket3, ybra+ybra3, yket+yket3, sybra, syket,  t* (Two-uu)*sval,                    i, h, j, h)  ! AppendRtLt

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendRL(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket, sval,                  tval, i, j, itype)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, sval,       FuncD(b,1)*tval, i, j, itype)
       If (jbra1 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket, Zero, -Sqrt2*tval/Dble(b+1), i, j, iLowerRL)
       If (jbra2 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket, Zero, -Sqrt2*tval/Dble(b+1), i, j, iLowerRL)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, sval,       FuncD(b,0)*tval, i, j, itype)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendRL(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket, sval,                  tval, i, j, itype)

       If (jbra1 /= 0 .And. jket0 /= 0) &
          Call AppendUpperR(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket,              t*(sval-FuncA(b,0,2)*tval), i, h, j)  ! AppendRLt
       If (jbra2 /= 0 .And. jket0 /= 0) &                                     
          Call AppendUpperR(jbra2, jket0, ybra+ybra2, yket+yket0, prbra, syket,              t*(sval+FuncA(b,2,0)*tval), i, h, j)  ! AppendRLt
       If (jbra3 /= 0 .And. jket1 /= 0) &
          Call AppendUpperR(jbra3, jket1, ybra+ybra3, yket+yket1, sybra, prket, t*(FuncA(b,2,1)*sval+FuncA(b,0,1)*tval), i, h, j)  ! AppendRLt
       If (jbra3 /= 0 .And. jket2 /= 0) &
          Call AppendUpperR(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket, t*(FuncA(b,0,1)*sval-FuncA(b,2,1)*tval), i, h, j)  ! AppendRLt

    Else If (deltab == 2) Then

       ! An RtLt segment with d'd=12 will always lead to a generator matrix
       ! element corresponding to the upper triangle of the CI Hamiltonian.
       If (jbra1 == jket2 .And. jket2 /= 0) &
          Call AppendTop(jket2, ybra+ybra1, yket+yket2, prbra, prket, tval, i, h, j, h)  ! AppendRtLt

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendRL(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket, Zero,                tval, i, j, itype)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendRL(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, Zero,   -FuncC(b, 1)*tval, i, j, itype)
       If (jbra1 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket, Zero, -FuncB(b,-1,1)*tval, i, j, iLowerRL)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendRL(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, Zero,   -FuncC(b,-1)*tval, i, j, itype)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendRL(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket, Zero,                tval, i, j, itype)

       If (jbra1 /= 0 .And. jket0 /= 0) &
          Call AppendUpperR(jbra1, jket0, ybra+ybra1, yket+yket0, prbra, syket,                tval, i, h, j)  ! AppendRLt
       If (jbra3 /= 0 .And. jket2 /= 0) &
          Call AppendUpperR(jbra3, jket2, ybra+ybra3, yket+yket2, sybra, prket, -FuncA(b,0,-1)*tval, i, h, j)  ! AppendRLt

    End If

    Return
  End Subroutine AppendRL



  Recursive Subroutine AppendRR(ibra, iket, ybra, yket, sybra, syket, sval, tval, i, j, itype)
    Implicit None
    Integer          :: ibra, iket, ybra, yket, sybra, syket
    Double Precision :: sval, tval
    Integer          :: i, j, itype
    ! End of dummy parameters.
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab, prbra, prket

    If (Dabs(sval) < Tiny .And. Dabs(tval) < Tiny)  Return

    If (DRT(iket)%k == kBorder) Then
       Call AddPartialLoop(ibra, iket, ybra, yket, sybra, syket, sval, tval, i, j, 0, itype)
       Return
    End If

    jbra0  =     DRT(ibra)%ju(0)
    jbra1  =     DRT(ibra)%ju(1)
    jbra2  =     DRT(ibra)%ju(2)
    jbra3  =     DRT(ibra)%ju(3)

    jket0  =     DRT(iket)%ju(0)
    jket1  =     DRT(iket)%ju(1)
    jket2  =     DRT(iket)%ju(2)
    jket3  =     DRT(iket)%ju(3)

    ybra0  =     DRT(ibra)%yu(0)
    ybra1  =     DRT(ibra)%yu(1)
    ybra2  =     DRT(ibra)%yu(2)
    ybra3  =     DRT(ibra)%yu(3)

    yket0  =     DRT(iket)%yu(0)
    yket1  =     DRT(iket)%yu(1)
    yket2  =     DRT(iket)%yu(2)
    yket3  =     DRT(iket)%yu(3)

    h      =     DRT(iket)%k + 1
    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b
    prbra  =     iProd(sybra, Irrep(h))
    prket  =     iProd(syket, Irrep(h))

    If (deltab == -2) Then

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendRR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket, Zero,              tval, i, j, itype)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendRR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, Zero,              tval, i, j, itype)
       If (jbra2 /= 0 .And. jket1 /= 0) &
          Call AppendRR(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket, Zero, FuncB(b,1,2)*tval, i, j, itype)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendRR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, Zero,   FuncD(b,1)*tval, i, j, itype)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendRR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket, Zero,              tval, i, j, itype)

       If (jbra0 /= 0 .And. jket1 /= 0) &
          Call AppendUpperR(jbra0, jket1, ybra+ybra0, yket+yket1, sybra, prket,               tval, i, h, j)  ! AppendRRt
       If (jbra2 /= 0 .And. jket3 /= 0) &
          Call AppendUpperR(jbra2, jket3, ybra+ybra2, yket+yket3, prbra, syket,  FuncA(b,3,2)*tval, i, h, j)  ! AppendRRt

    Else If (deltab == 0) Then

       ! An RtRt segment must only be attached above an RbR or RbRb,
       ! not an RRb segment, otherwise the resulting loop will
       ! correspond to the lower triangle of the CI Hamiltonian.
       If (jbra0 == jket3 .And. jket3 /= 0) Then
          If (itype == iLowerRbR) Then
             Call AppendTop(jket3, ybra+ybra0, yket+yket3, sybra, syket, Sqrt2*sval, i, h, j, h)  ! AppendRtRt
          Else If (itype == iLowerRbRb) Then
             ! t is the product of the RtRt singlet coupling value of Sqrt2 with Pt5.
             Call AppendTop(jket3, ybra+ybra0, yket+yket3, sybra, syket,     t*sval, i, h, j, h)  ! AppendRtRt
          End If
       End If

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendRR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,  sval,              tval, i, j, itype)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendRR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, -sval,  -FuncD(b,1)*tval, i, j, itype)
       If (jbra1 /= 0 .And. jket2 /= 0) &
          Call AppendRR(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket,  Zero, FuncB(b,1,2)*tval, i, j, itype)
       If (jbra2 /= 0 .And. jket1 /= 0) &
          Call AppendRR(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket,  Zero, FuncB(b,0,1)*tval, i, j, itype)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendRR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, -sval,  -FuncD(b,0)*tval, i, j, itype)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendRR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,  sval,              tval, i, j, itype)

       If (jbra0 /= 0 .And. jket1 /= 0) &
          Call AppendUpperR(jbra0, jket1, ybra+ybra0, yket+yket1, sybra, prket, t*(FuncA(b,2,1)*sval-FuncA(b,0,1)*tval), i, h, j)  ! AppendRRt
       If (jbra0 /= 0 .And. jket2 /= 0) &
          Call AppendUpperR(jbra0, jket2, ybra+ybra0, yket+yket2, sybra, prket, t*(FuncA(b,0,1)*sval+FuncA(b,2,1)*tval), i, h, j)  ! AppendRRt
       If (jbra1 /= 0 .And. jket3 /= 0) &
          Call AppendUpperR(jbra1, jket3, ybra+ybra1, yket+yket3, prbra, syket,             -t*(sval+FuncA(b,0,2)*tval), i, h, j)  ! AppendRRt
       If (jbra2 /= 0 .And. jket3 /= 0) &
          Call AppendUpperR(jbra2, jket3, ybra+ybra2, yket+yket3, prbra, syket,             -t*(sval-FuncA(b,2,0)*tval), i, h, j)  ! AppendRRt

    Else If (deltab == 2) Then

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendRR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket, Zero,              tval, i, j, itype)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendRR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, Zero,   FuncD(b,0)*tval, i, j, itype)
       If (jbra1 /= 0 .And. jket2 /= 0) &
          Call AppendRR(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket, Zero, FuncB(b,0,1)*tval, i, j, itype)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendRR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, Zero,              tval, i, j, itype)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendRR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket, Zero,              tval, i, j, itype)

       If (jbra0 /= 0 .And. jket2 /= 0) &
          Call AppendUpperR(jbra0, jket2, ybra+ybra0, yket+yket2, sybra, prket,                tval, i, h, j)  ! AppendRRt
       If (jbra1 /= 0 .And. jket3 /= 0) &
          Call AppendUpperR(jbra1, jket3, ybra+ybra1, yket+yket3, prbra, syket,  FuncA(b,-1,0)*tval, i, h, j)  ! AppendRRt

    End If

    Return
  End Subroutine AppendRR


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the upper middle segments  -------
!------------------------------------------------------------------------------

  ! These subroutines would contain dummy calls only and have been removed.


!------------------------------------------------------------------------------
!------  Subroutines for the construction of the upper joining segments  ------
!------------------------------------------------------------------------------

  Recursive Subroutine AppendUpperR(ibra, iket, ybra, yket, sybra, syket, value, i, j, k)
    Implicit None
    Integer          :: ibra, iket, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: i, j, k
    ! End of dummy parameters.
    Integer          :: jbra0, jbra1, jbra2, jbra3
    Integer          :: jket0, jket1, jket2, jket3
    Integer          :: ybra0, ybra1, ybra2, ybra3
    Integer          :: yket0, yket1, yket2, yket3
    Integer          :: h, b, deltab, prbra, prket

    If (Dabs(value) < Tiny) Return

    If (DRT(iket)%k == kBorder) Then
       Call AddPartialLoop(ibra, iket, ybra, yket, sybra, syket, value, Zero, i, j, k, iLowerR)
       Return
    End If

    jbra0  =     DRT(ibra)%ju(0)
    jbra1  =     DRT(ibra)%ju(1)
    jbra2  =     DRT(ibra)%ju(2)
    jbra3  =     DRT(ibra)%ju(3)

    jket0  =     DRT(iket)%ju(0)
    jket1  =     DRT(iket)%ju(1)
    jket2  =     DRT(iket)%ju(2)
    jket3  =     DRT(iket)%ju(3)

    ybra0  =     DRT(ibra)%yu(0)
    ybra1  =     DRT(ibra)%yu(1)
    ybra2  =     DRT(ibra)%yu(2)
    ybra3  =     DRT(ibra)%yu(3)

    yket0  =     DRT(iket)%yu(0)
    yket1  =     DRT(iket)%yu(1)
    yket2  =     DRT(iket)%yu(2)
    yket3  =     DRT(iket)%yu(3)

    h      =     DRT(iket)%k + 1
    b      =     DRT(iket)%b
    deltab = b - DRT(ibra)%b
    prbra  =     iProd(sybra, Irrep(h))
    prket  =     iProd(syket, Irrep(h))

    If (deltab == -1) Then

       If (jbra0 == jket1 .And. jket1 /= 0) &
          Call AppendTop(jket1, ybra+ybra0, yket+yket1, sybra, prket,              value, i, j, k, h)  ! AppendRt
       If (jbra2 == jket3 .And. jket3 /= 0) &
          Call AppendTop(jket3, ybra+ybra2, yket+yket3, prbra, syket, FuncA(b,2,1)*value, i, j, k, h)  ! AppendRt

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendUpperR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,            value,           i, j, k)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendUpperR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket,           -value,           i, j, k)
       If (jbra2 /= 0 .And. jket1 /= 0) &
          Call AppendUpperR(jbra2, jket1, ybra+ybra2, yket+yket1, prbra, prket,            value/Dble(b+1), i, j, k)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendUpperR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket, FuncC(b,1)*value,           i, j, k)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendUpperR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,           -value,           i, j, k)

    Else If (deltab == 1) Then

       If (jbra0 == jket2 .And. jket2 /= 0) &
          Call AppendTop(jket2, ybra+ybra0, yket+yket2, sybra, prket,              value, i, j, k, h)  ! AppendRt
       If (jbra1 == jket3 .And. jket3 /= 0) &
          Call AppendTop(jket3, ybra+ybra1, yket+yket3, prbra, syket, FuncA(b,0,1)*value, i, j, k, h)  ! AppendRt

       If (jbra0 /= 0 .And. jket0 /= 0) &
          Call AppendUpperR(jbra0, jket0, ybra+ybra0, yket+yket0, sybra, syket,            value,           i, j, k)
       If (jbra1 /= 0 .And. jket1 /= 0) &
          Call AppendUpperR(jbra1, jket1, ybra+ybra1, yket+yket1, prbra, prket, FuncC(b,1)*value,           i, j, k)
       If (jbra1 /= 0 .And. jket2 /= 0) &
          Call AppendUpperR(jbra1, jket2, ybra+ybra1, yket+yket2, prbra, prket,           -value/Dble(b+1), i, j, k)
       If (jbra2 /= 0 .And. jket2 /= 0) &
          Call AppendUpperR(jbra2, jket2, ybra+ybra2, yket+yket2, prbra, prket,           -value,           i, j, k)
       If (jbra3 /= 0 .And. jket3 /= 0) &
          Call AppendUpperR(jbra3, jket3, ybra+ybra3, yket+yket3, sybra, syket,           -value,           i, j, k)

    End If

    Return
  End Subroutine AppendUpperR



!------------------------------------------------------------------------------
!------  Subroutines for the construction of the top segments  ----------------
!------------------------------------------------------------------------------

  Subroutine AppendTop(itop, ybra, yket, sybra, syket, value, i, j, k, l)
    Implicit None
    Integer          :: itop, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: i, j, k, l
    ! End of dummy parameters.

    If (Dabs(value) < Tiny)  Return

    Call TwoBodyLoopComplete(itop, ybra, yket, sybra, syket, value, i, j, k, l)

    Return
  End Subroutine AppendTop


!------------------------------------------------------------------------------
!------  Subroutines to store partial and complete loops into a buffer  -------
!------------------------------------------------------------------------------

  Subroutine AddPartialLoop(ibra, iket, ybra, yket, sybra, syket, sval, tval, i, j, k, itype)
    Implicit None
    Integer          :: ibra, iket, ybra, yket, sybra, syket
    Double Precision :: sval, tval
    Integer          :: i, j, k, itype
    ! End of dummy parameters.
    Type(PartialLoopBuffer), Pointer :: newbuf

    If (.Not. Associated(lists%firstL)) Then
       Allocate(lists%firstL)
       Nullify(lists%firstL%next)
       lists%firstL%ncount = 1
    Else If (lists%firstL%ncount >= loopBufSize) Then
       Allocate(newbuf)
       newbuf%next   => lists%firstL
       lists%firstL  => newbuf
       newbuf%ncount =  1
    Else
       lists%firstL%ncount = lists%firstL%ncount + 1
    End If

    lists%firstL%part(lists%firstL%ncount)%singlet = sval
    lists%firstL%part(lists%firstL%ncount)%triplet = tval
    lists%firstL%part(lists%firstL%ncount)%itype   = itype
    lists%firstL%part(lists%firstL%ncount)%itip    = ibot
    lists%firstL%part(lists%firstL%ncount)%ibra    = ibra
    lists%firstL%part(lists%firstL%ncount)%iket    = iket
    lists%firstL%part(lists%firstL%ncount)%lexbra  = ybra
    lists%firstL%part(lists%firstL%ncount)%lexket  = yket
    lists%firstL%part(lists%firstL%ncount)%sybra   = sybra
    lists%firstL%part(lists%firstL%ncount)%syket   = syket
    lists%firstL%part(lists%firstL%ncount)%i       = i
    lists%firstL%part(lists%firstL%ncount)%j       = j
    lists%firstL%part(lists%firstL%ncount)%k       = k
    lists%firstL%part(lists%firstL%ncount)%l       = 0
    lists%ncountL = lists%ncountL + 1

    If (itype == iLowerW    .Or. &
        itype == iLowerOneR .Or. &
        itype == iLowerOneL .Or. &
        itype == iLowerR) Then
       lists%firstL%part(lists%firstL%ncount)%l = CanonicalIndex(i,j)
    End If

    Return
  End Subroutine AddPartialLoop



  Subroutine OneBodyLoopComplete(itop, ybra, yket, sybra, syket, value, i, j)
    Implicit None
    Integer          :: itop, ybra, yket, sybra, syket
    Double Precision :: value
    Integer          :: i, j
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
    lists%first1%loop(lists%first1%ncount)%p      = i
    lists%first1%loop(lists%first1%ncount)%q      = j
    lists%ncount1 = lists%ncount1 + 1

    Return
  End Subroutine OneBodyLoopComplete



  Subroutine TwoBodyLoopComplete(itop, ybra, yket, sybra, syket, value, i, j, k, l)
    Implicit None
    Integer          :: itop, ybra, yket, sybra, syket
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
    !
    ! This reordering is done to get complete agreement of the orbital
    ! indices with the loop-driven implementation for debugging purposes.
    ! A disagreement may occur due to the mode of loop construction.
    ! For the top-down approach as in the original loop-driven
    ! algorithm, i is always the smallest orbital index.
    ! For the bottom-up approach as in this module (gugalower),
    ! l is always the largest index.
    !
    If (Min(i,j) <= k) Then
       lists%first2%loop(lists%first2%ncount)%p   = CanonicalIndex(i,j)
       lists%first2%loop(lists%first2%ncount)%q   = CanonicalIndex(k,l)
    Else
       lists%first2%loop(lists%first2%ncount)%p   = CanonicalIndex(k,l)
       lists%first2%loop(lists%first2%ncount)%q   = CanonicalIndex(i,j)
    End If
    lists%ncount2 = lists%ncount2 + 1

    Return
  End Subroutine TwoBodyLoopComplete


End Module gugalower
