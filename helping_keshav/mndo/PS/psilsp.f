C     ******************************************************************
      SUBROUTINE PSILSP(R,ARO,BRO,BD,IMOD,W,W1,W2)
C
C  This subroutine is not intended to be human-readable.
C  It was generated by a Mathematica program (makeints.m).
C  General version of 14 May 1997.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARO(9), BRO(9), BD(5)
      DIMENSION W(461),W1(461),W2(461)
      PARAMETER (CN1=0.5D0)
      PARAMETER (CN2=-0.5D0)
      PARAMETER (CN3=2.D0)
      PARAMETER (CN4=0.25D0)
      PARAMETER (CN5=1.4142135623730950488016887242096980785697D0)
      PARAMETER (CN6=-1.4142135623730950488016887242096980785697D0)
      PARAMETER (CN7=3.D0)
      PARAMETER (CN8=1.5D0)
      PARAMETER (CN9=-1.5D0)
      PARAMETER (CN10=0.75D0)
      PARAMETER (CN11=-0.25D0)
      T1 = R**2
      T18 = SQRT(T1 + (ARO(1) + BRO(1))**2)
      W(1) = 1/T18
      T34 = SQRT(T1 + (ARO(2) + BRO(1))**2)
      W(2) = 1/T34
      T19 = SQRT(T1 + (ARO(1) + BRO(2))**2)
      W(16) = 1/T19
      T35 = SQRT(T1 + (ARO(2) + BRO(2))**2)
      W(17) = 1/T35
      T3 = (ARO(1) + BRO(3))**2
      T2 = -R
      T20 = T2 + BD(1)
      T9 = T20**2
      T21 = SQRT(T3 + T9)
      T22 = T2 - BD(1)
      T8 = T22**2
      T23 = SQRT(T3 + T8)
      W(51) = CN2/T21 + CN1/T23
      T10 = (ARO(2) + BRO(3))**2
      T36 = SQRT(T10 + T9)
      T37 = SQRT(T10 + T8)
      W(52) = CN2/T36 + CN1/T37
      T11 = CN3*BD(2)**2
      T4 = (ARO(1) + BRO(4))**2
      T25 = SQRT(T1 + T11 + T4)
      T24 = SQRT(T1 + (ARO(1) + BRO(5))**2)
      T6 = 1/T24
      T26 = SQRT(T1 + T4)
      T7 = CN2/T26
      T5 = CN1/T25 + T6 + T7
      W(66) = T5
      T38 = SQRT(T1 + (ARO(2) + BRO(5))**2)
      T14 = 1/T38
      T12 = (ARO(2) + BRO(4))**2
      T40 = SQRT(T1 + T12)
      T16 = CN2/T40
      T39 = SQRT(T1 + T11 + T12)
      T13 = T14 + T16 + CN1/T39
      W(67) = T13
      W(87) = T5
      W(88) = T13
      T29 = T2 + CN5*BD(2)
      T15 = T29**2
      T30 = SQRT(T15 + T4)
      T31 = T2 + CN6*BD(2)
      T17 = T31**2
      T32 = SQRT(T17 + T4)
      W(124) = CN4/T30 + CN4/T32 + T6 + T7
      T43 = SQRT(T12 + T15)
      T44 = SQRT(T12 + T17)
      W(125) = T14 + T16 + CN4/T43 + CN4/T44
      IF( IMOD.LE.0 ) RETURN
C
C   Expressions for first derivatives
C
      T46 = T18**2
      T47 = T18*T46
      T48 = 1/T47
      W1(1) = -(R*T48)
      T77 = T34**2
      T78 = T34*T77
      T79 = 1/T78
      W1(2) = -(R*T79)
      T49 = T19**2
      T50 = T19*T49
      T51 = 1/T50
      W1(16) = -(R*T51)
      T80 = T35**2
      T81 = T35*T80
      T82 = 1/T81
      W1(17) = -(R*T82)
      T52 = T23**2
      T53 = T23*T52
      T54 = 1/T53
      T55 = T21**2
      T56 = T21*T55
      T57 = 1/T56
      W1(51) = CN1*T22*T54 + CN2*T20*T57
      T83 = T37**2
      T84 = T37*T83
      T85 = 1/T84
      T86 = T36**2
      T87 = T36*T86
      T88 = 1/T87
      W1(52) = CN1*T22*T85 + CN2*T20*T88
      T63 = T24**2
      T64 = T24*T63
      T65 = 1/T64
      T28 = -(R*T65)
      T58 = T26**2
      T59 = T26*T58
      T62 = 1/T59
      T33 = CN1*R*T62
      T60 = T25**2
      T61 = T25*T60
      T66 = 1/T61
      T27 = T28 + T33 + CN2*R*T66
      W1(66) = T27
      T94 = T38**2
      T95 = T38*T94
      T96 = 1/T95
      T42 = -(R*T96)
      T89 = T40**2
      T90 = T40*T89
      T93 = 1/T90
      T45 = CN1*R*T93
      T91 = T39**2
      T92 = T39*T91
      T97 = 1/T92
      T41 = T42 + T45 + CN2*R*T97
      W1(67) = T41
      W1(87) = T27
      W1(88) = T41
      T70 = T32**2
      T71 = T32*T70
      T72 = 1/T71
      T74 = T30**2
      T75 = T30*T74
      T76 = 1/T75
      W1(124) = T28 + T33 + CN4*T31*T72 + CN4*T29*T76
      T105 = T43**2
      T106 = T105*T43
      T107 = 1/T106
      T101 = T44**2
      T102 = T101*T44
      T103 = 1/T102
      W1(125) = CN4*T107*T29 + CN4*T103*T31 + T42 + T45
      IF( IMOD.LE.1 ) RETURN
C
C   Expressions for second derivatives
C
      W2(1) = CN7*T1/(T46*T47) - T48
      W2(2) = CN7*T1/(T77*T78) - T79
      W2(16) = CN7*T1/(T49*T50) - T51
      W2(17) = CN7*T1/(T80*T81) - T82
      W2(51) = CN2*T54 + CN1*T57 + CN8*T8/(T52*T53) + CN9*T9/(T55*T56)
      W2(52) = CN8*T8/(T83*T84) + CN2*T85 + CN1*T88 + CN9*T9/(T86*T87)
      T68 = CN9*T1/(T58*T59)
      T69 = CN1*T62
      T73 = CN7*T1/(T63*T64) - T65
      T67 = CN8*T1/(T60*T61) + CN2*T66 + T68 + T69 + T73
      W2(66) = T67
      T100 = CN1*T93
      T104 = CN7*T1/(T94*T95) - T96
      T99 = CN9*T1/(T89*T90)
      T98 = T100 + T104 + CN8*T1/(T91*T92) + CN2*T97 + T99
      W2(67) = T98
      W2(87) = T67
      W2(88) = T98
      W2(124) = T68 + T69 + CN10*T17/(T70*T71) + CN11*T72 + T73 + 
     1 CN10*T15/(T74*T75) + CN11*T76
      W2(125) = T100 + CN11*T103 + T104 + CN11*T107 + 
     1 CN10*T15/(T105*T106) + CN10*T17/(T101*T102) + T99
      RETURN
      END