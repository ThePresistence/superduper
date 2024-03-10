      SUBROUTINE SPGTO2 (I,ZB,RAB,S1,S2)
C     *
C     CORE-VALENCE OVERLAP INTEGRALS FOR GTOS.
C     *
C     CORE    ORBITALS: STO-3G GAUSSIANS WITH SCALING FACTOR ZB.
C     VALENCE ORBITALS: CURRENT GAUSSIANS FROM SHELL I.
C     *
C     I      SHELL NUMBER FOR VALENCE SHELL OF ATOM A AT ORIGIN.
C     ZB     SCALING FACTOR FOR  CORE SHELL OF ATOM B AT (0,0,RAB).
C     RAB    DISTANCE A-B (IN AU).
C     S1     INTEGRAL S(A)-S(B).
C     S2     INTEGRAL SIGMA(A)-S(B).
C     *
      USE LIMIT, ONLY: LMGP, LMGS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (SIXTY=60.0D0)
C     UNSCALED STO-3G BASIS FROM SUBROUTINE STONG.
      PARAMETER (EX1=2.227660584D+00)
      PARAMETER (CS1=1.543289673D-01)
      PARAMETER (EX2=4.057711562D-01)
      PARAMETER (CS2=5.353281423D-01)
      PARAMETER (EX3=1.098175104D-01)
      PARAMETER (CS3=4.446345422D-01)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./GAUSS1/ KSTART(LMGS),KNG(LMGS),KTYPE(LMGS),NSHELL,NBASIS
     ./GAUSS3/ EXX(LMGP),C1(LMGP),C2(LMGP),C3(LMGP)
      DIMENSION CC(3),EC(3)
C *** INITIALIZATION.
      S1     = ZERO
      S2     = ZERO
      ISTART = KSTART(I)
      INUMB  = KNG(I)
      ABSQ   = RAB*RAB
C     SCALE STO-3G BASIS FOR CORE ORBITAL.
      EC(1)  = EX1*ZB**2
      EC(2)  = EX2*ZB**2
      EC(3)  = EX3*ZB**2
      CC(1)  = CS1*(TWO*EC(1)/PI)**0.75D0
      CC(2)  = CS2*(TWO*EC(2)/PI)**0.75D0
      CC(3)  = CS3*(TWO*EC(3)/PI)**0.75D0
C *** BRANCH ACCORDING TO NATURE OF ATOM A.
      IF(KTYPE(I).GT.0) GO TO 50
C *** CASE 0 0. ATOM A IS HYDROGEN OR HELIUM.
      DO 40 II=1,INUMB
      IG     = ISTART+II-1
      A      = EXX(IG)
      CSA    = C1(IG)
      DO 30 JJ=1,3
      B      = EC(JJ)
      CSB    = CC(JJ)
      AB     = A*B
      G      = A+B
      E      = ONE/G
      XQQ    = AB*ABSQ*E
      IF(XQQ.GT.SIXTY) GO TO 30
      S00    = (PI*E)**1.5D0 * EXP(-XQQ)
      S1     = S1 + CSA*CSB*S00
   30 CONTINUE
   40 CONTINUE
      RETURN
C *** CASE 1 0. ATOM A IS FIRST-ROW ATOM (Li-Ne).
   50 DO 70 II=1,INUMB
      IG     = ISTART+II-1
      A      = EXX(IG)
      CSA    = C1(IG)
      CPA    = C2(IG)
      DO 60 JJ=1,3
      B      = EC(JJ)
      CSB    = CC(JJ)
      AB     = A*B
      G      = A+B
      E      = ONE/G
      XQQ    = AB*ABSQ*E
      IF(XQQ.GT.SIXTY) GO TO 60
      S00    = (PI*E)**1.5D0 * EXP(-XQQ)
      S30    =  B*RAB*S00*E
      S1     = S1 + CSA*CSB*S00
      S2     = S2 + CPA*CSB*S30
   60 CONTINUE
   70 CONTINUE
      RETURN
      END
