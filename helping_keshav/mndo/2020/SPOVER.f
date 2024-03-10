      SUBROUTINE SPOVER (I,J,RAB,S)
C     *
C     ANALYTICAL TWO-CENTER OVERLAP INTEGRALS FOR GTOs (SP BASIS).
C     SIMPLIFIED VERSION OF SUBROUTINE SPSTV (IN PP95).
C     *
C     NOTATION.  I=INPUT, O=OUTPUT.
C     I,J        SHELL NUMBERS (I).
C     RAB        DISTANCE BETWEEN CENTERS OF GAUSSIANS A AND B (I).
C     S          OVERLAP INTEGRALS IN LOCAL COORDINATES (O).
C     *
      USE LIMIT, ONLY: LMGP, LMGS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (SIXTY=60.0D0)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./GAUSS1/ KSTART(LMGS),KNG(LMGS),KTYPE(LMGS),NSHELL,NBASIS
     ./GAUSS3/ EXX(LMGP),C1(LMGP),C2(LMGP),C3(LMGP)
      DIMENSION S(14)
C *** INITIALIZATION.
      DO 10 K=1,5
      S(K)   = ZERO
   10 CONTINUE
      ISTART = KSTART(I)
      INUMB  = KNG(I)
      JSTART = KSTART(J)
      JNUMB  = KNG(J)
      ABZ    =-RAB
      RABSQ  = RAB*RAB
C *** AVAILABLE SHELL TYPES.
C     0 0  S S  :  KTYPE(I)=0  KTYPE(J)=0
C     0 1  S P  :  KTYPE(I)=0  KTYPE(J)=1
C     1 0  P S  :  KTYPE(I)=1  KTYPE(J)=0
C     1 1  P P  :  KTYPE(I)=1  KTYPE(J)=1
C *** CASE 0 0. ATOM PAIR H-H.
      IF(KTYPE(I).EQ.0 .AND. KTYPE(J).EQ.0) THEN
         DO 25 II=1,INUMB
         IG     = ISTART+II-1
         A      = EXX(IG)
         CSA    = C1(IG)
         DO 20 JJ=1,JNUMB
         JG     = JSTART+JJ-1
         B      = EXX(JG)
         AB     = A*B
         G      = A+B
         E      = ONE/G
         XQQ    = AB*RABSQ*E
         IF(XQQ.GT.SIXTY) GO TO 20
         PIE    = PI*E
         CSB    = C1(JG)
         S00    = PIE*SQRT(PIE)* EXP(-XQQ)
         S(1)   = S(1) + CSA*CSB*S00
   20    CONTINUE
   25    CONTINUE
C *** CASE 0 1. ATOM PAIR H-X.
      ELSE IF(KTYPE(I).EQ.0 .AND. KTYPE(J).EQ.1) THEN
         DO 35 II=1,INUMB
         IG     = ISTART+II-1
         A      = EXX(IG)
         CSA    = C1(IG)
         DO 30 JJ=1,JNUMB
         JG     = JSTART+JJ-1
         B      = EXX(JG)
         AB     = A*B
         G      = A+B
         E      = ONE/G
         XQQ    = AB*RABSQ*E
         IF(XQQ.GT.SIXTY) GO TO 30
         PIE    = PI*E
         S00    = PIE*SQRT(PIE)* EXP(-XQQ)
         S03    = A*ABZ*S00*E
         CSB    = C1(JG)
         CPB    = C2(JG)
         S(1)   = S(1) + CSA*CSB*S00
         S(2)   = S(2) + CSA*CPB*S03
   30    CONTINUE
   35    CONTINUE
C *** CASE 1 0. ATOM PAIR X-H.
      ELSE IF(KTYPE(I).EQ.1 .AND. KTYPE(J).EQ.0) THEN
         DO 45 II=1,INUMB
         IG     = ISTART+II-1
         A      = EXX(IG)
         CSA    = C1(IG)
         CPA    = C2(IG)
         DO 40 JJ=1,JNUMB
         JG     = JSTART+JJ-1
         B      = EXX(JG)
         AB     = A*B
         G      = A+B
         E      = ONE/G
         XQQ    = AB*RABSQ*E
         IF(XQQ.GT.SIXTY) GO TO 40
         PIE    = PI*E
         S00    = PIE*SQRT(PIE)* EXP(-XQQ)
         S30    =-B*ABZ*S00*E
         CSB    = C1(JG)
         S(1)   = S(1) + CSA*CSB*S00
         S(3)   = S(3) + CPA*CSB*S30
   40    CONTINUE
   45    CONTINUE
C *** CASE 1 1. ATOM PAIR X-Y.
      ELSE IF(KTYPE(I).EQ.1 .AND. KTYPE(J).EQ.1) THEN
         DO 55 II=1,INUMB
         IG     = ISTART+II-1
         A      = EXX(IG)
         CSA    = C1(IG)
         CPA    = C2(IG)
         DO 50 JJ=1,JNUMB
         JG     = JSTART+JJ-1
         B      = EXX(JG)
         AB     = A*B
         G      = A+B
         E      = ONE/G
         XQQ    = AB*RABSQ*E
         IF(XQQ.GT.SIXTY) GO TO 50
         PIE    = PI*E
         S00    = PIE*SQRT(PIE)* EXP(-XQQ)
         S03    = A*ABZ*S00*E
         S30    =-B*ABZ*S00*E
         S11    = S00*PT5*E
         S33    = S00*(PT5*E-AB*ABZ*ABZ*E*E)
         CSB    = C1(JG)
         CPB    = C2(JG)
         S(1)   = S(1) + CSA*CSB*S00
         S(2)   = S(2) + CSA*CPB*S03
         S(3)   = S(3) + CPA*CSB*S30
         S(4)   = S(4) + CPA*CPB*S33
         S(5)   = S(5) + CPA*CPB*S11
   50    CONTINUE
   55    CONTINUE
      ENDIF
      RETURN
      END
