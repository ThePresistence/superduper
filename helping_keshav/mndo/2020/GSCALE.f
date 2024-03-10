      SUBROUTINE GSCALE (R,RI,FKO,NI,NJ)
C     *
C     KLOPMAN-OHNO SCALING OF TWO-CENTER TWO-ELECTRON INTEGRALS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     R         INTERATOMIC DISTANCE, IN BOHR (I).
C     RI(22)    ANALYTICAL INTEGRALS (I).
C     RI(22)    SCALED ANALYTICAL INTEGRALS (O).
C     FKO       KLOPMAN-OHNO FACTOR (O).
C     NI,NJ     ATOMIC NUMBERS (I).
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DPARM / LORBS(LMZ)
     ./REP   / GSS(LMZ),REPDUM(LMZ*5)
C    ./INOPT2/ IN2(300)
      DIMENSION RI(22)
C *** INPUT VARIABLES.
C     NPRINT = IN2(72)
C *** INITIALIZATION.
      IORBS  = LORBS(NI)
      JORBS  = LORBS(NJ)
C *** KLOPMAN-OHNO FACTOR.
      POI    = PT5*EV/GSS(NI)
      POJ    = PT5*EV/GSS(NJ)
      AEE    = (POI+POJ)**2
      SEMI   = EV/SQRT(R*R+AEE)
      FKO    = SEMI/RI(1)
C *** SCALING OF THE INTEGRALS.
      IF(IORBS.EQ.1 .AND. JORBS.EQ.1) THEN
         RI(1)  = RI(1)*FKO
      ELSE IF(IORBS.EQ.1 .AND. JORBS.EQ.4) THEN
         RI(1)  = RI(1)*FKO
         RI(5)  = RI(5)*FKO
         RI(11) = RI(11)*FKO
         RI(12) = RI(12)*FKO
      ELSE IF(IORBS.EQ.4 .AND. JORBS.EQ.1) THEN
         DO 10 I=1,4
         RI(I)  = RI(I)*FKO
   10    CONTINUE
      ELSE IF(IORBS.EQ.4 .AND. JORBS.EQ.4) THEN
         DO 20 I=1,22
         RI(I)  = RI(I)*FKO
   20    CONTINUE
      ENDIF
      RETURN
      END