      SUBROUTINE SCREEN (PA,PB,Q,A2,A,B,LM4,LMA2,LMQ,NPS)
C     *
C     CALCULATION OF THE COSMO SCREENING CHARGES.
C     *
C     NOTATION.  I=INPUT, O=OUTPUT, S=SCRATCH.
C     PA(LM4)    DENSITY MATRIX FOR ALPHA ELECTRONS (I).
C     PB(LM4)    DENSITY MATRIX FOR BETA ELECTRONS (I).
C     Q(LMQ)     POPULATIONS OF ONE-CENTER PAIRS (S).
C     A2(LMA2)   INVERSE OF COSMO A-MATRIX, PACKED (I).
C     A(NPS,NPS) INVERSE OF COSMO A-MATRIX, SQUARE (I).
C     B(LMQ,NPS) RECTANGULAR COSMO B-MATRIX (I)
C     LMQ        TOTAL NUMBER OF SOLUTE CHARGES (I).
C     NPS        TOTAL NUMBERS OF SEGMENTS ON COSMO SURFACE (I).
C     *
C     OUTPUT IN COMMON BLOCK COSMO8.
C     QS(LMNPS)  SCREENING CHARGES (O).
C     *
      USE LIMIT, ONLY: LM1, LMI, LMNPS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./COSMO7/ JPC(LMI+LM1+1)
     ./COSMO8/ COSURF(3,LMNPS),AREASG(LMNPS),QS(LMNPS)
     ./INOPT2/ IN2(300)
      DIMENSION PA(LM4),PB(LM4),Q(LMQ)
      DIMENSION A2(LMA2),A(NPS,NPS),B(LMQ,NPS)
C *** READ COSMO MATRICES FROM FILE (IF NEEDED).
      CALL CODISK (A,B,LMA2,LMQ*NPS,0)
C *** INITIALIZATION.
      MODCSM = IN2(235)
C *** ONE-CENTER PAIR POPULATIONS AND CORE CHARGES.
      CALL PQDEN (PA,PB,Q,LM4,LMQ)
C *** INITIALIZE THE SCREENING CHARGES.
      DO 30 I=1,NPS
      QS(I)  = ZERO
   30 CONTINUE
C *** CALCULATE THE SCREENING CHARGES.
      DO 80 I=1,NPS
      POSI   = ZERO
      DO 40 J=1,LMQ
      POSI   = POSI+Q(J)*B(J,I)
   40 CONTINUE
      IF(MODCSM.GT.0) THEN
         DO 50 K=1,I
         QS(K)  = QS(K)+POSI*A2(K+JPC(I))
   50    CONTINUE
         DO 60 K=I+1,NPS
         QS(K)  = QS(K)+POSI*A2(I+JPC(K))
   60    CONTINUE
      ELSE
         DO 70 K=1,NPS
         QS(K)  = QS(K)+POSI*A(K,I)
   70    CONTINUE
      ENDIF
   80 CONTINUE
      RETURN
      END
