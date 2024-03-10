      SUBROUTINE DIELEN (EDIE,PA,PB,Q,D2,D,LM4,LMD2,LMQ)
C     *
C     CALCULATION OF THE DIELECTRIC SCREENING ENERGY.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     EDIE      DIELECTRIC SCREENING ENERGY (O).
C     PA(LM4)   DENSITY MATRIX FOR ALPHA ELECTRONS (I).
C     PB(LM4)   DENSITY MATRIX FOR BETA ELECTRONS (I).
C     Q(LMQ)    POPULATIONS OF ONE-CENTER PAIRS (S).
C     D2(LMD2)  PACKED DIELECTRIC OPERATOR MATRIX (I).
C     D(LMQ,*)  SQUARE DIELECTRIC OPERATOR MATRIX (I).
C     *
      USE LIMIT, ONLY: LM1, LMI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./COSMO7/ JPC(LMI+LM1+1)
     ./INOPT2/ IN2(300)
      DIMENSION PA(LM4),PB(LM4),Q(LMQ),D2(LMD2),D(LMQ,LMQ)
C *** INITIALIZATION.
      MODCSM = IN2(235)
C *** ONE-CENTER PAIR POPULATIONS AND CORE CHARGES.
      CALL PQDEN (PA,PB,Q,LM4,LMQ)
C *** COMPUTATION OF THE SCREENING ENERGY.
      EDIE   = 0.D0
      DO 50 I=1,LMQ
      QI     = Q(I)
      IF(MODCSM.GT.0) THEN
         ID     = JPC(I)
         DO 30 J=1,I-1
         EDIE   = EDIE-QI*D2(ID+J)*Q(J)
   30    CONTINUE
         EDIE   = EDIE-0.5D0*QI*D2(ID+I)*QI
      ELSE
         DO 40 J=1,I-1
         EDIE   = EDIE-QI*D(J,I)*Q(J)
   40    CONTINUE
         EDIE   = EDIE-0.5D0*QI*D(I,I)*QI
      ENDIF
   50 CONTINUE
      RETURN
      END
