      SUBROUTINE DMAT (A,B,C,D,LMA,LMB,LMC,LMD,LMQ,NPS)
C     *
C     COMPUTATION OF THE PRODUCT B * INVERSE(A) * TRANS(B).
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     A(LMA)    INVERSE A-MATRIX FOR SURFACE CHARGE INTERACTIONS (I).
C     B(LMB)    B-MATRIX FOR SOLUTE/SURFACE INTERACTIONS (O).
C     C(LMC)    SCRATCH SPACE (S).
C     D(LMD)    DIELECTRIC OPERATOR MATRIX (O).
C     LMQ       NUMBER OF SOLUTE CHARGES, LMQ=LM6+NUMAT (I).
C     NPS       NUMBER OF SURFACE CHARGES (I).
C     *
C     LOGICAL DIMENSIONS OF THE ARRAYS FOR MODCSM.LE.0:
C     A(NPS,NPS), B(LMQ,NPS), C(LMQ,NPS), D(LMQ,LMQ).
C     *
C     LOGICAL DIMENSIONS OF THE ARRAYS FOR MODCSM.GE.1:
C     A(NPSLIN), B(LMQ,NPS), C(NPS), D(LMQLIN).
C     *
      USE LIMIT, ONLY: LM1, LMI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./COSMO2/ NCOSMO,NPSG(4)
     ./COSMO3/ FEPSI,RDS,DISEX2,SRAD(LM1)
     ./COSMO7/ JPC(LMI+LM1+1)
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
      DIMENSION A(LMA),B(LMB),C(LMC),D(LMD)
C *** INITIALIZATION.
      NPRINT = IN2(72)
      MODCSM = IN2(235)
      LMQLIN = (LMQ*(LMQ+1))/2
      NPSLIN = (NPS*(NPS+1))/2
C *** CONVERSION FROM ANGSTROM-1 TO ELECTRON VOLTS (EV*A0).
C     INCLUSION OF EPSILON-DEPENDENT CORRECTION FACTOR (FEPSI).
C     AN ADDITIONAL FACTOR OF -0.5 IN THE ORIGINAL COSMO IMPLEMENTATION
C     HAS BEEN REMOVED TO BE CONSISTENT WITH THE PUBLISHED FORMULAS.
C     CORRESPONDINGLY, FACTORS OF -2 HAVE BEEN REMOVED IN OTHER ROUTINES
C     WHERE THE DIELECTRIC OPERATOR MATRIX IS USED.
      FAC    = EV*A0*FEPSI
COLD  FAC    = TWO*13.6058D0*0.5292D0*FEPSI
C *** STORE INVERSE A-MATRIX AS LOWER TRIANGLE.
      IF(MODCSM.GE.1) THEN
         CALL LINEAR (A,A,NPS,NPS,NPSLIN)
      ENDIF
C *** CALCULATE B-MATRIX.
      CALL BMAT3 (B,LMQ,NPS)
C *** LIBRARY CALLS TO COMPUTE D = B * INVERSE(A) * TRANS(B).
      IF(MODCSM.LE.0) THEN
         CALL ZZERO (LMQ,NPS,C,LMQ)
         CALL ZZERO (LMQ,LMQ,D,LMQ)
         CALL DGEMM ('N','N',LMQ,NPS,NPS,ONE,B,LMQ,A,NPS,ZERO,C,LMQ)
         CALL DGEMM ('N','T',LMQ,LMQ,NPS,FAC,C,LMQ,B,LMQ,ZERO,D,LMQ)
C *** FORTRAN CODE TO COMPUTE D = B * INVERSE(A) * TRANS(B).
      ELSE
         DO 70 I=1,LMQ
         DO 40 K=1,NPS
         SUM    = ZERO
         DO 20 L=1,K
         SUM    = SUM+B(I+(L-1)*LMQ)*A(JPC(K)+L)
   20    CONTINUE
         DO 30 L=K+1,NPS
         SUM    = SUM+B(I+(L-1)*LMQ)*A(JPC(L)+K)
   30    CONTINUE
         C(K)   = SUM
   40    CONTINUE
         DO 60 J=1,I
         DIJ    = ZERO
         DO 50 K=1,NPS
         DIJ    = DIJ+C(K)*B(J+(K-1)*LMQ)
   50    CONTINUE
         IJ     = JPC(I)+J
         D(IJ)  = FAC*DIJ
   60    CONTINUE
   70    CONTINUE
      ENDIF
C *** DEBUG PRINT OF DIELECTRIC OPERATOR MATRIX.
      IF(NCOSMO.EQ.1 .AND. NPRINT.GT.5) THEN
         NB6 = NBF(6)
         WRITE(NB6,500)
         IF(MODCSM.LE.0) THEN
            CALL MATPRT (D,LMQ,LMQ,LMQ,LMQ)
         ELSE
            CALL VECPRT (D,LMQLIN,LMQ)
         ENDIF
         WRITE(NB6,510) FEPSI
      ENDIF
      RETURN
  500 FORMAT(///1X,'DIELECTRIC OPERATOR MATRIX (EV).'/)
  510 FORMAT(// 1X,'THE MATRIX ELEMENTS CONTAIN THE EPSILON-DEPENDENT',
     1          1X,'COSMO CORRECTION FACTOR: FEPSI =',F8.5)
      END
