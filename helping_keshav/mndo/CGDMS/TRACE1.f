C     ******************************************************************
      DOUBLE PRECISION FUNCTION TRACE1 (A,LM2,N)
C     *
C     TRACE OF MATRIX A.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     A(LM2,*)  INPUT MATRIX (I).
C     LM2       LEADING DIMENSION (I).
C     N         NUMBER OF ORBITALS (I).
C     TRACE1    TRACE OF THE MATRIX (O).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LM2,N)
      TRACE1 = 0.0D0
      DO 10 I=1,N
      TRACE1 = TRACE1+A(I,I)
   10 CONTINUE
      RETURN
      END
