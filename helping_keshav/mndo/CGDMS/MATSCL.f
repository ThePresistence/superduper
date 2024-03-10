C     ******************************************************************
      SUBROUTINE MATSCL (A,B,FACTOR,LM2,N)
C     *
C     SCALE MATRIX BY A FACTOR.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     A(LM2,*)  INPUT MATRIX (I).
C     B(LM2,*)  OUTPUT MATRIX (O).
C     FACTOR    SCALE FACTOR.
C     LM2       LEADING DIMENSION (I).
C     N         NUMBER OF ORBITALS (I).
C     *
C     ARRAYS A AND B MAY OCCUPY THE SAME CORE SPACE.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LM2,N),B(LM2,N)
      DO 20 J=1,N
      DO 10 I=1,N
      B(I,J) = FACTOR*A(I,J)
   10 CONTINUE
   20 CONTINUE
      RETURN
      END
