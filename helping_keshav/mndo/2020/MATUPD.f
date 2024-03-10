C     ******************************************************************
      SUBROUTINE MATUPD (A,B,C,FACTOR,LM2,N)
C     *
C     MATRIX UPDATE: C = A + FACTOR*B.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     A(LM2,*)  FIRST  INPUT MATRIX (I).
C     B(LM2,*)  SECOND INPUT MATRIX (I).
C     C(LM2,*)  RESULTING MATRIX (O).
C     FACTOR    UPDATE FACTOR (I).
C     LM2       LEADING DIMENSION (I).
C     N         NUMBER OF ORBITALS (I).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LM2,N),B(LM2,N),C(LM2,N)
      DO 20 J=1,N
      DO 10 I=1,N
      C(I,J) = A(I,J) + FACTOR*B(I,J)
   10 CONTINUE
   20 CONTINUE
      RETURN
      END
