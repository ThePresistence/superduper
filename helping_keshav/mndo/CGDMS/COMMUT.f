C     ******************************************************************
      DOUBLE PRECISION FUNCTION COMMUT (A,LM2,N)
C     *
C     GIVEN: MATRIX PRODUCT A = B*C OF TWO SYMMETRIC MATRICES B AND C.
C     FIND : MAXIMUM ABSOLUTE VALUE OF COMMUTATOR MATRIX ELEMENTS.
C     NOTE : C*B = C(transposed)*B(transposed) = (B*C)transposed
C     HENCE: [B,C] = B*C - C*B = B*C - (B*C)transposed
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     A(LM2,*)  INPUT MATRIX (I).
C     N         NUMBER OF ORBITALS (I).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LM2,N)
      COMMUT = 0.0D0
      DO 20 J=2,N
      DO 10 I=1,J
      COMMUT = MAX(COMMUT,ABS((A(I,J)-A(J,I))))
   10 CONTINUE
   20 CONTINUE
      RETURN
      END
