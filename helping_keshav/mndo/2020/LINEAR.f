      SUBROUTINE LINEAR (A,B,N,LM2,LM4)
C     *
C     STORE SYMMETRIC SQUARE MATRIX A(LM2,LM2) AS LINEAR ARRAY B(LM4).
C     B(LM4) CORRESPONDS TO THE LOWER TRIANGLE IN ROW-WISE ORDER.
C     ELEMENTS UP TO AND INCLUDING THE N-TH ROW ARE COPIED.
C     A AND B MAY OCCUPY THE SAME CORE SPACE.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LM2,LM2),B(LM4)
      IJ     = 0
      DO 20 I=1,N
      DO 10 J=1,I
      IJ     = IJ+1
      B(IJ)  = A(J,I)
   10 CONTINUE
   20 CONTINUE
      RETURN
      END