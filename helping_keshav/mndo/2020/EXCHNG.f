      SUBROUTINE EXCHNG (A,B,C,D,E,F,N)
C     *
C     COPY A TO B, C TO D, AND ARRAY E(N) TO ARRAY F(N).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(N),F(N)
      B      = A
      D      = C
      DO 10 I=1,N
      F(I)   = E(I)
   10 CONTINUE
      RETURN
      END
