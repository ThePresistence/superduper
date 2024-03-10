      SUBROUTINE MATOUT(A,B,NR,NC,LM2,LM3)
C     *
C     PRINT EIGENVALUES B(LM3) AND EIGENVECTORS A(LM2,LM3).
C     THE OUTPUT CONTAINS NR ROWS AND NC COLUMNS OF THE MATRIX A.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /NBFILE/ NBF(20)
      DIMENSION  A(LM2,LM3),B(LM3)
      NB6    = NBF(6)
      KA     = 1
      KC     = 10
   10 KB     = MIN(KC,NC)
      WRITE(NB6,500) (I,I=KA,KB)
      WRITE(NB6,510) (B(I),I=KA,KB)
      WRITE(NB6,520)
      N      = 0
      DO 20 I=1,NR
      WRITE(NB6,530) I,(A(I,J),J=KA,KB)
      N      = N+1
      IF(N.LT.10) GO TO 20
      WRITE(NB6,520)
      N      = 0
   20 CONTINUE
      IF(KB.EQ.NC) RETURN
      KA     = KC+1
      KC     = KC+10
      GO TO 10
  500 FORMAT(// 9X,I5,9I12)
  510 FORMAT(/  5X,10F12.5)
  520 FORMAT(   1X)
  530 FORMAT(   1X,I4,10F12.5)
      END