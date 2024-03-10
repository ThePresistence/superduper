      SUBROUTINE MATPRT(A,NR,NC,LM2,LM3)
C     *
C     PRINT THE MATRIX A(LM2,LM3), NR ROWS AND NC COLUMNS.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /NBFILE/ NBF(20)
      DIMENSION  A(LM2,LM3)
      NB6    = NBF(6)
      KA     = 1
      KC     = 10
   10 KB     = MIN(KC,NC)
      WRITE(NB6,500) (I,I=KA,KB)
      WRITE(NB6,510)
      N      = 0
      DO 20 I=1,NR
      WRITE(NB6,520) I,(A(I,J),J=KA,KB)
      N      = N+1
      IF(N.LT.10) GO TO 20
      WRITE(NB6,510)
      N      = 0
   20 CONTINUE
      IF(KB.EQ.NC) RETURN
      KA     = KC+1
      KC     = KC+10
      GO TO 10
  500 FORMAT(// 9X,I5,9I12)
  510 FORMAT(   1X)
  520 FORMAT(   1X,I4,10F12.5)
      END
