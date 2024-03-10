      SUBROUTINE VECPRT (A,ISIZE,NUMB)
C     *
C     PRINT PACKED LINEAR ARRAY A(ISIZE) AS LOWER TRIANGLE OF A MATRIX.
C     THE OUTPUT CONTAINS NUMB ROWS AND COLUMNS.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*4 LINE(31)
      COMMON /NBFILE/ NBF(20)
      DIMENSION A(ISIZE)
      DATA LINE / 31*'----' /
      IF(NUMB.LE.0) RETURN
      NB6    = NBF(6)
      LIMIT  = (NUMB*(NUMB+1))/2
      NA     = 1
   20 M      = MIN((NUMB+1-NA),10)
      MA     = 3*M+1
      M      = NA+M-1
      WRITE(NB6,500) (I,I=NA,M)
      WRITE(NB6,510) (LINE(K), K=1,MA)
      DO 30 I=NA,NUMB
      K      = (I*(I-1))/2
      L      = MIN((K+M),(K+I) )
      K      = K+NA
      WRITE(NB6,520) I,(A(N), N=K,L)
   30 CONTINUE
      IF(L.GE.LIMIT) RETURN
      NA     = M+1
      GO TO 20
  500 FORMAT(// 2X,10(7X,I5))
  510 FORMAT(   1X,31A4)
  520 FORMAT(   1X,I3,1X,10F12.5)
      END
