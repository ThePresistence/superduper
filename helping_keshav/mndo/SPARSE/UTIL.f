      SUBROUTINE UTIL (F,NFS,N,IA,IB,CUTF)
      IMPLICIT NONE
      INTEGER :: NFS,N,IA,IB,K,J,I
      DOUBLE PRECISION :: CUTF,F(N,*)
C
C *** COUNT NUMBER OF NON-ZERO ELEMENTS IN COLUMNS IA-IB.
      DO 20 K=IA,IB
      J      = K-IA+1
      DO 10 I=1,N
      IF(ABS(F(I,J)).GT.CUTF) THEN
         NFS = NFS+1
      ENDIF
   10 CONTINUE
   20 CONTINUE
      RETURN
      END
