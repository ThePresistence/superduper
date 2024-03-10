      SUBROUTINE WWSTEP (C12,CC,WW,LM7,NB3,KMAX,LMAX)
C     *
C     CALCULATION OF A SET OF (IJ,AB) INTEGRALS.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C12(KMAX),CC(LM7),WW(KMAX)
      IF(KMAX.GT.LMAX) GO TO 20
C     AO INTEGRALS IN CC(LM7).
      KK     = 1-KMAX
      DO 10 NN=1,KMAX
      KK     = KK+KMAX
      WW(NN) = DDOT(KMAX,C12,1,CC(KK),1)
   10 CONTINUE
      RETURN
C     AO INTEGRALS ON FILE NB3.
   20 REWIND NB3
      READ(NB3) CC
      KK     = 1-KMAX
      LL     = 0
      DO 30 NN=1,KMAX
      LL     = LL+1
      IF(LL.GT.LMAX) THEN
         READ(NB3) CC
         KK  = 1-KMAX
         LL  = 1
      ENDIF
      KK     = KK+KMAX
      WW(NN) = DDOT(KMAX,C12,1,CC(KK),1)
   30 CONTINUE
      RETURN
      END