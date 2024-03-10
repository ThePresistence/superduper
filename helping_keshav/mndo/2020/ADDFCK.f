      SUBROUTINE ADDFCK (F,PA,PB,Q,D2,D,LM4,LMD2,LMQ)
C     *
C     ADD THE COSMO CONTRIBUTIONS TO THE FOCK MATRIX.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     F(LM4)    FOCK MATRIX (I,O).
C     PA(LM4)   DENSITY MATRIX FOR ALPHA ELECTRONS (I).
C     PB(LM4)   DENSITY MATRIX FOR BETA ELECTRONS (I).
C     Q(LMQ)    EFFECTIVE POPULATIONS OF ONE-CENTER PAIRS (S).
C     D2(LMD2)  PACKED DIELECTRIC OPERATOR MATRIX (I).
C     D(LMQ,*)  SQUARE DIELECTRIC OPERATOR MATRIX (I).
C     *
      USE LIMIT, ONLY: LM1, LMI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./COSMO7/ JPC(LMI+LM1+1)
     ./FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
     ./INOPT2/ IN2(300)
      DIMENSION F(LM4),PA(LM4),PB(LM4),Q(LMQ),D2(LMD2),D(LMQ,LMQ)
C *** INITIALIZATION.
      MODCSM = IN2(235)
      LM6    = LMQ-NUMAT
      DO 10 J=1,LM6
      JM     = IP(J)
      Q(J)   = PA(JM)+PB(JM)
      IF(IP1(J).NE.IP2(J)) Q(J)=Q(J)*2.0D0
   10 CONTINUE
C *** OUTER LOOP OVER ALL ONE CENTER-PAIRS.
      DO 40 I=1,LM6
      IM     = IP(I)
      FIM    = 0.D0
C *** INNER LOOP OVER ALL ONE CENTER-PAIRS.
C     CODE FOR PACKED DIELECTRIC OPERATOR MATRIX.
      IF(MODCSM.GT.0) THEN
         DO 20 J=1,LM6
         MAXIJ  = MAX(I,J)
         IJ     = JPC(MAXIJ)+MIN(I,J)
         FIM    = FIM - D2(IJ) * Q(J)
   20    CONTINUE
      ELSE
C     CODE FOR SQUARE DIELECTRIC OPERATOR MATRIX.
         DO 30 J=1,LM6
         FIM    = FIM - D(J,I) * Q(J)
   30    CONTINUE
      ENDIF
      F(IM)  = F(IM) + FIM
   40 CONTINUE
      RETURN
      END
