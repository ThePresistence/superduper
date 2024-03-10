      SUBROUTINE FOCK2(F,PA,PB,Q,W,LM4,LM6,LM9)
C     *
C     TWO-CENTER TWO-ELECTRON CONTRIBUTIONS TO THE FOCK MATRIX.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     F(LM4)    FOCK MATRIX (I,O).
C     PA(LM4)   RHF OR UHF-ALPHA DENSITY MATRIX (I).
C     PB(LM4)   UHF-BETA DENSITY MATRIX (I).
C     Q(LM6)    SCRATCH ARRAY FOR ONE-CENTER PAIR TERMS (S).
C     W(LM9)    TWO-ELECTRON INTEGRALS (I).
C     *
      USE LIMIT, ONLY: LM1, LMI, LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IEQJ,KEQL
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./FLAG4 / INTSUM
     ./FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
     ./INDEX / INDX(LMX)
     ./INDEXW/ NW(LM1)
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
      DIMENSION F(LM4),PA(LM4),PB(LM4),Q(LM6),W(LM9)
      DIMENSION IWW(4,4)
      DATA IWW/1,2,4,7,2,3,5,8,4,5,6,9,7,8,9,10/
C *** INPUT OPTION.
      IOP = IN2(2)
      IF(IOP.GT.0) GO TO 200
C *** MNDO.
      IF(INTSUM.GT.LM9) THEN
         NB2 = NBF(2)
         REWIND NB2
         READ(NB2) W
      ENDIF
      DO 10 KL=1,LM6
      Q(KL)  = (PA(IP(KL))+PB(IP(KL)))
      IF(IP1(KL).NE.IP2(KL)) Q(KL)=Q(KL)*TWO
   10 CONTINUE
      KK     = 0
      DO 160 II=2,NUMAT
      IA     = NFIRST(II)
      IB     = NLAST(II)
      IC     = INDX(IA)
      IORBS  = IB-IA+1
      IW     = INDX(IORBS)+IORBS
      IQ     = NW(II)-1
      DO 150 JJ=1,II-1
      JA     = NFIRST(JJ)
      JB     = NLAST(JJ)
      JORBS  = JB-JA+1
      JW     = INDX(JORBS)+JORBS
      JQ     = NW(JJ)-1
      IF(INTSUM.GT.LM9) THEN
         KW  = KK+IW*JW
         IF(KW.GT.LM9)  THEN
            READ(NB2) W
            KK = 0
         ENDIF
      ENDIF
      IF(JA.EQ.JB .AND. IA.EQ.IB) THEN
         IJ  = INDX(IA)+IA
         KL  = INDX(JA)+JA
         IK  = IJ-IA+JA
         A   = W(KK+1)
         F(IJ) = F(IJ)+A*Q(JQ+1)
         F(KL) = F(KL)+A*Q(IQ+1)
         F(IK) = F(IK)-A*PA(IK)
      ELSE IF(IA.EQ.IB .AND. JW.EQ.10) THEN
         PSS = Q(IQ+1)
         IPS = IP(IQ+1)
         DO 20 K=1,10
         F(IPS) = F(IPS)+W(KK+K)*Q(JQ+K)
   20    CONTINUE
CDIR$ IVDEP
C$DIR NO_RECURRENCE
*VDIR NODEP
*VOCL LOOP,NOVREC
         DO 30 K=1,10
         F(IP(JQ+K)) = F(IP(JQ+K))+W(KK+K)*PSS
   30    CONTINUE
         IS  = IC+JA
         IX  = IS+1
         IY  = IS+2
         IZ  = IS+3
         F(IS) = F(IS)-PA(IS)*W(KK+1)-PA(IX)*W(KK+2)
     1                -PA(IY)*W(KK+4)-PA(IZ)*W(KK+7)
         F(IX) = F(IX)-PA(IS)*W(KK+2)-PA(IX)*W(KK+3)
     1                -PA(IY)*W(KK+5)-PA(IZ)*W(KK+8)
         F(IY) = F(IY)-PA(IS)*W(KK+4)-PA(IX)*W(KK+5)
     1                -PA(IY)*W(KK+6)-PA(IZ)*W(KK+9)
         F(IZ) = F(IZ)-PA(IS)*W(KK+7)-PA(IX)*W(KK+8)
     1                -PA(IY)*W(KK+9)-PA(IZ)*W(KK+10)
      ELSE IF(IW.EQ.10 .AND. JA.EQ.JB) THEN
         PSS = Q(JQ+1)
         JPS = IP(JQ+1)
         DO 40 K=1,10
         F(JPS) = F(JPS)+W(KK+K)*Q(IQ+K)
   40    CONTINUE
CDIR$ IVDEP
C$DIR NO_RECURRENCE
*VDIR NODEP
*VOCL LOOP,NOVREC
         DO 50 K=1,10
         F(IP(IQ+K)) = F(IP(IQ+K))+W(KK+K)*PSS
   50    CONTINUE
         IS  = IC+JA
         IX  = INDX(IA+1)+JA
         IY  = INDX(IA+2)+JA
         IZ  = INDX(IA+3)+JA
         F(IS) = F(IS)-PA(IS)*W(KK+1)-PA(IX)*W(KK+2)
     1                -PA(IY)*W(KK+4)-PA(IZ)*W(KK+7)
         F(IX) = F(IX)-PA(IS)*W(KK+2)-PA(IX)*W(KK+3)
     1                -PA(IY)*W(KK+5)-PA(IZ)*W(KK+8)
         F(IY) = F(IY)-PA(IS)*W(KK+4)-PA(IX)*W(KK+5)
     1                -PA(IY)*W(KK+6)-PA(IZ)*W(KK+9)
         F(IZ) = F(IZ)-PA(IS)*W(KK+7)-PA(IX)*W(KK+8)
     1                -PA(IY)*W(KK+9)-PA(IZ)*W(KK+10)
      ELSE IF(IW.EQ.10 .AND. JW.EQ.10) THEN
         DO 80 I=1,10
         ID1 = KK+10*(I-1)
         ID2 = KK+I-10
         IPJ = IP(JQ+I)
         IPI = IP(IQ+I)
         DO 60 K=1,10
         F(IPJ) = F(IPJ)+W(ID2+10*K)*Q(IQ+K)
   60    CONTINUE
         DO 70 K=1,10
         F(IPI) = F(IPI)+W(ID1+K)*Q(JQ+K)
   70    CONTINUE
   80    CONTINUE
         DO 100 I=1,4
         IS  = INDX(IA+I-1)+JA
         IX  = IS+1
         IY  = IS+2
         IZ  = IS+3
         DO 90 K=1,4
         KS  = INDX(IA+K-1)+JA
         KX  = KS+1
         KY  = KS+2
         KZ  = KS+3
         IJK = KK+10*(IWW(I,K)-1)
         F(IS) = F(IS)-PA(KS)*W(IJK+1)-PA(KX)*W(IJK+2)
     1                -PA(KY)*W(IJK+4)-PA(KZ)*W(IJK+7)
         F(IX) = F(IX)-PA(KS)*W(IJK+2)-PA(KX)*W(IJK+3)
     1                -PA(KY)*W(IJK+5)-PA(KZ)*W(IJK+8)
         F(IY) = F(IY)-PA(KS)*W(IJK+4)-PA(KX)*W(IJK+5)
     1                -PA(KY)*W(IJK+6)-PA(KZ)*W(IJK+9)
         F(IZ) = F(IZ)-PA(KS)*W(IJK+7)-PA(KX)*W(IJK+8)
     1                -PA(KY)*W(IJK+9)-PA(KZ)*W(IJK+10)
   90    CONTINUE
  100    CONTINUE
C     GENERAL CODE - ALSO VALID FOR D ORBITALS.
      ELSE
         KKW    = 0
CDIRX NOVECTOR
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 140 I=IA,IB
         KA     = INDX(I)
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 130 J=IA,I
         KB     = INDX(J)
         IJ     = KA+J
         IEQJ   = I.EQ.J
         PIJ    = PA(IJ)+PB(IJ)
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 120 K=JA,JB
         KC     = INDX(K)
         IK     = KA+K
         JK     = KB+K
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 110 L=JA,K
         IL     = KA+L
         JL     = KB+L
         KL     = KC+L
         KEQL   = K.EQ.L
         KKW    = KKW+1
         A      = W(KK+KKW)
         PKL    = PA(KL)+PB(KL)
         IF(IEQJ.AND.KEQL) THEN
            F(IJ) =  F(IJ) + A*PKL
            F(KL) =  F(KL) + A*PIJ
            F(IK) =  F(IK) - A*PA(JL)
         ELSE IF(IEQJ) THEN
            F(KL) =  F(KL) + A*PIJ
            F(IJ) =  F(IJ) + A*PKL*TWO
            F(IK) =  F(IK) - A*PA(JL)
            F(IL) =  F(IL) - A*PA(JK)
         ELSE IF(KEQL) THEN
            F(IJ) =  F(IJ) + A*PKL
            F(KL) =  F(KL) + A*PIJ*TWO
            F(IK) =  F(IK) - A*PA(JL)
            F(JK) =  F(JK) - A*PA(IL)
         ELSE
            F(IJ) =  F(IJ) + A*PKL*TWO
            F(KL) =  F(KL) + A*PIJ*TWO
            F(IK) =  F(IK) - A*PA(JL)
            F(JK) =  F(JK) - A*PA(IL)
            F(IL) =  F(IL) - A*PA(JK)
            F(JL) =  F(JL) - A*PA(IK)
         ENDIF
  110    CONTINUE
  120    CONTINUE
  130    CONTINUE
  140    CONTINUE
CDIRX VECTOR
      ENDIF
      KK     = KK+JW*IW
  150 CONTINUE
  160 CONTINUE
      RETURN
C *** MINDO AND CNDO.
  200 KR     = 1
      DO 240 II=2,NUMAT
      IA     = NFIRST(II)
      IB     = NLAST(II)
      IMINUS = II-1
      DO 230 JJ=1,IMINUS
      KR     = KR+1
      A      = W(KR)
      JA     = NFIRST(JJ)
      JB     = NLAST(JJ)
      DO 220 I=IA,IB
      KK     = INDX(I)+I
      DO 210 K=JA,JB
      LL     = INDX(K)+K
      IK     = INDX(I)+K
      F(KK)  = F(KK)+A*(PA(LL)+PB(LL))
      F(LL)  = F(LL)+A*(PA(KK)+PB(KK))
      F(IK)  = F(IK)-A*PA(IK)
  210 CONTINUE
  220 CONTINUE
  230 CONTINUE
      KR     = KR+1
  240 CONTINUE
      RETURN
      END
