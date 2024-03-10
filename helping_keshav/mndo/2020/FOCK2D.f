      SUBROUTINE FOCK2D (EF,F,PA,PB,W,LM4,LM9,NFIRST,NLAST)
C     *
C     TWO-CENTER TWO-ELECTRON CONTRIBUTIONS TO THE FOCK MATRIX.
C     SPECIAL VERSION FOR COMPUTING THE TWO-ELECTRON CONTRIBUTION
C     OF A GIVEN ATOM PAIR TO THE GRADIENT.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     EF        ENERGY CONTRIBUTION (O).
C     F(LM4)    FOCK MATRIX CONTRIBUTION (S).
C     PA(LM4)   RHF OR UHF-ALPHA DENSITY MATRIX (I).
C     PB(LM4)   UHF-BETA DENSITY MATRIX (I).
C     W(LM9)    TWO-ELECTRON INTEGRALS (I).
C     NFIRST()  INDEX OF FIRST ORBITAL AT A GIVEN ATOM (I).
C     NLAST()   INDEX OF LAST  ORBITAL AT A GIVEN ATOM (I).
C     *
      USE LIMIT, ONLY: LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IEQJ,KEQL
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./INDEX / INDX(LMX)
     ./INOPT2/ IN2(300)
      DIMENSION F(LM4),PA(LM4),PB(LM4),W(LM9)
      DIMENSION NFIRST(2),NLAST(2)
      DIMENSION PP(32)
      DIMENSION IWW(4,4),IXH(10),IXX(10)
      DATA IWW/1,2,4,7,2,3,5,8,4,5,6,9,7,8,9,10/
      DATA IXH/3,5,6,8,9,10,12,13,14,15/
      DATA IXX/15,20,21,26,27,28,33,34,35,36/
C *** INITIALIZATION.
      IA     = NFIRST(2)
      IB     = NLAST(2)
      JB     = NLAST(1)
      IF(IN2(2).GT.0) GO TO 200
C *** MNDO.
      IF(IB.EQ.2) THEN
         F(1)  = W(1)*(PA(3)+PB(3))
         F(2)  =-W(1)* PA(2)*PT5
         F(3)  = W(1)*(PA(1)+PB(1))
         EF    = EF+F(1)*PA(1)+F(2)*PA(2)+F(3)*PA(3)
         RETURN
      ENDIF
      DO 5  I=1,LM4
      F(I)   = ZERO
    5 CONTINUE
      IF(IB.EQ.5 .AND. JB.EQ.1) THEN
         PSS = PA(1)+PB(1)
         DO 10 I=1,10
         F(1)  = F(1)+W(I)*(PA(IXH(I))+PB(IXH(I)))
   10    CONTINUE
CDIR$ IVDEP
C$DIR NO_RECURRENCE
*VDIR NODEP
*VOCL LOOP,NOVREC
         DO 20 I=1,10
         F(IXH(I)) = F(IXH(I))+W(I)*PSS
   20    CONTINUE
         PA2   = -PA(2)*PT5
         PA4   = -PA(4)*PT5
         PA7   = -PA(7)*PT5
         PA11  = -PA(11)*PT5
         F(2)  = W(1)*PA2+W(2)*PA4+W(4)*PA7+W(7)*PA11
         F(4)  = W(2)*PA2+W(3)*PA4+W(5)*PA7+W(8)*PA11
         F(7)  = W(4)*PA2+W(5)*PA4+W(6)*PA7+W(9)*PA11
         F(11) = W(7)*PA2+W(8)*PA4+W(9)*PA7+W(10)*PA11
      ELSE IF(IB.EQ.5 .AND. JB.EQ.4) THEN
         F15 = ZERO
         PSS = PA(15)+PB(15)
         DO 30 I=1,10
         F15   = F15+W(I)*(PA(I)+PB(I))
         F(I)  = F(I)+W(I)*PSS
   30    CONTINUE
         F(15) = F15
         PA11  = -PA(11)*PT5
         PA12  = -PA(12)*PT5
         PA13  = -PA(13)*PT5
         PA14  = -PA(14)*PT5
         F(11) = W(1)*PA11+W(2)*PA12+W(4)*PA13+W(7)*PA14
         F(12) = W(2)*PA11+W(3)*PA12+W(5)*PA13+W(8)*PA14
         F(13) = W(4)*PA11+W(5)*PA12+W(6)*PA13+W(9)*PA14
         F(14) = W(7)*PA11+W(8)*PA12+W(9)*PA13+W(10)*PA14
      ELSE IF(IB.EQ.8 .AND. JB.EQ.4) THEN
CDIRX SHORTLOOP
C$DIR MAX_TRIPS(64)
*VDIR LOOPCNT=64
*VOCL LOOP,REPEAT(63)
         DO 40 I=1,10
         PP(I)  = PA(IXX(I))+PB(IXX(I))
   40    CONTINUE
         DO 60 I=1,10
         ID1    = 10*(I-1)
         ID2    = I-10
         FI1    = ZERO
         FI2    = ZERO
         DO 50 J=1,10
         FI1    = FI1+W(ID2+10*J)*PP(J)
         FI2    = FI2+W(ID1+J)*(PA(J)+PB(J))
   50    CONTINUE
         F(I)   = FI1
         F(IXX(I)) = FI2
   60    CONTINUE
CDIRX SHORTLOOP
C$DIR MAX_TRIPS(64)
*VDIR LOOPCNT=64
*VOCL LOOP,REPEAT(63)
         DO 70 I=11,32
         PP(I)  = -PA(I)*PT5
   70    CONTINUE
CDIRX NOVECTOR
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 90 I=1,4
         ID     = INDX(I+4)
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 80 K=1,4
         KD     = INDX(K+4)
         IK     = 10*(IWW(I,K)-1)
         F(ID+1)= F(ID+1)+PP(KD+1)*W(IK+1)+PP(KD+2)*W(IK+2)
     1                   +PP(KD+3)*W(IK+4)+PP(KD+4)*W(IK+7)
         F(ID+2)= F(ID+2)+PP(KD+1)*W(IK+2)+PP(KD+2)*W(IK+3)
     1                   +PP(KD+3)*W(IK+5)+PP(KD+4)*W(IK+8)
         F(ID+3)= F(ID+3)+PP(KD+1)*W(IK+4)+PP(KD+2)*W(IK+5)
     1                   +PP(KD+3)*W(IK+6)+PP(KD+4)*W(IK+9)
         F(ID+4)= F(ID+4)+PP(KD+1)*W(IK+7)+PP(KD+2)*W(IK+8)
     1                   +PP(KD+3)*W(IK+9)+PP(KD+4)*W(IK+10)
   80    CONTINUE
   90    CONTINUE
C     GENERAL CODE - ALSO VALID FOR D ORBITALS.
      ELSE
         KK     = 0
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
         DO 120 K=1,JB
         KC     = INDX(K)
         IK     = KA+K
         JK     = KB+K
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 110 L=1,K
         IL     = KA+L
         JL     = KB+L
         KL     = KC+L
         KEQL   = K.EQ.L
         KK     = KK+1
         A      = W(KK)
         F(IJ)  = F(IJ) + A*(PA(KL)+PB(KL))
         F(KL)  = F(KL) + A*PIJ
         A      = A*PT5
         IF(IEQJ.AND.KEQL) THEN
            F(IK) = F(IK) - A*PA(JL)
         ELSE IF(IEQJ) THEN
            F(IK) = F(IK) - A*PA(JL)
            F(IL) = F(IL) - A*PA(JK)
         ELSE IF(KEQL) THEN
            F(IK) = F(IK) - A*PA(JL)
            F(JK) = F(JK) - A*PA(IL)
         ELSE
            F(IK) = F(IK) - A*PA(JL)
            F(JK) = F(JK) - A*PA(IL)
            F(IL) = F(IL) - A*PA(JK)
            F(JL) = F(JL) - A*PA(IK)
         ENDIF
  110    CONTINUE
  120    CONTINUE
  130    CONTINUE
  140    CONTINUE
      ENDIF
CDIRX VECTOR
      DO 150 I=1,LM4
      EF     = EF+F(I)*PA(I)
  150 CONTINUE
      RETURN
C *** MINDO AND CNDO.
  200 A      = W(2)
      DO 210 I=1,LM4
      F(I)   = ZERO
  210 CONTINUE
      DO 230 I=IA,IB
      KA     = INDX(I)
      KK     = KA+I
      DO 220 K=1,JB
      LL     = INDX(K+1)
      IK     = KA+K
      F(KK)  = F(KK)+A*(PA(LL)+PB(LL))
      F(LL)  = F(LL)+A*(PA(KK)+PB(KK))
      F(IK)  = F(IK)-A* PA(IK)*PT5
  220 CONTINUE
  230 CONTINUE
      DO 240 I=1,LM4
      EF     = EF+F(I)*PA(I)
  240 CONTINUE
      RETURN
      END
