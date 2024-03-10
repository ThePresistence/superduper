      SUBROUTINE FOCKX (F,PA,PB,Q,W,LM4,LM6)
C     *
C     TWO-ELECTRON CONTRIBUTIONS TO MNDO-TYPE FOCK MATRIX.
C     SCALAR CODE FOR TWO-CENTER EXCHANGE CONTRIBUTIONS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     F(LM4)    FOCK MATRIX (I,O).
C     PA(LM4)   RHF OR UHF-ALPHA DENSITY MATRIX (I).
C     PB(LM4)   UHF-BETA DENSITY MATRIX (I).
C     Q(LM6)    SCRATCH ARRAY FOR ONE-CENTER PAIR TERMS (S).
C     W(LM6,*)  TWO-ELECTRON INTEGRALS (I).
C     *
      USE LIMIT, ONLY: LM1, LMI, LME, LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
     ./FINDX2/ JP1(LME),JP2(LME),JP3(LME),JX(LM1),JXLAST
     ./INDEX / INDX(LMX)
     ./INDEXW/ NW(LM1)
     ./UHF   / UHF
      DIMENSION F(LM4),PA(LM4),PB(LM4),Q(LM6),W(LM6,LM6)
      DIMENSION IWW(4,4)
      DATA IWW/1,2,4,7,2,3,5,8,4,5,6,9,7,8,9,10/
C
C *** COULOMB CONTRIBUTIONS.
C     ONE-CENTER EXCHANGE CONTRIBUTIONS ARE IMPLICITLY INCLUDED FOR RHF.
C$OMP PARALLEL DO SCHEDULE(STATIC)
      DO 10 KL=1,LM6
      Q(KL)  = (PA(IP(KL))+PB(IP(KL)))
      IF(IP1(KL).NE.IP2(KL)) Q(KL)=Q(KL)*TWO
   10 CONTINUE
C$OMP END PARALLEL DO

C$OMP PARALLEL DO SCHEDULE(STATIC)
C$OMP.PRIVATE(SUM)
      DO 20 IJ=1,LM6
      SUM    = ZERO
CDIR$ NOBL
      DO 15 KL=1,LM6
      SUM    = SUM+Q(KL)*W(KL,IJ)
   15 CONTINUE
      F(IP(IJ)) = F(IP(IJ))+SUM
   20 CONTINUE
C$OMP END PARALLEL DO
C
C *** TWO-CENTER EXCHANGE CONTRIBUTIONS.
C     OFFDIAGONAL TWO-CENTER TERMS (IJ,KL).
C$OMP PARALLEL DO SCHEDULE(DYNAMIC)
C$OMP.PRIVATE(A,IA,IB,IC,IJ,IJK,IJS,IJW,IK,IL,IORBS,IS,IW,IX,IY,IZ,
C$OMP.JA,JB,JK,JL,JORBS,JW,KA,KB,KL,KLS,KLW,KS,KX,KY,KZ)
      DO 55 II=1,NUMAT
      IA     = NFIRST(II)
      IB     = NLAST(II)
      IC     = INDX(IA)
      IORBS  = IB-IA+1
      IW     = INDX(IORBS)+IORBS
      IJ     = NW(II)-1
      DO 50 JJ=1,II-1
      JA     = NFIRST(JJ)
      JB     = NLAST(JJ)
      JORBS  = JB-JA+1
      JW     = INDX(JORBS)+JORBS
      KL     = NW(JJ)-1
      IF(IW.EQ.1 .AND. JW.EQ.1) THEN
         IS  = IC+JA
         F(IS) = F(IS)-PA(IS)*W(IJ+1,KL+1)
      ELSE IF(IW.EQ.1 .AND. JW.EQ.10) THEN
         IS  = IC+JA
         IX  = IS+1
         IY  = IS+2
         IZ  = IS+3
         IJS = IJ+1
         F(IS) = F(IS)-PA(IS)*W(KL+1,IJS)-PA(IX)*W(KL+2,IJS)
     1                -PA(IY)*W(KL+4,IJS)-PA(IZ)*W(KL+7,IJS)
         F(IX) = F(IX)-PA(IS)*W(KL+2,IJS)-PA(IX)*W(KL+3,IJS)
     1                -PA(IY)*W(KL+5,IJS)-PA(IZ)*W(KL+8,IJS)
         F(IY) = F(IY)-PA(IS)*W(KL+4,IJS)-PA(IX)*W(KL+5,IJS)
     1                -PA(IY)*W(KL+6,IJS)-PA(IZ)*W(KL+9,IJS)
         F(IZ) = F(IZ)-PA(IS)*W(KL+7,IJS)-PA(IX)*W(KL+8,IJS)
     1                -PA(IY)*W(KL+9,IJS)-PA(IZ)*W(KL+10,IJS)
      ELSE IF(IW.EQ.10 .AND. JW.EQ.1) THEN
         IS  = IC+JA
         IX  = INDX(IA+1)+JA
         IY  = INDX(IA+2)+JA
         IZ  = INDX(IA+3)+JA
         KLS = KL+1
         F(IS) = F(IS)-PA(IS)*W(IJ+1,KLS)-PA(IX)*W(IJ+2,KLS)
     1                -PA(IY)*W(IJ+4,KLS)-PA(IZ)*W(IJ+7,KLS)
         F(IX) = F(IX)-PA(IS)*W(IJ+2,KLS)-PA(IX)*W(IJ+3,KLS)
     1                -PA(IY)*W(IJ+5,KLS)-PA(IZ)*W(IJ+8,KLS)
         F(IY) = F(IY)-PA(IS)*W(IJ+4,KLS)-PA(IX)*W(IJ+5,KLS)
     1                -PA(IY)*W(IJ+6,KLS)-PA(IZ)*W(IJ+9,KLS)
         F(IZ) = F(IZ)-PA(IS)*W(IJ+7,KLS)-PA(IX)*W(IJ+8,KLS)
     1                -PA(IY)*W(IJ+9,KLS)-PA(IZ)*W(IJ+10,KLS)
      ELSE IF(IW.EQ.10 .AND. JW.EQ.10) THEN
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 40 I=1,4
         IS  = INDX(IA+I-1)+JA
         IX  = IS+1
         IY  = IS+2
         IZ  = IS+3
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 30 K=1,4
         KS  = INDX(IA+K-1)+JA
         KX  = KS+1
         KY  = KS+2
         KZ  = KS+3
         IJK = IJ+IWW(I,K)
         F(IS) = F(IS)-PA(KS)*W(KL+1,IJK)-PA(KX)*W(KL+2,IJK)
     1                -PA(KY)*W(KL+4,IJK)-PA(KZ)*W(KL+7,IJK)
         F(IX) = F(IX)-PA(KS)*W(KL+2,IJK)-PA(KX)*W(KL+3,IJK)
     1                -PA(KY)*W(KL+5,IJK)-PA(KZ)*W(KL+8,IJK)
         F(IY) = F(IY)-PA(KS)*W(KL+4,IJK)-PA(KX)*W(KL+5,IJK)
     1                -PA(KY)*W(KL+6,IJK)-PA(KZ)*W(KL+9,IJK)
         F(IZ) = F(IZ)-PA(KS)*W(KL+7,IJK)-PA(KX)*W(KL+8,IJK)
     1                -PA(KY)*W(KL+9,IJK)-PA(KZ)*W(KL+10,IJK)
   30    CONTINUE
   40    CONTINUE
      ELSE
C     GENERAL CODE   - ALSO VALID FOR D ORBITALS.
C     FIRST  SECTION - CONTRIBUTIONS FROM (II,KL).
C     SECOND SECTION - CONTRIBUTIONS FROM (IJ,KL) WITH I.NE.J.
         DO 43 I=IA,IB
         KA     = INDX(I)
         IJW    = IJ+INDX(I-IA+2)
         KLW    = KL
         DO 42 K=JA,JB
         IK     = KA+K
         DO 41 L=JA,K
         IL     = KA+L
         KLW    = KLW+1
         A      = W(KLW,IJW)
         F(IK) =  F(IK) - A*PA(IL)
         IF(K.NE.L) F(IL) =  F(IL) - A*PA(IK)
   41    CONTINUE
   42    CONTINUE
   43    CONTINUE
         DO 47 I=IA+1,IB
         KA     = INDX(I)
         DO 46 J=IA,I-1
         KB     = INDX(J)
         IJW    = IJ+INDX(I-IA+1)+J-IA+1
         KLW    = KL
         DO 45 K=JA,JB
         IK     = KA+K
         JK     = KB+K
         DO 44 L=JA,K
         IL     = KA+L
         JL     = KB+L
         KLW    = KLW+1
         A      = W(KLW,IJW)
         F(IK)  = F(IK) - A*PA(JL)
         F(JK)  = F(JK) - A*PA(IL)
         IF(K.NE.L) THEN
            F(IL)  = F(IL) - A*PA(JK)
            F(JL)  = F(JL) - A*PA(IK)
         ENDIF
   44    CONTINUE
   45    CONTINUE
   46    CONTINUE
   47    CONTINUE
      ENDIF
   50 CONTINUE
   55 CONTINUE
C$OMP END PARALLEL DO
      IF(.NOT.UHF) RETURN
C *** ONE-CENTER EXCHANGE CONTRIBUTIONS FOR UHF.
C     OFFDIAGONAL ONE-CENTER TERMS (II,KK) FOR SP-BASIS.
CDIR$ IVDEP
C$DIR NO_RECURRENCE
*VDIR NODEP
*VOCL LOOP,NOVREC
      DO 60 I=1,NUMAT
      IA     = NFIRST(I)
      IORBS  = NLAST(I)-IA+1
      IF(IORBS.EQ.4) THEN
         IJ     = NW(I)-1
         IXS    = INDX(IA+1)+IA
         IYS    = INDX(IA+2)+IA
         IZS    = INDX(IA+3)+IA
         IYX    = IYS+1
         IZX    = IZS+1
         IZY    = IZS+2
         F(IXS) = F(IXS)-PA(IXS)*W(IJ+3,IJ+1)
         F(IYS) = F(IYS)-PA(IYS)*W(IJ+6,IJ+1)
         F(IZS) = F(IZS)-PA(IZS)*W(IJ+10,IJ+1)
         F(IYX) = F(IYX)-PA(IYX)*W(IJ+6,IJ+3)
         F(IZX) = F(IZX)-PA(IZX)*W(IJ+10,IJ+3)
         F(IZY) = F(IZY)-PA(IZY)*W(IJ+10,IJ+6)
      ENDIF
   60 CONTINUE
C     DIAGONAL ONE-CENTER TERMS (IJ,IJ), GENERAL CODE.
CDIR$ IVDEP
C$DIR NO_RECURRENCE
*VDIR NODEP
*VOCL LOOP,NOVREC
C$OMP PARALLEL DO SCHEDULE(STATIC)
      DO 70 IJ=1,LM6
      F(IP(IJ)) = F(IP(IJ)) - PA(IP(IJ))*W(IJ,IJ)
   70 CONTINUE
C$OMP END PARALLEL DO
      DO 80 IJ=1,LM6
      I      = IP1(IJ)
      J      = IP2(IJ)
      IF(I.NE.J) THEN
         II    = INDX(I)+I
         JJ    = INDX(J)+J
         F(II) = F(II)-PA(JJ)*W(IJ,IJ)
         F(JJ) = F(JJ)-PA(II)*W(IJ,IJ)
      ENDIF
   80 CONTINUE
C     ALL ONE-CENTER TERMS ARE NOW INCLUDED FOR AN SP-BASIS.
C     THE FOLLOWING CODE COVERS ANY OFFDIAGONAL TERMS FOR AN SPD-BASIS.
      DO 110 N=1,NUMAT
      IA     = NFIRST(N)
      IORBS  = NLAST(N)-IA+1
      IF(IORBS.LE.4) GO TO 110
      KLMIN  = JX(N)
      KLMAX  = KLMIN+IORBS*IORBS-1
      KSTART = JP1(KLMIN)
      DO 100 KL=KLMIN,KLMAX
      K      = JP1(KL)
      L      = JP2(KL)
C     ONLY ELEMENTS F(IK) WITH I.GE.K ARE NEEDED. THEREFORE,
C     THE LOOP OVER IJ CAN START AT IJMIN.GE.KLMIN.
      IJMIN  = KLMIN+IORBS*(K-KSTART)
      DO 90 IJ=IJMIN,KLMAX
C     DIAGONAL TERMS HAVE BEEN INCLUDED ABOVE ALSO FOR AN SPD-BASIS.
      IF(JP3(IJ).EQ.JP3(KL)) GO TO 90
      WIJKL  = W(JP3(IJ),JP3(KL))
      IF(WIJKL.NE.ZERO) THEN
         I   = JP1(IJ)
         J   = JP2(IJ)
         IK  = INDX(I)+K
         JL  = INDX(MAX(J,L))+MIN(J,L)
         F(IK) = F(IK)-PA(JL)*WIJKL
      ENDIF
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE
      RETURN
      END
