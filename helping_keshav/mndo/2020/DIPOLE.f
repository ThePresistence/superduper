      SUBROUTINE DIPOLE (PA,PB,LM4,PRT)
C     *
C     POPULATION ANALYSIS, AND CALCULATION OF DIPOLE MOMENT.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     PA(LM4)   RHF OR UHF-ALPHA DENSITY MATRIX (I).
C     PB(LM4)   UHF-BETA DENSITY MATRIX (I).
C     PRT       LOGICAL PRINTING FLAG (I).
C     *
      USE LIMIT, ONLY: LM1, LMX, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2 ELEMNT
      LOGICAL PRT
      COMMON
     ./AMASS / AMS(LM1)
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CHARGE/ QI(LM1)
     ./CHARGP/ POP(3,LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./CPEOPT/ NFLCPE,NPTCPE,NCPEZ
     ./CPEOUT/ CPEENE,CPEDIP(3)
     ./DIPOL / DM,DMX(3),DMP(3)
     ./DIPOL1/ HYF(LMZ),HYFPD(LMZ)
     ./ELEMTS/ ELEMNT(107)
     ./INDEX / INDX(LMX)
     ./INFOB / PC(3,LM1),PM(3),PR(3,3)
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./PARDER/ CORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
      DIMENSION PA(LM4),PB(LM4)
      DIMENSION CG(3),DM1(4),DM2(4),PAO(9)
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** INPUT OPTIONS.
      KHARGE = IN2(65)
      NPRINT = IN2(72)
C *** INITIALIZATION.
      DM     = ZERO
      DO 10 I=1,3
      CG(I)  = ZERO
      DMX(I) = ZERO
      DMP(I) = ZERO
   10 CONTINUE
      DO 20 I=1,4
      DM1(I) = ZERO
      DM2(I) = ZERO
   20 CONTINUE
      IF(PRT) THEN
         WRITE(NB6,500)
         ITMAX = 0
         DO 30 I=1,NUMAT
         ITMAX = MAX(ITMAX,NLAST(I)-NFIRST(I)+1)
   30    CONTINUE
         IF(ITMAX.LE.1) THEN
            WRITE(NB6,510)
         ELSE IF(ITMAX.LE.4) THEN
            WRITE(NB6,520)
         ELSE
            WRITE(NB6,525)
         ENDIF
      ENDIF
C *** DETERMINE CENTER-OF-GRAVITY COORDINATES FOR IONS.
      IF(KHARGE.NE.0) THEN
         AMSUM  = ZERO
         DO 35 I=1,NUMAT
         AMSI   = AMS(I)
         AMSUM  = AMSUM+AMSI
         CG(1)  = CG(1)+AMSI*COORD(1,I)
         CG(2)  = CG(2)+AMSI*COORD(2,I)
         CG(3)  = CG(3)+AMSI*COORD(3,I)
   35    CONTINUE
         RAMSUM = ONE/AMSUM
         CG(1)  = CG(1)*RAMSUM
         CG(2)  = CG(2)*RAMSUM
         CG(3)  = CG(3)*RAMSUM
      ENDIF
C *** LOOP OVER ALL ATOMS.
      QISUM  = ZERO
      DO 70 I=1,NUMAT
      NI     = NAT(I)
      IA     = NFIRST(I)
      IB     = NLAST(I)
      IORBS  = IB-IA+1
C *** CALCULATE ATOMIC POPULATIONS.
      L      = 0
      W      = ZERO
      DO 40 J=IA,IB
      K      = INDX(J)+J
      L      = L+1
      P      = PA(K)+PB(K)
      PAO(L) = P
      W      = W+P
   40 CONTINUE
      QI(I)  = CORE(NI)-W
      QISUM  = QISUM + QI(I)
C *** CALCULATE SUBSHELL POPULATIONS.
      POP(1,I) = PAO(1)
      POP(2,I) = ZERO
      POP(3,I) = ZERO
      IF(L.GT.1) THEN
         POP(2,I) = PAO(2)+PAO(3)+PAO(4)
         IF(L.GT.4) THEN
            POP(3,I) = PAO(5)+PAO(6)+PAO(7)+PAO(8)+PAO(9)
         ENDIF
      ENDIF
C *** PRINT ATOMIC AND ORBITAL POPULATIONS.
      IF(PRT) THEN
         IF(L.LE.4) THEN
            WRITE(NB6,530) ELEMNT(NI),I,QI(I),W,(PAO(J),J=1,L)
         ELSE
            WRITE(NB6,530) ELEMNT(NI),I,QI(I),W,(PAO(J),J=1,4),POP(3,I)
         ENDIF
      ENDIF
C *** CHARGE CONTRIBUTIONS TO DIPOLE MOMENT.
C     CONVERSION FACTOR : 1 e  =  4.803242D-10 esu
C     TAKEN FROM J.PHYS.CHEM.REF.DATA 2, 663 (1973).
      DO 50 J=1,3
      DM1(J) = DM1(J)+4.80324D0*QI(I)*(COORD(J,I)-CG(J))
   50 CONTINUE
C *** SP HYBRIDIZATION CONTRIBUTIONS TO DIPOLE MOMENT.
      IF(IORBS.EQ.1) GO TO 70
      DO 60 J=1,3
      L      = INDX(IA+J)+IA
      DM2(J) = DM2(J)-HYF(NI)*(PA(L)+PB(L))
   60 CONTINUE
C *** PD HYBRIDIZATION CONTRIBUTIONS TO DIPOLE MOMENT.
      IF(IORBS.EQ.4) GO TO 70
      XT     = ONE/SQRT(THREE)
C     REQUIRED SUMS OF DENSITY MATRIX ELEMENTS.
C     DX     = P(XZ,Z) + P(X2Y2,X) + P(XY,Y) - ONE/SQRT(THREE)*P(Z2,X)
C     DY     = P(YZ,Z) - P(X2Y2,Y) + P(XY,X) - ONE/SQRT(THREE)*P(Z2,Y)
C     DZ     = P(XZ,X) +             P(YZ,Y) + TWO/SQRT(THREE)*P(Z2,Z)
      LL     = INDX(IA+5) + IA+3
      DX     = PA(LL)+PB(LL)
      LL     = INDX(IA+4) + IA+1
      DX     = DX+PA(LL)+PB(LL)
      LL     = INDX(IA+8) + IA+2
      DX     = DX+PA(LL)+PB(LL)
      LL     = INDX(IA+6) + IA+1
      DX     = DX-XT*(PA(LL)+PB(LL))
      LL     = INDX(IA+7) + IA+3
      DY     = PA(LL)+PB(LL)
      LL     = INDX(IA+4) + IA+2
      DY     = DY-(PA(LL)+PB(LL))
      LL     = INDX(IA+8) + IA+1
      DY     = DY+PA(LL)+PB(LL)
      LL     = INDX(IA+6) + IA+2
      DY     = DY-XT*(PA(LL)+PB(LL))
      LL     = INDX(IA+5) + IA+1
      DZ     = PA(LL)+PB(LL)
      LL     = INDX(IA+7) + IA+2
      DZ     = DZ+PA(LL)+PB(LL)
      LL     = INDX(IA+6) + IA+3
      DZ     = DZ+TWO*XT*(PA(LL)+PB(LL))
      DM2(1) = DM2(1) - DX*HYFPD(NI)
      DM2(2) = DM2(2) - DY*HYFPD(NI)
      DM2(3) = DM2(3) - DZ*HYFPD(NI)
   70 CONTINUE
C *** CHECK FOR CORRECT OVERALL CHARGE.
      IF(PRT .AND. NINT(QISUM).NE.KHARGE) THEN
         WRITE(NB6,535) QISUM,KHARGE
      ENDIF
C *** CALCULATE DIPOLE MOMENT.
      IF(NFLCPE.GT.0) THEN
         DO 80 J=1,3
         DMX(J) = DM1(J)+DM2(J)+CPEDIP(J)
   80    CONTINUE
      ELSE
         DO 90 J=1,3
         DMX(J) = DM1(J)+DM2(J)
   90    CONTINUE
      END IF
      CPETDP = SQRT(CPEDIP(1)**2+CPEDIP(2)**2+CPEDIP(3)**2)
      DM1(4) = SQRT(DM1(1)**2+DM1(2)**2+DM1(3)**2)
      DM2(4) = SQRT(DM2(1)**2+DM2(2)**2+DM2(3)**2)
      DM     = SQRT(DMX(1)**2+DMX(2)**2+DMX(3)**2)
C *** PRINTING SECTION.
      IF(PRT) THEN
         WRITE(NB6,540)
         WRITE(NB6,550) (DM1(J),J=1,4)
         WRITE(NB6,560) (DM2(J),J=1,4)
         IF(NFLCPE.GT.0) THEN
            WRITE(6,565) (CPEDIP(J),J=1,3),CPETDP
         END IF
         WRITE(NB6,570) (DMX(J),J=1,3),DM
      ENDIF
C *** TRANSFORM DIPOLE MOMENT TO PRINCIPAL AXES.
C     COMPUTE MOMENTS OF INERTIA AND TRANSFORMATION MATRIX.
      LPRINT = -5
      IF(PRT .AND. NPRINT.GT.0) LPRINT = 0
      IF(PRT .AND. NPRINT.GT.1) LPRINT = NPRINT
      CALL INERT (AMS,COORD,PC,PM,PR,ITYPE,NUMAT,NAT,LPRINT)
C     PERFORM THE TRANSFORMATION.
C     CALL DGEMV ('T',3,3,ONE,PR,3,DMX,1,ONE,DMP,1)
      DO 110 J=1,3
      DO 100 I=1,3
      DMP(J) = DMP(J) + PR(I,J)*DMX(I)
  100 CONTINUE
  110 CONTINUE
C     PRINTING SECTION.
      IF(PRT .AND. NPRINT.GT.0) THEN
         WRITE(NB6,580)
         WRITE(NB6,590) (DMP(J),J=1,3),DM
      ENDIF
      RETURN
  500 FORMAT(///5X,'NET ATOMIC CHARGES AND ORBITAL POPULATIONS.')
  510 FORMAT(/  5X,'ATOM',4X,'I',8X,'CHARGE',6X,'DENSITY',8X,'S'/)
  520 FORMAT(/  5X,'ATOM',4X,'I',8X,'CHARGE',6X,'DENSITY',8X,'S',
     1          10X,'PX',10X,'PY',10X,'PZ'/)
  525 FORMAT(/  5X,'ATOM',4X,'I',8X,'CHARGE',6X,'DENSITY',8X,'S',
     1          10X,'PX',10X,'PY',10X,'PZ',10X,'DSUM'/)
  530 FORMAT(   6X,A2,3X,I3,3X,7F12.5)
  535 FORMAT(/  5X,'TOTAL CHARGE',F12.5,'     DEVIATES FROM INPUT:',I5)
  540 FORMAT(///5X,'DIPOLE',14X,'X',11X,'Y',11X,'Z',9X,'TOTAL')
  550 FORMAT(/  5X,'POINT-CHARGE',4F12.5)
  560 FORMAT(   5X,'HYBRID      ',4F12.5)
  565 FORMAT(   5X,'CPE RESPONSE',4F12.5)
  570 FORMAT(   5X,'SUM         ',4F12.5)
  580 FORMAT(// 5X,'DIPOLE',14X,'A',11X,'B',11X,'C',9X,'TOTAL')
  590 FORMAT(/  5X,'PRINCIPAL AXIS',F10.5,3F12.5)
      END
