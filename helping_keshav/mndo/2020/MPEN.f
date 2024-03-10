      SUBROUTINE MPEN (E2,EIG,G,W,IMOCI,IPROD,NSYM,LM3,LM9,LM10,NB2,
     1                 IA,IB,JB,IOUT2,NNTOT,NUMB)
C     *
C     PERTURBATION TREATMENT WITH ONE MAIN CONFIGURATION.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL DEBUG
      CHARACTER*4 NAME(4)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DENERG/ DENER
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
      DIMENSION EIG(LM3),G(LM10,LM10),W(LM9)
      DIMENSION IMOCI(LM3),IPROD(8,8),NSYM(LM3)
      DIMENSION IPERT(4)
      DATA NAME/'RSMP','RSEN','BWMP','BWEN'/
      DATA ONE5/1.5D0/
      DATA EBWLIM/1.0D-08/
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** INPUT OPTIONS.
      MPERT  = IN2(135)
      IF(MPERT.LE.0) THEN
         IPERT(1) = 0
         IPERT(2) = 0
         IPERT(3) = 0
         IPERT(4) = 1
      ELSE
         IPERT(1) = MOD(MPERT/1000,10)
         IPERT(2) = MOD(MPERT/100,10)
         IPERT(3) = MOD(MPERT/10,10)
         IPERT(4) = MOD(MPERT,10)
      ENDIF
C *** INITIALIZE ENERGIES.
      E2     = ZERO
      EB     = ZERO
      EC     = ZERO
      IF(IOUT2.GE.5) WRITE(NB6,620)
C *** LOOP OVER OPTIONS.
      DO 100 MOP=1,4
      IF(IPERT(MOP).EQ.0) GO TO 100
      EBW    = ZERO
      IBW    = 0
   10 IBW    = IBW+1
      DEBUG  = IOUT2.GE.5 .AND. IBW.EQ.1
      IF(DEBUG) WRITE(NB6,650)
C     READ INTEGRALS.
      IF(NNTOT.GT.LM9) THEN
         REWIND NB2
         READ(NB2) W
      ENDIF
C *** COMPUTE THE SECOND-ORDER ENERGIES.
      KK     = 0
      CC     = ZERO
      CCSUM  = ZERO
      EESUM  = ZERO
C     OUTER IJ-LOOP.
      DO 80 I=IA,IB
      II     = IMOCI(I)
      IS     = NSYM(II)
      DI     = EIG(II)
      DO 70 J=1,JB
      IJ     = IMOCI(J)
      JS     = NSYM(IJ)
      NSIJ   = IPROD(IS,JS)
      DIJ    = DI-EIG(IJ)
      GIJ    = G(I,J)
      GJI    = G(J,I)
C     INNER KL-LOOP.
      DO 60 K=IA,I
      IK     = IMOCI(K)
      KS     = NSYM(IK)
      DIJK   = DIJ+EIG(IK)
      GIK    = G(I,K)
      GJK    = G(J,K)
      GKI    = G(K,I)
      GKJ    = G(K,J)
      DO 50 L=1,J
      IL     = IMOCI(L)
      LS     = NSYM(IL)
      NSKL   = IPROD(KS,LS)
      IF(NSIJ.NE.NSKL) GO TO 50
      GIL    = G(I,L)
      GJL    = G(J,L)
      GKL    = G(K,L)
      GLI    = G(L,I)
      GLJ    = G(L,J)
      GLK    = G(L,K)
C     INTEGRALS B1 AND B2.
      KK     = KK+1
      IF(KK.GT.LM9) THEN
         READ(NB2) W
         KK  = 1
      ENDIF
      B1     = W(KK)
      BA     = B1*B1
      IF(I.GT.K .AND. J.GT.L) THEN
         KK  = KK+1
         IF(KK.GT.LM9) THEN
            READ(NB2) W
            KK = 1
         ENDIF
         B2  = W(KK)
         BB  = B2*B2
         BM  = (B1-B2)**2
      ENDIF
C     ENERGY DENOMINATOR.
      DIJKL  = DIJK-EIG(IL)
C     MP TREATMENT.
      IF(MOP.EQ.1 .OR. MOP.EQ.3) THEN
         EA  = DIJKL-EBW
         EE  = -BA/EA
         IF(I.GT.K .OR.  J.GT.L) EE=EE*TWO
         IF(I.GT.K .AND. J.GT.L) EE=EE-TWO*(BB+BM)/EA
         CC  = -EE/EA
C     EN TREATMENT.
      ELSE IF(MOP.EQ.2 .OR. MOP.EQ.4) THEN
         EN0    = DIJKL+GIK+GJL-GIJ-GIL-GKJ-GKL-EBW
         IF(I.EQ.K .OR. J.EQ.L) THEN
            EA  = EN0+GJI+GLK
            IF(I.GT.K) EA=EA+GKI
            IF(J.GT.L) EA=EA+GLJ
            EE  = -BA/EA
            IF(I.GT.K .OR. J.GT.L) EE=EE*TWO
            CC  = -EE/EA
         ELSE
            GS  = GJI+GLI+GJK+GLK
            EA  = EN0+GKI+GLJ+PT5*GS
            EB  = EN0-GKI-GLJ+ONE5*GS
            BP  = (B1+B2)**2
            EE  = -BP/EA-THREE*BM/EB
            CC  = BP/EA**2+THREE*BM/EB**2
         ENDIF
      ENDIF
C     SUMMATION.
      EESUM  = EESUM+EE
      CCSUM  = CCSUM+CC
C     DEBUG PRINTING.
      IF(DEBUG) THEN
         IF(I.EQ.K .OR. J.EQ.L) THEN
            B2 = ZERO
            EB = ZERO
            EC = ZERO
         ENDIF
         WRITE(NB6,630) II,IJ,IK,IL,B1,B2,EA,EB,EC,EE,CC
      ENDIF
   50 CONTINUE
   60 CONTINUE
   70 CONTINUE
   80 CONTINUE
C *** BW TREATMENT.
      IF(MOP.GT.2) THEN
         IF(IOUT2.GE.0) THEN
            IF(IBW.EQ.1) WRITE(NB6,650)
            WRITE(NB6,710) IBW,EESUM
         ENDIF
         IF(ABS(EESUM-EBW).GE.EBWLIM) THEN
            EBW = EESUM
            GO TO 10
         ENDIF
      ENDIF
C *** PRINT THE RESULTS.
      EEV    = EESUM*EV
      IF(IPERT(MOP).EQ.1) E2=EEV
      IF(IOUT2.LE.-5) GO TO 100
      EKCAL = EEV*EVCAL
      WRITE(NB6,510) EESUM,EEV,EKCAL,NAME(MOP)
      CC0   = ONE/(ONE+CCSUM)
      C0    = SQRT(CC0)
C     DAVIDSON AND POPLE CORRECTION FOR MOP=4.
      IF(MOP.EQ.4) THEN
         CCD   = ONE-CC0
         EDSUM = CCD*EESUM
         EDEV  = CCD*EEV
         EDCAL = CCD*EKCAL
         WRITE(NB6,511) EDSUM,EDEV,EDCAL,NAME(MOP)
         FN    = NUMB*2
         THET  = 2.0D0*ACOS(C0)
         DNOM  = 2.0D0/COS(THET)-2.0D0
         TANSQ = TAN(THET)**2
         TANSQ = SQRT(FN**2+2.0D0*FN*TANSQ)-FN
         TANSQ = TANSQ/DNOM-1.0D0
         EPSUM = TANSQ*EESUM
         EPEV  = TANSQ*EEV
         EPCAL = TANSQ*EKCAL
         WRITE(NB6,512) EPSUM,EPEV,EPCAL,NAME(MOP)
         HDCAL = DENER+EKCAL
         HDEV  = HDCAL/EVCAL
         HDSUM = HDEV/EV
         WRITE(NB6,520) HDSUM,HDEV,HDCAL,NAME(MOP)
         HCAL1 = HDCAL+EDCAL
         HDEV  = HCAL1/EVCAL
         HDSUM = HDEV/EV
         WRITE(NB6,521) HDSUM,HDEV,HCAL1,NAME(MOP)
         HCAL2 = HDCAL+EPCAL
         HDEV  = HCAL2/EVCAL
         HDSUM = HDEV/EV
         WRITE(NB6,522) HDSUM,HDEV,HCAL2,NAME(MOP)
      ENDIF
      WRITE(NB6,530) C0,CC0
  100 CONTINUE
      RETURN
  510 FORMAT (/ 1X,'SECOND-ORDER ENERGY =',F12.7,' A.U.=',F12.7,' EV=',
     1             F12.7,' KCAL/MOL',10X,'OPTION=',A4)
  511 FORMAT (  1X,'DAVIDSON CORRECTION =',F12.7,' A.U.=',F12.7,' EV=',
     1             F12.7,' KCAL/MOL',10X,'OPTION=',A4)
  512 FORMAT (  1X,'POPLE CORRECTION    =',F12.7,' A.U.=',F12.7,' EV=',
     1             F12.7,' KCAL/MOL',10X,'OPTION=',A4)
  520 FORMAT (/ 1X,'HEAT OF FORMATION   =',F12.7,' A.U.=',F12.7,' EV=',
     1             F12.7,' KCAL/MOL',10X,'OPTION=',A4)
  521 FORMAT (  1X,'HEAT OF FORMATION   =',F12.7,' A.U.=',F12.7,' EV=',
     1             F12.7,' KCAL/MOL',10X,'OPTION=',A4,'+DAVIDSON')
  522 FORMAT (  1X,'HEAT OF FORMATION   =',F12.7,' A.U.=',F12.7,' EV=',
     1             F12.7,' KCAL/MOL',10X,'OPTION=',A4,'+POPLE')
  530 FORMAT (/ 1X,'COEFFICIENT OF SCF CONFIGURATION       ',F12.7/
     1          1X,'SQUARE OF COEFFICIENT                  ',F12.7)
  620 FORMAT (//1X,'DEBUG PRINT FOR SECOND-ORDER TERMS.'//
     1          1X,'  II  IJ  IK  IL',5X,'B1',10X,'B2',10X,'EA',
     2         10X,'EB',10X,'EC',10X,'EE',10X,'CC')
  630 FORMAT (  1X,4I4,8F12.7)
  650 FORMAT (  1X)
  710 FORMAT (  1X,'CYCLE',I3, 5X,'ENERGY =',F12.7,' A.U.')
      END
