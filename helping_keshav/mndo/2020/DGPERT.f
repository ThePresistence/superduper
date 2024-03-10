      SUBROUTINE DGPERT (E2,EIG,G,W,IMOCI,IPROD,NSYM,LM3,LM9,LM10,NB2,
     1                   IA,IB,JB,IOUT2,NNTOT,CA,CB,DELT,KK,LL,IABSCI)
C     *
C     PERTURBATION TREATMENT WITH TWO MAIN CONFIGURATIONS.
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
C *** INITIALIZATION.
      E2     = ZERO
      I4     = 0
      KKS    = NSYM(KK)
      LLS    = NSYM(LL)
      DO 10 I=1,IB
      II     = IMOCI(I)
      IF(II.EQ.KK) KMO=I
      IF(II.EQ.LL) LMO=I
   10 CONTINUE
      CBA    = CB/CA
      IF(IOUT2.GE.5) WRITE(NB6,620)
C *** LOOP OVER OPTIONS.
      DO 600 MOP=1,4
      IF(IPERT(MOP).LE.0) GO TO 600
      EBW    = ZERO
      IBW    = 0
   20 IBW    = IBW+1
      DBW    = DELT-EBW
      DEBUG  = IOUT2.GE.5 .AND. IBW.EQ.1
      IF(DEBUG) WRITE(NB6,650)
C     READ MATRIX ELEMENTS.
      IF(NNTOT.GT.LM9) THEN
         REWIND NB2
         READ(NB2) W
      ENDIF
      NN     = 0
C *** SINGLE EXCITATIONS, I-LMO AND KMO-I.
      ESUM1  = ZERO
      DO 100 I=1,IB
      II     = IMOCI(I)
      IF(II.EQ.LL) GO TO 100
      IS     = NSYM(II)
      NSIK   = IPROD(IS,KKS)
      NSIL   = IPROD(IS,LLS)
      IF(NSIK.NE.1 .AND. I.GE.IA) GO TO 100
      IF(NSIL.NE.1 .AND. I.LT.IA) GO TO 100
      NN     = NN+1
      IF(NN.GT.LM9) THEN
         READ(NB2) W
         NN  = 1
      ENDIF
      BA     = W(NN)
      IF(I.LT.IA) THEN
         EA  = EIG(LL)-EIG(II)
         IF(MOP.EQ.2 .OR. MOP.EQ.4) EA=EA-G(LMO,I)+TWO*G(I,LMO)
      ELSE
         EA  = EIG(II)-EIG(KK)
         IF(MOP.EQ.2 .OR. MOP.EQ.4) EA=EA-G(I,KMO)+TWO*G(KMO,I)
      ENDIF
      EE     =-BA*BA/(EA+DBW)
      ESUM1  = ESUM1+EE
      IF(DEBUG) WRITE(NB6,900) II,I4,I4,I4,BA,EA,EE
  100 CONTINUE
      IF(DEBUG) WRITE(NB6,650)
C *** DOUBLE AND QUADRUPLE EXCITATIONS.
      ESUM2  = ZERO
      ESUM4  = ZERO
      DO 203 I=IA,IB
      II     = IMOCI(I)
      IS     = NSYM(II)
      DO 202 J=1,JB
      IJ     = IMOCI(J)
      JS     = NSYM(IJ)
      NSIJ   = IPROD(IS,JS)
      DO 201 K=IA,I
      IK     = IMOCI(K)
      KS     = NSYM(IK)
      DO 200 L=1,J
      IL     = IMOCI(L)
      LS     = NSYM(IL)
      NSKL   = IPROD(KS,LS)
      IF(NSIJ.NE.NSKL) GO TO 200
      IF(II.EQ.LL .AND. IJ.EQ.KK .AND. IK.EQ.LL .AND. IL.EQ.KK) GOTO 200
C     EXCITATIONS J,L-I,K.
      NN     = NN+1
      IF(NN.GT.LM9) THEN
         READ(NB2) W
         NN  = 1
      ENDIF
      BA     = W(NN)
      CALL ERGCON (EA,EB,EIG,IMOCI,LM3,J,L,I4,I4,I,K,I4,I4,MOP,
     1             G,LM10)
      EE     =-BA*BA/(EA+DBW)
      ESUM2  = ESUM2+EE
      IF(DEBUG) WRITE(NB6,900) II,IJ,IK,IL,BA,EA,EE
      IF(I.EQ.K .OR. J.EQ.L) GO TO 210
      NN     = NN+1
      IF(NN.GT.LM9) THEN
         READ(NB2) W
         NN  = 1
      ENDIF
      BB     = W(NN)
      EE     =-BB*BB/(EB+DBW)
      ESUM2  = ESUM2+EE
      IF(DEBUG) WRITE(NB6,900) II,IJ,IK,IL,BB,EB,EE
C     EXCITATIONS J,L,KMO,KMO-I,K,LMO,LMO.
  210 IF(II.EQ.LL .OR. IJ.EQ.KK) GO TO 200
      IF(IK.EQ.LL .OR. IL.EQ.KK) GO TO 200
      IF(IABSCI.EQ.3) GO TO 200
      CALL ERGCON (EA,EB,EIG,IMOCI,LM3,J,L,KMO,KMO,I,K,LMO,LMO,MOP,
     1             G,LM10)
      BA     = BA*CBA
      EE     =-BA*BA/(EA+DBW)
      ESUM4  = ESUM4+EE
      IF(DEBUG) WRITE(NB6,900) II,IJ,IK,IL,BA,EA,EE
      IF(I.EQ.K .OR. J.EQ.L) GO TO 200
      BB     = BB*CBA
      EE     =-BB*BB/(EB+DBW)
      ESUM4  = ESUM4+EE
      IF(DEBUG) WRITE(NB6,900) II,IJ,IK,IL,BB,EB,EE
  200 CONTINUE
  201 CONTINUE
  202 CONTINUE
  203 CONTINUE
      IF(DEBUG) WRITE(NB6,650)
C *** TRIPLE EXCITATIONS.
      ESUM3  = ZERO
      IF(IABSCI.EQ.3) GO TO 550
C     EXCITATIONS I,KMO,KMO-J,K,LMO  (TYPE 8B,9,10).
      DO 300 I=1,JB
      II     = IMOCI(I)
      IF(II.EQ.KK) GO TO 300
      IS     = NSYM(II)
      NSILL  = IPROD(IS,LLS)
      DO 310 J=IA,IB
      IJ     = IMOCI(J)
      IF(IJ.EQ.LL) GO TO 310
      JS     = NSYM(IJ)
      DO 320 K=IA,J
      IK     = IMOCI(K)
      IF(IK.EQ.LL) GO TO 320
      KS     = NSYM(IK)
      NSJK   = IPROD(JS,KS)
      IF(NSILL.NE.NSJK) GO TO 320
      NN     = NN+1
      IF(NN.GT.LM9) THEN
         READ(NB2) W
         NN  = 1
      ENDIF
      BA     = W(NN)
      CALL ERGCON (EA,EB,EIG,IMOCI,LM3,I,KMO,KMO,I4,J,K,LMO,I4,MOP,
     1             G,LM10)
      EE     =-BA*BA/(EA+DBW)
      ESUM3  = ESUM3+EE
      IF(DEBUG) WRITE(NB6,900) II,IJ,IK,I4,BA,EA,EE
      IF(J.EQ.K) GO TO 320
      NN     = NN+1
      IF(NN.GT.LM9) THEN
         READ(NB2) W
         NN  = 1
      ENDIF
      BB     = W(NN)
      EE     =-BB*BB/(EB+DBW)
      ESUM3  = ESUM3+EE
      IF(DEBUG) WRITE(NB6,900) II,IJ,IK,I4,BB,EB,EE
  320 CONTINUE
  310 CONTINUE
  300 CONTINUE
C     EXCITATIONS J,K,KMO-I,LMO,LMO  (TYPE 8A,12,13).
      DO 400 I=IA,IB
      II     = IMOCI(I)
      IF(II.EQ.LL) GO TO 400
      IS     = NSYM(II)
      NSIKK  = IPROD(IS,KKS)
      DO 410 J=1,JB
      IJ     = IMOCI(J)
      IF(IJ.EQ.KK) GO TO 410
      JS     = NSYM(IJ)
      DO 420 K=1,J
      IK     = IMOCI(K)
      IF(IK.EQ.KK) GO TO 420
      KS     = NSYM(IK)
      NSJK   = IPROD(JS,KS)
      IF(NSIKK.NE.NSJK) GO TO 420
      NN     = NN+1
      IF(NN.GT.LM9) THEN
         READ(NB2) W
         NN  = 1
      ENDIF
      BA     = W(NN)
      CALL ERGCON (EA,EB,EIG,IMOCI,LM3,J,K,KMO,I4,I,LMO,LMO,I4,MOP,
     1             G,LM10)
      EE     =-BA*BA/(EA+DBW)
      ESUM3  = ESUM3+EE
      IF(DEBUG) WRITE(NB6,900) II,IJ,IK,I4,BA,EA,EE
      IF(J.EQ.K) GO TO 420
      NN     = NN+1
      IF(NN.GT.LM9) THEN
         READ(NB2) W
         NN  = 1
      ENDIF
      BB     = W(NN)
      EE     =-BB*BB/(EB+DBW)
      ESUM3  = ESUM3+EE
      IF(DEBUG) WRITE(NB6,900) II,IJ,IK,I4,BB,EB,EE
  420 CONTINUE
  410 CONTINUE
  400 CONTINUE
C     EXCITATIONS I,KMO,KMO-J,LMO,LMO  (TYPE 8C).
      DO 500 I=1,JB
      II     = IMOCI(I)
      IF(II.EQ.KK) GO TO 500
      IS     = NSYM(II)
      DO 510 J=IA,IB
      IJ     = IMOCI(J)
      IF(IJ.EQ.LL) GO TO 510
      JS     = NSYM(IJ)
      NSIJ   = IPROD(IS,JS)
      IF(NSIJ.NE.1) GO TO 510
      NN     = NN+1
      IF(NN.GT.LM9) THEN
         READ(NB2) W
         NN  = 1
      ENDIF
      BA     = W(NN)
      CALL ERGCON (EA,EB,EIG,IMOCI,LM3,I,KMO,KMO,I4,J,LMO,LMO,I4,MOP,
     1             G,LM10)
      EE     =-BA*BA/(EA+DBW)
      ESUM3  = ESUM3+EE
      IF(DEBUG) WRITE(NB6,900) II,IJ,I4,I4,BA,EA,EE
  510 CONTINUE
  500 CONTINUE
C *** ENERGY CORRECTION.
  550 EESUM  = ESUM1+ESUM2+ESUM3+ESUM4
C *** BW TREATMENT.
      IF(MOP.GT.2) THEN
         IF(IOUT2.GE.0) THEN
            IF(IBW.EQ.1) WRITE(NB6,650)
            WRITE(NB6,670) IBW,EESUM,ESUM1,ESUM2,ESUM3,ESUM4
         ENDIF
         IF(ABS(EESUM-EBW).GE.EBWLIM) THEN
            EBW = EESUM
            GO TO 20
         ENDIF
      ENDIF
C *** PRINT THE RESULTS.
      ESUM   = EESUM-DELT
      EEV    = ESUM*EV
      EKCAL  = EEV*EVCAL
      DELTAU =-DELT
      DELTEV = DELTAU*EV
      DELTKC = DELTEV*EVCAL
      HDCAL  = DENER+EKCAL
      HDEV   = HDCAL/EVCAL
      HDSUM  = HDEV/EV
      IF(IPERT(MOP).EQ.1) E2=EEV
      IF(IOUT2.LE.-5) GO TO 600
      IOPCI  = IABSCI-2
      WRITE(NB6,920) DELTAU,DELTEV,DELTKC,KK,LL
      WRITE(NB6,910) ESUM,EEV,EKCAL,NAME(MOP),IOPCI
      WRITE(NB6,930) HDSUM,HDEV,HDCAL,NAME(MOP),IOPCI
  600 CONTINUE
      RETURN
  620 FORMAT (//1X,'DEBUG PRINT FOR SECOND-ORDER TERMS.'//
     1          1X,'  II  IJ  IK  IL',5X,'BA',10X,'EA',10X,'EE')
  650 FORMAT (  1X)
  670 FORMAT (  1X,'CYCLE',I3 ,5X,'ENERGIES (IN A.U.)',3X,5F12.7)
  900 FORMAT (  1X,4I4,6F12.7)
  910 FORMAT (/ 1X,'SECOND-ORDER ENERGY =',F12.7,' A.U.=',F12.7,' EV=',
     1             F12.7,' KCAL/MOL',10X,'OPTION=',A4,I1)
  920 FORMAT (/ 1X,'CI STABILIZATION    =',F12.7,' A.U.=',F12.7,' EV=',
     1             F12.7,' KCAL/MOL',10X,'MO-VO ',2I3)
  930 FORMAT (/ 1X,'HEAT OF FORMATION   =',F12.7,' A.U.=',F12.7,' EV=',
     1             F12.7,' KCAL/MOL',10X,'OPTION=',A4,I1)
      END
