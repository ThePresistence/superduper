      SUBROUTINE WORDER (CC,W,LM7,LM9,NB2,NB3,IOUT2)
C     *
C     ORDERING OF THE TWO-ELECTRON INTEGRALS (CASE IMODE.GT.0).
C     *
C     THE ONE-CENTER AND TWO-CENTER TWO-ELECTRON INTEGRALS ARE
C     TRANSFERRED TO THE ARRAY CC(LM7) AND STORED AS A SQUARE
C     MATRIX (LM6,LM6). IF THERE IS NOT ENOUGH MEMORY, THE
C     INTEGRALS ARE SAVED ON DISK (LREC RECORDS OF LENGTH LM7).
C     *
      USE LIMIT, ONLY: LM1, LMX, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DPARM4/ REPD(52,LMZ)
     ./DPARM5/ INTIJ(243),INTKL(243),INTREP(243),INTRF1(243),INTRF2(243)
     ./FLAG4 / INTSUM
     ./INDEX / INDX(LMX)
     ./INDEXW/ NW(LM1)
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./REP   / GSS(LMZ),GPP(LMZ),GSP(LMZ),GP2(LMZ),HSP(LMZ),HPP(LMZ)
      DIMENSION CC(LM7),W(LM9)
C     *
C     CONTROL VARIABLES.
C     *
C     NOTATION.
C     LM6    - NUMBER OF ONE-CENTER PAIRS PER MOLECULE.
C     LM7    - LENGTH OF BUFFER FOR TWO-ELECTRON INTEGRALS W(LM6,LM6).
C     KTOT   - NUMBER OF INTEGRALS PER BUFFER (MULTIPLE OF LM6).
C     LMAX   - NUMBER OF INTEGRAL COLUMNS (LENGTH LM6) PER BUFFER.
      ITLAST = NLAST(NUMAT)-NFIRST(NUMAT)+1
      LM6    = NW(NUMAT)+INDX(ITLAST)+ITLAST-1
      LMAX   = LM7/LM6
      KTOT   = LM6*LM6
      IF(LM6.GT.LMAX) KTOT=LM6*LMAX
C     INPUT OPTION.
      IOP    = IN2(2)
      IF(IOP.EQ.2) GO TO 500
C     *
C     SORTING OF INTEGRALS FOR MNDO AND MINDO.
C     *
C     LREC   - NUMBER OF RECORDS.
      LREC   = 1+(LM6-1)/LMAX
      IF(LREC.GT.1) REWIND NB3
C     LOOP OVER RECORDS.
      DO 400 L=1,LREC
      LA     = LMAX*(L-1)+1
      LB     = LMAX*L
      KP     = LM6*(LA-1)
C     INITIALIZE INTEGRALS.
      DO 10 KK=1,KTOT
      CC(KK) = ZERO
   10 CONTINUE
C     *
C     LOOP OVER ONE-CENTER INTEGRALS.
C     *
      DO 130 II=1,NUMAT
      IT     = NLAST(II)-NFIRST(II)+1
      IA     = NW(II)
      IF(IA.GT.LB) GO TO 140
      IB     = IA+INDX(IT)+IT-1
      IF(IB.LT.LA) GO TO 130
      NI     = NAT(II)
C     HYDROGEN.
      IF(IT.EQ.1) THEN
         KK  = LM6*(IA-1)+IA-KP
         CC(KK) = GSS(NI)
         GO TO 130
      ENDIF
C     ONE-CENTER INTEGRALS FOR SP-BASIS.
      DO 110 I=IA,IB
      IF(I.LT.LA .OR. I.GT.LB) GO TO 110
      KS     = LM6*(I-1)+IA-1-KP
      KGO    = I-IA+1
      KK     = KS+KGO
      IF(KGO.EQ.1) THEN
         CC(KK   )=GSS(NI)
         CC(KS+3 )=GSP(NI)
         CC(KS+6 )=GSP(NI)
         CC(KS+10)=GSP(NI)
      ELSE IF(KGO.EQ.2 .OR. KGO.EQ.4 .OR. KGO.EQ.7) THEN
         CC(KK   )=HSP(NI)
      ELSE IF(KGO.EQ.3 .OR. KGO.EQ.6 .OR. KGO.EQ.10) THEN
         CC(KS+1 )=GSP(NI)
         CC(KS+3 )=GP2(NI)
         CC(KS+6 )=GP2(NI)
         CC(KS+10)=GP2(NI)
         CC(KK   )=GPP(NI)
      ELSE
         CC(KK   )=HPP(NI)
      ENDIF
  110 CONTINUE
C     ONE-CENTER INTEGRALS FOR SPD-BASIS.
      IF(IT.LE.4) GO TO 130
      DO 120 INTG=1,243
      IJ     = INTIJ(INTG)+IA-1
      IF(IJ.LT.LA .OR. IJ.GT.LB) GO TO 120
      KS     = LM6*(IJ-1)+IA-1-KP
      KK     = KS+INTKL(INTG)
      CC(KK) = REPD(INTREP(INTG),NI)
  120 CONTINUE
  130 CONTINUE
  140 CONTINUE
C     *
C     LOOP OVER TWO-CENTER INTEGRALS, MNDO.
C     *
      IF(NUMAT.EQ.1) GO TO 400
      IF(IOP.GT.0) GO TO 300
      IF(INTSUM.GT.LM9) THEN
         REWIND NB2
         READ(NB2) W
      ENDIF
      NA     = 0
      DO 280 II=2,NUMAT
      IT     = NLAST(II)-NFIRST(II)+1
      IW     = INDX(IT)+IT
      IA     = NW(II)
      IB     = IA+IW-1
      DO 270 JJ=1,II-1
      JT     = NLAST(JJ)-NFIRST(JJ)+1
      JW     = INDX(JT)+JT
      JA     = NW(JJ)
      JB     = JA+JW-1
      NO     = IW*JW
      IF(INTSUM.GT.LM9) THEN
         KW  = NA+NO
         IF(KW.GT.LM9) THEN
            READ(NB2) W
            NA  = 0
         ENDIF
      ENDIF
C     CASE II.GT.JJ.
      IF(IA.GT.LB .OR. IB.LT.LA) GO TO 230
      DO 220 I=IA,IB
      IF(I.LT.LA .OR. I.GT.LB) GO TO 220
      KS     = LM6*(I-1)+JA-1-KP
      NS     = NA+JW*(I-IA)
      DO 210 J=1,JW
      KK     = KS+J
      NN     = NS+J
      CC(KK) = W(NN)
  210 CONTINUE
  220 CONTINUE
C     CASE II.LT.JJ.
  230 IF(JA.GT.LB .OR. JB.LT.LA) GO TO 260
      DO 250 J=JA,JB
      IF(J.LT.LA .OR. J.GT.LB) GO TO 250
      KS     = LM6*(J-1)+IA-1-KP
      NS     = NA+J-JA+1-JW
      DO 240 I=1,IW
      KK     = KS+I
      NN     = NS+JW*I
      CC(KK) = W(NN)
  240 CONTINUE
  250 CONTINUE
  260 NA     = NA+NO
  270 CONTINUE
  280 CONTINUE
C     WRITE INTEGRAL RECORD, IF NECESSARY.
      IF(LREC.GT.1) WRITE(NB3) CC
      GO TO 400
C     *
C     LOOP OVER TWO-CENTER INTEGRALS, MINDO.
C     *
  300 DO 380 II=2,NUMAT
      IT     = NLAST(II)-NFIRST(II)+1
      IW     = INDX(IT)+IT
      IA     = NW(II)
      IB     = IA+IW-1
      DO 370 JJ=1,II-1
      JT     = NLAST(JJ)-NFIRST(JJ)+1
      JW     = INDX(JT)+JT
      JA     = NW(JJ)
      JB     = JA+JW-1
      IJ     = INDX(II)+JJ
      WIJ    = W(IJ)
C     CASE II.GT.JJ.
      IF(IA.GT.LB .OR. IB.LT.LA) GO TO 330
      DO 320 I=1,IT
      IL     = IA-1+INDX(I+1)
      IF(IL.LT.LA .OR. IL.GT.LB) GO TO 320
      KS     = LM6*(IL-1)+JA-1-KP
      DO 310 J=1,JT
      KK     = KS+INDX(J+1)
      CC(KK) = WIJ
  310 CONTINUE
  320 CONTINUE
C     CASE II.LT.JJ.
  330 IF(JA.GT.LB .OR. JB.LT.LA) GO TO 370
      DO 350 J=1,JT
      JL     = JA-1+INDX(J+1)
      IF(JL.LT.LA .OR. JL.GT.LB) GO TO 350
      KS     = LM6*(JL-1)+IA-1-KP
      DO 340 I=1,IT
      KK     = KS+INDX(I+1)
      CC(KK) = WIJ
  340 CONTINUE
  350 CONTINUE
  370 CONTINUE
  380 CONTINUE
C     WRITE INTEGRAL RECORD, IF NECESSARY.
      IF(LREC.GT.1) WRITE(NB3) CC
  400 CONTINUE
      GO TO 600
C     *
C     SORTING OF INTEGRALS FOR CNDO.
C     INTEGRALS ARE REQUIRED TO REMAIN IN MEMORY.
C     *
  500 DO 520 II=1,NUMAT
      NI     = NAT(II)
      KK     = LM6*(II-1)+II
      CC(KK) = GSS(NI)
      IF(II.EQ.1) GO TO 520
      DO 510 JJ=1,II-1
      WIJ    = W(INDX(II)+JJ)
      KIJ    = LM6*(II-1)+JJ
      KJI    = LM6*(JJ-1)+II
      CC(KIJ)= WIJ
      CC(KJI)= WIJ
  510 CONTINUE
  520 CONTINUE
C     *
C     DEBUG PRINT.
C     *
  600 IF(IOUT2.GT.0) THEN
         NB6 = NBF(6)
         WRITE(NB6,700)
         WRITE(NB6,710) LM6
         IF(LREC.GT.1) WRITE(NB6,715) LREC,LMAX,KTOT
         IF(IOUT2.LT.6) RETURN
         IMAX   = LM6
         IF(LREC.GT.1) THEN
            IMAX = LMAX
            REWIND NB3
            READ(NB3) CC
         ENDIF
         WRITE(NB6,720)
         KK     = 1
         DO 610 I=1,IMAX
         KA     = KK
         KB     = KA+LM6-1
         WRITE(NB6,730) (CC(K),K=KA,KB)
         KK     = KA+LM6
  610    CONTINUE
      ENDIF
      RETURN
  700 FORMAT(///1X,'AO INTEGRALS IN NEW ORDER.'/)
  710 FORMAT(//1X,'THE AO INTEGRALS ARE STORED IN A MATRIX WITH',I4,
     1            ' ROWS AND COLUMNS.')
  715 FORMAT(  1X,'THERE ARE',I4,' RECORDS EACH OF WHICH CONTAINS',I4,
     1            ' COLUMNS AND',I6,' INTEGRALS.')
  720 FORMAT(//1X,'INTEGRALS IN THE FIRST RECORD.'/)
  730 FORMAT(  1X,10F10.5)
      END