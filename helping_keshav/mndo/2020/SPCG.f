      FUNCTION SPCG(C1,C2,C3,C4,W,LM2,LM9)
C     *
C     SPCG CALCULATES THE VALENCE-SHELL REPULSION INTEGRAL BETWEEN
C     ELECTRON 1 IN MOLECULAR ORBITALS C1,C2 AND
C     ELECTRON 2 IN MOLECULAR ORBITALS C3,C4.
C     *
      USE LIMIT, ONLY: LM1, LMI, LMX, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DPARM4/ REPD(52,LMZ)
     ./DPARM5/ INTIJ(243),INTKL(243),INTREP(243),INTRF1(243),INTRF2(243)
     ./FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
     ./FLAG4 / INTSUM
     ./INDEX / INDX(LMX)
     ./INDEXW/ NW(LM1)
     ./NBFILE/ NBF(20)
     ./INOPT2/ IN2(300)
     ./REP   / GSS(LMZ),GPP(LMZ),GSP(LMZ),GP2(LMZ),HSP(LMZ),HPP(LMZ)
      DIMENSION C1(LM2),C2(LM2),C3(LM2),C4(LM2),W(LM9)
C *** INITIALIZATION.
      IOP    = IN2(2)
      SPCG   = ZERO
C *** TWO-CENTRE TERMS.
      IF(NUMAT.EQ.1) GO TO 90
      KK     = 0
      IF(IOP.LE.0 .AND. INTSUM.GT.LM9) THEN
         NB2 = NBF(2)
         REWIND NB2
         READ(NB2) W
      ENDIF
      DO 80 II=2,NUMAT
      IA     = NFIRST(II)
      IB     = NLAST(II)
      IORBS  = IB-IA+1
      IW     = INDX(IORBS)+IORBS
      DO 70 JJ=1,II-1
      JA     = NFIRST(JJ)
      JB     = NLAST(JJ)
C     MINDO, CNDO, INDO.
      IF(IOP.GT.0) THEN
         IJ  = INDX(II)+JJ
CDIR$ NOVECTOR
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 20 I=IA,IB
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 10 J=JA,JB
         SPCG   = SPCG+W(IJ)*(C1(I)*C2(I)*C3(J)*C4(J)
     1                       +C1(J)*C2(J)*C3(I)*C4(I))
   10    CONTINUE
   20    CONTINUE
      ELSE
C     MNDO AND RELATED METHODS.
         IF(INTSUM.GT.LM9) THEN
            JORBS = JB-JA+1
            JW  = INDX(JORBS)+JORBS
            KW  = KK+IW*JW
            IF(KW.GT.LM9) THEN
               READ(NB2) W
               KK = 0
            ENDIF
         ENDIF
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 60 I=IA,IB
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 50 J=IA,I
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 40 K=JA,JB
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
         DO 30 L=JA,K
         PROD   = C1(I)*C2(J)*C3(K)*C4(L)
     1           +C1(K)*C2(L)*C3(I)*C4(J)
         IF(I.NE.J) PROD=PROD+C1(J)*C2(I)*C3(K)*C4(L)
     1                       +C1(K)*C2(L)*C3(J)*C4(I)
         IF(K.NE.L) PROD=PROD+C1(I)*C2(J)*C3(L)*C4(K)
     1                       +C1(L)*C2(K)*C3(I)*C4(J)
         IF(I.NE.J .AND. K.NE.L) PROD=PROD+C1(J)*C2(I)*C3(L)*C4(K)
     1                                    +C1(L)*C2(K)*C3(J)*C4(I)
         KK     = KK+1
         SPCG   = SPCG+PROD*W(KK)
   30    CONTINUE
   40    CONTINUE
   50    CONTINUE
   60    CONTINUE
      ENDIF
   70 CONTINUE
   80 CONTINUE
   90 CONTINUE
C *** ONE-CENTRE TERMS.
      DO 110 II=1,NUMAT
      IS     = NFIRST(II)
      IORBS  = NLAST(II)-IS+1
      NI     = NAT(II)
C     (SS/SS)
      SPCG   = SPCG+C1(IS)*C2(IS)*C3(IS)*C4(IS)*GSS(NI)
      IF(IORBS.LE.1) GO TO 110
      IX     = IS+1
      IY     = IS+2
      IZ     = IS+3
C     (PP/PP) FOR P=X,Y,Z
      SPCG   = SPCG+GPP(NI)*
     1             (C1(IX)*C2(IX)*C3(IX)*C4(IX)+
     2              C1(IY)*C2(IY)*C3(IY)*C4(IY)+
     3              C1(IZ)*C2(IZ)*C3(IZ)*C4(IZ) )
C     (SS/PP)+(PP/SS) FOR P=X,Y,Z
      SPCG   = SPCG+GSP(NI)*
     1             (C1(IS)*C2(IS)*C3(IX)*C4(IX)+
     2              C1(IS)*C2(IS)*C3(IY)*C4(IY)+
     3              C1(IS)*C2(IS)*C3(IZ)*C4(IZ)+
     4              C1(IX)*C2(IX)*C3(IS)*C4(IS)+
     5              C1(IY)*C2(IY)*C3(IS)*C4(IS)+
     6              C1(IZ)*C2(IZ)*C3(IS)*C4(IS) )
C     (PP/P*P*)+(P*P*/PP) FOR P.NE.P*=X,Y,Z
      SPCG   = SPCG+GP2(NI)*
     1             (C1(IX)*C2(IX)*C3(IY)*C4(IY)+
     2              C1(IX)*C2(IX)*C3(IZ)*C4(IZ)+
     3              C1(IY)*C2(IY)*C3(IZ)*C4(IZ)+
     4              C1(IY)*C2(IY)*C3(IX)*C4(IX)+
     5              C1(IZ)*C2(IZ)*C3(IX)*C4(IX)+
     6              C1(IZ)*C2(IZ)*C3(IY)*C4(IY) )
C     (SP/SP)+(SP/PS)+(PS/SP)+(PS/PS) FOR P=X,Y,Z
      SPCG   = SPCG+HSP(NI)*
     1             (C1(IS)*C2(IX)*C3(IX)*C4(IS)+
     2              C1(IS)*C2(IX)*C3(IS)*C4(IX)+
     3              C1(IX)*C2(IS)*C3(IS)*C4(IX)+
     4              C1(IX)*C2(IS)*C3(IX)*C4(IS)+
     5              C1(IS)*C2(IY)*C3(IY)*C4(IS)+
     6              C1(IS)*C2(IY)*C3(IS)*C4(IY)+
     7              C1(IY)*C2(IS)*C3(IS)*C4(IY)+
     8              C1(IY)*C2(IS)*C3(IY)*C4(IS)+
     9              C1(IS)*C2(IZ)*C3(IZ)*C4(IS)+
     A              C1(IS)*C2(IZ)*C3(IS)*C4(IZ)+
     B              C1(IZ)*C2(IS)*C3(IS)*C4(IZ)+
     C              C1(IZ)*C2(IS)*C3(IZ)*C4(IS) )
      SPCG   = SPCG+HPP(NI)*
     1             (C1(IX)*C2(IY)*C3(IX)*C4(IY)+
     2              C1(IX)*C2(IY)*C3(IY)*C4(IX)+
     3              C1(IY)*C2(IX)*C3(IY)*C4(IX)+
     4              C1(IY)*C2(IX)*C3(IX)*C4(IY)+
     5              C1(IX)*C2(IZ)*C3(IX)*C4(IZ)+
     6              C1(IX)*C2(IZ)*C3(IZ)*C4(IX)+
     7              C1(IZ)*C2(IX)*C3(IZ)*C4(IX)+
     8              C1(IZ)*C2(IX)*C3(IX)*C4(IZ)+
     9              C1(IY)*C2(IZ)*C3(IY)*C4(IZ)+
     A              C1(IY)*C2(IZ)*C3(IZ)*C4(IY)+
     B              C1(IZ)*C2(IY)*C3(IZ)*C4(IY)+
     C              C1(IZ)*C2(IY)*C3(IY)*C4(IZ) )
C     ONE-CENTER TWO-ELECTRON INTEGRALS INVOLVING D ORBITALS.
      IF(IORBS.LE.4) GO TO 110
      NWII   = NW(II)-1
      DO 100 INTG=1,243
      IJ     = INTIJ(INTG)+NWII
      KL     = INTKL(INTG)+NWII
      I      = IP1(IJ)
      J      = IP2(IJ)
      K      = IP1(KL)
      L      = IP2(KL)
      PROD   = C1(I)*C2(J)*C3(K)*C4(L)
      IF(I.NE.J) PROD=PROD+C1(J)*C2(I)*C3(K)*C4(L)
      IF(K.NE.L) PROD=PROD+C1(I)*C2(J)*C3(L)*C4(K)
      IF(I.NE.J .AND. K.NE.L) PROD=PROD+C1(J)*C2(I)*C3(L)*C4(K)
      SPCG   = SPCG+PROD*REPD(INTREP(INTG),NI)
  100 CONTINUE
  110 CONTINUE
      RETURN
      END
