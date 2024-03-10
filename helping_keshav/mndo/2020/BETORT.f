      SUBROUTINE BETORT (H,S,B,COR,HR,LM2,LM4,ICNTR)
C     *
C     BETORT COMPUTES ORTHOGONALIZATION CORRECTIONS TO THE RESONANCE
C     INTEGRALS (IN EV) AS DEFINED IN THE OM2 METHOD.
C     BETORT IS CALLED ONCE FOR A GIVEN GEOMETRY.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     H(LM4)    ORTHOGONALIZATION CORRECTIONS (O).
C     S(LM4)    OVERLAP INTEGRALS (I).
C     B(LM4)    RESONANCE INTEGRALS (I).
C     COR()     MATRIX OF LOCAL CORE-ELECTRON ATTRACTIONS (I).
C     HR(LM4)   RELATIVE MAGNITUDE OF ORTHOGONALIZATION CORRECTIONS
C               PER ATOM PAIR AT THE REFERENCE GEOMETRY (I).
C               USED ONLY FOR ICNTR.GT.0 (GRADIENT CALCULATION).
C     ICNTR     TYPE OF INTEGRAL CALCULATION.
C               = 0  FULL ENERGY EVALUATION.
C               = I  GRADIENT EVALUATION FOR ATOM I.
C               =-I  GRADIENT EVALUATION FOR ATOM I - NO TEST ON CUTG.
C     *
      USE LIMIT, ONLY: LM1, LMX, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (SMALLG=1.0D-06)
      PARAMETER (SMALLS=1.0D-12)
      PARAMETER (PT625=0.0625D0)
      LOGICAL CASE1,CASE2,CASE3
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./INDEX / INDX(LMX)
     ./INOPT2/ IN2(300)
     ./PAROMM/ BETPI(LMZ*10),FTS1(LMZ),FTS2(LMZ),ZSCOR(LMZ*4)
C    ./NBFILE/ NBF(20)
      DIMENSION H(LM4),S(LM4),B(LM4)
      DIMENSION COR(LM2,LM2)
      DIMENSION HR(LM4)
      DIMENSION SIK(16),BIK(16),SJK(16),BJK(16)
      DIMENSION TS1(16),TS2(16),HI(4),HJ(4),HK(4)
C *** CONTROL VARIABLES.
C     NPRINT = IN2(72)
      MODE   = ABS(ICNTR)
      CUTS   = SMALLS
      CUTG   = SMALLG
      ICUTS  = IN2(25)
      ICUTG  = IN2(26)
      IF(ICUTS.GT.0) CUTS = 10.0D0**(-ICUTS)
      IF(ICUTG.GT.0) CUTG = 10.0D0**(-ICUTG)
C *** LOOP OVER ATOMS (I).
C$OMP PARALLEL DO SCHEDULE(dynamic)
C$OMP.default(none)
C$OMP.SHARED(numat,nat,nfirst,nlast,icntr,icutg,indx,hr,cutg)
C$OMP.SHARED(zero,pt25,fts1,fts2)
C$OMP.SHARED(mode,icuts,cuts,s,b,cor,h)
C$OMP.
C$OMP.PRIVATE(I,NI,IA,iorbs,j,ij,nj,ja,jorbs,case1,case2,case3)
C$OMP.PRIVATE(lmax,l,ts1,ts2,ft1,ft2)
C$OMP.PRIVATE(k,ka,korbs,kfirst,klast,ikss,jkss)
C$OMP.PRIVATE(iks,sik,bik,ikx,iky,ikz)
C$OMP.PRIVATE(jks,sjk,bjk,jkx,jky,jkz)
C$OMP.PRIVATE(hi,hj,hk)
C$OMP.PRIVATE(is,ix,iy,iz)
      DO 110 I=2,NUMAT
      NI     = NAT(I)
      IA     = NFIRST(I)
      IORBS  = NLAST(I)-IA+1
C *** LOOP OVER ATOMS (J).
      DO 100 J=1,I-1
      IF(ICNTR.GT.0 .AND. ICUTG.GE.0) THEN
         IJ  = INDX(I)+J
         IF(HR(IJ).LT.CUTG) GO TO 100
      ENDIF
      NJ     = NAT(J)
      JA     = NFIRST(J)
      JORBS  = NLAST(J)-JA+1
      CASE1  = IORBS.EQ.1 .AND. JORBS.EQ.1
      CASE2  = IORBS.EQ.1 .AND. JORBS.GE.4
      CASE3  = IORBS.GE.4 .AND. JORBS.EQ.1
      IF(CASE1) THEN
         LMAX = 1
      ELSE IF(CASE2) THEN
         LMAX = 4
      ELSE
         LMAX = 16
      ENDIF
      DO 10 L=1,LMAX
      TS1(L) = ZERO
      TS2(L) = ZERO
   10 CONTINUE
C     SCALE FACTORS FOR PERTURBATION TERMS.
      FT1    = PT25 *(FTS1(NI)+FTS1(NJ))
      FT2    = PT625*(FTS2(NI)+FTS2(NJ))
C *** LOOP OVER ATOMS (K) TO INCLUDE THREE-CENTER TERMS.
      KFIRST = 1
      KLAST  = NUMAT
C     IDENTIFY THREE-CENTER TERMS THAT MAY CONTRIBUTE TO GRADIENT.
      IF(MODE.GT.0) THEN
         IF(MODE.NE.I .AND. MODE.NE.J) THEN
            KFIRST = MODE
            KLAST  = MODE
         ENDIF
      ENDIF
      DO 90 K=KFIRST,KLAST
      IF(K.EQ.I .OR. K.EQ.J) GO TO 90
C     INITIALIZATION.
      KA     = NFIRST(K)
      KORBS  = NLAST(K)-KA+1
C     CHECK FOR CUTOFF IN CASE OF SMALL OVERLAP.
      IF(ICUTS.GE.0) THEN
         IF(I.GT.K) THEN
            IKSS = INDX(IA)+KA
         ELSE
            IKSS = INDX(KA)+IA
         ENDIF
C        IF(ABS(S(IKSS)).LT.CUTS) GO TO 90
         IF(J.GT.K) THEN
            JKSS = INDX(JA)+KA
         ELSE
            JKSS = INDX(KA)+JA
         ENDIF
         IF((S(IKSS)*S(JKSS)).LT.CUTS) GO TO 90
C        IF(ABS(S(IKSS)*S(JKSS)).LT.CUTS) GO TO 90
      ENDIF
C *** PREPARATION FOR FIRST PERTURBATION SUM.
      IF(I.GT.K) THEN
         IKS      = INDX(IA)+KA
         IF(IORBS.EQ.1) THEN
            DO 30 L=1,KORBS
            SIK(L)   = S(IKS-1+L)
            BIK(L)   = B(IKS-1+L)
   30       CONTINUE
         ELSE
            IKX      = INDX(IA+1)+KA
            IKY      = INDX(IA+2)+KA
            IKZ      = INDX(IA+3)+KA
            DO 35 L=1,KORBS
            SIK(L)   = S(IKS-1+L)
            SIK(L+4) = S(IKX-1+L)
            SIK(L+8) = S(IKY-1+L)
            SIK(L+12)= S(IKZ-1+L)
            BIK(L)   = B(IKS-1+L)
            BIK(L+4) = B(IKX-1+L)
            BIK(L+8) = B(IKY-1+L)
            BIK(L+12)= B(IKZ-1+L)
   35       CONTINUE
         ENDIF
      ELSE
         IKS      = INDX(KA)+IA
         SIK(1)   = S(IKS)
         BIK(1)   = B(IKS)
         IF(IORBS.GE.4 .AND. KORBS.EQ.1) THEN
            SIK(5)   = S(IKS+1)
            SIK(9)   = S(IKS+2)
            SIK(13)  = S(IKS+3)
            BIK(5)   = B(IKS+1)
            BIK(9)   = B(IKS+2)
            BIK(13)  = B(IKS+3)
         ELSE IF(IORBS.EQ.1 .AND. KORBS.GE.4) THEN
            IKX      = INDX(KA+1)+IA
            IKY      = INDX(KA+2)+IA
            IKZ      = INDX(KA+3)+IA
            SIK(2)   = S(IKX)
            SIK(3)   = S(IKY)
            SIK(4)   = S(IKZ)
            BIK(2)   = B(IKX)
            BIK(3)   = B(IKY)
            BIK(4)   = B(IKZ)
         ELSE IF(IORBS.GE.4 .AND. KORBS.GE.4) THEN
            IKX      = INDX(KA+1)+IA
            IKY      = INDX(KA+2)+IA
            IKZ      = INDX(KA+3)+IA
            SIK(5)   = S(IKS+1)
            SIK(9)   = S(IKS+2)
            SIK(13)  = S(IKS+3)
            SIK(2)   = S(IKX)
            SIK(6)   = S(IKX+1)
            SIK(10)  = S(IKX+2)
            SIK(14)  = S(IKX+3)
            SIK(3)   = S(IKY)
            SIK(7)   = S(IKY+1)
            SIK(11)  = S(IKY+2)
            SIK(15)  = S(IKY+3)
            SIK(4)   = S(IKZ)
            SIK(8)   = S(IKZ+1)
            SIK(12)  = S(IKZ+2)
            SIK(16)  = S(IKZ+3)
            BIK(5)   = B(IKS+1)
            BIK(9)   = B(IKS+2)
            BIK(13)  = B(IKS+3)
            BIK(2)   = B(IKX)
            BIK(6)   = B(IKX+1)
            BIK(10)  = B(IKX+2)
            BIK(14)  = B(IKX+3)
            BIK(3)   = B(IKY)
            BIK(7)   = B(IKY+1)
            BIK(11)  = B(IKY+2)
            BIK(15)  = B(IKY+3)
            BIK(4)   = B(IKZ)
            BIK(8)   = B(IKZ+1)
            BIK(12)  = B(IKZ+2)
            BIK(16)  = B(IKZ+3)
         ENDIF
      ENDIF
      IF(J.GT.K) THEN
         JKS      = INDX(JA)+KA
         IF(JORBS.EQ.1) THEN
            DO 40 L=1,KORBS
            SJK(L)   = S(JKS-1+L)
            BJK(L)   = B(JKS-1+L)
   40       CONTINUE
         ELSE
            JKX      = INDX(JA+1)+KA
            JKY      = INDX(JA+2)+KA
            JKZ      = INDX(JA+3)+KA
            DO 45 L=1,KORBS
            SJK(L)   = S(JKS-1+L)
            SJK(L+4) = S(JKX-1+L)
            SJK(L+8) = S(JKY-1+L)
            SJK(L+12)= S(JKZ-1+L)
            BJK(L)   = B(JKS-1+L)
            BJK(L+4) = B(JKX-1+L)
            BJK(L+8) = B(JKY-1+L)
            BJK(L+12)= B(JKZ-1+L)
   45       CONTINUE
         ENDIF
      ELSE
         JKS      = INDX(KA)+JA
         SJK(1)   = S(JKS)
         BJK(1)   = B(JKS)
         IF(JORBS.GE.4 .AND. KORBS.EQ.1) THEN
            SJK(5)   = S(JKS+1)
            SJK(9)   = S(JKS+2)
            SJK(13)  = S(JKS+3)
            BJK(5)   = B(JKS+1)
            BJK(9)   = B(JKS+2)
            BJK(13)  = B(JKS+3)
         ELSE IF(JORBS.EQ.1 .AND. KORBS.GE.4) THEN
            JKX      = INDX(KA+1)+JA
            JKY      = INDX(KA+2)+JA
            JKZ      = INDX(KA+3)+JA
            SJK(2)   = S(JKX)
            SJK(3)   = S(JKY)
            SJK(4)   = S(JKZ)
            BJK(2)   = B(JKX)
            BJK(3)   = B(JKY)
            BJK(4)   = B(JKZ)
         ELSE IF(JORBS.GE.4 .AND. KORBS.GE.4) THEN
            JKX      = INDX(KA+1)+JA
            JKY      = INDX(KA+2)+JA
            JKZ      = INDX(KA+3)+JA
            SJK(5)   = S(JKS+1)
            SJK(9)   = S(JKS+2)
            SJK(13)  = S(JKS+3)
            SJK(2)   = S(JKX)
            SJK(6)   = S(JKX+1)
            SJK(10)  = S(JKX+2)
            SJK(14)  = S(JKX+3)
            SJK(3)   = S(JKY)
            SJK(7)   = S(JKY+1)
            SJK(11)  = S(JKY+2)
            SJK(15)  = S(JKY+3)
            SJK(4)   = S(JKZ)
            SJK(8)   = S(JKZ+1)
            SJK(12)  = S(JKZ+2)
            SJK(16)  = S(JKZ+3)
            BJK(5)   = B(JKS+1)
            BJK(9)   = B(JKS+2)
            BJK(13)  = B(JKS+3)
            BJK(2)   = B(JKX)
            BJK(6)   = B(JKX+1)
            BJK(10)  = B(JKX+2)
            BJK(14)  = B(JKX+3)
            BJK(3)   = B(JKY)
            BJK(7)   = B(JKY+1)
            BJK(11)  = B(JKY+2)
            BJK(15)  = B(JKY+3)
            BJK(4)   = B(JKZ)
            BJK(8)   = B(JKZ+1)
            BJK(12)  = B(JKZ+2)
            BJK(16)  = B(JKZ+3)
         ENDIF
      ENDIF
C *** PREPARATION FOR SECOND PERTURBATION SUM.
      HI(1)    = COR(IA,K)
      IF(IORBS.GE.4) THEN
         HI(2) = COR(IA+1,K)
         HI(3) = HI(2)
         HI(4) = HI(2)
      ENDIF
      HJ(1)    = COR(JA,K)
      IF(JORBS.GE.4) THEN
         HJ(2) = COR(JA+1,K)
         HJ(3) = HJ(2)
         HJ(4) = HJ(2)
      ENDIF
      HK(1)    = COR(KA,I)+COR(KA,J)
      IF(KORBS.GE.4) THEN
         HK(2) = COR(KA+1,I)+COR(KA+1,J)
         HK(3) = HK(2)
         HK(4) = HK(2)
      ENDIF
C
C *** FIRST AND SECOND PERTURBATION SUM.
C     S(I)-S(J)
      DO 50 L=1,KORBS
      TS1(1) = TS1(1)+(SIK(L)*BJK(L)+BIK(L)*SJK(L))
      TS2(1) = TS2(1)+(SIK(L)*SJK(L)*(HI(1)+HJ(1)-HK(L)))
   50 CONTINUE
      IF(CASE1) GO TO 90
C     S(I)-P(J)
      IF(CASE2) THEN
         DO 60 L=1,KORBS
         TS1(2)=TS1(2)+(SIK(L)*BJK(L+4)+BIK(L)*SJK(L+4))
         TS1(3)=TS1(3)+(SIK(L)*BJK(L+8)+BIK(L)*SJK(L+8))
         TS1(4)=TS1(4)+(SIK(L)*BJK(L+12)+BIK(L)*SJK(L+12))
         TS2(2)=TS2(2)+(SIK(L)*SJK(L+4)*(HI(1)+HJ(2)-HK(L)))
         TS2(3)=TS2(3)+(SIK(L)*SJK(L+8)*(HI(1)+HJ(3)-HK(L)))
         TS2(4)=TS2(4)+(SIK(L)*SJK(L+12)*(HI(1)+HJ(4)-HK(L)))
   60    CONTINUE
      ELSE IF(CASE3) THEN
C     P(I)-S(J)
         DO 70 L=1,KORBS
         TS1(5) = TS1(5) +(SIK(L+4)*BJK(L)+BIK(L+4)*SJK(L))
         TS1(9) = TS1(9) +(SIK(L+8)*BJK(L)+BIK(L+8)*SJK(L))
         TS1(13)= TS1(13)+(SIK(L+12)*BJK(L)+BIK(L+12)*SJK(L))
         TS2(5) = TS2(5) +(SIK(L+4)*SJK(L)*(HI(2)+HJ(1)-HK(L)))
         TS2(9) = TS2(9) +(SIK(L+8)*SJK(L)*(HI(3)+HJ(1)-HK(L)))
         TS2(13)= TS2(13)+(SIK(L+12)*SJK(L)*(HI(4)+HJ(1)-HK(L)))
   70    CONTINUE
      ELSE
C     S(I)-P(J), P(I)-S(J), P(I)-P(J).
         DO 80 L=1,KORBS
         TS1(2) = TS1(2) +(SIK(L)*BJK(L+4)+BIK(L)*SJK(L+4))
         TS1(3) = TS1(3) +(SIK(L)*BJK(L+8)+BIK(L)*SJK(L+8))
         TS1(4) = TS1(4) +(SIK(L)*BJK(L+12)+BIK(L)*SJK(L+12))
         TS1(5) = TS1(5) +(SIK(L+4)*BJK(L)+BIK(L+4)*SJK(L))
         TS1(6) = TS1(6) +(SIK(L+4)*BJK(L+4)+BIK(L+4)*SJK(L+4))
         TS1(7) = TS1(7) +(SIK(L+4)*BJK(L+8)+BIK(L+4)*SJK(L+8))
         TS1(8) = TS1(8) +(SIK(L+4)*BJK(L+12)+BIK(L+4)*SJK(L+12))
         TS1(9) = TS1(9) +(SIK(L+8)*BJK(L)+BIK(L+8)*SJK(L))
         TS1(10)= TS1(10)+(SIK(L+8)*BJK(L+4)+BIK(L+8)*SJK(L+4))
         TS1(11)= TS1(11)+(SIK(L+8)*BJK(L+8)+BIK(L+8)*SJK(L+8))
         TS1(12)= TS1(12)+(SIK(L+8)*BJK(L+12)+BIK(L+8)*SJK(L+12))
         TS1(13)= TS1(13)+(SIK(L+12)*BJK(L)+BIK(L+12)*SJK(L))
         TS1(14)= TS1(14)+(SIK(L+12)*BJK(L+4)+BIK(L+12)*SJK(L+4))
         TS1(15)= TS1(15)+(SIK(L+12)*BJK(L+8)+BIK(L+12)*SJK(L+8))
         TS1(16)= TS1(16)+(SIK(L+12)*BJK(L+12)+BIK(L+12)*SJK(L+12))
         TS2(2) = TS2(2) +(SIK(L)*SJK(L+4)*(HI(1)+HJ(2)-HK(L)))
         TS2(3) = TS2(3) +(SIK(L)*SJK(L+8)*(HI(1)+HJ(3)-HK(L)))
         TS2(4) = TS2(4) +(SIK(L)*SJK(L+12)*(HI(1)+HJ(4)-HK(L)))
         TS2(5) = TS2(5) +(SIK(L+4)*SJK(L)*(HI(2)+HJ(1)-HK(L)))
         TS2(6) = TS2(6) +(SIK(L+4)*SJK(L+4)*(HI(2)+HJ(2)-HK(L)))
         TS2(7) = TS2(7) +(SIK(L+4)*SJK(L+8)*(HI(2)+HJ(3)-HK(L)))
         TS2(8) = TS2(8) +(SIK(L+4)*SJK(L+12)*(HI(2)+HJ(4)-HK(L)))
         TS2(9) = TS2(9) +(SIK(L+8)*SJK(L)*(HI(3)+HJ(1)-HK(L)))
         TS2(10)= TS2(10)+(SIK(L+8)*SJK(L+4)*(HI(3)+HJ(2)-HK(L)))
         TS2(11)= TS2(11)+(SIK(L+8)*SJK(L+8)*(HI(3)+HJ(3)-HK(L)))
         TS2(12)= TS2(12)+(SIK(L+8)*SJK(L+12)*(HI(3)+HJ(4)-HK(L)))
         TS2(13)= TS2(13)+(SIK(L+12)*SJK(L)*(HI(4)+HJ(1)-HK(L)))
         TS2(14)= TS2(14)+(SIK(L+12)*SJK(L+4)*(HI(4)+HJ(2)-HK(L)))
         TS2(15)= TS2(15)+(SIK(L+12)*SJK(L+8)*(HI(4)+HJ(3)-HK(L)))
         TS2(16)= TS2(16)+(SIK(L+12)*SJK(L+12)*(HI(4)+HJ(4)-HK(L)))
   80    CONTINUE
      ENDIF
   90 CONTINUE
C
C *** DEBUG PRINT.
C     IF(NPRINT.GT.5) THEN
C        NB6 = NBF(6)
C        WRITE(NB6,614) FT1,FT2
C 614    FORMAT (1X,'BETORT W 614',2F15.10)
C        WRITE(NB6,616) (L,(FT1*TS1(L)),(FT2*TS2(L)),L=1,LMAX)
C 616    FORMAT(1X,'BETORT W 616',I5,2F15.10)
C     ENDIF
C
C *** ROTATED MATRIX ELEMENTS.
      IS      = INDX(IA)+JA
C     S(I)-S(J)
      IF(CASE1) THEN
         H(IS)   = -FT1*TS1(1)+FT2*TS2(1)
C     S(I)-S(J), S(I)-P(J)
      ELSE IF(CASE2) THEN
         H(IS)   = -FT1*TS1(1)+FT2*TS2(1)
         H(IS+1) = -FT1*TS1(2)+FT2*TS2(2)
         H(IS+2) = -FT1*TS1(3)+FT2*TS2(3)
         H(IS+3) = -FT1*TS1(4)+FT2*TS2(4)
      ELSE IF(CASE3) THEN
C     S(I)-S(J), P(I)-S(J)
         IX      = INDX(IA+1)+JA
         IY      = INDX(IA+2)+JA
         IZ      = INDX(IA+3)+JA
         H(IS)   = -FT1*TS1(1) +FT2*TS2(1)
         H(IX)   = -FT1*TS1(5) +FT2*TS2(5)
         H(IY)   = -FT1*TS1(9) +FT2*TS2(9)
         H(IZ)   = -FT1*TS1(13)+FT2*TS2(13)
      ELSE
C     S(I)-S(J), S(I)-P(J), P(I)-S(J), P(I)-P(J).
         IX      = INDX(IA+1)+JA
         IY      = INDX(IA+2)+JA
         IZ      = INDX(IA+3)+JA
         H(IS)   = -FT1*TS1(1) +FT2*TS2(1)
         H(IS+1) = -FT1*TS1(2) +FT2*TS2(2)
         H(IS+2) = -FT1*TS1(3) +FT2*TS2(3)
         H(IS+3) = -FT1*TS1(4) +FT2*TS2(4)
         H(IX)   = -FT1*TS1(5) +FT2*TS2(5)
         H(IX+1) = -FT1*TS1(6) +FT2*TS2(6)
         H(IX+2) = -FT1*TS1(7) +FT2*TS2(7)
         H(IX+3) = -FT1*TS1(8) +FT2*TS2(8)
         H(IY)   = -FT1*TS1(9) +FT2*TS2(9)
         H(IY+1) = -FT1*TS1(10)+FT2*TS2(10)
         H(IY+2) = -FT1*TS1(11)+FT2*TS2(11)
         H(IY+3) = -FT1*TS1(12)+FT2*TS2(12)
         H(IZ)   = -FT1*TS1(13)+FT2*TS2(13)
         H(IZ+1) = -FT1*TS1(14)+FT2*TS2(14)
         H(IZ+2) = -FT1*TS1(15)+FT2*TS2(15)
         H(IZ+3) = -FT1*TS1(16)+FT2*TS2(16)
      ENDIF
  100 CONTINUE
  110 CONTINUE
C$OMP END PARALLEL DO
      RETURN
      END
