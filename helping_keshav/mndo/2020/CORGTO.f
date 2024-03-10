      SUBROUTINE CORGTO (ISHELL,KTYP,R,VB)
C     *
C     BASIC ANALYTICAL CORE-ELECTRON ATTRACTION INTEGRALS FOR GTOS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     ISHELL    SHELL NUMBER OF ORBITALS AT ATOM A (I).
C     KTYP      SHELL TYPE OF ORBITALS AT ATOM A (O). 0 S, 1 SP.
C     R         DISTANCE A-B, IN AU (I).
C     VB        RESULTS, IN AU (O).
C     *
      USE LIMIT, ONLY: LMGP, LMGS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TWENTY=20.0D0)
      PARAMETER (ONEPT5=1.5D0)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./FMGTO1/ XLIM,XMAX
     ./GAUSS1/ KSTART(LMGS),KNG(LMGS),KTYPE(LMGS),NSHELL,NBASIS
     ./GAUSS3/ EXX(LMGP),C1(LMGP),C2(LMGP),C3(LMGP)
     ./MATRIX/ AA(400),BA(400),CA(400),
     .         AB(400),BB(400),CB(400),
     .         AC(400),BC(400),CC(400),AB2(2400)
      DIMENSION VB(4)
      DIMENSION FMT(3)
C *** INITIALIZATION.
      RR     = R*R
      ISTART = KSTART(ISHELL)
      INUMB  = KNG(ISHELL)
      KTYP   = KTYPE(ISHELL)
      IF(KTYP.GT.0) GO TO 30
C *** S-SHELL AT ATOM A.
      VB(1)  = ZERO
CDIR$ NOVECTOR
C$DIR SCALAR
      DO 20 II=1,INUMB
      IG     = ISTART+II-1
      A      = EXX(IG)
      CSA    = C1(IG)
C$DIR SCALAR
      DO 10 JJ=1,INUMB
      JG     = ISTART+JJ-1
      B      = EXX(JG)
      CSB    = C1(JG)
      G      = A+B
      X      = G*RR
C     CALL AUXG0
      IF(X.GT.XLIM) THEN
         F0  = PT5*SQRT(PI/X)
      ELSE IF(X.LT.XMAX) THEN
         QQ  = X*TWENTY
         TH  = QQ-AINT(QQ)
         I   = NINT(QQ-TH)
         TH  = QQ-I
         TH2 = TH *(TH-ONE)
         TH3 = TH2*(TH-TWO)
         TH4 = TH2*(TH+ONE)
         F0  = AA(I+1)+TH*BA(I+1)-TH3*CA(I+1)+TH4*CA(I+2)
      ELSE
         CALL FMTGEN (FMT,X,1,ICK)
         F0  = FMT(1)
      ENDIF
      V00    = F0/G
      VB(1)  = VB(1) + CSA*CSB*V00
   10 CONTINUE
   20 CONTINUE
      VB(1)  = VB(1)*TWO*PI
      RETURN
C *** SP-SHELL AT ATOM A.
   30 CONTINUE
      DO 40 K=1,4
      VB(K)  = ZERO
   40 CONTINUE
C$DIR SCALAR
      DO 60 II=1,INUMB
      IG     = ISTART+II-1
      A      = EXX(IG)
      CSA    = C1(IG)
      CPA    = C2(IG)
C$DIR SCALAR
      DO 50 JJ=1,INUMB
      JG     = ISTART+JJ-1
      B      = EXX(JG)
      CSB    = C1(JG)
      CPB    = C2(JG)
      G      = A+B
      X      = G*RR
C     CALL AUXG2
      IF(X.GT.XLIM) THEN
         F0  = PT5*SQRT(PI/X)
         F1  = PT5*F0/X
         F2  = ONEPT5*F1/X
      ELSE IF(X.LT.XMAX) THEN
         QQ  = X*TWENTY
         TH  = QQ-AINT(QQ)
         I   = NINT(QQ-TH)
         TH  = QQ-I
         TH2 = TH *(TH-ONE)
         TH3 = TH2*(TH-TWO)
         TH4 = TH2*(TH+ONE)
         F0  = AA(I+1)+TH*BA(I+1)-TH3*CA(I+1)+TH4*CA(I+2)
         F1  = AB(I+1)+TH*BB(I+1)-TH3*CB(I+1)+TH4*CB(I+2)
         F2  = AC(I+1)+TH*BC(I+1)-TH3*CC(I+1)+TH4*CC(I+2)
      ELSE
         CALL FMTGEN (FMT,X,3,ICK)
         F0  = FMT(1)
         F1  = FMT(2)
         F2  = FMT(3)
      ENDIF
      FF1    = TWO*PI/G
      FF2    = PI/(G*G)
      V00    = F0*FF1
      V03    = F1*FF1*R
      V11    = (F0-F1)*FF2
      V33    = V11 + F2*FF1*RR
      VB(1)  = VB(1) + CSA*CSB*V00
      VB(2)  = VB(2) + CSA*CPB*V03
      VB(3)  = VB(3) + CPA*CPB*V33
      VB(4)  = VB(4) + CPA*CPB*V11
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
