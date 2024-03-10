      SUBROUTINE DHCORE (H,W,LM4,LM9,I,J,NAT,NFIRST,NLAST,COORD,ENUCLR)
C     *
C     PARTIAL CORE HAMILTONIAN FOR GRADIENT CALCULATION.
C     GENERAL VERSION FOR ATOM PAIR I-J.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     H(LM4)    CORE HAMILTONIAN MATRIX FOR ATOM PAIR I-J (O).
C     W(LM9)    TWO-ELECTRON INTEGRALS FOR ATOM PAIR I-J (O).
C     I,J       ATOM NUMBERS FOR PAIR I-J IN MASTER LIST (I).
C     NAT(2)    ATOMIC NUMBER OF ATOMS I AND J (I).
C     NFIRST(2) NUMBER OF FIRST ORBITAL AT ATOMS I AND J (I).
C     NLAST(2)  NUMBER OF LAST  ORBITAL AT ATOMS I AND J (I).
C     COORD()   CARTESIAN COORDINATES FOR ATOMS I AND J (I).
C     ENUCLR    CORE-CORE REPULSION ENERGY FOR PAIR I-J (O).
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IZERO=0)
      LOGICAL NOTOMI
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./EXTRA1/ RI(22),CORE(10,2),REPT(10,2),WW(2025)
     ./EXTRA2/ SIJ(14),T(14),YY(675)
      COMMON
     ./INOPT2/ IN2(300)
     ./PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./PARMIN/ F0(18),VS(18),VP(18),BETA(55),ALPHA(55)
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ),BETAP(LMZ),
     .         ALP(LMZ)
      DIMENSION H(LM4),W(LM9)
      DIMENSION NAT(2),NFIRST(2),NLAST(2),COORD(3,2)
C$OMP THREADPRIVATE (/EXTRA1/, /EXTRA2/)
C *** CHECK FOR SPECIAL CASE OF H-H PAIR.
      IOP    = IN2(2)
      NOTOMI = IOP.NE.-5 .AND. IOP.NE.-6 .AND. IOP.NE.-8 .AND. IOP.NE.-9
     1                   .AND. IOP.NE.-22.AND. IOP.NE.-23
      IF(NAT(1).EQ.1 .AND. NAT(2).EQ.1 .AND. NLAST(2).EQ.2
     1               .AND. NOTOMI) THEN
         CALL HHPAIR (1,2,COORD,SIJSS,HIJ,WIJ,ENUCLR)
         H(1) =-WIJ
         H(2) = HIJ
         H(3) =-WIJ
         IF(IOP.LE.0) THEN
            W(1) = WIJ
         ELSE
            W(2) = WIJ
         ENDIF
         RETURN
      ENDIF
C *** INITIALIZE SOME VARIABLES.
      KR     = 0
      NI     = NAT(2)
      IA     = NFIRST(2)
      IB     = NLAST(2)
      NJ     = NAT(1)
      JA     = NFIRST(1)
      JB     = NLAST(1)
      IORBS  = IB-IA+1
      JORBS  = JB-JA+1
      ISHELL = I
      JSHELL = J
      DO 10 L=1,LM4
      H(L)   = ZERO
   10 CONTINUE
C *** ROTATION MATRIX AND INTERATOMIC DISTANCE R (AU).
      CALL ROTMAT (1,2,JORBS,IORBS,2,COORD,R,YY)
C     *
C *** SECTION FOR MNDO AND RELATED METHODS.
C     *
C *** TWO-ELECTRON INTEGRALS.
      IF(IOP.GT.0) GO TO 20
      IW     = (IORBS*(IORBS+1))/2
      JW     = (JORBS*(JORBS+1))/2
C     TWO-ELECTRON INTEGRALS IN LOCAL COORDINATES.
      IF(NOTOMI) THEN
C        COMPUTE AND STORE THE SEMIEMPIRICAL INTEGRALS.
         CALL REPP (NI,NJ,R,RI,CORE)
         IF(IORBS.GE.9 .OR. JORBS.GE.9) THEN
            CALL REPPD (NI,NJ,R,RI,CORE,W,IW,JW,1)
         ENDIF
      ELSE
C        COMPUTE AND STORE THE ANALYTICAL INTEGRALS.
         CALL REPGAU (R,RI,ISHELL,JSHELL,NI,NJ,IZERO)
         CALL CORDEF (RI,REPT,NI,NJ)
C        SCALE THE ANALYTICAL INTEGRALS.
         CALL GSCALE (R,RI,FKO,NI,NJ)
C        SAVE INTEGRALS TO BE USED FOR CORE-ELECTRON ATTRACTIONS.
         CALL CORDEF (RI,CORE,NI,NJ)
      ENDIF
C     TRANSFORM TWO-ELECTRON INTEGRALS TO MOLECULAR COORDINATES.
      IP     = (IA*(IA+1))/2
      JP     = (JA*(JA+1))/2
      IF(IORBS.LE.4 .AND. JORBS.LE.4) THEN
         CALL ROTATE (IW,JW,IP,JP,KR,RI,YY,W,W,IW+JW,LM9,1)
      ELSE
         CALL ROTD  (W,YY,IW,JW)
      ENDIF
C *** RESONANCE INTEGRALS.
      IF(NOTOMI) THEN
         CALL BETAIJ (NI,NJ,R,SIJ,T)
         CALL ROTBET (IA,JA,IORBS,JORBS,T,YY,H,LM4)
      ELSE IF(IOP.NE.-9) THEN
         CALL BETOM  (IZERO,ISHELL,JSHELL,NI,NJ,R,SIJ,T)
         CALL ROTBET (IA,JA,IORBS,JORBS,T,YY,H,LM4)
      ENDIF
C *** CORE-ELECTRON ATTRACTIONS.
      IF(NOTOMI) THEN
         CALL ROTCOH (IA,JA,IORBS,JORBS,IP,JP,CORE,YY,H,LM4)
      ELSE
         CALL COROM  (ISHELL,JSHELL,NI,NJ,R,FKO,CORE,REPT,SIJ,T)
         CALL ROTCOH (IA,JA,IORBS,JORBS,IP,JP,CORE,YY,H,LM4)
      ENDIF
C *** CORE-CORE REPULSIONS.
      IF(NOTOMI) THEN
         WIJ = RI(1)
         CALL CORREP (I,J,NI,NJ,R,WIJ,ENUCLR)
      ELSE
         ENUCLR = FKO*TORE(NI)*TORE(NJ)*EV/R
      ENDIF
      RETURN
C     *
C *** SECTION FOR MINDO/3 AND CNDO/2.
C     *
C *** RESONANCE INTEGRALS.
   20 CONTINUE
      CALL BETAIJ (NI,NJ,R,SIJ,T)
      CALL ROTBET (IA,JA,IORBS,JORBS,T,YY,H,LM4)
C *** COULOMB INTERACTIONS FOR MINDO/3 AND CNDO/2.
      IF(IOP.EQ.1) THEN
         RIJ = R*A0
         WIJ = W1/SQRT(RIJ**2+(W2/F0(NI)+W2/F0(NJ))**2)
         CALL BONTYP (NI,NJ,NBOND)
         SCALE = EXP(-ALPHA(NBOND)*RIJ)
         IF(NBOND.EQ.7 .OR. NBOND.EQ.11) SCALE=ALPHA(NBOND)*EXP(-RIJ)
         ENUCLR = TORE(NI)*TORE(NJ)*(WIJ+(W1/RIJ-WIJ)*SCALE)
      ELSE
         ZI  = ZS(NI)
         ZJ  = ZS(NJ)
         CALL CREP (NI,NJ,ZI,ZJ,R,WIJ)
         ENUCLR = TORE(NI)*TORE(NJ)*EV/R
      ENDIF
      W(2)   = WIJ
      DO 30 L=IA,IB
      M      = L*(L+1)/2
      H(M)   = H(M)-TORE(NJ)*WIJ
   30 CONTINUE
      DO 40 L=JA,JB
      M      = L*(L+1)/2
      H(M)   = H(M)-TORE(NI)*WIJ
   40 CONTINUE
      RETURN
      END