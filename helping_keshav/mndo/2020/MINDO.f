      SUBROUTINE MINDO (H,W,LM4,LM9)
C     *
C     FULL CORE HAMILTONIAN FOR MINDO/3 AND CNDO/2.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     H(LM4)    CORE HAMILTONIAN MATRIX (O).
C     W(LM9)    TWO-ELECTRON INTEGRALS (O).
C     *
C     PARALLELIZED VERSION (FOR SHARED MEMORY).
C     CRAY MICROTASKING DIRECTIVES ARE INCLUDED.
C     NOTE THAT THE ONE-CENTER PART OF H(LM4) AND ENUCLR ARE GUARDED.
C     *
      USE LIMIT, ONLY: LM1, LMX, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./ENERGT/ EE,ENUCLR,EAT,ATHEAT
     ./EXTRA2/ SIJ(14),T(14),YY(675)
     ./INDEX / INDX(LMX)
     ./INOPT2/ IN2(300)
     ./PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./PARMIN/ F0(18),VS(18),VP(18),BETA(55),ALPHA(55)
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ),BETAP(LMZ),
     .         ALP(LMZ)
      DIMENSION H(LM4),W(LM9)
C$OMP THREADPRIVATE (/EXTRA2/)
C *** INITIALIZE SOME VARIABLES.
      IOP    = IN2(2)
      IJPAIR = INDX(NUMAT)+NUMAT
CMIC$ PROCESS
      ENUCLR = ZERO
      DO 10 I=1,LM4
      H(I)   = ZERO
   10 CONTINUE
      DO 20 I=1,IJPAIR
      W(I)   = ZERO
   20 CONTINUE
C
C *** DIAGONAL ONE-CENTER TERMS.
      DO 40 I=1,NUMAT
      NI     = NAT(I)
      IA     = NFIRST(I)
      LL     = INDX(IA)+IA
      H(LL)  = USS(NI)
      IORBS  = NLAST(I)-IA+1
      IF(IORBS.GE.4) THEN
         DO 30 J=IA+1,IA+3
         LL     = INDX(J)+J
         H(LL)  = UPP(NI)
   30    CONTINUE
      ENDIF
   40 CONTINUE
CMIC$ END PROCESS
C
C *** LOOP OVER ATOM PAIRS FOR OFFDIAGONAL TWO-CENTER TERMS.
C     ATOMS I AND J ARE IDENTIFIED AT THE BEGINNING OF THE LOOP.
CMIC$ DO GLOBAL FOR 64
C$OMP PARALLEL DO SCHEDULE(DYNAMIC)
C$OMP.PRIVATE(EN,I,IA,IB,IORBS,J,JA,JB,JORBS,NBOND,NI,NJ,
C$OMP.R,RIJ,SCALE,WIJ,ZI,ZJ)
      DO 70 IJ=2,IJPAIR
      I      = NINT(SQRT(TWO*IJ))
      J      = IJ-INDX(I)
      IF(I.EQ.J) GO TO 70
      NI     = NAT(I)
      NJ     = NAT(J)
      IA     = NFIRST(I)
      JA     = NFIRST(J)
      IB     = NLAST(I)
      JB     = NLAST(J)
      IORBS  = NLAST(I)-IA+1
      JORBS  = NLAST(J)-JA+1
C *** DISTANCE R (AU), ROTATION MATRIX, AND RESONANCE INTEGRALS.
      CALL ROTMAT (J,I,JORBS,IORBS,NUMAT,COORD,R,YY)
      CALL BETAIJ (NI,NJ,R,SIJ,T)
      CALL ROTBET (IA,JA,IORBS,JORBS,T,YY,H,LM4)
C *** COULOMB INTERACTIONS.
      IF(IOP.EQ.1) THEN
         RIJ = R*A0
         WIJ = W1/SQRT(RIJ**2+(W2/F0(NI)+W2/F0(NJ))**2)
         CALL BONTYP (NI,NJ,NBOND)
         SCALE = EXP(-ALPHA(NBOND)*RIJ)
         IF(NBOND.EQ.7 .OR. NBOND.EQ.11) SCALE=ALPHA(NBOND)*EXP(-RIJ)
         EN  = TORE(NI)*TORE(NJ)*(WIJ+(W1/RIJ-WIJ)*SCALE)
      ELSE
         ZI  = ZS(NI)
         ZJ  = ZS(NJ)
         CALL CREP (NI,NJ,ZI,ZJ,R,WIJ)
         EN  = TORE(NI)*TORE(NJ)*EV/R
      ENDIF
      W(IJ)  = WIJ
CMIC$ GUARD
C$OMP CRITICAL
      DO 50 L=IA,IB
      M      = INDX(L)+L
      H(M)   = H(M)-TORE(NJ)*WIJ
   50 CONTINUE
      DO 60 L=JA,JB
      M      = INDX(L)+L
      H(M)   = H(M)-TORE(NI)*WIJ
   60 CONTINUE
      ENUCLR = ENUCLR+EN
C$OMP END CRITICAL
CMIC$ END GUARD
   70 CONTINUE
C$OMP END PARALLEL DO
      RETURN
      END
