      SUBROUTINE BETAIJ (NI,NJ,R,S,T)
C     *
C     OVERLAP AND RESONANCE INTEGRALS IN LOCAL COORDINATES.
C     *
C     NOTATION. I=INPUT,O=OUTPUT.
C     NI,NJ     ATOMIC NUMBERS (I).
C     R         INTERNUCLEAR DISTANCE, IN ATOMIC UNITS (I).
C     S(14)     LOCAL OVERLAP INTEGRALS, DIMENSIONLESS (O).
C     T(14)     LOCAL RESONANCE INTEGRALS, IN EV (O).
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL NCASE1,NCASE2
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DPARM / LORBS(LMZ)
     ./DPARM1/ UDD(LMZ),ZD(LMZ),BETAD(LMZ)
     ./INOPT2/ IN2(300)
     ./PARMIN/ F0(18),VS(18),VP(18),BETA(55),ALPHA(55)
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ),BETAP(LMZ),
     .         ALP(LMZ)
      DIMENSION S(14),T(14)
C *** COMPUTE OVERLAP INTEGRALS.
      CALL OVERLP (NI,NJ,R,S,0,0)
C *** INITIALIZATION.
      IOP    = IN2(2)
      IORBS  = LORBS(NI)
      JORBS  = LORBS(NJ)
      NCASE1 = IORBS.EQ.1 .AND. JORBS.EQ.1
      NCASE2 = IORBS.EQ.1 .OR.  JORBS.EQ.1
C *** MNDO SECTION.
      IF(IOP.LE.0) THEN
         T(1) = PT5*(BETAS(NI)+BETAS(NJ))*S(1)
         IF(NCASE1) RETURN
         T(2) = PT5*(BETAS(NI)+BETAP(NJ))*S(2)
         T(3) = PT5*(BETAP(NI)+BETAS(NJ))*S(3)
         IF(.NOT.NCASE2) THEN
            TT   = PT5*(BETAP(NI)+BETAP(NJ))
            T(4) = TT*S(4)
            T(5) = TT*S(5)
         ENDIF
         IF(IORBS.GE.9 .OR. JORBS.GE.9) THEN
            T(6)  = PT5*(BETAD(NI)+BETAS(NJ))*S(6)
            T(7)  = PT5*(BETAS(NI)+BETAD(NJ))*S(7)
            IF(NCASE2) RETURN
            T(8)  = PT5*(BETAD(NI)+BETAP(NJ))*S(8)
            T(9)  = PT5*(BETAP(NI)+BETAD(NJ))*S(9)
            T(10) = PT5*(BETAD(NI)+BETAP(NJ))*S(10)
            T(11) = PT5*(BETAP(NI)+BETAD(NJ))*S(11)
            IF(IORBS.GE.9 .AND. JORBS.GE.9) THEN
               TT    = PT5*(BETAD(NI)+BETAD(NJ))
               T(12) = TT*S(12)
               T(13) = TT*S(13)
               T(14) = TT*S(14)
            ENDIF
         ENDIF
         RETURN
      ENDIF
C *** MINDO SECTION.
      IF(IOP.EQ.1) THEN
         CALL BONTYP (NI,NJ,NBOND)
         BET  = BETA(NBOND)
         T(1) = BET*(VS(NI)+VS(NJ))*S(1)
         IF(NCASE1) RETURN
         T(2) = BET*(VS(NI)+VP(NJ))*S(2)
         T(3) = BET*(VP(NI)+VS(NJ))*S(3)
         IF(NCASE2) RETURN
         TT   = BET*(VP(NI)+VP(NJ))
         T(4) = TT*S(4)
         T(5) = TT*S(5)
         RETURN
      ENDIF
C *** CNDO SECTION.
      IF(IOP.GE.2) THEN
         BET  = PT5*(BETAS(NI)+BETAS(NJ))
         IF(NI.GT.10 .OR. NJ.GT.10) BET=BET*0.75D0
         T(1) = BET*S(1)
         IF(NCASE1) RETURN
         T(2) = BET*S(2)
         T(3) = BET*S(3)
         IF(NCASE2) RETURN
         T(4) = BET*S(4)
         T(5) = BET*S(5)
         RETURN
      ENDIF
      END
