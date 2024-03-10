      SUBROUTINE VALPOT (NI,NJ,CORE,S,T,VALPP,VALPP1)
C     *
C     VALENCE-VALENCE PSEUDOPOTENTIAL (IN EV).
C     EVALUATED FROM SECOND-ORDER FORMULAS
C     USING OVERLAP AND RESONANCE INTEGRALS.
C     *
C     NOTATION: I=INPUT, O=OUTPUT.
C     NI,NJ     ATOMIC NUMBERS (I)
C     CORE      LOCAL CORE-ELECTRON ATTRACTION INTEGRALS (I).
C     S         LOCAL OVERLAP INTEGRALS (I).
C     T         LOCAL RESONANCE INTEGRALS (I).
C     VALPP     PSEUDOPOTENTIAL, FIRST  TERM (O).
C     VALPP1    PSEUDOPOTENTIAL, SECOND TERM (O).
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DPARM / LORBS(LMZ)
     ./INOPT2/ IN2(300)
     ./PAROPT/ USS(LMZ),UPP(LMZ),PARDUM(LMZ*5)
      DIMENSION S(14),T(14)
      DIMENSION CORE(10,2),VALPP(10,2),VALPP1(10,2)
C *** INITIALIZATION.
      IOP    = IN2(2)
      IORBS  = LORBS(NI)
      JORBS  = LORBS(NJ)
      USI    = USS(NI) + CORE(1,1)
      USJ    = USS(NJ) + CORE(1,2)
      IF(IORBS.GE.4) THEN
         UPI = UPP(NI) + CORE(3,1)
         UPPI= UPP(NI) + CORE(4,1)
      ENDIF
      IF(JORBS.GE.4) THEN
         UPJ = UPP(NJ) + CORE(3,2)
         UPPJ= UPP(NJ) + CORE(4,2)
      ENDIF
C *** EVALUATION.
C     CASE H-H.
      IF(IORBS.EQ.1 .AND. JORBS.EQ.1) THEN
         VALPP(1,1) = -S(1)*T(1)
         VALPP(1,2) =  VALPP(1,1)
C     CASE X-H.
      ELSE IF(IORBS.GE.4 .AND. JORBS.EQ.1) THEN
         VALPP(1,1) = -S(1)*T(1)
         VALPP(1,2) = -S(1)*T(1)-S(3)*T(3)
         VALPP(3,1) = -S(3)*T(3)
         VALPP(2,1) = -PT5 * ( S(1)*T(3)+S(3)*T(1) )
         IF(IOP.NE.-8) THEN
           VALPP1(1,1)=  S(1)*S(1)*(USI-USJ)
           VALPP1(1,2)=  S(1)*S(1)*(USJ-USI) + S(3)*S(3)*(USJ-UPI)
           VALPP1(3,1)=  S(3)*S(3)*(UPI-USJ)
           VALPP1(2,1)=  PT5*(S(1)*S(3)*(USI+UPI-TWO*USJ))
         ENDIF
C     CASE H-X.
      ELSE IF(IORBS.EQ.1 .AND. JORBS.GE.4) THEN
         VALPP(1,1) = -S(1)*T(1)-S(2)*T(2)
         VALPP(1,2) = -S(1)*T(1)
         VALPP(3,2) = -S(2)*T(2)
         VALPP(2,2) = -PT5 * ( S(1)*T(2)+S(2)*T(1) )
         IF(IOP.NE.-8) THEN
           VALPP1(1,1)=  S(1)*S(1)*(USI-USJ) + S(2)*S(2)*(USI-UPJ)
           VALPP1(1,2)=  S(1)*S(1)*(USJ-USI)
           VALPP1(3,2)=  S(2)*S(2)*(UPJ-USI)
           VALPP1(2,2)=  PT5*(S(2)*S(1)*(USJ+UPJ-TWO*USI))
         ENDIF
         RETURN
C *** CASE X-Y.
      ELSE
         VALPP(1,1) = -S(1)*T(1)-S(2)*T(2)
         VALPP(1,2) = -S(1)*T(1)-S(3)*T(3)
         VALPP(3,1) = -S(3)*T(3)-S(4)*T(4)
         VALPP(3,2) = -S(2)*T(2)-S(4)*T(4)
         VALPP(4,1) = -S(5)*T(5)
         VALPP(4,2) =  VALPP(4,1)
         VALPP(2,1) = -PT5 * (S(1)*T(3)+S(3)*T(1)+S(2)*T(4)+S(4)*T(2))
         VALPP(2,2) = -PT5 * (S(1)*T(2)+S(2)*T(1)+S(3)*T(4)+S(4)*T(3))
         IF(IOP.NE.-8) THEN
           VALPP1(1,1)=  S(1)*S(1)*(USI-USJ) + S(2)*S(2)*(USI-UPJ)
           VALPP1(1,2)=  S(1)*S(1)*(USJ-USI) + S(3)*S(3)*(USJ-UPI)
           VALPP1(3,1)=  S(3)*S(3)*(UPI-USJ) + S(4)*S(4)*(UPI-UPJ)
           VALPP1(3,2)=  S(2)*S(2)*(UPJ-USI) + S(4)*S(4)*(UPJ-UPI)
           VALPP1(4,1)=  S(5)*S(5)*(UPPI-UPPJ)
           VALPP1(4,2)=  -VALPP1(4,1)
           VALPP1(2,1)=  PT5*(S(3)*S(1)*(USI+UPI-2*USJ)+S(2)*S(4)
     1                   *(USI+UPI-TWO*UPJ))
           VALPP1(2,2)=  PT5*(S(1)*S(2)*(USJ+UPJ-2*USI)+S(3)*S(4)
     1                   *(USJ+UPJ-TWO*UPI))
         ENDIF
      ENDIF
      RETURN
      END
