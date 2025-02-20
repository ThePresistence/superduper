      SUBROUTINE HHPAIR (J,I,COORD,SIJ,HIJ,WIJ,ENUCLR)
C     *
C     INTEGRALS FOR A HYDROGEN-HYDROGEN PAIR.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     J,I      ATOM PAIR I-J (I).
C     COORD    CARTESIAN COORDINATES IN ANGSTROM (I).
C     SIJ      TWO-CENTER OVERLAP INTEGRAL (O).
C     HIJ      TWO-CENTER RESONANCE INTEGRAL (O).
C     WIJ      TWO-CENTER TWO-ELECTRON INTEGRAL (SS,SS) (O).
C     ENUCLR   CONTRIBUTION TO CORE-CORE REPULSION (O).
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./MULTIP/ DD(6,0:LMZ),PO(9,0:LMZ)
     ./INOPT2/ IN2(300)
     ./PARMIN/ F0(18),VS(18),VP(18),BETA(55),ALPHA(55)
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ),BETAP(LMZ),
     .         ALP(LMZ)
      DIMENSION COORD(3,*)
C *** INPUT OPTION.
      IOP    = IN2(2)
      IPAROK = IN2(11)
C *** DISTANCE R (AU), AND RESONANCE INTEGRAL HIJ.
      R      = SQRT( (COORD(1,J)-COORD(1,I))**2
     1              +(COORD(2,J)-COORD(2,I))**2
     2              +(COORD(3,J)-COORD(3,I))**2 ) / A0
      ZR     = ZS(1)*R
      IF(ZR.LT.BIGEXP) THEN
         SIJ = EXP(-ZR)*(ONE+ZR+ZR*ZR/THREE)
         HIJ = BETAS(1)*SIJ
      ELSE
         SIJ = ZERO
         HIJ = ZERO
      ENDIF
C *** TWO-ELECTRON INTEGRAL WIJ AND CORE-CORE REPULSION ENUCLR.
      IF(IOP.LE.0) THEN
         WIJ    = EV/SQRT(R*R+FOUR*PO(1,1)**2)
         RIJ    = R*A0
         ENI    = EXP(-ALP(1)*RIJ)
         ENUCLR = WIJ*(ONE+ENI+ENI)
         IF(IOP.EQ.-2 .OR. IOP.EQ.-7 .OR. IOP.EQ.-12 .OR. IOP.EQ.-17)
     1      CALL REPAM1(1,1,RIJ,ENUCLR)
         IF(IPAROK.EQ.8) CALL PDDG(1,1,RIJ,ENUCLR)
      ELSE IF(IOP.EQ.1) THEN
         RIJ    = R*A0
         WIJ    = W1/SQRT(RIJ**2+FOUR*(W2/F0(1))**2)
         SCALE  = EXP(-ALPHA(1)*RIJ)
         ENUCLR = WIJ+(W1/RIJ-WIJ)*SCALE
      ELSE
         ZI     = ZS(1)
         CALL CREP (1,1,ZI,ZI,R,WIJ)
         ENUCLR = EV/R
      ENDIF
      RETURN
      END
