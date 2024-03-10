      SUBROUTINE CREP (NI,NJ,ZI,ZJ,RAB,W)
C     *
C     REPULSION INTEGRAL (SS,SS) FOR CNDO.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     NI,NJ     ATOMIC NUMBERS (I).
C     ZI,ZJ     ORBITAL EXPONENTS (I).
C     RAB       INTERNUCLEAR DISTANCE IN ATOMIC UNITS (I).
C     W         REPULSION INTEGRAL (SS,SS) IN EV (O).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
C *** INITIALIZE SOME VARIABLES.
      LI     = 1
      LJ     = 1
      IF(NI.GT.2) LI=2
      IF(NJ.GT.2) LJ=2
      IF(NI.GT.10 .OR. NJ.GT.10) STOP 'CREP'
C     SWITCH EXPONENTS, IF NECESSARY.
      IF(LI.LE.LJ) THEN
         NTYPE = (LJ*(LJ-1))/2+LI
         ZA  = ZI
         ZB  = ZJ
      ELSE
         NTYPE = (LI*(LI-1))/2+LJ
         ZA  = ZJ
         ZB  = ZI
      ENDIF
C     CALCULATE REDUCED VARIABLES.
      Z      = PT5*(ZA+ZB)
      T      = (ZA-ZB)/(ZA+ZB)
      P      = Z*RAB
      PA     = ZA*RAB
      EA     = EXP(-TWO*PA)
      IF(T.EQ.ZERO) GO TO 10
      V      = PT5*(T+ONE/T)
      PB     = ZB*RAB
      EB     = EXP(-TWO*PB)
      APV    = ONE+V
      AMV    = ONE-V
      APV2   = APV*APV
      AMV2   = AMV*AMV
C *** CALCULATE THE INTEGRALS.
C     (1S,1S/1S,1S)
      IF(NTYPE.EQ.1) THEN
         A   = Z/P*(ONE-AMV2*PT25*(TWO+V+PA)*EA-APV2*PT25*
     1              (TWO-V+PB)*EB)
         GO TO 20
      ENDIF
C     (1S,1S/2S,2S)
      AMV3   = AMV2*AMV
      PB2    = PB*PB
      PB3    = PB2*PB
      V2     = V*V
      V3     = V2*V
      IF(NTYPE.EQ.2) THEN
         A   = Z/P*(ONE-AMV3*0.0625D0*(ONE-5.0D0*V-FOUR*V2-TWO*V
     1              *PA)*EA-APV2*(0.0625D0*(15.0D0-22.0D0*V+15.0D0*V2
     2              -FOUR*V3)+0.375D0*(THREE-THREE*V+V2)*PB+PT25*
     3              (TWO-V)*PB2+PB3/12.0D0)*EB)
         GO TO 20
      ENDIF
C     (2S,2S/2S,2S)
      APV3   = APV2*APV
      PA2    = PA*PA
      PA3    = PA2*PA
      V4     = V3*V
      A      = Z/P*(ONE-AMV3*(0.0625D0*(8.0D0-V-27.0D0*V2-30.0D0*V3
     1              -10.0D0*V4)+0.03125D0*(11.0D0-19.0D0*V-44.0D0*V2
     2              -20.0D0*V3)*PA+0.0625D0*(ONE-5.0D0*V-FOUR*V2)*PA2
     3              -V*PA3/24.0D0)*EA-APV3*(0.0625D0*(8.0D0+V-27.0D0*V2
     4              +30.0D0*V3-10.0D0*V4)+0.03125D0*(11.0D0+19.0D0*V
     5              -44.0D0*V2+20.0D0*V3)*PB+0.0625D0*(ONE+5.0D0*V
     6              -FOUR*V2)*PB2+V*PB3/24.0D0)*EB)
      GO TO 20
C *** SIMPLIFIED CALCULATION OF THE INTEGRALS FOR T=0.
   10 CONTINUE
      P2     = P*P
      P3     = P2*P
      P4     = P3*P
      P5     = P4*P
      P6     = P5*P
      P7     = P6*P
      IF(NTYPE.EQ.1) THEN
         A   = Z/P*(ONE-(ONE+1.375D0*P+0.75D0*P2+P3/6.0D0)*EA)
      ELSE IF(NTYPE.EQ.2) THEN
         A   = Z/P*(ONE-(ONE+1.5625D0*P+1.125D0*P2+23.0D0*P3/48.0D0
     1              +0.125D0*P4+P5/60.0D0)*EA)
      ELSE
         A   = Z/P*(ONE-(ONE+419.0D0*P/256.0D0+163.0D0*P2/128.0D0
     1              +119.0D0*P3/192.0D0+5.0D0*P4/24.0D0+P5/20.0D0
     2              +P6/120.0D0+P7/1260.0D0)*EA)
      ENDIF
C *** CONVERT TO EV.
   20 W      = A*EV
      RETURN
      END
