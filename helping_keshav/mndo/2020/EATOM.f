      FUNCTION EATOM (NI,ISS,IPP,IDD)
C     *
C     TOTAL ENERGY FOR AN ATOM WITH THE ATOMIC NUMBER NI
C     WHICH IS IN A CONFIGURATION WITH THE OCCUPATION NUMBERS
C     IS,IP,ID FOR THE S,P,D ORBITALS.
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATORB / IOS(LMZ),IOP(LMZ),IOD(LMZ)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DPARM / LORBS(LMZ)
     ./DPARM1/ UDD(LMZ),ZD(LMZ),BETAD(LMZ)
     ./DPARM3/ F0DD(LMZ),F2DD(LMZ),F4DD(LMZ),F0SD(LMZ),G2SD(LMZ),
     .         F0PD(LMZ),F2PD(LMZ),G1PD(LMZ),G3PD(LMZ)
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ),BETAP(LMZ),
     .         ALP(LMZ)
     ./REP   / GSS(LMZ),GPP(LMZ),GSP(LMZ),GP2(LMZ),HSP(LMZ),HPP(LMZ)
C    ./NBFILE/ NBF(20)
C     PREFACTORS FOR EXCHANGE TERMS.
      DIMENSION  CP(6),CD2(10),CD4(10)
      DATA CP / 0.D0, -1.D0, -3.D0,  -1.D0,   0.D0,   0.D0/
      DATA CD2/ 0.D0,-58.D0,-93.D0,-105.D0,-175.D0,-105.D0,-93.D0,
     1               -58.D0,  0.D0,   0.D0/
      DATA CD4/ 0.D0,  5.D0,-30.D0,-105.D0,-175.D0,-105.D0,-30.D0,
     1                 5.D0,  0.D0,   0.D0/
C *** INITIALIZATION.
      EATOM  = ZERO
      IF(ISS.LT.0 .OR. IPP.LT.0 .OR. IDD.LT.0) RETURN
      IF((ISS+IPP+IDD).EQ.0) THEN
         IS  = IOS(NI)
         IP  = IOP(NI)
         ID  = IOD(NI)
      ELSE
         IS  = ISS
         IP  = IPP
         ID  = IDD
      ENDIF
      E      = ZERO
C *** CONTRIBUTION FROM S ELECTRONS.
      IF(IS.EQ.1)THEN
         E   = E + USS(NI)
      ELSE IF(IS.EQ.2) THEN
         E   = E + TWO*USS(NI) + GSS(NI)
      ENDIF
C *** CONTRIBUTION FROM P ELECTRONS.
      IF(IP.GE.1) THEN
         E   = E + IP*UPP(NI) +IP*(IP-1)*PT5*GP2(NI) + CP(IP)*HPP(NI)
         IF(IS.EQ.1) THEN
           E = E + IP*GSP(NI) - MIN(IP,3)*HSP(NI)
         ELSE IF(IS.EQ.2)THEN
           E = E + IP*(TWO*GSP(NI)-HSP(NI))
         ENDIF
      ENDIF
C *** CONTRIBUTION FROM D ELECTRONS.
      IF(LORBS(NI).GE.9 .AND. ID.GE.1) THEN
         ADD = F0DD(NI)-(14.D0/441.D0)*(F2DD(NI)+F4DD(NI))
         E   = E + ID*UDD(NI) + ID*(ID-1)*PT5*ADD
     1           + (CD2(ID)*F2DD(NI)+CD4(ID)*F4DD(NI))/441.D0
         IF(IS.EQ.1) THEN
           E = E + ID*F0SD(NI) - MIN(ID,5)*G2SD(NI)/5.D0
         ELSE IF(IS.EQ.2)THEN
           E = E + ID*(TWO*F0SD(NI)-G2SD(NI)/5.D0)
         ENDIF
         IF(IP.GE.1)THEN
            APD = F0PD(NI)-(ONE/15.D0)*G1PD(NI)-(THREE/70.D0)*G3PD(NI)
            E   = E + IP*ID*APD
         ENDIF
      ENDIF
      EATOM  = E
C *** DEBUG PRINT.
C     NB6 = NBF(6)
C     WRITE(NB6,500) IS,IP,ID,NI,E
C 500 FORMAT(1X,'IS=',I2,'   IP=',I2,'  ID=',I2,'  NI=',I2,
C    1          '     EISOL=',F16.8,'EV')
      RETURN
      END
