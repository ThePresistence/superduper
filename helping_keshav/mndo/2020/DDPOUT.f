      SUBROUTINE DDPOUT (NI,JPRINT)
C     *
C     CHECK AND PRINT CHARGE SEPARATIONS AND ADDITIVE TERMS USED
C     TO COMPUTE THE TWO-CENTER TWO-ELECTRON INTEGRALS IN MNDO/D.
C     *
C     CHARGE SEPARATIONS DD(6,LMZ) IN ATOMIC UNITS.
C     ADDITIVE TERMS     PO(9,LMZ) IN ATOMIC UNITS.
C     SEE EQUATIONS (12)-(16) OF TCA PAPER FOR DD.
C     SEE EQUATIONS (19)-(26) OF TCA PAPER FOR PO.
C     *
C     INDEX OF DD AND PO : SS 1, SP 2, PP 3, SD 4, PD 5, DD 6.
C     MULTIPOLE          :  L=0,  L=1,  L=2,  L=2,  L=1,  L=2.
C     SPECIAL INDEX OF PO: PP 7, DD 8.
C     MULTIPOLE          :  L=0,  L=0.
C     FOR ATOMIC CORE    : ADDITIVE TERM PO(9,NI)
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*7 TEXT(9)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DPARM / LORBS(LMZ)
     ./DPARM3/ F0DD(LMZ),F2DD(LMZ),F4DD(LMZ),F0SD(LMZ),G2SD(LMZ),
     .         F0PD(LMZ),F2PD(LMZ),G1PD(LMZ),G3PD(LMZ)
     ./MULTIP/ DD(6,0:LMZ),PO(9,0:LMZ)
     ./NBFILE/ NBF(20)
     ./REP   / GSS(LMZ),GPP(LMZ),GSP(LMZ),GP2(LMZ),HSP(LMZ),HPP(LMZ)
      DATA TEXT/ 'SS    0', 'SP    1', 'PP    2', 'SD    2', 'PD    1',
     1           'DD    2', 'PP    0', 'DD    0', 'CORE  0'/
C *** GENERAL FORMULAS FOR ONE-CENTER MULTIPOLE-MULTIPOLE INTERACTIONS.
C     F1      DIPOLE-DIPOLE.
C     F2      QUADRUPOLE-QUADRUPOLE (SQUARE).
      F1(D1,P1) = (ONE/P1-ONE/SQRT(D1*D1+P1*P1))/FOUR
      F2(D2,P2) = (ONE/P2-TWO/SQRT(D2*D2+P2*P2)
     1                   +ONE/SQRT(TWO*D2*D2+P2*P2))/8.0D0
C *** SECTION FOR HYDROGEN.
      IF(JPRINT.LT.5) RETURN
      IF(LORBS(NI).EQ.1) THEN
         GSSNI  = PT5*EV/PO(1,NI)
C *** SECTION FOR HEAVIER ELEMENTS WITHOUT D ORBITALS.
      ELSE IF(LORBS(NI).EQ.4) THEN
         GSSNI  = PT5*EV/PO(1,NI)
         HSPNI  = EV*F1(DD(2,NI),PO(2,NI))
         HPPNI  = EV*F2(DD(3,NI),PO(3,NI))
C *** SECTION FOR HEAVIER ELEMENTS WITH D ORBITALS.
C     NOTE SPECIAL CONVENTION FOR DD(4,NI) AND DD(6,NI)
C     WHICH ARE STORED WITH AN EXTRA FACTOR OF SQRT(TWO).
C     THIS EXTRA FACTOR IS REMOVED AGAIN BELOW TO OBTAIN
C     VALUES (DD4SQ2,DD6SQ2) DEFINED AS IN THE TCA PAPER.
      ELSE IF(LORBS(NI).EQ.9) THEN
         DD4SQ2 = DD(4,NI)/SQRT(TWO)
         DD6SQ2 = DD(6,NI)/SQRT(TWO)
         GSSNI  = PT5*EV/PO(1,NI)
         HSPNI  = EV*F1(DD(2,NI),PO(2,NI))
         HPPNI  = EV*F2(DD(3,NI),PO(3,NI))
         HSDNI  = EV*F2(DD4SQ2  ,PO(4,NI))
         HPDNI  = EV*F1(DD(5,NI),PO(5,NI))
         HDDNI  = EV*F2(DD6SQ2  ,PO(6,NI))
         FPPNI  = PT5*EV/PO(7,NI)
         FDDNI  = PT5*EV/PO(8,NI)
         HSDREF = G2SD(NI)*0.2D0
         HPDREF = G1PD(NI)*FOUR /15.0D0
         HDDREF = F2DD(NI)*THREE/49.0D0
         FDDREF = F0DD(NI)
      ENDIF
C *** PRINTING SECTION.
      NB6    = NBF(6)
      WRITE(NB6,500) NI
      IF(LORBS(NI).EQ.1) THEN
         WRITE(NB6,510) TEXT(1),ZERO    ,PO(1,NI),GSSNI,GSS(NI)
      ELSE IF(LORBS(NI).EQ.4) THEN
         WRITE(NB6,510) TEXT(1),ZERO    ,PO(1,NI),GSSNI,GSS(NI)
         WRITE(NB6,510) TEXT(2),DD(2,NI),PO(2,NI),HSPNI,HSP(NI)
         WRITE(NB6,510) TEXT(3),DD(3,NI),PO(3,NI),HPPNI,HPP(NI)
      ELSE IF(LORBS(NI).EQ.9) THEN
         WRITE(NB6,510) TEXT(1),ZERO    ,PO(1,NI),GSSNI,GSS(NI)
         WRITE(NB6,510) TEXT(2),DD(2,NI),PO(2,NI),HSPNI,HSP(NI)
         WRITE(NB6,510) TEXT(3),DD(3,NI),PO(3,NI),HPPNI,HPP(NI)
         WRITE(NB6,520) TEXT(4),DD(4,NI),PO(4,NI),HSDNI,HSDREF,DD4SQ2
         WRITE(NB6,510) TEXT(5),DD(5,NI),PO(5,NI),HPDNI,HPDREF
         WRITE(NB6,520) TEXT(6),DD(6,NI),PO(6,NI),HDDNI,HDDREF,DD6SQ2
         WRITE(NB6,510) TEXT(7),ZERO    ,PO(7,NI),FPPNI,GSS(NI)
         WRITE(NB6,510) TEXT(8),ZERO    ,PO(8,NI),FDDNI,FDDREF
         WRITE(NB6,510) TEXT(9),ZERO    ,PO(9,NI)
      ENDIF
      RETURN
  500 FORMAT(///1X,'CHARGE SEPARATIONS DD AND ADDITIVE TERMS PO ',
     1             '(IN ATOMIC UNITS) FOR ELEMENT',I3,
     2       /  1X,'CORRESPONDING ONE-CENTER INTEGRALS (IN EV), ',
     3             'CALCULATED VS REFERENCE.',
     4       // 1X,'TYPE  L       DD',13X,'PO',13X,'CALC',11X,'REF'/)
  510 FORMAT(   1X,A,4F15.8)
  520 FORMAT(   1X,A,5F15.8,' (DD WITHOUT FACTOR SQRT2)')
      END
