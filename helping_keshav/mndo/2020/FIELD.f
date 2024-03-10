      SUBROUTINE FIELD (A,LM5,ICALL,SCFCAL,NPRINT)
C     *
C     FINITE FIELD CALCULATION OF ELECTRIC RESPONSE PROPERTIES.
C     DIPOLE MOMENT, POLARIZABILITY, AND
C     FIRST AND SECOND HYPERPOLARIZABILITY.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     A(LM5)    AVAILABLE BUFFER (S).
C     ICALL     CONTROL AND ERROR FLAG (I,O).
C     SCFCAL    EXTERNAL ROUTINE FOR ENERGY EVALUATION (I).
C     NPRINT    PRINTING FLAG (I).
C     *
C     THE ZERO-FIELD DENSITY MATRIX IS SAVED ON FILE (CALL DENSAV).
C     STANDARD CASE: FOR EACH FINITE-FIELD CALCULATION, THE ZERO-FIELD
C     DENSITY MATRIX IS READ (CALL DENSAV) AND USED AS INITIAL GUESS.
C     EXTRAPOLATION: FOR HALF OF THE FINITE-FIELD CALCULATIONS (E.G. AT
C     FIELD +E(IC)), THE STANDARD PROCEDURE IS USED. FOR THE OTHER HALF
C     OF THE CALCULATIONS (E.G. AT FIELD -E(IC)), THE INITIAL GUESS IS
C     OBTAINED BY SUITABLE EXTRAPOLATION (CALL DENEXT) WHICH SPEEDS UP
C     SCF CONVERGENCE.
C     *
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL SCFCAL
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DIPOL / DD,DC(3),DMP(3)
     ./ENERGT/ EE,ENUCLR,EAT,ATHEAT
     ./FIELD2/ FIFI(3)
     ./FIELD3/ ALPHAD,ALPHAE,BETAD,BETAE,GAMMAD,GAMMAE
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./PARDER/ CORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
      DIMENSION A(LM5)
      DIMENSION DIPE(3),DIPD(3)
      DIMENSION AIID(3),AIIE(3),AIJD(6)
      DIMENSION BIIID(3),BIIIE(3),BIJJD(6),BIJJE(6)
      DIMENSION GIIIID(3),GIIIIE(3),GIIJJD(3),GIIJJE(6)
      DIMENSION D0(3),DM(3,3),DM2(3,3),DP(3,3),DP2(3,3)
      DIMENSION DMEME(3,3),DMEPE(3,3),DPEME(3,3),DPEPE(3,3)
      DIMENSION EMEME(3,3),EMEPE(3,3),EPEME(3,3),EPEPE(3,3)
      DIMENSION ENC(3)
      DIMENSION WP(3),WM(3)
C *** INPUT OPTIONS.
      IFLD1  = IN2(22)
      IFLD2  = IN2(23)
      KITSCF = IN2(71)
C *** INITIALIZE SOME CONVERSION FACTORS.
C     Dipole: a.u. to Debye
      AUTODB = 2.5416D0
C     Electric Field: a.u. to Volt/meter
      AUTOVM = 51.4257D0
C     Polarizability: a.u. to Angstrom**3
      AUTOA3 = A0*A0*A0
C     First  hyperpolarizability: a.u to e.s.u. (*10**(-30))
      BAUESU = 8.6571D-03
C     Second hyperpolarizability: a.u to e.s.u. (*10**(-36))
      GAUESU = 5.05116D-04
C *** INITIALIZE SOME VARIABLES FOR SCF CALCULATIONS.
      ICALL  = 40
      IN2(71)= 999
C     IF(IFAST  .EQ.0) IFAST  =-1
      IF(IN2(73).EQ.0) IN2(73)=-1
C *** SAVE THE REFERENCE DENSITY MATRIX AT ZERO FIELD.
      CALL DENSAV (-1,A,LM5,0)
C *** ENERGY AND DIPOLE MOMENT AT ZERO FIELD.
      W0     = (EE+ENUCLR)/EV
      D00    = DD/AUTODB
      DO 10 I=1,3
      D0(I)  = DC(I)/AUTODB
   10 CONTINUE
C *** DEFINE THE APPLIED ELECTRIC FIELD (IN ATOMIC UNITS).
      IF(IFLD2.NE.0) THEN
         EVAL = DBLE(IFLD1)*10.D0**(-IFLD2)
      ELSE
         EVAL = 1.0D-03
      ENDIF
      REVAL  = ONE/EVAL
C *** INTERACTION OF THE NUCLEI WITH THE FIELD (IN ATOMIC UNITS).
      DO 30 I=1,3
      EN     = ZERO
      DO 20 K=1,NUMAT
      ZN     = CORE(NAT(K))
      EN     = EN - EVAL*COORD(I,K)*ZN/A0
   20 CONTINUE
      ENC(I) = EN
   30 CONTINUE
C     *
C *** LOOP OVER ALL COMPONENTS OF THE APPLIED FIELD.
C *** ALL QUANTITIES IN ATOMIC UNITS.
      IAIJ   = 0
      DO 90 IC=1,3
      FIFI(1) = ZERO
      FIFI(2) = ZERO
      FIFI(3) = ZERO
C *** FIELD +E(IC)
      FIFI(IC) = EVAL
      CALL DENSAV (-1,A,LM5,1)
      CALL SCF (A,LM5,ICALL,SCFCAL)
      IF(ICALL.EQ.-1) RETURN
      WPE    = (EE+ENUCLR)/EV
      DPE    = DC(IC)/AUTODB
      DO 40 I=1,3
      DP(I,IC) = DC(I)/AUTODB
   40 CONTINUE
C *** FIELD -E(IC)
      FIFI(IC) = -EVAL
C     CALL DENSAV (-1,A,LM5,1)
      CALL DENEXT (-1,A,LM5)
      CALL SCF (A,LM5,ICALL,SCFCAL)
      IF(ICALL.EQ.-1) RETURN
      WME    = (EE+ENUCLR)/EV
      DME    = DC(IC)/AUTODB
      DO 50 I=1,3
      DM(I,IC) = DC(I)/AUTODB
   50 CONTINUE
C *** FIELD +2E(IC)
      FIFI(IC) = EVAL + EVAL
      CALL DENSAV (-1,A,LM5,1)
      CALL SCF (A,LM5,ICALL,SCFCAL)
      IF(ICALL.EQ.-1) RETURN
      WP2E   = (EE+ENUCLR)/EV
      DP2E   = DC(IC)/AUTODB
      DO 60 I=1,3
      DP2(I,IC) = DC(I)/AUTODB
   60 CONTINUE
C *** FIELD -2E(IC)
      FIFI(IC) = -EVAL-EVAL
C     CALL DENSAV (-1,A,LM5,1)
      CALL DENEXT (-1,A,LM5)
      CALL SCF (A,LM5,ICALL,SCFCAL)
      IF(ICALL.EQ.-1) RETURN
      WM2E   = (EE+ENUCLR)/EV
      DM2E   = DC(IC)/AUTODB
      DO 70 I=1,3
      DM2(I,IC) = DC(I)/AUTODB
   70 CONTINUE
C *** ADD FIELD-NUCLEI INTERACTION.
C     ENERGIES.
      EN     = ENC(IC)
      WPE    = WPE + EN
      WME    = WME - EN
      WP2E   = WP2E + EN + EN
      WM2E   = WM2E - EN - EN
      WP(IC) = WPE
      WM(IC) = WME
C     ENERGY DIFFERENCES AND SUMS.
      DWE    = WPE  - WME
      DW2E   = WP2E - WM2E
      SWE    = WPE  + WME
      SW2E   = WP2E + WM2E
C     DIPOLE MOMENT DIFFERENCES AND SUMS.
      DDE    = DPE  - DME
      DD2E   = DP2E - DM2E
      SDE    = DPE  + DME
      SD2E   = DP2E + DM2E
C *** Dipole moment
      DIPE(IC) = DW2E/12.D0 - (DWE + DWE)/THREE
      DIPE(IC) = DIPE(IC)*REVAL*AUTODB
      DIPD(IC) = (SDE + SDE)/THREE - SD2E/6.D0
      DIPD(IC) = DIPD(IC)*AUTODB
C *** Polarizability (diagonal part)
      AIIE(IC) = 5.D0*W0/TWO - FOUR*SWE/THREE + SW2E/12.D0
      AIIE(IC) = AIIE(IC)*REVAL*REVAL
      AIID(IC) = (DDE + DDE)/THREE - DD2E/12.D0
      AIID(IC) = AIID(IC)*REVAL
C *** First hyper-polarizability (diagonal part)
      BIIIE(IC) = DWE - PT5*DW2E
      BIIIE(IC) = BIIIE(IC)*REVAL*REVAL*REVAL
      BIIID(IC) = (SD2E - SDE)/THREE
      BIIID(IC) = BIIID(IC)*REVAL*REVAL
C *** Second hyper-polarizability (diagonal part)
      GIIIIE(IC) = FOUR*SWE - SW2E - 6.D0*W0
      GIIIIE(IC) = GIIIIE(IC)*REVAL*REVAL*REVAL*REVAL
      GIIIID(IC) = PT5*DD2E - DDE
      GIIIID(IC) = GIIIID(IC)*REVAL*REVAL*REVAL
C *** CALCULATION OF OFF-DIAGONAL TERMS.
      IF(IC.EQ.1) GO TO 90
      DO 80 JC=1,IC-1
C *** Order: yx, xy, zx, xz, zy, yz (controlled by IAIJ)
      IAIJ  = IAIJ + 1
C *** Calculation based on dipole moments.
C *** Polarizability Aij (lower triangle: yx, zx, zy)
C *** First hyperpolarizability Bijj (lower triangle: yxx, zxx, zyy)
      DIJD   = DP(IC,JC)  - DM(IC,JC)
      DIJD2  = DP2(IC,JC) - DM2(IC,JC)
      SIJD   = DP(IC,JC)  + DM(IC,JC)
      SIJD2  = DP2(IC,JC) + DM2(IC,JC)
      AIJ    = (DIJD + DIJD)/THREE - DIJD2/12.D0
      BIJJ   = (SIJD2 - SIJD)/THREE
      AIJD(IAIJ)  = AIJ*REVAL
      BIJJD(IAIJ) = BIJJ*REVAL*REVAL
C *** Polarizability Aij (upper triangle: xy, xz, yz)
C *** First hyperpolarizability Bijj (upper triangle: xyy, xzz, yzz)
      IAIJ   = IAIJ + 1
      DJID   = DP(JC,IC) - DM(JC,IC)
      DJID2  = DP2(JC,IC) - DM2(JC,IC)
      SJID   = DP(JC,IC) + DM(JC,IC)
      SJID2  = DP2(JC,IC) + DM2(JC,IC)
      AJI    = (DJID + DJID)/THREE - DJID2/12.D0
      BJII   = (SJID2 - SJID)/THREE
      AIJD(IAIJ)  = AJI*REVAL
      BIJJD(IAIJ) = BJII*REVAL*REVAL
C *** Interaction of nuclei with the field (a.u.)
      ENJ    = ENC(JC)
C *** FIELD +E(IC), +E(JC)
      FIFI(IC) = EVAL
      FIFI(JC) = EVAL
      FIFI(6-IC-JC) = ZERO
      CALL DENSAV (-1,A,LM5,1)
      CALL SCF (A,LM5,ICALL,SCFCAL)
      IF(ICALL.EQ.-1) RETURN
      WPEPE  = (EE+ENUCLR)/EV
      WPEPE  = WPEPE + EN + ENJ
      EPEPE(IC,JC) = WPEPE
      EPEPE(JC,IC) = WPEPE
      DPEPE(IC,JC) = DC(IC)/AUTODB
      DPEPE(JC,IC) = DC(JC)/AUTODB
C *** FIELD -E(IC), -E(JC)
      FIFI(IC) = -EVAL
      FIFI(JC) = -EVAL
C     CALL DENSAV (-1,A,LM5,1)
      CALL DENEXT (-1,A,LM5)
      CALL SCF (A,LM5,ICALL,SCFCAL)
      IF(ICALL.EQ.-1) RETURN
      WMEME  = (EE+ENUCLR)/EV
      WMEME  = WMEME - EN - ENJ
      EMEME(IC,JC) = WMEME
      EMEME(JC,IC) = WMEME
      DMEME(IC,JC) = DC(IC)/AUTODB
      DMEME(JC,IC) = DC(JC)/AUTODB
C *** FIELD +E(IC), -E(JC)
      FIFI(IC) = EVAL
      FIFI(JC) =-EVAL
      CALL DENSAV (-1,A,LM5,1)
      CALL SCF (A,LM5,ICALL,SCFCAL)
      IF(ICALL.EQ.-1) RETURN
      WPEME  = (EE+ENUCLR)/EV
      WPEME  = WPEME + EN - ENJ
      EPEME(IC,JC) = WPEME
      EMEPE(JC,IC) = WPEME
      DPEME(IC,JC) = DC(IC)/AUTODB
      DMEPE(JC,IC) = DC(JC)/AUTODB
C *** FIELD -E(IC), +E(JC)
      FIFI(IC) =-EVAL
      FIFI(JC) = EVAL
C     CALL DENSAV (-1,A,LM5,1)
      CALL DENEXT (-1,A,LM5)
      CALL SCF (A,LM5,ICALL,SCFCAL)
      IF(ICALL.EQ.-1) RETURN
      WMEPE  = (EE+ENUCLR)/EV
      WMEPE  = WMEPE - EN + ENJ
      EMEPE(IC,JC) = WMEPE
      EPEME(JC,IC) = WMEPE
      DMEPE(IC,JC) = DC(IC)/AUTODB
      DPEME(JC,IC) = DC(JC)/AUTODB
   80 CONTINUE
   90 CONTINUE
C     *
C *** FINAL CALCULATION OF OFF-DIAGONAL TERMS.
C     Giijj (dipole), and Bijj and Giijj (energy)
      IJ     = 0
      DO 110 I=1,2
      DEE    = DP(I,I) - DM(I,I)
      DO 100 J=I+1,3
      IJ     = IJ+1
C *** Giijj (dipole) (xxyy, xxzz, yyzz)
      DEIEJ  = DPEPE(I,J) - DMEPE(I,J) + DPEME(I,J) - DMEME(I,J)
      GIIJJD(IJ) = PT5*DEIEJ - DEE
      GIIJJD(IJ) = GIIJJD(IJ)*REVAL*REVAL*REVAL
C *** Giijj (energy) (xxyy, xxzz, yyzz)
      SEIEJ  = EPEPE(I,J) + EMEPE(I,J) + EPEME(I,J) + EMEME(I,J)
      SII    = WP(I) + WM(I)
      SJJ    = WP(J) + WM(J)
      GIIJJE(IJ) = TWO*(SII + SJJ) - SEIEJ - FOUR*W0
      GIIJJE(IJ) = GIIJJE(IJ)*REVAL*REVAL*REVAL*REVAL
C *** Bijj (energy) (yxx, xyy, zxx, xzz, zyy, yzz)
      DIJ    = -EPEPE(I,J) + EMEPE(I,J) - EPEME(I,J) + EMEME(I,J)
      DJI    = -EPEPE(J,I) + EMEPE(J,I) - EPEME(J,I) + EMEME(J,I)
      DII    = WP(I) - WM(I)
      DJJ    = WP(J) - WM(J)
      BIJJ   = PT5*DIJ + DII
      BJII   = PT5*DJI + DJJ
      BIJJE(2*IJ-1) = BJII*REVAL*REVAL*REVAL
      BIJJE(2*IJ)   = BIJJ*REVAL*REVAL*REVAL
  100 CONTINUE
  110 CONTINUE
      FIFI(1) = ZERO
      FIFI(2) = ZERO
      FIFI(3) = ZERO
C *** CALCULATE OBSERVABLE PROPERTIES
C *** Dipole
      ALPHAD = AUTOA3*(AIID(1) + AIID(2) + AIID(3))/THREE
      BETADX = BIIID(1) + BIJJD(2) + BIJJD(4)
      BETADY = BIIID(2) + BIJJD(1) + BIJJD(6)
      BETADZ = BIIID(3) + BIJJD(3) + BIJJD(5)
      BETAMU = BETADX*D0(1) + BETADY*D0(2) + BETADZ*D0(3)
      IF(ABS(D00).GT.1.0D-06) THEN
         BETAD = 0.6D0*BETAMU/D00*BAUESU
      ELSE
         BETAD = ZERO
      ENDIF
      GAMMA1 = GIIIID(1) + GIIIID(2) + GIIIID(3)
      GAMMA2 = GIIJJD(1) + GIIJJD(2) + GIIJJD(3)
      GAMMAD = 0.2D0*(GAMMA1 + GAMMA2 + GAMMA2)*GAUESU
C *** Energy
      ALPHAE = AUTOA3*(AIIE(1) + AIIE(2) + AIIE(3))/THREE
      BETADX = BIIIE(1) + BIJJE(2) + BIJJE(4)
      BETADY = BIIIE(2) + BIJJE(1) + BIJJE(6)
      BETADZ = BIIIE(3) + BIJJE(3) + BIJJE(5)
      BETAMU = BETADX*D0(1) + BETADY*D0(2) + BETADZ*D0(3)
      IF(ABS(D00).GT.1.0D0-6) THEN
         BETAE = 0.6D0*BETAMU/D00*BAUESU
      ELSE
         BETAE = ZERO
      ENDIF
      GAMMA1 = GIIIIE(1) + GIIIIE(2) + GIIIIE(3)
      GAMMA2 = GIIJJE(1) + GIIJJE(2) + GIIJJE(3)
      GAMMAE = 0.2D0*(GAMMA1 + GAMMA2 + GAMMA2)*GAUESU
C *** RESTORE THE INITIAL PARAMETERS.
      ICALL  = 0
      IN2(71)= KITSCF
C     IF(IFAST  .EQ.-1) IFAST  =0
      IF(IN2(73).EQ.-1) IN2(73)=0
C *** PRINTING SECTION.
      IF(NPRINT.GE.-5) THEN
         EVM    = EVAL*AUTOVM
         NB6    = NBF(6)
         WRITE(NB6,500) EVAL,EVM
         WRITE(NB6,510) DIPD,DIPE
         WRITE(NB6,600)
         WRITE(NB6,520)
         WRITE(NB6,610) AIID(1),AIJD(2),AIJD(4),AIIE(1)
         WRITE(NB6,620) AIJD(1),AIID(2),AIJD(6),AIIE(2)
         WRITE(NB6,630) AIJD(3),AIJD(5),AIID(3),AIIE(3)
         WRITE(NB6,640) ALPHAD,ALPHAE
         WRITE(NB6,700)
         WRITE(NB6,520)
         WRITE(NB6,710) BIIID(1),BIJJD(2),BIJJD(4),
     1                  BIIIE(1),BIJJE(2),BIJJE(4)
         WRITE(NB6,720) BIJJD(1),BIIID(2),BIJJD(6),
     1                  BIJJE(1),BIIIE(2),BIJJE(6)
         WRITE(NB6,730) BIJJD(3),BIJJD(5),BIIID(3),
     1                  BIJJE(3),BIJJE(5),BIIIE(3)
         WRITE(NB6,740) BETAD,BETAE
         WRITE(NB6,800)
         WRITE(NB6,520)
         WRITE(NB6,810) GIIIID(1),GIIJJD(1),GIIJJD(2),
     1                  GIIIIE(1),GIIJJE(1),GIIJJE(2)
         WRITE(NB6,820) GIIJJD(1),GIIIID(2),GIIJJD(3),
     1                  GIIJJE(1),GIIIIE(2),GIIJJE(3)
         WRITE(NB6,830) GIIJJD(2),GIIJJD(3),GIIIID(3),
     1                  GIIJJE(2),GIIJJE(3),GIIIIE(3)
         WRITE(NB6,840) GAMMAD,GAMMAE
      ENDIF
      RETURN
  500 FORMAT(///1X,'***** STATIC ELECTRIC PROPERTIES *****'
     .       // 1X,'ELECTRIC PROPERTIES ARE EVALUATED FROM FINITE ',
     .             'PERTURBATION THEORY',
     .       /  1X,'BY DIFFERENTIATION OF BOTH THE DIPOLE MOMENT AND ',
     .             'THE TOTAL ENERGY.',
     .       /  1X,'APPLIED ELECTRIC FIELD:',F10.6,' A.U. =',F10.6,
     .             ' VOLT/METER')
  510 FORMAT(// 1X,'DIPOLE MOMENT (DEBYES)',
     .       //11X,'X',12X,'Y',12X,'Z',
     .       /  1X,'D ',3(1X,F12.6),11X,'(EXPECTATION VALUES)',
     .       /  1X,'E ',3(1X,F12.6),11X,'(ENERGY DERIVATIVES)')
  520 FORMAT(/ 14X,'FROM DIPOLE CALCULATIONS',
     .         15X,'FROM ENERGY CALCULATIONS',
     .       //11X,'X',12X,'Y',12X,'Z',16X,'X',12X,'Y',12X,'Z')
  600 FORMAT(//1X,'POLARIZABILITY (A.U.)')
  610 FORMAT(  1X,'X ',3(1X,F12.6), 5X,F12.6)
  620 FORMAT(  1X,'Y ',3(1X,F12.6),18X,F12.6)
  630 FORMAT(  1X,'Z ',3(1X,F12.6),31X,F12.6)
  640 FORMAT(/ 1X,'MEAN POLARIZABILITY (A**3):',F14.6,' (DIPOLE)',
     .       / 1X,'MEAN POLARIZABILITY (A**3):',F14.6,' (ENERGY)')
  700 FORMAT(//1X,'FIRST HYPERPOLARIZABILITIES BIII AND BIJJ (A.U.)')
  710 FORMAT(  1X,'X ',3(1X,F12.6),4X,3(1X,F12.6))
  720 FORMAT(  1X,'Y ',3(1X,F12.6),4X,3(1X,F12.6))
  730 FORMAT(  1X,'Z ',3(1X,F12.6),4X,3(1X,F12.6))
  740 FORMAT(/ 1X,'FIRST HYPERPOLARIZABILITY (10D-30 ESU):',F19.6,1X,
     .            '(DIPOLE)',
     .       / 1X,'FIRST HYPERPOLARIZABILITY (10D-30 ESU):',F19.6,1X,
     .            '(ENERGY)')
  800 FORMAT(//1X,'SECOND HYPERPOLARIZABILITIES GIIII AND GIIJJ (A.U.)')
  810 FORMAT(  1X,'X ',3(1X,F12.4),4X,3(1X,F12.4))
  820 FORMAT(  1X,'Y ',3(1X,F12.4),4X,3(1X,F12.4))
  830 FORMAT(  1X,'Z ',3(1X,F12.4),4X,3(1X,F12.4))
  840 FORMAT(/ 1X,'SECOND HYPERPOLARIZABILITY (10D-36 ESU):',F18.6,1X,
     .            '(DIPOLE)',
     .       / 1X,'SECOND HYPERPOLARIZABILITY (10D-36 ESU):',F18.6,1X,
     .            '(ENERGY)')
      END
