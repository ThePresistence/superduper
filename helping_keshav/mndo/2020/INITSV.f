      SUBROUTINE INITSV (JPRINT)
C     *
C     INPUT AND INITIALIZATION OF COSMO VARIABLES.
C     *
      USE LIMIT, ONLY: LM1, LMZ, LMNPS, NPPA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL RADIUS,COSGEL
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./COSMO1/ EGAS,ESOLUT,ESOLV,EDISP,EDIE
     ./COSMO2/ NCOSMO,NPSG,NPSH,NPSTOT,LMNSET
     ./COSMO3/ FEPSI,RDS,DISEX2,SRAD(LM1)
     ./COSMO6/ CAVITA(LMZ),ATOMDI(LMZ),BONDDI(LMZ),BONDSP(LMZ),
     .         DELTA(LMZ),OMEGA(LMZ)
     ./INOPT1/ IN1(300)
     ./INOPT2/ IN2(300)
     ./INOPT3/ XN3(50)
     ./INOPT4/ XN4(50)
     ./MOPAC / IMOPAC
     ./NBFILE/ NBF(20)
      DIMENSION RBONDI(LMZ),REMS(LMZ),RVDW(LMZ)
C     VAN-DER-WAALS RADII (IN ANGSTROM) TAKEN FROM
C     A. BONDI, J.PHYS.CHEM. 68, 441 (1964), TABLE I.
      DATA RBONDI /
     1    1.20D0,   1.40D0,  -1.00D0,  -1.00D0,  -1.00D0,
     2    1.70D0,   1.55D0,   1.52D0,   1.47D0,   1.54D0,
     3   -1.00D0,  -1.00D0,  -1.00D0,   2.10D0,   1.80D0,
     4    1.80D0,   1.75D0,   1.88D0,  -1.00D0,  -1.00D0,
     5   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     6   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     7   -1.00D0,  -1.00D0,   1.85D0,   1.90D0,   1.85D0,
     8    2.02D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     9   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     A   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     B   -1.00D0,   2.06D0,   1.98D0,   2.16D0,  -1.00D0,
     C   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     D   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     E   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     F   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     G   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     H   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     I   -1.00D0 /
C     VAN-DER-WAALS RADII (IN ANGSTROM) TAKEN FROM
C     J. EMSLEY, THE ELEMENTS, OXFORD UNIVERSITY PRESS, OXFORD (1997).
C     http://www.shef.ac.uk/~chem/web-elements/
      DATA REMS /
     1    1.20D0,   1.22D0,  -1.00D0,  -1.00D0,   2.08D0,
     2    1.85D0,   1.54D0,   1.40D0,   1.35D0,   1.60D0,
     3    2.31D0,  -1.00D0,   2.05D0,   2.00D0,   1.90D0,
     4    1.85D0,   1.81D0,   1.91D0,   2.31D0,  -1.00D0,
     5   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     6   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     7   -1.00D0,  -1.00D0,   2.00D0,   2.00D0,   1.95D0,
     8    1.98D0,   2.44D0,  -1.00D0,  -1.00D0,  -1.00D0,
     9   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     A   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,   2.00D0,
     B    2.20D0,   2.20D0,   2.15D0,   2.16D0,   2.62D0,
     C   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     D   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     E   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     F   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     G   -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,  -1.00D0,
     H   -1.00D0,  -1.00D0,   2.40D0,  -1.00D0,  -1.00D0,
     I   -1.00D0 /
C     VAN-DER-WAALS RADII (IN ANGSTROM) TAKEN FROM
C     MOPAC93 AS IMPLEMENTED IN CAChe BY SAM COLE.
C     THESE ARE THE ORIGINAL VALUES USED BY ANDREAS KLAMT
C     IN THE CASE OF H, Li, C-F, Na, Al-Cl, K, AND Ca.
C     THE VALUES DIFFER FOR Br AND I (ORIGINALLY 1.80, 2.05).
C     DIMENSION RMOPAC(LMZ)
C     DATA RMOPAC /
C    1    1.08D0,  -1.00D0,   1.80D0,   1.45D0,   1.90D0,
C    2    1.53D0,   1.48D0,   1.36D0,   1.30D0,  -1.00D0,
C    3    2.30D0,   1.75D0,   2.05D0,   2.10D0,   1.75D0,
C    4    1.70D0,   1.65D0,  -1.00D0,   2.80D0,   2.75D0,
C    5    1.80D0,   1.80D0,   1.80D0,   1.80D0,   1.80D0,
C    6    1.80D0,   1.80D0,   1.80D0,   1.83D0,   1.92D0,
C    7    2.30D0,   2.19D0,   1.95D0,   1.90D0,   1.85D0,
C    8   -1.00D0,   2.50D0,   2.20D0,   2.00D0,   2.00D0,
C    9    2.00D0,   2.00D0,   2.00D0,   2.00D0,   2.00D0,
C    A    2.00D0,   2.03D0,   2.12D0,   2.40D0,   2.23D0,
C    B    2.10D0,   2.06D0,   1.98D0,  -1.00D0,   2.72D0,
C    C    2.40D0,   2.20D0,   2.20D0,   2.20D0,   2.20D0,
C    D    2.20D0,   2.20D0,   2.20D0,   2.20D0,   2.20D0,
C    E    2.20D0,   2.20D0,   2.20D0,   2.20D0,   2.20D0,
C    F    2.20D0,   2.20D0,   2.20D0,   2.20D0,   2.20D0,
C    G    2.20D0,   2.20D0,   2.20D0,   2.23D0,   2.32D0,
C    H    2.53D0,   2.35D0,   2.23D0,   2.19D0,   2.11D0,
C    I   -1.00D0 /
C *** FILE NUMBERS.
      NB5    = NBF(5)
      NB6    = NBF(6)
C *** GENERAL INPUT OPTIONS.
      IOP    = IN2(2)
      JOP    = IN2(3)
      IFORM  = IN2(5)
      IPAROK = IN2(11)
      ICOSMO = ABS(IN2(28))
      IGRAD  = IN2(199)
C *** READ SPECIAL INPUT OPTIONS FOR COSMO.
      IF(IMOPAC.EQ.0) THEN
         IF(ICOSMO.GT.2) THEN
C           READ(NB5,500) NSPA,NVDW,IPOT,NITRO,MODCSM,
C    1                  EPSI,RSOLV,DELSC,DISEX
            IF(IFORM.LE.0) THEN
               READ(NB5,500) (IN1(I),I=231,235),(XN3(J),J=16,19)
            ELSE
               READ(NB5,*)   (IN1(I),I=231,235),(XN3(J),J=16,19)
            ENDIF
         ELSE
            IN1(231:235) = 0
            XN3(16:19) = ZERO
         ENDIF
      ENDIF
      IN2(231:235) = IN1(231:235)
      XN4(16:19)   = XN3(16:19)
C *** IMPOSE GENERAL DEFAULT VALUES.
      IF(IN2(231).EQ.0)   IN2(231)= 42
      IF(XN4(16).EQ.ZERO) XN4(16) = 78.4D0
      IF(XN4(17).EQ.ZERO) XN4(17) =  1.0D0
      IF(XN4(18).EQ.ZERO) XN4(18) = XN4(17)
      IF(XN4(19).EQ.ZERO) XN4(19) =  2.0D0
      NSPA   = IN2(231)
      NVDW   = IN2(232)
      IPOT   = IN2(233)
      EPSI   = XN4(16)
      RSOLV  = XN4(17)
      DELSC  = XN4(18)
      DISEX  = XN4(19)
C *** INITIALIZE COSMO VARIABLES.
C     NCOSMO : COUNTER FOR COSMO CALCULATIONS.
C     FEPSI  : DIELECTRIC FACTOR.
C     RDS    : DISTANCE BETWEEN VAN-DER-WAALS SURFACE AND SAS.
C     EGAS   : HEAT OF FORMATION IN THE GAS PHASE.
      NCOSMO = 0
      FEPSI  = (EPSI-ONE)/(EPSI+PT5)
      RDS    = MAX(DELSC,0.1D0)
      EGAS   = ZERO
C *** DEFINE VAN-DER-WAALS RADII (IN ANGSTROM).
C     READ SPECIAL VALUES FOR VAN-DER-WAALS RADII (OPTIONAL).
      IF(NVDW.GT.0) THEN
         RVDW(1:LMZ) =-ONE
         DO 60 I=1,NVDW
         IF(IFORM.LE.0) THEN
            READ(NB5,510) IVDW,RIVDW
         ELSE
            READ(NB5,*)   IVDW,RIVDW
         ENDIF
         IF(IVDW.LE.LMZ) RVDW(IVDW)=RIVDW
   60    CONTINUE
C     USE EMSLEY VALUES FOR VAN-DER-WAALS RADII (OPTIONAL).
      ELSE IF(NVDW.EQ.-1) THEN
         RVDW(1:LMZ) = REMS(1:LMZ)
C     USE CAChe/MOPAC93 VALUES FOR VAN-DER-WAALS RADII (OPTIONAL).
      ELSE IF(NVDW.EQ.-2) THEN
         RVDW(1:LMZ) = REMS(1:LMZ)
C     USE BONDI VALUES FOR VAN-DER-WAALS RADII (DEFAULT).
      ELSE
         RVDW(1:LMZ) = RBONDI(1:LMZ)
      ENDIF
C *** CHECK FOR SIMPLE INPUT ERRORS.
      KSTOP  = 0
      IF(NSPA.LT.12 .OR. NSPA.GT.1082) THEN
         WRITE(NB6,515)
         WRITE(NB6,520)
         KSTOP = KSTOP+1
      ENDIF
      IF(IPOT.LT.0 .OR. IPOT.GT.5) THEN
         WRITE(NB6,520)
         KSTOP = KSTOP+1
      ENDIF
      IF(IPOT.GE.3 .AND. IPOT.LE.5) THEN
         IF(JOP.NE.-1 .AND. IGRAD.EQ.0) THEN
            WRITE(NB6,530)
            KSTOP = KSTOP+1
         ENDIF
         IF(IOP.NE.-2 .AND. IOP.NE.0 .AND. IPOT.GT.3) THEN
            WRITE(NB6,540)
            KSTOP = KSTOP+1
         ENDIF
      ENDIF
      IF(DELSC.LT.0.1D0) THEN
         WRITE(NB6,550)
C        KSTOP = KSTOP+1
      ENDIF
      IF(DELSC.GT.RSOLV+0.5D0) THEN
         WRITE(NB6,560)
         KSTOP = KSTOP+1
      ENDIF
C *** COMPUTE CONTROL VARIABLES.
      NHYD   = 0
      NATMAX = 0
      DO 100 I=1,NUMAT
      NATMAX = MAX(NAT(I),NATMAX)
      IF(NAT(I).EQ.1) NHYD=NHYD+1
  100 CONTINUE
C *** DEFINE PARAMETERS FOR SEMIEMPIRICAL ELECTROSTATIC POTENTIALS.
C     THESE PARAMETERS ARE AVAILABLE ONLY FOR MNDO AND AM1.
C     D. BAKOWIES AND W. THIEL, J.COMP.CHEM. 17, 87 (1996).
      DELTA(1:LMZ) = ZERO
      OMEGA(1:LMZ) = ZERO
      IF(IPOT.GE.3 .AND. IPOT.LE.6) THEN
         CALL INIPOT (DELTA,OMEGA,LMZ,IOP,IPOT)
         IF(IPAROK.GE.1 .AND. IPAROK.LE.3) CALL MODPOT
      ENDIF
C *** DEFINE CAVITATION AND DISPERSION TERMS FOR WATER IN AM1.
C     A. GELESSUS, PH.D. THESIS, UNIVERSITY OF ZURICH, 1997, PAGE 65.
C     THE FOLLOWING OPTIMIZED PARAMETERS ARE SPECIFIED IN THIS SECTION.
C     CAVITA(I) : ATOMIC CAVITATION PARAMETERS IN CAL*ANGSTROM**(-2)
C     ATOMDI(I) : ATOMIC SUPERDELOCALIZABILIES IN CAL*EV*ANGSTROM**(-2)
C     BONDDI(I) : BOND   SUPERDELOCALIZABILIES IN CAL*EV*ANGSTROM**(-2)
C     BONDSP(I) : BOND   SUPERDELOCALIZABILIES IN KCAL*EV*ANGSTROM**(-2)
C     BONDSP(I) : SPECIAL VALUES USED ONLY FOR N-H AND O-H BONDS.
C     NOTE DIFFERENT UNITS FOR BONDDI(I) AND BONDSP(I), FACTOR 1000.
      COSGEL = .FALSE.
      IF(ICOSMO.EQ.2 .OR. ICOSMO.EQ.4) THEN
         CAVITA(1:LMZ) = ZERO
         ATOMDI(1:LMZ) = ZERO
         BONDDI(1:LMZ) = ZERO
         BONDSP(1:LMZ) = ZERO
         IF(EPSI.GT.75.0D0 .AND. EPSI.LT.82.0D0) THEN
            IF(IOP.EQ.-2) THEN
               COSGEL = .TRUE.
               DO 160 I=1,NUMAT
               NI = NAT(I)
               IF((NI.NE.1).AND.(NI.NE.6).AND.(NI.NE.7).AND.
     1            (NI.NE.8).AND.(NI.NE.9).AND.(NI.NE.15).AND.
     2            (NI.NE.16).AND.(NI.NE.17).AND.(NI.NE.35).AND.
     3            (NI.NE.53)) COSGEL = .FALSE.
  160          CONTINUE
               IF(COSGEL) THEN
                  CAVITA(1) = 55.25D0
                  CAVITA(6) = 43.39D0
                  CAVITA(7) = 68.02D0
                  CAVITA(8) = 86.43D0
                  CAVITA(9) = 89.87D0
                  CAVITA(15)= 567.83D0
                  CAVITA(16)= 47.59D0
                  CAVITA(17)= 39.55D0
                  CAVITA(35)= 59.25D0
                  CAVITA(53)= 63.82D0
                  ATOMDI(1) = 108.18D0
                  ATOMDI(6) = 85.50D0
                  ATOMDI(7) = 74.99D0
                  ATOMDI(8) = 65.47D0
                  ATOMDI(9) = 53.10D0
                  ATOMDI(15)= 598.20D0
                  ATOMDI(16)= 68.97D0
                  ATOMDI(17)= 83.93D0
                  ATOMDI(35)= 101.91D0
                  ATOMDI(53)= 105.17D0
                  BONDDI(1) = 103.74D0
                  BONDDI(6) = 143.76D0
                  BONDDI(7) = 142.00D0
                  BONDDI(8) = 140.01D0
                  BONDDI(9) = 66.06D0
                  BONDDI(15)= 475.38D0
                  BONDDI(16)= 81.00D0
                  BONDDI(17)= 67.57D0
                  BONDDI(35)= 112.48D0
                  BONDDI(53)= 138.10D0
                  BONDSP(1) = 0.10374D0
                  BONDSP(7) = 27.993D0
                  BONDSP(8) = 150.417D0
               ELSE
                  WRITE(NB6,630)
               ENDIF
            ENDIF
         ELSE
            WRITE(NB6,640)
         ENDIF
      ENDIF
      IF(IOP.NE.-2 .AND. (ICOSMO.EQ.2 .OR. ICOSMO.EQ.4)) THEN
         WRITE(NB6,650)
      ENDIF
C *** DEFINE THE SOLVATION RADII SRAD (IN ANGSTROM).
C     FOR A GIVEN ATOM I, SRAD(I) IS THE SUM OF THE CORRESPONDING
C     VAN-DER-WAALS RADIUS RVDW(I) AND THE SOLVENT RADIUS RSOLV.
C     SRAD(I) THUS SPECIFIES THE SOLVENT-ACCESSIBLE SURFACE (SAS).
      DO 170 I=1,NUMAT
      NI = NAT(I)
      AVDW = RVDW(NI)
      IF(AVDW.LT.ZERO) THEN
         WRITE(NB6,660) NI
         KSTOP = KSTOP+1
      ELSE
         SRAD(I) = AVDW+RSOLV
      ENDIF
  170 CONTINUE
C *** PREPARE FOR THE CONSTRUCTION OF THE COSMO SURFACE.
C     THE NUMBER OF SEGMENTS (NPSG IN GENERAL, NPSH FOR HYDROGEN)
C     MUST BE OF THE FORM 10*(3**K)*(4**L)+2 (SEE DVFILL).
C     FOR A GIVEN INPUT VALUE OF THE NUMBER OF SEGMENTS (NSPA)
C     THE CODE DETERMINES THE NEAREST ALLOWED VALUE FOR NPSG
C     THAT IS SMALLER THAN OR EQUAL TO THE INPUT VALUE (NSPA).
C     THE CODE ALSO EVALUATES AN APPROPRIATE VALUE FOR NPSH.
C     THE FOLLOWING RESULTS ARE OBTAINED UP TO NSPA=1082:
C     FOR NSPA=  6-  31: I4=0, NPSH= 12, NPSG=  12 (K=0,L=0).
C     FOR NSPA= 32-  41: I4=0, NPSH= 12, NPSG=  32 (K=1,L=0).
C     FOR NSPA= 42-  91: I4=1, NPSH= 12, NPSG=  42 (K=0,L=1).
C     FOR NSPA= 92- 121: I4=1, NPSH= 32, NPSG=  92 (K=2,L=0).
C     FOR NSPA=122- 161: I4=1, NPSH= 42, NPSG= 122 (K=1,L=1).
C     FOR NSPA=162- 271: I4=2, NPSH= 42, NPSG= 162 (K=0,L=2).
C     FOR NSPA=272- 361: I4=2, NPSH= 92, NPSG= 272 (K=3,L=0).
C     FOR NSPA=362- 481: I4=2, NPSH=122, NPSG= 362 (K=2,L=1).
C     FOR NSPA=482- 641: I4=2, NPSH=162, NPSG= 482 (K=1,L=2).
C     FOR NSPA=642- 811: I4=3, NPSH=162, NPSG= 642 (K=0,L=3).
C     FOR NSPA=812-1081: I4=3, NPSH=272, NPSG= 812 (K=4,L=0).
C     FOR NSPA=    1082: I4=3, NPSH=362, NPSG=1082 (K=3,L=2).
      X0     = LOG(NSPA*0.1D0-0.199999D0)
      Z3     = LOG(THREE)
      Z4     = LOG(FOUR)
      I4     = INT(X0/Z4)
      NPSG   = 0
      DO 180 I=0,I4
      X      = X0-I*Z4
      N      = 3**INT(X/Z3)*4**I
      NPSG   = MAX(NPSG,N)
  180 CONTINUE
      NPSH   = NPSG/3
      IF(MOD(NPSG,3).NE.0) NPSH=NPSG/4
      NPSG   = 10*NPSG+2
      NPSH   = MAX(12,NPSH*10+2)
      NPSTOT = NHYD*NPSH + (NUMAT-NHYD)*NPSG
      LMNSET = NPPA*NUMAT
      DISEX2 = (FOUR*(1.5D0+RSOLV-RDS)*DISEX)**2/NSPA
C *** PRINTING SECTION.
      IF(JPRINT.GE.0) THEN
         WRITE(NB6,700)
         IF(ICOSMO.LE.2) THEN
            WRITE(NB6,710)
         ENDIF
         IF(ICOSMO.GT.2 .OR. JPRINT.GT.0) THEN
            WRITE(NB6,720) (IN2(I),I=231,235),(XN4(J),J=16,19)
         ENDIF
         WRITE(NB6,730)
         IF(IPOT.LT.3) THEN
            WRITE(NB6,740)
         ELSE IF(IPOT.EQ.3 .OR. IPOT.EQ.4) THEN
            WRITE(NB6,750)
         ELSE IF(IPOT.EQ.5) THEN
            WRITE(NB6,760)
         ENDIF
         DO 210 I=1,NATMAX
         RADIUS=.FALSE.
         DO 190 J=1,NUMAT
         IF(NAT(J).EQ.I) THEN
            RADIUS=.TRUE.
            GO TO 200
         ENDIF
  190    CONTINUE
  200    CONTINUE
         IF(RADIUS) THEN
            IF(IPOT.LT.3) THEN
               WRITE(NB6,770) I,RVDW(I)
            ELSE IF(IPOT.EQ.3 .OR. IPOT.EQ.4) THEN
               WRITE(NB6,780) I,RVDW(I),OMEGA(I)
            ELSE IF(IPOT.EQ.5) THEN
               WRITE(NB6,790) I,RVDW(I),OMEGA(I),DELTA(I)
            ENDIF
         ENDIF
  210    CONTINUE
         IF(COSGEL) THEN
            WRITE(NB6,800)
         ENDIF
         WRITE(NB6,810) NPSH,NPSG,NPSTOT
      ENDIF
C *** REDEFINE ICOSMO AND IN2(28) (IF NECESSARY).
C     IN2(28) WILL BE SET TO -2,-1, 1, OR 2.
      IF(ICOSMO.GT.4) THEN
         WRITE(NB6,670)
         KSTOP = KSTOP+1
      ELSE IF(ICOSMO.GT.2) THEN
         ICOSMO = ICOSMO-2
      ENDIF
C     IF(ICOSMO.EQ.2 .AND. .NOT.COSGEL) THEN
C        ICOSMO = 1
C     ENDIF
      IF(IN2(28).LT.0) THEN
         IN2(28) =-ICOSMO
      ELSE IF(IN2(28).GT.0) THEN
         IN2(28) = ICOSMO
      ENDIF
C *** ERROR EXIT.
      IF(KSTOP.GT.0) THEN
         STOP 'INITSV'
      ENDIF
      RETURN
  500 FORMAT(5I5,5X,5F10.5)
  510 FORMAT(I5,F10.5)
  515 FORMAT(/  1X,'INPUT OPTION NSPA MUST BE BETWEEN 12 AND 1082.'/)
  520 FORMAT(/  1X,'INPUT OPTION IPOT MUST BE BETWEEN 0 AND 5.'/)
  530 FORMAT(/  1X,'FOR SEMIEMPIRICAL ELECTROSTATIC POTENTIALS',
     1          1X,'(IPOT=3-5) THE ANALYTICAL COSMO GRADIENT IS NOT',
     2          1X,'YET AVAILABLE.',
     2       /  1X,'OPTION IGRAD=1 MAY BE USED TO DO THE CALCULATION',
     3          1X,'WITH SLOW NUMERICAL GRADIENT EVALUATIONS.'/)
  540 FORMAT(/  1X,'SPECIAL PARAMETERS FOR SEMIEMPIRICAL ELECTROSTATIC',
     1          1X,'POTENTIALS (IPOT=4,5) ARE DEFINED ONLY FOR MNDO',
     2          1X,'AND AM1.'/)
  550 FORMAT(/  1X,'DELSC TOO SMALL: SET TO 0.1.'/)
  560 FORMAT(/  1X,'DELSC UNREASONABLY LARGE.'/)
  630 FORMAT(/  1X,'WARNING WITH REGARD TO COSMO INPUT OPTIONS:',
     1       /  1X,'THE MOLECULE CONTAINS ATOMS WHICH DO NOT ALLOW AN',
     2          1X,'AM1 CALCULATION WITH DISPERSION AND CAVITATION.',
     3       /  1X,'ONLY ELECTROSTATIC CONTRIBUTIONS CAN BE COMPUTED.')
  640 FORMAT(/  1X,'WARNING WITH REGARD TO COSMO INPUT OPTIONS:',
     1       /  1X,'THE SOLVENT IS NOT WATER (DIFFERENT EPSILON).',
     2       /  1X,'ONLY ELECTROSTATIC CONTRIBUTIONS CAN BE COMPUTED.')
  650 FORMAT(/  1X,'WARNING WITH REGARD TO COSMO INPUT OPTIONS:',
     1       /  1X,'THE SEMIEMPIRICAL METHOD IS NOT AM1.',
     2       /  1X,'ONLY ELECTROSTATIC CONTRIBUTIONS CAN BE COMPUTED.')
  660 FORMAT(/  1X,'MISSING VAN DER WAALS RADIUS FOR ELEMENT NI=',I2/)
  670 FORMAT(/  1X,'THE ABSOLUTE VALUE OF ICOSMO MUST BE BETWEEN 1',
     1          1X,'AND 4.'/)
  700 FORMAT(///1X,'*** COSMO OPTIONS FOR THE MOLECULE ***')
  710 FORMAT(// 1X,'DEFAULT VALUES ARE USED FOR COSMO OPTIONS')
  720 FORMAT(// 1X,'NSPA   =',I7,  3X,'NVDW   =',I7,
     1          3X,'IPOT   =',I7,  3X,'NITRO  =',I7,
     2          3X,'MODCSM =',I7,
     3       /  1X,'EPSI   =',F7.3,3X,'RSOLV  =',F7.3,
     4          3X,'DELSC  =',F7.3,3X,'DISEX  =',F7.3)
  730 FORMAT(// 1X,'ATOMIC PARAMETERS FOR COSMO CALCULATION')
  740 FORMAT(/  3X,'ELEMENT',9X,'RADIUS (A)')
  750 FORMAT(/  3X,'ELEMENT',9X,'RADIUS (A)',8X,'OMEGA (A-1)')
  760 FORMAT(/  3X,'ELEMENT',9X,'RADIUS (A)',8X,'OMEGA (A-1)',
     1          8X,'DELTA (A)'/)
  770 FORMAT(   1X,I6,8X,F12.5)
  780 FORMAT(   1X,I6,2X,2(6X,F12.5))
  790 FORMAT(   1X,I6,2X,3(6X,F12.5))
  800 FORMAT(// 1X,'AN AM1-COSMO CALCULATION IN WATER IS PERFORMED',
     1       /  1X,'ELECTROSTATIC, DISPERSION AND CAVITATION TERMS',
     2          1X,'ARE CALCULATED')
  810 FORMAT(// 1X,'NUMBER OF SEGMENTS FOR HYDROGEN ATOMS      ',I10,
     1       /  1X,'NUMBER OF SEGMENTS FOR OTHER ATOMS         ',I10,
     2       /  1X,'MAXIMUM NUMBER OF SEGMENTS FOR THE MOLECULE',I10,
     3       /  1X,'HIDDEN SEGMENTS WILL BE REMOVED AUTOMATICALLY')
      END
