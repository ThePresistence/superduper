      SUBROUTINE START (LM5,LEN,ICALL)
C     *
C     INPUT SECTION FOR MOLECULAR DATA.
C     *
      USE LIMIT, ONLY: LM1, LMZ, LMM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CALFRG
      COMMON
     ./AMASS / AMS(LM1)
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./ELEMTS/ ELEMNT(107)
     ./INFOB / PC(3,LM1),PM(3),PR(3,3)
     ./INOPT2/ IN2(300)
     ./MOPAC / IMOPAC
     ./PAR2  / MOL
     ./PAR24 / NREF1(LMM),NREF2(LMM)
     ./PAR26 / JOPREF(LMM),JCOSMO(LMM)
     ./PAR41 / COORDA(3,LM1),LASTAT,NUMATA,NFRGSI,CALFRG
     ./REFATC/ REFCRD(3,LM1)
     ./SCPAR / SCRQ0(LMZ),SCRQ1(LMZ),SCRQ2(LMZ),SCNET(LMZ),SCNEB(LMZ)
C     *
C *** INPUT OF DATA.
C     INPUT OF MOLECULAR GEOMETRY.
      CALL INPUT (ICALL)
      IF(ICALL.EQ.-9) RETURN
      MOL    = MOL+1
C     GENERAL INPUT OPTIONS.
      IOP    = IN2(2)
      JOP    = IN2(3)
      IFORM  = IN2(5)
      INREFD = IN2(10)
      MMINP  = IN2(12)
      NMR    = IN2(13)
      JPRINT = IN2(42)
      KCI    = IN2(77)
      IMMDP  = IN2(30)
C     CHECK AND PRINT CURRENT PARAMETERS.
      IF(ICALL.GE.0) CALL ACHECK (ICALL)
      IF(ICALL.GE.0 .AND. JPRINT.GE.5) CALL PRTPAR (IOP,0,JPRINT)
C     INITIALIZE MMOK CORRECTION.
      CALL MMOKIN
C     INITIALIZE DISPERSION CORRECTION.
      IF(IMMDP.NE.-1) CALL MMDPI
C     CHECK FOR HYDROGEN BONDS (MNDO/H).
      CALL HBONDT (JPRINT)
C     INPUT OF CORRELATION DATA.
      IABSCI = ABS(KCI)
      IF(IABSCI.EQ.0) THEN
         CALL CHKLRO
      ELSE IF(IABSCI.EQ.1) THEN
         CALL SETUPS (IFORM,JPRINT)
      ELSE IF(IABSCI.GE.2 .AND. IABSCI.LE.4) THEN
         CALL INPERT (IFORM,JPRINT)
      ELSE IF(IABSCI.GE.5) THEN
         IF(IABSCI.EQ.5) CALL INGUGA (IFORM,JPRINT)
         IF(IABSCI.GT.5) CALL INES (IFORM,JPRINT)
         ICROSS = IN2(160)
         IF(ICROSS.EQ.1 .OR. ICROSS.EQ.2 .OR. ICROSS.EQ.7) THEN
            CALL INMULS (IFORM,JPRINT)
         ELSE IF(ICROSS.GE.3 .AND. ICROSS.LE.5) THEN
            CALL INCONI (IFORM,JPRINT)
         ELSE IF(ICROSS.EQ.6) THEN
            CALL INEXDY (IFORM,JPRINT)
         ENDIF
      ENDIF
C     INPUT OF EXPERIMENTAL REFERENCE DATA.
      IF(INREFD.NE.0) THEN
         CALL INREF  (JPRINT,KSTOP,MOL,IFORM)
         CALL OUTREF (JPRINT,KSTOP,MOL,0)
         IN2(3) = JOPREF(MOL)
         IN2(28)= JCOSMO(MOL)
      ENDIF
C     INPUT AND INITIALIZATION FOR COSMO.
C     ICOSMO = IN2(28)
      IF(IN2(28).NE.0) THEN
         IF(IN2(28).LE.4) THEN
            CALL INITSV (JPRINT)
         ELSE
C           SMOOTH COSMO INITIALIZATION. GMETRY NEEDS TO BE CALLED
C           IN CASE WE ARE READING CONSTRAINTS.
            CALL GMETRY (+1)
            REFCRD(1:3, 1:NUMAT) = COORD(1:3, 1:NUMAT)
            CALL INITSC(LM6,NUMAT,NAT,ELEMNT,SCRQ0,SCRQ1,SCRQ2,
     1                  SCNET,SCNEB,JPRINT)
         ENDIF
         IF(INREFD.NE.0) JCOSMO(MOL)=IN2(28)
      ENDIF
C     INPUT AND INITIALIZATION FOR EXTERNAL POINTS.
      IF(MMINP.GT.0) THEN
         CALL INQMMM (JPRINT)
      ENDIF
C     INPUT FOR NMR TREATMENT.
      IF(NMR.NE.0) THEN
         CALL INPNMR (JPRINT)
      ENDIF
C     CHECK FOR ANALYTICAL DERIVATIVES.
      CALL CHKANA (LEN)
C     *
C *** DYNAMIC MEMORY ALLOCATION.
      CALL DYNMEM (LEN,LM5)
C     *
C *** FURTHER INITIALIZATION AND PRINTING.
C     DEFINITION OF PAIR INDICES.
      CALL DYNINT
C     DEFINITION OF ATOMIC MASSES.
      IF(IMOPAC.LE.0) CALL ATMASS
C     PRINT MOMENTS OF INERTIA AND RELATED DATA.
      IF(JPRINT.GT.1 .AND. NUMAT.GT.1) THEN
         CALL INERT (AMS,COORD,PC,PM,PR,ITYPE,NUMAT,NAT,JPRINT)
      ENDIF
C     DEFINITION OF GAUSSIAN BASIS SET.
      IF((IOP.EQ.-5 .OR. IOP.EQ.-6 .OR. IOP.EQ.-8 .OR. IOP.EQ.-9
     1     .OR. IOP.EQ.-22.OR. IOP.EQ.-23).AND.ICALL.GE.0) THEN
         CALL GTOMIN (0,0,JPRINT)
         IF(IOP.EQ.-5) CALL ECPDEF (JPRINT)
      ELSE IF((ICROSS.EQ.2 .OR. ICROSS.EQ.6 .OR. ICROSS.EQ.7)
     2          .AND. ICALL.GE.0) THEN
C-AK     NON-ADIABATIC COUPLING MATRIX ELEMENTS ARE IMPLEMENTED
C-AK     USING GAUSSIAN BASIS SETS ONLY, AND CURRENTLY THE CODE
C-AK     DOES NOT SUPPORT DIFFERENT EXPONENTS FOR S AND P ORBITALS
C-AK     OF ONE ATOM. THEREFORE MNDO-TYPE METHODS ARE EXCLUDED.
C-AK     IF DIFFERENT EXPONENTS WILL BE SUPPORTED THE FUTURE,
C-AK     AN STO-6G BASIS MIGHT BE INITIALIZED HERE.
         WRITE(NB6,610) ICROSS, IOP
         STOP 'START'
      ENDIF
      IF(CALFRG) THEN
         IF(JOP.NE.-1) THEN
            ICALL = -9
            WRITE(NB6,600) JOP
            RETURN
         ENDIF
         IF(NFRGSI.EQ.1) THEN
            COORD               = ZERO
            COORD(1:3,1:LASTAT) = COORDA(1:3,1:LASTAT)
         ELSE IF(NFRGSI.EQ.2) THEN
            COORD               = ZERO
            COORD(1:3,1:NUMAT)  = COORDA(1:3,LASTAT+1:NUMATA)
            NREF2(MOL)          = NREF1(MOL)
         ENDIF
      ENDIF
      RETURN
  600 FORMAT(///1X,'ERROR IN SUBROUTINE START:'
     1       /  1X,'INVALID TYPE OF CALCULATION (JOP =',I5,
     2             ') FOR A FRAGMENT.'
     3       /  1X,'SINGLE-POINT CALCULATION (JOP=-1) EXPECTED.')
  610 FORMAT(///1X,'ERROR IN SUBROUTINE START:'
     1       /  1X,'ICROSS =',I2,' CANNOT BE USED WITH IOP =',I3,'.')
      END