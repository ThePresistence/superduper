      SUBROUTINE PARDEF (IOP,IPAROK,IEXBAS,NMR)
C     *
C     DEFINITION OF PARAMETERS IN THE WORKING ARRAYS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     IOP       CHOICE OF SEMIEMPIRICAL METHOD (I).
C     IPAROK    FLAG FOR NON-STANDARD PARAMETERS (I).
C     IEXBAS    FLAG FOR EXTRA POLARIZATION FUNCTIONS (I).
C     NMR       FLAG FOR NMR PARAMETERS (I).
C     *
C     INPUT COMMON BLOCKS (DEFINED IN BLOCKDATA ROUTINES).
C     THESE COMMON BLOCKS ARE NEVER CHANGED.
C     AM1       AM1 PARAMETERS.
C     AM1GAU    AM1 PARAMETERS.
C     AM1D      AM1/d PARAMETERS.
C     AMDGAU    AM1/d PARAMETERS.
C     CONSTF    CONVERSION FACTORS.
C     CONSTN    CONSTANT NUMBERS.
C     DNBND     MAIN QUANTUM NUMBERS.
C     DPARM7    DEFAULT FLAGS FOR SEPARATE TREATMENT OF F0SD (ETC).
C     EXPHAT    EXPERIMENTAL HEATS OF FORMATION OF THE ATOMS.
C     MNDOBL    MNDO PARAMETERS.
C     MNDOEL    MNDO PARAMETERS.
C     MNDOD     MNDO/d PARAMETERS.
C     OM1BL     OM1 PARAMETERS.
C     OM2BL     OM2 PARAMETERS.
C     OM3BL     OM3 PARAMETERS.
C     OM4BL     OM4 PARAMETERS.
C     ODM2BL    ODM2 PARAMETERS.
C     ODM3BL    ODM3 PARAMETERS.
C     PM3       PM3 PARAMETERS.
C     PM3GAU    PM3 PARAMETERS.
C     *
C     OUTPUT COMMON BLOCKS (DEFINED IN ROUTINES CALLED FROM PARDEF).
C     ABOND     BOND PARAMETERS FOR MNDO/d.
C     AMPGAU    GAUSSIAN CORE REPULSION TERMS FOR AM1 AND PM3.
C     DIPOL1    HYBRIDIZATION TERMS FOR DIPOLE MOMENTS.
C     DPARM     NUMBER OF BASIS ORBITALS PER ATOM.
C     DPARM1    STANDARD d ORBITAL PARAMETERS.
C     DPARM2    AUXILIARY EXPONENTS FOR INTEGRALS WITH d ORBITALS.
C     DPARM3    SLATER-CONDON PARAMETERS.
C     DPARM4    ONE-CENTER TWO-ELECTRON INTEGRALS (spd).
C     DPARM8    CURRENT FLAGS FOR SEPARATE TREATMENT OF F0SD (ETC).
C     MULTIP    CHARGE SEPARATIONS AND ADDITIVE TERMS.
C     PAR20     USER-DEFINED EXTERNAL PARAMETERS.
C     PARAVL    FLAGS FOR AVAILABILITY OF PARAMETERS.
C     PARDER    ENERGIES OF ATOMS (EHEAT,EISOL).
C     PARMIN    SPECIAL MINDO/3 PARAMETERS.
C     PAROMM    OM1/OM2 RELATED PARAMETERS.
C     PAROPT    STANDARD MNDO/AM1/PM3/etc PARAMETERS.
C     PSC       FACTORIALS AND BINOMIALS.
C     PSNPRD    NMR PARAMETERS.
C     REP       ONE-CENTER TWO-ELECTRON INTEGRALS (sp).
C     *
      USE LIMIT, ONLY: LMZ, LMP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     COMMON BLOCKS USED IN THIS ROUTINE.
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./INOPT2/ IN2(300)
     ./MULTIP/ DD(6,0:LMZ),PO(9,0:LMZ)
     ./PAR20 / XXPAR(LMP),NFPAR(LMP),NSPAR(LMP),NNPAR
     ./PARAVL/ IMPAR(LMZ)
C     COMMON BLOCKS USED IN LOWER-LEVEL ROUTINES.
C     COMMON
C    ./ABOND / ALPB(LMZ,LMZ),MALPB(LMZ)
C    ./AMPGAU/ GUESS1(LMZ*12),GSCAL(LMZ),IMP(LMZ)
C    ./DIPOL1/ HYF(LMZ*2)
C    ./DPARM / LORBS(LMZ)
C    ./DPARM1/ UDD(LMZ*3)
C    ./DPARM2/ ZSN(LMZ*3)
C    ./DPARM3/ F0DD(LMZ*9)
C    ./DPARM4/ REPD(52,LMZ)
C    ./DPARM7/ IF0SD(LMZ*4)
C    ./DPARM8/ JF0SD(LMZ*4)
C    ./PARDER/ TORE(LMZ*3)
C    ./PARMIN/ F0(164)
C    ./PAROMM/ OMM(LMZ,16)
C    ./PAROPT/ USS(LMZ*7)
C    ./PSC   / FBIN(930)
C    ./PSNPRD/ XBETAN(LMZ*18+3),IDEFN(LMZ+1)
C    ./REP   / GSS(LMZ*6)
      SAVE NCOUNT
      DATA NCOUNT/0/
C *** INPUT OPTIONS.
      ICOSMO = IN2(28)
C *** GENERAL INITIALIZATION.
      DO 10 I=1,LMZ
      IMPAR(I) = 0
   10 CONTINUE
      IOPPRT = 0
      NCOUNT = NCOUNT+1
      IF(NCOUNT.EQ.1) NNPAR=0
C *** INITIALIZATION FOR EXTERNAL POINT CHARGES.
      DO 20 I=1,6
      DD(I,0) = ZERO
   20 CONTINUE
      DO 30 I=1,9
      PO(I,0) = ZERO
   30 CONTINUE
C *** DEFAULT NUMBER OF ORBITALS PER ELEMENT.
      CALL NUMBAS
C *** DEFINITION OF STANDARD PARAMETERS.
      IF(IPAROK.EQ.-1 .AND. (IOP.EQ.-2 .OR. IOP.EQ.-7)) THEN
         CALL PARAM (0,0)
      ENDIF
      IF(IOP.LE.-10 .AND. IOP.NE.-22 .AND. IOP.NE.-23) THEN
         CALL FBINOM
         CALL PARAMD (IOP,IPAROK)
         CALL MLIG (IOP)
      ELSE IF(IOP.LE.0 .OR. IOP.EQ.-22 .OR. IOP.EQ.-23) THEN
         CALL PARAM (IOP,IPAROK)
      ELSE IF(IOP.EQ.1) THEN
         CALL MPARAM
      ELSE IF(IOP.EQ.2) THEN
         CALL CPARAM
      ELSE IF(IOP.EQ.5 .OR. IOP.EQ.6) THEN
         CALL DFTPAR
      ENDIF
C *** SPECIAL PARAMETERS FOR POLARIZATION FUNCTIONS.
      IF(IEXBAS.GT.0) THEN
         CALL FBINOM
         CALL PARAME
      ENDIF
C *** DEFAULT NMR PARAMETERS.
      CALL PARNMR (NMR)
C *** USER-DEFINED EXTERNAL PARAMETERS.
      IF(IPAROK.GE.1 .AND. IPAROK.LE.3) THEN
         IF(NCOUNT.EQ.1) CALL EXTERN (IPAROK,IOPPRT)
C        INITIALIZE THE SIMPLE DISPERSION PARAMETERS FOR SCOSMO.
         IF(ICOSMO.GE.5) CALL INISD
         CALL MODPAR (XXPAR,NFPAR,NSPAR,NNPAR)
C        THE INITIALIZATION OF THE TWO-BODY PARAMETERS
C        NEEDS TO BE INCLUDED AT THIS POINT.
         IF(NMR.GT.0) CALL PARNMR (999)
      ENDIF
      RETURN
      END