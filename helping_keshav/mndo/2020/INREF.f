      SUBROUTINE INREF (JPRINT,KSTOP,MOL,IFORM)
C     *
C     INPUT OF EXPERIMENTAL REFERENCE DATA.
C     *
C     NOTATION. I=INPUT,O=OUTPUT.
C     JPRINT   PRINTING FLAG (I).
C     KSTOP    ERROR FLAG (O).
C     MOL      NUMBER OF THE CURRENT MOLECULE (I).
C     IFORM    FLAG FOR INPUT FORMAT (I).
C     *
      USE LIMIT, ONLY: LM1,LMV,LMF,LMM,LMPR,LMMR,MXNICS,LMSTAT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*2  ELEMNT
      CHARACTER*49 FORM(25),FORMI
      CHARACTER*30 KOMEXP,KOMI,KOMI1
      CHARACTER*1  POPSPD(3)
      LOGICAL UHF
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./DFP   / X(LMV),NVAR
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./ELEMTS/ ELEMNT(107)
     ./INOPT1/ IN1(300)
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
     ./PAR4  / EXCIT(20),NEXCIT,IEXCIT(4,20)
     ./PAR21 / ID1(7,LMF),LITEXP(LMF)
     ./PAR22 / XEXP(LMF),ERREXP(LMF),WEXP(LMF)
     ./PAR23 / KOMEXP(LMF)
     ./PAR24 / NREF1(LMM),NREF2(LMM)
     ./PAR26 / JOPREF(LMM),JCOSMO(LMM)
     ./PAR44 / NCOEFS(LMF),IMOL(LMMR,LMF),ICOEFS(LMMR,LMF)
     ./PARM1 / A(3,LM1),NC(LM1),NB(LM1),NA(LM1),NN(LM1),NATOMS
     ./UHF   / UHF
      DIMENSION ID(7)
      DIMENSION IDZERO(LMPR)
      DIMENSION NEED(LMPR)
C *** FOR A GIVEN PROPERTY (I) THE IDENTIFIERS ID(KMIN)...ID(7) ARE NOT
C     NEEDED AND THEREFORE SET TO 0. KMIN IS DEFINED AS IDZERO(I).
      DATA IDZERO/ 8,6,7,8,5,8,5,5,5,5,
     1             8,8,8,5,5,6,5,6,4,5,
     2             4,4,6,6,5,4,8,6,6,6,
     3             6,5,7,6,7,7,7,7,7,8,
     4             5,4,5,5,5,4,4,8,8,0/
C *** MESSAGES FOR INPUT ERRORS.
      DATA FORM( 1)/' OUTSIDE THE ALLOWED RANGE: ID1'/
      DATA FORM( 2)/' OUTSIDE THE ALLOWED RANGE: ID4'/
      DATA FORM( 3)/' OUTSIDE THE ALLOWED RANGE: ID4 OR ID5'/
      DATA FORM( 4)/' OUTSIDE THE ALLOWED RANGE: ID4, ID5, OR ID6'/
      DATA FORM( 5)/' OUTSIDE THE ALLOWED RANGE: ID4, ID5, ID6, OR ID7'/
      DATA FORM( 6)/' XI MUST BE POSITIVE: BOND LENGTH'/
      DATA FORM( 7)/' XI MUST BE POSITIVE: BOND ANGLE'/
      DATA FORM( 8)/' XI MUST BE POSITIVE: IONIZATION POTENTIAL'/
      DATA FORM( 9)/' XI MUST BE POSITIVE: EXCITATION ENERGY'/
      DATA FORM(10)/' XI MUST BE POSITIVE: DIPOLE MOMENT'/
      DATA FORM(11)/' XI MUST BE POSITIVE: POLARISABILITY'/
      DATA FORM(12)/' XI MUST BE POSITIVE: VIBRATIONAL WAVENUMBER'/
      DATA FORM(13)/' XI MUST BE POSITIVE: POPULATION'/
      DATA FORM(14)/' XI MUST BE POSITIVE: S**2 VALUE'/
      DATA FORM(15)/' XI MUST BE POSITIVE: IP-DIFFERENCE'/
      DATA FORM(16)/' ERRI MUST NOT BE NEGATIVE (EXPERIMENTAL ERROR)'/
      DATA FORM(17)/' LITI MUST NOT BE NEGATIVE (LITERATURE REFERENCE)'/
      DATA FORM(18)/' AT MOST 20 EXCITATION ENERGIES ALLOWED'/
      DATA FORM(19)/' S**2 MUST ONLY BE USED IN UHF'/
      DATA FORM(20)/' ID4 MUST BE GREATER THAN ID5: IP-DIFFERENCE'/
      DATA FORM(21)/' WI MUST NOT BE NEGATIVE (GRADIENT WEIGHT)'/
      DATA FORM(22)/' NO VARIABLES TO BE OPTIMIZED (HENCE NO GRADIENT)'/
      DATA FORM(23)/' PROPERTY NOT CALCULATED'/
      DATA FORM(24)/' OUTSIDE THE ALLOWED RANGE: IMOL'/
      DATA FORM(25)/' XI MUST BE POSITIVE: DISTANCE RMSD'/
C *** LABELS FOR POPULATIONS.
      DATA POPSPD/'s','p','d'/
C *** FILE NUMBERS.
      NB5    = NBF(5)
      NB6    = NBF(6)
C *** FLAGS AND COUNTERS.
      NCOEFS = 0
      IMOL   = 0
      ICOEFS = 0
C *** CHECK ON ALLOWED RANGE OF MOL.
      IF(MOL.LT.1 .OR. MOL.GT.LMM) THEN
         WRITE(NB6,400) MOL,LMM
         STOP 'INREF'
      ENDIF
C *** INITIALIZATION.
      IF(MOL.EQ.1) THEN
         NR  = 0
      ELSE
         NR  = NREF2(MOL-1)
      ENDIF
      NREF1(MOL) = NR+1
      JOP    = IN1(3)
      KCI    = IN1(77)
      IATERG = IN2(119)
      IUVCD  = IN2(146)
      KSTOP  = 0
      NEXCIT = 0
      NIMAG  = 0
      DO 10 I=1,LMPR
      NEED(I)= 0
   10 CONTINUE
      IF(JPRINT.GE.5) WRITE(NB6,500)
C *** READ ONE LINE PER REFERENCE DATUM.
   20 CONTINUE
C     IF(IFORM.LE.0) THEN
         READ(NB5,510,END=200) ID,XI,ERRI,WI,LITI,KOMI
C     ELSE
C        READ(NB5,*,END=200) ID,XI,ERRI,WI,LITI,KOMI
C     ENDIF
      IF(ID(1).EQ.0) GO TO 300
C *** CHECK FOR INPUT ERRORS.
C     RANGE OF REFERENCE PROPERTIES.
      ID1ABS = ABS(ID(1))
      IF(ID1ABS.GT.LMPR) THEN
         KSTOP = KSTOP+1
         KCODE = 1
         GO TO 100
      ENDIF
C     ANY EXPERIMENTAL ERROR MUST BE POSITIVE.
      IF(ERRI.LT.ZERO) THEN
         KSTOP = KSTOP+1
         KCODE = 16
         GO TO 100
      ENDIF
C     ANY LITERATURE REFERENCE MUST BE POSITIVE.
      IF(LITI.LT.0) THEN
         KSTOP = KSTOP+1
         KCODE = 17
         GO TO 100
      ENDIF
C     HEAT OF FORMATION AND RELATIVE ENERGIES (ETC).
      IF(ID1ABS.EQ.1  .OR. (ID1ABS.GE.11 .AND. ID1ABS.LE.15) .OR.
     1   ID1ABS.EQ.41 .OR.  ID1ABS.EQ.43 .OR.  ID1ABS.EQ.45  .OR.
     2   ID1ABS.EQ.48 .OR.  ID1ABS.EQ.49) THEN
         IF(ID(4).NE.0) THEN
            JMAX = 7
            IF(ID1ABS.GE.14) JMAX=4
            DO 30 J=4,JMAX
            IF(ID1ABS.EQ.41 .AND. ID(J).GT.NUMAT) THEN
               KSTOP = KSTOP+1
               KCODE = 5
               GO TO 100
            ELSE IF((ID(J).GT.MOL .OR. ID(J).LE.(-MOL)) .AND.
     1              ID1ABS.NE.41) THEN
               KSTOP = KSTOP+1
               KCODE = 5
               GO TO 100
            ENDIF
   30       CONTINUE
         ENDIF
         IF(ID(4).EQ.0 .AND. NUMAT.EQ.1) THEN
            IF(ID(5).LT.0 .OR. ID(5).GT.3 .OR.
     1         ID(6).LT.0 .OR. ID(6).GT.6 .OR.
     2         ID(7).LT.0 .OR. ID(7).GT.10) THEN
               KSTOP = KSTOP+1
               KCODE = 5
               GO TO 100
            ENDIF
         ENDIF
C     REACTION ENERGIES.
      ELSE IF(ID1ABS.EQ.44) THEN
         IF(ID(4).EQ.0) THEN
            KSTOP = KSTOP+1
            KCODE = 2
            GO TO 100
         ENDIF
         NCOEFS(NR+1)   = 1
         IMOL(1,NR+1)   = MOL
         ICOEFS(1,NR+1) = ID(4)
   31    CONTINUE
         READ(NB5,570,END=200) I, J
         IF(I.NE.0) THEN
            IF(J.GT.MOL .OR. J.LE.(-MOL)) THEN
               KSTOP = KSTOP+1
               KCODE = 24
               GO TO 100
            ENDIF
            NCOEFS(NR+1)              = NCOEFS(NR+1)+1
            ICOEFS(NCOEFS(NR+1),NR+1) = I
            IMOL(NCOEFS(NR+1),NR+1)   = J
            GO TO 31
         ENDIF
C     BOND LENGTH.
      ELSE IF(ID1ABS.EQ.2) THEN
         IF(XI.LT.ZERO) THEN
            KSTOP = KSTOP+1
            KCODE = 6
            GO TO 100
         ENDIF
         IF(ID(4).LE.0 .OR. ID(4).GT.NATOMS .OR.
     1      ID(5).LE.0 .OR. ID(5).GT.NATOMS) THEN
            KSTOP = KSTOP+1
            KCODE = 3
            GO TO 100
         ENDIF
C     BOND ANGLE.
      ELSE IF(ID1ABS.EQ.3) THEN
         IF(XI.LT.ZERO) THEN
            KSTOP = KSTOP+1
            KCODE = 7
            GO TO 100
         ENDIF
         IF(ID(4).LE.0 .OR. ID(4).GT.NATOMS .OR.
     1      ID(5).LE.0 .OR. ID(5).GT.NATOMS .OR.
     2      ID(6).LE.0 .OR. ID(6).GT.NATOMS) THEN
            KSTOP = KSTOP+1
            KCODE = 4
            GO TO 100
         ENDIF
C     DIHEDRAL ANGLE.
      ELSE IF(ID1ABS.EQ.4) THEN
         IF(ID(4).LE.0 .OR. ID(4).GT.NATOMS .OR.
     1      ID(5).LE.0 .OR. ID(5).GT.NATOMS .OR.
     2      ID(6).LE.0 .OR. ID(6).GT.NATOMS .OR.
     3      ID(7).LE.0 .OR. ID(7).GT.NATOMS) THEN
            KSTOP = KSTOP+1
            KCODE = 5
            GO TO 100
         ENDIF
C     IONIZATION POTENTIAL.
      ELSE IF(ID1ABS.EQ.5) THEN
         IF(XI.LT.ZERO) THEN
            KSTOP = KSTOP+1
            KCODE = 8
            GO TO 100
         ENDIF
C     EXCITATION ENERGY.
      ELSE IF(ID1ABS.EQ.6) THEN
         IF(XI.LT.ZERO) THEN
            KSTOP = KSTOP+1
            KCODE = 9
            GO TO 100
         ENDIF
         IF(ID(4).LT.0 .OR.  ID(4).GT.NUMB  .OR.
     1      ID(5).LT.0 .OR.  ID(5).GT.NORBS .OR.
     2      ID(5).GT.0 .AND. ID(5).LE.NUMB) THEN
            KSTOP = KSTOP+1
            KCODE = 3
            GO TO 100
         ENDIF
         NEXCIT = NEXCIT+1
         IF(NEXCIT.GT.20) THEN
            KSTOP  = KSTOP+1
            KCODE  = 18
            NEXCIT = 20
            GO TO 100
         ENDIF
C     DIPOLE MOMENT.
      ELSE IF(ID1ABS.EQ.7 .OR. ID1ABS.EQ.25) THEN
         IF(XI.LT.ZERO) THEN
            KSTOP = KSTOP+1
            KCODE = 10
            GO TO 100
         ENDIF
         IF(ID(4).LT.0 .OR. ID(4).GT.3) THEN
            KSTOP = KSTOP+1
            KCODE = 2
            GO TO 100
         ENDIF
C     POLARISABILITY (ALPHA).
      ELSE IF(ID1ABS.EQ.8) THEN
         IF(XI.LT.ZERO) THEN
            KSTOP = KSTOP+1
            KCODE = 11
            GO TO 100
         ENDIF
C     VIBRATIONAL WAVENUMBER.
      ELSE IF(ID1ABS.EQ.16) THEN
         IF(XI.LT.ZERO) THEN
            KSTOP = KSTOP+1
            KCODE = 12
            GO TO 100
         ENDIF
         I3N = 3*NUMAT
         IF(ID(4).LT. 1 .OR. ID(4).GT.I3N .OR.
     1      ID(5).LT.-1 .OR. ID(5).GT.8) THEN
            KSTOP = KSTOP+1
            KCODE = 3
            GO TO 100
         ENDIF
C     ATOMIC CHARGE.
      ELSE IF(ID1ABS.EQ.17) THEN
         IF(ID(4).LT.1 .OR. ID(4).GT.NUMAT) THEN
            KSTOP = KSTOP+1
            KCODE = 2
            GO TO 100
         ENDIF
C     POPULATION.
      ELSE IF(ID1ABS.EQ.18) THEN
         IF(XI.LT.ZERO) THEN
            KSTOP = KSTOP+1
            KCODE = 13
            GO TO 100
         ENDIF
         IF(ID(4).LT.1 .OR. ID(4).GT.NUMAT .OR.
     1      ID(5).LT.1 .OR. ID(5).GT.3) THEN
            KSTOP = KSTOP+1
            KCODE = 3
            GO TO 100
         ENDIF
C     SPIN EXPECTATION VALUE.
      ELSE IF(ID1ABS.EQ.19) THEN
         IF(.NOT.UHF) THEN
            KSTOP = KSTOP+1
            KCODE = 19
            GO TO 100
         ENDIF
         IF(XI.LT.ZERO) THEN
            KSTOP = KSTOP+1
            KCODE = 14
            GO TO 100
         ENDIF
C     GRADIENT COMPONENTS OR NORMS.
      ELSE IF(ID1ABS.EQ.20) THEN
         IF(ID(4).LT.-1 .OR. ID(4).GT.NVAR) THEN
            KSTOP = KSTOP+1
            KCODE = 2
            GO TO 100
         ENDIF
         IF(NUMAT.LE.1 .OR. NVAR.LE.0) THEN
            KSTOP = KSTOP+1
            KCODE = 22
            GO TO 100
         ENDIF
C     GRADIENT NORMS.
      ELSE IF(ID1ABS.GE.21 .AND. ID1ABS.LE.22) THEN
         IF(WI.LT.ZERO) THEN
            KSTOP = KSTOP+1
            KCODE = 21
            GO TO 100
         ENDIF
         IF(NUMAT.LE.1 .OR. NVAR.LE.0) THEN
            KSTOP = KSTOP+1
            KCODE = 22
            GO TO 100
         ENDIF
C     HIGHER IONIZATION POTENTIAL.
      ELSE IF(ID1ABS.EQ.23) THEN
         IF(XI.LT.ZERO) THEN
            KSTOP = KSTOP+1
            KCODE = 8
            GO TO 100
         ENDIF
         IF(ID(4).LT. 0 .OR. ID(4).GT.NORBS .OR.
     1      ID(5).LT.-1 .OR. ID(5).GT.8) THEN
            KSTOP = KSTOP+1
            KCODE = 3
            GO TO 100
         ENDIF
C     DIFFERENCE OF IONIZATION POTENTIALS.
      ELSE IF(ID1ABS.EQ.24) THEN
         IF(XI.LT.ZERO) THEN
            KSTOP = KSTOP+1
            KCODE = 15
            GO TO 100
         ENDIF
         IF(ID(4).LT.0 .OR. ID(4).GT.NUMB .OR.
     1      ID(5).LT.0 .OR. ID(5).GT.NUMB) THEN
            KSTOP = KSTOP+1
            KCODE = 3
            GO TO 100
         ENDIF
         IF(ID(4).GT.0 .AND. ID(5).GT.0 .AND.
     1      ID(4).LT.ID(5)) THEN
            KSTOP = KSTOP+1
            KCODE = 20
            GO TO 100
         ENDIF
C     DIFFERENCE OF VIBRATIONAL WAVENUMBERS.
      ELSE IF(ID1ABS.EQ.27) THEN
         I3N = 3*NUMAT
         IF(ID(4).LT. 1 .OR. ID(4).GT.I3N .OR.
     1      ID(5).LT.-1 .OR. ID(5).GT.8   .OR.
     2      ID(6).LT. 1 .OR. ID(6).GT.I3N .OR.
     3      ID(7).LT.-1 .OR. ID(7).GT.8) THEN
            KSTOP = KSTOP+1
            KCODE = 5
            GO TO 100
         ENDIF
C     ADIABATIC SOLVATION ENERGIES.
      ELSE IF(ID1ABS.EQ.28) THEN
         IF(ID(4).GE.MOL .OR. ID(4).LE.(-MOL)) THEN
            KSTOP = KSTOP+1
            KCODE = 5
            GO TO 100
         ENDIF
         IF(ID(5).LT.0 .OR. ID(5).GT.4) THEN
            KSTOP = KSTOP+1
            KCODE = 3
            GO TO 100
         ENDIF
C     VERTICAL SOLVATION ENERGIES.
      ELSE IF(ID1ABS.EQ.29) THEN
         IF(ID(5).LT.-4 .OR. ID(5).GT.0) THEN
            KSTOP = KSTOP+1
            KCODE = 3
            GO TO 100
         ENDIF
C     NMR CHEMICAL SHIFTS.
      ELSE IF(ID1ABS.EQ.30 .OR. ID1ABS.EQ.31) THEN
         IF(ID(4).LE.0 .OR. ID(4).GT.NUMAT .OR.
     1      ID(5).LT.0 .OR. ID(5).GT.NUMAT) THEN
            KSTOP = KSTOP+1
            KCODE = 3
            GO TO 100
         ENDIF
C     NICS CHEMICAL SHIFTS.
      ELSE IF(ID1ABS.EQ.32) THEN
         IF(ID(4).LE.0 .OR. ID(4).GT.MXNICS) THEN
            KSTOP = KSTOP+1
            KCODE = 2
            GO TO 100
         ENDIF
C     GUGA-CI EXCITATION ENERGIES.
      ELSE IF(ID1ABS.EQ.33) THEN
         IF(ID(4).LE.0 .OR. ID(4).GT.LMSTAT .OR.
     1      ID(5).LT.0 .OR. ID(5).GT.8) THEN
            KSTOP = KSTOP+1
            KCODE = 3
            GO TO 100
         ENDIF
         IF(KCI.NE.5) THEN
            KSTOP = KSTOP+1
            KCODE = 23
            GO TO 100
         ENDIF
C     GUGA-CI ROTATIONAL STRENGTHS.
      ELSE IF(ID1ABS.EQ.34) THEN
         IF(ID(4).LE.0 .OR. ID(4).GT.LMSTAT .OR.
     1      ID(5).LT.0 .OR. ID(5).GT.8) THEN
            KSTOP = KSTOP+1
            KCODE = 3
            GO TO 100
         ENDIF
         IF(KCI.NE.5 .OR. IUVCD.LT.2) THEN
            KSTOP = KSTOP+1
            KCODE = 23
            GO TO 100
         ENDIF
C     GUGA-CI OSCILLATOR STRENGTHS.
      ELSE IF(ID1ABS.EQ.35) THEN
         IF(ID(4).LE.0 .OR. ID(4).GT.LMSTAT .OR.
     1      ID(5).LT.0 .OR. ID(5).GT.8 .OR.
     2      ID(6).LE.0 .OR. ID(6).GT.3) THEN
            KSTOP = KSTOP+1
            KCODE = 4
            GO TO 100
         ENDIF
         IF(KCI.NE.5 .OR. IUVCD.LT.2) THEN
            KSTOP = KSTOP+1
            KCODE = 23
            GO TO 100
         ENDIF
C     GUGA-CI DIPOLE AND TRANSITION MOMENT MAGNITUDES AND ANGLES.
      ELSE IF(ID1ABS.GE.36 .AND. ID1ABS.LE.38) THEN
         IF(ID(4).LE.0 .OR. ID(4).GT.LMSTAT .OR.
     1      ID(5).LT.0 .OR. ID(5).GT.8 .OR.
     2      ID(6).LT.0 .OR. ID(6).GT.3) THEN
            KSTOP = KSTOP+1
            KCODE = 4
            GO TO 100
         ENDIF
         IF(KCI.NE.5 .OR. IUVCD.LT.2) THEN
            KSTOP = KSTOP+1
            KCODE = 23
            GO TO 100
         ENDIF
C     GUGA-CI DIPOLE MOMENT AND COMPONENTS.
      ELSE IF(ID1ABS.EQ.39) THEN
         IF(ID(4).LE.0 .OR. ID(4).GT.LMSTAT .OR.
     1      ID(5).LT.0 .OR. ID(5).GT.8 .OR.
     2      ID(6).LT.0 .OR. ID(6).GT.3) THEN
            KSTOP = KSTOP+1
            KCODE = 4
            GO TO 100
         ENDIF
         IF(KCI.NE.5 .OR. IUVCD.LT.1) THEN
            KSTOP = KSTOP+1
            KCODE = 23
            GO TO 100
         ENDIF
C     ORBITAL ENERGIES AND THEIR DIFFERENCES.
      ELSE IF(ID1ABS.EQ.40) THEN
         IF(ID(4).LT. 0 .OR. ID(4).GT.NORBS .OR.
     1      ID(5).LT.-1 .OR. ID(5).GT.8 .OR.
     2      ID(6).LT. 0 .OR. ID(6).GT.NORBS .OR.
     3      ID(7).LT.-1 .OR. ID(7).GT.8) THEN
            KSTOP = KSTOP+1
            KCODE = 5
            GO TO 100
         ENDIF
C     DISTANCE RMSD.
      ELSE IF(ID1ABS.EQ.46) THEN
         IF(XI.LT.ZERO) THEN
            KSTOP = KSTOP+1
            KCODE = 25
            GO TO 100
         ENDIF
      ENDIF
C *** DEBUG PRINT OF CURRENT INPUT LINE.
      IF(JPRINT.GE.5) WRITE(NB6,520) ID,XI,ERRI,WI,LITI,KOMI(1:20)
C *** DEFINE DEFAULT VALUES FOR CERTAIN INPUT OPTIONS.
      IF(WI.EQ.ZERO) WI=ONE
C     DEFINE UNNEEDED IDENTIFIERS TO BE ZERO.
      KMIN   = IDZERO(ID1ABS)
      DO 40 J=KMIN,7
      ID(J)  = 0
   40 CONTINUE
C     HEAT OF FORMATION.
      IF(ID1ABS.EQ.1 .OR. (ID1ABS.GE.11 .AND. ID1ABS.LE.15)
     1   .OR. ID1ABS.EQ.48 .OR. ID1ABS.EQ.49) THEN
         DO 50 J=4,7
         IF(ID(J).LT.0) THEN
            ID(J) = MOL+ID(J)
         ENDIF
   50    CONTINUE
C     EXCITATION ENERGY.
      ELSE IF(ID1ABS.EQ.6) THEN
         IF(ID(4).EQ.0) ID(4)=NUMB
         IF(ID(5).EQ.0) ID(5)=NUMB+1
         IF(ID(6).EQ.0) ID(6)=1
         ID(7) = NEXCIT
         DO 60 J=1,4
         IEXCIT(J,NEXCIT) = ID(J+3)
   60    CONTINUE
         NMOS  = MAX(NMOS,ID(5))
C     ATOMIC CHARGE.
      ELSE IF(ID1ABS.EQ.17) THEN
         NI    = NAT(ID(4))
         KOMI  = ELEMNT(NI)
C     POPULATION.
      ELSE IF(ID1ABS.EQ.18) THEN
         NI    = NAT(ID(4))
         KOMI  = POPSPD(ID(5))
         KOMI(3:4) = 'at'
         KOMI(6:7) = ELEMNT(NI)
C     SPIN EXPECTATION VALUE.
      ELSE IF(ID1ABS.EQ.19) THEN
         SZ    = (NALPHA-NBETA)*PT5
         XI    = SZ*(SZ+ONE)
         ERRI  = ZERO
         LITI  = 0
         KOMI  = 'S**2'
C     GRADIENT COMPONENT AND GRADIENT NORM.
      ELSE IF(ID1ABS.GE.20 .AND. ID1ABS.LE.22) THEN
         XI    = ZERO
         ERRI  = ZERO
         LITI  = 0
         IF(ID1ABS.EQ.20) KOMI='GRADIENT'
         IF(ID1ABS.EQ.21) KOMI='GNORM'
         IF(ID1ABS.EQ.22) KOMI='CNORM'
C     HIGHER IONIZATION POTENTIALS.
      ELSE IF(ID1ABS.EQ.23) THEN
         IF(ID(4).EQ.0) THEN
            IF(ID(5).EQ.0) THEN
               ID(4) = NUMB
            ELSE
               ID(4) = 1
            ENDIF
         ENDIF
C     DIFFERENCE OF IONIZATION POTENTIALS.
      ELSE IF(ID1ABS.EQ.24) THEN
         IF(ID(4).EQ.0) ID(4) = NUMB
         IF(ID(5).EQ.0) ID(5) = NUMB-1
C     NUMBER OF IMAGINARY FREQUENCIES.
      ELSE IF(ID1ABS.EQ.26) THEN
         ERRI  = 0
         KOMI  = 'NIMAG'
         NIMAG = NINT(XI)
C     DIFFERENCE OF VIBRATIONAL WAVENUMBERS.
      ELSE IF(ID1ABS.EQ.27) THEN
         I3N   = 3*NUMAT
         IF(ID(4).EQ.0 .AND. ID(5).EQ.0) ID(4) = I3N
         IF(ID(6).EQ.0 .AND. ID(7).EQ.0) ID(6) = I3N-1
C     ADIABATIC SOLVATION ENERGIES.
      ELSE IF(ID1ABS.EQ.28) THEN
         IF(ID(4).LT.0) THEN
            ID(4) = MOL+ID(4)
         ENDIF
         IF(ID(5).EQ.0) THEN
            IF(IN1(28).GT.0 .AND. IN1(28).LE.4) THEN
               ID(5) = IN1(28)
            ELSE
               ID(5) = 1
            ENDIF
         ENDIF
         LCOSMO = ID(5)
C     VERTICAL SOLVATION ENERGIES.
      ELSE IF(ID1ABS.EQ.29) THEN
         IF(ID(5).EQ.0) THEN
            IF(IN1(28).GE.-4 .AND. IN1(28).LT.0) THEN
               ID(5) = IN1(28)
            ELSE
               ID(5) =-1
            ENDIF
         ENDIF
         LCOSMO = ID(5)
C     GUGA-CI EXCITATION ENERGIES.
      ELSE IF(ID1ABS.EQ.33) THEN
         IF(ID(6).LE.0) THEN
            ID(6) = MOL+ID(6)
         ENDIF
C     ORBITAL ENERGIES.
      ELSE IF(ID1ABS.EQ.40) THEN
         IF(ID(4).EQ.0 .AND. ID(6).EQ.0) THEN
            IF(ID(5).EQ.0) THEN
               ID(4) = NUMB
            ELSE
               ID(4) = 1
            ENDIF
         ENDIF
C     PROTON AFFINITY.
      ELSE IF(ID1ABS.EQ.43) THEN
         IF(ID(4).LT.0) THEN
            ID(4) = MOL+ID(4)
         ENDIF
C     REACTION ENERGY.
      ELSE IF(ID1ABS.EQ.44) THEN
         DO 80 I=1,NCOEFS(NR+1)
         IF(IMOL(I,NR+1).LT.0) THEN
            IMOL(I,NR+1) = MOL+IMOL(I,NR+1)
         ENDIF
   80    CONTINUE
C     ADIABATIC EXCITATION ENERGY.
      ELSE IF(ID1ABS.EQ.45) THEN
         IF(ID(4).LT.0) THEN
            ID(4) = MOL+ID(4)
         ENDIF
      ENDIF
C *** STORE INPUT DATA.
      CALL TXCOPY (KOMI,KOMI1,30)
      NR     = NR+1
      IF(NR.GT.LMF) THEN
         WRITE(NB6,410) LMF
         STOP 'INREF'
      ENDIF
      DO 70 I=1,7
      ID1 (I,NR) = ID(I)
   70 CONTINUE
      XEXP  (NR) = XI
      ERREXP(NR) = ERRI
      WEXP  (NR) = WI
      LITEXP(NR) = LITI
      KOMEXP(NR) = KOMI1
      NEED(ID1ABS) = 1
      GO TO 20
C *** TREATMENT OF INPUT ERRORS.
  100 CONTINUE
      IF(JPRINT.GT.0) THEN
         WRITE(NB6,520) ID,XI,ERRI,WI,LITI,KOMI(1:20)
         WRITE(NB6,530)
         FORMI  = FORM(KCODE)
         WRITE(NB6,540) FORMI
      ENDIF
      GO TO 20
C *** END-OF-FILE DURING INPUT.
  200 CONTINUE
      IF(JPRINT.GT.0) WRITE(NB6,550)
C *** GENERATE ADDITIONAL ENTRIES BY DEFAULT.
  300 CONTINUE
C     SPIN EXPECTATION VALUE.
      IF(UHF .AND. NEED(19).EQ.0 .AND. NR.LT.LMF) THEN
         SZ  = (NALPHA-NBETA)*PT5
         XI  = SZ*(SZ+ONE)
         NR  = NR+1
         DO 310 I=2,7
         ID1 (I,NR) = 0
  310    CONTINUE
         ID1 (1,NR) = 19
         XEXP  (NR) = XI
         ERREXP(NR) = ZERO
         WEXP  (NR) =-ONE
         LITEXP(NR) = 0
         KOMEXP(NR) = 'S**2'
         NEED(19)   = 1
      ENDIF
C     CARTESIAN GRADIENT NORM.
      IF(NEED(22).EQ.0 .AND. JOP.NE.-1 .AND. NR.LT.LMF) THEN
         NR  = NR+1
         DO 320 I=2,7
         ID1 (I,NR) = 0
  320    CONTINUE
         ID1 (1,NR) = 22
         XEXP  (NR) = ZERO
         ERREXP(NR) = ZERO
         WEXP  (NR) =-ONE
         LITEXP(NR) = 0
         KOMEXP(NR) = 'CNORM'
         NEED(22)   = 1
      ENDIF
C     NUMBER OF IMAGINARY FREQUENCIES.
      IF(NEED(26).EQ.0 .AND. NR.LT.LMF .AND.
     1  (NEED(16).EQ.1 .OR. NEED(27).EQ.1 .OR. JOP.EQ.3)) THEN
         NR  = NR+1
         DO 330 I=2,7
         ID1 (I,NR) = 0
  330    CONTINUE
         ID1 (1,NR) = 26
         XEXP  (NR) = ZERO
         ERREXP(NR) = ZERO
         WEXP  (NR) =-ONE
         LITEXP(NR) = 0
         KOMEXP(NR) = 'NIMAG'
         NEED(26)   = 1
      ENDIF
C *** SET MOLECULE-SPECIFIC OPTIONS (JOPREF,JCOSMO).
      JOPREF(MOL) = JOP
      IF((NEED(16).EQ.1 .OR. NEED(26).EQ.1 .OR. NEED(27).EQ.1)
     1   .AND. JOP.GE.0) THEN
         IF(NIMAG.LE.0) THEN
            JOPREF(MOL) = MAX(JOP,3)
         ELSE
            JOPREF(MOL) = MAX(JOP,4)
         ENDIF
      ENDIF
      JCOSMO(MOL) = IN1(28)
      IF(NEED(28).EQ.1 .OR. NEED(29).EQ.1) THEN
         JCOSMO(MOL) = LCOSMO
      ENDIF
C *** REQUEST THERMOCHEMICAL CALCULATIONS FOR THE ODM2 OR ODM3 METHOD
      IF(IATERG.EQ.-1 .AND. NEED(1).EQ.1 .AND. JOP.GE.0) THEN
         IF(NIMAG.LE.0) THEN
            JOPREF(MOL) = MAX(JOP,3)
         ELSE
            JOPREF(MOL) = MAX(JOP,4)
         ENDIF
      ENDIF
      IF(IATERG.EQ.-1 .AND. NEED(47).EQ.1 .AND. JOP.GE.0) THEN
         IF(NIMAG.LE.0) THEN
            JOPREF(MOL) = MAX(JOP,3)
         ELSE
            JOPREF(MOL) = MAX(JOP,4)
         ENDIF
      ENDIF
C *** REQUEST THERMOCHEMICAL CALCULATIONS FOR OLD METHODS
      IF(IATERG.EQ.1 .AND. NEED(42).EQ.1 .AND. JOP.GE.0) THEN
         IF(NIMAG.LE.0) THEN
            JOPREF(MOL) = MAX(JOP,3)
         ELSE
            JOPREF(MOL) = MAX(JOP,4)
         ENDIF
      ENDIF
C *** SET UPPER LIMIT AND EXIT.
      NREF2(MOL) = NR
      IF(JPRINT.GT.0 .AND. KSTOP.GT.0) WRITE(NB6,560) KSTOP
      RETURN
  400 FORMAT(///1X,'NUMBER OF CURRENT MOLECULE: MOL =',I5,
     1       /  1X,'OUTSIDE THE ALLOWED RANGE FROM 1-',I5,
     2       /  1X,'STOP IN SUBROUTINE INREF.'//)
  410 FORMAT(///1X,'TOO MANY REFERENCE DATA. MAXIMUM NUMBER =',I5,
     1       /  1X,'STOP IN SUBROUTINE INREF.'//)
  420 FORMAT(///1X,'TOO MANY REFERENCE REACTIONS. MAXIMUM NUMBER =',I5,
     1       /  1X,'STOP IN SUBROUTINE INREF.'//)
  430 FORMAT(///1X,'TOO MANY REFERENCE MOLECULES. MAXIMUM NUMBER =',I5,
     1       /  1X,'STOP IN SUBROUTINE INREF.'//)
  500 FORMAT(///1X,'INPUT FOR EXPERIMENTAL REFERENCE DATA.',
     1       /  1X,'COPY OF ORIGINAL INPUT IN SUBROUTINE INREF.',
     2       /  1X,'ID1 ID2 ID3 ID4 ID5 ID6 ID7     XI        ERRI',
     3          4X,'WI  LIT'/)
  510 FORMAT(   I2,6I3,2F10.5,F5.1,I5,A)
  520 FORMAT(   1X,I3,6I4,2F10.5,F5.1,I5,1X,A)
  530 FORMAT(   1X,'THE PRECEDING LINE OF INPUT IS IGNORED. REASON:')
  540 FORMAT(   1X,A)
  550 FORMAT(/  1X,'END-OF-FILE ENCOUNTERED IN SUBROUTINE INREF.')
  560 FORMAT(/  1X,'NUMBER OF LINES IGNORED BY SUBROUTINE INREF.',I5)
  570 FORMAT(   8X,2I3)
      END
