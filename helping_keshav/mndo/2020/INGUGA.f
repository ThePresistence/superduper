      SUBROUTINE INGUGA (IFORM,JPRINT)
C     *
C     INPUT OF OPTIONS FOR GUGA CONFIGURATION INTERACTION.
C     *
      USE LIMIT, ONLY: LM1, LMX, LMACT, LMREF, LMGRD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CICLO / CLOTHR
     ./CICONJ/ ICONJ(LM1)
     ./CIREFS/ ICIREF(LMACT,LMREF),JCIREF(LMACT,LMREF)
     ./CIMOS / IMOCI(LMX)
     ./CIMOSY/ JMOCI(LMX),MOCISY(LMX)
     ./CIROOT/ IROOTA(8)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./IJWORK/ IWORK(LMX),JWORK(LM1*3)
     ./INOPT1/ IN1(300)
     ./INOPT2/ IN2(300)
     ./MOPAC / IMOPAC
     ./NBFILE/ NBF(20)
     ./OCCFL / IMOCC,NOCCA,NOCCB,MSUB,MOSUMA,MOSUMB,MOCCA(8),MOCCB(8)
     ./OCCFL2/ DOMEGA,EFERMI,NFLOAT,NDOCC,NUMOCC
     ./OCCNM / OCCA(LMX),OCCB(LMX)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
      DIMENSION IN2SAV(300),ISWAP(2,5)
C *** FILE NUMBERS.
      NB5    = NBF(5)
      NB6    = NBF(6)
C *** GENERAL INPUT OPTIONS.
      IMULT  = IN2(66)
      IUHF   = IN2(70)
C
C *** IMPOSE DEFAULTS ON INPSCF OPTIONS RELEVANT TO GUGACI.
C     THIS NEEDS TO BE DONE AFTER INPUT HAS BEEN CALLED.
      IF(IUHF.EQ.-6) THEN
C        IN2(216)=NFLOAT.
         IF(IN2(216).LE.0 .OR. IN2(216).GT.NORBS) THEN
            IN2(216) = NORBS
         ENDIF
         NFLOAT  = IN2(216)
C        IN2(217)=NDOCC.
         IF(IN2(217).LE.0 .OR. IN2(217).GT.NALPHA+NBETA-NFLOAT) THEN
            IN2(217) = (NALPHA+NBETA-NFLOAT)/2
         ENDIF
         NDOCC   = IN2(217)
         NUMOCC  = NFLOAT + NDOCC
      ELSE
         NUMOCC  = NUMB
      ENDIF
C
C *** FIRST AND SECOND LINE: GENERAL OPTIONS.
      IF(IMOPAC.EQ.0) THEN
         IF(IFORM.LE.0) THEN
            READ(NB5,500) (IN1(I),I=131,150)
            READ(NB5,500) (IN1(I),I=151,170)
         ELSE
            READ(NB5,*)   (IN1(I),I=131,150)
            READ(NB5,*)   (IN1(I),I=151,170)
         ENDIF
      ENDIF
C
C *** CHECK FOR KEEPING PREVIOUS CI OPTIONS.
      IF(IN1(158).EQ.1) THEN
         IN2SAV(131:170) = IN2(131:170)
         WRITE(NB6,900) IN1(158)
      ENDIF
C
C *** SET CURRENT CI OPTIONS.
      IN2(131:170) = IN1(131:170)
C
C *** DEBUG PRINT.
      IF(JPRINT.GE.5) THEN
         WRITE(NB6,505)
         WRITE(NB6,520) IN2(131:170)
      ENDIF
C
C *** IMPOSE DEFAULT VALUES AND CHECK FOR SIMPLE ERRORS.
      KSTOP  = 0
C     ACTIVE OCCUPIED ORBITALS (ICI1).
      IF(IN2(131).LE.0) THEN
         IF(IUHF.EQ.-6) THEN
            IN2(131) = NFLOAT
         ELSE IF(IMULT.EQ.0 .OR. IMULT.EQ.2) THEN
            IN2(131) = 1
         ELSE
            IN2(131) = 2
         ENDIF
      ELSE IF(IN2(131).GT.NUMOCC) THEN
         IN2(131) = NUMOCC
         IF(JPRINT.GT.-5) WRITE(NB6,600) IN1(131),IN2(131)
      ENDIF
C     ACTIVE UNOCCUPIED ORBITALS (ICI2).
      IF(IN2(132).LE.0) THEN
         IF(IUHF.EQ.-6) THEN
            IN2(132) = 0
         ELSE IF(IMULT.EQ.0 .OR. IMULT.EQ.2) THEN
            IN2(132) = 1
         ELSE
            IN2(132) = 0
         ENDIF
      ELSE IF(IN2(132).GT.(NORBS-NUMOCC)) THEN
         IN2(132) = NORBS-NUMOCC
         IF(JPRINT.GT.-5) WRITE(NB6,610) IN1(132),IN2(132)
      ENDIF
C     TOTAL NUMBER OF ACTIVE ORBITALS (ICI1+ICI2=NACTIV).
      ICI1   = IN2(131)
      ICI2   = IN2(132)
      NACTIV = ICI1+ICI2
      IF(NACTIV.GT.LMACT) THEN
         KSTOP = KSTOP+1
         WRITE(NB6,620) NACTIV,LMACT
      ENDIF
C     NUMBER OF REFERENCE CONFIGURATIONS (NCIREF).
      IF(IN2(136).GT.LMREF) THEN
         KSTOP = KSTOP+1
         WRITE(NB6,630) IN2(136),LMREF
         IN2(136) = LMREF
      ENDIF
C     DEFINITION OF REFERENCE CONFIGURATIONS (NCIREF,MCIREF,LEVEXC).
      IF(IN2(136).EQ.0) THEN
         IN2(137) = 0
         IN2(138) = 0
         IF(JPRINT.GT.-5) THEN
            IF(IN1(137).GT.0) WRITE(NB6,640)
            IF(IN1(138).GT.0) WRITE(NB6,650)
         ENDIF
      ENDIF
C     DEFINITION OF EXCITATION LEVEL (NCIREF,LEVEXC).
      IF(IN2(136).GT.0) THEN
         IF(IN2(138).EQ.0) IN2(138) = 2
         IF(IN2(138).LT.0) IN2(138) = 0
      ENDIF
C     SPECIFIC ROOT OF INTEREST (LROOT).
      IF(IN2(140).EQ.0) THEN
         IN2(140) = 1
      ENDIF
C     NUMBER OF CI ROOTS (IROOT).
      IF(IN2(139).EQ.0) THEN
         IF(IN2(145).EQ.11) THEN
            IN2(139) = 1
         ELSE
            IN2(139) = MAX(1,IN2(140))
         ENDIF
      ENDIF
C     TOTAL CHARGE OF CI STATE (CICHG).
      IF(IN2(141).EQ.0) THEN
         IN2(141) = IN2(65)
      ELSE IF(IN2(141).EQ.9999) THEN
         IN2(141) = 0
      ENDIF
C     MULTIPLICITY (MULTCI).
      IF(IN2(142).EQ.0 .OR. IN2(142).LT.-1) THEN
         IF(IN2(141).NE.IN2(65)) THEN
            KSTOP = KSTOP+1
            WRITE(NB6,700)
         ENDIF
         IF(IMULT.GT.0) THEN
            IN2(142) = IMULT
         ELSE
            IN2(142) = 1
         ENDIF
      ENDIF
C     SYMMETRY (NCISYM).
      IF(IN2(139).LT.0 .AND. IN2(143).NE.0) THEN
         KSTOP = KSTOP+1
         WRITE(NB6,710)
      ENDIF
C     IN-CORE OR DIRECT CALCULATION (CIDIR).
C     THE DEFAULT IS SHAPE-DRIVEN IN-CORE CI.
      IF(IN2(144).EQ.0) THEN
         IN2(144) = 1
      ELSE IF(IN2(144).GT. 3) THEN
         IN2(144) =  3
      ELSE IF(IN2(144).LT.-3) THEN
         IN2(144) = -3
      ENDIF
C     METHOD OF DIAGONALIZATION (CIDIAG).
      IF(IN2(145).LT.0) THEN
         IN2(145) = 0
      ENDIF
      IF(IN2(145).EQ.11 .AND. IN2(139).NE.1) THEN
         KSTOP = KSTOP+1
         WRITE(NB6,720)
      ENDIF
C     THRESHOLD FOR PRINTING LEADING CONFIGURATIONS (CILEAD).
      IF(IN2(150).LE.0) THEN
         IN2(150) = 1000
      ENDIF
C     NUMBER OF ACTIVE OCCUPIED SPECIAL ORBITALS (JCI1).
      IF(IN2(151).LT.0) THEN
         IN2(151) = 0
      ELSE IF(IN2(151).GT.IN2(131)) THEN
         WRITE(NB6,605) IN2(151), IN2(131)
         IN2(151) = IN2(131)
      ENDIF
C     NUMBER OF ACTIVE VIRTUAL SPECIAL ORBITALS (JCI2).
      IF(IN2(152).LT.0) THEN
         IN2(152) = 0
      ELSE IF(IN2(152).GT.IN2(132)) THEN
         WRITE(NB6,615) IN2(152), IN2(132)
         IN2(152) = IN2(132)
      ENDIF
C     PI POPULATION THRESHOLD (PIPOP).
      IF(IN2(153).LT.0) THEN
         KSTOP = KSTOP + 1
         WRITE(NB6,725)
      ELSE IF(IN2(153).EQ.0) THEN
         IN2(153) = 4000
      ENDIF
C     THRESHOLD FOR AUTOMATIC SELECTION OF REFERENCES (CISELT).
      IF(IN2(155).LE.0) THEN
         IN2(155) = 85
      ENDIF
C     ORBITAL MAPPING IN CONSECUTIVE CI RUNS (IMOMAP)
C     OPTION IMOMAP=2 IS ONLY AVAILABLE IF IT IS SUPPORTED
C     BY THE CHOSEN OPTIMIZER.
C     IMOMAP=2 NOW SUPPORTED BY DYNAMICS (ICROSS=6)
      IF(IN2(156).EQ.2 .AND. IN2(160).NE.6) THEN
         IF(IN2(3).NE.0) THEN
            KSTOP = KSTOP + 1
            WRITE(NB6,750)
         ELSE IF(IN2(8).LT.0) THEN
            KSTOP = KSTOP + 1
            WRITE(NB6,760)
         ELSE IF(IN2(8).EQ.0 .AND. IN2(194).NE.0) THEN
            KSTOP = KSTOP + 1
            WRITE(NB6,770)
         ENDIF
      ENDIF
C     IMOMAP=3 IS INTENDED FOR USE WITH SINGLE-POINT CALCULATIONS,
C     E.G. AS PART OF THE INTERFACE TO CHEMSHELL.
      IF(IN2(156).EQ.3 .AND. IN2(3).GE.0) THEN
         KSTOP = KSTOP + 1
         WRITE(NB6,775)
      ENDIF
C     MAXIMUM NUMBER OF CI GRADIENTS TO BE COMPUTED (NCIGRD).
      IF(IN2(159).GT.LMGRD) THEN
         KSTOP = KSTOP+1
         WRITE(NB6,730) IN2(159),LMGRD
      ENDIF
C     ENFORCE EXPLICIT NCIGRD=2 FOR VARYING MULTIPLICITIES (MULTCI=-1)
      IF(IN2(142).EQ.-1 .AND. IN2(159).NE.2) THEN
         KSTOP = KSTOP+1
         WRITE(NB6,780)
      ENDIF
C     CONICAL INTERSECTIONS AND EXCITED STATE DYNAMICS (ICROSS).
      IF(IN2(160).LT.0 .OR. IN2(160).GT.7) THEN
         KSTOP = KSTOP+1
         WRITE(NB6,740) IN2(160)
      ENDIF
C     ONLY ALLOW VARYING MULTIPLICITIES (MULTCI=-1)
C     FOR ICROSS=1-5 CALCULATIONS
      IF(IN2(142).EQ.-1 .AND. (IN2(160).LT.1 .OR. IN2(160).GT.5)) THEN
         KSTOP = KSTOP+1
         WRITE(NB6,790)
      ENDIF
C     MINIMUM SUBSPACE DIMENSION IN DAVIDSON DIAGONALIZATION (MINDAV).
      IF(IN2(161).LT.0) THEN
         IN2(161) = 0
      ENDIF
C     MAXIMUM SUBSPACE DIMENSION IN DAVIDSON DIAGONALIZATION (MAXDAV).
      IF(IN2(162).LT.0) THEN
         IN2(162) = 0
      ENDIF
C     MAXIMUM NUMBER OF ITERATIONS IN DAVIDSON DIAGONALIZATION (KITDAV).
      IF(IN2(163).LE.0) THEN
         IF(IN2(145).EQ.12) THEN
            IN2(163) = 100
         ELSE
            IN2(163) = 0
         ENDIF
      ENDIF
C     CONVERGENCE CRITERION IN DAVIDSON DIAGONALIZATION (NRMDAV).
      IF(IN2(164).LE.0) THEN
         IN2(164) = 7
      ENDIF
C     MAXIMUM NUMBER OF FAILED MO MAPPING ATTEMPTS ALLOWED (MAXMAP).
      IF(IN2(165).LE.0) THEN
         IN2(165) = 5
      ENDIF
C     THRESHOLD FOR MO MAPPING OVERLAP CRITERION SUCCESS (MAPTHR).
      IF(IN2(166).LE.0) THEN
         IN2(166) = 90
      ENDIF
C     NUMBER OF CONJUGATED ATOMS (NCONJ)
      IF(IN2(167).LT.0) THEN
         IN2(167) = 0
      ENDIF
C
C *** THIRD AND FOURTH LINE: DEFINITION OF ACTIVE SPACE.
      MOVO   = IN2(134)
C     VALIDATION OF OPTION MOVO.
      IF(MOVO.LT.-5.OR.MOVO.EQ.2.OR.MOVO.GT.3) THEN
         KSTOP = KSTOP+1
         WRITE(NB6,625) MOVO
      ENDIF
C     DEFAULT ASSIGNMENT.
      IF(MOVO.LE.0) THEN
         DO 20 I=1,NACTIV
         JMOCI(I) = NUMOCC-ICI1+I
   20    CONTINUE
      ENDIF
      IF(MOVO.NE.3) THEN
         MOCISY(1:NACTIV) = 0
      ENDIF
C     READ NUMBERS OF ACTIVE ORBITALS.
      IF(MOVO.EQ.1.OR.MOVO.EQ.3) THEN
         IF(IFORM.LE.0) THEN
            READ(NB5,500) (JMOCI(I),I=1,NACTIV)
         ELSE
            READ(NB5,*)   (JMOCI(I),I=1,NACTIV)
         ENDIF
      ENDIF
C     READ SYMMETRY OF ACTIVE ORBITALS.
      IF(MOVO.EQ.3) THEN
         IF(IFORM.LE.0) THEN
            READ(NB5,500) (MOCISY(I),I=1,NACTIV)
         ELSE
            READ(NB5,*)   (MOCISY(I),I=1,NACTIV)
         ENDIF
      ENDIF
C
C *** DEFINE LAST RELEVANT ORBITAL.
      NMOS = NORBS
C
C *** FIFTH LINE: DEFINITION OF REFERENCE CONFIGURATIONS.
      MULSCF = MAX(1,IMULT)
      NCIREF = IN2(136)
      MCIREF = IN2(137)
      NCICHG = IN2(141)
      MULTCI = IN2(142)
C     THE OCCUPATION NUMBERS OF THE ACTIVE ORBITALS (IWORK) IN THE
C     SCF CONFIGURATION ARE DETERMINED THROUGH A CALL TO OCCSTD.
C     NELACT IS THE NUMBER OF ACTIVE ELECTRONS.
      CALL OCCSTD
      NELACT = IN2(65)-NCICHG
      DO 70 I=1,NACTIV
      IWORK(I) = NINT(OCCA(NUMOCC-ICI1+I)*TWO)
      NELACT = NELACT+IWORK(I)
   70 CONTINUE
C     IN THE CASE OF FULL CI, REDEFINE CONTROL VARIABLES (NCIREF,LEVEXC)
C     TO CONFORM TO THE CONVENTIONS IN THE GUGA-CI CODE.
      IF(IN2(136).EQ.0) THEN
         NCIREF   = 1
         IN2(136) = 1
         IN2(138) = NELACT
         IF(JPRINT.GT.-5) THEN
            WRITE(NB6,651)
            WRITE(NB6,652) IN2(138)
         ENDIF
      ENDIF
C     INITIALIZATION OF REFERENCE CONFIGURATIONS (JCIREF).
      IF(NCIREF.GT.0) THEN
         JCIREF(1:NACTIV,1:NCIREF) = 0
      ENDIF
C     DEFAULT ASSIGNMENTS.
      IF(NCIREF.GT.0 .AND. (MCIREF.EQ.0 .OR. MCIREF.EQ.3)) THEN
         NUMTWO = NELACT / 2
         NUMONE = NELACT - 2 * NUMTWO
         IF(NCIREF.EQ.1 .AND. MULSCF.EQ.MULTCI
     1                  .AND. NCICHG.EQ.IN2(65)) THEN
            DO I=1,NACTIV
               JCIREF(I,1) = IWORK(I)
            ENDDO
         ELSE IF(NCIREF.EQ.1 .AND. IN2(138).EQ.NELACT) THEN
            JCIREF(1:NUMTWO,       1) = 2
            JCIREF(NUMTWO+1:NACTIV,1) = 0
            IF(NUMONE.EQ.1)  JCIREF(NUMTWO+1,1) = 1
         ELSE IF(NCIREF.EQ.2 .AND. NUMONE.EQ.0) THEN
            JCIREF(1:NUMTWO,       1:2) = 2
            JCIREF(NUMTWO+1:NACTIV,1:2) = 0
            JCIREF(NUMTWO,           2) = 0
            JCIREF(NUMTWO+1,         2) = 2
         ELSE IF(NCIREF.EQ.3 .AND. NUMONE.EQ.0) THEN
            JCIREF(1:NUMTWO,       1:3) = 2
            JCIREF(NUMTWO+1:NACTIV,1:3) = 0
            JCIREF(NUMTWO:NUMTWO+1,  2) = 1
            JCIREF(NUMTWO,           3) = 0
            JCIREF(NUMTWO+1,         3) = 2
         ELSE IF(NCIREF.EQ.6 .AND. NUMONE.EQ.0) THEN
C     AN EXAMPLE FOR MORE COMPLICATED CASES.
            JCIREF(1:NUMTWO,       1:6) = 2
            JCIREF(NUMTWO+1:NACTIV,1:6) = 0
            JCIREF(NUMTWO,           2) = 0
            JCIREF(NUMTWO+1,         2) = 2
            JCIREF(NUMTWO-1,         3) = 0
            JCIREF(NUMTWO+1,         3) = 2
            JCIREF(NUMTWO-1,         4) = 1
            JCIREF(NUMTWO,           4) = 1
            JCIREF(NUMTWO+1,         4) = 2
            JCIREF(NUMTWO-1,         5) = 1
            JCIREF(NUMTWO+1,         5) = 1
            JCIREF(NUMTWO,           6) = 1
            JCIREF(NUMTWO+1,         6) = 1
         ELSE
            KSTOP = KSTOP+1
            WRITE(NB6,660) NCIREF,NCICHG,MULTCI
         ENDIF
      ENDIF
C     DIRECT INPUT OF REFERENCE CONFIGURATIONS.
      IF(NCIREF.GT.0 .AND. (MCIREF.EQ.1 .OR. MCIREF.EQ.4)) THEN
         DO 110 J=1,NCIREF
         IF(IFORM.LE.0) THEN
            READ(NB5,500) (JCIREF(I,J),I=1,NACTIV)
         ELSE
            READ(NB5,*)   (JCIREF(I,J),I=1,NACTIV)
         ENDIF
  110    CONTINUE
      ENDIF
C     INPUT OF EXCITATION INDICES.
      IF(NCIREF.GT.0 .AND. MCIREF.EQ.2) THEN
         IF(JPRINT.GE.5) THEN
            WRITE(NB6,550)
         ENDIF
         DO 140 J=1,NCIREF
         DO 120 I=1,NACTIV
         JCIREF(I,J) = IWORK(I)
  120    CONTINUE
         IF(IFORM.LE.0) THEN
            READ(NB5,500) (ISWAP(1,K),ISWAP(2,K),K=1,5)
         ELSE
            READ(NB5,*)   (ISWAP(1,K),ISWAP(2,K),K=1,5)
         ENDIF
         IF(JPRINT.GE.5) THEN
            WRITE(NB6,560) J,(ISWAP(1,K),ISWAP(2,K),K=1,5)
         ENDIF
         DO 130 K=1,5
         IF(ISWAP(1,K).EQ.0 .AND. ISWAP(2,K).EQ.0) GO TO 130
         I1 = ISWAP(1,K)-(NUMOCC-ICI1)
         I2 = ISWAP(2,K)-(NUMOCC-ICI1)
         IF(I1.GT.0    .AND. I1.LE.ICI1 .AND.
     1      I2.GT.ICI1 .AND. I2.LE.NACTIV) THEN
            JCIREF(I1,J) = JCIREF(I1,J) - 1
            JCIREF(I2,J) = JCIREF(I2,J) + 1
         ELSE
            KSTOP = KSTOP+1
            WRITE(NB6,670) J,K,ISWAP(1,K),ISWAP(2,K)
         ENDIF
  130    CONTINUE
  140    CONTINUE
      ENDIF
C     CHECK FOR SIMPLE ERRORS.
      KSTOP0 = KSTOP
      DO 160 J=1,NCIREF
      NELCHK = 0
      DO 150 I=1,NACTIV
      IF(JCIREF(I,J).LT.0 .OR. JCIREF(I,J).GT.2) THEN
         KSTOP = KSTOP+1
         WRITE(NB6,680) I,J,JCIREF(I,J)
      ELSE
         NELCHK = NELCHK+JCIREF(I,J)
      ENDIF
  150 CONTINUE
      IF(KSTOP0.EQ.0 .AND. NELCHK.NE.NELACT) THEN
         KSTOP = KSTOP+1
         WRITE(NB6,690) J,NELCHK,NELACT
      ENDIF
  160 CONTINUE
C
C *** SIXTH LINE: REQUESTED NUMBER OF CI ROOTS OF EACH IRREP.
      IROOTA = 0
      IROOT  = IN2(139)
      IF(IROOT.LT.0) THEN
         IF(IFORM.LE.0) THEN
            READ(NB5,500) IROOTA
         ELSE
            READ(NB5,*)   IROOTA
         ENDIF
      ENDIF
C
C *** SEVENTH LINE: INDICES OF CONJUGATED ATOMS.
      NCONJ = IN2(167)
      IF (MOVO .LT. 0 .AND. NCONJ .GT. 0) THEN
         IF(IFORM.LE.0) THEN
            READ(NB5,500) (ICONJ(I),I=1,NCONJ)
         ELSE
            READ(NB5,*)   (ICONJ(I),I=1,NCONJ)
         ENDIF
      ENDIF
C
C
C *** NEXT TWO LINES: INPUT FOR MULTI-SURFACE CALCULATIONS.
C     FIRST  LINE: CI STATES INVOLVED.
C     SECOND LINE: SPECIFIC OPTIONS AND CRITERIA.
C     INPUT DONE IN SEPARATE ROUTINES CALLED AFTER INGUGA.
C     ICROSS=1-2: CALL INMULS FOR SINGLE POINT MULTISTATE CALCULATION.
C     ICROSS=3-5: CALL INCONI FOR CONICAL INTERSECTIONS.
C     ICROSS=6: CALL INEXDY FOR EXCITED STATE DYNAMICS.
C
C
C *** CHECK WHETHER TO RESTORE PREVIOUS INPUT OPTIONS.
      IF(IN1(158).EQ.1) THEN
         IN2(131) = IN2SAV(131)
         IN2(132) = IN2SAV(132)
         IN2(134) = IN2SAV(134)
         IN2(136) = IN2SAV(136)
         IN2(137) = IN2SAV(137)
         IN2(138) = IN2SAV(138)
         IN2(141) = IN2SAV(141)
         IN2(142) = IN2SAV(142)
         IN2(143) = IN2SAV(143)
         IN2(144) = IN2SAV(144)
         IN2(151) = IN2SAV(151)
         IN2(152) = IN2SAV(152)
         IN2(153) = IN2SAV(153)
         IN2(155) = IN2SAV(155)
         IN2(168) = IN2SAV(168)
         IN2(169) = IN2SAV(169)
         IN2(170) = IN2SAV(170)
         WRITE(NB6,910)
      ENDIF
C *** PRINTING SECTION.
      IF(JPRINT.GT.-5) THEN
         WRITE(NB6,510)
         WRITE(NB6,520) IN2(131:170)
         WRITE(NB6,570)
         WRITE(NB6,580) (JMOCI(I),I=1,NACTIV)
         WRITE(NB6,585) (MOCISY(I),I=1,NACTIV)
      ENDIF
      IF(JPRINT.GE.0 .AND. NCIREF.GT.0) THEN
         WRITE(NB6,590)
         DO 170 J=1,NCIREF
         WRITE(NB6,560) J,(JCIREF(I,J),I=1,NACTIV)
  170    CONTINUE
      ENDIF
      IF(IROOT.LT.0) THEN
         WRITE(NB6,525) (I,I=1,8),IROOTA
      ENDIF
      IF(JPRINT.GE.0.AND.IN1(158).EQ.1) THEN
         WRITE(NB6,920) IN1(158)
         WRITE(NB6,570)
         WRITE(NB6,580) (IMOCI(I),I=1,IN2(131)+IN2(132))
         WRITE(NB6,590)
         DO 180 J=1,IN2(136)
         WRITE(NB6,560) J,(ICIREF(I,J),I=1,IN2(131)+IN2(132))
  180    CONTINUE
      ENDIF
C
C *** STOP IN CASE OF FATAL INPUT ERRORS.
      IF(KSTOP.GT.0) THEN
         WRITE(NB6,800) KSTOP
         STOP 'INGUGA'
      ENDIF
      RETURN
  500 FORMAT(20I4)
  505 FORMAT(///1X,'*** OPTIONS FOR THE GUGA-CI CALCULATIONS    ***',
     1       /  1X,'*** FROM INPUT - NO DEFAULTS APPLIED YET    ***')
  510 FORMAT(///1X,'*** OPTIONS FOR THE GUGA-CI CALCULATIONS    ***',
     1       /  1X,'*** DEFINED EITHER EXPLICITLY OR BY DEFAULT ***')
  520 FORMAT(/  1X,'ICI1   =',I5,5X,'ICI2   =',I5,5X,'IOUTCI =',I5,
     1          5X,'MOVO   =',I5,5X,'MPERT  =',I5,
     2       /  1X,'NCIREF =',I5,5X,'MCIREF =',I5,5X,'LEVEXC =',I5,
     3          5X,'IROOT  =',I5,5X,'LROOT  =',I5,
     4       /  1X,'CICHG  =',I5,5X,'MULTCI =',I5,5X,'NCISYM =',I5,
     5          5X,'CIDIR  =',I5,5X,'CIDIAG =',I5,
     6       /  1X,'IUVCD  =',I5,5X,'IMCD   =',I5,5X,'IPOP   =',I5,
     7          5X,'CIPLOT =',I5,5X,'CILEAD =',I5,
     8       /  1X,'JCI1   =',I5,5X,'JCI2   =',I5,5X,'PIPOP  =',I5,
     9          5X,'INATUR =',I5,5X,'CISELT =',I5,
     A       /  1X,'IMOMAP =',I5,5X,'GUG157 =',I5,5X,'KEEPCI =',I5,
     B          5X,'NCIGRD =',I5,5X,'ICROSS =',I5,
     C       /  1X,'MINDAV =',I5,5X,'MAXDAV =',I5,5X,'KITDAV =',I5,
     D          5X,'NRMDAV =',I5,5X,'MAXMAP =',I5,
     E       /  1X,'MAPTHR =',I5,5X,'NCONJ  =',I5,5X,'GUG168 =',I5,
     F          5X,'GUG169 =',I5,5X,'GUG170 =',I5)
  525 FORMAT(/  8X,8I4/1X,'IROOTA:',8I4/)
  550 FORMAT(/  1X,'INGUGA: ECHO INPUT FOR OPTION MCIREF=2.',
     1       /  1X,'EXCITATION INDICES TO DEFINE REFERENCES.'/)
  560 FORMAT(/  1X,'REF',I3,2X,20I4,/(9X,20I4))
C    1       /  1X,         8X,20I4,
C    2       /  1X,         8X,20I4)
  570 FORMAT(///1X,'ORBITALS INCLUDED IN THE ACTIVE CI SPACE.',
     1       /  1X,'THE NUMBERING REFERS TO THE SCF SOLUTION.'/)
  580 FORMAT(   1X,'MOS',5X,20I4,/(9X,20I4))
C    1       /  1X,      8X,20I4,
C    2       /  1X,      8X,20I4)
  585 FORMAT(/  1X,'MOSYM',3X,20I4,/(9X,20I4))
C    1       /  1X,        8X,20I4,
C    2       /  1X,        8X,20I4)
  590 FORMAT(///1X,'REFERENCE CONFIGURATIONS FOR GUGA-CI.',
     1       /  1X,'OCCUPATION NUMBERS FOR THE ACTIVE ORBITALS.')
  600 FORMAT(/  1X,'INGUGA: ICI1 RESET FROM',I4,' TO', I4)
  605 FORMAT(/  1X,'INGUGA: JCI1 RESET FROM',I4,' TO', I4)
  610 FORMAT(/  1X,'INGUGA: ICI2 RESET FROM',I4,' TO', I4)
  615 FORMAT(/  1X,'INGUGA: JCI2 RESET FROM',I4,' TO', I4)
  620 FORMAT(/  1X,'INGUGA: ERROR - TOO MANY ACTIVE ORBITALS.',
     1       /  1X,'        REQUESTED',I5,
     2       /  1X,'        MAXIMUM  ',I5)
  625 FORMAT(/  1X,'INGUGA: ERROR - INVALID OPTION MOVO (',I2,').')
  630 FORMAT(/  1X,'INGUGA: ERROR - TOO MANY REFERENCE CONFIGURATIONS.',
     1       /  1X,'        REQUESTED',I5,
     2       /  1X,'        MAXIMUM  ',I5)
  640 FORMAT(/  1X,'INGUGA: NCIREF=0 IMPLIES FULL CI, RESET MCIREF=0.')
  650 FORMAT(/  1X,'INGUGA: NCIREF=0 IMPLIES FULL CI, RESET LEVEXC=0.')
  651 FORMAT(/  1X,'INGUGA: FULL CI, RESET NCIREF = 1.')
  652 FORMAT(/  1X,'INGUGA: FULL CI, RESET LEVEXC =',I2,'.')
  660 FORMAT(/  1X,'INGUGA: AN AUTOMATIC GENERATION OF REFERENCE',
     1          1X,'CONFIGURATIONS IS NOT POSSIBLE.',
     2       /  1X,'        PLEASE CHECK THE FOLLOWING INPUT OPTIONS.',
     3       /  1X,'        NUMBER OF REFERENCES (NCIREF)',I5,
     4       /  1X,'        CHARGE FOR CI        (CICHG) ',I5,
     5       /  1X,'        MULTIPLICITY FOR CI  (MULTCI)',I5)
  670 FORMAT(/  1X,'INGUGA: EXCITATION INDICES OUT OF ALLOWED RANGE.',
     1       /  1X,'        REFERENCE CONFIGURATION ',I6,
     2       /  1X,'        EXCITATION NUMBER (PAIR)',I6,
     3       /  1X,'        EXCITATION INDICES      ',2I6)
  680 FORMAT(/  1X,'INGUGA: REFERENCE OCCUPANCY OUT OF RANGE.',
     1       /  1X,'        ORBITAL / CONFIGURATION ',2I6,
     2       /  1X,'        OCCUPATION NUMBER       ',I6)
  690 FORMAT(/  1X,'INGUGA: ERROR IN REFERENCE CONFIGURATION ',I6,
     1       /  1X,'        ACTUAL NUMBER OF ACTIVE ELECTRONS',I6,
     2       /  1X,'        TRUE   NUMBER OF ACTIVE ELECTRONS',I6)
  700 FORMAT(/  1X,'INGUGA: PLEASE SPECIFY MULTIPLICITY OF CI STATE')
  710 FORMAT(/  1X,'INGUGA: NCISYM AND NEGATIVE IROOT SPECIFIED',
     1       /  1X,'        SIMULTANEOUSLY.')
  720 FORMAT(/  1X,'INGUGA: CIDIAG.EQ.11 CAN CALCULATE ONLY ONE ROOT.')
  725 FORMAT(/  1X,'INGUGA: NEGATIVE PIPOP IS NO LONGER SUPPORTED.')
  730 FORMAT(/  1X,'INGUGA: ERROR - TOO MANY CI GRADIENTS (NCIGRD).',
     1       /  1X,'        REQUESTED',I5,
     2       /  1X,'        MAXIMUM  ',I5)
  740 FORMAT(/  1X,'INGUGA: ERROR - OPTION NOT AVAILABLE: ICROSS =',I3)
  750 FORMAT(/  1X,'INGUGA: ERROR - IMOMAP=2 ONLY AVAILABLE FOR JOP=0')
  760 FORMAT(/  1X,'INGUGA: ERROR - IMOMAP=2 ONLY AVAILABLE FOR IEF>=0')
  770 FORMAT(/  1X,'INGUGA: ERROR - IMOMAP=2 ONLY AVAILABLE FOR LSUB=0')
  775 FORMAT(/  1X,'INGUGA: ERROR - IMOMAP=3 ONLY AVAILABLE FOR JOP<0')
  780 FORMAT(/  1X,'INGUGA: ERROR - MULTCI=-1 REQUIRES NCIGRD=2')
  790 FORMAT(/  1X,'INGUGA: ERROR - MULTCI=-1 REQUIRES ICROSS=1-5')
  800 FORMAT(/  1X,'INGUGA: NUMBER OF FATAL INPUT ERRORS',I6/)
  900 FORMAT(/  1X,'INGUGA: PREVIOUS CI OPTIONS SAVED SINCE KEEPCI =',
     1             I3/1X,'******* EXPERIMENTAL OPTION,',
     2          1X,'USE AT OWN RISK. *******')
  910 FORMAT(/  1X,'INGUGA: PREVIOUS CI OPTIONS RESTORED.')
  920 FORMAT(///1X,'INGUGA: EFFECTIVE OPTIONS DUE TO KEEPCI =',I3)
      END
