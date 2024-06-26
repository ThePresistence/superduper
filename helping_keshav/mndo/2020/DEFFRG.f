      SUBROUTINE DEFFRG (IFRAG,COORDS,NUMATS,NATS,NFRAGS,NCHRGS,JPRINT)
C     *
C     EXTRACT FRAGMENT DATA FROM MOLECULAR DATA.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     IFRAG     NUMBER OF CURRENT FRAGMENT (I).
C     COORDS    CARTESIAN COORDINATES OF MOLECULE (I).
C     NUMATS    NUMBER OF ATOMS IN MOLECULE (I).
C     NATS      ATOMIC NUMBERS OF MOLECULE (I).
C     NFRAGS    FRAGMENT ASSIGNMENT OF ATOMS IN MOLECULE (I).
C     NCHRGS    FORMAL CHARGES OF ATOMS IN MOLECULE (I).
C     JPRINT    PRINTING OPTION (I).
C     *
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF
      CHARACTER*80 KTITLE,KOMENT
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DELEMT/ NELMD,NFOCK
     ./DPARM / LORBS(LMZ)
     ./ENERGT/ EE,ENUCLR,EAT,ATHEAT
     ./FLAG1 / KTITLE,KOMENT
      COMMON
     ./HALFE / IODD,JODD
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
     ./PARDER/ CORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./UHF   / UHF
      DIMENSION COORDS(3,NUMATS),NATS(NUMATS)
      DIMENSION NFRAGS(NUMATS),NCHRGS(NUMATS)
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** EXTRACT CARTESIAN COORDINATES AND ATOMIC NUMBERS.
      KHARGE = 0
      NUMAT  = 0
      DO 10 I=1,NUMATS
      IF(IFRAG.EQ.NFRAGS(I)) THEN
         KHARGE = KHARGE+NCHRGS(I)
         NUMAT  = NUMAT+1
         COORD(1,NUMAT) = COORDS(1,I)
         COORD(2,NUMAT) = COORDS(2,I)
         COORD(3,NUMAT) = COORDS(3,I)
         NAT(NUMAT)     = NATS(I)
      ENDIF
   10 CONTINUE
C *** DETERMINE THE FOLLOWING VARIABLES.
C     NFIRST(I) FIRST BASIS ORBITAL OF ATOM I.
C     NLAST(I)  LAST  BASIS ORBITAL OF ATOM I.
C     NEL       NUMBER OF ELECTRONS.
C     NELMD     NUMBER OF ATOMS WITH D ORBITALS.
C     NFOCK     STEP SIZE IN FOCK (5 SP, 10 SPD).
C     IMULT     MULTIPLICITY OF FRAGMENT.
      NELMD  = 0
      NFOCK  = 5
      NEL    = -KHARGE
      IB     = 0
      DO 20 I=1,NUMAT
      NFIRST(I) = IB+1
      NI     = NAT(I)
      NEL    = NEL+NINT(CORE(NI))
      IB     = IB+LORBS(NI)
      NLAST(I) = IB
      IF(LORBS(NI).GE.9) THEN
         NELMD = NELMD+1
         NFOCK = 10
      ENDIF
   20 CONTINUE
C     IMULT  = 0  FOR EVEN NUMBER OF ELECTRONS
C     IMULT  = 2  FOR ODD  NUMBER OF ELECTRONS
      IMULT  = 0
      IF(MOD(NEL,2).EQ.1) IMULT=2
C *** PRINT THE CARTESIAN COORDINATES.
      IF(JPRINT.GE.5) THEN
         WRITE(NB6,600) IFRAG,KOMENT,KTITLE
         WRITE(NB6,610)
         DO 30 I=1,NUMAT
         NI  = NAT(I)
         IA  = NFIRST(I)
         IB  = NLAST (I)
         WRITE(NB6,620) I,NI,(COORD(J,I),J=1,3),IA,IB
   30    CONTINUE
      ENDIF
C *** INTERNUCLEAR DISTANCES.
      KSTOP  = 0
      IF(JPRINT.GE.5) THEN
         WRITE(NB6,630)
         CALL RIJPRT (COORD,NAT,NUMAT,KSTOP)
      ELSE
         CALL RIJCHK (COORD,NAT,NUMAT,KSTOP)
      ENDIF
      IF(KSTOP.GT.0) THEN
         WRITE(NB6,640)
         STOP 'INPFRG'
      ENDIF
C *** SET UHF OPTION.
      UHF    = .FALSE.
C *** DETERMINE OCCUPATION OF THE MOLECULAR ORBITALS.
C     NEL      NUMBER OF ELECTRONS       (NEL = NALPHA+NBETA).
C     NALPHA   NUMBER OF ALPHA ELECTRONS (CORRECT FOR RHF AND UHF).
C     NBETA    NUMBER OF BETA  ELECTRONS (CORRECT FOR RHF AND UHF).
C     NBETA    EQUAL TO NUMBER OF DOUBLY OCCUPIED MOLECULAR ORBITALS
C              IN RHF CALCULATIONS WITH IMULT.NE.1.
C     NUMB     NUMBER OF HIGHEST OCCUPIED MOLECULAR ORBITAL.
      NORBS  = NLAST(NUMAT)
      IODD   = 0
      JODD   = 0
      NODD   = MAX(1,IMULT)-1
      NBETA  = (NEL-NODD)/2
      NALPHA = NBETA+NODD
      NUMB   = NALPHA
      NCLO   = NBETA
      NMOS   = NUMB
C     IMULT  = 0  ! FOR EVEN NUMBER OF ELECTRONS
C     IMULT  = 2  ! FOR ODD  NUMBER OF ELECTRONS
      IF(IMULT.EQ.2) THEN
        IODD = NUMB
      ENDIF
C *** COMPUTE SUMS OF ATOMIC ENERGIES AND HEATS OF FORMATION.
      ATHEAT = ZERO
      EAT    = ZERO
      DO 50 I=1,NUMAT
      NI     = NAT(I)
      ATHEAT = ATHEAT + EHEAT(NI)
      EAT    = EAT + EISOL(NI)
   50 CONTINUE
C *** PRINTING SECTION.
      IF(JPRINT.GE.5) THEN
         IUHF   = -1
         KITSCF = IN2(71)
         IF(IN2(67).GE.31 .AND. IN2(67).LE.39) KITSCF=IN2(67)-30
         IFAST  = 2
         IDIAG  = IN2(74)
C        IDIAG  = 2
         KCI    = 0
         NSTART =-1
         WRITE(NB6,700) KHARGE,NEL,NCLO
         IF(IODD .GT.0) WRITE(NB6,710) IODD
         IF(IMULT.GT.0) WRITE(NB6,720) IMULT
         IF(IMULT.GT.0) WRITE(NB6,730)
         WRITE(NB6,740)
         WRITE(NB6,750) KHARGE ,IMULT  ,IN2(67),IN2(68),IN2(69),
     1                  IUHF   ,KITSCF ,IN2(72),IFAST  ,IDIAG  ,
     2                  IN2(75),IN2(76),KCI    ,NSTART ,IN2(79)
         WRITE(NB6,760)
         WRITE(NB6,790)
      ENDIF
      RETURN
  600 FORMAT(///1X,'*** FRAGMENT',I5,//5X,A,/5X,A)
  610 FORMAT(// 5X,'INITIAL CARTESIAN COORDINATES (ANGSTROMS)',
     1       // 5X,'ATOM NO.',8X,'ATOMIC NO.',9X,'X-COORDINATE',
     2          8X,'Y-COORDINATE',8X,'Z-COORDINATE',
     3          7X,'ATOMIC ORBITALS'//)
  620 FORMAT(   7X,I3,16X,I2,14X,3(F10.6,10X),I3,' TO ',I3 )
  630 FORMAT(///5X,'INITIAL INTERATOMIC DISTANCES (ANGSTROMS)')
  640 FORMAT(// 1X,'THERE ARE AT LEAST TWO ATOMS WITH ZERO DISTANCE.',
     1       /  1X,'CHECK INPUT GEOMETRY. STOP.'//)
  700 FORMAT(///1X,'MOLECULAR CHARGE     ',I5,
     1       /  1X,'NUMBER OF ELECTRONS  ',I5,
     2       /  1X,'DOUBLY OCCUPIED MOS  ',I5)
  710 FORMAT(   1X,'ONE ELECTRON IN MO   ',I5)
  720 FORMAT(   1X,'MULTIPLICITY         ',I5)
  730 FORMAT(/  1X,'DOUBLET ASSUMED BY DEFAULT FOR ODD-ELECTRON SYSTEM')
  740 FORMAT(// 1X,'*** OPTIONS FOR THE MOLECULE UNDER STUDY    ***',
     1       /  1X,'*** DEFINED EITHER EXPLICITLY OR BY DEFAULT ***')
  750 FORMAT(/  1X,'KHARGE =',I5,5X,'IMULT  =',I5,5X,'KTRIAL =',I5,5X,
     1             'KGEOM  =',I5,5X,'IPUBO  =',I5,
     2       /  1X,'IUHF   =',I5,5X,'KITSCF =',I5,5X,'NPRINT =',I5,5X,
     3             'IFAST  =',I5,5X,'IDIAG  =',I5,
     4       /  1X,'KSYM   =',I5,5X,'NUMSYM =',I5,5X,'KCI    =',I5,5X,
     5             'NSTART =',I5,5X,'NSTEP  =',I5)
  760 FORMAT(/  1X,'RHF CALCULATION')
  790 FORMAT(/  1X,'FAST DIAGONALIZATIONS NOT ALLOWED IN SCF')
      END
