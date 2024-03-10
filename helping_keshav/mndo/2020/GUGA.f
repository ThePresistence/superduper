      SUBROUTINE GUGA (C,E,F,H,W,LM2,LM3,LM4,LM6,
     1                 NROOT,NOPRT,A,LM5,ICALL)
C     *
C     CONTROL ROUTINE FOR GUGA CI CALCULATION.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     C(LM2,LM3)     EIGENVECTORS (I).
C     CNO(LM2,LM3)   NATURAL ORBITAL MO COEFFICIENTS (O).
C     E(LM3)         EIGENVALUES (I).
C     F(LM4)         FOCK MATRIX (I).
C     H(LM4)         ONE-ELECTRON INTEGRALS (I).
C     W(LM6,LM6)     TWO-ELECTRON INTEGRALS (I).
C     NROOT          NUMBER OF AVAILABLE EIGENVECTORS (I).
C     NOPRT          PRINTING FLAG (I).
C     A(LM5)         AVAILABLE BUFFER (S).
C     ICALL          CONTROL AND ERROR FLAG (I,O).
C     SEE SUBROUTINE SCF FOR DEFINITION OF ICALL.
C     *
      USE LIMIT, ONLY: LM1, LMX, LMZ, LMACT, LMREF, LMPROP, LMSTAT,
     1                 LMGRD, MAXGPU
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL NOPRT,PRT
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CIPRP / CIPROP(LMPROP,LMSTAT),ICISYM(LMSTAT)
     ./CIREFS/ ICIREF(LMACT,LMREF),JCIREF(LMACT,LMREF)
     ./CIMOS / IMOCI(LMX)
     ./CIROOT/ IROOTA(8)
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./DENERG/ DENER
     ./DNBND / III(LMZ),IIID(LMZ)
     ./DPARM1/ UDD(LMZ),ZD(LMZ),BETAD(LMZ)
     ./ENERGT/ EE,ENUCLR,EAT,ATHEAT
     ./FLAG3 / KRESET,MPRINT,T2
     ./GPUPRP/ GPUNAM(MAXGPU),NVCAPA(MAXGPU),MIBGLO(MAXGPU),NUMGPU
     ./GRDORG/ IGRST(LMGRD),ISTATE,JSTATE
     ./GUGA1 / NCIO,NCIGAM
     ./GUGA2 / LCI1,LCI2,LCI3,LCI4
     ./GUGA3 / EEREF
     ./GUGA4 / XNAC(LMGRD,LMGRD)
     ./GUGA5 / MOLGUG,ICICAL,MULSAV
     ./IJWORK/ NSYM(LMX),ISYM(3,LM1)
     ./INOPT2/ IN2(300)
     ./MMCOM1/ EMM
     ./MMDP  / EMMDP
     ./CHCORR/ EHCORR
     ./NBFILE/ NBF(20)
     ./OCCFL / IMOCC,NOCCA,NOCCB,MSUB,MOSUMA,MOSUMB,MOCCA(8),MOCCB(8)
     ./OCCFL2/ DOMEGA,EFERMI,NFLOAT,NDOCC,NUMOCC
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
     ./PAR2  / MOL
     ./PARDER/ CORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ*3)
     ./SKIPA / ISKPA
      DIMENSION C(LM2,LM3),E(LM3),F(LM4),H(LM4),W(LM6,LM6)
      DIMENSION A(LM5),OCCR(LMX),NSYACT(LMX),OCCACT(LMX)
      CHARACTER( 4) LABELS(8)
      CHARACTER( 3) GRPLAB
      CHARACTER(80) GPUNAM
C     FLAG FOR CHANGE OF SIZE OF ACTIVE SPACE.
      LOGICAL LRESAS
C     INITIAL POSITION OF EACH ACTIVE ORBITAL.
      DIMENSION IMOPOS(LMACT)
C     SAVE INFORMATION ABOUT MOLECULAR SYMMETRY BECAUSE
C     MOSYM WILL BE CALLED ONLY ONCE FOR EACH GEOMETRY.
      SAVE NUMLAB, LABELS, GRPLAB, ISUB, ISUBSV, NSYACT, OCCACT
C     INITIALIZE INFORMATION ABOUT THE PREVIOUS GEOMETRY.
      DATA ISUBSV/0/,NSYACT/LMX*0/,OCCACT/LMX*0.0D0/
C *** A NON-ZERO MOLGUG INDICATES THAT THE DRT ETC. HAS BEEN
C     INITIALIZED FOR MOLECULE MOLGUG. AFTER CLEANING UP
C     MOLGUG IS RESET TO ZERO.
      DATA MOLGUG/0/
C *** COUNT CI FULL CALCULATIONS FOR THE CURRENT MOLECULE.
      DATA ICICAL/0/
C *** CHECK IF MULTIPLICITY CHANGES BETWEEN CALCULATIONS.
      DATA MULSAV/-1/
C *** FILE NUMBERS
      NB6    = NBF(6)
C *** INPUT OPTIONS.
      IATERG = IN2(119)
      ICI1   = IN2(131)
      ICI2   = IN2(132)
      IOUTCI = IN2(133)
      NCIREF = IN2(136)
      IROOT  = IN2(139)
      MULTCI = IN2(142)
      NCISYM = IN2(143)
      IUVCD  = IN2(146)
      IMCD   = IN2(147)
      IPOP   = IN2(148)
      INATUR = IN2(154)
      IMOMAP = IN2(156)
      KEEPCI = IN2(158)
      ICROSS = IN2(160)
      INOUT  = IN2(212)
      NACTIV = ICI1+ICI2
C *** CHECK INOUT.
      IF(INOUT.NE.0) THEN
         WRITE(NB6,600)
         STOP 'GUGA'
      ENDIF
C *** INITIALIZATION.
      IF(NOPRT .AND. IOUTCI.LE.0) THEN
         PRT      = .FALSE.
         IN2(133) = -11
      ELSE
         PRT      = .TRUE.
      ENDIF
      IF(ICALL.GT.10) THEN
         IUVCD  = -1
         IMCD   = -1
         IPOP   = -1
         INATUR = -1
      ENDIF
      LRESAS = .FALSE.
      IF(MULSAV.EQ.-1) MULSAV = MULTCI
      CALL GUGDAT(C,E,W,LM2,LM3,LM6,NROOT,OCCR,IMOPOS,
     1            NUMLAB,LABELS,GRPLAB,ISUB,MOLGUG,MOL,ICALL,NOPRT)
      IF(ICALL.EQ.-1) RETURN
C     NUMBER OF ACTIVE ORBITALS MAY HAVE BEEN CHANGED BY GUGDAT
      IF ((ICI1 .NE. IN2(131)) .OR. (ICI2 .NE. IN2(132))) THEN
         ICI1   = IN2(131)
         ICI2   = IN2(132)
         NACTIV = ICI1+ICI2
         LRESAS = .TRUE.
      ENDIF
C *** DEFINITION OF CALL NUMBER.
C     1: INITIALIZATION AND CALCULATION OF ALL COUPLING COEFFICIENTS
C        FOR CIDIR.EQ.1 OR PARTIAL LOOP VALUES FOR CIDIR.EQ.3.
C     2: CI CALCULATION FOR EACH SYMMETRY (INCLUDING COMPUTATION OF
C        COUPLING COEFFICIENTS FOR CIDIR.EQ.2) AND CHECKING OF
C        REFERENCE CONFIGURATIONS.
C     3: CALCULATE ONE- AND TWO-PARTICLE (TRANSITION) DENSITY
C        MATRICES NEEDED FOR THE ANALYTIC GRADIENT AND NON-ADIABATIC
C        COUPLING VECTORS (FULLY DIRECT ALGORITHM FOR CIDIR.EQ.2).
C     4: GENERATE NATURAL ORBITALS, PERFORM POPULATION ANALYSIS, AND
C        CALCULATE SPECTROSCOPIC PROPERTIES (FULLY DIRECT ALGORITHM
C        FOR CIDIR.EQ.2).
C     5: NUMERIC CALCULATION OF NON-ADIABATIC COUPLING ELEMENTS.
C     6: SAVE CI COEFFICIENTS.
C     7: DEALLOCATE ALL PREVIOUSLY ALLOCATED DYNAMIC MEMORY.
C
C     REMARKS:
C     DURING A GEOMETRY OPTIMIZATION OR FORCE CONSTANT ANALYSIS
C     A SINGLE-POINT ENERGY CALCULATION MAY BE FOLLOWED BY A REQUEST
C     TO COMPUTE THE GRADIENT. IN THIS CASE PSDRV IS CALLED DIRECTLY
C     WITHOUT ANOTHER INVOCATION OF SUBROUTINE GUGA. THEREFORE THE
C     ONE- AND TWO-PARTICLE DENSITY MATRICES FOR THE GRADIENT MUST
C     ALWAYS BE COMPUTED DURING A GEOMETRY OPTIMIZATION OR FORCE
C     CONSTANT ANALYSIS, REGARDLESS OF THE VALUE OF ICALL2.
C
C     COUPLING COEFFICIENTS OR PARTIAL LOOP VALUES MUST ONLY BE REUSED
C     IF POINT GROUP IS UNCHANGED OR SYMMETRY IS TURNED OFF (NCISYM.LT.0).
C     OTHERWISE THESE QUANTITIES WILL BECOME INVALID IF SYMMETRY IS
C     LOWERED OR ORBITALS OF DIFFERENT SYMMETRY HAVE SWITCHED. THIS
C     CONDITION IS DETECTED AND THE CORRESPONDING QUANTITIES WILL BE
C     DELETED AND RECOMPUTED.
C
C     ENERGY DERIVATIVES BY SLOW FINITE DIFFERENCES (ISKPA.EQ.1)
C     MAY NOT WORK IF THE STATE OF INTEREST IS SELECTED BY SYMMETRY
C     AND SYMMETRY WILL BE LOWERED BY THE DISPLACEMENTS.
C     IF THIS IS ATTEMPTED THE PROGRAM WILL ISSUE A WARNING.
C
      IF(IOUTCI.GE.5) THEN
         WRITE(NB6,530) ICALL, MOL, MOLGUG, ICICAL, MULSAV, MULTCI
      ENDIF
      ICALL1 = ICALL/10
      ICALL2 = ICALL-10*ICALL1
      IF(ISKPA.EQ.1 .AND.
     1   (NUMLAB.GT.1 .AND. (IROOT.LT.0 .OR. NCISYM.GT.0)) .AND.
     2   (ICALL1.EQ.2 .OR. ICALL1.EQ.3 .OR. ICALL2.GE.1)) THEN
C        ERROR HAS BEEN TURNED INTO A WARNING.
         WRITE(NB6,610)
C        STOP 'GUGA'
      ENDIF
C     IF THE SAVED CI INFORMATION IS OUT-DATED CLEAR IT.
   10 IF(MOLGUG.GT.0 .AND. (
     1   (MOLGUG.NE.MOL  .AND. KEEPCI.NE.1) .OR.
     2   (NCISYM.GE.0 .AND. ISUB.NE.ISUBSV) .OR.
     3   (NCISYM.GE.0 .AND.
     4    ANY(NSYM(IMOCI(1:NACTIV)).NE.NSYACT(1:NACTIV))) .OR.
     5   (ANY(OCCR(IMOCI(1:NACTIV)).NE.OCCACT(1:NACTIV))) .OR.
     6   (NCIREF.NE.IN2(136)) .OR.
     7   (MULTCI.NE.MULSAV) .OR.
     8   (IMOCC.GE.5) .OR.
     9   (LRESAS))) THEN
         MOLGUG = 0
         ICICAL = 0
         MULSAV = MULTCI
         ICALLG = 7
         CALL GUGACI (ICALL,ICALLG,LM2,LM3,LM4,LM6,NUMAT,NUMB,NCIO,
     1                NCIGAM,IUVCD,IMCD,IPOP,INATUR,ISTATE,JSTATE,
     2                NUMLAB,NUMGPU,NB6,IN2,IMOCI,ICIREF,IROOTA,EE,
     3                ENUCLR,C,F,W,OCCR,NSYM,COORD,NAT,NFIRST,NLAST,
     4                CORE,ZS,ZP,ZD,III,IIID,ICISYM,A(LCI1),A(LCI2),
     5                A(LCI3),A(LCI4),CIPROP,XNAC,LABELS,GRPLAB,NVCAPA)
      ENDIF
C     SAVE POINT GROUP AND ACTIVE MO INFORMATION.
      ISUBSV           = ISUB
      NSYACT(1:NACTIV) = NSYM(IMOCI(1:NACTIV))
      OCCACT(1:NACTIV) = OCCR(IMOCI(1:NACTIV))
C     RESTORE NUMBER OF REFERENCES.
      NCIREF = IN2(136)
C     IF THE DRT ETC. IS NOT INITIALIZED, DO SO.
      IF(MOLGUG.EQ.0) THEN
         MOLGUG = MOL
         ICALLG = 1
         CALL GUGACI (ICALL,ICALLG,LM2,LM3,LM4,LM6,NUMAT,NUMB,NCIO,
     1                NCIGAM,IUVCD,IMCD,IPOP,INATUR,ISTATE,JSTATE,
     2                NUMLAB,NUMGPU,NB6,IN2,IMOCI,ICIREF,IROOTA,EE,
     3                ENUCLR,C,F,W,OCCR,NSYM,COORD,NAT,NFIRST,NLAST,
     4                CORE,ZS,ZP,ZD,III,IIID,ICISYM,A(LCI1),A(LCI2),
     5                A(LCI3),A(LCI4),CIPROP,XNAC,LABELS,GRPLAB,NVCAPA)
      ENDIF
C     DO A CI IF A FULL CALCULATION IS REQUESTED OR
C     NO CI HAS BEEN DONE SINCE THE LAST CLEANUP.
      IF(ICALL2.LE.1 .OR. ICICAL.EQ.0 .OR. KRESET.EQ.1) THEN
         EE     = EEREF
         ICALLG = 2
         CALL GUGACI (ICALL,ICALLG,LM2,LM3,LM4,LM6,NUMAT,NUMB,NCIO,
     1                NCIGAM,IUVCD,IMCD,IPOP,INATUR,ISTATE,JSTATE,
     2                NUMLAB,NUMGPU,NB6,IN2,IMOCI,ICIREF,IROOTA,EE,
     3                ENUCLR,C,F,W,OCCR,NSYM,COORD,NAT,NFIRST,NLAST,
     4                CORE,ZS,ZP,ZD,III,IIID,ICISYM,A(LCI1),A(LCI2),
     5                A(LCI3),A(LCI4),CIPROP,XNAC,LABELS,GRPLAB,NVCAPA)
C        CALCULATE CI HEATS OF FORMATION FOR ALL STATES STORED IN CIPROP.
         IF(IATERG.EQ.1) THEN
            CIPROP(1,:) =
     &      (CIPROP(1,:)+EEREF+ENUCLR-EAT)*EVCAL+ATHEAT+EMM+EMMDP+EHCORR
         ELSEIF(IATERG.EQ.-1) THEN
            CIPROP(1,:) =
     &      (CIPROP(1,:)+EEREF+ENUCLR)*EVCAL+EMM+EMMDP+EHCORR
         ENDIF
C        RESET CISELT SO THAT NO MORE REFERENCE CONFIGURATONS WILL BE ADDED.
         IN2(155) = 0
C        COUNT FULL CALCULATIONS EXCEPT FOR A NUMERICAL GRADIENT.
         IF(KRESET.EQ.0) ICICAL = ICICAL + 1
C        IF NUMBER OF REFERENCES HAS CHANGED SAVE NEW REFERENCES AND REDO CI.
         IF(NCIREF.NE.IN2(136)) THEN
            DO IREF=NCIREF+1,IN2(136)
               DO IACT=1,NACTIV
                  JCIREF(IMOPOS(IACT),IREF) = ICIREF(IACT,IREF)
               ENDDO
            ENDDO
            GOTO 10
         ENDIF
      ENDIF
C     DURING A GEOMETRY OPTIMIZATION OR FORCE CONSTANT ANALYSIS, THE
C     ANALYTIC GRADIENT MIGHT BE REQUESTED WITHOUT ANOTHER CALL TO GUGA,
C     SO ALWAYS PREPARE THE CORRESPONDING DENSITIES IN THESE CASES.
      IF(ICALL1.GE.2 .OR. ICALL2.GE.1) THEN
         ICALLG = 3
         CALL GUGACI (ICALL,ICALLG,LM2,LM3,LM4,LM6,NUMAT,NUMB,NCIO,
     1                NCIGAM,IUVCD,IMCD,IPOP,INATUR,ISTATE,JSTATE,
     2                NUMLAB,NUMGPU,NB6,IN2,IMOCI,ICIREF,IROOTA,EE,
     3                ENUCLR,C,F,W,OCCR,NSYM,COORD,NAT,NFIRST,NLAST,
     4                CORE,ZS,ZP,ZD,III,IIID,ICISYM,A(LCI1),A(LCI2),
     5                A(LCI3),A(LCI4),CIPROP,XNAC,LABELS,GRPLAB,NVCAPA)
      ENDIF
C     NATURAL ORBITALS, POPULATION ANALYSIS, SPECTROSCOPIC PROPERTIES.
C     THIS CALL WILL ALSO COMPUTE THE AO-BASIS DENSITY MATRIX OF THE
C     CURRENT STATE NEEDED IN QM/MM CALCULATIONS.
      IF(ICALL1.LE.2 .AND. ICALL2.LE.1 .AND. KRESET.EQ.0) THEN
         ICALLG = 4
         CALL GUGACI (ICALL,ICALLG,LM2,LM3,LM4,LM6,NUMAT,NUMB,NCIO,
     1                NCIGAM,IUVCD,IMCD,IPOP,INATUR,ISTATE,JSTATE,
     2                NUMLAB,NUMGPU,NB6,IN2,IMOCI,ICIREF,IROOTA,EE,
     3                ENUCLR,C,F,W,OCCR,NSYM,COORD,NAT,NFIRST,NLAST,
     4                CORE,ZS,ZP,ZD,III,IIID,ICISYM,A(LCI1),A(LCI2),
     5                A(LCI3),A(LCI4),CIPROP,XNAC,LABELS,GRPLAB,NVCAPA)
      ENDIF
C     NUMERIC NON-ADIABATIC COUPLING ELEMENT.
      IF(MOD(ICROSS,10).EQ.6 .AND. KRESET.EQ.0 .AND. ICALL2.LE.1
     1   .AND. ICICAL.GE.2) THEN
         ICALLG = 5
         CALL GUGACI (ICALL,ICALLG,LM2,LM3,LM4,LM6,NUMAT,NUMB,NCIO,
     1                NCIGAM,IUVCD,IMCD,IPOP,INATUR,ISTATE,JSTATE,
     2                NUMLAB,NUMGPU,NB6,IN2,IMOCI,ICIREF,IROOTA,EE,
     3                ENUCLR,C,F,W,OCCR,NSYM,COORD,NAT,NFIRST,NLAST,
     4                CORE,ZS,ZP,ZD,III,IIID,ICISYM,A(LCI1),A(LCI2),
     5                A(LCI3),A(LCI4),CIPROP,XNAC,LABELS,GRPLAB,NVCAPA)
      ENDIF
C     SAVE CI COEFFICIENTS.
      IF((MOD(ICROSS,10).EQ.6 .OR. ICROSS.EQ.1 .OR. ICROSS.EQ.2 .OR.
     1    ICROSS.EQ.7) .AND. KRESET.EQ.0 .AND. ICALL2.LE.1) THEN
         ICALLG = 6
         CALL GUGACI (ICALL,ICALLG,LM2,LM3,LM4,LM6,NUMAT,NUMB,NCIO,
     1                NCIGAM,IUVCD,IMCD,IPOP,INATUR,ISTATE,JSTATE,
     2                NUMLAB,NUMGPU,NB6,IN2,IMOCI,ICIREF,IROOTA,EE,
     3                ENUCLR,C,F,W,OCCR,NSYM,COORD,NAT,NFIRST,NLAST,
     4                CORE,ZS,ZP,ZD,III,IIID,ICISYM,A(LCI1),A(LCI2),
     5                A(LCI3),A(LCI4),CIPROP,XNAC,LABELS,GRPLAB,NVCAPA)
      ENDIF
C *** RESTORE OUTPUT LEVEL.
      IN2(133) = IOUTCI
C *** PRINT CI HEAT OF FORMATION.
      IF(PRT .AND. ICALL.GE.0) THEN
         IF(IATERG.EQ.1) THEN
            DENER = ATHEAT+EVCAL*(EE+ENUCLR-EAT)+EMM+EMMDP+EHCORR
            WRITE(NB6,520) DENER
         ELSEIF(IATERG.EQ.-1) THEN
            DENER = EVCAL*(EE+ENUCLR)+EMM+EMMDP+EHCORR
            WRITE(NB6,525) DENER
         ENDIF
      ENDIF
      RETURN
  520 FORMAT(/  1X,'CI HEAT OF FORMATION ',F15.5,' KCAL/MOL'/)
  525 FORMAT(/  1X,'CI TOTAL ENERGY      ',F15.5,' KCAL/MOL'/)
  530 FORMAT(// 1X,'DEBUG PRINT OF FLAGS IN SUBROUTINE GUGA.',
     1       /  1X,'ICALL: ',I4/1X,'MOL:   ',I4/1X,'MOLGUG:',I4,
     2       /  1X,'ICICAL:',I4/1X,'MULSAV:',I4/1X,'MULTCI:',I4)
  600 FORMAT(   1X,'ERROR: INOUT MUST BE ZERO FOR GUGA-CI.')
  610 FORMAT(/  1X,'WARNING: ENERGY DERIVATIVES BY SLOW FINITE ',
     1             'DIFFERENCES MAY NOT WORK'/1X,'IF STATE OF ',
     2             'INTEREST IS SELECTED BY SYMMETRY.')
      END



      SUBROUTINE GUGDAT(C,E,W,LM2,LM3,LM6,NROOT,OCCR,IMOPOS,
     1                  NUMLAB,LABELS,GRPLAB,ISUB,
     2                  MOLGUG,MOL,ICALL,NOPRT)
C     *
C     PREPARE INPUT DATA FOR GUGA-CI CALCULATION.
C     *
      USE LIMIT, ONLY: LM1, LMX, LMACT, LMREF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL NOPRT
      COMMON
     ./CIREFS/ ICIREF(LMACT,LMREF),JCIREF(LMACT,LMREF)
     ./CIMAP / IFMAP,NEWMAP
     ./CIMOS / IMOCI(LMX)
     ./CIMOSY/ JMOCI(LMX),MOCISY(LMX)
     ./FLAG3 / KRESET,MPRINT,T2
     ./IJWORK/ NSYM(LMX),ISYM(3,LM1)
     ./HALFE / IODD(2)
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./OCCFL / IMOCC,NOCCA,NOCCB,MSUB,MOSUMA,MOSUMB,MOCCA(8),MOCCB(8)
     ./OCCFL2/ DOMEGA,EFERMI,NFLOAT,NDOCC,NUMOCC
     ./OCCNM / OCCA(LMX),OCCB(LMX)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
     ./SYMLAB/ NGROUP(7),IRREP(26)
     ./SYMLIM/ IRREPA(7),IRREPN(7)

      DIMENSION C(LM2,LM3),E(LM3),W(LM6,LM6)
C     OCCUPATION NUMBERS OF ALL ORBITALS, 0.LE.OCCR(I).LE.2;
C     FOR ORBITALS IN THE CI ACTIVE SPACE, THESE ARE THE REFERENCE
C     OCCUPATION NUMBERS INCLUDED IN THE GUGA-CI SEGMENT VALUES.
      DIMENSION OCCR(LMX)
C     INITIAL POSITION OF EACH ACTIVE ORBITAL.
      DIMENSION IMOPOS(LMACT)
C     NUMBER OF MOS WITH SYMMETRY I.
      DIMENSION NMOSYM(8)
C     NUMBERING ORBITALS OF SAME SYMMETRY.
      DIMENSION IMOSYM(LMX)
C     ORBITAL INDEX BY IMOSYM AND SYMMETRY.
      DIMENSION KMOCI(LMX,8)
C     CLASSIFYING MOS BY THEIR OCCUPANCIES IN THE REFERENCES.
      DIMENSION MOTYPE(LMX), MOTYSV(LMX)
      SAVE MOTYPE
C     KMOMAP ASSIGNS A CURRENT ORBITAL TO EACH SAVED ORBITAL.
      DIMENSION KMOMAP(LMX)
C     LMOMAP INDICATES WHICH CURRENT ORBITALS HAVE BEEN MAPPED.
      DIMENSION LMOMAP(LMX)
      LOGICAL LMOMAP
      DIMENSION LLMOMAP(LMX)
      LOGICAL LLMOMAP
C     LISTS OF UNMATCHED CURRENT AND SAVED ORBITALS.
      DIMENSION IFAILC(LMX), IFAILS(LMX)
C     SCRATCH ARRAY.
      DIMENSION ITMP(LMREF)
C     FOR IMOMAP=1/2, IF THE MO COEFFICIENTS HAVE BEEN SAVED
C     IN CMOSAV FOR THE FOLLOWING RUN, LMOSAV IS SET TO .TRUE..
C     FOR IMOMAP=3, IF THE MO COEFFICIENTS ARE AVAILABLE ON FILE,
C     THEY WILL BE READ INTO CMOSAV AND LMOSAV WILL BE SET TO .TRUE..
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CMOSAV
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DOTPRO
      LOGICAL LMOSAV
      SAVE CMOSAV, LMOSAV
      DATA LMOSAV/.FALSE./
      CHARACTER(4) IRREP,LABELS(8),RLABEL(8)
      CHARACTER(3) NGROUP,GRPLAB
      CHARACTER(1) STAR
C     ACTIVE ORBITAL MAPPING STATUS
      DIMENSION ADTPRO(LMX)
      CHARACTER(9) CMAPST(3)
      DATA CMAPST/'UNCHANGED','CHANGED','FAILED'/
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** INPUT OPTIONS.
      NPRINT = IN2(72)
      ICI1   = IN2(131)
      ICI2   = IN2(132)
      IOUTCI = IN2(133)
      MOVO   = IN2(134)
      NCIREF = IN2(136)
      MCIREF = IN2(137)
      NCISYM = IN2(143)
      IMOMAP = IN2(156)
      KEEPCI = IN2(158)
      MAPTHR = IN2(166)
      NACTIV = ICI1+ICI2
C *** INITIALIZATION.
      IF(NOPRT .AND. IOUTCI.LE.0) THEN
         IOUTMO   = -11
         IN2(133) = -11
      ELSE
         IOUTMO   = IOUTCI
      ENDIF
      IOUTPI = IOUTMO
      ICALL1 = ICALL/10
      ICALL2 = ICALL-10*ICALL1
      PMAPTH = DBLE(MAPTHR)/1.0D2
C *** ORBITAL OCCUPATION NUMBERS.
      IF(IMOCC.LT.5) THEN
         OCCR(1:NUMB)      = 2.D0
         OCCR(NUMB+1:NMOS) = 0.D0
         IF(IODD(1).NE.0) THEN
            OCCR(IODD(1))  = 1.D0
            IF(IODD(2).NE.0) OCCR(IODD(2)) = 1.D0
         ENDIF
      ELSE
         OCCR(1:NUMOCC)    = OCCA(1:NUMOCC) * 2.0D0
      ENDIF
C *** IF ONLY THE GRADIENT IS REQUESTED, ALL DATA HAS BEEN PREPARED
C     BEFORE THE LAST FULL CALCULATION, SO SIMPLY RETURN.
      IF(ICALL2.GE.2 .AND. KRESET.EQ.0) RETURN
C *** ASSIGNMENT OF MO SYMMETRY.
      IF(NPRINT.GE.IOUTCI) IOUTMO=-11
      CALL MOSYM (C,E,LM2,LM3,NROOT,NSYM,ISYM,ISUB,IOUTMO)
      IF(ISUB.GT.0 .AND. NCISYM.GE.0) THEN
         NUMLAB           = IRREPN(ISUB)
         ILAB             = IRREPA(ISUB)+1
         JLAB             = IRREPA(ISUB)+NUMLAB
         LABELS(1:NUMLAB) = IRREP(ILAB:JLAB)
         GRPLAB           = NGROUP(ISUB)
      ELSE
         NUMLAB           = 1
         LABELS(1)        = 'A   '
         GRPLAB           = 'C1 '
         NSYM             = 1
      END IF
      IF(NUMLAB.LT.8) LABELS(NUMLAB+1:8) = 'XXX '
C *** COUNT MOS WITH THE SAME SYMMETRY.
      NMOSYM = 0
      KMOCI  = 0
      DO I=1,NMOS
         NMOSYM(NSYM(I))           = NMOSYM(NSYM(I)) + 1
         IMOSYM(I)                 = NMOSYM(NSYM(I))
         KMOCI(IMOSYM(I), NSYM(I)) = I
      ENDDO
C *** FOR IMOMAP=3, READ SAVED MO DATA FROM FILE IF AVAILABLE
      IF (IMOMAP.EQ.3) THEN
         IF (NEWMAP) THEN
            LMOSAV=0
            IF (ALLOCATED(CMOSAV)) DEALLOCATE(CMOSAV)
         END IF
         IF (LMOSAV .OR. ALLOCATED(CMOSAV)) THEN
            WRITE(NB6,700)
            STOP 'GUGDAT'
         ENDIF
         ALLOCATE(CMOSAV(NORBS, NMOS))
         CALL MOMSAV(1, NORBS, NMOS, NACTIV, NCIREF,
     1               ICI1, ICI2, MOVO, MCIREF, NCISYM,
     2               MOTYPE, CMOSAV, IMOCI, ICIREF, IMOPOS,
     3               LMOSAV)
         IF (.NOT. LMOSAV) THEN
            WRITE(NB6,710)
            DEALLOCATE(CMOSAV)
         ENDIF
      ENDIF
C
C *** SKIP MANIPULATION OF ACTIVE ORBITAL NUMBERS AND REFERENCE
C     CONFIGURATIONS IF PRIOR INFORMATION IS TO BE KEPT.
C
      IF((MOLGUG.GT.0 .AND. KEEPCI.EQ.1) .OR.
     1   (MOLGUG.GT.0 .AND. MOLGUG.EQ.MOL .AND.
     2    IMOMAP.GE.1 .AND. LMOSAV) .OR.
     3   (IMOMAP.EQ.3 .AND. LMOSAV)) GOTO 10
C
C *** DEFINITION OF MOS TO BE INCLUDED IN CORRELATION TREATMENT.
      IF(MOVO.LT.0) THEN
         CALL PIMOCI (C,E,LM2,LM3,IOUTPI,.FALSE.)
      ELSE
         DO I=1,NACTIV
            IF(MOCISY(I).EQ.0) THEN
               IMOCI(I) = JMOCI(I)
            ELSE
               IMOCI(I) = KMOCI(JMOCI(I), MOCISY(I))
               IF(IMOCI(I) == 0) THEN
                  WRITE(NB6,610) LABELS(MOCISY(I)), MOCISY(I)
                  STOP 'GUGDAT'
               ENDIF
            ENDIF
         ENDDO
      ENDIF
C *** WHEN REFERENCE CONFIGURATIONS ARE AUTOMATICALLY GENERATED,
C     AN ASCENDING ORDER OF ACTIVE MO INDICES IS EXPECTED.
C     THUS IF ACTIVE ORBITALS HAVE BEEN SPECIFIED BY THE USER
C     AND THE REFERENCE CONFIGURATIONS WERE DETERMINED BY DEFAULT,
C     SORT ACTIVE MOS ACCORDING TO INCREASING INDEX.
C
      IF(MOVO.GT.0 .AND. (MCIREF.EQ.0 .OR. MCIREF.EQ.3)) THEN
         IF(IOUTCI.GE.5) THEN
            WRITE(NB6,'(/1X,A)')'UNSORTED ACTIVE MOS:'
            WRITE(NB6,'( 1X,30I4)')  IMOCI(1:NACTIV)
         ENDIF
         DO I=1,NACTIV-1
            K = I
            DO J=I+1,NACTIV
               IF(IMOCI(J).LT.IMOCI(K))  K = J
            ENDDO
            IF(K.NE.I) THEN
               IF(IOUTCI.GE.5)
     1           WRITE(NB6,'(1X,A,I3,A,I3,A)')
     2           'EXCHANGING ACTIVE MOS ',I,' AND ',K,'.'
               IHELP       = IMOCI(I)
               IMOCI(I)    = IMOCI(K)
               IMOCI(K)    = IHELP
            ENDIF
         ENDDO
         IF(IOUTCI.GE.5) THEN
            WRITE(NB6,'(/1X,A)')'SORTED ACTIVE MOS:'
            WRITE(NB6,'( 1X,30I4)')  IMOCI(1:NACTIV)
         ENDIF
      ENDIF
C *** INITIALIZE REFERENCE CONFIGURATIONS.
      ICIREF = JCIREF
C *** INITIALIZE ACTIVE MO POSITIONS.
      DO I=1,NACTIV
         IMOPOS(I) = I
      ENDDO
C *** DETERMINE ACTIVE MO TYPES.
C      -1: ORBITAL IS NOT IN THE ACTIVE SPACE.
C       0: ORBITAL IS UNOCCUPIED IN ALL REFERENCE CONFIGURATIONS.
C       1: ORBITAL HAS VARYING OCCUPANCIES IN THE REFERENCE CONFIGURATIONS.
C       2: ORBITAL IS DOUBLY OCCUPIED IN ALL REFERENCE CONFIGURATIONS.
C
      MOTYPE(1:NMOS) = -1
      DO I=1,NACTIV
         MOTYPE(IMOCI(I)) = ICIREF(I,1)
         IF(MOTYPE(IMOCI(I)).EQ.0 .OR. MOTYPE(IMOCI(I)).EQ.2) THEN
            DO J=2, NCIREF
               IF(ICIREF(I,J).NE.MOTYPE(IMOCI(I))) THEN
                  MOTYPE(IMOCI(I)) = 1
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      IF(IOUTCI.GE.5) THEN
         WRITE(NB6,'(/1X,A)')    'UNSORTED MO TYPES:'
         WRITE(NB6,'( 1X,30I4)') IMOCI(1:NACTIV)
         WRITE(NB6,'( 1X,30I4)') MOTYPE(IMOCI(1:NACTIV))
         WRITE(NB6,*)
      ENDIF
C *** SORT ACTIVE MOS ACCORDING TO DECREASING TYPE AND INCREASING NUMBER
C     SO THAT THE DRT WILL BECOME AS COMPACT AS POSSIBLE.
C
      DO I=1,NACTIV-1
         K = I
         DO J=I+1,NACTIV
            IF( MOTYPE(IMOCI(J)).GT.MOTYPE(IMOCI(K)) .OR.
     1         (MOTYPE(IMOCI(J)).EQ.MOTYPE(IMOCI(K)) .AND.
     2          IMOCI(J).LT.IMOCI(K)) )
     3         K = J
         ENDDO
         IF(K.NE.I) THEN
            IF(IOUTCI.GE.5)
     1         WRITE(NB6,'(1X,A,I3,A,I3,A)')
     2           'EXCHANGING ACTIVE MOS ',I,' AND ',K,'.'
            IHELP       = IMOCI(I)
            IMOCI(I)    = IMOCI(K)
            IMOCI(K)    = IHELP
            IHELP       = IMOPOS(I)
            IMOPOS(I)   = IMOPOS(K)
            IMOPOS(K)   = IHELP
            ITMP        = ICIREF(I,:)
            ICIREF(I,:) = ICIREF(K,:)
            ICIREF(K,:) = ITMP
         ENDIF
      ENDDO
      IF(IOUTCI.GE.5) THEN
         WRITE(NB6,'(/1X,A)')    'SORTED MO TYPES:'
         WRITE(NB6,'( 1X,30I4)') IMOCI(1:NACTIV)
         WRITE(NB6,'( 1X,30I4)') MOTYPE(IMOCI(1:NACTIV))
         WRITE(NB6,'( 1X,30I4)') IMOPOS(1:NACTIV)
      ENDIF
C *** MAP CURRENT MOS ON PREVIOUS ONES IF REQUESTED BY INPUT.
C Change by Spoerkel: IMOMAP.GE.1 to IMOMAP.NE.0
C                     if IMOMAP==-1, only the phase of orbitals is mapped
   10 IF(IMOMAP.NE.0) THEN
C        CHECK IF SAVED ORBITALS ARE STILL VALID.
         IF(LMOSAV .AND. MOL.GT.MOLGUG .AND. KEEPCI.NE.1
     1       .AND. IMOMAP.NE.3) THEN
            DEALLOCATE(CMOSAV)
            LMOSAV = .FALSE.
         ENDIF
         IF(LMOSAV) THEN
C           MAPPING IS ASSUMED SUCCESSFUL UNTIL PROVEN OTHERWISE
            IFMAP = 0
C           CONSISTENCY CHECK.
            IF(UBOUND(CMOSAV,1).NE.NORBS .OR.
     1         UBOUND(CMOSAV,2).NE.NMOS) THEN
               WRITE(NB6,620)
     1           UBOUND(CMOSAV,1), UBOUND(CMOSAV,2), NORBS, NMOS
               STOP
            ENDIF
C           MAP THE CURRENT ORBITALS ON THE PREVIOUS ONES.
            LMOMAP(1:NMOS) = .FALSE.
            LLMOMAP(1:NMOS) = .FALSE.
            NFAIL          = 0
            ALLOCATE(DOTPRO(LMX,LMX))
C           CALCULATE ALL DOTPRODUCTS
            DO I=1,NMOS
               DO J=1,NMOS
                  DOTPRO(I,J) = DDOT(NORBS, CMOSAV(1,I), 1, C(1,J), 1)
               ENDDO
            ENDDO
C           FIND BIGGEST OVERLAP AND ASSIGN
C           THIS MAPPING. REPEAT NMOS TIMES.
            DO K=1,NMOS
               dotmax=0.D0
               Imax=0
               Jmax=0
C              LOOP OVER SAVED ORBITALS. CYCLE IF ALREADY ASSIGNED.
               DO J=1,NMOS
                  IF(LMOMAP(J)) CYCLE
C                 LOOP OVER CURRENT ORBITALS. CYCLE IF ALREADY ASSIGNED.
                  DO I=1,NMOS
                     IF(LLMOMAP(I)) CYCLE
C                    SAVE BIGGEST OVERLAP.
                     if (dabs(DOTPRO(I,J)).GT.dabs(dotmax)) then
                        Imax=I
                        Jmax=J
                        dotmax=DOTPRO(I,J)
                     endif
                  ENDDO
               ENDDO
C              ASSIGN BIGGEST OVERLAP MAPPING.
               KMOMAP(Imax) = Jmax
               ADTPRO(Imax) = dotmax
               LMOMAP(Jmax) = .TRUE.
               LLMOMAP(Imax) = .TRUE.
               if (DABS(dotmax) .LT. 0.5D0) NFAIL=NFAIL+1
               IF(MOTYPE(Imax).NE.-1) THEN
                  IF(IMOMAP.GE.2.AND.DABS(dotmax).LT.PMAPTH)THEN
                     WRITE(NB6,580) Imax, dotmax
                     IFMAP = 1
C                          PREVENT AN UNNECESSARY CI CALCULATION
                     ICALL = -1
                  ENDIF
               ENDIF
C              REVERSE PHASE OF NEGATIVE OVERLAP.
               IF(dotmax.LT.0) C(:,Jmax) = -C(:,Jmax)
            ENDDO
C           IF A MAPPING OVERLAP IS LESS THAN 0.5, PRINT MAPPINGS.
            if (NFAIL.GT.0) THEN
               Write(*,'("WARNING: ",I5," MAPPINGS BELOW 0.5")') NFAIL
               Write(*,*)
               if (IOUTCI.GE.2) THEN
                  DO I=1,NMOS
                     Write(*,'("MAPPED ",2I5,F12.8)') I,KMOMAP(I),
     &                  DOTPRO(I,KMOMAP(I))
                  ENDDO
               END IF
            END IF
            NFAIL=0
            DEALLOCATE(DOTPRO)
            IF(IMOMAP.EQ.1 .AND. IOUTCI.GT.0) THEN
C              PRINT OUT MATCHED ACTIVE ORBITAL MAPPING INFORMATION
C              SIMILAR TO THE IMOMAP=2 PRINTOUT BUT UNMATCHED PRINTED LATER
               WRITE(NB6,545)
               DO I = 1, NMOS
                  IF(LMOMAP(I) .AND. MOTYPE(I).NE.-1) THEN
                     MAPSTA = 1
                     IF (KMOMAP(I).NE.I) MAPSTA = 2
                     WRITE(NB6,560) I, CMAPST(MAPSTA),
     1                    KMOMAP(I), ADTPRO(I)
                  ENDIF
               ENDDO
            ENDIF
            IF(IMOMAP.EQ.-1 .AND. IOUTCI.GT.0) THEN
C              PRINT OUT MATCHED ACTIVE ORBITAL PHASE MAPPING INFORMATION
C              SIMILAR TO THE IMOMAP=2 PRINTOUT BUT UNMATCHED PRINTED LATER
               WRITE(NB6,546)
               DO I = 1, NMOS
                  IF(LMOMAP(I) .AND. MOTYPE(I).NE.-1) THEN
                     MAPSTA = 1
                     IF (KMOMAP(I).NE.I) MAPSTA = 2
                     IF (ADTPRO(I).LT.0) THEN
                        WRITE(NB6,561) I, CMAPST(MAPSTA),
     1                    KMOMAP(I), "-1"
                     ELSE
                        WRITE(NB6,561) I, CMAPST(MAPSTA),
     1                    KMOMAP(I), " 1"
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
            IF(IMOMAP.GE.2 .AND. IOUTCI.GT.0) THEN
C              PRINT ACTIVE ORBITAL MAPPING INFORMATION
               WRITE(NB6,550)
               DO I=1, NMOS
                  IF(MOTYPE(I).NE.-1) THEN
                     MAPSTA=1
                     IF(DABS(ADTPRO(I)).LT.PMAPTH) THEN
                        MAPSTA=3
                     ELSE IF(KMOMAP(I).NE.I) THEN
                        MAPSTA=2
                     ENDIF
                     IF(DABS(ADTPRO(I)).GT.0.5D0) THEN
                        WRITE(NB6,560) I, CMAPST(MAPSTA),
     1                                 KMOMAP(I), ADTPRO(I)
                     ELSE
                        WRITE(NB6,570) I, CMAPST(MAPSTA),
     1                                 KMOMAP(I)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
C          IF IMOMAP=1 OR (IMOMAP>=2 AND MAPPING IS SUCCESSFUL)
C           UPDATE ORBITALS
            IF(IMOMAP.EQ.1 .OR. (IMOMAP.GE.2.AND.IFMAP.EQ.0)) THEN
C              UPDATE THE ACTIVE ORBITAL INDICES.
               IMOCI(1:NACTIV) = KMOMAP(IMOCI(1:NACTIV))
C              UPDATE MO TYPES.
               MOTYSV(1:NMOS)  = MOTYPE(1:NMOS)
               DO I=1,NMOS
                  MOTYPE(KMOMAP(I)) = MOTYSV(I)
               ENDDO
C              THROW AWAY THE OLD MOS.
               DEALLOCATE(CMOSAV)
               LMOSAV = .FALSE.
               IF (MOVO .EQ. -5) THEN
                  CALL PIMOCI (C,E,LM2,LM3,IOUTPI,.TRUE.)
               END IF

            ENDIF
            IF(IMOMAP.EQ.-1) THEN
               DEALLOCATE(CMOSAV)
               LMOSAV = .FALSE.
            ENDIF
         ENDIF
C        SAVE CURRENT MOS FOR NEXT CI RUN.
C        NOTE LMOSAV CAN ONLY BE TRUE IF IMOMAP>=2 AND MAPPING FAILED
         IF(.NOT. LMOSAV) THEN
            ALLOCATE(CMOSAV(NORBS, NMOS))
            CMOSAV = C(1:NORBS, 1:NMOS)
            LMOSAV = .TRUE.
            IF (IMOMAP.EQ.3) THEN
C              SAVE MOS TO FILE INSTEAD OF INTERNALLY.
               CALL MOMSAV(0, NORBS, NMOS, NACTIV, NCIREF,
     1               ICI1, ICI2, MOVO, MCIREF, NCISYM,
     2               MOTYPE, CMOSAV, IMOCI, ICIREF, IMOPOS,
     3               LMOSAV)
               DEALLOCATE(CMOSAV)
               LMOSAV = .FALSE.
            ENDIF
         ELSE
            WRITE(NB6,640)
         ENDIF
      ENDIF
C *** PRINTING SECTION.
      IF(ICALL.EQ.-1) RETURN
      IF(IOUTCI.GT.0) THEN
C        PRINT MO SYMMETRY AND OCCUPATION NUMBERS.
         WRITE(NB6,500)
         DO I=1,NMOS
            STAR = '-'
            IF(MOTYPE(I).GE.0) STAR = CHAR(MOTYPE(I) + 48)
            WRITE(NB6,510) I, E(I), LABELS(NSYM(I)), NSYM(I),
     1                     IMOSYM(I), OCCR(I), STAR
         ENDDO
C        ASSEMBLE RIGHT-JUSTIFIED SYMMETRY LABELS.
         DO I=1,NUMLAB
            DO J=4,1,-1
               IF(LABELS(I)(J:J).NE.' ') EXIT
            ENDDO
            RLABEL(I)        = ' '
            RLABEL(I)(5-J:4) = LABELS(I)(1:J)
         ENDDO
C        PRINT REFERENCE CONFIGURATIONS.
         WRITE(NB6,520)
         NBLK = (NACTIV+29)/30
         DO K=1,NBLK
            J1 = 30*K-29
            J2 = MIN(30*K, NACTIV)
            WRITE(NB6,'(/1X,30I4)') IMOCI(J1:J2)
            WRITE(NB6,'( 1X,30A4)') RLABEL(NSYM(IMOCI(J1:J2)))
            DO I=1,NCIREF
               WRITE(NB6,'(1X,30I4)') ICIREF(J1:J2, I)
            ENDDO
         ENDDO
      ENDIF
      RETURN
  500 FORMAT(///1X,'   MO     EIGENVALUE     ',
     1             'LABEL   NSYM  IMOSYM     OCC.    ACTIVE'/)
  510 FORMAT(   1X,I5,F15.5,A10,I6,I7,F12.5,A7)
  520 FORMAT(///1X,'REFERENCE CONFIGURATIONS')
  530 FORMAT(// 6X,'UNMATCHED ORBITALS (OVERLAP LESS THAN 0.5):',
     1       // 6X,'SAVED ORBITALS   CURRENT ORBITALS',
     2       /  6X,'--------------   ----------------',
     3       /  6X,'INDEX   ACTIVE   INDEX   ENERGY')
  540 FORMAT(   1X,I2,I7,I8,I9,F12.5)
  545 FORMAT(// 6X,'MATCHED ACTIVE ORBITALS (OVERLAP OVER 0.5):',
     1       // 6X,'SAVED     STATUS     CURRENT  OVERLAP',
     2       /  6X,'-------------------------------------')
  546 FORMAT(// 6X,'MATCHED ACTIVE ORBITALS PHASES:',
     1       // 6X,'SAVED     STATUS     CURRENT  PHASE',
     2       /  6X,'-------------------------------------')
  550 FORMAT(// 6X,'ACTIVE ORBITAL MAPPING:',
     1       // 6X,'SAVED     STATUS     CURRENT  OVERLAP',
     2       /  6X,'-------------------------------------')
  560 FORMAT(   6X,1X,I2,6X,A9,5X,I2,5X,F6.3)
  561 FORMAT(   6X,1X,I2,6X,A9,5X,I2,5X,A2)
  570 FORMAT(   6X,1X,I2,6X,A9,5X,I2,5X,' -----')
  580 FORMAT(/  1X,'MAPPING FAILED FOR ACTIVE ORBITAL ',I3,
     1             ' WITH OVERLAP ',F6.3)
  590 FORMAT(/  1X,'MAPPING FAILED FOR ACTIVE ORBITAL ',I3,
     1             ' WITH OVERLAP UNDER 0.5')
  610 FORMAT(/  1X,'ERROR: TOO FEW MOS WITH SYMMETRY ',A3,' (',I1,').')
  620 FORMAT(/  1X,'ERROR: INCONSISTENT MO COEFFICIENT MATRIX ',
     1             'DIMENSIONS.',
     2       /  1X,'PREVIOUS DIMENSIONS (ROWS, COLS):',2I6,
     3       /  1X,'CURRENT  DIMENSIONS (ROWS, COLS):',2I6)
  630 FORMAT(/  1X,'ERROR: TOO MANY UNMATCHED ORBITALS.')
  640 FORMAT(/  1X,'ERROR: ACTIVE ORBITAL MAPPING FAILED',
     1       /  1X,'CURRENT ORBITALS WILL NOT BE SAVED')
  700 FORMAT(/  1X,'ERROR: THERE SHOULD BE NO INTERNALLY SAVED MOS')
  710 FORMAT(/  1x,'WARNING: IMOMAP=3 MAPPING WILL NOT BE CARRIED OUT')
      END