      SUBROUTINE INQMMM (JPRINT)
C     *
C     INPUT OF EXTERNAL POINTS.
C     *
C     OPTIONS AVAILABLE IN THE STANDARD PROGRAM.
C     MMINP=1: CALCULATE PROPERTIES AT EXTERNAL POINTS.
C     MMINP=2: QM TREATMENT IN THE PRESENCE OF EXTERNAL POINT CHARGES.
C     *
C     THIS ROUTINE IS DESIGNED TO READ THE MINIMUM AMOUNT OF DATA
C     NEEDED TO TEST THE PROGRAM IN SUCH APPLICATIONS (MMINP=1,2).
C     IT IS NOT INTENDED AS A COMPLETE INPUT FOR QM/MM CALCULATIONS.
C     *
      USE LIMIT, ONLY: LM1, LMZ, LM1M
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./IJWORK/ LTEMP(LM1),LTEMP1(LM1*11)
     ./INOPT1/ IN1(300)
     ./INOPT2/ IN2(300)
     ./MOPAC / IMOPAC
     ./NBFILE/ NBF(20)
     ./PARM1 / A(3,LM1),NABC(LM1*4),NATOMS
     ./QMMM1 / COORDM(3,LM1M),CHARGM(LM1M)
     ./QMMM2 / LINK(LM1)
     ./QMMM3 / DELTAM(LMZ),OMEGAM(LMZ)
     ./QMMM6 / ISELCT(LM1+LM1M)
C    ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
C    ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
C *** FILE NUMBERS.
      NB5    = NBF(5)
      NB6    = NBF(6)
C *** INITIALIZATION.
      IOP    = IN2(2)
      IFORM  = IN2(5)
      IPAROK = IN2(11)
      MMINP  = IN2(12)
      INGEOM = IN2(54)
      KSTOP  = 0
      LINK(1:NUMAT) = 0
C *** CHECK OPTIONS.
      IF(MMINP.LE.0 .OR. MMINP.GT.2) THEN
         WRITE(NB6,500) MMINP
         KSTOP = KSTOP+1
      ENDIF
C *** READ CONTROL VARIABLES.
C     121-127: NUMATM,MMCOUP,MMPOT,MMLINK,NLINK,MMFILE,MCHARG,MMSKIP
      IF(IMOPAC.EQ.0) THEN
         IF(IFORM.LE.0) THEN
            READ(NB5,700) (IN1(I),I=120,127)
         ELSE
            READ(NB5,*)   (IN1(I),I=120,127)
         ENDIF
      ENDIF
      IN2(120:127) = IN1(120:127)
C *** CHECK CONTROL VARIABLES.
      IF(IN2(120).LE.0 .OR. IN2(120).GT.LM1M) THEN
         WRITE(NB6,510) IN2(120)
         KSTOP = KSTOP+1
      ENDIF
      IF(IN2(121).LT.0 .OR. IN2(121).GT.4) THEN
         WRITE(NB6,520) IN2(121)
         KSTOP = KSTOP+1
      ENDIF
      IF(IN2(122).LT.0 .OR. IN2(122).GT.8) THEN
         WRITE(NB6,530) IN2(122)
         KSTOP = KSTOP+1
      ENDIF
      IF(IN2(123).LT.0 .OR. IN2(123).GT.3) THEN
         WRITE(NB6,540) IN2(123)
         KSTOP = KSTOP+1
      ENDIF
      IF(IN2(124).LT.0 .OR. IN2(124).GT.NUMAT) THEN
         WRITE(NB6,550) IN2(124)
         KSTOP = KSTOP+1
      ENDIF
      IF(IN2(127).LT.0 .OR. IN2(127).GT.2) THEN
         WRITE(NB6,560) IN2(127)
         KSTOP = KSTOP+1
      ENDIF
C *** IMPOSE DEFAULT VALUES.
      IF(IN2(121).GE.2 .AND. IN2(122).EQ.0) IN2(122)=4
      IF(IN2(121).GT.0 .AND. IN2(123).EQ.0) IN2(123)=1
      IF(IN2(124).EQ.0 .AND. IN2(123).GT.0) IN2(123)=0
      NUMATM = IN2(120)
C     MMCOUP = IN2(121)
      MMPOT  = IN2(122)
C     MMLINK = IN2(123)
      NLINK  = IN2(124)
      MMFILE = IN2(125)
C     MCHARG = IN2(126)
      MMSKIP = IN2(127)
C *** READ INFORMATION ON LINK ATOMS.
C     LTEMP(I): NUMBERS OF LINK ATOMS IN FULL LIST (WITH DUMMY ATOMS).
C     LINK(I) : LABELS GENERATED (LINK(I)=1 IF ATOM I IS A LINK ATOM).
C     LINK(I) : REFERS TO LIST OF REAL ATOMS (WITHOUT DUMMY ATOMS).
      IF(NLINK.GT.0) THEN
         IF(IFORM.LE.0) THEN
            READ(NB5,700) (LTEMP(I),I=1,NLINK)
         ELSE
            READ(NB5,*)   (LTEMP(I),I=1,NLINK)
         ENDIF
         DO 30 I=1,NLINK
         LINKQM = LTEMP(I)
         IF(LINKQM.GT.0 .AND. LINKQM.LE.NATOMS) THEN
            LINK(LINKQM) = 1
         ELSE
            WRITE(NB6,570) LINKQM
            KSTOP = KSTOP+1
         ENDIF
   30    CONTINUE
C        REMOVE DUMMY ATOMS.
         IF(KSTOP.EQ.0 .AND. NUMAT.NE.NATOMS) THEN
            II  = 0
            DO 40 I=1,NATOMS
            IF(NAT(I).LT.99) THEN
               II  = II+1
               IF(LINK(I).EQ.1 .AND. II.LT.I) THEN
                  LINK(II) = LINK(I)
                  LINK(I)  = 0
               ENDIF
            ENDIF
   40       CONTINUE
         ENDIF
      ENDIF
C *** ERROR EXIT.
      IF(KSTOP.GT.0) THEN
         WRITE(NB6,580)
         STOP 'INQMMM'
      ENDIF
C *** INITIALIZATION FOR FLAGS THAT LABEL FIXED MM ATOMS.
C     DEFAULT: NO FIXED MM ATOMS - COMPUTE FULL GRADIENT.
      ISELCT(1:NUMAT+NUMATM) = 0
      SKIPX=0.0D0
      SKIPY=0.0D0
      SKIPZ=0.0D0
      SKIPRAD=0.0D0
      IF(MMSKIP.EQ.2) THEN
         IF(IFORM.EQ.0) THEN
            READ(NB5,'(3F12.7,F12.4)') SKIPX,SKIPY,SKIPZ,SKIPRAD
         ELSE
            READ(NB5,*) SKIPX,SKIPY,SKIPZ,SKIPRAD
         ENDIF
         WRITE(NB6,*)
         WRITE(NB6,'("Freeze sphere: ",4F12.5)')
     1                       SKIPX,SKIPY,SKIPZ,SKIPRAD
      ENDIF
C *** SPECIFY INPUT FILE FOR COORDINATES AND CHARGES.
      IF(MMFILE.EQ.0) THEN
         NB = NBF(5)
      ELSE
         NB = NBF(20)
      ENDIF
C *** READ CARTESIAN COORDINATES AND CHARGES OF EXTERNAL POINTS.
      IF(MMFILE.GE.0 .AND. INGEOM.NE.-1) THEN
C        IF-STATEMENT KEPT FOR BACKWARD COMPATIBILITY.
         IF(IFORM.EQ.0) THEN
            IF(IN2(127).NE.1) THEN
               DO 50 I=1,NUMATM
               READ(NB,710,ERR=100,END=100)
     1                    (COORDM(J,I),J=1,3),CHARGM(I)
               DIST=DSQRT((COORDM(1,I)-SKIPX)**2+
     1                    (COORDM(2,I)-SKIPY)**2+
     2                    (COORDM(3,I)-SKIPZ)**2)
               IF (IN2(127).EQ.2.AND.DIST.GT.SKIPRAD) ISELCT(NUMAT+I)=1
   50          CONTINUE
            ELSE
               DO 60 I=1,NUMATM
               READ(NB,710,ERR=100,END=100)
     1           (COORDM(J,I),J=1,3),CHARGM(I),ISELCT(NUMAT+I)
   60          CONTINUE
            ENDIF
         ELSE
            IF(IN2(127).NE.1) THEN
               DO 70 I=1,NUMATM
               READ(NB,*,ERR=100,END=100)
     1                    (COORDM(J,I),J=1,3),CHARGM(I)
               DIST=DSQRT((COORDM(1,I)-SKIPX)**2+
     1                    (COORDM(2,I)-SKIPY)**2+
     2                    (COORDM(3,I)-SKIPZ)**2)
               IF (IN2(127).EQ.2.AND.DIST.GT.SKIPRAD) ISELCT(NUMAT+I)=1
   70          CONTINUE
            ELSE
               DO 80 I=1,NUMATM
               READ(NB,*,ERR=100,END=100)
     1           (COORDM(J,I),J=1,3),CHARGM(I),ISELCT(NUMAT+I)
   80          CONTINUE
            ENDIF
         ENDIF
      ENDIF
C *** INITIALIZATION FOR ELECTROSTATIC POTENTIAL.
      IF(MMPOT.EQ.0 .OR. MMPOT.GE.3) THEN
         CALL INIPOT (DELTAM,OMEGAM,LMZ,IOP,MMPOT)
         IF(IPAROK.GE.1 .AND. IPAROK.LE.3) CALL MODPOT
      ENDIF
C *** GSBP: READ GAMMA AND OMEGA MATRICES
      IF(IN2(121).EQ.4) THEN
           CALL  INGSBP
      ENDIF
C *** PRINTING SECTION.
      IF(JPRINT.GE.0) THEN
         WRITE(NB6,600)
         WRITE(NB6,610) (IN2(I),I=120,128)
         IF(MMFILE.GE.0) THEN
            IF(NUMATM.LE.100 .OR. JPRINT.GE.5) THEN
               WRITE(NB6,620)
               DO 90 I=1,NUMATM
               WRITE(NB6,630) I,(COORDM(J,I),J=1,3),CHARGM(I),
     1                        ISELCT(NUMAT+I)
   90          CONTINUE
            ENDIF
         ENDIF
         IF(NLINK.GT.0) THEN
            WRITE(NB6,640) (LTEMP(I),I=1,NLINK)
         ENDIF
      ENDIF
      RETURN
C *** ERROR EXIT.
  100 CONTINUE
      STOP 'INQMMM'
  500 FORMAT(/  1X,'OPTION OUT OF RANGE: MMINP  =',I5)
  510 FORMAT(/  1X,'OPTION OUT OF RANGE: NUMATM =',I5)
  520 FORMAT(/  1X,'OPTION OUT OF RANGE: MMCOUP =',I5)
  530 FORMAT(/  1X,'OPTION OUT OF RANGE: MMPOT  =',I5)
  540 FORMAT(/  1X,'OPTION OUT OF RANGE: MMLINK =',I5)
  550 FORMAT(/  1X,'OPTION OUT OF RANGE: NLINK  =',I5)
  560 FORMAT(/  1X,'OPTION OUT OF RANGE: MMSKIP =',I5)
  570 FORMAT(/  1X,'OPTION OUT OF RANGE: LINKQM =',I5)
  580 FORMAT(/  1X,'FATAL INPUT ERROR. STOP.')
  600 FORMAT(///1X,'*** OPTIONS FOR EXTERNAL POINTS AND QM/MM   ***',
     1       /  1X,'*** DEFINED EITHER EXPLICITLY OR BY DEFAULT ***')
  610 FORMAT(/  1X,'NUMATM =',I5,5X,'MMCOUP =',I5,5X,'MMPOT  =',I5,
     1          5X,'MMLINK =',I5,5X,'NLINK  =',I5,
     2       /  1X,'MMFILE =',I5,5X,'MCHARG =',I5,5X,'MMSKIP =',I5)
  620 FORMAT(///1X,'CARTESIAN COORDINATES (ANGSTROM) AND CHARGES (E)',
     1          1X,'OF THE EXTERNAL POINTS.',
     2       // 5X,'I',10X,'X',11X,'Y',11X,'Z',9X,'CHARGE   SKIPPED'/)
  630 FORMAT(1X,I5,2X,4F12.5,I8)
  640 FORMAT(///1X,'NUMBERS OF LINK ATOMS IN THE LIST OF QM ATOMS.',
     1       //(1X,16I5))
  700 FORMAT(16I5)
  710 FORMAT(3F12.7,F8.4,I3)
      END
