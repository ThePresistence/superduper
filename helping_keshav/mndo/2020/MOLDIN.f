      SUBROUTINE MOLDIN (NB6,NB13,NUMAT,I3N)
C     *
C     READ MOLDEN DATA FROM FILE NB13 AND WRITE THEM TO FILE NB93.
C     FILE NB93 'molden.dat' WILL BE READ BY MOLDEN.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 AZ
      CHARACTER*20 AWHERE
      CHARACTER*40 ATITLE
      CHARACTER*2  ATOMG,ATOMCG
      LOGICAL LA, LPRINT
      PARAMETER (THRVIB=1.0D0)
C     ALLOCATABLE WORKING ARRAYS FOR SAVING MOLDEN DATA.
      DOUBLE PRECISION, ALLOCATABLE :: G(:)
      DOUBLE PRECISION, ALLOCATABLE :: CG(:)
      DOUBLE PRECISION, ALLOCATABLE :: GGEO(:)
      DOUBLE PRECISION, ALLOCATABLE :: CGGEO(:)
C     AUTOMATIC FORTRAN90 ARRAYS FOR LABELS AND LOGICAL FLAGS.
      DIMENSION ATOMG(NUMAT),ATOMCG(NUMAT),LA(I3N)
C
C     FILE HANDLING.
      NB93 = 93
      OPEN (NB93, FILE='molden.dat', STATUS='UNKNOWN')
      WRITE (NB93,'(a)') '[Molden Format] '
      REWIND (NB13)
C     INITIALIZATION.
      NG     = 0
      NCG    = 0
      NGGEO  = 0
      NCGGEO = 0
      NATOMS = 0
      LPRINT = .FALSE.
C
C *** FIRST PASS OVER FILE NB13: MAIN WORK DONE HERE.
C     ERROR CHECKS AND OUTPUT TO FILE NB93 FOR MOST KEYWORDS.
C     [GEOCONV] AND [GEOMETRIES] ARE HANDLED IN A SECOND PASS
C     BECAUSE THE NUMBER OF AVAILABLE GRADIENTS (NG, NCG) AND
C     GEOMETRIES (NGGEO, NCGGEO) MUST BE DETERMINED FIRST.
C
C     READ LINE FROM FILE NB13 WITH MOLDEN DATA.
   10 READ(NB13,'(A)',ERR=100,END=100) AZ
C
C *** [GEOCONV]
      IF(AZ(1:9).EQ.'[GEOCONV]') THEN
        LPRINT = .FALSE.
        IF(AZ(12:13).EQ.'CG') THEN
           NCG = NCG + 1
        ELSE IF(AZ(12:13).EQ.'G ') THEN
           NG  = NG + 1
        ENDIF
        GO TO 10
      ENDIF
C
C *** [GEOMETRIES]
      IF(AZ(1:12).EQ.'[GEOMETRIES]') THEN
         LPRINT = .FALSE.
         IF(AZ(15:16).EQ.'CG') THEN
            NCGGEO = NCGGEO + 1
            READ(AZ(23:),*) NUMATS
            IF(NATOMS.EQ.0) THEN
               NATOMS = NUMATS
            ELSE IF(NUMATS.NE.NATOMS) THEN
               WRITE(NB6,*) ' Inconsistent number of atoms '
            ENDIF
         ELSE IF(AZ(15:16).EQ.'G ') THEN
            NGGEO = NGGEO + 1
            READ(AZ(23:),*) NUMATS
            IF(NATOMS.EQ.0) THEN
               NATOMS = NUMATS
            ELSE IF(NUMATS.NE.NATOMS) THEN
               WRITE(NB6,*) ' Inconsistent number of atoms '
               WRITE(NB6,*) ' Please check fort.', NB13
            ENDIF
         ENDIF
         GO TO 10
      ENDIF
C
C *** [Atoms]
      IF(AZ(1:7).EQ.'[Atoms]') THEN
         WRITE(NB93,'(a)') AZ
         LPRINT = .TRUE.
         GO TO 10
      ENDIF
C
C *** [STO]
      IF(AZ(1:5).EQ.'[STO]') THEN
         WRITE(NB93,'(a)') AZ
         LPRINT = .TRUE.
         GO TO 10
      ENDIF
C
C *** [MO]
      IF(AZ(1:4).EQ.'[MO]') THEN
         LPRINT = .FALSE.
         IF(AZ(6:10).EQ.'     ') THEN
            WRITE(NB93,'(a)') AZ
            LPRINT = .TRUE.
        ENDIF
        GO TO 10
      ENDIF
C
C *** [FR-COORD]
      IF(AZ(1:10).EQ.'[FR-COORD]') THEN
        LPRINT = .FALSE.
        READ(AZ(11:),*) NATFR
        IF(NATFR.NE.NATOMS) THEN
           WRITE(NB6,*) ' Inconsistency detected (in [FR-COORD]):',
     1                  ' different number of atoms '
           WRITE(NB6,*) ' check fort.', NB13
           GO TO 999
        ENDIF
        WRITE(NB93,'(A)') AZ
        DO 20 I=1,NATFR
        READ(NB13,'(A)',END=190) AZ
        WRITE(NB93,'(A)') AZ
        IF(AZ(1:1).EQ.'[') THEN
           AWHERE = '[FR-COORD]'
           GO TO 990
        ENDIF
   20   CONTINUE
        GO TO 10
  190   AWHERE = '[FR-COORD]'
        GO TO 990
      ENDIF
C
C *** [FREQ]
      IF(AZ(1:6).EQ.'[FREQ]') THEN
        LPRINT = .FALSE.
        READ(AZ(7:),*) NFR
        IF(NFR.NE.3*NATFR) THEN
           WRITE(NB6,*) ' Inconsistency detected (in [FREQ]): ',
     1                  ' NFR.ne.3*NATFR ', NFR, 3*NATFR
        ENDIF
        WRITE(NB93,'(A)') AZ
        DO 30 I=1,NFR
        AWHERE = '[FREQ]'
        READ(NB13,*,ERR=290,END=290) FREQ
        IF(ABS(FREQ).LT.THRVIB) THEN
           LA(I) = .FALSE.
        ELSE
           LA(I) = .TRUE.
           WRITE(NB93,'(F20.10)') FREQ
        ENDIF
   30   CONTINUE
        GO TO 10
  290   AWHERE = '[FREQ]'
        GO TO 990
      ENDIF
C
C *** [FR-NORM-COORD]
      IF(AZ(1:15).EQ.'[FR-NORM-COORD]') THEN
         LPRINT = .FALSE.
         READ(AZ(16:),*) NNFR
         IF(NNFR.NE.NFR) THEN
            WRITE(NB6,*)
     1        ' Inconsistency detected (in [FR-NORM-COORD]): ',
     2        ' NNFR.ne.NFR ', NNFR, NFR
           WRITE(NB6,*) ' Please check fort.', NB13
           GO TO 999
         ENDIF
         WRITE(NB93,'(A)') AZ
         NVIB = 0
         DO 50 I=1,NNFR
            READ(NB13,'(A)',END=390) AZ
            IF(AZ(1:1).EQ.'[') THEN
               AWHERE = '[FR-NORM-COORD]'
               GO TO 990
            ENDIF
            IF(LA(I)) THEN
               NVIB = NVIB + 1
               WRITE(NB93,'(A,I4)') 'vibration  ', NVIB
            ENDIF
            DO 40 J=1,NNFR/3
            READ(NB13,'(A)',END=390) AZ
            IF(LA(I)) THEN
               WRITE(NB93,'(A)') AZ
            ENDIF
            IF(AZ(1:1) .EQ. '[') THEN
               AWHERE = '[FR-NORM-COORD]'
               GO TO 990
            ENDIF
   40   CONTINUE
   50   CONTINUE
        IF(NVIB+6.NE.NNFR .AND. NVIB+5.NE.NNFR) THEN
           WRITE(NB6,*) ' NUMBER OF MODES WITH ''ZERO'' FREQUENCY ',
     1        'INCORRECT - CHECK fort.', NB13, ' AND/OR ADJUST THRVIB '
           GO TO 999
        ENDIF
        GO TO 10
  390   AWHERE = '[FR-NORM-COORD]'
        GO TO 990
      ENDIF
C
C *** NO 'KEYWORD' FOUND.
      IF(LPRINT) THEN
         WRITE(NB93,'(A)') AZ
      ENDIF
      GO TO 10
C
C *** SECOND PASS OVER FILE NB13: REMAINING WORK DONE HERE.
  100 CONTINUE
C     ASSIGN ARRAY DIMENSIONS.
      NAG    = 3*NG
      NACG   = 3*NCG
      NAGEO  = 3*NATOMS*NGGEO
      NACGEO = 3*NATOMS*NCGGEO
      NAGSUM = NAG+NACG+NAGEO+NACGEO
      IF(NAGSUM.EQ.0) GO TO 999
C     ALLOCATE WORKING ARRAYS.
      IF(NAG.GT.0) ALLOCATE (G(NAG),STAT=IER1)
      IF(NACG.GT.0) ALLOCATE (CG(NACG),STAT=IER2)
      IF(NAGEO.GT.0) ALLOCATE (GGEO(NAGEO),STAT=IER3)
      IF(NACGEO.GT.0) ALLOCATE (CGGEO(NACGEO),STAT=IER4)
      IF(IER1.NE.0 .OR. IER2.NE.0 .OR. IER3.NE.0 .OR. IER4.NE.0) THEN
         WRITE(NB6,*) ' ALLOCATE failed in MOLDIN '
         WRITE(NB6,*) ' Needed ', NAG, NACG, NAGEO, NACGEO
      ENDIF
C     INITIALIZATION.
      IG     = 0
      ICG    = 0
      IGGEO  = 0
      ICGGEO = 0
      REWIND(NB13)
C
C *** READ ONE LINE FROM FILE NB13.
  110 READ(NB13,'(A)',END=200) AZ
C
C *** [GEOCONV]
      IF(AZ(1:9).EQ.'[GEOCONV]') THEN
         IF(AZ(12:13).EQ.'CG') THEN
            READ(AZ(14:),*) (CG(J+ICG),J=1,3)
            ICG = ICG + 3
         ELSE IF(AZ(12:13).EQ.'G ') THEN
            READ(AZ(14:),*) (G(J+IG),J=1,3)
            IG = IG + 3
        ENDIF
        GO TO 110
      ENDIF
C
C *** [GEOMETRIES]
      IF(AZ(1:12).EQ.'[GEOMETRIES]') THEN
        IF(AZ(15:16).EQ.'CG') THEN
           READ(NB13,'(A)') ATITLE
           DO 120 I=1,NATOMS
           READ(NB13,'(A)') AZ
           ATOMCG(I) = AZ(1:2)
           READ(AZ(6:),*) (CGGEO(J+ICGGEO),J=1,3)
           ICGGEO = ICGGEO + 3
  120      CONTINUE
        ELSE IF(AZ(15:16) .EQ. 'G ') THEN
           READ(NB13,'(A)') ATITLE
           DO 130 I=1,NATOMS
           READ(NB13,'(A)') AZ
           ATOMG(I) = AZ(1:2)
           READ(AZ(6:),*) (GGEO(J+IGGEO),J=1,3)
           IGGEO = IGGEO + 3
  130      CONTINUE
        ENDIF
        GO TO 110
      ENDIF
C
C *** NO RELEVANT 'KEYWORD' FOUND.
      GO TO 110
C
C *** SAVE REMAINING DATA ON FILE NB93.
  200 CONTINUE
      IF(ICG.GT.0) THEN
         NCG = ICG/3
         WRITE(NB93,'(A)') '[GEOCONV] '
         WRITE(NB93,'(A)') 'energy '
         WRITE(NB93,'(F20.10)') (CG(3*J-2),J=1,NCG)
         WRITE(NB93,'(A)') 'max-force '
         WRITE(NB93,'(F20.10)') (CG(3*J)  ,J=1,NCG)
         WRITE(NB93,'(A)') 'rms-force '
         WRITE(NB93,'(F20.10)') (CG(3*J-1),J=1,NCG)
      ENDIF
      IF(ICGGEO.GT.0 .AND. NATOMS.GT.0) THEN
         NCGGEO = ICGGEO/(3*NATOMS)
         WRITE(NB93,'(A)') '[GEOMETRIES] XYZ '
         DO 220 I=1,NCGGEO
         WRITE(NB93,'(I5)') NATOMS
         WRITE(NB93,'(A)') ATITLE
         DO 210 J=1,NATOMS
         JJ = 3*(NATOMS*(I-1)+J-1)
         WRITE(NB93,'(2X,A2,3F20.10)') ATOMCG(J), (CGGEO(K+JJ),K=1,3)
  210    CONTINUE
  220    CONTINUE
      ENDIF
C
C *** DEALLOCATE WORKING ARRAYS.
      IF(NAG.GT.0) DEALLOCATE (G,STAT=IERR)
      IF(NACG.GT.0) DEALLOCATE (CG,STAT=IERR)
      IF(NAGEO.GT.0) DEALLOCATE (GGEO,STAT=IERR)
      IF(NACGEO.GT.0) DEALLOCATE (CGGEO,STAT=IERR)
      RETURN
C *** ERROR SECTION.
  990 CONTINUE
      WRITE(NB6,*) ' Something went wrong !'
      WRITE(NB6,*) ' Too few lines in ', AWHERE, 'block  - ',
     1             ' please check fort.', NB13
  999 CONTINUE
      RETURN
      END
