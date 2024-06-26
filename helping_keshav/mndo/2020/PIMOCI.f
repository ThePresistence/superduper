      SUBROUTINE PIMOCI (C,E,LM2,LM3,IOUT2,PRTONLY)
C     *
C     DEFINITION OF THE IMOCI ARRAY FOR THE CORRELATION TREATMENT.
C     *
C     MOVO=-1,-2,-3: PI CORRELATION IN PLANAR MOLECULES.
C     PI MOS BUILT FROM PX (MOVO=-1), PY (MOVO=-2), OR PZ (MOVO=-3).
C     *
C     MOVO=-4: D-SHELL CORRELATION IN TRANSITION METAL COMPOUNDS.
C     MOS WITH HIGHEST D POPULATIONS INCLUDED IN CI TREATMENT.
C     *
C     MOVO=-5: PSEUDO PI CORRELATION IN NON-PLANAR MOLECULES.
C     EACH CONJUGATED ATOM HAS AN ASSOCIATED LOCAL PI PLANE,
C     DEFINED BY THE ATOM AND ITS TWO NEAREST CONJUGATED NEIGHBOURS.
C     *
      USE LIMIT, ONLY: LM1, LMX, LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CICLO / CLOTHR
     ./CICONJ/ ICONJ(LM1)
     ./CIMOS / IMOCI(LMX)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./GUGA1 / NCIO,NCIGAM
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
      LOGICAL   LISPI(LM3)
      LOGICAL   LHASPO(LM1)
      LOGICAL   PRTONLY
      LOGICAL   MAKEPI
      DIMENSION C(LM2,LM3),E(LM3)
      DIMENSION IND(LMX),PIPOP(LMX),PIPOPA(LM1),UNIVEC(3,LM1)
      DIMENSION IPINN(2,LM1)
      DIMENSION ISORT(LMX)
      CHARACTER(LEN=30) PICHAR
      CHARACTER(LEN=1) INDCHAR
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** INPUT OPTIONS.
      ICI1   = IN2(131)
      ICI2   = IN2(132)
      MOVO   = IN2(134)
      JCI1   = IN2(151)
      JCI2   = IN2(152)
      IPIPOP = IN2(153)
      NCONJ  = IN2(167)
      IF(IPIPOP.GT.0) THEN
         TPIPOP = 1.0D-04 * DBLE(IPIPOP)
      ELSE
         TPIPOP = 0.0D0
      ENDIF
C *** ONLY PRINT PI POPULATIONS AND LEAVE ROUTINE.
      IF(PRTONLY) THEN
      IF(IOUT2.GE.1) THEN
         NHASPO = 0
         DO 8 I=1, NUMAT
            LHASPO(I) = .FALSE.
    8    CONTINUE
         DO 9 I=1, NUMAT
            IORBS  = NLAST(I)-NFIRST(I)+1
            IF (IORBS.GE.4) THEN
               LHASPO(I) = .TRUE.
               NHASPO = NHASPO + 1
            ENDIF
    9    CONTINUE
         CALL FINDNN(IOUT2, LHASPO, IPINN)
         DO 11 I=1,NUMAT
            UNIVEC(1,I) = ZERO
            UNIVEC(2,I) = ZERO
            UNIVEC(3,I) = ZERO
   11    CONTINUE
         CALL NRMPIP(IOUT2, LHASPO, IPINN, UNIVEC)
         DO 31 I=1,NUMAT
            PIPOPA(I) = ZERO
   31    CONTINUE
         DO 51 J=1,LM3
            PIPOP(J) = ZERO
            DO 41 I=1,NUMAT
               IA       = NFIRST(I)
               PPI      = UNIVEC(1,I)*C(IA+1,J) +
     1                    UNIVEC(2,I)*C(IA+2,J) +
     2                    UNIVEC(3,I)*C(IA+3,J)
               PIPOP(J) = PIPOP(J)+PPI**2
               IF(J.LE.NUMB) PIPOPA(I) = PIPOPA(I)+TWO*PPI**2
   41       CONTINUE
   51    CONTINUE
         NB6 = NBF(6)
         DO 121 J=1,LM3
            IND(J) = 0
  121    CONTINUE
         DO 131 I=1,ICI1+ICI2
            IND(IMOCI(I)) = 1
  131    CONTINUE
         WRITE(NB6,*)
         WRITE(NB6,551)
         DO 141 J=1,LM3
            PICHAR=""
            DO 142 I=1,FLOOR(PIPOP(J)*31.0D0)
               PICHAR=TRIM(ADJUSTL(PICHAR))//"="
  142       CONTINUE
            ITEMPPIPOP=MIN(MAX(FLOOR(PIPOP(J)*31.0D0),1),30)
            IF (PIPOP(J).LT.TPIPOP) PICHAR(ITEMPPIPOP:ITEMPPIPOP)="|"
            INDCHAR="-"
            IF (IND(J).EQ.1) INDCHAR="1"
            WRITE(NB6,561) J,E(J),PIPOP(J),PICHAR,INDCHAR
  141    CONTINUE
         IF(IOUT2.GE.2) THEN
            WRITE(NB6,570)
            DO 151 I=1,NUMAT
               WRITE(NB6,561) I,PIPOPA(I)
  151       CONTINUE
         END IF
      ENDIF
      GOTO 160
      ENDIF
C *** INITIALIZATION.
      IF(MOVO.GE.0 .OR. MOVO.LT.-5) THEN
         WRITE(NB6,600) MOVO
         STOP 'PIMOCI'
      ENDIF
      TNOPI = 0.01D0
C *** DETERMINE WHICH ATOMS ARE CONSIDERED WHEN CALCULATING PI/D-SHELL
C     POPULATIONS (LHASPO)
C     NCONJ = 0: ALL ATOMS WITH P ORBITALS (MOVO=-1/-2/-3/-5) OR
C                D ORBITALS (MOVO=-4) ARE CONSIDERED
C     NCONJ > 0: READ IN CONSIDERED ATOMS FROM ICONJ ARRAY
      NHASPO = 0
      DO 5 I=1, NUMAT
         LHASPO(I) = .FALSE.
    5 CONTINUE
      IF (NCONJ.EQ.0) THEN
         DO 6 I=1, NUMAT
            IORBS  = NLAST(I)-NFIRST(I)+1
            IF ((MOVO.NE.-4 .AND. IORBS.GE.4) .OR.
     1          (MOVO.EQ.-4 .AND. IORBS.GE.9)) THEN
               LHASPO(I) = .TRUE.
               NHASPO = NHASPO + 1
            ENDIF
    6    CONTINUE
      ELSE
         DO 7 I=1, NCONJ
            IF (ICONJ(I).LT.1 .OR. ICONJ(I).GT.NUMAT) THEN
               WRITE(NB6, 1000) ICONJ(I)
               STOP 'PIMOCI'
            ELSE IF (LHASPO(ICONJ(I))) THEN
               WRITE(NB6, 1010) ICONJ(I)
               STOP 'PIMOCI'
            ENDIF
            IORBS  = NLAST(ICONJ(I))-NFIRST(ICONJ(I))+1
            IF (MOVO.NE.-4 .AND. IORBS.LT.4) THEN
               WRITE(NB6, 1020) ICONJ(I)
               STOP 'PIMOCI'
            ELSE IF (MOVO.EQ.-4 .AND. IORBS.LT.9) THEN
               WRITE(NB6, 1030) ICONJ(I)
               STOP 'PIMOCI'
            ENDIF
            LHASPO(ICONJ(I)) = .TRUE.
            NHASPO = NHASPO + 1
    7    CONTINUE
      ENDIF
C *** FOR NON-PLANAR PI CORRELATION, DETERMINE 2 NEAREST NEIGHBOURS
C     FOR EACH CONJUGATED ATOM.
C     THESE THREE POINTS DEFINE THE LOCAL PI PLANE FOR EACH ATOM.
      IF (MOVO .EQ. -5) THEN
C        GET NEAREST NEIGHBOURS.
C        (CURRENTLY RECALCULATED WITH EVERY CALL TO PIMOCI, BUT FOR
C         MNDO-ONLY OPTIMISATIONS IT COULD BE WORTH CALLING IT ONLY
C         ONCE DURING INITIALISATION AND SAVING THE RESULTS, UNDER THE
C         ASSUMPTION THAT THE NEAREST NEIGHBOURS ARRAY IS UNLIKELY TO
C         CHANGE.)
         IF (NHASPO .LT. 3) THEN
            WRITE(NB6, 1040) NHASPO
            STOP 'PIMOCI'
         ENDIF
         CALL FINDNN(IOUT2, LHASPO, IPINN)
      ENDIF
C *** FOR PI CORRELATION, CALCULATE UNIT VECTOR ALONG PI AXIS
C     FOR EACH ATOM.
      DO 10 I=1,NUMAT
      UNIVEC(1,I) = ZERO
      UNIVEC(2,I) = ZERO
      UNIVEC(3,I) = ZERO
   10 CONTINUE
      IF (MOVO .GE. -3) THEN
C        CASE MOVO=-1. PLANAR MOLECULE. X AXIS AS PI AXIS.
C        CASE MOVO=-2. PLANAR MOLECULE. Y AXIS AS PI AXIS.
C        CASE MOVO=-3. PLANAR MOLECULE. Z AXIS AS PI AXIS.
         IPI    = ABS(MOVO)
         DO 20 I=1,NUMAT
            UNIVEC(IPI,I) = ONE
   20    CONTINUE
      ELSE IF (MOVO .EQ. -5) THEN
C        CASE MOVO=-5. NON-PLANAR MOLECULE. VARIABLE PI AXIS.
C        GET UNIVEC, THE UNIT NORMAL VECTORS TO THE PI PLANES.
         CALL NRMPIP(IOUT2, LHASPO, IPINN, UNIVEC)
      ENDIF
C *** CALCULATE POPULATIONS FOR ALL MOLECULAR ORBITALS.
C     CASE MOVO=-1,-2,-3,-5: PI POPULATIONS.
C     CASE MOVO=-4: D-SHELL POPULATIONS.
      DO 30 I=1,NUMAT
      PIPOPA(I) = ZERO
   30 CONTINUE
      DO 50 J=1,LM3
      PIPOP(J) = ZERO
      DO 40 I=1,NUMAT
      IA     = NFIRST(I)
      IF((MOVO.GE.-3 .OR. MOVO.EQ.-5) .AND. LHASPO(I)) THEN
         PPI = UNIVEC(1,I)*C(IA+1,J)+UNIVEC(2,I)*C(IA+2,J)
     1        +UNIVEC(3,I)*C(IA+3,J)
         PIPOP(J) = PIPOP(J)+PPI**2
         IF(J.LE.NUMB) PIPOPA(I) = PIPOPA(I)+TWO*PPI**2
      ELSE IF(MOVO.EQ.-4 .AND. LHASPO(I)) THEN
         PDELTA = C(IA+4,J)**2+C(IA+5,J)**2+C(IA+6,J)**2
     1           +C(IA+7,J)**2+C(IA+8,J)**2
         PIPOP(J) = PIPOP(J)+PDELTA
         IF(J.LE.NUMB) PIPOPA(I) = PIPOPA(I)+TWO*PDELTA
      ENDIF
   40 CONTINUE
   50 CONTINUE
C *** STORE INDICES OF ACTIVE ORBITALS IN IMOCI.
C     FIND OCCUPIED PI OR D-LIKE MOS USING AO POPULATIONS.
      NORB1 = 0
      DO 60 I=NUMB,1,-1
      IF(PIPOP(I).GE.TPIPOP) THEN
         NORB1 = NORB1 + 1
         IMOCI(NORB1) = I
      ENDIF
   60 CONTINUE
C     DEBUG PRINT: OCCUPIED SPECIAL ORBITAL INDICES.
      IF(IOUT2.GE.5) THEN
         WRITE(NB6,500) IMOCI(1:NORB1)
      ENDIF
C     SET A FLAG FOR EACH OCCUPIED SPECIAL ORBITAL.
      LISPI                 = .FALSE.
      LISPI(IMOCI(1:NORB1)) = .TRUE.
C     IF MORE SPECIAL ORBITALS HAVE BEEN REQUESTED THAN HAVE BEEN
C     IDENTIFIED, PRINT WARNING AND REDUCE THRESHOLD UNTIL ENOUGH
C     SPECIAL ORBITALS ARE IDENTIFIED (DOWN TO A MINIMUM OF TNOPI).
      IF(JCI1.GT.NORB1) THEN
         WRITE(NB6,650)
         CALL INDEXX(NUMB, PIPOP, ISORT)
         DO 65 I=NUMB,1,-1
            IF (PIPOP(ISORT(I)).LT.TNOPI) EXIT
            IF (.NOT. LISPI(ISORT(I))) THEN
               NORB1 = NORB1 + 1
               IMOCI(NORB1) = ISORT(I)
               LISPI(ISORT(I)) = .TRUE.
            ENDIF
            IF (JCI1.EQ.NORB1) EXIT
   65    CONTINUE
         REDPI1=PIPOP(IMOCI(NORB1))
         IF(JCI1.GT.NORB1) THEN
            WRITE(NB6,610) NORB1, JCI1
            JCI1 = NORB1
         ENDIF
      ENDIF
C     IF THE VALUE OF JCI1 IS USER-DEFINED,
C     USE JCI1 OCCUPIED SPECIAL ORBITALS.
      IF(JCI1.GT.0) THEN
         NORB1 = JCI1
      ELSE
         JCI1  = NORB1
      ENDIF
C     ADD SIGMA ORBITALS IF ICI1.GT.JCI1.
      DO 70 I=NUMB,1,-1
         IF(NORB1.GE.ICI1) EXIT
         IF(.NOT.LISPI(I)) THEN
            NORB1 = NORB1 + 1
            IMOCI(NORB1) = I
         ENDIF
   70 CONTINUE
C     WARN USER IF TOO MANY ACTIVE ORBITALS HAVE BEEN REQUESTED.
      IF(NORB1.LT.ICI1) THEN
         WRITE(NB6,630) ICI1, NORB1
         ICI1 = NORB1
      ENDIF
C     FIND VIRTUAL PI OR D-LIKE MOS USING AO POPULATIONS.
      NORB2 = 0
      DO 80 I=NUMB+1,NMOS
      IF(PIPOP(I).GE.TPIPOP) THEN
         NORB2 = NORB2 + 1
         IMOCI(ICI1+NORB2) = I
      ENDIF
   80 CONTINUE
C     DEBUG PRINT: VIRTUAL SPECIAL ORBITAL INDICES.
      IF(IOUT2.GE.5) THEN
         WRITE(NB6,510) IMOCI(ICI1+1:ICI1+NORB2)
      ENDIF
C     SET A FLAG FOR EACH VIRTUAL SPECIAL ORBITAL.
      LISPI(IMOCI(ICI1+1:ICI1+NORB2)) = .TRUE.
C     IF MORE SPECIAL ORBITALS HAVE BEEN REQUESTED THAN HAVE BEEN
C     IDENTIFIED, PRINT WARNING AND REDUCE THRESHOLD UNTIL ENOUGH
C     SPECIAL ORBITALS ARE IDENTIFIED (DOWN TO A MINIUM OF TNOPI).
      IF(JCI2.GT.NORB2) THEN
         WRITE(NB6,660)
         NVIRT = NMOS-NUMB
         CALL INDEXX(NVIRT, PIPOP(NUMB+1), ISORT)
         DO 85 I=NVIRT,1,-1
            IF (PIPOP(NUMB+ISORT(I)).LT.TNOPI) EXIT
            IF (.NOT. LISPI(NUMB+ISORT(I))) THEN
               NORB2 = NORB2 + 1
               IMOCI(ICI1+NORB2) = NUMB+ISORT(I)
               LISPI(NUMB+ISORT(I)) = .TRUE.
            ENDIF
            IF (JCI2.EQ.NORB2) EXIT
   85    CONTINUE
         REDPI2=PIPOP(IMOCI(ICI1+NORB2))
         IF(JCI2.GT.NORB2) THEN
            WRITE(NB6,620) NORB2, JCI2
            JCI2 = NORB2
         ENDIF
      ENDIF
C     IF THE VALUE OF JCI2 IS USER-DEFINED,
C     USE JCI2 VIRTUAL SPECIAL ORBITALS.
      IF(JCI2.GT.0) THEN
         NORB2 = JCI2
      ELSE
         JCI2  = NORB2
      ENDIF
C     ADD SIGMA* ORBITALS IF ICI2.GT.JCI2.
      DO 90 I=NUMB+1,NMOS
         IF(NORB2.GE.ICI2) EXIT
         IF(.NOT.LISPI(I)) THEN
            NORB2 = NORB2 + 1
            IMOCI(ICI1+NORB2) = I
         ENDIF
   90 CONTINUE
C     WARN USER IF TOO MANY ACTIVE ORBITALS HAVE BEEN REQUESTED.
      IF(NORB2.LT.ICI2) THEN
         WRITE(NB6,640) ICI2, NORB2
         ICI2 = NORB2
      ENDIF
C *** FINALLY SORT ACTIVE MOS IN ASCENDING ORDER.
      LM10   = ICI1+ICI2
      DO 110 I=1,LM10-1
      K = I
      DO 100 J=I+1,LM10
      IF(IMOCI(J).LT.IMOCI(K))  K = J
  100 CONTINUE
      IF(K.NE.I) THEN
         J        = IMOCI(I)
         IMOCI(I) = IMOCI(K)
         IMOCI(K) = J
      ENDIF
  110 CONTINUE
C *** PRINTING SECTION.
      IF(IOUT2.GE.0) THEN
         NB6 = NBF(6)
         WRITE(NB6,520)
         WRITE(NB6,530) (IMOCI(I),I=1,ICI1)
         WRITE(NB6,540)
         WRITE(NB6,530) (IMOCI(I),I=ICI1+1,LM10)
         IF(IOUT2.GE.1) THEN
            DO 120 J=1,LM3
            IND(J) = 0
  120       CONTINUE
            DO 130 I=1,LM10
            IND(IMOCI(I)) = 1
  130       CONTINUE
            WRITE(NB6,551)
            DO 140 J=1,LM3
               PICHAR=""
               DO 143 I=1,FLOOR(PIPOP(J)*31.0D0)
                  PICHAR=TRIM(ADJUSTL(PICHAR))//"="
  143          CONTINUE
               IF (J.LE.NUMB) THEN
                  REDPI=REDPI1
               ELSE
                  REDPI=REDPI2
               END IF
               ITEMPPIPOP=MIN(MAX(FLOOR(REDPI*31.0D0),1),30)
               IF (PIPOP(J).LT.REDPI) PICHAR(ITEMPPIPOP:ITEMPPIPOP)=":"
               ITEMPPIPOP=MIN(MAX(FLOOR(TPIPOP*31.0D0),1),30)
               IF (PIPOP(J).LT.TPIPOP) PICHAR(ITEMPPIPOP:ITEMPPIPOP)="|"
               INDCHAR="-"
               IF (IND(J).EQ.1) INDCHAR="1"
               WRITE(NB6,561) J,E(J),PIPOP(J),PICHAR,INDCHAR
  140       CONTINUE
            IF(IOUT2.GE.2) THEN
               WRITE(NB6,570)
               DO 150 I=1,NUMAT
               WRITE(NB6,561) I,PIPOPA(I)
  150          CONTINUE
            ENDIF
         ENDIF
      ENDIF
  160 RETURN
  500 FORMAT (// 1X,'OCCUPIED SPECIAL ORBITALS:'/1X,30I4)
  510 FORMAT (// 1X,'VIRTUAL SPECIAL ORBITALS:'/1X,30I4)
  520 FORMAT (// 1X,'IMOCI VALUES FOR OCCUPIED ORBITALS'/)
  530 FORMAT (   1X,10I5)
  540 FORMAT (// 1X,'IMOCI VALUES FOR UNOCCUPIED ORBITALS'/)
  551 FORMAT (// 1X,'    I  EIGENVALUE    PIPOP  ',12X,'PI-GRAPH',12X,
     1                     ' FLAG'/)
  561 FORMAT (   1X,I5,2X,2F10.5,3X,A30,3X,A1)
  570 FORMAT (// 1X,' ATOM   PIPOP-SUM'/)
  600 FORMAT (/  1X,'UNEXPECTED VALUE FOR MOVO (',I2,').')
  610 FORMAT (// 1X,'WARNING: TOO MANY SPECIAL ORBITALS REQUESTED.'
     1        /  1X,'OCCUPIED SPECIAL ORBITALS FOUND:    ',I5
     2        /  1X,'OCCUPIED SPECIAL ORBITALS REQUESTED:',I5
     3        /  1X,'RESETTING JCI1.')
  620 FORMAT (// 1X,'WARNING: TOO MANY SPECIAL ORBITALS REQUESTED.'
     1        /  1X,'VIRTUAL SPECIAL ORBITALS FOUND:    ',I5
     2        /  1X,'VIRTUAL SPECIAL ORBITALS REQUESTED:',I5
     3        /  1X,'RESETTING JCI2.')
  630 FORMAT (// 1X,'WARNING: TOO MANY ACTIVE ORBITALS REQUESTED.'
     1        /  1X,'OCCUPIED ACTIVE ORBITALS REQUESTED:',I5
     2        /  1X,'OCCUPIED ORBITALS PRESENT:         ',I5
     3        /  1X,'RESETTING ICI1.')
  640 FORMAT (// 1X,'WARNING: TOO MANY ACTIVE ORBITALS REQUESTED.'
     1        /  1X,'VIRTUAL ACTIVE ORBITALS REQUESTED:',I5
     2        /  1X,'VIRTUAL ORBITALS PRESENT:         ',I5
     3        /  1X,'RESETTING ICI2.')
  650 FORMAT (/  1X,'WARNING: TOO FEW OCCUPIED ORBITALS ABOVE PIPOP.',
     1        /  1X,'         PIPOP THRESHOLD WILL BE REDUCED.')
  660 FORMAT (/  1X,'WARNING: TOO FEW VIRTUAL ORBITALS ABOVE PIPOP.',
     1        /  1X,'         PIPOP THRESHOLD WILL BE REDUCED.')
 1000 FORMAT (// 1X,'ERROR: ICONJ ATOM OUT OF RANGE:',I5)
 1010 FORMAT (// 1X,'ERROR: ICONJ ATOM REQUESTED TWICE:',I5)
 1020 FORMAT (// 1X,'ERROR: ICONJ ATOM HAS NO P ORBITALS:',I5)
 1030 FORMAT (// 1X,'ERROR: ICONJ ATOM HAS NO D ORBITALS:',I5)
 1040 FORMAT (// 1X,'ERROR: AT LEAST THREE CONJUGATED ATOMS'
     1        /  1X,'       REQUIRED FOR MOVO=-5. CURRENTLY:',I5)
      END


      SUBROUTINE FINDNN(IOUT2, LHASPO, IPINN)
C     *
C     FOR EVERY ATOM WITH P ORBITALS, FIND THE TWO NEAREST NEIGHBOURS
C     THAT ALSO HAVE P ORBITALS.
C     TOGETHER THE THESE THREE POINTS DETERMINE THE LOCAL PI PLANE
C     FOR A GIVEN ATOM, FOR USE WITH THE MOVO=-5 OPTION
C     (PI CORRELATION IN NON-PLANAR MOLECULES).
C     *
C     THIS ROUTINE USES A NAIVE LINEAR SEARCH ALGORITHM.
C     THIS IS FINE FOR SMALL SYSTEMS, BUT SOMETHING MORE SOPHISTICATED
C     MIGHT BE REQUIRED FOR LARGER SYSTEMS.
C     *
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./NBFILE/ NBF(20)
      LOGICAL LHASPO(LM1)
      DIMENSION IPINN(2,LM1)
      DIMENSION DPINN(2,LM1)
C *** FILE NUMBERS
      NB6    = NBF(6)
C *** LOOP OVER ATOMS
      DO 20 I=1, NUMAT
C        INITIALISATION - NO NEIGHBOURS FOUND YET
         IPINN(1,I) = 0
         IPINN(2,I) = 0
         DPINN(1,I) = ZERO
         DPINN(2,I) = ZERO
         IF (LHASPO(I)) THEN
C           START OF LINEAR SEARCH LOOP
            DO 10 J=1, NUMAT
               IF (J .EQ. I) CYCLE
               IF (.NOT. LHASPO(J)) CYCLE
               DIST = SQRT( (COORD(1,I)-COORD(1,J))**2
     1                     +(COORD(2,I)-COORD(2,J))**2
     2                     +(COORD(3,I)-COORD(3,J))**2 )
c              DETERMINE IF IT IS A NEAREST NEIGHBOUR
               IF (IPINN(1,I).EQ.0 .OR. DIST.LT.DPINN(1,I)) THEN
                  IPINN(2,I) = IPINN(1,I)
                  DPINN(2,I) = DPINN(1,I)
                  IPINN(1,I) = J
                  DPINN(1,I) = DIST
               ELSE IF (IPINN(2,I).EQ.0 .OR. DIST.LT.DPINN(2,I)) THEN
                  IPINN(2,I) = J
                  DPINN(2,I) = DIST
               ENDIF
 10         CONTINUE
            IF (IPINN(1,I).EQ.0) THEN
               WRITE(NB6, 1000) I
               STOP 'FINDNN IN PIMOCI'
            ELSE IF (IPINN(2,I).EQ.0) THEN
               WRITE(NB6, 1010) I
               STOP 'FINDNN IN PIMOCI'
            ENDIF
         ENDIF
 20   CONTINUE
C *** PRINTING SECTION
      IF (IOUT2 .GE. 2) THEN
         WRITE(NB6, 500)
         DO 100 I=1, NUMAT
            IF (LHASPO(I)) THEN
               WRITE(NB6, 510) I, IPINN(1,I), DPINN(1,I),
     1                            IPINN(2,I), DPINN(2,I)
            ENDIF
 100     CONTINUE
      ENDIF
      RETURN
 500  FORMAT (// 1X,'NON-PLANAR PI CORRELATION SETUP'
     1         / 1X,'NEAREST NEIGHBOURS FOR RELEVANT ATOMS'
     1        // 1X,'  ATOM    NN1    DIST      NN2    DIST'/)
 510  FORMAT(    1X,I5,2X,I5,F10.5,2X,I5,F10.5)
 1000 FORMAT( // 1X,'ERROR: NO NEAREST NEIGHBOUR FOUND FOR ATOM',I5)
 1010 FORMAT( // 1X,'ERROR: NO 2ND NEAREST NEIGHBOUR FOUND FOR ATOM',I5)
      END


      SUBROUTINE NRMPIP(IOUT2, LHASPO, IPINN, UNIVEC)
C     *
C     FIND THE UNIT NORMAL VECTORS FOR ALL PI PLANES.
C     *
C     PI PLANES ARE DEFINED FROM THE NEAREST NEIGHBOUR ARRAY (IPINN)
C     AND THE CURRENT GEOMETRY.
C     UNIT VECTORS ARE RETURNED IN THE UNIVEC ARRAY.
C     *
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./NBFILE/ NBF(20)
      LOGICAL LHASPO(LM1)
      DIMENSION IPINN(2,LM1)
      DIMENSION UNIVEC(3,LM1)
      DIMENSION VA(3), VB(3), VNRM(3)
C *** FILE NUMBERS
      NB6    = NBF(6)
C *** LOOP OVER ATOMS
      DO 40 I=1, NUMAT
         IF (LHASPO(I)) THEN
C           VECTORS DEFINING THE PLANE
            DO 10 L=1,3
               VA(L) = COORD(L,IPINN(1,I)) - COORD(L,I)
 10         CONTINUE
            DO 20 L=1,3
               VB(L) = COORD(L,IPINN(2,I)) - COORD(L,I)
 20         CONTINUE
C           NORMAL VECTOR
            VNRM(1) = VA(2)*VB(3) - VA(3)*VB(2)
            VNRM(2) = VA(3)*VB(1) - VA(1)*VB(3)
            VNRM(3) = VA(1)*VB(2) - VA(2)*VB(1)
C           NORMALISE
            DNORM = DNRM2(3, VNRM, 1)
C           STORE NORMALISED UNIT VECTOR IN UNIVEC
            DO 30 L=1,3
               UNIVEC(L,I) = (ONE/DNORM) * VNRM(L)
 30         CONTINUE
         ENDIF
 40   CONTINUE
C     *** PRINTING SECTION
      IF (IOUT2 .GE. 2) THEN
         WRITE(NB6, 500)
         DO 100 I=1, NUMAT
            IF (LHASPO(I)) THEN
               WRITE(NB6, 510) I, UNIVEC(1,I), UNIVEC(2,I), UNIVEC(3,I)
            ENDIF
 100     CONTINUE
      ENDIF
      RETURN
 500  FORMAT ( / 1X,'UNIT NORMAL VECTORS TO PI PLANES'
     1        // 1X,'  ATOM       X         Y         Z'/)
 510  FORMAT(    1X,I5,2X,3F10.5)
      END
