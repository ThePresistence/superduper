      SUBROUTINE INPATH (RC,LMR,LREACT,LTOTAL,A,NATOMS,IFORM,INT,BFACT)
C     *
C     INPUT OF DATA FOR REACTION PATH OR GRID.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL INT
      COMMON
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./PARM6 / LX,LY,LZ
      DIMENSION RC(LMR),A(3,NATOMS)
C *** FILE NUMBERS.
      NB5    = NBF(5)
      NB6    = NBF(6)
C *** CHECK FOR INTERPOLATED PATH.
      KGEOM  = IN2(68)
      IF(KGEOM.EQ.4) GO TO 100
C *** FIRST LINE OF INPUT.
      IF(IFORM.LE.0) THEN
         READ(NB5,500) L1,L2,L3,STEP
      ELSE
         READ(NB5,*)   L1,L2,L3,STEP
      ENDIF
C     CHECK FOR SIMPLE INPUT ERRORS.
      KSTOP  = 0
      IF(L1.LT.1 .OR. L1.GT.NATOMS) THEN
         KSTOP = KSTOP+1
         WRITE(NB6,600) L1,L2,L3,STEP
         WRITE(NB6,610) NATOMS
      ENDIF
      IF(L2.LT.1 .OR. L2.GT.3) THEN
         KSTOP = KSTOP+1
         IF(KSTOP.EQ.1) WRITE(NB6,600) L1,L2,L3,STEP
         WRITE(NB6,620)
      ENDIF
      IF(L3.LT.1 .OR. L3.GT.LMR-1) THEN
         KSTOP = KSTOP+1
         IF(KSTOP.EQ.1) WRITE(NB6,600) L1,L2,L3,STEP
         WRITE(NB6,630) LMR-1
      ENDIF
      IF(KSTOP.GT.0) STOP 'INPATH'
C     STORE THE INPUT DATA.
      LREACT = 3*(L1-1)+L2
      LTOTAL = L3+1
C     INITIAL POINT FROM GEOMETRY INPUT.
      RC(1)  = A(L2,L1)
C *** READ THE VALUES FOR THE REACTION COORDINATE.
      IF(STEP.EQ.0.D0) THEN
         IF(IFORM.LE.0) THEN
            READ(NB5,510) (RC(I),I=2,LTOTAL)
         ELSE
            READ(NB5,*)   (RC(I),I=2,LTOTAL)
         ENDIF
C        CONVERT TO RADIANS, IF NECESSARY.
         IF(INT .AND. L2.GT.1) THEN
            DO 10 I=2,LTOTAL
            RC(I)  = RC(I)*BFACT
   10       CONTINUE
         ENDIF
C *** CALCULATE THE VALUES FOR THE REACTION COORDINATE
C     FROM THE INITIAL POINT AND THE STEP SIZE.
      ELSE
         IF(INT .AND. L2.GT.1) STEP=STEP*BFACT
         DO 20 I=2,LTOTAL
         RC(I)  = RC(1)+(I-1)*STEP
   20    CONTINUE
      ENDIF
      RETURN
C *** INPUT FOR LINEAR INTERPOLATION BETWEEN TWO INPUT GEOMETRIES.
  100 CONTINUE
      IF(IFORM.LE.0) THEN
         READ(NB5,520) L1,LX,LY,LZ
      ELSE
         READ(NB5,*)   L1,LX,LY,LZ
      ENDIF
      IF(L1.LE.2) L1=11
      IF(LX.NE.-1) LX=1
      IF(LY.NE.-1) LY=1
      IF(LZ.NE.-1) LZ=1
C     CHECK FOR SIMPLE INPUT ERRORS.
      IF(L1.GT.LMR) THEN
         WRITE(NB6,700) L1
         WRITE(NB6,710) LMR
         STOP 'INPATH'
      ENDIF
C     CALCULATE THE FACTORS FOR LINEAR INTERPOLATION.
C     GEOMETRY (I) = (ONE-RC(I)) * GEOMETRY (1) + RC(I) * GEOMETRY (2)
      RC(1)  = 0.D0
      RC(L1) = 1.D0
      DELTA  = 1.D0/(L1-1)
      DO 110 I=2,L1-1
      RC(I)  = RC(I-1) + DELTA
  110 CONTINUE
      LREACT = L1
      LTOTAL = L1
      RETURN
  500 FORMAT(3I5,F10.5)
  510 FORMAT(8F10.5)
  520 FORMAT(4I5)
  600 FORMAT(///5X,'**** WRONG INPUT FOR REACTION PATH.',/5X,3I5,F10.5)
  610 FORMAT(/  5X,'FIRST VARIABLE OUTSIDE THE ALLOWED RANGE 1-',I3)
  620 FORMAT(/  5X,'SECOND VARIABLE OUTSIDE THE ALLOWED RANGE 1-3')
  630 FORMAT(/  5X,'THIRD VARIABLE OUTSIDE THE ALLOWED RANGE 1-',I3,
     1       /  5X,'TOO MANY POINTS ON PATH.')
  700 FORMAT(///5X,'**** WRONG INPUT FOR INTERPOLATION. L1 =',5X,I5)
  710 FORMAT(/  5X,'MAXIMUM NUMBER OF POINTS ON PATH: LMR =',5X,I5)
      END
