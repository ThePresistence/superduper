      SUBROUTINE NMRAVG (ID1,NR1,NR2,LMF,IDNMR,SHIFT,XSHIFT,NUMAT)
C     *
C     AVERAGE NMR CHEMICAL SHIFTS FOR EQUIVALENT ATOMS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     ID1(7,*)  DEFINITION OF REFERENCE DATA (I).
C     NR1       INDEX OF FIRST REFERENCE DATUM FOR A GIVEN MOLECULE (I).
C     NR2       INDEX OF LAST  REFERENCE DATUM FOR A GIVEN MOLECULE (I).
C     LMF       COLUMN DIMENSION OF ID1 ARRAY (I).
C     IDNMR     PROPERTY INDEX FOR NMR CHEMICAL SHIFTS (I).
C               = 30  REFERENCE DATA FROM THE LIQUID PHASE.
C               = 31  REFERENCE DATA FROM THE GAS PHASE.
C     SHIFT     NMR CHEMICAL SHIFTS, NOT AVERAGED (I).
C     XSHIFT    NMR CHEMICAL SHIFTS, AVERAGED (O).
C     NUMAT     NUMBER OF ATOMS (I).
C     *
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ID1(7,LMF)
      DIMENSION SHIFT(2,NUMAT),XSHIFT(2,NUMAT)
      DIMENSION KSHIFT(LM1)
      IF(IDNMR.LT.30 .OR. IDNMR.GT.31) RETURN
C *** INITIALIZATION.
      ISHIFT = IDNMR-29
      DO 10 I=1,NUMAT
      XSHIFT(ISHIFT,I) = SHIFT(ISHIFT,I)
   10 CONTINUE
C *** PREPROCESSING OF NMR CHEMICAL SHIFTS.
C     AVERAGE OVER EQUIVALENT ATOMS.
      NSHIFT = 0
      MSHIFT = 0
      DO 20 NR=NR1,NR2
      IF(ID1(1,NR).EQ.IDNMR) THEN
         NSHIFT = NSHIFT+1
         IF(ID1(5,NR).GT.0 .AND. ID1(5,NR).NE.ID1(4,NR)) MSHIFT=MSHIFT+1
      ENDIF
   20 CONTINUE
C     NSHIFT: NUMBER OF REFERENCE DATA WHICH ARE CHEMICAL SHIFTS.
C     MSHIFT: NUMBER OF ATOMS WHICH HAVE BEEN DECLARED EQUIVALENT.
      IF(NSHIFT.GT.0 .AND. MSHIFT.GT.0) THEN
         DO 30 I=1,NUMAT
         KSHIFT(I) = 0
   30    CONTINUE
C        CHECK FOR EQUIVALENCIES.
         DO 60 NR=NR1,NR2
         IF(ID1(1,NR).NE.IDNMR) GO TO 60
         ID4 = ID1(4,NR)
         ID5 = ID1(5,NR)
         IF(KSHIFT(ID4).GT.0) GO TO 60
C        DETERMINE NUMBER OF EQUIVALENT ATOMS (NEQUIV).
         NEQUIV = 1
         SHIFTS = SHIFT(ISHIFT,ID4)
         DO 40 NRX=NR+1,NR2
         IF(ID1(1,NRX).NE.IDNMR) GO TO 40
         ID4X = ID1(4,NRX)
         ID5X = ID1(5,NRX)
         IF(ID5X.EQ.ID5) THEN
            NEQUIV = NEQUIV+1
            SHIFTS = SHIFTS+SHIFT(ISHIFT,ID4X)
            KSHIFT(ID4X) = ID5
         ENDIF
   40    CONTINUE
C        AVERAGE COMPUTED CHEMICAL SHIFTS (SHIFTS).
C        IMPOSE EQUIVALENT VALUES.
         IF(NEQUIV.GT.1) THEN
            SHIFTS = SHIFTS/DBLE(NEQUIV)
            DO 50 NRX=NR,NR2
            IF(ID1(1,NRX).NE.IDNMR) GO TO 50
            ID4X = ID1(4,NRX)
            ID5X = ID1(5,NRX)
            IF(ID5X.EQ.ID5) THEN
               XSHIFT(ISHIFT,ID4X) = SHIFTS
            ENDIF
   50       CONTINUE
         ENDIF
   60    CONTINUE
      ENDIF
      RETURN
      END