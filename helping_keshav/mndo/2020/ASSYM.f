      SUBROUTINE ASSYM (NUMAT,AMS,COORD,ISUB,MSYM,IEL,ISYM,JPRINT)
*VOCL TOTAL,SCALAR
C     *
C     ASSIGN ABELIAN POINT GROUP (ISUB), DEFINE UP TO THREE (MSYM)
C     ASSOCIATED CHARACTERISTIC SYMMETRY ELEMENTS (IEL), AND DETERMINE
C     THE SYMMETRY RELATIONS BETWEEN THE ATOMS (ISYM).
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     NUMAT     NUMBER OF ATOMS (I).
C     AMS()     ATOMIC MASSES (I).
C     COORD()   CARTESIAN COORDINATES OF THE ATOMS (I).
C     ISUB      INTERNAL LABEL (1..7) FOR ABELIAN POINT GROUP (O).
C     MSYM      NUMBER OF CHARACTERISTIC SYMMETRY ELEMENTS (O).
C     IEL(3)    INTERNAL LABELS FOR THESE SYMMETRY ELEMENTS (O).
C     ISYM(3,*) PERMUTED TABLE OF ATOM NUMBERS AFTER APPLYING
C               THE CORRESPONDING SYMMETRY OPERATIONS (O).
C     JPRINT    PRINTING FLAG (I).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 NGROUP(7)
      CHARACTER*2 NEL(7)
      COMMON /NBFILE/ NBF(20)
      DIMENSION AMS(NUMAT),COORD(3,NUMAT)
      DIMENSION IEL(3),ISYM(3,NUMAT)
      DIMENSION CG(3)
      DATA NGROUP/'Cs ','C2 ','C2v','D2h','C2h','D2 ','Ci '/
      DATA NEL   /'YZ' ,'XZ' ,'XY' ,'X ' ,'Y ' ,'Z ' ,'I ' /
      DATA ZERO  /0.0D0/
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** COMPUTE CENTER-OF-GRAVITY COORDINATES CG(3)
C     WHICH ARE USED IN SUBROUTINE ASSEL.
      AM     = ZERO
      CG(1)  = ZERO
      CG(2)  = ZERO
      CG(3)  = ZERO
      DO 10 I=1,NUMAT
      AM     = AM+AMS(I)
      CG(1)  = CG(1)+AMS(I)*COORD(1,I)
      CG(2)  = CG(2)+AMS(I)*COORD(2,I)
      CG(3)  = CG(3)+AMS(I)*COORD(3,I)
   10 CONTINUE
      CG(1)  = CG(1)/AM
      CG(2)  = CG(2)/AM
      CG(3)  = CG(3)/AM
C *** ASSIGN POINT GROUP BY LOOPING OVER SYMMETRY ELEMENTS.
      ISUB   = 0
      MSYM   = 0
      DO 30 I=1,7
      MSYM1  = MSYM+1
      CALL ASSEL(I,NUMAT,AMS,COORD,CG,ISYM,MSYM1,IERROR)
      IF(IERROR.GT.0) GO TO 30
      MSYM   = MSYM1
      IEL(MSYM) = I
      IF(MSYM.EQ.3) GO TO 40
   30 CONTINUE
   40 CONTINUE
      MSYM1  = MSYM
      IF(MSYM.GT.0) LEL=IEL(MSYM)
      IF(MSYM.EQ.1) THEN
         ISUB = 1
         IF(LEL.GT.3) ISUB=2
         IF(LEL.EQ.7) ISUB=7
      ELSE IF(MSYM.EQ.2) THEN
         MSYM = 0
      ELSE IF(MSYM.EQ.3) THEN
         ISUB = 3
         IF(LEL.EQ.3) ISUB=4
         IF(LEL.EQ.7) ISUB=5
         IF(IEL(1).EQ.4) ISUB=6
      ENDIF
C *** PRINTING SECTION.
      IF(JPRINT.LT.0) RETURN
      IF(MSYM.EQ.0) THEN
         IF(JPRINT.GT.0) WRITE(NB6,500)
         IF(MSYM1.EQ.2)  WRITE(NB6,510)
      ELSE
         WRITE(NB6,520) NGROUP(ISUB)
         IF(JPRINT.GT.0) THEN
            WRITE(NB6,530)
            DO 50 I=1,MSYM
            M   = IEL(I)
            WRITE(NB6,540) NEL(M),(ISYM(I,J),J=1,NUMAT)
   50       CONTINUE
         ENDIF
      ENDIF
      RETURN
  500 FORMAT (//5X,'POINT GROUP C1. NO SYMMETRY.')
  510 FORMAT (/ 5X,'INCONSISTENT ASSIGNMENT OF POINT GROUP.'
     1        / 5X,'CALCULATION DONE WITHOUT SYMMETRY.')
  520 FORMAT (//5X,'POINT GROUP ',A3,' ASSIGNED.')
  530 FORMAT (//5X,'SYMMETRY ELEMENT AND SYMMETRY RELATIONS BETWEEN ',
     1             'THE ATOMS.',
     2        / 5X,'SYMMETRY OPERATIONS PERMUTE THE ATOMS (1..N) ',
     3             'AS FOLLOWS.'/)
  540 FORMAT (  5X,A2,2X,30I4,/(9X,30I4))
      END