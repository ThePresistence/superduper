      SUBROUTINE SYMNUM (AMS,PC,ITYPE,N,ISYM,NUMSYM,LPRINT)
C     *
C     SYMMETRY NUMBER FOR ROTATIONAL PARTITION FUNCTION.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 NM1(2),NM3(40),NM4(8)
      COMMON /NBFILE/ NBF(20)
      DIMENSION AMS(N),PC(3,N),ISYM(3,N)
      DIMENSION IEL(3)
      DATA NM1/'C0v','D0h'/
      DATA NM3/'C3 ','C4 ','C5 ','C6 ','C7 ','C8 ',
     1         'D3 ','D4 ','D5 ','D6 ','D7 ','D8 ',
     2         'C3v','C4v','C5v','C6v','C7v','C8v',
     3         'C3h','C4h','C5h','C6h','C7h','C8h',
     4         'D3h','D4h','D5h','D6h','D7h','D8h',
     5         'D2d','D3d','D4d','D5d','D6d','D7d',
     6         'D8d','S4 ','S6 ','S8 '/
      DATA NM4/'C1 ','Cs ','C2 ','C2v','D2h','C2h','D2 ','Ci '/
      DATA PI/3.141592653589793D0/
C     AN ATOM IS CONSIDERED TO LIE ON A SYMMETRY AXIS IF THE RELEVANT
C     COORDINATES DIFFER BY LESS THAN SMALL.
      DATA SMALL/1.0D-07/
C     TWO ATOMIC MASSES ARE CONSIDERED TO BE EQUAL IF THEY DIFFER BY
C     LESS THAN SMALLM.
      DATA SMALLM/1.0D-07/
C *** FILE NUMBERS.
      NB6    = NBF(6)
C     *
C     INITIALIZATION.
      IF(NUMSYM.GT.0) THEN
         IF(LPRINT.GE.0) WRITE(NB6,580) NUMSYM
         RETURN
      ENDIF
      NUMSYM = 1
C     *
C     DETERMINATION OF ABELIAN POINT GROUP.
      CALL ASSYM (N,AMS,PC,ISUB,MSYM,IEL,ISYM,-1)
C     *
C     LINEAR MOLECULE.
      IF(ITYPE.EQ.1) THEN
         IF(ISUB.EQ.4) NUMSYM=2
         GO TO 100
      ENDIF
C     *
C     SPHERICAL TOP.
      IF(ITYPE.EQ.2) THEN
         NUMSYM = 12
         IF(ISUB.EQ.4) NUMSYM=24
         GO TO 100
      ENDIF
C     *
C     ASYMETRIC TOP.
      IF(ITYPE.EQ.4) THEN
         IF(ISUB.EQ.2 .OR. ISUB.EQ.3 .OR. ISUB.EQ.5) NUMSYM=2
         IF(ISUB.EQ.4 .OR. ISUB.EQ.6) NUMSYM=4
         GO TO 100
      ENDIF
C     *
C     SYMMETRIC TOP.
C     BY CONVENTION THE PRINCIPAL-AXIS COORDINATES HAVE THE Z-AXIS
C     AS HIGHEST-ORDER SYMMETRY AXIS. WE FIRST DETERMINE THE ORDER
C     OF THIS AXIS WHICH MUST BE BETWEEN 2 AND 8.
      DO 30 I=8,2,-1
      PH     = 2.0D0*PI/I
      COSPH  = COS(PH)
      SINPH  = SIN(PH)
      DO 20 J=1,N
      T1     = PC(1,J)
      T2     = PC(2,J)
      IF(ABS(T1).LT.SMALL .AND. ABS(T2).LT.SMALL) GO TO 20
      TX     = COSPH*T1-SINPH*T2
      TY     = SINPH*T1+COSPH*T2
      TZ     = PC(3,J)
      AMSJ   = AMS(J)
      DO 10 K=1,N
      IF(ABS(AMS(K)-AMSJ).GT.SMALLM) GO TO 10
      IF(ABS(PC(1,K)-TX).GT.SMALL) GO TO 10
      IF(ABS(PC(2,K)-TY).GT.SMALL) GO TO 10
      IF(ABS(PC(3,K)-TZ).GT.SMALL) GO TO 10
      GO TO 20
   10 CONTINUE
      GO TO 30
   20 CONTINUE
      IORDER = I
      GO TO 40
   30 CONTINUE
      IF(LPRINT.GT.-5) WRITE(NB6,500)
      RETURN
C     CHECK WHETHER IORDER IS EVEN OR ODD.
   40 IODD   = IORDER-2*(IORDER/2)
C     CLASSIFY DEGENERATE POINT GROUPS BY IPG. NOTATION AS FOLLOWS.
C     1-6 CN, 7-12 DN, 13-18 CNV, 19-24 CNH, 25-30 DNH, 31-37 DND,
C     38-40 SN. FOR MORE DETAIL, SEE DATA STATEMENT FOR ARRAY NM3.
C     CLASSIFICATION BASED ON THE ABELIAN SUBGROUP (ISUB) AND THE
C     ORDER OF THE HIGHEST SYMMETRY AXIS (IORDER,IODD).
      IF(ISUB.EQ.0) THEN
         IPG = IORDER-2
         GO TO 50
      ENDIF
      IF(ISUB.EQ.1) THEN
         IF(IODD.EQ.1) IPG=IORDER+10
         IF(IODD.EQ.0) IPG=IORDER+16
         GO TO 50
      ENDIF
      IF(ISUB.EQ.2) THEN
         IF(IODD.EQ.0) IPG=IORDER-2
         IF(IODD.EQ.1) IPG=IORDER+4
         IF(IORDER.EQ.2) IPG=38
         GO TO 50
      ENDIF
      IF(ISUB.EQ.3) THEN
         IF(IODD.EQ.0) IPG=IORDER+10
         IF(IODD.EQ.1) IPG=IORDER+22
         IF(IORDER.EQ.2) IPG=31
         GO TO 50
      ENDIF
      IF(ISUB.EQ.4) THEN
         IPG = IORDER+22
         GO TO 50
      ENDIF
      IF(ISUB.EQ.5) THEN
         IF(IODD.EQ.0) IPG=IORDER+16
         IF(IODD.EQ.1) IPG=IORDER+29
         GO TO 50
      ENDIF
      IF(ISUB.EQ.6) THEN
         IPG = IORDER+4
         GO TO 50
      ENDIF
      IF(ISUB.EQ.7) IPG=39
C     DEFINE SYMMETRY NUMBER FOR SYMMETRIC TOPS.
   50 ICD    = 1
      IF(IPG.GT. 6 .AND. IPG.LT.13) ICD=2
      IF(IPG.GT.24 .AND. IPG.LT.38) ICD=2
      NUMSYM = IORDER*ICD
C     *
C     PRINTING SECTION.
  100 IF(LPRINT.LT.0) RETURN
      WRITE(NB6,510) NUMSYM
      IF(ITYPE.EQ.1) THEN
         WRITE(NB6,520) NM1(NUMSYM)
         RETURN
      ENDIF
      IF(ITYPE.EQ.2) THEN
         IF(NUMSYM.EQ.12) WRITE(NB6,530)
         IF(NUMSYM.EQ.24) WRITE(NB6,535)
         RETURN
      ENDIF
      IF(ITYPE.EQ.4) THEN
         WRITE(NB6,540) NM4(ISUB+1)
         RETURN
      ENDIF
      WRITE(NB6,540) NM3(IPG)
      IF(IPG.LT.7) THEN
         WRITE(NB6,550)
         WRITE(NB6,570)
      ENDIF
      IF(IPG.EQ.14 .OR. IPG.EQ.16 .OR. IPG.EQ.18) THEN
         WRITE(NB6,560)
         WRITE(NB6,570)
      ENDIF
      RETURN
  500 FORMAT(///5X,'CALCULATION OF SYMMETRY NUMBER FOR SYMMETRIC TOP.',
     1       /  5X,'THE ORDER OF THE Z-AXIS IS NOT BETWEEN 2 AND 8.',
     2       /  5X,'UNLESS THE ORDER IS HIGHER THAN 8, THERE MUST ',
     3       /  5X,'BE SOME INCONSISTENCY. THE SYMMETRY NUMBER IS ',
     4       /  5X,'ASSUMED TO BE 1 IN THE FOLLOWING CALCULATION.')
  510 FORMAT(///5X,'SYMMETRY NUMBER',I3,' ASSIGNED ',
     1             'FOR ROTATIONAL PARTITION FUNCTION.')
  520 FORMAT(/  5X,'LINEAR MOLECULE WITH POINT GROUP ',A3)
  530 FORMAT(/  5X,'TETRAHEDRAL POINT GROUP ASSIGNED.')
  535 FORMAT(/  5X,'OCTAHEDRAL POINT GROUP ASSIGNED.')
  540 FORMAT(/  5X,'POINT GROUP    ',A3,' ASSIGNED ',
     1             'USING PRINCIPAL-AXIS COORDINATES.')
  550 FORMAT(/  5X,'THE PRESENT PROGRAM DOES NOT DISTINGUISH BETWEEN ',
     1             'CN AND DN POINT GROUPS, IN GENERAL.')
  560 FORMAT(/  5X,'THE PRESENT PROGRAM DOES NOT DISTINGUISH BETWEEN ',
     1             'CNV AND DND POINT GROUPS (N=4,6,8 EVEN).')
  570 FORMAT(   5X,'IF THE POINT GROUP IS ACTUALLY OF D-TYPE, THE ',
     1             'ASSIGNED SYMMETRY NUMBER IS TOO SMALL BY A ',
     2             'FACTOR OF 2.')
  580 FORMAT(// 5X,'SYMMETRY NUMBER',I3,' FROM INPUT DATA.')
      END
