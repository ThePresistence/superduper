      SUBROUTINE PDBINP (NB5,NB6,JPRINT)
C     *
C     READ A PDB-COMPATIBLE INPUT FILE.
C     *
C     NOTATION. I=INPUT.
C     NB5       FILE NUMBER FOR PDB INPUT (I).
C     NB6       FILE NUMBER FOR STANDARD OUTPUT (I).
C     JPRINT    PRINTING FLAG (I).
C     *
C     PDB INPUT DATA ARE SAVED IN COMMON BLOCKS:
C     ATOMC, ATOMS, FRGMT1.
C     OTHER RELEVANT COMMON BLOCKS ARE ALSO FILLED:
C     DFP, PARM1, PARM2, PARM3, PARM4, PARM5.
C     *
      USE LIMIT, ONLY: LM1, LMG, LMR, LMS, LMV, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*6 RECORD
      CHARACTER*4 NAME
      CHARACTER*3 RES
      CHARACTER*2 ELEM
      CHARACTER*2 ELEMNT
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./DFP   / XX(LMV),NVAR
     ./ELEMTS/ ELEMNT(107)
     ./FRGMT1/ NFRAGS(LM1),NCHRGS(LM1)
     ./INOPT2/ IN2(300)
     ./PARM1 / A(3,LM1),NC(LM1),NB(LM1),NA(LM1),NN(LM1),NATOMS
     ./PARM2 / NSYM,LPAR(LMS),LNUM(LMS),LDEP(LMS)
     ./PARM3 / LOC(LMV),NV
     ./PARM4 / RC(LMR),LREACT,LTOTAL
     ./PARM5 / RC1(LMG),RC2(LMG),LGRID1,LTOT1,LGRID2,LTOT2
C     *
C *** INITIALIZATION.
      NUM    = 0
      IFAIL  = 0
      NCHRG  = 0
      NRES   = 0
C *** READ PDB FILE (NB5).
   10 CONTINUE
      READ(NB5,500,ERR=30,END=100) RECORD,I,NAME,RES,IRES,X,Y,Z,ELEM
      IF(RECORD(1:3).EQ.'END') GO TO 100
      IF(RECORD.NE.'ATOM  ' .AND. RECORD.NE.'HETATM') GO TO 10
C     ECHO INPUT LINES FOR DEBUG PRINT.
C     WRITE(NB6,505) RECORD,I,NAME,RES,IRES,X,Y,Z,ELEM
C *** IDENTIFY ELEMENT.
      IF(ELEM(1:1).EQ.' ') THEN
         ELEM(1:1) = ELEM(2:2)
         ELEM(2:2) = ' '
      ENDIF
      DO 20 N=1,LMZ
      IF(ELEM.EQ.ELEMNT(N)) THEN
         NUM = NUM + 1
         NNN = N
         GO TO 40
      ENDIF
   20 CONTINUE
C *** ERROR WHEN READING INPUT LINE.
C     ENCOUNTERING END-OF-RECORD IN LAST INPUT LINE IS NO ERROR.
   30 CONTINUE
      IFAIL  = IFAIL+1
      GO TO 10
C *** STORE ATOMIC NUMBER, CARTESIAN COORDINATES, AND RESIDUE NUMBER.
   40 CONTINUE
      NAT(NUM) = NNN
      COORD(1,NUM) = X
      COORD(2,NUM) = Y
      COORD(3,NUM) = Z
      NFRAGS(NUM)  = IRES
C *** ASSIGN FORMAL CHARGE OF RESIDUE (NCHRGS).
C     STORE ITS VALUE AT A LOCATION DEFINED BY ITS FIRST ATOM.
C     CHECK FOR FIRST ATOM (NUM1) OF RESIDUE.
      NCHRGS(NUM)  = 0
      IF(NUM.EQ.1) THEN
         NUM1 = 1
      ELSE IF(IRES.NE.NFRAGS(NUM-1)) THEN
         NUM1 = NUM
      ENDIF
C     ASSIGN POSITIVE CHARGE FOR ARGININE AND LYSINE.
C     ASSIGN NEGATIVE CHARGE FOR GLUTAMATE AND ASPARTATE.
      IF(NUM1.EQ.NUM) THEN
         IF(RES.EQ.'ARG' .OR. RES.EQ.'LYS') THEN
            NCHRGS(NUM) = 1
         ELSE IF(RES.EQ.'GLU' .OR. RES.EQ.'ASP') THEN
            NCHRGS(NUM) =-1
         ENDIF
      ENDIF
C     UPDATE TOTAL CHARGE AND MAXIMUM RESIDUE NUMBER.
      NCHRG  = NCHRG+NCHRGS(NUM)
      NRES   = MAX(IRES,NRES)
      GO TO 10
C *** INPUT LOOP FINISHED.
C     PDB CONVENTION: FIRST RESIDUE - AMINO TERMINUS (N).
C     PDB CONVENTION: LAST RESIDUE - CARBOXYL TERMINUS (C).
C     OUR CONVENTION: N TERMINUS GETS FORMALLY A POSITIVE CHARGE.
C     OUR CONVENTION: C TERMINUS GETS FORMALLY A NEGATIVE CHARGE.
  100 CONTINUE
      NCHRGS(1)    = NCHRGS(1)   +1
      NCHRGS(NUM1) = NCHRGS(NUM1)-1
C     *
C *** SAVE RELEVANT INFORMATION.
      NATOMS = NUM
      NUMAT  = NUM
      NSYM   = 0
      LTOTAL = 0
      LTOT1  = 0
      LTOT2  = 0
      NV     = 0
      DO 120 N=1,NUMAT
      NA(N)  = 0
      NB(N)  = 0
      NC(N)  = 0
      NN(N)  = NAT(N)
      DO 110 J=1,3
      NV     = NV+1
      LOC(NV)= 3*(N-1)+J
      XX(NV) = COORD(J,N)
      A(J,N) = COORD(J,N)
  110 CONTINUE
  120 CONTINUE
      NVAR   = NV
C     SET FLAG FOR CARTESIAN COORDINATES.
      IN2(4) = 1
C     *
C *** PRINTING SECTION.
      IF(JPRINT.GE.-5) THEN
         WRITE(NB6,510) NUM
         IF(IFAIL.GT.0 .AND. JPRINT.GE.5) WRITE(NB6,520) IFAIL
         WRITE(NB6,530) NRES
         WRITE(NB6,540) NCHRG,IN2(65)
      ENDIF
      IF(JPRINT.GE.2) THEN
         WRITE(NB6,550)
         DO 130 N=1,NUMAT
         WRITE(NB6,560) N,NAT(N),COORD(1,N),COORD(2,N),COORD(3,N),
     1                  NFRAGS(N),NCHRGS(N)
  130    CONTINUE
         WRITE(NB6,570)
      ENDIF
      RETURN
C     SOME UNNEEDED ENTRIES IN THE PDB RECORDS ARE SKIPPED.
C     FULL FORMAT OF PDB ATOM RECORDS IS GIVEN FOR DOCUMENTATION.
C 500 FORMAT(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,6X,A4,A2,A2)
  500 FORMAT(A6,I5,1X,A4,1X,A3,1X,1X,I4,1X,3X,3F8.3,6X,6X,6X,4X,A2,2X)
C 505 FORMAT(1X,A6,I5,1X,A4,1X,A3,1X,1X,I4,1X,3X,3F8.3,6X,6X,6X,4X,A2)
  510 FORMAT(// 1X,'INPUT FROM PDB ATOMS SECTION.',
     1       /  1X,'NUMBER OF LINES READ SUCCESSFULLY:',I5)
  520 FORMAT(   1X,'NUMBER OF LINES WITH READ FAILURE:',I5)
  530 FORMAT(   1X,'LARGEST PDB RESIDUE NUMBER       :',I5)
  540 FORMAT(   1X,'TOTAL CHARGE ASSIGNED            :',I5,
     1       /  1X,'MOLECULAR CHARGE FROM INPUT      :',I5)
  550 FORMAT(///1X,'INPUT DATA READ FROM PDB FILE.',
     1       // 1X,'    I    NAT      X         Y         Z   ',
     2             '  RESIDUE   CHARGE'/)
  560 FORMAT(   1X,I5,I6,3F10.3,I7,I10)
  570 FORMAT(   1X)
      END
