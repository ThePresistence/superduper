      SUBROUTINE DYNVB (LEN,LM5,LM5X,LM5R)
C     *
C     DYNAMIC MEMORY ALLOCATION FOR VALENCE BOND INTERFACE.
C     *
C     DEFINE BUFFER FOR THE OVERLAP MATRIX UNLESS THIS HAS ALREADY
C     BEEN DONE BEFORE (E.G. FOR OM2-TYPE METHODS).
C     *
C     NOTATION.
C     LEN  = AVAILABLE LENGTH OF BLANK COMMON AREA (INPUT VALUE,
C            NOT MODIFIED IN THIS ROUTINE).
C     LM5  = LAST ADDRESS OF BLANK COMMON AREA THAT HAS PREVIOUSLY
C            BEEN RESERVED (INPUT VALUE, NOT MODIFIED IN THIS ROUTINE).
C     LM5X = LAST ADDRESS OF BLANK COMMON AREA THAT IS ACTUALLY USED
C            BY THE ARRAYS DEFINED PRESENTLY (DETERMINED HERE).
C     LM5R = LAST ADDRESS OF BLANK COMMON AREA THAT NEEDS TO REMAIN
C            RESERVED LATER ON (DETERMINED IN THIS ROUTINE).
C     *
C     THE ARRAYS ALLOCATED IN THIS ROUTINE START AT ADDRESS LM5+1
C     AND EXTEND TO ADDRESS LM5X.
C     *
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OM2
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./INOPT2/ IN2(300)
     ./LIMITS/ LM2,LM3,LM4,LM6,LM7,LM8,LM9
     ./LMGRAD/ LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8
     ./NBFILE/ NBF(20)
C     *
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** INPUT VARIABLES.
      IOP    = IN2(2)
      JPRINT = IN2(42)
      IVBSE  = IN2(213)
C *** INITIALIZATION.
      LM5X   = 0
      LM5R   = 0
      OM2    = (IOP.EQ.-6 .OR. IOP.EQ.-8 .OR. IOP.EQ.-9
     1          .OR. IOP.EQ.-22 .OR. IOP.EQ.-23) .AND. NUMAT.GT.2
      IF(OM2 .OR. IVBSE.EQ.0) RETURN
C *** MEMORY ALLOCATION.
      LM5X   = LM5+LM4+LM2*LM2
      IF(LM5X.LE.LEN) THEN
         LG1 = LM5+1
         LG2 = LG1+LM4
      ELSE
         WRITE(NB6,600) LM5X-LM5,LEN-LM5
         STOP 'DYNVB'
      ENDIF
C *** PRINTING SECTION.
      IF(JPRINT.GE.2) THEN
         WRITE(NB6,610) LG1,LG1+LM4-1,LM4
         WRITE(NB6,620) LG2,LG2+LM2*LM2-1,LM2*LM2
      ENDIF
      RETURN
  600 FORMAT(///1X,'BUFFER REQUESTED FOR VALENCE BOND INTERFACE',I10,
     1       /  1X,'MAXIMUM AVAILABLE BUFFER (WORDS)           ',I10,
     2       /  1X,'PROGRAM WILL STOP.'//)
  610 FORMAT(///1X,'MEMORY ALLOCATION FOR VALENCE BOND INTERFACE',
     1       // 1X,'    START    FINAL   LENGTH    CONTENTS',
     2       // 1X,3I9,4X,'OVERLAP MATRIX')
  620 FORMAT(   1X,3I9,4X,'BUFFER FOR TRANSFORMATION')
      END
