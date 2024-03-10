      SUBROUTINE DYNDMS (LEN,LM5,LM5X,LM5R)
C     *
C     DYNAMIC MEMORY ALLOCATION FOR DENSITY MATRIX SEARCH.
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./INOPT2/ IN2(300)
     ./LIMITS/ LM2,LM3,LM4,LM6,LM7,LM8,LM9
     ./LMDMS / LD1,LD2,LD3,LD4,LD5
     ./NBFILE/ NBF(20)
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** INPUT VARIABLES.
      JPRINT = IN2(42)
      LINDMS = IN2(56)
C *** INITIALIZATION.
      LM5X   = 0
      LM5R   = 0
C     *
C *** SECTION FOR CONJUGATE GRADIENT DENSITY MATRIX SEARCH.
C     CONVENTIONAL CODE WITH FULL ARRAYS.
C     *
      IF(LINDMS.GT.0) THEN
         LENDMS = 5*LM2*LM3
         LM5X   = LM5+LENDMS
         IF(LM5X.GT.LEN) THEN
            WRITE(NB6,500) LEN,LM5X
            STOP 'DYNSCR'
         ENDIF
         LDL = LM2*LM3
         LD1 = LM5+1
         LD2 = LD1+LDL
         LD3 = LD2+LDL
         LD4 = LD3+LDL
         LD5 = LD4+LDL
         IF(JPRINT.GE.2) THEN
C           WRITE(NB6,510) LEN,LD1,LM5X
            WRITE(NB6,520) LD1,LD1+LDL-1,LDL,LD2,LD2+LDL-1,LDL,
     1                     LD3,LD3+LDL-1,LDL,LD4,LD4+LDL-1,LDL,
     2                     LD5,LD5+LDL-1,LDL
         ENDIF
      ENDIF
      RETURN
  500 FORMAT(///1X,'NOT ENOUGH MEMORY FOR IN-CORE CG-DMS PROCEDURE',
     1       /  1X,'MAXIMUM AVAILABLE MEMORY (WORDS)           ',I10,
     2       /  1X,'MEMORY REQUESTED (WORDS)                   ',I10,
     3       /  1X,'PROGRAM WILL STOP.')
C 510 FORMAT(///1X,'MEMORY ALLOCATION FOR IN-CORE CG-DMS PROCEDURE',
C    1       // 1X,'AVAILABLE LENGTH OF BLANK COMMON           ',I10,
C    2       /  1X,'ADDRESS OF FIRST AVAILABLE WORD            ',I10,
C    3       /  1X,'ADDRESS OF LAST  AVAILABLE WORD            ',I10)
  520 FORMAT(///1X,'MEMORY ALLOCATION FOR IN-CORE CG-DMS PROCEDURE',
     1       // 1X,'    START    FINAL   LENGTH    CONTENTS',
     2       // 1X,3I9,4X,'CG-DMS GRADIENT MATRIX',
     3       /  1X,3I9,4X,'CG-DMS SEARCH DIRECTION VECTOR',
     4       /  1X,3I9,4X,'FIRST SCRATCH ARRAY',
     5       /  1X,3I9,4X,'SECOND SCRATCH ARRAY',
     6       /  1X,3I9,4X,'THIRD SCRATCH ARRAY')
      END
