      SUBROUTINE DYNOM2 (LEN,LM5,LM5X,LM5R)
C     *
C     DYNAMIC MEMORY ALLOCATION FOR ORTHOGONALIZATION CORRECTIONS.
C     RELEVANT FOR OM2 AND OM3.
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
      INOUT  = IN2(212)
C *** INITIALIZATION.
      LM5X   = 0
      LM5R   = 0
      LSCMAX = LEN-LM5
      OM2    = (IOP.EQ.-6 .OR. IOP.EQ.-8 .OR. IOP.EQ.-9
     1          .OR.IOP.EQ.-22 .OR.IOP.EQ.-23) .AND. NUMAT.GT.2
C     *
C     THE DYNAMIC MEMORY ALLOCATION FOR THE SCF TREATMENT (IN DYNSCF)
C     DEFINES THE FOLLOWING ORDER OF ARRAYS (IN STANDARD NOTATION):
C     ARRAY: C(L22),E(LM2),W(LM9),F(L22),H(LM4),D(LM4),P(LM4),Q(LD)
C     START: LS1....LS2....LS9....LS6....LS7....LS5....LS8....LS3
C     *
C     THE DYNAMIC MEMORY ALLOCATION FOR ORTHOGONALIZATION CORRECTIONS
C     EMPLOYS THE FOLLOWING ARRAYS:
C     ARRAY: S(LM4),B(LM4),COR(),HG(LM4),SG(LM4),BG(LM4),CORG(),HREF().
C     START: LG1....LG2....LG3...LG4.....LG5.....LG6.....LG7....LG8
C     THE FIRST FOUR ARRAYS ARE USED FOR THE ENERGY   EVALUATION.
C     THE LAST  FIVE ARRAYS ARE USED FOR THE GRADIENT EVALUATION.
C     THE ARRAYS START AT ADDRESS LM5+1.
C     *
C     NOTE THAT THE FIRST FOUR ARRAYS ARE COMPUTED IN THE INTEGRAL
C     SECTION AT THE REFERENCE GEOMETRY, BEFORE THE SCF ITERATIONS.
C     THEY ARE USED AS INPUT FOR THE GRADIENT EVALUATION, AFTER THE
C     SCF ITERATIONS. HENCE, IF THE SCF SECTION MAKES USE OF THE
C     GENERAL SCRATCH SPACE (E.G. IN DIIS), THE FIRST FOUR ARRAYS
C     MUST REMAIN UNCHANGED (E.G. BY USING ONLY SCRATCH SPACE BEYOND
C     ADDRESS LG5). ALTERNATIVELY, THE FIRST FOUR ARRAYS MAY BE SAVED
C     ON FILE IN THE INTEGRAL SECTION AND READ FROM FILE AGAIN IN THE
C     GRADIENT SECTION (INPUT OPTION INOUT.GT.0), OR THEY MAY BE
C     RECOMPUTED IN THE GRADIENT SECTION (INPUT OPTION INOUT.LT.0).
C     *
      IF(OM2) THEN
         LMC    = LM2*NUMAT
         LMD    = NUMAT*(NUMAT+1)/2
         LENORT = 5*LM4+2*LMC+LMD
         IF(LENORT.LE.LSCMAX) THEN
            LG1 = LM5+1
         ELSE
            WRITE(NB6,600) LENORT,LSCMAX
            STOP 'DYNSCR'
         ENDIF
         LG2    = LG1+LM4
         LG3    = LG2+LM4
         LG4    = LG3+LMC
         LG5    = LG4+LM4
         LG6    = LG5+LM4
         LG7    = LG6+LM4
         LG8    = LG7+LMC
         LM5X   = LM5+LENORT
         IF(INOUT.EQ.0) THEN
           LM5R = LG5-1
         ENDIF
         IF(JPRINT.GE.2) THEN
C           WRITE(NB6,610) LEN,LG1,LM5X,LM5R
            WRITE(NB6,620) LG1,LG1+LM4-1,LM4,LG2,LG2+LM4-1,LM4,
     1                     LG3,LG3+LMC-1,LMC,LG4,LG4+LM4-1,LM4,
     2                     LG5,LG5+LM4-1,LM4,LG6,LG6+LM4-1,LM4,
     3                     LG7,LG7+LMC-1,LMC,LG8,LG8+LMD-1,LMD
         ENDIF
      ENDIF
      RETURN
  600 FORMAT(///1X,'BUFFER REQUESTED FOR ORTHOGONALIZATION     ',I10,
     1       /  1X,'MAXIMUM AVAILABLE BUFFER (WORDS)           ',I10,
     2       /  1X,'PROGRAM WILL STOP.'//)
C 610 FORMAT(///1X,'MEMORY ALLOCATION FOR ORTHOGONALIZATION (OM2)',
C    1       // 1X,'AVAILABLE LENGTH OF BLANK COMMON           ',I10,
C    2       /  1X,'ADDRESS OF FIRST AVAILABLE WORD            ',I10,
C    3       /  1X,'ADDRESS OF LAST  AVAILABLE WORD            ',I10,
C    4       /  1X,'ADDRESS OF LAST  RESERVED  WORD            ',I10)
  620 FORMAT(///1X,'MEMORY ALLOCATION FOR ORTHOGONALIZATION (OM2)',
     1       // 1X,'    START    FINAL   LENGTH    CONTENTS',
     2       // 1X,3I9,4X,'OVERLAP MATRIX, REFERENCE',
     3       /  1X,3I9,4X,'RESONANCE MATRIX, REFERENCE',
     4       /  1X,3I9,4X,'LOCAL CORE MATRIX, REFERENCE',
     5       /  1X,3I9,4X,'ORTHOGONALIZATION CORRECTIONS',
     6       /  1X,3I9,4X,'OVERLAP MATRIX, GRADIENT',
     7       /  1X,3I9,4X,'RESONANCE MATRIX, GRADIENT',
     8       /  1X,3I9,4X,'LOCAL CORE MATRIX, GRADIENT',
     9       /  1X,3I9,4X,'ATOM-PAIR BASED THRESHOLDS, GRADIENT')
      END
