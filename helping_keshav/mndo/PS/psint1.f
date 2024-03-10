C     ******************************************************************
C
C     Store one-center two-electron integrals in matrix form.
C
C     ******************************************************************
      SUBROUTINE PSINT1(NTYP,NPAIR,AINTS)
C
C   Transform one-center two-electron integrals from
C   compact storagem to the full matrix form. 
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NTYP   - Atom type
C      NPAIR  - Number of electron pairs on atom
C      AINTS  - Output array. Integrals which do not
C               exist for one-center case (core integrals)
C               and derivatives are left unchanged.
C
C   Accessed common blocks:
C
C      REP    - One-center integrals in SP basis
C      DPARM4 - One-center integrals in D  basis
C
C   Local storage:
C
C      58 DOUBLE PRECISION words + 795 INTEGER words
C
C   Module logic:
C
C      Integral symmetry relations are mostly lifted from DPARM5
C      common, although unlike the main code I have a neat general
C      version which works with either SP or SPD basis...
C
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (AU  =3.67511944138184490995957D-2)
      PARAMETER (ZERO=0D0)
C
      COMMON 
     ./REP   / GSS(LMZ),GPP(LMZ),GSP(LMZ),GP2(LMZ),HSP(LMZ),HPP(LMZ)
     ./DPARM4/ REPD(52,LMZ)
C
      DIMENSION AINTS(46,46)
      DIMENSION UNIQUE(58), IJ(265), KL(265), II(265)
      DATA 
     .IJ/ 2, 2, 2, 2, 4, 7,11, 4, 7,11, 4, 4, 7, 7,11,11, 3, 5, 8,
     .    6, 9,10,16,46,22,37,29, 2, 2, 2, 2, 2,26, 8,13,39,18,40,
     .   32,20,33, 3, 3, 3, 5, 5, 8, 8,11,23,29,11,16,16,46,46,22,
     .   22,37,37, 4, 4, 4, 7, 7, 7,11,11,16,46,22,37, 4, 7,11,11,
     .   24,25, 3, 5, 4, 6, 9,10,12,38,17,30,28,36, 9,10,21,43,35,
     .   45, 6, 9, 9,10, 4, 7,23,23,29,29, 4, 7,27,44, 4, 6,14, 5,
     .   34,10, 7,12,27, 7,12,38,17,30,23,29,23,22,37,28,36,17,30,
     .   23,23,16,46,27,44,12,38,23,23,26,18,32,26,26,24,25,24,25,
     .   20,33,13,39,24,24,40,25,14,25,16,46,22,37,29,22,37,29,29,
     .   16,46,29,29,21,43,22,35,45,12,38,17,17,30,13,39,18,14,40,
     .   32,20,33,13,39,18,40,32,20,20,33,18,31,31,31,19,19,19,15,
     .   15,41,41,41,28,36,21,43,45,28,28,36,16,16,46,46,22,22,22,
     .   37,37,37,22,35,27,44,34,37,12,30,14,33,32,15,34,36,21,43,
     .   35,34,45,37,27,13,40,39,14,27,44,21,45,43,34,16,46,42/
      DATA
     .KL/ 2, 4, 7,11, 2, 2, 2, 4, 7,11, 7,11,11, 4, 4, 7, 3, 5, 8,
     .    6, 9,10, 2, 2, 2, 2, 2,16,46,22,37,29, 8,26, 3, 5, 8, 3,
     .    8, 3, 5,13,40,20,39,33,18,32,23,11,11,29, 4, 7, 4, 7, 4,
     .   11, 7,11,16,46,22,16,46,37,22,37,11,11, 7, 4,37,22,16,46,
     .    3, 5,24,25,12,38,17,30, 4, 6, 9,10, 9,10,28,36, 9,10, 6,
     .    9,35,21,45,43,23,23, 4, 7, 4, 7,29,29, 4, 6,27,44, 5,14,
     .   10,34,12, 7, 7,27,12,38,17,30,23,23,29,23,23,17,30,28,36,
     .   22,37,23,23,12,38,27,44,16,46,26,26,26,18,32,24,25,20,33,
     .   24,25,24,25,13,40,24,39,25,14,16,46,22,37,29,29,29,22,37,
     .   29,29,16,46,17,30,12,38,17,22,35,21,45,43,13,39,18,14,40,
     .   32,20,33,20,33,32,20,18,13,40,39,15,31,19,41,31,19,41,18,
     .   15,31,19,41,28,36,28,36,28,21,45,43,22,37,22,37,16,46,37,
     .   16,46,22,27,44,22,35,30,12,37,34,33,14,15,32,36,34,21,43,
     .   35,34,45,27,37,40,13,14,39,27,44,45,21,34,43,46,16,42/
      DATA
     .II/ 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5,
     .    6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 9, 9, 9, 9,
     .    9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,11,11,12,12,12,12,12,
     .   12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,
     .   14,14,14,14,15,15,15,15,15,15,15,15,16,16,16,16,17,17,17,
     .   17,17,17,17,17,18,18,18,18,19,19,19,19,20,20,20,20,21,21,
     .   22,22,23,23,24,24,25,25,25,25,25,26,26,27,27,27,27,27,27,
     .   27,27,28,28,28,28,28,28,28,28,29,30,30,30,30,31,31,32,32,
     .   32,32,33,33,33,33,33,33,34,34,35,35,35,35,35,36,36,36,36,
     .   37,37,37,37,38,38,38,38,38,38,38,38,38,38,39,39,39,39,39,
     .   39,39,39,40,40,40,40,40,40,40,40,41,41,41,41,41,41,41,41,
     .   41,41,41,41,42,42,43,43,43,43,43,43,44,44,44,44,44,44,44,
     .   44,44,44,45,45,45,45,46,46,46,46,47,47,48,48,49,49,50,50,
     .   50,50,50,51,51,52,52,53,53,54,54,55,55,56,56,57,57,58/
C
      DO 100 J=1,NPAIR
          DO 110 I=1,NPAIR
              AINTS(1+I,1+J) = ZERO
  110     CONTINUE
  100 CONTINUE
C
C    Move integrals to the local array, so that we could treat
C    SP and SPD cases uniformly.
C
      UNIQUE(1) = AU * GSS(NTYP)
      IF( NPAIR.GT.1 ) THEN
          UNIQUE(2) = AU * GSP(NTYP)
          UNIQUE(3) = AU * GPP(NTYP)
          UNIQUE(4) = AU * GP2(NTYP)
          UNIQUE(5) = AU * HSP(NTYP)
          UNIQUE(6) = AU * HPP(NTYP)
          IF( NPAIR.GT.10 ) THEN
              DO 200 I=1,52
                  UNIQUE(6+I) = AU * REPD(I,NTYP)
  200         CONTINUE
              NZERO = 265
          ELSE
              NZERO = 22
          ENDIF
      ELSE
          NZERO = 1
      ENDIF
C
C    Now copy integrals to their destination
C
      DO 300 I=1,NZERO
          AINTS(IJ(I),KL(I)) = UNIQUE(II(I))
  300 CONTINUE
      RETURN
      END
