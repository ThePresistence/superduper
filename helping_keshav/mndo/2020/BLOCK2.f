      BLOCKDATA BLOCK2
C     *
C     INITIALIZATION FOR D ORBITALS.
C     *
C     COMMON BLOCKS: ATORB, DELEMT, DNBND, DPARM5.
C     *
C     SPECIAL CONVENTION.
C     ELEMENT 86 IS A CONNECTION ATOM FOR QM/MM TREATMENTS.
C     CORE CHARGE 1, ONE ELECTRON IN A 2S ORBITAL.
C     NON-STANDARD VALUES FOR IOS(86),IOP(86),III(86),IIID(86).
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATORB / IOS(LMZ),IOP(LMZ),IOD(LMZ)
     ./DELEMT/ NELMD,NFOCK
     ./DNBND / III(LMZ),IIID(LMZ)
     ./DPARM5/ INTIJ(243),INTKL(243),INTREP(243),INTRF1(243),INTRF2(243)
C     GROUND-STATE OCCUPATION NUMBERS IOS,IOP,IOD FOR S,P,D ORBITALS.
      DATA IOS/1,2,1,7*2,1,7*2,1,2, 2,2,2,1,2,2,2,2,1,2,6*2,1,3*2,1,1,
     1         2,1,1,0,1,7*2,1,22*2,1,1,6*2,1/
C    1         2,1,1,0,1,7*2,1,22*2,1,1,7*2/
      DATA IOP/4*0,1,2,3,4,5,6, 0,0,1,2,3,4,5,6,12*0,1,2,3,4,5,6,12*0,
     1             1,2,3,4,5,6,26*0,1,2,3,4,5,0/
C    1             1,2,3,4,5,6,26*0,1,2,3,4,5,6/
      DATA IOD/20*0,1,2,3,5,5,6,7,8,10, 0,6*0,2*0,1,2,4,5,5,7,8,3*10,
     1          6*0,0,0,1,14*0,2,3,4,5,6,7,9,10,0,6*0/
C     NELMD  -  NUMBER OF ATOMS WITH D ORBITALS.
C     NFOCK  -  STEP SIZE IN SUBROUTINE FOCK (5 FOR SP, 10 FOR SPD).
      DATA NELMD,NFOCK /0,5/
C     III -ARRAY -  PRINCIPAL QUANTUM NUMBERS OF SP-AO
C     IIID-ARRAY -  PRINCIPAL QUANTUM NUMBERS OF  D-AO
      DATA III  /2*1,8*2,8*3,18*4,18*5,31*6,2/
C     DATA III  /2*1,8*2,8*3,18*4,18*5,32*6/
      DATA IIID /30*3,18*4,32*5,5*6,3/
C     DATA IIID /30*3,18*4,32*5,6*6/
C     INDICES FOR ONE-CENTER TWO-ELECTRON INTEGRALS.
C     W(INTIJ(I),INTKL(I)) = REP(INTREP(I))
C     IN THE RAFFENETTI SCHEME, THE INTEGRALS INCLUDE AN EXTRA TERM
C     -0.25*(REP(INTRF1)+REP(INTRF2)).
      DATA INTIJ/
     1    1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4,
     2    4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 9,
     3    9, 9, 9,10,10,10,10,10,10,11,11,11,11,11,11,12,12,12,12,12,
     4   13,13,13,13,13,14,14,14,15,15,15,15,15,15,15,15,15,15,16,16,
     5   16,16,16,17,17,17,17,17,18,18,18,19,19,19,19,19,20,20,20,20,
     6   20,21,21,21,21,21,21,21,21,21,21,21,21,22,22,22,22,22,22,22,
     7   22,22,23,23,23,23,23,24,24,24,24,24,25,25,25,25,26,26,26,26,
     8   26,26,27,27,27,27,27,28,28,28,28,28,28,28,28,28,28,29,29,29,
     9   29,29,30,30,30,31,31,31,31,31,32,32,32,32,32,33,33,33,33,33,
     A   34,34,34,34,35,35,35,35,35,36,36,36,36,36,36,36,36,36,36,36,
     B   36,37,37,37,37,38,38,38,38,38,39,39,39,39,39,40,40,40,41,42,
     C   42,42,42,42,43,43,43,43,44,44,44,44,44,45,45,45,45,45,45,45,
     D   45,45,45/
      DATA INTKL/
     1   15,21,28,36,45,12,19,23,39,11,15,21,22,26,28,36,45,13,24,32,
     2   38,34,37,43,11,15,21,22,26,28,36,45,17,25,31,16,20,27,44,29,
     3   33,35,42,15,21,22,28,36,45, 3, 6,11,21,26,36, 2,12,19,23,39,
     4    4,13,24,32,38,14,17,31, 1, 3, 6,10,15,21,22,28,36,45, 8,16,
     5   20,27,44, 7,14,17,25,31,18,30,40, 2,12,19,23,39, 8,16,20,27,
     6   44, 1, 3, 6,10,11,15,21,22,26,28,36,45, 3, 6,10,15,21,22,28,
     7   36,45, 2,12,19,23,39, 4,13,24,32,38, 7,17,25,31, 3, 6,11,21,
     8   26,36, 8,16,20,27,44, 1, 3, 6,10,15,21,22,28,36,45, 9,29,33,
     9   35,42,18,30,40, 7,14,17,25,31, 4,13,24,32,38, 9,29,33,35,42,
     A    5,34,37,43, 9,29,33,35,42, 1, 3, 6,10,11,15,21,22,26,28,36,
     B   45, 5,34,37,43, 4,13,24,32,38, 2,12,19,23,39,18,30,40,41, 9,
     C   29,33,35,42, 5,34,37,43, 8,16,20,27,44, 1, 3, 6,10,15,21,22,
     D   28,36,45/
      DATA INTREP/
     1    1, 1, 1, 1, 1, 3, 3, 8, 3, 9, 6, 6,12,14,13, 7, 6,15, 8, 3,
     2    3,11, 9,14,17, 6, 7,12,18,13, 6, 6, 3, 2, 3, 9,11,10,11, 9,
     3   16,10,11, 7, 6, 4, 5, 6, 7, 9,17,19,32,22,40, 3,33,34,27,46,
     4   15,33,28,41,47,35,35,42, 1, 6, 6, 7,29,38,22,31,38,51, 9,19,
     5   32,21,32, 3,35,33,24,34,35,35,35, 3,34,33,26,34,11,32,44,37,
     6   49, 1, 6, 7, 6,32,38,29,21,39,30,38,38,12,12, 4,22,21,19,20,
     7   21,22, 8,27,26,25,27, 8,28,25,26,27, 2,24,23,24,14,18,22,39,
     8   48,45,10,21,37,36,37, 1,13,13, 5,31,30,20,29,30,31, 9,19,40,
     9   21,32,35,35,35, 3,42,34,24,33, 3,41,26,33,34,16,40,44,43,50,
     A   11,44,32,39,10,21,43,36,37, 1, 7, 6, 6,40,38,38,21,45,30,29,
     B   38, 9,32,19,22, 3,47,27,34,33, 3,46,34,27,33,35,35,35,52,11,
     C   32,50,37,44,14,39,22,48,11,32,49,37,44, 1, 6, 6, 7,51,38,22,
     D   31,38,29/
      DATA INTRF1/
     1   19,19,19,19,19, 3, 3, 8, 3, 3,33,33, 8,27,25,35,33,15, 8, 3,
     2    3,34, 3,27,15,33,35, 8,28,25,33,33, 3, 2, 3, 3,34,24,35, 3,
     3   41,26,35,35,33, 2,23,33,35, 3,15, 1,32,22,40, 3, 6,11,14, 0,
     4   15, 6,18,16, 0, 7,11,16,19,33,33,35,29,44,22,48,44,52, 3, 1,
     5   32,21,32, 3,11, 6,10,11, 7,11,11, 3,11, 6,10,11,34,32,38,37,
     6   50,19,33,35,33,32,44,29,21,37,36,44,44, 8, 8, 2,22,21, 1,20,
     7   21,22, 8,14,10,13,14, 8,18,13,10,14, 2,10, 5,10,27,28,22,37,
     8   31,43,24,21,37,30,39,19,25,25,23,48,36,20,29,36,48, 3, 1,40,
     9   21,32,11, 7,11, 3,16,11,10, 6, 3,16,10, 6,11,41,40,38,45,49,
     A   34,38,32,37,26,21,45,30,37,19,35,33,33,40,44,44,21,43,36,29,
     B   44, 3,32, 1,22, 3, 0,14,11, 6, 3, 0,11,14, 6,11,11, 7,51,35,
     C   32,49,37,38,27,37,22,31,35,32,50,39,38,19,33,33,35,52,44,22,
     D   48,44,29/
      DATA INTRF2/
     1   19,19,19,19,19, 9, 9,12, 9, 3,33,33, 8,27,25,35,33,17,12, 9,
     2    9,35, 3,27,15,33,35, 8,28,25,33,33, 9, 4, 9, 3,35,26,34, 3,
     3   42,24,34,35,33, 2,23,33,35, 3,15,19,32,22,40, 9,33,35,27,47,
     4   17,33,28,42,46,35,34,41,19,33,33,35,29,44,22,48,44,52, 3,19,
     5   32,21,32, 9,34,33,26,35,35,34,34, 9,35,33,24,35,35,32,44,39,
     6    0,19,33,35,33,32,44,29,21,37,36,44,44, 8, 8, 2,22,21,19,20,
     7   21,22,12,27,24,25,27,12,28,25,24,27, 4,26,23,26,27,28,22,37,
     8   48,43,26,21,39,36,37,19,25,25,23,48,36,20,29,36,48, 3,19,40,
     9   21,32,34,35,34, 9,41,35,26,33, 9,42,24,33,35,42,40,44,43, 0,
     A   35,44,32,37,24,21,43,36,39,19,35,33,33,40,44,44,21,43,36,29,
     B   44, 3,32,19,22, 9,46,27,35,33, 9,47,35,27,33,34,34,35,52,34,
     C   32, 0,39,44,27,37,22,48,34,32, 0,37,44,19,33,33,35,52,44,22,
     D   48,44,29/
      END