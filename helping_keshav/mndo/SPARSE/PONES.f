      SUBROUTINE PONES (P,PS,IPS,JPS,LM6,N)
C     *
C     EXTRACT ONE-CENTER MATRIX ELEMENTS P(LM6) FROM A SPARSE MATRIX
C     WHICH IS STORED IN CSR FORMAT (PS,IPS,JPS).
C     THE ORDER OF COLUMN INDICES IS ARBITRARY.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     P(LM6)    VECTOR OF ONE-CENTER MATRIX ELEMENTS (O).
C     PS(*)     NONZERO VALUES OF SPARSE MATRIX (I).
C     IPS(N+1)  POINTER ARRAY FOR SPARSE MATRIX IN CSR FORMAT (I).
C     JPS(*)    COLUMN INDICES FOR SPARSE MATRIX IN CSR FORMAT (I).
C     *
      USE LIMIT, ONLY: LM1
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
      INTEGER :: I,LM6,N,IJN,IJL,J,IJP,IA,IB,II,NFIRST,NLAST,
     +           NAT,L,NUMAT
      DOUBLE PRECISION :: P(LM6)
      DOUBLE PRECISION, DIMENSION (:), POINTER :: PS
      INTEGER, DIMENSION (:), POINTER :: IPS,JPS
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
C
C *** INITIALIZATION.
      DO 10 I=1,LM6
      P(I) = 0.0D0
   10 CONTINUE
C *** TRANSFER FROM SPARSE STORAGE TO LOCAL ARRAY.
      IJP = 0
      DO 50 II=1,NUMAT
      IA  = NFIRST(II)
      IB  = NLAST(II)
      DO 40 I=IA,IB
      IJN = I-IA+1
      IJL = 0
      DO 20 L=IPS(I),IPS(I+1)-1
      J   = JPS(L)
      IF(J.GE.IA .AND. J.LE.I) THEN
         P(IJP+J-IA+1) = PS(L)
         IJL = IJL+1
         IF(IJL.GE.IJN) GO TO 30
      ENDIF
   20 CONTINUE
   30 CONTINUE
      IJP = IJP+IJN
   40 CONTINUE
   50 CONTINUE
      RETURN
      END