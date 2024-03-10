      DOUBLE PRECISION FUNCTION DDOTS1 (A,IIA,JJA,N,MODE2)
C     *
C     GENERAL CODE FOR SCALAR PRODUCT A*A OF A SPARSE SYMMETRIC MATRIX.
C     THE MATRIX IS STORED IN THE CSR FORMAT.
C     THE ORDER OF COLUMN INDICES IS ARBITRARY.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     A(*)      NONZERO VALUES OF INPUT MATRIX (I).
C     IIA(N+1)  POINTER ARRAY FOR A IN CSR FORMAT (I).
C     JJA(*)    COLUMN INDICES FOR A IN CSR FORMAT (I).
C     N         DIMENSION OF THE MATRIX (I).
C     MODE2     RANGE OF SCALAR PRODUCT (I).
C               = 0 ALL ELEMENTS OF THE MATRIX.
C               = 1 ONLY UPPER TRIANGLE INCLUDING DIAGONAL.
C     DDOTS1    SCALAR PRODUCT (O).
C     *
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
C
      INTEGER NNA,I,JA,N,MODE2,K,IIA(N+1),JJA(IIA(N+1)-1)
      DOUBLE PRECISION DIASUM,A(IIA(N+1)-1)
C
C *** INITIALIZATION.
      DDOTS1 = 0.0D0
      NNA    = IIA(N+1)-1
      IF(NNA.EQ.0) RETURN
C *** COMPUTE SCALAR PRODUCT FOR FULL MATRIX.
      DDOTS1 = DOT_PRODUCT(A(1:NNA),A(1:NNA))
      IF(MODE2.LE.0) RETURN
C *** COMPUTE DIAGONAL SUM.
      DIASUM = 0.0D0
      DO 30 I=1,N
C     LOOP OVER COLUMN INDICES OF MATRIX A.
      DO 20 K=IIA(I),IIA(I+1)-1
      JA     = JJA(K)
      IF(JA.EQ.I) THEN
         DIASUM = DIASUM + A(K)**2
         GO TO 30
      ENDIF
   20 CONTINUE
   30 CONTINUE
C *** COMPUTE SCALAR PRODUCT FOR UPPER TRIANGLE.
      DDOTS1 = 0.5D0*(DDOTS1-DIASUM) + DIASUM
      RETURN
      END
