      SUBROUTINE TRACEAP (NO,A,IAA,KKA,Tr)
C------------------------------------------------------------------------
C     calculates trace of a matrix stored in compressed sparse row format
C
C     on input:
C                A,IAA,KKA - matrix stored in compressed sparse row format
C                NO        - dimension of equivalent square matrix
C     on output:
C                Tr        - trace of matrix on input
C------------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION Tr
      INTEGER I,J,NO
      DOUBLE PRECISION, DIMENSION (:), POINTER :: A
      INTEGER, DIMENSION (:), POINTER :: IAA,KKA
C
      Tr=0.0D0
      DO 134 I=1,NO
      DO 135 J=KKA(I),KKA(I+1)-1
      IF (IAA(J).EQ.I) THEN
         Tr=Tr+A(J)
         GOTO 134
      ENDIF
 135  CONTINUE
 134  CONTINUE
      RETURN
      END
