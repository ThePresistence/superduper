      SUBROUTINE TRACENP (NO,A,IAA,KKA,B,IAB,KKB,Y,Tr)
C--------------------------------------------------------------------------
C     calculates trace of not formed matrix product A*B
C     matrices stored in compressed sparse row format.
C--------------------------------------------------------------------------
C     on input:
C                 A,IAA,KKA - matrix stored in compressed sparse row format
C                 B,IAB,KKB - matrix stored in compressed sparse row format
C                 NO        - dimension of equivalent square matrix
C     on output:
C                 Tr        - trace of not formed matrix product A*B
C---------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NO,I,J
      LOGICAL(1) Z(NO)
      DOUBLE PRECISION Tr
      DOUBLE PRECISION, DIMENSION (:), POINTER :: A,B,Y
      INTEGER, DIMENSION (:), POINTER :: IAA,KKA,IAB,KKB
C
      Tr=0.0D0
C
      DO 235 I=1,NO
      Y(I)=0.
      Z(I)=.FALSE.
 235  CONTINUE
C
      DO 234 I=1,NO
C
      DO 237 J=KKB(I),KKB(I+1)-1
      Y(IAB(J))=B(J)
      Z(IAB(J))=.TRUE.
 237  CONTINUE
C
      DO 236 J=KKA(I),KKA(I+1)-1
      IF (Z(IAA(J))) THEN
         Tr=A(J)*Y(IAA(J))+Tr
      ENDIF
 236  CONTINUE
C
      DO 238 J=KKB(I),KKB(I+1)-1
      Y(IAB(J))=0.
      Z(IAB(J))=.FALSE.
 238  CONTINUE
C
 234  CONTINUE
C
      RETURN
      END
