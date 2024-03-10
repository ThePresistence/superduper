C     ******************************************************************
      SUBROUTINE MATCOP (A,B,LM2,N,MODE)
C     *
C     MATRIX COPY: B=A OR B=-A.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     A(LM2,*)  INPUT MATRIX (I).
C     B(LM2,*)  COPIED MATRIX (O).
C     LM2       LEADING DIMENSION (I).
C     N         NUMBER OF ORBITALS (I).
C     MODE      SIGN OF COPIED MATRIX (I).
C               = 1  POSITIVE (DEFAULT).
C               =-1  NEGATIVE.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LM2,N),B(LM2,N)
C     USE BLAS ROUTINE DCOPY FOR MODE=1 AND LM2=N.
      IF(MODE.GE.0 .AND. LM2.EQ.N) THEN
         CALL DCOPY (N*N,A,1,B,1)
      ELSE IF(MODE.GE.0) THEN
C     USE FORTRAN CODE OTHERWISE.
         DO 20 J=1,N
         DO 10 I=1,N
         B(I,J) = A(I,J)
   10    CONTINUE
   20    CONTINUE
      ELSE
         DO 40 J=1,N
         DO 30 I=1,N
         B(I,J) = -A(I,J)
   30    CONTINUE
   40    CONTINUE
      ENDIF
      RETURN
      END
