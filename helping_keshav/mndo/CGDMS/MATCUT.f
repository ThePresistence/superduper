C     ******************************************************************
      SUBROUTINE MATCUT (A,CUTOFF,LM2,N,SPARSE)
C     *
C     SET MATRIX ELEMENTS A(I,J) BELOW CUTOFF TO ZERO.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     A(LM2,*)  MATRIX (I,O).
C     CUTOFF    CUTOFF VALUE (I)
C     LM2       LEADING DIMENSION (I).
C     N         NUMBER OF ORBITALS (I).
C     SPARSE    PERCENTAGE OF ELEMENTS BELOW CUTOFF.
C     *
C     ARRAYS A AND B MAY OCCUPY THE SAME CORE SPACE.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      DIMENSION A(LM2,N)
      NCUT = 0
      DO 20 J=1,N
      DO 10 I=1,N
      IF(ABS(A(I,J)).LT.CUTOFF) THEN
         A(I,J) = ZERO
         NCUT   = NCUT+1
      ENDIF
   10 CONTINUE
   20 CONTINUE
      SPARSE = 100.0D0*DBLE(NCUT)/DBLE(N*N)
      RETURN
      END
