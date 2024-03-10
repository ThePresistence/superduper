C     ******************************************************************
      SUBROUTINE MATDEB (A,LM2,NR,NC)
C     *
C     DEBUG PRINT FOR MATRIX A.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     A(LM2,*)  INPUT MATRIX (I).
C     LM2       LEADING DIMENSION OF THE MATRIX (I).
C     NR        NUMBER OF ROWS PRINTED (I).
C     NC        NUMBER OF COLUMNS PRINTED (I).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      DIMENSION A(LM2,NC)
      NB6 = 6
C     PRINT ONE NONZERO MATRIX ELEMENT PER LINE.
      DO 20 J=1,NC
      DO 10 I=1,NR
      IF(A(I,J).NE.ZERO) WRITE(NB6,500) I,J,A(I,J)
   10 CONTINUE
   20 CONTINUE
      IF(NR.NE.NC) RETURN
C     COMPUTE AND PRINT MAXIMUM ASYMMETRY.
      ASSYM  = ZERO
      DO 40 J=2,NC
      DO 30 I=1,J
      ASSYM  = MAX(ASSYM,ABS((A(I,J)-A(J,I))))
   30 CONTINUE
   40 CONTINUE
      WRITE(NB6,510) ASSYM
      RETURN
  500 FORMAT( 1X,'MATDEB: ',2I5,G20.10)
  510 FORMAT(/1X,'MATDEB: MAXIMUM ASYMMETRY =',G20.10)
      END
