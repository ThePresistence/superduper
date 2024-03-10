C     ******************************************************************
      SUBROUTINE PONEP (P,PP,LM2,LM6)
C     *
C     EXTRACT ONE-CENTER DENSITY MATRIX ELEMENTS PP(LM6) FROM THE
C     FULL SQUARE DENSITY MATRIX P(LM2,LM2).
C     *
      USE LIMIT, ONLY: LMI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
      DIMENSION P(LM2,LM2),PP(LM6)
      DO 10 IJ=1,LM6
      I      = IP1(IJ)
      J      = IP2(IJ)
      PP(IJ) = P(I,J)
   10 CONTINUE
      RETURN
      END
