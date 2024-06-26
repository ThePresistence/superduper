      SUBROUTINE MMOKF (NNHCO)
C     *
C     MOLECULAR MECHANICS ENERGY CORRECTION FOR PEPTIDE BONDS.
C     INCLUDED FOR COMPATIBILITY WITH MOPAC(6.0).
C     CODE ADAPTED FROM MOPAC(6.0) WRITTEN BY J.J.P.STEWART.
C     *
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./MMCOM1/ EMM
     ./MMCOM2/ NHCO(4,1+2*LM1/3),HTYPE
      EMM    = 0.0D0
      IF(NNHCO.LE.0) RETURN
      DO 10 I=1,NNHCO
      CALL DIHED(NHCO(1,I),NHCO(2,I),NHCO(3,I),NHCO(4,I),ANGLE,SIN1)
      EMM    = EMM+HTYPE*SIN1**2
   10 CONTINUE
      RETURN
      END
