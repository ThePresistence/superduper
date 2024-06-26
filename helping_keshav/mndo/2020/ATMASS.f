      SUBROUTINE ATMASS
C     *
C     DEFINITION OF ATOMIC MASSES.
C     *
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./AMASS / AMS(LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ISTOPE/ CMS(LMZ),BMS(LMZ)
      DO 10 I=1,NUMAT
      AMS(I) = BMS(NAT(I))
   10 CONTINUE
      RETURN
      END
