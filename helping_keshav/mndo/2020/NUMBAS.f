      SUBROUTINE NUMBAS
C     *
C     INITIALIZE THE NUMBER OF BASIS FUNCTIONS PER ATOM.
C     THE DEFAULT VALUES WILL BE REDEFINED, IF NECESSARY.
C     SEE SUBROUTINES PARAMD AND PARAME.
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./DPARM / LORBS(LMZ)
      LORBS(1) = 1
      LORBS(2) = 1
      DO 10 I=3,LMZ
      LORBS(I) = 4
   10 CONTINUE
      LORBS(86)= 1
      RETURN
      END
