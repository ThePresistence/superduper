      SUBROUTINE PDDG (NI,NJ,R,ENUCLR)
C     *
C     MODIFICATION OF CORE REPULSION FUNCTION FOR MNDO AND PM3.
C     *
C     CONTRIBUTIONS FROM PAIRWISE DISTANCE DIRECTED GAUSSIANS.
C     M.P. REPASKY, J. CHANDRASEKHAR, W.L. JORGENSEN, JCC 2002
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./PDDG1 / DA(LMZ,2),PA(LMZ,2)
      DATA CUTOFF/25.0D0/
C *** INITIALIZATION.
      TORENI = TORE(NI)
      TORENJ = TORE(NJ)
      ADD    = ZERO
C *** GENERAL SECTION.
      DO 20 I=1,2
      DO 10 J=1,2
      XX     = 10.0D0*(R-DA(NI,I)-DA(NJ,J))**2
      IF(XX.LT.CUTOFF) THEN
         ADD = ADD + (TORENI*PA(NI,I)+TORENJ*PA(NJ,J))*EXP(-XX)
      ENDIF
   10 CONTINUE
   20 CONTINUE
      ENUCLR = ENUCLR + ADD/(TORENI+TORENJ)
C     WRITE(6,500) NI,NJ,R,ADD/(TORENI+TORENJ)
      RETURN
C 500 FORMAT(1X,'PDDG: NI, NJ, R(A), PDDG-TERM (EV)',2I5,2F20.10)
      END
