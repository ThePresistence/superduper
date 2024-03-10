C     ******************************************************************
      SUBROUTINE PSIDRE(IORBB,AINT1,AINT2,OUT2)
C
C  This subroutine is not intended to be human-readable.
C  It was generated automagically by Mathematica program
C  coded by: Serge Pachkovsky
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AINT1(461), AINT2(461), OUT2(461,6)
      COMMON /PSDROT/ D1(3), D2A(6), D2B(6)
      SAVE /PSDROT/
C
C
C   CORE and SS distributions on B
C
      DO 10 ID=1,6
          DO 1000 II = 1,30
              OUT2(II,ID) = AINT2(II)*D2A(ID) + AINT1(II)*D2B(ID)
 1000     CONTINUE
   10 CONTINUE
      IF( IORBB.LE.0 ) RETURN
C
C   SP and PP distributions on B
C
      DO 20 ID=1,6
          DO 1010 II = 31,138
              OUT2(II,ID) = AINT2(II)*D2A(ID) + AINT1(II)*D2B(ID)
 1010     CONTINUE
   20 CONTINUE
      IF( IORBB.LE.3 ) RETURN
C
C   SD, PD and DD distributions on B
C
      DO 30 ID=1,6
          DO 1020 II = 139,461
              OUT2(II,ID) = AINT2(II)*D2A(ID) + AINT1(II)*D2B(ID)
 1020     CONTINUE
   30 CONTINUE
      RETURN
      END
