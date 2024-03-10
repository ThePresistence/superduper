C     ******************************************************************
      SUBROUTINE PSIDRS(IORBB,AINT1,OUT1)
C
C  This subroutine is not intended to be human-readable.
C  It was generated automagically by Mathematica program
C  coded by: Serge Pachkovsky
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AINT1(461), OUT1(461,3)
      COMMON /PSDROT/ D1(3), D2A(6), D2B(6)
      SAVE /PSDROT/
C
C
C   CORE and SS distributions on B
C
      DO 10 ID=1,3
          DO 1000 II = 1,2
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1000     CONTINUE
          DO 1010 II = 16,17
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1010     CONTINUE
   10 CONTINUE
      IF( IORBB.LE.0 ) RETURN
C
C   SP and PP distributions on B
C
      DO 20 ID=1,3
          DO 1020 II = 51,52
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1020     CONTINUE
          DO 1030 II = 66,67
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1030     CONTINUE
          DO 1040 II = 87,88
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1040     CONTINUE
          DO 1050 II = 124,125
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1050     CONTINUE
   20 CONTINUE
      IF( IORBB.LE.3 ) RETURN
C
C   SD, PD and DD distributions on B
C
      DO 30 ID=1,3
          DO 1060 II = 155,156
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1060     CONTINUE
          DO 1070 II = 204,205
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1070     CONTINUE
          DO 1080 II = 249,250
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1080     CONTINUE
          DO 1090 II = 264,265
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1090     CONTINUE
          DO 1100 II = 309,310
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1100     CONTINUE
          DO 1110 II = 334,335
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1110     CONTINUE
          DO 1120 II = 367,368
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1120     CONTINUE
          DO 1130 II = 406,407
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1130     CONTINUE
          DO 1140 II = 447,448
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1140     CONTINUE
   30 CONTINUE
      RETURN
      END
