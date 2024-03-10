C     ******************************************************************
      SUBROUTINE PSIDRP(IORBB,AINT1,OUT1)
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
          DO 1000 II = 1,6
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1000     CONTINUE
          DO 1010 II = 16,21
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1010     CONTINUE
   10 CONTINUE
      IF( IORBB.LE.0 ) RETURN
C
C   SP and PP distributions on B
C
      DO 20 ID=1,3
          DO 1020 II = 31,32
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1020     CONTINUE
          DO 1030 II = 41,42
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1030     CONTINUE
          DO 1040 II = 51,56
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1040     CONTINUE
          DO 1050 II = 66,71
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1050     CONTINUE
              OUT1(83,ID) = AINT1(83)*D1(ID)
          DO 1060 II = 87,92
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1060     CONTINUE
          DO 1070 II = 104,105
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1070     CONTINUE
          DO 1080 II = 114,115
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1080     CONTINUE
          DO 1090 II = 124,129
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1090     CONTINUE
   20 CONTINUE
      IF( IORBB.LE.3 ) RETURN
C
C   SD, PD and DD distributions on B
C
      DO 30 ID=1,3
          DO 1100 II = 139,140
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1100     CONTINUE
          DO 1110 II = 145,146
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1110     CONTINUE
          DO 1120 II = 155,160
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1120     CONTINUE
          DO 1130 II = 170,171
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1130     CONTINUE
              OUT1(180,ID) = AINT1(180)*D1(ID)
          DO 1140 II = 184,185
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1140     CONTINUE
          DO 1150 II = 194,195
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1150     CONTINUE
          DO 1160 II = 204,209
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1160     CONTINUE
          DO 1170 II = 219,220
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1170     CONTINUE
          DO 1180 II = 229,230
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1180     CONTINUE
          DO 1190 II = 239,240
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1190     CONTINUE
          DO 1200 II = 249,254
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1200     CONTINUE
          DO 1210 II = 264,269
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1210     CONTINUE
          DO 1220 II = 279,280
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1220     CONTINUE
          DO 1230 II = 289,290
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1230     CONTINUE
          DO 1240 II = 299,300
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1240     CONTINUE
          DO 1250 II = 309,314
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1250     CONTINUE
          DO 1260 II = 324,325
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1260     CONTINUE
          DO 1270 II = 334,339
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1270     CONTINUE
          DO 1280 II = 351,352
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1280     CONTINUE
          DO 1290 II = 357,358
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1290     CONTINUE
          DO 1300 II = 367,372
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1300     CONTINUE
          DO 1310 II = 382,383
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1310     CONTINUE
              OUT1(392,ID) = AINT1(392)*D1(ID)
          DO 1320 II = 396,397
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1320     CONTINUE
          DO 1330 II = 406,411
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1330     CONTINUE
          DO 1340 II = 423,424
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1340     CONTINUE
              OUT1(433,ID) = AINT1(433)*D1(ID)
          DO 1350 II = 437,438
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1350     CONTINUE
          DO 1360 II = 447,452
              OUT1(II,ID) = AINT1(II)*D1(ID)
 1360     CONTINUE
   30 CONTINUE
      RETURN
      END
