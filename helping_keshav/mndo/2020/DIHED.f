      SUBROUTINE DIHED (I,J,K,L,ANGLE,SINTH)
C     *
C     CALCULATION OF DIHEDRAL ANGLE BETWEEN ATOMS I,J,K,L (IN DEGREES).
C     THE DIHEDRAL ANGLE IS DEFINED AS THE ANGLE BETWEEN THE VECTORS
C     1.  I-J AND K-L MEASURED CLOCKWISE        IN THE DIRECTION J TO K,
C     2.  K-L AND I-J MEASURED COUNTERCLOCKWISE IN THE DIRECTION K TO J.
C     CONVENTIONS ACCORDING TO KLYNE AND PRELOG.
C     *
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (SMALL=1.0D-06)
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
C     PUT K AT ORIGIN
      XI1    = COORD(1,I)-COORD(1,K)
      XJ1    = COORD(1,J)-COORD(1,K)
      XL1    = COORD(1,L)-COORD(1,K)
      YI1    = COORD(2,I)-COORD(2,K)
      YJ1    = COORD(2,J)-COORD(2,K)
      YL1    = COORD(2,L)-COORD(2,K)
      ZI1    = COORD(3,I)-COORD(3,K)
      ZJ1    = COORD(3,J)-COORD(3,K)
      ZL1    = COORD(3,L)-COORD(3,K)
C     ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
C     COORDINATES AFTER ROTATION.
C     I   (XI2,YI2,ZI1).
C     J   (0  ,YJ2,0  ).
C     K   (0  ,0  ,0  ).
C     L   (XL2,YL2,ZL1).
      DIST   = SQRT(XJ1**2+YJ1**2+ZJ1**2)
      COSA   = ZJ1/DIST
      DDD    = ONE-COSA**2
      IF(DDD.LE.ZERO) THEN
         YXDIST = ZERO
      ELSE
         YXDIST = DIST*SQRT(DDD)
      ENDIF
      IF(YXDIST.LE.SMALL) THEN
         XI2    = XI1
         XL2    = XL1
         YI2    = YI1
         YL2    = YL1
         COSTH  = COSA
         SINTH  = ZERO
      ELSE
         COSPH  = YJ1/YXDIST
         SINPH  = XJ1/YXDIST
         XI2    = XI1*COSPH-YI1*SINPH
C        XJ2    = XJ1*COSPH-YJ1*SINPH
         XL2    = XL1*COSPH-YL1*SINPH
         YI2    = XI1*SINPH+YI1*COSPH
         YJ2    = XJ1*SINPH+YJ1*COSPH
         YL2    = XL1*SINPH+YL1*COSPH
         COSTH  = COSA
         SINTH  = YJ2/DIST
      ENDIF
C     ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
C     COORDINATES AFTER ROTATION.
C     I   (XI2,YI3,UNKNOWN).
C     J   (0  ,0  ,DIST   ).
C     K   (0  ,0  ,0      ).
C     L   (XL2,YL3,UNKNOWN).
      YI3    = YI2*COSTH-ZI1*SINTH
      YL3    = YL2*COSTH-ZL1*SINTH
C     ANGLE BETWEEN PROJECTIONS OF I-J AND K-L INTO XY-PLANE.
C     ANGLE BETWEEN POINTS (XI2,YI3) - (0,0) - (XL2,YL3).
      AA     = XL2*XL2+YL3*YL3
      BB     = XI2*XI2+YI3*YI3
      IF(AA.LT.SMALL .OR. BB.LT.SMALL) THEN
         ANGLE  = ZERO
      ELSE
         ANORM  = ONE/SQRT(AA)
         BNORM  = ONE/SQRT(BB)
         A1     = XL2*ANORM
         A2     = YL3*ANORM
         B1     = XI2*BNORM
         B2     = YI3*BNORM
         SINTH  = A1*B2-A2*B1
         COSTH  = A1*B1+A2*B2
         IF(COSTH.GT. ONE) COSTH= ONE
         IF(COSTH.LT.-ONE) COSTH=-ONE
         ANGLE  = ACOS(COSTH)
         IF(SINTH.LT.ZERO) ANGLE=TWO*PI-ANGLE
         ANGLE  = ANGLE*AFACT
      ENDIF
      RETURN
      END