      SUBROUTINE DVFILL (NPPA,DIRVEC)
C     *
C     CONSTRUCTION OF A HOMOGENEOUS SET OF POINTS ON A UNIT SPHERE.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     NPPA      NUMBER OF POINTS TO BE GENERATED (I).
C               NPPA MUST BE OF THE FORM 10*(3**K)*(4**L)+2
C               ALLOWED VALUES ARE 12,32,42,92,122,162,...
C               DEFAULT VALUE FOR THE FINE GRID : 1082 (K=3,L=1).
C     DIRVEC    CARTESIAN COORDINATES OF THE POINTS GENERATED (O).
C               THESE COORDINATES (X,Y,Z) DEFINE THE DIRECTION VECTORS
C               FROM THE ORIGIN TO THE POINTS ON THE UNIT SPHERE.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./NBFILE/ NBF(20)
      DIMENSION DIRVEC(3,NPPA)
      INTEGER FSET(3,20),KSET(2,30)
      DATA KSET/ 1, 2, 1, 3, 1, 4, 1, 5, 1, 6,
     1          12,11,12,10,12, 9,12, 8,12, 7,
     2           2, 3, 3, 4, 4, 5, 5, 6, 6, 2,
     3           7, 8, 8, 9, 9,10,10,11,11, 7,
     4           2,7,7,3,3,8,8,4,4,9,9,5,5,10,10,6,6,11,11,2/
      DATA FSET/ 1, 2, 3, 1, 3, 4, 1, 4, 5, 1, 5, 6, 1, 6, 2,
     1          12,11,10,12,10, 9,12, 9, 8,12, 8, 7,12, 7,11,
     2           2, 3, 7, 3, 4, 8, 4, 5, 9, 5, 6,10, 6, 2,11,
     3           7, 8, 3, 8, 9, 4, 9,10, 5,10,11, 6,11, 7, 2/
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** FIRST 12 POINTS : CORNERS OF A REGULAR ICOSAHEDRON.
C     EQUIVALENT VIEW : FACE-CENTERS OF A REGULAR DODECAHEDRON.
C     CONVENTION  : POINT  1 LIES ON THE NEGATIVE X-AXIS.
C     CONVENTION  : POINT 12 LIES ON THE POSITIVE X-AXIS.
      DIRVEC(1,1) =-ONE
      DIRVEC(2,1) = ZERO
      DIRVEC(3,1) = ZERO
      ND = 1
      R  = SQRT(0.8D0)
      H  = SQRT(0.2D0)
      DO 11 I=-1,1,2
      DO 10 J= 1,5
      ND = ND+1
C     BETA = ONE + J*1.25663706D0 + (I+1)*0.3141593D0
      BETA = ONE + J*0.4D0*PI     + (I+1)*0.1D0*PI
      DIRVEC(2,ND) = R*COS(BETA)
      DIRVEC(3,ND) = R*SIN(BETA)
      DIRVEC(1,ND) = I*H
   10 CONTINUE
   11 CONTINUE
      DIRVEC(2,12) = ZERO
      DIRVEC(3,12) = ZERO
      DIRVEC(1,12) = ONE
      ND = 12
C *** CHECK ON THE VALIDITY OF CURRENT NPPA VALUE.
C     NPPA MUST BE OF THE FORM 10*(3**K)*(4**L)+2.
C     K AND L ARE DETERMINED DURING THIS CHECK.
      M = (NPPA-2)/10
      DO 20 K=0,10
      IF((M/3)*3.NE.M) GO TO 30
      M = M/3
   20 CONTINUE
   30 DO 40 L=0,10
      IF((M/4)*4.NE.M) GO TO 50
      M = M/4
   40 CONTINUE
   50 IF((10*3**K*4**L+2).NE.NPPA) THEN
         WRITE(NB6,500) NPPA
         STOP 'DVFILL'
      ENDIF
C *** CREATE NEW POINTS ON EACH EXISTING EDGE.
C     THERE ARE 30*(M-1) NEW POINTS, M=(2**L)*(3**KH).
C     NPPA =   12, K=0, L=0, KH=0, M=1,   0 NEW POINTS.
C     NPPA =   42, K=0, L=1, KH=0, M=2,  30 NEW POINTS.
C     NPPA = 1082, K=3, L=1, KH=1, M=6, 150 NEW POINTS.
      KH = K/2
      M  = 2**L*3**KH
      DO 71 I=1,30
      NA = KSET(1,I)
      NB = KSET(2,I)
      DO 70 J=1,M-1
      ND = ND+1
      DO 60 IX=1,3
      DIRVEC(IX,ND) = DIRVEC(IX,NA)*(M-J)+DIRVEC(IX,NB)*J
   60 CONTINUE
   70 CONTINUE
   71 CONTINUE
C     WRITE(NB6,900) NPPA,K,L,KH,M,ND
C *** CREATE POINTS WITHIN EACH EXISTING TRIANGLE.
C     THERE ARE 10*(M-2)*(M-1) NEW POINTS, M=(2**L)*(3**KH).
C     NPPA =   12, K=0, L=0, KH=0, M=1,   0 NEW POINTS.
C     NPPA =   42, K=0, L=1, KH=0, M=2,   0 NEW POINTS.
C     NPPA = 1082, K=3, L=1, KH=1, M=6, 200 NEW POINTS.
      DO 92 I=1,20
      NA = FSET(1,I)
      NB = FSET(2,I)
      NC = FSET(3,I)
      DO 91 J1=1,M-1
      DO 90 J2=1,M-J1-1
      ND = ND+1
      DO 80 IX=1,3
      DIRVEC(IX,ND) = DIRVEC(IX,NA)*(M-J1-J2)
     1               +DIRVEC(IX,NB)*J1+DIRVEC(IX,NC)*J2
   80 CONTINUE
   90 CONTINUE
   91 CONTINUE
   92 CONTINUE
C     WRITE(NB6,900) NPPA,K,L,KH,M,ND
C *** CREATE ADDITIONAL SUBGRIDS (FOR K ODD).
C     THERE ARE 10*M*(M+1) NEW POINTS, M=(2**L)*(3**KH).
C     NPPA = 1082, K=3, L=1, KH=1, M=6, 420 NEW POINTS.
      IF(K.EQ.2*KH) GO TO 140
      T = ONE/THREE
      DO 112 I=1,20
      NA = FSET(1,I)
      NB = FSET(2,I)
      NC = FSET(3,I)
      DO 111 J1=0,M-1
      DO 110 J2=0,M-J1-1
      ND = ND+1
      DO 100 IX=1,3
      DIRVEC(IX,ND) = DIRVEC(IX,NA)*(M-J1-J2-2*T)
     1               +DIRVEC(IX,NB)*(J1+T)+DIRVEC(IX,NC)*(J2+T)
  100 CONTINUE
  110 CONTINUE
  111 CONTINUE
  112 CONTINUE
C     WRITE(NB6,900) NPPA,K,L,KH,M,ND
C     THERE ARE 10*M*(M-1) NEW POINTS, M=(2**L)*(3**KH).
C     NPPA = 1082, K=3, L=2, KH=1, M=6, 300 NEW POINTS.
      T = TWO/THREE
      DO 132 I=1,20
      NA = FSET(1,I)
      NB = FSET(2,I)
      NC = FSET(3,I)
      DO 131 J1=0,M-2
      DO 130 J2=0,M-J1-2
      ND = ND+1
      DO 120 IX=1,3
      DIRVEC(IX,ND) = DIRVEC(IX,NA)*(M-J1-J2-2*T)
     1               +DIRVEC(IX,NB)*(J1+T)+DIRVEC(IX,NC)*(J2+T)
  120 CONTINUE
  130 CONTINUE
  131 CONTINUE
  132 CONTINUE
C     WRITE(NB6,900) NPPA,K,L,KH,M,ND
C *** NORMALIZE ALL VECTORS.
  140 CONTINUE
      DO 170 I=1,NPPA
      DIST = ZERO
      DO 150 IX=1,3
      DIST = DIST+DIRVEC(IX,I)**2
  150 CONTINUE
      DIST = ONE/SQRT(DIST)
      DO 160 IX=1,3
      DIRVEC(IX,I) = DIRVEC(IX,I)*DIST
  160 CONTINUE
  170 CONTINUE
      RETURN
  500 FORMAT(//1X,'VALUE OF NPPA NOT ALLOWED:',I10,
     1       / 1X,'NPPA MUST BE OF THE FORM : 10*(3**K)*(4**L)+2'/)
C 900 FORMAT(  1X,'DVFILL: NPPA,K,L,KH,M,ND',I7,4I4,I8)
      END