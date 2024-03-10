      SUBROUTINE GMETRY (JCXYZ)
C     *
C     CALCULATE CARTESIAN COORDINATES FROM INTERNAL COORDINATES.
C     *
C     NOTATION. I=INPUT,O=OUTPUT.
C     JCXYZ     OPTION FOR TREATMENT OF DUMMY ATOMS (I).
C               =-2 ONLY REMOVE DUMMY ATOMS FROM ARRAY (COORD).
C               = 0 CALCULATE COORDINATES AND REMOVE DUMMY ATOMS.
C               = 1 CALCULATE COORDINATES AND KEEP DUMMY ATOMS.
C     *
C     CONVENTION FOR THE CALLING ROUTINE: UPON EXIT, THERE MUST NOT
C     BE ANY DUMMY ATOMS IN THE COORDINATE ARRAY (COORD). HENCE, ANY
C     CALL WITH JCXYZ=1 MUST BE FOLLOWED BY A CALL WITH JCXYZ=-2.
C     *
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (SMALL =1.0D-07)
      PARAMETER (SMALLT=1.0D-10)
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./INOPT1/ IN1(300)
     ./NBFILE/ NBF(20)
     ./PARM1 / A(3,LM1),NA(LM1),NB(LM1),NC(LM1),NN(LM1),NATOMS
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** INPUT OPTIONS.
C     GMETRY MAY BE CALLED FROM READMO/GETGEO, HENCE USE ARRAY IN1.
      IGEOM  = IN1(4)
C *** CHECK FOR REMOVAL OF DUMMY ATOMS.
      IF(JCXYZ.EQ.-2) GO TO 50
C *** CHECK FOR CARTESIAN INPUT GEOMETRY.
      IF(IGEOM.GT.0) THEN
         DO 10 I=1,NATOMS
         COORD(1,I) = A(1,I)
         COORD(2,I) = A(2,I)
         COORD(3,I) = A(3,I)
   10    CONTINUE
         GO TO 50
      ENDIF
C *** DEFINE FIRST TWO ATOMS.
      COORD(1,1) = ZERO
      COORD(2,1) = ZERO
      COORD(3,1) = ZERO
      COORD(1,2) = A(1,2)
      COORD(2,2) = ZERO
      COORD(3,2) = ZERO
      IF(NATOMS.EQ.2) GO TO 50
C *** CHECK FOR LINEARITY.
      DO 20 I=3,NATOMS
      COSA   = COS(A(2,I))
      TEST   = ONE-ABS(COSA)
      MC     = NC(I)
      MB     = NB(I)
      XB     = COORD(1,MB)-COORD(1,MC)
      IF(TEST.GT.SMALL) GO TO 30
      COORD(1,I) = COORD(1,MC)+A(1,I)*COSA*(XB/ABS(XB))
      COORD(2,I) = ZERO
      COORD(3,I) = ZERO
   20 CONTINUE
      GO TO 50
C *** PUT FIRST NONLINEAR ATOM IN XY-PLANE WITH POSITIVE Y.
C     THIS IS USUALLY THE THIRD ATOM OF THE MOLECULE.
   30 COORD(1,I) = COORD(1,MC)+A(1,I)*COSA*(XB/ABS(XB))
      COORD(2,I) = A(1,I)*SIN(A(2,I))
      COORD(3,I) = ZERO
      IF(I.EQ.NATOMS) GO TO 50
C *** DEFINE THE REMAINING ATOMS.
      II     = I+1
      DO 40 I=II,NATOMS
      COSA   = COS(A(2,I))
      MB     = NB(I)
      MC     = NC(I)
      XB     = COORD(1,MB)-COORD(1,MC)
      YB     = COORD(2,MB)-COORD(2,MC)
      ZB     = COORD(3,MB)-COORD(3,MC)
      RBC    = ONE/SQRT(XB*XB+YB*YB+ZB*ZB)
      TEST   = ONE-ABS(COSA)
C     ATOMS MC, MB, AND (I) ARE COLLINEAR.
      IF(TEST.LT.SMALLT) THEN
         RBC = A(1,I)*RBC*COSA
         COORD(1,I) = COORD(1,MC)+XB*RBC
         COORD(2,I) = COORD(2,MC)+YB*RBC
         COORD(3,I) = COORD(3,MC)+ZB*RBC
         GO TO 40
      ENDIF
C     THE ATOMS ARE NOT COLLINEAR.
      MA     = NA(I)
      XA     = COORD(1,MA)-COORD(1,MC)
      YA     = COORD(2,MA)-COORD(2,MC)
      ZA     = COORD(3,MA)-COORD(3,MC)
C     IN EVALUATING THE COORDINATES OF ATOM I WITH RESPECT TO THE
C     REFERENCE ATOMS MC-MB-MA, WE USE THE FOLLOWING STRATEGY.
C     WE DEFINE A LOCAL COORDINATE SYSTEM WHICH IS ORIENTED SUCH THAT
C     MC=(0,0,0), MB(XLB,0,0), MA(XLA,YLA,0)
C     AND CALCULATE THE COORDINATES OF ATOM I IN THIS LOCAL SYSTEM.
C     THESE COORDINATES ARE THEN TRANSFORMED TO THE ORIGINAL SYSTEM
C     BY THREE ROTATIONS AROUND THE AXES X,Y,Z. THE ROTATION MATRIX
C     ELEMENTS ARE COMPUTED FIRST BY STANDARD METHODS.
      XYB    = SQRT(XB*XB+YB*YB)
C     ROTATION AROUND THE Z-AXIS TO MAKE YB VANISH AND XB POSITIVE.
C     SWITCH ATOMS MB AND MA IF MC-MB IS PARALLEL TO THE Z-AXIS.
      K      = -1
      IF(XYB.LE.SMALL) THEN
         XPA = ZA
         ZA  =-XA
         XA  = XPA
         XPB = ZB
         ZB  =-XB
         XB  = XPB
         XYB = SQRT(XB*XB+YB*YB)
         K   = 1
      ENDIF
      IF(XYB.LT.SMALL) THEN
         WRITE(NB6,500) I,MC,MB,MA
         STOP 'GMETRY'
      ENDIF
      RXYB   = ONE/XYB
      COSTH  = XB*RXYB
      SINTH  = YB*RXYB
      XPA    = XA*COSTH+YA*SINTH
      YPA    = YA*COSTH-XA*SINTH
C     LOCAL COORDINATES - MC(0,0,0), MB(UNDEF,0,ZB), MA(XPA,YPA,ZA).
C     ROTATION ABOUT THE Y-AXIS TO MAKE ZB VANISH.
      SINPH  = ZB*RBC
      COSPH  = SQRT(ABS(ONE-SINPH*SINPH))
C     XQA    = XPA*COSPH+ZA*SINPH
      ZQA    = ZA*COSPH-XPA*SINPH
C     LOCAL COORDINATES - MC(0,0,0), MB(UNDEF,0,0), MA(XQA,YPA,ZQA).
C     ROTATION ABOUT THE X-AXIS TO MAKE ZA VANISH AND YA POSITIVE.
      YZA    = SQRT(YPA**2+ZQA**2)
      IF(YZA.LT.SMALL) THEN
         WRITE(NB6,500) I,MC,MB,MA
         STOP 'GMETRY'
      ENDIF
      RYZA   = ONE/YZA
      COSKH  = YPA*RYZA
      SINKH  = ZQA*RYZA
C     LOCAL COORDINATES - MC(0,0,0), MB(UNDEF,0,0), MA(XQA,UNDEF,0).
C     ALL NONZERO LOCAL COORDINATES ARE POSITIVE.
C     THE COORDINATES OF ATOM I ARE EVALUATED IN THE NEW FRAME.
C     THE SIGN OF ZD CONFORMS TO THE KLYNE-PRELOG CONVENTION FOR
C     DIHEDRAL ANGLES, SEE QCPE BULLETIN 4, 97 (1984).
      SINA   = SIN(A(2,I))
      SIND   = SIN(A(3,I))
      COSD   = COS(A(3,I))
      XD     = A(1,I)*COSA
      YD     = A(1,I)*SINA*COSD
      ZD     =-A(1,I)*SINA*SIND
C     TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
      YPD    =  YD*COSKH- ZD*SINKH
      ZPD    =  ZD*COSKH+ YD*SINKH
      XPD    =  XD*COSPH-ZPD*SINPH
      ZQD    = ZPD*COSPH+ XD*SINPH
      XQD    = XPD*COSTH-YPD*SINTH
      YQD    = YPD*COSTH+XPD*SINTH
      IF(K.EQ.1) THEN
         XRD =-ZQD
         ZQD = XQD
         XQD = XRD
      ENDIF
      COORD(1,I) = XQD+COORD(1,MC)
      COORD(2,I) = YQD+COORD(2,MC)
      COORD(3,I) = ZQD+COORD(3,MC)
   40 CONTINUE
C *** NOW REMOVE THE DUMMY ATOM COORDINATES, IF ANY, FROM COORD.
   50 IF(NATOMS.EQ.NUMAT) RETURN
      IF(JCXYZ.EQ.1) RETURN
      J      = 0
      DO 70 I=1,NATOMS
      IF(NN(I).EQ.99) GO TO 70
      J      = J+1
      COORD(1,J) = COORD(1,I)
      COORD(2,J) = COORD(2,I)
      COORD(3,J) = COORD(3,I)
   70 CONTINUE
      RETURN
  500 FORMAT(///5X,'ERROR IN CALCULATING THE DIHEDRAL ANGLE OF ATOM',I4,
     1       /  5X,'THE REFERENCE ATOMS',3I4,'  ARE COLLINEAR'///)
      END
