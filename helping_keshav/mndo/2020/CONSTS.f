      SUBROUTINE CONSTS (NPS,NSET,LMNSET)
C     *
C     CONSTRUCT OR UPDATE THE SOLVENT-ACCESSIBLE SURFACE (SAS).
C     *
C     FIRST CALL  TO CONSTS (NCOSMO.EQ.1): A NEW SURFACE IS CONSTRUCTED.
C     LATER CALLS TO CONSTS (NCOSMO.GT.1): THE OLD SURFACE IS UPDATED.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     NPS       NUMBER OF SEGMENTS ON THE SURFACE (O).
C     NSET      INDEX ARRAY GENERATED DURING SURFACE CONSTRUCTION (O).
C     LMNSET    MAXIMUM DIMENSION OF INDEX ARRAY (I).
C     *
      USE LIMIT, ONLY: LM1, LMNPS, NPPA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL DOUT(NPPA)
      SAVE NN
      SAVE DIRSM,DIRSMH,DIRTM,OCOORD,XA,XSP
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./COSMO2/ NCOSMO,NPSG,NPSH,NPSTOT(2)
     ./COSMO3/ FEPSI,RDS,DISEX2,SRAD(LM1)
     ./COSMO4/ IATSP(LMNPS)
     ./COSMO5/ ATOMAR(LM1),AREA1,NPS1
     ./COSMO8/ COSURF(3,LMNPS),AREASG(LMNPS),QSG(LMNPS)
     ./COSMO9/ DIRVEC(3,NPPA),TM(3,3,LM1),NAR(LMNPS),NSETF(LMNPS)
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
      DIMENSION NSET(LMNSET)
      DIMENSION IPIV(LMNPS),NN(3,LM1)
      DIMENSION DIRSM(3,NPPA),DIRSMH(3,NPPA/3),DIRTM(3,NPPA)
      DIMENSION OCOORD(3,LM1)
      DIMENSION XA(3),XSP(3,LMNPS)
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** INITIALIZATION.
      JPRINT = IN2(42)
      NSPA   = IN2(231)
      NCOSMO = NCOSMO+1
      IF(LMNPS.LT.3*NUMAT) THEN
         WRITE(5,500)
         STOP 'CONSTS'
      ENDIF
C *** DEFINE INITIAL SEGMENTATIONS AND GRIDS ON THE COSMO SURFACE.
      IF(NCOSMO.EQ.1) THEN
C        GENERAL SEGMENTATION FOR ATOMS (DIRSM).
         CALL DVFILL (NPSG,DIRSM)
C        SMALLER SEGMENTATION FOR HYDROGEN (DIRSMH).
         CALL DVFILL (NPSH,DIRSMH)
C        FINE GRID ON THE UNIT SPHERE.
         CALL DVFILL (NPPA,DIRVEC)
      ENDIF
C *** INITIAL SEARCH FOR 3 NEAREST NEIGHBOURS (NN(K,I)) OF ATOM I.
      IF(NCOSMO.EQ.1) THEN
         DO 11 I=1,NUMAT
         NN1    = 0
         NN2    = 0
         NN3    = 0
         DIST1  = 1.0D20
         DIST2  = 1.0D20
         DIST3  = 1.0D20
         DO 10 J=1,NUMAT
         IF(J.EQ.I) GO TO 10
         DIST   = (COORD(1,I)-COORD(1,J))**2
         DIST   = DIST+(COORD(2,I)-COORD(2,J))**2
         DIST   = DIST+(COORD(3,I)-COORD(3,J))**2
         IF((DIST+0.05D0).LT.DIST3) THEN
            DIST3 = DIST
            NN3   = J
         ENDIF
         IF((DIST3+0.05D0).LT.DIST2) THEN
            DIST  = DIST2
            DIST2 = DIST3
            DIST3 = DIST
            NN3   = NN2
            NN2   = J
         ENDIF
         IF((DIST2+0.05D0).LT.DIST1) THEN
            DIST  = DIST1
            DIST1 = DIST2
            DIST2 = DIST
            NN2   = NN1
            NN1   = J
         ENDIF
   10    CONTINUE
         NN(1,I) = NN1
         NN(2,I) = NN2
         NN(3,I) = NN3
   11    CONTINUE
         IF(JPRINT.GT.5) THEN
            WRITE(NB6,510)
            DO 12 I=1,NUMAT
            WRITE(NB6,520) I,NN(1,I),NN(2,I),NN(3,I)
   12       CONTINUE
         ENDIF
      ENDIF
C *** RESTORE THE INFORMATION ON THE PREVIOUS SEGMENTATION IN AN UPDATE.
C     THE PREVIOUS COORDINATES (COSURF) OF THE SEGMENT CENTERS ARE
C     PROJECTED ONTO A UNIT SPHERE CENTERED AT ATOM (I).
C     IATSP(I) : NUMBER OF ATOM (IAT) BELONGING TO A SEGMENT (I).
C     SRAD(IAT): RADIUS OF SAS-SURFACE AT ATOM (IAT).
C     RI       : RADIUS OF THE SURFACE CONTAINING THE SEGMENT CENTERS.
C     COSURF() : DIRECTION VECTORS (COORDINATES) OF THE SEGMENT CENTERS.
C     NPS3     : INDEX FOR THE FIRST SEGMENT OVERALL (REDEFINED LATER).
      IF(NCOSMO.GT.1) THEN
         NPS3 = LMNPS-NPS
         DO 20 I=NPS,1,-1
         IAT = IATSP(I)
         RI  = SRAD(IAT)-RDS
         IATSP(NPS3+I) = IAT
         COSURF(1,NPS3+I) = (COSURF(1,I)-OCOORD(1,IAT))/RI
         COSURF(2,NPS3+I) = (COSURF(2,I)-OCOORD(2,IAT))/RI
         COSURF(3,NPS3+I) = (COSURF(3,I)-OCOORD(3,IAT))/RI
   20    CONTINUE
         NPS3 = NPS3+1
      ENDIF
C *** INITIALIZATION.
      SDIS  = ZERO
      INSET = 1
      NPS   = 0
      AREA  = ZERO
C     *
C *** DEFINE THE SURFACE FOR EACH ATOM.
C     *
      DO 220 I=1,NUMAT
C     DS IS THE MEAN SURFACE AREA OF A SEGMENT (IN UNITS OF PI).
      DS    = SQRT(FOUR/NSPA)
      IF(NAT(I).EQ.1) DS=TWO*DS
C     IT WOULD SEEM MORE LOGICAL TO REPLACE THE TWO PRECEDING LINES BY:
C     IF(NAT(I).EQ.1) THEN
C        DS = SQRT(FOUR/NSPH)
C     ELSE
C        DS = SQRT(FOUR/NSPG)
C     ENDIF
C     THE ORIGINAL CODE IS KEPT FOR THE SAKE OF CONSISTENCY.
      C2DS  = COS(TWO*DS)
      R     = SRAD(I)
      RI    = R-RDS
      XA(1) = COORD(1,I)
      XA(2) = COORD(2,I)
      XA(3) = COORD(3,I)
      OCOORD(1,I) = XA(1)
      OCOORD(2,I) = XA(2)
      OCOORD(3,I) = XA(3)
      NPS0  = NPS+1
C *** UPDATE: THE PREVIOUS COORDINATES (COSURF) OF THE SEGMENT CENTERS
C     OF ATOM (I) ARE ROTATED ON THE UNIT SPHERE USING THE PREVIOUS
C     TRANSFORMATION MATRIX (TM), SUCH THAT THE STANDARD ORIENTATION
C     IS RECOVERED (SEE DVFILL). HENCE, THE NEW COORDINATES (COSURF)
C     SHOULD BE EQUIVALENT TO THE STANDARD COORDINATES (DIRSM,DIRSMH)
C     EXCEPT FOR THE REMOVAL AND RENUMBERING OF SOME COORDINATES.
C     NOTE THAT THE ROTATIONAL TRANSFORMATION IN THE DO 50 LOOP IS
C     THE REVERSE OF THAT IN THE DO 70 LOOP.
C     NPS0  : INDEX FOR FIRST SEGMENT OF ATOM (I).
C     NPS2  : INDEX FOR FIRST SEGMENT OF ATOM (I) IN AN UPDATE.
C     NPS3  : INDEX FOR FIRST SEGMENT OF ATOM (I+1) (REDEFINED HERE).
      IF(NCOSMO.GT.1) THEN
         IF(NPS.GE.NPS3) STOP 'NPS.GT.NPS3'
         NPS2=NPS3
         DO 30 IPS=NPS2,LMNPS
         IF(IATSP(IPS).NE.I) GO TO 40
   30    CONTINUE
   40    NPS3=IPS
C        TRANSFORM COSURF ACCORDING TO TM(INV).
         DO 50 J=NPS2,NPS3-1
         XX1=COSURF(1,J)
         XX2=COSURF(2,J)
         XX3=COSURF(3,J)
         COSURF(1,J)=XX1*TM(1,1,I)+XX2*TM(1,2,I)+XX3*TM(1,3,I)
         COSURF(2,J)=XX1*TM(2,1,I)+XX2*TM(2,2,I)+XX3*TM(2,3,I)
         COSURF(3,J)=XX1*TM(3,1,I)+XX2*TM(3,2,I)+XX3*TM(3,3,I)
   50    CONTINUE
      ENDIF
C *** BUILD NEW TRANSFORMATION MATRIX (TM) FOR THE LOCAL COORDINATE
C     SYSTEM DEFINED BY THE THREE NEAREST NEIGHBOURS (NN1,NN2,NN3).
      NN1=NN(1,I)
      NN2=NN(2,I)
      NN3=NN(3,I)
      IF(NN1.EQ.0) THEN
         TM(1,1,I)=ONE
         TM(1,2,I)=ZERO
         TM(1,3,I)=ZERO
      ELSE
         DIST1=(XA(1)-COORD(1,NN1))**2
         DIST1=DIST1+(XA(2)-COORD(2,NN1))**2
         DIST1=DIST1+(XA(3)-COORD(3,NN1))**2
         DIST=ONE/SQRT(DIST1)
         TM(1,1,I)=(COORD(1,NN1)-XA(1))*DIST
         TM(1,2,I)=(COORD(2,NN1)-XA(2))*DIST
         TM(1,3,I)=(COORD(3,NN1)-XA(3))*DIST
      ENDIF
C     THE CODE UP TO NEXT DO-LOOP HAS BEEN CHANGED SLIGHTLY
C     COMPARED WITH OLDER VERSIONS.
   60 IF(NN2.EQ.0) THEN
         XX1=TM(1,3,I)
         XX2=TM(1,1,I)
         XX3=TM(1,2,I)
         SP=XX1*TM(1,1,I)+XX2*TM(1,2,I)+XX3*TM(1,3,I)
      ELSE
         DIST2=(XA(1)-COORD(1,NN2))**2
         DIST2=DIST2+(XA(2)-COORD(2,NN2))**2
         DIST2=DIST2+(XA(3)-COORD(3,NN2))**2
         DIST=ONE/SQRT(DIST2)
         XX1=(COORD(1,NN2)-XA(1))*DIST
         XX2=(COORD(2,NN2)-XA(2))*DIST
         XX3=(COORD(3,NN2)-XA(3))*DIST
         SP=XX1*TM(1,1,I)+XX2*TM(1,2,I)+XX3*TM(1,3,I)
         IF(SP*SP.GT.0.99D0) THEN
            NN2=NN3
            NN3=0
            DIST2=DIST3
            GO TO 60
         ENDIF
      ENDIF
C     COMPUTE THE REMAINING MATRIX ELEMENTS.
      SININV=ONE/SQRT(ONE-SP*SP)
      TM(2,1,I)=(XX1-SP*TM(1,1,I))*SININV
      TM(2,2,I)=(XX2-SP*TM(1,2,I))*SININV
      TM(2,3,I)=(XX3-SP*TM(1,3,I))*SININV
      TM(3,1,I)=TM(1,2,I)*TM(2,3,I)-TM(2,2,I)*TM(1,3,I)
      TM(3,2,I)=TM(1,3,I)*TM(2,1,I)-TM(2,3,I)*TM(1,1,I)
      TM(3,3,I)=TM(1,1,I)*TM(2,2,I)-TM(2,1,I)*TM(1,2,I)
C *** TRANSFORM THE DIRECTION VECTORS OF THE FINE GRID ACCORDING TO TM.
C     DIRVEC : DIRECTION VECTORS OF GRID POINTS, STANDARD ORIENTATION.
C     DIRTM  : DIRECTION VECTORS OF GRID POINTS, LOCAL COORDINATES.
      DO 70 J=1,NPPA
      XX1=DIRVEC(1,J)
      XX2=DIRVEC(2,J)
      XX3=DIRVEC(3,J)
      DIRTM(1,J)=XX1*TM(1,1,I)+XX2*TM(2,1,I)+XX3*TM(3,1,I)
      DIRTM(2,J)=XX1*TM(1,2,I)+XX2*TM(2,2,I)+XX3*TM(3,2,I)
      DIRTM(3,J)=XX1*TM(1,3,I)+XX2*TM(2,3,I)+XX3*TM(3,3,I)
   70 CONTINUE
C *** FIND THE POINTS OF THE BASIC GRID ON THE SAS.
C     NAREA:NUMBER OF POINTS ASSOCIATED WITH THE CURRENT ATOM.
      NAREA=0
      DO 90 J=1,NPPA
      DOUT(J)=.TRUE.
C     XX* : CARTESIAN COORDINATES OF GRID POINT (J) ON THE SAS.
      XX1 = XA(1) + DIRTM(1,J)*R
      XX2 = XA(2) + DIRTM(2,J)*R
      XX3 = XA(3) + DIRTM(3,J)*R
C     CHECK DISTANCE OF GRID POINT (J) TO ALL ATOMS.
C     DOUT(J) IS TRUE IF THE NEAREST ATOM IS NOT ATOM (I).
C     IN THIS CASE THE GRID POINT (J) IS DISCARDED.
      DO 80 K=1,NUMAT
      IF(K.EQ.I) GO TO 80
      DIST=(XX1-COORD(1,K))**2
      DIST=DIST+(XX2-COORD(2,K))**2
      DIST=DIST+(XX3-COORD(3,K))**2
      DIST=SQRT(DIST)-SRAD(K)
      IF(DIST.LT.ZERO) GO TO 90
   80 CONTINUE
C     ASSIGN GRID POINT (J) TO THE CURRENT ATOM (I), SET DOUT(J).
      NAREA=NAREA+1
      DOUT(J)=.FALSE.
   90 CONTINUE
C     ATOMAR(I):EXPOSED SURFACE FOR ATOM I.
      ATOMAR(I)=FOUR*PI*NAREA*RI*RI/NPPA
C     MOVE TO THE NEXT ATOM IF THE EXPOSED SURFACE IS ZERO.
      IF(NAREA.EQ.0) GO TO 220
      AREA=AREA+NAREA*RI*RI
C *** DEFINE THE POSITION OF THE SEGMENT CENTERS (COSURF)
C     BY A SUITABLE ROTATIONAL TRANSFORMATION ON THE UNIT SPHERE.
C     IN CASE OF AN UPDATE: ROTATE THE OLD DIRECTION VECTORS OF
C     THE SEGMENT CENTERS TO THE NEW LOCAL COORDINATE SYSTEM.
      IF(NCOSMO.GT.1) THEN
         DO 100 J=NPS2,NPS3-1
         NPS=NPS+1
         IATSP(NPS)=I
         XX1=COSURF(1,J)
         XX2=COSURF(2,J)
         XX3=COSURF(3,J)
         COSURF(1,NPS)=XX1*TM(1,1,I)+XX2*TM(2,1,I)+XX3*TM(3,1,I)
         COSURF(2,NPS)=XX1*TM(1,2,I)+XX2*TM(2,2,I)+XX3*TM(3,2,I)
         COSURF(3,NPS)=XX1*TM(1,3,I)+XX2*TM(2,3,I)+XX3*TM(3,3,I)
  100    CONTINUE
C     OTHERWISE: ROTATE THE INITIALLY GENERATED DIRECTION VECTORS OF
C     THE SEGMENT CENTERS TO THE NEW LOCAL COORDINATE SYSTEM.
C     INITIAL VECTORS FOR HYDROGEN   : DIRSMH (NPSH CENTERS).
C     INITIAL VECTORS FOR OTHER ATOMS: DIRSM  (NPSG CENTERS).
      ELSE
         IF(NAT(I).EQ.1) THEN
            IF((NPS+NPSH).GT.LMNPS) THEN
               WRITE(NB6,530)
               STOP 'CONSTS'
            ENDIF
            DO 110 J=1,NPSH
            NPS=NPS+1
            IATSP(NPS)=I
            XX1=DIRSMH(1,J)
            XX2=DIRSMH(2,J)
            XX3=DIRSMH(3,J)
            COSURF(1,NPS)=XX1*TM(1,1,I)+XX2*TM(2,1,I)+XX3*TM(3,1,I)
            COSURF(2,NPS)=XX1*TM(1,2,I)+XX2*TM(2,2,I)+XX3*TM(3,2,I)
            COSURF(3,NPS)=XX1*TM(1,3,I)+XX2*TM(2,3,I)+XX3*TM(3,3,I)
  110       CONTINUE
         ELSE
            IF((NPS+NPSG).GT.LMNPS) THEN
               WRITE(NB6,530)
               STOP 'CONSTS'
            ENDIF
            DO 115 J=1,NPSG
            NPS=NPS+1
            IATSP(NPS)=I
            XX1=DIRSM(1,J)
            XX2=DIRSM(2,J)
            XX3=DIRSM(3,J)
            COSURF(1,NPS)=XX1*TM(1,1,I)+XX2*TM(2,1,I)+XX3*TM(3,1,I)
            COSURF(2,NPS)=XX1*TM(1,2,I)+XX2*TM(2,2,I)+XX3*TM(3,2,I)
            COSURF(3,NPS)=XX1*TM(1,3,I)+XX2*TM(2,3,I)+XX3*TM(3,3,I)
  115       CONTINUE
         ENDIF
      ENDIF
C *** ASSIGN EACH VALID GRID POINT TO A SEGMENT (LOOP 120).
C     NPS0 : INDEX FOR FIRST SEGMENT OF CURRENT ATOM (I).
C     NPS  : INDEX FOR LAST  SEGMENT OF CURRENT ATOM (I).
  120 SDIS0=SDIS
      DO 130 IPS=NPS0,NPS
      NAR(IPS)=0
      XSP(1,IPS)=ZERO
      XSP(2,IPS)=ZERO
      XSP(3,IPS)=ZERO
  130 CONTINUE
C     LOOP OVER ALL POINTS OF THE FINE GRID (J).
      DO 150 J=1,NPPA
      IF(DOUT(J)) GO TO 150
      SPM=-ONE
      X1=DIRTM(1,J)
      X2=DIRTM(2,J)
      X3=DIRTM(3,J)
C     LOOP OVER THE SEGMENTS (IPS) OF THE CURRENT ATOM (I).
C     SP : SCALAR PRODUCT BETWEEN THE DIRECTION VECTORS OF THE CURRENT
C        : GRID POINT (J) AND A SEGMENT CENTER (IPS) WHICH DEFINES THE
C        : COSINE OF THE ANGLE BETWEEN THESE TWO VECTORS.
C     SPM: MAXIMUM COSINE VALUE CORRESPONDING TO THE MINIMUM ANGLE.
C     IPM: ASSOCIATE NUMBER OF THE SEGMENT WHICH IS CLOSEST TO THE
C        : CURRENT GRID POINT (J) ON THE UNIT SPHERE.
C     CHANGE OF IF-STATEMENT ON 14 MAY 1995 BY ACHIM GELESSUS TO OBTAIN
C     HARDWARE-INDEPENDENT RESULTS. ADOPTED BY ANDREAS KLAMT FOR COSMO.
      DO 140 IPS=NPS0,NPS
      SP=X1*COSURF(1,IPS)+X2*COSURF(2,IPS)+X3*COSURF(3,IPS)
      IF((SP-SPM).GE.1.0D-8) THEN
         SPM=SP+1.0D-5
         IPM=IPS
      ENDIF
  140 CONTINUE
C     IF THERE IS NO SEGMENT CENTER WITHIN AN ANGULAR DISTANCE OF C2DS,
C     THE GRID POINT IS TREATED AS A NEW SEGMENT CENTER. IN THIS CASE,
C     THE ASSIGNMENT PROCEDURE IS RESTARTED (GO TO 120).
      IF(SPM.LT.C2DS) THEN
         NPS=NPS+1
         IF(NPS.GT.LMNPS) THEN
            WRITE(NB6,530)
            STOP 'CONSTS'
         ENDIF
         COSURF(1,NPS)=DIRTM(1,J)
         COSURF(2,NPS)=DIRTM(2,J)
         COSURF(3,NPS)=DIRTM(3,J)
         IATSP(NPS)=I
         GO TO 120
      ENDIF
C     ASSIGN THE GRID POINT (J) TO SEGMENT (IPM) AND COUNT THEM (NAR).
C     COMPUTE THE AVERAGE POSITION (XSP) OF THE POINTS (J) ON THE FINE
C     GRID WHICH ARE ASSIGN TO A SEGMENT (IPM).
      NAR(IPM)=NAR(IPM)+1
      XSP(1,IPM)=XSP(1,IPM)+DIRTM(1,J)
      XSP(2,IPM)=XSP(2,IPM)+DIRTM(2,J)
      XSP(3,IPM)=XSP(3,IPM)+DIRTM(3,J)
  150 CONTINUE
      SDIS=ZERO
      IPS=NPS0-1
      IF(NPS.LT.IPS) GO TO 120
C *** DEFINE THE SEGMENT CENTERS (LOOP 160) AT THE CURRENT ATOM (I).
  160 IPS=IPS+1
C     IF SEGMENT (IPS) HAS NO GRID POINTS ASSIGNED, IT IS REMOVED
C     FROM THE LIST OF SEGMENTS, AND THE ASSIGNMENT PROCEDURE IS
C     RESTARTED (GO TO 120).
  170 IF(NAR(IPS).EQ.0) THEN
         NPS=NPS-1
         IF(NPS.LT.IPS) GO TO 120
         DO 180 JPS=IPS,NPS
         NAR(JPS)=NAR(JPS+1)
         XSP(1,JPS)=XSP(1,JPS+1)
         XSP(2,JPS)=XSP(2,JPS+1)
         XSP(3,JPS)=XSP(3,JPS+1)
  180    CONTINUE
         GO TO 170
      ENDIF
C     COMPUTE THE COORDINATES (COSURF) OF THE NEW SEGMENT CENTER (IPS)
C     IN THE LOCAL COORDINATE SYSTEM, BY SCALING THE AVERAGE COORDINATES
C     (XPS) OF THE ASSOCIATE GRID POINTS SUCH THAT THEY CORRESPOND TO
C     A POINT ON THE UNIT SPHERE.
      DIST=XSP(1,IPS)**2+XSP(2,IPS)**2+XSP(3,IPS)**2
      SDIS=SDIS+DIST
      DIST=ONE/SQRT(DIST)
      COSURF(1,IPS)=XSP(1,IPS)*DIST
      COSURF(2,IPS)=XSP(2,IPS)*DIST
      COSURF(3,IPS)=XSP(3,IPS)*DIST
      IF(IPS.LT.NPS) GO TO 160
C     REPEAT THE ASSIGNMENT PROCEDURE IF THE SEGMENTATION HAS NOT BEEN
C     STABLE (GO TO 120).
      IF(ABS(SDIS-SDIS0).GT.1.0D-5) GO TO 120
C *** COMPUTE THE COORDINATES (COSURF) OF ALL SEGMENT CENTERS (IPS)
C     OF THE CURRENT ATOM (I) IN THE MOLECULAR COORDINATE SYSTEM.
C     NSETF(IPS): INDEX OF FIRST GRID POINT OF A GIVEN SEGMENT (IPS).
C                 REFERS TO THE FULL LIST OF ALL VALID GRID POINTS.
C     NAR(IPS)  : NUMBER OF GRID POINTS IN A GIVEN SEGMENT (IPS).
C     XSP()     : DIRECTION VECTOR OF SEGMENT CENTER ON UNIT SPHERE.
C     COSURF()  : ACTUAL CARTESIAN COORDINATES OF SEGMENT CENTER.
      DO 190 IPS=NPS0,NPS
      NSETF(IPS)=INSET
      INSET=INSET+NAR(IPS)
      NAR(IPS)=0
      XSP(1,IPS)=COSURF(1,IPS)
      XSP(2,IPS)=COSURF(2,IPS)
      XSP(3,IPS)=COSURF(3,IPS)
      COSURF(1,IPS)=XA(1)+COSURF(1,IPS)*RI
      COSURF(2,IPS)=XA(2)+COSURF(2,IPS)*RI
      COSURF(3,IPS)=XA(3)+COSURF(3,IPS)*RI
  190 CONTINUE
C *** LOOP OVER ALL GRID POINTS (J) AT THE CURRENT ATOM (I).
C     IDENTIFY THE SEGMENT (IPM) ASSOCIATED WITH GRID POINT (J).
C     IDENTIFY THE RUNNING INDEX (NARA+1) OF (J) IN (IPM).
C     ASSIGN TO EACH VALID GRID POINT IN THE FULL LIST (NSET) OF
C     ALL GRID POINTS ON THE SAS THE INDEX (J) IN THE LIST OF
C     GRID POINTS AT THE CURRENT ATOM (I).
C     TECHNICAL NOTE: THE DIMENSION OF THE SCRATCH ARRAY NSET(NPPA*LM1)
C     IS RATHER LARGE, BUT THERE ARE LESS THAN NPPA*NUMAT ENTRIES.
      DO 210 J=1,NPPA
      IF(DOUT(J)) GO TO 210
      SPM=-ONE
      X1=DIRTM(1,J)
      X2=DIRTM(2,J)
      X3=DIRTM(3,J)
C     ANALOGOUS COMMENTS AS GIVEN FOR THE DO 140 LOOP (SEE ABOVE).
C     THE FOLLOWING 8 LINES INCLUDING THE IF-STATEMENT SHOULD BE
C     REDUNDANT SINCE PATHOLOGICAL CASES SHOULD HAVE BEEN REMOVED
C     AFTER THE DO 140 LOOP.
      DO 200 IPS=NPS0,NPS
      SP=X1*XSP(1,IPS)+X2*XSP(2,IPS)+X3*XSP(3,IPS)
      IF((SP-SPM).GE.1.D-8) THEN
         SPM=SP+1.0D-5
         IPM=IPS
      ENDIF
  200 CONTINUE
      IF(SPM.LT.C2DS) GO TO 210
      NARA=NAR(IPM)
      NSET(NSETF(IPM)+NARA)=J
      NAR(IPM)=NARA+1
  210 CONTINUE
  220 CONTINUE
C *** SAVE RESULTS IN COMMON BLOCK.
C     AREA1  : MOLECULAR SURFACE AREA (IN ANGSTROM**2).
C     NPS1   : TOTAL NUMBER OF SEGMENTS.
      AREA1  = AREA*FOUR*PI/NPPA
      NPS1   = NPS
C *** DEBUG PRINT.
      IF(NCOSMO.EQ.1 .AND. JPRINT.GE.0) THEN
         WRITE(NB6,540) NPS
         IF(JPRINT.GE.5) THEN
            WRITE(NB6,550)
            DO 290 I=1,NPS
            WRITE(NB6,560) I,(COSURF(J,I),J=1,3)
  290       CONTINUE
         ENDIF
      ENDIF
      RETURN
  500 FORMAT(// 1X,'CONSTS: LMNPS MUST BE INCREASED FOR THIS SYSTEM'/)
  510 FORMAT(/  1X,'NEAREST NEIGHBOURS FOR LOCAL TRANSFORMATIONS',
     1          1X,'DURING COSMO SURFACE CONSTRUCTION.',
     2       // 1X,'     I   NN1   NN2   NN3'/)
  520 FORMAT(   1X,4I6)
  530 FORMAT(// 1X,'CONSTS: NPS IS GREATER THAN LMNPS.',
     1       /  1X,'USE SMALLER INPUT VALUE FOR NSPA.'/)
  540 FORMAT(// 1X,'NUMBER OF SEGMENTS ON THE MOLECULAR SURFACE',I10)
  550 FORMAT(///1X,'COORDINATES OF SEGMENT CENTERS (ANGSTROM).'/)
  560 FORMAT(   1X,I5,3F20.10)
      END