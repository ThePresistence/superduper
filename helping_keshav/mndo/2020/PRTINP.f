      SUBROUTINE PRTINP (IGEOM,JPRINT)
C     *
C     PRINT INPUT DATA.
C     *
      USE LIMIT, ONLY: LMV, LMS, LMR, LMG
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 KTITLE,KOMENT
      CHARACTER TEXT(18)*60
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./FLAG1 / KTITLE,KOMENT
     ./INOPT2/ IN2(300)
     ./MOPAC / IMOPAC
     ./NBFILE/ NBF(20)
     ./PARM2 / NSYM,LOCPAR(LMS),IDEPFN(LMS),LOCDEP(LMS)
     ./PARM3 / LOC(LMV),NV
     ./PARM4 / RC(LMR),LREACT,LTOTAL
     ./PARM5 / RC1(LMG),RC2(LMG),LGRID1,LTOT1,LGRID2,LTOT2
     ./PARM6 / LPT,LPT1,LPT2
     ./PARM7 / DEPFAC
      DIMENSION KDLIST(33)
      DATA TEXT/
     1' BOND LENGTH    IS SET EQUAL TO THE REFERENCE BOND LENGTH   ',
     2' BOND ANGLE     IS SET EQUAL TO THE REFERENCE BOND ANGLE    ',
     3' DIHEDRAL ANGLE IS SET EQUAL TO THE REFERENCE DIHEDRAL ANGLE',
     4' DIHEDRAL ANGLE VARIES AS  90 DEGREES - REFERENCE DIHEDRAL  ',
     5' DIHEDRAL ANGLE VARIES AS  90 DEGREES + REFERENCE DIHEDRAL  ',
     6' DIHEDRAL ANGLE VARIES AS 120 DEGREES - REFERENCE DIHEDRAL  ',
     7' DIHEDRAL ANGLE VARIES AS 120 DEGREES + REFERENCE DIHEDRAL  ',
     8' DIHEDRAL ANGLE VARIES AS 180 DEGREES - REFERENCE DIHEDRAL  ',
     9' DIHEDRAL ANGLE VARIES AS 180 DEGREES + REFERENCE DIHEDRAL  ',
     1' DIHEDRAL ANGLE VARIES AS 240 DEGREES - REFERENCE DIHEDRAL  ',
     2' DIHEDRAL ANGLE VARIES AS 240 DEGREES + REFERENCE DIHEDRAL  ',
     3' DIHEDRAL ANGLE VARIES AS 270 DEGREES - REFERENCE DIHEDRAL  ',
     4' DIHEDRAL ANGLE VARIES AS 270 DEGREES - REFERENCE DIHEDRAL  ',
     5' DIHEDRAL ANGLE VARIES AS - REFERENCE DIHEDRAL              ',
     6' BOND LENGTH VARIES AS HALF THE REFERENCE BOND LENGTH       ',
     7' BOND ANGLE VARIES AS HALF THE REFERENCE BOND ANGLE         ',
     8' BOND ANGLE VARIES AS 180 DEGREES - REFERENCE BOND ANGLE    ',
     9' BOND LENGTH IS MULTIPLIED BY A USER-SUPPLIED FACTOR        '/
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** INPUT OPTIONS.
      KGEOM  = IN2(68)
C *** PRINT SYMMETRY DATA.
      IF(JPRINT.GE.0) THEN
         IF(NSYM.GT.0) THEN
            WRITE(NB6,100) KOMENT,KTITLE
            WRITE(NB6,110)
            DO 10 I=1,33
            KDLIST(I) = 0
   10       CONTINUE
            DO 20 I=1,NSYM
            IV2 = IDEPFN(I)
            KDLIST(IV2) = 1
            WRITE(NB6,120) LOCPAR(I),IV2,LOCDEP(I)
   20       CONTINUE
            IF(KDLIST(18).EQ.1) WRITE(NB6,130) DEPFAC
            IF(IMOPAC.EQ.1 .OR. JPRINT.GT.0) THEN
               WRITE(NB6,140)
               DO 30 J=1,18
               IF(KDLIST(J).EQ.1) WRITE(NB6,150) J,TEXT(J)
   30          CONTINUE
            ENDIF
         ENDIF
         WRITE(NB6,100) KOMENT,KTITLE
         WRITE(NB6,160)
      ENDIF
C *** IMPOSE SYMMETRY AND PRINT INITIAL GEOMETRY.
      CALL SYMTRY(JPRINT)
      IF(JPRINT.GE.0) WRITE(NB6,170) NV
C *** PRINT REACTION PATH DATA.
      IF(KGEOM.NE.4) THEN
         LPT    = 0
         LPT1   = 0
         LPT2   = 0
      ENDIF
      IF(JPRINT.GE.0) THEN
         IF(KGEOM.EQ.1 .AND. LTOTAL.GT.0) THEN
            L1  = (LREACT+2)/3
            L2  = LREACT-3*(L1-1)
            WRITE(NB6,180) L2,L1,LTOTAL
            IF(IGEOM.EQ.0 .AND. L2.GT.1) THEN
               WRITE(NB6,190) (RC(I)*AFACT,I=1,LTOTAL)
            ELSE
               WRITE(NB6,190) (RC(I),I=1,LTOTAL)
            ENDIF
         ENDIF
         IF(LTOT1.GT.0 .AND. LTOT2.GT.0) THEN
            L1  = (LGRID1+2)/3
            L2  = LGRID1-3*(L1-1)
            M1  = (LGRID2+2)/3
            M2  = LGRID2-3*(M1-1)
            WRITE(NB6,200) L2,L1,LTOT1,M2,M1,LTOT2
C           BUG FIX BY LASSE SPOERKEL.
C           USE CORRECT CONVERSION FACTORS IF
C           THE TYPES OF BOTH COORDINATES DIFFER.
            IF(IGEOM.EQ.0 .AND. L2.GT.1) THEN
               WRITE(NB6,210) (RC1(I)*AFACT,I=1,LTOT1)
            ELSE
               WRITE(NB6,210) (RC1(I),I=1,LTOT1)
            ENDIF
            IF(IGEOM.EQ.0 .AND. M2.GT.1) THEN
               WRITE(NB6,220) (RC2(I)*AFACT,I=1,LTOT2)
            ELSE
               WRITE(NB6,220) (RC2(I),I=1,LTOT2)
            ENDIF
         ENDIF
         IF(KGEOM.EQ.4 .AND. LTOTAL.GT.0) THEN
            WRITE(NB6,230) LTOTAL,LPT,LPT1,LPT2
            WRITE(NB6,240) (RC(I),I=1,LTOTAL)
         ENDIF
      ENDIF
      RETURN
  100 FORMAT(///1X,A/1X,A)
  110 FORMAT(// 5X,'SYMMETRY CONDITIONS',
     1       // 5X,'REFERENCE     SYMMETRY     DEPENDENT',
     2       /  5X,'  ATOM        RELATION       ATOM   '/)
  120 FORMAT(   5X,I5,9X,I5,8X,I5)
  130 FORMAT(// 5X,'FACTOR FOR SYMMETRY RELATION 18 IS',F20.10)
  140 FORMAT(// 5X,'DESCRIPTION OF THE SYMMETRY RELATIONS USED'/)
  150 FORMAT(   5X,I5,5X,A)
  160 FORMAT(// 5X,'INPUT GEOMETRY'/)
  170 FORMAT(///5X,'**********',//5X,'*  ',
     1             'VARIABLES TO BE OPTIMIZED, OF WHICH THERE ARE',I4)
  180 FORMAT(// 5X,'REACTION COORDINATE: VARIABLE',I2,' OF ATOM',I3,
     1             ', NUMBER OF POINTS:',I4)
  190 FORMAT(/  5X,'VALUES OF THE REACTION COORDINATE ALONG THE PATH.',
     1       / (5X,8F10.5))
  200 FORMAT(// 5X,'FIRST  REACTION COORDINATE: VARIABLE',I2,' OF ATOM',
     1              I3,', NUMBER OF POINTS:',I3,
     2       /  5X,'SECOND REACTION COORDINATE: VARIABLE',I2,' OF ATOM',
     3              I3,', NUMBER OF POINTS:',I3)
  210 FORMAT(/  5X,'VALUES OF THE FIRST REACTION COORDINATE.',
     1       / (5X,8F10.5))
  220 FORMAT(/  5X,'VALUES OF THE SECOND REACTION COORDINATE.',
     1       / (5X,8F10.5))
  230 FORMAT(// 5X,'LINEAR INTERPOLATION BETWEEN TWO INPUT GEOMETRIES:',
     1             ' TOTAL NUMBER OF POINTS:',I4,
     2       /  5X,'PHASE FACTORS XYZ:',3I4)
  240 FORMAT(/  5X,'FACTORS FOR INTERPOLATION ALONG THE PATH.',
     1       / (5X,10F10.4))
      END
