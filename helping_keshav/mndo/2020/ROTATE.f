      SUBROUTINE ROTATE(IW,JW,IP,JP,KR,RI,YY,W,WW,LM6,LM9,IMODE)
C     *
C     TWO-ELECTRON REPULSION INTEGRALS.
C     TRANSFORMATION FROM LOCAL TO MOLECULAR COORDINATES.
C     *
C     STORAGE OF THE TRANSFORMED INTEGRALS.
C     OPTION IMODE.LE.0: SQUARE ARRAY WW(LM6,LM6).
C     OPTION IMODE.GT.0: LINEAR ARRAY W(LM9).
C     OPTION IMODE.LE.0: CALLS FROM HCORE AND HCOREP.
C     OPTION IMODE.GT.0: CALLS FROM HCORE AND DHCORE.
C     A GIVEN CALL EITHER REFERS TO W(LM9) OR WW(LM6,LM6).
C     *
C     INPUT DATA.
C     IW,JW    NUMBER OF ONE-CENTER PAIRS AT ATOMS I AND J.
C     IP,JP    ADDRESS OF (SS,SS) IN WW(LM6,LM6).
C     KR+1     ADDRESS OF (SS,SS) IN W(LM9).
C     RI       LOCAL TWO-ELECTRON INTEGRALS.
C     YY       PRECOMBINED ROTATION MATRIX ELEMENTS.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(LM9),WW(LM6,LM6)
      DIMENSION RI(22)
      DIMENSION YY(15,45)
      DIMENSION SUM(6)
      DIMENSION T1(6),T2(6),T3(6),T4(6)
      DIMENSION T5(6),T6(6),T7(6),T8(6),T9(6)
      DIMENSION SSPB(9),PASS(9),PSPS(6)
      DIMENSION PSPP(6,3),PPPS(3,6),PPPP(6,6)
      DIMENSION IPP(3),JPP(3),ISS(6),JSS(6)
      DATA IPP/ 1, 3, 6/
      DATA JPP/10,30,60/
      DATA ISS/ 2, 4, 5, 7, 8, 9/
      DATA JSS/20,40,50,70,80,90/
      IF(IW.EQ.1 .AND. JW.EQ.1) GO TO 110
C *** TRANSFORM THE INTEGRALS.
      DO 10 I=1,6
      SUM(I) = YY(I,6)+YY(I,10)
   10 CONTINUE
C     INTEGRAL TYPES (SS,PS) AND (SS,PP).
      IF(JW.GT.1) THEN
         SSPB(1) = RI(5)*YY(1,2)
         SSPB(3) = RI(5)*YY(2,2)
         SSPB(6) = RI(5)*YY(3,2)
         DO 20 I=1,6
         SSPB(ISS(I)) = RI(11)*YY(I,3)+RI(12)*SUM(I)
   20    CONTINUE
      ENDIF
C     INTEGRAL TYPES (PS,SS) AND (PP,SS).
      IF(IW.GT.1) THEN
         PASS(1) = RI(2)*YY(1,2)
         PASS(3) = RI(2)*YY(2,2)
         PASS(6) = RI(2)*YY(3,2)
         DO 30 I=1,6
         PASS(ISS(I)) = RI(3)*YY(I,3)+RI(4)*SUM(I)
   30    CONTINUE
      ENDIF
      IF(IW.GT.1 .AND. JW.GT.1) THEN
C        INTEGRAL TYPE (PS,PS) AND AUXILIARY TERMS FOR (PS,PP).
         DO 40 I=1,6
         PSPS(I) = RI( 6)*YY(I,3)+RI( 7)*SUM(I)
         T1(I)   = RI(13)*YY(I,3)+RI(14)*SUM(I)
         T2(I)   = RI(15)*YY(I,8)
         T3(I)   = RI(15)*YY(I,5)
   40    CONTINUE
C        INTEGRAL TYPE (PS,PP).
         DO 60 I=1,3
         DO 50 J=1,6
         PSPP(J,I) = YY(I,2)*T1(J)+YY(I,7)*T2(J)+YY(I,4)*T3(J)
   50    CONTINUE
   60    CONTINUE
C        AUXILIARY TERMS FOR (PP,PS), AND (PP,PP).
         DO 70 I=1,6
         T1(I)   = RI( 8)*YY(I,3)+RI( 9)*SUM(I)
         T2(I)   = RI(10)*YY(I,8)
         T3(I)   = RI(10)*YY(I,5)
         T4(I)   = RI(16)*YY(I,3)+RI(17)*SUM(I)
         T5(I)   = RI(18)*YY(I,3)+RI(19)*YY(I,10)+RI(21)*YY(I,6)
         T6(I)   = RI(18)*YY(I,3)+RI(19)*YY(I,6) +RI(21)*YY(I,10)
         T7(I)   = RI(20)*YY(I,8)
         T8(I)   = RI(20)*YY(I,5)
         T9(I)   = RI(22)*YY(I,9)
   70    CONTINUE
C        INTEGRAL TYPES (PP,PS) AND (PP,PP).
         DO 100 I=1,6
         DO 80 J=1,3
         PPPS(J,I) = YY(J,2)*T1(I)+YY(J,7)*T2(I)+YY(J,4)*T3(I)
   80    CONTINUE
         DO 90 J=1,6
         PPPP(J,I) = YY(J,3)*T4(I)+YY(J,10)*T5(I)+YY(J,6)*T6(I)
     1              +YY(J,8)*T7(I)+YY(J,5) *T8(I)+YY(J,9)*T9(I)
   90    CONTINUE
  100    CONTINUE
      ENDIF
C *** STORE INTEGRALS.
  110 CONTINUE
      IF(IMODE.GT.0) THEN
         K    = KR+1
         W(K) = RI(1)
         IF(JW.GT.1) THEN
C           INTEGRAL TYPES (SS,PS) AND (SS,PP).
            DO 120 I=1,9
            W(K+I) = SSPB(I)
  120       CONTINUE
         ENDIF
         IF(IW.GT.1 .AND. JW.EQ.1) THEN
C           INTEGRAL TYPES (PS,SS) AND (PP,SS).
            DO 130 I=1,9
            W(K+I) = PASS(I)
  130       CONTINUE
         ENDIF
         IF(IW.GT.1 .AND. JW.GT.1) THEN
C           INTEGRAL TYPES (PS,SS) AND (PP,SS).
            DO 140 I=1,9
            W(K+I*10) = PASS(I)
  140       CONTINUE
C           INTEGRAL TYPE (PS,PS).
            W(K+11) = PSPS(1)
            W(K+13) = PSPS(2)
            W(K+16) = PSPS(4)
            W(K+31) = PSPS(2)
            W(K+33) = PSPS(3)
            W(K+36) = PSPS(5)
            W(K+61) = PSPS(4)
            W(K+63) = PSPS(5)
            W(K+66) = PSPS(6)
C           INTEGRAL TYPE (PS,PP).
            DO 160 I=1,3
            IJ  = K+JPP(I)
            DO 150 J=1,6
            W(IJ+ISS(J)) = PSPP(J,I)
  150       CONTINUE
  160       CONTINUE
C           INTEGRAL TYPES (PP,PS) AND (PP,PP).
            DO 190 I=1,6
            IJ  = K+JSS(I)
            DO 170 J=1,3
            W(IJ+IPP(J)) = PPPS(J,I)
  170       CONTINUE
            DO 180 J=1,6
            W(IJ+ISS(J)) = PPPP(J,I)
  180       CONTINUE
  190       CONTINUE
         ENDIF
      ELSE
         WW(IP,JP) = RI(1)
         IF(JW.GT.1) THEN
C           INTEGRAL TYPE (SS,PS) AND (SS,PP).
            DO 210 I=1,9
            WW(IP,JP+I) = SSPB(I)
  210       CONTINUE
         ENDIF
         IF(IW.GT.1) THEN
C           INTEGRAL TYPE (PS,SS) AND (PP,SS).
            DO 220 I=1,9
            WW(IP+I,JP) = PASS(I)
  220       CONTINUE
         ENDIF
         IF(IW.GT.1 .AND. JW.GT.1) THEN
C           INTEGRAL TYPE (PS,PS).
            WW(IP+1,JP+1) = PSPS(1)
            WW(IP+3,JP+1) = PSPS(2)
            WW(IP+6,JP+1) = PSPS(4)
            WW(IP+1,JP+3) = PSPS(2)
            WW(IP+3,JP+3) = PSPS(3)
            WW(IP+6,JP+3) = PSPS(5)
            WW(IP+1,JP+6) = PSPS(4)
            WW(IP+3,JP+6) = PSPS(5)
            WW(IP+6,JP+6) = PSPS(6)
C           INTEGRAL TYPE (PS,PP).
            DO 240 I=1,3
            IJ  = IP+IPP(I)
            DO 230 J=1,6
            WW(IJ,JP+ISS(J)) = PSPP(J,I)
  230       CONTINUE
  240       CONTINUE
C           INTEGRAL TYPES (PP,PS) AND (PP,PP).
            DO 270 I=1,6
            IJ  = IP+ISS(I)
            DO 250 J=1,3
            WW(IJ,JP+IPP(J)) = PPPS(J,I)
  250       CONTINUE
            DO 260 J=1,6
            WW(IJ,JP+ISS(J)) = PPPP(J,I)
  260       CONTINUE
  270       CONTINUE
         ENDIF
      ENDIF
      RETURN
      END
