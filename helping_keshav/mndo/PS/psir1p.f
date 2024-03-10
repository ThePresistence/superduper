C     ******************************************************************
      SUBROUTINE PSIR1P(IORBB,AINT,OUT,RMSP,RMPP)
C
C  This subroutine is not intended to be human-readable.
C  It was generated automagically by Mathematica program
C  coded by: Serge Pachkovsky
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AINT(461), OUT(46,46)
      DIMENSION RMSP(3,3), RMPP(6,6)
C     Next four statements commented out by WT: Not used.
C     RCO1(I0,IND) = AINT(I0)
C     RCO2(I0,IND) = 0
C     RSS1(I0,IND) = AINT(I0)
C     RSS2(I0,IND) = 0
      RSP1(I0,IND) = AINT(I0)*RMSP(3,IND)
      RSP2(I0,IND) = AINT(I0)*RMSP(1,IND)
      RSP3(I0,IND) = AINT(I0)*RMSP(2,IND)
      RSP4(I0,IND) = 0
      RPP1(I0,IND) = AINT(I0)*RMPP(1,IND) + AINT(1 + I0)*RMPP(3,IND) + 
     1 AINT(2 + I0)*RMPP(6,IND)
      RPP2(I0,IND) = AINT(I0)*RMPP(4,IND)
      RPP3(I0,IND) = AINT(I0)*RMPP(5,IND)
      RPP4(I0,IND) = AINT(I0)*RMPP(2,IND)
      RPP5(I0,IND) = AINT(I0)*RMPP(1,IND) + AINT(1 + I0)*RMPP(3,IND)
      RPP6(I0,IND) = 0
C
C   CORE and SS distributions on B
C
      DO 1000 II=1,3
          OUT(2 + II,1) = RSP1(3,II)
          OUT(2 + II,2) = RSP1(18,II)
 1000 CONTINUE
      DO 1010 II=1,6
          OUT(5 + II,1) = RPP1(4,II)
          OUT(5 + II,2) = RPP1(19,II)
 1010 CONTINUE
      IF( IORBB.LE.0 ) RETURN
C
C   SP and PP distributions on B
C
      DO 1020 II=1,3
          OUT(2 + II,5) = RSP1(53,II)
          OUT(2 + II,6) = RSP1(68,II)
          OUT(2 + II,8) = RSP1(89,II)
          OUT(2 + II,11) = RSP1(126,II)
 1020 CONTINUE
      DO 1030 II=1,3
          OUT(2 + II,3) = RSP2(31,II)
          OUT(2 + II,9) = RSP2(104,II)
 1030 CONTINUE
      DO 1040 II=1,3
          OUT(2 + II,4) = RSP3(41,II)
          OUT(2 + II,10) = RSP3(114,II)
 1040 CONTINUE
      DO 1050 II=1,3
          OUT(2 + II,7) = RSP4(0,II)
 1050 CONTINUE
      DO 1060 II=1,6
          OUT(5 + II,5) = RPP1(54,II)
          OUT(5 + II,6) = RPP1(69,II)
          OUT(5 + II,8) = RPP1(90,II)
          OUT(5 + II,11) = RPP1(127,II)
 1060 CONTINUE
      DO 1070 II=1,6
          OUT(5 + II,3) = RPP2(32,II)
          OUT(5 + II,9) = RPP2(105,II)
 1070 CONTINUE
      DO 1080 II=1,6
          OUT(5 + II,4) = RPP3(42,II)
          OUT(5 + II,10) = RPP3(115,II)
 1080 CONTINUE
      DO 1090 II=1,6
          OUT(5 + II,7) = RPP4(83,II)
 1090 CONTINUE
      IF( IORBB.LE.3 ) RETURN
C
C   SD, PD and DD distributions on B
C
      DO 1100 II=1,3
          OUT(2 + II,14) = RSP1(157,II)
          OUT(2 + II,20) = RSP1(206,II)
          OUT(2 + II,25) = RSP1(251,II)
          OUT(2 + II,27) = RSP1(266,II)
          OUT(2 + II,32) = RSP1(311,II)
          OUT(2 + II,34) = RSP1(336,II)
          OUT(2 + II,37) = RSP1(369,II)
          OUT(2 + II,41) = RSP1(408,II)
          OUT(2 + II,46) = RSP1(449,II)
 1100 CONTINUE
      DO 1110 II=1,3
          OUT(2 + II,13) = RSP2(145,II)
          OUT(2 + II,17) = RSP2(184,II)
          OUT(2 + II,22) = RSP2(219,II)
          OUT(2 + II,23) = RSP2(229,II)
          OUT(2 + II,30) = RSP2(299,II)
          OUT(2 + II,33) = RSP2(324,II)
          OUT(2 + II,36) = RSP2(357,II)
          OUT(2 + II,45) = RSP2(437,II)
 1110 CONTINUE
      DO 1120 II=1,3
          OUT(2 + II,15) = RSP3(170,II)
          OUT(2 + II,18) = RSP3(194,II)
          OUT(2 + II,24) = RSP3(239,II)
          OUT(2 + II,28) = RSP3(279,II)
          OUT(2 + II,29) = RSP3(289,II)
          OUT(2 + II,38) = RSP3(382,II)
          OUT(2 + II,40) = RSP3(396,II)
          OUT(2 + II,43) = RSP3(423,II)
 1120 CONTINUE
      DO 1130 II=1,3
          OUT(2 + II,12) = RSP4(0,II)
          OUT(2 + II,16) = RSP4(0,II)
          OUT(2 + II,19) = RSP4(0,II)
          OUT(2 + II,21) = RSP4(0,II)
          OUT(2 + II,26) = RSP4(0,II)
          OUT(2 + II,31) = RSP4(0,II)
          OUT(2 + II,35) = RSP4(0,II)
          OUT(2 + II,39) = RSP4(0,II)
          OUT(2 + II,42) = RSP4(0,II)
          OUT(2 + II,44) = RSP4(0,II)
 1130 CONTINUE
      DO 1140 II=1,6
          OUT(5 + II,14) = RPP1(158,II)
          OUT(5 + II,20) = RPP1(207,II)
          OUT(5 + II,25) = RPP1(252,II)
          OUT(5 + II,27) = RPP1(267,II)
          OUT(5 + II,32) = RPP1(312,II)
          OUT(5 + II,34) = RPP1(337,II)
          OUT(5 + II,37) = RPP1(370,II)
          OUT(5 + II,41) = RPP1(409,II)
          OUT(5 + II,46) = RPP1(450,II)
 1140 CONTINUE
      DO 1150 II=1,6
          OUT(5 + II,13) = RPP2(146,II)
          OUT(5 + II,17) = RPP2(185,II)
          OUT(5 + II,22) = RPP2(220,II)
          OUT(5 + II,23) = RPP2(230,II)
          OUT(5 + II,30) = RPP2(300,II)
          OUT(5 + II,33) = RPP2(325,II)
          OUT(5 + II,36) = RPP2(358,II)
          OUT(5 + II,45) = RPP2(438,II)
 1150 CONTINUE
      DO 1160 II=1,6
          OUT(5 + II,15) = RPP3(171,II)
          OUT(5 + II,18) = RPP3(195,II)
          OUT(5 + II,24) = RPP3(240,II)
          OUT(5 + II,28) = RPP3(280,II)
          OUT(5 + II,29) = RPP3(290,II)
          OUT(5 + II,38) = RPP3(383,II)
          OUT(5 + II,40) = RPP3(397,II)
          OUT(5 + II,43) = RPP3(424,II)
 1160 CONTINUE
      DO 1170 II=1,6
          OUT(5 + II,16) = RPP4(180,II)
          OUT(5 + II,39) = RPP4(392,II)
          OUT(5 + II,44) = RPP4(433,II)
 1170 CONTINUE
      DO 1180 II=1,6
          OUT(5 + II,12) = RPP5(139,II)
          OUT(5 + II,35) = RPP5(351,II)
 1180 CONTINUE
      DO 1190 II=1,6
          OUT(5 + II,19) = RPP6(0,II)
          OUT(5 + II,21) = RPP6(0,II)
          OUT(5 + II,26) = RPP6(0,II)
          OUT(5 + II,31) = RPP6(0,II)
          OUT(5 + II,42) = RPP6(0,II)
 1190 CONTINUE
      RETURN
      END
