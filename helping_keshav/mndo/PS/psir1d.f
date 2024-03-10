C     ******************************************************************
      SUBROUTINE PSIR1D(IORBB,AINT,OUT,RMSP,RMPP,RMSD,RMPD,RMDD)
C
C  This subroutine is not intended to be human-readable.
C  It was generated automagically by Mathematica program
C  coded by Serge Pachkovsky
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AINT(461), OUT(46,46)
      DIMENSION RMSP(3,3), RMPP(6,6)
      DIMENSION RMSD(5,5), RMPD(15,15), RMDD(15,15)
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
      RSD1(I0,IND) = AINT(I0)*RMSD(3,IND)
      RSD2(I0,IND) = AINT(I0)*RMSD(2,IND)
      RSD3(I0,IND) = AINT(I0)*RMSD(4,IND)
      RSD4(I0,IND) = AINT(I0)*RMSD(1,IND) + AINT(1 + I0)*RMSD(3,IND)
      RSD5(I0,IND) = AINT(I0)*RMSD(5,IND)
      RSD6(I0,IND) = AINT(I0)*RMSD(1,IND)
      RSD7(I0,IND) = 0
      RPD1(I0,IND) = AINT(I0)*RMPD(4,IND) + AINT(1 + I0)*RMPD(9,IND) + 
     1 AINT(2 + I0)*RMPD(11,IND)
      RPD2(I0,IND) = AINT(I0)*RMPD(1,IND) + AINT(1 + I0)*RMPD(6,IND) + 
     1 AINT(2 + I0)*RMPD(7,IND) + AINT(3 + I0)*RMPD(14,IND)
      RPD3(I0,IND) = AINT(I0)*RMPD(2,IND) + AINT(1 + I0)*RMPD(8,IND) + 
     1 AINT(2 + I0)*RMPD(12,IND) + AINT(3 + I0)*RMPD(13,IND)
      RPD4(I0,IND) = 0
      RDD1(I0,IND) = AINT(I0)*RMDD(1,IND) + AINT(1 + I0)*RMDD(3,IND) + 
     1 AINT(2 + I0)*RMDD(6,IND) + AINT(3 + I0)*RMDD(10,IND) + 
     2 AINT(4 + I0)*RMDD(15,IND)
      RDD2(I0,IND) = AINT(I0)*RMDD(2,IND) + AINT(1 + I0)*RMDD(5,IND) + 
     1 AINT(2 + I0)*RMDD(14,IND)
      RDD3(I0,IND) = AINT(I0)*RMDD(7,IND) + AINT(1 + I0)*RMDD(9,IND) + 
     1 AINT(2 + I0)*RMDD(12,IND)
      RDD4(I0,IND) = AINT(I0)*RMDD(1,IND) + AINT(1 + I0)*RMDD(3,IND) + 
     1 AINT(2 + I0)*RMDD(4,IND) + AINT(3 + I0)*RMDD(6,IND) + 
     2 AINT(4 + I0)*RMDD(10,IND) + AINT(5 + I0)*RMDD(15,IND)
      RDD5(I0,IND) = AINT(I0)*RMDD(8,IND) + AINT(1 + I0)*RMDD(13,IND)
      RDD6(I0,IND) = AINT(I0)*RMDD(3,IND) + AINT(1 + I0)*RMDD(4,IND) + 
     1 AINT(2 + I0)*RMDD(10,IND)
      RDD7(I0,IND) = 0
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
      DO 1020 II=1,5
          OUT(11 + II,1) = RSD1(7,II)
          OUT(11 + II,2) = RSD1(22,II)
 1020 CONTINUE
      DO 1030 II=1,15
          OUT(16 + II,1) = RPD1(8,II)
          OUT(16 + II,2) = RPD1(23,II)
 1030 CONTINUE
      DO 1040 II=1,15
          OUT(31 + II,1) = RDD1(11,II)
          OUT(31 + II,2) = RDD1(26,II)
 1040 CONTINUE
      IF( IORBB.LE.0 ) RETURN
C
C   SP and PP distributions on B
C
      DO 1050 II=1,3
          OUT(2 + II,5) = RSP1(53,II)
          OUT(2 + II,6) = RSP1(68,II)
          OUT(2 + II,8) = RSP1(89,II)
          OUT(2 + II,11) = RSP1(126,II)
 1050 CONTINUE
      DO 1060 II=1,3
          OUT(2 + II,3) = RSP2(31,II)
          OUT(2 + II,9) = RSP2(104,II)
 1060 CONTINUE
      DO 1070 II=1,3
          OUT(2 + II,4) = RSP3(41,II)
          OUT(2 + II,10) = RSP3(114,II)
 1070 CONTINUE
      DO 1080 II=1,3
          OUT(2 + II,7) = RSP4(0,II)
 1080 CONTINUE
      DO 1090 II=1,6
          OUT(5 + II,5) = RPP1(54,II)
          OUT(5 + II,6) = RPP1(69,II)
          OUT(5 + II,8) = RPP1(90,II)
          OUT(5 + II,11) = RPP1(127,II)
 1090 CONTINUE
      DO 1100 II=1,6
          OUT(5 + II,3) = RPP2(32,II)
          OUT(5 + II,9) = RPP2(105,II)
 1100 CONTINUE
      DO 1110 II=1,6
          OUT(5 + II,4) = RPP3(42,II)
          OUT(5 + II,10) = RPP3(115,II)
 1110 CONTINUE
      DO 1120 II=1,6
          OUT(5 + II,7) = RPP4(83,II)
 1120 CONTINUE
      DO 1130 II=1,5
          OUT(11 + II,5) = RSD1(57,II)
          OUT(11 + II,11) = RSD1(130,II)
 1130 CONTINUE
      DO 1140 II=1,5
          OUT(11 + II,3) = RSD2(33,II)
          OUT(11 + II,9) = RSD2(106,II)
 1140 CONTINUE
      DO 1150 II=1,5
          OUT(11 + II,4) = RSD3(43,II)
          OUT(11 + II,10) = RSD3(116,II)
 1150 CONTINUE
      DO 1160 II=1,5
          OUT(11 + II,6) = RSD4(72,II)
          OUT(11 + II,8) = RSD4(93,II)
 1160 CONTINUE
      DO 1170 II=1,5
          OUT(11 + II,7) = RSD5(84,II)
 1170 CONTINUE
      DO 1180 II=1,15
          OUT(16 + II,5) = RPD1(58,II)
          OUT(16 + II,6) = RPD1(74,II)
          OUT(16 + II,8) = RPD1(95,II)
          OUT(16 + II,11) = RPD1(131,II)
 1180 CONTINUE
      DO 1190 II=1,15
          OUT(16 + II,3) = RPD2(34,II)
          OUT(16 + II,9) = RPD2(107,II)
 1190 CONTINUE
      DO 1200 II=1,15
          OUT(16 + II,4) = RPD3(44,II)
          OUT(16 + II,10) = RPD3(117,II)
 1200 CONTINUE
      DO 1210 II=1,15
          OUT(16 + II,7) = RPD4(0,II)
 1210 CONTINUE
      DO 1220 II=1,15
          OUT(31 + II,5) = RDD1(61,II)
          OUT(31 + II,11) = RDD1(134,II)
 1220 CONTINUE
      DO 1230 II=1,15
          OUT(31 + II,3) = RDD2(38,II)
          OUT(31 + II,9) = RDD2(111,II)
 1230 CONTINUE
      DO 1240 II=1,15
          OUT(31 + II,4) = RDD3(48,II)
          OUT(31 + II,10) = RDD3(121,II)
 1240 CONTINUE
      DO 1250 II=1,15
          OUT(31 + II,6) = RDD4(77,II)
          OUT(31 + II,8) = RDD4(98,II)
 1250 CONTINUE
      DO 1260 II=1,15
          OUT(31 + II,7) = RDD5(85,II)
 1260 CONTINUE
      IF( IORBB.LE.3 ) RETURN
C
C   SD, PD and DD distributions on B
C
      DO 1270 II=1,3
          OUT(2 + II,14) = RSP1(157,II)
          OUT(2 + II,20) = RSP1(206,II)
          OUT(2 + II,25) = RSP1(251,II)
          OUT(2 + II,27) = RSP1(266,II)
          OUT(2 + II,32) = RSP1(311,II)
          OUT(2 + II,34) = RSP1(336,II)
          OUT(2 + II,37) = RSP1(369,II)
          OUT(2 + II,41) = RSP1(408,II)
          OUT(2 + II,46) = RSP1(449,II)
 1270 CONTINUE
      DO 1280 II=1,3
          OUT(2 + II,13) = RSP2(145,II)
          OUT(2 + II,17) = RSP2(184,II)
          OUT(2 + II,22) = RSP2(219,II)
          OUT(2 + II,23) = RSP2(229,II)
          OUT(2 + II,30) = RSP2(299,II)
          OUT(2 + II,33) = RSP2(324,II)
          OUT(2 + II,36) = RSP2(357,II)
          OUT(2 + II,45) = RSP2(437,II)
 1280 CONTINUE
      DO 1290 II=1,3
          OUT(2 + II,15) = RSP3(170,II)
          OUT(2 + II,18) = RSP3(194,II)
          OUT(2 + II,24) = RSP3(239,II)
          OUT(2 + II,28) = RSP3(279,II)
          OUT(2 + II,29) = RSP3(289,II)
          OUT(2 + II,38) = RSP3(382,II)
          OUT(2 + II,40) = RSP3(396,II)
          OUT(2 + II,43) = RSP3(423,II)
 1290 CONTINUE
      DO 1300 II=1,3
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
 1300 CONTINUE
      DO 1310 II=1,6
          OUT(5 + II,14) = RPP1(158,II)
          OUT(5 + II,20) = RPP1(207,II)
          OUT(5 + II,25) = RPP1(252,II)
          OUT(5 + II,27) = RPP1(267,II)
          OUT(5 + II,32) = RPP1(312,II)
          OUT(5 + II,34) = RPP1(337,II)
          OUT(5 + II,37) = RPP1(370,II)
          OUT(5 + II,41) = RPP1(409,II)
          OUT(5 + II,46) = RPP1(450,II)
 1310 CONTINUE
      DO 1320 II=1,6
          OUT(5 + II,13) = RPP2(146,II)
          OUT(5 + II,17) = RPP2(185,II)
          OUT(5 + II,22) = RPP2(220,II)
          OUT(5 + II,23) = RPP2(230,II)
          OUT(5 + II,30) = RPP2(300,II)
          OUT(5 + II,33) = RPP2(325,II)
          OUT(5 + II,36) = RPP2(358,II)
          OUT(5 + II,45) = RPP2(438,II)
 1320 CONTINUE
      DO 1330 II=1,6
          OUT(5 + II,15) = RPP3(171,II)
          OUT(5 + II,18) = RPP3(195,II)
          OUT(5 + II,24) = RPP3(240,II)
          OUT(5 + II,28) = RPP3(280,II)
          OUT(5 + II,29) = RPP3(290,II)
          OUT(5 + II,38) = RPP3(383,II)
          OUT(5 + II,40) = RPP3(397,II)
          OUT(5 + II,43) = RPP3(424,II)
 1330 CONTINUE
      DO 1340 II=1,6
          OUT(5 + II,16) = RPP4(180,II)
          OUT(5 + II,39) = RPP4(392,II)
          OUT(5 + II,44) = RPP4(433,II)
 1340 CONTINUE
      DO 1350 II=1,6
          OUT(5 + II,12) = RPP5(139,II)
          OUT(5 + II,35) = RPP5(351,II)
 1350 CONTINUE
      DO 1360 II=1,6
          OUT(5 + II,19) = RPP6(0,II)
          OUT(5 + II,21) = RPP6(0,II)
          OUT(5 + II,26) = RPP6(0,II)
          OUT(5 + II,31) = RPP6(0,II)
          OUT(5 + II,42) = RPP6(0,II)
 1360 CONTINUE
      DO 1370 II=1,5
          OUT(11 + II,14) = RSD1(161,II)
          OUT(11 + II,20) = RSD1(210,II)
          OUT(11 + II,25) = RSD1(255,II)
          OUT(11 + II,27) = RSD1(270,II)
          OUT(11 + II,32) = RSD1(315,II)
          OUT(11 + II,37) = RSD1(373,II)
          OUT(11 + II,46) = RSD1(453,II)
 1370 CONTINUE
      DO 1380 II=1,5
          OUT(11 + II,13) = RSD2(147,II)
          OUT(11 + II,17) = RSD2(186,II)
          OUT(11 + II,22) = RSD2(221,II)
          OUT(11 + II,23) = RSD2(231,II)
          OUT(11 + II,30) = RSD2(301,II)
          OUT(11 + II,33) = RSD2(326,II)
          OUT(11 + II,36) = RSD2(359,II)
          OUT(11 + II,45) = RSD2(439,II)
 1380 CONTINUE
      DO 1390 II=1,5
          OUT(11 + II,15) = RSD3(172,II)
          OUT(11 + II,18) = RSD3(196,II)
          OUT(11 + II,24) = RSD3(241,II)
          OUT(11 + II,28) = RSD3(281,II)
          OUT(11 + II,29) = RSD3(291,II)
          OUT(11 + II,38) = RSD3(384,II)
          OUT(11 + II,40) = RSD3(398,II)
          OUT(11 + II,43) = RSD3(425,II)
 1390 CONTINUE
      DO 1400 II=1,5
          OUT(11 + II,34) = RSD4(340,II)
          OUT(11 + II,41) = RSD4(412,II)
 1400 CONTINUE
      DO 1410 II=1,5
          OUT(11 + II,16) = RSD5(181,II)
          OUT(11 + II,39) = RSD5(393,II)
          OUT(11 + II,44) = RSD5(434,II)
 1410 CONTINUE
      DO 1420 II=1,5
          OUT(11 + II,12) = RSD6(141,II)
          OUT(11 + II,35) = RSD6(353,II)
 1420 CONTINUE
      DO 1430 II=1,5
          OUT(11 + II,19) = RSD7(0,II)
          OUT(11 + II,21) = RSD7(0,II)
          OUT(11 + II,26) = RSD7(0,II)
          OUT(11 + II,31) = RSD7(0,II)
          OUT(11 + II,42) = RSD7(0,II)
 1430 CONTINUE
      DO 1440 II=1,15
          OUT(16 + II,14) = RPD1(162,II)
          OUT(16 + II,20) = RPD1(211,II)
          OUT(16 + II,25) = RPD1(256,II)
          OUT(16 + II,27) = RPD1(271,II)
          OUT(16 + II,32) = RPD1(316,II)
          OUT(16 + II,34) = RPD1(342,II)
          OUT(16 + II,37) = RPD1(374,II)
          OUT(16 + II,41) = RPD1(414,II)
          OUT(16 + II,46) = RPD1(454,II)
 1440 CONTINUE
      DO 1450 II=1,15
          OUT(16 + II,13) = RPD2(148,II)
          OUT(16 + II,17) = RPD2(187,II)
          OUT(16 + II,22) = RPD2(222,II)
          OUT(16 + II,23) = RPD2(232,II)
          OUT(16 + II,30) = RPD2(302,II)
          OUT(16 + II,33) = RPD2(327,II)
          OUT(16 + II,36) = RPD2(360,II)
          OUT(16 + II,45) = RPD2(440,II)
 1450 CONTINUE
      DO 1460 II=1,15
          OUT(16 + II,15) = RPD3(173,II)
          OUT(16 + II,18) = RPD3(197,II)
          OUT(16 + II,24) = RPD3(242,II)
          OUT(16 + II,28) = RPD3(282,II)
          OUT(16 + II,29) = RPD3(292,II)
          OUT(16 + II,38) = RPD3(385,II)
          OUT(16 + II,40) = RPD3(399,II)
          OUT(16 + II,43) = RPD3(426,II)
 1460 CONTINUE
      DO 1470 II=1,15
          OUT(16 + II,12) = RPD4(0,II)
          OUT(16 + II,16) = RPD4(0,II)
          OUT(16 + II,19) = RPD4(0,II)
          OUT(16 + II,21) = RPD4(0,II)
          OUT(16 + II,26) = RPD4(0,II)
          OUT(16 + II,31) = RPD4(0,II)
          OUT(16 + II,35) = RPD4(0,II)
          OUT(16 + II,39) = RPD4(0,II)
          OUT(16 + II,42) = RPD4(0,II)
          OUT(16 + II,44) = RPD4(0,II)
 1470 CONTINUE
      DO 1480 II=1,15
          OUT(31 + II,14) = RDD1(165,II)
          OUT(31 + II,20) = RDD1(214,II)
          OUT(31 + II,25) = RDD1(259,II)
          OUT(31 + II,27) = RDD1(274,II)
          OUT(31 + II,32) = RDD1(319,II)
          OUT(31 + II,37) = RDD1(377,II)
          OUT(31 + II,46) = RDD1(457,II)
 1480 CONTINUE
      DO 1490 II=1,15
          OUT(31 + II,13) = RDD2(152,II)
          OUT(31 + II,17) = RDD2(191,II)
          OUT(31 + II,22) = RDD2(226,II)
          OUT(31 + II,23) = RDD2(236,II)
          OUT(31 + II,30) = RDD2(306,II)
          OUT(31 + II,33) = RDD2(331,II)
          OUT(31 + II,36) = RDD2(364,II)
          OUT(31 + II,45) = RDD2(444,II)
 1490 CONTINUE
      DO 1500 II=1,15
          OUT(31 + II,15) = RDD3(177,II)
          OUT(31 + II,18) = RDD3(201,II)
          OUT(31 + II,24) = RDD3(246,II)
          OUT(31 + II,28) = RDD3(286,II)
          OUT(31 + II,29) = RDD3(296,II)
          OUT(31 + II,38) = RDD3(389,II)
          OUT(31 + II,40) = RDD3(403,II)
          OUT(31 + II,43) = RDD3(430,II)
 1500 CONTINUE
      DO 1510 II=1,15
          OUT(31 + II,34) = RDD4(345,II)
          OUT(31 + II,41) = RDD4(417,II)
 1510 CONTINUE
      DO 1520 II=1,15
          OUT(31 + II,16) = RDD5(182,II)
          OUT(31 + II,39) = RDD5(394,II)
          OUT(31 + II,44) = RDD5(435,II)
 1520 CONTINUE
      DO 1530 II=1,15
          OUT(31 + II,12) = RDD6(142,II)
          OUT(31 + II,35) = RDD6(354,II)
 1530 CONTINUE
      DO 1540 II=1,15
          OUT(31 + II,19) = RDD7(0,II)
          OUT(31 + II,21) = RDD7(0,II)
          OUT(31 + II,26) = RDD7(0,II)
          OUT(31 + II,31) = RDD7(0,II)
          OUT(31 + II,42) = RDD7(0,II)
 1540 CONTINUE
      RETURN
      END
