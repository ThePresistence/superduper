C     ******************************************************************
      SUBROUTINE PSIR1S(IORBB,AINT,OUT)
C
C  This subroutine is not intended to be human-readable.
C  It was generated automagically by Mathematica program
C  coded by: Serge Pachkovsky
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AINT(461), OUT(46,46)
      RCO1(I0,IND) = AINT(I0)
      RCO2(I0,IND) = 0
      RSS1(I0,IND) = AINT(I0)
      RSS2(I0,IND) = 0
C
C   CORE and SS distributions on B
C
          OUT(1,1) = RCO1(1,1)
          OUT(1,2) = RCO1(16,1)
          OUT(2,1) = RSS1(2,1)
          OUT(2,2) = RSS1(17,1)
      IF( IORBB.LE.0 ) RETURN
C
C   SP and PP distributions on B
C
          OUT(1,5) = RCO1(51,1)
          OUT(1,6) = RCO1(66,1)
          OUT(1,8) = RCO1(87,1)
          OUT(1,11) = RCO1(124,1)
          OUT(1,3) = RCO2(0,1)
          OUT(1,4) = RCO2(0,1)
          OUT(1,7) = RCO2(0,1)
          OUT(1,9) = RCO2(0,1)
          OUT(1,10) = RCO2(0,1)
          OUT(2,5) = RSS1(52,1)
          OUT(2,6) = RSS1(67,1)
          OUT(2,8) = RSS1(88,1)
          OUT(2,11) = RSS1(125,1)
          OUT(2,3) = RSS2(0,1)
          OUT(2,4) = RSS2(0,1)
          OUT(2,7) = RSS2(0,1)
          OUT(2,9) = RSS2(0,1)
          OUT(2,10) = RSS2(0,1)
      IF( IORBB.LE.3 ) RETURN
C
C   SD, PD and DD distributions on B
C
          OUT(1,14) = RCO1(155,1)
          OUT(1,20) = RCO1(204,1)
          OUT(1,25) = RCO1(249,1)
          OUT(1,27) = RCO1(264,1)
          OUT(1,32) = RCO1(309,1)
          OUT(1,34) = RCO1(334,1)
          OUT(1,37) = RCO1(367,1)
          OUT(1,41) = RCO1(406,1)
          OUT(1,46) = RCO1(447,1)
          OUT(1,12) = RCO2(0,1)
          OUT(1,13) = RCO2(0,1)
          OUT(1,15) = RCO2(0,1)
          OUT(1,16) = RCO2(0,1)
          OUT(1,17) = RCO2(0,1)
          OUT(1,18) = RCO2(0,1)
          OUT(1,19) = RCO2(0,1)
          OUT(1,21) = RCO2(0,1)
          OUT(1,22) = RCO2(0,1)
          OUT(1,23) = RCO2(0,1)
          OUT(1,24) = RCO2(0,1)
          OUT(1,26) = RCO2(0,1)
          OUT(1,28) = RCO2(0,1)
          OUT(1,29) = RCO2(0,1)
          OUT(1,30) = RCO2(0,1)
          OUT(1,31) = RCO2(0,1)
          OUT(1,33) = RCO2(0,1)
          OUT(1,35) = RCO2(0,1)
          OUT(1,36) = RCO2(0,1)
          OUT(1,38) = RCO2(0,1)
          OUT(1,39) = RCO2(0,1)
          OUT(1,40) = RCO2(0,1)
          OUT(1,42) = RCO2(0,1)
          OUT(1,43) = RCO2(0,1)
          OUT(1,44) = RCO2(0,1)
          OUT(1,45) = RCO2(0,1)
          OUT(2,14) = RSS1(156,1)
          OUT(2,20) = RSS1(205,1)
          OUT(2,25) = RSS1(250,1)
          OUT(2,27) = RSS1(265,1)
          OUT(2,32) = RSS1(310,1)
          OUT(2,34) = RSS1(335,1)
          OUT(2,37) = RSS1(368,1)
          OUT(2,41) = RSS1(407,1)
          OUT(2,46) = RSS1(448,1)
          OUT(2,12) = RSS2(0,1)
          OUT(2,13) = RSS2(0,1)
          OUT(2,15) = RSS2(0,1)
          OUT(2,16) = RSS2(0,1)
          OUT(2,17) = RSS2(0,1)
          OUT(2,18) = RSS2(0,1)
          OUT(2,19) = RSS2(0,1)
          OUT(2,21) = RSS2(0,1)
          OUT(2,22) = RSS2(0,1)
          OUT(2,23) = RSS2(0,1)
          OUT(2,24) = RSS2(0,1)
          OUT(2,26) = RSS2(0,1)
          OUT(2,28) = RSS2(0,1)
          OUT(2,29) = RSS2(0,1)
          OUT(2,30) = RSS2(0,1)
          OUT(2,31) = RSS2(0,1)
          OUT(2,33) = RSS2(0,1)
          OUT(2,35) = RSS2(0,1)
          OUT(2,36) = RSS2(0,1)
          OUT(2,38) = RSS2(0,1)
          OUT(2,39) = RSS2(0,1)
          OUT(2,40) = RSS2(0,1)
          OUT(2,42) = RSS2(0,1)
          OUT(2,43) = RSS2(0,1)
          OUT(2,44) = RSS2(0,1)
          OUT(2,45) = RSS2(0,1)
      RETURN
      END
