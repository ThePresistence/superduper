      SUBROUTINE INEIGF (IMOPAC,JOP)
C     *
C     INPUT OF OPTIONS FOR EIGENVECTOR FOLLOWING.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./INOPT1/ IN1(300)
     ./INOPT2/ IN2(300)
     ./INOPT3/ XN3(50)
     ./INOPT4/ XN4(50)
     ./NBFILE/ NBF(20)
C *** FILE NUMBERS.
      NB5    = NBF(5)
C *** READ SPECIAL INPUT FOR EIGENVECTOR FOLLOWING OR HDLC OPTIMIZER.
      IF(IMOPAC.EQ.0) THEN
         IF(IN2(8).GT.2 .OR. IN2(8).LT.-2 .OR. IN2(61).EQ.1) THEN
            READ(NB5,500) (IN1(I),I=80,89),(XN3(J),J=1,6)
            XN3(7) = ZERO
         ELSE
            IN1(80:89) = 0
            XN3(1:7)   = ZERO
         ENDIF
      ENDIF
      IN2(80:89) = IN1(80:89)
      XN4(1:7)   = XN3(1:7)
C *** IMPOSE GENERAL DEFAULT VALUES FOR EIGENVECTOR FOLLOWING.
      IF(IN2(8).GT.0) THEN
         IF(XN4(1).EQ.ZERO) XN4(1)=0.2D0
         IF(XN4(2).EQ.ZERO) XN4(2)=0.001D0
         IF(XN4(5).EQ.ZERO) XN4(5)=4.0D0
         IF(XN4(6).EQ.ZERO) XN4(6)=0.8D0
         IF(XN4(7).EQ.ZERO) THEN
            XN4(7)=ONE
            IF(IN2(43).GT.1) XN4(7)=ONE/IN2(43)
            IF(IN2(89).GT.1) XN4(7)=ONE*IN2(89)
         ENDIF
         IF(JOP.EQ.0 .OR. JOP.EQ.3 .OR. JOP.EQ.5) THEN
            IF(IN2(80).LT.0) IN2(80)=0
            IF(IN2(82).EQ.0) IN2(82)=2
            IF(IN2(83).LT.0) IN2(83)=0
            IF(XN4(3).EQ.ZERO) XN4(3)=0.5D0
         ELSE IF(JOP.EQ.1 .OR. JOP.EQ.4 .OR. JOP.EQ.6) THEN
            IF(IN2(80).LE.0) IN2(80)=1
            IF(IN2(82).EQ.0) IN2(82)=1
            IF(IN2(83).LE.0) IN2(83)=1
            IF(XN4(3).EQ.ZERO) XN4(3)=0.3D0
         ENDIF
         IF(IN2(82).LT.0) IN2(82)=0
         IF(IN2(81).LE.0) IN2(81)=9999
C *** IMPOSE GENERAL DEFAULT VALUES FOR HDLC OPTIONS.
      ELSE IF(IN2(8).LT.0) THEN
         IF(XN4(1).EQ.ZERO) XN4(1)=0.02D0
         IF(XN4(2).EQ.ZERO) XN4(2)=0.0002D0
         IF(XN4(3).EQ.ZERO) XN4(3)=1.0D0
         IF(XN4(5).EQ.ZERO) XN4(5)=10.0D0
         IF(XN4(6).EQ.ZERO) XN4(6)=0.6D0
         IF(JOP.EQ.0 .OR. JOP.EQ.3 .OR. JOP.EQ.5) THEN
            IF(IN2(80).LT.0) IN2(80)=0
            IF(IN2(81).EQ.0) IN2(81)=-1
            IF(IN2(82).EQ.0) IN2(82)=2
            IF(IN2(83).LT.0) IN2(83)=0
         ELSE IF(JOP.EQ.1 .OR. JOP.EQ.4 .OR. JOP.EQ.6) THEN
            IF(IN2(80).LE.0) IN2(80)=1
            IF(IN2(81).LT.0) IN2(81)=0
            IF(IN2(82).EQ.0) IN2(82)=1
            IF(IN2(83).LE.0) IN2(83)=1
         ENDIF
         IF(IN2(82).LT.0) IN2(82)=0
      ENDIF
C *** IMPOSE NEWTON-RAPHSON DEFAULT VALUES.
      IF(IN2(8).EQ.2) THEN
         IN2(81) = 1
         IN2(83) = 1
         IN2(84) = 2
         IN2(86) = 0
         IN2(88) = 1
      ENDIF
      RETURN
  500 FORMAT(10I2,6F10.5)
      END
