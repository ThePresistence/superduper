      SUBROUTINE MMDPI
C     *
C     INITIALIZATION OF DISPERSION PARAMETERS.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./INOPT2/ IN2(300)
     ./OMXDC / VS6,VCD

      IOP=IN2(2)
      IMMDP=IN2(30)
      IF((IOP.EQ.-2.OR.IOP.EQ.-7).AND.IMMDP.GT.0) THEN
C     AM1-D AND PM3-D. PCCP, 2007, 9, 2362.
        VS6=1.4D0
        VCD=23.0D0
      ELSE IF(IMMDP.EQ.1) THEN
C     USING THE DAMPING FUNCTION IN JCP, 2001, 114, 5149 - 5155.
        IF    (IOP.EQ.-6) THEN
          VS6=0.8D0
          VCD=0.6D0
        ELSEIF(IOP.EQ.-8) THEN
          VS6=0.6D0
          VCD=1.4D0
        ENDIF
      ELSE IF(IMMDP.EQ.2) THEN
C     USING THE DAMPING FUNCTION IN JCC, 2006, 27, 1787 - 1799.
        IF    (IOP.EQ.-6) THEN
          VS6=0.5D0
          VCD=12.0D0
        ELSEIF(IOP.EQ.-8) THEN
          VS6=0.6D0
          VCD=18.6D0
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MMDPI
