C     ******************************************************************
      SUBROUTINE PSISTD(IORBB,AINT)
C
C  This subroutine is not intended to be human-readable.
C  It was generated automagically by Mathematica program
C  coded by: Serge Pachkovsky
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AINT(46,46)
      DIMENSION TEMP(46)
      COMMON /PSISTC/ IPERM(68)
      EXTERNAL PSISTB
C
      ICTR   = 1
   10 ICNT   = IPERM(ICTR)
      IFIRST = IPERM(ICTR+1)
      IF( ICNT.EQ.-1 ) THEN
          IF( IORBB.LE.IFIRST ) RETURN
          ICTR = ICTR + 2
          GOTO 10
      ENDIF
      ICTR   = ICTR + 2
      TEMP(1) = AINT(1,IFIRST)
      TEMP(2) = AINT(2,IFIRST)
      TEMP(3) = AINT(3,IFIRST)
      TEMP(5) = AINT(4,IFIRST)
      TEMP(8) = AINT(5,IFIRST)
      TEMP(4) = AINT(6,IFIRST)
      TEMP(6) = AINT(7,IFIRST)
      TEMP(7) = AINT(8,IFIRST)
      TEMP(9) = AINT(9,IFIRST)
      TEMP(10) = AINT(10,IFIRST)
      TEMP(11) = AINT(11,IFIRST)
      TEMP(12) = AINT(12,IFIRST)
      TEMP(17) = AINT(13,IFIRST)
      TEMP(23) = AINT(14,IFIRST)
      TEMP(30) = AINT(15,IFIRST)
      TEMP(38) = AINT(16,IFIRST)
      TEMP(13) = AINT(17,IFIRST)
      TEMP(14) = AINT(18,IFIRST)
      TEMP(15) = AINT(19,IFIRST)
      TEMP(18) = AINT(20,IFIRST)
      TEMP(19) = AINT(21,IFIRST)
      TEMP(20) = AINT(22,IFIRST)
      TEMP(24) = AINT(23,IFIRST)
      TEMP(25) = AINT(24,IFIRST)
      TEMP(26) = AINT(25,IFIRST)
      TEMP(31) = AINT(26,IFIRST)
      TEMP(32) = AINT(27,IFIRST)
      TEMP(33) = AINT(28,IFIRST)
      TEMP(39) = AINT(29,IFIRST)
      TEMP(40) = AINT(30,IFIRST)
      TEMP(41) = AINT(31,IFIRST)
      TEMP(16) = AINT(32,IFIRST)
      TEMP(21) = AINT(33,IFIRST)
      TEMP(22) = AINT(34,IFIRST)
      TEMP(27) = AINT(35,IFIRST)
      TEMP(28) = AINT(36,IFIRST)
      TEMP(29) = AINT(37,IFIRST)
      TEMP(34) = AINT(38,IFIRST)
      TEMP(35) = AINT(39,IFIRST)
      TEMP(36) = AINT(40,IFIRST)
      TEMP(37) = AINT(41,IFIRST)
      TEMP(42) = AINT(42,IFIRST)
      TEMP(43) = AINT(43,IFIRST)
      TEMP(44) = AINT(44,IFIRST)
      TEMP(45) = AINT(45,IFIRST)
      TEMP(46) = AINT(46,IFIRST)
      IF( ICNT.GT.0 ) THEN
          IDST = IFIRST
          DO 20 I = 1, ICNT
              ISRC = IPERM(ICTR)
              ICTR  = ICTR + 1
              AINT(1,IDST) = AINT(1,ISRC)
              AINT(2,IDST) = AINT(2,ISRC)
              AINT(3,IDST) = AINT(3,ISRC)
              AINT(5,IDST) = AINT(4,ISRC)
              AINT(8,IDST) = AINT(5,ISRC)
              AINT(4,IDST) = AINT(6,ISRC)
              AINT(6,IDST) = AINT(7,ISRC)
              AINT(7,IDST) = AINT(8,ISRC)
              AINT(9,IDST) = AINT(9,ISRC)
              AINT(10,IDST) = AINT(10,ISRC)
              AINT(11,IDST) = AINT(11,ISRC)
              AINT(12,IDST) = AINT(12,ISRC)
              AINT(17,IDST) = AINT(13,ISRC)
              AINT(23,IDST) = AINT(14,ISRC)
              AINT(30,IDST) = AINT(15,ISRC)
              AINT(38,IDST) = AINT(16,ISRC)
              AINT(13,IDST) = AINT(17,ISRC)
              AINT(14,IDST) = AINT(18,ISRC)
              AINT(15,IDST) = AINT(19,ISRC)
              AINT(18,IDST) = AINT(20,ISRC)
              AINT(19,IDST) = AINT(21,ISRC)
              AINT(20,IDST) = AINT(22,ISRC)
              AINT(24,IDST) = AINT(23,ISRC)
              AINT(25,IDST) = AINT(24,ISRC)
              AINT(26,IDST) = AINT(25,ISRC)
              AINT(31,IDST) = AINT(26,ISRC)
              AINT(32,IDST) = AINT(27,ISRC)
              AINT(33,IDST) = AINT(28,ISRC)
              AINT(39,IDST) = AINT(29,ISRC)
              AINT(40,IDST) = AINT(30,ISRC)
              AINT(41,IDST) = AINT(31,ISRC)
              AINT(16,IDST) = AINT(32,ISRC)
              AINT(21,IDST) = AINT(33,ISRC)
              AINT(22,IDST) = AINT(34,ISRC)
              AINT(27,IDST) = AINT(35,ISRC)
              AINT(28,IDST) = AINT(36,ISRC)
              AINT(29,IDST) = AINT(37,ISRC)
              AINT(34,IDST) = AINT(38,ISRC)
              AINT(35,IDST) = AINT(39,ISRC)
              AINT(36,IDST) = AINT(40,ISRC)
              AINT(37,IDST) = AINT(41,ISRC)
              AINT(42,IDST) = AINT(42,ISRC)
              AINT(43,IDST) = AINT(43,ISRC)
              AINT(44,IDST) = AINT(44,ISRC)
              AINT(45,IDST) = AINT(45,ISRC)
              AINT(46,IDST) = AINT(46,ISRC)
              IDST  = ISRC
   20     CONTINUE
          ILAST = ISRC
      ELSE
          ILAST = IFIRST
      ENDIF
      DO 30 I = 1, 46
          AINT(I,ILAST) = TEMP(I)
   30 CONTINUE
      GOTO 10
      END
