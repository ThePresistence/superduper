C     ******************************************************************
      SUBROUTINE PSISTP(IORBB,AINT)
C
C  This subroutine is not intended to be human-readable.
C  It was generated automagically by Mathematica program
C  coded by: Serge Pachkovsky
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AINT(46,46)
      DIMENSION TEMP(11)
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
              IDST  = ISRC
   20     CONTINUE
          ILAST = ISRC
      ELSE
          ILAST = IFIRST
      ENDIF
      DO 30 I = 1, 11
          AINT(I,ILAST) = TEMP(I)
   30 CONTINUE
      GOTO 10
      END
