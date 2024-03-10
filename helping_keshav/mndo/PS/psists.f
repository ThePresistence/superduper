C     ******************************************************************
      SUBROUTINE PSISTS(IORBB,AINT)
C
C  This subroutine is not intended to be human-readable.
C  It was generated automagically by Mathematica program
C  coded by: Serge Pachkovsky
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AINT(46,46)
      DIMENSION TEMP(2)
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
      IF( ICNT.GT.0 ) THEN
          IDST = IFIRST
          DO 20 I = 1, ICNT
              ISRC = IPERM(ICTR)
              ICTR  = ICTR + 1
              AINT(1,IDST) = AINT(1,ISRC)
              AINT(2,IDST) = AINT(2,ISRC)
              IDST  = ISRC
   20     CONTINUE
          ILAST = ISRC
      ELSE
          ILAST = IFIRST
      ENDIF
      DO 30 I = 1, 2
          AINT(I,ILAST) = TEMP(I)
   30 CONTINUE
      GOTO 10
      END
