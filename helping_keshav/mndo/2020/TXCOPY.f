      SUBROUTINE TXCOPY (TEXT,TEXT1,N)
C     *
C     COPY TEXT TO TEXT1 WITHOUT LEADING BLANKS.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) TEXT,TEXT1
      CHARACTER*1  BLANK
      DATA BLANK/' '/
      TEXT1  = BLANK
      IBLANK = 0
      DO 10 I=1,N
      IF(TEXT(I:I).NE.BLANK) THEN
         IBLANK = I
         GO TO 20
      ENDIF
   10 CONTINUE
   20 CONTINUE
      IF(IBLANK.GT.0) TEXT1 = TEXT(IBLANK:)
      RETURN
      END
