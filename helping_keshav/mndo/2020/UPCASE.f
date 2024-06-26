      SUBROUTINE UPPCAS(KEYWRD,LEN)
C     *
C     CONVERT LOWER CASE TO UPPER CASE IN A CHARACTER STRING.
C     ADAPTED FROM MOPAC(6.0) WRITTEN BY J.J.P.STEWART.
C     *
      CHARACTER*(*) KEYWRD
      ICAPA  = ICHAR('A')
      ILOWA  = ICHAR('a')
      ILOWZ  = ICHAR('z')
      DO 10 I=1,LEN
         ILINE=ICHAR(KEYWRD(I:I))
         IF(ILINE.GE.ILOWA.AND.ILINE.LE.ILOWZ) THEN
            KEYWRD(I:I)=CHAR(ILINE+ICAPA-ILOWA)
         ENDIF
   10 CONTINUE
      RETURN
      END
