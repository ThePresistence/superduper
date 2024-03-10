      SUBROUTINE SYMREF (XCALC,X,NSYM,N,ID4,ID5,ID5MOD)
C     *
C     FIND THE REQUESTED VALUE XCALC IN THE ARRAY X(N).
C     *
C     NOTATION.
C     XCALC    REQUESTED VALUE.
C     X(N)     ARRAY OF AVAILABLE VALUES.
C     NSYM(N)  ARRAY OF SYMMETRY LABELS (STANDARD CONVENTIONS).
C     N        NUMBER OF AVAILABLE DATA.
C     ID4      FIRST  IDENTIFIER (SEE INPUT DESCRIPTION).
C     ID5      SECOND IDENTIFIER (SEE INPUT DESCRIPTION).
C     ID5MOD   ASSIGNMENT VIA ID5: ORDER OF AVAILABLE VALUES.
C              =-1 DESCENDING ORDER (DEFAULT).
C              = 1 ASCENDING ORDER (DEFAULT).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),NSYM(N)
      XCALC = 0.D0
      IF(ID4.LE.0 .OR. ID4.GT.N) RETURN
      IF(ID5.EQ.0) THEN
         XCALC = X(ID4)
      ELSE IF(ID5.EQ.-1) THEN
         IF(ID5MOD.LE.0) THEN
            XCALC = X(N+1-ID4)
         ELSE
            XCALC = X(ID4)
         ENDIF
      ELSE IF(ID5.GT.0) THEN
         IF(ID5MOD.LE.0) THEN
            JA = N
            JB = 1
            JD =-1
         ELSE
            JA = 1
            JB = N
            JD = 1
         ENDIF
         NID5  = 0
         DO 10 J=JA,JB,JD
         IF(ID5.EQ.NSYM(J)) THEN
            NID5 = NID5+1
            IF(ID4.EQ.NID5) THEN
               XCALC = X(J)
               GO TO 20
            ENDIF
         ENDIF
   10    CONTINUE
   20    CONTINUE
      ENDIF
      RETURN
      END
