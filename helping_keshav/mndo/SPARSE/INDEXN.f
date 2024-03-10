C     ******************************************************************
      SUBROUTINE INDEXN (N,IA,IND)
C     *
C     INDEXES IA(N) SUCH THAT IA(IND(J)) IS IN ASCENDING ORDER.
C     THE INPUT ARRAY IA(N) AND N ARE NOT CHANGED.
C     THE OUTPUT ARRAY IND(N) CONTAINS THE INFORMATION FROM THE SORT.
C     *
C     ADAPTED FROM
C     W.H.PRESS, B.P.FLANNERY, S.A.TEUKOLSKY, AND V.T.VETTERLING,
C     NUMERICAL RECIPES, CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE, 1986,
C     CHAPTER 8.3, PAGE 233.
C     *
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
      INTEGER N,L,J,I,IA(N),IND(N),IR,IXT,IQ
      DO 10 J=1,N
      IND(J) = J
   10 CONTINUE
      L      = N/2+1
      IR     = N
   20 CONTINUE
      IF(L.GT.1) THEN
         L   = L-1
         IXT = IND(L)
         IQ  = IA(IXT)
      ELSE
         IXT = IND(IR)
         IQ  = IA(IXT)
         IND(IR) = IND(1)
         IR  = IR-1
         IF(IR.EQ.1) THEN
            IND(1) = IXT
            RETURN
         ENDIF
      ENDIF
      I      = L
      J      = L+L
   30 IF(J.LE.IR) THEN
         IF(J.LT.IR) THEN
            IF(IA(IND(J)).LT.IA(IND(J+1))) J=J+1
         ENDIF
         IF(IQ.LT.IA(IND(J))) THEN
            IND(I) = IND(J)
            I = J
            J = J+J
         ELSE
            J = IR+1
         ENDIF
      GO TO 30
      ENDIF
      IND(I) = IXT
      GO TO 20
      END
