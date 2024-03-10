      SUBROUTINE INDEXX (N,A,IND)
C     *
C     INDEXES A(N) SUCH THAT A(IND(J)) IS IN ASCENDING ORDER.
C     THE INPUT ARRAY A(N) AND N ARE NOT CHANGED.
C     THE OUTPUT ARRAY IND(N) CONTAINS THE INFORMATION FROM THE SORT.
C     *
C     ADAPTED FROM
C     W.H.PRESS, B.P.FLANNERY, S.A.TEUKOLSKY, AND V.T.VETTERLING,
C     NUMERICAL RECIPES, CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE, 1986,
C     CHAPTER 8.3, PAGE 233.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),IND(N)
      DO 10 J=1,N
      IND(J) = J
   10 CONTINUE
      L      = N/2+1
      IR     = N
   20 CONTINUE
      IF(L.GT.1) THEN
         L   = L-1
         IXT = IND(L)
         Q   = A(IXT)
      ELSE
         IXT = IND(IR)
         Q   = A(IXT)
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
            IF(A(IND(J)).LT.A(IND(J+1))) J=J+1
         ENDIF
         IF(Q.LT.A(IND(J))) THEN
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
