C     ******************************************************************
      SUBROUTINE ROTBP (IAR,JAC,IAOS,JAOS,T,YY,H,LMR,LMF)
C     *
C     THIS ROUTINE TRANSFORMS TWO-CENTER ONE-ELECTRON INTEGRALS FROM
C     LOCAL TO MOLECULAR COORDINATES, AND INCLUDES THEM IN H(LMR,LMF).
C     *
C     NOTATION.
C     IAR       ROW    INDEX OF S ORBITAL AT ATOM I IN H(LMR,LMF).
C     JAC       COLUMN INDEX OF S ORBITAL AT ATOM J IN H(LMR,LMF).
C     IAOS      NUMBER OF ATOMIC ORBITALS AT ATOM I.
C     JAOS      NUMBER OF ATOMIC ORBITALS AT ATOM J.
C     T         LOCAL TWO-CENTER ONE-ELECTRON INTEGRALS.
C     YY        PRECALCULATED COMBINATION OF ROTATION MATRIX ELEMENTS.
C     H         CORE HAMILTONIAN MATRIX.
C     LMR       ROW    DIMENSION OF CORE HAMILTONIAN MATRIX.
C     LMF       COLUMN DIMENSION OF CORE HAMILTONIAN MATRIX.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(LMR,LMF)
      DIMENSION T(14),YY(15,45)
      DIMENSION HDD(15)
      DIMENSION INDX(5)
      DATA INDX /0,1,3,6,10/
C *** SECTION FOR AN SP-BASIS.
C     S(I)-S(J)
      H(IAR,JAC)  = T(1)
      IF(IAOS.EQ.1 .AND. JAOS.EQ.1) RETURN
C     S(I)-P(J)
      IF(JAOS.GT.1) THEN
         H(IAR,JAC+1) = T(2)*YY(1,2)
         H(IAR,JAC+2) = T(2)*YY(2,2)
         H(IAR,JAC+3) = T(2)*YY(3,2)
      ENDIF
C     P(I)-S(J)
      IF(IAOS.GT.1) THEN
         H(IAR+1,JAC) = T(3)*YY(1,2)
         H(IAR+2,JAC) = T(3)*YY(2,2)
         H(IAR+3,JAC) = T(3)*YY(3,2)
C     P(I)-P(J).
         IF(JAOS.GT.1) THEN
            T45     = T(4)-T(5)
            H(IAR+1,JAC+1) = YY(1,3)*T45+T(5)
            H(IAR+1,JAC+2) = YY(2,3)*T45
            H(IAR+1,JAC+3) = YY(4,3)*T45
            H(IAR+2,JAC+1) = H(IAR+1,JAC+2)
            H(IAR+2,JAC+2) = YY(3,3)*T45+T(5)
            H(IAR+2,JAC+3) = YY(5,3)*T45
            H(IAR+3,JAC+1) = H(IAR+1,JAC+3)
            H(IAR+3,JAC+2) = H(IAR+2,JAC+3)
            H(IAR+3,JAC+3) = YY(6,3)*T45+T(5)
         ENDIF
      ENDIF
C *** SECTION INVOLVING D ORBITALS.
C     D(I)-S(J)
      IF(IAOS.EQ.9) THEN
         DO 10 I=1,5
         M      = IAR+3+I
         H(M,JAC) = T(6)*YY(I,11)
   10    CONTINUE
C     D(I)-P(J)
         IF(JAOS.GT.1) THEN
            IJ     = 0
            DO 30 I=1,5
            M      = IAR+3+I
            DO 20 J=1,3
            IJ     = IJ+1
            H(M,JAC+J) = T(8) * YY(IJ,12)
     1                  +T(10)*(YY(IJ,18)+YY(IJ,25))
   20       CONTINUE
   30       CONTINUE
         ENDIF
      ENDIF
C     S(I)-D(J)
      IF(JAOS.EQ.9) THEN
         DO 40 I=1,5
         H(IAR,JAC+3+I) = T(7)*YY(I,11)
   40    CONTINUE
         IF(IAOS.GT.1) THEN
C     P(I)-D(J)
            DO 60 I=1,3
            M      = IAR+I
            DO 50 J=1,5
            IJ     = 3*(J-1)+I
            H(M,JAC+3+J) = T(9) * YY(IJ,12)
     1                    +T(11)*(YY(IJ,18)+YY(IJ,25))
   50       CONTINUE
   60       CONTINUE
C     D(I)-D(J)
            IF(IAOS.EQ.9) THEN
               DO 70 I=1,15
               HDD(I) = T(12)* YY(I,15)
     1                 +T(13)*(YY(I,21)+YY(I,28))
     2                 +T(14)*(YY(I,36)+YY(I,45))
   70          CONTINUE
               DO 90 I=1,5
               M      = IAR+3+I
               DO 80 J=1,5
               IF(I.GE.J) THEN
                  IJ  = INDX(I)+J
               ELSE
                  IJ  = INDX(J)+I
               ENDIF
               H(M,JAC+3+J) = HDD(IJ)
   80          CONTINUE
   90          CONTINUE
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
