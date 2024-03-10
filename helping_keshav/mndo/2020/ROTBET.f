      SUBROUTINE ROTBET (IA,JA,IORBS,JORBS,T,YY,H,LM4)
C     *
C     THIS ROUTINE TRANSFORMS TWO-CENTER ONE-ELECTRON INTEGRALS FROM
C     LOCAL TO MOLECULAR COORDINATES, AND INCLUDES THEM IN H(LM4).
C     USEFUL FOR RESONANCE INTEGRALS AND FOR OVERLAP INTEGRALS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     IA        INDEX OF FIRST ORBITAL AT ATOM I (I).
C     JA        INDEX OF FIRST ORBITAL AT ATOM J (I).
C     IORBS     NUMBER OF ORBITALS AT ATOM I (I).
C     JORBS     NUMBER OF ORBITALS AT ATOM J (I).
C     T(14)     LOCAL TWO-CENTER ONE-ELECTRON INTEGRALS (I).
C     YY()      PRECOMPUTED COMBINATION OF ROTATION MATRIX ELEMENTS (I).
C     H(LM4)    ONE-ELECTRON MATRIX IN MOLECULAR COORDINATES (O).
C     *
      USE LIMIT, ONLY: LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./INDEX / INDX(LMX)
      DIMENSION H(LM4)
      DIMENSION T(14),YY(15,45)
      DIMENSION HDD(15)
C *** SECTION FOR AN SP-BASIS.
C     S(I)-S(J)
      IS     = INDX(IA)+JA
      H(IS)  = T(1)
      IF(IORBS.EQ.1 .AND. JORBS.EQ.1) RETURN
C     S(I)-P(J)
      IF(JORBS.GE.4) THEN
         H(IS+1) = T(2)*YY(1,2)
         H(IS+2) = T(2)*YY(2,2)
         H(IS+3) = T(2)*YY(3,2)
      ENDIF
C     P(I)-S(J)
      IF(IORBS.GE.4) THEN
         IX      = INDX(IA+1)+JA
         IY      = INDX(IA+2)+JA
         IZ      = INDX(IA+3)+JA
         H(IX)   = T(3)*YY(1,2)
         H(IY)   = T(3)*YY(2,2)
         H(IZ)   = T(3)*YY(3,2)
C     P(I)-P(J).
         IF(JORBS.GE.4) THEN
            T45     = T(4)-T(5)
            H(IX+1) = YY(1,3)*T45+T(5)
            H(IX+2) = YY(2,3)*T45
            H(IX+3) = YY(4,3)*T45
            H(IY+1) = H(IX+2)
            H(IY+2) = YY(3,3)*T45+T(5)
            H(IY+3) = YY(5,3)*T45
            H(IZ+1) = H(IX+3)
            H(IZ+2) = H(IY+3)
            H(IZ+3) = YY(6,3)*T45+T(5)
         ENDIF
      ENDIF
C *** SECTION INVOLVING D ORBITALS.
C     D(I)-S(J)
      IF(IORBS.GE.9) THEN
         DO 10 I=1,5
         M      = INDX(IA+3+I)+JA
         H(M)   = T(6)*YY(I,11)
   10    CONTINUE
C     D(I)-P(J)
         IF(JORBS.GE.4) THEN
            IJ     = 0
            DO 30 I=1,5
            M      = INDX(IA+3+I)+JA
            DO 20 J=1,3
            IJ     = IJ+1
            H(M+J) = T(8) * YY(IJ,12)
     1              +T(10)*(YY(IJ,18)+YY(IJ,25))
   20       CONTINUE
   30       CONTINUE
         ENDIF
      ENDIF
C     S(I)-D(J)
      IF(JORBS.GE.9) THEN
         M      = INDX(IA)+JA+3
         DO 40 I=1,5
         H(M+I) = T(7)*YY(I,11)
   40    CONTINUE
         IF(IORBS.GE.4) THEN
C     P(I)-D(J)
            DO 60 I=1,3
            M      = INDX(IA+I)+JA+3
            DO 50 J=1,5
            IJ     = 3*(J-1)+I
            H(M+J) = T(9) * YY(IJ,12)
     1              +T(11)*(YY(IJ,18)+YY(IJ,25))
   50       CONTINUE
   60       CONTINUE
C     D(I)-D(J)
            IF(IORBS.GE.9) THEN
               DO 70 I=1,15
               HDD(I) = T(12)* YY(I,15)
     1                 +T(13)*(YY(I,21)+YY(I,28))
     2                 +T(14)*(YY(I,36)+YY(I,45))
   70          CONTINUE
               DO 90 I=1,5
               M      = INDX(IA+3+I)+JA+3
               DO 80 J=1,5
               IF(I.GE.J) THEN
                  IJ  = INDX(I)+J
               ELSE
                  IJ  = INDX(J)+I
               ENDIF
               H(M+J) = HDD(IJ)
   80          CONTINUE
   90          CONTINUE
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
