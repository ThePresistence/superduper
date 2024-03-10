      SUBROUTINE ROTCOH (IA,JA,IORBS,JORBS,IP,JP,CORE,YY,H,LMH)
C     *
C     THIS ROUTINE TRANSFORMS THE CORE ELECTRON ATTRACTION INTEGRALS
C     FROM LOCAL TO MOLECULAR COORDINATES, AND INCLUDES THEM IN THE
C     CORE HAMILTONIAN.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     IA        INDEX OF FIRST BASIS ORBITAL AT ATOM I (I).
C     JA        INDEX OF FIRST BASIS ORBITAL AT ATOM J (I).
C     IORBS     NUMBER OF BASIS ORBITALS AT ATOM I (I).
C     JORBS     NUMBER OF BASIS ORBITALS AT ATOM J (I).
C     IP        INDEX FOR (S,S) OF ATOM I IN LINEAR ARRAY H (I,S).
C     JP        INDEX FOR (S,S) OF ATOM J IN LINEAR ARRAY H (I,S).
C     CORE()    LOCAL CORE ELECTRON ATTRACTION INTEGRALS (I).
C     YY()      PRECOMBINED ELEMENTS OF ROTATION MATRIX (I).
C     H(LMH)    CORE HAMILTONIAN MATRIX (O).
C     *
C     DEPENDING ON THE INPUT DATA IN THE ARGUMENT LIST, THE INTEGRALS
C     ARE INCLUDED EITHER IN THE FULL CORE HAMILTONIAN H(LM4) OR IN
C     THE ONE-CENTER PART H(LM6) OF THE CORE HAMILTONIAN.
C
C     ARGUMENT       FIRST CASE              SECOND CASE
C     IA             NFIRST(I)               1
C     JA             NFIRST(J)               1
C     IP             INDX(IA)+IA             NW(I)
C     JP             INDX(JA)+JA             NW(J)
C     LMH            LM4                     LM6
C     *
C     FOR IORBS.LE.0 OR JORBS.LE.0, THERE ARE NO BASIS FUNCTIONS AT
C     ATOMS I OR J, RESPECTIVELY, AND HENCE NO CONTRIBUTIONS TO THE
C     CORE HAMILTONIAN. THIS CASE MAY ARISE IN CALCULATIONS WITH
C     EXTERNAL POINT CHARGES.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(LMH)
      DIMENSION CORE(10,2),YY(15,45)
      DIMENSION HPP(6),HDP(15),HDD(15)
C$DIR SCALAR
*VDIR NOVECTOR
*VOCL LOOP,SCALAR
      DO 60 KK=1,2
      IF(KK.EQ.1) THEN
         IS  = IP
         K   = IA-1
         L   = IORBS
      ELSE
         IS  = JP
         K   = JA-1
         L   = JORBS
      ENDIF
      IF(L.LE.0) GO TO 60
C     S-S
      H(IS)  = H(IS)+CORE(1,KK)
      IF(L.EQ.1) GO TO 60
C     INTERMEDIATE RESULTS FOR P-P
      DO 10 I=1,6
      HPP(I) = CORE(3,KK)*YY(I,3)+CORE(4,KK)*(YY(I,6)+YY(I,10))
   10 CONTINUE
C     P-S
      IX     = IS+1+K
      IY     = IX+2+K
      IZ     = IY+3+K
      H(IX)  = H(IX)+CORE(2,KK)*YY(1,2)
      H(IY)  = H(IY)+CORE(2,KK)*YY(2,2)
      H(IZ)  = H(IZ)+CORE(2,KK)*YY(3,2)
C     P-P
      H(IX+1)= H(IX+1)+HPP(1)
      H(IY+1)= H(IY+1)+HPP(2)
      H(IY+2)= H(IY+2)+HPP(3)
      H(IZ+1)= H(IZ+1)+HPP(4)
      H(IZ+2)= H(IZ+2)+HPP(5)
      H(IZ+3)= H(IZ+3)+HPP(6)
C     INTERMEDIATE RESULTS FOR D-P AND D-D
      IF(L.EQ.4) GO TO 60
      DO 20 I=1,15
      HDP(I) = CORE(6,KK)*YY(I,12)+CORE( 8,KK)*(YY(I,18)+YY(I,25))
      HDD(I) = CORE(7,KK)*YY(I,15)+CORE( 9,KK)*(YY(I,21)+YY(I,28))
     1                            +CORE(10,KK)*(YY(I,36)+YY(I,45))
   20 CONTINUE
C     D-S
      IDP    = 0
      IDD    = 0
      ID     = IZ+3+K
      DO 50 I=5,9
      H(ID+1)= H(ID+1)+CORE(5,KK)*YY(I-4,11)
C     D-P
      DO 30 J=2,4
      IDP    = IDP+1
      H(ID+J)= H(ID+J)+HDP(IDP)
   30 CONTINUE
C     D-D
      DO 40 J=5,I
      IDD    = IDD+1
      H(ID+J)= H(ID+J)+HDD(IDD)
   40 CONTINUE
      ID     = ID+I+K
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
