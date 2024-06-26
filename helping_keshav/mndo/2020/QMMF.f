      SUBROUTINE QMMF (H,PA,PB,Q,LM4,LM6)
C     *
C     CALCULATE THE ELECTRIC FIELD (ESF) GENERATED BY THE QM ATOMS
C     AT THE EXTERNAL POINTS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     H(LM6)    ONE-CENTER ELEMENTS OF CORE HAMILTONIAN MATRIX (S).
C     PA(LM4)   DENSITY MATRIX, ALPHA ELECTRONS (I).
C     PB(LM4)   DENSITY MATRIX, BETA  ELECTRONS (I).
C     Q(LM6)    PRECOMBINED ONE-CENTER DENSITY MATRIX ELEMENTS (S).
C     *
C     UNITS.
C     COORDM()  CARTESIAN COORDINATES OF EXTERNAL POINTS, IN ANGSTROM.
C     EP(L)     ELECTROSTATIC POTENTIAL, IN EV/E = V.
C     ESF(K,M)  ELECTRIC FIELD, IN KCAL/(MOL*ANGSTROM*E).
C     *
C     THE ELECTRIC FIELD (ESF) IS THE NEGATIVE GRADIENT OF THE
C     ELECTROSTATIC POTENTIAL (EP). ESF IS COMPUTED FROM EP BY
C     FINITE DIFFERENCE USING DISPLACEMENTS +DQ AND -DQ.
C     SCHEMATICALLY FOR THE K-TH CARTESIAN COMPONENT:
C     ESF(K) = - (EP(+DQ)-EP(-DQ)) / (TWO*DQ) = (EP(2)-EP(1))*RDQ
C     FOR MORE DETAILS SEE COMMENTS IN SUBROUTINE QMINT.
C     REF: D. BAKOWIES AND W. THIEL, J.COMPUT.CHEM. 17, 87 (1996).
C     *
      USE LIMIT, ONLY: LMI, LM1M
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (DQ=1.0D-06)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
     ./INOPT2/ IN2(300)
     ./QMMM1 / COORDM(3,LM1M),CHARGM(LM1M)
     ./QMMM4 / ELPOT(LM1M),ESF(3,LM1M)
C    ./NBFILE/ NBF(20)
      DIMENSION H(LM6),PA(LM4),PB(LM4),Q(LM6)
      DIMENSION DEL(2),EP(2),P(3)
C *** FILE NUMBERS.
C     NB6    = NBF(6)
C *** INPUT OPTIONS.
      NUMATM = IN2(120)
      MMPOT  = IN2(122)
C *** INITIALIZATION.
      PTCHG  = ONE
      RDQ    = PT5/DQ
      DEL(1) = DQ
      DEL(2) =-DQ
C *** COMPUTE ONE-CENTER DENSITY MATRIX ELEMENTS (Q).
      DO 10 I=1,LM6
      Q(I)  = (PA(IP(I))+PB(IP(I)))
      IF(IP1(I).NE.IP2(I)) Q(I)=Q(I)*TWO
   10 CONTINUE
C *** LOOP OVER ALL EXTERNAL POINTS (M).
      DO 60 M=1,NUMATM
      P(1)   = COORDM(1,M)
      P(2)   = COORDM(2,M)
      P(3)   = COORDM(3,M)
C *** ELECTRIC FIELD.
C     LOOP OVER COMPONENTS (K) AND DISPLACEMENTS (L).
      DO 50 K=1,3
      DO 40 L=1,2
      P(K)   = COORDM(K,M) + DEL(L)
      DO 20 I=1,LM6
      H(I)   = ZERO
   20 CONTINUE
C     INTEGRALS FOR ELECTROSTATIC POTENTIAL, SEE EQ.(2) IN REF.
      IF(MMPOT.LE.2) THEN
         CALL QMINTC (H,LM6,ENUCQM,P,PTCHG,1,0)
      ELSE
         CALL QMINT  (H,LM6,ENUCQM,P,PTCHG,1,0)
      ENDIF
C     ONE-ELECTRON CONTRIBUTION TO ELECTROSTATIC POTENTIAL.
      ELECQM = DDOT(LM6,H,1,Q,1)
C     ELECTROSTATIC POTENTIAL.
      EP(L)  = ELECQM + ENUCQM
   40 CONTINUE
C     FINITE DIFFERENCE FORMULA FOR ELECTRIC FIELD.
      ESF(K,M) = (EP(2)-EP(1))*RDQ*EVCAL
      P(K)     = COORDM(K,M)
C     WRITE(NB6,500) M,K,EP(2),EP(1),ESF(K,M)
   50 CONTINUE
   60 CONTINUE
      RETURN
C 500 FORMAT(1X,'QMMF:M,K,EP(2),EP(1),ESF(K,M)',2I5,3F15.10)
      END
