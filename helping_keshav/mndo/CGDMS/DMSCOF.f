C     ******************************************************************
      SUBROUTINE DMSCOF (F,H,G,P,Q1,Q2,LM2,N,NPRINT,B,C,D)
C     *
C     COEFFICIENTS FROM ANALYTIC DMS LINE SEARCH.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     F(LM2,*)  FOCK MATRIX (I).
C     G(LM2,*)  GRADIENT MATRIX (I).
C     H(LM2,*)  SEARCH DIRECTION MATRIX (I).
C     P(LM2,*)  DENSITY MATRIX (I).
C     Q1(LM2,*) SCRATCH PRODUCT MATRIX (S), PRODUCT HF (O).
C     Q2(LM2,*) SCRATCH PRODUCT MATRIX (S), PRODUCT HH (O).
C     LM2       LEADING DIMENSION OF ALL MATRICES (I).
C     N         NUMBER OF ORBITALS (I).
C     NPRINT    PRINTING FLAG (I).
C     B         FIRST  COEFFICIENT (O).
C     C         SECOND COEFFICIENT (O).
C     D         THIRD  COEFFICIENT (O).
C     *
C     TRACE2(A,B,LM2,N) COMPUTES THE TRACE OF A*B.
C     TRACE3(A,B,LM2,N) COMPUTES THE TRACE OF A(transpose)*B.
C     BOTH ARE EQUIVALENT IF MATRIX A IS SYMMETRIC.
C     THE MORE EFFICIENT FUNCTION TRACE3 MAY BE USED IN THIS CASE.
C     THE CODE WILL THEN CONTAIN BOTH CALLS (ONE OF THEM COMMENTED OUT).
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (ONE=1.0D0)
      PARAMETER (TWO=2.0D0)
      PARAMETER (THREE=3.0D0)
      PARAMETER (FOUR=4.0D0)
      COMMON /NBFILE/ NBF(20)
      DIMENSION F(LM2,N),G(LM2,N),H(LM2,N),P(LM2,N)
      DIMENSION Q1(LM2,N),Q2(LM2,N)
C *** COEFFICIENT B = Tr(HG).
C     B   = TRACE2(H,G,LM2,N)
      B   = TRACE3(H,G,LM2,N)
C *** COMPUTE PRODUCTS Q1=HF AND Q2=HP BY MATRIX MULTIPLICATION.
      CALL DGEMM ('N','N',LM2,N,N,ONE,H,LM2,F,LM2,ZERO,Q1,LM2)
      CALL DGEMM ('N','N',LM2,N,N,ONE,H,LM2,P,LM2,ZERO,Q2,LM2)
C *** COEFFICIENT C = 3 Tr(HHF) - 2 Tr(PHHF) - 2 Tr(HPHF) - 2 Tr(HHPF)
C     COEFFICIENT C = 3 Tr(HHF) - 4 Tr(PHHF) - 2 Tr(HPHF)
C     VALID DUE TO THE EQUALITY:    Tr(PHHF) = Tr(HHPF)
C     COEFFICIENT C1 = 3 Tr(HHF).
C     C1 = THREE*TRACE2(H,Q1,LM2,N)
      C1 = THREE*TRACE3(H,Q1,LM2,N)
C     COEFFICIENT C2 = -4 Tr(PHHF) = -4 Tr[(HP)transpose(HF)]
C     PH = (HP)transpose IS VALID SINCE H AND P ARE SYMMETRIC.
      C2 = -FOUR*TRACE3(Q2,Q1,LM2,N)
C     COEFFICIENT C3 = -2 Tr(HPHF).
      C3 = -TWO*TRACE2(Q2,Q1,LM2,N)
      C   = C1+C2+C3
C *** COMPUTE PRODUCT Q2=HH BY MATRIX MULTIPLICATION.
      CALL DGEMM ('N','N',LM2,N,N,ONE,H,LM2,H,LM2,ZERO,Q2,LM2)
C *** COEFFICIENT D = -2 Tr(HHHF).
C     D   = -TWO*TRACE2(Q2,Q1,LM2,N)
      D   = -TWO*TRACE3(Q2,Q1,LM2,N)
C *** DEBUG PRINT - RESULTS.
      IF(NPRINT.GE.7) THEN
         NB6 = NBF(6)
         WRITE (NB6,500) B,C,D
         WRITE (NB6,510) C1,C2,C3
         IF(NPRINT.GT.9) THEN
            WRITE(NB6,520) C1/THREE,-C2/FOUR,-C3/TWO
            WRITE (NB6,530)
            CALL MATPRT (G,N,N,LM2,N)
            WRITE (NB6,540)
            CALL MATPRT (H,N,N,LM2,N)
            WRITE (NB6,550)
            CALL MATPRT (P,N,N,LM2,N)
            WRITE (NB6,560)
            CALL MATPRT (F,N,N,LM2,N)
         ENDIF
      ENDIF
      RETURN
  500 FORMAT(/  1X,'DMSCOF: COEFFICIENT B  =',G20.10,
     1       /  1X,'DMSCOF: COEFFICIENT C  =',G20.10,
     2       /  1X,'DMSCOF: COEFFICIENT D  =',G20.10)
  510 FORMAT(   1X,'DMSCOF: COEFFICIENT C1 =',G20.10,
     1       /  1X,'DMSCOF: COEFFICIENT C2 =',G20.10,
     2       /  1X,'DMSCOF: COEFFICIENT C3 =',G20.10)
  520 FORMAT(   1X,'DMSCOF: Trace(HHF)     =',G20.10,
     1       /  1X,'DMSCOF: Trace(PHHF)    =',G20.10,
     2       /  1X,'DMSCOF: Trace(HPHF)    =',G20.10,
     3       /  1X,'DMSCOF: Trace(HHPF)    = Trace(PHHF)')
  530 FORMAT(///1X,'DMSCOF: GRADIENT MATRIX G.'/)
  540 FORMAT(///1X,'DMSCOF: SEARCH DIRECTION MATRIX H.'/)
  550 FORMAT(///1X,'DMSCOF: DENSITY MATRIX P.'/)
  560 FORMAT(///1X,'DMSCOF: FOCK MATRIX F.'/)
C 570 FORMAT(///1X,'DMSCOF: PRODUCT HF.'/)
C 580 FORMAT(///1X,'DMSCOF: PRODUCT HH.'/)
C 590 FORMAT(///1X,'DMSCOF: PRODUCT HP.'/)
      END