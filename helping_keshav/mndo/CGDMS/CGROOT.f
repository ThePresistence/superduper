C     ******************************************************************
      SUBROUTINE CGROOT (P,H,F,Q1,Q2,Q3,LM2,N,XMU,B,C,D,XSQMAX,XI,FI,
     1                   DI,TR1,TR3,NITER,ICG,NCG,NPRINT,ICALL)
C     *
C     SELECTION OF PROPER ROOT FOR CONJUGATE GRADIENT UPDATE.
C     CHOOSE ROOT ASSOCIATED WITH LOWER VALUE OF THE FUNCTIONAL.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     P(LM2,*)  DENSITY MATRIX (I), UNCHANGED (O).
C     H(LM2,*)  SEARCH DIRECTION MATRIX (I), UNCHANGED (O).
C     F(LM2,*)  FOCK MATRIX (I), UNCHANGED (O).
C     Q1(LM2,*) SCRATCH ARRAY (S), UPDATED DENSITY MATRIX (O).
C     Q2(LM2,*) SCRATCH ARRAY (S).
C     Q3(LM2,*) SCRATCH ARRAY (S), PURIFIED DENSITY MATRIX (O).
C     LM2       LEADING DIMENSION OF ALL MATRICES (I).
C     N         NUMBER OF ORBITALS (I).
C     XMU       LAGRANGEAN MULTIPLIER (I).
C     B,C,D     COEFFICIENTS IN QUADRATIC EQUATION (I).
C     XSQMAX    CRITERION FOR ADOPTING SOLUTION OF LINEAR EQUATION (I).
C     XI        SELECTED ROOT (O).
C     FI        ASSOCIATED VALUE OF FUNCTIONAL (O).
C     DI        MAXIMUM DEVIATION FROM IDEMPOTENCY FOR Q1 (O).
C     TR1       TRACE OF UPDATED DENSITY MATRIX Q1 (O).
C     TR3       TRACE OF PURIFIED DENSITY MATRIX Q3 (O).
C     NITER     NUMBER OF CURRENT SCF ITERATION (I).
C     ICG       NUMBER OF CURRENT CG CYCLE (I).
C     NCG(*)    TOTAL NUMBER OF OPERATIONS DURING CG SEARCHES (I,O).
C               (2) PURIFICATIONS DONE.
C               (3) MATRIX MULTIPLICATIONS DONE.
C               (4) LINEAR ROOTS ACCEPTED.
C               (5) QUADRATIC ROOTS SELECTED ON PHYSICAL GROUNDS.
C               (6) QUADRATIC ROOTS SELECTED BY COMPARISON OF F.
C     NPRINT    PRINTING FLAG (I).
C     ICALL     ERROR FLAG (I,O), NORMALLY NOT CHANGED.
C               =-1 ON OUTPUT: FATAL ERROR IN CGROOT.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL NOTX1,NOTX2,BOTHOK
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./NBFILE/ NBF(20)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
      DIMENSION P(LM2,N),H(LM2,N),F(LM2,N)
      DIMENSION Q1(LM2,N),Q2(LM2,N),Q3(LM2,N)
      DIMENSION NCG(10)
C *** INITIALIZATION.
      NB6    = NBF(6)
      XEL    = DBLE(NALPHA+NBETA)
      DI     = ZERO
      TR1    = ZERO
      TR3    = ZERO
C     REDEFINE COEFFICIENTS.
      CC     = TWO*C
      DD     = THREE*D
C     LIMITS FOR DIAGONAL DENSITY MATRIX ELEMENTS.
C     THESE ELEMENTS SHOULD BE BETWEEN 0 AND 1 UNDER OUR CONVENTIONS.
C     MUCH LOOSER LIMITS ARE USED TO AVOID ROUNDING ERROR PROBLEMS.
      PPMIN  =-ONE
      PPMAX  = TWO  
C     *
C *** SOLUTION OF LINEAR EQUATION.
C     CHECK CRITERION FOR SOLVING LINEAR EQUATION (MSQ FLAG).
      MSQ    = 0
      IF(CC.NE.ZERO) THEN
         XSQ = (FOUR*B*DD/CC)/CC
         IF(ABS(XSQ).LT.XSQMAX) THEN
            MSQ = 1
         ENDIF
         IF(NPRINT.GE.7) THEN
            BCC = -B/CC
            WRITE(NB6,300) XSQ,BCC
         ENDIF
      ENDIF
C     SELECT UPDATE FROM LINEAR EQUATION.
      IF(DD.EQ.ZERO .OR. MSQ.GT.0) THEN
         XI     = -B/CC
         CALL MATUPD (P,H,Q1,XI,LM2,N)
         TR1    = TRACE1 (Q1,LM2,N)
         CALL PURIFY (Q1,Q2,Q3,LM2,N)
C        TRACE3 MAY BE USED SINCE Q3 IS SYMMETRIC.
C        FI     = TRACE2(Q3,F,LM2,N) + XMU*(TR1-XEL)
         FI     = TRACE3(Q3,F,LM2,N) + XMU*(TR1-XEL)
         NCG(2) = NCG(2)+1
         NCG(3) = NCG(3)+2
         NCG(4) = NCG(4)+1
         IF(NPRINT.GE.5) THEN
            TR3 = TRACE1 (Q3,LM2,N)
            CALL MATDEV (Q1,Q2,LM2,N,1,DI,DIRMS)
            IF(NPRINT.GE.7) THEN
               WRITE(NB6,310)
               WRITE(NB6,320) NITER,ICG,XI,FI
               WRITE(NB6,530) NITER,ICG,DI,DIRMS
               WRITE(NB6,540) NITER,ICG,TR1,TR3
               WRITE(NB6,550) XEL
            ENDIF
         ENDIF
         RETURN
      ENDIF
C     *
C *** SOLUTION OF QUADRATIC EQUATION.
C     ERROR EXIT: THERE IS NO REAL SOLUTION.
      TMP    = CC**2 - FOUR*B*DD
      IF(TMP.LT.ZERO) THEN
         WRITE(NB6,900) TMP
C        STOP 'CGROOT'
         ICALL = -1
         RETURN
      ENDIF
C     SELECT UPDATE FROM ROOTS X1 AND X2.
      X1 = (-CC+SQRT(TMP))/(TWO*DD)
      X2 = (-CC-SQRT(TMP))/(TWO*DD)
C *** CHECK THE TWO ROOTS WITH REGARD TO THE CRITERION THAT THE DIAGONAL
C     ELEMENTS OF THE UPDATED DENSITY MATRIX MUST BE BETWEEN 0 AND 2.
      CALL MATUPM (P,H,X1,LM2,N,UPMIN1,UPMAX1)
      CALL MATUPM (P,H,X2,LM2,N,UPMIN2,UPMAX2)
      NOTX1  = UPMIN1.LT.PPMIN .OR. UPMAX1.GT.PPMAX
      NOTX2  = UPMIN2.LT.PPMIN .OR. UPMAX2.GT.PPMAX
      IF(NOTX1 .AND. NOTX2) THEN
         WRITE(NB6,400)
         WRITE(NB6,410) NITER,ICG,X1,UPMIN1,UPMAX1
         WRITE(NB6,420) NITER,ICG,X2,UPMIN2,UPMAX2
         WRITE(NB6,430)
C        STOP 'CGROOT'
         ICALL = -1
         RETURN
      ELSE IF(NOTX1) THEN
         BOTHOK = .FALSE.
         NCG(5) = NCG(5)+1
         Y1     = X2
      ELSE IF(NOTX2) THEN
         BOTHOK = .FALSE.
         NCG(5) = NCG(5)+1
         Y1     = X1
      ELSE
         BOTHOK = .TRUE.
         NCG(6) = NCG(6)+1
         Y1     = X1
         Y2     = X2
      ENDIF
C *** EVALUATE FUNCTIONAL FOR THE FIRST ROOT (Y1).
      CALL MATUPD (P,H,Q1,Y1,LM2,N)
      TR11   = TRACE1(Q1,LM2,N) 
      CALL PURIFY (Q1,Q2,Q3,LM2,N)
C     TRACE3 MAY BE USED SINCE Q3 IS SYMMETRIC.
C     F1     = TRACE2(Q3,F,LM2,N) + XMU*(TR11-XEL)
      F1     = TRACE3(Q3,F,LM2,N) + XMU*(TR11-XEL)
      FI     = F1
      XI     = Y1
      TR1    = TR11
      NCG(2) = NCG(2)+1
      NCG(3) = NCG(3)+2
      IF(NPRINT.GE.5) THEN
         TR31 = TRACE1 (Q3,LM2,N)
         CALL MATDEV (Q1,Q2,LM2,N,1,DIDEM1,DIRMS1)
         TR3  = TR31
         DI   = DIDEM1
      ENDIF
C *** EVALUATE FUNCTIONAL FOR THE SECOND ROOT (Y2).
C     MAKE SURE THAT ARRAY Q3 IS NOT OVERWRITTEN.
      IF(BOTHOK) THEN 
         CALL MATUPD (P,H,Q1,Y2,LM2,N)
         TR12   = TRACE1(Q1,LM2,N)
         CALL PURIFY (Q1,Q2,P,LM2,N)
C        TRACE3 MAY BE USED SINCE P IS SYMMETRIC.
C        F2     = TRACE2(P,F,LM2,N) + XMU*(TR12-XEL)
         F2     = TRACE3(P,F,LM2,N) + XMU*(TR12-XEL)
         NCG(2) = NCG(2)+1
         NCG(3) = NCG(3)+2
         IF(NPRINT.GE.5) THEN
            TR32 = TRACE1 (Q3,LM2,N)
            CALL MATDEV (Q1,Q2,LM2,N,1,DIDEM2,DIRMS2)
         ENDIF
C        ADOPT CG UPDATE LEADING TO LOWER VALUE OF FUNCTIONAL (FI).
C        CASE F1 ABOVE F2. REDEFINE Q3 AND P.
         IF(F1.GT.F2) THEN
            XI  = Y2
            FI  = F2
            DI  = DIDEM2
            TR1 = TR12
            TR3 = TR32
            CALL MATCOP (P,Q3,LM2,N,1)
            CALL MATUPD (Q1,H,P,-Y2,LM2,N)
         ELSE
C        CASE F1 BELOW F2. REDEFINE P AND Q1.
            CALL MATUPD (Q1,H,P,-Y2,LM2,N)
            CALL MATUPD (P,H,Q1,Y1,LM2,N)
         ENDIF
      ENDIF
C *** DEBUG PRINT.
      IF(NPRINT.GE.7) THEN
         WRITE(NB6,400)
         WRITE(NB6,410) NITER,ICG,X1,UPMIN1,UPMAX1
         WRITE(NB6,420) NITER,ICG,X2,UPMIN2,UPMAX2
         IF(BOTHOK) THEN
            WRITE(NB6,500)
            WRITE(NB6,520) NITER,ICG,Y1,F1,Y2,F2
            WRITE(NB6,530) NITER,ICG,DIDEM1,DIRMS1,DIDEM2,DIRMS2
            WRITE(NB6,540) NITER,ICG,TR11,TR31,TR12,TR32
            WRITE(NB6,550) XEL
         ELSE
            WRITE(NB6,510)
            WRITE(NB6,520) NITER,ICG,Y1,F1
            WRITE(NB6,530) NITER,ICG,DIDEM1,DIRMS1
            WRITE(NB6,540) NITER,ICG,TR1,TR3
            WRITE(NB6,550) XEL
         ENDIF
      ENDIF
      RETURN
  300 FORMAT(   1X,'CGROOT: TERM 4*B*D/C/C =',G20.10,
     1       /  1X,'CGROOT: LINEAR ROOT XI =',G20.10)
  310 FORMAT(   1X,'CGROOT: NITER,ICG,XI,FI - LINEAR ROOT')
  320 FORMAT(   1X,'CGROOT: LINEAR',I5,I3,3X,4G20.10)
  400 FORMAT(   1X,'CGROOT: NITER,ICG,STEP,UPMIN,UPMAX')
  410 FORMAT(   1X,'CGROOT: UPD1  ',I5,I3,3X,3G20.10)
  420 FORMAT(   1X,'CGROOT: UPD2  ',I5,I3,3X,3G20.10)
  430 FORMAT(// 1X,'CGROOT: BOTH ROOTS ARE UNACCEPTABLE.',
     1       /  1X,'CGROOT: THE UPDATE WOULD GENERATE A DENSITY MATRIX',
     2          1X,'WITH UNPHYSICAL DIAGONAL ELEMENTS.')
  500 FORMAT(   1X,'CGROOT: NITER,ICG,X1,F1,X2,F2')
  510 FORMAT(   1X,'CGROOT: NITER,ICG,XI,FI - PRESELECTION')
  520 FORMAT(   1X,'CGROOT: QUADR ',I5,I3,3X,4G20.10)
  530 FORMAT(   1X,'CGROOT: IDEMP ',I5,I3,3X,4G20.10)
  540 FORMAT(   1X,'CGROOT: TRACE ',I5,I3,3X,4G20.10)
  550 FORMAT(   1X,'CGROOT: NUMBER OF ELECTRONS',G20.10)
  900 FORMAT(///5X,'ERROR IN CONJUGATE GRADIENT SEARCH.',
     1       /  5X,'STEP FOR DENSITY MATRIX CANNOT BE DEFINED.',
     2       /  5X,'NEGATIVE ARGUMENT FOR SQUARE ROOT:',G20.10/)
      END