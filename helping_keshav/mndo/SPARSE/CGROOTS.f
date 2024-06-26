      SUBROUTINE CGROOTS (PS,JPS,IPS,HS,JHS,IHS,FS,JFS,IFS,Q8,JQ8,
     +                    IQ8,N,XMU,B,C,D,XSQMAX,XI,FI,DI,TR1,TR3,
     +                    NITER,ICG,NCG,NPRINT,ICALL,FILT,
     +                    FILTALL,CUTM)
C     *
C     SELECTION OF PROPER ROOT FOR CONJUGATE GRADIENT UPDATE.
C     CHOOSE ROOT ASSOCIATED WITH LOWER VALUE OF THE FUNCTIONAL.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     PS(*)     DENSITY MATRIX (I), UNCHANGED (O).
C     HS(*)     SEARCH DIRECTION MATRIX (I), UNCHANGED (O).
C     FS(*)     FOCK MATRIX (I), UNCHANGED (O).
C     Q8(*)     SCRATCH ARRAY (S), PURIFIED DENSITY MATRIX (O).
C     IX(N+1)   POINTERS FOR SPARSE MATRICES IN CSR FORMAT (S).
C               X=FS,GS,HS,PS,Q8.
C     JX(*)     COLUMN INDICES FOR SPARSE MATRICES IN CSR FORMAT (S).
C               X=FS,GS,HS,PS,Q8.
C     N         NUMBER OF ORBITALS (I).
C     XMU       LAGRANGEAN MULTIPLIER (I).
C     B,C,D     COEFFICIENTS IN QUADRATIC EQUATION (I).
C     XSQMAX    CRITERION FOR ADOPTING SOLUTION OF LINEAR EQUATION (I).
C     XI        SELECTED ROOT (O).
C     FI        ASSOCIATED VALUE OF FUNCTIONAL (O).
C     DI        MAXIMUM DEVIATION FROM IDEMPOTENCY FOR Q1 (O).
C     TR1       TRACE OF UPDATED DENSITY MATRIX Q1 (O).
C     TR3       TRACE OF PURIFIED DENSITY MATRIX Q8 (O).
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
C               =-1 ON OUTPUT: FATAL ERROR IN CGROOTS.
C     FILT      FLAG FOR REMOVING SMALL ELEMENTS IN PRODUCT A*B (I).
C     FILTALL   FLAG FOR REMOVING SMALL ELEMENTS IN SUM A+B (I).
C     CUTM      CUTOFF FOR SMALL MATRIX ELEMENTS (I).
C     *
      USE module3
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
      INTERFACE
      SUBROUTINE PURIFYN (PS,JPS,IPS,N,DD1,DD2,NPR,FILT,FILTALL,CUTM)
      USE module3
      IMPLICIT NONE
      INTEGER :: N,NPR
      DOUBLE PRECISION :: DD1,DD2,CUTM
      DOUBLE PRECISION, DIMENSION (:), POINTER :: PS
      INTEGER, DIMENSION (:), POINTER :: IPS,JPS
      LOGICAL :: FILT,FILTALL
      END SUBROUTINE PURIFYN
      END INTERFACE
C
      LOGICAL NOTX1,NOTX2,BOTHOK
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./NBFILE/ NBF(20)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
      INTEGER :: N,NCG(10),ie,NITER,NB6,ICALL,N4,NPRINT,NUMB,NORBS,
     +           NMOS,NALPHA,NBETA,ICG,NBF,MSQ,I,J
      DOUBLE PRECISION :: B,CC,DD,XI,BCC,XMU,QQ,FI,TR1,
     +                    TWO,X1,X2,TMP,Y1,Y2,ONE,FOUR,CUTM,F1,TR11,
     +                    TR12,TR31,XEL,TR32,F2,XSQ,D,DIDEM1,DIRMS1,
     +                    PT5,PT25,ZERO,THREE,C,DI,TR3,DIDEM2,DIRMS2,
     +                    DIRMS,UPMIN1,UPMAX1,UPMIN2,UPMAX2,XSQMAX,
     +                    PPMIN,PPMAX
C
      INTEGER, ALLOCATABLE :: IW(:)
      DOUBLE PRECISION, DIMENSION (:), POINTER :: Q8,Q1,PS,FS,HS,WR
      INTEGER, DIMENSION (:), POINTER :: IQ8,JQ8,IQ1,JQ1,JPS,IPS,JFS,
     +                                   IFS,JHS,IHS
      LOGICAL :: FILT,FILTALL,GO1,GO2
      SAVE Q1,JQ1,IQ1
C
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
C     ALLOCATE GENERAL SCRATCH ARRAY.
      ALLOCATE (WR(N),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'CGROOTS','WR',N,0)
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
         ALLOCATE (IW(N),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGROOTS','IW',N,0)
C        PERFORM UPDATE (Q8 = PS + XI*HS).
         CALL aplbdgp (N,JPS,IPS,JHS,IHS,N4,IW)
         ALLOCATE (Q8(N4),JQ8(N4),IQ8(N+1),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGROOTS','Q8',N4,0)
         DO 10 I=1,IHS(N+1)-1
         HS(I)=XI*HS(I)
   10    CONTINUE
         CALL aplbp (N,PS,JPS,IPS,HS,JHS,IHS,Q8,JQ8,IQ8,N4,IW,ie)
            IF(ie.NE.0) CALL XERSPA (ie,'aplb','CGROOTS',1)
         DO 20 I=1,IHS(N+1)-1
         HS(I)=HS(I)/XI
   20    CONTINUE
C        REMOVE SMALL ELEMENTS (Q8).
         IF(FILTALL) THEN
            CALL filterp (N,1,CUTM,Q8,JQ8,IQ8,Q8,JQ8,IQ8,N4,ie)
              IF (ie.ne.0) CALL XERSPA (ie,'filter','CGROOTS',4)
         ENDIF
         DEALLOCATE (IW,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGROOTS','IW',N4,1)
C        COMPUTE TRACE (TR1) OF UPDATED DENSITY MATRIX.
         CALL TRACEAP (N,Q8,JQ8,IQ8,TR1)
C        PERFORM ONE McWEENY PURIFICATION (Q8).
         CALL PURIFYN (Q8,JQ8,IQ8,N,DI,DIRMS,NPRINT,FILT,FILTALL,CUTM)
C        EVALUATE FUNCTIONAL (FI).
         CALL TRACENP (N,Q8,JQ8,IQ8,FS,JFS,IFS,WR,FI)
         DEALLOCATE (WR,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGROOTS','WR',N,1)
         NULLIFY (WR)
         FI     = FI + XMU*(TR1-XEL)
         NCG(2) = NCG(2)+1
         NCG(3) = NCG(3)+2
         NCG(4) = NCG(4)+1
C        DEBUG PRINT.
         IF(NPRINT.GE.5) THEN
            CALL TRACEAP (N,Q8,JQ8,IQ8,TR3)
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
C        STOP 'CGROOTS'
         ICALL = -1
         RETURN
      ENDIF
C     SELECT UPDATE FROM ROOTS X1 AND X2.
      IF (CC.LT.ZERO) THEN
         QQ=-(CC-SQRT(TMP))/2
      ELSE
         QQ=-(CC+SQRT(TMP))/2
      ENDIF
      X1 = QQ/DD
      X2 = B/QQ
C
C *** CHECK THE TWO ROOTS WITH REGARD TO THE CRITERION THAT THE DIAGONAL
C     ELEMENTS OF THE UPDATED DENSITY MATRIX MUST BE REASONABLE.
      CALL SPAUPM (PS,HS,IPS,IHS,JPS,JHS,N,X1,X2,UPMIN1,UPMAX1,
     +             UPMIN2,UPMAX2)
      NOTX1  = UPMIN1.LT.PPMIN .OR. UPMAX1.GT.PPMAX
      NOTX2  = UPMIN2.LT.PPMIN .OR. UPMAX2.GT.PPMAX
      IF(NOTX1 .AND. NOTX2) THEN
         WRITE(NB6,400)
         WRITE(NB6,410) NITER,ICG,X1,UPMIN1,UPMAX1
         WRITE(NB6,420) NITER,ICG,X2,UPMIN2,UPMAX2
         WRITE(NB6,430)
C        STOP 'CGROOTS'
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
C
C *** EVALUATE FUNCTIONAL FOR THE FIRST ROOT (Y1).
      ALLOCATE (IW(N),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'CGROOTS','IW',N,0)
      CALL aplbdgp (N,JPS,IPS,JHS,IHS,N4,IW)
      ALLOCATE (Q8(N4),JQ8(N4),IQ8(N+1),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'CGROOTS','Q8',N4,0)
      DO 30 I=1,IHS(N+1)-1
      HS(I)=Y1*HS(I)
  30  CONTINUE
      CALL aplbp (N,PS,JPS,IPS,HS,JHS,IHS,Q8,JQ8,IQ8,N4,IW,ie)
         IF(ie.NE.0) CALL XERSPA (ie,'aplb','CGROOTS',2)
      DO 40 I=1,IHS(N+1)-1
      HS(I)=HS(I)/Y1
   40 CONTINUE
C     REMOVE SMALL ELEMENTS.
      IF(FILTALL) THEN
         CALL filterp (N,1,CUTM,Q8,JQ8,IQ8,Q8,JQ8,IQ8,N4,ie)
           IF (ie.ne.0) CALL XERSPA (ie,'filter','CGROOTS',4)
      ENDIF
C     COMPUTE TRACE (TR11) OF UPDATED DENSITY MATRIX.
      CALL TRACEAP (N,Q8,JQ8,IQ8,TR11)
C     PERFORM ONE McWEENY PURIFICATION (Q8).
      CALL PURIFYN (Q8,JQ8,IQ8,N,DIDEM1,DIRMS1,NPRINT,FILT,FILTALL,CUTM)
C     EVALUATE FUNCTIONAL (F1).
      CALL TRACENP (N,Q8,JQ8,IQ8,FS,JFS,IFS,WR,F1)
C     CHECK DIAGONAL DENSITY MATRIX ELEMENTS.
      GO2=.FALSE.
      DO 60 I=1,N
      DO 50 J=IQ8(I),IQ8(I+1)-1
      IF (JQ8(J).EQ.I) THEN
         IF ((Q8(J).LT.PPMIN).OR.(Q8(J).GT.PPMAX)) THEN
            GO2=.TRUE.
            GO TO 80
         ENDIF
         GO TO 60
      ENDIF
   50 CONTINUE
   60 CONTINUE
C     SAVE RESULTS FOR FIRST ROOT.
   80 F1     = F1 + XMU*(TR11-XEL)
      FI     = F1
      XI     = Y1
      TR1    = TR11
      NCG(2) = NCG(2)+1
      NCG(3) = NCG(3)+2
      DI   = DIDEM1
      IF(NPRINT.GE.5) THEN
         CALL TRACEAP (N,Q8,JQ8,IQ8,TR31)
         TR3  = TR31
      ENDIF
C
c *** EVALUATE FUNCTIONAL FOR THE SECOND ROOT (Y2).
C     MAKE SURE THAT ARRAY Q8 IS NOT OVERWRITTEN.
      IF(BOTHOK) THEN
         ALLOCATE (Q1(N4),JQ1(N4),IQ1(N+1),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGROOTS','Q1',N4,0)
         DO 90 I=1,IHS(N+1)-1
         HS(I)=Y2*HS(I)
   90    CONTINUE
         CALL aplbp (N,PS,JPS,IPS,HS,JHS,IHS,Q1,JQ1,IQ1,N4,IW,ie)
            IF(ie.NE.0) CALL XERSPA (ie,'aplb','CGROOTS',3)
         DO 100 I=1,IHS(N+1)-1
         HS(I)=HS(I)/Y2
  100    CONTINUE
C        REMOVE SMALL ELEMENTS.
         IF(FILTALL) THEN
            CALL filterp (N,1,CUTM,Q1,JQ1,IQ1,Q1,JQ1,IQ1,N4,ie)
              IF (ie.ne.0) CALL XERSPA (ie,'filter','CGROOTS',4)
         ENDIF
C        COMPUTE TRACE (TR12) OF UPDATED DENSITY MATRIX.
         CALL TRACEAP (N,Q1,JQ1,IQ1,TR12)
C        PERFORM ONE McWEENY PURIFICATION (Q1).
         CALL PURIFYN (Q1,JQ1,IQ1,N,DIDEM2,DIRMS2,NPRINT,FILT,FILTALL,
     +                 CUTM)
C        EVALUATE FUNCTIONAL (F2).
         CALL TRACENP (N,Q1,JQ1,IQ1,FS,JFS,IFS,WR,F2)
C        CHECK DIAGONAL DENSITY MATRIX ELEMENTS.
         GO1=.FALSE.
         DO 120 I=1,N
         DO 110 J=IQ1(I),IQ1(I+1)-1
         IF (JQ1(J).EQ.I) THEN
            IF((Q1(J).LT.ZERO).OR.(Q1(J).GT.TWO)) THEN
               IF(GO2) THEN
                  WRITE(NB6,430)
C                 STOP 'CGROOTS'
                  ICALL = -1
                  RETURN
               ENDIF
               GO TO 130
            ENDIF
            GO TO 120
         ENDIF
  110    CONTINUE
  120    CONTINUE
C        SAVE RESULTS FOR SECOND ROOT.
  130    F2     = F2 + XMU*(TR12-XEL)
         NCG(2) = NCG(2)+1
         NCG(3) = NCG(3)+2
         IF(NPRINT.GE.5) THEN
            CALL TRACEAP (N,Q8,JQ8,IQ8,TR32)
         ENDIF
C        ADOPT CG UPDATE LEADING TO LOWER VALUE OF FUNCTIONAL (FI).
C        CASE F1 ABOVE F2. REDEFINE Q8 AND P.
         IF ((F1.GT.F2).OR.(GO2)) THEN
            XI  = Y2
            FI  = F2
            DI  = DIDEM2
            TR1 = TR12
            TR3 = TR32
C           SAVE THE APPROPRIATE DENSITY MATRIX (Q8).
            N4=IQ1(N+1)-1
            DEALLOCATE (Q8,JQ8,STAT=ie)
               IF(ie.NE.0) CALL XERALL (ie,'CGROOTS','Q8',N4,1)
            NULLIFY (Q8,JQ8)
            ALLOCATE (Q8(N4),JQ8(N4),STAT=ie)
               IF(ie.NE.0) CALL XERALL (ie,'CGROOTS','Q8',N4,0)
            CALL copmatp (N,Q1,JQ1,IQ1,Q8,JQ8,IQ8,1)
         ENDIF
      ENDIF
C
C *** DEALLOCATE SCRATCH ARRAYS.
      DEALLOCATE (IW,WR,STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'CGROOTS','IW',N,1)
      NULLIFY (WR)
      IF (ASSOCIATED(Q1)) THEN
         DEALLOCATE (Q1,JQ1,IQ1,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGROOTS','Q1',N4,1)
         NULLIFY (Q1,JQ1,IQ1)
      ENDIF
C
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
