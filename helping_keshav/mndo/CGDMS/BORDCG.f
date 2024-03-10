      SUBROUTINE BORDCG (D,P,PN,PL,LM4,N,NITER,KEXT)
C     *
C     UPDATE OF DENSITY MATRIX.
C     OPTIONALLY COMBINED WITH EXTRAPOLATION OR DAMPING.
C     CGDMS VERSION USING PACKED LINEAR ARRAYS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     D(LM4)    DIFFERENCE DENSITY MATRIX (I,O).
C     P(LM4)    DENSITY MATRIX (I,O).
C     PN(LM4)   NEW DENSITY MATRIX (I).
C     PL        MAXIMUM CHANGE IN DIAGONAL ELEMENTS (O).
C     N         NUMBER OF BASIS FUNCTIONS (I).
C     NITER     NUMBER OF SCF ITERATION (I).
C     KEXT      TYPE OF SCF ITERATION (I,O).
C               =-1 CONVERGED DENSITY, NO MODIFICATION ALLOWED (I).
C               = 0 STANDARD CASE, NO EXTRAPOLATION OR DAMPING (O).
C               = 1 EXTRAPOLATION FOR DENSITY DONE (O).
C               = 2 DAMPING FOR DENSITY DONE (O).
C     *
C     COMMENTS ON OTHER AVAILABLE INPUT OPTIONS.
C     NSTEP .GT.0  -  ATTEMPT EXTRAPOLATION, IF CURRENTLY POSSIBLE.
C     NSTEP .LT.0  -  ATTEMPT DAMPING, IF CURRENTLY POSSIBLE.
C     NSTART.LT.0  -  NO EXTRAPOLATION OR DAMPING, PL STILL EVALUATED.
C     NITER .LT.0  -  NO EXTRAPOLATION OR DAMPING, PL NOT EVALUATED.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (YCRIT=0.06D0)
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
      DIMENSION D(LM4),P(LM4),PN(LM4)
      SAVE YL
C *** INPUT OPTIONS.
      NPRINT = IN2(72)
      NSTART = IN2(78)
      NSTEP  = IN2(79)
C     *
C *** CALCULATE MAXIMUM CHANGE (PL) IN DIAGONAL MATRIX ELEMENTS.
      IF(NITER.GE.0) THEN
         PL  = ZERO
         DO 80 I=1,N
         II  = I*(I+1)/2
         PL  = MAX(ABS(PN(II)-P(II)),PL)
   80    CONTINUE
      ENDIF
C     *
C *** SIMPLEST CASE: NO EXTRAPOLATION OR DAMPING.
C     *
C     COPY DENSITY MATRIX AND RETURN.
      IF(NITER.LT.0 .OR. NSTART.LT.0 .OR. KEXT.LT.0) THEN
         CALL DCOPY (LM4,PN,1,P,1)
         RETURN
      ENDIF
C     *
C *** SECTION FOR EXTRAPOLATION OR DAMPING.
C     *
C *** INITIALIZATION.
      IF(NITER.EQ.1 .AND. NSTEP.GT.0) THEN
         YL = ZERO
         DO 110 I=1,LM4
         D(I)   = ZERO
  110    CONTINUE
         YLAMB  = ZERO
         DEN1   = ZERO
         DEN2   = ZERO
      ENDIF
C     COMPUTE DOT PRODUCTS INVOLVING OLD DIFFERENCE DENSITY MATRIX.
      IF(NITER.GT.1 .AND. NSTEP.GT.0) THEN
         DEN1   = DDOT(LM4,D,1,D,1)
         YLAMB  = DDOT(LM4,D,1,PN,1)-DDOT(LM4,D,1,P,1)
      ENDIF
C *** COMPUTE NEW DIFFERENCE DENSITY MATRIX D(LM2,LM3).
      DO 120 I=1,LM4
      D(I) = PN(I)-P(I)
  120 CONTINUE
C *** EXTRAPOLATION.
      KEXT   = 0
      IF(NSTEP.GT.0) THEN
         FAC    = ONE
         DEN2   = DDOT(LM4,D,1,D,1)
         IF(NITER.GE.NSTART) THEN
            II  = NITER-NSTART
            II  = II-(II/NSTEP)*NSTEP
            IF(II.EQ.0) KEXT=1
         ENDIF
         IF(NITER.LT.2 .OR. DEN1.EQ.ZERO .OR. DEN2.EQ.ZERO) THEN
            KEXT = 0
         ELSE
            YLAMB = YLAMB/DEN1
            IF(ABS(YLAMB).GE.ONE) YLAMB=YLAMB*DEN1/DEN2
            IF(ABS(YLAMB-YL).GT.YCRIT) KEXT=0
            IF(YLAMB.EQ.ONE) KEXT=0
            YL = YLAMB
         ENDIF
         IF(KEXT.EQ.1) FAC = ONE/(ONE-YLAMB)-ONE
C *** DAMPING.
      ELSE IF(NSTEP.LT.0) THEN
         IF(NITER.GE.NSTART .AND. NITER.GT.1) THEN
            KEXT = 2
            FAC  = DBLE(MIN(-NSTEP,9))/10.0D0
         ENDIF
      ENDIF
C *** UPDATE OF DENSITY MATRIX.
      IF(KEXT.GT.0) THEN
         DO 130 I=1,LM4
         P(I)   = PN(I)+FAC*D(I)
  130    CONTINUE
      ELSE
         CALL DCOPY (LM4,PN,1,P,1)
      ENDIF
C *** DEBUG PRINT.
      IF(NPRINT.GE.5) THEN
         NB6 = NBF(6)
         WRITE(NB6,530)
         WRITE(NB6,540) NITER,KEXT,YL,YLAMB,DEN1,DEN2,FAC
         IF(NPRINT.GE.9 .AND. KEXT.GT.0) THEN
            WRITE(NB6,550) NITER,FAC
            CALL VECPRT (P,LM4,N)
         ENDIF
      ENDIF
      RETURN
  530 FORMAT(// 1X,'BORDCG: NITER,KEXT,YL,YLAMB,DEN1,DEN2,FAC')
  540 FORMAT(   1X,'BORDCG: ',I5,I3,5G15.5)
  550 FORMAT(// 1X,'BORDCG: MODIFIED DENSITY MATRIX, NITER =',I3,
     1          5X,'FAC =',F10.5/)
      END