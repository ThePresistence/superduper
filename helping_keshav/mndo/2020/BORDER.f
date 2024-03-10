      SUBROUTINE BORDER (C,D,P,PN,PL,LM2,LM3,LM4,N,NOCC,IODD,JODD,
     1                   NITER,KEXT,INOUT,NB1,ISPIN)
C     *
C     CALCULATION OF DENSITY MATRIX.
C     OPTIONALLY COMBINED WITH EXTRAPOLATION OR DAMPING.
C     GENERAL VERSION USING PACKED LINEAR ARRAYS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     C(LM2,*)  EIGENVECTORS (I).
C     D(LM4)    DIFFERENCE DENSITY MATRIX (I,O).
C     P(LM4)    DENSITY MATRIX (I,O).
C     PN(LM4)   NEW DENSITY MATRIX (S).
C     PL        MAXIMUM CHANGE IN DIAGONAL ELEMENTS (O).
C     N         NUMBER OF BASIS FUNCTIONS (I).
C     NOCC      NUMBER OF OCCUPIED ORBITALS (I).
C     IODD      INDEX OF FIRST  SINGLY OCCUPIED RHF-MO (I).
C     JODD      INDEX OF SECOND SINGLY OCCUPIED RHF-MO (I).
C     NITER     NUMBER OF SCF ITERATION (I).
C     KEXT      TYPE OF SCF ITERATION (I,O).
C               =-1 CONVERGED DENSITY, NO MODIFICATION ALLOWED (I).
C               = 0 STANDARD CASE, NO EXTRAPOLATION OR DAMPING (O).
C               = 1 EXTRAPOLATION FOR DENSITY DONE (O).
C               = 2 DAMPING FOR DENSITY DONE (O).
C     INOUT     FLAG FOR USE OF FILES (I).
C     NB1       FILE NUMBER (I). USED FOR INOUT=2.
C     ISPIN     TYPE OF DENSITY MATRIX (I).
C               = 1 RHF OR UHF-ALPHA DENSITY MATRIX (I).
C               = 2 UHF-BETA DENSITY MATRIX (I).
C     *
C     COMMENTS ON OTHER AVAILABLE INPUT OPTIONS.
C     NSTEP .GT.0  -  ATTEMPT EXTRAPOLATION, IF CURRENTLY POSSIBLE.
C     NSTEP .LT.0  -  ATTEMPT DAMPING, IF CURRENTLY POSSIBLE.
C     NSTART.LT.0  -  NO EXTRAPOLATION OR DAMPING, PL STILL EVALUATED.
C     NITER .LT.0  -  NO EXTRAPOLATION OR DAMPING, PL NOT EVALUATED.
C     IMOCC .LE.1  -  USE DEFAULT OCCUPATION NUMBERS (IMPLICITLY).
C     IMOCC .EQ.2  -  USE DEFAULT OCCUPATION NUMBERS (EXPLICITLY).
C     IMOCC .GT.2  -  USE SPECIAL OCCUPATION NUMBERS OCC(LMX).
C     IMOCC .GE.5  -  USE FERMI-TYPE OCCUPATION NUMBERS (EXPLICITLY).
C     *
      USE LIMIT, ONLY: LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TINY=1.0D-10)
      PARAMETER (YCRIT=0.06D0)
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./OCCFL / IMOCC,NOCCA(2),MSUB,MOSUMA,MOSUMB,MOCCA(8),MOCCB(8)
     ./OCCNM / OCC(LMX,2)
      DIMENSION C(LM2,LM3),D(LM4),P(LM4),PN(LM4)
      SAVE YL
C     ALLOCATABLE SCRATCH ARRAY FOR DIIS TREATMENT.
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: TEMP
C *** INPUT OPTIONS.
      NPRINT = IN2(72)
      NSTART = IN2(78)
      NSTEP  = IN2(79)
C     *
C *** CALCULATE DENSITY MATRIX.
C     *
C     MATRIX MULTIPLICATION PN = C*C(TRANSPOSE).
      CALL ZZERO(N,N,PN,LM2)
      IF(IMOCC.LE.4) THEN
         CALL DGEMM('N','T',N,N,NOCC,ONE,C,LM2,C,LM2,ZERO,PN,LM2)
         CALL LINEAR (PN,PN,N,LM2,LM4)
      ENDIF
C     CORRECTIONS FOR HALF-ELECTRON CASE.
      IF(IMOCC.LE.1) THEN
         IF(IODD.GT.0) THEN
            IJ = 0
            DO 20 I=1,N
            DO 10 J=1,I
            IJ = IJ+1
            PN(IJ) = PN(IJ)-PT5*C(I,IODD)*C(J,IODD)
   10       CONTINUE
   20       CONTINUE
            IF(JODD.GT.0) THEN
               IJ = 0
               DO 40 I=1,N
               DO 30 J=1,I
               IJ = IJ+1
               PN(IJ) = PN(IJ)-PT5*C(I,JODD)*C(J,JODD)
   30          CONTINUE
   40          CONTINUE
            ENDIF
         ENDIF
      ELSE IF(IMOCC.LE.4) THEN
C     CORRECTIONS FOR NON-STANDARD OCCUPATION NUMBERS.
         NMO    = NOCCA(ISPIN)
         DO 70 K=1,NMO
         OCCDIF = OCC(K,ISPIN)-ONE
         IF(ABS(OCCDIF).GT.TINY) THEN
            IJ = 0
            DO 60 I=1,N
            DO 50 J=1,I
            IJ = IJ+1
            PN(IJ) = PN(IJ) + C(I,K)*C(J,K)*OCCDIF
   50       CONTINUE
   60       CONTINUE
         ENDIF
   70    CONTINUE
      ENDIF
C     COMPUTE DENSITY MATRIX (P) USING EXPLICIT OCCUPATION NUMBERS.
      IF(IMOCC.GE.5) THEN
         NMO    = NOCCA(ISPIN)
         ALLOCATE (TEMP(LM2,LM2))
         TEMP(1:N,1:NMO) = C(1:N,1:NMO)
         DO 75 I=1,NMO
         CALL DSCAL (N,OCC(I,ISPIN),TEMP(1:N,I),1)
   75    CONTINUE
C        DENSITY MATRIX IS THE PRODUCT C*OCC*C(T) = TEMP*C(T)
         CALL DGEMM('N','T',N,N,NMO,ONE,TEMP,LM2,C,LM2,ZERO,PN,LM2)
         CALL LINEAR (PN,PN,N,LM2,LM4)
         DEALLOCATE (TEMP)
      ENDIF
C     RESTORE OLD DENSITY MATRIX (P) FROM FILE IF NECESSARY.
      IF(INOUT.GE.2 .AND. NSTART.LT.0) THEN
         IF(ISPIN.EQ.1) REWIND NB1
         READ(NB1) P
      ENDIF
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
C     RESTORE OLD DIFFERENCE DENSITY MATRIX (D) FROM FILE IF NECESSARY.
      IF(NITER.GT.1 .AND. NSTEP.GT.0) THEN
         IF(INOUT.GE.2) THEN
            IF(ISPIN.EQ.1) REWIND NB1
            READ(NB1) D
         ENDIF
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
      IF(NPRINT.GE.6) THEN
         NB6 = NBF(6)
         WRITE(NB6,530)
         WRITE(NB6,540) NITER,KEXT,YL,YLAMB,DEN1,DEN2,FAC
         IF(NPRINT.GE.9 .AND. KEXT.GT.0) THEN
            WRITE(NB6,550) NITER,FAC
            CALL VECPRT (P,LM4,N)
         ENDIF
      ENDIF
      RETURN
  530 FORMAT(// 1X,'BORDER: NITER,KEXT,YL,YLAMB,DEN1,DEN2,FAC')
  540 FORMAT(   1X,'BORDER: ',I5,I3,5G15.5)
  550 FORMAT(// 1X,'BORDER: MODIFIED DENSITY MATRIX, NITER =',I3,
     1          5X,'FAC =',F10.5/)
      END