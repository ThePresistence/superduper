      SUBROUTINE CONINT (ARRAY,LM5,ICALL,SCFCAL)
C     *
C     CONTROL ROUTINE FOR CONICAL INTERSECTION SEARCH.
C     *
C     FIRST PART : MULTIPLE CI RUNS TO GET NEEDED DATA.
C     SECOND PART: CIMINELLI AND BEARPARK SEARCH ALGORITHMS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     ARRAY     AVAILABLE BUFFER (S).
C     LM5       DIMENSION OF ARRAY (I).
C     ICALL     CONTROL AND ERROR FLAG (I,O).
C     SCFCAL    EXTERNAL ROUTINE FOR ENERGY EVALUATION (I).
C     *
C     OUTPUT OF RESULTS VIA COMMON BLOCK ERG.
C     E         TARGET FUNCTION TO BE MINIMIZED (O).
C     G         GRADIENT OF TARGET FUNCTION (O).
C     GNORM     ASSOCIATE GRADIENT NORM (O).
C     *
C     ICALL DEFINES THE TYPE OF CALCULATION IN SUBROUTINE SCF.
C     SEE SUBROUTINE SCF FOR THE LIST OF POSSIBLE VALUES.
C     SEE SUBROUTINE SCFMUL FOR SPECIAL CONVENTIONS.
C     *
C     ICALL  ALSO SERVES AS AN ERROR FLAG FOR SUBROUTINE SCFCAL.
C      -1    RETURNED FROM SCFCAL IF THERE IS NO SCF CONVERGENCE
C     *
C     REFERENCES FOR CONICAL INTERSECTION SEARCH.
C     (1) C. CIMINELLI, G. GRANUCCI, AND M. PERSICO, CHEM. EUR. J.
C         10, 2327-2341 (2004).
C     (2) M. J. BEARPARK, M. A. ROBB, AND H. B. SCHLEGEL, CHEM. PHYS.
C         LETT. 223, 269-274 (1994).
C     (3) R. IZZO AND M. KLESSINGER, J. COMPUT. CHEM. 21, 52-62 (2000).
C     (4) A. TONIOLO, M. BEN-NUN, AND T. J. MARTINEZ, J. PHYS. CHEM. A
C         106, 4679-4689 (2002).
C     (5) M. R. MANAA AND D. R. YARKONY, J. CHEM. PHYS. 99, 5251-5256
C         (1993).
C     (6) D. R. YARKONY, J. PHYS. CHEM. A 108, 3200-3205 (2004).
C     *
      USE LIMIT, ONLY: LM1, LMV, LM1M, LMPROP, LMSTAT, LMGRD, LMNAC, 
     1                 LMYL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     CHARACTER*(*) CIFILE
C
      EXTERNAL SCFCAL
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CGRAD / CG(3,LM1+LM1M)
     ./CGRADM/ CGX(3,LM1+LM1M,LMGRD),CNORMX(LMGRD)
     ./CIPRP / CIPROP(LMPROP,LMSTAT),ICISYM(LMSTAT)
     ./CNACV / CC(3,LM1+LM1M),CCNORM,CCSUM
     ./CONINF/ E0,E1,GF2NRM
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DFP   / X(LMV),NVAR
     ./ERG   / E,G(LMV),GNORM,CNORM
     ./GRDORG/ IGRST(LMGRD),ISTATE,JSTATE
     ./INOPT2/ IN2(300)
     ./INOPT4/ XN4(50)
     ./NBFILE/ NBF(20)
     ./YARCON/ YLGCV(LMYL),YLGCT(LMYL),IYLGC(5,LMYL)
     ./YARLAG/ YLVAL(LMYL),YLGRD(LMYL),YLD(LMYL),NYL
     ./YARUPD/ YARV(3,LM1+LM1M,LMYL),YARS(LMV,2),YARG(LMV,LMYL,2),IXTRAP
      DIMENSION ARRAY(LM5)
      DIMENSION ENERGX(LMGRD),F2(3*NUMAT)
      DIMENSION GF2(LMV)
      DIMENSION OFTMP(3,LM1+LM1M)
C      DIMENSION GCHECK(LMV)
      DATA DELTA/1.0D-05/
      DATA TINY/1.0D-10/
      DATA ICICAL/0/
      SAVE ICICAL
      SAVE EREF0,EREF1
C *** FILE NUMBERS.
      NB6    = NBF(6)
      NB7    = NBF(7)
C *** INPUT OPTIONS.
      NSAV13 = IN2(17)
      NSAV15 = IN2(18)
      IPRINT = IN2(38)
      MPRINT = IN2(41)
      NPRINT = IN2(72)
      NUMATM = IN2(120)
      KEEPCI = IN2(158)
      NCIGRD = IN2(159)
      ICROSS = IN2(160)
C *** INITIALIZATION.
      ICALL1 = ICALL/10
      ICALL2 = ICALL-10*ICALL1
      IREFE  = ISTATE
      I3N    = 3*NUMAT
      NUMALL = NUMAT+NUMATM
C     FOR USE WITH SCFSAV
      IF(ICALL2.EQ.0) THEN
         ICALLS = 0
      ELSE
C        CARTESIAN AND 'INTERNAL' GRADIENTS AVAILABLE
         ICALLS = 1
      ENDIF
C
C *** LOOP OVER ALL CI CALCULATIONS.
      CALL SCFMUL (ARRAY,LM5,ICALL,SCFCAL,ENERGX,CGX,CNORMX,NUMALL)
      IF(ICALL.EQ.-1) RETURN
      ICICAL = ICICAL+1
C
C *** CONICAL INTERACTION SEARCH.
      E0 = ENERGX(1)
      E1 = ENERGX(2)
      IF(IREFE.EQ.0) THEN
         EREF0 = E0
         EREF1 = E1
         IF(IPRINT.GT.-5) WRITE(NB6,600) EREF0,IGRST(1),EREF1,IGRST(2)
      ENDIF
C     CALPHA, CBETA CORRESPOND TO THE FOLLOWING:
C     C1, C2 FOR CIMINELLI
C     C3, C4 FOR BEARPARK
C     T1, T2 FOR YARKONY
      CALPHA = XN4(25)
      CBETA  = XN4(26)
      GAPCON = XN4(27)
C     IF(E1.LT.EREF1 .AND. DABS(E1-E0).LT.GAPCON) THEN
C        WRITE(CIFILE,610) ICICAL
C        OPEN(NB7,FILE=CIFILE(1:16),FORM='FORMATTED',STATUS="UNKNOWN")
C        CALL GEOSAV(NUMAT,5,NB7)
C        CLOSE(NB7)
C        IF(IPRINT.GT.0) WRITE(NB6,620) ICICAL,CALPHA,CBETA,E0,E1
C      ENDIF
C
C *** CIMINELLI ALGORITHM (REFERENCE 1).
C     THE TARGET FUNCTION E TO BE MINIMIZED IS THE AVERAGE (E1+E0)/2
C     OF THE TWO STATE ENERGIES PLUS A PENALTY FUNCTION THAT INCREASES
C     MONOTONOUSLY WITH THE ENERGY DIFFERENCE (E1-E0), EQ.(3) OF REF.1.
C     THE GRADIENT OF THE TARGET FUNCTION IS THE AVERAGE (G1+G0)/2 
C     OF THE CORRESPONDING TWO ENERGY GRADIENTS PLUS THE DERIVATIVE
C     GSCALE*(G1-G0) OF THE PENALTY TERM.
C     NOTE: G1 STORED IN CGX STARTING AT CGX(1,1,2).
C     NOTE: G0 STORED IN CGX STARTING AT CGX(1,1,1).
      IF(ICROSS.EQ.3) THEN
         WDIFF  = ONE+((E1-E0)*(E1-E0))/(CBETA*CBETA)
         E      = PT5*(E0+E1)+CALPHA*CBETA*CBETA*DLOG(WDIFF)
         IF(ICALL2.EQ.0) GO TO 100
C        GRADIENT OF THE TARGET FUNCTION.
         GSCALE = TWO*CALPHA*(E1-E0)/WDIFF
         CALL DSCAL (I3N,PT5-GSCALE,CGX,1)
         CALL DAXPY (I3N,PT5+GSCALE,CGX(1,1,2),1,CGX,1)
         CALL GEOGRD (CGX,G,NUMAT,NVAR)
C        GRADIENT NORMS
         CNORM  = SQRT(DDOT(3*NUMAT,CGX,1,CGX,1))
         GNORM  = SQRT(DDOT(NVAR,G,1,G,1))
C
C *** BEARPARK ALGORITHM (REFERENCES 2 AND 3).
      ELSE IF(ICROSS.EQ.4) THEN
         IF(ICALL2.EQ.0) GO TO 90
C        FIRST WE DETERMINE THE GRADIENT OF THE TARGET FUNCTION.
C        THE REQUIRED NONADIABATIC COUPLING VECTOR BETWEEN THE TWO
C        STATES ISTATE AND JSTATE IS DETERMINED BY CALLING SCF WITH
C        THE FOLLOWING SPECIAL CONVENTIONS.
C        (1) ICALL=ICALL2=3     (USED IN SCFCAL, LOCAL)
C        (2) IN2(160)=ICROSS+10 (USED IN PSDST1 AND PSOMX, RESET).
         ISTATE = IGRST(1)
         JSTATE = IGRST(2)
         ICALL2 = 3
         IN2(160) = ICROSS+10
         CALL SCF (ARRAY,LM5,ICALL2,SCFCAL)
         IN2(160) = ICROSS
         IF(ICALL2.EQ.-1) RETURN
         CALL DCOPY (3*NUMALL,CG,1,CC,1)
CWT      THE ORIGINAL CODE SCALED THE COUPLING VECTOR BY HBAKMA (H-BAR).
CWT      THIS IS NOT DONE HERE BECAUSE THE GRADIENT DIFFERENCE VECTOR
CWT      AND THE COUPLING VECTOR HAVE THE SAME DIMENSION, SEE EQ.(1)
CWT      AND EQ.(2) OF REF.2, AND ARE COMPUTED CORRESPONDINGLY (SCF).
CWT      THE COUPLING VECTOR CC WILL BE NORMALIZED ANYWAY (SEE BELOW).
CWT      CALL DSCAL (3*NUMALL,HBAKMA,CC,1)
C
C        COMPUTE GRADIENT DIFFERENCE VECTOR X1=G1-G0, EQ.(1) OF REF.2.
C        AFTER CALL DAXPY: G1-G0 STORED IN CGX STARTING AT CGX(1,1,1).
         CALL DSCAL (I3N,-ONE,CGX,1)
         CALL DAXPY (I3N,ONE,CGX(1,1,2),1,CGX,1)
C        NORMALIZE X1=G1-G0 AS INDICATED BY EQ.(4) OF REF.2.
C        NORMALIZATION IS ACHIEVED BY SCALING WITH THE INVERSE NORM.
C        (BUT THIS IS NOT ATTEMPTED IF NORM IS ZERO) 
         DNORM = DNRM2 (I3N,CGX,1)
         IF(DNORM.GT.TINY) THEN
            CALL DSCAL (I3N,ONE/DNORM,CGX,1)
         ENDIF
C
C        COMPUTE PROJECTION OF X2=CC ALONG CGX.
         DPROJ = DDOT (I3N,CC,1,CGX,1)
C        ORTHOGONALIZE CC WITH RESPECT TO CGX: CC - <CC|CGX> CGX
C        FIRST ORTHOGONALIZE AND THEN NORMALIZE CC (SEE REF.3).
         CALL DAXPY (I3N,-DPROJ,CGX,1,CC,1)
C        NORMALIZE NONADIABATIC COUPLING VECTOR CC, EQ.(2) OF REF.2.
         DNORM = DNRM2 (I3N,CC,1)
         IF(DNORM.GT.TINY) THEN
            CALL DSCAL (3*NUMAT,ONE/DNORM,CC,1)
         ENDIF
C
C        PROJECTION P OF THE UPPER-STATE GRADIENT G1 ONTO THE
C        N-2 ORTHOGONAL COMPLEMENT OF THE PLANE SPANNED BY X1,X2.
C        PROJECTION MATRIX ACCORDING TO EQ.(3) OF REF.3:
C        P = 1 - X1*X1(TRANSPOSED) - X2*X2(TRANSPOSED)
C        F2 = P*G1 = G1 - X1*X1(TRANSPOSED)*G1 - X2*X2(TRANSPOSED)*G1
C        THE RHS VECTOR MULTIPLICATIONS CONTAIN DOT PRODUCTS:
C        DDOT1 = X1(TRANSPOSED)*G1
C        DDOT2 = X2(TRANSPOSED)*G1
C        HENCE THE PROJECTED GRADIENT F2 IS THE LINEAR COMBINATION:
C        F2 = G1 - X1*DDOT1 - X2*DDOT2
C        NOTE: G1 STORED IN CGX STARTING AT CGX(1,1,2).
C        NOTE: X1 STORED IN CGX STARTING AT CGX(1,1,1).
C        NOTE: X2 STORED IN CC  STARTING AT CC(1,1).
CWT      NOTE: ORIGINAL CODE REPLACED - ERROR IN INNER LOOP.
         DDOT1 = DDOT(I3N,CGX,1,CGX(1,1,2),1)
         DDOT2 = DDOT(I3N,CC ,1,CGX(1,1,2),1)
         CALL DCOPY (I3N,CGX(1,1,2),1,F2,1)
         CALL DAXPY (I3N,-DDOT1,CGX,1,F2,1)
         CALL DAXPY (I3N,-DDOT2,CC ,1,F2,1)
C
C        CALCULATE PROJECTED GRADIENT USED FOR CONVERGENCE CRITERIA
C        (REF. 3) AND TESTING FOR HESSIAN RESETS (REF. 4)
         CALL GEOGRD (F2,GF2,NUMAT,NVAR)
         GF2NRM = DNRM2(NVAR,GF2,1)
C
C        AT THIS POINT WE HAVE ASSEMBLED THE PROJECTED GRADIENT F2,
C        SEE EQ.(5) OF REF.2 AND EQ.(2) OF REF.3.
C        THE CHOICE OF TARGET FUNCTION E AND ASSOCIATE GRADIENT G
C        IS SLIGHTLY DIFFERENT IN REFERENCES 2 AND 3.
C        WE FOLLOW THE MORE EXPLICIT PROCEDURE OF REFERENCE 3.
C        CALPHA = ALPHA0 AND CBETA = ALPHA1 IN EQ.(4) OF REF.3.
C        TARGET GRADIENT DEFINED ACCORDING TO EQ.(4) OF REF.3.
C        CARTESIAN TARGET GRADIENT STORED IN CGX STARTING AT CGX(1,1,1).
C        AFTER CALL GEOGRD: OUTPUT TARGET GRADIENT STORED IN G.
         F2FAC = CALPHA*(ONE-CBETA)
         CALL DSCAL (I3N,TWO*(E1-E0)*CALPHA*CBETA,CGX,1)
         CALL DAXPY (I3N,F2FAC,F2,1,CGX,1)
         CALL GEOGRD (CGX,G,NUMAT,NVAR)

C        GRADIENT NORMS
         CNORM  = SQRT(DDOT(3*NUMAT,CGX,1,CGX,1))
         GNORM  = SQRT(DDOT(NVAR,G,1,G,1))

C        THE 'ENERGY' IS DEFINED HERE AS THE ENERGY GAP
C        THIS QUANTITY IS NOT MINIMISED DIRECTLY (SEE REF. 4)
C        BUT IT IS USED AS A TEST FOR RESETTING THE HESSIAN IF IT
C        DROPS BY A LARGE AMOUNT.
C        NOTE: MUST FOLLOW CALL TO SCF() OR E WILL BE OVERWRITTEN 
  90     E = E1-E0
         
C *** YARKONY ALGORITHM
      ELSE IF (ICROSS.EQ.5) THEN

         IF(ICALL2.EQ.0) GO TO 95

C        FIRST DEFINE THE GRADIENT OF THE TARGET FUNCTION
C        RIGHT HAND SIDE OF EQ.(2.5) OF REF.5.

C        CALCULATE THE GRADIENT OF THE INTERSTATE COUPLING.
C        THIS QUANTITY IS RELATED TO THE NONADIABATIC COUPLING VECTOR
C        BY EQ.(2.10) OF REF.5. 
C        DUE TO THE APPROXIMATIONS MADE IN MNDO, THE SAME PROCEDURE 
C        AS THE BEARPARK ALGORTIHM CAN BE USED, BECAUSE
C        h (OF EQ.2.10) IS RETURNED, NOT f.
C        THIS IS NOT AN APPROXIMATION IN THE CONTEXT OF THE YARKONY 
C        ALGORITHM - h IS THE DESIRED QUANTITY
         ISTATE = IGRST(1)
         JSTATE = IGRST(2)
         ICALL2 = 3
         IN2(160) = ICROSS+10
         CALL SCF (ARRAY,LM5,ICALL2,SCFCAL)
         IN2(160) = ICROSS
         IF(ICALL2.EQ.-1) RETURN
         CALL DCOPY (3*NUMALL,CG,1,CC,1)

C        COMPUTE GRADIENT DIFFERENCE VECTOR = G0-G1
         CALL DSCAL (I3N,-ONE,CGX(1,1,2),1)
         CALL DAXPY (I3N,ONE,CGX,1,CGX(1,1,2),1)

C        START OF EXTRAPOLATABLE FUNCTION SECTION (SEE REF. 6).

C        USE (G0+G1)/2 INSTEAD OF G0
C        AS THE LATTER IS NOT WELL-BEHAVED ON THE SEAM
         CALL DAXPY (I3N,-PT5,CGX(1,1,2),1,CGX,1)

C        STORE (G0+G1)/2 FOR HESSIAN UPDATE IN UPDHES
         CALL GEOGRD (CGX(1,1,1),YARS(1,1),NUMAT,NVAR)

C        CHECK IF SEAM HAS BEEN LOCATED
C        IF SO, THE GRADIENT DIFFERENCE VECTOR AND INTERSTATE
C        COUPLING GRADIENT CAN BE ORTHOGONALISED
         IF(IXTRAP.EQ.0 .AND. ABS(E0-E1).LT.CALPHA) THEN
            IXTRAP = 1
            WRITE(NB6,625)
         ENDIF
C        TURN ORTHOGONALISATION ROUTINE OFF AGAIN IF THE
C        SEAM IS LOST
         IF(IXTRAP.GT.0 .AND. ABS(E0-E1).GT.CBETA) THEN
            IXTRAP = 0
            WRITE(NB6,626)
         ENDIF

C        START OF ORTHOGONALISATION ROUTINE
         IF(IXTRAP.GT.0) THEN

C           FIRST, HALVE THE GRADIENT DIFFERENCE VECTOR
C           AS IMPLIED BY COMPARING EQ. 7 AND EQ. 9 OF REF. 6.
            CALL DSCAL (I3N,PT5,CGX(1,1,2),1)
            CALL DSCAL (I3N,PT5,YARV(1,1,1),1)

C           CALCULATE ORTHOGONALISING ANGLE USING EQ. 9 OF REF. 6
            DSPG   = DDOT (I3N,CGX(1,1,2),1,CGX(1,1,2),1)
            DSPH   = DDOT (I3N,CC,1,CC,1)
            DSP2HG = TWO*DDOT(I3N,CC,1,CGX(1,1,2),1)
            TAN4TH = -DSP2HG/(DSPH-DSPG)
            TWOTHT = ATAN(TAN4TH)/TWO
            CTWO = COS(TWOTHT)
            STWO = SIN(TWOTHT)

C           ORTHOGONALISE THE VECTORS (EQ. 8)
            CALL DCOPY(I3N,CGX(1,1,2),1,OFTMP,1)
            CALL DSCAL(I3N,CTWO,OFTMP,1)
            CALL DAXPY(I3N,STWO,CC,1,OFTMP,1)
            CALL DSCAL(I3N,CTWO,CC,1)
            CALL DAXPY(I3N,-STWO,CGX(1,1,2),1,CC,1)
            CALL DCOPY(I3N,OFTMP,1,CGX(1,1,2),1)

C           EQ. 9 IS ONLY UNIQUE UP TO TRANSPOSITIONS AND SIGN CHANGES.
C           THESE MUST BE TESTED FOR AND SWAPPED AS APPROPRIATE OTHERWISE
C           THE ORTHOGONALISED VECTORS WILL NOT BE SLOWLY VARYING.

C           TEST FOR TRANSPOSITIONS AND SIGN CHANGES USING OVERLAP
C           CRITERION. THIS MAY NOT BE THE MOST EFFICIENT WAY OF TESTING,
C           BUT NO CLUES ARE GIVEN IN REF. 6 TO YARKONY'S OWN PROCEDURE.
            DNRMG  = SQRT(DDOT(I3N,CGX(1,1,2),1,CGX(1,1,2),1))
            DNRMH  = SQRT(DDOT(I3N,CC,1,CC,1))
            DNRMG0 = SQRT(DDOT(I3N,YARV(1,1,1),1,YARV(1,1,1),1))
            DNRMH0 = SQRT(DDOT(I3N,YARV(1,1,2),1,YARV(1,1,2),1))
            
            DSPGG0 = DDOT(I3N,CGX(1,1,2),1,YARV(1,1,1),1)
            DSPHG0 = DDOT(I3N,CC,1,YARV(1,1,1),1)
            DSPGH0 = DDOT(I3N,CGX(1,1,2),1,YARV(1,1,2),1)
            DSPHH0 = DDOT(I3N,CC,1,YARV(1,1,2),1)

            CGG0   = DSPGG0/(DNRMG*DNRMG0)
            CHG0   = DSPHG0/(DNRMH*DNRMG0)
            CGH0   = DSPGH0/(DNRMG*DNRMH0)
            CHH0   = DSPHH0/(DNRMH*DNRMH0)

            IF(ABS(CHG0).GT.ABS(CGG0)) THEN
C              TRANSPOSITION DETECTED
               ITRANS = 1
               SCG = CHG0
               SCH = CGH0
            ELSE
C              NO TRANSPOSITION DETECTED
               ITRANS = 0
               SCG = CGG0
               SCH = CHH0
            ENDIF

            IF(SCG.LT.0) THEN
C              SIGN CHANGE IN GRAD DIFF VECTOR
               SIGNG = -ONE
            ELSE
C              NO SIGN CHANGE IN GRAD DIFF VECTOR
               SIGNG = ONE
            ENDIF

            IF(SCH.LT.0) THEN
C              SIGN CHANGE IN GRADIENT OF INTERSTATE COUPLING
               SIGNH = -ONE 
            ELSE
C              NO SIGN CHANGE IN GRADIENT OF INTERSTATE COUPLING
               SIGNH = ONE
            ENDIF

C           STORE VECTORS IN YARV WITH APPROPRIATE TRANSPOSITIONS/
C           SIGN CHANGES
            IF(ITRANS.EQ.0) THEN
               CALL DSCAL(I3N,SIGNG,CGX(1,1,2),1)
               CALL DSCAL(I3N,SIGNH,CC,1)
               CALL DCOPY(I3N,CGX(1,1,2),1,YARV(1,1,1),1)               
               CALL DCOPY(I3N,CC,1,YARV(1,1,2),1)
            ELSE
               CALL DSCAL(I3N,SIGNG,CC,1)
               CALL DSCAL(I3N,SIGNH,CGX(1,1,2),1)
               CALL DCOPY(I3N,CC,1,YARV(1,1,1),1)
               CALL DCOPY(I3N,CGX(1,1,2),1,YARV(1,1,2),1)               
            ENDIF

C           DOUBLE GRADIENT DIFFERENCE VECTOR FOR USE IN EQ. 7 OF REF. 6
            CALL DSCAL (I3N,TWO,YARV(1,1,1),1)

C           CHECK ORTHOGONALITY
C           OTEST = DDOT(I3N,YARV(1,1,2),1,YARV(1,1,1),1)
C           WRITE(NB6,*) 'ORTHOGONALITY TEST: ',OTEST

         ELSE

C           IF WE ARE NOT ON THE SEAM,
C           COPY UNORTHOGONALISED VECTORS DIRECTLY TO YARV
            CALL DCOPY(I3N,CGX(1,1,2),1,YARV(1,1,1),1)
            CALL DCOPY(I3N,CC,1,YARV(1,1,2),1)

         ENDIF

C        END OF EXTRAPOLATABLE FUNCTION SECTION

C        CONVERT GRAD DIFF VECTOR AND GRADIENT OF THE INTERSTATE
C        COUPLING VECTOR TO A FORMAT SUITABLE FOR USE IN 
C        HESSIAN UPDATE ROUTINES UPDHES AND YARHES
         CALL GEOGRD (YARV(1,1,1),YARG(1,1,1),NUMAT,NVAR)
         CALL GEOGRD (YARV(1,1,2),YARG(1,2,1),NUMAT,NVAR)

C        ADD SCALED G0-G1 TO G0 STARTING AT CGX(1,1,1)
         CALL DAXPY (I3N,YLVAL(1),YARV(1,1,1),1,CGX,1)

C        ADD SCALED INTERSTATE COUPLING VECTOR
         CALL DAXPY (I3N,YLVAL(2),YARV(1,1,2),1,CGX,1)

C        CALCULATE GEOMETRICAL CONSTRAINT CONTRIBUTIONS
C        TO THE GRADIENT, IF ANY
         IF(NYL.GT.2) THEN
            DO 120 I=1,(NYL-2)
               IF(IYLGC(1,I).EQ.1) THEN
                  DO 130 J=1,3
                     IATM = IYLGC(2,I)
                     IATN = IYLGC(3,I)
                     WMN = (COORD(J,IATM)-COORD(J,IATN))
                     YARV(J,IATM,I+2) = TWO*WMN
                     YARV(J,IATN,I+2) = -TWO*WMN
 130              CONTINUE
               ELSE IF(IYLGC(1,I).EQ.2) THEN
C                 CALCULATE ANGLE GRADIENTS NUMERICALLY
                  CALL BANGLE(IYLGC(2,I),IYLGC(3,I),IYLGC(4,I),
     1                        ANGO)
                  ANGO = ANGO/AFACT
                  DO 150 K=2,4
                  DO 140 J=1,3
                     IATM = IYLGC(K,I)
                     COORD(J,IATM) = COORD(J,IATM) + DELTA
                     CALL BANGLE(IYLGC(2,I),IYLGC(3,I),IYLGC(4,I),
     1                        ANGN)
                     ANGN = ANGN/AFACT
                     WMN = ((ANGN*ANGN-ANGO*ANGO)/DELTA)
                     YARV(J,IATM,I+2) = WMN
                     COORD(J,IATM) = COORD(J,IATM) - DELTA
 140              CONTINUE
 150              CONTINUE
               ELSE IF(IYLGC(1,I).EQ.3) THEN
C                 CALCULATE DIHEDRAL GRADIENTS NUMERICALLY
                  CALL DIHED(IYLGC(2,I),IYLGC(3,I),IYLGC(4,I),
     1                       IYLGC(5,I),ANGO,SINANG)
                  IF((YLGCT(I)-180.0D0).GT.ANGO) THEN
                     ANGO = ANGO + 360.0D0
                  ENDIF
                  ANGO = ANGO/AFACT
                  DO 170 K=2,5
                  DO 160 J=1,3
                     IATM = IYLGC(K,I)
                     COORD(J,IATM) = COORD(J,IATM) + DELTA
                     CALL DIHED(IYLGC(2,I),IYLGC(3,I),IYLGC(4,I),
     1                        IYLGC(5,I),ANGN,SINANG)
                     IF((YLGCT(I)-180.0D0).GT.ANGN) THEN
                        ANGN = ANGN + 360.0D0
                     ENDIF
                     ANGN = ANGN/AFACT
                     WMN = ((ANGN*ANGN-ANGO*ANGO)/DELTA)
                     YARV(J,IATM,I+2) = WMN
                     COORD(J,IATM) = COORD(J,IATM) - DELTA
 160              CONTINUE
 170              CONTINUE
               ENDIF
 120        CONTINUE

C        CONVERT SCALED CONSTRAINT CONTRIBUTIONS
C        TO A FORMAT SUITABLE FOR UPDHES AND YARHES
            DO 190 I=3,NYL
               CALL GEOGRD (YARV(1,1,I),YARG(1,I,1),NUMAT,NVAR)
 190        CONTINUE

C        ADD SCALED CONSTRAINT CONTRIBUTIONS TO THE GRADIENT
            DO 200 I=3,NYL
               CALL DAXPY (I3N,YLVAL(I),YARV(1,1,I),1,CGX,1)
 200        CONTINUE
         ENDIF

C        COPY LAGRANGIAN GRADIENT INTO G USING APPROPRIATE COORDINATES 
         CALL GEOGRD (CGX,G,NUMAT,NVAR)

C        GRADIENTS CORRESPONDING TO GRAD. DIFFERENCE VECTOR
C        AND GRADIENT OF THE INTERSTATE COUPLING VECTOR. SEE EQ. (2.5) 
C        OF REF. 5
C
C        IN REF. 6, YARKONY DOES NOT INDICATE WHAT SHOULD BE DONE WITH
C        THE RESIDUAL ENERGY DIFFERENCE WHEN ORTHOGONALISATION IS 
C        SWITCHED ON. IN THE ORIGINAL MNDO IMPLEMENTATION (AS REPORTED
C        IN THE IMPLEMENTATION PAPER), IT WAS
C        TRANSPOSED (AND HALVED) AND ITS SIGN CHANGED TO BE CONSISTENT
C        WITH ANY SUCH CHANGES OF THE VECTORS. THIS WAS AN EMPIRICAL 
C        FINDING THAT IMPROVED CONVERGENCE BEHAVIOUR.
C
C        CURRENT IMPROVED IMPLEMENTATION (FOLLOWING DL-FIND):
C        If orthogonalisation is switched on, the orthogonalisation 
C        procedure should also be applied to the final two gradients:
C        - Halve the first gradient (energy difference), like g was 
C          halved.
C        - Apply the rotation matrix from Eq. 17, i.e.
C
C        |  cos(2th)  sin(2th) | | (1/2)dE |  =  |  (1/2).cos(2th).dE |
C        | -sin(2th)  cos(2th) | |    0    |     | -(1/2).sin(2th).dE |
C
C        - If g and h were transposed, or there were sign changes, 
C          apply the same corrections to these gradients.
C        - Double the first gradient again.
C
C        IN PRINCIPLE, A THRESHOLD IS NO LONGER REQUIRED FOR 
C        ORTHOGONALISATION AS THE CONSTRAINT EQUATIONS NOW HOLD FOR
C        ANY VALUE OF THE ENERGY GAP. HOWEVER, THE THRESHOLD IS STILL
C        USEFUL TO PREVENT G AND H INFORMATION ENTERING THE HESSIAN
C        TOO EARLY, WHICH COULD FORCE THE STRUCTURE TOWARDS AN 
C        UNWANTED HIGH ENERGY SEAM.
C 
         IF(IXTRAP.EQ.0) THEN
            YLGRD(1) = (E0-E1)
            YLGRD(2) = ZERO
         ELSE IF(ITRANS.EQ.1) THEN
            YLGRD(1) = -ONE * STWO * (E0-E1) * SIGNG
            YLGRD(2) = PT5 * CTWO * (E0-E1) * SIGNH
         ELSE
            YLGRD(1) = CTWO * (E0-E1) * SIGNG
            YLGRD(2) = -PT5 * STWO * (E0-E1) * SIGNH
         ENDIF

C        GRADIENTS CORRESPONDING TO GEOMETRICAL CONSTRAINTS
         IF(NYL.GT.2) THEN
            DO 210 I=1,(NYL-2)
C              CALCULATE CURRENT VALUE OF CONSTRAINED PARAMETER
               IF(IYLGC(1,I).EQ.1) THEN
                  J = IYLGC(2,I)
                  K = IYLGC(3,I)
                  YLGCV(I) = SQRT( (COORD(1,K)-COORD(1,J))**2
     1                            +(COORD(2,K)-COORD(2,J))**2
     2                            +(COORD(3,K)-COORD(3,J))**2 )
               ELSE IF(IYLGC(1,I).EQ.2) THEN
                  CALL BANGLE(IYLGC(2,I),IYLGC(3,I),IYLGC(4,I),
     1                        YLGCV(I))
               ELSE IF(IYLGC(1,I).EQ.3) THEN
                  CALL DIHED(IYLGC(2,I),IYLGC(3,I),IYLGC(4,I),
     1                       IYLGC(5,I),YLGCV(I),SINANG)
                  IF((YLGCT(I)-180.0D0).GT.YLGCV(I)) THEN
                     YLGCV(I) = YLGCV(I) + 360.0D0
                  ENDIF
               ENDIF            
C              CALCULATE GRADIENTS
               IF(IYLGC(1,I).EQ.1) THEN
                  YLGRD(I+2) = (YLGCV(I)*YLGCV(I))-(YLGCT(I)*YLGCT(I))
               ELSE IF(IYLGC(1,I).GE.2) THEN
                  ANGV = YLGCV(I)/AFACT
                  ANGT = YLGCT(I)/AFACT
                  YLGRD(I+2) = ANGV*ANGV - ANGT*ANGT
               ENDIF
 210        CONTINUE
         ENDIF

C        CALCULATE GRADIENT NORMS
         CNORM  = DDOT(3*NUMAT,CGX,1,CGX,1)
         GNORM  = DDOT(NVAR,G,1,G,1)
         YNORM  = DDOT(NYL,YLGRD,1,YLGRD,1)
         CNORM  = SQRT(CNORM+YNORM)
         GNORM  = SQRT(GNORM+YNORM)

C        CHECK OF G USING YARS AND YARG
C        GCHECK MUST BE DECLARED IF THIS IS UNCOMMENTED
C         CALL DCOPY(NVAR,YARS(1,1),1,GCHECK,1)
C         DO 220 I=1,NYL
C            CALL DAXPY(NVAR,YLVAL(I),YARG(1,I,1),1,GCHECK,1)
C 220     CONTINUE
C         GCHKN = DDOT(NVAR,GCHECK,1,GCHECK,1)
C         GCHKN = SQRT(GCHKN+YNORM)

C        THE 'ENERGY' IS DEFINED AS THE ENERGY GAP
C        FOR USE WITH TEST FOR RESETTING HESSIAN IN EIGF().
C        IT IS NOT OTHERWISE USED

 95      E = ABS(E1-E0)

      ENDIF
C
C *** SAVE 'ENERGY' AND 'GRADIENT' INFORMATION HERE
C     INSTEAD OF SCF.f, WHERE INDIVIDUAL STATE INFORMATION
C     WOULD BE SAVED/OVERWRITTEN REPEATEDLY.
C     
C     THIS ALLOWS CIMINELLI/BEARPARK 'GRADIENTS' TO BE USED
C     IN EXTERNAL OPTIMIZERS IF DESIRED, WITHOUT ANY CHANGES
C     NECESSARY TO THE FORMAT OF FILE NB15.
C     IN THE CASE OF YARKONY, ONLY THE GRADIENT OF THE LAGRANGIAN
C     CAN BE SAVED, WHICH IS NOT SUFFICIENT INFORMATION FOR A
C     YARKONY OPTIMIZATION.
C
C     HOWEVER, EXTERNAL OPTIMIZERS WOULD NORMALLY BE EXPECTED TO
C     INTERFACE VIA ICROSS=1/2.
C
      IF(ICALL2.NE.0) THEN
          CALL DCOPY(I3N,CGX,1,CG,1)         
      ENDIF

      CALL SCFSAV (NSAV13,NSAV15,ICALL1,ICALLS,ICALLS)
C
C *** PRINTING SECTION.
  100 CONTINUE
      IF(IPRINT.GE.0) THEN
         WRITE(NB6,630) ICICAL,E,E0,E1,(E1-E0)
         IF(ICALL2.NE.0) THEN
            IF(ICROSS.EQ.3) THEN
               WRITE(NB6,640) GNORM
               WRITE(NB6,650) CNORM
            ELSE IF(ICROSS.EQ.4) THEN
               WRITE(NB6,660) GNORM
               WRITE(NB6,670) CNORM
            ELSE IF(ICROSS.EQ.5) THEN
               WRITE(NB6,680) GNORM
               WRITE(NB6,690) CNORM
            ENDIF
         ENDIF
         IF(IPRINT.GE.2) THEN
            IF(ICROSS.EQ.3) THEN
               WRITE(NB6,700) E,CALPHA,CBETA
            ELSE IF(ICROSS.EQ.4) THEN
               WRITE(NB6,710) E,CALPHA,CBETA
            ELSE IF(ICROSS.EQ.5) THEN
               WRITE(NB6,720) E,CALPHA,CBETA
            ENDIF
         ENDIF
      ENDIF
      RETURN
  500 FORMAT (///1X,'TIME FOR MULTIPLE CI RUNS   ',F10.3,' SECONDS')
  510 FORMAT (// 1X,'I  STATE',7X,'ENERGY',10X,'CNORM',10X,'GNORM'/)
  520 FORMAT (   1X,I1,I5,3F15.5)
  600 FORMAT (// 1X,'REFERENCE ENERGIES FOR SEARCH TAKEN FROM ',
     1              'CALCULATION AT INITIAL GEOMETRY'
     2        /  1X,'EREF0 =',F12.5,' KCAL/MOL FOR CI STATE',I3,
     3        /  1X,'EREF1 =',F12.5,' KCAL/MOL FOR CI STATE',I3)
C 610 FORMAT ('ci_',I9,'.inp')
C 620 FORMAT (/ 1X,'NEW MOPAC-TYPE INPUT FILE WRITTEN (MODE=5) ',
C    1             'FOR THE CURRENT STEP WITH ICICAL =',I9,
C    2        / 1X,'CALPHA =',F5.1,3X,'CBETA =',F5.1,
C    3          3X,'E0 =',F10.1,3X,'E1=',F10.1)
  625 FORMAT (/ 1X,'CI SEAM FOUND - USING EXTRAPOLATABLE FUNCTIONS') 
  626 FORMAT (/ 1X,'CI SEAM LOST - EXTRAPOLATABLE FUNCTIONS OFF')  
  630 FORMAT (/ 1X,'CONINT: ICICAL =',I5,
     1          3X,'E =',F8.1,3X,'E0 =',F8.1,3X,'E1=',F8.1,
     2          3X,'E1-E0 =',F16.9)
  640 FORMAT (/ 1X,'CIMINELLI GRADIENT NORM          ',3X,F20.5)
  650 FORMAT (  1X,'CIMINELLI CARTESIAN GRADIENT NORM',3X,F20.5)
  660 FORMAT (/ 1X,'BEARPARK GRADIENT NORM           ',3X,F20.5)
  670 FORMAT (  1X,'BEARPARK CARTESIAN GRADIENT NORM ',3X,F20.5)
  680 FORMAT (/ 1X,'YARKONY GRADIENT NORM            ',3X,F20.5)
  690 FORMAT (  1X,'YARKONY CARTESIAN GRADIENT NORM  ',3X,F20.5)
  700 FORMAT (/ 1X,'TARGET FUNCTION FOR CIMINELLI SEARCH:',F15.5,
     1          1X,'KCAL/MOL,   COEFFICIENTS:',2F8.3)
  710 FORMAT (/ 1X,'TARGET FUNCTION FOR BEARPARK SEARCH:',F15.5,
     1          1X,'KCAL/MOL,   COEFFICIENTS:',2F8.3)
  720 FORMAT (/ 1X,'TARGET FUNCTION FOR YARKONY SEARCH:',F15.5,
     1          1X,'KCAL/MOL,   INITIAL LMS:',2F8.3)
      END
