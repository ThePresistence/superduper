      SUBROUTINE EINVIT(NM,N,D,E,E2,M,W,IND,Z,IERR,RV1,RV2,RV3,RV4,RV6)
C
C     AUTHORS-
C        THIS IS A MODIFICATION OF ROUTINE TINVIT FROM EISPACK EDITION 3
C        DATED AUGUST 1983.
C        TINVIT IS A TRANSLATION OF THE INVERSE ITERATION TECHNIQUE
C        IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C        HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C        THIS VERSION IS BY S. T. ELBERT (AMES LABORATORY-USDOE).
C        SOME MINOR MODIFICATIONS BY W. THIEL, JANUARY 1991.
C
C     PURPOSE -
C        THIS ROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
C        SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES.
C
C     METHOD -
C        INVERSE ITERATION.
C
C     ON ENTRY -
C        NM     - INTEGER
C                 MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C                 ARRAY PARAMETERS AS DECLARED IN THE CALLING ROUTINE
C                 DIMENSION STATEMENT.
C        N      - INTEGER
C        D      - W.P. REAL (N)
C                 CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C        E      - W.P. REAL (N)
C                 CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C                 IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C        E2     - W.P. REAL (N)
C                 CONTAINS THE SQUARES OF CORRESPONDING ELEMENTS OF E,
C                 WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
C                 E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
C                 THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE
C                 SUM OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST
C                 CONTAIN 0.0 IF THE EIGENVALUES ARE IN ASCENDING ORDER,
C                 OR 2.0 IF THE EIGENVALUES ARE IN DESCENDING ORDER.
C                 IF TQLRAT, BISECT, TRIDIB, OR IMTQLV
C                 HAS BEEN USED TO FIND THE EIGENVALUES, THEIR
C                 OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE.
C        M      - INTEGER
C                 THE NUMBER OF SPECIFIED EIGENVECTORS.
C        W      - W.P. REAL (M)
C                 CONTAINS THE M EIGENVALUES IN ASCENDING
C                 OR DESCENDING ORDER.
C        IND    - INTEGER (M)
C                 CONTAINS IN FIRST M POSITIONS THE SUBMATRIX INDICES
C                 ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C                 1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX
C                 FROM THE TOP, 2 FOR THOSE BELONGING TO THE SECOND
C                 SUBMATRIX, ETC.
C        IERR   - INTEGER (LOGICAL UNIT NUMBER)
C                 LOGICAL UNIT FOR ERROR MESSAGES
C
C     ON EXIT -
C        ALL INPUT ARRAYS ARE UNALTERED.
C        Z      - W.P. REAL (NM,M)
C                 CONTAINS THE ASSOCIATED SET OF ORTHONORMAL
C                 EIGENVECTORS. ANY VECTOR WHICH WHICH FAILS TO CONVERGE
C                 IS LEFT AS IS (BUT NORMALIZED) WHEN ITERATING STOPPED.
C        IERR   - INTEGER
C                  0      FOR NORMAL RETURN,
C                 -R      IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C                         EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS.
C                         (ONLY LAST FAILURE TO CONVERGE IS REPORTED)
C
C        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
C
C        RV1    - W.P. REAL (N)
C                 DIAGONAL ELEMENTS OF U FROM LU DECOMPOSITION
C        RV2    - W.P. REAL (N)
C                 SUPER(1)-DIAGONAL ELEMENTS OF U FROM LU DECOMPOSITION
C        RV3    - W.P. REAL (N)
C                 SUPER(2)-DIAGONAL ELEMENTS OF U FROM LU DECOMPOSITION
C        RV4    - W.P. REAL (N)
C                 ELEMENTS DEFINING L IN LU DECOMPOSITION
C        RV6    - W.P. REAL (N)
C                 APPROXIMATE EIGENVECTOR
C
C     DIFFERENCES FROM EISPACK 3 -  DUE TO S. T. ELBERT
C        EPS3 IS SCALED BY  EPSCAL  (ENHANCES CONVERGENCE, BUT
C           LOWERS ACCURACY)!
C        ONE MORE ITERATION (MINIMUM 2) IS PERFORMED AFTER CONVERGENCE
C           (ENHANCES ACCURACY)!
C        REPLACE LOOP WITH PYTHAG WITH SINGLE CALL TO DNRM2!
C        IF NOT CONVERGED, USE PERFORMANCE INDEX TO DECIDE ON ERROR
C           VALUE SETTING, BUT DO NOT STOP!
C        L.U. FOR ERROR MESSAGES PASSED THROUGH IERR
C        USE PARAMETER STATEMENTS AND GENERIC INTRINSIC FUNCTIONS
C        USE LEVEL 1 BLAS
C        USE IF-THEN-ELSE TO CLARIFY LOGIC
C        LOOP OVER SUBSPACES MADE INTO DO LOOP.
C        LOOP OVER INVERSE ITERATIONS MADE INTO DO LOOP
C        ZERO ONLY REQUIRED PORTIONS OF OUTPUT VECTOR
C
C        EXTERNAL ROUTINES -
C        EPSLON
C        BLAS(1)--DASUM, DAXPY, DDOT, DNRM2, DSCAL
C        INTRINSIC FUNCTIONS - ABS, MAX, SQRT
C
C     MODIFICATIONS BY WALTER THIEL, JANUARY 1991.
C        REMOVE BLAS(1) REFERENCES--DASUM,DAXPY,DDOT,DSCAL
C           AND REPLACE THEM BY FORTRAN CODE
C           LOOK FOR 'CBLAS' TO SEE THE REPLACEMENTS
C        COSMETIC CHANGES IN COMMENT CARDS ETC.
C
C     MODIFICATIONS BY WALTER THIEL, DECEMBER 1999.
C        USE CONTINUE STATEMENT TO END DO LOOPS
C        IN ORDER TO AVOID WARNINGS FROM F90 COMPILER.
C
C     MODIFICATIONS BY AXEL KOSLOWSKI, JULY 2014.
C        EXPLICITLY INITIALIZE SCRATCH VECTORS WITH ZERO.
C
      LOGICAL CONVGD
      INTEGER GROUP,I,IERR,ITS,J,JJ,M,N,NM,P,Q,R,S,SUBMAT,TAG
      INTEGER IND(M)
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),Z(NM,M)
      DOUBLE PRECISION RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      DOUBLE PRECISION ANORM,EPS2,EPS3,EPS4,NORM,ORDER,RHO,U,UK,V
      DOUBLE PRECISION X0,X1,XU
      DOUBLE PRECISION EPSCAL,GRPTOL,HUNDRD,ONE,TEN,ZERO
      DOUBLE PRECISION EPSLON, ESTPI1, DNRM2
CBLAS DOUBLE PRECISION EPSLON, ESTPI1, DASUM, DDOT, DNRM2
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, GRPTOL = 0.001D0)
      PARAMETER (EPSCAL = 0.5D0, HUNDRD = 100.0D0, TEN = 10.0D0)
C
  001 FORMAT(' EIGENVECTOR ROUTINE EINVIT DID NOT CONVERGE FOR VECTOR'
     *      ,I5,'.  NORM =',1PE10.2,' PERFORMANCE INDEX =',E10.2/
     *      ' (AN ERROR HALT WILL OCCUR IF THE PI IS GREATER THAN 100)')
C
C-----------------------------------------------------------------------
C
C-AK  RV3(N-1) IS ACCESSED WITHOUT PRIOR INITIALIZATION IN ONE EXAMPLE
C-AK  THAT I ANALYZED. TO BE ON THE SAFE SIDE, ALL SCRATCH VECTORS ARE
C-AK  NOW INITIALIZED WITH ZERO.
      RV1(1:N) = ZERO
      RV2(1:N) = ZERO
      RV3(1:N) = ZERO
      RV4(1:N) = ZERO
      RV6(1:N) = ZERO
      LUEMSG = IERR
      IERR = 0
      X0 = ZERO
      UK = ZERO
      NORM = ZERO
      EPS2 = ZERO
      EPS3 = ZERO
      EPS4 = ZERO
      GROUP = 0
      TAG = 0
      ORDER = ONE - E2(1)
      Q = 0
      DO 930 SUBMAT = 1, N
         P = Q + 1
C        .......... ESTABLISH AND PROCESS NEXT SUBMATRIX ..........
         DO 120 Q = P, N-1
            IF (E2(Q+1) .EQ. ZERO) GO TO 140
  120    CONTINUE
         Q = N
C        .......... FIND VECTORS BY INVERSE ITERATION ..........
  140    CONTINUE
         TAG = TAG + 1
         ANORM = ZERO
         S = 0
         DO 920 R = 1, M
            IF (IND(R) .NE. TAG) GO TO 920
            ITS = 1
            X1 = W(R)
            IF (S .NE. 0) GO TO 510
C        .......... CHECK FOR ISOLATED ROOT ..........
            XU = ONE
            IF (P .EQ. Q) THEN
               RV6(P) = ONE
               CONVGD = .TRUE.
               GO TO 860
            END IF
            NORM = ABS(D(P))
            DO 500 I = P+1, Q
               NORM = MAX( NORM, ABS(D(I)) + ABS(E(I)) )
  500       CONTINUE
C        .......... EPS2 IS THE CRITERION FOR GROUPING,
C                   EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                   ROOTS ARE MODIFIED BY EPS3,
C                   EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW .........
            EPS2 = GRPTOL * NORM
            EPS3 = EPSCAL * EPSLON(NORM)
            UK = Q - P + 1
            EPS4 = UK * EPS3
            UK = EPS4 / SQRT(UK)
            S = P
            GROUP = 0
            GO TO 520
C        .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
  510       IF (ABS(X1-X0) .GE. EPS2) THEN
C                 ROOTS ARE SEPARATE
               GROUP = 0
            ELSE
C                 ROOTS ARE CLOSE
               GROUP = GROUP + 1
               IF (ORDER * (X1 - X0) .LE. EPS3) X1 = X0 + ORDER * EPS3
            END IF
C        .......... ELIMINATION WITH INTERCHANGES AND
C                   INITIALIZATION OF VECTOR ..........
  520       CONTINUE
            U = D(P) - X1
            V = E(P+1)
            RV6(P) = UK
            DO 550 I = P+1, Q
               RV6(I) = UK
               IF (ABS(E(I)) .GT. ABS(U)) THEN
C                 EXCHANGE ROWS BEFORE ELIMINATION
C                  *** WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
C                      E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY .......
                  XU = U / E(I)
                  RV4(I) = XU
                  RV1(I-1) = E(I)
                  RV2(I-1) = D(I) - X1
                  IF(I.NE.Q) RV3(I-1) = E(I+1)
                  U = V - XU * RV2(I-1)
                  V = -XU * RV3(I-1)
               ELSE
C                 STRAIGHT ELIMINATION
                  XU = E(I) / U
                  RV4(I) = XU
                  RV1(I-1) = U
                  RV2(I-1) = V
                  RV3(I-1) = ZERO
                  U = D(I) - X1 - XU * V
                  IF(I.NE.Q) V = E(I+1)
               END IF
  550       CONTINUE
            IF (ABS(U) .LE. EPS3) U = EPS3
            RV1(Q) = U
            RV2(Q) = ZERO
            RV3(Q) = ZERO
C     ........ DO INVERSE ITERATIONS
            CONVGD = .FALSE.
            DO 800 ITS = 1, 5
               IF (ITS .EQ. 1) GO TO 600
C                    .......... FORWARD SUBSTITUTION ..........
                  IF (NORM .EQ. ZERO) THEN
                     RV6(S) = EPS4
                     S = S + 1
                     IF (S .GT. Q) S = P
                  ELSE
                     XU = EPS4 / NORM
CBLAS                CALL DSCAL (Q-P+1, XU, RV6(P), 1)
                     DO 560 I = P, Q
                     RV6(I) = RV6(I) * XU
  560                CONTINUE
                  END IF
C                     ... ELIMINATION OPERATIONS ON NEXT VECTOR
                  DO 590 I = P+1, Q
                     U = RV6(I)
C                         IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                         WAS PERFORMED EARLIER IN THE
C                         TRIANGULARIZATION PROCESS ..........
                     IF (RV1(I-1) .EQ. E(I)) THEN
                        U = RV6(I-1)
                        RV6(I-1) = RV6(I)
                     ELSE
                        U = RV6(I)
                     END IF
                     RV6(I) = U - RV4(I) * RV6(I-1)
  590             CONTINUE
  600          CONTINUE
C           .......... BACK SUBSTITUTION
               RV6(Q) = RV6(Q) / RV1(Q)
               V = U
               U = RV6(Q)
               NORM = ABS(U)
               DO 620 I = Q-1, P, -1
                  RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
                  V = U
                  U = RV6(I)
                  NORM = NORM + ABS(U)
  620          CONTINUE
               IF (GROUP .EQ. 0) GO TO 700
C                 ....... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                         MEMBERS OF GROUP ..........
                  J = R
                  DO 680 JJ = 1, GROUP
  630                J = J - 1
                     IF (IND(J) .NE. TAG) GO TO 630
CBLAS                CALL DAXPY(Q-P+1, -DDOT(Q-P+1,RV6(P),1,Z(P,J),1),
CBLAS*                          Z(P,J),1,RV6(P),1)
                     XU = ZERO
                     DO 640 I = P, Q
                     XU = XU + RV6(I) * Z(I,J)
  640                CONTINUE
                     DO 660 I = P, Q
                     RV6(I) = RV6(I) - XU * Z(I,J)
  660                CONTINUE
  680             CONTINUE
CBLAS             NORM = DASUM(Q-P+1, RV6(P), 1)
                  NORM = ZERO
                  DO 720 I = P, Q
                  NORM = NORM + ABS(RV6(I))
  720             CONTINUE
  700          CONTINUE
               IF (CONVGD) GO TO 840
               IF (NORM .GE. ONE) CONVGD = .TRUE.
  800       CONTINUE
C        .......... NORMALIZE SO THAT SUM OF SQUARES IS
C                   1 AND EXPAND TO FULL ORDER ..........
  840       CONTINUE
            XU = ONE / DNRM2(Q-P+1,RV6(P),1)
  860       CONTINUE
            DO 870 I = 1, P-1
               Z(I,R) = ZERO
  870       CONTINUE
            DO 890 I = P,Q
               Z(I,R) = RV6(I) * XU
  890       CONTINUE
            DO 900 I = Q+1, N
               Z(I,R) = ZERO
  900       CONTINUE
            IF (.NOT.CONVGD) THEN
               RHO = ESTPI1(Q-P+1,X1,D(P),E(P),Z(P,R),ANORM)
               IF (RHO .GE. TEN  .AND.  LUEMSG .GT. 0)
     *             WRITE(LUEMSG,001) R,NORM,RHO
C               *** SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
               IF (RHO .GT. HUNDRD) IERR = -R
            END IF
            X0 = X1
  920    CONTINUE
         IF (Q .EQ. N) GO TO 940
  930 CONTINUE
  940 CONTINUE
      RETURN
      END
