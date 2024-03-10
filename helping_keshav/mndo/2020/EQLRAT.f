      SUBROUTINE EQLRAT(N,DIAG,E,E2IN,D,IND,IERR,E2)
C
C     AUTHORS -
C        THIS IS A MODIFICATION OF ROUTINE TQLRAT FROM EISPACK EDITION 3
C        DATED AUGUST 1983.
C        TQLRAT IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
C        ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C        THIS VERSION IS BY S. T. ELBERT (AMES LABORATORY-USDOE).
C        SOME MINOR MODIFICATIONS BY W. THIEL, JANUARY 1991.
C
C     PURPOSE -
C        THIS ROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C        TRIDIAGONAL MATRIX
C
C     METHOD -
C        RATIONAL QL
C
C     ON ENTRY -
C        N      - INTEGER
C                 THE ORDER OF THE MATRIX.
C        D      - W.P. REAL (N)
C                 CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C        E2     - W.P. REAL (N)
C                 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF
C                 THE INPUT MATRIX IN ITS LAST N-1 POSITIONS.
C                 E2(1) IS ARBITRARY.
C
C      ON EXIT -
C        D      - W.P. REAL (N)
C                 CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C                 ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C                 ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C                 THE SMALLEST EIGENVALUES.
C        E2     - W.P. REAL (N)
C                 DESTROYED.
C        IERR   - INTEGER
C                 0          FOR NORMAL RETURN,
C                 J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                            DETERMINED AFTER 30 ITERATIONS.
C
C     DIFFERENCES FROM EISPACK 3 - DUE TO S. T. ELBERT
C        G=G+B INSTEAD OF IF(G.EQ.0) G=B ; B=B/64
C        F77 BACKWARD LOOPS INSTEAD OF F66 CONSTRUCT
C        GENERIC INTRINSIC FUNCTIONS
C        ARRARY  IND  ADDED FOR USE BY EINVIT
C
C     EXTERNAL ROUTINES -
C        EPSLON
C        INTRINSIC--ABS, SIGN, SQRT
C
C     MODIFICATIONS BY WALTER THIEL, JANUARY 1991.
C        REPLACE SOME GO-TO-CONSTRUCTS BY BLOCK-IF-CONSTRUCTS
C        COSMETIC CHANGES IN COMMENTS ETC.
C
C     MODIFICATIONS BY WALTER THIEL, DECEMBER 1999.
C        USE CONTINUE STATEMENT TO END DO LOOPS
C        IN ORDER TO AVOID WARNINGS FROM F90 COMPILER.
C
      INTEGER I,J,L,M,N,II,L1,IERR
      INTEGER IND(N)
      DOUBLE PRECISION D(N),E(N),E2(N),DIAG(N),E2IN(N)
      DOUBLE PRECISION B,C,F,G,H,P,R,S,T,EPSLON
      DOUBLE PRECISION SCALE,ZERO,ONE,TWO
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (SCALE= 1.0D0/64.0D0)
C
      IERR = 0
      D(1)=DIAG(1)
      IND(1) = 1
      K = 0
      ITAG = 0
      IF (N .EQ. 1) GO TO 1001
      DO 100 I = 2, N
         D(I)=DIAG(I)
         E2(I-1) = E2IN(I)
  100 CONTINUE
      F = ZERO
      T = ZERO
      B = EPSLON(ONE)
      C = B *B
      B = B * SCALE
      E2(N) = ZERO
C     .......... MAIN LOOP ...................
      DO 290 L = 1, N
         H = ABS(D(L)) + ABS(E(L))
         IF (T .LT. H) THEN
            T = H
            B = EPSLON(T)
            C = B * B
            B = B * SCALE
         ENDIF
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
         M = L - 1
  110    M = M + 1
         IF (E2(M) .GT. C) GO TO 110
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS AN EXIT
C                FROM THE LOOP ..........
         IF (M .GT. K) THEN
            IF (M .NE. N) E2IN(M+1) = ZERO
            K = M
            ITAG = ITAG + 1
         ENDIF
         IF (M .EQ. L) GO TO 210
C        ....... ITERATE .......
         DO 205 J = 1, 30
C              .......... FORM SHIFT ..........
            L1 = L + 1
            S = SQRT(E2(L))
            G = D(L)
            P = (D(L1) - G) / (TWO * S)
            R = SQRT(P*P+ONE)
            D(L) = S / (P + SIGN(R,P))
            H = G - D(L)
            DO 140 I = L1, N
                D(I) = D(I) - H
  140       CONTINUE
            F = F + H
C              .......... RATIONAL QL TRANSFORMATION ..........
            G = D(M) + B
            H = G
            S = ZERO
            DO 200 I = M-1,L,-1
               P = G * H
               R = P + E2(I)
               E2(I+1) = S * R
               S = E2(I) / R
               D(I+1) = H + S * (H + D(I))
               G = D(I) - E2(I) / G   + B
               H = G * P / R
  200       CONTINUE
            E2(L) = S * G
            D(L) = H
C              .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST
            IF (H .EQ. ZERO) GO TO 210
            IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
            E2(L) = H * E2(L)
            IF (E2(L) .EQ. ZERO) GO TO 210
  205    CONTINUE
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
      IERR = L
      GO TO 1001
C     .......... CONVERGED ............
  210    P = D(L) + F
C           .......... ORDER EIGENVALUES ..........
         I = 1
         IF (L .EQ. 1) GO TO 250
            IF (P .LT. D(1)) GO TO 230
               I = L
C           .......... LOOP TO FIND ORDERED POSITION
  220          I = I - 1
               IF (P .LT. D(I)) GO TO 220
               I = I + 1
               IF (I .EQ. L) GO TO 250
  230       CONTINUE
            DO 240 II = L, I+1, -1
               D(II) = D(II-1)
               IND(II) = IND(II-1)
  240       CONTINUE
  250    CONTINUE
         D(I) = P
         IND(I) = ITAG
  290 CONTINUE
 1001 RETURN
      END
