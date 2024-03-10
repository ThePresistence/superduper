      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  Purpose
*  =======
*
*     DASUM takes the sum of the absolute values.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DABS,MOD
*     ..
      DASUM = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
*        code for increment equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,6)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DTEMP + DABS(DX(I))
            END DO
            IF (N.LT.6) THEN
               DASUM = DTEMP
               RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,6
            DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I+1)) +
     $              DABS(DX(I+2)) + DABS(DX(I+3)) +
     $              DABS(DX(I+4)) + DABS(DX(I+5))
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DTEMP = DTEMP + DABS(DX(I))
         END DO
      END IF
      DASUM = DTEMP
      RETURN
      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     DAXPY constant times a vector plus a vector.
*     uses unrolled loops for increments equal to one.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DY(I) + DA*DX(I)
            END DO
         END IF
         IF (N.LT.4) RETURN
         MP1 = M + 1
         DO I = MP1,N,4
            DY(I) = DY(I) + DA*DX(I)
            DY(I+1) = DY(I+1) + DA*DX(I+1)
            DY(I+2) = DY(I+2) + DA*DX(I+2)
            DY(I+3) = DY(I+3) + DA*DX(I+3)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
         END DO
      END IF
      RETURN
      END
      DOUBLE PRECISION FUNCTION DCABS1(Z)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX Z
*     ..
*     ..
*  Purpose
*  =======
*
*  DCABS1 computes absolute value of a double complex number 
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,DIMAG
*
      DCABS1 = ABS(DBLE(Z)) + ABS(DIMAG(Z))
      RETURN
      END
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     DCOPY copies a vector, x, to a vector, y.
*     uses unrolled loops for increments equal to one.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,7)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DX(I)
            END DO
            IF (N.LT.7) RETURN
         END IF   
         MP1 = M + 1
         DO I = MP1,N,7
            DY(I) = DX(I)
            DY(I+1) = DX(I+1)
            DY(I+2) = DX(I+2)
            DY(I+3) = DX(I+3)
            DY(I+4) = DX(I+4)
            DY(I+5) = DX(I+5)
            DY(I+6) = DX(I+6)
         END DO
      ELSE      
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DY(IY) = DX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     DDOT forms the dot product of two vectors.
*     uses unrolled loops for increments equal to one.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      DDOT = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DTEMP + DX(I)*DY(I)
            END DO
            IF (N.LT.5) THEN
               DDOT=DTEMP
            RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     $            DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DTEMP + DX(IX)*DY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      DDOT = DTEMP
      RETURN
      END
      DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*)
*     ..
*
*  Purpose
*  =======
*
*  DNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DNRM2 := sqrt( x'*x )
*
*  Further Details
*  ===============
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to DLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
      INTEGER IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
*     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE IF (N.EQ.1) THEN
          NORM = ABS(X(1))
      ELSE
          SCALE = ZERO
          SSQ = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
*
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (X(IX).NE.ZERO) THEN
                  ABSXI = ABS(X(IX))
                  IF (SCALE.LT.ABSXI) THEN
                      SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                      SCALE = ABSXI
                  ELSE
                      SSQ = SSQ + (ABSXI/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
*
      DNRM2 = NORM
      RETURN
*
*     End of DNRM2.
*
      END
      SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
*     .. Scalar Arguments ..
      DOUBLE PRECISION C,S
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     DROT applies a plane rotation.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*       code for both increments equal to 1
*
         DO I = 1,N
            DTEMP = C*DX(I) + S*DY(I)
            DY(I) = C*DY(I) - S*DX(I)
            DX(I) = DTEMP
         END DO
      ELSE
*
*       code for unequal increments or equal increments not equal
*         to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = C*DX(IX) + S*DY(IY)
            DY(IY) = C*DY(IY) - S*DX(IX)
            DX(IX) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
      SUBROUTINE DROTG(DA,DB,C,S)
*     .. Scalar Arguments ..
      DOUBLE PRECISION C,DA,DB,S
*     ..
*
*  Purpose
*  =======
*
*     DROTG construct givens plane rotation.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION R,ROE,SCALE,Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DABS,DSIGN,DSQRT
*     ..
      ROE = DB
      IF (DABS(DA).GT.DABS(DB)) ROE = DA
      SCALE = DABS(DA) + DABS(DB)
      IF (SCALE.EQ.0.0d0) THEN
         C = 1.0d0
         S = 0.0d0
         R = 0.0d0
         Z = 0.0d0
      ELSE
         R = SCALE*DSQRT((DA/SCALE)**2+ (DB/SCALE)**2)
         R = DSIGN(1.0d0,ROE)*R
         C = DA/R
         S = DB/R
         Z = 1.0d0
         IF (DABS(DA).GT.DABS(DB)) Z = S
         IF (DABS(DB).GE.DABS(DA) .AND. C.NE.0.0d0) Z = 1.0d0/C
      END IF
      DA = R
      DB = Z
      RETURN
      END
      SUBROUTINE DROTM(N,DX,INCX,DY,INCY,DPARAM)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DPARAM(5),DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
*
*     (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
*     (DY**T)
*
*     DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
*     LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
*     WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
*
*     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
*
*       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
*     H=(          )    (          )    (          )    (          )
*       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
*     SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         number of elements in input vector(s)
*
*  DX     (input/output) DOUBLE PRECISION array, dimension N
*         double precision vector with N elements
*
*  INCX   (input) INTEGER
*         storage spacing between elements of DX
*
*  DY     (input/output) DOUBLE PRECISION array, dimension N
*         double precision vector with N elements
*
*  INCY   (input) INTEGER
*         storage spacing between elements of DY
*
*  DPARAM (input/output)  DOUBLE PRECISION array, dimension 5 
*     DPARAM(1)=DFLAG
*     DPARAM(2)=DH11
*     DPARAM(3)=DH21
*     DPARAM(4)=DH12
*     DPARAM(5)=DH22
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DFLAG,DH11,DH12,DH21,DH22,TWO,W,Z,ZERO
      INTEGER I,KX,KY,NSTEPS
*     ..
*     .. Data statements ..
      DATA ZERO,TWO/0.D0,2.D0/
*     ..
*
      DFLAG = DPARAM(1)
      IF (N.LE.0 .OR. (DFLAG+TWO.EQ.ZERO)) RETURN
      IF (INCX.EQ.INCY.AND.INCX.GT.0) THEN
*
         NSTEPS = N*INCX
         IF (DFLAG.LT.ZERO) THEN
            DH11 = DPARAM(2)
            DH12 = DPARAM(4)
            DH21 = DPARAM(3)
            DH22 = DPARAM(5)
            DO I = 1,NSTEPS,INCX
               W = DX(I)
               Z = DY(I)
               DX(I) = W*DH11 + Z*DH12
               DY(I) = W*DH21 + Z*DH22
            END DO
         ELSE IF (DFLAG.EQ.ZERO) THEN
            DH12 = DPARAM(4)
            DH21 = DPARAM(3)
            DO I = 1,NSTEPS,INCX
               W = DX(I)
               Z = DY(I)
               DX(I) = W + Z*DH12
               DY(I) = W*DH21 + Z
            END DO
         ELSE
            DH11 = DPARAM(2)
            DH22 = DPARAM(5)
            DO I = 1,NSTEPS,INCX
               W = DX(I)
               Z = DY(I)
               DX(I) = W*DH11 + Z
               DY(I) = -W + DH22*Z
            END DO
         END IF
      ELSE
         KX = 1
         KY = 1
         IF (INCX.LT.0) KX = 1 + (1-N)*INCX
         IF (INCY.LT.0) KY = 1 + (1-N)*INCY
*
         IF (DFLAG.LT.ZERO) THEN
            DH11 = DPARAM(2)
            DH12 = DPARAM(4)
            DH21 = DPARAM(3)
            DH22 = DPARAM(5)
            DO I = 1,N
               W = DX(KX)
               Z = DY(KY)
               DX(KX) = W*DH11 + Z*DH12
               DY(KY) = W*DH21 + Z*DH22
               KX = KX + INCX
               KY = KY + INCY
            END DO
         ELSE IF (DFLAG.EQ.ZERO) THEN
            DH12 = DPARAM(4)
            DH21 = DPARAM(3)
            DO I = 1,N
               W = DX(KX)
               Z = DY(KY)
               DX(KX) = W + Z*DH12
               DY(KY) = W*DH21 + Z
               KX = KX + INCX
               KY = KY + INCY
            END DO
         ELSE
             DH11 = DPARAM(2)
             DH22 = DPARAM(5)
             DO I = 1,N
                W = DX(KX)
                Z = DY(KY)
                DX(KX) = W*DH11 + Z
                DY(KY) = -W + DH22*Z
                KX = KX + INCX
                KY = KY + INCY
            END DO
         END IF
      END IF
      RETURN
      END
      SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DD1,DD2,DX1,DY1
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DPARAM(5)
*     ..
*
*  Purpose
*  =======
*
*     CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
*     THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)*DX1,DSQRT(DD2)*
*     DY2)**T.
*     WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
*
*     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
*
*       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
*     H=(          )    (          )    (          )    (          )
*       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
*     LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
*     RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE
*     VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
*
*     THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
*     INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
*     OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
*
*
*  Arguments
*  =========
*
*  DD1    (input/output) DOUBLE PRECISION
*
*  DD2    (input/output) DOUBLE PRECISION
*
*  DX1    (input/output) DOUBLE PRECISION
*
*  DY1    (input) DOUBLE PRECISION
*
*  DPARAM (input/output)  DOUBLE PRECISION array, dimension 5
*     DPARAM(1)=DFLAG
*     DPARAM(2)=DH11
*     DPARAM(3)=DH21
*     DPARAM(4)=DH12
*     DPARAM(5)=DH22
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DFLAG,DH11,DH12,DH21,DH22,DP1,DP2,DQ1,DQ2,DTEMP,
     $                 DU,GAM,GAMSQ,ONE,RGAMSQ,TWO,ZERO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DABS
*     ..
*     .. Data statements ..
*
      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/
      DATA GAM,GAMSQ,RGAMSQ/4096.D0,16777216.D0,5.9604645D-8/
*     ..

      IF (DD1.LT.ZERO) THEN
*        GO ZERO-H-D-AND-DX1..
         DFLAG = -ONE
         DH11 = ZERO
         DH12 = ZERO
         DH21 = ZERO
         DH22 = ZERO
*
         DD1 = ZERO
         DD2 = ZERO
         DX1 = ZERO
      ELSE
*        CASE-DD1-NONNEGATIVE
         DP2 = DD2*DY1
         IF (DP2.EQ.ZERO) THEN
            DFLAG = -TWO
            DPARAM(1) = DFLAG
            RETURN
         END IF 
*        REGULAR-CASE..
         DP1 = DD1*DX1
         DQ2 = DP2*DY1
         DQ1 = DP1*DX1
*
         IF (DABS(DQ1).GT.DABS(DQ2)) THEN
            DH21 = -DY1/DX1
            DH12 = DP2/DP1
*
            DU = ONE - DH12*DH21
*
           IF (DU.GT.ZERO) THEN
             DFLAG = ZERO
             DD1 = DD1/DU
             DD2 = DD2/DU
             DX1 = DX1*DU
           END IF
         ELSE

            IF (DQ2.LT.ZERO) THEN
*              GO ZERO-H-D-AND-DX1..
               DFLAG = -ONE
               DH11 = ZERO
               DH12 = ZERO
               DH21 = ZERO
               DH22 = ZERO
*
               DD1 = ZERO
               DD2 = ZERO
               DX1 = ZERO
            ELSE
               DFLAG = ONE
               DH11 = DP1/DP2
               DH22 = DX1/DY1
               DU = ONE + DH11*DH22
               DTEMP = DD2/DU
               DD2 = DD1/DU
               DD1 = DTEMP
               DX1 = DY1*DU
            END IF
         END IF

*     PROCEDURE..SCALE-CHECK
         IF (DD1.NE.ZERO) THEN
            DO WHILE ((DD1.LE.RGAMSQ) .OR. (DD1.GE.GAMSQ))
               IF (DFLAG.EQ.ZERO) THEN
                  DH11 = ONE
                  DH22 = ONE
                  DFLAG = -ONE
               ELSE
                  DH21 = -ONE
                  DH12 = ONE
                  DFLAG = -ONE
               END IF
               IF (DD1.LE.RGAMSQ) THEN
                  DD1 = DD1*GAM**2
                  DX1 = DX1/GAM
                  DH11 = DH11/GAM
                  DH12 = DH12/GAM
               ELSE
                  DD1 = DD1/GAM**2
                  DX1 = DX1*GAM
                  DH11 = DH11*GAM
                  DH12 = DH12*GAM
               END IF
            ENDDO
         END IF
  
         IF (DD2.NE.ZERO) THEN
            DO WHILE ( (DABS(DD2).LE.RGAMSQ) .OR. (DABS(DD2).GE.GAMSQ) )
               IF (DFLAG.EQ.ZERO) THEN
                  DH11 = ONE
                  DH22 = ONE
                  DFLAG = -ONE
               ELSE
                  DH21 = -ONE
                  DH12 = ONE
                  DFLAG = -ONE
               END IF
               IF (DABS(DD2).LE.RGAMSQ) THEN
                  DD2 = DD2*GAM**2
                  DH21 = DH21/GAM
                  DH22 = DH22/GAM
               ELSE
                  DD2 = DD2/GAM**2
                  DH21 = DH21*GAM
                  DH22 = DH22*GAM
               END IF      
            END DO
         END IF
     
      END IF

      IF (DFLAG.LT.ZERO) THEN
         DPARAM(2) = DH11
         DPARAM(3) = DH21
         DPARAM(4) = DH12
         DPARAM(5) = DH22
      ELSE IF (DFLAG.EQ.ZERO) THEN
         DPARAM(3) = DH21
         DPARAM(4) = DH12 
      ELSE
         DPARAM(2) = DH11
         DPARAM(5) = DH22
      END IF

  260 CONTINUE
      DPARAM(1) = DFLAG
      RETURN
      END
      
     
     
     
      SUBROUTINE DSCAL(N,DA,DX,INCX)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  Purpose
*  =======
*
*     DSCAL scales a vector by a constant.
*     uses unrolled loops for increment equal to one.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DX(I) = DA*DX(I)
            END DO
            IF (N.LT.5) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            DX(I) = DA*DX(I)
            DX(I+1) = DA*DX(I+1)
            DX(I+2) = DA*DX(I+2)
            DX(I+3) = DA*DX(I+3)
            DX(I+4) = DA*DX(I+4)
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DX(I) = DA*DX(I)
         END DO
      END IF
      RETURN
      END
      DOUBLE PRECISION FUNCTION DSDOT(N,SX,INCX,SY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      REAL SX(*),SY(*)
*     ..
*
*  AUTHORS
*  =======
*  Lawson, C. L., (JPL), Hanson, R. J., (SNLA), 
*  Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)
*
*  Purpose
*  =======
*  Compute the inner product of two vectors with extended
*  precision accumulation and result.
*
*  Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
*  DSDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
*  where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
*  defined in a similar way using INCY.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         number of elements in input vector(s)
*
*  SX     (input) REAL array, dimension(N)
*         single precision vector with N elements
*
*  INCX   (input) INTEGER
*          storage spacing between elements of SX
*
*  SY     (input) REAL array, dimension(N)
*         single precision vector with N elements
*
*  INCY   (input) INTEGER
*         storage spacing between elements of SY
*
*  DSDOT  (output) DOUBLE PRECISION
*         DSDOT  double precision dot product (zero if N.LE.0)
*
*  Further Details
*  ===============
*
*  REFERENCES
*      
*  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
*  Krogh, Basic linear algebra subprograms for Fortran
*  usage, Algorithm No. 539, Transactions on Mathematical
*  Software 5, 3 (September 1979), pp. 308-323.
*
*  REVISION HISTORY  (YYMMDD)
*
*  791001  DATE WRITTEN
*  890831  Modified array declarations.  (WRB)
*  890831  REVISION DATE from Version 3.2
*  891214  Prologue converted to Version 4.0 format.  (BAB)
*  920310  Corrected definition of LX in DESCRIPTION.  (WRB)
*  920501  Reformatted the REFERENCES section.  (WRB)
*  070118  Reformat to LAPACK style (JL)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,KX,KY,NS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DBLE
*     ..
      DSDOT = 0.0D0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.INCY .AND. INCX.GT.0) THEN
*
*     Code for equal, positive, non-unit increments.
*
         NS = N*INCX
         DO I = 1,NS,INCX
            DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I))
         END DO
      ELSE
*
*     Code for unequal or nonpositive increments.
*
         KX = 1
         KY = 1
         IF (INCX.LT.0) KX = 1 + (1-N)*INCX
         IF (INCY.LT.0) KY = 1 + (1-N)*INCY
         DO I = 1,N
            DSDOT = DSDOT + DBLE(SX(KX))*DBLE(SY(KY))
            KX = KX + INCX
            KY = KY + INCY
         END DO
      END IF
      RETURN
      END
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     interchanges two vectors.
*     uses unrolled loops for increments equal one.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*       code for both increments equal to 1
*
*
*       clean-up loop
*
         M = MOD(N,3)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DX(I)
               DX(I) = DY(I)
               DY(I) = DTEMP
            END DO
            IF (N.LT.3) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,3
            DTEMP = DX(I)
            DX(I) = DY(I)
            DY(I) = DTEMP
            DTEMP = DX(I+1)
            DX(I+1) = DY(I+1)
            DY(I+1) = DTEMP
            DTEMP = DX(I+2)
            DX(I+2) = DY(I+2)
            DY(I+2) = DTEMP
         END DO
      ELSE
*
*       code for unequal increments or equal increments not equal
*         to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DX(IX)
            DX(IX) = DY(IY)
            DY(IY) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
      DOUBLE PRECISION FUNCTION DZASUM(N,ZX,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX ZX(*)
*     ..
*
*  Purpose
*  =======
*
*     DZASUM takes the sum of the absolute values.
*
*  Further Details
*  ===============
*
*     jack dongarra, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION STEMP
      INTEGER I,NINCX
*     ..
*     .. External Functions ..
      DOUBLE PRECISION DCABS1
      EXTERNAL DCABS1
*     ..
      DZASUM = 0.0d0
      STEMP = 0.0d0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
         DO I = 1,N
            STEMP = STEMP + DCABS1(ZX(I))
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            STEMP = STEMP + DCABS1(ZX(I))
         END DO
      END IF
      DZASUM = STEMP
      RETURN
      END
      DOUBLE PRECISION FUNCTION DZNRM2(N,X,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX X(*)
*     ..
*
*  Purpose
*  =======
*
*  DZNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DZNRM2 := sqrt( x**H*x )
*
*  Further Details
*  ===============
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to ZLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION NORM,SCALE,SSQ,TEMP
      INTEGER IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,DIMAG,SQRT
*     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE
          SCALE = ZERO
          SSQ = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL ZLASSQ( N, X, INCX, SCALE, SSQ )
*
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (DBLE(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(DBLE(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
              IF (DIMAG(X(IX)).NE.ZERO) THEN
                  TEMP = ABS(DIMAG(X(IX)))
                  IF (SCALE.LT.TEMP) THEN
                      SSQ = ONE + SSQ* (SCALE/TEMP)**2
                      SCALE = TEMP
                  ELSE
                      SSQ = SSQ + (TEMP/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
*
      DZNRM2 = NORM
      RETURN
*
*     End of DZNRM2.
*
      END
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  Purpose
*  =======
*
*     IDAMAX finds the index of element having max. absolute value.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DABS
*     ..
      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
         DMAX = DABS(DX(1))
         DO I = 2,N
            IF (DABS(DX(I)).GT.DMAX) THEN
               IDAMAX = I
               DMAX = DABS(DX(I))
            END IF
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         IX = 1
         DMAX = DABS(DX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (DABS(DX(IX)).GT.DMAX) THEN
               IDAMAX = I
               DMAX = DABS(DX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
      END
