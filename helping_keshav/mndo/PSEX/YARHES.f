      SUBROUTINE YARHES (HESS,LMH,IPRINT)
C     *
C     UPDATE OF THE HESSIAN FOR YARKONY CI SEARCH.
C     NB: THE TOP LEFT SQUARE OF THE HESSIAN IS UPDATED IN UPDHES.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     HESS      HESSIAN MATRIX FOR NEWTON-RAPHSON OPTIMIZATION (I,O)
C     LMH       ALLOCATED DIMENSIONS OF HESSIAN (I)
C     IPRINT    PRINTING FLAG (I)
C     *
      USE LIMIT, ONLY: LM1,LMV,LM1M,LMPROP,LMSTAT,LMGRD,LMNAC,LMYL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DFP   / XPARAM(LMV),NVAR
     ./NBFILE/ NBF(20)
     ./YARCON/ YLGCV(LMYL),YLGCT(LMYL),IYLGC(5,LMYL)
     ./YARLAG/ YLVAL(LMYL),YLGRD(LMYL),YLD(LMYL),NYL
     ./YARUPD/ YARV(3,LM1+LM1M,LMYL),YARS(LMV,2),YARG(LMV,LMYL,2),IXTRAP
      DIMENSION HESS(LMH,LMH)
C *** FILE NUMBERS.
      NB6    = NBF(6)
C *** INITIALIZATION.
      I3N    = 3*NUMAT
      NUMALL = NUMAT+NUMATM
      NVARLM = NVAR+NYL

C *** UPDATE HESSIAN
      DO 10 I=1,NVAR
         DO 20 J=1,NYL
            HESS(I,NVAR+J) = YARG(I,J,1)
            HESS(NVAR+J,I) = YARG(I,J,1)
 20      CONTINUE
 10   CONTINUE

C     BOTTOM RIGHT SQUARE OF HESSIAN = 0 
      DO 30 I=1,NYL
         DO 40 J=1,NYL
            HESS(NVAR+I,NVAR+J) = ZERO
 40      CONTINUE
 30   CONTINUE

      RETURN
      END
