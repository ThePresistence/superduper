      SUBROUTINE DFPUPD (HS,JHS,IHS,GS,JGS,IGS,Q1,JQ1,IQ1,Q2,JQ2,IQ2,
     +                   Q9,JQ9,IQ9,N,FILT,FILTALL,CUTM,ICG,NPRINT)
C     *
C *** DAVIDON-FLETCHER-POWELL UPDATE OF SEARCH DIRECTION (H).
C     PRELIMINARY VERSION.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     HS(*)     NEW SEARCH DIRECTION MATRIX (O).
C     GS(*)     GRADIENT FOR OLD/NEW DENSITY MATRIX (I,O).
C     Q1(*)     OLD/NEW INVERSE HESSIAN (I,O).
C     Q2(*)     DIFFERENCE BETWEEN OLD AND NEW DENSITY MATRIX (I,S).
C     Q9(*)     GRADIENT FOR NEW DENSITY MATRIX (I,S).
C     IX(N+1)   POINTERS FOR SPARSE MATRICES IN CSR FORMAT (S).
C               X=GS,HS,Q1,Q8,Q9.
C     JX(*)     COLUMN INDICES FOR SPARSE MATRICES IN CSR FORMAT (S).
C               X=GS,HS,Q1,Q8,Q9.
C     N         NUMBER OF ORBITALS (I).
C     FILT      FLAG FOR REMOVING SMALL ELEMENTS IN PRODUCT A*B (I).
C     FILTALL   FLAG FOR REMOVING SMALL ELEMENTS IN SUM A+B (I).
C     CUTM      CUTOFF FOR SMALL MATRIX ELEMENTS (I).
C     ICG       NUMBER OF CURRENT CG CYCLE (I).
C     NPRINT    PRINTING FLAG (I).
C     *
      USE module3
      IMPLICIT NONE
      LOGICAL :: FILT,FILTALL
      COMMON 
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
C    ./NBFILE/ NBF(20)
      INTEGER :: ie,N,N4,NPRINT,ICG,N2,I
C     INTEGER :: NBF,NB6
      DOUBLE PRECISION :: ZERO,ONE,TWO,THREE,FOUR,PT5,PT25,AZN,BZN,CUTM
      INTEGER, ALLOCATABLE ::IW(:)
      DOUBLE PRECISION, DIMENSION (:), POINTER :: GS,HS,Q1,Q2,Q9,Q8,WR
      INTEGER, DIMENSION (:), POINTER :: JGS,IGS,JHS,IHS,JQ1,IQ1,
     +                                   JQ2,IQ2,JQ9,IQ9,JQ8,IQ8
C
C *** DFP: FIRST COMPONENT OF UPDATE TO INVERSE HESSIAN.
      ALLOCATE (IW(N),WR(N),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','IW',N,0)
C     SYMBOLIC ADDITION TO GET NUMBER OF ELEMENTS (N4).
      CALL aplbdgp (N,JGS,IGS,JQ9,IQ9,N4,IW)
C     COMPUTE DIFFERENCE (HS) OF NEW AND OLD GRADIENT MATRIX.
      ALLOCATE (HS(N4),JHS(N4),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','HS',N4,0)
      N2=IGS(N+1)-1
      DO 10 I=1,N2
      GS(I) = -GS(I)
   10 CONTINUE
      CALL aplbp (N,Q9,JQ9,IQ9,GS,JGS,IGS,HS,JHS,IHS,N4,IW,ie)
         IF(ie.NE.0) CALL XERSPA (ie,'aplb','DFPDMS',1)
C     REMOVE SMALL ELEMENTS (HS).
      IF(FILTALL) THEN
         CALL filterp (N,1,CUTM,HS,JHS,IHS,HS,JHS,IHS,N4,ie)
           IF (ie.ne.0) CALL XERSPA (ie,'filter','DFPDMS',1)
      ENDIF
C     COPY NEW GRADIENT (Q9->GS).
      N4=IQ9(N+1)-1
      IF (IGS(N+1)-1.LT.N4) THEN
         DEALLOCATE (GS,JGS,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','GS',N4,1)
         NULLIFY (GS,JGS)
         ALLOCATE (GS(N4),JGS(N4),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','GS',N4,0)
      ENDIF
      CALL copmatp (N,Q9,JQ9,IQ9,GS,JGS,IGS,1)
C
C     TRACE (AZN) OF PRODUCT (Q2*HS) OF DIFFERENCE MATRICES.
C     Q2: DIFFERENCE OF NEW AND OLD DENSITY MATRIX.
C     HS: DIFFERENCE OF NEW AND OLD GRADIENT MATRIX.
      CALL TRACENP (N,Q2,JQ2,IQ2,HS,JHS,IHS,WR,AZN)
C     SYMBOLIC MULTIPLICATION (Q2*Q2) TO GET NUMBER OF ELEMENTS (N4).
      CALL amubdgp (N,JQ2,IQ2,JQ2,IQ2,N4,IW)
      IF (IQ9(N+1)-1.LT.N4) THEN
         DEALLOCATE (Q9,JQ9,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','Q9',N4,1)
         NULLIFY (Q9,JQ9)
         ALLOCATE (Q9(N4),JQ9(N4),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','Q9',N4,0)
      ENDIF
C     ACTUAL MULTIPLICATION (Q9=Q2*Q2).
      CALL amubp (N,Q2,JQ2,IQ2,Q2,JQ2,IQ2,Q9,JQ9,IQ9,N4,IW,ie)
         IF(ie.NE.0) CALL XERSPA (ie,'amub','DFPDMS',1)
C     DEFINE FIRST COMPONENT OF UPDATE (Q9) TO INVERSE HESSIAN.
      DO 20 I=1,N4
      Q9(I)=Q9(I)/AZN
   20 CONTINUE
C     REMOVE SMALL ELEMENTS (Q9).
      IF(FILT) THEN
         CALL filterp (N,1,CUTM,Q9,JQ9,IQ9,Q9,JQ9,IQ9,N4,ie)
            IF(ie.ne.0) CALL XERSPA (ie,'filter','DFPDMS',1)
      ENDIF
      ALLOCATE (IQ8(N+1),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','IQ8',N,0)
C
C *** DFP: SECOND COMPONENT OF UPDATE TO INVERSE HESSIAN.
C
C     THE ITERATIONS REQUIRE THE MATRIX PRODUCT OF PREVIOUS
C     INVERSE HESSIAN (Q1) AND GRADIENT DIFFERENCE MATRIX (HS).
C     IN THE FIRST ITERATION, THE INVERSE HESSIAN IS A UNIT MATRIX.
C     THEREFORE SIMPLY COPY GRADIENT DIFFERENCE MATRIX (HS->Q8).
      IF(ICG.EQ.1) THEN
         N4=IHS(N+1)-1
         ALLOCATE (Q8(N4),JQ8(N4),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','Q8',N4,0)
         CALL copmatp (N,HS,JHS,IHS,Q8,JQ8,IQ8,1)
         BZN=DOT_PRODUCT(Q8,Q8)
      ELSE
C     LATER CG ITERATIONS REQUIRE FULL MATRIX MULTIPLICATION.
C     SYMBOLIC MULTIPLICATION TO GET NUMBER OF ELEMENTS (N4).
         CALL amubdgp (N,JQ1,IQ1,JHS,IHS,N4,IW)
C        ACTUAL MULTIPLICATION (Q8=Q1*HS).
         ALLOCATE (Q8(N4),JQ8(N4),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','Q8',N4,0)
         CALL amubp (N,Q1,JQ1,IQ1,HS,JHS,IHS,Q8,JQ8,IQ8,N4,IW,ie)
            IF(ie.NE.0) CALL XERSPA (ie,'amub','DFPDMS',1)
C        REMOVE SMALL ELEMENTS (Q8).
         IF(FILT) THEN
            CALL filterp (N,1,CUTM,Q8,JQ8,IQ8,Q8,JQ8,IQ8,N4,ie)
               IF(ie.ne.0) CALL XERSPA (ie,'filter','DFPDMS',1)
         ENDIF
C        TRACE (BZN) OF PRODUCT (Q8*HS).
         CALL TRACENP (N,Q8,JQ8,IQ8,HS,JHS,IHS,WR,BZN)
      ENDIF
C     NEXT OPERATION: MULTIPLY PRODUCT MATRIX (Q8) AND TRANSPOSE
C     OF GRADIENT DIFFERENCE MATRIX (HS). THE CODE MAKES USE OF
C     THE SYMMETRY OF HS AND MULTIPLIES BY HS (NOT TRANSPOSED).
C     SYMBOLIC MULTIPLICATION TO GET NUMBER OF ELEMENTS (N4).
      CALL amubdgp (N,JQ8,IQ8,JHS,IHS,N4,IW)
      IF (IQ2(N+1)-1.LT.N4) THEN
         DEALLOCATE (Q2,JQ2,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','Q2',N4,1)
         NULLIFY (Q2,JQ2)
         ALLOCATE (Q2(N4),JQ2(N4),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','Q2',N4,0)
      ENDIF
C     ACTUAL MULTIPLICATION (Q2=Q8*HS).
      CALL amubp (N,Q8,JQ8,IQ8,HS,JHS,IHS,Q2,JQ2,IQ2,N4,IW,ie)
         IF(ie.NE.0) CALL XERSPA (ie,'amub','DFPDMS',1)
C     REMOVE SMALL ELEMENTS (Q2).
      IF(FILT) THEN
         CALL filterp (N,1,CUTM,Q2,JQ2,IQ2,Q2,JQ2,IQ2,N4,ie)
           IF (ie.ne.0) CALL XERSPA (ie,'filter','DFPDMS',1)
      ENDIF
      DEALLOCATE (HS,JHS,Q8,JQ8,IQ8,WR,STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','HS',N4,1)
      NULLIFY (HS,JHS,Q8,JQ8,IQ8,WR)
C     NEXT OPERATION: MULTIPLY PRODUCT MATRIX (Q2) BY CURRECT
C     INVERSE HESSIAN (Q1) FROM THE RIGHT.
C     FIRST ITERATION: UNIT MATRIX, HENCE COPY (Q2->HS).
      IF(ICG.EQ.1) THEN
         N4=IQ2(N+1)-1
         ALLOCATE (HS(N4),JHS(N4),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','HS',N4,0)
         CALL copmatp (N,Q2,JQ2,IQ2,HS,JHS,IHS,1)
      ELSE
C     LATER CG ITERATIONS REQUIRE FULL MATRIX MULTIPLICATION.
C     SYMBOLIC MULTIPLICATION TO GET NUMBER OF ELEMENTS (N4).
         CALL amubdgp (N,JQ2,IQ2,JQ1,IQ1,N4,IW)
C        ACTUAL MULTIPLICATION (HS=Q2*Q1).
         ALLOCATE (HS(N4),JHS(N4),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','HS',N4,0)
         CALL amubp (N,Q2,JQ2,IQ2,Q1,JQ1,IQ1,HS,JHS,IHS,N4,IW,ie)
            IF(ie.NE.0) CALL XERSPA (ie,'amub','DFPDMS',1)
      ENDIF
C     DEFINE SECOND COMPONENT OF UPDATE (HS) TO INVERSE HESSIAN.
      DO 30 I=1,N4
      HS(I)=-HS(I)/BZN
   30 CONTINUE
C     REMOVE SMALL ELEMENTS (HS).
      IF(FILT) THEN
         CALL filterp (N,1,CUTM,HS,JHS,IHS,HS,JHS,IHS,N4,ie)
            IF(ie.ne.0) CALL XERSPA (ie,'filter','DFPDMS',1)
      ENDIF
C
C *** DFP: UPDATE THE INVERSE HESSIAN (Q1+Q9+HS).
C
C     SYMBOLIC ADDITION TO GET NUMBER OF ELEMENTS (N4).
      CALL aplbdgp (N,JQ9,IQ9,JHS,IHS,N4,IW)
      IF (IQ2(N+1)-1.LT.N4) THEN
         DEALLOCATE (Q2,JQ2,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','Q2',N4,1)
         NULLIFY (Q2,JQ2)
         ALLOCATE (Q2(N4),JQ2(N4),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','Q2',N4,0)
      ENDIF
C     ADD TWO COMPONENTS OF UPDATE (Q2=Q9+HS).
      CALL aplbp (N,Q9,JQ9,IQ9,HS,JHS,IHS,Q2,JQ2,IQ2,N4,IW,ie)
         IF(ie.NE.0) CALL XERSPA (ie,'aplb','DFPDMS',1)
      DEALLOCATE (Q9,JQ9,IQ9,STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','Q9',N4,1)
      NULLIFY (Q9,JQ9,IQ9)
C     REMOVE SMALL ELEMENTS (Q2).
      IF(FILTALL) THEN
         CALL filterp (N,1,CUTM,Q2,JQ2,IQ2,Q2,JQ2,IQ2,N4,ie)
            IF(ie.ne.0) CALL XERSPA (ie,'filter','DFPDMS',1)
      ENDIF
C     COPY PREVIOUS INVERSE HESSIAN (Q1->HS).
      N4=IQ1(N+1)-1
      IF (IHS(N+1)-1.LT.N4) THEN
         DEALLOCATE (HS,JHS,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','HS',N4,1)
         NULLIFY (HS,JHS)
         ALLOCATE (HS(N4),JHS(N4),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','HS',N4,0)
      ENDIF
      CALL copmatp (N,Q1,JQ1,IQ1,HS,JHS,IHS,1)
C     COMPLETE THE UPDATE OF THE INVERSE HESSIAN.
C     SYMBOLIC ADDITION TO GET NUMBER OF ELEMENTS (N4).
      CALL aplbdgp (N,JQ2,IQ2,JHS,IHS,N4,IW)
      IF (IQ1(N+1)-1.LT.N4) THEN
         DEALLOCATE (Q1,JQ1,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','Q1',N4,1)
         NULLIFY (Q1,JQ1)
         ALLOCATE (Q1(N4),JQ1(N4),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','Q1',N4,0)
      ENDIF
C     PERFORM THE FINAL UPDATE (Q1=Q2+HS).
      CALL aplbp (N,Q2,JQ2,IQ2,HS,JHS,IHS,Q1,JQ1,IQ1,N4,IW,ie)
         IF(ie.NE.0) CALL XERSPA (ie,'aplb','DFPDMS',1)
      DEALLOCATE (Q2,JQ2,IQ2,STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','Q2',N4,1)
      NULLIFY (Q2,JQ2,IQ2)
C     REMOVE SMALL ELEMENTS (Q1).
      IF(FILTALL) THEN
         CALL filterp (N,1,CUTM,Q1,JQ1,IQ1,Q1,JQ1,IQ1,N4,ie)
            IF(ie.ne.0) CALL XERSPA (ie,'filter','DFPDMS',1)
      ENDIF
C
C *** DFP: DEFINE NEW SEARCH DIRECTION (HS=-Q1*GS).
C
C     SYMBOLIC MULTIPLICATION TO GET NUMBER OF ELEMENTS (N4).
      CALL amubdgp (N,JQ1,IQ1,JGS,IGS,N4,IW)
      IF (IHS(N+1)-1.LT.N4) THEN
         DEALLOCATE (HS,JHS,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','HS',N4,1)
         NULLIFY (HS,JHS)
         ALLOCATE (HS(N4),JHS(N4),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','HS',N4,0)
      ENDIF
C     ACTUAL MULTIPLICATION OF UPDATED INVERSE HESSIAN (Q1) AND
C     GRADIENT MATRIX (GS).
      CALL amubp (N,Q1,JQ1,IQ1,GS,JGS,IGS,HS,JHS,IHS,N4,IW,ie)
         IF(ie.NE.0) CALL XERSPA (ie,'amub','DFPDMS',1)
C     REMOVE SMALL ELEMENTS (HS).
      IF(FILT) THEN
         CALL filterp (N,1,CUTM,HS,JHS,IHS,HS,JHS,IHS,N4,ie)
            IF(ie.ne.0) CALL XERSPA (ie,'filter','DFPDMS',1)
      ENDIF
C     INCLUDE NEGATIVE SIGN.
      N4=IHS(N+1)-1
      DO 40 I=1,N4
      HS(I)=-HS(I)
   40 CONTINUE
      DEALLOCATE (IW,STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'DFPDMS','IW',N4,1)
      RETURN
      END
