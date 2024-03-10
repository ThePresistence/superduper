      SUBROUTINE PURIFYN (PS,JPS,IPS,N,DD1,DD2,NPR,FILT,FILTALL,CUTM)
C     *
C     McWEENY PURIFICATION TRANSFORMATION: R = 3 PP - 2 PPP.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     PS(*)     INPUT DENSITY MATRIX (I), PURIFIED DENSITY MATRIX (O).
C     IX(N+1)   POINTERS FOR SPARSE MATRICES IN CSR FORMAT.
C               X=PS (I,S,O).
C     JX(*)     COLUMN INDICES FOR SPARSE MATRICES IN CSR FORMAT.
C               X=PS (I,S,O).
C     N         NUMBER OF ORBITALS (I).
C     DD1       MAXIMUM DEVIATION BETWEEN P AND PP ON DIAGONAL (O).
C     DD1       RMS DEVIATION BETWEEN P AND PP ON DIAGONAL (O).
C     NPR       PRINTING FLAG (I).
C     FILT      FLAG FOR REMOVING SMALL ELEMENTS IN PRODUCT A*B (I).
C     FILTALL   FLAG FOR REMOVING SMALL ELEMENTS IN SUM A+B (I).
C     CUTM      CUTOFF FOR SMALL MATRIX ELEMENTS (I).
C     *
      USE module3
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
      INTERFACE
      SUBROUTINE MATDEVS (A,JA,IA,B,JB,IB,N,MODE,DEVMAX,DEVRMS)
      IMPLICIT NONE
      INTEGER :: N,MODE
      DOUBLE PRECISION :: DEVMAX,DEVRMS
      DOUBLE PRECISION, DIMENSION (:), POINTER :: B,A
      INTEGER, DIMENSION (:), POINTER :: JB,IB,JA,IA
      END SUBROUTINE MATDEVS
      END INTERFACE
C
      INTEGER :: N,ie,I,N4,N5,NPR,NB6,NBF
      DOUBLE PRECISION :: TWO,THREE,DD1,DD2,CUTM
      PARAMETER (TWO=2.0D0)
      PARAMETER (THREE=3.0D0)
      INTEGER, ALLOCATABLE :: IW(:)
      DOUBLE PRECISION, DIMENSION (:), POINTER :: PS,QS,Q5
      INTEGER, DIMENSION (:), POINTER :: IPS,JPS,IQS,JQS,JQ5,IQ5
      LOGICAL :: FILT,FILTALL
C     DOUBLE PRECISION :: ZTIME,TT1A,TT1B
C     COMMON /ZTIMES/ ZTIME(20)
      COMMON /NBFILE/ NBF(20)
C
C *** PRODUCT QS = PS*PS = PP
      ALLOCATE (IW(N+1),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'PURIFYN','IW',N,0)
      CALL amubdgp (N,JPS,IPS,JPS,IPS,N4,IW)
      ALLOCATE (QS(N4),JQS(N4),IQS(N+1),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'PURIFYN','QS',N4,0)
C     CALL CPUSEC (TT1A)
      CALL amubp (N,PS,JPS,IPS,PS,JPS,IPS,QS,JQS,IQS,N4,IW,ie)
C     CALL CPUSEC (TT1B)
C     ZTIME(8) = ZTIME(8) + TT1B - TT1A
C     WRITE(6,*) ' PURIFYN: PP ',IPS(N+1)-1,IPS(N+1)-1,N4,TT1B-TT1A
         IF(ie.NE.0) CALL XERSPA (ie,'amub','PURIFYN',1)
C     REMOVE SMALL ELEMENTS.
      IF(FILT) THEN
         CALL filterp (N,1,CUTM,QS,JQS,IQS,QS,JQS,IQS,N4,ie)
           IF (ie.ne.0) CALL XERSPA (ie,'filter','PURIFYN',1)
      ENDIF
C     MAXIMUM DEVIATIONS FOR DEBUG PRINT.
      IF (NPR.GE.5) THEN
         CALL MATDEVS (PS,JPS,IPS,QS,JQS,IQS,N,1,DD1,DD2)
      ENDIF
C
C *** PRODUCT Q5 = PS*QS = PPP
      CALL amubdgp (N,JPS,IPS,JQS,IQS,N4,IW)
      ALLOCATE (Q5(N4),JQ5(N4),IQ5(N+1),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'PURIFYN','Q5',N4,0)
C     CALL CPUSEC (TT1A)
      CALL amubp (N,PS,JPS,IPS,QS,JQS,IQS,Q5,JQ5,IQ5,N4,IW,ie)
C     CALL CPUSEC (TT1B)
C     ZTIME(9) = ZTIME(9) + TT1B - TT1A
C     WRITE(6,*) ' PURIFYN: PPP',IPS(N+1)-1,IQS(N+1)-1,N4,TT1B-TT1A
         IF(ie.NE.0) CALL XERSPA (ie,'amub','PURIFYN',2)
C     REMOVE SMALL ELEMENTS.
      IF(FILT) THEN
         CALL filterp (N,1,CUTM,Q5,JQ5,IQ5,Q5,JQ5,IQ5,N4,ie)
           IF (ie.ne.0) CALL XERSPA (ie,'filter','PURIFYN',1)
      ENDIF
C
C *** PURIFIED MATRIX PS = 3*QS - 2*Q5 = 3 PP - 2 PPP
      DEALLOCATE (PS,JPS,STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'PURIFYN','PS',N4,1)
      NULLIFY (PS,JPS)
      CALL aplbdgp (N,JQS,IQS,JQ5,IQ5,N4,IW)
      DO 10 I=1,IQS(N+1)-1
      QS(I)=THREE*QS(I)
   10 CONTINUE
      DO 20 I=1,IQ5(N+1)-1
      Q5(I)=-TWO*Q5(I)
   20 CONTINUE
      N5 = IPS(N+1)-1
      ALLOCATE (PS(N4),JPS(N4),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'PURIFYN','PS',N4,0)
      CALL aplbp (N,QS,JQS,IQS,Q5,JQ5,IQ5,PS,JPS,IPS,N4,IW,ie)
         IF(ie.NE.0) CALL XERSPA (ie,'amub','PURIFYN',2)
C     DEBUG PRINT.
      IF(NPR.GE.7) THEN
         NB6 = NBF(6)
         WRITE(NB6,600) N5,IPS(N+1)-1
         WRITE(NB6,610) IQS(N+1)-1,IQ5(N+1)-1
      ENDIF
      DEALLOCATE (IW,QS,JQS,IQS,Q5,JQ5,IQ5,STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'PURIFYN','IW',N,1)
      NULLIFY (QS,JQS,IQS,Q5,JQ5,IQ5)
C     REMOVE SMALL ELEMENTS.
      IF(FILTALL) THEN
         CALL filterp (N,1,CUTM,PS,JPS,IPS,PS,JPS,IPS,N4,ie)
           IF (ie.ne.0) CALL XERSPA (ie,'filter','PURIFYN',4)
      ENDIF
      RETURN
  600 FORMAT(/  1X,'PURIFYN: NONZERO ELEMENTS IN P , R   :',2I10)
  610 FORMAT(/  1X,'PURIFYN: NONZERO ELEMENTS IN PP, PPP :',2I10)
      END
