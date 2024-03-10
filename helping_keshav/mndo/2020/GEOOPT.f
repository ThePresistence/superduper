      SUBROUTINE GEOOPT (A,LM5,ICALL,SCFCAL)
C     *
C     DRIVER ROUTINE FOR GEOMETRY OPTIMIZATIONS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     A(LM5)    AVAILABLE BUFFER (S).
C     ICALL     CONTROL AND ERROR FLAG (O,S).
C     SCFCAL    EXTERNAL ROUTINE FOR ENERGY EVALUATION (I).
C     *
C     ICALL=0   INITIAL VALUE UPON INPUT.
C     ICALL>0   INTERNAL CONTROL IN WORKING ROUTINES (S).
C     ICALL<0   ERROR FLAG (O).
C     ICALL=-1  NO SCF CONVERGENCE.
C     ICALL=-3  NO CONVERGENCE IN GEOMETRY OPTIMIZATION.
C     ICALL=-8  TIME LIMIT.
C     ICALL=10  SUCCESSFUL OPTIMIZATION.
C     *
      USE LIMIT, ONLY: LM1, LMV, LMYL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL SCFCAL
      COMMON
     ./DFP   / X(LMV),NVAR
     ./INOPT2/ IN2(300)
     ./LMFOR / LF1,LF2,LF3,LF4,LF5,LF6,LF7,LF8,LFA
     ./NBFILE/ NBF(20)
     ./OPTCRT/ NSUCC
     ./OVERLY/ IOV,JOV,KOV,LOV
     ./PARM1 / AA(LMV),NA(LM1,4),NATOMS
     ./STATCG/ NCGTOT(20)
     ./YARLAG/ YLVAL(LMYL),YLGRD(LMYL),YLD(LMYL),NYL
      DIMENSION A(LM5)
C *** FILE NUMBERS.
      NB4    = NBF(4)
      NB6    = NBF(6)
      NB7    = NBF(7)
      NB9    = NBF(9)
C *** INITIALIZATION.
      IOP    = IN2(2)
      JOP    = IN2(3)
      IEF    = IN2(8)
      NSAV7  = IN2(14)
      NSAV9  = IN2(16)
      IPRINT = IN2(38)
      ICROSS = IN2(160)
      KOV    = 0
      LOV    = 0
C *** CHECK FOR THE SPECIAL CASE OF A RESTART.
      IF(JOP.GT.2 .AND. IN2(37).GT.0) THEN
         REWIND(NB4)
         READ(NB4) IDREAD
         IF(IDREAD.EQ.2) RETURN
      ENDIF
C *** CALL WORKING ROUTINES.
C     HDLC OPTIMIZER FOR AN ENERGY MINIMUM OR TRANSITION STATE.
      IF(JOP.GE.0 .AND. JOP.LE.4 .AND. JOP.NE.2 .AND. IEF.LT.0) THEN
C        CALL GMETRY (+1)
         CALL HDLOPT (A,LM5,ICALL,SCFCAL)
         CALL PRTHDL
         ID4 = 4
         REWIND(NB4)
         WRITE(NB4) ID4
C     OPTIMIZATION FOR AN ENERGY MINIMUM.
      ELSE IF(JOP.EQ.0 .OR. JOP.EQ.3) THEN
         IF(IEF.EQ.0) THEN
            CALL FLEPO (A(LF1),LFA,A,LM5,ICALL,SCFCAL)
         ELSE IF(IEF.GT.0) THEN
            NIMAG = 0
            IF(ICROSS.EQ.5) THEN
               NVARLM = NVAR+NYL
            ELSE
               NVARLM = NVAR
            ENDIF
            CALL EIGF (A(LF1),A(LF2),A(LF8),NVARLM,NIMAG,
     1                 A,LM5,ICALL,SCFCAL)
         ENDIF
         IF(NSUCC.GE.10) ICALL=-3
C     OPTIMIZATION FOR A TRANSITION STATE.
      ELSE IF(JOP.EQ.1 .OR. JOP.EQ.4) THEN
         IF(IEF.EQ.0) THEN
            CALL TSOPT (A(LF1),A(LF2),A(LF8),NVAR,A,LM5,ICALL,
     1                  SCFCAL)
         ELSE IF(IEF.GT.0) THEN
            NIMAG = 1
            CALL EIGF (A(LF1),A(LF2),A(LF8),NVAR,NIMAG,
     1                 A,LM5,ICALL,SCFCAL)
         ENDIF
         IF(NSUCC.GE.10) ICALL=-3
      ELSE
         RETURN
      ENDIF
C *** COUNT UNSUCCESSFUL OPTIMIZATIONS.
      IF(ICALL.EQ.-3) NCGTOT(11) = NCGTOT(11)+1
      IF(ICALL.EQ.-1) NCGTOT(12) = NCGTOT(12)+1
C *** RESET OPTIONS (MIDDLE,IFAST) AND SAVE RESULTS.
      IF(IN2(37).GT.0) IN2(37)=0
      IF(IN2(73).LT.0) IN2(73)=0
      IF(ICALL.EQ.-1) RETURN
      IF((IOP.EQ.-3 .OR. IOP.EQ.-13) .AND. IPRINT.GT.-5) CALL HBONDT(0)
      IF(NSAV7.GT.0) CALL GEOSAV (NATOMS,NSAV7,NB7)
      IF(NSAV9.GT.0) CALL PDBSAV (NB9)
      IF(ICALL.EQ.-8) RETURN
C *** ERROR EXIT: UNSUCCESSFUL OPTIMIZATION.
      IF(ICALL.EQ.-3 .AND. IPRINT.GT.-5) THEN
         WRITE(NB6,500)
         IF(JOP.GT.2) WRITE(NB6,510)
      ENDIF
      RETURN
  500 FORMAT(///1X,'UNSUCCESSFUL GEOMETRY OPTIMIZATION.')
  510 FORMAT(/  1X,'FORCE CONSTANTS WILL NOT BE COMPUTED.')
      END