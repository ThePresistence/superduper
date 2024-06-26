      SUBROUTINE TSSAV (Q,R,LMQ,MODE,MIDDLE)
C     *
C     TSOPT INFORMATION IS SAVED (MODE=0) OR READ (MODE=1,2).
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     Q         Q MATRIX IN NON-LINEAR LEAST-SQUARES ALGORITHM (S).
C     R         R MATRIX IN NON-LINEAR LEAST-SQUARES ALGORITHM (S).
C     LMQ       DIMENSION OF Q AND R MATRICES (I).
C     MODE      CHOICE BETWEEN WRITE/READ, SEE ABOVE (I).
C     MIDDLE    NO ACTION FOR MIDDLE.LT.0 (I).
C     *
      USE LIMIT, ONLY: LMV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ID4=1)
      COMMON
     ./CYCLES/ ICYC(2)
     ./DFP   / X(LMV),N
     ./DFPGO / GLAST(LMV),XLAST(LMV)
     ./ERG   / ENERGY,GRAD(LMV+2)
     ./FLAG2 / SECADD,TIME1
     ./FLPOCM/ G(LMV),P(LMV),ALPHA(2)
      COMMON
     ./NBFILE/ NBF(20)
     ./OPCOM1/ XD(LMV),GD(LMV),GNORM(8)
     ./OPCOM2/ IREPET(4)
     ./OPTRST/ CNCADD(7)
     ./OVERLY/ IOV,JOV,KOV,LOV
     ./PARM6 / LPT(3)
      DIMENSION Q(LMQ,LMQ),R(LMQ,LMQ)
C *** FILE NUMBERS.
      NB4    = NBF(4)
      NB6    = NBF(6)
C *** INITIALIZATION.
      IF(MIDDLE.LT.0) RETURN
      REWIND NB4
      IF(MODE.GT.0) GO TO 10
C *** MODE=0. SAVE TSOPT INFORMATION ON FILE NB4.
      CALL CPUSEC(TIME2)
      SECTOT = SECADD+TIME2-TIME1
      IF((KOV.EQ.0 .AND. LOV.EQ.1) .OR. LOV.GT.4) THEN
         WRITE(NB4) ID4
         WRITE(NB4) ((Q(J,I),J=1,N),I=1,N)
         WRITE(NB4) ((R(J,I),J=1,N),I=1,N)
      ELSE
         READ (NB4) IDREAD
         READ (NB4) ((Q(J,I),J=1,N),I=1,N)
         READ (NB4) ((R(J,I),J=1,N),I=1,N)
      ENDIF
      WRITE(NB4) ICYC,IREPET,IOV,JOV,KOV,LOV,LPT
      WRITE(NB4) SECTOT,ALPHA,GNORM,CNCADD,ENERGY
      WRITE(NB4) (X(I),I=1,N)
      WRITE(NB4) (G(I),I=1,N)
      WRITE(NB4) (P(I),I=1,N)
      WRITE(NB4) (XLAST(I),I=1,N)
      WRITE(NB4) (GLAST(I),I=1,N)
      WRITE(NB4) (XD(I),I=1,N)
      WRITE(NB4) (GD(I),I=1,N)
      RETURN
C *** MODE=1,2. READ TSOPT INFORMATION FROM FILE NB4.
   10 CONTINUE
      READ(NB4,ERR=100,END=200) IDREAD
      IF(IDREAD.NE.ID4) GO TO 300
      READ(NB4,ERR=100,END=200) ((Q(J,I),J=1,N),I=1,N)
      READ(NB4,ERR=100,END=200) ((R(J,I),J=1,N),I=1,N)
      IF(MODE.GT.1) RETURN
      READ(NB4,ERR=100,END=200) ICYC,IREPET,IOV,JOV,KOV,LOV,LPT
      READ(NB4,ERR=100,END=200) SECADD,ALPHA,GNORM,CNCADD,ENERGY
      READ(NB4,ERR=100,END=200) (X(I),I=1,N)
      READ(NB4,ERR=100,END=200) (G(I),I=1,N)
      READ(NB4,ERR=100,END=200) (P(I),I=1,N)
      READ(NB4,ERR=100,END=200) (XLAST(I),I=1,N)
      READ(NB4,ERR=100,END=200) (GLAST(I),I=1,N)
      READ(NB4,ERR=100,END=200) (XD(I),I=1,N)
      READ(NB4,ERR=100,END=200) (GD(I),I=1,N)
      RETURN
C     ERRORS EXITS.
  100 WRITE(NB6,500)
      STOP 'TSSAV'
  200 WRITE(NB6,510)
      STOP 'TSSAV'
  300 WRITE(NB6,500)
      WRITE(NB6,520) IDREAD,ID4
      STOP 'TSSAV'
  500 FORMAT(///1X,'ERROR WHEN READING RESTART FILE',
     1       /  1X,'FOR TRANSITION STATE SEARCH',
     2       /  1X,'IN SUBROUTINE TSSAV.'/)
  510 FORMAT(///1X,'END-OF-FILE ENCOUNTERED WHEN READING RESTART FILE',
     1       /  1X,'FOR TRANSITION STATE SEARCH',
     2       /  1X,'IN SUBROUTINE TSSAV.'/)
  520 FORMAT(   1X,'WRONG TYPE OF RESTART FILE: IDREAD =',I3,
     1       /  1X,'EXPECTED VALUE OF LABEL   : IDREAD =',I3,/)
      END
