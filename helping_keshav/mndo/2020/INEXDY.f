      SUBROUTINE INEXDY (IFORM,JPRINT)
C     *
C     INPUT OF SPECIFIC OPTIONS FOR SURFACE HOPPING.
C     *
      USE LIMIT, ONLY: LMGRD,LMNAC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./GRDORG/ IGRST(LMGRD),ISTATE,JSTATE
     ./INOPT2/ IN2(300)
     ./INOPT3/ XN3(50)
     ./INOPT4/ XN4(50)
C     ./NACOPT/ ICURST,NESTEP,ITULLY
     ./NBFILE/ NBF(20)
C *** FILE NUMBERS.
      NB5    = NBF(5)
      NB6    = NBF(6)
C *** INPUT VARIABLES.
      NCIGRD = IN2(159)
      ICROSS = IN2(160)
C *** READ INPUT AND INITIALIZE ARRAYS.
C *** DEFINITION OF STATES
      IF(NCIGRD.GE.1 .AND. NCIGRD.LE.LMGRD) THEN
         IF(IFORM.LE.0) THEN
            READ(NB5,500) (IGRST(I),I=1,NCIGRD)
         ELSE
            READ(NB5,*)   (IGRST(I),I=1,NCIGRD)
         ENDIF
      ELSE
         IN2(159) = 2
         IGRST(1) = 1
         IGRST(2) = 2
      ENDIF
C *** PRINTING SECTION.
      IF(JPRINT.GT.-5) THEN
         WRITE(NB6,100)
         WRITE(NB6,120) IN2(159),(IGRST(I),I=1,IN2(159))
         IF(NCIGRD.LE.0 .OR. NCIGRD.GT.LMGRD) WRITE(NB6,130)
      ENDIF
      RETURN
  100 FORMAT(///1X,'*** OPTIONS FOR SURFACE HOPPING ALGORITHM   ***',
     1       /  1X,'*** DEFINED EITHER EXPLICITLY OR BY DEFAULT ***')
  120 FORMAT(// 1X,'NUMBER OF CI STATES CONSIDERED:',I3,
     1       /  1X,'LABELS OF THESE CI STATES:     ',10I3)
  130 FORMAT(/  1X,'OPTIONS AND CI STATES SELECTED BY DEFAULT')
  500 FORMAT(20I4)
      END