      SUBROUTINE ERRSAV (NSAV15, ICALL)
C     *
C     SAVE SCF ERROR MESSAGES TO FILE NSAV15
C     TO INFORM EXTERNAL PROGRAMS OF CALCULATION FAILURE.
C     *
C     NSAV15    CONTROL OVER WHICH DATA ARE SAVED ON FILE NB15 (I).
C               = 0      DO NOT SAVE ANYTHING.
C               = 1-4,6  FORMATTED ERROR MESSAGE.
C               = 5,7    UNFORMATTED ERROR MESSAGE.
C               = 9      FORMATTED, BUT DO NOT REWIND FILE NB15.
C     ICALL     INTERNAL ERROR CODE (I).
C               = -1, IFMAP=0  SCF CONVERGENCE FAILURE.
C               = -1, IFMAP=1  MO MAPPING FAILURE.
C     *
C     CONVENTION FOR EXTERNAL (OUTPUT) ERROR CODES:
C
C       0:   UNKNOWN FAILURE
C       1:   SCF CONVERGENCE FAILURE
C       2:   MO MAPPING FAILURE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CIMAP / IFMAP
     ./NBFILE/ NBF(20)
C
      IF(NSAV15.LE.0) RETURN
C     INPUT OPTIONS.
      NB15   = NBF(15)

      IF (NSAV15.EQ.5 .OR. NSAV15.EQ.7) THEN
         CALL errsav_bin
         RETURN
      END IF

      IF(NSAV15.LT.9) REWIND NB15
C *** WRITE APPROPRIATE ERROR MESSAGE.
      IF(ICALL.EQ.-1) THEN
         IF(IFMAP.EQ.0) THEN
            WRITE(NB15, 100) 1
            WRITE(NB15, 110)
         ELSE
            WRITE(NB15, 100) 2
            WRITE(NB15, 120)
         ENDIF
      ELSE
         WRITE(NB15, 100) 0
         WRITE(NB15, 130)
      ENDIF
      RETURN
  100 FORMAT(' ERROR: CODE =',I5)
  110 FORMAT(' SCF CONVERGENCE FAILURE')
  120 FORMAT(' MO MAPPING FAILURE')
  130 FORMAT(' UNKNOWN ERROR')

      CONTAINS

      SUBROUTINE errsav_bin
       !---------------------------------------------------------------------
       ! Write unformatted (binary) file for data transfer to ChemShell
       ! Special format conventions for error message:
       ! I:    Length of content vector (= number of sections in file = 1)
       ! I(1): -1, special value to signify that an error has occurred:
       ! I:    Error code (integer)
       !---------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER(4),    PARAMETER      :: Nsec = 1 ! # of sections in file
        CHARACTER(32), PARAMETER      :: fname = 'fort.15'
        INTEGER(4)                    :: items(Nsec) = 0
        INTEGER(4)                    :: ierr
        !--------------------------------------------------------------------
        OPEN(UNIT=NB15, FILE=TRIM(fname), STATUS='REPLACE',
     $        FORM='UNFORMATTED', ACTION='WRITE')

        ! Construct and write content vector
        items(1) = -1           ! Signifies error message

        WRITE(NB15) Nsec
        WRITE(NB15) items

        ! Write appropriate error code
        IF(ICALL.EQ.-1) THEN
           IF(IFMAP.EQ.0) THEN
              ierr = 1
           ELSE
              ierr = 2
           ENDIF
        ELSE
           ierr = 0
        ENDIF

        WRITE(NB15) ierr

        CLOSE(NB15)
        RETURN
      END SUBROUTINE errsav_bin

      END SUBROUTINE ERRSAV
