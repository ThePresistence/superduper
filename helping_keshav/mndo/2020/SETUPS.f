      SUBROUTINE SETUPS (IFORM,JPRINT)
C     *
C     INPUT FOR MINIMAL CONFIGURATION INTERACTION.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./CCIS  / KMO,LMO,NC,LROOT
     ./INOPT2/ IN2(300)
     ./MOPAC / IMOPAC
     ./NBFILE/ NBF(20)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
C *** FILE NUMBERS.
      NB5    = NBF(5)
      NB6    = NBF(6)
C *** INPUT OPTIONS.
      IMULT  = IN2(66)
C *** INITIALIZATION AND INPUT.
      IF(IMOPAC.LE.0) THEN
         IF(IFORM.LE.0) THEN
            READ(NB5,100) KMO,LMO,NC,LROOT
         ELSE
            READ(NB5,*)   KMO,LMO,NC,LROOT
         ENDIF
      ENDIF
      IF(LROOT.EQ.0) LROOT=1
C *** DEFAULT VALUES FOR 2*2 OR 3*3 CLOSED-SHELL CI.
      IF(IMULT.EQ.0) THEN
         IF(KMO.LE.0 .OR. KMO.GT.NUMB) KMO=NUMB
         IF(LMO.LE.NUMB .OR. LMO.GT.NORBS) LMO=NUMB+1
         IF(NC.NE.2) NC=3
         NMOS   = MAX(NMOS,LMO)
C *** DEFAULT VALUES FOR 3*3 SINGLET CI (HALF ELECTRON METHOD).
      ELSE IF(IMULT.EQ.1) THEN
         KMO    = NUMB-1
         LMO    = NUMB
         NC     = 3
         NMOS   = MAX(NMOS,LMO)
C *** DEFAULT VALUES FOR 2*2 DOUBLET CI (HALF ELECTRON METHOD).
      ELSE IF(IMULT.EQ.2) THEN
         IF(KMO.LE.0 .OR. KMO.GT.NUMB) KMO=NUMB
         IF(LMO.LT.NUMB .OR. LMO.GT.NORBS) LMO=NUMB+1
         NC     = 2
         NMOS   = MAX(NMOS,LMO)
C *** ERROR EXIT.
      ELSE
         WRITE(NB6,120) IMULT
         STOP 'SETUPS'
      ENDIF
C *** PRINT CI DATA.
      IF(JPRINT.GE.0) THEN
         WRITE(NB6,110) KMO,LMO,NC,LROOT
      ENDIF
      RETURN
  100 FORMAT(10I5)
  110 FORMAT(// 1X,'INPUT DATA FOR CONFIGURATION INTERACTION',
     1       /  1X,'EXCITATION FROM ORBITAL       ',I5,
     2       /  1X,'EXCITATION INTO ORBITAL       ',I5,
     3       /  1X,'NUMBER OF CONFIGURATIONS      ',I5,
     4       /  1X,'ROOT FOR GEOMETRY OPTIMIZATION',I5)
  120 FORMAT(// 1X,'CORRELATION TREATMENT ARE NOT AVAILABLE ',
     1             'FOR STATE WITH IMULT =',I2//)
      END