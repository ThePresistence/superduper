      SUBROUTINE CHKLRO
C     *
C     MAKE SURE NO EXCITED STATE IS SELECTED IF KCI IS UNSET.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./INOPT1/ IN1(300)
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
      NB6   = NBF(6)
      KCI   = IN2(77)
      LROOT = IN1(140)
      IF(KCI.EQ.0.AND.LROOT.NE.0) THEN
         WRITE(NB6,500)
         WRITE(NB6,510)
         STOP 'CHKLRO'
      ENDIF
  500 FORMAT(/1X,'LROOT SELECTED BUT KCI UNSET.')
  510 FORMAT( 1X,'USE KCI TO SPECIFY AN EXCITED-STATE METHOD.')
      END