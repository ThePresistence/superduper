      SUBROUTINE HDLOPT
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,600)
      STOP 'HDLOPT'
      RETURN
  600 FORMAT(//1X,'HDLC OPTIMIZER NOT INCLUDED ',
     1            'IN THIS DISTRIBUTION. SORRY.')
      END
