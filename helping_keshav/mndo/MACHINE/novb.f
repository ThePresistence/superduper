      SUBROUTINE VBFACE
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,600)
      STOP 'VBFACE'
      RETURN
  600 FORMAT(//1X,'SEMIEMPIRICAL VALENCE-BOND MODULE NOT INCLUDED ',
     1            'IN THIS DISTRIBUTION. SORRY.')
      END
