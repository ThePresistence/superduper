      SUBROUTINE DFTB
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,600)
      STOP 'DFTB'
      RETURN
  600 FORMAT(//1X,'SEMIEMPIRICAL DFTB MODULE NOT INCLUDED ',
     1            'IN THIS DISTRIBUTION. SORRY.')    
      END
      SUBROUTINE DFTPAR
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,610)
      STOP 'DFTB'
      RETURN
  610 FORMAT(//1X,'SEMIEMPIRICAL DFTB MODULE NOT INCLUDED ',
     1            'IN THIS DISTRIBUTION. SORRY.')
      END