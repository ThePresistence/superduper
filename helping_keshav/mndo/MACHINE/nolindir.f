      SUBROUTINE SCFCG
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,600)
      STOP 'SCFCG'
      RETURN
  600 FORMAT(//1X,'FULL-MATRIX DENSITY MATRIX SEARCH NOT INCLUDED ',
     1            'IN THIS DISTRIBUTION. SORRY.')
      END



      SUBROUTINE ITERCG
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,600)
      STOP 'ITERCG'
      RETURN
  600 FORMAT(//1X,'FULL-MATRIX DENSITY MATRIX SEARCH NOT INCLUDED ',
     1            'IN THIS DISTRIBUTION. SORRY.')
      END



      SUBROUTINE DIRCG
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,600)
      STOP 'DIRCG'
      RETURN
  600 FORMAT(//1X,'FULL-MATRIX DENSITY MATRIX SEARCH NOT INCLUDED ',
     1            'IN THIS DISTRIBUTION. SORRY.')
      END
