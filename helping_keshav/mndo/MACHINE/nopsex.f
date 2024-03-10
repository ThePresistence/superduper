      SUBROUTINE SCROSS
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,600)
      STOP 'SCROSS'
      RETURN
  600 FORMAT(//1X,'SURFACE CROSSING NOT INCLUDED ',
     1            'IN THIS DISTRIBUTION. SORRY.')
      END


      SUBROUTINE CONINT
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,600)
      STOP 'CONINT'
      RETURN
  600 FORMAT(//1X,'CONICAL INTERSECTION SEARCH NOT INCLUDED ',
     1            'IN THIS DISTRIBUTION. SORRY.')
      END


      SUBROUTINE YARHES
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,600)
      STOP 'YARHES'
      RETURN
  600 FORMAT(//1X,'CONICAL INTERSECTION SEARCH NOT INCLUDED ',
     1            'IN THIS DISTRIBUTION. SORRY.')
      END


      SUBROUTINE BPARKD
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,600)
      STOP 'BPARKD'
      RETURN
  600 FORMAT(//1X,'CONICAL INTERSECTION SEARCH NOT INCLUDED ',
     1            'IN THIS DISTRIBUTION. SORRY.')
      END


      SUBROUTINE MULSTA
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,600)
      STOP 'MULSTA'
      RETURN
  600 FORMAT(//1X,'SINGLE-POINT CALCULATION OF MULTIPLE CI STATES NOT',
     1            ' INCLUDED IN THIS DISTRIBUTION. SORRY.')
      END


      SUBROUTINE BO_DYNAM
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,600)
      STOP 'BO_DYNAM'
      RETURN
  600 FORMAT(//1X,'BORN-OPPENHEIMER MOLECULAR DYNAMICS NOT',
     1            ' INCLUDED IN THIS DISTRIBUTION. SORRY.')
      END


      SUBROUTINE NACANA
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,600)
      STOP 'NACANA'
      RETURN
  600 FORMAT(//1X,'FULL NON-ADIABATIC COUPLING ELEMENTS NOT',
     1            ' INCLUDED IN THIS DISTRIBUTION. SORRY.')
      END


      SUBROUTINE NACNUM
      COMMON
     ./NBFILE/ NBF(20)
      NB6 = NBF(6)
      WRITE(NB6,600)
      STOP 'NACNUM'
      RETURN
  600 FORMAT(//1X,'FULL NON-ADIABATIC COUPLING ELEMENTS NOT',
     1            ' INCLUDED IN THIS DISTRIBUTION. SORRY.')
      END