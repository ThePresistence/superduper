MODULE CONST

! Physical constants. (Formerly common block CONSTF.)
  REAL(KIND=8), SAVE      :: A0, AFACT, EV, EVCAL, PI, W1, W2, BIGEXP
  REAL(KIND=8), SAVE      :: AU2DBY, AU2WAV, CFAC

! Numerical constants. (Formerly common block CONSTN.)
  REAL(KIND=8), PARAMETER :: ZERO  = 0.0D0
  REAL(KIND=8), PARAMETER :: ONE   = 1.0D0
  REAL(KIND=8), PARAMETER :: TWO   = 2.0D0
  REAL(KIND=8), PARAMETER :: THREE = 3.0D0
  REAL(KIND=8), PARAMETER :: FOUR  = 4.0D0
  REAL(KIND=8), PARAMETER :: PT5   = 0.5D0
  REAL(KIND=8), PARAMETER :: PT25  = 0.25D0

  INTEGER,      PARAMETER :: IZERO = 0
  INTEGER,      PARAMETER :: IONE  = 1

CONTAINS

  SUBROUTINE OLDCF
    A0     =  0.529167D0
    AFACT  = 57.29577951308232D0
    AU2WAV =  2.1947463D5
    CFAC   =  4.803242D0
    EV     = 27.21D0
    EVCAL  = 23.061D0
    PI     =  3.141592653589793D0
    W1     = 14.399D0
    W2     =  7.1995D0
    BIGEXP = 50.0D0
!   Derived constants.
    AU2DBY = A0 * CFAC
  END SUBROUTINE OLDCF

END MODULE CONST
