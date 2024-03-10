      SUBROUTINE GENSW( R,LO,HI,SW )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION R,LO,HI,WID,X,X2,X3,SW
      IF ( R .GE. HI ) THEN
         SW = 1.0D0
      ELSE IF ( R .LE. LO ) THEN
         SW = 0.0D0
      ELSE
         WID = HI-LO
         IF ( WID .LE. 0.00001D0 ) THEN
            WRITE(6,*)"WID TOO LOW IN GENSW",LO,HI,WID
            STOP
         END IF
C         WRITE(6,*)R,LO,HI,WID
         X = (R-LO)/WID
         X2=X*X
         X3=X2*X
         SW = X3*(10.0D0-15.0D0*X+6.0D0*X2)
      END IF
      END SUBROUTINE
