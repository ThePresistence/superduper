      module module1
      INTERFACE
      SUBROUTINE DMSGRDS (FS,IFS,JFS,GS,IGS,JGS,PS,IPS,JPS,XMU,N,MPDIAG, &
                          NPRINT,NCG,FILT,FILTALL,CUTM,PRECON)
      USE module3
      IMPLICIT NONE
      INTEGER :: N,NCG(10),NPRINT,MPDIAG
      DOUBLE PRECISION :: XMU,CUTM
      DOUBLE PRECISION, DIMENSION (:), POINTER :: GS,FS,PS
      INTEGER, DIMENSION (:), POINTER :: IGS,JGS,IFS,JFS,IPS,JPS
      LOGICAL :: FILT,FILTALL,PRECON
      END SUBROUTINE DMSGRDS
!
      SUBROUTINE DFPUPD (HS,JHS,IHS,GS,JGS,IGS,Q1,JQ1,IQ1,Q2,JQ2,IQ2, &
                         Q9,JQ9,IQ9,N,FILT,FILTALL,CUTM,ICG,NPRINT)
      IMPLICIT NONE
      LOGICAL :: FILT,FILTALL
      INTEGER :: N,NPRINT,ICG
      DOUBLE PRECISION :: CUTM
      DOUBLE PRECISION, DIMENSION (:), POINTER :: GS,HS,Q1,Q2,Q9
      INTEGER, DIMENSION (:), POINTER :: JGS,IGS,JHS,IHS,JQ1,IQ1,JQ2,IQ2,JQ9,IQ9
      END SUBROUTINE DFPUPD 
!
      SUBROUTINE CGUPDS (Q8,JQ8,IQ8,HS,JHS,IHS,GS,JGS,IGS,N,FILTALL, &
                         CUTM,NITER,ICG,NPRINT)
      USE module3
      IMPLICIT NONE
      LOGICAL :: FILTALL
      INTEGER :: N,NPRINT,ICG,NITER
      DOUBLE PRECISION :: CUTM
      DOUBLE PRECISION, DIMENSION (:), POINTER :: Q8,GS,HS
      INTEGER, DIMENSION (:), POINTER :: IQ8,JQ8,JGS,IGS,JHS,IHS
      END SUBROUTINE CGUPDS
!
      SUBROUTINE DMSCOFS (FS,IFS,JFS,HS,IHS,JHS,GS,IGS,JGS,PS,IPS, &
                        JPS,N,NPRINT,B,C,D,FILT,CUTM)
      USE module3
      IMPLICIT NONE
      INTEGER :: N,NPRINT
      DOUBLE PRECISION :: C,B,D,CUTM
      INTEGER, DIMENSION (:), POINTER :: IFS,JFS,IHS,JHS,IGS,JGS,IPS,JPS
      DOUBLE PRECISION, DIMENSION (:), POINTER :: FS,HS,GS,PS
      LOGICAL :: FILT
      END SUBROUTINE DMSCOFS
      END INTERFACE
      end module module1
