      SUBROUTINE GUGAFU (IND,NCIGAM,L,M,IP,IQ,IR,IS)
C     *
C     Converts packed index to individual contributions.
C     *
C     Adapted from subroutine PSCIFU written by Serge Pachkovsky.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     IND       Packed index (I).
C     NCIGAM    Size of two-particle density matrix (I).
C     L,M       Indices of the CI matrix element (O), L.LE.M.
C     IP,IQ     Indices for MOs involved (O).
C     IR,IS     IP.LE.IQ, IR.LE.IS, (IP,IQ).LE.(IR,IS)
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TWO=2.D0)
      PARAMETER (C3O4=0.75D0)
C
      INDX   = MOD(IND-1,NCIGAM) + 1
      INDY   = (IND-1)/NCIGAM
      M      = NINT(SQRT(TWO*INDY-C3O4))
      L      = INDY - (M*(M-1))/2
      IRS    = NINT(SQRT(TWO*INDX-C3O4))
      IPQ    = INDX - (IRS*(IRS-1))/2
      IS     = NINT(SQRT(TWO*IRS-C3O4))
      IR     = IRS - (IS*(IS-1))/2
      IQ     = NINT(SQRT(TWO*IPQ-C3O4))
      IP     = IPQ - (IQ*(IQ-1))/2
      RETURN
      END
