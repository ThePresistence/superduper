C     ******************************************************************
C
C     Auxiliary routines for storage of CPHF K matrix.
C
C     ******************************************************************
C
      FUNCTION IPSBPE(ISTART)
C
C   Compute ending row of a K matrix in a block packed form
C   which would still fit LKATMP buffer, given starting row.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ISTART - Starting row of a block
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGB2 - Parameters of "response" computation.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      All computations are done in double precision to
C      sidestep problem of overflowing default INTEGER
C      type with intermediate values. This could give
C      us answers which are a bit off the point, so we'll
C      adjust them to the right direction using integer
C      arifmetics
C
C   Bugs:
C
C      This is probably too close to brute force solution.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON 
     ./PSDOPT/ DSTORE, DSTEP,  DECONV, DPCONV, DCPHF, DPREC, DCDIFF,
     .         DSHIFT, DBASCR, DNCOFF,
     .         IUMIX,  IURHS,  IUK,    IURES,
     .         IPRINT, IENRG,  ICORE,  IDISK,  IMIX,   IDENS,
     .         INDSYM, IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC,
     .         IROWS,  IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP,
     .         IKMODE, ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT
     ./PSDGB2/ NPAIR, NPAIR2, NOCCA, NVACA, NOCCB, NVACB,
     .         IQSZA, IQSZB,  IQSZ,  ICPV1, ICPVL, NCPVRS
     ./PSPRT / NB6
      SAVE /PSDOPT/, /PSDGB2/, /PSPRT/
C
      DET  = DBLE(ISTART-1)**2 + 4.0D0 * DBLE(IQSZ) * DBLE(IROWS)
      IEND = INT((ISTART-1 + SQRT(DET))/2)
C
      IHAVE = IQSZ*IROWS
      IUSED = IEND*(IEND-ISTART+1)
      IF( IUSED.GT.IHAVE ) THEN
          DO 100 I=IEND,ISTART,-1
              IUSED = I*(I-ISTART+1)
              IF( IUSED.LE.IHAVE ) THEN
                  IPSBPE = I
                  RETURN
              ENDIF
  100     CONTINUE
          WRITE(NB6,10000)
          STOP 'IPSBPE'
      ELSE 
          DO 200 I=IEND+1,IQSZ
              IUSED = I*(I-ISTART+1)
              IF( IUSED.GT.IHAVE ) THEN
                  IPSBPE = I - 1
                  RETURN
              ENDIF
  200     CONTINUE
          IPSBPE = IQSZ
          RETURN
      ENDIF
C
10000 FORMAT(' CANNOT FIT SINGLE ROW IN A BUFFER. PROGRAM WILL STOP.')
      END
C
      SUBROUTINE PSBPSZ(ISIZE,IKSTRP,ISTAT)
C
C   Compute storage requirements for CPHF K matrix held in
C   block packed form.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ISIZE  - On successful return, contains number of words
C               needed for K storage
C      IKSTRP - Number of strips in the CPHF K matrix
C      ISTAT  - Nonzero if arifmetic overflow occured during
C               computation, ISIZE have undefined value.
C
C   Accessed common blocks:
C
C      PSDGB2 - Parameters of "response" computation
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-AK  PARAMETER (IHUGE=2147483647)
C-AK  PARAMETER (IHUGE=9223372036854775807)
      PARAMETER (IHUGE=HUGE(1))
C
      COMMON 
     ./PSDGB2/ NPAIR, NPAIR2, NOCCA, NVACA, NOCCB, NVACB,
     .         IQSZA, IQSZB,  IQSZ,  ICPV1, ICPVL, NCPVRS
      SAVE /PSDGB2/
C
      ISTART = 1
      ISIZE  = 0
      ISTAT  = 0
      IKSTRP = 0
   10 CONTINUE
          IEND   = IPSBPE(ISTART)
          IDELTA = IEND*(IEND-ISTART+1)
          IF( IHUGE-ISIZE.LE.IDELTA ) THEN
              ISTAT = 1
              RETURN
          ENDIF
          ISIZE  = ISIZE + IDELTA
          ISTART = IEND + 1
          IKSTRP = IKSTRP + 1
          IF( ISTART.LE.IQSZ ) GOTO 10
      RETURN
C
      END
C
