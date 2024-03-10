C     ******************************************************************
C
C     Modules for automatic selection of options for analytical second
C     derivatives: profile matching and execution time estimation.
C
C     ******************************************************************
      SUBROUTINE PSDMPR(IPROF,INFO)
C
C   Try to match profile IPROF against filled entries in options
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IPROF  - Number of profile to match
C      INFO   - On return, zero if match succeded, 1 if
C               some options do not match, 2 if end of
C               profile table detected, 3 if options
C               incompatibility was detected by PSDCPR.
C
C   Accessed common blocks:
C
C      PSDGBL - Scratch common for global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C      PSDPRF - Computation profiles.
C      PSPRT  - Printing unit.
C
C   Modified common blocks:
C
C      PSDOPT - Computation options. 
C
C   Local storage:
C
C      None.
C
C   Module logic:
C
C      See "Options:" section in PSDRV for description of
C      matching rules and default rules for computing parameters.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IIOPCN=24)
      PARAMETER (IPRFCN=69)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ONE=1.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON 
     ./PSDOPT/ DSTORE, DSTEP,  DECONV, DPCONV, DCPHF, DPREC, DCDIFF,
     .         DSHIFT, DBASCR, DNCOFF,
     .         IUMIX,  IURHS,  IUK,    IURES,
     .         IPRINT, IENRG,  ICORE,  IDISK,  IMIX,   IDENS,
     .         INDSYM, IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC,
     .         IROWS,  IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP,
     .         IKMODE, ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSDGB2/ NPAIR, NPAIR2, NOCCA, NVACA, NOCCB, NVACB, 
     .         IQSZA, IQSZB,  IQSZ,  ICPV1, ICPVL, NCPVRS
     ./PSDPRF/ IPROPT(IIOPCN,IPRFCN)
      SAVE /PSDOPT/, /PSDPRF/, /PSDGBL/, /PSDGB2/
      DIMENSION IOPTS(IIOPCN)
      EQUIVALENCE (IOPTS(1),IPRINT)
C
      INFO = 1
C
C    Setting IPROPT(1) to -666 is a good way to sift out unused
C    entries.
C
      IT = IPROPT(1,IPROF)
      IF( IT.EQ.-666) THEN
C
C         End-of-profile-table was found, skip options matching
C         and try to substitute wildcard options using usual rules.
C
          INFO = 2
      ELSE
          DO 200 IOPT=1,IIOPCN
              IF( IOPT.NE.3 .AND. IOPT.NE.4 .AND. IOPT.NE.6 ) THEN
                  IT = IPROPT(IOPT,IPROF)
                  IF( IOPTS(IOPT).EQ.0 ) THEN
                      IOPTS(IOPT) = IT
                  ELSE IF( IT.NE.0 .AND. IT.NE.IOPTS(IOPT) ) THEN
                      RETURN
                  ENDIF
              ENDIF
  200     CONTINUE
C
          IT = IPROPT(6,IPROF)
          IF( IDENS .EQ.0 .OR. IDENS.EQ.2 .AND. IT.GE.3 ) THEN
              IDENS  = IT
          ELSE IF( IT.NE.0 .AND. IT.NE.IDENS ) THEN
              RETURN
          ENDIF
C
          INFO = 0
      ENDIF
C
C   Process "wildcard" arguments left after matching
C
      IF( IMIX  .EQ.0   ) IMIX   = 3
      IF( IDENS .EQ.0   ) IDENS  = 7
      IF( IDENS .EQ.2   ) IDENS  = 7
      IF( INDSYM.EQ.0   ) INDSYM = 1
      IF( IQSWAP.EQ.0   ) IQSWAP = 2
      IF( INRHS .LT.0   ) THEN
          INRHS  = INT( (DBLE(NCPVRS)+2**(-INRHS-1)-1) / 2**(-INRHS-1) )
          IF( INRHS.LE.0 ) INRHS = 1
      ENDIF
      IF( INRHS.LE.0      ) INRHS = 3
      IF( INRHS.GT.NCPVRS ) INRHS = NCPVRS
      IF( ISOLVE.EQ.0   ) THEN
          IF( INRHS.EQ.1 ) THEN
              IF( IPRECT.EQ.1 ) THEN
                  ISOLVE = 1
              ELSE
                  ISOLVE = 6
              ENDIF
          ELSE
              IF( IPRECT.EQ.1 ) THEN
                  ISOLVE = 2
              ELSE
                  ISOLVE = 6
              ENDIF
          ENDIF
      ENDIF
      IF( IPRECT.EQ.0 .AND. IDENS.NE.1 .AND. IDENS.NE.3  ) THEN
          IF( IDENS.EQ.4 .OR. IDENS.EQ.5 ) THEN
              IPRECT = 4
          ELSE IF( IDENS.GT.5 ) THEN
              IPRECT = 3
          ELSE
              IF( ISOLVE.EQ.1 .OR. ISOLVE.EQ.2 ) THEN
                  IPRECT = 1
              ELSE
                  IPRECT = 2
              ENDIF
          ENDIF
      ENDIF
      IF( IAVEIT.LE.0   ) THEN
          DAVEIT = MAX(LOG10(DBLE(IQSZ))-2.1D0*LOG10(DCPHF)-3.0D0,6.0D0)
          IF( IAVEIT.LT.0 ) DAVEIT = DAVEIT * ABS(IAVEIT)
      ELSE
          DAVEIT = DBLE(IAVEIT)
      ENDIF
      IF( ISOLVE.EQ.4 .OR. ISOLVE.EQ.5 ) DAVEIT = DAVEIT * 5.0D0
      IF( ISOLVE.GT.1 ) THEN
          DAVEIT = DAVEIT*HALF*(ONE+ONE/SQRT(DBLE(INRHS)))
          IF( (ISOLVE.EQ.2 .OR. ISOLVE.EQ.3) .AND. INRHS.GE.3 ) THEN
              DAVEIT = 0.9D0 * DAVEIT
          ENDIF
          DAVEIT = MAX(DAVEIT,5.0D0)
      ENDIF
      IAVEIT = NINT(DAVEIT)
      IF( IMAXIT.LE.0   ) THEN
C         IMAXIT = NINT(10.0D0 + 1.5D0*IAVEIT)
C        The original expression is miscompiled by MIPSPro6 f77 at high
C        optimization levels.
          IMAXIT = 10 + (3*IAVEIT+1)/2
      ENDIF
      IF( IKRSAV.GT.0 ) THEN
          IKRSAV = IKRSAV - 1
      ELSE
          IF( ISOLVE.GE.4 .AND. IAVEIT.GT.2*INRHS ) THEN
              IKRSAV = IAVEIT - INRHS
          ELSE
              IKRSAV = 0
          ENDIF
      ENDIF
      IF( IKRVEC.LE.0   ) THEN
          IF( ISOLVE.EQ.1 .OR. ISOLVE.EQ.2 .OR. ISOLVE.EQ.3 ) THEN
              IF( ISOLVE.EQ.1 ) THEN
                  ITEMP  = MAX(INRHS*(IAVEIT+1),IMAXIT+1)
              ELSE
                  ITEMP  = MAX(INRHS*(IAVEIT+1),(INRHS*(IMAXIT+1))/2,
     .                     IMAXIT+1) + IKRSAV
              ENDIF
              IF( IKRVEC.EQ.0 ) THEN
                  IKRVEC = ITEMP
              ELSE
                  IKRVEC = NINT( (DBLE(ITEMP)*IKRVEC)/10 )
                  IF( IKRVEC.LE.0 ) IKRVEC = 1
              ENDIF
          ENDIF
          IF( ISOLVE.EQ.4 .OR. ISOLVE.EQ.5 .OR. ISOLVE.EQ.6 ) THEN
              IKRVEC = INRHS + IKRSAV
          ENDIF
      ENDIF
      IF( IROWS .LE.0   ) THEN
          IROWS  =  NINT( DBLE(IQSZ) / 2**(-IROWS) )
          IF( IROWS.LE.0 ) IROWS = 1
      ENDIF
      IF( IDSTRP.EQ.0   ) THEN
          IF( IDENS.LT.8 ) THEN
C
C             Make strip arena to occupy as much space as 
C             Krylov vectors did, because it is already here...
C
              IDSTRP = IKRVEC * IQSZ / NROS
          ENDIF
C
C         Default value of 10 should make sure that stage 4 does
C         not take excessive amount of time compared with stage 3.
C         YMMV on different architectures, though.
C
          IDSTRP = MAX(IDSTRP,10)
          IDSTRP = MIN(NCPVRS,IDSTRP)
      ENDIF
C
      IF( IHLST .EQ.0 ) IHLST  = 2
      IF( IHLWRP.EQ.0 ) IHLWRP = 1
      IF( IKMODE.EQ.0 ) IKMODE = 1
C
      IF( DSHIFT.EQ.ZERO ) THEN
          IF( IPRECT.EQ.2 ) THEN
              DSHIFT = 0.15D0
          ELSE IF( IPRECT.EQ.4 ) THEN
              DSHIFT = 0.90D0
          ELSE
              DSHIFT = 1.0D0
          ENDIF
      ENDIF
C
      IF( DBASCR.LE.ZERO ) THEN
          IF( ISOLVE.LE.3 ) THEN
              DBASCR = 2.0D-01
          ENDIF
          IF( ISOLVE.GE.4 ) THEN
              DBASCR = 1.0D-10
          ENDIF
      ENDIF
      IF( DNCOFF.EQ.ZERO ) THEN
          DNCOFF = 1.D-6
      ENDIF
C
      CALL PSDCPR(IPROF,.FALSE.,.FALSE.,IOK)
      IF( IOK.NE.0 ) INFO = 3
C
      RETURN
      END
C
      SUBROUTINE PSDCPR(IPROF,ILPR,ILDIE,INFO)
C
C   Check for inconsistent computation options
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IPROF  - Number of the profile to being processed.
C      ILPR   - .TRUE. if error message should be printed.
C      ILDIE  - .TRUE. if PSDCPR should print error message
C               and exit upon encountering inconsistent settings.
C      INFO   - On return, zero if options are consistent.
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGBL - Scratch common for global computation options.
C      FLAG4  - Integrals storage mode in the calling process.
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
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, LDIE, LPRINT, ILDIE, ILPR
      LOGICAL  DODIP, LIMAG, DOPTCH
      COMMON 
     ./PSDOPT/ DSTORE, DSTEP,  DECONV, DPCONV, DCPHF, DPREC, DCDIFF,
     .         DSHIFT, DBASCR, DNCOFF,
     .         IUMIX,  IURHS,  IUK,    IURES,
     .         IPRINT, IENRG,  ICORE,  IDISK,  IMIX,   IDENS,
     .         INDSYM, IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC,
     .         IROWS,  IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP,
     .         IKMODE, ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSPRT / NB6
      COMMON
     ./INOPT2/ IN2(300)
      SAVE /PSDOPT/, /PSDGBL/, /PSPRT /
C
      IMODE  = IN2(211)
      INFO   = 1
      LPRINT = ILPR
      LDIE   = ILDIE
      LPRINT = LPRINT .OR. LDIE .OR. IPRINT.GT.1
C
C    Final verification of constraints
C
      IF( IDENS.EQ.1 .AND. IMODE.GT.1 ) THEN
          IF(LPRINT) WRITE(NB6,10050) IPROF
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IDENS.NE.8 .AND. IQSWAP.EQ.3 ) THEN
          IF(LPRINT) WRITE(NB6,10060) IPROF, IDENS, IQSWAP
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IDENS.EQ.8 .AND. ISOLVE.EQ.1 .AND. INRHS.NE.1 ) THEN
          IF(LPRINT) WRITE(NB6,10070) IPROF, IDENS
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IENRG.NE.0 .AND. IENRG.NE.1 .AND. IENRG.NE. 4321 ) THEN
          IF(LPRINT) WRITE(NB6,10080) IPROF, IENRG
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IMIX.LT.1 .OR. IMIX.GT.4 ) THEN
          IF(LPRINT) WRITE(NB6,10090) IPROF, IMIX
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IMIX.EQ.2 .AND. IDENS.EQ.8 ) THEN
          IF(LPRINT) WRITE(NB6,10091) IPROF
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( ( IMIX.EQ.4 .AND. METHOD.NE.-1 ) .OR.
     .    ( IMIX.NE.4 .AND. METHOD.EQ.-1 ) ) THEN
          IF(LPRINT) WRITE(NB6,10092) IPROF, IMIX
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      END IF
      IF( IDENS.NE.1 .AND. .NOT.(IDENS.GE.3 .AND. IDENS.LE.8) ) THEN
          IF(LPRINT) WRITE(NB6,10100) IPROF, IDENS
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( INDSYM.NE.1 .AND. INDSYM.NE.2 ) THEN
          IF(LPRINT) WRITE(NB6,10110) IPROF, INDSYM
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IQSWAP.NE.-1.AND..NOT.(IQSWAP.GE.1 .AND .IQSWAP.LE.3) ) THEN
          IF(LPRINT) WRITE(NB6,10120) IPROF, IQSWAP
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IAVEIT.LE.1 ) THEN
          IF(LPRINT) WRITE(NB6,10130) IPROF, IAVEIT
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IMAXIT.LT.IAVEIT ) THEN
          IF(LPRINT) WRITE(NB6,10140) IPROF, IMAXIT
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IDENS.NE.1 .AND. IDENS.NE.3 .AND. ( IKRVEC.LT.INRHS .OR. 
     .        ( (ISOLVE.EQ.1 .OR. ISOLVE.EQ.2 .OR. ISOLVE.EQ.3) .AND. 
     .        IKRVEC.LT.IMAXIT+1 ) ) ) THEN
          IF(LPRINT) WRITE(NB6,10150) IPROF, IKRVEC
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IROWS.LT.1 ) THEN
          IF(LPRINT) WRITE(NB6,10160) IPROF, IROWS
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IDENS.GE.4 .AND. (IPRECT.LT.1 .OR. IPRECT.GT.6) ) THEN
          IF(LPRINT) WRITE(NB6,10180) IPROF, IPRECT
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IPRECT.GT.1 .AND. IDENS.EQ.3 ) THEN
          IF(LPRINT) WRITE(NB6,10185) IPROF
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IPRECT.EQ.4 .AND. IDENS.NE.4 .AND. IDENS.NE.5 ) THEN
          IF(LPRINT) WRITE(NB6,10190) IPROF
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IPRECT.EQ.1 .AND. ISOLVE.GE.3 ) THEN
          IF(LPRINT) WRITE(NB6,10192) IPROF, ISOLVE
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IDSTRP.LT.1 ) THEN
          IF(LPRINT) WRITE(NB6,10195) IPROF, IDSTRP
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IKMODE.NE.1 .AND. IKMODE.NE.2 ) THEN
          IF(LPRINT) WRITE(NB6,10197) IPROF, IKMODE
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( HALFEL.OR.DOCI ) THEN
          IF( IDENS.EQ.1 ) THEN
              IF(LPRINT) WRITE(NB6,10600) IPROF
              IF(LDIE) STOP 'PSDCPR'
              RETURN
          ENDIF
          IF( IDENS.GE.8 ) THEN
              IF(LPRINT) WRITE(NB6,10605) IPROF
              IF(LDIE) STOP 'PSDCPR'
              RETURN
          ENDIF
          IF( IHLST.NE.1 .AND. IHLST.NE.2 ) THEN
              IF(LPRINT) WRITE(NB6,10610) IPROF, IHLST
              IF(LDIE) STOP 'PSDCPR'
              RETURN
          ENDIF
          IF( IHLWRP.NE.1 .AND. IHLWRP.NE.2 ) THEN
              IF(LPRINT) WRITE(NB6,10620) IPROF, IHLWRP
              IF(LDIE) STOP 'PSDCPR'
              RETURN
          ENDIF
      ENDIF
      IF( ISOLVE.LT.1 .OR. ISOLVE.GT.6 ) THEN
          IF(LPRINT) WRITE(NB6,10630) IPROF, ISOLVE
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IKRSAV.LT.0 .OR. IKRSAV.GT.IKRVEC ) THEN
          IF(LPRINT) WRITE(NB6,10640) IPROF, IKRSAV, 0, IKRVEC
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( INRHS.LE.0 ) THEN
          IF(LPRINT) WRITE(NB6,10670) IPROF, INRHS
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
      IF( IDENS.EQ.1 .AND. NPTCHG.NE.0 ) THEN
          IF(LPRINT) WRITE(NB6,10680) 
          IF(LDIE) STOP 'PSDCPR'
          RETURN
      ENDIF
*     IF( IKRSAV.NE.0 .AND. (ISOLVE.EQ.4 .OR. ISOLVE.EQ.5 
*    .                                   .OR. ISOLVE.EQ.6 ) ) THEN
*         IF(LPRINT) WRITE(NB6,10650) IPROF, IKRSAV, ISOLVE
*         IF(LDIE) STOP 'PSDCPR'
*         RETURN
*     ENDIF
C
C    NMR-related restrictions
C
      IF( MODE.EQ.10 ) THEN
          IF( IDENS.EQ.1 ) THEN
              IF(LPRINT) WRITE(NB6,10800)
              IF(LDIE) STOP 'PSDCPR'
          ENDIF
          IF( IDENS.EQ.8 ) THEN
              IF(LPRINT) WRITE(NB6,10810)
              IF(LDIE) STOP 'PSDCPR'
          ENDIF
          IF( NMRLEV.LT.1 .OR. NMRLEV.GT.3 ) THEN
              IF(LPRINT) WRITE(NB6,10820) NMRLEV
              IF(LDIE) STOP 'PSDCPR'
          ENDIF
          IF( NMRLEV.EQ.3 .AND. INTCTL.NE.1   .AND. INTCTL.NE.2
     .                    .AND. .NOT. (INTCTL.GE.11 .AND. INTCTL.LE.19)
     .                    .AND. INTCTL.NE.101 .AND. INTCTL.NE.102 ) THEN
              IF(LPRINT) WRITE(NB6,10830) INTCTL
              IF(LDIE) STOP 'PSDCPR'
          ENDIF
      ENDIF
      INFO = 0
      RETURN
10050 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' DISK-BASED INTERGRALS ARE NOT SUPPORTED BY ',
     .       'EXTRAPOLATION MODULE' )
10060 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' IDENS = ', I4, ' IS NOT COMPATIBLE WITH IQSWAP = ', I4)
10070 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' IDENS = ', I4, ' CAN HANDLE ONLY ONE RIGHT-HAND SIDE ',
     .       'AT TIME.')
10080 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' INVALID IENRG VALUE ', I7 )
10090 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' INVALID IMIX VALUE ', I7 )
10091 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' IMIX=2 IS INCOMPATIBLE WITH IDENS=8')
10092 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' IMIX=4 MUST APPEAR WITH METHOD=-1')
10100 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' INVALID IDENS VALUE ', I7 )
10110 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' INVALID INDSYM VALUE ', I7 )
10120 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' INVALID IQSWAP VALUE ', I7 )
10130 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' INVALID IAVEIT VALUE ', I7 )
10140 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' INVALID IMAXIT VALUE ', I7 )
10150 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' INVALID IKRVEC VALUE ', I7 )
10160 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' INVALID IROWS VALUE ', I7 )
10180 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' INVALID IPRECT VALUE ', I7 )
10185 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' PRECONDITIONING IS ONLY AVAILABLE FOR ITERATIVE SOLVER.')
10190 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' IPRECT = 4 IS ONLY AVAILABLE FOR IDENS = 4,5')
10192 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' ISOLVE = ', I2, ' REQUIRES SYMMETRIZING CONDITIONER.' )
10195 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' INVALID IDSTRP VALUE ', I7 )
10197 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' INVALID IKMODE VALUE ', I7 )
10600 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' SEMI-NUMERIC TREATMENT IS NOT AVAILABLE FOR ',
     .       'HALF-ELECTRON OR CI DERIVATIVES.' )
10605 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' SEPARATE CPHF PHASE IS NECESSARY FOR ',
     .       'HALF-ELECTRON OR CI DERIVATIVES.' )
10610 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' IHLST VALUE (', I5, ') IS INVALID.' )
10620 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' IHLWRP VALUE (', I5, ') IS INVALID.' )
10630 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' ISOLVE VALUE (', I5, ') IS INVALID.' )
10640 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' IKRSAV VALUE (', I5, ') IS OUT OF BOUNDS.',
     .       '  SHOULD BE FROM ', I5, ' TO ', I5, '.')
10650 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' IKRSAV VALUES OTHER THAN 0 (', I5, 
     .       ') ARE INCOMPATIBLE WITH SOLVER ', I5 )
10660 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' IDENS VALUE ', I5, ' IS INVALID.')
10670 FORMAT(' FOR THE PROFILE ', I3, ':',
     .       ' INRHS VALUE ', I5, ' IS INVALID.')
10680 FORMAT(' IDENS = 1 IS NOT VALID IN THE PRESENCE OF ',
     .       'POINT CHARGES.')
10800 FORMAT(' NUMERIC DIFFERENTIATION OF DENSITY MATRIX NOT',
     .       ' IMPLEMENTED FOR CHEMICAL SHIFTS COMPUTATION.' )
10810 FORMAT(' SEPARATE CPHF PHASE IS REQUIRED FOR CHEMICAL SHIFTS',
     .       ' COMPUTATION.')
10820 FORMAT(' NMRLEV VALUE ',I4,' IS INCORRECT, SHOULD BE EITHER 1',
     .       ', 2 OR 3.')
10830 FORMAT(' INTCTL VALUE ',I4,' IS INCORRECT, SHOULD BE 1,',
     .       ' 2, 11, 12, 13, 14, 15, 101, OR 102.')
C
      END
C
      BLOCK DATA PSDMPF
C
C   Initialization of computation profiles and timing constants.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C      PSDOPT - Default computation options
C      PSDPRF - Computation profiles.
C      PSDTMC - Constants for estimation of execution times
C               of various profiles. See PSDEXT for discussion.
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C      ICORE is initialized to the maximum addressable size of
C      scratch array for second derivatives. This might be inappropriate
C      for some configurations which do not have that much physical
C      memory available.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IIOPCN=24)
      PARAMETER (IPRFCN=69)
      PARAMETER (IDCORE=1500000)
      PARAMETER (NTIMEC=33)
C
      COMMON 
     ./PSDOPT/ DSTORE, DSTEP,  DECONV, DPCONV, DCPHF, DPREC, DCDIFF,
     .         DSHIFT, DBASCR, DNCOFF,
     .         IUMIX,  IURHS,  IUK,    IURES,
     .         IPRINT, IENRG,  ICORE,  IDISK,  IMIX,   IDENS,
     .         INDSYM, IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC,
     .         IROWS,  IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP,
     .         IKMODE, ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT
     ./PSDPRF/ IPROPT(IIOPCN,IPRFCN)
     ./PSDTMC/ ATIME(NTIMEC)
      SAVE /PSDPRF/, /PSDTMC/, /PSDOPT/
C
      DATA ICORE/IDCORE/
      DATA IDENS/0/
C
C                   I  I  I  I     I    I  I  I     I  I  I  I  N  I
C                 I N  Q  A  M  I  K  I P  N  D  I  H  K  S  K  M  N
C               I D D  S  V  A  N  R  R R  C  S  H  L  M  O  R  R  T
C               M E S  W  E  X  R  C  O E  P  T  L  W  O  L  S  L  C
C               I N Y  A  I  I  H  E  W C  U  R  S  R  D  V  A  E  T
C               X S M  P  T  T  S  C  S T  S  P  T  P  E  E  V  V  L
C
      DATA ((IPROPT(I,J),I=1,IIOPCN),J=1,19)/
     1 0,0,0,0, 1,1,0,-1, 0, 0, 0, 0, 0,0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0,
     2 0,0,0,0, 3,1,0,-1, 0, 0, 0, 0, 0,0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0,
     3 0,0,0,0, 1,3,0,-1, 0, 0,-1, 0, 0,0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0,
     4 0,0,0,0, 1,3,0,-1, 0, 0,-1, 0, 0,0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0,
     5 0,0,0,0, 3,3,0,-1, 0, 0,-1, 0, 0,0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0,
     6 0,0,0,0, 3,3,0,-1, 0, 0,-1, 0, 0,0, 0, 0, 2, 2, 1, 0, 0, 0, 0, 0,
     7 0,0,0,0, 3,3,0,-1, 0, 0,-1, 0, 0,0, 0, 0, 2, 1, 2, 0, 0, 0, 0, 0,
     8 0,0,0,0, 3,3,0,-1, 0, 0,-1, 0, 0,0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0,
     9 0,0,0,0, 1,4,0,-1, 0, 0,-1, 0, 0,0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0,
     A 0,0,0,0, 1,4,0,-1, 0, 0,-1, 0, 0,0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0,
     1 0,0,0,0, 3,4,0, 2, 0, 0,-1, 0, 0,0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0,
     2 0,0,0,0, 3,4,0, 2, 0, 0,-1, 0, 0,0, 0, 0, 2, 2, 1, 0, 0, 0, 0, 0,
     3 0,0,0,0, 3,4,0, 2, 0, 0,-1, 0, 0,0, 0, 0, 2, 1, 2, 0, 0, 0, 0, 0,
     4 0,0,0,0, 3,4,0, 2, 0, 0,-1, 0, 0,0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0,
     5 0,0,0,0, 3,4,0, 2, 0, 0,12, 0, 0,0, 0, 0, 2, 0, 1, 3, 0, 0, 0, 0,
     6 0,0,0,0, 3,4,0, 2, 0, 0,12, 0, 0,0, 0, 0, 2, 0, 1, 6, 0, 0, 0, 0,
     7 0,0,0,0, 3,4,0, 2, 0, 0,12, 0, 0,0, 0, 0, 2, 0, 2, 3, 0, 0, 0, 0,
     8 0,0,0,0, 3,4,0, 2, 0, 0,12, 0, 0,0, 0, 0, 2, 0, 2, 6, 0, 0, 0, 0,
     9 0,0,0,0, 3,4,0, 2, 0, 0, 1, 0, 0,0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0/
      DATA ((IPROPT(I,J),I=1,IIOPCN),J=20,38)/
     B 0,0,0,0, 3,4,0, 2, 0, 0, 1, 0, 0,0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0,
     1 0,0,0,0, 1,5,0,-1, 0, 0,-1, 0, 0,0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0,
     2 0,0,0,0, 1,5,0,-1, 0, 0,-1, 0, 0,0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0,
     3 0,0,0,0, 3,5,0, 2, 0, 0,-1, 0, 8,0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0,
     4 0,0,0,0, 3,5,0, 2, 0, 0,-1, 0, 8,0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0,
     5 0,0,0,0, 3,5,0, 2, 0, 0,12, 0, 8,0, 0, 0, 2, 0, 1, 3, 0, 0, 0, 0,
     6 0,0,0,0, 3,5,0, 2, 0, 0,12, 0, 8,0, 0, 0, 2, 0, 1, 6, 0, 0, 0, 0,
     7 0,0,0,0, 3,5,0, 2, 0, 0,12, 0, 8,0, 0, 0, 2, 0, 2, 3, 0, 0, 0, 0,
     8 0,0,0,0, 3,5,0, 2, 0, 0,12, 0, 8,0, 0, 0, 2, 0, 2, 6, 0, 0, 0, 0,
     9 0,0,0,0, 3,6,0, 2, 0, 0,-1, 0, 8,0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0,
     C 0,0,0,0, 3,6,0, 2, 0, 0,-1, 0, 8,0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0,
     1 0,0,0,0, 1,7,0,-1, 0, 0,-1, 0, 0,0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0,
     2 0,0,0,0, 1,7,0,-1, 0, 0,12, 0, 0,0, 0, 0, 2, 0, 0, 3, 0, 0, 0, 0,
     3 0,0,0,0, 1,7,0,-1, 0, 0,12, 0, 0,0, 0, 0, 2, 0, 0, 6, 0, 0, 0, 0,
     4 0,0,0,0, 1,7,0,-1, 0, 0, 6, 0, 0,0, 0, 0, 2, 0, 0, 3, 0, 0, 0, 0,
     5 0,0,0,0, 1,7,0,-1, 0, 0, 6, 0, 0,0, 0, 0, 2, 0, 0, 6, 0, 0, 0, 0,
     6 0,0,0,0, 1,7,0,-1, 0, 0, 1, 0, 0,0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0,
     7 0,0,0,0, 3,7,0, 2, 0, 0,-1, 0, 0,0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
     8 0,0,0,0, 3,7,0, 2, 0, 0,12, 0, 0,0, 0, 0, 1, 1, 0, 3, 0, 0, 0, 0/
      DATA ((IPROPT(I,J),I=1,IIOPCN),J=39,57)/
     9 0,0,0,0, 3,7,0, 2, 0, 0,12, 0, 0,0, 0, 0, 1, 1, 0, 6, 0, 0, 0, 0,
     D 0,0,0,0, 3,7,0, 2, 0, 0, 6, 0, 0,0, 0, 0, 1, 1, 0, 3, 0, 0, 0, 0,
     1 0,0,0,0, 3,7,0, 2, 0, 0, 6, 0, 0,0, 0, 0, 1, 1, 0, 6, 0, 0, 0, 0,
     2 0,0,0,0, 3,7,0, 2, 0, 0, 1, 0, 0,0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
     3 0,0,0,0, 3,7,0, 2, 0, 0,-1, 0, 0,0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0,
     4 0,0,0,0, 3,7,0, 2, 0, 0,12, 0, 0,0, 0, 0, 2, 1, 0, 3, 0, 0, 0, 0,
     5 0,0,0,0, 3,7,0, 2, 0, 0,12, 0, 0,0, 0, 0, 2, 1, 0, 6, 0, 0, 0, 0,
     6 0,0,0,0, 3,7,0, 2, 0, 0, 6, 0, 0,0, 0, 0, 2, 1, 0, 3, 0, 0, 0, 0,
     7 0,0,0,0, 3,7,0, 2, 0, 0, 6, 0, 0,0, 0, 0, 2, 1, 0, 6, 0, 0, 0, 0,
     8 0,0,0,0, 3,7,0, 2, 0, 0, 1, 0, 0,0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0,
     9 0,0,0,0, 3,7,0, 2, 0, 0,-1, 0, 0,0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0,
     E 0,0,0,0, 3,7,0, 2, 0, 0,12, 0, 0,0, 0, 0, 1, 2, 0, 3, 0, 0, 0, 0,
     1 0,0,0,0, 3,7,0, 2, 0, 0,12, 0, 0,0, 0, 0, 1, 2, 0, 6, 0, 0, 0, 0,
     2 0,0,0,0, 3,7,0, 2, 0, 0, 6, 0, 0,0, 0, 0, 1, 2, 0, 3, 0, 0, 0, 0,
     3 0,0,0,0, 3,7,0, 2, 0, 0, 6, 0, 0,0, 0, 0, 1, 2, 0, 6, 0, 0, 0, 0,
     4 0,0,0,0, 3,7,0, 2, 0, 0, 1, 0, 0,0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0,
     5 0,0,0,0, 3,7,0, 2, 0, 0,-1, 0, 0,0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0,
     6 0,0,0,0, 3,7,0, 2, 0, 0,12, 0, 0,0, 0, 0, 2, 2, 0, 3, 0, 0, 0, 0,
     7 0,0,0,0, 3,7,0, 2, 0, 0,12, 0, 0,0, 0, 0, 2, 2, 0, 6, 0, 0, 0, 0/
      DATA ((IPROPT(I,J),I=1,IIOPCN),J=58,IPRFCN)/
     8 0,0,0,0, 3,7,0, 2, 0, 0, 6, 0, 0,0, 0, 0, 2, 2, 0, 3, 0, 0, 0, 0,
     9 0,0,0,0, 3,7,0, 2, 0, 0, 6, 0, 0,0, 0, 0, 2, 2, 0, 6, 0, 0, 0, 0,
     F 0,0,0,0, 3,7,0, 2, 0, 0, 1, 0, 0,0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0,
     1 0,0,0,0, 3,8,0, 3, 0, 0,-1, 0, 0,0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0,
     2 0,0,0,0, 3,8,0, 3, 0, 0,12, 0, 0,0, 0, 0, 2, 0, 0, 3, 0, 0, 0, 0,
     3 0,0,0,0, 3,8,0, 3, 0, 0,12, 0, 0,0, 0, 0, 2, 0, 0, 6, 0, 0, 0, 0,
     4 0,0,0,0, 3,8,0, 3, 0, 0, 6, 0, 0,0, 0, 0, 2, 0, 0, 3, 0, 0, 0, 0,
     5 0,0,0,0, 3,8,0, 3, 0, 0, 6, 0, 0,0, 0, 0, 2, 0, 0, 6, 0, 0, 0, 0,
     6 0,0,0,0, 3,8,0, 3, 0, 0, 1, 0, 0,0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0,
     7 0,0,0,0, 4,7,2,-1, 0, 0, 0, 0, 0,0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0,
     8 -666,
     9   0,0,0, 0,0,0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     A -666,1849781345,1954112361,1818321952,1919252073,
     B 1952544361,1936025888,1886217580,1701672046,
     G 1684365856,0544826707,1735532645,1632659555,
     1 1869312886,2037058080,0,0,0,0,0,0,0,0,0/
C
C
C    These values are sppropriate for IBM RS6K/370. Relative
C    magnitudes are probably Ok for all Power machines.
C
C     DATA ATIME/
*    .  .76611155886032547779D-08,  .72710922324223257601D+07,
*    .  .11838547921624482939D-02,  .10879423525578336012D+04,
*    .  .36693498702836045595D+04,  .11284974193280705723D+01,
*    .  .10942932743984086219D+01,  .11000930079563893127D+01,
*    .  .18045458680442425958D+02,  .10000000000000000000D+01,
*    .  .87087510558423062434D+01,  .46460086712684054433D+01,
*    .  .38703005011540798641D+01,  .48600651406926660414D+01,
*    .  .22982385289542750506D+01,  .72167227110144818170D+01,
*    .  .12785491791060922040D+02,  .10919607665286454878D+01,
*    .  .17877096766771857772D+00,  .10954642435439927439D+01,
*    .  .11144453943689307529D+01,  .17034605826276407736D+01,
*    .  .77128152607415927378D+01,  .37812295168295122494D+01,
*    .  .46230404867890708687D+01,  .31738482945008712477D+01,
*    .  .31829045706817566952D+01,  .20104626912501231573D+02,
*    .  .19778086852612055679D+02,  .43865873314535024363D+01,
*    .  .58328612304692240009D+00,  .25160642376600033110D+01,
*    .  .14780469585878728278D+01/ 
C
C    These values are a good bet for SGI Indigo^2 with 100MHz IP22 board
C    (i.e. R4000/R4010 with 1M of secondary cache)
C
C     DATA ATIME/
*    .  .33897714362650199794D-07,  .48221066906416945858D+06,
*    .  .11838547921624482939D-02,  .86239698853449283433D+03,
*    .  .29085617469372659798D+04,  .89425126038847924104D+00,
*    .  .86725551632268471636D+00,  .84986873090337755610D+00,
*    .  .71145222609253027102D+01,  .10000000000000000000D+01,
*    .  .44679280794511679886D+01,  .72414292337728203286D+00,
*    .  .79245669702216703367D+00,  .30384014339932354787D+01,
*    .  .00000000000000000000D+00,  .56781561715572674842D+01,
*    .  .10047597689654859110D+02,  .86427519162607813197D+00,
*    .  .00000000000000000000D+00,  .86702763512536207280D+00,
*    .  .88169808796822335406D+00,  .88800547055933065632D+00,
*    .  .53191001261263561872D+01,  .18167091552453020764D+01,
*    .  .13927983973520423611D+01,  .38706379820380441004D+01,
*    .  .15436511512210506769D+01,  .71411149696758116079D+01,
*    .  .13979900591983339453D+02,  .27860667024716860851D+01,
*    .  .11907431212852809921D+02,  .39170919714717911120D+00,
*    .  .71809574054424074596D+00/ 
C
C    These are good for a single CPU of the SGI Power Challenge.
C    (75 MHz IP21 - R8000/R8010 with 4M of the secondary cache)
C
      DATA ATIME/
     . 0.55067698234993090000D-08, 0.44644311645886106000D+07,
     . 0.11838547921624483000D-02, 0.53086244357101873000D+03,
     . 0.17904112237750444000D+04, 0.55059594112638544000D+00,
     . 0.53390536936663580000D+00, 0.53572030456025188000D+00,
     . 0.59131228426634935000D-02, 0.10000000000000000000D+01,
     . 0.30917082251120809000D+00, 0.51381284070551880000D+00,
     . 0.26914696244563974000D+01, 0.27320855337937466000D+01,
     . 0.13145550430588191000D+01, 0.83915553424762446000D+00,
     . 0.62716229053164216000D+01, 0.53237443969540443000D+00,
     . 0.79411856279178228000D+00, 0.53407424357631117000D+00,
     . 0.54317406256215039000D+00, 0.38588794925264516000D+01,
     . 0.36200478626588999000D+01, 0.52666906111748819000D+01,
     . 0.23209708430448512000D+01, 0.17156617872445636000D-01,
     . 0.64723136094532319000D+00, 0.00000000000000000000D+00,
     . 0.99972588205148671000D+01, 0.55065693209445576000D+00,
     . 0.40081886418433079000D+01, 0.00000000000000000000D+00,
     . 0.00000000000000000000D+00/ 
C
C    Generic values that might be used on other machines without
C    optimizing the constants in the timing model. This ensures
C    that the program will run formally, but certainly not in an
C    optimum manner.
C
C    The first entry should be taken as the inverse of the flop
C    rate (FLOPs/sec) for the LINPACK benchmark on the target
C    machine. The second entry measures the startup time and
C    may be taken to be zero (or, as a slightly better estimate,
C    as 0.1 times the flop rate). The remaining entries are all
C    set equal to one (implying that the machine runs with the
C    same speed for all subtasks).
C
C     DATA ATIME/
*    . .1D-07, 0.0D0, 31*1.0D0/
      END
