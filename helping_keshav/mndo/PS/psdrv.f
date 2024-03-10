C     ******************************************************************
C
C     Driver routines for analytical derivative code
C
C     ******************************************************************
      SUBROUTINE PSDRV(IMOD,DERIV,FORCE,DIPDRV,LDF,A,LDA,ROALP,
     .                 ROBET,CALP,CBET,LDC,EALP,EBET,ONEDEN,TWODEN)
C
C   Compute cartesian energy derivatives (analytically) and second 
C   derivatives (analytically or in mixed mode, depending on control
C   flags) for MNDO/d-type methods. Optionally recompute total energy
C   of molecule to do consistency check with the rest of the code.
C
C   Coded by: Serge Pachkovsky
C
C   Limitations:
C
C      Only MNDO, AM1, PM3 & MNDO/d are currently supported. Non-invariant
C      flavour of MNDO, molecular mechanics corrections etc. are not.
C
C      If second derivatives are computed in mixed mode, two-electron
C      integrals should be in-core.
C
C   Parameters:
C
C      IMOD   - Computation type
C           1 = Compute first derivatives of energy
C           2 = Compute first and second derivatives of energy
C         102 = Compute first and second derivatives of energy and
C               the first derivatives of the dipole moment.
C      DERIV  - Output derivatives, should be at least 3*NATOMS
C               Derivatives are in the atomic units.
C      FORCE  - Output force constants, should be at least 3*NATOMS
C               in second dimension. Not accessed if IMOD.LE.1
C               Force constants are in the atomic units.
C      DIPDRV - Output: Derivatives of the dipole moment. Not accessed
C               if IMOD.NE.102. Dipole moment derivatives are computed
C               in atomic units, multiply by 10^11 e c (in SI units)
C               to convert into Debye/Angstrom.
C      LDF    - Leading dimension of FORCE array, should be at least
C               3*NATOMS. Ignored if IMOD.LE.1
C      A      - Scratch array used by SCFCAL, not accessed if IMOD.LE.1
C      LDA    - Dimension of A array, ignored if IMOD.LE.1
C      ROALP  - Alpha density matrix in packed format
C      ROBET  - Beta density matrix in packed format (not used if 
C               UHF flag not set)
C      CALP   - Alpha orbital coefficients. Not used if IMOD.LE.1
C      CBET   - Beta orbital coefficients (not used if IMOD.LE.1 or
C               UHF flag not set)
C      LDC    - Leading dimension of CALP, CBET matrices
C      EALP   - Energy of alpha molecular orbitals, EV. Not used if
C               IMOD.LE.1
C      EBET   - Energy of beta molecular orbitals, EV. Not used if
C               UHF flag not set or IMOD.LE.1
C      ONEDEN - Precalculated one-particle density matrix.
C               Used only for GUGA-CI (KCI.EQ.5).
C      TWODEN - Precalculated two-particle density matrix.
C               Used only for GUGA-CI (KCI.EQ.5).
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. Although some of the
C               options could be modified by strategy routines,
C               they are restored to the original values on return.
C      PSDYNM - Handles for the dynamic memory blocks
C      PSPRT  - Printing unit
C      PSPRTF - Debug output tuning flags
C
C   Modified common blocks:
C
C      PSDGBL - Scratch common for global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C
C   Local storage:
C
C      IROPCN (=10) DOUBLE PRECISION cells and IIOPCN (=24) INTEGER
C      cells are used to save and restore computation options
C      modified inside PSDRV.
C
C   Other notes:
C
C      PSDGBL, PSDGB2 and PSPRTF are SAVED so that PSDSTR could be 
C      safely separated from PSDRV and be called at the setup phase 
C      instead of computation phase. PSDYNM is OK to leave unSAVE'ed,
C      because it will be filled afresh after call to PSDRV.
C
C      Introducing extra sequence points in memory allocation might
C      eliminate some N^2 arrays. Generally, I am too lax with most
C      N^2 memory items.
C
C   Quirks:
C
C      Some of the parameters are offensive in Russian. This is
C      unintentional.
C
C   Bugs:
C
C      For the call to DSPSV, DOUBLE PRECISION argument is
C      passed for the (INTEGER) IPIV formal parameter. This
C      is violation of standard and could cause problems
C      with some Fortran compilers. Unfortunately, I can't
C      figure a way to cast it to INTEGER while remaining
C      compliant to F77
C
C      Memory and disk allocation routines assume 32-bit integers,
C      and won't attempt to allocate more then 2**31-1 DOUBLE
C      PRECISION words of memory or disk space, even if IDISK and
C      ICORE are set to sufficiently high values. Setting IHUGE
C      to the appropriate value should probably help.
C
C      Passing identical arguments for EALP and EBET, or for CALP
C      and CBET in RHF case will make the code non-conforming to 
C      the F77 standard. This is unlikely to cause any problems,
C      since EBET and CBET are not examined in this case.
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IIOPCN=24)
      PARAMETER (IROPCN=10)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
      PARAMETER (MAXCIO=LMACT)
C
      PARAMETER (SONE  =-1.0D0)
      PARAMETER (ZERO  = 0.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, ISRES, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      COMMON 
     ./INOPT2/ IN2(300)
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
      COMMON
     ./PSDYNM/ LDXPA,  LDXPB,  LDXPDA, LDXPAT, LDXPBT, LDPA,
     .         LDPB,   LHIDAT, LCCNTA, LCCNTB, LPCNTA, LPCNTB, 
     .         LAI2T1, LAI2T2, LAI2T3, LAI2T4, LAI2T5, LKS,
     .         LQ,     LQTMP,  LXTMP,  LIDSPS, LCPITM, LKA, 
     .         LKATMP, LA2MTM, LA2MT1, LCOND,  LO2ATM, LO2AT1,
     .         LCPIT1, LCPAYA, LCPAYB, LCPAY,  LCPATA, LCPATB,
     .         LSHIFT, LCNTMP, LHLO,   LHLQ,   LHLZM,  LHLZA,
     .         LHLQEX, LC1YY,  LC1YR,  LC1RR,  LC1RB,  LC1YB,
     .         LC1BB,  LC1RHS, LC1TMA, LC1TMX, LC1X,   LC1BAS,
     .         LC1RSP, LC1TMP, LC1SNG, LC1SVT, LC1SCR, LHA0,
     .         LHAB,   LH0B,   LCIFI,  LCIFC,  LCIEI,  LCIEC,
     .         LCIVEC, LCIORB, LCIGAM, LCIDLT, LCI2ES, LCIH,
     .         LCIVC1, LCI2TM, LCIINT
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
     ./PSDSTT/ STIMES(-1:9), XTIMES(-1:9)
      COMMON
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN),
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN),
     .         IACTI(MAXACT), IACTJ(MAXACT)
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSPRTF/, /PSDSTT/, /PSHALF/, 
     .     /PSDCIP/, /PSPRT /
C
      DIMENSION DERIV(*), FORCE(LDF,*), A(LDA)
      DIMENSION ROALP(*), ROBET(*), CALP(LDC,*), CBET(LDC,*)
      DIMENSION EALP(*), EBET(*)
      DIMENSION ONEDEN(*), TWODEN(*)
      DIMENSION DIPDRV(LDF,3)
C
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: DUMP
      DIMENSION ROPT(IROPCN), ROPTS(IROPCN), IUNITS(4), 
     .          IOPT(IIOPCN), IOPTS(IIOPCN)
      EQUIVALENCE (DSTORE,ROPT(1)),(IUMIX,IUNITS(1)),(IPRINT,IOPT(1))
C
      DIMENSION DUMMY(1)
C
C    Initialize timing
C
      CALL PSTIME(-2)
C
C    Save options before strategy had time to modify them.
C
      DO 10 I=1,IROPCN
          ROPTS(I) = ROPT(I)
   10 CONTINUE
      DO 20 I=1,IIOPCN
          IOPTS(I) = IOPT(I)
   20 CONTINUE
C
C   Do the following:
C    1. Fetch necessary constants from main program common blocks
C       and initialize vital control flags.
C    2. Initialize tables for fractional occupations and half-electron.
C    3. Complete initialization of exact parameters of computation.
C
      CALL PSDCNS(IMOD)
      IF(DORESP) THEN
           CALL PSHINI(EALP,EBET)
C
C          Confusion danger!
C          After the previous call, EALP and EBET, as well as occupation 
C          numbers (OCCA and OCCB) might be in a different order. If 
C          DORESP flag is still set, they have to be restored by calls
C          to PSUSRX before returning to the calling module.
C
          IF(DORESP) THEN
              CALL PSSRTM(CALP,CBET,LDC)
          ENDIF
      ENDIF
      IF(DOCI) THEN
C
C        First stage of CI initialization. Everything what can be done
C        without dynamic memory allocation is done here. In particular,
C        sizes of all temporaries should be determined on return.
C
          CALL PSCIN1
      ENDIF
      CALL PSDSTR(LDF,LDC,EALP,EBET)
C
C    Zero out derivatives
C
      ENERGY = ZERO
      DO 100 I=1,NVARS
          DERIV(I) = ZERO
  100 CONTINUE
      IF( MODE.GT.1 ) THEN
          DO 120 I=1,NVARS
              DO 110 J=1,NVARS
                  FORCE(I,J) = ZERO
  110         CONTINUE
  120     CONTINUE
          IF(DODIP) THEN
              DO 130 J=1,3
                  DO 128 I=1,NVARS
                      DIPDRV(I,J) = ZERO
  128             CONTINUE
  130         CONTINUE
          ENDIF
      ENDIF
C
C    Check for early return for single particle
C
      IF( NVARS.LE.3 ) GOTO 999
C
C    Memory allocation stage
C
      CALL PSDADM(IMEMST,LDC,ITOP,IDSK,.FALSE.)
      IF( IMEMST.NE.0 ) THEN
          WRITE(NB6,11058)
          STOP 'PSDRV'
      ENDIF
      IF( IPRINT.GE.0) THEN
          IF (ITOP.NE.0) WRITE(NB6,11553) ITOP
          IF (IDSK.NE.0) WRITE(NB6,11554) IDSK
      ENDIF
C
C    Allocate scratch array and determine starting time
C
      ALLOCATE (DUMP(ITOP))
      CALL PSTIME(0)
C
C    Stage 1 - Generate all integrals and store necessary
C              intermediate results
C
      CALL PSMSTG(DUMP,1)
C
C    At this point, dynamic memory is already available, so that
C    we can finish CI initialization.
C
      IF(DOCI) THEN
          CALL PSCIN2(CALP,LDC,EALP,DUMP,ONEDEN,TWODEN)
      ENDIF
      IF( IHLST.EQ.2 ) THEN
          IF(HALFEL) CALL PSHFLO(CALP(1,IOPEN1),LDC,NUMOPN,DUMP)
          IF(DOCI)   CALL PSHFLO(DUMP(IPSMOF(LCIORB)),NORBS,NCIO,DUMP)
      ENDIF
C    This is a bit of a convoluted call. It all comes from the lack
C    of pointers in Fortran 77...
      IF(DOCI) THEN
          CALL PSDST1(ENERGY,DERIV,FORCE,LDF,ROALP,ROBET,
     .                DUMP(IPSMOF(LCIORB)),NORBS,DUMP)
      ELSE IF(HALFEL) THEN
          CALL PSDST1(ENERGY,DERIV,FORCE,LDF,ROALP,ROBET,
     .                CALP(1,IOPEN1),LDC,DUMP)
      ELSE
          CALL PSDST1(ENERGY,DERIV,FORCE,LDF,ROALP,ROBET,CALP,LDC,DUMP)
      ENDIF
C
      IF(DODIP) THEN
          CALL PSDDPS(DIPDRV,LDF,ROALP,ROBET)
      ENDIF
      IF( IENRG.GE.1 .AND. IPRINT.GE.1 ) THEN
          WRITE(NB6,11082) ENERGY
      ENDIF
      CALL PSTIME(1)
C
C    Create direct access file for right-hand sides of CPHF equations.
C
      IF( DORESP .AND. IDENS.NE.1 .AND. IQSWAP.LT.3 
     .                            .AND. IQSWAP.GE.0 ) THEN
          IF( IPRINT.GE.5 ) THEN
              WRITE(NB6,13490) NINT(IQSZ*DSTORE)
          ENDIF
          OPEN(UNIT=IURHS,STATUS='SCRATCH',ACCESS='DIRECT',
     .         FORM='UNFORMATTED',RECL=NINT(IQSZ*DSTORE))
      ENDIF
C
      IF( DORESP .AND. IQSWAP.NE.3 .AND. MODE.GE.2 
     .                             .AND. IDENS.NE.1 ) THEN
C
C     Stage 2 -  Build right-hand sides of CPHF equations. This is not
C                necessary for the half-electron first derivatives, which
C                are taken care of by Z-vector approach.
C
C
          CALL PSMSTG(DUMP,2)
          CALL PSDS1A(CALP,CBET,LDC,DUMP)
          CALL PSTIME(2)
      ENDIF
C
      IF( DORESP .AND. LPPASS ) THEN
          IF( IMIX.NE.3 ) THEN
              DO 150 I=1,NVARS
                  WRITE(NB6,10581) I
                  CALL PSDPPM(NORBS,DUMP(IPSMOF(LDXPA)+NROS*(I-1)))
                  IF( UHF ) THEN
                      WRITE(NB6,10730) I
                      CALL PSDPPM(NORBS,DUMP(IPSMOF(LDXPB)+NROS*(I-1)))
                  ENDIF
  150         CONTINUE
          ENDIF
          IF( IDENS.GE.3 .AND. IQSWAP.EQ.-1 .AND. MODE.GE.2 ) THEN
              WRITE(NB6,10024)
              CALL PSDPGM(IQSZ,NCPVRS,DUMP(IPSMOF(LQ)),IQSZ)
          ENDIF
          WRITE(NB6,11030)
          CALL PSDPDR(NVARS,DERIV)
          IF( MODE.GT.1 ) THEN
              WRITE(NB6,11040)
              CALL PSDPFC(NVARS,FORCE,LDF)
              IF(DODIP) THEN
                  WRITE(NB6,11050)
                  CALL PSDPGM(NVARS,3,DIPDRV,LDF)
              ENDIF
          ENDIF
      ENDIF
C
      IF( DORESP .AND. IDENS.GE.3 .AND. IDENS.LE.7 
     .                            .AND. IMIX.EQ.2 ) THEN
C
C         CPHF equations will be solved on a separate pass.
C         Derivatives of the Fock matrix should be swapped out
C         before we can process with stages 3 to 8.
C
          OPEN(IUMIX,STATUS='SCRATCH',ACCESS='SEQUENTIAL',
     .         FORM='UNFORMATTED')
          REWIND(IUMIX)
          CALL PSDWS(IUMIX,DUMP(IPSMOF(LDXPA)),NVARS*NROS)
          IF( UHF ) THEN
              CALL PSDWS(IUMIX,DUMP(IPSMOF(LDXPB)),NVARS*NROS)
          ENDIF
      ENDIF
C
      IF(DORESP) THEN
          IF( HALFEL.OR.DOCI ) THEN
C
C   Stage 3 - Compute right-hand side for the Z-vector response 
C             problem and store it along with othen right-hand
C             CPHF vectors.
C
              CALL PSMSTG(DUMP,3)
              CALL PSHLQ(CALP,LDC,EALP,DUMP)
              IF(IN2(77).EQ.1) THEN
C                How do you like this beauty?
C                (I'd be glad to pass DUMP, but we want to change type of
C                 IFI and IEI coefficients from dp to int)
                  CALL PSCICK(EALP,DUMP(IPSMOF(LCIFI)),
     .                     DUMP(IPSMOF(LCIFC)),DUMP(IPSMOF(LCIEI)),
     .                     DUMP(IPSMOF(LCIEC)),DUMP(IPSMOF(LCIVEC)),
     .                     DUMP(IPSMOF(LCIH)),DUMP(IPSMOF(LCIVC1)),
     .                     DUMP(IPSMOF(LCIINT)))
              ENDIF
              CALL PSTIME(3)
          ENDIF
C
          IF( IPRECT.GT.1 .AND. (IPRECT.NE.4 .AND. (IDENS.EQ.4 .OR. 
     .                            IDENS.EQ.5) .OR. IDENS.GT.5) ) THEN
C
C   Stage 4 - Build shift vector unless using the "best" shift which 
C             will be computed during integrals transformation.
C
              CALL PSMSTG(DUMP,4)
              CALL PSDS2S(CALP,CBET,LDC,EALP,EBET,DUMP)
              CALL PSTIME(4)
          ENDIF
          IF( IDENS.GE.3 .AND. IDENS.LT.6 ) THEN
C
C   Stage 5 - Explicitly construct CPHF K matrix
C
              CALL PSMSTG(DUMP,5)
              IF( IDENS.EQ.5 .AND. IROWS.LT.IQSZ ) THEN
                  OPEN(IUK,STATUS='SCRATCH',ACCESS='SEQUENTIAL',
     .                     FORM='UNFORMATTED')
                  REWIND(IUK)
              ENDIF
              CALL PSDS2P(CALP,CBET,LDC,EALP,EBET,DUMP)
              CALL PSTIME(5)
              IF( IDENS.EQ.3 .OR. IDENS.EQ.4 ) THEN
C
C                 Have CPHF K matrix in-core
C
                  IF( LPPASS ) THEN
                      WRITE(NB6,11100)
                      CALL PSDPPM(IQSZ,DUMP(IPSMOF(LKS)))
                  ENDIF
              ENDIF
          ENDIF
C
          IF( IDENS.GT.3 ) THEN
C
C   Stage 6 - Compute preconditioner vector using the shift
C             vector computed on stages 4 or 5. We don't need
C             this for the direct solver, which works with the
C             non-conditioned system.
C
              CALL PSMSTG(DUMP,6)
              IF( LPCPHF .AND. IPRECT.NE.1 ) THEN
                  WRITE(NB6,11105)
                  CALL PSXPRT(DUMP(IPSMOF(LSHIFT)))
              ENDIF
              IF( IPRECT.NE.1 ) THEN
                  CALL PSDS3P(EALP,EBET,DUMP(IPSMOF(LCOND)),
     .                                  DUMP(IPSMOF(LSHIFT)))
              ELSE
C
C                 IPRECT=1 is a special case, since it's shift
C                 vector is (implicitly) zero.
C
                  CALL PSDS3P(EALP,EBET,DUMP(IPSMOF(LCOND)),DUMMY)
              ENDIF
              CALL PSTIME(6)
          ENDIF
C
C   Stage 7 - solution of the CPHF equations in a separate phase
C
          IF( IDENS.GE.3 .AND. IDENS.LT.8 ) THEN
              CALL PSMSTG(DUMP,7)
              IF( IQSWAP.EQ.1 ) THEN
C
C                 Bring right-hand CPHF vectors back to memory if they
C                 were swapped out and close scratch file - we do not 
C                 need it anymore
C
                  IQ = IPSMOF(LQ)
                  DO 300 IX=ICPV1,ICPVL
                      CALL PSDRD(IURHS,IX+(1-ICPV1),
     .                           DUMP(IQ+(IX-ICPV1)*IQSZ),IQSZ)
  300             CONTINUE
                  CLOSE(UNIT=IURHS)
                  IF( LPPASS ) THEN
                      WRITE(NB6,11120)
                      CALL PSDPGM(IQSZ,NCPVRS,DUMP(IQ),IQSZ)
                  ENDIF
              ENDIF
              IF( IDENS.EQ.3 ) THEN
C
C                 This one is simple: call LAPACK routine to solve
C                 CPHF equations, and we have completed.
C
                  IF( IQSWAP.GE.2 ) THEN
                      WRITE(NB6,17086)
                      STOP 'PSDRV'
                  ENDIF
C
C                 We favor iterative CPHF solution, which is formulated
C                 for the problem A*X + B = 0. DSPSV solves A*X = B, so
C                 we have to multiply right-hand side by -1 in this case.
C
                  CALL DSCAL(IQSZ*NCPVRS,SONE,DUMP(IPSMOF(LQ)),1)
                  CALL DSPSV('U',IQSZ,NCPVRS,DUMP(IPSMOF(LKS)),
     .                       DUMP(IPSMOF(LIDSPS)),DUMP(IPSMOF(LQ)),
     .                       IQSZ,INFO)
                  IF( INFO.NE.0 ) THEN
                      WRITE(NB6,17455) INFO
                      STOP 'PSDRV'
                  ENDIF
C
C                 We had solved symmetrized system, instead of one we were
C                 interested in. All "real" CPHF vectors (but not Z-vectors)
C                 have to be rescaled to produce the answer we need.
C                 (Since scaling will have to be re-done then computing
C                 derivatives of the density matrix, we can as well leave it
C                 out...)
C
*                 IQ = IPSMOF(LQ)
*                 DO 400 IX=1,ICPVL
*                     CALL PSDUNO(DUMP(IQ+(IX-ICPV1)*IQSZ))
* 400             CONTINUE
              ELSE
C
C                 Iterative solver, either old-fascioned (PSDST3)
C                 or hi-tech (PSLINS)
C
                  IF( ISOLVE.EQ.1 ) THEN
                      CALL PSDST3(CALP,CBET,LDC,EALP,EBET,DUMP)
                  ELSE
                      CALL PSLINS(CALP,CBET,LDC,EALP,EBET,DUMP)
                  ENDIF
                  IF( IDENS.EQ.5 .AND. IROWS.LT.IQSZ ) THEN
                      CLOSE(UNIT=IUK)
                  ENDIF
              ENDIF
              CALL PSTIME(7)
              IF( LPPASS ) THEN
                  IF(IDENS.GE.3.AND.IDENS.LE.7.AND.IQSWAP.LT.2) THEN
C
C                     CPHF solutions are in core, print 'em
C
                      WRITE(NB6,11110)
                      CALL PSDPGM(IQSZ,NCPVRS,DUMP(IPSMOF(LQ)),IQSZ)
                  ENDIF
              ENDIF
          ENDIF
C
          IF( HALFEL.OR.DOCI ) THEN
C
C   Stage 8 - Construct AO-basis representation of the Z-vector now.
C
              CALL PSMSTG(DUMP,8)
              CALL PSHMKZ(CALP,LDC,DUMP)
              CALL PSTIME(8)
          ENDIF
C
C   Stage 9 - Summation of the "response" contributions.
C
          CALL PSMSTG(DUMP,9)
C
C         Bring mixed derivatives back to memory if they were swapped 
C         out during stages 2 and 3. After that, it's Ok to close 
C         scratch file, since we won't need it anymore
C
          IF( IMIX.EQ.2 ) THEN
              REWIND(IUMIX)
              CALL PSDRS(IUMIX,DUMP(IPSMOF(LDXPA)),NVARS*NROS)
              IF( UHF ) THEN
                  CALL PSDRS(IUMIX,DUMP(IPSMOF(LDXPB)),NVARS*NROS)
              ENDIF
              CLOSE(UNIT=IUMIX)
          ENDIF
C
          CALL PSDST4(DERIV,FORCE,DIPDRV,LDF,ROALP,ROBET,DUMP,A,LDA,
     .                    CALP,CBET,LDC,EALP,EBET)
          CALL CPUSEC(T2)
          CALL PSTIME(9)
          IF( IDENS.GE.2 .AND. IQSWAP.EQ.2 ) THEN
C
C             Delete scratch file with CPHF equations solutions
C             before finishing
C
              CLOSE(UNIT=IURHS)
          ENDIF
C
C         Do a final on the integrity of memory blocks
C
          CALL PSMSTG(DUMP,0)
      ENDIF
C
C         Deallocate scratch array
C
      DEALLOCATE (DUMP)
  999 CONTINUE
      IF(DORESP) THEN
C
C         Order of orbital coefficients, occupation numbers and 
C         orbital energies have to be restored now.
C
          CALL PSUSRT(EALP,EBET,CALP,CBET,LDC)
      ENDIF
      IF( IPRINT.GE.1 ) THEN
C
C         Report computed derivatives
C
          IF( NVARS.LT.20 .OR. IPRINT.GE.5 ) THEN
              WRITE(NB6,11146)
              CALL PSDPDR(NVARS,DERIV)
          ENDIF

          IF( MODE.GT.1 ) THEN
              IF( NVARS.LT.20 .OR. IPRINT.GE.5 ) THEN
                  WRITE(NB6,11238)
                  CALL PSDPFC(NVARS,FORCE,LDF)
                  IF(DODIP) THEN
                      WRITE(NB6,11055)
                      CALL PSDPGM(NVARS,3,DIPDRV,LDF)
                  ENDIF
              ENDIF
          ENDIF
      ENDIF
C
      CALL PSTIME(-1)
C
C    Special code: geometry optimization in the presence of the offset
C    force. This should help to do "quasi-ab-initio" geometry 
C    optimizations. 
C
C     CALL PSDROF(NVARS,DERIV,INFO)
C     IF( INFO.GT.0 .AND. 
C    .        ((IPRINT.GE.1.AND.NVARS.LT.20) .OR. IPRINT.GE.5) ) THEN
C         WRITE(NB6,12010)
C         CALL PSDPDR(NVARS,DERIV)
C     ENDIF
C
C    If binary dumps requested, produce it.
C
      IF( IURES.GT.0 ) THEN
          INQUIRE(UNIT=IURES,OPENED=ISRES)
          IF( .NOT.ISRES ) THEN
              OPEN(UNIT=IURES,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     .                        FORM='UNFORMATTED')
              REWIND(IURES)
          ENDIF
          WRITE(IURES) MODE, NVARS, XTIMES(-1)
          WRITE(IURES) (DERIV(I),I=1,NVARS)
          IF( MODE.GE.2 ) THEN
              WRITE(IURES) ((FORCE(I,J),I=1,NVARS),J=1,NVARS)
          ENDIF
          CALL PSBDMP
      ENDIF
C
C    Clean up and exit.
C
      DO 1010 I=1,IROPCN
          ROPT(I) = ROPTS(I)
 1010 CONTINUE
      DO 1020 I=1,IIOPCN
          IOPT(I) = IOPTS(I)
 1020 CONTINUE
      RETURN
11082 FORMAT(' VARIABLE PART OF HEAT OF FORMATION = ',G20.14)
11146 FORMAT(' COMPUTED DERIVATIVES: '/)
11238 FORMAT(' COMPUTED FORCE CONSTANTS: '/)
11553 FORMAT(' MEMORY REQUIRED    :',I10,' WORDS')
11554 FORMAT(' DISK SPACE REQUIRED:',I10,' WORDS')
11058 FORMAT(' NOT ENOUGH DYNAMIC MEMORY OR DISK SPACE TO DO',
     .       ' ANALYTICAL DERIVATIVES',
     .      /' PROGRAM WILL STOP.')
17086 FORMAT(' IQSWAP = 2 IS INCOMPATIBLE WITH LAPACK SOLVER.',
     .      /' PROGRAM WILL STOP.' )
17455 FORMAT(' DSPSV FAILED WITH ERROR CODE ', I8,
     .      /' PROGRAM WILL STOP.' )
13490 FORMAT(' OPENING DIRECT ACCESS FILE FOR RIGHT-HAND SIDES,',
     .       ' RECL = ', I9 )
C
C    Debugging formats
C
10581 FORMAT(' MIXED DERIVATIVE INVOLVING VARIABLE ', I4, 
     .       ' AND ALPHA DENSITY:'/ )
10730 FORMAT(' MIXED DERIVATIVE INVOLVING VARIABLE ', I4, 
     .       ' AND BETA DENSITY:'/ )
10024 FORMAT(' RIGHT-HAND SIDE OF CPHF EQUATIONS FOR VARIABLE (NOT',
     .       ' INCLUDING Z-VECTOR)')
11030 FORMAT(' FIRST DERIVATIVES AFTER STAGE 1:'/)
11040 FORMAT(' SECOND DERIVATIVES AFTER STAGE 1:'/)
11050 FORMAT(' DIPOLE MOMENT DERIVATIVES AFTER STAGE 1 (AU):'/)
11055 FORMAT(' DIPOLE MOMENT DERIVATIVES (AU):'/)
11100 FORMAT(' CPHF K MATRIX (PACKED):'/)
11105 FORMAT(' CPHF SHIFT VECTOR:'/)
11110 FORMAT(' CPHF SOLUTION VECTORS:'/)
11120 FORMAT(' CPHF RIGHT-HAND VECTORS AFTER RELOAD:'/)
12010 FORMAT(' GRADIENTS WITH OFFSET FORCES APPLIED:'/)
      END
C
      SUBROUTINE PSDCNS(IMOD)
C
C   Control options interface routine between derivatives
C   code and main program. This will also set debugging
C   output control flags.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD   - Prepare for:
C           1 = First derivatives computation
C           2 = Second derivatives computation
C          10 = NMR chemical shifts computation.
C         102 = Second derivatives + dipole moment derivatives
C
C   Accessed common blocks:
C
C      ATOMS  - Number of atoms, number of orbitals
C      UHF    - UHF/RHF/Halfelectron
C      HALFE  - RHF/Halfel
C      INOPT2 - Input options
C      ORBITS - Number of alpha and beta electrons in
C               system (note that NORBST in ORBITS
C               corresponds to NORBS elsewhere in program)
C      FLAG5  - Total charge of the system
C      SETCI  - CI flag
C      PSDOPT - Computation options (only output control word,
C               IPRINT) is used.
C      PSMMFL - Options related to point charges.
C
C   Modified common blocks:
C
C      PSDGBL - Global computation parameters.
C      PSDGB2 - Some of the parameters relevant for
C               "response" computation.
C      PSPRT  - Printing unit
C      PSPRTF - Debug output tuning flags.
C      PSO3PR - Precision control for three-center integrals.
C
C   Local storage:
C
C   Module logic:
C
C      This is just a conversion routine, which allows us to
C      localize dependence of the derivatives code on the layout
C      and meaning of variables we can't control.
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL UHF, GUHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
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
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
     ./PSO3PR/ EPSL, EPSR
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSPRTF/, /PSO3PR/, /PSPRT /
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./UHF   / GUHF
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./ORBITS/ NUMB,NORBST,NMOS,NALPHA,NBETA
C    ./PSMMFL/ MMINP,NUMATM,MMCOUP,MMPOT,MMLINK,NLINK,MFILE,MCHARG
C
      NB6    = NBF(6)
      IOP    = IN2(2)
      MMINP  = IN2(12)
      ICHRGE = IN2(65)
      IMULT  = IN2(66)
      KCI    = IN2(77)
      NUMATM = IN2(120)
      MMCOUP = IN2(121)
C     MMPOT  = IN2(122)
C     MMLINK = IN2(123)
C     NLINK  = IN2(124)
C     MFILE  = IN2(125)
C     MCHARG = IN2(126)
C
      MODE   = MOD(IMOD,100)
      DODIP  = .FALSE.
      MODEFL = IMOD/100
      IF( MOD(MODEFL,2) .GE. 1 ) DODIP  = .TRUE.
      IF( MOD(MODEFL,4) .GE. 2 ) DOPTCH = .TRUE.
      IF( MODE.LT.0 ) THEN
          WRITE(NB6,10000) MODE
          STOP 'PSDCNS'
      ENDIF
      IF( MODE.GE.10 ) THEN
          LIMAG = .TRUE.
      ELSE
          LIMAG = .FALSE.
      ENDIF
      IF( DODIP .AND. (MODE.LT.2 .OR. MODE.GE.10) ) THEN
          WRITE(NB6,10005) IMOD
          STOP 'PSDCNS'
      ENDIF
      LPHALF = .FALSE.
      LPCPHF = .FALSE.
      LPPASS = .FALSE.
      LPMMGR = .FALSE.
      LPINTS = .FALSE.
      LPCPHK = .FALSE.
      LPSUMM = .FALSE.
      LPNUME = .FALSE.
      LPA2MO = .FALSE.
      LPSOLV = .FALSE.
      LPNMR  = .FALSE.
      LPCI   = .FALSE.
      IF( IPRINT.GE.16 ) THEN
          IF( MOD(IPRINT,   32).GE.   16 ) LPHALF = .TRUE.
          IF( MOD(IPRINT,   64).GE.   32 ) LPCPHF = .TRUE.
          IF( MOD(IPRINT,  128).GE.   64 ) LPPASS = .TRUE.
          IF( MOD(IPRINT,  256).GE.  128 ) LPMMGR = .TRUE.
          IF( MOD(IPRINT,  512).GE.  256 ) LPINTS = .TRUE.
          IF( MOD(IPRINT, 1024).GE.  512 ) LPCPHK = .TRUE.
          IF( MOD(IPRINT, 2048).GE. 1024 ) LPSUMM = .TRUE.
          IF( MOD(IPRINT, 4096).GE. 2048 ) LPNUME = .TRUE.
          IF( MOD(IPRINT, 8192).GE. 4096 ) LPA2MO = .TRUE.
          IF( MOD(IPRINT,16384).GE. 8192 ) LPSOLV = .TRUE.
          IF( MOD(IPRINT,32768).GE.16384 ) LPNMR  = .TRUE.
          IF( MOD(IPRINT,65536).GE.32768 ) LPCI   = .TRUE.
      ENDIF
C    Minimal set of control variables
      NATOM  = NUMAT
      NORBS  = NLAST(NATOM) 
      NROS   = ( NORBS * (NORBS+1) ) / 2
      NPTCHG = 0
C    The magic incantation for enabling point charge contributions is
C    copied over from SCFCAL, and *should* be identical to that one.
      IF( MMINP.EQ.2 .AND. MMCOUP.GE.2 ) NPTCHG = NUMATM
      NVARS  = 3 * NUMAT
      IF(DOPTCH) NVARS = NVARS + 3*NPTCHG
      UHF    = GUHF
      DOCI   = KCI.EQ.1 .OR. KCI.EQ.5
      HALFEL = IMULT.GT.0 .AND. .NOT. (UHF.OR.DOCI)
      IF( UHF .AND. DOCI ) THEN
          WRITE(NB6,10006)
          STOP 'PSDCNS'
      ENDIF
      IF( IMULT.NE.0 .AND. MODE.EQ.10 ) THEN
          WRITE(NB6,10008)
          STOP 'PSDCNS'
      ENDIF
      IF( (HALFEL.OR.DOCI) .AND. MODE.GT.1 ) THEN
          WRITE(NB6,10010)
          STOP 'PSDCNS'
      ENDIF
      DORESP = .FALSE.
      IF( HALFEL .OR. DOCI .OR. MODE.GT.1 ) THEN
          DORESP = .TRUE.
      ENDIF
      IF( IMIX.EQ.4 ) THEN
          IF (MODE.NE.1) THEN
             WRITE (NB6,"('PSDCNS: HYBRID DERIVATIVES ARE SUPPORTED "//
     .                  "ONLY FOR GRADIENTS.')")
             STOP 'PSDCNS'
          END IF
          METHOD = -1
      ELSE
          IF( IOP.EQ.0 .OR. IOP.EQ.-1 .OR. IOP.EQ.-4 
     .                                .OR. IOP.EQ.-10 ) THEN
              METHOD = 0
          ELSE IF( IOP.EQ.-2 .OR. IOP.EQ.-12 ) THEN
              METHOD = 1
          ELSE IF( IOP.EQ.-7 ) THEN
              METHOD = 2
          ELSE IF( IOP.EQ.-20 .OR. IOP.EQ.-21 ) THEN
              METHOD = 3
          ELSE
              WRITE(NB6,10020) IOP
              STOP 'PSDCNS'
          ENDIF
      ENDIF
      IF( MODE.EQ.0 .AND. METHOD.NE.3 ) THEN
          WRITE(NB6,10022)
          STOP 'PSDCNS'
      ENDIF
      IF( IPRINT.GE.5 ) THEN
          WRITE(NB6,11100) MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, UHF,
     .                   HALFEL, DOCI, METHOD, ICHRGE, DORESP,
     .                   DODIP, LIMAG, DOPTCH
      ENDIF
      NOCCA  = NALPHA
      NVACA  = NORBS - NOCCA
      IF( UHF ) THEN
          NOCCB  = NBETA
          NVACB  = NORBS - NOCCB
      ELSE
          NOCCB  = 0
          NVACB  = 0
      ENDIF
      EPSL   = 2D-6
      EPSR   = 1D-6
      RETURN
10000 FORMAT(/' GARBAGE MODE PARAMETER ', I5, ' PASSED TO PSDRV'
     .       /' JOB WILL BE TERMINATED' )
10005 FORMAT(/' DIPOLE MOMENT DERIVATIVES CANNOT BE COMPUTED WITHOUT',
     .        ' SECOND ENERGY DERIVATIVES.'
     .       /' IMOD VALUE ', I5, ' IS HENCE INVALID.')
10006 FORMAT(/' UNRECTRICTED CI IS NOT SUPPORTED.')
10008 FORMAT(/' NMR CHEMICAL SHIFTS ARE NOT SUPPORTED FOR OPEN-SHELL',
     .        ' SYSTEMS.')
10010 FORMAT(/' SECOND ANALYTIC DERIVATIVES ARE NOT AVAILABLE FOR ',
     .        'HALF-ELECTRON METHOD AND CI'
     .       /' JOB WILL BE TERMINATED' )
10020 FORMAT(/' ANALYTIC DERIVATIVES ARE NOT AVAILABLE FOR METHOD ',
     .        'SELECTED BY IOP = ', I5
     .       /' JOB WILL BE TERMINATED' )
10022 FORMAT(/' ONE-ELECTRON HAMILTONIAN CONSTRUCTION IS ONLY SUPPORTED'
     .       ,' FOR OM1/OM2 METHODS.')
11100 FORMAT( ' MODE   = ', I5, ' NATOM  = ', I5, ' NPTCHG = ', I5,
     .        ' NVARS  = ', I5, ' NORBS  = ', I5, ' NROS   = ', I5/
     .        ' UHF    = ', L5, ' HALFEL = ', L5, ' DOCI   = ', L5,
     .        ' METHOD = ', I5, ' ICHRGE = ', I5, ' DORESP = ', L5/
     .        ' DODIP  = ', L5, ' LIMAG  = ', L5, ' DOPTCH = ', L5 )
      END
C
      SUBROUTINE PSDSTR(LDF,LDC,EALP,EBET)
C
C   Strategy dispatcher for PSDRV: complete filling global computation
C   parameters (begun by PSDCNS) and determine optimal computation path 
C   for PSDRV.  Do minimal check for the inconsistencies in parameters 
C   which could cause garbage output.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      LDF    - Leading dimension of (output) force
C               constants array, check for consistency
C      LDC    - Leading dimension of (input) orbital
C               coefficients array, check for consistency
C      EALP   - Alpha orbital energies, checked for
C               sufficient occ-vac separation. Not used
C               if DORESP is FALSE or semi-numeric mode selected.
C      EBET   - Beta orbital energies, checked for
C               sufficient occ-vac separation. Not used
C               if DORESP is FALSE, IDENS.EQ.1, or .NOT.UHF
C 
C   Accessed common blocks:
C
C      PSDPRF - To force linker to barf if I'll change number
C               of profiles in psdmpr.f, but not here.
C      PSOCC  - Occupation number blocks.
C      PSPRT  - Printing unit.
C      ATOMS  - Number of atoms, number of orbitals
C      LIMITS - Leading array dimensions, to make sure
C               we would be able to deal with situation
C               later on
C      SCRT   - Convergence criteria used during during
C               previous SCF computation
C
C   Modified common blocks:
C
C      PSDOPT - Detailed computation options compatible with
C               values already filled in at the time of call
C      PSDGBL - Global computation parameters
C      PSDGB2 - Global computation parameters used by response
C               quantities computations.
C      PSDGB3 - Global computation parameters for the execution
C               time estimations.
C
C   Local storage:
C
C      IROPCN (=10) DOUBLE PRECISION cells and IIOPCN (=24) INTEGER
C      cells to save/restore original computation options between
C      calls to profile matching routine.
C
C   Module logic:
C
C      This module tries to match set of standard profiles against
C      partially specified computation options and memory/disk use
C      limits, selecting profile wich gives lowest expectation of
C      execution time. If two profiles have the same expectation of
C      execution time, one which have smaller memory and disk storage
C      requirements (in that order) is selected.
C
C   Bugs:
C
C      Passing garbage parameters could cause jobs termination
C      in PSDSTR. Unfortunately, Fortran do not provide exception
C      handling, so that I can't deal with these errors gracefully.
C      Stopping now is better then facing unsupported parameter
C      combination after couple of hours of number-crunching, though.
C
C      Selection of profiles is based upon fitted execution-time
C      prediction, so it could select inferior method if timed for
C      wastly different system or compiler.
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IIOPCN=24)
      PARAMETER (IROPCN=10)
      PARAMETER (IPRFCN=69)
C-AK  PARAMETER (HUGE=1D30)
      PARAMETER (DHUGE=HUGE(1.D0))
C-AK  PARAMETER (IHUGE=2147483647)
C-AK  PARAMETER (IHUGE=9223372036854775807)
      PARAMETER (IHUGE=HUGE(1))
      PARAMETER (ZERO=0.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL CPRREJ
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
     ./PSDGB3/ NOPAIR, INPAIR
     ./PSPRT / NB6
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSDGB3/, /PSPRT /
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./LIMITS/ LM2,LM3,LM4,LM6,LM7,LM8,LM9
     ./SCRT  / SCFCRT,PLCRT,SCFCRS(2)
C
      DIMENSION EALP(*), EBET(*)
      DIMENSION ROPT(IROPCN), ROPTS(IROPCN),
     .          IOPT(IIOPCN), IOPTS(IIOPCN)
      DIMENSION CPRREJ(IPRFCN)
      EQUIVALENCE (DSTORE,ROPT(1)),(IPRINT,IOPT(1))
C
C    Compute sum of the number of orbital pairs on each atom
C    and sum of the squared number of orbital pairs. These are
C    necessary to compute various limits involving number of
C    two-electron integrals.
C
      NPAIR  = 0
      NPAIR2 = 0
      NOPAIR = 0
      INPAIR = 0
      DO 500 I=1,NATOM
          IORBS  = NLAST(I) - NFIRST(I) + 1
          IPAIR  = (IORBS*(IORBS+1))/2
          NPAIR  = NPAIR + IPAIR
          NPAIR2 = NPAIR2 + IPAIR**2
          NOPAIR = NOPAIR + IORBS*IPAIR
          INPAIR = INPAIR + IPAIR*I
  500 CONTINUE
C
      IF( DSTEP .EQ.ZERO ) DSTEP  = 2.0D-4
C
C    This is all that's necessary for first derivatives
C
      IF(.NOT.DORESP) RETURN
C
C    Number of CPHF equations which need to be solved
C
      ICPV1 = 1
      ICPVL = 0
      IF( HALFEL.OR.DOCI ) ICPV1 = 0
      IF( MODE.EQ. 2 ) ICPVL = NVARS
      IF( MODE.EQ.10 ) ICPVL = 3
      NCPVRS = MAX(0,ICPVL-ICPV1+1)
C
      IF( IPRINT.GE.5 ) THEN
          WRITE(NB6,11130) NOCCA, NVACA, NOCCB, NVACB, IQSZA, IQSZB,
     .                     IQSZ, NPAIR, NPAIR2, ICPV1, ICPVL
          WRITE(NB6,11110) LDF, LDC, LM2, LM3
      ENDIF
      IF( MODE.EQ.2 .AND. LDF.LT.NVARS ) THEN
C
C    If second derivatives were requested, enough space should be
C    available...
C
          WRITE(NB6,10030) LDF, NVARS
          STOP 'PSDSTR'
      ENDIF
      IF( NORBS.GT.LDC .OR. LM3.NE.NORBS ) THEN
C
C    Orbital coefficients matrix should be at least as large in the
C    first dimension as there are basis orbitals.
C
C    There should be as much eigenvalues/eigenvectors as there
C    are orbitals - we have no use for occupied-only solution.
C
          WRITE(NB6,10040) 
          STOP 'PSDSTR'
      ENDIF
C
      IF( IDENS .EQ.0 ) IDENS  = 2
C
      IF( IUMIX .EQ.0 ) IUMIX  = 95
      IF( IURHS .EQ.0 ) IURHS  = 96
      IF( IUK   .EQ.0 ) IUK    = 97
      IF( DSTORE.EQ.0 ) DSTORE =  8
C
      IF( DPREC .EQ.ZERO ) THEN
          IF( DOCI.OR.HALFEL ) THEN
              DPREC = 1.0D-6
          ELSE
              DPREC = 1.0D-5
          ENDIF
      ENDIF
      IF( DCDIFF.EQ.ZERO ) DCDIFF = 1.0D-1
      IF( MAX(SCFCRT,PLCRT).GT.DPREC ) THEN
          WRITE(NB6,10300) DPREC
      ENDIF
C
C         Pick parameters which will give desired precision.
C
      IF( DCPHF .EQ.ZERO ) DCPHF  = DPREC 
      IF( DECONV.EQ.ZERO ) DECONV = MAX(DPREC*DSTEP,SCFCRT)
      IF( DPCONV.EQ.ZERO ) DPCONV = MAX(DPREC*DSTEP,PLCRT)
      IF( IPRINT.GE.5 ) THEN
        WRITE(NB6,11119)
        WRITE(NB6,11120) IENRG,  ICORE,  IDISK,  IMIX,   IDENS,  INDSYM,
     .                   IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC, IROWS,
     .                   IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP, IKMODE,
     .                   ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT,
     .                   IUMIX,  IURHS,  IUK,    IURES, 
     .                   DSTORE, DSTEP,  DECONV, DPCONV, DCPHF,  DPREC,
     .                   DCDIFF, DSHIFT, DBASCR, DNCOFF
      ENDIF
C
C    Match options against available profiles, select fastest
C    possible match.
C
      DO 900 I=1,IIOPCN
          IOPTS(I) = IOPT(I)
  900 CONTINUE
      DO 902 I=1,IROPCN
          ROPTS(I) = ROPT(I)
  902 CONTINUE
C
C    "Best overall" profile
C
      IBEST  = 0
      TBEST  = DHUGE
      IMBEST = IHUGE
      IDBEST = IHUGE
C
C    "Fastest overall" profile
C
      IFAST  = 0
      TFAST  = DHUGE
      IMFAST = IHUGE
      IDFAST = IHUGE
C
C    "Smallest overall" profile
C
      ITINY  = 0
      TTINY  = DHUGE
      IMTINY = IHUGE
      IDTINY = IHUGE
C
      DO 1000 I=1,IPRFCN
          CALL PSDMPR(I,INFO)
          CPRREJ(I) = INFO.EQ.3
          IF( INFO.EQ.0 .OR. (INFO.EQ.2 .AND. IBEST.EQ.0) ) THEN
              IF( IPRINT.GE.5 ) THEN
                  WRITE(NB6,11300) I
              ENDIF
              CALL PSDADM(INFO,LDC,ITOP,IDSK,.TRUE.)
              IF( IPRINT.GE.5 ) THEN
                  WRITE(NB6,11310) I, INFO, ITOP, IDSK
              ENDIF
              IF( INFO.EQ.0 .OR. INFO.EQ.1 ) THEN
                  CALL PSDEXT(TIME)
                  IF( IPRINT.GE.5 ) THEN
                      WRITE(NB6,11320) I, TIME
                  ENDIF
                  IF( TIME.GT.DHUGE ) TIME = DHUGE
                  IF( TIME.LT.TFAST .OR.
     .                TIME.EQ.TFAST .AND. ITOP.LT.IMFAST .OR.
     .                TIME.EQ.TFAST .AND. ITOP.EQ.IMFAST .AND. 
     .                                    IDSK.LT.IDFAST ) THEN
                      IFAST  = I
                      TFAST  = TIME
                      IMFAST = ITOP
                      IDFAST = IDSK
                  ENDIF
                  IF( ITOP.LT.IMTINY .OR.
     .                ITOP.EQ.IMTINY .AND. TIME.LT.TTINY .OR.
     .                ITOP.EQ.IMTINY .AND. TIME.EQ.TTINY .AND.
     .                                     IDSK.LT.IDTINY ) THEN
                      ITINY  = I
                      TTINY  = TIME
                      IMTINY = ITOP
                      IDTINY = IDSK
                  ENDIF
                  IF( INFO.EQ.0 ) THEN
                      IF( TIME.LT.TBEST .OR.
     .                    TIME.EQ.TBEST .AND. ITOP.LT.IMBEST .OR.
     .                    TIME.EQ.TBEST .AND. ITOP.EQ.IMBEST .AND. 
     .                                        IDSK.LT.IDBEST ) THEN
                          IBEST  = I
                          TBEST  = TIME
                          IMBEST = ITOP
                          IDBEST = IDSK
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF
          DO 990 J=1,IIOPCN
              IOPT(J) = IOPTS(J)
  990     CONTINUE
          DO 992 J=1,IROPCN
              ROPT(J) = ROPTS(J)
  992     CONTINUE
 1000 CONTINUE
      IF( IPRINT.GE.1 .AND. IFAST.NE.0 ) THEN
          WRITE(NB6,10500) 'FASTEST', IFAST, TFAST, IMFAST, IDFAST
      ENDIF
      IF( IPRINT.GE.1 .OR. (IBEST.EQ.0 .AND. ITINY.NE.0) ) THEN
          WRITE(NB6,10500) 'SMALLEST', ITINY, TTINY, IMTINY, IDTINY
      ENDIF
      IF( IBEST.EQ.0 ) THEN
          WRITE(NB6,10200)
          DO 1100 I=1,IPRFCN
              IF( CPRREJ(I) ) THEN
                  CALL PSDMPR(I,ISTAT)
                  CALL PSDCPR(I,.TRUE.,.FALSE.,ISTAT)
                  DO 1090 J=1,IIOPCN
                      IOPT(J) = IOPTS(J)
 1090             CONTINUE
                  DO 1092 J=1,IROPCN
                      ROPT(J) = ROPTS(J)
 1092             CONTINUE
              ENDIF
 1100     CONTINUE
          STOP 'PSDSTR'
      ENDIF
      CALL PSDMPR(IBEST,INFO)
C
C    Call to PSDEXT will compute correct per-stage execution times
C
      CALL PSDEXT(TTEMP)
      IF( IPRINT.GE.0 ) THEN
          WRITE(NB6,10250) IBEST, TBEST
          IF( IDENS.GT.3 ) WRITE(NB6,10260) IAVEIT
      ENDIF
C
      IF( IPRINT.GE.0 ) THEN
          IF( IDENS.GE.3 ) THEN
              IF( DCPHF .GT.DPREC ) THEN
                  WRITE(NB6,10310) DPREC
              ENDIF
          ENDIF
          IF( IDENS.EQ.1 .AND. MAX(DECONV,DPCONV).GT.DPREC*DSTEP ) THEN
              WRITE(NB6,10320) DPREC
          ENDIF
      ENDIF
C
C    Check energy separation of occupation blocks. SCFCRT/DCPHF is the
C    minimal separation which is still acceptable to get DCPHF precision.
C    This call will check for near-degeneracy and block overlap as well,
C    which could potentially result in more severe consequences than 
C    loss of precision, so it should not be conditioned on IPRINT.
C
      CALL PSCHKE(EALP,EBET,SCFCRT/DCPHF)
C
      IF( IPRINT.GE.5 ) THEN
        WRITE(NB6,11120) IENRG,  ICORE,  IDISK,  IMIX,   IDENS,  INDSYM,
     .                   IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC, IROWS,
     .                   IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP, IKMODE,
     .                   ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT,
     .                   IUMIX,  IURHS,  IUK,    IURES, 
     .                   DSTORE, DSTEP,  DECONV, DPCONV, DCPHF,  DPREC,
     .                   DCDIFF, DSHIFT, DBASCR, DNCOFF
      ENDIF
C
C    Die if conflicting options were selected
C
      CALL PSDCPR(0,.TRUE.,.TRUE.,ISTAT)
      RETURN
10000 FORMAT(/' GARBAGE MODE PARAMETER ', I5, ' PASSED TO PSDRV'
     .       /' JOB WILL BE TERMINATED' )
10010 FORMAT(/' ANALYTIC DERIVATIVES ARE NOT AVAILABLE FOR ',
     .        'HALF-ELECTRON METHOD'
     .       /' JOB WILL BE TERMINATED' )
10020 FORMAT(/' ANALYTIC DERIVATIVES ARE NOT AVAILABLE FOR METHOD ',
     .        'SELECTED BY IOP = ', I5
     .       /' JOB WILL BE TERMINATED' )
10030 FORMAT(/' FORCE CONSTANT MATRIX IS TOO SMALL TO HOLD ALL ',
     .        'VALUES: HAVE ', I4, ' NEED ', I4
     .       /' JOB WILL BE TERMINATED' )
10040 FORMAT(/' UNSUPPORTED MEMORY MODE IN CALLING MODULE'
     .       /' JOB WILL BE TERMINATED' )
10200 FORMAT(/' CANNOT FIND ROUTE WHICH SATISFIES ALL GIVEN ',
     .        'CONSTRAINTS.'
     .       /' THIS WAS CAUSED BY THE INSUFFICIENT AMOUNT ',
     .        ' MEMORY OR DISK SPACE ALLOCATED TO THE PROGRAM,'
     .       /' OR ONE OF THE OPTIONS MISMATCHES REPORTED BELOW.')
10250 FORMAT(/' SELECTED PROFILE ', I3, ' WITH EXPECTED EXECUTION',
     .        ' TIME ', F17.3, ' SECONDS ' )
10260 FORMAT( ' ESTIMATED NUMBER OF CPHF ITERATIONS PER VARIABLE: ',I5)
10300 FORMAT(/' WARNING: PRECISION OF SCF MAY BE INSUFFICIENT',
     .       /' TO GIVE THE DESIRED ACCURACY OF ', G14.5,
     .        ' FOR THE ANALYTICAL DERIVATIVES')
10310 FORMAT(/' WARNING: PRECISION OF CPHF PROCESS MAY BE INSUFFICIENT',
     .       /' TO GIVE THE DESIRED ACCURACY OF ', G14.5,
     .        ' FOR THE ANALYTICAL DERIVATIVES')
10320 FORMAT(/' WARNING: PRECISION OF NUMERIC STEPS MAY BE INSUFFICIENT'
     .       /' TO GIVE THE DESIRED ACCURACY OF ', G14.5,
     .        ' FOR THE ANALYTICAL DERIVATIVES')
10500 FORMAT(/' ', A8, ' AVAILABLE PROFILE (', I3, ') WITH EXPECTED',
     .        ' EXECUTION TIME OF ', F12.2, ' SECONDS'
     .       /' REQUIRES ', I12, ' WORDS OF MEMORY AND ', I12,
     .        ' WORDS OF DISK SPACE' )
C
C    Formats for debugging output - don't waste your time
C    to prettify them.
C
11100 FORMAT( ' MODE   = ', I5, ' NATOM  = ', I5, ' NVARS  = ', I5,
     .        ' NORBS  = ', I5, ' NROS   = ', I5/
     .        ' UHF    = ', L5, ' HALFEL = ', L5, ' METHOD = ', I5 )
11110 FORMAT( ' LDF    = ', I5, ' LDC    = ', I5, ' LM2    = ', I5,
     .        ' LM3    = ', I5 )
11119 FORMAT( ' LOOKING FOR THE BEST MATCH FOR OPTIONS:' )
11120 FORMAT( ' IENRG  = ', I5/
     .        ' ICORE  = ', I9, ' IDISK  = ', I9/
     .        ' IMIX   = ', I5, ' IDENS  = ', I5, ' INDSYM = ', I5,
     .        ' IQSWAP = ', I5, ' IAVEIT = ', I5/
     .        ' IMAXIT = ', I5, ' INRHS  = ', I5, ' IKRVEC = ', I5, 
     .        ' IROWS  = ', I5, ' IPRECT = ', I5/
     .        ' INCPUS = ', I5, ' IDSTRP = ', I5, ' IHLST  = ', I5,
     .        ' IHLWRP = ', I5, ' IKMODE = ', I5/
     .        ' ISOLVE = ', I5, ' IKRSAV = ', I5, ' NMRLEV = ', I5,
     .        ' INTCTL = ', I5, ' ICIOPT = ', I5/
     .        ' IUMIX  = ', I5, ' IURHS  = ', I5, ' IUK    = ', I5,
     .        ' IURES  = ', I5/
     .        ' DSTORE = ', G20.14, ' DSTEP  = ', G20.14, 
     .        ' DECONV = ', G20.14/
     .        ' DPCONV = ', G20.14, ' DCPHF  = ', G20.14,
     .        ' DPREC  = ', G20.14/ 
     .        ' DCDIFF = ', G20.14, ' DSHIFT = ', G20.14,
     .        ' DBASCR = ', G20.14/ 
     .        ' DNCOFF = ', G20.14) 
11130 FORMAT( ' NOCCA  = ', I5, ' NVACA  = ', I5, ' NOCCB  = ', I5,
     .        ' NVACB  = ', I5/
     .        ' IQSZA  = ', I5, ' IQSZB  = ', I5, ' IQSZ   = ', I5/
     .        ' NPAIR  = ', I5/
     .        ' NPAIR2 = ', I5, ' ICPV1  = ', I5, ' ICPVL  = ', I5 )
11300 FORMAT( ' PROFILE ', I3, ' MATCHES VALUES GIVEN' )
11310 FORMAT( ' PSDADM FOR PROFILE ', I3, ' COMPLETED WITH STATUS ',
     .        I4, 
     .        ' AND ITOP = ', I10, ', IDSK = ', I10 )
11320 FORMAT( ' TIME ESTIMATION FOR PROFILE ', I3, ' IS ', G12.5 )
      END
C
      SUBROUTINE PSDWS(IUNIT,A,ISIZE)
C
C   Write variable-sized array to sequential file
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IUNIT  - Unit to write array to. It should be already
C               connected for unformatted sequential access
C      A      - Array to write
C      ISIZE  - Size of the array
C 
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Implied-DO construct could have a significant memory
C      overhead in naively implemented compiler and is better
C      avoided if large amounts of data have to be transferred.
C      At the same time, Fortran does not allow array I/O for
C      arrays of unspecified size, so we have to pass it to
C      a very simple routine as a dummy argument.
C
C      Also, if particular Fortran implementation would have
C      unreasonable limitations on record length for the 
C      unformatted file, we would be able to patch for it in
C      an uniform manner
C
C   Notes:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(ISIZE)
C
      WRITE(IUNIT) A
C
      RETURN
      END
C
      SUBROUTINE PSDRS(IUNIT,A,ISIZE)
C
C   Read variable-sized array from sequential file
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IUNIT  - Unit to write array to. It should be already
C               connected for unformatted sequential access
C      A      - Array to fill
C      ISIZE  - Size of the array
C 
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      See PSDWS for motivation
C
C   Notes:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(ISIZE)
C
      READ(IUNIT) A
C
      RETURN
      END
C
      SUBROUTINE PSDWD(IUNIT,IREC,A,ISIZE)
C
C   Write variable-sized array to direct access file
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IUNIT  - Unit to write array to. It should be already
C               connected for unformatted direct access
C      IREC   - Record number
C      A      - Array to write
C      ISIZE  - Size of the array
C 
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Implied-DO construct could have a significant memory
C      overhead in naively implemented compiler and is better
C      avoided if large amounts of data have to be transferred.
C      At the same time, Fortran does not allow array I/O for
C      arrays of unspecified size, so we have to pass it to
C      a very simple routine as a dummy argument.
C
C   Notes:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(ISIZE)
C
      WRITE(IUNIT,REC=IREC) A
C
      RETURN
      END
C
      SUBROUTINE PSDRD(IUNIT,IREC,A,ISIZE)
C
C   Read variable-sized array from direct access file
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IUNIT  - Unit to read array from. It should be already
C               connected for unformatted direct access
C      IREC   - Record number
C      A      - Array to read to
C      ISIZE  - Size of the array
C 
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Implied-DO construct could have a significant memory
C      overhead in naively implemented compiler and is better
C      avoided if large amounts of data have to be transferred.
C      At the same time, Fortran does not allow array I/O for
C      arrays of unspecified size, so we have to pass it to
C      a very simple routine as a dummy argument.
C
C   Notes:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(ISIZE)
C
      READ(IUNIT,REC=IREC) A
C
      RETURN
      END
C
      SUBROUTINE PSBDMP
C
C   Produce binary dump of execution statistics for derivatives
C   module.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C   Accessed common blocks:
C
C      PSDOPT, 
C      PSDGBL,
C      PSDGB2,
C      PSDGB3,
C      PSOCC,
C      PSHALF,
C      PSDSTT  - are dumped to provide data point for
C                fitting timing model.
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
      PARAMETER (MAXGRP=10)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
C
      PARAMETER (IPSSIG=452789862)
C
      LOGICAL UHF, HALFEL, DOCI, ISRES, DORESP, DODIP, LIMAG, DOPTCH
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
     ./PSDGB3/ NOPAIR, INPAIR
     ./PSOCC / DOCCA(MAXGRP), DOCCB(MAXGRP), NMGRPA, NMGRPB,
     .         IOPN1A, IOPN1B, IOPNLA, IOPNLB, IMAXBL,
     .         IG1STA(MAXGRP), IGCNTA(MAXGRP),
     .         IG1STB(MAXGRP), IGCNTB(MAXGRP),
     .         IGBASA(MAXGRP,MAXGRP), IGBASB(MAXGRP,MAXGRP)
      COMMON
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN),
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN),
     .         IACTI(MAXACT), IACTJ(MAXACT)
     ./PSDSTT/ STIMES(-1:9), XTIMES(-1:9)
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSDSTT/
      SAVE /PSDGB3/, /PSOCC /, /PSHALF/
C
      IF( IURES.EQ.0 ) RETURN
C
      INQUIRE(UNIT=IURES+1,OPENED=ISRES)
      IF( .NOT.ISRES ) THEN
          OPEN(UNIT=IURES+1,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     .                    FORM='UNFORMATTED')
          REWIND(IURES+1)
      ENDIF
C
      WRITE(IURES+1) IPSSIG, 3
      WRITE(IURES+1)
     .     DSTORE, DSTEP,  DECONV, DPCONV, DCPHF, DPREC, DCDIFF,
     .     DSHIFT, DBASCR, DNCOFF,
     .     IUMIX,  IURHS,  IUK,    IURES,
     .     IPRINT, IENRG,  ICORE,  IDISK,  IMIX,   IDENS,
     .     INDSYM, IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC,
     .     IROWS,  IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP,
     .     IKMODE, ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT
      WRITE(IURES+1)
     .     MODE, NATOM, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .     UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG
      WRITE(IURES+1)
     .     NPAIR, NPAIR2, NOCCA, NVACA, NOCCB, NVACB, 
     .     IQSZA, IQSZB,  IQSZ,  ICPV1, ICPVL, NCPVRS
      WRITE(IURES+1)
     .     NOPAIR, INPAIR
      WRITE(IURES+1)
     .     DOCCA, DOCCB, NMGRPA, NMGRPB,
     .     IOPN1A, IOPN1B, IOPNLA, IOPNLB, IMAXBL,
     .     IG1STA, IGCNTA, IG1STB, IGCNTB, IGBASA, IGBASB
      WRITE(IURES+1)
     .     HLFH, HLFG, NUMOIJ, IOI, IOJ,
     .     IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID,
     .     IACTI, IACTJ
      WRITE(IURES+1)
     .     XTIMES

      RETURN
      END
C
      SUBROUTINE PSTIME(IPOINT)
C
C   Record and optionally report per-stage execution time.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IPOINT - Execution checkpoint to record time for.
C               -2 is used to specify start time
C               -1 is used to specify total time
C
C   Accessed common blocks:
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
      PARAMETER (ZERO=0.D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      CHARACTER*40 STEP
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
     ./PSDSTT/ STIMES(-1:9), XTIMES(-1:9)
     ./PSPRT / NB6
      SAVE /PSDOPT/, /PSDGBL/, /PSDSTT/, /PSPRT /, T1, TSTART
      DIMENSION STEP(-1:9)
C
      DATA STEP( 0) /' SETTING UP TOOK                      '/
      DATA STEP( 1) /' INTEGRALS TOOK                       '/
      DATA STEP( 2) /' RIGHT-HAND SIDES TOOK                '/
      DATA STEP( 3) /' Z-VECTOR RHS TOOK                    '/
      DATA STEP( 4) /' SHIFT TOOK                           '/
      DATA STEP( 5) /' TRANSFORMATION TOOK                  '/
      DATA STEP( 6) /' CONDITIONER TOOK                     '/
      DATA STEP( 7) /' CPHF SOLUTION TOOK                   '/
      DATA STEP( 8) /' Z-VECTOR MO TO AO TOOK               '/
      DATA STEP( 9) /' SUMMATION TOOK                       '/
      DATA STEP(-1) /' TOTAL TIME FOR ANALYTICAL DERIVATIVES'/
C
      IF( IPOINT.LT.-2 .OR. IPOINT.GT.9 ) THEN
          WRITE(NB6,10000) IPOINT
          STOP 'PSTIME'
      ENDIF
C
      IF(IPOINT.EQ.-2) THEN
          CALL CPUSEC(TSTART)
          T1 = TSTART
          DO 100 I=-1,9
          XTIMES(I) = ZERO
  100     CONTINUE
          RETURN
      ENDIF
C
      CALL CPUSEC(T2)
      IF( IPOINT.NE.-1 ) THEN
          XTIMES(IPOINT) = T2 - T1
          T1 = T2
      ELSE
          XTIMES(IPOINT) = T2 - TSTART
      ENDIF
      IF( IPRINT.GE.1 .OR. 
     .        IPOINT.EQ.-1 .AND. IPRINT.GE.0 .AND. MODE.GE.2 ) THEN
          WRITE(NB6,10100) STEP(IPOINT), XTIMES(IPOINT), STIMES(IPOINT)
      ENDIF
C
      RETURN
10000 FORMAT(' INVALID EXECUTION CHECKPOINT ', I5, ' IN PSTIME')
10100 FORMAT(1X,A40,F14.3,' SECONDS (EXPECTED ',F14.3,')')
      END
C
      LOGICAL FUNCTION LHAVPS()
C
C   This function returns a flag indicating
C   that the PS directory is present.
C
      LHAVPS = .TRUE.
      RETURN
      END
