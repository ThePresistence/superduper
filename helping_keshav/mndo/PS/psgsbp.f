C     ******************************************************************
C
C     Formerly: Driver routines for analytical derivative code
C     Now: CPHF equation solver for GSBP      
C
C     ******************************************************************
      SUBROUTINE PSGSBP(IMOD,DERIV,FORCE,DIPDRV,LDF,A,LDA,ROALP,
     .                 ROBET,CALP,CBET,LDC,EALP,EBET,CRGDRV)
C   PSGSBP uses parts of the original PSDRV routine to set up and solve
C   the CPHF equation. The Mulliken charge derivatives are returned in
C   the CRGDRV array.
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
     ./ATOMS / NUMAT 
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSPRTF/, /PSDSTT/, /PSHALF/, 
     .     /PSDCIP/, /PSPRT /
C
      DIMENSION DERIV(*), FORCE(LDF,*), A(LDA)
      DIMENSION ROALP(*), ROBET(*), CALP(LDC,*), CBET(LDC,*)
      DIMENSION EALP(*), EBET(*)
      DIMENSION DIPDRV(LDF,3)
C
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: DUMP
      DIMENSION ROPT(IROPCN), ROPTS(IROPCN), IUNITS(4), 
     .          IOPT(IIOPCN), IOPTS(IIOPCN)
      EQUIVALENCE (DSTORE,ROPT(1)),(IUMIX,IUNITS(1)),(IPRINT,IOPT(1))
C
      DIMENSION DUMMY(1)
      DIMENSION CRGDRV(LDF,NUMAT)
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
C    This is a bit of a convoluted call. It all comes from the lack
C    of pointers in Fortran 77...
          CALL PSDST1(ENERGY,DERIV,FORCE,LDF,ROALP,ROBET,CALP,LDC,DUMP)
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
C          IF( HALFEL.OR.DOCI ) THEN
C
C   Stage 3 - Compute right-hand side for the Z-vector response 
C             problem and store it along with othen right-hand
C             CPHF vectors.
C
C              CALL PSMSTG(DUMP,3)
C              CALL PSHLQ(CALP,LDC,EALP,DUMP)
C              IF(IN2(77).EQ.1) THEN
C                How do you like this beauty?
C                (I'd be glad to pass DUMP, but we want to change type of
C                 IFI and IEI coefficients from dp to int)
C                  CALL PSCICK(EALP,DUMP(IPSMOF(LCIFI)),
C     .                     DUMP(IPSMOF(LCIFC)),DUMP(IPSMOF(LCIEI)),
C     .                     DUMP(IPSMOF(LCIEC)),DUMP(IPSMOF(LCIVEC)),
C     .                     DUMP(IPSMOF(LCIH)),DUMP(IPSMOF(LCIVC1)),
C     .                     DUMP(IPSMOF(LCIINT)))
C              ENDIF
C              CALL PSTIME(3)
C          ENDIF
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
C          IF( HALFEL.OR.DOCI ) THEN
C
C   Stage 8 - Construct AO-basis representation of the Z-vector now.
C
C              CALL PSMSTG(DUMP,8)
C              CALL PSHMKZ(CALP,LDC,DUMP)
C              CALL PSTIME(8)
C          ENDIF
C
C   Stage 9 - Summation of the "response" contributions.
C
          CALL PSMSTG(DUMP,9)

C         Bring mixed derivatives back to memory if they were swapped 
C         out during stages 2 and 3. After that, it's Ok to close 
C         scratch file, since we won't need it anymore

          IF( IMIX.EQ.2 ) THEN
              REWIND(IUMIX)
              CALL PSDRS(IUMIX,DUMP(IPSMOF(LDXPA)),NVARS*NROS)
              IF( UHF ) THEN
                  CALL PSDRS(IUMIX,DUMP(IPSMOF(LDXPB)),NVARS*NROS)
              ENDIF
              CLOSE(UNIT=IUMIX)
          ENDIF
C TB      CALL MODIFIED VERSION OF PSDST4 ROUTINE TO CALCULATE
C         MULLIKEN CHARGE DERIVATIVES EXCLUSIVELY
          CALL PSDS4M(DERIV,FORCE,DIPDRV,LDF,ROALP,ROBET,DUMP,A,LDA,
     .                    CALP,CBET,LDC,EALP,EBET,CRGDRV)
          CALL CPUSEC(T2)
          CALL PSTIME(9)
          IF( IDENS.GE.2 .AND. IQSWAP.EQ.2 ) THEN

C             Delete scratch file with CPHF equations solutions
C             before finishing

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
C
C
C
      SUBROUTINE PSDS4M(DERIV,FORCE,DIPDRV,LDF,ROALP,ROBET,DUMP,A,LDA,
     .                  CALP,CBET,LDC,EALP,EBET,CRGDRV)
C
C
C   Module logic:
C   MODIFIED VERSION OF PSDST4.F ROUTINE TO CALCULATE MULLIKEN CHARGE
C   DERIVATIVES (STORED IN CRGDRV) EXCLUSIVELY   
C
      USE LIMIT, ONLY: LM1, LM1M, LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ONE=1.0D0)
      PARAMETER (TWO=2.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
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
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./INDEX / INDX(LMX) 
      SAVE /PSDOPT/, /PSDGBL/, /PSPRTF/, /PSPRT /
C
      COMMON 
     ./INOPT2/ IN2(300)
     ./QMMM6 / ISELCT(LM1+LM1M)
C
      DIMENSION DERIV(*), FORCE(LDF,*), DUMP(*), ROALP(*), ROBET(*)
      DIMENSION CALP(LDC,*), CBET(LDC,*), A(LDA), EALP(*), EBET(*)
      DIMENSION DIPDRV(LDF,3), DIPD(3)
      DIMENSION DUMMY(1)
      DIMENSION CRGDRV(LDF,NUMAT)
C
C    Get dynamic memory offsets
C
      IF( MODE.GE.2 ) THEN
          IDPA = IPSMOF(LDPA)
          IF(UHF) THEN
              IDPB = IPSMOF(LDPB)
          ELSE
              IDPB = IDPA
          ENDIF
      ENDIF
      IF( IMIX.EQ.4 ) THEN
          IDXPAT = IPSMOF(LDXPAT)
          IF(UHF) THEN
            IDXPBT = IDXPAT + NROS
          ELSE
            IDXPBT = IDXPAT
          END IF
          IEA = IDXPAT
          IEB = IDXPBT
      ELSE IF( IMIX.EQ.3 ) THEN
          IDXPDA = IPSMOF(LDXPDA)
          IDXPAT = IPSMOF(LDXPAT)
          IEA    = IDXPAT
          IF(UHF) THEN
              IDXPBT = IPSMOF(LDXPBT)
              IEB    = IDXPBT
          ELSE
              IDXPBT = IDXPAT
          ENDIF
      ELSE
          IDXPA = IPSMOF(LDXPA)
          IF(UHF) IDXPB = IPSMOF(LDXPB)
      ENDIF
      IF( HALFEL.OR.DOCI ) THEN
          IHLZA = IPSMOF(LHLZA)
      ENDIF
C
      IVARTP = NVARS
      IF( MODE.LE.1 ) IVARTP = 1
      DO 1000 IP=1,IVARTP,IDSTRP
          IF( MODE.GE.2 ) THEN
              IPTOP = MIN(IVARTP,IP+IDSTRP-1)
              IPL   = IPTOP - IP + 1
              IF(LPSUMM) THEN
                  WRITE(NB6,11000) IP, IPTOP
              ENDIF
              DO 200 IPS=IP,IPTOP
                  IPOFF = (IPS-IP)*NROS
                  IF( IDENS.EQ.1 ) THEN
                      CALL PSDXTR(IPS,IPOFF,A,LDA,DUMP,ROALP,ROBET,
     .                            CALP,CBET,LDC)
                  ELSE
                      CALL PSCPHD(IPS,IPOFF,DUMP,CALP,CBET,LDC,EALP,
     .                            EBET)
                  ENDIF
                  IF(DODIP) THEN
                      CALL PSDDPR(DIPD,DUMP(IDPA+IPOFF),
     .                                 DUMP(IDPB+IPOFF))
                      DIPDRV(IPS,1) = DIPDRV(IPS,1) + DIPD(1)
                      DIPDRV(IPS,2) = DIPDRV(IPS,2) + DIPD(2)
                      DIPDRV(IPS,3) = DIPDRV(IPS,3) + DIPD(3)
                  ENDIF
                  CALL PSDSCP(NORBS,NROS,DUMP(IDPA+IPOFF))
                  IF(UHF) THEN
                      CALL PSDSCP(NORBS,NROS,DUMP(IDPB+IPOFF))
                  ELSE
                      CALL DSCAL(NROS,TWO,DUMP(IDPA+IPOFF),1)
                  ENDIF
C                  CALL PSDPPM(NORBS,DUMP(IDPA+IPOFF))
                  DO 300 IA=1,NUMAT
                    JA = NFIRST(IA)
                    JB = NLAST(IA)
                    CRGD = 0.D0
                    DO 400 JORB=JA,JB
                      LL = (JORB*(JORB+1)/2) -1
                      CRGD = CRGD - DUMP(IDPA+IPOFF+LL)
                      IF(UHF) THEN
                        CRGD = CRGD - DUMP(IDPB+IPOFF+LL)
                      ENDIF
  400               CONTINUE
                    CRGDRV(IPS,IA) = CRGD
  300             CONTINUE

  200         CONTINUE
          ENDIF
C
  800     CONTINUE
 1000 CONTINUE
C
      RETURN
11000 FORMAT(' COMPUTING DERIVATIVE OF DENSITY MATRIX ON VARIABLE(S) ',
     .        I4, ' - ', I4 )
11100 FORMAT(' CONTRIBUTION TO VARIABLE ', I4, ' IS: ',
     .        F14.7, ' (ALPHA)' )
11105 FORMAT(' CONTRIBUTION TO VARIABLE ', I4, ' IS: ',
     .        F14.7, ' (BETA)' )
11110 FORMAT(' Z-VECTOR CONTRIBUTION TO VAR ', I4, ' IS ', F20.14 )
      END
