C     ******************************************************************
C
C     Extrapolation of core Hamiltonian and two-center integrals.
C     Numerical evaluation of density matrix derivatives.
C
C     ******************************************************************
      SUBROUTINE PSDXTR(IVAR,IPOFF,A,LDA,DUMP,ROALP,ROBET,
     .                             CALP,CBET,LDC)
C
C   Compute density matrix derivative by numeric differentiation along
C   variable IVAR. Use data stored on the stage 1 to extrapolate H
C   and two-center integrals.
C
C   Coded by: Serge Pachkovsky
C
C   Limitations:
C
C      Two-center integrals are assumed to be in core, either in
C      linear or in matrix form.
C
C   Parameters:
C
C      IVAR   - Number of variable to compute density derivative for.
C      IPOFF  - Number of elements of DP array(s) to skip.
C      A      - Scratch array used by PSCFCA (simplified version
C               of SCFCAL routine)
C      LDA    - Dimension of A array
C      DUMP   - Base array for dynamic data
C      ROALP  - Alpha density matrix in packed format.
C               ROALP (and ROBET in UHF case) is actually located
C               inside A array and is changed by calls to PSCFCA 
C      ROBET  - Beta density matrix in packed format (not used if 
C               UHF flag not set)
C      CALP   - Alpha orbital coefficients. Actually located in
C               A array, and are changed by call to PSCFCA
C      CBET   - Beta orbital coefficients. Not used if UHF flag
C               not set.
C      LDC    - Leading dimension of the CALP, CBET arrays. Although
C               CALP and CBET are used as strictly linear arrays
C               here, we still need it to compute size of 'em.
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. Although some of the
C               options could be modified by strategy routines,
C               they are restored to the original values on return.
C      PSDYNM - Handles for the dynamic memory blocks
C      PSDGBL - Scratch common for global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C      PSEXTR - Used to pass information to the extrapolation 
C               routine
C      PSPRT  - Printing unit
C      PSPRTF - Debug output control flags
C
C   Modified common blocks:
C
C      ATOMC  - Coordinates of atoms. All changed are restored 
C               on exit.
C      FDIAG  - Allow fast diagonalization during first iteration,
C               because density matrix should be very good already.
C               All changes are restored on exit.
C      SCRT   - SCF convergence criteria. All changed are restored
C               on exit.
C
C   Local storage:
C
C      None.
C
C   Module logic:
C
C      PSDXTR changes value of cartesian variable IVAR by DSTEP
C      angstroms and do full SCF computation to get density matrix
C      at new geometry. If symmetric steps are requested (INDSYM=1),
C      computation is repeated with variable changed by -DSTEP
C      starting from density matrix extrapolated from density matrix
C      obtained on previous step. Then, density matrix derivative
C      as computed by finite-difference formula.
C
C      Calls to the SCF procedure are made through alternative
C      interface (PSCFCA instead of SCFCAL), which passes base address
C      of dynamic memory block to integrals extrapolation routine,
C      PSCORE, which computes H matrix and two-center integrals from
C      the values stored during stage 1 computation. If you do not like
C      the idea of calling alternate program, SCFCAL would do almost
C      the same job (extrapolation of integrals is in the first order)
C      wasting more time...
C
C      For derivatives computed with symmetric steps, second derivative
C      is also computed, so that appropriate warning could be issued
C      if selected step is too large. To save on a scratch array, second
C      derivative is computed as (2*(PL-P0) + (PH-PL))/DSTEP**2, there
C      PH, PL and P0 are densities computed during up-step, down-step
C      and at central geometries, respectively. This way, second 
C      derivatives end up in ROALP (and ROBET, in UHF case)
C
C   Bugs:
C
C      Passing aliased arguments (portions of A are aliased by ROALP,
C      ROBET, CALP and CBET) to a function which modifies them strictly
C      speaking is undefined behaivior under ISO/IEC 1539. I choose to
C      to ignore this violation of standard because:
C
C          1. It makes life much easier.
C          2. Probably it does not matter for array arguments of 
C             variable size (at least, I can't imagine sensible Fortran
C             processor where it would matter).
C          3. The rest of the program do it all the time.
C
C      Then LDC is not equal to NORBS, this module assumes that extra
C      components of CALP and CBET are not changed by calls to SCFCAL/
C      PSCFCA. Otherwise, checks on the change of orbital state will
C      break down horribly.
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO = 0.0D0)
      PARAMETER ( ONE = 1.0D0)
      PARAMETER (VONE =-1.0D0)
      PARAMETER ( TWO = 2.0D0)
      PARAMETER (BOHR = 0.529167D0)
      PARAMETER (SMALL= 1.0D-2)
C
      LOGICAL UHF, HALFEL, DOCI, WATCHC, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
C
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
     ./PSEXTR/ STEP, IATOM, INDX, IHIDAT
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./INOPT2/ IN2(300)
     ./SCRT  / SCFCRT,PLCRT,SCFCRS(2)
      SAVE /PSDOPT/, /PSDGBL/, /PSPRTF/, /PSPRT /
C
      DIMENSION A(LDA), DUMP(*), ROALP(*), ROBET(*)
      DIMENSION CALP(*), CBET(*)
C
      IPA    = IPSMOF(LPCNTA)
      IDPA   = IPSMOF(LDPA) + IPOFF
      ICA    = IPSMOF(LCCNTA)
      IF(UHF) THEN
          IPB    = IPSMOF(LPCNTB)
          IDPB   = IPSMOF(LDPB) + IPOFF
          ICB    = IPSMOF(LCCNTB)
      ENDIF
      WATCHC = .FALSE.
      IF( DCDIFF.GT.ZERO .AND. IPRINT.GT.-5 ) WATCHC = .TRUE.
      IATOM  = (IVAR-1)/3+1
      INDX   = MOD(IVAR-1,3)+1
      IHIDAT = IPSMOF(LHIDAT)
      ICSIZE = NORBS*LDC
C     IOLDFT = IFAST
      IOLDFT = IN2(73)
      OLDSCF = SCFCRT
      OLDPL  = PLCRT
      SCFCRT = DECONV
      PLCRT  = DPCONV
C     IF( IFAST  .EQ. 0 ) IFAST   = -1
      IF( IN2(73).EQ. 0 ) IN2(73) = -1
      IF( IPRINT.GE.5 ) THEN
          WRITE(NB6,11020) IVAR, IATOM, INDX
          WRITE(NB6,11030) SCFCRT, PLCRT, IN2(73).EQ.-1
      ENDIF
      OLDX   = COORD(INDX,IATOM)
C
C    Up-step
C
      CALL DCOPY(NROS,ROALP,1,DUMP(IPA),1)
      CALL DCOPY(ICSIZE,CALP,1,DUMP(ICA),1)
      IF(UHF) THEN
          CALL DCOPY(NROS,ROBET,1,DUMP(IPB),1)
          CALL DCOPY(ICSIZE,CBET,1,DUMP(ICB),1)
      ENDIF
      COORD(INDX,IATOM) = OLDX + DSTEP
      STEP = DSTEP/BOHR
      CALL PSCFCA(A,LDA,DUMP)
C     CALL SCFCAL(A,LDA)
      CALL DCOPY(NROS,ROALP,1,DUMP(IDPA),1)
      IF(UHF) CALL DCOPY(NROS,ROBET,1,DUMP(IDPB),1)
      IF( INDSYM.EQ.1 ) THEN
C
C         Down-step
C
          CALL DSCAL(NROS,VONE,ROALP,1)
          CALL DAXPY(NROS, TWO,DUMP(IPA),1,ROALP,1)
          IF(WATCHC) THEN
              CALL DAXPY(ICSIZE,VONE,DUMP(ICA),1,CALP,1)
              DCMAXA = CALP(IDAMAX(ICSIZE,CALP,1))
          ENDIF
          CALL DCOPY(ICSIZE,DUMP(ICA),1,CALP,1)
          IF(UHF) THEN
              CALL DSCAL(NROS,VONE,ROBET,1)
              CALL DAXPY(NROS, TWO,DUMP(IPB),1,ROBET,1)
              IF(WATCHC) THEN
                  CALL DAXPY(ICSIZE,VONE,DUMP(ICB),1,CBET,1)
                  DCMAXB = CBET(IDAMAX(ICSIZE,CBET,1))
                  DCMAX  = MAX(DCMAXA,DCMAXB)
              ENDIF
              CALL DCOPY(ICSIZE,DUMP(ICB),1,CBET,1)
          ELSE
              IF(WATCHC) THEN
                  DCMAX = DCMAXA
              ENDIF
          ENDIF
          IF(WATCHC) THEN
              IF( IPRINT.GT.5 ) THEN
                  WRITE(NB6,10020) 'UPSTEP', IVAR, DCMAX
              ENDIF
              IF( DCMAX.GT.DCDIFF ) THEN
                  WRITE(NB6,10030) 'UPSTEP', IVAR, DCDIFF
              ENDIF
          ENDIF
          COORD(INDX,IATOM) = OLDX - DSTEP
          STEP = -DSTEP/BOHR
          CALL PSCFCA(A,LDA,DUMP)
C         CALL SCFCAL(A,LDA)
          CALL DAXPY(NROS,VONE,ROALP,1,DUMP(IDPA),1)
          CALL DAXPY(NROS,VONE,DUMP(IPA),1,ROALP,1)
          CALL DSCAL(NROS, TWO,ROALP,1)
          CALL DAXPY(NROS, ONE,DUMP(IDPA),1,ROALP,1)
          CALL DSCAL(NROS,(BOHR/DSTEP)**2,ROALP,1)
          CALL DSCAL(NROS,BOHR/(TWO*DSTEP),DUMP(IDPA),1)
          D1MAXA = DUMP(IDPA+IDAMAX(NROS,DUMP(IDPA),1)-1)
          D2MAXA = ROALP(IDAMAX(NROS,ROALP,1))
          IF(UHF) THEN
              CALL DAXPY(NROS,VONE,ROBET,1,DUMP(IDPB),1)
              CALL DAXPY(NROS,VONE,DUMP(IPB),1,ROBET,1)
              CALL DSCAL(NROS, TWO,ROBET,1)
              CALL DAXPY(NROS, ONE,DUMP(IDPB),1,ROBET,1)
              CALL DSCAL(NROS,(BOHR/DSTEP)**2,ROBET,1)
              CALL DSCAL(NROS,BOHR/(TWO*DSTEP),DUMP(IDPB),1)
              D1MAXB = DUMP(IDPB+IDAMAX(NROS,DUMP(IDPB),1)-1)
              D2MAXB = ROBET(IDAMAX(NROS,ROBET,1))
              D1MAX  = MAX(D1MAXA,D1MAXB)
              D2MAX  = MAX(D2MAXA,D2MAXB)
          ELSE
              D1MAX = D1MAXA
              D2MAX = D2MAXA
          ENDIF
          IF( IPRINT.GT.5 ) THEN
              WRITE(NB6,10000) IVAR, D1MAX, D2MAX
          ENDIF
          IF( IPRINT.GE.-1 .AND. D2MAX*DSTEP*BOHR.GE.SMALL*D1MAX ) THEN
              WRITE(NB6,10010) IVAR, DSTEP
          ENDIF
          IF(WATCHC) THEN
              CALL DAXPY(ICSIZE,VONE,DUMP(ICA),1,CALP,1)
              DCMAXA = CALP(IDAMAX(ICSIZE,CALP,1))
              IF(UHF) THEN
                  CALL DAXPY(ICSIZE,VONE,DUMP(ICB),1,CBET,1)
                  DCMAXB = CBET(IDAMAX(ICSIZE,CBET,1))
                  DCMAX = MAX(DCMAXA,DCMAXB)
              ELSE
                  DCMAX = DCMAXA
              ENDIF
              IF( IPRINT.GT.5 ) THEN
                  WRITE(NB6,10020) 'DOWNSTEP', IVAR, DCMAX
              ENDIF
              IF( DCMAX.GT.DCDIFF ) THEN
                  WRITE(NB6,10030) 'DOWNSTEP', IVAR, DCDIFF
              ENDIF
          ENDIF
      ELSE
C
C         No down-step
C
          CALL DAXPY(NROS,VONE,DUMP(IPA),1,DUMP(IDPA),1)
          CALL DSCAL(NROS,BOHR/DSTEP,DUMP(IDPA),1)
          IF(UHF) THEN
              CALL DAXPY(NROS,VONE,DUMP(IPB),1,DUMP(IDPB),1)
              CALL DSCAL(NROS,BOHR/DSTEP,DUMP(IDPB),1)
          ENDIF
          IF(WATCHC) THEN
              CALL DAXPY(ICSIZE,VONE,DUMP(ICA),1,CALP,1)
              DCMAXA = CALP(IDAMAX(ICSIZE,CALP,1))
              IF(UHF) THEN
                  CALL DAXPY(ICSIZE,VONE,DUMP(ICB),1,CBET,1)
                  DCMAXB = CALP(IDAMAX(ICSIZE,CBET,1))
                  DCMAX = MAX(DCMAXA,DCMAXB)
              ELSE
                  DCMAX = DCMAXA
              ENDIF
              IF( IPRINT.GT.5 ) THEN
                  WRITE(NB6,10020) 'UPSTEP', IVAR, DCMAX
              ENDIF
              IF( DCMAX.GT.DCDIFF ) THEN
                  WRITE(NB6,10030) 'UPSTEP', IVAR, DCDIFF
              ENDIF
          ENDIF
      ENDIF
      CALL DCOPY(NROS,DUMP(IPA),1,ROALP,1)
      CALL DCOPY(ICSIZE,DUMP(ICA),1,CALP,1)
      IF(UHF) THEN
          CALL DCOPY(NROS,DUMP(IPB),1,ROBET,1)
          CALL DCOPY(ICSIZE,DUMP(ICB),1,CBET,1)
      ENDIF
      COORD(INDX,IATOM) = OLDX
C     IFAST   = IOLDFT
      IN2(73) = IOLDFT
      SCFCRT  = OLDSCF
      PLCRT   = OLDPL
C
      IF(LPNUME) THEN
          WRITE(NB6,11000) IVAR
          CALL PSDPPM(NORBS,DUMP(IDPA))
          IF(UHF) THEN
              WRITE(NB6,11010) IVAR
              CALL PSDPPM(NORBS,DUMP(IDPB))
          ENDIF
      ENDIF
C
      RETURN
10000 FORMAT(' FOR VARIABLE ', I4, ' MAXIMUM DENSITY DERIVATIVES ARE:'
     .      /' FIRST: ', E14.7, ' SECOND: ', E14.7 )
10010 FORMAT(' HIGH-ORDER CONTRIBUTION TO DENSITY MATRIX DERIVATIVE IS',
     .       ' SIGNIFICANT FOR VARIABLE ', I4
     .      /' STEP OF ', E14.7, ' ANGSTROM IS PROBABLY TOO LARGE' )
10020 FORMAT(' MAXIMUM ORBITAL COEFFICIENT CHANGE DURING ', A8, ' ON ',
     .       'VARIABLE ', I4, ' WAS ', E14.7 )
10030 FORMAT(' MAXIMUM CHANGE OF ORBITAL COEFFICIENT DURING ', A8, 
     .       ' ON VARIABLE ', I4
     .      /' EXCEEDS ', E14.7, '. CHANGE OF ELECTRONIC STATE ',
     .       'SUSPECTED.' )
11000 FORMAT(' DERIVATIVE OF ALPHA DENSITY MATRIX ON VARIABLE ', I4, 
     .       ':')
11010 FORMAT(' DERIVATIVE OF BETA DENSITY MATRIX ON VARIABLE ', I4, 
     .       ':')
11020 FORMAT(' COMPUTING NUMERICAL DENSITY DERIVATIVE ON VARIABLE ',
     .       I4, ' (ATOM ', I3, ', VAR ', I1, ')' )
11030 FORMAT(' SCF CONVERGENCE CRITERIA ARE: ', G14.7, ' AND ',
     .       G14.7, '. FAST DIAGONALISATION ALLOWED: ', L1 )
      END
C
      SUBROUTINE PSCORE (H,W,LDH,LDW,DUMP)
C
C   Extrapolate H matrix and two-electron integrals after small
C   change in nuclear coordinates (Coordinate index and change in
C   coordinate are passed in PSEXTR common)
C
C   Coded by: Serge Pachkovsky
C
C   Limitations:
C
C      Two-center integrals are assumed to be in core, either in
C      linear or in matrix form.
C
C      Steps should be along one of axes.
C
C   Parameters:
C
C      H      - Hamiltonian
C      W      - Two-electron integrals. If in matrix storage
C               form, one-center integrals are NOT computed.
C      LDH    - Dimension of H matrix
C      LDW    - Dimension of W matrix
C      DUMP   - Base array for dynamic data
C
C   Accessed common blocks:
C
C      PSDGBL - Scratch common for global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C      PSEXTR - Used to pass information from differentiator
C               routine
C      PSPRT  - Printing unit
C      PSPRTF - Debug output tuning flags
C      FLAG4  - Required to determine integrals storage mode
C      ATOMS  - Orbital numbers
C      PAROPT,
C      DPARM1 - USS, UPP, UDD
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      2277 DOUBLE PRECISION cells.
C
C   Module logic:
C
C      One-center two-electon integrals are just skipped (they are
C      assumed to be computed by previous SCF)
C      Two-center two-electron integrals are either extrapolated
C      from central values using first derivatives (if position of
C      one of the atoms had changed) or just copied.
C      H matrix is completely recomputed from stored and extrapolated
C      values.
C
C   Speedups possible:
C
C      Inlining calls to DCOPY and DAXPY might help.
C
C   Bugs:
C
C      I assume that CPU cache is big enough to hold entire ATEMP,
C      otherwise copying of integrals in block mode will exhibit
C      miserable performance.
C
C      There is no provision for catching too large steps which
C      would break extrapolation accuracy.
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (ONE =1.0D0)
      LOGICAL UHF, HALFEL, DOCI, MATRIX, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
C
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSDGB2/ NPAIR, NPAIR2, NOCCA, NVACA, NOCCB, NVACB,
     .         IQSZA, IQSZB,  IQSZ,  ICPV1, ICPVL, NCPVRS
     ./PSEXTR/ STEP, IATOM, INDX, IHIDAT
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./INOPT2/ IN2(300)
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ),BETAP(LMZ),
     .         ALP(LMZ)
     ./DPARM1/ UDD(LMZ),ZD(LMZ),BETAD(LMZ)
      SAVE /PSDGBL/, /PSDGB2/, /PSPRTF/, /PSPRT /
C
      DIMENSION H(LDH),W(LDW),DUMP(*)
C
      DIMENSION ATEMP(45*45+45+45+9*9)
      DIMENSION AH(9*9)
C
C     MATRIX = IMODE.LE.0
      MATRIX = IN2(211).LE.0
      IIN    = IHIDAT
      IOUT   = 1
      IPAIRA = 0
      DO 1000 NA=1,NATOM
          NFRSTA = NFIRST(NA)
          NORBA  = NLAST(NA) - NFIRST(NA) + 1
          NPAIRA = ( NORBA * (NORBA+1) ) / 2
          NTYPA  = NAT(NA)
          IPAIRB = 0
          AH(1)  = USS(NTYPA)
          IF( NORBA.GT.1 ) THEN
              U      = UPP(NTYPA)
              AH( 2) = ZERO
              AH( 3) = U
              AH( 4) = ZERO
              AH( 5) = ZERO
              AH( 6) = U
              AH( 7) = ZERO
              AH( 8) = ZERO
              AH( 9) = ZERO
              AH(10) = U
              IF( NORBA.GT.4 ) THEN
                  U    = UDD(NTYPA)
                  IH   = 10
                  DO 50 I=4,8
                      DO 48 J=1,I
                          AH(IH+J) = ZERO
   48                 CONTINUE
                      IH     = IH + I + 1
                      AH(IH) = U
   50             CONTINUE
              ENDIF
          ENDIF 
          DO 800 NB=1,NA-1
              NFRSTB = NFIRST(NB)
              NORBB  = NLAST(NB) - NFIRST(NB) + 1
              NPAIRB = ( NORBB * (NORBB+1) ) / 2
              NINTS  = NPAIRA*NPAIRB
              NOVER  = NORBA*NORBB
              NBLOCK = NINTS + NOVER + NPAIRA + NPAIRB
C
              CALL DCOPY(NBLOCK,DUMP(IIN),1,ATEMP,1)
              IF( NA.EQ.IATOM .OR. NB.EQ.IATOM ) THEN
                  IIN = IIN + INDX * NBLOCK
                  IF( NB.EQ.IATOM ) THEN
                      CALL DAXPY(NBLOCK,-STEP,DUMP(IIN),1,ATEMP,1)
                  ELSE
                      CALL DAXPY(NBLOCK, STEP,DUMP(IIN),1,ATEMP,1)
                  ENDIF
                  IIN = IIN + (4-INDX) * NBLOCK
              ELSE
                  IIN = IIN + 4 * NBLOCK
              ENDIF
C
C             All requisite integrals are computed and stored in ATEMP
C             C00L :) Now store two-electron integrals to W array
C
              IF(MATRIX) THEN
                  IW = IPAIRB*NPAIR + IPAIRA
                  IA = 0
                  DO 100 J=1,NPAIRB
                      DO 98 I=1,NPAIRA
                          W(IW+I) = ATEMP(IA+I)
   98                 CONTINUE
                      IW = IW + NPAIR
                      IA = IA + NPAIRA
  100             CONTINUE
                  IW = IPAIRA*NPAIR + IPAIRB
                  DO 104 I=1,NPAIRA
                      IA = I
                      DO 102 J=1,NPAIRB
                          W(IW+J) = ATEMP(IA)
                          IA = IA + NPAIRA
  102                 CONTINUE
                      IW = IW + NPAIR
  104             CONTINUE
              ELSE
                  CALL DCOPY(NINTS,ATEMP,1,W(IOUT),1)
                  IOUT = IOUT + NINTS
              ENDIF
C
C             Store off-diagonal block to the H matrix
C
              IH = (NFRSTA*(NFRSTA+1))/2 - NFRSTA + NFRSTB - 1
              IA = NINTS + NPAIRA + NPAIRB 
              DO 200 I=1,NORBA
                  IAT = IA + I
                  DO 198 J=1,NORBB
                      H( IH+J ) = ATEMP( IAT )
                      IAT       = IAT + NORBA
  198             CONTINUE
                  IH = IH + NFRSTA + I - 1
  200         CONTINUE
C
C             Add core contributions to the diagonal block of H matrix
C             at the atom B, accumulate contributions to the diagonal
C             block on A
C
              IH = (NFRSTB*(NFRSTB+1))/2-1
              IA = NINTS + NPAIRA
              DO 300 I=1,NORBB
                  DO 298 J=1,I
                      H(IH+J) = H(IH+J) + ATEMP(IA+J)
  298             CONTINUE
                  IH = IH+NFRSTB+I-1
                  IA = IA + I
  300         CONTINUE
              CALL DAXPY(NPAIRA,ONE,ATEMP(NINTS+1),1,AH,1)
C
              IPAIRB = IPAIRB + NPAIRB
  800     CONTINUE
C
C         Fill diagonal block of H at the atom A
C
          IH = (NFRSTA*(NFRSTA+1))/2-1
          IA = 0
          DO 900 I=1,NORBA
              DO 898 J=1,I
                  H(IH+J) = AH(IA+J)
  898         CONTINUE
              IH = IH+NFRSTA+I-1
              IA = IA + I
  900     CONTINUE
C
          IPAIRA = IPAIRA + NPAIRA
 1000 CONTINUE
C
      IF(LPNUME) THEN
          WRITE(NB6,11000) 
          CALL PSDPPM(NORBS,H)
          IF(MATRIX) THEN
              WRITE(NB6,11100)
              CALL PSDPGM(NPAIR,NPAIR,W,NPAIR)
          ELSE
              WRITE(NB6,11200)
              IW = 0
              DO 2000 NA=1,NATOM
                  NFRSTA = NFIRST(NA)
                  NORBA  = NLAST(NA) - NFIRST(NA) + 1
                  NPAIRA = ( NORBA * (NORBA+1) ) / 2
                  DO 1900 NB=1,NA-1
                      NFRSTB = NFIRST(NB)
                      NORBB  = NLAST(NB) - NFIRST(NB) + 1
                      NPAIRB = ( NORBB * (NORBB+1) ) / 2
                      WRITE(NB6,11300) NA, NB
                      CALL PSDPGM(NPAIRA,NPAIRB,W(IW),NPAIRA)
 1900             CONTINUE
                  IW = IW + NPAIRA*NPAIRB
 2000         CONTINUE
          ENDIF
      ENDIF
C
      RETURN
11000 FORMAT(' EXTRAPOLATED MATRIX H IS:'/)
11100 FORMAT(' EXTRAPOLATED TWO-ELECTRON INTEGRALS (MATRIX) ARE:'/)
11200 FORMAT(' EXTRAPOLATED TWO-ELECTRON INTEGRALS (ARRAY) ARE:')
11300 FORMAT(' ATOM PAIR ', I4, ',', I4, ':' )
      END
