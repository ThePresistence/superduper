C
C     Driver routines for analytical gradients of half-electron and CI
C     wavefunctions, where fully analytical gradients of the integrals
C     are not available.
C
C     ******************************************************************
      SUBROUTINE PSOMX(IMOD,DERIV,A,LDA,ROALP,ROBET,CALP,CBET,LDC,
     .                 EALP,EBET,ONEDEN,TWODEN)
C
C   Compute cartesian energy gradients for half-electron and CI 
C   wavefunctions, using numerical differentiation for the derivatives
C   of the integrals, and analytical Z-vector formalism for the 
C   wavefunction response.
C
C   Coded by: Serguei Patchkovskii
C
C   Limitations:
C
C      No support for higher derivatives or imaginary perturbations.
C      No support of I/O of intermediate quantities - only in-core
C      solution is implemented.
C
C   Parameters:
C
C      IMOD   - Computation type
C           1 = Compute first derivatives of energy for all atoms
C         201 = Compute first and second derivatives of energy 
C               for all atoms and point charges
C      DERIV  - Output derivatives, should be at least 3*NATOMS
C               Derivatives are in the atomic units.
C      ROALP  - Alpha density matrix in packed format
C      ROBET  - Beta density matrix in packed format (not used if 
C               UHF flag not set)
C      CALP   - Alpha orbital coefficients.
C      CBET   - Beta orbital coefficients (not used if 
C               UHF flag not set)
C      LDC    - Leading dimension of CALP, CBET matrices
C      EALP   - Energy of alpha molecular orbitals, EV. 
C      EBET   - Energy of beta molecular orbitals, EV. Not used if
C               UHF flag not set
C      ONEDEN - Precalculated one-particle density matrix diagonal.
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
C   Module logic:
C
C      Energy gradients of non-variational wavefunctions, like HE
C      or CI methods, with respect to a parameter a, are given by:
C
C         a      a       a                a    
C        E  = EFF + Tr( H  P ) + 0.5*Tr( G P ) 
C
C                          a                      a
C                 + Tr( F P  ) + Tr( (mu nu|la si) \Gamma )
C
C      where F is the Fock matrix, F^a is its static derivative 
C      with respect to the parameter a, P is the density matrix,
C      P^a is the first-order density matrix, (mu nu|la si)^a is
C      a derivative of a two-electron repulsion integral, 
C      \Gamma is the HE/CI part of the two-particle density
C      matrix. EFF is the core repulsion energy (i.e. closed-form
C      energy part, which does not involve density matrix), 
C      H is the one-electron Hamiltonian matrix, and G=(F-H).
C
C      The first-order density matrix can, in principle, be 
C      obtained by solving CPHF equations, with F^a used as
C      a perturbation. In practice, this route is inefficient,
C      so that this part of the derivative is reformulated in 
C      terms of the Z-vector (or the relaxation contribution 
C      to the density matrix, if you prefer) giving:
C
C               a           a    
C        Tr( F P  ) =  Tr( F Z ) 
C
C      The Z-matrix term term is evaluated in PSDST4.
C
C      All remaining terms are evaluated in PSOMS1.
C
C      The rest of this routine takes care of building the
C      Z vector and other auxiliary quantities.
C
C   Other notes:
C
C      See PSDRV for more information.
C
C   Quirks:
C
C   Bugs:
C
C      This routine is very inefficient for calculation of the
C      energy gradients of variational SCF wavefunctions.
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
      DIMENSION DERIV(*)
      DIMENSION ROALP(*), ROBET(*), CALP(LDC,*), CBET(LDC,*)
      DIMENSION EALP(*), EBET(*)
      DIMENSION ONEDEN(*), TWODEN(*)
      DIMENSION A(*)
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
      ROPTS(:) = ROPT(:)
      IOPTS(:) = IOPT(:)
C
C   IMOD can only ask for gradients, possibly on the point charges
C   as well.
C
      IF (IMOD.NE.1 .AND. IMOD.NE.201) THEN
        WRITE (NB6,"('PSOMX: IMOD HAD AN INVALID VALUE = ',I6)") IMOD
        STOP 'PSOMX'
      END IF
C
C    OMx gradients are supported for a small subset of computational
C    paths. Before calling the general strategy driver, we shall
C    choose a compatible set of options, overriding user input if 
C    necessary.
C
C    At the moment, we must use numerical differentiation of the 
C    Fock matrix, AO CPHF with a separate solution phase, in-core
C    storage of the CPHF variables, and the AO route for the Z-vector
C    reduction.
C
      IMIX   =  4  ! Numerical differentiation of the Fock matrix
      IDENS  =  7  ! AO CPHF, separate solution
      IQSWAP = -1  ! all in-core
      IHLWRP =  2  ! AO Z-vector route
      IHLST  =  2  ! In-core MO integral contractions
      IF (ANY( (IOPT(:)/=IOPTS(:)) .AND. (IOPTS(:)/=0)) ) THEN
        IF (IPRINT.GE.0) THEN
          WRITE (NB6,"('PSOMX: OVERRIDING USER INPUT FOR ANALYTICAL"//
     .               " GRADIENT OPTIONS')") 
        END IF
      END IF
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
C       Confusion danger!
C       After the previous call, EALP and EBET, as well as occupation 
C       numbers (OCCA and OCCB) may be in a different order. If 
C       DORESP flag is still set, they have to be restored by calls
C       to PSUSRX before returning to the calling module.
C
        IF(DORESP) THEN
          CALL PSSRTM(CALP,CBET,LDC)
        ENDIF
      ENDIF
      IF(DOCI) THEN
C
C       First stage of CI initialization. Everything what can be done
C       without dynamic memory allocation is done here. In particular,
C       sizes of all temporaries should be determined on return.
C
        CALL PSCIN1
      ENDIF
C
C    Call general strategy driver, to select the remaining parameters
C
      CALL PSDSTR(0,LDC,EALP,EBET)
C
C    Zero out derivatives
C
      ENERGY = ZERO
      DERIV(1:NVARS) = ZERO
C
C    Check for early return for single particle
C
      IF( NVARS.LE.3 ) GOTO 999
C
C    Memory allocation stage
C
      CALL PSDADM(IMEMST,LDC,ITOP,IDSK,.FALSE.)
      IF( IMEMST.NE.0 ) THEN
        WRITE(NB6,"('PSOMX: INSUFFICIENT MEMORY OR DISK FOR "//
     .            "ANALYTICAL DERIVATIVES.')")
        STOP 'PSOMX - MEMORY'
      ENDIF
      IF( IPRINT.GE.0) THEN
        IF (ITOP.NE.0) 
     .    WRITE(NB6,"(' MEMORY REQUIRED    :',I10,' WORDS')") ITOP
        IF (IDSK.NE.0) 
     .    WRITE(NB6,"(' DISK SPACE REQUIRED:',I10,' WORDS')") IDSK
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
      IF(HALFEL) CALL PSHFLO(CALP(1,IOPEN1),LDC,NUMOPN,DUMP)
      IF(DOCI)   CALL PSHFLO(DUMP(IPSMOF(LCIORB)),NORBS,NCIO,DUMP)
C
C    Calculate static part of the derivatives
C
      CALL PSOMS1(DERIV,ROALP,ROBET,DUMP,A,LDA)
C
      CALL PSTIME(1)
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
C           How do you like this beauty?
C          (I'd be glad to pass DUMP, but we want to change type of
C             IFI and IEI coefficients from dp to int)
            CALL PSCICK(EALP,DUMP(IPSMOF(LCIFI)),
     .               DUMP(IPSMOF(LCIFC)),DUMP(IPSMOF(LCIEI)),
     .               DUMP(IPSMOF(LCIEC)),DUMP(IPSMOF(LCIVEC)),
     .               DUMP(IPSMOF(LCIH)),DUMP(IPSMOF(LCIVC1)),
     .               DUMP(IPSMOF(LCIINT)))
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
C
        IF( IDENS.GT.3 ) THEN
C
C   Stage 6 - Compute preconditioner vector using the shift
C             vector computed on stages 4 or 5. We don't need
C             this for the direct solver, which works with the
C             non-conditioned system.
C
          CALL PSMSTG(DUMP,6)
          IF( IPRECT.NE.1 ) THEN
            CALL PSDS3P(EALP,EBET,DUMP(IPSMOF(LCOND)),
     .                            DUMP(IPSMOF(LSHIFT)))
          ELSE
C
C              IPRECT=1 is a special case, since it's shift
C              vector is (implicitly) zero.
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
C
C           Iterative solver, either old-fashioned (PSDST3)
C           or hi-tech (PSLINS)
C
            IF( ISOLVE.EQ.1 ) THEN
              CALL PSDST3(CALP,CBET,LDC,EALP,EBET,DUMP)
            ELSE
              CALL PSLINS(CALP,CBET,LDC,EALP,EBET,DUMP)
            ENDIF
          CALL PSTIME(7)
        END IF
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
        CALL PSDST4(DERIV,FORCE,DIPDRV,1,ROALP,ROBET,DUMP,A,LDA,
     .              CALP,CBET,LDC,EALP,EBET)
        CALL CPUSEC(T2)
        CALL PSTIME(9)
C
C      Do a final on the integrity of memory blocks
C
        CALL PSMSTG(DUMP,0)
      ENDIF
C
C    Deallocate scratch array
C
      DEALLOCATE (DUMP)
C
  999 CONTINUE
      IF(DORESP) THEN
C
C       Order of orbital coefficients, occupation numbers and 
C       orbital energies have to be restored now.
C
        CALL PSUSRT(EALP,EBET,CALP,CBET,LDC)
      ENDIF
      IF( IPRINT.GE.1 ) THEN
C
C       Report computed derivatives
C
        IF( NVARS.LT.20 .OR. IPRINT.GE.5 ) THEN
          WRITE(NB6,"(' COMPUTED DERIVATIVES: '/)")
          CALL PSDPDR(NVARS,DERIV)
        ENDIF
      ENDIF
C
      CALL PSTIME(-1)
C
C    Clean up and exit.
C
      ROPT(:) = ROPTS(:)
      IOPT(:) = IOPTS(:)
      RETURN
      END
C
      SUBROUTINE PSOMS1(DERIV,ROALP,ROBET,DUMP,A,LDA)
C
C   Compute static part of the derivatives for methods where
C   analytical gradients of the integrals are not available
C
C   Origin: Modification of PSDST1
C
C   Coded by: S. Patchkovskii
C
C   Parameters:
C
C      DERIV  - Output first derivatives. For RHF and UHF case
C               they are completely determined on stage 1. In
C               half-electron method, additional terms have to
C               be computed later.
C      ROALP  - Alpha density matrix in packed format
C      ROBET  - Beta density matrix in packed format (not used if 
C               UHF flag not set)
C      DUMP   - Scratch array. 
C      A      - General buffer according to MNDO conventions.
C             - Contains input density matrices A(LS8) and A(LU8).
C             - Provides scratch arrays for core Hamiltonian A(LS7),
C             - for Fock matrices A(LS6) and A(LU6), and for
C             - two-center two-electron integrals A(LS9).
C             - This array is passed from SCFCAL.
C      LDA    - Dimension of the array A(LDA) needed by COMPHF
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGBL - Global computation parameters
C      PSDYNM - Dynamic memory handles
C      PSPRT  - Printing unit.
C      PSPRTF - Debug output tuning flags
C      PSHALF - Half-electron parameters
C      ATOMS  - Starting orbital number, number of orbitals
C               and atom types.
C      INOPT2 - Input options.
C      QMMM1  - Point charges and their positions
C      QMMM2  - Indicates "link" atoms, which do not interact with
C               the MM part of the system.
C      LIMITS - Size of 2-e matrix (LM9)
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C     This routine is numerically less inefficient than the standard
C     finite-difference routines for the gradient due to the need to
C     evaluate and store the Fock matrix derivatives and two-electron
C     derivatives.
C
      USE LIMIT, ONLY: LEN, LM1, LMZ, LM1M, LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
      PARAMETER (MAXCIO=LMACT)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ONE=1.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      LOGICAL ORBDER
      LOGICAL DOCHRM, DONAC
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
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      COMMON
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN),
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN),
     .         IACTI(MAXACT), IACTJ(MAXACT)
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
     ./PSPRT / NB6
      SAVE /PSDOPT/, /PSDGBL/, /PSPRTF/, /PSHALF/, /PSDCIP/, /PSPRT /
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./INOPT2/ IN2(300)
     ./QMMM1 / COORDM(3,LM1M),CHARGM(LM1M)
     ./QMMM2 / LINK(LM1)
     ./QMMM6 / ISELCT(LM1+LM1M)
     ./LIMITS/ LM2,LM3,LM4,LM6,LM7,LM8,LM9
C
      DIMENSION DERIV(*), DUMP(*), ROALP(*), ROBET(*)
      DIMENSION A(LDA)
C
      DIMENSION AINTS(45,45), ADERIV(45,45)
      DIMENSION XJ(9,9), XK(9,9)
C
      EXTERNAL TR_PACK
C
C    Local dynamic memory. We need to keep track of:
C
C     1. Fock matrix derivatives (FDERIV)
C        The first index gives orbitals, in the packed format.
C        The second index runs up to 1 (.NOT.UHF) or 2 (UHF).
C     2. Derivatives of the 2-electron integrals (ADERIV)
C
      ALLOCATABLE FDERIV(:,:), TWOINT(:), TWODER(:)
C
      NSPIN = 1
      IF (UHF) NSPIN = 2
      ORBDER = HALFEL.OR.DOCI
C
CWT   DOCHRM : Flag for skipping computation for fixed MM atoms
C     DONAC  : Flag for computation of nonadiabatic couplings
      DOCHRM = IN2(127).GT.0
      DONAC  = IN2(160).GT.10
C
C    Prepare for a half-electron or CI calculation
C
      IF (ORBDER) THEN
        IHLO  = IPSMOF(LHLO)
        IHLOA = IHLO
        IF(HALFEL) THEN
            NOIJ   = NUMOIJ
            NACT   = NUMOPN
        ENDIF
        IF(DOCI) THEN
            NOIJ   = (NCIO*(NCIO+1))/2
            NACT   = NCIO
            ICI2ES = IPSMOF(LCI2ES)
            ICI2TM = IPSMOF(LCI2TM)
            ICIGAM = IPSMOF(LCIGAM)
        ENDIF
        IF( IHLST.NE.2 ) THEN
          WRITE (NB6,"('PSOMS1: IHLST = ',I3,' IS NOT SUPPORTED')") 
     .           IHLST
          STOP 'PSOMS1 - IHLST'
        END IF
      END IF
C
      ALLOCATE (FDERIV(NROS,NSPIN),TWOINT(LM9),TWODER(LM9),STAT=IALLOC)
      IF (IALLOC.NE.0) THEN
        WRITE (NB6,"('PSOMS1: ALLOCATION FAILURE. NROS = ',I8,"//
     .             "' NSPIN = ',I2,' LM9 = ',I8,' IALLOC = ',I6)") 
     .         NROS, NSPIN, LM9
        STOP 'PSOMS1 - MEMORY'
      END IF
C
C    Calculate and store two-electron integrals.
C
      IF (DORESP) THEN
        CALL PSOMFG(1,HALF,0,0,A,LDA,ROALP,ROBET,
     1              FDERIV,NROS,TWOINT,DSCAL)
        IAI2T1 = IPSMOF(LAI2T1)
      END IF
C
      IF (LPINTS) THEN
        WRITE (NB6,"(' SPIN-ALPHA DENSITY MATRIX: ')")
        CALL PSDPPM(NORBS,ROALP)
        IF (UHF) THEN
          WRITE (NB6,"(' SPIN-BETA DENSITY MATRIX: ')")
          CALL PSDPPM(NORBS,ROBET)
        END IF
      END IF
C
C    Loop over all atoms and point charges. We need here:
C
C    1. Derivatives of the one-electron part of the Hamiltonian
C    2. Derivatives of the two-electron integrals
C
      ISTL = 1
      CENTRES: DO ICNTR=1,NATOM+NPTCHG
C
CWT     Loop only over QM atoms if there are no MM atoms (.NOT.DOPTCH)
C       or if nonadiabatic couplings are computed (DONAC)
        IF ((.NOT.DOPTCH.OR.DONAC) .AND. (ICNTR.GT.NATOM)) EXIT CENTRES
C       Skip gradient calculation for fixed MM atoms
        IF (DOCHRM .AND. ISELCT(ICNTR).NE.0) CYCLE CENTRES
C
        IA = 3*(ICNTR-1)
        COMPONENTS: DO IC=1,3
          CALL PSOMFG(1,HALF,ICNTR,IC,A,LDA,ROALP,ROBET,
     .                FDERIV,NROS,TWODER,DSCAL)
C
CWT      Skip SCF derivative evaluation in NAC calculation
         IF (.NOT.DONAC) THEN
          IF (LPINTS) THEN
            WRITE (NB6,"(/' DERIVATIVE OF CORE REPULSION ENERGY, "//
     .                 "ICNTR = ',I6,' IC = ',I3,' G = ',G18.9)") 
     .             ICNTR, IC, DSCAL
            WRITE (NB6,"(/' SPIN-ALPHA H+0.5*G MATRIX DERIVATIVE, "//
     .                 "ICNTR = ',I6,' IC = ',I3)") ICNTR, IC
            CALL PSDPPM(NORBS,FDERIV(:,1))
            IF (UHF) THEN
              WRITE (NB6,"(/' SPIN-BETA H+0.5*G MATRIX DERIVATIVE, "//
     .                   "ICNTR = ',I6,' IC = ',I3)") ICNTR, IC
              CALL PSDPPM(NORBS,FDERIV(:,2))
            END IF
          END IF ! LPINTS
C
C        Take the derivative of the variational component
C        first. This is just the trace of the P F product.
C
          DR = TR_PACK(NROS,ROALP,FDERIV(:,1))
          IF (UHF) THEN
            DR = DR + TR_PACK(NROS,ROBET,FDERIV(:,2))
          ELSE
            DR = DR * 2
          END IF
          DERIV(IA+IC) = DERIV(IA+IC) - DR + DSCAL
C
          IF (LPINTS) THEN
            WRITE (NB6,"('CENTRE = ',I6,' IC = ',I1,"//
     .                 "' VAR. FORCE = ',F18.9)") ICNTR, IC, DR
          END IF
         END IF ! .NOT.DONAC
C
          IF (.NOT.(DORESP.OR.ORBDER) .OR. ICNTR.GT.NATOM) THEN
            CYCLE COMPONENTS
          END IF
C
C        The remainder of the loop is executed only when CPHF equations
C        will need to be solved, and only for the atoms - not the 
C        point charges
C
          NA     = ICNTR
          NFRSTA = NFIRST(NA)
          NORBA  = NLAST(NA) - NFIRST(NA) + 1
          NPAIRA = ( NORBA * (NORBA+1) ) / 2
C
          IHLOB  = IHLO
          IF (ORBDER .AND. IC==1) THEN
            IHLOAW = IHLOA
            IHLOA  = IHLOA + NOIJ*NPAIRA
          END IF
C
          ISTR = 1
          ATOM2: DO NB=1,NA
            IB     = 3*(NB-1)
            NFRSTB = NFIRST(NB)
            NORBB  = NLAST(NB) - NFIRST(NB) + 1
            NPAIRB = ( NORBB * (NORBB+1) ) / 2
C
            IF(ORBDER) THEN
              IHLOBW = IHLOB
              IHLOB  = IHLOB + NOIJ*NPAIRB
            ENDIF
C
C          Extract two-electron integrals from the MNDO-standard matrix
C
            CALL PSOMX2(TWOINT,LM9,ISTL,NPAIRA,ISTR,NPAIRB,AINTS,45)
C
            IF (IC==1) THEN
              IF (LPINTS) THEN
                WRITE(NB6,"(' TWO-CENTER INTEGRALS FOR: ',I6,',',I6)") 
     .                NA, NB
                CALL PSDPGM(NPAIRA,NPAIRB,AINTS,45)
              END IF ! LPINTS
C
C            Store 2E integrals for use by CPHF routines. 
C            This should only be done once per pair!
C
C
              IF (DORESP) THEN
                DO JB=1,NPAIRB
                  DUMP(IAI2T1:IAI2T1+NPAIRA-1) = AINTS(1:NPAIRA,JB)
                  IAI2T1 = IAI2T1 + NPAIRA
                END DO
              END IF
            END IF ! IC==1
C
C          There are no one-centre contributions for static orbital
C          terms.
C
            IF (NA==NB) CYCLE ATOM2
C
            IF (ORBDER) THEN
C
C             Fetch 2-electron integral derivatives for this atom block
C
              CALL PSOMX2(TWODER,LM9,ISTL,NPAIRA,ISTR,NPAIRB,ADERIV,45)
C
              IF (LPINTS) THEN
                WRITE(NB6,"(' TWO-CENTER INTEGRALS FOR: ',I6,',',I6,"//
     .                    "' IC = ',I3)") NA, NB, IC
                CALL PSDPGM(NPAIRA,NPAIRB,ADERIV,45)
              END IF ! LPINTS
C
              DR = 0
C
              IF (HALFEL) THEN
                CALL PSHSJK(NPAIRA,NPAIRB,ADERIV,45,DUMP(IHLOAW),
     .                      DUMP(IHLOBW),XJ,XK,9)
                DO K=1,NUMOPN
                  DO J=1,K-1
                    DR = DR + HLFH(J,K)*XJ(J,K) + HLFG(J,K)*XK(J,K)
                  END DO
                  DR = DR + HLFH(K,K)*XJ(K,K)
                END DO
C
                IF (LPHALF.AND.LPINTS) THEN
                  WRITE(NB6,"(' ATOM PAIR ',I4,',',I4)") NA, NB
                  WRITE(NB6,"('     XJ(',I2,'):')") IC
                  CALL PSDPSU(NUMOPN,XJ(1,1),9)
                  WRITE(NB6,"('     XK(',I2,'):')") IC
                  CALL PSDPSU(NUMOPN,XK(1,1),9)
                END IF ! LPHALF.AND.LPINTS
C
              END IF ! HALFEL
C
              IF (DOCI) THEN
                CALL PSZRMB(NPAIRA,NOIJ,DUMP(ICI2TM),NPAIRA)
                CALL PSZRMB(NOIJ,NOIJ,DUMP(ICI2ES),NOIJ)
                CALL DGEMM('N','N',NPAIRA,NOIJ,NPAIRB,ONE,
     .                     ADERIV,45,DUMP(IHLOBW),NPAIRB,
     .                     ZERO,DUMP(ICI2TM),NPAIRA)
                CALL DSYR2K('U','T',NOIJ,NPAIRA,ONE,DUMP(IHLOAW),
     .                      NPAIRA,DUMP(ICI2TM),NPAIRA,ZERO,
     .                      DUMP(ICI2ES),NOIJ)
                CALL PSU2PS(NOIJ,DUMP(ICI2ES),NOIJ)
                IF(LPCI.AND.LPINTS) THEN
                  WRITE(NB6,"('     IJKL(',I2,'):')") IC
                  CALL PSDPPM(NOIJ,DUMP(ICI2ES))
                ENDIF
                CALL PSSCTN(NACT,DUMP(ICI2ES))
                DR = DDOT(NCIGAM,DUMP(ICI2ES),1,DUMP(ICIGAM),1)
              ENDIF ! DOCI
C
C            Update gradients
C
              DERIV(IA+IC) = DERIV(IA+IC) - DR
              DERIV(IB+IC) = DERIV(IB+IC) + DR
C
              IF (LPINTS) THEN
                WRITE (NB6,"(' ATOMS = ',I6,', ',I6,' IC = ',I1,"//
     .                     "' ORB. FORCE = ',F18.9)") NA, NB, IC, DR
              END IF ! LPINTS
C
            END IF ! ORBDER
C
C          End of non-variational W.F. terms
C
            ISTR = ISTR + NPAIRB
          END DO ATOM2
C
C        End of CPHF-related terms
C
        END DO COMPONENTS
C
        ISTL = ISTL + NPAIRA
      END DO CENTRES
C
C    Release local memory
C
      DEALLOCATE (FDERIV,TWOINT,TWODER)
C
      RETURN
C
C    Debugging formats
C
      END
C
      SUBROUTINE PSOMFG(IMODE,GWGT,ICNTR,IC,A,LDA,ROALP,ROBET,
     .                  FDERIV,LDF,ADERIV,DSCAL)
C
C   Compute static derivatives of the Fock matrix and, optionally,
C   of the two-electron integrals, using numerical differentiation.
C
C   Coded by: S. Patchkovskii
C
C   Parameters:
C
C      IMODE  - Input:  0 = calculate Fock matrix derivative
C                       1 = ... also 2-electron integral derivative
C      GWGT   - Input:  Fraction of the 2-electron integral matrix,
C                       to be added to the one-electron terms.
C                         0.0 will give one-electron Hamiltonian matrix.
C                         1.0 will give the Fock matrix.
C                       GWGT has no effect for ICNTR corresponding to
C                       a point charge.
C      ICNTR  - Input:  Index of the atom or point charge.
C                             0 means unperturbed geometry
C                             1 to NATOM are atoms
C                       NATOM+1 to NATOM+NPTCHG-1 are point charges
C      IC     - Input:  Cartesian component of the derivative
C      A      - Input:  General buffer according to MNDO conventions,
C                       contains density matrices at correct location.
C      LDA    - Input:  Required dimension of the buffer A(LDA) in
C                       COMPHF (LDA=LM5).
C      ROALP  - Input:  Alpha density matrix in packed format
C      ROBET  - Input:  Beta density matrix in packed format (not used 
C                       if UHF flag not set)
C      FDERIV - Output: Derivative of the Fock matrix.
C                       1st dim: packed orbital pair index
C                       2nd dim: spin component
C                       Values are in Hartree
C      LDF    - Input:  Leading dimension of the FDERIV matrix
C      ADERIV - Output: Derivatives of the 2-electron integrals.
C                       Values are in Hartree
C               Not used if IMODE=0
C      DSCAL  - Output: Derivative of the "scalar" contributions
C                       to the total energy, such as core repulsion
C                       and force field corrections.
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options
C      PSDGBL - Global computation parameters
C      PSPRT  - Printing unit
C      ATOMS  - Orbital indices
C      ATOMC  - Atomic coordinates
C      QMMM1  - Point charges and their positions
C
C   Local storage:
C
C   Module logic:
C
C     We do standard two-sided numerical derivative, using
C     DSTEP as out Cartesian displacement.
C
C     For point charges, we can skip 2-electron integral part,
C     and differentiate just the 2-electron integrals instead.
C
C   Bugs:
C
C
      USE LIMIT, ONLY: LM1, LMZ, LM1M, LEN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (EV =27.2100D0)
      PARAMETER (BOHR=0.529167D0)
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
     ./PSPRT / NB6
      SAVE /PSDOPT/, /PSDGBL/, /PSPRT /
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ATOMC / COORD(3,LM1)
     ./QMMM1 / COORDM(3,LM1M),CHARGM(LM1M)
     ./LIMITS/ LM2,LM3,LM4,LM6,LM7,LM8,LM9
     ./LMSCF / LS1,LS2,LS3,LS4,LS5,LS6,LS7,LS8,LS9
     ./LMUHF / LU1,LU2,LU3,LU4,LU5,LU6,LU7,LU8
C
      DIMENSION ROALP(*), ROBET(*)
      DIMENSION FDERIV(LDF,*), ADERIV(*)
      DIMENSION A(LDA)
C
C    0. Sanity
C
      IF (LM4.NE.NROS) THEN
        WRITE (NB6,"('PSOMFG: SANITY CHECK: LM4 = ',I8,' /= ',I8,"//
     .             "' = NROS')") LM4, NROS
        STOP 'PSOMFG - LM4'
      END IF
      IF (DSTEP.LE.0) THEN
        WRITE (NB6,"('PSOMFG: SANITY CHECK: DSTEP = ',G12.6,' <= 0')") 
     .         DSTEP
        STOP 'PSOMFG - DSTEP'
      END IF
      IF (LM6**2 /= LM9) THEN
        WRITE (NB6,"('PSOMFG: 2-ELECTRON INTEGRALS MATRIX STORAGE IS"//
     .             " INCONSISTENT: LM6**2 = ',I8,' /= ',I8,' = LM9')") 
     .         LM6**2, LM9
        STOP 'PSOMFG - WSTORE'
      END IF
C
C    1. Do we need 2-electron terms ?
C
      ICALL = 0
      IF (ICNTR.LE.NATOM) ICALL = 1
C
C    2. Increment the coordinate
C
      IF (ICNTR.GT.0) CALL STEP(DSTEP)
C
      CALL COMPHF(A,LDA,ICNTR,ESCALM)
      IF (IMODE.NE.0 .AND. ICNTR.LE.0) CALL WSTORE(A(LS9),LM6,1)
C
C       Store Fock matrix and 2-electron integrals
C
      IF (ICALL.EQ.0) THEN
        FDERIV(1:NROS,1) = A(LS7:LS7+LM4-1)
        IF (UHF) THEN
          FDERIV(1:NROS,2) = FDERIV(1:NROS,1)
        END IF
        IF (IMODE.NE.0) THEN
          ADERIV(1:LM9) = 0
        END IF
      ELSE
        FDERIV(1:NROS,1) = A(LS7:LS7+LM4-1) + GWGT*A(LS6:LS6+LM4-1)
        IF (UHF) THEN
          FDERIV(1:NROS,2) = A(LS7:LS7+LM4-1) + GWGT*A(LU6:LU6+LM4-1)
        END IF
        IF (IMODE.NE.0) THEN
          ADERIV(1:LM9) = A(LS9:LS9+LM9-1)
        END IF
      END IF
C
      IF (ICNTR.LE.0) THEN
        SCALE = 1/(EV)
        FDERIV(1:NROS,1) = FDERIV(1:NROS,1) * SCALE
        IF (UHF) THEN
          FDERIV(1:NROS,2) = FDERIV(1:NROS,2) * SCALE
        END IF
        IF (IMODE.NE.0) THEN
          ADERIV(1:LM9) = ADERIV(1:LM9) * SCALE
        END IF
        DSCAL = 0
        RETURN
      END IF
C
C    3. Decrement the coordinate
C
      CALL STEP(-2*DSTEP)
C
      CALL COMPHF(A,LDA,ICNTR,ESCALP)
C
C    4. Calculate the derivatives
C
      SCALE = 1/(2*DSTEP*EV)
C
      IF (ICALL.EQ.0) THEN
        FDERIV(1:NROS,1) = SCALE*(FDERIV(1:NROS,1)-A(LS7:LS7+LM4-1))
        IF (UHF) THEN
          FDERIV(1:NROS,2) = FDERIV(1:NROS,1)
        END IF
        IF (IMODE.NE.0) THEN
          ADERIV(1:LM9) = 0
        END IF
      ELSE
        FDERIV(1:NROS,1) = SCALE*(
     .      FDERIV(1:NROS,1)-A(LS7:LS7+LM4-1)-GWGT*A(LS6:LS6+LM4-1))
        IF (UHF) THEN
          FDERIV(1:NROS,2) = SCALE*(
     .      FDERIV(1:NROS,2)-A(LS7:LS7+LM4-1)-GWGT*A(LU6:LU6+LM4-1))
        END IF
        IF (IMODE.NE.0) THEN
          ADERIV(1:LM9) = SCALE*(ADERIV(1:LM9)-A(LS9:LS9+LM9-1))
        END IF
      END IF
      DSCAL = SCALE*(ESCALP-ESCALM)
C
C    5. Clean up, debug print, and leave.
C
      CALL STEP(DSTEP)

      IF(IPRINT.GT.5) THEN
         WRITE (NB6,"(//' ICNTR: ',I4)") ICNTR
         WRITE (NB6,"(//' PSOMFG: FOCK MATRIX DERIVATIVES')")
         WRITE (NB6,"(1X,10F10.6)") FDERIV(1:NROS,1)
         WRITE (NB6,"(//' PSOMFG: TWO-ELECTRON DERIVATIVES')")
         WRITE (NB6,"(1X,10F10.6)") ADERIV(1:LM9)
      ENDIF

      RETURN
C
      CONTAINS 
C
        SUBROUTINE STEP(DX)
          IF (ICNTR.LE.0 .OR. ICNTR.GT.NATOM+NPTCHG) THEN
            WRITE (NB6,"('PSOMFG: BAD CENTRE INDEX ',I6,'. VALID"//
     .                 " RANGE IS 1 TO ',I6)") ICNTR, NATOM+NPTCHG
            STOP 'PSOMFG - CENTRE'
          END IF
C
          IF (ICNTR.LE.NATOM) THEN
            COORD (IC,ICNTR      ) = COORD (IC,ICNTR      ) + DX*BOHR
          ELSE
            COORDM(IC,ICNTR-NATOM) = COORDM(IC,ICNTR-NATOM) + DX*BOHR
          END IF
        END SUBROUTINE STEP
C
      END SUBROUTINE PSOMFG
C
      DOUBLE PRECISION FUNCTION TR_PACK(LEN,A,B)
C
C   Compute trace of the product of two symmetric matrices in
C   packed format.
C
C   Coded by: S. Patchkovskii
C
C   Parameters:
C
C      LEN    - Input: Size of the arrays A and B
C      A      - Input: Matrix A, in the packed format
C      B      - Input: Matrix B, in the packed format
C
C   Result:
C
C      TR_PACK = Tr(A B)
C
C   Module logic:
C
C     This is -almost- a scalar product of the two vectors.
C
C     The only problem is that off-diagonal elements must be
C     counted twice. In fact, it's better to do it the other
C     way around: count everything twice, then subtract the
C     diagonal part of the product. This way, we can use
C     (presumably highly optimizer) scalar product for 
C     O(LEN) operations, and Fortran code for the remaining
C     O(SQRT(LEN)) operations.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(LEN), B(LEN)
C
      TR_PACK = 2*DOT_PRODUCT(A,B)
C
      IROW = 0
      IPOS = 0
      DIAGONAL: DO 
        IROW = IROW + 1
        IPOS = IPOS + IROW
        IF (IPOS>LEN) EXIT DIAGONAL
C
        TR_PACK = TR_PACK - A(IPOS)*B(IPOS)
      END DO DIAGONAL
C
      END FUNCTION TR_PACK
C
      SUBROUTINE PSOMX2(W,LM9,ISTL,NL,ISTR,NR,AINTS,LDA)
C
C   Extract a diatomic subblock from MNDO 2-electron integral matrix
C
C   Coded by: S. Patchkovskii
C
C   Parameters:
C
C      W      - Input:  2-electron integrals, as returned by
C                       subroutine PSOMFG. Must be in the matrix
C                       form.
C      LM9    - Input:  Total Size of the W array
C      ISTL   - Input:  Starting index for the left-hand atom
C      NL     - Input:  Number of pairs for the left-hand atom
C      ISTR   - Input:  Starting index for the right-hand atom
C      NR     - Input:  Number of pairs for the right-hand atom
C      AINTS  - Output: Diatomic integral block
C      LDA    - Input:  Leading dimension of AINTS
C
C   Accessed common blocks:
C
C      ATOMS  - Atom indices
C      INOPT2 - IN2(211) gives the value of IMODE
C
C   Local storage:
C
C   Module logic:
C
C     We do standard two-sided numerical derivative, using
C     DSTEP as out Cartesian displacement.
C
C     For point charges, we can skip 2-electron integral part,
C     and differentiate just the 2-electron integrals instead.
C
C   Bugs:
C
C     Only IMODE=0 is supported at the moment. I have no clue
C     how to handle one-centre integrals otherwise.
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON
     ./INOPT2/ IN2(300)
     ./PSPRT / NB6
     ./PSDGB2/ NPAIR, NPAIR2, NOCCA, NVACA, NOCCB, NVACB,
     .         IQSZA, IQSZB,  IQSZ,  ICPV1, ICPVL, NCPVRS
C
      DIMENSION W(NPAIR,NPAIR), AINTS(LDA,*)
C
      IF (IN2(211).GT.0) THEN
        WRITE (NB6,"('PSOMX2: IMODE = ',I6,' IS NOT SUPPORTED.')") 
     .         IN2(211)
        STOP 'PSOMX2 - IMODE'
      END IF
C
      IF (LM9/=NPAIR**2) THEN
        WRITE (NB6,"('PSOMX2: NPAIR**2 = ',I8,' /= ',I8,' = LM9')") 
     .         NPAIR**2, LM9
        STOP 'PSOMX2 - LM9'
      END IF
C
      IF (ISTL<=0 .OR. ISTL+NL-1 > NPAIR .OR. 
     .    ISTR<=0 .OR. ISTR+NR-1 > NPAIR ) THEN
        WRITE (NB6,"()") ISTL, NL, ISTR, NR, 1, NPAIR
        STOP 'PSOMX2 - INDEX'
      END IF
C
      AINTS(1:NL,1:NR) = W(ISTL:ISTL+NL-1,ISTR:ISTR+NR-1)
C
      END SUBROUTINE PSOMX2
