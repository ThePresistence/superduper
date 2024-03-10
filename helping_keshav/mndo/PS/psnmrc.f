C     ******************************************************************
C
C     Control routines for computation of NMR chemical shifts.
C
C     ******************************************************************
      SUBROUTINE PSNMR(INC,SIGMA,ROALP,CALP,LDC,EALP)
C
C   Compute NMR dipole shielding tensor at MNDO/SCF level.
C
C   Coded by: Serge Pachkovsky
C
C   Limitations:
C
C      Only closed-shell systems are supported.
C
C   Parameters:
C
C      INC    - Per-atom logical array. True if the corresponding
C               atom is to be included in shielding computation.
C      SIGMA  - Output shielding tensors.
C      ROALP  - Alpha density matrix in packed format
C      CALP   - Alpha orbital coefficients. Not used if IMOD.LE.1
C      LDC    - Leading dimension of CALP matrix
C      EALP   - Energy of alpha molecular orbitals, EV.
C
C   Options: (COMMOND /PSDOPT/)
C
C      See PSDRV for the description of options.
C
C   Accessed common blocks:
C
C      ATOMS  - Atom types table
C      ELEMNT - Element names
C      PSDOPT - Computation options. All modifications are (or at least
C               should be ;-) restored on return.
C      PSDGBL - Global computation parameters (these are modified by
C               called routines)
C      PSDYNM - Dynamic memory handles (modified by called routines)
C      PSPRT  - Printing unit.
C      PSPRTF - Debugging output control flags (modified by called
C               routines).
C
C   Modified common blocks:
C
C      Generally speaking, all PS* common blocks (except for PSDOPT)
C      can potentially be modified.
C
C   Local storage:
C
C      Some 18 DOUBLE PRECISION words plus some 48 INTEGER words.
C
C   Module logic:
C
C      This a specialized version of PSDRV, so that some (hopefully)
C      helpful comments may be gleaned from there.
C
C      All parts specific to the nature of perturbation are in 
C      psnmr.f (i.e., H0B, HA0 and HAB partial derivatives code).
C      The rest of the program is (hopefully) a generic treatment
C      of a purely imaginary perturbation.
C
C   Bugs:
C
C      Execution time estimations do not work at present.
C
      USE LIMIT, ONLY: LM1, LMZ, MXNICS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IIOPCN=24)
      PARAMETER (IROPCN=10)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (SONE  =-1.0D0)
C
      LOGICAL INC, ISRES
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      CHARACTER*2 ELEMNT
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ELEMTS/ ELEMNT(107)
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
     ./PSNICS/ COORDN(3,MXNICS), NNICS
      COMMON
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ), 
     .         IBETA
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSPRTF/, /PSDSTT/, /PSNICS/, 
     .     /PSNPAR/, /PSPRT /
      EXTERNAL PS3JBD
C
      DIMENSION INC(*), SIGMA(3,3,*)
      DIMENSION ROALP(*), CALP(LDC,*), EALP(*)
C
      DIMENSION ROPT(IROPCN), ROPTS(IROPCN), IUNITS(4), 
     .          IOPT(IIOPCN), IOPTS(IIOPCN)
      EQUIVALENCE (DSTORE,ROPT(1)),(IUMIX,IUNITS(1)),(IPRINT,IOPT(1))
C
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: DUMP
      DIMENSION DUMMY(1), SPRNCP(3)
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
      CALL PSDCNS(10)
      IF(DORESP) THEN
C
C          Calling PSHINI here will allow us to handle both implicit
C          and explicit occupation numbers in a uniform manner.
C
           CALL PSHINI(EALP,DUMMY)
C
C          After the previous call, EALP and EBET, as well as occupation
C          numbers (OCCA and OCCB) might be in a different order. If 
C          DORESP flag is still set, they have to be restored by calls
C          to PSUSRX before returning to the calling module.
C
          IF(DORESP) THEN
C             Sort orbitals
              CALL PSSRTM(CALP,DUMMY,LDC)
          ENDIF
      ENDIF
      CALL PSDSTR(0,LDC,EALP,DUMMY)
      CALL PSNIPR
      IF( IPRINT.GT.-1 ) WRITE(NB6,10100)
      IF( NMRLEV.EQ.1 .AND. IPRINT.GT.-1 ) WRITE(NB6,10101)
      IF( NMRLEV.EQ.2 .AND. IPRINT.GT.-1 ) WRITE(NB6,10102)
      IF( NMRLEV.EQ.3 .AND. IPRINT.GT.-1 ) WRITE(NB6,10103)
      IF( NNICS.GT.0  .AND. NMRLEV.NE.3  ) WRITE(NB6,10104)
      IF( INDSYM.EQ.2 .AND. IPRINT.GT.-1 ) WRITE(NB6,10105)
      IF( INDSYM.NE.2 .AND. IPRINT.GT.-1 ) WRITE(NB6,10106)
      IF( NMRLEV.EQ.3 .AND. IPRINT.GT.-1 ) THEN
         IF( MOD(INTCTL,100).EQ.1   ) WRITE(NB6,10121)
         IF( MOD(INTCTL,100).EQ.2   ) WRITE(NB6,10122)
         IF( MOD(INTCTL,200).GE.100 ) WRITE(NB6,10123)
         IF( MOD(INTCTL,20) .GT.10  ) WRITE(NB6,10124) MOD(INTCTL,20)-10
      ENDIF
      IF( IPRINT.GT.-1 ) WRITE(NB6,10100)
C
C    Memory allocation stage
C
      CALL PSDADM(IMEMST,LDC,ITOP,IDSK,.FALSE.)
      IF( IMEMST.NE.0 ) THEN
          WRITE(NB6,11058)
          STOP 'PSNMR'
      ENDIF
      IF( IPRINT.GE.0 ) THEN
          IF (ITOP.NE.0) WRITE(NB6,11553) ITOP
          IF (IDSK.NE.0) WRITE(NB6,11554) IDSK
      ENDIF
C
C    Allocate scratch array and determine starting time
C
      ALLOCATE (DUMP(ITOP))
      CALL PSTIME(0)
C
C    A lot can be skipped if we know that density matrix won't
C    change....
C
      IF(DORESP) THEN
C
C    Stage 1 - Generate and store integrals necessary for CPHF iterative
C              part.
C
          CALL PSMSTG(DUMP,1)
          CALL PSNST1(DUMP)
          CALL PSTIME(1)
C
C    Create direct access file for right-hand sides of CPHF equations.
C    Case of IQSWAP.EQ.3 (immediately consumed CPHF solutions, IDENS=8)
C    can never happen for NMR chemical shifts.
C
          IF( IQSWAP.GE.0 ) THEN
              IF( IPRINT.GE.5 ) THEN
                  WRITE(NB6,13490) NINT(IQSZ*DSTORE)
              ENDIF
              OPEN(UNIT=IURHS,STATUS='SCRATCH',ACCESS='DIRECT',
     .             FORM='UNFORMATTED',RECL=NINT(IQSZ*DSTORE))
          ENDIF
C
C     Stage 2 -  Build right-hand sides of CPHF equations.
C
          CALL PSMSTG(DUMP,2)
          CALL PSNS1A(CALP,LDC,DUMP)
          CALL PSTIME(2)
C
          IF( LPPASS .AND. IQSWAP.EQ.-1 ) THEN
              WRITE(NB6,10024)
              CALL PSDPGM(IQSZ,NCPVRS,DUMP(IPSMOF(LQ)),IQSZ)
          ENDIF
C
C     Stage 3 is dedicated to half-electron derivatives and is never
C     present in the NMR case.
C
          IF( IPRECT.GT.1 .AND. (IPRECT.NE.4 .AND. (IDENS.EQ.4 .OR. 
     .                            IDENS.EQ.5) .OR. IDENS.GT.5) ) THEN
C
C   Stage 4 - Build shift vector unless using the "best" shift which 
C             will be computed during integrals transformation.
C
              CALL PSMSTG(DUMP,4)
              CALL PSDS2S(CALP,DUMMY,LDC,EALP,DUMMY,DUMP)
              CALL PSTIME(4)
          ENDIF
          IF( IDENS.GE.3 .AND. IDENS.LT.6 ) THEN
C
C   Stage 5 - Explicitly construct CPHF K matrix
C
              CALL PSMSTG(DUMP,5)
              IF( IDENS.EQ.5 .AND. IROWS.LT.IQSZ ) THEN
                  OPEN(IUK,STATUS='SCRATCH',ACCESS='SEQUENTIAL',
     .                 FORM='UNFORMATTED')
                  REWIND(IUK)
              ENDIF
              CALL PSDS2P(CALP,DUMMY,LDC,EALP,DUMMY,DUMP)
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
                  CALL PSDS3P(EALP,DUMMY,DUMP(IPSMOF(LCOND)),
     .                              DUMP(IPSMOF(LSHIFT)))
              ELSE
C
C                 IPRECT=1 is a special case, since its shift
C                 vector is (implicitly) zero.
C
                  CALL PSDS3P(EALP,DUMMY,DUMP(IPSMOF(LCOND)),DUMMY)
              ENDIF
              CALL PSTIME(6)
          ENDIF
C
C   Stage 7 - Solution of the CPHF equations in a separate phase.
C             This stage is unconditionally required for NMR
C             chemical shifts.
C
          CALL PSMSTG(DUMP,7)
          IF( IQSWAP.EQ.1 ) THEN
C
C             Bring right-hand CPHF vectors back to memory if they
C             were swapped out and close scratch file - we do not 
C             need it anymore.
C
              IQ = IPSMOF(LQ)
              DO 300 IX=ICPV1,ICPVL
                  CALL PSDRD(IURHS,IX+(1-ICPV1),
     .                   DUMP(IQ+(IX-ICPV1)*IQSZ),IQSZ)
  300         CONTINUE
              CLOSE(UNIT=IURHS)
              IF( LPPASS ) THEN
                  WRITE(NB6,11120)
                  CALL PSDPGM(IQSZ,NCPVRS,DUMP(IQ),IQSZ)
              ENDIF
          ENDIF
          IF( IDENS.EQ.3 ) THEN
C
C         This one is simple: call LAPACK routine to solve
C         CPHF equations, and we have completed.
C
              IF( IQSWAP.GE.2 ) THEN
                  WRITE(NB6,17086)
                  STOP 'PSNMR'
              ENDIF
C
C        We favor iterative CPHF solution, which is formulated
C        for the problem A*X + B = 0. DSPSV solves A*X = B, so
C        we have to multiply right-hand side by -1 in this case.
C
              CALL DSCAL(IQSZ*NCPVRS,SONE,DUMP(IPSMOF(LQ)),1)
              CALL DSPSV('U',IQSZ,NCPVRS,DUMP(IPSMOF(LKS)),
     .               DUMP(IPSMOF(LIDSPS)),DUMP(IPSMOF(LQ)),
     .               IQSZ,INFO)
              IF( INFO.NE.0 ) THEN
                  WRITE(NB6,17455) INFO
                  STOP 'PSDRV'
              ENDIF
          ELSE
C
C         Iterative solver, either old-fascioned (PSDST3)
C         or hi-tech (PSLINS)
C
              IF( ISOLVE.EQ.1 ) THEN
                  CALL PSDST3(CALP,DUMMY,LDC,EALP,DUMMY,DUMP)
              ELSE
                  CALL PSLINS(CALP,DUMMY,LDC,EALP,DUMMY,DUMP)
              ENDIF
              IF( IDENS.EQ.5 .AND. IROWS.LT.IQSZ ) THEN
                  CLOSE(UNIT=IUK)
              ENDIF
          ENDIF
          CALL PSTIME(7)
          IF( LPPASS ) THEN
              IF(IDENS.GE.3.AND.IDENS.LE.7.AND.IQSWAP.LT.2) THEN
C
C             CPHF solutions are in core, print 'em
C
                  WRITE(NB6,11110)
                  CALL PSDPGM(IQSZ,NCPVRS,DUMP(IPSMOF(LQ)),IQSZ)
              ENDIF
          ENDIF
C
C   Stage 8 - Compute first-order density matrix in AO basis and
C             discare CPHF solutions, which we don't need anymore.
C
          CALL PSMSTG(DUMP,8)
          CALL PSNMKD(CALP,LDC,DUMP)
          IF( IDENS.GE.2 .AND. IQSWAP.EQ.2 ) THEN
              CLOSE(UNIT=IURHS)
          ENDIF
C        Prescreen three-center integrals.
          IF( NMRLEV.GT.2 ) THEN
              IF(DORESP) THEN
                  CALL PSNSCR(ROALP,DUMP(IPSMOF(LDPA)))
              ELSE
                  CALL PSNSCR(ROALP,DUMMY)
              ENDIF
          ENDIF
          CALL PSTIME(8)
          IF( LPPASS ) THEN
C
C         Print our imaginary density matrices. Printout is
C         actually transposed!
C
              WRITE(NB6,11115) 'X'
              CALL PSDPPM(NORBS,DUMP(IPSMOF(LDPA)+0*NROS))
              WRITE(NB6,11115) 'Y'
              CALL PSDPPM(NORBS,DUMP(IPSMOF(LDPA)+1*NROS))
              WRITE(NB6,11115) 'Z'
              CALL PSDPPM(NORBS,DUMP(IPSMOF(LDPA)+2*NROS))
          ENDIF
      ENDIF
C
C   Stage 9 - Summation of the "response" contributions.
C
      CALL PSMSTG(DUMP,9)
      DO 900 IAT=1,NATOM
          IF(INC(IAT)) THEN
              CALL PSNSG(IAT,SIGMA(1,1,IAT),ROALP,DUMP)
              IF( IPRINT.GE.1 ) THEN
                  WRITE(NB6,10200) IAT, ELEMNT(NAT(IAT))
                  CALL PSDPGM(3,3,SIGMA(1,1,IAT),3)
              ENDIF
              CALL PSNANL(SIGMA(1,1,IAT),AVE,ANIS,ASYM,SPRNCP)
              IF( IPRINT.GT.0 ) THEN
                  WRITE(NB6,10220) IAT, ELEMNT(NAT(IAT)), 
     .                           AVE, ANIS, ASYM, SPRNCP
              ENDIF
          ENDIF
  900 CONTINUE
      DO 910 IAT=1,NNICS
          CALL PSNSG(IAT+NATOM,SIGMA(1,1,IAT+NATOM),ROALP,DUMP)
          IF( IPRINT.GE.1 ) THEN
              WRITE(NB6,10201) IAT, (COORDN(I,IAT), I=1,3)
              CALL PSDPGM(3,3,SIGMA(1,1,IAT+NATOM),3)
          ENDIF
          CALL PSNANL(SIGMA(1,1,IAT+NATOM),AVE,ANIS,ASYM,SPRNCP)
          IF( IPRINT.GT.0 ) THEN
              WRITE(NB6,10221) IAT, (COORDN(I,IAT), I=1,3), 
     .                       AVE, ANIS, ASYM, SPRNCP
          ENDIF
  910 CONTINUE
C
C     Print summary of NMR results.
C
      IF( IPRINT.GE.-1 ) THEN
          WRITE(NB6,10300)
          DO 920 IAT=1,NATOM
          IF(INC(IAT)) THEN
             CALL PSNANL(SIGMA(1,1,IAT),AVE,ANIS,ASYM,SPRNCP)
             WRITE(NB6,10302) IAT, ELEMNT(NAT(IAT)), 
     .            (COORD(I,IAT), I=1,3), AVE, ANIS, ASYM, SPRNCP
          ENDIF
  920     CONTINUE
          DO 930 IAT=1,NNICS
          CALL PSNANL(SIGMA(1,1,IAT+NATOM),AVE,ANIS,ASYM,SPRNCP)
          WRITE(NB6,10302) IAT, 'NICS',
     .         (COORDN(I,IAT), I=1,3), AVE, ANIS, ASYM, SPRNCP
  930     CONTINUE
          WRITE(NB6,10100)
      ENDIF
      CALL PSTIME(9)
C
C         Do a final check on the integrity of memory blocks.
C
      CALL PSMSTG(DUMP,0)
C
C         Deallocate scratch array
C
      DEALLOCATE (DUMP)

      IF(DORESP) THEN
C
C        Order of orbital coefficients, occupation numbers and 
C        orbital energies have to be restored now.
C
          CALL PSUSRT(EALP,DUMMY,CALP,DUMMY,LDC)
      ENDIF
C
C     Print computation time for NMR section.
      IPRSAV = IPRINT
      IPRINT = -1
      CALL PSTIME(-1)
      IPRINT = IPRSAV
      IF( IPRINT.GE.0 ) WRITE(NB6,10400) XTIMES(-1)
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
          WRITE(IURES) MODE, NATOM, XTIMES(-1)
          WRITE(IURES) (INC(K),((SIGMA(I,J,K),I=1,3),J=1,3),K=1,NVARS)
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
C    IBETA is reset to the default MNDO resonance function on return.
      IBETA = 1
      RETURN
10100 FORMAT()
10101 FORMAT(' CHEMICAL SHIFTS ARE EVALUATED USING ONE-CENTER TERMS.')
10102 FORMAT(' CHEMICAL SHIFTS ARE EVALUATED USING TWO-CENTER TERMS.')
10103 FORMAT(' CHEMICAL SHIFTS ARE EVALUATED USING THREE-CENTER TERMS.')
10104 FORMAT(' WARNING: NICS DOES NOT MAKE MUCH SENSE WITHOUT ',
     .       'THREE-CENTER TERMS.')
10105 FORMAT(' THE HERMITICITY OF THE NMR MATRICES WILL BE USED TO ',
     .       'REDUCE NUMBER OF INTEGRALS.')
10106 FORMAT(' THE HERMITICITY OF THE NMR MATRICES WILL NOT BE USED.')
10121 FORMAT(' THREE-CENTER INTEGRALS ARE COMPUTED WITH INCOMPLETE ',
     .       'GAMMA FUNCTION EXPANSION.')
10122 FORMAT(' THREE-CENTER INTEGRALS ARE COMPUTED WITH MODIFIED ',
     .       'BESSEL FUNCTION EXPANSION.')
10123 FORMAT(' VARIABLE-ACCURACY THREE-CENTER INTEGRALS ARE USED.')
10124 FORMAT(' THREE-CENTER INTEGRALS ARE COMPUTED WITH STO-',I1,
     .       'G EXPANSIONS.')
10200 FORMAT(' SHIELDING TENSOR FOR NUCLEUS ',I4,' (',A2,')',
     .       ' IS (P.P.M.):')
10201 FORMAT(' NICS TENSOR AT CENTER ',I4,' AT ',3(F10.4,1X),
     .       'IS (P.P.M.):')
10220 FORMAT(' FOR NUCLEUS ',I4,' (',A2,') AVERAGE SHIELDING ',
     .       'PARAMETERS ARE (P.P.M.):'
     .      /'    ISOTROPIC SHIELDING:  ',F15.4
     .      /'    SHIELDING ANISOTROPY: ',F15.4
     .      /'    SHIELDING ASYMMETRY:  ',F15.4
     .      /'    PRINCIPAL COMPONENTS: ',3(F15.4,1X))
10221 FORMAT(' FOR CENTER ',I4,' AT ',3(F10.4,1X),
     .       'AVERAGE NICS IS (P.P.M.):'
     .      /'    ISOTROPIC SHIELDING:  ',F15.4
     .      /'    SHIELDING ANISOTROPY: ',F15.4
     .      /'    SHIELDING ASYMMETRY:  ',F15.4
     .      /'    PRINCIPAL COMPONENTS: ',3(F15.4,1X))
10300 FORMAT(///
     .       ' SUMMARY OF NMR RESULTS (P.P.M.) ',
     .     //' ATOM TYPE       COORDINATES (ANGSTROM)     ',
     .       '     ISOTROPIC  SHIELDING  SHIELDING ',
     .       '      PRINCIPAL TENSOR ELEMENTS      ',
     .      /'                 X          Y          Z    ',
     .       '     SHIELDING ANISOTROPY  ASYMMETRY ',
     .       '    SIGMA11    SIGMA22    SIGMA33'/)
10302 FORMAT(1X,I4,1X,A4,3(1X,F10.4),3X,3(1X,F10.3),2X,3(1X,F10.3))
10400 FORMAT(' TIME FOR NMR COMPUTATION',F14.3,' SECONDS')
11553 FORMAT(' MEMORY REQUIRED    :',I10,' WORDS')
11554 FORMAT(' DISK SPACE REQUIRED:',I10,' WORDS')
11058 FORMAT(' NOT ENOUGH DYNAMIC MEMORY OR DISK SPACE TO DO ',
     .       'ANALYTICAL DERIVATIVES.'
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
10024 FORMAT(' RIGHT HAND OF CPHF EQUATIONS FOR VARIABLE (NOT ',
     .       'INCLUDING Z-VECTOR)')
11100 FORMAT(' CPHF K MATRIX (PACKED):'/)
11105 FORMAT(' CPHF SHIFT VECTOR:'/)
11110 FORMAT(' CPHF SOLUTION VECTORS:'/)
11115 FORMAT(' IMAGINARY DENSITY MATRIX DERIVATIVES ON VARIABLE ',A,':')
11120 FORMAT(' CPHF RIGHT-HAND VECTORS AFTER RELOAD:'/)
      END
C
      SUBROUTINE PSNIPR
C
C   Initialize NMR-related atomic parameters.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C   Accessed common blocks:
C
C      PAROPT,
C      DNBND,
C      DPARM,
C      DPARM1 - General atomic parameters
C      PSNPRD - "True" NMR parameters, read-only.
C
C   Modified common blocks:
C
C      PSNPAR - NMR-specific atomic parameters.
C
C   More detailed description of common blocks PSNPRD and PSNPAR:
C
C      PSNPRD - Predefined NMR parameters (read-only).
C      PSNPAR - Copies of the NMR parameters (working arrays).
C
C       A. Control flags
C
C          IDEFN - Element-specific flag for NMR parameters.
C                - .le.0  Standard parameters are used.
C                - .ge.1  NMR-specific parameters were defined.
C                - .ge.2  Resonance SCF parameters were defined.
C                - .ge.3  One-center SCF parameters were defined.
C
C       B. NMR-specific parameters used only in NMR routines.
C          Changing any of these does not affect "normal" properties.
C
C          BETAN - Beta value used for dipole part of the
C                  HA0 derivative, in eV.
C          ZETAN - Orbital exponents used for all NMR integrals
C                  except R_\mu x R_\nu H^0_{\mu\nu} in HA0,
C                  in atomic units.
C          CN    - Contraction coefficients for Slater AOs used
C                  in basis functions expansion.
C          SCALN - Scaling coefficient for shielding tensor.
C          OFFN  - Offset for chemical shift.
C          NN    - Main quantum number of primitive Slater AOs
C                  in basis funcion expansions.
C          MAXN  - Total number of primitive Slater AOs in basis
C                  expansion. Default 1.
C          MAXLA - Largest orbital moment in the basis set expansion
C                  at a given atom.
C          IBETA - Expression for the SCF resonance integrals.
C                  1 = Normal MNDO ( (1/2)( beta1 + beta2 ) S12 )
C                  Other values of IBETA are no longer supported.
C
C          Note  : Definitions refer to PSNPAR, analogous definitions
C                  for XBETAN, XZETAN, XSCALN, and XOFFN in PSNPRD.
C          Note  : The code allows to expand a basis function in terms
C                  of up to 3 Slater AOs. This is not used in practice.
C                  Basis functions are single Slater AOs.
C
C       C. SCF parameters used in the computation of the wave function.
C          Changing any of these will affect "normal" properties.
C
C          XZETAS - Orbital exponents for resonance integrals, in au.
C          XBETAS - Resonance integral scaling coefficients, in eV.
C          XUS    - One-center one-electron energies, in eV.
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ONE=1.D0)
      COMMON
     ./DNBND / III(LMZ),IIID(LMZ)
     ./DPARM / LORBS(LMZ)
     ./DPARM1/ UDD(LMZ),ZD(LMZ),BETAD(LMZ)
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ),BETAP(LMZ),
     .         ALP(LMZ)
      COMMON
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZETAN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
     ./PSNPRD/ XBETAN(LMZ,3),XZETAN(LMZ,3),XZETAS(LMZ,3),XBETAS(LMZ,3),
     .         XUS(LMZ,3),XSCALN(0:LMZ),XOFFN(0:LMZ,2),IDEFN(0:LMZ)
      SAVE /PSNPAR/
      EXTERNAL PSOV1C
C     Internal flag for resonance integrals.
      IBETA  = 1
C     NICS scale is a special case.
      IF( IDEFN(0).LE.0 ) THEN
          SCALN(0) = ONE
      ELSE
          SCALN(0) = XSCALN(0)
      ENDIF
C     All the remaining parameters are generic.
      DO 30 I=1,LMZ
      CN   (1,0,I) = ONE
      CN   (1,1,I) = ONE
      CN   (1,2,I) = ONE
      NN   (1,0,I) = III  (I)
      NN   (1,1,I) = III  (I)
      NN   (1,2,I) = IIID (I)
      MAXN (0,I)   = 1
      MAXN (1,I)   = 1
      MAXN (2,I)   = 1
      IF( LORBS(I).LE.1 ) THEN
         MAXLA(I) = 0
      ELSE IF( LORBS(I).LE.4 ) THEN
         MAXLA(I) = 1 
      ELSE
         MAXLA(I) = 2
      ENDIF
      IF( IDEFN(I).LE.0 ) THEN
         BETAN(0,I)   = BETAS(I)
         BETAN(1,I)   = BETAP(I)
         BETAN(2,I)   = BETAD(I)
         ZETAN(1,0,I) = ZS   (I)
         ZETAN(1,1,I) = ZP   (I)
         ZETAN(1,2,I) = ZD   (I)
         SCALN(I)     = ONE
      ELSE
         BETAN(0,I)   = XBETAN(I,1)
         BETAN(1,I)   = XBETAN(I,2)
         BETAN(2,I)   = XBETAN(I,3)
         ZETAN(1,0,I) = XZETAN(I,1)
         ZETAN(1,1,I) = XZETAN(I,2)
         ZETAN(1,2,I) = XZETAN(I,3)
         SCALN(I)     = XSCALN(I)
      ENDIF
      DO 20 L=0,MAXL
      DO 10 J=2,MAXC
      ZETAN(J,L,I) = ZETAN(1,L,I)
      CN   (J,L,I) = CN(1,L,I)
      NN   (J,L,I) = NN(1,L,I)
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSNMKD(CALP,LDC,DUMP)
C
C  Convert CPHF solution for imaginary perturbations into density
C  matrix derivatives.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients
C      LDC    - Leading dimention of CALP
C      DUMP   - Base of dynamic memory items.
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGBL - Global computation parameters
C      PSDGB2 - Some more computation parameters
C      PSDYNM - Dynamic memory handles
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      None.
C
C   Module logic:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/
C
      DIMENSION CALP(LDC,*), DUMP(*), DUMMY(1)
C
      DO 500 IP=ICPV1,ICPVL
C
C    Get the CPHF solution vector
C
          IF( IQSWAP.EQ.2 ) THEN
              IXTMP = IPSMOF(LXTMP)
              CALL PSDRD(IURHS,IP+(1-ICPV1),DUMP(IXTMP),IQSZ)
              IALP  = IXTMP
          ELSE
              IQ    = IPSMOF(LQ)
              IALP  = IQ + (IP-ICPV1)*IQSZ
          ENDIF
          IDPA = IPSMOF(LDPA) + (IP-ICPV1)*NROS
C
C    Check solutions for very large coupling coefficiients
C
          IF( IPRINT.GE.0 ) THEN
              CALL PSXCHB(IP,DUMP(IALP))
          ENDIF
C
C    Convert CPHF solution into derivatives of the density matrices
C
          CALL PSXM2P(CALP,DUMMY,LDC,DUMP(IALP),DUMP(IDPA),DUMMY,DUMP)
C
  500 CONTINUE
C
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION PSNTRC(T)
C
C   Return one-third of the trace of a 3x3 matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      T      - 3x3 matrix, naturally
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (THIRD=.333333333333333333333333333333333333D0)
      DIMENSION T(3,3)
C
      PSNTRC = THIRD*(T(1,1)+T(2,2)+T(3,3))
C
      RETURN
      END
C
      SUBROUTINE PSNSG(IA,SIGMA,ROALP,DUMP)
C
C   Compute NMR shielding tensor given zeroth and first order density
C   matrices.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Atom to compute shielding at
C      SIGMA  - Output shielding tensor times 10^-6, atomic units.
C      ROALP  - Zeroth order density matrix
C      DUMP   - Base of dynamic memory items.
C
C   Accessed common blocks:
C
C      ATOMS  - Atom types
C      PSDOPT - Computation options. 
C      PSDGBL - Global computation parameters
C      PSDYNM - Dynamic memory handles
C      PSNPAR - NMR parameters
C      PSPRT  - Printing unit
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      Some 72 DOUBLE PRECISION cells.
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO= 0.D0)
      PARAMETER (ONE = 1.D0)
      PARAMETER (VLIGHT=137.0359895D0)
      PARAMETER (SCALE=2.D6/(VLIGHT**2))
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON 
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
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
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
     ./PSPRT / NB6
      SAVE /PSDOPT/, /PSDGBL/, /PSNPAR/, /PSPRT /
C
      EXTERNAL PSNTRC
      DOUBLE PRECISION PSNTRC
      DIMENSION ROALP(*), DUMP(*), SIGMA(3,3)
      DIMENSION SD3(3,3), SD1(3,3), SD2(3,3), SD(3,3)
      DIMENSION SP3(3,3), SP1(3,3), SP2(3,3), SP(3,3)
      DIMENSION ST3(3,3), ST1(3,3), ST2(3,3)
C
      IHAB = IPSMOF(LHAB)
      IH0B = IPSMOF(LH0B)
      IF(DORESP) IDPA = IPSMOF(LDPA)
C
      IF( IA.LE.NATOM ) THEN
          ITYPA = NAT(IA)
      ELSE
          ITYPA = 0
      ENDIF
      XSCALE = SCALE * SCALN(ITYPA)
C
C    Compute diamagnetic part of the shielding tensor.
C
      CALL PSNHAD(IA,DUMP(IHAB),NORBS,NORBS,ROALP,NROS)
      IF( IPRINT.GE.1 .AND. IA.LE.NATOM ) WRITE(NB6,10108) IA
      IF( IPRINT.GE.1 .AND. IA.GT.NATOM ) WRITE(NB6,10109) IA-NATOM
      IF( (NMRLEV.GE.1 .AND. IPRINT.GE.1) .OR. NMRLEV.EQ.1 ) THEN
          DO 200 IY=1,3
              DO 190 IX=1,3
                  SD1(IX,IY) = PSNML1(IA,DUMP(IHAB+(3*(IY-1)+(IX-1))*
     .                                (NORBS**2)),NORBS,ROALP)
  190         CONTINUE
  200     CONTINUE
          CALL DSCAL(9,XSCALE,SD1,1)
          IF( IPRINT.GE.2 ) THEN
              WRITE(NB6,10110)
              CALL PSDPGM(3,3,SD1,3)
          ENDIF
          IF( IPRINT.GE.1 ) WRITE(NB6,10111) PSNTRC(SD1)
          CALL DCOPY(9,SD1,1,SD,1)
      ENDIF
      IF( (NMRLEV.GE.2 .AND. IPRINT.GE.1) .OR. NMRLEV.EQ.2 ) THEN
          DO 250 IY=1,3
              DO 240 IX=1,3
                  SD2(IX,IY) = PSNML2(IA,DUMP(IHAB+(3*(IY-1)+(IX-1))*
     .                                (NORBS**2)),NORBS,ROALP)
  240         CONTINUE
  250     CONTINUE
          CALL DSCAL(9,XSCALE,SD2,1)
          IF( IPRINT.GE.2 ) THEN
              WRITE(NB6,10120)
              CALL PSDPGM(3,3,SD2,3)
          ENDIF
          IF( IPRINT.GE.1 ) WRITE(NB6,10121) PSNTRC(SD2), 
     .                                     PSNTRC(SD2) - PSNTRC(SD1)
          CALL DCOPY(9,SD2,1,SD,1)
      ENDIF
      IF( (NMRLEV.GE.3 .AND. IPRINT.GE.1) .OR. NMRLEV.EQ.3 ) THEN
          DO 300 IY=1,9
              CALL PSU2PS(NORBS,DUMP(IHAB+(IY-1)*(NORBS**2)),NORBS)
              CALL PSDSCP(NORBS,NROS,DUMP(IHAB+(IY-1)*(NORBS**2)))
  300     CONTINUE
          CALL PSZRMB(3,3,SD3,3)
          CALL DGEMV('T',NROS,9,ONE,DUMP(IHAB),NORBS**2,ROALP,1,ZERO,
     .               SD3,1)
          CALL DSCAL(9,XSCALE,SD3,1)
          IF( IPRINT.GE.2 ) THEN
              WRITE(NB6,10130)
              CALL PSDPGM(3,3,SD3,3)
          ENDIF
          IF( IPRINT.GE.1 ) WRITE(NB6,10131) PSNTRC(SD3), 
     .                                     PSNTRC(SD3) - PSNTRC(SD2)
          CALL DCOPY(9,SD3,1,SD,1)
      ENDIF
C
C    Compute paramagnetic part of the shielding tensor. Parts dependent
C    on the first-order density should be omitted if DORESP is not set,
C    since the first-order density is zero in this case.
C
      CALL PSNH0B(IA,DUMP(IH0B),NORBS,NORBS,ROALP,DUMP(IDPA),NROS)
      CALL PSNHAP(IA,DUMP(IHAB),DUMP(IH0B),NORBS,NORBS)
      IF( IPRINT.GE.1 .AND. IA.LE.NATOM ) WRITE(NB6,10208) IA
      IF( IPRINT.GE.1 .AND. IA.GT.NATOM ) WRITE(NB6,10209) IA - NATOM
      IF( (NMRLEV.GE.1 .AND. IPRINT.GE.1) .OR. NMRLEV.EQ.1 ) THEN
          DO 400 IY=1,3
              DO 390 IX=1,3
                  SP1(IX,IY) = PSNML1(IA,DUMP(IHAB+(3*(IY-1)+(IX-1))*
     .                                (NORBS**2)),NORBS,ROALP) 
                  IF(DORESP) 
     .            SP1(IX,IY) = SP1(IX,IY)
     .                       + PSNML1(IA,DUMP(IH0B+(IY-1)*(NORBS**2)),
     .                                NORBS,DUMP(IDPA+(IX-1)*NROS))
  390         CONTINUE
  400     CONTINUE
          CALL DSCAL(9,XSCALE,SP1,1)
          IF( IPRINT.GE.2 ) THEN
              WRITE(NB6,10210)
              CALL PSDPGM(3,3,SP1,3)
          ENDIF
          IF( IPRINT.GE.1 ) WRITE(NB6,10211) PSNTRC(SP1)
          CALL DCOPY(9,SP1,1,SP,1)
      ENDIF
      IF( (NMRLEV.GE.2 .AND. IPRINT.GE.1) .OR. NMRLEV.EQ.2 ) THEN
          DO 500 IY=1,3
              DO 490 IX=1,3
                  SP2(IX,IY) = PSNML2(IA,DUMP(IHAB+(3*(IY-1)+(IX-1))*
     .                                (NORBS**2)),NORBS,ROALP) 
                  IF(DORESP)
     .            SP2(IX,IY) = SP2(IX,IY)
     .                       + PSNML2(IA,DUMP(IH0B+(IY-1)*(NORBS**2)),
     .                                NORBS,DUMP(IDPA+(IX-1)*NROS))
  490         CONTINUE
  500     CONTINUE
          CALL DSCAL(9,XSCALE,SP2,1)
          IF( IPRINT.GE.2 ) THEN
              WRITE(NB6,10220)
              CALL PSDPGM(3,3,SP2,3)
          ENDIF
          IF( IPRINT.GE.1 ) WRITE(NB6,10221) PSNTRC(SP2),
     .                                     PSNTRC(SP2) - PSNTRC(SP1)
          CALL DCOPY(9,SP2,1,SP,1)
      ENDIF
      IF( (NMRLEV.GE.3 .AND. IPRINT.GE.1) .OR. NMRLEV.EQ.3 ) THEN
          DO 600 IY=1,3
              CALL PSU2PS(NORBS,DUMP(IH0B+(IY-1)*(NORBS**2)),NORBS)
              CALL PSDSCP(NORBS,NROS,DUMP(IH0B+(IY-1)*(NORBS**2)))
  600     CONTINUE
          DO 700 IY=1,9
              CALL PSU2PS(NORBS,DUMP(IHAB+(IY-1)*(NORBS**2)),NORBS)
              CALL PSDSCP(NORBS,NROS,DUMP(IHAB+(IY-1)*(NORBS**2)))
  700     CONTINUE
          CALL PSZRMB(3,3,SP3,3)
          CALL DGEMV('T',NROS,9,ONE,DUMP(IHAB),NORBS**2,ROALP,1,ZERO,
     .               SP3,1)
          IF(DORESP)
     .    CALL DGEMM('T','N',3,3,NROS, ONE,DUMP(IDPA),NROS,DUMP(IH0B),
     .               NORBS**2,ONE,SP3,3)
          CALL DSCAL(9,XSCALE,SP3,1)
          IF( IPRINT.GE.2 ) THEN
              WRITE(NB6,10230)
              CALL PSDPGM(3,3,SP3,3)
          ENDIF
          IF( IPRINT.GE.1 ) WRITE(NB6,10231) PSNTRC(SP3),
     .                                     PSNTRC(SP3) - PSNTRC(SP2)
          CALL DCOPY(9,SP3,1,SP,1)
      ENDIF
C
C    Compute total shielding tensor and return.
C
      IF( IPRINT.GE.1 .AND. IA.LE.NATOM ) WRITE(NB6,10308) IA
      IF( IPRINT.GE.1 .AND. IA.GT.NATOM ) WRITE(NB6,10309) IA-NATOM
      IF( (NMRLEV.GE.1 .AND.IPRINT.GE.1) .OR. NMRLEV.EQ.1 ) THEN
          CALL DCOPY(9,SD1,1,ST1,1)
          CALL DAXPY(9,ONE,SP1,1,ST1,1)
          IF( IPRINT.GE.2 ) THEN
              WRITE(NB6,10310)
              CALL PSDPGM(3,3,ST1,3)
          ENDIF
          IF( IPRINT.GE.1 ) WRITE(NB6,10311) PSNTRC(ST1)
          CALL DCOPY(9,ST1,1,SIGMA,1)
      ENDIF
      IF( (NMRLEV.GE.2 .AND.IPRINT.GE.1) .OR. NMRLEV.EQ.2 ) THEN
          CALL DCOPY(9,SD2,1,ST2,1)
          CALL DAXPY(9,ONE,SP2,1,ST2,1)
          IF( IPRINT.GE.2 ) THEN
              WRITE(NB6,10320)
              CALL PSDPGM(3,3,ST2,3)
          ENDIF
          IF( IPRINT.GE.1 ) WRITE(NB6,10321) PSNTRC(ST2),
     .                                     PSNTRC(ST2) - PSNTRC(ST1)
          CALL DCOPY(9,ST2,1,SIGMA,1)
      ENDIF
      IF( (NMRLEV.GE.3 .AND.IPRINT.GE.2) .OR. NMRLEV.EQ.3 ) THEN
          CALL DCOPY(9,SD3,1,ST3,1)
          CALL DAXPY(9,ONE,SP3,1,ST3,1)
          IF( IPRINT.GE.2 ) THEN
              WRITE(NB6,10330) IA
              CALL PSDPGM(3,3,ST3,3)
          ENDIF
          IF( IPRINT.GE.1 ) WRITE(NB6,10331) PSNTRC(ST3),
     .                                     PSNTRC(ST3) - PSNTRC(ST2)
          CALL DCOPY(9,ST3,1,SIGMA,1)
      ENDIF
      RETURN
C
10108 FORMAT(/,' DIAMAGNETIC SHIELDING PARAMETERS FOR THE ATOM ',I4)
10109 FORMAT(' DIAMAGNETIC SHIELDING PARAMETERS FOR THE NICS CENTER ',
     .       I4)
10110 FORMAT('       1-CENTER DIAMAGNETIC SHIELDING TENSOR:')
10120 FORMAT('     1,2-CENTER DIAMAGNETIC SHIELDING TENSOR:')
10130 FORMAT('   1,2,3-CENTER DIAMAGNETIC SHIELDING TENSOR:')
10111 FORMAT('       1-CENTER ISOTROPIC DIAMAGNETIC SHIELDING ',F10.3)
10121 FORMAT('     1,2-CENTER ISOTROPIC DIAMAGNETIC SHIELDING ',F10.3,
     .       ' (+ ',F10.3,')')
10131 FORMAT('   1,2,3-CENTER ISOTROPIC DIAMAGNETIC SHIELDING ',F10.3,
     .       ' (+ ',F10.3,')')
C
10208 FORMAT(' PARAMAGNETIC SHIELDING PARAMETERS FOR THE ATOM ',I4)
10209 FORMAT(' PARAMAGNETIC SHIELDING PARAMETERS FOR THE NICS CENTER ',
     .       I4)
10210 FORMAT('       1-CENTER PARAMAGNETIC SHIELDING TENSOR:')
10220 FORMAT('     1,2-CENTER PARAMAGNETIC SHIELDING TENSOR:')
10230 FORMAT('   1,2,3-CENTER PARAMAGNETIC SHIELDING TENSOR:')
10211 FORMAT('       1-CENTER ISOTROPIC PARAMAGNETIC SHIELDING ',F10.3)
10221 FORMAT('     1,2-CENTER ISOTROPIC PARAMAGNETIC SHIELDING ',F10.3,
     .       ' (+ ',F10.3,')')
10231 FORMAT('   1,2,3-CENTER ISOTROPIC PARAMAGNETIC SHIELDING ',F10.3,
     .       ' (+ ',F10.3,')')
C
10308 FORMAT('  SHIELDING PARAMETERS FOR THE ATOM ',I4)
10309 FORMAT('  SHIELDING PARAMETERS FOR THE NICS CENTER ',I4)
10310 FORMAT('       1-CENTER SHIELDING TENSOR:')
10320 FORMAT('     1,2-CENTER SHIELDING TENSOR:')
10330 FORMAT('   1,2,3-CENTER SHIELDING TENSOR:')
10311 FORMAT('       1-CENTER ISOTROPIC SHIELDING ',F10.3)
10321 FORMAT('     1,2-CENTER ISOTROPIC SHIELDING ',F10.3,
     .       ' (+ ',F10.3,')')
10331 FORMAT('   1,2,3-CENTER ISOTROPIC SHIELDING ',F10.3,
     .       ' (+ ',F10.3,')')
      END
C
      SUBROUTINE PSNANL(SIGMA,AVE,ANIS,ASYM,DIAG)
C
C   Compute average first-order shielding parameters from the
C   shielding tensor.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      SIGMA  - Input shielding tensor.
C      AVE    - Output average shielding
C      ANIS   - Output shielding anisotropy
C      ASYM   - Output shielding asymmetry
C      DIAG   - Output principal components of the shielding tensor
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      Some 27 DOUBLE PRECISION cells.
C
C   Module logic:
C
C      The input tensor is symmetrized and diagonalized. Average
C      shielding parameters are then computed from the principal
C      values.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (THIRD=0.3333333333333333333333333333333333333333333D0)
      PARAMETER (SMALL=1.D-8)
C
      DIMENSION SIGMA(3,3)
      DIMENSION SS(6), VECT(3,3), DIAG(3), TEMP(9)
C
C    Extract symmetric part of the shielding tensor
C
      SS(1) = SIGMA(1,1)
      SS(2) = HALF*(SIGMA(1,2)+SIGMA(2,1))
      SS(3) = SIGMA(2,2)
      SS(4) = HALF*(SIGMA(1,3)+SIGMA(3,1))
      SS(5) = HALF*(SIGMA(2,3)+SIGMA(3,2))
      SS(6) = SIGMA(3,3)
C
C    Diagonalize symmetric part and compute average shielding 
C    parameters.
C
      CALL TDIAG(SS,VECT,DIAG,TEMP,6,3,3,0)
C
      AVE  = THIRD*(DIAG(1)+DIAG(2)+DIAG(3))
      ANIS = DIAG(3) - HALF*(DIAG(1)+DIAG(2))
      IF( (DIAG(3)-AVE).GT.SMALL) THEN
          ASYM = (DIAG(2)-DIAG(1))/(DIAG(3)-AVE)
      ELSE
          ASYM = ZERO
      ENDIF
      RETURN
      END
C
