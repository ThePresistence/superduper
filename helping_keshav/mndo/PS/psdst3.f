C     ******************************************************************
C
C     Iterative solution of CPHF equations.
C
C     ******************************************************************
      SUBROUTINE PSDS3P(EALP,EBET,COND,SHIFT)
C
C   Compute postconditioner matrix to use with iterative
C   CPHF solver.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      EALP   - Alpha orbital energies
C      EBET   - Beta orbital energies (not used if UHF flag
C               not set.
C      COND   - Post-conditioner matrix
C      SHIFT  - Shift vector. Used only if IPRECT.GT.1
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C      PSDGB2 - "Response" computation parameters
C      PSPRT  - Printing unit.
C      PSPRTF - Debug output tuning flags
C      PSOCC  - Occupation groups description
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
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSDGB2/ NPAIR, NPAIR2, NOCCA, NVACA, NOCCB, NVACB,
     .         IQSZA, IQSZB,  IQSZ,  ICPV1, ICPVL, NCPVRS
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
     ./PSOCC / DOCCA(MAXGRP), DOCCB(MAXGRP), NMGRPA, NMGRPB,
     .         IOPN1A, IOPN1B, IOPNLA, IOPNLB, IMAXBL,
     .         IG1STA(MAXGRP), IGCNTA(MAXGRP),
     .         IG1STB(MAXGRP), IGCNTB(MAXGRP),
     .         IGBASA(MAXGRP,MAXGRP), IGBASB(MAXGRP,MAXGRP)
      SAVE /PSDGBL/, /PSDGB2/, /PSPRTF/, /PSOCC/, /PSPRT /
C
      DIMENSION EALP(*), EBET(*), COND(*), SHIFT(*)
C
      CALL PSDS3R(EALP,COND,SHIFT,NMGRPA,DOCCA,IG1STA,IGCNTA,
     .            IGBASA,MAXGRP)
      IF(UHF) THEN
          CALL PSDS3R(EBET,COND(1+IQSZA),SHIFT(1+IQSZA),NMGRPB,DOCCB,
     .                IG1STB,IGCNTB,IGBASB,MAXGRP)
      ENDIF
      IF(LPCPHF) THEN
          WRITE(NB6,11010) 
          CALL PSXPRT(COND)
      ENDIF
C
      RETURN
11010 FORMAT(' POSTCONDITIONER SUPERVECTOR: '/)
      END
C
      SUBROUTINE PSDS3R(E,COND,SHIFT,NMGRP,DOCC,IG1ST,IGCNT,IGBAS,LDI)
C
C   Compute single-spin part of the postconditioner matrix for use 
C   with iterative CPHF solver.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      E      - Orbital energies
C      COND   - Single-spin part of the post-conditioner matrix
C      SHIFT  - Single-spin part of the shift vector. Used only 
C               if IPRECT.GT.1
C      NMGRP  - Number of occupation groups
C      DOCC   - Occupation numbers
C      IG1ST  - Indices of the first orbital in the occupation
C               group
C      IGCNT  - Number of orbitals in the occupation groups
C      IGBAS  - Base indices of the orbital pairs in the 
C               supervectors
C      LDI    - Leading dimension of IGBAS
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
C
      DIMENSION E(*), COND(*), SHIFT(*), DOCC(NMGRP), IG1ST(NMGRP)
      DIMENSION IGCNT(NMGRP), IGBAS(LDI,*)
C
      DO 1000 I1=1,NMGRP
          IND1   = IG1ST(I1)
          N1     = IGCNT(I1)
          O1     = DOCC(I1)
          DO 900 I2=I1+1,NMGRP
              IND2   = IG1ST(I2)
              N2     = IGCNT(I2)
              O2     = DOCC(I2)
              IBLOCK = IGBAS(I1,I2)
              CALL PSDS3Q(N2,N1,O2,O1,E(IND2),E(IND1),COND(IBLOCK),
     .                    SHIFT(IBLOCK))
  900     CONTINUE
 1000 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSDS3Q(N2,N1,O2,O1,E2,E1,COND,SHIFT)
C
C   Compute single occupation pair block of the single spin part 
C   of the postconditioner matrix
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N2     - Number of MOs with smaller occupation number
C      N1     - Number of MOs with larger occupation number
C      O2     - Smaller occupation number
C      O1     - Larger occupation number
C      E2     - Orbital energies of the MOs with smaller occupation
C               number
C      E1     - Orbital energies of the MOs with larger occupation
C               number
C      COND   - (Part of) post-conditioner matrix
C      SHIFT  - (Part of) diagonal shift vector, used only for
C               IPRECT > 1
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options
C      PSPRT  - Printing unit.
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
      PARAMETER (EV = 27.21D0)
      PARAMETER (AU = 0.0367511944138184490995957D0)
      PARAMETER (ONE=1.0D0)
      PARAMETER (ZERO=0.0D0)
C
      COMMON
     ./PSDOPT/ DSTORE, DSTEP,  DECONV, DPCONV, DCPHF, DPREC, DCDIFF,
     .         DSHIFT, DBASCR, DNCOFF,
     .         IUMIX,  IURHS,  IUK,    IURES,
     .         IPRINT, IENRG,  ICORE,  IDISK,  IMIX,   IDENS,
     .         INDSYM, IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC,
     .         IROWS,  IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP,
     .         IKMODE, ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT
     ./PSPRT / NB6
      SAVE /PSDOPT/, /PSPRT /
C
      DIMENSION E2(*), E1(*), COND(N2,N1), SHIFT(N2,N1)
C
      IF( IPRECT.LE.1 ) THEN
          EVOCC = EV * (O1 - O2)
          DO 200 IORB1=1,N1
              DO 100 IORB2=1,N2
                  COND(IORB2,IORB1) = EVOCC / (E2(IORB2) - E1(IORB1))
  100         CONTINUE
  200     CONTINUE
      ELSE
          AUOCC = AU / (O1 - O2)
          DO 400 IORB1=1,N1
              DO 300 IORB2=1,N2
                  T = (E2(IORB2) - E1(IORB1))*AUOCC + SHIFT(IORB2,IORB1)
                  IF( T.LE.ZERO ) THEN
                      WRITE(NB6,10100) IORB2, IORB1, T
                      STOP 'PSDS3Q'
                  ENDIF
                  COND(IORB2,IORB1) = ONE / SQRT( T )
  300         CONTINUE
  400     CONTINUE
      ENDIF
C
      RETURN
10100 FORMAT(' SQUARE OF PRECONDITIONER MATRIX ELEMENT (', I4, ',',
     .       I4, ' IS NON-POSITIVE (', G20.14, ')'
     .      /' PROGRAM WILL STOP.' )
      END
C
      SUBROUTINE PSDS3C(COND,X,LDX,ICNT)
C
C   Apply postconditioned matrix COND to solution vector(s) X
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      COND   - Post-conditioner matrix
C      X      - Array of solution vectors
C      LDX    - Leading dimension of X
C      ICNT   - Number of vectors to transform
C
C   Accessed common blocks:
C
C      PSDGB2 - "Response" computation parameters.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Second version of PSDS3C (for ICNT>1) tries to
C      improve cache (or vector registers) use by 
C      splitting update of multiply vectors in strips.
C
C   Bugs:
C
C      Unfortunately, there is no BLAS call to do it.
C      VECLEN should be a multiply of 4.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IVECL=64)
C
      COMMON 
     ./PSDGB2/ NPAIR, NPAIR2, NOCCA, NVACA, NOCCB, NVACB,
     .         IQSZA, IQSZB,  IQSZ,  ICPV1, ICPVL, NCPVRS
      SAVE /PSDGB2/
C
      DIMENSION X(LDX,*), COND(*)
C
      IF( ICNT.EQ.1 ) THEN
          DO 100 I=1,IQSZ
              X(I,1) = X(I,1) * COND(I)
  100     CONTINUE
      ELSE
          IFULL = (IQSZ/IVECL)*IVECL
          DO 400 IBL=1,IFULL,IVECL
              DO 300 IVEC=1,ICNT
                  DO 200 I=IBL,IBL+IVECL-1,4
                      X(I+0,IVEC) = X(I+0,IVEC) * COND(I+0)
                      X(I+1,IVEC) = X(I+1,IVEC) * COND(I+1)
                      X(I+2,IVEC) = X(I+2,IVEC) * COND(I+2)
                      X(I+3,IVEC) = X(I+3,IVEC) * COND(I+3)
  200             CONTINUE
  300         CONTINUE
  400     CONTINUE
          DO 600 IVEC=1,ICNT
              DO 500 I=IFULL+1,IQSZ,1
                  X(I,IVEC) = X(I,IVEC) * COND(I)
  500         CONTINUE
  600     CONTINUE
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSDST3(CALP,CBET,LDC,EALP,EBET,DUMP)
C
C   Solve first-order CPHF equations iteratively. This could
C   solve either all equations (then called through this
C   primary entry point) or single equation (if called through
C   alternative entry point PSDS3O which have one additional
C   parameter, number of variable)
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients
C      CBET   - Beta orbital coefficients (not used if
C               UHF flag not set)
C      LDC    - Leading dimension of CALP and CBET
C      EALP   - Alpha orbital energies
C               Neither EALP not EBET are actually used,
C               but they are required parameters for
C               PSDS2P which is called for IDENS=6 (direct
C               MO-basis CPHF)
C      EBET   - Beta orbital energies
C      DUMP   - Base of dynamic memory items
C
C      IIEQ   - Variable to solve CPHF equations for, duplicated
C               in IEQ after entry. IIEQ cannot be used directly
C               since Fortran does not require short-cirquiting
C               logical expressions, and IIEQ is not defined along
C               the primary entry point.
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGB2 - "Response" computation parameters.
C      PSDYNM - Dynamic memory handles
C      PSPRT  - Printing unit.
C      PSPRTF - Debug output tuning flags
C
C   Modified common blocks:
C
C      PSCPST - Iterative CPHF solution statistics
C
C   Local storage:
C
C      3*LM1 (=1500) INTEGER cells are used for variable
C      scheduling array.
C
C   Module logic:
C
C      This is actually a dispatch routine, actual work is
C      done is PSDS3S. PSDST3 will first set up parameters
C      for "optimistic" solution pass, which assumes all
C      equations in each set of INRHS equations could be
C      solved with IAVEIT iterations. If some of equations
C      remain unsolved, PSDST3 requests extra pass with
C      expected average number of iterations set to the
C      doubled maximum number of iterations for yet unsolved
C      variables. This is repeated untill IMAXIT is reached
C      or all equations are solved.
C
C   Bugs:
C
C      Scheduling array should be dynamically allocated.
C
C      Handling of alternative entry point PSDS3O is a kludge.
C      Unnecessary table scan introduced is unlike to become a
C      bottleneck, though, at it is N**2 at most, with N**3 and
C      N**4 steps in the controlled loop.
C     
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXVRS=3*LM1)
C
      PARAMETER (ZERO=0.0D0)
C
      LOGICAL SINGLE
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
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
     ./PSCPST/ NVECS,  NITERS, NWASTE, NMAXVE
     ./PSAXTM/ AXTIME
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSDOPT/, /PSDGB2/, /PSPRTF/, /PSAXTM/, /PSCPST/, /PSPRT /
C
      DIMENSION CALP(LDC,*), CBET(LDC,*), EALP(*), EBET(*), DUMP(*)
      DIMENSION ISCD(0:MAXVRS)
C
      SINGLE = .FALSE.
      IEQ    = 0
      GOTO 1
          ENTRY PSDS3O(CALP,CBET,LDC,EALP,EBET,DUMP,IIEQ)
          SINGLE = .TRUE.
          IEQ    = IIEQ
    1 CONTINUE
      IF( IPRINT.GE.0 .AND. (.NOT.SINGLE .OR. IEQ.EQ.ICPV1) ) THEN
          WRITE(NB6,10200)
      ENDIF
      IF( NCPVRS.GT.MAXVRS ) THEN
          WRITE(NB6,10000) NCPVRS, MAXVRS
          STOP 'PSDST3'
      ENDIF
      IF( IKRVEC.LT.INRHS ) THEN
          WRITE(NB6,10010) IKRVEC, INRHS
          STOP 'PSDST3'
      ENDIF
C
      IF( .NOT.SINGLE .OR. IEQ.EQ.1 ) THEN
          NVECS  = 0
          NITERS = 0
          NWASTE = 0
          NMAXVE = 0
          AXTIME = ZERO
      ENDIF
C
C    Initialise scheduling vector
C
      IF( SINGLE ) THEN
          DO 98 I=ICPV1,ICPVL
              ISCD(I) = 1
   98     CONTINUE
          ISCD(IEQ) = 0
          INFO      = IEQ - 1
      ELSE
          DO 100 I=ICPV1,ICPVL
              ISCD(I) = 0
  100     CONTINUE
          INFO      = 0
      ENDIF
C
C    Precompute some memory control parameters
C
      ICPITM = IPSMOF(LCPITM)
      IDOT   = ICPITM
      IDOTSZ = INRHS*(IMAXIT+1)*(IMAXIT+3)
      IKR    = IDOT + IDOTSZ
C
C    Do "optimistic" solution pass, which expects that IKRVEC of 
C    temporary vectors generated is enougth to solve INRHS equations
C
      NRHS  = INRHS
C
 1000 CONTINUE
      IF(LPCPHF) THEN
          WRITE(NB6,11100) NRHS
      ENDIF
      CALL PSDS3S(CALP,CBET,LDC,EALP,EBET,DUMP,ISCD,NRHS,IMAXIT,
     .            DUMP(IDOT),IQSZ,DUMP(IKR),INFO)
      IF(LPCPHF) THEN
          WRITE(NB6,11110) INFO, NCPVRS
      ENDIF
      IF( INFO.EQ.NCPVRS ) THEN
C
C         We've got solutions for all desired vectors,
C         report statictics and bail out
C
          IF( IPRINT.GE.0 ) THEN
              WRITE(NB6,10100) NCPVRS, NVECS, NITERS, NWASTE, 
     .                       DBLE(NVECS-NWASTE)/DBLE(NCPVRS), NMAXVE,
     .                       DBLE(NITERS)/DBLE(NCPVRS)
              IF( IPRINT.GE.2 ) WRITE(NB6,10110) AXTIME
          ENDIF
C
C         Save actual number of iterations - this might be interesting
C         for time-estimation routine statistics collection.
C
          IAVEIT = NINT(DBLE(NVECS-NWASTE)/DBLE(NCPVRS))
          IMAXIT = NMAXVE
          RETURN
      ENDIF
      IF( SINGLE ) THEN
          IF( INFO.NE.IEQ ) THEN
              WRITE(NB6,10021) 1-ISCD(IEQ), IEQ, NCPVRS
              STOP 'PSDST3'
          ENDIF
          RETURN
      ENDIF
C
C         Do additional "pessimistic" pass(es), doubling expected
C         number of iterations each time until IMAXIT is reached.
C
          ITERS = 0
          DO 2000 I=ICPV1,ICPVL
              IF( ISCD(I).LT.0 ) ITERS = MAX(ITERS,1-ISCD(I))
 2000     CONTINUE
          IF( ITERS.GE.IMAXIT ) THEN
              WRITE(NB6,10020) ITERS, NCPVRS-INFO, NCPVRS
              STOP 'PSDST3'
          ENDIF
          IF( IPRINT.GE.1 ) THEN
              WRITE(NB6,11000) NCPVRS - INFO, ITERS
          ENDIF
          ITERS = MIN(2*ITERS,IMAXIT)
          NRHS = MIN(INRHS,IKRVEC/(ITERS+1))
          GOTO 1000
C
10000 FORMAT(' STATIC ARRAY OVERFLOWED IN PSDST3. COMPILE-TIME LIMIT: ',
     .       I4, ' PASSED VALUE: ', I5
     .      /' PROGRAM WILL STOP.' )
10010 FORMAT(' CONFLICTING PARAMETERS: IKRVEC (', I5, ') SMALLER THAN',
     .       ' INRHS (', I5, ')',
     .      /' PROGRAM WILL STOP.' )
10020 FORMAT(' FAILED TO ACHIEVE REQUESTED CPHF CONVERGENCE IN ', I5,
     .       ' ITERATIONS'
     .      /' ON ', I4, ' VARIABLES OUT OF ', I4, 
     .      /' PROGRAM WILL STOP.' )
10021 FORMAT(' FAILED TO ACHIEVE REQUESTED CPHF CONVERGENCE IN ', I5,
     .       ' ITERATIONS'
     .      /' ON VARIABLE ', I4, ' (TOTAL NUMBER OF VARIABLES ', I4, 
     .       ')'
     .      /' PROGRAM WILL STOP.' )
10100 FORMAT(' ', I5, ' CPHF EQUATIONS WERE SOLVED.', 
     .      /' ', I7, ' KRYLOV  VECTORS WERE GENERATED IN ', I7, 
     .       ' PASSES.'
     .      /' ', I7, ' VECTORS WERE DISCARDED DUE TO MEMORY ',
     .       'CONSTRAINTS.'
     .      /' AVERAGE NUMBER OF USEFUL VECTORS PER VARIABLE: ', F14.7,
     .      /' MAXIMUM NUMBER OF USEFUL VECTORS PER VARIABLE: ', I7,
     .      /' AVERAGE NUMBER OF ITERATIONS PER VARIABLE: ', F14.7 )
10110 FORMAT(' TIME SPENT EVALUATING THE RESPONSE VECTORS: ', F14.3, 
     .       ' SEC.')
10200 FORMAT(' USING PLAIN ITERATIVE CPHF SOLVER ' )
11100 FORMAT(' DOING CPHF SOLUTION PASS WITH NRHS = ', I5 )
C
11000 FORMAT(' CPHF FOR ', I4, ' VARIABLE(S) FAILED TO CONVERGE WITH',
     .       ' MAXIMUM ITERATION COUNT OF ', I4 
     .      /' RUNNING EXTRA CPHF PASS' )
11110 FORMAT(' ', I5, ' EQUATIONS SOLVED OUT OF ', I5 )
      END
C
      SUBROUTINE PSDS3S(CALP,CBET,LDC,EALP,EBET,DUMP,ISCD,NRHS,MAXIT,
     .                  DOTS,IVSZ,VECS,INFO)
C
C   Run single solution pass over all CPHF right-hand sides yet unsolved
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients
C      CBET   - Beta orbital coefficients (not used if
C               UHF flag not set)
C      LDC    - Leading dimension of CALP and CBET
C      EALP   - Alpha orbital energies
C               Neither EALP not EBET are actually used,
C               but they are required parameters for
C               PSDS2P which is called for IDENS=6 (direct
C               MO-basis CPHF)
C      EBET   - Beta orbital energies
C      DUMP   - Base of dynamic memory items
C      ISCD   - State of CPHF computation for corresponding
C               variable
C      NRHS   - Number of right-hand sides to solve
C               simultaneously. This is NOT equal to INRHS!
C      MAXIT  - Maximum number of iterations allowed. This
C               is equal to IMAXIT and was renamed to avoid
C               name clash
C      DOTS   - Dot products for the vectors from Krylov
C               subspaces generated for each variable
C      IVSZ   - Leading dimension of VECS and RHS. This is
C               equal to IQSZ.
C      VECS   - Temporary space for Krylov vectors generated
C      INFO   - Number of variables which were solved so far
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGB2 - "Response" computation parameters.
C      PSDYNM - Dynamic memory handles
C      PSPRT  - Printing unit.
C      PSPRTF - Debug print controls
C
C   Modified common blocks:
C
C      PSCPST - Iterative CPHF solution statistics.
C
C   Local storage:
C
C      3*LM1 (=1500) DOUBLE PRECISION cells and 3*LM1 INTEGER
C      cells for the list of selected variables and error 
C      estimations.
C
C   Module logic:
C
C      This module solves systems of CPHF linear equations by
C      constructing locally-orthogonalised Krylov subspace of
C      the matrix 1-KA, where KA is conditioned CPHF matrix KS,
C      and solving linear system of the reduced dimention. 
C      System of equations closely follows method suggested by
C      Pople et al, but with additional renormalization of basis
C      vectors to norm 1. Unfortunately, methods using non-
C      orthogonalized or locally-orthogonalized basis vectors
C      (which require half as much operations as globally-
C      orthogonal method of Pople et al do) suffer from numeric 
C      instability for large number of iterations and can't 
C      be used reliably.
C
C      (1-A)*X = B(0)
C
C           n         ~
C      X = Sum alp(i)*B(i)
C          i=0
C
C      ~
C      B(i)     = gam(i)*B(i)
C
C      gam(i)   = 1.0/|B(i)|
C
C                   ~       i           ~
C      B(i+1)   = A*B(i) - Sum bet(j,i)*B(j)
C                          j=0
C
C      bet(j,i) = <B(j)|A|B(i)>
C
C                             n
C      (1-bet(0,0))*alp(0) - Sum alp(i)*bet(0,i) = gam(0)**-1
C                            i=1
C
C      - alp(k-1)/gam(k) + (1-bet(k,k))*alp(k) - 
C
C                n
C             - Sum  alp(i)*bet(k,i) = 0, k = 1...n
C              i=k+1
C
C      |R| = |alp(n)/gam(n+1)| (provided that alp(i) satisfy
C                               equations above)
C
C      Iterations are continued either until |R| falls below
C      DCPHF threshold or storage for basis vectors is exhausted. If
C      more then one equations are solved simultaneosly, equations
C      with largest estimation of residual error might be discarded
C      to free space for remaining equations.
C
C
C      Set of DOTS matrices are used to store needed quantities.
C
C      DOTS( i,j) store corresponding bet(i,j)
C                     ~    ~
C      DOTS(-1,i) is <B(i)|B(i)>
C
C      DOTS(-2,i) is gam(i)
C
C      Note that vectors are actually computed for -A matrix, so
C      sign of correlation coefficients as well as resulting solution 
C      coefficients have to be adjusted to account for the fact.
C      Although this is not elegant, it introduces negligibly
C      small overhead.
C      
C
C   Speedups possible:
C
C      Using basis vectors originated from iterations on other
C      variables to build reduced system might be useful. According
C      to MJS Dewar, DA Liotard, J.Mol.Struct.(Theochem) 206(1990)
C      123-133 this might reduce number of basis vectors generated
C      during solution phase by a factor of 2.
C
C   Bugs:
C
C      Array of variables indices should be dynamically allocated.
C
C      Array for reduced linear system and solution vector should
C      be dynamically allocated.
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXVRS=3*LM1)
      PARAMETER (ZERO  = 0.0D0)
      PARAMETER (ONE   = 1.0D0)
      PARAMETER (SONE  =-1.0D0)
C
      LOGICAL SHUTUP, FULCND
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
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
     ./PSCPST/ NVECS,  NITERS, NWASTE, NMAXVE
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSDOPT/, /PSDGB2/, /PSPRTF/, /PSCPST/, /PSPRT /
C
      DIMENSION CALP(LDC,*), CBET(LDC,*), EALP(*), EBET(*), DUMP(*)
      DIMENSION ISCD(0:MAXVRS), DOTS(-2:MAXIT,0:MAXIT,NRHS)
      DIMENSION VECS(IVSZ,*)
C
      DIMENSION IVARS(MAXVRS), ERRS(MAXVRS), SHUTUP(MAXVRS)
C
      IF( IQSWAP.LT.2 ) IQ    = IPSMOF(LQ)
      IF( IQSWAP.EQ.3 ) IXTMP = IPSMOF(LXTMP)
      ICOND = IPSMOF(LCOND)
      ITEMP = IPSMOF(LCPIT1)
      ITA   = ITEMP
      ITX   = ITA + IMAXIT**2
      IPIV  = ITX + IMAXIT
      FULCND = .FALSE.
      IF( IPRECT.GT.1 ) FULCND = .TRUE.
C
      INEXTV = ICPV1
      IPASS  = 0
   10 CONTINUE
      IF(INEXTV.GT.ICPVL) THEN
          RETURN
      ENDIF
      IPASS  = IPASS + 1
C
C    First, pick up the variables we are going to solve for and
C    move them to the adjacent elements of temporary arrays.
C    Then right-hand sides are out-of-core, we lose nothing.
C    Then right-hand sides are in core, it's the price we have
C    to pay for been able to use Level 3 BLAS for doing matrix-
C    vector products
C
          IVARC = 0
          DO 100 I=INEXTV,ICPVL
              IF( ISCD(I).LE.0 ) THEN
                  ISCD(I)      = 0
                  IVARC        = IVARC + 1
                  IVARS(IVARC) = I
                  IF( IQSWAP.EQ.3 ) THEN
                      CALL PSDS1X(CALP,CBET,LDC,DUMP,I,VECS(1,IVARC))
                  ELSE IF( IQSWAP.EQ.2 ) THEN
                      CALL PSDRD(IURHS,I+(1-ICPV1),VECS(1,IVARC),IQSZ)
                  ELSE
                      CALL DCOPY(IQSZ,DUMP(IQ+(I-ICPV1)*IQSZ),1,
     .                           VECS(1,IVARC),1)
                  ENDIF
                  IF( IVARC.GE.NRHS ) GOTO 101
              ENDIF
  100     CONTINUE
  101     CONTINUE
          IF( IVARC.EQ.0 ) THEN
              RETURN
          ENDIF
C
C         Apply conditioner to all selected vectors at once
C
          CALL PSDS3C(DUMP(ICOND),VECS(1,1),IQSZ,IVARC)
          INEXTV = I+1
          IF( IPRINT.GE.1 ) THEN
              WRITE(NB6,10000) IPASS, INEXTV, INFO
          ENDIF
C
C         Zeroth iteration is special, do it separately
C
          DO 200 I=1,IVARC
              DOTS(-2,0,I) = DNRM2(IQSZ,VECS(1,I),1)
              IF( DOTS(-2,0,I).NE.ZERO ) THEN
                  CALL DSCAL(IQSZ,ONE/DOTS(-2,0,I),VECS(1,I),1)
              ENDIF
              DOTS(-1,0,I) = ONE
              ERRS(I)      = DOTS(-2,0,I)
              SHUTUP(I)    = .FALSE.
  200     CONTINUE
C
          DO 1000 ITER=1,MAXIT+1
C
C    Check for convergence
C
              NIVARC = 0
              DO 250 IX=1,IVARC
                  CALL PSDS3T(ITER-1,DOTS(-2,0,IX),MAXIT,DUMP(ITA),
     .                        DUMP(ITX),DUMP(IPIV),ERRS(IX),SHUTUP(IX))
                  IF(LPCPHF) THEN
                      WRITE(NB6,11200) ITER-1, IX, IVARS(IX), ERRS(IX)
                      CALL PSDPGM(1,ITER-1,DUMP(ITX),1)
                  ENDIF
                  IF( ERRS(IX).LE.DCPHF ) THEN
C
C                     Convergence achieved, compute solution vector
C                     and retire solved variable
C
                      INFO = INFO + 1
                      ISCD(IVARS(IX)) = ITER
                      IF( ITER-1.GT.NMAXVE ) NMAXVE = ITER-1
C
C                     Adjust sign of solution coefficients (see note above)
C
                      CALL DSCAL(ITER/2,SONE,DUMP(ITX),2)
                      CALL PSDS3X(ITER-1,DUMP(ITX),VECS(1,IX),
     .                            IQSZ*IVARC,IQSZ)
                      IF(FULCND) THEN
C
C                         Problem was conditioned, convert vector to the
C                         true solution. I can't condition all vectors at
C                         once, because this vector might be out-of-memory
C                         by the time next one is available.
C
                          CALL PSDS3C(DUMP(ICOND),VECS(1,IX),IQSZ,1)
                      ENDIF
*                     IF( IVARS(IX).GE.1 ) THEN
C
C                         All solutions which do not correspond to Z-vector(s)
C                         are to be casted back to the unsymmetrical set of
C                         variables.
C
*                         CALL PSDUNO(VECS(1,IX))
*                     ENDIF
                      IF(LPCPHF) THEN
                          WRITE(NB6,11300) IVARS(IX), ITER-1
                          CALL PSXPRT(VECS(1,IX))
                      ENDIF
                      IF( IQSWAP.EQ.3 ) THEN
                          CALL DCOPY(IQSZ,VECS(1,IX),1,DUMP(IXTMP),1)
                      ELSE IF( IQSWAP.EQ.2 ) THEN
                          CALL PSDWD(IURHS,IVARS(IX)+(1-ICPV1),
     .                               VECS(1,IX),IQSZ)
                      ELSE
                          CALL DCOPY(IQSZ,VECS(1,IX),1,
     .                               DUMP(IQ+(IVARS(IX)-ICPV1)*IQSZ),1)
                      ENDIF
                  ELSE
                      NIVARC = NIVARC + 1
                  ENDIF
  250         CONTINUE
              IF( NIVARC.EQ.0 ) GOTO 10
              IF( ITER.EQ.MAXIT+1) GOTO 1001
C
              IF( (ITER+1)*NIVARC.GE.IKRVEC ) THEN
C
C                 We have not enough temporary space to compute
C                 next generation of vectors for all source vectors.
C                 Mark few vectors for retirement and continue.
C
                  IDROP  = NIVARC - IKRVEC/(ITER+1)
                  NWASTE = NWASTE + IDROP*(ITER-1)
C
C                 Mark vectors with largest error estimation for 
C                 retirement. This code is extremly naive, but it's
C                 only N^2 - next step is at least N^4, so we probably
C                 shouldn't bother...
C
                  DO 280 I=1,IDROP
                      J = IDAMAX(IVARC,ERRS,1)
                      ISCD(IVARS(J)) = -ITER
                      IF(LPCPHF) THEN
                          WRITE(NB6,11600) J, IVARS(J), ERRS(J)
                      ENDIF
                      ERRS(J)        = ZERO
  280             CONTINUE
                  NIVARC = NIVARC - IDROP
              ENDIF
              IF( NIVARC.NE.IVARC ) THEN
C
C             Discard vectors retired on this step
C
                  IIN  = 1
                  IOUT = 1
                  DO 350 I=1,ITER
                      DO 348 IX=1,IVARC
                          IF( ABS(ISCD(IVARS(IX))).NE.ITER ) THEN
                              IF( IIN.NE.IOUT ) THEN
                                  CALL DCOPY(IQSZ,VECS(1,IIN),1,
     .                                            VECS(1,IOUT),1)
                              ENDIF
                              IOUT = IOUT + 1
                          ENDIF
                          IIN = IIN + 1
  348                 CONTINUE
  350             CONTINUE
                  IOUT = 1
                  DO 360 IX=1,IVARC
                      IF( ABS(ISCD(IVARS(IX))).NE.ITER ) THEN
                          IF( IOUT.NE.IX ) THEN
                              IVARS(IOUT)  = IVARS(IX)
                              ERRS (IOUT)  = ERRS (IX)
                              SHUTUP(IOUT) = SHUTUP(IX)
                              DO 356 K=0,ITER
                                  DO 354 J=-2,ITER
                                      DOTS(J,K,IOUT)=DOTS(J,K,IX)
  354                             CONTINUE
  356                         CONTINUE
                          ENDIF
                          IOUT = IOUT + 1
                      ENDIF
  360             CONTINUE
                  IVARC = NIVARC
              ENDIF
C
C             Compute next generation of Krylov vectors. This
C             step is _very_ expensive, the rest of the module
C             does next to nothing compared with it.
C
              IPREV = (ITER-1)*IVARC+1
              INEXT = ITER*IVARC+1
              CALL PSDS3L(CALP,CBET,LDC,EALP,EBET,DUMP,IVARC,
     .                    VECS(1,IPREV),VECS(1,INEXT),IQSZ)
              NITERS = NITERS + 1
              NVECS  = NVECS  + IVARC
C              
              DO 370 IX=1,IVARC
C
C                 Compute self-correlation coefficients and orthogonalise
C                 new vectors with respect to single previous generation.
C
                  INEW = INEXT+IX-1
C                Compute dot products needed to orthogonalize
C                new vector with respect to previous ones
                  CALL PSZRMB(ITER+0,1,DOTS(0,ITER-1,IX),
     .                                (MAXIT+1)*(MAXIT+3))
                  CALL DGEMM('T','N',ITER+0,1,IQSZ,ONE,VECS(1,IX),
     .                       IQSZ*IVARC,VECS(1,INEW),IQSZ,ZERO,
     .                       DOTS(0,ITER-1,IX),(MAXIT+1)*(MAXIT+3))
C
                  DO 364 I=0,ITER-1
                      CALL DAXPY(IQSZ,-DOTS(I,ITER-1,IX)/DOTS(-1,I,IX),
     .                           VECS(1,I*IVARC+IX),1,VECS(1,INEW),1)
  364             CONTINUE
C
C                 Adjust sign of dot products (see note above)
C
                  DO 366 I=ITER-1,0,-2
                      DOTS(I,ITER-1,IX) = -DOTS(I,ITER-1,IX)
  366             CONTINUE
C
C                 Renormalize new basis vector to one, save old norm
C
                  DOTS(-2,ITER,IX)   = DNRM2(IQSZ,VECS(1,INEW),1)
                  IF( DOTS(-2,ITER,IX).NE.ZERO ) THEN
                      CALL DSCAL(IQSZ,ONE/DOTS(-2,ITER,IX),
     .                                VECS(1,INEW),1)
                  ENDIF
C
C                 Norm of the new basis vector is always one.
C
                  DOTS(-1,ITER,IX) = ONE
  370         CONTINUE
              IF(LPCPHF) THEN
                  WRITE(NB6,11000) ITER, IVARC
                  DO 500 I=1,IVARC
                      WRITE(NB6,11010) I, IVARS(I)
                      CALL PSDPSU(ITER,DOTS(0,0,I),MAXIT+3)
                      WRITE(NB6,11011)
                      CALL PSDPGM(1,ITER,DOTS(-1,0,I),MAXIT+3)
                      WRITE(NB6,11012)
                      CALL PSDPGM(1,ITER+1,DOTS(-2,0,I),MAXIT+3)
  500             CONTINUE
              ENDIF
C
 1000     CONTINUE
 1001     CONTINUE
C
      IF( MAXIT.LT.IMAXIT ) GOTO 10
C
C     Any variables remaining had failed to converge with the maximum
C     iteration count allowed, mark 'em. Also there is no point in
C     continuing: we won't be able to converge these variables anyway
C
      DO 1100 IX=1,IVARC
          IF( ISCD(IVARS(IX)).EQ.0 ) THEN
              ISCD(IVARS(IX)) = -IMAXIT-1
          ENDIF
 1100 CONTINUE
      RETURN
C
10000 FORMAT(' CPHF PASS ', I5, ' NEXT VARIABLE ', I5, 
     .       ', VARIABLES SOLVED SO FAR ', I5 )
11000 FORMAT(' CPHF ITERATION ', I4, ' NUMBER OF VARIABLES ', I4/)
11010 FORMAT(' REDUCED BASIS CORRELATION VECTOR ', I4, 
     .       ' VARIABLE ', I4, ' IS:'/)
11011 FORMAT(' BASIS VECTORS NORMS ARE:'/)
11012 FORMAT(' SCALING COEFFICIENTS ARE:'/)
11100 FORMAT(' STATIC ARRAY OVERFLOW IN PSDS3S.'
     .      /' PROGRAM WILL STOP.' )
11200 FORMAT(' ON ITERATION ', I4, ' VECTOR ', I4, ' VARIABLE ',
     .       I4, ' (ERR = ', G16.7, ') IS:'/)
11300 FORMAT(' CPHF SOLUTION VECTOR ON VARIABLE ', I4, ' COMPUTED IN ',
     .       I4, ' ITERATIONS:'/)
11500 FORMAT(' PROGRAM BUG: ATTEMPT TO RETIRE VECTOR ',
     .       I4, ' VARIABLE ', I4, ' STATUS ', I4,
     .       ' WHICH WAS ALREADY RETIRED.'
     .      /' PROGRAM WILL STOP.' )
11600 FORMAT(' RETIRING VECTOR ', I4, ' VARIABLE ', I4, ' ERROR ',
     .       'ESTIMATION ', G14.7 )
      END
C
      SUBROUTINE PSDS3T(ITER,DOTS,MAXIT,A,X,APIV,ERR,SHUTUP)
C
C   Solve CPHF equations in subspace of Krylov vectors and
C   estimate residual error of solution.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ITER   - Number of iterations
C      DOTS   - Dot products of Krylov vectors over conditioned
C               matrix K
C      MAXIT  - Dimension of DOTS
C      A      - Scratch space for linear system
C      X      - Solution vector
C      APIV   - Scratch space for linear solver. Actually
C               it should be integer, and is used as such by
C               solver...
C      ERR    - Error estimation
C      SHUTUP - .TRUE. if no more warnings on partial loss of
C               precision should be printed.
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options
C      PSPRT  - Printing unit.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C      Passing DOUBLE PRECISION argument for IPIV parameter
C      of DGESV is a violation of standard.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (ONE =1.0D0)
C
      LOGICAL SHUTUP
C
      COMMON
     ./PSDOPT/ DSTORE, DSTEP,  DECONV, DPCONV, DCPHF, DPREC, DCDIFF,
     .         DSHIFT, DBASCR, DNCOFF,
     .         IUMIX,  IURHS,  IUK,    IURES,
     .         IPRINT, IENRG,  ICORE,  IDISK,  IMIX,   IDENS,
     .         INDSYM, IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC,
     .         IROWS,  IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP,
     .         IKMODE, ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT
     ./PSPRT / NB6
      SAVE /PSDOPT/, /PSPRT /
      DIMENSION DOTS(-2:MAXIT,0:MAXIT), A(0:ITER-1,0:ITER-1)
      DIMENSION X(0:*), APIV(*)
C
      IF( ITER.EQ.0 ) THEN
          ERR = DOTS(-2,0)
          RETURN
      ENDIF
C
      X(0) = DOTS(-2,0)
      DO 50 I=1,ITER-1
          X(I) = ZERO
   50 CONTINUE
C
      DO 80 I=0,ITER-1
          DO 60 J=0,I-2
              A(I,J) = ZERO
   60     CONTINUE
          IF( I.GT.0 ) THEN
              A(I,I-1) = -DOTS(-2,I)
          ENDIF
          A(I,I) = ONE - DOTS(I,I)
          DO 70 J=I+1,ITER-1
              A(I,J) = -DOTS(I,J)
   70     CONTINUE
   80 CONTINUE
C
      CALL DGESV(ITER,1,A,ITER,APIV,X,ITER,INFO)
      IF( INFO.LT.0 ) THEN
          WRITE(NB6,11000) INFO
          STOP 'PSDS3S'
      ENDIF
      IF( INFO.GE.1 ) THEN
          IF( IPRINT.GE.-1 ) THEN
              WRITE(NB6,11010) INFO
          ENDIF
          DO 150 I=0,ITER-1
              X(I) = ZERO
  150     CONTINUE
      ENDIF
C
C    Compute error vector in reduced basis
C
      XERR = ABS(X(ITER-1))*DOTS(-2,ITER)*SQRT(DOTS(-1,ITER))
C
      IF( XERR.GT.ERR .AND. IPRINT.GE.0 .AND. .NOT.SHUTUP ) THEN
          WRITE(NB6,11020) ITER, ERR, XERR
          IF( IPRINT.LE.1 ) THEN
              SHUTUP = .TRUE.
          ENDIF
      ENDIF
      ERR = XERR
C
      RETURN
11000 FORMAT(' DGESV FAILED WITH ERROR CODE ', I5, 'IN PSDS3S.',
     .      /' PROGRAM WILL STOP.')
11010 FORMAT(' WARNING: REDUCED SPACE MATRIX IS SINGULAR (INFO=',I4,
     .       ').'
     .      /' THIS MAY INDICATE TOO STRINGENT CPHF',
     .       ' CONVERGENCE CRITERIA.'
     .      /' CORRESPONDING SOLUTION VECTOR WILL BE ZEROED.')
11020 FORMAT('    WARNING: RESIDUAL ERROR INCREASED ON ITERATION ', I4, 
     .       ' FROM ', G14.7, ' TO ', G14.7 )
      END
C
C
      SUBROUTINE PSDS3X(ITER,X,VECS,LDV,IQSZ)
C
C   Compute CPHF solution vector from reduced space solution
C   and Krylov vectors
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ITER   - Number of iterations
C      X      - Solution vector
C      VECS   - Krylov vectors
C      LDV    - Leading dimension of VECS
C      IQSZ   - Length of single vector
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
C
      DIMENSION X(*), VECS(LDV,*)
C
      IF( ITER.EQ.0 ) THEN
          CALL PSZRV(IQSZ,VECS)
      ELSE
          CALL DSCAL(IQSZ,X(1),VECS(1,1),1)
          DO 100 I=2,ITER
              CALL DAXPY(IQSZ,X(I),VECS(1,I),1,VECS(1,1),1)
  100     CONTINUE
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSDS3L(CALP,CBET,LDC,EALP,EBET,DUMP,IVARC,AIN,OUT,LDA)
C
C   Do conditioning if necessary and call PSDS3M to compute new
C   generation of Krylov vectors. 
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients
C      CBET   - Beta orbital coefficients (not used if
C               UHF flag not set)
C      LDC    - Leading dimension of CALP and CBET
C      EALP   - Alpha orbital energies
C               Neither EALP not EBET are actually used,
C               but they are required parameters for
C               PSDS2P which is called for IDENS=6 (direct
C               MO-basis CPHF)
C      EBET   - Beta orbital energies
C      DUMP   - Base of dynamic memory items
C      IVARC  - Number of vectors to process simultaneously
C      AIN    - Input vectors from previous generation
C      OUT    - New vectors
C      LDA    - Leading dimension of AIN and OUT. Assumed to
C               be equal to IQSZ in a lot of places here!
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGB2 - "Response" computation parameters.
C      PSDYNM - Dynamic memory handles
C
C   Modified common blocks:
C
C      PSAXTM - Cumulative time taken by PSDS3L, updated if IPRINT.GE.2
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
      PARAMETER (SONE=-1.0D0)
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
     ./PSAXTM/ AXTIME
      SAVE /PSDOPT/, /PSDGB2/, /PSAXTM/
C
      DIMENSION CALP(LDC,*), CBET(LDC,*), EALP(*), EBET(*), DUMP(*)
      DIMENSION AIN(LDA,*), OUT(LDA,*)
C
      IF( IPRINT.GE.2 ) THEN
          CALL CPUSEC(AXSTRT)
      ENDIF
      ICOND  = IPSMOF(LCOND)
      IF( IPRECT.LT.2 ) THEN
          CALL PSDS3M(CALP,CBET,LDC,EALP,EBET,DUMP,IVARC,AIN,OUT,LDA)
      ELSE
          ICNTMP = IPSMOF(LCNTMP)
          IF( IDENS.GE.7 ) ISHIFT = IPSMOF(LSHIFT)
C
C         Convert vectors from conditioned space into normal MO-basis
C         needed to compute matrix-vector products.
C
          CALL DCOPY(IQSZ*IVARC,AIN,1,DUMP(ICNTMP),1)
          CALL PSDS3C(DUMP(ICOND),DUMP(ICNTMP),IQSZ,IVARC)
C
C         Compute matrix-vector products
C
          CALL PSDS3M(CALP,CBET,LDC,EALP,EBET,DUMP,IVARC,DUMP(ICNTMP),
     .            OUT,LDA)
C
C         For AO-basis methods, add diagonal shift correction
C
          IF( IDENS.GE.7 ) THEN
              CALL PSDS3C(DUMP(ISHIFT),DUMP(ICNTMP),IQSZ,IVARC)
              CALL DAXPY(IQSZ*IVARC,SONE,DUMP(ICNTMP),1,OUT,1)
          ENDIF
      ENDIF
C
C     Transform vectors back into conditioned basis
C
      CALL PSDS3C(DUMP(ICOND),OUT,IQSZ,IVARC)
      IF( IPRINT.GE.2 ) THEN
          CALL CPUSEC(AXEND)
          AXTIME = AXTIME + (AXEND-AXSTRT)
      ENDIF
      RETURN
      END
C
      SUBROUTINE PSDS3M(CALP,CBET,LDC,EALP,EBET,DUMP,IVARC,AIN,OUT,LDA)
C
C   Generate new set of Krylov vectors
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients
C      CBET   - Beta orbital coefficients (not used if
C               UHF flag not set)
C      LDC    - Leading dimension of CALP and CBET
C      EALP   - Alpha orbital energies
C               Neither EALP not EBET are actually used,
C               but they are required parameters for
C               PSDS2P which is called for IDENS=6 (direct
C               MO-basis CPHF)
C      EBET   - Beta orbital energies
C      DUMP   - Base of dynamic memory items
C      IVARC  - Number of vectors to process simultaneously
C      AIN    - Input vectors from previous generation
C      OUT    - New vectors
C      LDA    - Leading dimension of AIN and OUT. Assumed to
C               be equal to IQSZ in a lot of places here!
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGB2 - "Response" computation parameters.
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
C      This is dispatch routine for different methods of
C      computing next generation of Krylov vectors. Two
C      basic families are: MO-basis methods (IDENS=4,5,6)
C      and AO-basis methods (IDENS=7,8). Methods from the
C      same group should produce exactly the same vectors.
C
C      Actual computation is two-step. First, previous
C      generation vectors are multiplied from the left
C      by symmetric matrix of transformed ERIs. Then, 
C      vectors are multiplied by post-conditioned vector
C      (inverse differences in orbital energies) to 
C      simulate multiplication by unsymmetric scaled
C      matrix by Pople et al. As NVARS is usually much
C      less than IQSZ, this is both faster and requires
C      much less memory.
C
C   Bugs:
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ONE =1.0D0)
      PARAMETER (ZERO=0.0D0)
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
      SAVE /PSDOPT/, /PSDGB2/
C
      DIMENSION CALP(LDC,*), CBET(LDC,*), EALP(*), EBET(*), DUMP(*)
      DIMENSION AIN(LDA,*), OUT(LDA,*)
C
C    Compute products
C
      IF( IDENS.EQ.4 ) THEN
C
C    Straight version, just call BLAS function for each input vector
C    in turn, and that's all.
C
          IKS    = IPSMOF(LKS)
          DO 100 I=1,IVARC
              CALL PSZRVB(IQSZ,OUT(1,I))
              CALL DSPMV('U',IQSZ,ONE,DUMP(IKS),AIN(1,I),1,
     .                   ZERO,OUT(1,I),1)
  100     CONTINUE
      ELSE IF( IDENS.EQ.5 ) THEN
C
C    Out-of-core block packed matrix version. First, clear all output
C    vectors.
C
          CALL PSZRM(IQSZ,IVARC,OUT,LDA)
C
C    Now, read K matrix in block by block and update products
C
          IKATMP = IPSMOF(LKATMP)
          ISTART = 1
          IF( IROWS.LT.IQSZ ) THEN
              REWIND(IUK)
          ENDIF
  200     CONTINUE
              IEND   = IPSBPE(ISTART)
              ISIZE  = IEND*(IEND-ISTART+1)
              IF( IROWS.LT.IQSZ ) THEN
                  CALL PSDRS(IUK,DUMP(IKATMP),ISIZE)
              ENDIF
              CALL PSDS3N(ISTART,IEND,DUMP(IKATMP),IVARC,AIN,OUT,LDA)
              ISTART = IEND + 1
              IF( ISTART.LE.IQSZ ) GOTO 200
      ELSE IF( IDENS.EQ.6 ) THEN
C
C    Direct version. First, zero output vectors.
C
          CALL PSZRM(IQSZ,IVARC,OUT,LDA)
C
C    Now, build K matrix block by block and update products
C
          CALL PSDS2D(CALP,CBET,LDC,EALP,EBET,DUMP,IVARC,AIN,OUT,LDA)
      ELSE IF( IDENS.EQ.7 .OR. IDENS.EQ.8 ) THEN
C
C    AO basis version (this  should be a big winner for large systems -
C    it have smaller scaling power!)
C
          DO 800 I=1,IVARC
              CALL PSDS3A(CALP,CBET,LDC,DUMP,AIN(1,I),OUT(1,I))
  800     CONTINUE
      ELSE
          STOP 'PSDS3M'
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSDS3N(ISTART,IEND,AKATMP,IVARC,AIN,OUT,LDA)
C
C   Update new set of Krylov vectors using block of CPHF K
C   matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ISTART - Starting row of CPHF K matrix in block AKATMP
C      IEND   - Ending row of CPHF K matrix in block AKATMP
C      AKATMP - Block of CPHF K matrix
C      IVARC  - Number of vectors to process simultaneously
C      AIN    - Input vectors from previous generation
C      OUT    - New vectors
C      LDA    - Leading dimension of AIN and OUT. Assumed to
C               be equal to IQSZ in a lot of places here!
C
C   Accessed common blocks:
C
C      None.
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
C                     IROWS
C                    <----->
C
C            / - - - - - - - - \             v v v  ^
C            | .     y y y y   | ^           v v v  |
C            |   .   y y y y   | | ICOLS     v v v  |
C            |     . y y y y   | v           v v v  |
C         ^  | x x x t         |          X  v v v  | IQSZ
C   IROWS |  | x x x t t       |             v v v  |
C         |  | x x x t t t     |             v v v  |
C         v  | x x x t t t t   |             v v v  |
C            |               . |             v v v  v
C            \ - - - - - - - - /
C
C              <--->                         <---->
C              ICOLS                          IVARC
C
C   Bugs:
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ONE = 1.0D0)
C
      DIMENSION AKATMP(IEND,*), AIN(LDA,*), OUT(LDA,*)
C
      IROWS = IEND - ISTART + 1
      ICOLS = IEND - IROWS
      IF( ICOLS.GT.0 ) THEN
C
C         Our fragment of K matrix have rectangular part
C
          CALL DGEMM('N','N',ICOLS,IVARC,IROWS,ONE,AKATMP,IEND,
     .                       AIN(ISTART,1),LDA,ONE,OUT,LDA)
          CALL DGEMM('T','N',IROWS,IVARC,ICOLS,ONE,AKATMP,IEND,
     .                       AIN,LDA,ONE,OUT(ISTART,1),LDA)
      ENDIF
C
C     Process triangular part of K
C
      CALL DSYMM('L','U',IROWS,IVARC,ONE,AKATMP(ISTART,1),IEND,
     .                   AIN(ISTART,1),LDA,ONE,OUT(ISTART,1),LDA)
C
      RETURN
      END
C
      SUBROUTINE PSDS3A(CALP,CBET,LDC,DUMP,AIN,OUT)
C
C   Generate new Krylov vector by transformation to AO basis
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients
C      CBET   - Beta orbital coefficients (not used if
C               UHF flag not set)
C      LDC    - Leading dimension of CALP and CBET
C      DUMP   - Base of dynamic memory items
C      AIN    - Input Krylov vector
C      OUT    - Generated vector
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C      PSDYNM - Dynamic memory handles
C      PSPRTF - Debug output control flags
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
C      First, Krylov vectors are tarnsformed from AO basis into
C      MO basis. Then, "response" matrix (matrices in UHF case)
C      are constructed. Finally, response matrices are transformed
C      back into MO basis.
C
C      Construction of response matrix is formally equivalent to
C      formation of Fock matrices with AO form of Krylov vector
C      used instead of density matrix. We use separate version,
C      because a) it doesn't matter, this is N^3 step while 
C      transformations are N^4. b) I could use unpacked version
C      of "density" this way and scale it in any way I please.
C      c) my routine is as fast as standard one, at least on RS/6k
C
C   Bugs:
C
C      Calling this code with IQSZ.EQ.0 would be disastrous.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
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
C    ./PSPRT / NB6
C     SAVE /PSPRT /
      SAVE /PSDGBL/
C
      DIMENSION CALP(LDC,*), CBET(LDC,*), DUMP(*), AIN(*), OUT(*)
C
C    Dynamic memory items
C
      ICPAYA = IPSMOF(LCPAYA)
      ICPAY  = IPSMOF(LCPAY)
      ICPATA = IPSMOF(LCPATA)
      IF(UHF) THEN
          ICPAYB = IPSMOF(LCPAYB)
          ICPATB = IPSMOF(LCPATB)
      ELSE
          ICPAYB = ICPAYA
          ICPATB = ICPATA
      ENDIF
      IAI2T1 = IPSMOF(LAI2T1)
C
C    Convert Krylov vector to the AO basis. Upper triangles of
C    symmetrised version of matrices are computed. 
C
      CALL PSXM2A(CALP,CBET,LDC,AIN,DUMP(ICPAYA),DUMP(ICPAYB),
     .            NORBS,DUMP)
C
C    Complete construction of generalized density matrices and
C    compute response matrix in AO basis.
C
      CALL PSRESP(DUMP(ICPAYA),DUMP(ICPAYB),DUMP(ICPAY),DUMP(IAI2T1),
     .            DUMP(ICPATA),DUMP(ICPATB),NORBS,.TRUE.)
C
C    Transform "response" matrices back to the MO basis, forming
C    alpha and beta parts of new Krylov vector, respectively.
C
      CALL PSXU2M(CALP,CBET,LDC,DUMP(ICPATA),DUMP(ICPATB),NORBS,
     .            OUT,DUMP)
C     IF(LPSUMM) THEN
C         WRITE(NB6,11200)
C         CALL PSXPRT(OUT)
C     ENDIF
C
      RETURN
11200 FORMAT(' NEW KRYLOV VECTOR IS:'/)
      END
