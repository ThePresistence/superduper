C     ******************************************************************
C
C     Solve CPHF equations iteratively, sharing information on
C     solving for different right-hand sides.
C
C     ******************************************************************
      SUBROUTINE PSLINS(CALP,CBET,LDC,EALP,EBET,DUMP)
*     ENTRY PSLIN1 below
C
C   Prepare for iterative solution of CPHF equations with
C   information sharing between different RHS.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients
C      CBET   - Beta orbital coefficients (not used if UHF flag
C               not set)
C      LDC    - Leading dimensions of CALP and CBET
C      EALP   - Alpha orbital energies
C      EBET   - Beta orbital energies (not used if UHF flag
C               not set)
C      DUMP   - Dynamic memory arena
C      IEQ    - Variable to compute solution for (valid only if
C               called through the auxiliary entry point PSLIN1)
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options.
C      PSDGB2 - "Response" computation parameters.
C      PSDYNM - Dynamic memory handles
C      PSPRT  - Printing unit
C      PSPRTF - Debug output tuning flags
C
C   Modified common blocks:
C
C      PSDLNS - CPHF equations solution statistics. The following
C               information is collected:
C          DLSPAS - Calls to PSLISA
C          DLSRHS - Linear equations solved
C          DLSDOT - Dot products of the length IQSZ taken
C          DLSDAX - DAXPYs of the length IQSZ taken
C          DLSAUX - Auxiliary linear problems solved
C          DLSRES - Residuals computed
C          DLSORT - Orthogonalization problems solved
C          DLSSVD - SVD calls made
C          DLSSVW - Cumulative complexity of SVD calls
C          DLSCMP - Number of basis space compactions made
C          DLSITR - Iterations made
C          DLSVEC - Number of response vectors computed
C          DLSSVP - Number of basis vectors removed by SVD
C          DLSSZP - Number of basis vectors removed by total 
C                   basis size constraints
C          DLSSHP - Number of basis vectors removed by shared
C                   basis size constraints
C          DLSMXB - Maximum number of basis vectors used in
C                   solution
C          DLSMXI - Maximum number of iterations made in solving 
C                   a single batch of right-hand-sides.
C          DLSDFO - Number of defective (partial) orthogonalizations
C                   of basis vectors.
C
C          Counts are reset by the call to PSLINS (or the
C          first call to PSLIN1). Although these are essentially
C          integer counts, the range of the default INTEGER type
C          is insufficient to hold corresponding values on 32-bit
C          systems, so that DOUBLE PRECISION is used instead.
C
C   Local storage:
C
C   Module logic:
C
C      PSLINS groups right-hand sides of CPHF equations into
C      batches, does necessary preconditioning and submits
C      batches to PSLISA, which actually solves them.
C
C   Bugs:
C
C      Calls to PSLIN1 should only be made with IEQ changing
C      from ICPV1 to ICPVL with step 1, or we'll loose internal
C      state information. At the best, it will result in inefficiency...
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D0)
C
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      LOGICAL SINGLE
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
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
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
     ./PSDLNS/ DLSPAS, DLSRHS, DLSDOT, DLSDAX, DLSAUX, DLSRES,
     .         DLSORT, DLSSVD, DLSSVW, DLSCMP, DLSITR, DLSVEC,
     .         DLSSVP, DLSSZP, DLSSHP, DLSMXB, DLSMXI, DLSDFO
     ./PSAXTM/ AXTIME
      SAVE /PSDOPT/, /PSDGB2/, /PSPRTF/, /PSDLNS/, /PSAXTM/
      SAVE /PSPRT /
C
      DIMENSION CALP(LDC,*),CBET(LDC,*),EALP(*),EBET(*),DUMP(*)
C
      SAVE IBASSZ
C
C     Jump to main entry point - solve for all RHS
C
      IXV1   = ICPV1
      IXVL   = ICPVL
      SINGLE = .FALSE.
      GOTO 50
C
C     Special entry point - either fetch solution already computed,
C     or process a new batch of RHS's
C
      ENTRY PSLIN1(CALP,CBET,LDC,EALP,EBET,DUMP,IEQ)
          IF( IQSWAP.EQ.3 ) THEN
              IXTMP = IPSMOF(LXTMP)
          ENDIF
          IXV1   = ((IEQ-ICPV1)/INRHS)*INRHS + ICPV1
          IF( IEQ.GT.IXV1 ) THEN
C
C             Solution was already computed, copy it to the expected place.
C
              IF( IQSWAP.EQ.3 ) THEN
                  IRHS = IPSMOF(LC1RHS)
                  CALL DCOPY(IQSZ,DUMP(IRHS+(IEQ-IXV1)*IQSZ),1,
     .                            DUMP(IXTMP),1)
              ENDIF
              RETURN
          ENDIF
          IXVL   = MIN(ICPVL,IXV1+INRHS-1)
          SINGLE = .TRUE.
   50 CONTINUE
      IYY   = IPSMOF(LC1YY)
      IYR   = IPSMOF(LC1YR)
      IRR   = IPSMOF(LC1RR)
      IRB   = IPSMOF(LC1RB)
      IYB   = IPSMOF(LC1YB)
      IBB   = IPSMOF(LC1BB)
      IF( IQSWAP.GT.1 ) THEN
          IRHS = IPSMOF(LC1RHS)
      ELSE
          IQ   = IPSMOF(LQ)
      ENDIF
      ITMA  = IPSMOF(LC1TMA)
      ITMX  = IPSMOF(LC1TMX)
      IX    = IPSMOF(LC1X)
      IBAS  = IPSMOF(LC1BAS)
      IRSP  = IPSMOF(LC1RSP)
      ITMP  = IPSMOF(LC1TMP)
      ISNG  = IPSMOF(LC1SNG)
      ISVT  = IPSMOF(LC1SVT)
      ISCR  = IPSMOF(LC1SCR)
C
      ICOND = IPSMOF(LCOND)
C
      IF( IXV1.EQ.ICPV1 ) THEN
C
C         Clear static counts
C
          IBASSZ = 0
C
C         Initialize statistics
C
          DLSPAS = 0
          DLSRHS = 0
          DLSDOT = 0
          DLSDAX = 0
          DLSAUX = 0
          DLSRES = 0
          DLSORT = 0
          DLSSVD = 0
          DLSSVW = 0
          DLSCMP = 0
          DLSITR = 0
          DLSVEC = 0
          DLSSVP = 0
          DLSSZP = 0
          DLSSHP = 0
          DLSMXB = 0
          DLSMXI = 0
          DLSDFO = 0
          AXTIME = ZERO
      ENDIF
C
      DO 1000 IV1ST=IXV1,IXVL,INRHS
          IVLAST = MIN(IXVL,IV1ST+INRHS-1)
          IVCNT  = IVLAST - IV1ST + 1
          IF(LPCPHF) THEN
              WRITE(NB6,11000) IV1ST, IVLAST
          ENDIF
C
C         Get necessary RHS vectors and precondition 'em
C
          IF( IQSWAP.LE.1 ) THEN
              IRHS = IQ + (IV1ST-ICPV1)*IQSZ
          ELSE IF( IQSWAP.EQ.2 ) THEN
              DO 100 IV=IV1ST,IVLAST
                  CALL PSDRD(IURHS,IV+(1-ICPV1),
     .                             DUMP(IRHS+(IV-IV1ST)*IQSZ),IQSZ)
  100         CONTINUE
          ELSE IF( IQSWAP.EQ.3 ) THEN
              DO 110 IV=IV1ST,IVLAST
                  CALL PSDS1X(CALP,CBET,LDC,DUMP,IV,
     .                             DUMP(IRHS+(IV-IV1ST)*IQSZ))
  110         CONTINUE
          ELSE
              WRITE(NB6,10000) IQSWAP
              STOP 'PSLINS'
          ENDIF
          IF(LPCPHF) THEN
              WRITE(NB6,11100)
              CALL PSDPGM(IQSZ,IVCNT,DUMP(IRHS),IQSZ)
          ENDIF
          CALL PSDS3C(DUMP(ICOND),DUMP(IRHS),IQSZ,IVCNT)
          IF(LPCPHF) THEN
              WRITE(NB6,11200)
              CALL PSDPGM(IQSZ,IVCNT,DUMP(IRHS),IQSZ)
          ENDIF
C
C         Compute solutions (Golly! That's a *call*! :-)
C
          CALL PSLISA(CALP,CBET,LDC,EALP,EBET,DUMP,IV1ST,
     .                IVCNT,DUMP(IRHS),IBASSZ,DUMP(IBAS),DUMP(IRSP),
     .                IQSZ,DUMP(IYY),DUMP(IYR),DUMP(IYB),INRHS,
     .                DUMP(IRR),DUMP(IRB),DUMP(IBB),IKRVEC,
     .                DUMP(ITMA),DUMP(ITMX),DUMP(ITMP),
     .                DUMP(ISNG),DUMP(ISVT),DUMP(IX),DUMP(ISCR),INFO)
          IF( INFO.NE.0 ) THEN
              WRITE(NB6,10100) INFO
              STOP 'PSLINS'
          ENDIF
          IF(LPCPHF) THEN
              WRITE(NB6,11300)
              CALL PSDPGM(IQSZ,IVCNT,DUMP(IRHS),IQSZ)
          ENDIF
C
C         If necessary, postcondition RHS's, then store 'em
C
          IF( IPRECT.GT.1 ) THEN
              CALL PSDS3C(DUMP(ICOND),DUMP(IRHS),IQSZ,IVCNT)
              IF(LPCPHF) THEN
                  WRITE(NB6,11400)
                  CALL PSDPGM(IQSZ,IVCNT,DUMP(IRHS),IQSZ)
              ENDIF
          ENDIF
          IF( IQSWAP.EQ.2 ) THEN
              DO 900 IV=IV1ST,IVLAST
                  CALL PSDWD(IURHS,IV+(1-ICPV1),
     .                             DUMP(IRHS+(IV-IV1ST)*IQSZ),IQSZ)
  900         CONTINUE
          ENDIF
 1000 CONTINUE
C
      IF( SINGLE .AND. IQSWAP.EQ.3 ) THEN
          CALL DCOPY(IQSZ,DUMP(IRHS),1,DUMP(IXTMP),1)
      ENDIF
      IF( IVLAST.EQ.ICPVL ) THEN
          IF( IPRINT.GE.0 ) THEN
              WRITE(NB6,10200) 
     .             DLSPAS, DLSRHS, DLSITR, DLSMXI, DLSVEC,
     .             DLSVEC/NCPVRS,
     .             DLSSVP, DLSSZP, DLSMXB
          ENDIF
          IF( IPRINT.GE.1 ) THEN
              WRITE(NB6,10300) 
     .             DLSDOT, DLSDAX, DLSAUX, DLSRES,
     .             DLSORT, DLSDFO, DLSSVD, DLSSVW, DLSCMP, 
     .             DLSSHP
          ENDIF
          IF( IPRINT.GE.2 ) WRITE(NB6,10310) AXTIME
          IAVEIT = NINT(DBLE(DLSVEC)/NCPVRS)
          IMAXIT = (NINT(DLSMXB)+INRHS-1)/INRHS
      ENDIF
C
      RETURN
10000 FORMAT(' IQSWAP = ', I5, ' IS NOT EXPECTED IN PSLINS.')
10100 FORMAT(' ITERATIVE LINEAR SOLVER FAILED WITH ERROR CODE ', I5 )
10200 FORMAT(' ITERATIVE SOLVER STATISTICS:'/
     .'     CALLS TO SOLVER                                    ', F14.0/
     .'     LINEAR EQUATIONS SOLVED                            ', F14.0/
     .'     ITERATIONS MADE                                    ', F14.0/
     .'     MAXIMUM ITERATIONS MADE FOR A BATCH OF RHS         ', F14.0/
     .'     RESPONSE VECTORS COMPUTED                          ', F14.0/
     .'     AVERAGE RESPONSE VECTORS PER VARIABLE              ', F14.7/
     .'     BASIS VECTORS REJECTED BY SVD                      ', F14.0/
     .'     BASIS VECTORS REMOVED DUE TO MEMORY CONSTRAINTS    ', F14.0/
     .'     MAXIMUM NUMBER OF BASIS VECTORS USED               ', F14.0)
10300 FORMAT(/
     .'     DDOTs COMPUTED                                     ', F14.0/
     .'     DAXPYs COMPUTED                                    ', F14.0/
     .'     AUXILIARY LINEAR SYSTEMS SOLVED                    ', F14.0/
     .'     RESIDUALS CONPUTED                                 ', F14.0/
     .'     ORTHOGONALIZATION PROBLEMS SOLVED                  ', F14.0/
     .'     DEFECTIVE ORTHOGONALIZATIONS ENCOUNTERED           ', F14.0/
     .'     SVD CALLS MADE                                     ', F14.0/
     .'     SVD CUMULATIVE COMPLEXITY                          ', F14.0/
     .'     BASIS SPACE COMPACTIONS                            ', F14.0/
     .'     BASIS VECTORS REMOVED BY SHARED BASIS CONSTRAINTS  ', F14.0)
10310 FORMAT(/
     .'     TIME SPENT EVALUATING THE RESPONSE VECTORS (SEC.)  ', F14.3)
C
11000 FORMAT(' SOLVING CPHF EQUATIONS FOR VARIABLES ', I4, '-', I4 )
11100 FORMAT(' UNPRECONDITIONED RHS SET IS: ' )
11200 FORMAT(' PRECONDITIONED RHS SET IS: ' )
11300 FORMAT(' RAW SOLUTION IS: ' )
11400 FORMAT(' CONDITIONED SOLUTION IS: ' )
      END
C
      SUBROUTINE PSLISA(CALP,CBET,LDC,EALP,EBET,DUMP,IV1ST,
     .           IVCNT,RHS,IBASSZ,BAS,RSP,LDX,YY,YR,YB,LDY,
     .           RR,RB,BB,LDR,TMA,TMX,TMP,SNG,SVT,X,SCR,INFO)
C
C   Solve CPHF equations iteratively (hi-tech solver with information 
C   sharing)
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients
C      CBET   - Beta orbital coefficients (not used if UHF flag
C               not set)
C      LDC    - Leading dimensions of CALP and CBET
C      EALP   - Alpha orbital energies
C      EBET   - Beta orbital energies (not used if UHF flag
C               not set)
C      DUMP   - Dynamic memory arena
C      IV1ST  - Index of the first variable (used for reporting only)
C      IVCNT  - Number of right-hand sides in this batch
C      RHS    - Input: right-hand sides
C               Output: solution vectors
C      IBASSZ - Number of basis vectors in BAS/RSP (should be
C               preserved across calls to PSLISA)
C      BAS    - Basis vectors (IKRSAV first elements should be 
C               preserved between calls to PSLISA)
C      RSP    - Response vectors (IKRSAV first elements should be 
C               preserved between calls to PSLISA)
C      LDX    - Leading dimension of RHS, BAS and RSP
C      YY     - Dot products of the current residuals.
C      YR     - Dot products of RHS vectors with response vectors
C               First index - RHS, second index - response vector
C      YB     - Dot products of RHS vectors with basis vectors
C               First index - RHS, second index - basis vector
C      LDY    - Leading dimension of YY, YR and YB
C      RR     - Dot products of response vectors, should be preserved
C               across calls
C      RB     - Dot products of response vectors with basis vectors
C               First index - response vector, second index - basis 
C               vector. Should be preserved across calls.
C      BB     - Dot products of basis vectors. Should be preserved
C               across calls.
C      LDR    - Leading dimension of RR, RB and BB
C      TMA    - Scratch matrix for auxiliary linear problems,
C               LDRxLDR
C      TMX    - Scratch matrix for RHS of auxiliary linear problems,
C               LDRxLDY
C      TMP    - Scratch space for conventional linear solver
C      SNG    - Scratch for singular values of the basis set vectors
C      SVT    - Scratch for SVD
C      X      - Scratch for solution/residual vectors
C      SCR    - Cumulative scores of basis vector performance. Should
C               be preserved across calls.
C      INFO   - Status.
C               0 = solver had converged to the desired accuracy
C
C      SMALL (constant parameter defined below) determines the number of
C      rows retained in LSQ solution of quasi-CG problems.
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options.
C      PSDGB2 - "Response" computation parameters.
C      PSPRT  - Printing unit
C      PSPRTF - Debug output tuning flags
C
C   Modified common blocks:
C
C      PSDLNS - CPHF equations solution statistics
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C      Then exactly one basis vector is used in restart, inefficient
C      code will be used.
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXVRS=3*LM1+2)
      PARAMETER (MAXPRM=3*((MAXVRS+1)/2)+2)
C
      PARAMETER (ZERO= 0.0D0)
      PARAMETER (ONE = 1.0D0)
      PARAMETER (TWO = 2.0D0)
      PARAMETER (SONE=-1.0D0)
C
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      LOGICAL CGLIKE, MINRES, ORTRES, GENCG
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
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
     ./PSDLNS/ DLSPAS, DLSRHS, DLSDOT, DLSDAX, DLSAUX, DLSRES,
     .         DLSORT, DLSSVD, DLSSVW, DLSCMP, DLSITR, DLSVEC,
     .         DLSSVP, DLSSZP, DLSSHP, DLSMXB, DLSMXI, DLSDFO
      SAVE /PSDOPT/, /PSDGB2/, /PSPRTF/, /PSDLNS/, /PSPRT /
C
      DIMENSION CALP(LDC,*), CBET(LDC,*), EALP(*), EBET(*), DUMP(*)
      DIMENSION RHS(LDX,IVCNT), BAS(LDX,*), RSP(LDX,*), X(LDX,*)
      DIMENSION YY(LDY,*), YR(LDY,*), YB(LDY,*)
      DIMENSION RR(LDR,*), RB(LDR,*), BB(LDR,*)
      DIMENSION TMA(LDR,*), TMX(LDR,*), TMP(*)
      DIMENSION SNG(*), SVT(*), SCR(*)
C
      DIMENSION IVECID(MAXVRS), IPERMS(MAXPRM), DUMMY(1)
C
C    Initialize control variables
C
      CGLIKE = .FALSE.
      MINRES = .FALSE.
      ORTRES = .FALSE.
      GENCG  = .FALSE.
      IF( ISOLVE.GE.4 .AND. ISOLVE.LE.6 ) CGLIKE = .TRUE.
      IF( ISOLVE.EQ.2 .OR.  ISOLVE.EQ.4 ) MINRES = .TRUE.
      IF( ISOLVE.EQ.3 .OR.  ISOLVE.EQ.5 ) ORTRES = .TRUE.
      IF( ISOLVE.EQ.6 ) GENCG = .TRUE.
C
      IF( IPRINT.GE.1 .AND. IV1ST.EQ.ICPV1 ) THEN
          IF( ISOLVE.EQ.2 .AND. IKRSAV.EQ.0 ) 
     .                            WRITE(NB6,10700) IVCNT
          IF( ISOLVE.EQ.2 .AND. IKRSAV.GT.0 ) 
     .                            WRITE(NB6,10702) IVCNT, IKRSAV
          IF( ISOLVE.EQ.3 .AND. IKRSAV.EQ.0 )
     .                            WRITE(NB6,10704) IVCNT
          IF( ISOLVE.EQ.3 .AND. IKRSAV.GT.0 )
     .                            WRITE(NB6,10706) IVCNT, IKRSAV
          IF( ISOLVE.EQ.4 .AND. IKRSAV.EQ.0 ) 
     .                            WRITE(NB6,10708) IVCNT
          IF( ISOLVE.EQ.4 .AND. IKRSAV.GT.0 ) 
     .                            WRITE(NB6,10710) IVCNT, IKRSAV
          IF( ISOLVE.EQ.5 .AND. IKRSAV.EQ.0 ) 
     .                            WRITE(NB6,10712) IVCNT
          IF( ISOLVE.EQ.5 .AND. IKRSAV.GT.0 ) 
     .                            WRITE(NB6,10714) IVCNT, IKRSAV
          IF( ISOLVE.EQ.6 .AND. IKRSAV.EQ.0 ) 
     .                            WRITE(NB6,10716) IVCNT
          IF( ISOLVE.EQ.6 .AND. IKRSAV.GT.0 ) 
     .                            WRITE(NB6,10718) IVCNT, IKRSAV
      ENDIF
C
      ISVTSZ = 1024 + 3*INRHS**2 + MAX(3*INRHS+IQSZ,5*IQSZ-4,2*INRHS**2,
     .             INRHS**2+INRHS*IQSZ+INRHS)
      INFO   = 0 
      ILEFT  = IVCNT
      ITER   = 0
      DLSPAS = DLSPAS + 1
C
C    IVECID contains the original index assigned to the vector
C    which had occupied slot I in the RHS array then PSLISA was 
C    called.
C
      DO 50 I=1,IVCNT
          IVECID(I) = I
   50 CONTINUE
C
C    If CG-like methods are used, initialize current solution
C    vectors with zeros.
C
      IF(CGLIKE) THEN
          CALL PSZRM(IQSZ,IVCNT,X,LDX)
      ENDIF
C
C    Compute norms of right-hand side vectors, immediately retire
C    right-hand sides for which zero is an acceptable solution.
C
      IIN  = 1
C
      DO 100 I=1,IVCNT
          YY(IIN,IIN) = DDOT(IQSZ,RHS(1,IIN),1,RHS(1,IIN),1)
          ERR   = SQRT(YY(IIN,IIN))
          IF( ERR.LE.DCPHF ) THEN
C
C             Move vector from ITOP place to current one, and 
C             exchange indices, since zero is an acceptable 
C             solution for the current vector.
C
              IF(LPSOLV) THEN
                  WRITE(NB6,11000) IVECID(IIN), ERR
              ENDIF
              CALL DCOPY(IQSZ,RHS(1,ILEFT),1,RHS(1,IIN),1)
              IF( .NOT.CGLIKE ) THEN
C
C                For CG-like methods, current solution is already 
C                zero, so we don't have to bother.
C
                  CALL PSZRV(IQSZ,RHS(1,ILEFT))
              ENDIF
              ITMP          = IVECID(IIN)
              IVECID(IIN)   = IVECID(ILEFT)
              IVECID(ILEFT) = ITMP
              ILEFT         = ILEFT - 1
          ELSE
              IIN           = IIN + 1
          ENDIF
  100 CONTINUE
      IF(LPSOLV) THEN
          WRITE(NB6,11100) (YY(I,I), I=1,ILEFT)
      ENDIF
      IF( ILEFT.EQ.0 ) GOTO 9000
C
C    Restart: compute dot products with existing basis 
C    and response vectors. This step is not meaninful for
C    CG-like methods, which do not keep information between
C    batches.
C
      IF( IBASSZ.GT.0 ) THEN
          IF(CGLIKE) THEN
              WRITE(NB6,10050) IBASSZ
              STOP 'PSLISA'
          ENDIF
          DLSDOT = DLSDOT + IBASSZ*ILEFT
          IF( MINRES.OR.ORTRES ) THEN
              CALL PSZRMB(ILEFT,IBASSZ,YR,LDY)
          ENDIF
          IF(MINRES) THEN
              IF( ILEFT.EQ.1 ) THEN
                  CALL DGEMV('T',IQSZ,IBASSZ,ONE,RSP,LDX,
     .                                RHS,1,ZERO,YR,LDY)
              ELSE
                  CALL DGEMM('T','N',ILEFT,IBASSZ,IQSZ,ONE,RHS,LDX,
     .                               RSP,LDX,ZERO,YR,LDY)
              ENDIF
          ENDIF
          IF(ORTRES) THEN
              IF( ILEFT.EQ.1 ) THEN
                  CALL DGEMV('T',IQSZ,IBASSZ,ONE,BAS,LDX,
     .                                RHS,1,ZERO,YB,LDY)
              ELSE
                  CALL DGEMM('T','N',ILEFT,IBASSZ,IQSZ,ONE,RHS,LDX,
     .                               BAS,LDX,ZERO,YB,LDY)
              ENDIF
          ENDIF
      ENDIF
      IF(LPSOLV) THEN
          WRITE(NB6,11200) IBASSZ
          IF( IBASSZ.GT.0 ) THEN
              IF(ORTRES) THEN
                  WRITE(NB6,11300)
                  CALL PSDPGM(ILEFT,IBASSZ,YB,LDY)
              ENDIF
              IF(MINRES) THEN
                  WRITE(NB6,11400)
                  CALL PSDPGM(ILEFT,IBASSZ,YR,LDY)
              ENDIF
          ENDIF
      ENDIF
C
C    Prepare auxiliary linear system(s) on solution coefficients
C
  300 CONTINUE
          IF( IBASSZ.GT.0 ) THEN
              IF(MINRES) THEN
C
C                 Minimize norm of the residual (MINRES)
C
                  DO 320 I=1,IBASSZ
                      DO 318 J=1,IBASSZ
                          TMA(J,I) = RR(J,I)
  318                 CONTINUE
  320             CONTINUE
                  DO 330 I=1,ILEFT
                      DO 328 J=1,IBASSZ
                          TMX(J,I) = YR(I,J)
  328                 CONTINUE
  330             CONTINUE
                  IF(LPSOLV) THEN
                      WRITE(NB6,11500)
                  ENDIF
              ELSE IF(ORTRES) THEN
C
C                 Make residual orthogonal to the basis vectors (NLSSQ)
C
                  DO 350 I=1,IBASSZ
                      DO 348 J=1,IBASSZ
                          TMA(J,I) = RB(J,I)
  348                 CONTINUE
  350             CONTINUE
                  DO 360 I=1,ILEFT
                      DO 358 J=1,IBASSZ
                          TMX(J,I) = YB(I,J)
  358                 CONTINUE
  360             CONTINUE
                  IF(LPSOLV) THEN
                      WRITE(NB6,11600)
                  ENDIF
              ELSE IF(GENCG) THEN
C
C                 Make residual orthogonal to previous one (CG)
C
                  DO 380 I=1,IBASSZ
                      DO 378 J=1,IBASSZ
                          TMA(J,I) = RB(I,J) + RB(J,I)
  378                 CONTINUE
  380             CONTINUE
                  DO 390 I=1,ILEFT
                      DO 388 J=1,IBASSZ
                          TMX(J,I) = TWO*YB(I,J)
  388                 CONTINUE
  390             CONTINUE
              ELSE
                  WRITE(NB6,10000) ISOLVE
                  STOP 'PSLISA'
              ENDIF
              IF(LPSOLV) THEN
                  WRITE(NB6,11700)
                  CALL PSDPGM(IBASSZ,IBASSZ,TMA,LDR)
                  WRITE(NB6,11800)
                  CALL PSDPGM(IBASSZ,ILEFT,TMX,LDR)
              ENDIF
C
C             Solve auxiliary linear systems, compute residuals
C             and, if convergence was achieved, solution vectors.
C             Although norms of the residuals can in principle
C             be estimated without computing residuals, the
C             expression is numerically unstable then RHS and
C             basis vectors are unrelated. Additionally, in
C             most cases, we'll need residuals shortly afterwards,
C             so there is no big harm in computing them now.
C
              CALL DGESV(IBASSZ,ILEFT,TMA,LDR,TMP,TMX,LDR,INFO)
              DLSAUX = DLSAUX + 1
              IF( INFO.NE.0 ) THEN
                  WRITE(NB6,10100) INFO
                  STOP 'PSLISA'
              ENDIF
              IF(LPSOLV) THEN
                  WRITE(NB6,11900)
                  CALL PSDPGM(IBASSZ,ILEFT,TMX,LDR)
              ENDIF
C
C                Recompute basis vectors performance scores.
C
              DO 398 J=1,IBASSZ
                  SCR(J) = ZERO
                  DO 396 I=1,ILEFT
                      SCR(J) = SCR(J) + TMX(J,I)**2
  396             CONTINUE
  398         CONTINUE
C
C            Computing residues... 
C
              DLSDAX = DLSDAX + IBASSZ*ILEFT
              DLSRES = DLSRES + ILEFT
              IF(CGLIKE) THEN
C
C                If we are using CG-like method, residues can be updated
C                in-place. Additionally, we'll have to update current
C                solution vectors.
C
                  DLSDAX = DLSDAX + IBASSZ*ILEFT
                  DLSRES = DLSRES + ILEFT
                  IF( ILEFT.EQ.1 ) THEN
                      CALL DGEMV('N',IQSZ,IBASSZ,SONE,RSP,LDX,
     .                               TMX,1,ONE,RHS,1)
                      CALL DGEMV('N',IQSZ,IBASSZ,SONE,BAS,LDX,
     .                               TMX,1,ONE,X,1)
                  ELSE
                      CALL DGEMM('N','N',IQSZ,ILEFT,IBASSZ,SONE,RSP,
     .                           LDX,TMX,LDR,ONE,RHS,LDX)
                      CALL DGEMM('N','N',IQSZ,ILEFT,IBASSZ,SONE,BAS,
     .                           LDX,TMX,LDR,ONE,X,LDX)
                  ENDIF
              ELSE
                  DO 400 I=1,ILEFT
                      CALL DCOPY(IQSZ,RHS(1,I),1,BAS(1,IBASSZ+I),1)
  400             CONTINUE
                  IF( ILEFT.EQ.1 ) THEN
                      CALL DGEMV('N',IQSZ,IBASSZ,SONE,RSP,LDX,
     .                               TMX,1,ONE,BAS(1,IBASSZ+1),1)
                  ELSE
                      CALL DGEMM('N','N',IQSZ,ILEFT,IBASSZ,SONE,RSP,
     .                           LDX,TMX,LDR,ONE,BAS(1,IBASSZ+1),LDX)
                  ENDIF
              ENDIF
C
C            Testing for convergence...
C
              IIN = 1
              DO 420 I=1,ILEFT
                  INDO = IBASSZ + IIN
                  IF( LPSOLV.AND.LPCPHF ) THEN
                      WRITE(NB6,12100) IVECID(IIN)
                      IF(CGLIKE) THEN
                          CALL PSXPRT(RHS(1,IIN))
                      ELSE
                          CALL PSXPRT(BAS(1,INDO))
                      ENDIF
                  ENDIF
                  IF(CGLIKE) THEN
                      ERR = DNRM2(IQSZ,RHS(1,IIN),1)
                  ELSE
                      ERR = DNRM2(IQSZ,BAS(1,INDO),1)
                  ENDIF
                  DLSDOT = DLSDOT + 1
                  IF(LPSOLV) THEN
                      WRITE(NB6,12200) IVECID(IIN), ERR
                  ENDIF
                  IF( ERR.LE.DCPHF ) THEN
C
C                    If non-CG method is used, compute solution.
C                    Also exchange the current equation with the
C                    last one, so that yet unsolved equations are
C                    contigous.
C
                      DLSRHS = DLSRHS + 1
                      IF( IPRINT.GE.1 ) THEN
                         WRITE(NB6,10600) IV1ST+IVECID(IIN)-1, ITER, ERR
                      ENDIF
                      IF( ILEFT.NE.IIN ) THEN
                          CALL DCOPY(IQSZ,RHS(1,ILEFT),1,RHS(1,IIN),1)
                      ENDIF
                      IF( .NOT.CGLIKE ) THEN
                          DLSDAX = DLSDAX + IBASSZ
                          CALL PSZRVB(IQSZ,RHS(1,ILEFT))
                          CALL DGEMV('N',IQSZ,IBASSZ,SONE,BAS,LDX,
     .                               TMX(1,IIN),1,ZERO,RHS(1,ILEFT),1)
                      ENDIF
                      IF(CGLIKE) THEN
                          CALL DCOPY(IQSZ,X(1,IIN),1,RHS(1,ILEFT),1)
                          CALL DCOPY(IQSZ,X(1,ILEFT),1,X(1,IIN),1)
                          CALL DCOPY(IQSZ,RHS(1,ILEFT),1,X(1,ILEFT),1)
                      ENDIF
                      CALL DCOPY(IBASSZ,TMX(1,ILEFT),1,TMX(1,IIN),1)
                      CALL DCOPY(IBASSZ,YB(ILEFT,1),LDY,YB(IIN,1),LDY)
                      CALL DCOPY(IBASSZ,YR(ILEFT,1),LDY,YR(IIN,1),LDY)
                      IF( .NOT.CGLIKE ) THEN
                          CALL DCOPY(IQSZ,BAS(1,IBASSZ+ILEFT),1,
     .                                    BAS(1,INDO),1)
                      ENDIF
                      ITMP          = IVECID(ILEFT)
                      IVECID(ILEFT) = IVECID(IIN)
                      IVECID(IIN)   = ITMP
                      IF(LPSOLV) THEN
                          WRITE(NB6,12300) IVECID(ILEFT)
                          IF(LPCPHF) THEN
                              IF(CGLIKE) THEN
                                  CALL PSXPRT(X(1,ILEFT))
                              ELSE
                                  CALL PSXPRT(RHS(1,ILEFT))
                              ENDIF
                          ENDIF
                      ENDIF
                      ILEFT     = ILEFT - 1
                  ELSE
                      IIN = IIN + 1
                  ENDIF
  420         CONTINUE
              IF( ILEFT.EQ.0 ) GOTO 9000
          ELSE
C
C             Then there is no basis vectors yet, residuals are
C             equal to RHS vectors.
C
              DO 500 I=1,ILEFT
                  CALL DCOPY(IQSZ,RHS(1,I),1,BAS(1,I),1)
  500         CONTINUE
          ENDIF
C
C        If maximum allowed number of iterations was made, 
C        indicate unsuccessful termination.
C
          IF( ITER.GE.IMAXIT ) THEN
              IF( IPRINT.GE.0 ) THEN
                  WRITE(NB6,12400) ITER, ILEFT, IVCNT
              ENDIF
              INFO = ILEFT
              GOTO 9000
          ENDIF
C
C        If we are extending basis set, orthogonalize new vectors
C        with respect to existing ones (MINRES) or with respect to
C        response vectors (GENCG). ORTRES do not need orthogonalization,
C        since new residuals are already orthogonal to basis vectors.
C
          INEWBS = ILEFT
          IF( IBASSZ.GT.0 .AND. .NOT.ORTRES ) THEN
              IF(CGLIKE) THEN
                   DO 680 I=1,ILEFT
                       CALL DCOPY(IQSZ,RHS(1,I),1,BAS(1,IBASSZ+I),1)
  680              CONTINUE
              ENDIF
              IF( MINRES.OR.GENCG ) THEN
                  CALL PSZRMB(IBASSZ,INEWBS,TMX,LDR)
              ENDIF
              IF(MINRES) THEN
                  DO 700 I=1,IBASSZ
                      DO 690 J=1,IBASSZ
                          TMA(J,I) = BB(J,I)
  690                 CONTINUE
  700             CONTINUE
                  DLSDOT = DLSDOT + IBASSZ*INEWBS
                  IF( INEWBS.EQ.1 ) THEN
                      CALL DGEMV('T',IQSZ,IBASSZ,SONE,BAS,LDX,
     .                               BAS(1,IBASSZ+1),1,ZERO,TMX,1)
                  ELSE
                      CALL DGEMM('T','N',IBASSZ,INEWBS,IQSZ,SONE,BAS,
     .                           LDX,BAS(1,IBASSZ+1),LDX,ZERO,TMX,LDR)
                  ENDIF
              ENDIF
C
              IF(GENCG) THEN
                  DO 750 I=1,IBASSZ
                      DO 740 J=1,IBASSZ
                          TMA(J,I) = RB(J,I)
  740                 CONTINUE
  750             CONTINUE
                  DLSDOT = DLSDOT + IBASSZ*INEWBS
                  IF( INEWBS.EQ.1 ) THEN
                      CALL DGEMV('T',IQSZ,IBASSZ,SONE,RSP,LDX,
     .                               BAS(1,IBASSZ+1),1,ZERO,TMX,1)
                  ELSE
                      CALL DGEMM('T','N',IBASSZ,INEWBS,IQSZ,SONE,RSP,
     .                           LDX,BAS(1,IBASSZ+1),LDX,ZERO,TMX,LDR)
                  ENDIF
              ENDIF
C
              IF(LPSOLV) THEN
                  WRITE(NB6,12900)
                  CALL PSDPGM(IBASSZ,IBASSZ,TMA,LDR)
                  WRITE(NB6,13000)
                  CALL PSDPGM(IBASSZ,INEWBS,TMX,LDR)
              ENDIF
              CALL DGESV(IBASSZ,INEWBS,TMA,LDR,TMP,TMX,LDR,INFO)
              DLSORT = DLSORT + 1
              IF( INFO.NE.0 ) THEN
                  WRITE(NB6,10400) INFO
                  STOP 'PSLISA'
              ENDIF
C
C                Compute composite basis vectors performance scores
C
              DO 770 J=1,IBASSZ
                  DO 768 I=1,INEWBS
                      SCR(J) = SCR(J) + TMX(J,I)**2
  768             CONTINUE
  770         CONTINUE
C
              DLSDAX = DLSDAX + IBASSZ*INEWBS
              IF( INEWBS.EQ.1 ) THEN
                  CALL DGEMV('N',IQSZ,IBASSZ,ONE,BAS,LDX,TMX,1,
     .                           ONE,BAS(1,IBASSZ+1),1)
              ELSE
                  CALL DGEMM('N','N',IQSZ,INEWBS,IBASSZ,ONE,BAS,
     .                       LDX,TMX,LDR,ONE,BAS(1,IBASSZ+1),LDX)
              ENDIF
              IF( LPSOLV.AND.LPCPHF ) THEN
                  WRITE(NB6,13200)
                  CALL PSDPGM(IQSZ,INEWBS,BAS(1,IBASSZ+1),LDX)
              ENDIF
          ENDIF
C
C        For CG-like methods, (part of) old basis vectors can be 
C        discarded now.
C
          IF( CGLIKE .AND. IBASSZ.GT.0 ) THEN
C
C            If there is an orthogonalization trail, preserve the
C            basis vectors with best scores.
C
              IDEST = 1
              ISRC  = IBASSZ + 1
              IF( IKRSAV.GT.0 ) THEN
                  IF( IBASSZ.GT.IKRSAV ) THEN
                      IDELV  = IBASSZ - IKRSAV
                      DLSCMP = DLSCMP + 1
                      DLSSZP = DLSSZP + IDELV
                      CALL PSLTRM(IVCNT,IBASSZ,IDELV,SCR,BAS,
     .                            RSP,LDX,YR,YB,LDY,BB,RB,RR,LDR)
                      IBASSZ = IKRSAV
                  ENDIF
                  IDEST  = IBASSZ + 1
              ELSE
                  IBASSZ = 0
              ENDIF
C
              IF( .NOT.ORTRES ) THEN
                  DO 800 I=1,ILEFT
                      CALL DCOPY(IQSZ,BAS(1,ISRC+I-1),1,
     .                                BAS(1,IDEST+I-1),1)
  800             CONTINUE
              ENDIF
              IF(ORTRES) THEN
                  DO 810 I=1,ILEFT
                      CALL DCOPY(IQSZ,RHS(1,I),1,BAS(1,IDEST+I-1),1)
  810             CONTINUE
              ENDIF
          ENDIF
C
C        Renormalize basis vectors before extracting optimal basis.
C        This is absolute necessity for CG, where computing conjugate
C        direction might result in wild variations in the length
C        of direction vector.
C
          IF( INEWBS.GT.0 ) THEN
              DLSDOT = DLSDOT + INEWBS
              IOUT   = 1
              DO 850 I=1,INEWBS
                  DNRM = DNRM2(IQSZ,BAS(1,IBASSZ+I),1)
                  IF(LPSOLV) THEN
                      WRITE(NB6,13300) I, DNRM
                  ENDIF
                  IF( DNRM.GT.ZERO ) THEN
                      DLSDAX = DLSDAX + 1
                      CALL DSCAL(IQSZ,ONE/DNRM,BAS(1,IBASSZ+I),1)
                      IF( I.NE.IOUT ) THEN
                          CALL DCOPY(IQSZ,BAS(1,IBASSZ+I),1,
     .                                    BAS(1,IBASSZ+IOUT),1)
                      ENDIF
                      IOUT = IOUT + 1
                  ELSE
                      INEWBS = INEWBS - 1
                  ENDIF
  850         CONTINUE
          ENDIF
C
C        If number of vectors not yet converged is more than 1,
C        do SVD on new residual vectors to extract optimal basis.
C        When all vectors are already orthogonal with respect
C        to existing basis, optimal basis should remain so.
C
          IF( INEWBS.GT.1 ) THEN
              DLSSVD = DLSSVD + 1
              DLSSVW = DLSSVW + INEWBS**2
              CALL DGESVD('O','N',IQSZ,INEWBS,BAS(1,IBASSZ+1),LDX,SNG,
     .                            DUMMY,1,DUMMY,1,SVT,ISVTSZ,INFO)
              IF( INFO.NE.0 ) THEN
                  WRITE(NB6,10200) INFO
                  STOP 'PSLISA'
              ENDIF
              IF(LPSOLV) THEN
                  WRITE(NB6,12500) (SNG(I),I=1,INEWBS)
                  WRITE(NB6,12600) NINT(SVT(1)), ISVTSZ
                  IF(LPCPHF) THEN
                      WRITE(NB6,12700)
                      CALL PSDPGM(IQSZ,INEWBS,BAS(1,IBASSZ+1),LDX)
                  ENDIF
              ENDIF
              IF( NINT(SVT(1)).GT.ISVTSZ .AND. IPRINT.GE.-1 ) THEN
                  WRITE(NB6,10300) NINT(SVT(1)), ISVTSZ
              ENDIF
C
C             Delete vectors with singular values too small
C             This loop should start from 1, so that we'll
C             catch the case of basis vectors vanishing completely.
C
              DO 900 I=1,INEWBS
                  IF( SNG(I).LE.DBASCR*SNG(1) ) THEN
                      INEWBS = I - 1
                      GOTO 910
                  ENDIF
  900         CONTINUE
  910         CONTINUE
              IF(LPSOLV) THEN
                  WRITE(NB6,12800) ILEFT - INEWBS
              ENDIF
              DLSSVP = DLSSVP + ILEFT - INEWBS
          ENDIF
C
          IF( INEWBS.LE.0 ) THEN
              WRITE(NB6,10500)
              STOP 'PSLISA'
          ENDIF
C
C        If the total number of basis vectors is too large, kill ones
C        with the worst performance scores. 
C
          IF( IBASSZ+INEWBS.GT.IKRVEC ) THEN
              IF(CGLIKE) THEN
                  WRITE(NB6,10550)
                  STOP 'PSLISA'
              ENDIF
              IDELV  = IBASSZ+INEWBS-IKRVEC
              DLSCMP = DLSCMP + 1
              DLSSZP = DLSSZP + IDELV
              CALL PSLTRM(ILEFT,IBASSZ,IDELV,SCR,BAS,RSP,LDX,
     .                                 YR,YB,LDY,BB,RB,RR,LDR)
C
C            Push down new basis vectors
C
              IOUT = IBASSZ - IDELV + 1
              DO 1200 IIN=IBASSZ+1,IBASSZ+INEWBS
                  IF( IIN.NE.IOUT ) THEN
                      CALL DCOPY(IQSZ,BAS(1,IIN),1,BAS(1,IOUT),1)
                  ENDIF
                  IOUT = IOUT + 1
 1200         CONTINUE
              IBASSZ = IBASSZ - IDELV
          ENDIF
C
C        Update dot products of current residual vectors (YY)
C
*         IF(.FALSE.) THEN
*             DLSDOT = DLSDOT + (ILEFT*(ILEFT+1))/2
*             INEWV  = IBASSZ+1
*             IF( ILEFT.EQ.1 ) THEN
*                 YY(1,1) = DDOT(IQSZ,RHS,1,RHS,1)
*             ELSE
*                 CALL DSYRK('U','T',ILEFT,IQSZ,ONE,RHS,LDX,ZERO,YY,LDY)
*             ENDIF
C
C            Clone symmetric parts
C
*             DO 1300 I=1,ILEFT
*                 DO 1290 J=I+1,ILEFT
*                     YY(J,I) = YY(I,J)
*1290             CONTINUE
*1300         CONTINUE
*         ENDIF
C
C        Update dot products of basis vectors. (BB and YB goes here)
C
          IF(MINRES) THEN
              DLSDOT = DLSDOT + INEWBS*IBASSZ 
     .                        + DBLE((INEWBS*(INEWBS+1))/2)
              INEWV  = IBASSZ+1
              CALL PSZRMB(IBASSZ,INEWBS,BB(1,INEWV),LDR)
              CALL PSZRMB(INEWBS,INEWBS,BB(INEWV,INEWV),LDR)
              IF( INEWBS.EQ.1 ) THEN
                  CALL DGEMV('T',IQSZ,INEWV,ONE,BAS,LDX,BAS(1,INEWV),1,
     .                                ZERO,BB(1,INEWV),1)
              ELSE
                  CALL DGEMM('T','N',IBASSZ,INEWBS,IQSZ,ONE,BAS,LDX,
     .                           BAS(1,INEWV),LDX,ZERO,BB(1,INEWV),LDR)
                  CALL DSYRK('U','T',INEWBS,IQSZ,ONE,BAS(1,INEWV),LDX,
     .                           ZERO,BB(INEWV,INEWV),LDR)
              ENDIF
C
C            Clone symmetric parts
C
              DO 1500 I=1,INEWBS
                  INEWV = IBASSZ + I
                  DO 1490 J=1,INEWV-1
                      BB(INEWV,J) = BB(J,INEWV)
 1490             CONTINUE
 1500         CONTINUE
          ENDIF
          IF( ORTRES.OR.GENCG ) THEN
C
C            For CG-like methods, all (Y,B) products are to
C            be recomputed, since Y were updated.
C
              IF(CGLIKE) THEN
                  INEWST = 1
                  INEWL  = IBASSZ+INEWBS
              ELSE
                  INEWST = IBASSZ+1
                  INEWL  = INEWBS
              ENDIF
C
              DLSDOT = DLSDOT + INEWL*ILEFT
              CALL PSZRMB(ILEFT,INEWL,YB(1,INEWST),LDY)
              IF( ILEFT.EQ.1 ) THEN
                  CALL DGEMV('T',IQSZ,INEWL,ONE,BAS(1,INEWST),
     .                       LDX,RHS,1,ZERO,YB(1,INEWST),LDY)
              ELSE
                  CALL DGEMM('T','N',ILEFT,INEWL,IQSZ,ONE,RHS,LDX,
     .                      BAS(1,INEWST),LDX,ZERO,YB(1,INEWST),LDY)
              ENDIF
          ENDIF
          IF(LPSOLV) THEN
              IF(MINRES) THEN
                  WRITE(NB6,13500)
                  CALL PSDPGM(IBASSZ+INEWBS,IBASSZ+INEWBS,BB,LDR)
              ENDIF
              IF( ORTRES.OR.GENCG ) THEN
                  WRITE(NB6,13600)
                  CALL PSDPGM(IVCNT,IBASSZ+INEWBS,YB,LDY)
              ENDIF
          ENDIF
C
C        For all new basis vectors, compute corresponding response
C        vectors.
C
          CALL PSDS3L(CALP,CBET,LDC,EALP,EBET,DUMP,INEWBS,
     .                    BAS(1,IBASSZ+1),RSP(1,IBASSZ+1),LDX)
          DLSITR = DLSITR + 1
          DLSVEC = DLSVEC + INEWBS
C
C        PSDS3L computes R'=(A-1)B, convert 'em to R by adding basis
C        vectors
C
          DLSDAX = DLSDAX + INEWBS
          DO 1700 I=1,INEWBS
              CALL DAXPY(IQSZ,ONE,BAS(1,IBASSZ+I),1,RSP(1,IBASSZ+I),1)
 1700     CONTINUE
C
          IF( LPSOLV.AND.LPCPHF ) THEN
              WRITE(NB6,13700)
              CALL PSDPGM(IQSZ,INEWBS,RSP(1,IBASSZ+1),LDX)
          ENDIF
C
C        Update dot products of response vectors. (RB, RR and YR goes here)
C
          IF(MINRES) THEN
              DLSDOT = DLSDOT + INEWBS*IBASSZ 
     .                        + DBLE((INEWBS*(INEWBS+1))/2)
              INEWV  = IBASSZ+1
              CALL PSZRMB(IBASSZ,INEWBS,RR(1,INEWV),LDR)
              CALL PSZRMB(INEWBS,INEWBS,RR(INEWV,INEWV),LDR)
              IF( INEWBS.EQ.1 ) THEN
                  CALL DGEMV('T',IQSZ,INEWV,ONE,RSP,LDX,RSP(1,INEWV),1,
     .                                ZERO,RR(1,INEWV),1)
              ELSE
                  CALL DGEMM('T','N',IBASSZ,INEWBS,IQSZ,ONE,RSP,LDX,
     .                           RSP(1,INEWV),LDX,ZERO,RR(1,INEWV),LDR)
                  CALL DSYRK('U','T',INEWBS,IQSZ,ONE,RSP(1,INEWV),LDX,
     .                           ZERO,RR(INEWV,INEWV),LDR)
              ENDIF
C
C            Clone symmetric parts
C
              DO 1800 I=1,INEWBS
                  INEWV = IBASSZ + I
                  DO 1790 J=1,INEWV-1
                      RR(INEWV,J) = RR(J,INEWV)
 1790             CONTINUE
 1800         CONTINUE
          ENDIF
C
          IF( ORTRES.OR.GENCG ) THEN
              DLSDOT = DLSDOT + (2*IBASSZ+INEWBS)*INEWBS
              INEWV  = IBASSZ + 1
              CALL PSZRMB(INEWBS,IBASSZ+INEWBS,RB(INEWV,1),LDR)
              CALL PSZRMB(IBASSZ,INEWBS,RB(1,INEWV),LDR)
              IF( INEWBS.EQ.1 ) THEN
                  CALL DGEMV('T',IQSZ,IBASSZ+1,ONE,BAS,LDX,
     .                            RSP(1,INEWV),1,ZERO,RB(INEWV,1),LDR)
                  CALL DGEMV('T',IQSZ,IBASSZ,ONE,RSP,LDX,
     .                            BAS(1,INEWV),1,ZERO,RB(1,INEWV),1)
              ELSE
                  CALL DGEMM('T','N',INEWBS,IBASSZ+INEWBS,IQSZ,ONE,
     .                               RSP(1,INEWV),LDX,BAS,LDX,ZERO,
     .                               RB(INEWV,1),LDR)
                  CALL DGEMM('T','N',IBASSZ,INEWBS,IQSZ,ONE,RSP,LDX,
     .                               BAS(1,INEWV),LDX,ZERO,
     .                               RB(1,INEWV),LDR)
              ENDIF
          ENDIF
C    
          IF( MINRES.OR.GENCG ) THEN
C
C            For CG-like methods, all (Y,R) products are to
C            be recomputed, since Y were updated.
C
              IF(CGLIKE) THEN
                  INEWST = 1
                  INEWL  = IBASSZ+INEWBS
              ELSE
                  INEWST = IBASSZ+1
                  INEWL  = INEWBS
              ENDIF
C
              DLSDOT = DLSDOT + INEWL*ILEFT
              CALL PSZRMB(ILEFT,INEWL,YR(1,INEWST),LDY)
              IF( ILEFT.EQ.1 ) THEN
                  CALL DGEMV('T',IQSZ,INEWL,ONE,RSP(1,INEWST),
     .                       LDX,RHS,1,ZERO,YR(1,INEWST),LDY)
              ELSE
                  CALL DGEMM('T','N',ILEFT,INEWL,IQSZ,ONE,RHS,LDX,
     .                      RSP(1,INEWST),LDX,ZERO,YR(1,INEWST),LDY)
              ENDIF
          ENDIF
          IBASSZ = IBASSZ + INEWBS
          IF( IBASSZ.GT.DLSMXB ) DLSMXB = IBASSZ
          IF(LPSOLV) THEN
              IF( ORTRES.OR.GENCG ) THEN
                  WRITE(NB6,13800)
                  CALL PSDPGM(IBASSZ,IBASSZ,RB,LDR)
              ENDIF
              IF(MINRES) THEN
                  WRITE(NB6,13900)
                  CALL PSDPGM(IBASSZ,IBASSZ,RR,LDR)
              ENDIF
              IF( MINRES.OR.GENCG ) THEN
                  WRITE(NB6,14000)
                  CALL PSDPGM(ILEFT,IBASSZ,YR,LDY)
              ENDIF
*             IF(.FALSE.) THEN
*                 WRITE(NB6,14050)
*                 CALL PSDPGM(ILEFT,ILEFT,YY,LDY)
*             ENDIF
          ENDIF
C
C        Repeat.
C
          ITER = ITER + 1
          IF( ITER.GT.DLSMXI ) DLSMXI = ITER
          GOTO 300
C
C    Trim number of basis vectors to IKRSAV using precomputed 
C    performance index.
C
 9000 CONTINUE
      IF( .NOT.CGLIKE .AND. IBASSZ.GT.IKRSAV ) THEN
          DLSCMP = DLSCMP + 1
          DLSSHP = DLSSHP + IBASSZ - IKRSAV
          CALL PSLTRM(IVCNT,IBASSZ,IBASSZ-IKRSAV,SCR,BAS,RSP,LDX,
     .                       YR,YB,LDY,BB,RB,RR,LDR)
          IBASSZ = IKRSAV
      ENDIF
      IF(CGLIKE) THEN
          IBASSZ = 0
      ENDIF
C
C    Put vectors in the places they actually belong.
C    This is done in two steps: First, build minimal set of
C    loops, then actually rearrange the elements. Since 
C    PSMKPT builds the list for the "come-from" problem, and
C    our ordiring table is "go-to", reverce ordering call
C    will be used.
C
      IF(LPSOLV) THEN
          WRITE(NB6,14100) (IVECID(I),I=1,IVCNT)
      ENDIF
      IF(CGLIKE) THEN
          DO 9100 I=1,IVCNT
              CALL DCOPY(IQSZ,X(1,I),1,RHS(1,IVECID(I)),1)
 9100     CONTINUE
      ELSE
          CALL PSMKPT(IVCNT,IVECID,IPERMS)
          IF(LPSOLV) THEN
              WRITE(NB6,14200) (IPERMS(I),I=2,IPERMS(1)+1)
          ENDIF
          CALL PSUSRV(IPERMS,RHS,LDX,X(1,1),IQSZ)
      ENDIF
      RETURN
C
10000 FORMAT(' UNSUPORTED SOLVER TYPE ', I5, ' IN PSLISA')
10050 FORMAT(' INITIAL BASIS IS NON-EMPTY FOR SOLVER ', I5 )
10100 FORMAT(' DGESV/DGELSS FAILED FOR AUXILIARY SYSTEM WITH ERROR',
     .       ' CODE ', I5, ' IN PSLISA' )
10200 FORMAT(' DGESVD FAILED WITH ERROR CODE ', I5, ' IN PSLISA' )
10300 FORMAT(' WARNING: SUBOPTIMAL ROUTE USED BY DGESVD. BUFFER ',
     .       'REQUIRED TO USE THE OPTIMAL '/' ROUTE IS ', I8, 
     .       ' WORDS, WHILE ONLY ', I8, ' WORDS WERE ALLOCATED.' )
10400 FORMAT(' DGESV FAILED FOR ORTHOGONALIZATION PROBLEM WITH ',
     .       'ERROR CODE ', I5, ' IN PSLISA' )
10500 FORMAT(' BASIS SPACE COLLAPSED IN PSLISA.' )
10550 FORMAT(' VECTOR SPACE LIMIT BUMPED FOR CG-LIKE METHOD IN PSLISA')
10600 FORMAT(' CPHF SOLUTION FOR VARIABLE ', I5, ' CONVERGED AFTER ', 
     .       I5, ' ITERATIONS (ERR = ', G12.6, ')' )
10700 FORMAT(' USING GLOBALLY ORTHOGONAL MINRES SOLVER. ',
     .       I4, ' EQUATIONS WILL BE SOLVED AT ONCE ' )
10702 FORMAT(' USING GLOBALLY ORTHOGONAL MINRES SOLVER. ',
     .            I4, ' EQUATIONS WILL BE SOLVED AT ONCE. '/
     .       ' ', I4, ' BASIS VECTORS WILL BE SHARED BETWEEN PASSES.')
10704 FORMAT(' USING GLOBALLY ORTHOGONAL ORTRES SOLVER. ',
     .       I4, ' EQUATIONS WILL BE SOLVED AT ONCE ' )
10706 FORMAT(' USING GLOBALLY ORTHOGONAL ORTRES SOLVER. ',
     .            I4, ' EQUATIONS WILL BE SOLVED AT ONCE. '/
     .       ' ', I4, ' BASIS VECTORS WILL BE SHARED BETWEEN PASSES.')
10708 FORMAT(' USING LOCALLY ORTHOGONAL MINRES SOLVER. ',
     .       I4, ' EQUATIONS WILL BE SOLVED AT ONCE ' )
10710 FORMAT(' USING LOCALLY ORTHOGONAL MINRES SOLVER. ',
     .            I4, ' EQUATIONS WILL BE SOLVED AT ONCE. '/
     .       ' ', I4, ' BASIS VECTORS WILL BE USED AS A ',
     .            'STABILIZATION TRAIL.')
10712 FORMAT(' USING LOCALLY ORTHOGONAL ORTRES SOLVER. ',
     .       I4, ' EQUATIONS WILL BE SOLVED AT ONCE ' )
10714 FORMAT(' USING LOCALLY ORTHOGONAL ORTRES SOLVER. ',
     .            I4, ' EQUATIONS WILL BE SOLVED AT ONCE. '/
     .       ' ', I4, ' BASIS VECTORS WILL BE USED AS A ',
     .            'STABILIZATION TRAIL.')
10716 FORMAT(' USING LOCALLY ORTHOGONAL GENCG SOLVER. ',
     .       I4, ' EQUATIONS WILL BE SOLVED AT ONCE ' )
10718 FORMAT(' USING LOCALLY ORTHOGONAL GENCG SOLVER. ',
     .            I4, ' EQUATIONS WILL BE SOLVED AT ONCE. '/
     .       ' ', I4, ' BASIS VECTORS WILL BE USED AS A ',
     .            'STABILIZATION TRAIL.')
C
11000 FORMAT(' ZERO IS AN ACCEPTABLE SOLUTION FOR RHS ', I5, 
     .       ', RETIRING IT WITH ERROR NORM ', G10.4 )
11100 FORMAT(' (Y,Y) VECTOR IS:'/ (5(' ',G14.7)) )
11200 FORMAT(' RESTART WITH ', I5, ' BASIS VECTORS' )
11300 FORMAT(' (Y,B) FOR RESTART BASIS IS:' )
11400 FORMAT(' (Y,R) FOR RESTART BASIS IS:' )
11500 FORMAT(' MINIMIZING RESIDUAL' )
11600 FORMAT(' SOLVING FOR ORTHOGONAL RESIDUAL' )
11700 FORMAT(' AUXILIARY LINEAR SYSTEM MATRIX IS:' )
11800 FORMAT(' AUXILIARY LINEAR SYSTEM RHS IS:' )
11900 FORMAT(' SOLUTION COEFFICIENTS ARE:' )
11910 FORMAT(' EFFECTIVE RANK OF THE PROBLEM IS ', I5 )
12000 FORMAT(' CURRENT SOLUTION FOR RHS ', I5, ' IS:' )
12100 FORMAT(' CURRENT RESIDUAL FOR RHS ', I5, ' IS:' )
12200 FORMAT(' RESIDUAL NORM FOR RHS ', I5, ' IS ', G14.7 )
12300 FORMAT(' SOLUTION FOR RHS ', I5, ' CONVERGED' )
12400 FORMAT(' ITERATIVE SOLVER FAILED TO CONVERGE IN ', I5, 
     .       ' ITERATIONS ON ', I5, ' VARIABLES OUT OF ', I5 )
12500 FORMAT(' SINGULAR VALUES OF THE OPTIMAL BASIS ARE:'/
     .       (5(' ',G14.7)) )
12600 FORMAT(' OPTIMAL SCRATCH SIZE FOR DGESVD IS ', I6, ' WORDS. ', 
     .       I6, ' WORDS WERE AVAILABLE.' )
12700 FORMAT(' OPTIMAL BASIS SET IS:' )
12800 FORMAT(' ', I5, ' BASIS VECTORS REJECTED SINCE SIGMA IS SMALL' )
12900 FORMAT(' ORTHOGONALIZATION AUXILIARY SYSTEM IS:' )
13000 FORMAT(' ORTHOGONALIZATION AUXILIARY SYSTEM RHS IS:' )
13090 FORMAT(' ORTHOGONALIZING COEFFICIENTS CLIPPED' )
13100 FORMAT(' ORTHOGONALIZING COEFFICIENTS ARE:' )
13200 FORMAT(' ORTHOGONALIZED NEW VECTORS ARE:' )
13300 FORMAT(' NORM OF NON-SCALED VECTOR ', I5, ' IS ', G14.7 )
13400 FORMAT(' ', I5, ' VECTORS VANISHED AFTER ORTHOGONALIZATION' )
13500 FORMAT(' UPDATED (B,B) IS:' )
13600 FORMAT(' UPDATED (Y,B) IS:' )
13700 FORMAT(' NEW RESPONSE VECTORS ARE:' )
13800 FORMAT(' UPDATED (R,B) IS:' )
13900 FORMAT(' UPDATED (R,R) IS:' )
14000 FORMAT(' UPDATED (Y,R) IS:' )
14050 FORMAT(' UPDATED (Y,Y) IS:' )
14100 FORMAT(' SOLUTION VECTORS ARE IN ORDER:'/(12(1X,I6)))
14200 FORMAT(' SORTING PERMUTATION TABLE IS:'/(12(1X,I6)))
      END
C
      SUBROUTINE PSLTRM(IVCNT,IBASSZ,IDELV,SCR,BAS,RSP,LDX,
     .                               YR,YB,LDY,BB,RB,RR,LDR)
C
C   Remove part of the basis vectors with the worst performance
C   scores.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IVCNT  - Number of valid rows in YR and YB
C      IBASSZ - Original size of the basis set.
C      IDELV  - Number of basis vectors to delete.
C      BAS    - Basis vectors.
C      RSP    - Response vectors
C      LDX    - Leading dimension of BAS and RSP
C      YR     - Dot products of RHS vectors with response vectors
C               First index - RHS, second index - response vector
C      YB     - Dot products of RHS vectors with basis vectors
C               First index - RHS, second index - basis vector
C      LDY    - Leading dimension of YR and YB
C      BB     - Dot products of basis vectors. 
C      RB     - Dot products of response vectors with basis vectors
C               First index - response vector, second index - basis 
C               vector. 
C      RR     - Dot products of response vectors
C      LDR    - Leading dimension of RR, RB and BB
C
C   Accessed common blocks:
C
C      PSDGB2 - "Response" computation parameters.
C      PSPRT  - Printing unit
C      PSPRTF - Debug output tuning flags
C
C   Modified common blocks:
C
C      PSDLNS - CPHF equations solution statistics
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (SONE=-1.0D0)
      PARAMETER (STWO=-2.0D0)
C
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
C
      COMMON
     ./PSDGB2/ NPAIR, NPAIR2, NOCCA, NVACA, NOCCB, NVACB,
     .         IQSZA, IQSZB,  IQSZ,  ICPV1, ICPVL, NCPVRS
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSDGB2/, /PSPRTF/, /PSPRT /
C
      DIMENSION BAS(LDX,*), RSP(LDX,*), YR(LDY,*), YB(LDY,*)
      DIMENSION RR(LDR,*), RB(LDR,*), BB(LDR,*), SCR(*)
C
C    Shortcut.
C
      IF(LPSOLV) THEN
          WRITE(NB6,11000) IDELV, IBASSZ
      ENDIF
      IF( IDELV.GE.IBASSZ ) RETURN
C
C    Select few smallest scores for deletion. We do not care
C    about the performance of the code below, since actual
C    compaction of vectors will be much more expensive.
C
      IREST = IBASSZ - IDELV
      DO 100 I=1,IDELV
          DO 70 J=1,IBASSZ
              IF( SCR(J).GE.SONE ) THEN
                  SCRMIN = SCR(J)
                  JMIN   = J
                  GOTO 80
              ENDIF
   70     CONTINUE
   80     CONTINUE
          DO 90 J=JMIN+1,IBASSZ
              IF( SCR(J).GE.SONE .AND. SCR(J).LT.SCRMIN ) THEN
                  SCRMIN = SCR(J)
                  JMIN   = J
              ENDIF
   90     CONTINUE
          IF(LPSOLV) THEN
              WRITE(NB6,11100) JMIN, SCR(JMIN)
          ENDIF
          SCR(JMIN) = STWO
  100 CONTINUE
C
C    Remove vectors marked for deletion, compressing columns of
C    dot product matrices in the progress.
C
      IOUT = 1
      DO 200 IIN=1,IBASSZ
          IF( SCR(IIN).GE.SONE ) THEN
              IF( IIN.NE.IOUT ) THEN
                  CALL DCOPY(IQSZ,BAS(1,IIN),1,BAS(1,IOUT),1)
                  CALL DCOPY(IQSZ,RSP(1,IIN),1,RSP(1,IOUT),1)
                  CALL DCOPY(IVCNT,YR(1,IIN),1,YR(1,IOUT),1)
                  CALL DCOPY(IVCNT,YB(1,IIN),1,YB(1,IOUT),1)
                  CALL DCOPY(IBASSZ,BB(1,IIN),1,BB(1,IOUT),1)
                  CALL DCOPY(IBASSZ,RB(1,IIN),1,RB(1,IOUT),1)
                  CALL DCOPY(IBASSZ,RR(1,IIN),1,RR(1,IOUT),1)
              ENDIF
              IOUT = IOUT + 1
          ENDIF
  200 CONTINUE
C
C    Compact columns of the dot product matrices and scores table.
C
      IOUT = 1
      DO 300 IIN=1,IBASSZ
          IF( SCR(IIN).GE.SONE ) THEN
              IF( IIN.NE.IOUT ) THEN
                  CALL DCOPY(IREST,BB(IIN,1),LDR,BB(IOUT,1),LDR)
                  CALL DCOPY(IREST,RB(IIN,1),LDR,RB(IOUT,1),LDR)
                  CALL DCOPY(IREST,RR(IIN,1),LDR,RR(IOUT,1),LDR)
                  SCR(IOUT) = SCR(IIN)
              ENDIF
              IOUT = IOUT + 1
          ENDIF
  300 CONTINUE
*     IF(LPSOLV) THEN
*         WRITE(NB6,11200)
*         CALL PSDPGM(IVCNT,IREST,YR,LDY)
*         WRITE(NB6,11300)
*         CALL PSDPGM(IVCNT,IREST,YB,LDY)
*         WRITE(NB6,11400)
*         CALL PSDPGM(IREST,IREST,BB,LDR)
*         WRITE(NB6,11500)
*         CALL PSDPGM(IREST,IREST,RB,LDR)
*         WRITE(NB6,11600)
*         CALL PSDPGM(IREST,IREST,RR,LDR)
*         IF(LPCPHF) THEN
*             WRITE(NB6,11700)
*             CALL PSDPGM(IQSZ,IREST,BAS,LDX)
*             WRITE(NB6,11800)
*             CALL PSDPGM(IQSZ,IREST,RSP,LDX)
*         ENDIF
*     ENDIF
      RETURN
11000 FORMAT(' DROPPING ', I5, ' BASIS VECTORS OUT OF ', I5 )
11100 FORMAT(' VECTOR ', I5, ' (SCORE ', G14.7,') MARKED FOR DELETION' )
11200 FORMAT(' (Y,R) AFTER COMPACTION IS:' )
11300 FORMAT(' (Y,B) AFTER COMPACTION IS:' )
11400 FORMAT(' (B,B) AFTER COMPACTION IS:' )
11500 FORMAT(' (R,B) AFTER COMPACTION IS:' )
11600 FORMAT(' (R,R) AFTER COMPACTION IS:' )
11700 FORMAT(' BASIS VECTORS AFTER COMPACTION ARE:' )
11800 FORMAT(' RESPONSE VECTORS AFTER COMPACTION ARE:' )
      END
