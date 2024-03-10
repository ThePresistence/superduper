C     ******************************************************************
C
C     Computation of right-hand side vectors for CPHF.
C
C     ******************************************************************
      SUBROUTINE PSDS1A(CALP,CBET,LDC,DUMP)
C
C   Compute right-hand sides of CPHF equations from the data
C   precomputed by PSDST1. This is N^4 step, but, fortunately,
C   it consitst entirely of BLAS calls and should be relatively
C   efficient. Alternative interface, PSDS1X might be used to
C   compute just one solution vector.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients.
C      CBET   - Beta orbital coefficients (Not used if UHF flag not set)
C      LDC    - Leading dimension of CALP, CBET matrices
C      DUMP   - Base of dynamic memory items
C
C      IEQ    - Number of variable to compute right-hand side for,
C               used by alternative entry point only.
C      VEC    - Place to store resulting vector to, used by alternative
C               entry point only.
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. Although some of the
C               options could be modified by strategy routines,
C               they are restored to the original values on return.
C      PSDYNM - Handles for the dynamic memory blocks
C      PSDGBL - Global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C      PSPRTF - Debug output tuning flags
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C      Single-variable entry, PSDS1X, is a kludge.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL UHF, HALFEL, DOCI, SINGLE, DORESP, DODIP, LIMAG, DOPTCH
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
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSPRTF/, /PSPRT /
C
      DIMENSION DUMP(*), CALP(LDC,*), CBET(LDC,*), VEC(*)
C
      IMINVR = ICPV1
      IMAXVR = ICPVL
      SINGLE = .FALSE.
      GOTO 1
          ENTRY PSDS1X(CALP,CBET,LDC,DUMP,IEQ,VEC)
          IMINVR = IEQ
          IMAXVR = IEQ
          SINGLE = .TRUE.
    1 CONTINUE
C
C    Get dynamic memory offsets
C
      IF( IMIX.EQ.3 ) THEN
          IDXPDA = IPSMOF(LDXPDA)
          IDXPAT = IPSMOF(LDXPAT)
          IEA    = IDXPAT
          IF(UHF) THEN
              IDXPBT = IPSMOF(LDXPBT)
              IEB    = IDXPBT
          ELSE
              IDXPBT = IDXPAT
              IEB    = IEA
          ENDIF
      ELSE
          IDXPA = IPSMOF(LDXPA)
          IF(UHF) IDXPB = IPSMOF(LDXPB)
      ENDIF
      IF( IQSWAP.LT.0 ) THEN
          IQ     = IPSMOF(LQ)
      ELSE
          IQTMP  = IPSMOF(LQTMP)
      ENDIF
      IF(LPCPHF) THEN
          WRITE(NB6,9980)
          CALL PSDPGM(NORBS,NORBS,CALP,LDC)
          IF(UHF) THEN
              WRITE(NB6,9990)
              CALL PSDPGM(NORBS,NORBS,CBET,LDC)
          ENDIF
      ENDIF
C
      DO 800 IX=IMINVR,IMAXVR
          IF( IX.EQ.0 ) THEN
C
C         Right-hand side for the Z-vector is a special case, and
C         is handled elsewhere.
C
          ELSE
              IF(IMIX.EQ.3) THEN
                  CALL PSDBMX(IX,DUMP(IDXPDA),DUMP(IDXPAT),DUMP(IDXPBT))
              ELSE
                  IEA = IDXPA + (IX-1)*NROS
                  IF(UHF) THEN
                      IEB = IDXPB + (IX-1)*NROS
                  ELSE
                      IEB = IEA
                  ENDIF
              ENDIF
              IF( IQSWAP.LT.0 ) THEN
                  IQTMP = IQ + (IX-ICPV1) * IQSZ
              ENDIF
              CALL PSXA2M(CALP,CBET,LDC,DUMP(IEA),DUMP(IEB),
     .                    DUMP(IQTMP),DUMP)
              IF(LPCPHF) THEN
                  WRITE(NB6,10000) IX
                  CALL PSXPRT(DUMP(IQTMP))
              ENDIF
          ENDIF
          IF(SINGLE) THEN
              CALL DCOPY(IQSZ,DUMP(IQTMP),1,VEC,1)
          ELSE IF( IQSWAP.GE.0 ) THEN
C
C             Swap out transformed vector
C
              CALL PSDWD(IURHS,IX+(1-ICPV1),DUMP(IQTMP),IQSZ)
          ENDIF
  800 CONTINUE
C
      RETURN
 9980 FORMAT(' ORBITAL COEFFICIENT MATRIX (ALPHA):')
 9990 FORMAT(' ORBITAL COEFFICIENT MATRIX (BETA):')
10000 FORMAT(' CPHF RIGHT-HAND VECTOR ', I4, ' IN MO BASIS:' )
      END
C
      SUBROUTINE PSNS1A(CALP,LDC,DUMP)
C
C   Compute magnetic field derivatives of one-electron hamiltonian
C   matrix and transform them into MO basis for use as right-hand
C   sides of CPHF equations.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients.
C      LDC    - Leading dimension of CALP matrix
C      DUMP   - Base of dynamic memory items
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options.
C      PSDYNM - Handles for the dynamic memory blocks
C      PSDGBL - Global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C      PSPRTF - Debug output tuning flags
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C        a0                                   a0
C      -H   is actually computed rather than H   itself, since our 
C      iterative solver expects it this way.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (SONE=-1.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
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
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSPRTF/, /PSPRT /
C
      DIMENSION DUMP(*), CALP(LDC,*)
C
C    Get dynamic memory offsets
C
      IF( IQSWAP.LT.0 ) THEN
          IQ     = IPSMOF(LQ)
      ELSE
          IQTMP  = IPSMOF(LQTMP)
      ENDIF
      IHA0 = IPSMOF(LHA0)
C
      CALL PSNHA0(DUMP(IHA0),NORBS,NORBS)
C
      DO 800 IX=ICPV1,ICPVL
          IF( IQSWAP.LT.0 ) THEN
              IQTMP = IQ + (IX-ICPV1) * IQSZ
          ENDIF
          IEA = IHA0 + (IX-ICPV1)*(NORBS**2)
          CALL PSNU2M(CALP,LDC,DUMP(IEA),NORBS,DUMP(IQTMP),DUMP)
C
C        Change sign of the right-hand side, since we are going to
C        solve AX + B = 0 problem.
C
          CALL DSCAL(IQSZ,SONE,DUMP(IQTMP),1)
          IF(LPCPHF) THEN
              WRITE(NB6,10000) IX
              CALL PSXPRT(DUMP(IQTMP))
          ENDIF
          IF( IQSWAP.GE.0 ) THEN
              CALL PSDWD(IURHS,IX+(1-ICPV1),DUMP(IQTMP),IQSZ)
          ENDIF
  800 CONTINUE
C
      RETURN
 9980 FORMAT(' ORBITAL COEFFICIENTS MATRIX (ALPHA):')
10000 FORMAT(' CPHF RIGHT-HAND VECTOR ', I4, ' IN MO BASIS:' )
      END
