C     ******************************************************************
C
C     Incorporate contributions from the density matrix derivatives
C     to the matrix of second derivatives.
C
C     ******************************************************************
      SUBROUTINE PSDST4(DERIV,FORCE,DIPDRV,LDF,ROALP,ROBET,DUMP,A,LDA,
     .                  CALP,CBET,LDC,EALP,EBET)
C
C   Combine contributions from density matrix dependence on molecular
C   coordinates into derivatives
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      DERIV  - First derivatives, adjusted in half-electron case.
C      FORCE  - Input/Output force constants. Not used if MODE.LT.2
C      DIPDRV - Input/Output derivatives of the dipole moment. Not used
C               if DODIP is not set.
C      LDF    - Leading dimension of FORCE array.
C      ROALP  - Alpha density matrix in packed format
C      ROBET  - Beta density matrix in packed format (not used if 
C               UHF flag not set)
C      DUMP   - Scratch array
C      A      - Array parameter which should be passed to SCFCAL 
C      LDA    - Dimension of A
C      CALP   - Alpha orbital coefficients
C      CBET   - Beta orbital coefficients, not used if UHF flag
C               not set.
C      LDC    - Leading dimension of CALP and CBET arrays.
C      EALP   - Orbital energies, might be needed by IDENS=8
C      EBET   - Beta orbital energies, might be needed by IDENS=8
C               then UHF is .TRUE.
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGBL - Global computation parameters
C      PSDYNM - Dynamic memory handles
C      PSPRT  - Printing unit
C      PSPRTF - Debug output tuning flags
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
C      PSDST4 computes derivatives of density matrices with respect
C      to each nuclear variable, either by finite difference (IDENS
C      .EQ.1) or from results of CPHF computation and combines them
C      with mixed coordinate-density derivatives of energy to form
C      indirect contributions to the second derivative.
C
C      Mixed derivatives can be passed in memory from previous
C      stages (IMIX=1,2) or recomputed for each new derivative of
C      density matrix from the stored contributions (IMIX=3).
C      The last possibility (IMIX=4) is to calculate derivatives
C      numerically, by finite-difference differentiation of the
C      Fock matrix.
C
C      Prior to summing, off-diagonal elements of density matrix
C      derivative are multiplied by factor of 2 to account for 
C      symmetry of density matrix.
C
CWT    There is an option to skip the calculation of the contribution
C      to the gradient for fixed MM atoms in the case of IMIX=4.
C      The general input option IN2(127) must be positive, and the
C      flag ISELCT(IAT) must be nonzero for a given MM atom.
C
C   Bugs:
C
C      When running in block mode, PSDST4 will not produce debugging
C      output on "response" contributions computed.
C
C      Treatment of first derivatives of the half-electron correction
C      is ugly. It is necessary if we do not want to do an extra
C      pass over Fock matrix derivatives, though.
C
      USE LIMIT, ONLY: LM1, LM1M
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
  200         CONTINUE
          ENDIF
          DO 800 IX=1,NVARS
              IC  = 1+MOD(IX-1,3)
              IAT = 1+(IX-IC)/3
              IF(IMIX.EQ.4) THEN
CWT             Skip calculation only for fixed MM atoms
                IF(IN2(127).LE.0 .OR. ISELCT(IAT).EQ.0) THEN
                  CALL PSOMFG(0,ONE,IAT,IC,A,LDA,ROALP,ROBET,
     .                        DUMP(IDXPAT),NROS,DUMMY,DUMMY)
                ENDIF
              ELSE IF(IMIX.EQ.3) THEN
                  CALL PSDBMX(IX,DUMP(IDXPDA),DUMP(IDXPAT),
     .                        DUMP(IDXPBT))
              ELSE
                  IEA = IDXPA + (IX-1)*NROS
                  IF(UHF) IEB = IDXPB + (IX-1)*NROS
              ENDIF
C
C             Contributions to second derivatives
C
              IF( MODE.GE.2 ) THEN
                  IF( IPTOP.EQ.IP ) THEN
                      XA = DDOT(NROS,DUMP(IEA),1,DUMP(IDPA),1)
                      FORCE(IX,IP) = FORCE(IX,IP) + XA
                      IF(LPSUMM) WRITE(NB6,11100) IX, XA
                      IF(UHF) THEN
                          XB = DDOT(NROS,DUMP(IEB),1,DUMP(IDPB),1)
                          FORCE(IX,IP) = FORCE(IX,IP) + XB
                          IF(LPSUMM) WRITE(NB6,11105) IX, XB
                      ENDIF
                  ELSE
                      CALL DGEMV('T',NROS,IPL,ONE,DUMP(IDPA),NROS,
     .                           DUMP(IEA),1,ONE,FORCE(IX,IP),LDF)
                      IF(UHF) THEN
                          CALL DGEMV('T',NROS,IPL,ONE,DUMP(IDPB),NROS,
     .                               DUMP(IEB),1,ONE,FORCE(IX,IP),LDF)
                      ENDIF
                  ENDIF
              ENDIF
C
C             Half-electron contributions to the first derivatives
C
              IF( (HALFEL.OR.DOCI) .AND. IP.EQ.1 ) THEN
CWT             Skip calculation only for fixed MM atoms
                IF(IN2(127).LE.0 .OR. ISELCT(IAT).EQ.0) THEN
                  XA = DDOT(NROS,DUMP(IEA),1,DUMP(IHLZA),1)
                  DERIV(IX) = DERIV(IX) + XA
                  IF(LPHALF) THEN
                      WRITE(NB6,11110) IX, XA
                  ENDIF
                ENDIF
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
C
      SUBROUTINE PSDBMX(IX,DXPDAT,DXPATM,DXPBTM)
C
C   Compute mixed density matrix-coordinate derivative for
C   the single coordinate from the array of precomputed data.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IX     - Coordinate to compute derivative by
C      DXPDAT - Precomputed pair derivatives
C      DXPATM - Mixed derivative on alpha density
C      DXPBTM - Mixed derivative on beta density (not used
C               if UHF flag not set)
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C      PSPRT  - Printing unit
C      PSPRTF - Debug output tuning flags
C      ATOMS  - Orbital indices
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      252 DOUBLE PRECISION cells.
C
C   Module logic:
C
C     See section of PSDST1 dealing with IMIX.EQ.1
C
C   Speedups possible:
C
C     Inlining calls to PSDAPA, PSDAPD, PSDSPA and PSDSPD with
C     simultaneous instantiation for NORBA,NORBB=1,4,9 might
C     help.
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      SAVE /PSDGBL/, /PSPRTF/, /PSPRT /
C
C*$*INLINE ROUTINE (PSDAPA,PSDAPD,PSDSPA,PSDSPD,PSU2UL)
CVD$R EXPAND (PSDAPA,PSDAPD,PSDSPA,PSDSPD,PSU2UL)
CVD$R SEARCH (psdst1.f)
C
      DIMENSION DXPDAT(*), DXPATM(*), DXPBTM(*)
      DIMENSION EAA(45), EBB(45), EABALP(9,9), EABBET(9,9)
C
      IATOM = (IX-1)/3+1
      IPTCH = IATOM - NATOM
      INDX  = MOD(IX-1,3)
      DO 100 I=1,NROS
          DXPATM(I) = ZERO
  100 CONTINUE
      IF(UHF) THEN
          DO 102 I=1,NROS
              DXPBTM(I) = ZERO
  102     CONTINUE
          N = 2
      ELSE
          N = 1
      ENDIF
      IIN = 0
C
      IF(DOPTCH) THEN
          MAXPT = NPTCHG
      ELSE
          MAXPT = MIN(1,NPTCHG)
          IF( IPTCH.GT.0 ) THEN
              WRITE(NB6,11100) IX, INDX, IPTCH
              STOP 'PSDBMX'
          ENDIF
      ENDIF
C
C-AK  DO 2000 NA=2,NATOM
      DO 2000 NA=1,NATOM
          NFRSTA = NFIRST(NA)
          NORBA  = NLAST(NA) - NFIRST(NA) + 1
          NPAIRA = ( NORBA * (NORBA+1) ) / 2
          DO 998 NB=1,NA-1
              NFRSTB = NFIRST(NB)
              NORBB  = NLAST(NB) - NFIRST(NB) + 1
              NPAIRB = ( NORBB * (NORBB+1) ) / 2
              NBLOCK = NPAIRA+NPAIRB+NORBA*NORBB*N
C
              IF( NA.EQ.IATOM .OR. NB.EQ.IATOM ) THEN
                  IIN = IIN + INDX*NBLOCK
                  DO 300 I=1,NPAIRA
                      EAA(I) = DXPDAT(IIN+I)
  300             CONTINUE
                  IIN = IIN + NPAIRA
                  DO 304 I=1,NPAIRB
                      EBB(I) = DXPDAT(IIN+I)
  304             CONTINUE
                  IIN = IIN + NPAIRB
                  DO 310 K=1,NORBB
                      DO 308 I=1,NORBA
                          EABALP(I,K) = DXPDAT(IIN+I)
  308                 CONTINUE
                      IIN = IIN + NORBA
  310             CONTINUE
                  IF(UHF) THEN
                      DO 316 K=1,NORBB
                          DO 314 I=1,NORBA
                              EABBET(I,K) = DXPDAT(IIN+I)
  314                     CONTINUE
                          IIN = IIN + NORBA
  316                 CONTINUE
                  ENDIF
                  IF( NA.EQ.IATOM ) THEN
                      CALL PSDAPD(DXPATM,NFRSTA,NORBA,EAA)
                      CALL PSDAPD(DXPATM,NFRSTB,NORBB,EBB)
                      CALL PSDAPA(DXPATM,NFRSTA,NFRSTB,NORBA,
     .                            NORBB,EABALP)
                      IF(UHF) THEN
                          CALL PSDAPD(DXPBTM,NFRSTA,NORBA,EAA)
                          CALL PSDAPD(DXPBTM,NFRSTB,NORBB,EBB)
                          CALL PSDAPA(DXPBTM,NFRSTA,NFRSTB,NORBA,
     .                                NORBB,EABBET)
                      ENDIF
                  ENDIF
                  IF( NB.EQ.IATOM ) THEN
                      CALL PSDSPD(DXPATM,NFRSTA,NORBA,EAA)
                      CALL PSDSPD(DXPATM,NFRSTB,NORBB,EBB)
                      CALL PSDSPA(DXPATM,NFRSTA,NFRSTB,NORBA,
     .                            NORBB,EABALP)
                      IF(UHF) THEN
                          CALL PSDSPD(DXPBTM,NFRSTA,NORBA,EAA)
                          CALL PSDSPD(DXPBTM,NFRSTB,NORBB,EBB)
                          CALL PSDSPA(DXPBTM,NFRSTA,NFRSTB,NORBA,
     .                                NORBB,EABBET)
                      ENDIF
                  ENDIF
                  IIN = IIN + (2-INDX)*NBLOCK
              ELSE
                  IIN = IIN + 3*NBLOCK
              ENDIF
  998     CONTINUE
C        Contributions from point charges, if any. All are of the same
C        size, which simplified our life a good deal.
          NBLOCK = NPAIRA
          DO 1900 NB=1,MAXPT
              IF( NA.EQ.IATOM .OR. NB.EQ.IPTCH ) THEN
                  IIN = IIN + INDX*NBLOCK
                  DO 1300 I=1,NPAIRA
                      EAA(I) = DXPDAT(IIN+I)
 1300             CONTINUE
                  IF( NA.EQ.IATOM ) THEN
                      CALL PSDAPD(DXPATM,NFRSTA,NORBA,EAA)
                      IF(UHF) THEN
                          CALL PSDAPD(DXPBTM,NFRSTA,NORBA,EAA)
                      ENDIF
                  ENDIF
                  IF( NB.EQ.IPTCH ) THEN
                      CALL PSDSPD(DXPATM,NFRSTA,NORBA,EAA)
                      IF(UHF) THEN
                          CALL PSDSPD(DXPBTM,NFRSTA,NORBA,EAA)
                      ENDIF
                  ENDIF
                  IIN = IIN + (3-INDX)*NBLOCK
              ELSE
                  IIN = IIN + 3*NBLOCK
              ENDIF
 1900     CONTINUE
 2000 CONTINUE
C
      IF(LPINTS) THEN
          WRITE(NB6,11000) IX
          CALL PSDPPM(NORBS,DXPATM)
          IF(UHF) THEN
              WRITE(NB6,11010) IX
              CALL PSDPPM(NORBS,DXPBTM)
          ENDIF
      ENDIF
C
      RETURN
11000 FORMAT(' RECOMPUTED ALPHA MIXED DERIVATIVES MATRIX ON ',
     .       'VARIABLE ', I4/)
11010 FORMAT(' RECOMPUTED BETA MIXED DERIVATIVES MATRIX ON ',
     .       'VARIABLE ', I4/)
11100 FORMAT(' CATASTROPHE IN PSDBMX: CALLED WITH IX = ',I5,
     .       ' (POINT CHARGE ',I5,', COMPONENT ',I1,
     .       '), BUT DOPTCH IS FALSE.')
      END
C
      SUBROUTINE PSCPHD(IP,IPOFF,DUMP,CALP,CBET,LDC,EALP,EBET)
C
C   Compute density matrix derivatives from CPHF solution vector
C   for the variable IP. Depending on options, this could actually
C   cause CPHF equations to be solved instead of fetching solution
C   vector from auxiliary storage.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IP     - Solution vector number
C      IPOFF  - Elements of DP array(s) to skip
C      DUMP   - Base of dynamic memory items
C      CALP   - Alpha orbital coefficients
C      CBET   - Beta orbital coefficients, not used if UHF flag
C               not set.
C      LDC    - Leading dimension of CALP and CBET arrays.
C      EALP   - Alpha orbital energies (used only if IDENS=8)
C      EBET   - Beta orbital energies (used only if IDENS=8 and
C               UHF is .TRUE.)
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGBL - Global computation parameters
C      PSDGB2 - "Response" computation parameters.
C      PSDYNM - Dynamic memory handles
C      PSPRT  - Printing unit
C      PSPRTF - Debug output tuning flags
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
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      DIMENSION CALP(LDC,*), CBET(LDC,*), DUMP(*), EALP(*), EBET(*)
C
C    If necessary, solve single CPHF equation
C
      IF( IDENS.EQ.8 ) THEN
          IF( ISOLVE.EQ.1 ) THEN
              CALL PSDS3O(CALP,CBET,LDC,EALP,EBET,DUMP,IP)
          ELSE
              CALL PSLIN1(CALP,CBET,LDC,EALP,EBET,DUMP,IP)
          ENDIF
      ENDIF
C
C    Get the CPHF solution vector
C
      IF( IQSWAP.EQ.3 ) THEN
          IXTMP = IPSMOF(LXTMP)
          IALP  = IXTMP
      ELSE IF( IQSWAP.EQ.2 ) THEN
          IXTMP = IPSMOF(LXTMP)
          CALL PSDRD(IURHS,IP+(1-ICPV1),DUMP(IXTMP),IQSZ)
          IALP  = IXTMP
      ELSE
          IQ    = IPSMOF(LQ)
          IALP  = IQ + (IP-ICPV1)*IQSZ
      ENDIF
      IDPA = IPSMOF(LDPA) + IPOFF
      IF(UHF) THEN
          IDPB = IPSMOF(LDPB) + IPOFF
      ELSE
          IDPB = IDPA
      ENDIF
C
C    Check solutions for very large coupling coefficiients
C
      IF( IPRINT.GE.0 ) THEN
          CALL PSXCHB(IP,DUMP(IALP))
      ENDIF
C
C    Convert CPHF solution into derivatives of the density matrices
C
      CALL PSXM2P(CALP,CBET,LDC,DUMP(IALP),DUMP(IDPA),DUMP(IDPB),DUMP)
C
      IF(LPNUME) THEN
          WRITE(NB6,11000) IP
          CALL PSDPPM(NORBS,DUMP(IDPA))
          IF(UHF) THEN
              WRITE(NB6,11010) IP
              CALL PSDPPM(NORBS,DUMP(IDPB))
          ENDIF
      ENDIF
C
      RETURN
11000 FORMAT(' ALPHA DENSITY MATRIX DERIVATIVE ON VARIABLE ', I5 /)
11010 FORMAT(' BETA DENSITY MATRIX DERIVATIVE ON VARIABLE ', I5 /)
      END
C
