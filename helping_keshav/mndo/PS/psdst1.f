C     ******************************************************************
C
C     Compute intermediate two-electron derivative terms.
C
C     ******************************************************************
      SUBROUTINE PSDST1(ENERGY,DERIV,FORCE,LDF,ROALP,ROBET,
     .                  CALP,LDC,DUMP)
C
C   Compute partial cartesian energy derivatives for MNDO/d-type
C   methods and store intermediate results for computation of
C   density-dependent contributions.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ENERGY - Part of the total energy which depends on nuclear
C               coordinates.
C      DERIV  - Output first derivatives. For RHF and UHF case
C               they are completely determined on stage 1. In
C               half-electron method, additional terms have to
C               be computed later.
C      FORCE  - Output force constants, should be at least 3*NATOMS
C               in second dimension. Not accessed if IMOD.LE.1
C      LDF    - Leading dimension of FORCE array, should be at least
C               3*NATOMS. Ignored if IMOD.LE.1
C      ROALP  - Alpha density matrix in packed format
C      ROBET  - Beta density matrix in packed format (not used if 
C               UHF flag not set)
C      CALP   - Alpha orbital coefficients positioned over the
C               first active MO, not used unless one of HALFEL or 
C               DOCI is set and IHLST is not 2.
C      LDC    - Leading dimension of CALP matrix.
C      DUMP   - Scratch array. 
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGBL - Global computation parameters
C      PSEXRO - Mixed derivatives with respect to coordinate/density
C      PSDYNM - Dynamic memory handles
C      PSPRT  - Printing unit.
C      PSPRTF - Debug output tuning flags
C      PSHALF - Half-electron parameters
C      ATOMS  - Starting orbital number, number of orbitals
C               and atom types.
C      INOPT2 - Input options.
C      PARDER - Core charges
C      QMM1   - Point charges and their positions
C      QMM2   - Indicates "link" atoms, which do not interact with
C               the MM part of the system.
C
C   Modified common blocks:
C
C      PSPAIR - Pair parameters
C      PSDENS - Density matrices
C
C   Local storage:
C
C      22030 DOUBLE PRECISION cells for two-electron repulsion
C      integrals, overlap integrals, their derivatives and storage
C      of pair force contributions.
C
C      90*MAXACT + 1620 (5670 in the present version) DOUBLE PRECISION
C      cells are used for quantities related to the half-electron
C      derivatives.
C
C   Module logic:
C
C   Bugs:
C
C      Arrays in COMMON /PSDENS/ should be exactly in
C      order ALP, BET, {sum} - PSDPAA and PSDPAB expect
C      to find them at these places.
C
C      This is a *HUGE* function. Splitting it into several independent
C      fragments might be a good idea. Unfortunately, it is integral-
C      driven, and has a great deal of overlap between code used in
C      different modes.
C
CWT    When using this routine to compute the CI contribution to
C      nonadiabatic couplings (see further comments below), the
C      calculation of the SCF contribution could be skipped.
C
      USE LIMIT, ONLY: LM1, LMZ, LM1M, LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
      PARAMETER (MAXCIO=LMACT)
      PARAMETER (EV =27.2100D0)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (ONE =1.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      LOGICAL ORBDER
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
     ./PSPAIR/ ZCOREA, ZCOREB,
     .         NA, NB, NTYPA, NTYPB, NFRSTA, NORBA, NORBB, 
     .         NFRSTB, NPAIRA, NPAIRB, N1CNTA
     ./PSDENS/ PAAALP(45), PAABET(45), ZAA, PAA(45),
     .         PBBALP(45), PBBBET(45), ZBB, PBB(45),
     .         PABALP(9,9), PABBET(9,9), PAB(9,9)
     ./PSEXRO/ EAA(45,3), EBB(45,3), EABALP(9,9,3), EABBET(9,9,3)
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
     ./PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./QMMM1 / COORDM(3,LM1M),CHARGM(LM1M)
     ./QMMM2 / LINK(LM1)
C
      DIMENSION DERIV(*), FORCE(LDF,*), DUMP(*), ROALP(*), ROBET(*)
      DIMENSION CALP(LDC,*)
C
C    Magic to force ftnchek to shut up :-0)
C
      DIMENSION PAAX(136), PBBX(136), PABX(243)
      EQUIVALENCE (PAAX(1),PAAALP(1)),(PBBX(1),PBBALP(1))
      EQUIVALENCE (PABX(1),PABALP(1,1))
C
C    Local temporary storage
C
      DIMENSION AINTS(46,46,10), AOVER(9,9,10)
      DIMENSION FPSUM(10), FPALP(10), FPBET(10), FCORE(10), FSUM(10)
      DIMENSION FHALF(10), XJ(9,9,10), XK(9,9,10)
C
C    Some very local constants
C        I1STX is the index of the first derivative to process
C        ILASTX is the index of the last derivative to process
C    Overlap integrals are special: even if we are not going
C    to compute heat of formation, we could still need scaled
C    overlap integrals for extrapolation of H matrix in numeric
C    density derivatives computation. Hence,
C        I1STO and
C        ILASTO
C
      IF( MODE.LE.1 ) THEN
          ILASTX =  4
      ELSE
          ILASTX = 10
      ENDIF
      IF( IENRG.EQ.1 ) THEN
          I1STX  =  1
      ELSE
          I1STX  =  2
      ENDIF
      I1STO  = I1STX
      ILASTO = ILASTX
      IF( MODE.GE.2 .AND. IDENS.EQ.1 ) THEN
          I1STO = 1
      ENDIF
C
C    Some "pointers" for response-related temporaries
C
      IF(DORESP) THEN
          IF( IMIX.EQ.3 ) THEN
              IDXPDA = IPSMOF(LDXPDA) - 1
          ELSE
              IDXPA = IPSMOF(LDXPA) 
              CALL PSZRM(NROS,NVARS,DUMP(IDXPA),NROS)
              IDXPA = IDXPA - 4 * NROS
              IF( UHF ) THEN
                  IDXPB = IPSMOF(LDXPB) 
                  CALL PSZRM(NROS,NVARS,DUMP(IDXPB),NROS)
                  IDXPB = IDXPB - 4 * NROS
              ENDIF
          ENDIF
          IF( IDENS.EQ.1 ) THEN
              IHIDAT = IPSMOF(LHIDAT) - 1
          ELSE
              IAI2T1 = IPSMOF(LAI2T1) - 1
          ENDIF
      ENDIF
C
C    Half-electron and CI- data.
C
      ORBDER = HALFEL.OR.DOCI
      IF(ORBDER) THEN
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
              IHLOAW = IHLOA
              IHLOB  = IHLOA + NOIJ*45
              IHLOBW = IHLOB
          ENDIF
      ENDIF
C
C    Loop over all atom pairs
C
      DO 9000 NA=1,NATOM
          NFRSTA = NFIRST(NA)
          NORBA  = NLAST(NA) - NFIRST(NA) + 1
          NPAIRA = ( NORBA * (NORBA+1) ) / 2
          NTYPA  = NAT(NA)
          ZCOREA = TORE( NTYPA )
          IA     = NA*3-2
          IF( NORBA.EQ.1 ) THEN
              N1CNTA =   1
          ELSE IF( NORBA.EQ.4 ) THEN
              N1CNTA =  22
          ELSE IF( NORBA.EQ.9 ) THEN
              N1CNTA = 265
          ELSE
              WRITE(NB6,10000) NA, NORBA
              STOP 'PSDST1'
          ENDIF
C
C         Fetch density matrix fragment around atom A to the
C         local array
C
          CALL PSDPAA(ROALP,ROBET,NFRSTA,NORBA,PAAX)
C
C         Fetch (or recompute) contracted orbital coefficients needed
C         for half-electron case.
C
          IF(ORBDER) THEN
              IF( IHLST.EQ.2 ) THEN
                  IHLOAW = IHLOA
                  IHLOA  = IHLOA + NOIJ*NPAIRA
                  IHLOB  = IHLO
              ELSE
                  CALL PSHSGO(NA,NORBA,NACT,CALP(NFRSTA,1),LDC,
     .                        DUMP(IHLOA))
              ENDIF
          ENDIF
          DO 2090 NB=1,NA-1
              NFRSTB = NFIRST(NB)
              NORBB  = NLAST(NB) - NFIRST(NB) + 1
              NPAIRB = ( NORBB * (NORBB+1) ) / 2
              NTYPB  = NAT(NB)
              ZCOREB = TORE( NTYPB )
              IB     = NB*3-2
              IF( NORBB.NE.1 .AND. NORBB.NE.4 .AND. NORBB.NE.9 ) THEN
                  WRITE(NB6,10000) NB, NORBB
                  STOP 'PSDST1'
              ENDIF
              IF(LPINTS) THEN
                   WRITE(NB6,10100) NA, NB
                   WRITE(NB6,10110) ZCOREA, ZCOREB, NA, NB, NTYPA, 
     .                            NTYPB, NFRSTA, NORBA, NORBB, NFRSTB,
     .                            NPAIRA, NPAIRB, N1CNTA
              ENDIF
C
C             Fetch density matrix block around A and cross-block
C             for A-B pair.
C
              CALL PSDPAA(ROALP,ROBET,NFRSTB,NORBB,PBBX)
              CALL PSDPAB(ROALP,ROBET,NFRSTA,NFRSTB,NORBA,NORBB,PABX)
C
C             Fetch (or recompute) contracted orbital coefficients needed
C             for half-electron case.
C
              IF(ORBDER) THEN
                  IF( IHLST.EQ.2 ) THEN
                      IHLOBW = IHLOB
                      IHLOB  = IHLOB + NOIJ*NPAIRB
                  ELSE
                      CALL PSHSGO(NB,NORBB,NACT,CALP(NFRSTB,1),LDC,
     .                            DUMP(IHLOB))
                  ENDIF
              ENDIF
C
C             Compute necessary two-center integrals. Scale overlap
C             integrals by corresponding betas.
C
              CALL PSINTS(MODE,NA,NB,AINTS,AINTS(1,1,2),AINTS(1,1,5))
              CALL PSOVER(MODE,NA,NB,AOVER)
              DO 510 I=I1STO,ILASTO
                  CALL PSDABT(AOVER(1,1,I))
  510         CONTINUE
C
C             Compute pair energy, pair force and pair force constants
C
              DO 610 I=I1STX,ILASTX
                  CALL PSDFPS(AINTS(1,1,I),AOVER(1,1,I),FPSUM(I))
                  CALL PSDFPX(AINTS(1,1,I),PABALP,FPALP(I))
  610         CONTINUE
              IF( UHF ) THEN
                  DO 620 I=I1STX,ILASTX
                      CALL PSDFPX(AINTS(1,1,I),PABBET,FPBET(I))
  620             CONTINUE
              ELSE
                  DO 630 I=I1STX,ILASTX
                      FPBET(I) = FPALP(I)
  630             CONTINUE
              ENDIF
C
C             Compute core repulsion term and it's derivatives, then
C             sum all contributions to pair force.
C
              CALL PSDFCR(AINTS,FCORE)
              DO 700 I=I1STX,ILASTX
                  FSUM(I) = FPSUM(I) + FPALP(I) + FPBET(I) + FCORE(I)
  700         CONTINUE
C
C            Half-electron contribution
C
              IF(HALFEL) THEN
                  DO 750 I=I1STX,ILASTX
                      CALL PSHSJK(NPAIRA,NPAIRB,AINTS(2,2,I),46,
     .                            DUMP(IHLOAW),DUMP(IHLOBW),
     .                            XJ(1,1,I),XK(1,1,I),9)
                      FHALF(I) = ZERO
                      DO 730 K=1,NUMOPN
                          DO 720 J=1,K-1
                              FHALF(I) = FHALF(I) + HLFH(J,K)*XJ(J,K,I) 
     .                                            + HLFG(J,K)*XK(J,K,I)
  720                     CONTINUE
                          FHALF(I) = FHALF(I) + HLFH(K,K)*XJ(K,K,I)
  730                 CONTINUE
                      IF( IENRG.NE.4321 ) FSUM(I) = FSUM(I) + FHALF(I)
  750             CONTINUE
                  IF(LPHALF.AND.LPINTS) THEN
                      WRITE(NB6,10390) NA, NB
                      DO 770 I=I1STX,ILASTX
                          WRITE(NB6,10394) I
                          CALL PSDPSU(NUMOPN,XJ(1,1,I),9)
                          WRITE(NB6,10396) I
                          CALL PSDPSU(NUMOPN,XK(1,1,I),9)
  770                 CONTINUE
                  ENDIF
              ENDIF
C
C            CI contribution. Derivatives of all ERIs over MOs in active
C            space have to be computed here.
C
CWT          CI energy derivatives require the two-electron density
C            matrix of a given CI state as input in DUMP(ICIGAM).
C            The CI contribution to nonadiabatic couplings between two
C            CI states can be computed by providing the corresponding
C            two-electron transition density as input in DUMP(ICIGAM).
C
              IF(DOCI) THEN
                  DO 780 I=I1STX,ILASTX
                      CALL PSZRMB(NPAIRA,NOIJ,DUMP(ICI2TM),NPAIRA)
                      CALL PSZRMB(NOIJ,NOIJ,DUMP(ICI2ES),NOIJ)
                      CALL DGEMM('N','N',NPAIRA,NOIJ,NPAIRB,ONE,
     .                           AINTS(2,2,I),46,DUMP(IHLOBW),NPAIRB,
     .                           ZERO,DUMP(ICI2TM),NPAIRA)
                      CALL DSYR2K('U','T',NOIJ,NPAIRA,ONE,DUMP(IHLOAW),
     .                            NPAIRA,DUMP(ICI2TM),NPAIRA,ZERO,
     .                            DUMP(ICI2ES),NOIJ)
                      CALL PSU2PS(NOIJ,DUMP(ICI2ES),NOIJ)
                      IF(LPCI.AND.LPINTS) THEN
                          WRITE(NB6,10398) I
                          CALL PSDPPM(NOIJ,DUMP(ICI2ES))
                      ENDIF
                      CALL PSSCTN(NACT,DUMP(ICI2ES))
                      FHALF(I) = DDOT(NCIGAM,DUMP(ICI2ES),1,
     .                                       DUMP(ICIGAM),1)
                      IF( IENRG.NE.4321 ) FSUM(I) = FSUM(I) + FHALF(I)
  780             CONTINUE
              ENDIF
C
C            Debugging output
C
              IF(LPINTS) THEN
                  WRITE(NB6,10400)
                  CALL PSDPGM(NPAIRA,NPAIRB,AINTS(2,2,1),46)
                  DO 790 I=2,5
                      WRITE(NB6,10405) I - 1
                      CALL PSDPGM(NPAIRA,NPAIRB,AINTS(2,2,I),46)
  790             CONTINUE
                  WRITE(NB6,10410)
                  CALL PSDPGM(NORBA,NORBB,AOVER,9)
                  WRITE(NB6,10300)
                  CALL PSDPPM(NORBA,PAAALP)
                  IF( UHF ) THEN
                      WRITE(NB6,10310)
                      CALL PSDPPM(NORBA,PAABET)
                  ENDIF
                  WRITE(NB6,10320)
                  CALL PSDPPM(NORBA,PAA)
                  WRITE(NB6,10330)
                  CALL PSDPPM(NORBB,PBBALP)
                  IF( UHF ) THEN
                      WRITE(NB6,10340)
                      CALL PSDPPM(NORBB,PBBBET)
                  ENDIF
                  WRITE(NB6,10350)
                  CALL PSDPPM(NORBB,PBB)
                  WRITE(NB6,10360)
                  CALL PSDPGM(NORBA,NORBB,PABALP,9)
                  IF( UHF ) THEN
                      WRITE(NB6,10370)
                      CALL PSDPGM(NORBA,NORBB,PABBET,9)
                  ENDIF
                  WRITE(NB6,10380)
                  CALL PSDPGM(NORBA,NORBB,PAB,9)
                  WRITE(NB6,10200)
                  CALL PSDPDR(10,FPSUM)
                  WRITE(NB6,10210)
                  CALL PSDPDR(10,FPALP)
                  IF( UHF ) THEN
                      WRITE(NB6,10220)
                      CALL PSDPDR(10,FPBET)
                  ENDIF
                  WRITE(NB6,10230)
                  CALL PSDPDR(10,FCORE)
                  IF(HALFEL) THEN
                      WRITE(NB6,10235)
                      CALL PSDPDR(10,FHALF)
                  ENDIF
                  WRITE(NB6,10240)
                  CALL PSDPDR(10,FSUM)
              ENDIF
C
CWT           Check for computation of nonadiabiatic couplings (NAC).
C             Only the CI contributions (FHALF) are needed for NAC.
C             Static part computed here, response computed in psdst4.
C             See equation 8 in Theor. Chem. Acc. 2005, 114, xxx.
C
              NACTSK = IN2(160)
              IF(NACTSK.GT.10) THEN
                 FSUM(I1STX:ILASTX) = FHALF(I1STX:ILASTX)
              ENDIF
C
C             Add pair contributions to corresponding derivatives
C
              IF( IENRG.EQ.1 ) THEN
                  ENERGY = ENERGY + FSUM(1)
              ENDIF
              IF( MODE.GE.1 ) THEN
                  DO 800 I=0,2
                      DERIV(IA+I) = DERIV(IA+I) - FSUM(2+I)
                      DERIV(IB+I) = DERIV(IB+I) + FSUM(2+I)
  800             CONTINUE
              ENDIF
              IF( MODE.GE.2 ) THEN
                  DO 810 I=0,2
                      FORCE(IA+I,IA+I)=FORCE(IA+I,IA+I)+FSUM(5+I)
                      FORCE(IB+I,IB+I)=FORCE(IB+I,IB+I)+FSUM(5+I)
                      FORCE(IA+I,IB+I)=FORCE(IA+I,IB+I)-FSUM(5+I)
                      FORCE(IB+I,IA+I)=FORCE(IB+I,IA+I)-FSUM(5+I)
  810             CONTINUE
                  DO 820 I=0,2
                      I1 = MAX(I-1,0)
                      I2 = MIN(I+1,2)
                      FORCE(IA+I1,IA+I2)=FORCE(IA+I1,IA+I2)+FSUM(8+I)
                      FORCE(IA+I2,IA+I1)=FORCE(IA+I2,IA+I1)+FSUM(8+I)
                      FORCE(IB+I1,IB+I2)=FORCE(IB+I1,IB+I2)+FSUM(8+I)
                      FORCE(IB+I2,IB+I1)=FORCE(IB+I2,IB+I1)+FSUM(8+I)
                      FORCE(IA+I1,IB+I2)=FORCE(IA+I1,IB+I2)-FSUM(8+I)
                      FORCE(IA+I2,IB+I1)=FORCE(IA+I2,IB+I1)-FSUM(8+I)
                      FORCE(IB+I1,IA+I2)=FORCE(IB+I1,IA+I2)-FSUM(8+I)
                      FORCE(IB+I2,IA+I1)=FORCE(IB+I2,IA+I1)-FSUM(8+I)
  820             CONTINUE
              ENDIF
C
C             Do whichever processing is needed for this atom pair
C             during phase 1, and pick next pair.
C
              IF(DORESP) THEN
C
C                 Compute mixed cartesian-density matrix derivatives.
C
                  CALL PSDERO(AINTS,AOVER)
                  IF( IMIX.EQ.3 ) THEN
C
C                     Store computed components of mixed derivatives
C
                      DO 918 J=1,3
                          DO 900 I=1,NPAIRA
                              DUMP(IDXPDA+I) = EAA(I,J)
  900                     CONTINUE
                          IDXPDA = IDXPDA + NPAIRA
                          DO 904 I=1,NPAIRB
                              DUMP(IDXPDA+I) = EBB(I,J)
  904                     CONTINUE
                          IDXPDA = IDXPDA + NPAIRB
                          DO 910 K=1,NORBB
                              DO 908 I=1,NORBA
                                  DUMP(IDXPDA+I) = EABALP(I,K,J)
  908                         CONTINUE
                              IDXPDA = IDXPDA + NORBA
  910                     CONTINUE
                          IF( UHF ) THEN
                              DO 916 K=1,NORBB
                                  DO 914 I=1,NORBA
                                      DUMP(IDXPDA+I) = EABBET(I,K,J)
  914                             CONTINUE
                                  IDXPDA = IDXPDA + NORBA
  916                         CONTINUE
                          ENDIF
  918                 CONTINUE
                  ELSE
C
C                     Combine computed components into complete matrix
C
                      DO 930 I=1,3
                          J = IDXPA+(3*NA+I)*NROS
                          CALL PSDAPD(DUMP(J),NFRSTA,NORBA,EAA(1,I))
                          CALL PSDAPD(DUMP(J),NFRSTB,NORBB,EBB(1,I))
                          CALL PSDAPA(DUMP(J),NFRSTA,NFRSTB,NORBA,
     .                                NORBB,EABALP(1,1,I))
                          J = IDXPA+(3*NB+I)*NROS
                          CALL PSDSPD(DUMP(J),NFRSTA,NORBA,EAA(1,I))
                          CALL PSDSPD(DUMP(J),NFRSTB,NORBB,EBB(1,I))
                          CALL PSDSPA(DUMP(J),NFRSTA,NFRSTB,NORBA,
     .                                NORBB,EABALP(1,1,I))
                          IF( UHF ) THEN
                              J = IDXPB+(3*NA+I)*NROS
                              CALL PSDAPD(DUMP(J),NFRSTA,NORBA,EAA(1,I))
                              CALL PSDAPD(DUMP(J),NFRSTB,NORBB,EBB(1,I))
                              CALL PSDAPA(DUMP(J),NFRSTA,NFRSTB,NORBA,
     .                                    NORBB,EABBET(1,1,I))
                              J = IDXPB+(3*NB+I)*NROS
                              CALL PSDSPD(DUMP(J),NFRSTA,NORBA,EAA(1,I))
                              CALL PSDSPD(DUMP(J),NFRSTB,NORBB,EBB(1,I))
                              CALL PSDSPA(DUMP(J),NFRSTA,NFRSTB,NORBA,
     .                                    NORBB, EABBET(1,1,I))
                          ENDIF
  930                 CONTINUE
                  ENDIF
                  IF( IDENS.EQ.1 ) THEN
C
C                     Store integrals and derivatives needed to extrapolate
C                     integrals and hamiltonian after small change in nuclear
C                     coordinates (needed by IDENS=1)
C
                      DO 970 K=1,4
                          DO 954 J=1,NPAIRB
                               DO 952 I=1,NPAIRA
                                   DUMP(IHIDAT+I) = EV*AINTS(I+1,J+1,K)
  952                          CONTINUE
                               IHIDAT = IHIDAT + NPAIRA
  954                     CONTINUE
                          DO 958 I=1,NPAIRA
                              DUMP(IHIDAT+I) = EV * AINTS(I+1,1,K) * 
     .                                         (-ZCOREB)
  958                     CONTINUE
                          IHIDAT = IHIDAT + NPAIRA
                          DO 962 J=1,NPAIRB
                              DUMP(IHIDAT+J) = EV * AINTS(1,J+1,K) *
     .                                         (-ZCOREA)
  962                     CONTINUE
                          IHIDAT = IHIDAT + NPAIRB
                          DO 968 J=1,NORBB
                              DO 966 I=1,NORBA
                                  DUMP(IHIDAT+I) = EV * AOVER(I,J,K)
  966                         CONTINUE
                              IHIDAT = IHIDAT + NORBA
  968                     CONTINUE
  970                 CONTINUE
                  ELSE
C
C                     Store two-center integrals needed by CPHF
C
                      DO 1010 J=1,NPAIRB
                          DO 1008 I=1,NPAIRA
                              DUMP(IAI2T1+I) = AINTS(1+I,1+J,1)
 1008                     CONTINUE
                          IAI2T1 = IAI2T1 + NPAIRA
 1010                 CONTINUE
C
                  ENDIF
                  IF(LPINTS) THEN
                      WRITE(NB6,10500)
                      CALL PSDPGM(NPAIRA,3,EAA,45)
                      WRITE(NB6,10510)
                      CALL PSDPGM(NPAIRB,3,EBB,45)
                      WRITE(NB6,10520)
                      CALL PSDPGM(81,3,EABALP,81)
                      IF( UHF ) THEN
                          WRITE(NB6,10530)
                          CALL PSDPGM(81,3,EABBET,81)
                      ENDIF
                  ENDIF
              ENDIF
 2090     CONTINUE
          IF( DORESP .AND. IDENS.GE.3 ) THEN
C
C         Compute and store one-center integrals if necessary for CPHF
C
              CALL PSINT1(NTYPA,NPAIRA,AINTS)
              DO 2098 J=1,NPAIRA
                   DO 2096 I=1,NPAIRA
                       DUMP(IAI2T1+I) = AINTS(1+I,1+J,1)
 2096              CONTINUE
                   IAI2T1 = IAI2T1 + NPAIRA
 2098         CONTINUE
              IF(LPINTS) THEN
                  WRITE(NB6,10600) NA
                  CALL PSDPGM(NPAIRA,NPAIRA,AINTS(2,2,1),46)
              ENDIF
          ENDIF
C        Contributions from point charges, only applicable if THERE ARE
C        some point charges :-)
          IF( NPTCHG.EQ.0 ) GOTO 9000
          IF( DORESP .AND. IMIX.EQ.3 ) THEN
C            Clear workspace for accumulation of pair contributions to Fock
C            matrix derivatives.
              LEN = NPAIRA*3
              IF(DOPTCH) LEN = LEN * NPTCHG
              DUMP(IDXPDA+1) = ZERO
              CALL DCOPY(LEN-1,DUMP(IDXPDA+1),0,DUMP(IDXPDA+2),1)
              ITXPDA = IDXPDA
              IDXPDA = IDXPDA + LEN
          ENDIF
          MMLINK = IN2(123)
          IF( MMLINK.NE.1 .OR. LINK(NA).LE.0 ) THEN
C           This is not a link atom, so it feels all point charges all right.
C           All centers in the point charge loop have the same properties,
C           except for charge. The 4090 loop is actually a cloned 2090,
C           with all unncessary things cut out.
             NTYPB  = 0
             NORBB  = 0
             NFRSTB = 0
             NPAIRB = 0
             DO 4090 NB=1,NPTCHG
                ZCOREB = CHARGM(NB)
                IB     = (NATOM+NB)*3-2
                IF(LPINTS) THEN
                     WRITE(NB6,10102) NA, NB
                     WRITE(NB6,10110) ZCOREA, ZCOREB, NA, NB, NTYPA, 
     .                              NTYPB, NFRSTA, NORBA, NORBB, NFRSTB, 
     .                              NPAIRA, NPAIRB, N1CNTA
                ENDIF
C              Compute core attraction integrals and their derivatives.
C              Negative NB parameter indicates point charge.
                CALL PSINTS(MODE,NA,-NB,AINTS,AINTS(1,1,2),AINTS(1,1,5))
C              Compute pair energy, pair force and pair force constants
C              Point charges do not affect two-electron terms, so there 
C              are no spin-dependent components in this case.
                DO 2610 I=I1STX,ILASTX
                    FPSUM(I) = -ZCOREB*DDOT(NPAIRA,PAA,1,AINTS(2,1,I),1)
 2610           CONTINUE
C              Compute core repulsion term and it's derivatives, then
C              sum all contributions to pair force.
                CALL PSDFCR(AINTS,FCORE)
                DO 2700 I=I1STX,ILASTX
                    FSUM(I) = FPSUM(I) + FCORE(I)
 2700           CONTINUE
C              Debugging output
                IF(LPINTS) THEN
                    WRITE(NB6,10200)
                    CALL PSDPDR(10,FPSUM)
                    WRITE(NB6,10230)
                    CALL PSDPDR(10,FCORE)
                    WRITE(NB6,10240)
                    CALL PSDPDR(10,FSUM)
                ENDIF
C              Add pair contributions to corresponding derivatives
                IF( IENRG.EQ.1 ) THEN
                    ENERGY = ENERGY + FSUM(1)
                ENDIF
                IF( MODE.GE.1 ) THEN
                    DO 2800 I=0,2
                        DERIV(IA+I) = DERIV(IA+I) - FSUM(2+I)
                        IF( .NOT.DOPTCH ) GOTO 2800
                        DERIV(IB+I) = DERIV(IB+I) + FSUM(2+I)
 2800               CONTINUE
                ENDIF
                IF( MODE.GE.2 ) THEN
                    DO 2810 I=0,2
                        FORCE(IA+I,IA+I)=FORCE(IA+I,IA+I)+FSUM(5+I)
                        IF( .NOT.DOPTCH ) GOTO 2810
                        FORCE(IB+I,IB+I)=FORCE(IB+I,IB+I)+FSUM(5+I)
                        FORCE(IA+I,IB+I)=FORCE(IA+I,IB+I)-FSUM(5+I)
                        FORCE(IB+I,IA+I)=FORCE(IB+I,IA+I)-FSUM(5+I)
 2810               CONTINUE
                    DO 2820 I=0,2
                       I1 = MAX(I-1,0)
                       I2 = MIN(I+1,2)
                       FORCE(IA+I1,IA+I2)=FORCE(IA+I1,IA+I2)+FSUM(8+I)
                       FORCE(IA+I2,IA+I1)=FORCE(IA+I2,IA+I1)+FSUM(8+I)
                       IF( .NOT.DOPTCH ) GOTO 2820
                       FORCE(IB+I1,IB+I2)=FORCE(IB+I1,IB+I2)+FSUM(8+I)
                       FORCE(IB+I2,IB+I1)=FORCE(IB+I2,IB+I1)+FSUM(8+I)
                       FORCE(IA+I1,IB+I2)=FORCE(IA+I1,IB+I2)-FSUM(8+I)
                       FORCE(IA+I2,IB+I1)=FORCE(IA+I2,IB+I1)-FSUM(8+I)
                       FORCE(IB+I1,IA+I2)=FORCE(IB+I1,IA+I2)-FSUM(8+I)
                       FORCE(IB+I2,IA+I1)=FORCE(IB+I2,IA+I1)-FSUM(8+I)
 2820               CONTINUE
                ENDIF
C              Do whichever processing is needed for this atom pair
C              during phase 1, and pick next pair.
                IF(DORESP) THEN
C                  Compute mixed cartesian-density matrix derivatives.
                    CALL PSDERO(AINTS,AOVER)
                    IF( IMIX.EQ.3 ) THEN
C                     Store computed components of mixed derivatives.
C                     If we are not interested in gradients on point charges,
C                     this will combine all contributions due to those.
                       DO 2918 J=1,3
                          DO 2900 I=1,NPAIRA
                             DUMP(ITXPDA+I) = DUMP(ITXPDA+I) + EAA(I,J)
 2900                     CONTINUE
                          ITXPDA = ITXPDA + NPAIRA
 2918                  CONTINUE
                       IF( .NOT.DOPTCH ) ITXPDA = ITXPDA - 3*NPAIRA
                    ELSE
C                     Combine computed components to form complete matrix
                       DO 2930 I=1,3
                          J = IDXPA+(3*NA+I)*NROS
                          CALL PSDAPD(DUMP(J),NFRSTA,NORBA,EAA(1,I))
                          IF(DOPTCH) THEN
                             J = IDXPA+(3*(NB+NATOM)+I)*NROS
                             CALL PSDSPD(DUMP(J),NFRSTA,NORBA,EAA(1,I))
                          ENDIF
                          IF(UHF) THEN
                             J = IDXPB+(3*NA+I)*NROS
                             CALL PSDAPD(DUMP(J),NFRSTA,NORBA,EAA(1,I))
                             IF(DOPTCH) THEN
                                J = IDXPB+(3*(NB+NATOM)+I)*NROS
                                CALL PSDSPD(DUMP(J),NFRSTA,NORBA,
     .                                      EAA(1,I))
                             ENDIF
                          ENDIF
 2930                  CONTINUE
                    ENDIF
                    IF(LPINTS) THEN
                        WRITE(NB6,10500)
                        CALL PSDPGM(NPAIRA,3,EAA,45)
                    ENDIF
                ENDIF
 4090       CONTINUE
          ENDIF
 9000 CONTINUE
C
      RETURN
10000 FORMAT(/' ATOM ', I5, ' HAS ', I5, ' ORBITALS, WHICH CANNOT BE',
     .        ' THE CASE FOR S, P OR D BASIS SET.'
     .       /' STOP NOW.' )
C
C    Debugging formats
C
10100 FORMAT(' PROCESSING PAIR ', I5, ', ', I5 )
10102 FORMAT(' PROCESSING PAIR ', I5, ' (ATOM), ', I5, '(CORE) ')
10110 FORMAT(' ZCOREA = ', F16.7, ' ZCOREB = ', F16.7/
     .       ' NA     = ', I5, ' NB     = ', I5, ' NTYPA  = ', I5,
     .       ' NTYPB  = ', I5, ' NFRSTA = ', I5/
     .       ' NORBA  = ', I5, ' NORBB  = ', I5, ' NFRSTB = ', I5,
     .       ' NPAIRA = ', I5, ' NPAIRB = ', I5/
     .       ' N1CNTA = ', I5 )
10200 FORMAT(' CONTRIBUTION FROM TOTAL DENSITY:' )
10210 FORMAT(' CONTRIBUTION FROM ALPHA DENSITY:' )
10220 FORMAT(' CONTRIBUTION FROM BETA  DENSITY:' )
10230 FORMAT(' CONTRIBUTION FROM CORE INTERACTION:' )
10235 FORMAT(' CONTRIBUTION FROM MO INTEGRALS:' )
10240 FORMAT(' COMPLETE PAIR SUM:' )
10300 FORMAT(' AA ALPHA DENSITY:' )
10310 FORMAT(' AA BETA  DENSITY:' )
10320 FORMAT(' AA DENSITY:' )
10330 FORMAT(' BB ALPHA DENSITY:' )
10340 FORMAT(' BB BETA  DENSITY:' )
10350 FORMAT(' BB DENSITY:' )
10360 FORMAT(' AB ALPHA DENSITY:' )
10370 FORMAT(' AB BETA  DENSITY:' )
10380 FORMAT(' AB DENSITY:' )
10390 FORMAT(' CONTRIBUTIONS FROM PAIR ', I4, ',', I4, ' TO OPEN M.O. ',
     .       ' AND DERIVATIVES:' )
10394 FORMAT('     XJ(',I2,'):')
10396 FORMAT('     XK(',I2,'):')
10398 FORMAT('     IJKL(',I2,'):')
10400 FORMAT(' TWO-CENTER INTEGRALS:' )
10405 FORMAT(' TWO-CENTER INTEGRALS DERIVATIVES, VAR ', I2, ':' )
10410 FORMAT(' RESONANCE INTEGRALS:' )
10500 FORMAT(' EAA:' )
10510 FORMAT(' EBB:' )
10520 FORMAT(' EABALP:' )
10530 FORMAT(' EABBET:' )
10600 FORMAT(' ONE-CENTER TWO-ELECTRON INTEGRALS FOR ATOM ', I5, ':' )
      END
C
      SUBROUTINE PSDABT(OVER)
C
C   Convert overlap integrals into resonance integrals by
C   multiplying them with beta paramaters.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      OVER   - Input: Overlap integrals
C               Output: Resonance integrals
C
C   Accessed common blocks:
C
C      PSPAIR - Pair papameters
C      PAROPT,
C      DPARM1 - Beta parameters
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C     This might be made to run much faster by moving
C     scaling by beta parameters to PSOVER.
C
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (HALFAU=0.183755972069092245497978684D-1)
C
      COMMON 
     ./PSPAIR/ ZCOREA, ZCOREB,
     .         NA, NB, NTYPA, NTYPB, NFRSTA, NORBA, NORBB, 
     .         NFRSTB, NPAIRA, NPAIRB, N1CNTA
      COMMON
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ),BETAP(LMZ),
     .         ALP(LMZ)
     ./DPARM1/ UDD(LMZ),ZD(LMZ),BETAD(LMZ)
C
      DIMENSION OVER(9,9)
C S-S
      SCALE = HALFAU * (BETAS(NTYPA) + BETAS(NTYPB))
      OVER(1,1) = OVER(1,1) * SCALE
      IF( NORBA.GT.1 ) THEN
C P-S
          SCALE = HALFAU * (BETAP(NTYPA) + BETAS(NTYPB))
          OVER(2,1) = OVER(2,1) * SCALE
          OVER(3,1) = OVER(3,1) * SCALE
          OVER(4,1) = OVER(4,1) * SCALE
          IF( NORBA.GT.4 ) THEN
C D-S
              SCALE = HALFAU * (BETAD(NTYPA) + BETAS(NTYPB))
              OVER(5,1) = OVER(5,1) * SCALE
              OVER(6,1) = OVER(6,1) * SCALE
              OVER(7,1) = OVER(7,1) * SCALE
              OVER(8,1) = OVER(8,1) * SCALE
              OVER(9,1) = OVER(9,1) * SCALE
          ENDIF
      ENDIF
      IF( NORBB.LE.1 ) RETURN
C S-P
      SCALE = HALFAU * (BETAS(NTYPA) + BETAP(NTYPB))
      OVER(1,2) = OVER(1,2) * SCALE
      OVER(1,3) = OVER(1,3) * SCALE
      OVER(1,4) = OVER(1,4) * SCALE
      IF( NORBA.GT.1 ) THEN
C P-P
          SCALE = HALFAU * (BETAP(NTYPA) + BETAP(NTYPB))
          DO 100 I=2,4
              OVER(2,I) = OVER(2,I) * SCALE
              OVER(3,I) = OVER(3,I) * SCALE
              OVER(4,I) = OVER(4,I) * SCALE
  100     CONTINUE
          IF( NORBA.GT.4 ) THEN
C D-P
              SCALE = HALFAU * (BETAD(NTYPA) + BETAP(NTYPB))
              DO 200 I=2,4
                  OVER(5,I) = OVER(5,I) * SCALE
                  OVER(6,I) = OVER(6,I) * SCALE
                  OVER(7,I) = OVER(7,I) * SCALE
                  OVER(8,I) = OVER(8,I) * SCALE
                  OVER(9,I) = OVER(9,I) * SCALE
  200         CONTINUE
          ENDIF
      ENDIF
      IF( NORBB.LE.4 ) RETURN
C S-D
      SCALE = HALFAU * (BETAS(NTYPA) + BETAD(NTYPB))
      OVER(1,5) = OVER(1,5) * SCALE
      OVER(1,6) = OVER(1,6) * SCALE
      OVER(1,7) = OVER(1,7) * SCALE
      OVER(1,8) = OVER(1,8) * SCALE
      OVER(1,9) = OVER(1,9) * SCALE
      IF( NORBA.GT.1 ) THEN
C P-D
          SCALE = HALFAU * (BETAP(NTYPA) + BETAD(NTYPB))
          DO 300 I=5,9
              OVER(2,I) = OVER(2,I) * SCALE
              OVER(3,I) = OVER(3,I) * SCALE
              OVER(4,I) = OVER(4,I) * SCALE
  300     CONTINUE
          IF( NORBA.GT.4 ) THEN
C D-D
              SCALE = HALFAU * (BETAD(NTYPA) + BETAD(NTYPB))
              DO 400 I=5,9
                  OVER(5,I) = OVER(5,I) * SCALE
                  OVER(6,I) = OVER(6,I) * SCALE
                  OVER(7,I) = OVER(7,I) * SCALE
                  OVER(8,I) = OVER(8,I) * SCALE
                  OVER(9,I) = OVER(9,I) * SCALE
  400         CONTINUE
          ENDIF
      ENDIF
      RETURN
      END
C
      SUBROUTINE PSDPAA(ROALP,ROBET,NFIRST,ICOUNT,P)
C
C   Fetch diagonal block of density matrix to the local 
C   storage array, multiplying all off-diagonal elements
C   of density matrix by 2. Also store diagonal entries of
C   total density matrix into separate array. This should 
C   improve performance for large computations by increasing 
C   cache hit ratio and shouldn't hurt anybody otherwise, 
C   according to test benchmarks.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ROALP  - Packed alpha density matrix
C      ROBET  - Packed beta density matrix
C               Net accessed if UHF flag not set
C      NFIRST - First orbital of the square orbital
C               block to fetch density for
C      ICOUNT - Number of orbitals
C      P      - Storage place
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C
C   Modified common blocks:
C
C      PSDENS - Modified indirectly through passed 
C               parameter P
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TWO=2.0D0)
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
C
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      SAVE /PSDGBL/
C
      DIMENSION ROALP(*), ROBET(*), P(3*45+1)
C
      IFIRST = ( NFIRST * (NFIRST+1) ) / 2
      IOUT   = 0
      IIN    = IFIRST-1
C
      IF( UHF ) THEN
          DO 120 I = 1, ICOUNT
              DO 110 J=1,I-1
                  P( 0+IOUT+J) = TWO * ROALP(IIN+J)
                  P(45+IOUT+J) = TWO * ROBET(IIN+J)
                  P(91+IOUT+J) = P(0+IOUT+J) + P(45+IOUT+J)
  110         CONTINUE
              P( 0+IOUT+I) = ROALP(IIN+I)
              P(45+IOUT+I) = ROBET(IIN+I)
              P(91+IOUT+I) = P(0+IOUT+I) + P(45+IOUT+I)
              IOUT = IOUT + I
              IIN  = IIN  + I + NFIRST-1
  120     CONTINUE
      ELSE
          DO 220 I = 1, ICOUNT
              DO 210 J=1,I-1
                  P( 0+IOUT+J) = TWO * ROALP(IIN+J)
                  P(91+IOUT+J) = TWO * P(0+IOUT+J)
  210         CONTINUE
              P( 0+IOUT+I) = ROALP(IIN+I)
              P(91+IOUT+I) = TWO * P(0+IOUT+I)
              IOUT = IOUT + I
              IIN  = IIN  + I + NFIRST-1
  220     CONTINUE
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSDPAB(ROALP,ROBET,NFRSTA,NFRSTB,ICNTA,ICNTB,P)
C
C   Fetch off-diagonal block of density matrix to the local 
C   storage array.  This should improve performance for 
C   large computations by increasing cache hit ratio and 
C   shouldn't hurt anybody otherwise, according to test 
C   benchmarks.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ROALP  - Packed alpha density matrix
C      ROBET  - Packed beta density matrix
C               Net accessed if UHF flag not set
C      NFRSTA - First orbital on atom A
C      NFRSTB - First orbital on atom B
C               This should be less than NFRSTA
C               (Strictly speaking, NFRSTB+ICNTB
C               should be less than or equal to NFRSTA)
C      ICNTA  - Number of orbitals on atom A
C      ICNTB  - Number of orbitals on atom B
C      P      - Storage place
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C
C   Modified common blocks:
C
C      PSDENS - Modified indirectly through passed 
C               parameter P
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TWO=2.0D0)
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
C
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSPRT / NB6
      SAVE /PSDGBL/, /PSPRT /
C
      DIMENSION ROALP(*), ROBET(*), P(9,9,3)
C
      IF( NFRSTA.LE.NFRSTB ) THEN
          WRITE(NB6,10000)
          STOP 'PSDPAB'
      ENDIF
      IFIRST = ( NFRSTA * (NFRSTA+1) ) / 2 - NFRSTA + NFRSTB
      IIN    = IFIRST-1
C
      IF( UHF ) THEN
          DO 120 I = 1, ICNTA
              DO 110 J=1,ICNTB
                  P(I,J,1) = ROALP(IIN+J)
                  P(I,J,2) = ROBET(IIN+J)
                  P(I,J,3) = P(I,J,1) + P(I,J,2)
  110         CONTINUE
              IIN  = IIN  + I + NFRSTA-1
  120     CONTINUE
      ELSE
          DO 220 I = 1, ICNTA
              DO 210 J=1,ICNTB
                  P(I,J,1) = ROALP(IIN+J)
                  P(I,J,3) = TWO * P(I,J,1)
  210         CONTINUE
              IIN  = IIN  + I + NFRSTA-1
  220     CONTINUE
      ENDIF
C
      RETURN
10000 FORMAT(' PROGRAMMING ERROR IN PSDST1' )
      END
C
      SUBROUTINE PSDFPS(AINTS,AOVER,F)
C
C   Compute component of derivative of pair-energy which 
C   depend on total density but not on individual alpha and 
C   beta densities.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      AINTS  - Two-center, two-electron integrals and core
C               interaction integrals
C      AOVER  - Resonance integrals
C      F      - Output sum value
C
C   Accessed common blocks:
C
C      PSPAIR - Pair papameters
C      PSDENS - Density matrices
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
      PARAMETER (TWO=2.0D0)
C
      COMMON 
     ./PSPAIR/ ZCOREA, ZCOREB,
     .         NA, NB, NTYPA, NTYPB, NFRSTA, NORBA, NORBB, 
     .         NFRSTB, NPAIRA, NPAIRB, N1CNTA
     ./PSDENS/ PAAALP(45), PAABET(45), ZAA, PAA(45),
     .         PBBALP(45), PBBBET(45), ZBB, PBB(45),
     .         PABALP(9,9), PABBET(9,9), PAB(9,9)
C
      DIMENSION AINTS(46,46), AOVER(9,9)
C
      F = 0D0
C
      DO 120 IB = 1, NPAIRB
          DO 110 IA = 1, NPAIRA
              F = F + AINTS(1+IA,1+IB) * PAA(IA) * PBB(IB)
  110     CONTINUE
  120 CONTINUE
C
      DO 210 IA = 1, NPAIRA
          F = F - ZCOREB * PAA(IA) * AINTS(1+IA,1)
  210 CONTINUE
C
      DO 310 IB = 1, NPAIRB
          F = F - ZCOREA * PBB(IB) * AINTS(1,1+IB)
  310 CONTINUE
C
      TF = 0D0
      DO 420 IB = 1, NORBB
C         IFORT 12.1 WILL BREAK THE CODE HERE WHEN VECTORIZING.
CDEC$     NOVECTOR
          DO 410 IA = 1, NORBA
              TF = TF + AOVER(IA,IB) * PAB(IA,IB)
  410     CONTINUE
  420 CONTINUE
      F = F + TWO * TF
C
      RETURN
      END
C
      SUBROUTINE PSDFPX(AINTS,P,F)
C
C   Compute component of derivative of pair-energy which 
C   depend on individual alpha and beta densities
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      AINTS  - Two-center, two-electron integrals and core
C               interaction integrals (core interaction
C               integrals are not used)
C      P      - Block of density matrix shared by A and B
C      F      - Output sum value
C
C   Accessed common blocks:
C
C      PSPAIR - Pair papameters
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
      PARAMETER (TWO=2.0D0)
      PARAMETER (HALF=0.5D0)
C
      COMMON 
     ./PSPAIR/ ZCOREA, ZCOREB,
     .         NA, NB, NTYPA, NTYPB, NFRSTA, NORBA, NORBB, 
     .         NFRSTB, NPAIRA, NPAIRB, N1CNTA
C
      DIMENSION AINTS(46,46), P(9,9)
      INTEGER MU, NU, LA, SI
C
      F = 0D0
C
      IINTB = 2
      DO 420 LA = 1, NORBB
          DO 210 SI = 1, LA - 1
C
              IINTA = 2
              DO 120 MU = 1, NORBA
                  DO 110 NU = 1, MU - 1
                      F = F + AINTS(IINTA,IINTB) * (P(MU,LA)* 
     .                        P(NU,SI) + P(MU,SI)*P(NU,LA))
                      IINTA = IINTA + 1
  110             CONTINUE
                  F = F + AINTS(IINTA,IINTB)*P(MU,LA)*P(MU,SI)
                  IINTA = IINTA + 1
  120         CONTINUE
C
              IINTB = IINTB + 1
  210     CONTINUE
C
          IINTA = 2
          DO 320 MU = 1, NORBA
              DO 310 NU = 1, MU - 1
                  F = F + AINTS(IINTA,IINTB)*P(MU,LA)*P(NU,LA)
                  IINTA = IINTA + 1
  310         CONTINUE
              F = F + HALF * AINTS(IINTA,IINTB)*(P(MU,LA)**2)
              IINTA = IINTA + 1
  320     CONTINUE
C
          IINTB = IINTB + 1
  420 CONTINUE
      F = -TWO * F
C
      RETURN
      END
C
      SUBROUTINE PSDFCR(AINTS,F)
C
C   Compute core repulsion contribution and it's derivatives
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      AINTS  - Two-center integrals (only (1,1,X) entry
C               is used)
C      F      - Core repulsion energy w/ derivatives
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C      PSDIST - Interatomic distance, computed by PSINTS
C      PSPAIR - Pair papameters
C      PAROPT - Core repulsion parameters
C      PSDROT - Derivatives projection coefficients.
C      AMPGAU - Additional parameters for core repulsion 
C               in AM1 and PM3
C      ABOND  - Additional parameters for core repulsion 
C               in MNDO/d
C      INOPT2 - Current values of input flags:
C               IN2(2)=IOP defines method that is used.
C               IN2(122)=MMPOT defines interaction potential
C               for point charges
C      QMM3   - Core repulsion parameters for point charges
C
C   Modified common blocks:
C
C   Local storage:
C
C      5*MAXGAU (=40) DOUBLE PRECISION cells are used for
C      differentiating gaussians in AM1 and PM3 core repulsion
C      formulas.
C
C   Module logic:
C
C   Bugs:
C
C      Naive method with double rotation is employed for
C      accounting for functional contribution to derivative.
C
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (BOHR=0.529167D0)
      PARAMETER (AUINEV=0.03675119441381844909959573686D0)
      PARAMETER (ZERO =0.0D0)
      PARAMETER (ONE  =1.0D0)
      PARAMETER (TWO  =2.0D0)
      PARAMETER (IELH = 1)
      PARAMETER (IELB = 5)
      PARAMETER (IELC = 6)
      PARAMETER (IELN = 7)
      PARAMETER (IELO = 8)
      PARAMETER (IELF = 9)
      PARAMETER (IELCL=17)
      PARAMETER (IELBR=35)
      PARAMETER (IELI =53)
      PARAMETER (MAXGAU=8)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSDIST/ DIST
     ./PSPAIR/ ZCOREA, ZCOREB,
     .         NA, NB, NTYPA, NTYPB, NFRSTA, NORBA, NORBB, 
     .         NFRSTB, NPAIRA, NPAIRB, N1CNTA
     ./PSDROT/ D1(3), D2A(6), D2B(6)
     ./PSPRT / NB6
      COMMON
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ),BETAP(LMZ),
     .         ALP(LMZ)
     ./AMPGAU/ FN1(LMZ,4),FN2(LMZ,4),FN3(LMZ,4),GSCAL(LMZ),IMP(LMZ)
     ./ABOND / ALPB(LMZ,LMZ),MALPB(LMZ)
     ./INOPT2/ IN2(300)
     ./QMMM3 / DELTAM(LMZ),OMEGAM(LMZ)
      SAVE /PSDROT/, /PSDIST/, /PSDGBL/, /PSPRT /
C
      DIMENSION AINTS(46,46,10), F(10)
C
      DIMENSION FUNC(10), AINT(10)
      DIMENSION A(MAXGAU), B(MAXGAU), C(MAXGAU), G(MAXGAU), T(MAXGAU)
C
      IF( METHOD.LT.0 .OR. METHOD.GT.2 ) THEN
          WRITE(NB6,10000) METHOD
          STOP 'PSDFCR'
      ENDIF
C     Input options.
      IOP   = IN2(2)
      MMPOT = IN2(122)
C     Common terms for all methods
      XDIST = DIST * BOHR
      ZZ    = ZCOREA * ZCOREB
C     ITYPA and ITYPB are equivalent to NTYPA and NTYPB, but can be modified
C     in certain circumstances.
      ITYPA = NTYPA
      ITYPB = NTYPB
C
C     Assign standard atomic parameters.
      IF(NTYPA.GT.0) ALPI = ALP(NTYPA)
      IF(NTYPB.GT.0) ALPJ = ALP(NTYPB)
C     Bond parameters for MNDO/d
      IF(IOP.EQ.-10 .AND. NTYPA.GT.0 .AND. NTYPB.GT.0) THEN
         MALPNI = MALPB(NTYPA)
         MALPNJ = MALPB(NTYPB)
         IF(MALPNI.GT.0 .AND. NTYPB.LE.MALPNI) THEN
            IF(ALPB(NTYPB,NTYPA).GT.ZERO) ALPI = ALPB(NTYPB,NTYPA)
         ENDIF
         IF(MALPNJ.GT.0 .AND. NTYPA.LE.MALPNJ) THEN
            IF(ALPB(NTYPA,NTYPB).GT.ZERO) ALPJ = ALPB(NTYPA,NTYPB)
         ENDIF
      ENDIF
C     Distinguish between atom-atom and atom-charge case.
      IF( NTYPB.NE.0 ) THEN
C         Second atom is not a pure core, use normal parameters.
          EXPA  = EXP( - ALPI * XDIST )
          EXPB  = EXP( - ALPJ * XDIST )
      ELSE
C         Second atom is a pure core, get special parameters. Shifting
C         exponent's argument by a constant is equivalent to multiplying
C         this term by a constant, so that derivatives below are still
C         valid.
          ALPI  = OMEGAM(NTYPA)
          EXPA  = EXP( - ALPI * ( XDIST - DELTAM(NTYPA) ) )
          ALPJ  = ZERO
          EXPB  = ZERO
C         In two special cases, point charge has an interaction potential
C         of its own. MMPOT.EQ.7 uses a very sharp repulsive core, while
C         MMPOT.EQ.8 uses hydrogen "core" parameter, whichever it happens
C         to be. In the later case, we'll have to add hydrogen's Gaussians
C         for AM1 and PM3, too.
          IF( MMPOT.EQ.7 ) THEN
              ALPI  = 5.D0
              EXPB  = EXP( - ALPI * XDIST )
          ENDIF
          IF( MMPOT.EQ.8 ) THEN
              ALPI  = ALP(IELH)
              EXPB  = EXP( - ALPI * XDIST )
              ITYPB = IELH
          ENDIF
      ENDIF
C     End of the special treatment
C
      IF(NTYPA.EQ.IELH.AND.(NTYPB.EQ.IELN.OR.NTYPB.EQ.IELO)) THEN
          VAL = 1 + EXPA + XDIST*EXPB
      ELSE IF((NTYPA.EQ.IELN.OR.NTYPA.EQ.IELO).AND.NTYPB.EQ.IELH) THEN
          VAL = 1 + XDIST*EXPA + EXPB
      ELSE
          VAL = 1 + EXPA + EXPB
      ENDIF
      VAL = VAL * ZZ
      FUNC(1) = VAL
      AINT(1) = AINTS(1,1,1)
      IF( MODE.GE.1 ) THEN
          IF(NTYPA.EQ.IELH.AND.(NTYPB.EQ.IELN.OR.NTYPB.EQ.IELO)) THEN
              DRV = -ALPI*EXPA - (ALPJ*XDIST-ONE)*EXPB
          ELSE IF((NTYPA.EQ.IELN.OR.NTYPA.EQ.IELO).AND.
     .                                            NTYPB.EQ.IELH) THEN
              DRV = -(ALPI*XDIST-ONE)*EXPA - ALPJ*EXPB
          ELSE
              DRV = -ALPI*EXPA - ALPJ*EXPB
          ENDIF
          DRV = DRV * ZZ * BOHR
          DO 110 I=1,3
              FUNC(1+I) = D1(I) * DRV
              AINT(1+I) = AINTS(1,1,1+I)
  110     CONTINUE
      ENDIF
      IF( MODE.GE.2 ) THEN
          IF(NTYPA.EQ.IELH.AND.(NTYPB.EQ.IELN.OR.NTYPB.EQ.IELO)) THEN
              FRC = (ALPI**2)*EXPA + ((ALPJ**2)*XDIST-TWO*ALPJ)*EXPB
          ELSE IF((NTYPA.EQ.IELN.OR.NTYPA.EQ.IELO).AND.
     .                                            NTYPB.EQ.IELH) THEN
              FRC = ((ALPI**2)*XDIST-TWO*ALPI)*EXPA + (ALPJ**2)*EXPB
          ELSE
              FRC = (ALPI**2)*EXPA + (ALPJ**2)*EXPB
          ENDIF
          FRC = FRC * ZZ * (BOHR**2)
          DO 120 I=1,6
              FUNC(4+I) = D2A(I) * FRC + D2B(I) * DRV
              AINT(4+I) = AINTS(1,1,4+I)
  120     CONTINUE
      ENDIF
C
C    Multiply (c,c) integrals with derivatives by scaling function
C
      CALL PSXML(MODE,AINT,FUNC,F)
      IF( METHOD.EQ.0 ) RETURN
      IF( NTYPB.EQ.0 .AND. MMPOT.NE.7 .AND. MMPOT.NE.8 ) RETURN
C
C     AM1 and PM3 need additional terms here.
C
      IGCNT = 0
      IF(METHOD.EQ.1.AND.(ITYPA.EQ.IELB.OR.ITYPB.EQ.IELB)) THEN
C
C         Pair parameters for boron in AM1
C         This approach _is_ ugly, I would prefer to add dummy
C         indices to the element parameters table and just substitute
C         these indices for boron, but I have no control over 
C         format of parameters, unfortunately.
C
          ITYPB = ITYPA+ITYPB-IELB
          ITYPA = IELB
          IF( ITYPB.EQ.IELH ) THEN
              A(1)  = 0.412253D0 
              A(2)  =-0.149917D0
              B(1)  = 10.0D0
              B(2)  =  6.0D0
              C(1)  = 0.832586D0
              C(2)  = 1.186220D0
              IGCNT = 2
          ELSE IF( ITYPB.EQ.IELC ) THEN
              A(1)  = 0.261751D0
              A(2)  = 0.050275D0
              B(1)  = 8.0D0
              B(2)  = 5.0D0
              C(1)  = 1.063995D0
              C(2)  = 1.936492D0
              IGCNT = 2
          ELSE IF( ITYPB.EQ.IELF  .OR. ITYPB.EQ.IELCL .OR. 
     .             ITYPB.EQ.IELBR .OR. ITYPB.EQ.IELI ) THEN
              A(1)  = 0.359244D0
              A(2)  = 0.074729D0 
              B(1)  = 9.0D0
              B(2)  = 9.0D0
              C(1)  = 0.819351D0
              C(2)  = 1.574414D0
              IGCNT = 2
          ELSE
              A(1)  = 0.182613D0
              A(2)  = 0.118587D0
              A(3)  =-0.073280D0
              B(1)  = 6.0D0
              B(2)  = 6.0D0
              B(3)  = 5.0D0
              C(1)  = 0.727592D0
              C(2)  = 1.466639D0
              C(3)  = 1.570975D0
              IGCNT = 3
          ENDIF
      ELSE
C
C         General case parameters for first atom
C
          DO 500 I=1,IMP(ITYPA)
              A(I) = FN1(ITYPA,I)
              B(I) = FN2(ITYPA,I)
              C(I) = FN3(ITYPA,I)
  500     CONTINUE
          IGCNT = IMP(ITYPA)
      ENDIF
      IF( ITYPB.EQ.ITYPA ) THEN
C
C         Just double A coefficients - atom B have the same
C         gaussians. It would be an error to consult FN1-FN3
C         now, because second element could as well be boron...
C
          DO 600 I=1,IGCNT
              A(I) = TWO * A(I)
  600     CONTINUE
      ELSE
C
C         Element B is different from A, fetch additional
C         gaussians.
C
          IF( ITYPB.NE.0 ) THEN
              DO 700 I=1,IMP(ITYPB)
                  A(IGCNT+I) = FN1(ITYPB,I)
                  B(IGCNT+I) = FN2(ITYPB,I)
                  C(IGCNT+I) = FN3(ITYPB,I)
  700         CONTINUE
              IGCNT = IGCNT + IMP(ITYPB)
          ENDIF
      ENDIF
      IF( IGCNT.GT.MAXGAU ) THEN
          WRITE(NB6,10010)
          STOP 'PSDFCR'
      ENDIF
C
C     Compute additional term
C
      RXDIST = ONE / XDIST
      GG     = ONE
      IF(ITYPA.GT.0)  GG = GG * GSCAL(ITYPA)
      IF(ITYPB.GT.0)  GG = GG * GSCAL(ITYPB)
      DO 800 I=1,IGCNT
          A(I) = A(I) * AUINEV
          G(I) = A(I) * EXP( -B(I) * (XDIST-C(I))**2 ) * RXDIST * ZZ *
     .           GG
          F(1) = F(1) + G(I)
  800 CONTINUE
      IF( MODE.GE.1 ) THEN
          DRV    = ZERO
          XDIST2 = XDIST**2
          DO 850 I=1,IGCNT
              T(I) = -ONE + TWO*B(I)*C(I)*XDIST - TWO*B(I)*XDIST2
              DRV = DRV + G(I)*T(I)*RXDIST
  850     CONTINUE
          DRV  = DRV * BOHR
          F(2) = F(2) + DRV*D1(1)
          F(3) = F(3) + DRV*D1(2)
          F(4) = F(4) + DRV*D1(3)
      ENDIF
      IF( MODE.GE.2 ) THEN
          FRC    = ZERO
          RXDST2 = RXDIST ** 2
          DO 900 I=1,IGCNT
              FRC = FRC + RXDST2*TWO*G(I)*( -T(I) - B(I)*XDIST2 + 
     .                  TWO * (B(I)*XDIST*(C(I)-XDIST))**2)
  900     CONTINUE
          FRC   = FRC * BOHR**2
          F( 5) = F( 5) + FRC*D2A(1) + DRV*D2B(1)
          F( 6) = F( 6) + FRC*D2A(2) + DRV*D2B(2)
          F( 7) = F( 7) + FRC*D2A(3) + DRV*D2B(3)
          F( 8) = F( 8) + FRC*D2A(4) + DRV*D2B(4)
          F( 9) = F( 9) + FRC*D2A(5) + DRV*D2B(5)
          F(10) = F(10) + FRC*D2A(6) + DRV*D2B(6)
      ENDIF

      RETURN
10000 FORMAT(' CORE REPULSION FUNCTION NOT DEFINED FOR METHOD ', I4
     .      /' PROGRAM WILL STOP.' )
10010 FORMAT(' TOO MANY GAUSSIANS IN CORE REPULSION FUNCTION',
     .      /' PROGRAM WILL STOP.' )
      END
C
      SUBROUTINE PSDERO(AINTS,AOVER)
C
C   Compute mixed derivatives of energy with respect to density
C   matrix elements and atomic coordinates
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      AINTS  - Two-center electron integrals
C      AOVER  - Resonance integrals
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C      PSPAIR - Pair papameters
C      PSDENS - Density matrices
C
C   Modified common blocks:
C
C      PSEXRO - Mixed derivatives with respect to coordinate/density
C
C   Local storage:
C
C   Module logic:
C
C      Unfortunately, I see no way of doing spin-dependent portion
C      with a BLAS call...
C
C   Speedups possible:
C
C      Inlining calls to DGEMV and PSDERX shoudl probably speed up
C      this module on most machines.
C
C   Bugs:
C
C      Strictly speaking, computed quantities are defined only if
C      P is simultaneously subjected to idempotency condition. 
C      Equivalent well-defined quantities (which are computed later
C      as the right-hand sides of CPHF equations) are much more 
C      expensive to recompute with N^2 limit on memory use.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ONE =1D0)
      PARAMETER (ZERO=0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSPAIR/ ZCOREA, ZCOREB,
     .         NA, NB, NTYPA, NTYPB, NFRSTA, NORBA, NORBB, 
     .         NFRSTB, NPAIRA, NPAIRB, N1CNTA
     ./PSDENS/ PAAALP(45), PAABET(45), ZAA, PAA(45),
     .         PBBALP(45), PBBBET(45), ZBB, PBB(45),
     .         PABALP(9,9), PABBET(9,9), PAB(9,9)
     ./PSEXRO/ EAA(45,3), EBB(45,3), EABALP(9,9,3), EABBET(9,9,3)
      SAVE /PSDGBL/
C
C*$*INLINE ROUTINE (PSZRV,DGEMV,PSDERX)
CVD$R EXPAND(PSZRV,DGEMV,PSDERX)
CVD$R SEARCH(psmutl.f,dblas2.f)
C
      DIMENSION AINTS(46,46,10), AOVER(9,9,10)
C
      DIMENSION PAAX(46), PBBX(46)
      EQUIVALENCE (ZAA,PAAX(1)),(ZBB,PBBX(1))
C
      ZAA = -ZCOREA
      ZBB = -ZCOREB
      CALL PSZRMB(NPAIRA,3,EAA,45)
      CALL PSZRMB(NPAIRB,3,EBB,45)
      DO 100 I=1,3
          CALL DGEMV('N',NPAIRA,NPAIRB+1,ONE,AINTS(2,1,1+I),46,
     .               PBBX,1,ZERO,EAA(1,I),1)
          CALL DGEMV('T',NPAIRA+1,NPAIRB,ONE,AINTS(1,2,1+I),46,
     .               PAAX,1,ZERO,EBB(1,I),1)
  100 CONTINUE
C
      CALL PSDERX(AINTS,AOVER,PABALP,EABALP)
      IF( UHF ) THEN
          CALL PSDERX(AINTS,AOVER,PABBET,EABBET)
      ENDIF
      RETURN
      END
C
      SUBROUTINE PSDERX(AINTS,AOVER,PAB,EAB)
C
C   Compute mixed derivatives of energy with respect to the density
C   matrix elements of the single spin and atomic coordinates
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      AINTS  - Two-center electron integrals
C      AOVER  - Resonance integrals
C      PAB    - Partial density matrix
C      EAB    - Output derivatives
C
C   Accessed common blocks:
C
C      PSPAIR - Pair papameters
C
C   Modified common blocks:
C
C      PSEXRO - Modified through passed parameter EAB
C
C   Local storage:
C
C   Module logic:
C
C   Speedups possible:
C
C      Unrolling specific versions for {NORBA,NORBB} = {1,4,9}
C      might improve performance on vector and heavily pipelined
C      machines.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON 
     ./PSPAIR/ ZCOREA, ZCOREB,
     .         NA, NB, NTYPA, NTYPB, NFRSTA, NORBA, NORBB, 
     .         NFRSTB, NPAIRA, NPAIRB, N1CNTA
C
      DIMENSION AINTS(46,46,10), AOVER(9,9,10)
      DIMENSION PAB(9,9), EAB(9,9,3)
      INTEGER MU, NU, LA, SI
C    Early return, if possible
      IF( NORBA.EQ.0 .OR. NORBA.EQ.0 ) RETURN
      DO 110 MU=1,NORBA
          DO 100 LA=1,NORBB
              EAB(MU,LA,1) = AOVER(MU,LA,2)
              EAB(MU,LA,2) = AOVER(MU,LA,3)
              EAB(MU,LA,3) = AOVER(MU,LA,4)
  100     CONTINUE
  110 CONTINUE
      IINTB = 2
      DO 210 LA = 1, NORBB
          DO 160 SI = 1, LA-1
              IINTA = 2
              DO 150 MU = 1, NORBA
                  DO 130 NU = 1, MU-1
                      DO 120 I=1,3
                          AI    = AINTS(IINTA,IINTB,I+1)
                          EAB(MU,LA,I) = EAB(MU,LA,I) - AI*PAB(NU,SI)
                          EAB(NU,LA,I) = EAB(NU,LA,I) - AI*PAB(MU,SI)
                          EAB(MU,SI,I) = EAB(MU,SI,I) - AI*PAB(NU,LA)
                          EAB(NU,SI,I) = EAB(NU,SI,I) - AI*PAB(MU,LA)
  120                 CONTINUE
                      IINTA = IINTA + 1
  130             CONTINUE
                  DO 140 I=1,3
                      AI           = AINTS(IINTA,IINTB,I+1)
                      EAB(MU,LA,I) = EAB(MU,LA,I) - AI*PAB(MU,SI)
                      EAB(MU,SI,I) = EAB(MU,SI,I) - AI*PAB(MU,LA)
  140             CONTINUE
                  IINTA = IINTA + 1
  150         CONTINUE
          IINTB = IINTB + 1
  160     CONTINUE
          IINTA = 2
          DO 200 MU = 1, NORBA
              DO 180 NU = 1, MU-1
                  DO 170 I=1,3
                      AI           = AINTS(IINTA,IINTB,I+1)
                      EAB(MU,LA,I) = EAB(MU,LA,I)-AI*PAB(NU,LA)
                      EAB(NU,LA,I) = EAB(NU,LA,I)-AI*PAB(MU,LA)
  170             CONTINUE
                  IINTA = IINTA + 1
  180         CONTINUE
              DO 190 I=1,3
                  AI           = AINTS(IINTA,IINTB,I+1)
                  EAB(MU,LA,I) = EAB(MU,LA,I)-AI*PAB(MU,LA)
  190         CONTINUE
              IINTA = IINTA + 1
  200     CONTINUE
          IINTB = IINTB + 1
  210 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSDAPD(A,NFIRST,ICOUNT,X)
C
C   Add pair contribution to diagonal block of mixed 
C   coordinate-density derivative held in packed form
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      A      - Matrix contaning mixed derivative
C      NFIRST - First orbital of the square orbital
C               block to add to whole matrix
C      ICOUNT - Number of orbitals
C      P      - Pair contribution
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      This is essentially reverse of PSDPAA, see it for 
C      comments.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(*), X(45)
C
      IFIRST = ( NFIRST * (NFIRST+1) ) / 2
      IIN    = 0
      IOUT   = IFIRST-1
C
      DO 220 I = 1, ICOUNT
          DO 210 J=1,I
              A(IOUT+J) = A(IOUT+J) + X(IIN+J)
  210     CONTINUE
          IIN  = IIN  + I
          IOUT = IOUT + I + NFIRST-1
  220 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSDSPD(A,NFIRST,ICOUNT,X)
C
C   Subtract pair contribution from the diagonal block of mixed 
C   coordinate-density derivative held in packed form
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      A      - Matrix contaning mixed derivative
C      NFIRST - First orbital of the square orbital
C               block to subtract from whole matrix
C      ICOUNT - Number of orbitals
C      P      - Pair contribution
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      This is essentially reverse of PSDPAA, see it for 
C      comments.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(*), X(45)
C
      IFIRST = ( NFIRST * (NFIRST+1) ) / 2
      IIN    = 0
      IOUT   = IFIRST-1
C
      DO 220 I = 1, ICOUNT
          DO 210 J=1,I
              A(IOUT+J) = A(IOUT+J) - X(IIN+J)
  210     CONTINUE
          IIN  = IIN  + I
          IOUT = IOUT + I + NFIRST-1
  220 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSDAPA(A,NFRSTA,NFRSTB,ICNTA,ICNTB,X)
C
C   Add pair contribution to the off-diagonal block of the
C   mixed coordinate-density derivative held in packed storage
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      A      - Packed mixed derivatives matrix
C      NFRSTA - First orbital on atom A
C      NFRSTB - First orbital on atom B
C               This should be less than NFRSTA
C               (Strictly speaking, NFRSTB+ICNTB
C               should be less than or equal to NFRSTA)
C      ICNTA  - Number of orbitals on atom A
C      ICNTB  - Number of orbitals on atom B
C      X      - What to add
C
C   Accessed common blocks:
C
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
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION A(*), X(9,9)
C
      IF( NFRSTA.LE.NFRSTB ) THEN
          WRITE(NB6,10000)
          STOP 'PSDAPA'
      ENDIF
      IFIRST = ( NFRSTA * (NFRSTA+1) ) / 2 - NFRSTA + NFRSTB
      IOUT   = IFIRST-1
C
      DO 220 I=1,ICNTA
          DO 210 J=1,ICNTB
              A(IOUT+J) = A(IOUT+J) + X(I,J)
  210     CONTINUE
          IOUT = IOUT + I + NFRSTA-1
  220 CONTINUE
C
      RETURN
10000 FORMAT(' PROGRAMMING ERROR IN PSDAPA' )
      END
C 
      SUBROUTINE PSDSPA(A,NFRSTA,NFRSTB,ICNTA,ICNTB,X)
C
C   Subtract pair contribution from the off-diagonal block of the
C   mixed coordinate-density derivative held in packed storage
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      A      - Packed mixed derivatives matrix
C      NFRSTA - First orbital on atom A
C      NFRSTB - First orbital on atom B
C               This should be less than NFRSTA
C               (Strictly speaking, NFRSTB+ICNTB
C               should be less than or equal to NFRSTA)
C      ICNTA  - Number of orbitals on atom A
C      ICNTB  - Number of orbitals on atom B
C      X      - What to subtract
C
C   Accessed common blocks:
C
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
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION A(*), X(9,9)
C
      IF( NFRSTA.LE.NFRSTB ) THEN
          WRITE(NB6,10000)
          STOP 'PSDAPA'
      ENDIF
      IFIRST = ( NFRSTA * (NFRSTA+1) ) / 2 - NFRSTA + NFRSTB
      IOUT   = IFIRST-1
C
      DO 220 I=1,ICNTA
          DO 210 J=1,ICNTB
              A(IOUT+J) = A(IOUT+J) - X(I,J)
  210     CONTINUE
          IOUT = IOUT + I + NFRSTA-1
  220 CONTINUE
C
      RETURN
10000 FORMAT(' PROGRAMMING ERROR IN PSDAPA' )
      END
C 
