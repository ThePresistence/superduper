C    ******************************************************************
C
C     Dynamic memory allocation for analytical derivative code.
C
C     ******************************************************************
      SUBROUTINE PSDADM(ISTAT,LDC,ITOP,IDSK,TRIAL)
C
C   Allocate dynamic memory blocks according to values set
C   in /PSDOPT/, /PSDGBL/ and /PSDGB2/. Set ISTAT to non-zero
C   value if not enough memory is available.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ISTAT  - Output: completion status
C           0 = Memory requirements could be met
C           1 = Memory limit reached
C           2 = Actual memory requirements could not be
C               computed due to limited range of INTEGER
C               type. Treat same as ISTAT=1.
C      LDC    = Leading dimension of orbital coeffitients
C               matrix, which could be different from the
C               number of orbitals to optimize cache use
C      ITOP   = On successful return, set to the maximum
C               number of memory words required during all
C               stages
C      IDSK   = On successful return, set to the number
C               of disk words required during all stages.
C      TRIAL  = .TRUE. if called from strategy dispatcher,
C               and should produce minumum of output
C 
C   Accessed common blocks:
C
C      PSDOPT - Computation options 
C      PSDGBL - Global computation parameters
C      PSDGB2 - Global computation parameters used by response
C               computations.
C      PSPRT  - Printing unit
C      PSPRTF - Debug output tuning flags
C      PSHALF - Half-electron constants
C      PSOCC  - Orbital occupation numbers
C
C   Modified common blocks:
C
C      PSDYNM - IDs of dynamic memory blocks are stored here
C
C   Local storage:
C
C      21 INTEGER cells for block lifetime records
C
C   Module logic:
C
C      This module does not actually touch any memory it
C      allocates. It could be called with different computation
C      options to see how much memory this computation mode
C      would require and (probably) re-adjust them.
C      To help in locating references to unallocated memory
C      blocks, unused memory handles are initialized with sequential
C      small integers.
C
C      See separate documentation for more details.
C
C   Bugs:
C
C     See notes on IHUGE in PSMMGR.
C     Trapping of overflow during computation of memory requirements
C     could fail allocation even if result is still withing allowed
C     range of values for INTEGER type. This should be of no concern
C     unless you really expect to be able to allocate arrays of
C     IHUGE DOUBLE PRECISION variables.
C
C     Forbidden combinations of options could result in chunks of
C     memory being needlessly allocated.
C
C   Notes:
C
C     Setting MEMLAX parameter at the call to PSMINI to the non-zero
C     value might result in less memory-compaction overhead at the
C     expence of increased memory use.
C
C     Conditionals are *intentionally* non-fused in the body of the
C     routine, since this should pose no problem for a decent compiler,
C     but eases subsequent code modifications substantially, as well
C     as simplifying debugging.
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NUMSTG=9)
      PARAMETER (NUMMEM=75)
C
C-AK  PARAMETER (IHUGE=2147483647)
C-AK  PARAMETER (IHUGE=9223372036854775807)
      PARAMETER (IHUGE=HUGE(1))
      PARAMETER (MAXGRP=10)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
      PARAMETER (MAXCIO=LMACT)
C
      LOGICAL UHF, HALFEL, DOCI, TRIAL, DORESP, DODIP, LIMAG, DOPTCH
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
      COMMON
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
     ./PSOCC / DOCCA(MAXGRP), DOCCB(MAXGRP), NMGRPA, NMGRPB,
     .         IOPN1A, IOPN1B, IOPNLA, IOPNLB, IMAXBL,
     .         IG1STA(MAXGRP), IGCNTA(MAXGRP),
     .         IG1STB(MAXGRP), IGCNTB(MAXGRP),
     .         IGBASA(MAXGRP,MAXGRP), IGBASB(MAXGRP,MAXGRP)
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN), 
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN), 
     .         IACTI(MAXACT), IACTJ(MAXACT)
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSPRTF/, /PSHALF/, /PSOCC/,
     .     /PSDCIP/, /PSPRT /
C
C     Do little horrid equivalences, so that we can allocate blocks
C     uniformly after collecting all the lifetime data piecewise.
C     Zeroth element of each column contain desired size of the
C     block, elements 1 to NUMSTG are the active stages.
C
      DIMENSION LALL(NUMMEM), IACTBL(0:NUMSTG,NUMMEM)
      DIMENSION IDXPA (0:NUMSTG), IDXPB (0:NUMSTG), IDXPDA(0:NUMSTG)
      DIMENSION IDXPAT(0:NUMSTG), IDXPBT(0:NUMSTG), IDPA  (0:NUMSTG)
      DIMENSION IDPB  (0:NUMSTG), IHIDAT(0:NUMSTG), ICCNTA(0:NUMSTG)
      DIMENSION ICCNTB(0:NUMSTG), IPCNTA(0:NUMSTG), IPCNTB(0:NUMSTG)
      DIMENSION IAI2T1(0:NUMSTG), IAI2T2(0:NUMSTG), IAI2T3(0:NUMSTG)
      DIMENSION IAI2T4(0:NUMSTG), IAI2T5(0:NUMSTG), IKS   (0:NUMSTG)
      DIMENSION IQ    (0:NUMSTG), IQTMP (0:NUMSTG), IXTMP (0:NUMSTG)
      DIMENSION IIDSPS(0:NUMSTG), ICPITM(0:NUMSTG), IKA   (0:NUMSTG)
      DIMENSION IKATMP(0:NUMSTG), IA2MTM(0:NUMSTG), IA2MT1(0:NUMSTG)
      DIMENSION ICOND (0:NUMSTG), IO2ATM(0:NUMSTG), IO2AT1(0:NUMSTG)
      DIMENSION ICPIT1(0:NUMSTG), ICPAYA(0:NUMSTG), ICPAYB(0:NUMSTG)
      DIMENSION ICPAY (0:NUMSTG), ICPATA(0:NUMSTG), ICPATB(0:NUMSTG)
      DIMENSION ISHIFT(0:NUMSTG), ICNTMP(0:NUMSTG), IHLO  (0:NUMSTG)
      DIMENSION IHLQ  (0:NUMSTG), IHLZM (0:NUMSTG), IHLZA (0:NUMSTG)
      DIMENSION IHLQEX(0:NUMSTG), IC1YY (0:NUMSTG), IC1YR (0:NUMSTG)
      DIMENSION IC1RR (0:NUMSTG), IC1RB (0:NUMSTG), IC1YB (0:NUMSTG)
      DIMENSION IC1BB (0:NUMSTG), IC1RHS(0:NUMSTG), IC1TMA(0:NUMSTG)
      DIMENSION IC1TMX(0:NUMSTG), IC1X  (0:NUMSTG), IC1BAS(0:NUMSTG)
      DIMENSION IC1RSP(0:NUMSTG), IC1TMP(0:NUMSTG), IC1SNG(0:NUMSTG)
      DIMENSION IC1SVT(0:NUMSTG), IC1SCR(0:NUMSTG), IHA0  (0:NUMSTG)
      DIMENSION IHAB  (0:NUMSTG), IH0B  (0:NUMSTG), ICIFI (0:NUMSTG)
      DIMENSION ICIFC (0:NUMSTG), ICIEI (0:NUMSTG), ICIEC (0:NUMSTG)
      DIMENSION ICIVEC(0:NUMSTG), ICIORB(0:NUMSTG), ICIGAM(0:NUMSTG)
      DIMENSION ICIDLT(0:NUMSTG), ICI2ES(0:NUMSTG), ICIH  (0:NUMSTG)
      DIMENSION ICIVC1(0:NUMSTG), ICI2TM(0:NUMSTG), ICIINT(0:NUMSTG)
C
      EQUIVALENCE (LALL(1),LDXPA)
      EQUIVALENCE (IDXPA (0),IACTBL(0, 1)), (IDXPB (0),IACTBL(0, 2))
      EQUIVALENCE (IDXPDA(0),IACTBL(0, 3)), (IDXPAT(0),IACTBL(0, 4))
      EQUIVALENCE (IDXPBT(0),IACTBL(0, 5)), (IDPA  (0),IACTBL(0, 6))
      EQUIVALENCE (IDPB  (0),IACTBL(0, 7)), (IHIDAT(0),IACTBL(0, 8))
      EQUIVALENCE (ICCNTA(0),IACTBL(0, 9)), (ICCNTB(0),IACTBL(0,10))
      EQUIVALENCE (IPCNTA(0),IACTBL(0,11)), (IPCNTB(0),IACTBL(0,12))
      EQUIVALENCE (IAI2T1(0),IACTBL(0,13)), (IAI2T2(0),IACTBL(0,14))
      EQUIVALENCE (IAI2T3(0),IACTBL(0,15)), (IAI2T4(0),IACTBL(0,16))
      EQUIVALENCE (IAI2T5(0),IACTBL(0,17)), (IKS   (0),IACTBL(0,18))
      EQUIVALENCE (IQ    (0),IACTBL(0,19)), (IQTMP (0),IACTBL(0,20))
      EQUIVALENCE (IXTMP (0),IACTBL(0,21)), (IIDSPS(0),IACTBL(0,22))
      EQUIVALENCE (ICPITM(0),IACTBL(0,23)), (IKA   (0),IACTBL(0,24))
      EQUIVALENCE (IKATMP(0),IACTBL(0,25)), (IA2MTM(0),IACTBL(0,26))
      EQUIVALENCE (IA2MT1(0),IACTBL(0,27)), (ICOND (0),IACTBL(0,28))
      EQUIVALENCE (IO2ATM(0),IACTBL(0,29)), (IO2AT1(0),IACTBL(0,30))
      EQUIVALENCE (ICPIT1(0),IACTBL(0,31)), (ICPAYA(0),IACTBL(0,32))
      EQUIVALENCE (ICPAYB(0),IACTBL(0,33)), (ICPAY (0),IACTBL(0,34))
      EQUIVALENCE (ICPATA(0),IACTBL(0,35)), (ICPATB(0),IACTBL(0,36))
      EQUIVALENCE (ISHIFT(0),IACTBL(0,37)), (ICNTMP(0),IACTBL(0,38))
      EQUIVALENCE (IHLO  (0),IACTBL(0,39)), (IHLQ  (0),IACTBL(0,40))
      EQUIVALENCE (IHLZM (0),IACTBL(0,41)), (IHLZA (0),IACTBL(0,42))
      EQUIVALENCE (IHLQEX(0),IACTBL(0,43)), (IC1YY (0),IACTBL(0,44)) 
      EQUIVALENCE (IC1YR (0),IACTBL(0,45)), (IC1RR (0),IACTBL(0,46)) 
      EQUIVALENCE (IC1RB (0),IACTBL(0,47)), (IC1YB (0),IACTBL(0,48))
      EQUIVALENCE (IC1BB (0),IACTBL(0,49)), (IC1RHS(0),IACTBL(0,50)) 
      EQUIVALENCE (IC1TMA(0),IACTBL(0,51)), (IC1TMX(0),IACTBL(0,52)) 
      EQUIVALENCE (IC1X  (0),IACTBL(0,53)), (IC1BAS(0),IACTBL(0,54))
      EQUIVALENCE (IC1RSP(0),IACTBL(0,55)), (IC1TMP(0),IACTBL(0,56)) 
      EQUIVALENCE (IC1SNG(0),IACTBL(0,57)), (IC1SVT(0),IACTBL(0,58)) 
      EQUIVALENCE (IC1SCR(0),IACTBL(0,59)), (IHA0  (0),IACTBL(0,60))
      EQUIVALENCE (IHAB  (0),IACTBL(0,61)), (IH0B  (0),IACTBL(0,62))
      EQUIVALENCE (ICIFI (0),IACTBL(0,63)), (ICIFC (0),IACTBL(0,64))
      EQUIVALENCE (ICIEI (0),IACTBL(0,65)), (ICIEC (0),IACTBL(0,66))
      EQUIVALENCE (ICIVEC(0),IACTBL(0,67)), (ICIORB(0),IACTBL(0,68))
      EQUIVALENCE (ICIGAM(0),IACTBL(0,69)), (ICIDLT(0),IACTBL(0,70))
      EQUIVALENCE (ICI2ES(0),IACTBL(0,71)), (ICIH  (0),IACTBL(0,72))
      EQUIVALENCE (ICIVC1(0),IACTBL(0,73)), (ICI2TM(0),IACTBL(0,74))
      EQUIVALENCE (ICIINT(0),IACTBL(0,75))
C
C     Make all blocks unallocated and clear list of active stages.
C
      DO 10 I=1,NUMMEM
          LALL(I) = -I
          DO 5 J=0,NUMSTG
              IACTBL(J,I) = 0
    5     CONTINUE
   10 CONTINUE
C
C     Initialize memory allocator
C
      ISTAT = 0
      ITOP  = 0
      IDSK  = 0
      IF( LPMMGR ) THEN
          CALL PSMINI(0,5)
      ELSE IF( IPRINT.LT.0 ) THEN
          CALL PSMINI(0,-1)
      ELSE
          CALL PSMINI(0,0)
      ENDIF
      IF(DORESP) THEN
C
C    Check for possible arifmetic overflows and return allocation failure
C
          IF( MODE.NE.10 .AND. (IMIX.EQ.1 .OR. IMIX.EQ.2) ) THEN
              IF( IHUGE/NROS .LE. NVARS ) THEN
                  IF(LPMMGR) WRITE(NB6,11000) DBLE(NROS)*NVARS
                  ISTAT = 2
                  RETURN
              ENDIF
          ENDIF
          IF( MODE.EQ.2 .AND. IHUGE/NROS .LE. IDSTRP ) THEN
              IF(LPMMGR) WRITE(NB6,11050) DBLE(NROS)*IDSTRP
              ISTAT = 2
              RETURN
          ENDIF
          IF( IDENS.EQ.3 .OR. IDENS.EQ.4 .OR. IDENS.EQ.5 ) THEN
              IF( IHUGE/(IQSZ+1)-1 .LE. IQSZ ) THEN
                  IF(LPMMGR) WRITE(NB6,11200) (DBLE(IQSZ)*(IQSZ+1))
                  ISTAT = 2
                  RETURN
              ENDIF
          ENDIF
          IF( IDENS.GE.3 ) THEN
              IF( IHUGE/IQSZ .LE. NCPVRS ) THEN
                  IF(LPMMGR) WRITE(NB6,11300) DBLE(IQSZ)*NVARS
                  ISTAT = 2
                  RETURN
              ENDIF
          ENDIF
          IF( DOPTCH .AND. IHUGE/(3*NPAIR) .LE. NPTCHG ) THEN
              IF(LPMMGR) WRITE(NB6,11310) DBLE(3*NPAIR)*NPTCHG
              ISTAT = 2
              RETURN
          ENDIF
C
C    Compute disk space requirements (disk temporaries may be used
C    for fully-analytic force constants, half-electron derivatives 
C    and NMR chemical shifts).
C
          IF( IDENS.GE.3 ) THEN
              IF( MODE.NE.10 .AND. IMIX.EQ.2 ) THEN
                  IDSK = IDSK + NVARS*NROS
                  IF(UHF) THEN
                      IF( IHUGE-IDSK.LE.NVARS*NROS ) THEN
                          IF(LPMMGR) WRITE(NB6,11400) DBLE(NROS)*NVARS
                          ISTAT = 2
                          RETURN
                      ENDIF
                      IDSK = IDSK + NVARS*NROS
                  ENDIF
              ENDIF
              IF( IQSWAP.NE.-1 .AND. IDENS.NE.8 ) THEN
                  IF( IHUGE-IDSK.LE.NCPVRS*IQSZ ) THEN
                      IF(LPMMGR) WRITE(NB6,11500) DBLE(IQSZ)*NCPVRS
                      ISTAT = 2
                      RETURN
                  ENDIF
                  IDSK = IDSK + NCPVRS*IQSZ
              ENDIF
              IF( IDENS.EQ.5 .AND. IROWS.LT.IQSZ ) THEN
                  CALL PSBPSZ(ISIZE,IKSTRP,ISTAT)
                  IF( ISTAT.NE.0 ) THEN
                      IF(LPMMGR) WRITE(NB6,11550)
                      ISTAT = 2
                      RETURN
                  ENDIF
                  IF( IPRINT.GE.1 .AND. .NOT. TRIAL ) THEN
                      ISIZEU = (IQSZ*(IQSZ+1))/2
                      WRITE(NB6,11560) ISIZE - ISIZEU, 
     .                               DBLE(ISIZE-ISIZEU)*100D0/ISIZEU
                  ENDIF
                  IF( IHUGE-IDSK.LE.ISIZE ) THEN
                      IF(LPMMGR) WRITE(NB6,11600) DBLE(ISIZE)
                      ISTAT = 2
                      RETURN
                  ENDIF
                  IDSK = IDSK + ISIZE
              ENDIF
          ENDIF
      ENDIF
C
C    We'll have to compute second derivatives and therefore need
C    some temporaries. If you think the following code is horrible,
C    please suggest better way of doing it...
C
C LHLO
      IF( HALFEL.OR.DOCI ) THEN
          IF(HALFEL) IHLO(0) = NUMOIJ
          IF(DOCI)   IHLO(0) = (NCIO*(NCIO+1))/2
          IF( IHLST.EQ.2 ) IHLO(0) = IHLO(0) * NPAIR
C        45 is actually 9*(9+1)/2, i.e. number of orbital pairs in MNDO/d
          IF( IHLST.EQ.1 ) IHLO(0) = IHLO(0) * 45 * 2
          IHLO(1) = 1
      ENDIF
C
      IF(DORESP) THEN
C LDXPA, LDXPB
          IF( MODE.NE.10 .AND. IMIX.NE.3 ) THEN
              IDXPA(0) = NROS*NVARS
              IDXPA(1) = 1
              IDXPA(9) = 1
              IF( MODE.GE.2 .AND. IQSWAP.NE.3 .AND. IDENS.NE.1 ) 
     .            IDXPA(2) = 1
              IF( IMIX.EQ.1 ) THEN
                  DO 100 I=2,8
                      IDXPA(I) = 1
  100             CONTINUE
              ENDIF
              IF(UHF) THEN
                  DO 110 I=0,NUMSTG
                      IDXPB(I) = IDXPA(I)
  110             CONTINUE
              ENDIF
          ENDIF
C LDXPDA
          IF( MODE.NE.10 .AND. IMIX.EQ.3 ) THEN
              IF(.NOT.UHF) THEN
                  IDXPDA(0) = 3*(NORBS*(NORBS+1))/2 + 3*(NATOM-2)*NPAIR
              ELSE
                  IDXPDA(0) = 6*(NORBS*(NORBS+1))/2 + 3*(NATOM-3)*NPAIR
              ENDIF
C            Storage for contributions from point charges
              IF( NPTCHG.GT.0 ) THEN
                  IF(DOPTCH) THEN
                      IDXPDA(0) = IDXPDA(0) + NPAIR*NPTCHG*3
                  ELSE
                      IDXPDA(0) = IDXPDA(0) + NPAIR*3
                  ENDIF
              ENDIF
              DO 150 I=1,NUMSTG
                  IDXPDA(I) = 1
  150         CONTINUE
          ENDIF
C LDXPAT, LDXPBT
          IF( MODE.NE.10 .AND. (IMIX.EQ.3 .OR. IMIX.EQ.4) ) THEN
              IDXPAT(0) = NROS
              IDXPAT(1) = 1
              IF( MODE.GE.2 ) THEN
                  IF( IDENS.NE.1 .AND. IQSWAP.NE.3  ) IDXPAT(2) = 1
              ENDIF
              IDXPAT(9) = 1
              IF(UHF) THEN
                  IF (IMIX.EQ.3) THEN
                      IDXPBT = IDXPAT
                  ELSE
                      IDXPAT = 2 * IDXPAT
                  END IF
              ENDIF
          ENDIF
C LDPA, LDPB - real perturbation
          IF( MODE.EQ.2 ) THEN
              IDPA(0) = NROS*IDSTRP
              IDPA(9) = 1
              IF(UHF) THEN
                  DO 250 I=0,NUMSTG
                      IDPB(I) = IDPA(I)
  250             CONTINUE
              ENDIF
          ENDIF
C LDPA, LDPB - imaginary perturbation, only RHF case is currently supported
          IF( MODE.EQ.10 ) THEN
              IDPA(0) = NROS*3
              IDPA(8) = 1
              IDPA(9) = 1
          ENDIF
C LHIDAT, LCCNTA, LCCNTB, LPCNTA, LPCNTB
          IF( IDENS.EQ.1 ) THEN
              IHIDAT(0) = 2*(NPAIR*(NPAIR+2*NATOM-4)-NPAIR2+
     .                    NORBS*(NORBS+1))
              DO 300 I=1,NUMSTG
                  IHIDAT(I) = 1
  300         CONTINUE
              ICCNTA(0) = NORBS*LDC
              ICCNTA(9) = 1
              IPCNTA(0) = NROS
              IPCNTA(9) = 1
              IF(UHF) THEN
                  DO 310 I=0,NUMSTG
                       ICCNTB(I) = ICCNTA(I)
                       IPCNTB(I) = IPCNTA(I)
  310             CONTINUE
              ENDIF
          ENDIF
C LAI2T1
          IF( IDENS.NE.1 ) THEN
              IAI2T1(0) = (NPAIR**2 + NPAIR2)/2
              IAI2T1(1) = 1
              IF( IDENS.EQ.7 ) THEN
                  DO 350 I=2,7
                      IAI2T1(I) = 1
  350         CONTINUE
              ENDIF
              IF( IDENS.EQ.8 ) THEN
                  DO 360 I=2,9
                      IAI2T1(I) = 1
  360         CONTINUE
              ENDIF
              IF( IDENS.EQ.3 .OR. IDENS.EQ.4 .OR. IDENS.EQ.5 ) THEN
                  DO 370 I=2,5
                      IAI2T1(I) = 1
  370         CONTINUE
              ENDIF
              IF( IDENS.EQ.6 ) THEN
                  DO 380 I=2,7
                      IAI2T1(I) = 1
  380         CONTINUE
              ENDIF
              IF( (IPRECT.EQ.3 .OR. IPRECT.EQ.5 .OR. IPRECT.EQ.6) .AND. 
     .                                               IDENS.GE.4 ) THEN
                  DO 390 I=2,4
                      IAI2T1(I) = 1
  390         CONTINUE
              ENDIF
              IF( HALFEL.OR.DOCI ) THEN
                  IAI2T1(2) = 1
                  IAI2T1(3) = 1
              ENDIF
          ENDIF
C LAI2T2, LAI2T3, LAI2T4, LAI2T5
          IAI2T2(0) = NORBS*NPAIR
          IAI2T3(0) = NPAIR*NORBS
          IAI2T4(0) = 0
          IAI2T5(0) = 0
          IF(HALFEL) THEN
              IAI2T4(0) = NUMOPN
              IAI2T5(0) = NUMOPN*NORBS
          ENDIF
          IF(DOCI) THEN
              IAI2T4(0) = NCIO
              IAI2T5(0) = NCIO*NORBS
          ENDIF
          DO 410 I=1,NMGRPA-1
              IAI2T4(0) = MAX(IAI2T4(0),IGCNTA(I))
              DO 400 J=I+1,NMGRPA
                  IAI2T5(0) = MAX(IAI2T5(0),IGCNTA(I)*IGCNTA(J))
  400         CONTINUE
  410     CONTINUE
          IF(UHF) THEN
              DO 430 I=1,NMGRPB-1
                  IAI2T4(0) = MAX(IAI2T4(0),IGCNTB(I))
                  DO 420 J=I+1,NMGRPB
                      IAI2T5(0) = MAX(IAI2T5(0),IGCNTB(I)*IGCNTB(J))
  420             CONTINUE
  430         CONTINUE
          ENDIF
          IAI2T4(0) = IAI2T4(0)*NORBS
          IF( IKMODE.EQ.1 .AND. (IDENS.EQ.3 .OR. IDENS.EQ.4 
     .                                      .OR. IDENS.EQ.5) ) THEN
              IAI2T2(5) = 1
              IAI2T3(5) = 1
              IAI2T4(5) = 1
              IAI2T5(5) = 1
          ENDIF
          IF( IKMODE.EQ.1 .AND. IDENS.EQ.6 ) THEN
              IAI2T2(7) = 1
              IAI2T3(7) = 1
              IAI2T4(7) = 1
              IAI2T5(7) = 1
          ENDIF
          IF( HALFEL.OR.DOCI ) THEN
              IAI2T2(3) = 1
              IAI2T3(3) = 1
              IAI2T4(3) = 1
              IAI2T5(3) = 1
          ENDIF
          IF( (IPRECT.EQ.3 .OR. IPRECT.EQ.5 .OR. IPRECT.EQ.6) 
     .                                      .AND. IDENS.GE.4 ) THEN
              IAI2T2(4) = 1
              IAI2T3(4) = 1
              IAI2T4(4) = 1
              IAI2T5(4) = 1
          ENDIF
C LKS
          IF( IDENS.EQ.3 .OR. IDENS.EQ.4 ) THEN
              IKS(0) = (IQSZ*(IQSZ+1))/2
              IKS(5) = 1
              IKS(6) = 1
              IKS(7) = 1
          ENDIF
C LQ 
          IF( IDENS.NE.1 .AND. IQSWAP.LT.2 ) THEN
              IQ(0) = IQSZ*NCPVRS
              IQ(2) = 1
              DO 500 I=7,9
                  IQ(I) = 1
  500         CONTINUE
              IF( MODE.EQ.10 ) IQ(9) = 0
              IF( IQSWAP.EQ.-1 ) THEN
                  DO 510 I=3,6
                      IQ(I) = 1
  510             CONTINUE
              ENDIF
          ENDIF
C LQTMP
          IF( IQSWAP.GE.0 .AND. IDENS.GE.2 ) THEN
              IQTMP(0) = IQSZ
              IF( IQSWAP.NE.3 ) IQTMP(2) = 1
              IF( IDENS .EQ.8 ) IQTMP(9) = 1
              IF( HALFEL.OR.DOCI ) IQTMP(3) = 1
          ENDIF
C LXTMP
          IF( IQSWAP.GE.2 .AND. IDENS.GE.2 ) THEN
              IXTMP(0) = IQSZ
              IF( MODE.NE.10 ) IXTMP(9) = 1
              IF( MODE.EQ.10 ) IXTMP(8) = 1
          ENDIF
C LIDSPS
          IF( IDENS.EQ.3 ) THEN
C-AK          IIDSPS(0) = (IQSZ+1)/2
              IIDSPS(0) = IQSZ
              IIDSPS(7) = 1
          ENDIF
C LCPITM, LCPIT1
          IF( IDENS.GE.4 .AND. ISOLVE.EQ.1 ) THEN
              ICPITM(0) = INRHS*(IMAXIT+1)*(IMAXIT+3) + IKRVEC*IQSZ
C-AK          ICPIT1(0) = IMAXIT**2 + (3*IMAXIT+1)/2
              ICPIT1(0) = 2*IMAXIT**2 + 3*IMAXIT
              IF( IDENS.EQ.8 ) THEN
                  ICPITM(9) = 1
                  ICPIT1(9) = 1
              ELSE
                  ICPITM(7) = 1
                  ICPIT1(7) = 1
              ENDIF
          ENDIF
C LC1YY, LC1YR, LC1RR, LC1RB, LC1YB, LC1BB, LC1RHS, LC1TMA, LC1TMX, 
C LC1X, LC1BAS, LC1RSP, LC1TMP, LC1SNG, LC1SVT, LC1SCR 
          IF( IDENS.GE.4 .AND. ISOLVE.GE.2 ) THEN
              IC1YY (0) = INRHS**2
              IC1YR (0) = INRHS*IKRVEC
              IC1RR (0) = IKRVEC**2
              IC1RB (0) = IKRVEC**2
              IC1YB (0) = INRHS*IKRVEC
              IC1BB (0) = IKRVEC**2
              IC1TMA(0) = IKRVEC**2
              IC1TMX(0) = INRHS*IKRVEC
              IF( ISOLVE.EQ.2 .OR. ISOLVE.EQ.3 ) THEN
                  IC1X  (0) = IQSZ
              ENDIF
              IF( ISOLVE.EQ.4 .OR. ISOLVE.EQ.5 .OR. ISOLVE.EQ.6 ) THEN
                  IC1X  (0) = IQSZ*INRHS
              ENDIF
              IF( ISOLVE.EQ.2 .OR. ISOLVE.EQ.3 .OR. ISOLVE.EQ.4 
     .                                         .OR. ISOLVE.EQ.6 ) THEN
                  IC1BAS(0) = IQSZ*(IKRVEC+INRHS)
              ENDIF
              IF( ISOLVE.EQ.5 ) THEN
                  IC1BAS(0) = IQSZ*IKRVEC
              ENDIF
              IC1RSP(0) = IQSZ*IKRVEC
C-AK          IC1TMP(0) = (IKRVEC+1)/2
              IC1TMP(0) = IKRVEC
              IC1SNG(0) = INRHS
              IC1SVT(0) = 1024 + 3*INRHS**2 + MAX(3*INRHS+IQSZ,5*IQSZ-4,
     .                        2*INRHS**2,INRHS**2+INRHS*IQSZ+INRHS)
              IC1SCR(0) = IKRVEC
              IF( IQSWAP.GT.1 ) THEN
                  IC1RHS(0) = IQSZ*INRHS
              ENDIF
              IF( IDENS.EQ.8 ) THEN
                  ISTAGE = 9 
              ELSE
                  ISTAGE = 7
              ENDIF
              IC1YY (ISTAGE) = 1
              IC1YR (ISTAGE) = 1
              IC1RR (ISTAGE) = 1
              IC1RB (ISTAGE) = 1
              IC1YB (ISTAGE) = 1
              IC1BB (ISTAGE) = 1
              IC1RHS(ISTAGE) = 1
              IC1TMA(ISTAGE) = 1
              IC1TMX(ISTAGE) = 1
              IC1X  (ISTAGE) = 1
              IC1BAS(ISTAGE) = 1
              IC1RSP(ISTAGE) = 1
              IC1TMP(ISTAGE) = 1
              IC1SNG(ISTAGE) = 1
              IC1SVT(ISTAGE) = 1
              IC1SCR(ISTAGE) = 1
          ENDIF
C LKATMP
          IF( IDENS.EQ.5 .OR. IDENS.EQ.6 ) THEN
              IKATMP(0) = IROWS*IQSZ
              IKATMP(7) = 1
              IF( IDENS.EQ.5 ) THEN
                  IKATMP(5) = 1
                  IF( IROWS.GE.IQSZ ) IKATMP(6) = 1
              ENDIF
          ENDIF
C LA2MTM
          IF( IDENS.GE.2 ) THEN
              IA2MTM(0) = IMAXBL*NORBS
              IF( IKMODE.EQ.2 ) THEN
                  IF( IDENS.EQ.3 .OR. IDENS.EQ.4 .OR. IDENS.EQ.5 ) 
     .                             IA2MTM(5) = 1
                  IF( IDENS.EQ.6 ) IA2MTM(7) = 1
              ENDIF
              IF( IQSWAP.NE.3 ) IA2MTM(2) = 1
              IF( IDENS .EQ.7 ) IA2MTM(7) = 1
              IF( IDENS .EQ.8 ) IA2MTM(9) = 1
              IF( IHLWRP.EQ.2 .AND. HALFEL.OR.DOCI ) IA2MTM(3) = 1
          ENDIF
C LA2MT1
          IF( IDENS.GE.2 .AND. MODE.GE.2 ) THEN
              IA2MT1(0) = NORBS**2
              IF( IQSWAP.NE.3 ) IA2MT1(2) = 1
              IF( IQSWAP.EQ.3 ) IA2MT1(9) = 1
          ENDIF
C LCOND
          IF( IDENS.GE.4 ) THEN
              ICOND(0) = IQSZ
              ICOND(6) = 1
              ICOND(7) = 1
              IF( IDENS.EQ.8 ) THEN
                  ICOND(8) = 1
                  ICOND(9) = 1
              ENDIF
          ENDIF
C LO2ATM
          IF( IDENS.NE.1 ) THEN
              IO2ATM(0) = IMAXBL*NORBS
              IF( IDENS.EQ.7 ) IO2ATM(7) = 1
              IF( IDENS.EQ.8 ) IO2ATM(9) = 1
              IF( HALFEL.OR.DOCI ) IO2ATM(8) = 1
              IF( MODE.EQ.2 ) IO2ATM(9) = 1
              IF( MODE.EQ.10 ) IO2ATM(8) = 1
          ENDIF
C LO2AT1
          IF( IDENS.NE.1 .AND. MODE.GE.2 ) THEN
              IO2AT1(0) = NORBS**2
              IF( MODE.EQ.2  ) IO2AT1(9) = 1
              IF( MODE.EQ.10 ) IO2AT1(8) = 1
          ENDIF
C LCPAYA, LCPAYB, LCPAY, LCPATA, LCPATB
          ICPAYA(0) = NORBS**2
          IF( IDENS.EQ.7 ) ICPAYA(7) = 1
          IF( IDENS.EQ.8 ) ICPAYA(9) = 1
          IF( IHLWRP.EQ.2 .AND. HALFEL.OR.DOCI ) ICPAYA(3) = 1
          IF( IKMODE.EQ.2 ) THEN
              IF( IDENS.EQ.3 .OR. IDENS.EQ.4 .OR. IDENS.EQ.5 ) 
     .            ICPAYA(5) = 1
              IF( IDENS.EQ.6 ) ICPAYA(7) = 1
          ENDIF
          DO 700 I=0,NUMSTG
              ICPAY (I) = ICPAYA(I)
              ICPATA(I) = ICPAYA(I)
  700     CONTINUE
          IF(UHF) THEN
              DO 710 I=0,NUMSTG
                  ICPAYB(I) = ICPAYA(I)
                  ICPATB(I) = ICPAYA(I)
  710         CONTINUE
          ENDIF
C LSHIFT
          IF( IPRECT.GE.2 .AND. IDENS.GE.4 ) THEN
              ISHIFT(0) = IQSZ
              IF( IPRECT.NE.4 ) ISHIFT(4) = 1
              ISHIFT(5) = 1
              ISHIFT(6) = 1
              IF( IDENS.EQ.6 .OR. IDENS.EQ.7 ) ISHIFT(7) = 1
              IF( IDENS.EQ.8 ) THEN
                  DO 750 I=7,9
                      ISHIFT(I) = 1
  750             CONTINUE
              ENDIF
          ENDIF
C LCNTMP
          IF( IPRECT.GE.2 .AND. IDENS.GE.4 ) THEN
              ICNTMP(0) = IQSZ * INRHS
              IF( IDENS.NE.8 ) ICNTMP(7) = 1
              IF( IDENS.EQ.8 ) ICNTMP(9) = 1
          ENDIF
C LHLQ
          IF( HALFEL.OR.DOCI ) THEN
              IF(HALFEL) IHLQ(0) = IQSZ + NMRACT
              IF(DOCI)   IHLQ(0) = IQSZ + NCIO*NORBS
              DO 800 I=3,8
                  IHLQ(I) = 1
  800         CONTINUE
          ENDIF
C LHLZM
          IF( HALFEL.OR.DOCI ) THEN
              IHLZM(0) = NORBS**2
              IHLZM(8) = 1
          ENDIF
C LHLZA
          IF( HALFEL.OR.DOCI ) THEN
              IHLZA(0) = NROS
              IHLZA(8) = 1
              IHLZA(9) = 1
          ENDIF
C LHLQEX
          IF( HALFEL.OR.DOCI ) THEN
              IF(HALFEL) IHLQEX(0) = NUMOPN*NORBS
              IF(DOCI)   IHLQEX(0) = NCIO*NORBS
              IHLQEX(3) = 1
          ENDIF
      ENDIF
C
C    NMR-specific quantities
C
      IF( MODE.EQ.10 ) THEN
          IHA0(0) = 3*(NORBS**2)
          IHA0(2) = 1
          IHAB(0) = 9*(NORBS**2)
          IHAB(9) = 1
          IH0B(0) = 3*(NORBS**2)
          IH0B(9) = 1
      ENDIF
C
C    CI-specific quantities. Some quantities are also shared
C    with half-electron data, see above.
C
      IF(DOCI) THEN
C-AK      ICIFI (0) = (NCIFS+1)/2
          ICIFI (0) = NCIFS
          ICIFC (0) = NCIFS
C-AK      ICIEI (0) = (NCIES+1)/2
          ICIEI (0) = NCIES
          ICIEC (0) = NCIES
          DO 1190 I=1,3
              ICIFI(I) = 1
              ICIFC(I) = 1
              ICIEI(I) = 1
              ICIEC(I) = 1
 1190     CONTINUE
          ICIVEC(0) = NCIC
          ICIVEC(1) = 1
          ICIVEC(2) = 1
          ICIVEC(3) = 1
          ICIORB(0) = NCIO*NORBS
          DO 1200 I=1,3
              ICIORB(I) = 1
 1200     CONTINUE
          ICIGAM(0) = NCIGAM
          ICIGAM(1) = 1
          ICIGAM(2) = 1
          ICIGAM(3) = 1
          ICIDLT(0) = NCIO
          DO 1300 I=1,8
              ICIDLT(I) = 1
 1300     CONTINUE
          ICI2ES(0) = ((NCIO*(NCIO+1))/2)**2
          ICI2ES(1) = 1
C-AK      ICIH  (0) = NCIC**2
          ICIH  (0) = 1
          ICIH  (3) = 1
          ICIVC1(0) = NCIC
          ICIVC1(3) = 1
          ICI2TM(0) = ((NCIO*(NCIO+1))/2)*45
          ICI2TM(1) = 1
          ICIINT(0) = NCIGAM
          ICIINT(3) = 1
      ENDIF
C
C    Perform actual memory allocation
C
      DO 9000 I=1,NUMMEM
          IF( IACTBL(0,I).GT.0 ) THEN
              DO 8990 J=1,NUMSTG
                  IF( IACTBL(J,I).NE.0 ) THEN
                      CALL PSMALC(IACTBL(0,I),NUMSTG,IACTBL(1,I),
     .                            LALL(I))
                      GOTO 8991
                  ENDIF
 8990         CONTINUE
 8991         CONTINUE
          ENDIF
 9000 CONTINUE
C
      IF( IPRINT.GE.5 .AND. .NOT. TRIAL ) THEN
          WRITE(NB6,10000) (LALL(I), I=1,NUMMEM)
          CALL PSMSTS
      ENDIF
C
C    Check whether we'll fit in available memory, or not
C
      ITOP = IPSMUS()
      IF( ITOP.GT.ICORE ) ISTAT = 1
      IF( IDISK.NE.0 .AND. IDSK.GT.IDISK ) ISTAT = 1
      RETURN
10000 FORMAT( ' MEMORY BLOCK(S) ID(S) ARE: '/
     .' LDXPA  = ',I8,' LDXPB  = ',I8,' LDXPDA = ',I8,' LDXPAT = ',I8,
     .' LDXPBT = ',I8/' LDPA   = ',I8,' LDPB   = ',I8,' LHIDAT = ',I8,
     .' LCCNTA = ',I8,' LCCNTB = ',I8/' LPCNTA = ',I8,' LPCNTB = ',I8,
     .' LAI2T1 = ',I8,' LAI2T2 = ',I8,' LAI2T3 = ',I8/' LAI2T4 = ',I8,
     .' LAI2T5 = ',I8,' LKS    = ',I8,' LQ     = ',I8,' LQTMP  = ',I8/
     .' LXTMP  = ',I8,' LIDSPS = ',I8,' LCPITM = ',I8,' LKA    = ',I8,
     .' LKATMP = ',I8/' LA2MTM = ',I8,' LA2MT1 = ',I8,' LCOND  = ',I8,
     .' LO2ATM = ',I8,' LO2AT1 = ',I8/' LCPIT1 = ',I8,' LCPAYA = ',I8,
     .' LCPAYB = ',I8,' LCPAY  = ',I8,' LCPATA = ',I8/' LCPATB = ',I8,
     .' LSHIFT = ',I8,' LCNTMP = ',I8,' LHLO   = ',I8,' LHLQ   = ',I8/
     .' LHLZM  = ',I8,' LHLZA  = ',I8,' LHLQEX = ',I8,' LC1YY  = ',I8,
     .' LC1YR  = ',I8/' LC1RR  = ',I8,' LC1RB  = ',I8,' LC1YBX = ',I8,
     .' LC1BB  = ',I8,' LC1RHS = ',I8/' LC1TMA = ',I8,' LC1TMX = ',I8,
     .' LC1X   = ',I8,' LC1BAS = ',I8,' LC1RSP = ',I8/' LC1TMP = ',I8,
     .' LC1SNG = ',I8,' LC1SVT = ',I8,' LC1SCR = ',I8,' LHA0   = ',I8/
     .' LHAB   = ',I8,' LH0B   = ',I8,' LCIFI  = ',I8,' LCIFC  = ',I8,
     .' LCIEI  = ',I8/' LCIEC  = ',I8,' LCIVEC = ',I8,' LCIORB = ',I8,
     .' LCIGAM = ',I8,' LCIDLT = ',I8/' LCI2ES = ',I8,' LCIH   = ',I8,
     .' LCIVC1 = ',I8,' LCI2TM = ',I8,' LCIINT = ',I8)

C
C    Formats used for debugging printouts
C
11000 FORMAT( ' NROS*NVARS (', E14.7, ') WON''T FIT IN INTEGER' )
11050 FORMAT( ' NROS*IDSTRP (', E14.7, ') WON''T FIT IN INTEGER' )
11100 FORMAT( ' IQSZ**2 (', E14.7, ') WON''T FIT IN INTEGER' )
11200 FORMAT( ' IQSZ*(IQSZ+1) (', E14.7, ') WON''T FIT IN INTEGER' )
11300 FORMAT( ' IQSZ*NCPVRS (', E14.7, ') WON''T FIT IN INTEGER' )
11310 FORMAT( ' NPAIR*NPTCHG*3 (', E14.7, ') WON''T FIT IN INTEGER' )
11400 FORMAT( ' BETA PART OF IUMIX SIZE (', E14.7, ') WOULD OVERFLOW',
     .        ' IDSK' )
11500 FORMAT( ' IURHS SIZE (', E14.7, ') WOULD OVERFLOW IDSK' )
11550 FORMAT( ' SIZE OF IUK FILE WOULD OVERFLOW INTEGER' )
11560 FORMAT( ' ', I8, ' WORDS (', G14.7, '%) WOULD BE WASTED IN ',
     /        'IUK FILE' )
11600 FORMAT( ' IUK SIZE (', E14.7, ') WOULD OVERFLOW IDSK' )
      END
