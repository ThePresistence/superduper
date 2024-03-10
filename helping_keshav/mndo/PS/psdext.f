C     ******************************************************************
C
C     Execution time model for the analytical derivatives.
C
C     ******************************************************************
      SUBROUTINE PSDEXT(XTIME)
C
C   Estimate computation time for the current set of computation
C   options using fitted formulas. 
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      XTIME  - Estimated total execution time. Detailed per-stage
C               execution time is in PSDSTT common.
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGBL - Scratch common for global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C      PSDTMC - Scaling constants for estimation of execution
C               time.
C      PSOCC  - Occupation numbers description
C      PSHALF - Open shells description
C      PSDCIP - CI parameters
C      PSDGB3 - Global computation parameters needed for the execution
C               time estimations only
C
C   Modified common blocks:
C
C      PSDSTT - Expected per-stage execution times, computed 
C               on return.
C
C   Local storage:
C
C   Module logic:
C
C      This is essentially a simulation of PSDRV with linear
C      fitted approximations to the excution times. See 
C      function body for the description of the actual terms
C      computed.
C
C      Meaning of the ATIME constants is (loosely) as follows:
C      (all parameters except ATIME(1) are in BLAS3 flops per
C       operation)
C
C        1   "Machine speed" in seconds per BLAS3 flop
C        2   Constant setup costs
C        4   Computation of the single two-center integral along 
C            with first derivatives
C        5   Computation of the single two-center integral along 
C            with first and second derivatives
C        6   Precomputing single element of Oij (half-electron)
C        7   Recomputing single element of Oij (half-electron)
C        8   Static part of the first derivative of the 
C            half-electron correction
C        9   Recomputing single set of the Fock matrix derivatives 
C            in A.O. basis 
C       10   BLAS3 flop, arbitrary fixed at 1.0
C       11   Spin-independant part of the Fock matrix
C       12   Spin-dependant part of the Fock matrix
C       13   Single integral transform over first index
C       14   Single integral transform over second index
C       15   Single integral transform over third index
C       16   Single integral transform over forth index
C       17   DSYR2 performance on matrices ca. NORBS**2
C       18   Shift vector cropping
C       19   I/O transfer penalty (flops/dble)
C       20   Sqrt costs
C       21   Shift->conditioner costs (ca. division costs)
C       22   Decomposition part of the direct solver
C       23   Substitution part of the direct solver
C       24   DSYMV performance on large vectors
C       25   DSYMM/DGEMM performance on narrow matrices
C       26   DDOT performance
C       27   Conditioning costs (ca. mul + 3 mem)
C       28   DGEMM performance on flat matrices (iterative solver)
C       29   Numeric density matrix derivatives costs
C       30   Scaling coefficient for partial AO->MO transformations.
C       31   Penalty for fetching single DOUBLE PRECISION variable
C            from the main memory.
C       32   Set of DAXPYs implemented via calls to DGEMM/DGEMV
C       33   Costs of SVD
C
C   Bugs:
C
C      Some auxiliary times are computed even if never used due
C      to the combination of options.
C
C      Timing for the CPHF K matrix formation are approximate
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NTIMEC=33)
      PARAMETER (MAXGRP=10)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXCIO=LMACT)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
C
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (TWO=2.0D0)
      PARAMETER (THREE=3.0D0)
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
     ./PSDGB3/ NOPAIR, INPAIR
     ./PSOCC / DOCCA(MAXGRP), DOCCB(MAXGRP), NMGRPA, NMGRPB, 
     .         IOPN1A, IOPN1B, IOPNLA, IOPNLB, IMAXBL,
     .         IG1STA(MAXGRP), IGCNTA(MAXGRP), 
     .         IG1STB(MAXGRP), IGCNTB(MAXGRP),
     .         IGBASA(MAXGRP,MAXGRP), IGBASB(MAXGRP,MAXGRP)
      COMMON
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN), 
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN),
     .         IACTI(MAXACT), IACTJ(MAXACT)
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
     ./PSDTMC/ ATIME(NTIMEC)
     ./PSDSTT/ STIMES(-1:9), XTIMES(-1:9)
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSDTMC/, /PSDSTT/
      SAVE /PSOCC /, /PSHALF/, /PSDGB3/, /PSDCIP/
C
C   Zero out expected execution times, since were is a chance we'll
C   never visit some of the stages.
C
      STIMES(0) = ATIME(2)
      DO 100 I=1,9
          STIMES(I) = ZERO
  100 CONTINUE
C
C   Some time estimations of the general interest...
C
      DNPAIR = DBLE(NPAIR)
      DINPAI = DBLE(INPAIR)
      DNPAI2 = DBLE(NPAIR2)
      DNINTS = HALF*DBLE(DNPAIR**2+DNPAI2)
C
      IF(DORESP) THEN
          DNVARS = DBLE(NVARS)
          DNORBS = DBLE(NORBS)
          DNORB2 = DNORBS**2
          DNATOM = DBLE(NATOM)
          DIQSZ  = DBLE(IQSZ)
          DIQSZA = DBLE(IQSZA)
          DIQSZB = DBLE(IQSZB)
          DNROS  = DBLE(NROS)
C
C   Time for A.O.->M.O. transformation
C
          TA2MA = ATIME(10)*2*DNORBS*
     .                  (DNORBS*(DNORBS-IGCNTA(NMGRPA)) + DIQSZA)
          IF(UHF) THEN
              TA2MB = ATIME(10)*2*DNORBS*
     .                      (DNORBS*(DNORBS-IGCNTB(NMGRPB)) + DIQSZB)
          ELSE
              TA2MB = ZERO
          ENDIF
          TA2M = TA2MA + TA2MB
C
C   Time for M.O.->A.O. transformation
C
          TM2AA = ZERO
          DO 170 I=1,NMGRPA-1
              TM2AA = TM2AA + (NMGRPA-I)*IGCNTA(I)
  170     CONTINUE
          TM2AA = TA2MA + ATIME(10)*2*DNORB2*TM2AA
          TM2AB = ZERO
          IF(UHF) THEN
              DO 180 I=1,NMGRPB-1
                  TM2AB = TM2AB + (NMGRPB-I)*IGCNTB(I)
  180         CONTINUE
              TM2AB = TA2MB + ATIME(10)*2*DNORB2*TM2AB
          ENDIF
          TM2A = TM2AA + TM2AB
C
C    Computation of "Fock" matrix
C
          ADTIM  = ATIME(12)*( (2*DNPAIR-DNORBS)**2 + (2*DNPAIR-DNORBS)
     .      + 2*DNPAI2 - DNORBS - 4*NOPAIR + (1+DNATOM)*DNPAIR - DINPAI)
          TFOCKA = ATIME(11)*2*DNPAIR*(DNPAIR+DNATOM ) + ADTIM
          TFOCK  = TFOCKA
          IF(UHF) TFOCK = TFOCK + ADTIM
C
C    Per-stage integrals transformation times
C
          TINT1  = ATIME(13)*2*(2*DNPAIR-DNORBS)*DNPAIR
          TINT2  = ATIME(14)*2*DNORBS*(DNPAIR+DNATOM)
          TINT3  = ATIME(15)*(4*DNPAIR-2*DNORBS)
          TINT4  = ATIME(16)*2*DNORBS
C
C    DSYR2 cost
C
          TDSYR2 = ATIME(17)*2*DNORBS*(DNORBS+2)
C
C    CPHF K formation cost
C
          IF( IDENS.GE.3 .AND. IDENS.LE.6 ) THEN
               IF( IKMODE.EQ.1 ) THEN
                   TK = PSDEKM(TINT1,TINT2,TINT3,TINT4)
               ELSE
                   TK = PSDEKA(TDSYR2,TFOCKA,TFOCK,TA2MA,ATIME(30))
               ENDIF
          ENDIF
          IF( IDENS.EQ.5 .OR. IDENS.EQ.6 ) THEN
              CALL PSBPSZ(IKSIZE,IKSTRP,IDUMMY)
          ELSE
              IKSIZE = 0
              IKSTRP = 0
          ENDIF
C
C    Fock matrix derivatives recomputation cost
C
          TMIX = ATIME(9)*DNVARS*(DNORB2 + DNPAIR)
          IF(UHF) TMIX = 2*TMIX
      ENDIF
C
C    Stage 1 - Two-center integrals & directly related quantities
C
      IF( HALFEL.OR.DOCI ) THEN
          IF(HALFEL) THEN
              NOIJ = NUMOIJ
              NACT = NUMOPN
          ENDIF
          IF(DOCI) THEN
              NOIJ = (NCIO*(NCIO+1))/2
              NACT = NCIO
          ENDIF
          IF( IHLST.EQ.2 ) THEN
              STIMES(1) = STIMES(1) + ATIME(6)*3*NOIJ*DNPAIR
          ELSE
              STIMES(1) = STIMES(1) + ATIME(7)*3*NOIJ*DINPAI
          ENDIF
          IF(HALFEL) THEN
              ADTIM = NUMOIJ*DNINTS + NUMOPN*(3*NUMOPN-1)*DNPAIR
          ENDIF
          IF(DOCI) THEN
              ADTIM = 4*(NOIJ**2)*DNINTS + NCIGAM*DNATOM*(DNATOM+1)/2
          ENDIF
          ADTIM = ATIME(8)*ADTIM
      ENDIF
      IF( MODE.EQ.1 ) THEN
          STIMES(1) = STIMES(1) + ATIME(4)*DNINTS
          IF( HALFEL.OR.DOCI ) STIMES(1) = STIMES(1) + 3*ADTIM
      ELSE
          STIMES(1) = STIMES(1) + ATIME(5)*DNINTS
          IF( HALFEL.OR.DOCI ) STIMES(1) = STIMES(1) + 9*ADTIM
      ENDIF
      IF( DORESP .AND. IMIX.LE.2 ) THEN
          STIMES(1) = STIMES(1) + TMIX
      ENDIF
C
C    Stage 2 - CPHF equations right-hand sides
C
      IF( DORESP .AND. IQSWAP.NE.3 .AND. MODE.GE. 2 
     .           .AND. IDENS.NE.1 ) THEN
          IF( IMIX.EQ.3 ) THEN
              STIMES(2) = STIMES(2) + TMIX
          ENDIF
          STIMES(2) = STIMES(2) + DNVARS*TA2M
          IF( IQSWAP.GE.0 ) THEN
              STIMES(2) = STIMES(2) + ATIME(19)*DNVARS*DIQSZ
          ENDIF
      ENDIF
C
C     Swapping out of the Fock matrix derivatives
C
      IF(DORESP .AND. IDENS.GE.3 .AND. IDENS.LE.7 .AND. IMIX.EQ.2) THEN
          STIMES(2) = STIMES(2) + ATIME(19)*DNVARS*DNROS
          IF(UHF) STIMES(2) = STIMES(2) + ATIME(19)*DNVARS*DNROS
      ENDIF
C
C    Stage 3 - Half-electron RHS
C
      IF( DORESP .AND. (HALFEL.OR.DOCI) ) THEN
          STIMES(3) = STIMES(3) + NACT*TINT1 + NOIJ*TINT2 + 
     .                NACT*(2*NACT+1)*TINT3 + 
     .                NACT*(2*NACT+1)*DNORBS*TINT4
          IF(HALFEL) IACTX = NACT*(NACT-1)
          IF(DOCI)   IACTX = NACT*NORBS
          IF( IHLWRP.EQ.1 ) THEN
              STIMES(3) = STIMES(3) + (NACT-1)*TINT1 + (NACT-1)*
     .                    DNORBS*TINT2 + (THREE/TWO)*IACTX*( DNORBS - 
     .                    IGCNTA(NMGRPA) )*TINT3 + (THREE/TWO)*IACTX
     .                    *IQSZA*TINT4
          ELSE
              IF( IACTX .GT. 0 ) THEN
                  STIMES(3) = STIMES(3) + IACTX*TDSYR2 + TFOCKA + TA2MA
              ENDIF
          ENDIF
          IF( IQSWAP.GT.0 ) THEN
              STIMES(3) = STIMES(3) + ATIME(19)*DIQSZ
          ENDIF
      ENDIF
C
C    Stage 4 - Shift vector
C
      IF( DORESP .AND. IPRECT.GT.1 .AND. (IPRECT.NE.4 .AND. (IDENS.EQ.4
     .                         .OR. IDENS.EQ.5) .OR. IDENS.GT.5) ) THEN
          STIMES(4) = STIMES(4) + ATIME(18) * DIQSZ
          IF( IPRECT.EQ.6 ) THEN
              STIMES(4) = STIMES(4) + 
     .                    (DNORBS - IGCNTA(NMGRPA)) * (TINT1+TINT2) +
     .                    DIQSZA * (TINT3+TINT4)
              IF(UHF) STIMES(4) = STIMES(4) + 
     .                    (DNORBS - IGCNTB(NMGRPB)) * (TINT1+TINT2) +
     .                    DIQSZB * (TINT3+TINT4)
          ENDIF
          IF( IPRECT.EQ.3 ) THEN
              STIMES(4) = STIMES(4) + 
     .                NMGRPA*(NMGRPA+1)*(TINT1+2*TINT2+2*TINT3+2*TINT4)
              IF(UHF) STIMES(4) = STIMES(4) + 
     .                NMGRPB*(NMGRPB+1)*(TINT1+2*TINT2+2*TINT3+2*TINT4)
          ENDIF
          IF( IPRECT.EQ.5 ) THEN
              STIMES(4) = STIMES(4) + (THREE/TWO)*
     .                NMGRPA*(NMGRPA+1)*(TINT1+2*TINT2+2*TINT3+2*TINT4)
              IF(UHF) STIMES(4) = STIMES(4) + (THREE/TWO)*
     .                NMGRPB*(NMGRPB+1)*(TINT1+2*TINT2+2*TINT3+2*TINT4)
          ENDIF
      ENDIF
C
C    Stage 5 - CPHF K matrix
C
      IF( DORESP .AND. IDENS.GE.3 .AND. IDENS.LT.6 ) THEN
           STIMES(5) = STIMES(5) + TK
           IF( IDENS.EQ.5 .AND. IROWS.LT.IQSZ ) THEN
               STIMES(5) = STIMES(5) + IKSIZE*ATIME(19)
           ENDIF
      ENDIF
C
C    Stage 6 - Preconditioner vector
C
      IF( DORESP .AND. IDENS.GT.3 ) THEN
          STIMES(6) = STIMES(6) + ATIME(21)*IQSZ
          IF( IPRECT.GT.1 ) THEN
              STIMES(6) = STIMES(6) + ATIME(20)*IQSZ
          ENDIF
      ENDIF
C
C    Stage 7 - CPHF solution on a separate phase
C
      IF( DORESP .AND. IDENS.GE.3 .AND. IDENS.LT.8 ) THEN
          STIMES(7) = STIMES(7) + 
     .                    PSDEIS(TK,TM2A,TFOCK,TA2M,IKSIZE,IKSTRP)
      ENDIF
C
C    Stage 8 - Half-electron RHS M.O.->A.O.
C
      IF( DORESP .AND. (HALFEL.OR.DOCI) ) THEN
          STIMES(8) = STIMES(8) + TM2A + NMRACT*TDSYR2
      ENDIF
C
C    Stage 9 - summation of the response contributions
C
      IF( DORESP ) THEN
          IF( IMIX.EQ.2 ) THEN
              STIMES(9) = STIMES(9) + ATIME(19)*NVARS*DNROS
          ENDIF
          IF( IDENS.EQ.8 ) THEN
              IF( MODE.GE.2 ) THEN
                  IF( IMIX.EQ.3 ) THEN
                      STIMES(9) = STIMES(9) + TMIX
                  ENDIF
                  STIMES(9) = STIMES(9) + DNVARS*TA2M
                  IF( IQSWAP.GE.0 ) THEN
                      STIMES(9) = STIMES(9) + ATIME(19)*DNVARS*DIQSZ
                  ENDIF
              ENDIF
              STIMES(9) = STIMES(9) + 
     .                        PSDEIS(TK,TM2A,TFOCK,TA2M,IKSIZE,IKSTRP)
          ENDIF
          IF( HALFEL.OR.DOCI ) THEN
              STIMES(9) = STIMES(9) + ATIME(26)*2*DNROS*NVARS
              IF( IMIX.EQ.3 ) THEN
                  STIMES(9) = STIMES(9) + TMIX
              ENDIF
          ENDIF
          IF( MODE.GE.2 ) THEN
              IF( IDENS.EQ.1 ) THEN
C
C                Warning, DECONV is < 1.0 !
C
                  IF(UHF) THEN
                      STIMES(9) = STIMES(9) - TWO*
     .                        ATIME(29)*DNVARS*DNORBS**3*LOG10(DECONV)
                  ELSE
                      STIMES(9) = STIMES(9) -
     .                        ATIME(29)*DNVARS*DNORBS**3*LOG10(DECONV)
                  ENDIF
              ELSE
                  STIMES(9) = STIMES(9) + NCPVRS*TM2A
              ENDIF
              STIMES(9) = STIMES(9) + ATIME(26)*2*DNROS*NVARS**2
              IF( IMIX.EQ.3 ) THEN
                  STIMES(9) = STIMES(9) + 
     .                        TMIX * DBLE((NVARS+IDSTRP-1)/IDSTRP)
                  IF( HALFEL.OR.DOCI ) THEN
                      STIMES(9) = STIMES(9) - TMIX
                  ENDIF
              ENDIF
          ENDIF
      ENDIF
C
C   Scale execution times by the machine speed and return.
C
      XTIME = ZERO
      DO 9000 I=0,9
          STIMES(I) = STIMES(I) * ATIME(1)
          XTIME = XTIME + STIMES(I)
 9000 CONTINUE
      STIMES(-1) = XTIME
C
      RETURN
      END
C
      FUNCTION PSDEKM(TINT1,TINT2,TINT3,TINT4)
C
C   Estimate construction time of the CPHF K matrix using explicit
C   transformation of integrals
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      TINT1  - Time required for the single first-index transform
C      TINT2  - Time required for the single second-index transform
C      TINT3  - Time required for the single third-index transform
C      TINT4  - Time required for the single forth-index transform
C
C   Accessed common blocks:
C
C      PSDGBL - Scratch common for global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C      PSOCC  - Occupation numbers description
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C     Timings are approximate.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXGRP=10)
C
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ONE=1.0D0)
      PARAMETER (TWO=2.0D0)
      PARAMETER (THREE=3.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSDGB2/ NPAIR, NPAIR2, NOCCA, NVACA, NOCCB, NVACB, 
     .         IQSZA, IQSZB,  IQSZ,  ICPV1, ICPVL, NCPVRS
     ./PSOCC / DOCCA(MAXGRP), DOCCB(MAXGRP), NMGRPA, NMGRPB, 
     .         IOPN1A, IOPN1B, IOPNLA, IOPNLB, IMAXBL,
     .         IG1STA(MAXGRP), IGCNTA(MAXGRP), 
     .         IG1STB(MAXGRP), IGCNTB(MAXGRP),
     .         IGBASA(MAXGRP,MAXGRP), IGBASB(MAXGRP,MAXGRP)
      SAVE /PSDGBL/, /PSDGB2/, /PSOCC /
C
      D1ST   = ZERO
      D2ND   = ZERO
      D3RD   = ZERO
      D4TH   = ZERO
C
      DNORBS = DBLE(NORBS)
      DIQSZA = DBLE(IQSZA)
      DIQSZB = DBLE(IQSZB)
C
C     Alpha-Alpha block
C
      DO 100 I=1,NMGRPA-1
          D1ST = D1ST + IGCNTA(I)*(NMGRPA-I)
  100 CONTINUE
C
      IF( NMGRPA.GE.2 ) THEN
          D2ND = D2ND + HALF*IGCNTA(1)*(ONE+IGCNTA(1)+TWO*IGCNTA(2))
          IF( NMGRPA.GE.3 ) THEN
              D2ND = D2ND + IGCNTA(1)*(NMGRPA-2)*(IGCNTA(1)+IGCNTA(2))
              DO 200 I=3,NMGRPA
                  D2ND = D2ND + IGCNTA(1)*(NMGRPA+1-I)*IGCNTA(I)
  200         CONTINUE
              DO 300 I=2,NMGRPA-1
                  D2ND = D2ND + DNORBS*(NMGRPA-I)*IGCNTA(I)
  300         CONTINUE
          ENDIF
      ENDIF
C
      D3RD = D3RD + HALF*DIQSZA*(THREE*DNORBS-
     .                           TWO*IGCNTA(1)-IGCNTA(NMGRPA))
      D4TH = D4TH + THREE*HALF*DIQSZA*(DIQSZA+1) - DIQSZA
C
      IF(UHF) THEN
C
          D3RD = D3RD + DIQSZB*(DNORBS-IGCNTA(NMGRPA))
          D4TH = D4TH + DIQSZB*DIQSZA
C
          DO 1100 I=1,NMGRPB-1
              D1ST = D1ST + IGCNTB(I)*(NMGRPB-I)
 1100     CONTINUE
C
          IF( NMGRPB.GE.2 ) THEN
              D2ND = D2ND + HALF*IGCNTB(1)*(ONE+IGCNTB(1)+TWO*IGCNTB(2))
              IF( NMGRPB.GE.3 ) THEN
                  D2ND = D2ND + IGCNTB(1)*(NMGRPB-2)*
     .                              (IGCNTB(1)+IGCNTB(2))
                  DO 1200 I=3,NMGRPB
                          D2ND = D2ND + IGCNTB(1)*(NMGRPB+1-I)*IGCNTB(I)
 1200             CONTINUE
                  DO 1300 I=2,NMGRPB-1
                      D2ND = D2ND + DNORBS*(NMGRPB-I)*IGCNTB(I)
 1300             CONTINUE
              ENDIF
          ENDIF
C
          D3RD = D3RD + HALF*DIQSZB*(THREE*DNORBS-
     .                           TWO*IGCNTB(1)-IGCNTB(NMGRPB))
          D4TH = D4TH + THREE*HALF*DIQSZB*(DIQSZB+1) - DIQSZB
C
      ENDIF
C
      PSDEKM = TINT1*D1ST + TINT2*D2ND + TINT3*D3RD + TINT4*D4TH
C
      RETURN
      END
C
      FUNCTION PSDEKA(TDSYR2,TFOCKA,TFOCK,TA2MA,A30)
C
C   Estimate construction time of the CPHF K matrix using AO
C   basis approach
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      TDSYR2 - Time for the single invocation of DSYR2 on the
C               NORBSxNORBS matrix
C      TFOCKA - Time for construction of the alpha Fock matrix
C      TFOCK  - Time for construction of the both Fock matrices
C      TA2MA  - Time for complete A.O.->M.O. transformation of
C               alpha part
C      A10    - Time for the single flop in DSYMM and DGEMM
C
C   Accessed common blocks:
C
C      PSDGBL - Scratch common for global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C      PSOCC  - Occupation numbers description
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C     Timings for the partial transforms are approximate.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXGRP=10)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSDGB2/ NPAIR, NPAIR2, NOCCA, NVACA, NOCCB, NVACB, 
     .         IQSZA, IQSZB,  IQSZ,  ICPV1, ICPVL, NCPVRS
     ./PSOCC / DOCCA(MAXGRP), DOCCB(MAXGRP), NMGRPA, NMGRPB, 
     .         IOPN1A, IOPN1B, IOPNLA, IOPNLB, IMAXBL,
     .         IG1STA(MAXGRP), IGCNTA(MAXGRP), 
     .         IG1STB(MAXGRP), IGCNTB(MAXGRP),
     .         IGBASA(MAXGRP,MAXGRP), IGBASB(MAXGRP,MAXGRP)
      SAVE /PSDGBL/, /PSDGB2/, /PSOCC /
C
C    Construction of RHS
C
      PSDEKA = TDSYR2*IQSZ
C
C    "Fock" matrices formation
C
      PSDEKA = PSDEKA + IQSZA * TFOCKA
      IF(UHF) THEN
          PSDEKA = PSDEKA + IQSZB * TFOCK
      ENDIF
C
C    Complete A.O.->M.O. transforms
C
      IF(UHF) THEN
          PSDEKA = PSDEKA + IQSZB * TA2MA
      ENDIF
C
C    Partial A.O.->M.O. transforms
C
      PSDEKA = PSDEKA + A30*PSDEKF(NORBS,IQSZA,NMGRPA,IG1STA,IGCNTA)
      IF(UHF) THEN
          PSDEKA = PSDEKA + A30*PSDEKF(NORBS,IQSZB,NMGRPB,IG1STB,IGCNTB)
      ENDIF
      RETURN
      END
C
      FUNCTION PSDEKF(NORBS,IQSZ,NMGRP,IG1ST,IGCNT)
C
C   Estimate number of FLOPS required to perform partial transformations
C   of "Fock" matrices required in the construction of the triangular
C   fragment of the CPHF K matrix for the single spin.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NORBS  - Number of M.O.s
C      IQSZ   - Size of a spin-specific part of the solution vector
C      NMGRP  - Number of occupation number groups
C      IG1ST  - Indices of the first orbital in the group.
C      IGCNT  - Number of orbitals in the group
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
      PARAMETER (HALF=0.5D0)
      PARAMETER (TWO =2.0D0)
      DIMENSION IGCNT(NMGRP), IG1ST(NMGRP)
C
      DNORBS = DBLE(NORBS)
      DIQSZ  = DBLE(IQSZ)
      PSDEKF = DIQSZ + HALF*DIQSZ*(DIQSZ+1)
      DO 900 I=1,NMGRP-1
          DN     = DBLE(IGCNT(I))
          DNSF   = DBLE(IG1ST(I+1)-1)
          PSDEKF = PSDEKF + DNORBS*DN*(DNORBS-DNSF)*DNSF
  900 CONTINUE
      PSDEKF = PSDEKF * TWO*DNORBS
C
      RETURN
      END
C
      FUNCTION PSDEIS(TK,TM2A,TFOCK,TA2M,IKSIZE,IKSTRP)
C
C   Estimate computation time for the solution phase of the CPHF
C   process.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      TK     - Time for the CPHF K matrix computation
C      TM2A   - Time for M.O.->A.O. transform
C      TFOCK  - Time for the Fock matrix construction
C      TA2M   - Time for A.O.->M.O. transform
C      IKSIZE - Size of the CPHF K matrix on disk (valid only
C               then IDENS is 5)
C      IKSTRP - Number of strips in the block-packed format
C               of the CPHF K matrix (valid only if IDENS is
C               5 or 6)
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGBL - Scratch common for global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C      PSDTMC - Scaling constants for estimation of execution
C               time.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      This is essentially a simulation of PSDRV with linear
C      fitted approximations to the excution times. See 
C      function body for the description of the actual terms
C      computed.
C
C      All computations have to be done in DOUBLE PRECISION, since
C      the range of INTEGER is not sufficient to hold counts which
C      are generated during time estimations for large systems with
C      esoteric profiles.
C
C   Bugs:
C
C      Some auxiliary times are computed even if never used due
C      to the combination of options.
C
C      Timing for the CPHF K matrix formation are approximate
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NTIMEC=33)
C
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
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
     ./PSDTMC/ ATIME(NTIMEC)
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSDTMC/
C
      PSDEIS = ZERO
      DIQSZ  = DBLE(IQSZ)
C
      IF( IQSWAP.EQ.1 ) THEN
          PSDEIS = PSDEIS + ATIME(19)*DIQSZ*NCPVRS
      ENDIF
      IF( IDENS.EQ.3 ) THEN
          PSDEIS = PSDEIS + ATIME(22)*DIQSZ**3
     .                    + ATIME(23)*DIQSZ**2*NCPVRS
      ENDIF
      IF( IDENS.NE.3 ) THEN
          NPASS  = (NCPVRS+INRHS-1)/INRHS
          NVECS  = NCPVRS*IAVEIT
          DNCPVR = DBLE(NCPVRS)
          DAVEIT = DBLE(NCPVRS)
          DNPASS = DBLE(NPASS)
          DNVECS = DBLE(NVECS)
          DNRHS  = DBLE(INRHS)
          DKRSAV = DBLE(IKRSAV)
C
C         RHS I/O penalty
C
          IF( IQSWAP.EQ.2 ) THEN
              PSDEIS = PSDEIS + ATIME(19)*2*DIQSZ*DNCPVR
          ENDIF
C
C         CPHF K fetching/recomputation penalty
C
          IF( IDENS.EQ.5 .AND. IROWS.LT.IQSZ ) THEN
              PSDEIS = PSDEIS + ATIME(19)*DBLE(IKSIZE)*DNPASS*DAVEIT
          ENDIF
          IF( IDENS.EQ.6 ) THEN
              PSDEIS = PSDEIS + TK*DNPASS*DAVEIT
          ENDIF
C
C         Computing (Symm)K X products
C
          IF( IDENS.EQ.4 ) THEN
              PSDEIS = PSDEIS + ATIME(24)*DNVECS*2*DIQSZ*(DIQSZ+1)
          ENDIF
          IF( IDENS.EQ.5 .OR. IDENS.EQ.6 ) THEN
              PSDEIS = PSDEIS + ATIME(25)*DNVECS*2*DIQSZ*(DIQSZ+1)
          ENDIF
          IF( IDENS.GE.7 ) THEN
              PSDEIS = PSDEIS + DNVECS*(TM2A + TFOCK + TA2M)
          ENDIF
C
C         Cache memory penalties
C
          IF( IDENS.EQ.5 .OR. IDENS.EQ.6 ) THEN
              PSDEIS = PSDEIS + ATIME(31)*( DNVECS*DIQSZ*DBLE(IKSTRP) + 
     .                          DAVEIT*DNPASS*HALF*DIQSZ*(DIQSZ+1) )
          ENDIF
C
C         Conditioning
C
          ICND = NCPVRS + NVECS
          IF( IPRECT.GT.1 ) THEN
              ICND = ICND + NCPVRS + NVECS
              IF( IDENS.GE.7 ) THEN
                  ICND = ICND + NVECS
              ENDIF
          ENDIF
          IF( HALFEL ) THEN
              ICND = ICND - 1
          ENDIF
          PSDEIS = PSDEIS + ATIME(27)*DIQSZ*ICND
          IF( ISOLVE.EQ.1 ) THEN
C
C             Simple solver, dominated by dot products.
C
              DIDDOT  = DNCPVR*(DAVEIT**2+9*DAVEIT+1)
              DIDAXP  = 0
              DISVD   = 0 
          ENDIF
          IF( ISOLVE.EQ.2 ) THEN
C
C             Minimum residual solver.
C
              DIDDOT  = DNPASS*DNRHS*( (1+DAVEIT)*DBLE( 2 + 5*DNRHS + 
     .                        3*DKRSAV + 3*DAVEIT/2 ) + DKRSAV )
              DIDAXP  = DNPASS*( DNRHS**2*DAVEIT*(2+DAVEIT) + 2*DNRHS*
     .                        (1+DAVEIT) + DKRSAV*(1+DAVEIT+DNRHS) )
              DISVD   = DNPASS*DNRHS**2*(1+DAVEIT)
          ENDIF
          IF( ISOLVE.EQ.3 ) THEN
C
C             Orthogonal residual solver (globally orthogonal CG)
C
              DIDDOT  = DNPASS*DNRHS*( DKRSAV + (1+DAVEIT)*( 1 + 
     .                         4*DNRHS + 2*DKRSAV + DNRHS*DAVEIT ) )
              DIDAXP  = DNPASS*DNRHS*( DKRSAV - DNRHS + (1+DAVEIT)*
     .                         (2 + DNRHS + DKRSAV + (DNRHS*DAVEIT)/2 ))
              DISVD   = DNPASS*DNRHS**2*(1+DAVEIT)
          ENDIF
          IF( ISOLVE.EQ.4 ) THEN
C
C             Minimal residue with locally orthogonal basis
C
              DIDDOT  = DNPASS*DNRHS*(1+DAVEIT)*(2+4*DKRSAV+3*DNRHS)
              DIDAXP  = DNPASS*DNRHS*(1+DAVEIT)*(1+3*DKRSAV+3*DNRHS)
              DISVD   = DNPASS*DNRHS**2*(1+DAVEIT)
          ENDIF
          IF( ISOLVE.EQ.5 ) THEN
C
C             Orthogonal residual with locally orthogonal basis
C
              DIDDOT  = DNPASS*DNRHS*(1+DAVEIT)*(1+3*DKRSAV+2*DNRHS)
              DIDAXP  = DNPASS*DNRHS*(1+DAVEIT)*(1+2*DKRSAV+2*DNRHS)
              DISVD   = DNPASS*DNRHS**2*(1+DAVEIT)
          ENDIF
          IF( ISOLVE.EQ.6 ) THEN
C
C             Generalized CG
C
              DIDDOT  = DNPASS*DNRHS*(1+DAVEIT)*(1+5*DKRSAV+4*DNRHS)
              DIDAXP  = DNPASS*DNRHS*(1+DAVEIT)*(1+3*DKRSAV+3*DNRHS)
              DISVD   = DNPASS*DNRHS**2*(1+DAVEIT)
          ENDIF
          IF( INRHS.EQ.1 ) DISVD = 0
          PSDEIS = PSDEIS + ATIME(28)*DIQSZ*DIDDOT 
          PSDEIS = PSDEIS + ATIME(32)*DIQSZ*DIDAXP 
          PSDEIS = PSDEIS + ATIME(33)*DIQSZ*DISVD
      ENDIF
C
      RETURN
      END
