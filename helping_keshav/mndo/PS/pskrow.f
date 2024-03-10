C     ******************************************************************
C
C     Compute a single row of the CPHF K matrix in the AO basis.
C
C     ******************************************************************
      SUBROUTINE PSKROW(DUMP,CALP,CBET,LDC,INBET,MOI,MOJ,OUTBET,I3,
     .                  MOK,I4,MOL,ROW)
C
C   Compute arbitrary part of the arbitrary vector of the CPHF K
C   matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      DUMP   - Base of the dynamic memory items
C      CALP   - Alpha orbital coefficients
C      CBET   - Beta orbital coefficients (not used if
C               UHF flag not set)
C      LDC    - Leading dimension of CALP and CBET
C      INBET  - .TRUE. if MOI and MOJ correspond to
C               beta M.O.s
C      MOI,
C      MOJ    - Orbitals identifying row of the K matrix
C               to compute, MOI <= MOJ
C      OUTBET - .TRUE. if complete alpha portion of the
C               K matrix row is desired, with MOL and
C               MOK specifying offsets in the beta part.
C      I3     - Occupation block MOK belongs to
C      I4     - Occupation block MOL belongs to, I4 > I3
C      MOK,
C      MOL    - Orbitals identifying ending point of the
C               row to compute, MOL > MOK
C      ROW    - Place for the output row.
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C      PSDGB2 - Response computation parameters
C      PSDYNM - Dynamic memory handles
C      PSOCC  - Occupation number groups
C      PSPRT  - Printing unit
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
C      This is the slightly glorifyed version of PSDS3A,
C      see it for comments. The basic idea is to feed
C      vectors containg zeros in all positions except
C      (MOI,MOJ) to the PSDS3A, which will compute it's
C      product by the CPHF K matrix, giving desired row.
C
C   Bugs:
C
C      Zeroing entire CPAYA and CPAYB is an overkill.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXGRP=10)
C
      PARAMETER (ONE  = 1.0D0)
      PARAMETER (SONE =-1.0D0)
C
      LOGICAL INBET, OUTBET
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSDGB2/ NPAIR, NPAIR2, NOCCA, NVACA, NOCCB, NVACB,
     .         IQSZA, IQSZB,  IQSZ,  ICPV1, ICPVL, NCPVRS
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
     ./PSOCC / DOCCA(MAXGRP), DOCCB(MAXGRP), NMGRPA, NMGRPB,
     .         IOPN1A, IOPN1B, IOPNLA, IOPNLB, IMAXBL,
     .         IG1STA(MAXGRP), IGCNTA(MAXGRP),
     .         IG1STB(MAXGRP), IGCNTB(MAXGRP),
     .         IGBASA(MAXGRP,MAXGRP), IGBASB(MAXGRP,MAXGRP)
     ./PSPRT / NB6
      SAVE /PSDGBL/, /PSDGB2/, /PSOCC/, /PSPRT /
C
      DIMENSION CALP(LDC,*), CBET(LDC,*), DUMP(*), ROW(*)
C
C    Sanity checks
C
      IF( (INBET.OR.OUTBET) .AND. .NOT.UHF ) THEN
          WRITE(NB6,10100)
          STOP 'PSKROW'
      ENDIF
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
      IA2MTM = IPSMOF(LA2MTM)
C
C    Construct AO representation of the special RHS vector.
C    The procedure is equivalent to PSXM2A in this special
C    case, but is much cheaper.
C
      CALL PSZRM(NORBS,NORBS,DUMP(ICPAYA),NORBS)
      IF(UHF) THEN
          CALL PSZRM(NORBS,NORBS,DUMP(ICPAYB),NORBS)
      ENDIF
      IF( .NOT.LIMAG ) THEN
C
C    Real perturbation
C
          IF( .NOT.INBET ) THEN
              CALL DSYR2('U',NORBS,ONE,CALP(1,MOI),1,CALP(1,MOJ),1,
     .               DUMP(ICPAYA),NORBS)
          ELSE
              CALL DSYR2('U',NORBS,ONE,CBET(1,MOI),1,CBET(1,MOJ),1,
     .               DUMP(ICPAYB),NORBS)
          ENDIF
      ELSE
C
C    Imaginary perturbation, only alpha part may be present. Unfortunately,
C    BLAS does not make provisions for handling antisymmetric matrices, so
C    that calls to DGER below make two times more work than is necessary.
C
          IF(INBET) THEN
              WRITE(NB6,10105)
              STOP 'PSKROW'
          ENDIF
C
          CALL DGER(NORBS,NORBS, ONE,CALP(1,MOI),1,CALP(1,MOJ),1,
     .           DUMP(ICPAYA),NORBS)
          CALL DGER(NORBS,NORBS,SONE,CALP(1,MOJ),1,CALP(1,MOI),1,
     .           DUMP(ICPAYA),NORBS)
      ENDIF
C
C     Compute "response" matrices
C
      CALL PSRESP(DUMP(ICPAYA),DUMP(ICPAYB),DUMP(ICPAY),DUMP(IAI2T1),
     .            DUMP(ICPATA),DUMP(ICPATB),NORBS,OUTBET)
C
C    Transform "response" matrices back to the MO basis. This is a 
C    little trickier than usual PSXU2M, since we do not want to
C    compute more elements of the output vector then is necessary.
C
      IF(OUTBET) THEN
C
C         Alpha part should be computed entirely, with fractional
C         beta part.
C
          IF( IQSZA.NE.0 ) THEN
              CALL PSXAM1(NMGRPA,IG1STA,IGCNTA,IGBASA,MAXGRP,CALP,LDC,
     .                    DUMP(ICPATA),NORBS,ROW,DUMP(IA2MTM),NORBS)
          ENDIF
          IF( IQSZB.NE.0 ) THEN
              IF( MOK.EQ.IG1STB(NMGRPB-1)+IGCNTB(NMGRPB-1)-1 .AND.
     .            MOL.EQ.IG1STB(NMGRPB  )+IGCNTB(NMGRPB  )-1 ) THEN
C
C                 Beta part happen to be complete, so more efficient
C                 treatment is available.
C
                  CALL PSXAM1(NMGRPB,IG1STB,IGCNTB,IGBASB,MAXGRP,CBET,
     .                        LDC,DUMP(ICPATB),NORBS,ROW(1+IQSZA),
     .                        DUMP(IA2MTM),NORBS)
              ELSE
                  CALL PSXAMP(NMGRPB,IG1STB,IGCNTB,IGBASB,MAXGRP,CBET,
     .                        LDC,DUMP(ICPATB),NORBS,ROW(1+IQSZA),
     .                        DUMP(IA2MTM),NORBS,I3,MOK,I4,MOL)
              ENDIF
          ENDIF
      ELSE
C
C        Partial alpha part is desired
C
          IF( IQSZA.NE.0 ) THEN
              IF( MOK.EQ.IG1STA(NMGRPA-1)+IGCNTA(NMGRPA-1)-1 .AND.
     .            MOL.EQ.IG1STA(NMGRPA  )+IGCNTA(NMGRPA  )-1 ) THEN
C
C                 Alpha part happen to be complete, so more efficient
C                 treatment is available.
C
                  CALL PSXAM1(NMGRPA,IG1STA,IGCNTA,IGBASA,MAXGRP,
     .                        CALP,LDC,DUMP(ICPATA),NORBS,ROW,
     .                        DUMP(IA2MTM),NORBS)
              ELSE
                  CALL PSXAMP(NMGRPA,IG1STA,IGCNTA,IGBASA,MAXGRP,
     .                        CALP,LDC,DUMP(ICPATA),NORBS,ROW,
     .                        DUMP(IA2MTM),NORBS,I3,MOK,I4,MOL)
              ENDIF
          ENDIF
      ENDIF
C
      RETURN
C
10100 FORMAT(' BETA HALF OF THE RHF VECTOR REQUESTED FROM PSKROW...' )
10105 FORMAT(' UHF IS NOT SUPPORTED IN PSKROW FOR IMAGINARY',
     .       ' PERTURBATIONS')
C
11000 FORMAT(' KRYLOV VECTOR IS:'/)
11010 FORMAT(' YA IS:'/)
11015 FORMAT(' YB IS:'/)
11020 FORMAT(' Y IS:'/)
11100 FORMAT(' ALPHA "RESPONSE" MATRIX IS:'/)
11105 FORMAT(' BETA "RESPONSE" MATRIX IS:'/)
11200 FORMAT(' NEW KRYLOV VECTOR IS:'/)
      END
