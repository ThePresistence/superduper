C     ******************************************************************
C
C     Explicit formation of linear equations for CPHF.
C     Computation of the CPHF K matrix. RHF and UHF cases.
C
C     ******************************************************************
C
      SUBROUTINE PSDS2P(CALP,CBET,LDC,EALP,EBET,DUMP)
C
C   Originally, PSDS2D was implemented as an alternative entry point
C   to subroutine PSDS2P. This caused problems with the Intel Fortran
C   compiler probably because LDA appears in the dimension statement
C   of AIN and OUT but is undefined when calling PSDS2P. Making PSDS2P
C   and PSDS2D separate subroutines, PSDS2P passing dummy arguments to
C   PSDS2D, fixes this problem.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CALP(LDC,*), CBET(LDC,*), EALP(*), EBET(*), DUMP(*)
      DIMENSION AIN(1,1), OUT(1,1)
      IVARC = 0
      LDA   = 1
      CALL PSDS2D(CALP,CBET,LDC,EALP,EBET,DUMP,IVARC,AIN,OUT,LDA)
      RETURN
      END
C
      SUBROUTINE PSDS2D(CALP,CBET,LDC,EALP,EBET,DUMP,IVARC,AIN,OUT,LDA)
C
C   Compute symmetric packed version of CPHF K matrix.
C   Resulting matrix could either be stored in memory
C   (IDENS=3,4,5), written to the file IUK (IDENS=5) or
C   used on-the fly by iterative solver (IDENS=6)
C
C   Calls from iterative solver should use alternative
C   entry point PSDS2D, which accepts additional parameters
C   needed by matrix product update routine.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients
C      CBET   - Beta orbital coefficients, not used if UHF flag
C               not set.
C      LDC    - Leading dimension of CALP and CBET arrays.
C      EALP   - Alpha orbital energies
C      EBET   - Beta orbital energies (not used if UHF flag
C               not set.
C      DUMP   - Base of dynamic memory items
C
C      IVARC  - Number of variables processed in parallel. Used
C               by direct entry point PSDS2D only
C      AIN    - Input vector(s). Used by direct entry point only.
C      OUT    - Output vector(s). Used by direct entry point only.
C      LDA    - Leading dimension of AIN and OUT. Used by direct entry
C               point only.
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGBL - Global computation parameters
C      PSDGB2 - Parameters of "response" computation
C      PSDYNM - Dynamic memory handles
C      PSPRTF - Debug output tuning flags
C      PSOCC  - Occupation groups description
C
C   Modified common blocks:
C
C      PS2XCH - Variables shared between PSDS2P and it's
C               workhorses. Here is brief description of
C               variables in PS2XCH:
C          DIAGWT  = 2.0 for UHF, 4.0 for RHF
C          SCALLJ  = 0.5 for UHF, 0.75 for RHF
C          WIKJL   = -1 always.
C          WILJK   = -1 for real perturbations, 1 for imaginary ones.
C          MOI,MOJ - Compound row index of the CPHF K matrix
C          I1      - Occupation group number of the MOI orbital
C          IND1    - Index of the first orbital in the I1 group
C          LAST1   - Index of the last orbital in the I1 group
C          N1      - Number of orbitals in the I1 group
C          O1      - Occupation number of the I1 group
C          I2, IND2, LAST2, N2, O2
C                  - Same for the occupation group of MOJ orbital
C          ISTART  - Number of first row in the KATMP buffer,
C                    valid if IDENS is 5 or 6
C          IEND    - Number of the last row in the KATMP buffer,
C          ICURR   - Current row in the KATMP buffer
C          IOLDST, IOLDEN
C                  - Values of ISTART and IEND prior to the call to
C                    PSDS2L, valid if IDENS is 6
C          IROW    - Memory arena offset of the current subrow
C          ISHIFT  - Position of the shift vector for the current
C                    row (memory arena index), valid if USESHF is
C                    .TRUE.
C          IKATMP  - Base offset of the KATMP buffer, valid if
C                    IDENS is 5 or 6
C          SETSHF  - .TRUE. if the shift vector should be determined
C          USESHF  - .TRUE. if the shift vector should be applied
C                    (always true if SETSHF is set)
C          FLUSH   - .TRUE. if KATMP wrapped after the last row.
C                    Valid if IDENS is 6.
C
C   Local storage:
C
C   Module logic:
C
C      Then IKMODE is 1, explicit transformation of integrals
C      is used. Then IKMODE is 2, Fock-type procedure is used.
C
C      This code takes advantage of the fact that second-inder
C      transformed ERIs in NDDO are only N^2 in size if we fix
C      the first index. Therefore, by precomputing and storing
C      all needed two-index transformed ERIs, we could compute
C      CPHF K matrix without either redundant work or sorting
C      pass. See program text for more detailed comments.
C
C      Subroutines PSDS2I, PSDS2J, PSDS2L and PSDS2K are very
C      tightly coupled with PSDS2P/PSDS2D, and could have been
C      represented much better by internal subroutines.
C
C   Bugs:
C
C      First-, and second- index transforms are performed
C      more then once for the same values of indices if number
C      of occupation groups is more than 2.
C
C      This subroutine depends on the exact ordering of the
C      occupation group blocks in the solution vector.
C
C      Passing around common variables in PS2XCH is ugly, and
C      can inhibit some optimizations. Inlining all PSDS2* routines
C      called from PSDS2P is probably a way to go.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXGRP=10)
C
      PARAMETER (ONE = 1.0D0)
      PARAMETER (SONE=-1.0D0)
      PARAMETER (TWO = 2.0D0)
      PARAMETER (FOUR= 4.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      LOGICAL SETSHF, USESHF, FLUSH
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
     ./PS2XCH/ DIAGWT, SCALLJ, O1, O2, WIKJL, WILJK,
     .         I1, IND1, N1, LAST1, MOI, I2, IND2, LAST2, N2, MOJ, 
     .         ISTART, IEND, ICURR, IOLDST, IOLDEN,
     .         IROW, ISHIFT, IKATMP, SETSHF, USESHF, FLUSH
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSPRTF/,/PSOCC/, /PSPRT /
C
      DIMENSION CALP(LDC,*), CBET(LDC,*), EALP(*), EBET(*), DUMP(*)
      DIMENSION AIN(LDA,*), OUT(LDA,*)
C
      FLUSH = .FALSE.
      WIKJL = SONE
      WILJK = SONE
      IF(LIMAG) WILJK = ONE
C
      IF( IKMODE.EQ.1 ) THEN
          IINTS = IPSMOF(LAI2T5)
          IF(UHF) THEN
              DIAGWT = TWO
          ELSE
              DIAGWT = FOUR
          ENDIF
          SCALLJ = (DIAGWT-ONE)/DIAGWT
      ENDIF
      IF( IDENS.EQ.3 .OR. IDENS.EQ.4 ) THEN
          IKS    = IPSMOF(LKS)
          IROW   = IKS
      ELSE IF( IDENS.EQ.5 .OR. IDENS.EQ.6 ) THEN
          IKATMP = IPSMOF(LKATMP)
          IROW   = IKATMP
          ISTART = 1
          IEND   = IPSBPE(ISTART)
          ICURR  = ISTART
      ELSE
          STOP 'PSDS2P'
      ENDIF
      SETSHF = .FALSE.
      USESHF = .FALSE.
      IF( IPRECT.GT.1 ) THEN
          ISHIFT = IPSMOF(LSHIFT)
          USESHF = .TRUE.
          IF( IPRECT.EQ.4 ) THEN
              SETSHF = .TRUE.
          ENDIF
      ENDIF
C
C     First quarter of K matrix is the same for UHF and RHF
C     modulo DIAGWT parameter.
C
      DO 2000 I1=1,NMGRPA-1
          IND1  = IG1STA(I1)
          N1    = IGCNTA(I1)
          O1    = DOCCA(I1)
          LAST1 = IND1 + N1 - 1
          DO 1900 I2=I1+1,NMGRPA
              IND2  = IG1STA(I2)
              N2    = IGCNTA(I2)
              O2    = DOCCA(I2)
              LAST2 = IND2 + N2 - 1
              DO 1800 MOI=IND1,LAST1
                  IF( IKMODE.EQ.1 ) THEN
C
C                     Perform all necessary first and second transforms
C
                      CALL PST1ST(DUMP,CALP(1,MOI))
                      CALL PSDS2I(DUMP,CALP,LDC,NMGRPA,IG1STA,IGCNTA)
                  ENDIF
C
                  DO 1700 MOJ=IND2,LAST2
C
C                     Compute single row of the K matrix
C
                      IROWSZ = IGBASA(I1,I2) + (MOI-IND1)*N2 + MOJ-IND2
                      IF( IKMODE.EQ.1 ) THEN
                          CALL PSDS2J(DUMP,CALP,LDC,NMGRPA,IG1STA,
     .                               IGCNTA,IGBASA,MAXGRP,DUMP(IINTS),
     .                               DUMP(IROW))
                      ELSE
                          CALL PSKROW(DUMP,CALP,CBET,LDC,.FALSE.,MOI,
     .                               MOJ,.FALSE.,I1,MOI,I2,MOJ,
     .                               DUMP(IROW))
                      ENDIF
C
C                         Advance to the next row
C
                      CALL PSDS2L(DUMP,EALP,IROWSZ)
                      IF(LPCPHK.AND.(IDENS.EQ.3.OR.IDENS.EQ.4)) THEN
C
C                         Report CPHF K matrix so far.
C
                          WRITE(NB6,11000) MOI, MOJ
                          CALL PSDPPM(IROWSZ,DUMP(IKS))
                      ENDIF
                      IF(FLUSH) THEN
                          CALL PSDS3N(IOLDST,IOLDEN,DUMP(IKATMP),IVARC,
     .                                AIN,OUT,LDA)
                          FLUSH = .FALSE.
                      ENDIF
 1700             CONTINUE
 1800         CONTINUE
 1900     CONTINUE
 2000 CONTINUE
C
      IF(.NOT.UHF) RETURN
      IF(LPCPHK) THEN
          WRITE(NB6,11050)
      ENDIF
C
C     UHF beta-alpha and beta-beta parts of CPHF matrix K
C
      DO 4000 I1=1,NMGRPB-1
          IND1  = IG1STB(I1)
          N1    = IGCNTB(I1)
          O1    = DOCCB(I1)
          LAST1 = IND1 + N1 - 1
          DO 3900 I2=I1+1,NMGRPB
              IND2  = IG1STB(I2)
              N2    = IGCNTB(I2)
              O2    = DOCCB(I2)
              LAST2 = IND2 + N2 - 1
              DO 3800 MOI=IND1,LAST1
                  IF( IKMODE.EQ.1 ) THEN
C
C                     Perform all necessary first and second transforms
C
                      CALL PST1ST(DUMP,CBET(1,MOI))
                      CALL PSDS2I(DUMP,CBET,LDC,NMGRPB,IG1STB,IGCNTB)
                  ENDIF
C
                  DO 3700 MOJ=IND2,LAST2
                      IROWSZ = IGBASB(I1,I2) + (MOI-IND1)*N2 + MOJ-IND2
                      IF( IKMODE.EQ.1 ) THEN
C
C                         beta-alpha part
C
                          CALL PSDS2K(DUMP,CALP,LDC,NMGRPA,IG1STA,
     .                      IGCNTA,IGBASA,MAXGRP,DUMP(IINTS),DUMP(IROW))
                          IROW = IROW + IQSZA
C
C                         beta-beta part
C
                          CALL PSDS2J(DUMP,CBET,LDC,NMGRPB,IG1STB,
     .                      IGCNTB,IGBASB,MAXGRP,DUMP(IINTS),DUMP(IROW))
                      ELSE
                          CALL PSKROW(DUMP,CALP,CBET,LDC,.TRUE.,MOI,
     .                              MOJ,.TRUE.,I1,MOI,I2,MOJ,DUMP(IROW))
                          IROW = IROW + IQSZA
                      ENDIF
C
C                         Advance to the next row
C
                      CALL PSDS2L(DUMP,EBET,IROWSZ)
                      IF(LPCPHK.AND.(IDENS.EQ.3.OR.IDENS.EQ.4)) THEN
C
C                         Report CPHF K matrix so far.
C
                          WRITE(NB6,11100) MOI, MOJ
                          CALL PSDPPM(IQSZA+IROWSZ,DUMP(IKS))
                      ENDIF
                      IF(FLUSH) THEN
                          CALL PSDS3N(IOLDST,IOLDEN,DUMP(IKATMP),IVARC,
     .                                AIN,OUT,LDA)
                          FLUSH = .FALSE.
                      ENDIF
 3700             CONTINUE
 3800         CONTINUE
 3900     CONTINUE
 4000 CONTINUE
C
      RETURN
11000 FORMAT(' I = ', I7, ' J = ', I7, ' (ALPHA). CPHF K MATRIX SO',
     .       ' FAR:'/)
11050 FORMAT(' STARTING BETA PART OF CPHF MATRIX K ')
11100 FORMAT(' I = ', I7, ' J = ', I7, ' (BETA). CPHF K MATRIX SO',
     .       ' FAR:'/)
      END
C
      SUBROUTINE PSDS2I(DUMP,C,LDC,NMGRP,IG1ST,IGCNT)
C
C   Precompute second-index transformed integrals for the formation
C   of K matrix. This routine should only be called by PSDS2P/PSDS2D.
C   
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      DUMP   - Base of dynamic memory items
C      C      - Orbital coefficients
C      LDC    - Leading dimension of C matrix
C      NMGRP  - Number of orbitals occupation blocks
C      IG1ST  - Base orbitals of the occuptation block
C      IGCNT  - Number of orbitals in the occupation block
C
C   Accessed common blocks:
C
C      PS2XCH - Variables shared between PSDS2P and it's
C               workhorses.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      a. Computation of (IJ,KL) integrals requires transformation
C         by all J from block (I2)
C      b. Computation of (IK,JL) integrals requires transformation
C         by all K from all blocks below block containing I, plus
C         all K<I from the same block.
C      c. Computation of (IL,JK) integrals requires transformation 
C         by Ls from groups (2) to (I2) if I1 is equal to 1, or
C         from (2) to (NMGRP) if I1 is not equal to 1.
C
C      Integrals transformed by the second index J are placed to the
C      slot J to simplify retrieving.
C
C   Bugs:
C
C      At least one occupation group is to be present.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL SETSHF, USESHF, FLUSH
      COMMON 
     ./PS2XCH/ DIAGWT, SCALLJ, O1, O2, WIKJL, WILJK,
     .         I1, IND1, N1, LAST1, MOI, I2, IND2, LAST2, N2, MOJ, 
     .         ISTART, IEND, ICURR, IOLDST, IOLDEN,
     .         IROW, ISHIFT, IKATMP, SETSHF, USESHF, FLUSH
C
      DIMENSION DUMP(*), C(LDC,*), IG1ST(NMGRP), IGCNT(NMGRP)
C
      IFIRST = IG1ST(1)
      ILAST  = IFIRST + IGCNT(1) - 1
      IF( I2.EQ.2 ) ILAST = MOI
      DO 100 J=IFIRST,ILAST
          CALL PST2ND(DUMP,C(1,J),J)
  100 CONTINUE
C
      IGLAST = NMGRP
      IF( I1.EQ.1 ) IGLAST = I2
      DO 300 IX=2,IGLAST
          IFIRST = IG1ST(IX)
          ILAST  = IFIRST + IGCNT(IX) - 1
          DO 200 J=IFIRST,ILAST
              CALL PST2ND(DUMP,C(1,J),J)
  200     CONTINUE
  300 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSDS2J(DUMP,C,LDC,NMGRP,IG1ST,IGCNT,IGBAS,LDI,
     .                  AINTS,ROW)
C
C   Compute alpha-alpha (or beta-beta) part of the single row of
C   the CPHF K matrix. This routine should only be called by 
C   PSDS2P/PSDS2D.
C   
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      DUMP   - Base of dynamic memory items
C      C      - Orbital coefficients
C      LDC    - Leading dimension of C matrix
C      NMGRP  - Number of orbitals occupation blocks
C      IG1ST  - Base orbitals of the occuptation block
C      IGCNT  - Number of orbitals in the occupation block
C      IGBAS  - Base indices of the orbital pairs blocks
C      LDI    - Leading dimension of IGBAS
C      AINTS  - Hiding place of the fully-transformed integrals
C      ROW    - Row of the CPHF K matrix to fill
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C      PS2XCH - Variables shared between PSDS2P and it's
C               workhorses.
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C      Loop controls are overly complicated by the need to support
C      multiply occupation groups.
C
C      IROWSZ should probably be set elswhere, but here it is so easy...
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL SETSHF, USESHF, FLUSH
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PS2XCH/ DIAGWT, SCALLJ, O1, O2, WIKJL, WILJK,
     .         I1, IND1, N1, LAST1, MOI, I2, IND2, LAST2, N2, MOJ, 
     .         ISTART, IEND, ICURR, IOLDST, IOLDEN,
     .         IROW, ISHIFT, IKATMP, SETSHF, USESHF, FLUSH
      SAVE /PSDGBL/
C
      DIMENSION DUMP(*), C(LDC,*), IG1ST(NMGRP), IGCNT(NMGRP)
      DIMENSION IGBAS(LDI,*), AINTS(*), ROW(*)
C
      IF(LIMAG) THEN
C
C     Zero out the new row, so that we can reuse code later on. This
C     is slightly less efficient than just replacing DAXPY's with
C     DCOPY's later on, but what the heck!
C
          N3LAST = MOI-IND1+1
          N4LAST = MOJ-IND2+1
          IROWSZ = IGBAS(I1,I2) + (MOI-IND1)*IGCNT(I2) + (MOJ-IND2)
          CALL PSZRV(IROWSZ,ROW)
      ELSE
C
C     Compute (I,J,K,L) integrals. Orbital blocks with (I3,I4)<(I1,I2)
C     are all rectangular and are computed with a single call. Orbital block
C     (I3,I4) = (I1,I2) have to be split into largest rectangular
C     subblock and residual vector. Additionally, elements with L.EQ.J
C     are to be rescaled, so that we can avoid computing corresponding
C     (I,L,J,K) integrals, which are exactly the same in that case
C
C     Since all (IJKL) integrals cancel for imaginary perturbation,
C     we do not need to compute 'em in the later case.
C
          N3LAST = MOI-IND1+1
          N4LAST = MOJ-IND2+1
          DO 1000 I3=1,I1
              IND3   = IG1ST(I3)
              N3     = IGCNT(I3)
              LAST3  = IND3 + N3 - 1
              I4LAST = NMGRP
              IF( I3.EQ.I1 ) I4LAST = I2 - 1
              IF( I3.LT.I4LAST ) THEN
                  CALL PST3RD(DUMP,MOJ,C(1,IND3),LDC,N3)
                  DO 900 I4=I3+1,I4LAST
                      IND4   = IG1ST(I4)
                      N4     = IGCNT(I4)
                      LAST4  = IND4 + N4 - 1
                      IBLOCK = IGBAS(I3,I4)
                      CALL PST4TH(DUMP,C(1,IND4),LDC,N3,N4)
                      CALL DCOPY(N3*N4,AINTS,1,ROW(IBLOCK),1)
                      IF( I4.EQ.I2 ) THEN
                          CALL DSCAL(N3,SCALLJ,ROW(IBLOCK+N4LAST-1),N4)
                      ENDIF
  900             CONTINUE
              ENDIF
 1000     CONTINUE
C
C     The last block is clipped by the diagonal of K matrix.
C     If third transform was already computed for all orbitals,
C     we'll use that, otherwise we'll do separate transforms 
C     for rectangular block and residual vector
C
          IND4   = IG1ST(I2)
          N4     = IGCNT(I2)
          LAST4  = IND4 + N4 - 1
          N3     = N3LAST
          IROWSZ = (N3-1)*N4 + N4LAST
          IBLOCK = IGBAS(I1,I2)
          IF( I1.LT.I4LAST ) THEN
              CALL PST4TH(DUMP,C(1,IND4),LDC,N3,N4)
              CALL DCOPY(IROWSZ,AINTS,1,ROW(IBLOCK),1)
          ELSE
              N3FULL = N3 - 1 + N4LAST/N4
              IF( N3FULL.GE.1 ) THEN
                   CALL PST3RD(DUMP,MOJ,C(1,IND3),LDC,N3FULL)
                   CALL PST4TH(DUMP,C(1,IND4),LDC,N3FULL,N4)
                   CALL DCOPY(N3FULL*N4,AINTS,1,ROW(IBLOCK),1)
              ENDIF
              IF( N3FULL.NE.N3 ) THEN
                   CALL PST3RD(DUMP,MOJ,C(1,MOI),LDC,1)
                   CALL PST4TH(DUMP,C(1,IND4),LDC,1,N4LAST)
                   CALL DCOPY(N4LAST,AINTS,1,ROW(IBLOCK+N3FULL*N4),1)
              ENDIF
          ENDIF
          CALL DSCAL(N3,SCALLJ,ROW(IBLOCK+N4LAST-1),N4)
          IROWSZ = IBLOCK + IROWSZ - 1
          CALL DSCAL(IROWSZ,DIAGWT,ROW,1)
      ENDIF
C
C     (I,K,J,L) integrals. Processing order is a bit unusual to
C     save on third transforms.
C
C     Transform full rows of blocks
C
      DO 1300 I3=1,I1-1
          IND3   = IG1ST(I3)
          N3     = IGCNT(I3)
          LAST3  = IND3 + N3 - 1
          DO 1200 MOK=IND3,LAST3
              CALL PST3RD(DUMP,MOK,C(1,MOJ),LDC,1)
              DO 1100 I4=I3+1,NMGRP
                  IND4   = IG1ST(I4)
                  N4     = IGCNT(I4)
                  IBLOCK = IGBAS(I3,I4) + N4*(MOK-IND3)
                  CALL PST4TH(DUMP,C(1,IND4),LDC,1,N4)
                  CALL DAXPY(N4,WIKJL,AINTS,1,ROW(IBLOCK),1)
 1100         CONTINUE
 1200     CONTINUE
 1300 CONTINUE
C
C     Block row with I3 = I1. All transforms with MOK<MOI
C     and I4<=I2 are necessary
C
      I3=I1
      IND3   = IG1ST(I3)
      N3     = IGCNT(I3)
      LAST3  = IND3 + N3 - 1
      IF( I3+1 .LE. I2 ) THEN
          DO 1500 MOK=IND3,MOI-1
              CALL PST3RD(DUMP,MOK,C(1,MOJ),LDC,1)
              DO 1400 I4=I3+1,I2
                  IND4   = IG1ST(I4)
                  N4     = IGCNT(I4)
                  IBLOCK = IGBAS(I3,I4) + N4*(MOK-IND3)
                  CALL PST4TH(DUMP,C(1,IND4),LDC,1,N4)
                  CALL DAXPY(N4,WIKJL,AINTS,1,ROW(IBLOCK),1)
 1400         CONTINUE
 1500     CONTINUE
      ENDIF
C
C     MOK=MOI, rows with I4<I2 will enter completely
C
      MOK = MOI
      CALL PST3RD(DUMP,MOK,C(1,MOJ),LDC,1)
      DO 1600 I4=I3+1,I2-1
          IND4   = IG1ST(I4)
          N4     = IGCNT(I4)
          IBLOCK = IGBAS(I3,I4) + N4*(MOK-IND3)
          CALL PST4TH(DUMP,C(1,IND4),LDC,1,N4)
          CALL DAXPY(N4,WIKJL,AINTS,1,ROW(IBLOCK),1)
 1600 CONTINUE
C
C     MOK=MOI, residual row for I4=I2 will enter
C
      I4     = I2
      IND4   = IG1ST(I4)
      N4     = IGCNT(I4)
      IBLOCK = IGBAS(I3,I4) + N4*(MOK-IND3)
      CALL PST4TH(DUMP,C(1,IND4),LDC,1,N4LAST)
      CALL DAXPY(N4LAST,WIKJL,AINTS,1,ROW(IBLOCK),1)
C
C     MOK>MOI, only rows for I4<I2 will enter
C
      IF( I3+1 .LE. I2-1 ) THEN
          DO 1800 MOK=MOI+1,LAST3
              CALL PST3RD(DUMP,MOK,C(1,MOJ),LDC,1)
              DO 1700 I4=I3+1,I2-1
                  IND4   = IG1ST(I4)
                  N4     = IGCNT(I4)
                  IBLOCK = IGBAS(I3,I4) + N4*(MOK-IND3)
                  CALL PST4TH(DUMP,C(1,IND4),LDC,1,N4)
                  CALL DAXPY(N4,WIKJL,AINTS,1,ROW(IBLOCK),1)
 1700         CONTINUE
 1800     CONTINUE
      ENDIF
C
C     (I,L,J,K) integrals. Only ones with L.NE.J are to
C     be computed. It saves us computation of NVACA*NOCCA*(NOCCA+1)/2
C     integrals, which is absolutely insignificant for large systems...
C
C     For imaginary perturbations, we however have to compute all of 'em,
C     since we had dropped (IJKL) integrals altogether.
C
      I4LAST = NMGRP
      IF( I1.EQ.1 ) I4LAST = I2
      DO 3000 I4=2,I4LAST
          IND4  = IG1ST(I4)
          N4    = IGCNT(I4)
          LAST4 = IND4 + N4 - 1
C
C         Warning: proper check is I1.EQ.1 .AND. MOI.EQ.IG1ST(1),
C         but it happens to be equivalent to this one.
C
          IF( MOI.EQ.1 .AND. I4.EQ.I2 ) THEN
              LAST4 = IND4 + N4LAST - 1
          ENDIF
          DO 2900 MOL=IND4,LAST4
              IF( MOL.NE.MOJ .OR. LIMAG ) THEN
                  CALL PST3RD(DUMP,MOL,C(1,MOJ),LDC,1)
                  I3LAST = MIN(I4-1,I1)
                  DO 2800 I3=1,I3LAST
                      IND3   = IG1ST(I3)
                      N3     = IGCNT(I3)
                      IBLOCK = IGBAS(I3,I4) + (MOL-IND4)
                      IF( I3.EQ.I1 .AND. I4.EQ.I2 ) THEN
                          IF( MOL.LE.MOJ ) THEN
                              N3 = N3LAST
                          ELSE
                              N3 = N3LAST - 1
                          ENDIF
                      ENDIF
                      LAST3 = IND3 + N3 - 1
                      CALL PST4TH(DUMP,C(1,IND3),LDC,1,N3)
                      CALL DAXPY(N3,WILJK,AINTS,1,ROW(IBLOCK),N4)
 2800             CONTINUE
              ENDIF
 2900     CONTINUE
 3000 CONTINUE
C
C             We have completed with this row.
C
      RETURN
      END
C
      SUBROUTINE PSDS2K(DUMP,C,LDC,NMGRP,IG1ST,IGCNT,IGBAS,LDI,
     .                  AINTS,ROW)
C
C   Compute beta-alpha part of the single row of the CPHF K matrix.
C   This routine should only be called by PSDS2P/PSDS2D.
C   
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      DUMP   - Base of dynamic memory items
C      C      - Orbital coefficients
C      LDC    - Leading dimension of C matrix
C      NMGRP  - Number of orbitals occupation blocks
C      IG1ST  - Base orbitals of the occuptation block
C      IGCNT  - Number of orbitals in the occupation block
C      IGBAS  - Base indices of the orbital pairs blocks
C      LDI    - Leading dimension of IGBAS
C      AINTS  - Hiding place of the fully-transformed integrals
C      ROW    - Row of the CPHF K matrix to fill
C
C   Accessed common blocks:
C
C      PS2XCH - Variables shared between PSDS2P and it's
C               workhorses.
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
      LOGICAL SETSHF, USESHF, FLUSH
      COMMON 
     ./PS2XCH/ DIAGWT, SCALLJ, O1, O2, WIKJL, WILJK,
     .         I1, IND1, N1, LAST1, MOI, I2, IND2, LAST2, N2, MOJ, 
     .         ISTART, IEND, ICURR, IOLDST, IOLDEN,
     .         IROW, ISHIFT, IKATMP, SETSHF, USESHF, FLUSH
C
      DIMENSION DUMP(*), C(LDC,*), IG1ST(NMGRP), IGCNT(NMGRP)
      DIMENSION IGBAS(LDI,*), AINTS(*), ROW(*)
C
C     (Ibeta,Jbeta,Kalpha,Lalpha) integrals blocks are always 
C     rectangular, and can therefore be handled efficiently.
C
      DO 1000 I3=1,NMGRP-1
          IND3   = IG1ST(I3)
          N3     = IGCNT(I3)
          CALL PST3RD(DUMP,MOJ,C(1,IND3),LDC,N3)
          DO 900 I4=I3+1,NMGRP
              IND4   = IG1ST(I4)
              N4     = IGCNT(I4)
              IBLOCK = IGBAS(I3,I4)
              CALL PST4TH(DUMP,C(1,IND4),LDC,N3,N4)
              CALL DCOPY(N3*N4,AINTS,1,ROW(IBLOCK),1)
              CALL DSCAL(N3*N4,TWO,ROW(IBLOCK),1)
  900     CONTINUE
 1000 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSDS2L(DUMP,E,IROWSZ)
C
C   Finish processing of the row of the CPHF K matrix.
C   This routine should only be called by PSDS2P/PSDS2D.
C   
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      DUMP   - Base of dynamic memory items
C      E      - Orbital energies
C      IROWSZ - Size of the current row
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options
C      PSDGB2 - Response-related quantities
C
C   Modified common blocks:
C
C      PS2XCH - Variables shared between PSDS2P and it's
C               workhorses.
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (AU  = 3.67511944138184490995957368D-2)
C
      LOGICAL SETSHF, USESHF, FLUSH
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
     ./PS2XCH/ DIAGWT, SCALLJ, O1, O2, WIKJL, WILJK,
     .         I1, IND1, N1, LAST1, MOI, I2, IND2, LAST2, N2, MOJ, 
     .         ISTART, IEND, ICURR, IOLDST, IOLDEN,
     .         IROW, ISHIFT, IKATMP, SETSHF, USESHF, FLUSH
      SAVE /PSDOPT/, /PSDGB2/
C
      DIMENSION DUMP(*), E(*)
C
      IROW = IROW+IROWSZ
      IF(SETSHF) THEN
          DUMP(ISHIFT) = DUMP(IROW-1) * DSHIFT
      ENDIF
      IF(USESHF) THEN
          DUMP(IROW-1) = DUMP(IROW-1) - DUMP(ISHIFT)
          ISHIFT = ISHIFT + 1
      ENDIF
      IF( IDENS.EQ.3 ) THEN
C
C         Direct solver needs complete K matrix, so stuff
C         E(vac) - E(occ) to the diagonal and then we are
C         ready to switch to the next row
C
          GAMMA = (E(MOJ) - E(MOI))/(O1 - O2)
          DUMP(IROW-1) = DUMP(IROW-1) + GAMMA * AU 
      ELSE IF( IDENS.EQ.4 ) THEN

C         Iterative solver would work Ok with this version of
C         K matrix provided we apply proper post-conditioner.

      ELSE IF( IDENS.EQ.5 .OR. IDENS.EQ.6 ) THEN
          IF( ICURR.EQ.IEND ) THEN
              IF( IDENS.EQ.5 ) THEN
C
C                 We do not need to flush K matrix if there is enough
C                 space to hold it entirely in memory.
C
                  IF( IROWS.LT.IQSZ ) THEN
                      ISIZE = IEND*(IEND-ISTART+1)
                      CALL PSDWS(IUK,DUMP(IKATMP),ISIZE)
                  ENDIF
              ELSE
                  IOLDST = ISTART
                  IOLDEN = IEND
                  FLUSH  = .TRUE.
              ENDIF
              ISTART = IEND + 1
              IEND   = IPSBPE(ISTART)
              ICURR  = ISTART
              IROW   = IKATMP
          ELSE
C
C             IROW have to be advanced over "slack" entries in the
C             blocked-packed storage format.
C
              IROW  = IROW + (IEND-ICURR)
              ICURR = ICURR + 1
          ENDIF
      ELSE
          STOP 'PSDS2L'
      ENDIF
      RETURN
      END
C
      SUBROUTINE PSDS2S(CALP,CBET,LDC,EALP,EBET,DUMP)
C
C   Computes shift vector used in iterative CPHF solution.
C   This function might be called for IPRECT=2,3,5,6
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients
C      CBET   - Beta orbital coefficients, not used if UHF flag
C               not set.
C      LDC    - Leading dimension of CALP and CBET arrays.
C      EALP   - Alpha orbital energies
C      EBET   - Beta orbital energies (not used if UHF flag
C               not set.
C      DUMP   - Base of dynamic memory items
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGBL - Global computation parameters
C      PSDGB2 - Parameters of "response" computation
C      PSDYNM - Dynamic memory handles
C      PSOCC  - Occupation number groups
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C      Then there is only one occupied or virtual orbital,
C      some of the integrals will be computed needlessly.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXGRP=10)
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
     ./PSOCC / DOCCA(MAXGRP), DOCCB(MAXGRP), NMGRPA, NMGRPB, 
     .         IOPN1A, IOPN1B, IOPNLA, IOPNLB, IMAXBL,
     .         IG1STA(MAXGRP), IGCNTA(MAXGRP), 
     .         IG1STB(MAXGRP), IGCNTB(MAXGRP),
     .         IGBASA(MAXGRP,MAXGRP), IGBASB(MAXGRP,MAXGRP)
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
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSOCC/, /PSPRT /
C
      DIMENSION CALP(LDC,*), CBET(LDC,*), EALP(*), EBET(*), DUMP(*)
C
      ISHIFT = IPSMOF(LSHIFT)
      IF( IPRECT.EQ.2 ) THEN
C
C         Use constant shift supplied by user
C
          DUMP(ISHIFT) = -DSHIFT
          IF( IPRINT.GE.1 ) THEN
              WRITE(NB6,10200) DUMP(ISHIFT)
          ENDIF
          CALL DCOPY(IQSZ-1,DUMP(ISHIFT),0,DUMP(ISHIFT+1),1)
      ELSE IF( IPRECT.EQ.3 .OR. IPRECT.EQ.5 ) THEN
C
C         Use interpolation expression
C
          CALL PSDS2R(NMGRPA,IG1STA,IGCNTA,IGBASA,MAXGRP,DUMP(ISHIFT),
     .                CALP,LDC,EALP,DUMP)
          IF(UHF) THEN
              CALL PSDS2R(NMGRPB,IG1STB,IGCNTB,IGBASB,MAXGRP,
     .                    DUMP(ISHIFT+IQSZA),CBET,LDC,EBET,DUMP)
          ENDIF
      ELSE IF( IPRECT.EQ.6 ) THEN
C
C         Approximate diagonal values of K by corresponding coulomb
C         integrals. This is usually a very good (though expensive)
C         approximation.
C
          CALL PSDS2X(NMGRPA,IG1STA,IGCNTA,IGBASA,MAXGRP,DUMP(ISHIFT),
     .                CALP,LDC,DUMP)
          IF(UHF) THEN
              CALL PSDS2X(NMGRPB,IG1STB,IGCNTB,IGBASB,MAXGRP,
     .                DUMP(ISHIFT+IQSZA),CBET,LDC,DUMP)
          ENDIF
      ELSE
          WRITE(NB6,11000) IPRECT
          STOP 'PSDS2S'
      ENDIF
C
C     Make sure that selected shift values won't bring us too 
C     close to singularity.
C
      CALL PSDS2U(NMGRPA,IG1STA,IGCNTA,IGBASA,MAXGRP,DOCCA,EALP,
     .            DUMP(ISHIFT))
      IF(UHF) THEN
          CALL PSDS2U(NMGRPB,IG1STB,IGCNTB,IGBASB,MAXGRP,DOCCB,EBET,
     .                DUMP(ISHIFT+IQSZA))
      ENDIF
C
      RETURN
11000 FORMAT(' PSDS2S CALLED WITH ILLEGAL VALUE OF IPRECT (', I7, ')'
     .      /' PROGRAM WILL STOP.' )
10200 FORMAT(' USING CONSTANT SHIFT VALUE OF ', G20.14, ' A.U.' )
      END
C
      SUBROUTINE PSDS2R(NMGRP,IG1ST,IGCNT,IGBAS,LDI,SHIFT,C,LDC,E,DUMP)
C
C   Compute single spin part of the CPHF shift vector using 
C   interpolation formula. This function might be called for 
C   IPRECT=3 or 5
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NMGRP  - Number of occupation groups
C      IG1ST  - Indices of the first orbital in a group
C      IGCNT  - Number of orbital in a group
C      IGBAS  - Base indices of orbital groups
C      LDI    - Leading dimension of IGBAS
C      SHIFT  - Computed shift matrix
C      C      - Orbital coefficients
C      LDC    - Leading dimension of C matrix
C      E      - Orbital energies
C      DUMP   - Base of dynamic memory pool
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C   Module logic:
C
C         Shift value for each active coupling block is
C         estimated from K(First1,First2,First1,First2),
C         K(Last1,First2,Last1,First2) and (in the case
C         of bi-linear interpolation) K(Last1,Last2,Last1,Last2)
C
C   Bugs:
C
C      Then there is only one occupied or virtual orbital,
C      some of the integrals will be computed needlessly.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (ONE=1.0D0)
      PARAMETER (TWO=2.0D0)
      PARAMETER (THREE=3.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (THIRD=ONE/THREE)
      PARAMETER (SMALL=1.0D-2)
      PARAMETER (AU=3.67511944138184490995957368D-2)
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
      DIMENSION IG1ST(NMGRP), IGCNT(NMGRP), IGBAS(LDI,*), SHIFT(*)
      DIMENSION C(LDC,*), E(*), DUMP(*)
C
      DO 1000 I1=1,NMGRP-1
          IND1  = IG1ST(I1)
          N1    = IGCNT(I1)
          LAST1 = IND1 + N1 - 1
          DO 900 I2=I1+1,NMGRP
              IND2  = IG1ST(I2)
              N2    = IGCNT(I2)
              LAST2 = IND2 + N2 - 1
              IOFF  = IGBAS(I1,I2)
C
              CALL PSDS2Q(C(1,IND1 ),C(1,IND2),LDC,DUMP,DKFF)
              CALL PSDS2Q(C(1,LAST1),C(1,IND2),LDC,DUMP,DKLF)
              IF( IPRECT.EQ.3 ) THEN
                  SHIFT(IOFF) = DSHIFT * (DKFF + DKLF) * HALF
                  IF( IPRINT.GE.1 ) THEN
                      WRITE(NB6,10200) SHIFT(IOFF)
                  ENDIF
                  CALL DCOPY(N1*N2-1,SHIFT(IOFF),0,SHIFT(IOFF+1),1)
              ELSE
                  CALL PSDS2Q(C(1,LAST1),C(1,LAST2),LDC,DUMP,DKLL)
C
C                 Compute coefficients of bi-linear interpolation
C
                  EOF = AU*E(IND1 )
                  EOL = AU*E(LAST1)
                  EVF = AU*E(IND2 )
                  EVL = AU*E(LAST2)
                  DO  = EOL - EOF
                  DV  = EVL - EVF
                  FO  = ZERO
                  FV  = ZERO
                  IF( DO.GT.SMALL ) FO = (DKLF - DKFF) / DO
                  IF( DV.GT.SMALL ) FV = (DKLL - DKLF) / DV

C                 If one of orbital energies separations is too small,
C                 we'll get averaged F0 gracefully using this formula
C                 instead of straightforward DLFF - EOF*F0 - EVF*FV
C
                  F0 = THIRD * ( DKFF + DKLL + DKLF - 
     .                     FO*(EOF + TWO*EOL) - FV*(TWO*EVF + EVL) )
                  F0 = F0 * DSHIFT
                  FO = FO * DSHIFT
                  FV = FV * DSHIFT
                  IF( IPRINT.GE.1 ) THEN
                      WRITE(NB6,10300) F0, FO, FV
                  ENDIF
C
                  FO = FO * AU
                  FV = FV * AU
                  CALL PSDS2T(N1,N2,E(IND1),E(IND2),F0,FO,FV,
     .                        SHIFT(IOFF))
              ENDIF
C
  900     CONTINUE
 1000 CONTINUE
C
      RETURN
10200 FORMAT(' USING CONSTANT SHIFT VALUE OF ', G20.14, ' A.U.' )
10300 FORMAT(' USING BI-LINEAR INTERPOLATION FOR SHIFT VALUE WITH',
     .       ' COEFFICIENTS (A.U.):'
     .      /1X, G20.14, ' (BASE) ', G20.14, ' (OCC) ', G20.14, 
     .       ' (VAC) ' )
      END
C
      SUBROUTINE PSDS2X(NMGRP,IG1ST,IGCNT,IGBAS,LDI,SHIFT,C,LDC,DUMP)
C
C   Compute single spin part of the CPHF shift vector using Coulomb 
C   integrals as an approximation to the diagonal values of K. This 
C   function might be called for IPRECT=6
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NMGRP  - Number of occupation groups
C      IG1ST  - Indices of the first orbital in a group
C      IGCNT  - Number of orbital in a group
C      IGBAS  - Base indices of orbital groups
C      LDI    - Leading dimension of IGBAS
C      SHIFT  - Computed shift matrix
C      C      - Orbital coefficients
C      LDC    - Leading dimension of C matrix
C      DUMP   - Base of dynamic memory pool
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDYNM - Dynamic memory handles
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C   Module logic:
C
C      Since (ij,ij) in the expression K(ijij) = A*(ij,ij) - (ii,jj)
C      is typically an order of magnitude smaller then (ii,jj),
C      very good shift value can be obtained by using (ii,jj) alone.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON 
     ./PSDOPT/ DSTORE, DSTEP,  DECONV, DPCONV, DCPHF, DPREC, DCDIFF,
     .         DSHIFT, DBASCR, DNCOFF,
     .         IUMIX,  IURHS,  IUK,    IURES,
     .         IPRINT, IENRG,  ICORE,  IDISK,  IMIX,   IDENS,
     .         INDSYM, IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC,
     .         IROWS,  IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP,
     .         IKMODE, ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT
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
      SAVE /PSDOPT/
C
      DIMENSION IG1ST(NMGRP), IGCNT(NMGRP), IGBAS(LDI,*), SHIFT(*)
      DIMENSION C(LDC,*), DUMP(*)
C
      IPLACE = IPSMOF(LAI2T5)
C
      DO 1000 I1=1,NMGRP-1
          IND1  = IG1ST(I1)
          N1    = IGCNT(I1)
          LAST1 = IND1 + N1 - 1
          DO 900 I2=I1+1,NMGRP
              IND2  = IG1ST(I2)
              N2    = IGCNT(I2)
              LAST2 = IND2 + N2 - 1
              IOFF  = IGBAS(I1,I2)
C
              DO 400 IOCC=IND1,LAST1
                  CALL PST1ST(DUMP,C(1,IOCC))
                  CALL PST2ND(DUMP,C(1,IOCC),1)
                  DO 300 IVAC=IND2,LAST2
                      CALL PST3RD(DUMP,1,C(1,IVAC),LDC,1)
                      CALL PST4TH(DUMP,C(1,IVAC),LDC,1,1)
                      SHIFT(IOFF) = -DSHIFT * DUMP(IPLACE)
                      IOFF = IOFF + 1
  300             CONTINUE
  400         CONTINUE
C
  900     CONTINUE
 1000 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSDS2Q(CI,CJ,LDC,DUMP,X)
C
C   Compute single diagonal element of the CPHF K matrix for use
C   in interpolated shift value computation. 
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CI     - Orbital coefficients of the orbital with larger
C               occupation number
C      CJ     - Orbital coefficients of the orbital with smaller
C               occupation number
C      LDC    - Leading dimension of CI and CJ
C      DUMP   - Base of dynamic memory pool
C      X      - Computed value, filled on return
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGBL - Global computation parameters
C      PSDGB2 - Parameters of "response" computation
C      PSDYNM - Dynamic memory handles
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C   Module logic:
C
C      K(i,j,i,j) = 3(ij,ij) - (ii,jj) (RHF)
C      K(i,j,i,j) =  (ij,ij) - (ii,jj) (UHF or NMR)
C
C   Bugs:
C
C      This module is *horribly* inefficient if more then a few
C      K values are to be computed.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (THREE=3.0D0)
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
      SAVE /PSDGBL/
C
      DIMENSION CI(LDC,1), CJ(LDC,1), DUMP(*)
C
      IPLACE = IPSMOF(LAI2T5)
      CALL PST1ST(DUMP,CI)
      CALL PST2ND(DUMP,CI,1)
      CALL PST2ND(DUMP,CJ,2)
      CALL PST3RD(DUMP,1,CJ,LDC,1)
      CALL PST4TH(DUMP,CJ,LDC,1,1)
      AIIJJ = DUMP(IPLACE)
C
C     Doing transformation in this order might slightly improve cache
C     use efficiency.
C
      CALL PST3RD(DUMP,2,CJ,LDC,1)
      CALL PST4TH(DUMP,CI,LDC,1,1)
      AIJIJ = DUMP(IPLACE)
C
      IF( UHF .OR. LIMAG ) THEN
          X = AIJIJ - AIIJJ
      ELSE
          X = THREE * AIJIJ - AIIJJ
      ENDIF
      RETURN
      END
C
      SUBROUTINE PSDS2T(N1,N2,E1,E2,F0,FO,FV,SHIFT)
C
C   Compute shift vector for a single occupation block pair of
C   the single spin using bi-linear interpolation.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N1     - Number of MOs in the first block
C      N2     - Number of MOs in the second block
C      E1     - Eigenvalues for the orbitals of the first block
C      E2     - Eigenvalues for the orbitals of the second block
C      F0,
C      FO,
C      FV     - Coefficients of bi-linear interpolation.
C               Note that F0 is in A.U.s, but FO and FV
C               are in EVs
C      SHIFT  - (Part of) shift matrix to fill
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
C   Module logic:
C
C   Bugs:
C
C     This function computes N1*(N2-1) products needlessly.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION E1(N1), E2(N2), SHIFT(N2,N1)
C
      DO 200 I=1,N1
          FX = F0 + FO*E1(I)
          DO 100 J=1,N2
              SHIFT(J,I) = FX + FV*E2(J)
  100     CONTINUE
  200 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSDS2U(NMGRP,IG1ST,IGCNT,IGBAS,LDI,DOCC,E,SHIFT)
C
C   Crop shift values which are too large for a single spin part of the
C   shift matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NMGRP  - Number of occupation-number groups
C      IG1ST  - First orbital of the occupation group
C      IGCNT  - Number of orbitals in the occupation group
C      IGBAS  - Starting indices of the occupation group pairs
C      LDI    - Leading dimension of the IGBAS
C      DOCC   - Occupation numbers of occupation groups
C      E      - Eigenvalues for one spin
C      SHIFT  - (Part of) shift matrix
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
C   Module logic:
C
C      CPHF Gamma matrix values are computed as:
C
C      Gamma(i,j) = (E(j) - E(i))/(OCC(j) - OCC(i))
C
C      Shift value should not bring diagonal part of K matrix
C      too close to singularity, so our shifts are cropped by
C      0.7 of corresponding Gamma.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (AU=3.67511944138184490995957368D-2)
      PARAMETER (PART=-0.7D0*AU)
C
      DIMENSION IG1ST(NMGRP), IGCNT(NMGRP), DOCC(NMGRP), IGBAS(LDI,*)
      DIMENSION E(*), SHIFT(*)
C
      DO 1000 I1=1,NMGRP-1
          IND1  = IG1ST(I1)
          N1    = IGCNT(I1)
          O1    = DOCC(I1)
          DO 900 I2=I1+1,NMGRP
              IND2  = IG1ST(I2)
              N2    = IGCNT(I2)
              O2    = DOCC(I2)
              IOFF  = IGBAS(I1,I2)
*?????
              GPART = PART / ( O1 - O2 )
C
              DO 200 I=1,N1
                  DO 100 J=1,N2
                      DELTA       = (E(IND2+J-1) - E(IND1+I-1)) * GPART
                      SHIFT(IOFF) = MAX(DELTA,SHIFT(IOFF))
                      IOFF        = IOFF + 1
  100             CONTINUE
  200         CONTINUE
C
  900     CONTINUE
 1000 CONTINUE
C
      RETURN
      END
