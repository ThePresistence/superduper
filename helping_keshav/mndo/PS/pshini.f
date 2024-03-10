C     ******************************************************************
C
C     Initialization for half-electron derivative calculations.
C     Definition of occupancies and other organisational work.
C
C     ******************************************************************
      SUBROUTINE PSHINI(EALP,EBET)
C
C   Explicit occupancies or half-electron:
C       Classify orbitals into occupation blocks and degenerate
C       shells. Compute permutation table for orbital coefficients
C       reordering.
C
C   Implicit occupancies:
C       Initialize permutation table and occupation blocks for 
C       "do-nothing" case, fill occupation numbers table with
C       explicit values for routines which might need them.
C
C   Half-electron computation:
C       Additionally, compute coefficients of the half-electron
C       correction.
C
C   In all cases, EALP, EBET, OCCA and OCCB are sorted into 
C   "occupancies-blocked" order and have to be restored by
C   caller if their values are necessary for other purposes.
C   Orbital coefficients are *not* sorted here.
C
C   Coded by: Serge Pachkovsky
C
C   Limitations:
C
C      Half-electron correction is only supported for RHF computations.
C
C      Only first derivatives of the half-electron correction are
C      supported at the moment.
C
C      This module might be called before PSDSTR, hence values in PSDOPT
C      (with the exception of IPRINT) and PSDGB2 (with the exception of
C      NOCC[AB] and NVAC[AB]) commons are not available at the moment of 
C      call. Neither are dynamic memory blocks. Unfortunately, this order 
C      can't be reversed without major surgery of PSDSTR
C
C      See also note on EIGTOL and OCCTOL in "Bugs"
C
C   Parameters:
C
C      EALP   - Energy of alpha molecular orbitals, EV. 
C               Are sorted in the descending order of the occupation
C               numbers on return.
C      EBET   - Energy of beta molecular orbitals, EV. 
C               Are sorted in the descending order of the occupation
C               numbers on return.
C
C   Accessed common blocks:
C
C      PSDOPT - Computation parameters
C      PSDGB2 - More system parameters (only NOCC[AB] and NVAC[AB]
C               are valid)
C      PSPRT  - Printing unit
C      PSPRTF - Debugging output control flags
C
C      OCCFL  - Orbitals occupation control flags
C      HALFE  - Half-electron paramaters
C
C   Modified common blocks:
C
C      PSDGBL - System parameters. HALFEL and DORESP flags might
C               be reset if generalised half-electron correction
C               vanish. NOCC[AB] and NVAC[AB], which are not
C               defined consistently by the main code in the case
C               of ther open shell, are redefined as number of orbitals
C               in the occupation groups with occupancies 1.0 and 0.0,
C               respectively.
C      PSOCC  - Occupation number groups description, filled on
C               return as following:
C          DOCCA(MAXGRP)  - Occupation numbers of alpha MOs,
C               one per occupation block. Values are sorted in 
C               the descending order, range is 1.0 to 0.0
C          DOCCB(MAXGRP)  - Occupation numbers of beta MOs.
C          NMGRPA         - Number of alpha occupation blocks.
C          NMGRPB         - Number of beta occupation blocks.
C          IOPN1A         - Index of the first alpha open block,
C                           or zero
C          IOPN1B         - Index of the first beta open block.
C                           or zero
C          IOPNLA         - Index of the last alpha open block.
C                           or zero
C          IOPNLB         - Index of the last beta open block.
C                           or zero
C          IMAXBL         - Length of the longest orbital block.
C          IG1STA(MAXGRP) - First orbital of each alpha occupation
C               block. This index refers to the orbital order after
C               the call to PSSRT[XV], not the original orbital
C               order, since that is not guaranteed to keep occupation
C               blocks together.
C          IGCNTA(MAXGRP) - Number of orbitals in each alpha occupation
C               block.
C          IG1STB(MAXGRP) - First orbitals of beta occupation blocks.
C          IGCNTB(MAXGRP) - Number of orbitals in beta occupation blocks.
C          IGBASA(MAXGRP,MAXGRP) - Base index of the corresponding alpha
C               occupation number group in the CPHF solution vector. All 
C               routines accessing CPHF solution vector should use these 
C               indices rather then make assumptions on the exact structure 
C               wherever at all possible.
C          IGBASB(MAXGRP,MAXGRP) - Base indices of beta groups, 1-based.
C               (i.e., IQSZA have to be added to it to get "real" offset,
C               in the compound vector)
C      PSHALF  - Values related to the half-electron correction, filled
C               on return as following:
C          HLFH(MAXOPN,MAXOPN)   - Weight coefficients of the J(I,J)
C               coulomb integrals in the half-electron correction. Only
C               values with I<=J actually enter energy expression on the
C               half-electron correction, the rest is selected so that
C               HLFH is symmetrical. I and J are relative to the IOPEN1.
C          HLFG(MAXOPN,MAXOPN)   - Weight coefficients of the K(I,J)
C               exchange integrals in the half-electron correction. Only
C               values with I<J actually enter energy expression, the rest
C               is defined so that HLFG is symmetrical and have the same
C               diagonal as the HLFH. I and J are relative to the IOPEN1.
C          NUMOIJ                - Number of Oij matrices which will be
C               needed to compute static part of the half-electron deriv.
C          IOI(MAXACT),
C          IOJ(MAXACT)           - Correspondence table between serial
C               numbers of needed Oij matrices and corresponding (i,j)
C               pairs. i and j are based on IOPEN1.
C          IOPEN1                - Index of the first open MO. This value
C               refer to the orbitals after the call to PSSRT[XV], since
C               original order does not gurantee continuity of the open
C               shells.
C          NUMOPN                - Number of open MOs.
C          NMRACT                - Number of redundant orbital transformation
C               indices which correspond to the non-degenerate pairs of MOs
C               and thus enter the expression for the first derivative of the
C               half-electron correction.
C          NMNACT                - Number of non-redundant orbital transformations
C               which enter the expression for the first derivative of the 
C               half-electron correction.
C          IACTID(MAXOPN,MAXOPN) - Correspondence table between open MO pairs
C               and active redundant variables, MO indices are relative to the
C               IOPEN1. Positive values mark redundant variables, negative
C               values - nonredundant ones.
C          IACTI(MAXACT),
C          IACTJ(MAXACT)         - Correspondence table between active redundant
C               variables and open MO pairs. MO indices are relative to the IOPEN1,
C               I < J.
C      PSORD1  - Correspondence tables for occupation-number ordering of the MOs.
C                See IPSSRI for descriptions.
C      PSORDR  - Permuting tables for the occupation-number ordering of the MOs.
C                See PSMKPT for the description of the table's structure.
C          IPERMA(MAXPRM) - Permutation table for alpha MOs
C          IPERMB(MAXPRM) - Permutation table for beta MOs
C
C      OCCNM   - Explicit occupation numbers. Are sorted in the descending order
C          on return.
C
C   Local storage:
C
C      LMX (4,500 at the moment) INTEGER cells are used for intermediate 
C      correspondence table of MOs.
C
C   Speedups possible:
C
C   Bugs:
C
C      Maximum numbers of different occupation number groups, open shells
C      and open orbitals are compilation-time constants. Although they could
C      have been converted into dynamic parameters, they were not due to the
C      following reasons: 
C          a. Number of different occupation groups and open shells in
C         most practical computations is small and should never exceed
C         built-in limit.
C          b. Validity of half-electron approximation for the case of
C         very large open shells is in doubt.
C          c. Program design is greatly simplified by these assumptions.
C         So probably is object code efficiency.
C
C      Ordering array and permutation tables are of the constant size.
C      Since the relevant constant is included anyway to get access to
C      the occupation numbers and arrays are O(N) only, this is probably 
C      no big deal.
C
C      This function might construct different formulas for the half-
C      electron correction depending on values of OCCTOL (in PSOCCI) and
C      EIGTOL (in PSHLFC and PSHLFG) which are selected somethat arbitrary. 
C      Present values should result in the same expression for the 
C      generalized half-electron correction as the code in forhe/hecor 
C      from MNDO93 for all reasonable inputs. If EIGTOL is set too high,
C      derivatives of the half-electron correction will be imprecise
C      due to neglect of significant redundant rotation coefficients.
C      If EIGTOL is too low, CPHF might become unstable due to inclusion
C      of nearly linear-dependent variables.
C
      USE LIMIT, ONLY: LM1, LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MAXGRP=10)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
      PARAMETER (MAXPRM=2+3*((LMX+1)/2))
C
      PARAMETER (ONE =1.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (SMALL=1.0D-8)
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
     ./PSDOPT/ DSTORE, DSTEP,  DECONV, DPCONV, DCPHF, DPREC, DCDIFF,
     .         DSHIFT, DBASCR, DNCOFF,
     .         IUMIX,  IURHS,  IUK,    IURES,
     .         IPRINT, IENRG,  ICORE,  IDISK,  IMIX,   IDENS,
     .         INDSYM, IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC,
     .         IROWS,  IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP,
     .         IKMODE, ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT
      COMMON
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN), 
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN), 
     .         IACTI(MAXACT), IACTJ(MAXACT)
     ./PSORDR/ IPERMA(MAXPRM), IPERMB(MAXPRM)
     ./PSORD1/ IORDEA(LMX), IORDEB(LMX)
      SAVE /PSDGBL/, /PSDGB2/, /PSPRTF/, /PSOCC /, /PSHALF/, 
     .     /PSORDR/, /PSORD1/, /PSDOPT/, /PSPRT /
C
      COMMON
     ./OCCFL / IMOCC,NOCCAX,NOCCBX,MSUB,MOSUMA,MOSUMB,MOCCA(8),MOCCB(8)
     ./OCCNM / OCCA(LMX),OCCB(LMX)
     ./HALFE / IOD,JOD
     ./INOPT2/ IN2(300)
C
      DIMENSION EALP(*), EBET(*)
      DIMENSION IORDER(LMX)
      EXTERNAL IPSSRI
C
      IF( UHF .AND. (HALFEL.OR.DOCI) ) THEN
          WRITE(NB6,10000)
          STOP 'PSHINI'
      ENDIF
C
      IF( HALFEL.AND.DOCI ) THEN
          WRITE(NB6,10002)
          STOP 'PSHINI'
      ENDIF
C
      IMULT = IN2(66)
C
      IF( IMOCC.LE.1 ) THEN
C
C         Implicit occupation numbers, OCCA and OCCB are not initialized
C         by calling module. So that we can use common-case code, they 
C         are initialized here. 
C
          DO 100 I=1,NOCCA
              OCCA(I) = ONE
  100     CONTINUE
          DO 110 I=NOCCA+1,NORBS
              OCCA(I) = ZERO
  110     CONTINUE
C
          IF(UHF) THEN
              DO 200 I=1,NOCCB
                  OCCB(I) = ONE
  200         CONTINUE
              DO 210 I=NOCCB+1,NORBS
                  OCCB(I) = ZERO
  210         CONTINUE
          ENDIF
          IF( HALFEL.OR.DOCI ) THEN
C             We can have at most two open orbitals in the case of
C             implicit occupation numbers.
              IF(IOD.NE.0) OCCA(IOD) = HALF
              IF(JOD.NE.0) OCCA(JOD) = HALF
          ENDIF
      ENDIF
C
C    Do a little sanity check on the occupation numbers
C    It is Ok to have fractional occupations in CI.
C
      IF( .NOT.(HALFEL.OR.DOCI) ) THEN
          DO 300 I=1,NORBS
              IF( OCCA(I).GT.SMALL .AND. OCCA(I).LT.ONE-SMALL ) THEN
                  WRITE(NB6,10200) I, 'ALPHA', OCCA(I)
                  STOP 'PSHINI'
              ENDIF
  300     CONTINUE
          IF(UHF) THEN
              DO 310 I=1,NORBS
                  IF( OCCB(I).GT.SMALL .AND. OCCB(I).LT.ONE-SMALL ) THEN
                      WRITE(NB6,10200) I, 'BETA', OCCB(I)
                      STOP 'PSHINI'
                  ENDIF
  310         CONTINUE
          ENDIF
      ENDIF
C
      IF(LPHALF) THEN
          WRITE(NB6,11000) IMOCC
          WRITE(NB6,11010) 'ALPHA', ( OCCA(I),I=1,NORBS )
          IF(UHF) THEN
              WRITE(NB6,11010) 'BETA ', ( OCCB(I),I=1,NORBS )
          ENDIF
      ENDIF
C
C     Create sorting permutations table(s) and sort occupation
C     numbers and orbital energies, then fill occupation groups
C     info tables in /PSOCC/
C
      CALL PSSRTO(NORBS,EALP,OCCA,IORDEA)
      IF(LPHALF) THEN
          WRITE(NB6,11020) 'ALPHA', (IORDEA(I),I=1,NORBS)
      ENDIF
C
C     If half-electron computation was requested, IOD and JOD
C     indices have to be converted now.
C
      IF( HALFEL .AND. IMOCC.LT.4 ) THEN
          NIOD = IPSSRI(IOD,IORDEA)
          NJOD = IPSSRI(JOD,IORDEA)
          IF(LPHALF) THEN
              WRITE(NB6,11025) NIOD, NJOD
          ENDIF
      ENDIF
C
C    Correspondence table is destroyed by PSMKPT, and have to be
C    copied to a temporary.
C
      DO 400 I=1,NORBS
          IORDER(I) = IORDEA(I)
  400 CONTINUE
      CALL PSMKPT(NORBS,IORDER,IPERMA)
      IF(LPHALF) THEN
          WRITE(NB6,11030) 'ALPHA', (IPERMA(I),I=2,IPERMA(1)+1)
      ENDIF
      CALL PSSRTX(IPERMA,EALP)
      CALL PSSRTX(IPERMA,OCCA)
      IMAXBL = 0
      CALL PSOCCI(NORBS,OCCA,NMGRPA,DOCCA,IG1STA,IGCNTA,MAXGRP,
     .            NOCCA,NVACA,IOPN1A,IOPNLA,IMAXBL)
      CALL PSOCCB(NMGRPA,IGCNTA,IGBASA,MAXGRP,IQSZA)
      IF(LPHALF) THEN
          WRITE(NB6,11040) 'ALPHA', (OCCA(I),EALP(I),I=1,NORBS)
          WRITE(NB6,11050) 'ALPHA', 
     .                   (DOCCA(I),IG1STA(I),IGCNTA(I),I=1,NMGRPA)
          WRITE(NB6,11053) 'ALPHA'
          CALL PSDPGI(NMGRPA,NMGRPB,IGBASA,MAXGRP)
      ENDIF
      IF(UHF) THEN
          CALL PSSRTO(NORBS,EBET,OCCB,IORDEB)
          IF(LPHALF) THEN
              WRITE(NB6,11020) 'BETA ', (IORDEB(I),I=1,NORBS)
          ENDIF
          DO 410 I=1,NORBS
              IORDER(I) = IORDEB(I)
  410     CONTINUE
          CALL PSMKPT(NORBS,IORDER,IPERMB)
          IF(LPHALF) THEN
              WRITE(NB6,11030) 'BETA ', (IPERMB(I),I=2,IPERMB(1)+1)
          ENDIF
          CALL PSSRTX(IPERMB,EBET)
          CALL PSSRTX(IPERMB,OCCB)
          CALL PSOCCI(NORBS,OCCB,NMGRPB,DOCCB,IG1STB,IGCNTB,MAXGRP,
     .                NOCCB,NVACB,IOPN1B,IOPNLB,IMAXBL)
          CALL PSOCCB(NMGRPB,IGCNTB,IGBASB,MAXGRP,IQSZB)
          IF(LPHALF) THEN
              WRITE(NB6,11040) 'BETA ', (OCCB(I),EBET(I),I=1,NORBS)
              WRITE(NB6,11050) 'BETA ', 
     .                       (DOCCB(I),IG1STB(I),IGCNTB(I),I=1,NMGRPB)
              WRITE(NB6,11053) 'BETA'
              CALL PSDPGI(NMGRPA,NMGRPB,IGBASA,MAXGRP)
          ENDIF
      ELSE
C
C         .NOT.UHF, fake some parameters with "empty" values.
C
          NMGRPB    = 0
          IPERMB(1) = 0
          IQSZB     = 0
      ENDIF
      IQSZ = IQSZA + IQSZB
C
      IF( IQSZ.EQ.0 ) THEN
C
C         Orbital structure of the system won't change no matter what.
C         Response computation is not necessary.
C
          DORESP = .FALSE.
          IF( IPRINT.GE.1 ) WRITE(NB6,10100)
      ENDIF
C
C     Verify occupation blocks energy separation
C     This call was moved to PSDSTR, since it performs checks closely
C     related to the sufficiency of CPHF precision performed there.
C
*     CALL PSCHKE(EALP,EBET)
      IF(HALFEL) THEN
C
C         If half-electron computation was requested, construct 
C         half-electron correction matrices H and G. Classical and
C         generalized corrections are handled differently, since
C         generalized correction is formulated for maximum-spin
C         configuration only.
C
          IF( IMOCC.LE.3 ) THEN
              CALL PSHLFC(IMULT,NIOD,NJOD,EALP)
          ELSE
              CALL PSHLFG(EALP)
          ENDIF
C
C         Count number of required OIJ
C
          NUMOIJ = (NUMOPN*(NUMOPN+1))/2
          IF(LPHALF) THEN
              WRITE(NB6,11055) IOPEN1
              WRITE(NB6,11060)
              CALL PSDPGM(NUMOPN,NUMOPN,HLFH,MAXOPN)
              WRITE(NB6,11070)
              CALL PSDPGM(NUMOPN,NUMOPN,HLFG,MAXOPN)
              WRITE(NB6,11075) NMNACT, NMRACT
              WRITE(NB6,11080) 
              CALL PSDPGI(NUMOPN,NUMOPN,IACTID,MAXOPN)
              WRITE(NB6,11090) (IACTI(I),IACTJ(I),I=1,NMRACT)
              WRITE(NB6,11100) NUMOIJ
*             WRITE(NB6,11110) (IOI(I),IOJ(I),I=1,NUMOIJ)
          ENDIF
      ENDIF
C
C     If response computation was found to be unnecessary, we should
C     undo changes to the orbital energies and occupation numbers.
C
      IF( .NOT.DORESP ) THEN
          CALL PSUSRX(IPERMA,EALP)
          CALL PSUSRX(IPERMA,OCCA)
          IF(UHF) THEN
              CALL PSUSRX(IPERMB,EBET)
              CALL PSUSRX(IPERMB,OCCB)
          ENDIF
      ENDIF
C
C     Uuufff! ;-)
C
      RETURN
C
10000 FORMAT(' HALF-ELECTRON CORRECTION AND CI ARE NOT SUPPORTED ',
     .       'IN UHF MODE')
10002 FORMAT(' BOTH CI AND HALF-ELECTRON MODES ARE SELECTED IN PSHINI.')
10100 FORMAT(' ORBITAL STRUCTURE OF THE SYSTEM WON''T CHANGE, '/
     .       ' HENCE RESPONSE COMPUTATION IS UNNECESSARY.' )
10200 FORMAT(' HALF-ELECTRON FLAG IS NOT SET, BUT THE OCCUPATION ',
     .       ' NUMBER OF THE ', I5, '-TH ', A, ' ORBITAL IS ', F12.8/
     .       ' COMPUTATION IS MEANINGLESS ANS WILL BE ABANDONED.' )
C
C   Debugging formats:
C
11000 FORMAT(' IMOCC = ', I4)
11010 FORMAT(' ', A, ' OCCUPATION NUMBERS ARE:'/ (8(1X,F10.6)))
11020 FORMAT(' ', A, ' INDICES AFTER SORTING BY OCCUPANCIES-BLOCKS:'/
     .       (12(1X,I5)))
11025 FORMAT(' TRANSLATED HALF-ELECTRON ORBITALS ARE: ', I5, 
     .       ' AND ', I5)
11030 FORMAT(' ', A, ' ORBITALS PERMUTATION TABLE: '/(12(1X,I5)))
11040 FORMAT(' SORTED ', A, ' ORBITAL OCCUPATIONS AND ENERGIES: '/ 
     .      (5(4X,F5.3,1X,F10.6)))
11050 FORMAT(' ', A, ' ORBITAL OCCUPATION BLOCKS ARE: '/
     .       (5X,F5.3,' AT ',I5,' LEN ',I5))
11053 FORMAT(' ', A, ' ORBITAL BLOCKS STARTING INDICES ARE: '/)
11055 FORMAT(' OPEN ORBITALS START AT ORBITAL ', I4)
11060 FORMAT(' HALF-ELECTRON CORRECTION H MATRIX: ')
11070 FORMAT(' HALF-ELECTRON CORRECTION G MATRIX: ')
11075 FORMAT(' ', I4, ' NON-REDUNDANT AND ', I4, ' REDUNDANT COUPLING',
     .       ' CONSTANTS ENTER THE HALF-ELECTRON DERIVATIVE' )
11080 FORMAT(' ACTIVE REDUNDANT VARIABLES ARE: ')
11090 FORMAT(' ACTIVE REDUNDANT VARIABLES CORRESPOND TO PAIRS: '/
     .       (5(1X,I4,',',I4,3X)))
11100 FORMAT(' NUMBER OF REQUIRED OIJ IS: ', I5 )
11110 FORMAT(' OIJ PAIRS ARE: '/ (5(1X,I4,',',I4,3X)))
      END
C
      SUBROUTINE PSSRTO(NORBS,E,OCC,IORDER)
C
C   Create correspondence table for the "occupancies-block" sorting
C   of MOs.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NORBS  - Number of orbitals to sort
C      E      - Energy of molecular orbitals sorted in the ascending
C               order
C      OCC    - Occupation numbers of the molecular orbitals
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Table of orbital indices (IORDER) is sorted by two keys.
C      Primary key is occupation number (sort in the descending
C      order), with orbital energy serving as a secondary key
C      (sort in the ascending order). Since most of the elements
C      are already in place, straight insertion algorithm is
C      used, with the following modification: next key is first
C      compared with the previous key, and if elements are in
C      order, sorting continues with the next element. Since 
C      very small number of elements is usually out of order,
C      this results in O(N) process for most practical inputs.
C      Worst-case behaiviour is O(N^2), which is still acceptable.
C
C   Bugs:
C
C      Algorithm used (straight insertion) is only optimal if the
C      table is almost ordered, with unordered values clustered
C      in the middle of the range. Fortunately, this is usually
C      the case for half-electron code and reasonable RHF excited
C      states. Even in the worst case, sorting is only O(N^2) and
C      should be masked by O(N^3) and O(N^4) steps.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ORDRED
C
      DIMENSION E(NORBS), OCC(NORBS), IORDER(NORBS)
C
C     ORDRED(I,J) returns .TRUE. if orbitals I and J are in
C     that order then sorted by "occupancies-block" criteria.
C
      ORDRED(I,J) = OCC(J).LT.OCC(I) .OR.
     .             (OCC(I).EQ.OCC(J) .AND. E(J).GE.E(I))
C
C     Initialize ordering table.
C
      DO 100 I=1,NORBS
          IORDER(I) = I
  100 CONTINUE
C
      DO 1000 KEY=2,NORBS
          IKEY  = IORDER(KEY)
          IF( .NOT. ORDRED(IORDER(KEY-1),IKEY) ) THEN
C
C             New element is not in order, find a place for it by 
C             bisections.
C
              ILEFT  = 0
              IRIGHT = KEY-1
C
C             Insertion could happen either somewhere in between 
C             ILEFT and IRIGHT.
C
  200         CONTINUE
                  IF( IRIGHT - ILEFT .LE. 1 ) GOTO 300
                  IMIDDL = (ILEFT + IRIGHT)/2
                  IF( ORDRED(IORDER(IMIDDL),IKEY) ) THEN
                      ILEFT  = IMIDDL
                  ELSE
                      IRIGHT = IMIDDL
                  ENDIF
                  GOTO 200
  300         CONTINUE
C
C             Shift elements with indices from ILEFT+1 to KEY-1
C             to the right.
C
              DO 400 I=KEY-1,ILEFT+1,-1
                  IORDER(I+1) = IORDER(I)
  400         CONTINUE
C
C             Insert key value at the opened position
C
              IORDER(ILEFT+1) = IKEY
          ENDIF
 1000 CONTINUE
C
      RETURN
      END
C
      INTEGER FUNCTION IPSSRI(IND,IORDER)
C
C   Find new position of the molecular orbital IND in the
C   correspondence table.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IND    - Pre-sort index of the orbital
C      IORDER - Ordering table
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Links starting at the index 'IND' are chased until one
C      leading to the target 'IND' is found. Since IORDER
C      describes permutation, this is guaranteed to happen.
C
C   Bugs:
C
C      Although this function takes constant time for most
C      practical inputs, worst-case execution time is linear 
C      on the size of the IORDER table, which is probably
C      excessive.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION IORDER(*)
C
C     Early return, if possible:
C
      IF( IND.EQ.0 ) THEN
          IPSSRI = 0
          RETURN
      ENDIF
C
      IPLACE = IND
  100 CONTINUE
          ITARGT = IORDER(IPLACE)
          IF( ITARGT.EQ.IND ) THEN
              IPSSRI = IPLACE
              RETURN
          ENDIF
          IPLACE = ITARGT
          GOTO 100
C
      END
C
      SUBROUTINE PSMKPT(NORBS,IORDER,IPERM)
C
C   Create minimal permutations set from the ordering table.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NORBS  - Number of entries in the ordering table.
C      IORDER - Ordering table, destroyed on return.
C      IPERM  - Permutations table, filled on return.
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      All non-degenerate cycles present in the ordering
C      table are extracted as follows:
C       a. Find first element of the IORDER which is not
C      a member of any of the loops already constructed.
C       b. Count number of elements in the loop by following
C      links.
C       c. If length of the loop is greater then one, output
C      corresponding permutation in the backward direction.
C       d. Clear members of the loop in the ordering table.
C      Steps c and d obviously can be omitted for the loops
C      of the length 1, since these elements are already in
C      place.
C
C      This algorithm is O(N) for all possible cases, since
C      each member of the ordering table will be manipulated
C      at most four times and at the least once.
C
C      Created permutation table is in the following format:
C
C          number of entries in the encoding, not counting the
C              final zero
C          loop length, loop members in "comefrom" order
C          ...
C          0
C
C      Since loops of the length zero can never occure, loop
C      length field set to zero means end-of-permutation table.
C
C   Bugs:
C
C      IPERM is silently assumed to be large enough to hold
C      constructed set of permutations. 3*((NORBS+1)/2)+2 
C      should be sufficient in all cases, but is probably
C      pure waste of memory most of the time...
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION IORDER(NORBS), IPERM(*)
C
      IOUT=2
C
      DO 1000 IFIRST=1,NORBS
          INEXT = IORDER(IFIRST)
          IF( INEXT.NE.0 .AND. INEXT.NE.IFIRST ) THEN
C
C             Count & store elements
C
              LENGTH         = 2
              IPERM(IOUT+1)  = IFIRST
              IORDER(IFIRST) = 0
              IOU2           = IOUT + 2
  100         CONTINUE
                  IPERM(IOU2)   = INEXT
                  IOU2          = IOU2 + 1
                  ITMP          = IORDER(INEXT)
                  IORDER(INEXT) = 0
                  IF( ITMP.EQ.IFIRST ) GOTO 200
                  LENGTH        = LENGTH + 1
                  INEXT         = ITMP
                  GOTO 100
  200         CONTINUE
              IPERM(IOUT) = LENGTH
              IOUT        = IOU2
          ENDIF
 1000 CONTINUE
C
      IPERM(IOUT) = 0
      IPERM(1)    = IOUT - 2
      RETURN
      END
C
      SUBROUTINE PSSRTX(IPERM,X)
C
C   Permute array of scalars in the forward direction according 
C   to the control array in IPERM.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IPERM  - Permutations table
C      X      - Array to permute
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Permutations are stored in the "comefrom" order.
C      See PSMKPT for the description.
C
C   Bugs:
C
C      This function is non-parallelizeable in the present form.
C      It need not be. Since all loops in the ordering table are
C      independent, adding table containing starting indices of
C      each loop will allow efficient parallelization on the
C      shared-memory machine.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION IPERM(*), X(*)
C
C     Early return if possible
C
      IF( IPERM(1).EQ.0 ) RETURN
C
      IN = 2
  100 CONTINUE
          LENGTH = IPERM(IN)
          IF( LENGTH.EQ.0 ) RETURN
          IFIRST = IPERM(IN+1)
          XSAVE  = X(IFIRST)
          IPREV  = IFIRST
          DO 200 I=2,LENGTH
              INEXT    = IPERM(IN+I)
              X(IPREV) = X(INEXT)
              IPREV    = INEXT
  200     CONTINUE
          X(IPREV) = XSAVE
          IN = IN + LENGTH + 1
          GOTO 100
C
      END
C
      SUBROUTINE PSUSRX(IPERM,X)
C
C   Permute array of scalars in the backward direction according 
C   to the control array in IPERM (thus reversing the effect of
C   the call to the PSSRTX)
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IPERM  - Permutations table
C      X      - Array to permute
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Permutations are stored in the "comefrom" order,
C      therefore they have to be traced in the reverce
C      order. Since each loop is prefixed with it's
C      length, this is easy deal. See PSMKPT for the
C      description of the permutation table.
C
C   Bugs:
C
C      This function is non-vectorizable in the present form.
C      It need not be, see PSSRTX.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION IPERM(*), X(*)
C
C     Early return if possible
C
      IF( IPERM(1).EQ.0 ) RETURN
C
      IN = 2
  100 CONTINUE
          LENGTH = IPERM(IN)
          IF( LENGTH.EQ.0 ) RETURN
          ILAST  = IPERM(IN+LENGTH)
          XSAVE  = X(ILAST)
          IPREV  = ILAST
          DO 200 I=LENGTH-1,1,-1
              INEXT    = IPERM(IN+I)
              X(IPREV) = X(INEXT)
              IPREV    = INEXT
  200     CONTINUE
          X(IPREV) = XSAVE
          IN = IN + LENGTH + 1
          GOTO 100
C
      END
C
      SUBROUTINE PSSRTV(IPERM,X,LDX,TEMP,LEN)
C
C   Permute array of vectors in the forward direction according 
C   to the control array in IPERM.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IPERM  - Permutations table.
C      X      - Array to permute.
C      LDX    - Leading dimension of the array X.
C      TEMP   - Scratch space to use during copying, should
C               be at least LEN elements long.
C      LEN    - Length of a single vector.
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      See PSSRTX.
C
C   Bugs:
C
C      See PSSRTX.
C
C      Calls to DCOPY is probably not the fastest way of doing
C      moves at least on some architectures.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION IPERM(*), X(LDX,*), TEMP(LEN)
C
C     Early return if possible
C
      IF( IPERM(1).EQ.0 ) RETURN
C
      IN = 2
  100 CONTINUE
          LENGTH = IPERM(IN)
          IF( LENGTH.EQ.0 ) RETURN
          IFIRST = IPERM(IN+1)
          CALL DCOPY(LEN,X(1,IFIRST),1,TEMP,1)
          IPREV  = IFIRST
          DO 200 I=2,LENGTH
              INEXT    = IPERM(IN+I)
              CALL DCOPY(LEN,X(1,INEXT),1,X(1,IPREV),1)
              IPREV    = INEXT
  200     CONTINUE
          CALL DCOPY(LEN,TEMP,1,X(1,IPREV),1)
          IN = IN + LENGTH + 1
          GOTO 100
C
      END
C
      SUBROUTINE PSUSRV(IPERM,X,LDX,TEMP,LEN)
C
C   Permute array of vectors in the reverce direction according 
C   to the control array in IPERM, thus undoing the previous
C   call to PSSRTV.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IPERM  - Permutations table.
C      X      - Array to permute.
C      LDX    - Leading dimension of the array X.
C      TEMP   - Scratch space to use during copying, should
C               be at least LEN elements long.
C      LEN    - Length of a single vector.
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      See PSUSRX.
C
C   Bugs:
C
C      See PSUSRX and PSSRTV.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION IPERM(*), X(LDX,*), TEMP(LEN)
C
C     Early return if possible
C
      IF( IPERM(1).EQ.0 ) RETURN
C
      IN = 2
  100 CONTINUE
          LENGTH = IPERM(IN)
          IF( LENGTH.EQ.0 ) RETURN
          ILAST  = IPERM(IN+LENGTH)
          CALL DCOPY(LEN,X(1,ILAST),1,TEMP,1)
          IPREV  = ILAST
          DO 200 I=LENGTH-1,1,-1
              INEXT    = IPERM(IN+I)
              CALL DCOPY(LEN,X(1,INEXT),1,X(1,IPREV),1)
              IPREV    = INEXT
  200     CONTINUE
          CALL DCOPY(LEN,TEMP,1,X(1,IPREV),1)
          IN = IN + LENGTH + 1
          GOTO 100
C
      END
C
      SUBROUTINE PSOCCI(NORBS,OCC,NMGRP,DOCC,IG1ST,IGCNT,MAXGRP,
     .                  NOCC,NVAC,IOPN1,IOPNL,IMAXBL)
C
C   Classify orbitals sorted by "occupanicies-blocks" into occupation
C   number groups.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NORBS  - Number of orbitals
C      OCC    - Occupation numbers
C      NMGRP  - Number of occupation blocks, filled on return.
C      DOCC   - Block occupation numbers, filled on return.
C      IG1ST  - Index of the first orbital in the block, filled on return.
C      IGCNT  - Number of orbitals in the block, filled on return.
C      MAXGRP - Size of the groups description arrays
C      NOCC   - Number of fully-occupied MOs, filled on return.
C      NVAC   - Number of fully-vacant MOs, filled on return.
C      IOPN1  - Index of the first open block, filled on return.
C      IOPNL  - Index of the last open block, filled on return.
C      IMAXBL - Length of the longest orbital block not counting
C               the last one so far.
C
C   Accessed common blocks:
C
C      PSPRT  - Printing unit
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Blocks of orbitals with occupation numbers between first and
C      last orbitals of the block differing by no more then OCCTOL
C      are extracted.
C
C   Bugs:
C
C       This function assumes that NORBS is at least 1.
C
C       OCCTOL value is arbitrary. In principle, occupation numbers
C       in MNDO93 code are constructed such that comparison on
C       equality should work Ok, so any OCCTOL small enough not
C       to make different occupation numbers equal should work Ok.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (OCCTOL=1.0D-8)
      PARAMETER (ONE=1.0D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION OCC(NORBS), DOCC(MAXGRP), IG1ST(MAXGRP), IGCNT(MAXGRP)
C
      NMGRP  = 0
      IOPN1  = 0
      IOPNL  = 0
      ISTART = 1
      AOCC   = OCC(1)
      DO 100 I=2,NORBS
          IF( ABS(OCC(I)-AOCC).GT.OCCTOL ) THEN
              NMGRP = NMGRP + 1
              IF( NMGRP.GT.MAXGRP ) THEN
                  WRITE(NB6,10000)
                  STOP 'PSOCCI'
              ENDIF
              DOCC(NMGRP)  = AOCC
              IG1ST(NMGRP) = ISTART
              IGCNT(NMGRP) = I - ISTART
              ISTART       = I
              AOCC         = OCC(I)
          ENDIF
  100 CONTINUE
      NMGRP = NMGRP + 1
      IF( NMGRP.GT.MAXGRP ) THEN
           WRITE(NB6,10000)
           STOP 'PSOCCI'
      ENDIF
      DOCC(NMGRP)  = AOCC
      IG1ST(NMGRP) = ISTART
      IGCNT(NMGRP) = NORBS + 1 - ISTART
C
C     Final touch: fill all the additional variables.
C
      IF( ABS(DOCC(1)-ONE).LT.OCCTOL ) THEN
          NOCC  = IGCNT(1)
          IOPN1 = 2
      ELSE
          NOCC  = 0
          IOPN1 = 1
      ENDIF
      IF( ABS(DOCC(NMGRP)).LT.OCCTOL ) THEN
          NVAC  = IGCNT(NMGRP)
          IOPNL = NMGRP - 1
      ELSE
          NVAC  = 0
          IOPNL = NMGRP
      ENDIF
      IF( IOPNL.LT.IOPN1 ) THEN
          IOPN1 = 0
          IOPNL = 0
      ENDIF
C
      DO 500 I=1,NMGRP-1
          IMAXBL = MAX(IMAXBL,IGCNT(I))
  500 CONTINUE
C
      RETURN
10000 FORMAT(' STATIC OCCUPATION CLASSIFICATION TABLE OVERFLOW.'/
     .       ' PROGRAM WILL STOP.' )
      END
C
      SUBROUTINE PSOCCB(NMGRP,IGCNT,IGBAS,LDI,IQSZ)
C
C   Compute starting indices of occupation number blocks in the CPHF
C   solution vector in the standard packing order.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NMGRP  - Number of occupation blocks
C      IGCNT  - Number of orbitals in the block
C      IGBAS  - Base indices of orbital blocks
C      LDI    - Leading dimension of IGBAS
C      IQSZ   - Size of the CPHF solution vector for this spin
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
      DIMENSION IGCNT(NMGRP), IGBAS(LDI,*)
C
      IX = 1
      DO 200 I1=1,NMGRP-1
          DO 100 I2=I1+1,NMGRP
              IGBAS(I1,I2) = IX
              IGBAS(I2,I1) = IX
              IX=IX+IGCNT(I1)*IGCNT(I2)
  100     CONTINUE
  200 CONTINUE
      IQSZ = IX-1
C
      RETURN
      END
C
      SUBROUTINE PSCHKE(EALP,EBET,DIFLOW)
C
C   Verify sufficient occupation blocks energy separation.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      EALP   - Energy of alpha molecular orbitals, EV. 
C      EBET   - Energy of beta molecular orbitals, EV.
C      DIFLOW - Minimal separation of energies of blocks
C               sufficient to ensure precision of CPHF.
C
C   Accessed common blocks:
C
C      PSDGBL - System parameters
C      PSOCC  - Occupation blocks description.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Speedups possible:
C
C   Module logic:
C
C      PSCHKF is called do do the actual verification of alpha
C      (and probably beta) orbital blocks.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MAXGRP=10)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSOCC / DOCCA(MAXGRP), DOCCB(MAXGRP), NMGRPA, NMGRPB, 
     .         IOPN1A, IOPN1B, IOPNLA, IOPNLB, IMAXBL,
     .         IG1STA(MAXGRP), IGCNTA(MAXGRP), 
     .         IG1STB(MAXGRP), IGCNTB(MAXGRP),
     .         IGBASA(MAXGRP,MAXGRP), IGBASB(MAXGRP,MAXGRP)
      SAVE /PSDGBL/, /PSOCC/
C
      DIMENSION EALP(*), EBET(*)
C
      CALL PSCHKF('ALPHA',NMGRPA,DOCCA,IG1STA,IGCNTA,EALP,DIFLOW)
      IF(UHF) THEN
          CALL PSCHKF('BETA ',NMGRPB,DOCCB,IG1STB,IGCNTB,EBET,DIFLOW)
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSCHKF(TAG,NMGRP,DOCC,IG1ST,IGCNT,E,DIFLOW)
C
C   Verify occupation blocks energy separation for a single spin.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      TAG    - Name of the spin block
C      NMGRP  - Number of occupation blocks
C      DOCC   - Occupation numbers
C      IG1ST  - First orbital of the block
C      IGCNT  - Number of orbitals in the block
C      E      - Orbital energies
C      LOWDIF - Minimal separation of energies of blocks
C               sufficient to ensure precision of CPHF.
C
C   Accessed common blocks:
C
C      PSDOPT - Control options (only IPRINT is valid at the moment
C               of call)
C      PSPRT  - Printing unit
C
C   Modified common blocks:
C
C   Local storage:
C
C   Speedups possible:
C
C   Module logic:
C
C      Orbital energies of the orbitals belonging to the different
C      occupation groups are verified:
C        a. not to overlap
C        b. to be separated by a sufficiently large gap.
C
C      Since case (a) is usually an error (and could cause abnormal
C      termination of the CPHF solution process), warning is always
C      issued.
C
C      (b) will often cause CPHF convergence problems and is usually
C      indicative of the wrong multiplicity attributed to the molecule.
C      Warning will be issued if IPRINT is 0 or higher.
C
C   Bugs:
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (EIGTOL=0.3D0)
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
      CHARACTER*(*) TAG
      DIMENSION DOCC(NMGRP), IG1ST(NMGRP), IGCNT(NMGRP), E(*)
C
      DO 200 I=2,NMGRP
          E1L = E(IG1ST(I))
          E1H = E(IG1ST(I)+IGCNT(I)-1)
          DO 100 J=1,I-1
              E2L = E(IG1ST(J))
              E2H = E(IG1ST(J)+IGCNT(J)-1)
C
              IF( E2L.LE.E1H .AND. E2H.GE.E1L ) THEN
                  WRITE(NB6,10000) TAG, DOCC(J), DOCC(I)
              ELSE
                  SEPMIN = MIN( ABS(E1H-E2L), ABS(E1L-E2H) )
                  IF( SEPMIN.LE.EIGTOL ) THEN
                      IF( IPRINT.GT.-1 ) THEN
                          WRITE(NB6,10010) TAG, DOCC(J), DOCC(I)
                      ENDIF
                  ENDIF
                  IF( SEPMIN.LT.DIFLOW ) THEN
                      IF( IPRINT.GT.-1 ) THEN
                          WRITE(NB6,10020) SEPMIN, TAG, DOCC(J), DOCC(I) 
                      ENDIF
                  ENDIF
              ENDIF
C
  100     CONTINUE
  200 CONTINUE
C
      RETURN
10000 FORMAT(/' ENERGIES OF THE ', A, ' OCCUPATION BLOCKS WITH ',
     .        'OCCUPATIONS ', F5.3, ' AND ', F5.3, ' OVERLAP'/
     .        ' VALIDITY OF CPHF SOLUTION IS IN DOUBT '/)
10010 FORMAT( ' ', A, ' ORBITAL BLOCKS WITH OCCUPATION NUMBERS ',
     .        F5.3, ' AND ', F5.3, ' ARE NEARLY DEGENERATE ' )
10020 FORMAT(/' ', G14.7, ' EV SEPARATION BETWEEN ', A, ' ORBITAL',
     .        ' BLOCKS WITH '
     .       /' OCCUPATION NUMBERS ', F5.3, ' AND ', F5.3,
     .        ' MAY BE INSUFFICIENT TO GIVE THE DESIRED CPHF PRECISION')
      END
C
      SUBROUTINE PSHLFC(IMULT,IOD,JOD,E)
C
C   Verify consistency of orbital occupation information and half-electron
C   information. Compute coefficients of the half-electron correction for
C   the classical case(s).
C
C   Coded by: Serge Pachkovsky
C
C   Limitations:
C
C      The only configurations supported are excited singlet, doublet
C      and triplet. 
C
C   Parameters:
C
C      IMULT  - Multiplicity of the system
C      IOD    - Index of the first MO carrying half-electron
C      JOD    - Index of the second MO carrying half-electron or zero
C      E      - Orbital energies
C
C   Accessed common blocks:
C
C      PSOCC  - Occupation blocks description.
C      PSPRT  - Printing unit
C
C   Modified common blocks:
C
C      PSHALF - Half-electron correction parameters, filled on return.
C
C   Local storage:
C
C   Module logic:
C
C      Occupation blocks are verifyed to agree with one of the
C      standard half-electron cases. That is: 
C        a. There could be at most three occupation blocks.
C           Only one of them can have fractional occupation number of 0.5
C        b. Number of orbitals in the fractional occupation block
C           should be either 1 (if JOD is zero and IOD is not) or 2
C           otherwise. IOD and JOD should be the indices of these orbitals.
C        c. For the multiplicity 2, JOD should be zero.
C
C      If these conditions are satisfied, HLFH and HLFG matrices are filled.
C      Then multplicity is either 1 or 3, and orbitals are non-degenerate,
C      corresponding coupling coefficient is added to the active redundant
C      set.
C
C   Gotchas:
C
C      IOD and JOD are *not* the same with the values in the HALFE common
C      of the MNDO93. 
C
C      CPHF treatment of half-electron correction for multiplicities 
C      1 and 3 is quite different depending on degeneracy or non-degeneracy
C      of orbitals.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (OCCTOL=1.0D-8)
      PARAMETER (EIGTOL=0.3D0)
      PARAMETER (MAXGRP=10)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
C
      PARAMETER (ONE =1.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ZERO=0.0D0)
C
      COMMON
     ./PSOCC / DOCCA(MAXGRP), DOCCB(MAXGRP), NMGRPA, NMGRPB, 
     .         IOPN1A, IOPN1B, IOPNLA, IOPNLB, IMAXBL,
     .         IG1STA(MAXGRP), IGCNTA(MAXGRP), 
     .         IG1STB(MAXGRP), IGCNTB(MAXGRP),
     .         IGBASA(MAXGRP,MAXGRP), IGBASB(MAXGRP,MAXGRP)
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN), 
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN), 
     .         IACTI(MAXACT), IACTJ(MAXACT)
     ./PSPRT / NB6
      SAVE /PSOCC/, /PSHALF/ , /PSPRT /
C
      DIMENSION E(*)
C
C     Check the consistency of the IOD, JOD and the orbitals occupations
C     blocks. There could be at most three occupation blocks, but only
C     one can have fractional occupation of 0.5. Block with fractional
C     occupations should consist entirely of IOD and JOD orbitals.
C
      IBLFRC = 2
      IF( ABS(DOCCA(1)-ONE).GT.OCCTOL ) THEN
          IBLFRC = 1
      ENDIF
      IF( ABS(DOCCA(IBLFRC)-HALF).GT.OCCTOL ) THEN
          WRITE(NB6,10010) (DOCCA(I),I=1,NMGRPA)
          STOP 'PSHLFC'
      ENDIF
      IF( NMGRPA.EQ.IBLFRC+1 ) THEN
          IF( ABS(DOCCA(NMGRPA)).GT.OCCTOL ) THEN
              WRITE(NB6,10010) (DOCCA(I),I=1,NMGRPA)
              STOP 'PSHLFC'
          ENDIF
      ELSEIF( NMGRPA.NE.IBLFRC ) THEN
          WRITE(NB6,10010) (DOCCA(I),I=1,NMGRPA)
          STOP 'PSHLFC'
      ENDIF
C
      IF( IGCNTA(IBLFRC).EQ.1 ) THEN
          IF( IG1STA(IBLFRC).NE.IOD .OR. JOD.NE.0 ) THEN
              WRITE(NB6,10020) IOD, JOD, IG1STA(IBLFRC), 0
              STOP 'PSHLFC'
          ENDIF
          IF( IMULT.NE.2 ) THEN
              WRITE(NB6,10025) IMULT
              STOP 'PSHLFC'
          ENDIF
      ELSEIF( IGCNTA(IBLFRC).EQ.2 ) THEN
          IF( IG1STA(IBLFRC).NE.IOD .OR. JOD.NE.IOD+1) THEN
             WRITE(NB6,10020) IOD, JOD, IG1STA(IBLFRC), IG1STA(IBLFRC)+1
             STOP 'PSHLFC'
          ENDIF
          IF( IMULT.NE.1 .AND. IMULT.NE.3 ) THEN
              WRITE(NB6,10027) IMULT
              STOP 'PSHLFC'
          ENDIF
      ELSE
          WRITE(NB6,10030) IGCNTA(IBLFRC)
          STOP 'PSHLFC'
      ENDIF
C
C     Half-electron data is consistent, proceed with building correction
C
      IOPEN1 = IOD
      IF( IMULT.EQ.1 ) THEN
C
C         Excited singlet
C
          NUMOPN      = 2
          HLFH(1,1)   = -0.25D0
          HLFH(2,2)   = -0.25D0
          HLFH(1,2)   = ZERO
          HLFH(2,1)   = ZERO
          HLFG(1,1)   = -0.25D0
          HLFG(2,2)   = -0.25D0
          HLFG(1,2)   = 1.5D0
          HLFG(2,1)   = 1.5D0
      ELSEIF( IMULT.EQ.2 ) THEN
C
C         Doublet
C
          NUMOPN      = 1
          HLFH(1,1)   = -0.25D0
          HLFG(1,1)   = -0.25D0
          NMRACT      = 0
          IACTID(1,1) = 0
      ELSEIF( IMULT.EQ.3 ) THEN
C
C         Triplet
C
          NUMOPN      = 2
          HLFH(1,1)   = -0.25D0
          HLFH(2,2)   = -0.25D0
          HLFH(1,2)   = ZERO
          HLFH(2,1)   = ZERO
          HLFG(1,1)   = -0.25D0
          HLFG(2,2)   = -0.25D0
          HLFG(1,2)   = -0.5D0
          HLFG(2,1)   = -0.5D0
      ELSE
          WRITE(NB6,10090) IMULT
          STOP 'PSHLFC'
      ENDIF
      IF( IMULT.EQ.1 .OR. IMULT.EQ.3 ) THEN
          IF( ABS(E(IOD)-E(JOD)).LT.EIGTOL ) THEN
C
C             Orbitals are degenerate
C
              NMRACT      = 0
              IACTID(1,1) = 0
              IACTID(2,1) = 0
              IACTID(1,2) = 0
              IACTID(2,2) = 0
          ELSE
C
C             Orbitals are non-degenerate
C
              NMRACT      = 1
              IACTID(1,1) = 0
              IACTID(2,1) = 1
              IACTID(1,2) = 1
              IACTID(2,2) = 0
              IACTI(1)    = 1
              IACTJ(1)    = 2
          ENDIF
      ENDIF
C
C     Non-redundant variables do not enter classical-case
C     expressions.
C
      NMNACT = 0
      RETURN
10010 FORMAT(' STRUCTURE OF ORBITAL OCCUPATION BLOCKS IS INCOMPATIBLE',
     .       ' WITH CLASSICAL'/
     .       ' HALF-ELECTRON CORRECTION. OCCUPANCIES ARE:'/
     .       (12(1X,F5.3)))
10020 FORMAT(' HALF-ELECTRON ORBITALS ARE ', I4, ' AND ', I4, ' WHILE ',
     .       I4, ' AND ', I4/
     .       ' ARE EXPECTED FROM OCCUPATION NUMBERS')
10025 FORMAT(' MULTIPLICITY OF ', I4, ' IS NOT FEASIBLE FOR SINGLE ',
     .       'OPEN SHELL')
10027 FORMAT(' MULTIPLICITY OF ', I4, ' IS NOT FEASIBLE FOR TWO ',
     .       'OPEN SHELLS')
10030 FORMAT(' NUMBER OF ORBITALS CARRYING HALF-ELECTRON (',I4,') IS ',
     .       ' NEITHER 1 NOR 2')
10090 FORMAT(' CLASSICAL HALF-ELECTRON CORRECTION IS NOT DEFINED FOR ',
     .       'THE MULTIPLICITY ', I4 )
      END
C
      SUBROUTINE PSHLFG(E)
C
C   Compute coefficients of the generalized half-electron correction.
C
C   Coded by: Serge Pachkovsky
C
C   Limitations:
C
C      Only configurations with the highest spin are supported.
C
C   Parameters:
C
C      E      - Orbital energies.
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options (only IPRINT is valid)
C      PSPRTF - Debug output control flags
C      PSOCC  - Occupation number blocks
C      PSHALF - Half-electron parameters, filled on return
C      PSPRT  - Printing unit
C
C   Modified common blocks:
C
C      PSDGBL - Global computation parameters (HALFEL and
C               DORESP flags might be reset if correction
C               vanish).
C
C   Local storage:
C
C      4*MAXOPN (40 at the moment) DOUBLE PRECISION cells for
C      occupation numbers of the open shells.
C
C      2*MAXOPN (20 at the moment) INTEGER cells for the
C      description of the open orbitals structure.
C
C   Module logic:
C
C      Orbitals are considered to be in the same shell if separation
C      of the orbital energies of the two consequitive orbitals
C      is less than or equal to EIGTOL. Therefore, block of M orbitals
C      could in principle span up to M*EIGTOL EVs.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MAXGRP=10)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
      PARAMETER (EIGTOL=0.3D0)
      PARAMETER (OCCTOL=1.0D-8)
C
      PARAMETER (ZERO=0.0D0)
      PARAMETER (ONE=1.0D0)
      PARAMETER (TWO=2.0D0)
      PARAMETER (FOUR=4.0D0)
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
      COMMON
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      SAVE /PSDOPT/, /PSPRTF/, /PSOCC/, /PSHALF/, /PSDGBL/, /PSPRT /
C
      DIMENSION E(*)
      DIMENSION ISH1ST(MAXOPN), ISHM(MAXOPN), DSHN(MAXOPN)
      DIMENSION DSHNA(MAXOPN), DSHNB(MAXOPN), DSHO(MAXOPN)
C
C     Split occupation blocks with fractional occupancies into 
C     (quasy-)degenerate open shells.
C
      IOPEN1 = 0
      ISHCNT = 0
      DO 200 IBLK=1,NMGRPA
          GROCC = DOCCA(IBLK)
          IF( ABS(GROCC-ONE).GT.OCCTOL .AND. ABS(GROCC).GT.OCCTOL ) THEN
              IF( IOPEN1.EQ.0 ) IOPEN1 = IG1STA(IBLK)
              IBASE = IG1STA(IBLK)
              I1    = 1
              DO 100 I=2,IGCNTA(IBLK)
                  IF( ABS(E(IBASE+I-2)-E(IBASE+I-1)).GT.EIGTOL ) THEN
                      ISHCNT = ISHCNT + 1
                      IF( ISHCNT.GT.MAXOPN ) THEN
                          WRITE(NB6,10000)
                          STOP 'PSHLFG'
                      ENDIF
                      ISH1ST(ISHCNT) = I1 + IBASE - IOPEN1
                      ISHM(ISHCNT)   = I - I1
                      DSHO(ISHCNT)   = GROCC
                      DSHN(ISHCNT)   = TWO * ISHM(ISHCNT) * GROCC
                      I1 = I
                  ENDIF
  100         CONTINUE
              ISHCNT = ISHCNT + 1
              IF( ISHCNT.GT.MAXOPN ) THEN
                  WRITE(NB6,10000)
                  STOP 'PSHLFG'
              ENDIF
              ISH1ST(ISHCNT) = I1 + IBASE - IOPEN1
              ISHM(ISHCNT)   = IGCNTA(IBLK) + 1 - I1
              DSHO(ISHCNT)   = GROCC
              DSHN(ISHCNT)   = TWO * ISHM(ISHCNT) * GROCC
          ENDIF
  200 CONTINUE
C
C    Determine alpha occupations of the degenerate shells and
C    total number of the open orbitals. Warn if the number of
C    electrons on the open shell is not natural.
C
      NUMOPN = 0
      DO 400 I=1,ISHCNT
          NUMOPN = NUMOPN + ISHM(I)
          DSHNA(I) = MIN(DBLE(ISHM(I)),DSHN(I))
          DSHNB(I) = DSHN(I) - DSHNA(I)
          IF( ABS( NINT(DSHN(I)) - DSHN(I) ).GT.OCCTOL ) THEN
              IF( IPRINT.GT.-5 ) THEN
                  WRITE(NB6,10010) E(IOPEN1+ISH1ST(I)-1), DSHN(I)
              ENDIF
          ENDIF
  400 CONTINUE
      IF(LPHALF) THEN
          WRITE(NB6,11000) NUMOPN,(ISH1ST(I),ISHM(I),DSHN(I),DSHNA(I),
     .                           DSHNB(I),I=1,ISHCNT)
      ENDIF
      IF( NUMOPN.GT.MAXOPN ) THEN
          WRITE(NB6,10020)
          STOP 'PSHLFG'
      ENDIF
C
C    Clear H, G and IACTID matrices
C
      DO 600 I=1,NUMOPN
         DO 500 J=1,NUMOPN
            HLFH(J,I)   = ZERO
            HLFG(J,I)   = ZERO
            IACTID(J,I) = 0
  500    CONTINUE
  600 CONTINUE
C
C    Fill H and G matrices (i.e., weights of coulomb and exchange
C    integrals in the half-electron correction.
C
      DO 1300 ISHA=1,ISHCNT
C
C         Intra-shell contributions go first
C
          IA  = ISH1ST(ISHA)
          AM  = ISHM(ISHA)
          AN  = DSHN(ISHA)
          ANA = DSHNA(ISHA)
          ANB = DSHNB(ISHA)
          HI  = -(AN/AM)**2/FOUR + ANB/AM
          IF( ISHM(ISHA).GT.1 ) THEN
              SA = -TWO/(AM*(AM-ONE))
              HJ = SA*(AN*(TWO*AM-AN)/(TWO*AM) - (ANA-ANB)/TWO)
              HK = SA*(((ANA-ANB)/TWO)**2 - (AN*(TWO*AM-AN))/(FOUR*AM))
          ENDIF
          DO 800 I=1,ISHM(ISHA)
              HLFH(IA+I-1,IA+I-1) = HI
              HLFG(IA+I-1,IA+I-1) = HI
              DO 700 J=1,I-1
                  HLFH(IA+I-1,IA+J-1) = HJ
                  HLFH(IA+J-1,IA+I-1) = HJ
                  HLFG(IA+I-1,IA+J-1) = HK
                  HLFG(IA+J-1,IA+I-1) = HK
  700         CONTINUE
  800     CONTINUE
C
C         Now, compute inter-shell contributions
C
          DO 1200 ISHB=1,ISHA-1
              IB  = ISH1ST(ISHB)
              BM  = ISHM(ISHB)
              BNA = DSHNA(ISHB)
              BNB = DSHNB(ISHB)
              HKX = -(ANA-ANB)*(BNA-BNB)/(TWO*AM*BM)
              DO 1100 I=1,ISHM(ISHA)
                  DO 1000 J=1,ISHM(ISHB)
                      HLFG(IA+I-1,IB+J-1) = HKX
                      HLFG(IB+J-1,IA+I-1) = HKX
 1000             CONTINUE
 1100         CONTINUE
 1200     CONTINUE
 1300 CONTINUE
C
C    Finally, fill active redundant and non-redundant variables matrix.
C
      NMRACT = 0
      NMNACT = 0
      DO 2000 ISHA=2,ISHCNT
          IA = ISH1ST(ISHA)
          MA = ISHM(ISHA)
          DO 1900 ISHB=1,ISHA-1
              IB = ISH1ST(ISHB)
              MB = ISHM(ISHB)
              DO 1800 I=1,MA
                  DO 1700 J=1,MB
                      IF( ABS(DSHO(ISHA)-DSHO(ISHB)).LE.OCCTOL ) THEN
C
C                         Occupations are equal, hence variable is redundant
C
                          NMRACT = NMRACT+1
                          IACTID(IA+I-1,IB+J-1) = NMRACT
                          IACTID(IB+J-1,IA+I-1) = NMRACT
                          IACTI(NMRACT)         = IB+J-1
                          IACTJ(NMRACT)         = IA+I-1
                      ELSE
C
C                         Different occupations
C
                          NMNACT = NMNACT + 1
                          IACTID(IA+I-1,IB+J-1) = -NMNACT
                          IACTID(IB+J-1,IA+I-1) = -NMNACT
                      ENDIF
 1700             CONTINUE
 1800         CONTINUE
 1900     CONTINUE
 2000 CONTINUE
C
C    If half-electron correction had degenerated to nothing, clear
C    half-electron and "response" flags, since this is actually
C    closed-shell computation.
C
      IF( NUMOPN.EQ.0 ) THEN
          IF( IPRINT.GT.-5 ) WRITE(NB6,10900)
          HALFEL = .FALSE.
          IF( MODE.LE.1 ) DORESP = .FALSE.
      ENDIF
      RETURN
10000 FORMAT(' MAXIMUM NUMBER OF OPEN SHELLS EXCEEDED.' )
10010 FORMAT(' NUMBER OF ELECTRONS IN THE OPEN SHELL WITH ENERGY ',
     .       G12.6, ' EV (', F8.6, ') IS NOT NATURAL' )
10020 FORMAT(' MAXIMUM NUMBER OF OPEN ORBITALS EXCEEDED.' )
10900 FORMAT(' COMPUTATION IS ACTUALLY CLOSED-SHELL, NO CPHF STEP ',
     .       'WILL BE PERFORMED' )
C
C    Debugging formats:
C
11000 FORMAT(' THERE ARE ', I4, ' OPEN ORBITALS WITH THE FOLLOWING ',
     .       'SHELL STRUCTURE:'/
     .       (' SHELL AT ', I5, ' (', I3, '-FOLD DEGENERATE) HAS ',
     .       F6.4, ' ELECTRONS (', F6.4, ' ALPHA, ', F6.4, ' BETA)') )
      END
C
      SUBROUTINE PSSRTM(CALP,CBET,LDC)
C
C   Sort molecular orbitals into occupation-block order.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients. Modified.
C      CBET   - Beta orbital coefficients. Modified.
C      LDC    - Leading dimension of the CALP and CBET.
C
C   Accessed common blocks:
C
C      PSDGBL - System parameters.
C      PSORDR - Permuting tables for the occupation-number ordering of the MOs.
C
C   Modified common blocks:
C
C   Local storage:
C
C      LMX (4,500 at the moment) DOUBLE PRECISION cells are used for temporary
C      storage of single MO vector.
C
C   Speedups possible:
C
C   Bugs:
C
C      Local storage of the temporary MO vector is ugly. This can be fixed by
C      doing multiply passes over CALP and CBET with reasonably small buffer,
C      or by allocating buffer dynamically.
C
      USE LIMIT, ONLY: LM1, LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MAXPRM=2+3*((LMX+1)/2))
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSORDR/ IPERMA(MAXPRM), IPERMB(MAXPRM)
      SAVE /PSDGBL/, /PSORDR/
C
      DIMENSION CALP(LDC,*), CBET(LDC,*)
      DIMENSION TEMP(LMX)
C
      CALL PSSRTV(IPERMA,CALP,LDC,TEMP,NORBS)
      IF(UHF) THEN
          CALL PSSRTV(IPERMB,CBET,LDC,TEMP,NORBS)
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSUSRT(EALP,EBET,CALP,CBET,LDC)
C
C   Restore ordiring of the MO coefficients, orbital energies and 
C   occupation numbers.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      EALP   - Energy of alpha molecular orbitals, EV. Modified.
C      EBET   - Energy of beta molecular orbitals, EV. Modified.
C      CALP   - Alpha orbital coefficients. Modified.
C      CBET   - Beta orbital coefficients. Modified.
C      LDC    - Leading dimension of the CALP and CBET.
C
C   Accessed common blocks:
C
C      PSDGBL - System parameters.
C      PSORDR - Permuting tables for the occupation-number ordering of the MOs.
C
C   Modified common blocks:
C
C      OCCNM  - Explicit occupation numbers. 
C
C   Local storage:
C
C      LMX (4,500 at the moment) DOUBLE PRECISION cells are used for temporary
C      storage of single MO vector.
C
C   Speedups possible:
C
C   Bugs:
C
C      Local storage of the temporary MO vector is ugly. This can be fixed by
C      doing multiply passes over CALP and CBET with reasonably small buffer,
C      or by allocating buffer dynamically.
C
      USE LIMIT, ONLY: LM1, LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MAXPRM=2+3*((LMX+1)/2))
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSORDR/ IPERMA(MAXPRM), IPERMB(MAXPRM)
      SAVE /PSDGBL/, /PSORDR/
C
      COMMON
     ./OCCNM / OCCA(LMX),OCCB(LMX)
C
      DIMENSION EALP(*), EBET(*), CALP(LDC,*), CBET(LDC,*)
      DIMENSION TEMP(LMX)
C
      CALL PSUSRX(IPERMA,EALP)
      CALL PSUSRX(IPERMA,OCCA)
      CALL PSUSRV(IPERMA,CALP,LDC,TEMP,NORBS)
      IF(UHF) THEN
          CALL PSUSRX(IPERMB,EBET)
          CALL PSUSRX(IPERMB,OCCB)
          CALL PSUSRV(IPERMB,CBET,LDC,TEMP,NORBS)
      ENDIF
C
      RETURN
      END
C
      BLOCK DATA PSORDI
C
C   Initialize permutation table to random values, to ensure program
C   crash if uninitialized data is used.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C   Accessed common blocks:
C
C      PSORDR - Permuting tables for the occupation-number ordering of the MOs.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Speedups possible:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MAXPRM=2+3*((LMX+1)/2))
C
      COMMON
     ./PSORDR/ IPERMA(MAXPRM), IPERMB(MAXPRM)
      SAVE /PSORDR/
C
      DATA IPERMA( 1)/-1047166573/, IPERMB( 1)/-1047166573/
      DATA IPERMA( 2)/-1465677437/, IPERMB( 2)/-1465677437/
      DATA IPERMA( 3)/-2121056124/, IPERMB( 3)/-2121056124/
      DATA IPERMA( 4)/-2053535323/, IPERMB( 4)/-2053535323/
      DATA IPERMA( 5)/-2119988827/, IPERMB( 5)/-2119988827/
      DATA IPERMA( 6)/-2052964215/, IPERMB( 6)/-2052964215/
      DATA IPERMA( 7)/-1802005627/, IPERMB( 7)/-1802005627/
      DATA IPERMA( 8)/-1803184733/, IPERMB( 8)/-1803184733/
      DATA IPERMA( 9)/-2054930302/, IPERMB( 9)/-2054930302/
      DATA IPERMA(10)/-1472142715/, IPERMB(10)/-1472142715/
      DATA IPERMA(11)/-1719171776/, IPERMB(11)/-1719171776/
      DATA IPERMA(12)/-0679378040/, IPERMB(12)/-0679378040/
      DATA IPERMA(13)/-1835620958/, IPERMB(13)/-1835620958/
      DATA IPERMA(14)/-1834467291/, IPERMB(14)/-1834467291/
C
      END
