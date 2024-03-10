C     ******************************************************************
C
C     Utilities for handling CPHF solution vector.
C
C     ******************************************************************
      SUBROUTINE PSXPRT(X)
C
C   Print quantity in the basis of CPHF-active MOs.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      X      - Quantity to print.
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C      PSOCC  - Occupation numbers classification
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
      SAVE /PSDGBL/, /PSDGB2/, /PSOCC/
C
      DIMENSION X(*)
C
      CALL PSXPR1('ALPHA',X,NMGRPA,IGCNTA,DOCCA,IGBASA,MAXGRP)
      IF(UHF) THEN
          CALL PSXPR1('BETA ',X(IQSZA+1),NMGRPB,IGCNTB,DOCCB,IGBASB,
     .                MAXGRP)
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSXPR1(ID,X,NMGRP,IGCNT,DOCC,IGBAS,LDI)
C
C   Print single spin part of quantity in the basis of CPHF-active MOs.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ID     - Spin name
C      X      - Quantity to print.
C      NMGRP  - Number of occupation number groups
C      IGCNT  - Number of MOs in each group
C      DOCC   - Occupation numbers
C      IGBAS  - Base indices of orbital blocks
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
      CHARACTER*(*) ID
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION X(*), IGCNT(NMGRP), DOCC(NMGRP), IGBAS(LDI,*)
C
      DO 600 I1=1,NMGRP-1
          N1 = IGCNT(I1)
          O1 = DOCC(I1)
          DO 300 I2=I1+1,NMGRP
              N2 = IGCNT(I2)
              O2 = DOCC(I2)
              WRITE(NB6,10000) ID, I1, O1, I2, O2
              CALL PSDPGM(N2,N1,X(IGBAS(I1,I2)),N2)
  300     CONTINUE
  600 CONTINUE
C
      RETURN
10000 FORMAT('   ', A, ' OCCUPATION BLOCKS ', I2, ' (', F5.3, 
     .       ') AND ', I2, ' (', F5.3, '):' )
      END
C
      SUBROUTINE PSXCHB(IP,X)
C
C   Report coupling coefficients exceeding the specified threshold.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IP     - CPHF solution number
C      X      - CPHF solution vector.
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C      PSOCC  - Occupation numbers classification
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
      PARAMETER (XLARGE=8.0D0)
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
      SAVE /PSDGBL/, /PSDGB2/, /PSOCC/
C
      DIMENSION X(*)
C
      CALL PSXCH1(IP,'ALPHA',X,XLARGE,NMGRPA,IGCNTA,IGBASA,
     .            MAXGRP,IG1STA)
      IF(UHF) THEN
          CALL PSXCH1(IP,'BETA ',X(IQSZA+1),XLARGE,NMGRPB,IGCNTB,
     .                IGBASB,MAXGRP,IG1STB)
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSXCH1(IP,ID,X,XLARGE,NMGRP,IGCNT,IGBAS,LDI,IG1ST)
C
C   Print values of the single-spin part of the CPHF equations solution
C   exceeding the specified threshold.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IP     - CPHF solution number
C      ID     - Spin name
C      X      - CPHF solution
C      XLARGE - Threshold
C      NMGRP  - Number of occupation number groups
C      IGCNT  - Number of MOs in each group
C      IGBAS  - Base indices of orbital blocks
C      LDI    - Leading dimension of IGBAS
C      IG1ST  - Starting indices of the MO blocks
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
C      In the present form, PSXCH1 is not parallelizeable. Since
C      I need a small number of the large coupling coefficients 
C      in no particular order, loops below can be easily done 
C      in parallel as long as all operations in the IF statement
C      are protected by a lock.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D0)
C
      CHARACTER*(*) ID
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION X(*), IGCNT(NMGRP), IGBAS(LDI,*), IG1ST(*)
C
      ICOUNT = 0
      VALMAX = ZERO
      DO 600 I1=1,NMGRP-1
          M1 = IG1ST(I1)
          N1 = IGCNT(I1)
          DO 300 I2=I1+1,NMGRP
              M2 = IG1ST(I2)
              N2 = IGCNT(I2)
              IPOS = IGBAS(I1,I2)
              DO 200 IOCC=1,N1
                  DO 100 IVAC=1,N2
                      IF( ABS(X(IPOS)).GE.XLARGE ) THEN
                          ICOUNT = ICOUNT + 1
                          IF( ABS(X(IPOS)).GE.VALMAX ) THEN
                              VALMAX = ABS(X(IPOS))
                              IOCCM  = M1+IOCC-1
                              IVACM  = M2+IVAC-1
                          ENDIF
                      ENDIF
                      IPOS = IPOS + 1
  100             CONTINUE
  200         CONTINUE
  300     CONTINUE
  600 CONTINUE
C
      IF( ICOUNT.GT.0 ) THEN
          WRITE(NB6,10000) ICOUNT, ID, IP, XLARGE, VALMAX, IOCCM, IVACM
      ENDIF
C
      RETURN
10000 FORMAT(/' ', I4, ' COUPLING COEFFICIENTS FROM THE ', A, 
     .        ' PART OF THE CPHF SOLUTION ', I4, ' EXCEED ', G9.3
     .       /'  THE LARGEST VALUE (', G15.7, ') CORRESPONDS TO THE ', 
     .        I5, '->', I5, ' MO TRANSITION.' )
      END
C
C
      SUBROUTINE PSXA2M(CALP,CBET,LDC,XAOA,XAOB,XMO,DUMP)
C
C   Transform response quantity from the AO basis (symmetric packed
C   format) into MO basis.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients. 
C      CBET   - Beta orbital coefficients (not used if UHF flag not set)
C      LDC    - Leading dimension of CALP, CBET matrices
C      XAOA   - Alpha part of the response quantity, AO basis,
C               packed format.
C      XAOB   - Beta part of the response quantity, AO basis
C               (not used if UHF flag not set)
C      XMO    - Response quantity in the MO basis, computed on return
C      DUMP   - Base of dynamic memory items
C
C   Accessed common blocks:
C
C      PSDYNM - Handles for the dynamic memory blocks
C      PSDGBL - Global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C      PSOCC  - Occupation numbers classification
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Transformation is to be computed over all pairs of MO blocks
C      with unequal occupation numbers. In the closed shell and UHF
C      case, this degenerates to occupied-virtual blocks, but it need
C      not be.
C
C      Transformation is two-stage, first over orbitals of the block
C      of orbitals with higher occupation number (N1 orbitals), then 
C      over block with smaller occupation numbers (N2 orbitals).
C
C      1st stage:
C
C            N1              NORBS               N1 
C         <----->           <----->           <----->
C        ^                 ^                 ^
C        |                 |                 |
C        |                 |                 |
C  NORBS |  TMP    = NORBS |   QA   x  NORBS |   C
C        |                 |                 |
C        |                 |                 |
C        v                 v                 v
C
C      2nd stage:
C
C            N1              NORBS               N1 
C         <----->           <----->           <----->
C        ^                 ^                 ^
C        |                 |                 |
C        |                 |    T            |
C   N2   |   QM    =   N2  |   C    x  NORBS |  TMP
C        |                 |                 |
C        |                 |                 |
C        v                 v                 v
C
C      Unfortunately, BLAS Level 3 have no support for multiplying
C      matrices in packed storage, therefore QA is first unpacked
C      into the upper triangle of a general matrix.
C
C      This module does not attempt to compute CPHF vector 0 
C      (right-hand side for the half-electron correction), since
C      it is little more complicated.
C
C   Bugs:
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
      SAVE /PSDGBL/, /PSDGB2/, /PSOCC/
C
      DIMENSION CALP(LDC,*), CBET(LDC,*), XAOA(*), XAOB(*)
      DIMENSION XMO(*), DUMP(*)
C
      IF(LIMAG) THEN
          WRITE(10000)
          STOP 'PSXA2M'
      ENDIF
C
      ITMP = IPSMOF(LA2MTM)
      ITM1 = IPSMOF(LA2MT1)
C
      IF( IQSZA.GT.0 ) THEN
          CALL PSPK2U(NORBS,XAOA,DUMP(ITM1),NORBS)
          CALL PSXAM1(NMGRPA,IG1STA,IGCNTA,IGBASA,MAXGRP,CALP,LDC,
     .                DUMP(ITM1),NORBS,XMO,DUMP(ITMP),NORBS)
      ENDIF
      IF(UHF .AND. IQSZB.GT.0 ) THEN
          CALL PSPK2U(NORBS,XAOB,DUMP(ITM1),NORBS)
          CALL PSXAM1(NMGRPB,IG1STB,IGCNTB,IGBASB,MAXGRP,CBET,LDC,
     .                DUMP(ITM1),NORBS,XMO(1+IQSZA),DUMP(ITMP),NORBS)
      ENDIF
C
      RETURN
10000 FORMAT(' PSXA2M IS NOT IMPLEMENTED IN IMAGINARY MODE ')
      END
C
      SUBROUTINE PSXU2M(CALP,CBET,LDC,XAOA,XAOB,LDX,XMO,DUMP)
C
C   Transform response quantity from the AO basis (stored in the
C   upper triangle of the square matrix) into MO basis.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients. 
C      CBET   - Beta orbital coefficients (not used if UHF flag not set)
C      LDC    - Leading dimension of CALP, CBET matrices
C      XAOA   - Alpha part of the response quantity, AO basis,
C               upper triangle of the square matrix is used.
C      XAOB   - Beta part of the response quantity, AO basis
C               (not used if UHF flag not set)
C      LDX    - Leading dimension of XAOA and XAOB matrices.
C      XMO    - Response quantity in the MO basis, computed on return
C      DUMP   - Base of dynamic memory items
C
C   Accessed common blocks:
C
C      PSDYNM - Handles for the dynamic memory blocks
C      PSDGBL - Global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C      PSOCC  - Occupation numbers classification
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      See PSXA2M
C
C   Bugs:
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
      SAVE /PSDGBL/, /PSDGB2/, /PSOCC/
C
      DIMENSION CALP(LDC,*), CBET(LDC,*), XAOA(LDX,*), XAOB(LDX,*)
      DIMENSION XMO(*), DUMP(*)
C
      IF(LIMAG) THEN
          CALL PSNU2M(CALP,LDC,XAOA,LDX,XMO,DUMP)
          RETURN
      ENDIF
C
      ITMP = IPSMOF(LA2MTM)
C
      IF( IQSZA.GT.0 ) THEN
          CALL PSXAM1(NMGRPA,IG1STA,IGCNTA,IGBASA,MAXGRP,CALP,LDC,
     .                XAOA,LDX,XMO,DUMP(ITMP),NORBS)
      ENDIF
      IF(UHF .AND. IQSZB.GT.0 ) THEN
          CALL PSXAM1(NMGRPB,IG1STB,IGCNTB,IGBASB,MAXGRP,CBET,LDC,
     .                XAOB,LDX,XMO(1+IQSZA),DUMP(ITMP),NORBS)
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSXAM1(NMGRP,IG1ST,IGCNT,IGBAS,LDI,C,LDC,AIN,LDA,OUT,
     .                  TEMP,NORB)
C
C   Transform single-spin quantity in the AO basis into the basis
C   of CPHF-active MOs.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NMGRP  - Number of orbital occupation groups
C      IG1ST  - Index of the first MO in the group
C      IGCNT  - Number of MOs in the group
C      IGBAS  - Base indices of orbital blocks
C      LDI    - Leading dimension of IGBAS array
C      C      - Orbital coefficients
C      LDC    - Leading dimention of the C matrix
C      AIN    - Quantity to transform (A.O. basis)
C      OUT    - Transformed quatity
C      TEMP   - Scratch space
C      NORBS  - Number of AOs.
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
      PARAMETER (ONE=1.0D0)
      PARAMETER (ZERO=0.0D0)
C
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSDGBL/, /PSPRTF/, /PSPRT /
C
      DIMENSION IG1ST(NMGRP), IGCNT(NMGRP), C(LDC,NORB) 
      DIMENSION AIN(LDA,NORB), OUT(*), TEMP(NORB,*)
      DIMENSION IGBAS(LDI,*)
C
      IF(LIMAG) THEN
          CALL PSNAM1(NMGRP,IG1ST,IGCNT,IGBAS,LDI,C,LDC,AIN,LDA,OUT,
     .                TEMP,NORB)
          RETURN
      ENDIF
C
      IF(LPSUMM) THEN
          WRITE(NB6,11000) 
          CALL PSDPSU(NORB,AIN,LDA)
      ENDIF
      DO 500 I1=1,NMGRP-1
          N1 = IGCNT(I1)
          CALL PSZRMB(NORB,N1,TEMP,NORB)
          CALL DSYMM('L','U',NORB,N1,ONE,AIN,LDA,C(1,IG1ST(I1)),
     .               LDC,ZERO,TEMP,NORB)
          IF(LPSUMM) THEN
              WRITE(NB6,11100) I1
              CALL PSDPGM(NORB,N1,TEMP,NORB)
          ENDIF
          DO 300 I2=I1+1,NMGRP
              N2 = IGCNT(I2)
              CALL PSZRMB(N2,N1,OUT(IGBAS(I1,I2)),N2)
              CALL DGEMM('T','N',N2,N1,NORB,ONE,C(1,IG1ST(I2)),LDC,
     .                   TEMP,NORB,ZERO,OUT(IGBAS(I1,I2)),N2)
  300     CONTINUE
  500 CONTINUE
C
      RETURN
11000 FORMAT(' AO QUANTITY TO TRANSFORM IS:' )
11100 FORMAT(' PARTIALLY TRANSFORMED SUB-BLOCK ', I2, ' IS:' )
      END
C
      SUBROUTINE PSNU2M(CALP,LDC,XAOA,LDX,XMO,DUMP)
C
C   Transform antisymmetric matrix in AO basis (supplied as an upper 
C   triangle of square matrix) into (MO) basis of CPHF-active rotations.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients. 
C      LDC    - Leading dimension of CALP matrix
C      XAOA   - Input: upper triangle of the antysimmetric matrix
C               Output: whole matrix is filled.
C      LDX    - Leading dimension of XAOA matrix.
C      XMO    - Output: quantity in MO basis.
C      DUMP   - Base of dynamic memory items
C
C   Accessed common blocks:
C
C      PSDYNM - Handles for the dynamic memory blocks
C      PSDGBL - Global computation options.
C      PSDGB2 - Global computation options needed only if response
C               quantities are computed.
C      PSOCC  - Occupation numbers classification
C      PSPRT  - Printing unit
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      See PSXA2M
C
C   Bugs:
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
      DIMENSION CALP(LDC,*), XAOA(LDX,*)
      DIMENSION XMO(*), DUMP(*)
C
      ITMP = IPSMOF(LA2MTM)
C
      IF( IQSZA.GT.0 ) THEN
          CALL PSNAM1(NMGRPA,IG1STA,IGCNTA,IGBASA,MAXGRP,CALP,LDC,
     .                XAOA,LDX,XMO,DUMP(ITMP),NORBS)
      ENDIF
      IF(UHF) THEN
          WRITE(NB6,11000)
          STOP 'PSNU2M'
      ENDIF
C
      RETURN
11000 FORMAT(' PSNU2M: TRANSFORMATION OF ANTISYMMETRIC MATRICES IS',
     .       ' NOT SUPPORTED IN UHF MODE' )
      END
C
      SUBROUTINE PSNAM1(NMGRP,IG1ST,IGCNT,IGBAS,LDI,C,LDC,AIN,LDA,OUT,
     .                  TEMP,NORBS)
C
C   Transform single-spin part of antisymmetric quantity in the AO basis 
C   into the basis of CPHF-active MOs.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NMGRP  - Number of orbital occupation groups
C      IG1ST  - Index of the first MO in the group
C      IGCNT  - Number of MOs in the group
C      IGBAS  - Base indices of orbital blocks
C      LDI    - Leading dimension of IGBAS array
C      C      - Orbital coefficients
C      LDC    - Leading dimention of the C matrix
C      AIN    - Input: AO quantity to transform, upper triangle filled
C               Output: whole matrix filled.
C      OUT    - Transformed quatity
C      TEMP   - Scratch space
C      NORBS  - Number of AOs.
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
C      First part of transformation makes two times more flops than is
C      strictly necessary, so as to map onto BLAS 3 routines. Depending 
C      on the relative quality of BLAS and compiler-generated code this 
C      may or may not be a win.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ONE=1.0D0)
      PARAMETER (ZERO=0.0D0)
C
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      COMMON
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSPRTF/, /PSPRT /
C
      DIMENSION IG1ST(NMGRP), IGCNT(NMGRP), C(LDC,NORBS) 
      DIMENSION AIN(LDA,NORBS), OUT(*), TEMP(NORBS,*)
      DIMENSION IGBAS(LDI,*)
C
      CALL PSU2AS(AIN,LDA,NORBS)
      IF(LPSUMM) THEN
          WRITE(NB6,11000) 
          CALL PSDPGM(NORBS,NORBS,AIN,LDA)
      ENDIF
      DO 500 I1=1,NMGRP-1
          N1 = IGCNT(I1)
          CALL PSZRMB(NORBS,N1,TEMP,NORBS)
          CALL DGEMM('N','N',NORBS,N1,NORBS,ONE,AIN,LDA,
     .                       C(1,IG1ST(I1)),LDC,ZERO,TEMP,NORBS)
          IF(LPSUMM) THEN
              WRITE(NB6,11100) I1
              CALL PSDPGM(NORBS,N1,TEMP,NORBS)
          ENDIF
          DO 300 I2=I1+1,NMGRP
              N2 = IGCNT(I2)
              CALL PSZRMB(N2,N1,OUT(IGBAS(I1,I2)),N2)
              CALL DGEMM('T','N',N2,N1,NORBS,ONE,C(1,IG1ST(I2)),LDC,
     .                   TEMP,NORBS,ZERO,OUT(IGBAS(I1,I2)),N2)
  300     CONTINUE
  500 CONTINUE
C
      RETURN
11000 FORMAT(' ANTISYMMETRIC AO QUANTITY TO TRANSFORM IS:' )
11100 FORMAT(' PARTIALLY TRANSFORMED SUB-BLOCK ', I2, ' IS:' )
      END
C
      SUBROUTINE PSXAMP(NMGRP,IG1ST,IGCNT,IGBAS,LDI,C,LDC,AIN,LDA,OUT,
     .                  TEMP,NORB,I3,MOK,I4,MOL)
C
C   Partially transform single-spin quantity in the AO basis into the 
C   basis of CPHF-active MOs.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NMGRP  - Number of orbital occupation groups
C      IG1ST  - Index of the first MO in the group
C      IGCNT  - Number of MOs in the group
C      IGBAS  - Base indices of orbital blocks
C      LDI    - Leading dimension of IGBAS array
C      C      - Orbital coefficients
C      LDC    - Leading dimention of the C matrix
C      AIN    - Quantity to transform (A.O. basis)
C      OUT    - Transformed quatity
C      TEMP   - Scratch space
C      NORBS  - Number of AOs.
C      I3     - Occupation block MOK belongs to
C      I4     - Occupation block MOL belongs to, I4 > I3
C      MOK,
C      MOL    - Orbitals identifying ending point of the
C               row to compute, MOL > MOK
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
      PARAMETER (ONE=1.0D0)
      PARAMETER (ZERO=0.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      SAVE /PSDGBL/
C
      DIMENSION IG1ST(NMGRP), IGCNT(NMGRP), C(LDC,NORB) 
      DIMENSION AIN(LDA,NORB), OUT(*), TEMP(NORB,*)
      DIMENSION IGBAS(LDI,*)
C
      IF(LIMAG) THEN
          CALL PSNAMP(NMGRP,IG1ST,IGCNT,IGBAS,LDI,C,LDC,AIN,LDA,OUT,
     .                TEMP,NORB,I3,MOK,I4,MOL)
          RETURN
      ENDIF
C
C    Process complete rows of blocks
C
      DO 200 I1=1,I3-1
          N1 = IGCNT(I1)
          CALL PSZRMB(NORB,N1,TEMP,NORB)
          CALL DSYMM('L','U',NORB,N1,ONE,AIN,LDA,C(1,IG1ST(I1)),
     .               LDC,ZERO,TEMP,NORB)
          DO 100 I2=I1+1,NMGRP
              N2 = IGCNT(I2)
              CALL PSZRMB(N2,N1,OUT(IGBAS(I1,I2)),N2)
              CALL DGEMM('T','N',N2,N1,NORB,ONE,C(1,IG1ST(I2)),LDC,
     .                   TEMP,NORB,ZERO,OUT(IGBAS(I1,I2)),N2)
  100     CONTINUE
  200 CONTINUE
C
C    Process complete blocks on the last row
C
      I1 = I3
      N1 = IGCNT(I1)
      CALL PSZRMB(NORB,N1,TEMP,NORB)
      CALL DSYMM('L','U',NORB,N1,ONE,AIN,LDA,C(1,IG1ST(I1)),
     .           LDC,ZERO,TEMP,NORB)
      DO 300 I2=I1+1,I4-1
          N2 = IGCNT(I2)
          CALL PSZRMB(N2,N1,OUT(IGBAS(I1,I2)),N2)
          CALL DGEMM('T','N',N2,N1,NORB,ONE,C(1,IG1ST(I2)),LDC,
     .               TEMP,NORB,ZERO,OUT(IGBAS(I1,I2)),N2)
  300 CONTINUE
C
C    Process the last, possibly incomplete block
C
      I2     = I4
      N2     = IGCNT(I2)
      N1LAST = MOK - IG1ST(I1) + 1
      N2LAST = MOL - IG1ST(I2) + 1
      N1FULL = N1LAST - 1 + N2LAST/N2
C
C    Largest rectangular subblock
C
      CALL PSZRMB(N2,N1FULL,OUT(IGBAS(I1,I2)),N2)
      CALL DGEMM('T','N',N2,N1FULL,NORB,ONE,C(1,IG1ST(I2)),LDC,
     .           TEMP,NORB,ZERO,OUT(IGBAS(I1,I2)),N2)
      IF( N1FULL.NE.N1LAST ) THEN
C        Remaining subrow
          CALL PSZRVB(N2LAST,OUT(IGBAS(I1,I2)+N1FULL*N2))
          CALL DGEMV('T',NORB,N2LAST,ONE,C(1,IG1ST(I2)),LDC,
     .            TEMP(1,N1LAST),1,ZERO,OUT(IGBAS(I1,I2)+N1FULL*N2),1)
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSNAMP(NMGRP,IG1ST,IGCNT,IGBAS,LDI,C,LDC,AIN,LDA,OUT,
     .                  TEMP,NORBS,I3,MOK,I4,MOL)
C
C   Partially transform single-spin antisymmetric quantity in the AO 
C   basis into the basis of CPHF-active MOs.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NMGRP  - Number of orbital occupation groups
C      IG1ST  - Index of the first MO in the group
C      IGCNT  - Number of MOs in the group
C      IGBAS  - Base indices of orbital blocks
C      LDI    - Leading dimension of IGBAS array
C      C      - Orbital coefficients
C      LDC    - Leading dimention of the C matrix
C      AIN    - Quantity to transform (A.O. basis)
C               Upper triangle should be provided on entry,
C               Matrix is expanded by symmetry on return.
C      OUT    - Transformed quatity
C      TEMP   - Scratch space
C      NORBS  - Number of AOs.
C      I3     - Occupation block MOK belongs to
C      I4     - Occupation block MOL belongs to, I4 > I3
C      MOK,
C      MOL    - Orbitals identifying ending point of the
C               row to compute, MOL > MOK
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
C      Half of flops done on the first transformation stage
C      is actually not necessary.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ONE=1.0D0)
      PARAMETER (ZERO=0.0D0)
C
      DIMENSION IG1ST(NMGRP), IGCNT(NMGRP), C(LDC,NORBS) 
      DIMENSION AIN(LDA,NORBS), OUT(*), TEMP(NORBS,*)
      DIMENSION IGBAS(LDI,*)
C
C    Complement matrix by lower triangle
C
      CALL PSU2AS(AIN,LDA,NORBS)
C
C    Process complete rows of blocks
C
      DO 200 I1=1,I3-1
          N1 = IGCNT(I1)
          CALL PSZRMB(NORBS,N1,TEMP,NORBS)
          CALL DGEMM('N','N',NORBS,N1,NORBS,ONE,AIN,LDA,
     .                       C(1,IG1ST(I1)),LDC,ZERO,TEMP,NORBS)
          DO 100 I2=I1+1,NMGRP
              N2 = IGCNT(I2)
              CALL PSZRMB(N2,N1,OUT(IGBAS(I1,I2)),N2)
              CALL DGEMM('T','N',N2,N1,NORBS,ONE,C(1,IG1ST(I2)),LDC,
     .                   TEMP,NORBS,ZERO,OUT(IGBAS(I1,I2)),N2)
  100     CONTINUE
  200 CONTINUE
C
C    Process complete blocks on the last row
C
      I1 = I3
      N1 = IGCNT(I1)
      CALL PSZRMB(NORBS,N1,TEMP,NORBS)
      CALL DGEMM('N','N',NORBS,N1,NORBS,ONE,AIN,LDA,
     .                   C(1,IG1ST(I1)),LDC,ZERO,TEMP,NORBS)
      DO 300 I2=I1+1,I4-1
          N2 = IGCNT(I2)
          CALL PSZRMB(N2,N1,OUT(IGBAS(I1,I2)),N2)
          CALL DGEMM('T','N',N2,N1,NORBS,ONE,C(1,IG1ST(I2)),LDC,
     .               TEMP,NORBS,ZERO,OUT(IGBAS(I1,I2)),N2)
  300 CONTINUE
C
C    Process the last, possibly incomplete block
C
      I2     = I4
      N2     = IGCNT(I2)
      N1LAST = MOK - IG1ST(I1) + 1
      N2LAST = MOL - IG1ST(I2) + 1
      N1FULL = N1LAST - 1 + N2LAST/N2
C    Largest rectangular subblock
      CALL PSZRMB(N2,N1FULL,OUT(IGBAS(I1,I2)),N2)
      CALL DGEMM('T','N',N2,N1FULL,NORBS,ONE,C(1,IG1ST(I2)),LDC,
     .           TEMP,NORBS,ZERO,OUT(IGBAS(I1,I2)),N2)
      IF( N1FULL.NE.N1LAST ) THEN
C        Remaining subrow
          CALL PSZRVB(N2LAST,OUT(IGBAS(I1,I2)+N1FULL*N2))
          CALL DGEMV('T',NORBS,N2LAST,ONE,C(1,IG1ST(I2)),LDC,
     .            TEMP(1,N1LAST),1,ZERO,OUT(IGBAS(I1,I2)+N1FULL*N2),1)
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSXM2A(CALP,CBET,LDC,XMO,XAOA,XAOB,LDX,DUMP)
C
C   Transform CPHF solution vector into AO basis (upper triangle
C   of square matrix)
C
C   If CPHF solution was already scaled by differences in 
C   occupation numbers, derivatives of the density matrix
C   will be computed.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients
C      CBET   - Beta orbital coefficients (not used if
C               UHF flag not set)
C      LDC    - Leading dimension of CALP and CBET
C      XMO    - CPHF solution vector in MO basis
C      XAOA   - Alpha part of the CPHF solution vector in AO
C               basis, upper triangle is computed.
C      XAOB   - Beta part of the CPHF solution vector in AO
C               basis (not used if UHF flag not set)
C      LDX    - Leading dimension of XAOA and XAOB.
C      DUMP   - Base of dynamic memory items
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C      PSDGB2 - "Response" computation parameters.
C      PSDYNM - Dynamic memory handles
C      PSOCC  - Occupation groups
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
      SAVE /PSDGBL/, /PSDGB2/, /PSOCC/
C
      DIMENSION CALP(LDC,*), CBET(LDC,*), DUMP(*)
      DIMENSION XMO(*), XAOA(LDX,*), XAOB(LDX,*)
C
C    Dynamic memory items
C
      IO2ATM = IPSMOF(LO2ATM)
C
      CALL PSXMA1(CALP,LDC,XMO,XAOA,NORBS,NMGRPA,IG1STA,
     .            IGCNTA,IGBASA,MAXGRP,DUMP(IO2ATM))
      IF(UHF) THEN
          CALL PSXMA1(CBET,LDC,XMO(1+IQSZA),XAOB,NORBS,NMGRPB,
     .                IG1STB,IGCNTB,IGBASB,MAXGRP,DUMP(IO2ATM))
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSXM2P(CALP,CBET,LDC,XMO,XAOA,XAOB,DUMP)
C
C   Transform CPHF solution vector into AO basis (Packed matrix)
C
C   If CPHF solution was already scaled by differences in 
C   occupation numbers, derivatives of the density matrix
C   will be computed.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients
C      CBET   - Beta orbital coefficients (not used if
C               UHF flag not set)
C      LDC    - Leading dimension of CALP and CBET
C      XMO    - CPHF solution vector in MO basis
C      XAOA   - Alpha part of the CPHF solution vector in AO
C               basis, in packed form.
C      XAOB   - Beta part of the CPHF solution vector in AO
C               basis (not used if UHF flag not set)
C      DUMP   - Base of dynamic memory items
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C      PSDGB2 - "Response" computation parameters.
C      PSDYNM - Dynamic memory handles
C      PSOCC  - Occupation groups
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
      SAVE /PSDGBL/, /PSDGB2/, /PSOCC/
C
      DIMENSION CALP(LDC,*), CBET(LDC,*), DUMP(*)
      DIMENSION XMO(*), XAOA(*), XAOB(*)
C
C    Dynamic memory items
C
      IO2ATM = IPSMOF(LO2ATM)
      IO2AT1 = IPSMOF(LO2AT1)
C
      CALL PSXMA1(CALP,LDC,XMO,DUMP(IO2AT1),NORBS,NMGRPA,IG1STA,
     .            IGCNTA,IGBASA,MAXGRP,DUMP(IO2ATM))
      CALL PSU2PK(NORBS,DUMP(IO2AT1),NORBS,XAOA)
      IF(UHF) THEN
          CALL PSXMA1(CBET,LDC,XMO(1+IQSZA),DUMP(IO2AT1),NORBS,NMGRPB,
     .                IG1STB,IGCNTB,IGBASB,MAXGRP,DUMP(IO2ATM))
          CALL PSU2PK(NORBS,DUMP(IO2AT1),NORBS,XAOB)
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSXMA1(C,LDC,X,Y,LDY,NMGRP,IG1ST,IGCNT,IGBAS,LDI,TEMP)
C
C   Tranform single-spin part of the CPHF solution vector into 
C   generalized density matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      C      - Orbital coefficients
C      LDC    - Leading dimension of C
C      X      - CPHF solution vector (already scaled by differences
C               in orbital occupations)
C      Y      - Generalized density
C      LDY    - Leading dimension of Y
C      NMGRP  - Number of occupation groups
C      IG1ST  - First orbital in the occupation group
C      IGCNT  - Number of orbitals in the occupation group
C      IGBAS  - Base indices of the occupation group pairs
C      LDI    - Leading dimension of IGBAS
C      TEMP   - Scratch space
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C      PSPRT  - Printing unit
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
      PARAMETER (ZERO=0.0D0)
      PARAMETER (ONE=1.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      COMMON
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      SAVE /PSDGBL/, /PSPRTF/, /PSPRT /
C
      DIMENSION C(LDC,*), X(*), Y(LDY,*), IG1ST(NMGRP), IGCNT(NMGRP)
      DIMENSION IGBAS(LDI,*), TEMP(*)
C
      IF(LIMAG) THEN
          CALL PSNMA1(C,LDC,X,Y,LDY,NMGRP,IG1ST,IGCNT,IGBAS,LDI,TEMP)
          RETURN
      ENDIF
C
      CALL PSZRM(NORBS,NORBS,Y,LDY)
      DO 1000 I1=1,NMGRP-1
          IND1  = IG1ST(I1)
          N1    = IGCNT(I1)
          DO 900 I2=I1+1,NMGRP
              IND2  = IG1ST(I2)
              N2    = IGCNT(I2)
              IBASE = IGBAS(I1,I2)
              CALL PSZRMB(NORBS,N1,TEMP,NORBS)
              CALL DGEMM ('N','N',NORBS,N1,N2,ONE,C(1,IND2),LDC,
     .                    X(IBASE),N2,ZERO,TEMP,NORBS)
              CALL DSYR2K('U','N',NORBS,N1,ONE,TEMP,NORBS,
     .                    C(1,IND1),LDC,ONE,Y,LDY)
              IF(LPSUMM) THEN
                  WRITE(NB6,11000) I1, I2
                  CALL PSDPGM(NORBS,N1,TEMP,NORBS)
                  WRITE(NB6,11100) I1, I2
                  CALL PSDPSU(NORBS,Y,LDY)
              ENDIF
C
  900     CONTINUE
 1000 CONTINUE
C
      RETURN
11000 FORMAT(' INTERMEDIATE QUANTITY FOR MO BLOCK ',I2,',',I2,':')
11100 FORMAT(' AO QUANTITY AFTER PROCESSING MO BLOCK ',I2,',',I2,':')
      END
C
      SUBROUTINE PSNMA1(C,LDC,X,Y,LDY,NMGRP,IG1ST,IGCNT,IGBAS,LDI,TEMP)
C
C   Tranform single-spin part of the CPHF solution vector into 
C   generalized density matrix in the case of imaginary perturbation.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      C      - Orbital coefficients
C      LDC    - Leading dimension of C
C      X      - CPHF solution vector (already scaled by differences
C               in orbital occupations)
C      Y      - Generalized density
C      LDY    - Leading dimension of Y
C      NMGRP  - Number of occupation groups
C      IG1ST  - First orbital in the occupation group
C      IGCNT  - Number of orbitals in the occupation group
C      IGBAS  - Base indices of the occupation group pairs
C      LDI    - Leading dimension of IGBAS
C      TEMP   - Scratch space
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C      PSPRT  - Printing unit
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
      PARAMETER ( ONE= 1.0D0)
      PARAMETER (ZERO= 0.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      COMMON
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      SAVE /PSDGBL/, /PSPRTF/, /PSPRT /
C
      DIMENSION C(LDC,*), X(*), Y(LDY,*), IG1ST(NMGRP), IGCNT(NMGRP)
      DIMENSION IGBAS(LDI,*), TEMP(*)
C
      CALL PSZRM(NORBS,NORBS,Y,LDY)
C
      DO 1000 I1=1,NMGRP-1
          IND1  = IG1ST(I1)
          N1    = IGCNT(I1)
          DO 900 I2=I1+1,NMGRP
              IND2  = IG1ST(I2)
              N2    = IGCNT(I2)
              IBASE = IGBAS(I1,I2)
              CALL PSZRMB(NORBS,N1,TEMP,NORBS)
              CALL DGEMM('N','N',NORBS,N1,N2,ONE,C(1,IND2),LDC,
     .                           X(IBASE),N2,ZERO,TEMP,NORBS)
              CALL DGEMM('N','T',NORBS,NORBS,N1,ONE,C(1,IND1),LDC,
     .                           TEMP,NORBS,ONE,Y,LDY)
C
*             CALL DGEMM('N','N',NORBS,N1,N2,ONE,C(1,IND2),LDC,X(IBASE),N2,ZERO,TEMP,NORBS)
*             CALL DGEMM('N','T',NORBS,NORBS,N1, ONE,TEMP,NORBS,C(1,IND1),LDC,ONE,Y,LDY)
*             CALL DGEMM('N','T',NORBS,NORBS,N1,SONE,C(1,IND1),LDC,TEMP,NORBS,ONE,Y,LDY)
C
              IF(LPSUMM) THEN
                  WRITE(NB6,11000) I1, I2
                  CALL PSDPGM(NORBS,N1,TEMP,NORBS)
                  WRITE(NB6,11100) I1, I2
                  CALL PSDPGM(NORBS,NORBS,Y,LDY)
              ENDIF
  900     CONTINUE
 1000 CONTINUE
C
C    Antisymmetrize resulting matrix, and we are done:
C
      DO 2000 I2=1,NORBS
          Y(I2,I2) = ZERO
          DO 1900 I1=1,I2-1
              AVE      = Y(I1,I2) - Y(I2,I1)
              Y(I1,I2) =  AVE
              Y(I2,I1) = -AVE
 1900     CONTINUE
 2000 CONTINUE
C
      RETURN
11000 FORMAT(' INTERMEDIATE QUANTITY FOR MO BLOCK ',I2,',',I2,':')
11100 FORMAT(' AO QUANTITY AFTER PROCESSING MO BLOCK ',I2,',',I2,':')
      END
C
