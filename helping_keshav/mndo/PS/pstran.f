C     ******************************************************************
C
C     Transform two-electron integrals into the MO basis.
C
C     ******************************************************************
      SUBROUTINE PST1ST(DUMP,CI)
C
C   Transform integrals by the first index, using integrals
C   computed by PSDST1, and store intermediate (i,mu,la,si) 
C   integrals in the safe place. 
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      DUMP   - Base array for temporary values
C      CI     - Orbital coefficients to transform by
C
C   Accessed common blocks:
C
C      PSDYNM - Handles for the dynamic memory blocks
C      PSDGBL - Scratch common for global computation options.
C      PSPRT  - Printing unit
C      PSPRTF - Debug output tuning flags
C      ATOMS  - Starting orbital number, number of orbitals
C               and atom types.
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
C     This function formally is N^3 time consumer, and is cache-
C     -friendly (Ok, I hope it is), so it is unlikely to be a 
C     good targed for performance optimization.
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
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
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      DIMENSION DUMP(*), CI(*)
      SAVE /PSDGBL/, /PSPRTF/, /PSPRT /
C
C*$*INLINE ROUTINE (PST1BU,PST1BT)
CVD$R EXPAND (PST1BU,PST1BT)
C
      IAI2T1 = IPSMOF(LAI2T1)
      IAI2T2 = IPSMOF(LAI2T2)
C
      DO 2000 NA=1,NATOM
          NFRSTA = NFIRST(NA)
          NORBA  = NLAST(NA) - NFIRST(NA) + 1
          NPAIRA = ( NORBA * (NORBA+1) ) / 2
          DO 1500 NB=1,NA-1
              NFRSTB = NFIRST(NB)
              NORBB  = NLAST(NB) - NFIRST(NB) + 1
              NPAIRB = ( NORBB * (NORBB+1) ) / 2
C
C             General diagonal block
C
              CALL PST1BU(NORBA,NPAIRB,CI(NFRSTA),
     .                    DUMP(IAI2T1),NPAIRA,DUMP(IAI2T2))
              IAI2T2 = IAI2T2 + NORBA*NPAIRB
              CALL PST1BT(NORBB,NPAIRA,CI(NFRSTB),
     .                    DUMP(IAI2T1),NPAIRA,DUMP(IAI2T2))
              IAI2T2 = IAI2T2 + NORBB*NPAIRA
              IAI2T1 = IAI2T1 + NPAIRA*NPAIRB
 1500     CONTINUE
C
C         Diagonal blocks are special
C
          CALL PST1BU(NORBA,NPAIRA,CI(NFRSTA),
     .                DUMP(IAI2T1),NPAIRA,DUMP(IAI2T2))
          IAI2T2 = IAI2T2 + NORBA*NPAIRA
          IAI2T1 = IAI2T1 + NPAIRA*NPAIRA
 2000 CONTINUE
      IF(LPA2MO) THEN
          WRITE(NB6,11010)
          CALL PSDPDR(IAI2T1-IPSMOF(LAI2T1),DUMP(IPSMOF(LAI2T1)))
          WRITE(NB6,11000) 
          CALL PSDPDR(IAI2T2-IPSMOF(LAI2T2),DUMP(IPSMOF(LAI2T2)))
      ENDIF
      RETURN
11000 FORMAT(' FIRST-INDEX TRANSFORMED INTEGRALS ARE:'/)
11010 FORMAT(' UNTRANSFORMED INTEGRALS ARE:'/)
      END
C
      SUBROUTINE PST1BU(NORBA,NPAIRB,CA,W,LDW,OUT)
C
C   Transform pair block of two-electron integrals
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NORBA  - Number of orbitals on first atom
C      NPAIRB - Number of electron pairs on second atom
C      CA     - Orbital coefficients for atom A functions
C      W      - Integrals block
C      LDW    - Leading dimension of array W
C      OUT    - Place for transformed integrals.
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
C      None.
C
C   Module logic:
C
C   Bugs:
C
C     This module assumes that first-level cache is big enough to
C     hold atom-atom block of input integrals plus transformed
C     integrals, or the performance will suffer badly.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
C
C*$*ASSERT NOARGUMENT ALIASING
C
      DIMENSION CA(*), W(LDW,*), OUT(NORBA,*)
C
      DO 100 KL=1,NPAIRB
          DO 98 I=1,NORBA
              OUT(I,KL) = ZERO
   98     CONTINUE
  100 CONTINUE
      IJ = 1
      DO 400 I=1,NORBA
          DO 300 J=1,I-1
              CAI = CA(I)
              CAJ = CA(J)
              DO 200 KL=1,NPAIRB
                  OUT(J,KL) = OUT(J,KL) + W(IJ,KL)*CAI
                  OUT(I,KL) = OUT(I,KL) + W(IJ,KL)*CAJ
  200         CONTINUE
              IJ = IJ + 1
  300     CONTINUE
          CAI = CA(I)
          DO 350 KL=1,NPAIRB
              OUT(I,KL) = OUT(I,KL) + W(IJ,KL)*CAI
  350     CONTINUE
          IJ = IJ + 1
  400 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PST1BT(NORBB,NPAIRA,CB,W,LDW,OUT)
C
C   Transform pair block of two-electron integrals
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NORBB  - Number of orbitals on second atom
C      NPAIRA - Number of electron pairs on first atom
C      CB     - Orbital coefficients for atom B functions
C      W      - Integrals block
C      LDW    - Leading dimension of array W
C      OUT    - Place for transformed integrals.
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
C      None.
C
C   Module logic:
C
C   Bugs:
C
C     This module assumes that first-level cache is big
C     enough to hold entire pair block of partially
C     transformed integrals, or the preformance will
C     suffer badly. Untransformed integrals better fit
C     the cache, too, but that's less critical, as they
C     are accesed sequentially.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
C
C*$*ASSERT NOARGUMENT ALIASING
C
      DIMENSION CB(*), W(LDW,*), OUT(NORBB,*)
C
      DO 100 IJ=1,NPAIRA
          DO 98 K=1,NORBB
              OUT(K,IJ) = ZERO
   98     CONTINUE
  100 CONTINUE
      KL = 1
      DO 400 K=1,NORBB
          DO 300 L=1,K-1
              CBL = CB(L)
              CBK = CB(K)
              DO 200 IJ=1,NPAIRA
                  OUT(K,IJ) = OUT(K,IJ) + W(IJ,KL)*CBL
                  OUT(L,IJ) = OUT(L,IJ) + W(IJ,KL)*CBK
  200         CONTINUE
              KL = KL + 1
  300     CONTINUE
          CBK = CB(K)
          DO 350 IJ=1,NPAIRA
              OUT(K,IJ) = OUT(K,IJ) + W(IJ,KL)*CBK
  350     CONTINUE
          KL = KL + 1
  400 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PST2ND(DUMP,CJ,J)
C
C   Transform integrals by set of the second indices, using 
C   results of PST2ND, and store resulting (i,j,la,si) integrals
C   in a slot J of a safe place. 
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      DUMP   - Base array for temporary values
C      CJ     - Orbital coefficients to transform by
C      J      - Slot number to store results to
C
C   Accessed common blocks:
C
C      PSDYNM - Handles for the dynamic memory blocks
C      PSDGBL - Scratch common for global computation options.
C      PSDGB2 - "Response" computation parameters.
C      PSPRT  - Printing unit
C      PSPRTF - Debug output tuning flags
C      ATOMS  - Starting orbital number, number of orbitals
C               and atom types.
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
C      Transformed blocks might be used several times by
C      CPHF matrix K contruction routines, so that all
C      integrals needed for CPHF might be generated 
C      sequentially - hence the orbitals block.
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ONE =1.0D0)
      PARAMETER (ZERO=0.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
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
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      SAVE /PSDGBL/, /PSDGB2/, /PSPRTF/, /PSPRT /
      DIMENSION DUMP(*), CJ(*)
C
C*$*INLINE ROUTINE (DGEMV)
CVD$R EXPAND (DGEMV)
CVD$R SEARCH (dblas2.f)
C
      IAI2T2 = IPSMOF(LAI2T2)
      IAI2T3 = IPSMOF(LAI2T3) + (J-1)*NPAIR
C
      I3A    = IAI2T3
      DO 2000 NA=1,NATOM
          NFRSTA = NFIRST(NA)
          NORBA  = NLAST(NA) - NFIRST(NA) + 1
          NPAIRA = ( NORBA * (NORBA+1) ) / 2
          I3B    = IAI2T3 
C
          DO 100 I=0,NPAIRA-1
              DUMP(I3A+I) = ZERO
  100     CONTINUE
C
          DO 1500 NB=1,NA-1
              NFRSTB = NFIRST(NB)
              NORBB  = NLAST(NB) - NFIRST(NB) + 1
              NPAIRB = ( NORBB * (NORBB+1) ) / 2
C
              CALL DGEMV('T',NORBA,NPAIRB,ONE,DUMP(IAI2T2),NORBA,
     .                   CJ(NFRSTA),1,ONE,DUMP(I3B),1)
              IAI2T2 = IAI2T2 + NORBA*NPAIRB
              CALL DGEMV('T',NORBB,NPAIRA,ONE,DUMP(IAI2T2),NORBB,
     .                   CJ(NFRSTB),1,ONE,DUMP(I3A),1)
              IAI2T2 = IAI2T2 + NORBB*NPAIRA
              I3B = I3B + NPAIRB
 1500     CONTINUE
          CALL DGEMV('T',NORBA,NPAIRA,ONE,DUMP(IAI2T2),NORBA,
     .               CJ(NFRSTA),1,ONE,DUMP(I3A),1)
          IAI2T2 = IAI2T2 + NORBA*NPAIRA
          I3A = I3A + NPAIRA
 2000 CONTINUE
C
      IF(LPA2MO) THEN
          WRITE(NB6,11000) J
          CALL PSDPDR(I3A-IAI2T3, DUMP(IAI2T3))
      ENDIF
C
      RETURN
11000 FORMAT(' SECOND-INDEX TRANSFORMED INTEGRALS (SLOT ', I5, 
     .       ') ARE:'/)
      END
C
      SUBROUTINE PST3RD(DUMP,J,CK,LDC,ICNTK)
C
C   Prepare to transform integrals by the third index.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      DUMP   - Base array for temporary values
C      J      - Index of the orbital transformed on
C               the previous step, should be berween
C               1 and ICNTJ
C      CK     - Orbital coefficients to transform by
C      LDC    - Leading dimension of CK matrix
C      ICNTK  - Number of columns of CK matrix to
C               process
C
C   Accessed common blocks:
C
C      PSDYNM - Handles for the dynamic memory blocks
C      PSDGBL - Scratch common for global computation options.
C      PSDGB2 - "Response" computation parameters.
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
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
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSDGBL/, /PSDGB2/, /PSPRTF/, /PSPRT /
      DIMENSION DUMP(*), CK(LDC,*)
C
      IAI2T3 = IPSMOF(LAI2T3)
      IAI2T4 = IPSMOF(LAI2T4)
      CALL PST3DO(DUMP(IAI2T3+(J-1)*NPAIR),DUMP(IAI2T4),NORBS,
     .            CK,LDC,ICNTK)
      IF(LPA2MO) THEN
          WRITE(NB6,11000) J
          CALL PSDPGM(NORBS,ICNTK,DUMP(IAI2T4),NORBS)
      ENDIF
      RETURN
11000 FORMAT(' THIRD INDEX TRANSFORMED INTEGRALS (J=', I4, 
     .       ') BLOCK IS:'/)
      END
C
      SUBROUTINE PST3DO(AI2,AI3,LDA,CK,LDC,ICNTK)
C
C   Transform integrals by the third index using the
C   results computed by PST2ND, and store intermediate
C   (i,j,k,nu) integrals in a safe place.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      AI2    - Integrals transformed by two indices
C      AI3    - Storage place for three-indices transformed
C               integrals
C      LDA    - Leading dimension of AI3 array
C      CK     - Orbital coefficients to transform by
C      LDC    - Leading dimension of CK matrix
C      ICNTK  - Number of columns of CK matrix to
C               process
C
C   Accessed common blocks:
C
C      PSDGBL - Scratch common for global computation options.
C      ATOMS  - Starting orbital number, number of orbitals
C               and atom types.
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
C   Spedups possible:
C
C      Instantiating inner loop for different number of orbitals
C      might be a good idea. So (though less likely) might be
C      moving of K loop inside innermost loop for large ICNTK.
C
C   Bugs:
C
C      As it is now, this function optimized cache performance
C      at the expence of repeated loop management overhead.
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      SAVE /PSDGBL/
C
      DIMENSION AI2(*), AI3(LDA,*), CK(LDC,*)
C
      DO 2100 K=1,ICNTK
          IIN = 1
          DO 2000 NA=1,NATOM
              NFRSTA = NFIRST(NA)
              NLASTA = NLAST (NA)
C
              DO 1000 I=NFRSTA,NLASTA
                  AI3(I,K) = ZERO
                  DO 800 J=NFRSTA,I-1
                      AI3(I,K) = AI3(I,K) + AI2(IIN)*CK(J,K)
                      AI3(J,K) = AI3(J,K) + AI2(IIN)*CK(I,K)
                      IIN = IIN + 1
  800             CONTINUE
                  AI3(I,K) = AI3(I,K) + AI2(IIN)*CK(I,K)
                  IIN = IIN + 1
 1000         CONTINUE
 2000     CONTINUE
 2100 CONTINUE
      RETURN
      END
C
      SUBROUTINE PST4TH(DUMP,CL,LDC,ICNTK,ICNTL)
C
C   Transform integrals by the forth index using values
C   computed by PST3RD, and store fullly transformed
C   two-electron integrals in the MO basis in a safe
C   place.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      DUMP   - Base array for temporary values
C      CL     - Orbital coefficients to transform by
C      LDC    - Leading dimension of CL matrix
C      ICNTK  - Number of columns of CK processed by
C               the previous call to PST3RD
C      ICNTL  - Number of columns of CL to process.
C
C   Accessed common blocks:
C
C      PSDYNM - Handles for the dynamic memory blocks
C      PSDGBL - Scratch common for global computation options.
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ONE =1.0D0)
      PARAMETER (ZERO=0.0D0)
C
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
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
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSDGBL/, /PSPRTF/, /PSPRT /
      DIMENSION DUMP(*), CL(LDC,*)
C
      IAI2T4 = IPSMOF(LAI2T4)
      IAI2T5 = IPSMOF(LAI2T5)
C
      IF( ICNTL.GT.0 .AND. ICNTK.GT.0 ) THEN
          CALL PSZRMB(ICNTL,ICNTK,DUMP(IAI2T5),ICNTL)
          CALL DGEMM('T','N',ICNTL,ICNTK,NORBS,ONE,CL,LDC,DUMP(IAI2T4),
     .               NORBS,ZERO,DUMP(IAI2T5),ICNTL)
      ENDIF
C
      IF(LPA2MO) THEN
          WRITE(NB6,11000) ICNTK, ICNTL
          CALL PSDPGM(ICNTL,ICNTK,DUMP(IAI2T5),ICNTL)
      ENDIF
      RETURN
11000 FORMAT(' FULLY-TRANSFORMED INTEGRALS BLOCK IS:',2I12/)
      END
C
