C     ******************************************************************
C
C     NMR-related computational routines.
C     Control routines are collected in psnmrc.f.
C
C     ******************************************************************
      SUBROUTINE PSNST1(DUMP)
C
C   Compute and store two-center two-electron repulsion integrals
C   for use by CPHF routines.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      DUMP   - Scratch array. 
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C      PSDYNM - Dynamic memory handles
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
C      2116 DOUBLE PRECISION cells for two-electron repulsion
C      integrals.
C
C   Module logic:
C
C      This a bare bones version of PSDST1.
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
      SAVE /PSDGBL/, /PSPRTF/, /PSPRT /
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
C
      DIMENSION DUMP(*)
      DIMENSION AINTS(46,46), DUMMY(1)
C
      IAI2T1 = IPSMOF(LAI2T1) - 1
C
C    Loop over all atom pairs
C
      DO 2100 NA=1,NATOM
          NORBA  = NLAST(NA) - NFIRST(NA) + 1
          NPAIRA = ( NORBA * (NORBA+1) ) / 2
          NTYPA  = NAT(NA)
C
          DO 2090 NB=1,NA-1
              NORBB  = NLAST(NB) - NFIRST(NB) + 1
              NPAIRB = ( NORBB * (NORBB+1) ) / 2
              CALL PSINTS(0,NA,NB,AINTS,DUMMY,DUMMY)
              IF(LPINTS) THEN
                  WRITE(NB6,10100) NA, NB
                  WRITE(NB6,10110) NORBA, NORBB, NPAIRA, NPAIRB
                  WRITE(NB6,10400)
                  CALL PSDPGM(NPAIRA,NPAIRB,AINTS(2,2),46)
              ENDIF
              DO 1010 J=1,NPAIRB
                  DO 1008 I=1,NPAIRA
                      DUMP(IAI2T1+I) = AINTS(1+I,1+J)
 1008             CONTINUE
                  IAI2T1 = IAI2T1 + NPAIRA
 1010         CONTINUE
 2090     CONTINUE
          CALL PSINT1(NTYPA,NPAIRA,AINTS)
          DO 2098 J=1,NPAIRA
               DO 2096 I=1,NPAIRA
                   DUMP(IAI2T1+I) = AINTS(1+I,1+J)
 2096          CONTINUE
               IAI2T1 = IAI2T1 + NPAIRA
 2098     CONTINUE
          IF(LPINTS) THEN
              WRITE(NB6,10600) NA
              CALL PSDPGM(NPAIRA,NPAIRA,AINTS(2,2),46)
          ENDIF
 2100 CONTINUE
C
      RETURN
C
C    Debugging formats
C
10100 FORMAT(' PROCESSING PAIR ', I5, ', ', I5 )
10110 FORMAT(' NORBA  = ', I5, ' NORBB  = ', I5,
     .       ' NPAIRA = ', I5, ' NPAIRB = ', I5)
10400 FORMAT(' TWO-CENTER TWO-ELECTRON INTEGRALS:' )
10600 FORMAT(' ONE-CENTER TWO-ELECTRON INTEGRALS FOR ATOM ', I5, ':' )
      END
C
      SUBROUTINE PSNHA0(HA,LDH1,LDH2)
C
C   Compute (upper triangle of) static first-order magnetic field 
C   derivatives of the hamiltonian matrix (intended for use as an 
C   RHS in magnetic field CPHF computation). Derivatives are in
C   atomic units multiplied by the speed of light.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      HA     - Output derivatives of the one-electron Hamiltonian
C      LDH1   - First leading dimension of the HA matrix
C      LDH2   - Second leading dimension of the HA matrix
C
C   Accessed common blocks:
C
C      ATOMS  - Range of orbitals carried by an atom
C      ATOMC  - Coordinates of atoms (don't you wish it was in a.u.s, boy?)
C      PSDOPT - We need INDSYM from here.
C      PSDGBL - Number of atoms, for consistency with derivatives code.
C      PSPRT  - Printing unit.
C      PSPRTF - Debugging output control flags.
C
C   Local storage:
C
C      7*(MAXL+1)**4 (ca. 570) DOUBLE PRECISION cells.
C
C   Module logic:
C
C   Bugs:
C
C      The code is not particularly cache-friendly.
C      Expressions were originally written for lower-triangle part
C      of HA and then hacked for upper triangle.
C
C      Lower triangle will be destroyed.
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      LOGICAL DEBUG
      PARAMETER (MAXL=2)
      PARAMETER (BOHR=0.529167D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (PT25=0.25D0)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ATOMC / COORD(3,LM1)
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
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSDGBL/, /PSPRTF/, /PSDOPT/, /PSPRT /
C
      DIMENSION HA(LDH1,LDH2,3)
C
      DIMENSION DH((MAXL+1)**2,(MAXL+1)**2),
     .          DL((MAXL+1)**2,(MAXL+1)**2,3), 
     .          DLT((MAXL+1)**2,(MAXL+1)**2,3),
     .          DD((MAXL+1)**2,(MAXL+1)**2,3)
C
      DEBUG = LPINTS .AND. LPNMR
C
      DO 5000 IA=1,NATOM
          NORBA = NLAST(IA)-NFIRST(IA)+1
          IBASA = NFIRST(IA)-1
          XA    = COORD(1,IA)/BOHR
          YA    = COORD(2,IA)/BOHR
          ZA    = COORD(3,IA)/BOHR
          IBMIN = 1
          IF( INDSYM.EQ.2 ) IBMIN = IA + 1
          DO 4000 IB=IBMIN,NATOM
              IF( IB.EQ.IA ) GOTO 3999
C
C            Two-center contributions, i.e. dipole moment, H and orbital
C            moment contributions.
C
              NORBB = NLAST(IB)-NFIRST(IB)+1
              IBASB = NFIRST(IB)-1
              XB    = COORD(1,IB)/BOHR
              YB    = COORD(2,IB)/BOHR
              ZB    = COORD(3,IB)/BOHR
              DX    = XA-XB
              DY    = YA-YB
              DZ    = ZA-ZB
              RX    = YA*ZB - YB*ZA
              RY    = XB*ZA - XA*ZB
              RZ    = XA*YB - YA*XB
C            Compute necessary integrals. PSNA2H and PSNA2I are
C            the alternate forms of the H^{00} contributions.
C            PSNA2H uses MNDO zetas and betas, while PSNA2I uses
C            NMR parameters. Using PSNA2I without making similar
C            changes to the zero-th order Hamiltonian in the SCF
C            part will cause loss of translational invariance.
              CALL PSNA2H(IA,IB,-DX,-DY,-DZ,DH)
C             CALL PSNA2I(IA,IB,-DX,-DY,-DZ,DH)
C            PSNA2J is somewhat slower than the *H and *I routines,
C            but handles different types of resonance integrals (H and I
C            can only do MNDO-type resonance integrals)
C             CALL PSNA2J(IA,IB,-DX,-DY,-DZ,DH)
              CALL PSNA2L(IA,IB,-DX,-DY,-DZ,DL,DLT)
              CALL PSNA2D(IA,IB,-DX,-DY,-DZ,DD)
              IF( DEBUG ) THEN
                  WRITE(NB6,11000) IA, IB
                  WRITE(NB6,11010)
                  CALL PSDPGM(NORBA,NORBB,DH(1,1),((MAXL+1)**2))
                  WRITE(NB6,11020) 'X'
                  CALL PSDPGM(NORBA,NORBB,DL (1,1,1),((MAXL+1)**2))
                  CALL PSDPGM(NORBB,NORBA,DLT(1,1,1),((MAXL+1)**2))
                  WRITE(NB6,11020) 'Y'
                  CALL PSDPGM(NORBA,NORBB,DL (1,1,2),((MAXL+1)**2))
                  CALL PSDPGM(NORBB,NORBA,DLT(1,1,2),((MAXL+1)**2))
                  WRITE(NB6,11020) 'Z'
                  CALL PSDPGM(NORBA,NORBB,DL (1,1,3),((MAXL+1)**2))
                  CALL PSDPGM(NORBB,NORBA,DLT(1,1,3),((MAXL+1)**2))
                  WRITE(NB6,11030) 'X'
                  CALL PSDPGM(NORBA,NORBB,DD(1,1,1),((MAXL+1)**2))
                  WRITE(NB6,11030) 'Y'
                  CALL PSDPGM(NORBA,NORBB,DD(1,1,2),((MAXL+1)**2))
                  WRITE(NB6,11030) 'Z'
                  CALL PSDPGM(NORBA,NORBB,DD(1,1,3),((MAXL+1)**2))
              ENDIF
C
C            Add contributions
C
              DO 2000 M1=1,NORBA
                  IH1 = IBASA + M1
                  DO 1900 M2=1,NORBB
                      IH2 = IBASB + M2
                      HA(IH1,IH2,1) = HALF*(
     .                  - DZ * DD(M1,M2,2) + DY * DD(M1,M2,3)
     .                  + RX * DH(M1,M2)
     .                  - HALF*(DL(M1,M2,1)-DLT(M2,M1,1))
     .                  )
                      HA(IH1,IH2,2) = HALF*(
     .                  - DX * DD(M1,M2,3) + DZ * DD(M1,M2,1)
     .                  + RY * DH(M1,M2)
     .                  - HALF*(DL(M1,M2,2)-DLT(M2,M1,2))
     .                  )
                      HA(IH1,IH2,3) = HALF*(
     .                  - DY * DD(M1,M2,1) + DX * DD(M1,M2,2)
     .                  + RZ * DH(M1,M2)
     .                  - HALF*(DL(M1,M2,3)-DLT(M2,M1,3))
     .                  )
 1900             CONTINUE
 2000         CONTINUE
 3999         CONTINUE
 4000     CONTINUE
C
C   One-center contributions, only orbital moment part is non-zero
C
          CALL PSNHA1(IA,DL)
          DO 4300 M2=1,NORBA
              IH2 = IBASA + M2
              DO 4200 M1=1,NORBA
                  IH1 = IBASA + M1
                  HA(IH1,IH2,1) = -(DL(M1,M2,1)-DL(M2,M1,1)) * PT25
                  HA(IH1,IH2,2) = -(DL(M1,M2,2)-DL(M2,M1,2)) * PT25
                  HA(IH1,IH2,3) = -(DL(M1,M2,3)-DL(M2,M1,3)) * PT25
 4200         CONTINUE
 4300     CONTINUE
 5000 CONTINUE
C
      IF( INDSYM.NE.2 ) THEN
C        Antisymmetrize HA0 matrices.
          CALL PSNMAT('HA0','X',' ',HA(1,1,1),LDH1)
          CALL PSNMAT('HA0','Y',' ',HA(1,1,2),LDH1)
          CALL PSNMAT('HA0','Z',' ',HA(1,1,3),LDH1)
      ENDIF
      IF(LPNMR) THEN
          WRITE(NB6,11140) 'X'
          CALL PSDPAU(NORBS,HA(1,1,1),LDH1)
          WRITE(NB6,11140) 'Y'
          CALL PSDPAU(NORBS,HA(1,1,2),LDH1)
          WRITE(NB6,11140) 'Z'
          CALL PSDPAU(NORBS,HA(1,1,3),LDH1)
      ENDIF
      RETURN
C
11000 FORMAT(' PAIR ', I5, ',', I5, ':')
11010 FORMAT(' CORE HAMILTONIAN BLOCK IS: ')
11020 FORMAT(' ORBITAL MOMENT BLOCK ',A,' IS: ')
11030 FORMAT(' DIPOLE MOMENT BLOCK ',A,' IS: ')
11040 FORMAT(' CURRENT H DERIVATIVE ',A,' IS: ')
11140 FORMAT(' FINAL H MAGNETIC FIELD DERIVATIVE (AO BASIS) ',A,' IS: ')
11050 FORMAT(' ATOM ', I5, ':' )
      END
C
      SUBROUTINE PSNA2H(IA,IB,DX,DY,DZ,DH)
C
C   Compute "normal" a two-center block of one particle hamiltonian 
C   matrix for HA0 (MNDO version).
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Index of the first atom.
C      IB     - Index of the second atom.
C      DX,DY,DZ
C             - Relative coordinates of the atom B (a.u.s)
C      DH     - Output H_{\mu\nu} integrals.
C
C   Accessed common blocks:
C
C      PAROPT - Orbital exponents and beta parameters for S- and P-type
C               orbitals
C      DPARM1 - Orbital exponents and beta parameters for D-type orbitals
C      ATOMS  - Orbital numbers and atomic types
C      DNBND  - Major quantum numbers
C
C   Local storage:
C
C      3*(MAXL+1)**4 + 2*(MAXL+1)**2 + 6*(MAXL+1) DOUBLE PRECISION
C      cell, i.e. 423 doublewords for MAXL=2.
C
C   Module logic:
C
C   Bugs:
C
C      Orbital moments are limited to 2 by arbitrary dimentioned
C      temporaries, while the code (except the fragment extracting
C      atomic parameters from MNDO94's COMMON blocks) can in 
C      principle handle arbitrary L's.
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (HALFAU=0.183755972069092245497978684D-1)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./DNBND / III(LMZ),IIID(LMZ)
     ./DPARM1/ UDD(LMZ),ZD(LMZ),BETAD(LMZ)
     ./PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ),BETAP(LMZ),
     .         ALP(LMZ)
      COMMON
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
     ./PSPRT / NB6
      SAVE /PSNPAR/, /PSPRT /
C
      DIMENSION DH((MAXL+1)**2,(MAXL+1)**2)
C
      DIMENSION ZA(0:MAXL), ZB(0:MAXL), BETAA(0:MAXL), BETAB(0:MAXL),
     .          NA(0:MAXL), NB(0:MAXL)
C   Emergency stop
      IF( IBETA.NE.1 ) STOP 'PSNA2H'
C   Gather parameters
      NORBA = NLAST(IA)-NFIRST(IA)+1
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPA = NAT(IA)
      ITYPB = NAT(IB)
C
      ZA   (0) = ZS   (ITYPA)
      BETAA(0) = BETAS(ITYPA)
      NA   (0) = III  (ITYPA)
      IF( NORBA.GT.1 ) THEN
          ZA   (1) = ZP   (ITYPA)
          BETAA(1) = BETAP(ITYPA)
          NA   (1) = III  (ITYPA)
      ENDIF
      IF( NORBA.GT.4 ) THEN
          ZA   (2) = ZD   (ITYPA)
          BETAA(2) = BETAD(ITYPA)
          NA   (2) = IIID (ITYPA)
      ENDIF
      IF( NORBA.GT.9 ) THEN
          WRITE(NB6,10000) NORBA
          STOP 'PSNA2H'
      ENDIF
C
      ZB   (0) = ZS   (ITYPB)
      BETAB(0) = BETAS(ITYPB)
      NB   (0) = III  (ITYPB)
      IF( NORBB.GT.1 ) THEN
          ZB   (1) = ZP   (ITYPB)
          BETAB(1) = BETAP(ITYPB)
          NB   (1) = III  (ITYPB)
      ENDIF
      IF( NORBB.GT.4 ) THEN
          ZB   (2) = ZD   (ITYPB)
          BETAB(2) = BETAD(ITYPB)
          NB   (2) = IIID (ITYPB)
      ENDIF
      IF( NORBB.GT.9 ) THEN
          WRITE(NB6,10000) NORBB
          STOP 'PSNA2H'
      ENDIF
C   Compute integrals
      LB     = 0
      IBASB  = 1
  100 CONTINUE
          LA    = 0
          IBASA = 1
  200     CONTINUE
              CALL PSOVX(NA(LA),LA,ZA(LA),NB(LB),LB,ZB(LB),DX,DY,DZ,
     .                   DH(IBASA,IBASB),(MAXL+1)**2)
              SC  = (BETAA(LA)+BETAB(LB))*HALFAU
              DO 500 M2=IBASB,IBASB+2*LB
                  DO 400 M1=IBASA,IBASA+2*LA
                      DH(M1,M2) = SC*DH(M1,M2)
  400             CONTINUE
  500         CONTINUE
          IBASA = IBASA + 2*LA+1
          LA    = LA    + 1
          IF( IBASA.LT.NORBA ) GOTO 200
      IBASB  = IBASB  + 2*LB+1
      LB     = LB     + 1
      IF( IBASB.LT.NORBB ) GOTO 100
C    Reordering of the integrals
      CALL PSSX94(NORBA,NORBB,DH,(MAXL+1)**2)
C
      RETURN
10000 FORMAT(' NUMBER OF AOS IS TOO BIG (',I3,') IN PSNA2H')
      END
C
      SUBROUTINE PSNA2L(IA,IB,DX,DY,DZ,DL,DLT)
C
C   Compute <\mu|L|\nu> and <\nu|L|\mu> integrals for a two-center
C   block of HA0.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Index of the first atom
C      IB     - Index of the second atom
C      DX,DY,DZ
C             - Relative coordinates of the atom B (a.u.s)
C      DL     - Output <\mu|L^\nu|\nu> integrals
C      DLT    - Output <\nu|L^\mu|\mu> integrals
C
C   Accessed common blocks:
C
C      ATOMS  - Orbital numbers and atomic types
C      PSNPAR - NMR-related orbital exponents and beta parameters
C
C   Local storage:
C
C      3*(MAXL+1)**4 + 2*(MAXL+1)**2 + 6*(MAXL+1) DOUBLE PRECISION
C      cell, i.e. 423 doublewords for MAXL=2.
C
C   Module logic:
C
C   Bugs:
C
C      Orbital moments are limited to 2 by arbitrary dimentioned
C      temporaries, while the code (except the fragment extracting
C      atomic parameters from MNDO94's COMMON blocks)can in 
C      principle handle arbitrary L's.
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION DL((MAXL+1)**2,(MAXL+1)**2,3), 
     .          DLT((MAXL+1)**2,(MAXL+1)**2,3)
C
      DIMENSION SN(2*MAXL+1,2*MAXL+1), ST(2*MAXL+1,2*MAXL+1)
C
C   Gather necessary parameters
C
      NORBA = NLAST(IA)-NFIRST(IA)+1
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPA = NAT(IA)
      ITYPB = NAT(IB)
C
C   Compute integrals. Warning: ST serves a double purpose
C   here, first as a temporary for elementary overlap integrals
C   and second as a transposed set of complete overlap integrals.
C   Since L acts only on the angular part of the function, it
C   is Ok to contract elementary integrals in the inner loop.
C
      DO 800 LA=0,MAXLA(ITYPA)
          IBASA = LA**2 + 1
          DO 700 LB=0,MAXLA(ITYPB)
              IBASB = LB**2 + 1
              DO 200 MB=1,2*LB+1
                  DO 100 MA=1,2*LA+1
                      SN(MA,MB) = ZERO
  100             CONTINUE
  200         CONTINUE
              DO 600 ILA=1,MAXN(LA,ITYPA)
                  NA = NN(ILA,LA,ITYPA)
                  ZA = ZN(ILA,LA,ITYPA)
                  CA = CN(ILA,LA,ITYPA)
                  DO 500 ILB=1,MAXN(LB,ITYPB)
                      NB = NN(ILB,LB,ITYPB)
                      ZB = ZN(ILB,LB,ITYPB)
                      CB = CN(ILB,LB,ITYPB)
                      CALL PSOVX(NA,LA,ZA,NB,LB,ZB,DX,DY,DZ,ST,2*MAXL+1)
                      DO 400 MB=1,2*LB+1
                          DO 300 MA=1,2*LA+1
                              SN(MA,MB) = SN(MA,MB) + CA*CB*ST(MA,MB)
  300                     CONTINUE
  400                 CONTINUE
  500             CONTINUE
  600         CONTINUE
              DO 620 MB=1,2*LB+1
                  DO 610 MA=1,2*LA+1
                      ST(MB,MA) = SN(MA,MB)
  610             CONTINUE
  620         CONTINUE
              CALL PSOLX(2*LA+1,LB,SN,2*MAXL+1,DL (IBASA,IBASB,1),
     .                                         (MAXL+1)**2,(MAXL+1)**2)
              CALL PSOLX(2*LB+1,LA,ST,2*MAXL+1,DLT(IBASB,IBASA,1),
     .                                         (MAXL+1)**2,(MAXL+1)**2)
  700     CONTINUE
  800 CONTINUE
C
C    Reordering of the integrals
C
      CALL PSSX94(NORBA,NORBB,DL (1,1,1),(MAXL+1)**2)
      CALL PSSX94(NORBA,NORBB,DL (1,1,2),(MAXL+1)**2)
      CALL PSSX94(NORBA,NORBB,DL (1,1,3),(MAXL+1)**2)
      CALL PSSX94(NORBB,NORBA,DLT(1,1,1),(MAXL+1)**2)
      CALL PSSX94(NORBB,NORBA,DLT(1,1,2),(MAXL+1)**2)
      CALL PSSX94(NORBB,NORBA,DLT(1,1,3),(MAXL+1)**2)
C
      RETURN
      END
C
      SUBROUTINE PSNA2D(IA,IB,DX,DY,DZ,DD)
C
C   Compute two-center \betan_{\mu\nu} <\mu|r^\nu|\nu> integrals for HA0
C   matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Index of the first atom
C      IB     - Index of the second atom
C      DX,DY,DZ
C             - Relative position of the atom B (a.u.s).
C      DD     - Output \beta <\mu|r^\nu|\nu> integrals
C
C   Accessed common blocks:
C
C      ATOMS  - Orbital numbers and atomic types
C      PSNPAR - NMR-related orbital exponents and beta parameters
C
C   Local storage:
C
C      3*(MAXL+1)**4 + 2*(MAXL+1)**2 + 6*(MAXL+1) DOUBLE PRECISION
C      cell, i.e. 423 doublewords for MAXL=2.
C
C   Module logic:
C
C   Bugs:
C
C      Orbital moments are limited to 2 by arbitrary dimentioned
C      temporaries, while the code (except the fragment extracting
C      atomic parameters from MNDO94's COMMON blocks)can in 
C      principle handle arbitrary L's.
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (HALFAU=0.183755972069092245497978684D-1)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION DD((MAXL+1)**2,(MAXL+1)**2,3)
C
      DIMENSION DDT(2*MAXL+1,2*MAXL+1,3), DDA(2*MAXL+1,2*MAXL+1,3)
C
C   Gather necessary parameters
C
      NORBA = NLAST(IA)-NFIRST(IA)+1
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPA = NAT(IA)
      ITYPB = NAT(IB)
C
C   Compute integrals
C
      DO 800 LA=0,MAXLA(ITYPA)
          IBASA = LA**2
          DO 700 LB=0,MAXLA(ITYPB)
              IBASB = LB**2
              DO 200 IX=1,3
                  DO 190 MB=1,2*LB+1
                      DO 180 MA=1,2*LA+1
                          DDA(MA,MB,IX) = ZERO
  180                 CONTINUE
  190             CONTINUE
  200         CONTINUE
              DO 600 ILA=1,MAXN(LA,ITYPA)
                  NA = NN(ILA,LA,ITYPA)
                  ZA = ZN(ILA,LA,ITYPA)
                  CA = CN(ILA,LA,ITYPA)
                  DO 500 ILB=1,MAXN(LB,ITYPB)
                      NB = NN(ILB,LB,ITYPB)
                      ZB = ZN(ILB,LB,ITYPB)
                      CB = CN(ILB,LB,ITYPB)
                      CALL PSODX(NA,LA,ZA,NB,LB,ZB,DX,DY,DZ,DDT,
     .                                             2*MAXL+1,2*MAXL+1)
                      DO 400 IX=1,3
                          DO 390 MB=1,2*LB+1
                              DO 380 MA=1,2*LA+1
                                  DDA(MA,MB,IX) = DDA(MA,MB,IX) 
     .                                          + CA*CB*DDT(MA,MB,IX)
  380                         CONTINUE
  390                     CONTINUE
  400                 CONTINUE
  500             CONTINUE
  600         CONTINUE
              SCN = HALFAU*(BETAN(LA,ITYPA) + BETAN(LB,ITYPB))
              DO 650 IX=1,3
                  DO 640 MB=1,2*LB+1
                      DO 630 MA=1,2*LA+1
                          DD(IBASA+MA,IBASB+MB,IX) = SCN*DDA(MA,MB,IX)
  630                 CONTINUE
  640             CONTINUE
  650         CONTINUE
  700     CONTINUE
  800 CONTINUE
C
C    Reordering of the integrals
C
      CALL PSSX94(NORBA,NORBB,DD(1,1,1),(MAXL+1)**2)
      CALL PSSX94(NORBA,NORBB,DD(1,1,2),(MAXL+1)**2)
      CALL PSSX94(NORBA,NORBB,DD(1,1,3),(MAXL+1)**2)
C
      RETURN
      END
C
      SUBROUTINE PSNHA1(IA,DL)
C
C   Compute one-center one-electron integrals necessary for the HA0
C   matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Index of the first atom
C      DL     - Output <\mu|L|\nu> integrals
C
C   Accessed common blocks:
C
C      ATOMS  - Orbital numbers and atomic types
C
C   Local storage:
C
C      (MAXL+1)**4 DOUBLE PRECISION cells
C
C   Module logic:
C
C      This is a degenerate version of PSNHA2, and shares all its
C      logic and shortcomings.
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (ONE=1.D0)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
C
      DIMENSION DL((MAXL+1)**2,(MAXL+1)**2,3)
      DIMENSION S00((MAXL+1)**2,(MAXL+1)**2)
C
      NORB = NLAST(IA)-NFIRST(IA)+1
C
C   Compute overlap integrals (i.e., initialize identity matrix ;-)
C   Since we make sure our integrals are normalized, this is Ok
C   for multiple-zeta case as well.
C
      DO 90 M2=1,NORB
          DO 70 M1=1,M2-1
              S00(M1,M2) = ZERO
   70     CONTINUE
          S00(M2,M2) = ONE
          DO 80 M1=M2+1,NORB
              S00(M1,M2) = ZERO
   80     CONTINUE
   90 CONTINUE
C
C   Compute integrals
C
      L     = 0
      IBAS  = 1
  100 CONTINUE
          CALL PSOLX(NORB,L,S00(1,IBAS),(MAXL+1)**2,DL(1,IBAS,1),
     .          (MAXL+1)**2,(MAXL+1)**2)
      IBAS  = IBAS  + 2*L+1
      L     = L     + 1
      IF( IBAS.LT.NORB ) GOTO 100
C
      CALL PSSX94(NORB,NORB,DL(1,1,1),(MAXL+1)**2)
      CALL PSSX94(NORB,NORB,DL(1,1,2),(MAXL+1)**2)
      CALL PSSX94(NORB,NORB,DL(1,1,3),(MAXL+1)**2)
C
      RETURN
      END
C
      SUBROUTINE PSNHAP(IA,HAB,H0B,LDH1,LDH2)
C
C   Compute upper triangle of each of the nine mixed one-electron
C   hamiltonian derivatives (paramagnetic part). Derivatives are in 
C   units c^-2, where c is speed of light.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Atom posessing the magnetic moment
C               IA>NATOM means NICS.
C      HAB    - Output HAB matrices (upper triangle filled)
C               First and second indices are AOs.
C               Third index is a magnetic field component (X,Y,Z)
C               Forth index is a magnetic dipole component (X,Y,Z)
C      H0B    - Input H0B matrices, as computed by PSNH0B
C               (i.e., with lower triangle being antisymmetric
C               version of the upper one).
C      LDH1,
C      LDH2   - Leading dimensions of HAB and H0B matrices.
C
C   Accessed common blocks:
C
C      ATOMS  - Range of orbitals carried by an atom
C      ATOMC  - Coordinates of atoms (don't you wish it was in a.u.s, boy?)
C      PSDOPT - Computation options, in particular NMRLEV.
C      PSDGBL - Number of atoms, for consistency with derivatives code.
C      PSPRT  - Printing unit.
C      PSPRTF - Debugging output controls.
C
C   Local storage:
C
C      Nope.
C
C   Module logic:
C
C   Bugs:
C
C      Selection of the atomic blocks at the initialization step is less
C      efficients than it might have been.
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      CHARACTER*1 VARN
      PARAMETER (BOHR=0.529167D0)
      PARAMETER (HALF=0.5D0)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ATOMC / COORD(3,LM1)
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
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSDGBL/, /PSPRTF/, /PSDOPT/, /PSPRT /
      DIMENSION VARN(3), RR(3)
      DIMENSION HAB(LDH1,LDH2,3,3), H0B(LDH1,LDH2,3)
      DATA VARN/'X','Y','Z'/
C
C    1                       L
C    - [R_\mu x R_\nu] <\mu|---|\nu> contribution is completely
C    2                      r^3
C     determined by H0B matrix, so that we can evaluate it first
C     (and initialize HAB at the same time as an added bonus)
C
      DO 900 IB=1,NATOM
          BX    = COORD(1,IB)
          BY    = COORD(2,IB)
          BZ    = COORD(3,IB)
          DO 890 IC=1,NATOM
              IF( NMRLEV.EQ.1 .AND. (IB.NE.IA .OR.  IC.NE.IA) .OR.
     .            NMRLEV.EQ.2 .AND. 
     .              (IB.NE.IA .AND. IC.NE.IA .AND. IB.NE.IC) ) GOTO 889
              CX    = COORD(1,IC)
              CY    = COORD(2,IC)
              CZ    = COORD(3,IC)
C            Actually, the expression is (we are omitting speed of light
C            on both sides, since it is taken care of elsewhere):
C                a,b                           0,b
C            H(P)      = - (1/2) (R  x R  )  H
C                mu,nu             mu   nu a   mu,nu
C
              RR(1) = (CY*BZ - CZ*BY)/(BOHR**2)
              RR(2) = (CZ*BX - CX*BZ)/(BOHR**2)
              RR(3) = (CX*BY - CY*BX)/(BOHR**2)
              DO 600 IHC=NFIRST(IC),NLAST(IC)
                  DO 590 IHB=NFIRST(IB),NLAST(IB)
                      HAB(IHB,IHC,1,1) = HALF*RR(1)*H0B(IHB,IHC,1)
                      HAB(IHB,IHC,1,2) = HALF*RR(1)*H0B(IHB,IHC,2)
                      HAB(IHB,IHC,1,3) = HALF*RR(1)*H0B(IHB,IHC,3)
                      HAB(IHB,IHC,2,1) = HALF*RR(2)*H0B(IHB,IHC,1)
                      HAB(IHB,IHC,2,2) = HALF*RR(2)*H0B(IHB,IHC,2)
                      HAB(IHB,IHC,2,3) = HALF*RR(2)*H0B(IHB,IHC,3)
                      HAB(IHB,IHC,3,1) = HALF*RR(3)*H0B(IHB,IHC,1)
                      HAB(IHB,IHC,3,2) = HALF*RR(3)*H0B(IHB,IHC,2)
                      HAB(IHB,IHC,3,3) = HALF*RR(3)*H0B(IHB,IHC,3)
  590             CONTINUE
  600         CONTINUE
  889         CONTINUE
  890     CONTINUE
  900 CONTINUE
C
C    Symmetrize result (this is not necessary, since HAB(P) defined
C    in this way is guaranteed to be symmetric as long as H0B is
C    antisymmetric, but who cares!).
C
      DO 3500 IY=1,3
          DO 3400 IX=1,3
              CALL PSNMSY(IA,'HAB(P)',VARN(IX),VARN(IY),
     .                    HAB(1,1,IX,IY),LDH1)
 3400     CONTINUE
 3500 CONTINUE
C
      IF(LPNMR) THEN
          DO 5000 IX=1,3
              DO 4900 IY=1,3
                  IF(IA.LT.NATOM) THEN
                      WRITE(NB6,11200) IA, VARN(IX), VARN(IY)
                  ELSE
                      WRITE(NB6,11202) IA-NATOM, VARN(IX), VARN(IY)
                  ENDIF
                  CALL PSDPAU(NORBS,HAB(1,1,IX,IY),LDH1)
 4900         CONTINUE
 5000     CONTINUE
      ENDIF
      RETURN
11200 FORMAT(' SYMMETRIZED HAB(P) MATRIX FOR ATOM ',I4,' (COMPONENT ',
     .       A,A,') IS:')
11202 FORMAT(' SYMMETRIZED HAB(P) MATRIX FOR NICS CENTER ',I4,
     .       ' (COMPONENT ',A,A,') IS:')
      END
C
      SUBROUTINE PSNHAD(IA,HAB,LDH1,LDH2,RO,LDR)
C
C   Compute upper triangle of each of the nine diamagnetic parts of
C   the mixed one-electron hamiltonian derivatives. Derivatives are 
C   in units c^-2, where c is speed of light.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Atom posessing the magnetic moment
C               IA>NATOM means NICS.
C      HAB    - Output HAB matrices (upper triangle filled)
C               First and second indices are AOs.
C               Third index is a magnetic field component (X,Y,Z)
C               Forth index is a magnetic dipole component (X,Y,Z)
C      LDH1,
C      LDH2   - Leading dimensions of HAB matrices.
C      RO     - Alpha density matrix.
C      LDR    - Leading dimension of ROALP.
C
C   Accessed common blocks:
C
C      ATOMS  - Range of orbitals carried by an atom
C      ATOMC  - Coordinates of atoms (don't you wish it was in a.u.s, boy?)
C      PSDOPT - Computation options, in particular NMRLEV.
C      PSDGBL - Number of atoms, for consistency with derivatives code.
C      PSPRT  - Printing unit.
C      PSPRTF - Debugging output controls.
C
C   Modified common blocks:
C
C      PSO3PR - Three-center integrals precision control. (EPSR is used here).
C
C   Local storage:
C
C      Some 9*(MAXL+1)**4 (ca. 730) DOUBLE PRECISION words.
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, MXNICS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI, VARINT
      LOGICAL PSNINC
      CHARACTER*1 VARN
      PARAMETER (MAXL=2)
      PARAMETER (BOHR=0.529167D0)
      PARAMETER (HALF=0.5D0)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ATOMC / COORD(3,LM1)
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
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
     ./PSNICS/ COORDN(3,MXNICS), NNICS
     ./PSO3PR/ EPSL, EPSR
      SAVE /PSDGBL/, /PSPRTF/, /PSDOPT/, /PSNICS/, /PSO3PR/, /PSPRT /
      EXTERNAL PSNDMX
      DIMENSION VARN(3), RR(3), RB(3), RC(3)
      DIMENSION HAB(LDH1,LDH2,3,3), RO(LDR)
      DIMENSION DD((MAXL+1)**2,3,3,(MAXL+1)**2)
      DIMENSION DL((MAXL+1)**2,3,3,(MAXL+1)**2)
      DATA VARN/'X','Y','Z'/
C
      VARINT = MOD(INTCTL,200).GE.100
      IGO    = MAX(MOD(INTCTL,20)-10,0)
C
      DO 100 IY=1,3
          DO 90 IX=1,3
              CALL PSNZR(IA,HAB(1,1,IX,IY),LDH1)
   90     CONTINUE
  100 CONTINUE
C
C    One-center contributions, centered at the IA atom.
C    NICS does not have one-center components.
C
      IF( IA.LE.NATOM ) THEN
          XA = COORD(1,IA)
          YA = COORD(2,IA)
          ZA = COORD(3,IA)
          IBASA = NFIRST(IA)-1
          NORBA = NLAST(IA)-NFIRST(IA)+1
          CALL PSHAB1(IA,DD)
          DO 1200 M2=1,NORBA
              IH2 = IBASA + M2
              DO 1190 M1=1,NORBA
                  IH1 = IBASA + M1
                  DO 1180 IY=1,3
                      DO 1170 IX=1,3
                          HAB(IH1,IH2,IX,IY) = 
     .                        HAB(IH1,IH2,IX,IY) + HALF*DD(M1,IX,IY,M2)
 1170                 CONTINUE
 1180             CONTINUE
 1190         CONTINUE
 1200     CONTINUE
      ELSE
          XA = COORDN(1,IA-NATOM)
          YA = COORDN(2,IA-NATOM)
          ZA = COORDN(3,IA-NATOM)
      ENDIF
C
C    Two-center contributions
C
      IF( NMRLEV.GT.1 ) THEN
          DO 1600 IB=1,NATOM
              IF( IB.EQ.IA ) GOTO 1599
              IBASB  = NFIRST(IB)-1
              NORBB  = NLAST(IB)-NFIRST(IB)+1
              RR(1)  = (XA-COORD(1,IB))/BOHR
              RR(2)  = (YA-COORD(2,IB))/BOHR
              RR(3)  = (ZA-COORD(3,IB))/BOHR
C
C            <B|A|A> contributions. (and <A|A|B> ones if INDSYM is 2).
C            PSAB2P is, strictly speaking, a paramagnetic contribution,
C            but it is needed here to make all matrices Hermitian.
C            Neither is present for NICS
C
              IF( IA.LE.NATOM ) THEN
                  CALL PSAB2P(IA,IB,RR(1),RR(2),RR(3),DL)
                  CALL PSAB2D(IA,IB,RR(1),RR(2),RR(3),DD)
                  CALL PSAB2M(NORBB,NORBA,RR,DL)
                  DO 1300 M2=1,NORBA
                      IH2 = IBASA + M2
                      DO 1290 M1=1,NORBB
                          IH1 = IBASB + M1
                          DO 1280 IY=1,3
                            DO 1270 IX=1,3
                              VL=HALF*(DD(M1,IX,IY,M2)+DL(M1,IX,IY,M2))
                              HAB(IH1,IH2,IX,IY)=HAB(IH1,IH2,IX,IY)+VL
                              IF( INDSYM.EQ.2 )
     .                        HAB(IH2,IH1,IX,IY)=HAB(IH2,IH1,IX,IY)+VL
 1270                       CONTINUE
 1280                     CONTINUE
 1290                 CONTINUE
 1300             CONTINUE
                  IF( INDSYM.NE.2 ) THEN
C
C            <A|A|B> contributions. The total contribution should be symmetric 
C            with <B|A|A> contribution, but the individual integrals need not,
C            so that we can use it as one more correctness check.
C
                      CALL PSAB2Q(IA,IB,RR(1),RR(2),RR(3),DL)
                      CALL PSAB2E(IA,IB,RR(1),RR(2),RR(3),DD)
C                    Note that PSAB2M wants to have right-minus-left distance,
C                    while we are feeding it with left-minus-right. Therefore
C                    all DL terms will have to be negated before the incorporation
C                    into HAB.
                      CALL PSAB2M(NORBA,NORBB,RR,DL)
                      DO 1400 M2=1,NORBB
                          IH2 = IBASB + M2
                          DO 1390 M1=1,NORBA
                              IH1 = IBASA + M1
                              DO 1380 IY=1,3
                                  DO 1370 IX=1,3
                                      VL = HALF*(DD(M1,IX,IY,M2)
     .                                          -DL(M1,IX,IY,M2))
                                      HAB(IH1,IH2,IX,IY) = 
     .                                      HAB(IH1,IH2,IX,IY) + VL
 1370                             CONTINUE
 1380                         CONTINUE
 1390                     CONTINUE
 1400                 CONTINUE
                  ENDIF
              ENDIF
*            End of the part excluded for NICS.
C
C            <B|A|B> contributions. NICS do have these contributions.
C
              CALL PSAB2B(IB,RR(1),RR(2),RR(3),DD)

              DO 1500 M2=1,NORBB
                  IH2 = IBASB + M2
                  DO 1490 M1=1,NORBB
                      IH1 = IBASB + M1
                      DO 1480 IY=1,3
                          DO 1470 IX=1,3
                              VL = HALF*DD(M1,IX,IY,M2)
                              HAB(IH1,IH2,IX,IY)=HAB(IH1,IH2,IX,IY)+VL
 1470                     CONTINUE
 1480                 CONTINUE
 1490             CONTINUE
 1500         CONTINUE
 1599         CONTINUE
 1600     CONTINUE
      ENDIF
C
C    Three-center contributions. All of those are present for NICS.
C
      IF( NMRLEV.GT.2 ) THEN
          DO 2700 IB=1,NATOM
              IF( IB.EQ.IA ) GOTO 2699
              IBASB = NFIRST(IB)-1
              NORBB = NLAST(IB)-NFIRST(IB)+1
              RB(1) = (COORD(1,IB)-XA)/BOHR
              RB(2) = (COORD(2,IB)-YA)/BOHR
              RB(3) = (COORD(3,IB)-ZA)/BOHR
              ICMAX = NATOM
              IF( INDSYM.EQ.2 ) ICMAX = IB-1
              DO 2600 IC=1,ICMAX
                  IF( IC.EQ.IB .OR. IC.EQ.IA ) GOTO 2599
                  IF( .NOT.PSNINC(IB,IC) ) GOTO 2599
                  IBASC = NFIRST(IC)-1
                  NORBC = NLAST(IC)-NFIRST(IC)+1
                  RC(1) = (COORD(1,IC)-XA)/BOHR
                  RC(2) = (COORD(2,IC)-YA)/BOHR
                  RC(3) = (COORD(3,IC)-ZA)/BOHR
                  RR(1) = RC(1) - RB(1)
                  RR(2) = RC(2) - RB(2)
                  RR(3) = RC(3) - RB(3)
C
C                Variable integrals precision control
C
                  IF(VARINT) THEN
                      ROMAX = PSNDMX(RO,IBASB+1,NORBB,IBASC+1,NORBC)
                      ROMAX = 5D0*MAX(ROMAX,DNCOFF)
                      EPSR  = MAX(1D-6,DNCOFF/ROMAX)
                      EPSL  = MAX(2D-6,DNCOFF/ROMAX)
                  ENDIF
C
C                <B|A|C> contributions. PSAB3P is, strictly speaking,
C                a paramagnetic contribution, but it have to be here
C                to make HAB(D) Hermitian.
C
                  IF( IGO.EQ.0 ) THEN
                      CALL PSAB3D(IB,IC,RB,RC,DD)
                      CALL PSAB3P(IB,IC,RB,RC,DL)
                  ELSE
                      CALL PSA3DG(IGO,IB,IC,RB,RC,DD)
                      CALL PSA3PG(IGO,IB,IC,RB,RC,DL)
                  ENDIF
                  CALL PSAB2M(NORBB,NORBC,RR,DL)
C
                  DO 2500 M2=1,NORBC
                      IH2 = IBASC + M2
                      DO 2490 M1=1,NORBB
                          IH1 = IBASB + M1
                          DO 2480 IY=1,3
                              DO 2470 IX=1,3
                                  VL = HALF*(DD(M1,IX,IY,M2)+
     .                                       DL(M1,IX,IY,M2))
                                  HAB(IH1,IH2,IX,IY) = 
     .                                HAB(IH1,IH2,IX,IY) + VL
                                  IF( INDSYM.EQ.2 ) HAB(IH2,IH1,IX,IY) = 
     .                                HAB(IH2,IH1,IX,IY) + VL
 2470                         CONTINUE
 2480                     CONTINUE
 2490                 CONTINUE
 2500             CONTINUE
 2599             CONTINUE
 2600         CONTINUE
 2699         CONTINUE
 2700     CONTINUE
      ENDIF
      IF(LPNMR) THEN
          DO 3000 IX=1,3
              DO 2900 IY=1,3
                  IF( IA.LE.NATOM ) THEN
                      WRITE(NB6,11198) IA, VARN(IX), VARN(IY)
                  ELSE
                      WRITE(NB6,11199) IA-NATOM, VARN(IX), VARN(IY)
                  ENDIF
                  CALL PSDPGM(NORBS,NORBS,HAB(1,1,IX,IY),LDH1)
 2900         CONTINUE
 3000     CONTINUE
      ENDIF
C
C    Symmetrize result
C
      DO 3500 IY=1,3
          DO 3400 IX=1,3
              CALL PSNMSY(IA,'HAB(D)',VARN(IX),VARN(IY),
     .                    HAB(1,1,IX,IY),LDH1)
 3400     CONTINUE
 3500 CONTINUE
C
      IF(LPNMR) THEN
          DO 5000 IX=1,3
              DO 4900 IY=1,3
                  IF( IA.LE.NATOM ) THEN
                      WRITE(NB6,11200) IA, VARN(IX), VARN(IY)
                  ELSE
                      WRITE(NB6,11201) IA-NATOM, VARN(IX), VARN(IY)
                  ENDIF
                  CALL PSDPAU(NORBS,HAB(1,1,IX,IY),LDH1)
 4900         CONTINUE
 5000     CONTINUE
      ENDIF
      RETURN
11198 FORMAT(' RAW HAB(D) MATRIX FOR ATOM ',I4,' (COMPONENT ',A,A,
     .       ') IS:')
11199 FORMAT(' RAW HAB(D) MATRIX FOR NICS POINT ',I4,' (COMPONENT ',A,A,
     .       ') IS:')
11200 FORMAT(' SYMMETRIZED HAB(D) MATRIX FOR ATOM ',I4,
     .       ' (COMPONENT ',A,A,') IS:')
11201 FORMAT(' SYMMETRIZED HAB(D) MATRIX FOR NICS POINT ',I4,
     .       ' (COMPONENT ',A,A,') IS:')
      END
C
      SUBROUTINE PSAB2M(NORB1,NORB2,R12,DL)
C
C   Evaluate [<\mu r|r^-3 L_b|\nu>,R_\nu - R_\mu]_a products.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NORB1  - Number of basis functions on the left-hand atom
C      NORB2  - Number of basis functions on the right-hand atom
C      R12    - Position of the right atom relative to the left atom
C      DL     - Input: <\mu r_\mu|r^-3 L|\nu> integrals, four-index array.
C                   First index: Orbitals at the left center.
C                   Second index: Projections of the r_\mu operator.
C                   Third index: Projections of the L operator.
C                   Forth index: Orbitals at the right center.
C               Output: <\mu [r_\mu,R_{12}]|r^-3 L\nu> integrals.
C                   First index: Orbitals at the left center.
C                   Second index: Projections of the [r_\mu,R_{12}] operator.
C                   Third index: Projections of the L operator.
C                   Forth index: Orbitals at the right center.
C               All operator projections are in the x,y,z order.
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
C
      DIMENSION R12(3)
      DIMENSION DL((MAXL+1)**2,3,3,(MAXL+1)**2)
C
      DO 500 M2=1,NORB2
          DO 400 IY=1,3
              DO 300 M1=1,NORB1
                  VX = DL(M1,2,IY,M2)*R12(3) - DL(M1,3,IY,M2)*R12(2)
                  VY = DL(M1,3,IY,M2)*R12(1) - DL(M1,1,IY,M2)*R12(3)
                  VZ = DL(M1,1,IY,M2)*R12(2) - DL(M1,2,IY,M2)*R12(1)
                  DL(M1,1,IY,M2) = VX
                  DL(M1,2,IY,M2) = VY
                  DL(M1,3,IY,M2) = VZ
  300         CONTINUE
  400     CONTINUE
  500 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSHAB1(IA,DD)
C
C   Compute one-center one-electron integrals necessary for the HAB
C   matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Index of the first atom
C      DD     - Output <\mu|r^-3 (\delta_{ab} r^2 - r_b r_a|\nu> integrals
C
C   Accessed common blocks:
C
C      ATOMS  - Atoms identities and range of orbitals
C      PSNPAR - NMR parameters
C
C   Local storage:
C
C      Some 9*((MAXL+2)**2+2)*(MAXL+1)**2 (ca. 1460) DOUBLE PRECISION
C      cells.
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION DD((MAXL+1)**2,3,3,(MAXL+1)**2)
C
      DIMENSION AX (2*MAXL+1,2*(MAXL+2)+1,3)
      DIMENSION DDT(2*MAXL+1,3,3,2*MAXL+1)
      EXTERNAL PSR31C
C
      NORB = NLAST(IA)-NFIRST(IA)+1
      ITYP = NAT(IA)
C
      DO 1000 L1=0,MAXLA(ITYP)
          IBAS1 = L1**2
          DO 900 L2=0,MAXLA(ITYP)
              IBAS2 = L2**2
              DO 400 M2=1,2*L2+1
                  DO 300 IY=1,3
                      DO 200 IX=1,3
                          DO 100 M1=1,2*L1+1
                              DD(IBAS1+M1,IX,IY,IBAS2+M2) = ZERO
  100                     CONTINUE
  200                 CONTINUE
  300             CONTINUE
  400         CONTINUE
              DO 800 IL2=1,MAXN(L2,ITYP)
                  N2  = NN(IL2,L2,ITYP)
                  Z2  = ZN(IL2,L2,ITYP)
                  C2  = CN(IL2,L2,ITYP)
C
C                Evaluate auxiliary integral. Left-hand function is contracted
C                at this point, while the right-hand one is still in the
C                elementary form.
C
                  VAL = ZERO
                  DO 500 IL1=1,MAXN(L1,ITYP)
                      N1  = NN(IL1,L1,ITYP)
                      Z1  = ZN(IL1,L1,ITYP)
                      C1  = CN(IL1,L1,ITYP)
                      VAL = VAL + C1*PSR31C(N1,Z1,N2+2,Z2)
  500             CONTINUE
C
C                Do angular part of the integration. Only integrals with equal
C                L and M on both sides can survive.
C
                  DO 540 M1=1,2*L1+1
                      DO 510 M2=1,2*(L2-2)+1
                          AX(M1,M2,1) = ZERO
  510                 CONTINUE
                      DO 520 M2=1,2*(L2+0)+1
                          AX(M1,M2,2) = ZERO
  520                 CONTINUE
                      DO 530 M2=1,2*(L2+2)+1
                          AX(M1,M2,3) = ZERO
  530                 CONTINUE
  540             CONTINUE
                  IF( L2-2.EQ.L1 ) THEN
                      DO 550 M=1,2*L1+1
                          AX(M,M,1) = VAL
  550                 CONTINUE
                  ENDIF
                  IF( L2+0.EQ.L1 ) THEN
                      DO 560 M=1,2*L1+1
                          AX(M,M,2) = VAL
  560                 CONTINUE
                  ENDIF
                  IF( L2+2.EQ.L1 ) THEN
                      DO 570 M=1,2*L1+1
                          AX(M,M,3) = VAL
  570                 CONTINUE
                  ENDIF
C
C                Evaluate quadrupole moment integrals and contract elementary
C                integrals by the second index.
C
                  CALL PSOQS(2*L1+1,N2,L2,Z2,AX,2*MAXL+1,2*(MAXL+2)+1,
     .                                          DDT,2*MAXL+1)
                  DO 700 M2=1,2*L2+1
                      DO 690 IY=1,3
                          DO 680 IX=1,3
                              DO 670 M1=1,2*L1+1
                                  DD(IBAS1+M1,IX,IY,IBAS2+M2) = 
     .                                DD(IBAS1+M1,IX,IY,IBAS2+M2) + 
     .                                C2*DDT(M1,IX,IY,M2)
  670                         CONTINUE
  680                     CONTINUE
  690                 CONTINUE
  700             CONTINUE
  800         CONTINUE
  900     CONTINUE
 1000 CONTINUE
C
      DO 1900 IY=1,3
          DO 1800 IX=1,3
              CALL PSSX94(NORB,NORB,DD(1,IX,IY,1),9*((MAXL+1)**2))
 1800     CONTINUE
 1900 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSAB2P(IA,IB,RX,RY,RZ,DL)
C
C   Compute two-center dipole-orbital moment contributions of the <B|A|A> 
C   type to the HAB matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Index of the right-hand atom (one carrying operator)
C      IB     - Index of the left-hand atom
C      RX,RY,RZ
C             - Relative coordinates of the atom A (right-hand one).
C               Coordinates are in Bohrs.
C      DL     - Output <r_\mu \mu|r^-3 L^\nu|\nu> integrals. Second
C               index runs over dipole moment projections, while third
C               index denotes orbital moment components.
C
C   Accessed common blocks:
C
C      ATOMS  - Atoms identities and range of orbitals
C      PSNPAR - NMR parameters
C
C   Local storage:
C
C      Some 28 (MAXL+1)**4 (ca. 2300) DOUBLE PRECISION cells.
C
C   Module logic:
C
C      is sorta missin', m'am.
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION DL((MAXL+1)**2,3,3,(MAXL+1)**2)
C
      DIMENSION RR3 (2*MAXL+1,2*MAXL+1,3)
      DIMENSION UR3 (2*MAXL+1,2*(MAXL+1)+1)
      DIMENSION DR3 (2*MAXL+1,2*(MAXL-1)+1)
      DIMENSION TR3 (2*(MAXL+1)+1,2*MAXL+1)
      DIMENSION TRR3(2*MAXL+1,2*MAXL+1,3)
      DIMENSION RR3L(2*MAXL+1,2*MAXL+1,3,3)
C
      NORBA = NLAST(IA)-NFIRST(IA)+1
      ITYPA = NAT(IA)
C
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPB = NAT(IB)
C
C    Zero out terms with LA = 0
C
      DO 20 IY=1,3
          DO 15 IX=1,3
              DO 10 MB=1,NORBB
                  DL(MB,IX,IY,1) = ZERO
   10         CONTINUE
   15     CONTINUE
   20 CONTINUE
C 
      DO 800 LB=0,MAXLA(ITYPB)
          IBASB = LB**2
          DO 790 LA=1,MAXLA(ITYPA)
              IBASA = LA**2
C
C            Evaluate all B-dipole auxiliary integrals over contractions
C            at B and A. Auxiliary integrals are in RR3 in the right order
C            (i.e. B at the left) at the end of this loop.
C
              DO 50 IX=1,3
                  DO 40 MB=1,2*LB+1
                      DO 30 MA=1,2*LA+1
                          RR3(MB,MA,IX) = ZERO
   30                 CONTINUE
   40             CONTINUE
   50         CONTINUE
              DO 600 ILB=1,MAXN(LB,ITYPB)
                  NB  = NN(ILB,LB,ITYPB)
                  ZB  = ZN(ILB,LB,ITYPB)
                  SCB = CN(ILB,LB,ITYPB)
C
C                Evaluate all auxiliary integrals on contracted basis function
C                at A. Transposed integrals are collected in UR3 and (if LB is
C                at least 1) in DR3.
C
                  DO 100 MB=1,2*(LB+1)+1
                      DO 90 MA=1,2*LA+1
                          UR3(MA,MB) = ZERO
   90                 CONTINUE
  100             CONTINUE
                  IF( LB.GE.1 ) THEN
                      DO 150 MB=1,2*(LB-1)+1
                          DO 140 MA=1,2*LA+1
                              DR3(MA,MB) = ZERO
  140                     CONTINUE
  150                 CONTINUE
                  ENDIF
                  DO 400 ILA=1,MAXN(LA,ITYPA)
                      NA  = NN(ILA,LA,ITYPA)
                      ZA  = ZN(ILA,LA,ITYPA)
                      SCA = CN(ILA,LA,ITYPA)
                      CALL PSORX(NB+1,LB+1,ZB,NA,LA,ZA,RX,RY,RZ,3,
     .                           TR3,2*(MAXL+1)+1)
                      DO 200 MB=1,2*(LB+1)+1
                          DO 190 MA=1,2*LA+1
                              UR3(MA,MB) = UR3(MA,MB) + SCA*TR3(MB,MA)
  190                     CONTINUE
  200                 CONTINUE
                      IF( LB.GE.1 ) THEN
                          CALL PSORX(NB+1,LB-1,ZB,NA,LA,ZA,RX,RY,RZ,3,
     .                               TR3,2*(MAXL+1)+1)
                          DO 300 MB=1,2*(LB-1)+1
                              DO 290 MA=1,2*LA+1
                                  DR3(MA,MB) = DR3(MA,MB) 
     .                                         + SCA*TR3(MB,MA)
  290                         CONTINUE
  300                     CONTINUE
                      ENDIF
  400             CONTINUE
C
C                Make (B-)dipole moment integrals out of auxiliary integrals
C                for each primitive Slater AO in contraction at center B and
C                accumulate contributions transposing them once again.
C
                  CALL PSODS(2*LA+1,NB,LB,ZB,DR3,2*MAXL+1,UR3,2*MAXL+1,
     .                       TRR3,2*MAXL+1,2*MAXL+1)
                  DO 500 IX=1,3
                      DO 490 MB=1,2*LB+1
                          DO 480 MA=1,2*LA+1
                              RR3(MB,MA,IX) = RR3(MB,MA,IX) 
     .                                      + SCB*TRR3(MA,MB,IX)
  480                     CONTINUE
  490                 CONTINUE
  500             CONTINUE
  600         CONTINUE
C
C            Apply orbital moment operator to the auxiliary B-dipole integrals.
C            All contractions are already carried out at this point.
C
              CALL PSOLX(2*LB+1,LA,RR3(1,1,1),2*MAXL+1,RR3L(1,1,1,1),
     .                   2*MAXL+1,3*(2*MAXL+1))
              CALL PSOLX(2*LB+1,LA,RR3(1,1,2),2*MAXL+1,RR3L(1,1,2,1),
     .                   2*MAXL+1,3*(2*MAXL+1))
              CALL PSOLX(2*LB+1,LA,RR3(1,1,3),2*MAXL+1,RR3L(1,1,3,1),
     .                   2*MAXL+1,3*(2*MAXL+1))
C
C            Stuff computed integrals into output array, and we have finished
C            with this (LB,LA) pair.
C
              DO 700 MA=1,2*LA+1
                  DO 690 IY=1,3
                      DO 680 IX=1,3
                          DO 670 MB=1,2*LB+1
                              DL(IBASB+MB,IX,IY,IBASA+MA) = 
     .                            RR3L(MB,MA,IX,IY)
  670                     CONTINUE
  680                 CONTINUE
  690             CONTINUE
  700         CONTINUE
C
  790     CONTINUE
  800 CONTINUE
C
C        Put all integrals into the MNDO-94 order.
C
      DO 900 IY=1,3
          DO 890 IX=1,3
              CALL PSSX94(NORBB,NORBA,DL(1,IX,IY,1),9*((MAXL+1)**2))
  890     CONTINUE
  900 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSAB2Q(IA,IB,RX,RY,RZ,DL)
C
C   Compute two-center dipole-orbital moment contributions of the <A|A|B> 
C   type to the HAB matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C                                                       ^
C      IA     - Index of the left-hand atom (carrying x L r^-3 operator)
C      IB     - Index of the right-hand atom
C      RX,RY,RZ
C             - Relative coordinates of the atom A (left-hand one) in
C               atomic units.
C                                 ^A
C                          A      L       B
C      DL     - Output < mu  r | ---- | nu > integrals.
C                             A    3
C                                 r
C                                  A
C               First index:  enumerates left-hand orbitals.
C               Second index, 1 to 3: components of the r_A operator,
C                   in the X,Y,Z order.
C               Third index: 1 to 3: components of the L^A operator,
C                   in the X,Y,Z order.
C               Forth index: enumerated right-hand orbitals.
C
C   Accessed common blocks:
C
C      ATOMS  - Atoms identities and range of orbitals
C      PSNPAR - NMR parameters
C
C   Local storage:
C
C      9*(MAXL+1)**2 DOUBLE PRECISION storage units, which
C      translates to some 270 DP words. Not excessive...
C
C   Module logic:
C
C   Bugs:
C
C      I hope there are none...
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION DL ((MAXL+1)**2,3,3,(MAXL+1)**2)
      DIMENSION DLT(2*MAXL+1,3,3,2*MAXL+1)
C
      NORBA = NLAST(IA)-NFIRST(IA)+1
      ITYPA = NAT(IA)
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPB = NAT(IB)
C    Zero out the total integral.
      DO 400 MB=1,NORBB
          DO 300 IY=1,3
              DO 200 IX=1,3
                  DO 100 MA=1,NORBA
                      DL(MA,IX,IY,MB) = ZERO
  100             CONTINUE
  200         CONTINUE
  300     CONTINUE
  400 CONTINUE
C
      DO 3000 LA=0,MAXLA(ITYPA)
          IBASA = LA**2
          DO 2900 LB=0,MAXLA(ITYPB)
              IBASB = LB**2
              DO 2000 JA=1,MAXN(LA,ITYPA)
                  NA = NN(JA,LA,ITYPA)
                  ZA = ZN(JA,LA,ITYPA)
                  CA = CN(JA,LA,ITYPA)
                  DO 1000 JB=1,MAXN(LB,ITYPB)
                      NB = NN(JB,LB,ITYPB)
                      ZB = ZN(JB,LB,ITYPB)
                      CB = CN(JB,LB,ITYPB)
                      CALL PSRLA3(NA,LA,ZA,NB,LB,ZB,RX,RY,RZ,
     .                            DLT,2*MAXL+1)
                      DO 900 MB=1,2*LB+1
                          DO 800 IY=1,3
                              DO 700 IX=1,3
                                  DO 600 MA=1,2*LA+1
                                      DL(IBASA+MA,IX,IY,IBASB+MB) = 
     .                                      DL(IBASA+MA,IX,IY,IBASB+MB)
     .                                    + CA*CB*DLT(MA,IX,IY,MB)
  600                             CONTINUE
  700                         CONTINUE
  800                     CONTINUE
  900                 CONTINUE
 1000             CONTINUE
 2000         CONTINUE
 2900     CONTINUE
 3000 CONTINUE
C
C        Put all integrals into the MNDO-94 order.
C
      DO 3900 IY=1,3
          DO 3890 IX=1,3
              CALL PSSX94(NORBA,NORBB,DL(1,IX,IY,1),9*((MAXL+1)**2))
 3890     CONTINUE
 3900 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSAB3P(IB,IC,RB,RC,DL)
C
C   Compute two-center dipole-orbital moment contributions of the <B|A|C> 
C   type. (Paramagnetic part of the HAB matrix). The actual integrals
C   are:
C                  A
C                 L
C        B  B      b       C
C     <mu  r  | ------ | nu >
C           a       3
C                  A
C                 r
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IB     - Index of the left-hand atom
C      IC     - Index of the right-hand atom
C      RB     - Coordinates of the left-hand atom relative to the operator
C               center (atomic units).
C      RC     - Coordinates of the right-hand atom relative to the operator
C               center (atomic units).
C      DL     - Output <r_\mu \mu|r^-3 L^\nu|\nu> integrals. 
C               First index: Orbitals at the left-hand center (B) in the
C                   MNDO94 order.
C               Second index: Components of the radius-vector at B.
C               Third index: Projections of the orbital moment operator
C                   (a.k.a. magnetic moment components).
C               Forth index: Orbitals at the right-hand center (C) in the
C                   MNDO94 order.
C
C   Accessed common blocks:
C
C      ATOMS  - Atoms identities and range of orbitals
C      PSNPAR - NMR parameters
C
C   Local storage:
C
C      21*(2*MAXL+1)**2 DOUBLE PRECISION storage units. With MAXL=2,
C      it's some 530 DP words.
C
C   Module logic:
C
C      Integrals are evaluated using reduction to N+1,L+/-1 orbitals
C      on the left-hand center.
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION DL((MAXL+1)**2,3,3,(MAXL+1)**2)
      DIMENSION RB(3), RC(3)
C
      DIMENSION ATDN (2*MAXL+1,3,2*(MAXL-1)+1)
      DIMENSION ATDNT(2*(MAXL-1)+1,3,2*MAXL+1)
      DIMENSION ATUP (2*MAXL+1,3,2*(MAXL+1)+1)
      DIMENSION ATUPT(2*(MAXL+1)+1,3,2*MAXL+1)
      DIMENSION AT   (2*MAXL+1,3,2*MAXL+1,3)
C
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPB = NAT(IB)
      NORBC = NLAST(IC)-NFIRST(IC)+1
      ITYPC = NAT(IC)
C
      DO 200 MC=1,NORBC
          DO 190 IY=1,3
              DO 180 IX=1,3
                  DO 170 MB=1,NORBB
                      DL(MB,IX,IY,MC) = ZERO
  170             CONTINUE
  180         CONTINUE
  190     CONTINUE
  200 CONTINUE
C
      DO 2000 LB=0,MAXLA(ITYPB)
          IBB = LB**2
          DO 1900 JB=1,MAXN(LB,ITYPB)
              NB = NN(JB,LB,ITYPB)
              ZB = ZN(JB,LB,ITYPB)
              CB = CN(JB,LB,ITYPB)
              DO 1000 LC=0,MAXLA(ITYPC)
                  IBC = LC**2
C                ATDN and ATUP are indexed by orbital on C, projection 
C                of L and orbital on B. ATDN are for orbital moments
C                one less than LB. ATUP are for orbital moments one more
C                than LB.
                  DO 400 MB=1,2*(LB-1)+1
                      DO 395 J=1,3
                          DO 390 MC=1,2*LC+1
                              ATDN(MC,J,MB) = ZERO
  390                     CONTINUE
  395                 CONTINUE
  400             CONTINUE
                  DO 500 MB=1,2*(LB+1)+1
                      DO 495 J=1,3
                          DO 490 MC=1,2*LC+1
                              ATUP(MC,J,MB) = ZERO
  490                     CONTINUE
  495                 CONTINUE
  500             CONTINUE
C                It is Ok to accumulate contraction on C, since dipole
C                reduction depends only on the quantum number on B.
                  DO 900 JC=1,MAXN(LC,ITYPC)
                      NC = NN(JC,LC,ITYPC)
                      ZC = ZN(JC,LC,ITYPC)
                      CC = CN(JC,LC,ITYPC)
C                    B is on the left, C is on the right. We might have
C                    transposed the centers, but then we would have to
C                    change the sign of the integrals.
                      IF( LB.GT.0 )
     .                CALL PS3LR3(NB+1,LB-1,ZB,RB,NC,LC,ZC,RC,
     .                            ATDNT,2*(MAXL-1)+1)
                      CALL PS3LR3(NB+1,LB+1,ZB,RB,NC,LC,ZC,RC,
     .                            ATUPT,2*(MAXL+1)+1)
C                    Accumulate contractions in C, placing B on the
C                    right-hand side.
                      DO 600 MB=1,2*(LB-1)+1
                          DO 595 J=1,3
                              DO 590 MC=1,2*LC+1
                                  ATDN(MC,J,MB) = ATDN(MC,J,MB) 
     .                                          + CC*ATDNT(MB,J,MC)
  590                         CONTINUE
  595                     CONTINUE
  600                 CONTINUE
                      DO 700 MB=1,2*(LB+1)+1
                          DO 695 J=1,3
                              DO 690 MC=1,2*LC+1
                                  ATUP(MC,J,MB) = ATUP(MC,J,MB) 
     .                                          + CC*ATUPT(MB,J,MC)
  690                         CONTINUE
  695                     CONTINUE
  700                 CONTINUE
  900             CONTINUE
C                AT is indexed by orbital on C, projection of L, orbital
C                on B and projection of the radius-vector at B.
                  DO 950 J=1,3
                      CALL PSODS(2*LC+1,NB,LB,ZB,
     .                           ATDN(1,J,1),3*(2*MAXL+1),
     .                           ATUP(1,J,1),3*(2*MAXL+1),
     .                           AT(1,J,1,1),3*(2*MAXL+1),2*MAXL+1)
  950             CONTINUE
C                Finally, DL is indexed by orbital on B, projection of the
C                radius-vector on B, projection of L and orbital on C.
C                Note that DL contains all the orbitals, while AT has
C                orbitals with the same LB and LC only.
                  DO 990 MC=1,2*LC+1
                      DO 980 IY=1,3
                          DO 970 IX=1,3
                              DO 960 MB=1,2*LB+1
                                  DL(IBB+MB,IX,IY,IBC+MC) = 
     .                                DL(IBB+MB,IX,IY,IBC+MC) 
     .                                    + CB*AT(MC,IY,MB,IX)
  960                         CONTINUE
  970                     CONTINUE
  980                 CONTINUE
  990             CONTINUE
 1000         CONTINUE
 1900     CONTINUE
 2000 CONTINUE
C
C    Put all integrals into the MNDO-94 order.
C
      DO 2900 IY=1,3
          DO 2890 IX=1,3
              CALL PSSX94(NORBB,NORBC,DL(1,IX,IY,1),9*((MAXL+1)**2))
 2890     CONTINUE
 2900 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSA3PG(IGO,IB,IC,RB,RC,DL)
C
C   Compute two-center dipole-orbital moment contributions of the <B|A|C> 
C   type (paramagnetic part of the HAB matrix) using Gaussian expansions.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IGO    - Order of gaussian expansion.
C      IB     - Index of the left-hand atom
C      IC     - Index of the right-hand atom
C      RB     - Coordinates of the left-hand atom relative to the operator
C               center (atomic units).
C      RC     - Coordinates of the right-hand atom relative to the operator
C               center (atomic units).
C      DL     - Output <r_\mu \mu|r^-3 L^\nu|\nu> integrals. 
C               First index: Orbitals at the left-hand center (B) in the
C                   MNDO94 order.
C               Second index: Components of the radius-vector at B.
C               Third index: Projections of the orbital moment operator
C                   (a.k.a. magnetic moment components).
C               Forth index: Orbitals at the right-hand center (C) in the
C                   MNDO94 order.
C
C   Accessed common blocks:
C
C      ATOMS  - Atoms identities and range of orbitals
C      PSNPAR - NMR parameters
C
C   Local storage:
C
C      9*(2*MAXL+1)**2 DOUBLE PRECISION words. With MAXL=2, it's some 225
C      words, or 1.8K.
C
C   Module logic:
C
C      This is gaussian version of PSAB3P, see that for comments.
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION DL((MAXL+1)**2,3,3,(MAXL+1)**2)
      DIMENSION RB(3), RC(3)
C
      DIMENSION DT(2*MAXL+1,3,3,2*MAXL+1)
C
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPB = NAT(IB)
      NORBC = NLAST(IC)-NFIRST(IC)+1
      ITYPC = NAT(IC)
C    Clear output
      DO 200 I2=1,NORBC
        DO 190 IY=1,3
          DO 180 IX=1,3
            DO 170 I1=1,NORBB
                DL(I1,IX,IY,I2) = ZERO
  170       CONTINUE
  180     CONTINUE
  190   CONTINUE
  200 CONTINUE
C
      DO 2000 LB=0,MAXLA(ITYPB)
          IBASB = LB**2
          DO 1900 JB=1,MAXN(LB,ITYPB)
              NB = NN(JB,LB,ITYPB)
              ZB = ZN(JB,LB,ITYPB)
              CB = CN(JB,LB,ITYPB)
              DO 1000 LC=0,MAXLA(ITYPC)
                  IBASC = LC**2
                  DO 600 JC=1,MAXN(LC,ITYPC)
                      NC = NN(JC,LC,ITYPC)
                      ZC = ZN(JC,LC,ITYPC)
                      CC = CN(JC,LC,ITYPC)
                      CALL PSGRL3(IGO,NB,LB,ZB,RB,NC,LC,ZC,RC,
     .                             DT,2*MAXL+1)
                      SC = CB*CC
                      DO 400 MC=1,2*LC+1
                        DO 390 IY=1,3
                          DO 380 IX=1,3
                            DO 370 MB=1,2*LB+1
                                I1 = IBASB+MB
                                I2 = IBASC+MC
                                DL(I1,IX,IY,I2) = DL(I1,IX,IY,I2)
     .                                          + SC*DT(MB,IX,IY,MC)
  370                       CONTINUE
  380                     CONTINUE
  390                   CONTINUE
  400                 CONTINUE
  600             CONTINUE
 1000         CONTINUE
 1900     CONTINUE
 2000 CONTINUE
C
C    Put all integrals into the MNDO-94 order.
C
      DO 2900 IY=1,3
          DO 2890 IX=1,3
              CALL PSSX94(NORBB,NORBC,DL(1,IX,IY,1),9*((MAXL+1)**2))
 2890     CONTINUE
 2900 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSAB2D(IA,IB,RX,RY,RZ,DD)
C
C   Compute two-center quadrupole contributions of the <B|A|A> type 
C   to the HAB matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Index of the right-hand atom (one carrying operator)
C      IB     - Index of the left-hand atom
C      RX,RY,RZ
C             - Relative coordinates of the atom A (right-hand one).
C               Coordinates are in Bohrs.
C      DD     - Output <\mu|r^-3 (\delta_{ab} r^2 - r_b r_a|\nu> integrals
C
C   Accessed common blocks:
C
C      ATOMS  - Atoms identities and range of orbitals
C      PSNPAR - NMR parameters
C
C   Local storage:
C
C      Next to nothing.
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION DD((MAXL+1)**2,3,3,(MAXL+1)**2)
      DIMENSION DDT(2*MAXL+1,3,3,2*MAXL+1)
C
      NORBA = NLAST(IA)-NFIRST(IA)+1
      ITYPA = NAT(IA)
C
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPB = NAT(IB)
C
C   MU is on the atom B, NU is on the atom A
C
      DO 800 LA=0,MAXLA(ITYPA)
          IBASA = LA**2
          DO 700 LB=0,MAXLA(ITYPB)
              IBASB = LB**2
              DO 200 MA=1,2*LA+1
                  DO 190 IY=1,3
                      DO 180 IX=1,3
                          DO 170 MB=1,2*LB+1
                              DD(IBASB+MB,IX,IY,IBASA+MA) = ZERO
  170                     CONTINUE
  180                 CONTINUE
  190             CONTINUE
  200         CONTINUE
              DO 600 ILA=1,MAXN(LA,ITYPA)
                  NA  = NN(ILA,LA,ITYPA)
                  ZA  = ZN(ILA,LA,ITYPA)
                  SCA = CN(ILA,LA,ITYPA)
                  DO 500 ILB=1,MAXN(LB,ITYPB)
                      NB  = NN(ILB,LB,ITYPB)
                      ZB  = ZN(ILB,LB,ITYPB)
                      SCB = CN(ILB,LB,ITYPB)
                      CALL PSOQR3(NB,LB,ZB,NA,LA,ZA,RX,RY,RZ,
     .                            DDT,2*MAXL+1)
                      DO 400 MA=1,2*LA+1
                          DO 390 IY=1,3
                              DO 380 IX=1,3
                                  DO 370 MB=1,2*LB+1
                                      DD(IBASB+MB,IX,IY,IBASA+MA) = DD(
     .                                       IBASB+MB,IX,IY,IBASA+MA) + 
     .                                       SCA*SCB*DDT(MB,IX,IY,MA)
  370                             CONTINUE
  380                         CONTINUE
  390                     CONTINUE
  400                 CONTINUE
  500             CONTINUE
  600         CONTINUE
  700     CONTINUE
  800 CONTINUE
C
C        Put all integrals into the MNDO-94 order.
C
      DO 1900 IY=1,3
          DO 1800 IX=1,3
              CALL PSSX94(NORBB,NORBA,DD(1,IX,IY,1),9*((MAXL+1)**2))
 1800     CONTINUE
 1900 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSABDR(NORBA,NORBB,DD)
C                        ^
C   Convert <\mu|X_a Y_b O |\nu> integrals into
C   <\mu|/\_{ab} X . Y - X_a Y_b |\nu> integrals.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NORBA  - Number of orbitals on the left-hand center.
C      NORBB  - Number of orbitals on the right-hand center.
C      DD     - Input: <\mu|r_a r_b O|\nu> integrals.
C                   First index: orbitals on the left-hand center.
C                   Second index: components of the first vector
C                       operator.
C                   Third index: components of the second vector
C                       operator.
C                   Forth index: orbitals at the right-hand center.
C               Output integrals, four-index array. 
C                   Second and third index are the components of the
C                   combined quadrupole operator. First and forth
C                   indices are unchanged.
C
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C      Nope.
C
C   Module logic:
C
C      Nope.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (ZERO=0.D0)
C
      DIMENSION DD((MAXL+1)**2,3,3,(MAXL+1)**2)
C
      DO 1500 MB=1,NORBB
          DO 1400 MA=1,NORBA
              TOT = ZERO
              DO 1200 IY=1,3
                  TOT = TOT + DD(MA,IY,IY,MB)
                  DO 1190 IX=1,3
                      DD(MA,IX,IY,MB) = -DD(MA,IX,IY,MB)
 1190             CONTINUE
 1200         CONTINUE
              DO 1300 IX=1,3
                  DD(MA,IX,IX,MB) = DD(MA,IX,IX,MB) + TOT
 1300         CONTINUE
 1400     CONTINUE
 1500 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSAB2E(IA,IB,RX,RY,RZ,DD)
C
C   Compute two-center quadrupole contributions of the <A|A|B> type 
C   to the HAB matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Index of the left-hand atom.
C      IB     - Index of the right-hand atom.
C      RX,RY,RZ
C             - Relative coordinates of the atom A (left-hand one) in
C               atomic units.
C      DD     - Output integrals.
C               First index: orbitals on the left-hand atom (MNDO94 order)
C               Second index: components of the left-center dipole, X,Y,Z
C               Third index: components of the right-center dipole, X,Y,Z
C               Forth index: orbitals on the right-hand atom (MNDO94 order)
C
C   Accessed common blocks:
C
C      ATOMS  - Atoms identities and range of orbitals
C      PSNPAR - NMR parameters
C
C   Local storage:
C
C      Forget it!
C
C   Module logic:
C
C      We are talking about these integrals:
C                       A   B    A  B
C                 /\   r . r  - r  r
C           A       ab           a  b        B
C      < \mu   | ---------------------- | \nu  >
C                           3
C                          A
C                         r 
C
C      We reduce 'em to < \mu^A r^A | r^A^-3 | \nu^B r^B >....
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION DD ((MAXL+1)**2,3,3,(MAXL+1)**2)
      DIMENSION DDT(2*MAXL+1,3,3,2*MAXL+1)
C
      NORBA = NLAST(IA)-NFIRST(IA)+1
      ITYPA = NAT(IA)
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPB = NAT(IB)
C
      DO 500 MB=1,NORBB
          DO 490 IY=1,3
              DO 480 IX=1,3
                  DO 470 MA=1,NORBA
                      DD(MA,IX,IY,MB) = ZERO
  470             CONTINUE
  480         CONTINUE
  490     CONTINUE
  500 CONTINUE
C
      DO 2000 LB=0,MAXLA(ITYPB)
          IBASB = LB**2
          DO 1900 LA=0,MAXLA(ITYPA)
              IBASA = LA**2
              DO 1800 JB=1,MAXN(LB,ITYPB)
                  NB = NN(JB,LB,ITYPB)
                  ZB = ZN(JB,LB,ITYPB)
                  CB = CN(JB,LB,ITYPB)
                  DO 1700 JA=1,MAXN(LA,ITYPA)
                      NA = NN(JA,LA,ITYPA)
                      ZA = ZN(JA,LA,ITYPA)
                      CA = CN(JA,LA,ITYPA)
                      CALL PSRAB3(NA,LA,ZA,NB,LB,ZB,RX,RY,RZ,
     .                            DDT,2*MAXL+1)
                      DO 1000 MB=1,2*LB+1
                          DO 990 IY=1,3
                              DO 980 IX=1,3
                                  DO 970 MA=1,2*LA+1
                                      DD(IBASA+MA,IX,IY,IBASB+MB) = 
     .                                    DD(IBASA+MA,IX,IY,IBASB+MB) + 
     .                                        CA*CB*DDT(MA,IX,IY,MB)
  970                             CONTINUE
  980                         CONTINUE
  990                     CONTINUE
 1000                 CONTINUE
 1700             CONTINUE
 1800         CONTINUE
 1900     CONTINUE
 2000 CONTINUE
C    Reduce all integrals to the quadrupole form
      CALL PSABDR(NORBA,NORBB,DD)
C    Put all integrals into the MNDO-94 order.
      DO 3900 IY=1,3
          DO 3800 IX=1,3
              CALL PSSX94(NORBA,NORBB,DD(1,IX,IY,1),9*((MAXL+1)**2))
 3800     CONTINUE
 3900 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSAB2B(IB,RX,RY,RZ,DD)
C
C   Compute two-center quadrupole contributions of the <B|A|B> type 
C   to the HAB(D) matrix, i.e.
C            \delta_{ab} r^A . r^B - r^B_b r^A_a
C   <\mu^B | ----------------------------------- | \nu^B >
C                               3
C                            r^A
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IB     - Index of the atom with orbitals.
C      RX,RY,RZ
C             - Relative coordinates of the atom A (i.e., of the operator)
C               Coordinates are in Bohrs.
C      DD     - Output integrals. 
C               First index: left-hand orbitals in the MNDO94 order.
C               Second index: components of the radius-vector centered
C                   at the operator (a.k.a. projections of the magnetic field),
C                   in the X,Y,Z order.
C               Third index: components of the radius-vector centered at the
C                   orbitals center (a.k.a. orientations of the magnetic moment),
C                   in the X,Y,Z order.
C               Forth index: right-hand orbitals in the MNDO94 order.
C
C   Accessed common blocks:
C
C      ATOMS  - Atoms identities and range of orbitals
C      PSNPAR - NMR parameters
C
C   Local storage:
C
C      Some 9*(2*MAXL+1)*(MAXL+1)**2 (ca. 400) DOUBLE PRECISION
C      cells.
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION DD((MAXL+1)**2,3,3,(MAXL+1)**2)
C
      DIMENSION TD(2*MAXL+1,3,2*(MAXL-1)+1)
      DIMENSION TU(2*MAXL+1,3,2*(MAXL+1)+1)
      DIMENSION TX(2*MAXL+1,2*MAXL+1,3)
C
      NORB = NLAST(IB)-NFIRST(IB)+1
      ITYP = NAT(IB)
C
C                     r^B_b r^A_a
C    Compute <\mu | -------------- | \nu> auxiliary integrals
C                            3
C                         r^A
C
      DO 900 L1=0,MAXLA(ITYP)
          IBAS1 = L1**2
          DO 800 L2=0,MAXLA(ITYP)
              IBAS2 = L2**2
              DO 200 M2=1,2*L2+1
                  DO 190 IY=1,3
                      DO 180 IX=1,3
                          DO 170 M1=1,2*L1+1
                              DD(IBAS1+M1,IX,IY,IBAS2+M2) = ZERO
  170                     CONTINUE
  180                 CONTINUE
  190             CONTINUE
  200         CONTINUE
              DO 700 I2=1,MAXN(L2,ITYP)
                  N2 = NN(I2,L2,ITYP)
                  Z2 = ZN(I2,L2,ITYP)
                  C2 = CN(I2,L2,ITYP)
                  DO 600 I1=1,MAXN(L1,ITYP)
                      N1 = NN(I1,L1,ITYP)
                      Z1 = ZN(I1,L1,ITYP)
                      C1 = CN(I1,L1,ITYP)
                      IF( L2.GE.1 )
     .                CALL PSOPR3(N1,L1,Z1,N2+1,L2-1,Z2,RX,RY,RZ,
     .                            TD,2*MAXL+1)
                      CALL PSOPR3(N1,L1,Z1,N2+1,L2+1,Z2,RX,RY,RZ,
     .                            TU,2*MAXL+1)
                      DO 500 IX=1,3
                          CALL PSODS(2*L1+1,N2,L2,Z2,TD(1,IX,1),
     .                             3*(2*MAXL+1),TU(1,IX,1),3*(2*MAXL+1),
     .                             TX,2*MAXL+1,2*MAXL+1)
                          DO 400 IY=1,3
                              DO 390 M2=1,2*L2+1
                                  IM2 = IBAS2 + M2
                                  DO 380 M1=1,2*L1+1
                                      IM1 = IBAS1 + M1
                                      DD(IM1,IX,IY,IM2) = DD(IM1,IX,IY,
     .                                    IM2) + C1*C2*TX(M1,M2,IY)
  380                             CONTINUE
  390                         CONTINUE
  400                     CONTINUE
  500                 CONTINUE
  600             CONTINUE
  700         CONTINUE
  800     CONTINUE
  900 CONTINUE
C    Convert integrals to the quadrupole form
      CALL PSABDR(NORB,NORB,DD)
C    Put all integrals into the MNDO-94 order.
      DO 1900 IY=1,3
          DO 1800 IX=1,3
              CALL PSSX94(NORB,NORB,DD(1,IX,IY,1),9*((MAXL+1)**2))
 1800     CONTINUE
 1900 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSAB3D(IB,IC,RB,RC,DD)
C
C   Compute two-center quadrupole contributions of the <B|A|C> type 
C   to the HAB(D) matrix, i.e.
C                   A   C    C  A
C             /\   r . r  - r  r
C               ab           b  a
C    <\mu  | ---------------------- | \nu  >
C        B               3               C
C                       A
C                      r 
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IB     - Index of the left-hand atom.
C      IC     - Index of the right-hand atom.
C      RB     - Position of the left-hand atom relative to the operator
C               center, atomic units.
C      RC     - Position of the right-hand atom relative to the operator
C               center, atomic units.
C      DD     - Output integrals, four-index array. 
C               First index: orbitals on the left-hand center.
C               Second index: operator projection (a.k.a. magnetic
C                   field component).
C               Third index: projection of the radius-vector at C (a.k.a.
C                   magnetic dipole component).
C               Forth index: orbitals at the right-hand center.
C
C   Accessed common blocks:
C
C      ATOMS  - Atoms identities and range of orbitals
C      PSNPAR - NMR parameters
C
C   Local storage:
C
C      Some 9*(2*MAXL+1)*(MAXL+1)**2 (ca. 400) DOUBLE PRECISION storage
C      units.
C
C   Module logic:
C
C      This is a replica of PSAB2B.
C
C   Possible optimizations:
C
C      Calls PSODS can be moved out of the inner loop by introducing
C      intermediate temporary.
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION RB(3), RC(3)
      DIMENSION DD((MAXL+1)**2,3,3,(MAXL+1)**2)
C
      DIMENSION TD(2*MAXL+1,3,2*(MAXL-1)+1)
      DIMENSION TU(2*MAXL+1,3,2*(MAXL+1)+1)
      DIMENSION TX(2*MAXL+1,2*MAXL+1,3)
C
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPB = NAT(IB)
      NORBC = NLAST(IC)-NFIRST(IC)+1
      ITYPC = NAT(IC)
C                      C   A
C                     r   r
C                B     b   a         C
C    Compute <\mu | ----------- | \nu > auxiliary integrals
C                        3  
C                       A 
C                      r
C
      DO 900 LB=0,MAXLA(ITYPB)
          IBASB = LB**2
          DO 800 LC=0,MAXLA(ITYPC)
              IBASC = LC**2
              DO 200 MC=1,2*LC+1
                  DO 190 IY=1,3
                      DO 180 IX=1,3
                          DO 170 MB=1,2*LB+1
                              DD(IBASB+MB,IX,IY,IBASC+MC) = ZERO
  170                     CONTINUE
  180                 CONTINUE
  190             CONTINUE
  200         CONTINUE
              DO 700 JC=1,MAXN(LC,ITYPC)
                  NC = NN(JC,LC,ITYPC)
                  ZC = ZN(JC,LC,ITYPC)
                  CC = CN(JC,LC,ITYPC)
                  DO 600 JB=1,MAXN(LB,ITYPB)
                      NB = NN(JB,LB,ITYPB)
                      ZB = ZN(JB,LB,ITYPB)
                      CB = CN(JB,LB,ITYPB)
                      IF( LC.GT.0 )
     .                CALL PS3PR3(NB,LB,ZB,RB,NC+1,LC-1,ZC,RC,
     .                            TD,2*MAXL+1)
                      CALL PS3PR3(NB,LB,ZB,RB,NC+1,LC+1,ZC,RC,
     .                            TU,2*MAXL+1)
                      DO 500 IX=1,3
                          CALL PSODS(2*LB+1,NC,LC,ZC,
     .                               TD(1,IX,1),3*(2*MAXL+1),
     .                               TU(1,IX,1),3*(2*MAXL+1),
     .                               TX,2*MAXL+1,2*MAXL+1)
                          DO 400 IY=1,3
                              DO 390 MC=1,2*LC+1
                                  IMC = IBASC + MC
                                  DO 380 MB=1,2*LB+1
                                      IMB = IBASB + MB
                                      DD(IMB,IX,IY,IMC) = DD(IMB,IX,IY,
     .                                    IMC) + CB*CC*TX(MB,MC,IY)
  380                             CONTINUE
  390                         CONTINUE
  400                     CONTINUE
  500                 CONTINUE
  600             CONTINUE
  700         CONTINUE
  800     CONTINUE
  900 CONTINUE
C    Convert integrals to the quadrupole form:
      CALL PSABDR(NORBB,NORBC,DD)
C    Put all integrals into the MNDO-94 order.
      DO 1900 IY=1,3
          DO 1800 IX=1,3
              CALL PSSX94(NORBB,NORBC,DD(1,IX,IY,1),9*((MAXL+1)**2))
 1800     CONTINUE
 1900 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSA3DG(IGO,IB,IC,RB,RC,DD)
C
C   Compute two-center quadrupole contributions of the <B|A|C> type 
C   to the HAB(D) matrix using Gaussian expansion.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IGO    - Order of Gaussian expansion to use.
C      IB     - Index of the left-hand atom.
C      IC     - Index of the right-hand atom.
C      RB     - Position of the left-hand atom relative to the operator
C               center, atomic units.
C      RC     - Position of the right-hand atom relative to the operator
C               center, atomic units.
C      DD     - Output integrals, four-index array. 
C               First index: orbitals on the left-hand center.
C               Second index: operator projection (a.k.a. magnetic
C                   field component).
C               Third index: projection of the radius-vector at C (a.k.a.
C                   magnetic dipole component).
C               Forth index: orbitals at the right-hand center.
C
C   Accessed common blocks:
C
C      ATOMS  - Atoms identities and range of orbitals
C      PSNPAR - NMR parameters
C
C   Local storage:
C
C      9*(2*MAXL+1)**2 DOUBLE PRECISION words. With MAXL=2, its some 225 DP
C      words, or 1.8K.
C
C   Module logic:
C
C      This is an alternative version of PSAB3D. See it for comments.
C
C   Possible optimizations:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION RB(3), RC(3)
      DIMENSION DD((MAXL+1)**2,3,3,(MAXL+1)**2)
C
      DIMENSION DT(2*MAXL+1,3,3,2*MAXL+1)
C
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPB = NAT(IB)
      NORBC = NLAST(IC)-NFIRST(IC)+1
      ITYPC = NAT(IC)
C    Zero output array
      DO 100 I2=1,NORBC
        DO 90 IY=1,3
          DO 80 IX=1,3
            DO 70 I1=1,NORBB
                DD(I1,IX,IY,I2) = ZERO
   70       CONTINUE
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
C
      DO 900 LB=0,MAXLA(ITYPB)
          IBASB = LB**2
          DO 800 LC=0,MAXLA(ITYPC)
              IBASC = LC**2
              DO 700 JC=1,MAXN(LC,ITYPC)
                  NC = NN(JC,LC,ITYPC)
                  ZC = ZN(JC,LC,ITYPC)
                  CC = CN(JC,LC,ITYPC)
                  DO 600 JB=1,MAXN(LB,ITYPB)
                      NB = NN(JB,LB,ITYPB)
                      ZB = ZN(JB,LB,ITYPB)
                      CB = CN(JB,LB,ITYPB)
                      CALL PSGRR3(IGO,NB,LB,ZB,RB,NC,LC,ZC,RC,
     .                             DT,2*MAXL+1)
                      SC = CB*CC
                      DO 400 MC=1,2*LC+1
                        DO 390 IY=1,3
                          DO 380 IX=1,3
                            DO 370 MB=1,2*LB+1
                                I1 = IBASB+MB
                                I2 = IBASC+MC
                                DD(I1,IX,IY,I2) = DD(I1,IX,IY,I2) 
     .                                          + SC*DT(MB,IX,IY,MC)
  370                       CONTINUE
  380                     CONTINUE
  390                   CONTINUE
  400                 CONTINUE
  600             CONTINUE
  700         CONTINUE
  800     CONTINUE
  900 CONTINUE
C    Convert integrals to the quadrupole form:
      CALL PSABDR(NORBB,NORBC,DD)
C    Put all integrals into the MNDO-94 order.
      DO 1900 IY=1,3
          DO 1800 IX=1,3
              CALL PSSX94(NORBB,NORBC,DD(1,IX,IY,1),9*((MAXL+1)**2))
 1800     CONTINUE
 1900 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSNH0B(IA,H0B,LDH1,LDH2,RO,RO1,LDR)
C
C   Compute upper triangle of each of the three magnetic moment
C   hamiltonian derivatives for a given atom. Derivatives are
C   in units of i/c, where i is Sqrt[-1] and c is speed of light.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Atom posessing the magnetic moment.
C               IA>NATOM means NICS.
C      H0B    - Output H0B matrices. Values are provided
C               in the upper triangles of H0B. Lower
C               triangle contains data for use by PSNHAB
C               (see below).
C      LDH1,
C      LDH2   - Leading dimensions of the H0B matrix.
C      RO     - Zero-order density matrix.
C      RO1    - First-order density matrices.
C      LDR    - Leading dimension of RO1.
C
C   Accessed common blocks:
C
C      ATOMS  - Range of orbitals carried by an atom
C      ATOMC  - Coordinates of atoms (don't you wish it was in a.u.s, boy?)
C      PSDOPT - Computation options, in particular NMRLEV.
C      PSDGBL - Number of atoms, for consistency with derivatives code.
C      PSPRTF - Debugging output control.
C
C   Modified common blocks:
C
C      PSO3PR - Precision control of three-center integrals (EPSL is used).
C
C   Local storage:
C
C      Some 3*(MAXL+1)**4 (ca. 240) DOUBLE PRECISION cells.
C
C   Module logic:
C                                           ^
C                                           L
C      Matrix elements of the operator <- ------> (note the minus sign!)
C                                              3
C                                        (r-R )
C                                            A
C      are evaluated here.
C
C      Although in the ideal world, H0B should be automatically antisymmetric
C      (since it is an imaginary part of a Hermitial operator matrix),
C      we compute the whole matrix as a precaution against lost of precision
C      in integral routines and antisymmetrise the result.
C
C      Output matrices actually contain *whole* matrix rather than just
C      an upper triangle, since this is the behaiviour required by PSNHAB.
C      However, this may change in the future, so that no other routines
C      should rely on lower triangle being correctly initialized.
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, MXNICS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI, VARINT
      LOGICAL PSNINC
      CHARACTER*1 VARN
      PARAMETER (BOHR=0.529167D0)
      PARAMETER (MAXL=2)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ATOMC / COORD(3,LM1)
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
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
     ./PSNICS/ COORDN(3,MXNICS), NNICS
     ./PSO3PR/ EPSL, EPSR
      SAVE /PSDGBL/, /PSPRTF/, /PSDOPT/, /PSNICS/, /PSO3PR/, /PSPRT /
      DIMENSION VARN(3), RB(3), RC(3)
      DIMENSION H0B(LDH1,LDH2,3), RO(LDR), RO1(LDR,3)
      DIMENSION DX((MAXL+1)**2,(MAXL+1)**2,3)
      DATA VARN/'X','Y','Z'/
C
      VARINT = MOD(INTCTL,200).GE.100
      IGO    = MAX(MOD(INTCTL,20)-10,0)
C
C    Clear arrays first, so that we can accumulate contributions
C    without thinking about previous state.
C
      CALL PSNZR(IA,H0B(1,1,1),LDH1)
      CALL PSNZR(IA,H0B(1,1,2),LDH1)
      CALL PSNZR(IA,H0B(1,1,3),LDH1)
C
C    One-center contributions, centered at the IA atom.
C    No one-center contributions to NICS.
C
      IF( IA.LE.NATOM ) THEN
          IBASA = NFIRST(IA)-1
          NORBA = NLAST(IA)-NFIRST(IA)+1
          XA = COORD(1,IA)
          YA = COORD(2,IA)
          ZA = COORD(3,IA)
          CALL PSNHB1(IA,DX)
          DO 600 M2=1,NORBA
              IH2 = IBASA+M2
              DO 550 M1=1,NORBA
                  IH1 = IBASA+M1
                  H0B(IH1,IH2,1) = H0B(IH1,IH2,1) - DX(M1,M2,1)
                  H0B(IH1,IH2,2) = H0B(IH1,IH2,2) - DX(M1,M2,2)
                  H0B(IH1,IH2,3) = H0B(IH1,IH2,3) - DX(M1,M2,3)
  550         CONTINUE
  600     CONTINUE
      ELSE
          XA = COORDN(1,IA-NATOM)
          YA = COORDN(2,IA-NATOM)
          ZA = COORDN(3,IA-NATOM)
      ENDIF
C
C    Two-center contributions
C
      IF( NMRLEV.GT.1 ) THEN
          DO 1600 IB=1,NATOM
              IF( IB.EQ.IA ) GOTO 1599
              IBASB  = NFIRST(IB)-1
              NORBB  = NLAST(IB)-NFIRST(IB)+1
              RX     = (COORD(1,IB)-XA)/BOHR
              RY     = (COORD(2,IB)-YA)/BOHR
              RZ     = (COORD(3,IB)-ZA)/BOHR
C
C        Two-center contributions of the (B|A|B) type.
C
              CALL PSNB2A(IB,RX,RY,RZ,DX)
              DO 1100 M2=1,NORBB
                  IH2 = IBASB + M2
                  DO 1090 M1=1,NORBB
                      IH1 = IBASB + M1
                      H0B(IH1,IH2,1) = H0B(IH1,IH2,1) - DX(M1,M2,1)
                      H0B(IH1,IH2,2) = H0B(IH1,IH2,2) - DX(M1,M2,2)
                      H0B(IH1,IH2,3) = H0B(IH1,IH2,3) - DX(M1,M2,3)
 1090             CONTINUE
 1100         CONTINUE
C        Two-center contribution of the (B|A|A) type and (A|A|B) type.
C        The later are just a transposed conjugated version of the
C        first. We are going for the <B|A|A> integrals rather than for
C        the more obvious <A|A|B> once since <B|A|A> are a bit easier
C        to evaluate. Neither of those is present for NICS.
              IF( IA.LE.NATOM ) THEN
                  CALL PSNB2C(IB,IA,-RX,-RY,-RZ,DX)
                  DO 1200 M2=1,NORBA
                      IH2 = IBASA + M2
                      DO 1190 M1=1,NORBB
                          IH1 = IBASB + M1
C                        Note that IH1 is an orbital at the atom B, while
C                        IH2 is an orbital at the atom A, so that we are
C                        adding this new term to the correct location with
C                        the right sign.
                          H0B(IH1,IH2,1) = H0B(IH1,IH2,1) - DX(M1,M2,1)
                          H0B(IH1,IH2,2) = H0B(IH1,IH2,2) - DX(M1,M2,2)
                          H0B(IH1,IH2,3) = H0B(IH1,IH2,3) - DX(M1,M2,3)
 1190                 CONTINUE
 1200             CONTINUE
                  IF( INDSYM.EQ.2 ) THEN
                      DO 1300 M2=1,NORBA
                          IH2 = IBASA + M2
                          DO 1290 M1=1,NORBB
                              IH1 = IBASB + M1
                              H0B(IH2,IH1,1)=H0B(IH2,IH1,1)+DX(M1,M2,1)
                              H0B(IH2,IH1,2)=H0B(IH2,IH1,2)+DX(M1,M2,2)
                              H0B(IH2,IH1,3)=H0B(IH2,IH1,3)+DX(M1,M2,3)
 1290                     CONTINUE
 1300                 CONTINUE
                  ENDIF
                  IF( INDSYM.NE.2 ) THEN
C                    Two-center contributions of the <A|A|B> type should be
C                    evaluated in the painful way rather than by hermiticity.
                      CALL PSNB2D(IA,IB,RX,RY,RZ,DX)
                      DO 1400 M2=1,NORBB
                          IH2 = IBASB + M2
                          DO 1390 M1=1,NORBA
                              IH1 = IBASA + M1
                              H0B(IH1,IH2,1)=H0B(IH1,IH2,1)-DX(M1,M2,1)
                              H0B(IH1,IH2,2)=H0B(IH1,IH2,2)-DX(M1,M2,2)
                              H0B(IH1,IH2,3)=H0B(IH1,IH2,3)-DX(M1,M2,3)
 1390                     CONTINUE
 1400                 CONTINUE
                  ENDIF
              ENDIF
*            End of the part not present for NICS
 1599         CONTINUE
 1600     CONTINUE
      ENDIF
C
C    Three-center contributions. All should be retained for NICS.
C
      IF( NMRLEV.GT.2 ) THEN
          DO 3000 IB=1,NATOM
              IF( IB.EQ.IA ) GOTO 2999
              IBASB  = NFIRST(IB)-1
              NORBB  = NLAST(IB)-NFIRST(IB)+1
              RB(1)  = (COORD(1,IB)-XA)/BOHR
              RB(2)  = (COORD(2,IB)-YA)/BOHR
              RB(3)  = (COORD(3,IB)-ZA)/BOHR
              ICMAX  = NATOM
              IF( INDSYM.EQ.2 ) ICMAX = IB - 1
              DO 2800 IC=1,ICMAX
                  IF( IC.EQ.IA .OR. IC.EQ.IB ) GOTO 2799
                  IF( .NOT. PSNINC(IB,IC) ) GOTO 2799
                  IBASC  = NFIRST(IC)-1
                  NORBC  = NLAST(IC)-NFIRST(IC)+1
                  RC(1)  = (COORD(1,IC)-XA)/BOHR
                  RC(2)  = (COORD(2,IC)-YA)/BOHR
                  RC(3)  = (COORD(3,IC)-ZA)/BOHR
C
C                Variable integrals precision control
C
                  IF(VARINT) THEN
                      ROMAX = MAX(
     .                    PSNDMX(RO,      IBASB+1,NORBB,IBASC+1,NORBC),
     .                    PSNDMX(RO1(1,1),IBASB+1,NORBB,IBASC+1,NORBC),
     .                    PSNDMX(RO1(1,2),IBASB+1,NORBB,IBASC+1,NORBC),
     .                    PSNDMX(RO1(1,3),IBASB+1,NORBB,IBASC+1,NORBC))
                      ROMAX = 5D0*MAX(ROMAX,DNCOFF)
                      EPSR  = MAX(1D-6,DNCOFF/ROMAX)
                      EPSL  = MAX(2D-6,DNCOFF/ROMAX)
                  ENDIF
C
                  IF( IGO.EQ.0 ) THEN
                      CALL PSNB3(IB,IC,RB,RC,DX)
                  ELSE
                      CALL PSNB3G(IGO,IB,IC,RB,RC,DX)
                  ENDIF
                  DO 2000 M2=1,NORBC
                      IH2 = IBASC + M2
                      DO 1900 M1=1,NORBB
                          IH1 = IBASB + M1
                          H0B(IH1,IH2,1) = H0B(IH1,IH2,1) - DX(M1,M2,1)
                          H0B(IH1,IH2,2) = H0B(IH1,IH2,2) - DX(M1,M2,2)
                          H0B(IH1,IH2,3) = H0B(IH1,IH2,3) - DX(M1,M2,3)
 1900                 CONTINUE
 2000             CONTINUE
                  IF( INDSYM.EQ.2 ) THEN
                      DO 2300 M2=1,NORBC
                          IH2 = IBASC + M2
                          DO 2100 M1=1,NORBB
                              IH1 = IBASB + M1
                              H0B(IH2,IH1,1)=H0B(IH2,IH1,1)+DX(M1,M2,1)
                              H0B(IH2,IH1,2)=H0B(IH2,IH1,2)+DX(M1,M2,2)
                              H0B(IH2,IH1,3)=H0B(IH2,IH1,3)+DX(M1,M2,3)
 2100                     CONTINUE
 2300                 CONTINUE
                  ENDIF
 2799             CONTINUE
 2800         CONTINUE
 2999         CONTINUE
 3000     CONTINUE
      ENDIF
C
      IF(LPNMR) THEN
          DO 3300 IX=1,3
              IF( IA.LE.NATOM ) THEN
                  WRITE(NB6,11198) IA, VARN(IX)
              ELSE
                  WRITE(NB6,11199) IA-NATOM, VARN(IX)
              ENDIF
              CALL PSDPGM(NORBS,NORBS,H0B(1,1,IX),LDH1)
 3300     CONTINUE
      ENDIF
C
C    Antisymmetrize result
C
      DO 3400 IX=1,3
          CALL PSNMAS(IA,'H0B',VARN(IX),' ',H0B(1,1,IX),LDH1)
 3400 CONTINUE
C
      IF(LPNMR) THEN
          DO 3500 IX=1,3
              IF( IA.LE.NATOM ) THEN
                  WRITE(NB6,11200) IA, VARN(IX)
              ELSE
                  WRITE(NB6,11201) IA-NATOM, VARN(IX)
              ENDIF
              CALL PSDPGM(NORBS,NORBS,H0B(1,1,IX),LDH1)
 3500     CONTINUE
      ENDIF
C
      RETURN
11198 FORMAT(' RAW H0B MATRIX FOR ATOM ',I4,' (COMPONENT ',A,') IS:')
11199 FORMAT(' RAW H0B MATRIX FOR NICS POSITION ',I4,' (COMPONENT ',A,
     .       ') IS:')
11200 FORMAT(' SYMMETRIZED H0B MATRIX FOR ATOM ',I4,' (COMPONENT ',A,
     .       ') IS:')
11201 FORMAT(' SYMMETRIZED H0B MATRIX FOR NICS POSITION ',I4,
     .       ' (COMPONENT ',A,') IS:')
      END
C
      SUBROUTINE PSNHB1(IA,DX)
C
C   Compute one-center one-electron integrals necessary for the H0B
C   matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Index of the first atom
C      DL     - Output <\mu|r^-3 L|\nu> integrals
C
C   Accessed common blocks:
C
C      ATOMS  - Atoms identities and orbital ranges.
C      PSNPAR - NMR parameters
C
C   Local storage:
C
C      (MAXL+1)**4 (ca. 80) DOUBLE PRECISION cells.
C
C   Module logic:
C
C   Bugs:
C
C      Not mine, fortunately :-)
C
C      Replacing call with PSR31C by PSR31B will cause us
C      to reproduce (incorrect) results by Wu, You and Dai.
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION DX((MAXL+1)**2,(MAXL+1)**2,3)
C
      DIMENSION X00(2*MAXL+1,2*MAXL+1)
      EXTERNAL PSR31C
C
      NORB = NLAST(IA)-NFIRST(IA)+1
      ITYP = NAT(IA)
C
C    Zero out output integrals: only diagional blocks of the same
C    L will contain non-zero elements, and the rest we don't need.
C
      DO 60 IX=1,3
          DO 50 M2=1,NORB
              DO 40 M1=1,NORB
                  DX(M1,M2,IX) = ZERO
   40         CONTINUE
   50     CONTINUE
   60 CONTINUE
C
C    Compute hyperfine coupling integrals. They are non-zero only for
C    the orbital pairs with the same (l,m) values. Since the integral
C    corresponding to l=0 will vanish due to the orbital moment
C    operator, we do not need to compute corresponding coupling
C    integral (which is nice, as it is singular ;-)
C
      DO 800 L=1,MAXLA(ITYP)
          DVAL = ZERO
          IBAS = L**2 + 1
          DO 300 I1=1,MAXN(L,ITYP)
              N1 = NN(I1,L,ITYP)
              Z1 = ZN(I1,L,ITYP)
              C1 = CN(I1,L,ITYP)
              DO 200 I2=1,MAXN(L,ITYP)
                  N2   = NN(I2,L,ITYP)
                  Z2   = ZN(I2,L,ITYP)
                  C2   = CN(I2,L,ITYP)
                  DVAL = DVAL + C1*C2*PSR31C(N1,Z1,N2,Z2)
  200         CONTINUE
  300     CONTINUE
          DO 400 M1=1,2*L+1
              DO 390 M2=1,2*L+1
                  X00(M2,M1) = ZERO
  390         CONTINUE
              X00(M1,M1) = DVAL
  400     CONTINUE
          CALL PSOLX(2*L+1,L,X00,2*MAXL+1,DX(IBAS,IBAS,1),
     .                           (MAXL+1)**2,(MAXL+1)**2)
  800 CONTINUE
C
      CALL PSSX94(NORB,NORB,DX(1,1,1),(MAXL+1)**2)
      CALL PSSX94(NORB,NORB,DX(1,1,2),(MAXL+1)**2)
      CALL PSSX94(NORB,NORB,DX(1,1,3),(MAXL+1)**2)
C
      RETURN
      END
C
      SUBROUTINE PSNB2A(IB,RX,RY,RZ,DX)
C
C   Compute two-center one-electron integrals of the (B|A|B)
C   type necessary for H0B matrices
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IB     - Index of the atom carrying orbitals
C      RX,RY,RZ
C             - Relative position of the atom B, atomic units.
C      DX     - Output <\mu|r^-3 L|\nu> integrals.
C               First index: Left-hand orbitals, in the MNDO94 order.
C               Second index: Right-hand orbitals, in the MNDO94 order.
C               Third index: Components of the orbital moment operator,
C                   in the X,Y,Z order.
C
C   Accessed common blocks:
C
C      PSNPAR - NMR-specific parameters
C      ATOMS  - AO ranges
C
C   Local storage:
C
C      Some 3*(MAXL+1)**4 (ca. 240) DOUBLE PRECISION cells.
C
C   Module logic:
C
C   Bugs:
C
C      Algorithm used by PSLBR3 is accurate only to approximately
C      10^-10 atomic units.
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION DX((MAXL+1)**2,(MAXL+1)**2,3)
      DIMENSION DT(2*MAXL+1,3,2*MAXL+1)
C
C    Gather parameters
C
      NORB = NLAST(IB)-NFIRST(IB)+1
      ITYP = NAT(IB)
C
C    Compute integrals
C
      DO 800 L2=0,MAXLA(ITYP)
          IBAS2 = L2**2
          DO 700 L1=0,MAXLA(ITYP)
              IBAS1 = L1**2
              DO 100 M2=1,2*L2+1
                  DO 90 IX=1,3
                      DO 80 M1=1,2*L1+1
                          DX(IBAS1+M1,IBAS2+M2,IX) = ZERO
   80                 CONTINUE
   90             CONTINUE
  100         CONTINUE
              DO 600 I2=1,MAXN(L2,ITYP)
                  N2 = NN(I2,L2,ITYP)
                  Z2 = ZN(I2,L2,ITYP)
                  C2 = CN(I2,L2,ITYP)
                  DO 500 I1=1,MAXN(L1,ITYP)
                      N1 = NN(I1,L1,ITYP)
                      Z1 = ZN(I1,L1,ITYP)
                      C1 = CN(I1,L1,ITYP)
                      CALL PSLBR3(N1,L1,Z1,N2,L2,Z2,-RX,-RY,-RZ,
     .                                     DT,2*MAXL+1)
                      DO 400 M2=1,2*L2+1
                          DO 300 IX=1,3
                              DO 200 M1=1,2*L1+1
                                  DX(IBAS1+M1,IBAS2+M2,IX) = 
     .                                DX(IBAS1+M1,IBAS2+M2,IX) 
     .                                + C1*C2*DT(M1,IX,M2)
  200                         CONTINUE
  300                     CONTINUE
  400                 CONTINUE
  500             CONTINUE
  600         CONTINUE
  700     CONTINUE
  800 CONTINUE
C
      CALL PSSX94(NORB,NORB,DX(1,1,1),(MAXL+1)**2)
      CALL PSSX94(NORB,NORB,DX(1,1,2),(MAXL+1)**2)
      CALL PSSX94(NORB,NORB,DX(1,1,3),(MAXL+1)**2)
C
      RETURN
      END
C
      SUBROUTINE PSNB2C(IA,IB,RX,RY,RZ,DX)
C
C   Compute two-center one-electron integrals of the (A|B|B)
C   type necessary for H0B matrices
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Index of the atom with the left set of the orbitals.
C      IB     - Index of the atom carrying right set of the orbitals
C               and the magnetic moment.
C      RX,RY,RZ
C             - Relative position of the atom right-hand atom (B), 
C               atomic units.
C      DL     - Output <\mu|r^-3 L|\nu> integrals. Orbitals of the
C               atom A run over first index.
C
C   Accessed common blocks:
C
C      ATOMS  - Atoms identities and orbital ranges
C      PSNPAR - NMR parameters.
C
C   Local storage:
C
C      Some (MAXL+1)**4 (ca. 80) DOUBLE PRECISION cells.
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION DX((MAXL+1)**2,(MAXL+1)**2,3)
      DIMENSION ORT(2*MAXL+1,2*MAXL+1)
      DIMENSION OR (2*MAXL+1,2*MAXL+1)
C
C    Gather parameters
C
      NORBA = NLAST(IA)-NFIRST(IA)+1
      ITYPA = NAT(IA)
C
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPB = NAT(IB)
C
C    Integrals with L=0 at the atom B are always zero:
C
      DO 50 IORBA=1,NORBA
          DX(IORBA,1,1) = ZERO
          DX(IORBA,1,2) = ZERO
          DX(IORBA,1,3) = ZERO
   50 CONTINUE
C
C    Compute auxiliary r^-3 integrals and act by the orbital
C    moment operator on 'em.
C
      DO 800 LB=1,MAXLA(ITYPB)
          IBASB = LB**2 + 1
          DO 700 LA=0,MAXLA(ITYPA)
              IBASA = LA**2 + 1
              DO 200 MB=1,2*LB+1
                  DO 100 MA=1,2*LA+1
                      OR(MA,MB) = ZERO
  100             CONTINUE
  200         CONTINUE
              DO 600 ILB=1,MAXN(LB,ITYPB)
                  NB  = NN(ILB,LB,ITYPB)
                  ZB  = ZN(ILB,LB,ITYPB)
                  CB  = CN(ILB,LB,ITYPB)
                  DO 500 ILA=1,MAXN(LA,ITYPA)
                      NA  = NN(ILA,LA,ITYPA)
                      ZA  = ZN(ILA,LA,ITYPA)
                      CA  = CN(ILA,LA,ITYPA)
                      CALL PSORX(NA,LA,ZA,NB,LB,ZB,RX,RY,RZ,3,
     .                                    ORT,2*MAXL+1)
                      DO 400 MB=1,2*LB+1
                          DO 300 MA=1,2*LA+1
                              OR(MA,MB) = OR(MA,MB) 
     .                                  + CA*CB*ORT(MA,MB)
  300                     CONTINUE
  400                 CONTINUE
  500             CONTINUE
  600         CONTINUE
              CALL PSOLX(2*LA+1,LB,OR,2*MAXL+1,DX(IBASA,IBASB,1),
     .                                (MAXL+1)**2,(MAXL+1)**2)
  700     CONTINUE
  800 CONTINUE
C
      CALL PSSX94(NORBA,NORBB,DX(1,1,1),(MAXL+1)**2)
      CALL PSSX94(NORBA,NORBB,DX(1,1,2),(MAXL+1)**2)
      CALL PSSX94(NORBA,NORBB,DX(1,1,3),(MAXL+1)**2)
C
      RETURN
      END
C
      SUBROUTINE PSNB2D(IA,IB,RX,RY,RZ,DX)
C
C   Compute two-center one-electron integrals of the <A|A|B>
C   type necessary for H0B matrices.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Index of the atom with the left set of the orbitals
C               and operator.
C      IB     - Index of the atom carrying right set of the orbitals.
C      RX,RY,RZ
C             - Relative position of the right-hand atom (B),
C               atomic units.
C      DX     - Output <\mu^A|r^A^-3 L^A|\nu^B> integrals.
C               First index: Orbitals at the left-hand atom (A) in
C                   the MNDO94 order.
C               Second index: Orbitals at the right-hand atom (B) in
C                   the MNDO94 order.
C               Third index: Components of the orbital moment operator
C                   in the X,Y,Z order.
C 
C   Accessed common blocks:
C
C      ATOMS  - Atoms identities and orbital ranges
C      PSNPAR - NMR parameters.
C
C   Local storage:
C
C      A little.
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION DX((MAXL+1)**2,(MAXL+1)**2,3)
      DIMENSION DT(2*MAXL+1,3,2*MAXL+1)
C
      NORBA = NLAST(IA)-NFIRST(IA)+1
      ITYPA = NAT(IA)
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPB = NAT(IB)
      DO 500 IX=1,3
          DO 490 MB=1,NORBB
              DO 480 MA=1,NORBA
                  DX(MA,MB,IX) = ZERO
  480         CONTINUE
  490     CONTINUE
  500 CONTINUE
C
      DO 3000 LB=0,MAXLA(ITYPB)
          IBASB = LB**2
          DO 2900 LA=0,MAXLA(ITYPA)
              IBASA = LA**2
              DO 2800 JB=1,MAXN(LB,ITYPB)
                  NB = NN(JB,LB,ITYPB)
                  ZB = ZN(JB,LB,ITYPB)
                  CB = CN(JB,LB,ITYPB)
                  DO 2700 JA=1,MAXN(LA,ITYPA)
                      NA = NN(JA,LA,ITYPA)
                      ZA = ZN(JA,LA,ITYPA)
                      CA = CN(JA,LA,ITYPA)
                      CALL PSLA3(NA,LA,ZA,NB,LB,ZB,-RX,-RY,-RZ,
     .                           DT,2*MAXL+1)
                      DO 1000 MB=1,2*LB+1
                          DO 990 IX=1,3
                              DO 980 MA=1,2*LA+1
                                  DX(IBASA+MA,IBASB+MB,IX) = 
     .                                DX(IBASA+MA,IBASB+MB,IX) + 
     .                                CA*CB*DT(MA,IX,MB)
  980                         CONTINUE
  990                     CONTINUE
 1000                 CONTINUE
 2700             CONTINUE
 2800         CONTINUE
 2900     CONTINUE
 3000 CONTINUE
C
      CALL PSSX94(NORBA,NORBB,DX(1,1,1),(MAXL+1)**2)
      CALL PSSX94(NORBA,NORBB,DX(1,1,2),(MAXL+1)**2)
      CALL PSSX94(NORBA,NORBB,DX(1,1,3),(MAXL+1)**2)
C
      RETURN
      END
C
      SUBROUTINE PSNB3(IB,IC,RB,RC,DX)
C
C   Compute three-center one-electron integrals of the (B|A|C)
C   type necessary for H0B matrices ( (B|A|B) integrals can
C   also be handled by this code. (B|A|A) integrals can not).
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IB     - Index of the atom carrying left-hand orbitals.
C      IC     - Index of the atom carrying right-hand orbitals.
C      RB     - Coordinates of the left-hand atom relative to the
C               operator center, atomic units.
C      RC     - Coordinates of the right-hand atom relative to the
C               operator center, atomic units.
C      DX     - Output <\mu|r^-3 L|\nu> integrals, three-dimensional
C               array.
C
C   Accessed common blocks:
C
C      PSNPAR - NMR-specific parameters
C      ATOMS  - AO ranges
C
C   Local storage:
C
C      A litle bit...
C
C   Module logic:
C
C   Accuracy:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION RB(3), RC(3)
      DIMENSION DX((MAXL+1)**2,(MAXL+1)**2,3)
      DIMENSION DT(2*MAXL+1,3,2*MAXL+1)
C
C    Gather parameters
C
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPB = NAT(IB)
      NORBC = NLAST(IC)-NFIRST(IC)+1
      ITYPC = NAT(IC)
C
C    Compute integrals
C
      DO 800 LC=0,MAXLA(ITYPC)
          IBASC = LC**2
          DO 700 LB=0,MAXLA(ITYPB)
              IBASB = LB**2
              DO 100 MC=1,2*LC+1
                  DO 90 IX=1,3
                      DO 80 MB=1,2*LB+1
                          DX(IBASB+MB,IBASC+MC,IX) = ZERO
   80                 CONTINUE
   90             CONTINUE
  100         CONTINUE
              DO 600 JC=1,MAXN(LC,ITYPC)
                  NC = NN(JC,LC,ITYPC)
                  ZC = ZN(JC,LC,ITYPC)
                  CC = CN(JC,LC,ITYPC)
                  DO 500 JB=1,MAXN(LB,ITYPB)
                      NB = NN(JB,LB,ITYPB)
                      ZB = ZN(JB,LB,ITYPB)
                      CB = CN(JB,LB,ITYPB)
                      CALL PS3LR3(NB,LB,ZB,RB,NC,LC,ZC,RC,DT,2*MAXL+1)
                      DO 400 MC=1,2*LC+1
                          DO 300 IX=1,3
                              DO 200 MB=1,2*LB+1
                                  DX(IBASB+MB,IBASC+MC,IX) = 
     .                                DX(IBASB+MB,IBASC+MC,IX) 
     .                                + CB*CC*DT(MB,IX,MC)
  200                         CONTINUE
  300                     CONTINUE
  400                 CONTINUE
  500             CONTINUE
  600         CONTINUE
  700     CONTINUE
  800 CONTINUE
C
      CALL PSSX94(NORBB,NORBC,DX(1,1,1),(MAXL+1)**2)
      CALL PSSX94(NORBB,NORBC,DX(1,1,2),(MAXL+1)**2)
      CALL PSSX94(NORBB,NORBC,DX(1,1,3),(MAXL+1)**2)
C
      RETURN
      END
C
      SUBROUTINE PSNB3G(IGO,IB,IC,RB,RC,DX)
C
C   Compute three-center one-electron integrals of the (B|A|C)
C   type necessary for H0B matrices using Gaussian expansion.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IGO    - Order of the Gaussian expansion.
C      IB     - Index of the atom carrying left-hand orbitals.
C      IC     - Index of the atom carrying right-hand orbitals.
C      RB     - Coordinates of the left-hand atom relative to the
C               operator center, atomic units.
C      RC     - Coordinates of the right-hand atom relative to the
C               operator center, atomic units.
C      DX     - Output <\mu|r^-3 L|\nu> integrals, three-dimensional
C               array.
C
C   Accessed common blocks:
C
C      PSNPAR - NMR-specific parameters
C      ATOMS  - AO ranges
C
C   Local storage:
C
C   Module logic:
C
C      This is a version of PSNB3, see it for comments.
C
C   Accuracy:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      PARAMETER (ZERO=0.D0)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      DIMENSION RB(3), RC(3)
      DIMENSION DX((MAXL+1)**2,(MAXL+1)**2,3)
      DIMENSION DT(2*MAXL+1,3,2*MAXL+1)
C    Gather parameters
      NORBB = NLAST(IB)-NFIRST(IB)+1
      ITYPB = NAT(IB)
      NORBC = NLAST(IC)-NFIRST(IC)+1
      ITYPC = NAT(IC)
C    Zero output array
      DO 100 IX=1,3
        DO 90 I2=1,NORBC
          DO 80 I1=1,NORBB
              DX(I1,I2,IX) = ZERO
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
C
      DO 800 LC=0,MAXLA(ITYPC)
          IBASC = LC**2
          DO 700 LB=0,MAXLA(ITYPB)
              IBASB = LB**2
              DO 600 JC=1,MAXN(LC,ITYPC)
                  NC = NN(JC,LC,ITYPC)
                  ZC = ZN(JC,LC,ITYPC)
                  CC = CN(JC,LC,ITYPC)
                  DO 500 JB=1,MAXN(LB,ITYPB)
                      NB = NN(JB,LB,ITYPB)
                      ZB = ZN(JB,LB,ITYPB)
                      CB = CN(JB,LB,ITYPB)
                      CALL PSGLR3(IGO,NB,LB,ZB,RB,NC,LC,ZC,RC,
     .                            DT,2*MAXL+1)
                      SC = CB*CC
                      DO 400 MC=1,2*LC+1
                        DO 300 IX=1,3
                          DO 200 MB=1,2*LB+1
                              I1 = IBASB+MB
                              I2 = IBASC+MC
                              DX(I1,I2,IX) = DX(I1,I2,IX) 
     .                                     + SC*DT(MB,IX,MC)
  200                     CONTINUE
  300                   CONTINUE
  400                 CONTINUE
  500             CONTINUE
  600         CONTINUE
  700     CONTINUE
  800 CONTINUE
C
      CALL PSSX94(NORBB,NORBC,DX(1,1,1),(MAXL+1)**2)
      CALL PSSX94(NORBB,NORBC,DX(1,1,2),(MAXL+1)**2)
      CALL PSSX94(NORBB,NORBC,DX(1,1,3),(MAXL+1)**2)
C
      RETURN
      END
C
      SUBROUTINE PSNZR(IA,X,LDX)
C
C   Zero out magnetic interactions matrix elements.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Atom posessing the magnetic moment
C               IA.GT.NATOM means NICS, essentially no atom.
C      X      - Matrix to zero out
C      LDX    - Leading dimension of the X matrix
C
C   Accessed common blocks:
C
C      ATOMS  - Range of orbitals carried by an atom
C      PSDOPT - Computation options, in particular NMRLEV.
C      PSDGBL - Number of atoms, for consistency with derivatives code.
C      PSPRT  - Printing unit.
C
C   Local storage:
C
C   Module logic:
C
C      Depending on the NMRLEV, one, two or three-center terms are
C      zeroed. This subroutine is necessary to avoid needless O(N^4)
C      clearing steps then NMRLEV is not 3.
C
C   Bugs:
C
C      Code for NMRLEV=2 clears orbitals block centered at IA somewhat
C      too intensely ;-)
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      PARAMETER (ZERO=0.D0)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
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
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSDGBL/, /PSDOPT/, /PSPRTF/, /PSPRT /
      DIMENSION X(LDX,NORBS)
C
      IF( .NOT.LPNMR .AND. NMRLEV.EQ.1 ) THEN
          IF( IA.LE.NATOM ) THEN
              DO 200 M2=NFIRST(IA),NLAST(IA)
                  DO 190 M1=NFIRST(IA),NLAST(IA)
                      X(M1,M2) = ZERO
  190             CONTINUE
  200         CONTINUE
          ENDIF
      ELSE IF( .NOT.LPNMR .AND. NMRLEV.EQ.2 ) THEN
          DO 300 IB=1,NATOM
              DO 290 M2=NFIRST(IB),NLAST(IB)
                  DO 280 M1=NFIRST(IB),NLAST(IB)
                      X(M1,M2) = ZERO
  280             CONTINUE
  290         CONTINUE
  300     CONTINUE
          IF( IA.LE.NATOM ) THEN
              DO 400 M2=1,NORBS
                  DO 390 M1=NFIRST(IA),NLAST(IA)
                      X(M1,M2) = ZERO
  390             CONTINUE
  400         CONTINUE
              DO 500 M2=NFIRST(IA),NLAST(IA)
                  DO 490 M1=1,NORBS
                      X(M1,M2) = ZERO
  490             CONTINUE
  500         CONTINUE
          ENDIF
      ELSE IF( LPNMR .OR. NMRLEV.EQ.3 ) THEN
          DO 900 M2=1,NORBS
              DO 890 M1=1,NORBS
                  X(M1,M2) = ZERO
  890         CONTINUE
  900     CONTINUE
      ELSE
          WRITE(NB6,11000) NMRLEV
          STOP 'PSNZR'
      ENDIF
      RETURN
C
11000 FORMAT(' NMRLEV IS STRANGE (',I5,') IN PSNZR.')
      END
C
      DOUBLE PRECISION FUNCTION PSNML1(IA,H,LDH,RO)
C
C   Compute Tr(H RO) using one-center terms only. Both matrices
C   should be either symmetric or antisymmetric.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Atom posessing the magnetic moment.
C               IA.GT.NATOM means NICS
C      H      - A hamiltonial matrix derivative, upper half
C               of the matrix filled.
C      LDH    - Leading dimension of the H matrix
C      RO     - Density matrix in the packed form (upper half
C               is stored column-wise).
C
C   Accessed common blocks:
C
C      ATOMS  - Range of orbitals carried by an atom
C      PSDGBL - Number of atoms, for consistency with derivatives code.
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      PARAMETER (TWO =2.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ZERO=0.0D0)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      COMMON
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      SAVE /PSDGBL/
      DIMENSION H(LDH,NORBS), RO(*)
C
      IF( IA.GT.NATOM ) THEN
          PSNML1 = ZERO
          RETURN
      ENDIF
C
      IBASA = NFIRST(IA) - 1
      NORBA = NLAST(IA) - NFIRST(IA) + 1
      IRO   = ((IBASA+2)*(IBASA+1))/2
C
      SUM   = ZERO
      DO 200 M2=1,NORBA
          DO 100 M1=1,M2-1
              SUM = SUM + H(IBASA+M1,IBASA+M2) * RO(IRO)
              IRO = IRO + 1
  100     CONTINUE
          SUM = SUM + HALF * H(IBASA+M2,IBASA+M2) * RO(IRO)
          IRO = IRO + IBASA + 1
  200 CONTINUE
      PSNML1 = TWO * SUM
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION PSNML2(IA,H,LDH,RO)
C
C   Compute Tr(H RO) using one- and two-center terms only. Both matrices
C   should be either symmetric or antisymmetric.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Atom posessing the magnetic moment
C               IA.GT.NATOM means NICS.
C      H      - A hamiltonial matrix derivative, upper half
C               of the matrix filled.
C      LDH    - Leading dimension of the H matrix
C      RO     - Density matrix in the packed form (upper half
C               is stored column-wise).
C
C   Accessed common blocks:
C
C      ATOMS  - Range of orbitals carried by an atom
C      PSDGBL - Number of atoms, for consistency with derivatives code.
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      PARAMETER (TWO =2.0D0)
      PARAMETER (ZERO=0.0D0)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      COMMON
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      SAVE /PSDGBL/
      DIMENSION H(LDH,NORBS), RO(*)
C
C    First, sum all one-center terms
C
      SUM1 = ZERO
      DO 200 IB=1,NATOM
          SUM1 = SUM1 + PSNML1(IB,H,LDH,RO)
  200 CONTINUE
      IF( IA.GT.NATOM ) THEN
          PSNML2 = SUM1
          RETURN
      ENDIF
C
C    Then, add two-center terms involving atom A
C    All this hubub should be skipped if A is not an atom
C
      IBASA = NFIRST(IA) - 1
      ITOPA = NLAST(IA) + 1
      NORBA = NLAST(IA) - NFIRST(IA) + 1
C
C    ... column band from the top of the matrix to the A-centered orbitals...
C
      IRO   = ((IBASA+1)*IBASA)/2 + 1
      SUM2  = ZERO
      DO 300 M2=1,NORBA
          DO 290 M1=1,IBASA
              SUM2 = SUM2 + H(M1,IBASA+M2)*RO(IRO)
              IRO  = IRO + 1
  290     CONTINUE
          IRO = IRO + M2
  300 CONTINUE
C
C    ... row band right of the A-centered orbitals...
C
      IRO   = (ITOPA*(ITOPA-1))/2+1+IBASA
      DO 400 M2=ITOPA,NORBS
          DO 390 M1=1,NORBA
              SUM2 = SUM2 + H(IBASA+M1,M2)*RO(IRO)
              IRO  = IRO + 1
  390     CONTINUE
          IRO = IRO + (M2-NORBA)
  400 CONTINUE
C
      PSNML2 = SUM1 + TWO*SUM2
      RETURN
      END
C
      SUBROUTINE PSNMSY(IA,TAG,LBL1,LBL2,X,LDX)
C
C   Symmetrize NMR-related matrices.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Atom posessing the magnetic moment
C               IA.GT.NATOM is NICS.
C      TAG    - Name of the quantity (used to issue diagnotic/debugging
C               messages)
C      LBL1,
C      LBL2   - Two additional quantity identifiers.
C      X      - (Input/Outputi) Property to symmetrize
C      LDX    - Leading dimension of X
C
C   Accessed common blocks:
C
C      ATOMS  - Range of orbitals carried by an atom
C      PSDGBL - Number of atoms, for consistency with derivatives code.
C      PSDOPT - NMRLEV comes from here.
C
C   Local storage:
C
C   Module logic:
C
C      If maximum symmetrization deviation exceeds the one-, two- or
C      three-center threshold, warning is issued.
C
C   Bugs:
C
C      One-center parts are processed more than once for NMRLEV>1.
C      Likewise, two-center parts are processed twice for NMRLEV=3.
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH, VARINT
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (TOL1C=1D-8)
      PARAMETER (TOL2C=1D-6)
      PARAMETER (TOL3C=3D-4)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      COMMON
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSDOPT/ DSTORE, DSTEP,  DECONV, DPCONV, DCPHF, DPREC, DCDIFF,
     .         DSHIFT, DBASCR, DNCOFF,
     .         IUMIX,  IURHS,  IUK,    IURES,
     .         IPRINT, IENRG,  ICORE,  IDISK,  IMIX,   IDENS,
     .         INDSYM, IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC,
     .         IROWS,  IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP,
     .         IKMODE, ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT
     ./PSPRT / NB6
      SAVE /PSDGBL/, /PSDOPT/, /PSPRT /
      CHARACTER*(*) TAG, LBL1, LBL2
      DIMENSION X(LDX,NORBS)
C
      IF( NMRLEV.LT.1 .OR. NMRLEV.GT.3 ) THEN
          WRITE(NB6,11000) NMRLEV 
          STOP 'PSNMSY'
      ENDIF
      DEV1C  = ZERO
      DEV2C  = ZERO
      DEV3C  = ZERO
      VARINT = MOD(INTCTL,200).GE.100
C
      IF( NMRLEV.GE.1 .AND. IA.LE.NATOM ) THEN
          DO 1900 M2=NFIRST(IA),NLAST(IA)
              DO 1800 M1=NFIRST(IA),M2-1
                  AVE      = HALF*(X(M1,M2) + X(M2,M1))
                  XDEV     =   ABS(X(M1,M2) - X(M2,M1))
                  X(M1,M2) = AVE
                  X(M2,M1) = AVE
                  DEV1C    = MAX(DEV1C,XDEV)
 1800         CONTINUE
 1900     CONTINUE
      ENDIF
      IF( NMRLEV.GE.2 ) THEN
          DO 2200 IB=1,NATOM
              DO 2100 M2=NFIRST(IB),NLAST(IB)
                  DO 2000 M1=NFIRST(IB),M2-1
                      AVE      = HALF*(X(M1,M2) + X(M2,M1))
                      XDEV     =   ABS(X(M1,M2) - X(M2,M1))
                      X(M1,M2) = AVE
                      X(M2,M1) = AVE
                      DEV2C    = MAX(DEV2C,XDEV)
 2000             CONTINUE
 2100         CONTINUE
 2200     CONTINUE
          IF( IA.LE.NATOM ) THEN
              DO 2400 M2=NFIRST(IA),NLAST(IA)
                  DO 2300 M1=1,NFIRST(IA)-1
                      AVE      = HALF*(X(M1,M2) + X(M2,M1))
                      XDEV     =   ABS(X(M1,M2) - X(M2,M1))
                      X(M1,M2) = AVE
                      X(M2,M1) = AVE
                      DEV2C    = MAX(DEV2C,XDEV)
 2300             CONTINUE
 2400         CONTINUE
              DO 2600 M2=NLAST(IA)+1,NORBS
                  DO 2500 M1=NFIRST(IA),NLAST(IA)
                      AVE      = HALF*(X(M1,M2) + X(M2,M1))
                      XDEV     =   ABS(X(M1,M2) - X(M2,M1))
                      X(M1,M2) = AVE
                      X(M2,M1) = AVE
                      DEV2C    = MAX(DEV2C,XDEV)
 2500             CONTINUE
 2600         CONTINUE
          ENDIF
      ENDIF
      IF( NMRLEV.GE.3 .AND. .NOT.VARINT ) THEN
          DO 3900 M2=1,NORBS
              DO 3800 M1=1,M2-1
                  AVE      = HALF*(X(M1,M2) + X(M2,M1))
                  XDEV     =   ABS(X(M1,M2) - X(M2,M1))
                  X(M1,M2) = AVE
                  X(M2,M1) = AVE
                  DEV3C    = MAX(DEV3C,XDEV)
 3800         CONTINUE
 3900     CONTINUE
      ENDIF
C
      IF( IPRINT.GE.2 ) 
     .        WRITE(NB6,11100) LBL1, LBL2, TAG, DEV1C, DEV2C, DEV3C
      IF( DEV1C.GT.TOL1C .AND. IPRINT.GE.-1 ) 
     .        WRITE(NB6,11201) LBL1, LBL2, TAG, DEV1C
      IF( DEV2C.GT.TOL2C .AND. IPRINT.GE.-1 ) 
     .        WRITE(NB6,11202) LBL1, LBL2, TAG, DEV2C
      IF( DEV3C.GT.TOL3C .AND. IPRINT.GE.-1 ) 
     .        WRITE(NB6,11203) LBL1, LBL2, TAG, DEV3C
C
      RETURN
11000 FORMAT(' NMRLEV IS STRANGE (',I4,') IN PSNMSY')
11100 FORMAT(' SYMMETRIZATION DEVIATIONS FOR THE ',A,A,' COMPONENT OF ',
     .       A,' ARE:'/' 1C = ',G15.8,' 2C = ',G15.8,' 3C = ',G15.8)
11201 FORMAT(' WARNING: 1-CENTER SYMMETRIZATION DEVIATION FOR THE '
     .       ,A,A,' COMPONENT OF ',A,' IS VERY LARGE: ',G15.8)
11202 FORMAT(' WARNING: 2-CENTER SYMMETRIZATION DEVIATION FOR THE '
     .       ,A,A,' COMPONENT OF ',A,' IS VERY LARGE: ',G15.8)
11203 FORMAT(' WARNING: 3-CENTER SYMMETRIZATION DEVIATION FOR THE '
     .       ,A,A,' COMPONENT OF ',A,' IS VERY LARGE: ',G15.8)
      END
C
      SUBROUTINE PSNMAS(IA,TAG,LBL1,LBL2,X,LDX)
C
C   Antisymmetrize NMR-related matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA     - Atom posessing the magnetic moment
C               IA.GT.NATOM is NICS
C      TAG    - Name of the quantity (used to issue diagnotic/debugging
C               messages)
C      LBL1,
C      LBL2   - Two additional quantity identifiers.
C      X      - (Input/Outputi) Property to antisymmetrize
C      LDX    - Leading dimension of X
C
C   Accessed common blocks:
C
C      ATOMS  - Range of orbitals carried by an atom
C      PSDGBL - Number of atoms, for consistency with derivatives code.
C      PSDOPT - NMRLEV comes from here.
C      PSPRT  - Printing unit.
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH, VARINT
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (TOL1C=1D-8)
      PARAMETER (TOL2C=1D-6)
      PARAMETER (TOL3C=3D-4)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      COMMON
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSDOPT/ DSTORE, DSTEP,  DECONV, DPCONV, DCPHF, DPREC, DCDIFF,
     .         DSHIFT, DBASCR, DNCOFF,
     .         IUMIX,  IURHS,  IUK,    IURES,
     .         IPRINT, IENRG,  ICORE,  IDISK,  IMIX,   IDENS,
     .         INDSYM, IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC,
     .         IROWS,  IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP,
     .         IKMODE, ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT
     ./PSPRT / NB6
      SAVE /PSDGBL/, /PSDOPT/, /PSPRT /
      CHARACTER*(*) TAG, LBL1, LBL2
      DIMENSION X(LDX,NORBS)
C
      IF( NMRLEV.LT.1 .OR. NMRLEV.GT.3 ) THEN
          WRITE(NB6,11000) NMRLEV
          STOP 'PSNMAS'
      ENDIF
      DEV1C  = ZERO
      DEV2C  = ZERO
      DEV3C  = ZERO
      VARINT = MOD(INTCTL,200).GE.100
C
      IF( NMRLEV.GE.1 .AND. IA.LE.NATOM ) THEN
          DO 1900 M2=NFIRST(IA),NLAST(IA)
              DO 1800 M1=NFIRST(IA),M2-1
                  AVE      = HALF*(X(M1,M2) - X(M2,M1))
                  XDEV     =   ABS(X(M1,M2) + X(M2,M1))
                  X(M1,M2) = AVE
                  X(M2,M1) =-AVE
                  DEV1C    = MAX(DEV1C,XDEV)
 1800         CONTINUE
 1900     CONTINUE
      ENDIF
      IF( NMRLEV.GE.2 ) THEN
          DO 2200 IB=1,NATOM
              DO 2100 M2=NFIRST(IB),NLAST(IB)
                  DO 2000 M1=NFIRST(IB),M2-1
                      AVE      = HALF*(X(M1,M2) - X(M2,M1))
                      XDEV     =   ABS(X(M1,M2) + X(M2,M1))
                      X(M1,M2) = AVE
                      X(M2,M1) =-AVE
                      DEV2C    = MAX(DEV2C,XDEV)
 2000             CONTINUE
 2100         CONTINUE
 2200     CONTINUE
          IF( IA.LE.NATOM ) THEN
              DO 2400 M2=NFIRST(IA),NLAST(IA)
                  DO 2300 M1=1,NFIRST(IA)-1
                      AVE      = HALF*(X(M1,M2) - X(M2,M1))
                      XDEV     =   ABS(X(M1,M2) + X(M2,M1))
                      X(M1,M2) = AVE
                      X(M2,M1) =-AVE
                      DEV2C    = MAX(DEV2C,XDEV)
 2300             CONTINUE
 2400         CONTINUE
              DO 2600 M2=NLAST(IA)+1,NORBS
                  DO 2500 M1=NFIRST(IA),NLAST(IA)
                      AVE      = HALF*(X(M1,M2) - X(M2,M1))
                      XDEV     =   ABS(X(M1,M2) + X(M2,M1))
                      X(M1,M2) = AVE
                      X(M2,M1) =-AVE
                      DEV2C    = MAX(DEV2C,XDEV)
 2500             CONTINUE
 2600         CONTINUE
          ENDIF
      ENDIF
      IF( NMRLEV.GE.3 .AND. .NOT.VARINT ) THEN
          DO 3900 M2=1,NORBS
              DO 3800 M1=1,M2-1
                  AVE      = HALF*(X(M1,M2) - X(M2,M1))
                  XDEV     =   ABS(X(M1,M2) + X(M2,M1))
                  X(M1,M2) = AVE
                  X(M2,M1) =-AVE
                  DEV3C    = MAX(DEV3C,XDEV)
 3800         CONTINUE
 3900     CONTINUE
      ENDIF
C
      IF( IPRINT.GE.2 ) 
     .        WRITE(NB6,11100) LBL1, LBL2, TAG, DEV1C, DEV2C, DEV3C
      IF( DEV1C.GT.TOL1C .AND. IPRINT.GE.-1 ) 
     .        WRITE(NB6,11201) LBL1, LBL2, TAG, DEV1C
      IF( DEV2C.GT.TOL2C .AND. IPRINT.GE.-1 ) 
     .        WRITE(NB6,11202) LBL1, LBL2, TAG, DEV2C
      IF( DEV3C.GT.TOL3C .AND. IPRINT.GE.-1 ) 
     .        WRITE(NB6,11203) LBL1, LBL2, TAG, DEV3C
C
      RETURN
11000 FORMAT(' NMRLEV IS STRANGE (',I4,') IN PSNMAS')
11100 FORMAT(' ANTISYMMETRIZATION DEVIATIONS FOR THE ',A,A,
     .       ' COMPONENT OF ',
     .       A,' ARE:'/' 1C = ',G15.8,' 2C = ',G15.8,' 3C = ',G15.8)
11201 FORMAT(' WARNING: 1-CENTER ANTISYMMETRIZATION DEVIATION FOR THE '
     .       ,A,A,' COMPONENT OF ',A,' IS VERY LARGE: ',G15.8)
11202 FORMAT(' WARNING: 2-CENTER ANTISYMMETRIZATION DEVIATION FOR THE '
     .       ,A,A,' COMPONENT OF ',A,' IS VERY LARGE: ',G15.8)
11203 FORMAT(' WARNING: 3-CENTER ANTISYMMETRIZATION DEVIATION FOR THE '
     .       ,A,A,' COMPONENT OF ',A,' IS VERY LARGE: ',G15.8)
      END
C
      SUBROUTINE PSNMAT(TAG,LBL1,LBL2,X,LDX)
C
C   Antisymmetrize NMR-related matrix, for the quantities without
C   atoms posessing special properties.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      TAG    - Name of the quantity (used to issue diagnotic/debugging
C               messages)
C      LBL1,
C      LBL2   - Two additional quantity identifiers.
C      X      - (Input/Outputi) Property to antisymmetrize
C      LDX    - Leading dimension of X
C
C   Accessed common blocks:
C
C      ATOMS  - Range of orbitals carried by an atom
C      PSDGBL - Number of atoms, for consistency with derivatives code.
C      PSDOPT - NMRLEV comes from here.
C      PSPRT  - Printing unit.
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH, VARINT
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (TOL1C=1D-12)
      PARAMETER (TOL2C=1D-10)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      COMMON
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSDOPT/ DSTORE, DSTEP,  DECONV, DPCONV, DCPHF, DPREC, DCDIFF,
     .         DSHIFT, DBASCR, DNCOFF,
     .         IUMIX,  IURHS,  IUK,    IURES,
     .         IPRINT, IENRG,  ICORE,  IDISK,  IMIX,   IDENS,
     .         INDSYM, IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC,
     .         IROWS,  IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP,
     .         IKMODE, ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT
     ./PSPRT / NB6
      SAVE /PSDGBL/, /PSDOPT/, /PSPRT /
      CHARACTER*(*) TAG, LBL1, LBL2
      DIMENSION X(LDX,NORBS)
C
      IF( NMRLEV.LT.1 .OR. NMRLEV.GT.3 ) THEN
          WRITE(NB6,11000) NMRLEV
          STOP 'PSNMAT'
      ENDIF
      DEV1C  = ZERO
      DEV2C  = ZERO
      VARINT = MOD(INTCTL,200).GE.100
C
      DO 2200 IB=1,NATOM
          DO 2100 M2=NFIRST(IB),NLAST(IB)
              DO 2000 M1=NFIRST(IB),M2-1
                  AVE      = HALF*(X(M1,M2) - X(M2,M1))
                  XDEV     =   ABS(X(M1,M2) + X(M2,M1))
                  X(M1,M2) = AVE
                  X(M2,M1) =-AVE
                  DEV1C    = MAX(DEV1C,XDEV)
 2000         CONTINUE
 2100     CONTINUE
 2200 CONTINUE
C
      IF( .NOT.VARINT ) THEN
          DO 3900 M2=1,NORBS
              DO 3800 M1=1,M2-1
                  AVE      = HALF*(X(M1,M2) - X(M2,M1))
                  XDEV     =   ABS(X(M1,M2) + X(M2,M1))
                  X(M1,M2) = AVE
                  X(M2,M1) =-AVE
                  DEV2C    = MAX(DEV2C,XDEV)
 3800         CONTINUE
 3900     CONTINUE
      ENDIF
C
      IF( IPRINT.GE.2 ) 
     .        WRITE(NB6,11100) LBL1, LBL2, TAG, DEV1C, DEV2C
      IF( DEV1C.GT.TOL1C .AND. IPRINT.GE.-1 ) 
     .        WRITE(NB6,11201) LBL1, LBL2, TAG, DEV1C
      IF( DEV2C.GT.TOL2C .AND. IPRINT.GE.-1 ) 
     .        WRITE(NB6,11202) LBL1, LBL2, TAG, DEV2C
C
      RETURN
11000 FORMAT(' NMRLEV IS STRANGE (',I4,') IN PSNMAT')
11100 FORMAT(' ANTISYMMETRIZATION DEVIATIONS FOR THE ',A,A,
     .       ' COMPONENT OF ',
     .       A,' ARE:'/' 1C = ',G15.8,' 2C = ',G15.8)
11201 FORMAT(' WARNING: 1-CENTER ANTISYMMETRIZATION DEVIATION FOR THE '
     .       ,A,A,' COMPONENT OF ',A,' IS VERY LARGE: ',G15.8)
11202 FORMAT(' WARNING: 2-CENTER ANTISYMMETRIZATION DEVIATION FOR THE '
     .       ,A,A,' COMPONENT OF ',A,' IS VERY LARGE: ',G15.8)
      END
