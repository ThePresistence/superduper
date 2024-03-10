C     ******************************************************************
C
C     Code related to the half-electron derivatives.
C
C     ******************************************************************
      SUBROUTINE PSHFLO(C,LDC,NUMOPN,DUMP)
C
C   Fill matrices of contracted orbital coefficients needed to compute
C   the static part of the half-electron correction derivative.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      C      - Alpha orbital coefficients, positioned over first
C               active orbital.
C      LDC    - Leading dimension of the C matrix.
C      NUMOPN - Number of active orbitals.
C      DUMP   - Scratch array.
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C      PSDYNM - Dynamic memory handles
C      ATOMS  - Starting orbital number, number of orbitals
C               and atom types.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
C
      DIMENSION C(LDC,*), DUMP(*)
C
      IHLO   = IPSMOF(LHLO)
      NUMOIJ = (NUMOPN*(NUMOPN+1))/2
C
      DO 500 I=1,NATOM
          NFRST = NFIRST(I)
          NORB  = NLAST(I) - NFRST + 1
          NPAIR = ( NORB * (NORB+1) ) / 2
          CALL PSHSGO(I,NORB,NUMOPN,C(NFRST,1),LDC,DUMP(IHLO))
          IHLO  = IHLO + NPAIR * NUMOIJ
  500 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSHSGO(IATOM,NORB,NUMOPN,C,LDC,O)
C
C   Compute atom block of the matrices of contracted orbital coefficients.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IATOM  - Number of the atom to construct contracted coefficients on
C      NORB   - Number of AOs on the atom in question
C      NUMOPN - Number of open orbitals
C      C      - Alpha orbital coefficients positioned over atom block
C               on the first index and over open orbitals on the second.
C      LDC    - Leading dimension of the C matrix.
C      O      - Linear array for the contracted orbital coefficients.
C               should have at least ( NORB*(NORB+1)*NUMOPN*(NUMOPN+1) )/4
C               elements.
C
C   Accessed common blocks:
C
C      PSPRT  - Printing unit
C      PSPRTF - Debug output tuning flags
C
C   Modified common blocks:
C
C   Local storage:
C
C   Speedups possible:
C
C      Unrolling this code with specific values of NORB (which can be
C      1, 4 or 9 in the present version) might be a good idea.
C
C   Module logic:
C
C      Contracted orbital coefficients are computed according to the
C      expressions:
C
C         O_{\mu\mu}^{ij} = C_{\mu i}C_{\mu j} and
C         O_{\mu\nu}^{ij} = C_{\mu i}C_{\nu j} + C_{\mu j}C_{\nu i}
C
C      Coefficients are stored in the order \mu <= \nu on the inner
C      index, i <= j on the outer one.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      COMMON 
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSPRTF/, /PSPRT /
C
      DIMENSION C(LDC,*), O(*)
      IOUT = 0
C
      DO 900 J=1,NUMOPN
          DO 800 I=1,J
              DO 700 NU=1,NORB
                  DO 600 MU=1,NU-1
                      O(IOUT+MU) = C(MU,I)*C(NU,J) + C(MU,J)*C(NU,I)
  600             CONTINUE
                  O(IOUT+NU) = C(NU,I)*C(NU,J)
                  IOUT       = IOUT + NU
  700         CONTINUE
  800     CONTINUE
  900 CONTINUE
C
      IF(LPHALF) THEN
          IPAIRS = ( NORB*(NORB+1) ) / 2
          IOUT   = 1
          WRITE(NB6,11000) IATOM
          DO 1400 J=1,NUMOPN
              DO 1300 I=1,J
                  WRITE(NB6,11010) I, J
                  CALL PSDPPM(NORB,O(IOUT))
                  IOUT = IOUT + IPAIRS
 1300         CONTINUE
 1400     CONTINUE
      ENDIF
C
      RETURN
11000 FORMAT(' CONTRACTED ORBITAL COEFFICIENT MATRICES FOR ATOM ', 
     .       I4, ' ARE:')
11010 FORMAT(' FOR THE ', I3, ',', I3, ' PAIR:')
      END
C
      SUBROUTINE PSHSJK(NPAIRA,NPAIRB,AINTS,LDA,OIJA,OIJB,XJ,XK,LDX)
C
C   Compute pair contribution to the static part of the derivatives
C   of two-electron integrals in the MO basis.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NPAIRA - Number of orbital pairs on the first atom
C      NPAIRB - Number of orbital pairs on the second atom
C      AINTS  - Integral derivatives, A orbital pair number is
C               the first index, matrix should be positioned
C               so that AINTS(1,1) correspond to (ss,ss) integral
C      LDA    - Leading dimension of the AINTS matrix
C      OIJA   - Contracted orbital coefficients on the A atom
C      OIJB   - Contracted orbital coefficients on the B atom
C      XJ     - Contributions to the Jij integral derivatives
C               (elements with i<=j are filled)
C      XK     - Contributions to the Kij integral derivatives
C               (elements with i<j are filled)
C      LDX    - Leading dimension of the XJ and XK matrices
C
C   Accessed common blocks:
C
C      PSHALF - Half-electron parameters
C      PSPRT  - Printing unit
C
C   Modified common blocks:
C
C   Local storage:
C
C      45*MAXACT (2025 in the present version) DOUBLE PRECISION cells
C      used for contractions on the second index.
C
C   Module logic:
C
C      Pair contributions to J' and K' are computed according
C      to the expressions:
C
C      XJij =     Sum   Oii(mu,nu)   Sum    Ojj(la,si) (mu nu,la si)' +
C              mu <= nu            la <= si
C
C           +     Sum   Ojj(mu,nu)   Sum    Oii(la,si) (mu nu,la si)' +
C              mu <= nu            la <= si
C
C      XKij = 2 * Sum   Oij(mu,nu)   Sum    Oij(la,si) (mu nu,la si)'
C              mu <= nu            la <= si
C
C      This is accomplished in two stages. First, entire matrix of
C      derivative integrals is contracted by the second coefficient
C      with the fragment of Oij belonging to B, giving TEMP(mu nu,i j).
C      Then, elements with i != j are contracted with Oij on atom A,
C      giving contributions to XK. Elements with i = j are contracted
C      with all k = l, producing contributions to XJ.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
C
      PARAMETER (TWO=2.0D0)
      PARAMETER (ONE=1.0D0)
      PARAMETER (ZERO=0.0D0)
C
      COMMON 
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN),
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN),
     .         IACTI(MAXACT), IACTJ(MAXACT)
     ./PSPRT / NB6
      SAVE /PSHALF/, /PSPRT /
C
      DIMENSION AINTS(LDA,*), OIJA(NPAIRA,*), OIJB(NPAIRB,*)
      DIMENSION XJ(LDX,*), XK(LDX,*)
      DIMENSION TEMP(45,MAXACT)
C
      IF( NPAIRA.GT.45 .OR. NUMOIJ.GT.MAXACT ) THEN
          WRITE(NB6,10000) NPAIRA, NUMOIJ, 45, MAXACT
          STOP 'PSHSJK'
      ENDIF
C
C     Contract by coefficients at the second atom:
C
C            NUMOIJ        NPAIRB
C          <------->    <---------->         NUMOIJ
C         ^            ^ ..........        <------->
C         |            |                  ^ .
C         |            |                  | .
C         |            |                  | .
C  NPAIRA |   TEMP   = |    AINTS      X  | . OIJB
C         |            |                  | .
C         |            |                  | .
C         |            |                  v .
C         v            v
C
      CALL PSZRMB(NPAIRA,NUMOIJ,TEMP,45)
      CALL DGEMM('N','N',NPAIRA,NUMOIJ,NPAIRB,ONE,AINTS,LDA,OIJB,
     .           NPAIRB,ZERO,TEMP,45)
C
C     Contract by coefficients at the first atom.
C
      IJ = 1
      DO 800 J=1,NUMOPN
C
C         (i,j) need only (i,j) to contribute to XKij
C
          DO 400 I=1,J-1
              XK(I,J) = TWO * DDOT(NPAIRA,TEMP(1,IJ),1,OIJA(1,IJ),1)
              IJ = IJ + 1
  400     CONTINUE
C
C         (i,i) elements have to be contracted with all (j,j) ones
C               (producing contributions to XJij)
C
          II = 1
          DO 500 I=1,J-1
C
C             XJ were already initialized by the DO-600 loop.
C             "used before set" warning issued by ftnchek is incorrect.
C
              XJ(I,J) = XJ(I,J) + 
     .                  DDOT(NPAIRA,TEMP(1,IJ),1,OIJA(1,II),1)
              II = II + I + 1
  500     CONTINUE
C
C         XJ(I,I) previously used here is equivalent to XJ(J,J), since
C         I is equal to J after the previous DO statement was executed
C         (even if the iteration count was zero), but the later form
C         is clearer...
C
          XJ(J,J) = TWO * DDOT(NPAIRA,TEMP(1,IJ),1,OIJA(1,II),1)
          II = II + J + 1
          DO 600 I=J+1,NUMOPN
              XJ(J,I) = DDOT(NPAIRA,TEMP(1,IJ),1,OIJA(1,II),1)
              II = II + I + 1
  600     CONTINUE
          IJ = IJ + 1
  800 CONTINUE
C
      RETURN
10000 FORMAT(' STATIC TEMPORARY OVERFLOW IN PSHSJK. NEED ',
     .       I4, ',', I4, ', HAVE ', I4, ',', I4 )
      END
C
      SUBROUTINE PSHLQ(CALP,LDC,EALP,DUMP)
C
C   Compute contraction coefficients for the half-electron
C   and CI derivatives.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients
C      LDC    - Leading dimension of CALP
C      EALP   - Alpha orbital energies
C      DUMP   - Base of the dynamic memory pool
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options 
C      PSDGBL - Global computation parameters
C      PSDGB2 - Some of the parameters relevant for
C               "response" computation.
C      PSPRT  - Printing unit
C      PSPRTF - Debug output tuning flags
C      PSHALF - Half-electron parameters
C      PSDYNM - Dynamic memory handles
C      PSDCIP - CI parameters
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
      PARAMETER (MAXCIO=LMACT)
C
      PARAMETER (ONE=1.0D0)
C
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
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
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN),
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN),
     .         IACTI(MAXACT), IACTJ(MAXACT)
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
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSPRTF/, /PSHALF/, /PSDCIP/
      SAVE /PSPRT /
C
      DIMENSION CALP(LDC,*), EALP(*), DUMP(*)
C    Interesting parameters
      IF(HALFEL) THEN
          NACT   = NUMOPN
      ELSE IF(DOCI) THEN
          NACT   = NCIO
          ICIDLT = IPSMOF(LCIDLT)
      ELSE
          WRITE(NB6,11000) 
          STOP 'PSHLQ'
      ENDIF
C    First, compute redundant Q-vector, QEX
      IHLQ   = IPSMOF(LHLQ)
      IHLQEX = IPSMOF(LHLQEX)
      IF(HALFEL) CALL PSHLQE(CALP,LDC,DUMP(IHLQEX),NORBS,DUMP)
      IF(DOCI)   CALL PSCIQE(CALP,LDC,DUMP(IHLQEX),NORBS,DUMP)
      IF( LPHALF.OR.LPCI ) THEN
          WRITE(NB6,11010)
          CALL PSDPGM(NORBS,NACT,DUMP(IHLQEX),NORBS)
      ENDIF
C    Wrap QEX onto minimal set of non-redundant and redundant 
C    variables
      IF(HALFEL) CALL PSHLQW(DUMP(IHLQEX),NORBS,DUMP(IHLQ))
      IF(DOCI)   CALL PSCIQW(DUMP(IHLQEX),NORBS,DUMP(IHLQ),DUMP(ICIDLT))
C
      IF( LPHALF.OR.LPCI ) THEN
          WRITE(NB6,11100)
          CALL PSXPRT(DUMP(IHLQ))
          WRITE(NB6,11110)
          IF(HALFEL) WRITE(NB6,11120) 
     .            (IACTI(I),IACTJ(I),DUMP(IHLQ+IQSZ+I-1),I=1,NMRACT)
          IF(DOCI) CALL PSDPGM(NORBS,NCIO,DUMP(IHLQ+IQSZ),NORBS)
      ENDIF
C    Convert contraction coefficients into right-hand side of
C    non-redundant response problem. 
      CALL PSDUNO(DUMP(IHLQ))
      IF(HALFEL) CALL PSHLUG(DUMP(IHLQ+IQSZ),EALP)
      IF(DOCI)   CALL PSCIUG(DUMP(IHLQ+IQSZ),EALP)
      IF( IQSWAP.LT.0 ) THEN
          IQTMP  = IPSMOF(LQ) - ICPV1 * IQSZ
      ELSE
          IQTMP  = IPSMOF(LQTMP)
      ENDIF
      IF( IHLWRP.EQ.1 ) THEN
          CALL PSHLWX(DUMP(IHLQ+IQSZ),DUMP(IQTMP),CALP,LDC,DUMP)
      ELSE
          CALL PSHLWA(DUMP(IHLQ+IQSZ),DUMP(IQTMP),CALP,LDC,DUMP)
      ENDIF
      CALL DAXPY(IQSZ,ONE,DUMP(IHLQ),1,DUMP(IQTMP),1)
      IF(LPHALF) THEN
          WRITE(NB6,11200)
          CALL PSXPRT(DUMP(IQTMP))
      ENDIF
      IF( IQSWAP.GT.0 ) THEN
          CALL PSDWD(IURHS,1-ICPV1,DUMP(IQTMP),IQSZ)
      ENDIF
C
      RETURN
11000 FORMAT(' PSHLQ CALLED WITH NEITHER HALFEL NOR DOCI SET.')
11010 FORMAT(' REDUNDANT CONTRACTION VECTOR (QEX):' )
11100 FORMAT(' NON-REDUNDANT CONTRACTION VECTOR (Q):' )
11110 FORMAT(' ADDITIONAL REDUNDANT CONTRACTIONS:' )
11120 FORMAT(3(1X,I3,',',I3,': ',F10.6,5X))
11200 FORMAT(' RIGHT-HAND SIDE FOR THE Z-VECTOR:' )
      END
C
      SUBROUTINE PSDUNO(X)
C
C   Scale non-redundant CPHF variables by inverce difference
C   in occupation numbers.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      X      - Variables to scale
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation options.
C      PSDGB2 - Global computation options needed only if response
C      PSOCC  - Occupation number groups description
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
      SAVE /PSOCC/, /PSDGBL/, /PSDGB2/
C
      DIMENSION X(*)
C
      CALL PSDUNP(X,NMGRPA,IGCNTA,DOCCA,IGBASA,MAXGRP)
      IF(UHF) THEN
          CALL PSDUNP(X(IQSZA+1),NMGRPB,IGCNTB,DOCCB,IGBASB,MAXGRP)
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSDUNP(X,NMGRP,IGCNT,DOCC,IGBAS,LDI)
C
C   Scale non-redundant CPHF variables for the single spin 
C   by inverce difference in occupation numbers.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      X      - Variables to scale
C      NMGRP  - Number of occupation number groups
C      IGCNT  - Sizes of the occupation number groups
C      DOCC   - Occupation numbers
C      IGBAS  - Starting indices of orbital blocks
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
      PARAMETER (ONE=1.0D0)
C
      DIMENSION X(*), IGCNT(NMGRP), DOCC(NMGRP), IGBAS(LDI,*)
C
      DO 500 I1=1,NMGRP-1
          N1 = IGCNT(I1)
          O1 = DOCC(I1)
          DO 300 I2=I1+1,NMGRP
              N2  = IGCNT(I2)
              O2  = DOCC(I2)
              ISZ = N1*N2
              SCL = ONE/(O1-O2)
              IF( SCL.NE.ONE ) CALL DSCAL(ISZ,SCL,X(IGBAS(I1,I2)),1)
  300     CONTINUE
  500 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSHLUG(XR,EALP)
C
C   Scale redundant CPHF variables by the difference in orbital
C   energies.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      XR     - Variables to scale, only variables listed in
C               IACTI/IACTJ are present
C      EALP   - Orbital energies, EVs
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options
C      PSOCC  - Occupation number groups description
C      PSHALF - Half-electron correction description
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
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
C
      PARAMETER (AU   = 3.67511944138184490995957368D-2)
      PARAMETER (TOL  = 1.0D-4)
      PARAMETER (ZERO = 0.0D0)
C
      COMMON 
     ./PSDOPT/ DSTORE, DSTEP,  DECONV, DPCONV, DCPHF, DPREC, DCDIFF,
     .         DSHIFT, DBASCR, DNCOFF,
     .         IUMIX,  IURHS,  IUK,    IURES,
     .         IPRINT, IENRG,  ICORE,  IDISK,  IMIX,   IDENS,
     .         INDSYM, IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC,
     .         IROWS,  IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP,
     .         IKMODE, ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN),
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN),
     .         IACTI(MAXACT), IACTJ(MAXACT)
     ./PSPRT / NB6
      SAVE /PSDOPT/, /PSHALF/, /PSPRT /
C
      DIMENSION XR(*), EALP(*)
C
      DO 700 K=1,NMRACT
          I  = IACTI(K) + IOPEN1 - 1
          J  = IACTJ(K) + IOPEN1 - 1
          DE = EALP(I) - EALP(J)
          IF( ABS(DE).LE.TOL ) THEN
C
C             Bad Thing (tm) happened - we are required to determine
C             coupling coefficient for (nearly) degenerate orbitals,
C             which is numerically unstable. All we can do is to
C             zero out contribution and try to continue...
C
              IF( IPRINT.GT.-5 ) THEN
                  WRITE(NB6,10000) I, J, TOL
              ENDIF
              XR(K) = ZERO
          ELSE
      
          ENDIF
          XR(K) = XR(K) / (AU*DE)
  700 CONTINUE
C
      RETURN
10000 FORMAT(' WARNING: ENERGY SEPARATION FOR ACTIVE REDUNDANT PAIR ',
     .       I4, ',', I4, ' IS LESS THAN ', G10.5, ' EV'/
     .       ' ATTEMPTING TO EXCLUDE PAIR AND CONTINUE - DERIVATIVES ',
     .       'ARE IN DOUBT.' )
      END
C
      SUBROUTINE PSHLQE(CALP,LDC,QEX,LDQ,DUMP)
C
C   Compute redundant contraction coefficients for 
C   the half-electron derivatives.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients
C      LDC    - Leading dimension of CALP
C      QEX    - Place for contraction coefficients
C      LDQ    - Leading dimension of QEX
C      DUMP   - Base of the dynamic memory pool
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C      PSHALF - Half-electron parameters
C      PSDYNM - Dynamic memory handles
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      QEX is computed according to the expression:
C
C                  open                open
C      QEX(p,i) = 2 Sum Hij (jj,ip) + 2 Sum Gij (ji,jp)
C                    j                   j
C
c      Actual code is a litle bit complicated by the desire
C      to use integral symmetry on the first and second
C      indices.
C
C   Bugs:
C
C      Only symmetry on the first two indices is used.
C      Hence, some integrals are computed more then once.
C      Hij and Gij are expected to be symmetrical.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
C
      PARAMETER (TWO=2.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN),
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN),
     .         IACTI(MAXACT), IACTJ(MAXACT)
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
      SAVE /PSDGBL/, /PSHALF/
C
      DIMENSION CALP(LDC,*), QEX(LDQ,*), DUMP(*)
C
C     DUMP(IINTS) will contain computed MO integrals
C
      IINTS = IPSMOF(LAI2T5)
C
C     Zero output array
C
      CALL PSZRM(NORBS,NUMOPN,QEX,LDQ)
C
      DO 1000 I=1,NUMOPN
          CALL PST1ST(DUMP,CALP(1,IOPEN1+I-1))
          DO 100 J=1,I
              CALL PST2ND(DUMP,CALP(1,IOPEN1+J-1),J)
  100     CONTINUE
C
C         Contributions from Hij, integrals are in the form (jj,ip).
C
          CALL PST3RD(DUMP,I,CALP(1,IOPEN1),LDC,NUMOPN)
          CALL PST4TH(DUMP,CALP,LDC,NUMOPN,NORBS)
          DO 200 K=1,NUMOPN
              CALL DAXPY(NORBS,TWO*HLFH(I,K),DUMP(IINTS+NORBS*(K-1)),1,
     .                   QEX(1,K),1)
  200     CONTINUE
C
C         Diagonal contributions from Gij, integrals are in the form (ii,ip)
C
          CALL DAXPY(NORBS,TWO*HLFG(I,I),DUMP(IINTS+NORBS*(I-1)),1,
     .               QEX(1,I),1)
C
C         Contributions from Gij, integrals are in form (ij,jp)
C         These will contribute both to QEX(*,I) and to QEX(*,J)
C
          DO 300 J=1,I-1
C
C             (ij,ip)
C
              CALL PST3RD(DUMP,J,CALP(1,IOPEN1+I-1),LDC,1)
              CALL PST4TH(DUMP,CALP,LDC,1,NORBS)
              CALL DAXPY(NORBS,TWO*HLFG(I,J),DUMP(IINTS),1,QEX(1,J),1)
C
C             (ji,jp)
C
              CALL PST3RD(DUMP,J,CALP(1,IOPEN1+J-1),LDC,1)
              CALL PST4TH(DUMP,CALP,LDC,1,NORBS)
              CALL DAXPY(NORBS,TWO*HLFG(I,J),DUMP(IINTS),1,QEX(1,I),1)
  300     CONTINUE
 1000 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSHLQW(QEX,LDQ,Q)
C
C   Convert redundant set of half-electron 
C   contraction coefficients into minimal one.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      QEX    - Redundant set of contraction coefficients
C      LDQ    - Leading dimension of the QEX matrix
C      Q      - Space for non-redundant coefficients
C
C   Accessed common blocks:
C
C      PSDGB2 - Response-related parameters
C      PSOCC  - Occupation blocks description
C      PSHALF - Half-electron parameters
C      PSDCIP - CI parameters
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Before contraction coefficients could be used for
C      solving response equations for Z-vector, they have
C      to be converted into standard active variable blocks
C      format. This is accomplished as following:
C
C      1. For closed-open blocks, QEX is copied unchanged
C      2. For open-virtual blocks, QEX is copied with sign reversed.
C      3. For open-open block, Q(i,j) is computed as QEX(i,j) - QEX(j,i)
C
C      Additionally, contraction coefficients for active
C      redundant CPHF variables, which are not covered
C      by basic format, are included (in the order
C      specified by IACTI/IACTJ) at the end of Q.
C
C   Bugs:
C
C      UHF is not supported.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXGRP=10)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
C
      PARAMETER (ONE=1.0D0)
      PARAMETER (SONE=-1.0D0)
C
      COMMON 
     ./PSDGB2/ NPAIR, NPAIR2, NOCCA, NVACA, NOCCB, NVACB,
     .         IQSZA, IQSZB,  IQSZ,  ICPV1, ICPVL, NCPVRS
     ./PSOCC / DOCCA(MAXGRP), DOCCB(MAXGRP), NMGRPA, NMGRPB,
     .         IOPN1A, IOPN1B, IOPNLA, IOPNLB, IMAXBL,
     .         IG1STA(MAXGRP), IGCNTA(MAXGRP),
     .         IG1STB(MAXGRP), IGCNTB(MAXGRP),
     .         IGBASA(MAXGRP,MAXGRP), IGBASB(MAXGRP,MAXGRP)
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN),
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN),
     .         IACTI(MAXACT), IACTJ(MAXACT)
      SAVE /PSHALF/, /PSOCC/, /PSDGB2/
C
      DIMENSION QEX(LDQ,*), Q(*)
C
C     Zero out Q, since most elements will never be touched otherwise.
C
      CALL PSZRV(IQSZ,Q)
      IOPENL = IOPEN1+NUMOPN-1
C
C     Coefficients for non-redundant CPHF variables
C
      DO 600 I1=1,NMGRPA-1
          IND1 = IG1STA(I1)
          N1   = IGCNTA(I1)
          DO 500 I2=I1+1,NMGRPA
              IND2 = IG1STA(I2)
              N2   = IGCNTA(I2)
C
C             It is well-known (to me ;-) that IND1 is always .LT. IND2
C             Therefore, the following formal test:
C             IF( IND1.LE.IOPENL .AND. IND2.GE.IOPEN1 .AND. IND2.LE.IOPENL ) THEN
C             can actually be written as:
C
              IF( IND2.GE.IOPEN1 .AND. IND2.LE.IOPENL ) THEN
C
C                 (!virtual)-open pair, Q(open,!virtual) += QEX(!virtual,open)
C                 !Virtual: J1(IND1)
C                 Open:     J2(IND2)
C
                  IQX = IGBASA(I1,I2)
                  DO 100 J1=1,N1
                      CALL DAXPY(N2,ONE,QEX(IND1+J1-1,IND2+1-IOPEN1),
     .                           LDQ,Q(IQX),1)
                      IQX = IQX + N2
  100             CONTINUE
              ENDIF
C
C             IF( IND1.GE.IOPEN1 .AND. IND1.LE.IOPENL .AND. IND2.GE.IOPEN1 ) THEN
C             Once again, IND2 is always .GT. IND1, hence test above reduces to:
C
              IF( IND1.GE.IOPEN1 .AND. IND1.LE.IOPENL ) THEN
C
C                 open-(!closed) pair, Q(!closed,open) -= QEX(!closed,open)
C                 Open:     J1(IND1)
C                 !Closed:  J2(IND2)
C
                  IQX = IGBASA(I1,I2)
                  DO 200 J1=1,N1
                      CALL DAXPY(N2,SONE,QEX(IND2,IND1+J1-IOPEN1),1,
     .                           Q(IQX),1)
                      IQX = IQX + N2
  200             CONTINUE
              ENDIF
  500     CONTINUE
  600 CONTINUE
C
C     And now for active redundant CPHF variables
C
      IQX = IQSZ + 1
      DO 900 IJ=1,NMRACT
          I      = IACTI(IJ)
          J      = IACTJ(IJ)
          Q(IQX) = QEX(IOPEN1+I-1,J) - QEX(IOPEN1+J-1,I)
          IQX    = IQX + 1
  900 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSHLWX(QRED,QADD,CALP,LDC,DUMP)
C
C   Distribute contributions to the half-electron correction
C   or CI derivative from active redundant CPHF variables onto 
C   non-redundant variables using explicit formation of
C   Wr matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      QRED   - Contraction coefficients for active redundant
C               CPHF variables, in the order specified by
C               IACTI/IACTJ. Coefficients should be already
C               scaled by inverse gamma matrix.
C      QADD   - Coefficients of non-redundant varibles, filled
C               on return.
C      CALP   - Alpha orbital coefficients
C      LDC    - Leading dimension of CALP
C      DUMP   - Base of the synamic memory pool
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options
C      PSDGBL - System parameters
C      PSDGB2 - Response-related parameters
C      PSOCC  - Occupation blocks description
C      PSHALF - Half-electron parameters
C      PSDYNM - Dynamic memory handles
C      PSDCIP - CI parameters.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Since integrals transformation is the most expensive
C      part of the process, programming objective was to
C      minimize number of redundant integral transformation.
C
C   Bugs:
C
C      UHF is not supported.
C
C      Due to a lot of conditionals, code is not particularly 
C      efficient.
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
      PARAMETER (MAXGRP=10)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
C
      PARAMETER (ZERO=0.0D0)
      PARAMETER (SONE=-1.0D0)
      PARAMETER (FOUR=4.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL HAVENZ
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
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN),
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN),
     .         IACTI(MAXACT), IACTJ(MAXACT)
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
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
      SAVE /PSHALF/, /PSOCC/, /PSDGBL/, /PSDGB2/, /PSDCIP/
C    Stracture of QRED is different depending on computation mode.
C    In half-electron mode, it is a vector of the dimension NMRACT.
C    In CI mode, it is a matrix of the dimension (NORBS,NCIO).
      DIMENSION QRED(*), QADD(*), CALP(LDC,*), DUMP(*)
      IPACK(I,J) = (I-1)*NORBS + J
C    DUMP(IINTS) will contain computed MO integrals
      IINTS = IPSMOF(LAI2T5)
C    Zero out additional dependence coefficients
      CALL PSZRV(IQSZ,QADD)
C    Initialization is slightly different for the CI and half-electron cases
      IF(HALFEL) THEN
          NACT   = NUMOPN
          JMAX   = NUMOPN
      ENDIF
      IF(DOCI) THEN
          NACT   = NCIO
          JMIN   = 1
          JMAX   = NORBS
      ENDIF
C    Loop over redundant coefficients. 
      DO 1000 I=1,NACT
          IF(HALFEL) THEN
              JMIN  = I + 1
              IREAL = IOPEN1 + (I-1)
          ENDIF
          IF(DOCI) THEN
              IREAL = INDCIO(I)
          ENDIF
C        Transformation have to be performed only then
C        there is an active redundant varible with such
C        first index. The actual check is different for
C        half-electron and CI cases.
          HAVENZ = .FALSE.
          IF(HALFEL) THEN
              DO 100 J=JMIN,JMAX
                  IF( IACTID(I,J).GT.0 ) THEN
                      IF( QRED(IACTID(I,J)).NE.ZERO ) HAVENZ = .TRUE.
                  ENDIF
  100         CONTINUE
          ENDIF
          IF(DOCI) THEN
              DO 120 J=JMIN,JMAX
                  IF( QRED(IPACK(I,J)).NE.ZERO ) HAVENZ = .TRUE.
  120         CONTINUE
          ENDIF
          IF( .NOT.HAVENZ ) GOTO 1000
C        Do transformation over first two indices.
          CALL PST1ST(DUMP,CALP(1,IREAL))
          DO 180 JREAL=1,NORBS
              CALL PST2ND(DUMP,CALP(1,JREAL),JREAL)
  180     CONTINUE
C
          DO 900 J=JMIN,JMAX
              IF(HALFEL) THEN
                  IJ    = IACTID(I,J)
                  JREAL = IOPEN1 + (J-1)
                  IF( IJ.LE.0 ) GOTO 900
              ENDIF
              IF(DOCI) THEN
                  IJ    = IPACK(I,J)
                  JREAL = J
              ENDIF
              SCALE = QRED(IJ)
              IF( SCALE.EQ.ZERO ) GOTO 900
C            (ij,kl) contributions first. Although we can't
C            see it (;-)) QADD is actually accessed linearly.
              DO 300 I1=1,NMGRPA-1
                  IND1 = IG1STA(I1)
                  N1   = IGCNTA(I1)
                  CALL PST3RD(DUMP,JREAL,CALP(1,IND1),LDC,N1)
                  DO 290 I2=I1+1,NMGRPA
                      IND2 = IG1STA(I2)
                      N2   = IGCNTA(I2)
                      CALL PST4TH(DUMP,CALP(1,IND2),LDC,N1,N2)
                      CALL DAXPY(N1*N2,FOUR*SCALE,DUMP(IINTS),1,
     .                           QADD(IGBASA(I1,I2)),1)
  290             CONTINUE
  300         CONTINUE
C            (ik,jl) and (il,jk) are a little tricky, since 
C            third-index transforms for all but first and last
C            blocks can be shared. 
              DO 600 I1=1,NMGRPA
                  IND1 = IG1STA(I1)
                  N1   = IGCNTA(I1)
                  DO 500 KREAL=IND1,IND1+N1-1
                      CALL PST3RD(DUMP,KREAL,CALP(1,JREAL),LDC,1)
C                    (il,jk) contribution to the (I2,I1) block
                      DO 350 I2=1,I1-1
                          IND2  = IG1STA(I2)
                          N2    = IGCNTA(I2)
                          IQADD = IGBASA(I1,I2) + (KREAL-IND1)
                          CALL PST4TH(DUMP,CALP(1,IND2),LDC,1,N2)
                          CALL DAXPY(N2,SONE*SCALE,DUMP(IINTS),1,
     .                               QADD(IQADD),N1)
  350                 CONTINUE
C                    (ik,jl) contribution to the (I1,I2) block
                      DO 400 I2=I1+1,NMGRPA
                          IND2  = IG1STA(I2)
                          N2    = IGCNTA(I2)
                          IQADD = IGBASA(I1,I2) + (KREAL-IND1)*N2
                          CALL PST4TH(DUMP,CALP(1,IND2),LDC,1,N2)
                          CALL DAXPY(N2,SONE*SCALE,DUMP(IINTS),1,
     .                               QADD(IQADD),1)
  400                 CONTINUE
  500             CONTINUE
  600         CONTINUE
C            Done with this (i,j) pair.
  900     CONTINUE
 1000 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSHLWA(QRED,QADD,CALP,LDC,DUMP)
C
C   Distribute contributions to the half-electron correction
C   and CI derivative from active redundant CPHF variables onto 
C   non-redundant variables using AO basis construction of 
C   Wx matrix
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      QRED   - Contraction coefficients for active redundant
C               CPHF variables, in the order specified by
C               IACTI/IACTJ. Coefficients should be already
C               scaled by inverse gamma matrix.
C      QADD   - Coefficients of non-redundant varibles, filled
C               on return.
C      CALP   - Alpha orbital coefficients
C      LDC    - Leading dimension of CALP
C      DUMP   - Base of the synamic memory pool
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options
C      PSDGBL - System parameters
C      PSDGB2 - Response-related parameters
C      PSHALF - Half-electron parameters
C      PSDYNM - Dynamic memory handles
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      This module is a variation on the subject of PSKROW.
C
C   Bugs:
C
C      UHF is not supported.
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
C
      PARAMETER (ZERO=0.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, HAVERD, DODIP, LIMAG, DOPTCH
      LOGICAL HAVENZ
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSDGB2/ NPAIR, NPAIR2, NOCCA, NVACA, NOCCB, NVACB,
     .         IQSZA, IQSZB,  IQSZ,  ICPV1, ICPVL, NCPVRS
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN),
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN),
     .         IACTI(MAXACT), IACTJ(MAXACT)
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
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
      SAVE /PSHALF/, /PSDGBL/, /PSDGB2/, /PSDCIP/
C
      DIMENSION QRED(*), QADD(*), CALP(LDC,*), DUMP(*)
      IPACK(I,J) = (I-1)*NORBS + J
C    Dynamic memory items
      ICPAYA = IPSMOF(LCPAYA)
      ICPAY  = IPSMOF(LCPAY)
      ICPATA = IPSMOF(LCPATA)
      IAI2T1 = IPSMOF(LAI2T1)
C
      CALL PSZRM(NORBS,NORBS,DUMP(ICPAYA),NORBS)
      HAVERD = .FALSE.
C    Method-specific initialization
      IF(HALFEL) THEN
          NACT = NUMOPN
          JMAX = NUMOPN
      ENDIF
      IF(DOCI) THEN
          NACT = NCIO
          JMIN = 1
          JMAX = NORBS
      ENDIF
C
      DO 1000 I=1,NACT
          IF(HALFEL) THEN
              JMIN  = I + 1
              IREAL = IOPEN1 + (I-1)
          ENDIF
          IF(DOCI) THEN
              IREAL = INDCIO(I)
          ENDIF
C        Transformation have to be performed only then
C        there is an active redundant varible with such
C        first index.
          HAVENZ = .FALSE.
          IF(HALFEL) THEN
              DO 100 J=JMIN,JMAX
                  IF( IACTID(I,J).GT.0 ) THEN
                      IF( QRED(IACTID(I,J)).NE.ZERO ) HAVENZ = .TRUE.
                  ENDIF
  100         CONTINUE
          ENDIF
          IF(DOCI) THEN
              DO 120 J=JMIN,JMAX
                  IF( QRED(IPACK(I,J)).NE.ZERO ) HAVENZ = .TRUE.
  120         CONTINUE
          ENDIF
          IF( .NOT.HAVENZ ) GOTO 1000
C
          DO 900 J=JMIN,JMAX
              IF(HALFEL) THEN
                  IJ    = IACTID(I,J)
                  JREAL = IOPEN1 + (J-1)
                  IF( IJ.LE.0 ) GOTO 900
              ENDIF
              IF(DOCI) THEN
                  IJ    = IPACK(I,J)
                  JREAL = J
              ENDIF
              SCALE = QRED(IJ)
              IF( SCALE.EQ.ZERO ) GOTO 900
              HAVERD = .TRUE.
              CALL DSYR2('U',NORBS,SCALE,CALP(1,IREAL),1,CALP(1,JREAL),
     .                   1,DUMP(ICPAYA),NORBS)
  900     CONTINUE
 1000 CONTINUE
C
      IF(HAVERD) THEN
          CALL PSRESP(DUMP(ICPAYA),DUMP(ICPAYA),DUMP(ICPAY),
     .                DUMP(IAI2T1),DUMP(ICPATA),DUMP(ICPATA),NORBS,
     .                .FALSE.)
          CALL PSXU2M(CALP,CALP,LDC,DUMP(ICPATA),DUMP(ICPATA),NORBS,
     .                QADD,DUMP)
      ELSE
          CALL PSZRV(IQSZ,QADD)
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSHMKZ(C,LDC,DUMP)
C
C   Compute AO representation of the Z-vector for the first
C   half-electron and CI derivatives.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      C      - Alpha orbital coefficients.
C      LDC    - Leading dimension of the C matrix.
C      DUMP   - Scratch array.
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options
C      PSDGBL - Global computation parameters
C      PSDGB2 - Response global computation parameters
C      PSDYNM - Dynamic memory handles
C      PSHALF - Half-electron parameters
C      PSPRT  - Printing unit
C      PSPRTF - Debugging output control flags
C      PSDCIP - CI parameters.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      MO representation of the half-lectron Z-vector consists
C      of non-redundant and redundant contributions, which can
C      be converted into AO representation according to the
C      expressions:
C
C      Z      =   Sum  Z    ( C     C     + C     C    ), mu != nu
C       mu nu     i,j   i j    mu i  nu j    nu i  mu j
C
C      Z      =   Sum  Z    C     C     , mu = nu
C       mu mu           i j  mu i  mu j
C
C      For the non-redundant part of the Z-vector, these expressions
C      are equivalent to the standard MO->AO transformation, with the
C      exception of diagonal elements, which are counted twice. 
C
C   Bugs:
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXOPN=10)
      PARAMETER (MAXACT=(MAXOPN*(MAXOPN+1))/2)
      PARAMETER (MAXCIO=LMACT)
C
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (IZERO=0)
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
      COMMON
     ./PSHALF/ HLFH(MAXOPN,MAXOPN), HLFG(MAXOPN,MAXOPN),
     .         NUMOIJ, IOI(MAXACT), IOJ(MAXACT),
     .         IOPEN1, NUMOPN, NMRACT, NMNACT, IACTID(MAXOPN,MAXOPN),
     .         IACTI(MAXACT), IACTJ(MAXACT)
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
      SAVE /PSDOPT/, /PSDGBL/, /PSDGB2/, /PSHALF/, /PSPRTF/, /PSDCIP/
      SAVE /PSPRT /
      DIMENSION C(LDC,*), DUMP(*), DUMMY(1)
      IPACK(I,J) = (I-1)*NORBS + J
C    Dynamic memory handles...
C    Since RHS for the Z-vector equations is no longer needed after
C    this point, we can safely use it to store non-redundant part of
C    the solution vector.
      IHLQ   = IPSMOF(LHLQ)
      IHLZA  = IPSMOF(LHLZA)
      IHLZM  = IPSMOF(LHLZM)
C    Method-specific parameters
      IF(HALFEL) NACT = NUMOPN
      IF(DOCI)   NACT = NCIO
C    Merge non-redundant part of the Z-vector in the Z-vector RHS
      IF( IQSWAP.GE.2 ) THEN
          CALL PSDRD(IURHS,1-ICPV1,DUMP(IHLQ),IQSZ)
      ELSE
          IQ = IPSMOF(LQ)
          CALL DCOPY(IQSZ,DUMP(IQ+(0-ICPV1)*IQSZ),1,DUMP(IHLQ),1)
      ENDIF
C
      IF(UHF) THEN
          WRITE(NB6,10100)
          STOP 'PSHMKZ'
      ENDIF
      IF( LPHALF.OR.LPCI ) THEN
          WRITE(NB6,11100) 
          CALL PSXPRT(DUMP(IHLQ))
          WRITE(NB6,11200)
          IF(HALFEL) WRITE(NB6,11202) (IACTI(I), IACTJ(I), 
     .                              DUMP(IHLQ+IQSZ+(I-1)), I=1,NMRACT)
          IF(DOCI) CALL PSDPGM(NORBS,NACT,DUMP(IHLQ+IQSZ),NORBS)
      ENDIF
C    Check magnitude of the coupling coefficients
      IF( IPRINT.GE.-1 ) THEN
          CALL PSXCHB(IZERO,DUMP(IHLQ))
      ENDIF
C    Transform Z-vector into AO basis. First, transform 
C    non-redundant part....
      CALL PSXM2A(C,C,LDC,DUMP(IHLQ),DUMP(IHLZM),DUMMY,NORBS,DUMP)
      IF( LPHALF.OR.LPCI ) THEN
          WRITE(NB6,11300)
          CALL PSDPSU(NORBS,DUMP(IHLZM),NORBS)
      ENDIF
C    ... then redundant part. We could have fetched vectors into
C    temporary storage and used DSYR2K, but since the number of
C    redundant variables should be small, it probably does not
C    matter.
      IF(DOCI) THEN
         JMAX = NORBS
      ENDIF
      DO 600 I=1,NACT
          IF(HALFEL) THEN
              IREAL = IOPEN1 + (I-1)
              JMAX  = I-1
          ENDIF
          IF(DOCI) THEN
              IREAL = INDCIO(I)
          ENDIF
          DO 500 J=1,JMAX
              IF(HALFEL) THEN
                  IJ    = IACTID(J,I)
                  JREAL = IOPEN1 + (J-1)
              ENDIF
              IF(DOCI) THEN
                  IJ    = IPACK(I,J)
                  JREAL = J
              ENDIF
              SCALE = DUMP(IHLQ+IQSZ+(IJ-1))
              IF( SCALE.EQ.ZERO ) GOTO 500
              CALL DSYR2('U',NORBS,SCALE,C(1,IREAL),1,C(1,JREAL),1,
     .                   DUMP(IHLZM),NORBS)
  500     CONTINUE
  600 CONTINUE
C    Diagonal elements were counted twice, correct it.
      CALL DSCAL(NORBS,HALF,DUMP(IHLZM),NORBS+1)
C    Pack computed representation of Z-vector
      CALL PSU2PK(NORBS,DUMP(IHLZM),NORBS,DUMP(IHLZA))
      IF( LPHALF.OR.LPCI ) THEN
          WRITE(NB6,11600)
          CALL PSDPPM(NORBS,DUMP(IHLZA))
      ENDIF
C
      RETURN
10100 FORMAT(' UHF IS NOT SUPPORTED IN HALF-ELECTRON MODE ')
C
11100 FORMAT(' M.O. REPRESENTATION OF THE HALF-ELECTRON Z-VECTOR:'/)
11200 FORMAT(' REDUNDANT PART:')
11202 FORMAT((2(2X,I4),5X,G15.9))
11300 FORMAT(' NON-REDUNDANT A.O. PART: '/)
11600 FORMAT(' A.O. REPRESENTATION OF THE HALF-ELECTRON Z-VECTOR:'/)
      END
