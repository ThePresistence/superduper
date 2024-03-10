C     ******************************************************************
C
C     Code related to the small-CI derivatives.
C
C     ******************************************************************
C
C     This code is *NOT* seriously tuned for optimal performance,
C     since it is expected to have a negligible effect on overall
C     execution time.
C
      SUBROUTINE PSCIN1
C
C   First part of the CI derivatives initialization.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      None.
C
C   Accessed common blocks:
C
C      PSDGBL,
C      PSDGB2 - Various parameters.
C      PSPRT  - Printing unit.
C      PSPRTF - Printing flags.
C      PSDOPT - Control options.
C      PSORD1 - MO ordering array.
C      CCIS   - CI parameters from the main program (default cases).
C      HALFE  - Multiplicity.
C      OCCNM  - Occupation numbers. Note that these are expected
C               to be sorted by PSHINI, while MOs in HALFE and CCIS
C               are in the original numbering system.
C
C   Modified common blocks:
C
C      PSDCIP - CI parameters.
C         NCIO   - Number of MOs involved in CI.
C         NCIC   - Number of configurations involved in CI.
C         ICIMOD - CI version. Small numbers are the standard 
C                  cases from the main code.
C            1 = 2x2 closed-shell CI
C            2 = 3x3 closed-shell CI
C            3 = 3x3 open-shell singlet CI
C            4 = 2x2 doublet (double->single) CI
C            5 = 2x2 doublet (single->vacant) CI
C         NCIGAM - Size of the two-particle density matrix
C                  over CI-active MOs.
C         NCIFS  - Number of F coefficients.
C         NCIES  - Number of E coefficients.
C         INDCIO - Array of MAXCIO elements, containing indices
C                  of the CI-active MOs in orbital coefficients
C                  matrix.
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMX, LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
C-AK  PARAMETER (IHUGE=2147483647)
C-AK  PARAMETER (IHUGE=9223372036854775807)
      PARAMETER (IHUGE=HUGE(1))
C
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ONE =1.0D0)
      PARAMETER (EPS =1.0D-10)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL ISONE, ISHALF, ISZERO
      EXTERNAL IPSSRI
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
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
     ./PSORD1/ IORDEA(LMX), IORDEB(LMX)
     ./PSPRT / NB6
C
      COMMON
     ./CCIS  / K,L,NC,LROOT
     ./CIMOS / IMOCI(LMX)
     ./INOPT2/ IN2(300)
     ./OCCNM / OCCA(LMX),OCCB(LMX)
      SAVE /PSDGBL/, /PSDOPT/, /PSDCIP/, /PSORD1/, /PSPRT /
      CHARACTER*20 MNAMES(6)
      DATA MNAMES/' 2x2 CLOSED SHELL ',' 3x3 CLOSED SHELL ',
     .            ' 3x3 OPEN SINGLET ',' 2x2 DOUBLET (D->S) ',
     .            ' 2x2 DOUBLET (S->V) ',' GUGA-CI '/
C
      ISONE (X) = ABS(X-ONE ).LT.EPS
      ISHALF(X) = ABS(X-HALF).LT.EPS
      ISZERO(X) = ABS(X-ZERO).LT.EPS
C
      IF( UHF.OR.HALFEL ) THEN
          WRITE(NB6,10000)
          STOP 'PSINC1'
      ENDIF
C
C    Check for type of CI treatment.
C    iabsci = 1 : Minimal CI.
C    iabsci = 5 : General GUGA-CI.
C
      IABSCI = ABS(IN2(77))
C
C    Decide which of the standard cases we are in for.
C    This is a replica of CIS() from the main program,
C    with the addition of reasonably savage (he-he ;-) 
C    consistency testing.
C
      IF(IABSCI.EQ.1) THEN
         IMULT = IN2(66)
         KTRUE = K
         LTRUE = L
         IF(DORESP) THEN
             KTRUE = IPSSRI(KTRUE,IORDEA)
             LTRUE = IPSSRI(LTRUE,IORDEA)
         ENDIF
         OCCK = OCCA(KTRUE)
         OCCL = OCCA(LTRUE)
         IF( IMULT.EQ.0 ) THEN
C           Closed-shell CI. K should be singly occupied, L should be vacant.
             IF( .NOT.( ISONE(OCCK).AND.ISZERO(OCCL) ) ) THEN
                 WRITE(NB6,10010) OCCK, OCCL
                 STOP 'PSCIN1'
             ENDIF
             IF( NC.EQ.2 ) THEN
                 ICIMOD = 1
                 NCIO   = 2
                 NCIC   = 2
                 NCIFS  = 5
                 NCIES  = 2
             ELSEIF( NC.EQ.3 ) THEN
                 ICIMOD = 2
                 NCIO   = 2
                 NCIC   = 3
                 NCIFS  = 9
                 NCIES  = 4
             ELSE
                 WRITE(NB6,10020) NC
                 STOP 'PSCIN1'
             ENDIF
         ELSE IF( IMULT.EQ.1 ) THEN
C           Open-shell singlet CI. Both K and L should be half-occupied.
             IF( .NOT.( ISHALF(OCCK).AND.ISHALF(OCCL) ) ) THEN
                 WRITE(NB6,10030) OCCK, OCCL
                 STOP 'PSCIN1'
             ENDIF
             ICIMOD = 3
             NCIO   = 2
             NCIC   = 3
             NCIFS  = 16
             NCIES  = 4
         ELSE IF( IMULT.EQ.2 ) THEN
C           Doublet CI. 
             IF( ISONE(OCCK).AND.ISHALF(OCCL) ) THEN
                 ICIMOD = 4
                 NCIO   = 2
                 NCIC   = 2
                 NCIFS  = 5
                 NCIES  = 2
             ELSE IF( ISHALF(OCCK).AND.ISZERO(OCCL) ) THEN
                 ICIMOD = 5
                 NCIO   = 2
                 NCIC   = 2
                 NCIFS  = 5
                 NCIES  = 2
             ELSE
                 WRITE(NB6,10040) OCCK, OCCL
                 STOP 'PSCIN1'
             ENDIF
         ELSE
             WRITE(NB6,10050) IMULT
             STOP 'PSCIN1'
         ENDIF
         INDCIO(1) = KTRUE
         INDCIO(2) = LTRUE
C
C    Initialization for general GUGA-CI.
C    The relevant variables are obtained from input in INGUGA (IN2)
C    and from the initial call to GUGACI in DYNGUG.
C
      ELSE IF(IABSCI.EQ.5) THEN
         ICIMOD = 6
         NCIO   = IN2(131)+IN2(132)
         NCIC   = 1
         NCIGAM = (NCIO*(NCIO+1)*(NCIO**2+NCIO+2))/8
         NCIFS  = 1
         NCIES  = 1
         DO 10 I=1,NCIO
         INDCIO(I) = IMOCI(I)
   10    CONTINUE
      ELSE
         WRITE(NB6,10060) IABSCI
         STOP 'PSCIN1'
      ENDIF
C
C    From this point, code is Ok for general case as well
C
      DNCIO  = DBLE(NCIO)
      TEMP   = DNCIO*(DNCIO+1)*(DNCIO**2+DNCIO+2)
      IF( TEMP.GE.DBLE(IHUGE) ) THEN
          WRITE(NB6,10200) TEMP/8
          STOP 'PSCIN1'
      ENDIF
      NCIGAM = (NCIO*(NCIO+1)*(NCIO**2+NCIO+2))/8
C     TEMP   = DBLE(NCIGAM)*DBLE(NCIC)*DBLE(NCIC+1)
C     IF( TEMP.GE.DBLE(IHUGE) ) THEN
C         WRITE(NB6,10210) TEMP/2
C         STOP 'PSCIN1'
C     ENDIF
C
      IF( IPRINT.GE.1 ) THEN
          WRITE(NB6,11500) NCIO, NCIC, ICIMOD, MNAMES(ICIMOD), 
     .                   NCIGAM, NCIFS, NCIES
          IF( IPRINT.GE.5 ) THEN
              WRITE(NB6,11510) ( INDCIO(I), I=1,NCIO )
          ENDIF
      ENDIF
      RETURN
10000 FORMAT(' PSINC1 CALLED WITH UHF OR HALFEL EQUAL TO .TRUE.')
10010 FORMAT(' WRONG OCCUPATION NUMBERS IN MINIMAL CLOSED SHELL CI ',
     .       'IN PSINC1: ',F10.6,' AND ',F10.6)
10020 FORMAT(' WRONG NUMBER OF CONFIGURATIONS IN MINIMAL CLOSED ',
     .       'SHELL CI IN PSINC1: ',I5)
10030 FORMAT(' WRONG OCCUPATION NUMBERS IN MINIMAL OPEN SHELL SINGLET',
     .       ' CI IN PSINC1: ',F10.6,' AND ',F10.6)
10040 FORMAT(' WRONG OCCUPATION NUMBERS IN MINIMAL DOUBLET CI IN ',
     .       'PSINC1: ',F10.6,' AND ',F10.6)
10050 FORMAT(' WRONG MULTIPLICITY IN MINIMAL CI IN PSINC1: ',I5)
10060 FORMAT(' ILLEGAL INPUT VALUE FOR IABSCI IN PSINC1: ',I5)
10200 FORMAT(' INTEGER CANNOT HOLD THE SIZE OF THE TWO-PARTICLE ',
     .       'DENSITY (',F15.0,') IN PSINC1.')
10210 FORMAT(' INTEGER CANNOT HOLD LARGEST PACKED F COEFFICIENT INDEX ('
     .       ,G15.9,') IN PSINC1.')
11500 FORMAT(' CI DERIVATIVES ARE COMPUTED WITH: '
     .      /' NCIO   = ',I8,' NCIC  = ',I8,' ICIMOD = ',I8,' (',A,') '
     .      /' NCIGAM = ',I8,' NCIFS = ',I8,' NCIES  = ',I8)
11510 FORMAT(' CI-ACTIVE MOLECULAR ORBITALS ARE:'/(10(1X,I7)))
      END
C
      SUBROUTINE PSCIN2(CALP,LDC,EALP,DUMP,ONEDEN,TWODEN)
C
C   Second (and last) part of the CI derivatives initialization.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      CALP   - Alpha orbital coefficients.
C      LDC    - Leading dimension of CALP.
C      EALP   - Alpha orbital energies.
C      DUMP   - Memory arena.
C      ONEDEN - Precalculated one-particle density matrix.
C      TWODEN - Precalculated two-particle density matrix.
C
C      Note: ONEDEN,TWODEN used only for GUGA-CI.
C
C   Accessed common blocks:
C
C      PSDGBL - Number of orbitals is needed from here.
C      PSDYNM - Dynamic memory IDs.
C      PSDCIP - Various CI parameters.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C      We are passing (correctly-sized) DOUBLE PRECISION chunks as
C      INTEGER arrays.
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
      PARAMETER (EPS=1.D-12)
      PARAMETER (ONE=1.D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      EXTERNAL IPSMOF
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
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
      SAVE /PSDGBL/, /PSDCIP/, /PSPRTF/, /PSPRT /
      DIMENSION CALP(LDC,*), EALP(*), DUMP(*)
      DIMENSION ONEDEN(*), TWODEN(*)
C    Get dynamic memory locations
C     UNCOMMENT THE FOLLOWING TWO LINES TO OBTAIN
C     DEBUG PRINT OF CI-RELATED DENSITIES (AK).
C     CALL PRDEN1(ONEDEN)
C     CALL PRDEN2(TWODEN)
      ICIFI  = IPSMOF(LCIFI)
      ICIFC  = IPSMOF(LCIFC)
      ICIEI  = IPSMOF(LCIEI)
      ICIEC  = IPSMOF(LCIEC)
      ICIVEC = IPSMOF(LCIVEC)
      ICIORB = IPSMOF(LCIORB)
      ICIGAM = IPSMOF(LCIGAM)
      ICIDLT = IPSMOF(LCIDLT)
C    First step: construct F and E coefficients for CI, and fetch 
C    configurations vector.
      IF( ICIMOD.GE.1 .AND. ICIMOD.LE.5 ) THEN
          CALL PSCIDF(DUMP(ICIFI),DUMP(ICIFC),DUMP(ICIEI),
     .                DUMP(ICIEC),DUMP(ICIVEC))
      ELSE IF( ICIMOD.EQ.6 ) THEN
C         MODGUG = 2
C         CALL PSGUGA(DUMP(ICIFI),DUMP(ICIFC),DUMP(ICIEI),
C    .                DUMP(ICIEC),DUMP(ICIVEC),CIVEC,EGEN,FGEN,
C    .                XEGEN,XFGEN,LMEPS,LMFPS,MODGUG)
C         IF( LPCI ) THEN
C            CALL PSGUGP(DUMP(ICIFI),DUMP(ICIFC),DUMP(ICIEI),
C    .                   DUMP(ICIEC),DUMP(ICIVEC))
C         ENDIF
      ELSE
          WRITE(NB6,10100) ICIMOD
          STOP 'PSCIN2'
      ENDIF
C    From this point, we don't care about the specific CI case anymore.
C    Check norm of the CI state vector
C     SNORM = DNRM2(NCIC,DUMP(ICIVEC),1)
C     IF( ABS(SNORM-ONE).GE.EPS ) THEN
C         WRITE(NB6,10200) SNORM
C         STOP 'PSCIN2'
C     ENDIF
C    Fetch orbital coefficients for CI-active space into a private array, 
C    so what we can treat them as contiguous in all cases.
      DO 1000 ITO=1,NCIO
          IFROM = INDCIO(ITO)
          CALL DCOPY(NORBS,CALP(1,IFROM),1,DUMP(ICIORB+(ITO-1)*NORBS),1)
 1000 CONTINUE
      IF( ICIMOD.NE.6 ) THEN
C        Compute two-particle density matrix over CI-active MOs
          CALL PSCIGA(DUMP(ICIFI),DUMP(ICIFC),DUMP(ICIVEC),DUMP(ICIGAM))
C        ... and Lagrangian multipliers for the same MOs.
          CALL PSCILG(DUMP(ICIEI),DUMP(ICIEC),DUMP(ICIVEC),DUMP(ICIDLT))
      ELSE
C        Fetch these matrices from GUGA-CI calculation.
          CALL PSGUG1(ONEDEN,DUMP(ICIDLT))
          CALL PSGUG2(TWODEN,DUMP(ICIGAM))
      ENDIF
C    This is about it.
      IF(LPCI) THEN
          WRITE(NB6,11010)
          CALL PSDPPM((NCIO*(NCIO+1))/2,DUMP(ICIGAM))
          WRITE(NB6,11020)
          CALL PSDPDR(NCIO,DUMP(ICIDLT))
          WRITE(NB6,11030)
          CALL PSDPGM(NORBS,NCIO,DUMP(ICIORB),NORBS)
      ENDIF
      RETURN
10100 FORMAT(' ICIMOD = ',I5,' IS NOT IMPLEMENTED IN PSCIN2.')
10200 FORMAT(' NORM OF THE CI STATE VECTOR (',G20.14,') IS NOT UNIT.')
11010 FORMAT(/' TWO-PARTICLE DENSITY MATRIX: '/)
11020 FORMAT(/' LAGRANGIAN MULTIPLIERS: '/)
11030 FORMAT(/' CI-ACTIVE MOS: '/)
      END
C
      SUBROUTINE PSCIDF(IF,F,IE,E,VEC)
C
C   Sets up coefficients and configurations vector for the standard
C   CI cases.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IF     - (Output) Packed indices of the F coefficients.
C      F      - (Output) F coefficients.
C      IE     - (Output) Packed indices of the E coefficients.
C      E      - (Output) E coefficients.
C      VEC    - (Output) Configurations vector.
C
C   Accessed common blocks:
C
C      PSDCIP - Various CI parameters.
C      CCIS   - CI state we are interested in.
C      XCI    - CI configurations vectors.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      We write matrix elements of the CI Hamiltomian as:
C
C            active  lm             active  lm
C      H   =   Sum  f     (pq|rs) +   Sum  e   eps
C       lm    pqrs   pqrs              r    r     r
C
C      where (pq|rs) are ERIs over molecular orbitals, and
C      eps are the orbital energies of the MO in the CI-active
C      space.
C
C      f and e coefficient matrices are sparse, and are
C      represented by a (packed) index plus a coefficient
C      value.
C
C   Bugs:
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL IPSCIF
      PARAMETER (MAXCIO=LMACT)
      PARAMETER (EIGHTH=0.125D0)
      PARAMETER (RSQ32 =0.17677669529663688110021109052621225982121D0)
      PARAMETER (FOURTH=0.250D0)
      PARAMETER (RSQ8  =0.35355339059327376220042218105242451964242D0)
      PARAMETER (C3O8  =0.375D0)
      PARAMETER (HALF  =0.500D0)
      PARAMETER (ONE   =1.000D0)
      PARAMETER (TWO   =2.000D0)
C
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      COMMON
     ./CCIS  / K,L,NC,LROOT
     ./XCI   / EO(6),ENGYCI(3),VECTCI(9),ECI(3),BKCAL(3)
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSDCIP/, /PSPRTF/, /PSPRT /
      DIMENSION IF(NCIFS), F(NCIFS), IE(NCIES), E(NCIES), VEC(NCIC)
C    First, fetch configurations vector. 
C    Thiel's code declares VECTCI as 3 by 3 matrix, but actually
C    handles it as either 3 by 3 or 2 by 2 matrix, depending on NCIC.
      DO 100 I=1,NCIC
          VEC(I) = VECTCI((LROOT-1)*NCIC+I)
  100 CONTINUE
      IF(LPCI) THEN
          WRITE(NB6,11050) LROOT, (VEC(I), I=1,NCIC)
      ENDIF
C    Now, produce F and E coefficients.
      GOTO (1001,2001,3001,4001,5001), ICIMOD
          WRITE(NB6,10100) ICIMOD
          STOP 'PSCIDF'
 1001 CONTINUE
C        2x2 closed shell CI
          IF( 1) =  IPSCIF(1,2,1,2,1,2)
          F ( 1) =  FOURTH
          IF( 2) =  IPSCIF(2,2,1,1,2,2)
          F ( 2) = -TWO
          IF( 3) =  IPSCIF(2,2,1,2,1,2)
          F ( 3) =  HALF
          IF( 4) =  IPSCIF(2,2,1,1,1,1)
          F ( 4) =  ONE
          IF( 5) =  IPSCIF(2,2,2,2,2,2)
          F ( 5) =  ONE
          IE( 1) =  IPSCIF(2,2,1,1,1,1)
          E ( 1) = -TWO
          IE( 2) =  IPSCIF(2,2,2,1,1,1)
          E ( 2) =  TWO
      GOTO 9001
 2001 CONTINUE
C        3x3 closed shell CI, continues above
          IF( 6) =  IPSCIF(2,3,1,1,1,2)
          F ( 6) =  RSQ8
          IF( 7) =  IPSCIF(2,3,1,2,2,2)
          F ( 7) = -RSQ8
          IF( 8) =  IPSCIF(3,3,1,1,2,2)
          F ( 8) = -HALF
          IF( 9) =  IPSCIF(3,3,1,2,1,2)
          F ( 9) =  HALF
          IE( 3) =  IPSCIF(3,3,1,1,1,1)
          E ( 3) = -ONE
          IE( 4) =  IPSCIF(3,3,2,1,1,1)
          E ( 4) =  ONE
      GOTO 1001
 3001 CONTINUE
C        3x3 open-shell CI
          IF( 1) =  IPSCIF(1,1,1,1,1,1)
          F ( 1) = -FOURTH
          IF( 2) =  IPSCIF(1,1,2,2,2,2)
          F ( 2) = -FOURTH
          IF( 3) =  IPSCIF(1,1,1,2,1,2)
          F ( 3) =  C3O8
          IF( 4) =  IPSCIF(1,2,1,1,1,2)
          F ( 4) =  RSQ32
          IF( 5) =  IPSCIF(1,2,1,2,2,2)
          F ( 5) = -RSQ32
          IF( 6) =  IPSCIF(2,2,1,1,1,1)
          F ( 6) =  FOURTH
          IF( 7) =  IPSCIF(2,2,2,2,2,2)
          F ( 7) =  FOURTH
          IF( 8) =  IPSCIF(2,2,1,1,2,2)
          F ( 8) = -HALF
          IF( 9) =  IPSCIF(2,2,1,2,1,2)
          F ( 9) =  EIGHTH
          IF(10) =  IPSCIF(1,3,1,1,1,2)
          F (10) = -RSQ32
          IF(11) =  IPSCIF(1,3,1,2,2,2)
          F (11) = +RSQ32
          IF(12) =  IPSCIF(2,3,1,2,1,2)
          F (12) =  FOURTH
          IF(13) =  IPSCIF(3,3,1,1,1,1)
          F (13) =  FOURTH
          IF(14) =  IPSCIF(3,3,2,2,2,2)
          F (14) =  FOURTH
          IF(15) =  IPSCIF(3,3,1,1,2,2)
          F (15) = -HALF
          IF(16) =  IPSCIF(3,3,1,2,1,2)
          F (16) =  EIGHTH
          IE( 1) =  IPSCIF(2,2,1,1,1,1)
          E ( 1) =  ONE
          IE( 2) =  IPSCIF(2,2,2,1,1,1)
          E ( 2) = -ONE
          IE( 3) =  IPSCIF(3,3,1,1,1,1)
          E ( 3) = -ONE
          IE( 4) =  IPSCIF(3,3,2,1,1,1)
          E ( 4) =  ONE
      GOTO 9001
 4001 CONTINUE
C        2x2 doublet CI, occ->open
          IF( 1) =  IPSCIF(1,1,2,2,2,2)
          F ( 1) = -FOURTH
          IF( 2) =  IPSCIF(1,2,1,2,2,2)
          F ( 2) =  EIGHTH
          IF( 3) =  IPSCIF(2,2,2,2,2,2)
          F ( 3) =  FOURTH
          IF( 4) =  IPSCIF(2,2,1,1,2,2)
          F ( 4) = -HALF
          IF( 5) =  IPSCIF(2,2,1,2,1,2)
          F ( 5) =  EIGHTH
          IE( 1) =  IPSCIF(2,2,1,1,1,1)
          E ( 1) = -ONE
          IE( 2) =  IPSCIF(2,2,2,1,1,1)
          E ( 2) =  ONE
      GOTO 9001
 5001 CONTINUE
C        2x2 doublet CI, open->vac
          IF( 1) =  IPSCIF(1,1,1,1,1,1)
          F ( 1) = -FOURTH
          IF( 2) =  IPSCIF(1,2,1,1,1,2)
          F ( 2) = -EIGHTH
          IF( 3) =  IPSCIF(2,2,1,1,1,1)
          F ( 3) =  FOURTH
          IF( 4) =  IPSCIF(2,2,1,1,2,2)
          F ( 4) = -HALF
          IF( 5) =  IPSCIF(2,2,1,2,1,2)
          F ( 5) =  EIGHTH
          IE( 1) =  IPSCIF(2,2,1,1,1,1)
          E ( 1) = -ONE
          IE( 2) =  IPSCIF(2,2,2,1,1,1)
          E ( 2) =  ONE
      GOTO 9001
C
 9001 CONTINUE
      RETURN
10100 FORMAT(' ICIMOD = ',I5,', NOT UNDERSTOOD IN PSCIDF.')
11050 FORMAT(' CI STATE VECTOR FOR THE STATE ',I4,' IS:'/(5(1X,F16.9)))
      END
C
      SUBROUTINE PSGUGA (IFCOF,FCOF,IECOF,ECOF,VEC,CIVEC,
     1                   EGEN,FGEN,IEGEN,IFGEN,LMEPS,LMFPS,MODGUG)
C     *
C     CONVERSION FROM GUGA-CI GENERATORS TO CI-COEFFICIENTS.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     IFCOF(LMFPS)    PACKED INDICES FOR TWO-ELECTRON COEFFICIENTS (O).
C     FCOF(LMFPS)     TWO-ELECTRON COEFFICIENTS (O).
C     IECOF(LMEPS)    PACKED INDICES FOR ONE-ELECTRON COEFFICIENTS (O).
C     ECOF(LMEPS)     ONE-ELECTRON COEFFICIENTS (O).
C     VEC(LMVEC)      CI VECTOR FOR STATE OF INTEREST (O).
C     CIVEC(LMVEC)    CI VECTOR FOR STATE OF INTEREST (I).
C     EGEN(LMEGEN)    ONE-ELECTRON GENERATORS (I).
C     FGEN(LMFGEN)    TWO-ELECTRON GENERATORS (I).
C     IEGEN(LMEGEN)   PACKED INDICES FOR ONE-ELECTRON GENERATORS (I).
C     IFGEN(LMFGEN)   PACKED INDICES FOR TWO-ELECTRON GENERATORS (I).
C     LMEPS           DIMENSION OF ARRAYS IECOF AND ECOF (O).
C     LMFPS           DIMENSION OF ARRAYS IECOF AND ECOF (O).
C     MODGUG          MODE OF OPERATION (I).
C                     =1 DETERMINE ARRAY SIZES FOR ONE-ELECTRON (LMEPS)
C                        AND TWO-ELECTRON (LMFPS) COEFFICIENTS.
C                     =2 COMPUTE THESE CI COEFFICIENTS (ECOF,FCOF) AND
C                        THE ASSOCIATED PACKED INDICES (IECOF,IFCOF).
C     *
C
C     Module logic:
C
C      We write matrix elements of the CI Hamiltomian as:
C
C            active  lm             active  lm
C      H   =   Sum  f     (pq|rs) +   Sum  e   eps
C       lm    pqrs   pqrs              r    r     r
C
C      where (pq|rs) are ERIs over molecular orbitals, and
C      eps are the orbital energies of the MO in the CI-active
C      space.
C
C      f and e coefficient matrices are sparse, and are
C      represented by a (packed) index plus a coefficient
C      value.
C
      USE LIMIT, ONLY: LM1, LMX, LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
      PARAMETER (MAXTMP=(MAXCIO*(MAXCIO+1))/2)
      PARAMETER (C3O4=0.75D0)
      PARAMETER (EIGHT=8.0D0)
      PARAMETER (TINY=1.0D-08)
      COMMON
     ./CIMOS / IMOCI(LMX)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./GUGA1 / LMEGEN,LMFGEN,LMVEC
     ./INOPT2/ IN2(300)
     ./HALFE / IODD(2)
      DIMENSION EGEN(LMEGEN),FGEN(LMFGEN)
      DIMENSION IEGEN(LMEGEN),IFGEN(LMFGEN)
      DIMENSION ECOF(*),FCOF(*) 
      DIMENSION IECOF(*),IFCOF(*) 
      DIMENSION VEC(LMVEC),CIVEC(LMVEC)
      DIMENSION E(MAXTMP),ETMP(MAXTMP),INDI(MAXTMP),INDJ(MAXTMP)
      DIMENSION FTMP(MAXTMP),ITMP(MAXTMP)
      DIMENSION OCC(MAXCIO)
C *** INITIALIZATION.
C     JPRINT = IN2(42)
      ICI1   = IN2(131)
      ICI2   = IN2(132)
      NCIO   = ICI1+ICI2
      NCIONE = (NCIO*(NCIO+1))/2
      NCIGAM = (NCIONE*(NCIONE+1))/2
      NCILOW = (ICI1*(ICI1+1))/2
C
C *** ORBITAL OCCUPATION NUMBERS
      DO 5 I=1,ICI1
      OCC(I) = TWO
 5    CONTINUE
      DO 6 I=1,NCIO
      IF(IMOCI(I).EQ.IODD(1) .OR. IMOCI(I).EQ.IODD(2))
     .   OCC(I) = ONE
 6    CONTINUE
C
C *** NONZERO ONE-ELECTRON COEFFICIENTS.
      IEPS = 0
C     STORE DIAGONAL ONE-ELECTRON GENERATORS FOR GIVEN MU-NU.
      DO 10 K=1,NCIO
      E(K) = ZERO
   10 CONTINUE
      MUNU = (IEGEN(1)-1) / NCIGAM
      DO 20 NE=1,LMEGEN
      MUNUE = (IEGEN(NE)-1) / NCIGAM
      IF(MUNUE.NE.MUNU) THEN
         DO 40 K=1,NCIO
         ETEMP = E(K)
C-AK     IF(K.LE.ICI1 .AND. MU.EQ.NU) ETEMP = ETEMP-TWO
         IF(K.LE.ICI1 .AND. MU.EQ.NU) ETEMP = ETEMP-OCC(K)
         IF(ABS(ETEMP).GT.TINY) THEN
            IEPS = IEPS+1
            IF(MODGUG.EQ.2) THEN
               IECOF(IEPS) = IPSCIF (MU,NU,1,1,K,K)
               ECOF (IEPS) = ETEMP
            ENDIF
         ENDIF
         E(K) = ZERO
   40    CONTINUE
         IF (MUNUE.EQ.0) GO TO 60
         MUNU = MUNUE
      ENDIF
      CALL PSCIFU (IEGEN(NE),MU,NU,I,J,K,L)
      IF(K.EQ.L) E(K) = EGEN(NE)
   20 CONTINUE
   50 CONTINUE
   60 CONTINUE
C     SAVE NUMBER OF NONZERO ONE-ELECTRON CI COEFFICIENTS.
      IF(MODGUG.EQ.1) THEN
         LMEPS = IEPS
      ENDIF
C
C *** NONZERO TWO-ELECTRON COEFFICIENTS.
      IFPS  = 0
      NFPS  = 0
      NOLDE = 1
      NOLDF = 1
C *** PROCESS TWO-ELECTRON GUGA GENERATORS FOR GIVEN MU-NU.
      MUNU = (IFGEN(1)-1) / NCIGAM
      DO 100 NF1=1,LMFGEN
      MUNUF = (IFGEN(NF1)-1) / NCIGAM
      IF(MUNUF.NE.MUNU) THEN
         NU   = NINT(SQRT(TWO*MUNU-C3O4))
         MU   = MUNU - (NU*(NU-1))/2
         IF(MU.EQ.NU) THEN
            DO 120 IJ=1,NCILOW
            E(IJ) = ZERO
  120       CONTINUE
            NUME = 0
            DO 130 NE=NOLDE,LMEGEN-1
            CALL PSCIFU (IEGEN(NE),MUE,NUE,I,J,K,L)
            IF(MUE.LT.MU .OR. (MUE.EQ.MU .AND. NUE.LT.NU)) GO TO 130
            IF(MUE.NE.MU .OR. NUE.NE.NU) GO TO 140
            KL = K + (L*(L-1))/2
            IF(KL.LE.NCILOW) THEN
               E(KL)      = EGEN(NE)
            ELSE
               NUME       = NUME + 1
               ETMP(NUME) = EGEN(NE)
               INDI(NUME) = K
               INDJ(NUME) = L
            ENDIF
  130       CONTINUE
  140       NOLDE = NE
C ***       INCLUDE ADDITIONAL CONTRIBUTIONS.
C-AK        LOWER ACTIVE ORBITALS FIRST.
            IJ = 0
            DO 205 I=1,ICI1
            DO 200 J=1,I
            IJ = IJ + 1
C           EXTRA TERMS INVOLVING ONE-ELECTRON GENERATORS.
            IF(ABS(E(IJ)).GT.TINY) THEN
               DO 170 K=1,ICI1
               ITEMP = IPSCIF (MU,NU,I,J,K,K)
               DO 150 NF=1,IFPS
               IF(ITEMP.EQ.ITMP(NF)) THEN
C-AK              FTMP(NF) = FTMP(NF) - E(IJ)*TWO
                  FTMP(NF) = FTMP(NF) - E(IJ)*OCC(K)
                  GO TO 155
               ENDIF
  150          CONTINUE
               IFPS = IFPS+1
C-AK           FTMP(IFPS) = -E(IJ)*TWO
               FTMP(IFPS) = -E(IJ)*OCC(K)
               ITMP(IFPS) = ITEMP
  155          CONTINUE
               ITEMP = IPSCIF (MU,NU,I,K,K,J)
               DO 160 NF=1,IFPS
               IF(ITEMP.EQ.ITMP(NF)) THEN
C-AK              FTMP(NF) = FTMP(NF) + E(IJ)
                  FTMP(NF) = FTMP(NF) + E(IJ)*OCC(K)*PT5
                  GO TO 165
               ENDIF
  160          CONTINUE
               IFPS = IFPS+1
C-AK           FTMP(IFPS) = E(IJ)
               FTMP(IFPS) = E(IJ)*OCC(K)*PT5
               ITMP(IFPS) = ITEMP
  165          CONTINUE
  170          CONTINUE
            ENDIF
C           EXTRA COULOMB AND EXCHANGE TERMS.
            ITEMP = IPSCIF(MU,NU,I,I,J,J)
            IF(I.EQ.J) THEN
C-AK           FTEMP = ONE
               FTEMP = OCC(I)*OCC(I)*PT25
            ELSE
C-AK           FTEMP = FOUR
               FTEMP = OCC(I)*OCC(J)
            ENDIF
            DO 180 NF=1,IFPS
            IF(ITEMP.EQ.ITMP(NF)) THEN
               FTMP(NF) = FTMP(NF) + FTEMP
               GO TO 185
            ENDIF
  180       CONTINUE
            IFPS = IFPS+1
            FTMP(IFPS) = FTEMP
            ITMP(IFPS) = ITEMP
  185       CONTINUE
            IF(I.NE.J) THEN
               ITEMP = IPSCIF(MU,NU,I,J,I,J)
               DO 190 NF=1,IFPS
               IF(ITEMP.EQ.ITMP(NF)) THEN
C-AK              FTMP(NF) = FTMP(NF) - TWO
                  FTMP(NF) = FTMP(NF) - PT5*OCC(I)*OCC(J)
                  GO TO 195
               ENDIF
  190          CONTINUE
               IFPS = IFPS+1
C-AK           FTMP(IFPS) = -TWO
               FTMP(IFPS) = -PT5*OCC(I)*OCC(J)
               ITMP(IFPS) = ITEMP
  195          CONTINUE
            ENDIF
  200       CONTINUE
  205       CONTINUE
C-AK        HIGHER ACTIVE ORBITALS.
            DO 2205 IJ=1,NUME
            EIJ = ETMP(IJ)
            I   = INDI(IJ)
            J   = INDJ(IJ)
C           EXTRA TERMS INVOLVING ONE-ELECTRON GENERATORS.
            DO 2170 K=1,ICI1
            ITEMP = IPSCIF (MU,NU,I,J,K,K)
            DO 2150 NF=1,IFPS
            IF(ITEMP.EQ.ITMP(NF)) THEN
C-AK           FTMP(NF) = FTMP(NF) - E(IJ)*TWO
               FTMP(NF) = FTMP(NF) - EIJ*OCC(K)
               GO TO 2155
            ENDIF
 2150       CONTINUE
            IFPS = IFPS+1
C-AK        FTMP(IFPS) = -E(IJ)*TWO
            FTMP(IFPS) = -EIJ*OCC(K)
            ITMP(IFPS) = ITEMP
 2155       CONTINUE
            ITEMP = IPSCIF (MU,NU,I,K,K,J)
            DO 2160 NF=1,IFPS
            IF(ITEMP.EQ.ITMP(NF)) THEN
C-AK           FTMP(NF) = FTMP(NF) + E(IJ)
               FTMP(NF) = FTMP(NF) + EIJ*OCC(K)*PT5
               GO TO 2165
            ENDIF
 2160       CONTINUE
            IFPS = IFPS+1
C-AK        FTMP(IFPS) = E(IJ)
            FTMP(IFPS) = EIJ*OCC(K)*PT5
            ITMP(IFPS) = ITEMP
 2165       CONTINUE
 2170       CONTINUE
C           NO EXTRA COULOMB AND EXCHANGE TERMS HERE.
 2205       CONTINUE
         ELSE
            NUME = 0
            DO 1130 NE=NOLDE,LMEGEN-1
            CALL PSCIFU (IEGEN(NE),MUE,NUE,I,J,K,L)
            IF(MUE.LT.MU .OR. (MUE.EQ.MU .AND. NUE.LT.NU)) GO TO 1130
            IF(MUE.NE.MU .OR. NUE.NE.NU) GO TO 1140
C-AK        KL     = K + (L*(L-1))/2
C-AK        E(KL)  = EGEN(NE)
            NUME       = NUME + 1
            E(NUME)    = EGEN(NE)
            INDI(NUME) = K
            INDJ(NUME) = L
 1130       CONTINUE
 1140       NOLDE = NE
C ***       INCLUDE ADDITIONAL CONTRIBUTIONS.
C           LOOP OVER ACTIVE ORBITALS.
            DO 1205 IJ=1,NUME
            I = INDI(IJ)
            J = INDJ(IJ)
C           EXTRA TERMS INVOLVING ONE-ELECTRON GENERATORS.
            DO 1170 K=1,ICI1
            ITEMP = IPSCIF (MU,NU,I,J,K,K)
            DO 1150 NF=1,IFPS
            IF(ITEMP.EQ.ITMP(NF)) THEN
C-AK           FTMP(NF) = FTMP(NF) - E(IJ)*TWO
               FTMP(NF) = FTMP(NF) - E(IJ)*OCC(K)
               GO TO 1155
            ENDIF
 1150       CONTINUE
            IFPS = IFPS+1
C-AK        FTMP(IFPS) = -E(IJ)*TWO
            FTMP(IFPS) = -E(IJ)*OCC(K)
            ITMP(IFPS) = ITEMP
 1155       CONTINUE
            ITEMP = IPSCIF (MU,NU,I,K,K,J)
            DO 1160 NF=1,IFPS
            IF(ITEMP.EQ.ITMP(NF)) THEN
C-AK           FTMP(NF) = FTMP(NF) + E(IJ)
               FTMP(NF) = FTMP(NF) + E(IJ)*OCC(K)*PT5
               GO TO 1165
            ENDIF
 1160       CONTINUE
            IFPS = IFPS+1
C-AK        FTMP(IFPS) = E(IJ)
            FTMP(IFPS) = E(IJ)*OCC(K)*PT5
            ITMP(IFPS) = ITEMP
 1165       CONTINUE
 1170       CONTINUE
C-AK        NO EXTRA COULOMB AND EXCHANGE TERMS HERE.
 1205       CONTINUE
         ENDIF
C ***    REMOVE ZERO COEFFICIENTS.
         NEW = 0
         DO 210 NF=1,IFPS
         IF(ABS(FTMP(NF)).GT.TINY) THEN
            NEW = NEW+1
            FTMP(NEW) = FTMP(NF)
            ITMP(NEW) = ITMP(NF)
         ENDIF
  210    CONTINUE
C ***    SAVE NONZERO COEFFICIENTS.
         IF(MODGUG.EQ.2) THEN
            DO 220 NF=1,NEW
            FCOF(NFPS+NF)  = FTMP(NF)
            IFCOF(NFPS+NF) = ITMP(NF)
  220       CONTINUE
         ENDIF
         NFPS = NFPS+NEW
         IF(MUNUF.EQ.0) GO TO 240
         MUNU = MUNUF
         IFPS = 0
      ENDIF
C-AK  CALL PSCIFU (IFGEN(NF1),MU,NU,I,J,K,L)
      IFPS = IFPS+1
      FTMP(IFPS) = FGEN(NF1)
C     FTMP(IFPS) = FGEN(NF1)*PT5
      ITMP(IFPS) = IFGEN(NF1)
  100 CONTINUE
  230 CONTINUE
  240 CONTINUE
C *** SAVE NUMBER OF NONZERO ONE-ELECTRON CI COEFFICIENTS.
      IF(MODGUG.EQ.1) THEN
         LMFPS = NFPS
      ENDIF
C *** DIVIDE TWO-ELECTRON COEFFICIENTS BY THEIR SYMMETRY NUMBER
C     TO CONFORM WITH THE CONVENTIONS IN THE GRADIENT CODE.
      IF(MODGUG.EQ.2) THEN
         DO 250 NF=1,NFPS
         CALL PSCIFU (IFCOF(NF),MUF,NUF,I,J,K,L)
         IF(I.NE.J) THEN
            IF(K.NE.L) THEN
               IF(I.NE.K .OR. J.NE.L) THEN
C                 (IJ,KL)
                  FCOF(NF) = FCOF(NF)/EIGHT
               ELSE
C                 (IJ,IJ)
                  FCOF(NF) = FCOF(NF)/FOUR 
               ENDIF
            ELSE
C              (IJ,KK)
               FCOF(NF) = FCOF(NF)/FOUR
            ENDIF
         ELSE
            IF(K.NE.L) THEN
C              (II,KL)
               FCOF(NF) = FCOF(NF)/FOUR
            ELSE
               IF(I.NE.K) THEN
C                 (II,KK)
                  FCOF(NF) = FCOF(NF)/TWO
C                 (II,II) UNCHANGED
               ENDIF  
            ENDIF  
         ENDIF  
  250    CONTINUE
      ENDIF  
C
C *** COPY CI VECTOR.
      IF(MODGUG.EQ.2) THEN
         DO 300 I=1,LMVEC
         VEC(I) = CIVEC(I)
  300    CONTINUE
      ENDIF
      RETURN
      END
C
      SUBROUTINE PSGUGP(IF,F,IE,E,VEC)
C
C   Prints coefficients and configurations vector for GUGA-CI.
C
C   Parameters:
C
C      IF     - (Output) Packed indices of the F coefficients.
C      F      - (Output) F coefficients.
C      IE     - (Output) Packed indices of the E coefficients.
C      E      - (Output) E coefficients.
C      VEC    - (Output) Configurations vector.
C
C   Accessed common blocks:
C
C      PSDCIP - Various CI parameters.
C      CCIS   - CI state we are interested in.
C      XCI    - CI configurations vectors.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      We write matrix elements of the CI Hamiltomian as:
C
C            active  lm             active  lm
C      H   =   Sum  f     (pq|rs) +   Sum  e   eps
C       lm    pqrs   pqrs              r    r     r
C
C      where (pq|rs) are ERIs over molecular orbitals, and
C      eps are the orbital energies of the MO in the CI-active
C      space.
C
C      f and e coefficient matrices are sparse, and are
C      represented by a (packed) index plus a coefficient
C      value.
C
C   Bugs:
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
C
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      COMMON
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSDCIP/, /PSPRTF/, /PSPRT /
      DIMENSION IF(NCIFS), F(NCIFS), IE(NCIES), E(NCIES), VEC(NCIC)
C     PRINT ONE-ELECTRON COEFFICIENTS.
      WRITE(NB6,500) NCIES
      DO 10 NE=1,NCIES
      CALL PSCIFU (IE(NE),MU,NU,I,J,K,L)
      WRITE(NB6,600) NE,E(NE),IE(NE),MU,NU,I,J,K,L
   10 CONTINUE
C     PRINT TWO-ELECTRON COEFFICIENTS.
      WRITE(NB6,510) NCIFS
      DO 20 NF=1,NCIFS
      CALL PSCIFU (IF(NF),MU,NU,I,J,K,L)
      WRITE(NB6,610) NF,F(NF),IF(NF),MU,NU,I,J,K,L
   20 CONTINUE
C     PRINT CI-VECTOR.
      WRITE(NB6,520) NCIC
      WRITE(NB6,530) (VEC(I),I=1,NCIC)
      RETURN
  500 FORMAT(//,' PSGUGP: ONE-ELECTRON COFFICIENTS.',I10/)
  510 FORMAT(//,' PSGUGP: TWO-ELECTRON COFFICIENTS.',I10/)
  520 FORMAT(//,' PSGUGP: CI-EIGENVECTORS, DIMENSION:',I10/)
  530 FORMAT(1X,5F20.10)
  600 FORMAT(' E',I8,F10.5,I10,6I5)
  610 FORMAT(' F',I8,F10.5,I10,6I5)
      END
C
      INTEGER FUNCTION IPSCIF(L,M,IP,IQ,IR,IS)
C
C   Computes a packing index for CI F coefficients.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L,M    - Indices of the CI matrix element.
C      IP,IQ,
C      IR,IS  - MO indices.
C
C   Accessed common blocks:
C
C      PSDCIP - We need size of the two-particle density 
C               from here.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C      Integer overflow is not expected to happen.
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
C
      COMMON
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
     ./PSPRT / NB6
      SAVE /PSDCIP/, /PSPRT /
      IPACK(I,J) = (J*(J-1))/2 + I
C    Parameter checking.
      IF( L .LE.0 .OR. L .GT.NCIC .OR. M .LE.0 .OR. M .GT.NCIC .OR. 
     .    IP.LE.0 .OR. IP.GT.NCIO .OR. IQ.LE.0 .OR. IQ.GT.NCIO .OR. 
     .    IR.LE.0 .OR. IR.GT.NCIO .OR. IS.LE.0 .OR. IS.GT.NCIO ) THEN
          WRITE(NB6,10100) L,M,IP,IQ,IR,IS
          STOP 'IPSCIF'
      ENDIF
C
      I1     = IPACK(MIN(IP,IQ),MAX(IP,IQ))
      I2     = IPACK(MIN(IR,IS),MAX(IR,IS))
      INDX   = IPACK(MIN(I1,I2),MAX(I1,I2))
      INDY   = IPACK(MIN(L,M),  MAX(L,M))
      IPSCIF = INDY*NCIGAM + INDX
      RETURN
10100 FORMAT(' IPSCIF CALLED WITH WRONG ARGUMENTS: ',6(1X,I7))
      END
C
      SUBROUTINE PSCIFU(IND,L,M,IP,IQ,IR,IS)
C
C   Converts packing index to individual contributions.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IND    - Packing index
C      L,M    - (Output) indices of the CI matrix element.
C               L <= M
C      IP,IQ,
C      IR,IS  - (Output) MO indices. IP <= IQ, IR <= IS,
C               (IP,IQ) <= (IR,IS)
C
C   Accessed common blocks:
C
C      PSDCIP - We need size of the two-particle density 
C               from here.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
C      Integer overflow is not expected to happen.
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
      PARAMETER (TWO=2.D0)
      PARAMETER (C3O4=0.75D0)
C
      COMMON
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
      SAVE /PSDCIP/
C
      INDX   = MOD(IND-1,NCIGAM) + 1
      INDY   = (IND-1)/NCIGAM
      M      = NINT(SQRT(TWO*INDY-C3O4))
      L      = INDY - (M*(M-1))/2
      IRS    = NINT(SQRT(TWO*INDX-C3O4))
      IPQ    = INDX - (IRS*(IRS-1))/2
      IS     = NINT(SQRT(TWO*IRS-C3O4))
      IR     = IRS - (IS*(IS-1))/2
      IQ     = NINT(SQRT(TWO*IPQ-C3O4))
      IP     = IPQ - (IQ*(IQ-1))/2
      RETURN
      END
C
      INTEGER FUNCTION IPSCIG(IP,IQ,IR,IS)
C
C   Computes a packing index for the two-particle density matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IP,IQ,
C      IR,IS  - MO indices.
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
C      Integer overflow is not expected to happen.
C
      USE LIMIT, ONLY: LMX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /INDEX / INDX(LMX)
      IPACK(I,J) = INDX(J) + I
      I1     = IPACK(MIN(IP,IQ),MAX(IP,IQ))
      I2     = IPACK(MIN(IR,IS),MAX(IR,IS))
      I3     = MIN(I1,I2)
      I4     = MAX(I1,I2)
      IPSCIG = (I4*(I4-1))/2 + I3
C
C     OLD CODE
C     IPACK(I,J) = (J*(J-1))/2 + I
C     I1     = IPACK(MIN(IP,IQ),MAX(IP,IQ))
C     I2     = IPACK(MIN(IR,IS),MAX(IR,IS))
C     IPSCIG = IPACK(MIN(I1,I2),MAX(I1,I2))
C
      RETURN
      END
C
C
      SUBROUTINE PSCIGA(IF,F,VEC,GAM)
C
C   Computes two-particle density matrix for a given CI state.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IF     - Packed indices of the F coefficients.
C      F      - F coefficients.
C      VEC    - CI state vector.
C      GAM    - (Output) two-particle density, in the packed form.
C
C   Accessed common blocks:
C
C      PSDCIP - Various CI parameters.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Two-particle density matrix is given by:
C
C               active       jk
C      GAM    =   Sum  A  A f
C         pqrs    j,k   j  k pqrs
C
C   Possible optimizations:
C
C      Precomputing A(J) A(K) products will simplify the inner
C      loop considerably.
C
C   Bugs:
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
      PARAMETER (TWO =2.D0)
C
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      COMMON
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSDCIP/, /PSPRTF/, /PSPRT /
      DIMENSION IF(NCIFS), F(NCIFS), VEC(NCIC), GAM(NCIGAM)
      EXTERNAL IPSCIG
C    Debugging output.
      IF(LPCI) THEN
          WRITE(NB6,11010)
          WRITE(NB6,11020) (IF(I), F(I), I=1,NCIFS)
      ENDIF
C    Zero out two-particle density first.
      CALL PSZRV(NCIGAM,GAM)
C
      DO 1000 I=1,NCIFS
          CALL PSCIFU(IF(I),L,M,IP,IQ,IR,IS)
          IPQRS = IPSCIG(IP,IQ,IR,IS)
          TEMP  = F(I)*VEC(L)*VEC(M)
          IF( L.NE.M ) TEMP = TEMP*TWO
          GAM(IPQRS) = GAM(IPQRS) + TEMP
 1000 CONTINUE
C
      RETURN
11010 FORMAT(/' WEIGHTS OF TWO-ELECTRON INTEGRALS IN CI HAMILTONIAN: '/)
11020 FORMAT((4(1X,I10,1X,G14.8)))
      END
C
      SUBROUTINE PSCILG(IE,E,VEC,DLT)
C
C   Computes lagrangian vector for a given CI state.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IE     - Packed indices of the E coefficients.
C      E      - F coefficients.
C      VEC    - CI state vector.
C      DLT    - (Output) Lagrangian for the CI-active orbitals.
C
C   Accessed common blocks:
C
C      PSDCIP - Various CI parameters.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Lagrangian multipliers are given by:
C
C             active       jk
C       /\  =   Sum  A  A e
C         r     j,k   j  k r
C
C   Possible optimizations:
C
C      Precomputing A(J) A(K) products will simplify the inner
C      loop considerably.
C
C   Bugs:
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
      PARAMETER (TWO   =2.D0)
C
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      COMMON
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
      SAVE /PSDCIP/, /PSPRTF/, /PSPRT /
      DIMENSION IE(NCIFS), E(NCIFS), VEC(NCIC), DLT(NCIO)
C    Debugging output.
      IF(LPCI) THEN
          WRITE(NB6,11010)
          WRITE(NB6,11020) (IE(I), E(I), I=1,NCIES)
      ENDIF
C    Zero out lagrangian multipliers first.
      CALL PSZRV(NCIO,DLT)
C
      DO 1000 I=1,NCIES
          CALL PSCIFU(IE(I),L,M,IP,IQ,IR,IS)
          TEMP  = E(I)*VEC(L)*VEC(M)
          IF( L.NE.M ) TEMP = TEMP*TWO
          DLT(IS) = DLT(IS) + TEMP
 1000 CONTINUE
C
      RETURN
11010 FORMAT(/' WEIGHTS OF ORBITAL ENERGIES IN CI HAMILTONIAN: '/)
11020 FORMAT((4(1X,I10,1X,G14.8)))
      END
C
      SUBROUTINE PSGUG1(ONEDEN,DLT)
C
C   Coded by: Axel Koslowski
C
C   Parameters:
C
C      ONEDEN - (Input)  One-particle density matrix for the
C                        CI-active orbitals.
C      DLT    - (Output) One-particle density matrix diagonal
C                        for the CI-active orbitals.
C
C   Accessed common blocks:
C
C      PSDCIP - Various CI parameters.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Copy the matrix diagonal into the vector.
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
      COMMON
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
      SAVE /PSDCIP/
      DIMENSION ONEDEN(LMACT,NCIO), DLT(NCIO)
      DO I=1,NCIO
         DLT(I) = ONEDEN(I,I)
      END DO
      RETURN
      END
C
      SUBROUTINE PSGUG2(TWODEN,GAM)
C
C   Coded by: Axel Koslowski
C
C   Parameters:
C
C      TWODEN - (Input)  two-particle density, in the packed form.
C      GAM    - (Output) two-particle density, in the packed form.
C
C   Accessed common blocks:
C
C      PSDCIP - Various CI parameters.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Divide matrix elements by their symmetry number while copying.
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
      COMMON
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
      SAVE /PSDCIP/
      DIMENSION TWODEN(NCIGAM), GAM(NCIGAM)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WEIGHT
      ALLOCATE(WEIGHT(NCIGAM))
      WEIGHT = 0.D0
      DO I=1,NCIO
         DO J=1,NCIO
            DO K=1,NCIO
               DO L=1,NCIO
                  IJKL = IPSCIG(I,J,K,L)
                  WEIGHT(IJKL) = WEIGHT(IJKL) + 1.D0
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      GAM = TWODEN / WEIGHT
      DEALLOCATE(WEIGHT)
      RETURN
      END
C
      SUBROUTINE PRDEN1(ONEDEN)
C
C   Coded by: Axel Koslowski
C
C   Parameters:
C
C      ONEDEN - (Input)  One-particle density matrix for the
C                        CI-active orbitals.
C
C   Accessed common blocks:
C
C      PSDCIP - Various CI parameters.
C      PSPRT  - Contains unit number of standard output.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Print the one-particle density matrix for the CI-active orbitals.
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
      COMMON
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
     ./PSPRT / NB6
      SAVE /PSDCIP/
      DIMENSION ONEDEN(LMACT,NCIO)
      WRITE(NB6,500)
C     Loop over column blocks.
      DO K=1,NCIO,10
         L = MIN(K+9, NCIO)
C        Print column indices.
         WRITE(NB6,510) (J, J=K,L)
C        Loop over lines and print.
         DO I=1,NCIO
            WRITE(NB6,520) I, ONEDEN(I, K:L)
         END DO
      END DO
      RETURN
  500 FORMAT(//1X,'ONE-PARTICLE DENSITY MATRIX PASSED TO PS MODULE.')
  510 FORMAT(/ 2X,10I12)
  520 FORMAT(  1X,I4,10F12.6)
      END
C
      SUBROUTINE PRDEN2(TWODEN)
C
C   Coded by: Axel Koslowski
C
C   Parameters:
C
C      TWODEN - (Input)  Two-particle density matrix for the
C                        CI-active orbitals.
C
C   Accessed common blocks:
C
C      PSDCIP - Various CI parameters.
C      PSPRT  - Contains unit number of standard output.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Print the lower triangle of the two-particle density matrix
C      for the CI-active orbitals.
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
      COMMON
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
     ./PSPRT / NB6
      SAVE /PSDCIP/
      DIMENSION TWODEN(NCIGAM), IROW(10), ICOL(10)
      WRITE(NB6,500)
C     Matrix dimension.
      NDIM = MATDIM(NCIGAM)
C     Loop over column blocks.
      DO K=1,NDIM,10
         L = MIN(K+9, NDIM)
C        Split column indices.
         DO I=K,L
            CALL SPLIT1(I, IROW(I-K+1), ICOL(I-K+1))
         END DO
C        Print column index components.
         WRITE(NB6,510) (IROW(I), I=1, L-K+1)
         WRITE(NB6,520) (ICOL(I), I=1, L-K+1)
C        Linear index of leftmost diagonal element.
         LINEAR = K*(K+1)/2
C        Print line by line.
         DO I=K, NDIM
            M = MIN(I-K, 9)
            CALL SPLIT1(I, IR, IC)
            WRITE(NB6,530) IR, IC, (TWODEN(LINEAR+J), J=0, M)
            LINEAR = LINEAR + I
         END DO
      END DO
      RETURN
  500 FORMAT(//1X,'TWO-PARTICLE DENSITY MATRIX PASSED TO PS MODULE.')
  510 FORMAT(/ 8X,10I12)
  520 FORMAT(  8X,10I12)
  530 FORMAT(  1X,2I5,10F12.6)
      END
C
      SUBROUTINE SPLIT1(N, I, J)
C
C   Coded by: Axel Koslowski
C
C   Parameters:
C
C      N      - (Input)  Linear index.
C      I      - (Output) Row index.
C      J      - (Output) Column index.
C
C   Accessed common blocks:
C
C      PSPRT  - Contains unit number of standard output.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Split linear index of the lower triangular matrix
C      into row and column indices (all indices 1-based).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./PSPRT / NB6
      ILOWER = CEILING(0.5D0 * SQRT(8.0D0 * N + 1.0D0) - 0.5D0)
      IUPPER = FLOOR  (0.5D0 * SQRT(8.0D0 * N - 7.0D0) + 0.5D0)
      IF(ILOWER.NE.IUPPER) THEN
         WRITE(NB6,10000) ILOWER, IUPPER, N
         STOP 'SPLIT1'
      ENDIF
      I = ILOWER
      J = N - I * (I-1) / 2
      RETURN
10000 FORMAT(1X,'CONTRADICTING ROW INDICES',I4,' AND',I4,' FOR N=',I8)
      END
C
      INTEGER FUNCTION MATDIM(N)
C
C   Coded by: Axel Koslowski
C
C   Parameters:
C
C      N      - (Input)  Linear size.
C
C   Accessed common blocks:
C
C      PSPRT  - Contains unit number of standard output.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Calculate dimension of triangular matrix of linear size N.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./PSPRT / NB6
      M = NINT(0.5D0 * SQRT(8.0D0 * N + 1.0D0) - 0.5D0)
      IF(2*N.NE.M*(M+1)) THEN
         WRITE(NB6,10000) N
         STOP 'MATDIM'
      ENDIF
      MATDIM = M
      RETURN
10000 FORMAT(1X,'NOT A TRIANGULAR MATRIX FOR N=',I8)
      END
C
      SUBROUTINE PSCIQE(CALP,LDC,QEX,LDQ,DUMP)
C
C   Compute redundant right-hand side for the CI Z-vector CPHF.
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
C      PSDGBL - Global computation parameters.
C      PSDYNM - Dynamic memory handles.
C      PSDCIP - CI parameters.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      QEX is computed according to the expression:
C
C               Active
C       q   = 4   Sum   Gamma     (rs|qi)
C        ip       rsq        pqrs
C
C               Active Active
C           = 4   Sum    Sum  Gamma     (2 - delta  ) (rs|qi)
C                  q    r<=s       pqrs           rs
C
C      The code is largerly snatched from PSHLQE.
C      We also tuck away all ERIs over CI-active MOs, so as
C      to run a check on the CI state vector later.
C
C   Bugs:
C
C      Some of the integrals are evaluated twice.
C
C      Contraction with two-particle density matrix doen't
C      use symmetry to a full extent.
C
C      A lot of the ERIs over CI-active MOs are stored twice.
C
      USE LIMIT, ONLY: LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXCIO=LMACT)
C
      PARAMETER (ZERO=0.0D0)
      PARAMETER (TWO =2.0D0)
      PARAMETER (FOUR=4.0D0)
C
      EXTERNAL IPSCIG
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
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
      SAVE /PSDGBL/, /PSDCIP/
C
      DIMENSION CALP(LDC,NORBS), QEX(LDQ,NCIO), DUMP(*)
C
      IINTS  = IPSMOF(LAI2T5)
      ICIORB = IPSMOF(LCIORB)
      ICIGAM = IPSMOF(LCIGAM)
      ICIINT = IPSMOF(LCIINT)
C    Zero output array
      CALL PSZRM(NORBS,NCIO,QEX,LDQ)
C
      DO 1000 IS=1,NCIO
          CALL PST1ST(DUMP,DUMP(ICIORB+(IS-1)*NORBS))
          DO 900 IR=1,IS
              SC = FOUR
              IF( IR.NE.IS ) SC = SC * TWO
              CALL PST2ND(DUMP,DUMP(ICIORB+(IR-1)*NORBS),1)
              CALL PST3RD(DUMP,1,DUMP(ICIORB),NORBS,NCIO)
              CALL PST4TH(DUMP,CALP,LDC,NCIO,NORBS)
              DO 400 IP=1,NCIO
                  IPTRUE = INDCIO(IP)
                  DO 300 IQ=1,IP
                      IPQRS  = IPSCIG(IP,IQ,IR,IS)
                      DUMP(ICIINT+(IPQRS-1)) = 
     .                    DUMP(IINTS+(IQ-1)*NORBS+(IPTRUE-1))
  300             CONTINUE
  400         CONTINUE
              DO 800 IP=1,NCIO
                  DO 700 I=1,NORBS
                      SUM = ZERO
                      DO 600 IQ=1,NCIO
                          IPQRS = IPSCIG(IP,IQ,IR,IS)
                          SUM = SUM + DUMP(ICIGAM+(IPQRS-1))
     .                              * DUMP(IINTS+(IQ-1)*NORBS+(I-1))
  600                 CONTINUE
                      QEX(I,IP) = QEX(I,IP) + SC*SUM
  700             CONTINUE
  800         CONTINUE
  900     CONTINUE
 1000 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSCIQW(QEX,LDQ,Q,DLT)
C
C   Convert redundant set of CI derivatives contraction coefficients 
C   into (quasi-) minimal one.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      QEX    - Redundant set of contraction coefficients.
C               Destroyed on return.
C      LDQ    - Leading dimension of the QEX matrix.
C      Q      - (Output) non-redundant coefficients.
C      DLT    - Lagrangian multipliers.
C
C   Accessed common blocks:
C
C      PSDGBL - General computation parameters
C      PSDGB2 - Response-related parameters
C      PSOCC  - Occupation blocks description
C      PSDCIP - CI parameters
C
C   Modified common blocks:
C
C   Local storage:
C
C      LMX INTEGER's.
C
C   Module logic:
C
C      The minimal set of the contraction coefficients is given by:
C
C      Q   = /\             if p is CI-active
C       pp     p
C
C      Q   = q              if i is CI-inactive, p>=i is CI-active
C       ip    ip
C
C      Q   = -q             if i is CI-inactive, p<i is CI-active
C       pi     ip
C
C      Q   = q  - q         if p<r are CI-active
C       pr    pr   rp
C
C      Non-redundant coefficients are in the standard occupation block
C      order. Redundant coefficients is a rectangular matrix of the
C      dimension (NORBS,NCIO). Entries in the redundant part corresponding
C      to non-redundant variables are zero, so there is no double counting.
C      (p,r) entries in the redundant part are stored twice, as (p,r)
C      and (r,p). However, only one of the entries (with p<r) is non-zero.
C
C      QEX also served as a scoreboard, since all processed non-redundant
C      coefficients are zeroed out.
C
C   Bugs:
C
C      Several of the "redundant" entries (all of which are zero) need
C      not be here. The results are still correct, though.
C
C      The code crawls through the entire Q unnecessarily.
C
      USE LIMIT, ONLY: LM1, LMX, LMACT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXGRP=10)
      PARAMETER (MAXCIO=LMACT)
C
      PARAMETER (ZERO=0.0D0)
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
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
      SAVE /PSOCC/, /PSDGBL/, /PSDGB2/, /PSDCIP/
      DIMENSION QEX(LDQ,NCIO), Q(*), DLT(NCIO)
      INTEGER IBACK(LMX)
C    Initialize CI activity vector
      DO 100 I=1,NORBS
          IBACK(I) = 0
  100 CONTINUE
      DO 110 I=1,NCIO
          IBACK(INDCIO(I)) = I
  110 CONTINUE
C    Coefficients for non-redundant CPHF variables. This one is a
C    bit convoluted since CI-active orbitals need not be continuous.
      DO 600 IBL1=1,NMGRPA-1
          IND1 = IG1STA(IBL1)
          N1   = IGCNTA(IBL1)
          DO 500 IBL2=IBL1+1,NMGRPA
              IND2 = IG1STA(IBL2)
              N2   = IGCNTA(IBL2)
              IPT  = IGBASA(IBL1,IBL2)
              DO 400 I1=1,N1
                  IORB1 = IND1 + (I1-1)
                  ICI1  = IBACK(IORB1)
                  DO 300 I2=1,N2
C                    By construction, IORB1 is always less than IORB2
                      IORB2  = IND2 + (I2-1)
                      ICI2   = IBACK(IORB2)
                      Q(IPT) = ZERO
                      IF( ICI1.NE.0 ) THEN
                          Q(IPT) = Q(IPT) - QEX(IORB2,ICI1)
                          QEX(IORB2,ICI1) = ZERO
                      ENDIF
                      IF( ICI2.NE.0 ) THEN
                          Q(IPT) = Q(IPT) + QEX(IORB1,ICI2)
                          QEX(IORB1,ICI2) = ZERO
                      ENDIF
                      IPT = IPT + 1
  300             CONTINUE
  400         CONTINUE
  500     CONTINUE
  600 CONTINUE
C    And now for active redundant CPHF variables. First, copy
C    all the remaining coefficients to the redundant part.
      IRED = IQSZ + 1
      CALL DCOPY(NCIO*NORBS,QEX,1,Q(IRED),1)
C    Adjust signs of the coefficients with i>p and replace
C    diagonal entries with corresponding Lagrangian multipliers.
      DO 900 I=1,NORBS
          DO 800 IP=1,NCIO
              IPREAL = INDCIO(IP)
              IF( I.GT.IPREAL ) THEN
                  Q(IRED+(IP-1)*NORBS+(I-1)) = 
     .                               -Q(IRED+(IP-1)*NORBS+(I-1))
              ENDIF
              IF( I.EQ.IPREAL ) THEN
                  Q(IRED+(IP-1)*NORBS+(I-1)) = -DLT(IP)
              ENDIF
  800     CONTINUE
  900 CONTINUE
C    Collapce active-active entries. It doesn't really matter which
C    one we are using, since both are treated alike.
      DO 1100 IP=1,NCIO
          IPREAL = INDCIO(IP)
          DO 1000 IQ=1,IP-1
              IQREAL = INDCIO(IQ)
              Q(IRED+(IP-1)*NORBS+(IQREAL-1)) = 
     .              Q(IRED+(IP-1)*NORBS+(IQREAL-1)) 
     .            + Q(IRED+(IQ-1)*NORBS+(IPREAL-1))
              Q(IRED+(IQ-1)*NORBS+(IPREAL-1)) = ZERO
 1000     CONTINUE
 1100 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSCICK(EALP,IF,F,IE,E,VEC,HCI,VEC1,AINT)
C
C   Make sure CI state vector is an eigenvector of the CI hamiltonian.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      EALP   - Orbital energies.
C      IF     - Control indices for the F coefficients.
C      F      - F coefficients (weights of ERIs in CI Hamiltonian)
C      IE     - Control indices for the E coefficients.
C      E      - E coefficients (weights of orbital energies in
C               CI Hamiltonian).
C      VEC    - CI state vector.
C      HCI    - Scratch pad for the CI Hamiltonian.
C      VEC1   - HCI*VEC.
C      AINT   - Packed ERIs over CI-active MOs.
C
C   Accessed common blocks:
C
C      PSDOPT - IPRINT is needed from here
C      PSDGBL - General computation parameters
C      PSDGB2 - Response-related parameters
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
      PARAMETER (AU   = 3.67511944138184490995957368D-2)
      PARAMETER (MAXCIO=LMACT)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (ONE =1.0D0)
      PARAMETER (EPSA=1.D-12)
      PARAMETER (EPSB=1.D-04)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
     ./PSDOPT/ DSTORE, DSTEP,  DECONV, DPCONV, DCPHF, DPREC, DCDIFF,
     .         DSHIFT, DBASCR, DNCOFF,
     .         IUMIX,  IURHS,  IUK,    IURES,
     .         IPRINT, IENRG,  ICORE,  IDISK,  IMIX,   IDENS,
     .         INDSYM, IQSWAP, IAVEIT, IMAXIT, INRHS,  IKRVEC,
     .         IROWS,  IPRECT, INCPUS, IDSTRP, IHLST,  IHLWRP,
     .         IKMODE, ISOLVE, IKRSAV, NMRLEV, INTCTL, ICIOPT
      SAVE /PSDGBL/, /PSDCIP/, /PSPRTF/, /PSDOPT/, /PSPRT /
      DIMENSION EALP(NORBS), IF(NCIFS), F(NCIFS), IE(NCIES), E(NCIES), 
     .          VEC(NCIC), HCI(NCIC,NCIC), VEC1(NCIC), AINT(NCIGAM)
      EXTERNAL IPSCIG
C    Multiply ERIs by their symmetry numbers
      CALL PSSCTN(NCIO,AINT)
C    Zero out the H*V product vector
      VEC1(1:NCIC) = 0.D0
C    Form the product H*V
C      ... ERI contributions go first
      DO 1000 I=1,NCIFS
          CALL PSCIFU(IF(I),L,M,IP,IQ,IR,IS)
          IPQRS    = IPSCIG(IP,IQ,IR,IS)
C         HCI(L,M) = HCI(L,M) + F(I)*AINT(IPQRS)
          HELEM    = F(I)*AINT(IPQRS)
          VEC1(L)  = VEC1(L) + HELEM * VEC(M)
          IF(L.NE.M)  VEC1(M) = VEC1(M) + HELEM * VEC(L)
 1000 CONTINUE
C      ... now, the orbital energies contributions
      DO 2000 I=1,NCIES
          CALL PSCIFU(IE(I),L,M,IP,IQ,IR,IS)
          ISREAL = INDCIO(IS)
C         HCI(L,M) = HCI(L,M) + E(I)*EALP(ISREAL)*AU
          HELEM    = E(I)*EALP(ISREAL)*AU
          VEC1(L)  = VEC1(L) + HELEM * VEC(M)
          IF(L.NE.M)  VEC1(M) = VEC1(M) + HELEM * VEC(L)
 2000 CONTINUE
C    The product vector H*V is now in shape
      EVAL = DDOT(NCIC,VEC,1,VEC1,1)
      CALL DAXPY(NCIC,-EVAL,VEC,1,VEC1,1)
      RESD = DNRM2(NCIC,VEC1,1)
      IF( IPRINT.GE.2 ) WRITE(NB6,11020) EVAL, RESD
      IF(LPCI) THEN
          WRITE(NB6,11030)(I,VEC(I),VEC1(I),I=1,NCIC)
      ENDIF
      IF( RESD.GE.EPSB ) THEN
          WRITE(NB6,11040) 'FATAL', RESD
          STOP 'PSCICK'
      ELSE IF( RESD.GE.EPSA .AND. IPRINT.GT.-5 ) THEN
          WRITE(NB6,11040) 'WARNING', RESD
      ENDIF
C    This is pure state all right.
      RETURN
11010 FORMAT(/' THE CI HAMILTONIAN MATRIX IS:'/)
11020 FORMAT(' CI STATE IS ',G15.8,' A.U. ABOVE REFERENCE POINT.'
     .      /' |H|V> - <V|H|V>V RESIDUE IS ',G15.8)
11030 FORMAT(/' CONFIGURATION ',5X,' WEIGHT ',4X,4X,' RESIDUE '/
     .       (5X,I10,           2X,G15.9,        2X,G15.9))
11040 FORMAT(1X,A,': CI STATE VECTOR DEVIATES FROM A PURE STATE BY ',
     .       G15.9,' IN PSCICK.')
      END
C
      SUBROUTINE PSCIUG(XR,EALP)
C
C   Scale redundant CPHF variables by the difference in orbital
C   energies. CI version.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      XR     - Variables to scale, two-dimensional array.
C               First dimension: MO index, 1 to NORBS
C               Second dimension: CI-active MO index, 1 to NCIO.
C               "diagonal" entries are lagrangian multipliers, and
C               should be left alone for now.
C      EALP   - Orbital energies, EVs
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options.
C      PSDGBL - Number of MOs is taken from here.
C      PSDCIP - CI parameters.
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
      PARAMETER (MAXCIO=LMACT)
C
      PARAMETER (AU   = 3.67511944138184490995957368D-2)
      PARAMETER (TOL  = 1.0D-5)
      PARAMETER (QTOL = 3.0D-8)
      PARAMETER (ZERO = 0.0D0)
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
C
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
     ./PSDCIP/ NCIO, NCIC, ICIMOD, NCIGAM, NCIFS, NCIES, INDCIO(MAXCIO)
     ./PSPRT / NB6
      SAVE /PSDOPT/, /PSDGBL/, /PSDCIP/, /PSPRT /
C
      DIMENSION XR(NORBS,NCIO), EALP(*)
C
      DO 900 J=1,NCIO
          JREAL = INDCIO(J)
          DO 800 IREAL=1,NORBS
C            Skip Lagrangian multipliers
              IF( IREAL.EQ.JREAL ) GOTO 800
C            If a particular coupling coefficient is inactive, we
C            don't want to sound an alarm over it.
              IF( XR(IREAL,J).NE.ZERO ) THEN
                  DE = EALP(IREAL) - EALP(JREAL)
                  IF( JREAL.LT.IREAL ) DE = -DE
                  IF( ABS(DE).LE.TOL ) THEN
                      IF( ABS(XR(IREAL,J)) .GE. QTOL ) THEN
                          WRITE(NB6,10010)
     @                      IREAL, JREAL, TOL, XR(IREAL,J), QTOL
                          STOP 'PSCIUG'
                      ENDIF
                      IF(IPRINT.GE.1)
     @                  WRITE(NB6,10020) IREAL, JREAL, XR(IREAL,J)
                      XR(IREAL,J) = ZERO
                  ELSE
                      IF(IPRINT.GE.5)
     @                  WRITE(NB6,10030) IREAL, JREAL, XR(IREAL,J)
                      XR(IREAL,J) = XR(IREAL,J) / (AU*DE)
                  ENDIF
              ENDIF
  800     CONTINUE
  900 CONTINUE
C
      RETURN
10000 FORMAT(' WARNING: ENERGY SEPARATION FOR ACTIVE REDUNDANT PAIR ',
     .       I4, ',', I4, ' IS LESS THAN ', G10.5, ' EV'/
     .       ' ATTEMPTING TO EXCLUDE PAIR AND CONTINUE - DERIVATIVES ',
     .       'ARE IN DOUBT.' )
10010 FORMAT(' ERROR: ENERGY SEPARATION FOR ACTIVE REDUNDANT PAIR ',
     .       I4, ',', I4, ' IS LESS THAN ', G10.5, ' EV,'/
     .       ' BUT Q COEFFICIENT IS ', E12.5,
     .       ' WHICH IS GREATER THAN ', G10.5, '.')
10020 FORMAT(' Q COEFFICIENT FOR ACTIVE REDUNDANT PAIR ',
     .       I4, ',', I4, ' OF ', E12.5, ' HAS BEEN NEGLECTED.')
10030 FORMAT(' Q COEFFICIENT FOR ACTIVE REDUNDANT PAIR ',
     .       I4, ',', I4, ' IS ', E12.5, '.')
      END
