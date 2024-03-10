C     ******************************************************************
C
C     NMR 3-center integrals over Slater AOs: Low-level routines.
C
C     ******************************************************************
C
C     Evaluation of the expansion coefficients for Sharma's expansion 
C     at the displaced center (Phys. Rev. A, 13, 517 (1976)).
C
      SUBROUTINE PS3SBT(L,MAXM,BETA,LDQ,LDK,LDP)
C
C   Compute block of values for the Sharma's beta coefficients.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L      - Order of the coefficents' orbital moment.
C      MAXM   - Largest (or smallest) projection of the orbital moment. 
C               Projections in the range 0 to MAXM (or MAXM to 0) will
C               be evaluated. 
C      BETA   - (Output) beta coefficeints, four-dimensional
C               array (see the expression below).
C                 First index:  IQ, 0 to L - M - 2*IP
C                 Second index: IK, 0 to L - M - 2*IP
C                 Third index:  IP, 0 to MIN(L,(L-M)/2)
C                 Forth index:  0 to ABS(MAXM)
C               Additionally, only IQ+IK not exceeding L - M - 2*IP will
C               be evaluated. BETA(*,*,*,0) is also used as a temporary 
C               holding place.
C       LDQ   - First leading dimension of BETA minus 1
C       LDK   - Second leading dimension of BETA minus 1
C       LDP   - Third leading dimension of BETA minus 1
C               
C   Accessed common blocks:
C
C       Nope.
C
C   Local storage:
C
C       5*MAXL+3 DOUBLE PRECISION storage units.
C
C   Module logic:
C
C       Beta coefficients are given by:                             P            P-L
C                            (2 L+1) (L-M)!  1/2       1        (-1) (2 L-2 P)! 4
C       BETA(L,M)(P,Q,K) = ( -------------- )    -------------- --------------------
C                                (L+M)!          (L-M-2 P-Q-K)!     P!Q!K!(L-P)!
C
C   Possible optimizations:
C
C   Bugs:
C
C      Maximum order of the coefficients is limited by MAXL.
C      This routine need more precision than REAL*8 can provide, and hence
C      uses non-standard REAL*16. This routine *will* compile and work
C      with REAL*8, but resulting coefficients will be imprecise.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=40)
C    Comment next line if your compiler does not support REAL*16
      REAL*16 BETA, FACT, RFACT, XFACT, ONE, FOUR, SCP, SCK, SCM, FOURI
C
      PARAMETER (ONE=1)
      PARAMETER (FOUR=4)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION FACT(0:2*MAXL), RFACT(0:2*MAXL), XFACT(0:MAXL)
      DIMENSION BETA(0:LDQ,0:LDK,0:LDP,0:*)
C    A bit of error-checking
      IF( ABS(MAXM).GT.L ) THEN
          WRITE(NB6,10010) MAXM, L
          STOP 'PS3SBT'
      ENDIF
      IF( L.GT.MAXL ) THEN
          WRITE(NB6,10020) L, MAXL
          STOP 'PS3SBT'
      ENDIF
      MAXQ = L + (ABS(MAXM)-MAXM)/2
      IF( LDQ.LT.MAXQ ) THEN
          WRITE(NB6,10030) LDQ, MAXQ
          STOP 'PS3SBT'
      ENDIF
      IF( LDK.LT.MAXQ ) THEN
          WRITE(NB6,10040) LDK, MAXQ
          STOP 'PS3SBT'
      ENDIF
      IF( LDP.LT.L ) THEN
          WRITE(NB6,10050) LDP, L
          STOP 'PS3SBT'
      ENDIF
C    Precompute factorials and inverse factorials
      FACT (0) = ONE
      RFACT(0) = ONE
      DO 1000 I=1,2*L
          FACT (I) = I*FACT(I-1)
          RFACT(I) = ONE/FACT(I)
 1000 CONTINUE
C    Precompute (2*I)!/I!/4**I
      FOURI = ONE
      DO 1050 I=0,L
          XFACT(I) = FACT(2*I)*RFACT(I)*FOURI
          FOURI    = FOURI/FOUR
 1050 CONTINUE
C    Compute parts of the expression independent of M, and store 'em
C    in the BETA(*,*,*,0).
      DO 2000 IP=0,MIN(L,MAXQ/2)
          SCP = XFACT(L-IP)*RFACT(IP)
          IF( MOD(IP,2).NE.0 ) SCP = -SCP
          MAXQT = MAXQ-2*IP
          DO 1900 IK=0,MAXQT
              SCK = SCP*RFACT(IK)
              DO 1800 IQ=0,MAXQT-IK
                  BETA(IQ,IK,IP,0) = SCK*RFACT(IQ)
 1800         CONTINUE
 1900     CONTINUE
 2000 CONTINUE
C    Do the M-dependent part
      DO 3000 MABS=ABS(MAXM),0,-1
          M   = SIGN(MABS,MAXM)
          SCM = SQRT((2*L+1)*FACT(L-M)*RFACT(L+M))
          DO 2900 IP=0,MIN(L,(L-M)/2)
              MAXQT = L-M-2*IP
              DO 2800 IK=0,MAXQT
                  DO 2700 IQ=0,MAXQT-IK
                      BETA(IQ,IK,IP,MABS) = BETA(IQ,IK,IP,0)*SCM*
     .                    RFACT(L-M-2*IP-IQ-IK)
 2700             CONTINUE
 2800         CONTINUE
 2900     CONTINUE
 3000 CONTINUE
C    Does that look perverse enough?
      RETURN
10010 FORMAT(' MAXM = ',I5,' IS INCOMPATIBLE WITH L = ',I5,' IN PS3SBT')
10020 FORMAT(' L = ',I5,' IS TOO BIG, LIMIT IS ',I5,' IN PS3SBT')
10030 FORMAT(' FACT WILL OVERFLOW ON FIRST DIMENSION IN PS3SBT. ',
     .       'HAVE = ',I5,' NEED = ',I5)
10040 FORMAT(' FACT WILL OVERFLOW ON SECOND DIMENSION IN PS3SBT. ',
     .       'HAVE = ',I5,' NEED = ',I5)
10050 FORMAT(' FACT WILL OVERFLOW ON THIRD DIMENSION IN PS3SBT. ',
     .       'HAVE = ',I5,' NEED = ',I5)
      END
C
      SUBROUTINE PS3SBN(LS,L,BN,LDS,LDN)
C
C   Compute block of values for the Sharma's BN coefficients.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      LS     - Expansion order.
C      L      - Orbital moment of the expanded orbitals
C      BN     - (Output) beta coefficients, three-dimensional
C               array. See expression below.
C                 First index, NU, 0 to L+LS
C                 Second index, IS, 0 to L+LS
C                 Third index, M, 0 to L
C               Only values with IS+NU.LE.L+LS are computed.
C      LDS    - First leading dimension of BN minus 1
C      LDN    - Second leading dimension of BN minus 1
C               
C   Accessed common blocks:
C
C   Local storage:
C
C      (MAXL+1)*((MAXLS+1)**3+(MAXL+1)*(2*MAXL+1)**2) DOUBLE PRECISION
C      words. With MAXLS=30, MAXL=3 and 8-byte DP word, it's some 941Kb.
C      This code wasn't intended for small memories anyway, wasn't it?
C
C   Precision:
C
C      Coefficients are good to almost machine precision with low L and
C      LS. At L=3 and LS=20, relative loos of precision is about 11 digits,
C      so that this routine should preferably be used with real*16 
C      temporaries.
C
C   Module logic:
C
C      This code is *not* intended as a direct replacement of the PS3SHD
C      precomputed coefficents. Low-order coefficients should be taken
C      from BLOCK DATA's, where they are computed to machine precision.
C      High-order coefficents, where some loss of significance is acceptable,
C      might be computed here.
C
C   Possible optimizations:
C
C   Bugs:
C
C      Maximum value of LS is limited by MAXLS.
C      Maximum value of L is limits by MAXL.
C      The REAL*16 temporaries are necessary to get meaningful results with
C      large L and LS.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXLS=40)
      PARAMETER (MAXL =3)
C    Comment next line if your compiler has trouble with REAL*16
      REAL*16 ZERO, TWO, HALF, BETALS, BETAL, SUM, SUMP, SUMQ
C
      PARAMETER (ZERO =0)
      PARAMETER (DZERO=0.D0)
      PARAMETER (TWO  =2)
      PARAMETER (HALF =1/TWO)
C
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION BETALS(0:MAXLS,0:MAXLS,0:MAXLS,0:MAXL)
      DIMENSION BETAL (0:2*MAXL,0:2*MAXL,0:MAXL,0:MAXL)
C
      DIMENSION BN(0:LDS,0:LDN,0:L)
C    Parameter verification
      IF( L.LT.0 .OR. L.GT.MAXL ) THEN
          WRITE(NB6,10010) L, MAXL
          STOP 'PS3SBN'
      ENDIF
      IF( LS.LT.0 .OR. LS.GT.MAXLS ) THEN
          WRITE(NB6,10020) LS, MAXLS
          STOP 'PS3SBN'
      ENDIF
      IF( LDS.LT.L+LS ) THEN
          WRITE(NB6,10030) LDS, L+LS
          STOP 'PS3SBN'
      ENDIF
      IF( LDN.LT.L+LS ) THEN
          WRITE(NB6,10040) LDN, L+LS
          STOP 'PS3SBN'
      ENDIF
C    Summation limits
      MAXM = MIN(L,LS)
      MAXO = L + LS
C    Auxiliary beta values
      CALL PS3SBT(LS,MAXM,BETALS,MAXLS,MAXLS,MAXLS)
      CALL PS3SBT(L,-MAXM,BETAL,2*MAXL,2*MAXL,MAXL)
C    Summation itself, a faithful reproduction of Sharma's 
C    expression for Bnu.
      DO 3000 M=0,MAXM
          DO 2900 IS=0,MAXO
              DO 2800 NU=0,MAXO-IS
                  IPPMAX = MIN(LS,NU)
                  SUM    = ZERO
                  DO 2700 IPP=0,IPPMAX
                      IQPMIN = MAX(0,IS-L-M)
                      DO 2600 IQP=IQPMIN,IS
                          SUMP  = ZERO
                          IQMIN = MAX(0,IPP+IQP+NU-LS+M)
                          IQMAX = MIN(NU-IPP,L+M-IS+IQP)
                          DO 2500 IQ=IQMIN,IQMAX
                              SUMQ  = ZERO
                              IPMAX = MIN(IS-IQP,L+M-IS+IQP-IQ)
                              DO 2400 IP=0,IPMAX
                                  SUMQ = SUMQ + BETAL(IS-IP-IQP,IQ,IP,M)
 2400                         CONTINUE
                              IF( MOD(IQ,2).NE.0 ) SUMQ = -SUMQ
                              SUMP=SUMP+BETALS(IQP,NU-IQ-IPP,IPP,M)*SUMQ
 2500                     CONTINUE
                          IF( MOD(IQP,2).NE.0 ) SUMP = -SUMP
                          SUM = SUM + SUMP
 2600                 CONTINUE
 2700             CONTINUE
                  IF( MOD(L,2).NE.0 ) SUM = -SUM
C                Truncation from REAL*16 to DOUBLE PRECISION is intentional
C                here: after this point, there is no possibility of numerical
C                instabilities arising.
                  BN(NU,IS,M) = DBLE(HALF*SUM)
 2800         CONTINUE
 2900     CONTINUE
 3000 CONTINUE
C    Zero out dummies
      DO 4000 M=MAXM+1,L
          DO 3900 IS=0,MAXO
              DO 3800 NU=0,MAXO-IS
                  BN(NU,IS,M) = DZERO
 3800         CONTINUE
 3900     CONTINUE
 4000 CONTINUE
C
      RETURN
10010 FORMAT(' L = ',I5,' IS ILLEGAL IN PS3SBN. VALID RANGE IS 0 TO ',
     .       I5)
10020 FORMAT(' LS = ',I5,' IS ILLEGAL IN PS3SBN. VALID RANGE IS 0 TO ',
     .       I5)
10030 FORMAT(' BN() WILL OVERFLOW ON THE FIRST DIMENSION IN PS3SBN.'/
     .       ' HAVE = ',I5,' WANT = ',I5)
10040 FORMAT(' BN() WILL OVERFLOW ON THE SECOND DIMENSION IN PS3SBN.'/
     .       ' HAVE = ',I5,' WANT = ',I5)
      END
