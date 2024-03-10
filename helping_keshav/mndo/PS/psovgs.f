C     ******************************************************************
C
C     NMR 3-center spin-orbit integrals over STO-nG functions.
C
C     ******************************************************************
C
C     This is an extremely naive implementation. Any proper current
C     integrals package should be able to do much better, but I wanted
C     to see how difficult it is to calculate  these integrals in a
C     straightforward manner.
C
C     User-callable top level routines defined here are:
C
C      PSGLR3 - Evaluates three-center one-electron r^-3 L integrals
C               over Slater functions by Gaussian expansion.
C      PSGRL3 - Evaluates three-center one-electron r^a r^-3 L integrals
C               over Slater functions by Gaussian expansion.
C      PSGRR3 - Evaluates three-center one-electron {x,y,z} r^-3 r^c
C               integrals over Slater functions by Gaussian expansion.
C
C   All integrals are computed over elementary nests, containing orbitals
C   in the form:
C
C                      Nx   Ny   Nz             -Alpha  r**2
C      /  / Pmax         p    p    p \  Nrad          q       \
C      |  |  Sum   C  X    Y    Z    |   Sum   E              |
C      \  \  p=1    p                /   q=1                  /
C
C   All orbitals in the nest share the same exponential part, but may
C   have different radial parts. Nests are characterized by the following
C   variables:
C
C      NO                   = Number of orbitals in the nest
C      NR                   = Number of radial components
C      NA                   = Maximum order of angular component
C      CS(IX,IY,IZ,IR,IORB) = C coefficients
C      ALP(IR)              = Alpha coefficients
C      MAXA                 = Top indices for the first three dimensions of CS
C      MAXR                 = Top index for the fourth dimension of CS.
C
C   The following elementary operation over orbital nests are defined by the
C   PSOVxy functions:
C
C      MK  = Creates a nest corresponding to |NL?> block of Slater functions
C      ZR  = Zeros out angular components
C      RS  = Applies r operator to the angular part of the nest
C      RV  = Computes vector product of the operator r and three-component
C            block of functions.
C      DX  = Applies gradient operator to the nest
C      TR  = Translates a nest
C      PR  = Computes product of the angular parts of two nests
C      DM  = Prints a verbose description of the nest
C      GP  = Computes position and exponent of a product of two
C            elementary gaussians.
C
      SUBROUTINE PSOGJ(MAXO,ALPHA,GAMMA,AINTJ,LDA)
C
C   Computes a set of elementary J integrals using one-dimentional
C   numerical integration.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      MAXO   - Maximum value of O parameter (minimum is 1).
C      ALPHA  - Value of the alpha parameter. Should be strictly positive.
C      GAMMA  - Value of the gamma parameter. May be positive, negative,
C               or zero.
C      AINTJ  - (Output) Two-dinetional array, place for integrals
C               First dimension: p parameter value, 0 to o/2
C               Second dimension: o parameter value, 1 to MAXO
C      LDA    - Leading dimension of the AINTJ array.
C
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C      Basically, (2*IOMAX+3)*IORDMX DP words. With IOMAX=10 and IORDMX=120,
C      it is about 2800 words, or 22K.
C
C   Module logic:
C
C      J integral is given by:
C                                   2
C                             -gamma   inf                2
C                             -------   /   o-1  - alpha r
C        J   (alpha,gamma) =  4 alpha   |  r    E           beta    (gamma r) d r
C         o,p                E          /                       o-2p
C                                       0
C
C      beta is second exponential integral, given by eq. 5.1.6 of Abramowitz
C      & Stegun and implemented by PSVXBG in Slater integrals package.
C
C      Analytical expression for this integral (in terms of error function
C      Erf) is known, but is fairly complicated and numerically unstable for
C      larger values of o and p parameters.
C
C      Numerical integration is performed with Laguerre formula on the interval
C      r = Abs[gamma]/(2*alpha) to infinity, using variable substitution:
C
C             Abs[gamma] + Sqrt[ gamma^2 + 4*alpha*t ]
C        r = ------------------------------------------
C                          2*alpha
C
C      Integration on the interval 0 to Abs[gamma]/(2*alpha) is performed by
C      Gaussian formula applied to five intervals separated by points:
C
C              Sqrt[ 2*n+1 +/- Sqrt[ 8*n+1 ] ]
C        Ai = ---------------------------------, where n is selected as
C                        2 Sqrt[alpha]
C
C        n = MAXO/Sqrt[2]
C
C      Some of the three gaussian partitions may be not present.
C
C   Precision:
C
C      For |gamma|/alpha ratio not exceeding 20, results are good to machine
C      precision. For ratios up to 100, accuracy may detiriorate to about
C      7 digits. For higher ratious, evaluation of beta function will
C      fail for large argument values.
C
C   Possible optimizations:
C
C      Better selection of integration points should help.
C      Most of the time (some 80% on R4000) is spent evaluating beta functions and
C      exponents.
C
C   Bugs:
C
C      Gaussian integration is clearly suboptimal for this integrand, especially
C      for large gamma/alpha ratious.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IORDL=10,IORDG=10,IOMAX=10)
      PARAMETER (MAXGP=7)
      PARAMETER (IORDMX=(MAXGP+4)*IORDG+IORDL)
      PARAMETER (HUGEV=300.D0)
C
      PARAMETER (ZERO =0.D0)
      PARAMETER (RSQ2 =0.7071067811865475244008443621048490392848D0)
      PARAMETER (HALF =0.5D0)
      PARAMETER (ONE  =1.D0)
      PARAMETER (TWO  =2.D0)
      PARAMETER (FOUR =4.D0)
      PARAMETER (EIGHT=8.D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION ARGS(IORDMX), WGTS(IORDMX)
      DIMENSION BETAS(0:IOMAX,IORDMX), EXPVAL(1:IOMAX,IORDMX)
      DIMENSION NPTS(MAXGP+1), I1ST(MAXGP+1), BORDER(MAXGP+1)
C
      DIMENSION AINTJ(0:LDA-1,MAXO)
C    Range checks
      IF( MAXO.LT.1 .OR. MAXO.GT.IOMAX ) THEN
          WRITE(NB6,11010) MAXO
          STOP 'PSOGJ'
      ENDIF
      IF( ALPHA.LE.ZERO ) THEN
          WRITE(NB6,11020) ALPHA
          STOP 'PSOGJ'
      ENDIF
      IF( LDA-1.LT.MAXO/2 ) THEN
          WRITE(NB6,11030) LDA
          STOP 'PSOGJ'
      ENDIF
C
      AGAMMA = ABS(GAMMA)
      GAMMA2 = GAMMA**2
C    Select integration partitions
      BLGR1  = HALF*AGAMMA/ALPHA
      BLGR2  = AGAMMA/ALPHA
      DN     = MAXO*RSQ2
      TMP1   = TWO*DN+ONE
      TMP2   = SQRT(EIGHT*DN+1)
      BGLOW  = ZERO
      IF( TMP2.LT.TMP1 ) BGLOW = SQRT( (TMP1-TMP2)/(FOUR*ALPHA) )
      BGHIGH = SQRT( (TMP1+TMP2)/(FOUR*ALPHA) )
      BGPEAK = BGHIGH - BGLOW
C
C    Potential partitions are:
C     A  0              to BGLOW-BGPEAK/2      Gauss
C     B  BGLOW-BGPEAK/2 to BGLOW               Gauss
C     C  BGLOW          to BGHIGH              Gauss
C     D  BGHIGH         to BGHIGH + BGPEAK     Gauss
C     E  BGHIGH+BGPEAK  to BLGR1/2             Gauss
C     F  BLGR1/2        to BLGR1 + BLGR1/2     Gauss
C     G  BLGR1/2        to BLGR2               Gauss
C     H  BLGR           to Infinity            Laguerre
C
      BORDER(1) = ZERO
      NPTS  (1) = IORDG
      BORDER(2) = BGLOW-BGPEAK*HALF
      NPTS  (2) = IORDG
      BORDER(3) = BGLOW
      NPTS  (3) = IORDG*2
      BORDER(4) = BGHIGH
      NPTS  (4) = IORDG
      BORDER(5) = BGHIGH+BGPEAK
      NPTS  (5) = IORDG
      BORDER(6) = BLGR1*HALF
      NPTS  (6) = IORDG*4
      BORDER(7) = BLGR1 + HALF*BLGR1
      NPTS  (7) = IORDG
      BORDER(8) = BLGR2
*     write(*,*)' borders: '
*     write(*,'((1X,F20.12))') border
C    Drop non-existent partitions
      NPART = 0
      DO 100 I=1,MAXGP
          IF( BORDER(I+1).GT.BORDER(NPART+1) ) THEN
              NPART           = NPART + 1
              NPTS(NPART)     = NPTS(I)
              BORDER(NPART+1) = BORDER(I+1)
          ENDIF
  100 CONTINUE
C    Assign integration orders and starting indices
      I1ST(1) = 1
      DO 200 I=1,NPART
          I1ST(I+1) = I1ST(I) + NPTS(I)
  200 CONTINUE
C    Get integration abscissas and weights for Gaussian integration
      DO 300 I=1,NPART
          CALL PSIGA(NPTS(I),NPTS(I),BORDER(I),BORDER(I+1),
     .               ARGS(I1ST(I)),WGTS(I1ST(I)))
  300 CONTINUE
C    Get integration abscissas and weights for Laguerre integration
C    Abscissas are for the temporary variable T at this point
      TLOW = ALPHA*BORDER(NPART+1)**2 - AGAMMA*BORDER(NPART+1)
      CALL PSILGA(IORDL,IORDL,ONE,TLOW,
     .            ARGS(I1ST(NPART+1)),WGTS(I1ST(NPART+1)))
C    Compute values of the R at the Laguerre integration points and
C    incorporate dr/dt into integration weights.
      R2ALP = ONE/(TWO*ALPHA)
      DO 400 I=I1ST(NPART+1),I1ST(NPART+1)+IORDL-1
          T = ARGS(I)
          TEMP = GAMMA2 + FOUR*ALPHA*T
          IF( TEMP.LT.ZERO ) TEMP = ZERO
          R = R2ALP*( AGAMMA + SQRT(TEMP) )
          ARGS(I) = R
          WGTS(I) = WGTS(I) / ( TWO*ALPHA*R - AGAMMA )
  400 CONTINUE
      NPART       = NPART + 1
      NPTS(NPART) = IORDL
      NPOINT      = I1ST(NPART) + IORDL - 1
      IF( NPOINT.GT.IORDMX ) THEN
          WRITE(NB6,11090) NPOINT, IORDMX
          STOP 'PSOGJ'
      ENDIF
C    From this point on, all integration points can be treated alike.
*     write(*,*)' npoint = ', npoint
*     write(*,'((1x,i3,2(2x,f20.12)))')(i,args(i),wgts(i),i=1,npoint)
C    Very large GAMMA*R arguments can kill PSVXBG, so we want
C    to drop these points altogether. We can only hope these points
C    are in the zero already...
      DO 500 IPT=1,NPOINT
          IF( MAX(ABS(GAMMA*ARGS(IPT)),ALPHA*ARGS(IPT)**2).GT.HUGEV ) 
     .                                                             THEN
              NPOINT = IPT-1
              GOTO 501
          ENDIF
  500 CONTINUE
  501 CONTINUE
C    Compute prefactors
      SCALE = EXP(-GAMMA2/(FOUR*ALPHA))
      DO 1000 IPT=1,NPOINT
          R = ARGS(IPT)
C        EXPVAL contains r^(o-1) Exp[-alpha r^2] values
C        multiplied by 2.0 (so as to compensate for beta values)
C        and by the integration weights.
C        There is no danger of division by zero, since Gaussian
C        integration points are always inside the interval.
          TEMP = TWO*EXP(-ALPHA*R**2)*WGTS(IPT)
          EXPVAL(1,IPT) = TEMP
          DO 800 IO=2,MAXO
              TEMP = TEMP*R
              EXPVAL(IO,IPT) = TEMP
  800     CONTINUE
C        BETAS contains values of the exponential integral beta at 
C        integration points multiplied by 0.5.
          CALL PSVXBG(GAMMA*R,MAXO,BETAS(0,IPT))
          DO 900 I=0,MAXO
              BETAS(I,IPT) = SCALE*BETAS(I,IPT)
  900     CONTINUE
 1000 CONTINUE
C    Accumulate integrals. This is not the best possible arrangement
C    of summation indices from cache line reuse point of view. However,
C    temporaries here should be small enough to fit in caches, so it
C    probably is not important.
      DO 6000 IO=1,MAXO
          DO 5000 IP=0,IO/2
              TEMP = ZERO
              DO 4000 IPT=1,NPOINT
                  TEMP = TEMP + EXPVAL(IO,IPT)*BETAS(IO-2*IP,IPT)
*         write(*,'(1X,I3,1X,I3,1X,I3,2X,F20.15,2X,F35.15)')
*    .    io, ip, ipt, args(ipt), EXPVAL(IO,IPT)*BETAS(IO-2*IP,IPT)
 4000          CONTINUE
               AINTJ(IP,IO) = TEMP
 5000     CONTINUE
 6000 CONTINUE
C    Here we are.
      RETURN
11010 FORMAT(' MAXO = ',I5,' IS INVALID IN PSOGJ.')
11020 FORMAT(' ALPHA = ',G15.8,' IS NOT POSITIVE IN PSOGJ.')
11030 FORMAT(' LDA = ',I5,' IS TOO SMALL IN PSOGJ.')
11090 FORMAT(' CATASATROPHE IN PSOGJ. WANT ',I4,' PTS, CAN HANDLE ONLY '
     .       ,I4)
      END
C
      SUBROUTINE PSOGIP(MAXO,ALPHA,A,AINTI,LDA)
C
C   Computes a set of I' integrals (local coordinate system).
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      MAXO   - Maximum value of I+J+K. Minimum is 1.
C      ALPHA  - Value of the orbital exponent. Should be strictly positive.
C      A      - Displacement of the operator center on Z axis. May be either 
C               positive or negative.
C      AINTI  - (Output) Three-dinetional array, place for integrals
C               First dimension: i parameter value, from 0 to MAXO
C               Second dimension: j parameter value, from 0 to MAXO
C               Third dimension: k parameter value, 0 to MAXO
C               Additionally, 1 <= i+j+k <= maxo. Only even values of
C               i and j result in non-zero integrals.
C      LDA    - First two leading dimensions of the AINTI array.
C
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C      About 0.25*(IOMAX+1)*(IOMAX+5) DP words. With IOMAX=10, it's only
C      40 or so words.
C
C   Module logic:
C
C      I' integral is given by:
C                                     2                
C                             -alpha a  (i-1)!! (j-1)!!
C      I'   (alpha,a) = 2 Pi E          ---------------
C        ijk                                (i+j)!!    
C
C                                                           i+j
C               Inf                       2  1              ---
C                /       i+j+k-1  -alpha r   /      k     2  2   2 alpha r t
C                | d r  r        E           | d t t  (1-t )    E
C                /                           /
C                0                          -1
C
C      It is evaluated as:
C                                            i+j             i+j
C                                            ---   i+j       --- - p
C                            (i-1)!! (j-1)!!  2  / --- \      2
C      I'   (alpha,a) = 2 Pi --------------- Sum |  2  | (-1)      J       (alpha,-2 alpha a)
C        ijk                     (i+j)!!     p=0 \  p  /            i+j+k,p
C
C      Integrals J are evaluated by PSOGJ.
C
C   Precision:
C
C      Is determined by J integrals. For moderate ALPHA (<5) and A (>0.1 <5)
C      values and small MAXS, integrals are good to machine precision.
C
C   Possible optimizations:
C
C      Probably are not worth the time. Evaluation of J integrals is much more
C      expensive than anything we do here.
C
C   Bugs:
C
C      It is not really necessary to have separate X and Y components of this
C      integral.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IOMAX=10)
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (STWO=-2.D0)
      PARAMETER (TWOPI=6.2831853071795864769252867665590057683943388D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION AINTJ(0:IOMAX/2,IOMAX)
      DIMENSION FACT(0:IOMAX), DFACT(-1:IOMAX), CP(0:IOMAX/2)
C
      DIMENSION AINTI(0:LDA-1,0:LDA-1,0:MAXO)
C    Range checks
      IF( MAXO.GT.IOMAX ) THEN
          WRITE(NB6,11010) MAXO, IOMAX
          STOP 'PSOGIP'
      ENDIF
      IF( LDA.LT.MAXO+1 ) THEN
          WRITE(NB6,11020) LDA, MAXO+1
          STOP 'PSOGIP'
      ENDIF
C    Compute J integrals
      CALL PSOGJ(MAXO,ALPHA,STWO*A*ALPHA,AINTJ,IOMAX/2+1)
C    Compute factorials and double factorials for this run
      FACT(0) = ONE
      DO 100 I=1,MAXO/2
          FACT(I) = FACT(I-1)*I
  100 CONTINUE
      DFACT(0) = ONE
      DO 200 I=2,MAXO,2
          DFACT(I) = DFACT(I-2)*I
  200 CONTINUE
      DFACT(-1) = ONE
      DO 300 I=1,MAXO,2
          DFACT(I) = DFACT(I-2)*I
  300 CONTINUE
C    This loop nest is in the wrong order for cache reuse in AINTI.
C    It shouldn't matter, since everything is going to fit in anyway.
      DO 4000 I=0,MAXO,2
          DO 3000 J=0,MAXO-I,2
              IJ2 = (I+J)/2
              SCALE = TWOPI*DFACT(I-1)*DFACT(J-1)*FACT(IJ2)/DFACT(I+J)
C            Compute distribution coefficients for (i+j)/2, multiplied
C            by overall scaling for this IJ pair.
              DO 500 IP=0,IJ2
                  CP(IP) = SCALE/( FACT(IP)*FACT(IJ2-IP) )
                  IF(MOD(IJ2-IP,2).NE.0) CP(IP) = -CP(IP)
  500         CONTINUE
C            Evaluate integrals
              DO 2000 K=0,MAXO-I-J
                  IJK = I+J+K
                  IF( IJK.LE.0 ) GOTO 2000
                  SUM = ZERO
                  DO 1000 IP=0,IJ2
                      SUM = SUM + CP(IP)*AINTJ(IP,IJK)
 1000             CONTINUE
                  AINTI(I,J,K) = SUM
 2000         CONTINUE
 3000     CONTINUE
 4000 CONTINUE
C    Clear zero and undefined integrals
      AINTI(0,0,0) = ZERO
      DO 5000 K=0,MAXO
          DO 4700 J=1,MAXO-K,2
              DO 4600 I=0,MAXO-(J+K)
                  AINTI(I,J,K) = ZERO
 4600         CONTINUE
 4700     CONTINUE
          DO 4900 J=0,MAXO-K,2
              DO 4800 I=1,MAXO-(J+K),2
                  AINTI(I,J,K) = ZERO
 4800         CONTINUE
 4900     CONTINUE
 5000 CONTINUE
C    Done!
      RETURN
11010 FORMAT(' MAXO (=',I5,') IS TOO BIG IN PSOGIP. LIMIT IS ',I5)
11020 FORMAT(' LDA (=',I5,') IS TOO SMALL IN PSOGIP. NEED AT LEAST ',I5)
      END
C
      SUBROUTINE PSOGIR(MAXO,A,B,ROT,LDR)
C
C   Computes 2D rotation coefficients for X**I Y**J factors.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      MAXO   - Maximum value of I+J.
C      A,B    - Transformation coefficients, see below.
C      ROT    - (Output) three-dimentional array of rotation coefficients.
C               First index, 0 to MAXO: IP index. Only even IP coefficients
C                   are actually computed.
C               Second index, 0 to MAXO: I index.
C               Third index, 0 to MAXO: J index.
C               Only indices satisfying inequalities: I+J <= MAXO and
C               IP <= I+J are computed.
C      LDR    - First two leading dimensions of the ROT matrix.
C
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C      About (IOMAX+1)*(IOMAX+5) DP words. With IOMAX=10, it's about 
C      170 words, or 1.3Kb.
C
C   Module logic:
C
C      If X and Y are assumed to transform according to:
C
C         X = A*U + B*V
C         Y =-B*U + A*V
C
C      X**I Y**J transforms according to:
C
C          I  J        I,J        I'  I+J-I'
C         X  Y  = Sum C          U   V
C                  I'  I',I+J-I'
C
C      where coefficients C are given by:
C
C          I,J             / I \ /  J \     I'-p  J-I'+2p  I+I'-2p
C         C          = Sum |   | |    | (-1)     A        B
C          I',I+J-I'    p  \ p / \I'-p/
C
C   Possible optimizations:
C
C      In pronciple, we could use A**2 + B**2 = 1 identity to simplify
C      expressions for C coefficients further.
C
C      Most of the time is spent on evaluation of the addresses inside
C      inner loop. Smart compiler might definitly help here...
C
C      In the first rotation, only even I's are necessary, so having
C      a separate routine in this case can shave a few more clocks.
C
C   Bugs:
C
C      As is, code works only for even I' (IP) values. If a general 
C      version is desired, enable code lines marked with CANYIP
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IOMAX=10)
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C    Local storage
      DIMENSION FACT(0:IOMAX), RFACT(0:IOMAX), BINOM(0:IOMAX,0:IOMAX)
      DIMENSION AN(0:IOMAX), BN(0:IOMAX)
C    Parameters
      DIMENSION ROT(0:LDR-1,0:LDR-1,0:IOMAX)
C    Range checks
      IF( MAXO.GT.IOMAX ) THEN
          WRITE(NB6,11010) MAXO, IOMAX
          STOP 'PSOGIR'
      ENDIF
      IF( LDR.LT.MAXO+1 ) THEN
          WRITE(NB6,11020) LDR, MAXO+1
      ENDIF
C    Compute factorials and binomials
      FACT (0)   = ONE
      RFACT(0)   = ONE
      BINOM(0,0) = ONE
      DO 200 I=1,MAXO
          FACT (I) = FACT(I-1) * I
          RFACT(I) = ONE/FACT(I)
          DO 100 J=0,I
              BINOM(J,I) = FACT(I)*RFACT(J)*RFACT(I-J)
  100     CONTINUE
  200 CONTINUE
C    Compute powers of A and B parameters
      AN(0) = ONE
      BN(0) = ONE
      DO 300 I=1,MAXO
          AN(I) = AN(I-1)*A
          BN(I) = BN(I-1)*B
  300 CONTINUE
C    Compute transformation coefficients
      DO 4000 J=0,MAXO
          DO 3000 I=0,MAXO-J
CANYIP        DO 2000 IP=0,I+J
              DO 2000 IP=0,I+J,2
                  SUM    = ZERO
                  ITLOW  = MAX(0,IP-J)
                  ITHIGH = MIN(I,IP)
                  IF( ITHIGH.GE.ITLOW ) THEN
                      JIP = J-IP
                      IIP = I+IP
                      DO 1000 IT0=ITLOW,ITHIGH-1,2
                          IT1 = IT0 + 1
                          SUM = SUM + BINOM(IT0,I)*BINOM(IP-IT0,J)*
     .                                    AN(JIP+2*IT0)*BN(IIP-2*IT0)
     .                              - BINOM(IT1,I)*BINOM(IP-IT1,J)*
     .                                    AN(JIP+2*IT1)*BN(IIP-2*IT1)
 1000                 CONTINUE
                      IF( MOD(ITHIGH-ITLOW,2).EQ.0 ) THEN
                          IT0 = ITHIGH
                          SUM = SUM + BINOM(IT0,I)*BINOM(IP-IT0,J)*
     .                                    AN(JIP+2*IT0)*BN(IIP-2*IT0)
                      ENDIF
CANYIP                IF( MOD(ITLOW+IP,2).NE.0 ) SUM = -SUM
                      IF( MOD(ITLOW,2).NE.0 ) SUM = -SUM
                  ENDIF
                  ROT(IP,I,J) = SUM
 2000         CONTINUE
 3000     CONTINUE
 4000 CONTINUE
C
      RETURN
11010 FORMAT(' MAXO (=',I5,') IS TOO BIG IN PSOGIR. LIMIT IS ',I5)
11020 FORMAT(' LDA (=',I5,') IS TOO SMALL IN PSOGIR. WANT AT LEAST ',I5)
      END
C
      SUBROUTINE PSOGI(MAXO,ALPHA,R,AINTI,LDA,ISZERO)
C
C   Computes a set of I integrals (molecular coordinate system).
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      MAXO   - Maximum value of I+J+K. Minimum is 1.
C      ALPHA  - Value of the orbital exponent. Should be strictly positive.
C      R      - Coordinates of the exponent center.
C      AINTI  - (Output) Three-dimetional array, place for integrals
C               First dimension: i parameter value, 0 to MAXO
C               Second dimension: j parameter value, 0 to MAXO
C               Third dimension: k parameter value, 0 to MAXO
C               Additionally, 1 <= i+j+k <= maxo
C      LDA    - First two leading dimensions of the AINTI array.
C      ISZERO - (Output) Set to .TRUE. if all elementary integrals
C               are smaller than tolerance.
C
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C      About 3*(IOMAX+1)**3 DP words. With IOMAX=10, it is about 4000
C      DP words, or some 32K.
C
C   Module logic:
C
C      I is given by real space integral:
C
C                  i j k           ->  -> 2
C               / x y z   - alpha (r - R )   ->
C         I   = | ------ E                 d r
C          ijk  /    3
C                   r
C
C      We evaluate it in local coordinate system pointing towards R,
C      and rotate into molecular coordinate system. Two-stage rotation
C      is used. First rotation (going molecular->local) is counterclock-
C      wise on OZ axis by ALPHA, bringing vector R into OZ'Y' plane. 
C      Second rotation is counterclockwise along OX' axis, making OZ''
C      axis to coincide with vector R.
C
C   Precision:
C
C      Is determined by precision of I' integrals in local coordinate
C      system.
C
C   Possible optimizations:
C
C      Going for spherical harmonics instead of cartesian components
C      should improve staling for large MAXO dramatically.
C
C      For small MAXO values (MAXO<=6, which covers everything we are
C      going to need for MNDO/d shieldings), execution time is dominated
C      by evaluation of beta integrals (PSVXBG) and exponentiations.
C
C      For larger MAXO values (MAXO<=10), evaluation of beta integrals
C      is about as expensive as rotation coefficients for cartesian
C      products (PSOGIR), with PSOGI and PSOGJ coming next.
C
C      For still larger MAXO values, execution time should be dominated
C      by PSOGIR, which scales as O(MAXO**4) and contains fairly 
C      complicated address experessions.
C
C   Bugs:
C
C      This routine scales as O(MAXO**4), and is impractical for large MAXO
C      values.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IOMAX=10)
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (EPS=1.D-18)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C    Temporaries
      DIMENSION AP(0:IOMAX,0:IOMAX,0:IOMAX)
      DIMENSION AX(0:IOMAX,0:IOMAX,0:IOMAX)
      DIMENSION ROT(0:IOMAX,0:IOMAX,0:IOMAX)
C    Parameters
      DIMENSION R(3), AINTI(0:LDA-1,0:LDA-1,0:MAXO)
      LOGICAL ISZERO
C    Range checks
      IF( MAXO.GT.IOMAX ) THEN
          WRITE(NB6,11010) MAXO, IOMAX
          STOP 'PSOGI'
      ENDIF
      IF( LDA.LT.MAXO+1 ) THEN
          WRITE(NB6,11020) LDA, MAXO+1
          STOP 'PSOGI'
      ENDIF
C    Evaluate distance and Euler angles for the transformation
C    into local coordinate system.
      A      = SQRT(R(1)**2+R(2)**2+R(3)**2)
      RXY    = SQRT(R(1)**2+R(2)**2)
      IF( A.GT.ZERO ) THEN
          SINALP = RXY/A
          COSALP = R(3)/A
      ELSE
          SINALP = ZERO
          COSALP = ONE
      ENDIF
      IF( RXY.GT.ZERO ) THEN
          SINBET = R(1)/RXY
          COSBET = R(2)/RXY
      ELSE
          SINBET = ZERO
          COSBET = ONE
      ENDIF
C    Evaluate integrals in the local coordinate system
      CALL PSOGIP(MAXO,ALPHA,A,AP,IOMAX+1)
C    Align local OZ axis with molecular OZ axis.
C    First, we'd need transformation matrix for Y**J*Z**K products
      CALL PSOGIR(MAXO,COSALP,SINALP,ROT,IOMAX+1)
C    We know local system integrals to be zero for odd powers
C    of X and Y, so we don't need to process all of them now.
C    (0,0,0) integral is technically undefined, but it won't do
C    any harm to evaluate it anyway.
      DO 1000 I=0,MAXO,2
          DO 900 K=0,MAXO-I
              DO 800 J=0,MAXO-(I+K)
                  JK  = J+K
                  SUM = ZERO
                  DO 700 JP=0,JK,2
                      SUM = SUM + ROT(JP,J,K)*AP(I,JP,JK-JP)
  700             CONTINUE
                  AX(I,J,K) = SUM
  800         CONTINUE
  900     CONTINUE
 1000 CONTINUE
C    Now, align local OX axis with molecular OX axis.
C    All integrals with odd powers of X in AX are still zero, which we can use.
      CALL PSOGIR(MAXO,COSBET,SINBET,ROT,IOMAX+1)
      ISZERO = .TRUE.
      DO 2000 K=0,MAXO
          DO 1900 J=0,MAXO-K
              DO 1800 I=0,MAXO-(J+K)
                  IJ  = I+J
                  SUM = ZERO
                  DO 1700 IP=0,IJ,2
                      SUM = SUM + ROT(IP,I,J)*AX(IP,IJ-IP,K)
 1700             CONTINUE
                  AINTI(I,J,K) = SUM
                  IF( ABS(SUM).GT.EPS ) ISZERO = .FALSE.
 1800         CONTINUE
 1900     CONTINUE
 2000 CONTINUE
C    Got it, dude?
      RETURN
11010 FORMAT(' MAXO (=',I5,') IS TOO BIG IN PSOGI. LIMIT IS ',I5)
11020 FORMAT(' LDA (=',I5,') IS TOO SMALL IN PSOGI. WANT AT LEAST ',I5)
      END
C
      SUBROUTINE PSOGMK(IGO,N,L,ZETA,NO,NR,NA,CS,ALP,MAXA,MAXR)
C
C   Generates STO-nG approximation for a shell of Slater functions.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IGO    - "n" in STO-nG.
C      N      - Major quantum number of Slater orbitals.
C      L      - Orbital quantum number of Slater orbitals.
C      ZETA   - Orbital exponent of Slater orbitals.
C      NO,NR,CS,ALP
C             - (Output) as described on top.
C      MAXA,MAXR
C             - Dimensioning parameters.
C
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C   Module logic:
C
C      This is a no-brainer function. General expressions *are*
C      available, at least for the angular part. However, they
C      are so complicated you wouldn't want to see 'em.
C
C      Orbitals in the output nest are in the order of increasing
C      "magnetic" quantum number.
C
C   Possible optimizations:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXGO=9,MAXPOW=5)
C
      PARAMETER (ZERO=0.D0)
C    Weight coefficients in the angular part
      PARAMETER (C1 =0.282094791773878143474039725780D0)
      PARAMETER (C2 =0.488602511902919921586384622838D0)
      PARAMETER (C3 =1.092548430592079070543385705803D0)
      PARAMETER (C4 =0.315391565252520006030893690296D0)
      PARAMETER (C5 =0.630783130505040012061787380591D0)
      PARAMETER (C6 =0.546274215296039535271692852901D0)
      PARAMETER (C7 =1.77013076977993053103683083262 D0)
      PARAMETER (C8 =0.590043589926643510345610277541D0)
      PARAMETER (C9 =2.89061144264055405538938015440 D0)
      PARAMETER (C10=0.457045799464465736158020696917D0)
      PARAMETER (C11=1.82818319785786294463208278767 D0)
      PARAMETER (C12=1.119528997770346174243187739597D0)
      PARAMETER (C13=0.746352665180230782828791826398D0)
      PARAMETER (C14=1.445305721320277027694690077199D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      EXTERNAL PSSNRM
      DIMENSION CS(0:MAXA,0:MAXA,0:MAXA,MAXR,*), ALP(MAXR)
C    Local temporary storage
      DIMENSION DQ(MAXGO)
C    (unnormalized) Gaussian expansions for (unnormalized) radial parts 
C    of Slater orbitals.
      DIMENSION GS(2,MAXGO,MAXGO,MAXPOW)
C    STO-1G 4S and 5S functions are *really* a joke in this case
C    1S, EPS =   0.0106812112657
      DATA (GS(I,1,1,1),I=1,2)/0.2709498090630D0, 0.464162509257D0/
C    2S, EPS =   0.0047118381051
      DATA (GS(I,1,1,2),I=1,2)/0.1012151084320D0, 0.391391278820D0/
C    3S, EPS =   0.0939041334213
      DATA (GS(I,1,1,3),I=1,2)/0.0529688175686D0, 0.656048317184D0/
C    4S, EPS =   3.25423759743
      DATA (GS(I,1,1,4),I=1,2)/0.0326460027360D0, 1.685966035680D0/
C    5S, EPS = 124.111496729
      DATA (GS(I,1,1,5),I=1,2)/0.0221691293830D0, 5.892152008050D0/
C    STO-2G
C    1S, EPS =   0.000789524224233
      DATA (GS(I,1,2,1),I=1,2)/0.8518186635030D0, 0.481012932267D0/
      DATA (GS(I,2,2,1),I=1,2)/0.1516232926810D0, 0.208059273742D0/
C    2S, EPS =   0.00364128056619 
*     DATA (GS(I,1,2,2),I=1,2)/5.0306586068800D0,-0.285757794358D0/
*     DATA (GS(I,2,2,2),I=1,2)/0.1024780007890D0, 0.397317946389D0/
C    2S, EPS =   0.0020032820933  
      DATA (GS(I,1,2,2),I=1,2)/0.1292278610780D0, 0.351846309146D0/
      DATA (GS(I,2,2,2),I=1,2)/0.0490858420488D0, 0.0650784650196D0/
C    3S, EPS =   0.00429800442761 
      DATA (GS(I,1,2,3),I=1,2)/0.6694095821980D0,-0.678063086956D0/
      DATA (GS(I,2,2,3),I=1,2)/0.0583713509450D0, 0.747852137075D0/
C    4S, EPS =   0.0137652491836
      DATA (GS(I,1,2,4),I=1,2)/0.2441785452670D0,-2.372501634110D0/
      DATA (GS(I,2,2,4),I=1,2)/0.0405109766372D0, 2.321661155990D0/
C    5S, EPS =   0.253174778095
      DATA (GS(I,1,2,5),I=1,2)/0.1213425653960D0,-11.18240813640D0/
      DATA (GS(I,2,2,5),I=1,2)/0.0313315214358D0, 10.35350937180D0/
C    E**-X, STO-3G
C    1S, EPS =   0.0000826250458967
      DATA (GS(I,1,3,1),I=1,2)/2.2276605842800D0, 0.355424585054D0/
      DATA (GS(I,2,3,1),I=1,2)/0.4057711561990D0, 0.343751184840D0/
      DATA (GS(I,3,3,1),I=1,2)/0.1098175103870D0, 0.107132229614D0/
C    2S, EPS =   0.0000514945588858
      DATA (GS(I,1,3,2),I=1,2)/2.5815783977100D0,-0.267113527980D0/
      DATA (GS(I,2,3,2),I=1,2)/0.1567622104300D0, 0.324889834023D0/
      DATA (GS(I,3,3,2),I=1,2)/0.0601833227159D0, 0.121807125858D0/
C    3S, EPS =   0.000454031990443
      DATA (GS(I,1,3,3),I=1,2)/0.5641487708760D0,-0.695268461731D0/
      DATA (GS(I,2,3,3),I=1,2)/0.0692442139135D0, 0.696608558492D0/
      DATA (GS(I,3,3,3),I=1,2)/0.0326952909674D0, 0.104204127770D0/
C    4S, EPS =   0.00132827548691
      DATA (GS(I,1,3,4),I=1,2)/0.2267938753100D0,-2.467638561850D0/
      DATA (GS(I,2,3,4),I=1,2)/0.0444817801874D0, 2.294792000980D0/
      DATA (GS(I,3,3,4),I=1,2)/0.0219529466442D0, 0.160684738911D0/
C    5S, EPS =   0.0728662906453
      DATA (GS(I,1,3,5),I=1,2)/0.1080198457610D0,-13.2598240755D0/
      DATA (GS(I,2,3,5),I=1,2)/0.0440811938188D0, 7.64000075806D0/
      DATA (GS(I,3,3,5),I=1,2)/0.0261081181021D0, 4.93625013507D0/
C    E**-X, STO-4G
C    1S, EPS =   0.0000109405425991
      DATA (GS(I,1,4,1),I=1,2)/5.216844533830D0, 0.2474661705050D0/
      DATA (GS(I,2,4,1),I=1,2)/0.954618276012D0, 0.3173636863600D0/
      DATA (GS(I,3,4,1),I=1,2)/0.265203410217D0, 0.2487491794820D0/
      DATA (GS(I,4,4,1),I=1,2)/0.088018627744D0, 0.0595294511502D0/
C    2S, EPS =   0.0000407125756891
*     DATA (GS(I,1,4,2),I=1,2)/2.484931632810D0,-0.2665419421510D0/
*     DATA (GS(I,2,4,2),I=1,2)/0.169226630379D0, 0.2864294179390D0/
*     DATA (GS(I,3,4,2),I=1,2)/0.0723839099766D0,0.1537637592600D0/
*     DATA (GS(I,4,4,2),I=1,2)/0.0334952950572D0,0.00983705948426D0/
C    2S, EPS =   0.0000202300400994
      DATA (GS(I,1,4,2),I=1,2)/11.61525550940D0,-0.1649747320930D0/
      DATA (GS(I,2,4,2),I=1,2)/2.000243111040D0,-0.2013734194590D0/
      DATA (GS(I,3,4,2),I=1,2)/0.160728068746D0, 0.3224443297920D0/
      DATA (GS(I,4,4,2),I=1,2)/0.0612574453247D0,0.1285093183690D0/
C    3S, EPS =   0.00000970551565816
      DATA (GS(I,1,4,3),I=1,2)/1.513265590940D0,-0.2694220373600D0/
      DATA (GS(I,2,4,3),I=1,2)/0.426249750774D0,-0.5451192090540D0/
      DATA (GS(I,3,4,3),I=1,2)/0.0764332086315D0,0.6548914676720D0/
      DATA (GS(I,4,4,3),I=1,2)/0.0376054506346D0,0.1836806351340D0/
C    4S, EPS =   0.0000662328705884
      DATA (GS(I,1,4,4),I=1,2)/0.324221283272D0,-1.0795778820700D0/
      DATA (GS(I,2,4,4),I=1,2)/0.166321717702D0,-1.6614958898700D0/
      DATA (GS(I,3,4,4),I=1,2)/0.0508109745068D0,2.1378627401300D0/
      DATA (GS(I,4,4,4),I=1,2)/0.0282906660023D0,0.5440577852970D0/
C    5S, EPS =   0.000959629615938
      DATA (GS(I,1,4,5),I=1,2)/0.8602284252050D0,1.04839787428D0/
      DATA (GS(I,2,4,5),I=1,2)/0.1189050199600D0,-12.073269348D0/
      DATA (GS(I,3,4,5),I=1,2)/0.0344607617627D0, 10.032225085D0/
      DATA (GS(I,4,4,5),I=1,2)/0.0197479879636D0,0.972000867009D0/
C    E**-X, STO-5G
C    1S, EPS =   1.72100730614D-6
      DATA (GS(I,1,5,1),I=1,2)/11.305636955600D0, 0.1724420895100D0/
      DATA (GS(I,2,5,1),I=1,2)/2.0717281782700D0, 0.2476774170580D0/
      DATA (GS(I,3,5,1),I=1,2)/0.5786484833230D0, 0.2780943503510D0/
      DATA (GS(I,4,5,1),I=1,2)/0.1975724572610D0, 0.1806503130170D0/
      DATA (GS(I,5,5,1),I=1,2)/0.0744527174597D0, 0.0348527400685D0/
C    2S, EPS =   2.02567314643D-6
      DATA (GS(I,1,5,2),I=1,2)/8.9849568621400D0,-0.181263322990D0/
      DATA (GS(I,2,5,2),I=1,2)/1.6737106356300D0,-0.183064196615D0/
      DATA (GS(I,3,5,2),I=1,2)/0.1944726668480D0, 0.236966576685D0/
      DATA (GS(I,4,5,2),I=1,2)/0.0880634563420D0, 0.193848767705D0/
      DATA (GS(I,5,5,2),I=1,2)/0.0424906852183D0, 0.0301550780537D0/
C    3S, EPS =   7.35634915183D-6
*     DATA (GS(I,1,5,3),I=1,2)/1.4666573856200D0,-0.278907270394D0/
*     DATA (GS(I,2,5,3),I=1,2)/0.4171278125700D0,-0.540203709451D0/
*     DATA (GS(I,3,5,3),I=1,2)/0.0790325670632D0, 0.616527228794D0/
*     DATA (GS(I,4,5,3),I=1,2)/0.0407967270196D0, 0.224286935260D0/
*     DATA (GS(I,5,5,3),I=1,2)/0.0189131313634D0, 0.00355595202125D0/
C    3S, EPS =   4.48963930259D-6
      DATA (GS(I,1,5,3),I=1,2)/4.2758779138700D0,-0.0698507234898D0/
      DATA (GS(I,2,5,3),I=1,2)/1.1324094328700D0,-0.2741897960980D0/
      DATA (GS(I,3,5,3),I=1,2)/0.4016256968180D0,-0.4950030686140D0/
      DATA (GS(I,4,5,3),I=1,2)/0.0773237062036D0, 0.6518953128230D0/
      DATA (GS(I,5,5,3),I=1,2)/0.0380070862703D0, 0.1920997391840D0/
C    4S, EPS =   4.74649780005D-6
      DATA (GS(I,1,5,4),I=1,2)/2.9802637828600D0, 0.0769914618562D0/
      DATA (GS(I,2,5,4),I=1,2)/0.3792228833390D0,-0.7927430993770D0/
      DATA (GS(I,3,5,4),I=1,2)/0.1789717223970D0,-1.9394126566400D0/
      DATA (GS(I,4,5,4),I=1,2)/0.0500211036037D0, 2.1419964760100D0/
      DATA (GS(I,5,5,4),I=1,2)/0.0278936168081D0, 0.5041041911000D0/
C    5S, EPS =   1.25477944549D-4
      DATA (GS(I,1,5,5),I=1,2)/0.7403763257480D0, 1.1675860976500D0/
      DATA (GS(I,2,5,5),I=1,2)/0.1367990863440D0,-7.4094063608900D0/
      DATA (GS(I,3,5,5),I=1,2)/0.0913530177928D0,-5.6534566099000D0/
      DATA (GS(I,4,5,5),I=1,2)/0.0372690731451D0, 9.7834354925400D0/
      DATA (GS(I,5,5,5),I=1,2)/0.0224149083571D0, 2.0609517973000D0/
C    E**-X, STO-6G
C    1S, EPS =   3.09302211981D-7
      DATA (GS(I,1,6,1),I=1,2)/23.103031488000D0, 0.1219838776490D0/
      DATA (GS(I,2,6,1),I=1,2)/4.2359155341200D0, 0.1841126016890D0/
      DATA (GS(I,3,6,1),I=1,2)/1.1850565186700D0, 0.2418174734420D0/
      DATA (GS(I,4,6,1),I=1,2)/0.4070988981800D0, 0.2385731276410D0/
      DATA (GS(I,5,6,1),I=1,2)/0.1580884151150D0, 0.1319064303470D0/
      DATA (GS(I,6,6,1),I=1,2)/0.0651095395441D0, 0.0212215075497D0/
C    2S, EPS =   2.45442858687D-7
      DATA (GS(I,1,6,2),I=1,2)/27.684962411900D0,-0.1096251845430D0/
      DATA (GS(I,2,6,2),I=1,2)/5.0771406271200D0,-0.1529697404040D0/
      DATA (GS(I,3,6,2),I=1,2)/1.4267860497100D0,-0.1471118083680D0/
      DATA (GS(I,4,6,2),I=1,2)/0.2040335728690D0, 0.2222710715220D0/
      DATA (GS(I,5,6,2),I=1,2)/0.0926029839916D0, 0.2064584404710D0/
      DATA (GS(I,6,6,2),I=1,2)/0.0441618397753D0, 0.0361066313284D0/
C    3S, EPS =   2.28726313371D-7
      DATA (GS(I,1,6,3),I=1,2)/3.27303193838D0,-0.0987953547155D0/
      DATA (GS(I,2,6,3),I=1,2)/0.920061131072D0,-0.317443161556D0/
      DATA (GS(I,3,6,3),I=1,2)/0.359334976485D0,-0.441582685005D0/
      DATA (GS(I,4,6,3),I=1,2)/0.0863668699103D0,0.528344542293D0/
      DATA (GS(I,5,6,3),I=1,2)/0.047973738124D0,0.308056433582D0/
      DATA (GS(I,6,6,3),I=1,2)/0.0272474114422D0,0.0290285892346D0/
C    4S, EPS =   
*     DATA (GS(I,1,6,4),I=1,2)/
*     DATA (GS(I,2,6,4),I=1,2)/
*     DATA (GS(I,3,6,4),I=1,2)/
*     DATA (GS(I,4,6,4),I=1,2)/
*     DATA (GS(I,5,6,4),I=1,2)/
*     DATA (GS(I,6,6,4),I=1,2)/
C    5S, EPS =   
*     DATA (GS(I,1,6,5),I=1,2)/
*     DATA (GS(I,2,6,5),I=1,2)/
*     DATA (GS(I,3,6,5),I=1,2)/
*     DATA (GS(I,4,6,5),I=1,2)/
*     DATA (GS(I,5,6,5),I=1,2)/
*     DATA (GS(I,6,6,5),I=1,2)/
C    E**-X, STO-7G
C    1S, EPS =   6.18716500055D-8
      DATA (GS(I,1,7,1),I=1,2)/45.0587089151D0,0.0878033213896D0/
      DATA (GS(I,2,7,1),I=1,2)/8.26349795846D0,0.135628748319D0/
      DATA (GS(I,3,7,1),I=1,2)/2.31329257354D0,0.191668575027D0/
      DATA (GS(I,4,7,1),I=1,2)/0.796146258704D0,0.22956857207D0/
      DATA (GS(I,5,7,1),I=1,2)/0.311276302655D0,0.201809727848D0/
      DATA (GS(I,6,7,1),I=1,2)/0.13228069847D0,0.0969446140791D0/
      DATA (GS(I,7,7,1),I=1,2)/0.0582363437708D0,0.013328771963D0/
C    2S, EPS =   1.3299556292D-7
      DATA (GS(I,1,7,2),I=1,2)/81.5049432456D0,-0.065075419859D0/
      DATA (GS(I,2,7,2),I=1,2)/14.8946374315D0,-0.0994996943919D0/
      DATA (GS(I,3,7,2),I=1,2)/4.09885347465D0,-0.135862432653D0/
      DATA (GS(I,4,7,2),I=1,2)/1.35155375512D0,-0.133223664523D0/
      DATA (GS(I,5,7,2),I=1,2)/0.206216296352D0,0.219478418392D0/
      DATA (GS(I,6,7,2),I=1,2)/0.0934934222218D0,0.209090776494D0/
      DATA (GS(I,7,7,2),I=1,2)/0.0444570499317D0,0.0372389578986D0/
C    2S, EPS =   1.57468233D-7
*     DATA (GS(I,1,7,2),I=1,2)/26.1983624959D0,-0.112482603409D0/
*     DATA (GS(I,2,7,2),I=1,2)/4.81709750318D0,-0.155040133949D0/
*     DATA (GS(I,3,7,2),I=1,2)/1.37291259802D0,-0.142634250733D0/
*     DATA (GS(I,4,7,2),I=1,2)/0.217217012254D0,0.190435390815D0/
*     DATA (GS(I,5,7,2),I=1,2)/0.106700768561D0,0.203391900241D0/
*     DATA (GS(I,6,7,2),I=1,2)/0.0565923318777D0,0.0684373454168D0/
*     DATA (GS(I,7,7,2),I=1,2)/0.0308042411091D0,0.00456674915222D0/
C    3S, EPS =   1.4734347876D-8
      DATA (GS(I,1,7,3),I=1,2)/6.84125034972D0,-0.0350159615443D0/
      DATA (GS(I,2,7,3),I=1,2)/1.91431548903D0,-0.136856642905D0/
      DATA (GS(I,3,7,3),I=1,2)/0.737374606655D0,-0.317170368775D0/
      DATA (GS(I,4,7,3),I=1,2)/0.333110061453D0,-0.382779616751D0/
      DATA (GS(I,5,7,3),I=1,2)/0.0895626825518D0,0.494805185372D0/
      DATA (GS(I,6,7,3),I=1,2)/0.0504060819399D0,0.339717892671D0/
      DATA (GS(I,7,7,3),I=1,2)/0.0286214242612D0,0.039823078346D0/
C    4S, EPS =   
*     DATA (GS(I,1,7,4),I=1,2)/
*     DATA (GS(I,2,7,4),I=1,2)/
*     DATA (GS(I,3,7,4),I=1,2)/
*     DATA (GS(I,4,7,4),I=1,2)/
*     DATA (GS(I,5,7,4),I=1,2)/
*     DATA (GS(I,6,7,4),I=1,2)/
*     DATA (GS(I,7,7,4),I=1,2)/
C    5S, EPS =   
*     DATA (GS(I,1,7,5),I=1,2)/
*     DATA (GS(I,2,7,5),I=1,2)/
*     DATA (GS(I,3,7,5),I=1,2)/
*     DATA (GS(I,4,7,5),I=1,2)/
*     DATA (GS(I,5,7,5),I=1,2)/
*     DATA (GS(I,6,7,5),I=1,2)/
*     DATA (GS(I,7,7,5),I=1,2)/
C    E**-X, STO-8G
C    1S, EPS =   1.35230468646D-8
      DATA (GS(I,1,8,1),I=1,2)/84.5781563304D0,0.0642516313276D0/
      DATA (GS(I,2,8,1),I=1,2)/15.5129625354D0,0.10038562421D0/
      DATA (GS(I,3,8,1),I=1,2)/4.34394734098D0,0.147063430328D0/
      DATA (GS(I,4,8,1),I=1,2)/1.49607112612D0,0.193752560198D0/
      DATA (GS(I,5,8,1),I=1,2)/0.586147458716D0,0.213469026225D0/
      DATA (GS(I,6,8,1),I=1,2)/0.250952839465D0,0.169178581227D0/
      DATA (GS(I,7,8,1),I=1,2)/0.114110939146D0,0.0717426601419D0/
      DATA (GS(I,8,8,1),I=1,2)/0.0529406321961D0,0.00858707209406D0/
C    2S, EPS =   1.68541637977D-8
      DATA (GS(I,1,8,2),I=1,2)/67.5836539437D0,-0.071293815016D0/
      DATA (GS(I,2,8,2),I=1,2)/12.3999093931D0,-0.107234321624D0/
      DATA (GS(I,3,8,2),I=1,2)/3.4853335123D0,-0.137685897833D0/
      DATA (GS(I,4,8,2),I=1,2)/1.2207233102D0,-0.118468054236D0/
      DATA (GS(I,5,8,2),I=1,2)/0.231516728766D0,0.166176496234D0/
      DATA (GS(I,6,8,2),I=1,2)/0.11723920827D0,0.206229869312D0/
      DATA (GS(I,7,8,2),I=1,2)/0.0627971189206D0,0.0889860052105D0/
      DATA (GS(I,8,8,2),I=1,2)/0.0340608326311D0,0.00859570484963D0/
C    3S, EPS =   7.0096570638D-9
      DATA (GS(I,1,8,3),I=1,2)/14.7377008331D0,-0.011463199645D0/
      DATA (GS(I,2,8,3),I=1,2)/4.1003253106D0,-0.0496440281041D0/
      DATA (GS(I,3,8,3),I=1,2)/1.55379292122D0,-0.146835941687D0/
      DATA (GS(I,4,8,3),I=1,2)/0.680196113775D0,-0.306080516366D0/
      DATA (GS(I,5,8,3),I=1,2)/0.325513824178D0,-0.361560317306D0/
      DATA (GS(I,6,8,3),I=1,2)/0.0902011524562D0,0.489297795095D0/
      DATA (GS(I,7,8,3),I=1,2)/0.0508040922698D0,0.345551393291D0/
      DATA (GS(I,8,8,3),I=1,2)/0.0288076643661D0,0.0415329558169D0/
C    3S, EPS =   1.07020278787D-8
*     DATA (GS(I,1,8,3),I=1,2)/6.63127796359D0,-0.0366093366878D0/
*     DATA (GS(I,2,8,3),I=1,2)/1.85776081414D0,-0.14188110559D0/
*     DATA (GS(I,3,8,3),I=1,2)/0.718444079828D0,-0.321988854015D0/
*     DATA (GS(I,4,8,3),I=1,2)/0.327614313783D0,-0.374395660119D0/
*     DATA (GS(I,5,8,3),I=1,2)/0.0917455120225D0,0.461211297566D0/
*     DATA (GS(I,6,8,3),I=1,2)/0.0533076578293D0,0.352417489392D0/
*     DATA (GS(I,7,8,3),I=1,2)/0.0317405794153D0,0.0630246274357D0/
*     DATA (GS(I,8,8,3),I=1,2)/0.0176120188384D0,0.000866470114259D0/
C    4S, EPS =   
*     DATA (GS(I,1,8,4),I=1,2)/
*     DATA (GS(I,2,8,4),I=1,2)/
*     DATA (GS(I,3,8,4),I=1,2)/
*     DATA (GS(I,4,8,4),I=1,2)/
*     DATA (GS(I,5,8,4),I=1,2)/
*     DATA (GS(I,6,8,4),I=1,2)/
*     DATA (GS(I,7,8,4),I=1,2)/
*     DATA (GS(I,8,8,4),I=1,2)/
C    5S, EPS =   
*     DATA (GS(I,1,8,5),I=1,2)/
*     DATA (GS(I,2,8,5),I=1,2)/
*     DATA (GS(I,3,8,5),I=1,2)/
*     DATA (GS(I,4,8,5),I=1,2)/
*     DATA (GS(I,5,8,5),I=1,2)/
*     DATA (GS(I,6,8,5),I=1,2)/
*     DATA (GS(I,7,8,5),I=1,2)/
*     DATA (GS(I,8,8,5),I=1,2)/
C    E**-X, STO-9G
C    1S, EPS =   3.18562489233D-9
      DATA (GS(I,1,9,1),I=1,2)/153.728624497D0,0.0477206905702D0/
      DATA (GS(I,2,9,1),I=1,2)/28.1979359285D0,0.0749947953185D0/
      DATA (GS(I,3,9,1),I=1,2)/7.89709680352D0,0.111913605321D0/
      DATA (GS(I,4,9,1),I=1,2)/2.72064070449D0,0.154827763712D0/
      DATA (GS(I,5,9,1),I=1,2)/1.06675720199D0,0.19149125807D0/
      DATA (GS(I,6,9,1),I=1,2)/0.457779546063D0,0.195508945403D0/
      DATA (GS(I,7,9,1),I=1,2)/0.209837187888D0,0.141014529442D0/
      DATA (GS(I,8,9,1),I=1,2)/0.100627425452D0,0.0534591749349D0/
      DATA (GS(I,9,9,1),I=1,2)/0.0487180150656D0,0.00565173408344D0/
C    2S, EPS =   2.77808884798D-9
      DATA (GS(I,1,9,2),I=1,2)/162.52006988D0,-0.046274483804D0/
      DATA (GS(I,2,9,2),I=1,2)/29.8089984575D0,-0.0717507122962D0/
      DATA (GS(I,3,9,2),I=1,2)/8.35081369071D0,-0.102416186932D0/
      DATA (GS(I,4,9,2),I=1,2)/2.88348272961D0,-0.124405650264D0/
      DATA (GS(I,5,9,2),I=1,2)/1.14030781225D0,-0.10393128104D0/
      DATA (GS(I,6,9,2),I=1,2)/0.238520383516D0,0.156228993047D0/
      DATA (GS(I,7,9,2),I=1,2)/0.121621719219D0,0.208360630747D0/
      DATA (GS(I,8,9,2),I=1,2)/0.0649005936477D0,0.0969036859144D0/
      DATA (GS(I,9,9,2),I=1,2)/0.0349352459479D0,0.0100577540405D0/
C    3S, EPS =   6.7098376694D-10
      DATA (GS(I,1,9,3),I=1,2)/12.5194657203D0,-0.0145421135843D0/
      DATA (GS(I,2,9,3),I=1,2)/3.50407264575D0,-0.0611537719673D0/
      DATA (GS(I,3,9,3),I=1,2)/1.35035980411D0,-0.168928602286D0/
      DATA (GS(I,4,9,3),I=1,2)/0.610812875422D0,-0.314634026746D0/
      DATA (GS(I,5,9,3),I=1,2)/0.306840656597D0,-0.325313063027D0/
      DATA (GS(I,6,9,3),I=1,2)/0.0961843396953D0,0.40532108712D0/
      DATA (GS(I,7,9,3),I=1,2)/0.0580690564327D0,0.372146371189D0/
      DATA (GS(I,8,9,3),I=1,2)/0.0360307483371D0,0.102559470245D0/
      DATA (GS(I,9,9,3),I=1,2)/0.0222722687361D0,0.00556566186546D0/
C    4S, EPS =   
*     DATA (GS(I,1,9,4),I=1,2)/
*     DATA (GS(I,2,9,4),I=1,2)/
*     DATA (GS(I,3,9,4),I=1,2)/
*     DATA (GS(I,4,9,4),I=1,2)/
*     DATA (GS(I,5,9,4),I=1,2)/
*     DATA (GS(I,6,9,4),I=1,2)/
*     DATA (GS(I,7,9,4),I=1,2)/
*     DATA (GS(I,8,9,4),I=1,2)/
*     DATA (GS(I,9,9,4),I=1,2)/
C    5S, EPS =   
*     DATA (GS(I,1,9,5),I=1,2)/
*     DATA (GS(I,2,9,5),I=1,2)/
*     DATA (GS(I,3,9,5),I=1,2)/
*     DATA (GS(I,4,9,5),I=1,2)/
*     DATA (GS(I,5,9,5),I=1,2)/
*     DATA (GS(I,6,9,5),I=1,2)/
*     DATA (GS(I,7,9,5),I=1,2)/
*     DATA (GS(I,8,9,5),I=1,2)/
*     DATA (GS(I,9,9,5),I=1,2)/
C    Output scalars
      NO = 2*L+1
      NR = IGO
      NA = L
C    Simple range checks
      IF( NA.GT.MAXA .OR. NR.GT.MAXR ) THEN
          WRITE(NB6,11010) NA, NR, MAXA, MAXR
          STOP 'PSOGMK'
      ENDIF
C    Generate radial part.
      IRPOW = N-L
      IF( IRPOW.LT.1 .OR. IRPOW.GT.MAXPOW .OR. 
     .                    IGO.LT.1 .OR. IGO.GT.MAXGO ) THEN
          WRITE(NB6,11020) IGO, N, L
          STOP 'PSOGMK'
      ENDIF
C    All Gaussian exponents should be non-zero, and should
C    be arranged in the decreasing order.
      DO 800 I=1,IGO
          IF( GS(1,I,IGO,IRPOW).LE.ZERO ) THEN
              WRITE(NB6,11030) I, IGO, IRPOW
              STOP 'PSOGMK'
          ENDIF
  800 CONTINUE
      DO 850 I=2,IGO
          IF( GS(1,I-1,IGO,IRPOW).LE.GS(1,I,IGO,IRPOW) ) THEN
              WRITE(NB6,11040) I-1, I, IGO, IRPOW
              STOP 'PSOGMK'
          ENDIF
  850 CONTINUE
C    Expansion seems to be valid, process it.
      SCA = ZETA**2
      SCD = PSSNRM(N,ZETA)/ZETA**(N-L-1)
      DO 900 I=1,IGO
          ALP(I) = SCA*GS(1,I,IGO,IRPOW)
          DQ (I) = SCD*GS(2,I,IGO,IRPOW)
  900 CONTINUE
C    Generate angular part
C    Most of the coefficients are zero
      DO 2000 IO=1,NO
          DO 1900 IZ=0,NA
              DO 1800 IY=0,NA-IZ
                  DO 1700 IX=0,NA-(IY+IZ)
                      CS(IX,IY,IZ,1,IO) = ZERO
 1700             CONTINUE
 1800         CONTINUE
 1900     CONTINUE
 2000 CONTINUE
C    Fill in non-zero angular coefficients for the first radial
C    component. Rest will be cloned aftewards.
      IF( L.EQ.0 ) THEN
          CS(0,0,0,1,1) = C1
      ELSE IF( L.EQ.1 ) THEN
C        M=-1 (y)
          CS(0,1,0,1,1) = C2
C        M= 0 (z)
          CS(0,0,1,1,2) = C2
C        M=+1 (x)
          CS(1,0,0,1,3) = C2
      ELSE IF( L.EQ.2 ) THEN
C        M=-2 (xy)
          CS(1,1,0,1,1) = C3
C        M=-1 (yz)
          CS(0,1,1,1,2) = C3
C        M= 0 (z^2)
          CS(2,0,0,1,3) =-C4
          CS(0,2,0,1,3) =-C4
          CS(0,0,2,1,3) = C5
C        M=+1 (xz)
          CS(1,0,1,1,4) = C3
C        M=+2 (x^2-y^2)
          CS(2,0,0,1,5) = C6
          CS(0,2,0,1,5) =-C6
      ELSE IF( L.EQ.3 ) THEN
C        M=-3 (xy(x-y))
          CS(2,1,0,1,1) = C7
          CS(0,3,0,1,1) =-C8
C        M=-2 (xyz)
          CS(1,1,1,1,2) = C9
C        M=-1 (y z^2)
          CS(2,1,0,1,3) =-C10
          CS(0,3,0,1,3) =-C10
          CS(0,1,2,1,3) = C11
C        M= 0 (z^3)
          CS(2,0,1,1,4) =-C12
          CS(0,2,1,1,4) =-C12
          CS(0,0,3,1,4) = C13
C        M=+1 (x z^2)
          CS(3,0,0,1,5) =-C10
          CS(1,2,0,1,5) =-C10
          CS(1,0,2,1,5) = C11
C        M=+2 (z(x^2-y^2))
          CS(2,0,1,1,6) = C14
          CS(0,2,1,1,6) =-C14
C        M=+3 (x(x^2-y^2))
          CS(3,0,0,1,7) = C8
          CS(1,2,0,1,7) =-C7
      ELSE
          WRITE(NB6,11050) L
          STOP 'PSOGMK'
      ENDIF
C    Replicate C coefficients over radial components
      DO 5000 IO=1,NO
          DO 4900 IR=NR,1,-1
              DO 4800 IZ=0,NA
                  DO 4700 IY=0,NA-IZ
                      DO 4600 IX=0,NA-(IY+IZ)
                          CS(IX,IY,IZ,IR,IO) = DQ(IR)*CS(IX,IY,IZ,1,IO)
 4600                 CONTINUE
 4700             CONTINUE
 4800         CONTINUE
 4900     CONTINUE
 5000 CONTINUE
C    Done!
      RETURN
11010 FORMAT(' NA = ',I4,' OR NR = ',I4,' IS TOO BIG FOR MAXA = ',I4,
     .       ' MAXR = ',I4,' IN PSOGMK.')
11020 FORMAT(' STO-',I2,'G EXPANSION IS NOT AVAILABLE IN PSOGMK ',
     .       'FOR N = ',I4,' L = ',I4)
11030 FORMAT(' GAUSSIAN EXPONENT FOR COMPONENT ',I2,' OF STO-',I1,
     .       'G EXPANSION WITH N-L= ',I2,' IS NON-POSITIVE IN PSOGMK')
11040 FORMAT(' COMPONENTS ',I2,' AND ',I2,' OF STO-',I1,
     .       'G EXPANSION WITH N-L= ',I2,' ARE OUT OF ORDER IN PSOGMK')
11050 FORMAT(' ANGULAR PART IS NOT AVAILABLE FOR L = ',I4,' IN PSOGMK')
      END
C
      SUBROUTINE PSOGZR(NO,NR,NA,CS,MAXA,MAXR)
C
C   Zeros out angular components.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      See top.
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Possible optimizations:
C
C      Static array dimensions and cloning should improve performance.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0)
C
      DIMENSION CS(0:MAXA,0:MAXA,0:MAXA,MAXR,NO)
C
      DO 800 IO=1,NO
          DO 700 IR=1,NR
              DO 600 IZ=0,NA
                  DO 500 IY=0,NA-IZ
                      DO 400 IX=0,NA-(IZ+IY)
                          CS(IX,IY,IZ,IR,IO) = ZERO
  400                 CONTINUE
  500             CONTINUE
  600         CONTINUE
  700     CONTINUE
  800 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSOGRS(NOI,NRI,NAI,CSI,ALPI,MAXAI,MAXRI,
     .                  NOO,NRO,NAO,CSO,ALPO,MAXAO,MAXRO)
C
C   Computes a product of r operator and angular part of a nest.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      *I     - Input orbitals nest
C      *O     - Output orbitals nest
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      On output, 3*NOI output orbitals are ordered as in (NOI,3) array.
C
C   Possible optimizations:
C
C      Compiled-in array bounds and cloning should improve performance.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION CSI(0:MAXAI,0:MAXAI,0:MAXAI,MAXRI,NOI),   ALPI(NRI)
      DIMENSION CSO(0:MAXAO,0:MAXAO,0:MAXAO,MAXRO,NOI,3), ALPO(MAXRO)
C    Output scalars are range checks
      NOO = 3*NOI
      NRO = NRI
      NAO = NAI+1
      IF( NRO.GT.MAXRO .OR. NAO.GT.MAXAO ) THEN
          WRITE(NB6,11010) NRO, NAO, MAXRO, MAXAO
          STOP 'PSOGRS'
      ENDIF
C    Orbital exponents won't change
      DO 100 IR=1,NRI
          ALPO(IR) = ALPI(IR)
  100 CONTINUE
C    Zeroing out all output array is overkill, but it is easier this way
      CALL PSOGZR(NOO,NRO,NAO,CSO,MAXAO,MAXRO)
C
      DO 1800 IO=1,NOI
          DO 1700 IR=1,NRI
              DO 1600 IZ=0,NAI
                  DO 1500 IY=0,NAI-IZ
                      DO 1400 IX=0,NAI-(IZ+IY)
                          VAL = CSI(IX,IY,IZ,IR,IO)
                          CSO(IX+1,IY+0,IZ+0,IR,IO,1) = VAL
                          CSO(IX+0,IY+1,IZ+0,IR,IO,2) = VAL
                          CSO(IX+0,IY+0,IZ+1,IR,IO,3) = VAL
 1400                 CONTINUE
 1500             CONTINUE
 1600         CONTINUE
 1700     CONTINUE
 1800 CONTINUE
C
      RETURN
11010 FORMAT(' NRO = ',I4,' OR NAO = ',I4,' IS TOO BIG FOR MAXRO = ',I4,
     .       ' MAXAO = ',I4,' IN PS OGRESS')
      END
C
      SUBROUTINE PSOGDX(NOI,NRI,NAI,CSI,ALPI,MAXAI,MAXRI,
     .                  NOO,NRO,NAO,CSO,ALPO,MAXAO,MAXRO)
C
C   Applies gradient operator to an orbitals nest.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      *I     - Input orbitals nest.
C      *O     - Output orbitals nest.
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Output orbitals are ordered as in (NOI,3) array.
C
C   Possible optimizations:
C
C      Compiled-in array dimensions and cloning should improve performance.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TWO =2.D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C    Parameters
      DIMENSION CSI(0:MAXAI,0:MAXAI,0:MAXAI,MAXRI,NOI),   ALPI(NRI)
      DIMENSION CSO(0:MAXAO,0:MAXAO,0:MAXAO,MAXRO,NOI,3), ALPO(MAXRO)
C    Output scalars are range checks
      NOO = 3*NOI
      NRO = NRI
      NAO = NAI+1
      IF( NRO.GT.MAXRO .OR. NAO.GT.MAXAO ) THEN
          WRITE(NB6,11010) NRO, NAO, MAXRO, MAXAO
          STOP 'PSOGDX'
      ENDIF
C    Orbital exponents won't change
      DO 100 IR=1,NRI
          ALPO(IR) = ALPI(IR)
  100 CONTINUE
C    Zero output
      CALL PSOGZR(NOO,NRO,NAO,CSO,MAXAO,MAXRO)
C
      DO 1800 IO=1,NOI
          DO 1700 IR=1,NRI
              SC = -TWO*ALPI(IR)
              DO 1600 IZ=0,NAI
                  DO 1500 IY=0,NAI-IZ
C                    Radial component
                      DO 1400 IX=0,NAI-(IZ+IY)
                          CS  = CSI(IX,IY,IZ,IR,IO)
                          VAL = SC * CS
C                        Radial component
                          CSO(IX+1,IY+0,IZ+0,IR,IO,1) = VAL
                          CSO(IX+0,IY+1,IZ+0,IR,IO,2) = VAL
                          CSO(IX+0,IY+0,IZ+1,IR,IO,3) = VAL
C                        Angular component
                          IF(IX.GT.0) CSO(IX-1,IY,IZ,IR,IO,1) = 
     .                                CSO(IX-1,IY,IZ,IR,IO,1) + IX*CS
                          IF(IY.GT.0) CSO(IX,IY-1,IZ,IR,IO,2) = 
     .                                CSO(IX,IY-1,IZ,IR,IO,2) + IY*CS
                          IF(IZ.GT.0) CSO(IX,IY,IZ-1,IR,IO,3) = 
     .                                CSO(IX,IY,IZ-1,IR,IO,3) + IZ*CS
 1400                 CONTINUE
 1500             CONTINUE
 1600         CONTINUE
 1700     CONTINUE
 1800 CONTINUE
C
      RETURN
11010 FORMAT(' NRO = ',I4,' OR NAO = ',I4,' IS TOO BIG FOR MAXRO = ',I4,
     .       ' MAXAO = ',I4,' IN PSOGDX')
      END
C
      SUBROUTINE PSOGTR(R,
     .                  NOI,NRI,NAI,CSI,ALPI,MAXAI,MAXRI,
     .                  NOO,NRO,NAO,CSO,ALPO,MAXAO,MAXRO)
C
C   Translates angular part of the orbitals nest to coordinate
C   origin.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      R      - Original position of the orbitals nest, vector
C               of the length 3.
C      *I     - Input orbitals nest.
C      *O     - Output orbitals nest.
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Output orbitals are not in the canonical form anymore, since
C      only angular part is translated. Radial part still refers to
C      the original position.
C
C      This function should perform reasonably well only if CSI is
C      fairly sparse.
C
C   Possible optimizations:
C
C      Static array bounds and cloning should improve performance.
C
C   Bugs:
C
C      Scaling of the algorithm we use here is frightful.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IOMAX=10)
      PARAMETER (EPS=1.D-15)
      PARAMETER (ONE=1.D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C    Parameters
      DIMENSION R(3)
      DIMENSION CSI(0:MAXAI,0:MAXAI,0:MAXAI,MAXRI,NOI), ALPI(NRI)
      DIMENSION CSO(0:MAXAO,0:MAXAO,0:MAXAO,MAXRO,NOI), ALPO(MAXRO)
C    Locals
      DIMENSION FACT(0:IOMAX), RFACT(0:IOMAX), BINOM(0:IOMAX,0:IOMAX)
      DIMENSION RPOW(3,0:IOMAX), CX(0:IOMAX,0:IOMAX)
      DIMENSION CY(0:IOMAX,0:IOMAX), CZ(0:IOMAX,0:IOMAX)
C    Output scalars are range checks
      NOO = NOI
      NRO = NRI
      NAO = NAI
      IF( NRO.GT.MAXRO .OR. NAO.GT.MAXAO ) THEN
          WRITE(NB6,11010) NRO, NAO, MAXRO, MAXAO
          STOP 'PSOGTR'
      ENDIF
      IF( NAI.GT.IOMAX ) THEN
          WRITE(NB6,11020) NAI, IOMAX
          STOP 'PSOGTR'
      ENDIF
C    Radial exponents wont change
      DO 100 IR=1,NRI
          ALPO(IR) = ALPI(IR)
  100 CONTINUE
C    Factorials and binomial coefficients
      FACT (0) = ONE
      RFACT(0) = ONE
      DO 200 I=1,NAI
          FACT (I) = I*FACT(I-1)
          RFACT(I) = ONE/FACT(I)
  200 CONTINUE
      DO 300 J=0,NAI
          DO 290 I=0,J
              BINOM(I,J) = FACT(J)*RFACT(I)*RFACT(J-I)
  2900    CONTINUE
  300 CONTINUE
C    Powers of coordinates of the orbital nest
      RPOW(1,0) = ONE
      RPOW(2,0) = ONE
      RPOW(3,0) = ONE
      DO 400 I=1,NAI
          RPOW(1,I) = -R(1)*RPOW(1,I-1)
          RPOW(2,I) = -R(2)*RPOW(2,I-1)
          RPOW(3,I) = -R(3)*RPOW(3,I-1)
  400 CONTINUE
C    Individual translation coefficients for powers of cartesian
C    coordinates
      DO 500 N=0,NAI
          DO 490 I=0,N
              CX(I,N) = BINOM(I,N)*RPOW(1,N-I)
              CY(I,N) = BINOM(I,N)*RPOW(2,N-I)
              CZ(I,N) = BINOM(I,N)*RPOW(3,N-I)
  490     CONTINUE
  500 CONTINUE
C    Zero output
      CALL PSOGZR(NOO,NRO,NAO,CSO,MAXAO,MAXRO)
C    Nest translation.
      DO 2900 IO=1,NOI
        DO 2800 IR=1,NRI
          DO 2700 IZ=0,NAI
            DO 2600 IY=0,NAI-IZ
              DO 2500 IX=0,NAI-(IY+IZ)
                  SC0  = CSI(IX,IY,IZ,IR,IO)
                  IF( ABS(SC0).GT.EPS ) THEN
                      DO 1900 IZO=0,IZ
                        SC1 = SC0*CZ(IZO,IZ)
                        DO 1800 IYO=0,IY
                          SC2 = SC1*CY(IYO,IY)
                          DO 1700 IXO=0,IX
                              CSO(IXO,IYO,IZO,IR,IO) = 
     .                        CSO(IXO,IYO,IZO,IR,IO) + SC2*CX(IXO,IX)
 1700                     CONTINUE
 1800                   CONTINUE
 1900                 CONTINUE
                  ENDIF
 2500         CONTINUE
 2600       CONTINUE
 2700     CONTINUE
 2800   CONTINUE
 2900 CONTINUE
C
      RETURN
11010 FORMAT(' NRO = ',I4,' OR NAO = ',I4,' IS TOO BIG FOR MAXRO = ',I4,
     .       ' MAXAO = ',I4,' IN PSOGTR')
11020 FORMAT(' NAI = ',I4,' IS TOO BIG IN PSOGTR. LIMIT IS ',I4)
      END
C
      SUBROUTINE PSOGPR(NA1,CS1,NA2,CS2,NAO,CSO)
C
C   Computes a product of the angular parts of two gaussian
C   components.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NA1    - Maximum order of the first orbital
C      CS1    - Coefficients of the first orbital
C      NA2    - Maximum order of the second orbital
C      CS2    - Coefficients of the second orbital
C      NAO    - (Output) maximum order of the product orbital
C      CSO    - (Output) coefficient of the  product orbital
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Possible optimizations:
c
C      This routine is already tuned fairly heavily. The only reasonable
C      optimization at this point is inlining of the whole beast.
C
C   Bugs:
C
C      Dimensions of of the coefficient arrays are static to reduce
C      address computation overhead.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0)
      PARAMETER (MAXAI=6,MAXAO=12)
C    Arguments
      DIMENSION CS1(0:MAXAI,0:MAXAI,0:MAXAI)
      DIMENSION CS2(0:MAXAI,0:MAXAI,0:MAXAI)
      DIMENSION CSO(0:MAXAO,0:MAXAO,0:MAXAO)
C    Output orders and range checks
      NAO = NA1 + NA2
C    Zero output coefficients
      DO 190 IZO=0,NAO
          DO 180 IYO=0,NAO-IZO
              DO 170 IXO=0,NAO-(IZO+IYO)
                  CSO(IXO,IYO,IZO) = ZERO
  170         CONTINUE
  180     CONTINUE
  190 CONTINUE
C    Evaluate products. These hideous nested loops are expensive, so 
C    we'd unroll several copies of the three nested ones.
      IF( NA2.GT.4 ) GOTO 9000
      ICASE = MIN(NA1+1,6)
      GOTO (1000,1010,1020,1030,1040,1090,
     .      2000,2010,2020,2030,2090,2090,
     .      3000,3010,3020,3090,3090,3090,
     .      4000,4010,4090,4090,4090,4090,
     .      5000,5090,5090,5090,5090,5090), NA2*6 + ICASE
C    We shouldn't get here
      STOP 'PSOGPR/??'
C    NA1 = 0, NA2 = 0
 1000 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      GOTO 9999
C    NA1 = 1, NA2 = 0
 1010 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      CSO(0,0,1) = CS1(0,0,1)*CS2(0,0,0)
      CSO(0,1,0) = CS1(0,1,0)*CS2(0,0,0)
      CSO(1,0,0) = CS1(1,0,0)*CS2(0,0,0)
      GOTO 9999
C    NA1 = 2, NA2 = 0
 1020 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      CSO(0,0,1) = CS1(0,0,1)*CS2(0,0,0)
      CSO(0,0,2) = CS1(0,0,2)*CS2(0,0,0)
      CSO(0,1,0) = CS1(0,1,0)*CS2(0,0,0)
      CSO(0,1,1) = CS1(0,1,1)*CS2(0,0,0)
      CSO(0,2,0) = CS1(0,2,0)*CS2(0,0,0)
      CSO(1,0,0) = CS1(1,0,0)*CS2(0,0,0)
      CSO(1,0,1) = CS1(1,0,1)*CS2(0,0,0)
      CSO(1,1,0) = CS1(1,1,0)*CS2(0,0,0)
      CSO(2,0,0) = CS1(2,0,0)*CS2(0,0,0)
      GOTO 9999
C    NA1 = 3, NA2 = 0
 1030 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      CSO(0,0,1) = CS1(0,0,1)*CS2(0,0,0)
      CSO(0,0,2) = CS1(0,0,2)*CS2(0,0,0)
      CSO(0,0,3) = CS1(0,0,3)*CS2(0,0,0)
      CSO(0,1,0) = CS1(0,1,0)*CS2(0,0,0)
      CSO(0,1,1) = CS1(0,1,1)*CS2(0,0,0)
      CSO(0,1,2) = CS1(0,1,2)*CS2(0,0,0)
      CSO(0,2,0) = CS1(0,2,0)*CS2(0,0,0)
      CSO(0,2,1) = CS1(0,2,1)*CS2(0,0,0)
      CSO(0,3,0) = CS1(0,3,0)*CS2(0,0,0)
      CSO(1,0,0) = CS1(1,0,0)*CS2(0,0,0)
      CSO(1,0,1) = CS1(1,0,1)*CS2(0,0,0)
      CSO(1,0,2) = CS1(1,0,2)*CS2(0,0,0)
      CSO(1,1,0) = CS1(1,1,0)*CS2(0,0,0)
      CSO(1,1,1) = CS1(1,1,1)*CS2(0,0,0)
      CSO(1,2,0) = CS1(1,2,0)*CS2(0,0,0)
      CSO(2,0,0) = CS1(2,0,0)*CS2(0,0,0)
      CSO(2,0,1) = CS1(2,0,1)*CS2(0,0,0)
      CSO(2,1,0) = CS1(2,1,0)*CS2(0,0,0)
      CSO(3,0,0) = CS1(3,0,0)*CS2(0,0,0)
      GOTO 9999
C    NA1 = 4, NA2 = 0
 1040 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      CSO(0,0,1) = CS1(0,0,1)*CS2(0,0,0)
      CSO(0,0,2) = CS1(0,0,2)*CS2(0,0,0)
      CSO(0,0,3) = CS1(0,0,3)*CS2(0,0,0)
      CSO(0,0,4) = CS1(0,0,4)*CS2(0,0,0)
      CSO(0,1,0) = CS1(0,1,0)*CS2(0,0,0)
      CSO(0,1,1) = CS1(0,1,1)*CS2(0,0,0)
      CSO(0,1,2) = CS1(0,1,2)*CS2(0,0,0)
      CSO(0,1,3) = CS1(0,1,3)*CS2(0,0,0)
      CSO(0,2,0) = CS1(0,2,0)*CS2(0,0,0)
      CSO(0,2,1) = CS1(0,2,1)*CS2(0,0,0)
      CSO(0,2,2) = CS1(0,2,2)*CS2(0,0,0)
      CSO(0,3,0) = CS1(0,3,0)*CS2(0,0,0)
      CSO(0,3,1) = CS1(0,3,1)*CS2(0,0,0)
      CSO(0,4,0) = CS1(0,4,0)*CS2(0,0,0)
      CSO(1,0,0) = CS1(1,0,0)*CS2(0,0,0)
      CSO(1,0,1) = CS1(1,0,1)*CS2(0,0,0)
      CSO(1,0,2) = CS1(1,0,2)*CS2(0,0,0)
      CSO(1,0,3) = CS1(1,0,3)*CS2(0,0,0)
      CSO(1,1,0) = CS1(1,1,0)*CS2(0,0,0)
      CSO(1,1,1) = CS1(1,1,1)*CS2(0,0,0)
      CSO(1,1,2) = CS1(1,1,2)*CS2(0,0,0)
      CSO(1,2,0) = CS1(1,2,0)*CS2(0,0,0)
      CSO(1,2,1) = CS1(1,2,1)*CS2(0,0,0)
      CSO(1,3,0) = CS1(1,3,0)*CS2(0,0,0)
      CSO(2,0,0) = CS1(2,0,0)*CS2(0,0,0)
      CSO(2,0,1) = CS1(2,0,1)*CS2(0,0,0)
      CSO(2,0,2) = CS1(2,0,2)*CS2(0,0,0)
      CSO(2,1,0) = CS1(2,1,0)*CS2(0,0,0)
      CSO(2,1,1) = CS1(2,1,1)*CS2(0,0,0)
      CSO(2,2,0) = CS1(2,2,0)*CS2(0,0,0)
      CSO(3,0,0) = CS1(3,0,0)*CS2(0,0,0)
      CSO(3,0,1) = CS1(3,0,1)*CS2(0,0,0)
      CSO(3,1,0) = CS1(3,1,0)*CS2(0,0,0)
      CSO(4,0,0) = CS1(4,0,0)*CS2(0,0,0)
      GOTO 9999
C    NA2 = 0
 1090 CONTINUE
      DO 1490 IZ1=0,NA1
        DO 1480 IY1=0,NA1-IZ1
          DO 1470 IX1=0,NA1-(IZ1+IY1)
             C1 = CS1(IX1,IY1,IZ1)
             CSO(IX1+0,IY1+0,IZ1+0)=CSO(IX1+0,IY1+0,IZ1+0)+C1*CS2(0,0,0)
 1470     CONTINUE
 1480   CONTINUE
 1490 CONTINUE
      GOTO 9999
C    NA1 = 0, NA2 = 1
 2000 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      CSO(0,0,1) = CS1(0,0,0)*CS2(0,0,1)
      CSO(0,1,0) = CS1(0,0,0)*CS2(0,1,0)
      CSO(1,0,0) = CS1(0,0,0)*CS2(1,0,0)
      GOTO 9999
C    NA1 = 1, NA2 = 1
 2010 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      CSO(0,0,1) = CS1(0,0,0)*CS2(0,0,1) + CS1(0,0,1)*CS2(0,0,0)
      CSO(0,0,2) = CS1(0,0,1)*CS2(0,0,1)
      CSO(0,1,0) = CS1(0,0,0)*CS2(0,1,0) + CS1(0,1,0)*CS2(0,0,0)
      CSO(0,1,1) = CS1(0,0,1)*CS2(0,1,0) + CS1(0,1,0)*CS2(0,0,1)
      CSO(0,2,0) = CS1(0,1,0)*CS2(0,1,0)
      CSO(1,0,0) = CS1(0,0,0)*CS2(1,0,0) + CS1(1,0,0)*CS2(0,0,0)
      CSO(1,0,1) = CS1(0,0,1)*CS2(1,0,0) + CS1(1,0,0)*CS2(0,0,1)
      CSO(1,1,0) = CS1(0,1,0)*CS2(1,0,0) + CS1(1,0,0)*CS2(0,1,0)
      CSO(2,0,0) = CS1(1,0,0)*CS2(1,0,0)
      GOTO 9999
C    NA1 = 2, NA2 = 1
 2020 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      CSO(0,0,1) = CS1(0,0,0)*CS2(0,0,1) + CS1(0,0,1)*CS2(0,0,0)
      CSO(0,0,2) = CS1(0,0,1)*CS2(0,0,1) + CS1(0,0,2)*CS2(0,0,0)
      CSO(0,0,3) = CS1(0,0,2)*CS2(0,0,1)
      CSO(0,1,0) = CS1(0,0,0)*CS2(0,1,0) + CS1(0,1,0)*CS2(0,0,0)
      CSO(0,1,1) = CS1(0,0,1)*CS2(0,1,0) + CS1(0,1,0)*CS2(0,0,1)
     .           + CS1(0,1,1)*CS2(0,0,0)
      CSO(0,1,2) = CS1(0,0,2)*CS2(0,1,0) + CS1(0,1,1)*CS2(0,0,1)
      CSO(0,2,0) = CS1(0,1,0)*CS2(0,1,0) + CS1(0,2,0)*CS2(0,0,0)
      CSO(0,2,1) = CS1(0,1,1)*CS2(0,1,0) + CS1(0,2,0)*CS2(0,0,1)
      CSO(0,3,0) = CS1(0,2,0)*CS2(0,1,0)
      CSO(1,0,0) = CS1(0,0,0)*CS2(1,0,0) + CS1(1,0,0)*CS2(0,0,0)
      CSO(1,0,1) = CS1(0,0,1)*CS2(1,0,0) + CS1(1,0,0)*CS2(0,0,1)
     .           + CS1(1,0,1)*CS2(0,0,0)
      CSO(1,0,2) = CS1(0,0,2)*CS2(1,0,0) + CS1(1,0,1)*CS2(0,0,1)
      CSO(1,1,0) = CS1(0,1,0)*CS2(1,0,0) + CS1(1,0,0)*CS2(0,1,0)
     .           + CS1(1,1,0)*CS2(0,0,0)
      CSO(1,1,1) = CS1(0,1,1)*CS2(1,0,0) + CS1(1,0,1)*CS2(0,1,0)
     .           + CS1(1,1,0)*CS2(0,0,1)
      CSO(1,2,0) = CS1(0,2,0)*CS2(1,0,0) + CS1(1,1,0)*CS2(0,1,0)
      CSO(2,0,0) = CS1(1,0,0)*CS2(1,0,0) + CS1(2,0,0)*CS2(0,0,0)
      CSO(2,0,1) = CS1(1,0,1)*CS2(1,0,0) + CS1(2,0,0)*CS2(0,0,1)
      CSO(2,1,0) = CS1(1,1,0)*CS2(1,0,0) + CS1(2,0,0)*CS2(0,1,0)
      CSO(3,0,0) = CS1(2,0,0)*CS2(1,0,0)
      GOTO 9999
C    NA1 = 3, NA2 = 1
 2030 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      CSO(0,0,1) = CS1(0,0,0)*CS2(0,0,1) + CS1(0,0,1)*CS2(0,0,0)
      CSO(0,0,2) = CS1(0,0,1)*CS2(0,0,1) + CS1(0,0,2)*CS2(0,0,0)
      CSO(0,0,3) = CS1(0,0,2)*CS2(0,0,1) + CS1(0,0,3)*CS2(0,0,0)
      CSO(0,0,4) = CS1(0,0,3)*CS2(0,0,1)
      CSO(0,1,0) = CS1(0,0,0)*CS2(0,1,0) + CS1(0,1,0)*CS2(0,0,0)
      CSO(0,1,1) = CS1(0,0,1)*CS2(0,1,0) + CS1(0,1,0)*CS2(0,0,1)
     .           + CS1(0,1,1)*CS2(0,0,0)
      CSO(0,1,2) = CS1(0,0,2)*CS2(0,1,0) + CS1(0,1,1)*CS2(0,0,1)
     .           + CS1(0,1,2)*CS2(0,0,0)
      CSO(0,1,3) = CS1(0,0,3)*CS2(0,1,0) + CS1(0,1,2)*CS2(0,0,1)
      CSO(0,2,0) = CS1(0,1,0)*CS2(0,1,0) + CS1(0,2,0)*CS2(0,0,0)
      CSO(0,2,1) = CS1(0,1,1)*CS2(0,1,0) + CS1(0,2,0)*CS2(0,0,1)
     .           + CS1(0,2,1)*CS2(0,0,0)
      CSO(0,2,2) = CS1(0,1,2)*CS2(0,1,0) + CS1(0,2,1)*CS2(0,0,1)
      CSO(0,3,0) = CS1(0,2,0)*CS2(0,1,0) + CS1(0,3,0)*CS2(0,0,0)
      CSO(0,3,1) = CS1(0,2,1)*CS2(0,1,0) + CS1(0,3,0)*CS2(0,0,1)
      CSO(0,4,0) = CS1(0,3,0)*CS2(0,1,0)
      CSO(1,0,0) = CS1(0,0,0)*CS2(1,0,0) + CS1(1,0,0)*CS2(0,0,0)
      CSO(1,0,1) = CS1(0,0,1)*CS2(1,0,0) + CS1(1,0,0)*CS2(0,0,1)
     .           + CS1(1,0,1)*CS2(0,0,0)
      CSO(1,0,2) = CS1(0,0,2)*CS2(1,0,0) + CS1(1,0,1)*CS2(0,0,1)
     .           + CS1(1,0,2)*CS2(0,0,0)
      CSO(1,0,3) = CS1(0,0,3)*CS2(1,0,0) + CS1(1,0,2)*CS2(0,0,1)
      CSO(1,1,0) = CS1(0,1,0)*CS2(1,0,0) + CS1(1,0,0)*CS2(0,1,0)
     .           + CS1(1,1,0)*CS2(0,0,0)
      CSO(1,1,1) = CS1(0,1,1)*CS2(1,0,0) + CS1(1,0,1)*CS2(0,1,0)
     .           + CS1(1,1,0)*CS2(0,0,1) + CS1(1,1,1)*CS2(0,0,0)
      CSO(1,1,2) = CS1(0,1,2)*CS2(1,0,0) + CS1(1,0,2)*CS2(0,1,0)
     .           + CS1(1,1,1)*CS2(0,0,1)
      CSO(1,2,0) = CS1(0,2,0)*CS2(1,0,0) + CS1(1,1,0)*CS2(0,1,0)
     .           + CS1(1,2,0)*CS2(0,0,0)
      CSO(1,2,1) = CS1(0,2,1)*CS2(1,0,0) + CS1(1,1,1)*CS2(0,1,0)
     .           + CS1(1,2,0)*CS2(0,0,1)
      CSO(1,3,0) = CS1(0,3,0)*CS2(1,0,0) + CS1(1,2,0)*CS2(0,1,0)
      CSO(2,0,0) = CS1(1,0,0)*CS2(1,0,0) + CS1(2,0,0)*CS2(0,0,0)
      CSO(2,0,1) = CS1(1,0,1)*CS2(1,0,0) + CS1(2,0,0)*CS2(0,0,1)
     .           + CS1(2,0,1)*CS2(0,0,0)
      CSO(2,0,2) = CS1(1,0,2)*CS2(1,0,0) + CS1(2,0,1)*CS2(0,0,1)
      CSO(2,1,0) = CS1(1,1,0)*CS2(1,0,0) + CS1(2,0,0)*CS2(0,1,0)
     .           + CS1(2,1,0)*CS2(0,0,0)
      CSO(2,1,1) = CS1(1,1,1)*CS2(1,0,0) + CS1(2,0,1)*CS2(0,1,0)
     .           + CS1(2,1,0)*CS2(0,0,1)
      CSO(2,2,0) = CS1(1,2,0)*CS2(1,0,0) + CS1(2,1,0)*CS2(0,1,0)
      CSO(3,0,0) = CS1(2,0,0)*CS2(1,0,0) + CS1(3,0,0)*CS2(0,0,0)
      CSO(3,0,1) = CS1(2,0,1)*CS2(1,0,0) + CS1(3,0,0)*CS2(0,0,1)
      CSO(3,1,0) = CS1(2,1,0)*CS2(1,0,0) + CS1(3,0,0)*CS2(0,1,0)
      CSO(4,0,0) = CS1(3,0,0)*CS2(1,0,0)
      GOTO 9999
C    NA2 = 1
 2090 CONTINUE
      DO 2490 IZ1=0,NA1
        DO 2480 IY1=0,NA1-IZ1
          DO 2470 IX1=0,NA1-(IZ1+IY1)
             C1    = CS1(IX1,IY1,IZ1)
             CSO(IX1+0,IY1+0,IZ1+0)=CSO(IX1+0,IY1+0,IZ1+0)+C1*CS2(0,0,0)
             CSO(IX1+1,IY1+0,IZ1+0)=CSO(IX1+1,IY1+0,IZ1+0)+C1*CS2(1,0,0)
             CSO(IX1+0,IY1+1,IZ1+0)=CSO(IX1+0,IY1+1,IZ1+0)+C1*CS2(0,1,0)
             CSO(IX1+0,IY1+0,IZ1+1)=CSO(IX1+0,IY1+0,IZ1+1)+C1*CS2(0,0,1)
 2470     CONTINUE
 2480   CONTINUE
 2490 CONTINUE
      GOTO 9999
C    NA1 = 0, NA2 = 2
 3000 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      CSO(0,0,1) = CS1(0,0,0)*CS2(0,0,1)
      CSO(0,0,2) = CS1(0,0,0)*CS2(0,0,2)
      CSO(0,1,0) = CS1(0,0,0)*CS2(0,1,0)
      CSO(0,1,1) = CS1(0,0,0)*CS2(0,1,1)
      CSO(0,2,0) = CS1(0,0,0)*CS2(0,2,0)
      CSO(1,0,0) = CS1(0,0,0)*CS2(1,0,0)
      CSO(1,0,1) = CS1(0,0,0)*CS2(1,0,1)
      CSO(1,1,0) = CS1(0,0,0)*CS2(1,1,0)
      CSO(2,0,0) = CS1(0,0,0)*CS2(2,0,0)
      GOTO 9999
C    NA1 = 1, NA2 = 2
 3010 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      CSO(0,0,1) = CS1(0,0,0)*CS2(0,0,1) + CS1(0,0,1)*CS2(0,0,0)
      CSO(0,0,2) = CS1(0,0,0)*CS2(0,0,2) + CS1(0,0,1)*CS2(0,0,1)
      CSO(0,0,3) = CS1(0,0,1)*CS2(0,0,2)
      CSO(0,1,0) = CS1(0,0,0)*CS2(0,1,0) + CS1(0,1,0)*CS2(0,0,0)
      CSO(0,1,1) = CS1(0,0,0)*CS2(0,1,1) + CS1(0,0,1)*CS2(0,1,0)
     .           + CS1(0,1,0)*CS2(0,0,1)
      CSO(0,1,2) = CS1(0,0,1)*CS2(0,1,1) + CS1(0,1,0)*CS2(0,0,2)
      CSO(0,2,0) = CS1(0,0,0)*CS2(0,2,0) + CS1(0,1,0)*CS2(0,1,0)
      CSO(0,2,1) = CS1(0,0,1)*CS2(0,2,0) + CS1(0,1,0)*CS2(0,1,1)
      CSO(0,3,0) = CS1(0,1,0)*CS2(0,2,0)
      CSO(1,0,0) = CS1(0,0,0)*CS2(1,0,0) + CS1(1,0,0)*CS2(0,0,0)
      CSO(1,0,1) = CS1(0,0,0)*CS2(1,0,1) + CS1(0,0,1)*CS2(1,0,0)
     .           + CS1(1,0,0)*CS2(0,0,1)
      CSO(1,0,2) = CS1(0,0,1)*CS2(1,0,1) + CS1(1,0,0)*CS2(0,0,2)
      CSO(1,1,0) = CS1(0,0,0)*CS2(1,1,0) + CS1(0,1,0)*CS2(1,0,0)
     .           + CS1(1,0,0)*CS2(0,1,0)
      CSO(1,1,1) = CS1(0,0,1)*CS2(1,1,0) + CS1(0,1,0)*CS2(1,0,1)
     .           + CS1(1,0,0)*CS2(0,1,1)
      CSO(1,2,0) = CS1(0,1,0)*CS2(1,1,0) + CS1(1,0,0)*CS2(0,2,0)
      CSO(2,0,0) = CS1(0,0,0)*CS2(2,0,0) + CS1(1,0,0)*CS2(1,0,0)
      CSO(2,0,1) = CS1(0,0,1)*CS2(2,0,0) + CS1(1,0,0)*CS2(1,0,1)
      CSO(2,1,0) = CS1(0,1,0)*CS2(2,0,0) + CS1(1,0,0)*CS2(1,1,0)
      CSO(3,0,0) = CS1(1,0,0)*CS2(2,0,0)
      GOTO 9999
C    NA1 = 2, NA2 = 2
 3020 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      CSO(0,0,1) = CS1(0,0,0)*CS2(0,0,1) + CS1(0,0,1)*CS2(0,0,0)
      CSO(0,0,2) = CS1(0,0,0)*CS2(0,0,2) + CS1(0,0,1)*CS2(0,0,1)
     .           + CS1(0,0,2)*CS2(0,0,0)
      CSO(0,0,3) = CS1(0,0,1)*CS2(0,0,2) + CS1(0,0,2)*CS2(0,0,1)
      CSO(0,0,4) = CS1(0,0,2)*CS2(0,0,2)
      CSO(0,1,0) = CS1(0,0,0)*CS2(0,1,0) + CS1(0,1,0)*CS2(0,0,0)
      CSO(0,1,1) = CS1(0,0,0)*CS2(0,1,1) + CS1(0,0,1)*CS2(0,1,0)
     .           + CS1(0,1,0)*CS2(0,0,1) + CS1(0,1,1)*CS2(0,0,0)
      CSO(0,1,2) = CS1(0,0,1)*CS2(0,1,1) + CS1(0,0,2)*CS2(0,1,0)
     .           + CS1(0,1,0)*CS2(0,0,2) + CS1(0,1,1)*CS2(0,0,1)
      CSO(0,1,3) = CS1(0,0,2)*CS2(0,1,1) + CS1(0,1,1)*CS2(0,0,2)
      CSO(0,2,0) = CS1(0,0,0)*CS2(0,2,0) + CS1(0,1,0)*CS2(0,1,0)
     .           + CS1(0,2,0)*CS2(0,0,0)
      CSO(0,2,1) = CS1(0,0,1)*CS2(0,2,0) + CS1(0,1,0)*CS2(0,1,1)
     .           + CS1(0,1,1)*CS2(0,1,0) + CS1(0,2,0)*CS2(0,0,1)
      CSO(0,2,2) = CS1(0,0,2)*CS2(0,2,0) + CS1(0,1,1)*CS2(0,1,1)
     .           + CS1(0,2,0)*CS2(0,0,2)
      CSO(0,3,0) = CS1(0,1,0)*CS2(0,2,0) + CS1(0,2,0)*CS2(0,1,0)
      CSO(0,3,1) = CS1(0,1,1)*CS2(0,2,0) + CS1(0,2,0)*CS2(0,1,1)
      CSO(0,4,0) = CS1(0,2,0)*CS2(0,2,0)
      CSO(1,0,0) = CS1(0,0,0)*CS2(1,0,0) + CS1(1,0,0)*CS2(0,0,0)
      CSO(1,0,1) = CS1(0,0,0)*CS2(1,0,1) + CS1(0,0,1)*CS2(1,0,0)
     .           + CS1(1,0,0)*CS2(0,0,1) + CS1(1,0,1)*CS2(0,0,0)
      CSO(1,0,2) = CS1(0,0,1)*CS2(1,0,1) + CS1(0,0,2)*CS2(1,0,0)
     .           + CS1(1,0,0)*CS2(0,0,2) + CS1(1,0,1)*CS2(0,0,1)
      CSO(1,0,3) = CS1(0,0,2)*CS2(1,0,1) + CS1(1,0,1)*CS2(0,0,2)
      CSO(1,1,0) = CS1(0,0,0)*CS2(1,1,0) + CS1(0,1,0)*CS2(1,0,0)
     .           + CS1(1,0,0)*CS2(0,1,0) + CS1(1,1,0)*CS2(0,0,0)
      CSO(1,1,1) = CS1(0,0,1)*CS2(1,1,0) + CS1(0,1,0)*CS2(1,0,1)
     .           + CS1(0,1,1)*CS2(1,0,0) + CS1(1,0,0)*CS2(0,1,1)
     .           + CS1(1,0,1)*CS2(0,1,0) + CS1(1,1,0)*CS2(0,0,1)
      CSO(1,1,2) = CS1(0,0,2)*CS2(1,1,0) + CS1(0,1,1)*CS2(1,0,1)
     .           + CS1(1,0,1)*CS2(0,1,1) + CS1(1,1,0)*CS2(0,0,2)
      CSO(1,2,0) = CS1(0,1,0)*CS2(1,1,0) + CS1(0,2,0)*CS2(1,0,0)
     .           + CS1(1,0,0)*CS2(0,2,0) + CS1(1,1,0)*CS2(0,1,0)
      CSO(1,2,1) = CS1(0,1,1)*CS2(1,1,0) + CS1(0,2,0)*CS2(1,0,1)
     .           + CS1(1,0,1)*CS2(0,2,0) + CS1(1,1,0)*CS2(0,1,1)
      CSO(1,3,0) = CS1(0,2,0)*CS2(1,1,0) + CS1(1,1,0)*CS2(0,2,0)
      CSO(2,0,0) = CS1(0,0,0)*CS2(2,0,0) + CS1(1,0,0)*CS2(1,0,0)
     .           + CS1(2,0,0)*CS2(0,0,0)
      CSO(2,0,1) = CS1(0,0,1)*CS2(2,0,0) + CS1(1,0,0)*CS2(1,0,1)
     .           + CS1(1,0,1)*CS2(1,0,0) + CS1(2,0,0)*CS2(0,0,1)
      CSO(2,0,2) = CS1(0,0,2)*CS2(2,0,0) + CS1(1,0,1)*CS2(1,0,1)
     .           + CS1(2,0,0)*CS2(0,0,2)
      CSO(2,1,0) = CS1(0,1,0)*CS2(2,0,0) + CS1(1,0,0)*CS2(1,1,0)
     .           + CS1(1,1,0)*CS2(1,0,0) + CS1(2,0,0)*CS2(0,1,0)
      CSO(2,1,1) = CS1(0,1,1)*CS2(2,0,0) + CS1(1,0,1)*CS2(1,1,0)
     .           + CS1(1,1,0)*CS2(1,0,1) + CS1(2,0,0)*CS2(0,1,1)
      CSO(2,2,0) = CS1(0,2,0)*CS2(2,0,0) + CS1(1,1,0)*CS2(1,1,0)
     .           + CS1(2,0,0)*CS2(0,2,0)
      CSO(3,0,0) = CS1(1,0,0)*CS2(2,0,0) + CS1(2,0,0)*CS2(1,0,0)
      CSO(3,0,1) = CS1(1,0,1)*CS2(2,0,0) + CS1(2,0,0)*CS2(1,0,1)
      CSO(3,1,0) = CS1(1,1,0)*CS2(2,0,0) + CS1(2,0,0)*CS2(1,1,0)
      CSO(4,0,0) = CS1(2,0,0)*CS2(2,0,0)
      GOTO 9999
C    NA2 = 2
 3090 CONTINUE
      DO 3490 IZ1=0,NA1
        DO 3480 IY1=0,NA1-IZ1
          DO 3470 IX1=0,NA1-(IZ1+IY1)
             C1    = CS1(IX1,IY1,IZ1)
             CSO(IX1+0,IY1+0,IZ1+0)=CSO(IX1+0,IY1+0,IZ1+0)+C1*CS2(0,0,0)
             CSO(IX1+1,IY1+0,IZ1+0)=CSO(IX1+1,IY1+0,IZ1+0)+C1*CS2(1,0,0)
             CSO(IX1+2,IY1+0,IZ1+0)=CSO(IX1+2,IY1+0,IZ1+0)+C1*CS2(2,0,0)
             CSO(IX1+0,IY1+1,IZ1+0)=CSO(IX1+0,IY1+1,IZ1+0)+C1*CS2(0,1,0)
             CSO(IX1+1,IY1+1,IZ1+0)=CSO(IX1+1,IY1+1,IZ1+0)+C1*CS2(1,1,0)
             CSO(IX1+0,IY1+2,IZ1+0)=CSO(IX1+0,IY1+2,IZ1+0)+C1*CS2(0,2,0)
             CSO(IX1+0,IY1+0,IZ1+1)=CSO(IX1+0,IY1+0,IZ1+1)+C1*CS2(0,0,1)
             CSO(IX1+1,IY1+0,IZ1+1)=CSO(IX1+1,IY1+0,IZ1+1)+C1*CS2(1,0,1)
             CSO(IX1+0,IY1+1,IZ1+1)=CSO(IX1+0,IY1+1,IZ1+1)+C1*CS2(0,1,1)
             CSO(IX1+0,IY1+0,IZ1+2)=CSO(IX1+0,IY1+0,IZ1+2)+C1*CS2(0,0,2)
 3470     CONTINUE
 3480   CONTINUE
 3490 CONTINUE
      GOTO 9999
C    NA1 = 0, NA2 = 3
 4000 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      CSO(0,0,1) = CS1(0,0,0)*CS2(0,0,1)
      CSO(0,0,2) = CS1(0,0,0)*CS2(0,0,2)
      CSO(0,0,3) = CS1(0,0,0)*CS2(0,0,3)
      CSO(0,1,0) = CS1(0,0,0)*CS2(0,1,0)
      CSO(0,1,1) = CS1(0,0,0)*CS2(0,1,1)
      CSO(0,1,2) = CS1(0,0,0)*CS2(0,1,2)
      CSO(0,2,0) = CS1(0,0,0)*CS2(0,2,0)
      CSO(0,2,1) = CS1(0,0,0)*CS2(0,2,1)
      CSO(0,3,0) = CS1(0,0,0)*CS2(0,3,0)
      CSO(1,0,0) = CS1(0,0,0)*CS2(1,0,0)
      CSO(1,0,1) = CS1(0,0,0)*CS2(1,0,1)
      CSO(1,0,2) = CS1(0,0,0)*CS2(1,0,2)
      CSO(1,1,0) = CS1(0,0,0)*CS2(1,1,0)
      CSO(1,1,1) = CS1(0,0,0)*CS2(1,1,1)
      CSO(1,2,0) = CS1(0,0,0)*CS2(1,2,0)
      CSO(2,0,0) = CS1(0,0,0)*CS2(2,0,0)
      CSO(2,0,1) = CS1(0,0,0)*CS2(2,0,1)
      CSO(2,1,0) = CS1(0,0,0)*CS2(2,1,0)
      CSO(3,0,0) = CS1(0,0,0)*CS2(3,0,0)
      GOTO 9999
C    NA1 = 1, NA2 = 3
 4010 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      CSO(0,0,1) = CS1(0,0,0)*CS2(0,0,1) + CS1(0,0,1)*CS2(0,0,0)
      CSO(0,0,2) = CS1(0,0,0)*CS2(0,0,2) + CS1(0,0,1)*CS2(0,0,1)
      CSO(0,0,3) = CS1(0,0,0)*CS2(0,0,3) + CS1(0,0,1)*CS2(0,0,2)
      CSO(0,0,4) = CS1(0,0,1)*CS2(0,0,3)
      CSO(0,1,0) = CS1(0,0,0)*CS2(0,1,0) + CS1(0,1,0)*CS2(0,0,0)
      CSO(0,1,1) = CS1(0,0,0)*CS2(0,1,1) + CS1(0,0,1)*CS2(0,1,0)
     .           + CS1(0,1,0)*CS2(0,0,1)
      CSO(0,1,2) = CS1(0,0,0)*CS2(0,1,2) + CS1(0,0,1)*CS2(0,1,1)
     .           + CS1(0,1,0)*CS2(0,0,2)
      CSO(0,1,3) = CS1(0,0,1)*CS2(0,1,2) + CS1(0,1,0)*CS2(0,0,3)
      CSO(0,2,0) = CS1(0,0,0)*CS2(0,2,0) + CS1(0,1,0)*CS2(0,1,0)
      CSO(0,2,1) = CS1(0,0,0)*CS2(0,2,1) + CS1(0,0,1)*CS2(0,2,0)
     .           + CS1(0,1,0)*CS2(0,1,1)
      CSO(0,2,2) = CS1(0,0,1)*CS2(0,2,1) + CS1(0,1,0)*CS2(0,1,2)
      CSO(0,3,0) = CS1(0,0,0)*CS2(0,3,0) + CS1(0,1,0)*CS2(0,2,0)
      CSO(0,3,1) = CS1(0,0,1)*CS2(0,3,0) + CS1(0,1,0)*CS2(0,2,1)
      CSO(0,4,0) = CS1(0,1,0)*CS2(0,3,0)
      CSO(1,0,0) = CS1(0,0,0)*CS2(1,0,0) + CS1(1,0,0)*CS2(0,0,0)
      CSO(1,0,1) = CS1(0,0,0)*CS2(1,0,1) + CS1(0,0,1)*CS2(1,0,0)
     .           + CS1(1,0,0)*CS2(0,0,1)
      CSO(1,0,2) = CS1(0,0,0)*CS2(1,0,2) + CS1(0,0,1)*CS2(1,0,1)
     .           + CS1(1,0,0)*CS2(0,0,2)
      CSO(1,0,3) = CS1(0,0,1)*CS2(1,0,2) + CS1(1,0,0)*CS2(0,0,3)
      CSO(1,1,0) = CS1(0,0,0)*CS2(1,1,0) + CS1(0,1,0)*CS2(1,0,0)
     .           + CS1(1,0,0)*CS2(0,1,0)
      CSO(1,1,1) = CS1(0,0,0)*CS2(1,1,1) + CS1(0,0,1)*CS2(1,1,0)
     .           + CS1(0,1,0)*CS2(1,0,1) + CS1(1,0,0)*CS2(0,1,1)
      CSO(1,1,2) = CS1(0,0,1)*CS2(1,1,1) + CS1(0,1,0)*CS2(1,0,2)
     .           + CS1(1,0,0)*CS2(0,1,2)
      CSO(1,2,0) = CS1(0,0,0)*CS2(1,2,0) + CS1(0,1,0)*CS2(1,1,0)
     .           + CS1(1,0,0)*CS2(0,2,0)
      CSO(1,2,1) = CS1(0,0,1)*CS2(1,2,0) + CS1(0,1,0)*CS2(1,1,1)
     .           + CS1(1,0,0)*CS2(0,2,1)
      CSO(1,3,0) = CS1(0,1,0)*CS2(1,2,0) + CS1(1,0,0)*CS2(0,3,0)
      CSO(2,0,0) = CS1(0,0,0)*CS2(2,0,0) + CS1(1,0,0)*CS2(1,0,0)
      CSO(2,0,1) = CS1(0,0,0)*CS2(2,0,1) + CS1(0,0,1)*CS2(2,0,0)
     .           + CS1(1,0,0)*CS2(1,0,1)
      CSO(2,0,2) = CS1(0,0,1)*CS2(2,0,1) + CS1(1,0,0)*CS2(1,0,2)
      CSO(2,1,0) = CS1(0,0,0)*CS2(2,1,0) + CS1(0,1,0)*CS2(2,0,0)
     .           + CS1(1,0,0)*CS2(1,1,0)
      CSO(2,1,1) = CS1(0,0,1)*CS2(2,1,0) + CS1(0,1,0)*CS2(2,0,1)
     .           + CS1(1,0,0)*CS2(1,1,1)
      CSO(2,2,0) = CS1(0,1,0)*CS2(2,1,0) + CS1(1,0,0)*CS2(1,2,0)
      CSO(3,0,0) = CS1(0,0,0)*CS2(3,0,0) + CS1(1,0,0)*CS2(2,0,0)
      CSO(3,0,1) = CS1(0,0,1)*CS2(3,0,0) + CS1(1,0,0)*CS2(2,0,1)
      CSO(3,1,0) = CS1(0,1,0)*CS2(3,0,0) + CS1(1,0,0)*CS2(2,1,0)
      CSO(4,0,0) = CS1(1,0,0)*CS2(3,0,0)
      GOTO 9999
C    NA2 = 3
 4090 CONTINUE
      DO 4490 IZ1=0,NA1
        DO 4480 IY1=0,NA1-IZ1
          DO 4470 IX1=0,NA1-(IZ1+IY1)
             C1    = CS1(IX1,IY1,IZ1)
             CSO(IX1+0,IY1+0,IZ1+0)=CSO(IX1+0,IY1+0,IZ1+0)+C1*CS2(0,0,0)
             CSO(IX1+1,IY1+0,IZ1+0)=CSO(IX1+1,IY1+0,IZ1+0)+C1*CS2(1,0,0)
             CSO(IX1+2,IY1+0,IZ1+0)=CSO(IX1+2,IY1+0,IZ1+0)+C1*CS2(2,0,0)
             CSO(IX1+3,IY1+0,IZ1+0)=CSO(IX1+3,IY1+0,IZ1+0)+C1*CS2(3,0,0)
             CSO(IX1+0,IY1+1,IZ1+0)=CSO(IX1+0,IY1+1,IZ1+0)+C1*CS2(0,1,0)
             CSO(IX1+1,IY1+1,IZ1+0)=CSO(IX1+1,IY1+1,IZ1+0)+C1*CS2(1,1,0)
             CSO(IX1+2,IY1+1,IZ1+0)=CSO(IX1+2,IY1+1,IZ1+0)+C1*CS2(2,1,0)
             CSO(IX1+0,IY1+2,IZ1+0)=CSO(IX1+0,IY1+2,IZ1+0)+C1*CS2(0,2,0)
             CSO(IX1+1,IY1+2,IZ1+0)=CSO(IX1+1,IY1+2,IZ1+0)+C1*CS2(1,2,0)
             CSO(IX1+0,IY1+3,IZ1+0)=CSO(IX1+0,IY1+3,IZ1+0)+C1*CS2(0,3,0)
             CSO(IX1+0,IY1+0,IZ1+1)=CSO(IX1+0,IY1+0,IZ1+1)+C1*CS2(0,0,1)
             CSO(IX1+1,IY1+0,IZ1+1)=CSO(IX1+1,IY1+0,IZ1+1)+C1*CS2(1,0,1)
             CSO(IX1+2,IY1+0,IZ1+1)=CSO(IX1+2,IY1+0,IZ1+1)+C1*CS2(2,0,1)
             CSO(IX1+0,IY1+1,IZ1+1)=CSO(IX1+0,IY1+1,IZ1+1)+C1*CS2(0,1,1)
             CSO(IX1+1,IY1+1,IZ1+1)=CSO(IX1+1,IY1+1,IZ1+1)+C1*CS2(1,1,1)
             CSO(IX1+0,IY1+2,IZ1+1)=CSO(IX1+0,IY1+2,IZ1+1)+C1*CS2(0,2,1)
             CSO(IX1+0,IY1+0,IZ1+2)=CSO(IX1+0,IY1+0,IZ1+2)+C1*CS2(0,0,2)
             CSO(IX1+1,IY1+0,IZ1+2)=CSO(IX1+1,IY1+0,IZ1+2)+C1*CS2(1,0,2)
             CSO(IX1+0,IY1+1,IZ1+2)=CSO(IX1+0,IY1+1,IZ1+2)+C1*CS2(0,1,2)
             CSO(IX1+0,IY1+0,IZ1+3)=CSO(IX1+0,IY1+0,IZ1+3)+C1*CS2(0,0,3)
 4470     CONTINUE
 4480   CONTINUE
 4490 CONTINUE
      GOTO 9999
C    NA1 = 0, NA2 = 4
 5000 CONTINUE
      CSO(0,0,0) = CS1(0,0,0)*CS2(0,0,0)
      CSO(0,0,1) = CS1(0,0,0)*CS2(0,0,1)
      CSO(0,0,2) = CS1(0,0,0)*CS2(0,0,2)
      CSO(0,0,3) = CS1(0,0,0)*CS2(0,0,3)
      CSO(0,0,4) = CS1(0,0,0)*CS2(0,0,4)
      CSO(0,1,0) = CS1(0,0,0)*CS2(0,1,0)
      CSO(0,1,1) = CS1(0,0,0)*CS2(0,1,1)
      CSO(0,1,2) = CS1(0,0,0)*CS2(0,1,2)
      CSO(0,1,3) = CS1(0,0,0)*CS2(0,1,3)
      CSO(0,2,0) = CS1(0,0,0)*CS2(0,2,0)
      CSO(0,2,1) = CS1(0,0,0)*CS2(0,2,1)
      CSO(0,2,2) = CS1(0,0,0)*CS2(0,2,2)
      CSO(0,3,0) = CS1(0,0,0)*CS2(0,3,0)
      CSO(0,3,1) = CS1(0,0,0)*CS2(0,3,1)
      CSO(0,4,0) = CS1(0,0,0)*CS2(0,4,0)
      CSO(1,0,0) = CS1(0,0,0)*CS2(1,0,0)
      CSO(1,0,1) = CS1(0,0,0)*CS2(1,0,1)
      CSO(1,0,2) = CS1(0,0,0)*CS2(1,0,2)
      CSO(1,0,3) = CS1(0,0,0)*CS2(1,0,3)
      CSO(1,1,0) = CS1(0,0,0)*CS2(1,1,0)
      CSO(1,1,1) = CS1(0,0,0)*CS2(1,1,1)
      CSO(1,1,2) = CS1(0,0,0)*CS2(1,1,2)
      CSO(1,2,0) = CS1(0,0,0)*CS2(1,2,0)
      CSO(1,2,1) = CS1(0,0,0)*CS2(1,2,1)
      CSO(1,3,0) = CS1(0,0,0)*CS2(1,3,0)
      CSO(2,0,0) = CS1(0,0,0)*CS2(2,0,0)
      CSO(2,0,1) = CS1(0,0,0)*CS2(2,0,1)
      CSO(2,0,2) = CS1(0,0,0)*CS2(2,0,2)
      CSO(2,1,0) = CS1(0,0,0)*CS2(2,1,0)
      CSO(2,1,1) = CS1(0,0,0)*CS2(2,1,1)
      CSO(2,2,0) = CS1(0,0,0)*CS2(2,2,0)
      CSO(3,0,0) = CS1(0,0,0)*CS2(3,0,0)
      CSO(3,0,1) = CS1(0,0,0)*CS2(3,0,1)
      CSO(3,1,0) = CS1(0,0,0)*CS2(3,1,0)
      CSO(4,0,0) = CS1(0,0,0)*CS2(4,0,0)
      GOTO 9999
C    NA2 = 4
 5090 CONTINUE
      DO 5490 IZ1=0,NA1
        DO 5480 IY1=0,NA1-IZ1
          DO 5470 IX1=0,NA1-(IZ1+IY1)
             C1    = CS1(IX1,IY1,IZ1)
             CSO(IX1+0,IY1+0,IZ1+0)=CSO(IX1+0,IY1+0,IZ1+0)+C1*CS2(0,0,0)
             CSO(IX1+1,IY1+0,IZ1+0)=CSO(IX1+1,IY1+0,IZ1+0)+C1*CS2(1,0,0)
             CSO(IX1+2,IY1+0,IZ1+0)=CSO(IX1+2,IY1+0,IZ1+0)+C1*CS2(2,0,0)
             CSO(IX1+3,IY1+0,IZ1+0)=CSO(IX1+3,IY1+0,IZ1+0)+C1*CS2(3,0,0)
             CSO(IX1+4,IY1+0,IZ1+0)=CSO(IX1+4,IY1+0,IZ1+0)+C1*CS2(4,0,0)
             CSO(IX1+0,IY1+1,IZ1+0)=CSO(IX1+0,IY1+1,IZ1+0)+C1*CS2(0,1,0)
             CSO(IX1+1,IY1+1,IZ1+0)=CSO(IX1+1,IY1+1,IZ1+0)+C1*CS2(1,1,0)
             CSO(IX1+2,IY1+1,IZ1+0)=CSO(IX1+2,IY1+1,IZ1+0)+C1*CS2(2,1,0)
             CSO(IX1+3,IY1+1,IZ1+0)=CSO(IX1+3,IY1+1,IZ1+0)+C1*CS2(3,1,0)
             CSO(IX1+0,IY1+2,IZ1+0)=CSO(IX1+0,IY1+2,IZ1+0)+C1*CS2(0,2,0)
             CSO(IX1+1,IY1+2,IZ1+0)=CSO(IX1+1,IY1+2,IZ1+0)+C1*CS2(1,2,0)
             CSO(IX1+2,IY1+2,IZ1+0)=CSO(IX1+2,IY1+2,IZ1+0)+C1*CS2(2,2,0)
             CSO(IX1+0,IY1+3,IZ1+0)=CSO(IX1+0,IY1+3,IZ1+0)+C1*CS2(0,3,0)
             CSO(IX1+1,IY1+3,IZ1+0)=CSO(IX1+1,IY1+3,IZ1+0)+C1*CS2(1,3,0)
             CSO(IX1+0,IY1+4,IZ1+0)=CSO(IX1+0,IY1+4,IZ1+0)+C1*CS2(0,4,0)
             CSO(IX1+0,IY1+0,IZ1+1)=CSO(IX1+0,IY1+0,IZ1+1)+C1*CS2(0,0,1)
             CSO(IX1+1,IY1+0,IZ1+1)=CSO(IX1+1,IY1+0,IZ1+1)+C1*CS2(1,0,1)
             CSO(IX1+2,IY1+0,IZ1+1)=CSO(IX1+2,IY1+0,IZ1+1)+C1*CS2(2,0,1)
             CSO(IX1+3,IY1+0,IZ1+1)=CSO(IX1+3,IY1+0,IZ1+1)+C1*CS2(3,0,1)
             CSO(IX1+0,IY1+1,IZ1+1)=CSO(IX1+0,IY1+1,IZ1+1)+C1*CS2(0,1,1)
             CSO(IX1+1,IY1+1,IZ1+1)=CSO(IX1+1,IY1+1,IZ1+1)+C1*CS2(1,1,1)
             CSO(IX1+2,IY1+1,IZ1+1)=CSO(IX1+2,IY1+1,IZ1+1)+C1*CS2(2,1,1)
             CSO(IX1+0,IY1+2,IZ1+1)=CSO(IX1+0,IY1+2,IZ1+1)+C1*CS2(0,2,1)
             CSO(IX1+1,IY1+2,IZ1+1)=CSO(IX1+1,IY1+2,IZ1+1)+C1*CS2(1,2,1)
             CSO(IX1+0,IY1+3,IZ1+1)=CSO(IX1+0,IY1+3,IZ1+1)+C1*CS2(0,3,1)
             CSO(IX1+0,IY1+0,IZ1+2)=CSO(IX1+0,IY1+0,IZ1+2)+C1*CS2(0,0,2)
             CSO(IX1+1,IY1+0,IZ1+2)=CSO(IX1+1,IY1+0,IZ1+2)+C1*CS2(1,0,2)
             CSO(IX1+2,IY1+0,IZ1+2)=CSO(IX1+2,IY1+0,IZ1+2)+C1*CS2(2,0,2)
             CSO(IX1+0,IY1+1,IZ1+2)=CSO(IX1+0,IY1+1,IZ1+2)+C1*CS2(0,1,2)
             CSO(IX1+1,IY1+1,IZ1+2)=CSO(IX1+1,IY1+1,IZ1+2)+C1*CS2(1,1,2)
             CSO(IX1+0,IY1+2,IZ1+2)=CSO(IX1+0,IY1+2,IZ1+2)+C1*CS2(0,2,2)
             CSO(IX1+0,IY1+0,IZ1+3)=CSO(IX1+0,IY1+0,IZ1+3)+C1*CS2(0,0,3)
             CSO(IX1+1,IY1+0,IZ1+3)=CSO(IX1+1,IY1+0,IZ1+3)+C1*CS2(1,0,3)
             CSO(IX1+0,IY1+1,IZ1+3)=CSO(IX1+0,IY1+1,IZ1+3)+C1*CS2(0,1,3)
             CSO(IX1+0,IY1+0,IZ1+4)=CSO(IX1+0,IY1+0,IZ1+4)+C1*CS2(0,0,4)
 5470     CONTINUE
 5480   CONTINUE
 5490 CONTINUE
      GOTO 9999
C    Generic loop
 9000 CONTINUE
      DO 9490 IZ1=0,NA1
        DO 9480 IY1=0,NA1-IZ1
          DO 9470 IX1=0,NA1-(IZ1+IY1)
              C1 = CS1(IX1,IY1,IZ1)
              DO 9390 IZ2=0,NA2
                DO 9380 IY2=0,NA2-IZ2
                  DO 9370 IX2=0,NA2-(IZ2+IY2)
                      C2  = CS2(IX2,IY2,IZ2)
                      IXO = IX1 + IX2
                      IYO = IY1 + IY2
                      IZO = IZ1 + IZ2
                      CSO(IXO,IYO,IZO) = CSO(IXO,IYO,IZO) + C1*C2
 9370             CONTINUE
 9380           CONTINUE
 9390         CONTINUE
 9470     CONTINUE
 9480   CONTINUE
 9490 CONTINUE
      GOTO 9999
C
 9999 CONTINUE
      RETURN
11010 FORMAT(' NAO = ',I4,' IS TOO BIG IN PSOGPR. LIMIT IS ',I4)
      END
C
      SUBROUTINE PSOGRV(NOI,NRI,NAI,CSI,ALPI,MAXAI,MAXRI,
     .                  NOO,NRO,NAO,CSO,ALPO,MAXAO,MAXRO)
C
C   Computes vector product of R operator with three-component 
C   orbitals.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      *I     - Input orbitals nest. Number of input orbitals
C               should be divisible by 3.
C      *O     - Output orbitals nest.
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Input orbitals should be in the order (IORB,IX).
C      Output is prouced in the same order.
C
C   Possible optimizations:
C
C      Compiled-in array dimensions and cloning should improve
C      performance.
C
C   Bugs:
C
C      Scaling of the algorithm we use here is frightful.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (EPS=1.D-14)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C    Parameters
      DIMENSION CSI(0:MAXAI,0:MAXAI,0:MAXAI,MAXRI,NOI/3,3), ALPI(NRI)
      DIMENSION CSO(0:MAXAO,0:MAXAO,0:MAXAO,MAXRO,NOI/3,3), ALPO(MAXRO)
C    Scalar output parameters and size checks.
      IF( 3*(NOI/3) .NE. NOI ) THEN
          WRITE(NB6,11010) NOI
          STOP 'PSOGRV'
      ENDIF
      NOO = NOI
      NRO = NRI
      NAO = NAI + 1
      IF( NRO.GT.MAXRO .OR. NAO.GT.MAXAO ) THEN
          WRITE(NB6,11020) NRO, NAO, MAXRO, MAXAO
          STOP 'PSOGRV'
      ENDIF
C    Radial part doesn't change
      DO 100 IR=1,NRI
          ALPO(IR) = ALPI(IR)
  100 CONTINUE
C    Clear output coefficients
      CALL PSOGZR(NOO,NRO,NAO,CSO,MAXAO,MAXRO)
C
      DO 1900 IO=1,NOI/3
        DO 1800 IR=1,NRI
          DO 1700 IZ=0,NAI
            DO 1600 IY=0,NAI-IZ
              DO 1500 IX=0,NAI-(IZ+IY)
                  CX = CSI(IX,IY,IZ,IR,IO,1)
                  CY = CSI(IX,IY,IZ,IR,IO,2)
                  CZ = CSI(IX,IY,IZ,IR,IO,3)
                  CSO(IX,IY+1,IZ,IR,IO,1) = CSO(IX,IY+1,IZ,IR,IO,1) + CZ
                  CSO(IX,IY,IZ+1,IR,IO,1) = CSO(IX,IY,IZ+1,IR,IO,1) - CY
                  CSO(IX,IY,IZ+1,IR,IO,2) = CSO(IX,IY,IZ+1,IR,IO,2) + CX
                  CSO(IX+1,IY,IZ,IR,IO,2) = CSO(IX+1,IY,IZ,IR,IO,2) - CZ
                  CSO(IX+1,IY,IZ,IR,IO,3) = CSO(IX+1,IY,IZ,IR,IO,3) + CY
                  CSO(IX,IY+1,IZ,IR,IO,3) = CSO(IX,IY+1,IZ,IR,IO,3) - CX
 1500         CONTINUE
 1600       CONTINUE
 1700     CONTINUE
 1800   CONTINUE
 1900 CONTINUE
C    For some cases of three-component orbitals, highest angular moment
C    contributions are always zero, so that we can weed them out now.
C    In principle, orbitals we are working with now are *always* this
C    way, but with heffalumps one never knows...
      DO 3000 IC=1,3
        DO 2900 IO=1,NOO/3
          DO 2800 IR=1,NRO
            DO 2700 IZ=0,NAO
              DO 2600 IY=0,NAO-IZ
                  IX = NAO - (IZ+IY)
                  IF( ABS(CSO(IX,IY,IZ,IR,IO,IC)).GT.EPS ) GOTO 4000
 2600         CONTINUE
 2700       CONTINUE
 2800     CONTINUE
 2900   CONTINUE
 3000 CONTINUE
C    Highest angular components are all zero, so drop 'em
      NAO = NAO - 1
 4000 CONTINUE
C
      RETURN
11010 FORMAT(' NOI = ',I4,' IS NOT DIVISIBLE BY 3 IN PSOGRV')
11020 FORMAT(' NRO = ',I4,' OR NAO = ',I4,' IS TOO BIG FOR MAXRO = ',I4,
     .       ' MAXAO = ',I4,' IN PSOGRV')
      END
C
      SUBROUTINE PSOGGP(ALP1,R1,ALP2,R2,ALPX,RX,SCG)
C
C   Computes product of two elementary gaussians.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ALP1   - Orbital exponent of the first gaussian
C      R1     - Coordinates of the first gaussian
C      ALP2   - Orbital exponent of the second gaussian
C      R2     - Coordinates of the second gaussian
C      ALPX   - (Output) exponent of the product gaussian
C      RX     - (Output) coordinates of the product gaussian
C      SCG    - (Output) overall scaling coefficient
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Possible optimizations:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      DIMENSION R1(3), R2(3), RX(3)
C
      ALPX = ALP1 + ALP2
      RA   = ONE/ALPX
      C1   = ALP1 * RA
      C2   = ALP2 * RA
      R12  = ZERO
      DO 100 I=1,3
          RX(I) = C1*R1(I) + C2*R2(I)
          R12   = R12 + (R1(I)-R2(I))**2
  100 CONTINUE
      SCG = EXP( -ALP1*ALP2*RA*R12 )
C
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION PSOGCN(NA,CS,AI)
C
C   Contracts elementary integrals with orbital expansion
C   coefficients.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NA     - Maximum order of components
C      CS     - Coefficients, three-dimensional array
C      AI     - Elementary integrals, three-dimensional array
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Possible optimizations:
C
C      Cloning for few more special cases might help. So can inlining.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXAO=12)
      PARAMETER (ZERO=0.D0)
C
      DIMENSION CS(0:MAXAO,0:MAXAO,0:MAXAO)
      DIMENSION AI(0:MAXAO,0:MAXAO,0:MAXAO)
C
      GOTO (0010,1010,2010,3010,4010,5010), NA+1
      GOTO 9010
C    NA = 0
 0010 CONTINUE
      PSOGCN = CS(0,0,0)*AI(0,0,0)
      RETURN
C    NA = 1
 1010 CONTINUE
      PSOGCN = CS(0,0,0)*AI(0,0,0) + CS(1,0,0)*AI(1,0,0) + 
     .         CS(0,1,0)*AI(0,1,0) + CS(0,0,1)*AI(0,0,1)
      RETURN
C    NA = 2
 2010 CONTINUE
      PSOGCN =                 CS(0,0,0)*AI(0,0,0) + CS(1,0,0)*AI(1,0,0) 
     . + CS(2,0,0)*AI(2,0,0) + CS(0,1,0)*AI(0,1,0) + CS(1,1,0)*AI(1,1,0)
     . + CS(0,2,0)*AI(0,2,0) + CS(0,0,1)*AI(0,0,1) + CS(1,0,1)*AI(1,0,1)
     . + CS(0,1,1)*AI(0,1,1) + CS(0,0,2)*AI(0,0,2)
      RETURN
C    NA = 3
 3010 CONTINUE
      PSOGCN =                 CS(0,0,0)*AI(0,0,0) + CS(1,0,0)*AI(1,0,0)
     . + CS(2,0,0)*AI(2,0,0) + CS(3,0,0)*AI(3,0,0) + CS(0,1,0)*AI(0,1,0)
     . + CS(1,1,0)*AI(1,1,0) + CS(2,1,0)*AI(2,1,0) + CS(0,2,0)*AI(0,2,0)
     . + CS(1,2,0)*AI(1,2,0) + CS(0,3,0)*AI(0,3,0) + CS(0,0,1)*AI(0,0,1)
     . + CS(1,0,1)*AI(1,0,1) + CS(2,0,1)*AI(2,0,1) + CS(0,1,1)*AI(0,1,1)
     . + CS(1,1,1)*AI(1,1,1) + CS(0,2,1)*AI(0,2,1) + CS(0,0,2)*AI(0,0,2)
     . + CS(1,0,2)*AI(1,0,2) + CS(0,1,2)*AI(0,1,2) + CS(0,0,3)*AI(0,0,3)
      RETURN
C    NA = 4
 4010 CONTINUE
      PSOGCN =                 CS(0,0,0)*AI(0,0,0) + CS(1,0,0)*AI(1,0,0)
     . + CS(2,0,0)*AI(2,0,0) + CS(3,0,0)*AI(3,0,0) + CS(4,0,0)*AI(4,0,0)
     . + CS(0,1,0)*AI(0,1,0) + CS(1,1,0)*AI(1,1,0) + CS(2,1,0)*AI(2,1,0)
     . + CS(3,1,0)*AI(3,1,0) + CS(0,2,0)*AI(0,2,0) + CS(1,2,0)*AI(1,2,0)
     . + CS(2,2,0)*AI(2,2,0) + CS(0,3,0)*AI(0,3,0) + CS(1,3,0)*AI(1,3,0)
     . + CS(0,4,0)*AI(0,4,0) + CS(0,0,1)*AI(0,0,1) + CS(1,0,1)*AI(1,0,1)
     . + CS(2,0,1)*AI(2,0,1) + CS(3,0,1)*AI(3,0,1) + CS(0,1,1)*AI(0,1,1)
     . + CS(1,1,1)*AI(1,1,1) + CS(2,1,1)*AI(2,1,1) + CS(0,2,1)*AI(0,2,1)
     . + CS(1,2,1)*AI(1,2,1) + CS(0,3,1)*AI(0,3,1) + CS(0,0,2)*AI(0,0,2)
     . + CS(1,0,2)*AI(1,0,2) + CS(2,0,2)*AI(2,0,2) + CS(0,1,2)*AI(0,1,2)
     . + CS(1,1,2)*AI(1,1,2) + CS(0,2,2)*AI(0,2,2) + CS(0,0,3)*AI(0,0,3)
     . + CS(1,0,3)*AI(1,0,3) + CS(0,1,3)*AI(0,1,3) + CS(0,0,4)*AI(0,0,4)
      RETURN
C    NA = 5
 5010 CONTINUE
      PSOGCN =                 CS(0,0,0)*AI(0,0,0) + CS(1,0,0)*AI(1,0,0)
     . + CS(2,0,0)*AI(2,0,0) + CS(3,0,0)*AI(3,0,0) + CS(4,0,0)*AI(4,0,0)
     . + CS(5,0,0)*AI(5,0,0) + CS(0,1,0)*AI(0,1,0) + CS(1,1,0)*AI(1,1,0)
     . + CS(2,1,0)*AI(2,1,0) + CS(3,1,0)*AI(3,1,0) + CS(4,1,0)*AI(4,1,0)
     . + CS(0,2,0)*AI(0,2,0) + CS(1,2,0)*AI(1,2,0) + CS(2,2,0)*AI(2,2,0)
     . + CS(3,2,0)*AI(3,2,0) + CS(0,3,0)*AI(0,3,0) + CS(1,3,0)*AI(1,3,0)
     . + CS(2,3,0)*AI(2,3,0) + CS(0,4,0)*AI(0,4,0) + CS(1,4,0)*AI(1,4,0)
     . + CS(0,5,0)*AI(0,5,0) + CS(0,0,1)*AI(0,0,1) + CS(1,0,1)*AI(1,0,1)
     . + CS(2,0,1)*AI(2,0,1) + CS(3,0,1)*AI(3,0,1) + CS(4,0,1)*AI(4,0,1)
     . + CS(0,1,1)*AI(0,1,1) + CS(1,1,1)*AI(1,1,1) + CS(2,1,1)*AI(2,1,1)
     . + CS(3,1,1)*AI(3,1,1) + CS(0,2,1)*AI(0,2,1) + CS(1,2,1)*AI(1,2,1)
     . + CS(2,2,1)*AI(2,2,1) + CS(0,3,1)*AI(0,3,1) + CS(1,3,1)*AI(1,3,1)
     . + CS(0,4,1)*AI(0,4,1) + CS(0,0,2)*AI(0,0,2) + CS(1,0,2)*AI(1,0,2)
     . + CS(2,0,2)*AI(2,0,2) + CS(3,0,2)*AI(3,0,2) + CS(0,1,2)*AI(0,1,2)
     . + CS(1,1,2)*AI(1,1,2) + CS(2,1,2)*AI(2,1,2) + CS(0,2,2)*AI(0,2,2)
     . + CS(1,2,2)*AI(1,2,2) + CS(0,3,2)*AI(0,3,2) + CS(0,0,3)*AI(0,0,3)
     . + CS(1,0,3)*AI(1,0,3) + CS(2,0,3)*AI(2,0,3) + CS(0,1,3)*AI(0,1,3)
     . + CS(1,1,3)*AI(1,1,3) + CS(0,2,3)*AI(0,2,3) + CS(0,0,4)*AI(0,0,4)
     . + CS(1,0,4)*AI(1,0,4) + CS(0,1,4)*AI(0,1,4) + CS(0,0,5)*AI(0,0,5)
      RETURN
C    Common case
 9010 CONTINUE
      PSOGCN = ZERO
      DO 9190 IZ=0,NA
        DO 9180 IY=0,NA-IZ
          DO 9170 IX=0,NA-(IZ+IY)
              PSOGCN = PSOGCN + CS(IX,IY,IZ)*AI(IX,IY,IZ)
 9170     CONTINUE
 9180   CONTINUE
 9190 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSGLR3(IGO,N1,L1,Z1,R1,N2,L2,Z2,R2,AINT,LDA)
C
C   Compute three-center matrix elements of the operator r^-3 L
C   using Gaussian expansion. See also PS3LR3.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IGO    - Gaussian expansion order to use.
C      N1     - Major quantum number of the STO block on left-hand center
C      L1     - Orbital quantum number of the STO block on left-hand center
C      Z1     - Orbital exponent on left-hand center
C      R1     - Array of X,Y,Z coordinates of left-hand center relative
C               to the operator center.
C      N2,L2,Z2,R2
C             - Attributes of the STO block on right-hand center
C      AINT   - (Output) integrals over STOs, three-index array
C               First index: 1 to 2*L1+1, magnetic quantum number of left-hand
C                            orbitals in the increasing order.
C               Second index: 1 to 3, projections of the orbital moment,
C                            in the X,Y,Z order.
C               Third index: 1 to 2*L2+1, magnetic quantum number of right-hand
C                            orbitals in the increasing order.
C      LDA    - Leading dimension of the AINT matrix.
C
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C      Basically, 12*(MAXAI+1)**3*MAXGO (2*MAXL+1) + 2*(MAXAO+1)**3 DOUBLE
C      PRECISION words. With MAXL=4,MAXGO=5,MAXAI=6 and MAXAO=12, this
C      amounts to some 190,000 DP words, or 1.5Mb. Quite a bit of this
C      can be shaved off by reusing coefficient arrays.
C
C   Module logic:
C
C      We are following a literal reduction path to the elementary integrals
C      given by PSOGI. There is nothing smart or complicated in here.
C
C   Precision:
C
C      With respect to gaussian basis, results should be good essentially
C      to machine precision, provided that orbital moment is sufficiently
C      small (<=3). Compared to exact integrals over Slater orbitals, accuracy
C      is determined by the STO-nG expansion used. With STO-5G and in 
C      "chemical" region (all intercenter distances from 1 to 6 Angstrom,
C      orbital exponents from 0.9 to 4.2 au^-1), deviations from exact STO
C      integrals are at most few percent of the largest integral in the nest,
C      and usually few tenths of a percent. At larger distances and with
C      higher exponents, STO-nG integrals fall off much more rapidly than STO
C      integrals do. At small distances and with smaller exponents, accuracy
C      is satisfactory for Slater parameters without gradient discontinuity 
C      at origin, but is absolutely erratic for 1S, 2S, 2P, 3P, 3D, 4D, etc.
C      orbitals. Accuracy of the results compared to exact STO integrals
C      gets progressively worse for smaller-order expansions. However, if
C      STO-3G expansions are used for 1S, 2S, 2P, 3S and 3P orbitals in "chemical"
C      region, 90% of integral values are within 10% of STO-5G results, so
C      STO-3G should be adequate for the core MNDO elements. Smaller-order
C      expansions (STO-1G and STO-2G) should never be used except in testing.
C
C   Possible optimizations:
C
C      Breakdown of the execution time depends very strongly on the parameters
C      values. For small orbital quantum numbers (L<=1), execution time with STO-5G
C      expansions on R4000 is dominated by EXP (26%), PSVXBG (15%), PSOGJ (15%) and 
C      PSOGPR (12%). All other routines contribute less than 4% each. For L<=2 STO-5G,
C      PSOGPR becomes a dominant contributor (34%), followed by EXP (15%), PSVXBG
C      (10%) and PSOGCN (7%). All the rest is below 4%. For still larger orbital
C      quantum numbers, PSOGPR and PSOGCN dominate absolutely due to unfavorable
C      scaling of these routines.
C
C      PSOGPR (multiplication of angular components of two expansion terms) and
C      PSOGCN (contraction of elementary integrals with expansion coefficients)
C      are already heavily tuned, so micro-optimizations should be worthless.
C      Effort spent in these routines can potentially be decreased by a factor
C      of n (from STO-nG ;-) for PSGLR3 and PSGRL3 to n**2 for PSGRR2 by using
C      the relationship between coefficients in angular components corresponding
C      to different gaussians. As long as no differential operator acts of the
C      orbital nest, ratio of corresponding coefficients from different gaussian
C      terms remains constant, so that neither products nor contactions need to be
C      evaluated. For L<=2, however, there is preciously little incentive to
C      implement this.
C
C      As it is now, STO-5G expansions for core MNDO elements (L<=1) run about
C      ten times faster than corresponding STO integrals. For L=2, where STO
C      integrals vary wildly in the execution time, speedup is between 1 and 10,
C      depending on parameters values.
C
C   Bugs:
C
C      It is not possible to adjust MAXAI and MAXAO without modifying PSOGPR
C      and PSOGCN as well.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXGO=9,MAXL=4)
      PARAMETER (MAXAI=6,MAXAO=12)
C
      PARAMETER (MAXORB=2*MAXL+1)
      PARAMETER (MXORBD=3*MAXORB)
C
      PARAMETER (ZERO=0.D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C    Parameters
      DIMENSION R1(3), R2(3), AINT(LDA,3,2*L2+1)
C    Temporaries. 
      LOGICAL ISZERO
      DIMENSION CS1 (0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MAXORB), ALP1 (MAXGO)
      DIMENSION CS1T(0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MAXORB), ALP1T(MAXGO)
      DIMENSION CS2 (0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MAXORB), ALP2 (MAXGO)
      DIMENSION CSD (0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MXORBD), ALPD (MAXGO)
      DIMENSION CSDT(0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MXORBD), ALPDT(MAXGO)
      DIMENSION CSRV(0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MXORBD), ALPRV(MAXGO)
      DIMENSION RX(3)
      DIMENSION AINTI(0:MAXAO,0:MAXAO,0:MAXAO)
      DIMENSION CSAP (0:MAXAO,0:MAXAO,0:MAXAO)
C
      EXTERNAL PSOGCN
C    Range checks
      IF( MAX(L1,L2).GT.MAXL ) THEN
          WRITE(NB6,11010) L1,L2,MAXL
          STOP 'PSGLR3'
      ENDIF
      IF( IGO.GT.MAXGO ) THEN
          WRITE(NB6,11020) IGO, MAXGO
          STOP 'PSGLR3'
      ENDIF
      IF( LDA.LT.2*L1+1 ) THEN
          WRITE(NB6,11020) LDA, 2*L1+1
          STOP 'PSGLR3'
      ENDIF
C    Expand Slater orbitals with Gaussian functions of the given order.
      CALL PSOGMK(IGO,N1,L1,Z1,NO1,NR1,NA1,CS1,ALP1,MAXAI,MAXGO)
      CALL PSOGMK(IGO,N2,L2,Z2,NO2,NR2,NA2,CS2,ALP2,MAXAI,MAXGO)
C    Evaluate gradient of the right-hand (second) orbital.
      CALL PSOGDX(NO2,NR2,NA2,CS2,ALP2,MAXAI,MAXGO,
     .            NOD,NRD,NAD,CSD,ALPD,MAXAI,MAXGO)
C    Translate radial parts of the orbitals.
      CALL PSOGTR(R1,NO1,NR1,NA1,CS1,ALP1,MAXAI,MAXGO,
     .            NO1T,NR1T,NA1T,CS1T,ALP1T,MAXAI,MAXGO)
      CALL PSOGTR(R2,NOD,NRD,NAD,CSD,ALPD,MAXAI,MAXGO,
     .            NODT,NRDT,NADT,CSDT,ALPDT,MAXAI,MAXGO)
C    Evaluate vector product of the rhs orbitals with r.
      CALL PSOGRV(NODT,NRDT,NADT,CSDT,ALPDT,MAXAI,MAXGO,
     .            NORV,NRRV,NARV,CSRV,ALPRV,MAXAI,MAXGO)
C    Clear output array
      DO 1900 IO2=1,NO2
        DO 1800 ILX=1,3
          DO 1700 IO1=1,NO1
              AINT(IO1,ILX,IO2) = ZERO
 1700     CONTINUE
 1800   CONTINUE
 1900 CONTINUE
C    Loop over all gaussian components on both sides.
      MAXO = NA1T + NARV
      DO 8000 IR2=1,NRRV
        DO 7000 IR1=1,NR1T
C          Compute product of gaussians 
            CALL PSOGGP(ALP1T(IR1),R1,ALPRV(IR2),R2,ALPX,RX,SCG)
C          Evaluate elementary integrals with the product gaussian
            CALL PSOGI(MAXO,ALPX,RX,AINTI,MAXAO+1,ISZERO)
C          If all elementary integrals are very small, there is no
C          point in summing 'em up
            IF(ISZERO) GOTO 6999
C          Loop over all orbitals
            DO 6000 IORV=1,NORV
C              Compute output orbital indices
                ILX = (IORV-1)/NO2 + 1
                IO2 = MOD(IORV-1,NO2) + 1
                DO 5000 IO1=1,NO1T
C                  Compute product of the angular parts
                    CALL PSOGPR(NA1T,CS1T(0,0,0,IR1,IO1),NARV,
     .                      CSRV(0,0,0,IR2,IORV),MAXOX,CSAP)
C                  Evaluate integral
                    VAL = SCG*PSOGCN(MAXO,AINTI,CSAP)
                    AINT(IO1,ILX,IO2) = AINT(IO1,ILX,IO2) + VAL
 5000           CONTINUE
 6000       CONTINUE
 6999       CONTINUE
 7000   CONTINUE
 8000 CONTINUE
C
      RETURN
11010 FORMAT(' L1 (',I4,') OR L2 (',I4,') IS TOO LARGE IN PSGLR3.',
     .       ' LIMIT IS ',I4)
11020 FORMAT(' GAUSSIAN EXPANSION ORDER (',I4,') IS TOO BIG IN PSGLR3.',
     .       ' LIMIT IS ',I4)
11030 FORMAT(' LDA (',I4,') IS TOO SMALL IN PSGLR3. SHOULD BE AT LEAST '
     .       ,I4)
      END
C
      SUBROUTINE PSGRL3(IGO,N1,L1,Z1,R1,N2,L2,Z2,R2,AINT,LDA)
C
C   Compute three-center matrix elements of the operator r^a r^-3 L
C   using Gaussian expansion. See also PS3RL3.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IGO    - Gaussian expansion order to use.
C      N1     - Major quantum number of the STO block on left-hand center
C      L1     - Orbital quantum number of the STO block on left-hand center
C      Z1     - Orbital exponent on left-hand center
C      R1     - Array of X,Y,Z coordinates of left-hand center relative
C               to the operator center.
C      N2,L2,Z2,R2
C             - Attributes of the STO block on right-hand center
C      AINT   - (Output) integrals over STOs, three-index array
C               First index: 1 to 2*L1+1, magnetic quantum number of left-hand
C                            orbitals in the increasing order.
C               Second index: 1 to 3, components of the radius-vector at the
C                            left-hand center, X,Y,Z order.
C               Third index: 1 to 3, projections of the orbital moment,
C                            in the X,Y,Z order.
C               Fourth index: 1 to 2*L2+1, magnetic quantum number of right-hand
C                            orbitals in the increasing order.
C      LDA    - Leading dimension of the AINT matrix.
C
C   Accessed common blocks:
C
C   Local storage:
C
C      Basically, 17*(MAXAI+1)**3*MAXGO*(2*MAXL+1) + 2*(MAXAO+1)**3 DOUBLE
C      PRECIUSION words. With MAXAI=6,MAXAO=12,MAXGO=5 and MAXL=4, it boils 
C      down to about 267,000 DP words, or some 2Mb. Most of those can be
C      eliminated by reusing arrays for expansion coefficients.
C
C   Module logic:
C
C      This is mostly a replica of PSGLR3. See it for general discussion of
C      accuracy, bugs and optimizations.
C
C   Precision:
C
C   Possible optimizations:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXGO=9,MAXL=4)
      PARAMETER (MAXAI=6,MAXAO=12)
C
      PARAMETER (MAXORB=2*MAXL+1)
      PARAMETER (MXORBD=3*MAXORB)
C
      PARAMETER (ZERO=0.D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C    Parameters
      DIMENSION R1(3), R2(3), AINT(LDA,3,3,2*L2+1)
C    Temporaries. 
      LOGICAL ISZERO
      DIMENSION CS1 (0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MAXORB), ALP1 (MAXGO)
      DIMENSION CS1R(0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MXORBD), ALP1R(MAXGO)
      DIMENSION CS1T(0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MXORBD), ALP1T(MAXGO)
      DIMENSION CS2 (0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MAXORB), ALP2 (MAXGO)
      DIMENSION CSD (0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MXORBD), ALPD (MAXGO)
      DIMENSION CSDT(0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MXORBD), ALPDT(MAXGO)
      DIMENSION CSRV(0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MXORBD), ALPRV(MAXGO)
      DIMENSION RX(3)
      DIMENSION AINTI(0:MAXAO,0:MAXAO,0:MAXAO)
      DIMENSION CSAP (0:MAXAO,0:MAXAO,0:MAXAO)
C
      EXTERNAL PSOGCN
C    Range checks
      IF( MAX(L1,L2).GT.MAXL ) THEN
          WRITE(NB6,11010) L1,L2,MAXL
          STOP 'PSGRL3'
      ENDIF
      IF( IGO.GT.MAXGO ) THEN
          WRITE(NB6,11020) IGO, MAXGO
          STOP 'PSGRL3'
      ENDIF
      IF( LDA.LT.2*L1+1 ) THEN
          WRITE(NB6,11020) LDA, 2*L1+1
          STOP 'PSGRL3'
      ENDIF
C    Expand Slater orbitals with Gaussian functions of the given order.
      CALL PSOGMK(IGO,N1,L1,Z1,NO1,NR1,NA1,CS1,ALP1,MAXAI,MAXGO)
      CALL PSOGMK(IGO,N2,L2,Z2,NO2,NR2,NA2,CS2,ALP2,MAXAI,MAXGO)
C    Multiply left-hand orbital by radius-vector
      CALL PSOGRS(NO1,NR1,NA1,CS1,ALP1,MAXAI,MAXGO,
     .            NO1R,NR1R,NA1R,CS1R,ALP1R,MAXAI,MAXGO)
C    Evaluate gradient of the right-hand (second) orbital.
      CALL PSOGDX(NO2,NR2,NA2,CS2,ALP2,MAXAI,MAXGO,
     .            NOD,NRD,NAD,CSD,ALPD,MAXAI,MAXGO)
C    Translate radial parts of the orbitals.
      CALL PSOGTR(R1,NO1R,NR1R,NA1R,CS1R,ALP1R,MAXAI,MAXGO,
     .            NO1T,NR1T,NA1T,CS1T,ALP1T,MAXAI,MAXGO)
      CALL PSOGTR(R2,NOD,NRD,NAD,CSD,ALPD,MAXAI,MAXGO,
     .            NODT,NRDT,NADT,CSDT,ALPDT,MAXAI,MAXGO)
C    Evaluate vector product of the rhs orbitals with r.
      CALL PSOGRV(NODT,NRDT,NADT,CSDT,ALPDT,MAXAI,MAXGO,
     .            NORV,NRRV,NARV,CSRV,ALPRV,MAXAI,MAXGO)
C    Clear output array
      DO 1900 IO2=1,NO2
        DO 1800 ILX=1,3
          DO 1700 IRX=1,3
              DO 1600 IO1=1,NO1
                  AINT(IO1,IRX,ILX,IO2) = ZERO
 1600         CONTINUE
 1700     CONTINUE
 1800   CONTINUE
 1900 CONTINUE
C    Loop over all gaussian components on both sides.
      MAXO = NA1T + NARV
      DO 8000 IR2=1,NRRV
        DO 7000 IR1=1,NR1T
C          Compute product of gaussians 
            CALL PSOGGP(ALP1T(IR1),R1,ALPRV(IR2),R2,ALPX,RX,SCG)
C          Evaluate elementary integrals with the product gaussian
            CALL PSOGI(MAXO,ALPX,RX,AINTI,MAXAO+1,ISZERO)
C          If elementary integrals are zero, skip this GTO pair
            IF(ISZERO) GOTO 6999
C          Loop over all orbitals
            DO 6000 IORV=1,NORV
C              Compute output orbital indices
                ILX = (IORV-1)/NO2 + 1
                IO2 = MOD(IORV-1,NO2) + 1
                DO 5000 IO1T=1,NO1T
                   IRX = (IO1T-1)/NO1 + 1
                   IO1 = MOD(IO1T-1,NO1) + 1
C                  Compute product of the angular parts
                    CALL PSOGPR(NA1T,CS1T(0,0,0,IR1,IO1T),NARV,
     .                      CSRV(0,0,0,IR2,IORV),MAXOX,CSAP)
C                  Evaluate integral
                    VAL = SCG*PSOGCN(MAXO,AINTI,CSAP)
                    AINT(IO1,IRX,ILX,IO2) = AINT(IO1,IRX,ILX,IO2) + VAL
 5000           CONTINUE
 6000       CONTINUE
 6999       CONTINUE
 7000   CONTINUE
 8000 CONTINUE
C
      RETURN
11010 FORMAT(' L1 (',I4,') OR L2 (',I4,') IS TOO LARGE IN PSGRL3.',
     .       ' LIMIT IS ',I4)
11020 FORMAT(' GAUSSIAN EXPANSION ORDER (',I4,') IS TOO BIG IN PSGRL3.',
     .       ' LIMIT IS ',I4)
11030 FORMAT(' LDA (',I4,') IS TOO SMALL IN PSGRL3. SHOULD BE AT LEAST '
     .       ,I4)
      END
C
      SUBROUTINE PSGRR3(IGO,N1,L1,Z1,R1,N2,L2,Z2,R2,AINT,LDA)
C
C   Compute three-center matrix elements of the operator {x,y,z} r^-3 r^c
C   using Gaussian expansion. See also PS3RR3.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IGO    - Gaussian expansion order to use.
C      N1     - Major quantum number of the STO block on left-hand center
C      L1     - Orbital quantum number of the STO block on left-hand center
C      Z1     - Orbital exponent on left-hand center
C      R1     - Array of X,Y,Z coordinates of left-hand center relative
C               to the operator center.
C      N2,L2,Z2,R2
C             - Attributes of the STO block on right-hand center
C      AINT   - Output integrals, four-index array.
C               First index, 1 to 2*L1+1: Orbitals at the left-hand center,
C                   order of increasing magnetic quantum number.
C               Second index, 1 to 3: Components of the dipole moment at the
C                   operator center, XYZ order.
C               Third index, 1 to 3: Components of the dipole moment at the
C                   right-hand center, XYZ order.
C               Forth index, 1 to 2*L2+1: Orbitals at the right-hand center,
C                   order of increasing magnetic quantum number.
C      LDA    - Leading dimension of the AINT matrix.
C
C   Accessed common blocks:
C
C   Local storage:
C
C      Basically, 12*(MAXAI+1)**3*MAXGO (2*MAXL+1) + 2*(MAXAO+1)**3 DOUBLE
C      PRECISION words. With MAXL=4,MAXGO=5,MAXAI=6 and MAXAO=12, this
C      amounts to some 190,000 DP words, or 1.5Mb. Quite a bit of this
C      can be shaved off by reusing coefficient arrays.
C
C   Module logic:
C
C      This is mostly a replica of PSGLR3. See it for general discussion of
C      accuracy, bugs and optimizations.
C
C   Precision:
C
C   Possible optimizations:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXGO=9,MAXL=4)
      PARAMETER (MAXAI=6,MAXAO=12)
C
      PARAMETER (MAXORB=2*MAXL+1)
      PARAMETER (MXORBD=3*MAXORB)
C
      PARAMETER (ZERO=0.D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C    Parameters
      DIMENSION R1(3), R2(3), AINT(LDA,3,3,2*L2+1)
C    Temporaries. 
      LOGICAL ISZERO
      DIMENSION CS1 (0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MAXORB), ALP1 (MAXGO)
      DIMENSION CS1T(0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MAXORB), ALP1T(MAXGO)
      DIMENSION CS1R(0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MXORBD), ALP1R(MAXGO)
      DIMENSION CS2 (0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MAXORB), ALP2 (MAXGO)
      DIMENSION CS2R(0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MXORBD), ALP2R(MAXGO)
      DIMENSION CS2T(0:MAXAI,0:MAXAI,0:MAXAI,MAXGO,MXORBD), ALP2T(MAXGO)
      DIMENSION RX(3)
      DIMENSION AINTI(0:MAXAO,0:MAXAO,0:MAXAO)
      DIMENSION CSAP (0:MAXAO,0:MAXAO,0:MAXAO)
C
      EXTERNAL PSOGCN
C    Range checks
      IF( MAX(L1,L2).GT.MAXL ) THEN
          WRITE(NB6,11010) L1,L2,MAXL
          STOP 'PSGRR3'
      ENDIF
      IF( IGO.GT.MAXGO ) THEN
          WRITE(NB6,11020) IGO, MAXGO
          STOP 'PSGRR3'
      ENDIF
      IF( LDA.LT.2*L1+1 ) THEN
          WRITE(NB6,11020) LDA, 2*L1+1
          STOP 'PSGRL3'
      ENDIF
C    Expand Slater orbitals with Gaussian functions of the given order.
      CALL PSOGMK(IGO,N1,L1,Z1,NO1,NR1,NA1,CS1,ALP1,MAXAI,MAXGO)
      CALL PSOGMK(IGO,N2,L2,Z2,NO2,NR2,NA2,CS2,ALP2,MAXAI,MAXGO)
C    Multiply right-hand orbital by radius-vector
      CALL PSOGRS(NO2,NR2,NA2,CS2,ALP2,MAXAI,MAXGO,
     .            NO2R,NR2R,NA2R,CS2R,ALP2R,MAXAI,MAXGO)
C    Translate radial parts of the orbitals.
      CALL PSOGTR(R1,NO1,NR1,NA1,CS1,ALP1,MAXAI,MAXGO,
     .            NO1T,NR1T,NA1T,CS1T,ALP1T,MAXAI,MAXGO)
      CALL PSOGTR(R2,NO2R,NR2R,NA2R,CS2R,ALP2R,MAXAI,MAXGO,
     .            NO2T,NR2T,NA2T,CS2T,ALP2T,MAXAI,MAXGO)
C    Multiply translated left-hand orbital by radius-vector
C    (right-hand one would do as well, but we are trying to strike
C    balance between orbital orders on both sides)
      CALL PSOGRS(NO1T,NR1T,NA1T,CS1T,ALP1T,MAXAI,MAXGO,
     .            NO1R,NR1R,NA1R,CS1R,ALP1R,MAXAI,MAXGO)
C    Clear output array
      DO 1900 IO2=1,NO2
        DO 1800 IRX=1,3
          DO 1700 IRO=1,3
              DO 1600 IO1=1,NO1
                  AINT(IO1,IRO,IRX,IO2) = ZERO
 1600         CONTINUE
 1700     CONTINUE
 1800   CONTINUE
 1900 CONTINUE
C    Loop over all gaussian components on both sides.
      MAXO = NA1R + NA2T
      DO 8000 IR2=1,NR2T
        DO 7000 IR1=1,NR1R
C          Compute product of gaussians 
            CALL PSOGGP(ALP1R(IR1),R1,ALP2T(IR2),R2,ALPX,RX,SCG)
C          Evaluate elementary integrals with the product gaussian
            CALL PSOGI(MAXO,ALPX,RX,AINTI,MAXAO+1,ISZERO)
C          If elementary integrals are zero, skip this GTO pair
            IF(ISZERO) GOTO 6999
C          Loop over all orbitals
            DO 6000 IO2T=1,NO2T
C              Compute output orbital indices
                IRX = (IO2T-1)/NO2 + 1
                IO2 = MOD(IO2T-1,NO2) + 1
                DO 5000 IO1R=1,NO1R
                   IRO = (IO1R-1)/NO1 + 1
                   IO1 = MOD(IO1R-1,NO1) + 1
C                  Compute product of the angular parts
                    CALL PSOGPR(NA1R,CS1R(0,0,0,IR1,IO1R),NA2T,
     .                      CS2T(0,0,0,IR2,IO2T),MAXOX,CSAP)
C                  Evaluate integral
                    VAL = SCG*PSOGCN(MAXO,AINTI,CSAP)
                    AINT(IO1,IRO,IRX,IO2) = AINT(IO1,IRO,IRX,IO2) + VAL
 5000           CONTINUE
 6000       CONTINUE
 6999       CONTINUE
 7000   CONTINUE
 8000 CONTINUE
C
      RETURN
11010 FORMAT(' L1 (',I4,') OR L2 (',I4,') IS TOO LARGE IN PSGRL3.',
     .       ' LIMIT IS ',I4)
11020 FORMAT(' GAUSSIAN EXPANSION ORDER (',I4,') IS TOO BIG IN PSGRL3.',
     .       ' LIMIT IS ',I4)
11030 FORMAT(' LDA (',I4,') IS TOO SMALL IN PSGRL3. SHOULD BE AT LEAST '
     .       ,I4)
      END
