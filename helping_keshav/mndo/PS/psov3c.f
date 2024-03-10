C     ******************************************************************
C
C     NMR 3-center integrals over Slater AOs: Top-level routines
C     for elementary 3-center NMR integrals using Slater AOs.
C
C     ******************************************************************
C
C     Incomplete gamma function expansion of Slater AOs at displaced
C     centers and implementation of spin-orbit and spin-dipolar
C     integrals.
C
C     Computation of various 3-center integrals using modified Sharma
C     expansion at the displaced center (Phys. Rev. A, 13, 517 (1976)).
C
C     This code may call low-level routines in any other source files.
C
C     User-callable routines include:
C     PS3LR3  - Computes three-center integrals over r^-3 L
C     PS3PR3  - Computes three-center integrals over {x,y,z} r^-3
C
      SUBROUTINE PS3GME(NMAX,Y,EN)
C
C   Compute batch of values for the function 1-Exp[-y]*Sum[y^j/j!,{j,0,n}]
C   needed in the evaluation of the generalized incomplete gamma function
C   of integer order (also see PSGMI). Computed values are expected to
C   have relative accuracy of EPS, or ca. one digit less than machine
C   accuracy, whichever is less.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NMAX   - Maximum order if EN required (minimum is always zero)
C      Y      - Argument value
C      EN     - Array for the function values, at least NMAX+1 elements
C
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C   Module logic:
C
C      If EN(N,Y) is larger than threshold, direct expression is used.
C      Otherwise (i.e. if cancellation of terms is severe), EN is
C      computed as a sum of the infinite continued fraction:
C                                      1
C      Exp[-y] y^(n+1)/n! ----------------------------, where
C                                        a1
C                           b0 + --------------------
C                                           a2
C                                 b1 + ------------
C                                       b2 + ....
C
C      bi = n+1+i
C      ai = -(n+1+(i-1)/2) y   if i odd
C      ai = (i/2) y            if i even
C
C      Both continuous fraction and summation method come from
C      Abramowitz and Stegun, poketbook edition.
C
C   Possible optimizations:
C
C   Bugs:
C
C      FP range overflow may occur in the continous fraction portion
C      for some valid arguments. Rescaling continous fraction should
C      take care of that at the cost of lower performance.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (THRS=0.1D0)
      PARAMETER (EPS=2.D-16)
C
      DIMENSION EN(0:NMAX)
C
C    First, let's try to use direct expression for as much terms as
C    we can.
C
      EXPY   = EXP(-Y)
      SUM    = ONE
      TERM   = ONE
      DO 1000 N=0,NMAX
          EN(N)  = ONE - EXPY * SUM
          IF( EN(N).LT.THRS ) GOTO 2000
          TERM   = TERM*Y/(N+1)
          SUM    = SUM + TERM
 1000 CONTINUE
      RETURN
 2000 CONTINUE
C
C    Do continous fraction expansion for one argument at a time,
C    temporarily store continous fractions in EN
C
      NMIN = N
      DO 3000 N=NMIN,NMAX
          AM2   = ONE
          AM1   = N+1
          BM2   = ZERO
          BM1   = ONE
          BI    = N+1
          AIE   = ZERO
          AIO   = -N*Y
C        In principle, we are after the AM?/BM? ratio - as soon as
C        it stabilizes, we'd got that we want. However, division in
C        the loop is expensive, so that we'll rewrite the conditional
C        to use only multiplications. We will also unroll the loop
C        once, and evaluate conditional once per two cycles...
 2200     CONTINUE
C            I is always odd on this iteration
              BI  = BI + ONE
              AIO = AIO - Y
              AM0 = BI*AM1 + AIO*AM2
              BM0 = BI*BM1 + AIO*BM2
              AM2 = AM1
              AM1 = AM0
              BM2 = BM1
              BM1 = BM0
C            I is always even on this iteration
              BI  = BI + ONE
              AIE = AIE + Y
              AM0 = BI*AM1 + AIE*AM2
              BM0 = BI*BM1 + AIE*BM2
              IF( ABS(AM0*BM1-AM1*BM0).LE.EPS*AM0*BM1 ) GOTO 2400
              AM2 = AM1
              AM1 = AM0
              BM2 = BM1
              BM1 = BM0
              GOTO 2200
 2400     CONTINUE
          EN(N) = AM0/BM0
 3000 CONTINUE
C
C    Rescale computed continous fractions. The loop below
C    *can't* use Y**N, or overflow may occur for valid
C    inputs.
C
      YPOWN = Y*EXPY
      DO 3500 N=1,NMIN
          YPOWN = YPOWN*Y/N
 3500 CONTINUE
      DO 4000 N=NMIN,NMAX
          EN(N) = YPOWN/EN(N)
          YPOWN = YPOWN*Y/(N+1)
 4000 CONTINUE
      END
C
      SUBROUTINE PS3GMI(KMIN,KMAX,NX,XVAL,Y,FUN,LDF)
C
C   Compute batch of values for the generalized incomplete gamma
C   function of positive integer order:
C
C      Gamma[k,x,x+y] = Integrate[ t^(k-1) Exp[-t], {t,x,x+y} ]
C
C   Result is expected to be good to (at least) machine accuracy
C   sans one digit.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      KMIN   - Minimum order of Gamma required (at least 1)
C      KMAX   - Maximum order of Gamma required
C      NX     - Number of X argument values
C      XVAL   - X argument values (each X >= 0)
C      Y      - Y argument values (Y > 0)
C      FUN    - Two-dimensional array for function values
C               First index: order of Gamma function, KMIN to KMAX
C               Second index: values of argument X, at least NX
C      LDF    - Leading dimension of the FUN array.
C      
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C      MAXO DOUBLE PRECISION words (ca. 1.6K)
C
C   Module logic:
C
C      Generalized incomplete gamma function is computed according
C      to the (finite) series expansion:
C
C      gamma[k,x,x+y] = Exp[-x] (k-1)! Sum[ x^i/i! Ex[k-1-i,y], {i,0,k-1} ],
C
C      where Ex[n,y] factors are provides by PS3GME above. The
C      expression is based on Abramowitz & Stegun (although not
C      directly taken from there).
C
C   Possible optimizations:
C
C   Bugs:
C
C      Ordex of the gamma functions computed is (arbitrarily) fixed
C      at MAXO (=100).
C
C      Computation can overflow or underflow for some valid inputs
C      and represntable outputs.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      PARAMETER (MAXO=200)
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
C
      DIMENSION XVAL(NX), FUN(KMIN:LDF+KMIN-1,NX)
      DIMENSION EX(0:MAXO-1)
C
C    Do a little error-checking
C
      IF( KMIN.LT.1 .OR. KMIN.GT.KMAX ) THEN
          WRITE(NB6,11000) KMIN
          STOP 'PS3GMI'
      ENDIF
      IF( KMAX.GT.MAXO ) THEN
          WRITE(NB6,11010) KMAX, MAXO
          STOP 'PS3GMI'
      ENDIF
      IF( NX.LE.0 ) THEN
          WRITE(NB6,11020) NX
          STOP 'PS3GMI'
      ENDIF
      IF( LDF.LT.KMAX-KMIN+1 ) THEN
          WRITE(NB6,11030) LDF, KMAX-KMIN+1
          STOP 'PS3GMI'
      ENDIF
      IF( Y.LT.ZERO ) THEN
          WRITE(NB6,11040) Y
          STOP 'PS3GMI'
      ENDIF
C
      CALL PS3GME(KMAX-1,Y,EX)
      DO 3000 IX=1,NX
          DO 1000 K=KMIN,KMAX
              FUN(K,IX) = ZERO
 1000     CONTINUE
          X    = XVAL(IX)
          IF( X.LT.ZERO ) THEN
              WRITE(NB6,11050) IX, X
              STOP 'PS3GMI'
          ENDIF
          SCLI = EXP(-X)
          TERM = ONE
C
          DO 2000 I=1,KMIN-1
              DO 1500 K=KMIN,KMAX
                  FUN(K,IX) = FUN(K,IX) + EX(K-I)*TERM
 1500         CONTINUE
              SCLI = SCLI*I
              TERM = TERM*X/I
 2000     CONTINUE
C
          DO 2600 I=KMIN,KMAX
              DO 2100 K=I,KMAX
                  FUN(K,IX) = FUN(K,IX) + EX(K-I)*TERM
 2100         CONTINUE
              FUN(I,IX) = FUN(I,IX) * SCLI
              SCLI = SCLI*I
              TERM = TERM*X/I
 2600     CONTINUE
C
 3000 CONTINUE
C
      RETURN
11000 FORMAT(' KMIN VALUE (',I5,') IS EITHER NON-POSITIVE OR GREATER ',
     .       'THAN KMAX IN PS3GMI.')
11010 FORMAT(' KMAX VALUE (',I5,') EXCEEDS COMPILE-TIME LIMIT OF ',
     .       I5,' IN PS3GMI.')
11020 FORMAT(' NUMBER OF X ARGUMENTS (',I5,
     .       ') IS NON-POSITIVE IN PS3GMI.')
11030 FORMAT(' LEADING DIMENSION OF THE OUTPUT ARRAY (',I5,
     .       ') IS TOO SMALL IN PS3GMI. AT LEAST ',I5,' IS NECESSARY.')
11040 FORMAT(' EXTENT ARGUMENT (',G20.15,') IS NEGATIVE IN PS3GMI.')
11050 FORMAT(' ',I5,'-TH BASE ARGUMENT (',G20.15,
     .       ') IS NEGATIVE IN PS3GMI.')
      END
C
      SUBROUTINE PS3SHB(MAXS,MAXM,NR,XS,CS,LDC,BPOLY,LDB1,LDB2)
C
C   Compute values of the auxiliary polynomial in Sharma's expansion.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      MAXS   - Maximum order of the polynomial.
C      MAXM   - Maximum magnetic quantum number.
C      NR     - Number of radial points to evaluate polynomial at.
C      XS     - (Array of) characteristic variable value ((r/a)^2) 
C               at each radial point, NR values.
C      CS     - Polynomial coefficients, three-index array:
C               1: (0:MAXS-IS) Summation variable NU
C               2: (0:MAXS) Index IS
C               3: (0:MAXM) Magnetic quantum number
C      LDC    - Two leading dimensions of the CS matrix
C      BPOLY  - Output polynomials, three-index array
C               1: (0:MAXS) Index IS
C               2: (0:MAXM) Magnetic quantum number
C               3: (1:NR) Radial point
C      LDB1   - First leading dimension of BPOLY
C      LDB2   - Second leading dimension of BPOLY
C
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C   Module logic:
C
C      Polynomials are computed using the Gorner expression.
C      Loops are arranged so that the number of read cache 
C      misses will be minimized provided that a single set
C      of coefficients will fit in cache.
C
C   Possible optimizations:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0)
      DIMENSION XS(NR), CS(0:LDC-1,0:LDC-1,0:MAXM)
      DIMENSION BPOLY(0:LDB1-1,0:LDB2-1,NR)
C
      DO 5000 M=0,MAXM
          DO 4000 IS=0,MAXS
              DO 3000 IR=1,NR
                  X   = XS(IR)
                  VAL = ZERO
                  DO 2000 NU=MAXS-IS,0,-1
                      VAL = VAL*X+CS(NU,IS,M)
 2000             CONTINUE
                  BPOLY(IS,M,IR) = VAL
 3000         CONTINUE
 4000     CONTINUE
 5000 CONTINUE
      END
C
      SUBROUTINE PS3SHP(N,L,ZETA,LS,NR,RS,FUN,LDF)
C
C   Compute values of the Sharma's expansion of a given order
C   at a set of points. Special case of A.EQ.ZERO
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N,L    - Major and orbital quantum number of the orbital
C               being expanded. L should be smaller than N.
C      ZETA   - Orbital exponent.
C      LS     - Order of the expansion term.
C      NR     - Number of radial points where function should be
C               computed
C      RS     - Array of radial arguments. Arguments should be
C               in an increasing order. (More strictly speaking,
C               it is sufficient to have all R<A preceed all R>=A)
C      FUN    - Two-dimensional array for computed expansion
C               function values.
C               First dimension:
C                  1 = value of the expansion function
C                  2 = bound for the error in the expansion function
C                      expansion is exact in this particular case,
C                      so the error estimation is zero.
C               Second dimension: absolute value of the magnetic 
C                  quantum number M, in the order 0 to L (i.e., at 
C                  least L + 1 entries). 
C               Third dimension: radial argument values (i.e.,
C                  at least NR entries)
C      LDF    - Second leading dimension of the FUN matrix.
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
C      Expansion with respect to the natural origin is equal to the function
C      itself if LS is equal to L, and is zero otherwise.
C
C   Possible optimizations:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0)
C
      DIMENSION RS(NR), FUN(2,0:LDF-1,NR)
C
      DO 200 I=1,NR
          VAL = ZERO
          IF( LS.EQ.L ) VAL = (RS(I)**N)*EXP(-ZETA*RS(I))
          DO 190 M=0,L
              FUN(1,M,I) = VAL
              FUN(2,M,I) = ZERO
  190     CONTINUE
  200 CONTINUE
      RETURN
      END
C
      SUBROUTINE PS3SHQ(N,L,ZETA,LS,A,NR,RS,GAMMAS,LDG,KHAVE,FUN,LDF)
C
C   Compute values of the Sharma's expansion of a given order
C   at a set of points. This function is not intended to be directly
C   callable, call PS3SHR instead.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N,L    - Major and orbital quantum number of the orbital
C               being expanded. L should be smaller than N.
C      ZETA   - Orbital exponent.
C      LS     - Order of the expansion term
C      A      - Separation of the centers. A is expected to be non-zero.
C      NR     - Number of radial points where function should be
C               computed
C      RS     - Array of radial arguments. Arguments should be
C               in an increasing order. (More strictly speaking,
C               it is sufficient to have all R<A preceed all R>=A)
C      GAMMAS - Input/Output array of the incomplete gamma functions.
C               First index: order of the gamma function, one-based.
C               Second index: argument value, one-to one correspondence
C                   with arguments.
C      LDG    - Leading dimension of the GAMMAS array.
C      KHAVE  - Input/Output: maximum order of the gamma function in
C                   the GAMMAS, should not exceed LDG.
C      FUN    - Two-dimensional array for computed expansion
C               function values.
C               First dimension:
C                  1 = value of the expansion function
C                  2 = bound for the error in the expansion function
C               Second dimension: absolute value of the magnetic 
C                  quantum number M, in the order 0 to L (i.e., at 
C                  least L + 1 entries). 
C               Third dimension: radial argument values (i.e.,
C                  at least NR entries)
C      LDF    - Second leading dimension of the FUN matrix.
C      
C   Accessed common blocks:
C
C      PS3SHD,
C      PS3SHE - Values of the expansion coefficients, defined in
C               psov3d.f's PS3SHS.
C      PSPRT  - Printing unit.
C
C   Local storage:
C
C      MAXNR*(1+(MAXL+1)*(MAXL+1+MAXLS)) DOUBLE PRECISION storage
C      units. This is some 6800 DP words.
C
C   Module logic:
C
C      The code is an (almost) straightforward implementation of the eq.
C      (21b) of Sharma's paper.
C
C      EPS should be set to reflect accuracy of the machine arithmetics.
C      1D-16 should be Ok for IEEE doubles.
C
C   Possible optimizations:
C
C      Fine-tuning KINC1 and KINCX might improve performance slightly.
C
C   Bugs:
C
C      L is limited to MAXL (=3). LS is limited to MAXLS (=15).
C      These limits are not entirely arbitrary, since the analytical
C      expressions for expansion coefficients are numerically unstable
C      and have to be precalculated at higher precision than is available
C      using DOUBLE PRECISION on most architectures. Limit of MAXLS
C      should not present any real problem, since the series is too
C      numerically unstable beyond that point.
C
C      Number of radial points is limited to MAXNR (=70). This limit
C      is arbitrary.
C
C      Major quantum number is limited to MAXN (=6). This limit is
C      arbitrary as well.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXN=6)
      PARAMETER (MAXL=3)
      PARAMETER (MAXLS=40)
      PARAMETER (MAXNR=80)
      PARAMETER (LCOEFF=274700)
      PARAMETER (KINC1=12)
      PARAMETER (KINCX=6)
      PARAMETER (EPS=1.D-16)
      PARAMETER (ZERO=0.D0)
      PARAMETER (TWO=2.D0)
      EXTERNAL PS3SHS
C
      COMMON
     ./PS3SHD/ COEFF(LCOEFF)
     ./PS3SHE/ IBASE(0:MAXLS,0:MAXL)
     ./PSPRT / NB6
      SAVE /PS3SHD/, /PS3SHE/, /PSPRT /
C
      DIMENSION COEFX(((MAXL+MAXLS+1)**2)*(MAXL+1))
      DIMENSION RS(NR), FUN(2,0:LDF-1,NR), GAMMAS(LDG,NR)
C
      DIMENSION BASR(MAXNR)
      DIMENSION BPOLY(0:MAXL+MAXLS,0:MAXL,MAXNR)
C    (A little of?) range checking
      IF( L.LT.0 .OR. L.GT.MAXL ) THEN
          WRITE(NB6,11000) L
          STOP 'PS3SHQ'
      ENDIF
      IF( LS.LT.0 .OR. LS.GT.MAXLS ) THEN
          WRITE(NB6,11010) LS
          STOP 'PS3SHQ'
      ENDIF
      IF( NR.LE.0 .OR. NR.GT.MAXNR ) THEN
          WRITE(NB6,11020) NR
          STOP 'PS3SHQ'
      ENDIF
      IF( LDF.LT.L+1 ) THEN
          WRITE(NB6,11030) LDF, L+1
          STOP 'PS3SHQ'
      ENDIF
      IF( IBASE(0,0).EQ.0 ) THEN
          WRITE(NB6,11040)
          STOP 'PS3SHQ'
      ENDIF
      IF( N.LE.0 .OR. N.GT.MAXN ) THEN
          WRITE(NB6,11050) N
          STOP 'PS3SHQ'
      ENDIF
      IF( L.GE.N ) THEN
          WRITE(NB6,11060) L, N
          STOP 'PS3SHQ'
      ENDIF
      IF( A.EQ.ZERO ) THEN
          WRITE(NB6,11062)
          STOP 'PS3SHQ'
      ENDIF
C
C    Evaluate gamma functions we'll need. Ones with
C    R < A have to be computed using order batching,
C    while R >= A allows radial batching as well.
      ABSA = ABS(A)
      MAXS = L+LS
      KMIN = N-L+1
      KMAX = KMIN+2*MAXS
      IF( KMAX.GT.KHAVE ) THEN
C        We had not encountered gammas that big before, so we'll
C        have to do them the hard way.
          KMIN = MAX(KHAVE+1,KMIN)
          IF( KMAX.GT.LDG ) THEN
              WRITE(NB6,11067) KMAX, LDG
              STOP 'PS3SHQ'
          ENDIF
          IF( KHAVE.EQ.0 ) THEN
             KHAVE = MIN(LDG,KMAX+KINC1)
          ELSE
             KHAVE = MIN(LDG,KMAX+KINCX)
          ENDIF
C
          DO 1000 I=1,NR
              R = RS(I)
              IF( R.LE.ZERO ) THEN
                  WRITE(NB6,11070) I, R
                  STOP 'PS3SHQ'
              ENDIF
              IF( R.LT.ABSA ) THEN
                  Y       = TWO*R*ZETA
                  BASR(I) = (ABSA-R)*ZETA
                  CALL PS3GMI(KMIN,KHAVE,1,BASR(I),Y,GAMMAS(KMIN,I),LDG)
              ELSE
                  GOTO 1100
              ENDIF
 1000     CONTINUE
          GOTO 1300
C        Rest of radial points can be handled as a single batch,
C        but we'd better make sure the calling guy is not cheating
C        first...
 1100     CONTINUE
          N1RBIG = I
          DO 1200 I=N1RBIG,NR
              R = RS(I)
              IF( R.LT.ABSA ) THEN
                  WRITE(NB6,11080) I, R
                  STOP 'PS3SHQ'
              ENDIF
              BASR(I) = (R - ABSA)*ZETA
 1200     CONTINUE
          CALL PS3GMI(KMIN,KHAVE,NR-N1RBIG+1,BASR(N1RBIG),TWO*ABSA*ZETA,
     .                GAMMAS(KMIN,N1RBIG),LDG)
 1300     CONTINUE
      ENDIF
C    Evaluate b_{\nu} polynomials for the entire batch of
C    integrals.
      DO 2000 I=1,NR
          BASR(I) = (RS(I)/ABSA)**2
 2000 CONTINUE
      IF( IBASE(LS,L).NE.0 ) THEN
C        Are coefficients already here? If not, compute and store'em
          IF( IBASE(LS,L).LT.0 ) THEN
              IBASE(LS,L) = ABS(IBASE(LS,L))
              CALL PS3SBN(LS,L,COEFF(IBASE(LS,L)),LS+L,LS+L)
          ENDIF
C        Coefficients already here, take an easy way out
          CALL PS3SHB(MAXS,L,NR,BASR,COEFF(IBASE(LS,L)),LS+L+1,
     .                BPOLY,MAXL+MAXLS+1,MAXL+1)
      ELSE
C        Compute coefficients in the local buffer
          CALL PS3SBN(LS,L,COEFX,LS+L,LS+L)
          CALL PS3SHB(MAXS,L,NR,BASR,COEFX,LS+L+1,
     .                BPOLY,MAXL+MAXLS+1,MAXL+1)
      ENDIF
C    Scale b_{\nu} values with the part of the prefactor independent
C    of R
      SCALEX = ((ZETA*ABSA)**(L-1))/(ZETA**N)
      IF( A.LT.ZERO .AND. MOD(MAXS,2).EQ.1 ) SCALEX = -SCALEX
      DO 3000 IS=0,MAXS
          SCALE = SCALEX/(ZETA*ABSA)**(2*IS)
          DO 2900 IR=1,NR
              DO 2800 M=0,L
                  BPOLY(IS,M,IR) = SCALE*BPOLY(IS,M,IR)
 2800         CONTINUE
 2900     CONTINUE
 3000 CONTINUE
C    Compute values of the expansion function and error
C    estimation.
      DO 6000 IR=1,NR
          R = RS(IR)
          SCALE = (ABSA/R)**LS
          DO 5000 M=0,L
              VAL  = ZERO
              AVAL = ZERO
              DO 4000 IS=0,MAXS
                  TERM = BPOLY(IS,M,IR)*GAMMAS(N-L+1+2*IS,IR)
                  VAL  = VAL + TERM
                  AVAL = AVAL + ABS(TERM)
 4000         CONTINUE
              FUN(1,M,IR) = SCALE*VAL
              FUN(2,M,IR) = EPS*SCALE*AVAL
 5000     CONTINUE
 6000 CONTINUE
      RETURN
11000 FORMAT(' L = ',I5,' IS NOT SUPPORTED IN PS3SHQ.')
11010 FORMAT(' LS = ',I5,' IS NOT SUPPORTED IN PS3SHQ.')
11020 FORMAT(' NR = ',I5,' IS INVALID IN PS3SHQ.')
11030 FORMAT(' LDF = ',I5,' IS INSUFFICIENT IN PS3SHQ.',
     .      /' AT LEAST ',I5,' IS NECESSARY.' )
11040 FORMAT(' BLOCK DATA PS3SHS WAS OMITTED AT LINK TIME. ',
     .       'PS3SHQ CANNOT PROCEED.')
11050 FORMAT(' MAJOR QUANTUM NUMBER N (',I5,') IS INVALID IN PS3SHQ.')
11060 FORMAT(' ORBITAL QUANTUM NUMBER L (',I5,
     .       ') IS INVALID FOR N = ',I5)
11062 FORMAT(' CENTER SEPARATION A IS ZERO IN PS3SHQ.')
11065 FORMAT(' ZERO INTERNUCLEAR DISTANCE FOR PS3SHQ EXPANSION.')
11067 FORMAT(' MAXIMUM OREDER OF GAMMA FUNCTION IS ',I5,
     .       ' IN PS3SHQ, WHILE THERE IS ONLY ',I5,' SLOTS.')
11070 FORMAT(1X,I5,'-TH RADIAL ARGUMENT IS NON-POSITIVE (',G25.15,
     .       ') IN PS3SHQ.')
11080 FORMAT(1X,I5,'-TH RADIAL ARGUMENT (',G25.15,
     .       ') IS SMALLER THAN THE INTERCENTER DISTANCE IN PS3SHQ.')
      END
C
      SUBROUTINE PS3SHR(N,L,ZETA,LS,A,NR,RS,FUN,LDF,ISAVE,LISAVE,
     .                  DSAVE,LDSAVE)
C
C   Compute values of the Sharma's expansion of a given order
C   at a set of points. This is a high-level driver for PS3SHQ.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N,L    - Major and orbital quantum number of the orbital
C               being expanded. L should be smaller than N.
C      ZETA   - Orbital exponent.
C      LS     - Order of the expansion term
C      A      - Separation of the centers
C      NR     - Number of radial points where function should be
C               computed
C      RS     - Array of radial arguments. Arguments should be
C               in an increasing order. (More strictly speaking,
C               it is sufficient to have all R<A preceed all R>=A)
C      FUN    - Two-dimensional array for computed expansion
C               function values.
C               First dimension: absolute value of the magnetic 
C                  quantum number M, in the order 0 to L (i.e., at 
C                  least L + 1 entries). 
C               Second dimension: radial argument values (i.e.,
C                  at least NR entries)
C      LDF    - Leading dimension of the FUN matrix.
C      ISAVE  - Integer persistent store, should be preserved
C               between calls to PS3SHR with the same ZETA's, A's
C               and radial coordinates. First element of the ISAVE
C               is expected to be set to -1 on the first call of the
C               series.
C      LISAVE - Length of the ISAVE array.
C      DSAVE  - Double precision persistent store.
C      LDSAVE - Length of the DSAVE array.
C      
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C      2*(MAXL+1)*MAXNR DOUBLE PRECISION storage units. This is
C      some 700 DP words.
C
C   Module logic:
C
C      This is a caching driver for the true Sharma's expansion
C      evaluator, PS3SHQ. Our job is:
C
C      1. Check cache for the presence of expansion function
C         of a given order evaluated previously, and return that
C         if we can.
C      2. If we can't, we should arrange for a call to PS3SHQ
C         reusing previously evaluated incomplete gamma functions
C         (which are over 70% of the effort in the uncached version).
C      3. Scan PS3SHQ's values for the loss of precision and discard
C         bad entries.
C      4. Store computed values in the cache. Cache is organized as
C         a FIFO. Current length (3) is sufficient to handle integrals
C         with dipole on operator center.
C
C      Internally, expansion function values are evaluated with an
C      error bar attached. However, after eliminating values with
C      too low relative precision, error bars are discarded and only
C      function values are available to the caller.
C
C      We also have quite a bit of persistent storage in ISAVE and DSAVE,
C      which are organized as follows:
C
C       Index  Object            Comments
C
C       IBSIGN ISIGN             Signature word. Any value but -1 or ISG is an error.
C       IBHAVE KHAVE             Highest order of gamma function present in DGAMMA.
C       IBICH  ICH               Next cache slot to be used on store, range is
C                                    0 to LCACHE - 1.
C       IBLSMX LSMX              Largest LS seen in this batch.
C       IBLSCH LSCH(LCACHE)      LS values for the cached expansion function sets.
C       IBLOW  IGOODL(0:MAXLS)   Index of the smallest argument value which is still
C                                    sufficiently numerically stable at given LS.
C                                    Zero if unknown.
C       IBHIGH IGOODH(0:MAXLS)   Index of the largest numerically stable argument value.
C
C       IBGAMM DGAMMA(LDG,MAXNR) Values of the auxiliary incomplete gamma function, i.e.
C                                Gamma[n,zeta*Abs[r-a],zeta*Abs[r+a]]
C                                    First index: order of the gamma function
C                                    Second index: argument value.
C       IBCH   CACHE(0:MAXL,MAXNR,LCACHE)
C                                Cached values of the expansion function.
C                                    First index: M parameter of the cached values.
C                                    Second index: Argument value.
C                                    Third index: Cache line.
C
C   Possible optimizations:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXN=6)
      PARAMETER (MAXL=3)
      PARAMETER (MAXLS=40)
      PARAMETER (MAXNR=80)
      PARAMETER (LDG=MAXN+2*(MAXL+MAXLS))
      PARAMETER (PLOSS =1.D-2)
      PARAMETER (ZERO=0.D0)
C    Number of cached function values sets.
      PARAMETER (LCACHE=3)
C    Our cache layout
      PARAMETER (IBSIGN=1)
      PARAMETER (IBHAVE=2)
      PARAMETER (IBICH =3)
      PARAMETER (IBLSMX=4)
      PARAMETER (IBLSCH=5)
      PARAMETER (IBLOW =IBLSCH+LCACHE)
      PARAMETER (IBHIGH=IBLOW+(MAXLS+1))
      PARAMETER (IICNT =IBHIGH+(MAXLS+1)-1)
C
      PARAMETER (IBGAMM=1)
      PARAMETER (IBCH  =IBGAMM+LDG*MAXNR)
      PARAMETER (IDCNT =IBCH+LCACHE*(MAXL+1)*MAXNR)
C
      EXTERNAL PS3SHS
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
      DIMENSION ISAVE(LISAVE), DSAVE(LDSAVE)
      DIMENSION RS(NR), FUN(0:LDF-1,NR)
C
      DIMENSION TFUN(2,0:MAXL,MAXNR)
C
      DATA ISG/'68'/
C    Check and initialize cache area
      IF( LISAVE.LT.IICNT .OR. LDSAVE.LT.IDCNT ) THEN
          WRITE(NB6,10010) IICNT, IDCNT, LISAVE, LDSAVE
          STOP 'PS3SHR'
      ENDIF
      IF( ISAVE(IBSIGN).NE.ISG ) THEN
C        Initialize cache area
          IF( ISAVE(IBSIGN).NE.-1 ) THEN
              WRITE(NB6,10020) ISAVE(IBSIGN)
              STOP 'PS3SHR'
          ENDIF
          ISAVE(IBHAVE) =  0
          ISAVE(IBLSMX) = -1
C        The upper limit is not in error, there are MAXLS+1 possible 
C        expansion orders.
          DO 90 I=0,MAXLS
              ISAVE(IBLOW +I) = 0
              ISAVE(IBHIGH+I) = 0
   90     CONTINUE
          ISAVE(IBICH ) =  0
          DO 100 I=0,LCACHE-1
              ISAVE(IBLSCH+I) = -1
  100     CONTINUE
          ISAVE(IBSIGN) = ISG
      ENDIF
C    (A little of?) Range checking
      IF( L.LT.0 .OR. L.GT.MAXL ) THEN
          WRITE(NB6,11000) L
          STOP 'PS3SHR'
      ENDIF
      IF( NR.LE.0 .OR. NR.GT.MAXNR ) THEN
          WRITE(NB6,11020) NR
          STOP 'PS3SHR'
      ENDIF
      IF( LDF.LT.L+1 ) THEN
          WRITE(NB6,11030) LDF, L+1
          STOP 'PS3SHR'
      ENDIF
      IF( N.LE.0 .OR. N.GT.MAXN ) THEN
          WRITE(NB6,11050) N
          STOP 'PS3SHR'
      ENDIF
C    Lookup our little cache first
      DO 1000 I=0,LCACHE-1
          IF( ISAVE(IBLSCH+I).EQ.LS ) THEN
C        Hurra, we'd got it in the cache!
              DO 900 IR=1,NR
                  DO 800 M=0,L
C                    DSAVE is adressed as a linearized three-dimensional array:
C                    first dimension is 0:MAXL, indexed by M
C                    second dimension in 1:MAXNR, indexed by IR
C                    third dimension is 0:LCACHE-1, indexed by I
                      FUN(M,IR) = DSAVE(IBCH + M + (IR-1)*(MAXL+1) + 
     .                                         I*(MAXL+1)*MAXNR)
  800             CONTINUE
  900         CONTINUE
              RETURN
          ENDIF
 1000 CONTINUE
C    Not in the cache yet, lets find expected stability region first.
C    We know that the stability region can only shrink with increasing
C    expansion order.
      IGOODL = 1
      IGOODH = NR
      DO 2000 I=MIN(LS,ISAVE(IBLSMX)),0,-1
          IF( ISAVE(IBLOW+I).NE.0 ) THEN
              IGOODL = ISAVE(IBLOW +I)
              IGOODH = ISAVE(IBHIGH+I)
              GOTO 2001
          ENDIF
 2000 CONTINUE
 2001 CONTINUE
C    Zero out omitted parts
      DO 2500 I=1,IGOODL-1
          DO 2490 M=0,L
              TFUN(1,M,I) = ZERO
 2490     CONTINUE
 2500 CONTINUE
      DO 2600 I=IGOODH+1,NR
          DO 2590 M=0,L
              TFUN(1,M,I) = ZERO
 2590     CONTINUE
 2600 CONTINUE
C    Perform the evalustion.
      IF( IGOODH.GE.IGOODL ) THEN
          IF( A.EQ.ZERO ) THEN
              CALL PS3SHP(N,L,ZETA,LS,NR,RS(IGOODL),
     .                    TFUN(1,0,IGOODL),MAXL+1)
          ELSE IF( MOD(INTCTL,100).EQ.1 ) THEN
              CALL PS3SHQ(N,L,ZETA,LS,A,IGOODH-IGOODL+1,RS(IGOODL),
     .            DSAVE(IBGAMM+(IGOODL-1)*LDG),LDG,ISAVE(IBHAVE),
     .            TFUN(1,0,IGOODL),MAXL+1)
          ELSE IF( MOD(INTCTL,100).EQ.2 ) THEN
              CALL PS3BA (N,L,ZETA,LS,A,IGOODH-IGOODL+1,RS(IGOODL),
     .            DSAVE(IBGAMM+(IGOODL-1)*LDG),LDG,ISAVE(IBHAVE),
     .            TFUN(1,0,IGOODL),MAXL+1)
          ELSE
              WRITE(NB6,11055) INTCTL
              STOP 'PS3SHR'
          ENDIF
      ENDIF
C    Screen computed values for instabilities
      NGOODL = IGOODH + 1
      NGOODH = IGOODL - 1
      DO 3000 I=IGOODL,IGOODH
          DO 2900 M=0,L
              IF( ABS(PLOSS*TFUN(1,M,I)).LT.ABS(TFUN(2,M,I)) ) THEN
                  TFUN(1,M,I) = ZERO
              ELSE
                  NGOODL = MIN(I,NGOODL)
                  NGOODH = MAX(I,NGOODH)
              ENDIF
 2900     CONTINUE
 3000 CONTINUE
C
      IF( NGOODL.GT.NGOODH ) THEN
*???? This is not a good idea to shut TLOSS message off, but *I* already know it
*???? is here, so what the deal?
          WRITE(NB6,11060) N,L,ZETA,LS,A
          NGOODL = IGOODL
          NGOODH = NGOODH
      ENDIF
C    Record new range for the future use
      ISAVE(IBLSMX   ) = MAX(LS,ISAVE(IBLSMX))
      ISAVE(IBLOW +LS) = NGOODL
      ISAVE(IBHIGH+LS) = NGOODH
C    Save function value into a cache slot and user buffer
      I = ISAVE(IBICH)
      DO 4900 IR=1,NR
          DO 4800 M=0,L
              FUN(M,IR) = TFUN(1,M,IR)
              DSAVE(IBCH + M + (IR-1)*(MAXL+1) + I*(MAXL+1)*MAXNR) 
     .                  = TFUN(1,M,IR)
 4800     CONTINUE
 4900 CONTINUE
C    Save cache id and advance free slot number
      ISAVE(IBLSCH+I) = LS
      I = I + 1
      IF( I.GE.LCACHE ) I = 0
      ISAVE(IBICH) = I
C    Done!
      RETURN
10010 FORMAT(' AT LEAST ',I5,' INTEGER AND ',I5,' DOUBLE PRECISION ',
     .       'WORDS ARE NEEDED IN PS3SHR. ONLY ',I5,' AND ',I5,
     .       ' ARE AVAILABLE.')
10020 FORMAT(' CALLER OF PS3SHR IS CLUELESS: SIGNATURE: ',I9)
11000 FORMAT(' L = ',I5,' IS NOT SUPPORTED IN PS3SHR.')
11020 FORMAT(' NR = ',I5,' IS INVALID IN PS3SHR.')
11030 FORMAT(' LDF = ',I5,' IS INSUFFICIENT IN PS3SHR.',
     .      /' AT LEAST ',I5,' IS NECESSARY.' )
11050 FORMAT(' MAJOR QUANTUM NUMBER N (',I5,') IS INVALID IN PS3SHR.')
11055 FORMAT(' INTCTL VALUE (',I5,') IS INVALID IN PS3SHR.')
11060 FORMAT(' TOTAL LOSS OF PRECISION IN PS3SHR.'/' N = ',I3,' L = ',
     .       I3,' ZETA = ',G18.12,' LS = ',I3,' A = ',G18.12)
      END
C
      SUBROUTINE PS3SCR(LS1,LS2,F,N1,L1,ALP1,A1,N2,L2,ALP2,A2,AINT,LDA,
     .                  ISAVE,LISAVE,DSAVE,LDSAVE)
C
C   Compute block of integrals over product of two Sharma's expansions
C   and arbitrary radial function.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      LS1    - Order of the first expansion term.
C      LS2    - Order of the second expansion term.
C      F      - External subroutine used to compute radial function.
C               F is expected to take a double precision argument (r)
C               and should return double precision value of the
C               radial function at the given point.
C      N1     - Major quantum number of the first expanded orbital.
C      L1     - Orbital quantum number of the first expanded orbital.
C      ALP1   - Orbital exponent of the first expanded orbital.
C      A1     - Displacement of the expansion center from the first
C               expanded orbital.
C      N2,L2,ALP2,A2
C             - Same as N1,L1,ALP1,A1 for the second expanded orbital.
C      AINT   - (Output) place for the computed integrals.
C               First dimension - magnetic quantum number of the first
C                   expanded orbital, 0 to L1.
C               Second dimension - magnetic quantum number of the second
C                   expanded orbital, 0 to L2.
C               Integrals are symmetric with respect to the signs of the
C               magnetic quantum numbers.
C      LDA    - Leading dimension of the AINT matrix.
C      ISAVE  - Temporary array preserved between calls to PS3SCR with
C               the same F,N*,L*,ALP* and A* parameters. First element of 
C               the ISAVE array should contain -1 if it is the very
C               first call in the series.
C      LISAVE - Length of the ISAVE array.
C      DSAVE  - Temporary array preserved between calls to PS3SCR, on the
C               same conditions as the ISAVE array.
C      LDSAVE - Length of the LDSAVE array.
C      
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C      ca. 2*MAXPTS*(MAXL+2) DOUBLE PRECISION units, i.e. some 600 DP
C      words for MAXL=3 and MAXPTS=60
C
C   Module logic:
C
C      The numerical integration used here is straightforward. The semi-
C      infinite region is divided into up to seven integration domains
C      (given below assuming A1<=A2. Second case is left as an excersize
C      to an interested reader ;-). R1L is defined as a distance at which
C      exponent ALP1 falls to BRD1. R1H is a distance at which EXP(-ALP1*R)
C      falls to BRD2. Obviously, BRD1 should be greater than BRD2.
C
C      a. Near-zero partition, 0 to A1 - R1H. IOZERO-points Gauss' formula
C         is used here. If R1H is >= A1, this partition is dropped.
C      b. Low arguments border partion, A1 - R1H to A1 - R1L. IOBORD-points
C         Gauss' formula is used here. If R1L is >= A1, this partition is
C         dropped.
C      c. First center partition, A1 - R1L to A1 + R1L, IOCENT-points Gauss'
C         integration.
C      d. Inter-center partition, A1 + R1L to A2 - R2L, IOSPAC-points Gauss'
C         integration.
C      e. Second-center partition, A2 - R2L to A2 + R2L, IOCENT-points Gauss'
C         integration. If partitions c and e overlap, single partition from
C         A1 - R1L to A2 + R2L is used instead of c, d and e for 2*IOCENT + 
C         IOSPAC-points Gauss' integration.
C      f. High-argument border partition, A2 + R2L to A2 + R2H. IOBORD-points
C         Gauss' integration.
C      g. Exponential tail partition, A2 + R2H to infinity. IOTAIL-points
C         Laguerre's integration with the exponent (ALP1+ALP2) is used here.
C
C      Sharma's expansion integrated here is numerically unstable for both
C      small and large arguments, especially at high expansion orders. 
C      Therefore, function values which are accurate to less than PLOSS
C      are discarded before integration is performed. For all interesting
C      integrands, this should not deteriorate integral precision significantly
C      since subintegrals in these regions are very small.
C
C   Accuracy:
C
C      Computed integrals are accurate to approximately 10^-5 *absolute*
C      once integrals are normalized. Accuracy can be cranked up to
C      approximately 10^-7 to 10^-8, but after that point numerical
C      instabilities in the Sharma's expansion at higher orders become
C      too pronounced.
C
C   Possible optimizations:
C
C      Guessing which function values are going to loose precision and
C      just omitting 'em should increase performance for higher expansion
C      terms.
C
C   Bugs:
C
C      Convergence of the integrals should be monitored, rather than
C      relied upon.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=3)
C
      PARAMETER (ZERO=0.D0)
C
      PARAMETER (BRD1  =3.D-1)
      PARAMETER (BRD2  =3.D-3)
      PARAMETER (IOZERO= 6)
      PARAMETER (IOBORD=12)
      PARAMETER (IOCENT=12)
      PARAMETER (IOSPAC=12)
      PARAMETER (IOTAIL=12)
C
      PARAMETER (MAXBLK=6)
      PARAMETER (MAXPTS=IOZERO+2*IOBORD+2*IOCENT+IOSPAC+IOTAIL)
C    Save area requirements for the auxiliary routines.
      PARAMETER (IISHC =90)
      PARAMETER (IDSHC =8350)
C    Layour of the save area.
      PARAMETER (ISIGN =1)
      PARAMETER (IBIO  =2)
      PARAMETER (IISHR1=3)
      PARAMETER (IISHR2=IISHR1+IISHC)
      PARAMETER (IICNT =IISHR2+IISHC-1)
C
      PARAMETER (IBARGS=1)
      PARAMETER (IBWGTS=IBARGS+MAXPTS)
      PARAMETER (IDSHR1=IBWGTS+MAXPTS)
      PARAMETER (IDSHR2=IDSHR1+IDSHC)
      PARAMETER (IDCNT =IDSHR2+IDSHC-1)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      EXTERNAL F
C
      DIMENSION AINT(0:LDA-1,0:L2)
      DIMENSION ISAVE(LISAVE), DSAVE(LDSAVE)
C
      DIMENSION VAL1(0:MAXL,MAXPTS), VAL2(0:MAXL,MAXPTS)
      DIMENSION PART(0:MAXBLK), IGO(MAXBLK)
      DATA ISG/'PS'/
C    Do we have enough space in the restart arrays?
      IF( LISAVE.LT.IICNT .OR. LDSAVE.LT.IDCNT ) THEN
          WRITE(NB6,10010) IICNT, IDCNT, LISAVE, LDSAVE
          STOP 'PS3SCR'
      ENDIF
C    Do we need to do initialization the hard way?
      IF( ISAVE(ISIGN).NE.ISG ) THEN
C        Were we called by someone with a clue?
          IF( ISAVE(ISIGN).NE.-1 ) THEN
              WRITE(NB6,10020) ISAVE(ISIGN)
              STOP 'PS3SCR'
          ENDIF
C        Make sure PS3SHR won't complain that we don't have a clue!
          ISAVE(IISHR1) = -1
          ISAVE(IISHR2) = -1
C        Prepare partitions for the gaussian integration. 
          R1L  = -LOG(BRD1)/ALP1
          R1H  = -LOG(BRD2)/ALP1
          R2L  = -LOG(BRD1)/ALP2
          R2H  = -LOG(BRD2)/ALP2
          PART(0) = ZERO
          IPARTS  = 0
C        Is there a near-zero partition?
          IF( MIN(A1-R1H,A2-R2H).GT.ZERO ) THEN
              IPARTS       = IPARTS + 1
              PART(IPARTS) = MIN(A1-R1H,A2-R2H)
              IGO (IPARTS) = IOZERO
          ENDIF
C        Is there a lower border partition?
          IF( MIN(A1-R1L,A2-R2L).GT.ZERO ) THEN
              IPARTS       = IPARTS + 1
              PART(IPARTS) = MIN(A1-R1L,A2-R2L)
              IGO (IPARTS) = IOBORD
          ENDIF
C        Are the near-center partions separated?
          IF( A1+R1L.LT.A2-R2L .OR. A2+R2L.LT.A1-R1L ) THEN
              IPARTS       = IPARTS + 1
              PART(IPARTS) = MIN(A1+R1L,A2+R2L)
              IGO (IPARTS) = IOCENT
              IPARTS       = IPARTS + 1
              PART(IPARTS) = MAX(A1-R1L,A2-R2L)
              IGO (IPARTS) = IOSPAC
              IPARTS       = IPARTS + 1
              PART(IPARTS) = MAX(A1+R1L,A2+R2L)
              IGO (IPARTS) = IOCENT
          ELSE
C        If they are not, combine near-center partions and space
C        partition. 
              IPARTS       = IPARTS + 1
              PART(IPARTS) = MAX(A1+R1L,A2+R2L)
              IGO (IPARTS) = 2*IOCENT + IOSPAC
          ENDIF
C        Add final gauss partition, and the rest should
C        be covered Ok by Legendre formula.
          IPARTS       = IPARTS + 1
          PART(IPARTS) = MAX(A1+R1H,A2+R2H)
          IGO (IPARTS) = IOBORD
C        Compute total number of the integration points
          IO = IOTAIL
          DO 100 I=1,IPARTS
              IO = IO + IGO(I)
  100     CONTINUE
          IF( IO.GT.MAXPTS ) THEN
C        This should never happen...
              WRITE(NB6,11010) IO
              STOP 'PS3SCR'
          ENDIF
C        Put the number of sample points away. Also set integration limits
C        to include all of the points.
          ISAVE(IBIO ) = IO
C    Determine integration abscissas and weights
          IBAS = 0
          DO 300 I=1,IPARTS
              CALL PSIGA(IGO(I),IGO(I),PART(I-1),PART(I),
     .                   DSAVE(IBARGS+IBAS),DSAVE(IBWGTS+IBAS))
              IBAS = IBAS + IGO(I)
  300     CONTINUE
          CALL PSILGA(IOTAIL,IOTAIL,ALP1+ALP2,PART(IPARTS),
     .                DSAVE(IBARGS+IBAS),DSAVE(IBWGTS+IBAS))
C
C    Evaluate modifier function at grid points
C
          DO 500 I=0,IO-1
              DSAVE(IBWGTS+I) = DSAVE(IBWGTS+I) * F(DSAVE(IBARGS+I))
  500     CONTINUE
          ISAVE(ISIGN) = ISG
      ENDIF
C    Get the desired range of the argument points.
      IO    = ISAVE(IBIO)
C    Evaluate integrands. It is also a duty of PS3SHR to weed out
C    points with less than PLOSS relative precision.
      CALL PS3SHR(N1,L1,ALP1,LS1,A1,IO,DSAVE(IBARGS),
     .            VAL1,MAXL+1,ISAVE(IISHR1),IISHC,DSAVE(IDSHR1),IDSHC)
      CALL PS3SHR(N2,L2,ALP2,LS2,A2,IO,DSAVE(IBARGS),
     .            VAL2,MAXL+1,ISAVE(IISHR2),IISHC,DSAVE(IDSHR2),IDSHC)
C    Evaluate integrals and error estimations.
C    4XXX-5XXX loops would have been better off with
C    BLAS level 3, but repacking of both VAL1 and VAL2
C    would have been necessary in this case.
      DO 5000 M2=0,L2
          DO 4900 M1=0,L1
              SUM = ZERO
              DO 4800 I=1,IO
                  STERM = DSAVE(IBWGTS+I-1)*VAL1(M1,I)*VAL2(M2,I)
                  SUM   = SUM + STERM
 4800         CONTINUE
              AINT(M1,M2) = SUM
 4900     CONTINUE
 5000 CONTINUE
C
10010 FORMAT(' UNABLE TO SAVE ',I5,' INTEGER AND ',I5,
     .       ' DOUBLE PRECISION VARIABLES.'/' ONLY ',I5,' INTEGER AND '
     .       ,I5,' DOUBLE OPRECISION WORDS ARE AVAILABLE TO PS3SCR.')
10020 FORMAT(' CLUELESS CALLER IN PS3SCR: SIGNATURE IS: ',I9)
11010 FORMAT(' TOTAL INTEGRATION ORDER IS TOO HIGH IN PS3SCR (',I5,').')
11020 FORMAT(' WARNING: GRID IS TOO ROUGH FOR THE FIRST ',
     .       'INTEGRAND IN PS3SCR.')
11030 FORMAT(' WARNING: GRID IS TOO ROUGH FOR THE SECOND ',
     .       'INTEGRAND IN PS3SCR.')
      END
C
      SUBROUTINE PS3RMY(L,LP,ALPHA,BETA,GAMMA,D,LDG)
C
C   Compute rotation matrix for the complex spherical harmonics of a
C   given order.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L      - Order of the spherical harmonics
C      LP     - Maximum magnetic quantum number needed in the base
C               (transformed) harmonics.
C      ALPHA,
C      BETA,
C      GAMMA  - Euler angles
C      D      - (Output) complex transformation matrix for rotated->base
C                        transform.
C               First index, (-LP:LP): magnetic quantum number of a
C                   spherical harmonic in the base coordinates.
C               Second index, (-L:L): magneric quantum number of a
C                   spherical harmonic in the rotated coordinates.
C      LDG    - Leading dimension of the rotation matrix
C      
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C      8*(LMAX+.5)*(LMAX+2) DOUBLE PRECISION storage units. With
C      LMAX = 15 it's about 2100 DP words.
C
C   Module logic:
C
C      Expression is taken from M.Weissbluth, "Atoms and Molecules",
C      Academic Press, N.Y., 1978. There is a lot of overlap with
C      3J expressions (see psovx.f), but we won't exploit it. The
C      coefficient for the m->mp transformation is (I is Sqrt(-1)):
C
C                        la                             1/2
C      l+min(m,-mp)  (-1)   ((l+m)!(l-m)!(l+mp)!(l-mp)!)
C           Sum      -------------------------------------- x
C     la=max(0,m-mp) la! (l+m-la)! (l-mp-la)! (mp-m+la)!
C
C         I mp alpha              2l+m-mp-2la             mp-m+2la  I m gamma
C      x E           (Cos(beta/2))           (Sin(beta/2))         E
C
C      The code below might look somewhat different, but it
C      actully implements the expression above ;-)
C      
C
C   Possible optimizations:
C
C      Factorials and inverse factorials can be pre-computed.
C
C      Symmetry relations between magnetic quantum numbers 
C      should be used.
C
C   Bugs:
C
C      The code is not strictly standard, since the DOUBLE COMPLEX
C      type and intrinsics are used.
C
C      Maximum orbital moment is (arbitrarily) limited to LMAX (=15).
C
C      Complex angular factors can become numerically unstable for 
C      high Ls (there should be no problem with the real angular 
C      factors), but at least we should trap in this case...
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LMAX=40)
      PARAMETER (ONE =1.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (EPS =1.0D-13)
C
      DOUBLE COMPLEX D, EALP, EGAM, DCONJG, DCMPLX
      INTRINSIC DCONJG,DCMPLX
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION D(-LP:LDG-1-LP,-L:L)
      DIMENSION FACT(0:2*LMAX), RFACT(0:2*LMAX)
      DIMENSION EALP(-LMAX:LMAX), EGAM(-LMAX:LMAX)
      DIMENSION COSBET(0:2*LMAX), SINBET(0:2*LMAX)
      DIMENSION A(0:2*LMAX,0:2*LMAX), B(0:2*LMAX,0:2*LMAX)
C
      IF( L.LT.0 .OR. L.GT.LMAX ) THEN
          WRITE(NB6,11010) L
          STOP 'PS3RMY'
      ENDIF
      IF( LP.LT.0 .OR. LP.GT.L ) THEN
          WRITE(NB6,11015) LP
          STOP 'PS3RMY'
      ENDIF
C    Prepare factorials and inverse factorials
      FACT  (0) = ONE
      RFACT (0) = ONE
      DO 500 I=1,2*L
          FACT  (I) = FACT(I-1)*I
          RFACT (I) = ONE/FACT(I)
  500 CONTINUE
C    Prepare complex angular factors
      EALP  (0) = ONE
      IF( LP.GT.0 ) THEN
          EALP( 1) = DCMPLX(COS(ALPHA),SIN(ALPHA))
          EALP(-1) = DCONJG(EALP(1))
          DO 600 I=2,LP
              EALP( I) = EALP(I-1)*EALP(1)
              EALP(-I) = DCONJG(EALP(I))
  600     CONTINUE
      ENDIF
      EGAM  (0) = ONE
      IF( L.GT.0 ) THEN
          EGAM( 1) = DCMPLX(COS(GAMMA),SIN(GAMMA))
          EGAM(-1) = DCONJG(EGAM(1))
          DO 700 I=2,L
              EGAM( I) = EGAM(I-1)*EGAM(1)
              EGAM(-I) = DCONJG(EGAM(I))
  700     CONTINUE
      ENDIF
C    Trap precision loss in complex angular factors
      IF( ABS(ABS(EALP(LP))-ONE).GT.EPS .OR.
     .    ABS(ABS(EGAM(L ))-ONE).GT.EPS ) THEN
          WRITE(NB6,11020) EPS
          STOP 'PS3RMY'
      ENDIF
C    Prepare real angular factors
      COSBET(0) = ONE
      SINBET(0) = ONE
      IF( L.GT.0 ) THEN
          COSBET(1) = COS(HALF*BETA)
          SINBET(1) = SIN(HALF*BETA)
          DO 800 I=2,2*L
              COSBET(I) = COSBET(I-1)*COSBET(1)
              SINBET(I) = SINBET(I-1)*SINBET(1)
  800     CONTINUE
      ENDIF
C    Wrap factors dependent on magnetic numbers only into
C    complex angular factors
      EALP(0) = EALP(0)*FACT(L)
      EGAM(0) = EGAM(0)*FACT(L)
      DO 900 I=1,LP
          SCALE    = SQRT(FACT(L+I)*FACT(L-I))
          EALP( I) = EALP( I)*SCALE
          EALP(-I) = EALP(-I)*SCALE
          EGAM( I) = EGAM( I)*SCALE
          EGAM(-I) = EGAM(-I)*SCALE
  900 CONTINUE
      DO 1000 I=LP+1,L
          SCALE    = SQRT(FACT(L+I)*FACT(L-I))
          EGAM( I) = EGAM( I)*SCALE
          EGAM(-I) = EGAM(-I)*SCALE
 1000 CONTINUE
C    Precompute scaled real angular factors
      DO 1200 J=0,2*L
          DO 1100 I=J,2*L
              SCALE  = RFACT(J)*RFACT(I-J)
              A(I,J) = COSBET(I)*SCALE
              B(I,J) = SINBET(I)*SCALE
 1100     CONTINUE
 1200 CONTINUE
C    Put the sign factor into B
      DO 1400 J=1,2*L,2
          DO 1300 I=J,2*L
              B(I,J) = -B(I,J)
 1300     CONTINUE
 1400 CONTINUE
C    The easy part - loop over all magnetic quantum numbers
C    and compute our answers
      DO 2000 M=-L,L
          DO 1900 MP=-LP,LP
              SUM = ZERO
              DO 1800 LA=MAX(0,M-MP),L+MIN(M,-MP)
                  SUM = SUM + A(2*L+M-MP-2*LA,L+M-LA)*B(MP-M+2*LA,LA)
 1800         CONTINUE
              D(MP,M) = SUM*EALP(MP)*EGAM(M)
 1900     CONTINUE
 2000 CONTINUE
      RETURN
11010 FORMAT(' ORBITAL MOMENT ',I5,' IS NOT SUPORTED IN PS3RMY.')
11015 FORMAT(' INVALID AUXILIARY MOMENT ',I5,' PASSED TO PS3RMY.')
11020 FORMAT(' PRECISION LOSS IN PS3RMY EXCEEDS ',G10.4,' IN PS3RMY.')
      END
C
      SUBROUTINE PS3EUL(R,A,EUL)
C
C   Evaluate distance and Euler angles for a given pair of atoms.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      R      - Relative position of the *rotation center* (i.e.,
C               direction of R is the direction of the rotated Z 
C               axis using Sharma's conventions), three components.
C      A      - (Output) Distance to operator center. A is always
C               non-negative and is measured in the *negative*
C               direction on the rotated Z axis.
C      EUL    - (Output) Euler angles for *a* rotation aligning
C               Z axis with R.
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
C      Euler angles for a transformation which aligns Z axis with a
C      specified direction vector are computed. The definition of
C      the Euler angles used here is as follows:
C       a. Right rotation on Z axis by gamma=eul(3), XYZ->X'Y'Z'
C       b. Right rotation on Y' axis by beta=eul(2), X'Y'Z'->X''Y''Z''
C       c. Right rotation on Z'' axis by alpha=eul(1), X''Y''Z''->X'''Y'''Z'''
C
C      Note that Euler angles are not uniquely defined by a direction
C      vector.
C
C   Possible optimizations:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0)
      DIMENSION R(3), EUL(3)
C
      A      = SQRT(R(1)**2+R(2)**2+R(3)**2)
      IF( A.EQ.ZERO ) THEN
          EUL(1) = ZERO
          EUL(2) = ZERO
          EUL(3) = ZERO
          RETURN 
      ENDIF
      EUL(3) = ZERO
      IF( R(1).EQ.ZERO .AND. R(2).EQ.ZERO ) THEN
          EUL(1) = ZERO
      ELSE
          EUL(1) = -ATAN2(R(2),R(1))
      ENDIF
      EUL(2) = -ACOS(R(3)/A)
C
      RETURN
      END
C
      SUBROUTINE PS3YLY(L,YLY,LDY)
C                         ^
C   Compute one-center <Y|L|Y> integrals, where Y is a complex spherical
C   harmonic and L the "physical" orbital moment operator (i.e., with
C   the -i factor included).
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L      - Order of spherical harmonics
C      YLY    - (Output) integrals, three-dimensional array
C               First dimension, -L to L: magnetic quantum numbers
C                   of the left harmonic in the increasing order
C               Second dimension, 1 to 3: projections of the orbital
C                   moment operator, in the order X,Y,Z. For X and Z,
C                   real part of the integral is given. For Y, imaginary
C                   part is given.
C               Third dimension, -L to L: magnetic quantum numbers of
C                   the right harmonic in the increasing order.
C      YLY    - Leading dimension of the YLY matrix.
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
C      Nope ;-)
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      DIMENSION YLY(-L:LDY-L-1,3,-L:L)
C
      DO 400 M2=-L,L
          SUP = HALF*SQRT(DBLE(L*(L+1)-M2*(M2+1)))
          SDN = HALF*SQRT(DBLE(L*(L+1)-M2*(M2-1)))
C        Lx
          DO 100 M1=-L,L
              YLY(M1,1,M2) = ZERO
  100     CONTINUE
          IF( M2.GT.-L ) YLY(M2-1,1,M2) = SDN
          IF( M2.LT. L ) YLY(M2+1,1,M2) = SUP
C        Ly
          DO 200 M1=-L,L
              YLY(M1,2,M2) = ZERO
  200     CONTINUE
          IF( M2.GT.-L ) YLY(M2-1,2,M2) = SDN
          IF( M2.LT. L ) YLY(M2+1,2,M2) =-SUP
C        Lz
          DO 300 M1=-L,L
              YLY(M1,3,M2) = ZERO
  300     CONTINUE
          YLY(M2,3,M2) = M2
  400 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PS3YDY(L1,L2,YYY,LDY)
C
C   Compute one-center <Y|Y_{1,m}|Y> integrals.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L1     - Order of the left-hand spherical harmonics.
C      L2     - Order of the right-hand spherical harmonics.
C      YYY    - (Output) integrals, three-dimensional array
C               First dimension, -L1 to L1: magnetic quantum numbers
C                   of the left harmonic in the increasing order
C               Second dimension, -1 to 1: magnetic quantum number
C                   of the operator harmonic.
C               Third dimension, -L2 to L2: magnetic quantum numbers of
C                   the right harmonic in the increasing order.
C      LDY    - Leading dimension of the YLY matrix.
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
C      Special-cased Gaunt formula:
C
C                                                       / L1 L2  1 \
C      <Y1|Y|Y2> = (-1)^(Lmax+M1) Sqrt[ (3/4 Pi) Lmax ] |          |
C                                                       \-M1 M2 MO /
C      where Lmax = Max[ L1, L2 ] and <> is Wigner's 3J symbol.
C
C   Accuracy:
C
C      Computed integrals are good to machine accuracy if L1 and L2
C      are not too large (up to 20).
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (DMAGIC=0.2387324146378430036533256450587715430516895D0)
      EXTERNAL PS3J
      DIMENSION YYY(-L1:LDY-L1-1,-1:1,-L2:L2)
C    Zero out everything first, and think later ;-)
      DO 1000 M2=-L2,L2
          DO 900 MO=-1,1
              DO 800 M1=-L1,L1
                  YYY(M1,MO,M2) = ZERO
  800         CONTINUE
  900     CONTINUE
 1000 CONTINUE
C    Catch zero-only cases
      IF( ABS(L1-L2).NE.1 ) RETURN
C
      LMAX  = MAX(L1,L2)
      SCALE = SQRT(DMAGIC*LMAX)
      IF( MOD(LMAX,2).EQ.1 ) SCALE = -SCALE
C
      DO 2000 M2=-L2,L2
          DO 1900 MO=-1,1
              M1 = MO + M2
              IF( ABS(M1).LE.L1 ) THEN
                  TERM = SCALE * PS3J(L1,L2,1,-M1,M2,MO)
                  IF( MOD(ABS(M1),2).EQ.1 ) TERM = -TERM
                  YYY(M1,MO,M2) = TERM
              ENDIF
 1900     CONTINUE
 2000 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PS3Y1C(ICNT1,ICNT2,AINT,LDA)
C
C   Convert integrals over complex order-1 spherical harmonics
C   into cartesian components integrals.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ICNT1  - Number of left-hand orbitals
C      ICNT2  - Number of right-hand orbitals
C      AINT   - Input: integrals over Y_{1,m}, three-dimensional DOUBLE
C                      COMPLEX array
C               Output: integrals over x/r, y/r, z/r, three-dimensional array
C               First dimension, 1 to ICNT1: left-hand orbitals.
C               Second dimension, -1 to 1: On input m value of the operator.
C                   On output, x, y and z components of the integrals in this
C                   order.
C               Third dimension, 1 to ICNT2: right-hand orbitals.
C      LDA    - Leading dimension of the LDA matrix.
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
C      X/R =   C1*( Y_{1,-1} - Y_{1,1} )
C      Y/R = I*C1*( Y_{1,-1} + Y_{1,1} )
C      Z/R = C2*Y_{1,0}
C
C      C1 = Sqrt[2 Pi/3]
C      C2 = Sqrt[4 Pi/3]
C      I  = Sqrt[-1]
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0)
      PARAMETER (C1=1.4472025091165353187260292545815915357442D0)
      PARAMETER (C2=2.0466534158929769769591032497785297214148D0)
      DOUBLE COMPLEX AINT, DCMPLX, X, Y, Z
      INTRINSIC DCMPLX
      DIMENSION AINT(LDA,-1:1,ICNT2)
C
      DO 2000 M2=1,ICNT2
          DO 1000 M1=1,ICNT1
              X = C1*             ( AINT(M1,-1,M2) - AINT(M1,1,M2) )
              Y = DCMPLX(ZERO,C1)*( AINT(M1,-1,M2) + AINT(M1,1,M2) )
              Z = C2*               AINT(M1, 0,M2)
              AINT(M1,-1,M2) = X
              AINT(M1, 0,M2) =-Y
              AINT(M1, 1,M2) = Z
 1000     CONTINUE
 2000 CONTINUE
C
      END
C
      SUBROUTINE PS3RTY(LS1,LS2,NC,L1I,EULER1,L2I,EULER2,YM,LDM,YL,LDL)
C
C   Rotate one-center integrals over complex spherical harmonics
C   into double-local coordinate system.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      LS1    - Order of the left-hand spherical harmonics
C      LS2    - Order of the right-hand spherical harmonics
C      NC     - Number of integral's components
C      L1I    - Maximum desired magnetic quantum number on
C               left-hand harmonics
C      EULER1 - Euler angles for molecular->left-hand local rotation
C      L2I    - Maximum desired magnetic quantum number on
C               right-hand harmonics
C      EULER2 - Euler angles for molecular->right-hand local rotation
C      YM     - Integrals in the molecular coordinate system, three-
C               dimentional *DOUBLE PRECISION* array:
C               First index, -LS1 to LS1: left-hand magnetic quantum number.
C               Second index, 1 to NC: components of the operator.
C               Third index, -LS2 to LS2: right-hand magnetic quantum number.
C      LDM    - Leading dimension of the YM matrix.
C      YL     - Integrals in the double-local coordinate system, three-
C               dimensional *DOUBLE COMPLEX* array:
C               First index, -L1 to L1: left-hand magnetic quantum number.
C               Second index, 1 to NC: component of the operator.
C               Third index, -L2 to L2: right-hand magnetic quantum number.
C      LDL    - Leading dimension of the YL matrix.
C      
C   Accessed common blocks:
C
C   Local storage:
C
C      6*(2*MAXL+1)*(2*MAXLS+1) DOUBLE PRECISION units. ca. 1300 DP
C      words with MAXL=3 and MAXLS=15
C
C   Module logic:
C
C      Two-stage transformation is used, first on left-hand side,
C      than on right-hand side.
C
C   Possible optimizations:
C
C      Using BLAS3 should increase performance for large LS, L1 and
C      L2 values. With present limits on LS, L1 and L2 this is unlikely
C      to help unless data cache is *very* small (say less than 4Kb).
C
C   Bugs:
C
C      The code is not strictly standard, since the DOUBLE COMPLEX
C      type and intrinsics are used.
C
C      Maximum order of spherical harmonics is fixed at MAXLS (=15).
C      Maximum magnetic quantum number of the double-local harmonics
C      is fixed at MAXL (=3).
C
C      Decent performance can only be expected if everything fits
C      into cache.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=3)
      PARAMETER (MAXLS=40)
      PARAMETER (ZERO=0.D0)
C
      DOUBLE COMPLEX RM1, RM2, YH, YL, SUM, DCONJG
      INTRINSIC DCONJG
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION EULER1(3), EULER2(3)
      DIMENSION YM(-LS1:LDM-LS1-1,NC,-LS2:LS2)
      DIMENSION YL(-L1I:LDL-L1I-1,NC,-L2I:L2I)
      DIMENSION RM1(-MAXL:MAXL,-MAXLS:MAXLS)
      DIMENSION RM2(-MAXL:MAXL,-MAXLS:MAXLS)
      DIMENSION YH (-MAXL:MAXL,-MAXLS:MAXLS)
C    Rangle checking...
      IF( LS1.LT.0 .OR. LS1.GT.MAXLS ) THEN
          WRITE(NB6,11010) LS1
          STOP 'PS3RTY'
      ENDIF
      IF( LS2.LT.0 .OR. LS2.GT.MAXLS ) THEN
          WRITE(NB6,11015) LS2
          STOP 'PS3RTY'
      ENDIF
      IF( L1I.LT.0 .OR. L1I.GT.MAXL .OR. 
     .    L2I.LT.0 .OR. L2I.GT.MAXL ) THEN
          WRITE(NB6,11020) L1I, L2I
          STOP 'PS3RTY'
      ENDIF
C    There will be no terms with ABS(M1)>LS, likewise for M2
      L1 = MIN(LS1,L1I)
      L2 = MIN(LS2,L2I)
      DO 500 I=1,NC
          DO 450 M1=L1+1,L1I
              DO 440 M2=-L2I,L2I
                  YL( M1,I,M2) = ZERO
                  YL(-M1,I,M2) = ZERO
  440         CONTINUE
  450     CONTINUE
          DO 490 M2=L2+1,L2I
              DO 480 M1=-L1,L1
                  YL(M1,I, M2) = ZERO
                  YL(M1,I,-M2) = ZERO
  480         CONTINUE
  490     CONTINUE
  500 CONTINUE
C    Compute rotation matrices
      CALL PS3RMY(LS1,L1,-EULER1(3),-EULER1(2),-EULER1(1),
     .            RM1(-L1,-LS1),2*MAXL+1)
      CALL PS3RMY(LS2,L2,-EULER2(3),-EULER2(2),-EULER2(1),
     .            RM2(-L2,-LS2),2*MAXL+1)
C    Loop over all operator components
      DO 3000 I=1,NC
C        First half-transformation. This is GEMM-type operation, actually,
C        but BLAS3 does not provide mixed-type routines, unfortunately.
          DO 1000 M2=-LS2,LS2
              DO 900 M1=-L1,L1
                  SUM = ZERO
                  DO 800 MP=-LS1,LS1
                      SUM = SUM + DCONJG(RM1(M1,MP))*YM(MP,I,M2)
  800             CONTINUE
                  YH(M1,M2) = SUM
  900         CONTINUE
 1000     CONTINUE
C        Second half-transformation. We should use BLAS3 here, but it'll
C        mean supporting complex routines, which I do not want to do...
          DO 2000 M2=-L2,L2
              DO 1900 M1=-L1,L1
                  SUM = ZERO
                  DO 1800 MP=-LS2,LS2
                      SUM = SUM + RM2(M2,MP)*YH(M1,MP)
 1800             CONTINUE
                  YL(M1,I,M2) = SUM
 1900         CONTINUE
 2000     CONTINUE
 3000 CONTINUE
C
      RETURN
11010 FORMAT(' INVALID LS1 (',I5,') PASSED TO PS3RTY.')
11015 FORMAT(' INVALID LS2 (',I5,') PASSED TO PS3RTY.')
11020 FORMAT(' INVALID L1 (',I5,') OR L2 (',I5,') PASSED TO PS3RTY.')
      END
C
      SUBROUTINE PS3RTM(L1,EUL1,L2,EUL2,NC,AL,LDL,AM,LDM)
C
C   Rotate integrals from double-local coordinate system into molecular
C   coordinate system and revert to real spherical harmonics.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L1     - Orbital quantum number of left-hand functions.
C      EUL1   - Euler angles for the molecular->left local transformation.
C      L2     - Orbital quantum number of right-hand functions.
C      EUL2   - Euler angles for the molecular->right local transformation.
C      NC     - Number of the operator compoments.
C      AL     - Integrals over complex functions in the double-local system.
C               First index, -L1 to L1: magnetic quantum number of the
C                   left-hand orbitals.
C               Second index, 1 to NC: operator components.
C               Third index, -L2 to L2: magnetic quantum number of the
C                   right-hand orbitals.
C      LDL    - Leading dimension of the AL matrix.
C      AM     - (Output) Integrals over real Slater AOs in the molecular system.
C               First index, -L1 to L1: magnetic quantum number of the left-hand
C                   orbitals
C               Second index, 1 to NC: operator components.
C               Third index, -L2 to L2: magnetic quantum number of the
C                   right-hand orbitals.
C      LDM    - Leading dimension of the AM matrix.
C      
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C      8*(2*MAXL+1)**2 DOUBLE PRECISION words, i.e. ca. 400 DP
C      words with MAXL=3.
C
C   Module logic:
C
C   Possible optimizations:
C
C      If orbital moments were larger, it might have made sence
C      to use BLAS3 to do rotations. They are not, and it doesn't, 
C      though.
C
C   Bugs:
C
C      The code is not strictly standard, since the DOUBLE COMPLEX
C      type and intrinsics are used.
C
C      Orbital moment of the orbitals is fixed at MAXL (=3).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=3)
      PARAMETER (ZERO=0.D0)
      PARAMETER (EPS =1.D-7)
      PARAMETER (RSQ2=0.70710678118654752440084436210484903928483594D0)
      DOUBLE COMPLEX AL, RM1, RM2, HT, FT, SUM, DIF, DCMPLX, DCONJG
      INTRINSIC DIMAG, DCMPLX, DCONJG
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION AL(-L1:LDL-L1-1,NC,-L2:L2)
      DIMENSION AM(-L1:LDM-L1-1,NC,-L2:L2)
      DIMENSION EUL1(3), EUL2(3)
      DIMENSION RM1(-MAXL:MAXL,-MAXL:MAXL), RM2(-MAXL:MAXL,-MAXL:MAXL)
      DIMENSION HT(-MAXL:MAXL,-MAXL:MAXL), FT(-MAXL:MAXL,-MAXL:MAXL)
C    Range checking...
      IF( L1.LT.0 .OR. L1.GT.MAXL .OR. L2.LT.0 .OR. L2.GT.MAXL ) THEN
          WRITE(NB6,11010) L1, L2
          STOP 'PS3RTM'
      ENDIF
C    Evaluate rotation matrices
      CALL PS3RMY(L1,L1,EUL1(1),EUL1(2),EUL1(3),RM1(-L1,-L1),2*MAXL+1)
      CALL PS3RMY(L2,L2,EUL2(1),EUL2(2),EUL2(3),RM2(-L2,-L2),2*MAXL+1)
C    Transform operator components one by one
      DO 6000 I=1,NC
C        Rotate left orbital ("ZGEMM")
          DO 1000 M2=-L2,L2
              DO 900 M1=-L1,L1
                  SUM = ZERO
                  DO 800 MC=-L1,L1
                      SUM = SUM + DCONJG(RM1(M1,MC))*AL(MC,I,M2)
  800             CONTINUE
                  HT(M1,M2) = SUM
  900         CONTINUE
 1000     CONTINUE
C        Rotate right orbital ("ZGEMM" once again)
          DO 2000 M2=-L2,L2
              DO 1900 M1=-L1,L1
                  SUM = ZERO
                  DO 1800 MC=-L2,L2
                      SUM = SUM + RM2(M2,MC)*HT(M1,MC)
 1800             CONTINUE
                  FT(M1,M2) = SUM
 1900         CONTINUE
 2000     CONTINUE
C        Convert left orbital to real harmonics. Hote
C        that the expressions here are different from
C        ones used for right orbitals due to the complex
C        conjugation of the left integrand.
          DO 3000 M2=-L2,L2
              DO 2800 M1=1,L1,2
                  SUM = RSQ2*(FT(-M1,M2)+FT(M1,M2))
                  DIF = RSQ2*(FT(-M1,M2)-FT(M1,M2))
                  FT(-M1,M2) = DCMPLX(-DIMAG(SUM),DBLE(SUM))
                  FT( M1,M2) = DIF
 2800         CONTINUE
              DO 2900 M1=2,L1,2
                  SUM = RSQ2*(FT(-M1,M2)+FT(M1,M2))
                  DIF = RSQ2*(FT(-M1,M2)-FT(M1,M2))
                  FT(-M1,M2) = DCMPLX(-DIMAG(DIF),DBLE(DIF))
                  FT( M1,M2) = SUM
 2900         CONTINUE
 3000     CONTINUE
C        Convert right orbital to real harmonics
          DO 4000 M1=-L1,L1
              DO 3800 M2=1,L2,2
                  SUM = RSQ2*(FT(M1,-M2)+FT(M1,M2))
                  DIF = RSQ2*(FT(M1,-M2)-FT(M1,M2))
                  FT(M1,-M2) = DCMPLX(DIMAG(SUM),-DBLE(SUM))
                  FT(M1, M2) = DIF
 3800         CONTINUE
              DO 3900 M2=2,L2,2
                  SUM = RSQ2*(FT(M1,-M2)+FT(M1,M2))
                  DIF = RSQ2*(FT(M1,-M2)-FT(M1,M2))
                  FT(M1,-M2) = DCMPLX(DIMAG(DIF),-DBLE(DIF))
                  FT(M1, M2) = SUM
 3900         CONTINUE
 4000     CONTINUE
C        Stuff integrals into output buffer, making sure they
C        are really real!
          DO 5000 M2=-L2,L2
              DO 4900 M1=-L1,L1
                  IF( ABS(DIMAG(FT(M1,M2))).GT.EPS ) THEN
                      WRITE(NB6,11020) I, M1, M2, FT(M1,M2)
                      STOP 'PS3RTM'
                  ENDIF
                  AM(M1,I,M2) = DBLE(FT(M1,M2))
 4900         CONTINUE
 5000     CONTINUE
 6000 CONTINUE
C
      RETURN
11010 FORMAT(' EITHER L1 (',I5,') OR L2 (',I5,') IS INVALID IN PS3RTM.')
11020 FORMAT(1X,I3,'-TH COMPONENT OF THE (',I3,',',I3,') INTEGRAL ',
     .       'REMAINED COMPLEX AFTER REDUCTION IN PS3RTM:'/1X,2G20.12)
      END
C
      DOUBLE PRECISION FUNCTION PS3AIA(X0,X1,X2)
C
C   Compute Aitken's approximation to a series limit.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      X0,X1,X2 - Three successive iterates of the series
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
C      Series is assumed to converge approximately geometrically,
C      and Aitken's correction is evaluated according to expr 3.9.7
C      in Abramowitz & Stegun.
C
C   Possible optimizations:
C
C   Bugs:
C
      DOUBLE PRECISION X0, X1, X2, D2X, DELTA, EPS
      PARAMETER (EPS=1D-2)
      D2X    = X0 - 2*X1 + X2
      DELTA  = (X0-X1)**2
      PS3AIA = X2
      IF( EPS*ABS(X0*D2X).GT.DELTA ) PS3AIA = X0 - DELTA/D2X
      RETURN
      END
C
      SUBROUTINE PS3ACA(ITER,L1,NC,L2,AX,ACCEL,MAXL,MAXC)
C
C   Store integrals for later extrapolation.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ITER   - Iteration number. Used to avoid copying uninitialized
C               data.
C      L1     - Orbital quantum number of the left-hand function.
C      NC     - Number of components in the integral
C      L2     - Orbital quantum number of the right-hand function.
C      AX     - Most recentry evaluated integrals. Three-index array.
C               First index: -L1 to L1, left-hand function
C               Second index: 1 to NC, integral's components
C               Third index: -L2 to L2, right-hand function
C      ACCEL  - Storage area, four-index array.
C               First index: -L1 to L1, left-hand function
C               Second index: 1 to NC, integral's components
C               Third index: -L2 to L2, right-hand function
C               Fourth index: 1 to 3, last three series sums,
C                             1 is the last sum, 2 is next to last, etc.
C      MAXL   - Maximum orbital quantum number allowed
C      MAXC   - Maximum number of integral components allowed
C      
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Accuracy:
C
C   Possible optimizations:
C
C   Bugs:
C
C      None. "Used before set" warning issued by ftnchek is in error.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX AX, ACCEL
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION AX   (-MAXL:MAXL,MAXC,-MAXL:MAXL)
      DIMENSION ACCEL(-MAXL:MAXL,MAXC,-MAXL:MAXL,3)
C    A litle check is won't harm a good girl...
      IF( L1.LT.0 .OR. L1.GT.MAXL .OR. L2.LT.0 .OR. L2.GT.MAXL 
     .                            .OR. NC.LT.0 .OR. NC.GT.MAXC ) THEN
          WRITE(NB6,11000) L1, L2, NC, MAXL, MAXC
          STOP 'PS3ACA'
      ENDIF
C    Move old enries back a bit
      DO 1000 I=MIN(ITER-1,2),1,-1
          DO 900 M2=-L2,L2
              DO 800 IC=1,NC
                  DO 700 M1=-L1,L1
                      ACCEL(M1,IC,M2,I+1) = ACCEL(M1,IC,M2,I)
  700             CONTINUE
  800         CONTINUE
  900     CONTINUE
 1000 CONTINUE
C    Put in newcomers
      DO 1900 M2=-L2,L2
          DO 1800 IC=1,NC
              DO 1700 M1=-L1,L1
                  ACCEL(M1,IC,M2,1) = AX(M1,IC,M2)
 1700         CONTINUE
 1800     CONTINUE
 1900 CONTINUE
C
      RETURN
11000 FORMAT(' BAD PS3ACA PARAMETERS: L1 = ',I4,' L2 = ',I4,
     .       ' NC = ',I4,' MAXL = ',I4,' MAXC = ',I4)
      END
C
      SUBROUTINE PS3ACL(L1,NC,L2,AX,ACCEL,MAXL,MAXC)
C
C   Extrapolate three last terms in the series.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L1     - Orbital quantum number of the left-hand function.
C      NC     - Number of components in the integral
C      L2     - Orbital quantum number of the right-hand function.
C      AX     - (Output) Extrapolated integrals. Three-index array.
C               First index: -L1 to L1, left-hand function
C               Second index: 1 to NC, integral's components
C               Third index: -L2 to L2, right-hand function
C      ACCEL  - Storage area, four-index array.
C               First index: -L1 to L1, left-hand function
C               Second index: 1 to NC, integral's components
C               Third index: -L2 to L2, right-hand function
C               Fourth index: 1 to 3, last three series sums,
C                             1 is the last sum, 2 is next to last, etc.
C      MAXL   - Maximum orbital quantum number allowed
C      MAXC   - Maximum number of integral components allowed
C      
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C   Accuracy:
C
C   Possible optimizations:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX AX, ACCEL, X0, X1, X2
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION AX   (-MAXL:MAXL,MAXC,-MAXL:MAXL)
      DIMENSION ACCEL(-MAXL:MAXL,MAXC,-MAXL:MAXL,3)
C
      INTRINSIC DIMAG, DCMPLX
      EXTERNAL PS3AIA
C    A litle check is won't harm a good girl...
      IF( L1.LT.0 .OR. L1.GT.MAXL .OR. L2.LT.0 .OR. L2.GT.MAXL 
     .                            .OR. NC.LT.0 .OR. NC.GT.MAXC ) THEN
          WRITE(NB6,11000) L1, L2, NC, MAXL, MAXC
          STOP 'PS3ACL'
      ENDIF
C
      DO 1000 M2=-L2,L2
          DO 900 IC=1,NC
              DO 800 M1=-L1,L1
                  X0 = ACCEL(M1,IC,M2,3)
                  X1 = ACCEL(M1,IC,M2,2)
                  X2 = ACCEL(M1,IC,M2,1)
                  AX(M1,IC,M2) = DCMPLX(
     .                           PS3AIA(DBLE(X0),DBLE(X1),DBLE(X2)),
     .                           PS3AIA(DIMAG(X0),DIMAG(X1),DIMAG(X2)))
  800         CONTINUE
  900     CONTINUE
 1000 CONTINUE
C
      RETURN
11000 FORMAT(' BAD PS3ACL PARAMETERS: L1 = ',I4,' L2 = ',I4,
     .       ' NC = ',I4,' MAXL = ',I4,' MAXC = ',I4)
      END
C
      SUBROUTINE PS3LR3(N1,L1,ALP1,R1,N2,L2,ALP2,R2,AINT,LDA)
C                                                           ^
C   Compute three-center integrals of the type <\phi_1|r^-3 L|\phi_2>
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N1     - Major quantum number of the STO block on left-hand center
C      L1     - Orbital quantum number of the STO block on left-hand center
C      ALP1   - Orbital exponent on left-hand center
C      R1     - Array of X,Y,Z coordinates of left-hand center relative
C               to the operator center.
C      N2,L2,ALP2,R2
C             - Attributes of the STO block on right-hand center
C      AINT   - (Output) integrals over STOs, three-index array
C               First index: -L1 to L1, magnetic quantum number of left-hand
C                            orbitals in the increasing order.
C               Second index: 1 to 3, projections of the orbital moment,
C                            in the X,Y,Z order.
C               Third index: -L2 to L2, magnetic quantum number of right-hand
C                            orbitals in the increasing order.
C      LDA    - Leading dimension of the AINT matrix.
C      
C   Accessed common blocks:
C
C      PSO3PR - Contains desired precision of the integrals group.
C      PSPRT  - Printing unit.
C
C   Local storage:
C
C      49*MAXL*(MAXL+1)+MAXL+12*MAXLS*(MAXLS+1) DOUBLE PRECISION storage
C      units. With MAXL=3 and MAXLS=15, it comes out as some 3500 DP words.
C      Add LDSAVE (=9200) DP words and LISAVE (=110) I words to that, so
C      that we can provide persistent storage to radial integration routines,
C      and it becomes quite sizeable at some 100K bytes.
C
C   Module logic:
C
C      Both left- and right- hand orbitals are expanded in terms of Sharma's
C      functions at the operator center. Each expansion term is in turn 
C      separated into the radial and angular parts which are then evaluated
C      separately.
C
C      See comments in the source and Sharma's paper on the expansion functions
C      for a (little) more details.
C
C   Accuracy:
C
C      10**-3 a.u. absolute, although in cases of the equal distances to the
C      operator center it might go up to 2*10**-3 a.u. It might be possible 
C      to do 10**-5 a.u. for strictly three-center integrals, but that's about it.
C
C      Alghough we do prepare for Aitken's delta-square convergence
C      acceleration, it seems to hurt more than it helps, so we don't
C      actually do it.
C
C      Accuracy of the code is limited by the accuracy of integrals over
C      Sharma's expansions and thus, ultimately, by the numerical 
C      instabilities in the expansion.
C
C   Possible optimizations:
C
C   Bugs:
C
C      The code is not strictly standard, since the DOUBLE COMPLEX
C      type and intrinsics are used.
C
C      Orbital moment of the orbitals is fixed at MAXL (=3).
C      Maximum expansion order is fixed at MAXLS (=15).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL = 3)
      PARAMETER (MAXLS=40)
      PARAMETER (ZERO =0.D0)
C
      PARAMETER (LISAVE=185)
      PARAMETER (LDSAVE=16900)
C
      COMMON
     ./PSO3PR/ EPSL, EPSR
     ./PSPRT / NB6
      SAVE /PSO3PR/, /PSPRT /
C
      EXTERNAL PSSNRM, PS3RM3, PS3AIA, PSO3BD
      DOUBLE PRECISION PS3RM3, PSSNRM
      DOUBLE COMPLEX YLYL, DCMPLX, ALINT, TERM, ACCEL
      LOGICAL CONV0, CONV1, CONV2, CONV3
      INTRINSIC DCMPLX, DBLE, DIMAG
      DIMENSION R1(3), R2(3), AINT(-L1:LDA-L1-1,3,-L2:L2)
      DIMENSION EULER1(3), EULER2(3)
      DIMENSION SHINT(0:MAXL,0:MAXL)
      DIMENSION ALINT (-MAXL:MAXL,3,-MAXL:MAXL)
      DIMENSION ACCEL (-MAXL:MAXL,3,-MAXL:MAXL,3)
      DIMENSION YLYM(-MAXLS:MAXLS,3,-MAXLS:MAXLS)
      DIMENSION YLYL(-MAXL:MAXL,3,-MAXL:MAXL)
      DIMENSION ISAVE(LISAVE), DSAVE(LDSAVE)
C    Range checking
      IF( L1.LT.0 .OR. L1.GT.MAXL .OR. L2.LT.0 .OR. L2.GT.MAXL ) THEN
          WRITE(NB6,11010) L1, L2
          STOP 'PS3LR3'
      ENDIF
C    Prepare scaling factors, inter-center distances and euler angles
      SCALE = PSSNRM(N1,ALP1)*PSSNRM(N2,ALP2)
      AEPS  = EPSL/SCALE
      CALL PS3EUL(R1,A1,EULER1)
      CALL PS3EUL(R2,A2,EULER2)
C    Zero out auxiliary integrals in the double-local coordinate system
      DO 400 M1=-L1,L1
          DO 390 I=1,3
              DO 380 M2=-L2,L2
                  ALINT(M1,I,M2) = ZERO
  380         CONTINUE
  390     CONTINUE
  400 CONTINUE
      CONV3 = .FALSE.
      CONV2 = .FALSE.
      CONV1 = .FALSE.
      CONV0 = .FALSE.
C    Initialize persistent store
      ISAVE(1) = -1
C    Loop over expansion orders until integrals converge
      DO 7000 LS=1,MAXLS
C        Evaluate the radial part of the expansion integrals
          CALL PS3SCR(LS,LS,PS3RM3,N1,L1,ALP1,A1,N2,L2,ALP2,A2,
     .                SHINT,MAXL+1,ISAVE,LISAVE,DSAVE,LDSAVE)
C        Evaluate angular part in the molecular coordinate system
          CALL PS3YLY(LS,YLYM(-LS,1,-LS),2*MAXLS+1)
C        Rotate angular part into double-local coordinate system
          CALL PS3RTY(LS,LS,3,L1,EULER1,L2,EULER2,YLYM(-LS,1,-LS),
     .                2*MAXLS+1,YLYL(-L1,1,-L2),2*MAXL+1)
C        Compensate Ly component of the rotated angular part, which was
C        purely imaginary in the first place. At the same time, the
C        operator we are interested in differs by a factor (0,-1) from 
C        that we had computed in PS3YLY, so that we have to multiply
C        X and Z components by (0,-1) and leave Y component alone.
C        If you are still following, you seem to be pretty sharp! ;-)
          DO 1000 M2=-L2,L2
              DO 900 M1=-L1,L1
                  YLYL(M1,1,M2) = DCMPLX(DIMAG(YLYL(M1,1,M2)),
     .                                   -DBLE(YLYL(M1,1,M2)))
                  YLYL(M1,2,M2) = -YLYL(M1,2,M2)
                  YLYL(M1,3,M2) = DCMPLX(DIMAG(YLYL(M1,3,M2)),
     .                                   -DBLE(YLYL(M1,3,M2)))
  900         CONTINUE
 1000     CONTINUE
C        Combine the radial and angular parts together and update
C        total integral in double-local system.
          CONV3 = CONV2
          CONV2 = CONV1
          CONV1 = CONV0
          CONV0 = .TRUE.
          DO 2000 M2=-L2,L2
              DO 1900 I=1,3
                  DO 1800 M1=-L1,L1
                      TERM = SHINT(ABS(M1),ABS(M2))*YLYL(M1,I,M2)
                      ALINT (M1,I,M2) = ALINT (M1,I,M2) + TERM
                      IF( ABS(TERM).GT.AEPS ) CONV0 = .FALSE.
 1800             CONTINUE
 1900         CONTINUE
 2000     CONTINUE
          CALL PS3ACA(LS,L1,3,L2,ALINT,ACCEL,MAXL,3)
C        Check convergence if we already had a chance to update all
C        components.
          IF( LS.GT.MAX(L1,L2)+1 .AND. CONV0 .AND. CONV1
     .                           .AND. CONV2 .AND. CONV3 ) GOTO 8000
 7000 CONTINUE
C    Failure!
      WRITE(NB6,11050) EPSL, MAXLS
      WRITE(NB6,11052) N1, L1, ALP1, R1, N2, L2, ALP2, R2
*     STOP 'PS3LR3'
C    Convergence achieved in double-local coordinates
 8000 CONTINUE
C    Compute extrapolated integrals ...
      CALL PS3ACL(L1,3,L2,ALINT,ACCEL,MAXL,3)
C    ... revert to the molecular coordinates and real harmonics
      CALL PS3RTM(L1,EULER1,L2,EULER2,3,ALINT(-L1,1,-L2),2*MAXL+1,
     .                                  AINT(-L1,1,-L2),LDA)
C    Now, normalize integrals...
      DO 9000 M2=-L2,L2
          DO 8900 I=1,3
              DO 8800 M1=-L1,L1
                  AINT(M1,I,M2) = SCALE*AINT(M1,I,M2)
 8800         CONTINUE
 8900     CONTINUE
 9000 CONTINUE
C    ... and away we go!
      RETURN
11010 FORMAT(' EITHER L1 (',I5,') OR L2 (',I5,') IS INVALID IN PS3LR3.')
11050 FORMAT(' UNABLE TO CONVERGE INTEGRALS TO ',E10.4,' USING ',I3,
     .       ' EXPANSION TERMS IN PS3LR3.')
11052 FORMAT(' N1 = ',I3,' L1 = ',I2,' ALP1 = ',F12.7,' R1 = ',3F15.7/
     .       ' N2 = ',I3,' L2 = ',I2,' ALP2 = ',F12.7,' R2 = ',3F15.7)
      END
C
      DOUBLE PRECISION FUNCTION PS3RM3(X)
      DOUBLE PRECISION X
      PS3RM3 = X**(-3)
      RETURN
      END
C
      SUBROUTINE PS3PR3(N1,L1,ALP1,R1,N2,L2,ALP2,R2,AINT,LDA)
C
C   Compute three-center integrals of the type <\phi_1|{x,y,z}/r^3|\phi_2>
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N1     - Major quantum number of the STO block on left-hand center
C      L1     - Orbital quantum number of the STO block on left-hand center
C      ALP1   - Orbital exponent on left-hand center
C      R1     - Array of X,Y,Z coordinates of left-hand center relative
C               to the operator center.
C      N2,L2,ALP2,R2
C             - Attributes of the STO block on right-hand center
C      AINT   - (Output) integrals over STOs, three-index array
C               First index: -L1 to L1, magnetic quantum number of left-hand
C                            orbitals in the increasing order.
C               Second index: 1 to 3, operator harmonics in the X,Y,Z order.
C               Third index: -L2 to L2, magnetic quantum number of right-hand
C                            orbitals in the increasing order.
C      LDA    - Leading dimension of the AINT matrix.
C      
C   Accessed common blocks:
C
C      PSO3PR - Contains desired precision of the integrals group.
C      PSPRT  - Printing unit.
C
C   Local storage:
C
C      12*MAXLS*(MAXLS+1) + 49*MAXL*(MAXL+1) + MAXL + 22 DOUBLE PRECISION
C      storage units, i.e. ca. 2800 DP words with MAXL=3 and MAXLS=14.
C      Add LDSAVE (=9200) DP words and LISAVE (=110) I words to that, so
C      that we can provide persistent storage to radial integration routines,
C      and it becomes quite sizeable at some 100K bytes.
C
C   Module logic:
C
C      Both left- and right- hand orbitals are expanded in terms of Sharma's
C      functions at the operator center. Each expansion term is in turn 
C      separated into the radial and angular parts which are then evaluated
C      separately.
C
C      The code is a modification of the PS3LR3, two major differences being:
C      1. Radial part of the operator function has r^-2 rather than r^-3
C         dependence.
C      2. Expansion orders for the left- and right- hand side orbitals should
C         differ by one for angular part to be non-zero.
C
C      Up to the final transformation, we are evaluating integrals over
C      Y_{1,m} r^-2 rather than x r^-2.
C
C   Accuracy:
C
C      10**-3 a.u. absolute in three-center cases. Might be as bad as
C      2*10**-3 a.u. in cases of equal distances to the operator center.
C
C   Possible optimizations:
C
C      Reusing radial function values rather than recomputing them when the
C      expansion order is encountered more than once might increase
C      performance by a factor of two, provided that radial function
C      evaluation is a time-limiting step.
C
C   Bugs:
C
C      The code is not strictly standard, since the DOUBLE COMPLEX
C      type and intrinsics are used.
C
C      Orbital moment of the orbitals is fixed at MAXL (=3).
C      Maximum expansion order is fixed at MAXLS (=14).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL = 3)
      PARAMETER (MAXLS=40)
      PARAMETER (ZERO =0.D0)
C
      PARAMETER (LISAVE=185)
      PARAMETER (LDSAVE=16900)
C
      COMMON
     ./PSO3PR/ EPSL, EPSR
     ./PSPRT / NB6
      SAVE /PSO3PR/, /PSPRT / 
C
      EXTERNAL PSSNRM, PS3RM2, PS3AIA, PSO3BD
      DOUBLE PRECISION PS3RM2, PSSNRM
      DOUBLE COMPLEX YDYL, DCMPLX, ALINT, TERM, ACCEL
      LOGICAL CONV0, CONV1, CONV2, CONV3
      INTRINSIC DCMPLX, DBLE, DIMAG
      DIMENSION R1(3), R2(3), AINT(-L1:LDA-L1-1,3,-L2:L2)
      DIMENSION EULER1(3), EULER2(3)
      DIMENSION SHINT(0:MAXL,0:MAXL)
      DIMENSION ALINT (-MAXL:MAXL,3,-MAXL:MAXL)
      DIMENSION ACCEL (-MAXL:MAXL,3,-MAXL:MAXL,3)
      DIMENSION YDYM(-MAXLS:MAXLS,3,-MAXLS:MAXLS)
      DIMENSION YDYL(-MAXL:MAXL,3,-MAXL:MAXL)
      DIMENSION ISAVE(LISAVE), DSAVE(LDSAVE)
C    Range checking
      IF( L1.LT.0 .OR. L1.GT.MAXL .OR. L2.LT.0 .OR. L2.GT.MAXL ) THEN
          WRITE(NB6,11010) L1, L2
          STOP 'PS3PR3'
      ENDIF
C    Prepare scaling factors, inter-center distances and euler angles
      SCALE = PSSNRM(N1,ALP1)*PSSNRM(N2,ALP2)
      AEPS  = EPSR/SCALE
      CALL PS3EUL(R1,A1,EULER1)
      CALL PS3EUL(R2,A2,EULER2)
C    Zero out auxiliary integrals in the double-local coordinate system
      DO 400 M1=-L1,L1
          DO 390 I=1,3
              DO 380 M2=-L2,L2
                  ALINT(M1,I,M2) = ZERO
  380         CONTINUE
  390     CONTINUE
  400 CONTINUE
      CONV3 = .FALSE.
      CONV2 = .FALSE.
      CONV1 = .FALSE.
      CONV0 = .FALSE.
C    Initialize persistent store
      ISAVE(1) = -1
C    Loop over expansion orders until integrals converge
      DO 7000 LS2=0,MAXLS
          CONV3 = CONV2
          CONV2 = CONV1
          CONV1 = CONV0
          CONV0 = .TRUE.
          DO 6900 LS1=ABS(LS2-1),MIN(LS2+1,MAXLS),2
C            Evaluate the radial part of the expansion integrals
              CALL PS3SCR(LS1,LS2,PS3RM2,N1,L1,ALP1,A1,N2,L2,ALP2,A2,
     .                           SHINT,MAXL+1,ISAVE,LISAVE,DSAVE,LDSAVE)
C            Evaluate angular part in the molecular coordinate system
              CALL PS3YDY(LS1,LS2,YDYM(-LS1,1,-LS2),2*MAXLS+1)
C            Rotate angular part into double-local coordinate system
              CALL PS3RTY(LS1,LS2,3,L1,EULER1,L2,EULER2,
     .                    YDYM(-LS1,1,-LS2),2*MAXLS+1,
     .                    YDYL(-L1,1,-L2),2*MAXL+1)
C            Combine the radial and angular parts together and update
C            total integral in double-local system.
              DO 2000 M2=-L2,L2
                  DO 1900 I=1,3
                      DO 1800 M1=-L1,L1
                          TERM = SHINT(ABS(M1),ABS(M2))*YDYL(M1,I,M2)
                          ALINT (M1,I,M2) = ALINT (M1,I,M2) + TERM
                          IF( ABS(TERM).GT.AEPS ) CONV0 = .FALSE.
 1800                 CONTINUE
 1900             CONTINUE
 2000         CONTINUE
 6900     CONTINUE
          CALL PS3ACA(LS2+1,L1,3,L2,ALINT,ACCEL,MAXL,3)
C        Check convergence if we already had a chance to update all
C        components.
          IF( LS2.GT.MAX(L1,L2)+1 .AND. CONV0 .AND. CONV1
     .                            .AND. CONV2 .AND. CONV3 ) GOTO 8000
 7000 CONTINUE
C    Failure!
      WRITE(NB6,11050) EPSR, MAXLS
      WRITE(NB6,11052) N1, L1, ALP1, R1, N2, L2, ALP2, R2
*     STOP 'PS3PR3'
C    Convergence achieved in double-local coordinates
 8000 CONTINUE
C    Compute extrapolated integrals
      CALL PS3ACL(L1,3,L2,ALINT,ACCEL,MAXL,3)
C    .... convert integrals over complex operator into real thing
      CALL PS3Y1C(2*L1+1,2*L2+1,ALINT(-L1,1,-L2),2*MAXL+1)
C    Second, revert to the molecular coordinates and real harmonics
      CALL PS3RTM(L1,EULER1,L2,EULER2,3,ALINT(-L1,1,-L2),2*MAXL+1,
     .                                  AINT(-L1,1,-L2),LDA)
C    Now, normalize integrals...
      DO 9000 M2=-L2,L2
          DO 8900 I=1,3
              DO 8800 M1=-L1,L1
                  AINT(M1,I,M2) = SCALE*AINT(M1,I,M2)
 8800         CONTINUE
 8900     CONTINUE
 9000 CONTINUE
C    ... and away we go!
      RETURN
11010 FORMAT(' EITHER L1 (',I5,') OR L2 (',I5,') IS INVALID IN PS3PR3.')
11050 FORMAT(' UNABLE TO CONVERGE INTEGRALS TO ',E10.4,' USING ',I3,
     .       ' EXPANSION TERMS IN PS3PR3.')
11052 FORMAT(' N1 = ',I3,' L1 = ',I2,' ALP1 = ',F12.7,' R1 = ',3F15.7/
     .       ' N2 = ',I3,' L2 = ',I2,' ALP2 = ',F12.7,' R2 = ',3F15.7)
      END
C
      DOUBLE PRECISION FUNCTION PS3RM2(X)
      DOUBLE PRECISION X
      PS3RM2 = X**(-2)
      RETURN
      END
C
      BLOCK DATA PSO3BD
C
C   Define default precision parameters for 3-center integrals.
C
C   Coded by: Serge Pachkovsky
C
C   Defined common blocks:
C
C      PSO3PR - Precision control for three-center integrals.
C         EPSL - desired precision of three-center integrals containing L operator
C         EPSR - desired precision of three-center integrals containing R operator
C
C      These values are re-initialized in PSDCNS.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON
     ./PSO3PR/ EPSL, EPSR
      SAVE /PSO3PR/
      DATA EPSL/2D-6/
      DATA EPSR/1D-6/
C
      END
C
