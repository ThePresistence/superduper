C     ******************************************************************
C
C     Compute overlap integrals for arbitrary quantum numbers.
C     Compute orbital and dipole moment integrals for L up to 3.
C
C     ******************************************************************
C
C     User-callable routines include:
C
C     PSSNRM  - Computes nurmalization coefficient for a Slater AO
C     PSOVX   - Computes overlap integrals
C     PSORX   - Computes two-center integrals of the (potentially)
C               arbitrary negative power of the distance.
C     PSSLM   - Computes values of the solid spherical harmonics of
C               potentially unlimited order.
C     PSOLX   - Computes orbital moment and related integrals
C     PSODS   - Computes arbitrary <\mu|O x|\nu> integrals
C     PSOQS   - Computes arbitrary integrals of type
C               <\mu|O \delta_{ab} r^2 - r_a r_b|\nu>.
C     PSODX   - Computes dipole moment integrals
C     PSOVP   - Computes two-center integrals with orbitals centered
C               at the same center.
C     PSOPR3  - Computes <mu|r^-3|nu> integrals with mu and nu on
C               one center and r^-3 on the other.
C     PSOQR3  - Computes magnetic dipole interaction integrals
C     PSOLB   - Computes <mu|O L |nu> integrals with L and nu centered
C               on different atoms.
C     PSLBR3  - Computes <mu|L/r^3|nu> integrals with mu and nu centered
C               on one atom and the operator on the other one.
C     PSLA3   - Computes <mu^A|L^A/r^3|nu^B> integrals.
C     PSRLA3  - Computes <mu^A r^A|L^A/r^3|nu^B> integrals.
C     PSRAB3  - Computes <mu^A r^A|(r^A)^-3|nu^B r^B> integrals.
C     PSR31C  - Computes one-center <\mu|r^-3|\nu> integrals
C     PSOV1C  - Computes one-center overlap integrals
C     PSSX94  - Converts integrals in the increasing-M order into
C               ad hoc order used in MNDO94.
C
      SUBROUTINE PSOVX(N1I,L1I,ALP1I,N2I,L2I,ALP2I,XI,YI,ZI,OVR,LDO)
C
C   Compute overlap integrals over block of Slater AOs.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N1I    - Major quantum number of the first AO
C      L1I    - Orbital quantum number of the first AO
C      ALP1I  - Orbital exponent of the first AO
C      N2I,L2I,ALP2I 
C             - Same for the second atom
C      XI,YI,ZI
C             - Relative coordinates of the second atom.
C      OVR    - Array for overlap integrals, at least 
C               (2*L1+1)*(2*L2+1) elements.
C      LDO    - Leading dimension of the OVR array.
C
C   Accessed common blocks:
C
C     PS3JDT  - Constants used by auxiliary functions.
C
C   Local storage:
C
C      ca. (MAXJ+1)**3 + 6 (MAXJ+1)**2 + 3 (MAXJ+1) DOUBLE PRECISION
C      cells. With MAXJ of 30, this amounts to 35650 cells, which is
C      large but probably not excessive yet.
C
C   Private functions:
C
C      All PSVX* routines correspond to final equations of Talman paper.
C      PS3J  computes needed 3J symbols
C      PSSLM computes (real) spherical harmonics.
C
C   Module logic:
C
C     J.D. Talman, Phys. Rev. A, 48, 243 (1993), with modifications
C     from W. Hierse and P.M. Oppeneer, Int. J. Quant. Chem. 52, 1249
C     (1994). Throughout the code, *1-variables correspond to non-primed
C     quantities of the paper, while *2-variables are primed ones.
C
C   Possible optimizations:
C
C     The execution time profile (on SGI Indigo^2) is relatively flat.
C     Two obvious candidates for the optimization are integer-integer
C     and double-integer exponentiation operations, which can be
C     eliminated in most places by careful loop unrolling. Another
C     good candidates for optimization are ps3j and psvxi. Caching
C     coefficients computed by psvxb, as suggested by Talman, is
C     probably worthless for low orbital moments we are interested in.
C
C     In the present state, the code is significantly slower than
C     Pople's sequence of implemented in MNDO94 for S-S overlaps,
C     approximately the same for P-P overlaps and somewhat faster
C     for D-D overlaps.
C
C   Bugs:
C
C     Sum N1+N2+L1+L2 should not exceed MAXJ, which means we
C     can do any type of overlap up to N=6, which is probably
C     adequate for now.
C
C     Expressions are singular at the zero interatomic distance.
C
C     The code is not strictly standard-conforming since the
C     non-standard DOUBLE COMPLEX type and associated DIMAG and
C     DCMPLX intrinsics are used.
C
C     If ALP1 and ALP2 are very dissimilar (i.e., with ALP1/ALP2
C     ratiou below 0.01 or above 100.), relative precision at
C     large interatomic separations is very low. Absolute precision
C     is still Ok provided that neither ALP1 nor ALP2 is zero.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXJ=30)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ONE=1.0D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      LOGICAL SWAP
      EXTERNAL PSSNRM
      DIMENSION OVR(-L1I:LDO-L1I-1,-L2I:L2I)
      DIMENSION BASF(0:MAXJ), BASG(0:MAXJ), DJ(0:MAXJ,0:MAXJ)
      DIMENSION Q((MAXJ+1)**2), F((MAXJ+1)**3), A(0:MAXJ)
      DIMENSION OVL(-MAXJ:MAXJ,-MAXJ:MAXJ)
C
      IF( N1I+N2I+L1I+L2I .GT. MAXJ ) THEN
          WRITE(NB6,11000) N1I, L1I, N2I, L2I
          STOP 'PSOVX'
      ENDIF
C
C    Our integral code can handle only L1.LE.L2 case. The second
C    possible case is handled by exchanging atoms and replacing
C    the computed integrals with complex conjugates.
C
      IF( L1I.LE.L2I ) THEN
          SWAP = .FALSE.
          L1   = L1I
          N1   = N1I
          ALP1 = ALP1I
          L2   = L2I
          N2   = N2I
          ALP2 = ALP2I
          X    = -XI
          Y    = -YI
          Z    = -ZI
      ELSE
          SWAP = .TRUE.
          L1   = L2I
          N1   = N2I
          ALP1 = ALP2I
          L2   = L1I
          N2   = N1I
          ALP2 = ALP1I
          X    = XI
          Y    = YI
          Z    = ZI
      ENDIF
C
      R = SQRT(X**2+Y**2+Z**2)
      X = X/R
      Y = Y/R
      Z = Z/R
C
C    Compute elemnetary basis functions and fold them into
C    J integrals and Q integrals
C
      CALL PSVXBF(HALF*(ALP1+ALP2)*R,N1+N2+L1+L2,BASF)
      CALL PSVXBG(HALF*(ALP2-ALP1)*R,N1+N2+L1+L2,BASG)
      CALL PSVXJ(N2+L1,N1-L1,L1+L2,R,BASF,BASG,DJ,MAXJ+1)
      CALL PSVXQ(L1,N2,L2,R,DJ,MAXJ+1,Q)
C
C    Compute F coefficients and fold Q integrals into A(/\)'s
C    Expressions for MINL and MINLB parameters of PSVXF are
C    simplified due to the inequality L2>=L1
C
      CALL PSVXF(L1,L2,L2-L1,L2-L1,F)
      CALL PSZRVB(2*L1+1,A(L2-L1))
      CALL DGEMV('T',(L1+1)*(2*L1+1),L1+1,ONE,F,2*(L1+1)*(2*L1+1),
     .               Q,1,ZERO,A(L2-L1),2)
C
C    Nomalize A(/\)
C
*     SCALE = SQRT((ALP1**(2*N1+1))*(ALP2**(2*N2+1))*RFACT(2*N1)*
*    .           RFACT(2*N2))*POW2(N1+N2+1)
      SCALE = PSSNRM(N1,ALP1)*PSSNRM(N2,ALP2)
      DO 1000 I=L2-L1,L2+L1,2
          A(I) = A(I) * SCALE
 1000 CONTINUE
C
C    Finally, convert local-system integrals A(/\) into 
C    overlap integrals, and normalize result. Expression
C    for the MINLB parameter of PSVXI is simplified due
C    to the inequality L2>=L1
C
      CALL PSVXI(L1,L2,L2-L1,X,Y,Z,A(L2-L1),OVL)
C
C    Store computed integrals, swapping centers in progress if necessary.
C
      IF( .NOT.SWAP ) THEN
          DO 1200 M2=-L2,L2
              DO 1100 M1=-L1,L1
                  OVR(M1,M2) = OVL(M1,M2)
 1100         CONTINUE
 1200     CONTINUE
      ELSE
          DO 1400 M1=-L1,L1
              DO 1300 M2=-L2,L2
                  OVR(M2,M1) = OVL(M1,M2)
 1300         CONTINUE
 1400     CONTINUE
      ENDIF
C
      RETURN
11000 FORMAT(' OVERLAP INTEGRAL REQUESTED FOR (',I3,',',I3,') (',
     .       I3,',',I3,') PAIR')
      END
C
      SUBROUTINE PSORX(N1,L1,ALP1,N2,L2,ALP2,XI,YI,ZI,IPOWN,OR,LDO)
C
C   Compute <\phi_a|r_b^-n|\phi_b> two-center integrals over Slater AOs.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N1I    - Major quantum number of the first AO
C      L1I    - Orbital quantum number of the first AO
C      ALP1I  - Orbital exponent of the first AO.
C      N2I    - Major quantum number of the second AO
C      L2I    - Orbital quantum number of the second AO
C      ALP2I  - Orbital exponent of the second AO. ALP2I
C               can be zero (Radial normalization by first
C               center AO will not be performed in this case).
C      XI,YI,ZI
C             - Relative coordinates of the second atom.
C      IPOWN  - Negative distance operator power, should be >= 1
C               (cases with IPOWN <= 0 are handled properly, but
C               using PSORX is wasteful in this case).
C               Operator is centered at the second atom.
C      OR     - Array for integrals, at least 
C               (2*L1+1)*(2*L2+1) elements.
C      LDO    - Leading dimension of the OVR array.
C
C   Accessed common blocks:
C
C     PS3JDT  - Constants used by auxiliary functions.
C
C   Local storage:
C
C     *Approximately* (MAXJ+1)**3 + 15/2 (MAXJ+1)**2 + 7 (MAXJ+1) +
C     3 (MAXSER+1), i.e. ca. 39,000 DOUBLE PRECISION cells with
C     MAXJ=30 and MAXSER=500.
C
C   Module logic:
C
C     J.D. Talman, Phys. Rev. A, 48, 243 (1993), with modifications
C     from W. Hierse and P.M. Oppeneer, Int. J. Quant. Chem. 52, 1249
C     (1994) and different expansion a-la R.M. Pitzer et al., J. Chem.
C     Phys. 37, 267 (1962) for Q_{\lambda L} integrals leading
C     to negative-M elemental integrals. Significant portions of the
C     code are carried over from PSOVX, which shares a good part of
C     auxiliary functions with PSORX.
C
C     Throughout the code, *1-variables correspond to non-primed
C     quantities of the paper, while *2-variables are primed ones.
C
C   Possible optimizations:
C
C   Bugs:
C
C     All problems of PSOVX apply here.
C
C     Expressions are singular at the zero interatomic distance.
C     Severe cancellation problems occure when ALP1 is large,
C     ALP2 is zero and distance is relatively large.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXJ=30)
      PARAMETER (MAXSER=500)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ONE=1.0D0)
      PARAMETER (TWO=2.0D0)
      PARAMETER (EPS=1.0D-16)
      PARAMETER (SMALL=1.D-20)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      EXTERNAL PSSNRM
      DIMENSION OR(-L1:LDO-L1-1,-L2:L2)
      DIMENSION QA((MAXJ+1)**2), QB((MAXJ+1)**2)
      DIMENSION QAG(0:MAXSER+MAXJ), QAS(0:MAXSER), QAI(0:MAXSER+MAXJ)
      DIMENSION QBF(-MAXJ:MAXJ), QBG(0:MAXJ), QBFT(1:MAXJ)
      DIMENSION QBK((MAXJ+1)**2), QBD(0:MAXJ/2,0:MAXJ)
      DIMENSION F((MAXJ+1)**3), A(0:MAXJ)
      DIMENSION OVL(-MAXJ:MAXJ,-MAXJ:MAXJ)
C    QA/QB, QBD, F, A
      IF( L1+L2.GT.MAXJ .OR.
C    QAG
     .    N2+L2+2*L1+1-IPOWN.GT.MAXSER .OR.
C    QBF
     .    N2-L2-IPOWN.LT.-MAXJ .OR. N1+N2-IPOWN.GT.MAXJ .OR.
C    QBG
     .    N1+L1+2*L2.GT.MAXJ .OR.
C    QBK
     .    (L2+L2+1)**2.GT.(MAXJ+1)**2 ) THEN
              WRITE(NB6,11000) N1,L1,N2,L2,IPOWN
              STOP 'PSORX'
      ENDIF
C
      R    = SQRT(XI**2+YI**2+ZI**2)
      X    = -XI/R
      Y    = -YI/R
      Z    = -ZI/R
      BET1 = ALP1*R
      BET2 = ALP2*R
      MINL = MAX(0,L2-L1)
      MAXL = L1+L2
C
C    First part of the integration: Compute Q^a by a series expansion
C
      CALL PSZRV((L1+1)*(MAXL-MINL+1),QA)
      QSCALE = HALF*(R**(N1+N2+1-IPOWN))
      EPSX   = EPS/QSCALE
      NTERMS = IPSRXL(EPSX,BET1,N1-L1)
      IF( NTERMS.EQ.10**6 ) THEN
C
C        Special case: overflow will occure
C
          WRITE(NB6,11400)
          NTERMS = MAXSER
      ELSE IF( NTERMS.GT.MAXSER ) THEN
C
C        Precision have to be lowered to fit into available space
C
          NTERM = NTERMS
          EPSY  = EPSX
  200     CONTINUE
              EPSY   = TWO*EPSY
              NTERMS = IPSRXL(EPSY,BET1,N1-L1)
              IF( NTERMS.GT.MAXSER ) GOTO 200
          WRITE(NB6,11500) EPSY/EPSX, NTERM
          EPSX  = EPSY
      ENDIF
C
C        Compute necessary values of the elementary integral and
C        series expansion of the integrand
C
      CALL PSXQAG(NTERMS+1+N2+L2+2*L1-IPOWN,BET2,QAG)
      CALL PSXQAS(NTERMS,N1-L1,BET1,QAS)
C
C        For each required L value, perform inner integration on it and, 
C        for each required value of LA, contract with the elementary integrals.
C        The ugly index expression in QA below is actually simulated
C        two-dimensional array.
C
      SCALE = EXP(-BET1)
      DO 1000 L=MINL,MAXL
          CALL PSXQAI(L,NTERMS,QAS,QAI)
          DO 900 LA=ABS(L-L2),MIN(L1,L+L2)
              IF( N2+LA-IPOWN+L+1.LT.0 ) THEN
                  WRITE(NB6,11100) LA, L
                  STOP 'PSORX'
              ENDIF
              QA(1+LA+(L-MINL)*(L1+1)) = SCALE*DDOT(NTERMS/2,
     .                   QAG(N2+LA-IPOWN+L+1),2,QAI(L+1),2)
  900     CONTINUE
 1000 CONTINUE
C
C        Compute elementary integrals for second part of the Q integral.
C        Elementary integrals include both usual A(N) and B(N) integrals
C        (which are lifted from PSOVX) and exponential integrals, which
C        are computed by a SLATEC function DEXINT. Result of DEXINT have
C        to be swapped, since they are generated in the opposite order
C        from that we want.
C
      IF( N1+N2-IPOWN.GE.0 ) THEN
          CALL PSVXBF(BET1+BET2,N1+N2-IPOWN,QBF(0))
      ENDIF
      IF( N2-L2-IPOWN.LT.0 ) THEN
          CALL DEXINT(BET1+BET2,1,1,IPOWN+L2-N2,MAX(D1MACH(4),EPSX),
     .                          QBFT,NZ,IERR)
          IF( IERR.NE.0 ) THEN
              WRITE(NB6,11200) IERR
              STOP 'PSORX'
          ENDIF
          IF( NZ.NE.0 ) THEN
              WRITE(NB6,11300) NZ
          ENDIF
          DO 2000 I=1,IPOWN+L2-N2
              QBF(-I) = QBFT(I)
 2000     CONTINUE
      ENDIF
      CALL PSVXBG(BET1,N1+L1+2*L2,QBG)
C
C        Compute intermediate K integrals and convert them to QB 
C        integrals. At the same time, accumulate total integral
C        value in the QA array.
C
      CALL PSXQBK(N2-L2-IPOWN,N2+L1-IPOWN,N1-L1,L1+L2,QBF,QBG,QBK)
      DO 3000 L=MINL,MAXL
          CALL PSVXD(L,QBD,MAXJ/2+1)
          DO 2900 LA=MAX(0,L-L2),L1
              QB(1+LA+(L-MINL)*(L1+1)) = PSXQBS(L,QBD,MAXJ/2+1,QBK,
     .                 N2-L2-IPOWN,N2+L1-IPOWN,L1+L2,N2+LA-L-IPOWN)
              QA(1+LA+(L-MINL)*(L1+1)) = QSCALE * 
     .             (QA(1+LA+(L-MINL)*(L1+1)) + 
     .              QB(1+LA+(L-MINL)*(L1+1)))
 2900     CONTINUE
 3000 CONTINUE
C
C    Compute F coefficients and fold Q integrals into A(/\)'s
C
      MINLB = ABS(L2-L1)
      MAXLB = L2+L1
      CALL PSVXF(L1,L2,MINL,MINLB,F)
      CALL PSZRVB(MAXLB-MINLB+2,A(MINLB))
      CALL DGEMV('T',(L1+1)*(MAXL-MINL+1),(MAXLB-MINLB)/2+1,ONE,F,
     .               2*(L1+1)*(MAXL-MINL+1),QA,1,ZERO,A(MINLB),2)
C
C    Nomalize A(/\) by prefactors of the orbitals with non-zero
C    exponents.
C
      SCALE = ONE
      IF( ALP1.GT.SMALL ) 
     .    SCALE = SCALE*PSSNRM(N1,ALP1)
*    .    SCALE = SCALE*SQRT(ALP1**(2*N1+1)*RFACT(2*N1)*POW2(2*N1+1))
      IF( ALP2.GT.SMALL ) 
     .    SCALE = SCALE*PSSNRM(N2,ALP2)
*    .    SCALE = SCALE*SQRT(ALP2**(2*N2+1)*RFACT(2*N2)*POW2(2*N2+1))
      DO 4000 I=MINLB,MAXLB,2
          A(I) = A(I) * SCALE
 4000 CONTINUE
C
C    Finally, convert local-system integrals A(/\) into 
C    overlap integrals
C
      CALL PSVXI(L1,L2,MINLB,X,Y,Z,A(MINLB),OVL)
C
C    Store computed integrals
C
      DO 5200 M2=-L2,L2
          DO 5100 M1=-L1,L1
              OR(M1,M2) = OVL(M1,M2)
 5100     CONTINUE
 5200 CONTINUE
C
      RETURN
11000 FORMAT(' STATIC ARRAYS WILL OVERFLOW IN PSORX FOR N1 = ',I4,
     .       ' L1 = ',I4,' N2 = ',I4,' L2 = ',I4,' IPOWN = ',I4)
11100 FORMAT(' SINGULAR INTEGRAL ENCOUNTERED IN PSORX, LA = ',I4,
     .       ' L = ',I4)
11200 FORMAT(' ERROR ',I4,' ENCOUNTERED IN DEXINT, PSORX WILL STOP')
11300 FORMAT(' UNDERFLOW DETECTED FOR ',I4,
     .       ' EXPONENTIAL INTEGRALS IN PSORX')
11400 FORMAT(' OVERFLOW IS LIKELY IN PSORX, INTEGRALS ARE UNRELIABLE.')
11500 FORMAT(' INTEGRALS PRECISION IN PSORX LOWERED BY A FACTOR OF ',
     .      G10.3/' FULL PRECISION WOULD HAVE REQUIRED MAXSER .GE. ',I6)
      END
C
      DOUBLE PRECISION FUNCTION PSSNRM(N,ALP)
C
C   Compute normalization coefficient for a Slater AO
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N      - Major quantum number of the Slater AO.
C      ALP    - Orbital exponent of the Slater AO.
C
C   Accessed common blocks:
C
C      PS3JDT - Various interesting constants.
C
C   Local storage:
C
C      Nope.
C
C   Module logic:
C
C   Possible optimizations:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXJ=30)
      PARAMETER (MXJ=90)
C
      COMMON /PS3JDT/FACT(0:MXJ+1),POW2(0:MXJ+1),RINT(0:MXJ+1),
     .               RFACT(0:MXJ+1), RSFACT(0:MXJ+1), GHALF(0:MXJ+1)
     ./PSPRT / NB6
      SAVE /PSPRT /
C
      IF( N.LT.0 .OR. N.GT.MAXJ/2 ) THEN
          WRITE(NB6,11010) N
          STOP 'PSSNRM'
      ENDIF
      PSSNRM = SQRT(ALP**(2*N+1)*RFACT(2*N)*POW2(2*N+1))
      RETURN
11010 FORMAT(' MAJOR QUANTUM NUMBER ',I5,' IS INVALID IN PSSNRM.')
      END
C
      INTEGER FUNCTION IPSRXL(EPS,BET1,IPOW)
C
C   Estimate minimum number of terms in series expansion necessary to
C   converge Q^a integral to a given absolute precision.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      EPS    - Required absolute precision
C      BET1   - Exponent value
C      IPOW   - R power at the displaced center
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      The algorithm employed is a horrible mixture of bracketing
C      (which is too slow by itself) and damped iterations (which
C      tend to converge to the *raising* edge of the series
C      expansion if left to its own devices).
C
C   Bugs:
C
C      The code below is horrible, but this can't be fixed without
C      resorting to (non-F77) DO WHILE, which I won't do.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HAVEDN
      PARAMETER (TWOPI=6.2831853071795864769252867665590057683943388D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ONE=1.0D0)
      PARAMETER (TWO=2.0D0)
      PARAMETER (TEN=10.0D0)
      PARAMETER (BIGEXP=650D0)
C
C    Adjust absolute precision by the expected scaling factor
C
      IF( BET1.GE.BIGEXP ) THEN
          IPSRXL = 10**6
          RETURN
      ENDIF
C
      EPSX = EPS*(HALF**IPOW)*EXP(BET1)
      IF( BET1.GT.ONE ) EPSX = EPSX/BET1**IPOW
C
C    Select an initial approximation for X by stepping up
C    until iterator is reasonable
C
      X   = TWO*BET1+TEN
   50 CONTINUE
          TEMP = LOG(EPSX*SQRT(TWOPI*X))/X
          IF( TEMP.LT.ONE ) GOTO 90 
          X = X*TWO
          GOTO 50
C
C    Solve transcedental equation by simple iterations, taking care
C    not to drop below BET1, which will give us wrong solution.
C    Additionally, we'll maintain top water mark so that our solution
C    won't shot into far away region for small BET1 values.
C    Since we actually want to get an integer, we can stop iterate
C    once we are reasonably close to one.
C
   90 CONTINUE
      HAVEDN = .FALSE.
  100 CONTINUE
          XNEW = BET1*EXP(ONE-TEMP)
          IF( ABS(XNEW-X).LT.HALF ) GOTO 200
          IF( XNEW.LT.X ) THEN
              IF( .NOT.HAVEDN ) THEN
                  HAVEDN = .TRUE.
                  XDOWN  = X
              ELSE IF( X.LT.XDOWN ) THEN
                  XDOWN = X
              ENDIF
          ELSE
              IF( HAVEDN ) XNEW = XDOWN
          ENDIF
          XOLD = X
          X    = HALF*(XOLD+XNEW)
  120 CONTINUE
          TEMP = LOG(EPSX*SQRT(TWOPI*X))/X
          IF( TEMP.LT.ONE ) GOTO 100
          X = HALF*(XOLD+X)
          GOTO 120
C
C    Round to a nearest integer and bail out.
C
  200 CONTINUE
      IPSRXL = NINT(XNEW+IPOW)
      RETURN
      END
C
      SUBROUTINE PSXQAG(MAXH,BET2,QAG)
C                                             1  h  -bet2 t
C   Compute values of the elementary integral | t  E        dt
C                                             0
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      MAXH   - Largest power of H required (smallest is 0)
C      BET2   - Exponent value
C      QAG    - Table of function values
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      The integral is actually a renormalized incomplete gamma
C      function: h! \gamma*(h+1,bet2). The integral value for
C      the largest h can therefore be computed as sum of the
C      series:
C                             n
C           -bet2      h! bet2
C          E      Sum ---------
C                 n=0  (h+1+n)!
C
C      Function values for smaller h can be computed by a stable
C      downwards recurrence:
C                       1                    -bet2
C          F   (bet2) = - ( bet2 F (bet2) + E     )
C           h-1         h         h
C
C      Using a library function for usually implemented \gamma
C      function is *not* acceptable, since unlike \gamma* it is
C      singular at bet2 = 0.
C
C      Both expansion and recurrence above are from Abramowitz
C      and Steguns's Handbook of Mathematical Functions.
C
C   Bugs:
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (EPS=1.D-16)
      PARAMETER (ONE=1.D0)
      DIMENSION QAG(0:MAXH)
C
      EXPB = EXP(-BET2)
C
C    Converge series for F(MAXH) first
C
      XN    = ONE/DBLE(MAXH+1)
      SUM   = XN
      ITERM = MAXH + 2
  100 CONTINUE
          XN    = XN * BET2/DBLE(ITERM)
          ITERM = ITERM + 1
          SUM   = SUM + XN
          IF( XN.GE.EPS ) GOTO 100
      FH = SUM*EXPB
      QAG(MAXH) = FH
C
C    Apply recurrence relation to get remaining integrals
C
      DO 200 IH=MAXH,1,-1
          FH = ( BET2*FH + EXPB ) / IH
          QAG(IH-1) = FH
  200 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSXQAS(NTERMS,IPOW,BET1,QAS)
C
C   Compute series expansion for (1+v)^ipow * exp(-beta v)
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NTERMS - Number of terms required in a series expansion
C      IPOW   - Power of (1+v) in the expression
C      BET1   - Exponent parameter in the expression
C      QAS    - (Output) coefficients of the series expansion
C
C   Accessed common blocks:
C
C      PS3JDT - Factorial values.
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXJ=30)
      PARAMETER (MXJ=90)
      PARAMETER (ONE =1.D0)
C
      COMMON /PS3JDT/FACT(0:MXJ+1),POW2(0:MXJ+1),RINT(0:MXJ+1),
     .               RFACT(0:MXJ+1), RSFACT(0:MXJ+1), GHALF(0:MXJ+1)
     ./PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION QAS(0:NTERMS)
      DIMENSION COMB(0:MAXJ)
C
      IF( IPOW.GT.MAXJ ) THEN
          WRITE(NB6,11000) IPOW
          STOP 'PSXQAS'
      ENDIF
C
C    Precompute table of number of permutations and initialize
C    coefficients to zero
C
      DO 100 I=0,IPOW
          COMB(I) = FACT(IPOW)*RFACT(I)*RFACT(IPOW-I)
  100 CONTINUE
      CALL PSZRV(NTERMS+1,QAS(0))
C
C    Go over the powers of BET1 updating all dependent coefficients
C
      TEMP = ONE
      DO 500 J=0,NTERMS
          DO 400 I=0,MIN(IPOW,NTERMS-J)
              QAS(I+J) = QAS(I+J) + COMB(I) * TEMP
  400     CONTINUE
          TEMP = TEMP * (-BET1)/(J+1)
  500 CONTINUE
C
      RETURN
11000 FORMAT(' STATIC TABLE OVERFLOW IN PSXQAS')
      END
C
      SUBROUTINE PSXQAI(L,NTERMS,QAS,QAI)
C
C   Perform inner integration for the QA.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L      - Degree of the Legenre polynome convolved with the
C               series.
C      NTERMS - Number of terms in series expansion
C      QAS    - (Input) coefficients of the series expansion
C      QAI    - (Output) coefficients of the series expansion
C
C   Accessed common blocks:
C
C      PS3JDT - Factorial, half-integer gamma function values and
C               powers of two are all here.
C
C   Local storage:
C
C   Module logic:
C                                                     2   2
C                                             t  N   t - v -2 v
C       The subprogram computes the integral: | v P (----------) dv,
C                                            -t    L    2 t
C       where P is a Legendre polynomial of the L-th order, for all
C              L
C       terms in the input series expansion. The expressions for 
C       coefficients are somethat esoteric, but not *that* difficilt ;-)
C       First L+1 terms of the output series are always zero due
C       to the nature of Legendre polynomial.
C
C   Bugs:
C
C       For short expansions, some of the coefficients may be computed
C       unnecessarily.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXJ=30)
      PARAMETER (MXJ=90)
      PARAMETER (HALF =0.5D0)
      PARAMETER (ONEP5=1.5D0)
C
      COMMON /PS3JDT/FACT(0:MXJ+1),POW2(0:MXJ+1),RINT(0:MXJ+1),
     .               RFACT(0:MXJ+1), RSFACT(0:MXJ+1), GHALF(0:MXJ+1)
     ./PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION QAS(0:NTERMS), QAI(0:NTERMS+L+1)
      DIMENSION CEVEN(0:MAXJ/2+1), CODD(0:MAXJ/2+1)
C
      IF( L.GE.MAXJ ) THEN
          WRITE(NB6,11000)
          STOP 'PSXQAI'
      ENDIF
C
C    Compute initial values of coefficients, separately for even and
C    odd values of N
C
      DO 100 IS=0,L/2
          CEVEN(IS) = POW2(2*IS  )*FACT(IS)*FACT(L  -IS)*RFACT(2*IS  )
     .               *RFACT(L  -2*IS)*GHALF(IS  )*GHALF(L-IS)/GHALF(L+1)
  100 CONTINUE
      DO 200 IS=0,(L+1)/2-1
          CODD (IS) =-POW2(2*IS+1)*FACT(IS)*FACT(L-1-IS)*RFACT(2*IS+1)
     .               *RFACT(L-1-2*IS)*GHALF(IS+1)*GHALF(L-IS)/GHALF(L+1)
  200 CONTINUE
      CALL PSZRV(NTERMS+L+2,QAI(0))
C
C    Process even and odd terms simultaneously. In each case, update 
C    coefficients for use on the next iteration.
C
      DO 1000 N=0,NTERMS-1,2
          CE         = QAS(N)
          CO         = QAS(N+1)
          DO 600 IS=0,MIN(N+2,L+1)/2-1
              QAI(N+L+1-2*IS) = QAI(N+L+1-2*IS)+CE*CEVEN(IS)+CO*CODD(IS)
              TEMP      = DBLE(N/2+1)/
     .                    (DBLE(N/2+1-IS)*(DBLE(N/2+L-IS)+ONEP5))
              CEVEN(IS) = CEVEN(IS)*TEMP*(DBLE(N/2)+HALF)
              CODD (IS) = CODD (IS)*TEMP*(DBLE(N/2)+ONEP5)
  600     CONTINUE
          DO 800 IS=MIN(N+2,L+1)/2,MIN(N,L)/2
              QAI(N+L+1-2*IS) = QAI(N+L+1-2*IS)+CE*CEVEN(IS)
              TEMP      = DBLE(N/2+1)/
     .                    (DBLE(N/2+1-IS)*(DBLE(N/2+L-IS)+ONEP5))
              CEVEN(IS) = CEVEN(IS)*TEMP*(DBLE(N/2)+HALF)
  800     CONTINUE
 1000 CONTINUE
C
C    If necessary, process the last (even) coefficient in the input
C    series. 
C
      IF( MOD(NTERMS,2).EQ.0 ) THEN
          CE         = QAS(NTERMS)
          DO 1100 IS=0,MIN(NTERMS,L)/2
              QAI(NTERMS+L+1-2*IS) = QAI(NTERMS+L+1-2*IS)+CE*CEVEN(IS)
 1100     CONTINUE
      ENDIF
C
      RETURN
11000 FORMAT(' STATIC TABLE OVERFLOW IN PSXQAI')
      END
C
      SUBROUTINE PSXQBK(MINM,MAXM,N,MAXP,QBF,QBG,QBK)
C
C   Compute integrmediate integral K for QB.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      MINM   - Minimum value of the first parameter
C      MAXM   - Maximum value of the first parameter
C      N      - Value of the second parameter
C      MAXP   - Maximum value of the third parameter
C               (minimum is always zero).
C      QBF    - F basis functions (\alpha_n or E_i exponential
C               integrals, depending on the index value)
C      QBG    - G basis functions (\beta_n exponential
C               integrals)
C      QBK    - (Output) K integrals
C
C   Accessed common blocks:
C
C      PS3JDT - Factorial, inverse factorial and powers of two.
C
C   Local storage:
C
C   Module logic:
C
C                           Inf  M  -bet2 t    t+1   2  2 p N  -bet1 u
C      K(M,N,p) is actually  |  t  E        dt  |  (t -u ) u  E        du,
C                            1                 t-1
C      and is computed by a binomial expansion to elementary integrals.
C      Verifying the validity of the expression is left as an
C      exersise to an interested reader ;-)
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXJ=30)
      PARAMETER (MXJ=90)
      PARAMETER (ZERO=0.D0)
C
      COMMON /PS3JDT/FACT(0:MXJ+1),POW2(0:MXJ+1),RINT(0:MXJ+1),
     .               RFACT(0:MXJ+1), RSFACT(0:MXJ+1), GHALF(0:MXJ+1)
     ./PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION QBF(-MAXJ:MAXJ), QBG(0:MAXJ), QBK(MINM:MAXM,0:MAXP)
C
      DIMENSION A(0:MAXJ), B(0:MAXJ)
C
      IF( MAXP.GT.MAXJ .OR. N.GT.MAXJ ) THEN
          WRITE(NB6,11000)
          STOP 'PSXQBK'
      ENDIF
C
C    Precompute number of pernutations for all possible counts out of N,
C    since these do not actually change throughout the computation.
C
      DO 100 I=0,N
          A(I) = FACT(N)*RFACT(I)*RFACT(N-I)
  100 CONTINUE
C
      DO 2000 IP=0,MAXP
C
C        The expression below actually differs from the correct
C        one by a factor of two. However, since our QBG integrals
C        are two times smaller than they should have been, these
C        two errors cancel out.
C
          DO 1100 I=0,IP
              B(I) = FACT(IP)*RFACT(I)*RFACT(IP-I)*POW2(IP-I+1)
 1100     CONTINUE
          DO 1800 M=MINM,MAXM-IP
              SUM = ZERO
              DO 1500 I=0,IP
                  DO 1400 J=0,N
                      SUM = SUM + B(I)*A(J)*QBF(M+IP-I+N-J)*QBG(IP+I+J)
 1400             CONTINUE
 1500         CONTINUE
              IF( MOD(IP,2).EQ.1 ) SUM = -SUM
              QBK(M,IP) = SUM
 1800     CONTINUE
 2000 CONTINUE
      RETURN
11000 FORMAT(' STATIC TEMPORARY OVERFLOWED IN PSXQBK ')
      END
C
      DOUBLE PRECISION FUNCTION PSXQBS(L,QBD,LDD,QBK,MINM,MAXM,MAXP,
     .                                 IBASM)
C
C   Compute intermediate integral Q^b.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L      - Total orbital moment
C      QBD    - Transformation coefficients
C      LDD    - Leading dimension of the QBD matrix
C      QBK    - K integrals
C      MINM   - Minimum value of the first parameter of K integrals
C      MAXM   - Maximum value of the first parameter of K integrals
C      MAXP   - Maximum value of the third parameter of K integrals
C               (minimum is always zero).
C      IBASM  - Base shift of first parameter of K integrals
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
      PARAMETER (ZERO=0.D0)
C
      DIMENSION QBD(0:LDD-1,0:*), QBK(MINM:MAXM,0:MAXP)
C
      PSXQBS = ZERO
      DO 200 K=0,L/2
          DO 100 IP=0,L-2*K
              PSXQBS = PSXQBS + QBD(K,IP)*QBK(IBASM+2*K,IP)
  100     CONTINUE
  200 CONTINUE
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PS3J(J1,J2,J3,M1,M2,M3)
C
C                               / J1  J2  J3 \
C   Computes Wigner 3j simbols  |            |
C                               \ M1  M2  M3 /
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      J*,M*  - Integers.
C
C   Accessed common blocks:
C
C     PS3JDT  - Various useful constants.
C
C   Local storage:
C
C   Module logic:
C
C     3J simbols are computed according to Racah's expression
C     taken from Landau and Lifshits.
C
C   Precision:
C
C     3J symbols were tested up to J=40. Absolute precision is
C     always about 2e-16 (i.e., machine precision). Relative
C     accuracy may be down to 5e-13 for small coefficients.
C
C   Bugs:
C
C     Half-integral orbital moments are not supported.
C     Sum of the orbital moments should be less than MAXJ.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXJ=90)
      PARAMETER (ZERO = 0.D0)
      PARAMETER (ONE  = 1.D0)
      PARAMETER (TEN  =10.D0)
C    Note that IHUGE *should* be exactly divisible by 2, since we
C    might want to take square root of it.
      PARAMETER (IHUGEH=25)
      PARAMETER (IHUGE =2*IHUGEH)
      PARAMETER (ISMALL=-10*IHUGE)
      PARAMETER (HUGE  =10.D0**IHUGE)
      PARAMETER (SHUGE =10.D0**IHUGEH)
      PARAMETER (TINY  =ONE/HUGE)
      PARAMETER (STINY =ONE/SHUGE)
C
      COMMON /PS3JDT/FACT(0:MXJ+1),POW2(0:MXJ+1),RINT(0:MXJ+1),
     .               RFACT(0:MXJ+1), RSFACT(0:MXJ+1), GHALF(0:MXJ+1)
     ./PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION TEMP(10)
C
      IJ0 = J1+J2+J3+1
      IJ1  =   J1 + J2 - J3
      IJ2  =   J1 - J2 + J3
      IJ3  = - J1 + J2 + J3
      IM1A =   J1 - M1
      IM1B =   J1 + M1
      IM2A =   J2 - M2
      IM2B =   J2 + M2
      IM3A =   J3 - M3
      IM3B =   J3 + M3
      IF( IJ1.LT.0 .OR. IJ2.LT.0 .OR. IJ3.LT.0 .OR. 
     .         IM1A.LT.0 .OR. IM1B.LT.0 .OR. 
     .         IM2A.LT.0 .OR. IM2B.LT.0 .OR. 
     .         IM3A.LT.0 .OR. IM3B.LT.0 .OR. M1+M2+M3.NE.0 ) THEN
          PS3J = ZERO
          RETURN
      ENDIF
C
C    Now, coefficient is not zero, so that we'll have to compute it.
C
      IF( IJ0.GT.MXJ ) THEN
          WRITE(NB6,1000) IJ0, MXJ
          STOP 'PS3J'
      ENDIF
C
      IZT1  = J2-J3-M1
      IZT2  = J1+M2-J3
      IZT3  = J1-M1
      IZT4  = J2+M2
      IMINZ = MAX(0,IZT1,IZT2)
      IMAXZ = MIN(IJ1,IZT3,IZT4)
C    Multiply all the factors taking care of under/overflows
      TEMP( 1) = FACT(IJ1)
      TEMP( 2) = FACT(IJ2)
      TEMP( 3) = FACT(IJ3)
      TEMP( 4) = FACT(IM1A)
      TEMP( 5) = FACT(IM1B)
      TEMP( 6) = FACT(IM2A)
      TEMP( 7) = FACT(IM2B)
      TEMP( 8) = FACT(IM3A)
      TEMP( 9) = FACT(IM3B)
      TEMP(10) = RFACT(IJ0)
      SCALE = TEMP(1)
      IBALS = 0
      DO 50 I=2,10
          SCALE = SCALE * TEMP(I)
              IF( SCALE.GT.HUGE ) THEN
              SCALE = SCALE * TINY
              IBALS = IBALS + IHUGE
          ELSEIF( SCALE.LT.TINY ) THEN
              SCALE = SCALE * HUGE
              IBALS = IBALS - IHUGE
          ENDIF
   50 CONTINUE
      SCALE = SQRT(SCALE)
      IBALS = IBALS/2
      IF( MOD(IBALS,IHUGE).NE.0 ) THEN
          IF( IBALS.GT.0 ) THEN
              SCALE = SCALE*SHUGE
              IBALS = IBALS-IHUGE/2
          ELSE
              SCALE = SCALE*STINY
              IBALS = IBALS+IHUGE/2
          ENDIF
      ENDIF
C
      PS3J  = ZERO
      IPS3J = ISMALL
      DO 100 IZ=IMINZ,IMAXZ,1
          TEMP(1) = RFACT(IZ)
          TEMP(2) = RFACT(IJ1-IZ)
          TEMP(3) = RFACT(IZT3-IZ)
          TEMP(4) = RFACT(IZT4-IZ)
          TEMP(5) = RFACT(IZ-IZT1)
          TEMP(6) = RFACT(IZ-IZT2)
          PROD  = SCALE
          IBALP = IBALS
          DO 80 I=1,6
              PROD = PROD*TEMP(I)
              IF( PROD.LT.TINY ) THEN
                  PROD  = PROD * HUGE
                  IBALP = IBALP - IHUGE
              ENDIF
   80     CONTINUE
          IF( IBALP.GE.IPS3J) THEN
              IF( MOD(IZ,2).NE.0 ) PROD = -PROD
              IF( IBALP.GT.IPS3J ) THEN
                  PS3J  = PROD
                  IPS3J = IBALP
              ELSE
                  PS3J  = PS3J + PROD
              ENDIF
          ENDIF
  100 CONTINUE
      IF( MOD(J1-J2-M3,2).NE.0 ) PS3J = -PS3J
      PS3J = PS3J*(TEN**IPS3J)
C
      RETURN
 1000 FORMAT(' SUM OF THE ORBITAL MOMENTS IN PS3J (', I5, 
     .       ') EXCEEDS COMPILED-IN LIMIT OF ', I3)
      END
C
      SUBROUTINE PSSLM(L,X,Y,Z,YLM)
C
C   Computes values of (real) spherical harmonics.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L      - Orbital moment
C      X,Y,Z  - *Normalized* coordinates of the point to compute
C               spherical harmonics at.
C      Y      - Output array of spherical harmonics, 2*L+1 values 
C               are stored.
C
C   Accessed common blocks:
C
C     PS3JDT  - Various useful constants.
C
C   Local storage:
C
C     (3/2) MAXJ DOUBLE PRECISION cells.
C
C   Module logic:
C
C     Values of spherical harmonics are computed according to 
C     expressions given in H.B.Schlegel, M.J.Frisch, Int.J.Quant.Chem.,
C     54, 83 (1995). Spherical harmonics defined in this way differ by
C     a factor of (-1)**M from the expression in Mathematica or 
C     in Landau&Lifshits and by (-i)**M from Talman's paper.
C
C   Bugs:
C
C     If X**2+Y**2+Z**2 is not unity, results are undefined.
C     Spherical harmonics with L>MAXJ/2 are not supported.
C     If spherical harmonics for several L are desired, repeated
C     calls to PSYLM are (relatively) inefficient.
C
C     This subroutine is not standard-conforming, since the
C     DOUBLE COMPLEX type is not included in the standard.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXJ=30)
      PARAMETER (MXJ=90)
      PARAMETER (HALF=0.5D0)
      PARAMETER (TWOPI=6.2831853071795864769252867665590057683943388D0)
      PARAMETER (SQ2=1.414213562373095048801688724209698078569671875D0)
C
      COMMON /PS3JDT/FACT(0:MXJ+1),POW2(0:MXJ+1),RINT(0:MXJ+1),
     .               RFACT(0:MXJ+1), RSFACT(0:MXJ+1), GHALF(0:MXJ+1)
C
      DOUBLE COMPLEX DCMPLX, CTEMP
      DOUBLE PRECISION A(0:MAXJ/2),B(0:MAXJ), YLM(-L:L)
      INTRINSIC DCMPLX, DIMAG
C
*     IF(L.GE.MAXJ/2) THEN
*         WRITE(NB6,1000) L, MAXJ/2-1
*         STOP 'PSYLK'
*     ENDIF
C
      DO 100 I=0,L/2
          A(I) = ((-1)**I)*FACT(2*(L-I))*RFACT(I)*RFACT(L-I)
  100 CONTINUE
      DO 200 I=0,L-1
          B(I) = Z**(L-I)*RFACT(L-I)
  200 CONTINUE
      B(L) = RFACT(0)
      SCALE=(L+HALF)/(TWOPI*POW2(2*L))
C
      TEMP = A(0)*B(0)
      DO 300 K=1,L/2
          TEMP = TEMP + A(K)*B(2*K)
  300 CONTINUE
      YLM(0) = TEMP * SQRT(SCALE)
C
      DO 500 M=1,L
          TEMP = A(0)*B(M)
          DO 400 K=1,(L-M)/2
              TEMP = TEMP + A(K)*B(M+2*K)
  400     CONTINUE
          TEMP = TEMP * SQ2 * SQRT(SCALE*FACT(L-M)*RFACT(L+M))
          CTEMP   = TEMP * (DCMPLX(X, Y)**M)
          YLM( M) = DBLE(CTEMP)
          YLM(-M) = DIMAG(CTEMP)
  500 CONTINUE
      RETURN
*1000 FORMAT(' PSYLK CALLED FOR SPHERICAL HARMONICS WITH L = ', I5, 
*    .       ' WHILE COMPILE-TIME LIMIT WAS SET TO ', I3 )
      END
C
      DOUBLE PRECISION FUNCTION PSVXFI(L1,L2,L3)
C
C   Computes a special case of 3J symbol necessary in F coefficient
C   evaluation.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L1,L2,L3 - Integers
C
C   Accessed common blocks:
C
C     PS3JDT  - Various useful constants.
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXJ=90)
      PARAMETER (ZERO=0.D0)
C
      COMMON /PS3JDT/FACT(0:MXJ+1),POW2(0:MXJ+1),RINT(0:MXJ+1),
     .               RFACT(0:MXJ+1), RSFACT(0:MXJ+1), GHALF(0:MXJ+1)
C
      I0 =   L1 + L2 + L3
      I1 =   L1 + L2 - L3
      I2 = - L1 + L2 + L3
      I3 =   L1 - L2 + L3
C
      IF( I1.LT.0 .OR. I2.LT.0 .OR. I3.LT.0 .OR. MOD(I0,2).EQ.1 ) THEN
          PSVXFI = ZERO
          RETURN
      ENDIF
C
      PSVXFI = FACT(I0/2)*FACT(I1)*
     .         RFACT(I1/2)*RFACT(I2/2)*RFACT(I3/2)*RFACT(I0+1)
      RETURN
      END
C
      SUBROUTINE PSVXF(L1,L2,MINL,MINLB,F)
C
C   Computes F coefficients required for overlap integrals computation.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L1     - Orbital moment of the first set of orbitals
C      L2     - Orbital moment of the second set of orbitals.
C               L2 should be >= L1
C      MINL   - Should be equal to MAX(L2-L1,0)
C      MINLB  - Should be equal to MAX(L2-L1,MOD(L2+L1,2))
C               Nope, to ABS(L2-L1)
C      F      - Space for the computed F coefficients, at least
C               ((2*L1+1)**2)*(L+1) values.
C
C   Accessed common blocks:
C
C     PS3JDT  - Various useful constants.
C
C   Local storage:
C
C     (MAXJ+1)(2*MAXJ+1) DOUBLE PRECISION cells.
C
C   Module logic:
C
C     Non-zero F coefficients are computed in stored in such an order
C     so that a single call to DGEMV is sufficient to produce all
C     A(/\)'s from Q(lamba,L). Only /\'s with the same parity as
C     the L1+L2 sum are computed (since the complementary set of
C     /\'s is indeterminate ;-)
C
C     Results computed by this code agree with the values computed
C     according to eq. (2.28) of Talman's paper rather than with
C     eq. (2.33), which is an inequivalent rendition of (2.28)
C
C   Bugs:
C
C     If L2 < L1, all the hell will broke lose.
C     Half of the space required for F is wasted.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXJ=30)
      PARAMETER (FOURPI=12.566370614359172953850573533118011536788678D0)
C
C                 lambda      L          /\
C
      DIMENSION F(0:L1,MINL:L2+L1,MINLB:L2+L1)
      DIMENSION A(0:MAXJ,0:2*MAXJ)
C
      SCALE = SQRT(FOURPI*(2*L1+1)*(2*L2+1))
      DO 200 L=MINL,L2+L1
          DO 100 LA=0,L1
              A(LA,L) = SCALE*((-1)**(L1+LA))*(2*L+1)*PSVXFI(L2,L,LA)
  100     CONTINUE
  200 CONTINUE
C
      DO 500 LB=MINLB,L2+L1,2
          SCALE = SQRT(DBLE(2*LB+1))/PSVXFI(LB,L2,L1)
          DO 400 L=MINL,L2+L1
              DO 300 LA=0,L1
                  F(LA,L,LB) = SCALE*A(LA,L)*PSVXFI(L,LB,L1-LA)
  300         CONTINUE
  400     CONTINUE
  500 CONTINUE
      RETURN
      END
C
      SUBROUTINE PSVXBF(U,MAXS,F)
C
C   Computes F basis functions required for overlap integral computation.
C   (These are actually \alpha_n exponential integrals, but let's keep quiet
C    about it ;-)
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      U      - Value of parameter
C      MAXS   - Maximum order of the function
C      F      - Computed basis function values, at least MAXS+1 elements
C               should be provided.
C
C   Accessed common blocks:
C
C     PS3JDT  - Various useful constants.
C
C   Local storage:
C
C   Module logic:
C
C      Eq. (2.19) of Talman's paper. Order of computation is slightly
C      rearranged to minimize the risk of under- or over- flow and
C      avoid divisions wherever possible.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXJ=90)
      PARAMETER (ONE=1.0D0)
C
      COMMON /PS3JDT/FACT(0:MXJ+1),POW2(0:MXJ+1),RINT(0:MXJ+1),
     .               RFACT(0:MXJ+1), RSFACT(0:MXJ+1), GHALF(0:MXJ+1)
C
      DIMENSION F(0:MAXS)
C
*     IF( MAXS.GT.MAXJ ) THEN
*         WRITE(NB6,1000) MAXS, MAXJ
*         STOP 'PSVXBF'
*     ENDIF
C
      RU     = ONE/U
      USCALE = EXP(-U)*RU
      SUM    = ONE
      X      = ONE
      F(0)   = USCALE
C
      DO 100 I=1,MAXS
          X      = X * U * RINT(I)
          SUM    = SUM + X
          USCALE = USCALE * I * RU
          F(I)   = USCALE * SUM
  100 CONTINUE
C
      RETURN
*1000 FORMAT(' PSVXBF CALLED WITH MAXS = ', I5, 
*    .       ' WHILE COMPILE-TIME LIMIT IS ', I3 )
      END
C
      SUBROUTINE PSVXBG(V,MAXS,G)
C
C   Computes G basis functions required for overlap integral computation.
C   (These are actually \beta_n exponetial integrals)
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      V      - Value of parameter
C      MAXS   - Maximum order of the function
C      G      - Computed basis function values, at least MAXS+1 elements
C               should be provided.
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Eqs. (60)-(64) of W.Hierse and P.M. Oppeneer, Int.J. Quant.Chem, 52,
C      1249 (1994) with a suitable correction by (-1)**S to account for the
C      difference with Tallman's expression.
C
C      Half of the value of the integral:
C
C                  1
C                  /  s  -v t
C      beta (v) =  | t  E     d t  (See eq. 5.1.6 in Abramowitz & Stegun)
C          s       /
C                 -1
C
C      is computed using stable recurrences.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (ONE=1.0D0)
      PARAMETER (TWO=2.0D0)
      PARAMETER (THREE=3.0D0)
      PARAMETER (FOUR=4.0D0)
      PARAMETER (EPS=1.D-16)
      PARAMETER (TINYV=1.D-30)
      PARAMETER (HUGEV=300.D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION G(0:MAXS)
C    If V is very small, overflow may occur in general expressions below.
C    In this case, we use exact limits for v -> 0
      IF( ABS(V).LT.TINYV ) THEN
          DO 10 I=1,MAXS,2
              G(I) = ZERO
   10     CONTINUE
          DO 15 I=0,MAXS,2
              G(I) = ONE/DBLE(I+1)
   15     CONTINUE
          RETURN
      ENDIF
C    If V is very large, the results in not representable. In this case, we'll
C    have to die...
      IF( ABS(V).GT.HUGEV ) THEN
          WRITE(NB6,11010) MAXS, V
          STOP 'PSVXBG'
      ENDIF
C
      SINHV = SINH(V)
      COSHV = COSH(V)
      IF( ABS(V).GT.MAXS ) THEN
C
C         Upward recurrence is stable
C
          RV   = ONE/V
          G(0) = SINHV*RV
          DO 100 I=0,MAXS-2,2
              G(I+1) = RV*( -(I+1)*G(I+0) + COSHV )
              G(I+2) = RV*( -(I+2)*G(I+1) + SINHV )
  100     CONTINUE
          IF( MOD(MAXS,2).EQ.1 ) G(MAXS) = RV*(-MAXS*G(MAXS-1) + COSHV)
      ELSE
C
C         Downward recurrence is stable, compute G(MAXS) first and start
C         from here. We can't use RINT table at this point, since there
C         is no guarantee what the sum will converge in MAXJ-MAXS iterations.
C         Note that although V can be negative, the termination test below
C         is Ok, since TEMP is multiplied by V twice.
C
          DI   = MAXS+1
          TEMP = ONE/DI
          SUMA = TEMP
          SUMB = ZERO
          V2   = V*V
          V3   = V*V*V
C        The following loop is written to provide good performance when
C        divisions and branches are expensive.
  200     CONTINUE
              DI1  = DI + ONE
              DI2  = DI + TWO
              DI3  = DI + THREE
              DI4  = DI + FOUR
              TEMP = TEMP*V/( DI1*DI2*DI3*DI4 )
              
              SUMB = SUMB + TEMP*DI4*( DI2*DI3 + V2 )
              SUMA = SUMA + TEMP*V  *( DI3*DI4 + V2 )

              DI   = DI4
              TEMP = TEMP*V3
          IF( TEMP.GE.EPS ) GOTO 200
          IF( MOD(MAXS,2).EQ.1 ) THEN
              G(MAXS)   = SINHV*SUMA - COSHV*SUMB
              G(MAXS-1) = (-V*G(MAXS)+COSHV)/MAXS
          ELSE
              G(MAXS)   = COSHV*SUMA - SINHV*SUMB
          ENDIF
C
          DO 500 I=2*(MAXS/2),2,-2
              G(I-1) = (-V*G(I-0)+SINHV)/DBLE(I-0)
              G(I-2) = (-V*G(I-1)+COSHV)/DBLE(I-1)
  500     CONTINUE
      ENDIF
C
      DO 600 I=1,MAXS,2
          G(I) = -G(I)
  600 CONTINUE
C
      RETURN
11010 FORMAT(' NON-REPRESENTABLE VALUES MAY OCCURE IN ',
     .       'PSVXBG FOR MAXS = ',I4,' V = ',G15.7)
      END
C
      SUBROUTINE PSVXB(M,N,B)
C
C   Computes B coefficients required for overlap integrals computation.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      M,N    - Row of B (third rank tensor) to compute.
C      B      - Space for the computed coefficients, at least 
C               M+N+1 elements
C
C   Accessed common blocks:
C
C     PS3JDT  - Various useful constants.
C
C   Local storage:
C
C   Module logic:
C
C     Expr. (2.13) of Talman's paper.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXJ=90)
      PARAMETER (ZERO=0.0D0)
C
      COMMON /PS3JDT/FACT(0:MXJ+1),POW2(0:MXJ+1),RINT(0:MXJ+1),
     .               RFACT(0:MXJ+1), RSFACT(0:MXJ+1), GHALF(0:MXJ+1)
     ./PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION B(0:M+N)
C
*     IF( M+N.GT.MAXJ ) THEN
*         WRITE(NB6,1000) M+N, MAXJ
*         STOP 'PSVXB'
*     ENDIF
C
      SCALE = FACT(M)*FACT(N)*((-1)**N)/POW2(M+N+1)
      DO 200 MU=0,M+N
          SUM = ZERO
          DO 100 I=MAX(0,MU-N),MIN(M,MU)
              SUM = SUM + ((-1)**I)*RFACT(M-I)*RFACT(I)
     .                             *RFACT(N-MU+I)*RFACT(MU-I)
  100     CONTINUE
          B(MU) = SCALE * ((-1)**MU) * SUM
  200 CONTINUE
C
      RETURN
*1000 FORMAT(' PSVXB CALLED WITH M+N = ', I5, 
*    .       ' WHILE COMPILE-TIME LIMIT WAS SET TO ', I3 )
      END
C
      SUBROUTINE PSVXJ(MAXM,N,MAXP,R,F,G,DJ,LDJ)
C
C   Computes elementary J integrals required for overlap integral
C   evaluation.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      MAXM   - Maximum value of M parameter
C      N      - Value of N parameter
C      MAXP   - Maximum value of P parameter
C      R      - Intercenter distance
C      F      - Values of F basis functions, at least MAXM+N+MAXP+1
C               elements
C      G      - Values of G basis functions, at least MAXM+N+MAXP+1
C               elements
C      DJ     - Space for integrals, at least (MAXM+1)*(MAXP+1)
C               elements
C      LDJ    - Leading dimension of the J matrix.
C
C   Accessed common blocks:
C
C   Local storage:
C
C      MAXJ DOUBLE PRECISION cells for a row of B tensor.
C
C   Module logic:
C
C      Expr. (2.16) of Talman's paper.
C
C   Bugs:
C
C      Some of the integrals are computed unnecessarily sometimes.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXJ=30)
      PARAMETER (ZERO=0.0D0)
*     COMMON /PSPRT / NB6
*     SAVE /PSPRT /
C
      DIMENSION F(0:MAXM+N+MAXP), G(0:MAXM+N+MAXP), DJ(0:LDJ-1,0:MAXP)
      DIMENSION BMN(0:MAXJ)
C
*     IF( LDJ.LT.MAXM+1 ) THEN
*         WRITE(NB6,1000) LDJ, MAXM+1
*         STOP 'PSVXJ'
*     ENDIF
C
      DO 300 M=0,MAXM
          CALL PSVXB(M,N,BMN)
          DO 200 IP=0,MAXP
              SUM = ZERO
              DO 100 MU=0,M+N
                  SUM=SUM+BMN(MU)*F(MU+IP)*G(M+N+IP-MU)
  100         CONTINUE
              DJ(M,IP) = SUM*(R**(M+N+2*IP+2))
  200     CONTINUE
  300 CONTINUE
C
      RETURN
*1000 FORMAT(' PSVXJ CALLED WITH LDJ = ', I3, ' WHILE AT LEAST ', I3, 
*    .       ' IS REQUIRED ')
      END
C
      SUBROUTINE PSVXD(L,D,LDD)
C
C   Computes D coefficients required for overlap integrals computation.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L      - Set of coefficients to compute
C      D      - Space for coefficients, at least (L+1)*(L/2+1) elements.
C      LDD    - Leading dimension of D matrix
C
C   Accessed common blocks:
C
C     PS3JDT  - Various useful constants.
C
C   Local storage:
C
C   Module logic:
C
C     Eq. (2.9) of Talman's paper.
C
C   Bugs:
C
C     Half of the space allocated for D is wasted.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXJ=90)
C
      COMMON /PS3JDT/FACT(0:MXJ+1),POW2(0:MXJ+1),RINT(0:MXJ+1),
     .               RFACT(0:MXJ+1), RSFACT(0:MXJ+1), GHALF(0:MXJ+1)
*    ./PSPRT / NB6
*     SAVE /PSPRT /
C
      DIMENSION D(0:LDD-1,0:L)
C
*     IF( LDD.LE.L/2 ) THEN
*         WRITE(NB6,1000) LDD, L+1
*         STOP 'PSVXD'
*     ENDIF
*     IF( L.GT.MAXJ/2 ) THEN
*         WRITE(NB6,1010) L, MAXJ/2
*         STOP 'PSVXD'
*     ENDIF
C
      DO 200 K=0,L/2
          SCALE = ((-1)**K)*FACT(2*(L-K))*RFACT(K)*RFACT(L-K)/
     .            POW2(2*(L-K))
          DO 100 IP=0,L-2*K
              D(K,IP) = SCALE*RFACT(IP)*RFACT(L-2*K-IP)
  100     CONTINUE
  200 CONTINUE
C
      RETURN
*1000 FORMAT(' PSVXD CALLED WITH LDD = ', I3, ' WHILE AT LEAST ', I3, 
*    .       ' IS REQUIRED' )
*1010 FORMAT(' PSVXD CALLED WITH L = ', I5, 
*    .       ' WHILE COMPILE-TIME LIMIT WAS SET TO ', I3 )
      END
C
      SUBROUTINE PSVXQ(L1,N2,L2,R,DJ,LDJ,Q)
C
C   Transforms J integrals into Q integrals.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L1     - Orbital moment of the first center
C      N2     - Major quantum number of the second center
C      L2     - Orbital moment of the second center
C      R      - Intercenter distance
C      DJ     - Table of J integrals
C      LDJ    - Leading dimension of the DJ table
C      Q      - Space for computed Q integrals
C
C   Accessed common blocks:
C
C   Local storage:
C
C      (MAXJ/2+1)(MAXJ+1) DOUBLE PRECISION cells for the D coefficients.
C
C   Module logic:
C
C     Eq. (2.10) of Talman's paper. Summation over K was fixed to
C     start from MAX(0,(L-N2-LA+1)/2) instead of zero, as per 
C     eq. (2.10) itself.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXJ=30)
      PARAMETER (ZERO=0.0D0)
C
      DIMENSION DJ(0:LDJ-1,0:L1+L2),Q(0:L1,L2-L1:L2+L1)
      DIMENSION D(0:MAXJ/2,0:MAXJ)
C
      RR2 = R**(-2)
      DO 500 L=L2-L1,L2+L1
          CALL PSVXD(L,D,MAXJ/2+1)
          DO 300 LA=0,L1
              SUM   = ZERO
              DO 200 K=MAX(0,(L-N2-LA+1)/2),L/2
                  DO 100 IP=0,L-2*K
                      SUM = SUM + D(K,IP) * RR2**(K+IP) * 
     .                            DJ(N2+LA+2*K-L,IP)
  100             CONTINUE
  200         CONTINUE
              Q(LA,L) = SUM * (R**(L1-LA+L-1))
  300     CONTINUE
  500 CONTINUE
C
      END
C
      SUBROUTINE PSVXI(L1,L2,MINLB,X,Y,Z,A,OVR)
C
C   Transforms A(/\) integrals into overlap integrals.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      L1     - Orbital moment of the first center
C      L2     - Orbital moment of the second center
C      MINLB  - Should be equal to MAX(L2-L1,MOD(L2-L1,2))
C               Nope, to ABS(L2-L1)
C      X,Y,Z  - Relative coordinates of the second
C               center, normalized to 1.
C      A      - A(/\) integrals
C      OVR    - Space for desired overlap integrals
C
C   Accessed common blocks:
C
C   Local storage:
C
C      2*MAXJ+1 DOUBLE PRECISION cells for spherical harmonics.
C
C   Module logic:
C
C     Eq. (2.22) of Talman's paper, with conversion to real
C     harmonics and correction for a different definition of
C     the phase of spherical harmonics.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXJ=30)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (RSQ2=7.07106781186547524400844362104849039284835938D-1)
C
      DIMENSION A(MINLB:L2+L1)
      DIMENSION OVR(-MAXJ:MAXJ,-MAXJ:MAXJ)
      DIMENSION SL(-MAXJ:MAXJ)
C
      DO 200 M2=-L2,L2
          DO 100 M1=-L1,L1
              OVR(M1,M2) = ZERO
  100     CONTINUE
  200 CONTINUE
C
      DO 900 LB=MINLB,L2+L1,2
          CALL PSSLM(LB,X,Y,Z,SL(-LB))
          TEMP  = PS3J(L1,L2,LB,0,0,0)
          SCALE = A(LB) * TEMP
          OVR(0,0) = OVR(0,0) + SCALE * TEMP * SL(0)
C
          DO 300 M1=1,MIN(L1,LB)
              TEMP = SCALE * PS3J(L1,L2,LB,M1,0,-M1) * (-1)**M1
              OVR( M1,0) = OVR( M1,0) + TEMP * SL( M1)
              OVR(-M1,0) = OVR(-M1,0) + TEMP * SL(-M1)
  300     CONTINUE
C
          DO 400 M2=1,MIN(L2,LB)
              TEMP = SCALE * PS3J(L1,L2,LB,0,M2,-M2) * (-1)**M2
              OVR(0, M2) = OVR(0, M2) + TEMP * SL( M2)
              OVR(0,-M2) = OVR(0,-M2) + TEMP * SL(-M2)
  400     CONTINUE
C
          DO 600 M2=1,L2
              TEMP1 = RSQ2*SCALE*(-1)**M2
              DO 500 M1=MAX(1,M2-LB),MIN(L1,M2-1)
                  TEMP = TEMP1*PS3J(L1,L2,LB,M1,-M2,M2-M1)
                  TDP  = TEMP*SL(M2-M1)
                  TDM  = TEMP*SL(M1-M2)
                  OVR( M1, M2) = OVR( M1, M2) + TDP
                  OVR(-M1,-M2) = OVR(-M1,-M2) + TDP
                  OVR( M1,-M2) = OVR( M1,-M2) + TDM
                  OVR(-M1, M2) = OVR(-M1, M2) - TDM
  500         CONTINUE
              IF( M2.LE.L1 ) THEN
                  TDP  = SCALE*(-1)**M2*PS3J(L1,L2,LB,M2,-M2,0)*SL(0)
                  OVR( M1, M2) = OVR( M1, M2) + TDP
                  OVR(-M1,-M2) = OVR(-M1,-M2) + TDP
              ENDIF
              DO 550 M1=M2+1,MIN(L1,M2+LB)
                  TEMP = RSQ2*SCALE*PS3J(L1,L2,LB,M1,-M2,M2-M1)*
     .                        (-1)**M1
                  TDP  = TEMP*SL(M1-M2)
                  TDM  = TEMP*SL(M2-M1)
                  OVR( M1, M2) = OVR( M1, M2) + TDP
                  OVR(-M1,-M2) = OVR(-M1,-M2) + TDP
                  OVR( M1,-M2) = OVR( M1,-M2) - TDM
                  OVR(-M1, M2) = OVR(-M1, M2) + TDM
  550         CONTINUE
  600     CONTINUE
          DO 800 M2=1,L2
              DO 700 M1=1,MIN(L1,LB-M2)
                  TEMP = RSQ2*SCALE*(-1)**(M2+M1)*
     .                        PS3J(L1,L2,LB,M1,M2,-M1-M2)
                  TSP  = TEMP * SL(  M1+M2 )
                  TSM  = TEMP * SL(-(M1+M2))
                  OVR( M1, M2) = OVR( M1, M2) + TSP
                  OVR(-M1,-M2) = OVR(-M1,-M2) - TSP
                  OVR( M1,-M2) = OVR( M1,-M2) + TSM
                  OVR(-M1, M2) = OVR(-M1, M2) + TSM
  700         CONTINUE
  800     CONTINUE

  900 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSOLX(NORBS,L2,X,LDX,XL,LDL1,LDL2)
C
C                                      ^
C   Transform two-center integrals <mu|X|nu> into integrals
C       ^ ^              +
C   <mu|X L |nu> (sans i*h factor) for real Slater-type orbitals.
C          i
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NORBS  - Number of AOs at the first center
C      L2     - Orbital moment of the second set of orbitals
C                   ^
C      X      - <mu|X|nu> integrals, (2*L1+1)*(2*L2+1) entries
C      LDX    - Leading dimension of the X matrix
C                          ^ ^
C      XL     - Output <mu|X L |nu> integrals, 3*(2*L1+1)*(2*L2+1)
C                             i
C               entries.
C      LDL1   - First leading dimension of the XL matrix
C      LDL2   - Second leading dimension of the XL matrix
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Both input and output integrals should be in the order of
C      decreasing projection of the orbital moments.
C
C   Bugs:
C
C      Transformation is handled by explicit formulas for
C      different Ls, so L2 should not exceed 3.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IX=1)
      PARAMETER (IY=2)
      PARAMETER (IZ=3)
      PARAMETER (ZERO=0.D0)
      PARAMETER (TWO=2.0D0)
      PARAMETER (THREE=3.0D0)
      PARAMETER (SQ3=1.7320508075688772935274463415058723669428052538D0)
      PARAMETER (SQ1P5=1.22474487139158904909864203735294569598297374D0)
      PARAMETER (SQ2P5=1.58113883008418966599944677221635926685977757D0)
      PARAMETER (SQ6=2.4494897427831780981972840747058913919659474807D0)
C
      DIMENSION X(LDX,-L2:L2)
      DIMENSION XL(LDL1,-L2:LDL2-L2-1,3)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      GOTO (0010,1010,2010,3010), L2+1
      WRITE(NB6,11000) L2
      STOP 'PSOLX'
C
C  S-type orbitals, these are *really* easy ;-)
C
 0010 CONTINUE
          DO 900 I=1,3
              DO 100 M1=1,NORBS
                  XL(M1,0,I) = ZERO
  100         CONTINUE
  900     CONTINUE
          RETURN
C
C  P-type orbitals
C
 1010 CONTINUE
          DO 1100 M1=1,NORBS
              XL(M1,-1,IX) = -X(M1, 0)
              XL(M1,-1,IY) =  ZERO
              XL(M1,-1,IZ) =  X(M1, 1)
              XL(M1, 0,IX) =  X(M1,-1)
              XL(M1, 0,IY) = -X(M1, 1)
              XL(M1, 0,IZ) =  ZERO
              XL(M1, 1,IX) =  ZERO
              XL(M1, 1,IY) =  X(M1, 0)
              XL(M1, 1,IZ) = -X(M1,-1)
 1100     CONTINUE
          RETURN
C
C  D-type orbitals
C
 2010 CONTINUE
          DO 2100 M1=1,NORBS
              XL(M1,-2,IX) = -X(M1, 1)
              XL(M1,-2,IY) =  X(M1,-1)
              XL(M1,-2,IZ) =  X(M1, 2) * TWO
              XL(M1,-1,IX) = -X(M1, 0) * SQ3 - X(M1, 2)
              XL(M1,-1,IY) = -X(M1,-2)
              XL(M1,-1,IZ) =  X(M1, 1)
              XL(M1, 0,IX) =  X(M1,-1) * SQ3
              XL(M1, 0,IY) = -X(M1, 1) * SQ3
              XL(M1, 0,IZ) =  ZERO
              XL(M1, 1,IX) =  X(M1,-2)
              XL(M1, 1,IY) =  X(M1, 0) * SQ3 - X(M1, 2)
              XL(M1, 1,IZ) = -X(M1,-1)
              XL(M1, 2,IX) =  X(M1,-1)
              XL(M1, 2,IY) =  X(M1, 1)
              XL(M1, 2,IZ) = -X(M1,-2) * TWO
 2100     CONTINUE
          RETURN
C
C  F-type orbitals
C
 3010 CONTINUE
          DO 3100 M1=1,NORBS
              XL(M1,-3,IX) = -X(M1, 2) * SQ1P5
              XL(M1,-3,IY) =  X(M1,-2) * SQ1P5
              XL(M1,-3,IZ) =  X(M1, 3) * THREE
              XL(M1,-2,IX) = -X(M1, 1) * SQ2P5 - X(M1, 3) * SQ1P5
              XL(M1,-2,IY) = -X(M1,-3) * SQ1P5 + X(M1,-1) * SQ2P5
              XL(M1,-2,IZ) =  X(M1, 2) * TWO
              XL(M1,-1,IX) = -X(M1, 0) * SQ6   - X(M1, 2) * SQ2P5
              XL(M1,-1,IY) = -X(M1,-2) * SQ2P5
              XL(M1,-1,IZ) =  X(M1, 1)
              XL(M1, 0,IX) =  X(M1,-1) * SQ6
              XL(M1, 0,IY) = -X(M1, 1) * SQ6
              XL(M1, 0,IZ) =  ZERO
              XL(M1, 1,IX) =  X(M1,-2) * SQ2P5
              XL(M1, 1,IY) =  X(M1, 0) * SQ6   - X(M1, 2) * SQ2P5
              XL(M1, 1,IZ) = -X(M1,-1)
              XL(M1, 2,IX) =  X(M1,-3) * SQ1P5 + X(M1,-1) * SQ2P5
              XL(M1, 2,IY) =  X(M1, 1) * SQ2P5 - X(M1, 3) * SQ1P5
              XL(M1, 2,IZ) = -X(M1,-2) * TWO
              XL(M1, 3,IX) =  X(M1,-2) * SQ1P5
              XL(M1, 3,IY) =  X(M1, 2) * SQ1P5
              XL(M1, 3,IZ) = -X(M1,-3) * THREE
 3100     CONTINUE
          RETURN
C
11000 FORMAT(' ORBITAL MOMENT IN PSOLX (',I4,') NOT IN A VALID RANGE.')
      END
C
      SUBROUTINE PSODS(NORBS,N2,L2,ALP2,D,LDD,U,LDU,XD,LDD1,LDD2)
C                             
C   Compute integrals <n,l,m|O x|n',l',m'> using reduction to integrals
C   <n,l,m|O|n'+1,l'-1,?> and <n,l,m|O|n'+1,l'+1,?>.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NORBS  - Number of AOs at the first center
C      N2     - Major quantum number of STOs at the second center
C      L2     - Orbital moment of STOs at the second center
C      ALP2   - Orbital exponent of STOs at the second center
C      D      - <n,l,?|O|n'+1,l'-1,?> integrals. Not used if L2=0.
C      LDX    - Leading dimension of the X matrix
C      U      - <n,l,?|O|n'+1,l'+1,?> integrals.
C      LDY    - Leading dimension of the Y matrix
C      XD     - Output dipole moment integrals.
C      LDL1   - First leading dimension of the XL matrix
C      LDL2   - Second leading dimension of the XL matrix
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Both input and output integrals should be in the order of
C      decreasing projection of the orbital moments.
C
C   Bugs:
C
C      Transformation is handled by explicit formulas for
C      different Ls, so L2 should not exceed 3. Major quantum
C      numbers should not exceed 8.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXN=8)
      PARAMETER (IX=1)
      PARAMETER (IY=2)
      PARAMETER (IZ=3)
      PARAMETER (RSQ12 =2.88675134594812882254574390250978727823801D-1)
      PARAMETER (RSQ15 =2.58198889747161125678617693318826640722195D-1)
      PARAMETER (S9O140=2.53546276418554973252885498209783587789236D-1)
      PARAMETER (SQ4O63=2.51976315339484818143011024156120040543834D-1)
      PARAMETER (SQ5O84=2.43975018237133294838596159060025047398453D-1)
      PARAMETER (SQ2O35=2.39045721866878727993763435938624996969376D-1)
      PARAMETER (RSQ18 =2.35702260395515841466948120701616346428279D-1)
      PARAMETER (SQ3O56=2.31455024943137865391641694145999880632329D-1)
      PARAMETER (RSQ20 =2.23606797749978969640917366873127623544062D-1)
      PARAMETER (RSQ21 =2.18217890235992381266097485415619451856403D-1)
      PARAMETER (SQ3O70=2.07019667802706265338838059388139596044859D-1)
      PARAMETER (RSQ24 =2.04124145231931508183107006225490949330496D-1)
      PARAMETER (S5O126=1.99204768222398939994802863282187497474480D-1)
      PARAMETER (RSQ28 =1.88982236504613613607258268117090030407876D-1)
      PARAMETER (S5O168=1.72516389835588554449031716156782996704049D-1)
      PARAMETER (RSQ42 =1.54303349962091910261094462763999920421552D-1)
      PARAMETER (R6    =1.66666666666666666666666666666666666666667D-1)
      PARAMETER (S3O140=1.46385010942279976903157695436015028439072D-1)
      PARAMETER (RSQ60 =1.29099444873580562839308846659413320361097D-1)
      PARAMETER (RSQ168=7.71516749810459551305472313819999602107762D-2)
      PARAMETER (RSQ280=5.97614304667196819984408589846562492423440D-2)
      PARAMETER (RSQ504=4.45435403187373974474255801466255869256669D-2)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION D(LDD,-(L2-1):(L2-1)), U(LDU,-(L2+1):(L2+1))
      DIMENSION XD(LDD1,-L2:LDD2-L2-1,3)
      DIMENSION C(-1:MAXN)
      DATA C(-1)/0.0000000000000000000000000000000000000000000000000D0/
      DATA C( 0)/1.4142135623730950488016887242096980785696718753769D0/
      DATA C( 1)/3.4641016151377545870548926830117447338856105076208D0/
      DATA C( 2)/5.4772255750516611345696978280080213395274469499798D0/
      DATA C( 3)/7.4833147735478827711674974646330986035120396155575D0/
      DATA C( 4)/9.4868329805051379959966806332981556011586654179757D0/
      DATA C( 5)/11.489125293076057319701222936437858636440528915966D0/
      DATA C( 6)/13.490737563232041465550305611495544972455368100218D0/
      DATA C( 7)/15.491933384829667540717061599129598443331686821166D0/
      DATA C( 8)/17.492855684535901412622458632636749229564195004658D0/
C
      IF( LDD.LT.NORBS .OR. LDU.LT.NORBS .OR. LDD1.LT.NORBS 
     .                                   .OR. LDD2.LT.2*L2+1 ) THEN
          WRITE(NB6,11002) NORBS, L2, LDD, LDD1, LDD2
          STOP 'PSODS'
      ENDIF
      IF( N2.GT.MAXN ) THEN
          WRITE(NB6,11000) N2
          STOP 'PSODS'
      ENDIF
      SC = C(N2)/ALP2
      GOTO (0010,1010,2010,3010), L2+1
      WRITE(NB6,11010) L2
      STOP 'PSODS'
C
C   S-type orbitals
C
 0010 CONTINUE
          DO 0100 M1=1,NORBS
              XD(M1, 0,IX) =  U(M1, 1)*SC*RSQ12
              XD(M1, 0,IY) =  U(M1,-1)*SC*RSQ12
              XD(M1, 0,IZ) =  U(M1, 0)*SC*RSQ12
 0100     CONTINUE
          RETURN
C
C   P-type orbitals
C
 1010 CONTINUE
          DO 1100 M1=1,NORBS
              XD(M1,-1,IX) =  U(M1,-2)*SC*RSQ20
              XD(M1,-1,IY) = -U(M1, 0)*SC*RSQ60 -U(M1, 2)*SC*RSQ20  
     .                       +D(M1, 0)*SC*RSQ12
              XD(M1,-1,IZ) =  U(M1,-1)*SC*RSQ20
C
              XD(M1, 0,IX) =  U(M1, 1)*SC*RSQ20
              XD(M1, 0,IY) =  U(M1,-1)*SC*RSQ20
              XD(M1, 0,IZ) =  U(M1, 0)*SC*RSQ15                     
     .                       +D(M1, 0)*SC*RSQ12
C
              XD(M1, 1,IX) = -U(M1, 0)*SC*RSQ60 +U(M1, 2)*SC*RSQ20  
     .                       +D(M1, 0)*SC*RSQ12
              XD(M1, 1,IY) =  U(M1,-2)*SC*RSQ20
              XD(M1, 1,IZ) =  U(M1, 1)*SC*RSQ20
 1100     CONTINUE
          RETURN
C
C   D-type orbitals
C
 2010 CONTINUE
          DO 2100 M1=1,NORBS
              XD(M1,-2,IX) =  U(M1,-3)*SC*SQ3O56 -U(M1,-1)*SC*RSQ280 
     .                       +D(M1,-1)*SC*RSQ20
              XD(M1,-2,IY) = -U(M1, 1)*SC*RSQ280 -U(M1, 3)*SC*SQ3O56 
     .                       +D(M1, 1)*SC*RSQ20
              XD(M1,-2,IZ) =  U(M1,-2)*SC*RSQ28
C
              XD(M1,-1,IX) =  U(M1,-2)*SC*RSQ28
              XD(M1,-1,IY) = -U(M1, 0)*SC*S3O140 -U(M1, 2)*SC*RSQ28  
     .                       +D(M1, 0)*SC*RSQ20
              XD(M1,-1,IZ) =  U(M1,-1)*SC*SQ2O35                     
     .                       +D(M1,-1)*SC*RSQ20
C
              XD(M1, 0,IX) =  U(M1, 1)*SC*SQ3O70                     
     .                       -D(M1, 1)*SC*RSQ60
              XD(M1, 0,IY) =  U(M1,-1)*SC*SQ3O70                     
     .                       -D(M1,-1)*SC*RSQ60
              XD(M1, 0,IZ) =  U(M1, 0)*SC*S9O140                     
     .                       +D(M1, 0)*SC*RSQ15
C
              XD(M1, 1,IX) = -U(M1, 0)*SC*S3O140 +U(M1, 2)*SC*RSQ28  
     .                       +D(M1, 0)*SC*RSQ20
              XD(M1, 1,IY) =  U(M1,-2)*SC*RSQ28
              XD(M1, 1,IZ) =  U(M1, 1)*SC*SQ2O35                     
     .                       +D(M1, 1)*SC*RSQ20
C
              XD(M1, 2,IX) = -U(M1, 1)*SC*RSQ280 +U(M1, 3)*SC*SQ3O56 
     .                       +D(M1, 1)*SC*RSQ20
              XD(M1, 2,IY) =  U(M1,-3)*SC*SQ3O56 +U(M1,-1)*SC*RSQ280 
     .                       -D(M1,-1)*SC*RSQ20
              XD(M1, 2,IZ) =  U(M1, 2)*SC*RSQ28
 2100     CONTINUE
          RETURN
C
C    F-type orbitals
C
 3010 CONTINUE
          DO 3100 M1=1,NORBS
              XD(M1,-3,IX) =  U(M1,-4)*SC*RSQ18  -U(M1,-2)*SC*RSQ504 
     .                       +D(M1,-2)*SC*SQ3O56
              XD(M1,-3,IY) = -U(M1, 4)*SC*RSQ18  -U(M1, 2)*SC*RSQ504 
     .                       +D(M1, 2)*SC*SQ3O56
              XD(M1,-3,IZ) =  U(M1,-3)*SC*R6
C
              XD(M1,-2,IX) =  U(M1,-3)*SC*RSQ24  -U(M1,-1)*SC*RSQ168 
     .                       +D(M1,-1)*SC*RSQ28
              XD(M1,-2,IY) = -U(M1, 3)*SC*RSQ24  -U(M1, 1)*SC*RSQ168 
     .                       +D(M1, 1)*SC*RSQ28
              XD(M1,-2,IZ) =  U(M1,-2)*SC*RSQ21                      
     .                       +D(M1,-2)*SC*RSQ28
C
              XD(M1,-1,IX) =  U(M1,-2)*SC*S5O168                     
     .                       -D(M1,-2)*SC*RSQ280
              XD(M1,-1,IY) = -U(M1, 2)*SC*S5O168 -U(M1, 0)*SC*RSQ42  
     .                       +D(M1, 2)*SC*RSQ280 +D(M1, 0)*SC*SQ3O70
              XD(M1,-1,IZ) =  U(M1,-1)*SC*SQ5O84                     
     .                       +D(M1,-1)*SC*SQ2O35
C
              XD(M1, 0,IX) =  U(M1, 1)*SC*S5O126                     
     .                       -D(M1, 1)*SC*S3O140
              XD(M1, 0,IY) =  U(M1,-1)*SC*S5O126                     
     .                       -D(M1,-1)*SC*S3O140
              XD(M1, 0,IZ) =  U(M1, 0)*SC*SQ4O63                     
     .                       +D(M1, 0)*SC*S9O140
C
              XD(M1, 1,IX) =  U(M1, 2)*SC*S5O168 -U(M1, 0)*SC*RSQ42  
     .                       -D(M1, 2)*SC*RSQ280 +D(M1, 0)*SC*SQ3O70
              XD(M1, 1,IY) =  U(M1,-2)*SC*S5O168                     
     .                       -D(M1,-2)*SC*RSQ280
              XD(M1, 1,IZ) =  U(M1, 1)*SC*SQ5O84                     
     .                       +D(M1, 1)*SC*SQ2O35
C
              XD(M1, 2,IX) =  U(M1, 3)*SC*RSQ24  -U(M1, 1)*SC*RSQ168 
     .                       +D(M1, 1)*SC*RSQ28
              XD(M1, 2,IY) =  U(M1,-3)*SC*RSQ24  +U(M1,-1)*SC*RSQ168 
     .                       -D(M1,-1)*SC*RSQ28
              XD(M1, 2,IZ) =  U(M1, 2)*SC*RSQ21                      
     .                       +D(M1, 2)*SC*RSQ28
C
              XD(M1, 3,IX) =  U(M1, 4)*SC*RSQ18  -U(M1, 2)*SC*RSQ504 
     .                       +D(M1, 2)*SC*SQ3O56
              XD(M1, 3,IY) =  U(M1,-4)*SC*RSQ18  +U(M1,-2)*SC*RSQ504 
     .                       -D(M1,-2)*SC*SQ3O56
              XD(M1, 3,IZ) =  U(M1, 3)*SC*R6
 3100     CONTINUE
          RETURN

11000 FORMAT(' MAJOR QUANTUM NUMBER IS TOO LARGE IN PSODS: ',I3)
11002 FORMAT(' WRONG LEADING DIMENSIONS IN PSODS:'/
     .       ' NORBS = ',I3,' L2 = ',I3,' LDD = ',I3,' LDD1 = ',I3,
     .       ' LDD2 = ',I3)
11010 FORMAT(' ORBITAL QUANTUM NUMBER IS TOO LARGE IN PSODS: ',I3)
      END
C
      SUBROUTINE PSOQS(NORBS,N2,L2,ALP2,AX,LDA1,LDA2,XD,LDX)
C                             
C   Compute integrals <n,l,m|O r^2 \delta_{ab} - r_a r_b |n',l',m'> 
C   using reduction to integrals <n,l,m|O|n'+2,l'-2,?>, <n,l,m|O|n'+2,l',?>
C   and <n,l,m|O|n'+2,l+2,?>.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NORBS  - Number of AOs at the first center
C      N2     - Major quantum number of STOs at the second center
C      L2     - Orbital moment of STOs at the second center
C      ALP2   - Orbital exponent of STOs at the second center
C      AX     - Input auxiliary integrals.
C               First index enumerates left orbitals
C               Second index enumerates magnetic quantum numbers of the
C                   orbitals on the second center (in the increasing order)
C               Third index enumerates possible values of the orbital
C                   quantum number of the orbitals at the second center
C                   (1 = l'-2, 2 = l', 3 = l'+2)
C      LDA1   - First leading dimension of the AX array
C      LDA2   - Second leading dimension of the AX array
C      XD     - Output ontegrals.
C      LDX    - First leading dimension of the XD array.
C
C   Accessed common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Bulk of the code is computer-generated, do not look for logic...
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (C1=0.16666666666666666666666666666666666666667D0)
      PARAMETER (C2=0.03726779962499649494015289447885460392401D0)
      PARAMETER (C3=0.06454972243679028141965442332970666018055D0)
      PARAMETER (C4=0.07453559924999298988030578895770920784802D0)
      PARAMETER (C5=0.2D0)
      PARAMETER (C6=0.05175491695067656633470951484703489901121D0)
      PARAMETER (C7=0.0133630620956212192342276740439876760777D0)
      PARAMETER (C8=0.05D0)
      PARAMETER (C9=0.04225771273642582887548091636829726463154D0)
      PARAMETER (C10=0.1D0)
      PARAMETER (C11=0.0400891862868636577026830221319630282331D0)
      PARAMETER (C12=0.03273268353539885718991462281234291777846D0)
      PARAMETER (C13=0.0534522483824848769369106961759507043108D0)
      PARAMETER (C14=0.06546536707079771437982924562468583555692D0)
      PARAMETER (C15=0.14285714285714285714285714285714285714286D0)
      PARAMETER (C16=0.0545544725589980953165243713539048629641D0)
      PARAMETER (C17=0.02061965247105806301818388501792705198741D0)
      PARAMETER (C18=0.04123930494211612603636777003585410397483D0)
      PARAMETER (C19=0.009221388919541468774236346189958094311507D0)
      PARAMETER (C20=0.03571428571428571428571428571428571428571D0)
      PARAMETER (C21=0.03857583749052297756527361569099998010539D0)
      PARAMETER (C22=0.014580296087995107727364786158963639237893D0)
      PARAMETER (C23=0.21428571428571428571428571428571428571429D0)
      PARAMETER (C24=0.04374088826398532318209435847689091771368D0)
      PARAMETER (C25=0.03688555567816587509694538475983237724603D0)
      PARAMETER (C26=0.05832118435198043090945914463585455695157D0)
      PARAMETER (C27=0.19047619047619047619047619047619047619048D0)
      PARAMETER (C28=0.03194382824999699566298819526758966050629D0)
      PARAMETER (C29=0.05050762722761053745720316872177493137749D0)
      PARAMETER (C30=0.11904761904761904761904761904761904761905D0)
      PARAMETER (C31=0.06388765649999399132597639053517932101259D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION AX(LDA1,0:LDA2-1,3), XD(LDX,3,3,-L2:L2)
C
      GOTO (0010,1010,2010), L2+1
      WRITE(NB6,11010) L2
      STOP 'PSOQS'
C  L = 0
 0010 CONTINUE
      DO 0090 M1=1,NORBS
          XD(M1,1,1,0) = C1*AX(M1,0,2) + C2*AX(M1,2,3) - C3*AX(M1,4,3)
          XD(M1,2,1,0) = -(C3*AX(M1,0,3))
          XD(M1,3,1,0) = -(C3*AX(M1,3,3))
          XD(M1,1,2,0) = -(C3*AX(M1,0,3))
          XD(M1,2,2,0) = C1*AX(M1,0,2) + C2*AX(M1,2,3) + C3*AX(M1,4,3)
          XD(M1,3,2,0) = -(C3*AX(M1,1,3))
          XD(M1,1,3,0) = -(C3*AX(M1,3,3))
          XD(M1,2,3,0) = -(C3*AX(M1,1,3))
          XD(M1,3,3,0) = C1*AX(M1,0,2) - C4*AX(M1,2,3)
 0090 CONTINUE
      GOTO 9000
C  L = 1
 1010 CONTINUE
      DO 1090 M1=1,NORBS
          XD(M1,1,1,-1) = C5*AX(M1,0,2) - C6*AX(M1,0,3) + C7*AX(M1,2,3)
          XD(M1,2,1,-1) = -(C8*AX(M1,2,2)) + C7*AX(M1,4,3) + 
     .       C6*AX(M1,6,3)
          XD(M1,3,1,-1) = -(C9*AX(M1,1,3))
          XD(M1,1,2,-1) = -(C8*AX(M1,2,2)) + C7*AX(M1,4,3) + 
     .       C6*AX(M1,6,3)
          XD(M1,2,2,-1) = C10*AX(M1,0,2) + C6*AX(M1,0,3) + 
     .       C11*AX(M1,2,3)
          XD(M1,3,2,-1) = -(C8*AX(M1,1,2)) + C12*AX(M1,3,3) + 
     .       C9*AX(M1,5,3)
          XD(M1,1,3,-1) = -(C9*AX(M1,1,3))
          XD(M1,2,3,-1) = -(C8*AX(M1,1,2)) + C12*AX(M1,3,3) + 
     .       C9*AX(M1,5,3)
          XD(M1,3,3,-1) = C5*AX(M1,0,2) - C13*AX(M1,2,3)
          XD(M1,1,1,0) = C5*AX(M1,1,2) + C12*AX(M1,3,3) - C9*AX(M1,5,3)
          XD(M1,2,1,0) = -(C9*AX(M1,1,3))
          XD(M1,3,1,0) = -(C8*AX(M1,2,2)) - C13*AX(M1,4,3)
          XD(M1,1,2,0) = -(C9*AX(M1,1,3))
          XD(M1,2,2,0) = C5*AX(M1,1,2) + C12*AX(M1,3,3) + C9*AX(M1,5,3)
          XD(M1,3,2,0) = -(C8*AX(M1,0,2)) - C13*AX(M1,2,3)
          XD(M1,1,3,0) = -(C8*AX(M1,2,2)) - C13*AX(M1,4,3)
          XD(M1,2,3,0) = -(C8*AX(M1,0,2)) - C13*AX(M1,2,3)
          XD(M1,3,3,0) = C10*AX(M1,1,2) - C14*AX(M1,3,3)
          XD(M1,1,1,1) = C10*AX(M1,2,2) + C11*AX(M1,4,3) - C6*AX(M1,6,3)
          XD(M1,2,1,1) = -(C8*AX(M1,0,2)) - C6*AX(M1,0,3) + 
     .       C7*AX(M1,2,3)
          XD(M1,3,1,1) = -(C8*AX(M1,1,2)) + C12*AX(M1,3,3) - 
     .       C9*AX(M1,5,3)
          XD(M1,1,2,1) = -(C8*AX(M1,0,2)) - C6*AX(M1,0,3) + 
     .       C7*AX(M1,2,3)
          XD(M1,2,2,1) = C5*AX(M1,2,2) + C7*AX(M1,4,3) + C6*AX(M1,6,3)
          XD(M1,3,2,1) = -(C9*AX(M1,1,3))
          XD(M1,1,3,1) = -(C8*AX(M1,1,2)) + C12*AX(M1,3,3) - 
     .       C9*AX(M1,5,3)
          XD(M1,2,3,1) = -(C9*AX(M1,1,3))
          XD(M1,3,3,1) = C5*AX(M1,2,2) - C13*AX(M1,4,3)
 1090 CONTINUE
      GOTO 9000
C  L = 2
 2010 CONTINUE
      DO 2090 M1=1,NORBS
          XD(M1,1,1,-2) = C15*AX(M1,0,2) - C16*AX(M1,0,3) + 
     .       C17*AX(M1,2,3)
          XD(M1,2,1,-2) = -(C3*AX(M1,0,1)) + C18*AX(M1,2,2) - 
     .       C19*AX(M1,4,3) + C16*AX(M1,8,3)
          XD(M1,3,1,-2) = -(C20*AX(M1,1,2)) - C21*AX(M1,1,3) + 
     .       C22*AX(M1,3,3)
          XD(M1,1,2,-2) = -(C3*AX(M1,0,1)) + C18*AX(M1,2,2) - 
     .       C19*AX(M1,4,3) + C16*AX(M1,8,3)
          XD(M1,2,2,-2) = C15*AX(M1,0,2) + C16*AX(M1,0,3) + 
     .       C17*AX(M1,2,3)
          XD(M1,3,2,-2) = -(C20*AX(M1,3,2)) + C22*AX(M1,5,3) + 
     .       C21*AX(M1,7,3)
          XD(M1,1,3,-2) = -(C20*AX(M1,1,2)) - C21*AX(M1,1,3) + 
     .       C22*AX(M1,3,3)
          XD(M1,2,3,-2) = -(C20*AX(M1,3,2)) + C22*AX(M1,5,3) + 
     .       C21*AX(M1,7,3)
          XD(M1,3,3,-2) = C23*AX(M1,0,2) - C18*AX(M1,2,3)
          XD(M1,1,1,-1) = C23*AX(M1,1,2) - C21*AX(M1,1,3) + 
     .       C22*AX(M1,3,3)
          XD(M1,2,1,-1) = -(C20*AX(M1,3,2)) + C22*AX(M1,5,3) + 
     .       C21*AX(M1,7,3)
          XD(M1,3,1,-1) = -(C20*AX(M1,0,2)) - C18*AX(M1,2,3)
          XD(M1,1,2,-1) = -(C20*AX(M1,3,2)) + C22*AX(M1,5,3) + 
     .       C21*AX(M1,7,3)
          XD(M1,2,2,-1) = C15*AX(M1,1,2) + C21*AX(M1,1,3) + 
     .       C24*AX(M1,3,3)
          XD(M1,3,2,-1) = -(C3*AX(M1,0,1)) - C17*AX(M1,2,2) + 
     .       C20*AX(M1,4,2) + C25*AX(M1,4,3) + C18*AX(M1,6,3)
          XD(M1,1,3,-1) = -(C20*AX(M1,0,2)) - C18*AX(M1,2,3)
          XD(M1,2,3,-1) = -(C3*AX(M1,0,1)) - C17*AX(M1,2,2) + 
     .       C20*AX(M1,4,2) + C25*AX(M1,4,3) + C18*AX(M1,6,3)
          XD(M1,3,3,-1) = C15*AX(M1,1,2) - C26*AX(M1,3,3)
          XD(M1,1,1,0) = C2*AX(M1,0,1) + C27*AX(M1,2,2) + 
     .       C18*AX(M1,4,2) + C28*AX(M1,4,3) - C20*AX(M1,6,3)
          XD(M1,2,1,0) = C18*AX(M1,0,2) - C20*AX(M1,2,3)
          XD(M1,3,1,0) = -(C17*AX(M1,3,2)) - C29*AX(M1,5,3)
          XD(M1,1,2,0) = C18*AX(M1,0,2) - C20*AX(M1,2,3)
          XD(M1,2,2,0) = C2*AX(M1,0,1) + C27*AX(M1,2,2) - 
     .       C18*AX(M1,4,2) + C28*AX(M1,4,3) + C20*AX(M1,6,3)
          XD(M1,3,2,0) = -(C17*AX(M1,1,2)) - C29*AX(M1,3,3)
          XD(M1,1,3,0) = -(C17*AX(M1,3,2)) - C29*AX(M1,5,3)
          XD(M1,2,3,0) = -(C17*AX(M1,1,2)) - C29*AX(M1,3,3)
          XD(M1,3,3,0) = -(C4*AX(M1,0,1)) + C30*AX(M1,2,2) - 
     .       C31*AX(M1,4,3)
          XD(M1,1,1,1) = C15*AX(M1,3,2) + C24*AX(M1,5,3) - 
     .       C21*AX(M1,7,3)
          XD(M1,2,1,1) = -(C20*AX(M1,1,2)) - C21*AX(M1,1,3) + 
     .       C22*AX(M1,3,3)
          XD(M1,3,1,1) = -(C3*AX(M1,0,1)) - C17*AX(M1,2,2) - 
     .       C20*AX(M1,4,2) + C25*AX(M1,4,3) - C18*AX(M1,6,3)
          XD(M1,1,2,1) = -(C20*AX(M1,1,2)) - C21*AX(M1,1,3) + 
     .       C22*AX(M1,3,3)
          XD(M1,2,2,1) = C23*AX(M1,3,2) + C22*AX(M1,5,3) + 
     .       C21*AX(M1,7,3)
          XD(M1,3,2,1) = -(C20*AX(M1,0,2)) - C18*AX(M1,2,3)
          XD(M1,1,3,1) = -(C3*AX(M1,0,1)) - C17*AX(M1,2,2) - 
     .       C20*AX(M1,4,2) + C25*AX(M1,4,3) - C18*AX(M1,6,3)
          XD(M1,2,3,1) = -(C20*AX(M1,0,2)) - C18*AX(M1,2,3)
          XD(M1,3,3,1) = C15*AX(M1,3,2) - C26*AX(M1,5,3)
          XD(M1,1,1,2) = -(C3*AX(M1,0,1)) + C18*AX(M1,2,2) + 
     .       C15*AX(M1,4,2) - C19*AX(M1,4,3) + C17*AX(M1,6,3) - 
     .       C16*AX(M1,8,3)
          XD(M1,2,1,2) = -(C16*AX(M1,0,3))
          XD(M1,3,1,2) = -(C20*AX(M1,3,2)) + C22*AX(M1,5,3) - 
     .       C21*AX(M1,7,3)
          XD(M1,1,2,2) = -(C16*AX(M1,0,3))
          XD(M1,2,2,2) = C3*AX(M1,0,1) - C18*AX(M1,2,2) + 
     .       C15*AX(M1,4,2) + C19*AX(M1,4,3) + C17*AX(M1,6,3) + 
     .       C16*AX(M1,8,3)
          XD(M1,3,2,2) = C20*AX(M1,1,2) - C21*AX(M1,1,3) - 
     .       C22*AX(M1,3,3)
          XD(M1,1,3,2) = -(C20*AX(M1,3,2)) + C22*AX(M1,5,3) - 
     .       C21*AX(M1,7,3)
          XD(M1,2,3,2) = C20*AX(M1,1,2) - C21*AX(M1,1,3) - 
     .       C22*AX(M1,3,3)
          XD(M1,3,3,2) = C23*AX(M1,4,2) - C18*AX(M1,6,3)
 2090 CONTINUE
      GOTO 9000
C
C    Final scaling of the integrals.
C
 9000 CONTINUE
      SC = SQRT(DBLE(2*N2+1)*DBLE(2*N2+2)*DBLE(2*N2+3)*DBLE(2*N2+4))
     .         /ALP2**2
      DO 9500 M2=-L2,L2
          DO 9400 IY=1,3
              DO 9300 IX=1,3
                  DO 9200 M1=1,NORBS
                      XD(M1,IX,IY,M2) = SC*XD(M1,IX,IY,M2)
 9200             CONTINUE
 9300         CONTINUE
 9400     CONTINUE
 9500 CONTINUE
      RETURN
C
11010 FORMAT(' ORBITAL QUANTUM NUMBER IS TOO LARGE IN PSOQS: ',I3)
      END
C
      SUBROUTINE PSODX(N1,L1,ALP1,N2,L2,ALP2,X,Y,Z,DD,LD1,LD2)
C
C   Compute dipole moment integrals over block of Slater AOs.
C   (i.e., <\mu|r-R_\nu|\nu> integrals)
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N1     - Major quantum number of the first AO
C      L1     - Orbital quantum number of the first AO
C      ALP1   - Orbital exponent of the first AO
C      N2,L2,ALP2 
C             - Same for the second atom
C      X,Y,Z  - Relative coordinates of the second atom.
C      DD     - Computed dipole moment integrals
C      LD1    - First leading dimension of D
C      LD2    - Second leading dimension of D
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
      PARAMETER (MAXL=5)
C
      DIMENSION DD(-L1:LD1-L1-1,-L2:LD2-L2-1,3)
      DIMENSION U(-MAXL:MAXL,-(MAXL+1):(MAXL+1)),
     .          D(-MAXL:MAXL,-MAXL:MAXL)
C
      CALL PSOVX(N1,L1,ALP1,N2+1,L2+1,ALP2,
     .           X,Y,Z,U(-L1,-(L2+1)),2*MAXL+1)
      IF(L2.GE.1)
     .CALL PSOVX(N1,L1,ALP1,N2+1,L2-1,ALP2,
     .           X,Y,Z,D(-L1,-(L2-1)),2*MAXL+1)
      CALL PSODS(2*L1+1,N2,L2,ALP2,D(-L1,-(L2-1)),2*MAXL+1,
     .                             U(-L1,-(L2+1)),2*MAXL+1,DD,LD1,LD2)
C
      RETURN
      END
C
      SUBROUTINE PSOQR3(N1,L1,ALP1,N2,L2,ALP2,X,Y,Z,DD,LD)
C
C   Compute magnetic dipole interactiopn integrals over block of Slater 
C   AOs.  (i.e., <\mu|(\delta_{ab} r^2 - r_b r_a) r^-3|\nu> integrals,
C   with the operator centered at the right atom).
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N1     - Major quantum number of the first AO
C      L1     - Orbital quantum number of the first AO
C      ALP1   - Orbital exponent of the first AO
C      N2,L2,ALP2 
C             - Same for the second atom
C      X,Y,Z  - Relative coordinates of the second atom.
C      DD     - Computed quadrupole moment integrals.
C               Two middle indices are components of the
C               quadrupole moment operator.
C      LD     - First leading dimension of D
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
      PARAMETER (MAXL=5)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION DD(LD,3,3,2*L2+1)
      DIMENSION AX(2*MAXL+1,2*(MAXL+2)+1,3)
C
      IF( L1.GT.MAXL .OR. L2.GT.MAXL ) THEN
          WRITE(NB6,11000) N1,L1,N2,L2
          STOP 'PSOQR3'
      ENDIF
C
      IF(L2.GE.2)
     .CALL PSORX(N1,L1,ALP1,N2+2,L2-2,ALP2,X,Y,Z,3,AX(1,1,1),2*MAXL+1)
      CALL PSORX(N1,L1,ALP1,N2+2,L2+0,ALP2,X,Y,Z,3,AX(1,1,2),2*MAXL+1)
      CALL PSORX(N1,L1,ALP1,N2+2,L2+2,ALP2,X,Y,Z,3,AX(1,1,3),2*MAXL+1)
      CALL PSOQS(2*L1+1,N2,L2,ALP2,AX,2*MAXL+1,2*(MAXL+2)+1,DD,LD)
      RETURN
C
11000 FORMAT(' STATIC TEMPORARY WILL OVERFLOW IN PSOQR3. N1 = ',I5,
     .       ' L1 = ',I5,' N2 = ',I5,' L2 = ',I5)
      END
C
      SUBROUTINE PSOVP(N1I,L1I,ALP1I,N2I,L2I,ALP2I,AP,LDA,XP,LDX)
C                             
C   Compute integrals <n,l,m|O|n',l',m'> with |nlm> and |n'l'm'> 
C   centered on the same atom using reduction to integrals <1|O|NLM>.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N1I    - Major quantum number of the first set of orbitals
C      L1I    - Orbital quantum number of the first set of orbitals
C      ALP1I  - Orbital exponent of the first set of orbitals
C      N2I    - Major quantum number of the second set of orbitals
C      L2I    - Orbital quantum number of the second set of orbitals
C      ALP12I - Orbital exponent of the second set of orbitals
C      AP     - Input <1|O|NLM> integrals
C      LDA    - Leading dimensions of the AP array. At least MIN(L1,L2)+1
C               entries by the second index should be filled in, with
C               2*(L1+L2)+1 values for L=L1+L2 in the first set, 2*(L1+L2)-3
C               values for L=L1+L2-2 in the second one etc.
C      XP     - Output <n,l,m|O|n',l',m'> integrals
C      LDX    - Leading dimension of the XP array.
C
C   Accessed common blocks:
C
C      PS3JDT - Factorial falues.
C
C   Local storage:
C
C   Module logic:
C
C      Logic be damned, enjoy the code!
C
C   Bugs:
C
C      L1 and L2 should not exceed 4.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXJ=30)
      PARAMETER (MXJ=90)
      PARAMETER (MAXL=4)
      PARAMETER (EIGHT=8.D0)
      PARAMETER (C1=0.282094791773878143474039725780386292922D0)
      PARAMETER (C2=0.1261566261010080024123574761182841973297D0)
      PARAMETER (C3=0.2185096861184158141086771411605376805381D0)
      PARAMETER (C4=0.2523132522020160048247149522365683946594D0)
      PARAMETER (C5=0.05839917008190184933692082809716963379762D0)
      PARAMETER (C6=0.2261790131595402924439377371777627807359D0)
      PARAMETER (C7=0.14304816810266883076266069110973459451098D0)
      PARAMETER (C8=0.1846743909223718017185109504744334766161D0)
      PARAMETER (C9=0.202300659403420632118494952299941262858D0)
      PARAMETER (C10=0.2335966803276073973476833123886785351905D0)
      PARAMETER (C11=0.24776669508347606189051829894840077192D0)
      PARAMETER (C12=0.04029925596769687763914853225434089898886D0)
      PARAMETER (C13=0.2384136135044480512711011518495576575183D0)
      PARAMETER (C14=0.180223751572868574874796394454691710471D0)
      PARAMETER (C15=0.06371871843402753980496308473446723537122D0)
      PARAMETER (C16=0.168583882836183860098745793583284385715D0)
      PARAMETER (C17=0.15607834722743986722048367225752691467006D0)
      PARAMETER (C18=0.1611970238707875105565941290173635959554D0)
      PARAMETER (C19=0.0901118757864342874373981972273458552355D0)
      PARAMETER (C20=0.2207281154418226173383958487494222360492D0)
      PARAMETER (C21=0.2417955358061812658348911935260453939332D0)
      PARAMETER (C22=0.04352817137756816534838292351587129669398D0)
      PARAMETER (C23=0.2303294329808903195101630973545931220929D0)
      PARAMETER (C24=0.07539300438651343081464591239258759357863D0)
      PARAMETER (C25=0.1994711402007163389699730299671909342379D0)
      PARAMETER (C26=0.15078600877302686162929182478517518715727D0)
      PARAMETER (C27=0.1946639002730061644564027603238987793254D0)
      PARAMETER (C28=0.1628675039676399738621282076127823348586D0)
      PARAMETER (C29=0.21324361862292308080822623611371272443013D0)
      PARAMETER (C30=0.2462325212298290689580146006325779688215D0)
      PARAMETER (C31=0.01694331772935932078559337544044873695741D0)
      PARAMETER (C32=0.2455320005465369006554142777588680492204D0)
      PARAMETER (C33=0.09403159725795938115801324192679543097401D0)
      PARAMETER (C34=0.05357947514468781075417011319294576064346D0)
      PARAMETER (C35=0.1901882698155455593956144712399123887117D0)
      PARAMETER (C36=0.188063194515918762316026483853590861948D0)
      PARAMETER (C37=0.06562118739530951484050111379824072166593D0)
      PARAMETER (C38=0.14175796661021041992692750556397343297754D0)
      PARAMETER (C39=0.14567312407894387607245142744035845369206D0)
      PARAMETER (C40=0.04482780509623634505574903692798520626119D0)
      PARAMETER (C41=0.15528807203695278910892675531973236908691D0)
      PARAMETER (C42=0.14867700967939759287824719771375688253285D0)
      PARAMETER (C43=0.08300496597356404770566062579330275396765D0)
      PARAMETER (C44=0.1793112203849453802229961477119408250448D0)
      PARAMETER (C45=0.11516471649044515975508154867729656104646D0)
      PARAMETER (C46=0.1694331772935932078559337544044873695741D0)
      PARAMETER (C47=0.1736173425847553353833100038948186208218D0)
      PARAMETER (C48=0.05947080387175903715129887908550275301314D0)
      PARAMETER (C49=0.2143179005787512430166804527717830425739D0)
      PARAMETER (C50=0.12679217987703037293040964749327492580777D0)
      PARAMETER (C51=0.2102610435016800040205957935304736622162D0)
      PARAMETER (C52=0.2273184612433489533555389667103530941543D0)
      PARAMETER (C53=0.2396146972445646498119247933731831767992D0)
      PARAMETER (C54=0.1682088348013440032164766348243789297729D0)
      PARAMETER (C55=0.011854396693264042131995716747889298489663D0)
      PARAMETER (C56=0.2548005986729750384987625828534773006668D0)
      PARAMETER (C57=0.07693494321105767549291992521283262534237D0)
      PARAMETER (C58=0.02217754547654999394672057123525469498494D0)
      PARAMETER (C59=0.1801712311720526720645076450295192895553D0)
      PARAMETER (C60=0.09932258459927991550394415507297167297186D0)
      PARAMETER (C61=0.04435509095309998789344114247050938996987D0)
      PARAMETER (C62=0.12147141927603090677236274273186672740847D0)
      PARAMETER (C63=0.13325523051897816043265035168318554670251D0)
      PARAMETER (C64=0.117520066950600237457234241211003121483D0)
      PARAMETER (C65=0.10864734032983335922332467636165953431551D0)
      PARAMETER (C66=0.2035507268673356680442159526539985096833D0)
      PARAMETER (C67=0.07112638015958425279197430048733579093798D0)
      PARAMETER (C68=0.1881827135584985209340713696226017702959D0)
      PARAMETER (C69=0.1795148674924679094834798254966094591322D0)
      PARAMETER (C70=0.15171775404828512353615527844971850932983D0)
      PARAMETER (C71=0.0858932642904357574957170905569918013544D0)
      PARAMETER (C72=0.1629710104947500388349870145424893014733D0)
      PARAMETER (C73=0.10257992428141023399055990028377683378983D0)
      PARAMETER (C74=0.06785024228911188934807198421799950322777D0)
      PARAMETER (C75=0.1774203638123999515737645698820375598795D0)
      PARAMETER (C76=0.04441841017299272014421678389439518223417D0)
      PARAMETER (C77=0.1778159503989606319799357512183394773449D0)
      PARAMETER (C78=0.1717865285808715149914341811139836027088D0)
      PARAMETER (C79=0.02564498107035255849763997507094420844746D0)
      PARAMETER (C80=0.11468784191000727492032497828934927029973D0)
      PARAMETER (C81=0.2217754547654999394672057123525469498494D0)
      PARAMETER (C82=0.2370879338652808426399143349577859697933D0)
      PARAMETER (C83=0.15386988642211535098583985042566525068474D0)
      PARAMETER (C84=0.03472346851695106707666200077896372416437D0)
      PARAMETER (C85=0.2329321080554291836633901329795985536304D0)
      PARAMETER (C86=0.06014281168637758170145216340113182955448D0)
      PARAMETER (C87=0.2083408111017064024599720046737823449862D0)
      PARAMETER (C88=0.08505477996612625195615650333838405978652D0)
      PARAMETER (C89=0.1837393247068666415534031186350740206646D0)
      PARAMETER (C90=0.159122922870344267348877276697247165908D0)
      PARAMETER (C91=0.14731920032792214039324856665532082953225D0)
      PARAMETER (C92=0.1964256004372295205243314222070944393763D0)
      PARAMETER (C93=0.2250337956076888051675144754103721947025D0)
      PARAMETER (C94=0.2405712467455103268058086536045273182179D0)
      PARAMETER (C95=0.011246068151364836493531155471437635112813D0)
      PARAMETER (C96=0.2502092208297798582872911052372942560252D0)
      PARAMETER (C97=0.07508081669196244299712675725709114300607D0)
      PARAMETER (C98=0.02514697286604716871419349298197758792349D0)
      PARAMETER (C99=0.2042949733241056242379415056776767024939D0)
      PARAMETER (C100=0.1126212250379436644956901358856367145091D0)
      PARAMETER (C101=0.06159725209742943823820981850790253133033D0)
      PARAMETER (C102=0.190364615027111657760929073833122209575D0)
      PARAMETER (C103=0.06653263642964998184016171370576408495481D0)
      PARAMETER (C104=0.12623680191215221204241927494147718357323D0)
      PARAMETER (C105=0.14188940657039987929134879296138810424551D0)
      PARAMETER (C106=0.13306527285929996368032342741152816990962D0)
      PARAMETER (C107=0.03373820445409450948059346641431290533844D0)
      PARAMETER (C108=0.14445836099979991918003483812542584125215D0)
      PARAMETER (C109=0.14046334619025075648807008556781415021313D0)
      PARAMETER (C110=0.06361736841212909169929946521178703545452D0)
      PARAMETER (C111=0.174223338642198573343716142640687294166D0)
      PARAMETER (C112=0.13272538654977693031118370022084007681277D0)
      PARAMETER (C113=0.09409135677924926046703568481130088514796D0)
      PARAMETER (C114=0.1785257973347715166486725999209820685151D0)
      PARAMETER (C115=0.09029786540801834345661648357930909656558D0)
      PARAMETER (C116=0.168315735882869616056559033255302911431D0)
      PARAMETER (C117=0.04486937006121239990457722675738827053758D0)
      PARAMETER (C118=0.2103946698535870200706987915691286392887D0)
      PARAMETER (C119=0.10668957023937637918796145073100368640697D0)
      PARAMETER (C120=0.2293756838200145498406499565786985405995D0)
      PARAMETER (C121=0.1652827715004524673443686414985987229164D0)
      PARAMETER (C122=0.05734392095500363746016248914467463514986D0)
      PARAMETER (C123=0.2061438342970458179897210173367803232506D0)
      PARAMETER (C124=0.06553590966286129995447141616534244017127D0)
      PARAMETER (C125=0.2304758133153235118342646070337224147494D0)
      PARAMETER (C126=0.13926380803358026240325175935135268536395D0)
      PARAMETER (C127=0.2385651315454840938723729945442013829545D0)
      PARAMETER (C128=0.1638397741571532498861785404133561004282D0)
      PARAMETER (C129=0.004764501055633018011170696931004515577894D0)
      PARAMETER (C130=0.2610929189739005798101139280793238944347D0)
      PARAMETER (C131=0.0358357089316044897943677127660920851739D0)
      PARAMETER (C132=0.01782713056910152329557972014130260796931D0)
      PARAMETER (C133=0.1973676950528223941795953744707073569994D0)
      PARAMETER (C134=0.09814013073014567533402673462121237667676D0)
      PARAMETER (C135=0.02521136982901936654914078458909484747415D0)
      PARAMETER (C136=0.14482829338784049941088150392118297407461D0)
      PARAMETER (C137=0.10158468630934461193700448884180434514327D0)
      PARAMETER (C138=0.04366737347227062168452015597187265081384D0)
      PARAMETER (C139=0.10240906836221683967456163835995282469089D0)
      PARAMETER (C140=0.12669836397082431364991375055025294630241D0)
      PARAMETER (C141=0.13166880217999307958157924931993934603115D0)
      PARAMETER (C142=0.09764321548286813795732945522869486291367D0)
      PARAMETER (C143=0.1962802614602913506680534692424247533535D0)
      PARAMETER (C144=0.014293503166899054033512090793013546733683D0)
      PARAMETER (C145=0.1709254378021244678289201177915904423016D0)
      PARAMETER (C146=0.07741397910978242532385262555133610327359D0)
      PARAMETER (C147=0.03300943160428774883052365189043330629385D0)
      PARAMETER (C148=0.1896247477774432756157417125047014442674D0)
      PARAMETER (C149=0.11379365909044610667997832450950090820148D0)
      PARAMETER (C150=0.14991152592791810268952683346048945111274D0)
      PARAMETER (C151=0.08169417525676891267331438070369458194683D0)
      PARAMETER (C152=0.1773777095591540391098792507703541248234D0)
      PARAMETER (C153=0.1713274582033349745023163810977022298431D0)
      PARAMETER (C154=0.10135869117665945091993100044020235704193D0)
      PARAMETER (C155=0.08733474694454124336904031194374530162769D0)
      PARAMETER (C156=0.15126821897411619929484470753456908484491D0)
      PARAMETER (C157=0.09774990997707300453226947369489655936965D0)
      PARAMETER (C158=0.02533967279416486272998275011005058926048D0)
      PARAMETER (C159=0.1691230102307236611787910854701557889493D0)
      PARAMETER (C160=0.06542675382009711688935115641414158445117D0)
      PARAMETER (C161=0.03196124471779480672433711167249291015992D0)
      PARAMETER (C162=0.1060034565616538316654248761356196273065D0)
      PARAMETER (C163=0.11992922074233448215162146676839156089019D0)
      PARAMETER (C164=0.13408494503421884629375947243533521794775D0)
      PARAMETER (C165=0.06392248943558961344867422334498582031983D0)
      PARAMETER (C166=0.12117204378875551047481548675943118176291D0)
      PARAMETER (C167=0.02998230518558362053790536669209789022255D0)
      PARAMETER (C168=0.07518995256510773072529365821957833114334D0)
      PARAMETER (C169=0.10199021561163840488763362306906826659925D0)
      PARAMETER (C170=0.1826737292399166251678861113561525086901D0)
      PARAMETER (C171=0.009577496073872774735268533745381086996268D0)
      PARAMETER (C172=0.10468280611215538702296185026262653512187D0)
      PARAMETER (C173=0.2183368673613531084226007798593632540692D0)
      PARAMETER (C174=0.08601992077982424398839713685150897224529D0)
      PARAMETER (C175=0.08655145984911107891058249054245084411636D0)
      PARAMETER (C176=0.2077235036378665893853979773018820258793D0)
      PARAMETER (C177=0.14762266593163468970699898930017151190222D0)
      PARAMETER (C178=0.10386175181893329469269898865094101293964D0)
      PARAMETER (C179=0.1952864309657362759146589104573897258273D0)
      PARAMETER (C180=0.0226644923581418677528074717931262814665D0)
      PARAMETER (C181=0.2254973469742982149050547806268743852658D0)
      PARAMETER (C182=0.11508946712408130295103272266152988107563D0)
      PARAMETER (C183=0.2358307703785999812089237888203503850832D0)
      PARAMETER (C184=0.14837393116990470670385569807277287562805D0)
      PARAMETER (C185=0.003349134660815399703528695822930909149505D0)
      PARAMETER (C186=0.2686623971685098981932989129819681558348D0)
      PARAMETER (C187=0.02845055206383370111678972019493431637519D0)
      PARAMETER (C188=0.10652530598454139683635066567930671201251D0)
      PARAMETER (C189=0.007104587489308430698352114480886666805154D0)
      PARAMETER (C190=0.1899730028876868537694832496215015501926D0)
      PARAMETER (C191=0.04609516266306470236685292140674448294988D0)
      PARAMETER (C192=0.11909891275269986241726055437740116531125D0)
      PARAMETER (C193=0.01588634057818821138510278381773778257234D0)
      PARAMETER (C194=0.12977626716344185344935463571628007686852D0)
      PARAMETER (C195=0.07791504569418438405126931222614691744872D0)
      PARAMETER (C196=0.11556668879983993534402787050034067300172D0)
      PARAMETER (C197=0.13504547338363840666405030279590328852435D0)
      PARAMETER (C198=0.03042002622521378267054212143165458344031D0)
      PARAMETER (C199=0.08495850966932633015139986393132712776523D0)
      PARAMETER (C200=0.11018851433363497822957909433239914861094D0)
      PARAMETER (C201=0.14153970941131175620817776014328386655746D0)
      PARAMETER (C202=0.07451354221385644069468472870742704673801D0)
      PARAMETER (C203=0.1908521052363872750978983956353611063636D0)
      PARAMETER (C204=0.02679307728652319762822956658344727319604D0)
      PARAMETER (C205=0.1962032737101399204371897459159672455097D0)
      PARAMETER (C206=0.12091484627129322974635631082847084459456D0)
      PARAMETER (C207=0.15288035920378502309925754971208638040008D0)
      PARAMETER (C208=0.1597879589768120952545259985189600680188D0)
      PARAMETER (C209=0.03759394332882173700944515133190833556228D0)
      PARAMETER (C210=0.1699170193386526603027997278626542555305D0)
      PARAMETER (C211=0.11324517011742907267341549271644674120789D0)
      PARAMETER (C212=0.07076985470565587810408888007164193327873D0)
      PARAMETER (C213=0.04501515779454613555468343426530109617478D0)
      PARAMETER (C214=0.06354536231275284554041113527095113028937D0)
      PARAMETER (C215=0.13329387650758653868321569031508415484404D0)
      PARAMETER (C216=0.11687256854127657607690396833922037617308D0)
      PARAMETER (C217=0.02133791404787527583759229014620073728139D0)
      PARAMETER (C218=0.13604249303380869016099095334472176596225D0)
      PARAMETER (C219=0.12319450419485887647641963701580506266066D0)
      PARAMETER (C220=0.09377577050283119169880348304206545618614D0)
      PARAMETER (C221=0.1763312242634143581804722943381810811719D0)
      PARAMETER (C222=0.15647803635108535614234346107213874006355D0)
      PARAMETER (C223=0.11290962813509911256044282177356106217756D0)
      PARAMETER (C224=0.0836984547021396689428469516051695594384D0)
      PARAMETER (C225=0.09946422485031802977692960273241333527216D0)
      PARAMETER (C226=0.1609676485359571884256341886760601607492D0)
      PARAMETER (C227=0.06914274399459705355027938211011672442482D0)
      PARAMETER (C228=0.07288285156561854406341764563912003644508D0)
      PARAMETER (C229=0.10208478235945702492908047518062957026679D0)
      PARAMETER (C230=0.1879697166441086850472257566595416778114D0)
      PARAMETER (C231=0.1875515410056623833976069660841309123723D0)
      PARAMETER (C232=0.1681252256510402655829430902465426927324D0)
      PARAMETER (C233=0.007112638015958425279197430048733579093798D0)
      PARAMETER (C234=0.10307191714852290899486050866839016162528D0)
      PARAMETER (C235=0.06848055384720518368051114222241145772233D0)
      PARAMETER (C236=0.2224087680946349593914389734483289560128D0)
      PARAMETER (C237=0.2344394262570779792470087076051636404653D0)
      PARAMETER (C238=0.14225276031916850558394860097467158187596D0)
      PARAMETER (C239=0.13696110769441036736102228444482291544466D0)
C
      LOGICAL TRANSP
      COMMON /PS3JDT/FACT(0:MXJ+1),POW2(0:MXJ+1),RINT(0:MXJ+1),
     .               RFACT(0:MXJ+1), RSFACT(0:MXJ+1), GHALF(0:MXJ+1)
     ./PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION AP(0:LDA-1,0:*), XP(-L1I:LDX-L1I-1,-L2I:L2I)
      DIMENSION XT(-MAXL:MAXL,-MAXL:MAXL)
C
      IF( N1I+N2I-1.GT.MAXJ ) THEN
          WRITE(NB6,11000) N1I, L1I, N2I, L2I
          STOP 'PSOVP'
      ENDIF
      IF( L1I.GT.MAXL .OR. L2I.GT.MAXL ) THEN
          WRITE(NB6,11100) N1I, L1I, N2I, L2I
          STOP 'PSOVP'
      ENDIF
C
      IF( L1I.GE.L2I ) THEN
          TRANSP = .FALSE.
          N1     = N1I
          L1     = L1I
          ALP1   = ALP1I
          N2     = N2I
          L2     = L2I
          ALP2   = ALP2I
      ELSE
          TRANSP = .TRUE.
          N1     = N2I
          L1     = L2I
          ALP1   = ALP2I
          N2     = N1I
          L2     = L1I
          ALP2   = ALP1I
      ENDIF
C
      GOTO (1000, 
     .      1100,1110,
     .      1200,1210,1220,
     .      1300,1310,1320,1330,
     .      1400,1410,1420,1430,1440) (L1*(L1+1))/2 + L2 + 1
      WRITE(NB6,11100) N1I, L1I, N2I, L2I
      STOP 'PSOVP'
C  L1 = 0  L2 = 0
 1000 CONTINUE
      XT(0,0) = C1*AP(0,0)
      GOTO 2000
C  L1 = 1  L2 = 0
 1100 CONTINUE
      XT(-1,0) = C1*AP(0,0)
      XT(0,0) = C1*AP(1,0)
      XT(1,0) = C1*AP(2,0)
      GOTO 2000
C  L1 = 1  L2 = 1
 1110 CONTINUE
      XT(-1,-1) = C1*AP(0,1) - C2*AP(2,0) - C3*AP(4,0)
      XT(0,-1) = C3*AP(1,0)
      XT(1,-1) = C3*AP(0,0)
      XT(-1,0) = C3*AP(1,0)
      XT(0,0) = C1*AP(0,1) + C4*AP(2,0)
      XT(1,0) = C3*AP(3,0)
      XT(-1,1) = C3*AP(0,0)
      XT(0,1) = C3*AP(3,0)
      XT(1,1) = C1*AP(0,1) - C2*AP(2,0) + C3*AP(4,0)
      GOTO 2000
C  L1 = 2  L2 = 0
 1200 CONTINUE
      XT(-2,0) = C1*AP(0,0)
      XT(-1,0) = C1*AP(1,0)
      XT(0,0) = C1*AP(2,0)
      XT(1,0) = C1*AP(3,0)
      XT(2,0) = C1*AP(4,0)
      GOTO 2000
C  L1 = 2  L2 = 1
 1210 CONTINUE
      XT(-2,-1) = C3*AP(2,1) - C5*AP(4,0) - C6*AP(6,0)
      XT(-1,-1) = C3*AP(1,1) - C7*AP(3,0) - C8*AP(5,0)
      XT(0,-1) = -(C2*AP(0,1)) + C9*AP(2,0)
      XT(1,-1) = C8*AP(1,0)
      XT(2,-1) = C6*AP(0,0) - C3*AP(0,1) + C5*AP(2,0)
      XT(-2,0) = C8*AP(1,0)
      XT(-1,0) = C3*AP(0,1) + C10*AP(2,0)
      XT(0,0) = C4*AP(1,1) + C11*AP(3,0)
      XT(1,0) = C3*AP(2,1) + C10*AP(4,0)
      XT(2,0) = C8*AP(5,0)
      XT(-2,1) = C6*AP(0,0) + C3*AP(0,1) - C5*AP(2,0)
      XT(-1,1) = C8*AP(1,0)
      XT(0,1) = -(C2*AP(2,1)) + C9*AP(4,0)
      XT(1,1) = C3*AP(1,1) - C7*AP(3,0) + C8*AP(5,0)
      XT(2,1) = C3*AP(2,1) - C5*AP(4,0) + C6*AP(6,0)
      GOTO 2000
C  L1 = 2  L2 = 2
 1220 CONTINUE
      XT(-2,-2) = C1*AP(0,2) - C14*AP(2,1) + C12*AP(4,0) - C13*AP(8,0)
      XT(-1,-2) = C17*AP(3,1) - C15*AP(5,0) - C16*AP(7,0)
      XT(0,-2) = -(C14*AP(0,1)) + C17*AP(2,0)
      XT(1,-2) = C16*AP(1,0) + C17*AP(1,1) - C15*AP(3,0)
      XT(2,-2) = C13*AP(0,0)
      XT(-2,-1) = C17*AP(3,1) - C15*AP(5,0) - C16*AP(7,0)
      XT(-1,-1) = C1*AP(0,2) + C19*AP(2,1) - C18*AP(4,0) - 
     .   C17*AP(4,1) - C14*AP(6,0)
      XT(0,-1) = C19*AP(1,1) + C20*AP(3,0)
      XT(1,-1) = C17*AP(0,1) + C14*AP(2,0)
      XT(2,-1) = C16*AP(1,0) - C17*AP(1,1) + C15*AP(3,0)
      XT(-2,0) = -(C14*AP(0,1)) + C17*AP(2,0)
      XT(-1,0) = C19*AP(1,1) + C20*AP(3,0)
      XT(0,0) = C1*AP(0,2) + C14*AP(2,1) + C21*AP(4,0)
      XT(1,0) = C19*AP(3,1) + C20*AP(5,0)
      XT(2,0) = -(C14*AP(4,1)) + C17*AP(6,0)
      XT(-2,1) = C16*AP(1,0) + C17*AP(1,1) - C15*AP(3,0)
      XT(-1,1) = C17*AP(0,1) + C14*AP(2,0)
      XT(0,1) = C19*AP(3,1) + C20*AP(5,0)
      XT(1,1) = C1*AP(0,2) + C19*AP(2,1) - C18*AP(4,0) + C17*AP(4,1) + 
     .   C14*AP(6,0)
      XT(2,1) = C17*AP(3,1) - C15*AP(5,0) + C16*AP(7,0)
      XT(-2,2) = C13*AP(0,0)
      XT(-1,2) = C16*AP(1,0) - C17*AP(1,1) + C15*AP(3,0)
      XT(0,2) = -(C14*AP(4,1)) + C17*AP(6,0)
      XT(1,2) = C17*AP(3,1) - C15*AP(5,0) + C16*AP(7,0)
      XT(2,2) = C1*AP(0,2) - C14*AP(2,1) + C12*AP(4,0) + C13*AP(8,0)
      GOTO 2000
C  L1 = 3  L2 = 0
 1300 CONTINUE
      XT(-3,0) = C1*AP(0,0)
      XT(-2,0) = C1*AP(1,0)
      XT(-1,0) = C1*AP(2,0)
      XT(0,0) = C1*AP(3,0)
      XT(1,0) = C1*AP(4,0)
      XT(2,0) = C1*AP(5,0)
      XT(3,0) = C1*AP(6,0)
      GOTO 2000
C  L1 = 3  L2 = 1
 1310 CONTINUE
      XT(-3,-1) = C6*AP(4,1) - C22*AP(6,0) - C23*AP(8,0)
      XT(-2,-1) = C8*AP(3,1) - C24*AP(5,0) - C25*AP(7,0)
      XT(-1,-1) = C9*AP(2,1) - C26*AP(4,0) + C5*AP(4,1) - C16*AP(6,0)
      XT(0,-1) = -(C7*AP(1,1)) + C27*AP(3,0)
      XT(1,-1) = -(C5*AP(0,1)) + C16*AP(2,0)
      XT(2,-1) = C25*AP(1,0) - C8*AP(1,1) + C24*AP(3,0)
      XT(3,-1) = C23*AP(0,0) - C6*AP(0,1) + C22*AP(2,0)
      XT(-3,0) = C28*AP(1,0)
      XT(-2,0) = C8*AP(0,1) + C29*AP(2,0)
      XT(-1,0) = C10*AP(1,1) + C13*AP(3,0)
      XT(0,0) = C11*AP(2,1) + C30*AP(4,0)
      XT(1,0) = C10*AP(3,1) + C13*AP(5,0)
      XT(2,0) = C8*AP(4,1) + C29*AP(6,0)
      XT(3,0) = C28*AP(7,0)
      XT(-3,1) = C23*AP(0,0) + C6*AP(0,1) - C22*AP(2,0)
      XT(-2,1) = C25*AP(1,0) + C8*AP(1,1) - C24*AP(3,0)
      XT(-1,1) = -(C5*AP(0,1)) + C16*AP(2,0)
      XT(0,1) = -(C7*AP(3,1)) + C27*AP(5,0)
      XT(1,1) = C9*AP(2,1) - C26*AP(4,0) - C5*AP(4,1) + C16*AP(6,0)
      XT(2,1) = C8*AP(3,1) - C24*AP(5,0) + C25*AP(7,0)
      XT(3,1) = C6*AP(4,1) - C22*AP(6,0) + C23*AP(8,0)
      GOTO 2000
C  L1 = 3  L2 = 2
 1320 CONTINUE
      XT(-3,-2) = C6*AP(2,2) - C33*AP(4,1) + C31*AP(6,0) - C32*AP(10,0)
      XT(-2,-2) = C8*AP(1,2) - C36*AP(3,1) + C34*AP(5,0) - C35*AP(9,0)
      XT(-1,-2) = -(C5*AP(2,2)) + C39*AP(4,1) - C37*AP(6,0) + 
     .   C33*AP(6,1) - C38*AP(8,0)
      XT(0,-2) = -(C36*AP(1,1)) + C38*AP(3,0)
      XT(1,-2) = -(C33*AP(0,1)) - C5*AP(0,2) + C38*AP(2,0) + 
     .   C39*AP(2,1) - C37*AP(4,0)
      XT(2,-2) = C35*AP(1,0)
      XT(3,-2) = C32*AP(0,0) - C6*AP(0,2) + C33*AP(2,1) - C31*AP(4,0)
      XT(-3,-1) = C42*AP(5,1) - C40*AP(7,0) - C41*AP(9,0)
      XT(-2,-1) = C8*AP(2,2) + C45*AP(4,1) - C43*AP(6,0) - 
     .   C42*AP(6,1) - C44*AP(8,0)
      XT(-1,-1) = C10*AP(1,2) + C48*AP(3,1) - C46*AP(5,0) - 
     .   C45*AP(5,1) - C47*AP(7,0)
      XT(0,-1) = -(C7*AP(0,2)) + C48*AP(2,1) + C49*AP(4,0)
      XT(1,-1) = C45*AP(1,1) + C47*AP(3,0)
      XT(2,-1) = C42*AP(0,1) - C8*AP(0,2) + C44*AP(2,0) - C45*AP(2,1) + 
     .   C43*AP(4,0)
      XT(3,-1) = C41*AP(1,0) - C42*AP(1,1) + C40*AP(3,0)
      XT(-3,0) = -(C51*AP(0,1)) + C50*AP(2,0)
      XT(-2,0) = C35*AP(3,0)
      XT(-1,0) = C9*AP(0,2) + C2*AP(2,1) + C52*AP(4,0)
      XT(0,0) = C11*AP(1,2) + C54*AP(3,1) + C53*AP(5,0)
      XT(1,0) = C9*AP(2,2) + C2*AP(4,1) + C52*AP(6,0)
      XT(2,0) = C35*AP(7,0)
      XT(3,0) = -(C51*AP(6,1)) + C50*AP(8,0)
      XT(-3,1) = C41*AP(1,0) + C42*AP(1,1) - C40*AP(3,0)
      XT(-2,1) = C42*AP(0,1) + C8*AP(0,2) + C44*AP(2,0) + C45*AP(2,1) - 
     .   C43*AP(4,0)
      XT(-1,1) = C45*AP(1,1) + C47*AP(3,0)
      XT(0,1) = -(C7*AP(2,2)) + C48*AP(4,1) + C49*AP(6,0)
      XT(1,1) = C10*AP(1,2) + C48*AP(3,1) - C46*AP(5,0) + C45*AP(5,1) + 
     .   C47*AP(7,0)
      XT(2,1) = C8*AP(2,2) + C45*AP(4,1) - C43*AP(6,0) + C42*AP(6,1) + 
     .   C44*AP(8,0)
      XT(3,1) = C42*AP(5,1) - C40*AP(7,0) + C41*AP(9,0)
      XT(-3,2) = C32*AP(0,0) + C6*AP(0,2) - C33*AP(2,1) + C31*AP(4,0)
      XT(-2,2) = C35*AP(1,0)
      XT(-1,2) = -(C33*AP(0,1)) + C5*AP(0,2) + C38*AP(2,0) - 
     .   C39*AP(2,1) + C37*AP(4,0)
      XT(0,2) = -(C36*AP(5,1)) + C38*AP(7,0)
      XT(1,2) = -(C5*AP(2,2)) + C39*AP(4,1) - C37*AP(6,0) - 
     .   C33*AP(6,1) + C38*AP(8,0)
      XT(2,2) = C8*AP(1,2) - C36*AP(3,1) + C34*AP(5,0) + C35*AP(9,0)
      XT(3,2) = C6*AP(2,2) - C33*AP(4,1) + C31*AP(6,0) + C32*AP(10,0)
      GOTO 2000
C  L1 = 3  L2 = 3
 1330 CONTINUE
      XT(-3,-3) = C1*AP(0,3) - C51*AP(2,2) + C57*AP(4,1) - 
     .   C55*AP(6,0) - C56*AP(12,0)
      XT(-2,-3) = C42*AP(3,2) - C60*AP(5,1) + C58*AP(7,0) - C59*AP(11,0)
      XT(-1,-3) = -(C33*AP(4,2)) + C63*AP(6,1) - C61*AP(8,0) + 
     .   C64*AP(8,1) - C62*AP(10,0)
      XT(0,-3) = -(C66*AP(1,1)) + C65*AP(3,0)
      XT(1,-3) = -(C64*AP(0,1)) - C33*AP(0,2) + C62*AP(2,0) + 
     .   C63*AP(2,1) - C61*AP(4,0)
      XT(2,-3) = C59*AP(1,0) + C42*AP(1,2) - C60*AP(3,1) + C58*AP(5,0)
      XT(3,-3) = C56*AP(0,0)
      XT(-3,-2) = C42*AP(3,2) - C60*AP(5,1) + C58*AP(7,0) - C59*AP(11,0)
      XT(-2,-2) = C1*AP(0,3) - C69*AP(4,1) + C67*AP(6,0) - 
     .   C70*AP(8,1) - C68*AP(10,0)
      XT(-1,-2) = C45*AP(3,2) + C73*AP(5,1) - C71*AP(7,0) - 
     .   C74*AP(7,1) - C72*AP(9,0)
      XT(0,-2) = -(C36*AP(0,2)) - C76*AP(2,1) + C75*AP(4,0)
      XT(1,-2) = C74*AP(1,1) + C45*AP(1,2) + C72*AP(3,0) + 
     .   C73*AP(3,1) - C71*AP(5,0)
      XT(2,-2) = C70*AP(0,1) + C68*AP(2,0)
      XT(3,-2) = C59*AP(1,0) - C42*AP(1,2) + C60*AP(3,1) - C58*AP(5,0)
      XT(-3,-1) = -(C33*AP(4,2)) + C63*AP(6,1) - C61*AP(8,0) + 
     .   C64*AP(8,1) - C62*AP(10,0)
      XT(-2,-1) = C45*AP(3,2) + C73*AP(5,1) - C71*AP(7,0) - 
     .   C74*AP(7,1) - C72*AP(9,0)
      XT(-1,-1) = C1*AP(0,3) + C2*AP(2,2) + C79*AP(4,1) - C39*AP(4,2) - 
     .   C77*AP(6,0) - C80*AP(6,1) - C78*AP(8,0)
      XT(0,-1) = C48*AP(1,2) + C60*AP(3,1) + C81*AP(5,0)
      XT(1,-1) = C39*AP(0,2) + C80*AP(2,1) + C78*AP(4,0)
      XT(2,-1) = C74*AP(1,1) - C45*AP(1,2) + C72*AP(3,0) - 
     .   C73*AP(3,1) + C71*AP(5,0)
      XT(3,-1) = -(C64*AP(0,1)) + C33*AP(0,2) + C62*AP(2,0) - 
     .   C63*AP(2,1) + C61*AP(4,0)
      XT(-3,0) = -(C66*AP(1,1)) + C65*AP(3,0)
      XT(-2,0) = -(C36*AP(0,2)) - C76*AP(2,1) + C75*AP(4,0)
      XT(-1,0) = C48*AP(1,2) + C60*AP(3,1) + C81*AP(5,0)
      XT(0,0) = C1*AP(0,3) + C54*AP(2,2) + C83*AP(4,1) + C82*AP(6,0)
      XT(1,0) = C48*AP(3,2) + C60*AP(5,1) + C81*AP(7,0)
      XT(2,0) = -(C36*AP(4,2)) - C76*AP(6,1) + C75*AP(8,0)
      XT(3,0) = -(C66*AP(7,1)) + C65*AP(9,0)
      XT(-3,1) = -(C64*AP(0,1)) - C33*AP(0,2) + C62*AP(2,0) + 
     .   C63*AP(2,1) - C61*AP(4,0)
      XT(-2,1) = C74*AP(1,1) + C45*AP(1,2) + C72*AP(3,0) + 
     .   C73*AP(3,1) - C71*AP(5,0)
      XT(-1,1) = C39*AP(0,2) + C80*AP(2,1) + C78*AP(4,0)
      XT(0,1) = C48*AP(3,2) + C60*AP(5,1) + C81*AP(7,0)
      XT(1,1) = C1*AP(0,3) + C2*AP(2,2) + C79*AP(4,1) + C39*AP(4,2) - 
     .   C77*AP(6,0) + C80*AP(6,1) + C78*AP(8,0)
      XT(2,1) = C45*AP(3,2) + C73*AP(5,1) - C71*AP(7,0) + C74*AP(7,1) + 
     .   C72*AP(9,0)
      XT(3,1) = -(C33*AP(4,2)) + C63*AP(6,1) - C61*AP(8,0) - 
     .   C64*AP(8,1) + C62*AP(10,0)
      XT(-3,2) = C59*AP(1,0) + C42*AP(1,2) - C60*AP(3,1) + C58*AP(5,0)
      XT(-2,2) = C70*AP(0,1) + C68*AP(2,0)
      XT(-1,2) = C74*AP(1,1) - C45*AP(1,2) + C72*AP(3,0) - 
     .   C73*AP(3,1) + C71*AP(5,0)
      XT(0,2) = -(C36*AP(4,2)) - C76*AP(6,1) + C75*AP(8,0)
      XT(1,2) = C45*AP(3,2) + C73*AP(5,1) - C71*AP(7,0) + C74*AP(7,1) + 
     .   C72*AP(9,0)
      XT(2,2) = C1*AP(0,3) - C69*AP(4,1) + C67*AP(6,0) + C70*AP(8,1) + 
     .   C68*AP(10,0)
      XT(3,2) = C42*AP(3,2) - C60*AP(5,1) + C58*AP(7,0) + C59*AP(11,0)
      XT(-3,3) = C56*AP(0,0)
      XT(-2,3) = C59*AP(1,0) - C42*AP(1,2) + C60*AP(3,1) - C58*AP(5,0)
      XT(-1,3) = -(C64*AP(0,1)) + C33*AP(0,2) + C62*AP(2,0) - 
     .   C63*AP(2,1) + C61*AP(4,0)
      XT(0,3) = -(C66*AP(7,1)) + C65*AP(9,0)
      XT(1,3) = -(C33*AP(4,2)) + C63*AP(6,1) - C61*AP(8,0) - 
     .   C64*AP(8,1) + C62*AP(10,0)
      XT(2,3) = C42*AP(3,2) - C60*AP(5,1) + C58*AP(7,0) + C59*AP(11,0)
      XT(3,3) = C1*AP(0,3) - C51*AP(2,2) + C57*AP(4,1) - C55*AP(6,0) + 
     .   C56*AP(12,0)
      GOTO 2000
C  L1 = 4  L2 = 0
 1400 CONTINUE
      XT(-4,0) = C1*AP(0,0)
      XT(-3,0) = C1*AP(1,0)
      XT(-2,0) = C1*AP(2,0)
      XT(-1,0) = C1*AP(3,0)
      XT(0,0) = C1*AP(4,0)
      XT(1,0) = C1*AP(5,0)
      XT(2,0) = C1*AP(6,0)
      XT(3,0) = C1*AP(7,0)
      XT(4,0) = C1*AP(8,0)
      GOTO 2000
C  L1 = 4  L2 = 1
 1410 CONTINUE
      XT(-4,-1) = C23*AP(6,1) - C84*AP(8,0) - C85*AP(10,0)
      XT(-3,-1) = C25*AP(5,1) - C86*AP(7,0) - C87*AP(9,0)
      XT(-2,-1) = C16*AP(4,1) - C88*AP(6,0) + C22*AP(6,1) - C89*AP(8,0)
      XT(-1,-1) = C27*AP(3,1) - C41*AP(5,0) + C24*AP(5,1) - C90*AP(7,0)
      XT(0,-1) = -(C26*AP(2,1)) + C35*AP(4,0)
      XT(1,-1) = -(C24*AP(1,1)) + C90*AP(3,0)
      XT(2,-1) = -(C22*AP(0,1)) + C89*AP(2,0) - C16*AP(2,1) + 
     .   C88*AP(4,0)
      XT(3,-1) = C87*AP(1,0) - C25*AP(1,1) + C86*AP(3,0)
      XT(4,-1) = C85*AP(0,0) - C23*AP(0,1) + C84*AP(2,0)
      XT(-4,0) = C91*AP(1,0)
      XT(-3,0) = C28*AP(0,1) + C92*AP(2,0)
      XT(-2,0) = C29*AP(1,1) + C93*AP(3,0)
      XT(-1,0) = C13*AP(2,1) + C94*AP(4,0)
      XT(0,0) = C30*AP(3,1) + C32*AP(5,0)
      XT(1,0) = C13*AP(4,1) + C94*AP(6,0)
      XT(2,0) = C29*AP(5,1) + C93*AP(7,0)
      XT(3,0) = C28*AP(6,1) + C92*AP(8,0)
      XT(4,0) = C91*AP(9,0)
      XT(-4,1) = C85*AP(0,0) + C23*AP(0,1) - C84*AP(2,0)
      XT(-3,1) = C87*AP(1,0) + C25*AP(1,1) - C86*AP(3,0)
      XT(-2,1) = -(C22*AP(0,1)) + C89*AP(2,0) + C16*AP(2,1) - 
     .   C88*AP(4,0)
      XT(-1,1) = -(C24*AP(1,1)) + C90*AP(3,0)
      XT(0,1) = -(C26*AP(4,1)) + C35*AP(6,0)
      XT(1,1) = C27*AP(3,1) - C41*AP(5,0) - C24*AP(5,1) + C90*AP(7,0)
      XT(2,1) = C16*AP(4,1) - C88*AP(6,0) - C22*AP(6,1) + C89*AP(8,0)
      XT(3,1) = C25*AP(5,1) - C86*AP(7,0) + C87*AP(9,0)
      XT(4,1) = C23*AP(6,1) - C84*AP(8,0) + C85*AP(10,0)
      GOTO 2000
C  L1 = 4  L2 = 2
 1420 CONTINUE
      XT(-4,-2) = C13*AP(4,2) - C97*AP(6,1) + C95*AP(8,0) - C96*AP(12,0)
      XT(-3,-2) = C16*AP(3,2) - C100*AP(5,1) + C98*AP(7,0) - 
     .   C99*AP(11,0)
      XT(-2,-2) = C17*AP(2,2) - C102*AP(4,1) + C101*AP(6,0) + 
     .   C97*AP(8,1) - C72*AP(10,0)
      XT(-1,-2) = -(C15*AP(3,2)) + C105*AP(5,1) - C103*AP(7,0) + 
     .   C100*AP(7,1) - C104*AP(9,0)
      XT(0,-2) = C12*AP(0,2) - C102*AP(2,1) + C106*AP(4,0)
      XT(1,-2) = -(C100*AP(1,1)) - C15*AP(1,2) + C104*AP(3,0) + 
     .   C105*AP(3,1) - C103*AP(5,0)
      XT(2,-2) = -(C97*AP(0,1)) + C72*AP(2,0)
      XT(3,-2) = C99*AP(1,0) - C16*AP(1,2) + C100*AP(3,1) - C98*AP(5,0)
      XT(4,-2) = C96*AP(0,0) - C13*AP(0,2) + C97*AP(2,1) - C95*AP(4,0)
      XT(-4,-1) = C109*AP(7,1) - C107*AP(9,0) - C108*AP(11,0)
      XT(-3,-1) = C16*AP(4,2) + C112*AP(6,1) - C110*AP(8,0) - 
     .   C109*AP(8,1) - C111*AP(10,0)
      XT(-2,-1) = C14*AP(3,2) + C115*AP(5,1) - C113*AP(7,0) - 
     .   C112*AP(7,1) - C114*AP(9,0)
      XT(-1,-1) = C20*AP(2,2) + C117*AP(4,1) + C15*AP(4,2) - 
     .   C111*AP(6,0) - C115*AP(6,1) - C116*AP(8,0)
      XT(0,-1) = -(C18*AP(1,2)) + C117*AP(3,1) + C118*AP(5,0)
      XT(1,-1) = -(C15*AP(0,2)) + C115*AP(2,1) + C116*AP(4,0)
      XT(2,-1) = C112*AP(1,1) - C14*AP(1,2) + C114*AP(3,0) - 
     .   C115*AP(3,1) + C113*AP(5,0)
      XT(3,-1) = C109*AP(0,1) - C16*AP(0,2) + C111*AP(2,0) - 
     .   C112*AP(2,1) + C110*AP(4,0)
      XT(4,-1) = C108*AP(1,0) - C109*AP(1,1) + C107*AP(3,0)
      XT(-4,0) = -(C120*AP(0,1)) + C119*AP(2,0)
      XT(-3,0) = -(C122*AP(1,1)) + C121*AP(3,0)
      XT(-2,0) = C17*AP(0,2) + C124*AP(2,1) + C123*AP(4,0)
      XT(-1,0) = C20*AP(1,2) + C126*AP(3,1) + C125*AP(5,0)
      XT(0,0) = C21*AP(2,2) + C128*AP(4,1) + C127*AP(6,0)
      XT(1,0) = C20*AP(3,2) + C126*AP(5,1) + C125*AP(7,0)
      XT(2,0) = C17*AP(4,2) + C124*AP(6,1) + C123*AP(8,0)
      XT(3,0) = -(C122*AP(7,1)) + C121*AP(9,0)
      XT(4,0) = -(C120*AP(8,1)) + C119*AP(10,0)
      XT(-4,1) = C108*AP(1,0) + C109*AP(1,1) - C107*AP(3,0)
      XT(-3,1) = C109*AP(0,1) + C16*AP(0,2) + C111*AP(2,0) + 
     .   C112*AP(2,1) - C110*AP(4,0)
      XT(-2,1) = C112*AP(1,1) + C14*AP(1,2) + C114*AP(3,0) + 
     .   C115*AP(3,1) - C113*AP(5,0)
      XT(-1,1) = -(C15*AP(0,2)) + C115*AP(2,1) + C116*AP(4,0)
      XT(0,1) = -(C18*AP(3,2)) + C117*AP(5,1) + C118*AP(7,0)
      XT(1,1) = C20*AP(2,2) + C117*AP(4,1) - C15*AP(4,2) - 
     .   C111*AP(6,0) + C115*AP(6,1) + C116*AP(8,0)
      XT(2,1) = C14*AP(3,2) + C115*AP(5,1) - C113*AP(7,0) + 
     .   C112*AP(7,1) + C114*AP(9,0)
      XT(3,1) = C16*AP(4,2) + C112*AP(6,1) - C110*AP(8,0) + 
     .   C109*AP(8,1) + C111*AP(10,0)
      XT(4,1) = C109*AP(7,1) - C107*AP(9,0) + C108*AP(11,0)
      XT(-4,2) = C96*AP(0,0) + C13*AP(0,2) - C97*AP(2,1) + C95*AP(4,0)
      XT(-3,2) = C99*AP(1,0) + C16*AP(1,2) - C100*AP(3,1) + C98*AP(5,0)
      XT(-2,2) = -(C97*AP(0,1)) + C72*AP(2,0)
      XT(-1,2) = -(C100*AP(1,1)) + C15*AP(1,2) + C104*AP(3,0) - 
     .   C105*AP(3,1) + C103*AP(5,0)
      XT(0,2) = C12*AP(4,2) - C102*AP(6,1) + C106*AP(8,0)
      XT(1,2) = -(C15*AP(3,2)) + C105*AP(5,1) - C103*AP(7,0) - 
     .   C100*AP(7,1) + C104*AP(9,0)
      XT(2,2) = C17*AP(2,2) - C102*AP(4,1) + C101*AP(6,0) - 
     .   C97*AP(8,1) + C72*AP(10,0)
      XT(3,2) = C16*AP(3,2) - C100*AP(5,1) + C98*AP(7,0) + C99*AP(11,0)
      XT(4,2) = C13*AP(4,2) - C97*AP(6,1) + C95*AP(8,0) + C96*AP(12,0)
      GOTO 2000
C  L1 = 4  L2 = 3
 1430 CONTINUE
      XT(-4,-3) = C23*AP(2,3) - C64*AP(4,2) + C131*AP(6,1) - 
     .   C129*AP(8,0) - C130*AP(14,0)
      XT(-3,-3) = C28*AP(1,3) - C66*AP(3,2) + C134*AP(5,1) - 
     .   C132*AP(7,0) - C133*AP(13,0)
      XT(-2,-3) = -(C22*AP(2,3)) + C63*AP(4,2) - C137*AP(6,1) + 
     .   C135*AP(8,0) + C134*AP(10,1) - C136*AP(12,0)
      XT(-1,-3) = -(C60*AP(5,2)) + C140*AP(7,1) - C138*AP(9,0) + 
     .   C141*AP(9,1) - C139*AP(11,0)
      XT(0,-3) = C57*AP(0,2) - C143*AP(2,1) + C142*AP(4,0)
      XT(1,-3) = -(C141*AP(1,1)) - C60*AP(1,2) + C139*AP(3,0) + 
     .   C140*AP(3,1) - C138*AP(5,0)
      XT(2,-3) = -(C134*AP(0,1)) - C22*AP(0,3) + C136*AP(2,0) + 
     .   C63*AP(2,2) - C137*AP(4,1) + C135*AP(6,0)
      XT(3,-3) = C133*AP(1,0)
      XT(4,-3) = C130*AP(0,0) - C23*AP(0,3) + C64*AP(2,2) - 
     .   C131*AP(4,1) + C129*AP(6,0)
      XT(-4,-2) = C70*AP(5,2) - C146*AP(7,1) + C144*AP(9,0) - 
     .   C145*AP(13,0)
      XT(-3,-2) = C25*AP(2,3) + C74*AP(4,2) - C149*AP(6,1) + 
     .   C147*AP(8,0) - C150*AP(10,1) - C148*AP(12,0)
      XT(-2,-2) = C29*AP(1,3) - C76*AP(3,2) - C153*AP(5,1) + 
     .   C151*AP(7,0) - C154*AP(9,1) - C152*AP(11,0)
      XT(-1,-2) = -(C24*AP(2,3)) + C73*AP(4,2) + C157*AP(6,1) + 
     .   C60*AP(6,2) - C155*AP(8,0) - C158*AP(8,1) - C156*AP(10,0)
      XT(0,-2) = -(C69*AP(1,2)) - C160*AP(3,1) + C159*AP(5,0)
      XT(1,-2) = -(C60*AP(0,2)) - C24*AP(0,3) + C158*AP(2,1) + 
     .   C73*AP(2,2) + C156*AP(4,0) + C157*AP(4,1) - C155*AP(6,0)
      XT(2,-2) = C154*AP(1,1) + C152*AP(3,0)
      XT(3,-2) = C150*AP(0,1) - C25*AP(0,3) + C148*AP(2,0) - 
     .   C74*AP(2,2) + C149*AP(4,1) - C147*AP(6,0)
      XT(4,-2) = C145*AP(1,0) - C70*AP(1,2) + C146*AP(3,1) - 
     .   C144*AP(5,0)
      XT(-4,-1) = -(C64*AP(6,2)) + C163*AP(8,1) - C161*AP(10,0) + 
     .   C164*AP(10,1) - C162*AP(12,0)
      XT(-3,-1) = C74*AP(5,2) + C166*AP(7,1) - C165*AP(9,0) - 
     .   C167*AP(9,1) - C150*AP(11,0)
      XT(-2,-1) = C16*AP(2,3) + C80*AP(4,2) + C168*AP(6,1) - 
     .   C63*AP(6,2) - C142*AP(8,0) - C169*AP(8,1) - C159*AP(10,0)
      XT(-1,-1) = C13*AP(1,3) + C60*AP(3,2) + C171*AP(5,1) - 
     .   C73*AP(5,2) - C170*AP(7,0) - C172*AP(7,1) - C159*AP(9,0)
      XT(0,-1) = -(C26*AP(0,3)) + C79*AP(2,2) + C174*AP(4,1) + 
     .   C173*AP(6,0)
      XT(1,-1) = C73*AP(1,2) + C172*AP(3,1) + C159*AP(5,0)
      XT(2,-1) = C63*AP(0,2) - C16*AP(0,3) + C169*AP(2,1) - 
     .   C80*AP(2,2) + C159*AP(4,0) - C168*AP(4,1) + C142*AP(6,0)
      XT(3,-1) = C167*AP(1,1) - C74*AP(1,2) + C150*AP(3,0) - 
     .   C166*AP(3,1) + C165*AP(5,0)
      XT(4,-1) = -(C164*AP(0,1)) + C64*AP(0,2) + C162*AP(2,0) - 
     .   C163*AP(2,1) + C161*AP(4,0)
      XT(-4,0) = -(C176*AP(1,1)) + C175*AP(3,0)
      XT(-3,0) = -(C66*AP(0,2)) - C178*AP(2,1) + C177*AP(4,0)
      XT(-2,0) = -(C76*AP(1,2)) + C180*AP(3,1) + C179*AP(5,0)
      XT(-1,0) = C27*AP(0,3) + C60*AP(2,2) + C182*AP(4,1) + C181*AP(6,0)
      XT(0,0) = C30*AP(1,3) + C83*AP(3,2) + C184*AP(5,1) + C183*AP(7,0)
      XT(1,0) = C27*AP(2,3) + C60*AP(4,2) + C182*AP(6,1) + C181*AP(8,0)
      XT(2,0) = -(C76*AP(5,2)) + C180*AP(7,1) + C179*AP(9,0)
      XT(3,0) = -(C66*AP(6,2)) - C178*AP(8,1) + C177*AP(10,0)
      XT(4,0) = -(C176*AP(9,1)) + C175*AP(11,0)
      XT(-4,1) = -(C164*AP(0,1)) - C64*AP(0,2) + C162*AP(2,0) + 
     .   C163*AP(2,1) - C161*AP(4,0)
      XT(-3,1) = C167*AP(1,1) + C74*AP(1,2) + C150*AP(3,0) + 
     .   C166*AP(3,1) - C165*AP(5,0)
      XT(-2,1) = C63*AP(0,2) + C16*AP(0,3) + C169*AP(2,1) + 
     .   C80*AP(2,2) + C159*AP(4,0) + C168*AP(4,1) - C142*AP(6,0)
      XT(-1,1) = C73*AP(1,2) + C172*AP(3,1) + C159*AP(5,0)
      XT(0,1) = -(C26*AP(2,3)) + C79*AP(4,2) + C174*AP(6,1) + 
     .   C173*AP(8,0)
      XT(1,1) = C13*AP(1,3) + C60*AP(3,2) + C171*AP(5,1) + 
     .   C73*AP(5,2) - C170*AP(7,0) + C172*AP(7,1) + C159*AP(9,0)
      XT(2,1) = C16*AP(2,3) + C80*AP(4,2) + C168*AP(6,1) + 
     .   C63*AP(6,2) - C142*AP(8,0) + C169*AP(8,1) + C159*AP(10,0)
      XT(3,1) = C74*AP(5,2) + C166*AP(7,1) - C165*AP(9,0) + 
     .   C167*AP(9,1) + C150*AP(11,0)
      XT(4,1) = -(C64*AP(6,2)) + C163*AP(8,1) - C161*AP(10,0) - 
     .   C164*AP(10,1) + C162*AP(12,0)
      XT(-4,2) = C145*AP(1,0) + C70*AP(1,2) - C146*AP(3,1) + 
     .   C144*AP(5,0)
      XT(-3,2) = C150*AP(0,1) + C25*AP(0,3) + C148*AP(2,0) + 
     .   C74*AP(2,2) - C149*AP(4,1) + C147*AP(6,0)
      XT(-2,2) = C154*AP(1,1) + C152*AP(3,0)
      XT(-1,2) = -(C60*AP(0,2)) + C24*AP(0,3) + C158*AP(2,1) - 
     .   C73*AP(2,2) + C156*AP(4,0) - C157*AP(4,1) + C155*AP(6,0)
      XT(0,2) = -(C69*AP(5,2)) - C160*AP(7,1) + C159*AP(9,0)
      XT(1,2) = -(C24*AP(2,3)) + C73*AP(4,2) + C157*AP(6,1) - 
     .   C60*AP(6,2) - C155*AP(8,0) + C158*AP(8,1) + C156*AP(10,0)
      XT(2,2) = C29*AP(1,3) - C76*AP(3,2) - C153*AP(5,1) + 
     .   C151*AP(7,0) + C154*AP(9,1) + C152*AP(11,0)
      XT(3,2) = C25*AP(2,3) + C74*AP(4,2) - C149*AP(6,1) + 
     .   C147*AP(8,0) + C150*AP(10,1) + C148*AP(12,0)
      XT(4,2) = C70*AP(5,2) - C146*AP(7,1) + C144*AP(9,0) + 
     .   C145*AP(13,0)
      XT(-4,3) = C130*AP(0,0) + C23*AP(0,3) - C64*AP(2,2) + 
     .   C131*AP(4,1) - C129*AP(6,0)
      XT(-3,3) = C133*AP(1,0)
      XT(-2,3) = -(C134*AP(0,1)) + C22*AP(0,3) + C136*AP(2,0) - 
     .   C63*AP(2,2) + C137*AP(4,1) - C135*AP(6,0)
      XT(-1,3) = -(C141*AP(1,1)) + C60*AP(1,2) + C139*AP(3,0) - 
     .   C140*AP(3,1) + C138*AP(5,0)
      XT(0,3) = C57*AP(6,2) - C143*AP(8,1) + C142*AP(10,0)
      XT(1,3) = -(C60*AP(5,2)) + C140*AP(7,1) - C138*AP(9,0) - 
     .   C141*AP(9,1) + C139*AP(11,0)
      XT(2,3) = -(C22*AP(2,3)) + C63*AP(4,2) - C137*AP(6,1) + 
     .   C135*AP(8,0) - C134*AP(10,1) + C136*AP(12,0)
      XT(3,3) = C28*AP(1,3) - C66*AP(3,2) + C134*AP(5,1) - 
     .   C132*AP(7,0) + C133*AP(13,0)
      XT(4,3) = C23*AP(2,3) - C64*AP(4,2) + C131*AP(6,1) - 
     .   C129*AP(8,0) + C130*AP(14,0)
      GOTO 2000
C  L1 = 4  L2 = 4
 1440 CONTINUE
      XT(-4,-4) = C1*AP(0,4) - C120*AP(2,3) + C188*AP(4,2) - 
     .   C187*AP(6,1) + C185*AP(8,0) - C186*AP(16,0)
      XT(-3,-4) = C109*AP(3,3) - C192*AP(5,2) + C191*AP(7,1) - 
     .   C189*AP(9,0) - C190*AP(15,0)
      XT(-2,-4) = -(C97*AP(4,3)) + C197*AP(6,2) - C195*AP(8,1) + 
     .   C193*AP(10,0) + C196*AP(12,1) - C194*AP(14,0)
      XT(-1,-4) = -(C192*AP(7,2)) + C200*AP(9,1) - C198*AP(11,0) + 
     .   C201*AP(11,1) - C199*AP(13,0)
      XT(0,-4) = C188*AP(0,2) - C203*AP(2,1) + C202*AP(4,0)
      XT(1,-4) = -(C201*AP(1,1)) - C192*AP(1,2) + C199*AP(3,0) + 
     .   C200*AP(3,1) - C198*AP(5,0)
      XT(2,-4) = -(C196*AP(0,1)) - C97*AP(0,3) + C194*AP(2,0) + 
     .   C197*AP(2,2) - C195*AP(4,1) + C193*AP(6,0)
      XT(3,-4) = C190*AP(1,0) + C109*AP(1,3) - C192*AP(3,2) + 
     .   C191*AP(5,1) - C189*AP(7,0)
      XT(4,-4) = C186*AP(0,0)
      XT(-4,-3) = C109*AP(3,3) - C192*AP(5,2) + C191*AP(7,1) - 
     .   C189*AP(9,0) - C190*AP(15,0)
      XT(-3,-3) = C1*AP(0,4) - C122*AP(2,3) - C208*AP(4,2) + 
     .   C206*AP(6,1) - C204*AP(8,0) - C207*AP(12,1) - C205*AP(14,0)
      XT(-2,-3) = C112*AP(3,3) + C213*AP(5,2) - C211*AP(7,1) + 
     .   C209*AP(9,0) - C212*AP(11,1) - C210*AP(13,0)
      XT(-1,-3) = -(C100*AP(4,3)) + C213*AP(6,2) + C216*AP(8,1) + 
     .   C192*AP(8,2) - C214*AP(10,0) + C217*AP(10,1) - C215*AP(12,0)
      XT(0,-3) = -(C208*AP(1,2)) - C219*AP(3,1) + C218*AP(5,0)
      XT(1,-3) = -(C192*AP(0,2)) - C100*AP(0,3) - C217*AP(2,1) + 
     .   C213*AP(2,2) + C215*AP(4,0) + C216*AP(4,1) - C214*AP(6,0)
      XT(2,-3) = C212*AP(1,1) + C112*AP(1,3) + C210*AP(3,0) + 
     .   C213*AP(3,2) - C211*AP(5,1) + C209*AP(7,0)
      XT(3,-3) = C207*AP(0,1) + C205*AP(2,0)
      XT(4,-3) = C190*AP(1,0) - C109*AP(1,3) + C192*AP(3,2) - 
     .   C191*AP(5,1) + C189*AP(7,0)
      XT(-4,-2) = -(C97*AP(4,3)) + C197*AP(6,2) - C195*AP(8,1) + 
     .   C193*AP(10,0) + C196*AP(12,1) - C194*AP(14,0)
      XT(-3,-2) = C112*AP(3,3) + C213*AP(5,2) - C211*AP(7,1) + 
     .   C209*AP(9,0) - C212*AP(11,1) - C210*AP(13,0)
      XT(-2,-2) = C1*AP(0,4) + C124*AP(2,3) - C224*AP(4,2) - 
     .   C222*AP(6,1) + C220*AP(8,0) - C197*AP(8,2) - C223*AP(10,1) - 
     .   C221*AP(12,0)
      XT(-1,-2) = C115*AP(3,3) + C229*AP(5,2) + C227*AP(7,1) - 
     .   C213*AP(7,2) - C225*AP(9,0) - C228*AP(9,1) - C226*AP(11,0)
      XT(0,-2) = -(C102*AP(0,3)) - C224*AP(2,2) + C230*AP(6,0)
      XT(1,-2) = C213*AP(1,2) + C115*AP(1,3) + C228*AP(3,1) + 
     .   C229*AP(3,2) + C226*AP(5,0) + C227*AP(5,1) - C225*AP(7,0)
      XT(2,-2) = C197*AP(0,2) + C223*AP(2,1) + C221*AP(4,0)
      XT(3,-2) = C212*AP(1,1) - C112*AP(1,3) + C210*AP(3,0) - 
     .   C213*AP(3,2) + C211*AP(5,1) - C209*AP(7,0)
      XT(4,-2) = -(C196*AP(0,1)) + C97*AP(0,3) + C194*AP(2,0) - 
     .   C197*AP(2,2) + C195*AP(4,1) - C193*AP(6,0)
      XT(-4,-1) = -(C192*AP(7,2)) + C200*AP(9,1) - C198*AP(11,0) + 
     .   C201*AP(11,1) - C199*AP(13,0)
      XT(-3,-1) = -(C100*AP(4,3)) + C213*AP(6,2) + C216*AP(8,1) + 
     .   C192*AP(8,2) - C214*AP(10,0) + C217*AP(10,1) - C215*AP(12,0)
      XT(-2,-1) = C115*AP(3,3) + C229*AP(5,2) + C227*AP(7,1) - 
     .   C213*AP(7,2) - C225*AP(9,0) - C228*AP(9,1) - C226*AP(11,0)
      XT(-1,-1) = C1*AP(0,4) + C126*AP(2,3) + C235*AP(4,2) - 
     .   C105*AP(4,3) - C233*AP(6,1) - C229*AP(6,2) - C231*AP(8,0) - 
     .   C234*AP(8,1) - C232*AP(10,0)
      XT(0,-1) = C117*AP(1,3) + C235*AP(3,2) + C234*AP(5,1) + 
     .   C236*AP(7,0)
      XT(1,-1) = C105*AP(0,3) + C229*AP(2,2) + C234*AP(4,1) + 
     .   C232*AP(6,0)
      XT(2,-1) = C213*AP(1,2) - C115*AP(1,3) + C228*AP(3,1) - 
     .   C229*AP(3,2) + C226*AP(5,0) - C227*AP(5,1) + C225*AP(7,0)
      XT(3,-1) = -(C192*AP(0,2)) + C100*AP(0,3) - C217*AP(2,1) - 
     .   C213*AP(2,2) + C215*AP(4,0) - C216*AP(4,1) + C214*AP(6,0)
      XT(4,-1) = -(C201*AP(1,1)) + C192*AP(1,2) + C199*AP(3,0) - 
     .   C200*AP(3,1) + C198*AP(5,0)
      XT(-4,0) = C188*AP(0,2) - C203*AP(2,1) + C202*AP(4,0)
      XT(-3,0) = -(C208*AP(1,2)) - C219*AP(3,1) + C218*AP(5,0)
      XT(-2,0) = -(C102*AP(0,3)) - C224*AP(2,2) + C230*AP(6,0)
      XT(-1,0) = C117*AP(1,3) + C235*AP(3,2) + C234*AP(5,1) + 
     .   C236*AP(7,0)
      XT(0,0) = C1*AP(0,4) + C128*AP(2,3) + C239*AP(4,2) + 
     .   C238*AP(6,1) + C237*AP(8,0)
      XT(1,0) = C117*AP(3,3) + C235*AP(5,2) + C234*AP(7,1) + 
     .   C236*AP(9,0)
      XT(2,0) = -(C102*AP(4,3)) - C224*AP(6,2) + C230*AP(10,0)
      XT(3,0) = -(C208*AP(7,2)) - C219*AP(9,1) + C218*AP(11,0)
      XT(4,0) = C188*AP(8,2) - C203*AP(10,1) + C202*AP(12,0)
      XT(-4,1) = -(C201*AP(1,1)) - C192*AP(1,2) + C199*AP(3,0) + 
     .   C200*AP(3,1) - C198*AP(5,0)
      XT(-3,1) = -(C192*AP(0,2)) - C100*AP(0,3) - C217*AP(2,1) + 
     .   C213*AP(2,2) + C215*AP(4,0) + C216*AP(4,1) - C214*AP(6,0)
      XT(-2,1) = C213*AP(1,2) + C115*AP(1,3) + C228*AP(3,1) + 
     .   C229*AP(3,2) + C226*AP(5,0) + C227*AP(5,1) - C225*AP(7,0)
      XT(-1,1) = C105*AP(0,3) + C229*AP(2,2) + C234*AP(4,1) + 
     .   C232*AP(6,0)
      XT(0,1) = C117*AP(3,3) + C235*AP(5,2) + C234*AP(7,1) + 
     .   C236*AP(9,0)
      XT(1,1) = C1*AP(0,4) + C126*AP(2,3) + C235*AP(4,2) + 
     .   C105*AP(4,3) - C233*AP(6,1) + C229*AP(6,2) - C231*AP(8,0) + 
     .   C234*AP(8,1) + C232*AP(10,0)
      XT(2,1) = C115*AP(3,3) + C229*AP(5,2) + C227*AP(7,1) + 
     .   C213*AP(7,2) - C225*AP(9,0) + C228*AP(9,1) + C226*AP(11,0)
      XT(3,1) = -(C100*AP(4,3)) + C213*AP(6,2) + C216*AP(8,1) - 
     .   C192*AP(8,2) - C214*AP(10,0) - C217*AP(10,1) + C215*AP(12,0)
      XT(4,1) = -(C192*AP(7,2)) + C200*AP(9,1) - C198*AP(11,0) - 
     .   C201*AP(11,1) + C199*AP(13,0)
      XT(-4,2) = -(C196*AP(0,1)) - C97*AP(0,3) + C194*AP(2,0) + 
     .   C197*AP(2,2) - C195*AP(4,1) + C193*AP(6,0)
      XT(-3,2) = C212*AP(1,1) + C112*AP(1,3) + C210*AP(3,0) + 
     .   C213*AP(3,2) - C211*AP(5,1) + C209*AP(7,0)
      XT(-2,2) = C197*AP(0,2) + C223*AP(2,1) + C221*AP(4,0)
      XT(-1,2) = C213*AP(1,2) - C115*AP(1,3) + C228*AP(3,1) - 
     .   C229*AP(3,2) + C226*AP(5,0) - C227*AP(5,1) + C225*AP(7,0)
      XT(0,2) = -(C102*AP(4,3)) - C224*AP(6,2) + C230*AP(10,0)
      XT(1,2) = C115*AP(3,3) + C229*AP(5,2) + C227*AP(7,1) + 
     .   C213*AP(7,2) - C225*AP(9,0) + C228*AP(9,1) + C226*AP(11,0)
      XT(2,2) = C1*AP(0,4) + C124*AP(2,3) - C224*AP(4,2) - 
     .   C222*AP(6,1) + C220*AP(8,0) + C197*AP(8,2) + C223*AP(10,1) + 
     .   C221*AP(12,0)
      XT(3,2) = C112*AP(3,3) + C213*AP(5,2) - C211*AP(7,1) + 
     .   C209*AP(9,0) + C212*AP(11,1) + C210*AP(13,0)
      XT(4,2) = -(C97*AP(4,3)) + C197*AP(6,2) - C195*AP(8,1) + 
     .   C193*AP(10,0) - C196*AP(12,1) + C194*AP(14,0)
      XT(-4,3) = C190*AP(1,0) + C109*AP(1,3) - C192*AP(3,2) + 
     .   C191*AP(5,1) - C189*AP(7,0)
      XT(-3,3) = C207*AP(0,1) + C205*AP(2,0)
      XT(-2,3) = C212*AP(1,1) - C112*AP(1,3) + C210*AP(3,0) - 
     .   C213*AP(3,2) + C211*AP(5,1) - C209*AP(7,0)
      XT(-1,3) = -(C192*AP(0,2)) + C100*AP(0,3) - C217*AP(2,1) - 
     .   C213*AP(2,2) + C215*AP(4,0) - C216*AP(4,1) + C214*AP(6,0)
      XT(0,3) = -(C208*AP(7,2)) - C219*AP(9,1) + C218*AP(11,0)
      XT(1,3) = -(C100*AP(4,3)) + C213*AP(6,2) + C216*AP(8,1) - 
     .   C192*AP(8,2) - C214*AP(10,0) - C217*AP(10,1) + C215*AP(12,0)
      XT(2,3) = C112*AP(3,3) + C213*AP(5,2) - C211*AP(7,1) + 
     .   C209*AP(9,0) + C212*AP(11,1) + C210*AP(13,0)
      XT(3,3) = C1*AP(0,4) - C122*AP(2,3) - C208*AP(4,2) + 
     .   C206*AP(6,1) - C204*AP(8,0) + C207*AP(12,1) + C205*AP(14,0)
      XT(4,3) = C109*AP(3,3) - C192*AP(5,2) + C191*AP(7,1) - 
     .   C189*AP(9,0) + C190*AP(15,0)
      XT(-4,4) = C186*AP(0,0)
      XT(-3,4) = C190*AP(1,0) - C109*AP(1,3) + C192*AP(3,2) - 
     .   C191*AP(5,1) + C189*AP(7,0)
      XT(-2,4) = -(C196*AP(0,1)) + C97*AP(0,3) + C194*AP(2,0) - 
     .   C197*AP(2,2) + C195*AP(4,1) - C193*AP(6,0)
      XT(-1,4) = -(C201*AP(1,1)) + C192*AP(1,2) + C199*AP(3,0) - 
     .   C200*AP(3,1) + C198*AP(5,0)
      XT(0,4) = C188*AP(8,2) - C203*AP(10,1) + C202*AP(12,0)
      XT(1,4) = -(C192*AP(7,2)) + C200*AP(9,1) - C198*AP(11,0) - 
     .   C201*AP(11,1) + C199*AP(13,0)
      XT(2,4) = -(C97*AP(4,3)) + C197*AP(6,2) - C195*AP(8,1) + 
     .   C193*AP(10,0) - C196*AP(12,1) + C194*AP(14,0)
      XT(3,4) = C109*AP(3,3) - C192*AP(5,2) + C191*AP(7,1) - 
     .   C189*AP(9,0) + C190*AP(15,0)
      XT(4,4) = C1*AP(0,4) - C120*AP(2,3) + C188*AP(4,2) - 
     .   C187*AP(6,1) + C185*AP(8,0) + C186*AP(16,0)
      GOTO 2000
C
C    Store results in the output array, and away we go
C
 2000 CONTINUE
      SCALE = SQRT( EIGHT*ALP1**(2*N1+1)*ALP2**(2*N2+1)/
     .              (ALP1+ALP2)**(2*N1+2*N2-1)*
     .              FACT(2*N1+2*N2-2)*RFACT(2*N1)*RFACT(2*N2) )
      IF(TRANSP) THEN
          DO 2100 M1=-L1,L1
              DO 2098 M2=-L2,L2
                  XP(M2,M1) = SCALE*XT(M1,M2)
 2098         CONTINUE
 2100     CONTINUE
      ELSE
          DO 2200 M1=-L1,L1
              DO 2198 M2=-L2,L2
                  XP(M1,M2) = SCALE*XT(M1,M2)
 2198         CONTINUE
 2200     CONTINUE
      ENDIF
      RETURN
11000 FORMAT(' COMPILE-TIME LIMIT EXCEEDED IN PSOVP, N1 = ',I4,
     .       ' L1 = ',I4,' N2 = ',I4,' L2 = ',I4)
11100 FORMAT(' UNIMPLEMENTED COMBINATION OF L1,L2 IN PSOVP, N1 = ',I4,
     .       ' L1 = ',I4,' N2 = ',I4,' L2 = ',I4)
      END
C
      SUBROUTINE PSOPR3(N1,L1,ALP1,N2,L2,ALP2,X,Y,Z,DX,LDX)
C                             x
C   Compute integrals <n,l,m|---|n',l',m'> (as well as analogous
C                            r^3
C   y and z integrals, with both orbitals centered on the same atom
C   and operator located on another one.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N1     - Major quantum number of the first set of orbitals
C      L1     - Orbital quantum number of the first set of orbitals
C      ALP1   - Orbital exponent of the first set of orbitals
C      N2     - Major quantum number of the second set of orbitals
C      L2     - Orbital quantum number of the second set of orbitals
C      ALP12  - Orbital exponent of the second set of orbitals
C      X,Y,Z  - Relative coordinates of the atom with x/r^3 operator.
C      DX     - Output integrals. Second index enumerates three
C               possible operators (in the order X,Y,Z), first and 
C               third index denote different orbitals. The arrangement 
C               is somewhat strange, but natural if you come to think 
C               about it.
C      LDX    - Leading dimension of the matrix DX.
C
C   Accessed common blocks:
C
C   Local storage:
C
C      (4*MAXL+1)*3*(MAX+1) = 255 DOUBLE PRECISION cells, which
C      is quite modest compared to the requirements of the called
C      routines.
C
C   Module logic:
C
C   Bugs:
C
C      L1 and L2 should nor exceed MAXL.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=4)
      PARAMETER (ZERO=0.D0)
      PARAMETER (SCALE=2.046653415892976976959103249778529721414D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION DX(LDX,3,2*L2+1)
      DIMENSION OVS(4*MAXL+1,3,0:MAXL)
C
      IF( L1.GT.MAXL .OR. L2.GT.MAXL ) THEN
          WRITE(NB6,11000) N1, L1, N2, L2
          STOP 'PSOPR3'
      ENDIF
C
C    Compute auxiliary integrals over the component Slater 
C    functions (aka charge distributions).
C
      DO 400 I=0,MIN(L1,L2)
          CALL PSORX(N1+N2-1,L1+L2-2*I,ALP1+ALP2,2,1,ZERO,X,Y,Z,3,
     .               OVS(1,1,I),4*MAXL+1)
  400 CONTINUE
C
C    Convert 'em to integrals over products of Slater functions,
C    once for each component of the interaction operator.
C
      CALL PSOVP(N1,L1,ALP1,N2,L2,ALP2,OVS(1,3,0),3*(4*MAXL+1),
     .           DX(1,1,1),3*LDX)
      CALL PSOVP(N1,L1,ALP1,N2,L2,ALP2,OVS(1,1,0),3*(4*MAXL+1),
     .           DX(1,2,1),3*LDX)
      CALL PSOVP(N1,L1,ALP1,N2,L2,ALP2,OVS(1,2,0),3*(4*MAXL+1),
     .           DX(1,3,1),3*LDX)
C
C    At this point, our integrals are over operators S_{1M}/r^3,
C    (where S_{1M} are real solid harmonics) rather than over
C    x/r^3 etc. Rescaling will take care of that.
C
      DO 900 M2=1,2*L2+1
          DO 800 I=1,3
              DO 700 M1=1,2*L1+1
                  DX(M1,I,M2) = SCALE * DX(M1,I,M2)
  700         CONTINUE
  800     CONTINUE
  900 CONTINUE
      RETURN
11000 FORMAT(' STATIC TEMPORARY WILL OVERFLOW IN PSOPR3, N1 = ',I3,
     .       ' L1 = ',I3,' N2 = ',I3,' L2 = ',I3)
      END
C
      SUBROUTINE PSOLB(NORBS,N,L,ALP,X,Y,Z,AP,LDA1,LDA3,AX,LDX)
C                         ^ ^ 
C   Compute integrals <mu|O L |nlm> (with |nlm> centered on atom B!=A)
C                                    ^     
C   using reduction to integrals <mu|O r_A|abc> with r_A been the
C   radius-vector relative to the position of atom A and |abc> 
C   (centered on the atom B) equal to |n-1,l-1,p>, |n-2,l,m> and |n-1,l,m>
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NORBS  - Number of <mu| orbitals.
C      N      - Major quantum number of the right set of orbitals.
C      L      - Orbital quantum number of the right set of orbitals.
C               L should be in the range 0 to N-1
C      ALP    - Orbital exponent of the right set of orbitals.
C      X,Y,Z  - Coordinates of the atom carrying L operator (i.e.,
C               atom A) relative to the atom carrying the right
C               set of orbitals (i.e., B)
C      AP     - Input <mu|O r_A|abc> integrals, organized as follows:
C               First index:  NORBS <mu| orbitals
C               Second index: Components of r_a operator, in the order
C                             X,Y,Z.
C               Third index:  2*b+1 |abc> orbitals, in the order of
C                             increasing c
C               Forth index:  Sets of |abc> orbitals.
C                             1: |n-1,l-1,c>
C                             2: |n-1,l,c>
C                             3: |n-2,l,c> (Not accessed if L.EQ.N-1)
C      LDA1   - First leading dimension of the AP matrix
C      LDA3   - Third leading dimension of the AP matrix (second and
C               forth dimensions are always equal to 3)
C      AX     - Output <mu|O L|nlm> integrals, organized as follows:
C               First index:  NORBS <mu| orbitals
C               Second index: components of the L operator, in the
C                             order X,Y,Z
C               Third index:  2*l+1 |nlm> orbitals, in the order
C                             of increasing c.
C
C   Accessed common blocks:
C
C   Local storage:
C
C      9*(MAXL+1)*MAXNOR (2025) DOUBLE PRECISION cells.
C
C   Module logic:
C
C      Transformation is performed in two steps. First, projected
C      integrals <mu|O t_a \partial\over{u}|nlm> (with t and u
C      taking all values of the set (x,y,z)) are computed
C      in the temporary array OD. Then, <mu|O L|nlm> integrals
C      are constructed.
C
C   Bugs:
C
C      L can't exceed MAXL. Likewise, NORBS can't be larger than
C      MAXNOR.
C      Some of the intermediate quanties are computed unnecessarily.
C      Integrals are always zero for some combinations of L1 and L2,
C      but are computed nonevertheless.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=4)
      PARAMETER (MAXNOR=25)
      PARAMETER (C1=1.7320508075688772935274463415058723669428D0)
      PARAMETER (C2=2.2360679774997896964091736687312762354406D0)
      PARAMETER (C3=1.290994448735805628393088466594133203611D0)
      PARAMETER (C4=2.581988897471611256786176933188266407222D0)
      PARAMETER (C5=3.240370349203930115482983718043998328853D0)
      PARAMETER (C6=2.6457513110645905905016157536392604257103D0)
      PARAMETER (C7=0.8366600265340755479781720257851874893928D0)
      PARAMETER (C8=2.898275349237887714743732831433954344628D0)
      PARAMETER (C9=3.346640106136302191912688103140749957571D0)
      PARAMETER (C10=2.049390153191919676644207736104210398147D0)
      PARAMETER (C11=3.549647869859769625540396974936970229049D0)
      PARAMETER (C12=4.242640687119285146405066172629094235709D0)
      PARAMETER (C13=3.674234614174767147295926112058837087949D0)
      PARAMETER (C14=3.D0)
      PARAMETER (C15=0.801783725737273154053660442639260564662D0)
      PARAMETER (C16=3.105295017040593980082570890822093940673D0)
      PARAMETER (C17=3.927922024247862862789754737481150133415D0)
      PARAMETER (C18=1.388730149658827192349850164875999283794D0)
      PARAMETER (C19=3.585685828003180919906451539079374954541D0)
      PARAMETER (C20=4.391550328268399307094730863080450853172D0)
      PARAMETER (C21=2.777460299317654384699700329751998567588D0)
      PARAMETER (C22=4.535573676110726726574198434810160729789D0)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION AP(LDA1,3,0:LDA3-1,3), AX(LDX,3,-L:L)
      DIMENSION OD(MAXNOR,3,3,-MAXL:MAXL)
C
      IF( L.LT.0 .OR. L.GT.MAXL .OR. NORBS.GT.MAXNOR ) THEN
          WRITE(NB6,11000) NORBS, N, L
          STOP 'PSOLB'
      ENDIF
C
C    Copy |n-1,l,m> integrals into intermediate array changing
C    sign in the progress. Since all three cartesian components
C    depend of |n-1,l,?> and |n-2,l,?> in the same way, only
C    one of the compoments need to be processed now.
C
      DO 300 M=-L,L
          DO 200 IOL=1,3
              DO 100 I=1,NORBS
                  OD(I,IOL,1,M) = -AP(I,IOL,L+M,2)
  100         CONTINUE
  200     CONTINUE
  300 CONTINUE
C
C    If necessary, mix in a tumble of |n-2,l,m> integrals
C
      IF( N.GT.L+1 ) THEN
          SC = (N-L-1)/SQRT(DBLE(N-1)*(N-1.5D0))
          DO 1300 M=-L,L
              DO 1200 IOL=1,3
                  DO 1100 I=1,NORBS
                      OD(I,IOL,1,M) = OD(I,IOL,1,M) + SC*AP(I,IOL,L+M,3)
 1100             CONTINUE
 1200         CONTINUE
 1300     CONTINUE
      ENDIF
C
C    Convert to three coordinate contributions
C
      SCX = X * ALP
      SCY = Y * ALP
      SCZ = Z * ALP
      DO 2300 M=-L,L
          DO 2200 IOL=1,3
              DO 2100 I=1,NORBS
                  OD(I,IOL,3,M) = SCZ * OD(I,IOL,1,M)
                  OD(I,IOL,2,M) = SCY * OD(I,IOL,1,M)
                  OD(I,IOL,1,M) = SCX * OD(I,IOL,1,M)
 2100         CONTINUE
 2200     CONTINUE
 2300 CONTINUE
C
C    Mix in serivatives of solid harmonics. This one is computer-generated.
C    (and is not intended for faint-hearted ;-)
C
      GOTO (3000,3100,3200,3300,3400), L+1
          WRITE(NB6,11100) L
          STOP 'PSOLB'
C  L = 0
 3000 CONTINUE
      GOTO 4000
C  L = 1
 3100 CONTINUE
      DO 3150 IT=1,NORBS
          DO 3140 IOL=1,3
              OD(IT,IOL,2,-1) = C1*AP(IT,IOL,0,1) + OD(IT,IOL,2,-1)
              OD(IT,IOL,3,0) = C1*AP(IT,IOL,0,1) + OD(IT,IOL,3,0)
              OD(IT,IOL,1,1) = C1*AP(IT,IOL,0,1) + OD(IT,IOL,1,1)
 3140     CONTINUE
 3150 CONTINUE
      GOTO 4000
C  L = 2
 3200 CONTINUE
      DO 3250 IT=1,NORBS
          DO 3240 IOL=1,3
              OD(IT,IOL,1,-2) = C2*AP(IT,IOL,0,1) + OD(IT,IOL,1,-2)
              OD(IT,IOL,2,-2) = C2*AP(IT,IOL,2,1) + OD(IT,IOL,2,-2)
              OD(IT,IOL,2,-1) = C2*AP(IT,IOL,1,1) + OD(IT,IOL,2,-1)
              OD(IT,IOL,3,-1) = C2*AP(IT,IOL,0,1) + OD(IT,IOL,3,-1)
              OD(IT,IOL,1,0) = -C3*AP(IT,IOL,2,1) + OD(IT,IOL,1,0)
              OD(IT,IOL,2,0) = -C3*AP(IT,IOL,0,1) + OD(IT,IOL,2,0)
              OD(IT,IOL,3,0) = C4*AP(IT,IOL,1,1) + OD(IT,IOL,3,0)
              OD(IT,IOL,1,1) = C2*AP(IT,IOL,1,1) + OD(IT,IOL,1,1)
              OD(IT,IOL,3,1) = C2*AP(IT,IOL,2,1) + OD(IT,IOL,3,1)
              OD(IT,IOL,1,2) = C2*AP(IT,IOL,2,1) + OD(IT,IOL,1,2)
              OD(IT,IOL,2,2) = -C2*AP(IT,IOL,0,1) + OD(IT,IOL,2,2)
 3240     CONTINUE
 3250 CONTINUE
      GOTO 4000
C  L = 3
 3300 CONTINUE
      DO 3350 IT=1,NORBS
          DO 3340 IOL=1,3
              OD(IT,IOL,1,-3) = C5*AP(IT,IOL,0,1) + OD(IT,IOL,1,-3)
              OD(IT,IOL,2,-3) = C5*AP(IT,IOL,4,1) + OD(IT,IOL,2,-3)
              OD(IT,IOL,1,-2) = C6*AP(IT,IOL,1,1) + OD(IT,IOL,1,-2)
              OD(IT,IOL,2,-2) = C6*AP(IT,IOL,3,1) + OD(IT,IOL,2,-2)
              OD(IT,IOL,3,-2) = C6*AP(IT,IOL,0,1) + OD(IT,IOL,3,-2)
              OD(IT,IOL,1,-1) = -C7*AP(IT,IOL,0,1) + OD(IT,IOL,1,-1)
              OD(IT,IOL,2,-1) = C8*AP(IT,IOL,2,1) + C7*AP(IT,IOL,4,1) + 
     .             OD(IT,IOL,2,-1)
              OD(IT,IOL,3,-1) = C9*AP(IT,IOL,1,1) + OD(IT,IOL,3,-1)
              OD(IT,IOL,1,0) = -C10*AP(IT,IOL,3,1) + OD(IT,IOL,1,0)
              OD(IT,IOL,2,0) = -C10*AP(IT,IOL,1,1) + OD(IT,IOL,2,0)
              OD(IT,IOL,3,0) = C11*AP(IT,IOL,2,1) + OD(IT,IOL,3,0)
              OD(IT,IOL,1,1) = 
     .            C8*AP(IT,IOL,2,1) - C7*AP(IT,IOL,4,1) + OD(IT,IOL,1,1)
              OD(IT,IOL,2,1) = -C7*AP(IT,IOL,0,1) + OD(IT,IOL,2,1)
              OD(IT,IOL,3,1) = C9*AP(IT,IOL,3,1) + OD(IT,IOL,3,1)
              OD(IT,IOL,1,2) = C6*AP(IT,IOL,3,1) + OD(IT,IOL,1,2)
              OD(IT,IOL,2,2) = -C6*AP(IT,IOL,1,1) + OD(IT,IOL,2,2)
              OD(IT,IOL,3,2) = C6*AP(IT,IOL,4,1) + OD(IT,IOL,3,2)
              OD(IT,IOL,1,3) = C5*AP(IT,IOL,4,1) + OD(IT,IOL,1,3)
              OD(IT,IOL,2,3) = -C5*AP(IT,IOL,0,1) + OD(IT,IOL,2,3)
 3340     CONTINUE
 3350 CONTINUE
      GOTO 4000
C  L = 4
 3400 CONTINUE
      DO 3450 IT=1,NORBS
          DO 3440 IOL=1,3
              OD(IT,IOL,1,-4) = C12*AP(IT,IOL,0,1) + OD(IT,IOL,1,-4)
              OD(IT,IOL,2,-4) = C12*AP(IT,IOL,6,1) + OD(IT,IOL,2,-4)
              OD(IT,IOL,1,-3) = C13*AP(IT,IOL,1,1) + OD(IT,IOL,1,-3)
              OD(IT,IOL,2,-3) = C13*AP(IT,IOL,5,1) + OD(IT,IOL,2,-3)
              OD(IT,IOL,3,-3) = C14*AP(IT,IOL,0,1) + OD(IT,IOL,3,-3)
              OD(IT,IOL,1,-2) = 
     .             -C15*AP(IT,IOL,0,1) + C16*AP(IT,IOL,2,1) + 
     .             OD(IT,IOL,1,-2)
              OD(IT,IOL,2,-2) = 
     .            C16*AP(IT,IOL,4,1) + C15*AP(IT,IOL,6,1) + 
     .             OD(IT,IOL,2,-2)
              OD(IT,IOL,3,-2) = C17*AP(IT,IOL,1,1) + OD(IT,IOL,3,-2)
              OD(IT,IOL,1,-1) = -C18*AP(IT,IOL,1,1) + OD(IT,IOL,1,-1)
              OD(IT,IOL,2,-1) = 
     .            C19*AP(IT,IOL,3,1) + C18*AP(IT,IOL,5,1) + 
     .             OD(IT,IOL,2,-1)
              OD(IT,IOL,3,-1) = C20*AP(IT,IOL,2,1) + OD(IT,IOL,3,-1)
              OD(IT,IOL,1,0) = -C21*AP(IT,IOL,4,1) + OD(IT,IOL,1,0)
              OD(IT,IOL,2,0) = -C21*AP(IT,IOL,2,1) + OD(IT,IOL,2,0)
              OD(IT,IOL,3,0) = C22*AP(IT,IOL,3,1) + OD(IT,IOL,3,0)
              OD(IT,IOL,1,1) = 
     .            C19*AP(IT,IOL,3,1) - C18*AP(IT,IOL,5,1) + 
     .             OD(IT,IOL,1,1)
              OD(IT,IOL,2,1) = -C18*AP(IT,IOL,1,1) + OD(IT,IOL,2,1)
              OD(IT,IOL,3,1) = C20*AP(IT,IOL,4,1) + OD(IT,IOL,3,1)
              OD(IT,IOL,1,2) = 
     .            C16*AP(IT,IOL,4,1) - C15*AP(IT,IOL,6,1) + 
     .             OD(IT,IOL,1,2)
              OD(IT,IOL,2,2) = 
     .            -C15*AP(IT,IOL,0,1) - C16*AP(IT,IOL,2,1) + 
     .             OD(IT,IOL,2,2)
              OD(IT,IOL,3,2) = C17*AP(IT,IOL,5,1) + OD(IT,IOL,3,2)
              OD(IT,IOL,1,3) = C13*AP(IT,IOL,5,1) + OD(IT,IOL,1,3)
              OD(IT,IOL,2,3) = -C13*AP(IT,IOL,1,1) + OD(IT,IOL,2,3)
              OD(IT,IOL,3,3) = C14*AP(IT,IOL,6,1) + OD(IT,IOL,3,3)
              OD(IT,IOL,1,4) = C12*AP(IT,IOL,6,1) + OD(IT,IOL,1,4)
              OD(IT,IOL,2,4) = -C12*AP(IT,IOL,0,1) + OD(IT,IOL,2,4)
 3440     CONTINUE
 3450 CONTINUE
      GOTO 4000
C
C    Take a vector product, scaling result on the fly, and away we go.
C
 4000 CONTINUE
      SC = ALP/SQRT(DBLE(N)*(N-0.5D0))
      DO 4500 M=-L,L
          DO 4400 IT=1,NORBS
              AX(IT,1,M) = SC*(OD(IT,2,3,M) - OD(IT,3,2,M))
              AX(IT,2,M) = SC*(OD(IT,3,1,M) - OD(IT,1,3,M))
              AX(IT,3,M) = SC*(OD(IT,1,2,M) - OD(IT,2,1,M))
 4400     CONTINUE
 4500 CONTINUE
      RETURN
11000 FORMAT(' STATIC TEMPORARY WILL OVERFLOW IN PSOLB. NORBS = ',I6,
     .       ' N = ',I4,' L = ',I4)
11100 FORMAT(' ASSERTION FAILURE IN PSOLB. L = ', I4 )
      END
C
      SUBROUTINE PSLA3(NA,LA,ZA,NB,LB,ZB,RX,RY,RZ,DL,LDL)
C
C                    ^A
C               A    L       B
C   Compute < mu  | ---- | nu > two-center integrals for given n and l.
C                     3
C                    r
C                     A
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NA     - Major quantum number of the left-hand center.
C      LA     - Orbital quantum number of the left-hand center.
C      ZA     - Orbital exponent of the left-hand center.
C               Operator is on the left-hand atom.
C      NB     - Major quantum number of the right-hand center.
C      LB     - Orbital quantum number of the right-hand center.
C      ZB     - Orbital exponent of the right-hand center.
C      RX,RY,RZ
C             - Relative coordinates of the atom A (left-hand one) in
C               atomic units.
C      DL     - Output integrals.
C               First index:  enumerates left-hand orbitals (increasing
C                   magnetic quantum number order).
C               Second index: 1 to 3: components of the L^A operator,
C                   in the X,Y,Z order.
C               Third index: enumerated right-hand orbitals (increasing
C                   magnetic quantum number order).
C      LDL    - Leading dimension of the DL matrix.
C
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C      ca. 56*(MAXL+.5)**2 DOUBLE PRECISION storage units. With MAXL=2,
C      the it's about 400 DP words.
C
C   Module logic:
C
C      Integrals are reduced to elementary integrals over r^-3 operator
C      using PSOLB. Elementary integrals are evaluated by PSORX and PSODS.
C      This is a little scaled-back version of PSRLA3, which is full of
C      really helpful comments (if you have time to hack through it).
C
C   Bugs:
C
C      Bugs? Which bugs? We donna wanna any stinking bugs!
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=3)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION DL(LDL,3,2*LB+1)
C
      DIMENSION ELT(2*MAXL+1,2*(MAXL+1)+1,2)
      DIMENSION ERT(2*MAXL+1,2*MAXL+1,3)
      DIMENSION ER (2*MAXL+1,3,2*MAXL+1,3)
C
      DIMENSION ILA(2), INB(3), ILB(3)
      DATA ILA/-1, 1/ INB/-1,-1,-2/ ILB/-1, 0, 0/
C
      IF( LA.GT.MAXL .OR. LB.GT.MAXL ) THEN
          WRITE(NB6,11000) LA, LB
          STOP 'PSLA3'
      ENDIF
C
      DO 2000 JB=1,3
          NBT = NB + INB(JB)
          LBT = LB + ILB(JB)
          IF( LBT.GE.0 .AND. LBT.LE.NBT ) THEN
              DO 1000 JA=1,2
                  NAT = NA + 1
                  LAT = LA + ILA(JA)
C                Elementary integrals. Since PSORX wants operator to be on the
C                right side, the integrals are transposed - i.e., B is on the
C                left, A is on the right (this incidentally makes relative
C                coordinates' sign right).
                  IF( LAT.GE.0 .AND. LAT.LE.NAT )
     .            CALL PSORX(NBT,LBT,ZB,NAT,LAT,ZA,RX,RY,RZ,3,
     .                       ELT(1,1,JA),2*MAXL+1)
 1000         CONTINUE
C            Produce intermediate <mu^A r^A|...|...> integrals. At this point,
C            B is still on the left, A is on the right - PSODS wants it this way.
              CALL PSODS(2*LBT+1,NA,LA,ZA,ELT(1,1,1),2*MAXL+1,
     .                                    ELT(1,1,2),2*MAXL+1,
     .                   ERT,2*MAXL+1,2*MAXL+1)
C            Transpose dipole intermediates are reshuffle 'em into the order PSOLB
C            would like, that is: first index are A orbitals, second index are r^A
C            components, third index are B orbitals and forth index are B orbitals
C            sets, that is JB.
              DO 1500 IX=1,3
                  DO 1490 MA=1,2*LA+1
                      DO 1480 MB=1,2*LBT+1
                          ER(MA,IX,MB,JB) = ERT(MB,MA,IX)
 1480                 CONTINUE
 1490             CONTINUE
 1500         CONTINUE
          ENDIF
 2000 CONTINUE
C    Final shot - ask PSOLB to evaluate our integrals, and away we go!
      CALL PSOLB(2*LA+1,NB,LB,ZB,RX,RY,RZ,ER,2*MAXL+1,2*MAXL+1,DL,LDL)
C
      RETURN
11000 FORMAT(' UNSUPPORTED PAIR LA = ',I3,' LB = ',I3,
     .       ' WAS REQUESTED FROM PSLA3.')
      END
C
      SUBROUTINE PSRLA3(NA,LA,ZA,NB,LB,ZB,RX,RY,RZ,DL,LDL)
C
C                      ^A
C               A      L       B
C   Compute < mu  r | ---- | nu > two-center integrals for given n and l.
C                  A    3
C                      r
C                       A
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NA     - Major quantum number of the left-hand center.
C      LA     - Orbital quantum number of the left-hand center.
C      ZA     - Orbital exponent of the left-hand center.
C               Operator is on the left-hand atom.
C      NB     - Major quantum number of the right-hand center.
C      LB     - Orbital quantum number of the right-hand center.
C      ZB     - Orbital exponent of the right-hand center.
C      RX,RY,RZ
C             - Relative coordinates of the atom A (left-hand one) in
C               atomic units.
C      DL     - Output integrals.
C               First index:  enumerates left-hand orbitals (increasing
C                   magnetic quantum number order).
C               Second index, 1 to 3: components of the r_A operator,
C                   in the X,Y,Z order.
C               Third index: 1 to 3: components of the L^A operator,
C                   in the X,Y,Z order.
C               Forth index: enumerated right-hand orbitals (increasing
C                   magnetic quantum number order).
C      LDL    - Leading dimension of the DL matrix.
C
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C      ca. 144*(MAXL+1)**2 - 36 DOUBLE PRECISION storage units. With MAXL=2,
C      it's about 1300 DP words.
C
C   Module logic:
C
C      Integrals are reduced to elementary integrals over r^-3 operator
C      (evaluated by PSORX) using PSODS (for the dipole moment reduction)
C      and PSOLB (for the displaced orbital moment reduction). The complete
C      reduction graph is a bit involved:
C
C      Final result is produced by PSODS, which wants to see integrals
C          over <NA+1,LA+1,?| and <NA+1,LA-1,?| and the rest of the
C          integral as it is.
C      These integrals are produced by PSOLB, which wants to see
C          <NA+2,LA+2,?|, <NA+2,LA+0,?| and <NA+2,LA-2,?| on the left
C          side, ra^-3 as an operator and |NB-1,LB-1,?>, |NB-2,LB-0,?>
C          and |NB-2,LB-0,?> on the right-hand side.
C      Elementary integrals over r^-3 are evaluated by PSORX. All in all,
C          we might need up to nine elementary integrals blocks (with
C          orbital moments as high as LA+2 and LB) to get a single orbital
C          moment block of the final integrals.
C
C
C      The naming scheme for the intermediate quantities is as follows:
C
C      Elementary integrals are indexed by JNA and JNB:
C          JNA is 1 = |NA+2,LA-2,?>
C                 2 = |NA+2,LA+0,?>
C                 3 = |NA+2,LA+2,?>
C
C              or 1 = |NA+2,LA-1,?>
C                 2 = |NA+2,LA+1,?>
C
C          JNB is 1 = |NB-1,LB-1,?>
C                 2 = |NB-1,LB-0,?>
C                 3 = |NB-2,LB-0,?>
C
C      ELM is an array of the elementary r^-3 integrals. (B on left, A on right).
C      EDT is an array of the r_A r^-3 integrals. (B on left, A on right).
C      EDP is an array of the <A r_A|r^-3|B> integrals. (A is on the left, B is
C          on the right). Second index is a r_A component, forth index is JNB,
C          fith index is JNA (variant 2).
C      ELT is an array of the <A|L^A r^-3|B> integrals. (A is on the left, B is
C          on the right). Second index is the projection of the L^A operator.
C      EL  is an array of the <A|L^A r^-3|B> integrals. (B is on the left, A is
C          on the right). Third index is the JNA (variant 2), forth index is
C          the projection of the L^A operator.
C      DLT is an array of the <A r_A|L^A_j r^-3|B> integrals. (B is on the left,
C          A is on the right). Third index in the dipole moment component.
C
C   Bugs:
C
C      I hope there are none...
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
C
      DIMENSION DL(LDL,3,3,2*LB+1)
C
      DIMENSION ELM(2*MAXL+1,2*(MAXL+2)+1,3)
      DIMENSION EDT(2*MAXL+1,2*(MAXL+1)+1,3)
      DIMENSION EDP(2*(MAXL+1)+1,3,2*MAXL+1,3,2)
      DIMENSION ELT(2*(MAXL+1)+1,3,2*MAXL+1)
      DIMENSION EL (2*MAXL+1,2*(MAXL+1)+1,2,3)
      DIMENSION DLT(2*MAXL+1,2*MAXL+1,3)
      DIMENSION INA(3), ILA(3), INB(3), ILB(3)
      DATA INA/ 2, 2, 2/ ILA/-2, 0, 2/
      DATA INB/-1,-1,-2/ ILB/-1, 0, 0/
C
      IF( LA.GT.MAXL .OR. LB.GT.MAXL ) THEN
          WRITE(NB6,11000) LA, LB
          STOP 'PSRLA3'
      ENDIF
C    Elementary integrals. Orbitals belonging to atom A
C    are on the right-hand side, since PSORX wants it
C    this way.
      DO 500 JNB=1,3
          NBT = NB + INB(JNB)
          LBT = LB + ILB(JNB)
          IF( NBT.LT.0 .OR. LBT.LT.0 .OR. LBT.GT.NBT ) GOTO 499
*         write(*,*)' jnb = ',jnb,' nbt = ',nbt,' lbt = ',lbt
          DO 200 JNA=1,3
              NAT = NA + INA(JNA)
              LAT = LA + ILA(JNA)
              IF( NAT.GE.0 .AND. LAT.GE.0 .AND. LAT.LE.NAT ) THEN
*                 write(*,*)'jna = ',jna,' nat = ',nat,' lat = ',lat
                  CALL PSORX(NBT,LBT,ZB,NAT,LAT,ZA,RX,RY,RZ,3,
     .                       ELM(1,1,JNA),2*MAXL+1)
              ENDIF
  200     CONTINUE
C        Dipole contraction of the elementary integrals
C        B is still on the left, A is on the right.
          DO 400 JNA=1,2
              NAT = NA + INA(JNA) - 1
              LAT = LA + ILA(JNA) + 1
              IF( LAT.LT.0 ) GOTO 399
*             write(*,*)' jna = ',jna,' nat = ',nat,' lat = ',lat
              CALL PSODS(2*LBT+1,NAT,LAT,ZA,ELM(1,1,JNA  ),2*MAXL+1,
     .                                      ELM(1,1,JNA+1),2*MAXL+1,
     .                   EDT,2*MAXL+1,2*(MAXL+1)+1)
C            Tranposition of the dipole intermediates
C            Now, A is on the left, B is on the right.
              DO 300 IX=1,3
                  DO 290 MA=1,2*LAT+1
                      DO 280 MB=1,2*LBT+1
                          EDP(MA,IX,MB,JNB,JNA) = EDT(MB,MA,IX)
  280                 CONTINUE
  290             CONTINUE
  300         CONTINUE
  399         CONTINUE
  400     CONTINUE
  499     CONTINUE
  500 CONTINUE
C        Loop over <N+2,L-1,?| and <N+2,L+1,?| integrals on the left
          DO 900 JNA=1,2
              NAT = NA + INA(JNA) - 1
              LAT = LA + ILA(JNA) + 1
              IF( LAT.LT.0 ) GOTO 899
C            Displaced orbital moment contraction. A is on the left, B is
C            on the right.
*             write(*,*)' jna = ',jna,' nat = ',nat,' lat = ',lat
              CALL PSOLB(2*LAT+1,NB,LB,ZB,RX,RY,RZ,EDP(1,1,1,1,JNA),
     .                   2*(MAXL+1)+1,2*MAXL+1,ELT,2*(MAXL+1)+1)
C            Interchange A and B once again.
              DO 700 MB=1,2*LB+1
                  DO 690 IX=1,3
                      DO 680 MA=1,2*LAT+1
                          EL(MB,MA,JNA,IX) = ELT(MA,IX,MB)
  680                 CONTINUE
  690             CONTINUE
  700         CONTINUE
  899         CONTINUE
  900     CONTINUE
C        Loop over components of the orbital moment operator
          DO 1200 IX=1,3
C            The last dipole contaction. B is on the left, A is on the right.
              CALL PSODS(2*LB+1,NA,LA,ZA,EL(1,1,1,IX),2*MAXL+1,
     .                   EL(1,1,2,IX),2*MAXL+1,DLT,2*MAXL+1,2*MAXL+1)
C            Transpose and store computed integrals.
              DO 1100 IY=1,3
                  DO 1090 MA=1,2*LA+1
                      DO 1080 MB=1,2*LB+1
                          DL(MA,IY,IX,MB) = DLT(MB,MA,IY)
 1080                 CONTINUE
 1090             CONTINUE
 1100         CONTINUE
 1200     CONTINUE
C
      RETURN
11000 FORMAT(' UNSUPPORTED PAIR LA = ',I3,' LB = ',I3,
     .       ' WAS REQUESTED FROM PSRLA3.')
      END
C
      SUBROUTINE PSRAB3(NA,LA,ZA,NB,LB,ZB,X,Y,Z,DD,LDD)
C
C   Compute two-center integrals <mu^A r^A|(r^A)^-3|nu^B r^B>.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NA     - Major quantum number of the left-hand set of orbitals
C      LA     - Orbital quantum number of the left-hand set of orbitals
C      ZA     - Orbital exponent of the left-hand set of orbitals
C      NB     - Major quantum number of the right-hand set of orbitals
C      LB     - Orbital quantum number of the right-hand set of orbitals
C      ZB     - Orbital exponent of the right-hand set of orbitals
C      X,Y,Z  - Relative coordinates of the left-hand atom (and hence,
C               of the operator) in atomic units.
C      DD     - Output integrals.
C               First index: orbitals on the left-hand atom, increasing
C                   M order.
C               Second index: component of the left-hand radius-vector
C                   operator.
C               Third index: component of the right-hand radius-vector
C                   operator.
C               Forth index: orbitals on the right-hand atom, increasing
C                   M order.
C      LDD    - Leading dimension of the matrix DD.
C
C   Accessed common blocks:
C
C      Nope.
C
C   Local storage:
C
C      56*MAXL**2 + 108*MAXL + 48 DOUBLE PRECISION storage units. With
C      MAXL=3, it's about 880 DP words.
C
C   Module logic:
C
C      Double dipole moment reduction is used.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=3)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION DD(LDD,3,3,2*LB+1)
C
      DIMENSION EXT(2*(MAXL+1)+1,2*(MAXL+1)+1,2)
      DIMENSION EDT(2*(MAXL+1)+1,2*MAXL+1,3)
      DIMENSION ED (2*MAXL+1,2*(MAXL+1)+1,3,2)
      DIMENSION DDX(2*MAXL+1,2*MAXL+1,3)
C
      IF( LA.GT.MAXL .OR. LB.GT.MAXL ) THEN
          WRITE(NB6,11000) LA, LB
          STOP 'PSRAB3'
      ENDIF
C
      DO 2000 IB=1,2
          LBT = LB + 2*IB-3
          IF( LBT.LT.0 ) GOTO 1999
          DO 1000 IA=1,2
              LAT = LA + 2*IA-3
C            Evaluate elementary <NB+1,LB+/-1|...|NA+1,LA+/-1> integrals
C            The integrals are transposed at this point, since PSORX
C            works this way.
              IF( LAT.GE.0 )
     .        CALL PSORX(NB+1,LBT,ZB,NA+1,LAT,ZA,X,Y,Z,3,
     .                   EXT(1,1,IA),2*(MAXL+1)+1)
 1000     CONTINUE
C        Construct r^A dipole integrals. B is still on the left, A is
C        on the right.
          CALL PSODS(2*LBT+1,NA,LA,ZA,EXT(1,1,1),2*(MAXL+1)+1,
     .                                EXT(1,1,2),2*(MAXL+1)+1,
     .               EDT,2*(MAXL+1)+1,2*MAXL+1)
C        Transpose r^A dipole integrals. After this point, A is on the
C        left, B is on the right, third index is r_A component and forth
C        index is connected to LB increment.
          DO 1500 IX=1,3
              DO 1490 MA=1,2*LA+1
                  DO 1480 MB=1,2*LBT+1
                      ED(MA,MB,IX,IB) = EDT(MB,MA,IX)
 1480             CONTINUE
 1490         CONTINUE
 1500     CONTINUE
 1999     CONTINUE
 2000 CONTINUE
C    Combine |NB+1,LB-1> and |NB+1,LB+1> integrals to get r^B components.
      DO 3000 IX=1,3
          CALL PSODS(2*LA+1,NB,LB,ZB,ED(1,1,IX,1),2*MAXL+1,
     .                               ED(1,1,IX,2),2*MAXL+1,
     .               DDX,2*MAXL+1,2*MAXL+1)
C        Stuff final integrals into the output buffer.
          DO 2900 IY=1,3
              DO 2890 MB=1,2*LB+1
                  DO 2880 MA=1,2*LA+1
                      DD(MA,IX,IY,MB) = DDX(MA,MB,IY)
 2880             CONTINUE
 2890         CONTINUE
 2900     CONTINUE
 3000 CONTINUE
C    Whoa!
      RETURN
11000 FORMAT(' LA (',I5,') OR LB (',I5,') IS TOO LARGE IN PSRAB3.')
      END
C
      SUBROUTINE PSLBR3(N1,L1,ALP1,N2,L2,ALP2,X,Y,Z,DX,LDX)
C                             L
C   Compute integrals <n,l,m|---|n',l',m'> with both orbitals
C                            r^3
C   located on the same center (and operator on the other one)
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      N1     - Major quantum number of the first set of orbitals
C      L1     - Orbital quantum number of the first set of orbitals
C      ALP1   - Orbital exponent of the first set of orbitals
C      N2     - Major quantum number of the second set of orbitals
C      L2     - Orbital quantum number of the second set of orbitals
C      ALP12  - Orbital exponent of the second set of orbitals
C      X,Y,Z  - Relative coordinates of the atom with L/r^3 operator.
C      DX     - Output integrals. Second index enumerates three
C               possible L components (in the order X,Y,Z), first and 
C               third index denote different orbitals. The arrangement 
C               is somewhat strange, but natural if you come to think 
C               about it.
C      LDX    - Leading dimension of the matrix DX.
C
C   Accessed common blocks:
C
C   Local storage:
C
C      9*(2*MAXL+1)**2 (729) DOUBLE PRECISION cells.
C
C   Module logic:
C
C   Bugs:
C
C      L1 and L2 should not exceed MAXL.
C      Values computed by PSLBR3 are accurate to at most 10^-10 a.u.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=4)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION DX(-L1:LDX-L1-1,3,-L2:L2)
      DIMENSION AP(2*MAXL+1,3,2*MAXL+1,3)
C
      IF( L1.GT.MAXL .OR. L2.GT.MAXL ) THEN
          WRITE(NB6,11000) N1, L1, N2, L2
          STOP 'PSLBR3'
      ENDIF
C
C    First, compute auxiliary integrals over products of Slater
C    function. We want the following set of integrals:
C       <n1,l1,m1|x/r^3|n2-1,l2-1,m2> (only needed if L2>0)
C       <n1,l1,m1|x/r^3|n2-1,l2,m2>
C       <n1,l1,m1|x/r^3|n2-2,l2,m2> (only needed if N2>L2+1)
C    as well as corresponding y and z versions.
C
      IF( L2.GT.0 )
     .CALL PSOPR3(N1,L1,ALP1,N2-1,L2-1,ALP2,X,Y,Z,AP(1,1,1,1),2*MAXL+1)
      CALL PSOPR3(N1,L1,ALP1,N2-1,L2,  ALP2,X,Y,Z,AP(1,1,1,2),2*MAXL+1)
      IF( N2.GT.L2+1 )
     .CALL PSOPR3(N1,L1,ALP1,N2-2,L2,  ALP2,X,Y,Z,AP(1,1,1,3),2*MAXL+1)
C
C    Then, convert 'em to what our heart desire.
C
      CALL PSOLB(2*L1+1,N2,L2,ALP2,X,Y,Z,AP,2*MAXL+1,2*MAXL+1,DX,LDX)
C
      RETURN
11000 FORMAT(' STATIC TEMPORARY WILL OVERFLOW IN PSLBR3. N1 = ',I4,
     .       ' L1 = ',I4,' N2 = ',I4,' L2 = ',I4)
      END
C
      DOUBLE PRECISION FUNCTION PSR31C(NA,ALPA,NB,ALPB)
C                             
C   Computes one-center integrals <\mu|r^-3|\nu> over STOs.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NA     - Major quantum number of the first orbital
C      ALPA   - Orbital exponent at the first atom
C      NB     - Major quantum number of the second orbital
C      ALPB   - Orbital exponent at the second atom
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
      PARAMETER (MAXJ=30)
      PARAMETER (MXJ=90)
      PARAMETER (TWO=2.D0)
C
      COMMON /PS3JDT/FACT(0:MXJ+1),POW2(0:MXJ+1),RINT(0:MXJ+1),
     .               RFACT(0:MXJ+1), RSFACT(0:MXJ+1), GHALF(0:MXJ+1)
     ./PSPRT / NB6
      SAVE /PSPRT /
C
      IF( NA.LT.0 .OR. NA.GE.MAXJ/2 .OR. NB.LT.0 .OR. NB.GE.MAXJ/2 
     .            .OR. NA+NB.LT.3 ) THEN
          WRITE(NB6,11000) NA, NB
          STOP 'PSR31C'
      ENDIF
      PSR31C = FACT(NA+NB-3)*RSFACT(2*NA)*RSFACT(2*NB)
     .       * SQRT(((TWO*ALPA)**(2*NA+1))*((TWO*ALPB)**(2*NB+1)))
     .       / ((ALPA+ALPB)**(NA+NB-2))
      RETURN
11000 FORMAT(' MAJOR QUANTUM NUMBERS ',I3, ' AND ', I3, 
     .       ' IN PSR31C ARE OUTSIDE PERMISSIBLE BONDS.' )
      END
C
      DOUBLE PRECISION FUNCTION PSOV1C(NA,ALPA,NB,ALPB)
C                             
C   Computes one-center overlap integrals over STOs.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NA     - Major quantum number of the first orbital
C      ALPA   - Orbital exponent at the first atom
C      NB     - Major quantum number of the second orbital
C      ALPB   - Orbital exponent at the second atom
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
      PARAMETER (MAXJ=30)
      PARAMETER (MXJ=90)
      PARAMETER (TWO=2.D0)
C
      COMMON /PS3JDT/FACT(0:MXJ+1),POW2(0:MXJ+1),RINT(0:MXJ+1),
     .               RFACT(0:MXJ+1), RSFACT(0:MXJ+1), GHALF(0:MXJ+1)
     ./PSPRT / NB6
      SAVE /PSPRT /
C
      IF( NA.LT.1.OR.NA.GE.MAXJ/2.OR.NB.LT.1.OR.NB.GE.MAXJ/2 ) THEN
          WRITE(NB6,11000) NA, NB
          STOP 'PSOV1C'
      ENDIF
      PSOV1C = FACT(NA+NB)*RSFACT(2*NA)*RSFACT(2*NB)
     .       * SQRT(((TWO*ALPA)**(2*NA+1))*((TWO*ALPB)**(2*NB+1)))
     .       / ((ALPA+ALPB)**(NA+NB+1))
      RETURN
11000 FORMAT(' MAJOR QUANTUM NUMBERS ',I3, ' AND ', I3, 
     .       ' IN PSOV1C ARE OUTSIDE PERMISSIBLE BONDS.' )
      END
C
      SUBROUTINE PSSX94(NORBSA,NORBSB,F,LDF)
C
C   Sort one-electron integrals over AOs into the MNDO94 order.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      NORBSA - Number of AOs at the first center
C      NORBSB - Number of AOs at the second center
C      F      - Input: integrals in increasing-M order
C               Output: integrals in MNDO94 order
C      LDF    - Leading dimension of the F matrix
C
C   Accessed common blocks:
C
C   Local storage:
C
C      MAXORB DOUBLE PRECISION cells.
C
C   Module logic:
C
C      is obvious. ;-')
C
C   Bugs:
C
C      Some of the copying is, strictly speaking, unncesary.
C      If NORBSA or NORBSB is not 1,4 or 9, all the hell will
C      broke lose.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXORB=9)
      COMMON /PSPRT / NB6
      SAVE /PSPRT /
      DIMENSION F(LDF,NORBSB)
      DIMENSION TEMP(MAXORB)
C
      IF( NORBSA.GT.MAXORB .OR. NORBSB.GT.MAXORB ) THEN
          WRITE(NB6,11000) NORBSA, NORBSB
          STOP 'PSSX94'
      ENDIF
C
C   Swap orbitals at the second atom.
C
      IF( NORBSB.GT.1 ) THEN
          DO 100 I=1,NORBSA
              TEMP(I) = F(I,2)
  100     CONTINUE
          DO 110 I=1,NORBSA
              F(I,2)  = F(I,4)
  110     CONTINUE
          DO 120 I=1,NORBSA
              F(I,4)  = F(I,3)
  120     CONTINUE
          DO 130 I=1,NORBSA
              F(I,3)  = TEMP(I)
  130     CONTINUE
          IF( NORBSB.GT.4 ) THEN
              DO 200 I=1,NORBSA
                  TEMP(I) = F(I,5)
  200         CONTINUE
              DO 210 I=1,NORBSA
                  F(I,5)  = F(I,9)
  210         CONTINUE
              DO 220 I=1,NORBSA
                  F(I,9)  = TEMP(I)
  220         CONTINUE
              DO 230 I=1,NORBSA
                  TEMP(I) = F(I,6)
  230         CONTINUE
              DO 240 I=1,NORBSA
                  F(I,6)  = F(I,8)
  240         CONTINUE
              DO 250 I=1,NORBSA
                  F(I,8)  = TEMP(I)
  250         CONTINUE
          ENDIF
      ENDIF
C
C  Swap orbitals at the first atom.
C
      IF( NORBSA.GT.1 ) THEN
          DO 500 I=1,NORBSB
              TEMP(I) = F(2,I)
  500     CONTINUE
          DO 510 I=1,NORBSB
              F(2,I)  = F(4,I)
  510     CONTINUE
          DO 520 I=1,NORBSB
              F(4,I)  = F(3,I)
  520     CONTINUE
          DO 530 I=1,NORBSB
              F(3,I)  = TEMP(I)
  530     CONTINUE
          IF( NORBSA.GT.4 ) THEN
              DO 600 I=1,NORBSB
                  TEMP(I) = F(5,I)
  600         CONTINUE
              DO 610 I=1,NORBSB
                  F(5,I)  = F(9,I)
  610         CONTINUE
              DO 620 I=1,NORBSB
                  F(9,I)  = TEMP(I)
  620         CONTINUE
              DO 630 I=1,NORBSB
                  TEMP(I) = F(6,I)
  630         CONTINUE
              DO 640 I=1,NORBSB
                  F(6,I)  = F(8,I)
  640         CONTINUE
              DO 650 I=1,NORBSB
                  F(8,I)  = TEMP(I)
  650         CONTINUE
          ENDIF
      ENDIF
C
      RETURN
11000 FORMAT(' SIZE OF THE ATOMIC INTEGRALS BLOCK IS TOO BIG (',I3,
     .       ',',I3,') IN PSSX94.' )
      END
C
      BLOCK DATA PS3JBD
C
C   Initialized common block PS3JDT with values used by PS3J.
C
C   Coded by: Serge Pachkovsky, with the help of faithful Mathematica.
C
C   Defined common blocks:
C
C     PS3JDT - Various interesting constants, namely:
C         FACT   = n!
C         POW2   = 2^n
C         RINT   = 1/n
C         RFACT  = 1/n!
C         RSFACT = 1/Sqrt[n!]
C         GHALF  = Gamma[n+1/2]/Gamma[1/2]
C
C   Bugs:
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXJ=90)
C
      COMMON /PS3JDT/FACT(0:MXJ+1),POW2(0:MXJ+1),RINT(0:MXJ+1),
     .               RFACT(0:MXJ+1), RSFACT(0:MXJ+1), GHALF(0:MXJ+1)
C
      DATA FACT  (0) / 1.D+0/
      DATA POW2  (0) / 1.D+0/
      DATA RINT  (0) / 0.D+0/
      DATA RFACT (0) / 1.D+0/
      DATA RSFACT(0) / 1.D+0/
      DATA GHALF (0) / 1.D+0/
      DATA FACT  (1) / 1.D+0/
      DATA POW2  (1) / 2.D+0/
      DATA RINT  (1) / 1.D+0/
      DATA RFACT (1) / 1.D+0/
      DATA RSFACT(1) / 1.D+0/
      DATA GHALF (1) / 0.5D+0/
      DATA FACT  (2) / 2.D+0/
      DATA POW2  (2) / 4.D+0/
      DATA RINT  (2) / 0.5D+0/
      DATA RFACT (2) / 0.5D+0/
      DATA RSFACT(2) / 0.707106781186547524400844362104849D+0/
      DATA GHALF (2) / 0.75D+0/
      DATA FACT  (3) / 6.D+0/
      DATA POW2  (3) / 8.D+0/
      DATA RINT  (3) / 0.333333333333333333333333333333333D+0/
      DATA RFACT (3) / 0.166666666666666666666666666666667D+0/
      DATA RSFACT(3) / 0.408248290463863016366214012450982D+0/
      DATA GHALF (3) / 1.875D+0/
      DATA FACT  (4) / 24.D+0/
      DATA POW2  (4) / 16.D+0/
      DATA RINT  (4) / 0.25D+0/
      DATA RFACT (4) / 41.666666666666666666666666666666667D-3/
      DATA RSFACT(4) / 0.204124145231931508183107006225491D+0/
      DATA GHALF (4) / 6.5625D+0/
      DATA FACT  (5) / 0.12D+3/
      DATA POW2  (5) / 32.D+0/
      DATA RINT  (5) / 0.2D+0/
      DATA RFACT (5) / 8.333333333333333333333333333333333D-3/
      DATA RSFACT(5) / 91.287092917527685576161630466800356D-3/
      DATA GHALF (5) / 29.53125D+0/
      DATA FACT  (6) / 0.72D+3/
      DATA POW2  (6) / 64.D+0/
      DATA RINT  (6) / 0.166666666666666666666666666666667D+0/
      DATA RFACT (6) / 1.388888888888888888888888888888889D-3/
      DATA RSFACT(6) / 37.267799624996494940152894478854604D-3/
      DATA GHALF (6) / 0.162421875D+3/
      DATA FACT  (7) / 5.04D+3/
      DATA POW2  (7) / 0.128D+3/
      DATA RINT  (7) / 0.142857142857142857142857142857143D+0/
      DATA RFACT (7) / 0.198412698412698412698412698412698D-3/
      DATA RSFACT(7) / 14.085904245475276291826972122765755D-3/
      DATA GHALF (7) / 1.0557421875D+3/
      DATA FACT  (8) / 40.32D+3/
      DATA POW2  (8) / 0.256D+3/
      DATA RINT  (8) / 0.125D+0/
      DATA RFACT (8) / 24.801587301587301587301587301587302D-6/
      DATA RSFACT(8) / 4.980119205559973499870071582054687D-3/
      DATA GHALF (8) / 7.91806640625D+3/
      DATA FACT  (9) / 0.36288D+6/
      DATA POW2  (9) / 0.512D+3/
      DATA RINT  (9) / 0.111111111111111111111111111111111D+0/
      DATA RFACT (9) / 2.755731922398589065255731922398589D-6/
      DATA RSFACT(9) / 1.660039735186657833290023860684896D-3/
      DATA GHALF (9) / 67.303564453125D+3/
      DATA FACT  (10) / 3.6288D+6/
      DATA POW2  (10) / 1.024D+3/
      DATA RINT  (10) / 0.1D+0/
      DATA RFACT (10) / 0.275573192239858906525573192239859D-6/
      DATA RSFACT(10) / 0.524950656957260037797939633658583D-3/
      DATA GHALF (10) / 0.6393838623046875D+6/
      DATA FACT  (11) / 39.9168D+6/
      DATA POW2  (11) / 2.048D+3/
      DATA RINT  (11) / 90.909090909090909090909090909090909D-3/
      DATA RFACT (11) / 25.052108385441718775052108385441719D-9/
      DATA RSFACT(11) / 0.158278578416163817828398057493079D-3/
      DATA GHALF (11) / 6.71353055419921875D+6/
      DATA FACT  (12) / 0.4790016D+9/
      DATA POW2  (12) / 4.096D+3/
      DATA RINT  (12) / 83.333333333333333333333333333333333D-3/
      DATA RFACT (12) / 2.087675698786809897921009032120143D-9/
      DATA RSFACT(12) / 45.691089927761735304439466708873392D-6/
      DATA GHALF (12) / 77.205601373291015625D+6/
      DATA FACT  (13) / 6.2270208D+9/
      DATA POW2  (13) / 8.192D+3/
      DATA RINT  (13) / 76.923076923076923076923076923076923D-3/
      DATA RFACT (13) / 0.160590438368216145993923771701549D-9/
      DATA RSFACT(13) / 12.672428274337012241581782357088822D-6/
      DATA GHALF (13) / 0.9650700171661376953125D+9/
      DATA FACT  (14) / 87.1782912D+9/
      DATA POW2  (14) / 16.384D+3/
      DATA RINT  (14) / 71.428571428571428571428571428571429D-3/
      DATA RFACT (14) / 11.470745597729724713851697978682106D-12/
      DATA RSFACT(14) / 3.38684891864543091458842524506959D-6/
      DATA GHALF (14) / 13.02844523174285888671875D+9/
      DATA FACT  (15) / 1.307674368D+12/
      DATA POW2  (15) / 32.768D+3/
      DATA RINT  (15) / 66.666666666666666666666666666666667D-3/
      DATA RFACT (15) / 0.764716373181981647590113198578807D-12/
      DATA RSFACT(15) / 0.874480630535623497631514529580712D-6/
      DATA GHALF (15) / 0.188912455860271453857421875D+12/
      DATA FACT  (16) / 20.922789888D+12/
      DATA POW2  (16) / 65.536D+3/
      DATA RINT  (16) / 62.5D-3/
      DATA RFACT (16) / 47.79477332387385297438207491117544D-15/
      DATA RSFACT(16) / 0.218620157633905874407878632395178D-6/
      DATA GHALF (16) / 2.9281430658342075347900390625D+12/
      DATA FACT  (17) / 0.355687428096D+15/
      DATA POW2  (17) / 0.131072D+6/
      DATA RINT  (17) / 58.823529411764705882352941176470588D-3/
      DATA RFACT (17) / 2.81145725434552076319894558301032D-15/
      DATA RSFACT(17) / 53.023176577281002838698207596591649D-9/
      DATA GHALF (17) / 48.31436058626442432403564453125D+12/
      DATA FACT  (18) / 6.402373705728D+15/
      DATA POW2  (18) / 0.262144D+6/
      DATA RINT  (18) / 55.555555555555555555555555555555556D-3/
      DATA RFACT (18) / 0.156192069685862264622163643500573D-15/
      DATA RSFACT(18) / 12.497682572615703325361452975928255D-9/
      DATA GHALF (18) / 0.845501310259627425670623779296875D+15/
      DATA FACT  (19) / 0.121645100408832D+18/
      DATA POW2  (19) / 0.524288D+6/
      DATA RINT  (19) / 52.631578947368421052631578947368421D-3/
      DATA RFACT (19) / 8.220635246624329716955981236872281D-18/
      DATA RSFACT(19) / 2.867165019077961915862419911964499D-9/
      DATA GHALF (19) / 15.641774239803107374906539916992187D+15/
      DATA FACT  (20) / 2.43290200817664D+18/
      DATA POW2  (20) / 1.048576D+6/
      DATA RINT  (20) / 50.D-3/
      DATA RFACT (20) / 0.411031762331216485847799061843614D-18/
      DATA RSFACT(20) / 0.641117588536780424092550271881354D-9/
      DATA GHALF (20) / 0.305014597676160593810677528381348D+18/
      DATA FACT  (21) / 51.09094217170944D+18/
      DATA POW2  (21) / 2.097152D+6/
      DATA RINT  (21) / 47.619047619047619047619047619047619D-3/
      DATA RFACT (21) / 19.57294106339126123084757437350543D-21/
      DATA RSFACT(21) / 0.139903327563683277929154381570029D-9/
      DATA GHALF (21) / 6.252799252361292173118889331817627D+18/
      DATA FACT  (22) / 1.12400072777760768D+21/
      DATA POW2  (22) / 4.194304D+6/
      DATA RINT  (22) / 45.454545454545454545454545454545455D-3/
      DATA RFACT (22) / 0.889679139245057328674889744250247D-21/
      DATA RSFACT(22) / 29.827489657110893526440716366250928D-12/
      DATA GHALF (22) / 0.134435183925767781722056120634079D+21/
      DATA FACT  (23) / 25.85201673888497664D+21/
      DATA POW2  (23) / 8.388608D+6/
      DATA RINT  (23) / 43.478260869565217391304347826086957D-3/
      DATA RFACT (23) / 38.681701706306840377169119315228123D-24/
      DATA RSFACT(23) / 6.219461528645935790146741966488423D-12/
      DATA GHALF (23) / 3.024791638329775088746262714266777D+21/
      DATA FACT  (24) / 0.62044840173323943936D+24/
      DATA POW2  (24) / 16.777216D+6/
      DATA RINT  (24) / 41.666666666666666666666666666666667D-3/
      DATA RFACT (24) / 1.611737571096118349048713304801172D-24/
      DATA RSFACT(24) / 1.269542268337733743314267195169946D-12/
      DATA GHALF (24) / 71.08260350074971458553717378526926D+21/
      DATA FACT  (25) / 15.511210043330985984D+24/
      DATA POW2  (25) / 33.554432D+6/
      DATA RINT  (25) / 40.D-3/
      DATA RFACT (25) / 64.469502843844733961948532192046872D-27/
      DATA RSFACT(25) / 0.253908453667546748662853439033989D-12/
      DATA GHALF (25) / 1.741523785768368007345660757739097D+24/
      DATA FACT  (26) / 0.403291461126605635584D+27/
      DATA POW2  (26) / 67.108864D+6/
      DATA RINT  (26) / 38.461538461538461538461538461538462D-3/
      DATA RFACT (26) / 2.479596263224797460074943545847957D-27/
      DATA RSFACT(26) / 49.795544612191937146986240757490316D-15/
      DATA GHALF (26) / 44.40885653709338418731434932234697D+24/
      DATA FACT  (27) / 10.888869450418352160768D+27/
      DATA POW2  (27) / 0.134217728D+9/
      DATA RINT  (27) / 37.037037037037037037037037037037037D-3/
      DATA RFACT (27) / 91.836898637955461484257168364739134D-30/
      DATA RSFACT(27) / 9.583157028764344602559973083500462D-15/
      DATA GHALF (27) / 1.176834698232974680963830257042195D+27/
      DATA FACT  (28) / 0.304888344611713860501504D+30/
      DATA POW2  (28) / 0.268435456D+9/
      DATA RINT  (28) / 35.714285714285714285714285714285714D-3/
      DATA RFACT (28) / 3.279889237069837910152041727312112D-30/
      DATA RSFACT(28) / 1.811046448070793658116887233276123D-15/
      DATA GHALF (28) / 32.362954201406803726505332068660355D+27/
      DATA FACT  (29) / 8.841761993739701954543616D+30/
      DATA POW2  (29) / 0.536870912D+9/
      DATA RINT  (29) / 34.482758620689655172413793103448276D-3/
      DATA RFACT (29) / 0.113099628864477169315587645769383D-30/
      DATA RSFACT(29) / 0.336302882628854624148194940917598D-15/
      DATA GHALF (29) / 0.92234419474009390620540196395682D+30/
      DATA FACT  (30) / 0.26525285981219105863630848D+33/
      DATA POW2  (30) / 1.073741824D+9/
      DATA RINT  (30) / 33.333333333333333333333333333333333D-3/
      DATA RFACT (30) / 3.769987628815905643852921525646106D-33/
      DATA RSFACT(30) / 61.40022498994531896057241577568715D-18/
      DATA GHALF (30) / 27.209153744832770233059357936726193D+30/
      DATA FACT  (31) / 8.22283865417792281772556288D+33/
      DATA POW2  (31) / 2.147483648D+9/
      DATA RINT  (31) / 32.258064516129032258064516129032258D-3/
      DATA RFACT (31) / 0.121612504155351794962997468569229D-33/
      DATA RSFACT(31) / 11.027805953831060903294206084777897D-18/
      DATA GHALF (31) / 0.829879189217399492108310417070149D+33/
      DATA FACT  (32) / 0.26313083693369353016721801216D+36/
      DATA POW2  (32) / 4.294967296D+9/
      DATA RINT  (32) / 31.25D-3/
      DATA RFACT (32) / 3.800390754854743592593670892788413D-36/
      DATA RSFACT(32) / 1.949459092890831498350848732335029D-18/
      DATA GHALF (32) / 26.14119446034808400141177813770969D+33/
      DATA FACT  (33) / 8.68331761881188649551819440128D+36/
      DATA POW2  (33) / 8.589934592D+9/
      DATA RINT  (33) / 30.303030303030303030303030303030303D-3/
      DATA RFACT (33) / 0.115163356207719502805868814932982D-36/
      DATA RSFACT(33) / 0.339357269271956958923882999874471D-18/
      DATA GHALF (33) / 0.849588819961312730045882789475565D+36/
      DATA FACT  (34) / 0.29523279903960414084761860964352D+39/
      DATA POW2  (34) / 17.179869184D+9/
      DATA RINT  (34) / 29.411764705882352941176470588235294D-3/
      DATA RFACT (34) / 3.387157535521161847231435733323006D-39/
      DATA RSFACT(34) / 58.199291537966008536289637900019698D-21/
      DATA GHALF (34) / 28.461225468703976456537073447431425D+36/
      DATA FACT  (35) / 10.3331479663861449296666513375232D+39/
      DATA POW2  (35) / 34.359738368D+9/
      DATA RINT  (35) / 28.571428571428571428571428571428571D-3/
      DATA RFACT (35) / 96.775929586318909920898163809228749D-42/
      DATA RSFACT(35) / 9.837475773099464661313797707198607D-21/
      DATA GHALF (35) / 0.981912278670287187750529033936384D+39/
      DATA FACT  (36) / 0.371993326789901217467999448150835D+42/
      DATA POW2  (36) / 68.719476736D+9/
      DATA RINT  (36) / 27.777777777777777777777777777777778D-3/
      DATA RFACT (36) / 2.688220266286636386691615661367465D-42/
      DATA RSFACT(36) / 1.639579295516577443552299617866434D-21/
      DATA GHALF (36) / 34.857885892795195165143780704741638D+39/
      DATA FACT  (37) / 13.763753091226345046315979581580902D+42/
      DATA POW2  (37) / 0.137438953472D+12/
      DATA RINT  (37) / 27.027027027027027027027027027027027D-3/
      DATA RFACT (37) / 72.65460179153071315382745030722879D-45/
      DATA RSFACT(37) / 0.269545175789756462029253873660073D-21/
      DATA GHALF (37) / 1.27231283508702462352774799572307D+42/
      DATA FACT  (38) / 0.523022617466601111760007224100074D+45/
      DATA POW2  (38) / 0.274877906944D+12/
      DATA RINT  (38) / 26.315789473684210526315789473684211D-3/
      DATA RFACT (38) / 1.911963205040281925100722376506021D-45/
      DATA RSFACT(38) / 43.726001475555501850577269073286858D-24/
      DATA GHALF (38) / 47.711731315763423382290549839615117D+42/
      DATA FACT  (39) / 20.397882081197443358640281739902897D+45/
      DATA POW2  (39) / 0.549755813888D+12/
      DATA RINT  (39) / 25.641025641025641025641025641025641D-3/
      DATA RFACT (39) / 49.024697565135433976941599397590277D-48/
      DATA RSFACT(39) / 7.001763889559218346960519325076389D-24/
      DATA GHALF (39) / 1.836901655656891800218186168825182D+45/
      DATA FACT  (40) / 0.815915283247897734345611269596116D+48/
      DATA POW2  (40) / 1.099511627776D+12/
      DATA RINT  (40) / 25.D-3/
      DATA RFACT (40) / 1.225617439128385849423539984939757D-48/
      DATA RSFACT(40) / 1.107076076486338787609666749137715D-24/
      DATA GHALF (40) / 72.557615398447226108618353668594689D+45/
      DATA FACT  (41) / 33.452526613163807108170062053440752D+48/
      DATA POW2  (41) / 2.199023255552D+12/
      DATA RINT  (41) / 24.39024390243902439024390243902439D-3/
      DATA RFACT (41) / 29.89310827142404510789121914487212D-51/
      DATA RSFACT(41) / 0.172896235561749710143568141603804D-24/
      DATA GHALF (41) / 2.938583423637112657399043323578085D+48/
      DATA FACT  (42) / 1.405006117752879898543142606244512D+51/
      DATA POW2  (42) / 4.398046511104D+12/
      DATA RINT  (42) / 23.80952380952380952380952380952381D-3/
      DATA RFACT (42) / 0.711740673129143931140267122496955D-51/
      DATA RSFACT(42) / 26.678468343012946123487786371058121D-27/
      DATA GHALF (42) / 0.121951212080940175282060297928491D+51/
      DATA FACT  (43) / 60.415263063373835637355132068513998D+51/
      DATA POW2  (43) / 8.796093022208D+12/
      DATA RINT  (43) / 23.255813953488372093023255813953488D-3/
      DATA RFACT (43) / 16.552108677421951886982956337138494D-54/
      DATA RSFACT(43) / 4.068428280973126836837834154540446D-27/
      DATA GHALF (43) / 5.182926513439957449487562661960847D+51/
      DATA FACT  (44) / 2.658271574788448768043625811014616D+54/
      DATA POW2  (44) / 17.592186044416D+12/
      DATA RINT  (44) / 22.727272727272727272727272727272727D-3/
      DATA RFACT (44) / 0.376184288123226179249612644025875D-54/
      DATA RSFACT(44) / 0.613338640657203481275610018818555D-27/
      DATA GHALF (44) / 0.225457303334638149052708975795297D+54/
      DATA FACT  (45) / 0.119622220865480194561963161495658D+57/
      DATA POW2  (45) / 35.184372088832D+12/
      DATA RINT  (45) / 22.222222222222222222222222222222222D-3/
      DATA RFACT (45) / 8.359650847182803983324725422797219D-57/
      DATA RSFACT(45) / 91.431126249121551457559310940651915D-30/
      DATA GHALF (45) / 10.03284999839139763284554942289071D+54/
      DATA FACT  (46) / 5.502622159812088949850305428800255D+57/
      DATA POW2  (46) / 70.368744177664D+12/
      DATA RINT  (46) / 21.739130434782608695652173913043478D-3/
      DATA RFACT (46) / 0.18173154015614791268097229179994D-57/
      DATA RSFACT(46) / 13.480784107615844193475778943030503D-30/
      DATA GHALF (46) / 0.456494674926808592294472498741527D+57/
      DATA FACT  (47) / 0.258623241511168180642964355153612D+60/
      DATA POW2  (47) / 0.140737488355328D+15/
      DATA RINT  (47) / 21.276595744680851063829787234042553D-3/
      DATA RFACT (47) / 3.866628513960593886829197697871054D-60/
      DATA RSFACT(47) / 1.966374459242337160505823677022503D-30/
      DATA GHALF (47) / 21.22700238409659954169297119148102D+57/
      DATA FACT  (48) / 12.413915592536072670862289047373375D+60/
      DATA POW2  (48) / 0.281474976710656D+15/
      DATA RINT  (48) / 20.833333333333333333333333333333333D-3/
      DATA RFACT (48) / 80.554760707512372642274952038980296D-63/
      DATA RSFACT(48) / 0.283821705842792038987956419783863D-30/
      DATA GHALF (48) / 1.008282613244588478230416131595348D+60/
      DATA FACT  (49) / 0.608281864034267560872252163321295D+63/
      DATA POW2  (49) / 0.562949953421312D+15/
      DATA RINT  (49) / 20.408163265306122448979591836734694D-3/
      DATA RFACT (49) / 1.643974708316579033515815347734292D-63/
      DATA RSFACT(49) / 40.545957977541719855422345683408946D-33/
      DATA GHALF (49) / 48.901706742362541194175182382374399D+60/
      DATA FACT  (50) / 30.414093201713378043612608166064769D+63/
      DATA POW2  (50) / 1.125899906842624D+15/
      DATA RINT  (50) / 20.D-3/
      DATA RFACT (50) / 32.879494166331580670316306954685835D-66/
      DATA RSFACT(50) / 5.734064367124908781068510871766589D-33/
      DATA GHALF (50) / 2.420634483746945789111671527927533D+63/
      DATA FACT  (51) / 1.551118753287382280224243016469303D+66/
      DATA POW2  (51) / 2.251799813685248D+15/
      DATA RINT  (51) / 19.607843137254901960784313725490196D-3/
      DATA RFACT (51) / 0.644695964045717268045417783425212D-66/
      DATA RSFACT(51) / 0.8029296133819684122794446712267D-33/
      DATA GHALF (51) / 0.12224204142922076235013941216034D+66/
      DATA FACT  (52) / 80.658175170943878571660636856403767D+66/
      DATA POW2  (52) / 4.503599627370496D+15/
      DATA RINT  (52) / 19.230769230769230769230769230769231D-3/
      DATA RFACT (52) / 12.397999308571485923950341988946393D-69/
      DATA RSFACT(52) / 0.111346303524506308255441020745856D-33/
      DATA GHALF (52) / 6.295465133604869261032179726257531D+66/
      DATA FACT  (53) / 4.2748832840600255642980137533894D+69/
      DATA POW2  (53) / 9.007199254740992D+15/
      DATA RINT  (53) / 18.867924528301886792452830188679245D-3/
      DATA RFACT (53) / 0.23392451525606577215000645262163D-69/
      DATA RSFACT(53) / 15.294591045728086813794122672902405D-36/
      DATA GHALF (53) / 0.33051191951425563620418943562852D+69/
      DATA FACT  (54) / 0.230843697339241380472092742683028D+72/
      DATA POW2  (54) / 18.014398509481984D+15/
      DATA RINT  (54) / 18.518518518518518518518518518518519D-3/
      DATA RFACT (54) / 4.331935467704921706481600974474631D-72/
      DATA RSFACT(54) / 2.081330215920799461087372271423983D-36/
      DATA GHALF (54) / 17.68238769401267653692413480612584D+69/
      DATA FACT  (55) / 12.696403353658275925965100847566517D+72/
      DATA POW2  (55) / 36.028797018963968D+15/
      DATA RINT  (55) / 18.181818181818181818181818181818182D-3/
      DATA RFACT (55) / 78.762463049180394663301835899538741D-75/
      DATA RSFACT(55) / 0.280646509062878590197158160257672D-36/
      DATA GHALF (55) / 0.963690129323690871262365346933858D+72/
      DATA FACT  (56) / 0.710998587804863451854045647463725D+75/
      DATA POW2  (56) / 72.057594037927936D+15/
      DATA RINT  (56) / 17.857142857142857142857142857142857D-3/
      DATA RFACT (56) / 1.406472554449649904701818498206049D-75/
      DATA RSFACT(56) / 37.502967275265698696993942093458853D-39/
      DATA GHALF (56) / 53.484802177464843355061276754829134D+72/
      DATA FACT  (57) / 40.526919504877216755680601905432322D+75/
      DATA POW2  (57) / 0.144115188075855872D+18/
      DATA RINT  (57) / 17.543859649122807017543859649122807D-3/
      DATA RFACT (57) / 24.67495709560789306494418417905349D-78/
      DATA RSFACT(57) / 4.967389364204087944423771598065769D-39/
      DATA GHALF (57) / 3.021891323026763649560962136647846D+75/
      DATA FACT  (58) / 2.350561331282878571829474910515075D+78/
      DATA POW2  (58) / 0.288230376151711744D+18/
      DATA RINT  (58) / 17.241379310344827586206896551724138D-3/
      DATA RFACT (58) / 0.425430294751860225257658347914715D-78/
      DATA RSFACT(58) / 0.652250178038963991862751285636931D-39/
      DATA GHALF (58) / 0.173758751074038909849755322857251D+78/
      DATA FACT  (59) / 0.138683118545689835737939019720389D+81/
      DATA POW2  (59) / 0.576460752303423488D+18/
      DATA RINT  (59) / 16.949152542372881355932203389830508D-3/
      DATA RFACT (59) / 7.210682961895936021316243184995175D-81/
      DATA RSFACT(59) / 84.915740365941201757196487020509921D-42/
      DATA GHALF (59) / 10.164886937831276226210686387149192D+78/
      DATA FACT  (60) / 8.320987112741390144276341183223364D+81/
      DATA POW2  (60) / 1.152921504606846976D+18/
      DATA RINT  (60) / 16.666666666666666666666666666666667D-3/
      DATA RFACT (60) / 0.120178049364932267021937386416586D-81/
      DATA RSFACT(60) / 10.96257494227210594634142227287644D-42/
      DATA GHALF (60) / 0.604810772800960935459535840035377D+81/
      DATA FACT  (61) / 0.507580213877224798800856812176625D+84/
      DATA POW2  (61) / 2.305843009213693952D+18/
      DATA RINT  (61) / 16.393442622950819672131147540983607D-3/
      DATA RFACT (61) / 1.970131956802168311835039121583381D-84/
      DATA RSFACT(61) / 1.403613891639067638874853396034929D-42/
      DATA GHALF (61) / 36.591051754458136595301918322140305D+81/
      DATA FACT  (62) / 31.469973260387937525653122354950764D+84/
      DATA POW2  (62) / 4.611686018427387904D+18/
      DATA RINT  (62) / 16.129032258064516129032258064516129D-3/
      DATA RFACT (62) / 31.776321883905940513468372928764214D-87/
      DATA RSFACT(62) / 0.178259142497393217190701951267799D-42/
      DATA GHALF (62) / 2.250349682899175400611067976811629D+84/
      DATA FACT  (63) / 1.982608315404440064116146708361898D+87/
      DATA POW2  (63) / 9.223372036854775808D+18/
      DATA RINT  (63) / 15.873015873015873015873015873015873D-3/
      DATA RFACT (63) / 0.504386061649300643070926554424829D-87/
      DATA RSFACT(63) / 22.45854095103465627338217714496675D-45/
      DATA GHALF (63) / 0.140646855181198462538191748550727D+87/
      DATA FACT  (64) / 0.126886932185884164103433389335161D+90/
      DATA POW2  (64) / 18.446744073709551616D+18/
      DATA RINT  (64) / 15.625D-3/
      DATA RFACT (64) / 7.88103221327032254798322741288795D-90/
      DATA RSFACT(64) / 2.807317618879332034172772143120844D-45/
      DATA GHALF (64) / 8.931075304006102371175176032971152D+87/
      DATA FACT  (65) / 8.247650592082470666723170306785496D+90/
      DATA POW2  (65) / 36.893488147419103232D+18/
      DATA RINT  (65) / 15.384615384615384615384615384615385D-3/
      DATA RFACT (65) / 0.121246649434928039199741960198276D-90/
      DATA RSFACT(65) / 0.348204895765306607211894245708544D-45/
      DATA GHALF (65) / 0.576054357108393602940798854126639D+90/
      DATA FACT  (66) / 0.544344939077443064003729240247843D+93/
      DATA POW2  (66) / 73.786976294838206464D+18/
      DATA RINT  (66) / 15.151515151515151515151515151515152D-3/
      DATA RFACT (66) / 1.837070445983758169693060003004184D-93/
      DATA RSFACT(66) / 42.861059786054732709835951410726255D-48/
      DATA GHALF (66) / 37.731560390599780992622324945294872D+90/
      DATA FACT  (67) / 36.471110918188685288249859096605464D+93/
      DATA POW2  (67) / 0.147573952589676412928D+21/
      DATA RINT  (67) / 14.925373134328358208955223880597015D-3/
      DATA RFACT (67) / 27.418961880354599547657611985137077D-96/
      DATA RSFACT(67) / 5.236311858584685132338020811704675D-48/
      DATA GHALF (67) / 2.509148765974885436009384608862109D+93/
      DATA FACT  (68) / 2.480035542436830599600990418569172D+96/
      DATA POW2  (68) / 0.295147905179352825856D+21/
      DATA RINT  (68) / 14.705882352941176470588235294117647D-3/
      DATA RFACT (68) / 0.40322002765227352275967076448731D-96/
      DATA RSFACT(68) / 0.634996084753499502047251416515707D-48/
      DATA GHALF (68) / 0.169367541703304766930633461098192D+96/
      DATA FACT  (69) / 0.171122452428141311372468338881273D+99/
      DATA POW2  (69) / 0.590295810358705651712D+21/
      DATA RINT  (69) / 14.492753623188405797101449275362319D-3/
      DATA RFACT (69) / 5.843768516699616271879286441845072D-99/
      DATA RSFACT(69) / 76.444545369173438965398109206527089D-51/
      DATA GHALF (69) / 11.601676606676376534748392085226177D+96/
      DATA FACT  (70) / 11.978571669969891796072783721689099D+99/
      DATA POW2  (70) / 1.180591620717411303424D+21/
      DATA RINT  (70) / 14.285714285714285714285714285714286D-3/
      DATA RFACT (70) / 83.482407381423089598275520597786739D-102/
      DATA RSFACT(70) / 9.136870765279713071611847472341708D-51/
      DATA GHALF (70) / 0.806316524164008169165013249923219D+99/
      DATA FACT  (71) / 0.850478588567862317521167644239926D+102/
      DATA POW2  (71) / 2.361183241434822606848D+21/
      DATA RINT  (71) / 14.08450704225352112676056338028169D-3/
      DATA RFACT (71) / 1.175808554667930839412331276025165D-102/
      DATA RSFACT(71) / 1.084347063752159952060748414892483D-51/
      DATA GHALF (71) / 56.845314953562575926133434119586959D+99/
      DATA FACT  (72) / 61.234458376886086861524070385274673D+102/
      DATA POW2  (72) / 4.722366482869645213696D+21/
      DATA RINT  (72) / 13.888888888888888888888888888888889D-3/
      DATA RFACT (72) / 16.330674370387928325171267722571741D-105/
      DATA RSFACT(72) / 0.127791526989812310972180020647395D-51/
      DATA GHALF (72) / 4.064440019179724178718540539550468D+102/
      DATA FACT  (73) / 4.470115461512684340891257138125051D+105/
      DATA POW2  (73) / 9.444732965739290427392D+21/
      DATA RINT  (73) / 13.698630136986301369863013698630137D-3/
      DATA RFACT (73) / 0.223707868087505867468099557843449D-105/
      DATA RSFACT(73) / 14.956866920832914446610781617092691D-54/
      DATA GHALF (73) / 0.294671901390530002957094189117409D+105/
      DATA FACT  (74) / 0.330788544151938641225953028221254D+108/
      DATA POW2  (74) / 18.889465931478580854784D+21/
      DATA RINT  (74) / 13.513513513513513513513513513513514D-3/
      DATA RFACT (74) / 3.023079298479809019839183214100655D-108/
      DATA RSFACT(74) / 1.738700462552365291639887940252999D-54/
      DATA GHALF (74) / 21.658384752203955217346422900129554D+105/
      DATA FACT  (75) / 24.809140811395398091946477116594034D+108/
      DATA POW2  (75) / 37.778931862957161709568D+21/
      DATA RINT  (75) / 13.333333333333333333333333333333333D-3/
      DATA RFACT (75) / 40.307723979730786931189109521342073D-111/
      DATA RSFACT(75) / 0.20076783601894698643569112640719D-54/
      DATA GHALF (75) / 1.613549664039194663692308506059652D+108/
      DATA FACT  (76) / 1.885494701666050254987932260861147D+111/
      DATA POW2  (76) / 75.557863725914323419136D+21/
      DATA RINT  (76) / 13.157894736842105263157894736842105D-3/
      DATA RFACT (76) / 0.530364789206984038568277756859764D-111/
      DATA RSFACT(76) / 23.029650218945663321654351444779763D-57/
      DATA GHALF (76) / 0.121822999634959197108769292207504D+111/
      DATA FACT  (77) / 0.145183092028285869634070784086308D+114/
      DATA POW2  (77) / 0.151115727451828646838272D+24/
      DATA RINT  (77) / 12.987012987012987012987012987012987D-3/
      DATA RFACT (77) / 6.887854405285506994393217621555378D-114/
      DATA RSFACT(77) / 2.624472214614875093894811233170713D-57/
      DATA GHALF (77) / 9.319459472074378578820850853874034D+111/
      DATA FACT  (78) / 11.324281178206297831457521158732046D+114/
      DATA POW2  (78) / 0.302231454903657293676544D+24/
      DATA RINT  (78) / 12.820512820512820512820512820512821D-3/
      DATA RFACT (78) / 88.305825708788551210169456686607412D-117/
      DATA RSFACT(78) / 0.297162961535902975342504844851203D-57/
      DATA GHALF (78) / 0.722258109085764339858615941175238D+114/
      DATA FACT  (79) / 0.894618213078297528685144171539832D+117/
      DATA POW2  (79) / 0.604462909807314587353088D+24/
      DATA RINT  (79) / 12.658227848101265822784810126582278D-3/
      DATA RFACT (79) / 1.117795262136563939369233628944398D-117/
      DATA RSFACT(79) / 33.433445262738985921693964019925236D-60/
      DATA GHALF (79) / 56.697261563232500678901351382256153D+114/
      DATA FACT  (80) / 71.569457046263802294811533723186532D+117/
      DATA POW2  (80) / 1.208925819614629174706176D+24/
      DATA RINT  (80) / 12.5D-3/
      DATA RFACT (80) / 13.97244077670704924211542036180497D-120/
      DATA RSFACT(80) / 3.737972816475134459336537294838897D-60/
      DATA GHALF (80) / 4.507432294276983803972657434889364D+117/
      DATA FACT  (81) / 5.797126020747367985879734231578109D+120/
      DATA POW2  (81) / 2.417851639229258349412352D+24/
      DATA RINT  (81) / 12.345679012345679012345679012345679D-3/
      DATA RFACT (81) / 0.172499268848235175828585436565493D-120/
      DATA RSFACT(81) / 0.415330312941681606592948588315433D-60/
      DATA GHALF (81) / 0.362848299689297196219798923508594D+120/
      DATA FACT  (82) / 0.475364333701284174842138206989405D+123/
      DATA POW2  (82) / 4.835703278458516698824704D+24/
      DATA RINT  (82) / 12.195121951219512195121951219512195D-3/
      DATA RFACT (82) / 2.103649620100428973519334592262115D-123/
      DATA RSFACT(82) / 45.865560283293487458390724499218288D-63/
      DATA GHALF (82) / 29.572136424677721491913612265950396D+120/
      DATA FACT  (83) / 39.455239697206586511897471180120611D+123/
      DATA POW2  (83) / 9.671406556917033397649408D+24/
      DATA RINT  (83) / 12.048192771084337349397590361445783D-3/
      DATA RFACT (83) / 25.345176145788300885775115569423077D-126/
      DATA RSFACT(83) / 5.034399283508242047248520314755466D-63/
      DATA GHALF (83) / 2.439701255035912023082873011940908D+123/
      DATA FACT  (84) / 3.314240134565353266999387579130131D+126/
      DATA POW2  (84) / 19.342813113834066795298816D+24/
      DATA RINT  (84) / 11.904761904761904761904761904761905D-3/
      DATA RFACT (84) / 0.301728287449860724830656137731227D-126/
      DATA RSFACT(84) / 0.549297995126380126159648736241153D-63/
      DATA GHALF (84) / 0.203715054795498653927419896497066D+126/
      DATA FACT  (85) / 0.281710411438055027694947944226061D+129/
      DATA POW2  (85) / 38.685626227668133590597632D+24/
      DATA RINT  (85) / 11.764705882352941176470588235294118D-3/
      DATA RFACT (85) / 3.54974455823365558624301338507326D-129/
      DATA RSFACT(85) / 59.579732780817803750498329963296379D-66/
      DATA GHALF (85) / 17.213922130219636256866981254002059D+126/
      DATA FACT  (86) / 24.22709538367273238176552320344126D+129/
      DATA POW2  (86) / 77.371252455336267181195264D+24/
      DATA RINT  (86) / 11.627906976744186046511627906976744D-3/
      DATA RFACT (86) / 41.276099514344832398174574245037907D-132/
      DATA RSFACT(86) / 6.42464781247539143254974870760248D-66/
      DATA GHALF (86) / 1.471790342133778899962126897217176D+129/
      DATA FACT  (87) / 2.10775729837952771721360051869939D+132/
      DATA POW2  (87) / 0.154742504910672534362390528D+27/
      DATA RINT  (87) / 11.494252873563218390804597701149425D-3/
      DATA RFACT (87) / 0.474437925452239452852581313161355D-132/
      DATA RSFACT(87) / 0.688794545167308546138812092037018D-66/
      DATA GHALF (87) / 0.127309864594571874846723976609286D+132/
      DATA FACT  (88) / 0.185482642257398439114796845645546D+135/
      DATA POW2  (88) / 0.309485009821345068724781056D+27/
      DATA RINT  (88) / 11.363636363636363636363636363636364D-3/
      DATA RFACT (88) / 5.391340061957266509688424013197219D-135/
      DATA RSFACT(88) / 73.425745225753524589527625144305752D-69/
      DATA GHALF (88) / 11.139613152025039049088347953312501D+132/
      DATA FACT  (89) / 16.507955160908461081216919262453619D+135/
      DATA POW2  (89) / 0.618970019642690137449562112D+27/
      DATA RINT  (89) / 11.235955056179775280898876404494382D-3/
      DATA RFACT (89) / 60.576854628733331569532854080867627D-138/
      DATA RSFACT(89) / 7.78311342771858436504425946554958D-69/
      DATA GHALF (89) / 0.985855763954215955844318793868156D+135/
      DATA FACT  (90) / 1.485715964481761497309522733620826D+138/
      DATA POW2  (90) / 1.237940039285380274899124224D+27/
      DATA RINT  (90) / 11.111111111111111111111111111111111D-3/
      DATA RFACT (90) / 0.673076162541481461883698378676307D-138/
      DATA RSFACT(90) / 0.820412190634367318141187348977996D-69/
      DATA GHALF (90) / 88.234090873902328048066532051199994D+135/
      DATA FACT  (91) / 0.135200152767840296255166568759495D+141/
      DATA POW2  (91) / 2.475880078570760549798248448D+27/
      DATA RINT  (91) / 10.989010989010989010989010989010989D-3/
      DATA RFACT (91) / 7.396441346609686394326355809629747D-141/
      DATA RSFACT(91) / 86.002565930381905055482162392258429D-72/
      DATA GHALF (91) / 7.985185224088160688350021150633599D+138/
      END
