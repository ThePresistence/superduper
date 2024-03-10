C     ******************************************************************
C
C     Compute overlap integrals and their derivatives
C     in the molecular coordinate system.
C
C     ******************************************************************
      SUBROUTINE PSOVER(IMOD,IA,IB,WX)
C
C   Compute overlap integrals along with first and second 
C   derivatives in molecular coordinate system
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD   - Computation type
C           0 = Compute integrals
C           1 = Compute integrals and first derivatives
C           2 = Compute integrals with 1st and 2nd 
C               derivatives
C      IA     - Number of first atom
C      IB     - Number of second atom
C      WX     - Output integrals. Third index enumerates:
C               integrals (=1), gradients (=2,3,4), and
C               second derivatives (=5,6,7,8,9,10).
C
C   Accessed common blocks:
C
C      ATOMS  - Indices of first/last orbital on atom,
C               atom types.
C      ATOMC  - Cartesian coordinates of atom
C      DNBND  - Major quantum numbers.
C      PSROT0,
C      PSROT1,
C      PSROT2,
C      PSDROT - Rotation matrices, their derivatives
C               and derivatives projection coefficients.
C               Those four are accessed indirectly by slave
C               routines.
C      PSDIST - Interatomic distance, a.u.
C
C   Modified common blocks:
C
C      PSOVLC - Used to return transformation coefficients
C               from PSOVC. PSOVLC is reported as volatile by
C               ftnchek. This is fine, as it should never be
C               used by any routines higher than PSOVER in the
C               calling tree. This would have been a good 
C               candidate for module data, if Fortran-77 had it.
C
C   Local storage:
C
C      333 DOUBLE PRECISION cells (some of these indirectly
C      through worker routines). This could be reduced
C      by not storing S2 (we don't need it afterwards)
C      and by re-using space occupied by Onn arrays for
C      TEMP array. I don't think it's worth the effort,
C      see memory requirements of psints.f for comparison.
C
C   Module logic:
C
C          First, overlap integrals are computed. If
C      first derivatives are desired, overlap integrals
C      for major quantum numbers one higher are computed
C      on each atom. If second derivatives are necessary,
C      overlap integrals for major quantum numbers two
C      higher on each of the atoms and one higher in both
C      are also computed. 
C          Then, first and second derivatives of overlap 
C      integrals are computed according to well-knows 
C      transformation formulas.
C          Then, derivatives of overlap integrals are projected
C      into molecular coordinate system
C          Finally, integrals are rotated into molecular
C      coordinate system using one-stage transformation
C      with combined rotation matrices.
C
C      See separate documentation for detailed timings.

C   Speedups possible:
C
C       Rewrite overlap computation module to compute
C       all necessary overlap integrals at once
C
C   Bugs:
C
C      S and P orbitals on the same center MUST have the
C      same major quantum numbers. There is no restriction 
C      on D orbitals.
C
C      The following common blocks used by PSOVER are 
C      supposed to retain values from previous call to 
C      PSINTS with the same atom types: 
C
C         PSROT0, PSROT1, PSROT2, PSDROT and PSDIST
C
C      They should be explicitly SAVE'd in all modules
C      using them to be conformant to F77 standard
C
C      First and second derivatives could have been
C      computed faster by writing explicit formulas for
C      them with following CSE a-la done for integrals.
C      Still, I do not feel it's worth the effort, see
C      timing results above along with timing data from
C      psints.f
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ISS=1,ISZ=2,IZS=3,IZZ=4,IXX=5,IZ2S=6,ISZ2=7,IZ2Z=8)
      PARAMETER (IZZ2=9,IXZX=10,IXXZ=11,IZ2Z2=12,IXZXZ=13,IXYXY=14)
C
      COMMON 
     ./PSDIST/ DIST
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./DNBND / III(LMZ),IIID(LMZ)
     ./PSOVLC/ C00, C10, C01, DS1, D10, D01, D11, D20, D02
      SAVE /PSDIST/
C
      DIMENSION WX(9,9,*)
C
C   Local temporary storage
C
      DIMENSION O00(14), O10(14), O01(14), O11(14), O20(14), O02(14)
      DIMENSION SL(14,0:2)
      EQUIVALENCE (SL,O00)
C
      CS1(N)   = C00 * O00(N) + C10 * O10(N) + C01 * O01(N)
      CS2(N)   = DS1 * SL(N,1) + D10 * O10(N) + D01 * O01(N) +
     .           D11 * O11(N) + D20 * O20(N) + D02 * O02(N)
C
      IORBA = NLAST(IA) - NFIRST(IA)
      IORBB = NLAST(IB) - NFIRST(IB)
      ITYPA = NAT(IA)
      ITYPB = NAT(IB)
      NSPA  = III (ITYPA) 
      NDA   = IIID(ITYPA)
      NSPB  = III (ITYPB) 
      NDB   = IIID(ITYPB)
C
C   Compute necessary overlap integrals.
C
      CALL OVERLP( ITYPA, ITYPB, DIST, O00, 0, 0 )
      IF( IMOD.GE.1 ) THEN
          CALL OVERLP( ITYPA, ITYPB, DIST, O10, 1, 0 )
          CALL OVERLP( ITYPA, ITYPB, DIST, O01, 0, 1 )
          IF( IMOD.GE.2 ) THEN 
              CALL OVERLP( ITYPA, ITYPB, DIST, O11, 1, 1 )
              CALL OVERLP( ITYPA, ITYPB, DIST, O20, 2, 0 )
              CALL OVERLP( ITYPA, ITYPB, DIST, O02, 0, 2 )
          ENDIF
      ENDIF
C
C   Transform intermediate integrals into derivatives and 
C   project them to the molecular coordinate system. This
C   step is only necessary if we do have derivatives...
C
      IF( IMOD.GE.1 ) THEN
C        SP-SP overlaps
          CALL PSOVC( NSPA, NSPB, DIST, IMOD )
          SL(ISS,1) = CS1(ISS)
          IF( IORBA.GT.0 ) THEN
              SL(IZS,1) = CS1(IZS)
          ENDIF
          IF( IORBB.GT.0 ) THEN
              SL(ISZ,1) = CS1(ISZ)
          ENDIF
          IF( IORBA.GT.0 .AND. IORBB.GT.0 ) THEN
              SL(IZZ,1) = CS1(IZZ)
              SL(IXX,1) = CS1(IXX)
          ENDIF
          IF( IMOD.GE.2 ) THEN
              SL(ISS,2) = CS2(ISS)
              IF( IORBA.GT.0 ) THEN
                  SL(IZS,2) = CS2(IZS)
              ENDIF
              IF( IORBB.GT.0 ) THEN
                  SL(ISZ,2) = CS2(ISZ)
              ENDIF
              IF( IORBA.GT.0 .AND. IORBB.GT.0 ) THEN
                  SL(IZZ,2) = CS2(IZZ)
                  SL(IXX,2) = CS2(IXX)
              ENDIF
          ENDIF
          IF( IORBA.GT.3 ) THEN
C            D-SP overlaps
              IF( NDA.NE.NSPA ) CALL PSOVC( NDA, NSPB, DIST, IMOD )
              SL(IZ2S,1) = CS1(IZ2S)
              IF( IORBB.GT.0 ) THEN
                  SL(IZ2Z,1) = CS1(IZ2Z)
                  SL(IXZX,1) = CS1(IXZX)
              ENDIF
              IF( IMOD.GE.2 ) THEN
                  SL(IZ2S,2) = CS2(IZ2S)
                  IF( IORBB.GT.0 ) THEN
                      SL(IZ2Z,2) = CS2(IZ2Z)
                      SL(IXZX,2) = CS2(IXZX)
                  ENDIF
              ENDIF
              IF( IORBB.GT.3 ) THEN
C                D-D  overlaps
                  IF( NDB.NE.NSPB ) CALL PSOVC( NDA, NDB, DIST, IMOD )
                  SL(IZ2Z2,1) = CS1(IZ2Z2)
                  SL(IXZXZ,1) = CS1(IXZXZ)
                  SL(IXYXY,1) = CS1(IXYXY)
                  IF( IMOD.GE.2 ) THEN
                      SL(IZ2Z2,2) = CS2(IZ2Z2)
                      SL(IXZXZ,2) = CS2(IXZXZ)
                      SL(IXYXY,2) = CS2(IXYXY)
                  ENDIF
              ENDIF
          ENDIF
          IF( IORBB.GT.3 ) THEN
C            SP-D overlaps
              IF( NSPA.NE.NDA .OR. NSPB.NE.NDB ) 
     .                  CALL PSOVC( NSPA, NDB, DIST, IMOD )
              SL(ISZ2,1) = CS1(ISZ2)
              IF( IORBA.GT.0 ) THEN
                  SL(IZZ2,1) = CS1(IZZ2)
                  SL(IXXZ,1) = CS1(IXXZ)
              ENDIF
              IF( IMOD.GE.2 ) THEN
                  SL(ISZ2,2) = CS2(ISZ2)
                  IF( IORBA.GT.0 ) THEN
                      SL(IZZ2,2) = CS2(IZZ2)
                      SL(IXXZ,2) = CS2(IXXZ)
                  ENDIF
              ENDIF
          ENDIF
      ENDIF
C    Rotate integrals and derivatives to correct orientation
      CALL PSOVRT(IMOD,IORBA+1,IORBB+1,SL,WX)
      RETURN
      END
C
      SUBROUTINE PSOVRT(IMOD,IORBA,IORBB,SL,WX)
C
C   Rotate overlap (or resonance) integrals and derivatives from
C   the local to the molecular coordinate system.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IMOD   - Computation type
C           0 = Compute integrals
C           1 = Compute integrals and first derivatives
C           2 = Compute integrals with 1st and 2nd 
C               derivatives
C      IORBA  - Number of orbitals on the left-hand atom
C      IORBB  - Number of orbitals on the rught-hand atom
C      SL     - Non-zero integrals and their radial gradients in the
C               local coordinate system.
C      WX     - (Output) integrals. The last index enumerates integrals
C               themselves (=1) and their first (=2,3,4) and second
C               (=5,6,7,8,9,10) Cartesian derivatives.
C
C   Accessed common blocks:
C
C      PSROT0,PSROT1,PSROT2, and PSDROT contain various rotation
C               matrices and their derivatives
C
C   Modified common blocks:
C
C      Nope.
C
C   Local storage:
C
C      A little.
C
C   Module logic:
C
C      Overlap-like integrals and their radial derivatives are first
C      projected into molecular coordinate system, and then rotated
C      using chain rules. See PSIDR, PSRM, and PSGVML for comments.
C
C   Speedups possible:
C
C      Inlining subroutine calls and cloning loops is the only hope.
C      It's going to be messy, though...
C
C   Bugs:
C
C       The code is limited to D orbitals at most.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ISS=1,ISZ=2,IZS=3,IZZ=4,IXX=5,IZ2S=6,ISZ2=7,IZ2Z=8)
      PARAMETER (IZZ2=9,IXZX=10,IXXZ=11,IZ2Z2=12,IXZXZ=13,IXYXY=14)
C
      COMMON 
     ./PSROT0/ RM (3,3), RMD (5,5), RMPP(6,6), RMPD(15,15), RMDD(15,15)
     ./PSROT1/ RM1(3,3,3), RMD1(5,5,3), RMPP1(6,6,3), RMPD1(15,15,3), 
     .         RMDD1(15,15,3)
     ./PSROT2/ RM2(3,3,6), RMD2(5,5,6), RMPP2(6,6,6), RMPD2(15,15,6), 
     .         RMDD2(15,15,6)
     ./PSDROT/ D1(3), D2A(6), D2B(6)
      SAVE /PSROT0/, /PSROT1/, /PSROT2/, /PSDROT/
C    Parameters
      DIMENSION SL(14,0:2)
      DIMENSION WX(9,9,*)
C    Local temporaries
      DIMENSION P1(14,3), P2(14,6)
      DIMENSION TEMP(9,9)
      LOGICAL DORB
C    Functions for projecting radial gradients into molecular system,
C    see PSIDR for definitions.
      CP1(I,J) = D1 (J) * SL(I,1)
      CP2(I,J) = D2A(J) * SL(I,2) + D2B(J) * SL(I,1)
C    Project first derivatives onto molecular coordinate system axes.
C    The loops below are designed to be concise. Any optimizer worth
C    its salt should be able to move IF's outside the loop if it is
C    likely to improve performance.
      IF( IMOD.GE.1 ) THEN
          DO 1000 I=1,3
C            S terms on the right-hand side center
              P1(ISS,I) = CP1(ISS,I)
              IF( IORBA.GT.1 ) THEN
                  P1(IZS, I) = CP1(IZS,I)
                  IF( IORBA.GT.4 ) P1(IZ2S,I) = CP1(IZ2S,I)
              ENDIF
C            P terms on the right-hand center
              IF( IORBB.GT.1 ) THEN
                  P1(ISZ,I) = CP1(ISZ,I)
                  IF( IORBA.GT.1 ) THEN
                      P1(IZZ,I) = CP1(IZZ,I)
                      P1(IXX,I) = CP1(IXX,I)
                      IF( IORBA.GT.4 ) THEN
                          P1(IZ2Z,I) = CP1(IZ2Z,I)
                          P1(IXZX,I) = CP1(IXZX,I)
                      ENDIF
                  ENDIF
C                D terms on the right-hand center
                  IF( IORBB.GT.4 ) THEN
                      P1(ISZ2,I) = CP1(ISZ2,I)
                      IF( IORBA.GT.1 ) THEN
                          P1(IZZ2,I) = CP1(IZZ2,I)
                          P1(IXXZ,I) = CP1(IXXZ,I)
                          IF( IORBA.GT.4 ) THEN
                              P1(IZ2Z2,I) = CP1(IZ2Z2,I)
                              P1(IXZXZ,I) = CP1(IXZXZ,I)
                              P1(IXYXY,I) = CP1(IXYXY,I)
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
 1000     CONTINUE
C        Project second derivatives
          IF( IMOD.GE.2 ) THEN
              DO 2000 I=1,6
C                S terms on the right-hand center
                  P2(ISS,I) = CP2(ISS,I)
                  IF( IORBA.GT.1 ) THEN
                      P2(IZS,I) = CP2(IZS,I)
                      IF( IORBA.GT.4 ) P2(IZ2S,I) = CP2(IZ2S,I)
                  ENDIF
C                P terms on the right-hand center
                  IF( IORBB.GT.1 ) THEN
                      P2(ISZ,I) = CP2(ISZ,I)
                      IF( IORBA.GT.1 ) THEN
                          P2(IZZ,I) = CP2(IZZ,I)
                          P2(IXX,I) = CP2(IXX,I)
                          IF( IORBA.GT.4 ) THEN
                              P2(IZ2Z,I) = CP2(IZ2Z,I)
                              P2(IXZX,I) = CP2(IXZX,I)
                          ENDIF
                      ENDIF
C                    D terms on the right-hand center
                      IF( IORBB.GT.4 ) THEN
                          P2(ISZ2,I) = CP2(ISZ2,I)
                          IF( IORBA.GT.1 ) THEN
                              P2(IZZ2,I) = CP2(IZZ2,I)
                              P2(IXXZ,I) = CP2(IXXZ,I)
                              IF( IORBA.GT.4 ) THEN
                                  P2(IZ2Z2,I) = CP2(IZ2Z2,I)
                                  P2(IXZXZ,I) = CP2(IXZXZ,I)
                                  P2(IXYXY,I) = CP2(IXYXY,I)
                              ENDIF
                          ENDIF
                      ENDIF
                  ENDIF
 2000         CONTINUE
          ENDIF
      ENDIF
C   Rotate integrals and derivatives to correct orientation
C   First, copy SS integrals which need no further transformation
      WX(1,1,1) = SL(ISS,0)
      IF( IMOD.GE.1 ) THEN
          DO 3190 I=1,3
              WX(1,1,1+I) = P1(ISS,I)
 3190     CONTINUE
          IF( IMOD.GE.2 ) THEN
              DO 3200 I=1,6
                  WX(1,1,4+I) = P2(ISS,I)
 3200         CONTINUE
          ENDIF
      ENDIF
C   Now, bite the bullet and unroll yet another copy of PSGVML...
      IF( IORBA.GT.1 .OR. IORBB.GT.1 ) THEN
          IA = IORBA - 1
          IB = IORBB - 1
          CALL PSOVRO(IA,IB,SL(1,0),WX,RM,RMPP,RMD,RMPD,RMDD)
          IF( IMOD.GT.0 ) THEN
C            First derivatives
              DO 4210 I=1,3
                  CALL PSOVRO(IA,IB,SL(1,0),TEMP,RM1(1,1,I),
     .                        RMPP1(1,1,I),RMD1(1,1,I), 
     .                        RMPD1(1,1,I),RMDD1(1,1,I))
                  CALL PSOVRO(IA,IB,P1(1,I),WX(1,1,1+I),RM,
     .                        RMPP,RMD,RMPD,RMDD)
                  CALL PSOVAD(IA,IB,TEMP,WX(1,1,1+I))
 4210         CONTINUE
C            Second derivatives
              IF( IMOD.GT.1 ) THEN
C                Fxy*G and F*Gxy contributions
                  DO 4220 I=1,6
                      CALL PSOVRO(IA,IB,SL(1,0),TEMP,RM2(1,1,I),
     .                            RMPP2(1,1,I),RMD2(1,1,I), 
     .                            RMPD2(1,1,I),RMDD2(1,1,I))
                      CALL PSOVRO(IA,IB,P2(1,I),WX(1,1,4+I),RM,
     .                            RMPP,RMD,RMPD,RMDD)
                      CALL PSOVAD(IA,IB,TEMP,WX(1,1,4+I))
 4220             CONTINUE
C                Fx*Gx contributions - diagonal second derivatives
                  DO 4230 I=1,3
                      CALL PSOVRO(IA,IB,P1(1,I),TEMP,RM1(1,1,I),
     .                            RMPP1(1,1,I),RMD1(1,1,I),
     .                            RMPD1(1,1,I),RMDD1(1,1,I))
                      CALL PSOVA2(IA,IB,TEMP,WX(1,1,4+I))
 4230             CONTINUE
C               Fx*Gy and Fy*Gx contributions - off-diagonal second derivatives
C               I1,I2 pair successfully assumes values (1,2), (1,3), (2,3)
                  DO 4240 I=1,3
                      I1 = MAX(I-1,1)
                      I2 = MIN(I+1,3)
C                    Fx*Gy contribution
                      CALL PSOVRO(IA,IB,P1(1,I1),TEMP,
     .                            RM1(1,1,I2),RMPP1(1,1,I2), 
     .                            RMD1(1,1,I2), RMPD1(1,1,I2),
     .                            RMDD1(1,1,I2))
                      CALL PSOVAD(IA,IB,TEMP,WX(1,1,4+I+3))
C                    Fy*Gx contribution
                      CALL PSOVRO(IA,IB,P1(1,I2),TEMP,
     .                            RM1(1,1,I1),RMPP1(1,1,I1), 
     .                            RMD1(1,1,I1), RMPD1(1,1,I1),
     .                            RMDD1(1,1,I1))
                      CALL PSOVAD(IA,IB,TEMP,WX(1,1,4+I+3))
 4240             CONTINUE
              ENDIF
          ENDIF
C        Make copies of symmetric portions of PP and DD overlaps
          IF( IORBA.GT.1 .AND. IORBB.GT.1 ) THEN
              IF( IORBA.GT.4 .AND. IORBB.GT.4 ) THEN
                  DORB = .TRUE.
              ELSE
                  DORB = .FALSE.
              ENDIF
              CALL PSOVCO(DORB,WX(1,1,1))
              IF( IMOD.GE.1 ) THEN
                  DO 5250 I=1,3
                      CALL PSOVCO(DORB,WX(1,1,1+I))
 5250             CONTINUE
                  IF( IMOD.GE.2 ) THEN
                      DO 5260 I=1,6
                          CALL PSOVCO(DORB,WX(1,1,4+I))
 5260                 CONTINUE
                  ENDIF
              ENDIF
          ENDIF
      ENDIF
      RETURN
      END
C
      SUBROUTINE PSOVC( NA, NB, DIST, IMOD )
C
C   Compute transformation coefficients for assembling
C   overlap integrals derivatives from overlap integrals
C   of higher major quantum number.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./PSOVLC/ C00, C10, C01, DS1, D10, D01, D11, D20, D02
C
      COEFF(N)     = SQRT( ( N + 1 ) * ( N + 0.5D0 ) )
C
      C00 = ( NA + NB + 1 ) / DIST
      C10 = - COEFF( NA )   / DIST
      C01 = - COEFF( NB )   / DIST
      IF( IMOD.LE.1 ) RETURN
C
      DS1  = ( NA + NB ) / DIST
      D10 = ( NA + NB + 2 ) * C10 / DIST
      D01 = ( NA + NB + 2 ) * C01 / DIST
      D11 = 2.0D0 * C10 * C01
      D20 = -( C10 * COEFF( NA + 1 ) ) / DIST
      D02 = -( C01 * COEFF( NB + 1 ) ) / DIST
      RETURN
      END
C
      SUBROUTINE PSOVRO(IORBA,IORBB,OVL,OUT,RM,RMPP,RMD,RMPD,RMDD)
C
C   "Rotate" single component of overlap matrix/derivative
C   RMxx could be either rotation matrices themselves or
C   corresponding cartesian derivatives. SS overlap is not
C   processed here because it have (implicit) unit rotation
C   matrix and is completely determined by projection of
C   local derivatives.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ISZ   =  2)
      PARAMETER (IZS   =  3)
      PARAMETER (IZZ   =  4)
      PARAMETER (IXX   =  5)
      PARAMETER (IZ2S  =  6)
      PARAMETER (ISZ2  =  7)
      PARAMETER (IZ2Z  =  8)
      PARAMETER (IZZ2  =  9)
      PARAMETER (IXZX  = 10)
      PARAMETER (IXXZ  = 11)
      PARAMETER (IZ2Z2 = 12)
      PARAMETER (IXZXZ = 13)
      PARAMETER (IXYXY = 14)
C
      DIMENSION OVL(14), OUT(9,9)
      DIMENSION RM(3,3), RMD(5,5), RMPP(6,6)
      DIMENSION RMPD(15,15), RMDD(15,15)
C
      CPP(I) = (RMPP(1,I)+RMPP(3,I))*OVL(IXX)+RMPP(6,I)*OVL(IZZ)
      CDP(I) = (RMPD(4,I)+RMPD(11,I))*OVL(IXZX)+RMPD(9,I)*OVL(IZ2Z)
      CPD(I) = (RMPD(4,I)+RMPD(11,I))*OVL(IXXZ)+RMPD(9,I)*OVL(IZZ2)
      CDD(I) = (RMDD(1,I)+RMDD(15,I))*OVL(IXYXY)+
     .         (RMDD(3,I)+RMDD(10,I))*OVL(IXZXZ)+RMDD(6,I)*OVL(IZ2Z2)
C
C   S orbital on B, PD orbitals on A
C
      IF( IORBA.GT.0 ) THEN
          DO 10 I=2,4
              OUT(I,1) = OVL(IZS) * RM(3,I-1)
   10     CONTINUE
          IF( IORBA.GT.3 ) THEN
              DO 20 I=5,9
                  OUT(I,1) = OVL(IZ2S) * RMD(3,I-4)
   20         CONTINUE
          ENDIF
      ENDIF
      IF( IORBB.LE.0 ) RETURN
C
C   P orbital on B, SPD orbitals on A
C
      DO 30 I=2,4
          OUT(1,I) = OVL(ISZ) * RM(3,I-1)
   30 CONTINUE
      IF( IORBA.GT.0 ) THEN
          OUT(2,2) = CPP(1)
          OUT(3,2) = CPP(2)
          OUT(3,3) = CPP(3)
          OUT(4,2) = CPP(4)
          OUT(4,3) = CPP(5)
          OUT(4,4) = CPP(6)
          IF( IORBA.GT.3 ) THEN
              DO 40 I=2,4
                  OUT(5,I) = CDP(I- 1)
                  OUT(6,I) = CDP(I+ 2)
                  OUT(7,I) = CDP(I+ 5)
                  OUT(8,I) = CDP(I+ 8)
                  OUT(9,I) = CDP(I+11)
   40         CONTINUE
          ENDIF
      ENDIF
      IF( IORBB.LE.3 ) RETURN
C
C   D orbital on B, SPD orbitals on A
C
      DO 50 I=5,9
          OUT(1,I) = OVL(ISZ2) * RMD(3,I-4)
   50 CONTINUE
      IF( IORBA.GT.0 ) THEN
          DO 60 I=2,4
              OUT(I,5) = CPD(I- 1)
              OUT(I,6) = CPD(I+ 2)
              OUT(I,7) = CPD(I+ 5)
              OUT(I,8) = CPD(I+ 8)
              OUT(I,9) = CPD(I+11)
   60     CONTINUE
          IF( IORBA.GT.3 ) THEN
              OUT(5,5) = CDD( 1)
              OUT(6,5) = CDD( 2)
              OUT(6,6) = CDD( 3)
              OUT(7,5) = CDD( 4)
              OUT(7,6) = CDD( 5)
              OUT(7,7) = CDD( 6)
              OUT(8,5) = CDD( 7)
              OUT(8,6) = CDD( 8)
              OUT(8,7) = CDD( 9)
              OUT(8,8) = CDD(10)
              OUT(9,5) = CDD(11)
              OUT(9,6) = CDD(12)
              OUT(9,7) = CDD(13)
              OUT(9,8) = CDD(14)
              OUT(9,9) = CDD(15)
          ENDIF
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE PSOVAD(IORBA,IORBB,FROM,TO)
C
C   Add symmetry-unique portion of overlap matrix
C   It is quite probable that straight additon of
C   rectangular matrix would be faster on vector
C   processors, although it would process 3 elements
C   unnecessary for PP pair (13 for DD pair)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FROM(9,9), TO(9,9)
C
      IF( IORBA.GT.0 ) THEN
C
C         P-S block
C
          DO 80 I=2,IORBA+1
             TO(I,1) = TO(I,1) + FROM(I,1)
   80     CONTINUE
          IF( IORBA.GT.3 .AND. IORBB.GT.0 ) THEN
C
C             D-P block
C
              DO 100 J=2,4
                  DO 90 I=5,9
                      TO(I,J) = TO(I,J) + FROM(I,J)
   90             CONTINUE
  100         CONTINUE
          ENDIF
      ENDIF
      IF( IORBB.GT.0 ) THEN
C
C         S-P block
C
          DO 110 I=2,IORBB+1
             TO(1,I) = TO(1,I) + FROM(1,I)
  110     CONTINUE
          IF( IORBB.GT.3 .AND. IORBA.GT.0 ) THEN
C
C             P-D block
C
              DO 130 J=5,9
                  DO 120 I=2,4
                      TO(I,J) = TO(I,J) + FROM(I,J)
  120             CONTINUE
  130         CONTINUE
          ENDIF
      ENDIF
      IF( IORBB.GT.0 .AND. IORBA.GT.0 ) THEN
C
C         P-P block
C
          DO 150 J=2,4
              DO 140 I=J,4
                  TO(I,J) = TO(I,J) + FROM(I,J)
  140         CONTINUE
  150     CONTINUE
          IF( IORBA.GT.3 .AND. IORBB.GT.3 ) THEN
C
C             D-D block
C
              DO 170 J=5,9
                  DO 160 I=J,9
                      TO(I,J) = TO(I,J) + FROM(I,J)
  160             CONTINUE
  170         CONTINUE
          ENDIF
      ENDIF
      RETURN
      END
C
      SUBROUTINE PSOVA2(IORBA,IORBB,FROM,TO)
C
C   Add doubled symmetry-unique portion of overlap matrix.
C   See notes to PSOVAD.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FROM(9,9), TO(9,9)
      PARAMETER (TWO=2.0D0)
C
      IF( IORBA.GT.0 ) THEN
C
C         P-S block
C
          DO 80 I=2,IORBA+1
             TO(I,1) = TO(I,1) + TWO * FROM(I,1)
   80     CONTINUE
          IF( IORBA.GT.3 .AND. IORBB.GT.0 ) THEN
C
C             D-P block
C
              DO 100 J=2,4
                  DO 90 I=5,9
                      TO(I,J) = TO(I,J) + TWO * FROM(I,J)
   90             CONTINUE
  100         CONTINUE
          ENDIF
      ENDIF
      IF( IORBB.GT.0 ) THEN
C
C         S-P block
C
          DO 110 I=2,IORBB+1
             TO(1,I) = TO(1,I) + TWO * FROM(1,I)
  110     CONTINUE
          IF( IORBB.GT.3 .AND. IORBA.GT.0 ) THEN
C
C             P-D block
C
              DO 130 J=5,9
                  DO 120 I=2,4
                      TO(I,J) = TO(I,J) + TWO * FROM(I,J)
  120             CONTINUE
  130         CONTINUE
          ENDIF
      ENDIF
      IF( IORBB.GT.0 .AND. IORBA.GT.0 ) THEN
C
C         P-P block
C
          DO 150 J=2,4
              DO 140 I=J,4
                  TO(I,J) = TO(I,J) + TWO * FROM(I,J)
  140         CONTINUE
  150     CONTINUE
          IF( IORBA.GT.3 .AND. IORBB.GT.3 ) THEN
C
C             D-D block
C
              DO 170 J=5,9
                  DO 160 I=J,9
                      TO(I,J) = TO(I,J) + TWO * FROM(I,J)
  160             CONTINUE
  170         CONTINUE
          ENDIF
      ENDIF

      RETURN
      END
C
      SUBROUTINE PSOVCO(DORB,OVL)
C
C   Replicate portions of overlap matrix by symmetry
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION OVL(9,9)
      LOGICAL DORB
C
      OVL(2,3) = OVL(3,2)
      OVL(2,4) = OVL(4,2)
      OVL(3,4) = OVL(4,3)
      IF( DORB ) THEN
          OVL(5,6) = OVL(6,5)
          OVL(5,7) = OVL(7,5)
          OVL(5,8) = OVL(8,5)
          OVL(5,9) = OVL(9,5)
          OVL(6,7) = OVL(7,6)
          OVL(6,8) = OVL(8,6)
          OVL(6,9) = OVL(9,6)
          OVL(7,8) = OVL(8,7)
          OVL(7,9) = OVL(9,7)
          OVL(8,9) = OVL(9,8)
      ENDIF
C
      RETURN
      END
