C    ******************************************************************
C
C     Routines for computing response quantities in the AO basis.
C     Analogous to Fock matrices.
C
C     ******************************************************************
      SUBROUTINE PSRESP(YA,YB,Y,AINTS,TA,TB,LDC,DOBETA)
C
C   Generate response matrices from alpha and beta generalized
C   density matrices.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      YA     - Alpha part of the generalized density matrix,
C               held in an upper triangle of the square array.
C               Destroyed on return.
C      YB     - Beta part of the generalized density matrix,
C               not used unless UHF is set.
C      Y      - Scratch space for the total generalized density.
C      AINTS  - Precomputed two-electron integrals.
C      TA     - Alpha response matrix (filled on return)
C      TB     - Beta response matrix (filled on return)
C               Not used unless UHF and DOBETA are both set.
C      LDC    - Leading dimension of YA, YB, Y, TA and TB
C               matrices
C      DOBETA - .TRUE. if beta response matrix is to be 
C               computed. This flag is ignored unless UHF is
C               also set.
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C      PSPRT  - Printing unit
C      PSPRTF - Debug output control flags
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
C     Scaling code could be much more efficient.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (HALF = 0.5D0)
      PARAMETER (ONE  = 1.0D0)
      PARAMETER (SONE =-1.0D0)
      PARAMETER (TWO  = 2.0D0)
      PARAMETER (FOUR = 4.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DOBETA, DODIP, LIMAG, DOPTCH
      LOGICAL LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM
      LOGICAL LPNUME, LPA2MO, LPSOLV, LPNMR, LPCI
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSPRT / NB6
     ./PSPRTF/ LPHALF, LPCPHF, LPPASS, LPMMGR, LPINTS, LPCPHK, LPSUMM,
     .         LPNUME, LPA2MO, LPSOLV, LPNMR,  LPCI
      SAVE /PSDGBL/, /PSPRTF/, /PSPRT /
C
      DIMENSION YA(LDC,*), YB(LDC,*), Y(LDC,*), TA(LDC,*), TB(LDC,*)
      DIMENSION AINTS(*)
C
C     Rescale generalized density matrices and compute total density.
C     Total density is not needed for imaginary perturbations
C     since only exchange part survives after transformation into MO
C     basis in this case.
C
      IF( .NOT.LIMAG ) THEN
          IF(UHF) THEN
              DO 100 I=1,NORBS
                  CALL DCOPY(I,     YA(1,I),1,Y(1,I),1)
                  CALL DAXPY(I,ONE, YB(1,I),1,Y(1,I),1)
                  CALL DSCAL(I,TWO, Y (1,I),1)
  100         CONTINUE
          ELSE
              DO 150 I=1,NORBS
                  CALL DCOPY(I,     YA(1,I),1,Y(1,I),1)
                  CALL DSCAL(I,FOUR,Y (1,I),1)
  150         CONTINUE
          ENDIF
          CALL DSCAL(NORBS,HALF,Y,LDC+1)
      ELSE
C
C        We have to change sign (i.e., effectively take complex
C        conjugate) of the density matrix to compensate for the
C        incorrect expression implemented in PSFOCA
C
          DO 200 I=1,NORBS
              CALL DSCAL(I,SONE,YA(1,I),1)
  200     CONTINUE
      ENDIF
C
      IF(LPSUMM) THEN
          IF( .NOT.LIMAG ) THEN
              WRITE(NB6,11010) 
              CALL PSDPSU(NORBS,YA,LDC)
              IF(UHF) THEN
                  WRITE(NB6,11015) 
                  CALL PSDPSU(NORBS,YB,LDC)
              ENDIF
              WRITE(NB6,11020)
              CALL PSDPSU(NORBS,Y,LDC)
          ELSE
              WRITE(NB6,11010) 
              CALL PSDPAU(NORBS,YA,LDC)
          ENDIF
      ENDIF
C
C    Now, compute alpha and beta "response" matrices. This
C    is the most complicated part, but the least significant
C    one speed-wise: both previous and following transformations
C    have scaling power one higher the this one.
C
      IF( .NOT.LIMAG ) THEN
          CALL PSFOCK(YA,YB,Y,AINTS,TA,TB,LDC,DOBETA)
      ELSE
          CALL PSFOCA(YA,AINTS,TA,LDC)
      ENDIF
      IF(LPSUMM) THEN
          IF( .NOT.LIMAG ) THEN
              WRITE(NB6,11100)
              CALL PSDPSU(NORBS,TA,LDC)
              IF(UHF .AND. DOBETA) THEN
                  WRITE(NB6,11105)
                  CALL PSDPSU(NORBS,TB,LDC)
              ENDIF
          ELSE
              WRITE(NB6,11100)
              CALL PSDPAU(NORBS,TA,LDC)
          ENDIF
      ENDIF
C
      RETURN
11010 FORMAT(' YA IS:'/)
11015 FORMAT(' YB IS:'/)
11020 FORMAT(' Y IS:'/)
11100 FORMAT(' ALPHA "RESPONSE" MATRIX IS:'/)
11105 FORMAT(' BETA "RESPONSE" MATRIX IS:'/)
      END
C
      SUBROUTINE PSFOCK(YALP,YBET,Y,AI,TALP,TBET,LDX,DOBETA)
C
C   Generate "response" matrices from generalised density
C   matrices and stored two-electron integrals
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      YALP   - Generalized alpha density matrix, stored
C               as upper triangle of square matrix
C      YBET   - Generalized beta density matrix, not used
C               if UHF flag not set
C      Y      - Generalized total density matrix, with
C               diagonal elements halved. Not required for
C               imaginary perturbations
C      AI     - Stored two-center integrals
C      TALP   - Place for alpha response matrix.
C      TBET   - Place for beta response matrix, not used if
C               UHF flag not set or DOBETA is .FALSE.
C      LDX    - Leading dimension of YALP, YBET, Y, TALP and
C               TBET
C      DOBETA - .TRUE. if beta response matrix is to be 
C               computed. This flag is ignored unless UHF is
C               also set.
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      180 DOUBLE PRECISION cells for accumulation of contributions
C      to the "response" matrix and storing sections of "density"
C
C   Module logic:
C
C      Loop over atoms cut'n'paste-ed from PST1ST, see it for
C      comments, too.
C
C      Contributions to the diagonal block of orbitals at atom
C      with higher number are accumulated in the local array
C      instead of going directly to "response" matrix. This should
C      improve cache use. Apart from it, code is straightforward.
C
C      For imaginary density matrix only exchange
C      contributions are computed, since Coulomb part cancels out
C      after transformation into MO basis.
C
C   Speedups possible:
C
C      Replacing calls to DGEMM, PSFADX and PSFADC with inlined
C      versions tuned to the pairs of (NORBA,NORBB) might give
C      significant speedup.
C
C   Bugs:
C
C      This is formally equivalent to computing Fock matrices 
C      from density matrices then hamiltonian matrix is zero, 
C      so that probably code from the main program could be
C      hacked to conform to this case. Nonevertheless, I prefer 
C      to have my code as self-contained as possible, so that 
C      I'll give my own version of it. It should be of no 
C      serious consequences, because this function should not
C      give a significant contributuion to the computation time.
C      Besides, it is actually *faster* than the version supplied
C      with main program ;-)
C
C      Calling PSFADC for diagonal block computes NORBA*(NORBA-1)/2 
C      values unnecessarily, because they are already known by
C      symmetry
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ONE =1.0D0)
      PARAMETER (ZERO=0.0D0)
C
C*$*INLINE ROUTINE (PSZRV,DCOPY,DGEMV,PSFADX,PSFADC,PSU2UL)
CVD$R EXPAND(PSZRV,DCOPY,DGEMV,PSFADX,PSFADC,PSU2UL)
CVD$R SEARCH(psmutl.f,dblas1.f,dblas2.f)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DOBETA, BETA, DODIP, LIMAG
      LOGICAL DOPTCH
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      SAVE /PSDGBL/
C
      DIMENSION YALP(LDX,*),YBET(LDX,*),Y(LDX,*)
      DIMENSION TALP(LDX,*),TBET(LDX,*),AI(*)
C
      DIMENSION YAA(45), YBB(45), TAA(45), TBB(45)
C
      BETA = UHF .AND. DOBETA
      IAI2T1 = 1
      CALL PSZRM(NORBS,NORBS,TALP,LDX)
      IF(BETA) THEN
          CALL PSZRM(NORBS,NORBS,TBET,LDX)
      ENDIF
C
      DO 2000 NA=1,NATOM
          NFRSTA = NFIRST(NA)
          NORBA  = NLAST(NA) - NFIRST(NA) + 1
          NPAIRA = ( NORBA * (NORBA+1) ) / 2
          IND    = 1
          DO 100 J=0,NORBA-1
              DO 98 I=0,J
                 YAA(IND+I) = Y(NFRSTA+I,NFRSTA+J)
                 TAA(IND+I) = ZERO
   98         CONTINUE
              IND = IND + J + 1
  100     CONTINUE
          DO 1500 NB=1,NA-1
              NFRSTB = NFIRST(NB)
              NORBB  = NLAST(NB) - NFIRST(NB) + 1
              NPAIRB = ( NORBB * (NORBB+1) ) / 2
              IND    = 1
              DO 200 J=0,NORBB-1
                  DO 198 I=0,J
                     YBB(IND+I) = Y(NFRSTB+I,NFRSTB+J)
  198             CONTINUE
                  IND = IND + J + 1
  200         CONTINUE
C
C             General off-diagonal block.
C
              CALL DGEMV('N',NPAIRA,NPAIRB,ONE,AI(IAI2T1),NPAIRA,
     .                       YBB,1,ONE,TAA,1)
              CALL PSZRVB(NPAIRB,TBB)
              CALL DGEMV('T',NPAIRA,NPAIRB,ONE,AI(IAI2T1),NPAIRA,
     .                       YAA,1,ZERO,TBB,1)
              CALL PSFADX(TALP(NFRSTB,NFRSTB),LDX,NORBB,TBB)
              CALL PSFADC(TALP(NFRSTB,NFRSTA),LDX,NORBB,NORBA,
     .                       AI(IAI2T1),NPAIRA,YALP(NFRSTB,NFRSTA))
              IF(BETA) THEN
                  CALL PSFADX(TBET(NFRSTB,NFRSTB),LDX,NORBB,TBB)
                  CALL PSFADC(TBET(NFRSTB,NFRSTA),LDX,NORBB,NORBA,
     .                       AI(IAI2T1),NPAIRA,YBET(NFRSTB,NFRSTA))
              ENDIF
C
              IAI2T1 = IAI2T1 + NPAIRA*NPAIRB
 1500     CONTINUE
C
C         Diagonal blocks are special, because we have to
C         complement Y matrix to the square atom block first
C
          CALL DGEMV('N',NPAIRA,NPAIRA,ONE,AI(IAI2T1),NPAIRA,
     .                   YAA,1,ONE,TAA,1)
          CALL PSFADX(TALP(NFRSTA,NFRSTA),LDX,NORBA,TAA)
          CALL PSU2UL(YALP(NFRSTA,NFRSTA),LDX,NORBA)
          CALL PSFADC(TALP(NFRSTA,NFRSTA),LDX,NORBA,NORBA,
     .                   AI(IAI2T1),NPAIRA,YALP(NFRSTA,NFRSTA))
          IF(BETA) THEN
              CALL PSFADX(TBET(NFRSTA,NFRSTA),LDX,NORBA,TAA)
              CALL PSU2UL(YBET(NFRSTA,NFRSTA),LDX,NORBA)
              CALL PSFADC(TBET(NFRSTA,NFRSTA),LDX,NORBA,NORBA,
     .                       AI(IAI2T1),NPAIRA,YBET(NFRSTA,NFRSTA))
          ENDIF
C
          IAI2T1 = IAI2T1 + NPAIRA*NPAIRA
 2000 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSFOCA(YALP,AI,TALP,LDX)
C
C   Generate exchange part of response matrices from generalised 
C   density matrices and stored two-electron integrals for 
C   antisymmetric density (i.e., for the imaginary part of the 
C   density).
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      YALP   - Generalized alpha density matrix, stored
C               as upper triangle of square matrix
C      AI     - Stored two-center integrals
C      TALP   - Place for alpha response matrix.
C      LDX    - Leading dimension of YALP, YBET, Y, TALP and
C               TBET
C
C   Accessed common blocks:
C
C      PSDGBL - Global computation parameters
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
C      This is a specialized version of PSFOCK.
C
C   Bugs:
C
C      This function actually computes *complex conjugate* of the
C      Fock matrix. This is compensated for in the calling routine
C      by changing sign of the density matrix before calling PSFOCA.
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C*$*INLINE ROUTINE (PSFADC,PSU2AS)
CVD$R EXPAND(PSFADC,PSU2AS)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      SAVE /PSDGBL/
C
      DIMENSION YALP(LDX,*), TALP(LDX,*),AI(*)
C
      IAI2T1 = 1
      CALL PSZRM(NORBS,NORBS,TALP,LDX)
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
C             General off-diagonal block.
C
              CALL PSFADC(TALP(NFRSTB,NFRSTA),LDX,NORBB,NORBA,
     .                       AI(IAI2T1),NPAIRA,YALP(NFRSTB,NFRSTA))
C
              IAI2T1 = IAI2T1 + NPAIRA*NPAIRB
 1500     CONTINUE
C
C         Diagonal blocks are special, because we have to
C         complement Y matrix to the square atom block first
C
          CALL PSU2AS(YALP(NFRSTA,NFRSTA),LDX,NORBA)
          CALL PSFADC(TALP(NFRSTA,NFRSTA),LDX,NORBA,NORBA,
     .                   AI(IAI2T1),NPAIRA,YALP(NFRSTA,NFRSTA))
C
          IAI2T1 = IAI2T1 + NPAIRA*NPAIRA
 2000 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSFADX(T,LDX,NORB,TXX)
C
C   Add contribution to the upper diagonal block of 
C   the "response" matrix
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      T      - "Response" matrix positioned over atom-atom block
C               to update
C      LDX    - Leading dimension of T matrix
C      NORB   - Number of orbitals in block
C      TXX    - Contribution stored as linear array
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(LDX,*), TXX(*)
C
      IND = 0
      DO 100 J=1,NORB
          DO 98 I=1,J
              T(I,J) = T(I,J) + TXX(IND+I)
   98     CONTINUE
          IND = IND + J
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE PSFADC(T,LDX,NORBB,NORBA,AI,LDA,Y)
C
C   Compute and add contribution to the (usually) off-diagonal block
C   of "response" matrix.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      T      - "Response" matrix positioned over B-A block
C               to update
C      LDX    - Leading dimension of T and Y matrices
C      NORBA  - Number of orbitals on atom A
C      NORBB  - Number of orbitals on atom B
C      AI     - Two-electron integrals block A-B
C      LDA    - Leading dimension of AI
C      Y      - Generalised density matrix positioned over B-A block
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(LDX,*), AI(LDA,*), Y(LDX,*)
      INTEGER MU,NU,LA,SI
C
      INDA = 1
      DO 500 NU=1,NORBA
          DO 400 MU=1,NU-1
              INDB = 1
              DO 300 SI=1,NORBB
                  DO 200 LA=1,SI-1
                      AINT     = AI(INDA,INDB)
                      T(LA,MU) = T(LA,MU) - Y(SI,NU) * AINT
                      T(SI,MU) = T(SI,MU) - Y(LA,NU) * AINT
                      T(LA,NU) = T(LA,NU) - Y(SI,MU) * AINT
                      T(SI,NU) = T(SI,NU) - Y(LA,MU) * AINT
                      INDB = INDB + 1
  200             CONTINUE
C
C                 LA = SI
C
                  AINT     = AI(INDA,INDB)
                  T(SI,MU) = T(SI,MU) - Y(SI,NU) * AINT
                  T(SI,NU) = T(SI,NU) - Y(SI,MU) * AINT
                  INDB = INDB + 1
  300         CONTINUE
              INDA = INDA + 1
  400     CONTINUE
C
C         MU = NU
C
          INDB = 1
          DO 430 SI=1,NORBB
              DO 420 LA=1,SI-1
                  AINT     = AI(INDA,INDB)
                  T(LA,NU) = T(LA,NU) - Y(SI,NU) * AINT
                  T(SI,NU) = T(SI,NU) - Y(LA,NU) * AINT
                  INDB = INDB + 1
  420         CONTINUE
C
C             MU = NU, LA = SI
C
              AINT     = AI(INDA,INDB)
              T(SI,NU) = T(SI,NU) - Y(SI,NU) * AINT
              INDB = INDB + 1
  430     CONTINUE
          INDA = INDA + 1
  500 CONTINUE
C
      RETURN
      END
