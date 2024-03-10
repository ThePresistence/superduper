C     ******************************************************************
C
C     Computation of the dipole moment derivatives.
C
C     ******************************************************************
      SUBROUTINE PSDDPS(DIPDRV,LDF,ROALP,ROBET)
C
C   Compute static part of the dipole moment derivatives.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      DIPDRV - Input: should be zero. Output: Dipole moment derivatives.
C      LDF    - Leading dimension of the DIPDRV
C      ROALP  - Alpha density matrix
C      ROBET  - Beta density matrix (not used unless UHF is set)
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGBL - Global computation parameters
C      ATOMS  - Starting orbital number, number of orbitals
C               and atom types.
C      PARDER - Core charges
C      AMASS  - Atomic masses of the atoms in molecule.
C
C   Modified common blocks:
C
C   Local storage:
C
C   Module logic:
C
C      Static part of the dipole moment derivative is computed 
C      (assuming the center-of-mass definition for ions) as:
C
C                      Q m
C      d D                A     A
C     ----- = Core  - ------ - Sum P
C      d x        A     M       mu  mu,mu
C         A
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (TWO=2.0D0)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      SAVE /PSDGBL/
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./AMASS / AMS(LM1)
C
      DIMENSION DIPDRV(LDF,3), ROALP(*), ROBET(*)
C
C    If we have an ion, compute the Q/M ratio, so that we produce
C    the same (i.e., equally non-physical) output as Thiel's code do.
C
      IF( ICHRGE.NE.0 ) THEN
          AMASS = ZERO
          DO 100 I=1,NATOM
              AMASS = AMASS + AMS(I)
  100     CONTINUE
          QOVM = DBLE(ICHRGE)/AMASS
      ELSE
          QOVM = ZERO
      ENDIF
C
C    Static contributions are present only then the the nuclear
C    motion coincides with the dipole direction. All three non-zero
C    components are the same.
C
      II = 1
      DO 500 I=1,NATOM
C
C        Compute sum of the diagonal elements of the atomic block of the 
C        density matrix.
C
          QAT = ZERO
          IF(UHF) THEN
              DO 200 IORB=NFIRST(I),NLAST(I)
                  QAT = QAT + ROALP(II) + ROBET(II)
                  II  = II + IORB + 1
  200         CONTINUE
          ELSE
              DO 210 IORB=NFIRST(I),NLAST(I)
                  QAT = QAT + ROALP(II)
                  II  = II + IORB + 1
  210         CONTINUE
              QAT = TWO*QAT
          ENDIF
C
          DRV = TORE(NAT(I)) - QOVM*AMS(I) - QAT
          DIPDRV(3*(I-1)+1,1) = DRV
          DIPDRV(3*(I-1)+2,2) = DRV
          DIPDRV(3*(I-1)+3,3) = DRV
  500 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSDDPR(DIPD,DPA,DPB)
C
C   Compute response part of the dipole moment derivatives.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      DIPD   - Output: Response part of the moment derivatives for
C               the density matrix derivatives given on input.
C      DPA    - Derivative of the alpha density matrix.
C      DPA    - Derivative of the beta density matrix 
C               (not used unless UHF is set).
C
C   Accessed common blocks:
C
C      PSDOPT - Computation options. 
C      PSDGBL - Global computation parameters
C      ATOMS  - Starting orbital number, number of orbitals
C               and atom types.
C      ATOMC  - Coordinates of the atoms.
C      MULTIP - Multipole separation parameters (A.K.A.
C               dipole moment integrals ;-)
C
C   Modified common blocks:
C
C      PSDENS - Storage place for a density matrix elements derivatives.
C
C   Local storage:
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LM1, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      PARAMETER (SRSQ3=-0.577350269189625764509148780501957455D0)
      PARAMETER (TWORS3=1.154700538379251529018297561003914911D0)
      PARAMETER (RBOHR=1.0D0/0.529167D0)
C
      PARAMETER (ISPX=2)
      PARAMETER (ISPY=4)
      PARAMETER (ISPZ=7)
      PARAMETER (IZYZ=32)
      PARAMETER (IXXY=38)
      PARAMETER (IYZ2=24)
      PARAMETER (IYX2Y2=13)
      PARAMETER (IZZ2=25)
      PARAMETER (IXXZ=17)
      PARAMETER (IYYZ=31)
      PARAMETER (IZXZ=19)
      PARAMETER (IXZ2=23)
      PARAMETER (IXX2Y2=12)
      PARAMETER (IYXY=39)
C
      LOGICAL UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
      COMMON 
     ./PSDGBL/ MODE, NATOM, NPTCHG, NVARS, NORBS, NROS, METHOD, ICHRGE,
     .         UHF, HALFEL, DOCI, DORESP, DODIP, LIMAG, DOPTCH
     ./PSDENS/ PAAALP(45), PAABET(45), ZAA, PAA(45),
     .         PBBALP(45), PBBBET(45), ZBB, PBB(45),
     .         PABALP(9,9), PABBET(9,9), PAB(9,9)
      SAVE /PSDGBL/
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ATOMC / COORD(3,LM1)
     ./MULTIP/ DD(6,0:LMZ),PO(9,0:LMZ)
C
      DIMENSION DIPD(3), DPA(*), DPB(*)
      DIMENSION PAAX(136) 
      EQUIVALENCE (PAAX(1),PAAALP(1))
C
      DIPD(1) = ZERO
      DIPD(2) = ZERO
      DIPD(3) = ZERO
C
      DO 500 IATOM=1,NATOM
          NFRST = NFIRST(IATOM)
          NORB  = NLAST(IATOM) - NFRST + 1
          CALL PSDPAA(DPA,DPB,NFRST,NORB,PAAX)
C
C        Contributions from the atomic charges: contract the diagonal
C        elements of the density matrix derivative with the negated
C        coordinate of the corresponding atom.
C
          QATD = ZERO
          II   = 1
          DO 100 I=1,NORB
              QATD = QATD + PAA(II)
              II   = II + I + 1
  100     CONTINUE
          QATD = -RBOHR * QATD
          DIPD(1) = DIPD(1) + QATD*COORD(1,IATOM)
          DIPD(2) = DIPD(2) + QATD*COORD(2,IATOM)
          DIPD(3) = DIPD(3) + QATD*COORD(3,IATOM)
C
C        Derivatives of the hybrid contributions: collect density matrix
C        elements by multipole and contract with the corresponding
C        charge separation.
C
          IF( NORB.GT.1 ) THEN
C
C            There are SP hybrid contributions, these are easy:
C
              DDSP = -DD(2,NAT(IATOM))
              DIPD(1) = DIPD(1) + DDSP*PAA(ISPX)
              DIPD(2) = DIPD(2) + DDSP*PAA(ISPY)
              DIPD(3) = DIPD(3) + DDSP*PAA(ISPZ)
C
              IF( NORB.GT.4 ) THEN
C
C                PD hybrid contributions are present, these are slightly more tricky:
C
                  DDPD = -DD(5,NAT(IATOM))
                  QX   = PAA(IZXZ)+SRSQ3*PAA(IXZ2)+PAA(IXX2Y2)+PAA(IYXY)
                  QY   = PAA(IZYZ)+PAA(IXXY)+SRSQ3*PAA(IYZ2)-PAA(IYX2Y2)
                  QZ   = TWORS3*PAA(IZZ2)+PAA(IXXZ)+PAA(IYYZ)
                  DIPD(1) = DIPD(1) + DDPD*QX
                  DIPD(2) = DIPD(2) + DDPD*QY
                  DIPD(3) = DIPD(3) + DDPD*QZ
              ENDIF
          ENDIF
  500 CONTINUE
      RETURN
      END
