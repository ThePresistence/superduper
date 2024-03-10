C     ******************************************************************
C
C     NMR prescreening routines for 3-center integrals.
C
C     ******************************************************************
      LOGICAL FUNCTION PSNINC(IA,IB)
C
C   TRUE if three-center integrals involving atoms IA and IB should
C   be explicitly included in the computation.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IA,IB  - Atom indices.
C
C   Accessed common blocks:
C
C      PSNPRS - Results of the integral prescreening:
C          IATMX  - Largest atomic index prescreened
C          IPAIRS - Bit-packed array of atomic pairs
C      PSPRT  - Printing unit
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      Nope.
C
C   Module logic:
C
C   Bugs:
C
C      Integers are assumed to be big enough to hold INTBIT bits.
C
C      We are relying on BTEST intrinsic, which is not standard
C      in Fortran77 (but is standard in Fortran90).
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (INTBIT=32)
      PARAMETER (INTPRS=((LM1*(LM1+1))/2+INTBIT-1)/INTBIT)
      COMMON
     ./PSNPRS/ IATMX, IPAIRS(INTPRS)
     ./PSPRT / NB6
      SAVE /PSNPRS/, /PSPRT /
      LOGICAL BTEST
      INTRINSIC BTEST
C    Range checking...
      IF( IA.LE.0 .OR. IA.GT.IATMX .OR. IB.LE.0 .OR. IB.GT.IATMX ) THEN
          WRITE(NB6,11000) IA, IB, IATMX
          STOP 'PSNINC'
      ENDIF
C
      IF( IA.LE.IB ) THEN
          IBIT = (IB*(IB-1))/2 + IA
      ELSE
          IBIT = (IA*(IA-1))/2 + IB
      ENDIF
      ILOC = (IBIT-1)/INTBIT + 1
      IOFF = (IBIT-1) - (ILOC-1)*INTBIT
      PSNINC = BTEST(IPAIRS(ILOC),IOFF)
      RETURN
11000 FORMAT(' ATOM INDICES ',I5,' AND ',I5,' ARE NOT FROM 1 TO ',I5,
     .       ' IN PSNINC.')
      END
C
      DOUBLE PRECISION FUNCTION PSNDMX(RO,IBASAI,NORBAI,IBASBI,NORBBI)
C
C   Find the maximum of two-center block of the density matrix stored
C   in the packed form.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      RO     - Density matrix in the packed form
C      IBASAI - First orbital on the first center
C      NORBAI - Number of orbitals on the first center
C      IBASBI - First orbital on the second center
C      NORBBI - Number of orbitals on the second center
C
C   Accessed common blocks:
C
C      Nope.
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      Nope.
C
C   Module logic:
C
C      Lifted from PSDPAB.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RO(*)
C
      IF( IBASAI.GT.IBASBI ) THEN
          IBASA = IBASAI
          IBASB = IBASBI
          NORBA = NORBAI
          NORBB = NORBBI
      ELSE
          IBASA = IBASBI
          IBASB = IBASAI
          NORBA = NORBBI
          NORBB = NORBAI
      ENDIF
C
      IFIRST = (IBASA*(IBASA+1))/2-IBASA+IBASB
      IIN    = IFIRST-1
      PSNDMX = 0.D0
      DO 200 I=1,NORBA
          DO 100 J=1,NORBB
              PSNDMX = MAX(PSNDMX,ABS(RO(IIN+J)))
  100     CONTINUE
          IIN = IIN + I + IBASA-1
  200 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSNMNA(ITYP,NMAX,AMIN)
C
C   Determine largest N and smallest ALP for the atom type.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ITYP   - Atom type
C      MMAX   - (Output) Largest major quantum number
C      AMIN   - (Output) Smallest orbital exponent
C
C   Accessed common blocks:
C
C      PSNPAR - Major quantum numbers and orbital exponents.
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      Nope.
C
C   Module logic:
C
C   Bugs:
C
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXL=2)
      PARAMETER (MAXC=3)
      COMMON
     ./PSNPAR/ BETAN(0:MAXL,LMZ), ZN(MAXC,0:MAXL,LMZ), 
     .         CN(MAXC,0:MAXL,LMZ), SCALN(0:LMZ),
     .         NN(MAXC,0:MAXL,LMZ), MAXN(0:MAXL,LMZ), MAXLA(LMZ),
     .         IBETA
      SAVE /PSNPAR/
C
      NMAX = NN(1,0,ITYP)
      AMIN = ZN(1,0,ITYP)
      DO 400 L=0,MAXLA(ITYP)
          DO 300 I=1,MAXN(L,ITYP)
              NMAX = MAX(NMAX,NN(I,L,ITYP))
              AMIN = MIN(AMIN,ZN(I,L,ITYP))
  300     CONTINUE
  400 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PSNSCR(RO,RO1)
C
C   Pre-screen three-center NMR integrals using zero- and first-order
C   density matrices.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      RO     - Zero-order density matrix.
C      ROA    - Array of the three first-order density matrices.
C               Not used if DORESP flag is not set.
C
C   Accessed common blocks:
C
C      PSDOPT - We need IPRINT from here
C      PSDGBL - NROS and DORESP come from here
C      PSPRT  - Printing unit
C      ATOMS  - Atom types and orbital indices
C      ATOMC  - Coordinates (in Angstroms)
C
C   Modified common blocks:
C
C      PSNPRS - Results of the integral prescreening:
C          IATMX  - Largest atomic index prescreened
C          IPAIRS - Bit-packed array of atomic pairs
C
C   Local storage:
C
C      Nope.
C
C   Module logic:
C
C      All atom pairs which have product of the S-type overlap (computed
C      with the largest major quantum number and smallest orbital exponent 
C      among any of the orbitals) and largest of the density matrix elements
C      in excess of DNCOFF are marked as needed in the computation.
C
C   Bugs:
C
C      Integers are assumed to be big enough to hold INTBIT bits.
C
C      We are relying on BTEST and IBSET intrinsics, which are not standard
C      in Fortran77 (but is standard in Fortran90).
C
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (BOHR=0.529167D0)
      PARAMETER (INTBIT=32)
      PARAMETER (INTPRS=((LM1*(LM1+1))/2+INTBIT-1)/INTBIT)
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
     ./PSNPRS/ IATMX, IPAIRS(INTPRS)
     ./PSPRT / NB6
      SAVE /PSNPRS/, /PSDGBL/, /PSDOPT/, /PSPRT /
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ATOMC / COORD(3,LM1)
      LOGICAL BTEST
      INTRINSIC IBSET, BTEST
      DIMENSION RO(NROS), RO1(NROS,3)
      DIMENSION OV(1,1)
C    Range and sanity checking...
      DO 100 I=0,INTBIT-1
          IDUMMY = IBSET(0,I)
          IF( BTEST(0,I) .OR. .NOT.BTEST(IDUMMY,I) ) THEN
              WRITE(NB6,11010) I
              STOP 'PSNSCR'
          ENDIF
  100 CONTINUE
      IF( NATOM.GT.LM1 ) THEN
          WRITE(NB6,11020) NATOM, LM1
          STOP 'PSNSCR'
      ENDIF
C    Zero out pairs array
      DO 200 I=1,((NATOM*(NATOM+1))/2+INTBIT-1)/INTBIT
          IPAIRS(I) = 0
  200 CONTINUE
      IATMX = NATOM
C
      INEED = 0
      DO 4000 IA=2,NATOM
          IBASA = NFIRST(IA)
          NORBA = NLAST(IA)-NFIRST(IA)+1
          CALL PSNMNA(NAT(IA),NAMAX,ZAMIN)
          DO 3900 IB=1,IA-1
              IBASB = NFIRST(IB)
              NORBB = NFIRST(IB)-NFIRST(IB)+1
              ROMAX = PSNDMX(RO,IBASA,NORBA,IBASB,NORBB)
              IF(DORESP) THEN
                  DO 1000 IRO=1,3
                      ROMAX = MAX(ROMAX,PSNDMX(RO1(1,IRO),IBASA,NORBA,
     .                                                    IBASB,NORBB))
 1000             CONTINUE
              ENDIF
C            Overlap integral won't be larger than 1, anyway
              IF( ROMAX.LT.DNCOFF ) GOTO 3899
C            Can't decide from density alone, examine overlap integral 
C            for this atom pair
              RX = (COORD(1,IB) - COORD(1,IA))/BOHR
              RY = (COORD(2,IB) - COORD(2,IA))/BOHR
              RZ = (COORD(3,IB) - COORD(3,IA))/BOHR
              CALL PSNMNA(NAT(IB),NBMAX,ZBMIN)
              CALL PSOVX(NAMAX,0,ZAMIN,NBMAX,0,ZBMIN,RX,RY,RZ,OV,1)
              IF( ABS(OV(1,1))*ROMAX.LT.DNCOFF ) GOTO 3899
C            Mark this pair as necessary.
              INEED = INEED + 1
              IBIT = (IA*(IA-1))/2 + IB
              ILOC = (IBIT-1)/INTBIT + 1
              IOFF = (IBIT-1) - (ILOC-1)*INTBIT
              IPAIRS(ILOC) = IBSET(IPAIRS(ILOC),IOFF)

 3899         CONTINUE
 3900     CONTINUE
 4000 CONTINUE
      IF( IPRINT.GE.1 ) WRITE(NB6,10050)INEED,(NATOM*(NATOM-1))/2,DNCOFF
      RETURN
10050 FORMAT(' THREE-CENTER TERMS INCLUDED FOR ',I9,
     .       ' ATOM PAIRS (OUT OF ',I9,'). CUTOFF = ',G10.4)
11010 FORMAT(' CANNOT MANIPULATE BIT ',I3,' IN AN INTEGER WORD. ',
     .       'PSNSCR WILL WORK INCORRECTLY.')
11020 FORMAT(' NUMBER OF ATOMS (',I4,') IS ABOVE THE LEGAL MAXIMUM (',
     .       I4,') IN PSNSCR.')
      END
