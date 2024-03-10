      SUBROUTINE MPCSYB (EIG,LM2,NB16,MUL)
C     *
C     WRITE GENERAL SYBYL OUTPUT DATA ON FILE NB16.
C     ADAPTED FROM MOPAC(6.0) WRITTEN BY J.J.P.STEWART.
C     *
      USE LIMIT, ONLY: LM1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CHARGE/ QI(LM1)
     ./DIPOL / DM,DMX(3),DMP(3)
     ./INOPT2/ IN2(300)
     ./NBFILE/ NBF(20)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
     ./SCFPRT/ EBIND,ENERGY,EIP,EIP2,SZ,S2,ITSAVE
      DIMENSION EIG(LM2)
C     INPUT OPTIONS.
      KHARGE = IN2(65)
C     CHARGE FLAG AND NUMBER OF ATOMS.
      WRITE(NB16,500,ERR=20) MUL,NUMAT
C     CARTESIAN COORDINATES AND ATOMIC CHARGES (STANDARD ZDO).
      DO 10 I=1,NUMAT
      WRITE(NB16,510,ERR=20) (COORD(J,I),J=1,3),QI(I)
   10 CONTINUE
C     TWO HIGHEST AND TWO LOWEST ORBITAL ENERGIES.
      I1     = MAX(1,NUMB-1)
      I2     = MIN(NORBS,NUMB+2)
      WRITE(NB16,520,ERR=20) (EIG(J),J=I1,I2),NUMB
C     HEAT OF FORMATION AND IONIZATION POTENTIAL.
      WRITE(NB16,530,ERR=20) ENERGY,EIP
C     MOLECULAR CHARGE AND DIPOLE MOMENT.
      WRITE(NB16,540,ERR=20) KHARGE,DM
      RETURN
C *** ERROR EXIT.
   20 CONTINUE
      NB6 = NBF(6)
      WRITE(NB6,550)
      RETURN
  500 FORMAT(2I4)
  510 FORMAT(4F12.6)
  520 FORMAT(4F12.6,2X,I4,2X,'HOMOs,LUMOs,# of occupied MOs')
  530 FORMAT(2F12.6,4X,'HF and IP')
  540 FORMAT(I4,F10.3,'  Charge,Dipole Moment')
  550 FORMAT(//1X,'Error writing SYBYL MOPAC output')
      END
