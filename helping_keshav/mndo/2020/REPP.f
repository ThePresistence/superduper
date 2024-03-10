      SUBROUTINE REPP (NI,NJ,R,A,CORE)
C     *
C     CALCULATION OF THE TWO-CENTRE TWO-ELECTRON INTEGRALS AND THE
C     CORE-ELECTRON ATTRACTION INTEGRALS IN LOCAL COORDINATES.
C     *
C     NOTATION. I=INPUT, O=OUTPUT.
C     NI        ATOMIC NUMBER OF FIRST  ATOM (I).
C     NJ        ATOMIC NUMBER OF SECOND ATOM (I).
C     R         INTERATOMIC DISTANCE, IN ATOMIC UNITS (I).
C     A(22)     LOCAL TWO-CENTRE TWO-ELECTRON  INTEGRALS, IN EV (O).
C     CORE()    LOCAL TWO-CENTRE CORE-ELECTRON INTEGRALS, IN EV (O).
C     *
C     POINT-CHARGE MULTIPOLES FROM TCA 1977 ARE EMPLOYED (SP BASIS).
C     BY DEFAULT, PENETRATION INTEGRALS ARE NEGLECTED.
C     HOWEVER, IF THE ADDITIVE TERM PO(9,N) FOR THE CORE DIFFERS
C     FROM THE ADDITIVE TERM PO(1,N) FOR SS, THE CORE-ELECTRON
C     ATTRACTION INTEGRALS ARE EVALUATED EXPLICITLY USING PO(9,N).
C     IN THIS CASE, PO(1,N) AND PO(7,N) ARE THE ADDITIVE TERMS
C     FOR THE MONOPOLES OF SS AND PP, RESPECTIVELY.
C     *
C     SPECIAL CONVENTION: NI=0 OR NJ=0 DENOTES AN EXTERNAL POINT
C     CHARGE WITHOUT BASIS ORBITALS. THE CHARGE IS 1 ATOMIC UNIT.
C     THE VALUES OF DD(I,0) AND PO(I,0) ARE DEFINED TO BE ZERO.
C     *
      USE LIMIT, ONLY: LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (SMALL=1.0D-06)
      COMMON
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./DPARM / LORBS(LMZ)
     ./INOPT2/ IN2(300)
     ./MULTIP/ DD(6,0:LMZ),PO(9,0:LMZ)
     ./PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
      DIMENSION A(22),CORE(10,2)
      DIMENSION X(69),PXX(33),PXY(69)
      DATA PXX/ 1.0D0   ,-0.5D0   ,-0.5D0   , 0.5D0   , 0.25D0  ,
     1          0.25D0  , 0.5D0   , 0.25D0  ,-0.25D0  , 0.25D0  ,
     2         -0.25D0  ,-0.125D0 ,-0.125D0 , 0.50D0  , 0.125D0 ,
     3          0.125D0 ,-0.5D0   ,-0.25D0  ,-0.25D0  ,-0.25D0  ,
     4          0.25D0  ,-0.125D0 , 0.125D0 ,-0.125D0 , 0.125D0 ,
     5          0.125D0 , 0.25D0  , 0.0625D0, 0.0625D0,-0.25D0  ,
     6          0.25D0  , 0.25D0  ,-0.25D0  /
      DATA PXY/ 1.0D0   ,-0.5D0   ,-0.5D0   , 0.5D0   , 0.25D0  ,
     1          0.25D0  , 0.5D0   , 0.25D0  ,-0.25D0  , 0.25D0  ,
     2         -0.25D0  ,-0.125D0 ,-0.125D0 ,-0.50D0  ,-0.50D0  ,
     3          0.5D0   , 0.25D0  , 0.25D0  , 0.5D0   , 0.25D0  ,
     4         -0.25D0  ,-0.25D0  ,-0.125D0 ,-0.125D0 , 0.5D0   ,
     5         -0.5D0   , 0.25D0  , 0.25D0  ,-0.25D0  ,-0.25D0  ,
     6         -0.25D0  , 0.25D0  ,-0.25D0  , 0.25D0  ,-0.125D0 ,
     7          0.125D0 ,-0.125D0 , 0.125D0 ,-0.125D0 , 0.125D0 ,
     8         -0.125D0 , 0.125D0 , 0.125D0 , 0.125D0 , 0.25D0  ,
     9          0.125D0 , 0.125D0 , 0.125D0 , 0.125D0 , 0.0625D0,
     A          0.0625D0, 0.0625D0, 0.0625D0,-0.25D0  , 0.25D0  ,
     B          0.25D0  ,-0.25D0  ,-0.25D0  , 0.25D0  , 0.25D0  ,
     C         -0.25D0  , 0.125D0 ,-0.125D0 ,-0.125D0 , 0.125D0 ,
     D         -0.125D0 , 0.125D0 , 0.125D0 ,-0.125D0 /
C *** INITIALIZATION.
      IOP = IN2(2)
      IF(NI.GT.0) THEN
         IORBS  = LORBS(NI)
         TORENI = TORE(NI)
      ELSE
         IORBS  = 1
         TORENI = ONE
      ENDIF
      IF(NJ.GT.0) THEN
         JORBS  = LORBS(NJ)
         TORENJ = TORE(NJ)
      ELSE
         JORBS  = 1
         TORENJ = ONE
      ENDIF
      R2        = R*R
      AEE       = (PO(1,NI)+PO(1,NJ))**2
C *** HYDROGEN - HYDROGEN.
      IF(IORBS.LE.1 .AND. JORBS.LE.1) THEN
         EE     = ONE/SQRT(R2+AEE)
         A(1)   = EE*EV
         CORE(1,1) = -TORENJ*A(1)
         CORE(1,2) = -TORENI*A(1)
C *** HEAVY ATOM - HYDROGEN
      ELSE IF(IORBS.GE.4 .AND. JORBS.LE.1) THEN
         DA     = DD(2,NI)
         QA     = DD(3,NI)
         TWOQA  = QA+QA
         ADE    = (PO(2,NI)+PO(1,NJ))**2
         AQE    = (PO(3,NI)+PO(1,NJ))**2
         X(1)   = R2+AEE
         X(2)   = R2+AQE
         X(3)   = (R+DA)**2+ADE
         X(4)   = (R-DA)**2+ADE
         X(5)   = (R-TWOQA)**2+AQE
         X(6)   = (R+TWOQA)**2+AQE
         X(7)   = R2+TWOQA*TWOQA+AQE
         DO 10 I=1,7
         X(I)   = PXY(I)/SQRT(X(I))
   10    CONTINUE
         A(1)   =  X(1)*EV
         A(2)   = (X(3)+X(4))*EV
         A(3)   = (X(1)+X(2)+X(5)+X(6))*EV
         A(4)   = (X(1)+X(2)+X(7))*EV
         CORE(1,1) = -TORENJ * A(1)
         CORE(2,1) = -TORENJ * A(2)
         CORE(3,1) = -TORENJ * A(3)
         CORE(4,1) = -TORENJ * A(4)
         CORE(1,2) = -TORENI * A(1)
C *** HYDROGEN - HEAVY ATOM
      ELSE IF(IORBS.LE.1 .AND. JORBS.GE.4) THEN
         DB     = DD(2,NJ)
         QB     = DD(3,NJ)
         TWOQB  = QB+QB
         AED    = (PO(1,NI)+PO(2,NJ))**2
         AEQ    = (PO(1,NI)+PO(3,NJ))**2
         X(1)   = R2+AEE
         X(2)   = R2+AEQ
         X(3)   = (R-DB)**2+AED
         X(4)   = (R+DB)**2+AED
         X(5)   = (R-TWOQB)**2+AEQ
         X(6)   = (R+TWOQB)**2+AEQ
         X(7)   = R2+TWOQB*TWOQB+AEQ
         DO 20 I=1,7
         X(I)   = PXY(I)/SQRT(X(I))
   20    CONTINUE
         A(1)   =  X(1)*EV
         A(5)   = (X(3)+X(4))*EV
         A(11)  = (X(1)+X(2)+X(5)+X(6))*EV
         A(12)  = (X(1)+X(2)+X(7))*EV
         CORE(1,1) = -TORENJ * A(1)
         CORE(1,2) = -TORENI * A(1)
         CORE(2,2) = -TORENI * A(5)
         CORE(3,2) = -TORENI * A(11)
         CORE(4,2) = -TORENI * A(12)
C *** HEAVY ATOM - HEAVY ATOM
      ELSE
         DA     = DD(2,NI)
         QA     = DD(3,NI)
         TWOQA  = QA+QA
         TWOQA2 = TWOQA*TWOQA
         RPDA2  = (R+DA)**2
         RMDA2  = (R-DA)**2
         RP2QA2 = (R+TWOQA)**2
         RM2QA2 = (R-TWOQA)**2
         ADE    = (PO(2,NI)+PO(1,NJ))**2
         AQE    = (PO(3,NI)+PO(1,NJ))**2
         ADD    = (PO(2,NI)+PO(2,NJ))**2
         ADQ    = (PO(2,NI)+PO(3,NJ))**2
         AQQ    = (PO(3,NI)+PO(3,NJ))**2
         TWOQAQ = TWOQA2+AQQ
         X(1)   = R2+AEE
         X(2)   = R2+AQE
         X(3)   = RPDA2+ADE
         X(4)   = RMDA2+ADE
         X(5)   = RM2QA2+AQE
         X(6)   = RP2QA2+AQE
         X(7)   = R2+TWOQA2+AQE
         X(8)   = RPDA2+ADQ
         X(9)   = RMDA2+ADQ
         X(10)  = R2+AQQ
         X(11)  = R2+TWOQAQ
         X(12)  = RP2QA2+AQQ
         X(13)  = RM2QA2+AQQ
         IF(NI.NE.NJ) THEN
            DB     = DD(2,NJ)
            QB     = DD(3,NJ)
            TWOQB  = QB+QB
            TWOQB2 = TWOQB*TWOQB
            TWOQBQ = TWOQB2+AQQ
            RPDB2  = (R+DB)**2
            RMDB2  = (R-DB)**2
            RP2QB2 = (R+TWOQB)**2
            RM2QB2 = (R-TWOQB)**2
            AED    = (PO(1,NI)+PO(2,NJ))**2
            AEQ    = (PO(1,NI)+PO(3,NJ))**2
            AQD    = (PO(3,NI)+PO(2,NJ))**2
            X(14)  = R2+AEQ
            X(15)  = RMDB2+AED
            X(16)  = RPDB2+AED
            X(17)  = RM2QB2+AEQ
            X(18)  = RP2QB2+AEQ
            X(19)  = R2+TWOQB2+AEQ
            X(20)  = RMDB2+AQD
            X(21)  = RPDB2+AQD
            X(22)  = R2+TWOQBQ
            X(23)  = RP2QB2+AQQ
            X(24)  = RM2QB2+AQQ
            X(25)  = R2+(DA-DB)**2+ADD
            X(26)  = R2+(DA+DB)**2+ADD
            X(27)  = (R+DA-DB)**2+ADD
            X(28)  = (R-DA+DB)**2+ADD
            X(29)  = (R-DA-DB)**2+ADD
            X(30)  = (R+DA+DB)**2+ADD
            X(31)  = RPDA2+TWOQB2+ADQ
            X(32)  = RMDA2+TWOQB2+ADQ
            X(33)  = RMDB2+TWOQA2+AQD
            X(34)  = RPDB2+TWOQA2+AQD
            X(35)  = (R+DA-TWOQB)**2+ADQ
            X(36)  = (R-DA-TWOQB)**2+ADQ
            X(37)  = (R+DA+TWOQB)**2+ADQ
            X(38)  = (R-DA+TWOQB)**2+ADQ
            X(39)  = (R+TWOQA-DB)**2+AQD
            X(40)  = (R+TWOQA+DB)**2+AQD
            X(41)  = (R-TWOQA-DB)**2+AQD
            X(42)  = (R-TWOQA+DB)**2+AQD
            X(43)  = R2+FOUR*(QA-QB)**2+AQQ
            X(44)  = R2+FOUR*(QA+QB)**2+AQQ
            X(45)  = R2+TWOQA2+TWOQBQ
            X(46)  = RM2QB2+TWOQAQ
            X(47)  = RP2QB2+TWOQAQ
            X(48)  = RP2QA2+TWOQBQ
            X(49)  = RM2QA2+TWOQBQ
            X(50)  = (R+TWOQA-TWOQB)**2+AQQ
            X(51)  = (R+TWOQA+TWOQB)**2+AQQ
            X(52)  = (R-TWOQA-TWOQB)**2+AQQ
            X(53)  = (R-TWOQA+TWOQB)**2+AQQ
            X(54)  = (R-QB)**2+(DA-QB)**2+ADQ
            X(55)  = (R+QB)**2+(DA-QB)**2+ADQ
            X(56)  = (R-QB)**2+(DA+QB)**2+ADQ
            X(57)  = (R+QB)**2+(DA+QB)**2+ADQ
            X(58)  = (R+QA)**2+(QA-DB)**2+AQD
            X(59)  = (R-QA)**2+(QA-DB)**2+AQD
            X(60)  = (R+QA)**2+(QA+DB)**2+AQD
            X(61)  = (R-QA)**2+(QA+DB)**2+AQD
            QMADD  = (QA-QB)**2+AQQ
            QPADD  = (QA+QB)**2+AQQ
            X(62)  = (R+QA-QB)**2+QMADD
            X(63)  = (R+QA+QB)**2+QMADD
            X(64)  = (R-QA-QB)**2+QMADD
            X(65)  = (R-QA+QB)**2+QMADD
            X(66)  = (R+QA-QB)**2+QPADD
            X(67)  = (R+QA+QB)**2+QPADD
            X(68)  = (R-QA-QB)**2+QPADD
            X(69)  = (R-QA+QB)**2+QPADD
            DO 30 I=1,69
            X(I)   = PXY(I)/SQRT(X(I))
   30       CONTINUE
            EE     = X(1)
            DZE    = X(3) +X(4)
            QZZE   = X(2) +X(5) +X(6)
            QXXE   = X(2) +X(7)
            EDZ    = X(15)+X(16)
            EQZZ   = X(14)+X(17)+X(18)
            EQXX   = X(14)+X(19)
            DXDX   = X(25)+X(26)
            DZDZ   = X(27)+X(28)+X(29)+X(30)
            X89    = X(8) +X(9)
            X2021  = X(20)+X(21)
            DZQXX  = X89  +X(31)+X(32)
            QXXDZ  = X2021+X(33)+X(34)
            DZQZZ  = X89  +X(35)+X(36)+X(37)+X(38)
            QZZDZ  = X2021+X(39)+X(40)+X(41)+X(42)
            X1011  = X(10)+X(11)
            X1213  = X(12)+X(13)
            X2324  = X(23)+X(24)
            QXXQXX = X1011+X(22)+X(43)+X(44)
            QXXQYY = X1011+X(22)+X(45)
            QXXQZZ = X1011+X2324+X(46)+X(47)
            QZZQXX = X(10)+X1213+X(22)+X(48)+X(49)
            QZZQZZ = X(10)+X1213+X2324+X(50)+X(51)+X(52)+X(53)
            DXQXZ  = X(54)+X(55)+X(56)+X(57)
            QXZDX  = X(58)+X(59)+X(60)+X(61)
            QXZQXZ = ZERO
            DO 40 I=62,69
            QXZQXZ = QXZQXZ+X(I)
   40       CONTINUE
         ELSE
            TWODA  = DA+DA
            X(14)  = R2+ADD
            X(15)  = RP2QA2+TWOQAQ
            X(16)  = RM2QA2+TWOQAQ
            X(17)  = R2+TWODA**2+ADD
            X(18)  = (R-TWODA)**2+ADD
            X(19)  = (R+TWODA)**2+ADD
            X(20)  = RPDA2+TWOQA2+ADQ
            X(21)  = RMDA2+TWOQA2+ADQ
            X(22)  = (R+DA-TWOQA)**2+ADQ
            X(23)  = (R-DA-TWOQA)**2+ADQ
            X(24)  = (R+DA+TWOQA)**2+ADQ
            X(25)  = (R-DA+TWOQA)**2+ADQ
            X(26)  = R2+FOUR*TWOQA2+AQQ
            X(27)  = R2+TWOQA2+TWOQAQ
            X(28)  = (R+TWOQA+TWOQA)**2+AQQ
            X(29)  = (R-TWOQA-TWOQA)**2+AQQ
            RMQA2  = (R-QA)**2
            RPQA2  = (R+QA)**2
            DMADD  = (DA-QA)**2+ADQ
            DPADD  = (DA+QA)**2+ADQ
            X(30)  = RMQA2+DMADD
            X(31)  = RPQA2+DMADD
            X(32)  = RMQA2+DPADD
            X(33)  = RPQA2+DPADD
            DO 50 I=1,33
            X(I)   = PXX(I)/SQRT(X(I))
   50       CONTINUE
            EE     = X(1)
            DZE    = X(3) +X(4)
            QZZE   = X(2) +X(5) +X(6)
            QXXE   = X(2) +X(7)
            EDZ    =-DZE
            EQZZ   = QZZE
            EQXX   = QXXE
            DXDX   = X(14)+X(17)
            DZDZ   = X(14)+X(18)+X(19)
            X89    = X(8) +X(9)
            DZQXX  = X89  +X(20)+X(21)
            QXXDZ  =-DZQXX
            DZQZZ  = X89  +X(22)+X(23)+X(24)+X(25)
            QZZDZ  =-DZQZZ
            X1010  = X(10)+X(10)*PT5
            X1111  = X(11)+X(11)
            X1213  = X(12)+X(13)
            QXXQXX = X1010+X1111+X(26)
            QXXQYY = X1111+X(10)+X(27)
            QXXQZZ = X(10)+X(11)+X1213+X(15)+X(16)
            QZZQXX = QXXQZZ
            QZZQZZ = X1010+X1213+X1213+X(28)+X(29)
            DXQXZ  = X(30)+X(31)+X(32)+X(33)
            QXZDX  =-DXQXZ
            QXZQXZ = QXXQZZ
         ENDIF
         A(1)  = EE
         A(2)  = DZE
         A(3)  = EE + QZZE
         A(4)  = EE + QXXE
         A(5)  = EDZ
         A(6)  = DZDZ
         A(7)  = DXDX
         A(8)  = EDZ + QZZDZ
         A(9)  = EDZ + QXXDZ
         A(10) = QXZDX
         A(11) = EE  + EQZZ
         A(12) = EE  + EQXX
         A(13) = DZE + DZQZZ
         A(14) = DZE + DZQXX
         A(15) = DXQXZ
         A(16) = EE + EQZZ + QZZE + QZZQZZ
         A(17) = EE + EQZZ + QXXE + QXXQZZ
         A(18) = EE + EQXX + QZZE + QZZQXX
         A(19) = EE + EQXX + QXXE + QXXQXX
         A(20) = QXZQXZ
         A(21) = EE + EQXX + QXXE + QXXQYY
         A(22) = PT5*(A(19)-A(21))
C        ORIGINAL (PP,PP) EXCHANGE INTEGRALS FROM FIRST MNDO PAPER.
C        THE FORMULA IS NOT ROTATIONALLY INVARIANT AND THUS OBSOLETE.
C        THE OLD CODE IS DOCUMENTED IN THE FORM OF COMMENT LINES.
C        IF(IOP.EQ.-4) THEN
C           IF(NI.NE.NJ) THEN
C              A(22) = PT25/SQRT(R2+TWO*(QA-QB)**2+AQQ)
C    1                +PT25/SQRT(R2+TWO*(QA+QB)**2+AQQ)
C    2                -PT5 /SQRT(R2+TWO*(QA**2+QB**2)+AQQ)
C           ELSE
C              A(22) = X(10)+X1111
C    1                +PT25/SQRT(R2+TWO*TWOQA2+AQQ)
C           ENDIF
C        ENDIF
         DO 60 I=1,22
         A(I)   = A(I)*EV
   60    CONTINUE
         CORE(1,1) = -TORENJ * A(1)
         CORE(2,1) = -TORENJ * A(2)
         CORE(3,1) = -TORENJ * A(3)
         CORE(4,1) = -TORENJ * A(4)
         CORE(1,2) = -TORENI * A(1)
         CORE(2,2) = -TORENI * A(5)
         CORE(3,2) = -TORENI * A(11)
         CORE(4,2) = -TORENI * A(12)
      ENDIF
C *** CALCULATE THE NUCLEAR ATTRACTION INTEGRALS IN LOCAL COORDINATES
C     WITH A SEPARATE ADDITIVE TERM FOR THE CORE (SP BASIS).
C     OMIT THE CALCULATION FOR IDENTICAL ADDITIVE TERMS (SS=CORE).
C     THIS OPTION IS ONLY VALID FOR SP-TYPE INTEGRALS IN MNDO/d.
      IF(IOP.GT.-10) RETURN
      ACI    = PO(9,NI)
      ACJ    = PO(9,NJ)
C *** ELECTRONS AT ATOM A (NI) AND CORE OF ATOM B (NJ).
      IF(ABS(ACJ-PO(1,NJ)).GT.SMALL) THEN
         CORE(1,1) = -TORENJ*EV/SQRT(R2+(PO(1,NI)+ACJ)**2)
         IF(IORBS.GE.4) THEN
            DA     = DD(2,NI)
            QA     = DD(3,NI)
            AEJ    = (PO(7,NI)+ACJ)**2
            ADJ    = (PO(2,NI)+ACJ)**2
            AQJ    = (PO(3,NI)+ACJ)**2
            TWOQA  = QA+QA
            X(1)   = R2+AEJ
            X(2)   = R2+AQJ
            X(3)   = (R+DA)**2+ADJ
            X(4)   = (R-DA)**2+ADJ
            X(5)   = (R-TWOQA)**2+AQJ
            X(6)   = (R+TWOQA)**2+AQJ
            X(7)   = R2+TWOQA*TWOQA+AQJ
            DO 70 I=1,7
            X(I)   = PXY(I)/SQRT(X(I))
   70       CONTINUE
            TOREV  = TORENJ*EV
            CORE(2,1) = -TOREV * (X(3)+X(4))
            CORE(3,1) = -TOREV * (X(1)+X(2)+X(5)+X(6))
            CORE(4,1) = -TOREV * (X(1)+X(2)+X(7))
         ENDIF
      ENDIF
C *** ELECTRONS AT ATOM B (NJ) AND CORE OF ATOM A (NI).
      IF(ABS(ACI-PO(1,NI)).GT.SMALL) THEN
         CORE(1,2) = -TORENI*EV/SQRT(R2+(PO(1,NJ)+ACI)**2)
         IF(JORBS.GE.4) THEN
            DB     = DD(2,NJ)
            QB     = DD(3,NJ)
            AEI    = (PO(7,NJ)+ACI)**2
            ADI    = (PO(2,NJ)+ACI)**2
            AQI    = (PO(3,NJ)+ACI)**2
            TWOQB  = QB+QB
            X(1)   = R2+AEI
            X(2)   = R2+AQI
            X(3)   = (R+DB)**2+ADI
            X(4)   = (R-DB)**2+ADI
            X(5)   = (R-TWOQB)**2+AQI
            X(6)   = (R+TWOQB)**2+AQI
            X(7)   = R2+TWOQB*TWOQB+AQI
            DO 80 I=1,7
            X(I)   = PXY(I)/SQRT(X(I))
   80       CONTINUE
            TOREV  = TORENI*EV
            CORE(2,2) =  TOREV * (X(3)+X(4))
            CORE(3,2) = -TOREV * (X(1)+X(2)+X(5)+X(6))
            CORE(4,2) = -TOREV * (X(1)+X(2)+X(7))
         ENDIF
      ENDIF
      RETURN
      END
