      SUBROUTINE CPEFRM(Q,ETAI,NU,RES,N,NCPE,ENTRPY,DBGPRT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./CPEOPT/ NFLCPE,NPTCPE,NCPEZ
     ./CPEPAR/ CPEZ(LMZ),CPEED0(LMZ),CPEQ0(LMZ),ATPS(LMZ),
     .         ANPS(LMZ),ATPT(LMZ), ANPT(LMZ),CPESL(LMZ),CPEEM0(LMZ),
     .         CPEM0(LMZ),CPERS(LMZ),CPERR(LMZ),CMXCAT(LMZ),CMXAN(LMZ),
     .         CPEFT(LMZ),CPERC(LMZ),CPEQM(LMZ),CPEPS,CPEZA(LMZ),
     .         CPEZB(LMZ)
C     ------------------------------------------------
C     INPUT/OUTPUT VARIABLES
      INTEGER N,NCPE
      DOUBLE PRECISION ETA,NU,RES,Q
      DIMENSION ETAI(NCPE,NCPE),NU(NCPE),RES(NCPE),Q(N)
      LOGICAL DBGPRT
C     ------------------------------------------------
      DOUBLE PRECISION LUMON,HOMON
      DIMENSION LUMON(N),HOMON(N)
      DOUBLE PRECISION LUMODN,HOMODN
      DIMENSION LUMODN(N),HOMODN(N)
      DOUBLE PRECISION HOMO,LUMO
      DIMENSION HOMO(N),LUMO(N)
      DOUBLE PRECISION DVEC,DMAT
      DIMENSION DVEC(NCPE),DMAT(NCPE,N)
      DOUBLE PRECISION GAMMA,GAMMAI,OMEGA
      DIMENSION GAMMA(N,N),GAMMAI(N,N),OMEGA(N)
      DOUBLE PRECISION MONRES,ENTRPY
      DIMENSION MONRES(N)
      DOUBLE PRECISION BETA
      DIMENSION BETA(N)
C     TEMPORARY ARRAYS
      DOUBLE PRECISION TMPN,TMPN2,TMPNC,TMPNC2,TMPNC3
      DIMENSION TMPN(N),TMPN2(N),TMPNC(NCPE),TMPNC2(NCPE),TMPNC3(NCPE)
      DOUBLE PRECISION A,B
      DOUBLE PRECISION N0,N0V
      DIMENSION N0V(N)
      DOUBLE PRECISION NVAL,NFELEC
C     -------------------------------------------------
      INTEGER I,J,ICPE,NATI
C      LOGICAL DBGPRT
C      CALL PTM(NCPE,ETAI,NCPE,4,4,"ETAI",6,.FALSE.)
C      CALL PTV(NU,NCPE,4,"NU",6,.FALSE.)
C      CALL PTV(Q,N,N,"Q",6,.FALSE.)
C      IF ( NPTCPE.GT.0 ) THEN
C         DBGPRT = .TRUE.
C      ELSE
C         DBGPRT = .FALSE.
C      END IF
C      IF(DBGPRT)WRITE(6,'(A)')"ENTERED CPEFRM - CPE FERMI POPULATION"
      CALL CLRV(DVEC,NCPE)
      CALL CLRM(NCPE,DMAT,N)
      CALL CLRV(N0V,N)
      CALL CLRV(TMPN,N)
      CALL CLRV(TMPN2,N)
      CALL CLRV(BETA,N)
      J = 0
      DO 10 I=1,NCPE,4
         DVEC(I)   = ONE
         J=J+1
         DMAT(I,J) = ONE
 10   CONTINUE
      IF(DBGPRT) THEN
         WRITE(6,'(A3,A3,1X,6A7)')"I","Z",
     1        "N0","Nhomo","Nlumo",
     1        "Q","dNhomo","dNlumo"
      END IF
      NFELEC = ZERO
      DO 100 I=1,N
         NATI   =  NAT(I)
         NVAL   =  TORE(NATI) - Q(I)
         N0V(I) =  TORE(NATI) - CMXCAT(NATI)
C     NFELEC = THE # OF ELECTRONS TO BE DETERMINED
C     FROM THE FERMI POPULATION
         NFELEC  =  NFELEC + NVAL - N0V(I)
C ---- THE CHANGE IN THE NUMBER OF ELECTRONS ----
C        THAT REPRESENT THE HOMO AND LUMO
C     THE Q(I) GETS US TO THE NEUTRAL ATOM
C     AND THE CMXAN(NATI) GETS US TO THE ANION
         LUMODN(I)  = Q(I) - CMXAN(NATI)
C     THE Q(I) GETS US TO THE NEUTRAL ATOM
         HOMODN(I) = Q(I)
C -----------------------------------------------
C     HERE ARE THE OCCUPATIONS THAT THE RESULTING HOMOS
C     AND LUMOS PROVIDE.  IF THE POPULATION (0.0->1.0) OF
C     HOMO(I) AND LUMO(I) ARE PHOMO(I) and PLUMO(I), THEN
C     THE NUMBER OF ELECTRONS ON ATOM I IS
C     N(I) = N0V(I) + PHOMO(I)*CATN(I) + PLUMO(I)*ANN(I)
C     AND THE PARTIAL CHARGE IS
C     QPART(I) = TORE(NAT(I)) - N(I)
C     (NOTE THAT TORE(NAT(I)) IS THE NUMBER OF
C     (VALENCE) ELECTRONS OF THE NEUTRAL ATOM)
         HOMON(I) =  CMXCAT(NATI)
         LUMON(I) = -CMXAN(NATI)
         IF(DBGPRT) THEN
            WRITE(6,'(2I3,1X,6F7.2)')I,NATI,
     1           N0V(I),HOMON(I),LUMON(I),Q(I),HOMODN(I),LUMODN(I)
         END IF
 100  CONTINUE
C     GAMMA = TRANS(D).ETAI.D
      CALL YCSA(NCPE,ETAI,DMAT,N,GAMMA)
C      CALL PTM(NCPE,DMAT,N,NCPE,N,"DMAT",6,.FALSE.)
C      CALL PTM(N,GAMMA,N,N,N,"GAMMA",6,.FALSE.)
C     GAMMAI = INVERSE(GAMMA)
      CALL SVDINV(N,GAMMA,N,GAMMAI)
      CALL PTM(N,GAMMAI,N,N,N,"GAMMA INVERSE",6,.FALSE.)
C     OMEGA = OMEGAI.TRANSP(D).ETAI.NU
      CALL YSV(NCPE,ETAI,NU,TMPNC)
      CALL YCV(N,DMAT,NCPE,TMPNC,TMPN)
      CALL YSV(N,GAMMAI,TMPN,OMEGA)
      CALL PTV(OMEGA,N,N,"Gamma.D.Ei.m",6,.FALSE.)
      IF ( DBGPRT ) THEN
         WRITE(6,'(A3,A3,1X,2A15)')"I","Z","HOMO","LUMO"
      END IF
      DO 200 I=1,N
C         LUMO(I) = LUMODN(I) *GAMMAI(I,I)+OMEGA(I)
C         HOMO(I) = HOMODN(I) *GAMMAI(I,I)+OMEGA(I)
         LUMO(I) = OMEGA(I)
         HOMO(I) = OMEGA(I)
         DO 210 J=1,N
            LUMO(I) = LUMO(I)-LUMODN(I)*GAMMAI(I,J)
            HOMO(I) = HOMO(I)-HOMODN(I)*GAMMAI(I,J)
 210     CONTINUE
         IF ( DBGPRT ) THEN
            WRITE(6,'(2I3,1X,2F15.5)')I,NAT(I),HOMO(I),LUMO(I)
         END IF
 200  CONTINUE
      CALL FERMIL(NFELEC,N0V,Q,N,HOMO,LUMO,HOMON,LUMON,
     a     CPEFT,MONRES,ENTRPY,DBGPRT)
C      CALL PTV(MONRES,N,N,"MONRES",6,.FALSE.)
      CALL ADDVV(ONE,MONRES,N,ONE,TMPN,TMPN2)
      CALL YSV(N,GAMMAI,TMPN2,BETA)
      CALL YAV(NCPE,DMAT,N,BETA,TMPNC2)
      CALL YSV(NCPE,ETAI,TMPNC2,TMPNC3)
      CALL ADDVV(ONE,TMPNC3,NCPE,-ONE,TMPNC,RES)
C      CALL PTV(RES,NCPE,4,"RES",6,.FALSE.)
C      IF(DBGPRT)WRITE(6,'(A)')"EXITING CPEFRM"
      END SUBROUTINE
      SUBROUTINE FERMIL(
     a     NFELEC,N0V,Q,
     a     N,
     a     HOMO,LUMO,NHOMO,NLUMO,
     a     FT,MONRES,ENTRPY,DBGPRT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (LM1=600)
      PARAMETER (LMX=9*LM1)
      PARAMETER (LMZ=86)
      COMMON
     .     /ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     .     /PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)
      INTEGER N
      LOGICAL DBGPRT
      DOUBLE PRECISION NFELEC,HOMO,LUMO,NHOMO,NLUMO,MONRES,ENTRPY
      DIMENSION HOMO(N),LUMO(N),MONRES(N),NHOMO(N),NLUMO(N)
      DOUBLE PRECISION N0V,Q,FT
      DIMENSION N0V(N),Q(N),FT(LMZ)
      DOUBLE PRECISION LEV
      DIMENSION LEV(2*N)
      DOUBLE PRECISION NMAP
      DIMENSION NMAP(2*N)
      DOUBLE PRECISION NLO,NHI,NMID,LLO,LHI,LMID
      DOUBLE PRECISION NTOL,KT,PLEV,BK
      INTEGER I,J,IMAP,AMAP,M,ITER
      DIMENSION IMAP(2*N),AMAP(2*N)
      LOGICAL LISODD
      BK = 0.3166815213D-05
      TOL = 1.0D-8
      M = 2*N
      CALL CLRV(MONRES,N)
      ENTRPY = 0.0D0
C ----------------------------------------------------
C  CREATE SCRATCH ARRAYS FOR FAST LOOKUP OF DATA
C     COPY THE HOMO'S AND LUMO'S INTO A COMBINED ARRAY
      DO 10 I=1,N
         J = (I-1)*2+1
C     ODD INDICIES ARE HOMOS
         LEV(J)  = HOMO(I)
C     EVEN INDICIES ARE LUMOS
         LEV(J+1)= LUMO(I)
 10   CONTINUE
C     SORT THE ARRAY IN ASCENDING ORDER
      CALL SORTV(LEV,M,IMAP,"A")
C     IMAP IS NOW THE MAPPING BACK TO THE UNSORTED ARRAY
C     NOW CREATE AMAP TO BE THE MAPPING TO THE ATOM
C     i.e., IT IS THE MAPPING FROM THE I'th ENERGY LEVEL
C     TO THE AMAP(I)'th ATOM
C     NMAP IS LIKE IMAP, BUT IT CONTAINS THE FLOATS
C     CORRESPONDING TO THE NUMBER OF ELECTRONS THAT
C     THE I'th ENERGY LEVEL CONTRIBUTES AT FULL OCCUPATION
      DO 20 I=1,M
         IF ( LISODD(IMAP(I)) ) THEN
C     ODD INDICIES ARE HOMOS
            AMAP(I)=(IMAP(I)-1)/2+1
            NMAP(I)=NHOMO(AMAP(I))
         ELSE
C     EVEN INDICIES ARE LUMOS
            AMAP(I)=IMAP(I)/2
            NMAP(I)=NLUMO(AMAP(I))
         END IF
 20   CONTINUE
C ---------------------------------------------------
C     MAKE AN INITIAL GUESS AT THE FERMI LEVEL
C     TO START THE RECURSIVE BISECTION SEARCH
      IF ( DBGPRT ) THEN
         WRITE(6,'(A,F10.4)')"# OF ELECTRONS TO BE FERMI POPULATED=",
     1        NFELEC
      END IF
      IF ( NFELEC .LT. TOL ) RETURN
      NHI  = 0.0D0
      NMID = 0.0D0
      NLO  = 0.0D0
      DO 100 I=1,M
         NHI = NHI + NMAP(I)
         IF ( NHI .LT. NFELEC ) THEN
            NMID = NHI
            LMID = LEV(IMAP(I))
         END IF
 100  CONTINUE
      LLO = LEV(1)-0.1D0
      LHI = LEV(M)+0.1D0
C --------------------------------------------------
C     RECURSIVE BISECTION SEARCH
      ITER = 0
      DO 200
         ITER = ITER + 1
C         WRITE(6,'(A,6E15.5)')"NLO,NMID,NHI,LLO,LMID,LHI",
C     a        NLO,NMID,NHI,LLO,LMID,LHI
C     TAKE A BISECTION STEP
C         WRITE(6,'(A,2E15.5)')"VAL/TOL",ABS(NFELEC - NMID), TOL
C         IF ( ITER .LT. 3 ) THEN
C            LMID = WTAVG(NLO,LHI,NHI,LLO)
C         ELSE
         LMID = (LHI+LLO)/2.0D0
C         END IF
C ---------------------------------------------------
C     CALCULATE THE NUMBER OF ELECTRONS AND ENTROPY
C     AT THIS FERMI LEVEL
         NMID = 0.0D0
         ENTRPY = 0.0D0
         DO 250 I=1,M
C     OCCUPATION FRACTION (PROBABILITY)
            KT = BK*FT(NAT(AMAP(I)))
C            WRITE(6,*)KT,FT(NAT(AMAP(I)))
            PLEV = NMAP(I) / ( 1.0D0 + EXP((LEV(I)-LMID)/KT) )
C            PLEV = 1.D0 / ( 1.0D0 + EXP((LEV(I)-LMID)/KT) )
C            NMID = NMID + NMAP(I)*PLEV
            NMID = NMID + PLEV
C            WRITE(6,'(A,E15.5)')"   NMID +=",NMID
            IF ( PLEV .GT. 1.D-20 ) THEN
               ENTRPY = ENTRPY - KT*PLEV*LOG(PLEV)
            END IF
 250     CONTINUE
C --------------------------------------------------
         IF ( ITER .GT. 1 ) THEN
            WRITE(6,'(2(A,E12.4),A,2E12.4)')"TEST dN=",
     1           ABS(NFELEC - NMID),
     2           " TOL=",TOL,
     3           " WHERE N0,NT=",NFELEC,NMID
         END IF
         IF ( ABS(NFELEC - NMID) .LT. TOL ) THEN
            IF ( ITER .GT. 1 ) EXIT
         ELSE IF ( NFELEC .GT. NMID ) THEN
C     WE NEED TO RAISE THE FERMI LEVEL
C            WRITE(6,*)"RAISE FERMI"
            NLO = NMID
            LLO = LMID
         ELSE IF ( NFELEC .LT. NMID ) THEN
C     WE NEED TO LOWER THE FERMI LEVEL
            NHI = NMID
            LHI = LMID
         ELSE
            WRITE(6,'(A)')"RECURSIVE BISECTION ERROR IN CPEFERMI.f"
         END IF
 200  CONTINUE
C ------- WE HAVE CONVERGED ------------------------
C -----------------------------------------
C     USE THIS FERMI LEVEL TO CALCULATE THE
C     NUMBER OF ELECTRONS ON EACH ATOM
      IF ( DBGPRT ) THEN
         WRITE(6,'(A,I5,A)')"CONVERGED IN ",ITER," ITERATIONS"
         WRITE(6,'(2A6,3A12)')"Level","Atom","ProbOcc",
     1        "N(1.0)","N(ProbOcc)"
      END IF
      DO 300 I=1,M
         KT = BK*FT(NAT(AMAP(I)))
         PLEV = 1.0D0 / ( 1.0D0 + EXP((LEV(I)-LMID)/KT) )
C         PLEV = NMAP(I) / ( 1.0D0 + EXP((LEV(I)-LMID)/KT) )
         IF ( DBGPRT ) THEN
            WRITE(6,'(2I6,3F12.6)')I,AMAP(I),PLEV,NMAP(I),
     1           PLEV*NMAP(I)
C     a           MONRES(AMAP(I))+PLEV
         END IF
         MONRES(AMAP(I)) = MONRES(AMAP(I))+PLEV*NMAP(I)
 300  CONTINUE
      IF ( DBGPRT ) THEN
         WRITE(6,'(A)')"DESIRED ATOM PARTIAL CHARGES"
         DO 310 I=1,N
            WRITE(6,'(2I6,F12.6)')I,NAT(I),
     1           TORE(NAT(I))-MONRES(I)-N0V(I)
 310     CONTINUE
      END IF
      IF ( DBGPRT ) THEN
         WRITE(6,'(A)')"PARTIAL CHARGES NEEDED TO REACH DESIRED VALUES"
      END IF
      DO 350 I=1,N
c                    |--- NEEDED CHANGE IN CHARGE STATE ---|
c                    |------ DESIRED CHARGE STATE ----|
c                    |- # OF ELECTRONS -|
         MONRES(I) =  TORE(NAT(I))-MONRES(I)-N0V(I)-Q(I)
         IF ( DBGPRT ) THEN
            WRITE(6,'(2I6,F12.6)')I,NAT(I),MONRES(I)
         END IF
 350  CONTINUE
      END SUBROUTINE
c$$$      SUBROUTINE POPLEV(N,NLEV,ELEV,PLEV,NELECI,ENTRPY,KT)
c$$$C*********************************************************
c$$$C     GIVEN A SET OF ENERGY LEVELS AND NUMBER OF ELECTRONS,
c$$$C     POPULATE THE ENERGY LEVELS USING THE FERMI LEVEL THAT
c$$$C     REPRODUCES THE NUMBER OF ELECTRONS.
c$$$C*********************************************************
c$$$      INTEGER N,NLEV
c$$$      DOUBLE PRECISION ELEV,PLEV,ENTRPY,KT,NELECI
c$$$      DIMENSION ELEV(NLEV),PLEV(NLEV),NELECI(N)
c$$$C     WORK ARRAYS
c$$$      INTEGER I,J,K,IA,NITER
c$$$      INTEGER IELEV
c$$$      DIMENSION IELEV(NLEV)
c$$$      DOUBLE PRECISION MU,SLOPE,INTER
c$$$      DOUBLE PRECISION NELEC,NTEST,NHI,NLO,MUHI,MULO,DN
c$$$
c$$$C     THE NUMBER OF ELECTRONS IS THE NUMBER OF ATOMS
c$$$      NELEC = 1.000000000000000D0 * N
c$$$C     SORT THE ARRAY OF ENERGY LEVELS FROM MOST NEGATIVE
c$$$C     TO LEAST NEGATIVE AND STORE THE SORT MAPPING INTO IELEV
c$$$      CALL SORTV(ELEV,NLEV,IELEV)
c$$$
c$$$C     OUT INITIAL GUESS AT THE FERMI LEVEL IS THE FERMI
c$$$C     LEVEL AT ZERO TEMPERATURE (THE N'TH ENERGY LEVEL)
c$$$C     DETERMINE THE FERMI LEVEL AT ZERO TEMPERATURE
c$$$      MU = ELEV(N)
c$$$
c$$$C     USE THIS AS A STARTING POINT FOR A RECURSIVE BISECTION
c$$$C     THAT FINDS THE FERMI LEVEL NEEDED (AT KT) THAT
c$$$C     REPRODUCES THE NUMBER OF ELECTRONS IN OUR SYSTEM
c$$$
c$$$C     MU = ELEV(1) CORRESPONDS TO NELEC=0.5
c$$$C     MU = 0       CORRESPONDS TO NELEC=N+N*N
c$$$      NHI  = NELEC+NELEC**2
c$$$      NLO  = 0.5D0
c$$$      MUHI = 0.0D0
c$$$      MULO = ELEV(1)
c$$$
c$$$      NITER = 0
c$$$      DO 1000
c$$$C     NTEST IS OUR TRIAL NUMBER OF ELECTRONS
c$$$         NTEST = 0.0D0
c$$$
c$$$C     -----------------------------------------------
c$$$C     CALCULATE THE FERMI DISTRIBUTION OF EACH LEVEL
c$$$C     -----------------------------------------------
c$$$         DO 80 I=1,NLEV
c$$$            PLEV(I) = 1.0D0 / ( 1.0D0 + EXP((ELEV(I)-MU)/KT) )
c$$$ 80      CONTINUE
c$$$
c$$$C     ----------------------------------------------
c$$$C     POPULATE LEVELS AND COUNT ELECTRONS
c$$$C     ----------------------------------------------
c$$$C     START WITH EACH ATOM HAVING NO ELECTRONS
c$$$         CALL CLRV(NELECI,N)
c$$$         DO 100 I=1,NLEV
c$$$C     WHAT ATOM DOES THIS LEVEL CORRESPOND TO?
c$$$            CALL LEV2AT(I,NLEV,N,IELEV,IA)
c$$$C     ADD TO THE TOTAL NUMBER OF ELECTRONS
c$$$            NTEST = NTEST + PLEV(I)
c$$$C     CALCULATE THE NUMBER OF ELECTRONS ON ATOM I
c$$$            NELECI(IA) = NELECI(IA) + PLEV(I)
c$$$ 100     CONTINUE
c$$$
c$$$C     WE CAN NOW CALCULATE THE ENTROPY
c$$$C     THIS IS ONLY USEFUL IF WE ARE CONVERGED
c$$$         ENTRPY = 0.0D0
c$$$         DO 85 I=1,NLEV
c$$$            IF ( PLEV(I) .LT. 1.0D-20 ) CYCLE
c$$$            ENTRPY = ENTRPY - PLEV(I)*LOG(PLEV(I))
c$$$ 85      CONTINUE
c$$$
c$$$C     -----------------------------------------------
c$$$C     CHECK FOR CONVERGENCE
c$$$C     -----------------------------------------------
c$$$C     DID WE GET THE RIGHT NUMBER OF ELECTRONS?
c$$$         DN = NTEST-NELEC
c$$$C         WRITE(6,'(A,4F25.5)')"MU,NTEST,DN",MU,NTEST,NELEC,DN
c$$$         IF ( ABS( DN ) .LT. 1.0D-10 ) THEN
c$$$            EXIT
c$$$         ELSE IF ( DN .GT. 0.0D0 ) THEN
c$$$C     WE NEED TO LOWER MU BECAUSE WE HAVE TOO MANY ELECTRONS
c$$$            MUHI = MU
c$$$            NHI  = NTEST
c$$$            MU = (MULO+MU)/2.0D0
c$$$         ELSE
c$$$C     WE NEED TO RAISE MU BECAUSE WE DON'T HAVE ENOUGH ELECTRONS
c$$$            MULO = MU
c$$$            NLO  = NTEST
c$$$            MU = (MUHI+MU)/2.0D0
c$$$         END IF
c$$$ 1000 CONTINUE
c$$$      END SUBROUTINE
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$      SUBROUTINE LEV2AT(ILEV,NLEV,N,ATMAP,IA)
c$$$C*******************************************************
c$$$C     GIVES THE ATOM INDEX IA FROM THE LEVEL INDEX ILEV
c$$$C*******************************************************
c$$$      INTEGER ILEV,N,ATMAP,IA
c$$$      DIMENSION ATMAP(NLEV)
c$$$      IA = ATMAP(ILEV)
c$$$      IF ( IA .GT. N ) THEN
c$$$C     THIS IS A LUMO AND THE ATOM IS THE ROW OF THE LUMO MATRIX
c$$$         IF ( MOD((IA-N),N) .GT. 0 ) THEN
c$$$            IA = (IA-N)/N + 1
c$$$         ELSE
c$$$            IA = (IA-N)/N
c$$$         END IF
c$$$      END IF
c$$$      END SUBROUTINE
c$$$
c$$$
c$$$
c$$$
c$$$      SUBROUTINE SORTV(V,N,IMAP)
c$$$C*******************************************************
c$$$C     SORTS VECTOR V(N) IN ASCENDING ORDER
c$$$C     AND KEEPS THE INDEX MAPPING OF THE SORT IN IMAP(N)
c$$$C*******************************************************
c$$$      INTEGER N,IMAP
c$$$      DOUBLE PRECISION V
c$$$      DIMENSION IMAP(N),V(N)
c$$$      INTEGER I,IMIN
c$$$      DOUBLE PRECISION VT,VMIN
c$$$      DIMENSION VT(N)
c$$$      DO 10 I=1,N
c$$$         CALL MINLV(V,N,IMIN,VMIN)
c$$$         IMAP(I) = IMIN
c$$$         VT(I)   = VMIN
c$$$         V(IMIN)    = 1.0D+10
c$$$ 10   CONTINUE
c$$$      CALL CPYV(VT,N,V)
c$$$      END SUBROUTINE
c$$$
c$$$
c$$$      SUBROUTINE MINLV(V,N,I,VAL)
c$$$C********************************************
c$$$C     FINDS THE MINIMUM VALUE IN VECTOR V(N)
c$$$C     STORES THE MINIMUM VALUE IN VAL
c$$$C     AND STORES THE INDEX IN I
c$$$C********************************************
c$$$      INTEGER N,I
c$$$      DOUBLE PRECISION V,VAL
c$$$      DIMENSION V(N)
c$$$      INTEGER J
c$$$      VAL = 1.0D+10
c$$$      DO 10 J=1,N
c$$$         IF ( V(J).LT.VAL ) THEN
c$$$            I=J
c$$$            VAL = V(J)
c$$$         END IF
c$$$ 10   CONTINUE
c$$$      END SUBROUTINE
