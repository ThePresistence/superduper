      SUBROUTINE CGDMSS (FS,PS,IFS,IPS,JFS,JPS,N,MPDIAG,NITER,NCG,
     +                   PCGCRT,ICALL)
C     *
C     DENSITY MATRIX SEARCH USING SPARSE MATRICES.
C     *
C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
C     FS(*)     FOCK MATRIX (I), UNCHANGED.
C     PS(*)     DENSITY MATRIX (I,O).
C     IX(N+1)   POINTERS FOR SPARSE MATRICES IN CSR FORMAT (S).
C               X=FS,PS.
C     JX(*)     COLUMN INDICES FOR SPARSE MATRICES IN CSR FORMAT (S).
C               X=FS,PS.
C     N         NUMBER OF ORBITALS (I).
C     MPDIAG    TYPE OF DENSITY MATRIX (I,O).
C     NITER     NUMBER OF CURRENT SCF ITERATION (I).
C     NCG(*)    TOTAL NUMBER OF OPERATIONS DURING CG SEARCHES (I,O).
C               (1) CG CYCLES DONE.
C               (2) PURIFICATIONS DONE.
C               (3) MATRIX MULTIPLICATIONS DONE.
C               (4) LINEAR ROOTS ACCEPTED.
C               (5) QUADRATIC ROOTS SELECTED ON PHYSICAL GROUNDS.
C               (6) QUADRATIC ROOTS SELECTED BY COMPARISON OF F.
C               (7) CG SEARCHES CONVERGED.
C               (8) CG SEARCHES NOT CONVERGED.
C     PCGCRT    CRITERION ON CONVERGENCE OF DENSITY MATRIX (I).
C     ICALL     ERROR FLAG (I,O), NORMALLY NOT CHANGED.
C               =-1 ON OUTPUT: FATAL CGDMS ERROR.
C     *
      USE module3
      USE module1
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
      INTERFACE
      SUBROUTINE CGROOTS (PS,JPS,IPS,HS,JHS,IHS,FS,JFS,IFS,Q8,JQ8,
     +                    IQ8,N,XMU,B,C,D,XSQMAX,XI,FI,DI,TR1,TR3,
     +                    NITER,ICG,NCG,NPRINT,ICALL,FILT,
     +                    FILTALL,CUTM)
      USE module1
      USE module3
      IMPLICIT NONE
      INTEGER :: N,NCG(10),NITER,NPRINT,ICALL,ICG
      DOUBLE PRECISION :: B,XI,XMU,FI,TR1,CUTM,D,C,DI,TR3,XSQMAX
      DOUBLE PRECISION, DIMENSION (:), POINTER :: Q8,PS,FS,HS
      INTEGER, DIMENSION (:), POINTER :: IQ8,JQ8,JPS,IPS,JFS,IFS,
     +                                   JHS,IHS
      LOGICAL :: FILT,FILTALL
      END SUBROUTINE CGROOTS
C
      SUBROUTINE MATDEVS (A,JA,IA,B,JB,IB,N,MODE,DEVMAX,DEVRMS)
      USE module3
      IMPLICIT NONE
      INTEGER :: N,MODE
      DOUBLE PRECISION :: DEVMAX,DEVRMS
      DOUBLE PRECISION, DIMENSION (:), POINTER :: B,A
      INTEGER, DIMENSION (:), POINTER :: JB,IB,JA,IA
      END SUBROUTINE MATDEVS
C
      SUBROUTINE PURIFYP (PS,JPS,IPS,RS,JRS,IRS,N,DD1,DD2,NPR,
     +                    FILT,FILTALL,CUTM)
      USE module3
      IMPLICIT NONE
      INTEGER :: N,NPR
      DOUBLE PRECISION :: DD1,DD2,CUTM
      DOUBLE PRECISION, DIMENSION (:), POINTER :: RS,PS
      INTEGER, DIMENSION (:), POINTER :: IRS,JRS,IPS,JPS
      LOGICAL :: FILT,FILTALL
      END SUBROUTINE PURIFYP
      END INTERFACE
C
      LOGICAL :: GETPUR,FILT,FILTALL,PRECON
      COMMON
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./INOPT2/ IN2(300)
     ./INOPT4/ XN4(50)
     ./NBFILE/ NBF(20)
     ./ORBITS/ NUMB,NORBS,NMOS,NALPHA,NBETA
      INTEGER :: N,NITER,MAXCG,NPRINT,MAXPUR,MLROOT,N4,ICG,NALPHA,
     +           NBETA,IN2,MPURIF,MIDEMP,ICALL,I,NCG(10),IPUR,ie,
     +           NUMB,NORBS,NMOS,MCMAX,MCGUPD,IPSMAX,NBF,MPDIAG,J,
     +           NB6,MCGPRE,N2,MCUTM
      DOUBLE PRECISION :: ZERO,ONE,DMCMAX,DMCRMS,XEL,DOTGG,DPMAX,DPRMS,
     +                    DOTGQ,TR1,DOTQQ,DOTHQ,DOTHG,GAMMA,XMU,TRACEP,
     +                    PT5,FI,THREE,FOUR,PT25,PMCMAX,PCGCRT,TWO,B,D,
     +                    C,XSQMAX,PIDEMP,DIDMAX,DIDRMS,PPRMS,PASYMM,
     +                    PPMAX,PCGMAX,XI,TR3,CUTP,CUTM,XN4
      INTEGER, ALLOCATABLE ::IW(:)
      DOUBLE PRECISION, DIMENSION (:), POINTER :: Q8,GS,PS,HS,Q1,Q2,Q5,
     +                                            FS,WR
      INTEGER, DIMENSION (:), POINTER :: IQ8,JQ8,JGS,IGS,JPS,IPS,JHS,
     +                                   IHS,JQ1,IQ1,JQ2,IQ2,IW1,IW2,
     +                                   JFS,IFS
      SAVE Q1,JQ1,IQ1,Q2,JQ2,IQ2
C
C *** INITIALIZATION.
      XEL    = DBLE(NALPHA+NBETA)
      NB6    = NBF(6)
      NPRINT = IN2(72)
      MAXCG  = IN2(171)
      MAXPUR = IN2(172)
      MCMAX  = IN2(173)
      MIDEMP = IN2(174)
      MPURIF = IN2(175)
      MLROOT = IN2(176)
      MCGPRE = IN2(177)
      MCGUPD = IN2(178)
      MCUTM  = IN2(181)
      CUTP   = XN4(21)
C     CUTOFF FOR INTERMEDIATE MATRICES (CUTM).
C     FILT   = T : APPLY CUTOFF TO MATRIX PRODUCTS.
C     FILTALL= T : APPLY CUTOFF TO MATRIX SUMS.
      IF(MCUTM.GT.0) THEN
         CUTM = XN4(20)/DBLE(MCUTM)
         FILT = .TRUE.
         FILTALL = .TRUE.
      ELSE
         FILT = .FALSE.
         FILTALL = .FALSE.
      ENDIF
C     CRITERIA FOR PURIFICATION.
      IF(MCMAX.GT.0) THEN
         PMCMAX = 10.0D0**(-MCMAX)
      ELSE
         PMCMAX = PCGCRT/TWO
      ENDIF
      IF(MIDEMP.GT.0) THEN
         PIDEMP = 10.0D0**(-MIDEMP)
      ELSE
         PIDEMP = ONE
      ENDIF
C     CRITERIA FOR ROOT SEARCH.
      IF(MLROOT.GT.0) THEN
         XSQMAX = 10.0D0**(-MLROOT)
      ELSE
         XSQMAX = -ONE
      ENDIF
      DIDMAX = ZERO
C     FLAG FOR DIAGONAL PRECONDITIONING.
C     IF(MCGPRE.EQ.1 .AND. MCGUPD.NE.3 .AND. NITER.GE.2) THEN
      IF(MCGPRE.EQ.1 .AND. MCGUPD.NE.3) THEN
         PRECON = .TRUE.
      ELSE
         PRECON = .FALSE.
      ENDIF
C
C *** DEBUG PRINT.
      IF(NPRINT.GE.5) THEN
         WRITE(NB6,400)
         IF(NPRINT.GE.9) THEN
            WRITE(NB6,410)
            WRITE(NB6,420) (IN2(I),I=171,186)
            WRITE(NB6,430)
            WRITE(NB6,420) N,NITER,(NCG(I),I=1,8)
            WRITE(NB6,440) NITER
            CALL SPAPRT (FS,IFS,JFS,N)
            WRITE(NB6,450) NITER
            CALL SPAPRT (PS,IPS,JPS,N)
         ENDIF
      ENDIF
C
C *** EVALUATE INITIAL GRADIENT (GS).
      CALL DMSGRDS (FS,IFS,JFS,GS,IGS,JGS,PS,IPS,JPS,XMU,N,MPDIAG,
     +              NPRINT,NCG,FILT,FILTALL,CUTM,PRECON)
      MPDIAG = 0
      IF(NPRINT.GE.9) THEN
         WRITE(NB6,460) NITER
         CALL SPAPRT (GS,IGS,JGS,N)
      ENDIF
C
C *** DEFINE INITIAL STEEPEST-DESCENT SEARCH DIRECTION (HS=-GS).
      N4=IGS(N+1)-1
      ALLOCATE (HS(N4),JHS(N4),IHS(N+1),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','HS',N4,0)
      CALL copmatp (N,GS,JGS,IGS,HS,JHS,IHS,-1)
C
C *** INITIALIZE UNIT MATRIX (Q1) FOR DFP UPDATE.
C     NOTE: THIS OPTION IS NOT INCLUDED ELSEWHERE (CGDMS).
      IF(MCGUPD.EQ.3) THEN
         ALLOCATE (Q1(N),JQ1(N),IQ1(N+1),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGDMS','Q1',N,0)
         IQ1(N+1)=N+1
         DO 10 I=1,N
         IQ1(I)=I; JQ1(I)=I; Q1(I)=ONE
   10    CONTINUE
      ENDIF
C
C *** CONJUGATE GRADIENT ITERATIONS.
C     *
      DO 100 ICG=1,MAXCG
      NCG(1) = NCG(1)+1
C
C *** GET OPTIMUM COEFFICIENTS (B,C,D) FROM ANALYTIC LINE SEARCH.
      CALL DMSCOFS (FS,IFS,JFS,HS,IHS,JHS,GS,IGS,JGS,PS,IPS,JPS,N,
     +              NPRINT,B,C,D,FILT,CUTM)
      NCG(3) = NCG(3)+3
C     RETURN FOR CONVERGED INPUT DENSITY (E.G. H2 MOLECULE).
      IF(C.EQ.ZERO .AND. D.EQ.ZERO) THEN
         DEALLOCATE (HS,JHS,IHS,GS,JGS,IGS,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','HS',N4,1)
         NULLIFY (HS,JHS,IHS,GS,JGS,IGS)
         IF (ASSOCIATED(Q1)) THEN
            DEALLOCATE (Q1,JQ1,IQ1,STAT=ie)
               IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q1',N4,1)
            NULLIFY (Q1,JQ1,IQ1)
         ENDIF
         IF (ASSOCIATED(Q2)) THEN
            DEALLOCATE (Q2,JQ2,IQ2,STAT=ie)
               IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q2',N4,1)
            NULLIFY (Q2,JQ2,IQ2)
         ENDIF
         RETURN
      ENDIF
C
C *** SOLVE QUADRATIC EQUATION (SOLUTIONS X1 AND X2)
C     0   = B + 2*C*X + 3*D*X**2
C *** AND SELECT BEST UPDATE OF THE DENSITY MATRIX (Q5=P+XI*H).
C     ORIGINAL DENSITY MATRIX P, ORIGINAL SEARCH DIRECTION H.
C     UPDATED DENSITY MATRIX   , PURIFIED DENSITY MATRIX Q8.
      CALL CGROOTS (PS,JPS,IPS,HS,JHS,IHS,FS,JFS,IFS,Q8,JQ8,IQ8,N,
     +              XMU,B,C,D,XSQMAX,XI,FI,DIDMAX,TR1,TR3,NITER,
     +              ICG,NCG,NPRINT,ICALL,FILT,FILTALL,CUTM)
C     RETURN IN CASE OF ERROR.
      IF(ICALL.EQ.-1) THEN
         DEALLOCATE (HS,JHS,IHS,GS,JGS,IGS,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','HS',N4,1)
         NULLIFY (HS,JHS,IHS,GS,JGS,IGS)
         IF (ASSOCIATED(Q1)) THEN
            DEALLOCATE (Q1,JQ1,IQ1,STAT=ie)
               IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q1',N4,1)
            NULLIFY (Q1,JQ1,IQ1)
         ENDIF
         IF (ASSOCIATED(Q2)) THEN
            DEALLOCATE (Q2,JQ2,IQ2,STAT=ie)
               IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q2',N4,1)
            NULLIFY (Q2,JQ2,IQ2)
         ENDIF
         IF (ASSOCIATED(Q8)) THEN
            DEALLOCATE (Q8,JQ8,IQ8,STAT=ie)
               IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q8',N4,1)
            NULLIFY (Q8,JQ8,IQ8)
         ENDIF
         RETURN
      ENDIF
C *** COMPARE UPDATED MATRIX (Q8) AND ORIGINAL MATRIX (PS).
C     MAXIMUM AND RMS DEVIATIONS IN DIAGONAL DENSITY MATRIX ELEMENTS.
      CALL MATDEVS (Q8,JQ8,IQ8,PS,JPS,IPS,N,1,DPMAX,DPRMS)
C     DEBUG PRINT FOR UPDATED MATRIX (Q8).
      IF(NPRINT.GE.9) THEN
         WRITE(NB6,500)
         WRITE(NB6,510) NITER,ICG,XI,FI,DPMAX,DPRMS
         WRITE(NB6,520) NITER,ICG,TR1
         CALL SPAPRT (Q8,IQ8,JQ8,N)
      ENDIF
C
C *** McWEENY PURIFICATION TRANSFORMATIONS.
C     - DONE IF CG UPDATE INDICATES CONVERGENCE (DPMAX.LT.PCGCRT)
C     - OR   IF PURIFY OPTION IS TURNED ON (ICG.GE.MPURIF)
C     - OR   IF THIS IS THE LAST CG CYCLE (ICG.EQ.MAXCG).
C     MPURIF.GE.0 IS REQUIRED FOR THE LATTER TWO CASES.
C     MPURIF.LT.0 ALLOWS PURIFICATION ONLY CLOSE TO CONVERGENCE.
      GETPUR = DPMAX.LT.PCGCRT
      GETPUR = GETPUR .OR. (MPURIF.GE.0 .AND. ICG.GE.MPURIF)
      GETPUR = GETPUR .OR. (MPURIF.GE.0 .AND. ICG.EQ.MAXCG)
      IF(GETPUR) THEN
         DO 40 IPUR=1,MAXPUR
         NCG(2) = NCG(2)+1
         NCG(3) = NCG(3)+2
C        COPY UPDATED DENSITY MATRIX FROM CGROOTS (Q8->Q5).
         N4=IQ8(N+1)-1
         ALLOCATE (Q5(N4),IW1(N4),IW2(N+1),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q5',N4,0)
         CALL copmatp (N,Q8,JQ8,IQ8,Q5,IW1,IW2,1)
         DEALLOCATE (Q8,JQ8,IQ8,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q8',N4,1)
         NULLIFY (Q8,JQ8,IQ8)
C        MCWEENY PURIFICATION (Q5->Q8).
         CALL PURIFYP (Q5,IW1,IW2,Q8,JQ8,IQ8,N,DIDMAX,DIDRMS,NPRINT,
     +                 FILT,FILTALL,CUTM)
C
C        CHECK WHETHER PURIFIED DENSITY (Q8) IS PHYSICALLY MEANINGFUL.
C        IF DIAGONAL ELEMENTS OF Q8 ARE OUTSIDE THE ALLOWED RANGE, 
C        THE PREVIOUS DENSITY MATRIX (Q5) IS USED AS AN UPDATE.
C        NOTE: THIS TEST IS NOT INCLUDED ELSEWHERE (CGDMS).
         DO 30 I=1,N
         DO 20 J=IQ8(I),IQ8(I+1)-1
         IF(JQ8(J).EQ.I) THEN
            IF(Q8(J).LT.-ONE .OR. Q8(J).GT.TWO) THEN
               N4=IW2(N+1)-1
               DEALLOCATE (Q8,JQ8,STAT=ie)
                  IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q8',N4,1)
               NULLIFY (Q8,JQ8)
               ALLOCATE (Q8(N4),JQ8(N4),STAT=ie)
                  IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q8',N4,0)
               CALL copmatp (N,Q5,IW1,IW2,Q8,JQ8,IQ8,1)
               DEALLOCATE (Q5,IW1,IW2,STAT=ie)
                  IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q5',N4,1)
               NULLIFY (Q5,IW1,IW2)
               GO TO 60 
            ENDIF
            GO TO 30
         ENDIF
   20    CONTINUE
   30    CONTINUE
C
C        COMPARE PREVIOUS (Q5) AND PURIFIED (Q8) MATRIX.
         CALL MATDEVS (Q5,IW1,IW2,Q8,JQ8,IQ8,N,1,DMCMAX,DMCRMS)
         DEALLOCATE (Q5,IW1,IW2,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q5',N4,1)
         NULLIFY (Q5,IW1,IW2)
C        DEBUG PRINT.
         IF(NPRINT.GE.7) THEN
            IF(IPUR.EQ.1) WRITE(NB6,530)
            WRITE(NB6,540) NITER,ICG,IPUR,DIDMAX,DIDRMS,DMCMAX,DMCRMS
         ENDIF
         IF(DMCMAX.LT.PMCMAX .AND. DIDMAX.LT.PIDEMP) GO TO 50
   40    CONTINUE
   50    CONTINUE
         MPURIF = MAX(MPURIF,0)
      ENDIF
   60 CONTINUE
C
C *** COMPARE PURIFIED MATRIX (Q8) AND ORIGINAL MATRIX (PS).
C     MAXIMUM AND RMS DEVIATIONS IN DIAGONAL DENSITY MATRIX ELEMENTS.
      CALL MATDEVS (Q8,JQ8,IQ8,PS,JPS,IPS,N,1,PPMAX,PPRMS)
C     DEBUG PRINT FOR PURIFIED MATRIX (Q8).
      IF(NPRINT.GE.7) THEN
         ALLOCATE (WR(N),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','WR',N,0)
         CALL TRACENP (N,Q8,JQ8,IQ8,FS,JFS,IFS,WR,FI)
         DEALLOCATE (WR,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','WR',N,1)
         NULLIFY (WR)
         FI = FI + XMU*(TR1-XEL)
         WRITE(NB6,550)
         WRITE(NB6,510) NITER,ICG,XI,FI,PPMAX,PPRMS
         IF(NPRINT.GE.9) THEN
            CALL TRACEAP (N,Q8,JQ8,IQ8,TR3)
            WRITE(NB6,560) NITER,ICG,TR3
            CALL SPAPRT (Q8,IQ8,JQ8,N)
         ENDIF
      ENDIF
C
C *** DIFFERENCE DENSITY (Q2=Q8-PS) FOR DFP SEARCH DIRECTION UPDATE.
C     NOTE: THIS OPTION IS NOT INCLUDED ELSEWHERE (CGDMS).
      IF(MCGUPD.EQ.3) THEN
C        SYMBOLIC ADDITION TO GET NUMBER OF ELEMENTS (N4).
         ALLOCATE (IW(N),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','IW',N,0)
         CALL aplbdgp (N,JQ8,IQ8,JPS,IPS,N4,IW)
C        DIFFERENCE (Q2) BETWEEN PURIFIED (Q8) AND OLD (PS) DENSITY.
         ALLOCATE (Q2(N4),JQ2(N4),IQ2(N+1),STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q2',N4,0)
         N2=IPS(N+1)-1
         DO 32 I=1,N2
         PS(I)=-PS(I)
 32      CONTINUE
         CALL aplbp (N,Q8,JQ8,IQ8,PS,JPS,IPS,Q2,JQ2,IQ2,N4,IW,ie)
            IF(ie.NE.0) CALL XERSPA (ie,'aplb','CGDMSS',1)
C        REMOVE SMALL ELEMENTS (Q2).
         IF(FILTALL) THEN
            CALL filterp (N,1,CUTM,Q2,JQ2,IQ2,Q2,JQ2,IQ2,N4,ie)
              IF(ie.ne.0) CALL XERSPA (ie,'filter','CGDMSS',4)
         ENDIF
         DEALLOCATE (IW,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','IW',N,1)
      ENDIF
C
C *** DEFINE NEW DENSITY MATRIX (PS=Q8).
      N4=IQ8(N+1)-1
      DEALLOCATE (PS,JPS,STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','PS',N4,1)
      NULLIFY (PS,JPS)
      IF(NPRINT.GE.7) THEN
         WRITE(NB6,740) 100.D0*(ONE-DBLE(N4)/DBLE(N**2))
      ENDIF
      ALLOCATE (PS(N4),JPS(N4),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','PS',N4,0)
C     COPY THE DENSITY MATRIX (Q8->PS).
      CALL copmatp (N,Q8,JQ8,IQ8,PS,JPS,IPS,1)
      PCGMAX = PPMAX
      DEALLOCATE (Q8,JQ8,IQ8,STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q8',N4,1)
      NULLIFY (Q8,JQ8,IQ8)
C
C *** CHECK FOR CONVERGENCE.
      IF(PCGMAX.LT.PCGCRT) THEN
         NCG(7) = NCG(7)+1
         IF(NPRINT.GE.7) THEN
            WRITE(NB6,570) PCGMAX,PCGCRT
         ENDIF
         GO TO 110
      ELSE IF(ICG.EQ.MAXCG) THEN
         NCG(8) = NCG(8)+1
         IF(NPRINT.GE.7) THEN
            WRITE(NB6,580) ICG,PCGMAX,PCGCRT
         ENDIF
         GO TO 110
      ENDIF
C
C *** PREPARE NEXT CONJUGATE GRADIENT CYCLE.
C     EVALUATE GRADIENT FOR NEW DENSITY MATRIX (Q8).
      CALL DMSGRDS (FS,IFS,JFS,Q8,IQ8,JQ8,PS,IPS,JPS,XMU,N,MPDIAG,
     +              NPRINT,NCG,FILT,FILTALL,CUTM,.FALSE.)
C
C *** CONJUGATE GRADIENT SEARCH DIRECTION UPDATE (HS).
      IF(MCGUPD.NE.3) THEN
         CALL CGUPDS (Q8,JQ8,IQ8,HS,JHS,IHS,GS,JGS,IGS,N,FILTALL,CUTM,
     +                NITER,ICG,NPRINT)
      ELSE
C *** DAVIDON-FLETCHER-POWELL SEARCH DIRECTION UPDATE (HS).
C     NOTE: THIS OPTION IS NOT INCLUDED ELSEWHERE (CGDMS).
         DEALLOCATE (HS,JHS,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','HS',N,1)
         NULLIFY (HS,JHS)
         CALL DFPUPD (HS,JHS,IHS,GS,JGS,IGS,Q1,JQ1,IQ1,Q2,JQ2,IQ2,Q8,
     +                JQ8,IQ8,N,FILT,FILTALL,CUTM,ICG,NPRINT)
      ENDIF
  100 CONTINUE
  110 CONTINUE
C     *
C *** SYMMETRIZE THE DENSITY MATRIX FROM THE CG SEARCH.
C     THE FINITE NUMERICAL ACCURACY OF THE COMPUTATION CAUSES THE
C     DENSITY MATRIX TO BECOME UNSYMMETRIC (VERY SLIGHTLY) DURING
C     THE CG SEARCH. THIS NUMERICAL NOISE CAN ACCUMULATE DURING THE
C     SCF ITERATIONS AND MAY LEAD TO CONVERGENCE PROBLEMS OR OTHER
C     NUMERICAL TROUBLE. SYMMETRIZATION AVOIDS THESE PROBLEMS.
      IF (ASSOCIATED(Q1)) THEN
         DEALLOCATE (Q1,JQ1,IQ1,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q1',N4,1)
         NULLIFY (Q1,JQ1,IQ1)
      ENDIF
      IF (ASSOCIATED(Q2)) THEN
         DEALLOCATE (Q2,JQ2,IQ2,STAT=ie)
            IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','Q2',N4,1)
         NULLIFY (Q2,JQ2,IQ2)
      ENDIF
      N4=IPS(N+1)-1
      DEALLOCATE (GS,JGS,HS,JHS,STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','GS',N4,1)
      NULLIFY (GS,JGS,HS,JHS)
      ALLOCATE (GS(N4),JGS(N4),HS(N4),JHS(N4),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','GS',N4,0)
C     GENERATE TRANSPOSE (GS) OF DENSITY MATRIX (PS).
      CALL csrcscp (N,PS,JPS,IPS,GS,JGS,IGS)
      CALL copmatp (N,PS,JPS,IPS,HS,JHS,IHS,1)
C     SYMBOLIC ADDITION TO GET NUMBER OF ELEMENTS (N4).
      ALLOCATE (IW(N+1),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','IW',N,0)
      CALL aplbdgp (N,JHS,IHS,JGS,IGS,N4,IW)
      DEALLOCATE (PS,JPS,STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','PS',N,1)
      NULLIFY (PS,JPS)
      ALLOCATE (PS(N4),JPS(N4),STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','PS',N4,0)
C     DEBUG PRINT.
      IF(NPRINT.GE.5) THEN
         DO 120 I=1,IGS(N+1)-1
         GS(I) = -GS(I)
  120    CONTINUE
         CALL aplbp (N,HS,JHS,IHS,GS,JGS,IGS,PS,JPS,IPS,N4,IW,ie)
            IF(ie.NE.0) CALL XERSPA (ie,'aplb','CGDMSS',1)
         PASYMM = ZERO
         IPSMAX = IDAMAXP(N4,PS)
         PASYMM = ABS(PS(IPSMAX))
         CALL TRACEAP (N,HS,JHS,IHS,TRACEP)
         WRITE(NB6,610) NITER,PASYMM
         WRITE(NB6,620) TRACEP
         WRITE(NB6,630) DIDMAX
         WRITE(NB6,400)
         DO 130 I=1,IGS(N+1)-1
         GS(I) = -GS(I)
  130    CONTINUE
      ENDIF
C     ADD DENSITY MATRIX AND ITS TRANSPOSE TO YIELD NEW MATRIX (PS).
      CALL aplbp (N,HS,JHS,IHS,GS,JGS,IGS,PS,JPS,IPS,N4,IW,ie)
         IF(ie.NE.0) CALL XERSPA (ie,'aplb','CGDMSS',1)
      DEALLOCATE (IW,HS,JHS,GS,JGS,IGS,IHS,STAT=ie)
         IF(ie.NE.0) CALL XERALL (ie,'CGDMSS','IW',N,1)
      NULLIFY (HS,JHS,GS,JGS,IGS,IHS)
C     SCALE BY 0.5 TO GET SYMMETRIZED DENSITY MATRIX (PS).
      DO 140 I=1,IPS(N+1)-1
      PS(I) = PT5*PS(I)
  140 CONTINUE
C     REMOVE SMALL ELEMENTS.
      IF(FILT) THEN
         CALL filterp (N,1,CUTP,PS,JPS,IPS,PS,JPS,IPS,N4,ie)
           IF(ie.ne.0) CALL XERSPA (ie,'filter','CGDMSS',4)
      ENDIF
      RETURN
  400 FORMAT(   1X)
  410 FORMAT(   1X,'CGDMS :  MAXCP MAXPUR  IDMSP  IDEMP MPURIF',
     1          1X,        'MLROOT MCGPRE MCGUPD MPSCAL MCUTAU',
     2          1X,        ' MCUTM  MCUTF  MCUTP  MCUT1  MCUT2',
     3          1X,        ' MCUTR')
  420 FORMAT(   1X,'CGDMS :',16I7)
  430 FORMAT(   1X,'CGDMS :  NORBS  NITER NCG(1) NCG(2) NCG(3)',
     1          1X,        'NGC(4) NCG(5) NCG(6) NCG(7) NCG(8)')
  440 FORMAT(// 1X,'CGDMS : INPUT FOCK MATRIX, NITER =',I3/)
  450 FORMAT(// 1X,'CGDMS : INPUT DENSITY MATRIX, NITER =',I3/)
  460 FORMAT(// 1X,'CGDMS : INITIAL GRADIENT MATRIX, NITER =',I3/)
  500 FORMAT(/  1X,'CGDMS : NITER,ICG,XI,FI,DPMAX,DPRMS')
  510 FORMAT(   1X,'CGDMS : ',I5,I3,3X,5G20.10)
  520 FORMAT(// 1X,'CGDMS : DENSITY MATRIX AFTER CG UPDATE.',
     1       /  1X,'CGDMS : NITER =',I3,'   ICG =',I3,'   TRACE =',
     2                      G20.10/)
  530 FORMAT(/  1X,'CGDMS : NITER,ICG,IPUR,DIDMAX,DIDRMS,DMCMAX,DMCRMS')
  540 FORMAT(   1X,'CGDMS : ',I5,2I3,4G20.10)
  550 FORMAT(/  1X,'CGDMS : NITER,ICG,XI,FI,PPMAX,PPRMS')
  560 FORMAT(// 1X,'CGDMS : DENSITY MATRIX AFTER PURIFICATION.',
     1       /  1X,'CGDMS : NITER =',I3,'   ICG =',I3,'   TRACE =',
     2                      G20.10/)
  570 FORMAT(/  1X,'CGDMS : DENSITY MATRIX CONVERGED IN CG SEARCH.',
     1       /  1X,'CGDMS : PCGMAX =',G20.10,
     2       /  1X,'CGDMS : PCGCRT =',G20.10,/)
  580 FORMAT(/  1X,'CGDMS : DENSITY MATRIX NOT CONVERGED IN CG SEARCH.',
     1       /  1X,'CGDMS : NUMBER OF CG CYCLES =',I3,
     2       /  1X,'CGDMS : PCGMAX =',G20.10,
     3       /  1X,'CGDMS : PCGCRT =',G20.10,/)
  610 FORMAT(   1X,'CGDMS : CG CYCLES DONE. NITER =',I3,
     1       /  1X,'CGDMS : DENSITY MATRIX WILL BE SYMMETRIZED',
     2       /  1X,'CGDMS : MAXIMUM ASYMMETRY OF DENSITY MATRIX =',
     3                      G20.10)
  620 FORMAT(   1X,'CGDMS : TRACE OF THE CURRENT DENSITY MATRIX =',
     1                      G20.10)
  630 FORMAT(   1X,'CGDMS : MAXIMUM DEVIATION FROM IDEMPOTENCY  =',
     1                      G20.10)
  730 FORMAT(   1X,'CGDMS : UNPHYSICAL DIAGONAL DENSITY MATRIX ELEMENT',
     1       /  1X,'CGDMS : I =',I6,5X,'VALUE =',G20.10)
  740 FORMAT(   1X,'CGDMS : SPARSITY OF INTERMEDIATE DENSITY MATRIX =',
     1                      F10.5,' %')
      END
