      SUBROUTINE ee132 (A,LDA,P,F1,N1,R,F2,N2,H,ICNTR,ESCAL,NDEN,I1n)
      USE LIMIT, ONLY: LM1, LMX, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ENERGT/ EE,ENUCLR,EAT,ATHEAT
     ./LIMITS/ LM2,LM3,LM4,LM6,LM7,LM8,LM9
     ./LMGRAD/ LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8
     ./LMSCF / LS1,LS2,LS3,LS4,LS5,LS6,LS7,LS8,LS9
     ./LMUHF / LU1,LU2,LU3,LU4,LU5,LU6,LU7,LU8
      DIMENSION A(LDA),P(LM3*LM3,N1),F1(LM3*LM3,N1)
      DIMENSION H(LM3*LM3),R(LM3*LM3,N2),F2(LM3*LM3,N2)

C *** SIMPLIFIED CODE FOR GRADIENT CALCULATIONS WITH DISPLACED ATOM.
C *** INTEGRAL CALCULATION FOR SCF.
      CALL HCORE2 (A(LS7),A(LS9),A(LG1),A(LG2),A(LG3),A(LG4),
     1             LM2,LM4,LM6,LM9,ICNTR)

C *** CORE REPULSION ENERGY.
      ESCAL = ENUCLR

C *** TWO-ELECTRON CONTRIBUTIONS TO FOCK MATRIX FOR GIVEN DENSITY.
      CALL FOCK132(F1,F1(1,N1*(NDEN-1)+1),P,P(1,N1*(NDEN-1)+1),N1,
     $             F2,F2(1,N2*(NDEN-1)+1),R,R(1,N2*(NDEN-1)+1),N2,
     $             A(LS9),LM3,LM6,NDEN,I1n,ICNTR)

      CALL SQUARE(A(LS7),H,LM3,LM3,LM4)

      RETURN
      END


      SUBROUTINE FOCK132(F1a,F1b,Qa,Qb,NQ,F2a,F2b,Ra,Rb,NR,W,LM3,LM6,
     $                  NDen,I1n,ICNTR)
!     *
!     BUILD CIS FOCK-LIKE MATRIX.
!     *
!     REFERENCE OCCUPATION NUMBERS OF DOUBLY OCCUPIED (CLOSED-SHELL)
!     ORBITALS ARE SET TO TWO, AND ALL OTHERS ARE SET TO ZERO.
!     *
!     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
!     Ra/Rb(LM3,LM3)  DESITY MATRIX(I).
!     W(LM6,LM6)      TWO-ELECTRON INTEGRALS(I).
!     F2a/F2b(LM3,LM3)      CIS FOCK-LIKE MATRIX(O).

      USE LIMIT, ONLY: LM1, LMI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      COMMON /CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
      COMMON /FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
      COMMON /INDEXW/ NW(LM1)
      REAL*8 Qa(LM3,LM3,*),Qb(LM3,LM3,*),F1a(LM3,LM3,*),F1b(LM3,LM3,*),
     $       Ra(LM3,LM3,*),Rb(LM3,LM3,*),F2a(LM3,LM3,*),F2b(LM3,LM3,*),
     $       W(LM6,LM6)
      INTEGER NL,LM3,LM6,NDen

      CALL VecInit(F1a,LM3*LM3*NQ,0d0)
      CALL VecInit(F2a,LM3*LM3*NR,0d0)
      IF(NDen.EQ.2) THEN
        CALL VecInit(F1b,LM3*LM3*NQ,0d0)
        CALL VecInit(F2b,LM3*LM3*NR,0d0)
      ENDIF
!     IP1(I).GE.IP2(I)
      DO II=1,NUMAT
      NA = NW(II) - 1
      IA = NFIRST(II)
      IB = NLAST(II)
      DO IMU = IA,IB
      DO INU = IA,IMU  
      NA = NA + 1
         DO JJ=1,NUMAT
         IF((II.EQ.ICNTR .OR. JJ.EQ.ICNTR) .AND. II.NE.JJ) THEN
         !IF((II.EQ.ICNTR .OR. JJ.EQ.ICNTR)) THEN
         !write(6,*) "atoms:",II,JJ
         NB = NW(JJ) - 1
         JA = NFIRST(JJ)
         JB = NLAST(JJ)
         DO IRO = JA,JB
         DO ISI = JA,IRO
         NB = NB + 1
         IF(NB.GE.NA)  THEN   
            Xij = W(NA,NB)
            Yij = -W(NA,NB)
!           COULOMB PART.
!           (MU,NU | RO,SI) + (MU,NU | SI,RO)  
!           (NU,MU | RO,SI) + (NU,MU | SI,RO)
            DO M=1,NR
            IF(I1n.EQ.1) THEN
               Rho1 = Ra(IRO,ISI,M) + Rb(IRO,ISI,M)
               IF(IRO.NE.ISI)  Rho1 = Rho1 + Ra(ISI,IRO,M) 
     $                                     + Rb(ISI,IRO,M)
               Rho2 = Ra(IMU,INU,M) + Rb(IMU,INU,M)
               IF(IMU.NE.INU)  Rho2 = Rho2 + Ra(INU,IMU,M) 
     $                                     + Rb(INU,IMU,M)
             
               F2a(IMU,INU,M) = F2a(IMU,INU,M) + Rho1*Xij
               IF(IMU.NE.INU) F2a(INU,IMU,M) = F2a(INU,IMU,M) + Rho1*Xij
               IF(NA.NE.NB) THEN
                 F2a(IRO,ISI,M) = F2a(IRO,ISI,M) + Rho2*Xij
                 IF(IRO.NE.ISI)  
     $             F2a(ISI,IRO,M) = F2a(ISI,IRO,M) + Rho2*Xij
               ENDIF
            ENDIF

!           EXCHANGE PART.
!           (MU,NU | RO,SI)
            F2a(IMU,IRO,M) = F2a(IMU,IRO,M) + Ra(INU,ISI,M)*Yij
!            write(6,*) M,":",IMU,IRO,"--",INU,ISI
            IF(IMU.NE.IRO .OR. INU.NE.ISI)   
     $        F2a(IRO,IMU,M) = F2a(IRO,IMU,M) + Ra(ISI,INU,M)*Yij
!            IF(IMU.NE.IRO .OR. INU.NE.ISI) 
!     $        write(6,*) M,":",IRO,IMU,"--",ISI,INU 
            IF(IRO.NE.ISI) THEN
!             (MU,NU | SI,RO)
              F2a(IMU,ISI,M) = F2a(IMU,ISI,M) + Ra(INU,IRO,M)*Yij
!              write(6,*) M,":",IMU,ISI,"--",INU,IRO
              IF(ISI.NE.IMU .OR. INU.NE.IRO) 
     $          F2a(ISI,IMU,M) = F2a(ISI,IMU,M) + Ra(IRO,INU,M)*Yij
!              IF(IMU.NE.ISI .OR. INU.NE.IRO) 
!     $          write(6,*) M,":",ISI,IMU,"--",IRO,INU 
            ENDIF
            IF(IMU.NE.INU ) THEN
!             (NU,MU | RO,SI)
              IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
              F2a(INU,IRO,M) = F2a(INU,IRO,M) + Ra(IMU,ISI,M)*Yij
!              write(6,*) M,":",INU,IRO,"--",IMU,ISI
              IF(INU.NE.IRO .OR. IMU.NE.ISI) 
     $          F2a(IRO,INU,M) = F2a(IRO,INU,M) + Ra(ISI,IMU,M)*Yij
!              IF(INU.NE.IRO .OR. IMU.NE.ISI) 
!     $          write(6,*) M,":",IRO,INU,"--",ISI,IMU
              ENDIF
              IF(IRO.NE.ISI) THEN
!               (NU,MU | SI,RO)
                F2a(INU,ISI,M) = F2a(INU,ISI,M) + Ra(IMU,IRO,M)*Yij
!                write(6,*) M,":",INU,ISI,"--",IMU,IRO
                IF(INU.NE.ISI .OR. IMU.NE.IRO) 
     $            F2a(ISI,INU,M) = F2a(ISI,INU,M) + Ra(IRO,IMU,M)*Yij
!                IF(INU.NE.ISI .OR. IMU.NE.IRO) 
!     $            write(6,*) M,":",ISI,INU,"--",IRO,IMU 
              ENDIF
            ENDIF
            

            IF(NDen.EQ.2) THEN
!             COULOMB PART.
!             (MU,NU | RO,SI) + (MU,NU | SI,RO)  
!             (NU,MU | RO,SI) + (NU,MU | SI,RO)
              IF(I1n.EQ.1) THEN
              F2b(IMU,INU,M) = F2b(IMU,INU,M) + Rho1*Xij
              IF(IMU.NE.INU)  
     $          F2b(INU,IMU,M) = F2b(INU,IMU,M) + Rho1*Xij
              IF(NA.NE.NB) THEN
                F2b(IRO,ISI,M) = F2b(IRO,ISI,M) + Rho2*Xij
                IF(IRO.NE.ISI)  
     $            F2b(ISI,IRO,M) = F2b(ISI,IRO,M) + Rho2*Xij
              ENDIF
              ENDIF

!             EXCHANGE PART.
!             (MU,NU | RO,SI)
              F2b(IMU,IRO,M) = F2b(IMU,IRO,M) + Rb(INU,ISI,M)*Yij
              IF(IMU.NE.IRO .OR. INU.NE.ISI) 
     $          F2b(IRO,IMU,M) = F2b(IRO,IMU,M) + Rb(ISI,INU,M)*Yij
              IF(IRO.NE.ISI) THEN
!               (MU,NU | SI,RO)
                F2b(IMU,ISI,M) = F2b(IMU,ISI,M) + Rb(INU,IRO,M)*Yij
                IF(ISI.NE.IMU .OR. INU.NE.IRO) 
     $            F2b(ISI,IMU,M) = F2b(ISI,IMU,M) + Rb(IRO,INU,M)*Yij
              ENDIF
              IF(IMU.NE.INU ) THEN
!               (NU,MU | RO,SI)
                IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
                F2b(INU,IRO,M) = F2b(INU,IRO,M) + Rb(IMU,ISI,M)*Yij
                IF(INU.NE.IRO .OR. IMU.NE.ISI) 
     $            F2b(IRO,INU,M) = F2b(IRO,INU,M) + Rb(ISI,IMU,M)*Yij
                ENDIF
                IF(IRO.NE.ISI) THEN
!                 (NU,MU | SI,RO)
                  F2b(INU,ISI,M) = F2b(INU,ISI,M) + Rb(IMU,IRO,M)*Yij
                  IF(INU.NE.ISI .OR. IMU.NE.IRO) 
     $              F2b(ISI,INU,M) = F2b(ISI,INU,M) + Rb(IRO,IMU,M)*Yij
                ENDIF
              ENDIF

            ENDIF ! Beta Spin
            ENDDO ! M

            DO M=1,NQ
!            IF(I1n.EQ.1) THEN
               Rho1 = Qa(IRO,ISI,M) + Qb(IRO,ISI,M)
               IF(IRO.NE.ISI)  Rho1 = Rho1 + Qa(ISI,IRO,M) 
     $                                     + Qb(ISI,IRO,M)
               Rho2 = Qa(IMU,INU,M) + Qb(IMU,INU,M)
               IF(IMU.NE.INU)  Rho2 = Rho2 + Qa(INU,IMU,M) 
     $                                     + Qb(INU,IMU,M)
             
!               write(6,*) IMU,INU,"--",IRO,ISI
!               if(IMU.NE.INU) write(6,*) INU,IMU,"--",IRO,ISI
!               if(NA.NE.NB) THEN
!                 write(6,*) IRO,ISI,"--",IMU,INU
!                 if(IRO.NE.ISI) write(6,*) ISI,IRO,"--",IMU,INU
!               endif
               F1a(IMU,INU,M) = F1a(IMU,INU,M) + Rho1*Xij
               IF(IMU.NE.INU) F1a(INU,IMU,M) = F1a(INU,IMU,M) + Rho1*Xij
               IF(NA.NE.NB) THEN
                 F1a(IRO,ISI,M) = F1a(IRO,ISI,M) + Rho2*Xij
                 IF(IRO.NE.ISI)  
     $              F1a(ISI,IRO,M) = F1a(ISI,IRO,M) + Rho2*Xij
               ENDIF
!            ENDIF

!            Yij = 0d0
!           EXCHANGE PART.
!           (MU,NU | RO,SI)
            F1a(IMU,IRO,M) = F1a(IMU,IRO,M) + Qa(INU,ISI,M)*Yij
!            write(6,*) M,":",IMU,IRO,"--",INU,ISI
            IF(IMU.NE.IRO .OR. INU.NE.ISI)   
     $        F1a(IRO,IMU,M) = F1a(IRO,IMU,M) + Qa(ISI,INU,M)*Yij
!            IF(IMU.NE.IRO .OR. INU.NE.ISI) 
!     $        write(6,*) M,":",IRO,IMU,"--",ISI,INU 
            IF(IRO.NE.ISI) THEN
!             (MU,NU | SI,RO)
              F1a(IMU,ISI,M) = F1a(IMU,ISI,M) + Qa(INU,IRO,M)*Yij
!              write(6,*) M,":",IMU,ISI,"--",INU,IRO
              IF(ISI.NE.IMU .OR. INU.NE.IRO) 
     $          F1a(ISI,IMU,M) = F1a(ISI,IMU,M) + Qa(IRO,INU,M)*Yij
!              IF(IMU.NE.ISI .OR. INU.NE.IRO) 
!     $          write(6,*) M,":",ISI,IMU,"--",IRO,INU 
            ENDIF
            IF(IMU.NE.INU ) THEN
!             (NU,MU | RO,SI)
              IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
              F1a(INU,IRO,M) = F1a(INU,IRO,M) + Qa(IMU,ISI,M)*Yij
!              write(6,*) M,":",INU,IRO,"--",IMU,ISI
              IF(INU.NE.IRO .OR. IMU.NE.ISI) 
     $          F1a(IRO,INU,M) = F1a(IRO,INU,M) + Qa(ISI,IMU,M)*Yij
!              IF(INU.NE.IRO .OR. IMU.NE.ISI) 
!     $          write(6,*) M,":",IRO,INU,"--",ISI,IMU
              ENDIF
              IF(IRO.NE.ISI) THEN
!               (NU,MU | SI,RO)
                F1a(INU,ISI,M) = F1a(INU,ISI,M) + Qa(IMU,IRO,M)*Yij
!                write(6,*) M,":",INU,ISI,"--",IMU,IRO
                IF(INU.NE.ISI .OR. IMU.NE.IRO) 
     $            F1a(ISI,INU,M) = F1a(ISI,INU,M) + Qa(IRO,IMU,M)*Yij
!                IF(INU.NE.ISI .OR. IMU.NE.IRO) 
!     $            write(6,*) M,":",ISI,INU,"--",IRO,IMU 
              ENDIF
            ENDIF

            IF(NDen.EQ.2) THEN
!             COULOMB PART.
!             (MU,NU | RO,SI) + (MU,NU | SI,RO)  
!             (NU,MU | RO,SI) + (NU,MU | SI,RO)
              F1b(IMU,INU,M) = F1b(IMU,INU,M) + Rho1*Xij
              IF(IMU.NE.INU)  
     $          F1b(INU,IMU,M) = F1b(INU,IMU,M) + Rho1*Xij
              IF(NA.NE.NB) THEN
                F1b(IRO,ISI,M) = F1b(IRO,ISI,M) + Rho2*Xij
                IF(IRO.NE.ISI)  
     $            F1b(ISI,IRO,M) = F1b(ISI,IRO,M) + Rho2*Xij
              ENDIF

!             EXCHANGE PART.
!             (MU,NU | RO,SI)
              F1b(IMU,IRO,M) = F1b(IMU,IRO,M) + Qb(INU,ISI,M)*Yij
              IF(IMU.NE.IRO .OR. INU.NE.ISI) 
     $          F1b(IRO,IMU,M) = F1b(IRO,IMU,M) + Qb(ISI,INU,M)*Yij
              IF(IRO.NE.ISI) THEN
!               (MU,NU | SI,RO)
                F1b(IMU,ISI,M) = F1b(IMU,ISI,M) + Qb(INU,IRO,M)*Yij
                IF(ISI.NE.IMU .OR. INU.NE.IRO) 
     $            F1b(ISI,IMU,M) = F1b(ISI,IMU,M) + Qb(IRO,INU,M)*Yij
              ENDIF
              IF(IMU.NE.INU ) THEN
!               (NU,MU | RO,SI)
                IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
                F1b(INU,IRO,M) = F1b(INU,IRO,M) + Qb(IMU,ISI,M)*Yij
                IF(INU.NE.IRO .OR. IMU.NE.ISI) 
     $            F1b(IRO,INU,M) = F1b(IRO,INU,M) + Qb(ISI,IMU,M)*Yij
                ENDIF
                IF(IRO.NE.ISI) THEN
!                 (NU,MU | SI,RO)
                  F1b(INU,ISI,M) = F1b(INU,ISI,M) + Qb(IMU,IRO,M)*Yij
                  IF(INU.NE.ISI .OR. IMU.NE.IRO) 
     $              F1b(ISI,INU,M) = F1b(ISI,INU,M) + Qb(IRO,IMU,M)*Yij
                ENDIF
              ENDIF

            ENDIF ! Beta Spin
            ENDDO ! NQ
          ENDIF ! NA,NB

         ENDDO !JB
         ENDDO !JA
         ENDIF !ICNTR
         ENDDO !JJ
      ENDDO !IB
      ENDDO !IA
      ENDDO !II

      RETURN
      END

      SUBROUTINE makejkx (A,LDA,Pa,Pb,Ja,Jb,Ka,Kb,N1,H,ICNTR,ESCAL,
     $                    NDEN)
      USE LIMIT, ONLY: LM1, LMX, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ENERGT/ EE,ENUCLR,EAT,ATHEAT
     ./LIMITS/ LM2,LM3,LM4,LM6,LM7,LM8,LM9
     ./LMGRAD/ LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8
     ./LMSCF / LS1,LS2,LS3,LS4,LS5,LS6,LS7,LS8,LS9
     ./LMUHF / LU1,LU2,LU3,LU4,LU5,LU6,LU7,LU8
      DIMENSION A(LDA),Pa(LM3*LM3,N1),Pb(LM3*LM3,N1)
      DIMENSION Ja(LM3*LM3,N1),Jb(LM3*LM3,N1)
      DIMENSION Ka(LM3*LM3,N1),Kb(LM3*LM3,N1)
      DIMENSION H(LM3*LM3)

C *** SIMPLIFIED CODE FOR GRADIENT CALCULATIONS WITH DISPLACED ATOM.
C *** INTEGRAL CALCULATION FOR SCF.
      CALL HCORE2 (A(LS7),A(LS9),A(LG1),A(LG2),A(LG3),A(LG4),
     1             LM2,LM4,LM6,LM9,ICNTR)

C *** CORE REPULSION ENERGY.
      ESCAL = ENUCLR

C *** TWO-ELECTRON CONTRIBUTIONS TO FOCK MATRIX FOR GIVEN DENSITY.
      CALL FOCK132a(Pa,Pb,Ja,Jb,Ka,Kb,N1,A(LS9),LM3,LM6,NDEN,ICNTR)

      CALL SQUARE(A(LS7),H,LM3,LM3,LM4)

      RETURN
      END


      SUBROUTINE FOCK132a(Qa,Qb,Ja,Jb,Ka,Kb,NQ,W,LM3,LM6,NDen,ICNTR)
!     *
!     BUILD CIS FOCK-LIKE MATRIX.
!     *
!     REFERENCE OCCUPATION NUMBERS OF DOUBLY OCCUPIED (CLOSED-SHELL)
!     ORBITALS ARE SET TO TWO, AND ALL OTHERS ARE SET TO ZERO.
!     *
!     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
!     Ra/Rb(LM3,LM3)  DESITY MATRIX(I).
!     W(LM6,LM6)      TWO-ELECTRON INTEGRALS(I).
!     F2a/F2b(LM3,LM3)      CIS FOCK-LIKE MATRIX(O).

      USE LIMIT, ONLY: LM1, LMI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      COMMON /CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
      COMMON /FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
      COMMON /INDEXW/ NW(LM1)
      REAL*8 Qa(LM3,LM3,*),Qb(LM3,LM3,*),Ja(LM3,LM3,*),Jb(LM3,LM3,*),
     $       Ka(LM3,LM3,*),Kb(LM3,LM3,*),W(LM6,LM6)
      INTEGER NQ,LM3,LM6,NDen

      CALL VecInit(Ja,LM3*LM3*NQ,0d0)
      CALL VecInit(Ka,LM3*LM3*NQ,0d0)
      IF(NDen.EQ.2) THEN
        CALL VecInit(Jb,LM3*LM3*NQ,0d0)
        CALL VecInit(Kb,LM3*LM3*NQ,0d0)
      ENDIF
!     IP1(I).GE.IP2(I)
      DO II=1,NUMAT
      NA = NW(II) - 1
      IA = NFIRST(II)
      IB = NLAST(II)
      DO IMU = IA,IB
      DO INU = IA,IMU  
      NA = NA + 1
         DO JJ=1,NUMAT
         IF((II.EQ.ICNTR .OR. JJ.EQ.ICNTR) .AND. II.NE.JJ) THEN
         !IF((II.EQ.ICNTR .OR. JJ.EQ.ICNTR)) THEN
         !write(6,*) "atoms:",II,JJ
         NB = NW(JJ) - 1
         LA = NFIRST(JJ)
         LB = NLAST(JJ)
         DO IRO = LA,LB
         DO ISI = LA,IRO
         NB = NB + 1
         IF(NB.GE.NA)  THEN   
            Xij = W(NA,NB)
            Yij = -W(NA,NB)
!           COULOMB PART.
!           (MU,NU | RO,SI) + (MU,NU | SI,RO)  
!           (NU,MU | RO,SI) + (NU,MU | SI,RO)
            DO M=1,NQ
            Rho1 = Qa(IRO,ISI,M) !+ Qb(IRO,ISI,M)
            IF(IRO.NE.ISI)  Rho1 = Rho1 + Qa(ISI,IRO,M) 
!     $                                  + Qb(ISI,IRO,M)
            Rho2 = Qa(IMU,INU,M) !+ Qb(IMU,INU,M)
            IF(IMU.NE.INU)  Rho2 = Rho2 + Qa(INU,IMU,M) 
!     $                                  + Qb(INU,IMU,M)
            
!            write(6,*) IMU,INU,"--",IRO,ISI
!            if(IMU.NE.INU) write(6,*) INU,IMU,"--",IRO,ISI
!            if(NA.NE.NB) THEN
!              write(6,*) IRO,ISI,"--",IMU,INU
!              if(IRO.NE.ISI) write(6,*) ISI,IRO,"--",IMU,INU
!            endif
            Ja(IMU,INU,M) = Ja(IMU,INU,M) + Rho1*Xij
            IF(IMU.NE.INU) Ja(INU,IMU,M) = Ja(INU,IMU,M) + Rho1*Xij
            IF(NA.NE.NB) THEN
              Ja(IRO,ISI,M) = Ja(IRO,ISI,M) + Rho2*Xij
              IF(IRO.NE.ISI)  
     $           Ja(ISI,IRO,M) = Ja(ISI,IRO,M) + Rho2*Xij
            ENDIF

!           Yij = 0d0
!           EXCHANGE PART.
!           (MU,NU | RO,SI)
            Ka(IMU,IRO,M) = Ka(IMU,IRO,M) + Qa(INU,ISI,M)*Yij
!            write(6,*) M,":",IMU,IRO,"--",INU,ISI
            IF(IMU.NE.IRO .OR. INU.NE.ISI)   
     $        Ka(IRO,IMU,M) = Ka(IRO,IMU,M) + Qa(ISI,INU,M)*Yij
!            IF(IMU.NE.IRO .OR. INU.NE.ISI) 
!     $        write(6,*) M,":",IRO,IMU,"--",ISI,INU 
            IF(IRO.NE.ISI) THEN
!             (MU,NU | SI,RO)
              Ka(IMU,ISI,M) = Ka(IMU,ISI,M) + Qa(INU,IRO,M)*Yij
!              write(6,*) M,":",IMU,ISI,"--",INU,IRO
              IF(ISI.NE.IMU .OR. INU.NE.IRO) 
     $          Ka(ISI,IMU,M) = Ka(ISI,IMU,M) + Qa(IRO,INU,M)*Yij
!              IF(IMU.NE.ISI .OR. INU.NE.IRO) 
!     $          write(6,*) M,":",ISI,IMU,"--",IRO,INU 
            ENDIF
            IF(IMU.NE.INU ) THEN
!             (NU,MU | RO,SI)
              IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
              Ka(INU,IRO,M) = Ka(INU,IRO,M) + Qa(IMU,ISI,M)*Yij
!              write(6,*) M,":",INU,IRO,"--",IMU,ISI
              IF(INU.NE.IRO .OR. IMU.NE.ISI) 
     $          Ka(IRO,INU,M) = Ka(IRO,INU,M) + Qa(ISI,IMU,M)*Yij
!              IF(INU.NE.IRO .OR. IMU.NE.ISI) 
!     $          write(6,*) M,":",IRO,INU,"--",ISI,IMU
              ENDIF
              IF(IRO.NE.ISI) THEN
!               (NU,MU | SI,RO)
                Ka(INU,ISI,M) = Ka(INU,ISI,M) + Qa(IMU,IRO,M)*Yij
!                write(6,*) M,":",INU,ISI,"--",IMU,IRO
                IF(INU.NE.ISI .OR. IMU.NE.IRO) 
     $            Ka(ISI,INU,M) = Ka(ISI,INU,M) + Qa(IRO,IMU,M)*Yij
!                IF(INU.NE.ISI .OR. IMU.NE.IRO) 
!     $            write(6,*) M,":",ISI,INU,"--",IRO,IMU 
              ENDIF
            ENDIF

            IF(NDen.EQ.2) THEN
!             COULOMB PART.
!             (MU,NU | RO,SI) + (MU,NU | SI,RO)  
!             (NU,MU | RO,SI) + (NU,MU | SI,RO)
              Rho1 = Qb(IRO,ISI,M) !+ Qb(IRO,ISI,M)
              IF(IRO.NE.ISI)  Rho1 = Rho1 + Qb(ISI,IRO,M)
!     $                                    + Qb(ISI,IRO,M)
              Rho2 = Qb(IMU,INU,M) !+ Qb(IMU,INU,M)
              IF(IMU.NE.INU)  Rho2 = Rho2 + Qb(INU,IMU,M)
!     $                                    + Qb(INU,IMU,M)
              Jb(IMU,INU,M) = Jb(IMU,INU,M) + Rho1*Xij
              IF(IMU.NE.INU)  
     $          Jb(INU,IMU,M) = Jb(INU,IMU,M) + Rho1*Xij
              IF(NA.NE.NB) THEN
                Jb(IRO,ISI,M) = Jb(IRO,ISI,M) + Rho2*Xij
                IF(IRO.NE.ISI)  
     $            Jb(ISI,IRO,M) = Jb(ISI,IRO,M) + Rho2*Xij
              ENDIF

!             EXCHANGE PART.
!             (MU,NU | RO,SI)
              Kb(IMU,IRO,M) = Kb(IMU,IRO,M) + Qb(INU,ISI,M)*Yij
              IF(IMU.NE.IRO .OR. INU.NE.ISI) 
     $          Kb(IRO,IMU,M) = Kb(IRO,IMU,M) + Qb(ISI,INU,M)*Yij
              IF(IRO.NE.ISI) THEN
!               (MU,NU | SI,RO)
                Kb(IMU,ISI,M) = Kb(IMU,ISI,M) + Qb(INU,IRO,M)*Yij
                IF(ISI.NE.IMU .OR. INU.NE.IRO) 
     $            Kb(ISI,IMU,M) = Kb(ISI,IMU,M) + Qb(IRO,INU,M)*Yij
              ENDIF
              IF(IMU.NE.INU ) THEN
!               (NU,MU | RO,SI)
                IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
                Kb(INU,IRO,M) = Kb(INU,IRO,M) + Qb(IMU,ISI,M)*Yij
                IF(INU.NE.IRO .OR. IMU.NE.ISI) 
     $            Kb(IRO,INU,M) = Kb(IRO,INU,M) + Qb(ISI,IMU,M)*Yij
                ENDIF
                IF(IRO.NE.ISI) THEN
!                 (NU,MU | SI,RO)
                  Kb(INU,ISI,M) = Kb(INU,ISI,M) + Qb(IMU,IRO,M)*Yij
                  IF(INU.NE.ISI .OR. IMU.NE.IRO) 
     $              Kb(ISI,INU,M) = Kb(ISI,INU,M) + Qb(IRO,IMU,M)*Yij
                ENDIF
              ENDIF
            ENDIF ! Beta Spin
            ENDDO ! NQ
          ENDIF ! NA,NB

         ENDDO !LB
         ENDDO !LA
         ENDIF !ICNTR
         ENDDO !JJ
      ENDDO !IB
      ENDDO !IA
      ENDDO !II

      RETURN
      END

      SUBROUTINE FOCK132b(Qa,Qb,Ja,Jb,Ka,Kb,NQ,W,LM3,LM6,NDen,ICNTR)
!     *
!     BUILD CIS FOCK-LIKE MATRIX.
!     *
!     REFERENCE OCCUPATION NUMBERS OF DOUBLY OCCUPIED (CLOSED-SHELL)
!     ORBITALS ARE SET TO TWO, AND ALL OTHERS ARE SET TO ZERO.
!     *
!     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
!     Ra/Rb(LM3,LM3)  DESITY MATRIX(I).
!     W(LM6,LM6)      TWO-ELECTRON INTEGRALS(I).
!     F2a/F2b(LM3,LM3)      CIS FOCK-LIKE MATRIX(O).

      USE LIMIT, ONLY: LM1, LMI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON  /ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
      COMMON  /CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
      COMMON  /FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
      COMMON  /INDEXW/ NW(LM1)
      REAL*8  Qa(NQ,*),Qb(NQ,*),Ja(NQ,*),Jb(NQ,*),
     $        Ka(NQ,*),Kb(NQ,*),W(LM6,LM6)
      INTEGER NQ,LM3,LM6,NDen

      CALL VecInit(Ja,LM3*LM3*NQ,0d0)
      CALL VecInit(Ka,LM3*LM3*NQ,0d0)
      IF(NDen.EQ.2) THEN
        CALL VecInit(Jb,LM3*LM3*NQ,0d0)
        CALL VecInit(Kb,LM3*LM3*NQ,0d0)
      ENDIF

      CALL MATTRANS2(Qa,LM3*LM3,NQ)
      IF(NDEN.EQ.2) CALL MATTRANS2(Qb,LM3*LM3,NQ)

!     IP1(I).GE.IP2(I)
      DO II=1,NUMAT
      NA = NW(II) - 1
      IA = NFIRST(II)
      IB = NLAST(II)
      DO IMU = IA,IB
      DO INU = IA,IMU  
      NA = NA + 1
         DO JJ=1,NUMAT
         IF((II.EQ.ICNTR .OR. JJ.EQ.ICNTR) .AND. II.NE.JJ) THEN
         !IF((II.EQ.ICNTR .OR. JJ.EQ.ICNTR)) THEN
         !write(6,*) "atoms:",II,JJ
         NB = NW(JJ) - 1
         LA = NFIRST(JJ)
         LB = NLAST(JJ)
         DO IRO = LA,LB
         DO ISI = LA,IRO
         NB = NB + 1
         IF(NB.GE.NA)  THEN   
            IMN = IMU+(INU-1)*LM3
            INM = INU+(IMU-1)*LM3
            IRS = IRO+(ISI-1)*LM3
            ISR = ISI+(IRO-1)*LM3
               
            IMR = IMU+(IRO-1)*LM3 
            IMS = IMU+(ISI-1)*LM3
            INR = INU+(IRO-1)*LM3
            INS = INU+(ISI-1)*LM3 
            IRM = IRO+(IMU-1)*LM3
            ISM = ISI+(IMU-1)*LM3
            IRN = IRO+(INU-1)*LM3
            ISN = ISI+(INU-1)*LM3

            Xij = W(NA,NB)
            Yij = -W(NA,NB)
!           COULOMB PART.
!           (MU,NU | RO,SI) + (MU,NU | SI,RO)  
!           (NU,MU | RO,SI) + (NU,MU | SI,RO)
            DO M=1,NQ
            Rho1 = Qa(M,IRS) 
            IF(IRO.NE.ISI)  Rho1 = Rho1 + Qa(M,ISR) 
            Rho1x = Rho1*Xij
            Ja(M,IMN) = Ja(M,IMN) + Rho1x
            IF(IMU.NE.INU) Ja(M,INM) = Ja(M,INM) + Rho1x
            IF(NA.NE.NB) THEN
              Rho2 = Qa(M,IMN) 
              IF(IMU.NE.INU)  Rho2 = Rho2 + Qa(M,INM) 
              Rho2x = Rho2*Xij
              Ja(M,IRS) = Ja(M,IRS) + Rho2x
              IF(IRO.NE.ISI)  Ja(M,ISR) = Ja(M,ISR) + Rho2x
            ENDIF

!           EXCHANGE PART.
!           (MU,NU | RO,SI)
            Ka(M,IMR) = Ka(M,IMR) + Qa(M,INS)*Yij
            IF(IMU.NE.IRO .OR. INU.NE.ISI)   
     $        Ka(M,IRM) = Ka(M,IRM) + Qa(M,ISN)*Yij
            IF(IRO.NE.ISI) THEN
!             (MU,NU | SI,RO)
              Ka(M,IMS) = Ka(M,IMS) + Qa(M,INR)*Yij
              IF(ISI.NE.IMU .OR. INU.NE.IRO) 
     $          Ka(M,ISM) = Ka(M,ISM) + Qa(M,IRN)*Yij
            ENDIF
            IF(IMU.NE.INU ) THEN
!             (NU,MU | RO,SI)
              IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
              Ka(M,INR) = Ka(M,INR) + Qa(M,IMS)*Yij
              IF(INU.NE.IRO .OR. IMU.NE.ISI) 
     $          Ka(M,IRN) = Ka(M,IRN) + Qa(M,ISM)*Yij
              ENDIF
              IF(IRO.NE.ISI) THEN
!               (NU,MU | SI,RO)
                Ka(M,INS) = Ka(M,INS) + Qa(M,IMR)*Yij
                IF(INU.NE.ISI .OR. IMU.NE.IRO) 
     $            Ka(M,ISN) = Ka(M,ISN) + Qa(M,IRM)*Yij
              ENDIF
            ENDIF

            IF(NDen.EQ.2) THEN
!             COULOMB PART.
!             (MU,NU | RO,SI) + (MU,NU | SI,RO)  
!             (NU,MU | RO,SI) + (NU,MU | SI,RO)
              Rho1 = Qb(M,IRS) 
              IF(IRO.NE.ISI)  Rho1 = Rho1 + Qb(M,ISR)
              Jb(M,IMN) = Jb(M,IMN) + Rho1*Xij
              IF(IMU.NE.INU) Jb(M,INM) = Jb(M,INM) + Rho1*Xij
              IF(NA.NE.NB) THEN
                Rho2 = Qb(M,IMN) 
                IF(IMU.NE.INU)  Rho2 = Rho2 + Qb(M,INM)
                Jb(M,IRS) = Jb(M,IRS) + Rho2*Xij
                IF(IRO.NE.ISI) Jb(M,ISR) = Jb(M,ISR) + Rho2*Xij
              ENDIF

!             EXCHANGE PART.
!             (MU,NU | RO,SI)
              Kb(M,IMR) = Kb(M,IMR) + Qb(M,INS)*Yij
              IF(IMU.NE.IRO .OR. INU.NE.ISI) 
     $          Kb(M,IRM) = Kb(M,IRM) + Qb(M,ISN)*Yij
              IF(IRO.NE.ISI) THEN
!               (MU,NU | SI,RO)
                Kb(M,IMS) = Kb(M,IMS) + Qb(M,INR)*Yij
                IF(ISI.NE.IMU .OR. INU.NE.IRO) 
     $            Kb(M,ISM) = Kb(M,ISM) + Qb(M,IRN)*Yij
              ENDIF
              IF(IMU.NE.INU ) THEN
!               (NU,MU | RO,SI)
                IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
                Kb(M,INR) = Kb(M,INR) + Qb(M,IMS)*Yij
                IF(INU.NE.IRO .OR. IMU.NE.ISI) 
     $            Kb(M,IRN) = Kb(M,IRN) + Qb(M,ISM)*Yij
                ENDIF
                IF(IRO.NE.ISI) THEN
!                 (NU,MU | SI,RO)
                  Kb(M,INS) = Kb(M,INS) + Qb(M,IMR)*Yij
                  IF(INU.NE.ISI .OR. IMU.NE.IRO) 
     $              Kb(M,ISN) = Kb(M,ISN) + Qb(M,IRM)*Yij
                ENDIF
              ENDIF
            ENDIF ! Beta Spin
            ENDDO ! NQ
          ENDIF ! NA,NB

         ENDDO !LB
         ENDDO !LA
         ENDIF !ICNTR
         ENDDO !JJ
      ENDDO !IB
      ENDDO !IA
      ENDDO !II

      CALL MATTRANS2(Qa,NQ,LM3*LM3)
      IF(NDEN.EQ.2) CALL MATTRANS2(Qb,NQ,LM3*LM3)

      CALL MATTRANS2(Ja,NQ,LM3*LM3)
      CALL MATTRANS2(Ka,NQ,LM3*LM3)
      IF(NDEN.EQ.2) THEN
        CALL MATTRANS2(Jb,NQ,LM3*LM3)
        CALL MATTRANS2(Kb,NQ,LM3*LM3)
      ENDIF


      RETURN
      END

