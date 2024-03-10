      SUBROUTINE USE131 (A,LDA,P,FDRV,ICNTR,ESCAL,NDEN,I1n)
      USE LIMIT, ONLY: LM1, LMX, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ENERGT/ EE,ENUCLR,EAT,ATHEAT
     ./LIMITS/ LM2,LM3,LM4,LM6,LM7,LM8,LM9
     ./LMGRAD/ LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8
     ./LMSCF / LS1,LS2,LS3,LS4,LS5,LS6,LS7,LS8,LS9
     ./LMUHF / LU1,LU2,LU3,LU4,LU5,LU6,LU7,LU8
      DIMENSION A(LDA),P(LM3*LM3,*),FDRV(LM3*LM3,*)
      REAL*8, DIMENSION(:), ALLOCATABLE :: TMP
      ALLOCATE(TMP(LM3*LM3))

C *** SIMPLIFIED CODE FOR GRADIENT CALCULATIONS WITH DISPLACED ATOM.
C *** INTEGRAL CALCULATION FOR SCF.
      CALL HCORE2 (A(LS7),A(LS9),A(LG1),A(LG2),A(LG3),A(LG4),
     1             LM2,LM4,LM6,LM9,ICNTR)

C *** CORE REPULSION ENERGY.
      ESCAL = ENUCLR

C *** TWO-ELECTRON CONTRIBUTIONS TO FOCK MATRIX FOR GIVEN DENSITY.
      CALL FOCK131(FDRV(1,1),FDRV(1,NDEN),P(1,1),P(1,NDEN),
     $             FDRV(1,NDEN+1),FDRV(1,2*NDEN),P(1,NDEN+1),
     $             P(1,2*NDEN),A(LS9),LM3,LM6,NDEN,I1n,ICNTR)
      !CALL MATPRNT(FDRV,LM3,LM3,6)

      !call vecinit(tmp,LM3*LM3,0d0)
      !CALL FOCKX (tmp,A(LS8),A(LU8),FDRV,A(LS9),LM4,LM6)
      !call square(tmp,FDRV,LM3,LM3,LM4)
      !CALL MATPRNT(FDRV,LM3,LM3,6)

      CALL SQUARE(A(LS7),TMP,LM3,LM3,LM4)
      !CALL MATPRNT(TMP,LM3,LM3,6)
      !CALL MATPRNT(FDRV,LM3,LM3,6)
      !CALL MATPRNT(FDRV(1,NDEN+1),LM3,LM3,6)
      FDRV(1:LM3*LM3,1) = FDRV(1:LM3*LM3,1) + TMP(1:LM3*LM3)
      IF(NDEN.EQ.2) THEN 
        CALL SQUARE(A(LU7),TMP,LM3,LM3,LM4)
        FDRV(1:LM3*LM3,2) = FDRV(1:LM3*LM3,2) + TMP(1:LM3*LM3)
      ENDIF

      DEALLOCATE(TMP)

      RETURN
      END


      SUBROUTINE USE131a (A,LDA,P,FDRV,ICNTR,ESCAL,NDEN,I1n)
      USE LIMIT, ONLY: LM1, LMX, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./ENERGT/ EE,ENUCLR,EAT,ATHEAT
     ./LIMITS/ LM2,LM3,LM4,LM6,LM7,LM8,LM9
     ./LMGRAD/ LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8
     ./LMSCF / LS1,LS2,LS3,LS4,LS5,LS6,LS7,LS8,LS9
     ./LMUHF / LU1,LU2,LU3,LU4,LU5,LU6,LU7,LU8
      DIMENSION A(LDA),P(LM3*LM3,*),FDRV(LM3*LM3,*)

C *** SIMPLIFIED CODE FOR GRADIENT CALCULATIONS WITH DISPLACED ATOM.
C *** INTEGRAL CALCULATION FOR SCF.
      CALL HCORE2 (A(LS7),A(LS9),A(LG1),A(LG2),A(LG3),A(LG4),
     1             LM2,LM4,LM6,LM9,ICNTR)

C *** CORE REPULSION ENERGY.
      ESCAL = ENUCLR

C *** TWO-ELECTRON CONTRIBUTIONS TO FOCK MATRIX FOR GIVEN DENSITY.
      CALL FOCK131(FDRV(1,1),FDRV(1,NDEN),P(1,1),P(1,NDEN),
     $             FDRV(1,NDEN+1),FDRV(1,2*NDEN),P(1,NDEN+1),
     $             P(1,2*NDEN),A(LS9),LM3,LM6,NDEN,I1n,ICNTR)
      !CALL MATPRNT(FDRV,LM3,LM3,6)

      !call vecinit(tmp,LM3*LM3,0d0)
      !CALL FOCKX (tmp,A(LS8),A(LU8),FDRV,A(LS9),LM4,LM6)
      !call square(tmp,FDRV,LM3,LM3,LM4)
      !CALL MATPRNT(FDRV,LM3,LM3,6)

      CALL SQUARE(A(LS7),FDRV(1,2*NDEN+1),LM3,LM3,LM4)
      !CALL MATPRNT(TMP,LM3,LM3,6)
      !CALL MATPRNT(FDRV,LM3,LM3,6)
      !CALL MATPRNT(FDRV(1,NDEN+1),LM3,LM3,6)
      IF(NDEN.EQ.2) THEN 
        CALL SQUARE(A(LU7),FDRV(1,3*NDEN),LM3,LM3,LM4)
      ENDIF

      RETURN
      END



      SUBROUTINE FOCK131(F1a,F1b,Qa,Qb,F2a,F2b,Ra,Rb,W,LM3,LM6,
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
      REAL*8 Qa(LM3,LM3),Qb(LM3,LM3),F1a(LM3,LM3),F1b(LM3,LM3),
     $       Ra(LM3,LM3),Rb(LM3,LM3),F2a(LM3,LM3),F2b(LM3,LM3),
     $       W(LM6,LM6)
      INTEGER NL,LM3,LM6,NDen

      CALL VecInit(F1a,LM3*LM3,0d0)
      CALL VecInit(F2a,LM3*LM3,0d0)
      IF(NDen.EQ.2) THEN
        CALL VecInit(F1b,LM3*LM3,0d0)
        CALL VecInit(F2b,LM3*LM3,0d0)
      ENDIF
      !call matprnt(Qa,LM3,LM3,6)
      !call matprnt(Ra,LM3,LM3,6)
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
            IF(I1n.EQ.1) THEN
               Rho1 = Ra(IRO,ISI) + Rb(IRO,ISI)
               IF(IRO.NE.ISI)  Rho1 = Rho1 + Ra(ISI,IRO) + Rb(ISI,IRO)
               Rho2 = Ra(IMU,INU) + Rb(IMU,INU)
               IF(IMU.NE.INU)  Rho2 = Rho2 + Ra(INU,IMU) + Rb(INU,IMU)
             
               F2a(IMU,INU) = F2a(IMU,INU) + Rho1*Xij
               IF(IMU.NE.INU) F2a(INU,IMU) = F2a(INU,IMU) + Rho1*Xij
               IF(NA.NE.NB) THEN
                 F2a(IRO,ISI) = F2a(IRO,ISI) + Rho2*Xij
                 IF(IRO.NE.ISI)  F2a(ISI,IRO) = F2a(ISI,IRO) + Rho2*Xij
               ENDIF
            ENDIF

!           EXCHANGE PART.
!           (MU,NU | RO,SI)
            F2a(IMU,IRO) = F2a(IMU,IRO) + Ra(INU,ISI)*Yij
!            write(6,*) M,":",IMU,IRO,"--",INU,ISI
            IF(IMU.NE.IRO .OR. INU.NE.ISI)   
     $        F2a(IRO,IMU) = F2a(IRO,IMU) + Ra(ISI,INU)*Yij
!            IF(IMU.NE.IRO .OR. INU.NE.ISI) 
!     $        write(6,*) M,":",IRO,IMU,"--",ISI,INU 
            IF(IRO.NE.ISI) THEN
!             (MU,NU | SI,RO)
              F2a(IMU,ISI) = F2a(IMU,ISI) + Ra(INU,IRO)*Yij
!              write(6,*) M,":",IMU,ISI,"--",INU,IRO
              IF(ISI.NE.IMU .OR. INU.NE.IRO) 
     $          F2a(ISI,IMU) = F2a(ISI,IMU) + Ra(IRO,INU)*Yij
!              IF(IMU.NE.ISI .OR. INU.NE.IRO) 
!     $          write(6,*) M,":",ISI,IMU,"--",IRO,INU 
            ENDIF
            IF(IMU.NE.INU ) THEN
!             (NU,MU | RO,SI)
              IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
              F2a(INU,IRO) = F2a(INU,IRO) + Ra(IMU,ISI)*Yij
!              write(6,*) M,":",INU,IRO,"--",IMU,ISI
              IF(INU.NE.IRO .OR. IMU.NE.ISI) 
     $          F2a(IRO,INU) = F2a(IRO,INU) + Ra(ISI,IMU)*Yij
!              IF(INU.NE.IRO .OR. IMU.NE.ISI) 
!     $          write(6,*) M,":",IRO,INU,"--",ISI,IMU
              ENDIF
              IF(IRO.NE.ISI) THEN
!               (NU,MU | SI,RO)
                F2a(INU,ISI) = F2a(INU,ISI) + Ra(IMU,IRO)*Yij
!                write(6,*) M,":",INU,ISI,"--",IMU,IRO
                IF(INU.NE.ISI .OR. IMU.NE.IRO) 
     $            F2a(ISI,INU) = F2a(ISI,INU) + Ra(IRO,IMU)*Yij
!                IF(INU.NE.ISI .OR. IMU.NE.IRO) 
!     $            write(6,*) M,":",ISI,INU,"--",IRO,IMU 
              ENDIF
            ENDIF

            IF(NDen.EQ.2) THEN
!             COULOMB PART.
!             (MU,NU | RO,SI) + (MU,NU | SI,RO)  
!             (NU,MU | RO,SI) + (NU,MU | SI,RO)
              F2b(IMU,INU) = F2b(IMU,INU) + Rho1*Xij
              IF(IMU.NE.INU)  
     $          F2b(INU,IMU) = F2b(INU,IMU) + Rho1*Xij
              IF(NA.NE.NB) THEN
                F2b(IRO,ISI) = F2b(IRO,ISI) + Rho2*Xij
                IF(IRO.NE.ISI)  
     $            F2b(ISI,IRO) = F2b(ISI,IRO) + Rho2*Xij
              ENDIF

!             EXCHANGE PART.
!             (MU,NU | RO,SI)
              F2b(IMU,IRO) = F2b(IMU,IRO) + Rb(INU,ISI)*Yij
              IF(IMU.NE.IRO .OR. INU.NE.ISI) 
     $          F2b(IRO,IMU) = F2b(IRO,IMU) + Rb(ISI,INU)*Yij
              IF(IRO.NE.ISI) THEN
!               (MU,NU | SI,RO)
                F2b(IMU,ISI) = F2b(IMU,ISI) + Rb(INU,IRO)*Yij
                IF(ISI.NE.IMU .OR. INU.NE.IRO) 
     $            F2b(ISI,IMU) = F2b(ISI,IMU) + Rb(IRO,INU)*Yij
              ENDIF
              IF(IMU.NE.INU ) THEN
!               (NU,MU | RO,SI)
                IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
                F2b(INU,IRO) = F2b(INU,IRO) + Rb(IMU,ISI)*Yij
                IF(INU.NE.IRO .OR. IMU.NE.ISI) 
     $            F2b(IRO,INU) = F2b(IRO,INU) + Rb(ISI,IMU)*Yij
                ENDIF
                IF(IRO.NE.ISI) THEN
!                 (NU,MU | SI,RO)
                  F2b(INU,ISI) = F2b(INU,ISI) + Rb(IMU,IRO)*Yij
                  IF(INU.NE.ISI .OR. IMU.NE.IRO) 
     $              F2b(ISI,INU) = F2b(ISI,INU) + Rb(IRO,IMU)*Yij
                ENDIF
              ENDIF

            ENDIF ! Beta Spin

            IF(I1n.EQ.1) THEN
               Rho1 = Qa(IRO,ISI) + Qb(IRO,ISI)
               IF(IRO.NE.ISI)  Rho1 = Rho1 + Qa(ISI,IRO) + Qb(ISI,IRO)
               Rho2 = Qa(IMU,INU) + Qb(IMU,INU)
               IF(IMU.NE.INU)  Rho2 = Rho2 + Qa(INU,IMU) + Qb(INU,IMU)
             
!               write(6,*) IMU,INU,"--",IRO,ISI
!               if(IMU.NE.INU) write(6,*) INU,IMU,"--",IRO,ISI
!               if(NA.NE.NB) THEN
!                 write(6,*) IRO,ISI,"--",IMU,INU
!                 if(IRO.NE.ISI) write(6,*) ISI,IRO,"--",IMU,INU
!               endif
               F1a(IMU,INU) = F1a(IMU,INU) + Rho1*Xij
               IF(IMU.NE.INU) F1a(INU,IMU) = F1a(INU,IMU) + Rho1*Xij
               IF(NA.NE.NB) THEN
                 F1a(IRO,ISI) = F1a(IRO,ISI) + Rho2*Xij
                 IF(IRO.NE.ISI)  F1a(ISI,IRO) = F1a(ISI,IRO) + Rho2*Xij
               ENDIF
            ENDIF

!            Yij = 0d0
!           EXCHANGE PART.
!           (MU,NU | RO,SI)
            F1a(IMU,IRO) = F1a(IMU,IRO) + Qa(INU,ISI)*Yij
!            write(6,*) M,":",IMU,IRO,"--",INU,ISI
            IF(IMU.NE.IRO .OR. INU.NE.ISI)   
     $        F1a(IRO,IMU) = F1a(IRO,IMU) + Qa(ISI,INU)*Yij
!            IF(IMU.NE.IRO .OR. INU.NE.ISI) 
!     $        write(6,*) M,":",IRO,IMU,"--",ISI,INU 
            IF(IRO.NE.ISI) THEN
!             (MU,NU | SI,RO)
              F1a(IMU,ISI) = F1a(IMU,ISI) + Qa(INU,IRO)*Yij
!              write(6,*) M,":",IMU,ISI,"--",INU,IRO
              IF(ISI.NE.IMU .OR. INU.NE.IRO) 
     $          F1a(ISI,IMU) = F1a(ISI,IMU) + Qa(IRO,INU)*Yij
!              IF(IMU.NE.ISI .OR. INU.NE.IRO) 
!     $          write(6,*) M,":",ISI,IMU,"--",IRO,INU 
            ENDIF
            IF(IMU.NE.INU ) THEN
!             (NU,MU | RO,SI)
              IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
              F1a(INU,IRO) = F1a(INU,IRO) + Qa(IMU,ISI)*Yij
!              write(6,*) M,":",INU,IRO,"--",IMU,ISI
              IF(INU.NE.IRO .OR. IMU.NE.ISI) 
     $          F1a(IRO,INU) = F1a(IRO,INU) + Qa(ISI,IMU)*Yij
!              IF(INU.NE.IRO .OR. IMU.NE.ISI) 
!     $          write(6,*) M,":",IRO,INU,"--",ISI,IMU
              ENDIF
              IF(IRO.NE.ISI) THEN
!               (NU,MU | SI,RO)
                F1a(INU,ISI) = F1a(INU,ISI) + Qa(IMU,IRO)*Yij
!                write(6,*) M,":",INU,ISI,"--",IMU,IRO
                IF(INU.NE.ISI .OR. IMU.NE.IRO) 
     $            F1a(ISI,INU) = F1a(ISI,INU) + Qa(IRO,IMU)*Yij
!                IF(INU.NE.ISI .OR. IMU.NE.IRO) 
!     $            write(6,*) M,":",ISI,INU,"--",IRO,IMU 
              ENDIF
            ENDIF

            IF(NDen.EQ.2) THEN
!             COULOMB PART.
!             (MU,NU | RO,SI) + (MU,NU | SI,RO)  
!             (NU,MU | RO,SI) + (NU,MU | SI,RO)
              F1b(IMU,INU) = F1b(IMU,INU) + Rho1*Xij
              IF(IMU.NE.INU)  
     $          F1b(INU,IMU) = F1b(INU,IMU) + Rho1*Xij
              IF(NA.NE.NB) THEN
                F1b(IRO,ISI) = F1b(IRO,ISI) + Rho2*Xij
                IF(IRO.NE.ISI)  
     $            F1b(ISI,IRO) = F1b(ISI,IRO) + Rho2*Xij
              ENDIF

!             EXCHANGE PART.
!             (MU,NU | RO,SI)
              F1b(IMU,IRO) = F1b(IMU,IRO) + Qb(INU,ISI)*Yij
              IF(IMU.NE.IRO .OR. INU.NE.ISI) 
     $          F1b(IRO,IMU) = F1b(IRO,IMU) + Qb(ISI,INU)*Yij
              IF(IRO.NE.ISI) THEN
!               (MU,NU | SI,RO)
                F1b(IMU,ISI) = F1b(IMU,ISI) + Qb(INU,IRO)*Yij
                IF(ISI.NE.IMU .OR. INU.NE.IRO) 
     $            F1b(ISI,IMU) = F1b(ISI,IMU) + Qb(IRO,INU)*Yij
              ENDIF
              IF(IMU.NE.INU ) THEN
!               (NU,MU | RO,SI)
                IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
                F1b(INU,IRO) = F1b(INU,IRO) + Qb(IMU,ISI)*Yij
                IF(INU.NE.IRO .OR. IMU.NE.ISI) 
     $            F1b(IRO,INU) = F1b(IRO,INU) + Qb(ISI,IMU)*Yij
                ENDIF
                IF(IRO.NE.ISI) THEN
!                 (NU,MU | SI,RO)
                  F1b(INU,ISI) = F1b(INU,ISI) + Qb(IMU,IRO)*Yij
                  IF(INU.NE.ISI .OR. IMU.NE.IRO) 
     $              F1b(ISI,INU) = F1b(ISI,INU) + Qb(IRO,IMU)*Yij
                ENDIF
              ENDIF

            ENDIF ! Beta Spin
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

