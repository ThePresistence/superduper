!
!     CREATED BY JIE LIU ON 3/2016 
!
SUBROUTINE MAKEFAO(Fa,Fb,Pa,Pb,W,NL,LM3,LM6,NDEN,I1n,Cj,Ck)
!*
!BUILD CIS FOCK-LIKE MATRIX.
!*
!REFERENCE OCCUPATION NUMBERS OF DOUBLY OCCUPIED (CLOSED-SHELL)
!ORBITALS ARE SET TO TWO, AND ALL OTHERS ARE SET TO ZERO.
!*
!NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
!Pa/Pb(LM3,LM3)  DESITY MATRIX(I).
!W(LM6,LM6)      TWO-ELECTRON INTEGRALS(I).
!Fa/Fb(LM3,LM3)      CIS FOCK-LIKE MATRIX(O).
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  REAL*8  Pa(LM3,LM3,*),Pb(LM3,LM3,*),Fa(LM3,LM3,*),Fb(LM3,LM3,*),W(LM6,LM6)
  INTEGER NL,LM3,LM6,NDEN
  REAL*8  Cj,Ck

  IF(NL.GT.4) THEN
    CALL MAKEFAOx2(Fa,Fb,Pa,Pb,W,NL,LM3,LM6,NDEN,I1n,Cj,Ck)
  ELSE
    CALL MAKEFAOx1(Fa,Fb,Pa,Pb,W,NL,LM3,LM6,NDEN,I1n,Cj,Ck)
  ENDIF

END SUBROUTINE MAKEFAO

SUBROUTINE MAKEFAOx1(Fa,Fb,Pa,Pb,W,NL,LM3,LM6,NDEN,I1n,Cj,Ck)
!*
!BUILD CIS FOCK-LIKE MATRIX.
!*
!REFERENCE OCCUPATION NUMBERS OF DOUBLY OCCUPIED (CLOSED-SHELL)
!ORBITALS ARE SET TO TWO, AND ALL OTHERS ARE SET TO ZERO.
!*
!NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
!Pa/Pb(LM3,LM3)  DESITY MATRIX(I).
!W(LM6,LM6)      TWO-ELECTRON INTEGRALS(I).
!Fa/Fb(LM3,LM3)      CIS FOCK-LIKE MATRIX(O).
  USE LIMIT, ONLY: LMI
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  COMMON /CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
  COMMON /FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
  REAL*8 Pa(LM3,LM3,*),Pb(LM3,LM3,*),Fa(LM3,LM3,*),Fb(LM3,LM3,*),W(LM6,LM6)
  INTEGER NL,LM3,LM6,NDEN
  REAL*8 Cj,Ck

  CALL VecInit(Fa,LM3*LM3*NL,0d0)
  IF(NDEN.EQ.2) CALL VecInit(Fb,LM3*LM3*NL,0d0)
  !IP1(I).GE.IP2(I)
  DO I=1,LM6
     IMU = IP1(I)
     INU = IP2(I)
     DO J=I,LM6
        IRO = IP1(J)
        ISI = IP2(J)
        Xij = W(I,J)*Cj
        Yij = -W(I,J)*Ck
        DO M=1,NL 
           !COULOMB PART.
           !(MU,NU | RO,SI) + (MU,NU | SI,RO)  
           !(NU,MU | RO,SI) + (NU,MU | SI,RO)
           IF(I1n.EQ.1) THEN
              Rho1 = Pa(IRO,ISI,M) + Pb(IRO,ISI,M)
              IF(IRO.NE.ISI)  Rho1 = Rho1 + Pa(ISI,IRO,M) + Pb(ISI,IRO,M)
              Rho2 = Pa(IMU,INU,M) + Pb(IMU,INU,M)
              IF(IMU.NE.INU)  Rho2 = Rho2 + Pa(INU,IMU,M) + Pb(INU,IMU,M)
            
              Fa(IMU,INU,M) = Fa(IMU,INU,M) + Rho1*Xij
              IF(IMU.NE.INU) Fa(INU,IMU,M) = Fa(INU,IMU,M) + Rho1*Xij
              IF(I.NE.J) THEN
                Fa(IRO,ISI,M) = Fa(IRO,ISI,M) + Rho2*Xij
                IF(IRO.NE.ISI)  Fa(ISI,IRO,M) = Fa(ISI,IRO,M) + Rho2*Xij
              ENDIF
           ENDIF

           !EXCHANGE PART.
           !(MU,NU | RO,SI)
           Fa(IMU,IRO,M) = Fa(IMU,IRO,M) + Pa(INU,ISI,M)*Yij
           IF(IMU.NE.IRO .OR. INU.NE.ISI)  & 
            Fa(IRO,IMU,M) = Fa(IRO,IMU,M) + Pa(ISI,INU,M)*Yij
           IF(IRO.NE.ISI) THEN
             !(MU,NU | SI,RO)
             Fa(IMU,ISI,M) = Fa(IMU,ISI,M) + Pa(INU,IRO,M)*Yij
             IF(ISI.NE.IMU .OR. INU.NE.IRO) &
              Fa(ISI,IMU,M) = Fa(ISI,IMU,M) + Pa(IRO,INU,M)*Yij
           ENDIF
           IF(IMU.NE.INU ) THEN
             !(NU,MU | RO,SI)
             IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
             Fa(INU,IRO,M) = Fa(INU,IRO,M) + Pa(IMU,ISI,M)*Yij
             IF(INU.NE.IRO .OR. IMU.NE.ISI) &
              Fa(IRO,INU,M) = Fa(IRO,INU,M) + Pa(ISI,IMU,M)*Yij
             ENDIF
             IF(IRO.NE.ISI) THEN
               !(NU,MU | SI,RO)
               Fa(INU,ISI,M) = Fa(INU,ISI,M) + Pa(IMU,IRO,M)*Yij
               IF(INU.NE.ISI .OR. IMU.NE.IRO) &
                Fa(ISI,INU,M) = Fa(ISI,INU,M) + Pa(IRO,IMU,M)*Yij
             ENDIF
           ENDIF

           IF(NDEN.EQ.2) THEN
             !COULOMB PART.
             !(MU,NU | RO,SI) + (MU,NU | SI,RO)  
             !(NU,MU | RO,SI) + (NU,MU | SI,RO)
             IF(I1N.EQ.1) THEN
               Fb(IMU,INU,M) = Fb(IMU,INU,M) + Rho1*Xij
               IF(IMU.NE.INU)  &
                Fb(INU,IMU,M) = Fb(INU,IMU,M) + Rho1*Xij
               IF(I.NE.J) THEN
                 Fb(IRO,ISI,M) = Fb(IRO,ISI,M) + Rho2*Xij
                 IF(IRO.NE.ISI)  &
                  Fb(ISI,IRO,M) = Fb(ISI,IRO,M) + Rho2*Xij
               ENDIF
             ENDIF

             !EXCHANGE PART.
             !(MU,NU | RO,SI)
             Fb(IMU,IRO,M) = Fb(IMU,IRO,M) + Pb(INU,ISI,M)*Yij
             IF(IMU.NE.IRO .OR. INU.NE.ISI) &
              Fb(IRO,IMU,M) = Fb(IRO,IMU,M) + Pb(ISI,INU,M)*Yij
             IF(IRO.NE.ISI) THEN
               !(MU,NU | SI,RO)
               Fb(IMU,ISI,M) = Fb(IMU,ISI,M) + Pb(INU,IRO,M)*Yij
               IF(ISI.NE.IMU .OR. INU.NE.IRO) &
                Fb(ISI,IMU,M) = Fb(ISI,IMU,M) + Pb(IRO,INU,M)*Yij
             ENDIF
             IF(IMU.NE.INU ) THEN
               !(NU,MU | RO,SI)
               IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
               Fb(INU,IRO,M) = Fb(INU,IRO,M) + Pb(IMU,ISI,M)*Yij
               IF(INU.NE.IRO .OR. IMU.NE.ISI) &
                Fb(IRO,INU,M) = Fb(IRO,INU,M) + Pb(ISI,IMU,M)*Yij
               ENDIF
               IF(IRO.NE.ISI) THEN
                 !(NU,MU | SI,RO)
                 Fb(INU,ISI,M) = Fb(INU,ISI,M) + Pb(IMU,IRO,M)*Yij
                 IF(INU.NE.ISI .OR. IMU.NE.IRO) &
                  Fb(ISI,INU,M) = Fb(ISI,INU,M) + Pb(IRO,IMU,M)*Yij
               ENDIF
             ENDIF

           ENDIF ! Beta Spin

        ENDDO !M
     ENDDO !J
  ENDDO !I

  RETURN
END

SUBROUTINE MAKEJK2(Ja,Jb,Ka,Kb,Pa,Pb,W,NL,LM3,LM6,NDEN,I1n,Cj,Ck)
!*
!BUILD CIS FOCK-LIKE MATRIX.
!*
!REFERENCE OCCUPATION NUMBERS OF DOUBLY OCCUPIED (CLOSED-SHELL)
!ORBITALS ARE SET TO TWO, AND ALL OTHERS ARE SET TO ZERO.
!*
!NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
!Pa/Pb(LM3,LM3)  DESITY MATRIX(I).
!W(LM6,LM6)      TWO-ELECTRON INTEGRALS(I).
!Fa/Fb(LM3,LM3)      CIS FOCK-LIKE MATRIX(O).
  USE LIMIT, ONLY: LMI
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  COMMON /CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
  COMMON /FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
  REAL*8 Ja(LM3,LM3,*),Jb(LM3,LM3,*),Ka(LM3,LM3,*),Kb(LM3,LM3,*)
  REAL*8 Pa(LM3,LM3,*),Pb(LM3,LM3,*),W(LM6,LM6)
  INTEGER NL,LM3,LM6,NDEN
  REAL*8 Cj,Ck

  CALL VecInit(Ja,LM3*LM3*NL,0d0)
  CALL VecInit(Ka,LM3*LM3*NL,0d0)
  IF(NDEN.EQ.2) THEN
    CALL VecInit(Jb,LM3*LM3*NL,0d0)
    CALL VecInit(Kb,LM3*LM3*NL,0d0)
  ENDIF
  !IP1(I).GE.IP2(I)
  DO I=1,LM6
     IMU = IP1(I)
     INU = IP2(I)
     DO J=I,LM6
        IRO = IP1(J)
        ISI = IP2(J)
        Xij = W(I,J)*Cj
        Yij = -W(I,J)*Ck
        DO M=1,NL 
           !COULOMB PART.
           !(MU,NU | RO,SI) + (MU,NU | SI,RO)  
           !(NU,MU | RO,SI) + (NU,MU | SI,RO)
           Rho1 = Pa(IRO,ISI,M) !+ Pb(IRO,ISI,M)
           IF(IRO.NE.ISI)  Rho1 = Rho1 + Pa(ISI,IRO,M) !+ Pb(ISI,IRO,M)
           Rho2 = Pa(IMU,INU,M) !+ Pb(IMU,INU,M)
           IF(IMU.NE.INU)  Rho2 = Rho2 + Pa(INU,IMU,M) !+ Pb(INU,IMU,M)
           
           Ja(IMU,INU,M) = Ja(IMU,INU,M) + Rho1*Xij
           IF(IMU.NE.INU) Ja(INU,IMU,M) = Ja(INU,IMU,M) + Rho1*Xij
           IF(I.NE.J) THEN
             Ja(IRO,ISI,M) = Ja(IRO,ISI,M) + Rho2*Xij
             IF(IRO.NE.ISI)  Ja(ISI,IRO,M) = Ja(ISI,IRO,M) + Rho2*Xij
           ENDIF

           !EXCHANGE PART.
           !(MU,NU | RO,SI)
           Ka(IMU,IRO,M) = Ka(IMU,IRO,M) + Pa(INU,ISI,M)*Yij
           IF(IMU.NE.IRO .OR. INU.NE.ISI)  & 
            Ka(IRO,IMU,M) = Ka(IRO,IMU,M) + Pa(ISI,INU,M)*Yij
           IF(IRO.NE.ISI) THEN
             !(MU,NU | SI,RO)
             Ka(IMU,ISI,M) = Ka(IMU,ISI,M) + Pa(INU,IRO,M)*Yij
             IF(ISI.NE.IMU .OR. INU.NE.IRO) &
              Ka(ISI,IMU,M) = Ka(ISI,IMU,M) + Pa(IRO,INU,M)*Yij
           ENDIF
           IF(IMU.NE.INU ) THEN
             !(NU,MU | RO,SI)
             IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
             Ka(INU,IRO,M) = Ka(INU,IRO,M) + Pa(IMU,ISI,M)*Yij
             IF(INU.NE.IRO .OR. IMU.NE.ISI) &
              Ka(IRO,INU,M) = Ka(IRO,INU,M) + Pa(ISI,IMU,M)*Yij
             ENDIF
             IF(IRO.NE.ISI) THEN
               !(NU,MU | SI,RO)
               Ka(INU,ISI,M) = Ka(INU,ISI,M) + Pa(IMU,IRO,M)*Yij
               IF(INU.NE.ISI .OR. IMU.NE.IRO) &
                Ka(ISI,INU,M) = Ka(ISI,INU,M) + Pa(IRO,IMU,M)*Yij
             ENDIF
           ENDIF

           IF(NDEN.EQ.2) THEN
             !COULOMB PART.
             !(MU,NU | RO,SI) + (MU,NU | SI,RO)  
             !(NU,MU | RO,SI) + (NU,MU | SI,RO)
             Rho1 = Pb(IRO,ISI,M) !+ Pb(IRO,ISI,M)
             IF(IRO.NE.ISI)  Rho1 = Rho1 + Pb(ISI,IRO,M) !+ Pb(ISI,IRO,M)
             Rho2 = Pb(IMU,INU,M) !+ Pb(IMU,INU,M)
             IF(IMU.NE.INU)  Rho2 = Rho2 + Pb(INU,IMU,M) !+ Pb(INU,IMU,M)
             Jb(IMU,INU,M) = Jb(IMU,INU,M) + Rho1*Xij
             IF(IMU.NE.INU)  &
              Jb(INU,IMU,M) = Jb(INU,IMU,M) + Rho1*Xij
             IF(I.NE.J) THEN
               Jb(IRO,ISI,M) = Jb(IRO,ISI,M) + Rho2*Xij
               IF(IRO.NE.ISI)  &
                Jb(ISI,IRO,M) = Jb(ISI,IRO,M) + Rho2*Xij
             ENDIF

             !EXCHANGE PART.
             !(MU,NU | RO,SI)
             Kb(IMU,IRO,M) = Kb(IMU,IRO,M) + Pb(INU,ISI,M)*Yij
             IF(IMU.NE.IRO .OR. INU.NE.ISI) &
              Kb(IRO,IMU,M) = Kb(IRO,IMU,M) + Pb(ISI,INU,M)*Yij
             IF(IRO.NE.ISI) THEN
               !(MU,NU | SI,RO)
               Kb(IMU,ISI,M) = Kb(IMU,ISI,M) + Pb(INU,IRO,M)*Yij
               IF(ISI.NE.IMU .OR. INU.NE.IRO) &
                Kb(ISI,IMU,M) = Kb(ISI,IMU,M) + Pb(IRO,INU,M)*Yij
             ENDIF
             IF(IMU.NE.INU ) THEN
               !(NU,MU | RO,SI)
               IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
               Kb(INU,IRO,M) = Kb(INU,IRO,M) + Pb(IMU,ISI,M)*Yij
               IF(INU.NE.IRO .OR. IMU.NE.ISI) &
                Kb(IRO,INU,M) = Kb(IRO,INU,M) + Pb(ISI,IMU,M)*Yij
               ENDIF
               IF(IRO.NE.ISI) THEN
                 !(NU,MU | SI,RO)
                 Kb(INU,ISI,M) = Kb(INU,ISI,M) + Pb(IMU,IRO,M)*Yij
                 IF(INU.NE.ISI .OR. IMU.NE.IRO) &
                  Kb(ISI,INU,M) = Kb(ISI,INU,M) + Pb(IRO,IMU,M)*Yij
               ENDIF
             ENDIF

           ENDIF ! Beta Spin

        ENDDO !M
     ENDDO !J
  ENDDO !I

  RETURN
END

SUBROUTINE MAKEJK(Ja,Jb,Ka,Kb,Pa,Pb,W,NL,LM3,LM6,NDEN,I1n,Cj,Ck)
!*
!BUILD CIS FOCK-LIKE MATRIX.
!*
!REFERENCE OCCUPATION NUMBERS OF DOUBLY OCCUPIED (CLOSED-SHELL)
!ORBITALS ARE SET TO TWO, AND ALL OTHERS ARE SET TO ZERO.
!*
!NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
!Pa/Pb(LM3,LM3)  DESITY MATRIX(I).
!W(LM6,LM6)      TWO-ELECTRON INTEGRALS(I).
!Fa/Fb(LM3,LM3)  CIS FOCK-LIKE MATRIX(O).
  USE LIMIT, ONLY: LMI
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  COMMON /CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
  COMMON /FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
  REAL*8  Ja(NL,LM3,LM3),Jb(NL,LM3,LM3),Ka(NL,LM3,LM3),Kb(NL,LM3,LM3)
  REAL*8  Pa(NL,LM3,LM3),Pb(NL,LM3,LM3),W(LM6,LM6)
  INTEGER NL,LM3,LM6,NDEN
  REAL*8  Cj,Ck

  CALL VecInit(Ja,LM3*LM3*NL,0d0)
  CALL VecInit(Ka,LM3*LM3*NL,0d0)
  IF(NDEN.EQ.2) THEN
    CALL VecInit(Jb,LM3*LM3*NL,0d0)
    CALL VecInit(Kb,LM3*LM3*NL,0d0)
  ENDIF
  CALL MATTRANS2(Pa,LM3*LM3,NL)
  IF(NDEN.EQ.2) CALL MATTRANS2(Pb,LM3*LM3,NL)
  !IP1(I).GE.IP2(I)
  DO I=1,LM6
     IMU = IP1(I)
     INU = IP2(I)
     DO J=I,LM6
        IRO = IP1(J)
        ISI = IP2(J)
        Xij = W(I,J)*Cj
        Yij = -W(I,J)*Ck
        DO M=1,NL 
           !COULOMB PART.
           !(MU,NU | RO,SI) + (MU,NU | SI,RO)  
           !(NU,MU | RO,SI) + (NU,MU | SI,RO)
           Rho1 = Pa(M,IRO,ISI) !+ Pb(M,IRO,ISI)
           IF(IRO.NE.ISI)  Rho1 = Rho1 + Pa(M,ISI,IRO) !+ Pb(M,ISI,IRO)
           Rho1 = Rho1*Xij
           Ja(M,IMU,INU) = Ja(M,IMU,INU) + Rho1
           IF(IMU.NE.INU) Ja(M,INU,IMU) = Ja(M,INU,IMU) + Rho1

           IF(I.NE.J) THEN
             Rho2 = Pa(M,IMU,INU) !+ Pb(M,IMU,INU)
             IF(IMU.NE.INU)  Rho2 = Rho2 + Pa(M,INU,IMU) !+ Pb(M,INU,IMU)
             Rho2 = Rho2*Xij
             Ja(M,IRO,ISI) = Ja(M,IRO,ISI) + Rho2
             IF(IRO.NE.ISI)  Ja(M,ISI,IRO) = Ja(M,ISI,IRO) + Rho2
           ENDIF

           !EXCHANGE PART.
           !(MU,NU | RO,SI)
           Ka(M,IMU,IRO) = Ka(M,IMU,IRO) + Pa(M,INU,ISI)*Yij
           IF(IMU.NE.IRO .OR. INU.NE.ISI)  & 
            Ka(M,IRO,IMU) = Ka(M,IRO,IMU) + Pa(M,ISI,INU)*Yij
           IF(IRO.NE.ISI) THEN
             !(MU,NU | SI,RO)
             Ka(M,IMU,ISI) = Ka(M,IMU,ISI) + Pa(M,INU,IRO)*Yij
             IF(ISI.NE.IMU .OR. INU.NE.IRO) &
              Ka(M,ISI,IMU) = Ka(M,ISI,IMU) + Pa(M,IRO,INU)*Yij
           ENDIF
           IF(IMU.NE.INU ) THEN
             !(NU,MU | RO,SI)
             IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
             Ka(M,INU,IRO) = Ka(M,INU,IRO) + Pa(M,IMU,ISI)*Yij
             IF(INU.NE.IRO .OR. IMU.NE.ISI) &
              Ka(M,IRO,INU) = Ka(M,IRO,INU) + Pa(M,ISI,IMU)*Yij
             ENDIF
             IF(IRO.NE.ISI) THEN
               !(NU,MU | SI,RO)
               Ka(M,INU,ISI) = Ka(M,INU,ISI) + Pa(M,IMU,IRO)*Yij
               IF(INU.NE.ISI .OR. IMU.NE.IRO) &
                Ka(M,ISI,INU) = Ka(M,ISI,INU) + Pa(M,IRO,IMU)*Yij
             ENDIF
           ENDIF

           IF(NDEN.EQ.2) THEN
             !COULOMB PART.
             !(MU,NU | RO,SI) + (MU,NU | SI,RO)  
             !(NU,MU | RO,SI) + (NU,MU | SI,RO)
             Rho1 = Pb(M,IRO,ISI) !+ Pb(M,IRO,ISI)
             IF(IRO.NE.ISI)  Rho1 = Rho1 + Pb(M,ISI,IRO) !+ Pb(M,ISI,IRO)
             Rho2 = Pb(M,IMU,INU) !+ Pb(M,IMU,INU)
             IF(IMU.NE.INU)  Rho2 = Rho2 + Pb(M,INU,IMU) !+ Pb(M,INU,IMU)
             Jb(M,IMU,INU) = Jb(M,IMU,INU) + Rho1*Xij
             IF(IMU.NE.INU)  &
              Jb(M,INU,IMU) = Jb(M,INU,IMU) + Rho1*Xij
             IF(I.NE.J) THEN
               Jb(M,IRO,ISI) = Jb(M,IRO,ISI) + Rho2*Xij
               IF(IRO.NE.ISI)  &
                Jb(M,ISI,IRO) = Jb(M,ISI,IRO) + Rho2*Xij
             ENDIF

             !EXCHANGE PART.
             !(MU,NU | RO,SI)
             Kb(M,IMU,IRO) = Kb(M,IMU,IRO) + Pb(M,INU,ISI)*Yij
             IF(IMU.NE.IRO .OR. INU.NE.ISI) &
              Kb(M,IRO,IMU) = Kb(M,IRO,IMU) + Pb(M,ISI,INU)*Yij
             IF(IRO.NE.ISI) THEN
               !(MU,NU | SI,RO)
               Kb(M,IMU,ISI) = Kb(M,IMU,ISI) + Pb(M,INU,IRO)*Yij
               IF(ISI.NE.IMU .OR. INU.NE.IRO) &
                Kb(M,ISI,IMU) = Kb(M,ISI,IMU) + Pb(M,IRO,INU)*Yij
             ENDIF
             IF(IMU.NE.INU ) THEN
               !(NU,MU | RO,SI)
               IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
               Kb(M,INU,IRO) = Kb(M,INU,IRO) + Pb(M,IMU,ISI)*Yij
               IF(INU.NE.IRO .OR. IMU.NE.ISI) &
                Kb(M,IRO,INU) = Kb(M,IRO,INU) + Pb(M,ISI,IMU)*Yij
               ENDIF
               IF(IRO.NE.ISI) THEN
                 !(NU,MU | SI,RO)
                 Kb(M,INU,ISI) = Kb(M,INU,ISI) + Pb(M,IMU,IRO)*Yij
                 IF(INU.NE.ISI .OR. IMU.NE.IRO) &
                  Kb(M,ISI,INU) = Kb(M,ISI,INU) + Pb(M,IRO,IMU)*Yij
               ENDIF
             ENDIF

           ENDIF ! Beta Spin

        ENDDO !M
     ENDDO !J
  ENDDO !I
  CALL MATTRANS2(Ja,NL,LM3*LM3)
  CALL MATTRANS2(ka,NL,LM3*LM3)
  IF(NDEN.EQ.2) THEN
    CALL MATTRANS2(Jb,NL,LM3*LM3)
    CALL MATTRANS2(kb,NL,LM3*LM3)
  ENDIF

  RETURN
END

SUBROUTINE MAKEJK3(Ja,Jb,Ka,Kb,Pa,Pb,W,NL,LM3,LM6,NDEN,I1n,Cj,Ck)
!*
!BUILD CIS FOCK-LIKE MATRIX.
!*
!REFERENCE OCCUPATION NUMBERS OF DOUBLY OCCUPIED (CLOSED-SHELL)
!ORBITALS ARE SET TO TWO, AND ALL OTHERS ARE SET TO ZERO.
!*
!NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
!Pa/Pb(LM3,LM3)  DESITY MATRIX(I).
!W(LM6,LM6)      TWO-ELECTRON INTEGRALS(I).
!Fa/Fb(LM3,LM3)  CIS FOCK-LIKE MATRIX(O).
  USE LIMIT, ONLY: LMI
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  COMMON /CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
  COMMON /FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
  REAL*8  Ja(NL,LM3,LM3),Jb(NL,LM3,LM3),Ka(NL,LM3,LM3),Kb(NL,LM3,LM3)
  REAL*8  Pa(NL,LM3,LM3),Pb(NL,LM3,LM3),W(LM6,LM6)
  INTEGER NL,LM3,LM6,NDEN
  REAL*8  Cj,Ck

  CALL VecInit(Ja,LM3*LM3*NL,0d0)
  CALL VecInit(Ka,LM3*LM3*NL,0d0)
  IF(NDEN.EQ.2) THEN
    CALL VecInit(Jb,LM3*LM3*NL,0d0)
    CALL VecInit(Kb,LM3*LM3*NL,0d0)
  ENDIF
  CALL MATTRANS2(Pa,LM3*LM3,NL)
  IF(NDEN.EQ.2) CALL MATTRANS2(Pb,LM3*LM3,NL)
  !IP1(I).GE.IP2(I)
  DO I=1,LM6
     IMU = IP1(I)
     INU = IP2(I)
     DO J=I,LM6
        IRO = IP1(J)
        ISI = IP2(J)
        Xij = W(I,J)*Cj
        Yij = -W(I,J)*Ck
        DO M=1,NL 
           !COULOMB PART.
           !(MU,NU | RO,SI) + (MU,NU | SI,RO)  
           !(NU,MU | RO,SI) + (NU,MU | SI,RO)
           Rho1 = Pa(M,IRO,ISI) !+ Pb(M,IRO,ISI)
           IF(IRO.NE.ISI)  Rho1 = Rho1 + Pa(M,ISI,IRO) !+ Pb(M,ISI,IRO)
           Rho1 = Rho1*Xij
           Ja(M,IMU,INU) = Ja(M,IMU,INU) + Rho1
           IF(IMU.NE.INU) Ja(M,INU,IMU) = Ja(M,INU,IMU) + Rho1

           IF(I.NE.J) THEN
             Rho2 = Pa(M,IMU,INU) !+ Pb(M,IMU,INU)
             IF(IMU.NE.INU)  Rho2 = Rho2 + Pa(M,INU,IMU) !+ Pb(M,INU,IMU)
             Rho2 = Rho2*Xij
             Ja(M,IRO,ISI) = Ja(M,IRO,ISI) + Rho2
             IF(IRO.NE.ISI)  Ja(M,ISI,IRO) = Ja(M,ISI,IRO) + Rho2
           ENDIF

           !EXCHANGE PART.
           !(MU,NU | RO,SI)
           Ka(M,IMU,IRO) = Ka(M,IMU,IRO) + Pa(M,INU,ISI)*Yij
           IF(IMU.NE.IRO .OR. INU.NE.ISI)  & 
            Ka(M,IRO,IMU) = Ka(M,IRO,IMU) + Pa(M,ISI,INU)*Yij
           IF(IRO.NE.ISI) THEN
             !(MU,NU | SI,RO)
             Ka(M,IMU,ISI) = Ka(M,IMU,ISI) + Pa(M,INU,IRO)*Yij
             IF(ISI.NE.IMU .OR. INU.NE.IRO) &
              Ka(M,ISI,IMU) = Ka(M,ISI,IMU) + Pa(M,IRO,INU)*Yij
           ENDIF
           IF(IMU.NE.INU ) THEN
             !(NU,MU | RO,SI)
             IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
             Ka(M,INU,IRO) = Ka(M,INU,IRO) + Pa(M,IMU,ISI)*Yij
             IF(INU.NE.IRO .OR. IMU.NE.ISI) &
              Ka(M,IRO,INU) = Ka(M,IRO,INU) + Pa(M,ISI,IMU)*Yij
             ENDIF
             IF(IRO.NE.ISI) THEN
               !(NU,MU | SI,RO)
               Ka(M,INU,ISI) = Ka(M,INU,ISI) + Pa(M,IMU,IRO)*Yij
               IF(INU.NE.ISI .OR. IMU.NE.IRO) &
                Ka(M,ISI,INU) = Ka(M,ISI,INU) + Pa(M,IRO,IMU)*Yij
             ENDIF
           ENDIF

           IF(NDEN.EQ.2) THEN
             !COULOMB PART.
             !(MU,NU | RO,SI) + (MU,NU | SI,RO)  
             !(NU,MU | RO,SI) + (NU,MU | SI,RO)
             Rho1 = Pb(M,IRO,ISI) !+ Pb(M,IRO,ISI)
             IF(IRO.NE.ISI)  Rho1 = Rho1 + Pb(M,ISI,IRO) !+ Pb(M,ISI,IRO)
             Rho2 = Pb(M,IMU,INU) !+ Pb(M,IMU,INU)
             IF(IMU.NE.INU)  Rho2 = Rho2 + Pb(M,INU,IMU) !+ Pb(M,INU,IMU)
             Jb(M,IMU,INU) = Jb(M,IMU,INU) + Rho1*Xij
             IF(IMU.NE.INU)  &
              Jb(M,INU,IMU) = Jb(M,INU,IMU) + Rho1*Xij
             IF(I.NE.J) THEN
               Jb(M,IRO,ISI) = Jb(M,IRO,ISI) + Rho2*Xij
               IF(IRO.NE.ISI)  &
                Jb(M,ISI,IRO) = Jb(M,ISI,IRO) + Rho2*Xij
             ENDIF

             !EXCHANGE PART.
             !(MU,NU | RO,SI)
             Kb(M,IMU,IRO) = Kb(M,IMU,IRO) + Pb(M,INU,ISI)*Yij
             IF(IMU.NE.IRO .OR. INU.NE.ISI) &
              Kb(M,IRO,IMU) = Kb(M,IRO,IMU) + Pb(M,ISI,INU)*Yij
             IF(IRO.NE.ISI) THEN
               !(MU,NU | SI,RO)
               Kb(M,IMU,ISI) = Kb(M,IMU,ISI) + Pb(M,INU,IRO)*Yij
               IF(ISI.NE.IMU .OR. INU.NE.IRO) &
                Kb(M,ISI,IMU) = Kb(M,ISI,IMU) + Pb(M,IRO,INU)*Yij
             ENDIF
             IF(IMU.NE.INU ) THEN
               !(NU,MU | RO,SI)
               IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
               Kb(M,INU,IRO) = Kb(M,INU,IRO) + Pb(M,IMU,ISI)*Yij
               IF(INU.NE.IRO .OR. IMU.NE.ISI) &
                Kb(M,IRO,INU) = Kb(M,IRO,INU) + Pb(M,ISI,IMU)*Yij
               ENDIF
               IF(IRO.NE.ISI) THEN
                 !(NU,MU | SI,RO)
                 Kb(M,INU,ISI) = Kb(M,INU,ISI) + Pb(M,IMU,IRO)*Yij
                 IF(INU.NE.ISI .OR. IMU.NE.IRO) &
                  Kb(M,ISI,INU) = Kb(M,ISI,INU) + Pb(M,IRO,IMU)*Yij
               ENDIF
             ENDIF

           ENDIF ! Beta Spin

        ENDDO !M
     ENDDO !J
  ENDDO !I
  CALL MATTRANS2(Ja,NL,LM3*LM3)
  CALL MATTRANS2(ka,NL,LM3*LM3)
  IF(NDEN.EQ.2) THEN
    CALL MATTRANS2(Jb,NL,LM3*LM3)
    CALL MATTRANS2(kb,NL,LM3*LM3)
  ENDIF

  RETURN
END

SUBROUTINE MAKEFAOx2(Fa,Fb,Pa,Pb,W,NL,LM3,LM6,NDEN,I1n,Cj,Ck)
!*
!BUILD CIS FOCK-LIKE MATRIX (NL>4).
!*
!Pa/Pb(LM3,LM3)  DESITY MATRIX(I).
!W(LM6,LM6)      TWO-ELECTRON INTEGRALS(I).
!Fa/Fb(LM3,LM3)  CIS FOCK-LIKE MATRIX(O).
  USE LIMIT, ONLY: LMI
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  COMMON  /CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
  COMMON  /FINDX1/ IP(LMI),IP1(LMI),IP2(LMI)
  REAL*8  Pa(NL,*),Pb(NL,*),Fa(NL,*),Fb(NL,*),W(LM6,LM6)
  REAL*8  Cj,Ck
  INTEGER NL,LM3,LM6,NDEN

  CALL VecInit(Fa,LM3*LM3*NL,0d0)
  IF(NDEN.EQ.2) CALL VecInit(Fb,LM3*LM3*NL,0d0)

  CALL MATTRANS2(Pa,LM3*LM3,NL)
  IF(NDEN.EQ.2) CALL MATTRANS2(Pb,LM3*LM3,NL)

  !IP1(I).GE.IP2(I)
  DO I=1,LM6
     IMU = IP1(I)
     INU = IP2(I)
     IMN = IMU+(INU-1)*LM3
     INM = INU+(IMU-1)*LM3
     DO J=I,LM6
        IRO = IP1(J)
        ISI = IP2(J)
        Xij = W(I,J)*Cj
        Yij = -W(I,J)*Ck
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
        DO M=1,NL 
           !COULOMB PART.
           !(MU,NU | RO,SI) + (MU,NU | SI,RO)  
           !(NU,MU | RO,SI) + (NU,MU | SI,RO)
           IF(I1n.EQ.1) THEN
              Rho1 = Pa(M,IRS) + Pb(M,IRS)
              IF(IRO.NE.ISI)  Rho1 = Rho1 + Pa(M,ISR) + Pb(M,ISR)
              Rho1x = Rho1*Xij
              Fa(M,IMN) = Fa(M,IMN) + Rho1x
              IF(IMU.NE.INU) Fa(M,INM) = Fa(M,INM) + Rho1x

              IF(I.NE.J) THEN
                Rho2 = Pa(M,IMN) + Pb(M,IMN)
                IF(IMU.NE.INU)  Rho2 = Rho2 + Pa(M,INM) + Pb(M,INM)
                Rho2x = Rho2*Xij
                Fa(M,IRS) = Fa(M,IRS) + Rho2x
                IF(IRO.NE.ISI)  Fa(M,ISR) = Fa(M,ISR) + Rho2x
              ENDIF
           ENDIF

           !EXCHANGE PART.
           !(MU,NU | RO,SI)
           Fa(M,IMR) = Fa(M,IMR) + Pa(M,INS)*Yij
           IF(IMU.NE.IRO .OR. INU.NE.ISI)  & 
            Fa(M,IRM) = Fa(M,IRM) + Pa(M,ISN)*Yij
           IF(IRO.NE.ISI) THEN
             !(MU,NU | SI,RO)
             Fa(M,IMS) = Fa(M,IMS) + Pa(M,INR)*Yij
             IF(ISI.NE.IMU .OR. INU.NE.IRO) &
              Fa(M,ISM) = Fa(M,ISM) + Pa(M,IRN)*Yij
           ENDIF
           IF(IMU.NE.INU ) THEN
             !(NU,MU | RO,SI)
             IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
             Fa(M,INR) = Fa(M,INR) + Pa(M,IMS)*Yij
             IF(INU.NE.IRO .OR. IMU.NE.ISI) &
              Fa(M,IRN) = Fa(M,IRN) + Pa(M,ISM)*Yij
             ENDIF
             IF(IRO.NE.ISI) THEN
               !(NU,MU | SI,RO)
               Fa(M,INS) = Fa(M,INS) + Pa(M,IMR)*Yij
               IF(INU.NE.ISI .OR. IMU.NE.IRO) &
                Fa(M,ISN) = Fa(M,ISN) + Pa(M,IRM)*Yij
             ENDIF
           ENDIF

           IF(NDEN.EQ.2) THEN
             !COULOMB PART.
             !(MU,NU | RO,SI) + (MU,NU | SI,RO)  
             !(NU,MU | RO,SI) + (NU,MU | SI,RO)
             IF(I1N.EQ.1) THEN
               Fb(M,IMN) = Fb(M,IMN) + Rho1x
               IF(IMU.NE.INU)  &
                Fb(M,INM) = Fb(M,INM) + Rho1x
               IF(I.NE.J) THEN
                 Fb(M,IRS) = Fb(M,IRS) + Rho2x
                 IF(IRO.NE.ISI)  &
                  Fb(M,ISR) = Fb(M,ISR) + Rho2x
               ENDIF
             ENDIF

             !EXCHANGE PART.
             !(MU,NU | RO,SI)
             Fb(M,IMR) = Fb(M,IMR) + Pb(M,INS)*Yij
             IF(IMU.NE.IRO .OR. INU.NE.ISI) &
              Fb(M,IRM) = Fb(M,IRM) + Pb(M,ISN)*Yij
             IF(IRO.NE.ISI) THEN
               !(MU,NU | SI,RO)
               Fb(M,IMS) = Fb(M,IMS) + Pb(M,INR)*Yij
               IF(ISI.NE.IMU .OR. INU.NE.IRO) &
                Fb(M,ISM) = Fb(M,ISM) + Pb(M,IRN)*Yij
             ENDIF
             IF(IMU.NE.INU ) THEN
               !(NU,MU | RO,SI)
               IF(INU.NE.ISI .OR. IMU.NE.IRO) THEN
               Fb(M,INR) = Fb(M,INR) + Pb(M,IMS)*Yij
               IF(INU.NE.IRO .OR. IMU.NE.ISI) &
                Fb(M,IRN) = Fb(M,IRN) + Pb(M,ISM)*Yij
               ENDIF
               IF(IRO.NE.ISI) THEN
                 !(NU,MU | SI,RO)
                 Fb(M,INS) = Fb(M,INS) + Pb(M,IMR)*Yij
                 IF(INU.NE.ISI .OR. IMU.NE.IRO) &
                  Fb(M,ISN) = Fb(M,ISN) + Pb(M,IRM)*Yij
               ENDIF
             ENDIF

           ENDIF ! Beta Spin

        ENDDO !M
     ENDDO !J
  ENDDO !I

  CALL MATTRANS2(Fa,NL,LM3*LM3)
  CALL MATTRANS2(Pa,NL,LM3*LM3)
  IF(NDEN.EQ.2) THEN
    CALL MATTRANS2(Fb,NL,LM3*LM3)
    CALL MATTRANS2(Pb,NL,LM3*LM3)
  ENDIF

  RETURN
END
