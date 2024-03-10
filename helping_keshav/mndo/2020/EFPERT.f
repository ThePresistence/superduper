      SUBROUTINE EFPERT (H,LM4)
C     *
C     ADD ELECTRIC FIELD PERTURBATION MATRIX ELEMENTS TO THE CORE
C     HAMILTONIAN H(LM4).
C     *
      USE LIMIT, ONLY: LM1, LMX, LMZ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
     ./FIELD2/ FIFI(3)
     ./INDEX / INDX(LMX)
     ./MULTIP/ DD(6,0:LMZ),PO(9,0:LMZ)
      DIMENSION H(LM4)
C *** INITIALIZATION.
      RSQRT3 = ONE/SQRT(THREE)
C *** LOOP OVER ALL ATOMS.
      DO 40 I=1,NUMAT
      NI     = NAT(I)
      IA     = NFIRST(I)
      IORBS  = NLAST(I)-IA+1
      IS     = INDX(IA)+IA
      K      = IA-1
C *** PERTURBATION CONTRIBUTIONS TO DIAGONAL ELEMENTS (IN EV)
C     FIFI(IC)     COMPONENTS OF THE APPLIED FIELD (IN A.U.)
C     COORD(IC,I)  CARTESIAN COORDINATES OF ATOM I (IN ANGSTROM)
      HDIAG  = ZERO
      DO 10 IC=1,3
      HDIAG  = HDIAG + COORD(IC,I)*FIFI(IC)
   10 CONTINUE
      HDIAG  = HDIAG / A0*EV
C *** S-S DIAGONAL ELEMENT.
      H(IS)  = H(IS) + HDIAG
C *** P-P DIAGONAL ELEMENTS.
      IF(IORBS.GE.4) THEN
         DO 20 J=IA+1,IA+3
         LL     = INDX(J)+J
         H(LL)  = H(LL) + HDIAG
   20    CONTINUE
C *** S-P OFF-DIAGONAL ELEMENTS.
C        Rsp is the radial sp dipole integral.
C        H(s,px) = H(s,px) + E(x)*Rsp
C        H(s,py) = H(s,py) + E(y)*Rsp
C        H(s,pz) = H(s,pz) + E(z)*Rsp
         ISX   = IS+1+K
         ISY   = ISX+2+K
         ISZ   = ISY+3+K
         H(ISX) = H(ISX) + FIFI(1)*DD(2,NI) * EV
         H(ISY) = H(ISY) + FIFI(2)*DD(2,NI) * EV
         H(ISZ) = H(ISZ) + FIFI(3)*DD(2,NI) * EV
C *** D-D DIAGONAL ELEMENTS.
         IF(IORBS.GE.9) THEN
            DO 30 J=IA+4,IA+8
            LL     = INDX(J)+J
            H(LL)  = H(LL) + HDIAG
   30       CONTINUE
C *** P-D OFF-DIAGONAL ELEMENTS.
            IDX2Y2  = ISZ+3+K
            IDXZ    = IDX2Y2+5+K
            IDZ2    = IDXZ+6+K
            IDYZ    = IDZ2+7+K
            IDXY    = IDYZ+8+K
C           X-components
C           Rpd is the radial pd dipole integral.
C           H(px,dz2)    = H(px,dz2)   - E(x)*Rpd/sqrt(3)
C           H(px,dx2-y2) = H(px,dx2-y2)+ E(x)*Rpd
C           H(py,dxy)    = H(py,dxy)   + E(x)*Rpd
C           H(pz,dxz)    = H(pz,dxz)   + E(x)*Rpd
            XF  = FIFI(1)*DD(5,NI) * EV
            H(IDX2Y2+2)  = H(IDX2Y2+2) + XF
            H(IDXZ+4)    = H(IDXZ+4)   + XF
            H(IDZ2+2)    = H(IDZ2+2)   - XF*RSQRT3
            H(IDXY+3)    = H(IDXY+3)   + XF
C           Y-components
C           H(py,dz2)    = H(py,dz2)   - E(y)*Rpd/sqrt(3)
C           H(py,dx2-y2) = H(py,dx2-y2)- E(y)*Rpd
C           H(px,dxy)    = H(px,dxy)   + E(y)*Rpd
C           H(pz,dyz)    = H(pz,dyz)   + E(y)*Rpd
            YF  = FIFI(2)*DD(5,NI) * EV
            H(IDX2Y2+3)  = H(IDX2Y2+3) - YF
            H(IDZ2+3)    = H(IDZ2+3)   - YF*RSQRT3
            H(IDYZ+4)    = H(IDYZ+4)   + YF
            H(IDXY+2)    = H(IDXY+2)   + YF
C           Z-components
C           H(pz,dz2)    = H(pz,dz2)   + E(z)*Rpd*2/sqrt(3)
C           H(px,dxz)    = H(px,dxz)   + E(z)*Rpd
C           H(py,dyz)    = H(py,dyz)   + E(z)*Rpd
            ZF  = FIFI(3)*DD(5,NI) * EV
            H(IDXZ+2)    = H(IDXZ+2)   + ZF
            H(IDYZ+3)    = H(IDYZ+3)   + ZF
            H(IDZ2+4)    = H(IDZ2+4)   + (ZF+ZF)*RSQRT3
         ENDIF
      ENDIF
   40 CONTINUE
      RETURN
      END
