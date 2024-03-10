MODULE RBM

  USE DataTypes
  USE CoulombMod
  implicit none

  !  Proportionality constant used to get the minimum number of particles
  REAL(SP), PRIVATE, PARAMETER :: PCMNP=0.50_SP
  !   REAL(SP), PRIVATE, PARAMETER :: PCMNP=0.10_SP

  !------------------------------------------------------------------------------
  !  STRUCTURE FOR STORING THE INFORMATION ABOUT A SET OF PARTICLES
  !------------------------------------------------------------------------------
  TYPE :: PARTICLE
     REAL(SP), POINTER :: Crd(:,:) => Null()     ! List of particles
     !  -- Crd(1,i) is the charge for the i-th particle
     !  -- Crd(2,i) is X coordinate for the i-th particle
     !  -- Crd(3,i) is Y coordinate for the i-th particle
     !  -- Crd(4,i) is Z coordinate for the i-th particle
     !  -- Crd(5,i) is Gaussian exponent factor for the i-th particle
     REAL(SP), POINTER :: pot(:) =>NULL()
     REAL(SP), POINTER :: Gradn(:,:) => Null()      ! fields
     !  -- Grad(1,i) is the X component of the potential gradient at the i-th particle
     !  -- Grad(2,i) is the Y component of the potential gradient at the i-th particle
     !  -- Grad(3,i) is the Z component of the potential gradient at the i-th particle
     INTEGER(I4B), POINTER :: KEY(:)   => Null() ! "Key" used to restore particle order
     REAL(SP) :: COFC(3) = 0.0_sp  ! Center of charge
     REAL(SP) :: MOFI(3,3) = 0.0_sp! Moment of inertia matrix
     REAL(SP) :: SPL(3)   = 0.0_sp ! Vector perpENDicular to the splitting plane
     REAL(SP), POINTER :: QLM(:,:,:)   => Null() ! Coefficients of the multipole 
     ! expansions. The coefficients are stored as follows:
     !   Q^+,    "QLM(1,l,m)"   is "Q^+_{lm}"   (multipole expansion)
     !           "QLM(2,l,m)"   is "Q^+_{lm}^x" (X derivative)
     !           "QLM(3,l,m)"   is "Q^+_{lm}^y" (Y derivative)
     !           "QLM(4,l,m)"   is "Q^+_{lm}^z" (Z derivative)
     !   Q^-,    "QLM(1,l-m,l)" is "Q^-_{lm}"   (multipole expansion)
     !           "QLM(2,l-m,l)" is "Q^-_{lm}^x" (X derivative)
     !           "QLM(3,l-m,l)" is "Q^-_{lm}^y" (Y derivative)
     !           "QLM(4,l-m,l)" is "Q^-_{lm}^z" (Z derivative)
     LOGICAL(LGD) :: GETQLM  = .true. ! To decide whether to compute "QLM" or not
     LOGICAL(LGD) :: GETSPL = .true.  ! To decide whether to compute "COFC" and "MOFI" or not
     LOGICAL(LGD) :: SPLIT = .true.  ! To decide whether to split the set or not
     INTEGER(I4B) :: N1 = -1      ! Number of particles in the first "child"
  
  END TYPE PARTICLE

  PRIVATE
  PUBLIC :: RecursiveBisectionMethod, DirectCoulomb

  INTERFACE RecursiveBisectionMethod
     MODULE PROCEDURE TLRB1_NOGRAD,TLRB1_GRAD, TLRB2_NOGRAD,TLRB2_GRAD
  END INTERFACE
  INTERFACE DirectCoulomb
     MODULE PROCEDURE DirectCoulomb1, DirectCoulomb_GRAD,DirectCoulomb_NOGRAD
  END INTERFACE

CONTAINS

  !*******************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Approximate and FAST evaluation of particle-particle interactions (Coulomb  !
  ! potential and electric fields) by the RECURSIVE BISECTION method.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE TLRB1_GRAD(Crd,pot,Grad,WS,L)
    IMPLICIT NONE
    REAL(SP), INTENT(IN) ::     WS    ! Well-separatedness PARAMETER
    INTEGER(I4B), INTENT(IN) :: L     ! Order of the multipole expansion
    !------------------------------------------------------------------------------
    !  INPUT-OUTPUT:
    REAL(SP), INTENT(INOUT), TARGET :: Crd(:,:)    ! List of particles
    REAL(SP), INTENT(INOUT), TARGET :: Pot(:)
    REAL(SP), INTENT(INOUT), TARGET :: Grad(:,:)    ! fields
    INTEGER(I4B) :: i,np
    INTEGER(I4B), TARGET :: key(SIZE(Crd,2))
    TYPE(PARTICLE) :: set

    !  Initialization
    set%Crd=>Crd; set%pot => pot; set%key=>KEY
    set%getqlm=.TRUE.; set%getspl=.TRUE.; set%split=.TRUE.

    np=SIZE(Crd,2); key=(/ (i,i=1,np) /)
    !  Execute the RECURSIVE bisection method

    set%Gradn=>Grad
    CALL TLRB0_GRAD(WS,L,set,set)
    !  Restore the original ordering
    CALL RORDRB(1,np,key,Crd,pot,Grad)

  END SUBROUTINE TLRB1_GRAD

  SUBROUTINE TLRB1_NOGRAD(Crd,pot,WS,L)
    IMPLICIT NONE
    REAL(SP), INTENT(IN) ::     WS    ! Well-separatedness PARAMETER
    INTEGER(I4B), INTENT(IN) :: L     ! Order of the multipole expansion
    !------------------------------------------------------------------------------
    !  INPUT-OUTPUT:
    REAL(SP), INTENT(INOUT), TARGET :: Crd(:,:)    ! List of particles
    REAL(SP), INTENT(INOUT), TARGET :: Pot(:)
    INTEGER(I4B) :: i,np
    INTEGER(I4B), TARGET :: key(SIZE(Crd,2))
    TYPE(PARTICLE) :: set

    !  Initialization
    set%Crd=>Crd; set%pot => pot; set%key=>KEY
    set%getqlm=.TRUE.; set%getspl=.TRUE.; set%split=.TRUE.

    np=SIZE(Crd,2); key=(/ (i,i=1,np) /)
    !  Execute the RECURSIVE bisection method

    CALL TLRB0_NOGRAD(WS,L,set,set)
    !  Restore the original ordering
    CALL RORDRB(1,np,key,Crd,pot)
  END SUBROUTINE TLRB1_NOGRAD


  SUBROUTINE TLRB2_grad(Crd1,pot1,Grad1,Crd2,pot2,Grad2,WS,L)
    IMPLICIT NONE
    REAL(SP),  INTENT(IN) ::     WS    ! Well-separatedness PARAMETER
    INTEGER(I4B),  INTENT(IN) :: L     ! Order of the multipole expansion
    !------------------------------------------------------------------------------
    !  INPUT-OUTPUT:
    REAL(SP), INTENT(INOUT), TARGET :: Crd1(:,:), Crd2(:,:)    ! List of particles  
    REAL(SP), INTENT(INOUT), TARGET :: Pot1(:), Pot2(:)
    REAL(SP), INTENT(INOUT), TARGET :: Grad1(:,:), Grad2(:,:)    ! fields
    INTEGER(I4B) :: i,np1,np2
    INTEGER(I4B), TARGET :: key1(SIZE(Crd1,2))
    INTEGER(I4B), TARGET :: key2(SIZE(Crd2,2))
    TYPE(PARTICLE) :: set1, set2

    !  Initialization
    set1%Crd=>Crd1; set1%pot=>pot1; set1%key=>key1
    set2%Crd=>Crd2; set2%pot=>pot2; set2%key=>key2

    set1%getqlm=.TRUE.; set1%getspl=.TRUE.; set1%split=.TRUE.
    set2%getqlm=.TRUE.; set2%getspl=.TRUE.; set2%split=.TRUE.

    np1=SIZE(Crd1,2); key1=(/ (i,i=1,np1) /)
    np2=SIZE(Crd2,2); key2=(/ (i,i=1,np2) /)

    !  Execute the RECURSIVE bisection method
    set1%Gradn=>Grad1
    set2%Gradn=>Grad2
    CALL TLRB0_GRAD(WS,L,set1,set2)
    !  Restore the original ordering
    CALL RORDRB(1,np1,key1,Crd1,pot1,Grad1)
    CALL RORDRB(1,np2,key2,Crd2,pot2,Grad2)
  END SUBROUTINE TLRB2_GRAD

  SUBROUTINE TLRB2_NOGRAD(Crd1,pot1,Crd2,pot2,WS,L)
    IMPLICIT NONE
    REAL(SP),  INTENT(IN) ::     WS    ! Well-separatedness PARAMETER
    INTEGER(I4B),  INTENT(IN) :: L     ! Order of the multipole expansion
    !------------------------------------------------------------------------------
    !  INPUT-OUTPUT:
    REAL(SP), INTENT(INOUT), TARGET :: Crd1(:,:), Crd2(:,:)    ! List of particles  
    REAL(SP), INTENT(INOUT), TARGET :: Pot1(:), Pot2(:)
    INTEGER(I4B) :: i,np1,np2
    INTEGER(I4B), TARGET :: key1(SIZE(Crd1,2))
    INTEGER(I4B), TARGET :: key2(SIZE(Crd2,2))
    TYPE(PARTICLE) :: set1, set2

    !  Initialization
    set1%Crd=>Crd1; set1%pot=>pot1; set1%key=>key1
    set2%Crd=>Crd2; set2%pot=>pot2; set2%key=>key2

    set1%getqlm=.TRUE.; set1%getspl=.TRUE.; set1%split=.TRUE.
    set2%getqlm=.TRUE.; set2%getspl=.TRUE.; set2%split=.TRUE.

    np1=SIZE(Crd1,2); key1=(/ (i,i=1,np1) /)
    np2=SIZE(Crd2,2); key2=(/ (i,i=1,np2) /)

    CALL TLRB0_NOGRAD(WS,L,set1,set2)
    !  Restore the original ordering
    CALL RORDRB(1,np1,key1,Crd1,pot1)
    CALL RORDRB(1,np2,key2,Crd2,pot2)
  END SUBROUTINE TLRB2_NOGRAD
  


  !*******************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Approximate and FAST evaluation of particle-particle interactions (Coulomb  !
  ! potential and fields) by the RECURSIVE BISECTION method (auxiliary           !
  ! SUBROUTINE).                                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE TLRB0_NOGRAD(WS,L,SET1,SET2)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    !  INPUT:
    REAL(SP), INTENT(IN) ::     WS    ! Well-separatedness PARAMETER
    INTEGER(I4B), INTENT(IN) ::            L     ! Order of the multipole expansion
    !------------------------------------------------------------------------------
    !  INPUT-OUTPUT:
    TYPE(PARTICLE), INTENT(INOUT) :: SET1  ! First set of particles
    TYPE(PARTICLE), INTENT(INOUT) :: SET2  ! Second set of particles
    !------------------------------------------------------------------------------
    INTEGER(I4B) :: i,n1,n2,idec
    REAL(SP), TARGET :: qlm11(1,0:L,0:L),qlm12(1,0:L,0:L), &
         qlm21(1,0:L,0:L),qlm22(1,0:L,0:L)
    REAL(SP), POINTER :: Crd1(:,:),Pot1(:),Crd2(:,:),Pot2(:)
    LOGICAL(LGD) :: itself
    TYPE(PARTICLE) :: set11,set12,set21,set22

    !  Initialization
    Crd1=>SET1%Crd; Pot1=>SET1%Pot; Crd2=>SET2%Crd; Pot2=>SET2%Pot
    n1=SIZE(Crd1,2); n2=SIZE(Crd2,2)

    set11%getqlm=.TRUE.; set11%getspl=.TRUE.; set11%split=.TRUE.
    set12%getqlm=.TRUE.; set12%getspl=.TRUE.; set12%split=.TRUE.
    set21%getqlm=.TRUE.; set21%getspl=.TRUE.; set21%split=.TRUE.
    set22%getqlm=.TRUE.; set22%getspl=.TRUE.; set22%split=.TRUE.
    set11%qlm=>qlm11; set12%qlm=>qlm12; set21%qlm=>qlm21; set22%qlm=>qlm22

    itself=ASSOCIATED(Crd1,Crd2)

    !  Get decision (direct sum, multipole expansion, or splitting)
    idec=DMAKRB(WS,L,SET1,SET2)

    !  Execute the decision
    SELECT CASE (idec)

       !  DIRECT SUMMATION
    CASE (0)
       !WRITE(*,*) 'SUMMATION'
       CALL DirectCoulomb_NOGRAD(Crd1,Pot1,Crd2,Pot2,itself)
       !  MULTIPOLE EXPANSION
    CASE (1)

       !WRITE(*,*) 'MULTIPOLE EXPANSION'

       IF (SET1%GETQLM) THEN
          CALL GMEXRB_NOGRAD(L,SET1)
       END IF
       CALL EMEXRB_NOGRAD(L,SET1%QLM,SET1%COFC,SET2)

       IF (SET2%GETQLM) THEN
          CALL GMEXRB_NOGRAD(L,SET2)
       END IF
       CALL EMEXRB_NOGRAD(L,SET2%QLM,SET2%COFC,SET1)

       !  SPLITTING
    CASE (2)
       !  Split the first distribution
       IF (SET1%SPLIT) THEN

          SET1%SPL=LEIGRB(SET1%MOFI)
          CALL PSPLRB(SET1,set11,set12)

          set11%qlm=>qlm11; set12%qlm=>qlm12
       ELSE
          i=SET1%N1
          set11%Crd=>SET1%Crd(:,1:i); set12%Crd=>SET1%Crd(:,i+1:)
          set11%Pot=>SET1%Pot(1:i); set12%Pot=>SET1%Pot(i+1:)
          set11%KEY=>SET1%KEY(1:i); set12%KEY=>SET1%KEY(i+1:)
       END IF

       IF (itself) THEN

          !  Recursion (interaction of "SET1" with itself) 
          CALL TLRB0_NOGRAD(WS,L,set11,set11)     
          CALL TLRB0_NOGRAD(WS,L,set12,set12)     
          CALL TLRB0_NOGRAD(WS,L,set11,set12)
       ELSE
          !  Split the second distribution
          IF (SET2%SPLIT) THEN

             SET2%SPL=LEIGRB(SET2%MOFI)
             CALL PSPLRB(SET2,set21,set22)

             set21%qlm=>qlm21; set22%qlm=>qlm22
          ELSE
             i=SET2%N1
             set21%Crd=>SET2%Crd(:,1:i); set22%Crd=>SET2%Crd(:,i+1:)
             set21%Pot=>SET2%Pot(1:i); set22%Pot=>SET2%Pot(i+1:)
             set21%KEY=>SET2%KEY(1:i); set22%KEY=>SET2%KEY(i+1:)
          END IF

          !  Recursion (interaction "SET1" - "SET2") 
          CALL TLRB0_NOGRAD(WS,L,set11,set21) 
          CALL TLRB0_NOGRAD(WS,L,set12,set22)
          CALL TLRB0_NOGRAD(WS,L,set11,set22)
          CALL TLRB0_NOGRAD(WS,L,set12,set21)
       END IF

    END SELECT

  END SUBROUTINE TLRB0_NOGRAD

  !*******************************************************************************

  SUBROUTINE DirectCoulomb1(Crd1,pot1,Grad1)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    !  INPUT:
    REAL(SP), INTENT(IN) :: Crd1(:,:)
    REAL(SP), INTENT(INOUT) :: pot1(:)
    REAL(SP),optional, INTENT(INOUT) :: Grad1(:,:)

    IF(PRESENT(Grad1))THEN
       CALL DirectCoulomb_GRAD(Crd1,pot1,Grad1,Crd1,pot1,Grad1,.true.)
    ELSE
       CALL DirectCoulomb_NOGRAD(Crd1,pot1,Crd1,pot1,.true.)
    ENDIF

  END SUBROUTINE DirectCoulomb1
  !*******************************************************************************

  !*******************************************************************************

  SUBROUTINE DirectCoulomb_GRAD(Crd1,pot1,Grad1,Crd2,pot2,Grad2,SameSet)
    USE DataTypes
    Use CoulombMod
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    !  INPUT:

    REAL(SP), INTENT(IN) :: Crd1(:,:),Crd2(:,:)
    REAL(SP), INTENT(INOUT) :: pot1(:), pot2(:)
    REAL(SP), INTENT(INOUT) :: Grad1(:,:), Grad2(:,:)
    LOGICAL(LGD),INTENT(IN) :: SameSet
    LOGICAL(LGD) :: PointCharge1, PointCharge2

    IF ( SIZE(Crd1,1) < 5 ) THEN
       PointCharge1 = .TRUE.
    ELSE
       PointCharge1 = .FALSE.
    END IF
    IF ( SIZE(Crd2,1) < 5 ) THEN
       PointCharge2 = .TRUE.
    ELSE
       PointCharge2 = .FALSE.
    END IF

    IF ( PointCharge1 .AND. PointCharge2 ) THEN
       CALL pointS_pointS(Crd1(2:4,:),Crd1(1,:),           pot1,Grad1,&
            Crd2(2:4,:),Crd2(1,:),           pot2,Grad2,sameset)
    ELSE IF ( PointCharge1 ) THEN
       CALL pointS_gausS(Crd1(2:4,:),Crd1(1,:),          pot1,Grad1,&
            Crd2(2:4,:),Crd2(1,:),Crd2(5,:),pot2,Grad2,sameset)
    ELSE IF ( PointCharge2 ) THEN
       CALL gausS_pointS(Crd1(2:4,:),Crd1(1,:),Crd1(5,:),pot1,Grad1,&
            Crd2(2:4,:),Crd2(1,:),          pot2,Grad2,sameset)
    ELSE
       CALL gausS_gausS(Crd1(2:4,:),Crd1(1,:),Crd1(5,:),pot1,Grad1,&
            Crd2(2:4,:),Crd2(1,:),Crd2(5,:),pot2,grad2,sameset)
    END IF
 
  END SUBROUTINE DirectCoulomb_GRAD

  SUBROUTINE DirectCoulomb_NOGRAD(Crd1,pot1,Crd2,pot2,SameSet)
    USE DataTypes
    Use CoulombMod
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    !  INPUT:
    REAL(SP), INTENT(IN) :: Crd1(:,:),Crd2(:,:)
    REAL(SP), INTENT(INOUT) :: pot1(:), pot2(:)
    LOGICAL(LGD),INTENT(IN) :: SameSet
    LOGICAL(LGD) :: PointCharge1, PointCharge2

    IF ( SIZE(Crd1,1) < 5 ) THEN
       PointCharge1 = .TRUE.
    ELSE
       PointCharge1 = .FALSE.
    END IF
    IF ( SIZE(Crd2,1) < 5 ) THEN
       PointCharge2 = .TRUE.
    ELSE
       PointCharge2 = .FALSE.
    END IF

    IF ( PointCharge1 .AND. PointCharge2 ) THEN
       CALL pointS_pointS(Crd1(2:4,:),Crd1(1,:),pot1,Crd2(2:4,:),Crd2(1,:),pot2,sameset)
    ELSE IF ( PointCharge1 ) THEN
       CALL pointS_gausS(Crd1(2:4,:),Crd1(1,:),pot1,Crd2(2:4,:),Crd2(1,:),Crd2(5,:),pot2,sameset)
    ELSE IF ( PointCharge2 ) THEN
       CALL gausS_pointS(Crd1(2:4,:),Crd1(1,:),Crd1(5,:),pot1,Crd2(2:4,:),Crd2(1,:),pot2,sameset)
    ELSE
       CALL gausS_gausS(Crd1(2:4,:),Crd1(1,:),Crd1(5,:),pot1,Crd2(2:4,:),Crd2(1,:),Crd2(5,:),pot2,sameset)
    END IF

  END SUBROUTINE DirectCoulomb_NOGRAD


  !*******************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Compute the eigenvector with lowest eigenvalue of a 3 X 3 symmetric matrix  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION LEIGRB(C) RESULT(EV)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    !  INPUT:
    REAL(SP), INTENT(IN) :: C(3,3)   !  The symmetric matrix
    !------------------------------------------------------------------------------
    !  OUTPUT: 
    REAL(SP) :: EV(3)                !  The eigenvector
    !------------------------------------------------------------------------------
    REAL(SP), PARAMETER :: PI=3.14159265358979323846_SP
    REAL(SP) :: a,b,q,s,t,r

    !  Get the lowest root of the characteristic polynomial (a cubic equation)
    a=-(C(1,1)+C(2,2)+C(3,3))
    b=C(1,1)*C(2,2)-C(2,1)**2+C(1,1)*C(3,3)-C(3,1)**2+C(2,2)*C(3,3)-C(3,2)**2
    q=ABS(a*a-3*b)/9
    s=2*SQRT(q)
    t=-C(1,1)*C(2,2)*C(3,3)-C(2,1)*C(3,2)*C(3,1)- &
         C(2,1)*C(3,2)*C(3,1)+C(2,2)*C(3,1)**2+ &
         C(3,3)*C(2,1)**2+C(1,1)*C(3,2)**2
    t=(2*a**3-9*a*b+27*t)/(27*q*s); IF (t>1) t=1; IF (t<-1) t=-1
    t=aCOS(t)
    r=-MAX(s*COS(t/3)+a/3,s*COS((t+2*PI)/3)+a/3,s*COS((t-2*PI)/3)+a/3)
    !  Get the corresponding eigenvector
    EV(1)=(C(2,2)-r)*(C(3,3)-r)-C(3,2)**2
    EV(2)=(-C(2,1)*(C(3,3)-r)+C(3,1)*C(3,2))
    EV(3)=(-C(3,1)*(C(2,2)-r)+C(2,1)*C(3,2))
  END FUNCTION LEIGRB

  !*******************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Compute the splitting information of a set of particles: it is assumed that !
  ! all the charges are unity, THEN the total charge, the center of charge and   !
  ! the moment of inertia matrix of this positive distribution are computed      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE INFldB(SET)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    !  INPUT-OUTPUT:
    TYPE(PARTICLE), INTENT(INOUT) :: SET  !  Set of particles
    !------------------------------------------------------------------------------
    INTEGER(I4B) :: i,np
    REAL(SP) :: t0,tx,ty,tz,px,py,pz,ixx,iyy,izz,iyx,izx,izy

    ! Initialization
    SET%COFC=0.0_SP; SET%MOFI=0.0_SP; np=SIZE(SET%Crd,2)
    px=0.0_SP; py=0.0_SP; pz=0.0_SP
    ixx=0.0_SP; iyy=0.0_SP; izz=0.0_SP; iyx=0.0_SP; izx=0.0_SP; izy=0.0_SP
    ! Summation over all the particles
    DO i=1,np
       tx=SET%Crd(2,i); ty=SET%Crd(3,i); tz=SET%Crd(4,i)
       px=px+tx; py=py+ty; pz=pz+tz
       ixx=ixx+tx*tx; iyy=iyy+ty*ty; izz=izz+tz*tz
       iyx=iyx-ty*tx; izx=izx-tz*tx; izy=izy-tz*ty
    END DO
    ! Final RESULTs
    t0=1.0_SP/np; tx=px*t0; ty=py*t0; tz=pz*t0
    SET%COFC(1)=tx; SET%COFC(2)=ty; SET%COFC(3)=tz
    SET%MOFI(2,1)=iyx+px*ty
    SET%MOFI(3,1)=izx+px*tz
    SET%MOFI(3,2)=izy+py*tz
    SET%MOFI(1,1)=iyy+izz-py*ty-pz*tz
    SET%MOFI(2,2)=ixx+izz-px*tx-pz*tz
    SET%MOFI(3,3)=ixx+iyy-px*tx-py*ty
    SET%MOFI(1,2)=SET%MOFI(2,1)
    SET%MOFI(1,3)=SET%MOFI(3,1)
    SET%MOFI(2,3)=SET%MOFI(3,2)
    SET%GETSPL=.FALSE.
  END SUBROUTINE INFldB

  !*******************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Split a "parent" set of particles (including their potential and fields)    !
  ! into two "children" sets by using a splitting plane passing through the      !
  ! center of the "parent" set, and oriented perpENDicularly to one of the three !
  ! principal axis of the "parent" set (the one with the smallest moment of      !
  ! inertia)                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE PSPLRB(SET,SET1,SET2)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    !  INPUT-OUTPUT:
    TYPE(PARTICLE), INTENT(INOUT) :: SET ! Parent set of particles
    !------------------------------------------------------------------------------
    !  OUTPUT:
    TYPE(PARTICLE), INTENT(OUT) :: SET1  ! Firt child set of particles
    TYPE(PARTICLE), INTENT(OUT) :: SET2  ! Second child set of particles
    !------------------------------------------------------------------------------
    INTEGER(I4B) :: i,j,i1,np,np1
    REAL(SP) :: t0,t1,v1,v2,v3,xi,yi,zi,xj,yj,zj
    INTEGER(I4B), POINTER :: key(:)
    REAL(SP), POINTER :: Crd(:,:),Grad(:,:),Pot(:)
    REAL(SP) :: CrdTmp(SIZE(set%Crd,1)),GradTmp(3),potTmp

    KEY=>SET%KEY; Crd=>SET%Crd; Grad=>SET%Gradn; Pot=>SET%Pot
    v1=SET%SPL(1); v2=SET%SPL(2); v3=SET%SPL(3)
    t0=SET%COFC(1)*v1+SET%COFC(2)*v2+SET%COFC(3)*v3
    np=SIZE(Crd,2); i=0; j=np+1; np1=0
    !  Look for a particle that should be in the second "child"
    DO
       DO
          i=i+1; xi=Crd(2,i); yi=Crd(3,i); zi=Crd(4,i)
          t1=xi*v1+yi*v2+zi*v3
          IF (t1<=t0.OR.i==np) EXIT
       END DO
       !  Look for a particle that should be in the first "child"
       DO
          j=j-1; xj=Crd(2,j); yj=Crd(3,j); zj=Crd(4,j)
          t1=xj*v1+yj*v2+zj*v3
          IF (t1>t0) np1=j
          IF (t1>t0.OR.j==1) EXIT
       END DO
       !  Exchange the "i" and "j" particles (potential and fields too) and repeat
       IF (i<j) THEN
          CrdTmp = Crd(:,i)
          Crd(:,i) = Crd(:,j); Crd(:,j) = CrdTmp
          potTmp = Pot(i)
          Pot(i) = Pot(j); Pot(j) = potTmp
          i1 = KEY(i)
          KEY(i)   = KEY(j);   KEY(j)   = i1

          IF(associated(Grad))THEN
             GradTmp = Grad(:,i)
             Grad(:,i) = Grad(:,j); Grad(:,j) = GradTmp
          ENDIF
       ELSE
          EXIT
       END IF
    END DO
    ! Save the information
    SET1%Crd=>Crd(:,1:np1); SET2%Crd=>Crd(:,np1+1:)
    SET1%Pot=>Pot(1:np1); SET2%Pot=>Pot(np1+1:)
    SET1%KEY=>KEY(1:np1); SET2%KEY=>KEY(np1+1:)
    SET%SPLIT=.FALSE.; SET%N1=np1
    IF(associated(Grad))THEN
       SET1%Gradn=>Grad(:,1:np1); SET2%Gradn=>Grad(:,np1+1:)
    ENDIF
  END SUBROUTINE PSPLRB

  !*******************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  "Decision maker" SUBROUTINE. It decides IF two sets of particles should     c
  ! interact by direct summation, should interact by multipole expansion, or     c
  ! should be split                                                              c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION DMAKRB(WS,L,SET1,SET2) RESULT(IDEC)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    !  INPUT:
    REAL(SP), INTENT(IN) :: WS         ! Well-separatedness PARAMETER
    INTEGER(I4B), INTENT(IN) ::            L      ! Order of the multipole expansion
    TYPE(PARTICLE), INTENT(INOUT) :: SET1   ! First set of particles
    TYPE(PARTICLE), INTENT(INOUT) :: SET2   ! Second set of particles
    !------------------------------------------------------------------------------
    !  OUTPUT:
    INTEGER(I4B) :: IDEC                   ! THE decision:
    !                                     -- "IDEC=0" for direct summation
    !                                     -- "IDEC=1" for multipole expansion
    !                                     -- "IDEC=2" for splitting
    !------------------------------------------------------------------------------
    INTEGER(I4B) :: np1,np2,mnpart
    REAL(SP) :: r1,r2,r12

    np1=SIZE(SET1%Crd,2); np2=SIZE(SET2%Crd,2)
    mnpart=PCMNP*(L+2)**2  ! Get the minimum number of particles
    IF (ASSOCIATED(SET1%Crd,SET2%Crd)) THEN
       !  Interaction of "SET1" with itself (direct sum or splitting)
       IF (np1<=4*mnpart) THEN
          IDEC=0
       ELSE
          IDEC=2
!!! Jana
          SET1%SPLIT=.TRUE.
!!!
          IF (SET1%GETSPL) CALL INFldB(SET1)
       END IF
    ELSE 
       !  Interaction "SET1" - "SET2" (direct sum, multipole expansion or splitting) 
       IF (SET1%GETSPL) CALL INFldB(SET1)
       IF (SET2%GETSPL) CALL INFldB(SET2)
       !  SMALL sets of particles, direct summation
       IDEC=0
       IF (np1>mnpart.AND.np2>mnpart) THEN
          !  SPLIT (the children should have more particles than "mnpart") ...
          IF (MAX(np1,np2).gt.2*mnpart) THEN
             IDEC=2
!!! Jana
             SET1%SPLIT=.TRUE.
             SET2%SPLIT=.TRUE.
!!!
          END IF

          ! ... or interact via MULTIPOLES (well-separated sets)
          r12=(SET1%COFC(1)-SET2%COFC(1))**2 + &
               (SET1%COFC(2)-SET2%COFC(2))**2 + &
               (SET1%COFC(3)-SET2%COFC(3))**2
          r1=(SET1%MOFI(1,1)+SET1%MOFI(2,2)+SET1%MOFI(3,3))/(2*np1)
          r2=(SET2%MOFI(1,1)+SET2%MOFI(2,2)+SET2%MOFI(3,3))/(2*np2)
          IF (r12>=2*(WS+1)**2*(r1+r2)) IDEC=1
       END IF
    END IF
  END FUNCTION DMAKRB

  !*******************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Restore the list of particles and fields to their original ordering. A      !
  ! RECURSIVE implementation of QUICKSORT is used.                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE RORDRB(IP0,IPL,KEY,Crd,pot,Grad)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    !  INPUT:
    INTEGER(I4B), INTENT(IN) ::           IP0       ! First particle
    INTEGER(I4B), INTENT(IN) ::           IPL       ! Last particle
    !------------------------------------------------------------------------------
    !  INPUT-OUTPUT:
    INTEGER(I4B), INTENT(INOUT) :: KEY(:)           ! Key (for reordering)
    REAL(SP), INTENT(INOUT) :: Crd(:,:)   ! List of particles
    REAL(SP), INTENT(INOUT) :: Pot(:)     ! List of potential 
    REAL(SP), optional, INTENT(INOUT) :: Grad(:,:)   ! List of fields
    !------------------------------------------------------------------------------
    INTEGER(I4B) :: i,j,i1,np,np0
    REAL(SP) :: CrdTmp(SIZE(Crd,1)),GradTmp(3), PotTmp

    np=IPL+1-IP0
    ! For 1 particle, DO nothing
    IF (np==2) THEN
       ! 2 particles CASE ...
       IF (KEY(IP0)>KEY(IPL)) THEN
          !     ... exchange particles "IP0" and "IPL"
          i1 = KEY(IP0)
          KEY(IP0)   = KEY(IPL);   KEY(IPL)   = i1

          CrdTmp = Crd(:,IP0)
          Crd(:,IP0) = Crd(:,IPL); Crd(:,IPL) = CrdTmp
          PotTmp = Pot(IP0)
          Pot(IP0) = Pot(IPL); Pot(IPL) = PotTmp
          IF(Present(GRAD))THEN
             GradTmp = Grad(:,IP0)
             Grad(:,IP0) = Grad(:,IPL); Grad(:,IPL) = GradTmp
          ENDIF
      ENDif
    ELSEif (np>2) THEN
       ! More than 2 particles CASE ...
       np0=(IP0+IPL)/2  ! Middle of the set
       i=IP0-1; j=IPL+1
       !     ... look for a particle that should be in the second half
       DO
          DO
             i=i+1
             IF (KEY(i)>np0.OR.i==IPL) EXIT
          END DO
          !     ... look for a particle that should be in the first half
          DO
             j=j-1
             IF (KEY(j)<=np0.OR.j==IP0) EXIT
          END DO
          !     ... exchange the "i" and "j" particles and repeat
          IF (i<j) THEN
             i1 = KEY(IP0); KEY(IP0) = KEY(IPL); KEY(IPL) = i1
             CrdTmp = Crd(:,IP0); Crd(:,IP0) = Crd(:,IPL); Crd(:,IPL) = CrdTmp
             PotTmp = Pot(IP0); Pot(IP0) = Pot(IPL); Pot(IPL) = PotTmp
             IF(Present(GRAD))THEN
                GradTmp = Grad(:,IP0)
                Grad(:,IP0) = Grad(:,IPL); Grad(:,IPL) = GradTmp
             ENDIF
          ELSE
             EXIT
          END IF
       END DO
       !     ... recursion
       CALL RORDRB(IP0,np0,key,Crd,Pot,Grad)
       CALL RORDRB(np0+1,IPL,key,Crd,Pot,Grad)
    ENDif
  END SUBROUTINE RORDRB

  !*******************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Computes the multipole expansion of a set of particles and its derivatives  !
  ! with respect to X, Y and Z                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GMEXRB_NOGRAD(L,SET)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    !  INPUT:
    INTEGER(I4B), INTENT(IN) ::            L    ! Order of the multipole expansion
    !------------------------------------------------------------------------------
    !  INPUT-OUTPUT:
    TYPE(PARTICLE), INTENT(INOUT) :: SET  ! Set of particles
    !------------------------------------------------------------------------------
    INTEGER(I4B) :: i,m1,l1,mm1,mp1,lm1,lmm,np
    REAL(SP) :: x(1:3),r2,plm,plm0,a,b,c,s,o1(1:3),i1,i2
    REAL(SP), POINTER :: Crd(:,:),qlm(:,:,:)
    REAL(SP) :: qlm0(0:L,0:L)

    lm1=L-1; Crd=>SET%Crd; np=SIZE(Crd,2)
    o1(1:3)=SET%COFC(1:3)
    qlm0=0
    DO i=1,np
       !  Some initialization ...
       x(1:3)=o1(1:3)-Crd(2:4,i)
       r2=DOT_PRODUCT(x,x)
       !  Contributions to "q_{l,0}"
       c=Crd(1,i); plm=c; plm0=0; i1=1.0_SP; i2=1.0_SP
       DO l1=0,lm1
          qlm0(l1,0)=qlm0(l1,0)+plm
          b=(x(3)*i1*plm-r2*plm0)/i2; plm0=plm; plm=b
          i1=i1+2.0_SP; i2=i2+i1
       END DO
       qlm0(l1,0)=qlm0(l1,0)+plm
       !  Rest of contributions ...
       c=c/2; s=c*x(2); c=c*x(1)
       DO m1=1,L
          plm=1; plm0=0; i1=m1+m1+1; i2=i1
          DO l1=m1,lm1
             ! ... add contributions to "q_{l1,m1}" and "q_{l1,-m1}"
             lmm=l1-m1
             qlm0(l1,m1)=qlm0(l1,m1)+plm*c
             qlm0(lmm,l1)=qlm0(lmm,l1)+plm*s
             b=(x(3)*i1*plm-r2*plm0)/i2; plm0=plm; plm=b
             i1=i1+2.0_SP; i2=i2+i1
          END DO
          lmm=l1-m1
          qlm0(l1,m1)=qlm0(l1,m1)+plm*c
          qlm0(lmm,l1)=qlm0(lmm,l1)+plm*s
          ! ... compute "N_{m+1,m+1}"
          b=c; a=1.0_SP/(2*(m1+1))
          c=(b*x(1)-s*x(2))*a; s=(s*x(1)+b*x(2))*a
       END DO
    END DO
    ! Copy to the output matrix
    qlm=>SET%QLM
    qlm=0; qlm(1,:L,:L)=qlm0(:L,:L)
    ! Get the coefficients in final form and the derivatives (minus)
    a=qlm(1,0,0)
    DO l1=1,L
       lmm=l1-1
       a=-qlm(1,l1,1); qlm(1,l1,1)=-2*a
       b=-qlm(1,lmm,l1); qlm(1,lmm,l1)=-2*b
    END DO
    DO m1=2,L
       mp1=m1+1; mm1=m1-1
       DO l1=m1,L
          lmm=l1-m1
          a=-qlm(1,l1,m1); qlm(1,l1,m1)=-2*a
          b=-qlm(1,lmm,l1); qlm(1,lmm,l1)=-2*b
       END DO
    END DO
    SET%GETQLM=.FALSE.
  END SUBROUTINE GMEXRB_NOGRAD

  !*******************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Evaluate the interaction between a set of particles and several multipole   !
  ! expansions (ME) with the same order and origen, adding the contribution to   !
  ! potential and fields                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE EMEXRB_NOGRAD(L,QLM,ORIG,SET)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    !  INPUT:
    INTEGER(I4B), INTENT(IN) :: L                          ! Order of the ME
    REAL(SP), INTENT(IN) :: QLM(1,0:L,0:L)  ! Coefficients of the ME
    REAL(SP), INTENT(IN) :: ORIG(3)             ! Origen of the ME
    !------------------------------------------------------------------------------
    !  INPUT-OUTPUT:
    TYPE(PARTICLE), INTENT(INOUT) :: SET             ! Set of particles
    !------------------------------------------------------------------------------
    INTEGER(I4B) :: i,m1,l1,lmm,np,lp1
    REAL(SP) :: x(1:3),rm2,plm,plm0,c0,s0,c,s,p, &
         o1(1:3),i1,i2,t1
    REAL(SP), POINTER :: Crd(:,:),Pot(:)

    Crd=>SET%Crd; Pot=>SET%Pot
    o1(1:3)=ORIG(1:3)
    np=SIZE(Crd,2); lp1=L+1
    DO i=1,np
       ! Some initialization ...
       x(1:3)=o1(1:3)-Crd(2:4,i)
       rm2=1/DOT_PRODUCT(x,x)
       ! Contributions Fldom "q_{l,0}"
       !      c=Crd(1,i)*SQRT(rm2)
       c=SQRT(rm2) !to give potential and not energy -DMY
       p=0.0_SP
       plm=c; plm0=0; i1=1.0_SP; i2=0.0_SP
       DO l1=0,L
          p=p+QLM(1,l1,0)*plm
          c0=(x(3)*i1*plm-i2*plm0)*rm2; plm0=plm; plm=c0
          i2=i2+i1; i1=i1+2.0_SP
       END DO
       ! Rest of contributions ...
       c=c*rm2; s=c*x(2); c=c*x(1)
       DO m1=1,lp1
          plm=1; plm0=0; i1=m1+m1+1; i2=0.0_SP
          DO l1=m1,L
             ! ... add contributions Fldom "M_{l1,m1}"
             lmm=l1-m1; c0=c*plm; s0=s*plm
             p=p+QLM(1,l1,m1)*c0+QLM(1,lmm,l1)*s0
             ! ... compute "M_{l,m}"
             c0=(x(3)*i1*plm-i2*plm0)*rm2; plm0=plm; plm=c0
             i2=i2+i1; i1=i1+2.0_SP
          END DO
          lmm=l1-m1; c0=c*plm; s0=s*plm
          ! Compute "M_{m+1,m+1}"
          c0=c; t1=(m1+m1+1)*rm2
          c=(c0*x(1)-s*x(2))*t1; s=(s*x(1)+c0*x(2))*t1
       END DO
       Pot(i)=Pot(i)+p
    END DO
  END SUBROUTINE EMEXRB_NOGRAD

  RECURSIVE SUBROUTINE TLRB0_GRAD(WS,L,SET1,SET2)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    !  INPUT:
    REAL(SP), INTENT(IN) ::     WS    ! Well-separatedness PARAMETER
    INTEGER(I4B), INTENT(IN) ::            L     ! Order of the multipole expansion
    !------------------------------------------------------------------------------
    !  INPUT-OUTPUT:
    TYPE(PARTICLE), INTENT(INOUT) :: SET1  ! First set of particles
    TYPE(PARTICLE), INTENT(INOUT) :: SET2  ! Second set of particles
    !------------------------------------------------------------------------------
    INTEGER(I4B) :: i,n1,n2,idec
    REAL(SP), TARGET :: qlm11(4,0:L+1,0:L+1),qlm12(4,0:L+1,0:L+1), &
         qlm21(4,0:L+1,0:L+1),qlm22(4,0:L+1,0:L+1)
    REAL(SP), POINTER :: Crd1(:,:),Pot1(:),Grad1(:,:),Crd2(:,:),Pot2(:), Grad2(:,:)
    LOGICAL(LGD) :: itself
    TYPE(PARTICLE) :: set11,set12,set21,set22

    !  Initialization
    Crd1=>SET1%Crd; Pot1=>SET1%Pot; Grad1=>SET1%Gradn
    Crd2=>SET2%Crd; Pot2=>SET2%Pot; Grad2=>SET2%Gradn
    n1=SIZE(Crd1,2); n2=SIZE(Crd2,2)
    set11%getqlm=.TRUE.; set11%getspl=.TRUE.; set11%split=.TRUE.
    set12%getqlm=.TRUE.; set12%getspl=.TRUE.; set12%split=.TRUE.
    set21%getqlm=.TRUE.; set21%getspl=.TRUE.; set21%split=.TRUE.
    set22%getqlm=.TRUE.; set22%getspl=.TRUE.; set22%split=.TRUE.
    set11%qlm=>qlm11; set12%qlm=>qlm12; set21%qlm=>qlm21; set22%qlm=>qlm22
    itself=ASSOCIATED(Crd1,Crd2)

    !...............................................................................
    !  Get decision (direct sum, multipole expansion, or splitting)
    idec=DMAKRB(WS,L,SET1,SET2)
    !  Execute the decision
    SELECT CASE (idec)

       !  DIRECT SUMMATION
    CASE (0)
       CALL DirectCoulomb_GRAD(Crd1,pot1,Grad1,Crd2,pot2,Grad2,itself)
       !  MULTIPOLE EXPANSION
    CASE (1)
       IF (SET1%GETQLM) CALL GMEXRB_GRAD(L,SET1)
       CALL EMEXRB_GRAD(L,SET1%QLM,SET1%COFC,SET2)
       IF (SET2%GETQLM) CALL GMEXRB_GRAD(L,SET2)
       CALL EMEXRB_GRAD(L,SET2%QLM,SET2%COFC,SET1)

       !  SPLITTING
    CASE (2)
       !  Split the first distribution
       IF (SET1%SPLIT) THEN
          SET1%SPL=LEIGRB(SET1%MOFI)
          CALL PSPLRB(SET1,set11,set12)
          set11%qlm=>qlm11; set12%qlm=>qlm12

       ELSE
          i=SET1%N1
          set11%Crd=>SET1%Crd(:,1:i); set12%Crd=>SET1%Crd(:,i+1:)
          set11%Pot=>SET1%Pot(1:i);   set12%Pot=>SET1%Pot(i+1:)
          set11%Gradn=>SET1%Gradn(:,1:i); set12%Gradn=>SET1%Gradn(:,i+1:)
          set11%KEY=>SET1%KEY(1:i);   set12%KEY=>SET1%KEY(i+1:)
       END IF

       IF (itself) THEN
          !  Recursion (interaction of "SET1" with itself)
          CALL TLRB0_GRAD(WS,L,set11,set11)
          CALL TLRB0_GRAD(WS,L,set12,set12)
          CALL TLRB0_GRAD(WS,L,set11,set12)
       ELSE
          !  Split the second distribution
          IF (SET2%SPLIT) THEN
             SET2%SPL=LEIGRB(SET2%MOFI)
             CALL PSPLRB(SET2,set21,set22)
             set21%qlm=>qlm21; set22%qlm=>qlm22
          ELSE
             i=SET2%N1
             set21%Crd=>SET2%Crd(:,1:i); set22%Crd=>SET2%Crd(:,i+1:)
             set21%Pot=>SET2%Pot(1:i);   set22%Pot=>SET2%Pot(i+1:)
             set21%Gradn=>SET2%Gradn(:,1:i); set22%Gradn=>SET2%Gradn(:,i+1:)
             set21%KEY=>SET2%KEY(1:i);   set22%KEY=>SET2%KEY(i+1:)
          END IF
          !  Recursion (interaction "SET1" - "SET2")
          CALL TLRB0_GRAD(WS,L,set11,set21)
          CALL TLRB0_GRAD(WS,L,set12,set22)
          CALL TLRB0_GRAD(WS,L,set11,set22)
          CALL TLRB0_GRAD(WS,L,set12,set21)
       END IF
    END SELECT
  END SUBROUTINE TLRB0_GRAD

  !*******************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Computes the multipole expansion of a set of particles and its derivatives  !
  ! with respect to X, Y and Z                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GMEXRB_GRAD(L,SET)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    !  INPUT:
    INTEGER(I4B), INTENT(IN) ::            L    ! Order of the multipole expansion
    !------------------------------------------------------------------------------
    !  INPUT-OUTPUT:
    TYPE(PARTICLE), INTENT(INOUT) :: SET  ! Set of particles
    !------------------------------------------------------------------------------
    INTEGER(I4B) :: i,m1,l1,mm1,mp1,lm1,lp1,lmm,lmm2,np
    REAL(SP) :: x(1:3),r2,plm,plm0,a,b,c,s,o1(1:3),i1,i2
    REAL(SP), POINTER :: Crd(:,:),qlm(:,:,:)
    REAL(SP) :: qlm0(0:L,0:L)

    lm1=L-1; Crd=>SET%Crd; np=SIZE(Crd,2)
    o1(1:3)=SET%COFC(1:3)
    qlm0=0
    DO i=1,np
       !  Some initialization ...
       x(1:3)=o1(1:3)-Crd(2:4,i)
       r2=DOT_PRODUCT(x,x)
       !  Contributions to "q_{l,0}"
       c=Crd(1,i); plm=c; plm0=0; i1=1.0_SP; i2=1.0_SP
       DO l1=0,lm1
          qlm0(l1,0)=qlm0(l1,0)+plm
          b=(x(3)*i1*plm-r2*plm0)/i2; plm0=plm; plm=b
          i1=i1+2.0_SP; i2=i2+i1
       END DO
       qlm0(l1,0)=qlm0(l1,0)+plm
       !  Rest of contributions ...
       c=c/2; s=c*x(2); c=c*x(1)
       DO m1=1,L
          plm=1; plm0=0; i1=m1+m1+1; i2=i1
          DO l1=m1,lm1
             ! ... add contributions to "q_{l1,m1}" and "q_{l1,-m1}"
             lmm=l1-m1
             qlm0(l1,m1)=qlm0(l1,m1)+plm*c
             qlm0(lmm,l1)=qlm0(lmm,l1)+plm*s
             b=(x(3)*i1*plm-r2*plm0)/i2; plm0=plm; plm=b
             i1=i1+2.0_SP; i2=i2+i1
          END DO
          lmm=l1-m1
          qlm0(l1,m1)=qlm0(l1,m1)+plm*c
          qlm0(lmm,l1)=qlm0(lmm,l1)+plm*s
          ! ... compute "N_{m+1,m+1}"
          b=c; a=1.0_SP/(2*(m1+1))
          c=(b*x(1)-s*x(2))*a; s=(s*x(1)+b*x(2))*a
       END DO
    END DO
    ! Copy to the output matrix
    qlm=>SET%QLM
    qlm=0; qlm(1,:L,:L)=qlm0(:L,:L)
    ! Get the coefficients in final form and the derivatives (minus)
    a=qlm(1,0,0); qlm(4,1,0)=a; qlm(2,1,1)=-a; qlm(3,0,1)=-a
    DO l1=1,L
       lp1=l1+1; lmm=l1-1
       a=qlm(1,l1,0); qlm(4,lp1,0)=a; qlm(2,lp1,1)=-a; qlm(3,l1,lp1)=-a
       a=-qlm(1,l1,1); qlm(1,l1,1)=-2*a; qlm(4,lp1,1)=qlm(1,l1,1)
       b=-qlm(1,lmm,l1); qlm(1,lmm,l1)=-2*b; qlm(4,l1,lp1)=qlm(1,lmm,l1)
       qlm(3,lmm,lp1)=qlm(3,lmm,lp1)+a; qlm(2,lmm,lp1)=qlm(2,lmm,lp1)+b
       qlm(2,lp1,2)=qlm(2,lp1,2)+a; qlm(3,lp1,2)=qlm(3,lp1,2)-b
       qlm(2,lp1,0)=qlm(2,lp1,0)-a; qlm(3,lp1,0)=qlm(3,lp1,0)-b
    END DO
    DO m1=2,L
       mp1=m1+1; mm1=m1-1
       DO l1=m1,L
          lp1=l1+1; lmm=l1-m1; lmm2=l1+2-m1
          a=-qlm(1,l1,m1); qlm(1,l1,m1)=-2*a; qlm(4,lp1,m1)=qlm(1,l1,m1)
          b=-qlm(1,lmm,l1); qlm(1,lmm,l1)=-2*b; qlm(4,lp1-m1,lp1)=qlm(1,lmm,l1)
          qlm(3,lmm2,lp1)=qlm(3,lmm2,lp1)+a
          qlm(2,lmm2,lp1)=qlm(2,lmm2,lp1)-b
          qlm(3,lmm,lp1)=qlm(3,lmm,lp1)+a; qlm(2,lmm,lp1)=qlm(2,lmm,lp1)+b
          qlm(2,lp1,mp1)=qlm(2,lp1,mp1)+a; qlm(3,lp1,mp1)=qlm(3,lp1,mp1)-b
          qlm(2,lp1,mm1)=qlm(2,lp1,mm1)-a; qlm(3,lp1,mm1)=qlm(3,lp1,mm1)-b
       END DO
    END DO
    SET%GETQLM=.FALSE.
  END SUBROUTINE GMEXRB_GRAD

  !*******************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Evaluate the interaction between a set of particles and several multipole   !
  ! expansions (ME) with the same order and origen, adding the contribution to   !
  ! potential and fields                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE EMEXRB_GRAD(L,QLM,ORIG,SET)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    !  INPUT:
    INTEGER(I4B), INTENT(IN) :: L                          ! Order of the ME
    REAL(SP), INTENT(IN) :: QLM(4,0:L+1,0:L+1)  ! Coefficients of the ME
    REAL(SP), INTENT(IN) :: ORIG(3)             ! Origen of the ME
    !------------------------------------------------------------------------------
    !  INPUT-OUTPUT:
    TYPE(PARTICLE), INTENT(INOUT) :: SET             ! Set of particles
    !------------------------------------------------------------------------------
    INTEGER(I4B) :: i,m1,l1,lmm,np,lp1
    REAL(SP) :: x(1:3),rm2,plm,plm0,c0,s0,c,s,p,fx(1:3), &
         o1(1:3),i1,i2,t1
    REAL(SP), POINTER :: Crd(:,:),Pot(:),Grad(:,:)

    Crd=>SET%Crd; Pot=>SET%Pot; Grad=>SET%Gradn
    o1(1:3)=ORIG(1:3)
    np=SIZE(Crd,2); lp1=L+1
    DO i=1,np
       ! Some initialization ...
       x(1:3)=o1(1:3)-Crd(2:4,i)
       rm2=1/DOT_PRODUCT(x,x)
       ! Contributions Fldom "q_{l,0}"
       !      c=Crd(1,i)*SQRT(rm2)
       c=SQRT(rm2) !to give potential and not energy -DMY
       p=0.0_SP; fx(1:3)=0.0_SP
       plm=c; plm0=0; i1=1.0_SP; i2=0.0_SP
       DO l1=0,L
          p=p+QLM(1,l1,0)*plm
          fx(1:3)=fx(1:3)+QLM(2:4,l1,0)*plm
          c0=(x(3)*i1*plm-i2*plm0)*rm2; plm0=plm; plm=c0
          i2=i2+i1; i1=i1+2.0_SP
       END DO
       fx(1:3)=fx(1:3)+QLM(2:4,l1,0)*plm
       ! Rest of contributions ...
       c=c*rm2; s=c*x(2); c=c*x(1)
       DO m1=1,lp1
          plm=1; plm0=0; i1=m1+m1+1; i2=0.0_SP
          DO l1=m1,L
             ! ... add contributions Fldom "M_{l1,m1}"
             lmm=l1-m1; c0=c*plm; s0=s*plm
             p=p+QLM(1,l1,m1)*c0+QLM(1,lmm,l1)*s0
             fx(1:3)=fx(1:3)+QLM(2:4,l1,m1)*c0+QLM(2:4,lmm,l1)*s0
             ! ... compute "M_{l,m}"
             c0=(x(3)*i1*plm-i2*plm0)*rm2; plm0=plm; plm=c0
             i2=i2+i1; i1=i1+2.0_SP
          END DO
          lmm=l1-m1; c0=c*plm; s0=s*plm
          fx(1:3)=fx(1:3)+QLM(2:4,l1,m1)*c0+QLM(2:4,lmm,l1)*s0
          ! Compute "M_{m+1,m+1}"
          c0=c; t1=(m1+m1+1)*rm2
          c=(c0*x(1)-s*x(2))*t1; s=(s*x(1)+c0*x(2))*t1
       END DO
       Pot(i)=Pot(i)+p
       Grad(1:2,i)=Grad(1:2,i)-fx(1:2)
       Grad(3,i)=Grad(3,i)+fx(3)
    END DO
  END SUBROUTINE EMEXRB_GRAD


 subroutine print_particle(set)
   type(particle) :: set

   if (associated(set%crd)) then
      write(6,*)'Crd is ',size(set%crd,1),' by ',size(set%crd,2)
   else
      write(6,*)'Crd is null'
   end if

   if (associated(set%Gradn)) then
      write(6,*)'Grad is ',size(set%gradn,1),' by ',size(set%gradn,2)
   else
      write(6,*)'Grad is null!'
   end if

   if (associated(set%key)) then
      write(6,*)'key is ',size(set%key)
   else
      write(6,*)'key is null!'
   end if

   write(6,*)'cofc = ',set%cofc
   write(6,*)'mofi = ',set%mofi
   write(6,*)'spl = ',set%spl
   if (associated(set%qlm)) then
      write(6,*)'qlm is ',size(set%qlm,1),' by ',size(set%qlm,2),' by ',size(set%qlm,3)
   else
      write(6,*)'qlm is null!'
   end if
   write(6,*)'getqlm = ',set%getqlm
   write(6,*)'getspl = ',set%getspl
   write(6,*)'n1 = ',set%n1

 end subroutine print_particle

END MODULE RBM
