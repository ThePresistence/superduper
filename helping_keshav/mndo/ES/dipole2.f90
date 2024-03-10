

! created by jie liu on 10/16
! see GUGACi/property.f90 for more details
!
SUBROUTiNE CiSPROP(Perm,R,CiE,CiPROP,CiStateSym,SymLabel,icisym,nroots,LM3)
  USE LiMiT, ONLY: LM1,LMX,LMZ,LMACT,LMREF,LMPROP,LMSTAT,LMGRD,MAXGPU
  Use Const, Only: A0, AU2DBY, AU2WAV, EV, OLDCF
  iMPLiCiT REAL*8 (A-H,O-Z)
  COMMON /ATOMC / COORD(3,LM1)
  COMMON /ATOMS / NUMAT,NAT(LM1),NFiRST(LM1),NLAST(LM1)
  COMMON /DPARM1/ UDD(LMZ),ZD(LMZ),BETAD(LMZ)
  COMMON /DNBND / iii(LMZ),iiiD(LMZ)
  COMMON /PAROPT/ USS(LMZ),UPP(LMZ),ZS(LMZ),ZP(LMZ),BETAS(LMZ*3)
  COMMON /PARDER/ CORE(LMZ)
  COMMON /NBFiLE/ NBF(20)
  COMMON /iNOPT2/ iN2(300)
  REAL*8 Perm(LM3*LM3,*),R(LM3*LM3,*),CiE(*),CiPROP(LMPROP,*)
  iNTEGER CiStateSym(*),icisym(*)
  CHARACTER(Len=4) SymLabel(*)
  REAL*8, Dimension(3) :: RefDip,Rmoi,dip,mur,mup,mag
  REAL*8, Dimension(:), Allocatable :: murabs, murphi, murtheta
  REAL*8, Dimension(:), Allocatable :: mupabs, mupphi, muptheta
  REAL*8, Dimension(:), Allocatable :: magabs, magphi, magtheta
  REAL*8, Dimension(:,:), Allocatable :: XYZ,Sao,TMP,TMP2
  REAL*8, Dimension(:,:,:), Allocatable :: Rst, Delst, RxDelst
  REAL*8, Dimension(:,:,:), Allocatable :: Rao,Delao,RxDelao
  iNTEGER s,t
  LOGiCAL ECP 

  iuvcd = iN2(146)
  if(iuvcd.le.0) return

  CALL OLDCF
  iOP = iN2(2)
  nb6 = NBF(6)
  ECP = iOP == -5 .OR. iOP == -6 .OR. iOP == -8 .OR. iOP == -9
  NZ  = SiZE(ZS)
  NAO = LM3
  NATOM = NUMAT

  allocate(XYZ(3,NATOM))
  allocate(Sao(NAO,NAO))
  allocate(TMP(NAO,NAO))
  allocate(TMP2(NAO,NAO))
  allocate(Rao(NAO,NAO,3))
  allocate(Delao(NAO,NAO,3))
  allocate(RxDelao(NAO,NAO,3))
 
  XYZ = COORD / A0
  Call CalcPropint(Sao, Rao, Delao, RxDelao, NAT, NFiRST, NLAST, & 
    XYZ, ZS, ZP, ZD, iii, iii, iiiD, NZ, NATOM, NAO, nb6, ECP)
  
  ! S**(-1/2)
  Call MatDiagSquare(Sao,TMP,NAO,ierror)
  if(ierror.ne.0) then
    write(6,*) "overlap diagonalization fails"
    stop
  endif
  call vecinit(TMP2,NAO*NAO,0.d0)
  do i=1,NAO
    if(TMP(i,1).le.0.d0) then
      write(6,*) "eigenvalue of S is negative"
      stop
    endif
    TMP2(i,i) = 1.d0/sqrt(TMP(i,1))
  enddo
  call matmult(TMP2,Sao,TMP,NAO,NAO,NAO,NAO,NAO,NAO,3)
  call matmult(Sao,TMP,TMP2,NAO,NAO,NAO,NAO,NAO,NAO,1)

  do i=1,3
    call matmult(TMP2,Rao(1,1,i),TMP,NAO,NAO,NAO,NAO,NAO,NAO,1)
    call matmult(TMP,TMP2,Rao(1,1,i),NAO,NAO,NAO,NAO,NAO,NAO,1)
    call matmult(TMP2,Delao(1,1,i),TMP,NAO,NAO,NAO,NAO,NAO,NAO,1)
    call matmult(TMP,TMP2,Delao(1,1,i),NAO,NAO,NAO,NAO,NAO,NAO,1)
    call matmult(TMP2,RxDelao(1,1,i),TMP,NAO,NAO,NAO,NAO,NAO,NAO,1)
    call matmult(TMP,TMP2,RxDelao(1,1,i),NAO,NAO,NAO,NAO,NAO,NAO,1)
  enddo

  mstate = nroots
  nstate = nroots

  allocate(    Rst(mstate,2,3))
  allocate(  Delst(mstate,2,3))
  allocate(RxDelst(mstate,2,3))

      Rst = 0.D0
    Delst = 0.D0
  RxDelst = 0.D0

  do i=1,nroots
     call vecdot(Rst(i,1,1),Perm(1,i),Rao(1,1,1),NAO*NAO)
     call vecdot(Rst(i,1,2),Perm(1,i),Rao(1,1,2),NAO*NAO)
     call vecdot(Rst(i,1,3),Perm(1,i),Rao(1,1,3),NAO*NAO)
  enddo

  do i=1,nroots
    call vecdot(    Rst(i,2,1),R(1,i),    Rao(1,1,1),NAO*NAO) 
    call vecdot(    Rst(i,2,2),R(1,i),    Rao(1,1,2),NAO*NAO) 
    call vecdot(    Rst(i,2,3),R(1,i),    Rao(1,1,3),NAO*NAO) 
    call vecdot(  Delst(i,2,1),R(1,i),  Delao(1,1,1),NAO*NAO) 
    call vecdot(  Delst(i,2,2),R(1,i),  Delao(1,1,2),NAO*NAO) 
    call vecdot(  Delst(i,2,3),R(1,i),  Delao(1,1,3),NAO*NAO) 
    call vecdot(RxDelst(i,2,1),R(1,i),RxDelao(1,1,1),NAO*NAO) 
    call vecdot(RxDelst(i,2,2),R(1,i),RxDelao(1,1,2),NAO*NAO) 
    call vecdot(RxDelst(i,2,3),R(1,i),RxDelao(1,1,3),NAO*NAO) 
  enddo

  ! Contribution of the atom cores:
  RefDip = 0.D0
  do i=1, natom
     RefDip = RefDip + Core(NAT(i)) * XYZ(:,i)
  enddo

  deallocate(XYZ)
  deallocate(Sao)
  deallocate(Rao)
  deallocate(Delao)
  deallocate(RxDelao)
  deallocate(TMP)
  deallocate(TMP2)

  allocate(murabs(mstate))
  allocate(mupabs(mstate))
  allocate(magabs(mstate))
  allocate(murphi(mstate))
  allocate(mupphi(mstate))
  allocate(magphi(mstate))
  allocate(murtheta(mstate))
  allocate(muptheta(mstate))
  allocate(magtheta(mstate))

  do s=2,mstate
    deltaE      = (CiE(s) - CiE(1)) / eV
    mur         = AU2DBY*    Rst(s,2,:)
    mup         = AU2DBY*  Delst(s,2,:) / deltaE
    mag         =         -RxDelst(s,2,:)
    murabs(s)   = Sqrt(dot_Product(mur,mur))
    mupabs(s)   = Sqrt(dot_Product(mup,mup))
    magabs(s)   = Sqrt(dot_Product(mag,mag))
    murphi(s)   = atan2deg(mur(2), mur(1))
    mupphi(s)   = atan2deg(mup(2), mup(1))
    magphi(s)   = atan2deg(mag(2), mag(1))
    murtheta(s) = acos2deg(mur(3), murabs(s))
    muptheta(s) = acos2deg(mup(3), mupabs(s))
    magtheta(s) = acos2deg(mag(3), magabs(s))
    CiProp(7,s) = murabs(s)
    CiProp(8,s) = mupabs(s)
    if (murtheta(s) < 45.D0) Then
       CiProp(10,s) = murphi(s)
       CiProp(11,s) = mupphi(s)
       CiProp(13,s) = murtheta(s)
       CiProp(14,s) = muptheta(s)
    else if (murtheta(s) > 135.D0) Then
       CiProp(10,s) = atan2deg(-mur(2), -mur(1))
       CiProp(11,s) = atan2deg(-mup(2), -mup(1))
       CiProp(13,s) = 180.D0 - murtheta(s)
       CiProp(14,s) = 180.D0 - muptheta(s)
    else if (Abs(murphi(s)) <= 90.D0) Then
       CiProp(10,s) = murphi(s)
       CiProp(11,s) = mupphi(s)
       CiProp(13,s) = murtheta(s)
       CiProp(14,s) = muptheta(s)
    else
       CiProp(10,s) = atan2deg(-mur(2), -mur(1))
       CiProp(11,s) = atan2deg(-mup(2), -mup(1))
       CiProp(13,s) = acos2deg(-mur(3),  murabs(s))
       CiProp(14,s) = acos2deg(-mup(3),  mupabs(s))
    end if
  end do

  ! Print results:
  !---------------

  ! Print state dipole moments:
  if (iuvcd > 0) Then
     write(nb6,'(A//A)') ' State dipole moments:', &
       &   "   #    Symmetry       eV         cm-1     mu_r(x)/D   mu_r(y)/D   mu_r(z)/D    |mu_r|/D     phi/deg   theta/deg"
     do s=1,nroots
        dip    = (RefDip - Rst(s,1,:)) * AU2DBY
        absdip = Sqrt(dot_Product(dip,dip))
        phi    = atan2deg(dip(2), dip(1))
        theta  = acos2deg(dip(3), absdip)
        icisym(   s) = CiStateSym(s)
        CiProp( 6,s) = absdip
        CiProp( 9,s) = phi
        CiProp(12,s) = theta
        CiProp(15,s) = dip(1)
        CiProp(16,s) = dip(2)
        CiProp(17,s) = dip(3)
        write(nb6,'(i4,A9,A1,i1,A1,F12.6,F12.2,4F12.6,2F12.4)')    &
          &   s, SymLabel(CiStateSym(s)), '(', CiStateSym(s), ')', &
          &   CiE(s)-CiE(1), AU2WAV*(CiE(s)-CiE(1))/eV, dip, absdip, phi, theta
     end do
     write(nb6,'(//)')
  end if

  ! Print oscillator and rotational strengths for transitions originating from state 1:
  if (iuvcd > 1) Then
     write(nb6,'(A//A)') " Properties of transitions   1 -> #:", &
       &   "   #    Symmetry       eV         cm-1        nm        f_r         f_p         f_rp        R/DBM"
     do t=2,nroots
        deltaE = (CiE(t) - CiE(1)) / eV
        mup    = -AU2DBY * Delst(t,2,:) / deltaE
        fr     =  dot_Product(Rst(t,2,:),Rst(t,2,:))*deltaE/1.5D0
        fp     =  dot_Product(Delst(t,2,:),Delst(t,2,:))/(deltaE*1.5D0)
        frp    =  dot_Product(Rst(t,2,:),Delst(t,2,:))/1.5D0
        rot    = -dot_Product(mup,RxDelst(t,2,:))

           CiProp(2,t) = rot
           CiProp(3,t) = fr
           CiProp(4,t) = fp
           CiProp(5,t) = frp

        write(nb6,'(i4,A9,A1,i1,A1,F12.6,F12.2,F10.2,4F12.6)')     &
          &   t, SymLabel(CiStateSym(t)), '(', CiStateSym(t), ')', &
          &   CiE(t)-CiE(1), AU2WAV*deltaE, 1.D7/(AU2WAV*deltaE), fr, fp, frp, rot
     end do
     write(nb6,'(//)')
  end if

  ! Print transition moments:
  if (iuvcd > 1) then
    write(nb6,'(A,i4,A//A)') " Dipole-length electric dipole transition moments", 1, " -> #:", &
      &   "   #    Symmetry       eV         cm-1     mu_r(x)/D   mu_r(y)/D   mu_r(z)/D    |mu_r|/D     phi/deg   theta/deg"
    do s=2,nroots
    deltaE = (CiE(s) - CiE(1)) / eV
    mur    = AU2DBY * Rst(s,2,:)
    write(nb6,'(i4,A9,A1,i1,A1,F12.6,F12.2,4F12.6,2F12.4)')  &
      & s, SymLabel(CiStateSym(s)), '(', CiStateSym(s), ')', &
      & CiE(s)-CiE(1), AU2WAV*deltaE, mur, murabs(s), murphi(s), murtheta(s)
    enddo
    write(nb6,'(//)')

    write(nb6,'(A,i4,A//A)') " Dipole-velocity electric dipole transition moments", 1, " -> #:", &
      &   "   #    Symmetry       eV         cm-1     mu_p(x)/D   mu_p(y)/D   mu_p(z)/D    |mu_p|/D     phi/deg   theta/deg"
    do s=2,nroots
    deltaE = (CiE(s) - CiE(1)) / eV
    mup    = AU2DBY * Delst(s,2,:) / deltaE
    write(nb6,'(i4,A9,A1,i1,A1,F12.6,F12.2,4F12.6,2F12.4)')    &
      & s, SymLabel(CiStateSym(s)), '(', CiStateSym(s), ')', &
      & CiE(s)-CiE(1), AU2WAV*deltaE, mup, mupabs(s), mupphi(s), muptheta(s)
    enddo
    write(nb6,'(//)')

    write(nb6,'(A,i4,A//A)') " Magnetic dipole transition moments", 1, " -> #:", &
      &   "   #    Symmetry       eV         cm-1      m(x)/iBM    m(y)/iBM    m(z)/iBM      |m|/BM     phi/deg   theta/deg"
    do s=2,nroots
    deltaE = (CiE(s) - CiE(1)) / eV
    mag    = -RxDelst(s,2,:)
    write(nb6,'(i4,A9,A1,i1,A1,F12.6,F12.2,4F12.6,2F12.4)')    &
      & s, SymLabel(CiStateSym(s)), '(', CiStateSym(s), ')', &
      & CiE(s)-CiE(1), AU2WAV*deltaE, mag, magabs(s), magphi(s), magtheta(s)
    enddo
    write(nb6,'(//)')
  endif

  deallocate(murabs)
  deallocate(mupabs)
  deallocate(magabs)
  deallocate(murphi)
  deallocate(mupphi)
  deallocate(magphi)
  deallocate(murtheta)
  deallocate(muptheta)
  deallocate(magtheta)

  deallocate(Rst)
  deallocate(Delst)
  deallocate(RxDelst)
  return

  contains

  double Precision Function atan2deg(y,x)
    double Precision            :: y, x
    ! end of dummy parameters
    double Precision            :: r
    double Precision, Parameter :: Pi = 3.1415926535897932385D0
    r = Sqrt(x*x + y*y)
    if (dabs(r) > 1.D-6) Then
       atan2deg = 180.D0 * atan2(y,x) / Pi
    else
       atan2deg = 0.D0
    end if
    Return
  end Function atan2deg



  double Precision Function acos2deg(z,r)
    double Precision :: z, r
    ! end of dummy parameters
    double Precision, Parameter :: Pi = 3.1415926535897932385D0
    if (dabs(r) > 1.D-6) Then
       if (z > r) Then
          acos2deg = 0.D0
       else if (dabs(z) > r) Then
          acos2deg = 180.D0
       else
          acos2deg = 180.D0 * acos(z/r) / Pi
       end if
    else
       acos2deg = 0.D0
    end if
    Return
  end Function acos2deg

end

