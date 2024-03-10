!
! Calculation of spectroscopic properties (CD, UV) of a set of CI states.
!
! Written by Axel Koslowski in 2000-2005 at MPI Muelheim.
!

Module property

  Use Const, Only: A0, AU2DBY, AU2WAV, EV, OLDCF
  Use matutils

  Implicit None

  Private

  Public :: Properties, Population, NaturalOrb

Contains

  !
  ! The following subroutine is the public driver routine for the calculation of the
  ! spectroscopic properties (i.e. UV, CD, and in a future version also MCD) of the
  ! molecule:
  !
  Subroutine Properties(CIE, CIStateSym, SymLabel, Perm, Trans, CMO, NAT, CoAng,        &
                        NFirst, NLast, Core, zs, zp, zd, nsp, nd, Occ, iactive, CIProp, &
                        icisym, nPerm, mTrans, nTrans, NB6, iop, iuvcd, imcd, icall)
    Implicit None
    Double Precision, Dimension(:)     :: CIE               ! CI eigenvalues.
    Integer,          Dimension(:)     :: CIStateSym        ! Symmetry numbers of CI states.
    Character(Len=4), Dimension(:)     :: SymLabel          ! Symmetry labels.
    Double Precision, Dimension(:,:,:) :: Perm              ! Permanent one-electron densities.
    Double Precision, Dimension(:,:,:) :: Trans             ! One-electron transition densities.
    Double Precision, Dimension(:,:)   :: CMO               ! MO coefficients.
    Integer,          Dimension(:)     :: NAT               ! Nuclear charges of all atoms.
    Double Precision, Dimension(:,:)   :: CoAng             ! Atom coordinates (Angstroem), Dimension(3,:).
    Integer,          Dimension(:)     :: NFirst, NLast     ! AO index range of each atom.
    Double Precision, Dimension(:)     :: Core              ! Core charge (nucleus plus core electrons) of each element.
    Double Precision, Dimension(:)     :: zs, zp, zd        ! Atomic orbital exponents per element.
    Integer,          Dimension(:)     :: nsp, nd           ! Principle quantum numbers of AOs per element.
    Double Precision, Dimension(:)     :: Occ               ! Orbital occupation numbers in SCF reference.
    Integer,          Dimension(:)     :: iactive           ! Indices of the active orbitals.
    Double Precision, Dimension(:,:)   :: CIProp            ! Array to store the calculated properties.
                                                            ! The first index specifies the property,
                                                            ! the second one the state.
                                                            ! The codes for the first index are:
                                                            !  1: excitation energy (eV),
                                                            !  2: rotational strength (DBM),
                                                            !  3: dipole-length oscillator strength (AU),
                                                            !  4: dipole-velocity oscillator strength (AU),
                                                            !  5: mixed oscillator strength (AU),
                                                            !  6: permanent dipole moment magnitude (D),
                                                            !  7: dipole-length dipole transition moment magnitude (D),
                                                            !  8: dipole-velocity dipole transition moment magnitude (D),
                                                            !  9: permanent dipole moment angle phi (deg),
                                                            ! 10: dipole-length dipole transition moment angle phi (deg),
                                                            ! 11: dipole-velocity dipole transition moment angle phi (deg),
                                                            ! 12: permanent dipole moment angle theta (deg),
                                                            ! 13: dipole-length dipole transition moment angle theta (deg),
                                                            ! 14: dipole-velocity dipole transition moment angle theta (deg),
                                                            ! 15: permanent dipole moment x coordinate,
                                                            ! 16: permanent dipole moment y coordinate,
                                                            ! 17: permanent dipole moment z coordinate.
                                                            !
    Integer,          Dimension(:)   :: icisym              ! Symmetry numbers of the LMSTAT lowest CI states.
    Integer                          :: nPerm               ! Number of permanent one-electron densities.
    Integer                          :: mTrans              ! Highest state from which transitions may originate.
    Integer                          :: nTrans              ! Highest state into which transitions may occur.
    Integer                          :: NB6                 ! Unit number of stdout.
    Integer                          :: iop                 ! Semiempirical method, IN2(1).
    Integer                          :: iuvcd               ! Flag controlling UV and CD calculation:
                                                            !  0: no calculation,
                                                            !  1: permanent dipole moments of all states,
                                                            !  2: also spectroscopic properties for transitions,
                                                            !     originating from the ground state,
                                                            !  3: also transition moments between excited states,
                                                            !  4: debug print.
    Integer                          :: imcd                ! Flag controlling MCD calculation (ignored)
                                                            !  0: no calculation.
    Integer                          :: icall               ! Set to -1 to indicate error.
                                                            !
    ! end of dummy parameters
    Double Precision, Dimension(:,:),   Allocatable ::   Coord  ! Atom coordinates (Bohr), Dimension(3,:).
    Double Precision, Dimension(:,:),   Allocatable ::     Sao  ! AO basis operator matrices.
    Double Precision, Dimension(:,:,:), Allocatable ::     Rao
    Double Precision, Dimension(:,:,:), Allocatable ::   Delao
    Double Precision, Dimension(:,:,:), Allocatable :: RxDelao
    Double Precision, Dimension(:,:,:), Allocatable ::     Rmo  ! MO basis operator matrices.
    Double Precision, Dimension(:,:,:), Allocatable ::   Delmo
    Double Precision, Dimension(:,:,:), Allocatable :: RxDelmo
    Double Precision, Dimension(:,:,:), Allocatable ::     Rst  ! State basis operator matrices.
    Double Precision, Dimension(:,:,:), Allocatable ::   Delst
    Double Precision, Dimension(:,:,:), Allocatable :: RxDelst
    Double Precision, Dimension(:),     Allocatable ::  Seigen  ! Eigenvalues of overlap matrix.
    Double Precision, Dimension(:,:),   Allocatable ::  Deorth  ! MOs in non-orthogonal AO basis.
    Double Precision, Dimension(:,:),   Allocatable ::     tmp  ! Scratch matrix.
    Double Precision, Dimension(3)                  ::  RefDip  ! Electronic dipole moment of the reference function.
    Double Precision, Dimension(3)                  ::    Rmoi  ! Dipole moment contribution of occupied orbital i.
    Double Precision, Dimension(3)                  ::     dip  ! Current permanent dipole moment.
    Double Precision, Dimension(3)                  ::     mur  ! Current dipole-length electric transition moment.
    Double Precision, Dimension(3)                  ::     mup  ! Current dipole-velocity electric transition moment.
    Double Precision, Dimension(3)                  ::     mag  ! Current magnetic transition moment.
    Double Precision, Dimension(:,:),   Allocatable :: murabs, murphi, murtheta  ! Transition moment magnitudes and angles.
    Double Precision, Dimension(:,:),   Allocatable :: mupabs, mupphi, muptheta
    Double Precision, Dimension(:,:),   Allocatable :: magabs, magphi, magtheta
    Integer :: natom    ! number of atoms
    Integer :: nao      ! number of atomic orbitals
    Integer :: nmo      ! number of molecular orbitals
    Integer :: nocc     ! number of occupied molecular orbitals
    Integer :: nactive  ! number of active molecular orbitals
    Integer :: mstate   ! leading dimension of state basis operators
    Integer :: nstate   ! second dimension of state basis operators
    Double Precision :: x, deltaE, fr, fp, frp, rot, absdip, phi, theta
    Integer :: m, i, j, ij, k, s, t, ierr, LMZ
    Logical :: ecp

    !Double Precision, Parameter :: bohr   =  0.529177249D0
    !Double Precision, Parameter :: eV     =  0.272113961D2
    !Double Precision, Parameter :: kayser =  0.219474631D6
    !Double Precision, Parameter :: debye  =  0.254174776D1

    If (iuvcd <= 0) Return

    CALL OLDCF

    natom    = UBound(NAT,1)
    nao      = NLast(natom)
    nmo      = UBound(CMO,2)
    nocc     = UBound(Occ,1)
    nactive  = UBound(iactive,1)
    mstate   = Max(nPerm, mTrans)
    nstate   = Max(nPerm, nTrans)

    Write(NB6,'(//)')

    If (iuvcd >= 4) Then
       Do i=1, nPerm
          Write(NB6,'(A,I4,A/)') ' One-electron density for state', i, ':'
          Call PrintMatrix(Perm(:,:,i), NB6)
       End Do

       ij = 0
       Do i=1, mTrans
          Do j=i+1, nTrans
             ij = ij + 1
             Write(NB6,'(A,I4,A,I4,A/)') ' One-electron density for transition', i, ' ->', j, ':'
             Call PrintMatrix(Trans(:,:,ij), NB6)
          End Do
       End Do
    End If

    Allocate(Coord(3,natom))
    Allocate(Sao(nao,nao))
    Allocate(Rao(nao,nao,3))
    Allocate(Delao(nao,nao,3))
    Allocate(RxDelao(nao,nao,3))

    Coord = CoAng / A0


    ! Calculate AO basis integrals:
    !------------------------------

    ecp = iop == -5 .Or. iop == -6 .Or. iop == -8 .Or. iop == -9 .Or. iop == -22 .Or. iop == -23
    LMZ = Size(zs)
    Call CalcPropInt(Sao, Rao, Delao, RxDelao, NAT, NFirst, NLast, Coord, &
                     zs, zp, zd, nsp, nsp, nd, LMZ, natom, nao, NB6, ecp)

    ! Print matrices if requested:
    If (iuvcd >= 4) Then
       Write(NB6,'(A/)') ' AO basis overlap integrals:'
       Call PrintMatrix(Sao, NB6)
       Write(NB6,'(A/)') ' AO basis R component x:'
       Call PrintMatrix(Rao(:,:,1), NB6)
       Write(NB6,'(A/)') ' AO basis R component y:'
       Call PrintMatrix(Rao(:,:,2), NB6)
       Write(NB6,'(A/)') ' AO basis R component z:'
       Call PrintMatrix(Rao(:,:,3), NB6)
       Write(NB6,'(A/)') ' AO basis DEL component x:'
       Call PrintMatrix(Delao(:,:,1), NB6)
       Write(NB6,'(A/)') ' AO basis DEL component y:'
       Call PrintMatrix(Delao(:,:,2), NB6)
       Write(NB6,'(A/)') ' AO basis DEL component z:'
       Call PrintMatrix(Delao(:,:,3), NB6)
       Write(NB6,'(A/)') ' AO basis R x DEL component x:'
       Call PrintMatrix(RxDelao(:,:,1), NB6)
       Write(NB6,'(A/)') ' AO basis R x DEL component y:'
       Call PrintMatrix(RxDelao(:,:,2), NB6)
       Write(NB6,'(A/)') ' AO basis R x DEL component z:'
       Call PrintMatrix(RxDelao(:,:,3), NB6)
    End If

    ! Calculate S**(-0.5) matrix and deorthogonalize MOs:
    !----------------------------------------------------

    Allocate(Seigen(nao))
    Allocate(tmp(nao,nao))

    ! Diagonalize overlap matrix:
    Call TRED2(nao, nao, Sao, Seigen, tmp(:,1), Sao)
    Call TQL2 (nao, nao,      Seigen, tmp(:,1), Sao, ierr)

    If (ierr /= 0) Then
       Write(NB6,'(1X,A)') 'The diagonalization of the overlap matrix failed.'
       Deallocate(Coord)
       Deallocate(Sao)
       Deallocate(Rao)
       Deallocate(Delao)
       Deallocate(RxDelao)
       Deallocate(Seigen)
       Deallocate(tmp)
       icall = -1
       Return
    EndIf

    If (iuvcd >= 4) Then
       Write(NB6,'(A/)') ' Eigenvalues of overlap matrix:'
       Call PrintVector(Seigen, NB6)
    End If

    ! Form square root of inverse of all eigenvalues:
    Seigen = 1.D0 / Sqrt(Seigen)

    ! Generate S**(-0.5) matrix in tmp:
    Do i=1,nao
       Do j=i,nao
          x = 0.D0
          Do k=1,nao
             x = x + Sao(i,k) * Sao(j,k) * Seigen(k)
          End Do
          tmp(i,j) = x
          tmp(j,i) = x
       End Do
    End Do

    Deallocate(Sao)
    Deallocate(Seigen)

    ! Transform MOs:
    Allocate(Deorth(nao,nmo))
    Deorth = Matmul(tmp, CMO(1:nao,1:nmo))

    ! Print S**(-0.5) matrix and MO coefficients if requested:
    If (iuvcd >= 4) Then
       Write(NB6,'(A/)') ' S**(-1/2) matrix:'
       Call PrintMatrix(tmp, NB6)
       Write(NB6,'(A/)') ' MO coefficients in orthogonal AO basis:'
       Call PrintMatrix(CMO, NB6)
       Write(NB6,'(A/)') ' MO coefficients in deorthogonal AO basis:'
       Call PrintMatrix(Deorth, NB6)
    End If

    Deallocate(tmp)

    ! Transform integrals from AO to MO basis:
    !-----------------------------------------

    ! Reference dipole moment:

    RefDip = 0.D0

    ! Contribution of the atom cores:
    Do i=1, natom
       RefDip = RefDip + Core(NAT(i)) * Coord(:,i)
    End Do

    If (iuvcd >= 4) Then
       Write(NB6,'(A/)') ' Dipole moment of the atom cores (Debye):'
       Write(NB6,'(3F12.6//)') RefDip * AU2DBY
    End If

    ! Electronic contribution to the reference dipole moment:
    Do i=1, nocc
       Call TransformMV(Rmoi(1), Rao(:,:,1), Deorth(:,i))
       Call TransformMV(Rmoi(2), Rao(:,:,2), Deorth(:,i))
       Call TransformMV(Rmoi(3), Rao(:,:,3), Deorth(:,i))
       RefDip = RefDip - Rmoi * Occ(i)
    End Do

    ! Active orbitals:

    Allocate(Rmo(nactive,nactive,3))
    Call TransformMM(Rmo(:,:,1), Rao(:,:,1), Deorth(:,iactive))
    Call TransformMM(Rmo(:,:,2), Rao(:,:,2), Deorth(:,iactive))
    Call TransformMM(Rmo(:,:,3), Rao(:,:,3), Deorth(:,iactive))
    Deallocate(Rao)

    Allocate(Delmo(nactive,nactive,3))
    Call TransformMM(Delmo(:,:,1), Delao(:,:,1), Deorth(:,iactive))
    Call TransformMM(Delmo(:,:,2), Delao(:,:,2), Deorth(:,iactive))
    Call TransformMM(Delmo(:,:,3), Delao(:,:,3), Deorth(:,iactive))
    Deallocate(Delao)

    Allocate(RxDelmo(nactive,nactive,3))
    Call TransformMM(RxDelmo(:,:,1), RxDelao(:,:,1), Deorth(:,iactive))
    Call TransformMM(RxDelmo(:,:,2), RxDelao(:,:,2), Deorth(:,iactive))
    Call TransformMM(RxDelmo(:,:,3), RxDelao(:,:,3), Deorth(:,iactive))
    Deallocate(RxDelao)

    ! Print matrices if requested:
    If (iuvcd >= 4) Then
       Write(NB6,'(A/)') ' Dipole moment of the reference function (Debye):'
       Write(NB6,'(3F12.6//)') RefDip * AU2DBY
       Write(NB6,'(A/)') ' Active MO basis R component x:'
       Call PrintMatrix(Rmo(:,:,1), NB6)
       Write(NB6,'(A/)') ' Active MO basis R component y:'
       Call PrintMatrix(Rmo(:,:,2), NB6)
       Write(NB6,'(A/)') ' Active MO basis R component z:'
       Call PrintMatrix(Rmo(:,:,3), NB6)
       Write(NB6,'(A/)') ' Active MO basis DEL component x:'
       Call PrintMatrix(Delmo(:,:,1), NB6)
       Write(NB6,'(A/)') ' Active MO basis DEL component y:'
       Call PrintMatrix(Delmo(:,:,2), NB6)
       Write(NB6,'(A/)') ' Active MO basis DEL component z:'
       Call PrintMatrix(Delmo(:,:,3), NB6)
       Write(NB6,'(A/)') ' Active MO basis R x DEL component x:'
       Call PrintMatrix(RxDelmo(:,:,1), NB6)
       Write(NB6,'(A/)') ' Active MO basis R x DEL component y:'
       Call PrintMatrix(RxDelmo(:,:,2), NB6)
       Write(NB6,'(A/)') ' Active MO basis R x DEL component z:'
       Call PrintMatrix(RxDelmo(:,:,3), NB6)
    End If

    ! Calculate state basis integrals:
    !---------------------------------

    Allocate(    Rst(mstate,nstate,3))
    Allocate(  Delst(mstate,nstate,3))
    Allocate(RxDelst(mstate,nstate,3))

        Rst = 0.D0
      Delst = 0.D0
    RxDelst = 0.D0

    Do m=1,nPerm
       Rst(m,m,1) = Contract(Perm(:,:,m), Rmo(:,:,1))
       Rst(m,m,2) = Contract(Perm(:,:,m), Rmo(:,:,2))
       Rst(m,m,3) = Contract(Perm(:,:,m), Rmo(:,:,3))
    End Do

    m = 0
    Do s=1,mTrans
       Do t=s+1,nTrans
          m = m + 1
              Rst(s,t,1) = Contract(Trans(:,:,m),     Rmo(:,:,1))
              Rst(s,t,2) = Contract(Trans(:,:,m),     Rmo(:,:,2))
              Rst(s,t,3) = Contract(Trans(:,:,m),     Rmo(:,:,3))
            Delst(s,t,1) = Contract(Trans(:,:,m),   Delmo(:,:,1))
            Delst(s,t,2) = Contract(Trans(:,:,m),   Delmo(:,:,2))
            Delst(s,t,3) = Contract(Trans(:,:,m),   Delmo(:,:,3))
          RxDelst(s,t,1) = Contract(Trans(:,:,m), RxDelmo(:,:,1))
          RxDelst(s,t,2) = Contract(Trans(:,:,m), RxDelmo(:,:,2))
          RxDelst(s,t,3) = Contract(Trans(:,:,m), RxDelmo(:,:,3))
       End Do
    End Do

    Deallocate(Deorth)
    Deallocate(Rmo)
    Deallocate(Delmo)
    Deallocate(RxDelmo)

    ! Print matrices if requested:
    If (iuvcd >= 4) Then
       Write(NB6,'(A/)') ' State basis R component x:'
       Call PrintUpper(Rst(:,:,1), NB6)
       Write(NB6,'(A/)') ' State basis R component y:'
       Call PrintUpper(Rst(:,:,2), NB6)
       Write(NB6,'(A/)') ' State basis R component z:'
       Call PrintUpper(Rst(:,:,3), NB6)
       Write(NB6,'(A/)') ' State basis DEL component x:'
       Call PrintUpper(Delst(:,:,1), NB6)
       Write(NB6,'(A/)') ' State basis DEL component y:'
       Call PrintUpper(Delst(:,:,2), NB6)
       Write(NB6,'(A/)') ' State basis DEL component z:'
       Call PrintUpper(Delst(:,:,3), NB6)
       Write(NB6,'(A/)') ' State basis R x DEL component x:'
       Call PrintUpper(RxDelst(:,:,1), NB6)
       Write(NB6,'(A/)') ' State basis R x DEL component y:'
       Call PrintUpper(RxDelst(:,:,2), NB6)
       Write(NB6,'(A/)') ' State basis R x DEL component z:'
       Call PrintUpper(RxDelst(:,:,3), NB6)
    End If

    ! Calculate transition moment magnitudes and angles:
    !---------------------------------------------------

    Allocate(murabs(mstate,nstate))
    Allocate(mupabs(mstate,nstate))
    Allocate(magabs(mstate,nstate))
    Allocate(murphi(mstate,nstate))
    Allocate(mupphi(mstate,nstate))
    Allocate(magphi(mstate,nstate))
    Allocate(murtheta(mstate,nstate))
    Allocate(muptheta(mstate,nstate))
    Allocate(magtheta(mstate,nstate))

    m = 0
    Do s=1,mTrans
       Do t=s+1,nTrans
          m = m + 1
          deltaE        = (CIE(t) - CIE(s)) / eV
          mur           = -AU2DBY*    Rst(s,t,:)
          mup           = -AU2DBY*  Delst(s,t,:) / deltaE
          mag           =         RxDelst(s,t,:)
          murabs(s,t)   = Sqrt(Dot_Product(mur,mur))
          mupabs(s,t)   = Sqrt(Dot_Product(mup,mup))
          magabs(s,t)   = Sqrt(Dot_Product(mag,mag))
          murphi(s,t)   = atan2deg(mur(2), mur(1))
          mupphi(s,t)   = atan2deg(mup(2), mup(1))
          magphi(s,t)   = atan2deg(mag(2), mag(1))
          murtheta(s,t) = acos2deg(mur(3), murabs(s,t))
          muptheta(s,t) = acos2deg(mup(3), mupabs(s,t))
          magtheta(s,t) = acos2deg(mag(3), magabs(s,t))
          If (s == 1 .And. t <= UBound(CIProp,2)) Then
             CIProp(7,t) = murabs(s,t)
             CIProp(8,t) = mupabs(s,t)
             If (murtheta(s,t) < 45.D0) Then
                CIProp(10,t) = murphi(s,t)
                CIProp(11,t) = mupphi(s,t)
                CIProp(13,t) = murtheta(s,t)
                CIProp(14,t) = muptheta(s,t)
             Else If (murtheta(s,t) > 135.D0) Then
                CIProp(10,t) = atan2deg(-mur(2), -mur(1))
                CIProp(11,t) = atan2deg(-mup(2), -mup(1))
                CIProp(13,t) = 180.D0 - murtheta(s,t)
                CIProp(14,t) = 180.D0 - muptheta(s,t)
             Else If (Abs(murphi(s,t)) <= 90.D0) Then
                CIProp(10,t) = murphi(s,t)
                CIProp(11,t) = mupphi(s,t)
                CIProp(13,t) = murtheta(s,t)
                CIProp(14,t) = muptheta(s,t)
             Else
                CIProp(10,t) = atan2deg(-mur(2), -mur(1))
                CIProp(11,t) = atan2deg(-mup(2), -mup(1))
                CIProp(13,t) = acos2deg(-mur(3),  murabs(s,t))
                CIProp(14,t) = acos2deg(-mup(3),  mupabs(s,t))
             End If
          End If
       End Do
    End Do

    ! Print results:
    !---------------

    ! Print state dipole moments:
    If (nPerm > 0) Then
       Write(NB6,'(A//A)') ' State dipole moments:', &
         &   "   #    Symmetry       eV         cm-1     mu_r(x)/D   mu_r(y)/D   mu_r(z)/D    |mu_r|/D     phi/deg   theta/deg"
       Do s=1,nPerm
          dip    = (RefDip - Rst(s,s,:)) * AU2DBY
          absdip = Sqrt(Dot_Product(dip,dip))
          phi    = atan2deg(dip(2), dip(1))
          theta  = acos2deg(dip(3), absdip)
          If (s <= Ubound(CIProp,2)) Then
             icisym(   s) = CIStateSym(s)
             CIProp( 6,s) = absdip
             CIProp( 9,s) = phi
             CIProp(12,s) = theta
             CIProp(15,s) = dip(1)
             CIProp(16,s) = dip(2)
             CIProp(17,s) = dip(3)
          End If
          Write(NB6,'(I4,A9,A1,I1,A1,F12.6,F12.2,4F12.6,2F12.4)')    &
            &   s, SymLabel(CIStateSym(s)), '(', CIStateSym(s), ')', &
            &   CIE(s)-CIE(1), AU2WAV*(CIE(s)-CIE(1))/eV, dip, absdip, phi, theta
       End Do
       Write(NB6,'(//)')
    End If

    ! Print oscillator and rotational strengths for transitions originating from state 1:
    If (nTrans > 1) Then
       Write(NB6,'(A//A)') " Properties of transitions   1 -> #:", &
         &   "   #    Symmetry       eV         cm-1        nm        f_r         f_p         f_rp        R/DBM"
       Do t=2,nTrans
          deltaE = (CIE(t) - CIE(1)) / eV
          mup    = -AU2DBY * Delst(1,t,:) / deltaE
          fr     =  Dot_Product(Rst(1,t,:),Rst(1,t,:))*deltaE/1.5D0
          fp     =  Dot_Product(Delst(1,t,:),Delst(1,t,:))/(deltaE*1.5D0)
          frp    =  Dot_Product(Rst(1,t,:),Delst(1,t,:))/1.5D0
          rot    = -Dot_Product(mup,RxDelst(1,t,:))

          If (t <= UBound(CIProp,2)) Then
             CIProp(2,t) = rot
             CIProp(3,t) = fr
             CIProp(4,t) = fp
             CIProp(5,t) = frp
          End If

          Write(NB6,'(I4,A9,A1,I1,A1,F12.6,F12.2,F10.2,4F12.6)')     &
            &   t, SymLabel(CIStateSym(t)), '(', CIStateSym(t), ')', &
            &   CIE(t)-CIE(1), AU2WAV*deltaE, 1.D7/(AU2WAV*deltaE), fr, fp, frp, rot
       End Do
       Write(NB6,'(//)')
    End If

    ! Print transition moments:
    Do s=1,mTrans
       Write(NB6,'(A,I4,A//A)') " Dipole-length electric dipole transition moments", s, " -> #:", &
         &   "   #    Symmetry       eV         cm-1     mu_r(x)/D   mu_r(y)/D   mu_r(z)/D    |mu_r|/D     phi/deg   theta/deg"
       Do t=s+1,nTrans
          deltaE = (CIE(t) - CIE(s)) / eV
          mur    = -AU2DBY * Rst(s,t,:)
          Write(NB6,'(I4,A9,A1,I1,A1,F12.6,F12.2,4F12.6,2F12.4)')  &
            & t, SymLabel(CIStateSym(t)), '(', CIStateSym(t), ')', &
            & CIE(t)-CIE(s), AU2WAV*deltaE, mur, murabs(s,t), murphi(s,t), murtheta(s,t)
       End Do
       Write(NB6,'(//)')
       Write(NB6,'(A,I4,A//A)') " Dipole-velocity electric dipole transition moments", s, " -> #:", &
         &   "   #    Symmetry       eV         cm-1     mu_p(x)/D   mu_p(y)/D   mu_p(z)/D    |mu_p|/D     phi/deg   theta/deg"
       Do t=s+1,nTrans
          deltaE = (CIE(t) - CIE(s)) / eV
          mup    = -AU2DBY * Delst(s,t,:) / deltaE
          Write(NB6,'(I4,A9,A1,I1,A1,F12.6,F12.2,4F12.6,2F12.4)')    &
            &   t, SymLabel(CIStateSym(t)), '(', CIStateSym(t), ')', &
            &   CIE(t)-CIE(s), AU2WAV*deltaE, mup, mupabs(s,t), mupphi(s,t), muptheta(s,t)
       End Do
       Write(NB6,'(//)')
       Write(NB6,'(A,I4,A//A)') " Magnetic dipole transition moments", s, " -> #:", &
         &   "   #    Symmetry       eV         cm-1      m(x)/iBM    m(y)/iBM    m(z)/iBM      |m|/BM     phi/deg   theta/deg"
       Do t=s+1,nTrans
          deltaE = (CIE(t) - CIE(s)) / eV
          mag    = RxDelst(s,t,:)
          Write(NB6,'(I4,A9,A1,I1,A1,F12.6,F12.2,4F12.6,2F12.4)')    &
            &   t, SymLabel(CIStateSym(t)), '(', CIStateSym(t), ')', &
            &   CIE(t)-CIE(s), AU2WAV*deltaE, mag, magabs(s,t), magphi(s,t), magtheta(s,t)
       End Do
       Write(NB6,'(//)')
    End Do

    Deallocate(murabs)
    Deallocate(mupabs)
    Deallocate(magabs)
    Deallocate(murphi)
    Deallocate(mupphi)
    Deallocate(magphi)
    Deallocate(murtheta)
    Deallocate(muptheta)
    Deallocate(magtheta)

    Deallocate(Rst)
    Deallocate(Delst)
    Deallocate(RxDelst)
    Deallocate(Coord)

    Return
  End Subroutine Properties



  !
  ! The following subroutine is the public driver routine for performing a
  ! population analysis of state iState or all states calculated previously
  ! (as controlled by the ipop flag).
  !
  Subroutine Population(Perm, PCIAO, CMO, NAT, NFirst, NLast, Core, Occ, iactive, nPerm, NB6, iState, ipop, icall)
    Implicit None
    Double Precision, Dimension(:,:,:)            :: Perm           ! Permanent one-electron densities.
    Double Precision, Dimension(:,:)              :: PCIAO          ! AO-basis one-electron density matrix (iState).
    Double Precision, Dimension(:,:)              :: CMO            ! MO coefficients.
    Integer,          Dimension(:)                :: NAT            ! Nuclear charges of all atoms.
    Integer,          Dimension(:)                :: NFirst, NLast  ! AO index range of each atom.
    Double Precision, Dimension(:)                :: Core           ! Number of valence electrons of each atom.
    Double Precision, Dimension(:)                :: Occ            ! Orbital occupation numbers in SCF reference.
    Integer,          Dimension(:)                :: iactive        ! Indices of the active orbitals.
    Integer                                       :: nPerm          ! Number of permanent one-electron densities.
    Integer                                       :: NB6            ! Unit number of stdout.
    Integer                                       :: iState         ! State of interest (non-negative, actual index).
    Integer                                       :: ipop           ! Flag controlling population analysis:
                                                                    !  0: no such evaluation,
                                                                    !  1: population analysis for state iState,
                                                                    !  2: population analysis for all states
                                                                    !     calculated previously,
                                                                    !  3: debug print.
                                                                    ! -1: special option: calculate PCIAO only.
    Integer                                       :: icall          ! Set to -1 to indicate error.
    ! end of dummy parameters

    Double Precision, Dimension(:,:), Allocatable :: CACT
    Double Precision, Dimension(:,:), Allocatable :: PAO
    Double Precision, Dimension(:,:), Allocatable :: temp
    Double Precision, Dimension(:),   Allocatable :: Q
    Integer                                       :: natom, norb, nocc, nact
    Integer                                       :: kfirst, klast, i, j, k, mu, nu

    If (ipop == 1 .Or. ipop < 0) Then
       kfirst = iState
       klast  = iState
    Else
       kfirst = 1
       klast  = nPerm
    End If

    natom = UBound(NAT,1)
    norb  = UBound(CMO,2)
    nocc  = UBound(Occ,1)
    nact  = UBound(iactive,1)

    Allocate(CACT(norb,nact))
    Allocate(PAO(norb,norb))
    Allocate(temp(nact,norb))
    Allocate(Q(natom))

    ! Copy CMO to CACT:
    CACT = CMO(1:norb, iactive)

    Do k=kfirst, klast
       If (ipop > 0) Then
          Write(NB6,'(//1X,"Population analysis of state",I3,":"/)') k
       End If

       If (ipop >= 3) Then
          Write(NB6,'(1X,"Active MO basis one-electron density matrix:"/)')
          Call PrintMatrix(Perm(:,:,k), NB6)
       End If

       ! Perm -> PAO transformation:
       Call DGEMM('N', 'T', nact, norb, nact, 1.D0, Perm(1,1,k), nact, CACT(1,1), norb, 0.D0, temp(1,1), nact)
       Call DGEMM('N', 'N', norb, norb, nact, 1.D0, CACT(1,1),   norb, temp(1,1), nact, 0.D0, PAO (1,1), norb)

       ! Occupied orbital contributions to PAO:
       Do nu=1, norb
          Do mu=1, norb
             Do i=1, nocc
                PAO(mu,nu) = PAO(mu,nu) + CMO(mu,i) * CMO(nu,i) * Occ(i)
             End Do
          End Do
       End Do

       ! Pass PAO of state iState to calling routine:
       If (k == iState) PCIAO = PAO

       ! Print PAO:
       If (ipop >= 3) Then
          Write(NB6,'(1X,"AO basis one-electron density matrix:"/)')
          Call PrintMatrix(PAO, NB6)
       End If

       If (ipop > 0) Then

          ! Calculate net atomic charges:
          Do i=1, natom
             Q(i) = Core(NAT(i))
             Do j=NFirst(i), NLast(i)
                Q(i) = Q(i) - PAO(j,j)
             End Do
          End Do

          ! Print net atomic charges:
          Write(NB6,'(/1X,"Net atomic charges (atomic units):"/)')
          Do i=1, natom
             Write(NB6,'(1X,I6,F12.6)') i, Q(i)
          End Do

       End If
    End Do

    Deallocate(Q)
    Deallocate(temp)
    Deallocate(PAO)
    Deallocate(CACT)

    Return
  End Subroutine Population



  Subroutine NaturalOrb(Perm, Occ, iactive, nPerm, NB6, iState, inatur, icall)
    Implicit None
    Double Precision, Dimension(:,:,:)            :: Perm           ! Permanent one-electron densities (active MO basis).
    Double Precision, Dimension(:)                :: Occ            ! Orbital occupation numbers in the SCF reference.
    Integer,          Dimension(:)                :: iactive        ! Indices of the active orbitals.
    Integer                                       :: nPerm          ! Number of permanent one-electron densities.
    Integer                                       :: NB6            ! Unit number of stdout.
    Integer                                       :: iState         ! State of interest (non-negative, actual index).
    Integer                                       :: inatur         ! Flag controlling natural orbital calculation:
                                                                    !  0: no such evaluation,
                                                                    !  1: natural orbitals for state iState,
                                                                    !  2: natural orbitals for all states
                                                                    !     calculated previously,
                                                                    !  3: debug print.
    Integer                                       :: icall          ! Set to -1 to indicate error.
    ! end of dummy parameters

    Double Precision, Dimension(:),   Allocatable :: Psp, Enat, work
    Double Precision, Dimension(:,:), Allocatable :: Cnat
    Integer                                       :: nact, nnact, i, j, ij, k, kfirst, klast, info

    If (inatur <= 0)  Return

    If (inatur == 1) Then
       kfirst = iState
       klast  = iState
    Else
       kfirst = 1
       klast  = nPerm
    End If

    nact  = UBound(iactive,1)
    nnact = nact * (nact+1) / 2

    Allocate(Psp(nnact))
    Allocate(Enat(nact))
    Allocate(Cnat(nact,nact))
    Allocate(work(3 * nact))

    Do k=kfirst, klast
       Write(NB6,'(//1X,"Natural orbitals for state",I3,":"//)') k

       If (inatur >= 3) Then
          Write(NB6,'(1X,"Active MO basis one-electron density matrix (square storage):"/)')
          Call PrintMatrix(Perm(:,:,k), NB6)
       End If

       ij = 0
       Do i=1, nact
          Do j=1, i
             ij      = ij + 1
             Psp(ij) = Perm(i,j,k)
          End Do
          Psp(ij) = Psp(ij) + Occ(iactive(i))
       End Do

       Write(NB6,'(1X,"Active MO basis one-electron density matrix (packed storage):"/)')
       Call PrintPackedLower(Psp, nact, NB6)

       Call DSPEV('V', 'U', nact, Psp(1), Enat(1), Cnat(1,1), nact, work(1), info)

       If (info < 0) Then
          Write(NB6,'(1X,"DSPEV: Argument ",I1," has an illegal value.")') -info
          icall = -1
          Return
       Else If (info > 0) Then
          Write(NB6,'(1X,"DSPEV: The algorithm failed to converge.")')
          icall = -1
          Return
       End If

       Write(NB6,'(1X,"Density matrix eigenvalues:"/)')
       Call PrintVector(Enat, NB6)

       Write(NB6,'(1X,"Density matrix eigenvectors:"/)')
       Call PrintMatrix(Cnat, NB6)
    End Do

    Deallocate(work)
    Deallocate(Cnat)
    Deallocate(Enat)
    Deallocate(Psp)

    Return
  End Subroutine NaturalOrb



  ! C = BT * A * B
  Subroutine TransformMM(C, A, B)
    Double Precision, Dimension(:,:)                      :: C, A, B
    ! end of dummy parameters
    Double Precision, Dimension(UBound(A,1), UBound(C,1)) :: B2, AB
    Integer                                               :: na, nc

    na = UBound(A,1)
    nc = UBound(C,1)
    B2 = B(1:na, 1:nc)

    Call DGEMM('N', 'N', na, nc, na, 1.D0, A (1,1), na, B2(1,1), na, 0.D0, AB(1,1), na)
    Call DGEMM('T', 'N', nc, nc, na, 1.D0, B2(1,1), na, AB(1,1), na, 0.D0, C (1,1), nc)

    Return
  End Subroutine TransformMM



  ! c = bT * A * b
  Subroutine TransformMV(c, A, b)
    Double Precision                         :: c
    Double Precision, Dimension(:,:)         :: A
    Double Precision, Dimension(:)           :: b
    ! end of dummy parameters
    Double Precision, Dimension(UBound(A,1)) :: Ab
    Double Precision                         :: DDOT
    Integer                                  :: na

    na = UBound(A,1)

    Call DGEMV('N', na, na, 1.D0, A(1,1), na, b, 1, 0.D0, Ab, 1)
    c = DDOT(na, b, 1, Ab, 1)

    Return
  End Subroutine TransformMV



  Double Precision Function Contract(A, B)
    Double Precision, Dimension(:,:) :: A, B
    ! end of dummy parameter
    Integer                          :: n, i, j
    Double Precision                 :: x

    n = UBound(A,1)

    If (UBound(A,2)/=n .Or. UBound(B,1)/=n .Or. UBound(B,2)/=n) Then
       Stop 'property::Contract'
    End If

    x = 0.D0
    Do i=1, n
       Do j=1, n
          x = x + A(i,j) * B(i,j)
       End Do
    End Do

    Contract = x
    Return
  End Function Contract



  Double Precision Function atan2deg(y,x)
    Double Precision            :: y, x
    ! end of dummy parameters
    Double Precision            :: r
    Double Precision, Parameter :: PI = 3.1415926535897932385D0
    r = Sqrt(x*x + y*y)
    If (dabs(r) > 1.D-6) Then
       atan2deg = 180.D0 * atan2(y,x) / PI
    Else
       atan2deg = 0.D0
    End If
    Return
  End Function atan2deg



  Double Precision Function acos2deg(z,r)
    Double Precision :: z, r
    ! end of dummy parameters
    Double Precision, Parameter :: PI = 3.1415926535897932385D0
    If (dabs(r) > 1.D-6) Then
       If (z > r) Then
          acos2deg = 0.D0
       Else If (dabs(z) > r) Then
          acos2deg = 180.D0
       Else
          acos2deg = 180.D0 * acos(z/r) / PI
       End If
    Else
       acos2deg = 0.D0
    End If
    Return
  End Function acos2deg

End Module property
