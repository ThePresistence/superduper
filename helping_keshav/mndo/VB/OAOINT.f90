
Module OAOINT

  Private

  Public :: OAOTRA    ! Integral transformation from orthogonal to non-orthogonal AO basis.

Contains

  !
  ! Integral transformation from orthogonal to non-orthogonal AO basis.
  !

  ! Driver routine.
  Subroutine OAOTRA (ActiveCMO, FockMatrix, TwoElectronAOIntegrals,  &
             &       OneElectronMOIntegrals, TwoElectronMOIntegrals, &
             &       NFirst, NLast, PrintLevel, Stdout, Ivbovr)
    Implicit None
    Double Precision, Dimension(:,:), Intent(In)    :: ActiveCMO
    Double Precision, Dimension(:),   Intent(In)    :: FockMatrix
    Double Precision, Dimension(:,:), Intent(In)    :: TwoElectronAOIntegrals
    Double Precision, Dimension(:,:), Intent(Out)   :: OneElectronMOIntegrals
    Double Precision, Dimension(:,:), Intent(Out)   :: TwoElectronMOIntegrals
    Integer, Dimension(:),            Intent(In)    :: NFirst, NLast
    Integer,                          Intent(In)    :: PrintLevel
    Integer,                          Intent(In)    :: Stdout
    Integer,                          Intent(In)    :: Ivbovr
    ! end of dummy parameters
    Double Precision, Dimension(:,:), Allocatable   :: fock, CC
    Integer                                         :: nao, naoao, nact, ncione
    Double Precision                                :: tuser, tsys, twall

    If (Ivbovr < 2) Return

    If (PrintLevel >= 2) Then
       Write(Stdout,'(//1X,A)') 'Integral transformation (OAO - AO)...'
       Call GetTime(tuser, tsys, twall)
    End If

    nao    = NLast(UBound(NLast,1))
    naoao  = UBound(TwoElectronAOIntegrals,1)
    nact   = UBound(ActiveCMO,2)
    ncione = nact * (nact+1) / 2

!       Write(Stdout,500) nao,naoao,nact,ncione
!   500 FORMAT(/1X,'nao,naoao,nact,ncione',4i5)

    ! Form active MO basis one-electron integrals:
    Allocate(fock(nao,nao))
    Call FormSquare(fock, FockMatrix)
    Call SimilarityTransform(OneElectronMOIntegrals, fock, ActiveCMO)
    Deallocate(fock)


    ! Form active MO basis two-electron integrals:
    If (Ivbovr > 2) Then    
       Allocate(CC(naoao,ncione))
       Call FormMOProducts(ActiveCMO, CC, NFirst, NLast)
       Call SimilarityTransform(TwoElectronMOIntegrals, TwoElectronAOIntegrals, CC)
       Deallocate(CC)
    End If

    If (PrintLevel >= 2)  Call PrintTime(tuser, tsys, twall, Stdout)

    If (PrintLevel > 5) Then
       Write(Stdout,'(//A/)') ' Transformed AO basis one-electron integrals:'
       Call PrintMatrix(OneElectronMOIntegrals, Stdout)
       If (Ivbovr > 2) Then    
          Write(Stdout,'(//A/)') ' Transformed AO basis two-electron integrals:'
          Call PrintMatrix(TwoElectronMOIntegrals, Stdout)
       End If
    End If

    Return
  End Subroutine OAOTRA



  Subroutine SimilarityTransform(C, A, B)  ! C = BT * A * B
    Implicit None
    Double Precision, Dimension(:,:)              :: C
    Double Precision, Dimension(:,:)              :: A, B
    ! end of dummy parameters
    Double Precision, Dimension(:,:), Allocatable :: P
    Integer                                       :: k, l, ldb

    k   = UBound(A,1)
    l   = UBound(B,2)
    ldb = UBound(B,1)

    Allocate(P(k,l))
    Call DGEMM('N', 'N', k, l, k, 1.D0, A(1,1), k,   B(1,1), ldb, 0.D0, P,      k)
    Call DGEMM('T', 'N', l, l, k, 1.D0, B(1,1), ldb, P(1,1), k,   0.D0, C(1,1), l)
    Deallocate(P)

    Return
  End Subroutine SimilarityTransform



  Subroutine FormSquare(square, packed)
    Implicit None
    Double Precision, Dimension(:,:) :: square
    Double Precision, Dimension(:)   :: packed
    ! end of dummy parameters
    Integer                          :: i, j, ij

    ij = 0
    Do j=1, UBound(square,2)
       Do i=1, j
          ij = ij + 1
          square(i,j) = packed(ij)
          square(j,i) = packed(ij)
       End Do
    End Do

    Return
  End Subroutine FormSquare



  Subroutine FormMOProducts(CMO, CC, NFirst, NLast)
    Implicit None
    Double Precision, Dimension(:,:) :: CMO, CC
    Integer,          Dimension(:)   :: NFirst, NLast
    ! end of dummy parameters
    Integer                          :: i, j, ij

    ij = 0
    Do i=1, UBound(CMO,2)
       Do j=1, i
          ij = ij + 1
          Call FormMOProduct(CMO, CC, NFirst, NLast, i, j, ij)
       End Do
    End Do

    Return
  End Subroutine FormMOProducts



  Subroutine FormMOProduct(CMO, CC, NFirst, NLast, i, j, ij)
    Implicit None
    Double Precision, Dimension(:,:) :: CMO, CC
    Integer,          Dimension(:)   :: NFirst, NLast
    Integer                          :: i, j, ij
    ! end of dummy parameters
    Integer                          :: ia, mu, nu, munu

    munu = 0
    Do ia=1, UBound(NFirst,1)
       Do mu=NFirst(ia), NLast(ia)
          Do nu=NFirst(ia), mu-1
             munu = munu + 1
             CC(munu,ij) = CMO(mu,i) * CMO(nu,j) + CMO(nu,i) * CMO(mu,j)
          End Do
          munu = munu + 1
          CC(munu,ij) = CMO(mu,i) * CMO(mu,j)
       End Do
    End Do

    Return
  End Subroutine FormMOProduct



  Subroutine PrintTime(tuser, tsys, twall, Stdout)
    Double Precision :: tuser, tsys, twall
    Integer :: Stdout
    ! end of dummy parameters
    Double Precision :: tuser1, tsys1, twall1, tcpu

    tuser1 = tuser
    tsys1  = tsys
    twall1 = twall

    Call GetTime(tuser, tsys, twall)

    tuser1 = tuser  - tuser1
    tsys1  = tsys   - tsys1
    twall1 = twall  - twall1
    tcpu   = tuser1 + tsys1

    Write(Stdout,'(/A,F10.2,A)') ' User   CPU time: ', tuser1, ' s'
    Write(Stdout,'( A,F10.2,A)') ' System CPU time: ', tsys1,  ' s'
    Write(Stdout,'( A,F10.2,A)') ' Total  CPU time: ', tcpu,   ' s'
    Write(Stdout,'( A,F10.2,A)') ' Wall clock time: ', twall1, ' s'

    Return
  End Subroutine PrintTime

  !
  ! This subroutine prints a rectangular matrix on unit NB6:
  !
  Subroutine PrintMatrix(A, NB6)
    Implicit None
    Double Precision, Dimension(:,:) :: A
    Integer                          :: NB6
    ! end of dummy parameters
    Integer                          :: i, j, k
    Integer                          :: ni, nj, nk, nl, nm

    ni = UBound(A,1)
    nj = UBound(A,2)
    nk = (nj+9) / 10

    Do k=1,nk
       nl = 10*k-9
       nm = Min(nj, nl+9)
       Write(NB6,'(4X,10I12)') (j,j=nl,nm)
       Do i=1,ni
          Write(NB6,'(1X,I6,10F12.6)') i, (A(i,j),j=nl,nm)
       End Do
       Write(NB6,'(//)')
    End Do

    Return
  End Subroutine PrintMatrix


End Module OAOINT
